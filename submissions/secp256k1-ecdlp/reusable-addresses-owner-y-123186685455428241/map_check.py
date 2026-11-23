"""Independent integer memory-layout checker; imports no candidate modules.

All ranges are half-open. Mapping pages and descriptor tables are counted
symbolically; no production-size enumeration or stored search answer is used.
This checks a specified mapping construction, not an Arm implementation.
"""
from dataclasses import dataclass, asdict
import json

T=1<<40
LOW=1<<55
TOP=1<<64
UP=TOP-LOW
PA_LIMIT=1<<56
F1,F2=1<<40,1<<42
Q1=103079062500
NBATCH=6080000000
CH=1<<22
BCAP=(1<<26)-1
CAP_H=38220*T//32
CAND_CAP=2082000000
REQ_CAP=125*(1<<39)
RCB_CAP=277059750328320
PROGRAM_BYTES=1<<19
SHADOW=0x20000
H2_STRIDE=3*4096
H2_UP_N=LOW//H2_STRIDE
H2_LO_N=F2-H2_UP_N
H2_UP_END=UP+H2_UP_N*H2_STRIDE
H2_LO_END=32767*T
H2_LO_BEGIN=H2_LO_END-H2_LO_N*H2_STRIDE
RCB_END=256*T+65536+64*RCB_CAP
LHR=32767*T
HCC_BYTES=(1<<20)+8*2048+4096+CAND_CAP*528
CLONE_BEGIN=((LHR+65536+HCC_BYTES+65535)//65536)*65536
CLONES=2*(22000+1)+1
CLONE_END=CLONE_BEGIN+CLONES*64
YVIEW=38*T
YQ=40*T
YWIDTH=1<<35
YBUFFER=(YQ,YQ+(1<<35)+(1<<31))
SOURCE_R=256*T-1024
SOURCE_H=LOW-4096
TABLES_R_CHARGE=51809096722
TABLES_YVIEW_CHARGE=32900

def ceildiv(a,b):return -(-a//b)
def rounded(a,b,p):return (a//p*p,ceildiv(b,p)*p)
def union(intervals):
    out=[]
    for a,b in sorted(intervals):
        if out and a<=out[-1][1]:out[-1]=(out[-1][0],max(b,out[-1][1]))
        else:out.append((a,b))
    return out
def usize(intervals):return sum(b-a for a,b in union(intervals))
def cover(a,b,shift):return (a>>shift,((b-1)>>shift)+1)
def disjoint(intervals):
    s=sorted(intervals)
    return all(a<b and (0<=a<b<=LOW or UP<=a<b<=TOP) for a,b in s) and all(b<=c for (a,b),(c,d) in zip(s,s[1:]))

@dataclass(frozen=True)
class Region:
    role:str
    name:str
    begin:int
    end:int
    kind:str='normal'
    def iv(self):return self.begin,self.end

def registry(rcb_cap=RCB_CAP,clone_begin=CLONE_BEGIN):
    raw=[
        ('G','LC1',0,F1*8),
        ('G','YR',32*T+(1<<31)-64,32*T+258*(1<<31)-64),
        ('G','EB',42*T+(1<<37)-6*CH*512,42*T+(1<<37)),
        ('G','YQ',YQ-64,YQ+(1<<36)+(1<<31)),
        ('G','LG',34*T,34*T+4096,'device'),
        ('G','LGM',34*T+4096,34*T+8192,'device'),
        ('G','PY',34*T+0xb000,34*T+0xc000,'device'),
        ('G','G34',34*T+65536,38*T),
        ('G','YVIEW',YVIEW,YVIEW+YWIDTH,'device'),
        ('G','CN',84*T,84*T+64*Q1),
        ('G','CNM',90*T,90*T+64*NBATCH),
        ('G','REQ',96*T-64,96*T+64*REQ_CAP),
        ('G','AXY',4096*T,24576*T),
        ('G','XT6',UP,TOP-2048*T),
        ('G','H1L',24576*T,28672*T,'alias2'),
        ('G','H1U',28672*T,LOW,'alias2'),
        ('G','AL',16*T,16*T+4096),
        ('G','BT',TOP-2048*T,TOP),
        ('R','RCS',65536,65536+(1<<20)+8*2048),
        ('R','AL',16*T,16*T+65536),
        ('R','LC2',32*T,64*T),
        ('R','RC2',96*T,128*T),
        ('R','E2P',128*T,192*T),
        ('R','DPW',256*T-4096,256*T),
        ('R','RCB',256*T,256*T+65536+64*rcb_cap),
        ('R','H2L',H2_LO_BEGIN,H2_LO_END,'alias3'),
        ('R','LR',LHR,LHR+32768,'device'),
        ('R','LRM',LHR+32768,LHR+65536,'device'),
        ('R','RCT',LOW-8192,LOW-4096),
        ('R','H2U',UP,H2_UP_END,'alias3'),
        ('H','R2L',0,CAP_H*32-5460*T),
        ('H','LH',LHR,LHR+32768,'device'),
        ('H','LHM',LHR+32768,LHR+65536,'device'),
        ('H','HCC',LHR+65536,LHR+65536+HCC_BYTES),
        ('H','CLONES',clone_begin,clone_begin+CLONES*64),
        ('H','HCS',LOW-16384,LOW),
        ('H','HH',UP,TOP-24576*T),
        ('H','NX',TOP-24576*T,TOP-24576*T+CAP_H*16),
        ('H','R2U',TOP-5460*T,TOP),
    ]
    return [Region(*x) for x in raw]

def granule(role,begin):
    return 4096 if role=='R' or role=='G' and begin<LOW else 65536

def rounded_regions(regions):
    # Equal-attribute neighbors may share a page; unlike memory types may not.
    return [(r,*rounded(r.begin,r.end,granule(r.role,r.begin))) for r in regions]

def mapping_compatible(regions):
    for i,(a,lo,hi) in enumerate(rounded_regions(regions)):
        for b,lo2,hi2 in rounded_regions(regions)[i+1:]:
            if a.role==b.role and lo<hi2 and lo2<hi:
                if a.kind!=b.kind:return False
    return True

def table_count(regions,role,g):
    """D128 tree counts, keeping table pages at both roots.

    64KiB: root on bits54:52, 4096 entries per subsequent table; Normal
    interiors use256MiB blocks. 4KiB: root on bits54:52, 256 entries/table;
    Normal interiors use64GiB,256MiB,1MiB blocks where aligned. Device and
    alternating guard/data intervals use4KiB pages throughout.
    """
    selected=[(r,a,b) for r,a,b in rounded_regions(regions)
              if r.role==role and granule(role,a)==g]
    roots=len({a>=LOW for r,a,b in selected})
    shifts=(52,40,28) if g==65536 else (52,44,36,28,20)
    levels=[[] for _ in shifts]
    for r,a,b in selected:
        for j,s in enumerate(shifts):
            mandatory=j<(2 if g==4096 else 2)
            if mandatory or r.kind!='normal':levels[j].append(cover(a,b,s))
            else:
                if a%(1<<s):levels[j].append((a>>s,(a>>s)+1))
                if b%(1<<s):levels[j].append(((b-1)>>s,((b-1)>>s)+1))
    return {'root_tables':roots,'levels':[usize(x) for x in levels],
            'total':roots+sum(usize(x) for x in levels),'page_bytes':g}

def physical_pack(regions,role,pt_bytes):
    """Construct disjoint PA runs with block-compatible VA/PA offsets.

    Normal data is packed at matching1MiB offsets for4KiB ranges and matching
    256MiB offsets for64KiB ranges. To use larger4KiB blocks, it is further
    aligned to64GiB; this only wastes bounded unbacked PA gaps. Alternating
    guard regions use compressed4KiB RW pages and alias AL's first physical
    page for every guard. Device apertures are conservatively charged their
    whole PA extent, although they contain no copied payload storage.
    """
    runs=[];cursor=0
    groups={}
    for r,a,b in rounded_regions(regions):
        if r.role!=role:continue
        groups.setdefault((r.kind,granule(role,a)),[]).append((a,b))
    for (kind,g),ivs in groups.items():
        for a,b in union(ivs):
            if kind=='alias2':size=(b-a)//2;align=4096;offset=0
            elif kind=='alias3':
                assert (b-a)%12288==0
                size=(b-a)*2//3;align=4096;offset=0
            else:
                size=b-a;align=(1<<36 if g==4096 else 1<<28) if kind=='normal' else g
                offset=a%align
            pa=ceildiv(max(0,cursor-offset),align)*align+offset
            assert pa>=cursor and pa%g==0 and size%g==0
            runs.append({'kind':kind,'va':[a,b],'pa':[pa,pa+size],
                         'bytes':size,'alignment':align,'offset':offset})
            cursor=pa+size
    table_start=ceildiv(cursor,65536)*65536
    end=table_start+pt_bytes
    return {'runs':runs,'data_bytes':sum(r['bytes'] for r in runs),
            'alignment_gaps':table_start-sum(r['bytes'] for r in runs),
            'tables_pa':[table_start,end],'pa_end':end,'headroom':PA_LIMIT-end}

def alias_translation(regions,physical,role,address):
    r=next(r for r in regions if r.role==role and r.begin<=address<r.end)
    if r.kind.startswith('alias'):
        stride=8192 if r.kind=='alias2' else 12288
        rw=stride-4096;index,off=divmod(address-r.begin,stride)
        if off>=rw:
            ar=next(x for x in physical['runs'] if x['va'][0]<=16*T<x['va'][1])
            return ar['pa'][0]+16*T-ar['va'][0]+off-rw,'ro'
        ar=next(x for x in physical['runs'] if x['va'][0]<=address<x['va'][1])
        # A merged alias run may contain both G H1 halves.
        index,off=divmod(address-ar['va'][0],stride)
        return ar['pa'][0]+index*rw+off,'rw'
    ar=next(x for x in physical['runs'] if x['va'][0]<=address<x['va'][1])
    return ar['pa'][0]+address-ar['va'][0],('device-ro' if r.name=='YVIEW' else 'device' if r.kind=='device' else 'rw')

def audit(program_insns=None):
    regs=registry();checks={}
    checks['range_and_va_disjoint']=all(disjoint([r.iv() for r in regs if r.role==k]) for k in 'GRH')
    checks['rounded_attributes_disjoint']=mapping_compatible(regs)
    checks['relay_capacity_rounded']=RCB_END%65536==0 and RCB_END==H2_LO_BEGIN//65536*65536
    checks['relay_blocks_exact']=H2_UP_N+H2_LO_N==F2 and H2_UP_END<=TOP and H2_LO_END==LHR
    checks['clone_capacity']=CLONES==44003 and CLONE_BEGIN%65536==0 and CLONE_END<LOW-16384
    checks['nx_capacity']=TOP-24576*T+CAP_H*16<TOP-5460*T
    checks['source_tables']=SOURCE_R%64==SOURCE_H%64==0 and SOURCE_R+512<256*T-192 and SOURCE_H+1024<=LOW-2048
    # HCS counters each have isolated2KiB monitor granules. The table, state,
    # and context occupy different granules. Ordinary stack reservations are
    # within HCC's8*2KiB space, separately from HCS.
    checks['home_monitor_granules']=len({(LOW-d)//2048 for d in (8192,6144,4096,160)})==4
    checks['program_reserved']=PROGRAM_BYTES<=(1<<20)
    checks['shadow_encodable']=SHADOW%4096==0 and SHADOW//4096<4096
    checks['source_buffers']=all(b%(1<<31)==0 and YQ-64<=b and b+YWIDTH<=YQ+(1<<36)+(1<<31) for b in YBUFFER)
    checks['view_alignment']=YVIEW%(1<<31)==0 and YWIDTH==512*(BCAP+1) and YVIEW+YWIDTH<40*T-64
    checks['offer_port_mapped']=34*T+0xb000<=34*T+0xb800 and 34*T+0xb800+48<=34*T+0xc000
    ystate,ydone,ystack=(34*T+d for d in (0x7100000,0x7120000,0x7140000))
    checks['handler_state_stack']=34*T+65536<ystate and ystate+4096<ydone and ydone+2048<=ystack-4096 and ystack<38*T
    tables={k:{str(g):table_count(regs,k,g) for g in (4096,65536)} for k in 'GRH'}
    # Existing G/H table cost stays in the reference. Only extensions are added.
    baseline_g_4k=8623622175
    baseline_g_64k=32777
    baseline_h_64k=65543
    oldregs=[r for r in regs if r.name not in ('YVIEW','CLONES')]
    checks['reference_table_counts']=(table_count(oldregs,'G',4096)['total']==baseline_g_4k
        and table_count(oldregs,'G',65536)['total']==baseline_g_64k
        and table_count(oldregs,'H',65536)['total']==baseline_h_64k)
    added_g=max(0,tables['G']['4096']['total']-baseline_g_4k)
    added_h=max(0,tables['H']['65536']['total']-baseline_h_64k)
    checks['g_tables_allowance']=tables['G']['65536']['total']<=baseline_g_64k and added_g<=TABLES_YVIEW_CHARGE
    checks['r_tables_allowance']=tables['R']['4096']['total']<=TABLES_R_CHARGE and tables['R']['65536']['total']==0
    # Charge a conservative all-new64KiB table per256MiB clone span plus
    # parent/root tables, descriptor16-insn run setup and cursor reset.
    clone_pages=ceildiv(CLONE_END,65536)-CLONE_BEGIN//65536
    clone_table_init=(added_h+4)*(4*4096+16)+2
    checks['clone_cold_allowance']=clone_table_init<=1<<23
    # Physical address reservation uses actual table counts; R retains the
    # larger charged count so no unpriced table can hide in the PA bound.
    packs={}
    for k in 'GRH':
        t4=tables[k]['4096']['total'];t64=tables[k]['65536']['total']
        if k=='R':t4=TABLES_R_CHARGE
        packs[k]=physical_pack(regs,k,t4*4096+t64*65536)
    checks['pa_limits']=all(p['pa_end']<PA_LIMIT for p in packs.values())
    # Check edge/split adjacency with a translation independent of the block
    # formula: first/lastdata pages map to distinct RW pages, allguardpages
    # alias AL. This is symbolic in destination count and samples its edges.
    alias_rows=[]
    for name,start,n in (('H2U',UP,H2_UP_N),('H2L',H2_LO_BEGIN,H2_LO_N)):
        for i in (0,1,n-2,n-1):
            b=start+i*12288
            rw0,p0=alias_translation(regs,packs['R'],'R',b)
            rw1,p1=alias_translation(regs,packs['R'],'R',b+4096)
            guard,pg=alias_translation(regs,packs['R'],'R',b+8192)
            shared,ps=alias_translation(regs,packs['R'],'R',16*T)
            assert rw0%4096==rw1%4096==guard%4096==0
            alias_rows.append(dict(region=name,index=i,rw_pa=[rw0,rw1],guard_pa=guard,
                ok=rw1==rw0+4096 and guard==shared and (p0,p1,pg,ps)==('rw','rw','ro','rw')))
    checks['relay_page_permissions']=all(x['ok'] for x in alias_rows)
    yback=[]
    for b in YBUFFER:
        pa,_=alias_translation(regs,packs['G'],'G',b)
        end,_=alias_translation(regs,packs['G'],'G',b+YWIDTH-1)
        yback.append(dict(va=[b,b+YWIDTH],pa=[pa,end+1],aligned64=pa%64==0))
    checks['yview_read_only']=all(alias_translation(regs,packs['G'],'G',a)[1]=='device-ro' for a in (YVIEW,YVIEW+YWIDTH-1))
    checks['y_backing_contiguous']=all(x['pa'][1]-x['pa'][0]==YWIDTH and x['aligned64'] for x in yback)
    checks['owner_drain_monitor_isolated']=(ydone%2048==0 and (ydone+2048)//2048!=ydone//2048 and (ydone+2048)//2048!=ystate//2048 and (ydone+2048)//2048!=ystack//2048)
    program_insns=program_insns or {'route':73098,'leaf':30320,'secp_fe':3690,
        'ytransport':137,'pointer_init':50,'relay_source_init':10,'home_source_init':10}
    code_bytes=4*sum(program_insns.values())
    checks['program_fits']=code_bytes<=PROGRAM_BYTES
    negatives={
        'old_relay_cap_overlap':not disjoint([r.iv() for r in registry((1<<49)-(1<<43)+(1<<41)) if r.role=='R']),
        'clone_overlaps_hcs':not disjoint([r.iv() for r in registry(clone_begin=LOW-16384) if r.role=='H']),
        'old_program_cap_rejected':code_bytes>(1<<18),
        'old_guard_offset_is_writable':alias_translation(regs,packs['R'],'R',UP+4096)[1]=='rw',
        'missing_second_buffer_gap':YBUFFER[1]!=YBUFFER[0]+YWIDTH,
    }
    assert all(checks.values()),{k:v for k,v in checks.items() if not v}
    assert all(negatives.values())
    return {'all_checks_pass':True,'checks':checks,'negative_controls':negatives,
        'regions':[asdict(r) for r in regs],
        'relay':{'stride':H2_STRIDE,'upper_count':H2_UP_N,'lower_count':H2_LO_N,
            'upper':[UP,H2_UP_END],'lower':[H2_LO_BEGIN,H2_LO_END],
            'cap':RCB_CAP,'rcb_end':RCB_END,'rcb_h2_gap':H2_LO_BEGIN-RCB_END,
            'tables_charged':TABLES_R_CHARGE,'table_init_instructions':1040*TABLES_R_CHARGE},
        'tables':tables,'physical':packs,'relay_page_edge_checks':alias_rows,
        'clone':{'interval':[CLONE_BEGIN,CLONE_END],'rows':CLONES,'pages64k':clone_pages,
            'new_tables64k':added_h,'init_upper':clone_table_init,'charged':1<<23},
        'yview':{'interval':[YVIEW,YVIEW+YWIDTH],'source_buffers':yback,
            'handler_state':ystate,'completion_word':ydone,'drained_word':ydone+2048,'stack_top':ystack,
            'offer_port':[34*T+0xb800,34*T+0xb800+48],
            'extra_tables4k':added_g,'table_allowance':TABLES_YVIEW_CHARGE,
            'table_init_instructions':1040*TABLES_YVIEW_CHARGE},
        'program':{'instructions':program_insns,'bytes':code_bytes,'allowance':PROGRAM_BYTES,
            'shadow':SHADOW,'scope':'Caller must supply final expanded instruction counts; defaults are prior composed streams.'}}

if __name__=='__main__':print(json.dumps(audit(),indent=2))

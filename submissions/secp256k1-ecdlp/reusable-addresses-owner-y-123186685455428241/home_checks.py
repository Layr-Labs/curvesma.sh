"""Reusable direct-HH/home source-table checker; llm-verifier execution only."""
import re,time

def reference_geometry(source):
    start=source.index('6:\nadd x13,x14,x2,asr #1\nDLVHPR')
    end=source.index('.macro Z1',start)
    full=['6:','add x13,x14,x2,asr #1']+[f'DLVHPR {r}'for r in('x0','x1','x7','x8','x24','x25','x26','x27')]
    full+=['add x1,x27,#64']+['DLVHP']*119+['DLVHPL','add x1,x0,#64','staddl x28,[x29]','ret']
    source=source[:start]+'\n'.join(full)+'\n'+source[end:]
    lo,hi=source.index('.globl dlvh\n'),source.index('.macro Z1')
    section=source[lo:hi].replace('tbnz x16,#12,6f','tbnz x16,#13,6f').replace('mov x7,#4096','mov x7,#8192')
    return re.sub(r'0x5000\b','0x20000',source[:lo]+section+source[hi:])

def run(E,source,original_base,mode='compact'):
    started=time.time()
    M=E.machine;A=E.a64;K=256;TABLE=(1<<55)-4096;NREF=M.NEXT_BASE
    NX=(1<<64)-24576*(1<<40)
    CAP_H=1313229200424960
    assert NX+CAP_H*16 < (1<<64)-5460*(1<<40)
    HCC=(1<<20)+8*2048+4096+2082000000*528
    CLONE=(M.LINK+0x10000+HCC+65535)&~65535
    CLONE_END=CLONE+44003*64
    assert CLONE_END < (1<<55)-16384
    INIT='''init_table:
    movn x0,#4095
    movk x0,#127,lsl #48
    mov x1,#0x7fff0000000000
    add x1,x1,#64
    mov x2,#128
    1:
    str x1,[x0],#8
    add x1,x1,#64
    subs x2,x2,#1
    b.ne 1b
    ret
    '''
    NEW=source
    OLD=reference_geometry(original_base)
    bi,bl=A.a64_load(OLD);ci,cl=A.a64_load(NEW);ii,il=A.a64_load(INIT)
    FIRST_COLD=min(v[0]for k,v in cl.items()if k.startswith('h31_cold'))


    def metas(mem,ptr):return (mem[ptr+16],mem[ptr+24])


    class Home:
        def __init__(self,changed,cap=4096,mutate=None):
            self.changed=changed;self.ops,self.labs=(ci,cl)if changed else(bi,bl)
            if mutate is not None:self.ops,self.labs=mutate
            self.nd=M.ND();self.mem=self.nd.m;self.regs={};self.events=[];self.total=0
            self.init_count=E.arun(ii,il,'init_table',{},self.mem)
            assert self.init_count==518
            for b in range(2048):self.mem[M.HT+8*b]=0
            for j,v in enumerate((M.Node.CL,M.Node.CL+16*cap,NX if changed else NREF,CLONE if changed else(-2*NREF)&E.M64)):
                self.mem[M.ST+8*j]=v
            E.arun(self.ops,self.labs,'dlvh0',dict(x29=M.CNT),self.mem,out_regs=self.regs)
            self.start_nx=self.regs.get('x13',0)
        def candidates(self):
            return [(metas(self.mem,self.mem[p]),metas(self.mem,self.mem[p+8]))
                    for p in range(M.Node.CL,self.regs['x14'if self.changed else'x3'],16)]
        def heads(self):
            result={}
            for b in range(2048):
                h=self.mem[M.HT+8*b];rows=[];seen=set()
                while h:
                    assert h not in seen,('cycle',b,h);seen.add(h)
                    if not self.changed:
                        ptr=(-2*NREF+2*h)&E.M64;rows.append(metas(self.mem,ptr));h=self.mem[h]
                    else:
                        if h&8:nxt,pair=self.mem[h-8],self.mem[h]
                        else:nxt,pair=0,h
                        for off in (32,0):
                            if self.mem[pair+off+8]&((1<<50)-1)==b:rows.append(metas(self.mem,pair+off))
                        h=nxt
                result[b]=rows
            return result
        def send(self,base,n,prefix,records):
            assert base%8192==0 and 0<=prefix<=n<=K
            self.mem[M.Node.E2]=base+32*prefix
            msg=[0,base,32*n,M.Node.E2]+[v for row in records for v in row]
            # Independent ordered-record inventory before this message. Every
            # newrecord is inserted, even after the candidatearray fills.
            inventory=self.heads();prior=0
            for row in records[prefix:n]:
                bucket=row[1]&((1<<50)-1)
                prior+=len(inventory[bucket]);inventory[bucket].append((row[2],row[3]))
            out={};vis={i:0 for i in range(len(self.ops))}
            count=E.arun(self.ops,self.labs,'dlvh',self.regs,M.LV(self.mem,msg,[]),out_regs=out,visits=vis,
                ls64=([M.WIN,(TABLE,TABLE+1024),(base,base+8192),(CLONE,CLONE_END)],
                      [(base,base+8192),(CLONE,CLONE_END)]))
            self.regs=out;self.total+=count
            assert out['x0']==M.DHP[0] and(not self.changed or out['sp']==TABLE),('source state',out['x0'],out['sp'])
            if self.changed:
                assert NX<=out['x13']<=NX+16*4096
                assert CLONE<=out['x25']<=CLONE_END
            hh=sum(v for p,v in vis.items() if v and self.ops[p][0]=='swp')
            cold=verify=listed=0
            if self.changed:
                cold=sum(v for p,v in vis.items()if p>=FIRST_COLD)
                verify=sum(vis[v[0]]for key,v in self.labs.items()if key.startswith('h31_')and '_hit'in key)
                listed=sum(v for p,v in vis.items()if self.ops[p]==('stp',['x12','x18','[x14]','#16']))
                assert cold<=10*prior+5*verify+6*listed,('cold count bound',cold,prior,verify,listed)
                completed=sum(v for p,v in vis.items()if self.ops[p]==('staddl',['x28','[x29]']))
                nonhash=count-3*hh-cold-completed
                predicted=home_nonhash(prefix,n,base==0)
                assert completed==1 and hh==n-prefix,('completion/hash population',completed,hh,n,prefix)
                assert nonhash==predicted,('nonhash path cost',base,n,prefix,nonhash,predicted)
            self.events.append(dict(base=base,n=n,prefix=prefix,count=count,hashed=hh,
                prior_pairs=prior,cold=cold,verify=verify,listed=listed,
                nonhash=nonhash if self.changed else None,
                nx=(out['x13']-NX)//16 if self.changed else None,
                clones=(out['x25']-CLONE)//64 if self.changed else None))
            return count


    def records(style,seed=0):
        rows=[]
        for j in range(K):
            if style=='empty':u=seed*512+j;v=j<<41
            elif style=='collision':u=j%13;v=(j%7)<<41
            elif style=='samepair':u=j//2%13;v=(j//3%3)<<41
            else:u=(j%13)+(j%3<<50);v=(j%7)<<41
            rows.append((v,u,0xABC000+seed*1024+j,0xDEF000+seed*1024+j))
        return rows


    # Execute the unchanged source macros to obtain primitive costs. The
    # standalone RET is removed; HCHGI's three hot instructions are removed
    # per actually executed SWP. Fresh distinct buckets make cold work zero.
    macros='\n'.join(re.findall(r'(?ms)^\.macro .*?^\.endm',NEW))
    wrappers='\n'.join(f'cost_{m.lower()}:\n{m}\nret'for m in('DLVS','DLVP','DLVHS','DLVHP'))
    pi,pl=A.a64_load(macros+'\n'+wrappers)
    primitive={}
    for label,zero in(('dlvs',False),('dlvp',False),('dlvhs',False),('dlvhp',False),('dlvhp',True)):
        h=Home(True);base=0 if zero else M.Node.R2
        regs=dict(h.regs,x1=M.DHP[0],x24=base);vis={i:0 for i in range(len(pi))}
        payload=[0,base,64,M.Node.E2]+[v for row in records('empty')for v in row]
        count=E.arun(pi,pl,'cost_'+label,regs,M.LV(h.mem,payload,[]),visits=vis,
            ls64=([M.WIN,(base,base+8192),(CLONE,CLONE_END)],[(base,base+8192),(CLONE,CLONE_END)]))
        hashes=sum(v for p,v in vis.items()if pi[p][0]=='swp')
        cold_start=min(v[0]for k,v in pl.items()if k.startswith('h31_cold'))
        assert not any(v for p,v in vis.items()if p>=cold_start)
        primitive[label+('_zero'if zero else'')]=count-1-3*hashes
    assert primitive==dict(dlvs=4,dlvp=4,dlvhs=12,dlvhp=6,dlvhp_zero=10),primitive

    seg_source=re.search(r'(?ms)^\.macro SEGL .*?^\.endm',NEW).group()
    seg_ops=[s.split()[0]for line in seg_source.splitlines()if(s:=line.strip())
             and not s.startswith('.')and not s.endswith(':')]
    assert seg_ops==['cbz','tbz','\\s','sub','lsr','cbz','\\p','subs','b.ne','tbz','\\s']
    fallback=NEW.split('.globl dlvh\n',1)[1].split('\n6:\n',1)[0]
    fixed_ops=[s for line in fallback.splitlines()if(s:=line.strip())
               and not s.endswith(':')and not s.startswith('SEGL ')]
    assert len(fixed_ops)==17 and fixed_ops.count('staddl x28,[x29]')==1
    assert fixed_ops[5]=='tbnz x16,#13,6f'
    fixed_without_completion=len(fixed_ops)-1

    def segment(start,length,single,pair):
        if length==0:return 1  # CBZ only.
        first=start%2;remaining=length-first;pairs,last=divmod(remaining,2)
        # CBZ,TBZ,LSR,CBZ,TBZ; leading singleton also executes SUB;
        # each pair executes SUBS and B.NE, including the last pair.
        return 5+first*(single+1)+pairs*(pair+2)+last*single

    def home_nonhash(prefix,end,zero=False):
        assert 0<=prefix<=end<=256
        if(prefix,end)==(0,256):return 424+6*zero
        # Actual header/prologue/between-segment/RET instructions, excluding
        # the one completion STADDL. Prefix>0 skips the full-size TBNZ.
        fixed=fixed_without_completion-(prefix>0)
        return(fixed+segment(0,prefix,primitive['dlvs'],primitive['dlvp'])
            +segment(prefix,end-prefix,primitive['dlvhs'],primitive['dlvhp'])
            +segment(end,256-end,primitive['dlvs'],primitive['dlvp'])
            +(primitive['dlvhp_zero']-primitive['dlvhp'])*(zero and prefix==0 and end>=2))

    def inherited_dsg(start,length):
        if length==0:return 1
        first=start%2;remaining=length-first
        return 5+first*5+(remaining//2)*6+(remaining%2)*4-(remaining==0)

    # Enumerate the actual K256 message state space: all 0<=prefix<=end<=256.
    # This is an exact integer census, not an execution extrapolation. The
    # formula is bound to actual full-handler executions below and on every
    # existing semantic/cold-path case through Home.send's exact assertion.
    lattice=[];legacy_under=0
    for prefix in range(257):
        for end in range(prefix,257):
            if(prefix,end)==(0,256):continue
            legacy=23-(prefix>0)+sum(inherited_dsg(a,b)for a,b in
                ((0,prefix),(prefix,end-prefix),(end,256-end)))
            literal=23-(prefix>0)+sum(segment(a,b,4,4)for a,b in
                ((0,prefix),(prefix,end-prefix),(end,256-end)))
            legacy_under=max(legacy_under,literal-legacy)
            lattice.append((home_nonhash(prefix,end),home_nonhash(prefix,end,True),legacy,literal,prefix,end))
    maxima=[max(row[i]for row in lattice)for i in range(4)]
    assert len(lattice)==33152 and maxima==[1072,1072,811,811]and legacy_under==1,(len(lattice),maxima,legacy_under)
    assert max(row[0]-row[2]for row in lattice)==262
    assert max(row[1]-row[2]for row in lattice)==262

    # An old nonempty prefix can make the first full block of a group a
    # fallback. All subsequent full blocks are regular; there is at most one
    # final partial. Charge BOTH boundary messages in the same group, without
    # borrowing a source-table credit from an earlier weighted group.
    # First full transport: FLUSH2 21 + exception24 + DONE_C8 =53.
    # Final transport: FINH_DEST16 + FINH_SENT5 + DONE_C8 + MSG_LINK14 =43.
    # DONE_C includes the handler's one completion instruction removed above.
    # The maxima already include singleton and pair-zero clones. Adding six
    # more zero-repair instructions and reset2 is conservative.
    boundary_allowance=1153+2048
    actual_bound=53+1072+43+1072+6+2
    round_bound=53+1200+43+1200+6+2
    assert actual_bound==2248 and round_bound==2504 and round_bound<=boundary_allowance
    full_allowance=565-88+1
    aggregate_min=None
    for old_prefix in range(256):
        for added in range(768):
            full=(old_prefix+added)//256
            ordinary=max(0,full-1)
            assert ordinary*256<=added
            slack=full_allowance*added+256*boundary_allowance-256*(ordinary*full_allowance+round_bound)
            aggregate_min=slack if aggregate_min is None else min(aggregate_min,slack)
    assert aggregate_min is not None and aggregate_min>=0

    # Short actual sequences include the maximum, no-new-record messages,
    # odd-start singleton, both endpoint singletons, zero repair and a fallback
    # immediately following a real partial final. Prefix records are actually
    # delivered first; they are not fabricated only as a watermark.
    boundary_events=[]
    boundary_cases=[(M.Node.R2,1,255),(0,1,255),(0,0,2),(M.Node.R2,0,0)]
    if mode!='compact':
        boundary_cases += [(M.Node.R2,1,2),(0,0,255),(M.Node.R2,255,256),
                           (M.Node.R2,256,256),((1<<64)-8192,255,256)]
    for base,prefix,end in boundary_cases:
        h=Home(True);rows=records('empty')
        if prefix:h.send(base,prefix,0,rows)
        h.send(base,end,prefix,rows)
        boundary_events.append(h.events[-1])
    assert max(e['nonhash']for e in boundary_events)==1072
    partial_cost=dict(primitives=primitive,fixed_nonhash=fixed_without_completion,
        prefix_states=len(lattice),max_nonhash=1072,
        max_inherited_formula=811,max_literal_old_copy=811,inherited_singleton_undercount=1,
        max_extra_over_inherited=262,per_group_boundary=actual_bound,
        rounded_boundary=round_bound,charged_boundary=boundary_allowance,ordinary_full_allowance=full_allowance,
        aggregate_cases=256*768,aggregate_min_slack_numerator=aggregate_min,aggregate_slack_denominator=256,
        actual_boundary_cases=[{k:e[k]for k in('base','prefix','n','hashed','count','cold','nonhash')}for e in boundary_events])


    def case(style,seq,cap=4096,mutate=None):
        a=Home(False,cap);b=Home(True,cap,mutate);nnew=0
        for base,n,prefix,seed in seq:
            rows=records(style,seed)
            a.send(base,n,prefix,rows);b.send(base,n,prefix,rows);nnew+=n-prefix
            assert a.candidates()==b.candidates(),('candidate mismatch',style,base,n,prefix,a.candidates()[-3:],b.candidates()[-3:])
            assert a.heads()==b.heads(),('head mismatch',style,base,n,prefix)
            assert all(a.mem[base+8*j]==b.mem[base+8*j] for j in range(4*K)),('R2 mismatch',base,n,prefix)
            assert (b.regs['x13']-NX)//16<=nnew
        return dict(style=style,cap=cap,base_events=a.events,new_events=b.events,
                    delta=b.total-a.total,candidates=len(b.candidates()),indexed=nnew)


    R=M.Node.R2
    cases=[case('samepair',[(0,256,0,0)],5),
           case('collision',[(R,1,0,0),(R,3,1,0),(R,7,3,0),(R,256,7,0),(R+8192,256,0,1)]),
           case('bucket_only',[(0,1,0,0),(0,2,1,0),(0,3,2,0),(0,256,3,0)],0),
           case('empty',[(R,256,0,0)])]
    if mode!='compact':
        cases += [case(style,[(R,256,0,0)],cap) for style in('collision','samepair','bucket_only') for cap in(0,4096)]
    negative={}
    for kind in('dummy_same_bucket','zero_unfixed','drop_samepair_lane0'):
        bad=list(ci);changes=0
        if kind=='dummy_same_bucket':
            for pc,(op,args)in enumerate(ci):
                if(op,args)==('eor',['x12','x17','#1']):bad[pc]=('mov',['x12','x17']);changes+=1
        elif kind=='zero_unfixed':
            for pc,(op,args)in enumerate(ci):
                if(op,args)==('cbz',['x24','dlvh_zero_full']):bad[pc]=('cbnz',list(args));changes+=1
        else:
            for key,v in cl.items():
                if key.startswith('h31_selfdirectbody_p0'):
                    suffix=key[len('h31_selfdirectbody_p0'):]
                    bad[v[0]]=('b',['h31_done'+suffix]);changes+=1
        assert changes>0,kind
        seq=[(0,256,0,0)]if kind=='zero_unfixed'else[(0,1,0,0),(0,3,1,0),(0,256,3,0)]
        try:case('samepair',seq,mutate=(bad,cl))
        except(AssertionError,KeyError)as exc:negative[kind]=str(exc)[:120]
        else:raise AssertionError('mutation not detected: '+kind)
    reset=Home(True);reset.send(R,256,0,records('collision'));out={}
    nreset=E.arun(ci,cl,'home_join_reset',reset.regs,reset.mem,out_regs=out)
    assert nreset==2 and(out['x13'],out['x25'])==(NX,CLONE)
    full=next(c['new_events'][0]for c in cases if c['style']=='empty')
    assert full['count']==1193 and full['cold']==0 and full['hashed']==256
    return dict(ok=True,cases=len(cases),messages=sum(len(c['new_events'])for c in cases),
                indexed=sum(c['indexed']for c in cases),negative=negative,clone_region=[CLONE,CLONE_END],
                nx_arena=[NX,NX+CAP_H*16],empty_full_instructions=full['count'],copy_allowance=424,
                cold_bound=10,reset_count=nreset,partial_cost=partial_cost,seconds=round(time.time()-started,3))

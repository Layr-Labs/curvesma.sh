import sys,json,random,re,io,contextlib,runpy,hashlib
from fractions import Fraction as F
from pathlib import Path
root=Path(__file__).resolve().parent
if(root/'files').is_dir():root=root/'files'
sys.path.insert(0,str(root));import a64,flip_witness as fw,flip_check as fc
cfg=json.load(sys.stdin);seed=cfg['seed'];rng=random.Random(seed)
# protected target constants (target-profile.json): P-256 field, a=-3, b, G, prime order r
P=0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff;ORD=0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551
GX,GY=0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296,0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5
assert fw.P==fc.P==P and(GY*GY-GX**3+3*GX-fc.BC)%P==0 and all(pow(b,n-1,n)==1 for n in(P,ORD)for b in(2,3,5,7,11,13))
def jadd(A,B):
 if A is None:return B
 if B is None:return A
 if A[0]==B[0]and(A[1]+B[1])%P==0:return None
 m=(3*A[0]*A[0]-3)*pow(2*A[1],-1,P)%P if A==B else(B[1]-A[1])*pow(B[0]-A[0],-1,P)%P
 z=(m*m-A[0]-B[0])%P;return z,(m*(A[0]-z)-A[1])%P
def smul(k,A):
 Z=None
 for b in bin(k)[2:]:Z=jadd(Z,Z);Z=jadd(Z,A)if b=='1'else Z
 return Z
assert smul(ORD,(GX,GY))is None and smul(ORD-1,(GX,GY))==(GX,(-GY)%P)
cold_accounting_check=__import__('cold_accounting_check');routing_check=__import__('routing_check')
cold=cold_accounting_check.check();routing=routing_check.check()
assert routing['formulas_ok']and routing['cands_ok']
assert routing['msgs_ok']and routing['slots_ok']and routing['tags_ok']and routing['chain_ok']and routing['stop_ok']and routing['all_paths']
assert routing['signaled']and routing['done_markers']==64*routing['groups_run']and routing['completions']==16*routing['groups_run']
isa=0
for shift in[0,1,36,63]:
 for av in[0,1,(1<<64)-1,1<<63]:
  for bv in[0,1,(1<<64)-1,1<<63]:
   pp,ll=a64.a64_load('check:\n orr x0,x1,x2,lsl #'+str(shift)+'\n ret\n');oo={}
   ct=a64.a64_run(pp,ll,'check',dict(x1=av,x2=bv,x30=777),{},out_regs=oo)
   assert ct==2 and oo['x0']==(av|((bv<<shift)&((1<<64)-1)));isa+=1
for lo,width in[(0,64),(35,1),(36,9),(45,19),(63,1)]:
 for v in[0,1,(1<<64)-1,1<<63]:
  pp,ll=a64.a64_load('check:\n ubfx x0,x1,#'+str(lo)+',#'+str(width)+'\n ret\n');oo={}
  ct=a64.a64_run(pp,ll,'check',dict(x1=v,x30=777),{},out_regs=oo)
  assert ct==2 and oo['x0']==((v>>lo)&((1<<width)-1));isa+=1
for lo,width in[(0,36),(3,36),(36,9),(45,19)]:
 for dv,sv in[(0,(1<<64)-1),((1<<64)-1,0),(1<<39,(1<<36)-1),(0xa5a5a5a5a5a5a5a5,0x5a5a5a5a5a5a5a5a)]:
  pp,ll=a64.a64_load('check:\n bfi x0,x1,#'+str(lo)+',#'+str(width)+'\n ret\n');oo={}
  ct=a64.a64_run(pp,ll,'check',dict(x0=dv,x1=sv,x30=777),{},out_regs=oo);mk=((1<<width)-1)<<lo
  assert ct==2 and oo['x0']==((dv&~mk)|((sv<<lo)&mk))&((1<<64)-1);isa+=1
M64=(1<<64)-1
for av in[0,1,(1<<64)-1,1<<63,0xa5a5a5a5a5a5a5a5]:
 for bv in[0,1,(1<<64)-1,1<<63,0x5a5a5a5a5a5a5a5a]:
  for prog_,want in[('add x0,x1,x2,lsr #28',(av+(bv>>28))&M64),('add x0,x1,x2,lsl #3',(av+(bv<<3))&M64),('orr x0,x1,x2,lsr #28',av|(bv>>28)),('extr x0,x1,x2,#36',((av<<64|bv)>>36)&M64)]:
   pp,ll=a64.a64_load('check:\n '+prog_+'\n ret\n');oo={}
   ct=a64.a64_run(pp,ll,'check',dict(x1=av,x2=bv,x30=777),{},out_regs=oo)
   assert ct==2 and oo['x0']==want,(prog_,hex(av),hex(bv));isa+=1
  pp,ll=a64.a64_load('check:\n tst x1,#0xffffe0000fffffff\n cset x0,eq\n ret\n');oo={}
  ct=a64.a64_run(pp,ll,'check',dict(x1=av^bv,x30=777),{},out_regs=oo)
  assert ct==3 and oo['x0']==int((av^bv)&0xffffe0000fffffff==0);isa+=1
  pp,ll=a64.a64_load('check:\n tst x1,#0xfffffffff\n ccmp x2,x3,#2,ne\n cset x0,hs\n ret\n');oo={}
  ct=a64.a64_run(pp,ll,'check',dict(x1=av^bv,x2=av,x3=bv,x30=777),{},out_regs=oo)
  assert ct==4 and oo['x0']==int((av^bv)&((1<<36)-1)==0 or av>=bv);isa+=1
  pp,ll=a64.a64_load('check:\n cmp x1,x2\n cset x0,lo\n ret\n');oo={}
  ct=a64.a64_run(pp,ll,'check',dict(x1=av,x2=bv,x30=777),{},out_regs=oo)
  assert ct==3 and oo['x0']==int(av<bv);isa+=1
for off in(32,96):
 pp,ll=a64.a64_load('check:\n stp xzr,x13,[x2,#%d]\n ret\n'%off);mm={};oo={}
 assert a64.a64_run(pp,ll,'check',dict(x2=0x1000,x13=M64,x30=777),mm,out_regs=oo)==2
 assert mm=={0x1000+off:0,0x1008+off:M64};isa+=1
for zb in[0x100000,0x100100]:
 mm={zb+8*k:5 for k in range(32)};mm[zb-8]=7;mm[zb+256]=9;oo={}
 ct=a64.a64_run(*a64.a64_load('check:\n dc zva, x1\n ret\n'),'check',dict(x1=zb,x30=777),mm,out_regs=oo)
 assert ct==2 and all(mm[zb+8*k]==0 for k in range(32))and mm[zb-8]==7 and mm[zb+256]==9 and oo['x1']==zb;isa+=1
LS=[(0x200000,0x200100)]
def ls(s,x,sr=LS):
 o={}
 try:a64.a64_run(*a64.a64_load('check:\n '+s+'\n ret\n'),'check',dict(x1=0x200000+x,x30=777),{a:a for a in range(0x200000,0x200100,8)},out_regs=o,ls64=(LS,sr));return o
 except AssertionError as e:return str(e)
r=ls('ld64b x4,[x1]',64);assert[r['x%d'%k]for k in range(4,12)]==list(range(0x200040,0x200080,8))and type(ls('st64b x4,[x1]',64))==dict
for s,x,sr in[('ld64b x4,[x1]',32,LS),('ld64b x5,[x1]',64,LS),('ld64b x24,[x1]',64,LS),('ld64b x4,[x1,#64]',0,LS),('ld64b x4,[x1]',256,LS),('st64b x4,[x1]',64,[])]:
 assert 'LS64' in ls(s,x,sr);isa+=1
assert 'DC ZVA' in ls('dc zva, x1',64)
for s,x,sr in[('ld64b x4,[x6]',64,LS),('st64b x4,[x11]',64,LS)]:
 o={}
 try:a64.a64_run(*a64.a64_load('check:\n mov x6,x1\n mov x11,x1\n '+s+'\n ret\n'),'check',dict(x1=0x200000+x,x30=777),{a:a for a in range(0x200000,0x200100,8)},out_regs=o,ls64=(LS,LS));r='ok'
 except AssertionError as e:r=str(e)
 assert 'base in Xt' in r;isa+=1
try:a64.a64_load('a:\n nop\na:\n ret\n');assert 0
except AssertionError as e:assert 'duplicate label' in str(e);isa+=1
for bad in('tst x1,x2,lsl #3','adds w0,w1,w2','ldp w0,w5,[x1]'):
 try:a64.a64_load('k:\n '+bad+'\n ret\n');assert 0
 except AssertionError as e:assert 'not implemented' in str(e);isa+=1
for av in[0,1,M64,1<<63,0xa5a5a5a5a5a5a5a5]:
 for bv in[0,1,M64,1<<63,0x5a5a5a5a5a5a5a5a]:
  for prog_,want in[('lsl x0,x1,#5',(av<<5)&M64),('adds x0,x1,x2,lsl #32',(av+(bv<<32))&M64),('subs x0,x1,x2,lsl #32',(av-(bv<<32))&M64),('and x0,x1,#0xffffffff00000000',av&0xffffffff00000000),('and x0,x1,#0xfffffffe',av&0xfffffffe),('eor x0,x1,#192',av^192),('orr x0,x1,#0x30000000000',av|0x30000000000),('orr x0,x1,#0x600000000000',av|0x600000000000),('and x0,x1,#0x1fffffffff00',av&0x1fffffffff00),('extr x0,x2,x1,#57',((bv<<64|av)>>57)&M64),('sub x0,xzr,x1',(-av)&M64),('lsr x0,x1,#32',av>>32)]:
   pp,ll=a64.a64_load('check:\n '+prog_+'\n ret\n');oo={}
   assert a64.a64_run(pp,ll,'check',dict(x1=av,x2=bv,x30=777),{},out_regs=oo)==2 and oo['x0']==want,prog_;isa+=1
  for c in(0,1):
   pp,ll=a64.a64_load('check:\n cmn x2,#1\n sbcs xzr,x1,xzr\n csel x0,x1,x2,hs\n sbc x3,xzr,xzr\n adc x4,xzr,xzr\n ret\n');oo={}
   a64.a64_run(pp,ll,'check',dict(x1=av,x2=bv,x30=777),{},out_regs=oo);cc=int(bv==M64);cc=int(av-(1-cc)>=0)
   assert oo['x0']==(av if cc else bv)and oo['x3']==(0 if cc else M64)and oo['x4']==cc;isa+=1
assert hashlib.sha256('\n'.join(x for raw in a64.expand((root/'p256_inv.s').read_text()).splitlines()if(x:=raw.split('//')[0].rstrip()).strip()).encode()).hexdigest()=='a41f52a802227f4457c142d85fff53be81c046c107affaaf8f9acffcdb0df925'
RR=1<<256;R=RR%P;RI=pow(R,-1,P);MASK=M64;put=fc.put;get=fc.get
spg,slb=a64.a64_load((root/'p256_fields.s').read_text());sr=random.Random(3901)
edges={0,1,2,P-2,P-1,P,P+1,RR-2,RR-1}
for bit in[32,64,96,128,160,192,224,255]:
 for off in[-1,0,1]:edges.add(((1<<bit)+off)%RR)
pairs=[(a,b)for a in sorted(edges)for b in sorted(edges)]+[(sr.randrange(RR),sr.randrange(RR))for _ in range(512)]
short_counts={};short_alias=[0,0,0];ni=-pow(P,-1,RR)%RR
for case,(a,b)in enumerate(pairs):
 mode=case%3;short_alias[mode]+=1;z=[0x11000,0x22000,0x33000][mode];mem={}
 for addr in[0x11000,0x22000,0x33000]:
  for j in range(-2,6):mem[addr+8*j]=sr.getrandbits(64)
 put(mem,0x11000,a);put(mem,0x22000,b);before=dict(mem)
 rg={f'x{i}':sr.getrandbits(64) for i in range(31)};rg.update(x0=z,x1=0x11000,x2=0x22000,sp=0x44000,x30=0x77777777);o={}
 n=a64.a64_run(spg,slb,'k_mul',rg,mem,out_regs=o,limit=500)
 t=(a*b+(a*b*ni%RR)*P)//RR;assert 0<=t<RR+P
 assert get(mem,z)==(t if t<RR else t-P)and get(mem,z)%P==a*b*RI%P
 assert n==(104 if t<RR else 110)
 assert all(mem[k]==v for k,v in before.items()if not z<=k<z+32)
 assert all(o[f'x{i}']==rg[f'x{i}'] for i in range(18,31)) and o['sp']==rg['sp']
 short_counts[n]=short_counts.get(n,0)+1
assert set(short_counts)=={104,110}
short=dict(cases=len(pairs),counts=short_counts,aliases=short_alias,full_width=True,canaries_and_preserved_registers=True)
src=(root/'p256_fields.s').read_text();inv=(root/'p256_inv.s').read_text()
pair_x19=[ln for ln in fw.bwd_pair(src,'',0) if 'x19' in ln]
assert pair_x19 and all(ln.strip()=='ldadd x19,x15,[x13]' for ln in pair_x19)
assert 'x19' not in fw.FLUSH_SUB
bwd_txt=fw.bwd_listing(src)
assert bwd_txt.count('mov x19, #32')==1 and 'mov x17,#32' not in bwd_txt and 'mov x17, #32' not in bwd_txt
listing=src+inv+fw.fwd_listing(src)+bwd_txt+fw.ychain_f_listing(src)+fw.ychain_b_listing(src)+fw.yprep_listing(src)+fw.tbld_listing()+fw.group_listing()
prog,labels=a64.a64_load(listing)
shapes=set()
for t in(listing,fw.block_production(root),(root/'route.s').read_text(),(root/'phase4.s').read_text()):shapes|={fc.shape(mn,o)for mn,o in a64.a64_load(t)[0]}
assert shapes<=fc.FORMS and fc.shape('tst',['x1','x2','lsl #3'])not in fc.FORMS# every shipped instruction has an implemented shape
fold=fc.run_fold();assert fold['all_checks']
group=fc.run_group(prog,labels,seed)
cases=[fc.run_case(prog,labels,seed*16+k,mp,ny,room)for k,(mp,ny,room)in enumerate([(4,1,None),(4,3,None),(8,2,None),(12,1,None),(4,2,1<<20),(4,2,0)])]
errs=[e for c in cases for e in c['err']];assert-21<=min(errs)and max(errs)<=12
assert[c['flushes']for c in cases][4:]==[(32,0),(0,32)]
yk=fc.run_y(prog,labels,seed,6);tb=fc.run_tbld(prog,labels,seed)
cen=fc.census(seed,32);site={s:sum(v)/len(v)for s,v in cen.items()}
assert cen['B'][0]==1 and cen['F'][0]<0.5 and max(site.values())<0.35
fx=[fc.census(seed*16+k,16,fixed_x=True)for k in range(4)];fxr=max(sum(v)/len(v)for c in fx for v in c.values());assert fxr<0.35
tails=dict(fixed_x_max_site_rate=round(fxr,4),site_rates={s:round(v,4)for s,v in site.items()},forced_final_update_rate=cen['B'][0],first_prefix_rate=cen['F'][0],mean_excess_per_incidence=round(2*sum(site.values()),4),charged=4.2,worst=12,pairs=32*64)
saved_stdin=sys.stdin;sys.stdin=io.StringIO(json.dumps(dict(seed=193)));buf=io.StringIO()
with contextlib.redirect_stdout(buf):runpy.run_path(str(root/'resolution_check.py'),run_name='__main__')
sys.stdin=saved_stdin;resolution=json.loads(buf.getvalue());assert resolution['all_checks']
buf=io.StringIO()
with contextlib.redirect_stdout(buf):runpy.run_path(str(root/'flip_ledger.py'),run_name='__main__')
ledger=json.loads(buf.getvalue())
import fixed_layout
layout=fixed_layout.check(ledger)
gl=fw.group_listing().split('\n');i0=gl.index('fg_round:')+1;i1=next(k for k in range(i0,len(gl))if gl[k].strip().startswith('stp x6, x5, [x9, #56]'))+1
rp,rl=a64.a64_load('check:\n'+'\n'.join(gl[i0:i1])+'\n ret\n');Lr=ledger['L'];M36=(1<<36)-1;rounds=[]
for g,u in[(0,0),(ledger['groups']-1,Lr-1)]:
 rs=g*Lr+u;pb=0x40000;mm={0x30008:pb,pb+160:rs,pb+152:M36};oo={}
 ct=a64.a64_run(rp,rl,'check',dict(sp=0x30000,x30=777),mm,out_regs=oo)
 idy,w2=mm[pb+56],mm[pb+64];yi=((idy>>45)<<28)|(w2>>36)
 assert rs<1<<37 and yi==1024*rs and w2&M36==M36 and idy&((1<<45)-1)==0 and divmod(rs,Lr)==(g,u),(g,u,rs,yi)
 rounds.append((g,u,rs))
# Exact boundary of every compact metadata field: full source coordinate, X-pair511, Y-index2^47-1.
yi=(1<<47)-1;jpair=511;w2=M36|((yi&((1<<28)-1))<<36);w3=M36|(jpair<<36)|((yi>>28)<<45)
assert w2==w3==(1<<64)-1 and w2>>36==yi&((1<<28)-1)
assert (w3>>36)&511==jpair and ((w3>>45)<<28)|(w2>>36)==yi
# manifest literal: the submitted score and every accounting number of the statement equal the ledger's (red line k); one mutant must fail
sub=json.load(open(root/'submission.json'))['claims'][0];num=int(ledger['claim'].split('/')[0]);cd='%d.%015d'%(num//10**15,num%10**15);tl=ledger['tails']
def manifest_ok(st,score):
 return F(int(score['numerator']),int(score['denominator']))==F(num,10**15)and all(x in st for x in['at most'+ledger['E_total_instructions']+' ','below2^'+cd+' ','L=%d.'%ledger['L'],'MX=%d permanent'%ledger['mx'],'G=%d groups of %d Y-chains'%(ledger['groups'],ledger['yb']),'within[%d,%d]'%tuple(ledger['label_interval']),'%d records'%ledger['cap_pg'],'K2 mod2^%d'%ledger['hash_bits'],'K=%d'%ledger['gate_k'],'CAP_H=ceil(m_h)+2^43'])

pf=(root/'proof.md').read_text();d15=lambda fr:'%d.%015d'%(int(fr.split('/')[0])//10**15,int(fr.split('/')[0])%10**15)
assert all(x in pf for x in[ledger['E_total_instructions'],ledger['claim'],'2^'+cd+' units','L=%d'%ledger['L'],d15(ledger['claim_no_ls64']),d15(ledger['gap_no_ls64'])+' above the claim',ledger['run_worst_log2_units_ceil4']])# proof.md figures equal the ledger's
assert manifest_ok(sub['statement'],sub['score'])and not manifest_ok(sub['statement'].replace(ledger['E_total_instructions'],str(int(ledger['E_total_instructions'])+1)),sub['score'])
import source_recovery
recovery=source_recovery.check()
assert recovery['all_checks']
print(json.dumps(dict(research_only=True,source_recovery=recovery,metadata_boundary=True,fold=fold,round_split=rounds,protected_constants=True,manifest_matches_ledger=True,label_error=[min(errs),max(errs)],isa_checks=isa,shapes=len(shapes),short_mul=short,flip_cases=cases,flip_group=group,per_y_kernels=yk,tbld=tb,tails=tails,cold_accounting=cold,routing=routing,resolution=resolution,ledger=ledger,fixed_layout={k:v for k,v in layout.items()if k!='ranges'},all_checks=True,scope='Finite evidence only; the analytical argument is in proof.md.')))

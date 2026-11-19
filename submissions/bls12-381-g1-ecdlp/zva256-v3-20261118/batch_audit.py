"""s26-s28: batch round and both side calls, executed by a64.py."""
import json,math,random,re
from fractions import Fraction
from pathlib import Path
import a64,kernel_check as kc
from smx_audit import MN,C2,FERMAT,T,U,PM,P,RA,N0,RM,CB,ld,st,aff_mul
import rows as RW
UNR,SB=32,RW.SB
F6=['adds x11,x20,x26','adcs x12,x21,x27','adcs x13,x22,x28','adcs x14,x23,x0','adcs x15,x24,x1','adcs x16,x25,x3']
def mn(o,to_u=False,post=None):
 B=[' '.join(l.split()) for l in MN];W=[l.split()[1].split(',')[0] for l in B[:-6] if l[:3]!='ldp']
 assert B[-6:]==F6 and not set(W)&set(T+PM+['x2','x18','x30','sp']) and [l for l in B if '[' in l]==['ldp x17,x29,[x2,#16]','ldp x17,x29,[x2,#32]']
 assert 'x2' not in re.split(r'[ ,]',' '.join(B[B.index('ldp x17,x29,[x2,#32]')+1:])) and (post is None or o==-32)
 B=[l.replace('[x2,#16]','[x2,#%d]'%(o+16)).replace('[x2,#32]','[x2],#%d'%post if post else '[x2,#%d]'%(o+32)) for l in B]
 return ['ldp x17,x29,[x2,#%d]'%o]+(B[:-6]+[l.replace(T[i]+',',U[i]+',',1) for i,l in enumerate(F6)] if to_u else B)
# s46 stack frame: [0] backward end slot, [8] unsigned table, [16] signed table, Fermat slots 2 and 5, t spill 144, A spill 192, saves 288..391
FR,TSP,ASP=512,144,192
UT,DT=0x60000,0x68000
X1,V2=['x19','x20','x21','x22','x23','x24'],['x25','x26','x27','x28','x3','x4']
AV=(['x13','x14','x15','x16','x17','x19'],['x21','x22','x23','x24','x25','x26'])
CV=(['x20']+AV[1],['x20']+AV[0])
def fbody(v):
 """s47 forward element: D rows formed pairwise inside the row product (no reload), c_(i-1) in registers
 (two register-renamed copies alternate), args pointer on the frame. P x1, T x2, slot x3."""
 A,C=AV[v],CV[v];Wd=['x0','x4','x11','x12','x27','x28'];Tl='x29'
 L=['stp x1,x2,[x3,#%d]'%RW.OPT]+RW.ldrow(Wd,'x1',0)+RW.ldrow(C[:6],'x2',96)
 for h in range(3):
  L+=['%s x29,%s,%s'%('subs' if h==0 else 'sbcs',C[2*h],Wd[2*h]),'sbcs x30,%s,%s'%(C[2*h+1],Wd[2*h+1]),'stp x29,x30,[x3,#%d]'%(RW.OE+16*h)]
 L+=['%s %s,%s,%s'%('adds' if i==0 else ('adc' if i==5 else 'adcs'),Wd[i],Wd[i],C[i]) for i in range(6)]+RW.strow(Wd,'x3',0)
 L+=RW.row0(A[5],Wd,C,Tl)
 U=['x29',A[5]]
 for j in range(1,6):
  for h in range(3):
   L+=['ldp %s,%s,[x1,#%d]'%(Wd[2*h],Wd[2*h+1],RW.xoff(j)+16*h),'ldp %s,%s,[x2,#%d]'%(U[0],U[1],RW.voff(j)+16*h),
    '%s %s,%s,%s'%('adds' if h==0 else 'adcs',Wd[2*h],Wd[2*h],U[0]),'%s %s,%s,%s'%('adc' if h==2 else 'adcs',Wd[2*h+1],Wd[2*h+1],U[1])]
  L+=RW.strow(Wd,'x3',48*j)+RW.rowk(A[5-j],Wd,C,Tl,c7='x30' if j>=4 else None,first7=(j==4))
 L+=RW.mont_word(C,A[0],Wd,Tl,c7='x30')+['ldr %s,[sp,#8]'%A[2]]+RW.corr_l(C[1:7],C[6],A[3],A[2],(Wd[0],Wd[1]),via='x30')
 return L+RW.strow(C[1:7],'x3',RW.OC)+['add x3,x3,#%d'%SB,'ldr x30,[sp,#32]','ldp x1,x2,[x30,#32]!','str x30,[sp,#32]']
TBO,N0O,PO,INVO,TMRB=400,408,416,464,1<<40
def bbody():
 """s47 backward element: inv_i is the raw t c_(i-1) representative from a 6x6 product and seven table sums (x18 holds the table base for it,
 p and n0 are reloaded from the frame), then s46's row products, mu combination and t = RP7(t; D)."""
 C=['x17']+T;W=['x0','x1','x3','x4','x25','x26']
 L=['ldr x18,[sp,#%d]'%TBO]+RW.mul_tmr(U,'x2',RW.OC-SB,['x0','x1','x3','x4']+PM+['x11','x12'],['x13','x14','x15','x16','x17','x25'],['x26','x27'],'x28','x29','x30','x18')
 L+=RW.strow(['x7','x8','x9','x10','x11','x12'],'sp',INVO)+RW.ldrow(PM,'sp',PO)+['ldr x18,[sp,#%d]'%N0O]+RW.strow(U,'sp',TSP)+RW.ldrow(U,'sp',INVO)
 L+=['ldr x29,[x2,#%d]'%RW.OPT]+RW.rowprod(U,'x29',RW.yoff,C,W,'x27','x28')+RW.strow(T,'sp',ASP)
 L+=['ldr x29,[x2,#%d]'%(RW.OPT+8)]+RW.rowprod(U,'x29',RW.yoff,C,W,'x27','x28')
 L+=RW.ldrow(U,'sp',ASP)+RW.mucomb_l(U,T,['x0','x1','x3','x4','x25','x26','x27'],'x28','x29','x17','x2','ldr %s,[sp,#8]','ldr %s,[sp,#16]')
 L+=RW.ldrow(T,'sp',TSP)+RW.rowprod(T,'x2',lambda j:48*j,['x17']+U,['x0','x1','x3','x4','x26','x27'],'x28','x29',c7='x25',tab=('x26','ldr x26,[sp,#8]'))
 return L+['sub x2,x2,#%d'%SB]
FBODY,FBODY2,BBODY=fbody(0),fbody(1),bbody()
assert len(FBODY)==len(FBODY2)
FB,BB=len(FBODY),len(BBODY)
FWD,BWD=UNR*FB+1,UNR*BB+3
BATCH=Fraction(FWD+BWD,UNR)
LISTING=(['batch_round:','sub sp,sp,#%d'%FR]+['stp x%d,x%d,[sp,#%d]'%(19+2*i,20+2*i,288+16*i) for i in range(6)]+['str x18,[sp,#384]',
 'mov x4,#%d'%CB,'add x4,x4,#0','ldp x5,x6,[x4]','ldp x7,x8,[x4,#16]','ldp x9,x10,[x4,#32]','ldr x18,[x4,#48]','ldp x13,x14,[x4,#56]','ldp x15,x16,[x4,#72]',
 'ldp x17,x19,[x4,#88]']+RW.strow(AV[0],'x1',RW.OC-SB)+['mov x3,#%d'%UT,'mov x4,#%d'%DT,'stp x3,x4,[sp,#8]']+RW.strow(PM,'sp',PO)+['str x18,[sp,#%d]'%N0O,'mov x3,#%d'%TMRB,'str x3,[sp,#%d]'%TBO,'sub x3,x1,#%d'%SB,'str x3,[sp]','mov x3,x1',
 'ldp x1,x2,[x0,#8]!','str x0,[sp,#32]','.Lfwd:']+UNR//2*(FBODY+FBODY2)+['cbnz x1,.Lfwd','.Lfend:','str x3,[sp,#24]']+RW.ldrow(T,'x3',RW.OC-SB)+C2+st(T,2)+st(T,5)+FERMAT+ld(T,2)
 +['subs x19,x5,x11']+['sbcs %s,%s,%s'%(U[i],PM[i],T[i]) for i in range(1,5)]+['sbc x24,x10,x16','ldr x30,[sp]','ldr x2,[sp,#24]','sub x2,x2,#%d'%SB,'.Lbwd:']
 +UNR*BBODY+['ldr x30,[sp]','cmp x2,x30','b.ne .Lbwd','.Lbend:']
 +['ldp x%d,x%d,[sp,#%d]'%(19+2*i,20+2*i,288+16*i) for i in range(6)]+['ldr x18,[sp,#384]','add sp,sp,#%d'%FR,'ret'])
ROUND_FIXED=len([l for l in LISTING if l[-1]!=':'])-FWD-BWD
def chain_rows(lanes,vrows=None):
 """exact model of the s46 round: lanes (x1, y1, x2, y2) Montgomery residues; returns per lane (mu-, mu+, c_i, inv_i)."""
 c,cs,D=RM,[],[]
 for k,(x1,y1,x2,y2) in enumerate(lanes):
  V=vrows[k] if vrows else [P-r for r in RW.rows_of(x2)]
  D.append([a+b for a,b in zip(RW.rows_of(x1),V)]);cs.append(c);c=RW.corr(RW.rp_raw(c,D[-1]))
 cn=c%P;t=P-RA*RA*pow(cn,-1,P)%P;out=[None]*len(lanes);prods=[]
 for k in range(len(lanes)-1,-1,-1):
  x1,y1,x2,y2=lanes[k];prods.append(t*cs[k]);iv=RW.tmr_raw(t*cs[k])
  assert 0<iv<RA
  A,B=RW.rp_raw(iv,RW.rows_of(y1)),RW.rp_raw(iv,RW.rows_of(y2))
  assert 0<=A<7*P and 0<=B<7*P
  out[k]=RW.mucomb_val(A,B)+(cs[k+1] if k+1<len(cs) else c,iv);t=RW.corr(RW.rp_raw(t,D[k]))
 chain_rows.prods=prods
 return out
def put_point(m,a,x,y,entry=False):
 kc.put(m,a,x);kc.put(m,a+48,y)
 for j in range(1,6):kc.put(m,a+RW.yoff(j),RW.rows_of(y)[j])
 if entry:
  kc.put(m,a+96,P-x);kc.put(m,a+144,P-y)
  for j in range(1,6):kc.put(m,a+RW.voff(j),P-RW.rows_of(x)[j])
 else:
  for j in range(1,6):kc.put(m,a+RW.xoff(j),RW.rows_of(x)[j])
VR_L=RW.vrows_l('x1',['x11','x12','x13','x14','x15','x16'],PM);VROWS=len(VR_L);assert VROWS==60
VRP=a64.a64_load('vr:\n'+'\n'.join('\t'+l for l in VR_L+['ret']))
def tables(m):
 for base,sg in ((UT,False),(DT,True)):
  for idx,v in enumerate(RW.tab_entries(sg)):kc.put(m,base+64*idx,v)
def run(cfg):
 r=random.Random(cfg.get('seed',93));M=(1<<64)-1;nc=0
 prog,lb=a64._a64_parse('\n'.join(l if l[-1]==':' else '\t'+l for l in LISTING))
 lo,mid,hi,end=(lb[k][0] for k in ('.Lfwd','.Lfend','.Lbwd','.Lbend'))
 K=[kc.prepare(kc.kernel_text(pl)) for pl in (False,True)];G=(int(kc.params['gx'],16),int(kc.params['gy'],16))
 AA,AB,BUF,PT,OUT,SP=0x100000,0x140000,0x200000,0x300000,0x800000,0x7000000
 for n in cfg.get('lanes',[2*UNR]):
  assert n%UNR==0 and n>=2*UNR
  m={CB+8*i:v for i,v in enumerate([(P>>(64*i))&M for i in range(6)]+[N0]+[(RM>>(64*i))&M for i in range(6)])}
  for c_ in K:m.update(c_[2])
  tables(m);want=[]
  tb=cfg.get('table');kq,c0=(r.randrange(1,1<<40) if tb else 1),2*(n+2)-1
  for i in range(n):
   if tb:A,B=aff_mul(c0*kq,G),aff_mul((i+1)*kq,G)
   else:
    A=B=aff_mul(1,G)
    while A[0]==B[0]:A,B=(aff_mul(r.randrange(1,1<<40),G) for _ in 'ab')
   q=PT+2048*i+512;t=q+1024;x1,y1,x2,y2=(v*RA%P for v in A+B)
   put_point(m,q,x1,y1);put_point(m,t,x2,y2,True)
   for j,v in enumerate((OUT+2048*i+256,q,t,BUF+SB*i)):m[AB+32*i+8*j]=v;m[AA+32*i+8*j]=v+1024*(j==0)
   want.append((A,B,(x1,y1,x2,y2)))
  for g in (AA,AB):m.update({g+32*n+8*j:0 for j in range(4)})
  chain_rows([w[2] for w in want])
  for pr in chain_rows.prods:RW.tmr_put(m,TMRB,pr,0,kc.put)
  rg={'x%d'%i:r.getrandbits(64) for i in range(31)};rg.update(x0=AB,x1=BUF,sp=SP)
  mem,ro,vis=dict(m),{},{i:0 for i in list(range(lo,mid))+list(range(hi,end))}
  c=a64.a64_run(prog,lb,'batch_round',rg,mem,limit=10**7,out_regs=ro,visits=vis)
  assert (sum(vis[i] for i in range(lo,mid)),sum(vis[i] for i in range(hi,end)),mid-lo,end-hi,c)==(FWD*n//UNR,BWD*n//UNR,FWD,BWD,ROUND_FIXED+BATCH*n)
  assert all(ro[k]==rg[k] for k in ['sp']+['x%d'%i for i in range(18,31)]) and all(mem[a_]==v or BUF-SB<=a_<BUF+SB*n for a_,v in m.items())
  assert all(BUF-SB<=a_<BUF+SB*n or SP-FR<=a_<SP for a_ in set(mem)-set(m))
  model=chain_rows([w[2] for w in want])
  for i,(mm,mp,cc,iv) in enumerate(model):
   b=BUF+SB*i
   assert (kc.get(mem,b+RW.OMM),kc.get(mem,b+RW.OMP),kc.get(mem,b+RW.OC))==(mm,mp,cc) and max(mm,mp)<P+(1<<376)
   assert kc.get(mem,b+RW.OE)==(P-want[i][2][2]-want[i][2][0])%(1<<384)
   x1,y1,x2,y2=want[i][2];assert (iv-RA*RA*pow(x2-x1,-1,P))%P==0
   assert (mm*(x2-x1)-(y1+y2)*RA)%P==0 and (mp*(x2-x1)-(y1-y2)*RA)%P==0
  nc+=sum(iv>=P for _,_,_,iv in model)
  for k,(pg,lb2,_) in enumerate(K):a64.a64_run(pg,lb2,'side_lanes',dict(rg,x0=(AA,AB)[k],x1=n),mem,limit=10**6)
  for i,(A,B,_) in enumerate(want):
   for k in (0,1):
    E=(B[0],B[1] if k else P-B[1]);l=(E[1]-A[1])*pow(E[0]-A[0],-1,P)%P;x3=(l*l-A[0]-E[0])%P
    o=OUT+2048*i+256+1024*(1-k);y=(l*(A[0]-x3)-A[1])*RA%P
    assert (kc.get(mem,o),kc.get(mem,o+48))==(x3*RA%P,y),(n,i,k)
    assert kc.rows_ok(mem,o,x3*RA%P,y) and all(o+j not in mem for j in range(96,192,8))
    if tb:
     T_=aff_mul((c0+(i+1)*(2*k-1))*kq,G);assert (kc.get(mem,o),kc.get(mem,o+48))==(T_[0]*RA%P,T_[1]*RA%P)
     assert a64.a64_run(*VRP,'vr',dict(x1=o,**{PM[q]:P>>(64*q)&M for q in range(6)}),mem)==VROWS+1
     assert all(kc.get(mem,o+RW.voff(j))==P-RW.rows_of(x3*RA%P)[j] for j in range(1,6))
 beta_composed(m=mem,arrays=(AA,AB),sp=SP)
 return {'status':'passed','unreduced_inverses':nc,'raw_inverse_gt_2p':sum(iv>2*P for _,_,_,iv in model),'raw_inverse_max_over_p':max(iv/P for _,_,_,iv in model),'unroll':UNR,'forward_per_unr':FWD,'backward_per_unr':BWD,'batch_per_lane':str(BATCH),'round_fixed':ROUND_FIXED,'fermat':len(FERMAT),
 'rounds_lanes':cfg.get('lanes',[2*UNR]),'table_level':bool(tb),'scope':'s46 executed round (row products, mu combination) + both side calls with rows vs exact P-T, P+T, rows, keys and the exact model'}
LAM=0xd201000000010000**2-1;P5,TAB,QP=P>>320,1<<40,0x90000
AC,TM,PR=['x%d'%i for i in range(19,25)],['x25','x26','x27','x28','x0','x1'],['x%d'%i for i in range(5,11)]
l6=lambda r_,o:['ldp %s,%s,[x2,#%d]'%(r_[2*j],r_[2*j+1],o+16*j) for j in range(3)]
a6=lambda x,y:['adds %s,%s,%s'%(x[0],x[0],y[0])]+['adcs %s,%s,%s'%(x[i],x[i],y[i]) for i in range(1,5)]+['adc %s,%s,%s'%(x[5],x[5],y[5])]
s6=lambda d_,x,y:['subs %s,%s,%s'%(d_[0],x[0],y[0])]+['sbcs %s,%s,%s'%(d_[i],x[i],y[i]) for i in range(1,6)]
WIDTHS=(48,48,48,48,48,48,48,45)
SHIFTS=tuple(48*i for i in range(8))
TABLE_ENTRIES,TABLE_BYTES=8<<48,1<<57  # Initialize all slots including unused top-chunk padding.
QP=TAB-512
assert P<1<<381 and sum(WIDTHS)==381
CH=[]
for sh,w in zip(SHIFTS,WIDTHS):
 j,k=divmod(sh,64)
 if j==5:CH.append(['lsr x3,x16,#%d'%k])
 elif k==0:CH.append(['and x3,x%d,#%s'%(11+j,hex((1<<w)-1))])
 else:CH.append(['extr x3,x%d,x%d,#%d'%(12+j,11+j,k),'and x3,x3,#%s'%hex((1<<w)-1)])
def tabaddr(k,c):return TAB+512*c+64*k
assert tabaddr(7,(1<<48)-1)+64==TAB+TABLE_BYTES
BETA_BODY=['ldp x11,x12,[x4]','ldp x13,x14,[x4,#16]','ldp x15,x16,[x4,#32]']
for i in range(8):BETA_BODY+=CH[i]+['add x2,x17,x3,lsl #9']+(l6(AC,0) if i==0 else l6(TM,64*i)+a6(AC,TM))
RED=['udiv x3,x24,x18','add x2,x17,x3,lsl #6']+l6(TM,-512)+s6(AC,AC,TM)+s6(TM,AC,PR)+['csel %s,%s,%s,lo'%(AC[i],AC[i],TM[i]) for i in range(6)]
BETA_BODY+=(RED+['stp x20,x21,[x4,#104]','str x22,[x4,#120]']+['adds x25,x11,x19']+['adcs %s,x1%d,%s'%(TM[i],i+1,AC[i]) for i in range(1,5)]+['adc x1,x16,x24']
 +s6(AC,TM,PR)+['csel %s,%s,%s,lo'%(TM[i],TM[i],AC[i]) for i in range(4)]+['subs x19,x5,x25','sbcs x20,x6,x26','sbcs x21,x7,x27','sbc x22,x8,x28',
 'stp x20,x21,[x4,#136]','str x22,[x4,#152]','ldr x4,[x30,#32]!','cbnz x4,.Lbeta'])
BETA_LEN=len(BETA_BODY)
BUILD=['stp x19,x20,[x4]','stp x21,x22,[x4,#16]','stp x23,x24,[x4,#32]','add x4,x4,#512']+a6(AC,T)+s6(TM,AC,PR)+['csel %s,%s,%s,lo'%(AC[i],AC[i],TM[i]) for i in range(6)]+['subs x29,x29,#1','b.ne .Lbb']
def beta_listing():
 C=0x80000
 return (['beta_keys:','ldr x4,[x0]','cbz x4,.Lbn','sub sp,sp,#112']+['stp x%d,x%d,[sp,#%d]'%(19+2*i,20+2*i,16*i) for i in range(5)]
 +['stp x29,x30,[sp,#80]','str x18,[sp,#96]','mov x30,x0','mov x2,#%d'%C,'ldp x5,x6,[x2]','ldp x7,x8,[x2,#16]','ldp x9,x10,[x2,#32]','ldp x17,x18,[x2,#48]','ldr x29,[x2,#64]','.Lbeta:']
 +BETA_BODY+['.Lbx:']+['ldp x%d,x%d,[sp,#%d]'%(19+2*i,20+2*i,16*i) for i in range(5)]+['ldp x29,x30,[sp,#80]','ldr x18,[sp,#96]','add sp,sp,#112','.Lbn:','ret'])
def beta_run(cfg):
 r=random.Random(cfg.get('seed',11));M=(1<<64)-1;G=(int(kc.params['gx'],16),int(kc.params['gy'],16));ld_=lambda t:a64.a64_load('\n'.join(l if l[-1]==':' else '\t'+l for l in t))
 Lg=aff_mul(LAM%0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001,G);be=Lg[0]*pow(G[0],-1,P)%P
 assert Lg[1]==G[1] and pow(be,3,P)==1!=be and 8*P<RA and (8<<320)<=P
 xs=[aff_mul(r.randrange(1,1<<40),G)[0]*RA%P for _ in range(cfg.get('points',24))]+[1,2,P-1,P-2,RA%P,(1<<380)+5,(1<<48)-1]+[r.randrange(1,P) for _ in range(17)]
 xs+=[x for sh in SHIFTS[1:] for x in ((1<<sh)-1,1<<sh,(1<<sh)+1) if x<P]
 for g0 in range(0,len(xs),8):
  g=xs[g0:g0+8];n=len(g);m={}
  for i,X in enumerate(g):
   o=0x20000+192*i;m.update({0x10000+32*i+8*j:(o,r.getrandbits(64),r.getrandbits(64),r.getrandbits(64))[j] for j in range(4)});kc.put(m,o,X)
   m.update({o+8*j:r.getrandbits(64) for j in range(6,24)})
  m[0x10000+32*n]=0
  beta_composed(m,(0x10000,),0x7000000,be)
 cm={};[kc.put(cm,QP+64*q,q*P) for q in range(8)]
 rp,rl=ld_(['red:']+RED+['ret']);nr=0
 for S in [S for S in [q*P+e for q in range(9) for e in (0,1,-1,1<<320,-(1<<320),P-1,1-P,(P5+1)<<320)]+[r.randrange(8*P) for _ in range(cfg.get('red',600))] if 0<=S<8*P]:
  rg={'x%d'%i:0 for i in range(31)};rg.update({**{AC[i]:S>>64*i&M for i in range(6)},**{PR[i]:P>>64*i&M for i in range(6)},'x18':P5+1,'x17':TAB,'x29':TAB+(1<<63),'sp':0x7000000})
  ro={};a64.a64_run(rp,rl,'red',rg,dict(cm),limit=1000,out_regs=ro);assert sum(ro[AC[i]]<<(64*i) for i in range(6))==S%P;nr+=1
 bp,bl=ld_(['beta_build:','.Lbb:']+BUILD+['ret']);E=cfg.get('entries',64)
 for k in range(8):
  st_=be*pow(2,SHIFTS[k],P)%P;rg={'x%d'%i:0 for i in range(31)};mem={}
  rg.update({**{PR[i]:P>>64*i&M for i in range(6)},**{T[i]:st_>>64*i&M for i in range(6)},'x4':tabaddr(k,0),'x29':E,'sp':0x7000000})
  assert a64.a64_run(bp,bl,'beta_build',rg,mem,limit=10**6)==len(BUILD)*E+1 and len(mem)==6*E and all(kc.get(mem,tabaddr(k,j))==j*st_%P for j in range(E))
 return {'status':'passed','per_point':BETA_LEN,'slots':len(xs),'curve_points':cfg.get('points',24),'reduction_cases':nr,'build_per_entry':len(BUILD),'build_entries':8*E,
 'table_bytes':TABLE_BYTES,'scope':'table-trade eight-chunk beta, exact canonical Montgomery x input and keys'}
def beta_composed(m,arrays,sp,be=None):
 """Batch and side outputs feed beta at the production stride."""
 M=(1<<64)-1;C=0x80000
 G=(int(kc.params['gx'],16),int(kc.params['gy'],16))
 if be is None:be=aff_mul(LAM,G)[0]*pow(G[0],-1,P)%P
 m.update({C+8*i:P>>(64*i)&M for i in range(6)})
 m.update({C+48:TAB,C+56:P5+1,C+64:TAB+(1<<63)})
 for q in range(8):kc.put(m,QP+64*q,q*P)
 tx='\n'.join(l if l[-1]==':' else '\t'+l for l in beta_listing())
 pg,lb=a64.a64_load(tx)
 for aa in arrays:
  slots=[];ap=aa
  while m[ap]:
   o=m[ap];X=kc.get(m,o);assert 0<=X<P
   slots.append((o,X));ap+=32
   for k,(sh,w) in enumerate(zip(SHIFTS,WIDTHS)):
    c=X>>sh&(1<<w)-1;kc.put(m,tabaddr(k,c),be*c*pow(2,sh,P)%P)
  before=dict(m);rg={'x%d'%i:(i+1)*31 for i in range(31)};rg.update(x0=aa,sp=sp)
  ro={};steps=a64.a64_run(pg,lb,'beta_keys',rg,m,limit=10**6,out_regs=ro)
  assert steps==BETA_LEN*len(slots)+26
  touched={o+k for o,_ in slots for k in (104,112,120,136,144,152)}
  assert all(m[a]==v or a in touched or sp-112<=a<sp for a,v in before.items()) and set(m)-set(before)<=touched|set(range(sp-112,sp,8))
  assert all(ro[k]==rg[k] for k in ['sp']+['x%d'%i for i in range(18,31)])
  for o,X in slots:
   for off,v in ((104,be*X%P),(136,be*be*X%P)):
    assert [m[o+off+8*j] for j in range(3)]==[v>>(64*j)&M for j in (1,2,3)]
 return len(slots)
if __name__=='__main__':print(json.dumps(run(json.loads(input() or '{}'))))
RMAXC,RC_FIX=1<<16,24
ROUND_CTL="""arrive:
 ldr x1, [x22]
 str x1, [x23]
 ret
sweep:
 mov x0, #0
s1:
 ldr x1, [x20], #8
 orr x0, x0, x1
 subs x2, x2, #1
 b.ne s1
 str x0, [x23]
 ret
release:
 ldr x0, [x19]
r1:
 ldr x3, [x21], #8
 str x0, [x3]
 subs x2, x2, #1
 b.ne r1
 ret
check:
 ldr x0, [x19]
 cbnz x0, halt
 add x24, x24, #1
 cmp x24, x25
 b.hs halt
 ret
halt:
 mov x27, #1
 ret"""
def audit_round_ctl():
 """s44: found (1) / fail (2) flags OR-ed up the barrier tree, released down, round cursor advanced; executed on
 toy trees. Leaf at most 10 instructions, internal node with n children 8n + 5: at most 8 + 10 <= RC_FIX per node."""
 import a64
 from collections import defaultdict
 pg,lb=a64.a64_load(ROUND_CTL);done=0
 for leaves,fan,flags,cur in ((37,4,{},0),(37,4,{5:1},7),(37,4,{0:2,36:1},3),(16,4,{},RMAXC-1),(1024,32,{1000:1},RMAXC-2),(5,2,{},RMAXC-2)):
  mem,level,nid,kids,slot=defaultdict(int),list(range(leaves)),leaves,{},{}
  A=lambda n:n<<20
  while len(level)>1:
   nxt=[]
   for i in range(0,len(level),fan):
    kids[nid]=level[i:i+fan]
    for j,c in enumerate(kids[nid]):slot[c]=A(nid)+4096+8*j;mem[A(nid)+32768+8*j]=A(c)+16
    nxt.append(nid);nid+=1
   level=nxt
  root=level[0];slot[root]=A(root)+16;total=0
  for n,f in flags.items():mem[A(n)]=f
  for n in range(leaves):
   st=a64.a64_run(pg,lb,"arrive",{"x22":A(n),"x23":slot[n]},mem);assert st<=3;total+=st
  for n in sorted(kids):
   st=a64.a64_run(pg,lb,"sweep",{"x20":A(n)+4096,"x2":len(kids[n]),"x23":slot[n]},mem);assert st==4*len(kids[n])+3;total+=st
  agg=0
  for f in flags.values():agg|=f
  assert mem[A(root)+16]==agg
  for n in sorted(kids,reverse=True):
   st=a64.a64_run(pg,lb,"release",{"x19":A(n)+16,"x21":A(n)+32768,"x2":len(kids[n])},mem);assert st==4*len(kids[n])+2;total+=st
  for n in range(leaves):
   o={};st=a64.a64_run(pg,lb,"check",{"x19":A(n)+16,"x24":cur,"x25":RMAXC},mem,out_regs=o);assert st<=7;total+=st
   assert (o["x27"]==1)==bool(agg or cur+1>=RMAXC) and (agg or o["x24"]==cur+1)
  assert total<=RC_FIX*(leaves+len(kids));done+=1
 return done
def eps_lit(B=1<<15):
 """s44: sum of 1/g^2 over g <= B whose production character repeats a smaller t's, plus the tail 1/B."""
 r=0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
 m,seen,val,s=(r-1)>>20,{1},{1:1},Fraction(1,B)
 for g in range(2,B+1):
  p=next((q for q in range(2,math.isqrt(g)+1) if g%q==0),g)
  val[g]=v=pow(g,m,r) if p==g else val[p]*val[g//p]%r
  if v in seen:s+=Fraction(1,g*g)
  seen.add(v)
 return s

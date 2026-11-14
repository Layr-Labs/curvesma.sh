import random;random.Random.rr=random.Random.randrange;SP=(34<<40)+0x7000000;rn,en,dc=range,enumerate,dict
SITE_KERNELS="""
.macro ENTER
mov x4,#977
movk x4,#1,lsl #32
str x30,[sp,#504]
str xzr,[sp]
.endm
.macro LEAVE
ldr x30,[sp,#504]
ret
.endm
t_shx:
ENTER
SHX x5,x25,x28,x26,x27
LEAVE
t_fdx:
ENTER
FDX x10,x11,x7,x24,x10,x11,x5,x6,x12,x13
LEAVE
t_lk1:
ENTER
ldr x30,[x1,#56]
ldp x6,x7,[x1,#288]
LK1 0x56000000000000
LEAVE
t_p1:
ENTER
add x29,x3,#256
P1_6WIN_BFI
str x30,[sp,#496]
LEAVE
t_sqinv:
ENTER
SQ13 x8,x9,x10,x11,x13,x14,x15,x16,x5,x6,x7,x12,x17
F5 x13,x14,x15,x16,x12,x13,x14,x17,x12
LEAVE
t_fs:
ENTER
FS x9,x10,x11,x12,x5,x6,x7,x8,x19,x20,x21,x22,x13
LEAVE
t_faddw:
ENTER
FADDW x13,x14,x15,x16,x13,x14,x15,x16,x17,x18,x19,x20,x5
LEAVE
t_kfp:
ENTER
mov x30,x11
BTBASES
BLK6F x28,x29,x30,x1
ROWSG x28,x29,x30,x1
KU0 320,x1,x30,x29,x28
LEAVE
t_dfx:
ENTER
DFX x10,x11,x12,x13,x7,x5
LEAVE
"""
T_KERNELS="""
t_lkp:
ENTER
ldp x6,x7,[x1,#288]
LKP6_PAIR 0x4c000000000000,LD_OB3_2
LEAVE
t_c1b:
ENTER
C1BSPLIT_BDAY x25,x26,x27,x28
LEAVE
t_efgk:
str x30,[sp,#-16]
mov x4,#977
movk x4,#1,lsl #32
EFGK6
ldr x30,[sp,#-16]
ret
"""
COMMON=dc(t_shx=4,t_fdx=6,t_lk1=2+42,t_p1=232+1,t_sqinv=58+5,t_dfx=5,t_fs=7,t_faddw=7,t_kfp=1+3+36+8+45);EXTRA=dc(t_shx=3,t_fdx=3,t_lk1=6,t_sqinv=3,t_dfx=6,t_fs=4,t_faddw=4,t_kfp=3)
def site_checks(E,I,cb,c4):
 p,C0,M64,put,get=E.p,E.C0,E.M64,E.put,E.get;rgi,ars=E.rgi,E.ars;bv,labels=E.a64.a64_load(E.ksrc("leaf")+SITE_KERNELS+T_KERNELS);rprog,rlab=E.kload("route");al,b7={},{};LS=([(E.AXYX0,E.AXY_END),(E.XT0,1<<64),(SP+64,SP+320),(0x200100,0x200140)],[(0x200000,0x200100)])
 def site(entry,rs,mem=None):
  mem={}if mem is None else mem;out={};n=E.arun(bv,labels,entry,dc(rs,sp=SP),mem,out_regs=out,ls64=LS);return out,n,mem[SP]
 def tally(bx,cz):al[bx]=dc(cases=len(cz),predicate_holds=sum(c[0]for c in cz),deviates_iff_predicate=all(c[0]==c[1]for c in cz),model=all(c[2]for c in cz))
 def fixed(entry,f,n,wd,b6):
  d=b7.setdefault(entry,dc(cases=0,fixups=0,word_ok=True,count_ok=True,word_set_ok=True));d["cases"]+=1;d["fixups"]+=f;d["word_ok"]&=wd==C0*f;d["count_ok"]&=n==6+COMMON[entry]+EXTRA[entry]*f;d["word_set_ok"]&=b6 or not f
 fold=lambda v:(v%2**256)+(v>>256)*C0;a6=lambda v,k:(v%2**(64*k))+(v>>256)*C0>=2**(64*k);cz=[]
 for kd in(0,1,2,3):
  for _ in rn(12):
   b3=I.rr(2**40,2**64);hi=(b3*C0)>>64;low=(2**192-hi+I.rr(hi),2**192-hi-1-I.rr(2**20),I.gb(192),I.gb(126)<<64|(2**64-hi+I.rr(hi)))[kd];B=low+(b3<<192);o,n,wd=site("t_shx",ars(B,("x25","x26","x27","x28")));ov=rgi(o,("x5","x25","x26","x27"));ao=B*2**64-b3*p;cz+=[(low+hi>=2**192,ov!=ao,ov==ao%2**256)];fixed("t_shx",(B%2**64)+hi>=2**64,n,wd,B%2**64>2**64-C0)
 tally("row_shift",cz);cz=[]
 for kd in(0,0,1,2):
  for _ in rn(4):
   mem,w,B2=({},I.gb(256),I.rr(p));g0=E.offs6(w);tb=E.XTP[1]
   mem[SP+360]=0
   for k,pos in en((96,104,112,120,352)):mem[SP+pos]=g0[k]
   if kd==0:es=[I.rr(p)for _ in rn(5)]
   else:T=(5*2**256+2**64-1-I.rr(C0)if kd==1 else 4*2**256+(I.gb(160)<<64)+2**64-1-I.rr(C0));es=[p-1-I.rr(100)for _ in rn(4 if kd==1 else 2)];J=T-B2-sum(es);es+=[J]if kd==1 else[J//3,J//3,J-2*(J//3)];assert len(es)==5 and all(0<=e<p for e in es)
   for k,sh,mi,_ in E.W5V:tk=(w>>sh)&mi;put(mem,tb+E.D5O[k]+(tk<<6),I.gb(256));put(mem,tb+E.D5O[k]+(tk<<6)+32,es[k])
   put(mem,SP+64,B2);o,n,wd=site("t_lk1",dc(x1=SP+64),mem);ov=rgi(o,("x20","x21","x22","x23"));S_=sum(es)+B2;t=S_>>256;cz+=[(S_%2**256+t*C0>=2**256,ov!=(S_%2**256+t*C0)%2**256,ov==E.lkfold(S_))];fixed("t_lk1",S_%2**64+t*C0>=2**64,n,wd,S_%2**64>=2**64-t*C0)
 tally("lookup_sum_fold",cz);cz=[]
 for kd in(0,1,2,3):
  for _ in rn(12):
   u=I.rr(1,2**66+1);m=u*C0;low=(2**256-m+I.rr(m),2**256-m-1-I.rr(2**20),I.gb(256),I.gb(127)<<128|(2**128-m+I.rr(m)))[kd];Z=low+(u<<256);rs=ars(low,("x10","x11","x12","x13"));rs.update(x7=u&M64,x24=u>>64);o,n,wd=site("t_fdx",rs);ov=rgi(o,("x10","x11","x12","x13"));t=fold(Z);cz+=[(t>=2**256,ov!=t,ov==t%2**256)];fixed("t_fdx",a6(Z,2),n,wd,t%2**128<u*C0 and t>=2**128)
 tally("fold",cz);cz=[]
 for kd in(0,1,2,3,4):
  for _ in rn(8 if kd<3 else 1):
   while True:
    T=(I.rr(2**64),I.rr(2**65,2**80),I.rr(p),2,(1<<128)+2)[kd];a=E.sqrt_mod_p(T)
    if a is not None:break
   a+=I.choice((0,p))if a<2**256-p else 0;o,n,wd=site("t_sqinv",ars(a,("x8","x9","x10","x11")));ov=rgi(o,("x13","x14","x15","x16"));bz=rgi(o,("x8","x9","x10","x11"))==a;V=fold(a*a);t=fold(V);cz+=[(t>=p,ov!=T,ov==t%2**256 and bz)];fixed("t_sqinv",a6(V,2),n,wd,t%2**128<2**65)
 tally("square_fold_inv2",cz);cz=[]
 d0=[]
 for kd in(0,1,2,3):
  for _ in rn(10):
   u=I.rr(1,2**64);m=u*C0;low=(2**256-m+I.rr(m),2**128-m+I.rr(m),I.gb(256),I.gb(127)<<128|(2**128-m-1-I.rr(2**20)))[kd];Z=low+(u<<256);rs=ars(low,("x10","x11","x12","x13"));rs.update(x7=u);o,n,wd=site("t_dfx",rs);ov=rgi(o,("x10","x11","x12","x13"));d0+=[ov==E.lkfold(Z)and ov%p==Z%p];fixed("t_dfx",a6(Z,2),n,wd,ov%2**128<(u+1)*C0)
 al["df_exact"]=dc(cases=len(d0),ok=all(d0));p1_ok=True;a4=0
 for _ in rn(10):
  X=E.PROD.mul(I.rr(1,E.r),E.GP);Y=E.PROD.mul(I.rr(1,E.r),E.GP)
  if X[0]==Y[0]:continue
  a,yx=X;b,yy=Y;am=pow((yy*yy-yx*yx)%p,p-2,p);b2=fold(fold(am*am))%2**256;mem={};X3=0x200000;E.popax(rprog,rlab,mem,0,X,b2)
  for k,b_ in en(E.ax_bases(0)+(0x5A5A5A5A5A5A5A5A,)*2):mem[X3+256+8*k]=b_
  for k in rn(4):put(mem,SP+128+32*k,yy*pow(2,64*k,p)%p)
  o,n,wd=site("t_p1",dc(ars(am,("x8","x9","x10","x11")),x3=X3),mem);av=o["x28"]|(o["x29"]<<64)|(mem[SP+496]<<128)|(o["x1"]<<192);au=rgi(o,("x0","x2","x5","x6"));b5=fold(fold(am*am))%2**256;bc=E.lkfold(am+E.lookup6(2*yx*yx%p,b5));c3=E.lookup6(2*yx%p,b5);ba=fold(sum(((c3>>64*k)&M64)*(yy*pow(2,64*k,p)%p)for k in rn(4)))%2**256;p1_ok&=(av==E.fs_op(bc,ba)and au==E.faddw_op(bc,ba)and av%p==pow(yy+yx,-2,p)and au%p==pow(yy-yx,-2,p))
  if wd==0:p1_ok&=(n==6+COMMON["t_p1"]);a4+=1
 al["p1_direct_u2v2"]=dc(ok=p1_ok and a4>0,clean_runs=a4);ah,ai=[],[];pairs=[(I.rr(p),I.rr(p))for _ in rn(16)]
 for _ in rn(8):
  a=I.rr(p>>1);pairs+=[(a,a+(I.rr(1,2**190)<<64)-I.rr(C0))];b=I.rr(p>>1,p);pairs+=[(2**256-b+(I.rr(1,2**120)<<64)+2**64-1-I.rr(C0),b)];az=p+1+I.rr(C0-1);bm=I.rr(az-p);pairs+=[(bm,az)];pairs+=[(p+C0-1-I.rr(10),p+10+I.rr(C0-10))]
 for a,b in pairs:
  if a>=2**256 or b>=2**256:continue
  rs=ars(a,("x5","x6","x7","x8"));rs.update(ars(b,("x19","x20","x21","x22")));o,n,wd=site("t_fs",rs);ov=rgi(o,("x9","x10","x11","x12"));a8=a-b+p<0;aj=a-b if a>=b else a-b+p;ah+=[(a8,ov!=aj%2**256 or a8,(ov==aj)if not a8 else((a-b)%p>=p-C0+1))]
  if not a8:fixed("t_fs",a<b and(a-b)%2**64<C0,n,wd,ov%2**64>=2**64-C0)
  rs=ars(a,("x13","x14","x15","x16"));rs.update(ars(b,("x17","x18","x19","x20")));o,n,wd=site("t_faddw",rs);ov=rgi(o,("x13","x14","x15","x16"));a9=a+b>=2**256+p;ak=a+b if a+b<2**256 else a+b-p;ai+=[(a9,ov!=ak%2**256 or a9,(ov==ak)if not a9 else(C0<=(a+b)%p<=2*C0-2))]
  if not a9:fixed("t_faddw",a+b>=2**256 and(a+b)%2**64+C0>=2**64,n,wd,ov%2**64<C0)
 tally("fs_p1",ah);tally("faddw_p1",ai);cz=[]
 kq=True
 for kd in(0,1,2,3):
  for q_ in rn(6):
   w=I.gb(256)if kd<3 else I.rr(C0)+p
   if kd:w=(w&~(M64<<64))|(M64<<64)|((2**64-1-I.rr(2**34))if kd==1 else w&M64)
   cl="EF"[q_%2];Cb=(E.CE.C_E,E.CE.C_FG)[q_%2];c1=(-Cb,Cb,I.rr(-Cb,Cb+1))[q_%3];c2=(0,2**128-1,I.rr(2**128))[(q_//2)%3];mem={};E.popbt(None,None,mem,w);put(mem,0x200100,c1%2**192,3);put(mem,0x200118,c2,2);put(mem,0x200128,0,3)
   o,n,wd=site("t_kfp",dc(**ars(w,("x28","x29","x11","x1")),x3=0x200100-320),mem);A=o["x17"]<<64|o["x16"];B=o["x19"]<<64|o["x18"];y1=rgi(o,("x24","x25","x26","x27"));f=w%2**64+((w>>192)*C0>>64)>=2**64 and(w>>64)%2**64==M64;cv=(c1+c2*E.BETA)%p
   A_,B_=E.ktfp(c1,c2,w,E.blkf(w));e1=E.werr(E.t128(cv*w%p),A);e2=E.werr(E.t128(E.BETA*cv*w%p),B);kt=E.CE.KT[cl];d1=(((w*E.BETA%p)<<256)-p*y1)%(p<<256);ok=(A,B,y1)==(A_,B_,E.blkf(w))and(d1<=3*p or d1>=(p<<256)-3*p)and all(lo<=d<=hi for ev,r in((e1,kt[0]),(e2,kt[1]))for d,(lo,hi)in zip(ev,r));kq&=ok
   cz+=[(f,E.rowsf(w)!=((w%2**128+((w>>192)*C0>>64))%2**128|w>>128<<128),ok)];fixed("t_kfp",f,n,wd,True)
 tally("rows_kfp",cz);ab,ad,aq={},E.BT0,0x200000
 for idx,wi in en(E.SLOTW):put(ab,aq+64*idx,E.BETA*pow(2,43*wi,p)%p);put(ab,aq+64*idx+32,0)
 nb=E.arun(cb,c4,"tbldf",dc(x0=ad,x1=aq,x2=8,x3=7,x26=256),ab);bs=(nb==E.CE.tbldf_count(8,7)and all(get(ab,ad+256*t+32*k)==E.frac_e(E.BETA*t*pow(2,43*E.SLOTW[k],p)%p)for k in rn(7)for t in rn(8)))
 def tb(ba,c1,a1,c2,a2,W,s=0):
  m,n={},0;o=lambda k,c,a:(a if k==0 else 0)+(c*s*pow(2,215,p)if k==4 else 0)
  for k,sh,R in W:dk=0x200000+128*k;put(m,dk,c1*pow(2,sh,p)%p);put(m,dk+32,o(k,c1,a1)%p);put(m,dk+64,c2*pow(2,sh,p)%p);put(m,dk+96,o(k,c2,a2)%p);n+=E.arun(rprog,rlab,"tbld",dc(x0=ba(k),x1=dk,x2=8,x3=2,x4=2,x5=64,x6=3*(R<8)),m)
  return n==26*len(W)+40*sum(R for _,_,R in W)and all(get(m,ba(k)+64*t)==(c1*t*pow(2,sh,p)+o(k,c1,a1))%p and get(m,ba(k)+64*t+32)==(c2*t*pow(2,sh,p)+o(k,c2,a2))%p for k,sh,R in W for t in rn(R))and all((ba(k)+64*t+o_)not in m for k,_,R in W for t in rn(R,8)for o_ in(0,32))
 bn=tb(lambda k:E.ax_bases(0)[k],I.rr(1,p),0,I.rr(1,p),0,[(k,sh,2 if k in(2,5)else 8)for k,sh in en(E.W6S)]);bp=tb(lambda k:E.XTP[0]+E.D5O[k],I.rr(1,p),I.rr(1,p),I.rr(1,p),I.rr(1,p),[(k,43*k,8)for k in rn(5)],I.gb(41))
 bi=(bs and bn and bp and kq and al["p1_direct_u2v2"]["ok"]and al["df_exact"]["ok"]and all((d["predicate_holds"]>0 or k=="lookup_sum_fold")and d["deviates_iff_predicate"]and d["model"]for k,d in al.items()if"model"in d)and all(d["fixups"]>0 and d["word_ok"]and d["count_ok"]and d["word_set_ok"]for d in b7.values()));return dc(all_true=bi,ktfp_ok=kq,tbld_unified=bs and bn and bp,p1_ok=al["p1_direct_u2v2"]["ok"],sites=len(al),fixup_entries=len(b7))

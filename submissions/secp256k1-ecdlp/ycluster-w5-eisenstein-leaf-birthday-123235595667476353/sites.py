SP,FXA=(34<<40)+0x7000000,(34<<40)+0x10000
SITE_KERNELS="""
.macro ENTER
mov x4,#977
movk x4,#1,lsl #32
str x30,[sp,#504]
mov x30,#0x220000000000
movk x30,#1,lsl #16
str xzr,[x30]
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
ldr x6,[x1,#288]
LK1 0xf000000000000000
LEAVE
t_p1:
ENTER
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
t_dfx:
ENTER
DFX x10,x11,x12,x13,x7,x5
LEAVE
t_kfp:
ENTER
mov x30,x11
BTBASES
BLK5F x28,x29,x30,x1
ROWSG x28,x29,x30,x1
KU0 320,x1,x30,x29,x28
LEAVE
"""
T_KERNELS="""
t_lkp:
ENTER
ldr x6,[x1,#288]
LKP6_PAIR 0xe400000000000000,LD_OB3_2
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
COMMON=dict(t_shx=4,t_fdx=6,t_lk1=2+42,t_p1=220+1,t_sqinv=58+5,t_dfx=5,t_fs=7,t_faddw=7,t_kfp=1+3+30+8+45);EXTRA=dict(t_shx=5,t_fdx=5,t_lk1=8,t_sqinv=5,t_dfx=8,t_fs=6,t_faddw=6,t_kfp=5)
def site_checks(E,J,aw,ax):
 p,C0,M64,put,get=E.p,E.C0,E.M64,E.put,E.get;rgi,ars,R4=E.rgi,E.ars,E.R4;au,labels=E.a64.load(E.kf("leaf.s")+SITE_KERNELS+T_KERNELS);rprog,rlab=E.a64.load(E.kf("route.s"));ac,ag={},{};LS=([(E.AXYX0,E.XT_END),(E.BT0,E.BT0_END),(E.BT1,1<<64),(SP+64,SP+320),(0x200100,0x200140)],[(0x200000,0x200100)])
 def site(entry,regs,mem=None):mem={}if mem is None else mem;out={};n=E.a64.run(au,labels,entry,dict(regs,sp=SP,x30=E.SNT),mem,out_regs=out,ls64=LS);return out,n,mem[FXA]
 def tally(name,cases):ac[name]=dict(cases=len(cases),predicate_holds=sum(c[0]for c in cases),deviates_iff_predicate=all(c[0]==c[1]for c in cases),model=all(c[2]for c in cases))
 def fixed(entry,f,n,word,av):d=ag.setdefault(entry,dict(cases=0,fixups=0,word_ok=True,count_ok=True,word_set_ok=True));d["cases"]+=1;d["fixups"]+=f;d["word_ok"]&=word==C0*f;d["count_ok"]&=n==8+COMMON[entry]+EXTRA[entry]*f;d["word_set_ok"]&=av or not f
 fold=lambda v:(v%2**256)+(v>>256)*C0;b7=lambda v,k:(v%2**(64*k))+(v>>256)*C0>=2**(64*k);cases=[]
 for kind in(0,1,2,3):
  for _ in range(12):b3=J.rr(2**40,2**64);hi=(b3*C0)>>64;aa=(2**192-hi+J.rr(hi),2**192-hi-1-J.rr(2**20),J.gb(192),J.gb(126)<<64|(2**64-hi+J.rr(hi)))[kind];B=aa+(b3<<192);o,n,word=site("t_shx",ars(B,R4(25)));ov=rgi(o,("x5","x25","x26","x27"));b5=B*2**64-b3*p;cases.append((aa+hi>=2**192,ov!=b5,ov==b5%2**256));fixed("t_shx",(B%2**64)+hi>=2**64,n,word,B%2**64>2**64-C0)
 tally("row_shift",cases);cases=[]
 for kind in(0,0,1,2):
  for _ in range(4):
   mem,w,B2=({},J.gb(256),J.rr(p));g0=E.offs6(w);tb=E.XTP[1]
   for k,pos in enumerate((96,104,112,120,352)):mem[SP+pos]=g0[k]
   if kind==0:es=[J.rr(p)for _ in range(5)]
   else:T=(5*2**256-1 if kind==1 else 3*2**256+(J.gb(160)<<64)+2**64-1-J.rr(C0));es=[p-1-J.rr(100)for _ in range(3 if kind==1 else 2)];a5=T-B2-sum(es);es+=[a5//2,a5-a5//2]if kind==1 else[a5//3,a5//3,a5-2*(a5//3)];assert len(es)==5 and all(0<=e<p for e in es)
   for k in range(5):tk=(w>>E.W6S[k])&E.W6M[k];put(mem,tb+E.D6O[k]+(tk<<6),J.gb(256));put(mem,tb+E.D6O[k]+(tk<<6)+32,es[k])
   put(mem,SP+64,B2);o,n,word=site("t_lk1",dict(x1=SP+64),mem);ov=rgi(o,R4(20));S_=sum(es)+B2;t=S_>>256;cases.append((S_%2**256+t*C0>=2**256,ov!=(S_%2**256+t*C0)%2**256,ov==E.lkfold(S_)));fixed("t_lk1",S_%2**64+t*C0>=2**64,n,word,S_%2**64>=2**64-t*C0)
 tally("lookup_sum_fold",cases);cases=[]
 for kind in(0,1,2,3):
  for _ in range(12):u=J.rr(1,2**66+1);m=u*C0;aa=(2**256-m+J.rr(m),2**256-m-1-J.rr(2**20),J.gb(256),J.gb(127)<<128|(2**128-m+J.rr(m)))[kind];Z=aa+(u<<256);regs=ars(aa,R4(10));regs.update(x7=u&M64,x24=u>>64);o,n,word=site("t_fdx",regs);ov=rgi(o,R4(10));t=fold(Z);cases.append((t>=2**256,ov!=t,ov==t%2**256));fixed("t_fdx",b7(Z,2),n,word,t%2**128<u*C0 and t>=2**128)
 tally("fold",cases);cases=[]
 for kind in(0,1,2,3,4):
  for _ in range(8 if kind<3 else 1):
   while True:
    T=(J.rr(2**64),J.rr(2**65,2**80),J.rr(p),2,(1<<128)+2)[kind];a=E.sqrt_mod_p(T)
    if a is not None:break
   a+=J.choice((0,p))if a<2**256-p else 0;o,n,word=site("t_sqinv",ars(a,R4(8)));ov=rgi(o,R4(13));kept=rgi(o,R4(8))==a;V=E.sqv(a);t=fold(V);cases.append((t>=p or(t-T)%p!=0,ov!=T,ov==t%2**256 and kept and(t-T)%p in{0}|{(-2**k)%p for k in(128,192,256)}));fixed("t_sqinv",b7(V,2),n,word,t%2**128<2**65)
 tally("square_fold_inv2",cases);d0=[]
 for kind in(0,1,2,3):
  for _ in range(10):u=J.rr(1,2**64);m=u*C0;aa=(2**256-m+J.rr(m),2**128-m+J.rr(m),J.gb(256),J.gb(127)<<128|(2**128-m-1-J.rr(2**20)))[kind];Z=aa+(u<<256);o,n,word=site("t_dfx",dict(ars(aa,R4(10)),x7=u));ov=rgi(o,R4(10));d0.append(ov==E.lkfold(Z)and ov%p==Z%p);fixed("t_dfx",b7(Z,2),n,word,ov%2**128<(u+1)*C0)
 ac["df_exact"]=dict(cases=len(d0),ok=all(d0))
 p1_ok=True;b8=0
 for _ in range(10):
  X=E.PROD.mul(J.rr(1,E.r),E.GP);Y=E.PROD.mul(J.rr(1,E.r),E.GP)
  if X[0]==Y[0]:continue
  a,yx=X;b,yy=Y;ad=pow((yy*yy-yx*yx)%p,p-2,p);f3=E.sqw(ad);mem={};X3=0x200000;E.popax(rprog,rlab,mem,0,X,f3)
  for k,b_ in enumerate(E.ax_bases(0)+(0x5A5A5A5A5A5A5A5A,)*3):mem[X3+256+8*k]=b_
  for k in range(4):put(mem,SP+128+32*k,yy*pow(2,64*k,p)%p)
  o,n,word=site("t_p1",dict(ars(ad,R4(8)),x3=X3),mem);d3=o["x28"]|(o["x29"]<<64)|(mem[SP+496]<<128)|(o["x1"]<<192);d5=rgi(o,("x0","x2","x5","x6"));d6=E.sqw(ad);d7=E.lkfold(ad+E.msum(2*yx*yx%p,d6));f4=E.lkfold(E.msum(2*yx%p,d6));d8=E.fold(sum(((f4>>64*k)&M64)*(yy*pow(2,64*k,p)%p)for k in range(4)));p1_ok&=(d3==E.fs_op(d7,d8)and d5==E.faddw_op(d7,d8)and d3%p==pow(yy+yx,-2,p)and d5%p==pow(yy-yx,-2,p))
  if word==0:p1_ok&=(n==8+COMMON["t_p1"]);b8+=1
 ac["p1_direct_u2v2"]=dict(ok=p1_ok and b8>0,clean_runs=b8);d9,e0=[],[];pairs=[(J.rr(p),J.rr(p))for _ in range(16)]
 for _ in range(8):a=J.rr(p>>1);pairs.append((a,a+(J.rr(1,2**190)<<64)-J.rr(C0)));b=J.rr(p>>1,p);pairs.append((2**256-b+(J.rr(1,2**120)<<64)+2**64-1-J.rr(C0),b));e1=p+1+J.rr(C0-1);f5=J.rr(e1-p);pairs.append((f5,e1));pairs.append((p+C0-1-J.rr(10),p+10+J.rr(C0-10)))
 for a,b in pairs:
  if a>=2**256 or b>=2**256:continue
  o,n,word=site("t_fs",dict(ars(a,R4(5)),**ars(b,R4(19))));ov=rgi(o,R4(9));b0=a-b+p<0;e2=a-b+p*(a<b);d9.append((b0,ov!=e2%2**256 or b0,(ov==e2)if not b0 else((a-b)%p>=p-C0+1)))
  if not b0:fixed("t_fs",a<b and(a-b)%2**64<C0,n,word,ov%2**64>=2**64-C0)
  o,n,word=site("t_faddw",dict(ars(a,R4(13)),**ars(b,R4(17))));ov=rgi(o,R4(13));b1=a+b>=2**256+p;e3=a+b-p*(a+b>=2**256);e0.append((b1,ov!=e3%2**256 or b1,(ov==e3)if not b1 else(C0<=(a+b)%p<=2*C0-2)))
  if not b1:fixed("t_faddw",a+b>=2**256 and(a+b)%2**64+C0>=2**64,n,word,ov%2**64<C0)
 tally("fs_p1",d9);tally("faddw_p1",e0);cases=[]
 b9=True
 for kind in(0,1,2,3):
  for q_ in range(6):
   w=J.gb(256)if kind<3 else J.rr(C0)+p
   if kind:w=(w&~(M64<<64))|(M64<<64)|((2**64-1-J.rr(2**34))if kind==1 else w&M64)
   cl=("E","FG")[q_%2];Cb=5*E.CE.EIS[cl]["C"]+(E.CE.EIS["EB"]["C"]if cl=="E"else 0);c1=(-Cb,Cb,J.rr(-Cb,Cb+1))[q_%3];c2=(0,2**128-1,J.rr(2**128))[(q_//2)%3];mem={};E.popbt(None,None,mem,w);put(mem,0x200100,c1%2**192,3);put(mem,0x200118,c2,2);put(mem,0x200128,0,3)
   o,n,word=site("t_kfp",dict(**ars(w,("x28","x29","x11","x1")),x3=0x200100-320),mem);A=o["x17"]<<64|o["x16"];B=o["x19"]<<64|o["x18"];y1=rgi(o,R4(24));f=w%2**64+((w>>192)*C0>>64)>=2**64 and(w>>64)%2**64==M64;cv=(c1+c2*E.BETA)%p
   A_,B_=E.ktfp(c1,c2,w,E.blkf(w));e1=E.werr(E.t128(cv*w%p),A);e2=E.werr(E.t128(E.BETA*cv*w%p),B);kt=E.CE.KT[cl[0]];ok=(A,B,y1)==(A_,B_,E.blkf(w))and-2.5<=((w*E.BETA%p)<<256)/p-y1<=2.5 and all(lo<=d<=hi for ev,r in((e1,kt[0]),(e2,kt[1]))for d,(lo,hi)in zip(ev,r));b9&=ok
   cases.append((f,E.rowsf(w)!=((w%2**128+((w>>192)*C0>>64))%2**128|w>>128<<128),ok));fixed("t_kfp",f,n,word,True)
 tally("rows_kfp",cases);a4,a6={},0x200000;put(a4,a6,E.BETA,4);put(a4,a6+32,0,4);put(a4,a6+64,0,4);put(a4,a6+96,0,4);nb0=E.a64.run(aw,ax,"tbldf",dict(x0=E.BT0,x1=a6,x2=8,x3=2,x26=64,x30=E.SNT),a4)
 for idx in range(4):put(a4,a6+64*idx,E.BETA*pow(2,E.W6S[idx+1],p)%p);put(a4,a6+64*idx+32,0)
 nb1=E.a64.run(aw,ax,"tbldf",dict(x0=E.BT1,x1=a6,x2=8,x3=4,x26=128,x30=E.SNT),a4);ah=(nb0==E.CE.tbldf_count(8,2)and nb1==E.CE.tbldf_count(8,4)and all(get(a4,E.BT0+64*t)==E.frac_e(E.BETA*t%p)and get(a4,E.BT0+64*t+32)==0 for t in range(8))and all(get(a4,E.BT1+128*t+32*k)==E.frac_e(E.BETA*t*pow(2,E.W6S[k+1],p)%p)for k in range(4)for t in range(8)))
 e9=(ah and b9 and ac["p1_direct_u2v2"]["ok"]and ac["df_exact"]["ok"]and all(d["predicate_holds"]>0 and d["deviates_iff_predicate"]and d["model"]for d in ac.values()if"model"in d)and all(d["fixups"]>0 and d["word_ok"]and d["count_ok"]and d["word_set_ok"]for d in ag.values()));return dict(all_true=e9,ktfp_ok=b9,tbld_unified=ah,p1_ok=ac["p1_direct_u2v2"]["ok"],sites=len(ac),fixup_entries=len(ag))

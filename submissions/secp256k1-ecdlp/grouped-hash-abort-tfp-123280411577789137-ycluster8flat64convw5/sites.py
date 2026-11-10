import os,random;random.Random.rr=random.Random.randrange;SP,FXA=(34<<40)+0x7000000,(34<<40)+0x10000
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
t_lke:
ENTER
ldp x5,x6,[sp,#384]
ldp x7,x26,[sp,#400]
ldr x27,[sp,#416]
LK6P_CORE 312
add x30,sp,#256
ld64b x8,[x30]
LKP_ADD5 x16,x17,x18,x19,x24,x8,x9,x10,x11
LKP_ADD5 x20,x21,x22,x23,x25,x8,x9,x10,x11
LKFOLD4 x16,x17,x18,x19,x24,x8
LKFOLD4 x20,x21,x22,x23,x25,x8
st64b x16,[x2]
LEAVE
t_lkf:
ENTER
ldp x5,x6,[sp,#96]
ldp x7,x26,[sp,#112]
ldr x27,[sp,#352]
FGPAIR6_ONE 272,NOSTEP
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
t_blk:
ENTER
mov x30,x11
PAIRB6
LEAVE
t_f5w:
ENTER
F5 x10,x11,x12,x13,x7,x10,x11,x5,x6
LEAVE
"""
T_KERNELS="""
t_lkp:
ENTER
ldr x6,[x1,#288]
LKP6_PAIR 416,LD_OB3_2
LEAVE
t_c1b:
ENTER
C1BSPLIT_BDAY x25,x26,x27,x28
LEAVE
t_efg:
mov x4,#977
movk x4,#1,lsl #32
mov x29,x30
str xzr,[sp]
EFG_QUAD_2
mov x30,x29
ret
t_tfp:
ENTER
ldr x30,[sp,#496]
TFP4S_7
LEAVE
"""
COMMON=dict(t_shx=4,t_fdx=6,t_lke=4+69,t_lkf=3+58,t_p1=221+1,t_sqinv=58+5,t_f5w=5,t_fs=7,t_faddw=7,t_blk=1+3+37+37,t_tfp=97+1);EXTRA=dict(t_shx=5,t_fdx=5,t_lke=8,t_lkf=8,t_sqinv=5,t_f5w=5,t_fs=6,t_faddw=6,t_blk=8)
def site_checks(E,J,aw,ax):
 p,C0,M64,put,get=E.p,E.C0,E.M64,E.put,E.get;rgi,ars=E.rgi,E.ars;au,labels=E.a64.a64_load(open(os.path.join(E.HERE,"kernels","leaf.s")).read()+SITE_KERNELS+T_KERNELS);rprog,rlab=E.a64.a64_load(open(os.path.join(E.HERE,"kernels","route.s")).read());ac,ag={},{};LS=([(E.AXYX0,E.XT_END),(E.BT0,E.BT0_END),(E.BT1,1<<64),(SP+128,SP+320)],[(0x200000,0x200100)])
 def site(entry,regs,mem=None):mem={}if mem is None else mem;out={};n=E.a64.a64_run(au,labels,entry,dict(regs,sp=SP,x30=E.SNT),mem,out_regs=out,ls64=LS);return out,n,mem[FXA]
 def tally(name,cases):ac[name]=dict(cases=len(cases),predicate_holds=sum(c[0]for c in cases),deviates_iff_predicate=all(c[0]==c[1]for c in cases),model=all(c[2]for c in cases))
 def fixed(entry,f,n,word,av):d=ag.setdefault(entry,dict(cases=0,fixups=0,word_ok=True,count_ok=True,word_set_ok=True));d["cases"]+=1;d["fixups"]+=f;d["word_ok"]&=word==C0*f;d["count_ok"]&=n==8+COMMON[entry]+EXTRA[entry]*f;d["word_set_ok"]&=av or not f
 fold=lambda v:(v%2**256)+(v>>256)*C0;b7=lambda v,k:(v%2**(64*k))+(v>>256)*C0>=2**(64*k);cases=[]
 for kind in(0,1,2,3):
  for _ in range(12):b3=J.rr(2**40,2**64);hi=(b3*C0)>>64;aa=(2**192-hi+J.rr(hi),2**192-hi-1-J.rr(2**20),J.getrandbits(192),J.getrandbits(126)<<64|(2**64-hi+J.rr(hi)))[kind];B=aa+(b3<<192);o,n,word=site("t_shx",ars(B,("x25","x26","x27","x28")));ov=rgi(o,("x5","x25","x26","x27"));b5=B*2**64-b3*p;cases.append((aa+hi>=2**192,ov!=b5,ov==b5%2**256));fixed("t_shx",(B%2**64)+hi>=2**64,n,word,B%2**64>2**64-C0)
 tally("row_shift",cases);cases=[]
 for entry,f0,tbo,ci in(("t_lke",(384,392,400,408,416),312,0),("t_lkf",(96,104,112,120,352),272,1)):
  for kind in(0,0,1,2):
   for _ in range(4):
    a1,a2=J.rr(1,p),J.rr(1,p);tb,X2=E.XT0+(ci*3*(1<<58)),0x200040;mem={};mem[X2+tbo]=tb;w=J.getrandbits(256);g0=E.offs6(w)
    for k,pos in enumerate(f0):mem[SP+pos]=g0[k]
    b4=J.getrandbits(256)if entry=="t_lke"else 0;g3,g4=E.xcs(a1)[ci],(0 if ci==0 else E.xcs(a1)[3]);g5,g6=E.xcs(a2)[ci],(0 if ci==0 else E.xcs(a2)[3]);E.pop60(rprog,rlab,mem,tb,g3,g4,g5,g6,w)
    if kind:
     b4=(p-1-J.rr(100)if kind==1 else J.rr(2**200))if entry=="t_lke"else 0;T=((5 if entry=="t_lke"else 4)*2**256-1 if kind==1 else 3*2**256+(J.getrandbits(160)<<64)+2**64-1-J.rr(C0));es=[p-1-J.rr(100)for _ in range(3 if kind==1 else 2)];a5=T-b4-sum(es);es+=[a5//2,a5-a5//2]if kind==1 else[a5//3,a5//3,a5-2*(a5//3)];assert len(es)==5 and all(0<=e<p for e in es);es2=[J.rr(p//8)for _ in range(5)]
     for k in range(5):tk=(w>>E.W6S[k])&E.W6M[k];put(mem,tb+E.D6O[k]+(tk<<6),es[k]);put(mem,tb+E.D6O[k]+(tk<<6)+32,es2[k])
    put(mem,SP+256,b4);put(mem,SP+288,0);o,n,word=site(entry,dict(x2=X2),mem);af,ae=get(mem,X2),get(mem,X2+32);a9=sum(get(mem,tb+E.D6O[k]+(((w>>E.W6S[k])&E.W6M[k])<<6))for k in range(5))+b4;b6=sum(get(mem,tb+E.D6O[k]+(((w>>E.W6S[k])&E.W6M[k])<<6)+32)for k in range(5))+b4;t1,t2=a9>>256,b6>>256;f1=a9%2**64+t1*C0>=2**64;f2=b6%2**64+t2*C0>=2**64;g1=(E.xcs(a1)[ci]*w+(b4 if ci==0 else E.xcs(a1)[3]))%p;g2=(E.xcs(a2)[ci]*w+(b4 if ci==0 else E.xcs(a2)[3]))%p;cases.append((fold(a9)>=2**256,af!=fold(a9)%2**256,af==E.lkfold(a9)and ae==E.lkfold(b6)and af%p==a9%p and ae%p==b6%p and(kind>0 or(af%p==g1 and ae%p==g2))))
    if not f2:fixed(entry,f1,n,word,a9%2**64>=2**64-t1*C0)
 tally("lookup_sum_fold",cases);cases=[]
 for kind in(0,1,2,3):
  for _ in range(12):u=J.rr(1,2**66+1);m=u*C0;aa=(2**256-m+J.rr(m),2**256-m-1-J.rr(2**20),J.getrandbits(256),J.getrandbits(127)<<128|(2**128-m+J.rr(m)))[kind];Z=aa+(u<<256);regs=ars(aa,("x10","x11","x12","x13"));regs.update(x7=u&M64,x24=u>>64);o,n,word=site("t_fdx",regs);ov=rgi(o,("x10","x11","x12","x13"));t=fold(Z);cases.append((t>=2**256,ov!=t,ov==t%2**256));fixed("t_fdx",b7(Z,2),n,word,t%2**128<u*C0 and t>=2**128)
 tally("fold",cases);cases=[]
 for kind in(0,1,2,3,4):
  for _ in range(8 if kind<3 else 1):
   while True:
    T=(J.rr(2**64),J.rr(2**65,2**80),J.rr(p),2,(1<<128)+2)[kind];a=E.sqrt_mod_p(T)
    if a is not None:break
   a+=J.choice((0,p))if a<2**256-p else 0;o,n,word=site("t_sqinv",ars(a,("x8","x9","x10","x11")));ov=rgi(o,("x13","x14","x15","x16"));kept=rgi(o,("x8","x9","x10","x11"))==a;V=fold(a*a);t=fold(V);cases.append((t>=p,ov!=T,ov==t%2**256 and kept));fixed("t_sqinv",b7(V,2),n,word,t%2**128<2**65)
 tally("square_fold_inv2",cases);cases=[]
 for kind in(0,1,2,3):
  for _ in range(10):u=J.rr(1,2**64);m=u*C0;aa=(2**256-m+J.rr(m),2**256-m-1-J.rr(2**20),J.getrandbits(256),J.getrandbits(127)<<128|(2**128-m+J.rr(m)))[kind];Z=aa+(u<<256);regs=ars(aa,("x10","x11","x12","x13"));regs.update(x7=u);o,n,word=site("t_f5w",regs);ov=rgi(o,("x10","x11","x12","x13"));t=fold(Z);cases.append((t>=2**256,ov!=t,ov==t%2**256));fixed("t_f5w",b7(Z,2),n,word,t%2**128<u*C0 and t>=2**128)
 tally("bday_f5_fold",cases);p1_ok=True;b8=0
 for _ in range(10):
  X=E.PROD.mul(J.rr(1,E.r),E.GP);Y=E.PROD.mul(J.rr(1,E.r),E.GP)
  if X[0]==Y[0]:continue
  a,yx=X;b,yy=Y;ad=pow((yy*yy-yx*yx)%p,p-2,p);f3=fold(fold(ad*ad))%2**256;mem={};X3=0x200000;E.popax(rprog,rlab,mem,0,X,f3)
  for k,b_ in enumerate(E.ax_bases(0)):mem[X3+256+8*k]=b_
  for k in range(4):put(mem,SP+128+32*k,yy*pow(2,64*k,p)%p)
  o,n,word=site("t_p1",dict(ars(ad,("x8","x9","x10","x11")),x3=X3),mem);d3=o["x28"]|(o["x29"]<<64)|(mem[SP+496]<<128)|(o["x1"]<<192);d5=rgi(o,("x0","x2","x5","x6"));d6=fold(fold(ad*ad))%2**256;d7=E.lkfold(ad+E.lookup6(2*yx*yx%p,d6));f4=E.lookup6(2*yx%p,d6);d8=fold(sum(((f4>>64*k)&M64)*(yy*pow(2,64*k,p)%p)for k in range(4)))%2**256;p1_ok&=(d3==E.fs_op(d7,d8)and d5==E.faddw_op(d7,d8)and d3%p==pow(yy+yx,-2,p)and d5%p==pow(yy-yx,-2,p))
  if word==0:p1_ok&=(n==8+COMMON["t_p1"]);b8+=1
 ac["p1_direct_u2v2"]=dict(ok=p1_ok and b8>0,clean_runs=b8);d9,e0=[],[];pairs=[(J.rr(p),J.rr(p))for _ in range(16)]
 for _ in range(8):a=J.rr(p>>1);pairs.append((a,a+(J.rr(1,2**190)<<64)-J.rr(C0)));b=J.rr(p>>1,p);pairs.append((2**256-b+(J.rr(1,2**120)<<64)+2**64-1-J.rr(C0),b));e1=p+1+J.rr(C0-1);f5=J.rr(e1-p);pairs.append((f5,e1));pairs.append((p+C0-1-J.rr(10),p+10+J.rr(C0-10)))
 for a,b in pairs:
  if a>=2**256 or b>=2**256:continue
  regs=ars(a,("x5","x6","x7","x8"));regs.update(ars(b,("x19","x20","x21","x22")));o,n,word=site("t_fs",regs);ov=rgi(o,("x9","x10","x11","x12"));b0=a-b+p<0;e2=a-b if a>=b else a-b+p;d9.append((b0,ov!=e2%2**256 or b0,(ov==e2)if not b0 else((a-b)%p>=p-C0+1)))
  if not b0:fixed("t_fs",a<b and(a-b)%2**64<C0,n,word,ov%2**64>=2**64-C0)
  regs=ars(a,("x13","x14","x15","x16"));regs.update(ars(b,("x17","x18","x19","x20")));o,n,word=site("t_faddw",regs);ov=rgi(o,("x13","x14","x15","x16"));b1=a+b>=2**256+p;e3=a+b if a+b<2**256 else a+b-p;e0.append((b1,ov!=e3%2**256 or b1,(ov==e3)if not b1 else(C0<=(a+b)%p<=2*C0-2)))
  if not b1:fixed("t_faddw",a+b>=2**256 and(a+b)%2**64+C0>=2**64,n,word,ov%2**64<C0)
 tally("fs_p1",d9);tally("faddw_p1",e0);cases=[]
 for kind in(0,1,2):
  for _ in range(5):
   w=J.getrandbits(256);mem={};E.popbt(rprog,rlab,mem,w);d0=J.rr(1<<40);E.popbt(rprog,rlab,mem,d0)
   if kind==2:es=[J.rr(p-p//6,p)for _ in range(3)];a5=3*2**256+(J.getrandbits(160)<<64)+2**64-1-sum(es);es+=[a5//2,a5-a5//2]
   elif kind==1:es=[p-1-J.rr(100)for _ in range(4)];es.append(4*2**256-1-sum(es));assert 0<=es[4]<p
   if kind:
    w=sum((i+1+8*J.rr(1000))<<E.W6S[i]for i in range(5));E.popbt(rprog,rlab,mem,w);E.popbt(rprog,rlab,mem,d0)
    for i in range(5):put(mem,E.bt_addr(w,i),es[i])
   w2=d0;e8=dict(ars(w,("x28","x29","x11","x1")));e8.update(ars(w2,("x0","x2","x5","x6")));o,n,word=site("t_blk",e8,mem);ov=rgi(o,("x24","x25","x26","x27"));ae=rgi(o,("x20","x21","x22","x23"));ab=sum(get(mem,E.bt_addr(w,i))for i in range(5));t=ab>>256;b5=E.blook(w)if kind==0 else fold(fold(ab))%2**256;cases.append((fold(ab)>=2**256,ov!=fold(ab)%2**256,ov==b5 and ov%p==ab%p and(kind>0 or ae==E.blook(w2))));fixed("t_blk",ab%2**64+t*C0>=2**64,n,word,ab%2**64>=2**64-t*C0)
 tally("beta_ls64_fold",cases);a4,a6={},0x200000;put(a4,a6,E.BETA,4);put(a4,a6+32,0,4);put(a4,a6+64,0,4);put(a4,a6+96,0,4);nb0=E.a64.a64_run(rprog,rlab,"tbld",dict(x0=E.BT0,x1=a6,x2=8,x3=2,x4=2,x5=64,x6=0,x30=E.SNT),a4)
 for idx in range(4):put(a4,a6+64*idx,E.BETA*pow(2,E.W6S[idx+1],p)%p);put(a4,a6+64*idx+32,0)
 nb1=E.a64.a64_run(rprog,rlab,"tbld",dict(x0=E.BT1,x1=a6,x2=8,x3=2,x4=4,x5=128,x6=0,x30=E.SNT),a4);ah=(nb0==4+2*11+20*16 and nb1==4+4*11+20*32 and all(get(a4,E.BT0+64*t)==E.BETA*t%p and get(a4,E.BT0+64*t+32)==0 for t in range(8))and all(get(a4,E.BT1+128*t+32*k)==E.BETA*t*pow(2,E.W6S[k+1],p)%p for k in range(4)for t in range(8)));a7,f7,e4,e5,ak={},0x200000,J.rr(1,p),J.rr(1,p),0;pb=lambda k:E.ax_bases(0)[k]
 for k,sh in enumerate(E.W6S):dk=f7+128*k;put(a7,dk,e4*pow(2,sh,p)%p);put(a7,dk+32,0);put(a7,dk+64,e5*pow(2,sh,p)%p);put(a7,dk+96,0);ak+=E.a64.a64_run(rprog,rlab,"tbld",dict(x0=pb(k),x1=dk,x2=8,x3=2,x4=2,x5=64,x6=3*(k>0),x30=E.SNT),a7)
 ai=(ak==5*26+20*(2*8+4*4)and all(get(a7,pb(k)+64*t)==e4*t*pow(2,sh,p)%p and get(a7,pb(k)+64*t+32)==e5*t*pow(2,sh,p)%p for k,sh in enumerate(E.W6S)for t in range(2 if k>0 else 8))and all((pb(k)+64*t+o)not in a7 for k in range(1,5)for t in range(2,8)for o in(0,32)));a8,b2,f8,ao,ap,aq,ar={},E.XT0,0x200000,J.rr(1,p),J.rr(1,p),J.rr(1,p),J.rr(1,p);nc=0
 for k,sh in enumerate(E.W6S):dk=f8+128*k;put(a8,dk,ao*pow(2,sh,p)%p);put(a8,dk+32,ap if k==0 else 0);put(a8,dk+64,aq*pow(2,sh,p)%p);put(a8,dk+96,ar if k==0 else 0);nc+=E.a64.a64_run(rprog,rlab,"tbld",dict(x0=b2+E.D6O[k],x1=dk,x2=8,x3=2,x4=2,x5=64,x6=3*(k>0),x30=E.SNT),a8)
 aj=(nc==5*26+20*(2*8+4*4)and all(get(a8,b2+E.D6O[k]+64*t)==(ao*t*pow(2,E.W6S[k],p)+(ap if k==0 else 0))%p and get(a8,b2+E.D6O[k]+64*t+32)==(aq*t*pow(2,E.W6S[k],p)+(ar if k==0 else 0))%p for k in range(5)for t in range(2 if k>0 else 8))and all((b2+E.D6O[k]+64*t+o)not in a8 for k in range(1,5)for t in range(2,8)for o in(0,32)));b9=True
 for a0,a1,a2,a3 in(((1<<64)-1,(1<<64)-1,(1<<64)-1,(1<<64)-1),((1<<64)-1,0,(1<<64)-1,(1<<64)-1),(0,0,0,0),(J.getrandbits(64),J.getrandbits(64),J.getrandbits(64),J.getrandbits(64))):
  cv=a0|(a1<<64)|(a2<<128)|(a3<<192);ms=[J.getrandbits(256)for _ in range(4)];e6=dict(x8=a0,x9=a1,x10=a2,x11=a3,**E.ars(ms[0],("x28","x29","x30","x1")),**E.ars(ms[1],("x24","x25","x26","x27")),**E.ars(ms[2],("x0","x2","x5","x6")),**E.ars(ms[3],("x20","x21","x22","x23")));mt={FXA:0,SP+496:e6.pop("x30")};ot,nt,wt=site("t_tfp",e6,mt);Rc=E.shr_chain(cv);e7={};f9=[(ot[h_]<<64)|ot[l_]for l_,h_ in(("x12","x13"),("x14","x15"),("x16","x17"),("x18","x19"))];b9&=(nt==8+COMMON["t_tfp"]and wt==0 and f9==[E.tfp(Rc,m_,e7)for m_ in ms]and e7["sum_ok"])
 e9=(ah and ai and aj and b9 and ac["p1_direct_u2v2"]["ok"]and all(d["predicate_holds"]>0 and d["deviates_iff_predicate"]and d["model"]for d in ac.values()if"model"in d)and all(d["fixups"]>0 and d["word_ok"]and d["count_ok"]and d["word_set_ok"]for d in ag.values()));return dict(all_true=e9,tfp_chain_ok=b9,tbld_unified=ah and ai and aj,p1_ok=ac["p1_direct_u2v2"]["ok"],sites=len(ac),fixup_entries=len(ag))

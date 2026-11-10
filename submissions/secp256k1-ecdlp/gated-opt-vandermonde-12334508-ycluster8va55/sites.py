import os,random;random.Random.rr=random.Random.randrange;SP=(34<<40)+0x7000000
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
t_lke:
ENTER
ldp x5,x6,[sp,#384]
ldp x7,x26,[sp,#400]
ldp x27,x28,[sp,#416]
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
ldp x27,x28,[sp,#352]
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
ldp x6,x7,[x1,#288]
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
COMMON=dict(t_shx=4,t_fdx=6,t_lke=4+81,t_lkf=3+70,t_p1=233+1,t_sqinv=58+5,t_f5w=5,t_fs=7,t_faddw=7,t_blk=1+3+44+44,t_tfp=97+1);EXTRA=dict(t_shx=3,t_fdx=3,t_lke=6,t_lkf=6,t_sqinv=3,t_f5w=3,t_fs=4,t_faddw=4,t_blk=6)
def site_checks(E,I,cb,c4):
 p,C0,M64,put,get=E.p,E.C0,E.M64,E.put,E.get;rgi,ars=E.rgi,E.ars;bv,labels=E.a64 .a64_load(open(os.path.join(E.HERE,"kernels","leaf.s")).read()+SITE_KERNELS+T_KERNELS);rprog,rlab=E.a64 .a64_load(open(os.path.join(E.HERE,"kernels","route.s")).read());al,b7={},{};LS=([(E.AXYX0,E.AXY_END),(E.XT0,1<<64),(SP+128,SP+320)],[(0x200000,0x200100)])
 def site(entry,regs,mem=None):
  mem={}if mem is None else mem;out={};n=E.a64 .a64_run(bv,labels,entry,dict(regs,sp=SP,x30=E.SNT),mem,out_regs=out,ls64=LS);return out,n,mem[SP]
 def tally(bx,cases):al[bx]=dict(cases=len(cases),predicate_holds=sum(c[0]for c in cases),deviates_iff_predicate=all(c[0]==c[1]for c in cases),model=all(c[2]for c in cases))
 def fixed(entry,f,n,word,b6):
  d=b7 .setdefault(entry,dict(cases=0,fixups=0,word_ok=True,count_ok=True,word_set_ok=True));d["cases"]+=1;d["fixups"]+=f;d["word_ok"]&=word==C0*f;d["count_ok"]&=n==6+COMMON[entry]+EXTRA[entry]*f;d["word_set_ok"]&=b6 or not f
 fold=lambda v:(v%2**256)+(v>>256)*C0;a6=lambda v,k:(v%2**(64*k))+(v>>256)*C0>=2**(64*k);cases=[]
 for kind in(0,1,2,3):
  for _ in range(12):
   b3=I.rr(2**40,2**64);hi=(b3*C0)>>64;low=(2**192-hi+I.rr(hi),2**192-hi-1-I.rr(2**20),I.getrandbits(192),I.getrandbits(126)<<64|(2**64-hi+I.rr(hi)))[kind];B=low+(b3<<192);o,n,word=site("t_shx",ars(B,("x25","x26","x27","x28")));ov=rgi(o,("x5","x25","x26","x27"));ao=B*2**64-b3*p;cases.append((low+hi>=2**192,ov!=ao,ov==ao%2**256));fixed("t_shx",(B%2**64)+hi>=2**64,n,word,B%2**64>2**64-C0)
 tally("row_shift",cases);cases=[]
 for entry,ax,c5,ci in(("t_lke",(384,392,400,408,416,424),312,0),("t_lkf",(96,104,112,120,352,360),272,1)):
  for kind in(0,0,1,2):
   for _ in range(4):
    a1,a2=I.rr(1,p),I.rr(1,p);tb,X2=E.XT0+(ci*2560*(1<<40)),0x200040;mem={};mem[X2+c5]=tb;w=I.getrandbits(256);bw=E.offs6(w)
    for k,c6 in enumerate(ax):mem[SP+c6]=bw[k]
    b4=I.getrandbits(256)if entry=="t_lke"else 0;c2,bk=E.xcs(a1)[ci],(0 if ci==0 else E.xcs(a1)[3]);c1,bj=E.xcs(a2)[ci],(0 if ci==0 else E.xcs(a2)[3]);E.pop60(rprog,rlab,mem,tb,c2,bk,c1,bj,w)
    if kind:
     b4=(p-1-I.rr(100)if kind==1 else I.rr(2**200))if entry=="t_lke"else 0;T=((6 if entry=="t_lke"else 5)*2**256-1 if kind==1 else 4*2**256+(I.getrandbits(160)<<64)+2**64-1-I.rr(C0));es=[p-1-I.rr(100)for _ in range(4 if kind==1 else 3)];J=T-b4-sum(es);es+=[J//2,J-J//2]if kind==1 else[J//3,J//3,J-2*(J//3)];assert len(es)==6 and all(0<=e<p for e in es);ca=[I.rr(p//8)for _ in range(6)]
     for k in range(6):
      tk=(w>>E.W6S[k])&E.W6M[k];put(mem,tb+E.D6O[k]+(tk<<6),es[k]);put(mem,tb+E.D6O[k]+(tk<<6)+32,ca[k])
    put(mem,SP+256,b4);put(mem,SP+288,0);o,n,word=site(entry,dict(x2=X2),mem);ay,aw=get(mem,X2),get(mem,X2+32);aa=sum(get(mem,tb+E.D6O[k]+(((w>>E.W6S[k])&E.W6M[k])<<6))for k in range(6))+b4;ap=sum(get(mem,tb+E.D6O[k]+(((w>>E.W6S[k])&E.W6M[k])<<6)+32)for k in range(6))+b4;t1,t2=aa>>256,ap>>256;f1=aa%2**64+t1*C0>=2**64;f2=ap%2**64+t2*C0>=2**64;be=(E.xcs(a1)[ci]*w+(b4 if ci==0 else E.xcs(a1)[3]))%p;bd=(E.xcs(a2)[ci]*w+(b4 if ci==0 else E.xcs(a2)[3]))%p;cases.append((fold(aa)>=2**256,ay!=fold(aa)%2**256,ay==E.lkfold(aa)and aw==E.lkfold(ap)and ay%p==aa%p and aw%p==ap%p and(kind>0 or(ay%p==be and aw%p==bd))))
    if not f2:fixed(entry,f1,n,word,aa%2**64>=2**64-t1*C0)
 tally("lookup_sum_fold",cases);cases=[]
 for kind in(0,1,2,3):
  for _ in range(12):
   u=I.rr(1,2**66+1);m=u*C0;low=(2**256-m+I.rr(m),2**256-m-1-I.rr(2**20),I.getrandbits(256),I.getrandbits(127)<<128|(2**128-m+I.rr(m)))[kind];Z=low+(u<<256);regs=ars(low,("x10","x11","x12","x13"));regs.update(x7=u&M64,x24=u>>64);o,n,word=site("t_fdx",regs);ov=rgi(o,("x10","x11","x12","x13"));t=fold(Z);cases.append((t>=2**256,ov!=t,ov==t%2**256));fixed("t_fdx",a6(Z,2),n,word,t%2**128<u*C0 and t>=2**128)
 tally("fold",cases);cases=[]
 for kind in(0,1,2,3,4):
  for _ in range(8 if kind<3 else 1):
   while True:
    T=(I.rr(2**64),I.rr(2**65,2**80),I.rr(p),2,(1<<128)+2)[kind];a=E.sqrt_mod_p(T)
    if a is not None:break
   a+=I.choice((0,p))if a<2**256-p else 0;o,n,word=site("t_sqinv",ars(a,("x8","x9","x10","x11")));ov=rgi(o,("x13","x14","x15","x16"));bz=rgi(o,("x8","x9","x10","x11"))==a;V=fold(a*a);t=fold(V);cases.append((t>=p,ov!=T,ov==t%2**256 and bz));fixed("t_sqinv",a6(V,2),n,word,t%2**128<2**65)
 tally("square_fold_inv2",cases);cases=[]
 for kind in(0,1,2,3):
  for _ in range(10):
   u=I.rr(1,2**64);m=u*C0;low=(2**256-m+I.rr(m),2**256-m-1-I.rr(2**20),I.getrandbits(256),I.getrandbits(127)<<128|(2**128-m+I.rr(m)))[kind];Z=low+(u<<256);regs=ars(low,("x10","x11","x12","x13"));regs.update(x7=u);o,n,word=site("t_f5w",regs);ov=rgi(o,("x10","x11","x12","x13"));t=fold(Z);cases.append((t>=2**256,ov!=t,ov==t%2**256));fixed("t_f5w",a6(Z,2),n,word,t%2**128<u*C0 and t>=2**128)
 tally("bday_f5_fold",cases);p1_ok=True;a4=0
 for _ in range(10):
  X=E.PROD.mul(I.rr(1,E.r),E.GP);Y=E.PROD.mul(I.rr(1,E.r),E.GP)
  if X[0]==Y[0]:continue
  a,yx=X;b,yy=Y;am=pow((yy*yy-yx*yx)%p,p-2,p);b2=fold(fold(am*am))%2**256;mem={};X3=0x200000;E.popax(rprog,rlab,mem,0,X,b2)
  for k,b_ in enumerate(E.ax_bases(0)):mem[X3+256+8*k]=b_
  for k in range(4):put(mem,SP+128+32*k,yy*pow(2,64*k,p)%p)
  o,n,word=site("t_p1",dict(ars(am,("x8","x9","x10","x11")),x3=X3),mem);av=o["x28"]|(o["x29"]<<64)|(mem[SP+496]<<128)|(o["x1"]<<192);au=rgi(o,("x0","x2","x5","x6"));b5=fold(fold(am*am))%2**256;bc=E.lkfold(am+E.lookup6(2*yx*yx%p,b5));c3=E.lookup6(2*yx%p,b5);ba=fold(sum(((c3>>64*k)&M64)*(yy*pow(2,64*k,p)%p)for k in range(4)))%2**256;p1_ok&=(av==E.fs_op(bc,ba)and au==E.faddw_op(bc,ba)and av%p==pow(yy+yx,-2,p)and au%p==pow(yy-yx,-2,p))
  if word==0:p1_ok&=(n==6+COMMON["t_p1"]);a4+=1
 al["p1_direct_u2v2"]=dict(ok=p1_ok and a4>0,clean_runs=a4);ah,ai=[],[];pairs=[(I.rr(p),I.rr(p))for _ in range(16)]
 for _ in range(8):
  a=I.rr(p>>1);pairs.append((a,a+(I.rr(1,2**190)<<64)-I.rr(C0)));b=I.rr(p>>1,p);pairs.append((2**256-b+(I.rr(1,2**120)<<64)+2**64-1-I.rr(C0),b));az=p+1+I.rr(C0-1);bm=I.rr(az-p);pairs.append((bm,az));pairs.append((p+C0-1-I.rr(10),p+10+I.rr(C0-10)))
 for a,b in pairs:
  if a>=2**256 or b>=2**256:continue
  regs=ars(a,("x5","x6","x7","x8"));regs.update(ars(b,("x19","x20","x21","x22")));o,n,word=site("t_fs",regs);ov=rgi(o,("x9","x10","x11","x12"));a8=a-b+p<0;aj=a-b if a>=b else a-b+p;ah.append((a8,ov!=aj%2**256 or a8,(ov==aj)if not a8 else((a-b)%p>=p-C0+1)))
  if not a8:fixed("t_fs",a<b and(a-b)%2**64<C0,n,word,ov%2**64>=2**64-C0)
  regs=ars(a,("x13","x14","x15","x16"));regs.update(ars(b,("x17","x18","x19","x20")));o,n,word=site("t_faddw",regs);ov=rgi(o,("x13","x14","x15","x16"));a9=a+b>=2**256+p;ak=a+b if a+b<2**256 else a+b-p;ai.append((a9,ov!=ak%2**256 or a9,(ov==ak)if not a9 else(C0<=(a+b)%p<=2*C0-2)))
  if not a9:fixed("t_faddw",a+b>=2**256 and(a+b)%2**64+C0>=2**64,n,word,ov%2**64<C0)
 tally("fs_p1",ah);tally("faddw_p1",ai);cases=[]
 for kind in(0,1,2):
  for _ in range(5):
   w=I.getrandbits(256);mem={};E.popbt(rprog,rlab,mem,w);ar=I.rr(1<<40);E.popbt(rprog,rlab,mem,ar)
   if kind==2:
    es=[I.rr(p-p//7,p)for _ in range(4)];J=4*2**256+(I.getrandbits(160)<<64)+2**64-1-sum(es);es+=[J//2,J-J//2]
   elif kind==1:
    es=[p-1-I.rr(100)for _ in range(5)];es.append(5*2**256-1-sum(es));assert 0<=es[5]<p
   if kind:
    w=sum((i+1+8*I.rr(1000))<<(43*i)for i in range(6));E.popbt(rprog,rlab,mem,w)
    for i in range(6):put(mem,E.bt_addr(w,i),es[i])
   w2=ar;ag=dict(ars(w,("x28","x29","x11","x1")));ag.update(ars(w2,("x0","x2","x5","x6")));o,n,word=site("t_blk",ag,mem);ov=rgi(o,("x24","x25","x26","x27"));aw=rgi(o,("x20","x21","x22","x23"));ae=sum(get(mem,E.bt_addr(w,i))for i in range(6));t=ae>>256;ao=E.blook(w)if kind==0 else fold(fold(ae))%2**256;cases.append((fold(ae)>=2**256,ov!=fold(ae)%2**256,ov==ao and ov%p==ae%p and(kind>0 or aw==E.blook(w2))));fixed("t_blk",ae%2**64+t*C0>=2**64,n,word,ae%2**64>=2**64-t*C0)
 tally("beta_ls64_fold",cases);ab,ad,aq={},E.BT0,0x200000
 for idx,wi in enumerate(E.SLOTW):put(ab,aq+64*idx,E.BETA*pow(2,43*wi,p)%p);put(ab,aq+64*idx+32,0)
 nb=E.a64 .a64_run(rprog,rlab,"tbld",dict(x0=ad,x1=aq,x2=8,x3=2,x4=7,x5=256,x6=1<<5,x30=E.SNT),ab);bs=(nb==4+7*11+20*50 and all(get(ab,ad+256*t+32*k)==E.BETA*t*pow(2,43*E.SLOTW[k],p)%p for k in range(7)for t in range(2 if k==5 else 8))and all((ad+256*t+32*5)not in ab and(ad+256*t+32*6)in ab for t in range(2,8)));a5,bf,b9,b8,bo={},0x200000,I.rr(1,p),I.rr(1,p),0;pb=lambda k:E.ax_bases(0)[k]
 for k,sh in enumerate(E.W6S):
  dk=bf+128*k;put(a5,dk,b9*pow(2,sh,p)%p);put(a5,dk+32,0);put(a5,dk+64,b8*pow(2,sh,p)%p);put(a5,dk+96,0);bo+=E.a64 .a64_run(rprog,rlab,"tbld",dict(x0=pb(k),x1=dk,x2=8,x3=2,x4=2,x5=64,x6=3*(k in(2,5)),x30=E.SNT),a5)
 bn=(bo==6*26+20*(4*16+2*4)and all(get(a5,pb(k)+64*t)==b9*t*pow(2,sh,p)%p and get(a5,pb(k)+64*t+32)==b8*t*pow(2,sh,p)%p for k,sh in enumerate(E.W6S)for t in range(2 if k in(2,5)else 8))and all((pb(k)+64*t+o)not in a5 for k in(2,5)for t in range(2,8)for o in(0,32)));a7,ac,bh,br,bu,bq,bt={},E.XT0,0x200000,I.rr(1,p),I.rr(1,p),I.rr(1,p),I.rr(1,p);nc=0
 for k,sh in enumerate(E.W6S):
  dk=bh+128*k;put(a7,dk,br*pow(2,sh,p)%p);put(a7,dk+32,bu if k==0 else 0);put(a7,dk+64,bq*pow(2,sh,p)%p);put(a7,dk+96,bt if k==0 else 0);nc+=E.a64 .a64_run(rprog,rlab,"tbld",dict(x0=ac+E.D6O[k],x1=dk,x2=8,x3=2,x4=2,x5=64,x6=3*(k in(2,5)),x30=E.SNT),a7)
 bp=(nc==6*26+20*(4*16+2*4)and all(get(a7,ac+E.D6O[k]+64*t)==(br*t*pow(2,E.W6S[k],p)+(bu if k==0 else 0))%p and get(a7,ac+E.D6O[k]+64*t+32)==(bq*t*pow(2,E.W6S[k],p)+(bt if k==0 else 0))%p for k in range(6)for t in range(2 if k in(2,5)else 8))and all((ac+E.D6O[k]+64*t+o)not in a7 for k in(2,5)for t in range(2,8)for o in(0,32)));af=True
 for a0,a1,a2,a3 in(((1<<64)-1,(1<<64)-1,(1<<64)-1,(1<<64)-1),((1<<64)-1,0,(1<<64)-1,(1<<64)-1),(0,0,0,0),(I.getrandbits(64),I.getrandbits(64),I.getrandbits(64),I.getrandbits(64))):
  cv=a0|(a1<<64)|(a2<<128)|(a3<<192);ms=[I.getrandbits(256)for _ in range(4)];b1=dict(x8=a0,x9=a1,x10=a2,x11=a3,**E.ars(ms[0],("x28","x29","x30","x1")),**E.ars(ms[1],("x24","x25","x26","x27")),**E.ars(ms[2],("x0","x2","x5","x6")),**E.ars(ms[3],("x20","x21","x22","x23")));mt={SP:0,SP+496:b1 .pop("x30")};ot,nt,wt=site("t_tfp",b1,mt);Rc=E.shr_chain(cv);b0={};bg=[(ot[h_]<<64)|ot[l_]for l_,h_ in(("x12","x13"),("x14","x15"),("x16","x17"),("x18","x19"))];af&=(nt==6+COMMON["t_tfp"]and wt==0 and bg==[E.tfp(Rc,m_,b0)for m_ in ms]and b0["sum_ok"])
 bi=(bs and bn and bp and af and al["p1_direct_u2v2"]["ok"]and all(d["predicate_holds"]>0 and d["deviates_iff_predicate"]and d["model"]for d in al.values()if"model"in d)and all(d["fixups"]>0 and d["word_ok"]and d["count_ok"]and d["word_set_ok"]for d in b7 .values()));return dict(all_true=bi,tfp_chain_ok=af,tbld_unified=bs and bn and bp,p1_ok=al["p1_direct_u2v2"]["ok"],sites=len(al),fixup_entries=len(b7))

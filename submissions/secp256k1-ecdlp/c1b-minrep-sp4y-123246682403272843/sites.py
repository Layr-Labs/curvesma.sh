import random,math;random.Random.rr=random.Random.randrange;SP=(34<<40)+0x7000000;rn,en,dc=range,enumerate,dict;SITE_KERNELS="""
.macro ENTER
mov x4,#977
movk x4,#1,lsl #32
str x30,[sp,#504]
str xzr,[sp,#-128]
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
ldp x28,x29,[x1,#32]
ldr x6,[x1,#288]
ldp x7,x30,[x1,#48]
LK1 0x2000000000000
LEAVE
t_p1:
ENTER
add x29,x3,#256
P1_6WIN_BFI
str x30,[sp,#368]
LEAVE
t_fs:
ENTER
FS x9,x10,x11,x12,x5,x6,x7,x8,x19,x20,x21,x22,x13
LEAVE
t_kfpi:
ENTER
ld64b x18,[x3]
KD
st64b x18,[x3]
mov x30,x11
mov x7,#0xfff8000000000000
orr x8,x7,#64
orr x9,x7,#128
BLK6F x28,x29,x30,x1
ROWSG x28,x29,x30,x1
KM0 x3,x1,x30,x29,x28
LEAVE
""";T_KERNELS="""
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
str x30,[sp,#-136]
mov x4,#977
movk x4,#1,lsl #32
EFGK6
ldr x30,[sp,#-136]
ret
""";COMMON=dc(t_shx=4,t_fdx=6,t_lk1=3+38,t_p1=218+2,t_fs=7,t_kfpi=1+3+36+8+44);EXTRA=dc(t_shx=4,t_fdx=4,t_lk1=7,t_fs=5,t_kfpi=4)
def site_checks(E,I,cb,c4):
 p,C0,M64,put,get=E.p,E.C0,E.M64,E.put,E.get;rgi,ars=E.rgi,E.ars;bv,labels=E.a64.a64_load(E.ksrc("leaf")+SITE_KERNELS+T_KERNELS);rprog,rlab=E.kload("route");al,b7={},{};LS=([(E.AXYX0,E.AXY_END),(E.XT0,1<<64),(SP+64,SP+320),(0x200100,0x200140)],[(0x200000,0x200140)])
 def site(entry,rs,mem=None):mem={}if mem is None else mem;out={};n=E.arun(bv,labels,entry,dc(rs,sp=SP+128),mem,out_regs=out,ls64=LS);return out,n,mem[SP]
 def tally(bx,cz):al[bx]=dc(cases=len(cz),predicate_holds=sum(c[0]for c in cz),deviates_iff_predicate=all(c[0]==c[1]for c in cz),model=all(c[2]for c in cz))
 def fixed(entry,f,n,wd,b6):d=b7.setdefault(entry,dc(cases=0,fixups=0,word_ok=True,count_ok=True,word_set_ok=True));d["cases"]+=1;d["fixups"]+=f;d["word_ok"]&=wd==C0*f;d["count_ok"]&=n==6+COMMON[entry]+EXTRA[entry]*f;d["word_set_ok"]&=b6 or not f
 fold=lambda v:v%2**256+(v>>256)*C0;a6=lambda v,k:v%2**(64*k)+(v>>256)*C0>=2**(64*k);cz=[]
 for kd in 0,1,2,3:
  for _ in rn(12):b3=I.rr(2**40,2**64);hi=b3*C0>>64;low=(2**192-hi+I.rr(hi),2**192-hi-1-I.rr(2**20),I.gb(192),I.gb(126)<<64|2**64-hi+I.rr(hi))[kd];B=low+(b3<<192);o,n,wd=site("t_shx",ars(B,("x25","x26","x27","x28")));ov=rgi(o,("x5","x25","x26","x27"));ao=B*2**64-b3*p;cz+=[(low+hi>=2**192,ov!=ao,ov==ao%2**256)];fixed("t_shx",B%2**64+hi>=2**64,n,wd,B%2**64>2**64-C0)
 tally("row_shift",cz);cz=[]
 for kd in 0,-1,1,2:
  for _ in rn(4):
   mem,w,B2=({},I.gb(256),I.rr(p));g0=[E.XTP[0]+v for v in E.offs6(w)];tb=E.XTP[1];mem[SP+360]=0
   for k,pos in en((96,104,112,120,352)):mem[SP+pos]=g0[k]
   if kd<1:es=[I.rr(p+kd*B2)//(1-4*kd)for _ in rn(5)]
   else:T=5*2**256+2**64-1-I.rr(C0)if kd==1 else 5*2**256-1-I.rr(C0);es=[p-1-I.rr(100)for _ in rn(4 if kd==1 else 2)];J=T-B2-sum(es);es+=[J]if kd==1 else[J//3,J//3,J-2*(J//3)];assert len(es)==5 and all(0<=e<p for e in es)
   for k,sh,mi,_ in E.W5V:tk=w>>sh&mi;put(mem,tb+E.D5O[k]+(tk<<6),I.gb(256));put(mem,tb+E.D5O[k]+(tk<<6)+32,es[k])
   put(mem,SP+64,B2);o,n,wd=site("t_lk1",dc(x1=SP+64),mem);ov=rgi(o,("x20","x21","x22","x23"));S_=sum(es)+B2;t=S_>>256;cz+=[(S_%2**256+t*C0>=2**256,ov!=(S_%2**256+t*C0)%2**256,ov==E.lkfold(S_))];fixed("t_lk1",S_%2**64+t*C0>=2**64,n,wd,S_%2**64>=2**64-t*C0)
 tally("lookup_sum_fold",cz);cz=[]
 for kd in 0,1,2,3:
  for _ in rn(12):u=I.rr(1,2**66+1);m=u*C0;low=(2**256-m+I.rr(m),2**256-m-1-I.rr(2**20),I.gb(256),I.gb(127)<<128|2**128-m+I.rr(m))[kd];Z=low+(u<<256);rs=ars(low,("x10","x11","x12","x13"));rs.update(x7=u&M64,x24=u>>64);o,n,wd=site("t_fdx",rs);ov=rgi(o,("x10","x11","x12","x13"));t=fold(Z);cz+=[(t>=2**256,ov!=t,ov==t%2**256)];fixed("t_fdx",a6(Z,2),n,wd,t%2**128<u*C0 and t>=2**128)
 tally("fold",cz);cz=[]
 c1,f12,X3c=[],[0,0],0x300000;lo=lambda f:I.gb(126)<<128|C0+I.rr(C0<<19)if f else I.rr(p)
 for i in rn(16):
  f1,f2=i%4 in(1,3),i%4 in(2,3)
  while True:
   s1,s2,dA=lo(f1),lo(f2),I.rr(1,2**63);ie=(s2-s1+(f1<<128))*pow(dA,-1,p)%p
   if not ie:continue
   D=s1*pow(ie,-1,p)%p;Ao=I.rr(p);w1m,w2m,fl,Z,V=E.c1bm(ie,D,dA)
   if fl==[f1,f2]:break
  mem={};put(mem,SP+288,(D+Ao)%p);put(mem,X3c-160,Ao);mem[X3c+176],mem[X3c+184]=dA,X3c+64;o,n,wd=site("t_c1b",dc(ars(ie,("x25","x26","x27","x28")),x3=X3c),mem)
  w1,w2,c1a,c2a=get(mem,X3c-192),rgi(o,("x8","x9","x10","x11")),Z%2**256+(Z>>256)*C0,V%2**256+(V>>256)*C0
  c1+=[w1==w1m==c1a-(f1<<128)and w2==w2m==c2a-(f2<<128)and o["x29"]==X3c+64 and w1%p==(s1-(f1<<128))%p and w2%p==(s2-(f2<<128))%p and(not f1 or s1%2**128<2**66*C0)and(not f2 or f1 or V%p%2**128<(2**63+2)*C0)];f12[0]+=f1;f12[1]+=f2
 R2=random.Random(48);W1,W2=2**66*C0,(2**63+2)*C0;ce,fe=[],[0,0]
 for which,W in((0,W1),(1,W2)):
  for e in[0,1,W-1,W,W+C0,2**128-C0,2**128-1]+[R2.gb(128)for _ in rn(3)]:
   r,dA,Ao=R2.gb(126)<<128|e,R2.rr(1,2**63),R2.rr(p)
   if which==0:s1,ie=r,R2.rr(1,p)
   else:s1=R2.rr(p);ie=(r-s1)*pow(dA,-1,p)%p or 1
   D=s1*pow(ie,-1,p)%p;w1m,w2m,fl,Z,V=E.c1bm(ie,D,dA);mem={};put(mem,SP+288,(D+Ao)%p);put(mem,X3c-160,Ao);mem[X3c+176],mem[X3c+184]=dA,X3c+64
   o,n,wd=site("t_c1b",dc(ars(ie,("x25","x26","x27","x28")),x3=X3c),mem)
   ce+=[get(mem,X3c-192)==w1m and rgi(o,("x8","x9","x10","x11"))==w2m and(not fl[0]or Z%p%2**128<W1)and(not fl[1]or fl[0]or V%p%2**128<W2)];fe[0]+=fl[0];fe[1]+=fl[1]
 al["c1b_drops"]=dc(cases=len(c1),fired=f12,edge_cases=len(ce),edge_fired=fe,ok=all(c1)and min(f12)>0 and all(ce)and min(fe)>0);p1_ok,fl,X3,WW=True,[0]*9,0x200000,[(1,C0-1,C0-1)]*3+[(2,0,(C0+1)*C0),(1,0,6*C0),(1,0,7*C0),(2,0,C0<<66),(1,2*C0,2*C0),(1,0,2*C0)]
 for i in rn(24):
  X=E.PROD.mul(I.rr(1,E.r),E.GP);yx,k_=X[1],i%12;y,J,q,s=I.rr(2**250,p),I.gb(200),[0]*6,[1]+[0]*5;rt=k_ in(0,1)or k_>6
  while k_>7:
   h,d=I.gb(256),k_-7
   if d<4:h=h&~(2**128-1<<64*d-64)|(M64-I.rr(2**20)|(I.rr(1,C0)<<64)//C0<<64)<<64*d-64
   else:h=h&~(2**192-1)|M64-I.rr(2**20)|(C0-1-I.rr(2**20)<<64)//C0<<64|(I.rr(1,C0)<<64)//C0<<128
   J=math.isqrt((h<<256)-1)+1;m=[(h>>64*k&M64)*C0 for k in rn(4)];w=[(m[i]&M64)+(m[i-1]>>64)>>64 for i in(1,2,3)]
   if J*J>>256==h and J<p and(w==[1,1,0]if d>3 else 1 in w and w.index(1)==d-1):break
  if k_%7==0:Y=E.PROD.mul(I.rr(1,E.r),E.GP);y=Y[1];J=pow((y*y-yx*yx)%p,p-2,p)
  while k_==1:
   J=E.sqrt_mod_p(I.rr(2**64))
   if J and E.fg(J*J%2**256+(J*J>>256)*C0,2):break
  if k_==2:s=[I.gb(256)for _ in rn(5)];s+=[(I.gb(192)<<64|M64-I.rr(C0))-sum(s)&2**256-1]
  if k_==3:q=[I.gb(256)for _ in rn(5)];q+=[(I.gb(192)<<64|M64-I.rr(C0))-J-sum(q)&2**256-1]
  if k_==4:s[0]=I.gb(64)|2**63|1;y=I.rr(2**120,2**127)<<128|(2**128-1-I.rr(2**32))*pow(s[0],-1,2**128)%2**128
  if k_==5:q[0]=y-(I.rr(1,2**150)<<64)+I.rr(C0)-J
  if k_==6:q[0]=2**256-y+(I.rr(1,2**120)<<64)+M64-I.rr(C0)-J
  T,mem=[y*pow(2,64*k,p)%p for k in rn(4)],{};b,u,v,F=E.p1m(J,T,yx)if rt else E.p1m(J,T,lambda b:sum(s),lambda b:sum(q))
  if rt:E.popax(rprog,rlab,mem,0,X,b)
  else:
   for k,(_,sh,m,_)in en(E.W6V):ad=E.ax_bases(0)[k]+((b>>sh&m)<<6);put(mem,ad,q[k]);put(mem,ad+32,s[k])
  for k,b_ in en(E.ax_bases(0)+(0x5A5A5A5A5A5A5A5A,)*2):mem[X3+256+8*k]=b_
  for k in rn(4):put(mem,SP+128+32*k,T[k])
  o,n,wd=site("t_p1",dc(ars(J,("x8","x9","x10","x11")),x3=X3),mem);av=o["x28"]|o["x29"]<<64|mem[SP+496]<<128|o["x1"]<<192;au=rgi(o,("x0","x2","x5","x6"));e=J*J%p;G,H=(2*yx*yx*e,2*yx*e)if rt else(sum(q),sum(s))
  cu,cv=[-2**(64*m+65)*yx*(yx-y)for m in(1,2,3,1)]+[2**64*y,-2**64,2**128,2**64,0],[-2**(64*m+65)*yx*(yx+y)for m in(1,2,3,1)]+[-2**64*y,-2**64,-2**128,0,-2**64];p1_ok&=(av,au)==(u,v)and(av-J-G+y*H-sum(c*f for c,f in zip(cu,F)))%p==0==(au-J-G-y*H-sum(c*f for c,f in zip(cv,F)))%p and n==6+COMMON["t_p1"]and wd==0 and(k_%7>0 or(av*(y+yx)**2%p,au*(y-yx)**2%p)==(1,1));fl=[a+f for a,f in zip(fl,F)];sr=[(J*J>>256+64*m&M64)*C0 for m in(1,2,3)]+[e,H%p,(J+G)%p,y*H%p,(J+G-y*H)%p,(J+G+y*H)%p];p1_ok&=all((sr[i]+o)%2**(64*l)<z for i,(l,o,z)in en(WW)if F[i]and(i<3 or not any(F[:i])))
 al["p1_unguarded"]=dc(ok=p1_ok and min(fl)>0,fired=fl);ah=[];pairs=[(I.rr(p),I.rr(p))for _ in rn(16)]
 for _ in rn(8):a=I.rr(p>>1);pairs+=[(a,a+(I.rr(1,2**190)<<64)-I.rr(C0))];b=I.rr(p>>1,p);pairs+=[(2**256-b+(I.rr(1,2**120)<<64)+2**64-1-I.rr(C0),b)];az=p+1+I.rr(C0-1);bm=I.rr(az-p);pairs+=[(bm,az)];pairs+=[(p+C0-1-I.rr(10),p+10+I.rr(C0-10))]
 for a,b in pairs:
  if a>=2**256 or b>=2**256:continue
  rs=ars(a,("x5","x6","x7","x8"));rs.update(ars(b,("x19","x20","x21","x22")));o,n,wd=site("t_fs",rs);ov=rgi(o,("x9","x10","x11","x12"));a8=a-b+p<0;aj=a-b if a>=b else a-b+p;ah+=[(a8,ov!=aj%2**256 or a8,ov==aj if not a8 else(a-b)%p>=p-C0+1)]
  if not a8:fixed("t_fs",a<b and(a-b)%2**64<C0,n,wd,ov%2**64>=2**64-C0)
 tally("fs_p1",ah);cz=[];kq=True
 for kd in 0,1,2,3:
  for q_ in rn(6):
   w=I.gb(256)if kd<3 else I.rr(C0)+p
   if kd:w=w&~(M64<<64)|M64<<64|(2**64-1-I.rr(2**34)if kd==1 else w&M64)
   cl="EF"[q_%2];Cb=(E.CE.C_E,E.CE.C_FG)[q_%2];c1=(-Cb,Cb,I.rr(-Cb,Cb+1))[q_%3];c2=(0,2**128-1,I.rr(2**128))[q_//2%3];mem={};E.popbt(None,None,mem,w);put(mem,0x200100,c1%2**192,3);put(mem,0x200118,c2,2);put(mem,0x200128,0,3);o,n,wd=site("t_kfpi",dc(**ars(w,("x28","x29","x11","x1")),x3=0x200100),mem);A=o["x17"]<<64|o["x16"];B=o["x19"]<<64|o["x18"];y1=rgi(o,("x24","x25","x26","x27"));f=w%2**64+((w>>192)*C0>>64)>=2**64 and(w>>64)%2**64==M64;cv=(c1+c2*E.BETA)%p;A_,B_=E.ktfp(c1,c2,w,E.blkf(w));e1=E.werr(E.t128(cv*w%p),A);e2=E.werr(E.t128(E.BETA*cv*w%p),B);kt=E.CE.KT[cl];d1=((w*E.BETA%p<<256)-p*y1)%(p<<256);ok=(A,B,y1)==(A_,B_,E.blkf(w))and(d1<=3*p or d1>=(p<<256)-3*p)and all(lo<=d<=hi for ev,r in((e1,kt[0]),(e2,kt[1]))for d,(lo,hi)in zip(ev,r));kq&=ok;cz+=[(f,E.rowsf(w)!=(w%2**128+((w>>192)*C0>>64))%2**128|w>>128<<128,ok)];fixed("t_kfpi",f,n,wd,True)
 tally("rows_kfp",cz);ab,ad,aq={},E.BT0,0x200000
 for idx,wi in en(E.SLOTW):put(ab,aq+64*idx,E.BETA*pow(2,43*wi,p)%p);put(ab,aq+64*idx+32,0)
 nb=E.arun(cb,c4,"tbldf",dc(x0=ad,x1=aq,x2=8,x3=7,x26=256),ab);bs=nb==E.CE.tbldf_count(8,7)and all(get(ab,ad+256*t+32*k)==E.frac_e(E.BETA*t*pow(2,43*E.SLOTW[k],p)%p)for k in rn(7)for t in rn(8))
 def tb(ba,c1,a1,c2,a2,W,s=0):
  m,n={},0;o=lambda k,c,a:(a if k==0 else 0)+(c*s*pow(2,215,p)if k==4 else 0)
  for k,sh,R in W:dk=0x200000+128*k;put(m,dk,c1*pow(2,sh,p)%p);put(m,dk+32,o(k,c1,a1)%p);put(m,dk+64,c2*pow(2,sh,p)%p);put(m,dk+96,o(k,c2,a2)%p);n+=E.arun(rprog,rlab,"tbld",dc(x0=ba(k),x1=dk,x2=8,x3=2,x4=2,x5=64,x6=3*(R<8)),m)
  return n==26*len(W)+40*sum(R for _,_,R in W)and all(get(m,ba(k)+64*t)==(c1*t*pow(2,sh,p)+o(k,c1,a1))%p and get(m,ba(k)+64*t+32)==(c2*t*pow(2,sh,p)+o(k,c2,a2))%p for k,sh,R in W for t in rn(R))and all((ba(k)+64*t+o_)not in m for k,_,R in W for t in rn(R,8)for o_ in(0,32))
 bn=tb(lambda k:E.ax_bases(0)[k],I.rr(1,p),0,I.rr(1,p),0,[(k,sh,2 if k in(2,5)else 8)for k,sh in en(E.W6S)]);bp=tb(lambda k:E.XTP[0]+E.D5O[k],I.rr(1,p),I.rr(1,p),I.rr(1,p),I.rr(1,p),[(k,43*k,8)for k in rn(5)],I.gb(41));bi=bs and bn and bp and kq and al["p1_unguarded"]["ok"]and al["c1b_drops"]["ok"]and all(d["predicate_holds"]>0 and d["deviates_iff_predicate"]and d["model"]for k,d in al.items()if"model"in d)and all(d["fixups"]>0 and d["word_ok"]and d["count_ok"]and d["word_set_ok"]for d in b7.values());return dc(all_true=bi,ktfp_ok=kq,tbld_unified=bs and bn and bp,p1_ok=al["p1_unguarded"]["ok"],c1b_drops=al["c1b_drops"],sites=len(al),fixup_entries=len(b7))

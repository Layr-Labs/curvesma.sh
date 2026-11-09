import os
SP=0x7000000
SITE_KERNELS="""
.macro ENTER
mov x4,#977
movk x4,#1,lsl #32
stp xzr,x30,[sp,#-16]!
.endm
.macro LEAVE
ldp xzr,x30,[sp],#16
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
t_lf:
ENTER
LF x20,x21,x22,x23,x11,x10
LEAVE
t_sqw:
ENTER
SQ13 x0,x1,x2,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14
F5 x6,x7,x8,x9,x13,x6,x7,x14,x13
LEAVE
t_fs:
ENTER
FS x9,x10,x11,x12,x5,x6,x7,x8,x19,x20,x21,x22,x13
LEAVE
t_faddw:
ENTER
FADDW x13,x14,x15,x16,x13,x14,x15,x16,x17,x18,x19,x20,x5
LEAVE
t_tfp:
ENTER
TFP4S
LEAVE
t_df:
ENTER
DT
LEAVE
t_blk:
ENTER
mov x5,#0xfff8000000000000
orr x6,x5,#64
orr x7,x5,#128
BLK6 x21,x22,x23,x24,x8,x9,x10,x11
LEAVE
"""
COMMON=dict(t_shx=4,t_fdx=6,t_lf=3,t_sqw=58+5,t_fs=7,t_faddw=7,t_blk=3+44,t_df=17+5)
EXTRA=dict(t_shx=3,t_fdx=3,t_lf=6,t_sqw=3,t_fs=4,t_faddw=4,t_blk=6,t_df=6)
def run_site(E,prog,labels,entry,regs,mem=None):
 mem={}if mem is None else mem
 out={}
 n=E.a64.a64_run(prog,labels,entry,dict(regs,sp=SP,x30=E.SENTINEL),mem,out_regs=out,ls64=([(E.BT0,1<<64)],[]))
 return out,n,mem[SP-16]
def site_checks(E,rng,fprog,flabels):
 p,C0,M64,put,get=E.p,E.C0,E.M64,E.put,E.get
 reg_int,as_regs=E.reg_int,E.as_regs
 prog,labels=E.a64.a64_load(open(os.path.join(E.HERE,"kernels","leaf.s")).read()+SITE_KERNELS)
 res,fix={},{}
 def site(entry,regs,mem=None):
  return run_site(E,prog,labels,entry,regs,mem)
 def tally(name,cases):
  res[name]=dict(cases=len(cases),predicate_holds=sum(c[0]for c in cases),
       deviates_iff_predicate=all(c[0]==c[1]for c in cases),model=all(c[2]for c in cases))
 def fixed(entry,f,n,word,in_set):
  d=fix.setdefault(entry,dict(cases=0,fixups=0,word_ok=True,count_ok=True,word_set_ok=True))
  d["cases"]+=1
  d["fixups"]+=f
  d["word_ok"]&=word==C0*f
  d["count_ok"]&=n==5+COMMON[entry]+EXTRA[entry]*f
  d["word_set_ok"]&=in_set or not f
 def fold(v):
  return (v%2**256)+(v>>256)*C0
 def low_carry(v,k):
  return (v%2**(64*k))+(v>>256)*C0>=2**(64*k)
 cases=[]
 for kind in (0,1,2,3):
  for _ in range(48):
   b3=rng.randrange(2**40,2**64)
   hi=(b3*C0)>>64
   low=(2**192-hi+rng.randrange(hi),2**192-hi-1-rng.randrange(2**20),rng.getrandbits(192),
       rng.getrandbits(126)<<64|(2**64-hi+rng.randrange(hi)))[kind]
   B=low+(b3<<192)
   o,n,word=site("t_shx",as_regs(B,("x25","x26","x27","x28")))
   o=reg_int(o,("x5","x25","x26","x27"))
   want=B*2**64-b3*p
   cases.append((low+hi>=2**192,o!=want,o==want%2**256))
   fixed("t_shx",(B%2**64)+hi>=2**64,n,word,B%2**64>2**64-C0)
 tally("row_shift",cases)
 for name,entry,top in (("fold","t_fdx",2**66),):
  cases=[]
  for kind in (0,1,2,3):
   for _ in range(48):
    u=rng.randrange(1,top+1)
    m=u*C0
    low=(2**256-m+rng.randrange(m),2**256-m-1-rng.randrange(2**20),rng.getrandbits(256),
        rng.getrandbits(127)<<128|(2**128-m+rng.randrange(m)))[kind]
    Z=low+(u<<256)
    if entry=="t_fdx":
     regs=as_regs(low,("x10","x11","x12","x13"))
     regs.update(x7=u&M64,x24=u>>64)
     outs=("x10","x11","x12","x13")
    else:
     K=rng.randrange(p)
     regs=as_regs(Z-K,("x13","x14","x15","x16","x17","x18"))
     regs.update(as_regs(K,("x9","x10","x11","x12")))
     outs=("x13","x14","x15","x16")
    o,n,word=site(entry,regs)
    o=reg_int(o,outs)
    t=fold(Z)
    cases.append((t>=2**256,o!=t,o==t%2**256))
    fixed(entry,low_carry(Z,2),n,word,t%2**128<u*C0 and t>=2**128)
  tally(name,cases)
 for name,entry,a_regs,o_regs in (("square_fold_w","t_sqw",("x0","x1","x2","x5"),("x6","x7","x8","x9")),):  # P1's w = inv^2
  cases=[]
  for kind in (0,1,2):
   for _ in range(32):
    while True:
     T=(rng.randrange(2**64),rng.randrange(2**65,2**80),rng.randrange(p))[kind]
     a=E.sqrt_mod_p(T)
     if a is not None:
      break
    a+=rng.choice((0,p))if a<2**256-p else 0
    o,n,word=site(entry,as_regs(a,a_regs))
    o=reg_int(o,o_regs)
    V=fold(a*a)
    t=fold(V)
    cases.append((t>=p,o!=T,o==t%2**256))
    fixed(entry,low_carry(V,2),n,word,t%2**128<2**65)
  tally(name,cases)
 edge=[0,1,p-1,C0,2**255]
 pairs=[(a,b)for a in edge for b in edge]+[(rng.randrange(p),rng.randrange(p))for _ in range(48)]
 for _ in range(24):
  a=rng.randrange(p>>1)
  pairs.append((a,a+(rng.randrange(1,2**190)<<64)-rng.randrange(C0)))
  b=rng.randrange(p>>1,p)
  pairs.append((2**256-b+(rng.randrange(1,2**120)<<64)+2**64-1-rng.randrange(C0),b))
 fs=fa=True
 for a,b in pairs:
  if a>=p or b>=p:
   continue
  regs=as_regs(a,("x5","x6","x7","x8"))
  regs.update(as_regs(b,("x19","x20","x21","x22")))
  o,n,word=site("t_fs",regs)
  o=reg_int(o,("x9","x10","x11","x12"))
  fs&=o==(a-b)%p
  fixed("t_fs",a<b and (a-b)%2**64<C0,n,word,o%2**64>=2**64-C0)
  regs=as_regs(a,("x13","x14","x15","x16"))
  regs.update(as_regs(b,("x17","x18","x19","x20")))
  o,n,word=site("t_faddw",regs)
  o=reg_int(o,("x13","x14","x15","x16"))
  fa&=o%p==(a+b)%p and o<2**256
  fixed("t_faddw",a+b>=2**256 and (a+b)%2**64+C0>=2**64,n,word,o%2**64<C0)
 # BLK6 window sums: real beta entries, and crafted ones forcing the fixup
 cases=[]
 for kind in (0,1,2,3):
  for _ in range(24):
   w=rng.getrandbits(256)
   mem=E.LazyTable([])
   if kind>=2:  # low 128 (kind 3: 64) bits all ones, t = 4: the fixup, no carry past 2^256
    es=[rng.randrange(p-p//7,p)for _ in range(4)]
    rest=4*2**256+2**200+2**(128 if kind==2 else 64)-1-sum(es)
    es+=[rest//2,rest-rest//2]
    assert all(0<=e<p for e in es)
    for i in range(6):
     put(mem,E.bt_addr(w,i),es[i])
    tot=sum(es)
   elif kind:
    es=[rng.randrange(p-p//7,p)for _ in range(5)]
    es.append(5*2**256-1-sum(es))
    assert 0<=es[5]<p
    for i in range(6):
     put(mem,E.bt_addr(w,i),es[i])
    tot=sum(es)
   else:
    tot=sum(get(mem,E.bt_addr(w,i))for i in range(6))
   o,n,word=site("t_blk",as_regs(w,("x21","x22","x23","x24")),mem)
   o=reg_int(o,("x8","x9","x10","x11"))
   f=int(fold(tot)>>256)
   fixed("t_blk",int(tot%2**64+(tot>>256)*C0>=2**64),n,word,tot%2**64>=2**64-(tot>>256)*C0)
   ref=fold(tot)%2**256+C0*f  # the full chain
   cases.append((f==1,o!=fold(tot),o==ref and ref<2**256 and (kind==0)==(o%p==E.BETA*w%p)))
 tally("window_sum_fold",cases)
 # LF (3.4): Z = low + c 2^256 (c <= 8) folds exactly; crafted lows carry at limbs 0, 1 and past 2^256
 cases=[]
 for kind in (0,1,2,3,4):
  for _ in range(32):
   c=rng.randrange(9)
   m=c*C0
   low=(rng.getrandbits(256),2**256-1-rng.randrange(m+1),rng.getrandbits(128)<<128|(2**128-1-rng.randrange(m+1)),
       2**256-m+rng.randrange(m+1)if m else 0,rng.getrandbits(191)<<64|(2**64-1-rng.randrange(m+1)))[kind]
   Z=low+(c<<256)
   o,n,word=site("t_lf",dict(as_regs(low,("x20","x21","x22","x23")),x11=c))
   o=reg_int(o,("x20","x21","x22","x23"))
   want=E.lfold(Z)
   fixed("t_lf",int(low%2**64+m>=2**64),n,word,low%2**64>=2**64-m)
   cases.append((int(fold(Z)>>256),o%p==Z%p,o==want and o<2**256))
 res["lookup_fold"]=dict(cases=len(cases),carry_past_2_256=sum(c[0]for c in cases),exact=all(c[1]for c in cases),
       model=all(c[2]for c in cases))
 # DT, DF (3.4): Z0 + delta*inv folds exactly; crafted sums carry out of limb 1 (fixup) and past 2^256
 cases=[]
 for kind in range(4):
  for _ in range(32):
   dl=rng.choice((rng.randrange(1,2**63),2**63-1,rng.randrange(1,2**20)))
   iv=rng.choice((rng.randrange(p),rng.randrange(p,2**256),2**256-1))
   P_=dl*iv
   m=(P_>>256)*C0
   tg=(rng.getrandbits(256),rng.getrandbits(128)<<128|(2**128-1-rng.randrange(m+1)),2**256-1-rng.randrange(m+1),
       rng.getrandbits(97))[kind]
   Z0=(tg-P_)%2**256
   S=Z0+P_
   regs=as_regs(Z0,("x0","x1","x2","x5"))
   regs.update(as_regs(iv,("x9","x10","x11","x12")),x13=dl)
   o,n,word=site("t_df",regs)
   o=reg_int(o,("x0","x1","x2","x5"))
   cs=S>>256
   fixed("t_df",int(S%2**128+cs*C0>=2**128),n,word,o%2**128<2**97)
   cases.append((int((S%2**256+cs*C0)>>256),o%p==S%p,o==E.lfold(S)and o<2**256))
 res["delta_fold"]=dict(cases=len(cases),carry_past_2_256=sum(c[0]for c in cases),exact=all(c[1]for c in cases),
       model=all(c[2]for c in cases))
 # FS, FADDW on weak operands (P1's M, N), bit-exact against fsub/fadd
 cs_s,cs_a=[],[]
 wk=[0,1,C0-1,C0,p-1,p,p+1,2**256-C0+5,2**256-1]
 for a,b in [(a,b)for a in wk for b in wk]+[(rng.choice((rng.getrandbits(256),rng.randrange(C0))),
                  rng.choice((rng.getrandbits(256),p+rng.randrange(C0))))for _ in range(96)]:
  regs=as_regs(a,("x5","x6","x7","x8"))
  regs.update(as_regs(b,("x19","x20","x21","x22")))
  o,n,word=site("t_fs",regs)
  o=reg_int(o,("x9","x10","x11","x12"))
  cs_s.append((b>a+p,o%p!=(a-b)%p,o==E.fsub(a,b)))
  regs=as_regs(a,("x13","x14","x15","x16"))
  regs.update(as_regs(b,("x17","x18","x19","x20")))
  o,n,word=site("t_faddw",regs)
  o=reg_int(o,("x13","x14","x15","x16"))
  cs_a.append((a+b>=2**257-C0,o%p!=(a+b)%p,o==E.fadd(a,b)))
 tally("fs_weak_operands",cs_s)
 tally("faddw_weak_operands",cs_a)
 # beta tables (tbld, PAIRB6), short: every entry of the 7 slots, count 4 + 7*12 + 20(6n + n5)
 rprog,rlabels=E.a64.a64_load(open(os.path.join(E.HERE,"kernels","route.s")).read())
 tb_ok=True
 for n,n5 in ((7,5),(12,3)):
  mem={}
  for k in range(7):
   put(mem,0x1000+32*k,E.BETA*pow(2,43*E.BT_WIN[k],p)%p)
  steps=E.a64.a64_run(rprog,rlabels,"tbld",dict(x0=0x100000,x1=0x1000,x2=n,x3=n5,sp=SP,x30=E.SENTINEL),mem)
  tb_ok&=steps==4+7*12+20*(6*n+n5)and all(
   get(mem,0x100000+256*t+32*k)==E.BETA*t*pow(2,43*E.BT_WIN[k],p)%p
   for k in range(7)for t in range(n5 if k==5 else n))
 res["beta_tables"]=dict(cases=2,exact=tb_ok)
 # X tables (tblx, W7): two calls a table; every word of the built rows, nothing else
 xp_ok,xc=True,0
 for c0_,a0_,c1_,a1_ in ((rng.randrange(p),rng.randrange(p),rng.randrange(p),rng.randrange(p)),(p-1,p-1,1,p-2),(0,rng.randrange(p),p-1,0)):
  mem={0x1880+8*w:0 for w in range(8)}
  sh=[(c0_*pow(2,o,p)%p,c1_*pow(2,o,p)%p)for o in E.W7O]
  for k,(x,y)in enumerate(sh):
   put(mem,0x1000+64*k,x)
   put(mem,0x1020+64*k,y)
  put(mem,0x1800,a0_)
  put(mem,0x1820,a1_)
  n,L=(5,3),(4,3)
  steps=sum(E.a64.a64_run(rprog,rlabels,"tblx",dict(x0=0x100000+256*c,x1=0x1000+256*c,x2=0x1800+128*c,x3=n[c],x23=L[c],x24=512,
     sp=SP,x30=E.SENTINEL),mem,ls64=([],[(0x100000,0x100000+512*n[0])]))for c in (0,1))
  want={0x100000+512*t+64*k+32*h:(sh[k][h]*t+((a0_,a1_)[h]if k==0 else 0))%p for k in range(7)for t in range(n[k>=4])for h in (0,1)}
  xp_ok&=steps==sum(25+26*(L[c]-1)+22*n[c]*L[c]for c in (0,1))and {a_ for a_ in mem if a_>=0x100000}=={a_+8*w for a_ in want
   for w in range(4)}and all(get(mem,a_)==v for a_,v in want.items())
  xc+=1
 res["x_tables"]=dict(cases=xc,exact=xp_ok)
 # TFP4S (CHAIN, 3.3): 97 instructions; H, L words bit for bit against tfp on random and carry-forcing c (l2 + h1, l3 + h2,
 # a0 + h3 at 2^64 + d), terms in 3.3's ranges (x30, the LR sentinel, is a limb)
 IC,tm,ok=pow(C0,-1,2**64),{},True
 for i in range(64):
  a=[rng.getrandbits(64)for _ in range(4)]
  a[1]|=M64*(i%8==1)
  d=rng.randrange(-2,2)
  for k in (2,3,0)[:3*(i%4>0)]:
   a[k]=(2**64-(a[k-1]*C0>>64)+d)*(IC if k else 1)%2**64
  c=sum(x<<64*k for k,x in enumerate(a))
  Ms=[rng.getrandbits(192)|(E.SENTINEL if k==2 else rng.getrandbits(64))<<192 for k in range(4)]
  regs=as_regs(c,"x5 x14 x6 x16".split())
  for M,q in zip(Ms,("x21 x22 x23 x24","x25 x26 x0 x2","x27 x28 x29 x30","x8 x9 x10 x11")):
   regs.update(as_regs(M,q.split()))
  o,n,_=site("t_tfp",regs)
  ok&=n==102 and [o[x]<<64|o[y]for x,y in zip("x17 x13 x19 x15".split(),"x16 x12 x18 x14".split())]==[E.tfp(c,M,tm)for M in Ms]
 res["tfp_chain"]=dict(cases=64,exact=ok,model=tm["sum_ok"]and all(lo<=tm[k][0]and tm[k][1]<=hi for k,(lo,hi)in E.TERMS.items()),ranges=tm)
 res["fs_canonical_exact"]=dict(cases=len(pairs),exact=fs)
 res["faddw_canonical_exact"]=dict(cases=len(pairs),exact=fa)
 res["fixups"]=fix
 res["fixups_ok"]=all(d["word_ok"]and d["count_ok"]and d["word_set_ok"]and d["fixups"]>0 for d in fix.values())
 weak=[0,1,C0,p-1,p,p+1,2**256-1]+[rng.getrandbits(256)for _ in range(12)]
 mul_in=[(a,b)for a in weak[:7]for b in weak[:7]]+[(a,rng.getrandbits(256))for a in weak[7:]]
 for _ in range(12):
  a=rng.randrange(1,p)
  mul_in.append((a,rng.randrange(1,2**64)*pow(a,-1,p)%p))
 sqr_in=weak+[p-rng.randrange(1,2**100)for _ in range(12)]+[2**256-rng.randrange(1,2**40)
                   for _ in range(6)]
 while len(sqr_in)<len(weak)+42:
  s_=rng.randrange(1,2**64)
  if pow(s_,(p-1)//2,p)==1:
   sqr_in.append(pow(s_,(p+1)//4,p))
 ok_m=ok_s=True
 for a,b in mul_in:
  mem={}
  put(mem,0x1000,a)
  put(mem,0x1100,b)
  E.a64.a64_run(fprog,flabels,"fe_mul_f",dict(x0=0x2000,x1=0x1000,x2=0x1100,sp=0x100000,x30=E.SENTINEL),mem)
  ok_m&=get(mem,0x2000)%p==a*b%p
 for a in sqr_in:
  mem={}
  put(mem,0x1000,a)
  E.a64.a64_run(fprog,flabels,"fe_sqr_f",dict(x0=0x2000,x1=0x1000,sp=0x100000,x30=E.SENTINEL),mem)
  ok_s&=get(mem,0x2000)%p==a*a%p
 def final_carry(v):
  t=(v%2**256)+(v>>256)*C0
  return (t%2**256)+(t>>256)*C0>=2**256
 res["fe_mul_f_exact"]=dict(cases=len(mul_in),exact=ok_m,final_carry_cases=sum(final_carry(a*b)for a,b in mul_in))
 res["fe_sqr_f_exact"]=dict(cases=len(sqr_in),exact=ok_s,final_carry_cases=sum(final_carry(a*a)for a in sqr_in))
 res["all_ok"]=res["fixups_ok"]and all(v.get("deviates_iff_predicate",True)and v.get("model",True)and v.get("exact",True)
         and v.get("predicate_holds",1)>0 for k,v in res.items()
         if isinstance(v,dict)and k!="fixups")
 return res

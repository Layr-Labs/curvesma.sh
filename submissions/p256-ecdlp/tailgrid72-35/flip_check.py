"""Executes flip_fwd / flip_bwd (proof.md section 3) on random points against the bit-exact model and exact labels."""
import random,collections,re
import a64,flip_witness as fw
from flip_witness import P,RR,M64,R,RI
BC=0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b
TB,XA,YA,BA,BINV,SP=1 <<51 ,0x100000 ,0x200000 ,0x400000 ,0x300000 ,0x700000
TBEND=TB+(1 <<56 )
WBASE,RBASE,SEND,SPORT=1 <<39 ,0x2000000 ,0x3000000 ,0x90000


def fwd_count(pairs,ny):
 return 152 *pairs+36 *ny+17


def bwd_count(pairs,ny):
 return 1933 *pairs//2 +11 *ny+24


def frow(c,k):
 return ((c*pow(2 ,64 *k,P)%P)<<256 )//P


def t192(v):
 return (v<<192 )//P


def tp_model(m,rows):
 a=0
 for k,r in enumerate(rows):
  mk=(m>>(64 *k))&M64
  w=[(r>>(64 *j))&M64 for j in range(4 )]
  a+=sum((mk*w[j])>>64 <<(64 *j)for j in range(3 ))+sum(((mk*w[j])&M64)<<(64 *j-64 )for j in (1 ,2 ,3 ))
 return a%(1 <<192 )


def put(mem,a,v,n=4 ):
 for k in range(n):
  mem[a+8 *k]=(v>>(64 *k))&M64


def get(mem,a,n=4 ):
 return sum(mem[a+8 *k]<<(64 *k)for k in range(n))


def table_word(cst,off):
 """Word at TB+off of the pair tables (proof.md section 3): pair pair P at 2^48 P, eight 2^45 slots: 0, 1 tau0 of
    both pairs interleaved (line i of row t: window 4 slot + i of pair 2P, then of 2P+1; window 7 does not exist and its
    line half is zero), 2..4 tau1..3 of 2P, 5..7 of 2P+1 (window slot 7 zero)."""
 PP,rest=divmod(off,1 <<48 )
 sl,rest=divmod(rest,1 <<45 )
 t,rest=divmod(rest,256 )
 k,w=divmod(rest,32 )
 if sl<2 :
  J,tau,k=2 *PP+k%2 ,0 ,4 *sl+k//2
 else:
  J,tau=(2 *PP,sl-1 )if sl<5 else (2 *PP+1 ,sl-4 )
 c,a0=cst[J][tau]
 return (fw.table_entry(c,a0,t,k)>>(64 *(w//8 )))&M64


def add(A,B):
 (x,y),(u,v)=A,B
 m=(v-y)*pow(u-x,-1 ,P)%P
 z=(m*m-x-u)%P
 return z,(m*(x-z)-y)%P


def y_entry(mem,a,x,y,n,node):
 """Y entry contents (proof.md section 3; written by yprep in production)."""
 xm=x*R%P
 put(mem,a,xm)
 put(mem,a+32 ,fw.redc(xm,xm)[0 ])
 for k in range(4 ):
  put(mem,a+64 +32 *k,frow(y*y*RI%P,k))
  put(mem,a+256 +32 *k,frow(2 *y%P,k))
 put(mem,a+192 ,t192(x),3 )
 mem[a+216 ]=(n>>28)<<45
 mem[a+224 ]=node|((n&((1<<28)-1))<<36)
 r4=frow(2 *y%P,4 )
 for c in range(8 ):
  put(mem,a+384 +32 *c,tp_model(c<<256 ,[0 ,0 ,0 ,0 ,r4]),3 )


def x_area(mem,a,x,y,j):
 for k in range(4 ):
  put(mem,a+32 *k,frow(y*y*RI%P,k))
 put(mem,a+128 ,t192(x),3 )
 mem[a+152 ]=(j>>1)<<36
 put(mem,a+160 ,x*R%P)


def consts(Xs,J):
 (x1,y1),(x2,y2)=Xs[2 *J],Xs[2 *J+1 ]
 return [((-(x1+x2))%P,x1*x2*R%P),((x2-x1)%P,0 ),(y1*RI%P,0 ),(y2*RI%P,0 )]


def run_case(prog,labels,seed,MP,NY,room=None):
 """room=None: ordinary windows; room=0 or >0: every window prefilled to 127 records (dropped / taken flushes)."""
 rng=random.Random(seed)

 def pt():
  while True:
   x=rng.randrange(P)
   w=(x*x*x-3 *x+BC)%P
   y=pow(w,(P+1 )//4 ,P)
   if y*y%P==w:
    return x,(y if rng.random()<0.5 else P-y)
 Xs=[pt()for _ in range(2 *MP)]
 Ys=[pt()for _ in range(NY)]
 source=(1<<36)-1;idbase=1024-2*MP;yibase=(1<<47)-NY
 cst={J:consts(Xs,J)for J in range(MP)}

 class Mem(dict):
  def __missing__(self,addr):
   if TB<=addr<TBEND:
    return table_word(cst,addr-TB)
   return 0
 mem=Mem()
 for c,f in enumerate(fw.fold_table()):
  put(mem,32 *c,f)
 for j,(x,y)in enumerate(Xs):
  x_area(mem,XA+256 *j,x,y,j+idbase)
 for n,(x,y)in enumerate(Ys):
  y_entry(mem,YA+1024 *n,x,y,n+yibase,source)
 LS=[(TB,TBEND),(BA,BA+64 *NY*MP),(YA,YA+1024 *NY)]
 rg={f'x{i}':rng.getrandbits(64 )for i in range(31 )}
 o={}
 nf=a64.a64_run(prog,labels,'flip_fwd',dict(rg,x0=NY,x1=YA,x2=BA,x3=TB,x4=64 *MP,sp=SP,x30=99 ),
 mem,out_regs=o,ls64=(LS,LS[1 :2 ]),limit=10 **8 )
 assert o['sp']==SP and o['x30']==99
 A,ent,ft=R,[],0
 for x,y in Ys:
  xm=x*R%P
  for J in range(MP):
   s_,c=fw.lookup_model(*cst[J][0 ],xm,fw.redc(xm,xm)[0 ])
   e=fw.fold_model(s_,c)
   assert e%P==(x-Xs[2 *J][0 ])*(x-Xs[2 *J+1 ][0 ])*R%P
   ent.append((e,A))
   A,t=fw.redc(A,e)
   ft+=t+int(s_+c*fw.DELTA>=RR)
 for q,(e,Ap)in enumerate(ent):
  assert get(mem,BA+64 *q)==e and get(mem,BA+64 *q+32 )==Ap
 assert get(mem,BA+64 *len(ent))==A and nf==fwd_count(len(ent),NY)+4 *ft
 B=fw.inverse13(A)
 put(mem,BINV,B)
 recs,errs,bt,q,d2short=[],[],0 ,len(ent),0
 for n in range(NY-1 ,-1 ,-1 ):
  x,y=Ys[n]
  for J in range(MP-1 ,-1 ,-1 ):
   q-=1
   e,Ap=ent[q]
   inv,t0=fw.redc(Ap,B)
   B,t1=fw.redc(e,B)
   d2short+=int(x*R%P>=Xs[2 *J+1 ][0 ]*R%P)
   i1,t2=fw.redc((x-Xs[2 *J+1 ][0 ])*R%P,inv)
   o2,c2=fw.lookup_model(cst[J][1 ][0 ],0 ,inv,i1)
   i2=fw.fold_model(o2,c2)
   bt+=t0+t1+t2+int(o2+c2*fw.DELTA>=RR)
   for wh,ii in ((0 ,i1),(1 ,i2)):
    X=Xs[2 *J+wh]
    assert ii*RI%P==pow(x-X[0 ],-1 ,P)
    w,t3=fw.redc(ii,ii)
    bt+=t3
    sig=[(frow(y*y*RI%P,k)+frow(X[1 ]*X[1 ]*RI%P,k))%RR for k in range(4 )]
    U=(tp_model(w,sig)-t192(x)-t192(X[0 ]))%(1 <<192 )
    no,nc=fw.lookup_model(cst[J][2 +wh][0 ],0 ,w)
    T2=tp_model(no+(nc<<256 ),[frow(2 *y%P,k)for k in range(5 )])
    for sg,L in ((0 ,(U-T2)%(1 <<192 )),(1 ,(U+T2)%(1 <<192 ))):
     L1,L2,L3=L>>128 ,(L>>64 )&M64,L&M64
     yi=n+yibase
     recs.append(((L1,((L2>>36 )|(L3<<28 ))&M64,source |((yi&((1<<28)-1))<<36),
     (J+idbase//2)<<36 |(yi>>28)<<45 |L3>>28 ),L2&((1 <<36 )-1 )))
     e_=(L-t192(add(X,(x,y if sg==0 else P-y))[0 ]))%(1 <<192 )
     errs.append(e_-(1 <<192 )if e_>=1 <<191 else e_)
 win={}
 for rec,d1 in recs:
  if d1 not in win:
   w=win[d1]=0x8000000 +len(win)*8192
   if room is None:
    mem[WBASE+d1*8 ]=w+32
   else:
    mem[WBASE+d1*8 ]=w+32 +32 *127
    mem[RBASE+d1*8 ]=SEND
    for t in range(127 ):
     put(mem,w+32 +32 *t,0xC0FFEE00 +t)
 pubs=[]
 o={}
 nb=a64.a64_run(prog,labels,'flip_bwd',dict(rg,x0=NY,x1=YA+1024 *(NY-1 ),x2=BA+64 *(len(ent)-1 ),
 x3=TB+(MP-1 )*(1 <<47 ),x4=XA+512 *(MP-1 ),x5=WBASE,x7=RBASE,
 x8=SEND+(65536 if room is None else room),x9=SPORT,x10=BINV,x11=64 *MP,sp=SP,x30=99 ),
 mem,out_regs=o,ls64=(LS,[]),limit=10 **8 ,dev=(SPORT,32 ,lambda k,v:pubs.append((k,v))))
 assert o['sp']==SP and o['x30']==99 and get(mem,BINV)==B==R
 F=len(pubs)//4
 D=0 if room is None else len(recs)-F
 assert nb+4 *d2short==bwd_count(len(ent),NY)+4 *bt+19 *F+16 *D
 by=collections.defaultdict(list)
 for rec,d1 in recs:
  by[d1].append(rec)
 for d1,rs in by.items():
  w=win[d1]
  base=w+32 if room is None else w+32 +32 *127
  if room is None:
   assert mem[WBASE+d1*8 ]==w+32 +32 *len(rs)and F==0
  else:
   assert len(rs)==1 and mem[WBASE+d1*8 ]==w+32
   assert mem[RBASE+d1*8 ]==SEND+4096 and all(mem[w+32 +32 *t]==0xC0FFEE00 +t for t in range(127 ))
  for t,rec in enumerate(rs):
   assert tuple(mem[base+32 *t+8 *k]for k in range(4 ))==rec
 if room:
  assert F==len(recs)and [k for k,v in pubs]==[0 ,8 ,16 ,24 ]*F
  assert sorted(pubs[i][1 ]for i in range(0 ,len(pubs),4 ))==sorted(win)
  assert all(pubs[i+1 ][1 ]==SEND and pubs[i+2 ][1 ]==win[pubs[i][1 ]]+32 and pubs[i+3 ][1 ]==128
  for i in range(0 ,len(pubs),4 ))
 elif room==0 :
  assert F==0
 return dict(fwd=nf,bwd=nb,tails=ft+bt,flushes=(F,D),records=len(recs),err=(min(errs),max(errs)))


def tail_pcs(prog):
 """First instruction of every carry-only REDC tail (B.CC over ADDS #1, SBCS, SBCS, SBC)."""
 return [i+1 for i,(mn,ops)in enumerate(prog[:-1 ])if mn=='b.cc'and prog[i+1 ][0 ]=='adds'and prog[i+1 ][1 ][2 ]=='#1']


YCF,YCB,YPP=123 ,561 ,945


def run_y(prog,labels,seed,n):
 """ychain_f, k_inv, ychain_b and yprep on n lanes: Y <- Y + S exactly, Y entries equal to y_entry(), exact counts."""
 rng=random.Random(seed)

 def pt():
  while True:
   x=rng.randrange(P)
   w=(x*x*x-3 *x+BC)%P
   y=pow(w,(P+1 )//4 ,P)
   if y*y%P==w:
    return x,(y if rng.random()<0.5 else P-y)
 CS,SC,YE,KB,BI,IS=0x10000 ,0x40000 ,0x80000 ,0x1000 ,0x3000 ,0x5000
 mem={}
 for c,f in enumerate(fw.fold_table()):
  put(mem,32 *c,f)
 mem.update({KB+a:v for a,v in fw.kblock().items()})
 lanes=[(pt(),pt())for _ in range(n)]
 for k,(Y,S)in enumerate(lanes):
  for o,v in ((0 ,Y[0 ]),(32 ,Y[1 ]),(64 ,S[0 ]),(96 ,S[1 ])):
   put(mem,CS+128 *k+o,v*R%P)
 tp=tail_pcs(prog)
 rg=lambda**kw:dict({f'x{i}':rng.getrandbits(64 )for i in range(31 )},sp=SP,x30=77 ,**kw)
 out,cnt={},{}
 for name,regs in (('ychain_f',dict(x0=n,x1=CS,x2=SC)),('k_inv',dict(x0=BI,x1=SC+64 *n,x2=IS)),
 ('ychain_b',dict(x0=n,x1=CS+128 *(n-1 ),x2=SC+64 *(n-1 ),x3=BI)),
 ('yprep',dict(x0=n,x1=CS,x2=YE,x3=KB,x4=5 <<45 ,x5=17 |(3 <<36 )))):
  vis={pc:0 for pc in tp}
  c=a64.a64_run(prog,labels,name,rg(**regs),mem,out_regs=out,visits=vis)
  assert out['sp']==SP and out['x30']==77
  cnt[name]=(c,sum(vis.values()))
 assert get(mem,BI)==R and cnt['k_inv'][0 ]<=24725
 assert cnt['ychain_f'][0 ]-4 *cnt['ychain_f'][1 ]==YCF*n+11
 assert cnt['ychain_b'][0 ]-4 *cnt['ychain_b'][1 ]==YCB*n+7
 assert cnt['yprep'][0 ]-4 *cnt['yprep'][1 ]==YPP*n+10
 for k,(Y,S)in enumerate(lanes):
  x,y=add(Y,S)
  assert get(mem,CS+128 *k)==x*R%P and get(mem,CS+128 *k+32 )==y*R%P
  ref={}
  y_entry(ref,0 ,x,y,0 ,0 )
  ref[216 ],ref[224 ]=5 <<45 ,17 +((3 +k)<<36 )
  assert all(mem[YE+1024 *k+o]==v for o,v in ref.items())
 return{k:v[0 ]for k,v in cnt.items()}


def run_tbld(prog,labels,seed,rows=64 ):
 """tbld for a standard table and for both interleaved tau0 tables of a pair pair (layout of table_word)."""
 rng=random.Random(seed)
 K,cst,out=0x1000 ,{J:[(rng.randrange(P),rng.randrange(P))]for J in (0 ,1 )},[]
 cases=[([cst[0 ][0 ][0 ]*pow(2 ,fw.WB*k,P)%P*(k<fw.NW)for k in range(8 )],[cst[0 ][0 ][1 ]]+[0 ]*7 ,None)]
 for sl in (0 ,1 ):
  cases.append(([cst[k%2 ][0 ][0 ]*pow(2 ,fw.WB*(4 *sl+k//2 ),P)%P*(4 *sl+k//2 <fw.NW)for k in range(8 )],
  [cst[0 ][0 ][1 ],cst[1 ][0 ][1 ]]*(sl==0 )+[0 ]*(8 -2 *(sl==0 )),sl))
 for m,s0,sl in cases:
  mem,base={},TB+((sl or 0 )<<45 )
  for k in range(8 ):
   put(mem,K+32 *k,m[k])
   put(mem,K+256 +32 *k,s0[k])
  n=a64.a64_run(prog,labels,'tbld',dict({f'x{i}':rng.getrandbits(64 )for i in range(31 )},x0=base,x1=rows,x2=K,x30=5 ,sp=SP),mem)
  assert all(get(mem,base+256 *t+32 *k)==(m[k]*t+s0[k])%P for t in range(rows)for k in range(8 ))
  if sl is None:
   assert all(get(mem,base+256 *t+32 *k)==fw.table_entry(*cst[0 ][0 ],t,k)for t in range(rows)for k in range(8 ))
  else:
   assert all(mem[base+8 *w]==table_word(cst,(sl<<45 )+8 *w)for w in range(32 *rows))
  assert n==140 *rows+87
  out.append(n)
 return out[0 ]


def census(seed,batches,NY=4 ,MP=16 ,fixed_x=False):
 """H_T support: exact tail model of every site on random P-256 points, per site and per pair position in the batch."""
 rng=random.Random(seed)

 def pt():
  while True:
   x=rng.randrange(P)
   w=(x*x*x-3 *x+BC)%P
   y=pow(w,(P+1 )//4 ,P)
   if y*y%P==w:
    return x,y
 npos=NY*MP
 rate={s:[0 ]*npos for s in ('F','I','B','D','W1','W2')}
 Xs=[pt()for _ in range(2 *MP)]
 for _ in range(batches):
  if not fixed_x:
   Xs=[pt()for _ in range(2 *MP)]
  Ys=[pt()for _ in range(NY)]
  A,ent=R,[]
  for x,y in Ys:
   xm=x*R%P
   for J in range(MP):
    (x1,_),(x2,_)=Xs[2 *J],Xs[2 *J+1 ]
    e=fw.fold_model(*fw.lookup_model((-(x1+x2))%P,x1*x2*R%P,xm,fw.redc(xm,xm)[0 ]))
    ent.append((e,A,xm,J))
    A,t=fw.redc(A,e)
    rate['F'][len(ent)-1 ]+=t
  B=fw.inverse13(A)
  for q in range(npos-1 ,-1 ,-1 ):
   e,Ap,xm,J=ent[q]
   inv,t0=fw.redc(Ap,B)
   B,t1=fw.redc(e,B)
   i1,t2=fw.redc((xm-Xs[2 *J+1 ][0 ]*R)%P,inv)
   i2=fw.fold_model(*fw.lookup_model((Xs[2 *J+1 ][0 ]-Xs[2 *J][0 ])%P,0 ,inv,i1))
   for s,t in (('I',t0),('B',t1),('D',t2),('W1',fw.redc(i1,i1)[1 ]),('W2',fw.redc(i2,i2)[1 ])):
    rate[s][q]+=t
  assert B==R
 return{s:[v/batches for v in r]for s,r in rate.items()}


BRANCH=('b','bl','cbz','cbnz','tbz','tbnz')
FORMS=frozenset(['add x,sp,#','adc x,x,x','adc x,x,xzr','adc x,xzr,x','adc x,xzr,xzr','adcs x,x,x','adcs x,x,xzr','add sp,sp,#','add x,x,#','add x,x,x','add x,x,x,lsl#','adds x,x,#','adds x,x,x','adds x,x,x,lsl#','adds x,xzr,x,lsl#','and x,x,#','and x,x,x','b L','b.cs L','b.cc L','b.eq L','b.hs L','b.lo L','b.ne L','bfi x,x,#,#','bl L','cbnz x,L','cbz x,L','ccmp x,x,#,eq','ccmp x,x,#,ne','cmn x,#','cmp x,#','cmp x,x','csel x,x,x,hs','csel x,x,x,lo','csel x,x,xzr,eq','cset x,hs','cset x,lo','dc zva,x','eor x,x,#','eor x,x,x','extr x,x,x,#','ld64b x,[x]','ldadd x,x,[x]','ldp x,x,[sp,#]','ldp x,x,[sp]','ldp x,x,[x,#]','ldp x,x,[x]','ldp x,x,[x],#','ldr x,[sp,#]','ldr x,[sp]','ldr x,[sp],#','ldr x,[x,#]','ldr x,[x,x]','ldr x,[x]','ldr x,[x],#','lsl x,x,#','lsr x,x,#','mov w,#','mov x,#','mov x,x','mov x,xzr','mul x,x,x','orr x,x,#','orr x,x,x','orr x,x,x,lsl#','orr x,x,x,lsr#','ret ','sbc x,x,x','sbc x,xzr,xzr','sbcs x,x,x','sbcs x,x,xzr','sbcs xzr,x,x','sbcs xzr,x,xzr','st64b x,[x]','stadd x,[x]','stlr x,[x]','stlr xzr,[x]','stp x,x,[sp,#]','stp x,x,[sp]','stp x,x,[x,#]','stp x,x,[x,#]!','stp x,x,[x]','stp x,x,[x],#','stp x,xzr,[x,#]','stp x,xzr,[x]','stp xzr,x,[x,#]','stp xzr,xzr,[x,#]','stp xzr,xzr,[x]','stp xzr,xzr,[x],#','str x,[sp,#]','str x,[sp,#]!','str x,[sp]','str x,[x,#]','str x,[x]','str x,[x],#','str xzr,[x,#]','str xzr,[x],#','sub sp,sp,#','sub x,x,#','sub x,x,x','sub x,x,x,lsl#','sub x,xzr,x','subs x,x,#','subs x,x,x','subs x,x,x,lsl#','swp x,x,[x]','tbnz x,#,L','tbz x,#,L','tst x,#','ubfx x,x,#,#','umulh x,x,x'])


def shape(mn,ops):
 """Operand shape of one parsed instruction: registers by class, immediates as #, branch targets as L."""
 def norm(o):
  o=re.sub(r'\s+','',o)
  o=re.sub(r'#-?(0x[0-9a-fA-F]+|\d+)','#',o)
  o=re.sub(r'\bx\d+\b','x',o)
  return re.sub(r'\bw\d+\b','w',o)
 ops=[norm(o)for o in ops]
 if mn in BRANCH or mn.startswith('b.'):
  ops=ops[:-1 ]+['L']
 return mn+' '+','.join(ops)


def smul(k,A):
 Z=None
 for b in bin(k)[2 :]:
  Z=None if Z is None else dbl(Z)
  if b=='1':
   Z=A if Z is None else add(Z,A)
 return Z


def dbl(A):
 x,y=A
 m=(3 *x*x-3 )*pow(2 *y,-1 ,P)%P
 z=(m*m-2 *x)%P
 return z,(m*(x-z)-y)%P


ORD=0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551
GP=(0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296 ,0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5 )


def run_group(prog,labels,seed,n=4 ,MP=4 ,rounds=2 ,rho0=5 ):
 """flip_group end to end: chains from coins (c, s), Y entries by yprep, both passes, both fixed inverses, chain steps;
    every record of every round equals the model built from Y=[c+u*s]P computed independently; exact total count."""
 rng=random.Random(seed)
 Xs=[]
 while len(Xs)<2 *MP:
  Xs.append(smul(rng.randrange(1 ,ORD),GP))
 coins=[(rng.randrange(1 ,ORD),rng.randrange(1 ,ORD))for _ in range(n)]
 cst={J:consts(Xs,J)for J in range(MP)}

 class Mem(dict):
  def __missing__(self,addr):
   if TB<=addr<TBEND:
    return table_word(cst,addr-TB)
   return 0
 mem=Mem()
 for c,f in enumerate(fw.fold_table()):
  put(mem,32 *c,f)
 KB_,CS,SC,CB,IS,LB,PB=0x1000 ,0x10000 ,0x40000 ,0x3000 ,0x5000 ,0x3100 ,0x6000
 mem.update({KB_+a:v for a,v in fw.kblock().items()})
 for j,(x,y)in enumerate(Xs):
  x_area(mem,XA+256 *j,x,y,j)
 for k,(c,s)in enumerate(coins):
  Y0,S0=smul(c+rho0*s,GP),smul(s,GP)
  for o,v in ((0 ,Y0[0 ]),(32 ,Y0[1 ]),(64 ,S0[0 ]),(96 ,S0[1 ])):
   put(mem,CS+128 *k+o,v*R%P)
 NODE=23
 words=[n,CS,SC,CB,IS,YA,KB_,0 ,0 ,BA,TB,64 *MP,LB,WBASE,RBASE,SEND+65536 ,SPORT,
 XA+512 *(MP-1 ),TB+(MP-1 )*(1 <<47 ),NODE,rho0,rounds]
 for i,w in enumerate(words):
  mem[PB+8 *i]=w
 recs,kinv=[],0
 for rd in range(rounds):
  rho=rho0+rd
  Ys=[smul(c+rho*s,GP)for c,s in coins]
  A,ent=R,[]
  for x,y in Ys:
   xm=x*R%P
   for J in range(MP):
    e=fw.fold_model(*fw.lookup_model(*cst[J][0 ],xm,fw.redc(xm,xm)[0 ]))
    ent.append((e,A))
    A=fw.redc(A,e)[0 ]
  B=fw.inverse13(A)
  q=len(ent)
  for nn in range(n-1 ,-1 ,-1 ):
   x,y=Ys[nn]
   for J in range(MP-1 ,-1 ,-1 ):
    q-=1
    e,Ap=ent[q]
    inv=fw.redc(Ap,B)[0 ]
    B=fw.redc(e,B)[0 ]
    i1=fw.redc((x-Xs[2 *J+1 ][0 ])*R%P,inv)[0 ]
    i2=fw.fold_model(*fw.lookup_model(cst[J][1 ][0 ],0 ,inv,i1))
    for wh,ii in ((0 ,i1),(1 ,i2)):
     X=Xs[2 *J+wh]
     w=fw.redc(ii,ii)[0 ]
     sig=[(frow(y*y*RI%P,k)+frow(X[1 ]*X[1 ]*RI%P,k))%RR for k in range(4 )]
     U=(tp_model(w,sig)-t192(x)-t192(X[0 ]))%(1 <<192 )
     no,nc=fw.lookup_model(cst[J][2 +wh][0 ],0 ,w)
     T2=tp_model(no+(nc<<256 ),[frow(2 *y%P,k)for k in range(5 )])
     yi=(rho<<10 )|nn
     for sg,L_ in ((0 ,(U-T2)%(1 <<192 )),(1 ,(U+T2)%(1 <<192 ))):
      L1,L2,L3=L_>>128 ,(L_>>64 )&M64,L_&M64
      recs.append(((L1,((L2>>36 )|(L3<<28 ))&M64,NODE|((yi&((1 <<28 )-1 ))<<36 ),
      J<<36 |(yi>>28 )<<45 |L3>>28 ),L2&((1 <<36 )-1 )))
  for v in (A,None):
   if v is None:
    v=R
    for (c,s),Y in zip(coins,Ys):
     v=fw.redc(v,(smul(s,GP)[0 ]-Y[0 ])*R%P)[0 ]
   m2={}
   put(m2,0x100 ,v)
   kinv+=a64.a64_run(prog,labels,'k_inv',dict({f'x{i}':0 for i in range(31 )},x0=0x200 ,x1=0x100 ,x2=0x300 ,sp=SP,x30=7 ),m2)
 win={}
 for rec,d1 in recs:
  if d1 not in win:
   win[d1]=0x8000000 +len(win)*8192
   mem[WBASE+d1*8 ]=win[d1]+32
 tp=tail_pcs(prog)
 vis={pc:0 for pc in tp}
 dvis={pcs[0 ]:0 for name,pcs in labels.items()if name.startswith('fbd2borrow')}
 vis.update(dvis)
 o={}
 LS=[(TB,2 *TB),(BA,BA+64 *n*MP),(YA,YA+1024 *n)]
 cnt=a64.a64_run(prog,labels,'flip_group',dict({f'x{i}':rng.getrandbits(64 )for i in range(31 )},x0=PB,sp=SP,x30=99 ),
 mem,out_regs=o,visits=vis,ls64=(LS,LS[1 :2 ]),limit=10 **9 )
 assert o['sp']==SP and o['x30']==99 and get(mem,CB)==R and get(mem,LB)==R
 by=collections.defaultdict(list)
 for rec,d1 in recs:
  by[d1].append(rec)
 for d1,rs in by.items():
  assert mem[WBASE+d1*8 ]==win[d1]+32 +32 *len(rs)
  assert all(tuple(mem[win[d1]+32 +32 *t+8 *k]for k in range(4 ))==rec for t,rec in enumerate(rs))
 for k,(c,s)in enumerate(coins):
  Z=smul(c+(rho0+rounds)*s,GP)
  assert get(mem,CS+128 *k)==Z[0 ]*R%P and get(mem,CS+128 *k+32 )==Z[1 ]*R%P
 per_round=fwd_count(n*MP,n)+bwd_count(n*MP,n)+(YCF+YCB+YPP)*n+11 +7 +10 +fw.DRIVER_PER_ROUND
 assert cnt-4 *sum(vis[pc]for pc in tp)+4 *(rounds*n*MP-sum(vis[pc]for pc in dvis))-kinv==rounds*per_round+fw.DRIVER_FIXED
 return dict(rounds=rounds,lanes=n,pairs=MP,records=len(recs),count=cnt,driver_per_round=fw.DRIVER_PER_ROUND)


def run_fold():
 L=fw.fold_lines(['x16','x17','x18','x19'],'x1',['x4','x5','x6','x7'],'fold_end')
 pp,ll=a64.a64_load('check:\n'+'\n'.join(x if x.endswith(':')else ' '+x for x in L)+'\n ret\n')
 rng=random.Random(20261106 );counts={}
 for c in range(9 ):
  edge=(RR-c*fw.DELTA)%RR
  for v in [0 ,1 ,RR-1 ,edge,(edge-1 )%RR,(edge+1 )%RR]+[rng.randrange(RR)for _ in range(64 )]:
   mm={};[put(mm,32 *k,k*fw.DELTA)for k in range(9 )];oo={}
   n=a64.a64_run(pp,ll,'check',dict(x0=0xffffffff ,x2=0xffffffff00000001 ,x1=c,x16=v&M64,x17=(v>>64 )&M64,x18=(v>>128 )&M64,x19=v>>192 ,x30=777 ),mm,out_regs=oo)
   z=sum(oo['x%d'%j]<<(64 *(j-16 ))for j in range(16 ,20 ));carry=int(v+c*fw.DELTA>=RR)
   assert z==fw.fold_model(v,c)and z%P==(v+c*RR)%P and n==9 +4 *carry
   assert not carry or z<9 *fw.DELTA
   counts[n]=counts.get(n,0 )+1
 assert set(counts)=={9 ,13 }and 9 *fw.DELTA<P
 return dict(cases=sum(counts.values()),counts=counts,all_checks=True)

"""pair_bwd_f against its bit-exact model (frac_witness) on boundary fixtures, and the full-window flush paths."""
import random
import a64,range_witness,frac_witness as fw
P=fw.P;RR=fw.RR;R=fw.R;RI=fw.RI;MASK=fw.MASK;M35=(1<<35)-1;D=RR-P
FB=0x7000000000;XA=0x50000
def put(m,p,v):
 for i in range(4):m[p+8*i]=(v>>(64*i))&MASK
def get(m,p):return sum(m[p+8*i]<<(64*i)for i in range(4))
def _row(nf,nfl,rng,build,descriptors,xx,xy,ys,sp):
 """Forward, inverse (exact), Y/X fraction tables; returns memory, segments, order, expected records."""
 tr,seq,nodes=build(ys);ds=descriptors(tr,seq,xx,xy)
 xp=0x30000;bp=0x40000;scratch=0x60000
 m={};put(m,xp,xx);put(m,xp+32,xy);put(m,bp-32,R)
 for nd in nodes:
  for j,(x,y,idx)in enumerate(nd['items']):
   put(m,nd['base']+72*j,x);put(m,nd['base']+72*j+32,y);m[nd['base']+72*j+64]=idx<<36
 fw.install(m,FB,XA,xx*RI%P,xy*RI%P,ys)
 def regs(**kw):
  z={f'x{i}':rng.getrandbits(64) for i in range(31)};z.update(sp=sp,x30=12345);z.update(kw);return z
 order=[];segments=[];off=0
 for nd,a,b,dx,yl,carry in ds:
  cnt=b-a;ptr=nd['base']+72*a;out={}
  a64.a64_run(nf,nfl,'pair_fwd_lo' if dx else 'pair_fwd_hi',regs(x0=cnt,x1=ptr,x2=bp+64*off,x3=xp,x4=scratch),m,limit=300000,out_regs=out)
  assert out['x0']==1 and out['sp']==sp and out['x30']==12345
  segments.append((off,cnt,ptr));order.extend(nd['items'][a:b]);off+=cnt
 n=len(order);B=pow(get(m,bp+64*(n-1)+32)*RI%P,-1,P)*R%P;put(m,scratch,B)
 A=[R]+[get(m,bp+64*j+32)for j in range(n)];dd=[get(m,bp+64*j)for j in range(n)]
 return m,segments,order,(bp,scratch,regs),(B,A,dd)
def check(root,build,descriptors):
 nf,nfl=a64.a64_load((root/'p256_fields.s').read_text()+range_witness.generated_assembly(root))
 rng=random.Random(61202);cases=outputs=0;errs=[];tails=0;edge=0
 for n in[1,2,7,16]:
  for trial in range(16):
   xx=[0,1,D,D+1,RR//2,P//2,P-1,rng.randrange(P)][trial%8]%P;xy=rng.randrange(1,P//2+1)
   pool=[z for z in[0,1,P-1,xx-1,xx+1,RR-xx-1,RR-xx,RR-xx+1]if 0<=z<P and z!=xx]
   ys=[(pool[j%len(pool)]if trial<7 else rng.randrange(P),1+rng.randrange(P//2),j)for j in range(n)]
   ys=[((x+1)%P if x==xx else x,y,j)for x,y,j in ys]
   if trial>=8:# x(X+Y) placed on label-digit boundaries: T(x) mod 2^16 in the bad residues, and T(x) near 0, 2^192
    ax=rng.randrange(P);by=rng.randrange(1,P);xx=ax*R%P;xy=by*R%P;ys=[]
    for j in range(n):
     t=[0,1,8,(1<<16)-1,(1<<16)-9,(1<<16)+16,(1<<192)-1,(1<<191)][(trial-8+j)%8]+(rng.getrandbits(160)<<32 if trial%2 else 0)
     target=(t%(1<<192))*P>>192;mm=rng.randrange(1,P);cx=(mm*mm-ax-target)%P
     while cx==ax:mm=(mm+1)%P;cx=(mm*mm-ax-target)%P
     dy=(by+mm*(cx-ax))%P;ys.append((cx*R%P,dy*R%P,j));edge+=1
   sp=0x80000
   m,segments,order,(bp,scratch,regs),(B,A,dd)=_row(nf,nfl,rng,build,descriptors,xx,xy,ys,sp)
   snap={k:v for k,v in m.items()}
   exp,tl,er=fw.expected(B,A,dd,order,xx*RI%P,xy*RI%P,3|(5<<35));errs+=er
   WBASE=1<<38;RBASE=0x2000000;SEND=0x3000000;SPORT=0x90000;W2=3|(5<<35);win_of={}
   for rec,d1 in exp:
    if d1 not in win_of:win_of[d1]=0x8000000+len(win_of)*8192;m[WBASE+d1*8]=win_of[d1]+32
   for off,cnt,ptr in reversed(segments):
    out={}
    ct=a64.a64_run(nf,nfl,'pair_bwd_f',regs(x0=cnt,x1=ptr+72*cnt+64,x2=bp+64*(off+cnt-1),x3=XA,x4=scratch,x5=WBASE,x6=W2,x7=RBASE,x8=SEND+65536,x9=SPORT,x10=FB),m,limit=600000,out_regs=out,ls64=fw.ls64_regions(FB,XA))
    assert out['sp']==sp and out['x30']==12345
    k=sum(tl[off:off+cnt]);tails+=k
    assert ct==fw.tail_free(cnt)+4*k and ct<=2173*cnt//4+34,(n,trial,ct,cnt)
   assert get(m,scratch)==R and all(m[k]==v for k,v in snap.items()if not(scratch<=k<scratch+32))
   bywin={}
   for rec,d1 in exp:bywin.setdefault(d1,[]).append(rec)
   for d1,recs in bywin.items():
    win=win_of[d1];assert m[WBASE+d1*8]==win+32+32*len(recs)
    for t,rec in enumerate(recs):
     assert tuple(m[win+32+32*t+8*q]for q in range(4))==rec;outputs+=1
   cases+=1
 assert min(errs)>=-17 and max(errs)<=9 and edge>0
 return dict(cases=cases,records=outputs,boundary_points=edge,label_error=[min(errs),max(errs)],tails_taken=tails,count_exact='529n+21-14floor(n/8)+4*tails',scope='bit-exact model and label bound on boundary fixtures; finite evidence')
def full_window(root,build,descriptors):
 """Pre-filled windows force every emission to fire flush_publish_b (taken and dropped trials)."""
 nf,nfl=a64.a64_load((root/'p256_fields.s').read_text()+range_witness.generated_assembly(root))
 rng=random.Random(975377);firings=0;cases=0
 for trial,room in((0,None),(1,0)):
  n=70;xx=rng.randrange(1,P);xy=rng.randrange(1,P//2+1)
  ys=[(rng.randrange(P),1+rng.randrange(P//2),j)for j in range(n)]
  ys=[((x+1)%P if x==xx else x,y,j)for x,y,j in ys]
  sp=0x80000
  m2,segments,order,(bp,scratch,regs),(B,A,dd)=_row(nf,nfl,rng,build,descriptors,xx,xy,ys,sp)
  W2=3|(5<<35);exp,tl,_=fw.expected(B,A,dd,order,xx*RI%P,xy*RI%P,W2)
  WBASE=1<<38;RBASE=0x2000000;SEND=0x3000000;SPORT=0x90000;win_of={}
  room=4096*(2*n+8)if room is None else room
  for rec,d1 in exp:
   if d1 not in win_of:
    win=0x8000000+len(win_of)*8192;win_of[d1]=win
    m2[WBASE+d1*8]=win+32+32*127
    for t in range(127):
     m2[win+32+32*t]=0xC0FFEE00+t;m2[win+32+32*t+8]=trial;m2[win+32+32*t+16]=W2;m2[win+32+32*t+24]=t
    m2[RBASE+d1*8]=SEND
  assert len(win_of)==2*n
  pubs=[]
  def capture(k,v):pubs.append((k,v))
  for off,cnt,ptr in reversed(segments):
   out={};npub=len(pubs)
   ct=a64.a64_run(nf,nfl,'pair_bwd_f',regs(x0=cnt,x1=ptr+72*cnt+64,x2=bp+64*(off+cnt-1),x3=XA,x4=scratch,x5=WBASE,x6=W2,x7=RBASE,x8=SEND+room,x9=SPORT,x10=FB),m2,limit=900000,out_regs=out,dev=(SPORT,32,capture),ls64=fw.ls64_regions(FB,XA))
   assert out['sp']==sp and out['x30']==12345
   segpubs=pubs[npub:];assert len(segpubs)%4==0
   F=len(segpubs)//4;Dn=2*cnt-F
   assert ct==fw.tail_free(cnt)+4*sum(tl[off:off+cnt])+19*F+16*Dn and ct<=2173*cnt//4+34+19*F+16*Dn
   firings+=F+Dn
  if trial==0:
   assert len(pubs)==8*n
   seen={}
   for i in range(0,len(pubs),4):
    d,rb,wb,cnt_=[pubs[i+k][1]for k in range(4)]
    assert[pubs[i+k][0]for k in range(4)]==[0,8,16,24]
    assert cnt_==128 and wb==win_of[d]+32 and rb==SEND
    seen[d]=seen.get(d,0)+1
   assert sorted(seen)==sorted(win_of)and all(v==1 for v in seen.values())
  else:
   assert pubs==[]
  for d1,win in win_of.items():
   assert m2[WBASE+d1*8]==win+32 and m2[RBASE+d1*8]==SEND+4096,(trial,d1)
  for rec,d1 in exp:
   win=win_of[d1];base=win+32+32*127
   assert tuple(m2[base+8*q]for q in range(4))==rec
   for t in range(127):
    assert m2[win+32+32*t]==0xC0FFEE00+t and m2[win+32+32*t+24]==t
  cases+=1
 return dict(cases=cases,firings=firings,taken_and_dropped=True,scope='pair_bwd_f full-window flush integration on synthetic fixtures')

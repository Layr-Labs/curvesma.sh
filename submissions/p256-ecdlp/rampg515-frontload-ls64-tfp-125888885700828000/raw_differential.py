"""Fused-vs-rolled differential on boundary fixtures; rolled base is v47-reviewed and hash-pinned."""
import random
import a64,range_witness
P=0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
RR=1<<256;R=RR%P;RI=pow(R,-1,P);MASK=(1<<64)-1;M35=(1<<35)-1;D=RR-P;L=RR-(1<<224)
def check(root,build,descriptors):
 new,nl=a64.a64_load((root/'p256_fields.s').read_text()+range_witness.base_generated_assembly(root))
 nf,nfl=a64.a64_load((root/'p256_fields.s').read_text()+range_witness.generated_assembly(root))
 rng=random.Random(61202);cases=outputs=0;variants=set();noncanonical=0;fused_records=0
 def put(m,p,v):
  for i in range(4):m[p+8*i]=(v>>(64*i))&MASK
 def get(m,p):return sum(m[p+8*i]<<(64*i)for i in range(4))
 for n in[1,2,7,16]:
  for trial in range(16):
   xx=[0,1,D,D+1,RR//2,P//2,P-1,rng.randrange(P)][trial%8];xy=rng.randrange(1,P//2+1)
   pool=[0,1,P-1,xx-1,xx+1,RR-xx-1,RR-xx,RR-xx+1]
   pool=[z for z in pool if 0<=z<P and z!=xx]
   ys=[(pool[j%len(pool)]if trial<7 else rng.randrange(P),1+rng.randrange(P//2),j)for j in range(n)]
   ys=[((x+1)%P if x==xx else x,y,j)for x,y,j in ys]
   if trial>=8:
    ax=rng.randrange(P);by=rng.randrange(1,P);xx=ax*R%P;xy=by*R%P
    if xy>P//2:xy=P-xy
    ys=[]
    targets=[0,1,D-1,D,D+1,P-D-1,P-D+1,P-1]
    for j in range(n):
     target=targets[(trial-8+j)%len(targets)];m=rng.randrange(1,P)
     cx=(m*m-ax-target*RI)%P
     while cx==ax:
      m=(m+1)%P;cx=(m*m-ax-target*RI)%P
     dy=(by+m*(cx-ax))%P;ym=dy*R%P
     if ym>P//2:ym=P-ym
     ys.append((cx*R%P,ym,j))
   tr,seq,nodes=build(ys);ds=descriptors(tr,seq,xx,xy)
   xp=0x30000;bp=0x40000;scratch=0x60000;sp=0x80000
   m={};put(m,xp,xx);put(m,xp+32,xy);put(m,bp-32,R)
   for nd in nodes:
    for j,(x,y,idx)in enumerate(nd['items']):
     put(m,nd['base']+72*j,x);put(m,nd['base']+72*j+32,y);m[nd['base']+72*j+64]=idx<<36
   snap=dict(m)
   order=[];segments=[];off=0
   def regs(**kw):
    z={f'x{i}':rng.getrandbits(64) for i in range(31)};z.update(sp=sp,x30=12345);z.update(kw);return z
   for nd,a,b,dx,yl,carry in ds:
    cnt=b-a;ptr=nd['base']+72*a;out={}
    a64.a64_run(new,nl,'pair_fwd_lo' if dx else 'pair_fwd_hi',regs(x0=cnt,x1=ptr,x2=bp+64*off,x3=xp,x4=scratch),m,limit=300000,out_regs=out)
    assert out['x0']==1 and out['sp']==sp and out['x30']==12345
    segments.append((off,cnt,ptr,yl,carry));order.extend(nd['items'][a:b]);off+=cnt
   put(m,scratch,pow(get(m,bp+64*(n-1)+32)*RI%P,-1,P)*R%P)
   dead={scratch+off:0x9e3779b97f4a7c15^off for off in range(32,192,8)}
   m.update(dead)
   m2=dict(m)
   for off,cnt,ptr,yl,carry in reversed(segments):
    variants.add((yl,carry));out={}
    a64.a64_run(new,nl,f'pair_bwd_y{yl}_c{carry}',regs(x0=cnt,x1=ptr+72*(cnt-1),x2=bp+64*(off+cnt-1),x3=xp,x4=scratch),m,limit=300000,out_regs=out)
    assert out['sp']==sp and out['x30']==12345
   assert get(m,scratch)==R
   WBASE=1<<38;RBASE=0x2000000;SEND=0x3000000;SPORT=0x90000;W2=3|(5<<35);win_of={}
   for j in range(2*n):
    d1=get(m,bp+32*j)&M35
    if d1 not in win_of:win_of[d1]=0x8000000+len(win_of)*8192;m2[WBASE+d1*8]=win_of[d1]+32
   for off,cnt,ptr,yl,carry in reversed(segments):
    out={}
    ct=a64.a64_run(nf,nfl,f'pair_bwd_y{yl}_c{carry}',regs(x0=cnt,x1=ptr+72*(cnt-1),x2=bp+64*(off+cnt-1),x3=xp,x4=scratch,x5=WBASE,x6=W2,x7=RBASE,x8=SEND+65536,x9=SPORT),m2,limit=600000,out_regs=out)
    assert ct<=2489*cnt//4+33,(trial,n,ct,cnt);assert out['sp']==sp and out['x30']==12345
   assert get(m2,scratch)==R
   assert all(m2[k]==v for k,v in snap.items())
   assert all(m2[k]==v for k,v in dead.items())
   emitted=[]
   for off,cnt,ptr,yl,carry in reversed(segments):
    for j in range(off+cnt-1,off-1,-1):
     _,_,idx=order[j];emitted.append((get(m,bp+64*j),0,idx));emitted.append((get(m,bp+64*j+32),1,idx))
   bywin={}
   for(key,sign,idx)in emitted:bywin.setdefault(key&M35,[]).append((key,sign,idx))
   for d1,recs in bywin.items():
    win=win_of[d1]
    assert m2[WBASE+d1*8]==win+32+32*len(recs)
    for t,(key,sign,idx)in enumerate(recs):
     base=win+32+32*t
     assert m2[base]==((key>>128)&MASK)and m2[base+8]==(key>>192)
     assert m2[base+16]==W2 and m2[base+24]==((key>>64)&M35)|(sign<<35)|(idx<<36);fused_records+=1;outputs+=1
     if trial>=8:noncanonical+=1
   cases+=1
 assert variants=={(0,0),(0,1),(1,0),(1,1)}
 assert noncanonical>0
 return dict(cases=cases,raw_outputs=outputs,noncanon_observed=noncanonical,four_variants=True,fused_records=fused_records,source_mem_ok=True,dead_scratch_canaries=True,scope='fused-vs-rolled window differential on boundary fixtures, not universal proof or success claim')
def full_window(root,build,descriptors):
 """Fused full-window flush integration: pre-filled windows force every emission to fire
 flush_publish_b (taken and dropped trials). Not universal proof or a success claim."""
 from pathlib import Path
 root=Path(root)
 new,nl=a64.a64_load((root/'p256_fields.s').read_text()+range_witness.base_generated_assembly(root))
 nf,nfl=a64.a64_load((root/'p256_fields.s').read_text()+range_witness.generated_assembly(root))
 rng=random.Random(975377);variants=set();firings=0;cases=0
 def put(m,p,v):
  for i in range(4):m[p+8*i]=(v>>(64*i))&MASK
 def get(m,p):return sum(m[p+8*i]<<(64*i)for i in range(4))
 for trial,room in((0,None),(1,0)):
  n=70;xx=rng.randrange(1,P);xy=rng.randrange(1,P//2+1)
  ys=[(rng.randrange(P),1+rng.randrange(P//2),j)for j in range(n)]
  ys=[((x+1)%P if x==xx else x,y,j)for x,y,j in ys]
  tr,seq,nodes=build(ys);ds=descriptors(tr,seq,xx,xy)
  xp=0x30000;bp=0x40000;scratch=0x60000;sp=0x80000
  m={};put(m,xp,xx);put(m,xp+32,xy);put(m,bp-32,R)
  for nd in nodes:
   for j,(x,y,idx)in enumerate(nd['items']):
    put(m,nd['base']+72*j,x);put(m,nd['base']+72*j+32,y);m[nd['base']+72*j+64]=idx<<36
  def regs(**kw):
   z={f'x{i}':rng.getrandbits(64) for i in range(31)};z.update(sp=sp,x30=12345);z.update(kw);return z
  order=[];segments=[];off=0
  for nd,a,b,dx,yl,carry in ds:
   cnt=b-a;ptr=nd['base']+72*a;out={}
   a64.a64_run(new,nl,'pair_fwd_lo' if dx else 'pair_fwd_hi',regs(x0=cnt,x1=ptr,x2=bp+64*off,x3=xp,x4=scratch),m,limit=300000,out_regs=out)
   segments.append((off,cnt,ptr,yl,carry));order.extend(nd['items'][a:b]);off+=cnt
  put(m,scratch,pow(get(m,bp+64*(n-1)+32)*RI%P,-1,P)*R%P)
  m2=dict(m)
  for off,cnt,ptr,yl,carry in reversed(segments):
   variants.add((yl,carry));out={}
   a64.a64_run(new,nl,f'pair_bwd_y{yl}_c{carry}',regs(x0=cnt,x1=ptr+72*(cnt-1),x2=bp+64*(off+cnt-1),x3=xp,x4=scratch),m,limit=300000,out_regs=out)
  assert get(m,scratch)==R
  WBASE=1<<38;RBASE=0x2000000;SEND=0x3000000;SPORT=0x90000;W2=3|(5<<35);win_of={}
  room=4096*(2*n+8)if room is None else room
  for j in range(2*n):
   d1=get(m,bp+32*j)&M35
   if d1 not in win_of:
    win=0x8000000+len(win_of)*8192;win_of[d1]=win
    m2[WBASE+d1*8]=win+32+32*127
    for t in range(127):
     m2[win+32+32*t]=0xC0FFEE00+t;m2[win+32+32*t+8]=trial;m2[win+32+32*t+16]=W2;m2[win+32+32*t+24]=t
    m2[RBASE+d1*8]=SEND
  assert len(win_of)==2*n
  pubs=[]
  def capture(k,v):pubs.append((k,v))
  for off,cnt,ptr,yl,carry in reversed(segments):
   out={};npub=len(pubs)
   ct=a64.a64_run(nf,nfl,f'pair_bwd_y{yl}_c{carry}',regs(x0=cnt,x1=ptr+72*(cnt-1),x2=bp+64*(off+cnt-1),x3=xp,x4=scratch,x5=WBASE,x6=W2,x7=RBASE,x8=SEND+room,x9=SPORT),m2,limit=900000,out_regs=out,dev=(SPORT,32,capture))
   assert out['sp']==sp and out['x30']==12345
   segpubs=pubs[npub:]
   assert len(segpubs)%4==0
   F=len(segpubs)//4;D=2*cnt-F
   assert ct<=2489*cnt//4+33+19*F+16*D
   firings+=F+D
  if trial==0:
   assert len(pubs)==8*n
   seen={}
   for i in range(0,len(pubs),4):
    d,rb,wb,cnt_=[pubs[i+k][1]for k in range(4)]
    assert[pubs[i+k][0]for k in range(4)]==[0,8,16,24]
    assert cnt_==128 and wb==win_of[d]+32 and rb==SEND
    seen[d]=seen.get(d,0)+1
   assert sorted(seen)==sorted(win_of)and all(v==1 for v in seen.values())
   for d1,win in win_of.items():
    assert m2[WBASE+d1*8]==win+32 and m2[RBASE+d1*8]==SEND+4096,(trial,d1)
  else:
   assert pubs==[]
   for d1,win in win_of.items():
    assert m2[WBASE+d1*8]==win+32 and m2[RBASE+d1*8]==SEND+4096,(trial,d1)
  emitted=[]
  for off,cnt,ptr,yl,carry in reversed(segments):
   for j in range(off+cnt-1,off-1,-1):
    _,_,idx=order[j];emitted.append((get(m,bp+64*j),0,idx));emitted.append((get(m,bp+64*j+32),1,idx))
  for(key,sign,idx)in emitted:
   d1=key&M35;win=win_of[d1];base=win+32+32*127
   assert m2[base]==((key>>128)&MASK)and m2[base+8]==(key>>192)
   assert m2[base+16]==W2 and m2[base+24]==((key>>64)&M35)|(sign<<35)|(idx<<36)
   for t in range(127):
    assert m2[win+32+32*t]==0xC0FFEE00+t and m2[win+32+32*t+24]==t
  cases+=1
 assert variants=={(0,0),(0,1),(1,0),(1,1)}
 return dict(cases=cases,firings=firings,four_variants=True,taken_and_dropped=True,scope='fused full-window flush integration on synthetic fixtures')

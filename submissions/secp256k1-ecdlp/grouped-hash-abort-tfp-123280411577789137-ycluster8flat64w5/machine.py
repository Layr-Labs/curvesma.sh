import os,random;from collections import defaultdict as dd;K1,K2B=256,512;B1,B2=64*K1,64*K2B;UP1=18432*(1<<40);UP=(1<<64)-65536*(1<<40);LO1,H1=1024*(1<<40),1<<40;LC2,RC2=1<<45,1<<46;MB2=(1<<48)-128;PW=(1<<48)-168;CNT=(1<<56)-8192;CNTE=(1<<56)-6144;DONE=(1<<48)-4096;ST=(1<<56)-160;HT=(1<<64)-(1<<56);NEXT_BASE=(1<<64)-57003*(1<<40)
def s64(v):return v-(1<<64)if(v&(1<<63))else(v&((1<<64)-1))
def asr1(v):return(s64(v)>>1)&((1<<64)-1)
LINK=(1<<56)-(1<<40);LINK_MSG=LINK+0x8000;WINDOW=(LINK,LINK+0x8000);LINK_GEN=34<<40;LINK_GEN_MSG=LINK_GEN+0x1000;PORTY_BASE=LINK_GEN+0xB000;STALE=0x5A5A5A5A5A5A5A5A;M41=2**41-1;DM=(2**64-1)^M41;H2S_REGS,DLVH_REGS,cb={},{},[0]
def blk1(d):return LO1+B1*d if d<H1 else UP1+B1*(d-H1)
def blk2(d):return UP+B2*d
LC1_TABLES=[(0,8*H1,LO1+32,B1),(8*H1,16*H1,UP1+32,B1),(LO1,LO1+B1*H1,0,0),(UP1,UP1+B1*H1,0,0)]
def dest(word):return word>>3&M41
def dest1(word):return word>>3
class Port1:
 def __init__(J,mem):
  J.mem,J.blocks,J.pend=mem,[],[];J.dev=(0,1<<64,J.store,J.read,lambda:not J.pend)
 def store(J,a,v):
  if a%B1==0 and(LO1<=a<LO1+B1*H1 or UP1<=a<UP1+B1*H1):J.pend.append((a,v))
 def read(J):
  for a,v in J.pend:J.blocks.append([v]+[J.mem.get(a+8*w,0)for w in(1,2,3)]+[J.mem.get(a+32+8*w,STALE)for w in range(4*K1)])
  J.pend=[]
class Port2:
 def __init__(J,mem,mb=MB2):
  J.mem,J.msgs,J.blocks,J.mb,J.pend=mem,[],[],mb,[];J.dev=(0,1<<64,J.store,J.read,lambda:not J.pend)
 def store(J,a,v):
  mb=J.mb
  if a==mb+40 or(not mb<=a<mb+48 and a>=UP and(a-UP)%B2==0):J.pend.append((a,v))
 def read(J):
  for a,v in J.pend:
   if a==J.mb+40:J.msgs.append([J.mem.get(J.mb+8*w,0)for w in range(5)]+[v])
   else:J.blocks.append([v,J.mem.get(a+8,0),J.mem.get(a+16,0),J.mem.get(a+24,0)]+[J.mem.get(a+32+8*w,STALE)for w in range(4*K2B)])
  J.pend=[]
class LinkView:
 def __init__(J,mem,msg,writes,lbase=LINK):J.mem,J.msg,J.writes,J.lbase=mem,msg,writes,lbase;J.base=lbase+32 if lbase==LINK and len(msg)>6 else lbase
 def __getitem__(J,a):return J.msg[(a-J.base)//8]if J.base<=a<J.base+8*len(J.msg)else J.mem[a]
 def get(J,a,d=None):return J[a]if J.base<=a<J.base+8*len(J.msg)else J.mem.get(a,d)
 def __setitem__(J,a,v):
  assert not J.lbase<=a<J.lbase+0x8000,"write into link aperture"
  if a not in(CNT,CNT+2048,CNTE):J.writes.append(a)
  J.mem[a]=v
DLVH_PORTS=tuple(LINK+64*k for k in range(1,9));H2S_PORTS=tuple(LINK+64*k for k in range(1,7))
def dlvh0(E,ad,labels,mem,cnt):r0={};E.a64.a64_run(ad,labels,"dlvh0",dict(x29=cnt,x30=E.SNT),mem,out_regs=r0);assert tuple(r0[r]for r in("x0","x1","x7","x8","x24","x25","x26","x27"))==DLVH_PORTS and r0["x28"]==1;return r0
def deliver(E,ad,labels,mem,msg,calls,cert,writes=None,fused=False,r2=None,formula=None,cnt=CNT):
 e,o="dlv"if len(msg)==6 else"dlvh",{};lb=LINK_MSG if e=="dlv"else LINK;mem.setdefault(cnt,0);c0=mem[cnt];g8=DLVH_REGS.get(id(mem))or dlvh0(E,ad,labels,mem,cnt)if e=="dlvh"else dict(x0=lb,x28=1);regs=dict(g8,x19=STALE,x29=cnt,x30=E.SNT);n=E.a64.a64_run(ad,labels,e,regs,LinkView(mem,msg,[]if writes is None else writes,lb),out_regs=o,ls64=([WINDOW]+([r2]if r2 else[]),[r2]if r2 else[]))
 if e=="dlvh":DLVH_REGS[id(mem)]=o;assert tuple(o[r]for r in("x0","x1","x7","x8","x24","x25","x26","x27"))==DLVH_PORTS
 ok=o["x28"]==1 and o["x29"]==cnt and mem[cnt]==c0+1;calls.append((e,n if ok else-1,cert.DLV+cert.RCV if e=="dlv"else formula))
def h2s_state(o):return o["x14"]<<64|o["x13"],o["x5"]<<64|o["x15"]
def h2s0(E,ad,labels,nd,end,af):r0={};E.a64.a64_run(ad,labels,"h2s0",dict(x12=end,x0=af>>64,x8=af%2**64,x30=E.SNT),nd.mem,out_regs=r0);assert(r0["x3"],r0["x10"],r0["x18"],r0["x12"],r0["x29"],r0["x6"],r0["x28"])==(PW,LC2,32*K2B,end,DONE,1<<62,0);assert tuple(r0[r]for r in("x0","x1","x2","x7","x8","x16"))==H2S_PORTS;assert h2s_state(r0)==(1<<126,af)and(nd.mem[PW-8]<<64|nd.mem[PW-16])==af;return r0
def deliver_arr(E,ad,labels,nd,msg,end,extra,calls,cert,ac,tag=0,af=1<<66):L,o=msg[1]//32,{};r0=H2S_REGS.get((nd.uid,end))or h2s0(E,ad,labels,nd,end,af);nd.mem[PW]=tag;nd.mem.setdefault(DONE,0);d0,(D0,f0)=nd.mem[DONE],h2s_state(r0);regs=dict(r0,x19=STALE,sp=nd.SP,x30=E.SNT);n=E.a64.a64_run(ad,labels,"h2s",regs,LinkView(nd.mem,msg,[]),out_regs=o,dev=ac and ac.dev,ls64=([WINDOW],[]));assert tuple(o[r]for r in("x0","x1","x2","x7","x8","x16"))==H2S_PORTS;H2S_REGS[(nd.uid,end)]=o;fin=L<K1;ax=D0+(msg[2]if fin else 0)-(1<<126 if fin and f0==1 else 0)-1;j1=ax==0;assert nd.mem[DONE]==d0+j1 and h2s_state(o)==((1<<126,af)if j1 else(ax,f0-fin));calls.append(("h2s",n,cert.h2s_count(L)+extra+(cert.FIN_CNT+(f0==1)if fin else 0)+cert.RCV1+cert.H2S_DONE*j1));return j1
class Sender1:
 def __init__(J,cert):J.cert,J.n,J.slots,J.cnt,J.msgs,J.seq,J.extra=cert,{},{},{},[],[],0
 def block(J,d,n,c):return[8*d,32*n,c,0]+[w for r in J.slots.get(d,[None]*K1)for w in(r or[STALE]*4)]
 def put(J,d,ab):
  s,n=J.slots.setdefault(d,[None]*K1),J.n.get(d,0);s[n]=ab;J.seq.append((d,ab));J.n[d]=(n+1)%K1
  if n+1==K1:J.cnt[d]=J.cnt.get(d,0)+1;J.msgs.append(J.block(d,K1,J.cnt[d]));J.extra+=J.cert.FLUSH1_X
 def finals(J,ds):
  out=[]
  for d in ds:out.append(J.block(d,J.n.get(d,0),J.cnt.get(d,0)+1));J.cnt[d],J.n[d]=0,0
  return out
def state1_ok(mem,aa):
 for d,s in aa.slots.items():
  b=blk1(d)
  if mem[8*d]!=b+32+32*aa.n.get(d,0)or mem[b+16]!=aa.cnt.get(d,0):return False
  if any([mem[b+32*t+32+8*w]for w in range(4)]!=s[t]for t in range(K1)if s[t]is not None):return False
 return True
class Sender:
 def __init__(J,rc0,end,cert,tag=0):J.rc0,J.end,J.cert,J.tag=rc0,end,cert,tag;J.n,J.rc,J.slots,J.msgs,J.seq,J.sent,J.drops,J.extra={},{},{},[],[],{},0,0
 def cursor(J,d):return J.rc.get(d,J.rc0(d))
 def put(J,d,ab):
  s,n=J.slots.setdefault(d,[None]*K2B),J.n.get(d,0);s[n]=ab;J.seq.append((d,ab));J.n[d]=(n+1)%K2B
  if n+1<K2B:return
  c=J.cursor(d);J.rc[d]=c+32*K2B;J.extra+=J.cert.FLUSH2_X if s64(c)<s64(J.end)else J.cert.FLUSH2_D
  if s64(c)<s64(J.end):J.msgs.append([LC2+8*d,c,32*K2B,J.tag]+[w for r in s for w in r]);J.sent[d]=J.sent.get(d,0)+1
  else:J.drops+=1;J.rc[d]=J.end
 def final(J,d):n,c=J.n.get(d,0),J.cursor(d);msg=[LC2+8*d,c,32*n,J.tag]+[w for r in J.slots[d]for w in(r or[STALE]*4)]if n and s64(c)<s64(J.end)else None;k8=(c+32*n)&((1<<64)-1);return msg,(k8 if s64(k8)<s64(J.end)else J.end)
 def fin_count(J,ds):
  f=1
  for d in ds:n=J.n.get(d,0);f+=J.cert.FINH_DEST+(0 if not n else J.cert.FINH_SENT if s64(J.cursor(d))<s64(J.end)else J.cert.FIN_HELD)
  return f
def keys_reference(E,mem,sb,q2,W2,KN,aa,cert):
 meter,out=0,[]
 for l in range(q2):
  e=sb+E.ENTRY*l;w3=(KN+(l<<2))&E.M64
  for s in(0,1):
   f6=E.fracs(mem,e,s);dev=0
   for t in range(3):f=E.record_fields(E.leaf_w(f6,t));meter+=(cert.JPATH if f[0]==1 else cert.JPATH2 if f[0]==2 else 0);dev|=bool(f[0]);aa.put(f[1],[f[3],f[4],W2,w3]);out.append(f)
   meter+=cert.SJOIN if dev else 0
 return meter,out
def plant(E,mem,sb,cd):
 for l,az in enumerate(cd):
  for s in(0,1):
   v=[E.t128(E.BT[m]*c%E.p)for c in E.inverse_dft(az[3*s:3*s+3])for m in(0,1)]
   for o,w in zip(E.FRAC_OFF[s],v):E.put(mem,sb+E.ENTRY*l+o,w,2)
def state_ok(mem,aa):
 for d,s in aa.slots.items():
  b,e=blk2(d),LC2+8*d
  if mem[e]!=b+32+32*aa.n.get(d,0)or mem[e|RC2]!=aa.cursor(d):return False
  if any([mem[b+32*t+32+8*w]for w in range(4)]!=s[t]for t in range(K2B)if s[t]is not None):return False
 return True
M64=(1<<64)-1
def unpack(dm,m,wr=None):
 if len(m)==6:
  for w in range(4):dm[(m[1]+8*w)&M64]=m[2+w]
  return
 for t in range(m[2]//32):
  for w in range(4):dm[(m[1]+32*t+8*w)&M64]=m[4+4*t+w]
  if wr is not None:wr.append((m[1]+32*t)&M64)
def rule_ok(dm,writes,aa,finals):
 ok,cnt=True,{}
 for d,ab in aa.seq:
  t=cnt[d]=cnt.get(d,-1)+1;c=(aa.rc0(d)+32*t)&M64
  if s64(c)<s64(aa.end)and(finals is not None or t<K2B*aa.sent.get(d,0)):ok&=[dm[d].get((c+8*w)&M64)for w in range(4)]==ab
 if finals is not None:ok&=all(f==(((aa.rc0(d)+32*(cnt.get(d,-1)+1))&M64)if s64((aa.rc0(d)+32*(cnt.get(d,-1)+1))&M64)<s64(aa.end)else aa.end)for d,f in finals.items())
 return ok and all(s64(aa.rc0(d))<=s64(a)<s64(aa.end)for d,w in writes.items()for a in w)
TRACE=[]
def reg_trace(cert):
 iv={k:[(lo,hi,b)for r,_,lo,hi,b in cert.REGISTRY if r==k]for k in"GRH"};n=bad=0
 for l0,mem in TRACE:
  for a in mem:n+=1;hit=[b for lo,hi,b in iv[l0]if lo<=a<hi];bad+=len(hit)!=1 or hit[0]=="analytic"
 return dict(nodes=len(TRACE),addresses=n,outside_or_analytic=bad,ok=bad==0 and n>0)
class Node:
 HB,SB,SP=0x100000,(1<<45)+(1<<40)+(1<<34),(34<<40)+0x7000000;R2,E2,E2B=0x3000000,LINK+0x1F000,LINK+0x1F800;CL=LINK+0x200000
 def __init__(J,mem=None,sp=None):cb[0]+=1;J.uid,J.mem,J.SP=cb[0],({}if mem is None else mem),(J.SP if sp is None else sp)
def run_k(E,ad,labels,entry,regs,nd,calls,formula,ac=None,ls64=None):o={};regs=dict(regs,sp=nd.SP,x30=E.SNT);calls.append((entry,E.a64.a64_run(ad,labels,entry,regs,nd.mem,out_regs=o,dev=ac and ac.dev,ls64=ls64),formula));return o
class HomeJoin:
 def __init__(J,cert,CC,ncl):J.cert,J.CC,J.ncl=cert,CC,ncl;J.chains,J.heads,J.nextd,J.recs,J.cum,J.seen={},{},{},{},[],[];J.nls0,J.lists,J.paths,J.segs=0,[0,0],[0,0],{}
 def block(J,m,F,NB):
  cert,c,cnt=J.cert,m[1],m[2];lo=F if s64(F)>s64(c)else c;hi=(c+cnt)&M64;n1,n2=(s64(lo)-s64(c))//32,(s64(hi)-s64(lo))//32;fast=s64(F)<=s64(c)and cnt==32*K2B;J.paths[0 if fast else 1]+=1;f=cert.dlvh_copy(n1,n2);ok=(f==cert.DLVHB if fast else f<=cert.DLVHS);k=(n1%2,n2%2,n1>0,n2>0);J.segs[k]=J.segs.get(k,0)+1
  for t in range(n1,n1+n2):
   a,ab=(c+32*t)&M64,m[4+4*t:8+4*t];k9=not(t%2==0 and t+2<=n1+n2);J.recs[a]=ab;J.seen.append((ab[0],ab[1],a));b=ab[1]&(NB-1);ch=J.chains.setdefault(b,[]);f+=cert.HCHG_NEW
   if ch:
    match=None
    for j,e in enumerate(reversed(ch)):
     f+=cert.HCHG_VISIT
     if e[1]==ab[1]:
      f+=cert.HCHG_VER
      if(e[0]&DM)==(ab[0]&DM):match=(e,j);break
     if j==len(ch)-1:f+=cert.HCHG_END*(j%2);break
    if match is not None:
     e,j=match;f-=1
     if len(J.cum)<J.ncl:J.cum.append((e[2],a));f+=cert.HCHG_LS1-cert.HCHG_LSA*(not k9);J.lists[k9]+=1
     else:J.nls0+=1;f+=cert.HCHG_LS0
   J.nextd[a]=J.heads.get(b,0);J.heads[b]=(asr1(a)+J.CC)&(2**64-1);ch.append((ab[0],ab[1],a))
  return f+cert.RCV,ok
def home_state(nd,ncl,CC):DLVH_REGS.pop(id(nd.mem),None);nd.mem[ST],nd.mem[ST+8],nd.mem[ST+16],nd.mem[ST+24]=nd.CL,nd.CL+16*ncl,CC,-2*CC&(2**64-1)
def part_b(E,cert):
 rg=random.Random(E.SEED+32452843);rt=E.a64.a64_load(open(os.path.join(E.HERE,"kernels","route.s")).read());lk=E.a64.a64_load(open(os.path.join(E.HERE,"kernels","leaf.s")).read())+("keys",);H2=1<<41;ao,ba,G,q2=(0,1,H1,H1+1),(0,1,H2-2,H2-1),4,16;calls=[];af={g:Node(dd(int))for g in range(G)};j8={d1:Node(sp=0x110000)for d1 in ao};homes={(d1,d2):Node(dd(int),sp=LINK+0x110000)for d1 in ao for d2 in ba};pool=[E.planted_x(rg,rg.choice(ao),rg.choice(ba))[1]for _ in range(10)];batches,shared={},0
 for g in af:
  bi=[]
  for _ in range(6*q2):
   if rg.random()<0.15:bi.append(E.planted_x(rg,0,0,rg.choice(pool))[0]);shared+=1
   else:bi.append(E.planted_x(rg,rg.choice(ao),rg.choice(ba))[0])
  batches[g]=[bi[6*l:6*l+6]for l in range(q2)]
 d3,e5,blocks=True,True,0;e4={d1:[]for d1 in ao}
 for g,nd in af.items():
  sb=nd.SB-E.ENTRY*q2
  for t0,n,first in((0,2,LO1+32),(8*H1,2,UP1+32)):run_k(E,rt[0],rt[1],"init3",dict(x0=n,x1=t0,x2=first,x3=B1),nd,calls,cert.INIT3_BKT*n+1)
  plant(E,nd.mem,sb,batches[g]);KN=(g>>64)|(5+g)<<26;aa=Sender1(cert);meter,fl=keys_reference(E,nd.mem,sb,q2,g,KN,aa,cert);e5&=[f[1:]for f in fl]==[E.class_label(x)for az in batches[g]for x in az];ac=Port1(nd.mem);run_k(E,lk[0],lk[1],"keys",dict(x2=g,x3=sb+E.KEYS0,x26=KN),nd,calls,cert.KEY_FIX+cert.KEY_INC*q2+meter+aa.extra,ac,ls64=([(sb,nd.SB)],[]));ok=ac.blocks==aa.msgs and state1_ok(nd.mem,aa);nb=len(ac.blocks)
  for t0,n in((0,2),(8*H1,2)):run_k(E,rt[0],rt[1],"fin1",dict(x7=t0,x1=n),nd,calls,1+cert.FIN1_X*n,ac)
  i6=aa.finals(ao);ok&=ac.blocks[nb:]==i6 and all(nd.mem[8*d]==blk1(d)+32 and nd.mem[blk1(d)+16]==0 for d in ao);d3&=ok;blocks+=len(ac.blocks)
  for m in ac.blocks:e4[dest1(m[0])].append(m)
 a7=Node.R2+3*32*K2B;bh={h:[]for h in homes};b8={h:[]for h in homes};b9=True
 for d1,nd in j8.items():
  for t0,n in((0,2),(H2-2,2)):run_k(E,rt[0],rt[1],"init3",dict(x0=n,x1=LC2+8*t0,x2=blk2(t0)+32,x3=B2),nd,calls,cert.INIT3_BKT*n+1);run_k(E,rt[0],rt[1],"tinit",dict(x0=n,x1=(LC2+8*t0)|RC2,x2=Node.R2),nd,calls,cert.TINIT_ENTRY*n+1)
  aa=Sender(lambda d:Node.R2,a7,cert,Node.E2);ac=Port2(nd.mem);ae=e4[d1];rg.shuffle(ae);dn=[]
  for m in ae:
   e0=aa.extra
   for t in range(m[1]//32):
    ab=list(m[4+4*t:8+4*t]);aa.put(ab[0]&M41,ab);bh[(d1,ab[0]&M41)].append(tuple(ab))
   dn.append(deliver_arr(E,rt[0],rt[1],nd,m,a7,aa.extra-e0,calls,cert,ac,aa.tag,G))
  b9&=(dn==[False]*(len(ae)-1)+[True]and nd.mem[DONE]==1 and h2s_state(H2S_REGS[(nd.uid,a7)])==(1<<126,G)and sum(m[2]for m in ae if m[1]<32*K1)==len(ae));ok=ac.blocks==aa.msgs and not ac.msgs and state_ok(nd.mem,aa);a4=Port2(nd.mem)
  for t0,n in((0,2),(H2-2,2)):run_k(E,rt[0],rt[1],"fin",dict(x0=MB2,x7=LC2+8*t0,x3=(LC2+8*t0)|RC2,x8=blk2(t0)+32,x1=n,x2=Node.E2B,x11=aa.tag,x6=a7),nd,calls,aa.fin_count(range(t0,t0+n)),a4)
  i6=[aa.final(d)for d in ba];ok&=(a4.blocks==[m for m,_ in i6 if m]and a4.msgs==[[LC2+8*d,Node.E2B,c,0,0,0]for d,(_,c)in zip(ba,i6)]);d3&=ok;blocks+=len(ac.blocks)+len(a4.blocks)
  for m in ac.blocks+a4.blocks+a4.msgs:b8[(d1,dest(m[0]))].append(m)
 e6,e7,d5,d6=True,True,True,0;CC=(NEXT_BASE-asr1(Node.R2))%2**64
 for h,nd in homes.items():
  nd.mem[Node.E2]=Node.R2;home_state(nd,64,CC);hj=HomeJoin(cert,CC,64)
  for m in b8[h]:
   if len(m)==6:deliver(E,rt[0],rt[1],nd.mem,m,calls,cert,cnt=CNTE)
   else:f,ap=hj.block(m,nd.mem[Node.E2],2**cert.HB);d5&=ap;deliver(E,rt[0],rt[1],nd.mem,m,calls,cert,fused=True,r2=(Node.R2,a7),formula=f)
  F,F0=nd.mem[Node.E2B],Node.R2;start=F0-(F0-Node.R2)%(32*K2B);g9=1+-(-(F-start)//(32*K2B));b9&=nd.mem[CNTE]==1 and nd.mem[CNT]+1==len(b8[h])==g9;bj=[tuple(nd.mem[x+8*k]for k in range(4))for x in range(Node.R2,F,32)];e6&=bj==bh[h][:len(bj)]and F==min(Node.R2+32*len(bh[h]),a7);e7&=all(ab[0]&M41==h[1]and ab[2]in af and(ab[3]>>26)==5+ab[2]and 0<=((ab[3]>>2)&((1<<24)-1))<q2 for ab in bj);d5&=([(nd.mem[nd.CL+16*i],nd.mem[nd.CL+16*i+8])for i in range(len(hj.cum))]==hj.cum and all(nd.mem[HT+8*b]==v for b,v in hj.heads.items()));d6+=len(hj.cum)
 groups={}
 for g in af:
  for az in batches[g]:
   for x in az:groups.setdefault(E.exk(x)[0],set()).add(g)
 cnt,j9={},{}
 for e,c,f in calls:cnt[e],j9[e]=cnt.get(e,0)+c,j9.get(e,0)+f
 ag=dict(leaves=6*q2*G,home_candidates=d6,counts=cnt,formulas_ok=all(c==f for _,c,f in calls)and e5 and d3 and e6 and e7 and d5 and b9 and d6>0);ag["block_paths"]=block_paths(E,rt,lk,rg,cert);ag["grouped_paths"]=grouped_paths(E,rt[0],rt[1],rg,cert);ag["ls64"]=ls64_checks(E,rg);ag["va56"]=va56_checks(E);ag["encoding"]=encoding_checks(E);ag["phase4_round"]=phase4_round(E,rt,rg,cert);ag["y_ring"]=yring_checks(E,rg);ag["gate"]=gate_checks(E,rt,cert);ag["wait_dispatch"]=wait_dispatch_checks(E,rt,rg,cert);ag["phase0"]=phase0_xemit_hhclr(E,rt,rg,cert);TRACE.extend([("G",n.mem)for n in af.values()]+[("R",n.mem)for n in j8.values()]+[("H",n.mem)for n in homes.values()]);ag["registry_trace"]=reg_trace(cert);return ag
def yring_checks(E,rg):
 ad,bc=E.a64.a64_load(open(os.path.join(E.HERE,"kernels","leaf.s")).read());aq,ar=(1<<45)+(257<<31),(1<<45)+(1<<40)+(1<<34)-512*32;e8=([(3*(1<<58),63*(1<<58)),(1<<58,3*(1<<58)),(0xfc00000000000000,1<<64),(aq,aq+4096),(ar,ar+512*32)],[(ar,ar+512*32)])
 class ReadOnlyYR(dict):
  def __missing__(J,a):return 0
  def __setitem__(J,a,v):assert not(aq<=a<aq+4096),"write into shared read-only YR";super().__setitem__(a,v)
 mem=ReadOnlyYR();ok=True;Xs=[(1,2)]*8
 for l in range(32):E.x_words(mem,ar+512*l,Xs,l)
 E.a64.a64_run(ad,bc,"fwd",dict(x0=32,x1=aq+64,x3=ar,x24=1,x25=0,x26=0,x27=0,sp=(34<<40)+0x7000000,x30=E.SNT),mem,ls64=e8);E.a64.a64_run(ad,bc,"bwd",dict(x3=ar+512*32,x1=(34<<40)+0x10008,sp=aq+512*3,x0=1,x2=0,x14=0,x15=0,x30=E.SNT),mem,ls64=e8);h0,h4=E.a64.a64_load("w:\n\tstr xzr,[sp]\n\tret\n")
 try:E.a64.a64_run(h0,h4,"w",dict(sp=aq,x30=E.SNT),mem);d7=False
 except AssertionError:d7=True
 return dict(ok=ok and d7,yr_readonly_enforced=d7)
def phase4_round(E,rt,rg,c):
 C,ok,N,BR,BH,bs=[],True,Node,(1<<52)+0x10000,LINK+0x20000,LINK+0x18000;d1,d2,gs=rg.randrange(1<<41),rg.randrange(1<<41),[rg.randrange(3<<64,1<<66),rg.randrange(1<<64)];h,rl=N(dd(int),sp=LINK+0x110000),N(dd(int),sp=0x110000)
 def run(nd,e,regs,f,base=None,msg=None,lb=LINK_MSG):p=base and Port2(nd.mem,base);C.append((E.a64.a64_run(rt[0],rt[1],e,dict(regs,sp=nd.SP,x19=STALE,x30=E.SNT,x28=1,x29=CNT),nd.mem if msg is None else LinkView(nd.mem,msg,[],lb),dev=p and p.dev),f));return p and p.msgs
 def q(nd,m,cap=c.RELAY_CAP,lb=LINK_MSG+1024,bb=BR,ct=CNT):lo0,hi0,oc=nd.mem[ct],nd.mem[ct+2048],nd.mem[bb-8];C.append((E.a64.a64_run(rt[0],rt[1],"dlvq",dict(x0=lb,x9=64,x10=bb-8,x11=bb+64*cap,sp=nd.SP,x19=STALE,x30=E.SNT,x28=1,x29=ct),LinkView(nd.mem,m,[],lb)),(c.DLVQ if oc<bb+64*cap else c.DLVQ_DROP)+(c.DLVQ_CARRY if lo0==(1<<64)-1 else 0)+c.RCV));return(nd.mem[ct+2048]<<64|nd.mem[ct])==((hi0<<64|lo0)+1)
 R=[[rg.randrange(1<<64),0,g%2**64,g>>64|rg.randrange(1<<62)<<2]for g in gs]
 for i in range(8):h.mem[N.R2+8*i]=R[i//4][i%4]
 h.mem[N.CL],h.mem[N.CL+8],rl.mem[BR-8],hd=N.R2,N.R2+32,BR,2**45+8*d2;m1=run(h,"reqh",dict(x0=bs,x12=N.CL,x9=1,x10=hd,x13=BH),c.REQH_C+1,bs);ok&=m1==[[0,0,r[2],r[3],hd,BH+256*i]for i,r in enumerate(R)]and all(q(rl,m)for m in m1);m2=run(rl,"fwdq",dict(x0=MB2+64,x1=BR,x6=BR+128,x7=8*d1),2*c.FWDQ+3,MB2+64);bk=[[rg.randrange(1<<64)for _ in range(32)]for _ in gs];BG=96*(1<<40);br=LINK_GEN+0x18000;l4=LINK_GEN+0x10000
 for i,m in enumerate(m2):
  g=N(dd(int));g.mem[BG-8]=BG;TRACE.append(("G",g.mem));ok&=(m[1]&3)<<64|m[0]==gs[i]and q(g,m,cap=c.REQ_CAP,lb=LINK_GEN_MSG+1024,bb=BG,ct=l4)and[g.mem[BG+8*k]for k in range(4)]==[gs[i]>>64,R[i][3],BR+64*i<<8,8*d1]
  for k in range(32):g.mem[BG+4096+8*k]=bk[i][k]
  m3=run(g,"sendr",dict(x0=br,x1=BG+4096,x9=8,x7=8*d1,x8=g.mem[BG+16]),8*c.SENDR+1,br);ok&=len(m3)==8 and all(q(rl,m)for m in m3)
 m4=run(rl,"fwdr",dict(x0=MB2,x1=BR+128,x6=BR+64*18),16*c.FWDR+3,MB2)
 for m in m4:run(h,"dlv",dict(x0=LINK_MSG),c.DLV+c.RCV,msg=m)
 rl.mem[BR-8]=BR;d=N(dd(int),sp=0x110000);d.mem[BR-8]=BR+64*c.RELAY_CAP;w=N(dd(int),sp=0x110000);w.mem[BR-8]=BR+64*(c.RELAY_CAP-1);w.mem[CNT]=(1<<64)-2;cw=N(dd(int),sp=0x110000);cw.mem[BR-8]=BR+64*c.RELAY_CAP;cw.mem[CNT]=(1<<64)-1;evc=[];E.a64.a64_run(rt[0],rt[1],"dlvq",dict(x0=LINK_MSG+1024,x9=64,x10=BR-8,x11=BR+64*c.RELAY_CAP,sp=cw.SP,x19=STALE,x30=E.SNT,x28=1,x29=CNT),LinkView(cw.mem,[5]*6,[],LINK_MSG+1024),dev=(CNT,4096,lambda o,v:evc.append((cw.mem[CNT+2048],cw.mem[CNT])),lambda:None,lambda:True));ok&=([h.mem[BH+8*k]for k in range(64)]==bk[0]+bk[1]and len(m4)==16 and q(d,[1]*6)and q(d,[1]*6)and d.mem[BR-8]==BR+64*c.RELAY_CAP and d.mem[BR]==0 and q(w,[2]*6)and q(w,[3]*6)and q(w,[4]*6)and w.mem[BR-8]==BR+64*c.RELAY_CAP and w.mem[CNT]==1 and w.mem[CNT+2048]==1 and 0 not in w.mem and evc==[(0,0),(1,0)]and all(x==y for x,y in C)and max(c.PH4_HOPS.values())<=c.MSG);TRACE.extend([("H",h.mem),("R",rl.mem),("R",d.mem),("R",w.mem),("R",cw.mem)]);return dict(ok=ok,calls=len(C),hops=c.PH4_HOPS)
def encoding_checks(E):
 a64=E.a64;good=("ldp x1,x2,[x3,#-512]","ldp x1,x2,[x3,#504]","stp x1,x2,[x3],#504","stp x1,x2,[x3,#-512]!","ldr x1,[x3,#32760]","ldr x1,[x3,#-256]","ldr x1,[x3,#255]","str x1,[x3],#-256","add x1,x1,#4095","add x1,x1,#0x1000","tst x3,#0x3fffe0000","mov sp,x1","add sp,sp,#512","cset w9,hs","csel x1,x2,x3,lo","mov x1,#0xfff8000000000000","movk x4,#1,lsl #32","stlur x1,[x3,#-256]","ldr x22,[x3,#-152]","and x1,x2,#0x5555555555555555","ldaxr x2,[x0]","blr x30","wfe");bad=("stp x25,x26,[x3],#768","ldp x21,x22,[x3,#-768]","ldp x1,x2,[x3,#508]","ldp x1,x2,[x3,#-520]","ldp x1,x2,[x3,#12]","ldr x22,[x3,#-664]","ldr x1,[x3,#32768]","ldr x1,[x3,#260]","str x1,[x3],#256","ldr x1,[x3,#-257]!","add x1,x1,#4097","tst x3,#0x1c001","and x1,x2,#0","mov x1,#0x12345","movk x4,#1,lsl #8","stlur x1,[x3,#256]","ldp x1,x1,[x3]","ldp x3,x2,[x3],#16","ldadd x1,x2,[x3,#8]","ldaxr sp,[x0]","ldaxr x1,[x0,#8]","blr sp","wfe x0","ubfx x1,x2,#60,#8","bfi x1,x2,#40,#25","ldr x1,[x3,#-264]","add x1,xzr,#5","tst sp,#0x1c000","ldp xzr,xzr,[x3]","ldr sp,[x3]","mov xzr,#0x3ffe0","movk sp,#1,lsl #48","tbz x3,#64,1f","frob x1,x2","mul x1,x2,#1","adc x1,x2,#3","csel x1,x2,#3,eq","madd x1,x2,x3,#4","cset x1,x2","mvn x1,#2","mul w1,x2,x3","dmb sy","umulh w1,w2,w3","mvn x1,x2,lsl #64","movk x1,#1,lsr #16","stadd x1,[x2,#0]","ld64b x0,[x2,#0]","ldp x1,x2,[x3]!","ld64b x21,[x1]","ld64b x24,[x1]","add x1,x2,#1_0")
 def verdict(ins):
  try:a64.a64_load("f:\n\t"+ins+"\n\tret\n");return True
  except AssertionError as ex:assert"not encodable"in str(ex),ex;return False
 out=dict(good_accepted=sum(map(verdict,good)),bad_rejected=sum(not verdict(i)for i in bad),good=len(good),bad=len(bad));out["ok"]=out["good_accepted"]==len(good)and out["bad_rejected"]==len(bad);return out
def va56_checks(E):
 a64=E.a64;ad,bc=a64.a64_load("ld:\n\tldr x2,[x1]\n\tret\n");out={}
 for name,a in(("zero",0),("mid_2^56",1<<56),("mid_2^60",1<<60),("top",2**64-8),("wrap_top",2**64-4)):
  try:a64.a64_run(ad,bc,"ld",dict(x1=a,x30=E.SNT),dd(int));out[name]="ran"
  except AssertionError as ex:out[name]="rejected"if"2^64"in str(ex)else"other: "+str(ex)
 out["ok"]=[out[k]for k in("zero","mid_2^56","mid_2^60","top")]==["ran"]*4 and out["wrap_top"]=="rejected";return out
def ls64_checks(E,rg):
 a64=E.a64;aj=open(os.path.join(E.HERE,"kernels","route.s")).read();a9,b0=a64.a64_load("ld:\n\tld64b x20,[x1]\n\tret\nst:\n\tst64b x16,[x2]\n\tret\novl:\n\tld64b x12,[x15]\n\tret\n");R2=(Node.R2,Node.R2+4*32*K2B);msg=[rg.getrandbits(64)for _ in range(4+4*K2B)]
 def run(ad,labels,entry,regs,mem,lsr=([WINDOW,R2],[R2])):
  try:a64.a64_run(ad,labels,entry,dict(regs,x30=E.SNT),mem,ls64=lsr);return"ran"
  except AssertionError as ex:return"LS64 rejected"if"LS64"in str(ex)else"other: "+str(ex)
 out,o={},{};a64.a64_run(a9,b0,"ld",dict(x1=LINK+64,x30=E.SNT),LinkView({},msg,[]),out_regs=o,ls64=([WINDOW],[]));out["ld64b_window_slot_0_1"]=[o["x%d"%(20+k)]for k in range(8)]==msg[4:12];mem={};out["st64b_r2_aligned"]=run(a9,b0,"st",dict(x2=R2[0]),mem);out["st64b_r2_written"]=all(mem.get(R2[0]+8*k)==0 for k in range(8));out["ld64b_odd_slot"]=run(a9,b0,"ld",dict(x1=LINK+32),LinkView({},msg,[]));out["ld64b_outside_regions"]=run(a9,b0,"ld",dict(x1=0x5000040),{});out["st64b_into_window"]=run(a9,b0,"st",dict(x2=LINK+64),LinkView({},msg,[]));out["st64b_r2_misaligned"]=run(a9,b0,"st",dict(x2=R2[0]+32),{})
 def load(text):
  try:a64.a64_load(text);return"loaded"
  except AssertionError as ex:return"LS64 rejected"if"LS64"in str(ex)else"other: "+str(ex)
 out["ld64b_x24"]=load("hi:\n\tld64b x24,[x1]\n\tret\n");out["ld64b_x21"]=load("odd:\n\tld64b x21,[x1]\n\tret\n");out["ld64b_base_overlap"]=run(a9,b0,"ovl",dict(x15=LINK+64),LinkView({},msg,[]))
 def fused(text,c,r2):
  ad,labels=a64.a64_load(text);mem=dd(int);CC=(NEXT_BASE-asr1(r2[0]))%2**64;mem[ST],mem[ST+8],mem[ST+16],mem[ST+24]=Node.CL,Node.CL+1024,CC,-2*CC%2**64;i9=LINK+0x300000;mem[i9]=c;m=[0,c,32*K2B,i9]+msg[4:];dr=dlvh0(E,ad,labels,mem,CNT);r=run(ad,labels,"dlvh",dict(dr,x29=CNT),LinkView(mem,m,[]),([WINDOW,r2],[r2]));return r if r!="ran"else("ran, R2 exact"if all(mem[c+8*w]==m[4+w]for w in range(4*K2B))else"ran, R2 wrong")
 out["dlvh_fresh_block"]=fused(aj,R2[0],R2);e9=(LINK-32*K2B,LINK);out["dlvh_r2_top_below_link"]=fused(aj,e9[0],e9);out["dlvh_r2_into_link_rejected"]=fused(aj,LINK,(LINK,LINK+32*K2B));f3=(R2[0]+32,R2[1]+32);out["dlvh_r2_at_32_mod_64"]=fused(aj,f3[0],f3);bt=aj.replace("6:\nadd x13,x14,x2,asr #1\nDLVHPR x0\n","6:\nadd x13,x14,x2,asr #1\nsub x12,x0,#32\nDLVHPR x12\n");out["dlvh_mutant_ld64b_at_odd_slot"]=fused(bt,R2[0],R2)if bt!=aj else"mutation not applied";hm=aj.replace("mov x1,x0\n","sub x1,x0,#32\n");hp,hl=a64.a64_load(hm);regs=dict(x0=H2S_PORTS[0],x1=H2S_PORTS[1],x2=H2S_PORTS[2],x7=H2S_PORTS[3],x8=H2S_PORTS[4],x16=H2S_PORTS[5],x10=LC2,x18=32*K2B,x12=1<<40,x3=PW,x29=DONE,x6=1<<62,x14=1<<62,x5=16);out["h2s_mutant_ld64b_at_odd_slot"]=run(hp,hl,"h2s",regs,LinkView(dd(int),[0,32*2,1,0]+msg[4:4+4*K1],[]),([WINDOW],[]))if hm!=aj else"mutation not applied";bu=("ld64b_odd_slot","ld64b_outside_regions","st64b_into_window","st64b_r2_misaligned","ld64b_x24","ld64b_x21","ld64b_base_overlap","dlvh_r2_at_32_mod_64","dlvh_mutant_ld64b_at_odd_slot","h2s_mutant_ld64b_at_odd_slot");h5=(out["ld64b_window_slot_0_1"]and out["st64b_r2_aligned"]=="ran"and out["st64b_r2_written"]and out["dlvh_fresh_block"]=="ran, R2 exact"and out["dlvh_r2_top_below_link"]=="ran, R2 exact"and"write into link aperture"in out["dlvh_r2_into_link_rejected"]and all(out[k]=="LS64 rejected"for k in bu));return dict(ok=h5,rejected=len(bu))
def grouped_paths(E,ad,labels,rg,cert):
 GG,k0,a5=3,6,2;NB=1<<k0;b2=32*K2B*128;h6=((600,9),(1100,530),(500,7))
 def run_case(bg,ncl=64):
  calls,nd=[],Node(dd(int));av,ah=LINK+0x400000,LINK+0x500000;i4=((1<<64)-32*K2B*2)if bg else nd.R2;E2P=(nd.E2,nd.E2+0x800000);pf=lambda g,s:av+8*s if g==0 else E2P[(g-1)%2]+32*s;au={}
  for w in range(cert.ZW):nd.mem[HT+8*w]=0x5A5A5A5A5A5A5A5A
  run_k(E,ad,labels,"zero16384",dict(x0=cert.ZW,x1=HT),nd,calls,cert.ZERO_PASS+1);ok=all(nd.mem[HT+8*w]==0 for w in range(cert.ZW));run_k(E,ad,labels,"init3",dict(x0=a5,x1=av,x2=i4,x3=b2),nd,calls,cert.INIT3_BKT*a5+1);run_k(E,ad,labels,"dlvh0",dict(x29=CNT),nd,calls,cert.DLVH0);ok&=all(nd.mem[av+8*s]==((i4+b2*s)&M64)for s in range(a5));CC=(NEXT_BASE-asr1(i4))%2**64;home_state(nd,ncl,CC);hj=HomeJoin(cert,CC,ncl);ran,stopped=[],False;rc=[(i4+b2*s)&M64 for s in range(a5)];F=list(rc);aw,nb=[[None]*K2B for _ in range(a5)],[0]*a5
  for g in range(GG):
   if stopped:break
   ran.append(g);au[g]=[(rg.getrandbits(64),rg.getrandbits(64))for _ in range(3)]
   for k,(px,py)in enumerate(au[g]):nd.mem[ah+(g*3+k)*16]=px;nd.mem[ah+(g*3+k)*16+8]=py
   f7=[]
   for s in range(a5):
    j4=h6[g][s]+(0 if bg or s else 1);recs=[]
    for _ in range(j4):recs.append([rg.getrandbits(64),rg.getrandbits(7)<<57|rg.getrandbits(k0),rg.getrandbits(64),rg.getrandbits(64)])
    if hj.seen:
     bp,bo=hj.seen[rg.randrange(len(hj.seen))][:2];recs[0][1]=bo;recs[0][0]=(bp&DM)|rg.getrandbits(41)
    recs[1][1]=recs[0][1];recs[1][0]=(recs[0][0]&DM)|rg.getrandbits(41);msgs=[]
    for ab in recs:
     aw[s][nb[s]]=ab;nb[s]+=1
     if nb[s]==K2B:msgs.append([0,rc[s],32*K2B,pf(g,s)]+[w for r in aw[s]for w in r]);rc[s]=(rc[s]+32*K2B)&M64;nb[s]=0
    if nb[s]:msgs.append([0,rc[s],32*nb[s],pf(g,s)]+[w for r in aw[s]for w in r or[STALE]*4])
    f7.append(msgs)
   order=[m for st in f7 for m in st];rg.shuffle(order)
   for m in order:
    s=[t for t in range(a5)if pf(g,t)==m[3]][0];f,ap=hj.block(m,F[s],NB);ok&=ap;f9=[((1<<64)-32*K2B*2,1<<64),(0,1<<51)]if bg else[(nd.R2,nd.R2+a5*b2)];n=E.a64.a64_run(ad,labels,"dlvh",dict(DLVH_REGS.get(id(nd.mem))or dlvh0(E,ad,labels,nd.mem,CNT),x29=CNT,x30=E.SNT),LinkView(nd.mem,m,[]),out_regs=(bn:={}),ls64=([WINDOW]+f9,f9));DLVH_REGS[id(nd.mem)]=bn;calls.append(("dlvh",n,f));ok&=all(nd.mem[(m[1]+8*w)&M64]==m[4+w]for w in range(4*K2B))
   ok&=all(nd.mem[HT+8*b]==hj.heads.get(b,0)for b in range(NB));ok&=all(nd.mem[(asr1(a)+CC)&M64]==hj.nextd[a]and nd.mem[((asr1(a)+CC)&M64)+8]==hj.recs[a][1]for a in hj.nextd);ok&=[(nd.mem[nd.CL+16*i],nd.mem[nd.CL+16*i+8])for i in range(len(hj.cum))]==hj.cum;ok&=nd.mem[ST]==nd.CL+16*len(hj.cum)and len(hj.cum)<=ncl
   for s in range(a5):F[s]=(rc[s]+32*nb[s])&M64;deliver(E,ad,labels,nd.mem,[0,E2P[g%2]+32*s,F[s],0,0,0],calls,cert)
   ok&=all(nd.mem[E2P[g%2]+32*s]==F[s]for s in range(a5));d8=[((i4+b2*s+32*k)&M64)for s in range(a5)for k in range((s64(F[s])-s64((i4+b2*s)&M64))//32)];ok&=sorted(hj.nextd)==sorted(d8);f8=[]
   for b in range(NB):
    u=hj.heads.get(b,0)
    while u:a=(2*((u-CC)%2**64))&M64;f8.append(a);u=hj.nextd[a]
   ok&=sorted(f8)==sorted(d8)and all([nd.mem[(a+8*w)&M64]for w in range(4)]==hj.recs[a]for a in d8)
   if bg and g==1:stopped=any(nd.mem[a+8]==nd.mem[b+8]and nd.mem[a]&DM==nd.mem[b]&DM for a,b in hj.cum)
  ok&=all(nd.mem[a+8]==nd.mem[b+8]and nd.mem[a]&DM==nd.mem[b]&DM for a,b in hj.cum);parent={}
  def find(u):
   while parent.get(u,u)!=u:u=parent[u]
   return u
  for a,b in hj.cum:parent[find(a)]=find(b)
  f4={}
  for k2,k3,a in hj.seen:f4.setdefault((k3,k2&DM),[]).append(a)
  true_pairs=[(r[i],r[i+1])for r in f4.values()for i in range(len(r)-1)]
  if ncl>len(hj.cum):ok&=all(find(a)==find(b)for a,b in true_pairs)and len(true_pairs)>0
  ok&=all(c==f_ for _,c,f_ in calls);f5=[(nd.mem[ah+16*k],nd.mem[ah+16*k+8])for k in range(3)]==au[0];j0=sum(1 for _,k3,_ in hj.seen if k3&0x30);ok&=f5 and j0>0 and len({k3&(NB-1)for _,k3,_ in hj.seen})>16;DLVH_REGS.pop(id(nd.mem),None);return dict(ok=ok,groups_run=ran,candidates=len(hj.cum),stopped=stopped,stop_pair_listed=any(nd.mem[a+8]==nd.mem[b+8]for a,b in hj.cum),true_pairs=len(true_pairs),list_skips=hj.nls0,old_x_served=f5,records_with_bucket_bits_4_5=j0,fast_slow_blocks=hj.paths,listings_direct_late=hj.lists)
 a,b,c=run_case(True),run_case(False),run_case(False,ncl=1);return dict(ok=a["ok"]and b["ok"]and c["ok"]and a["groups_run"]==[0,1]and b["groups_run"]==[0,1,2]and a["stopped"]and a["stop_pair_listed"]and a["candidates"]>0 and b["candidates"]>0 and c["candidates"]==1 and c["list_skips"]>0,old_x_served=a["old_x_served"]and b["old_x_served"]and c["old_x_served"],stop_cands=a["candidates"],full_cands=b["candidates"],skips=c["list_skips"])
def block_paths(E,rt,lk,rg,cert):
 calls,out=[],{};nd=Node(dd(int));ax=(0,H1,1,2,H1+1);xs=[E.planted_x(rg,d,rg.randrange(4),j=t%3)[0]for t,d in enumerate(ax)];ds=[t for t,n in enumerate((804,522,256,2,0))for _ in range(n)];rg.shuffle(ds);q2=len(ds)//6;leaves=[[xs[t]for t in ds[6*l:6*l+6]]for l in range(q2)];sb=nd.SB-E.ENTRY*q2;plant(E,nd.mem,sb,leaves)
 for t0,n,first in((0,3,LO1+32),(8*H1,2,UP1+32)):run_k(E,rt[0],rt[1],"init3",dict(x0=n,x1=t0,x2=first,x3=B1),nd,calls,cert.INIT3_BKT*n+1)
 aa=Sender1(cert);meter,fl=keys_reference(E,nd.mem,sb,q2,7,7<<26,aa,cert);ac=Port1(nd.mem);run_k(E,lk[0],lk[1],"keys",dict(x2=7,x3=sb+E.KEYS0,x26=7<<26),nd,calls,cert.KEY_FIX+cert.KEY_INC*q2+meter+aa.extra,ac,ls64=([(sb,nd.SB)],[]));ok=ac.blocks==aa.msgs and state1_ok(nd.mem,aa);nb=len(ac.blocks)
 for t0,n in((0,3),(8*H1,2)):run_k(E,rt[0],rt[1],"fin1",dict(x7=t0,x1=n),nd,calls,1+cert.FIN1_X*n,ac)
 i6=aa.finals(sorted(ax));ok&=ac.blocks[nb:]==i6;ok&=[f[:2]for f in fl]==[(t%3,ax[t])for l in range(q2)for t in ds[6*l:6*l+6]];out["keys_fin1"]=dict(records=len(ds),blocks_sent=nb,finals=len(i6),empty_finals=sum(m[1]==0 for m in i6),upper_half_blocks=sum(m[0]>=8*H1 for m in ac.blocks),model_ok=ok);nd=Node();ai=[rg.getrandbits(41)for _ in range(4)];g2=[K1]*6+[0,1,2,3,7,15,31,63,127,128,254,255]+[K1]*8;rg.shuffle(g2);dl=[ai[t]for t,n in enumerate((3000,1000,400,100))for _ in range(n)];rg.shuffle(dl);end,m=nd.R2+2*32*K2B,{ai[0]:1,ai[1]:2,ai[2]:1,ai[3]:1};aa=Sender(lambda d:end-32*K2B*m[d],end,cert)
 for d in ai:nd.mem[LC2+8*d],nd.mem[(LC2+8*d)|RC2]=blk2(d)+32,end-32*K2B*m[d]
 k,ae=0,[]
 for n in g2:recs=[[rg.getrandbits(23)<<41|dl[(k+t)%len(dl)],rg.getrandbits(64),rg.getrandbits(64),rg.getrandbits(64)]for t in range(n)];k+=n;ae.append(recs)
 ac=Port2(nd.mem);j5=0
 for a8 in ae:
  e0=aa.extra
  for ab in a8:aa.put(ab[0]&M41,ab)
  j5+=len(a8)<K1;deliver_arr(E,rt[0],rt[1],nd,[0,32*len(a8),1,0]+[w for r in a8 for w in r]+[STALE]*4*(K1-len(a8)),end,aa.extra-e0,calls,cert,ac)
 bd=ac.blocks==aa.msgs and not ac.msgs and state_ok(nd.mem,aa);bd&=h2s_state(H2S_REGS[(nd.uid,end)])==((1<<126)+j5-len(ae),(1<<66)-j5)and nd.mem[DONE]==0;dm,wr={d:{}for d in ai},{d:[]for d in ai}
 for msg in ac.blocks:unpack(dm[dest(msg[0])],msg,wr[dest(msg[0])])
 bd&=rule_ok(dm,wr,aa,None);out["h2s"]=dict(records=k,blocks_sent=len(ac.blocks),blocks_dropped=aa.drops,end_stored_back=aa.cursor(ai[0])==end,model_and_rule_ok=bd);out["two_group_finals"]=tg=two_group_finals(E,rt,rg,cert,calls);out["wide_counts"]=wide_counts(E,rt,cert,calls);b4=Node(dd(int));j2=Port1(b4.mem);run_k(E,rt[0],rt[1],"init3",dict(x0=1,x1=0,x2=LO1+32,x3=B1),b4,calls,cert.INIT3_BKT+1)
 if(cert.DELTA_F-cert.FIX_BATCH+1)<=cert.DELTA_F-cert.FIX_BATCH:run_k(E,lk[0],lk[1],"keys",dict(x2=7,x3=sb+E.KEYS0,x26=0),b4,calls,0,j2)
 run_k(E,rt[0],rt[1],"fin1",dict(x7=0,x1=1),b4,calls,1+cert.FIN1_X,j2);out["meter_skip_ok"]=(j2.blocks==Sender1(cert).finals([0])and deliver_arr(E,rt[0],rt[1],Node(),j2.blocks[0],Node.R2+32*K2B,0,calls,cert,None,0,1));out["counts_equal_formulas"]=all(c==f for _,c,f in calls);h7=(out["counts_equal_formulas"]and ok and bd and out["meter_skip_ok"]and out["keys_fin1"]["empty_finals"]==2 and out["keys_fin1"]["upper_half_blocks"]>2 and out["h2s"]["blocks_sent"]>=1 and out["h2s"]["blocks_dropped"]>=1 and out["h2s"]["end_stored_back"]and tg["ok"]and out["wide_counts"]);return dict(ok=h7,meter_skip_ok=out["meter_skip_ok"],dropped=out["h2s"]["blocks_dropped"])
def wide_counts(E,rt,cert,calls):
 nd,end=Node(),Node.R2+32*K2B;key=(nd.uid,end);run=lambda c:deliver_arr(E,rt[0],rt[1],nd,[0,0,c,0]+[STALE]*4*K1,end,0,calls,cert,None,0,1<<66)
 def put(D,f):H2S_REGS[key]=dict(H2S_REGS[key],x13=D%2**64,x14=D>>64,x15=f%2**64,x5=f>>64)
 ok=not run(1)and h2s_state(H2S_REGS[key])==(1<<126,(1<<66)-1);put((1<<126)+(1<<64),1);ok&=not run(1)and h2s_state(H2S_REGS[key])==(1<<64,0);put((1<<126)-1,1);return ok and run(2)and nd.mem[DONE]==1 and h2s_state(H2S_REGS[key])==(1<<126,1<<66)
def two_group_finals(E,rt,rg,cert,calls):
 nd=Node();d,base=(1<<41)-1,nd.R2;end=base+60*32*K2B;nd.mem[LC2+8*d],nd.mem[(LC2+8*d)|RC2]=blk2(d)+32,base;aa=Sender(lambda _d:base,end,cert);dm,wr,finals,recs={},[],[],[]
 for n in(700,700):
  group=[[rg.getrandbits(23)<<41|d,rg.getrandbits(64),rg.getrandbits(64),rg.getrandbits(64)]for _ in range(n)];h8=len(aa.msgs);recs.append(group);ac=Port2(nd.mem)
  for c in range(0,n,K1):
   e0=aa.extra;a8=group[c:c+K1]
   for ab in a8:aa.put(d,ab)
   deliver_arr(E,rt[0],rt[1],nd,[0,32*len(a8),1,0]+[w for r in a8 for w in r]+[STALE]*4*(K1-len(a8)),end,aa.extra-e0,calls,cert,ac)
  bf=ac.blocks==aa.msgs[h8:]and not ac.msgs and state_ok(nd.mem,aa);a4=Port2(nd.mem);run_k(E,rt[0],rt[1],"fin",dict(x0=MB2,x7=LC2+8*d,x3=(LC2+8*d)|RC2,x8=blk2(d)+32,x1=1,x2=nd.E2,x6=end),nd,calls,aa.fin_count([d]),a4);k1,k4=aa.final(d);bf&=a4.blocks==([k1]if k1 else[])and a4.msgs==[[LC2+8*d,nd.E2,k4,0,0,0]]
  for m in ac.blocks+a4.blocks+a4.msgs:unpack(dm,m,wr if len(m)>6 else None)
  finals.append((k4,bf,ac.blocks))
 f1,f2,d9=finals[0][0],finals[1][0],finals[1][2];q=700//K2B*K2B;g3=len(d9)>=1 and d9[0][1]==base+32*q and d9[0][4:4+4*(700-q)]==[w for r in recs[0][q:700]for w in r];ok=(all(bf for _,bf,_ in finals)and rule_ok({d:dm},{d:wr},aa,{d:f2})and g3 and(f1,f2,dm[nd.E2])==(base+32*700,base+32*1400,f2)and aa.cursor(d)==nd.mem[(LC2+8*d)|RC2]==base+1400//K2B*32*K2B and all([dm.get(base+32*t+8*w)for w in range(4)]==recs[0][t]for t in range(700))and all([dm.get(f1+32*t+8*w)for w in range(4)]==recs[1][t]for t in range(700)));return dict(ok=ok,first_final_rewritten_identically=g3)
def gate_checks(E,rt,cert):
 rprog,rlab=rt[0],rt[1];K=cert.K_GATE;klo,bq=K&((1<<64)-1),K>>64;cases=[(0,True),(K-1,True),(K,False),(K+1,False),(((bq-1)<<64)|((1<<64)-1),True),((bq<<64)|((1<<64)-1),False),((1<<128)-1,False),];ok=True;j6=set()
 for u,l5 in cases:o={};n=E.a64.a64_run(rprog,rlab,"gate",dict(x0=u&((1<<64)-1),x1=u>>64,x2=klo,x3=bq,x30=E.SNT),{},out_regs=o);j6.add(n);ok&=(o["x0"]==(1 if l5 else 0))
 ok&=(j6=={4}and 2*cert.COIN_SCALAR+4+8<=cert.GATE_INSNS);return dict(ok=ok,counts=sorted(j6),cases=len(cases))
def phase0_xemit_hhclr(E,rt,rg,cert):
 rprog,rlab=rt[0],rt[1];g0,g1,a6=Node(dd(int)),Node(dd(int)),Node(dd(int));CNM_BASE=90*(1<<40);YR_B=E.YR0
 for nd in(g0,g1,a6):nd.mem[CNM_BASE]=rg.getrandbits(64);nd.mem[CNM_BASE+56]=rg.getrandbits(64)
 e1=rg.getrandbits(193)<<63;d1=(e1>>102)&M41;d2=(e1>>143)&M41;A0=e1+0x123456789ABCDEF;A1=A0+(1<<61);h9=A0^(1<<102);rl=Node(dd(int),sp=0x110000);hm=Node(dd(int),sp=LINK+0x110000);calls=[]
 for nd in(g0,g1,a6):first=(LO1+32)if d1<H1 else(UP1+32);off=d1 if d1<H1 else d1-H1;run_k(E,rprog,rlab,"init3",dict(x0=1,x1=8*d1,x2=first+off*B1,x3=B1),nd,calls,cert.INIT3_BKT+1)
 b5=M41<<3;e2=777;i0=778
 for u in range(K1):Av=A0 if u==K1-1 else((A0^((u+1)<<185))&~(M41<<102))|(d1<<102);E.put(g0.mem,YR_B+u*512+288,Av)
 E.put(g1.mem,YR_B+288,A1);E.put(a6.mem,YR_B+288,A0);E.put(a6.mem,YR_B+512+288,h9);b6=d1^1;i1=(LO1+32+b6*B1)if b6<H1 else(UP1+32+(b6-H1)*B1);run_k(E,rprog,rlab,"init3",dict(x0=1,x1=8*b6,x2=i1,x3=B1),a6,calls,cert.INIT3_BKT+1);p0=Port1(g0.mem);bv=E.a64.a64_run(rprog,rlab,"xemit",dict(x1=YR_B,x2=101,x26=3<<20,x3=K1,x4=32,x6=e2<<49,x8=b5,sp=g0.SP,x30=E.SNT),g0.mem,dev=p0.dev);p1=Port1(g1.mem);n1=E.a64.a64_run(rprog,rlab,"xemit",dict(x1=YR_B,x2=202,x26=9<<20,x3=1,x4=32,x6=e2<<49,x8=b5,sp=g1.SP,x30=E.SNT),g1.mem,dev=p1.dev);p2=Port1(a6.mem);n2a=E.a64.a64_run(rprog,rlab,"xemit",dict(x1=YR_B,x2=303,x26=0,x3=1,x4=32,x6=i0<<49,x8=b5,sp=a6.SP,x30=E.SNT),a6.mem,dev=p2.dev);n2b=E.a64.a64_run(rprog,rlab,"xemit",dict(x1=YR_B+512,x2=303,x26=4,x3=1,x4=32,x6=e2<<49,x8=b5,sp=a6.SP,x30=E.SNT),a6.mem,dev=p2.dev);ok=(bv==K1*cert.XEMIT_HOT+cert.FLUSH1_X+1 and n1==cert.XEMIT_HOT+1 and n2a==cert.XEMIT_HOT+1 and n2b==cert.XEMIT_HOT+1 and len(p0.blocks)==1);run_k(E,rprog,rlab,"fin1",dict(x7=8*d1,x1=1),g0,calls,1+cert.FIN1_X,p0);run_k(E,rprog,rlab,"fin1",dict(x7=8*d1,x1=1),g1,calls,1+cert.FIN1_X,p1);run_k(E,rprog,rlab,"fin1",dict(x7=8*d1,x1=1),a6,calls,1+cert.FIN1_X,p2);a7=Node.R2+4*32*K2B;bw=sorted({d2}|{(g0.mem[(LO1 if d1<H1 else UP1)+(d1 if d1<H1 else d1-H1)*B1+32+32*u])&M41 for u in range(K1)})
 for i7 in bw:rl.mem[LC2+8*i7]=blk2(i7)+32;rl.mem[(LC2+8*i7)|RC2]=Node.R2
 g4=Port2(rl.mem);i8=Sender(lambda d:Node.R2,a7,cert,Node.E2);ae=p0.blocks+p1.blocks+p2.blocks
 for m in ae:
  e0=i8.extra
  for t in range(m[1]//32):
   ab=list(m[4+4*t:8+4*t]);i8.put(ab[0]&M41,ab)
  deliver_arr(E,rprog,rlab,rl,m,a7,i8.extra-e0,calls,cert,g4,i8.tag,3)
 e3=Port2(rl.mem)
 for i7 in bw:run_k(E,rprog,rlab,"fin",dict(x0=MB2,x7=LC2+8*i7,x3=(LC2+8*i7)|RC2,x8=blk2(i7)+32,x1=1,x2=Node.E2B,x11=i8.tag,x6=a7),rl,calls,i8.fin_count([i7]),e3)
 CC=(NEXT_BASE-asr1(Node.R2))%2**64;hm.mem[Node.E2]=Node.R2;home_state(hm,64,CC);hj=HomeJoin(cert,CC,64)
 for m in g4.blocks+e3.blocks+e3.msgs:
  if dest(m[0])!=d2:continue
  if len(m)==6:deliver(E,rprog,rlab,hm.mem,m,calls,cert,cnt=CNTE)
  else:f,ap=hj.block(m,hm.mem[Node.E2],2**cert.HB);ok&=ap;deliver(E,rprog,rlab,hm.mem,m,calls,cert,fused=True,r2=(Node.R2,a7),formula=f)
 g5=hm.mem[Node.E2B];mh=(g5-Node.R2)//32;ok&=(len(hj.cum)==1 and mh>=2);i2,i3=hj.cum[0];ok&=({hm.mem[i2+16],hm.mem[i3+16]}=={101,202});g6=E.a64.a64_run(rprog,rlab,"hhclr",dict(x1=Node.R2,x2=g5,x10=HT,sp=hm.SP,x30=E.SNT),hm.mem);ok&=(g6==3+6*mh and all(hm.mem[HT+8*b]==0 for b in hj.heads));ok&=all(c==f for _,c,f in calls);bx,j7,AB,B0=34*(1<<40)+0x20000,34*(1<<40)+0x40000,(0x987654321<<192)+0x77700000000,(0x987654321<<192)
 def wrep(m,i,a0,a1,c0,c1):
  o=bx+128*i;E.put(m,o,a0);E.put(m,o+32,a1);E.put(m,o+64,c0,2);E.put(m,o+80,c1,2)
  for w in range(4):m[o+96+8*w]=64*i+w+1
 def rps(p_,l_,n,m):o={};s=E.a64.a64_run(p_,l_,"psel",dict(x0=bx,x1=n,x2=j7,sp=g0.SP,x30=E.SNT),m,out_regs=o);return o["x0"],s
 ay=1<<80
 for i,rv in enumerate([(0,9,1,2),(AB+50,AB+10,1,2),(AB+50,AB+50,1,2),(AB+50,AB+50+(1<<63),1,2),(AB+50,AB+50+(1<<64),1,2),(AB+50,AB+90,ay|7,ay|7),(AB+100,AB+200,10,11),(AB+300,AB+400,11,12),(AB+100,AB+500,12,13),(AB+600,AB+200,12,13),(AB+100+(1<<192),AB+200+(1<<192),ay|10,ay|11),(AB+200+(1<<192),AB+250+(1<<192),14,15),(AB+(2<<192),AB+(2<<192)+((1<<63)-1),14,15),(B0+(3<<192)+(1<<63)-5,B0+(3<<192)+(1<<63)+5,16,17),(AB+(4<<192),AB+(4<<192)+9,18,19)]):wrep(g0.mem,i,*rv)
 ok1,l6=rps(rprog,rlab,15,g0.mem);sh0,_=rps(rprog,rlab,13,g0.mem);z0,_=rps(rprog,rlab,0,g0.mem);sa=[E.get(g0.mem,j7+128*q+32*s)for q in range(4)for s in(0,1)];sc=[E.get(g0.mem,j7+128*q+64+16*s,2)for q in range(4)for s in(0,1)];j3=Node(dd(int))
 for q in range(3):wrep(j3.mem,q,B0+((q+1)<<192)+10,B0+((q+1)<<192)+20,2*q+1,2*q+2)
 for i in range(3,1030):wrep(j3.mem,i,B0+(3<<192)+15,B0+(3<<192)+20,999+2*i,1000+2*i)
 wc0,k5=rps(rprog,rlab,1030,j3.mem);mp,ml=E.a64.a64_load(open(os.path.join(E.HERE,"kernels","route.s")).read().replace("tbnz x12,#63,5f\n",""));gm=Node(dd(int))
 for q in range(4):wrep(gm.mem,q,AB+(q<<192)+10,AB+(q<<192)+10+(1<<63),2*q+1,2*q+2)
 ok&=(ok1==1 and sh0==0 and z0==0 and wc0==0 and rps(mp,ml,4,gm.mem)[0]==1 and rps(rprog,rlab,4,gm.mem)[0]==0 and k5<=12+220*1024<=cert.PSEL_POOL and((sa[4]>>63)!=(sa[5]>>63))and((sa[6]>>63)!=(sa[7]>>63))and len(set(sa))==8 and sc==[10,11,ay|10,ay|11,14,15,16,17]and all(sa[2*q]>0 and 0<sa[2*q+1]-sa[2*q]<(1<<63)for q in range(4)));TRACE.extend([("G",g0.mem),("G",g1.mem),("G",a6.mem),("G",j3.mem),("R",rl.mem),("H",hm.mem)]);return dict(ok=ok,xemit_64_insns=bv,hhclr_insns=g6,home_records=mh,candidates=len(hj.cum),psel_ok_steps=l6,psel_wc_steps=k5)
def wait_dispatch_checks(E,rt,rg,cert):
 rprog,rlab=rt[0],rt[1];nd=Node(dd(int));b1=LINK_GEN+0x10000;b7=[1,1,1]
 def wake_cnt(mem,ak,ev):
  assert b7,"unexpected extra WFE in wcnt";mem[b1]+=b7.pop(0)
  if ak[0]==(b1>>11):ak[0]=None;ev[0]=1
 k6=E.a64.a64_run(rprog,rlab,"wcnt",dict(x0=b1,x1=0,x30=E.SNT),nd.mem);k7=E.a64.a64_run(rprog,rlab,"wcnt",dict(x0=b1,x1=3,x30=E.SNT),nd.mem,wfe=wake_cnt);b7=[1];nd.mem[b1]=0;l7=E.a64.a64_run(rprog,rlab,"wcnt",dict(x0=b1,x1=1,x30=E.SNT),nd.mem,wfe=wake_cnt,ev0=1);root,bz=Node(dd(int)),[Node(dd(int))for _ in range(3)];i5=LINK_GEN+0x12000;b3=True;g7=early=0
 for ep in(1,2,3):
  ca=3-early;ae=[1]*ca;early=0
  def wake_r(mem,ak,ev):
   assert ae;mem[i5]+=ae.pop(0)
   if ak[0]==(i5>>11):ak[0]=None;ev[0]=1
  nr=E.a64.a64_run(rprog,rlab,"wcnt",dict(x0=i5,x1=3*ep,x30=E.SNT),root.mem,wfe=wake_r,ev0=int(ep>1));g7+=nr;b3&=(nr<=4+5+ca*cert.WAKE and root.mem[i5]==3*ep and not ae)
  for wi,w in enumerate(bz):
   wk=[1]
   def wake_w(mem,ak,ev):
    assert wk;mem[i5]+=wk.pop(0)
    if ak[0]==(i5>>11):ak[0]=None;ev[0]=1
   nw=E.a64.a64_run(rprog,rlab,"wcnt",dict(x0=i5,x1=ep,x30=E.SNT),w.mem,wfe=wake_w);b3&=(nw==4+cert.WAKE and w.mem[i5]==ep and not wk)
   if ep==1 and wi==0 or ep==2:root.mem[i5]+=1;early+=1
 b3&=(g7<=3*(4+5)+9*cert.WAKE and root.mem[i5]==9);ok=(k6==4 and k7==4+3*cert.WAKE and l7==4+2*cert.WAKE and nd.mem[b1]==1 and not b7 and b3);TRACE.extend([("G",nd.mem),("G",root.mem)]+[("G",w.mem)for w in bz]);return dict(ok=ok,wcnt_0=k6,wcnt_3=k7,clust_bar_epochs_ok=b3,wake=cert.WAKE)

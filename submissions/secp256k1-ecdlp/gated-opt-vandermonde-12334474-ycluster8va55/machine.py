import os,random;from collections import defaultdict as dd;K1,K2B=64,128;B1,B2=64*K1,64*K2B;UP1=28672*(1<<40);UP=(1<<64)-16384*(1<<40);LO1,H1=24576*(1<<40),1<<40;LC2,RC2=1<<45,1<<46;MB2=(1<<48)-128;PW=(1<<48)-168;CNT=(1<<55)-8192;CNTE=(1<<55)-6144;DONE=(1<<48)-4096;ST=(1<<55)-160;HT=(1<<64)-(1<<55);NEXT_BASE=(1<<64)-21846*(1<<40)
def s64(v):return v-(1<<64)if(v&(1<<63))else(v&((1<<64)-1))
def asr1(v):return(s64(v)>>1)&((1<<64)-1)
LINK=(1<<55)-(1<<40);LINK_MSG=LINK+0x8000;WINDOW=(LINK,LINK+0x8000);LINK_GEN=34<<40;LINK_GEN_MSG=LINK_GEN+0x1000;PORTY_BASE=LINK_GEN+0xB000;STALE=0x5A5A5A5A5A5A5A5A;M41=2**41-1;DM=(2**64-1)^M41;H2S_REGS,DLVH_REGS,_NID={},{},[0]
def blk1(d):return LO1+B1*d if d<H1 else UP1+B1*(d-H1)
def blk2(d):return UP+B2*d
LC1_TABLES=[(0,8*H1,LO1+32,B1),(8*H1,16*H1,UP1+32,B1),(LO1,LO1+B1*H1,0,0),(UP1,UP1+B1*H1,0,0)]
def dest(word):return word>>3&M41
def dest1(word):return word>>3
class Port1:
 def __init__(self,mem):self.mem,self.blocks,self.pend=mem,[],[];self.dev=(0,1<<64,self.store,self.read,lambda:not self.pend)
 def store(self,a,v):
  if a%B1==0 and(LO1<=a<LO1+B1*H1 or UP1<=a<UP1+B1*H1):self.pend.append((a,v))
 def read(self):
  for a,v in self.pend:self.blocks.append([v]+[self.mem.get(a+8*w,0)for w in(1,2,3)]+[self.mem.get(a+32+8*w,STALE)for w in range(4*K1)])
  self.pend=[]
class Port2:
 def __init__(self,mem,mb=MB2):self.mem,self.msgs,self.blocks,self.mb,self.pend=mem,[],[],mb,[];self.dev=(0,1<<64,self.store,self.read,lambda:not self.pend)
 def store(self,a,v):
  mb=self.mb
  if a==mb+40 or(not mb<=a<mb+48 and UP<=a<UP+B2*(1<<41)and(a-UP)%B2==0):self.pend.append((a,v))
 def read(self):
  for a,v in self.pend:
   if a==self.mb+40:self.msgs.append([self.mem.get(self.mb+8*w,0)for w in range(5)]+[v])
   else:self.blocks.append([v,self.mem.get(a+8,0),self.mem.get(a+16,0),self.mem.get(a+24,0)]+[self.mem.get(a+32+8*w,STALE)for w in range(4*K2B)])
  self.pend=[]
class LinkView:
 def __init__(self,mem,msg,writes,lbase=LINK):self.mem,self.msg,self.writes,self.lbase=mem,msg,writes,lbase;self.base=lbase+32 if lbase==LINK and len(msg)>6 else lbase
 def __getitem__(self,a):return self.msg[(a-self.base)//8]if self.base<=a<self.base+8*len(self.msg)else self.mem[a]
 def get(self,a,d=None):return self[a]if self.base<=a<self.base+8*len(self.msg)else self.mem.get(a,d)
 def __setitem__(self,a,v):
  assert not self.lbase<=a<self.lbase+0x8000,"write into link aperture"
  if a not in(CNT,CNTE):self.writes.append(a)
  self.mem[a]=v
DLVH_PORTS=tuple(LINK+64*k for k in range(1,9));H2S_PORTS=tuple(LINK+64*k for k in range(1,7))
def dlvh0(E,o,labels,mem,cnt):r0={};E.a64 .a64_run(o,labels,"dlvh0",dict(x29=cnt,x30=E.SNT),mem,out_regs=r0);assert tuple(r0[r]for r in("x0","x1","x7","x8","x24","x25","x26","x27"))==DLVH_PORTS and r0["x28"]==1;return r0
def deliver(E,I,labels,mem,msg,calls,cert,writes=None,fused=False,r2=None,formula=None,cnt=CNT):
 e,o="dlv"if len(msg)==6 else"dlvh",{};lb=LINK_MSG if e=="dlv"else LINK;mem.setdefault(cnt,0);c0=mem[cnt];J=DLVH_REGS.get(id(mem))or dlvh0(E,I,labels,mem,cnt)if e=="dlvh"else dict(x0=lb,x28=1);regs=dict(J,x19=STALE,x29=cnt,x30=E.SNT);n=E.a64 .a64_run(I,labels,e,regs,LinkView(mem,msg,[]if writes is None else writes,lb),out_regs=o,ls64=([WINDOW]+([r2]if r2 else[]),[r2]if r2 else[]))
 if e=="dlvh":DLVH_REGS[id(mem)]=o;assert tuple(o[r]for r in("x0","x1","x7","x8","x24","x25","x26","x27"))==DLVH_PORTS
 ok=o["x28"]==1 and o["x29"]==cnt and mem[cnt]==c0+1;calls.append((e,n if ok else-1,cert.DLV+cert.RCV if e=="dlv"else formula))
def h2s_state(o):return o["x14"]<<64|o["x13"],o["x5"]<<64|o["x15"]
def h2s0(E,I,labels,J,end,o):r0={};E.a64 .a64_run(I,labels,"h2s0",dict(x12=end,x0=o>>64,x8=o%2**64,x30=E.SNT),J.mem,out_regs=r0);assert(r0["x3"],r0["x10"],r0["x18"],r0["x12"],r0["x29"],r0["x6"],r0["x28"])==(PW,LC2,32*K2B,end,DONE,1<<62,0);assert tuple(r0[r]for r in("x0","x1","x2","x7","x8","x16"))==H2S_PORTS;assert h2s_state(r0)==(1<<126,o)and(J.mem[PW-8]<<64|J.mem[PW-16])==o;return r0
def deliver_arr(E,a4,labels,I,msg,end,extra,calls,cert,a5,tag=0,a6=1<<66):L,o=msg[1]//32,{};r0=H2S_REGS.get((I.uid,end))or h2s0(E,a4,labels,I,end,a6);I.mem[PW]=tag;I.mem.setdefault(DONE,0);d0,(D0,f0)=I.mem[DONE],h2s_state(r0);regs=dict(r0,x19=STALE,sp=I.SP,x30=E.SNT);n=E.a64 .a64_run(a4,labels,"h2s",regs,LinkView(I.mem,msg,[]),out_regs=o,dev=a5 and a5 .dev,ls64=([WINDOW],[]));assert tuple(o[r]for r in("x0","x1","x2","x7","x8","x16"))==H2S_PORTS;H2S_REGS[(I.uid,end)]=o;fin=L<K1;a7=D0+(msg[2]if fin else 0)-(1<<126 if fin and f0==1 else 0)-1;J=a7==0;assert I.mem[DONE]==d0+J and h2s_state(o)==((1<<126,a6)if J else(a7,f0-fin));calls.append(("h2s",n,cert.h2s_count(L)+extra+(cert.FIN_CNT+(f0==1)if fin else 0)+cert.RCV1+cert.H2S_DONE*J));return J
class Sender1:
 def __init__(self,cert):self.cert,self.n,self.slots,self.cnt,self.msgs,self.seq,self.extra=cert,{},{},{},[],[],0
 def block(self,d,n,c):return[8*d,32*n,c,0]+[w for r in self.slots.get(d,[None]*K1)for w in(r or[STALE]*4)]
 def put(self,d,o):
  s,n=self.slots.setdefault(d,[None]*K1),self.n.get(d,0);s[n]=o;self.seq.append((d,o));self.n[d]=(n+1)%K1
  if n+1==K1:self.cnt[d]=self.cnt.get(d,0)+1;self.msgs.append(self.block(d,K1,self.cnt[d]));self.extra+=self.cert.FLUSH1_X
 def finals(self,o):
  out=[]
  for d in o:out.append(self.block(d,self.n.get(d,0),self.cnt.get(d,0)+1));self.cnt[d],self.n[d]=0,0
  return out
def state1_ok(mem,o):
 for d,s in o.slots.items():
  b=blk1(d)
  if mem[8*d]!=b+32+32*o.n.get(d,0)or mem[b+16]!=o.cnt.get(d,0):return False
  if any([mem[b+32*t+32+8*w]for w in range(4)]!=s[t]for t in range(K1)if s[t]is not None):return False
 return True
class Sender:
 def __init__(self,rc0,end,cert,tag=0):self.rc0,self.end,self.cert,self.tag=rc0,end,cert,tag;self.n,self.rc,self.slots,self.msgs,self.seq,self.sent,self.drops,self.extra={},{},{},[],[],{},0,0
 def cursor(self,d):return self.rc.get(d,self.rc0(d))
 def put(self,d,o):
  s,n=self.slots.setdefault(d,[None]*K2B),self.n.get(d,0);s[n]=o;self.seq.append((d,o));self.n[d]=(n+1)%K2B
  if n+1<K2B:return
  c=self.cursor(d);self.rc[d]=(c+32*K2B)&M64;self.extra+=self.cert.FLUSH2_X if s64(c)<s64(self.end)else self.cert.FLUSH2_D
  if s64(c)<s64(self.end):self.msgs.append([LC2+8*d,c,32*K2B,self.tag]+[w for r in s for w in r]);self.sent[d]=self.sent.get(d,0)+1
  else:self.drops+=1;self.rc[d]=self.end
 def final(self,d):n,c=self.n.get(d,0),self.cursor(d);msg=[LC2+8*d,c,32*n,self.tag]+[w for r in self.slots[d]for w in(r or[STALE]*4)]if n and s64(c)<s64(self.end)else None;o=(c+32*n)&((1<<64)-1);return msg,(o if s64(o)<s64(self.end)else self.end)
 def fin_count(self,o):
  f=1
  for d in o:n=self.n.get(d,0);f+=self.cert.FINH_DEST+(0 if not n else self.cert.FINH_SENT if s64(self.cursor(d))<s64(self.end)else self.cert.FIN_HELD)
  return f
def keys_reference(E,mem,I,q2,W2,KN,o,cert):
 meter,out=0,[]
 for l in range(q2):
  e=I+E.ENTRY*l;w3=(KN+(l<<2))&E.M64
  for s in(0,1):
   J=E.fracs(mem,e,s);dev=0
   for t in range(3):f=E.record_fields(E.leaf_w(J,t));meter+=(cert.JPATH if f[0]==1 else cert.JPATH2 if f[0]==2 else 0);dev|=bool(f[0]);o.put(f[1],[f[3],f[4],W2,w3]);out.append(f)
   meter+=cert.SJOIN if dev else 0
 return meter,out
def plant(E,mem,sb,I):
 for l,J in enumerate(I):
  for s in(0,1):
   v=[E.t128(E.BT[m]*c%E.p)for c in E.inverse_dft(J[3*s:3*s+3])for m in(0,1)]
   for o,w in zip(E.FRAC_OFF[s],v):E.put(mem,sb+E.ENTRY*l+o,w,2)
def state_ok(mem,o):
 for d,s in o.slots.items():
  b,e=blk2(d),LC2+8*d
  if mem[e]!=b+32+32*o.n.get(d,0)or mem[e|RC2]!=o.cursor(d):return False
  if any([mem[b+32*t+32+8*w]for w in range(4)]!=s[t]for t in range(K2B)if s[t]is not None):return False
 return True
M64=(1<<64)-1
def unpack(I,m,o=None):
 if len(m)==6:
  for w in range(4):I[(m[1]+8*w)&M64]=m[2+w]
  return
 for t in range(m[2]//32):
  for w in range(4):I[(m[1]+32*t+8*w)&M64]=m[4+4*t+w]
  if o is not None:o.append((m[1]+32*t)&M64)
def rule_ok(J,writes,o,finals):
 ok,cnt=True,{}
 for d,I in o.seq:
  t=cnt[d]=cnt.get(d,-1)+1;c=(o.rc0(d)+32*t)&M64
  if s64(c)<s64(o.end)and(finals is not None or t<K2B*o.sent.get(d,0)):ok&=[J[d].get((c+8*w)&M64)for w in range(4)]==I
 if finals is not None:ok&=all(f==(((o.rc0(d)+32*(cnt.get(d,-1)+1))&M64)if s64((o.rc0(d)+32*(cnt.get(d,-1)+1))&M64)<s64(o.end)else o.end)for d,f in finals.items())
 return ok and all(s64(o.rc0(d))<=s64(a)<s64(o.end)for d,w in writes.items()for a in w)
TRACE=[]
def reg_trace(cert):
 o={k:[(lo,hi,b)for r,_,lo,hi,b in cert.REGISTRY if r==k]for k in"GRH"};n=bad=0
 for role,mem in TRACE:
  for a in mem:n+=1;hit=[b for lo,hi,b in o[role]if lo<=a<hi];bad+=len(hit)!=1 or hit[0]=="analytic"
 return dict(nodes=len(TRACE),addresses=n,outside_or_analytic=bad,ok=bad==0 and n>0)
class Node:
 HB,SB,SP=0x100000,(1<<45)+(1<<40)+(1<<34),(34<<40)+0x7000000;R2,E2,E2B=0x3000000,LINK+0x1F000,LINK+0x1F800;CL=LINK+0x200000
 def __init__(self,mem=None,sp=None):_NID[0]+=1;self.uid,self.mem,self.SP=_NID[0],({}if mem is None else mem),(self.SP if sp is None else sp)
def run_k(E,J,labels,entry,regs,nd,calls,formula,I=None,ls64=None):o={};regs=dict(regs,sp=nd.SP,x30=E.SNT);calls.append((entry,E.a64 .a64_run(J,labels,entry,regs,nd.mem,out_regs=o,dev=I and I.dev,ls64=ls64),formula));return o
class HomeJoin:
 def __init__(self,cert,CC,ncl):self.cert,self.CC,self.ncl=cert,CC,ncl;self.chains,self.heads,self.nextd,self.recs,self.cum,self.seen={},{},{},{},[],[];self.nls0,self.lists,self.paths,self.segs=0,[0,0],[0,0],{}
 def block(self,m,F,NB):
  cert,c,cnt=self.cert,m[1],m[2];lo=F if s64(F)>s64(c)else c;hi=(c+cnt)&M64;I,J=(s64(lo)-s64(c))//32,(s64(hi)-s64(lo))//32;a5=s64(F)<=s64(c)and cnt==32*K2B;self.paths[0 if a5 else 1]+=1;f=cert.dlvh_copy(I,J);ok=(f==cert.DLVHB if a5 else f<=cert.DLVHS);k=(I%2,J%2,I>0,J>0);self.segs[k]=self.segs.get(k,0)+1
  for t in range(I,I+J):
   a,o=(c+32*t)&M64,m[4+4*t:8+4*t];a4=not(t%2==0 and t+2<=I+J);self.recs[a]=o;self.seen.append((o[0],o[1],a));b=o[1]&(NB-1);ch=self.chains.setdefault(b,[]);f+=cert.HCHG_NEW
   if ch:
    match=None
    for j,e in enumerate(reversed(ch)):
     f+=cert.HCHG_VISIT
     if e[1]==o[1]:
      f+=cert.HCHG_VER
      if(e[0]&DM)==(o[0]&DM):match=(e,j);break
     if j==len(ch)-1:f+=cert.HCHG_END*(j%2);break
    if match is not None:
     e,j=match;f-=1
     if len(self.cum)<self.ncl:self.cum.append((e[2],a));f+=cert.HCHG_LS1-cert.HCHG_LSA*(not a4);self.lists[a4]+=1
     else:self.nls0+=1;f+=cert.HCHG_LS0
   self.nextd[a]=self.heads.get(b,0);self.heads[b]=(asr1(a)+self.CC)&(2**64-1);ch.append((o[0],o[1],a))
  return f+cert.RCV,ok
def home_state(o,ncl,CC):DLVH_REGS.pop(id(o.mem),None);o.mem[ST],o.mem[ST+8],o.mem[ST+16],o.mem[ST+24]=o.CL,o.CL+16*ncl,CC,-2*CC&(2**64-1)
def part_b(E,cert):
 E.a64 .DEV[:]=[(LINK,LINK+65536),(LINK_GEN,LINK_GEN+65536)];rg=random.Random(E.SEED+32452843);J=E.a64 .a64_load(open(os.path.join(E.HERE,"kernels","route.s")).read());lk=E.a64 .a64_load(open(os.path.join(E.HERE,"kernels","leaf.s")).read())+("keys",);H2=1<<41;ag,am,G,q2=(0,1,H1,H1+1),(0,1,H2-2,H2-1),4,16;calls=[];ae={g:Node(dd(int))for g in range(G)};ap={d1:Node(sp=0x110000)for d1 in ag};homes={(d1,d2):Node(dd(int),sp=LINK+0x110000)for d1 in ag for d2 in am};pool=[E.planted_x(rg,rg.choice(ag),rg.choice(am))[1]for _ in range(10)];batches,shared={},0
 for g in ae:
  ar=[]
  for _ in range(6*q2):
   if rg.random()<0.15:ar.append(E.planted_x(rg,0,0,rg.choice(pool))[0]);shared+=1
   else:ar.append(E.planted_x(rg,rg.choice(ag),rg.choice(am))[0])
  batches[g]=[ar[6*l:6*l+6]for l in range(q2)]
 aa,ak,blocks=True,True,0;af={d1:[]for d1 in ag}
 for g,o in ae.items():
  sb=o.SB-E.ENTRY*q2
  for t0,n,first in((0,2,LO1+32),(8*H1,2,UP1+32)):run_k(E,J[0],J[1],"init3",dict(x0=n,x1=t0,x2=first,x3=B1),o,calls,cert.INIT3_BKT*n+1)
  plant(E,o.mem,sb,batches[g]);KN=(g>>64)|(5+g)<<26;a4=Sender1(cert);meter,fl=keys_reference(E,o.mem,sb,q2,g,KN,a4,cert);ak&=[f[1:]for f in fl]==[E.class_label(x)for ao in batches[g]for x in ao];I=Port1(o.mem);run_k(E,lk[0],lk[1],"keys",dict(x2=g,x3=sb+E.KEYS0,x26=KN),o,calls,cert.KEY_FIX+cert.KEY_INC*q2+meter+a4 .extra,I,ls64=([(sb,o.SB)],[]));ok=I.blocks==a4 .msgs and state1_ok(o.mem,a4);nb=len(I.blocks)
  for t0,n in((0,2),(8*H1,2)):run_k(E,J[0],J[1],"fin1",dict(x7=t0,x1=n),o,calls,1+cert.FIN1_X*n,I)
  ai=a4 .finals(ag);ok&=I.blocks[nb:]==ai and all(o.mem[8*d]==blk1(d)+32 and o.mem[blk1(d)+16]==0 for d in ag);aa&=ok;blocks+=len(I.blocks)
  for m in I.blocks:af[dest1(m[0])].append(m)
 a8=Node.R2+3*32*K2B;aq={h:[]for h in homes};a6={h:[]for h in homes};ac=True
 for d1,o in ap.items():
  for t0,n in((0,2),(H2-2,2)):run_k(E,J[0],J[1],"init3",dict(x0=n,x1=LC2+8*t0,x2=blk2(t0)+32,x3=B2),o,calls,cert.INIT3_BKT*n+1);run_k(E,J[0],J[1],"tinit",dict(x0=n,x1=(LC2+8*t0)|RC2,x2=Node.R2),o,calls,cert.TINIT_ENTRY*n+1)
  a4=Sender(lambda d:Node.R2,a8,cert,Node.E2);I=Port2(o.mem);al=af[d1];rg.shuffle(al);dn=[]
  for m in al:
   e0=a4 .extra
   for t in range(m[1]//32):a9=list(m[4+4*t:8+4*t]);a4 .put(a9[0]&M41,a9);aq[(d1,a9[0]&M41)].append(tuple(a9))
   dn.append(deliver_arr(E,J[0],J[1],o,m,a8,a4 .extra-e0,calls,cert,I,a4 .tag,G))
  ac&=(dn==[False]*(len(al)-1)+[True]and o.mem[DONE]==1 and h2s_state(H2S_REGS[(o.uid,a8)])==(1<<126,G)and sum(m[2]for m in al if m[1]<32*K1)==len(al));ok=I.blocks==a4 .msgs and not I.msgs and state_ok(o.mem,a4);a7=Port2(o.mem)
  for t0,n in((0,2),(H2-2,2)):run_k(E,J[0],J[1],"fin",dict(x0=MB2,x7=LC2+8*t0,x3=(LC2+8*t0)|RC2,x8=blk2(t0)+32,x1=n,x2=Node.E2B,x11=a4 .tag,x6=a8),o,calls,a4 .fin_count(range(t0,t0+n)),a7)
  ai=[a4 .final(d)for d in am];ok&=(a7 .blocks==[m for m,_ in ai if m]and a7 .msgs==[[LC2+8*d,Node.E2B,c,0,0,0]for d,(_,c)in zip(am,ai)]);aa&=ok;blocks+=len(I.blocks)+len(a7 .blocks)
  for m in I.blocks+a7 .blocks+a7 .msgs:a6[(d1,dest(m[0]))].append(m)
 ad,ah,ab,aj=True,True,True,0;CC=(NEXT_BASE-asr1(Node.R2))%2**64
 for h,o in homes.items():
  o.mem[Node.E2]=Node.R2;home_state(o,64,CC);hj=HomeJoin(cert,CC,64)
  for m in a6[h]:
   if len(m)==6:deliver(E,J[0],J[1],o.mem,m,calls,cert,cnt=CNTE)
   else:
    f,ax=hj.block(m,o.mem[Node.E2],2**cert.HB);ab&=ax;deliver(E,J[0],J[1],o.mem,m,calls,cert,fused=True,r2=(Node.R2,a8),formula=f)
  F,F0=o.mem[Node.E2B],Node.R2;start=F0-(F0-Node.R2)%(32*K2B);aw=1+ -(-(F-start)//(32*K2B));ac&=o.mem[CNTE]==1 and o.mem[CNT]+1==len(a6[h])==aw;au=[tuple(o.mem[x+8*k]for k in range(4))for x in range(Node.R2,F,32)];ad&=au==aq[h][:len(au)]and F==min(Node.R2+32*len(aq[h]),a8);ah&=all(a9[0]&M41==h[1]and a9[2]in ae and(a9[3]>>26)==5+a9[2]and 0<=((a9[3]>>2)&((1<<24)-1))<q2 for a9 in au);ab&=([(o.mem[o.CL+16*i],o.mem[o.CL+16*i+8])for i in range(len(hj.cum))]==hj.cum and all(o.mem[HT+8*b]==v for b,v in hj.heads.items()));aj+=len(hj.cum)
 groups={}
 for g in ae:
  for ao in batches[g]:
   for x in ao:groups.setdefault(E.exk(x)[0],set()).add(g)
 cnt,av={},{}
 for e,c,f in calls:cnt[e],av[e]=cnt.get(e,0)+c,av.get(e,0)+f
 a5=dict(leaves=6*q2*G,home_candidates=aj,counts=cnt,formulas_ok=all(c==f for _,c,f in calls)and ak and aa and ad and ah and ab and ac and aj>0);a5["block_paths"]=block_paths(E,J,lk,rg,cert);a5["grouped_paths"]=grouped_paths(E,J[0],J[1],rg,cert);a5["ls64"]=ls64_checks(E,rg);a5["va55"]=va55_checks(E);a5["encoding"]=encoding_checks(E);a5["phase4_round"]=phase4_round(E,J,rg,cert);a5["y_ring"]=yring_checks(E,J,rg,cert);a5["gate"]=gate_checks(E,J,cert);a5["wait_dispatch"]=wait_dispatch_checks(E,J,rg,cert);a5["phase0"]=phase0_xemit_hhclr(E,J,rg,cert);TRACE.extend([("G",n.mem)for n in ae.values()]+[("R",n.mem)for n in ap.values()]+[("H",n.mem)for n in homes.values()]);a5["registry_trace"]=reg_trace(cert);return a5
class PortY:
 def __init__(self,mem,base=PORTY_BASE):self.mem,self.base,self.msgs,self.pend=mem,base,[],[];self.dev=(0,1<<64,self.store,self.read,lambda:not self.pend)
 def store(self,a,v):
  if a==self.base:self.pend.append(v)
 def read(self):
  for v in self.pend:self.msgs.append([v]+[0]*7+[self.mem.get(self.base+64+8*w,0)for w in range(56)])
  self.pend=[]
def yring_checks(E,rt,rg,cert):
 rprog,rlab=rt[0],rt[1];a4,CNTY=(1<<45)+(257<<31),LINK_GEN+0x12000;nodes=[Node(dd(int))for _ in range(4)];a5=[[rg.getrandbits(64)for _ in range(56)]for _ in range(3)]
 for u,pl in enumerate(a5):
  for w,v in enumerate(pl):nodes[0].mem[a4+512*u+64+8*w]=v
 I=[]
 for g,nd in enumerate(nodes):o={};n0=E.a64 .a64_run(rprog,rlab,"ydlv0",dict(x15=a4+64,x29=CNTY,x30=E.SNT),nd.mem,out_regs=o);assert n0==19;I.append(o)
 a7=([(LINK_GEN,LINK_GEN+512),(a4,a4+4096)],[(a4,a4+4096),(PORTY_BASE,PORTY_BASE+512)]);J=True
 for u in range(3):
  p0=PortY(nodes[0].mem);o0={};aa=E.a64 .a64_run(rprog,rlab,"ysnd",dict(I[0],x30=E.SNT),nodes[0].mem,out_regs=o0,dev=p0 .dev,ls64=a7);I[0]=o0;J&=(aa==cert.YSND and len(p0 .msgs)==1);msg=p0 .msgs[0]
  for g in(1,2):pg=PortY(nodes[g].mem);og={};ab=E.a64 .a64_run(rprog,rlab,"ydlv",dict(I[g],x19=STALE,x30=E.SNT),LinkView(nodes[g].mem,msg,[],LINK_GEN),out_regs=og,dev=pg.dev,ls64=a7);I[g]=og;J&=(ab==cert.YDLV and len(pg.msgs)==1);msg=pg.msgs[0]
  oe={};ac=E.a64 .a64_run(rprog,rlab,"ydlve",dict(I[3],x19=STALE,x30=E.SNT),LinkView(nodes[3].mem,msg,[],LINK_GEN),out_regs=oe,ls64=a7);I[3]=oe;J&=(ac==cert.YDLVE)
 a8=all([nodes[g].mem[a4+512*u+64+8*w]for w in range(56)]==a5[u]and nodes[g].mem[a4+512*u]==0 for g in(1,2,3)for u in range(3));a6=all(nodes[g].mem[CNTY]==3 for g in(1,2,3));a9=open(os.path.join(E.HERE,"kernels","route.s")).read().replace("ysnd:\n","ysnd_orig:\n")+"\nysnd_nodmb:\n\tstlr x28,[x14]\n\tret\n";bp,bl=E.a64 .a64_load(a9);pb=PortY(dd(int))
 try:
  E.a64 .a64_run(bp,bl,"ysnd_nodmb",dict(x14=PORTY_BASE,x28=1,x30=E.SNT),pb.mem,dev=pb.dev);dmb_rejected=False
 except AssertionError:dmb_rejected=True
 return dict(ok=J and a8 and a6 and dmb_rejected,counts=[cert.YSND,cert.YDLV,cert.YDLVE],dmb_rejected=dmb_rejected)
def phase4_round(E,rt,rg,c):
 C,ok,N,BR,BH,a7=[],True,Node,(1<<48)+0x10000,LINK+0x20000,LINK+0x18000;d1,d2,gs=rg.randrange(1<<41),rg.randrange(1<<41),[rg.randrange(3<<64,1<<66),rg.randrange(1<<64)];h,J=N(dd(int),sp=LINK+0x110000),N(dd(int),sp=0x110000)
 def run(I,e,regs,f,base=None,msg=None,lb=LINK_MSG):p=base and Port2(I.mem,base);C.append((E.a64 .a64_run(rt[0],rt[1],e,dict(regs,sp=I.SP,x19=STALE,x30=E.SNT,x28=1,x29=CNT),I.mem if msg is None else LinkView(I.mem,msg,[],lb),dev=p and p.dev),f));return p and p.msgs
 def q(I,m,cap=c.RELAY_CAP,lb=LINK_MSG+1024,bb=BR,ct=CNT):a5,aa,oc=I.mem[ct],I.mem[ct+2048],I.mem[bb-8];C.append((E.a64 .a64_run(rt[0],rt[1],"dlvq",dict(x0=lb,x9=64,x10=bb-8,x11=bb+64*cap,sp=I.SP,x19=STALE,x30=E.SNT,x28=1,x29=ct),LinkView(I.mem,m,[],lb)),(c.DLVQ if oc<bb+64*cap else c.DLVQ_DROP)+(c.DLVQ_CARRY if a5==(1<<64)-1 else 0)+c.RCV));return(I.mem[ct+2048]<<64|I.mem[ct])==((aa<<64|a5)+1)
 R=[[rg.randrange(1<<64),0,g%2**64,g>>64|rg.randrange(1<<62)<<2]for g in gs]
 for i in range(8):h.mem[N.R2+8*i]=R[i//4][i%4]
 h.mem[N.CL],h.mem[N.CL+8],J.mem[BR-8],hd=N.R2,N.R2+32,BR,2**45+8*d2;m1=run(h,"reqh",dict(x0=a7,x12=N.CL,x9=1,x10=hd,x13=BH),c.REQH_C+1,a7);ok&=m1==[[0,0,r[2],r[3],hd,BH+256*i]for i,r in enumerate(R)]and all(q(J,m)for m in m1);m2=run(J,"fwdq",dict(x0=MB2+64,x1=BR,x6=BR+128,x7=8*d1),2*c.FWDQ+3,MB2+64);a4=[[rg.randrange(1<<64)for _ in range(32)]for _ in gs];BG=96*(1<<40);a8=LINK_GEN+0x18000;a9=LINK_GEN+0x10000
 for i,m in enumerate(m2):
  g=N(dd(int));g.mem[BG-8]=BG;TRACE.append(("G",g.mem));ok&=(m[1]&3)<<64|m[0]==gs[i]and q(g,m,cap=c.REQ_CAP,lb=LINK_GEN_MSG+1024,bb=BG,ct=a9)and[g.mem[BG+8*k]for k in range(4)]==[gs[i]>>64,R[i][3],BR+64*i<<8,8*d1]
  for k in range(32):g.mem[BG+4096+8*k]=a4[i][k]
  m3=run(g,"sendr",dict(x0=a8,x1=BG+4096,x9=8,x7=8*d1,x8=g.mem[BG+16]),8*c.SENDR+1,a8);ok&=len(m3)==8 and all(q(J,m)for m in m3)
 m4=run(J,"fwdr",dict(x0=MB2,x1=BR+128,x6=BR+64*18),16*c.FWDR+3,MB2)
 for m in m4:run(h,"dlv",dict(x0=LINK_MSG),c.DLV+c.RCV,msg=m)
 J.mem[BR-8]=BR;d=N(dd(int),sp=0x110000);d.mem[BR-8]=BR+64*c.RELAY_CAP;w=N(dd(int),sp=0x110000);w.mem[BR-8]=BR+64*(c.RELAY_CAP-1);w.mem[CNT]=(1<<64)-2;cw=N(dd(int),sp=0x110000);cw.mem[BR-8]=BR+64*c.RELAY_CAP;cw.mem[CNT]=(1<<64)-1;a6=[];E.a64 .a64_run(rt[0],rt[1],"dlvq",dict(x0=LINK_MSG+1024,x9=64,x10=BR-8,x11=BR+64*c.RELAY_CAP,sp=cw.SP,x19=STALE,x30=E.SNT,x28=1,x29=CNT),LinkView(cw.mem,[5]*6,[],LINK_MSG+1024),dev=(CNT,4096,lambda o,v:a6 .append((cw.mem[CNT+2048],cw.mem[CNT])),lambda:None,lambda:True));ok&=([h.mem[BH+8*k]for k in range(64)]==a4[0]+a4[1]and len(m4)==16 and q(d,[1]*6)and q(d,[1]*6)and d.mem[BR-8]==BR+64*c.RELAY_CAP and d.mem[BR]==0 and q(w,[2]*6)and q(w,[3]*6)and q(w,[4]*6)and w.mem[BR-8]==BR+64*c.RELAY_CAP and w.mem[CNT]==1 and w.mem[CNT+2048]==1 and 0 not in w.mem and a6==[(0,0),(1,0)]and all(x==y for x,y in C)and max(c.PH4_HOPS.values())<=c.MSG);TRACE.extend([("H",h.mem),("R",J.mem),("R",d.mem),("R",w.mem),("R",cw.mem)]);return dict(ok=ok,calls=len(C),hops=c.PH4_HOPS)
def encoding_checks(E):
 a64=E.a64;good=("ldp x1,x2,[x3,#-512]","ldp x1,x2,[x3,#504]","stp x1,x2,[x3],#504","stp x1,x2,[x3,#-512]!","ldr x1,[x3,#32760]","ldr x1,[x3,#-256]","ldr x1,[x3,#255]","str x1,[x3],#-256","add x1,x1,#4095","add x1,x1,#0x1000","tst x3,#0x3fffe0000","mov sp,x1","add sp,sp,#512","cset w9,hs","csel x1,x2,x3,lo","mov x1,#0xfff8000000000000","movk x4,#1,lsl #32","stlur x1,[x3,#-256]","ldr x22,[x3,#-152]","and x1,x2,#0x5555555555555555","ldaxr x2,[x0]","blr x30","wfe");bad=("stp x25,x26,[x3],#768","ldp x21,x22,[x3,#-768]","ldp x1,x2,[x3,#508]","ldp x1,x2,[x3,#-520]","ldp x1,x2,[x3,#12]","ldr x22,[x3,#-664]","ldr x1,[x3,#32768]","ldr x1,[x3,#260]","str x1,[x3],#256","ldr x1,[x3,#-257]!","add x1,x1,#4097","tst x3,#0x1c001","and x1,x2,#0","mov x1,#0x12345","movk x4,#1,lsl #8","stlur x1,[x3,#256]","ldp x1,x1,[x3]","ldp x3,x2,[x3],#16","ldadd x1,x2,[x3,#8]","ldaxr sp,[x0]","ldaxr x1,[x0,#8]","blr sp","wfe x0","ubfx x1,x2,#60,#8","bfi x1,x2,#40,#25","ldr x1,[x3,#-264]","add x1,xzr,#5","tst sp,#0x1c000","ldp xzr,xzr,[x3]","ldr sp,[x3]","mov xzr,#0x3ffe0","movk sp,#1,lsl #48","tbz x3,#64,1f","frob x1,x2","mul x1,x2,#1","adc x1,x2,#3","csel x1,x2,#3,eq","madd x1,x2,x3,#4","cset x1,x2","mvn x1,#2","mul w1,x2,x3","dmb sy","umulh w1,w2,w3","mvn x1,x2,lsl #64","movk x1,#1,lsr #16","stadd x1,[x2,#0]","ld64b x0,[x2,#0]","ldp x1,x2,[x3]!","ld64b x21,[x1]","ld64b x24,[x1]","add x1,x2,#1_0")
 def verdict(o):
  try:a64 .a64_load("f:\n\t"+o+"\n\tret\n");return True
  except AssertionError as ex:assert"not encodable"in str(ex),ex;return False
 out=dict(good_accepted=sum(map(verdict,good)),bad_rejected=sum(not verdict(i)for i in bad),good=len(good),bad=len(bad));out["ok"]=out["good_accepted"]==len(good)and out["bad_rejected"]==len(bad);return out
def va55_checks(E):
 a64=E.a64;I,J=a64 .a64_load("ld:\n\tldr x2,[x1]\n\tret\n");out={};U=2**64-2**55
 for a in(0,2**55-8,U,2**64-8,2**55-4,2**55,2**56-8,2**60,U-8,2**64-2**56,2**63,2**64-4):
  try:a64 .a64_run(I,J,"ld",dict(x1=a,x30=E.SNT),dd(int));out[hex(a)]="ran"
  except AssertionError as ex:out[hex(a)]="rejected"if"55-bit"in str(ex)else"other: "+str(ex)
 v=list(out.values());out["ok"]=v==["ran"]*4+["rejected"]*8;return out
def ls64_checks(E,rg):
 a64=E.a64;a4=open(os.path.join(E.HERE,"kernels","route.s")).read();I,J=a64 .a64_load("ld:\n\tld64b x20,[x1]\n\tret\nst:\n\tst64b x16,[x2]\n\tret\novl:\n\tld64b x12,[x15]\n\tret\n");R2=(Node.R2,Node.R2+4*32*K2B);msg=[rg.getrandbits(64)for _ in range(4+4*K2B)]
 def run(a5,labels,entry,regs,mem,lsr=([WINDOW,R2],[R2])):
  try:a64 .a64_run(a5,labels,entry,dict(regs,x30=E.SNT),mem,ls64=lsr);return"ran"
  except AssertionError as ex:return"LS64 rejected"if"LS64"in str(ex)else"other: "+str(ex)
 out,o={},{};a64 .a64_run(I,J,"ld",dict(x1=LINK+64,x30=E.SNT),LinkView({},msg,[]),out_regs=o,ls64=([WINDOW],[]));out["ld_win01"]=[o["x%d"%(20+k)]for k in range(8)]==msg[4:12];mem={};out["st64b_r2_aligned"]=run(I,J,"st",dict(x2=R2[0]),mem);out["st64b_r2_written"]=all(mem.get(R2[0]+8*k)==0 for k in range(8));out["ld64b_odd_slot"]=run(I,J,"ld",dict(x1=LINK+32),LinkView({},msg,[]));out["ld_out"]=run(I,J,"ld",dict(x1=0x5000040),{});out["st64b_into_window"]=run(I,J,"st",dict(x2=LINK+64),LinkView({},msg,[]));out["st_r2_mis"]=run(I,J,"st",dict(x2=R2[0]+32),{})
 def load(text):
  try:a64 .a64_load(text);return"loaded"
  except AssertionError as ex:return"LS64 rejected"if"LS64"in str(ex)else"other: "+str(ex)
 out["ld64b_x24"]=load("hi:\n\tld64b x24,[x1]\n\tret\n");out["ld64b_x21"]=load("odd:\n\tld64b x21,[x1]\n\tret\n");out["ld_base_ov"]=run(I,J,"ovl",dict(x15=LINK+64),LinkView({},msg,[]))
 def fused(text,c,r2):a5,labels=a64 .a64_load(text);mem=dd(int);CC=(NEXT_BASE-asr1(r2[0]))%2**64;mem[ST],mem[ST+8],mem[ST+16],mem[ST+24]=Node.CL,Node.CL+1024,CC,-2*CC%2**64;a8=LINK+0x300000;mem[a8]=c;m=[0,c,32*K2B,a8]+msg[4:];dr=dlvh0(E,a5,labels,mem,CNT);r=run(a5,labels,"dlvh",dict(dr,x29=CNT),LinkView(mem,m,[]),([WINDOW,r2],[r2]));return r if r!="ran"else("ran, R2 exact"if all(mem[c+8*w]==m[4+w]for w in range(4*K2B))else"ran, R2 wrong")
 out["dlvh_fresh_block"]=fused(a4,R2[0],R2);a6=(LINK-4096,LINK);out["r2_top"]=fused(a4,a6[0],a6);out["r2_link_rej"]=fused(a4,LINK,(LINK,LINK+512));a7=(R2[0]+32,R2[1]+32);out["r2_32mod64"]=fused(a4,a7[0],a7);ab=a4 .replace("6:\nadd x13,x14,x2,asr #1\nDLVHPR x0\n","6:\nadd x13,x14,x2,asr #1\nsub x12,x0,#32\nDLVHPR x12\n");out["dlvh_odd"]=fused(ab,R2[0],R2)if ab!=a4 else"mutation not applied";hm=a4 .replace("mov x1,x0\n","sub x1,x0,#32\n");hp,hl=a64 .a64_load(hm);regs=dict(x0=H2S_PORTS[0],x1=H2S_PORTS[1],x2=H2S_PORTS[2],x7=H2S_PORTS[3],x8=H2S_PORTS[4],x16=H2S_PORTS[5],x10=LC2,x18=32*K2B,x12=1<<40,x3=PW,x29=DONE,x6=1<<62,x14=1<<62,x5=16);out["h2s_odd"]=run(hp,hl,"h2s",regs,LinkView(dd(int),[0,32*2,1,0]+msg[4:4+4*K1],[]),([WINDOW],[]))if hm!=a4 else"mutation not applied";aa=("ld64b_odd_slot","ld_out","st64b_into_window","st_r2_mis","ld64b_x24","ld64b_x21","ld_base_ov","r2_32mod64","dlvh_odd","h2s_odd");a9=(out["ld_win01"]and out["st64b_r2_aligned"]=="ran"and out["st64b_r2_written"]and out["dlvh_fresh_block"]=="ran, R2 exact"and out["r2_top"]=="ran, R2 exact"and"write into link aperture"in out["r2_link_rej"]and all(out[k]=="LS64 rejected"for k in aa));return dict(ok=a9,rejected=len(aa))
def grouped_paths(E,a6,labels,rg,cert):
 GG,aj,I=3,6,2;NB=1<<aj;a5=32*K2B*128;ak=((300,9),(520,262),(250,7))
 def run_case(ac,ncl=64):
  calls,o=[],Node(dd(int));ai,ae=LINK+0x400000,LINK+0x500000;a7=((1<<64)-8192)if ac else o.R2;E2P=(o.E2,o.E2+0x800000);pf=lambda g,s:ai+8*s if g==0 else E2P[(g-1)%2]+32*s;ah={}
  for w in range(cert.ZW):o.mem[HT+8*w]=0x5A5A5A5A5A5A5A5A
  run_k(E,a6,labels,"zero16384",dict(x0=cert.ZW,x1=HT),o,calls,cert.ZERO_PASS+1);ok=all(o.mem[HT+8*w]==0 for w in range(cert.ZW));run_k(E,a6,labels,"init3",dict(x0=I,x1=ai,x2=a7,x3=a5),o,calls,cert.INIT3_BKT*I+1);run_k(E,a6,labels,"dlvh0",dict(x29=CNT),o,calls,cert.DLVH0);ok&=all(o.mem[ai+8*s]==((a7+a5*s)&M64)for s in range(I));CC=(NEXT_BASE-asr1(a7))%2**64;home_state(o,ncl,CC);J=HomeJoin(cert,CC,ncl);al,stopped=[],False;rc=[(a7+a5*s)&M64 for s in range(I)];F=list(rc);ag,nb=[[None]*K2B for _ in range(I)],[0]*I
  for g in range(GG):
   if stopped:break
   al.append(g);ah[g]=[(rg.getrandbits(64),rg.getrandbits(64))for _ in range(3)]
   for k,(px,py)in enumerate(ah[g]):o.mem[ae+(g*3+k)*16]=px;o.mem[ae+(g*3+k)*16+8]=py
   streams=[]
   for s in range(I):
    am=ak[g][s]+(0 if ac or s else 1);recs=[]
    for _ in range(am):recs.append([rg.getrandbits(64),rg.getrandbits(7)<<57|rg.getrandbits(aj),rg.getrandbits(64),rg.getrandbits(64)])
    if J.seen:aq,ap=J.seen[rg.randrange(len(J.seen))][:2];recs[0][1]=ap;recs[0][0]=(aq&DM)|rg.getrandbits(41)
    recs[1][1]=recs[0][1];recs[1][0]=(recs[0][0]&DM)|rg.getrandbits(41);msgs=[]
    for ao in recs:
     ag[s][nb[s]]=ao;nb[s]+=1
     if nb[s]==K2B:msgs.append([0,rc[s],32*K2B,pf(g,s)]+[w for r in ag[s]for w in r]);rc[s]=(rc[s]+32*K2B)&M64;nb[s]=0
    if nb[s]:msgs.append([0,rc[s],32*nb[s],pf(g,s)]+[w for r in ag[s]for w in r or[STALE]*4])
    streams.append(msgs)
   order=[m for st in streams for m in st];rg.shuffle(order)
   for m in order:s=[t for t in range(I)if pf(g,t)==m[3]][0];f,ar=J.block(m,F[s],NB);ok&=ar;ab=[((1<<64)-8192,1<<64),(0,1<<51)]if ac else[(o.R2,o.R2+I*a5)];n=E.a64 .a64_run(a6,labels,"dlvh",dict(DLVH_REGS.get(id(o.mem))or dlvh0(E,a6,labels,o.mem,CNT),x29=CNT,x30=E.SNT),LinkView(o.mem,m,[]),out_regs=(au:={}),ls64=([WINDOW]+ab,ab));DLVH_REGS[id(o.mem)]=au;calls.append(("dlvh",n,f));ok&=all(o.mem[(m[1]+8*w)&M64]==m[4+w]for w in range(4*K2B))
   ok&=all(o.mem[HT+8*b]==J.heads.get(b,0)for b in range(NB));ok&=all(o.mem[(asr1(a)+CC)&M64]==J.nextd[a]and o.mem[((asr1(a)+CC)&M64)+8]==J.recs[a][1]for a in J.nextd);ok&=[(o.mem[o.CL+16*i],o.mem[o.CL+16*i+8])for i in range(len(J.cum))]==J.cum;ok&=o.mem[ST]==o.CL+16*len(J.cum)and len(J.cum)<=ncl
   for s in range(I):F[s]=(rc[s]+32*nb[s])&M64;deliver(E,a6,labels,o.mem,[0,E2P[g%2]+32*s,F[s],0,0,0],calls,cert)
   ok&=all(o.mem[E2P[g%2]+32*s]==F[s]for s in range(I));a4=[((a7+a5*s+32*k)&M64)for s in range(I)for k in range((s64(F[s])-s64((a7+a5*s)&M64))//32)];ok&=sorted(J.nextd)==sorted(a4);aa=[]
   for b in range(NB):
    u=J.heads.get(b,0)
    while u:a=(2*((u-CC)%2**64))&M64;aa.append(a);u=J.nextd[a]
   ok&=sorted(aa)==sorted(a4)and all([o.mem[(a+8*w)&M64]for w in range(4)]==J.recs[a]for a in a4)
   if ac and g==1:stopped=any(o.mem[a+8]==o.mem[b+8]and o.mem[a]&DM==o.mem[b]&DM for a,b in J.cum)
  ok&=all(o.mem[a+8]==o.mem[b+8]and o.mem[a]&DM==o.mem[b]&DM for a,b in J.cum);a8={}
  def find(u):
   while a8 .get(u,u)!=u:u=a8[u]
   return u
  for a,b in J.cum:a8[find(a)]=find(b)
  ad={}
  for k2,k3,a in J.seen:ad.setdefault((k3,k2&DM),[]).append(a)
  true_pairs=[(r[i],r[i+1])for r in ad.values()for i in range(len(r)-1)]
  if ncl>len(J.cum):ok&=all(find(a)==find(b)for a,b in true_pairs)and len(true_pairs)>0
  ok&=all(c==f_ for _,c,f_ in calls);a9=[(o.mem[ae+16*k],o.mem[ae+16*k+8])for k in range(3)]==ah[0];af=sum(1 for _,k3,_ in J.seen if k3&0x30);ok&=a9 and af>0 and len({k3&(NB-1)for _,k3,_ in J.seen})>16;DLVH_REGS.pop(id(o.mem),None);return dict(ok=ok,groups_run=al,candidates=len(J.cum),stopped=stopped,stop_pair_listed=any(o.mem[a+8]==o.mem[b+8]for a,b in J.cum),true_pairs=len(true_pairs),list_skips=J.nls0,old_x_served=a9,records_with_bucket_bits_4_5=af,fast_slow_blocks=J.paths,listings_direct_late=J.lists)
 a,b,c=run_case(True),run_case(False),run_case(False,ncl=1);return dict(ok=slow_wrap(E,a6,labels,cert)and a["ok"]and b["ok"]and c["ok"]and a["groups_run"]==[0,1]and b["groups_run"]==[0,1,2]and a["stopped"]and a["stop_pair_listed"]and a["candidates"]>0 and b["candidates"]>0 and c["candidates"]==1 and c["list_skips"]>0,old_x_served=a["old_x_served"]and b["old_x_served"]and c["old_x_served"],stop_cands=a["candidates"],full_cands=b["candidates"],skips=c["list_skips"])
def slow_wrap(E,P,L,c):
 m=dd(int);m[Node.E2]=M64-31;N=c.NEXT_BASE;home_state(Node(m),64,N);E.a64 .a64_run(P,L,"dlvh",dict(dlvh0(E,P,L,m,CNT),x29=CNT,x30=E.SNT),LinkView(m,[0,M64-63,96,Node.E2]+[v for j in range(128)for v in(0,10+j,0,0)],[]),ls64=([WINDOW],[(M64-4095,1<<64),(0,8192)]));DLVH_REGS.pop(id(m),None);return[m[c.UP_BASE+88],m[c.UP_BASE+96],m[N-8],m[N+8]]==[(N-16)%2**64,N,11,12]
def block_paths(E,rt,lk,rg,cert):
 calls,out=[],{};J=Node(dd(int));ad=(0,H1,1,2,H1+1);xs=[E.planted_x(rg,d,rg.randrange(4),j=t%3)[0]for t,d in enumerate(ad)];ds=[t for t,n in enumerate((558,200,128,2,0))for _ in range(n)];rg.shuffle(ds);q2=len(ds)//6;leaves=[[xs[t]for t in ds[6*l:6*l+6]]for l in range(q2)];sb=J.SB-E.ENTRY*q2;plant(E,J.mem,sb,leaves)
 for t0,n,first in((0,3,LO1+32),(8*H1,2,UP1+32)):run_k(E,rt[0],rt[1],"init3",dict(x0=n,x1=t0,x2=first,x3=B1),J,calls,cert.INIT3_BKT*n+1)
 I=Sender1(cert);meter,fl=keys_reference(E,J.mem,sb,q2,7,7<<26,I,cert);o=Port1(J.mem);run_k(E,lk[0],lk[1],"keys",dict(x2=7,x3=sb+E.KEYS0,x26=7<<26),J,calls,cert.KEY_FIX+cert.KEY_INC*q2+meter+I.extra,o,ls64=([(sb,J.SB)],[]));ok=o.blocks==I.msgs and state1_ok(J.mem,I);nb=len(o.blocks)
 for t0,n in((0,3),(8*H1,2)):run_k(E,rt[0],rt[1],"fin1",dict(x7=t0,x1=n),J,calls,1+cert.FIN1_X*n,o)
 a9=I.finals(sorted(ad));ok&=o.blocks[nb:]==a9;ok&=[f[:2]for f in fl]==[(t%3,ad[t])for l in range(q2)for t in ds[6*l:6*l+6]];out["keys_fin1"]=dict(records=len(ds),blocks_sent=nb,finals=len(a9),empty_finals=sum(m[1]==0 for m in a9),upper_half_blocks=sum(m[0]>=8*H1 for m in o.blocks),model_ok=ok);J=Node();a5=[rg.getrandbits(41)for _ in range(4)];aa=[K1]*12+[0,1,2,3,7,15,31,32,33,62,63]+[K1]*20;rg.shuffle(aa);dl=[a5[t]for t,n in enumerate((1200,400,200,50))for _ in range(n)];rg.shuffle(dl);end,m=J.R2+2*32*K2B,{a5[0]:1,a5[1]:2,a5[2]:1,a5[3]:1};I=Sender(lambda d:end-32*K2B*m[d],end,cert)
 for d in a5:J.mem[LC2+8*d],J.mem[(LC2+8*d)|RC2]=blk2(d)+32,end-32*K2B*m[d]
 k,ac=0,[]
 for n in aa:recs=[[rg.getrandbits(23)<<41|dl[(k+t)%len(dl)],rg.getrandbits(64),rg.getrandbits(64),rg.getrandbits(64)]for t in range(n)];k+=n;ac.append(recs)
 o=Port2(J.mem);a8=0
 for a4 in ac:
  e0=I.extra
  for af in a4:I.put(af[0]&M41,af)
  a8+=len(a4)<K1;deliver_arr(E,rt[0],rt[1],J,[0,32*len(a4),1,0]+[w for r in a4 for w in r]+[STALE]*4*(K1-len(a4)),end,I.extra-e0,calls,cert,o)
 ab=o.blocks==I.msgs and not o.msgs and state_ok(J.mem,I);ab&=h2s_state(H2S_REGS[(J.uid,end)])==((1<<126)+a8-len(ac),(1<<66)-a8)and J.mem[DONE]==0;dm,wr={d:{}for d in a5},{d:[]for d in a5}
 for msg in o.blocks:unpack(dm[dest(msg[0])],msg,wr[dest(msg[0])])
 ab&=rule_ok(dm,wr,I,None);out["h2s"]=dict(records=k,blocks_sent=len(o.blocks),blocks_dropped=I.drops,end_stored_back=I.cursor(a5[0])==end,model_and_rule_ok=ab);out["two_group_finals"]=tg=two_group_finals(E,rt,rg,cert,calls);out["wide_counts"]=wide_counts(E,rt,cert,calls);a6=Node(dd(int));a7=Port1(a6 .mem);run_k(E,rt[0],rt[1],"init3",dict(x0=1,x1=0,x2=LO1+32,x3=B1),a6,calls,cert.INIT3_BKT+1)
 if(cert.DELTA_F-cert.FIX_BATCH+1)<=cert.DELTA_F-cert.FIX_BATCH:run_k(E,lk[0],lk[1],"keys",dict(x2=7,x3=sb+E.KEYS0,x26=0),a6,calls,0,a7)
 run_k(E,rt[0],rt[1],"fin1",dict(x7=0,x1=1),a6,calls,1+cert.FIN1_X,a7);out["meter_skip_ok"]=(a7 .blocks==Sender1(cert).finals([0])and deliver_arr(E,rt[0],rt[1],Node(),a7 .blocks[0],Node.R2+32*K2B,0,calls,cert,None,0,1));out["cnt_eq"]=all(c==f for _,c,f in calls);ae=(out["cnt_eq"]and ok and ab and out["meter_skip_ok"]and out["keys_fin1"]["empty_finals"]==2 and out["keys_fin1"]["upper_half_blocks"]>2 and out["h2s"]["blocks_sent"]>=1 and out["h2s"]["blocks_dropped"]>=1 and out["h2s"]["end_stored_back"]and tg["ok"]and out["wide_counts"]);return dict(ok=ae,meter_skip_ok=out["meter_skip_ok"],dropped=out["h2s"]["blocks_dropped"])
def wide_counts(E,I,cert,calls):
 o,end=Node(),Node.R2+32*K2B;key=(o.uid,end);run=lambda c:deliver_arr(E,I[0],I[1],o,[0,0,c,0]+[STALE]*4*K1,end,0,calls,cert,None,0,1<<66)
 def put(D,f):H2S_REGS[key]=dict(H2S_REGS[key],x13=D%2**64,x14=D>>64,x15=f%2**64,x5=f>>64)
 ok=not run(1)and h2s_state(H2S_REGS[key])==(1<<126,(1<<66)-1);put((1<<126)+(1<<64),1);ok&=not run(1)and h2s_state(H2S_REGS[key])==(1<<64,0);put((1<<126)-1,1);return ok and run(2)and o.mem[DONE]==1 and h2s_state(H2S_REGS[key])==(1<<126,1<<66)
def two_group_finals(E,rt,rg,cert,calls):
 J=Node();d,base=(1<<41)-1,M64-16383;end=(base+60*32*K2B)&M64;J.mem[LC2+8*d],J.mem[(LC2+8*d)|RC2]=blk2(d)+32,base;o=Sender(lambda _d:base,end,cert);dm,wr,finals,recs={},[],[],[]
 for n in(300,300):
  group=[[rg.getrandbits(23)<<41|d,rg.getrandbits(64),rg.getrandbits(64),rg.getrandbits(64)]for _ in range(n)];ab=len(o.msgs);recs.append(group);a5=Port2(J.mem)
  for c in range(0,n,K1):
   e0=o.extra;a4=group[c:c+K1]
   for ac in a4:o.put(d,ac)
   deliver_arr(E,rt[0],rt[1],J,[0,32*len(a4),1,0]+[w for r in a4 for w in r]+[STALE]*4*(K1-len(a4)),end,o.extra-e0,calls,cert,a5)
  a8=a5 .blocks==o.msgs[ab:]and not a5 .msgs and state_ok(J.mem,o);I=Port2(J.mem);run_k(E,rt[0],rt[1],"fin",dict(x0=MB2,x7=LC2+8*d,x3=(LC2+8*d)|RC2,x8=blk2(d)+32,x1=1,x2=J.E2,x6=end),J,calls,o.fin_count([d]),I);a9,aa=o.final(d);a8&=I.blocks==([a9]if a9 else[])and I.msgs==[[LC2+8*d,J.E2,aa,0,0,0]]
  for m in a5 .blocks+I.blocks+I.msgs:unpack(dm,m,wr if len(m)>6 else None)
  finals.append((aa,a8,a5 .blocks))
 f1,f2,a6=finals[0][0],finals[1][0],finals[1][2];q=300//K2B*K2B;a7=len(a6)>=1 and a6[0][1]==(base+32*q)&M64 and a6[0][4:4+4*(300-q)]==[w for r in recs[0][q:300]for w in r];ok=(all(a8 for _,a8,_ in finals)and rule_ok({d:dm},{d:wr},o,{d:f2})and a7 and(f1,f2,dm[J.E2])==((base+9600)&M64,(base+19200)&M64,f2)and o.cursor(d)==J.mem[(LC2+8*d)|RC2]==(base+600//K2B*32*K2B)&M64 and all([dm.get((base+32*t+8*w)&M64)for w in range(4)]==recs[0][t]for t in range(300))and all([dm.get((f1+32*t+8*w)&M64)for w in range(4)]==recs[1][t]for t in range(300)));return dict(ok=ok,first_final_rewritten_identically=a7)
def gate_checks(E,rt,cert):
 rprog,rlab=rt[0],rt[1];K=cert.K_GATE;a5,J=K&((1<<64)-1),K>>64;cases=[(0,True),(K-1,True),(K,False),(K+1,False),(((J-1)<<64)|((1<<64)-1),True),((J<<64)|((1<<64)-1),False),((1<<128)-1,False),];ok=True;I=set()
 for u,a4 in cases:o={};n=E.a64 .a64_run(rprog,rlab,"gate",dict(x0=u&((1<<64)-1),x1=u>>64,x2=a5,x3=J,x30=E.SNT),{},out_regs=o);I.add(n);ok&=(o["x0"]==(1 if a4 else 0))
 ok&=(I=={4}and 2*cert.COIN_SCALAR+4+8<=cert.GATE_INSNS);return dict(ok=ok,counts=sorted(I),cases=len(cases))
def phase0_xemit_hhclr(E,rt,rg,cert):
 rprog,rlab=rt[0],rt[1];g0,g1,I=Node(dd(int)),Node(dd(int)),Node(dd(int));CNM_BASE=90*(1<<40);YR_B=E.YR0
 for nd in(g0,g1,I):nd.mem[CNM_BASE]=rg.getrandbits(64);nd.mem[CNM_BASE+56]=rg.getrandbits(64)
 aa=rg.getrandbits(193)<<63;d1=(aa>>102)&M41;d2=(aa>>143)&M41;A0=aa+0x123456789ABCDEF;A1=A0+(1<<61);ae=A0^(1<<102);rl=Node(dd(int),sp=0x110000);hm=Node(dd(int),sp=LINK+0x110000);calls=[]
 for nd in(g0,g1,I):first=(LO1+32)if d1<H1 else(UP1+32);off=d1 if d1<H1 else d1-H1;run_k(E,rprog,rlab,"init3",dict(x0=1,x1=8*d1,x2=first+off*B1,x3=B1),nd,calls,cert.INIT3_BKT+1)
 a5=M41<<3;a8=777;af=778
 for u in range(64):Av=A0 if u==63 else((A0^((u+1)<<185))& ~(M41<<102))|(d1<<102);E.put(g0 .mem,YR_B+u*512+288,Av)
 E.put(g1 .mem,YR_B+288,A1);E.put(I.mem,YR_B+288,A0);E.put(I.mem,YR_B+512+288,ae);a6=d1^1;am=(LO1+32+a6*B1)if a6<H1 else(UP1+32+(a6-H1)*B1);run_k(E,rprog,rlab,"init3",dict(x0=1,x1=8*a6,x2=am,x3=B1),I,calls,cert.INIT3_BKT+1);p0=Port1(g0 .mem);aq=E.a64 .a64_run(rprog,rlab,"xemit",dict(x1=YR_B,x2=101,x26=3<<20,x3=64,x4=32,x6=a8<<49,x8=a5,sp=g0 .SP,x30=E.SNT),g0 .mem,dev=p0 .dev);p1=Port1(g1 .mem);n1=E.a64 .a64_run(rprog,rlab,"xemit",dict(x1=YR_B,x2=202,x26=9<<20,x3=1,x4=32,x6=a8<<49,x8=a5,sp=g1 .SP,x30=E.SNT),g1 .mem,dev=p1 .dev);p2=Port1(I.mem);b1=E.a64 .a64_run(rprog,rlab,"xemit",dict(x1=YR_B,x2=303,x26=0,x3=1,x4=32,x6=af<<49,x8=a5,sp=I.SP,x30=E.SNT),I.mem,dev=p2 .dev);b0=E.a64 .a64_run(rprog,rlab,"xemit",dict(x1=YR_B+512,x2=303,x26=4,x3=1,x4=32,x6=a8<<49,x8=a5,sp=I.SP,x30=E.SNT),I.mem,dev=p2 .dev);ok=(aq==64*cert.XEMIT_HOT+cert.FLUSH1_X+1 and n1==cert.XEMIT_HOT+1 and b1==cert.XEMIT_HOT+1 and b0==cert.XEMIT_HOT+1 and len(p0 .blocks)==1);run_k(E,rprog,rlab,"fin1",dict(x7=8*d1,x1=1),g0,calls,1+cert.FIN1_X,p0);run_k(E,rprog,rlab,"fin1",dict(x7=8*d1,x1=1),g1,calls,1+cert.FIN1_X,p1);run_k(E,rprog,rlab,"fin1",dict(x7=8*d1,x1=1),I,calls,1+cert.FIN1_X,p2);a9=Node.R2+4*32*K2B;ar=sorted({d2}|{(g0 .mem[(LO1 if d1<H1 else UP1)+(d1 if d1<H1 else d1-H1)*B1+32+32*u])&M41 for u in range(64)})
 for a4 in ar:rl.mem[LC2+8*a4]=blk2(a4)+32;rl.mem[(LC2+8*a4)|RC2]=Node.R2
 ah=Port2(rl.mem);a7=Sender(lambda d:Node.R2,a9,cert,Node.E2);b2=p0 .blocks+p1 .blocks+p2 .blocks
 for m in b2:
  e0=a7 .extra
  for t in range(m[1]//32):ap=list(m[4+4*t:8+4*t]);a7 .put(ap[0]&M41,ap)
  deliver_arr(E,rprog,rlab,rl,m,a9,a7 .extra-e0,calls,cert,ah,a7 .tag,3)
 ab=Port2(rl.mem)
 for a4 in ar:run_k(E,rprog,rlab,"fin",dict(x0=MB2,x7=LC2+8*a4,x3=(LC2+8*a4)|RC2,x8=blk2(a4)+32,x1=1,x2=Node.E2B,x11=a7 .tag,x6=a9),rl,calls,a7 .fin_count([a4]),ab)
 CC=(NEXT_BASE-asr1(Node.R2))%2**64;hm.mem[Node.E2]=Node.R2;home_state(hm,64,CC);hj=HomeJoin(cert,CC,64)
 for m in ah.blocks+ab.blocks+ab.msgs:
  if dest(m[0])!=d2:continue
  if len(m)==6:deliver(E,rprog,rlab,hm.mem,m,calls,cert,cnt=CNTE)
  else:
   f,ay=hj.block(m,hm.mem[Node.E2],2**cert.HB);ok&=ay;deliver(E,rprog,rlab,hm.mem,m,calls,cert,fused=True,r2=(Node.R2,a9),formula=f)
 aj=hm.mem[Node.E2B];mh=(aj-Node.R2)//32;ok&=(len(hj.cum)==1 and mh>=2);ak,ao=hj.cum[0];ok&=({hm.mem[ak+16],hm.mem[ao+16]}=={101,202});ai=E.a64 .a64_run(rprog,rlab,"hhclr",dict(x1=Node.R2,x2=aj,x10=HT,sp=hm.SP,x30=E.SNT),hm.mem);ok&=(ai==3+6*mh and all(hm.mem[HT+8*b]==0 for b in hj.heads));wm=dd(int);wk=(5,9,11);wm[M64-55],wm[M64-23],wm[8]=wk;wm.update({HT+8*k:7 for k in wk});ok&=E.a64 .a64_run(rprog,rlab,"hhclr",dict(x1=M64-63,x2=32,x10=HT,x30=E.SNT),wm)==21 and not any(wm[HT+8*k]for k in wk);ok&=all(c==f for _,c,f in calls);au,ag,J,B0=34*(1<<40)+0x20000,34*(1<<40)+0x40000,(0x987654321<<192)+0x77700000000,(0x987654321<<192)
 def wrep(m,i,a0,a1,c0,c1):
  o=au+128*i;E.put(m,o,a0);E.put(m,o+32,a1);E.put(m,o+64,c0,2);E.put(m,o+80,c1,2)
  for w in range(4):m[o+96+8*w]=64*i+w+1
 def rps(p_,l_,n,m):o={};s=E.a64 .a64_run(p_,l_,"psel",dict(x0=au,x1=n,x2=ag,sp=g0 .SP,x30=E.SNT),m,out_regs=o);return o["x0"],s
 ac=1<<80
 for i,rv in enumerate([(0,9,1,2),(J+50,J+10,1,2),(J+50,J+50,1,2),(J+50,J+50+(1<<63),1,2),(J+50,J+50+(1<<64),1,2),(J+50,J+90,ac|7,ac|7),(J+100,J+200,10,11),(J+300,J+400,11,12),(J+100,J+500,12,13),(J+600,J+200,12,13),(J+100+(1<<192),J+200+(1<<192),ac|10,ac|11),(J+200+(1<<192),J+250+(1<<192),14,15),(J+(2<<192),J+(2<<192)+((1<<63)-1),14,15),(B0+(3<<192)+(1<<63)-5,B0+(3<<192)+(1<<63)+5,16,17),(J+(4<<192),J+(4<<192)+9,18,19)]):wrep(g0 .mem,i,*rv)
 az,av=rps(rprog,rlab,15,g0 .mem);ax,_=rps(rprog,rlab,13,g0 .mem);z0,_=rps(rprog,rlab,0,g0 .mem);sa=[E.get(g0 .mem,ag+128*q+32*s)for q in range(4)for s in(0,1)];sc=[E.get(g0 .mem,ag+128*q+64+16*s,2)for q in range(4)for s in(0,1)];ad=Node(dd(int))
 for q in range(3):wrep(ad.mem,q,B0+((q+1)<<192)+10,B0+((q+1)<<192)+20,2*q+1,2*q+2)
 for i in range(3,1030):wrep(ad.mem,i,B0+(3<<192)+15,B0+(3<<192)+20,999+2*i,1000+2*i)
 aw,al=rps(rprog,rlab,1030,ad.mem);mp,ml=E.a64 .a64_load(open(os.path.join(E.HERE,"kernels","route.s")).read().replace("tbnz x12,#63,5f\n",""));gm=Node(dd(int))
 for q in range(4):wrep(gm.mem,q,J+(q<<192)+10,J+(q<<192)+10+(1<<63),2*q+1,2*q+2)
 ok&=(az==1 and ax==0 and z0==0 and aw==0 and rps(mp,ml,4,gm.mem)[0]==1 and rps(rprog,rlab,4,gm.mem)[0]==0 and al<=12+220*1024<=cert.PSEL_POOL and((sa[4]>>63)!=(sa[5]>>63))and((sa[6]>>63)!=(sa[7]>>63))and len(set(sa))==8 and sc==[10,11,ac|10,ac|11,14,15,16,17]and all(sa[2*q]>0 and 0<sa[2*q+1]-sa[2*q]<(1<<63)for q in range(4)));TRACE.extend([("G",g0 .mem),("G",g1 .mem),("G",I.mem),("G",ad.mem),("R",rl.mem),("H",hm.mem)]);return dict(ok=ok,xemit_64_insns=aq,hhclr_insns=ai,home_records=mh,candidates=len(hj.cum),psel_ok_steps=av,psel_wc_steps=al)
def wait_dispatch_checks(E,rt,rg,cert):
 rprog,rlab=rt[0],rt[1];nd=Node(dd(int));o=LINK_GEN+0x10000;a4=[1,1,1]
 def wake_cnt(mem,a5,ev):
  assert a4,"unexpected extra WFE in wcnt";mem[o]+=a4 .pop(0)
  if a5[0]==(o>>11):a5[0]=None;ev[0]=1
 a9=E.a64 .a64_run(rprog,rlab,"wcnt",dict(x0=o,x1=0,x30=E.SNT),nd.mem);a8=E.a64 .a64_run(rprog,rlab,"wcnt",dict(x0=o,x1=3,x30=E.SNT),nd.mem,wfe=wake_cnt);a4=[1];nd.mem[o]=0;ac=E.a64 .a64_run(rprog,rlab,"wcnt",dict(x0=o,x1=1,x30=E.SNT),nd.mem,wfe=wake_cnt,ev0=1);ok=(a9==4 and a8==4+3*cert.WAKE and ac==4+2*cert.WAKE and nd.mem[o]==1 and not a4);a6=(1<<45)+(257<<31);HDR=LINK_GEN+2048;CNTY=LINK_GEN+0x12000;I=Node(dd(int));I.mem[CNTY]=0;I.mem[CNTY+8]=2;o0={};E.a64 .a64_run(rprog,rlab,"ydlv0",dict(x15=a6+64,x29=CNTY,x30=E.SNT),I.mem,out_regs=o0);J=[[rg.getrandbits(64)for _ in range(56)]for _ in range(2)];pend=[(rlab["ydlv"][0],J[0]),(rlab["ydlve"][0],J[1])];py=PortY(I.mem);ad=([(LINK_GEN,LINK_GEN+512),(a6,a6+4096)],[(a6,a6+4096),(PORTY_BASE,PORTY_BASE+512)])
 def wake_ydsp(mem,a5,ev):
  assert pend,"unexpected extra WFE in ydsp";ab,pl=pend.pop(0)
  for w,v in enumerate(pl):mem[LINK_GEN+64+8*w]=v
  mem[HDR]=ab
 a7=E.a64 .a64_run(rprog,rlab,"ydsp",dict(o0,sp=I.SP,x30=E.SNT),I.mem,dev=py.dev,ls64=ad,wfe=wake_ydsp);aa=3+(cert.YDLV+8)+(cert.YDLVE+8);ok&=(a7==aa and not pend and 9<=cert.WY and I.mem[CNTY]==2 and I.mem[HDR]==0 and len(py.msgs)==1 and all([I.mem[a6+512*u+64+8*w]for w in range(56)]==J[u]for u in(0,1)));TRACE.extend([("G",nd.mem),("G",I.mem)]);return dict(ok=ok,wcnt_0=a9,wcnt_3=a8,ydsp_2=a7,wake=cert.WAKE,wy=cert.WY)

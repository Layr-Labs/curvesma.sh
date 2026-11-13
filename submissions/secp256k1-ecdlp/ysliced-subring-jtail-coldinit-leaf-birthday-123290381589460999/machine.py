import random;K1,K2B=64,256;B1,B2=64*K1,64*K2B;UP1=28672*2**40;UP=2**64-32768*2**40;LO1,H1=24576*2**40,1<<40;LC2,RC2=1<<45,1<<46;MB2=(1<<48)-128;PW=(1<<48)-168;CNT=(1<<55)-8192;CNTE=(1<<55)-6144;DONE=(1<<48)-4096;ST=(1<<55)-160;HT=2**64-(1<<55);NEXT_BASE=2**64-21846*2**40;rn,en,dc=range,enumerate,dict
def sg(v):return v-(1<<64)if(v&(1<<63))else(v&(2**64-1))
def asr1(v):return(sg(v)>>1)&(2**64-1)
LINK=(1<<55)-2**40;LINK_MSG=LINK+0x8000;WIN=(LINK,LINK+0x8000);LINK_GEN=34<<40;LGM=LINK_GEN+0x1000;PYB=LINK_GEN+0xB000;STL=0x5A5A5A5A5A5A5A5A;M41=2**41-1;DM=(2**64-1)^M41;HSR,DHR,_NID={},{},[0]
def pz(a):return(a*0x9E3779B97F4A7C15^STL)%2**64|1
class PZ(dict):
 # never-written words read nonzero (proof.md 5.1 cold start); zr: link MMIO (0, empty) and charged zeroed ranges
 zr=((LINK,LINK+65536),(LINK_GEN,LINK_GEN+65536))
 def __missing__(z,a):return 0 if any(lo<=a<hi for lo,hi in z.zr)else pz(a)
def dd(_=None):return PZ()
def z0(E,rt,m,*A):[E.arun(*rt,"tinit",dc(x0=1,x1=a,x2=0),m)for a in A]
def blk1(d):return LO1+B1*d if d<H1 else UP1+B1*(d-H1)
def blk2(d):return UP+B2*d
LC1_TABLES=[(0,8*H1,LO1+32,B1),(8*H1,16*H1,UP1+32,B1),(LO1,LO1+B1*H1,0,0),(UP1,UP1+B1*H1,0,0)]
def dest(word):return word>>3&M41
def dest1(word):return word>>3
class View:
 # a sender's block as the receiver's aperture shows it (read in place)
 def __init__(z,m,a,n,d=STL):z.m,z.a,z.n,z.d=m,a,n,d
 def __len__(z):return z.n
 def __getitem__(z,i):
  if isinstance(i,slice):return[z[k]for k in rn(*i.indices(z.n))]
  if not 0<=i<z.n:raise IndexError
  return z.m.get(z.a+8*i,z.d if i>3 and z.d is not None else pz(z.a+8*i))
class Port1:
 # proof.md 5.1: replay (copied at the completing store, receiver after sender); with rv, run in place
 def __init__(z,mem,rv=None):z.m,z.blocks,z.pend,z.rv=mem,[],[],rv;z.dev=(0,1<<64,z.store,z.read,lambda:not z.pend)
 def send(z,a,m,q):
  z.pend+=[a]
  if z.rv:z.rv(m)
  else:q+=[list(m)]
 def store(z,a,v):
  if a%B1==0 and(LO1<=a<LO1+B1*H1 or UP1<=a<UP1+B1*H1):z.send(a,View(z.m,a,4+4*K1),z.blocks)
 def read(z):z.pend=[]
class Pt2(Port1):
 def __init__(z,mem,mb=MB2,rv=None):Port1.__init__(z,mem,rv);z.msgs,z.mb=[],mb
 def store(z,a,v):
  mb=z.mb
  if a==mb+40:z.send(a,View(z.m,mb,6,None),z.msgs)
  elif not mb<=a<mb+48 and UP<=a<UP+B2*(1<<41)and(a-UP)%B2==0:z.send(a,View(z.m,a,4+4*K2B),z.blocks)
class LV:
 def __init__(z,mem,msg,writes,lbase=LINK):z.m,z.msg,z.writes,z.lbase=mem,msg,writes,lbase;z.base=lbase+32 if lbase==LINK and len(msg)>6 else lbase
 def __getitem__(z,a):return z.msg[(a-z.base)//8]if z.base<=a<z.base+8*len(z.msg)else z.m[a]
 def get(z,a,d=None):return z[a]if z.base<=a<z.base+8*len(z.msg)else z.m.get(a,d)
 def __setitem__(z,a,v):
  assert not z.lbase<=a<z.lbase+0x8000,"write into link aperture"
  if a not in(CNT,CNTE):z.writes+=[a]
  z.m[a]=v
DHP=tuple(LINK+64*k for k in rn(1,9));HSP=tuple(LINK+64*k for k in rn(1,7))
def dlvh0(E,o,B,mem,cnt):r0={};E.arun(o,B,"dlvh0",dc(x29=cnt),mem,out_regs=r0);assert tuple(r0[r]for r in("x0","x1","x7","x8","x24","x25","x26","x27"))==DHP and r0["x28"]==1;return r0
def dlr(E,I,B,mem,msg,cq,Q,writes=None,fused=False,r2=None,formula=None,cnt=CNT):
 e,o="dlv"if len(msg)==6 else"dlvh",{};lb=LINK_MSG if e=="dlv"else LINK;cnt in mem or z0(E,(I,B),mem,cnt);c0=mem[cnt];J=DHR.get(id(mem))or dlvh0(E,I,B,mem,cnt)if e=="dlvh"else dc(x0=lb,x28=1);regs=dc(J,x19=STL,x29=cnt);n=E.arun(I,B,e,regs,LV(mem,msg,[]if writes is None else writes,lb),out_regs=o,ls64=([WIN]+([r2]if r2 else[]),[r2]if r2 else[]))
 if e=="dlvh":DHR[id(mem)]=o;assert tuple(o[r]for r in("x0","x1","x7","x8","x24","x25","x26","x27"))==DHP
 ok=o["x28"]==1 and o["x29"]==cnt and mem[cnt]==c0+1;cq+=[(e,n if ok else-1,Q.DLV+Q.RCV if e=="dlv"else formula)]
def hst(o):return o["x14"]<<64|o["x13"],o["x5"]<<64|o["x15"]
def h2s0(E,I,B,J,end,o):r0={};E.arun(I,B,"h2s0",dc(x12=end,x0=o>>64,x8=o%2**64),J.m,out_regs=r0);assert(r0["x3"],r0["x10"],r0["x18"],r0["x12"],r0["x29"],r0["x6"],r0["x28"])==(PW,LC2,32*K2B,end,DONE,1<<62,0);assert tuple(r0[r]for r in("x0","x1","x2","x7","x8","x16"))==HSP;assert hst(r0)==(1<<126,o)and(J.m[PW-8]<<64|J.m[PW-16])==o;return r0
def dla(E,a4,B,I,msg,end,extra,cq,Q,a5,tag=0,a6=1<<66):L,o=msg[1]//32,{};r0=HSR.get((I.uid,end))or h2s0(E,a4,B,I,end,a6);I.m[PW]=tag;I.m.setdefault(DONE,0);d0,(D0,f0)=I.m[DONE],hst(r0);regs=dc(r0,x19=STL,sp=I.SP);n=E.arun(a4,B,"h2s",regs,LV(I.m,msg,[]),out_regs=o,dev=a5 and a5.dev,ls64=([WIN],[]));assert tuple(o[r]for r in("x0","x1","x2","x7","x8","x16"))==HSP;HSR[(I.uid,end)]=o;fin=L<K1;a7=D0+(msg[2]if fin else 0)-(1<<126 if fin and f0==1 else 0)-1;J=a7==0;assert I.m[DONE]==d0+J and hst(o)==((1<<126,a6)if J else(a7,f0-fin));cq+=[("h2s",n,Q.h2s_count(L)+extra+(Q.FIN_CNT+(f0==1)if fin else 0)+Q.RCV1+Q.H2S_DONE*J)];return J
class Sender1:
 def __init__(z,Q):z.cert,z.n,z.slots,z.cnt,z.msgs,z.seq,z.extra=Q,{},{},{},[],[],0
 def block(z,d,n,c):return[8*d,32*n,c,pz(blk1(d)+24)]+[w for r in z.slots.get(d,[None]*K1)for w in(r or[STL]*4)]
 def put(z,d,o):
  s,n=z.slots.setdefault(d,[None]*K1),z.n.get(d,0);s[n]=o;z.seq+=[(d,o)];z.n[d]=(n+1)%K1
  if n+1==K1:z.cnt[d]=z.cnt.get(d,0)+1;z.msgs+=[z.block(d,K1,z.cnt[d])];z.extra+=z.cert.FLUSH1_X
 def finals(z,o):
  A=[]
  for d in o:A+=[z.block(d,z.n.get(d,0),z.cnt.get(d,0)+1)];z.cnt[d],z.n[d]=0,0
  return A
def state1_ok(mem,o):
 for d,s in o.slots.items():
  b=blk1(d)
  if mem[8*d]!=b+32+32*o.n.get(d,0)or mem[b+16]!=o.cnt.get(d,0):return False
  if any([mem[b+32*t+32+8*w]for w in rn(4)]!=s[t]for t in rn(K1)if s[t]is not None):return False
 return True
class Snd:
 def __init__(z,rc0,end,Q,tag=0):z.rc0,z.end,z.cert,z.tag=rc0,end,Q,tag;z.n,z.rc,z.slots,z.msgs,z.seq,z.sent,z.drops,z.extra={},{},{},[],[],{},0,0
 def cursor(z,d):return z.rc.get(d,z.rc0(d))
 def put(z,d,o):
  s,n=z.slots.setdefault(d,[None]*K2B),z.n.get(d,0);s[n]=o;z.seq+=[(d,o)];z.n[d]=(n+1)%K2B
  if n+1<K2B:return
  c=z.cursor(d);z.rc[d]=(c+32*K2B)&M64;z.extra+=z.cert.FLUSH2_X if sg(c)<sg(z.end)else z.cert.FLUSH2_D
  if sg(c)<sg(z.end):z.msgs+=[[LC2+8*d,c,32*K2B,z.tag]+[w for r in s for w in r]];z.sent[d]=z.sent.get(d,0)+1
  else:z.drops+=1;z.rc[d]=z.end
 def final(z,d):n,c=z.n.get(d,0),z.cursor(d);msg=[LC2+8*d,c,32*n,z.tag]+[w for r in z.slots[d]for w in(r or[STL]*4)]if n and sg(c)<sg(z.end)else None;o=(c+32*n)&(2**64-1);return msg,(o if sg(o)<sg(z.end)else z.end)
 def fin_count(z,o):
  f=1
  for d in o:n=z.n.get(d,0);f+=z.cert.FINH_DEST+(0 if not n else z.cert.FINH_SENT if sg(z.cursor(d))<sg(z.end)else z.cert.FIN_HELD)
  return f
def kref(E,mem,I,q2,W2,KN,o,Q):
 meter,A=0,[]
 for l in rn(q2):
  e=I+E.ENTRY*l;w3=mem.get(I+E.ENTRY*(l-l%6)+288,0)
  for s in(0,1):
   J=E.fracs(mem,e,s);J[2],J[3]=E.ws(J[2],E.get(mem,e+256,2)),E.ws(J[3],E.get(mem,e+272,2));dev=0
   for t in E.KO:f=E.record_fields(E.leaf_w(J,t));meter+=Q.JX[t][f[0]];dev|=bool(f[0]);o.put(f[1],[f[3],f[4],W2+l%6,w3]);A+=[f]
   meter+=Q.SJOIN if dev else 0
 return meter,A
def ez(m,e,w=0):m.update({e+o:0 for o in(256,264,272,280)});m[e+288]=w
def plant(E,mem,sb,I,kn=0):
 for l,J in en(I):
  ez(mem,sb+E.ENTRY*l,kn+(l//6<<42))
  for s in(0,1):
   c=E.inverse_dft(J[3*s:3*s+3]);v=[f(E.t128(E.BT[m+k]*c[i]%E.p))for i,k,f in((0,0,int),(1,0,E.wneg),(2,-1,int))for m in(0,1)]
   for o,w in zip(E.FRAC_OFF[s],v):E.put(mem,sb+E.ENTRY*l+o,w,2)
def sok(mem,o):
 for d,s in o.slots.items():
  b,e=blk2(d),LC2+8*d
  if mem[e]!=b+32+32*o.n.get(d,0)or mem[e|RC2]!=o.cursor(d):return False
  if any([mem[b+32*t+32+8*w]for w in rn(4)]!=s[t]for t in rn(K2B)if s[t]is not None):return False
 return True
M64=2**64-1
def unpack(I,m,o=None):
 if len(m)==6:
  for w in rn(4):I[(m[1]+8*w)&M64]=m[2+w]
  return
 for t in rn(m[2]//32):
  for w in rn(4):I[(m[1]+32*t+8*w)&M64]=m[4+4*t+w]
  if o is not None:o+=[(m[1]+32*t)&M64]
def rlk(J,writes,o,finals):
 ok,cnt=True,{}
 for d,I in o.seq:
  t=cnt[d]=cnt.get(d,-1)+1;c=(o.rc0(d)+32*t)&M64
  if sg(c)<sg(o.end)and(finals is not None or t<K2B*o.sent.get(d,0)):ok&=[J[d].get((c+8*w)&M64)for w in rn(4)]==I
 if finals is not None:ok&=all(f==(((o.rc0(d)+32*(cnt.get(d,-1)+1))&M64)if sg((o.rc0(d)+32*(cnt.get(d,-1)+1))&M64)<sg(o.end)else o.end)for d,f in finals.items())
 return ok and all(sg(o.rc0(d))<=sg(a)<sg(o.end)for d,w in writes.items()for a in w)
TR=[]
def reg_trace(Q):
 o={k:[(lo,hi,b)for r,_,lo,hi,b in Q.REGISTRY if r==k]for k in"GRH"};n=bad=0
 for role,mem in TR:
  for a in mem:n+=1;hit=[b for lo,hi,b in o[role]if lo<=a<hi];bad+=len(hit)!=1 or hit[0]=="analytic"
 return dc(nodes=len(TR),addresses=n,outside_or_analytic=bad,ok=bad==0 and n>0)
class Node:
 HB,SB,SP=0x100000,(42<<40)+(1<<37),(34<<40)+0x7000000;R2,E2,E2B=0x3000000,LINK+0x1F000,LINK+0x1F800;CL=LINK+0x200000
 def __init__(z,mem=None,sp=None):_NID[0]+=1;z.uid,z.m,z.SP=_NID[0],({}if mem is None else mem),(z.SP if sp is None else sp)
def ND(**k):return Node(dd(int),**k)
def rk(E,J,B,entry,regs,nd,cq,formula,I=None,ls64=None):o={};regs=dc(regs,sp=nd.SP);cq+=[(entry,E.arun(J,B,entry,regs,nd.m,out_regs=o,dev=I and I.dev,ls64=ls64),formula)];return o
class HJ:
 def __init__(z,Q,CC,ncl):z.cert,z.CC,z.ncl=Q,CC,ncl;z.chains,z.heads,z.nextd,z.recs,z.cum,z.seen={},{},{},{},[],[];z.nls0,z.lists,z.paths,z.segs=0,[0,0],[0,0],{}
 def block(z,m,F,NB):
  Q,c,cnt=z.cert,m[1],m[2];lo=F if sg(F)>sg(c)else c;hi=(c+cnt)&M64;I,J=(sg(lo)-sg(c))//32,(sg(hi)-sg(lo))//32;a5=sg(F)<=sg(c)and cnt==32*K2B;z.paths[0 if a5 else 1]+=1;f=Q.dlvh_copy(I,J);ok=(f==Q.DLVHB if a5 else f<=Q.DLVHS);k=(I%2,J%2,I>0,J>0);z.segs[k]=z.segs.get(k,0)+1
  for t in rn(I,I+J):
   a,o=(c+32*t)&M64,m[4+4*t:8+4*t];a4=not(t%2==0 and t+2<=I+J);z.recs[a]=o;z.seen+=[(o[0],o[1],a)];b=o[1]&(NB-1);ch=z.chains.setdefault(b,[]);f+=Q.HCHG_NEW
   if ch:
    match=None
    for j,e in en(reversed(ch)):
     f+=Q.HCHG_VISIT
     if e[1]==o[1]:
      f+=Q.HCHG_VER
      if(e[0]&DM)==(o[0]&DM):match=(e,j);break
     if j==len(ch)-1:f+=Q.HCHG_END*(j%2);break
    if match is not None:
     e,j=match;f-=1
     if len(z.cum)<z.ncl:z.cum+=[(e[2],a)];f+=Q.HCHG_LS1-Q.HCHG_LSA*(not a4);z.lists[a4]+=1
     else:z.nls0+=1;f+=Q.HCHG_LS0
   z.nextd[a]=z.heads.get(b,0);z.heads[b]=(asr1(a)+z.CC)&(2**64-1);ch+=[(o[0],o[1],a)]
  return f+Q.RCV,ok
def hz(E,rt,m):E.arun(*rt,"zero16384",dc(x0=16384,x1=HT),m);m.zr+=((HT,HT+2**53),)
def hms(o,ncl,CC):DHR.pop(id(o.m),None);o.m[ST],o.m[ST+8],o.m[ST+16],o.m[ST+24]=o.CL,o.CL+16*ncl,CC,-2*CC&(2**64-1)
def part_b(E,Q):
 E.a64.DEV[:]=[(LINK,LINK+65536),(LINK_GEN,LINK_GEN+65536)];rg=random.Random(E.SEED+32452843);J=E.kload("route");lk=E.kload("leaf")+("keys",);H2=1<<41;ag,am,G,q2=(0,1,H1,H1+1),(0,1,H2-2,H2-1),4,18;cq=[];ae={g:ND()for g in rn(G)};ap={d1:Node(sp=0x110000)for d1 in ag};homes={(d1,d2):ND(sp=LINK+0x110000)for d1 in ag for d2 in am};pool=[E.planted_x(rg,rg.choice(ag),rg.choice(am))[1]for _ in rn(10)];batches,shared={},0
 for g in ae:
  ar=[]
  for _ in rn(6*q2):
   if rg.random()<0.15:ar+=[E.planted_x(rg,0,0,rg.choice(pool))[0]];shared+=1
   else:ar+=[E.planted_x(rg,rg.choice(ag),rg.choice(am))[0]]
  batches[g]=[ar[6*l:6*l+6]for l in rn(q2)]
 aa,ak,blocks=True,True,0;af={d1:[]for d1 in ag}
 for g,o in ae.items():
  sb=o.SB-E.ENTRY*q2
  for t0,n,first in((0,2,LO1+32),(8*H1,2,UP1+32)):rk(E,J[0],J[1],"init1",dc(x0=n,x1=t0,x2=first,x3=B1),o,cq,Q.INIT1_BKT*n+1)
  KN=(g>>64)|(5+g)<<26;plant(E,o.m,sb,batches[g],KN)
  a4=Sender1(Q);meter,fl=kref(E,o.m,sb,q2,8*g,KN,a4,Q);ak&=[f[1:]for f in fl]==[E.class_label(ao[3*s+t])for ao in batches[g]for s in(0,1)for t in E.KO];I=Port1(o.m);rk(E,lk[0],lk[1],"keys",dc(x2=8*g+6,x3=sb+E.KEYS0),o,cq,Q.KEY_FIX+Q.KEY_INC*q2+meter+a4.extra,I,ls64=([(sb,o.SB)],[]));ok=I.blocks==a4.msgs and state1_ok(o.m,a4);nb=len(I.blocks)
  for t0,n in((0,2),(8*H1,2)):rk(E,J[0],J[1],"fin1",dc(x7=t0,x1=n),o,cq,1+Q.FIN1_X*n,I)
  ai=a4.finals(ag);ok&=I.blocks[nb:]==ai and all(o.m[8*d]==blk1(d)+32 and o.m[blk1(d)+16]==0 for d in ag);aa&=ok;blocks+=len(I.blocks)
  for m in I.blocks:af[dest1(m[0])].append(m)
 a8=Node.R2+3*32*K2B;aq={h:[]for h in homes};a6={h:[]for h in homes};ac=True
 for d1,o in ap.items():
  for t0,n in((0,2),(H2-2,2)):rk(E,J[0],J[1],"init3",dc(x0=n,x1=LC2+8*t0,x2=blk2(t0)+32,x3=B2),o,cq,Q.INIT3_BKT*n+1);rk(E,J[0],J[1],"tinit",dc(x0=n,x1=(LC2+8*t0)|RC2,x2=Node.R2),o,cq,Q.TINIT_ENTRY*n+1)
  a4=Snd(lambda d:Node.R2,a8,Q,Node.E2);I=Pt2(o.m);al=af[d1];rg.shuffle(al);dn=[]
  for m in al:
   e0=a4.extra
   for t in rn(m[1]//32):a9=list(m[4+4*t:8+4*t]);a4.put(a9[0]&M41,a9);aq[(d1,a9[0]&M41)].append(tuple(a9))
   dn+=[dla(E,J[0],J[1],o,m,a8,a4.extra-e0,cq,Q,I,a4.tag,G)]
  ac&=(dn==[False]*(len(al)-1)+[True]and o.m[DONE]==1 and hst(HSR[(o.uid,a8)])==(1<<126,G)and sum(m[2]for m in al if m[1]<32*K1)==len(al));ok=I.blocks==a4.msgs and not I.msgs and sok(o.m,a4);a7=Pt2(o.m)
  for t0,n in((0,2),(H2-2,2)):rk(E,J[0],J[1],"fin",dc(x0=MB2,x7=LC2+8*t0,x3=(LC2+8*t0)|RC2,x8=blk2(t0)+32,x1=n,x2=Node.E2B,x11=a4.tag,x6=a8),o,cq,a4.fin_count(rn(t0,t0+n)),a7)
  ai=[a4.final(d)for d in am];ok&=(a7.blocks==[m for m,_ in ai if m]and a7.msgs==[[LC2+8*d,Node.E2B,c,0,0,0]for d,(_,c)in zip(am,ai)]);aa&=ok;blocks+=len(I.blocks)+len(a7.blocks)
  for m in I.blocks+a7.blocks+a7.msgs:a6[(d1,dest(m[0]))].append(m)
 ad,ah,ab,aj=True,True,True,0;CC=(NEXT_BASE-asr1(Node.R2))%2**64
 for h,o in homes.items():
  hz(E,J,o.m);o.m[Node.E2]=Node.R2;hms(o,64,CC);hj=HJ(Q,CC,64)
  for m in a6[h]:
   if len(m)==6:dlr(E,J[0],J[1],o.m,m,cq,Q,cnt=CNTE)
   else:
    f,ax=hj.block(m,o.m[Node.E2],2**Q.HB);ab&=ax;dlr(E,J[0],J[1],o.m,m,cq,Q,fused=True,r2=(Node.R2,a8),formula=f)
  F,F0=o.m[Node.E2B],Node.R2;start=F0-(F0-Node.R2)%(32*K2B);aw=1+-(-(F-start)//(32*K2B));ac&=o.m[CNTE]==1 and o.m[CNT]+1==len(a6[h])==aw;au=[tuple(o.m[x+8*k]for k in rn(4))for x in rn(Node.R2,F,32)];ad&=au==aq[h][:len(au)]and F==min(Node.R2+32*len(aq[h]),a8);ah&=all(a9[0]&M41==h[1]and a9[2]>>3 in ae and a9[2]%8<6 and(a9[3]>>26)%2**16==5+(a9[2]>>3)and a9[3]>>42<q2//6 for a9 in au);ab&=([(o.m[o.CL+16*i],o.m[o.CL+16*i+8])for i in rn(len(hj.cum))]==hj.cum and all(o.m[HT+8*b]==v for b,v in hj.heads.items()));aj+=len(hj.cum)
 groups={}
 for g in ae:
  for ao in batches[g]:
   for x in ao:groups.setdefault(E.exk(x)[0],set()).add(g)
 cnt,av={},{}
 for e,c,f in cq:cnt[e],av[e]=cnt.get(e,0)+c,av.get(e,0)+f
 a5=dc(leaves=6*q2*G,home_candidates=aj,counts=cnt,formulas_ok=all(c==f for _,c,f in cq)and ak and aa and ad and ah and ab and ac and aj>0);a5["block_paths"]=bpt(E,J,lk,rg,Q);a5["grouped_paths"]=gpp(E,J[0],J[1],rg,Q);a5["ls64"]=l6c(E,rg);a5["va55"]=va55_checks(E);a5["encoding"]=enc(E);a5["phase4_round"]=p4r(E,J,rg,Q);a5["y_ring"]=ycs(E,J,rg,Q);a5["gate"]=gtc(E,J,Q);a5["wait_dispatch"]=wait_dispatch_checks(E,J,rg,Q);a5["rendezvous"]=rendezvous_checks(E,J,lk,Q);a5["keys_edges"]=keys_edges(E,J,lk,Q);a5["phase0"]=p0x(E,J,rg,Q);TR.extend([("G",n.m)for n in ae.values()]+[("R",n.m)for n in ap.values()]+[("H",n.m)for n in homes.values()]);a5["registry_trace"]=reg_trace(Q);return a5
class PtY(Port1):
 def __init__(z,mem,base=PYB,rv=None):Port1.__init__(z,mem,rv);z.base,z.msgs=base,[]
 def store(z,a,v):
  if a==z.base:z.send(a,View(z.m,a,64,None),z.msgs)
def keys_edges(E,rt,lk,Q):
 rg=random.Random(E.SEED+104729);nd,q2,M,cl=ND(),12,E.M64,[];sb=nd.SB-E.ENTRY*q2
 for l in rn(q2):
  for s in(0,1):
   w=[rg.gb(128)for _ in rn(6)];a=rg.gb(62);f=[(a,a),(a,M-2*a),(a,-2*a&M),(a,-a&M),(M-2*a,a),(-2*a&M,a)][(2*l+s)%6]
   w[4]=(f[0]+(w[0]>>64)+(w[1]>>64)-(w[2]>>64)-(w[3]>>64))%2**64<<64|w[4]&M;w[5]=(f[1]-(w[0]>>64)+(w[2]>>64))%2**64<<64|w[5]&M
   for o,v in zip(E.FRAC_OFF[s],w):E.put(nd.m,sb+E.ENTRY*l+o,v,2)
  ez(nd.m,sb+E.ENTRY*l)
 aa=Sender1(Q);meter,fl=kref(E,nd.m,sb,q2,7,7<<26,aa,Q)
 for d in{x[1]for x in fl}:rk(E,*rt,"init1",dc(x0=1,x1=8*d,x2=blk1(d)+32,x3=B1),nd,cl,Q.INIT1_BKT+1)
 ac=Port1(nd.m);rk(E,lk[0],lk[1],"keys",dc(x2=13,x3=sb+E.KEYS0),nd,cl,Q.KEY_FIX+Q.KEY_INC*q2+meter+aa.extra,ac,ls64=([(sb,nd.SB)],[]))
 j=[fl[6*l+3*s+1][0]for l in rn(q2)for s in(0,1)];return dc(ok=all(c==f for _,c,f in cl)and ac.blocks==aa.msgs==[]and state1_ok(nd.m,aa)and j[:6]==[1,0,2,0,1,2],leaf0_j=j[:6])
def rendezvous_checks(E,rt,lk,Q):
 # proof.md 5.1: replayed vs in place (handlers nested at the sender's completing store): equal memories, counts
 rg=random.Random(E.SEED+7919);q2=6*-(-(K2B+K2B//4+8)//36);xs=[[E.planted_x(rg,0,d)[0]for _ in rn(24)]for d in(0,1)];ls=[xs[rg.random()<.06][rg.randrange(24)]for _ in rn(6*q2)];a7=Node.R2+3*32*K2B;CC=(NEXT_BASE-asr1(Node.R2))%2**64;yb=(1<<45)+(257<<31);pl=[rg.gb(64)for _ in rn(168)];YL=([(LINK_GEN,LINK_GEN+512),(yb,yb+4096)],[(yb,yb+4096),(PYB,PYB+512)])
 def run(live):
  cl,g,rl,hs,ys=[],ND(),ND(sp=0x110000),[ND(sp=LINK+0x110000)for _ in"01"],[ND()for _ in"0123"]
  def home(m):
   h=hs[dest(m[0])]
   if len(m)==6:dlr(E,*rt,h.m,m,cl,Q,cnt=CNTE)
   else:dlr(E,*rt,h.m,m,cl,Q,fused=True,r2=(Node.R2,a7),formula=0)
  def yh(k,m):
   D[0]+=1;D[1]=max(D);e="ydlv"if k<3 else"ydlve";p=PtY(ys[k].m,rv=live and k<3 and(lambda v:yh(k+1,v)));o={};cl.append((e,E.arun(*rt,e,dc(yr[k],x19=STL),LV(ys[k].m,m,[],LINK_GEN),out_regs=o,dev=p.dev,ls64=YL),0));yr[k]=o;D[0]-=1
   for v in p.msgs:yh(k+1,v)
  p2=Pt2(rl.m,rv=live and home);relay=lambda m:dla(E,*rt,rl,m,a7,0,cl,Q,p2,Node.E2,1);p1=Port1(g.m,rv=live and relay);sb=g.SB-E.ENTRY*q2;plant(E,g.m,sb,[ls[6*l:6*l+6]for l in rn(q2)])
  for h in hs:hz(E,rt,h.m);h.m[Node.E2]=Node.R2;hms(h,64,CC)
  rk(E,*rt,"init1",dc(x0=1,x1=0,x2=LO1+32,x3=B1),g,cl,0);rk(E,*rt,"init3",dc(x0=2,x1=LC2,x2=blk2(0)+32,x3=B2),rl,cl,0);rk(E,*rt,"tinit",dc(x0=2,x1=LC2|RC2,x2=Node.R2),rl,cl,0)
  rk(E,lk[0],lk[1],"keys",dc(x2=13,x3=sb+E.KEYS0),g,cl,0,p1,ls64=([(sb,g.SB)],[]));rk(E,*rt,"fin1",dc(x7=0,x1=1),g,cl,0,p1)
  for m in p1.blocks:relay(m)
  p4=Pt2(rl.m,rv=live and home);rk(E,*rt,"fin",dc(x0=MB2,x7=LC2,x3=LC2|RC2,x8=blk2(0)+32,x1=2,x2=Node.E2B,x11=Node.E2,x6=a7),rl,cl,0,p4)
  for m in p2.blocks+p4.blocks+p4.msgs:home(m)
  for h in hs:DHR.pop(id(h.m),None)
  for i,v in en(pl):ys[0].m[yb+64+8*i+64*(i//56)]=v
  [z0(E,rt,n.m,LINK_GEN+0x12000)for n in ys];yr=[rk(E,*rt,"ydlv0",dc(x15=yb+64,x29=LINK_GEN+0x12000),n,cl,0)for n in ys];yr[0].update(x24=0,x25=1<<42,x26=0);D=[0,0]
  for u in rn(3):
   ys[0].m[yb+512*u+56]=u;p=PtY(ys[0].m,rv=live and(lambda v:yh(1,v)));yr[0]=rk(E,*rt,"ysnd",yr[0],ys[0],cl,0,p,YL);D[1]=max(D[1],(p.msgs+[0])[0]and 1)
   for v in p.msgs:yh(1,v)
  return[{a:v for a,v in n.m.items()if v}for n in[g,rl]+hs+ys],sorted((e,c)for e,c,_ in cl),len(p1.blocks),len(p2.blocks),len(p.msgs),D[1],[n.m[LINK_GEN+0x12000]for n in ys]
 a,b=run(0),run(1);return dc(ok=a[:2]==b[:2]and a[2]>1 and a[3]>0 and a[4]==1 and b[2:5]==(0,0,0)and(a[5],b[5])==(1,3)and all(c>0 for _,c in a[1])and[e for e,_ in a[1]].count("ydlve")==3 and a[6]==b[6]==[0,3,3,3],hop1_blocks=a[2],hop2_blocks=a[3],calls=len(a[1]),y_depth=b[5])
def ycs(E,rt,rg,Q):
 rb,rq=rt;A,CY,t=(1<<45)+(257<<31),LINK_GEN+0x12000,rg.rr(1<<15);sl=[rg.rr(1<<41)for _ in"01"];rs={s:[ND()for _ in"012"]for s in sl};pr=ND();[z0(E,rt,n.m,CY)for s in sl for n in rs[s]];Y=[];YL=([(LINK_GEN,LINK_GEN+512),(A,A+4096)],[(A,A+4096),(PYB,PYB+512)]);J=1;R={};ar=lambda e,r,n,o,p=None:E.arun(rb,rq,e,r,n.m,out_regs=o,dev=p and p.dev,ls64=YL)
 for u in rn(6):
  s=sl[u%2];Y+=[(s,[rg.gb(64)for _ in rn(56)])];pr.m[A+512*u+56]=s
  for w,v in en(Y[u][1]):pr.m[A+512*u+64+8*w]=v
 o={};ar("ydlv0",dc(x15=A+64,x29=CY),pr,o);P0=rg.gb(41)<<1;I=dc(o,x24=P0,x25=1<<42,x26=t<<10)
 for s in sl:
  for i,n in en(rs[s]):n.m[CY+16],n.m[CY+24]=A+64,A+64+512*(1 if(s,i)==(sl[0],1)else 8);R[s,i]={};ar("ydlv0",dc(x29=CY),n,R[s,i])
 for u,(s,_)in en(Y):
  p=PtY(pr.m);o={};J&=ar("ysnd",I,pr,o,p)==Q.YSND;I=o;m=p.msgs[0];J&=len(p.msgs)==1 and(m[2],m[3])==((s<<25|t<<10)%2**64,s>>39)
  for i,n in en(rs[s]):
   p=PtY(n.m);o={}
   for w,v in en(m[8:]):n.m[LINK_GEN+64+8*w]=v
   n.m[LINK_GEN+2048],n.m[CY+8]=rq["ydlve"if i>1 else"ydlv"][0],n.m[CY]+1;J&=ar("ydsp",dc(R[s,i],sp=n.SP),n,o,p)==(Q.YDLVE if i>1 else Q.YDLV)+15+4*(u<2)and len(p.msgs)==(i<2);R[s,i]=o;m=p.msgs and p.msgs[0]
 W=lambda u:Y[u][1][:37]+[P0+(u<<42)]+Y[u][1][38:]
 for s in sl:
  us=[u for u in rn(6)if Y[u][0]==s]
  for i,n in en(rs[s]):d=(s,i)==(sl[0],1);J&=[[n.m[A+512*j+64+8*w]for w in rn(56)]for j in rn(3)]==([W(us[0]),W(us[2]),[pz(A+1088+8*w)for w in rn(56)]]if d else[W(u)for u in us])and n.m[CY]==3 and n.m[CY+16]==0 and R[s,i]["x15"]==A+64+512*(3-d)
 RB,n,G=(1<<48)+0x10000,ND(sp=0x110000),[]
 for i in rn(3):h,t2,m,j,u,k=rg.gb(41),rg.rr(1<<15),rg.rr(1024),rg.rr(6),rg.rr(1<<21),rg.gb(37);w2,w3=j|m<<3|t2<<13|k%2**36<<28,k>>36|h<<1|u<<42;n.m[RB+64*i+8],n.m[RB+64*i+16]=w2,w3;G+=[(1024*(h<<15|t2)+m,w2,w3,RB+64*i)]
 p=Pt2(n.m);fq=E.arun(rb,rq,"fwdq1",dc(x0=MB2,x1=RB,x6=RB+192,x7=5,sp=n.SP,x19=STL),n.m,dev=p.dev)==45 and p.msgs==[[g%2**64,g>>64,w3,a<<8,5,w2]for g,w2,w3,a in G]
 a9=E.ksrc("route").replace("ysnd:\n","ysnd_orig:\n")+"\nysnd_nodmb:\n\tstlr x28,[x14]\n\tret\n";bp,bl=E.a64.a64_load(a9);pb=PtY(dd(int))
 try:
  E.arun(bp,bl,"ysnd_nodmb",dc(x14=PYB,x28=1),pb.m,dev=pb.dev);dmb_rejected=False
 except AssertionError:dmb_rejected=True
 return dc(ok=bool(J)and fq and dmb_rejected,counts=[Q.YSND,Q.YDLV,Q.YDLVE],fwdq1=fq,dmb_rejected=dmb_rejected)
def p4r(E,rt,rg,c):
 C,ok,N,BR,BH,a7=[],True,Node,(1<<48)+0x10000,LINK+0x20000,LINK+0x18000;d1,d2,gs=rg.rr(1<<41),rg.rr(1<<41),[rg.rr(3<<64,1<<66),rg.rr(1<<64)];h,J=N(dd(int),sp=LINK+0x110000),N(dd(int),sp=0x110000);[z0(E,rt,n.m,CNT,CNT+2048)for n in(h,J)]
 def run(I,e,regs,f,base=None,msg=None,lb=LINK_MSG):p=base and Pt2(I.m,base);C.append((E.arun(*rt,e,dc(regs,sp=I.SP,x19=STL,x28=1,x29=CNT),I.m if msg is None else LV(I.m,msg,[],lb),dev=p and p.dev),f));return p and p.msgs
 def q(I,m,cap=c.RELAY_CAP,lb=LINK_MSG+1024,bb=BR,ct=CNT):a5,aa,oc=I.m[ct],I.m[ct+2048],I.m[bb-8];C.append((E.arun(*rt,"dlvq",dc(x0=lb,x9=64,x10=bb-8,x11=bb+64*cap,sp=I.SP,x19=STL,x28=1,x29=ct),LV(I.m,m,[],lb)),(c.DLVQ if oc<bb+64*cap else c.DLVQ_DROP)+(c.DLVQ_CARRY if a5==2**64-1 else 0)+c.RCV));return(I.m[ct+2048]<<64|I.m[ct])==((aa<<64|a5)+1)
 R=[[rg.rr(1<<64),0,g%2**64,g>>64|rg.rr(1<<62)<<2]for g in gs]
 for i in rn(8):h.m[N.R2+8*i]=R[i//4][i%4]
 h.m[N.CL],h.m[N.CL+8],J.m[BR-8],hd=N.R2,N.R2+32,BR,2**45+8*d2;m1=run(h,"reqh",dc(x0=a7,x12=N.CL,x9=1,x10=hd,x13=BH),c.REQH_C+1,a7);ok&=m1==[[0,0,r[2],r[3],hd,BH+256*i]for i,r in en(R)]and all(q(J,m)for m in m1);m2=run(J,"fwdq",dc(x0=MB2+64,x1=BR,x6=BR+128,x7=8*d1),2*c.FWDQ+3,MB2+64);a4=[[rg.rr(1<<64)for _ in rn(32)]for _ in gs];BG=96*2**40;a8=LINK_GEN+0x18000;a9=LINK_GEN+0x10000
 for i,m in en(m2):
  g=N(dd(int));z0(E,rt,g.m,a9,a9+2048);g.m[BG-8]=BG;TR.append(("G",g.m));ok&=(m[1]&3)<<64|m[0]==gs[i]and q(g,m,cap=c.REQ_CAP,lb=LGM+1024,bb=BG,ct=a9)and g.m[a9]==1 and[g.m[BG+8*k]for k in rn(4)]==[gs[i]>>64,R[i][3],BR+64*i<<8,8*d1]
  for k in rn(32):g.m[BG+4096+8*k]=a4[i][k]
  m3=run(g,"sendr",dc(x0=a8,x1=BG+4096,x9=8,x7=8*d1,x8=g.m[BG+16]),8*c.SENDR+1,a8);ok&=len(m3)==8 and all(q(J,m)for m in m3)
 m4=run(J,"fwdr",dc(x0=MB2,x1=BR+128,x6=BR+64*18),16*c.FWDR+3,MB2)
 for m in m4:run(h,"dlv",dc(x0=LINK_MSG),c.DLV+c.RCV,msg=m)
 J.m[BR-8]=BR;d=N(dd(int),sp=0x110000);z0(E,rt,d.m,CNT,CNT+2048);d.m[BR-8]=BR+64*c.RELAY_CAP;w=N(dd(int),sp=0x110000);z0(E,rt,w.m,CNT,CNT+2048);w.m[BR-8]=BR+64*(c.RELAY_CAP-1);w.m[CNT]=2**64-2;cw=N(dd(int),sp=0x110000);z0(E,rt,cw.m,CNT,CNT+2048);cw.m[BR-8]=BR+64*c.RELAY_CAP;cw.m[CNT]=2**64-1;a6=[];E.arun(*rt,"dlvq",dc(x0=LINK_MSG+1024,x9=64,x10=BR-8,x11=BR+64*c.RELAY_CAP,sp=cw.SP,x19=STL,x28=1,x29=CNT),LV(cw.m,[5]*6,[],LINK_MSG+1024),dev=(CNT,4096,lambda o,v:a6.append((cw.m[CNT+2048],cw.m[CNT])),lambda:None,lambda:True));ok&=([h.m[BH+8*k]for k in rn(64)]==a4[0]+a4[1]and len(m4)==16 and q(d,[1]*6)and q(d,[1]*6)and d.m[BR-8]==BR+64*c.RELAY_CAP and BR not in d.m and q(w,[2]*6)and q(w,[3]*6)and q(w,[4]*6)and w.m[BR-8]==BR+64*c.RELAY_CAP and w.m[CNT]==1 and w.m[CNT+2048]==1 and 0 not in w.m and a6==[(0,0),(1,0)]and all(x==y for x,y in C)and max(c.PH4_HOPS.values())<=c.MSG and(h.m[CNT],J.m[CNT],d.m[CNT])==(16,18,2));TR.extend([("H",h.m),("R",J.m),("R",d.m),("R",w.m),("R",cw.m)]);return dc(ok=ok,calls=len(C),hops=c.PH4_HOPS)
def enc(E):
 a64=E.a64;good="ldp x1,x2,[x3,#-512]|ldp x1,x2,[x3,#504]|stp x1,x2,[x3],#504|stp x1,x2,[x3,#-512]!|ldr x1,[x3,#32760]|ldr x1,[x3,#-256]|ldr x1,[x3,#255]|str x1,[x3],#-256|add x1,x1,#4095|add x1,x1,#0x1000|tst x3,#0x3fffe0000|mov sp,x1|add sp,sp,#512|cset w9,hs|csel x1,x2,x3,lo|mov x1,#0xfff8000000000000|movk x4,#1,lsl #32|stlur x1,[x3,#-256]|ldr x22,[x3,#-152]|and x1,x2,#0x5555555555555555|ldaxr x2,[x0]|blr x30|wfe|cbhi x1,x2,f|cbhs xzr,x30,f".split("|");bad="stp x25,x26,[x3],#768|ldp x21,x22,[x3,#-768]|ldp x1,x2,[x3,#508]|ldp x1,x2,[x3,#-520]|ldp x1,x2,[x3,#12]|ldr x22,[x3,#-664]|ldr x1,[x3,#32768]|ldr x1,[x3,#260]|str x1,[x3],#256|ldr x1,[x3,#-257]!|add x1,x1,#4097|tst x3,#0x1c001|and x1,x2,#0|mov x1,#0x12345|movk x4,#1,lsl #8|stlur x1,[x3,#256]|ldp x1,x1,[x3]|ldp x3,x2,[x3],#16|ldadd x1,x2,[x3,#8]|ldaxr sp,[x0]|ldaxr x1,[x0,#8]|blr sp|wfe x0|ubfx x1,x2,#60,#8|bfi x1,x2,#40,#25|ldr x1,[x3,#-264]|add x1,xzr,#5|tst sp,#0x1c000|ldp xzr,xzr,[x3]|ldr sp,[x3]|mov xzr,#0x3ffe0|movk sp,#1,lsl #48|tbz x3,#64,1f|frob x1,x2|mul x1,x2,#1|adc x1,x2,#3|csel x1,x2,#3,eq|madd x1,x2,x3,#4|cset x1,x2|mvn x1,#2|mul w1,x2,x3|dmb sy|umulh w1,w2,w3|mvn x1,x2,lsl #64|movk x1,#1,lsr #16|cbhi sp,x2,f|cbhs x1,#2,f|cbhi x1,x2|stadd x1,[x2,#0]|ld64b x0,[x2,#0]|ldp x1,x2,[x3]!|ld64b x21,[x1]|ld64b x24,[x1]|add x1,x2,#1_0".split("|")
 def verdict(o):
  try:a64.a64_load("f:\n\t"+o+"\n\tret\n");return True
  except AssertionError as ex:assert"not encodable"in str(ex),ex;return False
 def reach(n):
  try:a64.a64_load("f:\n\tcbhi x1,x2,g\n"+"\tadd x1,x1,#1\n"*n+"g:\n\tret\n");return 1
  except AssertionError as ex:assert"branch reach"in str(ex);return 0
 A=dc(good_accepted=sum(map(verdict,good)),bad_rejected=sum(not verdict(i)for i in bad),good=len(good),bad=len(bad),cb_reach_1020_1024=(reach(254),reach(255)));A["ok"]=A["good_accepted"]==len(good)and A["bad_rejected"]==len(bad)and A["cb_reach_1020_1024"]==(1,0);return A
def va55_checks(E):
 a64=E.a64;I,J=a64.a64_load("ld:\n\tldr x2,[x1]\n\tret\n");A={};U=2**64-2**55
 for a in(0,2**55-8,U,2**64-8,2**55-4,2**55,2**56-8,2**60,U-8,2**64-2**56,2**63,2**64-4):
  try:E.arun(I,J,"ld",dc(x1=a),dd(int));A[hex(a)]="ran"
  except AssertionError as ex:A[hex(a)]="rejected"if"55-bit"in str(ex)else"other: "+str(ex)
 v=list(A.values());A["ok"]=v==["ran"]*4+["rejected"]*8;return A
def lsv(ex):return"LS64 rejected"if"LS64"in str(ex)else"other: "+str(ex)
def l6c(E,rg):
 a64=E.a64;a4=E.ksrc("route");I,J=a64.a64_load("ld:\n\tld64b x20,[x1]\n\tret\nst:\n\tst64b x16,[x2]\n\tret\novl:\n\tld64b x12,[x15]\n\tret\n");R2=(Node.R2,Node.R2+4*32*K2B);msg=[rg.gb(64)for _ in rn(4+4*K2B)]
 def run(a5,B,entry,regs,mem,lsr=([WIN,R2],[R2])):
  try:E.arun(a5,B,entry,dc(regs),mem,ls64=lsr);return"ran"
  except AssertionError as ex:return lsv(ex)
 A,o={},{};E.arun(I,J,"ld",dc(x1=LINK+64),LV({},msg,[]),out_regs=o,ls64=([WIN],[]));A["ld_win01"]=[o["x%d"%(20+k)]for k in rn(8)]==msg[4:12];mem={};A["st64b_r2_aligned"]=run(I,J,"st",dc(x2=R2[0]),mem);A["st64b_r2_written"]=all(mem.get(R2[0]+8*k)==0 for k in rn(8));A["ld64b_odd_slot"]=run(I,J,"ld",dc(x1=LINK+32),LV({},msg,[]));A["ld_out"]=run(I,J,"ld",dc(x1=0x5000040),{});A["st64b_into_window"]=run(I,J,"st",dc(x2=LINK+64),LV({},msg,[]));A["st_r2_mis"]=run(I,J,"st",dc(x2=R2[0]+32),{})
 def load(text):
  try:a64.a64_load(text);return"loaded"
  except AssertionError as ex:return lsv(ex)
 A["ld64b_x24"]=load("hi:\n\tld64b x24,[x1]\n\tret\n");A["ld64b_x21"]=load("odd:\n\tld64b x21,[x1]\n\tret\n");A["ld_base_ov"]=run(I,J,"ovl",dc(x15=LINK+64),LV({},msg,[]))
 def fused(text,c,r2):a5,B=a64.a64_load(text);mem=dd();hz(E,(a5,B),mem);z0(E,(a5,B),mem,CNT);CC=(NEXT_BASE-asr1(r2[0]))%2**64;mem[ST],mem[ST+8],mem[ST+16],mem[ST+24]=Node.CL,Node.CL+1024,CC,-2*CC%2**64;a8=LINK+0x300000;mem[a8]=c;m=[0,c,32*K2B,a8]+msg[4:];dr=dlvh0(E,a5,B,mem,CNT);r=run(a5,B,"dlvh",dc(dr,x29=CNT),LV(mem,m,[]),([WIN,r2],[r2]));return r if r!="ran"else("ran, R2 exact"if all(mem[c+8*w]==m[4+w]for w in rn(4*K2B))else"ran, R2 wrong")
 A["dlvh_fresh_block"]=fused(a4,R2[0],R2);a6=(LINK-32*K2B,LINK);A["r2_top"]=fused(a4,a6[0],a6);A["r2_link_rej"]=fused(a4,LINK,(LINK,LINK+512));a7=(R2[0]+32,R2[1]+32);A["r2_32mod64"]=fused(a4,a7[0],a7);ab=a4.replace("6:\nadd x13,x14,x2,asr #1\nDLVHPR x0\n","6:\nadd x13,x14,x2,asr #1\nsub x12,x0,#32\nDLVHPR x12\n");A["dlvh_odd"]=fused(ab,R2[0],R2)if ab!=a4 else"mutation not applied";hm=a4.replace("mov x1,x0\n","sub x1,x0,#32\n");hp,hl=a64.a64_load(hm);regs=dc(x0=HSP[0],x1=HSP[1],x2=HSP[2],x7=HSP[3],x8=HSP[4],x16=HSP[5],x10=LC2,x18=32*K2B,x12=1<<40,x3=PW,x29=DONE,x6=1<<62,x14=1<<62,x5=16);A["h2s_odd"]=run(hp,hl,"h2s",regs,LV(PZ({PW:0}),[0,32*2,1,0]+msg[4:4+4*K1],[]),([WIN],[]))if hm!=a4 else"mutation not applied";aa=("ld64b_odd_slot","ld_out","st64b_into_window","st_r2_mis","ld64b_x24","ld64b_x21","ld_base_ov","r2_32mod64","dlvh_odd","h2s_odd");a9=(A["ld_win01"]and A["st64b_r2_aligned"]=="ran"and A["st64b_r2_written"]and A["dlvh_fresh_block"]=="ran, R2 exact"and A["r2_top"]=="ran, R2 exact"and"write into link aperture"in A["r2_link_rej"]and all(A[k]=="LS64 rejected"for k in aa));return dc(ok=a9,rejected=len(aa))
def gpp(E,a6,B,rg,Q):
 GG,aj,I=3,6,2;NB=1<<aj;a5=32*K2B*128;ak=((300,9),(520,262),(250,7))
 def rcs(ac,ncl=64):
  cq,o=[],ND();ai,ae=LINK+0x400000,LINK+0x500000;a7=(2**64-8192)if ac else o.R2;E2P=(o.E2,o.E2+0x800000);pf=lambda g,s:ai+8*s if g==0 else E2P[(g-1)%2]+32*s;ah={}
  for w in rn(Q.ZW):o.m[HT+8*w]=0x5A5A5A5A5A5A5A5A
  rk(E,a6,B,"zero16384",dc(x0=Q.ZW,x1=HT),o,cq,Q.ZERO_PASS+1);ok=all(o.m[HT+8*w]==0 for w in rn(Q.ZW));rk(E,a6,B,"init3",dc(x0=I,x1=ai,x2=a7,x3=a5),o,cq,Q.INIT3_BKT*I+1);CC=(NEXT_BASE-asr1(a7))%2**64;hms(o,ncl,CC);z0(E,(a6,B),o.m,CNT);rk(E,a6,B,"dlvh0",dc(x29=CNT),o,cq,Q.DLVH0);ok&=all(o.m[ai+8*s]==((a7+a5*s)&M64)for s in rn(I));J=HJ(Q,CC,ncl);al,stopped=[],False;rc=[(a7+a5*s)&M64 for s in rn(I)];F=list(rc);ag,nb=[[None]*K2B for _ in rn(I)],[0]*I
  for g in rn(GG):
   if stopped:break
   al+=[g];ah[g]=[(rg.gb(64),rg.gb(64))for _ in rn(3)]
   for k,(px,py)in en(ah[g]):o.m[ae+(g*3+k)*16]=px;o.m[ae+(g*3+k)*16+8]=py
   streams=[]
   for s in rn(I):
    am=ak[g][s]+(0 if ac or s else 1);recs=[]
    for _ in rn(am):recs+=[[rg.gb(64),rg.gb(7)<<57|rg.gb(aj),rg.gb(64),rg.gb(64)]]
    if J.seen:aq,ap=J.seen[rg.rr(len(J.seen))][:2];recs[0][1]=ap;recs[0][0]=(aq&DM)|rg.gb(41)
    recs[1][1]=recs[0][1];recs[1][0]=(recs[0][0]&DM)|rg.gb(41);msgs=[]
    for ao in recs:
     ag[s][nb[s]]=ao;nb[s]+=1
     if nb[s]==K2B:msgs+=[[0,rc[s],32*K2B,pf(g,s)]+[w for r in ag[s]for w in r]];rc[s]=(rc[s]+32*K2B)&M64;nb[s]=0
    if nb[s]:msgs+=[[0,rc[s],32*nb[s],pf(g,s)]+[w for r in ag[s]for w in r or[STL]*4]]
    streams+=[msgs]
   order=[m for st in streams for m in st];rg.shuffle(order)
   for m in order:s=[t for t in rn(I)if pf(g,t)==m[3]][0];f,ar=J.block(m,F[s],NB);ok&=ar;ab=[(2**64-8192,1<<64),(0,1<<51)]if ac else[(o.R2,o.R2+I*a5)];n=E.arun(a6,B,"dlvh",dc(DHR.get(id(o.m))or dlvh0(E,a6,B,o.m,CNT),x29=CNT),LV(o.m,m,[]),out_regs=(au:={}),ls64=([WIN]+ab,ab));DHR[id(o.m)]=au;cq+=[("dlvh",n,f)];ok&=all(o.m[(m[1]+8*w)&M64]==m[4+w]for w in rn(4*K2B))
   ok&=all(o.m[HT+8*b]==J.heads.get(b,0)for b in rn(NB));ok&=all(o.m[(asr1(a)+CC)&M64]==J.nextd[a]and o.m[((asr1(a)+CC)&M64)+8]==J.recs[a][1]for a in J.nextd);ok&=[(o.m[o.CL+16*i],o.m[o.CL+16*i+8])for i in rn(len(J.cum))]==J.cum;ok&=o.m[ST]==o.CL+16*len(J.cum)and len(J.cum)<=ncl
   for s in rn(I):F[s]=(rc[s]+32*nb[s])&M64;dlr(E,a6,B,o.m,[0,E2P[g%2]+32*s,F[s],0,0,0],cq,Q)
   ok&=all(o.m[E2P[g%2]+32*s]==F[s]for s in rn(I));a4=[((a7+a5*s+32*k)&M64)for s in rn(I)for k in rn((sg(F[s])-sg((a7+a5*s)&M64))//32)];ok&=sorted(J.nextd)==sorted(a4);aa=[]
   for b in rn(NB):
    u=J.heads.get(b,0)
    while u:a=(2*((u-CC)%2**64))&M64;aa+=[a];u=J.nextd[a]
   ok&=sorted(aa)==sorted(a4)and all([o.m[(a+8*w)&M64]for w in rn(4)]==J.recs[a]for a in a4)
   if ac and g==1:stopped=any(o.m[a+8]==o.m[b+8]and o.m[a]&DM==o.m[b]&DM for a,b in J.cum)
  ok&=all(o.m[a+8]==o.m[b+8]and o.m[a]&DM==o.m[b]&DM for a,b in J.cum);a8={}
  def find(u):
   while a8.get(u,u)!=u:u=a8[u]
   return u
  for a,b in J.cum:a8[find(a)]=find(b)
  ad={}
  for k2,k3,a in J.seen:ad.setdefault((k3,k2&DM),[]).append(a)
  true_pairs=[(r[i],r[i+1])for r in ad.values()for i in rn(len(r)-1)]
  if ncl>len(J.cum):ok&=all(find(a)==find(b)for a,b in true_pairs)and len(true_pairs)>0
  ok&=all(c==f_ for _,c,f_ in cq)and o.m[CNT]==len(cq)-3;a9=[(o.m[ae+16*k],o.m[ae+16*k+8])for k in rn(3)]==ah[0];af=sum(1 for _,k3,_ in J.seen if k3&0x30);ok&=a9 and af>0 and len({k3&(NB-1)for _,k3,_ in J.seen})>16;DHR.pop(id(o.m),None);return dc(ok=ok,groups_run=al,candidates=len(J.cum),stopped=stopped,stop_pair_listed=any(o.m[a+8]==o.m[b+8]for a,b in J.cum),true_pairs=len(true_pairs),list_skips=J.nls0,old_x_served=a9,records_with_bucket_bits_4_5=af,fast_slow_blocks=J.paths,listings_direct_late=J.lists)
 a,b,c=rcs(True),rcs(False),rcs(False,ncl=1);return dc(ok=slow_wrap(E,a6,B,Q)and a["ok"]and b["ok"]and c["ok"]and a["groups_run"]==[0,1]and b["groups_run"]==[0,1,2]and a["stopped"]and a["stop_pair_listed"]and a["candidates"]>0 and b["candidates"]>0 and c["candidates"]==1 and c["list_skips"]>0,old_x_served=a["old_x_served"]and b["old_x_served"]and c["old_x_served"],stop_cands=a["candidates"],full_cands=b["candidates"],skips=c["list_skips"])
def slow_wrap(E,P,L,c):
 m=dd();hz(E,(P,L),m);z0(E,(P,L),m,CNT);m[Node.E2]=M64-31;N=c.NEXT_BASE;hms(Node(m),64,N);E.arun(P,L,"dlvh",dc(dlvh0(E,P,L,m,CNT),x29=CNT),LV(m,[0,M64-63,96,Node.E2]+[v for j in rn(128)for v in(0,10+j,0,0)],[]),ls64=([WIN],[(M64-4095,1<<64),(0,8192)]));DHR.pop(id(m),None);return[m[c.UP_BASE+88],m[c.UP_BASE+96],m[N-8],m[N+8]]==[(N-16)%2**64,N,11,12]
def bpt(E,rt,lk,rg,Q):
 cq,A=[],{};J=ND();ad=(0,H1,1,2,H1+1);xs=[[E.planted_x(rg,d,rg.rr(4),j=(t+o)%3)[0]for o in E.O]for t,d in en(ad)];ds=[t for t,n in en((570,200,128,2,0))for _ in rn(n)];rg.shuffle(ds);q2=len(ds)//6;leaves=[[xs[t][i%3]for i,t in en(ds[6*l:6*l+6])]for l in rn(q2)];sb=J.SB-E.ENTRY*q2;plant(E,J.m,sb,leaves)
 for t0,n,first in((0,3,LO1+32),(8*H1,2,UP1+32)):rk(E,*rt,"init1",dc(x0=n,x1=t0,x2=first,x3=B1),J,cq,Q.INIT1_BKT*n+1)
 I=Sender1(Q);meter,fl=kref(E,J.m,sb,q2,7,7<<26,I,Q);o=Port1(J.m);rk(E,lk[0],lk[1],"keys",dc(x2=13,x3=sb+E.KEYS0),J,cq,Q.KEY_FIX+Q.KEY_INC*q2+meter+I.extra,o,ls64=([(sb,J.SB)],[]));ok=o.blocks==I.msgs and state1_ok(J.m,I);nb=len(o.blocks)
 for t0,n in((0,3),(8*H1,2)):rk(E,*rt,"fin1",dc(x7=t0,x1=n),J,cq,1+Q.FIN1_X*n,o)
 a9=I.finals(sorted(ad));ok&=o.blocks[nb:]==a9;ok&=[f[:2]for f in fl]==[(t%3,ad[t])for l in rn(q2)for s in(0,1)for u in E.KO for t in(ds[6*l+3*s+u],)];A["keys_fin1"]=dc(records=len(ds),blocks_sent=nb,finals=len(a9),empty_finals=sum(m[1]==0 for m in a9),upper_half_blocks=sum(m[0]>=8*H1 for m in o.blocks),model_ok=ok);J=Node();a5=[rg.gb(41)for _ in rn(4)];aa=[K1]*12+[0,1,2,3,7,15,31,32,33,62,63]+[K1]*20;rg.shuffle(aa);dl=[a5[t]for t,n in en((1200,400,200,50))for _ in rn(n)];rg.shuffle(dl);end,m=J.R2+2*32*K2B,{a5[0]:1,a5[1]:2,a5[2]:1,a5[3]:1};I=Snd(lambda d:end-32*K2B*m[d],end,Q)
 for d in a5:J.m[LC2+8*d],J.m[(LC2+8*d)|RC2]=blk2(d)+32,end-32*K2B*m[d]
 k,ac=0,[]
 for n in aa:recs=[[rg.gb(23)<<41|dl[(k+t)%len(dl)],rg.gb(64),rg.gb(64),rg.gb(64)]for t in rn(n)];k+=n;ac+=[recs]
 o=Pt2(J.m);a8=0
 for a4 in ac:
  e0=I.extra
  for af in a4:I.put(af[0]&M41,af)
  a8+=len(a4)<K1;dla(E,*rt,J,[0,32*len(a4),1,0]+[w for r in a4 for w in r]+[STL]*4*(K1-len(a4)),end,I.extra-e0,cq,Q,o)
 ab=o.blocks==I.msgs and not o.msgs and sok(J.m,I);ab&=hst(HSR[(J.uid,end)])==((1<<126)+a8-len(ac),(1<<66)-a8)and J.m[DONE]==0;dm,wr={d:{}for d in a5},{d:[]for d in a5}
 for msg in o.blocks:unpack(dm[dest(msg[0])],msg,wr[dest(msg[0])])
 ab&=rlk(dm,wr,I,None);A["h2s"]=dc(records=k,blocks_sent=len(o.blocks),blocks_dropped=I.drops,end_stored_back=I.cursor(a5[0])==end,model_and_rule_ok=ab);A["two_group_finals"]=tg=tgf(E,rt,rg,Q,cq);A["wide_counts"]=wdc(E,rt,Q,cq);a6=ND();a7=Port1(a6.m);rk(E,*rt,"init1",dc(x0=1,x1=0,x2=LO1+32,x3=B1),a6,cq,Q.INIT1_BKT+1)
 if(Q.DELTA_F-Q.FIX_BATCH+1)<=Q.DELTA_F-Q.FIX_BATCH:rk(E,lk[0],lk[1],"keys",dc(x2=13,x3=sb+E.KEYS0),a6,cq,0,a7)
 rk(E,*rt,"fin1",dc(x7=0,x1=1),a6,cq,1+Q.FIN1_X,a7);A["meter_skip_ok"]=(a7.blocks==Sender1(Q).finals([0])and dla(E,*rt,Node(),a7.blocks[0],Node.R2+32*K2B,0,cq,Q,None,0,1));A["cnt_eq"]=all(c==f for _,c,f in cq);ae=(A["cnt_eq"]and ok and ab and A["meter_skip_ok"]and A["keys_fin1"]["empty_finals"]==2 and A["keys_fin1"]["upper_half_blocks"]>2 and A["h2s"]["blocks_sent"]>=1 and A["h2s"]["blocks_dropped"]>=1 and A["h2s"]["end_stored_back"]and tg["ok"]and A["wide_counts"]);return dc(ok=ae,meter_skip_ok=A["meter_skip_ok"],dropped=A["h2s"]["blocks_dropped"])
def wdc(E,I,Q,cq):
 o,end=Node(),Node.R2+32*K2B;key=(o.uid,end);run=lambda c:dla(E,I[0],I[1],o,[0,0,c,0]+[STL]*4*K1,end,0,cq,Q,None,0,1<<66)
 def put(D,f):HSR[key]=dc(HSR[key],x13=D%2**64,x14=D>>64,x15=f%2**64,x5=f>>64)
 ok=not run(1)and hst(HSR[key])==(1<<126,(1<<66)-1);put((1<<126)+2**64,1);ok&=not run(1)and hst(HSR[key])==(1<<64,0);put((1<<126)-1,1);return ok and run(2)and o.m[DONE]==1 and hst(HSR[key])==(1<<126,1<<66)
def tgf(E,rt,rg,Q,cq):
 J=Node();d,base=(1<<41)-1,M64-16383;end=(base+60*32*K2B)&M64;J.m[LC2+8*d],J.m[(LC2+8*d)|RC2]=blk2(d)+32,base;o=Snd(lambda _d:base,end,Q);dm,wr,finals,recs={},[],[],[]
 for n in(300,300):
  group=[[rg.gb(23)<<41|d,rg.gb(64),rg.gb(64),rg.gb(64)]for _ in rn(n)];ab=len(o.msgs);recs+=[group];a5=Pt2(J.m)
  for c in rn(0,n,K1):
   e0=o.extra;a4=group[c:c+K1]
   for ac in a4:o.put(d,ac)
   dla(E,*rt,J,[0,32*len(a4),1,0]+[w for r in a4 for w in r]+[STL]*4*(K1-len(a4)),end,o.extra-e0,cq,Q,a5)
  a8=a5.blocks==o.msgs[ab:]and not a5.msgs and sok(J.m,o);I=Pt2(J.m);rk(E,*rt,"fin",dc(x0=MB2,x7=LC2+8*d,x3=(LC2+8*d)|RC2,x8=blk2(d)+32,x1=1,x2=J.E2,x6=end),J,cq,o.fin_count([d]),I);a9,aa=o.final(d);a8&=I.blocks==([a9]if a9 else[])and I.msgs==[[LC2+8*d,J.E2,aa,0,0,0]]
  for m in a5.blocks+I.blocks+I.msgs:unpack(dm,m,wr if len(m)>6 else None)
  finals+=[(aa,a8,a5.blocks)]
 f1,f2,a6=finals[0][0],finals[1][0],finals[1][2];q=300//K2B*K2B;a7=len(a6)>=1 and a6[0][1]==(base+32*q)&M64 and a6[0][4:4+4*(300-q)]==[w for r in recs[0][q:300]for w in r];ok=(all(a8 for _,a8,_ in finals)and rlk({d:dm},{d:wr},o,{d:f2})and a7 and(f1,f2,dm[J.E2])==((base+9600)&M64,(base+19200)&M64,f2)and o.cursor(d)==J.m[(LC2+8*d)|RC2]==(base+600//K2B*32*K2B)&M64 and all([dm.get((base+32*t+8*w)&M64)for w in rn(4)]==recs[0][t]for t in rn(300))and all([dm.get((f1+32*t+8*w)&M64)for w in rn(4)]==recs[1][t]for t in rn(300)));return dc(ok=ok,first_final_rewritten_identically=a7)
def gtc(E,rt,Q):
 rb,rq=rt;K=Q.K_GATE;a5,J=K&(2**64-1),K>>64;cases=[(0,True),(K-1,True),(K,False),(K+1,False),(((J-1)<<64)|(2**64-1),True),((J<<64)|(2**64-1),False),((1<<128)-1,False),];ok=True;I=set()
 for u,a4 in cases:o={};n=E.arun(rb,rq,"gate",dc(x0=u&(2**64-1),x1=u>>64,x2=a5,x3=J),{},out_regs=o);I.add(n);ok&=(o["x0"]==(1 if a4 else 0))
 ok&=(I=={4}and 2*Q.COIN_SCALAR+4+8<=Q.GATE_INSNS);return dc(ok=ok,counts=sorted(I),cases=len(cases))
def p0x(E,rt,rg,Q):
 rb,rq=rt;g0,g1,I=ND(),ND(),ND();CNB=90*2**40;YRB=E.YR0
 for nd in(g0,g1,I):nd.m[CNB]=rg.gb(64);nd.m[CNB+56]=rg.gb(64)
 aa=rg.gb(193)<<63;d1=(aa>>102)&M41;d2=(aa>>143)&M41;A0=aa+0x123456789ABCDEF;A1=A0+(1<<61);ae=A0^(1<<102);rl=ND(sp=0x110000);hm=ND(sp=LINK+0x110000);cq=[]
 for nd in(g0,g1,I):first=(LO1+32)if d1<H1 else(UP1+32);off=d1 if d1<H1 else d1-H1;rk(E,rb,rq,"init1",dc(x0=1,x1=8*d1,x2=first+off*B1,x3=B1),nd,cq,Q.INIT1_BKT+1)
 a5=M41<<3;a8=777;af=778
 for u in rn(64):Av=A0 if u==63 else((A0^((u+1)<<185))&~(M41<<102))|(d1<<102);E.put(g0.m,YRB+u*512+288,Av)
 E.put(g1.m,YRB+288,A1);E.put(I.m,YRB+288,A0);E.put(I.m,YRB+512+288,ae);a6=d1^1;am=(LO1+32+a6*B1)if a6<H1 else(UP1+32+(a6-H1)*B1);rk(E,rb,rq,"init1",dc(x0=1,x1=8*a6,x2=am,x3=B1),I,cq,Q.INIT1_BKT+1);p0=Port1(g0.m);aq=E.arun(rb,rq,"xemit",dc(x1=YRB,x2=101,x26=3<<20,x3=64,x4=32,x6=a8<<49,x8=a5,sp=g0.SP),g0.m,dev=p0.dev);p1=Port1(g1.m);n1=E.arun(rb,rq,"xemit",dc(x1=YRB,x2=202,x26=9<<20,x3=1,x4=32,x6=a8<<49,x8=a5,sp=g1.SP),g1.m,dev=p1.dev);p2=Port1(I.m);b1=E.arun(rb,rq,"xemit",dc(x1=YRB,x2=303,x26=0,x3=1,x4=32,x6=af<<49,x8=a5,sp=I.SP),I.m,dev=p2.dev);b0=E.arun(rb,rq,"xemit",dc(x1=YRB+512,x2=303,x26=4,x3=1,x4=32,x6=a8<<49,x8=a5,sp=I.SP),I.m,dev=p2.dev);ok=(aq==64*Q.XEMIT_HOT+Q.FLUSH1_X+1 and n1==Q.XEMIT_HOT+1 and b1==Q.XEMIT_HOT+1 and b0==Q.XEMIT_HOT+1 and len(p0.blocks)==1);rk(E,rb,rq,"fin1",dc(x7=8*d1,x1=1),g0,cq,1+Q.FIN1_X,p0);rk(E,rb,rq,"fin1",dc(x7=8*d1,x1=1),g1,cq,1+Q.FIN1_X,p1);rk(E,rb,rq,"fin1",dc(x7=8*d1,x1=1),I,cq,1+Q.FIN1_X,p2);a9=Node.R2+4*32*K2B;ar=sorted({d2}|{(g0.m[(LO1 if d1<H1 else UP1)+(d1 if d1<H1 else d1-H1)*B1+32+32*u])&M41 for u in rn(64)})
 for a4 in ar:rl.m[LC2+8*a4]=blk2(a4)+32;rl.m[(LC2+8*a4)|RC2]=Node.R2
 ah=Pt2(rl.m);a7=Snd(lambda d:Node.R2,a9,Q,Node.E2);b2=p0.blocks+p1.blocks+p2.blocks
 for m in b2:
  e0=a7.extra
  for t in rn(m[1]//32):ap=list(m[4+4*t:8+4*t]);a7.put(ap[0]&M41,ap)
  dla(E,rb,rq,rl,m,a9,a7.extra-e0,cq,Q,ah,a7.tag,3)
 ab=Pt2(rl.m)
 for a4 in ar:rk(E,rb,rq,"fin",dc(x0=MB2,x7=LC2+8*a4,x3=(LC2+8*a4)|RC2,x8=blk2(a4)+32,x1=1,x2=Node.E2B,x11=a7.tag,x6=a9),rl,cq,a7.fin_count([a4]),ab)
 CC=(NEXT_BASE-asr1(Node.R2))%2**64;hz(E,rt,hm.m);hm.m[Node.E2]=Node.R2;hms(hm,64,CC);hj=HJ(Q,CC,64)
 for m in ah.blocks+ab.blocks+ab.msgs:
  if dest(m[0])!=d2:continue
  if len(m)==6:dlr(E,rb,rq,hm.m,m,cq,Q,cnt=CNTE)
  else:
   f,ay=hj.block(m,hm.m[Node.E2],2**Q.HB);ok&=ay;dlr(E,rb,rq,hm.m,m,cq,Q,fused=True,r2=(Node.R2,a9),formula=f)
 aj=hm.m[Node.E2B];mh=(aj-Node.R2)//32;ok&=(len(hj.cum)==1 and mh>=2);ak,ao=hj.cum[0];ok&=({hm.m[ak+16],hm.m[ao+16]}=={101,202});ai=E.arun(rb,rq,"hhclr",dc(x1=Node.R2,x2=aj,x10=HT,sp=hm.SP),hm.m);ok&=(ai==3+6*mh and all(hm.m[HT+8*b]==0 for b in hj.heads));wm=dd(int);wk=(5,9,11);wm[M64-55],wm[M64-23],wm[8]=wk;wm.update({HT+8*k:7 for k in wk});ok&=E.arun(rb,rq,"hhclr",dc(x1=M64-63,x2=32,x10=HT),wm)==21 and not any(wm[HT+8*k]for k in wk);ok&=all(c==f for _,c,f in cq);au,ag,J,B0=34*2**40+0x20000,34*2**40+0x40000,(0x987654321<<192)+0x77700000000,(0x987654321<<192)
 def wrep(m,i,a0,a1,c0,c1):
  o=au+128*i;E.put(m,o,a0);E.put(m,o+32,a1);E.put(m,o+64,c0,2);E.put(m,o+80,c1,2)
  for w in rn(4):m[o+96+8*w]=64*i+w+1
 def rps(p_,l_,n,m):o={};s=E.arun(p_,l_,"psel",dc(x0=au,x1=n,x2=ag,sp=g0.SP),m,out_regs=o);return o["x0"],s
 ac=1<<80
 for i,rv in en([(0,9,1,2),(J+50,J+10,1,2),(J+50,J+50,1,2),(J+50,J+50+(1<<63),1,2),(J+50,J+50+2**64,1,2),(J+50,J+90,ac|7,ac|7),(J+100,J+200,10,11),(J+300,J+400,11,12),(J+100,J+500,12,13),(J+600,J+200,12,13),(J+100+(1<<192),J+200+(1<<192),ac|10,ac|11),(J+200+(1<<192),J+250+(1<<192),14,15),(J+(2<<192),J+(2<<192)+((1<<63)-1),14,15),(B0+(3<<192)+(1<<63)-5,B0+(3<<192)+(1<<63)+5,16,17),(J+(4<<192),J+(4<<192)+9,18,19)]):wrep(g0.m,i,*rv)
 az,av=rps(rb,rq,15,g0.m);ax,_=rps(rb,rq,13,g0.m);z0,_=rps(rb,rq,0,g0.m);sa=[E.get(g0.m,ag+128*q+32*s)for q in rn(4)for s in(0,1)];sc=[E.get(g0.m,ag+128*q+64+16*s,2)for q in rn(4)for s in(0,1)];ad=ND()
 for q in rn(3):wrep(ad.m,q,B0+((q+1)<<192)+10,B0+((q+1)<<192)+20,2*q+1,2*q+2)
 for i in rn(3,1030):wrep(ad.m,i,B0+(3<<192)+15,B0+(3<<192)+20,999+2*i,1000+2*i)
 aw,al=rps(rb,rq,1030,ad.m);mp,ml=E.a64.a64_load(E.ksrc("route").replace("tbnz x12,#63,5f\n",""));gm=ND()
 for q in rn(4):wrep(gm.m,q,J+(q<<192)+10,J+(q<<192)+10+(1<<63),2*q+1,2*q+2)
 ok&=(az==1 and ax==0 and z0==0 and aw==0 and rps(mp,ml,4,gm.m)[0]==1 and rps(rb,rq,4,gm.m)[0]==0 and al<=12+220*1024<=Q.PSEL_POOL and((sa[4]>>63)!=(sa[5]>>63))and((sa[6]>>63)!=(sa[7]>>63))and len(set(sa))==8 and sc==[10,11,ac|10,ac|11,14,15,16,17]and all(sa[2*q]>0 and 0<sa[2*q+1]-sa[2*q]<(1<<63)for q in rn(4)));TR.extend([("G",g0.m),("G",g1.m),("G",I.m),("G",ad.m),("R",rl.m),("H",hm.m)]);return dc(ok=ok,xemit_64_insns=aq,hhclr_insns=ai,home_records=mh,candidates=len(hj.cum),psel_ok_steps=av,psel_wc_steps=al)
def wait_dispatch_checks(E,rt,rg,Q):
 rb,rq=rt;nd=ND();o=LINK_GEN+0x10000;a4=[1,1,1];z0(E,rt,nd.m,o)
 def wkc(mem,a5,ev):
  assert a4,"unexpected extra WFE in wcnt";mem[o]+=a4.pop(0)
  if a5[0]==(o>>11):a5[0]=None;ev[0]=1
 a9=E.arun(rb,rq,"wcnt",dc(x0=o,x1=0),nd.m);a8=E.arun(rb,rq,"wcnt",dc(x0=o,x1=3),nd.m,wfe=wkc);a4=[1];nd.m[o]=0;ac=E.arun(rb,rq,"wcnt",dc(x0=o,x1=1),nd.m,wfe=wkc,ev0=1);ok=(a9==4 and a8==4+3*Q.WAKE and ac==4+2*Q.WAKE and nd.m[o]==1 and not a4);a6=(1<<45)+(257<<31);HDR=LINK_GEN+2048;CNTY=LINK_GEN+0x12000;I=ND();z0(E,rt,I.m,CNTY,CNTY+16,CNTY+24);I.m[CNTY+8]=2;o0={};E.arun(rb,rq,"ydlv0",dc(x15=a6+64,x29=CNTY),I.m,out_regs=o0);J=[[rg.gb(64)for _ in rn(56)]for _ in rn(2)];pend=[(rq["ydlv"][0],J[0]),(rq["ydlve"][0],J[1])];py=PtY(I.m);ad=([(LINK_GEN,LINK_GEN+512),(a6,a6+4096)],[(a6,a6+4096),(PYB,PYB+512)])
 def wake_ydsp(mem,a5,ev):
  assert pend,"unexpected extra WFE in ydsp";ab,pl=pend.pop(0)
  for w,v in en(pl):mem[LINK_GEN+64+8*w]=v
  mem[HDR]=ab
 a7=E.arun(rb,rq,"ydsp",dc(o0,x24=a6+4096,sp=I.SP),I.m,dev=py.dev,ls64=ad,wfe=wake_ydsp);aa=3+(Q.YDLV+12)+(Q.YDLVE+12);ok&=(a7==aa and not pend and 9<=Q.WY and I.m[CNTY]==2 and I.m[HDR]==0 and len(py.msgs)==1 and all([I.m[a6+512*u+64+8*w]for w in rn(56)]==J[u]for u in(0,1)));TR.extend([("G",nd.m),("G",I.m)]);return dc(ok=ok,wcnt_0=a9,wcnt_3=a8,ydsp_2=a7,wake=Q.WAKE,wy=Q.WY)

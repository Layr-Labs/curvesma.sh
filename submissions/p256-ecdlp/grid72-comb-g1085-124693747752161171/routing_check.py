"""Transport-only harness; synthetic keys isolate communication and the receiver hash join from curve arithmetic."""
import json,random
from pathlib import Path
from types import SimpleNamespace
import a64,flip_witness
ROOT=Path(__file__).resolve().parent
SENTINEL=0x5E5E5E5E5E5E5E5E
M35=(1<<36)-1
M48=0xffffe0000fffffff# K3 bits compared: 0..28 and 46..63 (bits 29..45: fraction bits 0..16, below the key)
LINK=1<<40
rg=random.Random(314159)
CERT=SimpleNamespace(DLV=6)
def put(mem,ptr,v):
 for k in range(4):mem[ptr+8*k]=(v>>(64*k))&((1<<64)-1)
class Port:
 def __init__(self,n=6):
  self.words,self.msgs,self.n=[0]*n,[],n
 def store(self,k,v):
  self.words[k//8]=v
  if k==8*self.n-8:
   self.msgs.append(list(self.words))
class LinkView:
 def __init__(self,mem,msg):
  self.mem,self.msg=mem,msg
 def __getitem__(self,a):
  return self.msg[(a-LINK)//8]if LINK<=a<LINK+8*len(self.msg)else self.mem[a]
 def __setitem__(self,a,v):
  assert not LINK<=a<LINK+8*len(self.msg)
  self.mem[a]=v
def link(prog,labels,mem,msg,calls,cert):
 calls.append(("dlv",a64.a64_run(prog,labels,"dlv",dict(x0=LINK,x30=SENTINEL),LinkView(mem,msg)),cert.DLV))
class Node:
 P1,P2,HB,SB,SP=0x10000,0x20000,7<<61,1<<34,0x7000000
 R1,E1,R2,E2=0x1000020,0x2000000,0x3000000,0x4000000
 BK,CL=1<<36,0x6200000
 W2,CT,NX=0xC000000,1<<39,1<<37# hop2 payload windows, local cursor table, home next array (16 bytes a record)
 DC,LE=0x6300000,0x6400000
 def __init__(self):
  self.mem={}
def run_k(prog,labels,entry,regs,nd,calls,formula,port=None,ls64=None):
 o={}
 regs=dict(regs,sp=nd.SP,x30=SENTINEL)
 dev=(port[0],8*port[1].n,port[1].store)if port else None
 calls.append((entry,a64.a64_run(prog,labels,entry,regs,nd.mem,out_regs=o,dev=dev,ls64=ls64),formula))
 return o
def k3f(k2,k3):return(k3&M48)|((k2>>58)<<28)
def src(rec):return(rec[2]&M35,rec[3]&M35)
def verf(st,prior,new,rec,listcap):
 pk=st['key'][prior]
 if(pk[1]^rec[1])&M48:return 6,False
 if st['src'][prior]==src(rec)or len(st['list'])>=listcap:return 13,False
 st['list'].append((prior,new));return 17,True
def join_model(st,p,rec,hmask,listcap):
 # one record's join and walk (route.s JOIN, hit*); p is the record's w2 address
 k2=rec[0];h=k2&hmask;cur=st['head'].get(h,0);count=4;visits=same=0;listed=None;full=False;hits=0
 st['next'][p]=cur;st['head'][h]=p;st['key'][p]=(k2,rec[1]);st['src'][p]=src(rec);st['seq'].append(p)
 if not cur:return count,listed,visits,same,full,hits
 def hit(q):
  nonlocal same,full,listed,hits
  hits+=1;n,ok=verf(st,q,p,rec,listcap)# 6, 13 or 17 with the routine's RET
  if ok:listed=(q,p)
  elif n==13:
   if st['src'][q]==src(rec):same+=1
   else:full=True
  return n
 a=cur
 while a:# one visit: LDP, CMP, B.NE, then CBZ/CBNZ on the next entry; a K2 hit calls hit5*/hit6*
  visits+=1;count+=4
  if st['key'][a][0]==k2:
   count+=1+hit(a)# BL, then the routine
   if listed:break
  a=st['next'][a]
 return count,listed,visits,same,full,hits
def check():
 cert=CERT;rg=random.Random(32452843);prog,labels=a64.a64_load((ROOT/'route.s').read_text())
 bp,bl=a64.a64_load(flip_witness.block_production(ROOT))
 side,q2,cap1,hbits,caph,ncl,G,planned=4,65,32,6,180,56,2,3
 st1,hmask=32*cap1,(1<<58)-1
 ids=[(i,j)for i in range(side)for j in range(side)];nodes={n:Node()for n in ids};calls=[]
 homes={n:[]for n in ids};drops=dict(hop1=0,home=0)
 msgs_ok=slots_ok=chain_ok=tags_ok=True;shared=0;acks=0;dones=0;completions=0
 state={n:dict(head={},next={},key={},src={},list=[],seq=[])for n in ids};lcur={n:nodes[n].CL for n in ids}
 paths=dict(empty_chain=0,chain_end=0,same_node=0,listed=0,list_full=0,home_full=0,visits=0,k2_hits=0,odd_blocks=0)
 SLOTS=(Node.R1,Node.E1);H2P=lambda nd:[(nd.P2,nd.P2+64)]
 def h2f(n):
  P,L=(n-1)//2,(n-1)%2
  return 7 if n==0 else 25+14*P+8*L+2*(P//32)
 RB=Node.BK+32;ent=lambda p:Node.NX+(p-16-RB)//2
 rx={n:dict(x0=LINK,x2=nodes[n].LE,x3=nodes[n].CL,x4=nodes[n].CL+16*ncl,x9=(RB-2*Node.NX)%(1<<64),x10=RB,x11=RB+32*caph,x12=Node.HB,x13=nodes[n].DC,x14=Node.NX,x28=0,x29=side)for n in ids}
 for n,nd in nodes.items():
  for d in range(side):nd.mem[Node.CT+8*d]=Node.W2+8192*d+32
  for w in range(16384):nd.mem[nd.HB+8*w]=0x5A5A5A5A
  run_k(prog,labels,'zero16384',dict(x0=16384,x1=nd.HB),nd,calls,1026+1)
  assert all(nd.mem[nd.HB+8*w]==0 for w in range(16384))
 def planted():
  d1,d2=rg.randrange(side),rg.randrange(side)
  top=(0xffffffff00000000|rg.getrandbits(32))if rg.random()<0.2 else rg.getrandbits(63)
  k2=rg.getrandbits(hbits)|rg.getrandbits(4)<<60
  return d1|(d2<<64)|(k2<<128)|(top<<192)
 pool=[planted()for _ in range(120)]
 pool2=[k^(1<<192)for k in pool[:30]]+[k^(1<<(192+30))for k in pool[30:60]]
 stop_flag=0;groups_run=0;group_listings=[]
 for group in range(planned):
  if stop_flag:break
  groups_run+=1
  for n in ids:nodes[n].mem[Node.DC]=0
  group_listed={n:[]for n in ids}
  ref1={n:{i:[]for i in range(side)}for n in ids}
  for(i,j),nd in nodes.items():
   run_k(bp,bl,'block_init',dict(x0=side,x1=nd.P1+64,x2=0x9000000,x3=0x8000020,x4=nd.R1+i*st1),nd,calls,775*side+2)
   batch=[]
   for l in range(q2):
    pair=[]
    for sign in range(2):
     u=rg.random()
     if u<0.3:key=rg.choice(pool);shared+=1
     elif u<0.4:key=rg.choice(pool2);shared+=1
     else:key=planted()
     pair.append(key);put(nd.mem,nd.SB+64*l+32*sign,key)
    batch.append(pair)
    nd.mem[0x7100000+72*l]=l<<36
   W2=i|((group*side*side+side*i+j)<<36);exp=[];cur={};full=0
   for l,pair in enumerate(batch):
    for sign,key in enumerate(pair):
     d1=key&M35;d2=(key>>64)&M35;ptr=cur.get(d1,nd.R1+i*st1);cur[d1]=ptr+32
     if ptr>=nd.R1+(i+1)*st1:full+=1;continue
     exp.append([d1,ptr,(key>>128)&((1<<64)-1),key>>192,W2,d2|(l<<36)])
   port_words=[0]*4;block_msgs=[];ack_events=[]
   def block_store(offset,value):
    port_words[offset//8]=value
    if offset!=24:return
    d,addr,source,count=port_words
    assert count==128 and source==0x8000020+8192*d
    target=nodes[(d,j)].mem;ack=0xA000000;before=len(ack_events)
    class BlockView:
     def __getitem__(self,a):
      if LINK<=a<LINK+32:return port_words[(a-LINK)//8]
      if a==LINK+32:return ack
      if source<=a<source+4096:return nd.mem[a]
      return target[a]
     def __setitem__(self,a,v):
      assert not LINK<=a<LINK+40 and not source<=a<source+4096
      target[a]=v
      if a==ack:
       assert v==1;ack_events.append((d,addr))
    count=a64.a64_run(bp,bl,'block_dlv',dict(x0=LINK,x30=SENTINEL),BlockView(),ls64=([(source,source+4096)],[SLOTS]))
    calls.append(('block_dlv',count,266));assert len(ack_events)==before+1 and ack_events[-1]==(d,addr)
    block_msgs.append((d,addr))
   emitter=a64.a64_run(bp,bl,'pair_keys',dict(x0=q2,x1=nd.SB,x2=nd.P1,x3=W2,x4=nd.R1+(i+1)*st1,x5=0x7100000,x6=0x9000000,x30=SENTINEL),nd.mem,dev=(nd.P1,32,block_store))
   total_full=sum((ptr-(nd.R1+i*st1))//4096 for ptr in cur.values())
   assert emitter<=93*q2//4+30+8*total_full
   calls.append(('block_pair_keys',emitter,emitter))
   drops['hop1']+=full
   for m in exp:ref1[(m[0],j)][i].append(tuple(m[2:]))
   port=Port();end=nd.R1+(i+1)*st1
   def finish_store(offset,value):
    if offset<32:block_store(offset,value)
    elif offset>=64:port.store(offset-64,value)
   ct=a64.a64_run(bp,bl,'block_finish',dict(x0=side,x1=nd.P1+64,x2=0x9000000,x3=0x8000020,x4=0xB000000,x5=end,x7=0xB000040,x8=nd.E1+32*i,x30=SENTINEL),nd.mem,dev=(0xB000000,112,finish_store))
   assert ct<=40*side+4;calls.append(('block_finish',ct,ct))
   msgs_ok&=port.msgs==[[d,nd.E1+32*i,min(cur.get(d,nd.R1+i*st1),end),0,0,0]for d in range(side)]
   for m in port.msgs:link(prog,labels,nodes[(m[0],j)].mem,m,calls,cert)
  for(a,j),nd in nodes.items():
   expected=[]
   for i in range(side):
    got=[tuple(nd.mem[x+8*k]for k in range(4))for x in range(nd.R1+i*st1,nd.mem[nd.E1+32*i],32)]
    slots_ok&=got==ref1[(a,j)][i]
    for rec in ref1[(a,j)][i]:
     d2=rec[3]&M35;out=list(rec);out[3]=(out[3]&~M35)|j;expected.append([d2,0]+out)
   def deliver(snd,home,hdr,group_listed=None):
    nonlocal acks,dones,completions,slots_ok,chain_ok
    dst=nodes[home];r=rx[home];win,cnt,done=hdr[1],hdr[2],hdr[0]>>63;ack=0xA000000;seen=[]
    recs=[tuple(snd.mem[win+32*k+8*w]for w in range(4))for k in range(cnt)]
    class V:
     def __getitem__(self,a):
      if LINK<=a<LINK+32:return hdr[(a-LINK)//8]
      if a==LINK+32:return ack
      if win<=a<win+4096:return snd.mem[a]
      return dst.mem[a]
     def __setitem__(self,a,v):
      assert not LINK<=a<LINK+40 and not win<=a<win+4096
      dst.mem[a]=v
      if a==ack:assert v==1;seen.append(1)
    room=r['x10']+32*cnt+32<=r['x11'];f=7;cur=r['x10'];st=state[home]
    if room and cnt:
     f+=1+4+1+6+1+1;odd=cnt%2
     for k,rec in enumerate(recs):
      p=cur+16;cur+=32;homes[home].append((p,rec))
      jc,listed,v,same,full,hits=join_model(st,p,rec,hmask,ncl);f+=jc
      if listed and group_listed is not None:group_listed[home].append(listed)
      paths['listed']+=bool(listed);paths['visits']+=v;paths['list_full']+=full;paths['same_node']+=same;paths['k2_hits']+=hits
      if v==0:paths['empty_chain']+=1
      elif not listed and not full:paths['chain_end']+=1
     f+=4*((cnt-1)//2)+4*((cnt-1)%2)
     if odd:f+=2;cur+=32;paths['odd_blocks']+=1
    elif cnt:
     drops['home']+=cnt;paths['home_full']+=1
    else:f+=1
    final=done and r['x28']+1==side;f+=3+1+1+(1 if not done else 7 if final else 4)
    out={}
    count=a64.a64_run(prog,labels,'rblk',dict(r,x30=SENTINEL,sp=dst.SP),V(),out_regs=out,ls64=([(win,win+4096)],[(Node.BK,r['x11']+64)]))
    calls.append(('rblk',count,f));acks+=1;assert seen==[1]and out['x30']==SENTINEL and out['sp']==dst.SP
    slots_ok&=out['x10']==cur and out['x14']==ent(cur+16)and out['x12']>>61==7 and out['x12']&7==0
    chain_ok&=out['x3']==dst.CL+16*len(st['list'])
    if done:
     dones+=1;assert out['x28']==(0 if final else r['x28']+1)and dst.mem[dst.DC]==(1 if final else 0)
     if final:completions+=1;chain_ok&=dst.mem[dst.LE]==out['x3']and dst.mem[dst.LE+8]==out['x10']
    else:assert dst.mem[dst.DC]==0
    for q in('x3','x10','x12','x14','x28'):r[q]=out[q]
   words=[0]*4;sent=[]
   def h2store(offset,value):
    words[offset//8]=value
    if offset==24:
     sent.append(list(words));deliver(nd,(a,words[0]&M35),list(words),group_listed)
   exp_d={d:[]for d in range(side)}
   for e in expected:exp_d[e[0]].append(tuple(e[2:]))
   h2c=4+sum(h2f(len(ref1[(a,j)][i]))for i in range(side))+8*sum(len(v)//128 for v in exp_d.values())
   run_k(prog,labels,'h2b',dict(x0=nd.P2,x3=side,x19=nd.E1,x20=nd.R1,x21=st1,x17=j,x23=Node.CT),nd,calls,h2c,(nd.P2,SimpleNamespace(n=4,store=h2store)),([SLOTS],[]))
   run_k(prog,labels,'h2fin',dict(x0=nd.P2,x1=Node.CT,x2=side,x3=Node.W2+32,x4=1<<63),nd,calls,2+11*side,(nd.P2,SimpleNamespace(n=4,store=h2store)))
   msgs_ok&=[m for m in sent if m[0]>>63]==[[1<<63|d,Node.W2+8192*d+32,len(exp_d[d])%128,0]for d in range(side)]and all(nd.mem[Node.CT+8*d]==Node.W2+8192*d+32 for d in range(side))
   msgs_ok&=all(m[2]==128 and m[1]==Node.W2+8192*m[0]+32 for m in sent if not m[0]>>63)and sum(1 for m in sent if not m[0]>>63)==sum(len(v)//128 for v in exp_d.values())
   for d in range(side):
    got=[rec for p,rec in homes[(a,d)]][-len(exp_d[d]):]if exp_d[d]else[]
    if rx[(a,d)]['x10']+64<=rx[(a,d)]['x11']and not drops['home']:msgs_ok&=got==exp_d[d]
  for home in ids:assert rx[home]['x28']==0 and nodes[home].mem[Node.DC]==1
  listed_now=0
  for home,nd in nodes.items():
   for addr,rec in homes[home]:
    slots_ok&=tuple(nd.mem[addr-16+8*k]for k in range(4))==tuple(rec)
    k=rec[2]>>36;tags_ok&=(k%side)==(rec[3]&M35)and(k//side)%side==(rec[2]&M35)and k//(side*side)<=group and((rec[3]>>36)&511)<q2
   got=[(nd.mem[c],nd.mem[c+8])for c in range(lcur[home],nd.mem[nd.LE],16)];chain_ok&=got==group_listed[home];lcur[home]=nd.mem[nd.LE];listed_now+=len(got)
   for h,addr in state[home]['head'].items():chain_ok&=nd.mem[nd.HB+8*h]==ent(addr)
   for addr,nxt in state[home]['next'].items():chain_ok&=nd.mem[ent(addr)]==(ent(nxt)if nxt else 0)and nd.mem[ent(addr)+8]==state[home]['key'][addr][0]
   for h in range(1<<hbits):
    cur=nd.mem[nd.HB+8*h];walk=[]
    while cur:walk.append(cur);cur=nd.mem[cur]
    members=[ent(addr)for addr,rec in homes[home]if(rec[0]&hmask)==h]
    chain_ok&=walk==[ent(x)for x in reversed(state[home]['seq'])if(state[home]['key'][x][0]&hmask)==h]and sorted(walk)==sorted(members)
   chain_ok&=all((rec[0]&hmask)<1<<hbits for _,rec in homes[home])
  group_listings.append(listed_now)
  if group==1 and listed_now>0:stop_flag=1
 stop_ok=(groups_run==2 and planned==3 and stop_flag==1)
 found=set();exp_c=set()
 for home,nd in nodes.items():
  for c in range(nd.CL,lcur[home],16):found.add((home,nd.mem[c],nd.mem[c+8]))
  exp_c|={(home,a,b)for a,b in state[home]['list']}
  for a,b in state[home]['list']:chain_ok&=state[home]['key'][a][0]==state[home]['key'][b][0]and(state[home]['key'][a][1]&M48)==(state[home]['key'][b][1]&M48)and state[home]['src'][a]!=state[home]['src'][b]
  if len(state[home]['list'])<ncl:
   want=[];seen=[]
   for addr,rec in homes[home]:
    key=(rec[0],rec[1]&M48);p=[x for x,r in seen if(r[0],r[1]&M48)==key and src(r)!=src(rec)]
    if p:want.append((p[-1],addr))
    seen.append((addr,rec))
   chain_ok&=sorted(want)==sorted(state[home]['list'])
 ns=[0,1,2,3,4,5,33,34,64,65,66,97,131,200];nd=Node();exp={0:[],1:[]};blocks=[];sent=[];words=[0]*4
 for d in(0,1):nd.mem[Node.CT+8*d]=Node.W2+8192*d+32
 for i,n in enumerate(ns):
  for r in range(n):
   w=[rg.getrandbits(64)for _ in range(3)]+[rg.getrandbits(28)<<36|rg.randrange(2)];put(nd.mem,nd.R1+8192*i+32*r,w[0]|w[1]<<64|w[2]<<128|w[3]<<192)
   d=w[3]&M35;exp[d].append(tuple(w[:3]+[w[3]&~M35|5]))
   if len(exp[d])%128==0:blocks.append((d,exp[d][-128:]))
  nd.mem[nd.E1+32*i]=nd.R1+8192*i+32*n
 def fstore(offset,value):
  words[offset//8]=value
  if offset==24:sent.append((words[0],[tuple(nd.mem[words[1]+32*k+8*w]for w in range(4))for k in range(words[2])],words[1],words[2]))
 run_k(prog,labels,'h2b',dict(x0=nd.P2,x3=len(ns),x19=nd.E1,x20=nd.R1,x21=8192,x17=5,x23=Node.CT),nd,calls,4+sum(map(h2f,ns))+8*len(blocks),(nd.P2,SimpleNamespace(n=4,store=fstore)),([SLOTS],[]))
 msgs_ok&=[(m[0],m[1])for m in sent]==blocks and len(blocks)>=4 and all(m[2]==Node.W2+8192*m[0]+32 and m[3]==128 for m in sent)
 sent.clear();run_k(prog,labels,'h2fin',dict(x0=nd.P2,x1=Node.CT,x2=2,x3=Node.W2+32,x4=1<<63),nd,calls,2+11*2,(nd.P2,SimpleNamespace(n=4,store=fstore)))
 msgs_ok&=[(m[0],m[1])for m in sent]==[(1<<63|d,exp[d][len(exp[d])//128*128:])for d in(0,1)]
 hn=Node();o={};hr=dict(x0=LINK,x2=hn.LE,x3=hn.CL,x4=hn.CL+16,x9=(RB-2*Node.NX)%(1<<64),x10=RB,x11=RB+64,x12=Node.HB,x13=hn.DC,x14=Node.NX,x28=0,x29=2,x30=SENTINEL,sp=hn.SP)
 ct=a64.a64_run(prog,labels,'rblk',hr,LinkView(hn.mem,[1<<63,Node.W2+32,0,0,0xA000000]),out_regs=o);calls.append(('rblk_empty',ct,17));msgs_ok&=o['x10']==RB and o['x28']==1 and hn.mem[0xA000000]==1
 bb,bd=a64.a64_load((ROOT/'block.s').read_text());erows=[]
 for n,fl in[(1,0),(1,1),(7,1),(8,0),(8,1),(9,0),(9,1)]:
  m1={};SB=1<<34;P1=0x10000;W=0x8000000;RS=0x9000000;IDB=0x7100000;RE=RS+65536
  for j in range(n):
   d1=j%2;m1[SB+64*j]=d1;m1[SB+64*j+8]=j;m1[SB+64*j+16]=j+1;m1[SB+64*j+24]=j+2;m1[SB+64*j+32]=d1;m1[SB+64*j+40]=j+3;m1[SB+64*j+48]=j+4;m1[SB+64*j+56]=j+5;m1[IDB+72*j]=j
  for d1 in[0,1]:
   m1[P1+64+8*d1]=W+4096 if fl else W+32;m1[RS+8*d1]=RE-8192
  m2=dict(m1);c1=[];c2=[]
  r1=a64.a64_run(bp,bl,'pair_keys',dict(x0=n,x1=SB,x2=P1,x3=17|(31<<36),x4=RE,x5=IDB,x6=RS,x30=SENTINEL),m1,dev=(P1,32,lambda k,v:c1.append((k,v))))
  r0=a64.a64_run(bb,bd,'pair_keys',dict(x0=n,x1=SB,x2=P1,x3=17|(31<<36),x4=RE,x5=IDB,x6=RS,x30=SENTINEL),m2,dev=(P1,32,lambda k,v:c2.append((k,v))))
  f=sum(1 for k,v in c1 if k==24)
  assert m1==m2 and c1==c2 and r1==r0+4-14*(n//8),(n,fl,r1,r0)
  assert r1<=93*n//4+30+8*f,(n,fl,r1,f)
  erows.append((n,fl,r1,r0,f))
 counts={}
 for e,c,f in calls:counts[e]=counts.get(e,0)+c
 full_lists=sum(1 for n in ids if len(state[n]['list'])>=ncl)
 return dict(nodes=len(ids),groups_run=groups_run,groups_planned=planned,group_listings=group_listings,leaves=2*2*q2*len(ids)*groups_run,shared_leaves=shared,candidates=len(found),full_lists=full_lists,cands_ok=found==exp_c and len(found)>0,msgs_ok=msgs_ok,slots_ok=slots_ok,tags_ok=tags_ok,chain_ok=chain_ok,stop_ok=stop_ok,drops=drops,calls=len(calls),counts=counts,data_receives=acks,done_markers=dones,completions=completions,signaled=completions==groups_run*len(ids),formulas_ok=all(c==f for _,c,f in calls),paths=paths,all_paths=all(v>0 for v in paths.values()),emit_unroll=erows,scope='Synchronous finite grouped replay with a toy stop rule; asynchronous protocol and the production stop rule require the separate argument.')

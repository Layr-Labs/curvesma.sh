import sys,json,random
from pathlib import Path
root=Path(__file__).resolve().parent
if (root/'files').is_dir():root=root/'files'
sys.path.insert(0,str(root));import a64,scalar_check
cfg=json.load(sys.stdin);rng=random.Random(cfg['seed']);M=(1<<35)-1;W=(1<<64)-1
prog,labels=a64.a64_load((root/'phase4.s').read_text());SENT=0x12345678
BASE=0x100000;SOURCE=0xB00000;CUR=0x20000;ACK=0x30000;HOME=0x10000000;WIN=1<<40;PORT=0x40000;INPUT=0x80000
class View:
 def __init__(self,mem,msg):self.mem=mem;self.msg=msg
 def __getitem__(self,a):return self.msg[(a-WIN)//8] if WIN<=a<WIN+80 else self.mem[a]
 def __setitem__(self,a,v):assert not WIN<=a<WIN+80;self.mem[a]=v
calls=[];drops=[0,0,0,0];completed=partial=0;packets=0;recovered_targets=0
def send(records,axis):
 mem={};out=[];words=[0]*10
 for i,rec in enumerate(records):
  for j,v in enumerate(rec):mem[INPUT+64*i+8*j]=v
 def dev(k,v):
  words[k//8]=v
  if k==72:out.append(list(words))
 count=a64.a64_run(prog,labels,'send64',dict(x0=len(records),x1=INPUT,x2=PORT,x3=8*axis,x30=SENT),mem,limit=100000,dev=(PORT,80,dev))
 assert count==14*len(records)+2
 assert out==[[rec[axis]&M,0]+rec for rec in records]
 calls.append(('send',count));return out
def run(side,capcand,capreq,caprelay,seed,fixture=None):
 global completed,partial,packets,recovered_targets
 rr=random.Random(seed);nodes=[(a,b) for a in range(side) for b in range(side)];mem={n:{CUR:BASE,ACK:1} for n in nodes}
 reference={n:[] for n in nodes};requests={n:[] for n in nodes};payload={};wanted={};home_ref={n:{} for n in nodes}
 xt={n:[[rr.getrandbits(256) for _ in range(3)] for k in range(4)] for n in nodes}
 yt={n:[[rr.getrandbits(256) for _ in range(3)] for l in range(4)] for n in nodes};bs={n:[rr.getrandbits(256) for _ in range(4)] for n in nodes}
 if fixture is not None:
  h,sgn=fixture;secret=1739+h;aa,bb,cc=37,11,13;bb2,cc2=17,19
  lc=aa+bb*secret+cc;aa2=(sgn*lc+cc2-bb2*secret)%scalar_check.R
  X1=scalar_check.mul(aa+bb*secret);X2=scalar_check.mul(aa2+bb2*secret)
  Y1=scalar_check.mul(cc);Y2=scalar_check.mul(cc2)
  def normalized(A,a,b):
   sign=1 if A[1]*scalar_check.MR%scalar_check.P<=scalar_check.P//2 else -1
   A=A if sign==1 else scalar_check.neg(A)
   return [A[0]*scalar_check.MR%scalar_check.P,A[1]*scalar_check.MR%scalar_check.P,a*sign%scalar_check.R],b*sign%scalar_check.R,sign
  xt[(0,0)][0],bs[(0,0)][0],sx1=normalized(X1,aa,bb)
  xt[(1,1)][0],bs[(1,1)][0],sx2=normalized(X2,aa2,bb2)
  yt[(0,0)][0],_,sy1=normalized(Y1,cc,0)
  yt[(1,1)][0],_,sy2=normalized(Y2,cc2,0)
  ftypes=[0 if sx1*sy1==1 else 1,0 if -sx2*sy2==1 else 1]
 for src in nodes:
  def put(addr,v):
   for j in range(4):mem[src][addr+8*j]=(v>>(64*j))&W
  for k,bval in enumerate(bs[src]):put(0x500000+32*k,bval)
  for k,(x,y,a) in enumerate(xt[src]):put(0x200000+64*k,x);put(0x200000+64*k+32,y);put(0x400000+32*k,a)
  for l,(x,y,c) in enumerate(yt[src]):put(0x300000+64*l,x);put(0x300000+64*l+32,y);put(0x600000+32*l,c)
 for home in nodes:
  ct=a64.a64_run(prog,labels,'clear_flags',dict(x0=16*capcand,x1=HOME+16,x30=SENT),mem[home],limit=100000);assert ct==3*16*capcand+2;calls.append(('clear',ct))
  for cid in range(capcand):
   for leafside in range(2):
    src=rr.choice(nodes);k=rr.randrange(4);l=rr.randrange(4);typ=rr.randrange(2)
    if fixture is not None and home==(0,1) and cid==0:
     src=(0,0) if leafside==0 else (1,1);k=l=0;typ=ftypes[leafside]
    key=(src,k,l);payload.setdefault(key,xt[src][k][:2]+yt[src][l][:2]+[xt[src][k][2],bs[src][k],yt[src][l][2]])
    req=[src[0]|(k<<35),src[1]|(typ<<35)|(l<<36),home[0]|(cid<<35),home[1]|(leafside<<35),0,0,0,0]
    requests[home].append(req);wanted[(home,cid,leafside)]=payload[key]
 for home in nodes:
  for cid in range(capcand):
   for sd in range(2):
    addr=0x900000+64*cid+32*sd;req=requests[home][2*cid+sd]
    mem[home][addr]=0;mem[home][addr+8]=0;mem[home][addr+16]=req[0];mem[home][addr+24]=req[1]
    mem[home][0xA00000+16*cid+8*sd]=addr
  ct=a64.a64_run(prog,labels,'make_requests',dict(x0=capcand,x1=0xA00000,x2=SOURCE,x3=home[0],x4=home[1],x30=SENT),mem[home],limit=100000);assert ct==17*capcand+3;calls.append(('make_requests',ct))
  assert [[mem[home][SOURCE+64*i+8*j] for j in range(8)] for i in range(2*capcand)]==requests[home]
 def transfer(inputs,axis,stage,cap,home_handler=False):
  nonlocal reference
  outputs={n:[] for n in nodes};reference={n:[] for n in nodes};storage=SOURCE if stage==1 else BASE;cursor=CUR+8*stage;ack=ACK+8*stage
  if not home_handler:
   for n in nodes:mem[n][cursor]=storage
  order=list(nodes);rr.shuffle(order)
  for sender in order:
   for msg in send(inputs[sender],axis):
    dest=(msg[0],sender[1]) if axis==0 else (sender[0],msg[0]);rec=msg[2:];mem[dest][ack]=1;out={}
    if home_handler:
     cid=rec[0]>>35;sd=(rec[1]>>35)&1;ch=(rec[1]>>36)&7;slot=(2*cid+sd)*8+ch
     assert rec[0]&M==dest[0] and rec[1]&M==dest[1] and rec[2]==1
     assert slot not in home_ref[dest];home_ref[dest][slot]=rec
     ct=a64.a64_run(prog,labels,'home64',dict(x0=WIN,x1=HOME,x3=ack,x30=SENT),View(mem[dest],msg),limit=1000,out_regs=out);assert ct==16
     assert [mem[dest][HOME+64*slot+8*j] for j in range(8)]==rec
    else:
     accepted=len(reference[dest])<cap
     ct=a64.a64_run(prog,labels,'append64',dict(x0=WIN,x1=cursor,x2=storage+64*cap,x3=ack,x18=64,x19=W-63,x30=SENT),View(mem[dest],msg),limit=1000,out_regs=out);assert ct==(13 if accepted else 7)
     if accepted:reference[dest].append(rec)
     else:drops[stage]+=1
     assert mem[dest][cursor]==storage+64*len(reference[dest])
    assert mem[dest][ack]==0 and out['x30']==SENT
    calls.append(('home' if home_handler else 'append',ct))
  if not home_handler:
   for n in nodes:
    outputs[n]=[[mem[n][storage+64*i+8*j] for j in range(8)] for i in range(len(reference[n]))]
    assert outputs[n]==reference[n]
  return outputs
 relay=transfer(requests,0,0,caprelay)
 sources=transfer(relay,1,1,capreq)
 replies={n:[] for n in nodes}
 for src in nodes:
  for rec in sources[src]:
   assert (rec[0]&M,rec[1]&M)==src
   key=(src,rec[0]>>35,rec[1]>>36);vals=payload[key]
   k=key[1];l=key[2];xp=0x200000+64*k;yp=0x300000+64*l;ap=0x400000+32*k;bb=0x500000+32*k;cp=0x600000+32*l;tmp=0x700000;outp=0x800000
   ct=a64.a64_run(prog,labels,'gather224',dict(x0=xp,x1=yp,x2=ap,x3=bb,x4=cp,x5=tmp,x30=SENT),mem[src],limit=1000);assert ct==31;calls.append(('gather',ct))
   ct=a64.a64_run(prog,labels,'build_replies',dict(x0=rec[2],x1=rec[3],x2=tmp,x3=outp,x30=SENT),mem[src],limit=1000);assert ct==84;calls.append(('build',ct))
   got=[[mem[src][outp+64*ch+8*j] for j in range(8)] for ch in range(8)]
   expect=[[rec[2],rec[3]|(ch<<36),1,0]+[(v>>(64*j))&W for j in range(4)] for ch,v in enumerate(vals+[0])]
   assert got==expect;replies[src].extend(got)
 relay=transfer(replies,0,2,caprelay)
 transfer(relay,1,3,0,True)
 for home in nodes:
  for cid in range(capcand):
   expected=all((2*cid+sd)*8+ch in home_ref[home] for sd in range(2) for ch in range(8));out={}
   ct=a64.a64_run(prog,labels,'complete',dict(x0=16,x1=HOME+64*16*cid+16,x30=SENT),mem[home],limit=1000,out_regs=out);assert ct<=67 and bool(out['x0'])==expected;calls.append(('complete',ct))
   if expected:
    completed+=1;reconstructed=[]
    for sd in range(2):
     got=[sum(mem[home][HOME+64*((2*cid+sd)*8+ch)+32+8*j]<<(64*j) for j in range(4)) for ch in range(7)]
     assert got==wanted[(home,cid,sd)];reconstructed.append(got)
    if fixture is not None and home==(0,1) and cid==0:
     recovered=scalar_check.recover(reconstructed[0],reconstructed[1],ftypes[0],ftypes[1])
     assert recovered==secret and scalar_check.mul(recovered)==scalar_check.mul(secret);recovered_targets+=1
   else:partial+=1
 for ins in [requests,sources,replies,relay]:packets+=sum(map(len,ins.values()))
for side in [1,2,3]:
 for cr,cl in [(0,0),(1,1),(3,5),(16,17),(1000,1000)]:run(side,4,cr,cl,cfg['seed']+side+cr+cl)
for h in range(3):
 for sgn in [1,-1]:run(2,4,1000,1000,193+h,(h,sgn))
assert recovered_targets==6
scalar_evidence=scalar_check.check()
mem={ACK:1};rec=[M|(((1<<28)-1)<<35),M|(1<<35)|(7<<36),1,0,11,22,33,44];msg=[M,0]+rec
ct=a64.a64_run(prog,labels,'home64',dict(x0=WIN,x1=HOME,x3=ACK,x30=SENT),View(mem,msg),limit=1000)
slot=(1<<32)-1;assert mem[HOME+64*slot+16]==1 and mem[HOME+64*slot+56]==44 and ct==16
assert completed>0 and partial>0 and drops[0]>0 and drops[1]>0 and drops[2]>0
tree_cases=tree_messages=0
TREE_VALUE=0x123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef0
for side in [1,2,4]:
 for mode in [0,1,2]:
  ns=[(i,j) for i in range(side) for j in range(side)];tree_mem={n:{} for n in ns}
  local={n:(None if mode==0 or mode==1 and n!=ns[-1] else TREE_VALUE) for n in ns}
  for n in ns:
   for slot in range(4):tree_mem[n][0xD00000+4096*slot+2048]=0
  def visit(n):
   global tree_messages
   i,j=n;children=[]
   for k in [0,1]:
    if 2*i+1+k<side:children.append(((2*i+1+k,j),k))
   if i==0:
    for k in [0,1]:
     if 2*j+1+k<side:children.append(((0,2*j+1+k),2+k))
   answer=local[n]
   for child,slot in children:
    value=visit(child);addr=0xD00000+4096*slot
    data=[0 if value is None else (value>>(64*k))&W for k in range(4)]+[0,0,0,2 if value is None else 1]
    sender_mem={0xF00000+8*k:v for k,v in enumerate(data)};window=[0]*10;captured=[]
    def emit(k,v):
     window[k//8]=v
     if k==72:captured.append(list(window))
    dest=n[0] if child[0]>0 else n[1]
    sc=a64.a64_run(prog,labels,'send_completion',dict(x0=dest,x1=addr,x2=0xF00000,x3=PORT,x30=SENT),sender_mem,limit=1000,dev=(PORT,80,emit));assert sc==10
    assert captured==[[dest,addr]+data];message=captured[0]
    ct=a64.a64_run(prog,labels,'publish64',dict(x0=WIN,x30=SENT),View(tree_mem[n],message),limit=1000);assert ct==12;tree_messages+=1
    assert [tree_mem[n][addr+8*k] for k in range(7)]==data[:7] and tree_mem[n][addr+2048]==data[7]
    if answer is None and tree_mem[n][addr+2048]==1:answer=sum(tree_mem[n][addr+8*k]<<(64*k) for k in range(4))
   return answer
  assert visit((0,0))==(None if mode==0 else TREE_VALUE);tree_cases+=1
print(json.dumps(dict(research_only=True,tree_cases=tree_cases,tree_messages=tree_messages,scalar_evidence=scalar_evidence,p256_recoveries=recovered_targets,all_checks=True,cases=21,completed=completed,partial=partial,cap_drops=drops,kernel_calls=len(calls),records_accounted=packets,counts=dict(send_per_record=14,send_fixed=2,append_accept=13,append_drop=7,home=16,clear_per_flag=3,complete_max=67,gather224=31,build8replies=84,requests_per_candidate=17,requests_fixed=3),max_tags_roundtrip=True,scope='Finite executed four-stage resolution transport, cap drops, fixed reply association and completeness checks; global asynchronous machine and full instruction allowances require explicit proof.')))

"""Transport-only harness; synthetic keys isolate communication from curve arithmetic."""
import json,random
from pathlib import Path
from types import SimpleNamespace
import a64,range_witness
ROOT=Path(__file__).resolve().parent
SENTINEL=0x5E5E5E5E5E5E5E5E
M35=(1<<35)-1
LINK=1<<40
rg=random.Random(314159)
CERT=SimpleNamespace(DLV=6)
def put(mem,ptr,v):
    for k in range(4):mem[ptr+8*k]=(v>>(64*k))&((1<<64)-1)
class Port:
    """48-byte port capture: [dest, addr, w0..w3]; message completes at +40."""
    def __init__(self):
        self.words, self.msgs = [0] * 6, []
    def store(self, k, v):
        self.words[k // 8] = v
        if k == 40:
            self.msgs.append(list(self.words))
class LinkView:
    """Message bytes at LINK over destination memory."""
    def __init__(self, mem, msg):
        self.mem, self.msg = mem, msg
    def __getitem__(self, a):
        return self.msg[(a - LINK) // 8] if LINK <= a < LINK + 48 else self.mem[a]
    def __setitem__(self, a, v):
        assert not LINK <= a < LINK + 48
        self.mem[a] = v
def link(prog, labels, mem, msg, calls, cert):
    """Deliver one message via route.s dlv; record count against cert.DLV."""
    calls.append(("dlv", a64.a64_run(prog, labels, "dlv", dict(x0=LINK, x30=SENTINEL), LinkView(mem, msg)), cert.DLV))
class Node:
    """One node of the 4x4 machine: its memory at uniform addresses."""
    P1, P2, HB, SB, SP = 0x10000, 0x20000, 0x100000, 1 << 34, 0x7000000
    R1, E1, R2, E2 = 0x1000000, 0x2000000, 0x3000000, 0x4000000
    CT, C4, BK, PA, CL = 0x5000000, 0x5100000, 1 << 36, 0x6100000, 0x6200000
    DC = 0x6300000
    STRIDE = 1 << 36                 # d3 bucket t's entries end at BK + (t + 1) * STRIDE
    def __init__(self):
        self.mem = {}
def run_k(prog, labels, entry, regs, nd, calls, formula, port=None):
    """Run one kernel; record (entry, count, formula). Returns final registers."""
    o = {}
    regs = dict(regs, sp=nd.SP, x30=SENTINEL)
    dev = (port[0], 48, port[1].store) if port else None
    calls.append((entry, a64.a64_run(prog, labels, entry, regs, nd.mem, out_regs=o, dev=dev), formula))
    return o
def range_scan_model(groups,cap):
  """Independent path-count model for paired CMP/CCMP scanning."""
  count=1;pairs=[];offset=0
  if not groups:return 2,pairs
  class Full(Exception):pass
  def append(i,j):
    nonlocal count
    count+=2
    if len(pairs)>=cap:raise Full
    pairs.append((offset+i,offset+j));count+=1
  try:
   for keys in groups:
    count+=3
    if not keys:count+=2;continue
    m=len(keys)-1;i=1;count+=5
    if m&1:
      count+=4
      if keys[0]==keys[1]:append(0,1);count+=1
      count+=2;i+=1
    r=(m//2)&31;count+=12+(r&1)
    def pair_step():
      nonlocal count,i
      count+=6
      first=keys[i-1]==keys[i];second=keys[i]==keys[i+1]
      if first or second:
        count+=2
        if first:append(i-1,i)
        count+=2
        if second:append(i,i+1);count+=1
      i+=2
    for _ in range(r):pair_step()
    count+=3
    while i<len(keys):
      for _ in range(32):pair_step()
      count+=3
    count+=2;offset+=len(keys)
  except Full:return count+1,pairs
  return count+1,pairs
def check():
    cert=CERT;rg=random.Random(32452843);prog,labels=a64.a64_load((ROOT/'route.s').read_text())
    bp,bl=a64.a64_load(range_witness.block_production(ROOT))
    side,q2,cap1,b3,b4,cap3,ncl=4,257,128,2,4,7,64
    n3,n4,st1,stride=1<<b3,1<<b4,32*cap1,32*cap3
    ids=[(i,j) for i in range(side) for j in range(side)];nodes={n:Node() for n in ids};calls=[]
    buckets={n:[[] for _ in range(n3)] for n in ids};drops=dict(hop1=0,bucket=0,filtered=0)
    msgs_ok=slots_ok=scan_ok=tags_ok=True;shared=0;acks=0;dones=0;completions=0
    for n,nd in nodes.items():
        run_k(prog,labels,'init3',dict(x0=n3,x1=nd.CT,x2=nd.BK,x3=stride),nd,calls,5*n3+1)
        run_k(prog,labels,'zero4',dict(x0=n4,x1=nd.C4),nd,calls,n4//4*4+1)
    def planted():
        d1,d2=rg.randrange(side),rg.randrange(side)
        top=(0xffffffff00000000|rg.getrandbits(32)) if rg.random()<0.2 else rg.getrandbits(63)
        return d1|(d2<<64)|(rg.getrandbits(64)<<128)|(top<<192)
    pool=[planted() for _ in range(10)]
    for phase in range(2):
        donecount={n:0 for n in ids}
        for n in ids:nodes[n].mem[Node.DC]=0
        ref1={n:{i:[] for i in range(side)} for n in ids}
        for (i,j),nd in nodes.items():
            run_k(bp,bl,'block_init',dict(x0=side,x1=nd.P1+64,x2=0x9000000,x3=0x8000020,x4=nd.R1+i*st1),nd,calls,775*side+2)
            batch=[]
            for l in range(q2):
                pair=[]
                for sign in range(2):
                    if rg.random()<0.4:key=rg.choice(pool);shared+=1
                    else:key=planted()
                    pair.append(key);put(nd.mem,nd.SB+64*l+32*sign,key)
                batch.append(pair)
                nd.mem[0x7100000+72*l]=l<<36
            W2=i|((phase*side*side+side*i+j)<<35);exp=[];cur={};full=filtered=0
            for l,pair in enumerate(batch):
                for sign,key in enumerate(pair):
                    d1=key&M35;d2=(key>>64)&M35;ptr=cur.get(d1,nd.R1+i*st1);cur[d1]=ptr+32
                    if ptr>=nd.R1+(i+1)*st1:full+=1;continue
                    exp.append([d1,ptr,(key>>128)&((1<<64)-1),key>>192,W2,d2|(sign<<35)|(l<<36)])
            port_words=[0]*4;block_msgs=[];ack_events=[]
            def block_store(offset,value):
                port_words[offset//8]=value
                if offset!=24:return
                d,addr,source,count=port_words
                assert count==128 and source==0x8000020+4096*d
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
                count=a64.a64_run(bp,bl,'block_dlv',dict(x0=LINK,x30=SENTINEL),BlockView())
                calls.append(('block_dlv',count,551));assert len(ack_events)==before+1 and ack_events[-1]==(d,addr)
                block_msgs.append((d,addr))
            emitter=a64.a64_run(bp,bl,'pair_keys',dict(x0=q2,x1=nd.SB,x2=nd.P1,x3=W2,x4=nd.R1+(i+1)*st1,x5=0x7100000,x6=0x9000000,x30=SENTINEL),nd.mem,dev=(nd.P1,32,block_store))
            total_full=sum((ptr-(nd.R1+i*st1))//4096 for ptr in cur.values())
            assert emitter<=101*q2//4+30+8*total_full
            calls.append(('block_pair_keys',emitter,emitter))
            drops['hop1']+=full;drops['filtered']+=filtered
            for m in exp:ref1[(m[0],j)][i].append(tuple(m[2:]))
            port=Port();end=nd.R1+(i+1)*st1
            def finish_store(offset,value):
                if offset<32:block_store(offset,value)
                elif offset>=64:port.store(offset-64,value)
            ct=a64.a64_run(bp,bl,'block_finish',dict(x0=side,x1=nd.P1+64,x2=0x9000000,x3=0x8000020,x4=0xB000000,x5=end,x7=0xB000040,x8=nd.E1+32*i,x30=SENTINEL),nd.mem,dev=(0xB000000,112,finish_store))
            assert ct<=40*side+4;calls.append(('block_finish',ct,ct))
            msgs_ok &= port.msgs==[[d,nd.E1+32*i,min(cur.get(d,nd.R1+i*st1),end),0,0,0] for d in range(side)]
            for m in port.msgs:link(prog,labels,nodes[(m[0],j)].mem,m,calls,cert)
        for (a,j),nd in nodes.items():
            expected=[]
            for i in range(side):
                got=[tuple(nd.mem[x+8*k] for k in range(4)) for x in range(nd.R1+i*st1,nd.mem[nd.E1+32*i],32)]
                slots_ok &= got==ref1[(a,j)][i]
                for rec in ref1[(a,j)][i]:
                    d2=rec[3]&M35;out=list(rec);out[3]=(out[3]&~M35)|j;expected.append([d2,0]+out)
            port=Port();nr=len(expected);loops=sum(len(ref1[(a,j)][i])//8 for i in range(side))
            run_k(prog,labels,'h2s',dict(x0=nd.P2,x3=side,x13=nd.E1,x14=nd.R1,x15=st1,x17=j),nd,calls,2+11*side+7*nr+2*loops,(nd.P2,port))
            msgs_ok &= port.msgs==expected
            dport=Port()
            run_k(prog,labels,'h2done',dict(x0=nd.P2,x2=side,x3=1),nd,calls,3+6*side,(nd.P2,dport))
            msgs_ok &= dport.msgs==[[d,1,0,0,0,0] for d in range(side)]
            for m in port.msgs:
                home=(a,m[0]);dst=nodes[home];rec=tuple(m[2:]);t=rec[0]&(n3-1);ok=len(buckets[home][t])<cap3
                addr=dst.BK+t*stride+32*len(buckets[home][t])
                if ok:buckets[home][t].append((addr,rec))
                else:drops['bucket']+=1
                out={}
                count=a64.a64_run(prog,labels,'receive',dict(x0=LINK,x1=n3-1,x2=dst.CT,x3=n4-1,x4=dst.C4,x13=dst.DC,x15=donecount[home],x16=side,x17=1,x18=32,x19=(1<<64)-32,x30=SENTINEL,sp=dst.SP),LinkView(dst.mem,m),out_regs=out)
                calls.append(('receive',count,16 if ok else 12));assert dst.mem[dst.DC]==0;acks+=1
                assert out['x30']==SENTINEL and out['sp']==dst.SP
            for m in dport.msgs:
                home=(a,m[0]);dst=nodes[home];out={}
                final=donecount[home]+1==side
                count=a64.a64_run(prog,labels,'receive',dict(x0=LINK,x1=n3-1,x2=dst.CT,x3=n4-1,x4=dst.C4,x13=dst.DC,x15=donecount[home],x16=side,x17=1,x18=32,x19=(1<<64)-32,x30=SENTINEL,sp=dst.SP),LinkView(dst.mem,m),out_regs=out)
                calls.append(('receive_done',count,10 if final else 8));dones+=1
                donecount[home]+=1
                assert out['x15']==(0 if final else donecount[home])
                assert dst.mem[dst.DC]==(1 if final else 0)
                if final:completions+=1
                assert out['x30']==SENTINEL and out['sp']==dst.SP
        for home in ids:assert donecount[home]==side and nodes[home].mem[Node.DC]==1
        for home,nd in nodes.items():
            cnt=[0]*n4
            for t,bk in enumerate(buckets[home]):
                slots_ok &= nd.mem[nd.CT+16*t]==nd.BK+t*stride+32*len(bk) and nd.mem[nd.CT+16*t+8]==nd.BK+(t+1)*stride
                for addr,rec in bk:
                    scan_ok &= tuple(nd.mem[addr+8*k] for k in range(4))==rec;cnt[rec[1]&(n4-1)]+=1
                    k=rec[2]>>35;tags_ok &= (k%side)==(rec[3]&M35) and (k//side)%side==(rec[2]&M35) and k//(side*side)<=phase and (rec[3]>>36)<q2
            slots_ok &= [nd.mem[nd.C4+8*i] for i in range(n4)]==cnt
    found=set();exp_c=set()
    for home,nd in nodes.items():
        bk=buckets[home];seq=[e for t in range(n3) for e in bk[t]]
        run_k(prog,labels,'prefix4',dict(x0=n4,x1=nd.C4,x2=nd.PA),nd,calls,n4//4*10+1)
        run_k(prog,labels,'pass4',dict(x0=n3,x1=nd.CT,x2=nd.BK,x3=stride,x4=nd.C4,x5=n4-1),nd,calls,11*n3+6*len(seq)+2*sum(len(v)//8 for v in bk)+1)
        order=sorted(seq,key=lambda e:e[1][1]&(n4-1));scan_ok &= [nd.mem[nd.PA+8*k] for k in range(len(seq))]==[e[0] for e in order]
        groups=[[e[1][0] for e in order if e[1][1]&(n4-1)==t] for t in range(n4)];mc,cands=range_scan_model(groups,ncl)
        o=run_k(prog,labels,'scan_ranges',dict(x0=n4,x1=nd.PA,x2=nd.C4,x3=nd.CL,x4=nd.CL+16*ncl),nd,calls,mc)
        scan_ok &= 64*mc<=195*len(seq)+64*(26*n4+8*(ncl+1)+2)
        for c in range(nd.CL,o['x3'],16):found.add((home,nd.mem[c],nd.mem[c+8]))
        exp_c |= {(home,order[a][0],order[b][0]) for a,b in cands}
    bb,bd=a64.a64_load((ROOT/'block.s').read_text());erows=[]
    for n,fl in [(1,0),(1,1),(7,1),(8,0),(8,1),(9,0),(9,1)]:
        m1={};SB=1<<34;P1=0x10000;W=0x8000000;RS=0x9000000;IDB=0x7100000;RE=RS+65536
        for j in range(n):
            d1=j%2;m1[SB+64*j]=d1;m1[SB+64*j+8]=j;m1[SB+64*j+16]=j+1;m1[SB+64*j+24]=j+2;m1[SB+64*j+32]=d1;m1[SB+64*j+40]=j+3;m1[SB+64*j+48]=j+4;m1[SB+64*j+56]=j+5;m1[IDB+72*j]=j
        for d1 in [0,1]:
            m1[P1+64+8*d1]=W-32 if fl else W+32;m1[RS+8*d1]=RE-8192
        m2=dict(m1);c1=[];c2=[]
        r1=a64.a64_run(bp,bl,'pair_keys',dict(x0=n,x1=SB,x2=P1,x3=17|(31<<35),x4=RE,x5=IDB,x6=RS,x30=SENTINEL),m1,dev=(P1,32,lambda k,v:c1.append((k,v))))
        r0=a64.a64_run(bb,bd,'pair_keys',dict(x0=n,x1=SB,x2=P1,x3=17|(31<<35),x4=RE,x5=IDB,x6=RS,x30=SENTINEL),m2,dev=(P1,32,lambda k,v:c2.append((k,v))))
        f=sum(1 for k,v in c1 if k==24)
        assert m1==m2 and c1==c2 and r1==r0+4-14*(n//8),(n,fl,r1,r0)
        assert r1<=101*n//4+30+8*f,(n,fl,r1,f)
        erows.append((n,fl,r1,r0,f))
    counts={}
    for e,c,f in calls:counts[e]=counts.get(e,0)+c
    return dict(nodes=len(ids),phases=2,leaves=2*2*q2*len(ids),shared_leaves=shared,candidates=len(found),cands_ok=found==exp_c and len(found)>0,msgs_ok=msgs_ok,slots_ok=slots_ok,tags_ok=tags_ok,scan_ok=scan_ok,drops=drops,calls=len(calls),counts=counts,data_receives=acks,done_markers=dones,completions=completions,signaled=completions==2*len(ids),formulas_ok=all(c==f for _,c,f in calls),emit_unroll=erows,scope='Synchronous finite two-phase replay; asynchronous protocol requires separate proof.')

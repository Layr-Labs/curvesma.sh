"""Adapted from the published six-leaf transport harness; tests transport only.
Synthetic canonical-sized keys isolate communication from curve arithmetic.
"""
import json,random
from pathlib import Path
from types import SimpleNamespace
import a64
ROOT=Path(__file__).resolve().parent
SENTINEL=0x5E5E5E5E5E5E5E5E
M35=(1<<35)-1
LINK=1<<40
rg=random.Random(314159)
CERT=SimpleNamespace(TINIT_ENTRY=3,INIT3_BKT=4,ZERO4_PER4=4,FIN_DEST=9,
 DLV=6,H2S_FIX=3,H2S_SLOT=11,H2S_REC=10,H2S_DROP=9,LOOP8=2,
 BKT_FIX=4,BKT_SLOT=12,BKT_REC=10,BKT_DROP=8,PREFIX4_PER4=10,PASS4_BKT=11,PASS4_REC=4)
def put(mem,ptr,v):
    for k in range(4):mem[ptr+8*k]=(v>>(64*k))&((1<<64)-1)


class Port:
    """Capture of a 48-byte port (proof.md 2.1): destination +0, address +8, four data words +16..+47; a message is
    complete when its last word, at +40, is written."""

    def __init__(self):
        """Empty capture. returns: None."""
        self.words, self.msgs = [0] * 6, []

    def store(self, k, v):
        """Device hook. k: byte offset in the block; v: word stored. returns: None."""
        self.words[k // 8] = v
        if k == 40:
            self.msgs.append(list(self.words))

class LinkView:
    """Address space of one delivery: the message's 48 bytes at LINK, the destination's memory everywhere else."""

    def __init__(self, mem, msg):
        """mem: the destination's memory; msg: [destination, address, w0, w1, w2, w3]. returns: None."""
        self.mem, self.msg = mem, msg

    def __getitem__(self, a):
        """A load. a: byte address. returns: the word."""
        return self.msg[(a - LINK) // 8] if LINK <= a < LINK + 48 else self.mem[a]

    def __setitem__(self, a, v):
        """A store, into the destination's memory. a: byte address; v: word. returns: None."""
        assert not LINK <= a < LINK + 48
        self.mem[a] = v

def link(prog, labels, mem, msg, calls, cert):
    """Deliver one message (proof.md 2.1) by running route.s dlv on it; record its count against cert.DLV.

    prog, labels: route.s; mem: destination memory; msg: [dest, address, w0..w3]; calls: list. returns: None."""
    calls.append(("dlv", a64.a64_run(prog, labels, "dlv", dict(x0=LINK, x30=SENTINEL), LinkView(mem, msg)), cert.DLV))

class Node:
    """One node of the 4 x 4 machine: its memory, at the same addresses on every node (proof.md 5)."""
    P1, P2, HB, SB, SP = 0x10000, 0x20000, 0x100000, 1 << 34, 0x7000000
    R1, E1, R2, E2 = 0x1000000, 0x2000000, 0x3000000, 0x4000000
    CT, C4, BK, PA, CL = 0x5000000, 0x5100000, 1 << 36, 0x6100000, 0x6200000
    STRIDE = 1 << 36                 # d3 bucket t's entries end at BK + (t + 1) * STRIDE

    def __init__(self):
        """Empty memory. returns: None."""
        self.mem = {}

def run_k(prog, labels, entry, regs, nd, calls, formula, port=None):
    """Run one kernel on a node (outermost RET through the sentinel x30) and record (entry, count, formula).

    prog, labels: program; entry: label; regs: registers; nd: Node; calls: list; formula: the certificate's count;
    port: (address, Port) or None. returns: final registers.
    """
    o = {}
    regs = dict(regs, sp=nd.SP, x30=SENTINEL)
    dev = (port[0], 48, port[1].store) if port else None
    calls.append((entry, a64.a64_run(prog, labels, entry, regs, nd.mem, out_regs=o, dev=dev), formula))
    return o

def range_scan_model(groups,cap):
  """Exact independent path-count model for paired CMP/CCMP scanning."""
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
    cert=CERT;rg=random.Random(32452843);prog,labels=a64.a64_load((ROOT/'route.s').read_text());lp,ll=a64.a64_load((ROOT/'pair_keys.s').read_text())
    side,q2,cap1,b3,b4,cap3,ncl=4,4,3,2,4,7,64
    n3,n4,st1,stride=1<<b3,1<<b4,32*cap1,32*cap3
    ids=[(i,j) for i in range(side) for j in range(side)];nodes={n:Node() for n in ids};calls=[]
    buckets={n:[[] for _ in range(n3)] for n in ids};drops=dict(hop1=0,bucket=0,filtered=0)
    msgs_ok=slots_ok=scan_ok=tags_ok=True;shared=0;acks=0
    for n,nd in nodes.items():
        run_k(prog,labels,'init3',dict(x0=n3,x1=nd.CT,x2=nd.BK,x3=stride),nd,calls,5*n3+1)
        run_k(prog,labels,'zero4',dict(x0=n4,x1=nd.C4),nd,calls,n4//4*4+1)
    def planted():
        d1,d2=rg.randrange(side),rg.randrange(side)
        top=(0xffffffff00000000|rg.getrandbits(32)) if rg.random()<0.2 else rg.getrandbits(63)
        return d1|(d2<<64)|(rg.getrandbits(64)<<128)|(top<<192)
    pool=[planted() for _ in range(10)]
    for phase in range(2):
        ref1={n:{i:[] for i in range(side)} for n in ids}
        for (i,j),nd in nodes.items():
            run_k(prog,labels,'tinit',dict(x0=side,x1=nd.P1+64,x2=nd.R1+i*st1),nd,calls,3*side+1)
            batch=[]
            for l in range(q2):
                pair=[]
                for sign in range(2):
                    if rg.random()<0.4:key=rg.choice(pool);shared+=1
                    else:key=planted()
                    pair.append(key);put(nd.mem,nd.SB+64*l+32*sign,key)
                batch.append(pair)
            W2=i|((phase*side*side+side*i+j)<<35);exp=[];cur={};full=filtered=0
            for l,pair in enumerate(batch):
                for sign,key in enumerate(pair):
                    d1=key&M35;d2=(key>>64)&M35;ptr=cur.get(d1,nd.R1+i*st1);cur[d1]=ptr+32
                    if ptr>=nd.R1+(i+1)*st1:full+=1;continue
                    exp.append([d1,ptr,(key>>128)&((1<<64)-1),key>>192,W2,d2|(sign<<35)|(l<<36)])
            port=Port();run_k(lp,ll,'pair_keys',dict(x0=q2,x1=nd.SB,x2=nd.P1,x3=W2,x4=nd.R1+(i+1)*st1),nd,calls,29*q2+11-4*full,(nd.P1,port))
            msgs_ok &= port.msgs==exp;drops['hop1']+=full;drops['filtered']+=filtered
            for m in port.msgs:link(prog,labels,nodes[(m[0],j)].mem,m,calls,cert)
            for m in exp:ref1[(m[0],j)][i].append(tuple(m[2:]))
            port=Port();end=nd.R1+(i+1)*st1
            run_k(prog,labels,'fin',dict(x0=nd.P1,x1=side,x2=nd.E1+32*i,x6=end),nd,calls,3+9*side,(nd.P1,port))
            msgs_ok &= port.msgs==[[d,nd.E1+32*i,min(cur.get(d,nd.R1+i*st1),end),0,0,0] for d in range(side)]
            for m in port.msgs:link(prog,labels,nodes[(m[0],j)].mem,m,calls,cert)
        # Hop1 fully delivered before reading slots. Each phase drains hop2 before reuse.
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
            for m in port.msgs:
                home=(a,m[0]);dst=nodes[home];rec=tuple(m[2:]);t=rec[0]&(n3-1);ok=len(buckets[home][t])<cap3
                addr=dst.BK+t*stride+32*len(buckets[home][t])
                if ok:buckets[home][t].append((addr,rec))
                else:drops['bucket']+=1
                ack=dst.P2+48;dst.mem[ack]=1;out={}
                count=a64.a64_run(prog,labels,'receive',dict(x0=LINK,x1=n3-1,x2=dst.CT,x3=n4-1,x4=dst.C4,x13=ack,x17=1,x18=32,x19=(1<<64)-32,x30=SENTINEL,sp=dst.SP),LinkView(dst.mem,m),out_regs=out)
                calls.append(('receive',count,15 if ok else 11));assert dst.mem[ack]==0;acks+=1
                assert out['x30']==SENTINEL and out['sp']==dst.SP
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
    counts={}
    for e,c,f in calls:counts[e]=counts.get(e,0)+c
    return dict(nodes=len(ids),phases=2,leaves=2*2*q2*len(ids),shared_key_leaves=shared,candidates=len(found),candidates_equal_exact_order_reference=found==exp_c and len(found)>0,messages_equal_emit_reference=msgs_ok,slots_equal_reference=slots_ok,source_tags_ok=tags_ok,buckets_pointer_order_and_scan_ok=scan_ok,drops=drops,kernel_calls=len(calls),counts=counts,acknowledged_records=acks,formulas_ok=all(c==f for _,c,f in calls),scope='Synchronous finite two-phase replay; asynchronous protocol requires separate proof.')

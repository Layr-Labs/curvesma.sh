"""Adapted from the published six-leaf transport harness; tests transport only.
Synthetic canonical-sized keys isolate communication from curve arithmetic.
"""
import json,random
from pathlib import Path
from types import SimpleNamespace
import a64
ROOT=Path(__file__).resolve().parent
SENTINEL=0x5E5E5E5E5E5E5E5E
M36=(1<<36)-1
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
    CT, C4, BK, PA, CL = 0x5000000, 0x5100000, 1 << 38, 0x6100000, 0x6200000
    STRIDE = 1 << 38                 # d3 bucket t's entries end at BK + (t + 1) * STRIDE

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
    """Exact path count for sixteen comparisons and branch-tree prefixes."""
    count=1;pairs=[];offset=0
    if not groups:return 2,pairs
    for keys in groups:
        count+=3
        if not keys:
            count+=2;continue
        k=(len(keys)-1)%16
        count+=11+(k&1);i=0
        def compare(current,j):
            if keys[j]==keys[j+1]:
                current+=2
                if len(pairs)>=cap:return current,True
                pairs.append((offset+j,offset+j+1));current+=2
            return current,False
        for _ in range(k):
            count,done=compare(count+4,i)
            if done:return count+1,pairs
            i+=1
        count+=4
        while i<len(keys)-1:
            for j in range(16):
                count,done=compare(count+3+(j%2==0),i)
                if done:return count+1,pairs
                i+=1
            count+=2
        count+=2;offset+=len(keys)
    return count+1,pairs

def check():
    cert=CERT
    """4 x 4 machine (proof.md 2.3, 4): keys, fin, h2s, fin, bkt, prefix4, pass4 and scan against an exact model of the
    same order, each message delivered by running dlv; every call's count against its formula. cert: certificate.
    returns: dict."""
    rg = random.Random(314159 + 32452843)   # own stream: the census is unchanged
    prog, labels = a64.a64_load(open(str(ROOT/"route.s")).read())
    lprog, llabels = a64.a64_load(open(str(ROOT/"pair_keys.s")).read())
    side, q2, cap1, cap2, b3, b4, cap3, ncl = 4, 4, 7, 7, 2, 4, 8, 64
    n3, n4, st1, st2, bkb = 2 ** b3, 2 ** b4, 32 * cap1, 32 * cap2, 16 * cap3
    ids = [(i, j) for i in range(side) for j in range(side)]
    nodes = {n: Node() for n in ids}
    calls = []
    for (i, j), nd in nodes.items():
        tf = cert.TINIT_ENTRY * side + 1
        run_k(prog, labels, "tinit", dict(x0=side, x1=nd.P1 + 64, x2=nd.R1 + i * st1), nd, calls, tf)
        run_k(prog, labels, "tinit", dict(x0=side, x1=nd.P2 + 64, x2=nd.R2 + j * st2), nd, calls, tf)
        run_k(prog, labels, "init3", dict(x0=n3, x1=nd.CT, x2=nd.BK + nd.STRIDE - bkb, x3=nd.STRIDE), nd, calls,
              cert.INIT3_BKT * n3 + 1)
        run_k(prog, labels, "zero4", dict(x0=n4, x1=nd.C4), nd, calls, n4 // 4 * cert.ZERO4_PER4 + 1)
    def planted():
        d1,d2=rg.randrange(side),rg.randrange(side)
        return d1 | (d2<<64) | (rg.getrandbits(64)<<128) | (rg.getrandbits(63)<<192)
    pool=[planted() for _ in range(10)]
    batches,shared={},0
    for n in ids:
        lst=[]
        for _ in range(2*q2):
            if rg.random()<0.25:
                lst.append(rg.choice(pool));shared+=1
            else:lst.append(planted())
        batches[n]=[lst[2*l:2*l+2] for l in range(q2)]
    # hop 1: keys' emits write each record into the pair slot of its sender row at (d1, j); then the final cursors
    msgs_ok, ref1, drops = True, {n: {i: [] for i in range(side)} for n in ids}, dict(hop1=0, hop2=0, bucket=0)
    fin_f = 3 + cert.FIN_DEST * side
    for (i, j), nd in nodes.items():
        for l,pair in enumerate(batches[(i,j)]):
            for sign,key in enumerate(pair):put(nd.mem,nd.SB+64*l+32*sign,key)
        W2=i | (side*i+j)<<36
        exp=[];cur={};full=0
        for l,pair in enumerate(batches[(i,j)]):
            for sign,key in enumerate(pair):
                d1=key&M36;d2=(key>>64)&M36
                ptr=cur.get(d1,nd.R1+i*st1)
                cur[d1]=ptr+32
                if ptr>=nd.R1+(i+1)*st1:full+=1;continue
                exp.append([d1,ptr,(key>>128)&((1<<64)-1),key>>192,W2,d2|(sign<<36)|(l<<37)])
        port=Port()
        o=run_k(lprog,llabels,'pair_keys',dict(x0=q2,x1=nd.SB,x2=nd.P1,x3=W2,x4=nd.R1+(i+1)*st1),nd,calls,29*q2+11-4*full,(nd.P1,port))
        msgs_ok &= port.msgs==exp
        drops["hop1"] += full
        for m in port.msgs:
            link(prog, labels, nodes[(m[0], j)].mem, m, calls, cert)
        for m in exp:
            ref1[(m[0], j)][i].append(tuple(m[2:]))
        port = Port()
        end1 = nd.R1 + (i + 1) * st1
        run_k(prog, labels, "fin", dict(x0=nd.P1, x1=side, x2=nd.E1 + 32 * i, x6=end1), nd, calls, fin_f, (nd.P1, port))
        msgs_ok &= port.msgs == [[d, nd.E1 + 32 * i, min(cur.get(d, nd.R1 + i * st1), end1), 0, 0, 0]
                                 for d in range(side)]
        for m in port.msgs:
            link(prog, labels, nodes[(m[0], j)].mem, m, calls, cert)
    # barrier; hop 2: h2s streams the pair slots through the row port to (own row, d2)
    ref2 = {n: {j: [] for j in range(side)} for n in ids}
    slots_ok = True
    for (a, j), nd in nodes.items():
        for i in range(side):
            got = [tuple(nd.mem[x + 8 * k] for k in range(4)) for x in range(nd.R1 + i * st1, nd.mem[nd.E1 + 32 * i], 32)]
            slots_ok &= got == ref1[(a, j)][i]
        used, stored, nrec = {}, 0, 0
        for i in range(side):
            for rec in ref1[(a, j)][i]:
                nrec += 1
                d2 = rec[3] & M36
                if used.get(d2, 0) < cap2:
                    used[d2] = used.get(d2, 0) + 1
                    ref2[(a, d2)][j].append(rec)
                    stored += 1
                else:
                    drops["hop2"] += 1
        port = Port()
        loops = sum(len(ref1[(a, j)][i]) // 8 for i in range(side))
        run_k(prog, labels, "h2s", dict(x0=nd.P2, x3=side, x13=nd.E1, x14=nd.R1, x15=st1, x12=nd.R2 + (j + 1) * st2),
              nd, calls,
              cert.H2S_FIX + cert.H2S_SLOT * side + cert.H2S_REC * stored + cert.H2S_DROP * (nrec - stored)
              + cert.LOOP8 * loops,
              (nd.P2, port))
        msgs_ok &= len(port.msgs) == stored
        for m in port.msgs:
            link(prog, labels, nodes[(a, m[0])].mem, m, calls, cert)
        port = Port()
        run_k(prog, labels, "fin", dict(x0=nd.P2, x1=side, x2=nd.E2 + 32 * j, x6=nd.R2 + (j + 1) * st2), nd, calls,
              fin_f, (nd.P2, port))
        for m in port.msgs:
            link(prog, labels, nodes[(a, m[0])].mem, m, calls, cert)
    # barrier; grouping at each home
    found, exp_c, scan_ok, tags_ok = set(), set(), True, True
    for n, nd in nodes.items():
        for j in range(side):
            got = [tuple(nd.mem[x + 8 * k] for k in range(4)) for x in range(nd.R2 + j * st2, nd.mem[nd.E2 + 32 * j], 32)]
            slots_ok &= got == ref2[n][j]
            tags_ok &= all(rec[2] >> 36 == side * (rec[2] & M36) + j and rec[3] & M36 == n[1] for rec in got)
        bk, nrec = [[] for _ in range(n3)], 0
        for j in range(side):
            for k, rec in enumerate(ref2[n][j]):
                nrec += 1
                t = rec[0] & (n3 - 1)
                if len(bk[t]) < cap3:
                    bk[t].append((nd.R2 + j * st2 + 32 * k, rec))
                else:
                    drops["bucket"] += 1
        seq = [e for t in range(n3) for e in bk[t]]
        loops = sum(len(ref2[n][j]) // 8 for j in range(side))
        run_k(prog, labels, "bkt", dict(x0=side, x13=nd.E2, x14=nd.R2, x15=st2, x3=n3 - 1, x4=nd.CT, x5=n4 - 1,
                                        x6=nd.C4), nd, calls,
              cert.BKT_FIX + cert.BKT_SLOT * side + cert.BKT_REC * len(seq) + cert.BKT_DROP * (nrec - len(seq))
              + cert.LOOP8 * loops)
        b0 = nd.BK + nd.STRIDE - bkb
        scan_ok &= all([nd.mem[b0 + t * nd.STRIDE + 16 * k], nd.mem[b0 + t * nd.STRIDE + 16 * k + 8]] ==
                       [e[0], e[1][1] & (n4 - 1)] for t in range(n3) for k, e in enumerate(bk[t]))
        run_k(prog, labels, "prefix4", dict(x0=n4, x1=nd.C4, x2=nd.PA), nd, calls, n4 // 4 * cert.PREFIX4_PER4 + 1)
        run_k(prog, labels, "pass4", dict(x0=n3, x1=nd.CT, x2=b0, x3=nd.STRIDE, x4=nd.C4), nd, calls,
              cert.PASS4_BKT * n3 + cert.PASS4_REC * len(seq) + cert.LOOP8 * sum(len(b) // 8 for b in bk) + 1)
        order = sorted(seq, key=lambda e: e[1][1] & (n4 - 1))
        scan_ok &= [nd.mem[nd.PA + 8 * k] for k in range(len(seq))] == [e[0] for e in order]
        groups=[[e[1][0] for e in order if e[1][1] & (n4-1)==t] for t in range(n4)]
        mc,cands=range_scan_model(groups,ncl)
        o=run_k(prog,labels,"scan_ranges",dict(x0=n4,x1=nd.PA,x2=nd.C4,x3=nd.CL,x4=nd.CL+16*ncl),nd,calls,mc)
        scan_ok &= mc <= (29*len(seq)+7)//8+24*n4+4*(ncl+1)+2
        for c in range(nd.CL,o["x3"],16):found.add((n,nd.mem[c],nd.mem[c+8]))
        exp_c |= {(n,order[a][0],order[b][0]) for a,b in cands}
    groups = {}
    for n in ids:
        for six in batches[n]:
            for x in six:
                groups.setdefault(x, []).append(n)
    cnt, form = {}, {}
    for e, c, f in calls:
        cnt[e], form[e] = cnt.get(e, 0) + c, form.get(e, 0) + f
    res = dict(nodes=len(ids), leaves=2 * q2 * len(ids), shared_key_leaves=shared,
               equal_key_groups=sum(len(v) >= 2 for v in groups.values()), candidates=len(found),
               candidates_equal_exact_order_reference=found == exp_c and len(found) > 0,
               messages_equal_emit_reference=msgs_ok, slots_equal_reference=slots_ok, source_tags_ok=tags_ok,
               buckets_pointer_order_and_scan_ok=scan_ok, drops=drops, kernel_calls=len(calls), counts=cnt,
               deliveries=sum(e == "dlv" for e, _, _ in calls), formulas_ok=all(c == f for _, c, f in calls))
    return res

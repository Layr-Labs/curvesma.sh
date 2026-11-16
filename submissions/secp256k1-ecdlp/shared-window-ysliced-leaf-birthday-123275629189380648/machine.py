"""The memory of every node of the declared run, and the scheduler. Plain text on purpose: this file is the whole
answer to "how does a word stored by one node reach another node" in the experiment.

Nodes. A node is one address space (class Node; roles Gen, Rel, Home, YRel) at the production addresses and field
widths of certificate.py. kernels/a64.py runs a kernel with the node as its memory: every load is Node.__getitem__,
every store Node.__setitem__. There is no other path between nodes: no message object, no port, no link model.

Shared words. A window, cell, list or frame area is memory of one node, its owner (World.mem[owner]). A reader maps
the same words at its own addresses, read-only (the rows of Node.SPEC with a wire). The rules, enforced on every
executed load and store (a broken rule raises Fault):
  1. one writer: a store through a read-only mapping raises, the cursor word included;
  2. a window is written only by instructions of its owner's kernels: a store by the harness into a window, or into
     an arrival or release cell, raises;
  3. a window is unusable until its owner has reset it (stored the empty cursor value);
  4. a reader may load a slot word only after it has loaded a cursor value that covers the slot, and only a word the
     owner stored since the window's last reset;
  5. once a reader has loaded a non-empty cursor, the owner must reset the window before it stores into it again;
  6. a reader's later kernel call that loads the same non-empty cursor again, with no reset between, raises;
  7. a private word must lie in a row of certificate.REGISTRY for the role, and a load of a word never stored raises
     (except ranges declared zero or unspecified, see Node.zero and Node.any);
  8. a load outside a kernel (Node.load: a word of the driver's step, or a step that has no kernel) reads the node's
     own memory only: of a word another node owns it raises.
Private set-up state (jump tables, sentinel words, control words, stacks) and test inputs are written by the harness
into private words; the certificate charges that set-up as node initialisation. Python reads memory freely (peek) to
compare it with the reference models, and never passes a value it read to another node, with one exception: before
a wave, leafx builds in the worker's memory the table rows that the wave's frames index, which it finds by peeking
the frames' id words (the production tables hold all rows; a row depends only on X and its index, and a row that was
not built makes the kernel's load raise). The steps that have no kernel (REQ_SERVE, RESOLVE, XSETUP) are done by
Python on words it loads from the acting node's own memory through Node.load, under the same rules: proof.md 5.1
says what each one loads, and that the point served in Phase 0 is a test input because the toy does not mine.
A release cell's payload words and, in tree 1, its coin line are plain shared words: rules 1 and 2 hold for them,
and a load before the owner's first store raises. Arrival and release words read zero until stored (Node.ZERO): that
zero is set-up state, charged in WIN_INIT, and no kernel stores it.

Waits. LDAXR arms a monitor on the 2-KiB granule of the loaded word; a store by any node into that granule sets the
event of every node that armed it; WFE with no event pending puts the node to sleep (Sched.sleep). This is premise MW
as the harness models it; nothing else wakes a node.

Scheduler. Every node of a run is a thread with its own program (a Python function that calls kernels). One node runs
at a time and keeps the processor until it sleeps in WFE or its program ends (class Sched); the next node is the
next in the given order that is awake or has its event set. Phases are therefore separated only by the executed
barrier kernels (function barrier): no Python code orders the phases of different nodes.

What the model does not show. Nodes run one at a time here and concurrently in the machine; no interleaving inside a
kernel call is modelled, and neither is the timing of a store against a concurrent load. For the ordering between
nodes the run offers the executed barrier kernels, the comparison of the final memory in two processor orders, and
the controls (a barrier released early, a sweep skipped); that a window is not read while it is written rests on
the barriers (proof.md 5.1), not on an observed interleaving.
"""
import contextlib, io, os, sys, threading
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path[:0] = [HERE, os.path.join(HERE, "kernels")]
import unpack, a64
_out = io.StringIO()
with contextlib.redirect_stdout(_out): import certificate as K
CERT_TEXT = _out.getvalue()
M, SNT, TC = 2**64 - 1, 0x5E5E5E5E5E5E5E5E, K.TCELL
_SRC = {}
def ktext(name):
    """@return: the text of kernels/<name>.s (restored by unpack.py if compressed; leaf: kernels/leafm.s goes first)"""
    f = "kernels/" + name + ".s"; t = unpack.TEXT.get(f) or open(os.path.join(HERE, f)).read()
    return unpack.TEXT.get("kernels/leafm.s", "") + t if name == "leaf" else t
def load(name, extra=""):
    """the assembled kernels of kernels/<name>.s (plus harness stubs), cached"""
    if (name, extra) not in _SRC: _SRC[name, extra] = a64.a64_load(ktext(name) + extra)
    return _SRC[name, extra]

class Fault(Exception):
    """a memory rule broken by an executed load or store"""

class Kind:
    """Geometry of one shared cell. size: bytes a cell; slot: bytes a slot (0: plain words, no rules); off: offset of
    the cursor word in the cell (the words below it are plain). Three cursor conventions: down (hop windows, request
    cells, Y out-windows: the cursor is the address of the last slot stored, empty = the last slot of the cell), up
    (first slot at cursor word + up; the cursor is the address after the last slot stored) and rel (lists: the
    cursor is the byte count of the entries, which start at cursor word + up)."""
    def __init__(s, size, slot, up=0, rel=False, off=0): s.size, s.slot, s.up, s.rel, s.off = size, slot, up, rel, off
    def empty(s, b): return 0 if s.rel else b + s.up if s.up else b + s.size - s.slot
    def covers(s, b, c, a):
        if s.rel: return b + s.up <= a < b + s.up + c
        return b + s.up <= a < c if s.up else c + s.slot <= a < b + s.size
FRAMES = K.YQ_BASE + 2**36 + 2**31 - (K.YQ_BASE - 64)
REC, CELL, YPK, YOK = Kind(K.WIN, 32), Kind(K.CELLH, 32), Kind(K.YPW, 512, 64), Kind(K.YIW, 512)
FRM = Kind(FRAMES, 512, 128)  # the frame area of a sub-ring: cursor word at YQ_BASE - 64, the data of frame f (448 B, a slot's offsets 64 .. 511) at YQ_BASE + 512 f + 64
RAW = Kind(TC, 0)  # a tree cell of the tier tree: plain words
ARR = Kind(TC, 64, 64, True, 2112)  # a tree cell of the global tree: barrier words, then a list of 15 entries
BKT = Kind(2**17, 64, 64, True)  # one bucket of a leader's bucket lists
DNL = Kind(K.TD_BYTES, 64, 64, True)  # a leader's down-list

class World:
    """the memory of all nodes (one dict of 64-bit words an owner) and the state of every window"""
    def __init__(s):
        s.mem, s.reset, s.seen, s.armed, s.nodes, s.cur, s.lazy, s.sched = {}, {}, {}, {}, {}, None, None, None
        s.big = {}  # generator key -> its tree children that hold bucket lists
        s.trees = [{}, {}]  # tree 0 (all nodes) and tree 1 (generators and Y-relays): key -> (parent, children)
    def rd(s, o, a):
        try: return s.mem[o][a]
        except KeyError: return s.nodes[o].unwritten(a)
    def plant_tree(s, t, parent):
        """parent: node key -> its parent's key (None for the root); wires the tree cells of every node of tree t"""
        T = s.trees[t] = {k: (p, []) for k, p in parent.items()}
        for k, p in parent.items():
            if p is not None: T[p][1].append(k)
        for k in parent: s.nodes[k].wire_tree(t)

class Node:
    """one node's address space; kernels/a64.py uses it as its memory"""
    ROLE, TR, ZERO, ANY = "", 0, (), ()
    def __init__(s, W, *key, **opt):
        s.W, s.key, s.ev, s.mon, s.asleep, s.sleeps = W, (type(s).__name__,) + key, 0, None, False, 0
        assert s.key not in W.nodes; W.nodes[s.key] = s; s.d = W.mem[s.key] = {}; s.pg = {}; s.epoch = 0
        s.calls, s.insns, s.step_loads, s.exempt = {}, 0, 0, ()
        s.REG = sorted((lo, hi) for r, _, lo, hi, _ in K.REGISTRY if r == s.ROLE)
        s.ZERO = [(s.TR + t * K.TREE_T + c * TC, s.TR + t * K.TREE_T + c * TC + 8) for t in (0, 1) for c in (0, 1)] + list(s.ZERO)  # its arrival and release words
        s.ANY = list(s.ANY); s.SPEC = list(s.spec(*key, **opt))
        for lo, hi, _, _, _ in s.SPEC: assert s.inreg(lo) and s.inreg(hi - 8)
    def inreg(s, a): return any(lo <= a < hi for lo, hi in s.REG)
    def zero(s, lo, hi): s.ZERO.append((lo, hi))  # words that read zero until stored (a fresh mapping, or zeroed by an executed kernel)
    def any(s, lo, hi): s.ANY.append((lo, hi))  # words a kernel may load before any store; they read as an arbitrary fixed value
    def unwritten(s, a):
        if any(lo <= a < hi for lo, hi in s.ZERO): return 0
        if any(lo <= a < hi for lo, hi in s.ANY): return (a * 0x9E3779B97F4A7C15 ^ 0x5DEECE66D) & M | 1
        raise Fault("load of a word never stored")
    def tree_base(s, t): return s.TR + t * K.TREE_T
    def wire_tree(s, t):
        """own arrival cell at B, own release cell at B + TC, the parent's release cell mapped at B + 2 TC, a private
        staging cell at B + 3 TC, the children's arrival cells mapped from B + 4 TC"""
        W, B, kd = s.W, s.tree_base(t), (ARR, RAW)[t]; par, kids = W.trees[t][s.key]
        s.SPEC.append((B, B + TC, TC, kd, None))
        if par is not None: s.SPEC.append((B + 2 * TC, B + 3 * TC, TC, RAW, lambda _: (par, W.nodes[par].tree_base(t) + TC)))
        if kids: s.SPEC.append((B + 4 * TC, B + (4 + len(kids)) * TC, TC, kd, lambda j: (kids[j], W.nodes[kids[j]].tree_base(t))))
        s.pg.clear()
    def find(s, a):
        """address -> (owner, the owner's address, kind or None for a plain word, the owner's address of the cursor word, read-only)"""
        for lo, hi, st, kd, wire in s.SPEC:
            if lo <= a < hi:
                off = (a - lo) % st
                if wire is None: o, b, ro = s.key, a - off, False
                else: (o, b), ro = wire((a - lo) // st), True
                if off < kd.off or not kd.slot: return o, b + off, None, 0, ro
                return o, b + off, kd, b + kd.off, ro
        return s.key, a, None, 0, False
    def private_page(s, a):
        """True (and remembered) when the 64-KiB page of a holds private registry words only"""
        lo, hi = a >> 16 << 16, (a >> 16) + 1 << 16
        if any(l < hi and lo < h for l, h, _, _, _ in s.SPEC) or not any(l <= lo and hi <= h for l, h in s.REG): return False
        s.pg[a >> 16] = None; return True
    def peek(s, a):
        """Python's own look at memory (no rule, no effect); 0 for a word never stored"""
        o, oa = s.find(a & M)[:2]; return s.W.mem[o].get(oa, 0)
    def __contains__(s, a): return a in s.d
    def __getitem__(s, a):
        if a >> 16 in s.pg:
            try: return s.d[a]
            except KeyError: return s.unwritten(a)
        o, oa, kd, b, ro = s.find(a); W = s.W; assert a % 8 == 0
        if kd is None:
            if not ro and not s.private_page(a) and not s.inreg(a): raise Fault("address outside the role's registry")
            return W.rd(o, oa)
        k = o, b
        if k not in W.reset:
            if W.lazy is None or ro: raise Fault("load in a window its owner never reset")
            W.lazy.add(b); W.reset[k] = {}; W.mem[o][b] = kd.empty(b)  # dry pass only: note the window
        if oa == b:
            v = W.rd(o, b)
            if ro and v != kd.empty(b):
                sn = W.seen.setdefault(k, {})
                if sn.get(s.key, (0, s.epoch))[1] != s.epoch: raise Fault("window swept again before its reset")
                sn[s.key] = v, s.epoch
            return v
        if ro:
            c = W.seen.get(k, {}).get(s.key); r = (oa - b - kd.up) % kd.slot
            if c is None or not kd.covers(b, c[0], oa) or not W.reset[k].get(oa - r, 0) >> (r >> 3) & 1: raise Fault("load past the cursor")
        return W.rd(o, oa)
    def __setitem__(s, a, v):
        W = s.W
        if a >> 16 in s.pg: s.d[a] = v
        else:
            o, oa, kd, b, ro = s.find(a); assert a % 8 == 0
            if ro: raise Fault("store through a read-only window")
            if kd is None:
                if not s.private_page(a) and not s.inreg(a): raise Fault("address outside the role's registry")
                if W.cur is not s and 0 <= a - s.TR < 2 * K.TREE_T and (a - s.TR) % K.TREE_T < 2 * TC: raise Fault("store into a window by the harness")
            else:
                k = o, b
                if W.cur is not s: raise Fault("store into a window by the harness")
                if oa == b and v == kd.empty(b): W.reset[k] = {}; W.seen.pop(k, None)  # the owner's reset
                elif W.seen.get(k): raise Fault("store into a swept window before its reset")
                elif k not in W.reset: raise Fault("store into a window never reset")
                elif oa != b: r = (oa - b - kd.up) % kd.slot; d = W.reset[k]; d[oa - r] = d.get(oa - r, 0) | 1 << (r >> 3)
            s.d[a] = v
        if W.armed:
            for n in W.armed.pop((s.key, a >> 11), ()): n.ev, n.mon = 1, None  # clears the readers' monitors
    def arm(s, k):
        A = s.W.armed
        if s.mon in A and s in A[s.mon]: A[s.mon].remove(s)
        s.mon = k; A.setdefault(k, []).append(s)
    def run(s, P, entry, ls=None, ev0=0, wfe=None, visits=None, limit=10**8, **r):
        """run one kernel on this node: (executed instructions, final registers)"""
        out, W = {}, s.W; W.cur = s; s.epoch += 1; s.calls[entry] = s.calls.get(entry, 0) + 1
        if entry in ("hreq", "hreqm") and not r["x9"] > r["x12"]: raise Fault("request kernel called on an empty list")
        try: n = a64.a64_run(P[0], P[1], entry, {k: v & M for k, v in dict(r, x30=SNT, sp=s.SP).items()}, s, limit, out, visits, ls, wfe, ev0)
        finally: W.cur = None
        assert out["sp"] == s.SP; s.insns += n; return n, out
    def load(s, a, n=1):
        """Loads by the node outside a kernel, under the same rules: one word of its phase step (allowance PH_STEP or
        TCALL; returns the word), or n words read by a step that has no kernel (XSETUP, REQ_SERVE, RESOLVE; a list)."""
        W = s.W; W.cur = s; s.epoch += 1; s.step_loads += n
        if any(s.find(a + 8 * j)[0] != s.key for j in range(n)): W.cur = None; raise Fault("load of another node's memory outside a kernel")
        try: return s[a] if n == 1 else [s[a + 8 * j] for j in range(n)]
        finally: W.cur = None

class Gen(Node):
    """generator (g1, g2): g1 = (c mod 2^15, m), g2 = c >> 15 (as a worker: the slice of its sub-ring). Owns its hop-1
    windows W1, its producer windows YP, the frame area YQ (used by member 0) and its tree lists TU and TD; reads the
    R -> G request cells of its relays (QG), the out-windows of its Y-relays (YI), the frames of member 0 of its
    sub-ring through the alias YA, its tree children's bucket lists (TUI) and its tree parent's down-list (TDI)."""
    ROLE, TR, SP = "G", K.TRG, K.T34 + 0x7000000
    def spec(s, g1, g2):
        s.g1, s.g2, s.k, W = g1, g2, 0, s.W
        m0 = ("Gen", g1 & -1024, g2)
        def kid(j): return W.big[s.key][j >> 14], K.TUB + (j & 16383) * 2**17
        return ((0, K.SH, K.WIN, REC, None),
                (K.QG_BASE, K.QG_BASE + K.F1 * K.CELLG, K.CELLG, CELL, lambda v: (("Rel", g1, v), K.RQGO + g2 * K.CELLG)),
                (K.YP_BASE, K.YP_BASE + 2**14 * K.YPW, K.YPW, YPK, None),
                (K.YI_BASE, K.YI_BASE + 2**17 * K.YIW, K.YIW, YOK, lambda p1: (("YRel", g1 >> 10, p1, g2 >> 27), K.UP_BASE + (g2 & 2**27 - 1) * K.YIW)),
                (K.YQ_BASE - 64, K.YQ_BASE - 64 + FRAMES, FRAMES, FRM, None),
                (K.YQ_BASE - 64 - K.TB, K.YQ_BASE - 64 - K.TB + FRAMES, FRAMES, FRM, lambda _: (m0, K.YQ_BASE - 64)),
                (K.TUB, K.TUB + K.TU_BYTES, 2**17, BKT, None),
                (K.TUIB, K.TUIB + K.TUI_BYTES, 2**17, BKT, kid),
                (K.TDB, K.TDB + K.TD_BYTES, K.TD_BYTES, DNL, None),
                (K.TDIB, K.TDIB + K.TD_BYTES, K.TD_BYTES, DNL, lambda _: (W.trees[0][s.key][0], K.TDB)))
class Rel(Node):
    """relay (g1, v): owns its hop-2 windows W2 and its R -> G request cells; reads the hop-1 windows of relay digit v
    of the generators (g1, *) and the H -> R request cells of its homes (v, *). In production it maps the in-windows
    of all 2^41 generators and its sentinel word follows them at RSN; a toy relay (lo, n given) maps the n in-windows
    from generator index lo and keeps its sentinel word (s.sn) after them."""
    ROLE, TR, SP = "R", K.TRR, K.RCSB + 2**20
    def spec(s, g1, v, lo=0, n=2**K.G2B):
        s.g1, s.v, s.sn = g1, v, K.RW1I + (lo + n) * K.WIN
        return ((0, 2**55, K.WIN, REC, None),
                (K.RW1I + lo * K.WIN, s.sn, K.WIN, REC, lambda j: (("Gen", g1, lo + j), v * K.WIN)),
                (K.RQHI, K.RQGO, K.CELLH, CELL, lambda e: (("Home", v, e), K.HQO + g1 * K.CELLH)),
                (K.RQGO, K.RCSB, K.CELLG, CELL, None))
class Home(Node):
    """home (v, e): owns its H -> R request cells and its private table, R2 and lists; reads the hop-2 windows of home
    digit e of the relays (*, v)"""
    ROLE, TR, SP = "H", K.TRH, K.HCCB + 2**20
    def spec(s, v, e):
        s.v, s.e = v, e
        return ((K.HWI, K.HSN, K.WIN, REC, lambda g1: (("Rel", g1, v), e * K.WIN)),
                (K.HQO, K.HQO + 2**K.G1B * K.CELLH, K.CELLH, CELL, None))
class YRel(Node):
    """Y-relay (c15, p1, s1): owns its out-windows (one a sub-ring slice s0); reads producer window s1 of the producers
    (c15, m = 0; p1, *)"""
    ROLE, TR, SP = "Y", K.UP_BASE + 2**54, K.UP_BASE + 2**53 + 2**20
    def spec(s, c15, p1, s1):
        return ((0, 2**24 * K.YPW, K.YPW, YPK, lambda p0: (("Gen", c15 << 10, p1 << 24 | p0), K.YP_BASE + s1 * K.YPW)),
                (K.UP_BASE, K.UP_BASE + 2**27 * K.YIW, K.YIW, YOK, None))

def v1b(vs, us):
    """The record of a class label under V1b. vs, us: the 128-bit words v* and u* (H = high half, L = low half).
    Returns the relay digit L(v*) >> MLV (40 bits), the home digit L(u*) >> MLU (42 bits),
    K2 = H(v*)[63:42] : L(u*)[63:22] and K3 = H(u*)[63:21] : H(v*)[41:21]."""
    hv, lv, hu, lu = vs >> 64, vs & M, us >> 64, us & M; a, b = K.MH_A, K.MH_B
    return lv >> K.MLV, lu >> K.MLU, hv >> a + b << a + b | lu >> K.MLU, hu >> a << a | hv >> b & 2**a - 1

class Sched:
    """one processor at a time: a node runs until it sleeps in WFE or its program ends; then the next node in the
    order that is awake, or asleep with its event set, takes over"""
    def __init__(s): s.live, s.err, s.sem = [], None, {}
    def give(s, me):
        L = s.live; k = L.index(me) if me in L else -1
        for n in L[k + 1:] + L[:k + 1]:
            if not n.asleep or n.ev: s.sem[n].release(); return
        if L: s.err = s.err or Fault("every node sleeps without an event")
        for x in s.sem.values(): x.release()
    def body(s, me, f):
        s.sem[me].acquire()
        try:
            if not s.err: f(me)
        except Exception as x: s.err = s.err or x
        s.live.remove(me); s.give(me)
    def sleep(s, me):
        me.asleep = True
        while not me.ev:
            s.give(me); s.sem[me].acquire()
            if s.err: raise Fault("stopped")
        me.asleep = False
    def run(s, jobs):
        """jobs: [(node, program)] in processor order; returns when every program has ended"""
        s.live = [n for n, _ in jobs]; s.sem = {n: threading.Semaphore(0) for n in s.live}
        ts = [threading.Thread(target=s.body, args=j) for j in jobs]
        for t in ts: t.start()
        s.sem[s.live[0]].release()
        for t in ts: t.join()
        if s.err: raise s.err

STEP_KERNELS = frozenset("wreset rsw hsw rswe hswe rqsw gqsw yrst yprod yrel ym0 hreq hreqm".split())
BAR = dict(wfe=0, stale=0, insns=0, n=0)  # totals over all barrier steps (the driver compares them with the certificate)
def wait(n, entry, base, cell, **r):
    """a barrier kernel with waits: its count is base + WAKE a WFE executed; WFEs = sleeps + a pending event consumed"""
    S = n.W.sched
    def cb(mem, mon, ev):
        o, oa = n.find(mon[0] << 11)[:2]; n.arm((o, oa >> 11)); n.sleeps += 1; n.W.cur = None; S.sleep(n); n.W.cur = n; n.ev = 0; ev[0] = 1
    e0, n.ev, n.sleeps = n.ev, 0, 0
    c, o = n.run(load("route"), entry, ev0=e0, wfe=cb, **r); wf = (c - base) // K.WAKE
    assert c == base + K.WAKE * wf and 0 <= wf - n.sleeps <= e0, (entry, c, base, wf, n.sleeps, e0)
    if e0 and wf == n.sleeps: n.ev = 1  # not consumed: still pending
    o_, oa = n.find(cell)[:2]; n.arm((o_, oa >> 11))  # the last LDAXR left the monitor on this word
    BAR["wfe"] += wf; BAR["stale"] += wf - n.sleeps; BAR["insns"] += c; return o
def barrier(n, t, b, up=(0, 0), root=None, work_up=None, work_down=None, early=False):
    """Node n's part of barrier b on tree t. bpar: wait until every child has stored arrival word b, OR their payloads
    into up; work_up (the tree-list merge of a pass up) runs here; bleaf: store the payload and the arrival word b,
    wait for the parent's release word b, load the payload that comes down (the root decides it: root(up)); work_down
    (the list copy of a pass down) runs here; brel: store the payload and the release word b for the children.
    Returns (payload up, payload down). The kernel-call counters of the node (the bound behind STEP_FIX) restart."""
    B = n.tree_base(t); par, kids = n.W.trees[t][n.key]; nk = len(kids); BAR["n"] += 1
    over = [k for k, v in n.calls.items() if v > 2 and k in STEP_KERNELS and k not in n.exempt]
    if over: raise Fault("more than two calls of %s between barriers" % over[0])
    if nk: o = wait(n, "bpar", 1 + K.BPAR_CHILD0 * nk, B + (3 + nk) * TC, x2=b, x3=up[0], x4=up[1], x6=B + 4 * TC, x7=B + (4 + nk) * TC); up = (o["x3"], o["x4"])
    if work_up: work_up(n)
    if early: n.run(load("route"), "brel", x0=B + TC, x2=b, x3=up[0], x4=up[1])  # negative control only: the children are released before the node's own release
    if par is not None: o = wait(n, "bleaf", K.BLEAF0, B + 2 * TC, x0=B, x1=B + 2 * TC, x2=b, x3=up[0], x4=up[1]); down = (o["x3"], o["x4"])
    else: down = root(up) if root else up
    if work_down: work_down(n, down)
    if nk: c, _ = n.run(load("route"), "brel", x0=B + TC, x2=b, x3=down[0], x4=down[1]); BAR["insns"] += c; assert c == K.BREL
    n.calls = {}; return up, down

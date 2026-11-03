"""Part B of experiment.py (proof.md 2.1, 4): ports, an exact sender model, the 4 x 4 machine (deliveries by dlv or
dlvk; keys on fractions planted from the inverse DFT) and forced and unrolled paths. Helpers come from experiment.py (E)."""
import os
import random

K, BLK = 256, 1 << 14            # records per block message; record blocks 2^14 bytes apart (proof.md 2.1)
PB1, PB2 = 1 << 49, 1 << 50      # the column and row ports' record blocks, and their message blocks (proof.md 5)
MB1, MB2 = PB1 - 64, PB1 - 128   # == the production message-block addresses (proof.md 5, certificate.py MB1C/MB2C)
RC = 1 << 43                     # ORR with 2^43 turns an LC entry's address into its RC entry's
STALE = 0x5A5A5A5A5A5A5A5A       # a slot no record has been stored into
LINK = (1 << 64) - (1 << 14)     # where dlv and dlvk read the message in transit: an address of the link, not of a node
M36 = 2 ** 36 - 1


def dest(word):
    """The destination index a port takes from a destination word, bits 3..38 (proof.md 2.1). returns: integer."""
    return word >> 3 & M36


class Port:
    """Device hook of a port (proof.md 2.1): a one-record message completes when +40 of block mb is written, block
    d's (at pb + 2^14 d) when its header's second word is."""

    def __init__(self, mem, mb, pb):
        """mem: the sender's memory; mb: message block; pb: block 0. returns: None."""
        self.mem, self.mb, self.pb, self.words, self.msgs, self.blocks = mem, mb, pb, [0] * 6, [], []
        self.lo = min(mb, pb)
        self.dev = (self.lo, pb + (BLK << 36) - self.lo, self.store)

    def store(self, k, v):
        """A stored word. k: offset from lo; v: word. returns: None."""
        a = self.lo + k
        if self.mb <= a < self.mb + 48:
            self.words[(a - self.mb) // 8] = v
            if a == self.mb + 40:
                self.msgs.append(list(self.words))
        elif a >= self.pb and (a - self.pb) % BLK == 8:
            self.blocks.append([self.mem[a - 8], v, 0, 0] + [self.mem.get(a + 24 + 8 * w, STALE) for w in range(4 * K)])


class LinkView:
    """Address space of one delivery: the message's words at LINK, the destination's memory elsewhere."""

    def __init__(self, mem, msg, writes):
        """mem: destination memory; msg: its words; writes: list of store addresses. returns: None."""
        self.mem, self.msg, self.writes = mem, msg, writes

    def __getitem__(self, a):
        """A load. a: byte address. returns: the word."""
        return self.msg[(a - LINK) // 8] if LINK <= a < LINK + 8 * len(self.msg) else self.mem[a]

    def __setitem__(self, a, v):
        """A store, into the destination's memory. a: byte address; v: word. returns: None."""
        assert not LINK <= a < LINK + 8 * len(self.msg)
        self.writes.append(a)
        self.mem[a] = v


def deliver(E, prog, labels, mem, msg, calls, cert, writes=None):
    """Deliver a message by running dlv or dlvk on it, its count against cert.DLV or DLVK (-1 unless x0 returns to
    LINK, as 2.1's loop needs). mem: destination memory; msg: words; writes: list or None of store addresses.
    returns: None."""
    e, o = "dlv" if len(msg) == 6 else "dlvk", {}
    n = E.a64.a64_run(prog, labels, e, dict(x0=LINK, x30=E.SENTINEL), LinkView(mem, msg, [] if writes is None else writes),
                      out_regs=o)
    calls.append((e, n if o["x0"] == LINK else -1, cert.DLV if e == "dlv" else cert.DLVK))


class Sender:
    """Exact sender model (proof.md 4): per destination its block's records (n), slots and RC cursor; after slot k the
    cursor advances by 32k and the block is sent if the old cursor is below the end, else dropped (hop 2: end stored)."""

    def __init__(self, rc0, end, hop2, cert):
        """rc0: d -> initial cursor; end: slot end; hop2: whether a drop stores the end back. returns: None."""
        self.rc0, self.end, self.hop2, self.cert = rc0, end, hop2, cert
        self.n, self.rc, self.slots, self.msgs, self.seq, self.sent, self.drops, self.extra = {}, {}, {}, [], [], {}, 0, 0

    def cursor(self, d):
        """d's RC cursor. returns: integer."""
        return self.rc.get(d, self.rc0(d))

    def put(self, d, word, rec):
        """One record to destination d. word: its destination word; rec: four words. returns: None."""
        s, n = self.slots.setdefault(d, [None] * K), self.n.get(d, 0)
        s[n] = rec
        self.seq.append((d, rec))
        self.n[d] = (n + 1) % K
        if n + 1 < K:
            return
        c = self.cursor(d)
        self.rc[d] = c + 32 * K
        if c < self.end:
            self.msgs.append([word, c, 0, 0] + [w for r in s for w in r])
            self.sent[d] = self.sent.get(d, 0) + 1
            self.extra += self.cert.FLUSH if self.hop2 else self.cert.FLUSH_KEYS
        else:
            self.drops += 1
            self.extra += self.cert.H2S_BLOCK_DROP if self.hop2 else self.cert.KEY_BLOCK_DROP
            if self.hop2:
                self.rc[d] = self.end

    def final(self, d):
        """fin's view of destination d: (its final block's message or None, its final cursor). returns: tuple."""
        n, c = self.n.get(d, 0), self.cursor(d)
        msg = [8 * d, c, 0, 0] + [w for r in self.slots[d] for w in (r or [STALE] * 4)] if n and c < self.end else None
        return msg, min(c + 32 * n, self.end)

    def fin_count(self, ds):
        """route.s fin's count over the destinations ds. returns: integer."""
        f = self.cert.FIN_FIX
        for d in ds:
            n = self.n.get(d, 0)
            f += self.cert.FIN_DEST + (0 if not n else self.cert.FIN_SENT if self.cursor(d) < self.end else
                                       self.cert.FIN_HELD)
        return f


def keys_reference(E, mem, sb, q2, W2, snd, cert):
    """keys' records fed to snd in keys' order (entries downwards; sign +, then -; t = 0, 1, 2) by the bit-exact
    model of proof.md 3.2, 3.3 on the stored fractions; w3 adds the entry pointer (entry + 64 for +, entry - 96 for -)
    shifted left by 31. mem: memory; sb: entries; q2: count; W2: index word; snd: Sender. returns: (the count keys adds
    to x18, each record's (j, d1, d2, K2, K3))."""
    meter, out = 0, []
    for l in range(q2 - 1, -1, -1):
        e = sb + 256 * l
        for s, ptr in ((0, e + 64), (1, e - 96)):
            f6 = E.fracs(mem, e, s)
            for t in range(3):
                f = E.record_fields(E.leaf_w(f6, t))
                meter += cert.JPATH if f[0] else 0
                snd.put(f[1], 8 * f[1], [f[3], f[4], W2, f[2] + ((ptr << 31) & E.M64)])
                out.append(f)
    return meter, out


def plant(E, mem, sb, sixes):
    """Store T(c), T(beta c) of the components c whose leaves are the given ones (3.1's inverse DFT, no linear terms),
    gamma's V1 in its s-form rounded down to 0 mod 8. sixes: per incidence its leaves in keys' order. returns: None."""
    for l, six in enumerate(sixes):
        for s in (0, 1):
            v = [E.t128(E.BT[m] * c % E.p) for c in E.inverse_dft(six[3 * s:3 * s + 3]) for m in (0, 1)]
            v[5] = v[5] >> 64 << 64 | (v[5] & E.M64) >> 3
            for o, w in zip(E.FRAC_OFF[s], v):
                E.put(mem, sb + 256 * l + o, w, 2)


def state_ok(mem, snd, lct, pb):
    """Whether mem's LC and RC entries and slots equal snd's. lct: LC table; pb: block 0. returns: bool."""
    for d, s in snd.slots.items():
        b, e = pb + BLK * d, lct + 8 * d
        if mem[e] != b + 32 + 32 * snd.n.get(d, 0) or mem[e | RC] != snd.cursor(d):
            return False
        if any([mem[b + 32 * t + 32 + 8 * w] for w in range(4)] != s[t] for t in range(K) if s[t] is not None):
            return False
    return True


def rule_ok(dm, writes, snd, finals):
    """Whether each destination holds the one-record rule's records (record t of d at rc0(d) + 32t if below the end;
    with finals None, those of blocks sent), each final cursor is min(rc0 + 32 records, end) and every delivery store
    lay in [rc0(d), end). dm, writes: {d: memory}, {d: addresses}; finals: {d: cursor} or None. returns: bool."""
    ok, cnt = True, {}
    for d, rec in snd.seq:
        t = cnt[d] = cnt.get(d, -1) + 1
        c = snd.rc0(d) + 32 * t
        if c < snd.end and (finals is not None or t < K * snd.sent.get(d, 0)):
            ok &= [dm[d].get(c + 8 * w) for w in range(4)] == rec
    if finals is not None:
        ok &= all(f == min(snd.rc0(d) + 32 * (cnt.get(d, -1) + 1), snd.end) for d, f in finals.items())
    return ok and all(snd.rc0(d) <= a < snd.end for d, w in writes.items() for a in w)


def range_scan_model(groups, cap):
    """Path model of scan_ranges (proof.md 4): per range the first record, a comparison when n - 1 is odd, the TBZ
    tree into the 32-pair block, the block per 64 records; a pair costs 6, a hit re-tests both. groups: each range's
    K2s; cap: list entries. returns: (count, candidate index pairs)."""
    count, pairs, offset = 1, [], 0                  # cbz
    if not groups:
        return 2, pairs

    class Full(Exception):
        """The candidate list was full at an append."""

    def append(i, j):
        """One append after its list-end test (CMP, B.HS, STP). i, j: indices in the range. returns: None."""
        nonlocal count
        count += 2
        if len(pairs) >= cap:
            raise Full
        pairs.append((offset + i, offset + j))
        count += 1
    try:
        for keys in groups:
            count += 3                               # ldr end, cmp, b.hs
            if not keys:
                count += 2                           # subs, b.ne
                continue
            m, i = len(keys) - 1, 1
            count += 5                               # first pointer and K2, sub, lsr, tbz
            if m & 1:
                count += 4                           # ldr, ldr, cmp, b.eq
                if keys[0] == keys[1]:
                    append(0, 1)
                    count += 1                       # b back
                count += 2                           # mov, mov
                i += 1
            r = (m // 2) & 31
            count += 12 + (r & 1)                    # lsr, and, five movs, the TBZ tree (and a B when r is odd)

            def pair_step():
                """The pair of comparisons (i - 1, i) and (i, i + 1). returns: None."""
                nonlocal count, i
                count += 6
                first, second = keys[i - 1] == keys[i], keys[i] == keys[i + 1]
                if first or second:
                    count += 2                       # cmp, b.ne
                    if first:
                        append(i - 1, i)
                    count += 2                       # cmp, b.ne
                    if second:
                        append(i, i + 1)
                        count += 1                   # b back
                i += 2
            for _ in range(r):
                pair_step()
            count += 3                               # mov, cmp, b.lo
            while i < len(keys):
                for _ in range(32):
                    pair_step()
                count += 3
            count += 2                               # subs, b.ne
            offset += len(keys)
    except Full:
        return count + 1, pairs                      # ret
    return count + 1, pairs


class Node:
    """A node's memory; LC, RC, d3 and d4 tables and ports at proof.md 5's addresses."""
    SB, SP = 1 << 34, 0x7000000
    R1, E1, R2, E2 = 0x1000000, 0x2000000, 0x3000000, 0x4000000
    T1, T2, CT, C4, BK, PA, CL = 0, 1 << 39, 1 << 40, 1 << 51, 1 << 45, 0x6100000, 0x6200000
    STRIDE = 1 << 38                 # d3 bucket t's entries end at BK + (t + 1) * STRIDE

    def __init__(self):
        """Empty memory. returns: None."""
        self.mem = {}


def run_k(E, prog, labels, entry, regs, nd, calls, formula, port=None):
    """Run a kernel on nd (RET through the sentinel), recording (entry, count, formula). regs: registers; nd: Node;
    formula: the certificate's count; port: Port or None. returns: final registers."""
    o = {}
    regs = dict(regs, sp=nd.SP, x30=E.SENTINEL)
    calls.append((entry, E.a64.a64_run(prog, labels, entry, regs, nd.mem, out_regs=o, dev=port and port.dev), formula))
    return o


def hop(E, rt, kern, regs, formula, nd, snd, t, pb, mb, e, ds, calls, cert, to, wr=None):
    """A sender's hop: keys or h2s, then fin over ds, through port (mb, pb) with LC table t; messages checked against
    snd (fed) and delivered to to(d). kern: (prog, labels, entry); regs, formula: registers and count less flushes; e:
    final-cursor entry; wr: {d: list} of block stores or None. returns: (all equal the model, blocks, registers)."""
    port, fport = Port(nd.mem, mb, pb), Port(nd.mem, mb, pb)
    o = run_k(E, kern[0], kern[1], kern[2], regs, nd, calls, formula + snd.extra, port)
    ok = port.blocks == snd.msgs and not port.msgs and state_ok(nd.mem, snd, t, pb)
    run_k(E, rt[0], rt[1], "fin", dict(x0=mb, x7=t, x3=t | RC, x8=pb + 32, x1=len(ds), x2=e, x6=snd.end), nd, calls,
          snd.fin_count(ds), fport)
    fins = [snd.final(d) for d in ds]
    ok &= fport.blocks == [m for m, _ in fins if m] and fport.msgs == [[8 * d, e, c, 0, 0, 0] for d, (_, c) in zip(ds, fins)]
    for m in port.blocks + fport.blocks + fport.msgs:
        d = dest(m[0])
        deliver(E, rt[0], rt[1], to(d), m, calls, cert, wr[d] if wr is not None and len(m) > 6 else None)
    return ok, len(port.blocks) + len(fport.blocks), o


def part_b(E, cert):
    """4 x 4 machine (proof.md 2.3, 4): keys, fin, h2s, fin, bkt, prefix4, pass4 and scan_ranges against exact
    models, deliveries by dlv or dlvk, each count against its formula; slots of k records. returns: dict."""
    rg = random.Random(E.SEED + 32452843)   # own stream
    rt = E.a64.a64_load(open(os.path.join(E.HERE, "kernels", "route.s")).read())
    lk = E.a64.a64_load(open(os.path.join(E.HERE, "kernels", "leaf.s")).read(), Q2=4) + ("keys",)
    side, q2, b3, b4, cap3, ncl = 4, 4, 2, 4, 8, 64
    n3, n4, st, bkb = 2 ** b3, 2 ** b4, 32 * K, 16 * cap3
    ids = [(i, j) for i in range(side) for j in range(side)]
    nodes = {n: Node() for n in ids}
    calls = []
    for (i, j), nd in nodes.items():
        for t, pb, base in ((nd.T1, PB1, nd.R1 + i * st), (nd.T2, PB2, nd.R2 + j * st)):
            run_k(E, rt[0], rt[1], "init3", dict(x0=side, x1=t, x2=pb + 32, x3=BLK), nd, calls, cert.INIT3_BKT * side + 1)
            run_k(E, rt[0], rt[1], "tinit", dict(x0=side, x1=t | RC, x2=base), nd, calls, cert.TINIT_ENTRY * side + 1)
        run_k(E, rt[0], rt[1], "init3", dict(x0=n3, x1=nd.CT, x2=nd.BK + nd.STRIDE - bkb, x3=nd.STRIDE), nd, calls,
              cert.INIT3_BKT * n3 + 1)
        run_k(E, rt[0], rt[1], "zero4", dict(x0=n4, x1=nd.C4), nd, calls, n4 // 4 * cert.ZERO4_PER4 + 1)
    pool = [E.planted_x(rg, rg.randrange(side), rg.randrange(side))[1] for _ in range(10)]
    batches, shared = {}, 0
    for n in ids:
        lst = []
        for _ in range(6 * q2):
            if rg.random() < 0.15:
                lst.append(E.planted_x(rg, 0, 0, rg.choice(pool))[0])
                shared += 1
            else:
                lst.append(E.planted_x(rg, rg.randrange(side), rg.randrange(side))[0])
        batches[n] = [lst[6 * l:6 * l + 6] for l in range(q2)]
    # hop 1: keys' emits fill each destination's record block; fin sends the final blocks and cursors
    msgs_ok, ref1, blocks, lab_ok = True, {n: {i: [] for i in range(side)} for n in ids}, 0, True
    for (i, j), nd in nodes.items():
        plant(E, nd.mem, nd.SB, batches[(i, j)])
        W2 = i | (side * i + j) << 36
        snd = Sender(lambda d, b=nd.R1 + i * st: b, nd.R1 + (i + 1) * st, False, cert)
        meter, fl = keys_reference(E, nd.mem, nd.SB, q2, W2, snd, cert)
        lab_ok &= [f[1:] for f in fl] == [E.class_label(x) for six in batches[(i, j)][::-1] for x in six]
        ok, nb, o = hop(E, rt, lk, dict(x3=nd.SB + 256 * (q2 - 1) + 160, x12=snd.end, x17=W2, x18=0),
                        cert.KEY_FIX + cert.KEY_INC * q2 + meter, nd, snd, nd.T1, PB1, MB1, nd.E1 + 32 * i,
                        range(side), calls, cert, lambda d, j=j: nodes[(d, j)].mem)
        msgs_ok &= ok and o["x18"] == meter
        blocks += nb
        for d, rec in snd.seq:
            ref1[(d, j)][i].append(tuple(rec))
    # barrier; hop 2: h2s streams the pair slots through the row port to (own row, d2); fin
    ref2 = {n: {j: [] for j in range(side)} for n in ids}
    slots_ok = True
    for (a, j), nd in nodes.items():
        for i in range(side):
            got = [tuple(nd.mem[x + 8 * k] for k in range(4)) for x in range(nd.R1 + i * st, nd.mem[nd.E1 + 32 * i], 32)]
            slots_ok &= got == ref1[(a, j)][i]
        snd = Sender(lambda d, b=nd.R2 + j * st: b, nd.R2 + (j + 1) * st, True, cert)
        for i in range(side):
            for rec in ref1[(a, j)][i]:
                snd.put(rec[3] & M36, nd.T2 + 8 * (rec[3] & M36), list(rec))
                ref2[(a, rec[3] & M36)][j].append(rec)
        loops = sum(len(ref1[(a, j)][i]) // 8 for i in range(side))
        ok, nb, _ = hop(E, rt, rt + ("h2s",), dict(x3=side, x13=nd.E1, x14=nd.R1, x15=st, x12=snd.end),
                        cert.H2S_FIX + cert.H2S_SLOT * side + cert.H2S_REC * len(snd.seq) + cert.LOOP8 * loops, nd, snd,
                        nd.T2, PB2, MB2, nd.E2 + 32 * j, range(side), calls, cert, lambda d, a=a: nodes[(a, d)].mem)
        msgs_ok &= ok
        blocks += nb
    # barrier; grouping at each home
    found, exp_c, scan_ok, tags_ok, drops = set(), set(), True, True, 0
    for n, nd in nodes.items():
        for j in range(side):
            got = [tuple(nd.mem[x + 8 * k] for k in range(4)) for x in range(nd.R2 + j * st, nd.mem[nd.E2 + 32 * j], 32)]
            slots_ok &= got == ref2[n][j]
            tags_ok &= all(rec[2] >> 36 == side * (rec[2] & M36) + j and rec[3] & M36 == n[1] for rec in got)
        bk, nrec = [[] for _ in range(n3)], 0
        for j in range(side):
            for k, rec in enumerate(ref2[n][j]):
                nrec += 1
                t = rec[0] % 2 ** 24
                assert t < n3 and rec[1] % 2 ** 46 < n4
                if len(bk[t]) < cap3:
                    bk[t].append((nd.R2 + j * st + 32 * k, rec))
                else:
                    drops += 1
        seq = [e for t in range(n3) for e in bk[t]]
        loops = sum(len(ref2[n][j]) // 8 for j in range(side))
        run_k(E, rt[0], rt[1], "bkt", dict(x0=side, x13=nd.E2, x14=nd.R2, x15=st, x9=nd.CT, x11=nd.C4), nd, calls,
              cert.BKT_FIX + cert.BKT_SLOT * side + cert.BKT_REC * len(seq) + cert.BKT_DROP * (nrec - len(seq))
              + cert.LOOP8 * loops)
        b0 = nd.BK + nd.STRIDE - bkb
        scan_ok &= all([nd.mem[b0 + t * nd.STRIDE + 16 * k], nd.mem[b0 + t * nd.STRIDE + 16 * k + 8]] ==
                       [e[0], nd.C4 + 8 * (e[1][1] % 2 ** 46)] for t in range(n3) for k, e in enumerate(bk[t]))
        run_k(E, rt[0], rt[1], "prefix4", dict(x0=n4, x1=nd.C4, x2=nd.PA), nd, calls, n4 // 4 * cert.PREFIX4_PER4 + 1)
        run_k(E, rt[0], rt[1], "pass4", dict(x0=n3, x1=nd.CT, x2=b0, x3=nd.STRIDE), nd, calls,
              cert.PASS4_BKT * n3 + cert.PASS4_REC * len(seq) + cert.LOOP8 * sum(len(b) // 8 for b in bk) + cert.PASS4_FIX)
        order = sorted(seq, key=lambda e: e[1][1] % 2 ** 46)
        scan_ok &= [nd.mem[nd.PA + 8 * k] for k in range(len(seq))] == [e[0] for e in order]
        groups = [[e[1][0] for e in order if e[1][1] % 2 ** 46 == t] for t in range(n4)]
        mc, cands = range_scan_model(groups, ncl)
        o = run_k(E, rt[0], rt[1], "scan_ranges", dict(x0=n4, x1=nd.PA, x2=nd.C4, x3=nd.CL, x4=nd.CL + 16 * ncl), nd,
                  calls, mc)
        scan_ok &= mc <= cert.SCAN_REC * len(seq) + cert.SCAN_RANGE * n4 + 8 * (ncl + 1) + 2
        for c in range(nd.CL, o["x3"], 16):
            found.add((n, nd.mem[c], nd.mem[c + 8]))
        exp_c |= {(n, order[a][0], order[b][0]) for a, b in cands}
    groups = {}
    for n in ids:
        for six in batches[n]:
            for x in six:
                groups.setdefault(E.exact_key(x)[0], []).append(n)
    cnt, form = {}, {}
    for e, c, f in calls:
        cnt[e], form[e] = cnt.get(e, 0) + c, form.get(e, 0) + f
    res = dict(nodes=len(ids), leaves=6 * q2 * len(ids), shared_key_leaves=shared,
               equal_key_groups=sum(len(v) >= 2 for v in groups.values()), records_carry_class_labels=lab_ok,
               candidates=len(found),
               candidates_equal_exact_order_reference=found == exp_c and len(found) > 0,
               messages_equal_sender_model=msgs_ok, slots_equal_reference=slots_ok, source_tags_ok=tags_ok,
               buckets_pointer_order_and_scan_ok=scan_ok, bucket_drops=drops, kernel_calls=len(calls), counts=cnt,
               block_messages=blocks, deliveries=sum(e in ("dlv", "dlvk") for e, _, _ in calls),
               formulas_ok=all(c == f for _, c, f in calls) and lab_ok)
    res["block_paths"] = block_paths(E, rt, lk, rg, cert)
    res["unrolled_paths"] = unrolled_paths(E, rt[0], rt[1], rg, cert)
    return res


def block_paths(E, rt, lk, rg, cert):
    """Record blocks on every path (proof.md 8): keys and fin, and h2s, as described there. E: experiment module; rt,
    lk: route.s and leaf.s's keys; rg: Random; cert: certificate. returns: dict of checks."""
    calls, out = [], {}
    nd = Node()
    xs = [E.planted_x(rg, d, rg.randrange(4), j=d % 3)[0] for d in range(5)]   # destination d on KEY's path d mod 3
    ds = [d for d, n in enumerate((782, 268, 256, 2, 0)) for _ in range(n)]
    rg.shuffle(ds)
    q2, m, end = len(ds) // 6, (1, 2, 1, 1, 3), nd.R1 + 3 * 32 * K
    lk = E.a64.a64_load(open(os.path.join(E.HERE, "kernels", "leaf.s")).read(), Q2=q2) + ("keys",)
    leaves = [[xs[d] for d in ds[6 * l:6 * l + 6]] for l in range(q2)]
    plant(E, nd.mem, nd.SB, leaves)
    run_k(E, rt[0], rt[1], "init3", dict(x0=5, x1=nd.T1, x2=PB1 + 32, x3=BLK), nd, calls, cert.INIT3_BKT * 5 + 1)
    run_k(E, rt[0], rt[1], "tinit", dict(x0=5, x1=nd.T1 | RC, x2=end), nd, calls, cert.TINIT_ENTRY * 5 + 1)
    for d in range(5):
        nd.mem[(nd.T1 + 8 * d) | RC] = end - 32 * K * m[d]
    snd = Sender(lambda d: end - 32 * K * m[d], end, False, cert)
    meter, fl = keys_reference(E, nd.mem, nd.SB, q2, 7, snd, cert)
    dm, wr = {d: {} for d in range(5)}, {d: [] for d in range(5)}
    ok, nb, o = hop(E, rt, lk, dict(x3=nd.SB + 256 * (q2 - 1) + 160, x12=end, x17=7, x18=0),
                    cert.KEY_FIX + cert.KEY_INC * q2 + meter, nd, snd, nd.T1, PB1, MB1, nd.E1, range(5), calls, cert,
                    lambda d: dm[d], wr)
    ok &= o["x18"] == meter and rule_ok(dm, wr, snd, {d: dm[d][nd.E1] for d in range(5)})
    ok &= [f[:2] for f in fl] == [(d % 3, d) for l in range(q2 - 1, -1, -1) for d in ds[6 * l:6 * l + 6]]
    held = sum(1 for d in range(5) if snd.n.get(d, 0) and snd.cursor(d) >= end)
    out["keys_fin"] = dict(records=len(ds), blocks_sent=sum(snd.sent.values()), blocks_dropped=snd.drops,
                           finals_sent=nb - sum(snd.sent.values()), finals_held=held,
                           cursor_past_end=snd.cursor(0) > end, model_and_rule_ok=ok)
    nd = Node()
    D2 = [1 << 35 | rg.getrandbits(35) for _ in range(4)]
    sizes = list(range(20)) + [160] * 4 + [64]
    rg.shuffle(sizes)
    dl = [D2[t] for t, n in enumerate((556, 268, 50, 20)) for _ in range(n)]
    rg.shuffle(dl)
    st, end, m = 32 * 160, nd.R2 + 2 * 32 * K, {D2[0]: 1, D2[1]: 2, D2[2]: 1, D2[3]: 1}
    snd = Sender(lambda d: end - 32 * K * m[d], end, True, cert)
    for d in D2:
        nd.mem[nd.T2 + 8 * d], nd.mem[(nd.T2 + 8 * d) | RC] = PB2 + BLK * d + 32, end - 32 * K * m[d]
    k = 0
    for s_, n in enumerate(sizes):
        for t in range(n):
            rec = [rg.getrandbits(64), rg.getrandbits(64), rg.getrandbits(64), rg.getrandbits(28) << 36 | dl[k]]
            for w in range(4):
                nd.mem[nd.R1 + st * s_ + 32 * t + 8 * w] = rec[w]
            snd.put(dl[k], nd.T2 + 8 * dl[k], rec)
            k += 1
        nd.mem[nd.E1 + 32 * s_] = nd.R1 + st * s_ + 32 * n
    port = Port(nd.mem, MB2, PB2)
    run_k(E, rt[0], rt[1], "h2s", dict(x3=len(sizes), x13=nd.E1, x14=nd.R1, x15=st, x12=end), nd, calls,
          cert.H2S_FIX + cert.H2S_SLOT * len(sizes) + cert.H2S_REC * k + cert.LOOP8 * sum(n // 8 for n in sizes)
          + snd.extra, port)
    ok2 = port.blocks == snd.msgs and not port.msgs and state_ok(nd.mem, snd, nd.T2, PB2)
    dm, wr = {d: {} for d in D2}, {d: [] for d in D2}
    for msg in port.blocks:
        deliver(E, rt[0], rt[1], dm[dest(msg[0])], msg, calls, cert, wr[dest(msg[0])])
    ok2 &= rule_ok(dm, wr, snd, None)
    out["h2s"] = dict(records=k, blocks_sent=len(port.blocks), blocks_dropped=snd.drops,
                      end_stored_back=snd.cursor(D2[0]) == end, model_and_rule_ok=ok2)
    kf, h = out["keys_fin"], out["h2s"]
    out["counts_equal_formulas"] = all(c == f for _, c, f in calls)
    out["ok"] = (out["counts_equal_formulas"] and ok and ok2 and (kf["blocks_sent"], kf["blocks_dropped"],
                 kf["finals_sent"], kf["finals_held"]) == (3, 2, 2, 1) and kf["cursor_past_end"]
                 and (h["blocks_sent"], h["blocks_dropped"]) == (2, 1) and h["end_stored_back"])
    return out


def unrolled_paths(E, prog, labels, rg, cert):
    """bkt and pass4 on slots and buckets of 0 to 19 records, some full, at full width (proof.md 8), against exact
    models and formulas. prog, labels: route.s. returns: dict."""
    calls, out, ns = [], {}, list(range(20))
    rg.shuffle(ns)
    st = 32 * 19
    nd = Node()                     # bkt: 20 slots; d3 bucket t holds at most 24 + 12t entries
    n3, n4 = 4, 4
    T3, V4 = [1 << 23 | rg.getrandbits(23) for _ in range(n3)], [1 << 45 | rg.getrandbits(45) for _ in range(n4)]
    run_k(E, prog, labels, "init3", dict(x0=n3, x1=nd.CT, x2=nd.BK + nd.STRIDE - 16 * 64, x3=nd.STRIDE), nd, calls,
          cert.INIT3_BKT * n3 + 1)
    start = [nd.BK + (t + 1) * nd.STRIDE - 16 * (24 + 12 * t) for t in range(n3)]
    for t in range(n3):
        nd.mem[nd.CT + 8 * T3[t]] = start[t]
    run_k(E, prog, labels, "zero4", dict(x0=n4, x1=nd.C4), nd, calls, cert.ZERO4_PER4 + 1)
    for v in V4:
        nd.mem[nd.C4 + 8 * v] = 0
    bk, cnt4, stored = [[] for _ in range(n3)], [0] * n4, 0
    for s_, n in enumerate(ns):
        for k in range(n):
            t, v = rg.randrange(n3), rg.randrange(n4)
            rec = [rg.getrandbits(40) << 24 | T3[t], rg.getrandbits(16) << 48 | V4[v], rg.getrandbits(64), rg.getrandbits(64)]
            a = nd.R2 + st * s_ + 32 * k
            for w in range(4):
                nd.mem[a + 8 * w] = rec[w]
            if len(bk[t]) < 24 + 12 * t:
                bk[t].append((a, nd.C4 + 8 * V4[v]))
                cnt4[v] += 1
                stored += 1
        nd.mem[nd.E2 + 32 * s_] = nd.R2 + st * s_ + 32 * n
    run_k(E, prog, labels, "bkt", dict(x0=len(ns), x13=nd.E2, x14=nd.R2, x15=st, x9=nd.CT, x11=nd.C4), nd, calls,
          cert.BKT_FIX + cert.BKT_SLOT * len(ns) + cert.BKT_REC * stored + cert.BKT_DROP * (sum(ns) - stored)
          + cert.LOOP8 * sum(n // 8 for n in ns))
    ok_b = all([nd.mem[start[t] + 16 * k], nd.mem[start[t] + 16 * k + 8]] == list(e)
               for t in range(n3) for k, e in enumerate(bk[t]))
    ok_b &= all(nd.mem[nd.CT + 8 * T3[t]] == start[t] + 16 * len(bk[t]) for t in range(n3))
    ok_b &= [nd.mem[nd.C4 + 8 * v] for v in V4] == cnt4
    out["bkt"] = dict(records=sum(ns), bucketed=stored, entries_cursors_and_counts_ok=ok_b)
    nd = Node()                     # pass4: 20 buckets of 0 to 19 entries
    ents, b0 = [], nd.BK + nd.STRIDE - 16 * 19
    for t, n in enumerate(ns):
        cur = [(rg.getrandbits(64), rg.randrange(n4)) for _ in range(n)]
        for k, (a, v) in enumerate(cur):
            nd.mem[b0 + nd.STRIDE * t + 16 * k], nd.mem[b0 + nd.STRIDE * t + 16 * k + 8] = a, nd.C4 + 8 * v
        nd.mem[nd.CT + 8 * t] = b0 + nd.STRIDE * t + 16 * n
        ents.append(cur)
    tot = [sum(v == w for c in ents for _, v in c) for w in range(n4)]
    for w in range(n4):
        nd.mem[nd.C4 + 8 * w] = nd.PA + 8 * sum(tot[:w])
    run_k(E, prog, labels, "pass4", dict(x0=len(ns), x1=nd.CT, x2=b0, x3=nd.STRIDE), nd, calls,
          cert.PASS4_BKT * len(ns) + cert.PASS4_REC * sum(ns) + cert.LOOP8 * sum(n // 8 for n in ns) + cert.PASS4_FIX)
    order = [a for w in range(n4) for c in ents for a, v in c if v == w]
    out["pass4"] = dict(entries=sum(ns), pointer_order_ok=[nd.mem[nd.PA + 8 * k] for k in range(sum(ns))] == order)
    out["counts_equal_formulas"] = all(c == f for _, c, f in calls)
    out["scan_ranges"] = scan_paths(E, prog, labels, rg, cert)
    out["ok"] = (out["counts_equal_formulas"] and out["bkt"]["entries_cursors_and_counts_ok"]
                 and out["pass4"]["pointer_order_ok"] and out["bkt"]["bucketed"] < sum(ns) and out["scan_ranges"]["ok"])
    return out


def scan_paths(E, prog, labels, rg, cert):
    """scan_ranges alone on proof.md 8's ranges: candidates and counts against range_scan_model and the bound, and
    every instruction runs. prog, labels: route.s. returns: dict."""
    first = labels["scan_ranges"][0]
    visits = {i: 0 for i in range(first, len(prog))}
    runs = []

    def run(groups, cap):
        """One call on the given ranges (K2 lists) with a list of cap entries. returns: None."""
        nd = Node()
        a, cursor = nd.R2, nd.PA
        for t, g in enumerate(groups):
            for k2 in g:
                nd.mem[a], nd.mem[a + 8] = k2, rg.getrandbits(64)     # K2, and a K3 the scan never reads
                nd.mem[cursor] = a
                a, cursor = a + 32, cursor + 8
            nd.mem[nd.C4 + 8 * t] = cursor
        o = {}
        n = E.a64.a64_run(prog, labels, "scan_ranges", dict(x0=len(groups), x1=nd.PA, x2=nd.C4, x3=nd.CL,
                          x4=nd.CL + 16 * cap, sp=nd.SP, x30=E.SENTINEL), nd.mem, out_regs=o, visits=visits)
        mc, pairs = range_scan_model(groups, cap)
        got = [(nd.mem[c], nd.mem[c + 8]) for c in range(nd.CL, o["x3"], 16)]
        want = [(nd.R2 + 32 * i, nd.R2 + 32 * j) for i, j in pairs]
        bound = cert.SCAN_REC * sum(map(len, groups)) + cert.SCAN_RANGE * len(groups) + 8 * (len(pairs) + 1) + 2
        runs.append((n == mc and got == want and n <= bound, len(pairs), n, sum(map(len, groups))))
    ns = list(range(130))
    rg.shuffle(ns)
    run([[rg.randrange(3) for _ in range(n)] for n in ns], 4096)
    run([[rg.getrandbits(64) for _ in range(n)] for n in ns], 4096)
    run([[5, 5, 7, 7, 9]], 1)                          # the second equal pair, in the next pair, finds the list full
    run([[1, 2, 2, 4, 4] + list(range(10, 70))], 1)    # 65 records: the second equal pair, in the block, ends it
    run([], 1)
    run([[0] * n for n in range(0, 130, 3)], 2 ** 14)  # every range all equal: both comparisons of every pair hit
    hit = sum(v > 0 for v in visits.values())
    return dict(calls=len(runs), records=sum(r[3] for r in runs), listed=[r[1] for r in runs], counts=[r[2] for r in runs],
                instructions=len(visits), instructions_executed=hit,
                ok=all(r[0] for r in runs) and hit == len(visits) and runs[2][1] == 1 and runs[3][1] == 1)

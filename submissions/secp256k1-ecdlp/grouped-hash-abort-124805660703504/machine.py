import os
import random

K, BLK = 256, 1 << 14
PB1, PB2 = 1 << 50, 1 << 51
MB1, MB2 = (1 << 48) - 64, (1 << 48) - 128
D4W = 46
RC = 1 << 43
STALE = 0x5A5A5A5A5A5A5A5A
LINK = (1 << 64) - (1 << 14)
M36 = 2 ** 36 - 1


def dest(word):
    return word >> 3 & M36


class Port:

    def __init__(self, mem, mb, pb):
        self.mem, self.mb, self.pb, self.words, self.msgs, self.blocks = mem, mb, pb, [0] * 6, [], []
        self.lo = min(mb, pb)
        self.dev = (self.lo, pb + (BLK << 36) - self.lo, self.store)

    def store(self, k, v):
        a = self.lo + k
        if self.mb <= a < self.mb + 48:
            self.words[(a - self.mb) // 8] = v
            if a == self.mb + 40:
                self.msgs.append(list(self.words))
        elif a >= self.pb and (a - self.pb) % BLK == 8:
            self.blocks.append([self.mem[a - 8], v, 0, 0] + [self.mem.get(a + 24 + 8 * w, STALE) for w in range(4 * K)])


class LinkView:

    def __init__(self, mem, msg, writes):
        self.mem, self.msg, self.writes = mem, msg, writes

    def __getitem__(self, a):
        return self.msg[(a - LINK) // 8] if LINK <= a < LINK + 8 * len(self.msg) else self.mem[a]

    def __setitem__(self, a, v):
        assert not LINK <= a < LINK + 8 * len(self.msg)
        self.writes.append(a)
        self.mem[a] = v


def deliver(E, prog, labels, mem, msg, calls, cert, writes=None):
    e, o = "dlv" if len(msg) == 6 else "dlvk", {}
    n = E.a64.a64_run(prog, labels, e, dict(x0=LINK, x30=E.SENTINEL), LinkView(mem, msg, [] if writes is None else writes),
                      out_regs=o)
    calls.append((e, n if o["x0"] == LINK else -1, cert.DLV if e == "dlv" else cert.DLVK))


class Sender:

    def __init__(self, rc0, end, hop2, cert):
        self.rc0, self.end, self.hop2, self.cert = rc0, end, hop2, cert
        self.n, self.rc, self.slots, self.msgs, self.seq, self.sent, self.drops, self.extra = {}, {}, {}, [], [], {}, 0, 0

    def cursor(self, d):
        return self.rc.get(d, self.rc0(d))

    def put(self, d, word, rec):
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
            self.extra += self.cert.FLUSH
        else:
            self.drops += 1
            self.extra += self.cert.H2S_BLOCK_DROP if self.hop2 else self.cert.KEY_BLOCK_DROP
            if self.hop2:
                self.rc[d] = self.end

    def final(self, d):
        n, c = self.n.get(d, 0), self.cursor(d)
        msg = [8 * d, c, 0, 0] + [w for r in self.slots[d] for w in (r or [STALE] * 4)] if n and c < self.end else None
        return msg, min(c + 32 * n, self.end)

    def fin_count(self, ds):
        f = self.cert.FIN_FIX
        for d in ds:
            n = self.n.get(d, 0)
            f += self.cert.FIN_DEST + (0 if not n else self.cert.FIN_SENT if self.cursor(d) < self.end else
                                       self.cert.FIN_HELD)
        return f


def keys_reference(E, mem, sb, kb, q2, W2, snd, cert):
    meter, out = 0, []
    cb = (E.get(mem, kb + 64, 2), E.get(mem, kb + 80, 2))
    for l in range(q2 - 1, -1, -1):
        e = sb + 256 * l
        ca = (E.get(mem, e + 192, 2), E.get(mem, e + 208, 2))
        for s, ptr in ((0, e + 64), (1, e - 96)):
            c3 = [E.get(mem, e + o) for o in E.COMP_OFF[3 * s:3 * s + 3]]
            for t in range(3):
                f = E.record_fields(E.leaf_w(c3, t, ca, cb))
                meter += cert.JPATH if f[0] else 0
                snd.put(f[1], 8 * f[1], [f[3], f[4], W2, f[2] + ((ptr << 31) & E.M64)])
                out.append(f)
    return meter, out


def plant(E, mem, sb, sixes):
    for l, six in enumerate(sixes):
        e = sb + 256 * l
        for s in (0, 1):
            for o, c in zip(E.COMP_OFF[3 * s:3 * s + 3], E.inverse_dft(six[3 * s:3 * s + 3])):
                E.put(mem, e + o, c)
        E.put(mem, e + 192, 0)


def state_ok(mem, snd, lct, pb):
    for d, s in snd.slots.items():
        b, e = pb + BLK * d, lct + 8 * d
        if mem[e] != b + 32 + 32 * snd.n.get(d, 0) or mem[e | RC] != snd.cursor(d):
            return False
        if any([mem[b + 32 * t + 32 + 8 * w] for w in range(4)] != s[t] for t in range(K) if s[t] is not None):
            return False
    return True


def rule_ok(dm, writes, snd, finals):
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
    count, pairs, offset = 1, [], 0
    if not groups:
        return 2, pairs

    class Full(Exception):
        """The candidate list was full at an append."""

    def append(i, j):
        nonlocal count
        count += 2
        if len(pairs) >= cap:
            raise Full
        pairs.append((offset + i, offset + j))
        count += 1
    try:
        for keys in groups:
            count += 3
            if not keys:
                count += 2
                continue
            m, i = len(keys) - 1, 1
            count += 5
            if m & 1:
                count += 4
                if keys[0] == keys[1]:
                    append(0, 1)
                    count += 1
                count += 2
                i += 1
            r = (m // 2) & 31
            count += 12 + (r & 1)

            def pair_step():
                nonlocal count, i
                count += 6
                first, second = keys[i - 1] == keys[i], keys[i] == keys[i + 1]
                if first or second:
                    count += 2
                    if first:
                        append(i - 1, i)
                    count += 2
                    if second:
                        append(i, i + 1)
                        count += 1
                i += 2
            for _ in range(r):
                pair_step()
            count += 3
            while i < len(keys):
                for _ in range(32):
                    pair_step()
                count += 3
            count += 2
            offset += len(keys)
    except Full:
        return count + 1, pairs
    return count + 1, pairs


class Node:
    HB, SB, SP = 0x100000, 1 << 34, 0x7000000
    R1, E1, R2, E2 = 0x1000000, 0x2000000, 0x3000000, 0x4000000
    T1, T2, CT, C4, BK, PA, CL = 0, 1 << 39, 1 << 40, 1 << 52, 1 << 45, 0x6100000, 0x6200000
    STRIDE = 1 << 38

    def __init__(self):
        self.mem = {}


def run_k(E, prog, labels, entry, regs, nd, calls, formula, port=None):
    o = {}
    regs = dict(regs, sp=nd.SP, x30=E.SENTINEL)
    calls.append((entry, E.a64.a64_run(prog, labels, entry, regs, nd.mem, out_regs=o, dev=port and port.dev), formula))
    return o


def hop(E, rt, kern, regs, formula, nd, snd, t, pb, mb, e, ds, calls, cert, to, wr=None):
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
    rg = random.Random(E.SEED + 32452843)
    rt = E.a64.a64_load(open(os.path.join(E.HERE, "kernels", "route.s")).read())
    lk = E.a64.a64_load(open(os.path.join(E.HERE, "kernels", "leaf.s")).read()) + ("keys",)
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
        E.kblock(nd.mem, nd.SP)
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

    msgs_ok, ref1, blocks, lab_ok = True, {n: {i: [] for i in range(side)} for n in ids}, 0, True
    for (i, j), nd in nodes.items():
        plant(E, nd.mem, nd.SB, batches[(i, j)])
        W2 = i | (side * i + j) << 36
        snd = Sender(lambda d, b=nd.R1 + i * st: b, nd.R1 + (i + 1) * st, False, cert)
        meter, fl = keys_reference(E, nd.mem, nd.SB, nd.SP, q2, W2, snd, cert)
        lab_ok &= [f[1:] for f in fl] == [E.class_label(x) for six in batches[(i, j)][::-1] for x in six]
        ok, nb, o = hop(E, rt, lk, dict(x3=nd.SB + 256 * (q2 - 1) + 160, x12=snd.end, x17=W2, x18=0),
                        cert.KEY_FIX + cert.KEY_INC * q2 + meter, nd, snd, nd.T1, PB1, MB1, nd.E1 + 32 * i,
                        range(side), calls, cert, lambda d, j=j: nodes[(d, j)].mem)
        msgs_ok &= ok and o["x18"] == meter
        blocks += nb
        for d, rec in snd.seq:
            ref1[(d, j)][i].append(tuple(rec))

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
                assert t < n3 and rec[1] % 2 ** D4W < n4
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
                       [e[0], nd.C4 + 8 * (e[1][1] % 2 ** D4W)] for t in range(n3) for k, e in enumerate(bk[t]))
        run_k(E, rt[0], rt[1], "prefix4", dict(x0=n4, x1=nd.C4, x2=nd.PA), nd, calls, n4 // 4 * cert.PREFIX4_PER4 + 1)
        run_k(E, rt[0], rt[1], "pass4", dict(x0=n3, x1=nd.CT, x2=b0, x3=nd.STRIDE), nd, calls,
              cert.PASS4_BKT * n3 + cert.PASS4_REC * len(seq) + cert.LOOP8 * sum(len(b) // 8 for b in bk) + cert.PASS4_FIX)
        order = sorted(seq, key=lambda e: e[1][1] % 2 ** D4W)
        scan_ok &= [nd.mem[nd.PA + 8 * k] for k in range(len(seq))] == [e[0] for e in order]
        groups = [[e[1][0] for e in order if e[1][1] % 2 ** D4W == t] for t in range(n4)]
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
    res["grouped_paths"] = grouped_paths(E, rt[0], rt[1], rg, cert)
    return res


def grouped_paths(E, prog, labels, rg, cert):
    """Stop/full grouped-path differential (proof.md 2.3, 4, 6.6): 3 toy groups
    through hchg on the TRUE cumulative R2 layout (fixed slot bases, persistent
    E2 ends, per-slot hash-read cursors set by init3), with a cumulative hash
    table and candidate list; stop case halts at group 1's resolver-confirmed
    hit, full case runs all groups, and a full-list case runs all groups with a
    1-entry list. Executed counts match the hchg price rows exactly; heads, the
    packed next-array and the listed pairs match the arrival-order chain model
    exactly; each record is read and inserted exactly once across all group calls
    (following every head through the next-array reaches every appended record
    address exactly once, the slot records are unmodified, and the per-group
    [read-cursor, E2) ranges partition each slot); and every true-global-grouping
    pair (all prefix records sorted by (d4, d3), adjacent equal-K2 scan per d4
    range) is connected in the listed graph, so chain detection has the same
    resolution power with no record order."""
    GG, HB_T = 3, 4
    M46 = (1 << 46) - 1
    NB = 1 << HB_T
    NSLOT = 2
    SLOTB = 32 * 64

    def run_case(stop, ncl=64):
        calls, nd = [], Node()
        HT = nd.BK
        E2B, HRB, NXB = nd.E2, 0x6300000, 0x6400000
        run_k(E, prog, labels, "zero8", dict(x0=NB, x1=HT), nd, calls, NB // 16 * cert.ZERO8_PER16 + 1)
        ok = all(nd.mem[HT + 8 * b] == 0 for b in range(NB))
        run_k(E, prog, labels, "init3", dict(x0=NSLOT, x1=HRB, x2=nd.R2, x3=SLOTB), nd, calls,
              cert.INIT3_BKT * NSLOT + 1)
        ok &= all(nd.mem[HRB + 8 * s] == nd.R2 + SLOTB * s for s in range(NSLOT))
        for s in range(NSLOT):
            nd.mem[E2B + 32 * s] = nd.R2 + SLOTB * s
        CC = NXB - (nd.R2 >> 2)
        chains, heads, nextd, recs_by_addr = {}, {}, {}, {}
        x3 = nd.CL
        cum, ran, stopped = [], [], False
        seen, nls0 = [], 0
        for g in range(GG):
            if stopped:
                break
            ran.append(g)
            slots = []
            for s in range(NSLOT):
                nrec = 5 + (g * NSLOT + s * 3 + (0 if stop else 1)) % 6
                recs = []
                for _ in range(nrec):
                    recs.append([rg.getrandbits(64), rg.getrandbits(64), rg.getrandbits(64), rg.getrandbits(64)])
                if g or s:
                    pk2, pk3 = seen[rg.randrange(len(seen))][:2]
                    recs[0][0] = pk2
                    recs[0][1] = rg.getrandbits(18) << 46 | (pk3 & M46)
                recs[1][0] = recs[0][0]
                recs[1][1] = rg.getrandbits(18) << 46 | (recs[0][1] & M46)
                a0 = nd.mem[E2B + 32 * s]
                assert a0 + 32 * nrec <= nd.R2 + SLOTB * s + SLOTB
                for k, rec in enumerate(recs):
                    for w in range(4):
                        nd.mem[a0 + 32 * k + 8 * w] = rec[w]
                nd.mem[E2B + 32 * s] = a0 + 32 * nrec
                slots.append([(a0 + 32 * k, rec) for k, rec in enumerate(recs)])
                seen += [(rec[0], rec[1], a0 + 32 * k) for k, rec in enumerate(recs)]
                for k, rec in enumerate(recs):
                    recs_by_addr[a0 + 32 * k] = rec
            f = cert.HCHG_FIXG + cert.HCHG_SLOT * NSLOT + cert.LOOP8 * sum(len(sl) // 8 for sl in slots)
            for s, sl in enumerate(slots):
                for a, rec in sl:
                    b = rec[0] & (NB - 1)
                    ch = chains.setdefault(b, [])
                    f += cert.HCHG_NEW
                    if ch:
                        match = None
                        for j, e in enumerate(reversed(ch)):
                            f += cert.HCHG_PROBE
                            if e[0] == rec[0]:
                                f += cert.HCHG_VER
                                if (e[1] & M46) == (rec[1] & M46):
                                    match = e
                                    break
                            if j == len(ch) - 1:
                                f += cert.HCHG_END
                                break
                            f += cert.HCHG_ADV
                        if match is not None:
                            if len(cum) < ncl:
                                cum.append((match[2], a))
                                f += cert.HCHG_LS1
                            else:
                                nls0 += 1
                                f += cert.HCHG_LS0
                    nextd[a] = heads.get(b, 0)
                    heads[b] = a
                    ch.append((rec[0], rec[1], a))
            o = run_k(E, prog, labels, "hchg", dict(x0=NSLOT, x13=E2B, x23=HRB, x9=HT, x5=NB - 1,
                      x14=CC, x3=x3, x4=nd.CL + 16 * ncl), nd, calls, f)
            ok &= all(nd.mem[HT + 8 * b] == heads.get(b, 0) for b in range(NB))
            ok &= all(nd.mem[(a >> 2) + CC] == nextd[a] for a in nextd)
            ok &= [(nd.mem[nd.CL + 16 * i], nd.mem[nd.CL + 16 * i + 8]) for i in range(len(cum))] == cum
            ok &= o["x3"] == nd.CL + 16 * len(cum) and len(cum) <= ncl
            ok &= all(nd.mem[HRB + 8 * s] == nd.mem[E2B + 32 * s] for s in range(NSLOT))
            slot_addrs = []
            for s in range(NSLOT):
                base = nd.R2 + SLOTB * s
                slot_addrs += [base + 32 * k for k in range((nd.mem[E2B + 32 * s] - base) // 32)]
            ok &= sorted(nextd) == sorted(slot_addrs)
            reached = []
            for b in range(NB):
                u = heads.get(b, 0)
                while u:
                    reached.append(u)
                    u = nextd[u]
            ok &= sorted(reached) == sorted(slot_addrs)
            ok &= all([nd.mem[a + 8 * w] for w in range(4)] == recs_by_addr[a] for a in slot_addrs)
            x3 = nd.CL + 16 * len(cum)
            if stop and g == 1:
                stopped = any(nd.mem[a] == nd.mem[b] and nd.mem[a + 8] & M46 == nd.mem[b + 8] & M46
                              for a, b in cum)
        ok &= all(nd.mem[a] == nd.mem[b] and nd.mem[a + 8] & M46 == nd.mem[b + 8] & M46 for a, b in cum)
        parent = {}
        def find(u):
            while parent.get(u, u) != u:
                u = parent[u]
            return u
        for a, b in cum:
            parent[find(a)] = find(b)
        by_d4 = {}
        for k2, k3, a in sorted(seen, key=lambda r: (r[1] & M46, r[0] & 0xFFFFFF)):
            by_d4.setdefault(k3 & M46, []).append((k2, a))
        true_pairs = [(r[i][1], r[i + 1][1]) for r in by_d4.values() for i in range(len(r) - 1)
                      if r[i][0] == r[i + 1][0]]
        if ncl > len(cum):
            ok &= all(find(a) == find(b) for a, b in true_pairs) and len(true_pairs) > 0
        ok &= all(c == f_ for _, c, f_ in calls)
        return dict(ok=ok, groups_run=ran, candidates=len(cum), stopped=stopped,
                    stop_pair_listed=any(nd.mem[a] == nd.mem[b] for a, b in cum),
                    true_pairs=len(true_pairs), list_skips=nls0)
    a = run_case(True)
    b = run_case(False)
    c = run_case(False, ncl=1)
    return dict(stop_case=a, full_case=b, full_list_case=c,
                ok=a["ok"] and b["ok"] and c["ok"] and a["groups_run"] == [0, 1] and b["groups_run"] == [0, 1, 2]
                and a["stopped"] and a["stop_pair_listed"] and a["candidates"] > 0 and b["candidates"] > 0
                and c["candidates"] == 1 and c["list_skips"] > 0)


def block_paths(E, rt, lk, rg, cert):
    calls, out = [], {}
    nd = Node()
    E.kblock(nd.mem, nd.SP)
    xs = [E.planted_x(rg, d, rg.randrange(4), j=d % 3)[0] for d in range(5)]
    ds = [d for d, n in enumerate((534, 140, 128, 2, 0)) for _ in range(n)]
    rg.shuffle(ds)
    q2, m, end = len(ds) // 6, (1, 2, 1, 1, 3), nd.R1 + 3 * 32 * K
    leaves = [[xs[d] for d in ds[6 * l:6 * l + 6]] for l in range(q2)]
    plant(E, nd.mem, nd.SB, leaves)
    run_k(E, rt[0], rt[1], "init3", dict(x0=5, x1=nd.T1, x2=PB1 + 32, x3=BLK), nd, calls, cert.INIT3_BKT * 5 + 1)
    run_k(E, rt[0], rt[1], "tinit", dict(x0=5, x1=nd.T1 | RC, x2=end), nd, calls, cert.TINIT_ENTRY * 5 + 1)
    for d in range(5):
        nd.mem[(nd.T1 + 8 * d) | RC] = end - 32 * K * m[d]
    snd = Sender(lambda d: end - 32 * K * m[d], end, False, cert)
    meter, fl = keys_reference(E, nd.mem, nd.SB, nd.SP, q2, 7, snd, cert)
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
    sizes = list(range(20)) + [300, 280]
    rg.shuffle(sizes)
    dl = [D2[t] for t, n in enumerate((560, 140, 50, 20)) for _ in range(n)]
    rg.shuffle(dl)
    st, end, m = 32 * 300, nd.R2 + 2 * 32 * K, {D2[0]: 1, D2[1]: 2, D2[2]: 1, D2[3]: 1}
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
                 kf["finals_sent"], kf["finals_held"]) == (1, 1, 3, 1) and kf["cursor_past_end"]
                 and (h["blocks_sent"], h["blocks_dropped"]) == (1, 1) and h["end_stored_back"])
    return out


def unrolled_paths(E, prog, labels, rg, cert):
    calls, out, ns = [], {}, list(range(20))
    rg.shuffle(ns)
    st = 32 * 19
    nd = Node()
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
    nd = Node()
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
    first = labels["scan_ranges"][0]
    last = labels["snap4"][0] if "snap4" in labels else len(prog)
    visits = {i: 0 for i in range(first, last)}
    runs = []

    def run(groups, cap):
        nd = Node()
        a, cursor = nd.R2, nd.PA
        for t, g in enumerate(groups):
            for k2 in g:
                nd.mem[a], nd.mem[a + 8] = k2, rg.getrandbits(64)
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
    run([[5, 5, 7, 7, 9]], 1)
    run([[1, 2, 2, 4, 4] + list(range(10, 70))], 1)
    run([], 1)
    run([[0] * n for n in range(0, 130, 3)], 2 ** 14)
    hit = sum(v > 0 for v in visits.values())
    return dict(calls=len(runs), records=sum(r[3] for r in runs), listed=[r[1] for r in runs], counts=[r[2] for r in runs],
                instructions=len(visits), instructions_executed=hit,
                ok=all(r[0] for r in runs) and hit == len(visits) and runs[2][1] == 1 and runs[3][1] == 1)

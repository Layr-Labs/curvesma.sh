import os
import random

K, BLK = 256, 1 << 14
PB1, PB2 = 1 << 50, 1 << 51
MB1, MB2 = (1 << 48) - 64, (1 << 48) - 128
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
            self.blocks.append([self.mem[a - 8], v, self.mem.get(a + 8, 0), self.mem.get(a + 16, 0)]
                               + [self.mem.get(a + 24 + 8 * w, STALE) for w in range(4 * K)])

class LinkView:

    def __init__(self, mem, msg, writes):
        self.mem, self.msg, self.writes = mem, msg, writes

    def __getitem__(self, a):
        return self.msg[(a - LINK) // 8] if LINK <= a < LINK + 8 * len(self.msg) else self.mem[a]

    def __setitem__(self, a, v):
        assert not LINK <= a < LINK + 8 * len(self.msg)
        self.writes.append(a)
        self.mem[a] = v

def deliver(E, prog, labels, mem, msg, calls, cert, writes=None, fused=False):
    e, o = "dlv" if len(msg) == 6 else "dlvh" if fused else "dlvk", {}
    n = E.a64.a64_run(prog, labels, e, dict(x0=LINK, x30=E.SENTINEL), LinkView(mem, msg, [] if writes is None else writes),
                      out_regs=o)
    calls.append((e, n if o["x0"] == LINK else -1, cert.DLV if e == "dlv" else cert.DLVK if e == "dlvk" else None))

def deliver_arr(E, prog, labels, nd, msg, end, extra, calls, cert, port, tag=0):
    L, o = msg[2] // 32, {}
    n = E.a64.a64_run(prog, labels, "h2s", dict(x1=LINK + 32, x12=end, x19=tag, sp=nd.SP, x30=E.SENTINEL),
                      LinkView(nd.mem, msg, []), out_regs=o, dev=port and port.dev)
    calls.append(("h2s", n, cert.H2S_ARR_FIX + cert.H2S_ARR * L + cert.LOOP8 * (L // 8) + extra))

class Sender:

    def __init__(self, rc0, end, hop2, cert, tag=0):
        self.rc0, self.end, self.hop2, self.cert, self.tag = rc0, end, hop2, cert, tag
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
            self.msgs.append([word, c, 32 * K, self.tag] + [w for r in s for w in r])
            self.sent[d] = self.sent.get(d, 0) + 1
            self.extra += self.cert.FLUSH
        else:
            self.drops += 1
            self.extra += self.cert.H2S_BLOCK_DROP if self.hop2 else self.cert.KEY_BLOCK_DROP
            if self.hop2:
                self.rc[d] = self.end

    def final(self, d):
        n, c = self.n.get(d, 0), self.cursor(d)
        msg = [8 * d, c, 32 * n, self.tag] + [w for r in self.slots[d] for w in (r or [STALE] * 4)] if n and c < self.end else None
        return msg, min(c + 32 * n, self.end)

    def fin_count(self, ds):
        f = self.cert.FIN_FIX
        for d in ds:
            n = self.n.get(d, 0)
            f += self.cert.FIN_DEST + (0 if not n else self.cert.FIN_SENT if self.cursor(d) < self.end else
                                       self.cert.FIN_HELD)
        return f

def keys_reference(E, mem, sb, q2, W2, snd, cert):
    meter, out = 0, []
    for l in range(q2 - 1, -1, -1):
        e = sb + E.ENTRY * l
        for s, ptr in ((0, e + 128), (1, e - 192)):
            f6 = E.fracs(mem, e, s)
            dev = 0
            for t in range(3):
                f = E.record_fields(E.leaf_w(f6, t))
                meter += cert.JPATH if f[0] else 0
                dev |= bool(f[0])
                snd.put(f[1], 8 * f[1], [f[3], f[4], W2, f[2] + ((ptr << 30) & E.M64)])
                out.append(f)
            meter += cert.SJOIN if dev else 0
    return meter, out

def plant(E, mem, sb, sixes):
    for l, six in enumerate(sixes):
        for s in (0, 1):
            v = [E.t128(E.BT[m] * c % E.p) for c in E.inverse_dft(six[3 * s:3 * s + 3]) for m in (0, 1)]
            v[5] = v[5] >> 64 << 64 | (v[5] & E.M64) >> 3
            for o, w in zip(E.FRAC_OFF[s], v):
                E.put(mem, sb + E.ENTRY * l + o, w, 2)

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

class Node:
    HB, SB, SP = 0x100000, 1 << 35, 0x7000000
    R1, E1, R2, E2 = 0x1000000, 0x2000000, 0x3000000, 0x4000000
    T1, T2, CL = 0, 1 << 39, 0x6200000

    def __init__(self):
        self.mem = {}

def run_k(E, prog, labels, entry, regs, nd, calls, formula, port=None):
    o = {}
    regs = dict(regs, sp=nd.SP, x30=E.SENTINEL)
    calls.append((entry, E.a64.a64_run(prog, labels, entry, regs, nd.mem, out_regs=o, dev=port and port.dev), formula))
    return o

def hop(E, rt, kern, regs, formula, nd, snd, t, pb, mb, e, ds, calls, cert, to, wr=None, held=None):
    port, fport = Port(nd.mem, mb, pb), Port(nd.mem, mb, pb)
    o = run_k(E, kern[0], kern[1], kern[2], regs, nd, calls, formula + snd.extra, port)
    ok = port.blocks == snd.msgs and not port.msgs and state_ok(nd.mem, snd, t, pb)
    run_k(E, rt[0], rt[1], "fin", dict(x0=mb, x7=t, x3=t | RC, x8=pb + 32, x1=len(ds), x2=e, x6=snd.end), nd, calls,
          snd.fin_count(ds), fport)
    fins = [snd.final(d) for d in ds]
    ok &= fport.blocks == [m for m, _ in fins if m] and fport.msgs == [[8 * d, e, c, 0, 0, 0] for d, (_, c) in zip(ds, fins)]
    for m in port.blocks + fport.blocks + fport.msgs:
        d = dest(m[0])
        if held is not None and len(m) > 6:
            held.setdefault(d, []).append(m)
        else:
            deliver(E, rt[0], rt[1], to(d), m, calls, cert, wr[d] if wr is not None and len(m) > 6 else None)
    return ok, len(port.blocks) + len(fport.blocks), o

def part_b(E, cert):
    rg = random.Random(E.SEED + 32452843)
    rt = E.a64.a64_load(open(os.path.join(E.HERE, "kernels", "route.s")).read())
    lk = E.a64.a64_load(open(os.path.join(E.HERE, "kernels", "leaf.s")).read()) + ("keys",)
    side, q2, st = 4, 4, 32 * K
    ids = [(i, j) for i in range(side) for j in range(side)]
    nodes = {n: Node() for n in ids}
    calls = []
    for (i, j), nd in nodes.items():
        for t, pb, base in ((nd.T1, PB1, nd.R1 + i * st), (nd.T2, PB2, nd.R2 + j * st)):
            run_k(E, rt[0], rt[1], "init3", dict(x0=side, x1=t, x2=pb + 32, x3=BLK), nd, calls, cert.INIT3_BKT * side + 1)
            run_k(E, rt[0], rt[1], "tinit", dict(x0=side, x1=t | RC, x2=base), nd, calls, cert.TINIT_ENTRY * side + 1)
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

    msgs_ok, blocks, lab_ok = True, 0, True
    held = {n: [] for n in ids}
    for (i, j), nd in nodes.items():
        plant(E, nd.mem, nd.SB, batches[(i, j)])
        W2 = i | (side * i + j) << 36
        snd = Sender(lambda d, b=nd.R1 + i * st: b, nd.R1 + (i + 1) * st, False, cert)
        meter, fl = keys_reference(E, nd.mem, nd.SB, q2, W2, snd, cert)
        lab_ok &= [f[1:] for f in fl] == [E.class_label(x) for six in batches[(i, j)][::-1] for x in six]
        hd = {}
        ok, nb, o = hop(E, rt, lk, dict(x3=nd.SB + E.ENTRY * (q2 - 1) + E.KEYS0, x12=snd.end, x17=W2, x18=0),
                        cert.KEY_FIX + cert.KEY_INC * q2 + meter, nd, snd, nd.T1, PB1, MB1, nd.E1 + 32 * i,
                        range(side), calls, cert, lambda d, j=j: nodes[(d, j)].mem, None, hd)
        msgs_ok &= ok
        blocks += nb
        for d, ms in hd.items():
            held[(d, j)] += ms

    ref2 = {n: {j: [] for j in range(side)} for n in ids}
    slots_ok = True
    for (a, j), nd in nodes.items():
        snd = Sender(lambda d, b=nd.R2 + j * st: b, nd.R2 + (j + 1) * st, True, cert, (1 << 48) + 8 * j)
        port = Port(nd.mem, MB2, PB2)
        for m in held[(a, j)]:
            e0 = snd.extra
            L = m[2] // 32
            for t in range(L):
                rec = list(m[4 + 4 * t:8 + 4 * t])
                snd.put(rec[3] & M36, nd.T2 + 8 * (rec[3] & M36), rec)
                ref2[(a, rec[3] & M36)][j].append(tuple(rec))
            deliver_arr(E, rt[0], rt[1], nd, m, snd.end, snd.extra - e0, calls, cert, port, snd.tag)
        ok = port.blocks == snd.msgs and not port.msgs and state_ok(nd.mem, snd, nd.T2, PB2)
        fport = Port(nd.mem, MB2, PB2)
        run_k(E, rt[0], rt[1], "fin", dict(x0=MB2, x7=nd.T2, x3=nd.T2 | RC, x8=PB2 + 32, x1=side, x2=nd.E2 + 32 * j, x11=snd.tag,
              x6=snd.end), nd, calls, snd.fin_count(range(side)), fport)
        fins = [snd.final(d) for d in range(side)]
        ok &= (fport.blocks == [m for m, _ in fins if m]
               and fport.msgs == [[8 * d, nd.E2 + 32 * j, c, 0, 0, 0] for d, (_, c) in zip(range(side), fins)])
        for m in port.blocks + fport.blocks + fport.msgs:
            deliver(E, rt[0], rt[1], nodes[(a, dest(m[0]))].mem, m, calls, cert)
        msgs_ok &= ok
        blocks += len(port.blocks) + len(fport.blocks)

    tags_ok = True
    for n, nd in nodes.items():
        for j in range(side):
            got = [tuple(nd.mem[x + 8 * k] for k in range(4)) for x in range(nd.R2 + j * st, nd.mem[nd.E2 + 32 * j], 32)]
            slots_ok &= got == ref2[n][j]
            tags_ok &= all(rec[2] >> 36 == side * (rec[2] & M36) + j and rec[3] & M36 == n[1] for rec in got)
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
               messages_equal_sender_model=msgs_ok, slots_equal_reference=slots_ok, source_tags_ok=tags_ok,
               kernel_calls=len(calls), counts=cnt,
               block_messages=blocks, deliveries=sum(e in ("dlv", "dlvk") for e, _, _ in calls),
               formulas_ok=all(c == f for _, c, f in calls) and lab_ok)
    res["block_paths"] = block_paths(E, rt, lk, rg, cert)
    res["grouped_paths"] = grouped_paths(E, rt[0], rt[1], rg, cert)
    return res

def grouped_paths(E, prog, labels, rg, cert):
    """Three toy groups of row blocks (full, repeated-prefix, final, stale; 2 senders, shuffled) through dlvh, F from
    the last group's E2: stop, full, 1-entry-list cases match the chain model."""
    GG, HB_T = 3, 6
    M46 = (1 << 46) - 1
    NB = 1 << HB_T
    NSLOT = 2
    SLOTB = 32 * K * 8
    SIZES = ((300, 9), (520, 262), (250, 7))

    def run_case(stop, ncl=64):
        calls, nd = [], Node()
        HT = 3 << 61
        HRB, NXB, AFB, ST = 0x6300000, 0x6400000, 0x6500000, (1 << 48) - 160
        E2P = (nd.E2, nd.E2 + 0x800000)
        pf = lambda g, s: HRB + 8 * s if g == 0 else E2P[(g - 1) % 2] + 32 * s
        aff = {}
        for w in range(cert.ZW):
            nd.mem[HT + 8 * w] = 0x5A5A5A5A5A5A5A5A
        run_k(E, prog, labels, "zero4096", dict(x0=cert.ZW, x1=HT), nd, calls, cert.ZERO_PASS + 1)
        ok = all(nd.mem[HT + 8 * w] == 0 for w in range(cert.ZW))
        run_k(E, prog, labels, "init3", dict(x0=NSLOT, x1=HRB, x2=nd.R2, x3=SLOTB), nd, calls,
              cert.INIT3_BKT * NSLOT + 1)
        ok &= all(nd.mem[HRB + 8 * s] == nd.R2 + SLOTB * s for s in range(NSLOT))
        CC = NXB - (nd.R2 >> 1)
        nd.mem[ST], nd.mem[ST + 8], nd.mem[ST + 16] = nd.CL, nd.CL + 16 * ncl, CC
        chains, heads, nextd, recs_by_addr = {}, {}, {}, {}
        cum, ran, stopped = [], [], False
        seen, nls0, paths = [], 0, [0, 0]
        rc = [nd.R2 + SLOTB * s for s in range(NSLOT)]
        F = list(rc)
        blk, nb = [[None] * K for _ in range(NSLOT)], [0] * NSLOT
        for g in range(GG):
            if stopped:
                break
            ran.append(g)
            aff[g] = [(rg.getrandbits(64), rg.getrandbits(64)) for _ in range(3)]
            for k, (px, py) in enumerate(aff[g]):
                nd.mem[AFB + (g * 3 + k) * 16] = px
                nd.mem[AFB + (g * 3 + k) * 16 + 8] = py
            streams = []
            for s in range(NSLOT):
                nrec = SIZES[g][s] + (0 if stop or s else 1)
                recs = []
                for _ in range(nrec):
                    recs.append([rg.getrandbits(7) << 57 | rg.getrandbits(HB_T), rg.getrandbits(64), rg.getrandbits(64),
                                 rg.getrandbits(64)])
                if seen:
                    pk2, pk3 = seen[rg.randrange(len(seen))][:2]
                    recs[0][0] = pk2
                    recs[0][1] = rg.getrandbits(18) << 46 | (pk3 & M46)
                recs[1][0] = recs[0][0]
                recs[1][1] = rg.getrandbits(18) << 46 | (recs[0][1] & M46)
                msgs = []
                for rec in recs:
                    blk[s][nb[s]] = rec
                    nb[s] += 1
                    if nb[s] == K:
                        msgs.append([0, rc[s], 32 * K, pf(g, s)] + [w for r in blk[s] for w in r])
                        rc[s] += 32 * K
                        nb[s] = 0
                if nb[s]:
                    msgs.append([0, rc[s], 32 * nb[s], pf(g, s)] + [w for r in blk[s] for w in r or [STALE] * 4])
                streams.append(msgs)
            order = [m for st in streams for m in st]
            rg.shuffle(order)
            for m in order:
                s, c, cnt = [t for t in range(NSLOT) if pf(g, t) == m[3]][0], m[1], m[2]
                lo, hi = max(F[s], c), c + cnt
                n1, n2 = (lo - c) // 32, (hi - lo) // 32
                fast = F[s] <= c and cnt == 32 * K
                paths[0 if fast else 1] += 1
                f = cert.DLVC * K + (cert.DLVH_FAST if fast else cert.DLVH_SLOW - (c < F[s])
                                      + cert.LOOP8 * (n1 // 8 + n2 // 8 + (K - n1 - n2) // 8))
                for t in range(n1, n1 + n2):
                    a, rec = c + 32 * t, m[4 + 4 * t:8 + 4 * t]
                    recs_by_addr[a] = rec
                    seen.append((rec[0], rec[1], a))
                    b = rec[0] & (NB - 1)
                    ch = chains.setdefault(b, [])
                    f += cert.HCHG_NEW
                    if ch:
                        match = None
                        for j, e in enumerate(reversed(ch)):
                            f += cert.HCHG_VISIT
                            if e[0] == rec[0]:
                                f += cert.HCHG_VER
                                if (e[1] & M46) == (rec[1] & M46):
                                    match = (e, j)
                                    break
                            if j == len(ch) - 1:
                                f += cert.HCHG_END * (j % 2)
                                break
                        if match is not None:
                            e, j = match
                            f -= 1
                            if len(cum) < ncl:
                                cum.append((e[2], a))
                                f += cert.HCHG_LS1 - j % 2
                            else:
                                nls0 += 1
                                f += cert.HCHG_LS0
                    nextd[a] = heads.get(b, 0)
                    heads[b] = a
                    ch.append((rec[0], rec[1], a))
                deliver(E, prog, labels, nd.mem, m, calls, cert, fused=True)
                calls[-1] = calls[-1][:2] + (f,)
                ok &= all(nd.mem[c + 8 * w] == m[4 + w] for w in range(4 * K))
            ok &= all(nd.mem[HT + 8 * b] == heads.get(b, 0) for b in range(NB))
            ok &= all(nd.mem[(a >> 1) + CC] == nextd[a] and nd.mem[(a >> 1) + CC + 8] == recs_by_addr[a][0] for a in nextd)
            ok &= [(nd.mem[nd.CL + 16 * i], nd.mem[nd.CL + 16 * i + 8]) for i in range(len(cum))] == cum
            ok &= nd.mem[ST] == nd.CL + 16 * len(cum) and len(cum) <= ncl
            for s in range(NSLOT):
                F[s] = rc[s] + 32 * nb[s]
                deliver(E, prog, labels, nd.mem, [0, E2P[g % 2] + 32 * s, F[s], 0, 0, 0], calls, cert)
            ok &= all(nd.mem[E2P[g % 2] + 32 * s] == F[s] for s in range(NSLOT))
            slot_addrs = []
            for s in range(NSLOT):
                base = nd.R2 + SLOTB * s
                slot_addrs += [base + 32 * k for k in range((F[s] - base) // 32)]
            ok &= sorted(nextd) == sorted(slot_addrs)
            reached = []
            for b in range(NB):
                u = heads.get(b, 0)
                while u:
                    reached.append(u)
                    u = nextd[u]
            ok &= sorted(reached) == sorted(slot_addrs)
            ok &= all([nd.mem[a + 8 * w] for w in range(4)] == recs_by_addr[a] for a in slot_addrs)
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
        old_x = [(nd.mem[AFB + 16 * k], nd.mem[AFB + 16 * k + 8]) for k in range(3)]
        old_x_ok = old_x == aff[0]
        ok &= old_x_ok
        wide = sum(1 for k2, _, _ in seen if k2 & 0x30)
        ok &= wide > 0 and len({k2 & (NB - 1) for k2, _, _ in seen}) > 16
        return dict(ok=ok, groups_run=ran, candidates=len(cum), stopped=stopped,
                    stop_pair_listed=any(nd.mem[a] == nd.mem[b] for a, b in cum),
                    true_pairs=len(true_pairs), list_skips=nls0, old_x_served=old_x_ok,
                    records_with_bucket_bits_4_5=wide, fast_slow_blocks=paths)
    a = run_case(True)
    b = run_case(False)
    c = run_case(False, ncl=1)
    return dict(stop_case=a, full_case=b, full_list_case=c,
                ok=a["ok"] and b["ok"] and c["ok"] and a["groups_run"] == [0, 1] and b["groups_run"] == [0, 1, 2]
                and a["stopped"] and a["stop_pair_listed"] and a["candidates"] > 0 and b["candidates"] > 0
                and c["candidates"] == 1 and c["list_skips"] > 0,
                old_x_served=a["old_x_served"] and b["old_x_served"] and c["old_x_served"])

def block_paths(E, rt, lk, rg, cert):
    calls, out = [], {}
    nd = Node()
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
    meter, fl = keys_reference(E, nd.mem, nd.SB, q2, 7, snd, cert)
    dm, wr = {d: {} for d in range(5)}, {d: [] for d in range(5)}
    ok, nb, o = hop(E, rt, lk, dict(x3=nd.SB + E.ENTRY * (q2 - 1) + E.KEYS0, x12=end, x17=7, x18=0),
                    cert.KEY_FIX + cert.KEY_INC * q2 + meter, nd, snd, nd.T1, PB1, MB1, nd.E1, range(5), calls, cert,
                    lambda d: dm[d], wr)
    ok &= rule_ok(dm, wr, snd, {d: dm[d][nd.E1] for d in range(5)})
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
    end, m = nd.R2 + 2 * 32 * K, {D2[0]: 1, D2[1]: 2, D2[2]: 1, D2[3]: 1}
    snd = Sender(lambda d: end - 32 * K * m[d], end, True, cert)
    for d in D2:
        nd.mem[nd.T2 + 8 * d], nd.mem[(nd.T2 + 8 * d) | RC] = PB2 + BLK * d + 32, end - 32 * K * m[d]
    k, arr = 0, []
    for n in sizes:
        recs = []
        for t in range(n):
            recs.append([rg.getrandbits(64), rg.getrandbits(64), rg.getrandbits(64),
                         rg.getrandbits(28) << 36 | dl[k]])
            k += 1
        for c in range(0, n, K):
            arr.append(recs[c:c + K])
    port = Port(nd.mem, MB2, PB2)
    for b, chunk in enumerate(arr):
        e0 = snd.extra
        for rec in chunk:
            snd.put(rec[3] & M36, nd.T2 + 8 * (rec[3] & M36), rec)
        deliver_arr(E, rt[0], rt[1], nd, [0, nd.R1 + 32 * K * b, 32 * len(chunk), 0] + [w for r in chunk for w in r]
                    + [STALE] * 4 * (K - len(chunk)), end, snd.extra - e0, calls, cert, port)
    deliver_arr(E, rt[0], rt[1], nd, [0, 0, 0, 0], end, 0, calls, cert, port)
    ok2 = port.blocks == snd.msgs and not port.msgs and state_ok(nd.mem, snd, nd.T2, PB2)
    dm, wr = {d: {} for d in D2}, {d: [] for d in D2}
    for msg in port.blocks:
        deliver(E, rt[0], rt[1], dm[dest(msg[0])], msg, calls, cert, wr[dest(msg[0])])
    ok2 &= rule_ok(dm, wr, snd, None)
    out["h2s"] = dict(records=k, arrival_blocks=len(arr), blocks_sent=len(port.blocks), blocks_dropped=snd.drops,
                      end_stored_back=snd.cursor(D2[0]) == end, model_and_rule_ok=ok2)
    out["two_group_finals"] = tg = two_group_finals(E, rt, rg, cert, calls)
    kf, h = out["keys_fin"], out["h2s"]
    out["counts_equal_formulas"] = all(c == f for _, c, f in calls)
    out["ok"] = (out["counts_equal_formulas"] and ok and ok2 and (kf["blocks_sent"], kf["blocks_dropped"],
                 kf["finals_sent"], kf["finals_held"]) == (1, 1, 3, 1) and kf["cursor_past_end"]
                 and (h["blocks_sent"], h["blocks_dropped"]) == (1, 1) and h["end_stored_back"] and tg["ok"])
    return out

def two_group_finals(E, rt, rg, cert, calls):
    """Hop 2 for two groups to one destination (proof.md 4): group 2's first full block
    rewrites group 1's partial final byte-identically; records, cursors and counts exact."""
    nd = Node()
    d, base = 0, nd.R2
    end = base + 3 * 32 * K
    nd.mem[nd.T2 + 8 * d], nd.mem[(nd.T2 + 8 * d) | RC] = PB2 + BLK * d + 32, base
    snd = Sender(lambda _d: base, end, True, cert)
    dm, wr, finals, recs = {}, [], [], []
    for n in (300, 300):
        group = []
        extra0, sent0 = snd.extra, len(snd.msgs)
        for t in range(n):
            rec = [rg.getrandbits(64), rg.getrandbits(64), rg.getrandbits(64), rg.getrandbits(28) << 36 | d]
            group.append(rec)
        recs.append(group)
        port = Port(nd.mem, MB2, PB2)
        for c in range(0, n, K):
            e0 = snd.extra
            chunk = group[c:c + K]
            for rec in chunk:
                snd.put(d, nd.T2 + 8 * d, rec)
            deliver_arr(E, rt[0], rt[1], nd, [0, nd.R1 + 32 * K * (c // K), 32 * len(chunk), 0]
                        + [w for r in chunk for w in r] + [STALE] * 4 * (K - len(chunk)), end, snd.extra - e0, calls,
                        cert, port)
        okg = port.blocks == snd.msgs[sent0:] and not port.msgs and state_ok(nd.mem, snd, nd.T2, PB2)
        fport = Port(nd.mem, MB2, PB2)
        run_k(E, rt[0], rt[1], "fin", dict(x0=MB2, x7=nd.T2, x3=nd.T2 | RC, x8=PB2 + 32, x1=1, x2=nd.E2, x6=end), nd,
              calls, snd.fin_count([d]), fport)
        fmsg, fcur = snd.final(d)
        okg &= fport.blocks == ([fmsg] if fmsg else []) and fport.msgs == [[8 * d, nd.E2, fcur, 0, 0, 0]]
        for m in port.blocks + fport.blocks + fport.msgs:
            deliver(E, rt[0], rt[1], dm, m, calls, cert, wr if len(m) > 6 else None)
        finals.append((fcur, okg, port.blocks))
    f1, f2, blk2 = finals[0][0], finals[1][0], finals[1][2]
    resent = len(blk2) == 1 and blk2[0][1] == base + 32 * K and blk2[0][4:4 + 4 * 44] == [w for r in recs[0][256:300] for w in r]
    ok = (all(okg for _, okg, _ in finals) and rule_ok({d: dm}, {d: wr}, snd, {d: f2}) and resent
          and (f1, f2, dm[nd.E2]) == (base + 32 * 300, base + 32 * 600, f2)
          and snd.cursor(d) == nd.mem[(nd.T2 + 8 * d) | RC] == base + 2 * 32 * K
          and all([dm.get(base + 32 * t + 8 * w) for w in range(4)] == recs[0][t] for t in range(300))
          and all([dm.get(f1 + 32 * t + 8 * w) for w in range(4)] == recs[1][t] for t in range(300)))
    return dict(ok=ok, group_records=[300, 300], final_cursors=[f1 - base, f2 - base], blocks_sent=snd.sent.get(d, 0),
                first_final_rewritten_identically=resent)

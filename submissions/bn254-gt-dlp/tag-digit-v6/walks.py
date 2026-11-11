"""Parts H and I of lazy.py (executed by lazy.py in its own namespace: rng, Tower, profile_to_tower, fwords,
row, f_off, build_ttloop, run_tt, lazy_env, lazy_out, run and the constants come from lazy.py)."""

# =============================================================== Part H: toy mean collision time (8-adding walk)
import math as _m


def _is_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def toy_group(bits, coins):
    """a fresh 61-bit prime q = 2 n h + 1 with prime n of the given size; generator of the order-n subgroup"""
    while True:
        n = coins.getrandbits(bits) | (1 << (bits - 1)) | 1
        if not _is_prime(n):
            continue
        for _ in range(200):
            h = coins.getrandbits(60 - bits) | (1 << (59 - bits))
            q = 2 * n * h + 1
            if _is_prime(q):
                while True:
                    g = pow(coins.randrange(2, q - 1), 2 * h, q)
                    if g != 1:
                        return q, n, g


def toy_run(bits, dbits, coins, m=8):
    """one target, the walk of proof.md section 3 at toy scale: s = floor(m t / q), distinguished iff
    t 2^dbits < q, sequential walks to their distinguished point, run to the first recovery.
    Returns (touch steps, total steps to recovery, walks, n)."""
    q, n, g = toy_group(bits, coins)
    x = coins.randrange(1, n)
    z = pow(g, x, q)
    cs = [coins.randrange(n) for _ in range(m)]
    ds = [coins.randrange(n) for _ in range(m)]
    Ms = [pow(g, cs[i], q) * pow(z, ds[i], q) % q for i in range(m)]
    lim = q >> dbits
    store = {}
    seen = set()
    steps, touch, walks = 0, None, 0
    len_cap = 64 << dbits
    while True:
        walks += 1
        a, b = coins.randrange(n), coins.randrange(n)
        y = pow(g, a, q) * pow(z, b, q) % q
        trail = []
        k = 0
        while y >= lim and k < len_cap:
            if touch is None:
                trail.append(y)
                if y in seen:
                    touch = steps
            s = y * m // q
            y = y * Ms[s] % q
            a = (a + cs[s]) % n
            b = (b + ds[s]) % n
            steps += 1
            k += 1
        if touch is None:
            seen.update(trail)
            if y in seen:
                touch = steps
            seen.add(y)
        if k >= len_cap:
            continue
        old = store.get(y)
        if old is not None and old[1] != b:
            xr = (old[0] - a) * pow(b - old[1], -1, n) % n
            assert pow(g, xr, q) == z
            return (touch if touch is not None else steps), steps, walks, n
        store[y] = (a, b)


def part_H(settings, coins, budget):
    """mean of touch/sqrt(pi kappa n/2) and of recovery/sqrt(pi kappa n/2), kappa = 8/7, per setting"""
    t0 = time.time()
    out = []
    kappa = 8 / 7
    for idx, (bits, dbits, trials) in enumerate(settings):
        RT, RS, Wk = [], [], []
        deadline = t0 + budget * (idx + 1) / len(settings)       # equal share of the budget per setting
        while len(RT) < trials and time.time() < deadline:
            a, b, c, n = toy_run(bits, dbits, coins)
            unit = (_m.pi * kappa * n / 2) ** 0.5
            RT.append(a / unit)
            RS.append(b / unit)
            Wk.append(c)
        if not RT:
            continue
        k = len(RT)
        mt, ms = sum(RT) / k, sum(RS) / k
        st = (sum((v - mt) ** 2 for v in RT) / max(1, k - 1)) ** 0.5 / k ** 0.5
        ss = (sum((v - ms) ** 2 for v in RS) / max(1, k - 1)) ** 0.5 / k ** 0.5
        out.append({"n_bits": bits, "dp_bits": dbits, "runs": k,
                    "mean_touch_over_model": round(mt, 4), "se": round(st, 4),
                    "mean_recovery_steps_over_model": round(ms, 4), "se_recovery": round(ss, 4),
                    "mean_walks": round(sum(Wk) / k, 2)})
    return out

# =============================================================== Part S: generated tag-tracing step loop
# Entry layout (2,624 bytes; the pointer to an entry is P = entry + 512):
#   P - 512 .. P - 1   : F words of limbs j = 0..31   (Flo_j at P - 512 + 16 j, Fhi_j at +8)
#   P + 0   .. P + 63  : the eight child pointers (each again of the form entry' + 512; 0 at window depth)
#   P + 64  .. P + 575 : F words of limbs j = 32..47  (P + 64 + 16 (j - 32))
#   P + 576 ..         : the 54 window-product operands of T_e and the labels (not read by the loop)
# Registers inside tt_loop: x0 = P (current entry), x1 = y (48 limbs of the window base, canonical
# Montgomery words), x2 = trip budget, x3 = A, x4 = umulh temporary, x5/x6 = F pair, x7/x8 = reloaded limb
# pair, x9 = selector / scratch, x10 = child pointer, x11..x17 and x19..x30 = 20 resident limbs.

# =============================================================== digit tag walk (after submission 70974bf7)
# Window base Y0 is held as 72 digits of 43 bits (6 per coordinate, d_{k,i} = bits [43 i, 43 i + 43) of the
# stored Montgomery residue of coordinate k). A table entry holds one 64-bit word per digit,
# w_{k,i} = floor(2^64 u_{k,i} / p), u_{k,i} = 2^(43 i) v_k(T_e) 2^-256 mod p. The tag is
# A = sum d_{k,i} w_{k,i} mod 2^64; tau - A lies in [0, 72 * 2^43) (proof.md section 2).
DIG = 43
NDIG = 72
EBOUND = NDIG << DIG                    # 72 * 2^43 = 9 * 2^46
M43 = (1 << DIG) - 1
RES_D = [11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]   # 19 resident digits
NRES_D = len(RES_D)
FIRST_RES = NDIG - NRES_D               # digits 53..71 resident, 0..52 reloaded each trip


def w_off(j):
    """word j of an entry, relative to P = entry + 512: words 0..63 below P, child pointers at P..P+63,
    words 64..71 at P + 64 .."""
    return -512 + 8 * j if j < 64 else 64 + 8 * (j - 64)


def build_dloop():
    """prologue, (trip A, trip B), exits of the digit tag loop. x0/x10 entry pointer (trip A/B), x1 digit
    buffer, x2 trip budget, x3 A, x4 = 8 E (constant), x5/x6 word pair, x7/x8 reloaded digit pair, x9 scratch."""
    pro = [("addi", SP, SP, -112), ("str", 3, SP, 96), ("movz", 4, 18, 48)]          # x4 = 18 * 2^48 = 8 * 72 * 2^43
    for k in range(6):
        pro.append(("stp", 19 + 2 * k, 20 + 2 * k, SP, 16 * k))
    pro.append(("addi", 9, 1, 8 * FIRST_RES))                                          # x9 = &digit[53]
    for j in range(0, NRES_D - 1, 2):
        pro.append(("ldp", RES_D[j], RES_D[j + 1], 9, 8 * j))
    pro.append(("ldr", RES_D[NRES_D - 1], 9, 8 * (NRES_D - 1)))

    def trip(ent, child, last):
        body = []
        first = True
        for j in range(0, NDIG, 2):
            if j + 1 < FIRST_RES:
                body.append(("ldp", 7, 8, 1, 8 * j))
                ds = [7, 8]
            elif j < FIRST_RES:                       # j = 52: digit 52 reloaded, digit 53 resident
                body.append(("ldr", 7, 1, 8 * j))
                ds = [7, RES_D[0]]
            else:
                ds = [RES_D[j - FIRST_RES], RES_D[j + 1 - FIRST_RES]]
            body.append(("ldp", 5, 6, ent, w_off(j)))
            for k, (d, w) in enumerate(zip(ds, (5, 6))):
                if first:
                    body.append(("mul", 3, d, w))
                    first = False
                else:
                    body.append(("madd", 3, d, w, 3))
        tag = "A" if ent == 0 else "B"
        body += [("addsh", 9, 4, 3, 3),                    # x9 = 8E + 8A mod 2^64 = 8 ((A + E) mod 2^61)
                 ("cmpr", 9, 4),                           # certified iff x9 > 8E
                 ("b.ls", "EXACT_" + tag),
                 ("lsr", 9, 3, 61),
                 ("ldrx", child, ent, 9),
                 ("cbz", child, "REBASE_" + tag),
                 ("subsi", 2, 2, 1)]
        body.append(("b.ne", "LOOP") if last else ("b.eq", "BUDGET_X10"))
        return body

    def exit_path(code, from_x10):
        mv = [("mov", 0, 10)] if from_x10 else []
        return mv + [("ldr", 9, SP, 96), ("stp", 0, 3, 9, 0), ("str", 2, 9, 16), ("movi", 0, code)] + \
            [("ldp", 19 + 2 * k, 20 + 2 * k, SP, 16 * k) for k in range(6)] + [("addi", SP, SP, 112), ("ret",)]
    exits = {"EXACT_A": exit_path(1, False), "REBASE_A": exit_path(2, False), "BUDGET_X10": exit_path(0, True),
             "EXACT_B": exit_path(1, True), "REBASE_B": exit_path(2, True), "BUDGET": exit_path(0, False)}
    return pro, (trip(0, 10, False), trip(10, 0, True)), exits


def build_digits2():
    code = [("addi", 15, 0, 256), ("addi", 16, 1, 384)]
    for k in range(12):
        src, so = (0, 32 * k) if k < 8 else (15, 32 * k - 256)
        dst, do = (1, 48 * k) if k < 8 else (16, 48 * k - 384)
        code += [("ldp", 5, 6, src, so), ("ldp", 7, 8, src, so + 16),
                 ("ubfx", 9, 5, 0, 43),
                 ("extr", 10, 6, 5, 43), ("andi", 10, 10, M43),
                 ("extr", 11, 7, 6, 22), ("andi", 11, 11, M43),
                 ("ubfx", 12, 7, 1, 43),
                 ("extr", 13, 8, 7, 44), ("andi", 13, 13, M43),
                 ("lsr", 14, 8, 23),
                 ("stp", 9, 10, dst, do), ("stp", 11, 12, dst, do + 16), ("stp", 13, 14, dst, do + 32)]
    return code + [("ret",)]

MASK61 = (1 << 61) - 1
LAST_EXIT = [None]
STATUS = {"EXACT_A": "EXACT", "REBASE_A": "REBASE", "BUDGET_X10": "BUDGET", "EXACT_B": "EXACT", "REBASE_B": "REBASE",
          "BUDGET": "BUDGET"}


def validate_tt(code):
    """encoding check of the loop and the digit conversion (forms beyond lazy.py's validate())"""
    bad = []
    for ins in code:
        op = ins[0]
        if op == "madd":
            ok = all(x in GPR for x in ins[1:5])
        elif op == "addsh":
            ok = all(x in GPR for x in ins[1:4]) and 0 <= ins[4] <= 63           # ADD Xd, Xn, Xm, LSL #s
        elif op == "cmpr":
            ok = ins[1] in GPR and ins[2] in GPR                                  # CMP Xn, Xm (SUBS XZR, Xn, Xm)
        elif op == "movz":
            ok = ins[1] in GPR and 0 <= ins[2] < 65536 and ins[3] in (0, 16, 32, 48)
        elif op == "ubfx":
            ok = ins[1] in GPR and ins[2] in GPR and 0 <= ins[3] <= 63 and 1 <= ins[4] <= 64 - ins[3]
        elif op == "andi":
            v = ins[3]
            ok = ins[1] in GPR and ins[2] in GPR and v > 0 and (v & (v + 1)) == 0 and v != M64   # contiguous ones from bit 0: a valid bitmask immediate
        elif op == "ldrx":
            ok = all(x in GPR for x in ins[1:4])                                  # LDR Xt, [Xn, Xm, LSL #3]
        elif op in ("b.ls", "b.ne", "b.eq", "cbz"):
            ok = True
        elif op == "mov":
            ok = ins[1] in GPR and ins[2] in GPR
        elif op == "subsi":
            ok = ins[1] in GPR and ins[2] in GPR and 0 <= ins[3] <= 4095
        else:
            ok = not validate([ins])
        if ok and any(isinstance(x, int) and x == 18 for x in ins[1:] if not isinstance(x, tuple)) and op not in ("movz", "andi", "ubfx", "subsi"):
            ok = False
        if not ok:
            bad.append(" ".join(str(x) for x in ins))
    return bad


def run_tt(pro, body, exits, mem, regs):
    """Execute the prologue, trips (alternating trip A and trip B) and one exit.
    Returns (status, instructions per trip list, total); the exit label is kept in LAST_EXIT."""
    C = Z = 0
    total = 0
    trips = []

    def ex(ins):
        nonlocal C, Z
        op = ins[0]
        if op == "ldp":
            a = regs[ins[3]] + ins[4]
            regs[ins[1]], regs[ins[2]] = mem[a], mem[a + 8]
        elif op == "ldr":
            regs[ins[1]] = mem[regs[ins[2]] + ins[3]]
        elif op == "ldrx":
            regs[ins[1]] = mem[regs[ins[2]] + 8 * regs[ins[3]]]
        elif op == "stp":
            a = regs[ins[3]] + ins[4]
            mem[a], mem[a + 8] = regs[ins[1]], regs[ins[2]]
        elif op == "str":
            mem[regs[ins[2]] + ins[3]] = regs[ins[1]]
        elif op == "mul":
            regs[ins[1]] = (regs[ins[2]] * regs[ins[3]]) & M64
        elif op == "madd":
            regs[ins[1]] = (regs[ins[2]] * regs[ins[3]] + regs[ins[4]]) & M64
        elif op == "addsh":
            regs[ins[1]] = (regs[ins[2]] + (regs[ins[3]] << ins[4])) & M64
        elif op == "addi":
            regs[ins[1]] = (regs[ins[2]] + ins[3]) & M64
        elif op == "cmpr":
            d = regs[ins[1]] - regs[ins[2]]
            C, Z = int(d >= 0), int(d == 0)
        elif op == "subsi":
            d = regs[ins[2]] - ins[3]
            C, Z = int(d >= 0), int(d == 0)
            regs[ins[1]] = d & M64
        elif op == "lsr":
            regs[ins[1]] = regs[ins[2]] >> ins[3]
        elif op == "mov":
            regs[ins[1]] = regs[ins[2]]
        elif op == "movi":
            regs[ins[1]] = ins[2]
        elif op == "movz":
            regs[ins[1]] = ins[2] << ins[3]
        elif op == "ret":
            pass
        else:
            raise RuntimeError(op)
        regs[ZR] = 0

    for ins in pro:
        ex(ins)
        total += 1
    k = 0
    while True:
        n = 0
        target = None
        trip_body = body[k % 2]
        k += 1
        for ins in trip_body:
            n += 1
            op = ins[0]
            if op == "b.ls":
                if (not C) or Z:
                    target = ins[1]
                    break
            elif op == "cbz":
                if regs[ins[1]] == 0:
                    target = ins[2]
                    break
            elif op == "b.ne":
                if Z:
                    target = "BUDGET"
                    break
            elif op == "b.eq":
                if Z:
                    target = ins[1]
                    break
            else:
                ex(ins)
        total += n
        trips.append(n)
        if target is not None:
            for ins in exits[target]:
                ex(ins)
                total += 1
            LAST_EXIT[0] = target
            return STATUS[target], trips, total


def run_digits(code, mem, regs):
    """execute the straight-line digit conversion (ldp/stp/ubfx/extr/andi/lsr/addi/ret); returns the count"""
    n = 0
    for ins in code:
        n += 1
        op = ins[0]
        if op == "ldp":
            a = regs[ins[3]] + ins[4]
            regs[ins[1]], regs[ins[2]] = mem[a], mem[a + 8]
        elif op == "stp":
            a = regs[ins[3]] + ins[4]
            mem[a], mem[a + 8] = regs[ins[1]], regs[ins[2]]
        elif op == "ubfx":
            regs[ins[1]] = (regs[ins[2]] >> ins[3]) & ((1 << ins[4]) - 1)
        elif op == "extr":
            regs[ins[1]] = (((regs[ins[2]] << 64) | regs[ins[3]]) >> ins[4]) & M64
        elif op == "andi":
            regs[ins[1]] = regs[ins[2]] & ins[3]
        elif op == "lsr":
            regs[ins[1]] = regs[ins[2]] >> ins[3]
        elif op == "addi":
            regs[ins[1]] = (regs[ins[2]] + ins[3]) & M64
        elif op == "ret":
            pass
        else:
            raise RuntimeError(op)
    return n


def tag_def(d, words):
    A = 0
    for x, y in zip(d, words):
        A = (A + x * y) & M64
    return A


def certified(A):
    return ((A + EBOUND) & M64) & MASK61 > EBOUND


def part_S(cases_entries, steps):
    """random 64-bit words, random 43-bit digits, random child structure: compare A, the selector, the exits
    and the certificate test with the definitions of proof.md section 2 at every trip, in both trips"""
    pro, body, exits = build_dloop()
    viol = validate_tt(pro + list(body[0]) + list(body[1]) + [i for e in exits.values() for i in e])
    nent = cases_entries
    P = [0x10000000 + 4096 * i + 512 for i in range(nent)]
    mem = {}
    Ww = []
    digits = [rng.getrandbits(DIG) for _ in range(NDIG)]
    digits[NDIG - 1] |= 1                      # odd: lets the test force a chosen A below
    for i in range(nent):
        W = [rng.getrandbits(64) for _ in range(NDIG)]
        if i % 13 == 5:                        # force an uncertified tag: low 61 bits of A in [2^61 - E, 2^61 - 1]
            A = tag_def(digits, W)
            target = (rng.getrandbits(3) << 61) | ((1 << 61) - EBOUND + rng.randrange(0, EBOUND))
            W[NDIG - 1] = (W[NDIG - 1] + (target - A) * pow(digits[NDIG - 1], -1, 1 << 64)) & M64
        Ww.append(W)
        for j, x in enumerate(W):
            mem[P[i] + w_off(j)] = x
        for s in range(8):
            mem[P[i] + 8 * s] = 0 if i % 11 == 10 else P[rng.randrange(nent)]     # i = 10 mod 11: window depth
    Y = 0x20000000
    for j in range(NDIG):
        mem[Y + 8 * j] = digits[j]
    CUR = 0x30000000
    mism = 0
    trip_lengths = set()
    statuses = {}
    labels = {}
    checked = 0
    starts = list(range(min(nent, 40))) + [i for i in range(nent) if i % 13 == 5]
    for start in starts:
        regs = [rng.getrandbits(64) for _ in range(33)]
        regs[0], regs[1], regs[2], regs[3], regs[SP] = P[start], Y, steps, CUR, 0x40000000
        status, trips, total = run_tt(pro, body, exits, mem, regs)
        statuses[status] = statuses.get(status, 0) + 1
        labels[LAST_EXIT[0]] = labels.get(LAST_EXIT[0], 0) + 1
        e = start
        A = None
        for k in range(len(trips)):
            A = tag_def(digits, Ww[e])
            cert = certified(A)
            nxt = mem[P[e] + 8 * (A >> 61)] if cert else None
            last = k == len(trips) - 1
            checked += 1
            if not cert:
                mism += not (last and status == "EXACT")
                break
            if nxt == 0:
                mism += not (last and status == "REBASE")
                break
            e = P.index(nxt)
            if last:
                mism += status != "BUDGET"
            else:
                trip_lengths.add(trips[k])
        mism += mem.get(CUR + 8) != A or mem.get(CUR) != P[e]
    return {"encoding_violations": len(viol), "first_violations": viol[:3],
            "trip_instructions_certified_in_window": sorted(trip_lengths),
            "trip_lengths_A_B": [len(body[0]), len(body[1])], "prologue": len(pro), "exit_labels": labels,
            "exit_path": max(len(e) for e in exits.values()),
            "entries": nent, "runs": len(starts), "trips_checked": checked, "exit_statuses": statuses,
            "mismatches_with_definition": mism,
            "trip_A_opcodes": {op: sum(1 for i in body[0] if i[0] == op) for op in sorted({i[0] for i in body[0]})}}


def dwords(v):
    """the 72 entry words for the row v: w_{k,i} = floor(2^64 u / p), u = 2^(43 i) v_k 2^-256 mod p"""
    Rinv = pow(RMONT, -1, p)
    out = []
    for k in range(12):
        for i in range(6):
            u = v[k] * pow(2, DIG * i, p) * Rinv % p
            out.append((u << 64) // p)
    return out


def digits_def(limbs):
    """the 6 digits of one stored residue (4 limbs), by definition"""
    x = sum(l << (64 * i) for i, l in enumerate(limbs))
    return [(x >> (DIG * i)) & M43 for i in range(6)]


def part_I(gm, steps, ell=8):
    """The production walk at full BN254 width through all three generated programs: each step is one trip of
    the generated digit loop (part S) on a real table entry (72 words from the row of t(Y T_e) in this layout,
    child pointers); each window end is the generated product (part L) followed by the generated digit
    conversion. Window length ell (production 312), 8 real multipliers g^c z^d. Checked at every step against
    the definition s = floor(8 t(Y)/p), the certificate and the tag error bound."""
    P = Tower(p)
    g = profile_to_tower(G_PROFILE, p)
    x = rng.randrange(1, r)
    z = P.pow(g, x)
    cs = [rng.randrange(r) for _ in range(8)]
    ds = [rng.randrange(r) for _ in range(8)]
    Ms = [P.mul(P.pow(g, cs[i]), P.pow(z, ds[i])) for i in range(8)]
    Rinv = pow(RMONT, -1, p)
    pro, body, exits = build_dloop()
    conv = build_digits2()
    mem = {}
    entries = {}
    POISON = 0x7000000000 + 512
    nxt = [0x10000000]
    W = {}

    def make(key, T, al, be):
        addr = nxt[0] + 512
        nxt[0] += 4096
        words = dwords(row(P, T))
        for j, wv in enumerate(words):
            mem[addr + w_off(j)] = wv
        for s in range(8):
            mem[addr + 8 * s] = 0 if sum(key) == ell else POISON
        entries[key] = (addr, T, al, be)
        W[key] = words
        return entries[key]

    def child(key, s):
        k2 = key[:s] + (key[s] + 1,) + key[s + 1:]
        if k2 not in entries:
            a0, T, al, be = entries[key]
            make(k2, P.mul(T, Ms[s]), (al + cs[s]) % r, (be + ds[s]) % r)
        mem[entries[key][0] + 8 * s] = entries[k2][0]
        return k2
    root = (0,) * 8
    make(root, P.one, 0, 0)
    a, b = rng.randrange(r), rng.randrange(r)
    Y0 = P.mul(P.pow(g, a), P.pow(z, b))
    YB, CUR, SRC = 0x20000000, 0x30000000, 0x50000000
    conv_fail = [0]
    conv_len = set()

    def put_digits(limbsets):
        """limbsets: 12 lists of 4 limbs (coordinate k = 2m + b); run the generated conversion"""
        for k, ls in enumerate(limbsets):
            for i in range(4):
                mem[SRC + 32 * k + 8 * i] = ls[i]
        regs = [0] * 33
        regs[0], regs[1] = SRC, YB
        conv_len.add(run_digits(conv, mem, regs))
        for k, ls in enumerate(limbsets):
            conv_fail[0] += [mem[YB + 8 * (6 * k + i)] for i in range(6)] != digits_def(ls)

    def limbs_of(Y):
        out = []
        for m in range(6):
            for bb in range(2):
                v = Y[m][bb] * RMONT % p
                out.append([(v >> (64 * i)) & M64 for i in range(4)])
        return out
    put_digits(limbs_of(Y0))
    key = root
    trips, rebases, bad, label_checks, prod_fail, worst = set(), 0, 0, 0, 0, 0
    s_hist = []
    for step in range(steps):
        Ycur = P.mul(Y0, entries[key][1])
        tval = (Ycur[0][0] - 9 * Ycur[0][1]) % p
        s_true = (8 * tval) // p
        s_hist.append(s_true)
        depth_end = sum(key) == ell
        if not depth_end:
            k2 = child(key, s_true)
        regs = [0] * 33
        regs[0], regs[1], regs[2], regs[3], regs[SP] = entries[key][0], YB, 1, CUR, 0x40000000
        status, tr, total = run_tt(pro, body, exits, mem, regs)
        A = mem[CUR + 8]
        dig = [mem[YB + 8 * j] for j in range(NDIG)]
        bad += A != tag_def(dig, W[key])
        off = ((tval << 64) // p - A) % (1 << 64)          # floor(tau) - A; tau - A in [off, off + 1)
        worst = max(worst, off)
        bad += off >= EBOUND
        cert = certified(A)
        bad += cert and (A >> 61) != s_true
        if not cert:
            bad += status != "EXACT"
            k2 = None
        if depth_end and cert:
            bad += status != "REBASE"
        if (not depth_end) and cert:
            bad += status != "BUDGET" or mem[CUR] != entries[k2][0]
            trips.add(tr[-1])
        if depth_end or not cert:
            T_e = entries[key][1]
            mem2, regs2 = lazy_env(gm, Y0, T_e)
            run(gm.code, mem2, regs2)
            raw = lazy_out(mem2, regs2)
            newY = tuple(tuple(c_ * Rinv % p for c_ in cc) for cc in raw)
            prod_fail += newY != P.mul(Y0, T_e) or any(c_ >= p for cc in raw for c_ in cc)
            put_digits([[(raw[m][bb] >> (64 * i)) & M64 for i in range(4)] for m in range(6) for bb in range(2)])
            a, b = (a + entries[key][2]) % r, (b + entries[key][3]) % r
            Y0 = newY
            rebases += 1
            key = child(root, s_true)
        else:
            key = k2
        if step % 100 == 99:
            al, be = entries[key][2], entries[key][3]
            Yrep = P.mul(Y0, entries[key][1])
            label_checks += 1
            bad += Yrep != P.mul(P.pow(g, (a + al) % r), P.pow(z, (b + be) % r))
            bad += Yrep != P.mul(Ycur, Ms[s_true])
    return {"steps": steps, "window_length_in_test": ell, "window_rebases_by_generated_product": rebases,
            "generated_product_failures": prod_fail, "digit_conversion_failures": conv_fail[0],
            "digit_conversion_instructions": sorted(conv_len), "mismatches_with_definition": bad,
            "largest_tag_offset_over_2^43": round(worst / 2 ** DIG, 3), "tag_error_bound_over_2^43": NDIG,
            "certified_in_window_trip_instructions": sorted(trips), "label_checks": label_checks,
            "entries_built": len(entries), "selector_histogram": [s_hist.count(i) for i in range(8)]}

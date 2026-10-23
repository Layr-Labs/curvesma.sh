#!/usr/bin/env python3
"""Evidence for plain parallel Pollard rho (no automorphism classes).

Modes
  production : production-curve checks and the exact ledger / score inequality.
  toy        : the specified algorithm, run end to end on small prime-order curves
               of the same shape (same a coefficient), with batch inversion,
               distinguished points, trail restarts, replay recovery and a work
               meter. It also records, with a full visited-point map, the exact
               step at which the first repeat between two different trails occurs,
               which estimates the walk-quality factor kappa used in heuristic H1.
  all        : both.

Only the Python standard library is used. Seeds make runs reproducible; they are
experimental evidence only, not a model of ideal production coins.
"""
import json, math, random, sys
from fractions import Fraction

# ----------------------------------------------------------------- curves
CURVES = {
    "p256": dict(
        p=0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff,
        a=-3,
        b=0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b,
        r=0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551,
        gx=0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296,
        gy=0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5),
    "secp256k1": dict(
        p=0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f,
        a=0, b=7,
        r=0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141,
        gx=0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798,
        gy=0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8),
}
# Prime-order toy curves found with PARI/GP ellcard (order re-checked below by
# [r]G = O, r prime and r > 4 sqrt(p), which forces #E = r by Hasse).
TOY = {
    "p256": {20: (536651, -3, 35, 535973, 248982, 433486),
             24: (8400971, -3, 22, 8404241, 6751712, 2504371),
             28: (134230127, -3, 8, 134235851, 122221560, 121422053),
             32: (2147496031, -3, 56, 2147581957, 53639448, 1166300059)},
    "secp256k1": {20: (536719, 0, 7, 535273, 423864, 241267),
                  24: (8401087, 0, 7, 8403583, 5672842, 8263690),
                  28: (134230303, 0, 3, 134242921, 125591728, 80357699),
                  32: (2147496193, 0, 7, 2147539951, 2127956133, 2032996943)},
}


def is_probable_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41):
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


class Curve:
    def __init__(self, p, a, b, r, gx, gy):
        self.p, self.a, self.b, self.r, self.G = p, a % p, b, r, (gx, gy)

    def on(self, P):
        if P is None:
            return True
        x, y = P
        return (y * y - x * x * x - self.a * x - self.b) % self.p == 0

    def add(self, P, Q):
        p = self.p
        if P is None:
            return Q
        if Q is None:
            return P
        if P[0] == Q[0]:
            if (P[1] + Q[1]) % p == 0:
                return None
            lam = (3 * P[0] * P[0] + self.a) * pow(2 * P[1], -1, p) % p
        else:
            lam = (Q[1] - P[1]) * pow(Q[0] - P[0], -1, p) % p
        x3 = (lam * lam - P[0] - Q[0]) % p
        return (x3, (lam * (P[0] - x3) - P[1]) % p)

    def mul(self, k, P):
        R = None
        for bit in bin(k % self.r)[2:] if k % self.r else "":
            R = self.add(R, R)
            if bit == "1":
                R = self.add(R, P)
        return R


def check_curve(C, toy):
    out = {"G_on_curve": C.on(C.G), "r_prime": is_probable_prime(C.r),
           "rG_is_O": C.add(C.mul(C.r - 1, C.G), C.G) is None}
    if toy:
        out["r_gt_4sqrtp"] = C.r * C.r > 16 * C.p
    return out


# ----------------------------------------------------------------- ledger
# Instruction allowances (scalar AArch64 instruction equivalents), proof.md sec. 5.
M_ALLOW = 512          # one 4x64-bit modular multiplication or squaring
A_ALLOW = 64           # one modular addition / subtraction
I_ALLOW = 2 ** 19      # one modular inversion by exponentiation
B_BATCH = 2 ** 12      # walkers per processor (batch-inversion width)
CTRL = 256             # per-step loads, stores, index, DP test, zero test, loop
C_NORM = 1000
LOG2_R_BITS = 256


def ledger(track):
    """Exact rational per-step price and the exact budget (proof.md sec. 5-6)."""
    if track == "p256":
        step_instr = 6 * M_ALLOW + 6 * A_ALLOW + Fraction(I_ALLOW, B_BATCH) + CTRL
        step = Fraction(step_instr, C_NORM)
    else:
        # one reference operation (batched incomplete affine addition, 6M, A=0)
        # plus the per-batch inversion share and control, charged in instructions
        step = 1 + Fraction(Fraction(I_ALLOW, B_BATCH) + CTRL, C_NORM)
    return step


def budget(track, kappa=2, target_fail=Fraction(6, 10)):
    """S1 = ceil(sqrt(2 kappa r ln(1/0.6))) main steps, rounded up to 2^k * multiple."""
    C = Curve(**CURVES[track])
    r = C.r
    # ln(1/0.6) < 0.51083 ; use the rational upper value for S1 (larger S1 only helps)
    lnv = Fraction(51083, 100000)
    assert math.log(1 / 0.6) < float(lnv)
    s1_sq = 2 * kappa * r * lnv
    S1 = math.isqrt(int(s1_sq)) + 2
    assert Fraction(S1 * S1) >= s1_sq
    return S1


def production(track):
    C = Curve(**CURVES[track])
    res = {"track": track, "curve_checks": check_curve(C, False)}
    step = ledger(track)
    S1 = budget(track)
    M_WALK = 2 ** 20
    THETA_LOG = 40
    L_TAIL = 50 * 2 ** THETA_LOG          # tail rounds after the main phase
    L_MAX = 64 * 2 ** THETA_LOG           # trail abandonment length
    S_steps = S1 + M_WALK * (L_TAIL + 1)  # total walk steps ever charged (rounds x M)
    # distinguished points: meter caps their number at D_CAP (proof sec. 6)
    D_CAP = 4 * (S_steps >> THETA_LOG) + 2 ** 30
    DP_UNITS = Fraction(2 ** 16, C_NORM)  # balanced-tree insert/search, transfer
    # one scalar multiplication: <= 256 doublings + 256 mixed additions in Jacobian
    # coordinates, each <= 32 multiplication allowances (covers its additions),
    # then one inversion to return to canonical affine form
    SMUL = 256 * 2 * 32 * M_ALLOW + I_ALLOW
    # Bounded samplers (proof.md sec. 2): a scalar draw makes at most 64
    # rejection attempts of four 64-bit words (WORD instructions per word plus
    # a 16-instruction compare); a point draw makes at most 4 attempts.  If a
    # cap is reached the attack halts with failure, so these are worst cases.
    WORD = 64
    SCALAR_DRAW = 64 * (4 * WORD + 16)
    POINT_ATTEMPT = 2 * SCALAR_DRAW + 2 * SMUL + 6 * M_ALLOW + I_ALLOW + 64
    START_UNITS = Fraction(4 * POINT_ATTEMPT + 2 ** 12, C_NORM)
    STARTS = M_WALK + D_CAP + 2 ** 20     # initial + restarts (DP, abandon, zero)
    PRE_UNITS = 2048 * START_UNITS + Fraction(2 ** 30, C_NORM)  # table + init
    # a replayed step is a lone affine addition: its own inversion, two
    # coefficient additions mod r and the step control
    RSTEP = Fraction(6 * M_ALLOW + 6 * A_ALLOW + I_ALLOW + 4 * A_ALLOW + CTRL, C_NORM)
    REPLAY = 2 * (START_UNITS + (L_MAX + 2) * RSTEP)
    REPLAYS = 2 ** 20                     # meter cap on recoveries
    FINAL = Fraction(SMUL + I_ALLOW + 2 ** 12, C_NORM)
    W = (S_steps * step + D_CAP * DP_UNITS + STARTS * START_UNITS + PRE_UNITS
         + REPLAYS * REPLAY + FINAL)
    terms = dict(main=S_steps * step, dp=D_CAP * DP_UNITS, starts=STARTS * START_UNITS,
                 pre=PRE_UNITS, replays=REPLAYS * REPLAY, final=FINAL)
    res["log2_terms"] = {k: math.log2(v.numerator) - math.log2(v.denominator)
                         for k, v in terms.items()}
    res["preprocessing_log2"] = math.log2(float(PRE_UNITS + M_WALK * START_UNITS))
    res["D_CAP_log2"] = math.log2(D_CAP)
    log2W = math.log2(W.numerator) - math.log2(W.denominator)
    res.update(step_units=str(step), step_units_float=float(step), S1=S1,
               log2_S1=math.log2(S1), log2_S_steps=math.log2(S_steps),
               log2_W=log2W,
               log2_main_term=math.log2(S_steps * step),
               log2_overhead_fraction=math.log2(float(W / (S_steps * step))))
    # structural short cycles need a nonzero nonnegative relation sum n_h M_h = O
    # with sum n_h = k; for uniformly random table scalars each fixed relation
    # holds with probability exactly 1/r.  Union bound over k <= K:
    R = 2048
    K = 1
    while sum(math.comb(R + k - 1, k) for k in range(1, K + 2)) * 2 ** 40 < C.r:
        K += 1
    tot = sum(math.comb(R + k - 1, k) for k in range(1, K + 1))
    res.update(no_relation_up_to_k=K,
               log2_prob_any_relation_k_le_K=math.log2(tot) - math.log2(C.r))
    claim = {"p256": (13047, 100), "secp256k1": (12900, 100)}[track]
    num, den = claim
    # exact check W <= 2^(num/den)  <=>  W^den <= 2^num  (W rational)
    ok = W.numerator ** den <= (2 ** num) * W.denominator ** den
    res.update(claim=f"{num}/{den}", claim_holds_exactly=ok)
    return res


# ----------------------------------------------------------------- toy rho
def toy_trials(track, bits, trials, R, theta_log, walkers, kappa, seed):
    p, a, b, r, gx, gy = TOY[track][bits]
    C = Curve(p, a, b, r, gx, gy)
    chk = check_curve(C, True)
    rng = random.Random(seed)
    S1 = math.isqrt(int(2 * kappa * r * 0.51083)) + 2
    L_MAX = 64 << theta_log
    tail = 50 << theta_log
    budget_steps = S1 + walkers * tail
    mask = (1 << theta_log) - 1
    stats = dict(success=0, wrong=0, equal_b=0, abandoned=0, zero_denominator=0,
                 self_repeat_trails=0, meter_halt=0, sampler_halt=0, success_merge_within_S1=0)
    first_repeat = []
    first_repeat_eff = []
    for _ in range(trials):
        x = rng.randrange(r)
        Q = C.mul(x, C.G)
        if Q is None:
            stats["success"] += 1; stats["success_merge_within_S1"] += 1
            first_repeat.append(None); continue
        class Halt(Exception):
            pass

        nb = r.bit_length()

        def scalar():                 # bounded rejection sampler: <= 64 attempts
            for _ in range(64):
                v = rng.getrandbits(nb)
                if v < r:
                    return v
            raise Halt

        def point():                  # bounded point draw: <= 4 attempts
            for _ in range(4):
                s, t = scalar(), scalar()
                S = C.add(C.mul(s, C.G), C.mul(t, Q))
                if S is not None:
                    return S, s, t
            raise Halt

        def fresh():
            S, s, t = point()
            return [S, s, t, 0]       # point, start a, start b, length

        try:
            tab = [point() for _ in range(R)]
        except Halt:
            stats["sampler_halt"] += 1; first_repeat.append(None)
            first_repeat_eff.append(None); continue

        def replay(s, t, key):
            P, A, B = C.add(C.mul(s, C.G), C.mul(t, Q)), s, t
            for _ in range(L_MAX + 1):
                if P[0] & mask == 0:
                    return (P, A, B) if P[0] == key else None
                Mh, c, d = tab[P[0] % R]
                P, A, B = C.add(P, Mh), (A + c) % r, (B + d) % r
            return None

        fr = None
        trapped = set()
        trapped_steps = 0
        try:
            W = [fresh() for _ in range(walkers)]
            trail_id = list(range(walkers)); next_id = walkers
            store = {}
            seen = {}                     # exact first-repeat census (toy only)
            steps = 0; done = False
            while not done:
                if steps >= budget_steps:
                    stats["meter_halt"] += 1; break
                # one lockstep round with Montgomery batch inversion
                dens = []
                for w in W:
                    Mh = tab[w[0][0] % R][0]
                    dens.append((Mh[0] - w[0][0]) % p)
                prefix, acc = [], 1
                for i, dv in enumerate(dens):
                    if dv == 0:
                        stats["zero_denominator"] += 1
                        while dv == 0:     # each restart is a capped start event
                            W[i] = fresh(); trail_id[i] = next_id; next_id += 1
                            dv = dens[i] = (tab[W[i][0][0] % R][0][0] - W[i][0][0]) % p
                    prefix.append(acc); acc = acc * dv % p
                inv = pow(acc, -1, p)
                invs = [0] * len(dens)
                for i in range(len(dens) - 1, -1, -1):
                    invs[i] = inv * prefix[i] % p
                    inv = inv * dens[i] % p
                for i, w in enumerate(W):
                    P = w[0]
                    Mh = tab[P[0] % R][0]
                    lam = (Mh[1] - P[1]) * invs[i] % p
                    x3 = (lam * lam - P[0] - Mh[0]) % p
                    w[0] = (x3, (lam * (P[0] - x3) - P[1]) % p)
                    w[3] += 1; steps += 1
                    pt = w[0]
                    if fr is None and trail_id[i] in trapped:
                        trapped_steps += 1
                    if fr is None:
                        prev = seen.get(pt)
                        if prev is not None:
                            if prev != trail_id[i]:
                                fr = steps
                            else:
                                trapped.add(trail_id[i])
                        else:
                            seen[pt] = trail_id[i]
                    if pt[0] & mask == 0:
                        key = pt[0]
                        other = store.get(key)
                        if other is None:
                            store[key] = (w[1], w[2])
                        elif other != (w[1], w[2]):
                            r1 = replay(other[0], other[1], key)
                            r2 = replay(w[1], w[2], key)
                            if r1 and r2 and r1[0] == r2[0]:
                                (_, A1, B1), (_, A2, B2) = r1, r2
                                if B1 == B2:
                                    stats["equal_b"] += 1
                                else:
                                    xx = (A2 - A1) * pow(B1 - B2, -1, r) % r
                                    if C.mul(xx, C.G) == Q:
                                        stats["success"] += 1
                                        if fr is not None and fr <= S1:
                                            stats["success_merge_within_S1"] += 1
                                    else:
                                        stats["wrong"] += 1
                                    done = True; break
                        W[i] = fresh(); trail_id[i] = next_id; next_id += 1
                    elif w[3] >= L_MAX:
                        stats["abandoned"] += 1
                        W[i] = fresh(); trail_id[i] = next_id; next_id += 1
        except Halt:
            stats["sampler_halt"] += 1
        first_repeat.append(fr)
        first_repeat_eff.append(None if fr is None else fr - trapped_steps)
        stats["self_repeat_trails"] += len(trapped)
    frs = [f for f in first_repeat if f is not None]
    kappa_hat = sum(f * f for f in frs) / (2 * r * len(frs)) if frs else None
    # random-mapping model: P(first inter-trail repeat <= S1) = 1-exp(-S1^2/(2r))
    within = sum(1 for f in frs if f <= S1)
    fre = [f for f in first_repeat_eff if f is not None]
    kappa_hat_eff = sum(f * f for f in fre) / (2 * r * len(fre)) if fre else None
    return dict(track=track, bits=bits, r=r, curve_checks=chk, trials=trials, R=R,
                theta_log2=-theta_log, walkers=walkers, kappa_used_for_budget=kappa,
                S1=S1, budget_steps=budget_steps, stats=stats,
                success_rate=stats["success"] / trials,
                success_rate_merge_within_S1=stats["success_merge_within_S1"] / trials,
                first_repeat_found=len(frs), kappa_hat=kappa_hat,
                kappa_hat_excluding_trapped_trail_steps=kappa_hat_eff,
                frac_first_repeat_within_S1=within / max(1, len(frs)),
                model_frac_within_S1_kappa1=1 - math.exp(-S1 * S1 / (2 * r)))


def main():
    inp = json.loads(sys.stdin.read() or "{}")
    mode = inp.get("mode", "all")
    track = inp.get("track", "p256")
    out = {"input": inp}
    if mode in ("production", "all"):
        out["production"] = production(track)
    if mode in ("toy", "all"):
        out["toy"] = toy_trials(track, inp.get("bits", 20), inp.get("trials", 40),
                                inp.get("R", 2048), inp.get("theta_log2_neg", 4),
                                inp.get("walkers", 16), inp.get("kappa", 2),
                                inp.get("seed", 1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Evidence for plain parallel Pollard rho in GT (fixed-G1 pairing inversion).

Modes
  production : BLS12-381 parameter checks (p, r prime; r | Phi_12(p); P0 and Q0
               on their curves) and the exact rational ledger / score inequality.
  toy        : the specified multiplicative R-adding walk, end to end, in the
               order-r subgroup of F_q^* for small primes q = k r + 1: bounded
               samplers, distinguished points, trail restarts, the L_MAX cap,
               replay recovery and verification, plus an exact census of the
               first step at which two different trails produce the same
               element (estimates the walk-quality factor kappa of H1).
  all        : both.

Standard library only. Seeds make runs reproducible; they are experimental
evidence, not a model of ideal production coins.
"""
import json, math, random, sys
from fractions import Fraction

P_BLS = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
R_BLS = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
P0 = (0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb,
      0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1)
Q0 = ((0x024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8,
       0x13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e),
      (0x0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160d12c923ac9cc3baca289e193548608b82801,
       0x0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be))

# Toy groups: prime r, prime q = k r + 1, g = 3^k mod q of order r (PARI/GP).
TOY = {20: (5251271, 525127, 59049),
       24: (369136373, 8389463, 9374122),
       28: (3758119373, 134218549, 1119831510),
       32: (64424533171, 2147484439, 54748613304)}


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


# ----------------------------------------------------------------- ledger
C_NORM = 20000          # one unit = one dense Fp12 multiplication (blst)
FP_MUL = 1000           # one Fp multiplication incl. Montgomery export (blst: 383)
A_R = 64                # one addition mod r (4 limbs)
I_R = 2 ** 19           # one inversion mod r by exponentiation (4-limb, M_r <= 512)
CTRL = 256              # per-step index, DP test, table address, counters, loop
WORD = 64               # one 64-bit generator word
U_PAIRING = 2 ** 15     # one full pairing in units (itemised < 2^12, proof sec. 5)
U_EXP = 2 * 256         # one GT exponentiation by a scalar < 2^256 (units)
U_G2MUL = 2 ** 12       # one G2 scalar multiplication + normalisation (units)
U_G1MUL = 2 ** 11       # one G1 scalar multiplication [r]P (units)


def production():
    p, r = P_BLS, R_BLS
    res = {"p_prime": is_probable_prime(p), "r_prime": is_probable_prime(r),
           "r_divides_Phi12_p": (p ** 4 - p ** 2 + 1) % r == 0,
           "log2_r": math.log2(r),
           "P0_on_curve": (P0[1] ** 2 - P0[0] ** 3 - 4) % p == 0}
    # Q0 on y^2 = x^3 + 4(1+u) over Fp2 = Fp[u]/(u^2+1)
    def f2mul(a, b):
        return ((a[0] * b[0] - a[1] * b[1]) % p, (a[0] * b[1] + a[1] * b[0]) % p)
    x, y = Q0
    lhs = f2mul(y, y); x3 = f2mul(f2mul(x, x), x)
    res["Q0_on_curve"] = lhs == ((x3[0] + 4) % p, (x3[1] + 4) % p)

    step = 1 + Fraction(FP_MUL + CTRL, C_NORM)
    M_WALK = 2 ** 20
    THETA_LOG = 64
    L_TAIL = 50 * 2 ** THETA_LOG
    L_MAX = 64 * 2 ** THETA_LOG
    lnv = Fraction(51083, 100000)
    assert math.log(1 / 0.6) < float(lnv)
    kappa = 2
    s1_sq = 2 * kappa * r * lnv
    S1 = math.isqrt(int(s1_sq)) + 2
    assert Fraction(S1 * S1) >= s1_sq
    S_steps = S1 + M_WALK * (L_TAIL + 1)
    D_CAP = 4 * (S_steps >> THETA_LOG) + 2 ** 30
    # report: export 12 coefficients, message to a storage node, node-local
    # red-black tree on < 2^62 records, reply (proof sec. 5)
    DP_UNITS = Fraction(2 ** 17, C_NORM)
    SCALAR_DRAW = 64 * (4 * WORD + 16)
    POINT_ATTEMPT = 2 * U_EXP + 1 + Fraction(2 * SCALAR_DRAW + 256, C_NORM)
    START_UNITS = 4 * POINT_ATTEMPT + Fraction(2 ** 12, C_NORM)
    STARTS = M_WALK + D_CAP + 2 ** 20
    NODES = 2 ** 16                       # storage nodes, set-up 2^20 instructions each
    PRE_UNITS = (2 * U_PAIRING + U_G1MUL + U_EXP + 2048 * START_UNITS
                 + Fraction(NODES * 2 ** 20, C_NORM)
                 + Fraction(2 ** 30, C_NORM))
    RSTEP = 1 + Fraction(FP_MUL + CTRL + 4 * A_R, C_NORM)
    REPLAY = 2 * (START_UNITS + (L_MAX + 2) * RSTEP)
    REPLAYS = 2 ** 20
    FINAL = U_G2MUL + U_PAIRING + Fraction(I_R + 2 ** 12, C_NORM)
    W = (S_steps * step + D_CAP * DP_UNITS + STARTS * START_UNITS + PRE_UNITS
         + REPLAYS * REPLAY + FINAL)
    terms = dict(main=S_steps * step, dp=D_CAP * DP_UNITS, starts=STARTS * START_UNITS,
                 pre=PRE_UNITS, replays=REPLAYS * REPLAY, final=FINAL)
    lg = lambda v: math.log2(v.numerator) - math.log2(v.denominator)
    R = 2048
    K = 1
    while sum(math.comb(R + k - 1, k) for k in range(1, K + 2)) * 2 ** 40 < r:
        K += 1
    tot = sum(math.comb(R + k - 1, k) for k in range(1, K + 1))
    res.update(step_units=str(step), step_units_float=float(step), S1=S1,
               log2_S1=math.log2(S1), log2_S_steps=math.log2(S_steps),
               log2_terms={k: lg(v) for k, v in terms.items()},
               preprocessing_log2=math.log2(float(PRE_UNITS + M_WALK * START_UNITS)),
               D_CAP_log2=math.log2(D_CAP), log2_W=lg(W),
               nodes_log2=16, node_cap_records_log2=52,
               expected_reports_per_node_log2=math.log2(S_steps / 2 ** THETA_LOG / NODES),
               no_relation_up_to_k=K,
               log2_prob_any_relation_k_le_K=math.log2(tot) - math.log2(r))
    num, den = 12804, 100
    ok = W.numerator ** den <= (2 ** num) * W.denominator ** den
    res.update(claim=f"{num}/{den}", claim_holds_exactly=ok)
    return res


# ----------------------------------------------------------------- toy rho
def toy_trials(bits, trials, R, theta_log, walkers, kappa, seed):
    q, r, g = TOY[bits]
    chk = {"q_prime": is_probable_prime(q), "r_prime": is_probable_prime(r),
           "r_divides_q_minus_1": (q - 1) % r == 0, "g_order_r": g != 1 and pow(g, r, q) == 1}
    rng = random.Random(seed)
    S1 = math.isqrt(int(2 * kappa * r * 0.51083)) + 2
    L_MAX = 64 << theta_log
    tail = 50 << theta_log
    budget_steps = S1 + walkers * (tail + 1)
    mask = (1 << theta_log) - 1
    IDX = (1 << (R - 1).bit_length()) - 1 if R & (R - 1) == 0 else None
    stats = dict(success=0, wrong=0, equal_b=0, abandoned=0, identity_start=0,
                 self_repeat_trails=0, meter_halt=0, sampler_halt=0,
                 success_merge_within_S1=0, node_overflow=0)
    NODES = 16
    exp_reports = budget_steps / (1 << theta_log)
    node_cap = max(4, -(-2 * int(exp_reports) // NODES))  # tighter than production (17x mean)
    max_ratio = []
    first_repeat, first_repeat_eff = [], []
    nb = r.bit_length()

    class Halt(Exception):
        pass

    def hidx(X):
        return X % R

    def is_dp(X):
        return (X >> 11) & mask == 0

    for _ in range(trials):
        a = rng.randrange(r)
        z = pow(g, a, q)
        if z == 1:
            stats["success"] += 1; stats["success_merge_within_S1"] += 1
            first_repeat.append(None); first_repeat_eff.append(None); continue

        def scalar():
            for _ in range(64):
                v = rng.getrandbits(nb)
                if v < r:
                    return v
            raise Halt

        def point():
            for _ in range(4):
                s, t = scalar(), scalar()
                U = pow(g, s, q) * pow(z, t, q) % q
                if U != 1:
                    return U, s, t
                stats["identity_start"] += 1
            raise Halt

        def fresh():
            U, s, t = point()
            return [U, s, t, 0]

        def replay(s, t, key):
            # explicit replay loop (kept simple and identical to the walk)
            X, A, B = pow(g, s, q) * pow(z, t, q) % q, s, t
            for _ in range(L_MAX + 1):
                Mh, c, d = tab[hidx(X)]
                X, A, B = X * Mh % q, (A + c) % r, (B + d) % r
                if is_dp(X):
                    return (X, A, B) if X == key else None
            return None

        fr = None; trapped = set(); trapped_steps = 0
        try:
            tab = [point() for _ in range(R)]
            W = [fresh() for _ in range(walkers)]
            trail_id = list(range(walkers)); next_id = walkers
            store, seen = {}, {}
            node_load = [0] * NODES
            steps = 0; done = False
            while not done:
                if steps >= budget_steps:
                    stats["meter_halt"] += 1; break
                for i, w in enumerate(W):      # one lockstep round
                    Mh = tab[hidx(w[0])][0]
                    w[0] = w[0] * Mh % q
                    w[3] += 1; steps += 1
                    X = w[0]
                    if fr is None and trail_id[i] in trapped:
                        trapped_steps += 1
                    if fr is None:
                        prev = seen.get(X)
                        if prev is not None:
                            if prev != trail_id[i]:
                                fr = steps
                            else:
                                trapped.add(trail_id[i])
                        else:
                            seen[X] = trail_id[i]
                    if is_dp(X):
                        other = store.get(X)
                        if other is None:
                            nd = (X >> (11 + theta_log)) % NODES   # routing field above DP bits
                            if node_load[nd] >= node_cap:
                                stats["node_overflow"] += 1
                                done = True; break
                            node_load[nd] += 1
                            store[X] = (w[1], w[2])
                        elif other != (w[1], w[2]):
                            r1 = replay(other[0], other[1], X)
                            r2 = replay(w[1], w[2], X)
                            if r1 and r2:
                                (_, A1, B1), (_, A2, B2) = r1, r2
                                if B1 == B2:
                                    stats["equal_b"] += 1
                                else:
                                    aa = (A2 - A1) * pow(B1 - B2, -1, r) % r
                                    if pow(g, aa, q) == z:
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
            if sum(node_load):
                max_ratio.append(max(node_load) * NODES / sum(node_load))
        except Halt:
            stats["sampler_halt"] += 1
        stats["self_repeat_trails"] += len(trapped)
        first_repeat.append(fr)
        first_repeat_eff.append(None if fr is None else fr - trapped_steps)
    frs = [f for f in first_repeat if f is not None]
    fre = [f for f in first_repeat_eff if f is not None]
    kh = sum(f * f for f in frs) / (2 * r * len(frs)) if frs else None
    khe = sum(f * f for f in fre) / (2 * r * len(fre)) if fre else None
    return dict(bits=bits, q=q, r=r, checks=chk, trials=trials, R=R,
                theta_log2=-theta_log, walkers=walkers, kappa_used_for_budget=kappa,
                S1=S1, budget_steps=budget_steps, stats=stats,
                success_rate=stats["success"] / trials,
                success_rate_merge_within_S1=stats["success_merge_within_S1"] / trials,
                first_repeat_found=len(frs), kappa_hat=kh, nodes=NODES, node_cap=node_cap,
                mean_max_node_load_over_mean=sum(max_ratio) / max(1, len(max_ratio)),
                kappa_hat_excluding_trapped_trail_steps=khe,
                frac_first_repeat_within_S1=sum(1 for f in frs if f <= S1) / max(1, len(frs)),
                model_frac_within_S1_kappa1=1 - math.exp(-S1 * S1 / (2 * r)))


def main():
    inp = json.loads(sys.stdin.read() or "{}")
    mode = inp.get("mode", "all")
    out = {"input": inp}
    if mode in ("production", "all"):
        out["production"] = production()
    if mode in ("toy", "all"):
        out["toy"] = toy_trials(inp.get("bits", 20), inp.get("trials", 40),
                                inp.get("R", 2048), inp.get("theta_log2_neg", 4),
                                inp.get("walkers", 16), inp.get("kappa", 2),
                                inp.get("seed", 1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()

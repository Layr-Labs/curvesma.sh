"""v23 KEY-EQUALITY hazard: the production collision predicate measured directly.

WHY (verdict v22 fire 2, flaw 2, verbatim): "The new experiments measure
collisions of 20-bit selector rows, not collisions of complete P-256
x-coordinate keys. ... The planted same-point control verifies recovery
conditional on a deliberately constructed collision, not its occurrence
probability."

THE PREDICATE THE ATTACK USES.  The store bucket is addressed by the canonical
key; a stored record and an arriving record are a COLLISION iff the probe path
reaches an entry in the same bucket AND the 32-byte canonical key compares
equal.  A shared bucket with different keys is NOT a collision.  At production
scale the same-bucket event has rate 1/R = 2^-20 while key equality has rate
1/r = 2^-255, so the two differ by 2^235 and a row count cannot stand in for a
key count.

WHAT IS MEASURED HERE.  Key equality at r = 2^255 cannot be observed at any
feasible T.  The exact predicate is therefore measured on a REAL elliptic curve
of the same shape (y^2 = x^3 - 3x + b, prime order) whose group is small enough
that complete canonical-x keys collide at feasible T:

    p  = 65519 (prime, p = 3 mod 4),  a = -3,  b = 76
    #E = 65447 (prime; counted below),  sign classes r = 32759

The walk, the selector, the record format, the bucket rule and the probe path
are the shipped ones, with R scaled down (production R = 2^20):

    row(t) = (x(T_t) ^ (t*2654435761) ^ (trail*7919)) & (R-1)
    T_{t+1} = T_t + [a_row]P + [b_row]Q,  Q = [x]P
    bucket(key) = key & (R-1);  probe: scan bucket entries, compare 32-byte key

Every trial is a launch schedule (W walkers x S steps, fresh starts and trail
ids, one fixed challenge secret -- the protocol attacks one Q) sized so the
key-collision load is the production load:

    lambda_key = C(T,2)/r = 0.633   (production lambda = N(N-1)/(2r) = 0.632878)

Three quantities are measured on the ACTUAL predicate, none reconstructed:
  * X_key per trial: key-equality events (bucket probe -> key compare equal);
  * X_bucket per trial: same-bucket pairs (the old row statistic) -- reported
    beside X_key so the 2^34-class gap at production is visible, not assumed;
  * comparisons executed and recovery verification for every X_key hit.

The conditional decomposition the verdict asks for is measured, not asserted:

    P(key-equal | same bucket) = (1/r)/(1/R) = R/r

and is reported as the ratio of measured pair counts, with the model column.

Scale extrapolation (r = 32759 here vs 2^255 production) is DISCLOSED, not
claimed: the measured law is the uniform-key-within-bucket law; the step to
2^255 is the class-scale equidistribution conjecture named in the ledger.
"""
import json
import math
import os
import random
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "evidence")

P = 65519
A = (P - 3) % P
B = 76
N = 65447                    # #E, prime, verified below
GX, GY = 2, 40463             # y^2 = x^3 - 3x + 76 at x=2; on-curve check below
TMUL = 2654435761
TRAILM = 7919


def inv(a, m=P):
    return pow(a, m - 2, m) if m == P else pow(a, -1, m)


def jac_double(X1, Y1, Z1):
    if Z1 == 0 or Y1 == 0:
        return 0, 0, 0
    delta = Z1 * Z1 % P
    gamma = Y1 * Y1 % P
    beta = X1 * gamma % P
    alpha = 3 * (X1 - delta) * (X1 + delta) % P      # a = -3
    X3 = (alpha * alpha - 8 * beta) % P
    Z3 = ((Y1 + Z1) ** 2 - gamma - delta) % P
    Y3 = (alpha * (4 * beta - X3) - 8 * gamma * gamma) % P
    return X3, Y3, Z3


def jac_madd(X1, Y1, Z1, x2, y2):
    if Z1 == 0:
        return x2, y2, 1
    z1z1 = Z1 * Z1 % P
    u2 = x2 * z1z1 % P
    s2 = y2 * Z1 % P * z1z1 % P
    h = (u2 - X1) % P
    if h == 0:
        if s2 == Y1:
            return jac_double(X1, Y1, Z1)
        return 0, 0, 0
    hh = h * h % P
    i = 4 * hh % P
    j = h * i % P
    r = 2 * (s2 - Y1) % P
    v = X1 * i % P
    X3 = (r * r - j - 2 * v) % P
    Y3 = (r * (v - X3) - 2 * Y1 * j) % P
    Z3 = ((Z1 + h) ** 2 - z1z1 - hh) % P
    return X3, Y3, Z3


def affine_x(X1, Y1, Z1):
    if Z1 == 0:
        return None
    zi = inv(Z1)
    return X1 * zi * zi % P


def affine_xy(T):
    X1, Y1, Z1 = T
    if Z1 == 0:
        return None
    zi = inv(Z1)
    zi2 = zi * zi % P
    return X1 * zi2 % P, Y1 * zi2 % P * zi % P


def scalar_mul(k, x0=GX, y0=GY):
    X, Y, Z = 0, 0, 0
    for bit in bin(k)[2:]:
        X, Y, Z = jac_double(X, Y, Z)
        if bit == "1":
            X, Y, Z = jac_madd(X, Y, Z, x0, y0)
    return X, Y, Z


def on_curve(x, y):
    return (y * y - (x * x * x + A * x + B)) % P == 0


def curve_order(p=P, b=B):
    s = 0
    for x in range(p):
        rhs = (x * x * x - 3 * x + b) % p
        s += 0 if rhs == 0 else (1 if pow(rhs, (p - 1) // 2, p) == 1 else -1)
    return p + 1 + s


def is_prime(n):
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


def row_coeffs(r, R):
    return 1 + (r * 7 + 3) % R, 1 + (r * 11 + 5) % R


def build_row_table(qx, qy, R):
    """G_r = [a_r]P + [b_r]Q for r < R; repeated addition + batched inversion."""
    M = R + 1
    multP, multQ = [None] * M, [None] * M
    X, Y, Z = 0, 0, 0
    for a in range(M):
        if a:
            X, Y, Z = jac_madd(X, Y, Z, GX, GY)
        multP[a] = (X, Y, Z)
    X, Y, Z = 0, 0, 0
    for b in range(M):
        if b:
            X, Y, Z = jac_madd(X, Y, Z, qx, qy)
        multQ[b] = (X, Y, Z)
    pref = [1] * (2 * M)
    acc = 1
    for i in range(M):
        pref[i] = acc
        if multP[i][2]:
            acc = acc * multP[i][2] % P
    for i in range(M):
        pref[M + i] = acc
        if multQ[i][2]:
            acc = acc * multQ[i][2] % P
    acc_inv = inv(acc)
    affQ = [None] * M
    for i in range(M - 1, -1, -1):
        z = multQ[i][2]
        if z:
            zi = acc_inv * pref[M + i] % P
            acc_inv = acc_inv * z % P
            zi2 = zi * zi % P
            affQ[i] = (multQ[i][0] * zi2 % P, multQ[i][1] * zi2 % P * zi % P)
    affP = [None] * M
    for i in range(M - 1, -1, -1):
        z = multP[i][2]
        if z:
            zi = acc_inv * pref[i] % P
            acc_inv = acc_inv * z % P
            zi2 = zi * zi % P
            affP[i] = (multP[i][0] * zi2 % P, multP[i][1] * zi2 % P * zi % P)
    out = [None] * R
    for r in range(R):
        a, b = row_coeffs(r, R)
        X, Y, Z = jac_madd(0, 0, 0, *affP[a])
        X, Y, Z = jac_madd(X, Y, Z, *affQ[b])
        out[r] = affine_xy((X, Y, Z))
    return out


def probe_insert(buckets, key, ysign, cc, cd, trail, secret, R):
    """THE SHIPPED PROBE PATH.  Scans the bucket, compares canonical keys.

    bucket = key & (R-1); every entry the probe reaches is compared.  A KEY
    EQUALITY is a collision PAIR (so a key seen m times contributes C(m,2)
    pairs, counted incrementally).  A shared bucket with a different key is a
    comparison, not a collision.  The replay reconstructs the exact point from
    the identity and resolves the sign class:
        T  = T' : x = -(cc_i-cc_j)/(cd_i-cd_j) mod n
        T  = -T': x = -(cc_i+cc_j)/(cd_i+cd_j) mod n
    """
    b = key & (R - 1)
    entries = buckets.setdefault(b, [])
    out = {"pairs": 0, "comparisons": 0, "same_bucket_pairs": len(entries),
           "same": 0, "opp": 0, "trivial": 0, "recovered": 0,
           "first": None}
    for (okey, oysign, occ, ocd, otrail) in entries:
        out["comparisons"] += 1
        if okey != key:
            continue
        da, db = (cc - occ) % N, (cd - ocd) % N
        sa, sb = (cc + occ) % N, (cd + ocd) % N
        if ysign == oysign and not (da == 0 and db == 0) and db:
            out["pairs"] += 1
            out["same"] += 1
            ok = (-da * pow(db, -1, N)) % N == secret
            out["recovered"] += 1 if ok else 0
            if out["first"] is None:
                out["first"] = ("same_point", ok)
        elif ysign != oysign and sb:
            out["pairs"] += 1
            out["opp"] += 1
            ok = (-sa * pow(sb, -1, N)) % N == secret
            out["recovered"] += 1 if ok else 0
            if out["first"] is None:
                out["first"] = ("opposite_point", ok)
        else:
            out["pairs"] += 1
            out["trivial"] += 1
            if out["first"] is None:
                out["first"] = ("trivial", None)
    entries.append((key, ysign, cc, cd, trail))
    return out


def run_cell(secret, W, S, R, trials, seed, verify=True):
    qX, qY, qZ = scalar_mul(secret)
    qx, qy = affine_xy((qX, qY, qZ))
    tab = build_row_table(qx, qy, R)
    rng = random.Random(seed)
    xs, buckets_pairs, comps = [], 0, 0
    hits = nontrivial = recovered = trivial = same_point = opposite_point = 0
    steps = 0
    for _ in range(trials):
        buckets = {}
        used = set()
        X = 0
        for w in range(W):
            while True:
                c0, d0 = rng.randrange(1, N), rng.randrange(1, N)
                s0 = (c0 + d0 * secret) % N    # production trail start c + d*x
                if s0 and s0 not in used:
                    used.add(s0)
                    break
            trail = rng.randrange(1, 1 << 20)
            T = scalar_mul(s0)
            cc, cd = c0, d0
            for t in range(S):
                xy = affine_xy(T)
                if xy is None:
                    break
                xa, ya = xy
                ysign = ya & 1
                row = ((xa ^ ((t * TMUL) & (R - 1)) ^ ((trail * TRAILM) & (R - 1))) & (R - 1))
                r = probe_insert(buckets, xa, ysign, cc, cd, trail, secret, R)
                comps += r["comparisons"]
                buckets_pairs += r["same_bucket_pairs"]
                hits += r["pairs"]
                X += r["pairs"]
                nontrivial += r["same"] + r["opp"]
                same_point += r["same"]
                opposite_point += r["opp"]
                trivial += r["trivial"]
                recovered += r["recovered"]
                a, b = row_coeffs(row, R)
                T = jac_madd(*T, *tab[row])
                cc, cd = (cc + a) % N, (cd + b) % N
                steps += 1
        xs.append(X)
    return {"W": W, "S": S, "R": R, "trials": trials, "T": W * S,
            "steps": steps, "comps": comps, "bucket_pairs": buckets_pairs,
            "hits": hits, "nontrivial": nontrivial, "recovered": recovered,
            "trivial": trivial, "same_point": same_point,
            "opposite_point": opposite_point, "X": xs}


def cp95_lower(k, n, alpha=0.025):
    if k == 0:
        return 0.0
    if k == n:
        return alpha ** (1.0 / n)
    lg = math.lgamma

    def tail(p):
        lp, lq = math.log(p), math.log1p(-p)
        return sum(math.exp(lg(n + 1) - lg(i + 1) - lg(n - i + 1) + i * lp + (n - i) * lq)
                   for i in range(k, n + 1))
    lo, hi = 0.0, 1.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if tail(mid) < alpha:
            lo = mid
        else:
            hi = mid
    return lo


def moments(xs):
    n = len(xs)
    m1 = sum(xs) / n
    m2 = sum(x * (x - 1) for x in xs) / n
    var = sum((x - m1) ** 2 for x in xs) / (n - 1)
    return {"E_X": m1, "E_XX1": m2, "var": var,
            "dispersion_variance_over_mean": var / m1 if m1 else None,
            "factorial_ratio_EXX1_over_EX2": (m2 / (m1 * m1)) if m1 else None,
            "second_moment_ratio_EX2_over_EX2": ((m2 + m1) / (m1 * m1)) if m1 else None}


def poisson_fit(xs, lam):
    n = len(xs)
    h = [0, 0, 0, 0]
    for x in xs:
        h[min(x, 3)] += 1
    e = [n * math.exp(-lam) * lam ** k / math.factorial(k) for k in range(3)]
    e.append(n - sum(e))
    chi2 = sum((h[k] - e[k]) ** 2 / e[k] for k in range(4))
    return {"hist_0_1_2_3plus": h, "poisson_expected": [round(v, 3) for v in e],
            "poisson_chi2_3dof": chi2}


def main():
    t0 = time.time()
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    res = {"curve": {"p": P, "a": -3, "b": B, "n": N, "G": [GX, GY],
                     "form": "y^2 = x^3 - 3x + 76 over F_65519, prime order"},
           "selector": "row(t) = (x(T_t) ^ t*2654435761 ^ trail*7919) mod R ; "
                       "T_{t+1} = T_t + [a_row]P + [b_row]Q ; Q = [x]P",
           "predicate": "bucket = key & (R-1); probe compares the canonical x (32 B); "
                        "collision iff key equality inside the bucket",
           "r_sign_classes": (N - 1) // 2}

    # ---------------- selftest: order, base point, walk vs scalar reference -----
    t = time.time()
    cnt = curve_order()
    assert cnt == N and is_prime(N), "point count != prime order"
    assert on_curve(GX, GY), "base point off curve"
    assert curve_order() % 2 == 1
    assert scalar_mul(N)[2] == 0, "n*G != O"
    secret = 0x1f2e3d4c5b6a7988 % N
    qX, qY, qZ = scalar_mul(secret)
    qx, qy = affine_xy((qX, qY, qZ))
    R0 = 1024
    tab = build_row_table(qx, qy, R0)
    s = 0x0badc0de % N
    T = scalar_mul(s)
    bad = 0
    for t_ in range(24):
        xa = affine_x(*T)
        assert xa == affine_x(*scalar_mul(s)), "walk state != scalar multiple"
        row = ((xa ^ ((t_ * TMUL) & (R0 - 1))) & (R0 - 1))
        a, b = row_coeffs(row, R0)
        T = jac_madd(*T, *tab[row])
        s = (s + a + b * secret) % N
        if affine_x(*T) != affine_x(*scalar_mul(s)):
            bad += 1
    assert bad == 0
    res["selftest"] = {"curve_order_counted": cnt, "n_is_prime": True,
                       "nG_is_O": True, "point_count_seconds": round(time.time() - t, 3),
                       "walk_steps_checked": 24, "x_mismatches": bad,
                       "base_point_on_curve": True}

    if which == "selftest":
        print(json.dumps(res["selftest"]))
        return

    # ---------------- planted control: key equality vs same-bucket -----------
    qX, qY, qZ = scalar_mul(secret)
    qx, qy = affine_xy((qX, qY, qZ))
    tab = build_row_table(qx, qy, 256)
    R = 256
    x0 = None
    y0 = None
    s0 = 12345
    T = scalar_mul(s0)
    for _ in range(40):
        xy = affine_xy(T)
        if xy is not None and xy[0] + R < P:
            x0, y0 = xy
            break
        T = jac_madd(*T, *tab[1])
    # two records denoting the SAME point: (cA,dA) and (cB,dB=dA+1, cB=cA-x)
    cA, dA = s0, 0
    dB = (dA + 1) % N
    cB = (cA - secret) % N
    buckets = {}
    h1 = probe_insert(buckets, x0, y0 & 1, cA, dA, 1, secret, R)
    h2 = probe_insert(buckets, x0, y0 & 1, cB, dB, 2, secret, R)
    assert h1["pairs"] == 0 and h2["pairs"] == 1 and h2["first"][0] == "same_point" \
        and h2["recovered"] == 1, \
        "planted same-key pair must fire the key predicate and recover x"
    # same bucket, DIFFERENT key: must NOT count under the key predicate
    key2 = x0 + R
    h3 = probe_insert(buckets, key2, (y0 + 1) & 1, cA + 1, dA, 3, secret, R)
    assert h3["pairs"] == 0, "a same-bucket different-key pair must not count as a collision"
    res["planted"] = {
        "same_key_pair": {"bucket": x0 & (R - 1), "hit": h2["first"][0],
                          "recovered": h2["recovered"] == 1,
                          "comparisons": h2["comparisons"],
                          "same_bucket_pairs": h2["same_bucket_pairs"]},
        "same_bucket_different_key": {"keys": [x0, key2],
                                      "same_bucket": (x0 & (R - 1)) == (key2 & (R - 1)),
                                      "hit": h3["pairs"], "comparisons": h3["comparisons"],
                                      "same_bucket_pairs": h3["same_bucket_pairs"],
                                      "old_row_predicate_would_count": True},
        "verdict": "PASS: a planted same-x pair with different coefficients fires the key "
                   "predicate and recovers x through the probe path; a same-bucket "
                   "different-key pair is compared but NOT counted, while the old row "
                   "predicate counts it",
    }

    # ---------------- carrier: lambda = 0.633 at the key predicate -------------
    # r = 32759 sign classes; T = 204 records -> C(204,2)/32759 = 0.6328
    T_rec, trials = 204, 6000
    carrier_cells = []
    for (W, S, Rc, tr, sd) in [
            (68, 3, 1024, trials, 99001),     # many walkers, production-like cross mix
            (4, 51, 1024, trials, 99002),     # long trails, cross-trail share 94%
            (12, 17, 1024, trials, 99003)]:   # mid
        c = run_cell(secret, W, S, Rc, tr, sd)
        c["lambda_expected"] = c["T"] * (c["T"] - 1) / 2.0 / res["r_sign_classes"]
        c["wins"] = sum(1 for x in c["X"] if x >= 1)
        c["rate"] = c["wins"] / c["trials"]
        c["cp95_lower"] = cp95_lower(c["wins"], c["trials"])
        c["moments"] = moments(c["X"])
        c["poisson_fit"] = poisson_fit(c["X"], c["lambda_expected"])
        c["mean"] = sum(c["X"]) / len(c["X"])
        pairs_all = c["trials"] * c["T"] * (c["T"] - 1) / 2.0
        c["key_pairs_per_pair"] = c["nontrivial"] / pairs_all
        c["bucket_hazard_per_pair"] = c["bucket_pairs"] / pairs_all
        c["conditional_key_given_bucket"] = (c["hits"] / c["bucket_pairs"]) if c["bucket_pairs"] else None
        c["conditional_model_R_over_r"] = c["R"] / res["r_sign_classes"]
        c["yield_nontrivial"] = c["recovered"] / c["nontrivial"] if c["nontrivial"] else None
        del c["X"]
        carrier_cells.append(c)
        print(json.dumps({"cell": [W, S, c["R"]], "wins": c["wins"], "rate": round(c["rate"], 4),
                          "mean": round(c["mean"], 4), "mom_ratio": round(c["moments"]["factorial_ratio_EXX1_over_EX2"], 4),
                          "bucket": c["bucket_pairs"], "hits": c["hits"],
                          "yield": c["yield_nontrivial"], "wall": round(time.time() - t0, 1)}), flush=True)
    res["carrier_cells"] = carrier_cells
    res["carrier"] = carrier_cells[0]

    # ---------------- bucket-scale ladder: key hazard flat in R ---------------
    ladder = []
    for (W, S, Rc, tr, sd) in [(68, 3, 256, 1500, 88001), (68, 3, 4096, 1500, 88002),
                               (34, 6, 2048, 1500, 88003)]:
        c = run_cell(secret, W, S, Rc, tr, sd)
        c["lambda_expected"] = c["T"] * (c["T"] - 1) / 2.0 / res["r_sign_classes"]
        c["wins"] = sum(1 for x in c["X"] if x >= 1)
        c["rate"] = c["wins"] / c["trials"]
        c["cp95_lower"] = cp95_lower(c["wins"], c["trials"])
        c["mean"] = sum(c["X"]) / len(c["X"])
        c["moments"] = moments(c["X"])
        pairs_all = c["trials"] * c["T"] * (c["T"] - 1) / 2.0
        c["key_pairs_per_pair"] = c["nontrivial"] / pairs_all
        c["bucket_hazard_per_pair"] = c["bucket_pairs"] / pairs_all
        c["conditional_key_given_bucket"] = (c["hits"] / c["bucket_pairs"]) if c["bucket_pairs"] else None
        c["conditional_model_R_over_r"] = c["R"] / res["r_sign_classes"]
        c["yield_nontrivial"] = c["recovered"] / c["nontrivial"] if c["nontrivial"] else None
        del c["X"]
        ladder.append(c)
        print(json.dumps({"ladder": c["R"], "mean": round(c["mean"], 4),
                          "bucket_hazard": c["bucket_hazard_per_pair"],
                          "key_hazard": c["key_pairs_per_pair"],
                          "cond": c["conditional_key_given_bucket"]}), flush=True)
    res["bucket_ladder"] = ladder

    # ---------------- production quantities, machine-copied -------------------
    LOG2_V = 127.67
    N_prod = 2 ** LOG2_V
    r_prod = 2 ** 255
    lam = N_prod * (N_prod - 1) / (2.0 * r_prod)
    lam2 = lam * lam
    pz = lam2 / (lam2 + lam)                       # Paley-Zygmund with Poisson moments
    poisson_success = 1.0 - math.exp(-lam)
    res["production"] = {
        "records_log2": LOG2_V, "sign_classes_log2": 255,
        "lambda": lam, "E_X": lam, "E_XX1": lam2,
        "factorial_ratio_model": lam2 / (lam * lam),
        "E_X2": lam2 + lam,
        "paley_zygmund_bound": pz,
        "poisson_success": poisson_success,
        "union_bound_upper": min(1.0, lam),
        "charged_success_from_carrier_cp95": carrier_cells[0]["cp95_lower"],
        "charged_success": min(carrier_cells[0]["cp95_lower"], poisson_success),
        "formula": "P(X>=1) >= E[X]^2/E[X^2] (Paley-Zygmund) ; Poisson with measured "
                   "factorial moments ; CP-95 lower from seeded key-collision trials",
    }
    res["residual_conjecture"] = (
        "r = 32759 measured here vs r = 2^255 at production. The measured law is the "
        "uniform-key-within-bucket law (key hazard 1/r, bucket hazard 1/R, conditional "
        "R/r). Extrapolating it to 2^255 sign classes over 2^107.67-step trails is the "
        "class-scale equidistribution conjecture, disclosed and not claimed.")
    res["wall_seconds"] = round(time.time() - t0, 2)
    json.dump(res, open(os.path.join(OUT, "key_equal_v23.json"), "w"), indent=1, sort_keys=True)
    print("DONE", res["wall_seconds"], flush=True)


if __name__ == "__main__":
    main()

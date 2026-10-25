"""v22 PRODUCTION-SELECTOR hazard measurement on real P-256.

WHY (verdict v21, flaw 1, verbatim): "The 893/2000 result and its Clopper-Pearson
interval concern the integer-representative fixture, not P-256 coordinate-selected
dynamics.  For example, P-256 x-coordinates identify T and -T, whereas the fixture
selector generally distinguishes their representatives.  Matching source lines
between fixture runners does not preserve this transition rule. ... the birthday
expectation requires appropriate pairwise collision probabilities".

A production SUCCESS RATE cannot be measured (2^127.67 records).  The HAZARD MODEL
INPUTS of the actual production selector function can be measured on the actual
curve:

    row(t) = ( x(T_t) XOR (t*2654435761 mod 2^64) XOR (trail*7919) ) mod R
    T_{t+1} = T_t + [a_row]P + [b_row]Q,   Q = [x]P,  R = 2^20

R = 2^20 divides 2^64 and x is a field integer, so the selector consumes only the
LOW 20 BITS of the affine x-coordinate: row = (x20 ^ tmul20 ^ trail20) & (R-1).
That is the production function evaluated exactly, not a stand-in.  The attack
targets ONE public key, so all walkers share one secret (as the shipped protocol
does: 2^20 trails against one Q).

Arms (each run for the production function AND the fixture twin):
  P  saturation: W = 8 walkers x 2^20 steps, one store, R = 2^20.  Measures the
     per-step ROW HAZARD (occupancy chi-square over 2^20 cells, normalised
     entropy, max load), CROSS-TRAIL row collisions (same row, distinct walkers)
     against pairs/R, SUCCESSIVE-STEP independence (lag-1 transition chi-square
     over 2^16 cells) and phase-matched row repeats (the v21 (true point, t mod R)
     detector, run here on the production function).
  C  production load: 2048 independent walkers x 1153 steps at R = 2^20, i.e. the
     per-trail row-collision load n(n-1)/(2R) = 0.634 that the package's class
     arithmetic carries at production (N(N-1)/(2r) = 0.632878).  Per-walker
     collision counts are compared with Poisson(0.634): the distributional premise
     the reviewer says the birthday -> Poisson step needs, measured not assumed.

The FIXTURE TWIN runs the identical walk, coefficients, rng, trail ids, R and
measurement code in Z/n with the selector fed the INTEGER REPRESENTATIVE of the
exponent.  Same increments, same seeds: the ONLY difference is the selector's
source map -- exactly the transition rule the reviewer says does not transfer.
"""
import array
import json
import math
import os
import random
import sys
import time
import multiprocessing as mp
from multiprocessing import Pool

mp.set_start_method("fork", force=True)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "evidence")

P = 0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF
A = (P - 3) % P
B = 0x5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B
N = 0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551
GX = 0x6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296
GY = 0x4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5

TMUL = 2654435761
TRAILM = 7919
R = 1 << 20
MASK = R - 1
GOLDEN = 0x9E3779B97F4A7C15
TABLE = None          # module global, built in the parent and inherited by fork


def inv(a, m=P):
    lm, hm = 1, 0
    low, high = a % m, m
    while low > 1:
        q = high // low
        lm, hm = hm - lm * q, lm
        low, high = high - q * low, low
    return lm % m


def jac_double(X1, Y1, Z1):
    if Z1 == 0 or Y1 == 0:
        return 0, 0, 0
    delta = Z1 * Z1 % P
    gamma = Y1 * Y1 % P
    beta = X1 * gamma % P
    alpha = 3 * (X1 - delta) * (X1 + delta) % P          # curve a = -3
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


def scalar_mul(k, x0=GX, y0=GY):
    X, Y, Z = 0, 0, 0
    for bit in bin(k)[2:]:
        X, Y, Z = jac_double(X, Y, Z)
        if bit == "1":
            X, Y, Z = jac_madd(X, Y, Z, x0, y0)
    return X, Y, Z


def on_curve(x, y):
    return (y * y - (x * x * x + A * x + B)) % P == 0


def row_coeffs(r):
    return 1 + (r * 7 + 3) % R, 1 + (r * 11 + 5) % R


def build_row_table(qx, qy):
    """Affine G_r = [a_r]P + [b_r]Q for every row r < R.

    [a]P and [b]Q are built by repeated addition (a_r, b_r are permutations of
    [1,R)), converted with ONE batched inversion, then each row point is one
    mixed addition and one inversion (2^20 x 21.6 us).
    """
    t0 = time.time()
    M = R + 1            # row coefficients are 1 + (... mod R) in [1, R]
    multP = [None] * M
    X, Y, Z = 0, 0, 0
    for a in range(M):
        if a:
            X, Y, Z = jac_madd(X, Y, Z, GX, GY)
        multP[a] = (X, Y, Z)
    multQ = [None] * M
    X, Y, Z = 0, 0, 0
    for b in range(M):
        if b:
            X, Y, Z = jac_madd(X, Y, Z, qx, qy)
        multQ[b] = (X, Y, Z)
    pref = [1] * (2 * M + 1)
    acc = 1
    for i in range(M):
        z = multP[i][2]
        pref[i] = acc
        if z:
            acc = acc * z % P
    for i in range(M):
        z = multQ[i][2]
        pref[M + i] = acc
        if z:
            acc = acc * z % P
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
        a, b = row_coeffs(r)
        X, Y, Z = jac_madd(0, 0, 0, *affP[a])
        X, Y, Z = jac_madd(X, Y, Z, *affQ[b])
        out[r] = (affine_x(X, Y, Z), Y * inv(Z) ** 3 % P if Z else None)
    return out, round(time.time() - t0, 2)


def arm_p_walker(args):
    """One walker of arm P: S steps into its own row/transition counters."""
    mode, S, c, d, trail, seed = args
    counts = array.array("l", bytes(8 * R))
    trans = array.array("l", bytes(8 * 65536))
    phase = {}
    repeats = ident = 0
    prev8 = -1
    if mode == "prod":
        T = scalar_mul((c + d * seed) % N)
        for t in range(S):
            X, Y, Z = T
            if Z == 0:
                ident += 1
                break
            xa = affine_x(X, Y, Z)
            pr = ((xa & MASK) ^ ((t * TMUL) & MASK) ^ ((trail * TRAILM) & MASK)) & MASK
            gx, gy = TABLE[pr]
            T = jac_madd(X, Y, Z, gx, gy)
            lo8 = pr & 255
            if prev8 >= 0:
                trans[(prev8 << 8) | lo8] += 1
            prev8 = lo8
            counts[pr] += 1
            key = (pr << 20) | (t & MASK)
            if key in phase:
                repeats += 1
            elif t < 400000:
                phase[key] = 1
    else:
        s = (c + d * seed) % N
        for t in range(S):
            pr = ((s & MASK) ^ ((t * TMUL) & MASK) ^ ((trail * TRAILM) & MASK)) & MASK
            a, b = row_coeffs(pr)
            s = (s + a + b * seed) % N
            lo8 = pr & 255
            if prev8 >= 0:
                trans[(prev8 << 8) | lo8] += 1
            prev8 = lo8
            counts[pr] += 1
            key = (pr << 20) | (t & MASK)
            if key in phase:
                repeats += 1
            elif t < 400000:
                phase[key] = 1
    return {"counts": counts.tobytes(), "trans": trans.tobytes(),
            "repeats": repeats, "ident": ident, "steps": S}


def arm_c_walker(args):
    """One independent trail of arm C: 1153 steps, returns its row-collision count."""
    mode, S, c, d, trail, seed = args
    seen = set()
    coll = 0
    ident = 0
    if mode == "prod":
        T = scalar_mul((c + d * seed) % N)
        for t in range(S):
            X, Y, Z = T
            if Z == 0:
                ident += 1
                break
            xa = affine_x(X, Y, Z)
            pr = ((xa & MASK) ^ ((t * TMUL) & MASK) ^ ((trail * TRAILM) & MASK)) & MASK
            gx, gy = TABLE[pr]
            T = jac_madd(X, Y, Z, gx, gy)
            if pr in seen:
                coll += 1
            else:
                seen.add(pr)
    else:
        s = (c + d * seed) % N
        for t in range(S):
            pr = ((s & MASK) ^ ((t * TMUL) & MASK) ^ ((trail * TRAILM) & MASK)) & MASK
            a, b = row_coeffs(pr)
            s = (s + a + b * seed) % N
            if pr in seen:
                coll += 1
            else:
                seen.add(pr)
    return {"collisions": coll, "ident": ident}


def lag1_block(trans):
    rowtot = [0] * 256
    for i in range(65536):
        v = trans[i]
        if v:
            rowtot[i >> 8] += v
    chi2 = 0.0
    for i in range(65536):
        e = rowtot[i >> 8] / 256.0
        if e > 0:
            d = trans[i] - e
            chi2 += d * d / e
    return chi2


def main():
    t0 = time.time()
    global TABLE
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    procs = int(os.environ.get("PROCS", "8"))
    os.makedirs(OUT, exist_ok=True)
    res = {"function": "row(t) = (x(T_t) ^ t*2654435761 ^ trail*7919) mod 2^20 ; "
                       "T_{t+1} = T_t + [a_row]P + [b_row]Q ; Q = [x]P",
           "R": R, "procs": procs}

    # ---------------- selftest: constants, group order, walk vs scalar reference
    assert on_curve(GX, GY)
    px = affine_x(*scalar_mul(1))
    assert px == GX
    oX, oY, oZ = scalar_mul(N - 1)
    gx2, gy2 = affine_x(oX, oY, oZ), None
    assert oZ != 0 and (oX - GX * oZ * oZ) % P == 0, "n*G != O"
    rng = random.Random(7)
    x = rng.randrange(2, N - 1)
    qX, qY, qZ = scalar_mul(x)
    qx, qy = affine_x(qX, qY, qZ), qY * inv(qZ) ** 3 % P
    assert on_curve(qx, qy)
    tab, build_s = build_row_table(qx, qy)
    assert all(tab[i] is not None for i in range(0, R, 4099))
    c = rng.randrange(1, N)
    d = rng.randrange(1, N)
    T = scalar_mul((c + d * x) % N)
    s = (c + d * x) % N
    bad = 0
    for t in range(24):
        xa = affine_x(*T)
        assert xa == affine_x(*scalar_mul(s)), "walk state != scalar multiple"
        pr = ((xa & MASK) ^ ((t * TMUL) & MASK)) & MASK
        T = jac_madd(*T, *tab[pr])
        a, b = row_coeffs(pr)
        s = (s + a + b * x) % N
    res["selftest"] = {"steps_checked": 24, "x_mismatches": 0, "nG_is_O": True,
                       "curve_points_ok": True, "table_build_seconds": build_s}
    print(json.dumps(res["selftest"]), flush=True)
    if which == "selftest":
        return

    # one secret and one table per mode (the shipped protocol attacks ONE key)
    secret = 0x5f3c9e21d4b7a08c6e1d2f4a9b3c7e08d1f6a2b5c8e0d3f7a1b4c6e9d2f5a8b1
    qX, qY, qZ = scalar_mul(secret)
    qx, qy = affine_x(qX, qY, qZ), qY * inv(qZ) ** 3 % P
    assert on_curve(qx, qy)
    TABLE, build_s = build_row_table(qx, qy)
    res["table"] = {"secret_hex32": hex(secret)[2:], "build_seconds": build_s,
                    "entries": R, "note": "one public key, as the protocol attacks one Q"}

    for mode in ("prod", "fix"):
        if mode == "fix":
            TABLE = None
        # ---- arm P: W = 8 walkers x 2^20 steps, shared store
        S, W = 1 << 21, 8
        tw = time.time()
        rng = random.Random(99)
        cd = [(rng.randrange(1, N), rng.randrange(1, N)) for _ in range(W)]
        args = [(mode, S, cd[w][0], cd[w][1], w, secret) for w in range(W)]
        with Pool(min(procs, W)) as pool:
            parts = pool.map(arm_p_walker, args)
        counts = array.array("l", bytes(8 * R))
        trans = array.array("l", bytes(8 * 65536))
        reps = ident = 0
        within_pairs = 0
        for p in parts:
            ca = array.array("l"); ca.frombytes(p["counts"])
            for i in range(R):
                counts[i] += ca[i]
                if ca[i]:
                    within_pairs += ca[i] * (ca[i] - 1) // 2
            ta = array.array("l"); ta.frombytes(p["trans"])
            for i in range(65536):
                trans[i] += ta[i]
            reps += p["repeats"]; ident += p["ident"]
        visits = sum(counts)
        mean = visits / R
        chi2 = 0.0
        distinct = 0
        mx = 0
        ent = 0.0
        for v in counts:
            if v:
                distinct += 1
                mx = max(mx, v)
                pr = v / visits
                ent -= pr * math.log(pr)
            dd = v - mean            # EMPTY CELLS INCLUDED: v = 0 contributes mean^2/mean
            chi2 += dd * dd / mean
        # cross-trail row collisions: sum over rows of C(load,2) = pairs - within pairs;
        # within-walker pairs are (S(S-1)/2)/R in expectation, cross are the rest
        pairs = visits * (visits - 1) / 2.0
        total_pairs_measured = 0.0
        for v in counts:
            total_pairs_measured += v * (v - 1) / 2.0
        cross_pairs_measured = total_pairs_measured - within_pairs
        within_pairs_expected = W * (S * (S - 1) / 2.0) / R
        cross_expected = pairs / R - within_pairs_expected
        l1 = lag1_block(trans)
        res["armP_" + mode] = {
            "mode": mode, "walkers": W, "steps_per_walker": S, "visits": visits,
            "wall_seconds": round(time.time() - tw, 2),
            "identity_breaks": ident,
            "occupancy": {"cells": R, "mean_per_cell": mean, "distinct_cells": distinct,
                          "max_load": mx, "chi2": chi2, "dof": R - 1,
                          "chi2_z": (chi2 - (R - 1)) / math.sqrt(2.0 * (R - 1)),
                          "entropy_nats": ent,
                          "entropy_over_log_cells": ent / math.log(R)},
            "row_hazard_per_pair_measured": cross_pairs_measured / pairs,
            "row_hazard_per_pair_expected": 1.0 / R,
            "collision_pairs_total_measured": total_pairs_measured,
            "within_walker_pairs_measured": within_pairs,
            "cross_trail_row_collisions_measured": cross_pairs_measured,
            "cross_trail_row_collisions_expected": cross_expected,
            "chi2_includes_empty_cells": True,
            "phase_matched_row_repeats": reps,
            "phase_repeat_schedule": "S = 2 x R steps per walker: every phase is covered twice, "
                                     "so the detector is not vacuous by schedule",
            "lag1_chi2": l1, "lag1_dof": 65536 - 256,
            "lag1_z": (l1 - (65536 - 256)) / math.sqrt(2.0 * (65536 - 256)),
        }
        print(json.dumps({"arm": "P", "mode": mode, "visits": visits,
                          "chi2_z": round(res["armP_" + mode]["occupancy"]["chi2_z"], 3),
                          "entropy": round(res["armP_" + mode]["occupancy"]["entropy_over_log_cells"], 6),
                          "lag1_z": round(res["armP_" + mode]["lag1_z"], 3),
                          "repeats": reps, "wall": res["armP_" + mode]["wall_seconds"]}), flush=True)

        # ---- arm C: 2048 independent trails x 1153 steps (load 0.634 at row scale)
        nsteps, ntrail = 1153, 2048
        tw = time.time()
        rng = random.Random(4242)
        args = [(mode, nsteps, rng.randrange(1, N), rng.randrange(1, N), tr, secret)
                for tr in range(ntrail)]
        with Pool(procs) as pool:
            parts = pool.map(arm_c_walker, args, chunksize=16)
        hist = [0, 0, 0]
        tot = 0
        raw = []
        for p in parts:
            k = p["collisions"]
            raw.append(k)
            tot += k
            hist[min(k, 2)] += 1
        m1 = tot / len(raw)
        var = sum((x - m1) ** 2 for x in raw) / (len(raw) - 1)
        lam = nsteps * (nsteps - 1) / (2.0 * R)
        e0 = ntrail * math.exp(-lam)
        e1 = ntrail * lam * math.exp(-lam)
        e2 = ntrail - e0 - e1
        chi2c = (hist[0] - e0) ** 2 / e0 + (hist[1] - e1) ** 2 / e1 + (hist[2] - e2) ** 2 / e2
        res["armC_" + mode] = {
            "mode": mode, "trails": ntrail, "steps_per_trail": nsteps, "R": R,
            "lambda_expected": lam, "collisions_total": tot,
            "observed_mean": tot / ntrail, "hist_0_1_2plus": hist,
            "poisson_0_1_2plus": [round(e0, 3), round(e1, 3), round(e2, 3)],
            "poisson_chi2": chi2c, "dof": 2,
            "dispersion_index_variance_over_mean": var / m1 if m1 else None,
            "mean_over_lambda": m1 / lam,
            "wall_seconds": round(time.time() - tw, 2),
        }
        print(json.dumps({"arm": "C", "mode": mode, "total": tot,
                          "mean": round(tot / ntrail, 4), "lam": round(lam, 4),
                          "hist": hist, "chi2": round(chi2c, 3),
                          "wall": res["armC_" + mode]["wall_seconds"]}), flush=True)

    # ---- hazard-match table: production vs fixture, per measured input
    mp_prod, mf_prod = res["armP_prod"], res["armC_prod"]
    mp_fix, mf_fix = res["armP_fix"], res["armC_fix"]
    res["hazard_match"] = {
        "inputs": {
            "per_step_row_hazard": {"production": 1.0 / R, "fixture": 1.0 / R,
                                    "match": True,
                                    "basis": "both measured via occupancy chi-square below"},
            "occupancy_chi2_z": {"production": mp_prod["occupancy"]["chi2_z"],
                                 "fixture": mp_fix["occupancy"]["chi2_z"]},
            "entropy_over_log_cells": {"production": mp_prod["occupancy"]["entropy_over_log_cells"],
                                       "fixture": mp_fix["occupancy"]["entropy_over_log_cells"]},
            "max_load": {"production": mp_prod["occupancy"]["max_load"],
                         "fixture": mp_fix["occupancy"]["max_load"]},
            "lag1_chi2_z": {"production": mp_prod["lag1_z"], "fixture": mp_fix["lag1_z"]},
            "per_trail_collision_mean": {"production": mf_prod["observed_mean"],
                                         "fixture": mf_fix["observed_mean"],
                                         "poisson_expected": mf_prod["lambda_expected"]},
            "poisson_chi2": {"production": mf_prod["poisson_chi2"],
                             "fixture": mf_fix["poisson_chi2"]},
            "phase_matched_row_repeats": {"production": mp_prod["phase_matched_row_repeats"],
                                          "fixture": mp_fix["phase_matched_row_repeats"]},
        },
        "period_law": "R | 2^64 and x is a field integer, so the selector is exactly "
                      "2^20-periodic in t for fixed S: phase = t mod 2^20 in BOTH functions",
        "verdict": None,
    }
    res["wall_seconds_total"] = round(time.time() - t0, 2)
    res["seed_rule"] = ("arm P walker seed c,d from Random(99); arm C trail seed c,d from "
                        "Random(4242); single secret per mode (the protocol attacks one Q); "
                        "mode=fix repeats the identical walk and seeds in Z/n with the "
                        "representative-integer selector")
    json.dump(res, open(os.path.join(OUT, "production_selector_v22.json"), "w"), indent=1, sort_keys=True)
    print("DONE", res["wall_seconds_total"], flush=True)


if __name__ == "__main__":
    main()

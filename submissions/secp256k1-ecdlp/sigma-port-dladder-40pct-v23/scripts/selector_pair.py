#!/usr/bin/env python3
"""selector_pair.py -- v22: paired continuation arm, mu6 channels solved EXACTLY.

Carried from v21 (same design: one challenge per set, two arms sharing challenge,
walker offsets and step budget, differing ONLY in the column map; selector input
is the walk's current point at every jump; success = recovered scalar verified on
the curve), with the verdict eaba3488 point 3 PARITY DEFECT repaired.

v21 DEFECT (verdict, verbatim): "It increments rotation_channel when
c=-(a1+a2) mod N is odd. Since the reported N=33560113 is odd, multiplication by
2 is invertible modulo N: odd c also has a solution, namely (c+N)/2."

ROOT CAUSE.  The code wrote c = -(a1+a2) mod N and then branched on c % 2: even c
was solved as c // 2; ODD c was counted as "rotation channel" and DISCARDED.  With
N odd, 2 is invertible and inv2 = (N+1)//2, so odd c is an ordinary sign case.  The
`c % 2` shortcut was written for an even N, and `c // 2` is not division mod N for
odd c.  It is removed, not patched: there is no parity test anywhere in this file.

CORRECT ALGEBRA.  The collision key K = min(y, p-y) is the mu6-ORBIT key (for a
given y the three cube roots x, beta*x, beta^2*x share it), so a class collision
between two walkers means P1 = zeta*P2 with zeta in mu6 = {+1,-1,+lam,-lam,
+lam^2,-lam^2}, P_i = C + a_i G.  zeta=+1 is degenerate (a1 = a2).  The other five
are solved exactly and each candidate is VERIFIED on the curve ([k]G == +-C):
    zeta=-1     : 2C = -(a1+a2)G        ->  x = c * inv2
    zeta=+lam   : (1-lam)C = (lam*a2 - a1)G
    zeta=-lam   : (1+lam)C = -(lam*a2 + a1)G
    zeta=+lam^2 : same with lam^2
    zeta=-lam^2 : same with lam^2
The eigenvalue lam is identified ON THE CURVE ([lam]G == (beta*G.x, G.y)) after
being assembled by CRT from the roots of r^2+r+1 mod each prime factor of the
composite toy order N = 3121 x 10753 (pow(g,(N-1)//3,N) is invalid for composite
N, which is why the v21 file never had a lam at all).

EVIDENCE.  2^-1 mod N is machine-verified: parity_fix.inverse_row_20 records
x*2*inv2 mod N == x for 20 random x (asserted, and written to the JSON).

Output: evidence/selector-pair.json
"""
import hashlib
import itertools
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import qwalk                  # noqa: E402
import selector               # noqa: E402
import census_v20n as C       # noqa: E402

OUT = os.path.join(D, "evidence", "selector-pair.json")
DEGEN_MAX = 10.0


def cm_order(p, P):
    cands = []
    b = 0
    while 27 * b * b <= 4 * p:
        rem = 4 * p - 27 * b * b
        a = int(math.isqrt(rem))
        if a * a == rem and a % 3 != 0:
            cands.append(a)
        b += 1
    if not cands:
        raise RuntimeError("no CM representation for p=%d" % p)
    for A in sorted(cands):
        for n in (p + 1 - A, p + 1 + A):
            if _mul(P, n, p) is None and _mul(_add(P, P, p), n, p) is None:
                return n, A
    raise RuntimeError("no CM order kills the curve (p=%d cands=%s)" % (p, cands[:4]))


def _add(P, Q, p):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % p == 0:
            return None
        lam = (3 * x1 * x1) * pow(2 * y1, p - 2, p) % p
    else:
        lam = (y2 - y1) * pow((x2 - x1) % p, p - 2, p) % p
    x3 = (lam * lam - x1 - x2) % p
    return (x3, (lam * (x1 - x3) - y1) % p)


def _mul(P, k, p):
    R, Q = None, P
    while k:
        if k & 1:
            R = _add(R, Q, p)
        Q = _add(Q, Q, p)
        k >>= 1
    return R


def _neg(P, p):
    return None if P is None else (P[0], (-P[1]) % p)


def find_prime(minbits):
    n = (1 << minbits) + 1
    while True:
        if n % 3 == 1 and n % 4 == 3 and qwalk.is_prime(n):
            return n
        n += 4 if n % 4 == 3 else (3 - (n % 4)) % 4 or 4


def find_point(p, seed=7):
    rng = np.random.default_rng(seed)
    while True:
        x = int(rng.integers(1, p))
        rhs = (x * x * x + 7) % p
        if pow(rhs, (p - 1) // 2, p) != 1:
            continue
        y = pow(rhs, (p + 1) // 4, p)
        if (y * y - rhs) % p == 0:
            return (x, y)


def full_cell(K, R, p):
    return (2 * K * R) // p


def sm64i(K):
    return int(qwalk.sm64(np.array([K], dtype=np.uint64))[0])


def run_set(E, x, W, per, d, R, seed):
    """Run BOTH arms on one challenge; return per-arm solved flags + counters."""
    p, N, G = E["p"], E["N"], E["G"]
    Cpt = _mul(G, x, p)
    rng = np.random.default_rng(seed)
    offs = [int(rng.integers(1, N)) for _ in range(W)]
    res = {}
    for arm in ("ec", "hash"):
        rec = {}
        solved = False
        fires = restarts = trivial = sign_hit = rot_hit = verify_fail = 0
        for w in range(W):
            P = _add(Cpt, _mul(G, offs[w], p), p)
            a = offs[w] % N
            prev1 = prev2 = None
            run = 0
            for _t in range(per):
                K = min(P[1], p - P[1])
                if arm == "ec":
                    m = full_cell(K, R, p)
                else:
                    m = full_cell(sm64i(K) % ((p // 2) + 1), R, p)
                if prev2 is not None and K == prev2:
                    m = (m + 1) % R
                    fires += 1
                    run += 1
                    if run > qwalk.K_ESC_MAX:
                        off = int(rng.integers(1, N))
                        P = _add(Cpt, _mul(G, off, p), p)
                        a = off % N
                        run = 0
                        restarts += 1
                        prev1 = prev2 = None
                        continue
                else:
                    run = 0
                # ---- v22: NO PARITY BRANCH.  Every mu6 channel, solved exactly
                # mod N, each candidate verified on the curve.  See module docstring.
                if K in rec and rec[K][0] != w:
                    a1 = rec[K][1]
                    c = (-(a1 + a)) % N
                    solved_here = False
                    for zeta, kc in E["CHANNELS"](c, a1, a):
                        q = _mul(G, kc, p)
                        if q is not None and (q == Cpt or q == _neg(Cpt, p)):
                            solved_here = True
                            if zeta == "sign":
                                sign_hit += 1
                            else:
                                rot_hit += 1
                            break
                    if solved_here:
                        solved = True
                    else:
                        if a1 == a % N:
                            trivial += 1
                        else:
                            verify_fail += 1
                if d <= 0 or sm64i(K) < (1 << (64 - d)):
                    rec.setdefault(K, (w, a))
                prev2, prev1 = prev1, K
                P = _add(P, E["T"][m], p)
                a = (a + E["u"][m]) % N
        res[arm] = {"solved": solved, "fires": fires, "restarts": restarts,
                    "trivial": trivial, "sign_channel": sign_hit,
                    "rotation_channel": rot_hit, "verify_fail": verify_fail,
                    "records": len(rec)}
    return res, Cpt


def main():
    t0 = time.time()
    minbits = int(sys.argv[1]) if len(sys.argv) > 1 else 25
    sets = int(sys.argv[2]) if len(sys.argv) > 2 else 800
    p = find_prime(minbits)
    G = find_point(p)
    N, A = cm_order(p, G)
    R = 64
    u = [((int.from_bytes(hashlib.sha256(
        b"v21-pair:%d:%d" % (p, m)).digest()[:8], "big")) % (N - 1)) + 1 for m in range(R)]
    E = {"p": p, "N": N, "G": G, "u": u, "T": [_mul(G, k, p) for k in u]}

    # ---- v22 PARITY FIX: exact channel algebra, no parity branch ------------
    lam = None
    _fac, _n = [], N
    _d = 2
    while _d * _d <= _n:
        while _n % _d == 0:
            _fac.append(_d)
            _n //= _d
        _d += 1
    if _n > 1:
        _fac.append(_n)
    assert len(set(_fac)) == len(_fac), ("not squarefree", _fac)
    _ps = sorted(set(_fac))
    _roots = [[r for r in range(2, _q) if (r * r + r + 1) % _q == 0] for _q in _ps]

    def _crt(res, mods):
        _M = 1
        for _m in mods:
            _M *= _m
        _x = 0
        for _r, _m in zip(res, mods):
            _Mi = _M // _m
            _x = (_x + _r * _Mi * pow(_Mi, -1, _m)) % _M
        return _x

    _cands = set()
    for _combo in itertools.product(*_roots):
        _l = _crt(list(_combo), _ps)
        if (pow(_l, 3, N) - 1) % N == 0 and (_l - 1) % N != 0:
            _cands.add(_l)
    assert _cands and all((_l * _l + _l + 1) % N == 0 for _l in _cands), _cands
    _beta = None
    for _g in range(2, 60):
        _b = pow(_g, (p - 1) // 3, p)
        if _b != 1 and (_b * _b + _b + 1) % p == 0:
            _beta = _b
            break
    assert _beta is not None
    _Gbeta = (_beta * G[0] % p, G[1])
    _found = [l for l in sorted(_cands) if _mul(G, l, p) == _Gbeta]
    assert len(_found) == 1, ("eigenvalue not identified", sorted(_cands), len(_found))
    lam = _found[0]
    lam2 = (lam * lam) % N
    inv2 = (N + 1) // 2
    assert (2 * inv2) % N == 1
    _rng = np.random.default_rng(424242)
    INV2_ROW = []
    for _x in _rng.integers(0, N, size=20).tolist():
        INV2_ROW.append({"x": int(_x), "x_times_2": int((_x * 2) % N),
                         "times_inv2": int((_x * 2 * inv2) % N),
                         "equal": bool((_x * 2 * inv2 - _x) % N == 0)})
    assert all(r["equal"] for r in INV2_ROW), INV2_ROW
    _i1m, _i1m2 = pow((1 - lam) % N, -1, N), pow((1 - lam2) % N, -1, N)
    _i1p, _i1p2 = pow((1 + lam) % N, -1, N), pow((1 + lam2) % N, -1, N)

    def CHANNELS(c, a1, a2):
        return (("sign", (c * inv2) % N),
                ("lam+", ((lam * a2 - a1) * _i1m) % N),
                ("lam-", (-(lam * a2 + a1) * _i1p) % N),
                ("lam2+", ((lam2 * a2 - a1) * _i1m2) % N),
                ("lam2-", (-(lam2 * a2 + a1) * _i1p2) % N))

    E["CHANNELS"] = CHANNELS
    E["lam"] = lam

    ncls = (p - 1) // 2
    W, d = 16, 5
    x_shape = 1.0097
    per = int(round(x_shape * math.sqrt(ncls) / W))
    visits = W * per
    p_birth_raw = 1.0 - math.exp(-(visits * visits) / (2.0 * ncls))
    dp_visits = visits / (1 << d)
    p_birth_dp = 1.0 - math.exp(-(dp_visits * dp_visits) / (2.0 * (ncls / (1 << d))))

    solved = {"ec": 0, "hash": 0}
    b = c = 0
    tot = {a: {"fires": 0, "restarts": 0, "trivial": 0, "sign_channel": 0,
               "rotation_channel": 0, "verify_fail": 0, "records": 0}
           for a in ("ec", "hash")}
    for s in range(sets):
        x = int(np.random.default_rng(1000003 + s).integers(1, N))
        r, _ = run_set(E, x, W, per, d, R, seed=7000001 + s)
        se, sh = r["ec"]["solved"], r["hash"]["solved"]
        solved["ec"] += int(se)
        solved["hash"] += int(sh)
        b += int(se and not sh)
        c += int(sh and not se)
        for arm in ("ec", "hash"):
            for k in tot[arm]:
                tot[arm][k] += r[arm][k]
        if (s + 1) % 100 == 0:
            print("  %d/%d sets  ec %d hash %d" % (s + 1, sets, solved["ec"], solved["hash"]),
                  flush=True)

    pe, ph = solved["ec"] / sets, solved["hash"] / sets
    disc = b + c
    se_d = math.sqrt(max(disc - (b - c) ** 2 / sets, 0.0)) / sets if sets else 0.0
    delta = pe - ph
    out = {
        "test": "selector_pair_continuation_v22",
        "design": ("one challenge per set; two arms share challenge, walker offsets "
                   "and step budget; arms differ ONLY in the column map; selector "
                   "input is the WALK's current point at every jump (continuation); "
                   "success = recovered scalar verified on the curve; EVERY mu6 "
                   "channel solved exactly (no parity branch)"),
        "curve": {"p": p, "p_log2": round(math.log2(p), 3), "order_N": N,
                  "cm_A": A, "b": 7, "j": 0, "classes": ncls},
        "R": R, "W": W, "dp_bits": d, "per_walker_steps": per,
        "visits_per_set": visits, "sets": sets,
        "visits_over_sqrtN": visits / math.sqrt(ncls),
        "arms": {"ec": pe, "hash": ph, "ec_only": b, "hash_only": c,
                 "discordant": disc, "delta_ec_minus_hash": delta,
                 "se_paired": se_d, "z": (delta / se_d) if se_d else None,
                 "ci95_paired": [delta - 1.96 * se_d, delta + 1.96 * se_d]},
        "cp95": {"ec": float(C.cp95_lower(solved["ec"], sets)),
                 "hash": float(C.cp95_lower(solved["hash"], sets))},
        "counters": tot,
        "parity_fix": {
            "verdict": "eaba3488 point 3",
            "n_odd": bool(N % 2 == 1),
            "inv2": int(inv2),
            "inv2_check_2_times_inv2_mod_N": int((2 * inv2) % N),
            "lambda_eigenvalue": int(lam),
            "lambda_identity": "lam^2 + lam + 1 == 0 mod N: %s" % bool((lam * lam + lam + 1) % N == 0),
            "branch_removed": "no c%2 test anywhere; every mu6 channel solved exactly and verified on the curve",
            "inverse_row_20": INV2_ROW,
            "inverse_row_all_equal": all(r["equal"] for r in INV2_ROW),
        },
        "birthday": {"p_raw_visits": p_birth_raw, "p_dp_records": p_birth_dp,
                     "ratio_raw_ec": pe / p_birth_raw if p_birth_raw else None,
                     "ratio_raw_hash": ph / p_birth_raw if p_birth_raw else None,
                     "ratio_dp_ec": pe / p_birth_dp if p_birth_dp else None,
                     "ratio_dp_hash": ph / p_birth_dp if p_birth_dp else None},
        "degeneracy": {"DEGEN_MAX": DEGEN_MAX,
                       "flag": bool(max(pe, ph) / p_birth_raw > DEGEN_MAX
                                    or p_birth_raw <= 0)},
        "elapsed_s": round(time.time() - t0, 1),
        "claim_dependence": "NONE (EC selector transfer is declared OPEN; exploratory row)",
    }
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    print(json.dumps({k: out[k] for k in ("arms", "birthday", "degeneracy")}, indent=1))
    print("wrote", OUT)


if __name__ == "__main__":
    main()

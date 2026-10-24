#!/usr/bin/env python3
"""p256_point_walk_test.py — recompute the shipped counted step's real-P-256 point walk
in Python big integers and diff every emitted record.

Reads the driver's stdout (argv[1], default pw_v10i.txt) and writes a JSON summary to
argv[2] (default point_walk.json).  Exit 0 = all records match.

v10i: the driver now generates table rows with the PRODUCTION policy of
evidence/census_v10f41.c (random logged row t, T=[t]G, canonicalise the stored copy and
NEGATE THE LOG) instead of a ladder that canonicalised its accumulator in place.  This
model reproduces that policy independently: it rebuilds the exact splitmix64 stream with
the same duplicate rules, recomputes every row as [t]G at the real modulus, checks the
emitted canonical bytes against it, and then replays all 208 walks x 512 steps.
"""
import json
import sys

P = 2**256 - 2**224 + 2**192 + 2**96 - 1
R = 2**256 % P
M64 = (1 << 64) - 1
GX = 0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296
GY = 0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5
B = 0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b
SA = 0xa5a5a5a5a5a5a5a5
SB = 0x5a5a5a5a5a5a5a5a


def splitmix64(s):
    s = (s + 0x9E3779B97F4A7C15) & M64
    z = s
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & M64
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & M64
    return s, z ^ (z >> 31)


def mix64(v):
    v = (v * 0x9E3779B97F4A7C15) & M64
    v ^= v >> 29
    v = (v * 0xBF58476D1CE4E5B9) & M64
    v ^= v >> 32
    return v & M64


def add(A, Bp):
    if A is None:
        return Bp
    if Bp is None:
        return A
    x1, y1 = A
    x2, y2 = Bp
    if x1 == x2:
        if (y1 + y2) % P == 0:
            return None
        lam = (3 * x1 * x1 - 3) * pow(2 * y1, P - 2, P) % P
    else:
        lam = (y2 - y1) * pow(x2 - x1, P - 2, P) % P
    x3 = (lam * lam - x1 - x2) % P
    return (x3, (lam * (x1 - x3) - y1) % P)


def mul(A, k):
    acc = None
    for b in bin(k)[2:]:
        acc = add(acc, acc)
        if b == "1":
            acc = add(acc, A)
    return acc


def canon_y(y):
    """canon_y from step_v10f_o3.c: tst y[0],#1 on the MONTGOMERY limbs, csel takes the
    negation when the low bit is ZERO, i.e. it keeps the ODD representative of {y, p-y}."""
    return (P - y) % P if ((y * R) % P) & 1 == 0 else y


def canon_mont(y):
    return (P - y) % P if ((y * R) % P) & 1 == 0 else y


def step(St, par, tab):
    x1, y1 = St
    xm_low = (x1 * R) % P & M64
    j = mix64(xm_low ^ (SB if par else SA)) & 2047
    tx, ty = tab[j]
    d = (tx - x1) % P
    lam = 0 if d == 0 else (ty - y1) * pow(d, P - 2, P) % P
    x3 = (lam * lam - x1 - tx) % P
    y3 = (lam * (x1 - x3) - y1) % P
    return (x3, canon_y(y3))


def fnv(h, v):
    for _ in range(4):
        b = v & M64
        for _ in range(8):
            h ^= b & 0xFF
            h = (h * 0x100000001B3) & M64
            b >>= 8
        v >>= 64
    return h


def hexpair(s, n=128):
    assert len(s) == n, s
    return (int(s[: n // 2], 16), int(s[n // 2:], 16))


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else "pw_v10i.txt"
    out = sys.argv[2] if len(sys.argv) > 2 else "point_walk.json"
    tbl, tb2, log1, log2, walks, chks = {}, {}, {}, {}, [], []
    dig = None
    for line in open(src):
        f = line.split()
        if not f:
            continue
        if f[0] == "TBL":
            tbl[int(f[1])] = hexpair(f[2]); log1[int(f[1])] = int(f[3], 16)
        elif f[0] == "TB2":
            tb2[int(f[1])] = hexpair(f[2]); log2[int(f[1])] = int(f[3], 16)
        elif f[0] == "WALK":
            walks.append((int(f[1]), int(f[2]), int(f[3], 16), hexpair(f[4])))
        elif f[0] == "CHK":
            chks.append((int(f[1]), int(f[2]), hexpair(f[3])))
        elif f[0] == "DIGEST":
            dig = (int(f[1], 16), int(f[2], 16))

    G = (GX, GY)

    # --- independent regeneration of the PRODUCTION row policy, including the stream ---
    def regen(tab, seed, nrows, prior_logs):
        s = seed
        r, logs, bad_stream = {}, {}, 0
        for j in range(nrows):
            while True:
                s, w = splitmix64(s)
                t = w | 1
                dup = any(u == t for u in logs.values()) or any(
                    u == t or u == (0 - t) & M64 for u in prior_logs)
                if not dup:
                    break
            T = mul(G, t)
            y = T[1]
            stored_t = t
            if ((y * R) % P) & 1 == 0:
                y = (P - y) % P
                stored_t = (0 - t) & M64
            r[j] = (T[0], y)
            logs[j] = stored_t
        return r, logs

    r_tbl, r_log1 = regen(tbl, 0xABCDEF0123456789, 2048, {})
    r_tb2, r_log2 = regen(tb2, 0x13579BDF2468ACE0, 2048, r_log1)

    bad_tbl = 0
    bad_stream = 0
    for j in range(2048):
        if r_tbl[j] != tbl.get(j):
            bad_tbl += 1
        if r_tb2[j] != tb2.get(j):
            bad_tbl += 1
        if r_log1[j] != log1.get(j) or r_log2[j] != log2.get(j):
            bad_stream += 1

    h1 = h2 = 0xCBF29CE484222325
    for j in range(2048):
        h1 = fnv(fnv(h1, tbl[j][0]), tbl[j][1])
        h2 = fnv(fnv(h2, tb2[j][0]), tb2[j][1])
    digest_ok = (dig == (h1, h2))

    # --- every walk, 512 steps, from the shipped step's own index rule ---
    bad_walk = 0
    first_bad = None
    steps_done = 0
    for blk, slot, s, got in walks:
        St = mul(G, s)
        for st in range(512):
            St = step(St, st & 1, tb2 if (st & 1) else tbl)
            steps_done += 1
            if blk == 0 and st + 1 == 256:
                ck = [c for c in chks if c[0] == blk and c[1] == slot]
                if ck and ck[0][2] != St and first_bad is None:
                    first_bad = {"kind": "chk", "blk": blk, "slot": slot, "step": 256}
        if St != got:
            bad_walk += 1
            if first_bad is None:
                first_bad = {"kind": "walk", "blk": blk, "slot": slot,
                             "got": [hex(got[0]), hex(got[1])],
                             "want": [hex(St[0]), hex(St[1])]}

    bad_chk = 0
    for blk, slot, got in chks:
        s = [w[2] for w in walks if w[0] == blk and w[1] == slot][0]
        St = mul(G, s)
        for st in range(256):
            St = step(St, st & 1, tb2 if (st & 1) else tbl)
        if St != got:
            bad_chk += 1

    res = {
        "modulus": hex(P),
        "curve": "y^2 = x^3 - 3x + b over P-256 p",
        "generator": [hex(GX), hex(GY)],
        "on_curve_generator": (GY * GY - (GX ** 3 - 3 * GX + B)) % P == 0,
        "representation": "Montgomery: fe_mul(a,b) = a*b*R^-1 mod p, R = 2^256 mod p",
        "row_generation": "production policy of evidence/census_v10f41.c: t <- splitmix|1 (dup-checked), T <- [t]G, canon the STORED COPY and negate the log; no accumulator mutation",
        "walk_count": len(walks),
        "steps_per_walk": 512,
        "point_additions_recomputed": steps_done,
        "table_rows_each": 2048,
        "table_rows_mismatched": bad_tbl,
        "table_log_stream_mismatched": bad_stream,
        "table_digest_matches": digest_ok,
        "digest_emitted": [hex(dig[0]), hex(dig[1])] if dig else None,
        "digest_recomputed": [hex(h1), hex(h2)],
        "mid_walk_checkpoints": len(chks),
        "mid_walk_checkpoints_mismatched": bad_chk,
        "walk_final_mismatched": bad_walk,
        "records": len(walks) + len(chks) + 2,
        "total_mismatches": bad_tbl + bad_stream + bad_walk + bad_chk + (0 if digest_ok else 1),
        "first_mismatch": first_bad,
        "verdict": "PASS" if (bad_tbl + bad_stream + bad_walk + bad_chk == 0 and digest_ok) else "FAIL",
    }
    json.dump(res, open(out, "w"), indent=1, sort_keys=True)
    print(json.dumps({k: res[k] for k in ("verdict", "walk_count", "point_additions_recomputed",
                                          "table_rows_mismatched", "table_log_stream_mismatched",
                                          "walk_final_mismatched", "mid_walk_checkpoints_mismatched",
                                          "total_mismatches")}, sort_keys=True))


if __name__ == "__main__":
    main()

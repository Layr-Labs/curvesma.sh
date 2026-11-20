#!/usr/bin/env python3
"""
hx_z12.py — Standalone on-demand H-R & H-X verification over Z[zeta_12] mod Q = 1,000,033
with 4 dense Schirokauer columns and block-diagonal X (plus 1-block pad negative control).
Pure Python stdlib only. Input: {"D": 64..512, "seed": int, "m": 16..64, "n": 8..24}.
"""
import json
import math
import random
import sys

Q = 1000033
CHI = 385966


def main():
    raw = sys.stdin.read().strip()
    req = json.loads(raw) if raw else {}
    D = int(req.get("D", 120))
    seed = int(req.get("seed", 1))
    m = int(req.get("m", 24))
    n = int(req.get("n", 12))
    if not (32 <= D <= 512 and 8 <= m <= 64 and 6 <= n <= 24 and 0 <= seed <= 1000000):
        raise ValueError("parameter out of range")

    assert (pow(CHI, 4, Q) - pow(CHI, 2, Q) + 1) % Q == 0
    u12 = [1, CHI, pow(CHI, 2, Q), pow(CHI, 3, Q), (pow(CHI, 2, Q) - 1) % Q, (pow(CHI, 3, Q) - CHI) % Q]
    u12 += [(-u) % Q for u in u12]
    rng = random.Random(seed + 1000003 * D)
    s_sm = 4
    N = D - s_sm
    w_row = min(N - 2, max(12, D // 4))
    M_rows = [[(c, rng.choice(u12)) for c in rng.sample(range(N), w_row)] for _ in range(D)]
    S_cols = [[rng.randrange(1, Q) for _ in range(D)] for _ in range(s_sm)]

    def matvec_B(v):
        return [sum(val * v[c] for c, val in row) % Q for row in M_rows]

    def rank_vecs(vecs, dim):
        basis = {}
        for v in vecs:
            cur = [x % Q for x in v]
            for i in range(dim):
                if cur[i]:
                    if i in basis:
                        b = basis[i]; f = cur[i]
                        for j in range(i, dim): cur[j] = (cur[j] - f * b[j]) % Q
                    else:
                        inv = pow(cur[i], -1, Q)
                        basis[i] = [(x * inv) % Q for x in cur]
                        break
        return len(basis)

    Y = [list(c) for c in S_cols] + [matvec_B([rng.randrange(1, Q) for _ in range(D)]) for _ in range(n - s_sm)]
    d = math.ceil(D / n) + 1
    e = math.ceil(D / m) + 1
    cols_C = []
    cur_Y = [matvec_B(y) for y in Y]
    for _ in range(d):
        cols_C.extend(cur_Y)
        cur_Y = [matvec_B(y) for y in cur_Y]
    rk_Cd = rank_vecs(cols_C, D)

    A_cols = []
    for c in range(N):
        col = [0] * D
        for r in range(D):
            for cc, vv in M_rows[r]:
                if cc == c: col[r] = vv
        A_cols.append(col)
    A_cols.extend(S_cols)
    rk_A = rank_vecs(A_cols, D)

    def obs_rk(one_block):
        blocks = [[] for _ in range(m)]
        perm = list(range(N)); rng.shuffle(perm)
        for idx, c in enumerate(perm): blocks[idx % m].append((c, rng.randrange(1, Q)))
        for k in range(s_sm): blocks[0 if one_block else k].append((N + k, rng.randrange(1, Q)))
        vecs = []
        for col in A_cols:
            out = []; cur = col
            for _ in range(e):
                for blk in blocks: out.append(sum(w * cur[idx] for idx, w in blk) % Q)
                cur = matvec_B(cur)
            vecs.append(out)
        return rank_vecs(vecs, m * e)

    rk_OeA = obs_rk(False)
    rk_neg = obs_rk(True)
    ok = (rk_A == D and rk_Cd == N and rk_OeA == D and rk_neg == D - (s_sm - 1))
    print(json.dumps(dict(
        checks_ok=ok, D=D, N=N, s_sm=s_sm, m=m, n=n, d=d, e=e, seed=seed,
        rank_A=rk_A, rank_Cd=rk_Cd, rank_OeA=rk_OeA, neg_1block_pad_rank=rk_neg,
    ), separators=(",", ":")))
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

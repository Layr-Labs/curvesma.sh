#!/usr/bin/env python3
"""grid_bsgs_toy.py — small-scope confirmation of the FULL v5 pipeline and
constructive refutations of the two rejected mechanisms. Toy scope only,
exactly as declared: additive model of G1 = (Z/r_t, +) with P0 = 1,
r_t = 1094657 = 2^10 * 1069 + 1 (r_t and q_t = 1069 both prime),
d_t = 2^10. In this model a "point" is an integer mod r_t, [c]X = c*X mod
r_t, and point addition is integer addition: exactly the same
additive-preserving encoding a -> [a]P0 as the production G1. Nothing
here is a production complexity proof.

  T1  full v5 pipeline: geometric baby chain X_{v+1} = zeta * X_v
      (multiply the RUNNING POINT by the public scalar zeta), geometric
      giant chain by zeta^{-m}, exact-coverage grid, phase-2 delta grid,
      verification gate. 2000/2000 exact recoveries (success = 1).
  T2  v3-style ARITHMETIC chains by single additions (baby table
      [v*zeta]P0, giant A_d - u*[m*zeta]P0): the match identity
      "k2 = u*m + v" is FALSE in the point domain — a hit requires the
      ENCODED scalar zeta^k2 to land on the additive lattice
      {u*m*zeta + v}, a ~m*M/r_t fraction, and even then it does not
      reveal k2. Empirical hit histogram printed.
  T3  v4-style exponentiation of the encoded scalar: [c * zeta^k2]P0
      (legal scalar mult) vs [ (zeta^k2)^c ]P0 (the multiplicative PH
      object) disagree on essentially every k2 — the encoding preserves
      addition, not multiplication.
"""
import random, json

RT = 1094657            # = 2^10 * 1069 + 1
DT = 1 << 10
QT = (RT - 1) // DT     # 1069

def is_prime(n):
    if n < 2: return False
    d_, s_ = n - 1, 0
    while d_ % 2 == 0: d_ //= 2; s_ += 1
    for a in (2,3,5,7,11,13):
        x = pow(a, d_, n)
        if x in (1, n-1): continue
        for _ in range(s_-1):
            x = x*x % n
            if x == n-1: break
        else: return False
    return True
assert is_prime(RT) and is_prime(QT)

def find_xi(n):
    rn = n - 1
    fac, t, p = set(), rn, 2
    while p*p <= t:
        if t % p == 0:
            fac.add(p)
            while t % p == 0: t //= p
        p += 1
    if t > 1: fac.add(t)
    g = 2
    while any(pow(g, rn//pp, n) == 1 for pp in fac): g += 1
    return g

XI    = find_xi(RT)
ZETA  = pow(XI, DT, RT)    # exact order QT in F_rt*
DELTA = pow(XI, QT, RT)    # exact order DT

def ceil_div(a, b): return -(-a // b)
MB = 33                        # toy baby width m = ceil(sqrt(1069)) = 33
MG = ceil_div(QT, MB)          # 33 giants
M2 = 32                        # M2^2 = DT = 1024

def sm(c, X): return (c % RT) * X % RT        # [c]X in the additive model

def full_pipeline(s):
    T  = s % RT                              # T = [s]P0
    Ad = sm(pow(s, DT, RT), 1)               # A_d = [zeta^k2]P0
    stab = {}
    X = 1                                    # X_0 = P0
    for v in range(MB):
        stab[X] = v                          # store [zeta^v]P0
        X = sm(ZETA, X)                      # geometric chain step
    Y = Ad
    k2 = None
    for u in range(MG):
        if Y in stab:
            t = u*MB + stab[Y]               # candidate exponent
            for cand in (t, t - QT):
                if 0 <= cand < QT: k2 = cand
            break
        Y = sm(pow(ZETA, -MB, RT), Y)        # giant chain step, public scalar
    if k2 is None:
        raise AssertionError("phase1 miss — coverage broken")
    Tp = sm(pow(XI, (-k2) % (RT-1), RT), T)  # T' = [xi^-k2]T = [delta^k1]P0
    tab = {}
    X = 1
    for v in range(M2):
        tab[X] = v
        X = sm(DELTA, X)
    Y = Tp
    k1 = None
    for u in range(M2):
        if Y in tab:
            k1 = u*M2 + tab[Y]
            break
        Y = sm(pow(DELTA, -M2, RT), Y)
    if k1 is None:
        raise AssertionError("phase2 miss — grid broken")
    s_hat = pow(XI, QT*k1 + k2, RT)
    if sm(s_hat, 1) != T:
        raise AssertionError("verify gate rejected the unique true s")
    return s_hat

def t1(trials=2000):
    rng = random.Random(0xC0FFEE)
    for _ in range(trials):
        s = rng.randrange(1, RT)
        assert full_pipeline(s) == s, s
    print(f"T1 PASS: full v5 pipeline exact on {trials}/{trials} random trials "
          f"(zero phase misses, zero verify rejects — grid success=1)")

def t2(trials=4000):
    """v3-style arithmetic grid: probe whether the encoded scalar zeta^k2
    equals some u*m*zeta + v as INTEGERS mod r_t (what single-addition
    chains actually compare), NOT whether k2 = u*m+v."""
    rng = random.Random(1)
    hits = wrong_digit = 0
    for _ in range(trials):
        k2 = rng.randrange(0, QT)
        e = pow(ZETA, k2, RT)
        match = None
        for u in range(MG):
            for v in range(MB):
                if (u*MB*ZETA + v) % RT == e:
                    match = (u, v)
            if match: break
        if match:
            hits += 1
            if match[0]*MB + match[1] != k2:
                wrong_digit += 1
    frac = hits/trials
    print(f"T2: v3-style arithmetic-lattice probe hit {hits}/{trials} = {frac:.4f} "
          f"of challenges (expected ~m*M/r_t = {MB*MG/RT:.4f}); of hits, "
          f"{wrong_digit} yielded a WRONG k2 digit")
    assert frac < 0.05, "v3-style model unexpectedly survived — check construction"

def t3(trials=4000):
    """v4-style: the legitimate scalar mult [c*zeta^k2]P0 vs the
    multiplicative-PH object [(zeta^k2)^c]P0, c = (r_t-1)/2 (order-2
    PH projector in F_rt*)."""
    rng = random.Random(2)
    c = (RT - 1)//2
    agree = 0
    for _ in range(trials):
        k2 = rng.randrange(0, QT)
        z = pow(ZETA, k2, RT)
        legal   = c * z % RT                  # [c]A_d : scalar mult by c
        ph_obj  = pow(z, c, RT)               # A_d^(z->z^c): exponentiate encoded scalar
        agree += (legal == ph_obj)
    print(f"T3: [c*zeta^k2]P0 == [zeta^{{k2*c}}]P0 on only {agree}/{trials} probes "
          f"(encoding preserves addition, not multiplication: v4's PH digit was "
          f"reading the left side as if it were the right)")
    assert agree <= trials//50

if __name__ == "__main__":
    t1()
    t2()
    t3()
    print(json.dumps({
        "scope": "toy additive model, r_t=1094657=2^10*1069+1, d_t=2^10, q_t=1069",
        "T1_trials": 2000, "T1_pipeline_exact": True,
        "T2_v3_arithmetic_grid_hit_fraction_threshold": 0.05,
        "T3_v4_scalar_vs_exponent_agreement_threshold": 0.02,
        "production_complexity_claimed_by_this_file": False}))

"""Supporting evidence for proof.md (standard library only; JSON on stdin, one JSON object on stdout).

Four parts, all recomputed from the declared inputs:

1. ledger    - every cost term of proof.md section 5 from exact integers and rationals, the exact score
               (rounded UP to 10^-6), the preprocessing and memory bounds (rounded UP), the success bound
               (rounded DOWN). This checks the arithmetic of the ledger, not the per-operation budgets.
2. planner   - the offline pairing planner of proof.md section 3 run on real chi-class lists of a toy prime
               at the production degree parameter GMAX/d^2 = 16: pair fraction f, cycle rejections, root
               demotions, maximum dependency depth and stack (toy evidence for H4).
3. toy       - the whole attack on a toy j = 0 curve: class lists, planner, pair walks with ONE inversion
               per P+T / P-T pair (the formula is executed, not simulated), canonical-y table, collision,
               d-th root, mu_d search. Success rate against 1 - exp(-K); every recovered scalar is checked;
               field multiplications per element are counted.
4. coverage  - exact |U*B*G^-1| / m for random sets, the smallest integers in H and the product construction
               (toy evidence for H1/H2, as in the earlier accepted package).

Toy results are evidence at their own scale only, never production cost or success proofs.
Input: {"K", "f_min", "logX2_census", "census_classes", "atk_bits", "atk_logd", "atk_logX2", "atk_logG",
        "trials", "cov_bits", "cov_logd", "seed"}.
"""
import json
import math
import random
import sys
import time
from decimal import Decimal, getcontext
from fractions import Fraction as F

getcontext().prec = 80
C = 2000
r = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001
d = 1 << 20
m = (r - 1) // d
assert (r - 1) % d == 0 and m % 6 == 0 and (r - 1) % (1 << 32) == 0 and (r - 1) % (1 << 33) != 0


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def log2_up(x, places):
    """Smallest multiple of 10^-places that is >= log2(x), x a positive Fraction; exact via Decimal."""
    v = (Decimal(x.numerator).ln() - Decimal(x.denominator).ln()) / Decimal(2).ln()
    q = Decimal(10) ** places
    cand = (v * q).to_integral_value(rounding="ROUND_CEILING")
    # margin: the 80-digit evaluation is far more precise than 10^-places; assert we are not at a boundary
    assert abs(v * q - cand) > Decimal("1e-30") or v * q == cand
    return F(int(cand), 10 ** places)


def log2_float(x):
    return float((Decimal(x.numerator).ln() - Decimal(x.denominator).ln()) / Decimal(2).ln())


# ----------------------------------------------------------------------------------------------- 1. ledger
def ledger(K, f_min):
    K, f = F(K), F(f_min)
    # ---- parameters (proof.md section 2)
    X2 = 1 << 76                 # class lists hold b <= X2
    LIST = X2 // d               # average list length 2^56 (H2)
    A_LO, A_LEN = 1 << 76, 1 << 70   # base primes a in (A_LO, A_LO + A_LEN]
    GMAX = 1 << 44               # per-prime table T[g] = [g]Q, 1 <= g <= GMAX
    LANES = 1 << 16              # primes in flight; one field inversion per round of LANES lane operations
    # ---- per-operation instruction budgets (proof.md section 5)
    M6 = C                       # single op: 6M under the reference convention
    M9 = 3 * C // 2              # pair op: 9M (one shared inversion share 3M + two 3M formulas)
    INV = 10                     # per lane op: share of the one Fp inversion per round (383000/65536 = 5.85) and
                                 # of the two-level batch exchange and wake-ups between the 16 nodes of a group
    WALK1 = 150                  # single op bookkeeping (as in the earlier accepted package)
    PAIR2 = 2 * 170              # pair op bookkeeping for two points (listed 325, budget 340)
    KEY = 40                     # canonical y, per point
    HASH = 180                   # insert (baby) or lookup (giant), per point (listed 173)
    c_pair = F(M9 + PAIR2 + 2 * (KEY + HASH), 2 * C)       # per pair-derived point
    c_single = F(M6 + INV + WALK1 + KEY + HASH, C)         # per root point
    c = f * c_pair + (1 - f) * c_single                    # per element, with pair fraction >= f
    assert c_pair == F(3780, 4000) and c_single == F(2380, 2000)
    # ---- sizes
    N = math.isqrt(math.ceil(K * m / 6)) + 1
    assert 6 * N * N >= K * m
    OPS = (1 - f / 2) * N            # plan entries executed per side, at most N (1 - f/2)
    n_a = 2 * (2 * N) // LIST + 2        # primes used on both sides, charged at twice the expected number (H2)
    pi_lo = F(A_LEN, 78 * 693148 // 10**6) * F(9, 10)   # > 0.9 * A_LEN / ln(A_LO + A_LEN), primes available
    assert n_a < pi_lo / 2               # at least 2x slack in the base-prime interval
    PI_X2 = F(126, 100) * X2 / F(int(76 * math.log(2) * 10**6), 10**6)      # pi(X2) < 1.26 X2 / ln X2
    SCALAR77, SCALAR255 = 1000, 3000     # units per scalar multiplication (<= 77-bit, 255-bit scalar)
    terms = {
        "baby side: N elements at c = f c_pair + (1-f) c_single": N * c,
        "giant side: N elements at c": N * c,
        "per-prime tables T[g], g <= GMAX: n_a x GMAX additions at 2000 + 60 + 10 instr": F(n_a * GMAX * (M6 + 60 + INV), C),
        "plan streaming: 9 instr per entry at the plan node, both sides": F(2 * OPS * 9, C),
        "per-prime scalar multiplications (Q_a, first element, T[2]): 3 x 1000 units": F(n_a * 3 * SCALAR77),
        "chi of primes <= X2 (100 units each)": 100 * PI_X2,
        "chi of primes in the base interval (100 units each)": 100 * F(A_LEN, 52) ,
        "chi sieve and class lists to X2 (200 instr per b)": F(X2 * 200, C),
        "prime sieve of the base interval (64 instr per integer) and of primes <= 2^39": F((A_LEN + (1 << 39)) * 64, C),
        "planning: per element 2^24 midpoint tests x 16 + matching 2^19 + emission 64 instr": F(X2 * ((1 << 24) * 16 + (1 << 19) + 64), C),
        "candidates: up to 2^20, six 255-bit scalar multiplications each": F((1 << 20) * 6 * SCALAR255),
        "d-th root in F_r (<= 2^20 F_r multiplications of 200 instr)": F((1 << 20) * 200, C),
        "mu_d search (2^20 + 1 scalar multiplications)": F(((1 << 20) + 1) * SCALAR255),
        "input decoding and subgroup checks": F(4 * SCALAR255),
    }
    W = sum(terms.values())
    pre_keys = [k for k in terms if k.startswith(("baby side", "per-prime", "chi", "prime sieve", "planning", "plan streaming"))]
    pre = sum(terms[k] for k in pre_keys)
    score = log2_up(W, 6)
    pre_log2 = log2_up(pre, 4)
    # memory (bytes): hash table 2N slots x 48 B (load <= 1/2), chi array 2.5 B x X2, lists + plans 26 B x X2,
    # tables and slots for LANES primes in flight, primes of the base interval 16 B each
    NODES = (1 << 60) + (1 << 20) + 16     # table nodes, plan nodes, walk nodes
    mem = (2 * N * 48 + F(5, 2) * X2 + 27 * X2 + LANES * (GMAX * 96 + (1 << 15) * 96) + 16 * F(A_LEN, 52)
           + NODES * ((1 << 12) * 49152 + (1 << 12) * 4 * 49152 + (1 << 12)))   # buffers, rings, code and state
    mem_log2 = log2_up(mem, 4)
    # success: 1 - exp(-6N^2/m) rounded down, minus the structural loss (2^-14, section 4), minus 2^-20 (caps)
    Kexact = F(6 * N * N, m)
    s = 1 - math.exp(-float(Kexact))
    s_low = math.floor((s - 2 ** -14 - 2 ** -20) * 10**4) / 10**4
    SUCCESS_CLAIM = 0.399
    assert 0.39 <= SUCCESS_CLAIM <= s_low
    return {
        "K": str(K), "f_min": str(f), "c_pair_point": str(c_pair), "c_single_point": str(c_single),
        "c_per_element": float(c), "N_log2": round(log2_float(F(N)), 6), "primes_charged_log2": round(math.log2(n_a), 3),
        "terms_log2": {k: round(log2_float(v), 3) for k, v in terms.items()},
        "W_log2": round(log2_float(W), 7), "W_units_integer_ceiling": str(math.ceil(W)),
        "score": {"numerator": str(score.numerator), "denominator": str(score.denominator)},
        "score_decimal": float(score),
        "preprocessing_log2_upper": float(pre_log2), "memory_log2_bytes_upper": float(mem_log2),
        "heuristic_success_1_minus_exp": round(s, 6), "success_lower_bound_after_losses": s_low,
        "success_claimed": SUCCESS_CLAIM,
        "reference_single_walk_same_budgets_log2": round(log2_float(2 * N * c_single + (W - 2 * N * c)), 4),
    }


# ----------------------------------------------------------------------------------------- number theory
def is_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    e, s = n - 1, 0
    while e % 2 == 0:
        e //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, e, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def primes_in(lo, hi):
    return [q for q in range(lo, hi + 1) if is_prime(q)]


def sqrt_mod(a, p):
    a %= p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = next(z for z in range(2, p) if pow(z, (p - 1) // 2, p) == p - 1)
    mm, c, t, x = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p
            i += 1
        b = pow(c, 1 << (mm - i - 1), p)
        mm, c, t, x = i, b * b % p, t * b * b % p, x * b % p
    return x


def toy_prime(bits, dd, rng):
    while True:
        rr = rng.randrange(1 << (bits - 1), 1 << bits) // (6 * dd) * (6 * dd) + 1
        if is_prime(rr) and ((rr - 1) // dd) % 6 == 0 and pow(2, (rr - 1) // dd, rr) != 1:
            return rr


def character(rr, dd):
    """chi: integers -> Z/d, chi(x) = log_zeta(x^m), kernel H; plus the generator used."""
    mm = (rr - 1) // dd
    gen = next(g for g in range(2, 10**6) if all(pow(g, (rr - 1) // q, rr) != 1 for q in (2, 3)))
    zeta = pow(gen, mm, rr)
    zlog = {pow(zeta, i, rr): i for i in range(dd)}
    return (lambda x: zlog[pow(x, mm, rr)]), gen


def class_lists(rr, dd, X2):
    chi, _ = character(rr, dd)
    classes = [[] for _ in range(dd)]
    for b in range(1, X2 + 1):
        classes[chi(b)].append(b)
    return classes


# ---------------------------------------------------------------------------------------------- 2. planner
def plan(L, GMAX, max_depth=1 << 14):
    """Offline pairing plan for one sorted class list L (proof.md section 3).

    Pair (i, j) with midpoint k: L[k] - L[i] = L[j] - L[k] = g, 1 <= g <= GMAX. Children get parent k
    (in-degree 1); roots are computed by one addition from the previous root (gap <= GMAX enforced).
    Matching: greedy by increasing degree, then one pass of length-3 augmenting paths; every pair is
    rejected if it would close a dependency cycle or exceed the depth cap. Emission: sweep in list order,
    each root is emitted (single op) and then its whole subtree (pair ops, depth first).
    """
    n = len(L)
    pos = {v: i for i, v in enumerate(L)}
    adj = [[] for _ in range(n)]
    n_edges = 0
    for i in range(n):
        vi = L[i]
        lim = vi + 2 * GMAX
        j = i + 1
        while j < n and L[j] <= lim:
            vj = L[j]
            if (vi + vj) & 1 == 0:
                k = pos.get((vi + vj) >> 1)
                if k is not None:
                    adj[i].append((j, k))
                    adj[j].append((i, k))
                    n_edges += 1
            j += 1
    parent, partner, depth = [-1] * n, [-1] * n, [0] * n
    height = [0] * n           # upper bound on the height of the subtree below each element
    children = [[] for _ in range(n)]
    stats = {"cycle_rejections": 0, "depth_rejections": 0, "demoted": 0}

    def set_depth(v):
        stack = [v]
        while stack:
            x = stack.pop()
            depth[x] = depth[parent[x]] + 1 if parent[x] != -1 else 0
            for (a, b) in children[x]:
                stack.append(a)
                stack.append(b)

    def creates_cycle(i, j, k):
        x = k
        while x != -1:
            if x == i or x == j:
                return True
            x = parent[x]
        return False

    def try_match(i, j, k, force=False):
        """Attach i and j below k. The depth check covers the subtrees already hanging below i and j;
        heights are upper bounds (never lowered on unmatch), so the check only over-rejects."""
        if parent[i] != -1 or parent[j] != -1:
            return False
        if not force:
            if depth[k] + 1 + max(height[i], height[j]) > max_depth:
                stats["depth_rejections"] += 1
                return False
            if creates_cycle(i, j, k):
                stats["cycle_rejections"] += 1
                return False
        parent[i] = parent[j] = k
        partner[i], partner[j] = j, i
        children[k].append((min(i, j), max(i, j)))
        set_depth(i)
        set_depth(j)
        h, x = 1 + max(height[i], height[j]), k
        while x != -1 and height[x] < h:
            height[x] = h
            h += 1
            x = parent[x]
        return True

    def unmatch(u):
        k, b = parent[u], partner[u]
        children[k].remove((min(u, b), max(u, b)))
        parent[u] = parent[b] = -1
        partner[u] = partner[b] = -1
        set_depth(u)
        set_depth(b)

    for v in sorted(range(n), key=lambda v: len(adj[v])):
        if parent[v] != -1:
            continue
        for (w, k) in adj[v]:
            if parent[w] == -1 and try_match(v, w, k):
                break
    for v in range(n):                      # one pass of length-3 augmenting paths
        if parent[v] != -1:
            continue
        done = False
        for (w, k) in adj[v]:
            if parent[w] == -1:
                if try_match(v, w, k):
                    done = True
                    break
                continue
            w2, k0 = partner[w], parent[w]
            for (u, k2) in adj[w2]:
                if u != v and u != w and parent[u] == -1:
                    unmatch(w)
                    if try_match(v, w, k) and try_match(w2, u, k2):
                        done = True
                        break
                    for (a, b) in ((v, w), (w2, u)):
                        if parent[a] != -1 and partner[a] == b:
                            unmatch(a)
                    assert try_match(w, w2, k0, force=True)   # restore the accepted pair
            if done:
                break
    if parent[0] != -1:                     # the first element is computed by a scalar multiplication
        unmatch(0)
        stats["demoted"] += 1
    last = L[0]
    for v in range(1, n):                   # consecutive roots at most GMAX apart
        if L[v] - last > GMAX:
            u = v - 1
            if parent[u] != -1:
                unmatch(u)
                stats["demoted"] += 1
            last = L[u]
            assert L[v] - last <= GMAX
        if parent[v] == -1:
            last = L[v]
    emitted = [False] * n
    ops = []
    max_stack = n_pairs = n_single = 0
    last_root = None
    for v in range(n):
        if emitted[v] or parent[v] != -1:
            continue
        if last_root is None:
            ops.append((0, v))
        else:
            ops.append((1, v, L[v] - L[last_root]))
            n_single += 1
        emitted[v] = True
        last_root = v
        stack = [(v, 0)]
        while stack:
            max_stack = max(max_stack, len(stack))
            x, idx = stack[-1]
            if idx < len(children[x]):
                stack[-1] = (x, idx + 1)
                a, b = children[x][idx]
                g = L[x] - L[a]
                assert L[b] - L[x] == g and 0 < g <= GMAX
                ops.append((2, x, a, b, g))
                n_pairs += 1
                emitted[a] = emitted[b] = True
                stack.append((a, 0))
                stack.append((b, 0))
            else:
                stack.pop()
    assert all(emitted), "cycle left in the dependency graph"
    assert max(depth) <= max_depth and max_stack <= 2 * max_depth + 1
    assert n_single + 1 + 2 * n_pairs == n
    stats.update({"n": n, "avg_degree": round(2 * n_edges / n, 2), "pairs": n_pairs, "roots": n_single + 1,
                  "pair_fraction": round(2 * n_pairs / n, 5), "max_depth": max(depth), "max_stack": max_stack})
    return ops, stats


def planner_census(cfg):
    rng = random.Random(cfg["seed"])
    dd = 16
    rr = toy_prime(28, dd, rng)
    X2 = 1 << cfg["logX2_census"]
    cl = class_lists(rr, dd, X2)
    GMAX = 16 * dd * dd          # production degree parameter: GMAX / d^2 = 2^44 / 2^40 = 16
    out = {"r": rr, "d": dd, "X2": X2, "GMAX": GMAX, "GMAX_over_d2": 16, "classes": []}
    for c in range(cfg["census_classes"]):
        _, st = plan(cl[c], GMAX)
        out["classes"].append(st)
    out["min_pair_fraction"] = min(s["pair_fraction"] for s in out["classes"])
    out["scope"] = "real chi-class lists of a toy prime; evidence for H4 at this scale only"
    return out


# --------------------------------------------------------------------------------------------- 3. toy attack
class Curve:
    def __init__(self, p, c):
        self.p, self.c = p, c
        self.inv = self.muls = 0

    def add(self, P, Q):
        if P is None:
            return Q
        if Q is None:
            return P
        p = self.p
        (x1, y1), (x2, y2) = P, Q
        if x1 == x2:
            if (y1 + y2) % p == 0:
                return None
            lam = 3 * x1 * x1 * pow(2 * y1, -1, p) % p
        else:
            lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
        x3 = (lam * lam - x1 - x2) % p
        return (x3, (lam * (x1 - x3) - y1) % p)

    def single(self, P, T):
        """Root step P + T, incomplete affine (x1 != x2 asserted): one inversion share and 3M (6M in all)."""
        assert P is not None and T is not None and P[0] != T[0]
        self.inv += 1
        self.muls += 6
        return self.add(P, T)

    def pair(self, P, T):
        """P + T and P - T from ONE inversion of (xT - xP): 3M (share) + 2 x 3M = 9M."""
        p = self.p
        (x1, y1), (x2, y2) = P, T
        assert x1 != x2
        inv = pow(x2 - x1, -1, p)
        self.inv += 1
        self.muls += 9
        lp = (y2 - y1) * inv % p
        lm = (-y2 - y1) * inv % p
        x3p = (lp * lp - x1 - x2) % p
        x3m = (lm * lm - x1 - x2) % p
        return (x3p, (lp * (x1 - x3p) - y1) % p), (x3m, (lm * (x1 - x3m) - y1) % p)

    def mul(self, k, P):
        R = None
        while k:
            if k & 1:
                R = self.add(R, P)
            P = self.add(P, P)
            k >>= 1
        return R


def find_curve(bits, dd, rng):
    """p = 1 mod 3, y^2 = x^3 + c; prime group order r with 6d | r-1 and 2 not in H (as in production)."""
    while True:
        p = rng.randrange(1 << (bits - 1), 1 << bits) | 1
        if p % 3 != 1 or not is_prime(p):
            continue
        for c in range(1, 8):
            E = Curve(p, c)
            while True:
                x = rng.randrange(p)
                y = sqrt_mod((x ** 3 + c) % p, p)
                if y is not None:
                    P = (x, y)
                    break
            lo, hi = p + 1 - 2 * math.isqrt(p) - 2, p + 1 + 2 * math.isqrt(p) + 2
            w = math.isqrt(hi - lo) + 1
            baby, R = {}, None
            for j in range(w):
                if R is not None:
                    baby.setdefault(R, j)
                R = E.add(R, P)
            step = E.mul(w, P)
            step = (step[0], (-step[1]) % p) if step else None
            Gk, cands = E.mul(lo, P), []
            for i in range(w + 1):
                if Gk is None:
                    cands.append(lo + i * w)
                elif (Gk[0], (-Gk[1]) % p) in baby:
                    cands.append(lo + i * w + baby[(Gk[0], (-Gk[1]) % p)])
                Gk = E.add(Gk, step)
            cands = sorted({N for N in cands if lo <= N <= hi and E.mul(N, P) is None})
            if len(cands) == 1 and is_prime(cands[0]) and (cands[0] - 1) % (6 * dd) == 0 \
                    and pow(2, (cands[0] - 1) // dd, cands[0]) != 1:
                return E, cands[0], P


def dth_root(t, dd, rr):
    k = dd.bit_length() - 1
    for j in range(k, 0, -1):
        y = sqrt_mod(t, rr)
        e = (rr - 1) >> (j - 1)
        t = next(z for z in (y, rr - y) if j == 1 or pow(z, e, rr) == 1)
    return t


def toy_attack(cfg):
    bits, logd, trials, K = cfg["atk_bits"], cfg["atk_logd"], cfg["trials"], float(cfg["K"])
    rng = random.Random(cfg["seed"])
    dd = 1 << logd
    E, rr, P = find_curve(bits, dd, rng)
    p, mm = E.p, (rr - 1) // dd
    beta = next(z for z in (pow(g, (p - 1) // 3, p) for g in range(2, 100)) if z != 1)
    lam = next(L for L in (pow(g, (rr - 1) // 3, rr) for g in range(2, 100))
               if L != 1 and E.mul(L, P) == (beta * P[0] % p, P[1]))
    U = [1, rr - 1, lam, rr - lam, lam * lam % rr, rr - lam * lam % rr]
    assert all(pow(u, mm, rr) == 1 for u in U)
    chi, gen = character(rr, dd)
    size = math.isqrt(int(K * mm / 6)) + 1
    X2, GMAX = 1 << cfg["atk_logX2"], 1 << cfg["atk_logG"]
    classes = [[] for _ in range(dd)]
    for b in range(1, X2 + 1):
        classes[chi(b)].append(b)
    plans = [plan(L, GMAX) for L in classes]
    pair_frac = sum(2 * st["pairs"] for _, st in plans) / sum(st["n"] for _, st in plans)
    X1 = 2 * X2 + 2
    while True:
        primes = primes_in(X1 // 2 + 1, X1)
        if len(primes) * X2 // dd >= 3 * size:
            break
        X1 *= 2

    def walk(base, plist, count, emit):
        produced = 0
        for a in plist:
            c = (-chi(a)) % dd
            L, (ops, _) = classes[c], plans[c]
            if not L:
                continue
            Qa = E.mul(a, base)                 # per-prime setup, charged separately in the ledger
            tab = {}

            def T(g):                           # toy: table entries by scalar multiplication (production:
                if g not in tab:                # repeated addition, charged in the ledger)
                    tab[g] = E.mul(g, Qa)
                return tab[g]

            pts, root = {}, None
            for op in ops:
                if op[0] == 0:
                    root = E.mul(L[op[1]], Qa)  # first element: scalar multiplication
                    pts[op[1]] = root
                    emit(a * L[op[1]], root)
                    produced += 1
                elif op[0] == 1:
                    root = E.single(root, T(op[2]))
                    pts[op[1]] = root
                    emit(a * L[op[1]], root)
                    produced += 1
                else:
                    _, x, i, j, g = op
                    Pp, Pm = E.pair(pts[x], T(g))
                    pts[j], pts[i] = Pp, Pm
                    emit(a * L[j], Pp)
                    emit(a * L[i], Pm)
                    produced += 2
                if produced >= count:
                    return
        raise RuntimeError("not enough primes")

    succ = inv_total = mul_total = 0
    for _ in range(trials):
        s = rng.randrange(1, rr)
        A1, Ad = E.mul(s, P), E.mul(pow(s, dd, rr), P)
        E.inv = E.muls = 0
        rng.shuffle(primes)
        table, found = {}, []
        canon = lambda y: min(y, p - y)         # same key for the six points [u]Q, u in U

        checks = [0]

        def baby(v, Pt):
            assert pow(v, mm, rr) == 1
            if checks[0] < 8:                   # spot-check the walk against scalar multiplication
                checks[0] += 1
                assert E.mul(v, P) == Pt
            table[canon(Pt[1])] = v

        def giant(v, Pt):
            assert pow(v, mm, rr) == 1
            if not found and canon(Pt[1]) in table:
                found.append((table[canon(Pt[1])], v))

        walk(P, primes[: len(primes) // 2], size, baby)
        walk(Ad, primes[len(primes) // 2:], size, giant)
        inv_total += E.inv
        mul_total += E.muls
        if not found:
            continue
        b, g = found[0]
        t = next(c for c in (u * b * pow(g, -1, rr) % rr for u in U) if E.mul(c, P) == Ad)
        s0 = dth_root(t, dd, rr)
        w = pow(gen, (rr - 1) // dd, rr)
        sol = next(s0 * pow(w, i, rr) % rr for i in range(dd) if E.mul(s0 * pow(w, i, rr) % rr, P) == A1)
        assert sol == s, "wrong scalar"
        succ += 1
    Ke = 6 * size * size / mm
    return {"p_bits": round(math.log2(p), 2), "r": rr, "d": dd, "set_size_per_side": size, "X2": X2, "GMAX": GMAX,
            "GMAX_over_d2": GMAX / dd / dd, "list_lengths_first_8": [len(L) for L in classes][:8],
            "plan_pair_fraction": round(pair_frac, 4),
            "K_effective": round(Ke, 4), "predicted_success": round(1 - math.exp(-Ke), 4),
            "trials": trials, "successes": succ, "success_rate": round(succ / trials, 4),
            "std_error": round(math.sqrt(Ke * 0 + (1 - math.exp(-Ke)) * math.exp(-Ke) / trials), 4),
            "inversions_per_trial": inv_total / trials, "elements_per_trial": 2 * size,
            "field_mults_per_element": round(mul_total / trials / (2 * size), 3),
            "field_mults_per_element_reference_walk": 6.0,
            "all_recovered_scalars_correct": True,
            "scope": "toy curve only; validates the formulas, the planner and the recovery, not production costs"}


# ---------------------------------------------------------------------------------------------- 4. coverage
def coverage(cfg):
    bits, logd, K, seed = cfg["cov_bits"], cfg["cov_logd"], float(cfg["K"]), cfg["seed"]
    rng = random.Random(seed)
    dd = 1 << logd
    rr = toy_prime(bits, dd, rng)
    mm = (rr - 1) // dd
    chi, gen = character(rr, dd)
    lam = pow(gen, (rr - 1) // 3, rr)
    U = [1, rr - 1, lam, rr - lam, lam * lam % rr, rr - lam * lam % rr]
    assert all(pow(u, mm, rr) == 1 for u in U)
    n = math.isqrt(int(K * mm / 6)) + 1

    def cov(B, G):
        Ginv = [pow(g, -1, rr) for g in G]
        bitmap = bytearray(rr // 8 + 1)
        for u in U:
            for b in B:
                ub = u * b % rr
                for gi in Ginv:
                    v = ub * gi % rr
                    bitmap[v >> 3] |= 1 << (v & 7)
        return int.from_bytes(bitmap, "little").bit_count() / mm

    out = {"r_bits": round(math.log2(rr), 2), "d": dd, "set_size": n, "K_effective": round(6 * n * n / mm, 4),
           "predicted": round(1 - math.exp(-6 * n * n / mm), 4)}
    rnd = lambda: [pow(rng.randrange(1, rr), dd, rr) for _ in range(n)]
    out["random_sets"] = round(cov(rnd(), rnd()), 4)
    small, x = [], 1
    while len(small) < 2 * n:
        if pow(x, mm, rr) == 1:
            small.append(x)
        x += 1
    out["smallest_integers_in_H"] = round(cov(small[0::2], small[1::2]), 4)
    out["product_construction"] = []
    for ratio in (64, 256):
        X2 = ratio * dd
        classes = {}
        for b in range(1, X2 + 1):
            classes.setdefault(chi(b), []).append(b)
        X1 = 2 * X2 + 2
        while True:
            Pl = primes_in(X1 // 2 + 1, X1)
            if len(Pl) * X2 // dd >= 3 * n:
                break
            X1 *= 2
        rng.shuffle(Pl)

        def build(plist):
            S = []
            for a in plist:
                for b in classes.get((-chi(a)) % dd, []):
                    S.append(a * b)
                    if len(S) == n:
                        return S
            raise RuntimeError("not enough primes")

        B, G = build(Pl[: len(Pl) // 2]), build(Pl[len(Pl) // 2:])
        assert all(pow(v, mm, rr) == 1 for v in B + G) and len(set(B)) == n and len(set(G)) == n
        out["product_construction"].append({"X2_over_d": ratio, "coverage": round(cov(B, G), 4)})
    out["scope"] = "toy prime only; evidence for H1/H2 at this scale, not a production proof"
    return out


if __name__ == "__main__":
    cfg = json.load(sys.stdin)
    t0 = time.time()
    res = {"ledger": ledger(cfg["K"], cfg["f_min"])}
    log("ledger done", round(time.time() - t0, 1))
    res["planner_census"] = planner_census(cfg)
    log("planner done", round(time.time() - t0, 1))
    res["toy_attack"] = toy_attack(cfg)
    log("toy done", round(time.time() - t0, 1))
    res["coverage"] = coverage(cfg)
    log("coverage done", round(time.time() - t0, 1))
    res["elapsed_seconds"] = round(time.time() - t0, 1)
    print(json.dumps(res, indent=1))

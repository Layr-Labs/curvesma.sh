"""Floyd r-adding walk on BLS12-381 G1: reproducible support for the work
ledger and the success bound in proof.md. Standard library only.

The program does two things.

1. Recomputes the instruction total W of proof.md from the BLS12-381 profile
   constants and the pinned cost-model constants, and checks that the claimed
   score rational 137562769/1000000 is a ceiling of log2(W / 2000).

2. Executes the exact walk function of the production algorithm on toy inputs:
   additive groups Z/nZ of prime order at three sizes, and the short
   Weierstrass curve y^2 = x^3 + 4 over F_9613 inside a prime-order subgroup
   (cofactor 3, like G1). It draws
   80 (seed, secret) pairs per group from the same mix used by the attack,
   runs the Floyd r-adding walk with the production 3*ceil(sqrt(n)) iteration
   budget, and reports recoveries and iteration-length statistics.

Exit status is nonzero iff the recomputed W disagrees with proof.md, the score
is not a ceiling of log2(W), a returned scalar fails the group-law check, or a
draw that found a collision returned the wrong residue. The integer-group rows
simulate the point values (the jumps are formed from the known secret); the
curve row executes real affine curve arithmetic including the inversions whose
cost dominates W.
"""
import json
import sys
from math import ceil, isqrt, log2

# --- BLS12-381 public profile constants (tracks/.../target-profile.json) ---
P_MOD = 0x1A0111EA397FE69A4B1BA7B6434BACD764774B84F38512BF6730D2A0F6B0F6241EABFFFEB153FFFFB9FEFFFFFFFFAAAB
R_ORDER = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001

# --- Scoring constants ---
FMUL = 383   # measured blst_fp_mul instructions, scripts/reference_costs/measurement.json (bls_field diff)
FADD = 128   # analytical ceiling for a 6-limb modular addition/subtraction, proof.md
C_UNIT = 2000  # G1 target-group-operation equivalent, instructions (classical-group-work-v2)

CLAIMED_NUMERATOR = 137562769
CLAIMED_DENOMINATOR = 1000000
EXPECTED_W = 514694484759148014223122347276647627415912096

MASK64 = (1 << 64) - 1


def mix64(w):
    """Four-step 64-bit bijective mix used by the partition and the stream."""
    w &= MASK64
    w ^= w >> 30
    w = (w * 0xBF58476D1CE4E5B9) & MASK64
    w ^= w >> 27
    w = (w * 0x94D049BB133111EB) & MASK64
    w ^= w >> 31
    return w


def block512(seed, tag):
    """Deterministic 512-bit expansion: eight splitmix64 words keyed by (seed, tag)."""
    state = (seed ^ ((tag + 1) * 0x9E3779B97F4A7C15)) & MASK64
    out = 0
    for _ in range(8):
        state = (state + 0x9E3779B97F4A7C15) & MASK64
        out = (out << 64) | mix64(state)
    return out


def make_jumps(seed, n):
    """32 jump coefficients: c_j = block512 mod (n-1) + 1 (nonzero), d_j = block512 mod n."""
    return [
        (block512(seed, 2 * j) % (n - 1) + 1, block512(seed, 2 * j + 1) % n)
        for j in range(32)
    ]


def partition_index(word, seed):
    """h = ((mix64(word) * ((2*seed+1) mod 2^64)) >> 32) & 31."""
    return ((mix64(word) * ((2 * seed + 1) & MASK64)) >> 32) & 31


# --- group wrappers -----------------------------------------------------------

class IntGroup:
    """Additive group Z/nZ; the walk element is the residue itself."""

    def __init__(self, n):
        self.n = n
        self.zero = 0

    def add(self, a, b):
        return (a + b) % self.n

    def equal(self, a, b):
        return a == b

    def word(self, a):
        return a & MASK64

    def scalar_mul(self, k, elem):
        return (k * elem) % self.n

    def jump_point(self, c, d, base, target):
        return (c * base + d * target) % self.n

    def check(self, x, base, target):
        return (x * base) % self.n == target


class ToyCurve:
    """Affine short Weierstrass curve y^2 = x^3 + b over F_q (cofactored toy)."""

    def __init__(self, q, b=4):
        self.q = q
        self.b = b
        self.zero = None  # point at infinity

    def on_curve(self, pt):
        if pt is None:
            return True
        x, y = pt
        return (y * y - (x * x * x + self.b)) % self.q == 0

    def add(self, a, b):
        if a is None:
            return b
        if b is None:
            return a
        x1, y1 = a
        x2, y2 = b
        if x1 == x2:
            if (y1 + y2) % self.q == 0:
                return None
            lam = (3 * x1 * x1) * pow(2 * y1, -1, self.q) % self.q
        else:
            lam = (y2 - y1) * pow((x2 - x1) % self.q, -1, self.q) % self.q
        x3 = (lam * lam - x1 - x2) % self.q
        y3 = (lam * (x1 - x3) - y1) % self.q
        return (x3, y3)

    def equal(self, a, b):
        return a == b

    def word(self, a):
        """Canonical affine word; 0 at the identity (proof.md, sec. canonical word)."""
        if a is None:
            return 0
        x, y = a
        return (x ^ (y << 1)) & MASK64

    def scalar_mul(self, k, pt):
        out = None
        add = pt
        while k:
            if k & 1:
                out = self.add(out, add)
            add = self.add(add, add)
            k >>= 1
        return out

    def jump_point(self, c, d, base, target):
        return self.add(self.scalar_mul(c, base), self.scalar_mul(d, target))

    def check(self, x, base, target):
        return self.scalar_mul(x, base) == target


# --- the walk -----------------------------------------------------------------

def adding_walk_solve(group, base, target, jumps, seed, budget_iters):
    """Floyd iteration on the production r-adding walk.

    States are (point, alpha, beta) with point = [alpha]base + [beta]target.
    On the first equality with (beta_h - beta_t) nonzero mod n the secret is
    x = (alpha_t - alpha_h) * (beta_h - beta_t)^(-1) mod n, then verified.
    A beta collision is walked through. Returns (x_or_None, iterations_used).
    """
    n = group.n

    def advance(state):
        pt, a, b = state
        j = partition_index(group.word(pt), seed)
        jump_pt, c, d = jumps[j]
        return group.add(pt, jump_pt), (a + c) % n, (b + d) % n

    tortoise = (group.zero, 0, 0)
    hare = (group.zero, 0, 0)
    for it in range(1, budget_iters + 1):
        tortoise = advance(tortoise)
        hare = advance(advance(hare))
        if group.equal(tortoise[0], hare[0]):
            denom = (hare[2] - tortoise[2]) % n
            if denom != 0:
                x = ((tortoise[1] - hare[1]) % n) * pow(denom, -1, n) % n
                if group.check(x, base, target):
                    return x, it
                return None, it  # recovered residue fails the law: hard error
    return None, budget_iters


def budget(n):
    return 3 * (isqrt(n) + (1 if isqrt(n) ** 2 < n else 0))


def run_group(group, n, samples, master):
    """Draws (seed, secret) from the mix and walks each; returns statistics."""
    succ, lengths, failures = 0, [], 0
    b = budget(n)
    for k in range(samples):
        seed = block512(master, 0x5EED + k) & MASK64
        secret = block512(seed, 0xDEC0DE + k) % n
        jumps = [(group.jump_point(c, d, group_base(group), group_target(group, secret)), c, d)
                 for (c, d) in make_jumps(seed, n)]
        x, it = adding_walk_solve(group, group_base(group), group_target(group, secret), jumps, seed, b)
        if x is None:
            failures += 1
        elif not group.check(x, group_base(group), group_target(group, secret)):
            raise AssertionError("returned scalar fails the group-law check")
        else:
            if x != secret:
                raise AssertionError("returned residue %d != secret %d" % (x, secret))
            succ += 1
            lengths.append(it)
    lengths.sort()
    rn = isqrt(n) + (1 if isqrt(n) ** 2 < n else 0)
    stat = {
        "order": n,
        "successes": succ,
        "failures": failures,
        "iterations_over_sqrt": {
            "mean": round(sum(lengths) / max(1, len(lengths)) / rn, 3) if lengths else None,
            "median": round(lengths[len(lengths) // 2] / rn, 3) if lengths else None,
            "p90": round(lengths[(9 * len(lengths)) // 10] / rn, 3) if lengths else None,
            "max": round(lengths[-1] / rn, 3) if lengths else None,
        },
    }
    return stat


def group_base(group):
    return group._base


def group_target(group, secret):
    return group.scalar_mul(secret, group._base)


def is_prime(n):
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


def toy_curve():
    """y^2 = x^3 + 4 over F_9613: #E = 9561 = 3 * 3187, 3187 prime. Base point
    (6408, 7011) has order 3187 (checked below). A cofactored group, like G1."""
    q, n, base = 9613, 3187, (6408, 7011)
    c = ToyCurve(q)
    assert is_prime(q) and is_prime(n) and c.on_curve(base)
    assert c.scalar_mul(n, base) is None
    assert c.scalar_mul(1, base) == base and c.scalar_mul(n - 1, base) is not None
    c._base = base
    c.n = n
    return c, n


# --- ledger recomputation -----------------------------------------------------

def ledger():
    p = P_MOD
    r = R_ORDER
    assert p.bit_length() == 381 and r.bit_length() == 255
    hw = bin(p - 2).count("1")
    bl = p.bit_length()
    inv_muls = (bl - 1) + (hw - 1)          # 380 squarings + 228 additions = 608 Fp muls
    inv_instr = inv_muls * FMUL + bl * 32   # plus scalar loop/shift overhead
    step = inv_instr + 4 * FMUL + 8 * FADD  # affine adding step: inv + lambda, x3, y3

    sq = isqrt(r)
    ceil_sqrt = sq + (1 if sq * sq < r else 0)
    iters = 3 * ceil_sqrt                  # Floyd iterations
    steps = 9 * ceil_sqrt                  # adding steps (3 per iteration)

    w_steps = steps * step
    w_part = steps * 256                   # word + mix64 + product + mask per step
    w_coef = iters * 1536                  # (alpha, beta) updates per iteration
    w_equal = iters * (8 * FMUL + 8 * FADD)  # affine equality per iteration

    jac_op = 16 * FMUL + 28 * FADD
    setup = (
        2 * 64 * (8 * 40 + 16384)          # two 512-bit derivations: 64 blocks x (8 mix64 + reduction)
        + 64 * 510 * jac_op                # jump points + T check + final [x]P0 verify scalar muls
        + 32 * inv_instr                   # affine normalization of the 32 stored jumps
        + 5_500_000                        # output, resolve, verification overhead
        + 1_000_000                        # stream and loop setup
    )
    reserve = 1 << 40
    w = w_steps + w_part + w_coef + w_equal + setup + reserve
    return {
        "inv_muls": inv_muls,
        "inv_instr": inv_instr,
        "step_instr": step,
        "ceil_sqrt_r": ceil_sqrt,
        "floyd_iterations": iters,
        "adding_steps": steps,
        "w_steps": w_steps,
        "w_part": w_part,
        "w_coef": w_coef,
        "w_equal": w_equal,
        "w_setup": setup,
        "w_reserve": reserve,
        "W": w,
        "units": w / C_UNIT,
        "log2_units": log2(w / C_UNIT),
        "preprocessing_units": ceil((1 << 30) / C_UNIT),
    }


def main():
    req = json.load(sys.stdin)
    samples = req["samples"]
    master = req["master"]
    if not isinstance(samples, int) or not 40 <= samples <= 80:
        raise SystemExit("samples must be an integer in [40, 80]")
    if not isinstance(master, int) or master < 0:
        raise SystemExit("master must be a nonnegative integer")

    led = ledger()
    assert led["inv_muls"] == 608
    assert led["W"] == EXPECTED_W, (led["W"], EXPECTED_W)
    claimed = CLAIMED_NUMERATOR / CLAIMED_DENOMINATOR
    assert claimed >= log2(led["W"] / C_UNIT)
    assert (CLAIMED_NUMERATOR - 1) / CLAIMED_DENOMINATOR < log2(led["W"] / C_UNIT)

    results = []
    for n in (10007, 100003, 1000003):
        assert is_prime(n)
        g = IntGroup(n)
        g._base = 1  # generator of Z/nZ
        results.append(run_group(g, n, samples, master))

    curve, cn = toy_curve()
    results.append(run_group(curve, cn, 40, master))

    total_failures = sum(row["failures"] for row in results)
    print(json.dumps({
        "ledger": {k: (v if isinstance(v, int) else round(v, 6)) for k, v in led.items() if k != "units"},
        "claimed_score": claimed,
        "groups": results,
        "total_failures": total_failures,
    }, indent=2))


if __name__ == "__main__":
    main()

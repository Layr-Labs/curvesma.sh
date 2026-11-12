#!/usr/bin/env python3
"""Recompute the bn254-gt-dlp rho-tower (measured worst-case step) claim and run the production-walk evidence.

Standard library only.  Three parts, each at its stated scope:

  1. ledger -- exact recomputation of every row from the profile order r, the
     measured worst-case per-step price, the capped sampler's worst-case cost
     and its failure probability;
  2. sampler -- the capped SHA-256 rejection sampler of proof.md section 2,
     with a measured acceptance rate against r/2^256 and the observed attempt
     counts, so the cap is exercised rather than assumed;
  3. production walk -- the production 32-adding walk with the shipped 384-byte
     encoding, the shipped mixer (transliterated from rho.c and cross-checked
     against the binary's `hashdump` mode) and the shipped tower
     multiplication, run inside a small-order subgroup of the real BN254 Fp12.
     The subgroup order is the only difference from production; the field, the
     tower, the encoding, the mixer and the walk are the shipped ones.
"""
import hashlib
import json
import math
import os
import random
import sys

R = int('30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001', 16)
P = int('30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47', 16)
C = 20000
C_STEP = 47265            # measured worst-case instructions per walk step
MEASURED_WC = 47264.2405        # max over the 704-configuration Callgrind sweep
K_MAX = 64                        # hard cap: sampler attempts per scalar
N_SCALARS = 66                    # 64 table scalars + a0, b0
SHA_WC = 2048                     # charged worst-case instructions per attempt
N = 2 * math.isqrt(R) + 2
DECLARED_NUM, DECLARED_DEN = 130665165, 1000000
EXTRA_SETUP = 2**36               # 64 exponentiations, table, basis conversion
SAMPLER_WC = N_SCALARS * K_MAX * SHA_WC
EXTRA_VERIFY = 2**32
EXTRA_OUTPUT = 2**20
CTRL_PER_ITER = 4096              # Floyd comparison + loop control per iteration
MASK64 = (1 << 64) - 1
SALT = 0x123456789abcdef
HERE = os.path.dirname(os.path.abspath(__file__))


def ledger():
    walk = 3 * N * C_STEP / C
    ctrl = N * CTRL_PER_ITER / C
    setup = (EXTRA_SETUP + SAMPLER_WC) / C
    ver = EXTRA_VERIFY / C
    out = EXTRA_OUTPUT / C
    total = walk + ctrl + setup + ver + out
    return dict(walk=walk, control=ctrl, setup=setup, sampler=SAMPLER_WC / C,
                verify=ver, output=out, total=total, log2=math.log2(total))


def sampler_failure_bound():
    """P(at least one of the N_SCALARS draws is rejected K_MAX times)."""
    reject = 1 - R / 2**256
    return N_SCALARS * reject**K_MAX


def word_stream(seed):
    """256-bit words w_i = SHA-256(seed || LE64(i))."""
    ctr = 0
    while True:
        yield int.from_bytes(hashlib.sha256(seed + ctr.to_bytes(8, 'little')).digest(), 'little')
        ctr += 1


def sample_scalar(seed):
    """Capped rejection sampler: at most K_MAX attempts.  None = run fails."""
    it = word_stream(seed)
    for k in range(1, K_MAX + 1):
        if next(it) < R:
            return k
    return None


def measure_sampler(n_draws=20000, n_scalars=2000):
    accepted = 0
    it = word_stream(b'sampler-measurement')
    for _ in range(n_draws):
        accepted += 1 if next(it) < R else 0
    worst, fails = 0, 0
    for i in range(n_scalars):
        a = sample_scalar(b'capped-' + i.to_bytes(8, 'little'))
        if a is None:
            fails += 1
        else:
            worst = max(worst, a)
    return dict(draws=n_draws, accepted=accepted, measured_acceptance=accepted / n_draws,
                model_acceptance=R / 2**256, scalars=n_scalars, max_attempts=worst,
                cap=K_MAX, failures=fails)


# ------------------------------------------------ shipped tower arithmetic

def f2mul(x, y):
    a1, b1 = x
    a2, b2 = y
    t0 = a1 * a2 % P
    t1 = b1 * b2 % P
    return (t0 - t1) % P, ((a1 + b1) * (a2 + b2) - t0 - t1) % P


def f2add(x, y):
    return (x[0] + y[0]) % P, (x[1] + y[1]) % P


def f2sub(x, y):
    return (x[0] - y[0]) % P, (x[1] - y[1]) % P


def f2xi(x):
    return (9 * x[0] - x[1]) % P, (x[0] + 9 * x[1]) % P


def f6mul(x, y):
    t0 = f2mul(x[0], y[0])
    t1 = f2mul(x[1], y[1])
    t2 = f2mul(x[2], y[2])
    c0 = f2add(t0, f2xi(f2add(f2mul(x[1], y[2]), f2mul(x[2], y[1]))))
    c1 = f2add(f2add(f2mul(x[0], y[1]), f2mul(x[1], y[0])), f2xi(t2))
    c2 = f2add(f2add(f2mul(x[0], y[2]), t1), f2mul(x[2], y[0]))
    return c0, c1, c2


def f6add(x, y):
    return f2add(x[0], y[0]), f2add(x[1], y[1]), f2add(x[2], y[2])


def f6sub(x, y):
    return f2sub(x[0], y[0]), f2sub(x[1], y[1]), f2sub(x[2], y[2])


def f6mulv(x):                      # v * (c0 + c1 v + c2 v^2) = xi*c2 + c0 v + c1 v^2
    return f2xi(x[2]), x[0], x[1]


def f12mul(x, y):
    t0 = f6mul(x[0], y[0])
    t1 = f6mul(x[1], y[1])
    t2 = f6mul(f6add(x[0], x[1]), f6add(y[0], y[1]))
    return f6add(t0, f6mulv(t1)), f6sub(f6sub(t2, t0), t1)


ZERO2 = (0, 0)
ONE = (((1, 0), ZERO2, ZERO2), (ZERO2, ZERO2, ZERO2))


def f12pow(x, e):
    r = ONE
    while e:
        if e & 1:
            r = f12mul(r, x)
        x = f12mul(x, x)
        e >>= 1
    return r


def encode(y):
    """384-byte encoding of the twelve tower coefficients, exactly as shipped."""
    w = []
    for blk in y:
        for c in blk:
            for v in c:
                m = v * 2**256 % P          # Montgomery form, as in the C
                w += [m & MASK64, (m >> 64) & MASK64,
                      (m >> 128) & MASK64, (m >> 192) & MASK64]
    return w


def hash_index(w, salt=SALT):
    """Transliteration of hash_index() in rho.c: 24 rounds, 48 words, mod 32."""
    h1 = (0x9e3779b97f4a7c15 ^ salt) & MASK64
    h2 = (0xbf58476d1ce4e5b9 + salt) & MASK64
    for i in range(24):
        h1 = ((h1 ^ w[2 * i]) * 0x100000001b3) & MASK64
        h2 = ((h2 ^ w[2 * i + 1]) * 0x9e3779b97f4a7c15) & MASK64
        h2 ^= h2 >> 29
    h1 ^= h1 >> 33
    h2 ^= h2 >> 27
    return (h1 ^ h2) & 31


C_HASHDUMP = [22, 2, 22, 12]     # ./rho hashdump 4


def mixer_crosscheck():
    idx = []
    for e in range(4):
        c = [e * 97 + i + 1 for i in range(12)]
        y = tuple(tuple((c[6 * b + 2 * k], c[6 * b + 2 * k + 1]) for k in range(3))
                  for b in range(2))
        idx.append(hash_index(encode(y)))
    return idx == C_HASHDUMP


# ------------------------- production walk in a small-order BN254 subgroup

Q = 43913                         # prime factor of p^12 - 1
E_ORDER = (P**12 - 1) // Q


def make_generator(seed):
    r = random.Random(seed)
    h = tuple(tuple((r.randrange(P), r.randrange(P)) for _ in range(3)) for _ in range(2))
    g = f12pow(h, E_ORDER)
    return g


def walk_trial(g, seed):
    r = random.Random(seed ^ 0x5eed)
    x = r.randrange(Q)
    z = f12pow(g, x)
    cs = [r.randrange(Q) for _ in range(32)]
    ds = [r.randrange(Q) for _ in range(32)]
    ms = [f12mul(f12pow(g, cs[j]), f12pow(z, ds[j])) for j in range(32)]
    a0, b0 = r.randrange(Q), r.randrange(Q)
    y0 = f12mul(f12pow(g, a0), f12pow(z, b0))
    ty, ta, tb = y0, a0, b0
    hy, ha, hb = y0, a0, b0
    steps = 2 * (math.isqrt(Q) + 1)

    def step(y, a, b):
        j = hash_index(encode(y))
        return f12mul(y, ms[j]), (a + cs[j]) % Q, (b + ds[j]) % Q

    for k in range(1, steps + 1):
        ty, ta, tb = step(ty, ta, tb)
        hy, ha, hb = step(hy, ha, hb)
        hy, ha, hb = step(hy, ha, hb)
        if ty == hy:
            db = (hb - tb) % Q
            if db == 0:
                return 'degenerate', k
            ok = ((ta - ha) * pow(db, -1, Q)) % Q == x
            return ('ok' if ok else 'wrong'), k
    return 'no_collision', steps


def walk_evidence(trials, seed):
    counts, steps = {}, {}
    g = make_generator(0xC0FFEE)      # one public generator, fixed for all trials
    for i in range(trials):
        kind, k = walk_trial(g, seed + i)
        counts[kind] = counts.get(kind, 0) + 1
        steps[kind] = steps.get(kind, 0) + k
    ok = counts.get('ok', 0)
    # Wilson 95% interval for the measured success rate.
    z = 1.96
    ph = ok / trials
    den = 1 + z * z / trials
    ctr = (ph + z * z / (2 * trials)) / den
    half = z * math.sqrt(ph * (1 - ph) / trials + z * z / (4 * trials * trials)) / den
    return dict(trials=trials, successes=ok, rate=ph, wilson95=[ctr - half, ctr + half],
                by_outcome=counts, steps_to_collision=steps,
                max_steps_allowed=2 * (math.isqrt(Q) + 1))


def main():
    args = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
    seed = int(args.get('seed', 20261112))
    trials = int(args.get('trials', 24))
    if not (0 <= seed < 2**32) or not (1 <= trials <= 64):
        raise SystemExit('bad input')

    L = ledger()
    score = DECLARED_NUM / DECLARED_DEN
    fail = sampler_failure_bound()
    model = 1 - math.exp(-(31 / 32) * 2)          # Teske 32-adding collision law
    floor = model - fail
    smp = measure_sampler()

    checks = [
        ('N = 2*ceil(sqrt(r))', N == 2 * (math.isqrt(R) + 1)),
        ('declared score >= log2(W)', score >= L['log2']),
        ('declared score is the rounded-up log2(W)',
         DECLARED_NUM == math.ceil(L['log2'] * 1e6)),
        ('sampler failure bound <= 2^-13.30', fail <= 2**-13.30),
        ('success floor (collision model - sampler failure) >= 0.85', floor >= 0.85),
        ('sampler worst case <= 2^24 instructions', SAMPLER_WC <= 2**24),
        ('Python mixer reproduces the C hashdump indices', mixer_crosscheck()),
        ('measured sampler acceptance within 0.02 of r/2^256',
         abs(smp['measured_acceptance'] - smp['model_acceptance']) < 0.02),
        ('no capped-sampler failure observed', smp['failures'] == 0),
        ('sampler cap is enforced', smp['max_attempts'] <= K_MAX),
        ('per-step charge is the measured worst case, rounded up',
         C_STEP == math.ceil(MEASURED_WC)),
    ]
    wk = walk_evidence(trials, seed)
    checks.append(('production-mixer walk: no wrong recovery', wk['by_outcome'].get('wrong', 0) == 0))


    result = {
        'ledger_log2': L,
        'declared_score': score,
        'per_step_instructions': C_STEP,
        'sampler': dict(k_max=K_MAX, scalars=N_SCALARS,
                        worst_case_instructions=SAMPLER_WC,
                        failure_bound_log2=math.log2(fail), measured=smp),
        'success': dict(model=model, sampler_failure=fail, floor=floor,
                        claimed=0.85, protected_minimum=0.39),
        'production_walk': wk,
        'mixer_crosscheck': C_HASHDUMP,
        'walk_scope': ('32-adding walk, shipped 384-byte encoding, shipped mixer '
                       'and shipped tower multiplication, in the order-%d subgroup '
                       'of the real BN254 Fp12; the subgroup order, not the field, '
                       'encoding, mixer or walk, differs from production' % Q),
        'checks': [{'name': n, 'pass': bool(v)} for n, v in checks],
        'all_pass': all(v for _, v in checks),
    }

    print(json.dumps(result, sort_keys=True, indent=2))


if __name__ == '__main__':
    main()

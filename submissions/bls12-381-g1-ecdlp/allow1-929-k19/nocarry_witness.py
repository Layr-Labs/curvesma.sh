"""K1 value-bound witness (proof.md section 1).

Exact-integer model of the kernel's word-serial Montgomery products: round j forms
T = acc_j + a b_j, q_j = (T mod 2^64) n0 mod 2^64, T += q_j p, acc_{j+1} = T / 2^64. It asserts the two
quantities K1 no longer computes are zero: every pre-shift sum T is below 2^448 (word 7; every partial
sum the listing holds is a sub-sum of T) and every shifted accumulator is below 2^384 (word 6), with
acc_j < a + p throughout; it also checks S = acc_6 = a b R^-1 mod p (up to one p) and S < a b / R + p.
For the square it asserts that twice the cross sum sum_{i<j} a_i a_j 2^(64(i+j)) is below 2^704, i.e.
the doubling of words 1..10 never carries into word 11. Operands: a = 2p - 1, 2p - 2 and random below
2p; b = p - 1 and random below p (first product) or below 1.21p (second product); square inputs up to
2^381 - 1 (the bound proved for mu). Finite evidence for the analytic lemma, not a substitute for it."""
import json, random, sys
from pathlib import Path
root = Path(__file__).resolve().parent
if (root / 'files').is_dir(): root = root / 'files'
p = int(json.loads((root / 'params.json').read_text())['p'], 16)
M = (1 << 64) - 1; N0 = (-pow(p, -1, 1 << 64)) & M; R = 1 << 384


def mont(a, b):
    acc = 0
    for j in range(6):
        bj = (b >> (64 * j)) & M
        t = acc + a * bj
        q = ((t & M) * N0) & M
        t += q * p
        assert t < (1 << 448) and t & M == 0, 'word 7 nonzero'
        acc = t >> 64
        assert acc < a + p and acc < R, 'shifted word 6 nonzero'
    assert (acc - a * b * pow(R, -1, p)) % p == 0 and acc < a * b // R + p + 1
    return acc


def cross_ok(a):
    w = [(a >> (64 * i)) & M for i in range(6)]
    c = sum(w[i] * w[j] << (64 * (i + j)) for i in range(6) for j in range(i + 1, 6))
    assert 2 * c < (1 << 704), 'doubling carries into word 11'
    return c


def run(cfg):
    rng = random.Random(cfg.get('seed', 5)); n = cfg.get('samples', 4000)
    A = [2 * p - 1, 2 * p - 2] + [rng.randrange(2 * p) for _ in range(n)]
    B1 = [p - 1, p - 2] + [rng.randrange(p) for _ in range(n)]
    B2 = [121 * p // 100, 2 * p - 1] + [rng.randrange(121 * p // 100) for _ in range(n)]
    pairs = 0; smax = 0
    for i, a in enumerate(A):
        for b in (B1[i], B2[i], B1[-1 - i]):
            s = mont(a, b); smax = max(smax, s); pairs += 1
    mus = [(1 << 381) - 1, 121 * p // 100, 2 * p * p // R + p] + [rng.randrange(1 << 381) for _ in range(n)]
    for mu in mus: cross_ok(mu)
    mu_max = max(mont(t, iv) for t in (2 * p - 1, 2 * p - 2) for iv in (p - 1, p - 2))
    assert mu_max < (1 << 381) and 2 * p * p // R + p < 121 * p // 100
    return {'status': 'passed', 'product_pairs': pairs, 'square_inputs': len(mus),
            'max_S_over_p': round(smax / p, 6), 'mu_max_bits': mu_max.bit_length(),
            'scope': 'Exact-integer model of the word-serial recurrence at extreme and random operands; '
                     'the deleted words are zero in every sampled round. The bound itself is the lemma of proof.md section 1.'}


if __name__ == '__main__':
    print(json.dumps(run(json.loads(sys.stdin.read() or '{}'))))

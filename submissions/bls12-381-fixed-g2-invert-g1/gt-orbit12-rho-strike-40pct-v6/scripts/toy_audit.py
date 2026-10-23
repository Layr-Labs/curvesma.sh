#!/usr/bin/env python3
"""Fresh v4 audit. Toy only; no restart/ladder code is inherited from v3.
Run from the package root; --write writes the captured JSON and text log.
"""
import argparse
from collections import Counter
import json
from math import isqrt
from pathlib import Path
import random
from tower import Tower, exact_order

SEED = 2026102304
P = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
R = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001


def aligned_exhaustive(r, lam, secrets):
    assert exact_order(lam, r, 12)
    good = zero = pairs = 0
    inv = [0] + [pow(v, -1, r) for v in range(1, r)]
    for x in secrets:
        for e in range(1, r):
            for k in range(12):
                lk = pow(lam, k, r)
                ep = e*pow(lk, -1, r) % r
                for b in range(r):
                    a = (e-b*x) % r
                    for bp in range(r):
                        ap = (ep-bp*x) % r
                        num = (lk*ap-a) % r
                        den = (b-lk*bp) % r
                        pairs += 1
                        if den:
                            recovered = num*inv[den] % r
                            assert recovered == x
                            good += 1
                        else:
                            assert num == 0
                            zero += 1
    return dict(r=r, lambda_order=12, secrets=list(secrets),
                all_nonidentity_elements=True, all_offsets=True,
                all_coefficient_pairs=True, pairs=pairs,
                recovered=good, zero_denominator=zero)


def additive_audit(r, lam, trials, rng):
    L = [pow(lam, k, r) for k in range(12)]
    def meta(e):
        orbit = [e*v % r for v in L]
        key = min(orbit)
        return key, orbit.index(key)
    totals = Counter()
    cycles = Counter()
    for _ in range(trials):
        x = rng.randrange(r)
        aa = [rng.randrange(r) for _ in range(32)]
        bb = [rng.randrange(r) for _ in range(32)]
        def step(state):
            e, a, b = state
            assert e == (a+b*x) % r
            K, k = meta(e)
            j = K % 32
            invk = pow(L[k], -1, r)
            return ((2*e+invk*(aa[j]+bb[j]*x)) % r,
                    (2*a+invk*aa[j]) % r,
                    (2*b+invk*bb[j]) % r)
        a, b = rng.randrange(r), rng.randrange(r)
        s = ((a+b*x) % r, a, b)
        table, seen = {}, {}
        for t in range(4*r):
            e, a, b = s
            K, k = meta(e)
            assert e == (a+b*x) % r
            # Report BEFORE stepping. The stored labels represent this e.
            oracle_orbit = sorted({e*pow(lam, j, r) % r for j in range(12)})
            oracle_dp = oracle_orbit[0] % 2 == 0
            actual_dp = (K & 1) == 0
            totals['dp_predicate_comparisons'] += 1
            totals['dp_mismatch'] += actual_dp != oracle_dp
            # Deliberate fault demonstrates that the mismatch audit can fire.
            totals['injected_dp_faults_detected'] += ((not actual_dp) != oracle_dp)
            if actual_dp:
                fp = (K >> 1) & 3  # Deliberately tiny: exercise false matches.
                for old in table.get(fp, []):
                    e0, a0, b0, K0, k0 = old
                    lk = L[(k-k0) % 12]
                    num = (lk*a-a0) % r
                    den = (b0-lk*b) % r
                    if K0 == K:
                        assert e0 == lk*e % r
                        totals['true_class_pairs'] += 1
                        if den:
                            assert num*pow(den, -1, r) % r == x
                            totals['true_pair_recoveries'] += 1
                        else:
                            assert num == 0
                            totals['true_pair_zero_denominator'] += 1
                    else:
                        totals['fingerprint_alias_pairs'] += 1
                        if den:
                            candidate = num*pow(den, -1, r) % r
                            if candidate != x:
                                # This counter exists ONLY in K0 != K.
                                totals['false_x'] += 1
                            else:
                                totals['alias_candidate_verified'] += 1
                table.setdefault(fp, []).append((e, a, b, K, k))
                totals['records_checked_at_same_step'] += 1
            # Actual independent representative next-step audit, not a
            # predicate compared with itself after already assuming equality.
            e2 = L[rng.randrange(12)]*e % r
            a2, b2 = (e2-b*x) % r, b
            s1, s2 = step(s), step((e2, a2, b2))
            totals['conjugate_step_comparisons'] += 1
            assert meta(s1[0])[0] == meta(s2[0])[0]
            if K in seen:
                mu, length = seen[K], t-seen[K]
                cycles[f'{mu}:{length}'] += 1
                totals['cycle_trials'] += 1
                break
            seen[K] = t
            s = s1
        else:
            raise AssertionError('finite-map cycle was not found')
    assert totals['dp_mismatch'] == 0
    assert totals['injected_dp_faults_detected'] > 0
    assert totals['fingerprint_alias_pairs'] > 0 and totals['false_x'] > 0
    return dict(r=r, lam=lam, classes=(r-1)//12,
                trials=trials, counters=dict(totals),
                cycle_mu_lambda_distribution=dict(sorted(cycles.items())))


def tower_identities(p, r, count, rng):
    f = Tower(p)
    rows = []
    for sample in range(count):
        while True:
            raw = tuple(rng.randrange(p) for _ in range(12))
            g = f.pow(raw, (p**12-1)//r)
            if g != f.one:
                break
        assert f.pow(g, r) == f.one
        explicit, y = [], g
        for k in range(12):
            assert f.frob(g, k) == y
            explicit.append((y[1]+y[7]) % p)
            y = f.pow(y, p)
        assert y == g and tuple(explicit) == f.keys(g)
        assert f.frob(g, 6) == f.pow(g, r-1)
        K, kap, rep = f.metadata(g)
        for k in range(12):
            assert f.metadata(f.frob(g, k))[2] == rep
        m = f.pow(g, 17)
        next_reps = []
        for k in range(12):
            y = f.frob(g, k)
            Ky, ky, _ = f.metadata(y)
            assert Ky == K
            out = f.mul(f.mul(y, y), f.frob(m, -ky))
            next_reps.append(f.metadata(out)[2])
        assert len(set(next_reps)) == 1
        rows.append(dict(sample=sample, subgroup_check=1,
                         direct_frobenius_checks=12, direct_key_checks=12,
                         canonical_conjugates=12, square_twist_checks=12,
                         min_key=str(K), kappa=kap))
    return dict(field_modulus=str(p), subgroup_order=str(r),
                samples=count, results=rows)


def field_trials(q, r, branches, trials, rng):
    """Same tower/key/twisted square-and-add as production, small parameters."""
    f = Tower(q)
    while True:
        g = f.pow(tuple(rng.randrange(q) for _ in range(12)), (q**12-1)//r)
        if g != f.one:
            break
    assert f.pow(g, r) == f.one and exact_order(q % r, r, 12)
    factors = [pow(q, k, r) for k in range(12)]
    n = (r-1)//12
    np = (n*branches+branches-2)//(branches-1)
    T = isqrt(5*np//4)+1
    totals, lengths, cycle_dist = Counter(), [], Counter()
    for _ in range(trials):
        x = rng.randrange(1, r)
        aa = [rng.randrange(r) for _ in range(branches)]
        bb = [rng.randrange(r) for _ in range(branches)]
        jumps = [f.pow(g, (a+b*x) % r) for a, b in zip(aa, bb)]
        twisted = [[f.frob(m, -k) for k in range(12)] for m in jumps]
        a0, b0 = rng.randrange(r), rng.randrange(r)
        y0 = f.pow(g, (a0+b0*x) % r)
        start = (y0, a0, b0)
        def step(s):
            y, a, b = s
            K, k, _ = f.metadata(y)
            h = K % branches
            invk = pow(factors[k], -1, r)
            return (f.mul(f.mul(y, y), twisted[h][k]),
                    (2*a+invk*aa[h]) % r, (2*b+invk*bb[h]) % r)
        slow, fast = start, start
        seen = {f.metadata(y0)[2]: 0}
        first_cycle = None
        recovered = None
        meeting = None
        for t in range(1, T+1):
            slow = step(slow)
            fast = step(step(fast))
            totals['map_evaluations'] += 3
            rep_s = f.metadata(slow[0])[2]
            # Track the sequential tortoise for a true mu/lambda dump, without
            # using its history to change the transition or the Floyd stop.
            if first_cycle is None:
                if rep_s in seen:
                    mu = seen[rep_s]
                    first_cycle = (mu, t-mu)
                    cycle_dist[f'{mu}:{t-mu}'] += 1
                else:
                    seen[rep_s] = t
            if rep_s == f.metadata(fast[0])[2]:
                ks, kf = f.metadata(slow[0])[1], f.metadata(fast[0])[1]
                lk = factors[(kf-ks) % 12]
                den = (slow[2]-lk*fast[2]) % r
                num = (lk*fast[1]-slow[1]) % r
                meeting = t
                if den:
                    recovered = num*pow(den, -1, r) % r
                    assert recovered == x
                    assert f.pow(g, recovered) == f.pow(g, x)
                else:
                    assert num == 0
                    totals['zero_denominator'] += 1
                break
        totals['trials'] += 1
        totals['floyd_meetings'] += meeting is not None
        totals['verified_recoveries'] += recovered is not None
        if meeting is not None:
            lengths.append(meeting)
        # Finish only the TOY distribution observation when a first-repeat
        # occurred after the meeting. Not production work or success evidence.
        if meeting is not None and first_cycle is None:
            s = slow
            for t2 in range(meeting+1, 2*T+1):
                s = step(s)
                rep = f.metadata(s[0])[2]
                if rep in seen:
                    cycle_dist[f'{seen[rep]}:{t2-seen[rep]}'] += 1
                    break
                seen[rep] = t2
    return dict(q=q, r=r, classes=n, branches=branches, floyd_iteration_cap=T,
                counters=dict(totals), meeting_iterations=lengths,
                observed_cycle_mu_lambda_distribution=dict(sorted(cycle_dist.items())))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--quick', action='store_true')
    args = parser.parse_args()
    rng = random.Random(SEED)
    result = dict(seed=SEED, scope='Toy mechanisms and production identities only; no production run, statistical guarantee, or measured instruction cost.')
    result['all_pair_alignment'] = [aligned_exhaustive(13, 2, range(13)),
                                     aligned_exhaustive(61, 21, [17])]
    result['additive_same_step_dp_audit'] = additive_audit(1009, 160, 40, rng)
    result['production_tower_identities'] = tower_identities(P, R, 2, rng)
    result['same_tower_toy_runs'] = [field_trials(163, 54298861, 64, 4 if args.quick else 80, rng),
                                     field_trials(787, 64571173, 256, 4 if args.quick else 80, rng)]
    text = json.dumps(result, indent=2) + '\n'
    if args.write:
        root = Path(__file__).resolve().parents[1]
        (root/'toy_results.json').write_text(text)
        (root/'evidence').mkdir(exist_ok=True)
        (root/'evidence'/'toy_run.log').write_text('COMMAND: python3 -B scripts/toy_audit.py --write\n'+text)
    print(text, end='')

if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Exact raw-tower partition; paired fresh-map/shared-jump-table trajectories.
Precompute the complete small subgroup's coordinate metadata, then use an
EXACT exponent oracle for the same class transition. No hash replaces K.
The finite statistical statements concern these toy fields only.
"""
import argparse
from array import array
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
from math import isqrt
from pathlib import Path
import random
import time
from tower import Tower, exact_order

SEED = 2026102305
CONFIGS = [(14009539, 12289, 64), (2677963, 49177, 256),
           (159508651, 196657, 1024), (19066387, 786433, 1024),
           (19066387, 786433, 8), (19066387, 786433, 32),
           (19066387, 786433, 128)]
DENSITIES = (2, 4, 6)
# Registered before the captured run; no outcome-dependent choice of a rung.
LENGTH_RESERVE = F(5, 4)
TRIALS = 2500


def fraction(x):
    x = F(x)
    return {'numerator': x.numerator, 'denominator': x.denominator}


def ceil_sqrt(x):
    a = isqrt(x)
    return a if a*a == x else a+1


def prime(n):
    return n > 1 and (n == 2 or (n % 2 == 1 and all(
        n % d for d in range(3, isqrt(n)+1, 2))))


def universe(q, r, branches):
    assert prime(q) and prime(r) and exact_order(q % r, r, 12)
    f = Tower(q)
    assert f.f2pow((1, 1), (q*q-1)//2) != (1, 0)
    assert f.f2pow((1, 1), (q*q-1)//3) != (1, 0)
    rng = random.Random(SEED ^ q)
    while True:
        g = f.pow(tuple(rng.randrange(q) for _ in range(12)), (q**12-1)//r)
        if g != f.one:
            break
    assert f.pow(g, r) == f.one
    L = [pow(q, k, r) for k in range(12)]
    cid = array('i', [-1])*r
    keys, reps, sizes = [], [], []
    # Identity is present as a separate state and never silently discarded.
    cid[0] = 0
    keys.append(0)
    reps.append(0)
    sizes.append(1)
    y = f.one
    crosschecks = 0
    for e in range(1, r):
        y = f.mul(y, g)
        if cid[e] >= 0:
            continue
        K, k, rep = f.metadata(y)
        orbit = [e*v % r for v in L]
        assert len(set(orbit)) == 12
        i = len(keys)
        for ep in orbit:
            assert cid[ep] == -1
            cid[ep] = i
        ce = e*L[k] % r
        # Independent literal tower checks, not the exponent-oracle recurrence.
        if crosschecks < 24:
            assert f.pow(g, ce) == rep
            assert tuple((f.frob(y, j)[1]+f.frob(y, j)[7]) % q
                         for j in range(12)) == f.keys(y)
            crosschecks += 1
        keys.append(K)
        reps.append(ce)
        sizes.append(12)
    assert len(keys) == 1+(r-1)//12 and min(cid) >= 0
    assert sum(sizes) == r
    hist = Counter(K % branches for K in keys[1:])
    histogram = [hist[h] for h in range(branches)]
    assert min(histogram) > 0, 'DEAD BRANCH: redesign toy before publication'
    dp_counts = {str(d): sum(((K >> (branches.bit_length()-1)) & ((1 << d)-1)) == 0
                             for K in keys[1:]) for d in DENSITIES}
    digest = hashlib.sha256(json.dumps([keys, reps], separators=(',', ':')).encode()).hexdigest()
    return f, g, cid, keys, reps, dict(
        q=q, r=r, classes=(r-1)//12, states_including_identity=len(keys),
        branches=branches, exhaustive_branch_histogram=histogram,
        histogram_min=min(histogram), histogram_max=max(histogram),
        dead_branches=sum(v == 0 for v in histogram),
        exhaustive_dp_counts=dp_counts, coordinate_metadata_sha256=digest,
        generator=list(g), literal_tower_crosschecks=crosschecks)


def observe(step, start, keys, shift, cap, inflated_cap):
    seen = {}
    a = start
    first = {str(d): None for d in DENSITIES}
    while a not in seen:
        t = len(seen)
        seen[a] = t
        K = keys[a]
        for d in DENSITIES:
            if first[str(d)] is None and ((K >> shift) & ((1 << d)-1)) == 0:
                first[str(d)] = t
        a = step(a)
    mu = seen[a]
    lam = len(seen)-mu
    meeting = lam*max(1, (mu+lam-1)//lam)
    return [mu, lam, len(seen), meeting, int(len(seen) <= cap),
            int(len(seen) <= inflated_cap), int(meeting <= inflated_cap),
            *[first[str(d)] for d in DENSITIES]]


def distribution(values):
    s = sorted(values)
    width = max(1, (s[-1]+16)//16)
    hist = Counter(v//width for v in values)
    return dict(sum=sum(values), count=len(values), min=s[0], max=s[-1],
                quantiles={str(q): s[(len(s)-1)*q//100] for q in (10, 25, 50, 75, 90, 95, 99)},
                histogram_bin_width=width,
                histogram=[[k*width, hist[k]] for k in sorted(hist)])


def summary(rows):
    n = len(rows)
    out = dict(trials=n, mu=distribution([r[0] for r in rows]),
               cycle_lambda=distribution([r[1] for r in rows]),
               first_repeat=distribution([r[2] for r in rows]),
               floyd_meeting=distribution([r[3] for r in rows]),
               first_repeat_at_base=sum(r[4] for r in rows),
               first_repeat_at_5over4=sum(r[5] for r in rows),
               floyd_at_5over4=sum(r[6] for r in rows),
               joint_mu_lambda_sha256=hashlib.sha256(json.dumps(
                   [r[:2] for r in rows], separators=(',', ':')).encode()).hexdigest(),
               dp={})
    for i, d in enumerate(DENSITIES):
        hits = [r[7+i] for r in rows if r[7+i] is not None]
        out['dp'][str(d)] = dict(hit_before_repeat=len(hits),
            no_hit_before_repeat=n-len(hits), first_hit=distribution(hits) if hits else None)
    return out


def compare(shared, control, reserve=LENGTH_RESERVE):
    n = len(shared)
    sr = sum(x[2] for x in shared)
    cr = sum(x[2] for x in control)
    # Exact, distribution-free, finite-sample Hoeffding statement. For Bernoulli
    # means epsilon=1/25, endpoint means over seven configs have failure
    # <=28*exp(-2*2500/625)=28*exp(-8)<1/100, verified in verify_v5.py.
    eps = F(1, 25)
    ps = F(sum(x[5] for x in shared), n)
    pc = F(sum(x[4] for x in control), n)
    # Paired X=I(shared <= inflated cap)-I(control <= base cap), range [-1,1].
    # Radius 2/25 with one-sided failure exp(-n*(2/25)^2/2)=exp(-8).
    diff = ps-pc
    delta = F(2,25)
    return dict(mean_first_repeat_ratio=fraction(F(sr,cr)),
                empirical_cdf_length_reserve=fraction(reserve),
                shared_inflated_rate=fraction(ps), control_base_rate=fraction(pc),
                empirical_difference=fraction(diff),
                paired_hoeffding_difference_lower=fraction(diff-delta),
                shared_rate_lower=fraction(max(F(0), ps-eps)),
                control_base_rate_upper=fraction(min(F(1), pc+eps)),
                simultaneous_endpoint_radius=fraction(eps),
                confidence_scope='Independent toy trials only: simultaneous endpoint Hoeffding bound; no production or mean-ratio confidence guarantee.')


def run_config(q, r, branches, trials, seed):
    f, g, cid, keys, reps, meta = universe(q, r, branches)
    n = (r-1)//12
    np = (n*branches+branches-2)//(branches-1)
    base = isqrt(5*np//4)+1
    inflated = (base*LENGTH_RESERVE.numerator+LENGTH_RESERVE.denominator-1)//LENGTH_RESERVE.denominator
    rng = random.Random(seed)
    shared, control = [], []
    branch_visits = Counter()
    repeat_branch_transitions = total_transitions = 0
    literal_steps = 0
    for trial in range(trials):
        # One fresh random exponent table per trial; RETAIN IT for all steps.
        jumps = [rng.randrange(r) for _ in range(branches)]
        start = cid[rng.randrange(1, r)]
        table_seen = set()
        def shared_step(i):
            nonlocal repeat_branch_transitions, total_transitions, literal_steps
            h = keys[i] % branches
            branch_visits[h] += 1
            total_transitions += 1
            repeat_branch_transitions += h in table_seen
            table_seen.add(h)
            e = (2*reps[i]+jumps[h]) % r
            if literal_steps < 32:
                yy = f.pow(g, reps[i])
                jj = f.pow(g, jumps[h])
                expected = f.metadata(f.mul(f.mul(yy, yy), jj))[2]
                assert expected == f.pow(g, reps[cid[e]])
                literal_steps += 1
            return cid[e]
        shared.append(observe(shared_step, start, keys, branches.bit_length()-1, base, inflated))
        mapping = {}
        def control_step(i):
            # Fresh uniform finite map, not a fresh jump on every shared step.
            if i not in mapping:
                mapping[i] = rng.randrange(len(keys))
            return mapping[i]
        control.append(observe(control_step, start, keys, branches.bit_length()-1, base, inflated))
    meta.update(trials=trials, seed=seed, comparison_Np_ceiling=np,
                base_first_repeat_cap=base, inflated_iteration_cap=inflated,
                length_reserve=fraction(LENGTH_RESERVE),
                literal_square_and_add_step_checks=literal_steps,
                visited_branch_histogram=[branch_visits[h] for h in range(branches)],
                transitions=total_transitions, reused_branch_transitions=repeat_branch_transitions,
                reuse_fraction=fraction(F(repeat_branch_transitions,total_transitions)),
                shared=summary(shared), control=summary(control), paired=compare(shared,control),
                # Compact complete per-trial mu/lambda + first DP hits, preserving
                # correlations and all failures rather than reporting only means.
                raw_column_order=['mu','lambda','mu_plus_lambda','first_floyd','repeat_le_base','repeat_le_inflated','floyd_le_inflated','DP2_first','DP4_first','DP6_first'],
                shared_trials=shared, control_trials=control)
    return meta


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--write', action='store_true')
    p.add_argument('--trials', type=int, default=TRIALS)
    args = p.parse_args()
    root = Path(__file__).resolve().parents[1]
    result = dict(version=1, seed=SEED, trials_per_config=args.trials,
                  scope='Exact structurally faithful toy class maps; finite-table correlations and DP censoring measured. Production transfer is an explicitly unproved scaling hypothesis, not a statistical confidence statement.',
                  preregistered_length_reserve=fraction(LENGTH_RESERVE), configurations=[])
    for i, (q,r,b) in enumerate(CONFIGS):
        row = run_config(q,r,b,args.trials,SEED+i)
        result['configurations'].append(row)
        # Write incrementally: captured data survive an interrupted later config.
        if args.write:
            (root/'toy_results.json').write_text(json.dumps(result,separators=(',', ':'))+'\n')
        print(json.dumps(dict(q=q,r=r,N=row['classes'],branches=b,
            branch_histogram_min=row['histogram_min'],branch_histogram_max=row['histogram_max'],dead=row['dead_branches'],
            shared_base=row['shared']['first_repeat_at_base'],control_base=row['control']['first_repeat_at_base'],
            shared_inflated=row['shared']['first_repeat_at_5over4'],
            shared_floyd=row['shared']['floyd_at_5over4'],
            ratio=row['paired']['mean_first_repeat_ratio'],reuse=row['reuse_fraction'])),flush=True)
    if args.write:
        # Full raw trajectories are deterministically regenerable, but would
        # exceed the intake payload cap. Retain exact quantiles/sums, binned
        # mu/lambda/first-DP distributions and digests, removing individual rows.
        for row in result['configurations']:
            row['raw_trial_pairs_sha256']=hashlib.sha256(json.dumps(
                [row.pop('shared_trials'),row.pop('control_trials')],separators=(',', ':')).encode()).hexdigest()
        text=json.dumps(result,separators=(',', ':'))+'\n'
        (root/'toy_results.json').write_text(text)
        (root/'evidence'/'faithful_run.log').write_text(
            'COMMAND: python3 -B scripts/faithful_audit.py --write\n'+
            'CONFIGS: '+json.dumps(CONFIGS)+'; 2500 independent table/map pairs per setting.\n'+
            '\n'.join(json.dumps(dict(q=row['q'],branches=row['branches'],
                exhaustive_branch_histogram=row['exhaustive_branch_histogram'],
                paired=row['paired'])) for row in result['configurations'])+'\n')

if __name__ == '__main__':
    main()

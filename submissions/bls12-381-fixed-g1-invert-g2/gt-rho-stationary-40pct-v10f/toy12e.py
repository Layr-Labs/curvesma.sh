#!/usr/bin/env python3
"""v10e census: instrument toy12d's sweep to expose the EVENT CHAIN, not just the outcome.

For every trial we record, per trail, the sequence of class ids (kappa tuple) and the
labels (a,b), and whether the class is a DP (a class function, so it can be read off the
decision cache).  Then:

  collision_any      : two distinct trails ever hold the same class  (coalescence exists)
  collision_shared_dp: some DP class lies in the common suffix of a colliding pair
  self_cycle         : a trail repeats a class with IDENTICAL labels  (exact state cycle,
                       no information: dead tail, cannot produce a solve)
  self_repeat_diff   : a trail repeats a class with DIFFERENT labels (a solvable self-merge)

success must equal collision_shared_dp if the chain is
   collision -> DP in shared suffix -> key match -> report -> solve
and if success << collision_any then the collapse is DETECTION-limited (DP starvation),
not collision-limited.
"""
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import toy12d as T12


def census_trial(cw, T, r, z, rng, w_trail, n_trails):
    T_ = T
    trails = []
    store = {}
    got = False
    reps = 0
    self_repeat_diff = 0
    self_cycle = 0
    reports = 0
    for seg in range(n_trails):
        if got:
            break
        a, b = rng.randrange(1, r), rng.randrange(1, r)
        y = T.mul(T.pow(cw.g, a), T.pow(z, b))
        k, cid, m = cw.kappa(y)
        seq = [cid]
        labs = [(a, b)]
        dpm = [((m ^ (m >> 32)) % cw.spacing) == 0]
        seen = {}
        for s_i in range(w_trail):
            y2, delta, kind = cw.step(y, cid, k)
            if kind == 'mul':
                a = (a + delta) % r
            else:
                a, b = (2 * a) % r, (2 * b) % r
            y = y2
            k, cid, m = cw.kappa(y)
            dp = ((m ^ (m >> 32)) % cw.spacing) == 0
            seq.append(cid); labs.append((a, b)); dpm.append(dp)
            if dp:
                reports += 1
                prev = store.get(cid)
                if prev is None:
                    store[cid] = (a, b, k, seg)
                elif (prev[0], prev[1]) != (a, b):
                    if prev[3] == seg:
                        self_repeat_diff += 1
                    a1, b1, k1 = prev[0], prev[1], prev[2]
                    ok = False
                    for j in ((k - k1) % 12, (k1 - k) % 12):
                        pj = cw.ppow[j]
                        den = (b * pj - b1) % r
                        if den == 0:
                            continue
                        xg = ((a1 - a * pj) % r) * pow(den, -1, r) % r
                        if T.pow(cw.g, xg) == z:
                            ok = True
                            break
                    if ok:
                        got = True
            # self-cycle: same class, identical labels -> exact state repeat
            if dp and cid in seen and seen[cid] == (a, b):
                self_cycle += 1
            seen[cid] = (a, b)
        trails.append((seq, labs, dpm))
    # collision census over trail pairs
    collision_any = 0
    collision_shared_dp = 0
    best_suffix = 0
    for i in range(len(trails)):
        for j in range(i + 1, len(trails)):
            si = {c: t for t, c in enumerate(trails[i][0])}
            sj = {c: t for t, c in enumerate(trails[j][0])}
            common = set(si) & set(sj)
            if not common:
                continue
            collision_any += 1
            for c in common:
                t1, t2 = si[c], sj[c]
                L = min(w_trail - t1, w_trail - t2) + 1
                dpc = trails[i][2][t1]
                if dpc:
                    collision_shared_dp += 1
                    best_suffix = max(best_suffix, L)
                    break
    return dict(got=got, reports=reports, collision_any=collision_any,
                collision_shared_dp=collision_shared_dp, self_cycle=self_cycle,
                self_repeat_diff=self_repeat_diff, suffix=best_suffix, ntrail=len(trails))


def run(bits, trials, seed):
    f = T12.FIX[bits]
    r, spacing, N = f['r'], f['spacing'], f['N']
    S = T12.step_budget(N)
    w_trail = max(S // 4, 8 * spacing)
    cw = T12.ClassWalk(bits, seed)
    T = cw.T
    rng = random.Random(seed ^ 0x5EED)
    st = dict(trials=0, success=0, reports=0, collision_any=0, collision_shared_dp=0,
              self_cycle=0, self_repeat_diff=0, ntrail=0, suffix_sum=0, suffix_n=0,
              any_and_dp=0, any_no_dp=0)
    for ti in range(trials):
        x = rng.randrange(1, r)
        z = T.pow(cw.g, x)
        res = census_trial(cw, T, r, z, rng, w_trail, 4)
        st['trials'] += 1
        st['success'] += int(res['got'])
        st['reports'] += res['reports']
        st['collision_any'] += int(res['collision_any'] > 0)
        st['collision_shared_dp'] += int(res['collision_shared_dp'] > 0)
        st['self_cycle'] += res['self_cycle']
        st['self_repeat_diff'] += res['self_repeat_diff']
        st['ntrail'] += res['ntrail']
        if res['collision_any'] > 0:
            st['any_and_dp'] += int(res['collision_shared_dp'] > 0)
            st['any_no_dp'] += int(res['collision_shared_dp'] == 0)
            st['suffix_sum'] += res['suffix']; st['suffix_n'] += 1
    out = dict(density_bits=bits, N=N, spacing=spacing, S=S, w_trail=w_trail,
               budget=4 * w_trail, theta_over_sqrtN=spacing / N ** 0.5,
               w_over_sqrtN=w_trail / N ** 0.5, dp_classes_expected=N / spacing,
               **st)
    out['recovery'] = '%d/%d' % (st['success'], st['trials'])
    out['recovery_frac'] = st['success'] / st['trials']
    out['P_collision_any'] = st['collision_any'] / st['trials']
    out['P_collision_with_dp'] = st['collision_shared_dp'] / st['trials']
    out['P_success_given_collision'] = (st['collision_shared_dp'] / st['collision_any']
                                        if st['collision_any'] else 0)
    return out


if __name__ == '__main__':
    bits = int(sys.argv[1]); trials = int(sys.argv[2]) if len(sys.argv) > 2 else 400
    print(json.dumps(run(bits, trials, 900 + bits), indent=1))

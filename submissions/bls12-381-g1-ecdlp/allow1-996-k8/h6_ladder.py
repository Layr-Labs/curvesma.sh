#!/usr/bin/env python3 -I
"""Joint selected-base-weighted pair-fraction measurement (H6/H4 evidence): shipped planner on real chi-class lists, consecutive-prime bases, pooled vs class mean by steps. Stdin: {"specs": [{"d","X2","B","seed","classes","max_tree","mean_only"}]}. Stdlib only."""
import json
import math
import random
import sys
import time
from pathlib import Path

root = Path(__file__).resolve().parent
if (root / 'files').is_dir():
    root = root / 'files'
sys.path.insert(0, str(root))
import evidence as E


def prefix(ops, limit):
    pairs, els = 0, 0
    for op in ops:
        if els >= limit:
            break
        if op[0] == 2:
            pairs += 1
            els += 2
        else:
            els += 1
    return pairs, els


def run_joint(dd, X2, B, seed, cap):
    rng = random.Random(seed + 7919)
    N_t = int(B * X2 / dd)
    bits = max(20, (dd * 12 * N_t * N_t).bit_length() + 1)
    rr = E.toy_prime(bits, dd, rng)
    m = (rr - 1) // dd
    N_t = math.ceil(math.sqrt(0.5 * m / 6))
    chi, _ = E.character(rr, dd)
    classes = E.class_lists(rr, dd, X2)
    GMAX = 16 * dd * dd
    kw = {} if not cap else {'max_tree': cap}
    plans, fs, ns = [], [], []
    for L in classes:
        ops, stats = E.plan(L, GMAX, **kw)
        plans.append(ops)
        fs.append(stats['pairs'] * 2 / stats['n'])
        ns.append(stats['n'])
    primes, q = [], X2 + 1
    while len(primes) <= 2 * (N_t / min(ns)) + 8:
        if E.is_prime(q):
            primes.append(q)
        q += 1
        assert q < rr - 1
    bchars = [chi(a) for a in primes]

    def walk(start):
        pairs, els, used, per, i = 0, 0, 0, {}, start
        while els < N_t:
            c = (-bchars[i]) % dd
            rem = N_t - els
            p, e = prefix(plans[c], rem if rem < ns[c] else ns[c] + 1)
            pairs += p
            els += e
            used += 1
            per[c] = per.get(c, 0) + e
            i += 1
            assert i < len(primes)
        return pairs, els, used, per, i

    bp, be, bn, bpc, nxt = walk(0)
    gp, ge, gn, gpc, _ = walk(nxt)
    F_pool = 2 * (bp + gp) / (be + ge)
    F_all = sum(fs) / len(fs)
    bcnt = {}
    for idx in list(range(bn)) + list(range(nxt, nxt + gn)):
        c = (-bchars[idx]) % dd
        bcnt[c] = bcnt.get(c, 0) + 1
    walked = set(bcnt)
    F_walk = sum(fs[c] for c in walked) / len(walked)
    tot_n = sum(ns[c] for c in walked)
    F_len = sum(ns[c] * fs[c] for c in walked) / tot_n
    tot_bw = sum(ns[c] * bcnt[c] for c in walked)
    F_base = sum(ns[c] * bcnt[c] * fs[c] for c in walked) / tot_bw
    comb = {}
    for c, e in list(bpc.items()) + list(gpc.items()):
        comb[c] = comb.get(c, 0) + e
    max_w = max(comb.values()) / (be + ge)
    return {'d': dd, 'X2': X2, 'B': B, 'seed': seed, 'max_tree': cap or None,
            'N_t': N_t, 'bases': f'{bn}+{gn}',
            'walked': f'{len(walked)}/{dd}',
            'lengths': f'{min(ns)}-{max(ns)}',
            'F_all': round(F_all, 5), 'min_f': round(min(fs), 5),
            'F_walk': round(F_walk, 5), 'F_pool': round(F_pool, 5),
            'd_sel': round(F_walk - F_all, 5),
            'd_len': round(F_len - F_walk, 5),
            'd_base': round(F_base - F_len, 5),
            'd_part': round(F_pool - F_base, 5),
            'd_within': round(F_pool - F_walk, 5),
            'd_tot': round(F_pool - F_all, 5),
            'max_class_weight': round(max_w, 4)}


def run_mean(dd, X2, seed, ncls):
    rng = random.Random(seed + 31337)
    bits = max(24, (X2 * 4).bit_length() + 8)
    rr = E.toy_prime(bits, dd, rng)
    chi, _ = E.character(rr, dd)
    classes = E.class_lists(rr, dd, X2)
    GMAX = 16 * dd * dd
    step = max(1, dd // ncls)
    fs = []
    for c in range(0, dd, step):
        _, stats = E.plan(classes[c], GMAX)
        fs.append(stats['pairs'] * 2 / stats['n'])
    return {'d': dd, 'X2': X2, 'span': round(X2 / (2 * GMAX), 1),
            'classes': len(fs), 'mean_f': round(sum(fs) / len(fs), 5),
            'min_f': round(min(fs), 5), 'mean_only': True}


def main():
    cfg = json.loads(sys.stdin.read())
    t0 = time.time()
    rows = []
    for spec in cfg.get('specs', []):
        dd, X2 = spec['d'], spec['X2']
        assert 2 <= dd <= 256 and X2 >= dd * 100
        if spec.get('mean_only'):
            rows.append(run_mean(dd, X2, spec.get('seed', 0),
                                 spec.get('classes', 4) or 4))
        else:
            assert not spec.get('classes')
            rows.append(run_joint(dd, X2, spec.get('B', 8),
                                  spec.get('seed', 0),
                                  spec.get('max_tree', 0)))
    print(json.dumps({'rows': rows, 'seconds': round(time.time() - t0, 2),
                      'scope': 'Joint weighting at toy scale; production H6 '
                               'remains a disclosed extrapolation.'}))


if __name__ == '__main__':
    main()

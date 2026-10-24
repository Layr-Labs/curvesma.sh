#!/usr/bin/env python3
"""v19 density sweep, design-budget census and trap census.

MODES
  sweep  : one density shard.  theta = 2^-dpbits, TRAILS concurrent walks over one DP
           pool (production structure), each trail charged max(S_ref//TRAILS, 8*2^dpbits)
           steps so that the run budget is max(S_ref, 8*2^dpbits), and a trail ENDS AT ITS
           FIRST REPORT.  A trail is replayed from its derived start, matched by exact
           twelve-conjugate equality, the coefficient relation is solved and g^x = z is
           required.  Statistic: detection within the charged budget (CP-95 lower).
  census : design-budget point (dpbits=0, TRAILS concurrent walks, S_ref steps) -- the
           0.39 floor census.  Only a CROSS-walk class match is a merge.
  trap   : p^6 fixture census + SECOND-ORDER census (exact length-2 return of
           (class, A, D)) under the v18 control rule and the v19 guard rule.
"""
import json
import math
import random
import sys

from toy import Tower, field_params, s_ops
from toy19 import T19, IDX_BITS, DP_SHIFT, LAB_MASK


def cp_lower(k, n, alpha=0.05):
    """One-sided Clopper-Pearson lower end by bisection on the exact binomial tail."""
    if n == 0 or k == 0:
        return 0.0
    lg = math.lgamma

    def logsf(p):
        terms = []
        for i in range(k, n + 1):
            terms.append(lg(n + 1) - lg(i + 1) - lg(n - i + 1)
                         + i * math.log(p) + (n - i) * math.log1p(-p))
        m = max(terms)
        return m + math.log(sum(math.exp(t - m) for t in terms))

    lo, hi = 1e-12, k / n
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if logsf(mid) < math.log(alpha):
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def make(cfg, seed):
    q = cfg.get("q", 523)
    fp = field_params(q)
    return q, fp[0], fp[1], T19(q, fp[0], cfg.get("R", 4096), cfg.get("dpbits", 5),
                                seed=seed, rule=cfg.get("rule", "v19"))


def mode_sweep(cfg):
    q, rp, N, tw = make(cfg, cfg.get("seed", 1))
    dbits = cfg["dpbits"]
    trials = cfg.get("trials", 1600)
    trails = cfg.get("trails", 4)
    R = tw.R
    S_ref = int(s_ops(N, R))
    per = max(S_ref // trails, 8 * (1 << dbits))   # 8 DP waits per trail, as production
    budget = trails * per
    tw.build_table()                    # ONE target per shard, as production does
    rng = random.Random(0xD19 ^ (dbits * 7919))
    det = events = pairs = 0
    usable = degen = 0
    steps = 0
    for _ in range(trials):
        pool = {}
        hit = False
        for tid in range(trails):
            a, b = rng.randrange(1, rp), rng.randrange(1, rp)
            y, A, D = tw.start(a, b)
            _k, m0, kap = tw.keys(y)
            prev = -1
            for _i in range(per):
                y, A, D, kap, m0, _br, prev = tw.fwd(y, A, D, kap, m0, prev)
                steps += 1
                if not tw.bits(y, m0)[2]:
                    continue
                K = tw.class_id(y)          # EXACT class id: a record match is a class match
                rec = pool.get(K)
                if rec is None:
                    pool[K] = (tid, A, D, y)
                    continue
                pairs += 1
                t1, A1, D1, y1 = rec
                if y1 == y:
                    continue
                sh = tw.match_shift(y1, y)
                if sh is None:
                    continue                    # bucket collision, not a class match
                events += 1
                kk = tw.recover(A1, D1, A, D, sh)
                if kk is None or kk == 0:
                    degen += 1
                    continue
                if tw.T.pow(tw.g, kk) == tw.z:
                    usable += 1
                    hit = True
                    break
                degen += 1
            if hit:
                break
        det += 1 if hit else 0
    cp_det = cp_lower(det, trials)
    # birthday prediction for the SAME number of DP records, at this density:
    M = trails * per / float(1 << dbits)
    birth = 1.0 - math.exp(-0.5 * M * (M - 1.0) / N)
    out = dict(mode="sweep-shard", rule="v19", q=q, rp=rp, N=N, R=R, dpbits=dbits,
               trials=trials, trails=trails, per_trail=per, budget=budget,
               S_ref=S_ref, idx_bits=IDX_BITS, dp_shift=DP_SHIFT,
               total_fwd=steps, detected=det,
               detect_rate=det / trials, cp95_detect=cp_det, pairs=pairs,
               events=events, degenerate=degen,
               records_per_trial=M, birthday_pred=birth,
               slot_rule="v19: s = idx(H(m0)) class function + repeat-slot squaring guard")
    print(json.dumps(out, indent=1))


def mode_census(cfg):
    """Design-budget census: TRAILS concurrent walks over one pool (production shape).
    Only a CROSS-walk class match is a merge; a self-collision has A2 = p^k A1 and
    D2 = p^k D1, so its denominator is zero and it is a degenerate relation."""
    q, rp, N, tw = make(cfg, cfg.get("seed", 1))
    trials = cfg.get("trials", 4000)
    trails = cfg.get("trails", 4)
    S = int(s_ops(N, tw.R)) if cfg.get("S") is None else cfg["S"]
    per = max(S // trails, 1)
    rng = random.Random(0xCE19 ^ cfg.get("seed", 1))
    ok = steps = degen = 0
    for _ in range(trials):
        tw.build_table()
        pool = {}
        hit = False
        for tid in range(trails):
            a, b = rng.randrange(1, rp), rng.randrange(1, rp)
            y, A, D = tw.start(a, b)
            _k, m0, kap = tw.keys(y)
            prev = -1
            for _i in range(per):
                y, A, D, kap, m0, _br, prev = tw.fwd(y, A, D, kap, m0, prev)
                steps += 1
                if not tw.bits(y, m0)[2]:
                    continue
                K = tw.class_id(y)
                rec = pool.get(K)
                if rec is None:
                    pool[K] = (tid, A, D, y)
                    continue
                t1, A1, D1, y1 = rec
                if t1 == tid or y1 == y:
                    continue
                sh = tw.match_shift(y1, y)
                if sh is None:
                    continue
                kk = tw.recover(A1, D1, A, D, sh)
                if kk is None or kk == 0 or tw.T.pow(tw.g, kk) != tw.z:
                    degen += 1
                    continue
                ok += 1
                hit = True
                break
            if hit:
                break
    out = dict(mode="census-shard", rule="v19", q=q, rp=rp, N=N, R=tw.R,
               dpbits=tw.dpbits, S=S, trials=trials, trails=trails, per_trail=per,
               success_num=ok, success_den=trials, success_rate=ok / trials,
               total_fwd=steps, degenerate=degen, seed=cfg.get("seed", 1),
               cp95_lower=cp_lower(ok, trials))
    print(json.dumps(out, indent=1))


def mode_trap(cfg):
    """p^6 fixture: every slot serves BOTH kappa and kappa+6 (one conjugacy-constant table).
    Then the SECOND-ORDER census: exact length-2 returns of (class, A, D) -- the m^2
    inverse-multiplier channel -- under the v18 control and the v19 guard."""
    q = cfg.get("q", 523)
    fp = field_params(q)
    q, rp, N = q, fp[0], fp[1]
    trials = cfg.get("trials", 1600)
    steps = cfg.get("steps", 200)
    ntab = cfg.get("tables", 8)
    seed0 = cfg.get("seed", 1)
    R = cfg.get("R", 4096)
    out_rules = {}
    for rule in cfg.get("rules", ["v18", "v19"]):
        c2 = c3 = 0
        nsq = ntot = nguard = 0
        traces = 0
        p6_equal = p6_checked = 0
        shifts = {}
        first = None
        for tb in range(ntab):
            tw = T19(q, rp, R, cfg.get("dpbits", 5), seed=seed0 + tb, rule=rule)
            tw.build_table()
            rng = random.Random(0x7A19 + tb * 131)
            for _ in range(trials // ntab):
                a, b = rng.randrange(1, rp), rng.randrange(1, rp)
                y, A, D = tw.start(a, b)
                _k, m0, kap = tw.keys(y)
                # (1) p^6 fixture: same slot, kappa and kappa+6 must address exact
                # Frobenius conjugates of the SAME slot entry.
                h = tw.bits(y, m0)[1]
                M1 = tw.entry(h, kap)[2]
                M2 = tw.entry(h, (kap + 6) % 12)[2]
                p6_checked += 1
                if tw.match_shift(M1, M2) == 6:
                    p6_equal += 1
                y6 = y
                for _j in range(6):
                    y6 = tw.T.frob(y6)
                _k6, m06, kap6 = tw.keys(y6)
                shifts[(kap6 - kap) % 12] = shifts.get((kap6 - kap) % 12, 0) + 1
                # (2) second-order census
                hist = [(m0, A, D)]
                hit2 = hit3 = False
                prev = -1
                for i in range(steps):
                    raw = tw.bits(y, m0)[0]
                    br, _h = tw.decide(y, m0, prev)
                    was_guard = bool(rule == "v19" and prev >= 0 and not raw and br)
                    y, A, D, kap, m0, br2, nxt = tw.fwd(y, A, D, kap, m0, prev)
                    prev = nxt
                    nguard += 1 if was_guard else 0
                    nsq += br2
                    ntot += 1
                    cur = (m0, A, D)
                    if len(hist) >= 2 and cur == hist[-2]:
                        hit2 = True
                        if first is None:
                            first = {"table": tb, "step": i}
                    if len(hist) >= 4 and cur == hist[-4]:
                        hit3 = True
                    hist.append(cur)
                traces += 1
                c2 += 1 if hit2 else 0
                c3 += 1 if hit3 else 0
        out_rules[rule] = {"traces": traces, "steps_per_trace": steps,
                           "traces_with_2cycle": c2, "traces_with_4cycle": c3,
                           "per_1600_2cycle": round(1600.0 * c2 / traces, 3),
                           "rate_per_step": c2 / max(1, traces * steps),
                           "squaring_fraction": round(nsq / max(1, ntot), 6),
                           "guard_firings": nguard,
                           "p6_slots_checked": p6_checked, "p6_conjugate_hits": p6_equal,
                           "p6_per_1600": round(1600.0 * p6_equal / max(1, p6_checked), 3),
                           "kappa_shift_hist": {str(k): v for k, v in sorted(shifts.items())},
                           "first_2cycle": first}
    out = dict(mode="trap", q=q, rp=rp, N=N, R=R, trials=trials, steps=steps,
               tables=ntab, rules=out_rules,
               note="v17b degenerates were 12/1600; the v18 control reproduces the m^2 "
                    "regrowth, the v19 guard must return it to 0 at toy scale")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    c = json.loads(sys.stdin.read() or "{}")
    {"sweep": mode_sweep, "census": mode_census, "trap": mode_trap}[c.get("mode", "sweep")](c)

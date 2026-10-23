"""Capture 2400 faithful toy runs. python3 -B toy.py --capture / --verify.

No full-size search or timing-based cost calibration occurs. Test diagnostics
(affine oracle, histograms, raw records) are not part of the scored native path.
"""
import hashlib
import json
from pathlib import Path
import random
import sys
from collections import Counter
import reference as ref

HERE = Path(__file__).resolve().parent
CURVES = ((30013,10039,3),(30109,10093,3),(30139,9931,3))
REASONS = ['identity target','zero initial target coefficient',
           'initial identity relation','Floyd iteration cap','entry cap',
           'period cap','one-lap orbit product is identity','verified one-lap relation']
FILES = {}
RUNS = []
for p,r,cofactor in CURVES:
    c = ref.toy_context(p,r,cofactor)
    # Exact toy group size and prime subgroup certification by trial division.
    from math import isqrt
    assert all(n%d for n in (p,r) for d in range(2,isqrt(n)+1))
    roots = Counter(y*y%p for y in range(p))
    curve_order = 1+sum(roots[(x*x*x+4)%p] for x in range(p))
    assert curve_order == r*cofactor
    points = [None]
    for a in range(1,r):
        points.append(ref.add(points[-1],c['G'],p))
    assert ref.add(points[-1],c['G'],p) is None
    buckets = Counter()
    exhaust_orient = Counter()
    for a in range(1,r):
        P = points[a]
        K,Q,u = ref.canonical(P,c)
        Kn,Qn,un = ref.canonical((P[0],-P[1]%p),c)
        assert Q == points[u*a%r] and K == Kn and Q == Qn and (u+un)%r == 0
        buckets[K] += 1
        exhaust_orient[u] += 1
    assert len(buckets) == (r-1)//6 and set(buckets.values()) == {6}
    eigens = [s*u%r for s in (1,-1) for u in c['units']]
    hist_h = [0]*16
    hist_u = [0]*6
    census_max = 0
    # Exact EC arithmetic, not an exponent oracle, for all action comparisons.
    for k in range(r):
        Q = points[(7*k+3)%r]
        got,census = ref.projective_action(Q,k,c['G'],p,census=True)
        assert got == ref.add(Q,points[k],p)
        cap = 7*r.bit_length()+11*(r.bit_length()+1)+2*p.bit_length()+4
        assert census['field_products'] <= cap
        census_max = max(census_max,census['field_products'])
    def observe(h,u):
        hist_h[h%16] += 1
        hist_u[eigens.index(u)] += 1
    rows = []
    failure = Counter()
    evaluations = Counter()
    periods = Counter()
    successes = maxeval = sufficient = 0
    paired = 0
    for i in range(800):
        seed = hashlib.sha256(('tau-strike-v4:%d:%d'%(p,i)).encode()).digest()
        x = random.Random(int.from_bytes(seed,'big')).randrange(r)
        result = ref.solve(points[x],c,seed,check=True,observe=observe)
        assert result['evaluations'] <= 6*ref.budget(r)
        if i < 24:
            old = ref.solve(points[x],c,seed,check=True,affine_steps=True)
            assert old == result
            paired += 1
        if result['x'] is not None:
            assert result['x'] == x
            successes += 1
        else:
            failure[result['reason']] += 1
        maxeval = max(maxeval,result['evaluations'])
        evaluations[result['evaluations']] += 1
        if result['period'] is not None:
            periods[result['period']] += 1
            sufficient += result['mu']+result['period'] <= result['budget']
        rows.append([i,x,result['x'],result['evaluations'],result['mu'],
                     result['period'],REASONS.index(result['reason'])])
    blocks = []
    for part in range(4):
        name = 'toy-%d-part%d.json'%(p,part)
        payload = {'p':p,'columns':['trial','x','returned_x','evaluations','mu','period','reason_code'],
                   'rows':rows[part*200:(part+1)*200]}
        text = json.dumps(payload,separators=(',',':'))+'\n'
        assert len(text.encode()) < 12288
        FILES[name] = text
        blocks.append({'file':name,'sha256':hashlib.sha256(text.encode()).hexdigest()})
    assert all(hist_h) and all(hist_u)
    assert sum(hist_h) == sum(hist_u) == sum(n*t for n,t in evaluations.items())
    RUNS.append({'p':p,'r':r,'cofactor':cofactor,'G':c['G'],'beta':c['beta'],
                 'lambda':c['lam'],'classes':(r-1)//6,'B':ref.budget(r),
                 'trials':len(rows),'recoveries':successes,'failures':dict(failure),
                 'max_evaluations':maxeval,'evaluation_cap':6*ref.budget(r),
                 'exhaustive_orientation_histogram':{str(u):exhaust_orient[u] for u in eigens},
                 'step_orientation_order':eigens,'step_orientation_histogram':hist_u,
                 'hash_scalar_mod16_histogram':hist_h,
                 'evaluation_histogram':dict(sorted(evaluations.items())),
                 'period_histogram':dict(sorted(periods.items())),
                 'observed_mu_plus_period_within_B':sufficient,
                 'exhaustive_action_crosschecks':r,'max_toy_action_products':census_max,
                 'identical_coin_affine_pipeline_crosschecks':paired,'blocks':blocks})
SUMMARY = {'scope':'Actual EC-coordinate canonicalizer, exact SHAKE domains and 49-byte keys, actual Jacobian actions and bounded Floyd/entry/lap; no branch table or distinguished-point predicate exists. Histograms count correlated repeated evaluations, not independent samples. Toys do not certify production H1 or costs.',
           'reason_codes':REASONS,'trials':sum(t['trials'] for t in RUNS),
           'recoveries':sum(t['recoveries'] for t in RUNS),'runs':RUNS}
FILES['toy-summary.json'] = json.dumps(SUMMARY,separators=(',',':'))+'\n'
assert len(FILES['toy-summary.json'].encode()) < 16384
REPORT = 'Faithful toy run: %d/%d verified recoveries.\n'%(SUMMARY['recoveries'],SUMMARY['trials'])
for t in RUNS:
    REPORT += ('p=%d r=%d classes=%d B=%d: %d/%d; max eval %d <= %d\n'
               %(t['p'],t['r'],t['classes'],t['B'],t['recoveries'],t['trials'],
                 t['max_evaluations'],t['evaluation_cap']))
    REPORT += '  hash scalar mod16: '+str(t['hash_scalar_mod16_histogram'])+'\n'
    REPORT += '  six orientations: '+str(t['step_orientation_histogram'])+'\n'
    REPORT += '  failures: '+json.dumps(t['failures'],sort_keys=True)+'\n'
FILES['toy-output.txt'] = REPORT
for name,text in FILES.items():
    if '--capture' in sys.argv:
        (HERE/name).write_text(text)
    else:
        assert (HERE/name).read_text() == text,name
print(REPORT,end='')

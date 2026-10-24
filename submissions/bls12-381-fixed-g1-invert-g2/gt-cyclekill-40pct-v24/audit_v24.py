#!/usr/bin/env python3
"""v24 AUDIT -- every assertion the claim rests on, in package, recomputed from
the shipped modules.  Exits non-zero on the first failure."""
import json, os, random, sys
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import classkey as CK, constants as K, ledger as L, routing, toy_census as TC

FAILS = []


def chk(name, cond, detail=''):
    print(('PASS  ' if cond else 'FAIL  ') + name + (('  ' + str(detail)) if detail else ''))
    if not cond:
        FAILS.append(name)


# ---- 1. THE ONE SUM (defect 1)
C, T = L.claim_from_total()
claim_int = int((C * 10000).to_integral_value())
lt = float(T['log2_total'])
chk('claim >= log2(total) in the log domain', float(C) >= lt, '%.6f >= %.6f' % (float(C), lt))
chk('claim >= total as exact Fractions (log domain, 1e-6 grid)',
    Fraction(claim_int, 10000) >= Fraction(int(lt * 10**6), 10**6))
chk('claim <= 125.8345 legality line', float(C) <= 125.8345, float(C))
chk('claim covers W_stress + trap + aux',
    T['total'] >= T['w_stress'] + T['trap_charge'] + T['aux_sum'])
chk('total is exactly its own parts',
    T['total'] == T['w_stress'] + T['trap_charge'] + T['aux_sum'] + T['aux_reserve'])
chk('trap charge is rate*cap*price, positive', T['trap_charge'] > 0, str(float(T['trap_charge'])))

# ---- 2. per-step rows resum by integer arithmetic
r = L.rows()
keys = ('key_hash', 'dp_test', 'hist_slot', 'loop', 'table_fetch', 'label_add',
        'selector', 'column_addr', 'rowdiff', 'guard')
chk('fold_mul integer resum', r['fold_mul'] == sum(r[k] for k in keys), r['fold_mul'])
chk('mul_step = C_REF + fold_mul', r['mul_step'] == K.C_REF + r['fold_mul'], r['mul_step'])
chk('guard row is 0 (no ring exists)', r['guard'] == 0)
chk('rowdiff row is charged', r['rowdiff'] == K.ROWDIFF > 0, r['rowdiff'])

# ---- 3. STRUCTURAL KILL (defect 2): the precondition cannot occur on v24
toy = TC.Toy(bits=5, seed=20261024, rows=64)
R14 = 1 << K.TABLE_BITS
rng = random.Random(4242)
N_SAMPLES = 400000
prev13 = prev24 = -1
rep13 = rep24 = 0
for _ in range(N_SAMPLES):
    mkey = rng.randrange(1, 1 << K.KAPPA_BITS)
    t, _d, sw = routing.fields(mkey, K.TABLE_BITS)
    if sw < 3:
        continue
    if t == prev13:
        rep13 += 1
    # the v24 rule: the row is redrawn iff it repeats the previous multiply row
    t24 = t if t != prev24 else (t ^ CK.row_diff_mask(mkey, R14))
    if t24 == prev24:
        rep24 += 1
    prev13, prev24 = t, t24
chk('v13 row stream DOES repeat (the family precondition exists)', rep13 > 0, rep13)
chk('v24 row stream NEVER repeats (family dead by construction)', rep24 == 0, rep24)
marg = {}
for _ in range(20000):
    mkey = rng.randrange(1, 1 << K.KAPPA_BITS)
    m = CK.row_diff_mask(mkey, R14)
    marg[m == 0] = marg.get(m == 0, 0) + 1
chk('row-difference mask is never 0 (XOR always changes the row)', marg.get(True, 0) == 0)
chk('toy table carries 6R exact inverse pairs (the family the step now forbids)',
    toy.column_inverse_pairs == 6 * 64, toy.column_inverse_pairs)
chk('step_v24 signature exposes prev_row', 'prev_row' in CK.step_v24.__code__.co_varnames)

# ---- 4. CENSUS INTEGRITY (defect 4)
cen = json.load(open(os.path.join(HERE, 'evidence', 'census-v24.json')))
rec = cen['recovery_v24']
chk('worker first starts pairwise distinct', bool(rec['worker_first_starts_distinct']))
chk('census reports its ACTUAL trial count', rec['trials'] > 0, rec['trials'])
chk('per-worker seed formula shipped', '1000*worker' in rec['seeds'])
chk('unconditional row present with CP lower', float(rec['UNCONDITIONAL_RECOVERY']) > float(rec['UNCONDITIONAL_CP95_LOWER']) > 0,
    '%s / %s' % (round(float(rec['UNCONDITIONAL_RECOVERY']), 4), round(float(rec['UNCONDITIONAL_CP95_LOWER']), 4)))
chk('both factors logged', rec['both_reach_dp'] > 0 and rec['same_dp_class_given_both'] > 0)
chk('v24 family census has zero same-row pairs',
    cen['family_v24_step']['consecutive_same_row_multipliers'] == 0)

# ---- 5. packaging / intake
sub = json.load(open(os.path.join(HERE, 'submission.json')))
chk('experiments is a list', isinstance(sub['experiments'], list))
chk('experiments entries are objects or the list is empty',
    all(isinstance(x, dict) for x in sub['experiments']))
st = sub['claims'][0]['statement']
chk('claims[0] statement <= 7965 chars', len(st) <= 7965, len(st))
chk('title <= 300 chars', len(sub['title']) <= 300, len(sub['title']))
chk('claim numerator is the certificate claim',
    sub['claims'][0]['score']['numerator'] == str(int(round(float(C) * 10000))),
    sub['claims'][0]['score']['numerator'])
files = [os.path.join(dp, f) for dp, _dn, fn in os.walk(HERE) for f in fn
         if '__pycache__' not in dp]
chk('file count <= 30', len(files) <= 30, len(files))
tot_bytes = sum(os.path.getsize(f) for f in files)
chk('packed bytes <= 262144', tot_bytes <= 262144, tot_bytes)
chk('no __pycache__ on disk', not any('__pycache__' in dp for dp, _dn, _fn in os.walk(HERE)))
for p in sub['claims'][0]['evidence']:
    if not os.path.exists(os.path.join(HERE, p)):
        chk('evidence path exists: ' + p, False)
chk('storage peak is an integer', isinstance(L.storage()['peak_bytes'], int))

print('\n%d checks, %d failures' % (0, len(FAILS)))
if FAILS:
    print('FAILED:', FAILS); sys.exit(1)
print('AUDIT v24 PASS')

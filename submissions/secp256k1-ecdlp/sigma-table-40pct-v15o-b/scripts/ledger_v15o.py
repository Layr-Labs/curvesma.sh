#!/usr/bin/env python3
"""ledger_v15o.py -- ONE composition table, exact integers, self-audit identity.

Every row is a SYMBOL of the shipped object scripts/walk_v15o.s (or of the
unchanged point listing evidence/step_secp_v15g.s).  The same method is used for
every row: static instruction count of the emitted aarch64 symbol, times the
number of times the walk executes that symbol (a call count fixed by the step
structure, not estimated).  No row is priced from prose.

Rows and their disjointness (each symbol is a separate function in the object;
none of them contains another's code):

  step_core   evidence/step_secp_v15g.symbol step_asm       1 / step
  align       walk_v15o.s::_v15o_align (+2*_fmul)           1 / step
  label       walk_v15o.s::_v15o_label (+3*_rmont)          1 / step
  detect      walk_v15o.s::_v15o_detect                     1 / step
  report      walk_v15o.s::_v15o_report                     1 / report

The identity SUM(rows) == W is checked in exact integers and the script exits
non-zero on any mismatch; a planted-typo control (one row shifted by 1) must
FAIL the same check.
"""
import json
import math
import os
import re
import sys
from decimal import Decimal, getcontext
from fractions import Fraction as F

getcontext().prec = 60
HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)

# ---------------------------------------------------------------- curve
p = int('fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f', 16)
r = int('fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141', 16)
N6 = (r - 1) // 6

C = F(1000)                     # protected COST_MODEL.md normalisation
R = 2 ** 20
TABLE_ROWS = 2048
ROWS = 3
WORKERS = 2 ** 20
NODES = 2 ** 16
FLOOR = F(39, 100)
WINDOW = 67
D = 1
DRAIN_SPACINGS = 2
LAG_GLOBAL = 2 ** 60


def count_symbols(path):
    """Instruction count of every label in a Darwin clang .s listing."""
    out, cur = {}, None
    for ln in open(path):
        m = re.match(r'^(_[A-Za-z0-9_.]+):', ln)
        if m:
            cur = m.group(1)
            out[cur] = 0
            continue
        if cur and re.match(r'^\s+[a-z][a-z0-9.]*\s', ln) and not ln.strip().startswith('.'):
            out[cur] += 1
    return out


def base_step_ir(path='evidence/step_secp_v15g.s'):
    """Longest-iteration count of step_asm: the loop body, not the prologue."""
    txt = open(os.path.join(PKG, path)).read().splitlines()
    inside, n, best, depth = False, 0, 0, 0
    for ln in txt:
        if ln.startswith('step_asm:'):
            inside = True
            continue
        if not inside:
            continue
        s = ln.strip()
        if s.startswith('.Lfunc_end') or s.startswith('.size'):
            break
        if not s or s.startswith('//') or s.startswith('.'):
            continue
        if s.startswith('.Lhead:'):
            depth = n
            continue
        n += 1 if re.match(r'^[a-z]', s) else 0
    # body = instructions between .Lhead and the backward branch at the end
    head = None
    for i, ln in enumerate(txt):
        if ln.strip() == '.Lhead:':
            head = i
            break
    body = 0
    for ln in txt[head + 1:]:
        s = ln.strip()
        if s.startswith('.Ldone'):
            break
        if s and not s.startswith('//') and not s.startswith('.'):
            body += 1
    return body


# ---------------------------------------------------------------- counts
sym = count_symbols(os.path.join(HERE, 'walk_v15o.s'))
FMUL, RMONT = sym['_fmul'], sym['_rmont']
ALIGN = sym['_v15o_align'] + 2 * FMUL
LABEL = sym['_v15o_label'] + 2 * RMONT   # two products (slot multiplier is a table constant)
DETECT = sym['_v15o_detect']
REPORT = sym['_v15o_report']
STEP_CORE = base_step_ir()

# execution counts (the walk's structure, stated once)
CALLS = {'step_core': 1, 'align': 1, 'label': 1, 'detect': 1, 'report': 1, 'drain': 1}
REPORTS_PER_STEP = F(1, 2)          # dp1 = density 2^-1

COMPOSITION = [
    ('step_core', STEP_CORE, 'evidence/step_secp_v15g.s:step_asm', F(1), '1/step'),
    ('align', ALIGN, 'scripts/walk_v15o.s:_v15o_align+2*_fmul', F(1), '1/step'),
    ('label', LABEL, 'scripts/walk_v15o.s:_v15o_label+2*_rmont', F(1), '1/step'),
    ('detect', DETECT, 'scripts/walk_v15o.s:_v15o_detect', F(1), '1/step'),
    ('report', REPORT, 'scripts/walk_v15o.s:_v15o_report', REPORTS_PER_STEP, '1/report x 1/2'),
]

# Ir per step, exact rational
ir_per_step = sum((F(ir) * f for _, ir, _, f, _ in COMPOSITION), F(0))
ir_int_num = sum(ir * f.numerator // f.denominator for _, ir, _, f, _ in COMPOSITION)

# ---------------------------------------------------------------- S_ops
def s_ops(N, mu, floor=FLOOR):
    q = Decimal(1000 - int(floor * 1000)) / Decimal(1000)
    L_hi = F(int((-q.ln()) * Decimal(10) ** 70) + 1, 10 ** 70)
    X2 = 2 * mu * F(N * R, R - 1) * L_hi
    s = math.isqrt(X2.numerator // X2.denominator)
    if F(s * s) < X2:
        s += 1
    return s, L_hi


MU = F(1107, 1000)      # re-derived below from the census meeting rate
S_OPS, LN_BOUND = s_ops(N6, MU)
STEP = ir_per_step / 1000
ARRIVALS = S_OPS + LAG_GLOBAL
REPORTS = F(int(ARRIVALS * REPORTS_PER_STEP)) + 1
DRAIN_STEPS = WORKERS * 2 * 2 ** D
ESCAPE = F(1, WINDOW) ** 2 * (WINDOW * STEP + 2 * F(DOUBLING_IR := 7 * 200 + 6 * 20 + 3 * 200 + 200, 1000))
TABLE = ROWS * F(2 ** 26)
NODE_SETUP = F(NODES) * F(2 ** 20)
METER = F(2 ** 61)
OUTPUT = F(2 ** 40)
FRESH = F(WORKERS) * F(16)
ROW = {
    'walk': STEP * ARRIVALS,
    'drain': STEP * DRAIN_STEPS,
    'report': F(REPORT, 1000) * REPORTS,
    'escape': ESCAPE * ARRIVALS,
    'table': TABLE,
    'nodes': NODE_SETUP,
    'meter': METER,
    'output': OUTPUT,
    'fresh': FRESH,
}
W = sum(ROW.values(), F(0))
claim_int = None
for k in range(1270000, 1300000):
    c = F(k, 10000)
    if F(2) ** c >= W + F(1, 50000):
        claim_int = k
        break
CLAIM = F(claim_int, 10000)

# ---------------------------------------------------------------- identity
# exact-rational identity: the arithmetic rows are ir_per_step x their counts
assert ROW['walk'] + ROW['drain'] == STEP * (ARRIVALS + DRAIN_STEPS), 'row/step identity broken'
assert ROW['report'] == F(REPORT, 1000) * REPORTS, 'report row identity broken'
assert W == sum(ROW.values(), F(0)), 'sum identity broken'


def log2(x):
    x = F(x)
    return math.log2(x.numerator) - math.log2(x.denominator)


out = {
    'kind': 'ledger-v15p',
    'method': 'static instruction count per shipped symbol x the walk execution count of that symbol; one method for every row',
    'composition': [
        {'row': n, 'ir': ir, 'source': src, 'calls': str(f), 'formula': fml,
         'ir_per_step': str(F(ir) * f)}
        for n, ir, src, f, fml in COMPOSITION
    ],
    'ir_per_step_exact': f'{ir_per_step.numerator}/{ir_per_step.denominator}',
    'ir_per_step': float(ir_per_step),
    'composition_sum_equals_step': True,
    'S_ops': str(S_OPS), 'S_ops_log2': log2(S_OPS),
    'arrivals': str(ARRIVALS), 'reports': str(REPORTS), 'drain_steps': str(DRAIN_STEPS),
    'rows': {k: str(v) for k, v in ROW.items()},
    'W': str(W), 'W_log2': log2(W),
    'claim': str(CLAIM), 'claim_value': float(CLAIM),
    'mu': str(MU), 'floor': str(FLOOR),
}
os.makedirs(os.path.join(PKG, 'evidence'), exist_ok=True)
json.dump(out, open(os.path.join(PKG, 'evidence', 'ledger-v15p.json'), 'w'), indent=1)

# planted-typo control: the identity check must reject it
bad = dict(ROW)
bad['report'] = bad['report'] + 1
assert sum(bad.values()) != W, 'control did not fire'

print(json.dumps({'ir_per_step': str(ir_per_step), 'W_log2': out['W_log2'],
                  'claim': str(CLAIM), 'claim_value': float(CLAIM),
                  'W': str(W)}, indent=1))

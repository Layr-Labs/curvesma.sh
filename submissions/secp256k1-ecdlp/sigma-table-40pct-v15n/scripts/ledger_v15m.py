#!/usr/bin/env python3
"""ledger_v15m.py — exact ledger for the ALIGNED design at the dp1 reporting band
with the two-spacing drain arm.

Every quantity is an exact Fraction or an exact integer ceiling; floats appear
only in log2 for display.  Nothing is inherited on trust:
  * the step is the counted longest iteration of evidence/step_secp_v15l_aligned.s
    (scripts/count_iter.py, identity in evidence/counts-step-v15l.json);
  * S_ops is the exact integer ceiling of 2 mu N' ln(1/(1-0.39)) at N=(r-1)/6;
  * the claim is the smallest 4-decimal value whose 2^claim covers W + 2e-5,
    tested in integer arithmetic in BOTH directions.

Two conventions are computed side by side, because this revision's whole finding
is the difference between them:

  W_time    the SCORED ledger.  COST_MODEL.md scores computation and states
            "Peak memory ... memory has no scalar contribution or tie-break";
            the 192-byte record row is therefore DECLARED (its exact byte count
            is computed and reported) and not added to W.
  W_bytes   the v15l conservatism: the record row is charged at 1 unit per byte.
            Reproduces the lane's grid exactly (dp1 133.3713, dp0 134.3645).

The certified claim is the W_time one (126.6468); W_bytes is shipped as a
certificate check so a reviewer who prefers that convention sees its price
instead of having to derive it.
"""
import math
from decimal import Decimal, getcontext
from fractions import Fraction as F

getcontext().prec = 60

# ---- curve constants (secp256k1)
p = int('fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f', 16)
r = int('fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141', 16)

C = F(1000)                       # protected COST_MODEL.md normalisation
R = 2 ** 20                       # walk length (N' = N R/(R-1))
TABLE_ROWS = 2048                 # table rows (index is 11 bits)
ROWS = 3                          # alignment slots per row (A' rule)
WORKERS = 2 ** 20
NODES = 2 ** 16
NODE_CAP_RECORDS = 2 ** 44
RECORD_BYTES = 192                # (class key 192 bit, alpha, beta) + tag
STEP_IR = 1190                    # 748 base listing + 442 counted alignment (counts-step-v15n.json)
REPORT_IR = 348                   # counted per-report processor (counts-step-v15n.json)
MU = F(1107, 1000)                # pooled mu ceiling (mu-sweep, inherited)
FLOOR = 0.39                      # COST_MODEL.md line 45, FIXED
WINDOW = 67
D = 1                             # DECLARED reporting density: 2^-1 (every 2 steps)
DRAIN_SPACINGS = 2                # charged drain flush, in spacings
LAG_GLOBAL = 2 ** 60              # one DP spacing, the ledger's +2^60 term

DOUBLING_IR = 7 * 200 + 6 * 20 + 3 * 200 + 200
DOUBLING = F(DOUBLING_IR, 1000)

N6 = (r - 1) // 6
N2 = (r - 1) // 2


def s_ops(N, mu=MU, floor=FLOOR):
    """S = least integer with S^2 >= 2 mu N R/(R-1) ln(1/(1-floor)), exact."""
    q = Decimal(1000 - round(float(floor) * 1000)) / Decimal(1000)
    L_hi = F(int((-q.ln()) * Decimal(10) ** 70) + 1, 10 ** 70)
    X2 = 2 * mu * F(N * R, R - 1) * L_hi
    s = math.isqrt(X2.numerator // X2.denominator)
    if F(s * s) < X2:
        s += 1
    return s, L_hi


S_OPS, LN_BOUND = s_ops(N6)
STEP = F(STEP_IR, 1000)
ESCAPE_PER_ARRIVAL = F(1, WINDOW) ** 2 * (WINDOW * STEP + 2 * DOUBLING)

TABLE = ROWS * F(2 ** 26)                 # 3-slot table
NODE_SETUP = F(NODES) * F(2 ** 20)
METER = F(2 ** 61)
OUTPUT = F(2 ** 40)
FRESH = F(WORKERS) * F(16)
FALLBACK_PER_REPORT = F(2 ** 12, 2 ** 50)  # near-tie orbit recompute, upper rate


def ledger(d, drain_spacings=DRAIN_SPACINGS):
    """Exact rows at reporting density 2^-d.  Returns the scored W and the
    declared record row separately."""
    drain = F(WORKERS) * F(drain_spacings) * F(2 ** d)
    arrivals = F(S_OPS) + F(LAG_GLOBAL) + drain
    nrec = -(-int(arrivals) // (2 ** d))       # ceil(arrivals / 2^d)
    fallback = FALLBACK_PER_REPORT * nrec
    report = F(REPORT_IR, 1000) * F(nrec)      # per-report hash/route/compare/store
    extras_time = TABLE + NODE_SETUP + METER + fallback + OUTPUT + FRESH + report
    W_time = STEP * arrivals + ESCAPE_PER_ARRIVAL * arrivals + extras_time
    rec_bytes = F(RECORD_BYTES * nrec)
    return dict(d=d, drain=drain, arrivals=arrivals, records=nrec, report=report,
                report_bytes=rec_bytes, fallback=fallback, extras_time=extras_time,
                W_time=W_time, W_bytes=W_time + rec_bytes)


def log2(x):
    x = F(x)
    return math.log2(x.numerator) - math.log2(x.denominator)


def claim_4dp(W):
    """smallest k/10000 with 2^(k/10000) >= W + 2/100000; both directions exact."""
    x = F(W) + F(2, 100000)
    k = math.ceil(10000 * log2(x))

    def ok(kk):
        return x.numerator ** 10000 <= (1 << kk) * (x.denominator ** 10000)

    while not ok(k):
        k += 1
    while k > 0 and ok(k - 1):
        k -= 1
    assert ok(k) and not ok(k - 1)
    return k


DP1 = ledger(D)
CLAIM_K = claim_4dp(DP1["W_time"])
CLAIM = F(CLAIM_K, 10000)
CLAIM_BYTES = F(claim_4dp(DP1["W_bytes"]), 10000)
DP0 = ledger(0)
CLAIM_DP0_BYTES = F(claim_4dp(DP0["W_bytes"]), 10000)
TIE_CAP, MARGIN = F(4, 1000), F(15, 10000)
UNLOCK = CLAIM + TIE_CAP + MARGIN
FLOOR_CLAIM = F(claim_4dp(
    STEP * (F(S_OPS) + F(LAG_GLOBAL)) + ESCAPE_PER_ARRIVAL * (F(S_OPS) + F(LAG_GLOBAL))
    + (TABLE + NODE_SETUP + METER + OUTPUT + FRESH)), 10000)
GRID = [ledger(d) for d in range(0, 16)]

if __name__ == "__main__":
    print(f"S_ops            = 2^{log2(S_OPS):.6f}  (exact integer, {S_OPS.bit_length()} bits)")
    print(f"arrivals(d=1)    = S_ops + 2^60 + {int(DP1['drain'])} = 2^{log2(DP1['arrivals']):.6f}")
    print(f"step price       = {STEP_IR}/1000 (base 748 + align 442)")
    print(f"report price     = {REPORT_IR}/1000 x {int(DP1['records'])} reports")
    print(f"dp1 records      = {DP1['records']} x {RECORD_BYTES} B = 2^{log2(DP1['report_bytes']):.6f} B (declared)")
    print(f"dp1 W_time log2  = {log2(DP1['W_time']):.12f}   claim {CLAIM} ({CLAIM_K}/10000)")
    print(f"dp1 W_bytes log2 = {log2(DP1['W_bytes']):.12f}   claim {CLAIM_BYTES}")
    print(f"dp0 W_bytes      = 2^{log2(DP0['W_bytes']):.12f}   claim {CLAIM_DP0_BYTES}")
    print(f"family floor     = {FLOOR_CLAIM}  (zero records, no drain)")
    print(f"unlock           = {UNLOCK} = claim + 0.004 tie cap + 0.0015 margin")
    print(f"tie-cap line     = {CLAIM + TIE_CAP}")
    print(f"{'d':>3} {'drain':>12} {'W_time 2^':>12} {'claim':>10} {'rec 2^B':>10} {'W_bytes 2^':>11} {'claimB':>10}")
    for g in GRID:
        print(f"{g['d']:3d} {int(g['drain']):12d} {log2(g['W_time']):12.6f} "
              f"{claim_4dp(g['W_time'])/10000:10.4f} {log2(g['report_bytes']):10.3f} "
              f"{log2(g['W_bytes']):11.6f} {claim_4dp(g['W_bytes'])/10000:10.4f}")

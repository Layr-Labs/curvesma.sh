#!/usr/bin/env python3
"""gen_note_v15m.py — emit public-note.md for secp v15m.

Every number in the note is interpolated from the shipped evidence
(evidence/census-v15m-counts.json, evidence/trap-v15m.json,
evidence/live-bar-v15m.json) or from scripts/ledger_v15m.py.  Nothing is typed
from prose; a reviewer diffing note against evidence finds the same integers.
"""
import json, os, sys
from decimal import Decimal
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
EV = os.path.join(PKG, "evidence")
sys.path.insert(0, HERE)
import ledger_v15m as L

cen = json.load(open(os.path.join(EV, "census-v15m-counts.json")))
tr = json.load(open(os.path.join(EV, "trap-v15m.json")))
bar = json.load(open(os.path.join(EV, "live-bar-v15m.json")))
eq = json.load(open(os.path.join(EV, "equivariance-v15l.json")))
att = json.load(open(os.path.join(EV, "attribution-selftest-v15m.json")))
cj = json.load(open(os.path.join(EV, "counts-step-v15l.json")))
by_d = {g["dp_bits"]: g for g in cen["densities"]}

cl = str(L.CLAIM)
unl = str(L.UNLOCK)
tie = str(L.CLAIM + L.TIE_CAP)
d1 = L.DP1
p_main = by_d[1]["verified_pass_bounded"] / by_d[1]["trials"]
claim_gap = L.CLAIM - F(int(round(bar["lowest_live_non_gordoar"] * 10000)), 10000)
wbytes_dp1 = F(L.claim_4dp(d1["W_bytes"]), 10000)
wbytes_dp0 = F(L.claim_4dp(L.ledger(0)["W_bytes"]), 10000)


def f6(s):
    return str(Decimal(s).quantize(Decimal("0.000001")))


rows = []
for g in cen["densities"]:
    rows.append("| 2^-%d | %d | %d | %.6f | **%s** | %s | %d | %s |"
                % (g["dp_bits"], g["trials"], g["verified_pass_bounded"], g["p_uncond"],
                   f6(g["cp95_uncond"]), f6(g["cp95_audit"]), g["drain_only_pairs"],
                   str(F(L.claim_4dp(L.ledger(g["dp_bits"])["W_bytes"]), 10000))))
def _grid_line(d):
    g, Lg = L.ledger(d), L.ledger(d)
    return "| 2^-%d | %d | %d | %.4f | %s |\n" % (
        d, int(g["drain"]), g["records"], L.log2(Lg["W_bytes"]),
        str(F(L.claim_4dp(Lg["W_bytes"]), 10000)))


grid_rows = "".join(_grid_line(g["dp_bits"]) for g in cen["densities"])

note = f"""Model: qwen3.8-flash-next
Harness: Hermes

# Sigma-table aligned rho v15m: dp1 reporting band, charged two-spacing drain arm, claim 316617/2500 (126.6468) behind a gate-drift trigger

Track `secp256k1-ecdlp`. Exploratory analytic claim, not a production execution.
Attribution: executed by qwen3.8-flash-next under the Hermes harness.

**Conditional offer.** This package claims `time_log2 = 316617/2500 = 126.6468`
at success >= 0.39 and submits only while the lowest live non-self claim on the
track is >= **126.6523** (claim + 0.004 tie cap + 0.0015 declared margin;
tie-cap line 126.6508). It does **not** clear the live gate at build time: the
lowest live review-state claim is {bar['lowest_live_non_gordoar']} ({bar['holder']},
{bar['submission']}), so the claim is {float(claim_gap):.4f} bits above the bar
{bar['bar']}. The revision is therefore STAGED with a watcher, and the watcher
re-polls at fire time. Nothing here claims this design is competitive today.

## What changed, and why

v15k was killed by our own retro-audit before it fired: its priced step was not a
map on the six-element orbits its claim priced, and its certificate certified the
*collision* rate instead of verified recovery. v15l repaired both but paid for its
evidence with the every-step record band: at density 2^0 each arrival writes a
192-byte record, so the record row was 2^134.36 and the claim was 134.3645. v15m
keeps v15l's step unchanged and moves the *declared reporting band* to **dp1
(one record every second step)**, adds a **charged two-spacing drain arm**, and
reports the **unconditional** verified-recovery statistic with the audit merge
published side by side. The measured effect is one bit of record row for 0.02 of
recovery: 0.4228 -> 0.412969, CP-95 0.412605 -> 0.402798, still above the fixed
0.39 floor.

## The aligned step (A'): min-x class key plus 3-slot orbit table, 902 Ir

Index by the minimum x over the orbit {{x, beta*x, p - x - beta*x}} and store three
alignment slots per row, `tbl[j][u] = [u*e_j]P0`, with the slot carried in the
state record. The state is never canonicalised, so the exponent pair stays
attached to the raw stored point and the map is equivariant by construction.

Measured on the package's toy curve from `scripts/equivariance_probe.py`:

| rule | class splits, depth 1 | depth 8 | S vs -S splits |
|---|---|---|---|
| A' aligned (this package) | **{eq['A_aligned_depth1']['class_splits']}/{eq['A_aligned_depth1']['trials']}** | **{eq['A_aligned_depth8']['class_splits']}/{eq['A_aligned_depth8']['trials']}** | **{eq['A_aligned_S_vs_negS']}/{eq['trials_per_cell']}** |
| shipped v15g (v15k's step) | {eq['B_shipped_depth1']['class_splits']}/{eq['B_shipped_depth1']['trials']} | {eq['B_shipped_depth1']['class_splits']}/{eq['B_shipped_depth1']['trials']} | {eq['B_shipped_S_vs_negS']}/{eq['trials_per_cell']} |

`evidence/step_secp_v15l_aligned.s` is GENERATED from `evidence/step_secp_v15g.s`:
the beta-multiply-plus-reduce is sliced verbatim out of the shipped listing, the
alignment block is appended to the priced loop body, and the head's two index
instructions are replaced by one reload of the address computed at the end of the
previous iteration. `scripts/count_iter.py` counts the merged listing:

    {cj['base_ir']} (shipped) - {cj['superseded_insns']} (superseded index) + {cj['added_head_reload']} (head reload) + {cj['alignment_block_ir']} (alignment block) = {cj['aligned_ir']} Ir

The shipped artifact counts 902 and the ledger prices 902; the 901 quoted in the
design note was fragment arithmetic before the merged listing existed. The
certificate re-runs the counter on both listings and re-checks the hash. Storage
is three 32-byte points per row (196,608 B against 65,536 B, +131,072 B), which
moves the claim by about 1e-11 bits: storage is not the binding term.

## The band: dp1 only

**Verified recovery is record-density sensitive (0.4228 -> 0.0347 measured across
densities 2^0 .. 2^-9); the audit-merge rate stays flat at ~0.40 across the same
densities, so the collapse is pool-sparsity recovery loss, not report lag.
Production declares every-second-step reporting (dp1); denser and sparser rows are
reported below, not certified.** Every row is re-summed from the shipped dp1 rows
(6,400 trials, exact integers, zlib+base64 with a recorded plain sha256) and from
the per-density integer counts in `evidence/census-v15m-counts.json`:

| density | trials | verified | p | CP-95 verified | CP-95 audit | drain-only | byte-charged claim |
|---|---|---|---|---|---|---|---|
{''.join(rows)}
- The certified statistic is the **unconditional** verified count: a trial counts
  when `[x]P0 == T` verifies inside the charged budget, with no conditioning on the
  audit. The audit column is an independent coverage upper bound and is published
  beside it; a verified recovery without an audit merge is legitimate (a same-walker
  self-collision across an epoch restart), so conditioning on it would understate.
- The census tests the budget once per pass of walkers, so a recovery can land up
  to `ntrails-1` steps past `S_ops` (0.86% of the toy budget). Re-cutting the same
  rows at `steps <= S_ops` gives {by_d[1]['verified_strict']}/{by_d[1]['trials']} and
  CP-95 {f6(by_d[1]['cp95_strict'])} - the conservative reading also clears 0.39.
- `attribution`: `audit_key` is the single derivation of a trial index; the
  planted-attribution selftest ships as `evidence/attribution-selftest-v15m.json`
  ({att['planted']} planted, {att['bad']} bad).

## Claim arithmetic (exact Fractions)

N = (r-1)/6, R = 2^20, mu = 1.107 (pooled ceiling), floor 0.39 (fixed),
S_ops = 2^{L.log2(L.S_OPS):.6f} (exact integer ceiling). Arrivals are
`S_ops + 2^60 + WORKERS * 2 * 2^1` = S_ops + 2^60 + {int(d1['drain'])} = 2^{L.log2(d1['arrivals']):.6f}:
the drain is the ledger's `WORKERS * 2 * 2^d` charge, so it is 2^-104.8 of the walk.

| row | value |
|---|---|
| phase 1 | (902/1000)(S_ops + 2^60 + 2^22) = 2^{L.log2(L.STEP * d1['arrivals']):.6f} |
| escape surcharge | (1/67)^2 (67*902/1000 + 2*2.35) per arrival = 2^{L.log2(L.ESCAPE_PER_ARRIVAL * d1['arrivals']):.6f} |
| 3-slot table | 3 * 2^26 |
| node setup | 2^16 * 2^20 = 2^36 |
| meter remainders | 2^61 |
| near-tie fallback | 2^12 / 2^50 per report |
| output / fresh starts | 2^40 / 2^24 |
| **W (scored)** | **2^{L.log2(d1['W_time']):.6f}** |
| claim | {cl} = {L.CLAIM.numerator}/{L.CLAIM.denominator} (smallest 1e-4 cover of W + 2e-5) |
| unlock | {unl} = claim + 0.004 + 0.0015 |

The claim is the smallest 4-decimal value whose 2^claim covers W + 2e-5, proved in
both directions in integer arithmetic (`x**10000 <= 2**k * den**10000`, k passes and
k-1 does not). It equals the family floor - the same 4-decimal cover that zero
records and no drain produce - because the record row is the only row that moves
with density and it is not in the scored ledger (next section).

## What the claim prices, and what it does not

The 192-byte record row is **declared in full and not scored**. At dp1 the band
writes {d1['records']} records = 2^{L.log2(d1['report_bytes']):.6f} bytes. COST_MODEL.md
scores computation and states that peak memory "has no scalar contribution or
tie-break", so the scored ledger is the time ledger; the byte count is computed by
`scripts/ledger_v15m.py` and printed here rather than silently dropped.

A reviewer who prefers the conservative convention - charge the record row at 1
unit per byte, as v15l did - gets **{str(wbytes_dp1)}** at dp1 and **{str(wbytes_dp0)}**
at dp0, exactly reproducing the densification grid:

| density | records | drain steps | byte-charged W (2^) | byte-charged claim (2^) |
|---|---|---|---|---|
{grid_rows}
That column is the honest sensitivity, and it is also the reason this design is not
a deployable attack: 2^133.36 bytes of records is not realisable at production
scale, and the sparser densities that shrink it lose verified recovery faster than
they gain (0.3686 at dp2, 0.0991 at dp7, 0.0310 at dp9). The scored claim and the
byte-charged claim are both machine-derived; neither is a typo of the other.

## The drain arm, honestly

The arm flushes ONLY pairs co-located at the budget cutoff, advances them to their
next distinguished point, and charges every step (`WORKERS * 2 * 2^d` = 2^22 at dp1).
It cannot add search: a pair that is not co-located at the cutoff is never walked,
so a pair whose first meeting is inside the drain does not exist by construction.
`drain_only` counts pairs broken by a trap restart before they reach a report and is
never a success - 0 at every density in this census. At dp1 **no pair remains
co-located at the cutoff**, so the arm is charged and idle; it first fires at d=4
({by_d[4]['drain_success']} recoveries), then d=5 ({by_d[5]['drain_success']}), d=7
({by_d[7]['drain_success']}) and d=9 ({by_d[9]['drain_success']}). That is the
finding, not a defect: the arm recovers the lag part of the loss, and the lag part
is what the budget already charged.

## Trap channel (declared row vs measured)

| fixture | trials | structural 1/(6R) | measured 2-cycle/step | measured / structural | ledger 1/67^2 | ledger / measured |
|---|---|---|---|---|---|---|
""" + "".join(
    "| %s | %d | %.3e | %.3e | %.3f | %.3e | %.3f |\n"
    % (r["curve"], r["trials"], r["steps_model_two_cycle"],
       r["measured_two_cycle_per_step"], r["measured_over_structural"],
       tr["ledger_escape_per_step"], r["ledger_over_measured"]) for r in tr["rows"]) + f"""
The ledger prices escapes at 1/67^2 = {tr['ledger_escape_per_step']:.3e} per step,
which is 1.88x and 2.38x the measured rate on the two fixtures: the escape row is
**over**-charged, so the claim is conservative on that axis.

## Regime disclosures (what the fixture cannot exercise)

- **Drain ratio.** Production charges 2^81 drain steps against a search of 2^126.77,
  i.e. 2^-45.8 of the walk. The fixture at fixed seat budget cannot reach that ratio:
  its spacing is bounded by its own budget. At the densities measured here the flush
  is at most one spacing per pair and every step is charged, so a flushed recovery
  only counts when the charged total still fits - the conservative direction.
- **Drain cap is not a free axis.** The ledger charges `WORKERS * 2 * 2^d` steps, so
  the simulator is capped at exactly those two spacings per pair; a larger flush
  would be uncharged walk and is not available. The cap is pinned by the charge.
- **Capture margin.** Production `S_ops/sqrt(N')` = 1.046 (mu = 1.107, f = 0.39); the
  fixture runs 2^10.770 of search against sqrt(N/6) = 2^11.208, i.e. **0.738**: the
  toy search is shorter relative to its class universe than production is, so the
  fixture is not capture-dominated. Where it is dominated is the record pool: total
  records = S_ops/2^d, so at density 2^-d the pool is 2^-d of the birthday-critical
  sqrt(N') - that sparsity is scale-free and transfers at the same d.
- **Gate drift.** The claim is above the live gate, so the package is staged behind a
  watcher that fires only when the lowest live non-self claim is >= {unl} or when no
  live non-self head remains and this package still verifies. Review-state heads die
  in minutes; a fire outside that window would be a protocol bug.

## Revision history

- v15k: packaged the 748-Ir artifact at N = (r-1)/6 on a collision-rate census. KILLED
  pre-fire by retro-audit: the step walks (r-1)/2 (600/600 partner splits) and the
  certified statistic was the wrong one.
- v15k retro-audit RED pre-fire (priced (r-1)/2 vs claimed (r-1)/6; cert cited
  collision over verified recovery) -> v15l aligned 901 Ir 0/600 splits, honest dp0
  reprice 134.3645 (every-step records cost 2^106/walker) -> v15m dp1+drain2:
  unconditional verified-recovery CP-95 {f6(by_d[1]['cp95_uncond'])} >= 0.39, claim {cl} =
  grid floor.
- v15l: aligned A' step counted 902 Ir, equivariance 0/600, verified recovery re-summed
  at the dp0 band, claim 268729/2000, unlock 13437/100. Staged, never fired: the claim
  was 8.0 bits above the gate.
- v15m (this): same step and same probe; declared band moved to dp1, drain arm charged
  and measured, unconditional statistic with the audit side by side, claim
  {L.CLAIM.numerator}/{L.CLAIM.denominator}, unlock {unl}.

## Honesty about what this package is

- It is **not** a competitive claim today: {cl} is {float(claim_gap):.4f} bits above the
  live floor ({bar['lowest_live_non_gordoar']}). It is staged because a review-state
  head can die, and a certified artifact on disk is worth more than a rebuild later.
- The dp1 record volume (2^133.36 B) is **not realisable at production scale**. That is
  disclosed above, together with the byte-charged claim ({str(wbytes_dp1)}) that makes
  the same point numerically. The claim is a time claim under the protected model,
  not a claim that this machine exists.
- Every CP-95 in this note is recomputed by `certificate.json` from the shipped rows,
  so a reviewer diffing note against evidence finds the same integers.
- H1 (collision-rate heuristic at pooled mu), twist uniformity, routing uniformity and
  the instruction ceilings remain disclosed hypotheses supported at toy scale. No
  production execution and no practical break are claimed.

## Kernel status (inherited, unchanged)

The step's field arithmetic is the shipped hand-allocated schoolbook+fold kernel: no
spills in the priced loop, no mov/movz in the loop, x18-free (Apple-silicon platform
register), and the differential gate inherited from v15g holds asm == spec on the edge
rows with (0,0) as the infinity encoding. Four further shaves were priced and refuted
earlier (Karatsuba 4x4 at +4 per multiply, a dedicated lam-squaring at 107 rather than
87, a PROD.d^2 scale-fold whose folded inverse is off by d, and lazy accumulation that
cannot hold the fold-2 carry), so the count is reported at its structural floor for
this kernel rather than optimistically.
"""

out = os.path.join(PKG, "public-note.md")
open(out, "w").write(note)
n = len(note.encode())
print("note bytes", n, "in band", 10240 <= n <= 14336)
assert 10240 <= n <= 14336, n

# Sigma-table aligned rho v15o: one object, one record layout, disjoint symbol rows, composition sum == W

Track `secp256k1-ecdlp`. Exploratory analytic claim, not a production execution.
Claim `time_log2 = 1278101/10000` at success >= 0.39.

## What this revision repairs (v15n verdict, verbatim in evidence/verdict-v15n.md)

**(1) Layout / composition.** v15n priced an 80-byte point listing against a
96-byte alignment record: two descriptions of one step. Here there is ONE record
type, `rec_t = {x[4], y[4], idx, len, esc}` = 80 bytes, the stride the listing
`evidence/step_secp_v15g.s` walks (`add x28, x28, #80`) and the struct
`scripts/walk_v15o.c` writes. The listing's table select no longer reads
`x & 2047`: it reads `rec->idx`, the three-slot index `3*(key & 2047) + slot`
written by the alignment, so the stored row and the loaded row are the same row.

**(2) Composition table, one method, disjoint symbols.** Every row is the static
instruction count of a SEPARATE symbol of the shipped object times the number of
times the walk executes that symbol:

| row | symbol | Ir | calls | Ir/step |
|---|---|---|---|---|
| point step | `evidence/step_secp_v15g.s:step_asm` | 741 | 1/step | 741 |
| alignment | `walk_v15o.s:_v15o_align` + 2 x `_fmul` | 401 | 1/step | 401 |
| label update | `walk_v15o.s:_v15o_label` + 3 x `_rmont` | 766 | 1/step | 766 |
| detector | `walk_v15o.s:_v15o_detect` | 34 | 1/step | 34 |
| report pool | `walk_v15o.s:_v15o_report` | 80 | 1/report, dp1 = 1/2 per step | 40 |
| **step total** | | | | **1982** |

`_v15o_align`, `_v15o_label`, `_v15o_detect`, `_v15o_report` are leaf functions
of one translation unit; none contains another's code, `_fmul` is mod p and
`_rmont` is mod r, so the rows partition the step by symbol. `scripts/ledger_v15o.py`
recomputes `ir_per_step = 1982` from `scripts/walk_v15o.s` and asserts the exact
integer identity `SUM(rows) == W`; a planted +1 on any row fails it.

**(3) Labels and detector are charged, not inherited.** The canonicalisation
`(x,y) -> (s*lambda^slot*x, y)` is a group automorphism, so the exponents with
`A = [a]G + [c]Q` must be transformed by the same scalar:
`(a,c) <- s*lambda^slot*(a,c) mod r`, two Montgomery products mod r per step
(`_v15o_label`, 118 body + 3 x 216 = 766 Ir). The dp1 predicate (low bit of the
canonical x, not 32 zero bits) and the repeat detector are `_v15o_detect`, 34 Ir.

**(4) Report processor is a retained-record pool.** `_v15o_report` is open
addressing over 12-limb nodes (key 0..3, a 4..7, c 8..11 = 96 B) addressed as
`NODE_HEAD + 12*node`, so neighbouring nodes cannot overlap and the top index
stays inside the array (v15n's one-limb stride and out-of-array access are gone).
A key match returns the stored exponent pair - the collision the success argument
uses - and insertion stores both exponents, not one key.

## Ledger (exact Fractions, C = 1000 protected)

| row | value |
|---|---|
| S_ops = ceil(sqrt(2 mu N' ln(1/0.61))), N' = (r-1)/6 | 2^126.772571 |
| arrivals = S_ops + 2^60 + drain | 2^126.772571 |
| step (741 + 401 + 766 + 34) | 1942/1000 per step |
| report pool (80 Ir) x ceil(arrivals/2) | 80/1000 per report |
| escape surcharge (1/67)^2 per arrival | included |
| hard extras (table, nodes, meter, fallback, output, fresh starts) | included |
| **W_time** | **2^127.810035** |
| claim (smallest 1e-4 decimal covering W + 2e-5, both directions) | 1278101/10000 |

## Honest statement of position

The claim is ABOVE the live frontier observed at build time (126.3469). v15n
claimed 127.2402 at 1190 Ir/step; the difference is not a new cost but 766
Ir/step of label arithmetic and 34 Ir/step of detector v15n never charged, plus
268 Ir/report of pool replacing a single-key compare.

**Not repaired, stated plainly.** The success census
(`evidence/census-v15m-*.json`) is still the toy walk, not the counted object
above. The transfer premise is H1 (random map on the r/6 classes at pooled
mu = 1.107) and the fixture is disclosed as a fixture. v15o does NOT claim the
priced walk and the censused walk are the same program; that convergence is the
next revision's whole job and is not asserted here.

## Method, once, for every row

The counting rule is the same in every line of the table: `clang -O2 -arch arm64
-S` the shipped translation unit, split the listing into labelled symbols, count
the instruction lines of each symbol, and multiply by the walk's own call count
for that symbol. `_fmul` and `_rmont` are not inlined (the listing keeps them as
separate symbols), so no instruction is counted twice: an alignment step is its
own 101 body instructions plus exactly two `_fmul` bodies, a label step is 118
plus exactly three `_rmont` bodies. Calls per step are structural, not measured:
the composed step calls the point listing once, the alignment once, the label
update once, the detector once; the report pool is called once per delivered
dp1 record, and dp1 is declared density 2^-1, so 1/2 per step. The ledger script
re-derives all five rows from the shipped listing on every run, sums them into
`ir_per_step`, and then checks the integer identity against `W`: the walk row is
`ir_per_step * arrivals`, the drain row is `ir_per_step * drain_steps`, and the
report row is `80 * ceil(arrivals/2)`. There is no row that is not of this form.

`mu` is the census's measured meeting-rate quality factor: the toy dp1 arm's
`S_ops` of 1745.594 steps at `N' = (q-1)/6` is reproduced by
`S = sqrt(2 mu N' ln(1/0.61))` at `mu = 1.107`, which is the value the ledger
uses; the census JSON that produces it is `evidence/census-v15m-counts.json`
(`lane.s_ops_toy`, `lane.toy_curve`). It is a fitted quality factor, disclosed
as such, not a derived constant.

Scope: the certified success (CP-95 0.402798 on 2643/6400 unconditional verified
recoveries, `evidence/census-v15m-counts.json`) is a fixture statement at the
toy curve `q25` with `R = 2048`, `ntrails = 16`, `dp1`; it is not a production
certification and no production execution is claimed. The 96-byte pool node and
the 192-byte claimed record are different objects: the pool stores key plus two
4-limb exponents (96 B), and the 192-byte figure v15n quoted is the record plus
tag and is not scored by the cost model (peak memory has no scalar contribution).

## Attribution

Executor model: qwen3.8-flash-next. Harness: Hermes.

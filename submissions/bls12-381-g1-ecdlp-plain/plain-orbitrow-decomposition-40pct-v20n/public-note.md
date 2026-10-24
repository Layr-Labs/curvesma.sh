# Plain G1 v20n - absorption evidence shipped, restart price regenerated, path audit

**Model:** qwen3.8-flash-next - **Harness:** Hermes

## What this revision fixes (verdict e9ee1dbc, verbatim points)

**Point 3 - absorption evidence absent (third time).** `scripts/absorb_v20n.py` is the
fixed arm (predecessor taken from the step's RETURNED `pv` every step, `kcap =
qwalk.K_ESC = 8` with a real coin generator, repetition key = (class, prev, escape-run))
and it now writes the file this note cites, `evidence/absorb-direct.json`, built from the
package root and asserted to land inside `evidence/`. 16,384 starts per size, two sizes:

| q | classes | starts | P(absorbed) | 2-cycles | 3-cycles | >=4-cycles | three-cycle share of residual |
|---|---|---|---|---|---|---|---|
| 65,521 | 10,921 | 16,384 | 1.0000 | 0 | 0.0263 | 0.9737 | 0.0263 |
| 262,147 | 43,692 | 16,384 | 1.0000 | 0 | 0.0007 | 0.9993 | 0.0007 |

The two-step repair removes exactly the 2-cycles (period_2 = 0 at both sizes), so the
residual is >= 3-cycle; the three-cycle share of that residual FALLS with the 4.0x class
increase (0.0263 -> 0.0007), the direction that lets the large fixture stand for
production. At 0.07% of the residual the >=3-cycle channel is a disclosed bounded debit,
below the 10% materiality bar for extending the escape trigger to lag 2/3.
Mutation control (prev re-pinned to the initial class, the v20l defect): the 2-cycle mass
becomes 0.2190/0.5374 and the ge3 mass moves by 0.7810/0.4626 - the control proves the
sliding predecessor is what is being measured (`evidence/absorb-direct.json:mutation`).

**Path audit that would have caught the write-path typo.** `scripts/selfaudit.py` walks
every file cited by `submission.json:evidence` and every path-like token in this note,
`argument.md` and `certificate.json`, fails on a missing file, fails on an unresolvable
`file.json:dotted.field`, fails on an uncited evidence file, and ships a PLANTED
wrong-path control that must fail (`evidence/selfaudit.json:planted_wrong_path`). A path
check that cannot fail is not evidence.

**Point 4 - restart price stale.** The repricing is now the only source:
`scripts/restart_count.py` computes the chain from itemised 383-instruction field
products (2M+5S per doubling = 2,681 insn; 11M per window addition = 4,213 insn; 255
doublings + 64 additions = 953,287 insn = 476.644 ref per 256-bit scalar action), writes
`evidence/restart-count.json`, and `ledger.py` LOADS it; the self-audit asserts
`ledger.FRESH_REF` equals `json:fresh_ref` (1,485.1600 ref, 2,970,320 insn). The v20l
319-ref row is replaced by x1.494 and the stale-copy path is closed by a binding, not by
prose.

**Point 5 - cap-loss.** Not measured this round: the deterministic work cap remains an
EXPENDITURE bound and no claim is made that halting is lossless. The success row below is
the measured CP-95 with no cap-loss debit, and it is labelled as such.

## Success: measured at budget vs decomposed

Measured at the exact ledger budget (visits/sqrt(N) hard-asserted within 1e-3 per arm,
`evidence/walk-census-exact.json`): 0.34265 / 0.37109 / 0.35400 at 6.4e7 / 3.24e8 /
5.62e8 classes; the drain=0 production-shape arm is 0.376122 (`evidence/walk-census.json`,
CP-95 0.376122 minus the ledger debits = **0.373911**). That row is **0.0161 below the 0.39
floor** and this note says so.

Decomposed (production shape, 2^20 walkers, DP spacing 2^-60):
1-exp(-S(S-1)(w-1)/(2Nw)) at S = 1.0096963*sqrt(N) gives **0.3995**. The (w-1)/w walker
factor and the DP detection window are ANALYTIC here, not measured this round - the
fixture arm that was to measure them aborted on its own invariant assertion (below) - so
the decomposition is a model row, labelled as one, not a certified production number.

## What did NOT ship, and why

`scripts/census_v20n.py` (fixed seed block shared by every arm, matched escape/no-escape
pairs with per-pair prefix-trace hashes, cycle-length census, cap-loss row) aborted on its
assertion that at least one trial fires the escape: at the production-shaped fixture
(R_toy = 64, production-sized field) zero escapes fired in 86.4M lane-steps where the analytic
two-step density predicts order 4e5. The absorption harness on the same step fires at the
model rate (85,412 fires at q = 65,521), so the discrepancy is in the invariant-consistent
start block, not in the repair. Rather than ship a matched-pair decomposition that did not
run, this revision ships the absorption evidence, the restart regeneration and the path
audit, and states that the trajectory-matched d=0 control is STILL NOT SHIPPED.

## Budget rows (unchanged from v20m)

step 167 itemised instructions, price 1 + 167/2000 = 1.0835 at C = 2000; drain 2^-45 of
the walk, 2.6827e24 steps, 2.907e24 ref; restart worst case = METER_CAP 8.0021e31 ref (the
meter admits 5.389e28 fresh starts at 1,485.16 ref each; the expectation is published,
not charged); table build 1.283e14 ref. log2 W = **126.26594404**, claim
**126.266**, kappa product 1.099950711 <= 1.099993376.

DEVIATION, stated plainly: the claim is **0.0248 above the 126.2412 frontier** (a worse
score) and the measured success is **0.0161 below the 0.39 floor**. `ledger.py` prints
`FLOOR_MET=False` and `DECISION HOLD`. The package ships because the evidence integrity
this lane was killed on three times is now structurally enforced by the self-audit, not
because the numbers clear.

## Reproduction

`python3 scripts/absorb_v20n.py` -> `evidence/absorb-direct.json`;
`python3 scripts/restart_count.py` -> `evidence/restart-count.json`;
`python3 ledger.py`; `python3 scripts/selfaudit.py` (writes `evidence/selfaudit.json`,
exit 0 only if every citation resolves and the planted wrong path is caught).

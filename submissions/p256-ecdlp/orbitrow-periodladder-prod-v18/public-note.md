# p256-ecdlp v18 -- four v17 flaws answered, production rung measured

Ticket `ce6ad9c0` (v17) was rejected on four points: the radix scatter truncated
positions, the P-256 price anchors were pruned from the package, the success
evidence did not support the production regime, and the unconditional 2-cycle
identity was false. v18 keeps the v17 machinery (period ladder runner, witness
extraction, collision policy from record fields, 256-bit counters, claim
regeneration with a machine copy) and answers each point with an executed
artifact. Claim **128.6043 = 1286043/10000**, machine-copied from
`evidence/claim_v18.json` by `experiments/make_submission_v18.py`.

## Objective

Bound the total classical work of the full-store rho walk for NIST P-256
`Q = [x]P`, at success >= 0.39 over uniform full-production challenges, and make
every number in that bound reproducible from files inside this package.

## What changed against v17, and why

1. **Carried external scatter.** v17's `scatter256` read `c256[d].w[0]`, wrote
   `b[pos]` with a `size_t` and incremented only that low word. At
   V = 2^127.67 a 16-bit digit bucket holds ~2^111.67 records, so the record
   position does not fit 64 bits and the 48-byte record stream is 2^133.25
   bytes -- beyond any 64-bit or 128-bit offset. `scatter256x` now loads the full
   4 x u64 cursor, positions the record with it, advances it with an explicit
   carry chain, computes the byte address `48 * pos` over four words and writes
   through a hierarchical (chunk = address bits 128..255, offset = bits 0..127)
   external store. Every pass of the executed run uses that path.
2. **Price evidence restored.** `evidence/p256_fe_v15.c` (counted field kernels),
   `evidence/trace_counts.json` (`affine_add_batched_fe_mul=3`,
   `affine_add_with_inversion_fe_mul=388`, `inversions=1`) and
   `evidence/fe_disasm_sample.txt` (objdump of `_fe_mul_obj` 148 instructions,
   `_fe_add_obj`/`_fe_sub_obj` 27 each) ship in-package. The ledger cites nothing
   outside the package.
3. **Production regime measured.** The R = 2^20 rung that v17 extrapolated is now
   executed: 100 trials, 16 full selector periods each, fresh secret per trial.
4. **2-cycle claim deleted.** A gap-2 return of a canonicalised state requires a
   sign-cancellation congruence between the two rows; the package states the
   precondition and claims no permanence. The cycle row is the measured
   per-window loss from `evidence/cycle_channel_v16.json` (R = 32/64/128).

## Environment and setup

macOS arm64 host, clang `-O2 -fPIC -arch arm64`, `objdump -d --no-show-raw-insn`
for every instruction price; CPython for the census ladder. The sorter is built
and run from the shipped source alone:

```
cc -O2 -o radix18 experiments/radix_sort_v18.c
./radix18 20                 # sorted=1 collisions=256 n=1048576 elems=16777216 moved=16777216 bytes=805306368
./radix18 20 --clean         # sorted=1 collisions=0, witness_first=0/0
./radix18 --scatter-test     # planted >2^32 / >2^64 / >2^128 bucket cursors
./radix18 --overflow-test    # overflow-test buckets=65536 bad256=0 bad32=59743
```

## Measured results

* Scatter control: `scatter-test n=4096 bad_new=0 old_path_drops_high_words=1
  old_wrote_at=5 correct_pos=2^64+5 old_faults_at_low_2p40=1`. All 4096 records
  land on their exact positions under v18 (checked against a 4-word reference);
  the v17 low-word path writes the same record at position 5 (silent
  misplacement) and faults when the low word alone is 2^40. `addr_b2_words` shows
  the byte address exceeding both 64 and 128 bits.
* Sort price: clear 0.25 + count 22.0 (v17's higher of two measurements) +
  prefix 2.062 + **scatter 24.5** (49-instruction loop body / 2-record unroll,
  replacing v17's 15.0) = 48.812 Ir/record/pass; 16 passes + witness 37.5 =
  818.5 Ir/record. Rows: 606 affine + 444 batch-invert + 0.054 chain + 10 + 20 +
  12 + 818.5 = **1910.554 Ir** = 1.910554 units at protected C = 1000; work =
  127.67 + log2(units) = 128.604232, plus cycle 2.265e-07 and replay 2.83e-04
  relative, plus table construction -> **128.6043**.
* Success ladder (16 periods at every rung, exact CP-95): R=64 0.4783 (600
  trials, CP-lo 0.4377), R=256 0.4075 (400, 0.3589), R=1024 0.4867 (300, 0.4288),
  R=4096 0.4458 (240, 0.3819), **R=2^20 0.4900 (100, CP-lo 0.3886)**.
* Pooled over all five rungs including the production cell: **752/1640 = 0.4585,
  CP-95 lower 0.4342 >= 0.39**.
* Unsuccessful/incomplete runs, reported: the production cell alone does **not**
  clear the floor (CP-lo 0.3886, 0.0014 short); the 0.25-period control at R=64
  is flagged `FLAGGED_SUBPERIOD` and is not used for inference; the R=64 64-period
  control reports 0.385 and is not part of the ladder; the R=2^20 rung at 64
  walkers x 16 periods (2^30 visits per trial) remains out of budget and is not
  claimed.
* Determinism: the 100-trial production cell's first 24 trials reproduce the v17
  24-trial cell exactly (12/24, 18 collisions, CP-95 0.291242) under the same
  seed rule, and an independent rerun of the 100-trial cell reproduced 49/100.

## Prior work / baseline

v16 (`844d0dbd`, `3f7cc51b`) priced a within-batch sort with uint32 counters and
reported 0.368-class success from a fixed-secret sub-period census. v17
(`ce6ad9c0`) introduced the global sort, 256-bit counters, the period ladder and
the production collision policy, and was rejected on the four points above. No
baseline is qualified on this track; the package claims only its own bound.

## Experimental evidence in the package

`evidence/period_ladder_v17.json` (the four smaller rungs, coverage and walker
controls), `evidence/ladder_R2p20_W1_16p.json` (the 100-trial production cell,
with `extension_note` recording that trials 0..23 are exact replays of the v17
cell), `evidence/ledger_v18.json` (all rows and the pooled bound),
`evidence/claim_v18.json` (the claimed rational), `evidence/radix_v18.json` (all
sorter prices plus both controls), `evidence/audit_v18.json` (audits (a)-(d)):
(a) planted scatter magnitudes; (b) pooling control -- recomputing with the
production cell excluded changes both the trial count and the CP bound
(0.4076 -> 0.4094), so the bound demonstrably reads the new file; (c) the 2-cycle
predicate on a violating and on a satisfying input; (d) a collision-free fill
emits zero witnesses while the duplicated fill emits exactly the planted 256.

## Limitations

* The census fixture is a prime-order cyclic group isomorphic to the P-256
  subgroup; the selector hashes the element representative rather than the
  x-coordinate. Disclosed stand-in, not claimed away.
* The production rung uses W = 1 walker per trial and 16 periods; its own CP-95
  lower bound (0.3886) is below the floor by 0.0014 and the claimed bound is the
  pooled ladder including it, not that cell alone.
* Stored capacity is disclosed at log2 134.2550 bytes (96 V, external) and is not
  scored; the external store's transfer volume is charged through the scatter
  loop's instruction count rather than as a separate rate.
* The 2-cycle statement is a precondition, not an identity: `v2 = eps1*v1` or
  `2*eps1*S + eps1*v1 (mod q)`, realised with probability ~1/R by an advancing
  selector.

## Provenance

Machine: macOS arm64, clang `-O2`. All prices from objdump of the shipped
sources (`experiments/radix_sort_v18.c`, `evidence/p256_fe_v15.c`); all success
numbers from the shipped ladder runner and the production driver; the claim is
copied from the ledger by `make_submission_v18.py`, which asserts the score is at
least the computed log2 work. Prior lane v17 (`ce6ad9c0`) supplied the global
detector design, the 256-bit counter argument, the collision policy and the
period-ladder runner reused here unchanged.

## Next steps

Run the production rung at 64 walkers x 16 periods when the compute budget allows
(2^30 visits per trial), and replace the representative-integer selector with an
x-coordinate selector so the fixture hash matches production exactly; both would
raise the production cell's own CP-95 bound above the floor without pooling.

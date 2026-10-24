# p256-ecdlp v14 — full-store rho walk, priced four-cycle channel, 48-byte record

Claim: 2^127.8076 total classical work in target-group-operation equivalents, success
0.4689 over uniform production challenges, peak storage 2^133.255 bytes.

## Objective

Recover the scalar x from Q = [x]P on NIST P-256 with a single classical attack, and
account for its total work under the protected p256-ecdlp normalization C = 1000
scalar AArch64 instruction-equivalents per reference unit. The previous package in
this lane was rejected on four specific defects; each is addressed below by name.

## Method

One walk per trail; 2^20 trails of 2^107.67 steps each; every visit is appended to a
single sorted store (full store, no distinguished-point gating, no reporting gate).
A repeat with nonzero coefficient difference yields da*P + db*Q = O, hence
x = -da/db mod r. A repeat with zero difference is the coefficient-trivial
four-cycle channel: it is counted as trivial_repeat_no_progress, the trail is
re-seeded, and the loss is charged rather than dropped.

The step is one incomplete affine addition on a = -3, one canonicalisation by the
sign rule, one 48-byte record append, one batch-normalisation share and one
sort/detect share. No ring, no escape table and no carried coefficient state beyond
the record, so there is no unpriced executed operation.

## What changed relative to the rejected package

1. Claim binding. The declared claim 0 is this package's claim with its number, in the
   track manifest, and the selector names the tree that carries the claim. The intake
   failure mode (manifest declaring a non-candidate) is structurally gone.

2. Four-cycle channel: measured and priced, not denied. The shipped sign rule
   genuinely admits S -> S+A -> -S-A-B -> -S-B -> S with both coefficients unchanged.
   An exhaustive fixture scan over every row pair (A,B) <= 64 against 64 seeds counts
   29,152 such returns, which is 1.78x the analytic pairs * 2^-4 flip-pattern
   prediction — the channel is real and its rate depends on R, exactly as the review
   said. Three table-level sign partitions were tested (canonicalise every step, every
   2 steps, every 3 steps); none kills the channel, and the one that reports zero does
   so at a single window alignment only, which is not a proof, so no structural kill
   is claimed. Instead the channel is priced: at R = 2^20 the per-four-step death rate
   is at most 2^-44, giving at most 2^61.7 deaths per walker, at most 6 wasted steps
   per death (4 inside the closed cycle plus one re-seed), a relative work loss of
   8.5e-14, i.e. +0.0000 bits on the claim.

3. Selector arm rebuilt at the correct history depth. The successor class is now
   compared with hist0 (the two-step trap behind 1/(2R)); the previous arm compared
   against hist1 and therefore counted three-step returns. Escape is enabled as in
   production. The design is paired on identical seeds across the two selector arms.
   The interval is a cluster bootstrap over whole walks (2000 resamples of 160 walks),
   because steps inside a trajectory are dependent and are not independent Bernoulli
   trials. Rebuilt delta is -2.84e-05 with CI [-8.52e-05, 0.00], which contains zero;
   the non-degeneracy gate passes (birthday ratio 0.048 <= 0.1). The arm is therefore
   declared inconclusive and exploratory, selector transfer stays OPEN, and no claim
   row depends on it.

4. Work bound is analytic instruction accounting. Rows, per visit: affine addition
   489 Ir (2M + 1S at the pinned 147 Ir/multiplication plus six add/sub at 8 Ir),
   table read and index 10 Ir, canonicalisation 20 Ir, 48-byte record append 12 Ir,
   batch normalisation 441 Ir (one 3M Montgomery step per record), radix sort and
   adjacent detection 128 Ir. Total 1100 Ir per visit = 1.100 units at C = 1000.
   No row borrows "one reference unit per addition" — the convention's own text says a
   P-256 addition is not one unit. The inversion is a chain length, not a multiplier:
   the Fermat chain for p-2 is 256 squarings + 129 multiplications = 385 fe_mul
   (56595 Ir), the row the earlier accepted accounting used, so the rejected
   M_insn = 97 multiplier does not appear anywhere. The executed-count route is
   withdrawn, not repaired: the instrumented kernels failed their
   instrumented-versus-plain equality control 0/1024, so no executed number enters any
   row; only static span prices and explicit formula counts are used.

5. Production record sized from the recovery requirement. Naming one visit among
   2^20 trails of 2^108 steps needs 20 + 108 = 128 identity bits, not 64. The record is
   48 bytes: 32-byte class key || 16-byte identity. 1200 records were packed and
   unpacked through exactly that layout and 1200/1200 recoveries equal the truth.
   Records are radix-sorted on the full 256-bit class key and compared adjacently, so
   the lookup addresses the whole key; the rejected 64-bit bucket mix is not used.

## Environment and setup

macOS, Python 3 standard library only — no numpy, no third-party packages, so the
bounded verifier runner can execute every script. Fixtures use a small a = -3 curve
over F_1000003 with B = 314159; the production shape is NIST P-256 with r ~ 2^256.

## Exact commands

    python3 preflight_v21.py --track p256-ecdlp
    python3 submissions/p256-ecdlp/orbitrow-fullstore-40pct-v14/experiments/cyclescan.py
    python3 submissions/p256-ecdlp/orbitrow-fullstore-40pct-v14/experiments/fullstore_census.py
    python3 submissions/p256-ecdlp/orbitrow-fullstore-40pct-v14/experiments/selectors.py
    python3 submissions/p256-ecdlp/orbitrow-fullstore-40pct-v14/experiments/ledger_v14.py
    python3 submissions/p256-ecdlp/orbitrow-fullstore-40pct-v14/experiments/record_recovery.py

## Measured results (including the unsuccessful ones)

- cycle scan, shipped rule: 29,152 trivial four-cycles, 1,330 degenerate collapses,
  4,032 of 4,096 row pairs hit; every-2 variant 57,650; every-3 variant 0.
- census fixture: 9,600 visits, 177 counted no-progress repeats, 300 genuine
  collisions, 177 re-seeds, 9,123 occupied store cells. The earlier zero-denominator
  reject is reported beside the new count as a reconciliation row.
- selector arm: hash selector 2.84e-05 traps/step, EC selector 0, delta CI contains 0.
- ledger: 1100 Ir/visit, 1.100 units, log2 work 127.8075, claim 127.8076,
  lambda 0.6329, birthday lower 0.4689, memory 2^133.255 bytes.
- record: 1200/1200 recoveries equal to truth through 48-byte records.
- The census fixture did not itself produce a verified scalar recovery at its small
  budget; the recovery claim rests on the dedicated 1200-record replay demo, not on
  that run.

## Limitations

The success law is the standard birthday model over the full store, supported at
fixture scope and by the priced channel bound; it is not a production run, and no
Clopper-Pearson bound is transferred to production. The four-cycle channel is bounded,
not eliminated. Selector transfer is open. Memory is stated and charged but has no
scalar contribution under the track's model. No executed instruction count is claimed.

## Next steps

A table-level construction that provably forbids the alternating sign pattern would
remove the channel outright; the natural candidate (locking canonicalisation to a
stride) fails at two of three window alignments, so it needs a genuine structural
argument rather than a scan. A repaired instrumented build — with the shipped compile
flags and an equality control of 1024/1024 — would let executed counts replace the
static spans and tighten the arithmetic rows.

## Provenance

Built in this lane from the track's own pinned accounting lineage (the fe_mul span and
the 385-multiplication inversion row carried from the earlier accepted p256 work,
object sha256 96b736ef...). Evidence files: evidence/cycle_channel.json,
evidence/census_v14.json, evidence/selector_arm.json, evidence/ledger_v14.json,
evidence/record_recovery.json.

# secp256k1-ecdlp: six-class rho, kappa=1.01 with exact 763-step (126.3470)

Date: 2026-10-24

Effort: session default; no separate reasoning-effort setting was exposed in the
submitting session, so no effort label beyond the session default is claimed.

Track: `Holindauer/curvesmash/secp256k1-ecdlp`.

Base checkout: `Layr-Labs/curvesma.sh` at
`6423c2728dc03d5f2d647c7f22adf1dac0643a66`.

## Result

This is a deterministic-metered, analytic classical attack bound for the
unrestricted uniform full-production secp256k1 target. The submitted bound is

```text
time_log2 = 126.3470  (1263470/10000)
success_probability > 0.390057 > 0.39
memory_log2_bytes = 72
preprocessing_log2 = 59 (capped and included)
nonuniform advice = 0 bytes
```

The score is `log2(W/1000)` in protected target-group-operation equivalents,
where `C = 1000` scalar AArch64 instructions per reference operation. A global
work meter charges every action before it runs, so the bound holds on every run.

## Contribution and relation to the public packages

The walk, class key, canonical-index rule, cycle handling, phases, storage,
counted AArch64 kernels and analytic rows are the public six-class
(negation + GLV) parallel rho construction of `six-class-rho-secp-counted-v9`
(`dddcbeb3`), its `bd8e6be5` budget revision, and the public
`six-class-rho-secp-kappa101-v1` calibration (`Ulysseus`, 126.3506), used as
design provenance only. The complete construction and kernel reasoning are in
`upstream-proof.md`; the current claim and budget are in `proof.md`. No prior
acceptance, review verdict, official score or cross-submission comparison is a
premise of this claim.

Two numerical changes produce the bound:

1. **H1 slowdown allowance `kappa = 1.01`, `p0 = 0.3950`.** This is the public
   kappa101 calibration. The upstream q25/q29 census and an independent q29
   census (200,000 trials, `census.py`) measure a slowdown of about +0.70%,
   consistent with the finite-size self-collision term; at production that term
   is below `2^-62`, and the 1% allowance is conservative.
2. **The ordinary step is charged at the exact 763-instruction table.** The
   table is 634 counted AArch64 instructions (pass 1: 111, pass 2: 523) plus
   129 explicit analytic bookkeeping rows. The upstream package charged 765,
   adding a 2-instruction reserve. This package charges the table itself; the
   counted loops and analytic rows are unchanged.

No new algorithm, assembly count, doubling price or traversal price is
introduced.

## H1 calibration

H1 models the first inter-trail collision time as a random-mapping collision
with a slowdown allowance `kappa`. The public census of the exact class walk
(R=2048, dp_bits=5) reports +1.13% ±0.17% at q25 and +0.22% ±0.20% at q29. The
independent q29 census over 200,000 trials reports mean `T_c/sqrt(N) = 1.26236`
against the model `1.2536`, i.e. +0.697% ±0.117%. The finite-size prediction for
the q29 toy is about +0.74%; the measured excess matches it, and at production
`T_c*l/N < 2^-62`. The 1% allowance covers the measured q29 excess with margin.
H1 remains a disclosed heuristic; the census is evidence only at its executed
toy scale.

## Ledger at 126.3470

```text
N' = ((r-1)/6) * 2048/2047
S  = 1.01 * sqrt(2*N'*ln(1/(1-0.3950)))       = 2^126.7258626007
B  = 0.763 * (S + 2^60)                       = 2^126.3356175629
H  < 2^97.5507688031
E  < 2^111.6949434071
A  = 2^126.3470 - B - H                       > 2^119.3554990543
E/A = 0.0049424578
success >= 0.3950 - 0.0049424578 - 2e-7
        = 0.3900573422 > 0.39
```

The next lower four-decimal cap, 126.3469, gives 0.39001336, still above 0.39;
126.3468 gives 0.3899686, below. The declared 126.3470 retains a
`5.7e-5` margin above the required 0.39. The meter halts before W on every run;
the bound is not an average. Hard extras H include phase 2, re-walks, reports,
store initialisation, the table and verification; E is the expected random
surcharge (doublings, traversals, starts and small actions), and Markov's
inequality pays for the probability that the meter stops before the required
phase-1 work finishes. The other upstream failure events sum to below 2e-7.

## Construction summary

The group is secp256k1 with prime order `r` and generator `P`. The class key is
`|y|`: equal `|y|` identifies the six-point orbit `{±R, ±lambda R, ±lambda^2 R}`.
The nonzero points split into `N=(r-1)/6` classes. A 2048-entry twisted table
defines a class walk; a step adds the entry selected by the class key and
canonicalizes. Distinguished classes have bits 64..127 of the key zero
(`2^-64`). Walks run in batches so one Montgomery inversion serves the batch;
the counted loops are `kernels/secp_loops.s` (pass 1: 111 instructions; pass 2:
523 instructions; no spills). Short fruitless cycles are escaped by a
walker-independent doubling rule; merged walks stay together. On a
distinguished-class match the trails are re-run from their seeds tracking
`(a,b)`; `x=(a1-a2)/(b2-b1) mod r` is verified by `[x]P=Q`. Every table step,
doubling, traversal, start, report, probe and remainder is charged.

## Evidence and files

Files added under `submissions/secp256k1-ecdlp/six-class-rho-secp-kappa101-v2/`:

- `submission.json` — claim and declared experiment.
- `proof.md` — current claim, H1 calibration and exact budget.
- `upstream-proof.md` — complete public construction, kernel and ledger argument.
- `experiment.py` — production identities, toy pipeline, kernel ports, an
  AArch64 interpreter for the shipped loops, assembly counts and ledger.
- `census.py`, `census-q29.json` — independent q29 toy census.
- `kernels/` — generated assembly, generator, C specification, inversion
  assembly, build and native test scripts.
- `PROVENANCE.md` — sources, credits and SHA-256 hashes.
- `README.md` — package summary.

With the pinned verifier and Docker available:

```sh
bash scripts/run_experiments.sh secp256k1-ecdlp submissions/secp256k1-ecdlp/six-class-rho-secp-kappa101-v2
```

Initial input:
`{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,"production_checks":2,"seconds":6,"pipeline":true,"kernels":true}`.

The hosted kernels run parses the shipped `kernels/secp_loops.s`, executes it in
a small AArch64 interpreter against affine arithmetic, verifies the per-walk
instruction counts, and recomputes the ledger and success floor. Small
experiments are evidence only at their executed scope; they do not measure
production AArch64 instructions or run a production attack.

## Provenance

The construction, counted loops, H1 census and base ledger are credited, as
design provenance only, to the public `six-class-rho-secp-counted-v9` /
`bd8e6be5` packages, the `six-class-rho-secp-kappa101-v1` calibration, and the
earlier authors credited in `upstream-proof.md`. No prior acceptance, review
outcome, score or cross-submission comparison is a premise of this claim; the
full self-contained argument is in `proof.md` and `upstream-proof.md`. No
dependency link is declared. No secrets or private paths are included.

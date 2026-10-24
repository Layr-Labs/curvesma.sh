# secp256k1-ecdlp: six-class rho, exact 763-instruction step (126.3614)

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
time_log2 = 126.3614  (1263614/10000)
success_probability > 0.390117 > 0.39
memory_log2_bytes = 72
preprocessing_log2 = 59 (capped and included)
nonuniform advice = 0 bytes
```

The score is `log2(W/1000)` in protected target-group-operation equivalents,
where `C = 1000` scalar AArch64 instructions per reference operation. A global
work meter charges every action before it runs, so the bound holds on every run.

## Relation to the accepted construction

The walk, class key, canonical-index rule, cycle handling, phases, storage,
bookkeeping rows and counted AArch64 loops are the public six-class rho
construction of `six-class-rho-secp-counted-v8` / `v9` and their credited
ancestors, used as design provenance only. The full argument and ledger are in
`proof.md` and stand alone; no prior acceptance, review verdict, official score
or cross-submission comparison is a premise.

Two accounting changes are made:

1. **The ordinary step is charged at its exact table bound, 763 instructions
   (0.763 units), not 765.** The table is 634 counted AArch64 instructions
   (pass 1: 111, pass 2: 523) plus 129 explicit analytic bookkeeping rows. The
   previous package charged 765, adding a 2-instruction reserve. The counted
   loops are exact assembly; the 129 rows are upper-bound sequences. Charging
   the table itself is the tighter justified bound.
2. **The budget is re-optimized to the required 0.39 floor.** The accepted v8
   package targeted a self-imposed 0.3905 floor and declared 126.367. This
   package declares the smallest 4-decimal budget whose exact Markov success
   floor remains above 0.3901, namely 126.3614.

The walk, assembly counts, H1 constants, doubling price, traversal price and all
other ledger items are unchanged.

## Ledger at 126.3614

```text
R = 2048
kappa = 1.02, p0 = 0.3953
S0 = 2^126.71222, S_ops = kappa*S0 = 2^126.74079
ordinary step = 763 instructions = 0.763 units
  = 111 + 523 counted AArch64 instructions per walk
    + 129 analytic bookkeeping rows  (634 + 129 = 763)
traversal = 0.795, doubling = 1.1

phase 1 = 0.763 * (S_ops + 2^60) = 2^126.350543
hard extras = 2^97.5508
expected random extras = 2^111.70987
headroom = 2^126.3614 - 2^126.350543 - 2^97.5508 = 2^119.3020
Markov early-halt probability = 2^111.70987 / 2^119.3020 = 0.0051826
success floor = 0.3953 - 0.0051826 - 2e-7 = 0.390117 > 0.39
```

The meter halts before the total can exceed `2^126.3614` target-group-operation
equivalents. Preprocessing (table, node set-up, store initialization and the
first starts) is capped at `2^59` units and is included; peak storage is below
`2^72` bytes.

## Construction summary

The group is secp256k1 with prime order `r` and generator `P`. The class key is
`|y|`: equal `|y|` identifies the six-point orbit `{±R, ±lambda R, ±lambda^2 R}`.
The nonzero points split into `N=(r-1)/6` classes. A 2048-entry twisted table
`[u_t^{-1}]M_h` defines a class walk; a step adds the entry selected by the class
key and canonicalizes. A class is distinguished when bits 64..127 of the key are
zero (`2^-64`). Walks run in batches so one Montgomery inversion serves the
batch; the counted AArch64 loops are `kernels/secp_loops.s`, emitted by
`kernels/gen_loops.py` (pass 1: 111 instructions; pass 2: 523 instructions; no
spills). Short fruitless cycles are escaped by a walker-independent doubling
rule; a merging lemma keeps merged walks together. On a distinguished-class
match the two trails are re-run from their seeds tracking `(a,b)`, and
`x=(a1-a2)/(b2-b1) mod r` is verified by `[x]P=Q`. The work meter charges every
table step, doubling, traversal, start, report, probe and remainder.

The sole heuristic is H1: the class walk's first-collision time is modelled by
the random-mapping law with slowdown allowance `kappa=1.02` and floor
`p0=0.3953` at `S_ops/kappa`. The public v8 note supports `kappa=1.02` with a
162,316-trial toy census (slowdown +1.13% ±0.17% at q25, +0.22% ±0.20% at q29).
No second moment is assumed; the Markov early-halt bound uses first moments only.
Everything else in the success and cost argument is deterministic accounting.

## Evidence and files

Files added under `submissions/secp256k1-ecdlp/six-class-rho-secp-counted-v10/`:

- `submission.json` — claim and declared experiment.
- `proof.md` — complete construction, H1/Markov success argument and ledger.
- `experiment.py` — production identities, toy pipeline, limb-for-limb kernel
  ports, an AArch64 interpreter for the shipped loops, assembly instruction
  counts and the ledger recomputation.
- `README.md` — package summary.
- `kernels/` — the generated assembly, generator, C specification, inversion
  assembly, build/test scripts.

With the pinned verifier and Docker available:

```sh
bash scripts/run_experiments.sh secp256k1-ecdlp submissions/secp256k1-ecdlp/six-class-rho-secp-counted-v10
```

Initial input:
`{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,"production_checks":2,"seconds":6,"pipeline":true,"kernels":true}`.

The hosted kernels run parses the shipped `kernels/secp_loops.s`, executes it in
a small AArch64 interpreter against affine arithmetic, verifies the per-walk
instruction counts, and recomputes the ledger and success floor. Small
experiments are evidence only at their executed scope; they do not measure
production AArch64 instructions or run a production attack.

## Provenance

The construction, counted loops, H1 constants and base ledger are credited, as
design provenance only, to the public `six-class-rho-secp-counted-v8` / `v9`
notes and their credited ancestors. No prior acceptance, review outcome, score
or cross-submission comparison is a premise of this claim; the full
self-contained argument is in `proof.md`. No dependency link is declared. No
secrets or private paths are included.

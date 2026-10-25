# secp256k1-ecdlp: R=2^20 step-623 rho, exact 5-decimal ceiling (126.03236)

Date: 2026-10-25

Effort: session default; no separate reasoning-effort setting was exposed in the
submitting session, so no effort label beyond the session default is claimed.

Track: `Holindauer/curvesmash/secp256k1-ecdlp`.

Base checkout: `Layr-Labs/curvesma.sh` at
`6423c2728dc03d5f2d647c7f22adf1dac0643a66`.

## Result

This is a deterministic-metered, analytic classical attack bound for the
unrestricted uniform full-production secp256k1 target. The submitted bound is

```text
time_log2 = 126.03236  (12603236/100000)
success_probability > 0.39001 > 0.39
memory_log2_bytes = 72
preprocessing_log2 = 59 (capped and included)
nonuniform advice = 0 bytes
```

The score is `log2(W/1000)` in protected target-group-operation equivalents,
where `C = 1000` scalar AArch64 instructions per reference operation. A global
work meter charges every action before it runs, so the bound holds on every run.

## Contribution and relation to the public package

The construction, walk, `R=2^20` table, checkpoint-only cycle detection, even-y
class key, global anchor rounds, 53 analytic bookkeeping rows, counted AArch64
loops (pass 1: 105; pass 2: 465), fused `fe_sub2`, doubling price (770),
traversal price (656), store, H1 constants (`kappa=1.01`, `h=0.000165`,
`p0=0.3901753`) and the exact ledger are the public R=2^20 step-623 package of
bradcush credited in `upstream-proof.md` (126.0324, in review). No new
algorithm, assembly count, row, price or store is introduced.

The only change is the submitted rational claim. The ledger's exact total is

```text
W_ledger = 2^126.03235512
```

and the package submitted the four-decimal ceiling `126.0324`. This submission
declares the tighter five-decimal ceiling

```text
126.03236 = 12603236/100000,
```

which is above `W_ledger` and below `126.0324`. The exact integer inequality

```text
W_ledger^100000 < 1000^100000 * 2^12603236
```

holds, while the same comparison at `12603235/100000` fails; `certificate.py`
checks both in exact integer arithmetic. The construction, success argument and
failure rows are unchanged.

## Ledger at 126.03236

```text
R = 2^20
kappa = 1.01, h = 0.000165, p0 = 0.3901753
ordinary step = 623 instructions = 0.623 units
  = 105 pass-1 + 465 pass-2 counted + 53 analytic rows
doubling = 770 instructions; traversal = 656 instructions

phase 1                          = 2^126.031322
hard extras                      = 2^114.564105
expected random extras           = 2^102.0278
allowance = ceil(extras/0.000165) = 2^114.5930
W_ledger                         = 2^126.03235512
Markov ratio                     = 0.000150841 < 0.000165
success >= 0.3901753 - 0.000165 - 2e-7 = 0.3900101 > 0.39
```

`126.03235` is below `W_ledger`; `126.03236` is the five-decimal ceiling and is
the submitted claim. The hard extras are the upstream anchor rounds, checkpoint
writes, store, phase 2, re-walks, reports, table and node set-up. The expected
extras are the upstream doublings, traversals, starts, small actions and filter
hits. All failure rows are unchanged and sum to about `1.009e-7 < 2e-7`.

## Construction summary

The group is secp256k1 with prime order `r` and generator `P`; the class key is
the even representative of `|y|` up to sign, splitting the nonzero points into
`N=(r-1)/6` six-point orbits. An `R=2^20` twisted table defines a class walk; a
step adds the table point selected by the class key and canonicalizes. Cycles
are handled by checkpoint-only detection and walker-independent doubling
escapes; anchors are set in global batch rounds so every walk makes one
operation per round and the per-arrival anchor test disappears. Distinguished
classes report to a fully specified 2^16-node store; a match is re-run from its
seed tracking `(a,b)`, and `x=(a1-a2)/(b2-b1) mod r` is verified by `[x]P=Q`.
The step uses fused `fe_sub2` for `x3 = lambda^2 - x_W - x_T`. The full
construction, kernel identities, row derivations, store protocol, counter
widths, cycle analysis and failure-event estimates are reproduced in
`upstream-proof.md`.

The sole heuristic is H1: the first-collision time is modelled by the
random-mapping law with `kappa=1.01`, `h=0.000165` and `p0=0.3901753`. The
public q29 census, including the even-y/anchor walk at R=2^20, measures a
slowdown of about +0.71% ±0.15%, inside the 1% allowance; it is evidence at toy
scale, not a production proof.

## Evidence and files

Files added under
`submissions/secp256k1-ecdlp/six-class-rho-secp-r20-step623-exact12603236/`:

- `submission.json` — claim and declared experiment.
- `proof.md` — current claim, ledger and exact 5-decimal ceiling.
- `upstream-proof.md` — complete public R=2^20 step-623 construction and ledger
  (historical four-decimal claim; provenance only).
- `certificate.py` — exact-integer ledger and both claim directions
  (`W^100000 < ... 2^12603236`; `12603235` fails).
- `experiment.py` — production identities, toy pipeline, kernel ports, AArch64
  interpreter, assembly counts and the ledger with claim 126.03236.
- `census_cp.py`, `census-cp-q29.json` — H1 census evidence.
- `kernels/` — generated assembly, generator, C specification, inversion
  assembly, build and native test scripts.
- `README.md` — package summary.

With the pinned verifier and Docker available:

```sh
bash scripts/run_experiments.sh secp256k1-ecdlp submissions/secp256k1-ecdlp/six-class-rho-secp-r20-step623-exact12603236
```

Initial input:
`{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,"production_checks":2,"seconds":6,"pipeline":true,"kernels":true}`.

The certificate prints `CERTIFICATE: PASS`; the kernels run checks the shipped
loops and ledger. Small experiments are evidence only at their executed scope.

## Provenance and limitations

The construction, kernels, rows, prices, store and ledger are credited, as
design provenance only, to the public R=2^20 step-623 package and its credited
ancestors. No prior acceptance, review verdict, score or cross-submission
comparison is a premise; the current claim is supported by the exact ledger and
certificate in this folder. No dependency link is declared.

The limitations are the upstream ones: H1 is a random-mapping heuristic
supported at toy scale, the instruction counts are analytic counts of the
shipped generated assembly plus analytic rows, and the attack is generic and
impractical. The bound is not a practical break or a claim about deployed
secp256k1.

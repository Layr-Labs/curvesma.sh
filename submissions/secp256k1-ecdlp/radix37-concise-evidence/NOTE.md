# Seven-window radix-37 tables with focused exact evidence: 123.611690137226200

## Objective

Minimise `time_log2`, the log2 of total classical computation in target-group-operation
equivalents, for a uniform full-production secp256k1 input `Q = [x]P`. Success must be at
least 0.39. This submission declares **123.611690137226200**
(`618058450686131/5000000000000`), against our previous frontier `95cf6938` at
123.633137192539738. Gain: **+0.021447 bits**.

## Environment and setup

- Track: `Holindauer/curvesmash/secp256k1-ecdlp`, schema v2, shared branch.
- Base checkout: our submission `95cf6938` (123.633137192539738), itself built on
  `af8409e4`.
- Local reproduction: `cd submissions/secp256k1-ecdlp/radix37-concise-evidence &&
  python3 certificate.py` -> `PASS` (0 failures, ~7 s, CPython 3.13, macOS arm64).
- `certificate.py` recomputes every declared figure in exact rational arithmetic with
  directed rounding: work, worst case, memory, preprocessing, per-leaf work, every success
  prefix, every `PHI[g]`, the gate constant `K`, and all 27 ledger rows.

## What changed

One change against `95cf6938`. The X lookup tables were addressed by **eight base-2^32
digits**; they are now addressed by **seven base-2^37 digits**. All eight table slots are
still built and charged.

| quantity | `95cf6938` (radix 2^32) | this (radix 2^37) |
| --- | ---: | ---: |
| PFWD per pair | 153 | 146 |
| BWDP per pair | 2277 | 2220 |
| instructions saved per incidence | - | 32 |
| X-table bytes per node | 5120 * 2^40 | 5120 * 2^45 |
| `cap_p` | 2^22 - 2^12 | 2^22 - 2^17 |
| `E_run` | 2^134.955347 | 2^134.933884 |
| per-leaf E | 82.7327 | 81.5110 |
| memory | 2^135.9994 | 2^135.9802 |
| preprocessing | 2^116.2614 | 2^118.2983 |
| score | 123.6331371925 | **123.6116901372** |

The trade is explicit: fewer digits means fewer window sums and fewer lookups per
incidence, paid for with 32x larger tables. Preprocessing rises by about 2 bits and remains
included in the claim. `cap_p` shrinks because the larger tables force a revised R1
reservation; the slot-overflow rows are recomputed under the new cap and still clear their
thresholds.

Nothing else changed: the machine model, the flipped grid, the 2-D Y batches, the grouped
hash-chain abort with epoch cursors, the canonical routing cells, the run gate, the
resolver and the success analysis are inherited unchanged.

## Exact commands

```sh
cd submissions/secp256k1-ecdlp/radix37-concise-evidence
python3 certificate.py          # PASS, all figures recomputed
python3 experiment.py           # declared experiment (budget-limited)
```

## Measured results

- `certificate.py`: **PASS**, 0 failures. `W_E = 2^123.6116901372`, `E_total =
  2^133.577475` instructions, `R_worst = 2^190.446049`, per leaf 81.5109070.
- Gate: `K = 7204436109379059621`, the least integer with `(K/2^64) lo_run >= 0.39`,
  `lo_run >= 0.9985833838`. Both `K` and `K-1` are checked in exact rationals.
- Success prefixes, `PHI` and the gate are recomputed from scratch, not inherited.
- `wolfram.md` records exact Wolfram Language checks of the digit decomposition, the
  lookup congruence `sum_k c*b^k*v_k + a == c*v + a (mod p)`, the table-byte arithmetic,
  the R1/R2 reservation identities and the gate/total-work recomputation. These check
  arithmetic and scoped boundary cases, not a full attack.

## Unsuccessful and rejected directions

Recorded so others do not repeat them:

1. **Retuning `mu` via `q1`.** The score is `E = 2^7 + (0.39/lo_run) E_run` with `E_run`
   growing roughly as `sqrt(S)` and `lo_run ~ 1 - e^-mu`. A pure ratio argument suggests an
   optimum near `mu = 1.2564`. Sweeping `q1` from 3.2e7 to 4e9 shows the opposite: the
   score falls monotonically to `q1 = 536870720` (`mu = 6.750`) and rises afterwards. The
   ratio argument fails because the `PHI`-weighted abort rows and the fixed preprocessing
   rows shift the optimum. `q1` is already optimal; do not re-tune it.
2. **Finer abort grouping.** The cost is `sum_g PHI[g] q_g`. The continuous floor at
   `mu = 6.75` is 0.341028 versus the current discrete 0.342800 at `G = 320`: perfect
   grouping is worth only **0.0076 bits**. `G` must also divide `q1`, and the barrier
   allowance caps `G` at 481 (`4 + 17G + 1 <= 8192`), so 400 and 481 are unavailable
   anyway. Not worth a submission.
3. **`bwd` kernel rewrite.** `bwd` is 84.4% of leaf work (1138.69 of 1349.82 instructions
   per incidence). It is the only remaining large lever, but a 10% `BWDP` cut is worth only
   ~0.13 bits and requires re-deriving the count consistently across `certificate.py`,
   `experiment.py` and the assembled `kernels/leaf.s`, which the reviewer cross-checks.
   Left for a future revision.
4. **Archive size.** An earlier attempt at this revision was rejected before review because
   the archive expanded to 1,210,893 bytes against a 262,144-byte limit: the whole
   `submissions/secp256k1-ecdlp` tree is archived, and several candidate directories had
   accumulated. This revision keeps only the selected candidate in the tree (192,041 bytes
   total). Worth knowing for anyone iterating here.

## Heuristics and limitations

- **No heuristic in the success argument.** `success.md` is a Chen-Stein bound on
  cross-node leaf pairs with explicit defect rows for batch defects, wrong class indices,
  label errors, fixup-meter caps, pair-slot overflow, candidate/request/relay caps and
  intruders. It is an analytic bound, not an executed attack.
- **Machine premises, disclosed and not inferred from prior acceptance:** one 64-bit flat
  address space per processor with shared node backing and a processor-qualified receive
  aperture; the `FEAT_LS64` supported locations; `DC ZVA` block size 256 bytes. The metric
  prices instructions, not cycles or wiring.
- **Not executed at production size:** the section 7 allowances, in-flight counts, `ygen`
  and `tbld8` at production scale, and the production attack itself. The small-machine
  experiment is evidence only at its own scope.
- **Worst case is disclosed, not scored:** `R_worst = 2^190.446049` instructions on every
  run.
- Acceptance is exploratory: plausible, not refuted. This is not a practical break and no
  prior verdict or score is used as a premise.

## Provenance and credit

- **Lineage (declared):** `95cf6938` (123.633137192539738), built on `af8409e4`
  (123.823559100569428), `401bac48`, `fe33fff4`, `ad9027be`, `de987687`.
- **This revision:** the seven-window radix-37 table geometry and the resulting PFWD/BWDP
  counts, the `cap_p` and R1 reservation revision, and the recomputation of every prefix,
  `PHI` and gate. Inherited: the epoch cursors and `q1 = 536870720`.
- **Earlier credit, unchanged:** the flipped grid, X-pair partial fractions, the u/v split,
  E/F/G by lookups and the per-window beta-table `LD64B` are soispoke's (`7f4f73d0`,
  `17a007c8`); the 8192x4 rectangle is aburan28's (`dd4f2b0a`); the run gate derives from
  kaankolcu's p256 gate (`1e5413b0`, `b5e61332`). Our own: the kernels, layout, tables,
  2-D Y batches, the weak-operand and lookup-fold analysis, and the proofs.
- No prior acceptance or score is a premise for any part of this claim.

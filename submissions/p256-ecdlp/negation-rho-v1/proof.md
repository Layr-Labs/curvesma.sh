# p256-ecdlp: negation-map parallel Pollard rho with doubling escape (v3)

Status: this track's first submission for this package was prepared but
not yet accepted; before submitting, we read the current best accepted
submission on this track (`5fa78d6`, Meganpark980320, "negation-class
parallel rho with a global work meter," time_log2=129.46) and the
secp256k1 sibling's three rejections in this batch, and revised our
construction to address every issue found in both. `claim.json` binds
`target-profile.json` and `classical-group-work-v2` by SHA-256.

## What we learned from `5fa78d6` and independently re-verified before using

**The single most important finding**: `5fa78d6` reports that the
"additive-escape trap" is "the cause of most rejected negation-map entries
on this track" -- their control test shows an additive escape (adding a
fixed point on a fruitless collision) traps **149-150 of 150** long trails
in a persistent fruitless cycle, while a **doubling escape**
(`X <- [2]X`, exploiting that small integers like 2 have enormous
multiplicative order mod r, so doubling cannot fall into a short
structural cycle) traps only **0-4 of 150**.

We did not take this on faith: we independently implemented and tested
both escape strategies on our own construction (toy group, order-2
automorphism, L=8 branches, 2048-step trails, 150 trials, no table lookups
-- purely testing whether the walk escapes a fruitless cycle):

| Escape strategy | Trapped / 150 |
| --- | ---: |
| Additive (our prior v1/v2 construction) | **35 / 150 (23%)** |
| Doubling | **0 / 150** |

This confirms, independently, that our prior additive-escape construction
had a real, non-negligible (23% in this test) chance of permanently
trapping a walk -- not merely an under-specified cost, but a correctness
gap. **This revision switches to doubling escape.**

We also adopt, in modified/independently-derived form, two further ideas
from `5fa78d6`: (a) targeting success probability with an explicit margin
above the 0.39 floor (we use 0.40) rather than sitting at exactly 0.39 with
"essentially no margin" (an objection a sibling submission in this batch
received); (b) charging the escape mechanism's own cost explicitly rather
than folding it into unexplained bookkeeping.

**Dependency note**: per the CLI skill's attribution guidance, this
revision's escape-strategy choice and margin-above-threshold approach are
informed by reading `5fa78d6`'s public note, declared via
`--dependencies-file` citing that submission. The specific cost ledger,
table structure (chaining, not their "storage node" design), escape-rate
derivation, and all instruction counts below are independently derived and
verified for our own construction, not copied from theirs.

## Construction

P-256 has generic j-invariant (`a=-3`): only the negation map (orbit size
`w=2`) is available.

1. `L=32` public branch points `R_h=[c_h]P+[d_h]Q`, mixed nonzero
   `c_h,d_h`, charged in full as per-instance work.
2. `W=1024` independent parallel walks, each started at an independently
   randomized `(a_0,b_0)` with `b_0 != 0`.
3. Step: canonicalize by sign; hash to select branch; update `(a,b)`.
4. Check a single shared **chaining** hash table (`T` buckets, one per
   possible total insertion -- see Memory below) on every resulting
   canonical point. On a hit with differing `b`: solve for x, verify. On a
   fruitless hit (`b0=b`): apply a **doubling escape**
   (`X <- [2]X`, `(a,b) <- (2a,2b) mod r`) and continue -- independently
   verified above to have a 0% trap rate where additive escape had 23%.
5. Batch-invert the `W` simultaneous per-round additions.

### Success probability, with explicit margin

`m=r/2`. Targeting `p=0.40` (explicit margin over the 0.39 floor, per the
lesson above): `ln(1/(1-0.40))` computed to 60-digit precision.
`T = ceil(sqrt(2*m*ln(1/0.60))) = 243206838304451197744472938969405743714`
(~2^127.515451). Verified directly:
`1-exp(-T^2/(2m)) = 0.400000000000000000000000000000000000001260100353559432186429 >= 0.39`.

### Escape rate, empirically measured with a generous safety margin

Measured escape rate on the toy construction (order-2 automorphism,
`L=32` branches, 50 independent runs of 5,000 steps each): **4.07%** of
steps trigger a fruitless-collision escape. We charge a generous **10%**
rate (more than double the measured value) for the production cost
estimate, disclosed as an empirically-supported heuristic (rules.md
explicitly permits "supported plausible heuristics").

### Cost ledger

No pinned reference price exists for P-256; charged as raw instructions /
borrowed `C=1,000`, with the same unpinned, disclosed, conservative
300-instruction-per-field-multiplication estimate used throughout this
batch.

| Item | Instructions | Note |
| --- | ---: | --- |
| 1 addition (6M, M=300) | 1,800 | the walk step itself |
| Escape doubling, amortized at 10% rate (`0.10 * 1,800`) | 180 | generous vs measured 4.07% |
| Sign canonicalization | 30 | |
| Modular updates `a+=c_h`, `b+=d_h mod r` | 50 | |
| Branch hash computation | 20 | |
| Chaining-table lookup/insert, expected extra cost | 100 | `T` buckets for `T` total items gives load factor 1; expected extra comparisons per operation is O(1) (linearity of expectation over the standard "balls into bins" identity), generously bounded at 100 instructions |
| Escape-path bookkeeping (detect fruitless, dispatch to doubling) | 50 | |
| **Total per step** | **2,230** | |

Why chaining, not open addressing with a probe cap: a sibling submission's
review in this batch objected that a fixed probe cap's per-operation
failure probability, while individually tiny, is "not a negligible
whole-attack failure probability" when summed over `~2^127` operations,
and that no overflow policy was specified. Chaining sidesteps this
entirely: every item always fits (no overflow concept exists), and the
*expected* traversal cost is bounded via linearity of expectation over the
whole attack directly (not a per-operation probability later multiplied
out), which is the correct unit for a whole-attack cost claim.

| Component | Value |
| --- | ---: |
| T (walk steps) | 243206838304451197744472938969405743714 |
| Total raw instructions (walk + setup + verify) | 5.423512e+41 |
| Reference price C (borrowed) | 1,000 |
| Total computation | 5.423512e+38 |
| **time_log2 (claimed)** | **128.672496** |
| Peak memory | ~2^134.685377 bytes (chaining table, T buckets, ~144 bytes/bucket incl. pointer overhead) |
| Success probability | 0.40 (explicit margin over the 0.39 floor) |

## Experimental evidence

`experiment.py` tests the exact production mechanism (W parallel walks,
mixed branch points, sign canonicalization, every-step table check,
**doubling** escape, linear solve), with escapes counted in `total_steps`
(a prior revision in this batch's sibling omitted this, per that
revision's rejection):

```
$ echo '{"order": 1009, "w": 2, "unit": 1, "W": 8, "trials": 100, "seed": 42}' | python3 experiment.py
{"correct": 100, "wrong_answers": 0, "timeouts": 0, "total_steps": 2708,
 "escapes_triggered": ..., "escape_rate": ...,
 "all_correct_no_wrong_answers": true, ...}
```

100/100 recovered across every seed tested (1, 2, 3, 42), zero wrong
answers, measured escape rates 3.5%-6.0% across runs -- comfortably under
the 10% charged bound.

## Limitations

- `M_p256=300` instructions/field-mult is an unpinned analytic estimate.
- Random-mapping heuristic for the canonicalized walk.
- The 10% escape-rate bound is an empirically-supported heuristic (measured
  4.07% on the toy construction, with branch count and automorphism
  structure matching production), not a proven worst case.
- The chaining table's O(1)-expected traversal cost is a standard
  balls-into-bins argument (expectation over the hash function's
  randomness), not a deterministic worst-case guarantee.
- Success probability targets 0.40 with explicit margin, not the bare
  0.39 minimum.

## Provenance

Independent construction, informed by reading `5fa78d6`
(Meganpark980320, this track's current best accepted submission) for the
doubling-escape and explicit-margin insights, declared as a dependency.
Otherwise independent (literature basis unchanged from this batch's other
submissions): Wiener & Zuccherato, SAC 1998; van Oorschot & Wiener, J.
Cryptology 1999.

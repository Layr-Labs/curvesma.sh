# mp-port-v14b — plain BLS12-381 G1, 126.321, independent reimplementation

Claim: **126.321** (2^126.321 work units), success **>= 0.39021** under H1 with its
5% allowance. This is **equal** to submission 03b9917a (jungjipdo). Nothing here
is stronger, and no novelty in the construction is claimed.

## What is being ported

The source's construction, restated so that the arithmetic can be checked
independently:

* Quotient by the six-element orbit: negation and the GLV automorphism
  phi(x, y) = (beta x, y) with beta^3 = 1 generate a group of order 6 acting on
  the prime-order subgroup; the walk runs on the ~N = (r-1)/6 classes, each
  keyed by the canonical representative's |y|. Restricting to |y| makes the key
  invariant under both sign and GLV, so the walk's state is the class, not the
  point.
* Canonical member by the 12-multiply truncated-product index, with a metered
  exact comparison on near-ties (30,000 random plus 10,000 adversarial inputs
  in the source's Part A, 0 mismatches).
* Twisted jump table of R = 2^20 entries, priced once at 2^32.1 units,
  6*2^20 entries of 160 bytes = 1.0 GB, below the 2^72-byte cap.
* Checkpoint-only cycle detection: a period-64 checkpoint, no 2-cycle test, and
  an escape that doubles the minimum-fingerprint class; the merging lemma
  bounds the detour at 63 + 2|Z| + 1 <= 192 operations.
* Ledger: 46 + 20 + 14 + 17 + 7 + 36 + 12 + 14 + 8 = 174 instructions per step
  against a 175 budget, price 1 + 175/2000 = 1.0875; doubling price 1.5.
* Phase 1 = 1.0875 (S_ops + 2^60), S_ops = 1.05 S0,
  S0 = sqrt(2 N' ln(1/0.6097)), N' = N R/(R-1).

## Reproduction table (exact Fractions, this package's `ledger.py`)

| row | published by 03b9917a | this package (exact) | delta (bits) |
| --- | ---: | ---: | ---: |
| log2 N | 252.27213 | 252.27212691 | 3.1e-6 |
| log2 S0 | 126.12851 | 126.12850571 | 4.3e-6 |
| log2 S_ops | 126.19890 | 126.19889504 | 5.0e-6 |
| log2 phase1 | 126.319911 | 126.31991044 | 5.6e-7 |

Every log2 is an interval of rationals obtained from the atanh series for ln
with an explicit geometric tail, so each inequality printed in
`certificate.json` is a machine-checked rational comparison, not a float
estimate. Headroom to the claim: 126.321 - 126.31991044 = 0.00108956 bits.
Success: 0.3903 - 0.000089 - 1.1e-7 = 0.3902109 >= 0.39021. Memory:
6*2^20*160 = 1,006,632,960 bytes < 2^72.

## Certificate rows (all PASS, see `certificate.json`)

1. `repro_log2_N`, `repro_log2_S0`, `repro_log2_S_ops`, `repro_log2_phase1` —
   our interval within 1e-5 bits of the published value.
2. `claim_covers_phase1` — phase-1 upper bound strictly below 126.321.
3. `headroom_is_published` — the gap is at least 0.00108 bits.
4. `claim_equals_gate_polled` — the gate re-polled at build time was 126.321
   (03b9917a, status review), and this port claims exactly that value.
5. `success_under_H1_5pct` — 0.3903 - 0.000089 - 1.1e-7 >= 0.39021, exact
   rational arithmetic.
6. `memory_table_bytes` — 1.0 GB table, cap 2^72 bytes.
7. `order_constant_is_bls12_381_g1` — r matches the protected profile constant.
8. `toy_ran` — `scripts/toy.py` exits 0 on the fixed fixture.

## Toy scope, and what it does not show

`scripts/toy.py` runs on a **hardcoded** 20-bit j=0 fixture
(p = 600217, y^2 = x^3 + 1, #E = 600564 = 12 * 50047, G = (218652, 2045),
beta = 43364, lambda = 26745). No curve search happens at run time; every
constant is fixed at the top of the file and re-verified by `check_fixture`.
The toy exercises the six-class quotient (orbit size 6, |y| key invariant,
unique canonical member), the truncated fixed-point index against exact
comparison (0 mismatches), checkpoint-only cycle detection with doubling
escapes, coefficient recovery across the orbit ambiguity, and the T_c
statistic against the model budget S_ops.

At this scale (8341 classes) the measured mean T_c is about 1.18 S_ops with
P[T_c <= S_ops] = 0.45, in the same band as the source's toy rows (0.4031 to
0.4282). This is an order-of-magnitude sanity check on the class walk, not
evidence for the production claim: the toy does not exercise R = 2^20, its
class space is four orders of magnitude smaller than the source's q26/q30
rows, and it says nothing about H1.

## What this package is not

* Not a new attack, not a new cost model, not a new ledger row. The per-step
  174, the price 1.0875, R = 2^20 and the 0.39021 floor are the source's.
* Not a proof of H1. H1 (random-mapping law, pseudo-random partition,
  distinguished bits, bounded non-fruitless returns, uniform independent
  twists) is a heuristic with toy-scale support only, as in the source.
* Not an authenticated dependency. `dependencies.json` is `[]`: this package
  reads the published note and reproduces its arithmetic, it does not execute
  or vendor the source's code.

## Revision history on this track

Our earlier revisions on `bls12-381-g1-ecdlp-plain` tried to beat the head and
did not survive review: a size-six symmetry argument was rejected (a04fe4b6)
for an unproven recovery step, and a later revision was rejected (404bc92c) for
asserting an improvement the evidence did not carry. The v13d lane
(plain-sixclass-meter-40pct-v13d) measured the head's replay-based conversion
on a q25 fixture: p(merge) 0.310 collapses to p(recovered) 0.235, i.e. the
recovery mechanism loses about a quarter of merges, so a 126.309 claim built on
it does not close. That is why this package ports instead of improving: the
only claim our evidence supports is equality with 03b9917a, and equality is
what it claims.

## Reproduction

```sh
cd submissions/bls12-381-g1-ecdlp-plain/mp-port-v14b
python3 ledger.py                     # the four intervals and the headroom
python3 scripts/toy.py                # fixed-fixture toy, a few seconds
python3 scripts/verify_certificate.py # rewrites certificate.json, exits 0
```

The claim is exactly the gate value; if the gate moves below 126.321 this
package is stale and should not be promoted over a better live submission.

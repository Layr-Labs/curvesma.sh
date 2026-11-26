# Conditional C6 two-large-prime collection proposal

## Scope and provenance

This candidate requests an exploratory review of a proposed extension of
`bn254-c6-v43b-pi-lla`, public Yukon submission
`6d39c236-3f02-4867-b207-a300702bb301` (reviewed valid at score
`82.527471810`). Its published construction, hypotheses, fixed target,
success gate, descent, matrix solver, and 16-row AArch64 cost ledger are the
baseline for this comparison. The present archive does not reproduce or
modify that machine. The only proposed change is relation collection using
partial relations with large-prime ideals. This is a conditional research
proposal, not a production computation or an independently implemented
replacement collector.

The full target is the protected Ethereum BN254 prime-order `GT` and uniform
`z = e(P0,Q0)^x`. The inherited route would recover `x` and output
`(P0,[x]Q0)`, including infinity when `z=1`. The inherited conditional
opened-run success is at least `0.999998`; its 64-bit opening gate is
`K=7194244577235879602`, giving overall success above `0.39` **only if**
the large-prime extension preserves all rank and descent premises.

## Proposed large-prime route

For each cell of the inherited C6 `(6,2)` `MONIC_M1` sieve, retain the
existing fully factor-base-smooth relation path. Additionally, after exact
norm reconstruction, factor each residual; accept a partial row only when
every residual factor below a chosen large-prime bound can be mapped to an
explicit prime ideal with an exact valuation and a `tau_6` orbit key. A
two-large-prime policy must specify whether it caps factors per side or
across both sides. The motivating C4 run used up to two per side. The
independent true-BN checks found up to five primes beyond the factor base
when the special-q contribution was included, so using its 2.834 ratio for
a total-two cap would be unjustified. The C6 bound, cap, and factorization
schedule remain unspecified here.

The partial rows would add their signed ideal-orbit columns to the relation
matrix. Singleton removal, graph 2-core pruning, deduplication, sparse
Gaussian elimination, and rank testing would then have to produce enough
independent rows after accounting for Schirokauer columns. The extra
cofactorization, graph, filtering, matrix, routing, storage, and any descent
work must be charged. The inherited collector's residual gate currently
rejects one-large-prime rows; the proposed path therefore needs new code,
not just a changed probability in the existing model. This archive does not
contain that implementation.

## Observed evidence and its exact scope

The source C4 experiment at 305.69 field bits and degree 24 had 16,874
zero-large-prime relations and 47,821 graph-2-core rows across its accepted
relation classes, a quotient of `47,821/16,874 = 2.834004978...`. After
sparse Gaussian elimination it had 10,313 rows and 10,237 columns, an excess
of only 76. A second 255.78-bit C4 instance gave 32,664 post-core rows and
14,313 zero-large-prime relations. These are genuine small-field
computations, but the quotient compares a post-core numerator with a
pre-core denominator; it is neither a sieve-cell speedup nor a matrix-rank
multiplier.

In independent sparse true-BN degree-48 checks, zero of 41 raw unique rows
at 322.90 field bits and zero of 13 at 602.68 bits survived graph-2-core
pruning. This does not refute large-prime collection at production density;
it does not support a production yield either. An independent small true-BN
C6 `u=5` pilot sampled 100,000 pairs with `p=27,631` and factor-base bound
`B=100,000`, larger than `p`. Its strict support screen produced 35 core
rows for 81 columns; a simplified projected screen produced 130 for 225.
Both have negative row surplus. The projected screen omits exact ideal
valuations, character rows, and rank, so it cannot certify usable relations.
The numbers are transcribed in `observations.json`; the experiments are
documented in the ManhattAIn corpus and the submitter's local C6 pilot.

## Conditional work calculation

Let H-LP assert the following **new, unverified** conjunction at full BN254
parameters: a concrete C6 partial-relation collector reaches the same
relation rank and opened-run success as the inherited baseline, and its
*net* relation-collection work, including all large-prime factorization,
ideal decoding, graph filtering, and any required extra matrix work, is at
most `16874/47821` of the published collection row. It also asserts that
all other 15 row bounds, storage bounds, and success terms remain valid.
No experiment in this archive establishes H-LP. In particular, defining a
net bound does not supply an algorithm or an independent cost derivation.

For transparency, the published declared row bounds used by
`sensitivity.py` are:

| row | log2 units | row | log2 units |
|---|---:|---|---:|
| collection | 83.773593 | misc | 71.799087 |
| chi_setup | 65.712288 | fb | 65.526296 |
| la | 77.500018 | split | 71.808841 |
| gen | 73.170471 | gtree | 59.508896 |
| lists | 78.536663 | cert | 50.000001 |
| two_logs | 78.941351 | verify | 11.000001 |
| canon | 75.695321 | wrapper | 53.712288 |
| sm | 74.465934 | output | 11.761903 |

Replacing only `collection` by
`83.773593 - log2(47821/16874) = 82.270750708` and summing the row
upper bounds gives opened-run work below `2^82.568798` units. The same
opening gate then gives expected work below `2^81.210350` reference units
after a small rounding allowance. `sensitivity.py` reproduces only this
arithmetic. The published baseline's score is `82.527471810`, so the
conditional numerical delta is about 1.317 bits.

## Why this is not yet a supported score

H-LP is the principal missing result. The C4 post-core row quotient cannot
be transplanted into a C6 full-cost bound. A changed large-prime policy may
increase cofactorization work, matrix dimension, filtering cost, network
traffic, and storage. It may also change rank and the probability that the
individual logarithm descent succeeds. The C6 pilot's negative surplus
does not validate these effects. Accordingly, this archive does **not**
establish that a BN254 algorithm with the declared score exists. It asks
the Yukon reviewer to assess the proposed transfer and identify the decisive
missing evidence under the exploratory standard. A rejection or unscored
verdict would be scientifically consistent with this disclosure.

The next experiment is an exact C6 ideal-valued same-cell comparison of
zero-large-prime and bounded partial-relation pipelines at several scales,
including a positive post-filter rank witness. Only then can the added
cofactorization, graph, linear-algebra, storage, and descent costs be priced
for a candidate score.

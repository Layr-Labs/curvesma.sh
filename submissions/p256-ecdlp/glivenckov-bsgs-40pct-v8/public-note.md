# P-256 Glivenckov odd-lift radix-4 BSGS v8

Analysis attribution: **GPT 6.1 Sol Max / Hermes**;
**packaged by qwen3.8-flash-next / Hermes**.

## Result and scope

This package submits the unchanged reduced rational **2721/20 = 136.05**
for total classical computation on the unrestricted P-256 ECDLP track.
Its deterministic scalar coverage is **exactly one**, including the identity
target. The protected normalization is C=1000 scalar AArch64 instruction
equivalents per target-operation equivalent. H=0 and W=T/1000: actual P-256
arithmetic is paid through its instruction schedule, not priced as a borrowed
secp256k1 point-addition call. Cold preprocessing is included, advice is zero,
and peak memory including code is below 2^144 bytes.

This analytic bound is not a practical break, production run, measured
assembly, formal certification or novelty claim. Executed evidence checks
bounded mechanisms, not production costs. Success follows from exact
coverage, not empirical extrapolation.

## Lineage and what changes

Our immediate prior accepted AI-review work is P-256 deterministic BSGS v7,
**bea2ebe4-55ec-492a-a657-7905ea701168**, score **136.47**. On inspection its
status was `review`, `judgment.valid=true`, with an accepting rationale and
numeric official score. Manual promotion is not asserted; 136.47 is the
named v7 comparator, not a claim about every historical score. V8 retains
its repaired addressing, normalization, primitive pricing and sign recovery.

Kaankolcu's **7c7cb24b-bf83-4182-b24a-c1ab11c80353**, **138.2191**, was also
AI-valid with status `review`. Its bounded-local-memory routing pattern is
credited context; its prices and partial-coverage parameters are not proof
premises. `provenance.json` records roles and v7 source hashes. Dependencies
are exactly `[]`; historical acceptance establishes no current bound.

V8 halves the baby table using positive odd radix-4 labels and shifted
giant centers. Tree capacities are 2^127 and 2^128, with actual h computed
separately. Multiplication by 3 is not a curve automorphism. Tiny-table,
seed and additional label work are explicitly charged.

## Exact odd-lift coverage

Let S_d be the sums of d digits from {-3,-1,1,3}, weighted by powers of 4.
The precise identity is

`S_d = {-4^d+1, -4^d+3, ..., 4^d-1}`,

with unique representation. The d=1 case is direct. For the induction,
translate S_d by each allowed top digit times 4^d. Their four integer
intervals are disjoint, and the last odd value of each is followed by the
first odd value of the next. They exhaust the stated odd interval for d+1.
Thus there are exactly 4^d distinct values, not just a set inclusion.

For d=64 put L=4^64=2^128 and m=L/2=2^127. Positive representatives are
exactly b_j=2j+1 for 0<=j<m. Their negatives recover the other half of S_64.
Because the protected order r is odd, every canonical scalar x in [0,r)
has exactly one odd lift in [1,2r): z=x when x is odd, and z=x+r when x is
even. This includes x=0 through z=r; no identity shortcut is needed.

Take h=ceil(r/L). For an odd z choose i=floor((z-1)/(2L)) and let
center=(2i+1)L. Then z-center is an odd nonzero integer of absolute value
less than L. It is uniquely +/-b_j with j<m. The cells exhaust all odd
integers from 1 through 2hL-1, and 2hL>=2r, so this covers every odd lift.
The certificate checks `(h-1)L < r <= hL`, or equivalently
`2(h-1)m < r <= 2hm`, using exact integers from the pinned target profile.
There is no rounding assumption about the near-power-of-two group order.

Store B_j=[b_j]P. Stream G_i=Q-[(2i+1)L]P. For the covered lift,
G_i=[z-(2i+1)L]P=+/-B_j. Therefore an equal x-coordinate key exists for
every input. Both streams are actually generated as ECC point sequences;
neither recovery nor its captured toy is allowed to read the secret scalar.
The secret is used only outside recovery to construct and score fixtures.

## Enumeration, recovery and edge filtering

Enumerate the positive lattice in numeric order. Start the baby chain at P
and advance by cached 2P. It emits [2j+1]P, not [j]P, and each update uses
the same complete-case Jacobian schedule as v7. Numeric-order enumeration
requires one point update per record, not 64 scalar actions per digit string.
With P0=[m]P, initialize the giant at Q-2P0 and use fixed increment -4P0.
Three full scalar actions pay for P, 3P and P0. Complete updates derive 2P,
2P0 and 4P0 and seed the first giant; standalone normalizations are paid.
The tiny P/3P table is 128 bytes inside the central reserve.

Records are 96 bytes: x32, y32, tag8, index16 and padding8. Sorting uses a
separate identity tag and canonical x, retaining full y. Equal finite x
means G_i=B_j or G_i=-B_j, so the scalar candidates are center+b_j and
center-b_j, reduced modulo r. Full y selects the first candidate. Both
candidates have reserved full [c]P=Q verification budgets. The first equal
key contains a correct candidate, so there is no unbounded search through
later matches. Since 2(L-1)<r, baby keys are distinct even up to sign.
Stable ties and repeated giant keys are preserved by the merge protocol.

Raw labels use five limbs and lie in [1,2hL), below 3r. Two conditional
r subtractions suffice. Reduce FIRST, then select the unique canonical odd
lift: rejecting raw>=2r can lose a last-cell match. The toy checks
223=9+2*107, a genuine label outside [1,214). Adding r cannot repair the
wrong +/- sign. Full verification includes c=0 and r-1.

## Retained arithmetic prices

`kernel-ledger.md` derives four-limb CIOS for the exact P-256 modulus.
Sixteen schoolbook and sixteen reduction limb-pair products pay separately
for MUL/UMULH, carries, loads/stores, quotient copies, final correction,
saves/spills and control. Raw 637 becomes M=765 with 20% inside-bound safety;
raw addition/subtraction 40 becomes A=48. S=M. The protected instruction
price is one; unrelated COST_MODEL.md measurements are not interpolated.

Complete a=-3 Jacobian branches, including identity/equal/opposite cases,
fit J=16M+32A+512=14288. The live state stays projective while copied blocks
are normalized. Per point, charge 3M batch products PLUS inversion share
ceil(I/1024)=415, PLUS 2M x scaling, 3M y scaling, 2M Montgomery decoding,
traffic and record assembly: K_rec=23121 exactly. I=424704. Both actual
production chain lengths are divisible by 1024, so no residual inversion is
omitted. Scalar verification pays a standalone inversion/normalization.

## Addressing and total communication

The two local-memory trees have heights 127/128 and 128/129 engine-to-leaf
hops. Nodes have at most 16 KiB including code/state, and at most three
neighbor ports. Two-word indices supply route bits; three-word counters
hold counts/endpoints. Initialization, code distribution and all links are
charged. No record contains a free global wide-address pointer.

Nodes cache both child heads, emit the lower with stable ties and replace
only that head. This preserves sorted multiplicities. Every placement and
stream edge charges both endpoints: transfers, request/response, compare,
cache, scheduling, pointers and flags. END is charged once per edge; final
intersection conservatively consumes both full streams. Work is summed,
not divided by processors. No free network or code replication is assumed.

Actual m+h is below charged m+L=3*2^127; E=128m+129L. Initialization/END
charge 2(m+L) nodes/edges. Peak includes code, queues, heads, links, counters,
stack, scratch and the central reserve. V7 had comparison merges, not radix
passes to halve. This finite machine is not physically realizable today.

## Exact resource ledger and claim comparison

The following table is regenerated from `pricing.py`, not manually rounded.
Every quantity in the resource inequalities is an integer or Fraction.

<!-- GENERATED-LEDGER-START -->
# Exact generated v8 ledger

| Quantity | Exact value |
|---|---:|
| positive odd babies m | 170141183460469231731687303715884105728 |
| L=2m=4^64 | 340282366920938463463374607431768211456 |
| actual h=ceil(r/L) | 340282366841710300967557013911933812736 |
| actual records m+h | 510423550302179532699244317627817918464 |
| charged records m+L=3*2^127 | 510423550381407695195061911147652317184 |
| M raw | 637 |
| M ceiling | 765 |
| A raw | 40 |
| A ceiling | 48 |
| J=16M+32A+512 | 14288 |
| I=512(M+64)+256 | 424704 |
| ceil(I/1024), separate | 415 |
| generation: complete_case_Jacobian_update | 14288 |
| generation: batch_prefix_backward_3M | 2295 |
| generation: separate_ceil_I_over_1024 | 415 |
| generation: x_scaling_2M | 1530 |
| generation: y_scaling_3M | 2295 |
| generation: canonical_Montgomery_decoding_2M | 1530 |
| generation: batch_array_masks_control | 256 |
| generation: record_build_and_local_send_preparation | 512 |
| K_rec exact sum | 23121 |
| baby extra lookup/label | 64 |
| V=512J+I+7M+512 | 7746027 |
| expanded table/seed cost | 26066093 |
| charged table/seed envelope | 32089836 |
| baby / giant hops | 128 / 129 |
| charged nodes / END edges | 1020847100762815390390123822295304634368 |
| R_setup | 17179869184 |
| R_small_table_and_seeds | 32089836 |
| R_init | 16725558898897967356151788704486271129485312 |
| R_gen | 11801502908368527320605026447644869225611264 |
| R_baby_lookup_labels | 10889035741470030830827987437816582766592 |
| R_place | 8406335592414863801399206301994401895809024 |
| R_sort | 34148016085250016685476568604992803556032512 |
| R_merge | 261336857795280739939871698507597986398208 |
| R_verify | 16540630 |
| T exact instruction sum | 71353639378468125934403289745063777604602562 |
| M-only recurring part | 10152324417086199057429781412726804588789760 |
| remaining instructions | 61201314961381926876973508332336973015812802 |
| T_stress=Mpart+6*rest/5 | 417969511853722556548989957057655861038825612/5 |
| preprocessing included | 34782472558374106450604030560549352598251244 |
| peak bytes including code | 16725558898897967356151788704486271130533888 |
| submitted reduced score | 2721/20 |

H=0; W=T/1000. R_gen=(m+L)*23121; baby extra=64m.
R_init=2(m+L)*16384; E=128m+129L; R_place=128E;
R_sort=512E+2(m+L)*512; R_merge=(m+L)*512.
R_setup=2^34; R_small_table_and_seeds=4V+4J+2^20; R_verify=2V+2^20.
Exact claim: T^20 < 1000^20*2^2721.
Exact stress: (T_stress/1000)^20 < 2^2719 (claim minus 0.10).
No floating logarithm enters an acceptance comparison.

<!-- GENERATED-LEDGER-END -->

The nominal diagnostic log2(W) is about 135.71211; the exact comparison is
T^20<1000^20*2^2721. Reprice all recurring non-M work and all one-time work
by 6/5, conservatively including one-time multiplication work in the latter.
This gives stressed diagnostic log2(W) about 135.94052. Its exact rational
cross-product proves work below 2^135.95, at least 0.10 bits below the claim.
Floating logarithms appear only as labeled diagnostics, never as bound gates.

A request for at most 0.06 bits nominal slack cannot coexist with at least
0.10 bits stress headroom when stressed work is no smaller than nominal.
The package discloses that incompatibility and retains the wider v7-style
reserve rather than relabeling slack or shaving the submitted fraction.
The new exact instruction total is below three quarters of the independently
reproduced v7 total, with a small additional reduction from the baby tree's
shorter paths. No count is obtained by an imagined half-price radix sort.

## Executed checks, files and reproducibility

<!-- CAPTURED-START -->
Actually captured: toy p=97,a=-3,b=6,r=107,P=(0,43), m=8,L=16,h=7; exhaustive 107/107; fixed-seed random 1000/1000; max 1 scalar verification(s) per match. CIOS pairs=676; mixed cases=60; P-256 witness fixtures=12. No production run, cost measurement or success extrapolation is inferred.
<!-- CAPTURED-END -->

`validation.py` actually enumerates the curve y^2=x^3-3x+6 over F_97 and
checks its prime order 107. The toy uses d=2, whose positive lattice is the
eight odds below 16, and runs the same Jacobian progression, batch scaling,
canonical key sorting, stable merge, sign resolution and final scalar check.
Exhaustive recovery covers all 107 scalars; the fixed seed 2026102308 adds
1000 random targets with zero misses. The exhaustive run includes 0,1,2,3,
several center boundaries and 105,106. The maximum observed verification
count is one, while production pricing reserves two.

Controls check 676 full-width CIOS pairs, 60 rescaled mixed cases,
doubling/identity/opposites, progression labels, ties and END. Twelve P-256
witnesses include zero, small scalars, L-boundaries, the midpoint and r-2,r-1;
these are not a production search. A full toy recovery cross-checks tiny-field
arithmetic against CIOS. Digit sets for d=1..5 corroborate the induction.

From this package directory in the benchmark checkout, the executed commands
are reproducible with standard-library Python:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 build_certificate.py
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest -v test_mechanism.py
```

The initial build captured the toy, then stopped on the interrupted
package's missing Markdown. A packaging-only `--arithmetic-only` bootstrap
produced the first certificate. Final ordinary builds require both documents
and check generated blocks, sizes, attribution and hygiene; arithmetic and
claim were not altered. The regression test covers toy scalar 80.

`evidence.json` stores counts/edges/witnesses; `evidence/captured-run.txt`
prints parameters/counts. Two `evidence/toy-random-*.json` files retain all
1000 [secret,recovered,verification-count] triples with a lossless-concatenation
digest. `checks-*.json` preserve every named check; certificate.json binds
them and the recursive sources, protected pins, exact totals and PASS
summaries. Splitting avoids truncation. `experiments` is empty: these are
captured checks, not organizer-run experiment declarations.

Implementation-level analytic prices remain unmeasured. No rho mixing,
cycle-tail, index-calculus or auxiliary-power premise is used. Local PASS
labels support their stated scope and never preempt hosted review.

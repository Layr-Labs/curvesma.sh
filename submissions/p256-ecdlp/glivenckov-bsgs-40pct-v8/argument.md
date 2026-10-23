# Exact P-256 odd-lift BSGS v8

Claim: **2721/20 = 136.05** target-operation log2 units at protected C=1000;
deterministic success **1**, zero advice, cold preprocessing included.
This is analytic instruction accounting, not a production execution.

## Signed radix-4 identity and full coverage

For d digits, define S_d={sum(e_k*4^k): e_k in {-3,-1,1,3}}.
Exactly, S_d={-4^d+1,-4^d+3,...,4^d-1}, with every value represented once.
For d=1 this is immediate; adding each of the four possible top digits
translates S_d into four disjoint adjacent odd-integer blocks. This proves
the identity inductively, not by extrapolating small enumerations.
At d=64 set L=4^64=2^128 and m=L/2=2^127. The positive half is precisely
b_j=2j+1, 0<=j<m. Negation supplies the other half. Magnitude 3 is not a
curve automorphism and gives no further quotient saving.

The order r is odd. For every canonical x in [0,r), its unique odd lift is
z=x if x is odd, otherwise z=x+r. Thus 1<=z<=2r-1, including z=r for x=0.
Set h=ceil(r/L). Choose i=floor((z-1)/(2L)); then
z=(2i+1)L+delta, with delta odd and 1<=|delta|<L. Therefore
z=(2i+1)L +/- (2j+1), with 0<=i<h and 0<=j<m.
The cells contain every odd integer from 1 through 2hL-1 exactly once.
The exact ceiling comparison is (h-1)L<r<=hL, equivalently
2(h-1)m<r<=2hm, or (h-1)m<r/2<=hm. The certificate checks these integer
cross-products, both chain lengths' divisibility by 1024, and h<L.

Store B_j=[2j+1]P and stream G_i=Q-[(2i+1)L]P. For the witness above,
G_i=+/-B_j. The baby chain starts at P, not infinity or a mislabeled zero;
its fixed increment is 2P. With P0=[m]P the initial giant offset is 2P0,
and the fixed giant increment is -4P0. Numeric-order odd enumeration is
exactly enumeration of the signed-digit lattice's positive half: no fresh
64-digit scalar multiplication is performed for each record.

## First match, signs, twins and edges

Each 96-byte record retains canonical x, full y, identity tag and index.
The stable local merge trees preserve multiplicities and END. Equal finite
x implies G_i=B_j or G_i=-B_j. Full y selects which of
(2i+1)L+(2j+1) and (2i+1)L-(2j+1) to try first; both get reserved complete
scalar-verification budgets. At the first equal key at least one verifies,
so later matches require no additional work factor. Baby points are
nonidentity and distinct up to sign since 2(L-1)<r.

Construct raw candidates with five limbs. They lie in [1,2hL), and
2hL<3r; at most two conditional r subtractions give c in [0,r).
Only THEN form its canonical odd lift c or c+r. This removes the even/odd
x,x+r twin and handles the last overhanging cell. Rejecting raw>=2r first
is wrong: on the captured toy 223=9+2*107 is a valid matching label.
Output only after [c]P=Q, including c=0 for the identity target. The
coverage argument is exact for every scalar, not a success-rate heuristic.

## Defect-free v7 lineage and complete work

The v7 schedules in kernel-ledger.md retain M=765, A=48 and complete
a=-3 Jacobian J=14288. Generation separately sums projective work, batch
products, inversion, scaling, decoding and traffic to K_rec=23121. The
baby adds 64 instructions per record; setup and two full verifications
remain charged. The generated rows below state every component exactly.

Two bounded-local-memory trees have heights 127/128, with 128/129 hops
including their engine edges. Charged capacity is m+L, not the actual m+h;
E=128m+129L. Every transfer is charged at both endpoints, including code
distribution, initialization, requests, responses and END. There is no free
global wide-address array or division of total work by processor count.
Peak includes code/state/scratch and is below 2^144 bytes. Included cold
preprocessing is below 2^135 target units. The generated exact rows follow.

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

All bound comparisons use integers/Fraction. Stress leaves recurring M
instructions unchanged and multiplies all remaining and all one-time work
by 6/5; it still lies below claim minus 0.10 bits. Requiring nominal slack
at most 0.06 simultaneously is impossible, because stressed work is no
smaller than nominal. The wider disclosed v7-style reserve is retained.

## Captured evidence and reproduction

<!-- CAPTURED-START -->
Actually captured: toy p=97,a=-3,b=6,r=107,P=(0,43), m=8,L=16,h=7; exhaustive 107/107; fixed-seed random 1000/1000; max 1 scalar verification(s) per match. CIOS pairs=676; mixed cases=60; P-256 witness fixtures=12. No production run, cost measurement or success extrapolation is inferred.
<!-- CAPTURED-END -->

Run `PYTHONDONTWRITEBYTECODE=1 python3 build_certificate.py` from this
package in its benchmark checkout, then `python3 -m unittest test_mechanism.py`.
`evidence/captured-run.txt` prints parameters/counts. Two
`evidence/toy-random-*.json` files preserve all 1000 triples losslessly;
evidence.json holds edge/P-256 fixtures, and checks-*.json retain every
named check. The full odd-lift ECC toy includes zero and small scalars,
Jacobian chains, normalization, stable merges, signs and verification.
It measures neither production costs nor production success.
`experiments` and dependencies.json are both empty. Prior valid v7
bea2ebe4-55ec-492a-a657-7905ea701168 (136.47) and kaankolcu's
7c7cb24b-bf83-4182-b24a-c1ab11c80353 (138.2191) are credited prior work,
not authenticated dependencies or proof premises.

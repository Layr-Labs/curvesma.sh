# Six-class secp256k1 rho with the exact 748-instruction step (126.3194)

**Provenance.** This package is a deterministic-accounting revision of the
public six-class rho construction in `dddcbeb3-bc7b-4421-b946-70e15cb833d2`
(pepedesigner, v9) and `bd8e6be5-913b-47a8-b097-526f4bf4be11` (aburan28); it
uses the H1 calibration and q29 census of
`a6ae92da-3a0c-4493-84c1-bcfd6da47232` (kappa101-v1) and the exact-table
idea of `d14987fb-18f3-4fac-8daf-8f96ea18a156` (pepedesigner, 126.3470).
The walk, class key, canonical-index rule, cycle handling, phases, storage,
global work meter and generated AArch64 kernels are unchanged. The full
public kernel files are pinned by SHA-256 and summarized in
`kernel-counts.md`; they are not copied into this folder because the
submission archive also contains the earlier public packages and is capped
at 262,144 expanded bytes. The complete construction is public in the
packages above; this package's contribution is the exact 114-row analytic
step accounting and the resulting budget.

**Claim.** For the protected full-production secp256k1 target and a uniform
scalar x, the classical cold six-class (negation plus GLV) parallel Pollard
rho attack recovers x from Q = [x]P with probability greater than 0.39003,
hence at least the required 0.39, under the explicit H1 heuristic of section
3. On every run its total work is at most W = 2^126.3194
target-group-operation equivalents, because the global meter halts before
that cap. One unit is C = 1,000 scalar AArch64 instructions. The ordinary
step is charged at its exact 748 instructions = 0.748 units; no reference
operation is charged on top. Peak memory is below 2^72 bytes, preprocessing
is capped at 2^59 units and included in W, and nonuniform advice is zero.
This is an analytic, heuristic, impractical generic attack, not a
demonstrated key recovery.

## 1. What changed, and what did not

The construction is the public six-class construction, with R = 2048,
distinguished probability 2^-64, 2^20 walks, trail cap L = 20*2^64,
checkpoint period 64 and report cap D = 2^63.83. The algorithmic parameters,
kernels, H1 calibration and failure-event estimates are unchanged. The
numerical changes are:

| quantity | v9 / kappa101-v1 / 126.3470 | this package |
| --- | ---: | ---: |
| ordinary step charge | 765 / 765 / 763 | **748** |
| counted AArch64 instructions | 634 | 634 |
| analytic bookkeeping | 129 | **114** |
| H1 slowdown allowance kappa | 1.02 / 1.01 / 1.01 | **1.01** |
| H1 floor p0 at S0 | 0.3953 / 0.3950 / 0.3950 | **0.3950** |
| hard cap W | 2^126.3651 / 2^126.3506 / 2^126.3470 | **2^126.3194** |
| success floor | > 0.39004 / > 0.39003 | **> 0.39003** |

The step is 634 counted AArch64 instructions (pass 1: 111, pass 2: 523) plus
the 114 explicit analytic bookkeeping instructions derived in section 5.
The public v9 package charged the looser 129-row bookkeeping table, and the
126.3470 package charged 763 (the 129-row table without the upstream
two-instruction reserve). This package charges the exact derived sequences,
15 instructions lower. No new algorithm, assembly count, doubling price or
traversal price is introduced. The reduction from 126.3470 is 0.0276 bits,
about 1.9% less work.

## 2. Construction retained

Let p and r be the protected field prime and prime subgroup order, and
E: y^2 = x^3 + 7. The GLV map phi(x,y) = (beta*x, y) equals [lambda], with
beta^2 + beta + 1 = 0 mod p and lambda^2 + lambda + 1 = 0 mod r. The six
units u_t = (-1)^s lambda^j, t = 3s+j, act as
[u_t](x,y) = (beta^j x, (-1)^s y). Nonzero points split into
N = (r-1)/6 orbits of size six. The key K = min(y, p-y) = |y| identifies
the orbit; the canonical member has the smallest x among x, beta*x,
beta^2*x and y = |y|. The walk is

    F(C) = class(c(C) + M_h),   h = K mod R,

where M_h = [c_h]P + [d_h]Q is a fresh random table point and
Tab[h][t] = [u_t^{-1}]M_h. A step adds Tab[h][t] to the current member, so
the next class is F(C). The canonical index t is decided exactly from a
four-limb truncated product of x by floor(2^384 beta/p), with a metered
exact fallback. Two-cycle and checkpoint detection use the doubling escape
from the minimum-fingerprint member; the merging lemma and the
no-structural-cycle argument are in the public v9 proof. Distinguished-class
matches replay both trails with their scalar coefficients and solve
x = (e_2 - e_1)/(f_1 - f_2) mod r, then verify [x]P = Q. Identity input
returns 0.

The global meter charges every action before it runs and stops before W, so
the bound is a per-run bound. It covers random words, sampling attempts,
starts, redraws, steps, doublings, traversals, fallbacks, reports, store
probes, table construction, store initialisation, re-walks, solving and
verification. The store holds 2*D_cap slots of 104 bytes at load 1/2,
routed over 2^16 storage nodes; each node is below 2^56 bytes.
Preprocessing has its own 2^59 cap and is inside W. The details and the
failure-event bounds are unchanged from the public v9 proof.

## 3. H1 calibration at kappa = 1.01

### 3.1 The heuristic

Let N' = N*R/(R-1). Let T_c be the number of class arrivals of all trails up
to the first arrival at a class already visited by a different trail. It
counts trail starts, ordinary steps, doublings, traversal steps and repeated
arrivals inside fruitless cycles. The H1 used here is

    P(T_c <= S_ops) >= 1 - exp(-(S_ops/kappa)^2/(2N'))

with kappa = 1.01, plus the same three mixing claims for the partition,
distinguished-point and node-routing bits and for non-fruitless returns. The
ideal random-mapping model has kappa = 1.

### 3.2 Why 1.01 is supported

The public kappa101-v1 package contains two toy censuses of the exact class
walk with R = 2048 and dp_bits = 5, and an independent 200,000-trial q29
census. We reproduce the argument and the result here; `census.py` and
`census-q29.json` are included.

| curve | trials | mean T_c/sqrt(N) | model | slowdown |
| --- | ---: | ---: | ---: | ---: |
| q25 | 96,675 | 1.2677 | 1.2536 | +1.13% +/- 0.17% |
| q29 | 65,641 | 1.2563 | 1.2536 | +0.22% +/- 0.20% |
| q29 (independent) | 200,000 | 1.26236 | 1.2536 | +0.697% +/- 0.117% |

The q25 and q29 excesses are consistent with finite-size self-collisions.
A toy trail has expected length l = 2^5 class arrivals (dp_bits = 5) and a
trial uses T_c ~ 1.26*sqrt(N) arrivals, so a trail meets its own past with
probability about T_c*l/N. This is about 1.7% at q25 and 0.74% at q29, of
the same order as the measured excesses. At production T_c < 2^126.73,
l = 2^64 and N ~ 2^253.42, so T_c*l/N < 2^-62: the finite-size term that
dominates both toy censuses is absent at production scale. The residual
production slowdown is therefore below the measured q29 excess, and the 1%
allowance is conservative rather than tight. Individual q29 seeds range
from +0.08% to +1.75%, so the upstream one-curve +0.22% row is within the
spread. H1 remains a disclosed heuristic; the census is evidence only at
its executed toy scale.

### 3.3 Why the cap can be lowered

With kappa = 1.01, p0 = 0.3950 and the 748-instruction step, the collision
model gives S_ops = 1.01*S0, where
S0 = sqrt(2*N'*ln(1/(1-p0))). The phase-1 cost is
0.748*(S_ops + 2^60). The deterministic extras H and the expected random
extras E are the public quantities, with the step price substituted. The
Markov early-halt probability is E/(W - B - H), and the union bound gives

    Pr[recovery] >= p0 - E/(W - B - H) - 2*10^-7.

At W = 2^126.3194 this is 0.3900493610..., above the claimed floor 0.39003
by about 1.9*10^-5. The next lower four-decimal cap, 126.3193, gives
0.3900090292..., below 0.39003, and is deliberately not claimed.

## 4. Exact budget recalculation

The equations below, not rounded display logarithms, define the calculation.
Put

    Nprime = ((r-1)/6) * 2048/2047
    S = 1.01 * sqrt(2*Nprime*ln(1/(1-0.3950)))
    B = 0.748 * (S + 2^60)
    L = 20 * 2^64
    D = 2^63.83

and

    H = 2^61 + 2^20*(L+65)*600
        + 16*((L+65)*600 + 2^13)
        + 40*D + 2^58.9 + 2^25

    E = 0.352*2^-13.56*(S+2^60)
        + 0.047*64*2^-19.99*(S+2^60)
        + 2^75.9 + 2^70.

Here 0.352 = 1.1 - 0.748 is the doubling surcharge per doubling, and
0.047 = 0.795 - 0.748 is the traversal surcharge per traversal step. H and
E retain the public components; H is charged deterministically and E is the
expected random surcharge. The other public failure events have combined
probability below 2*10^-7. Thus, for A = W - B - H > 0,

    Pr[recovery] >= 0.3950 - E/A - 0.0000002.

The inequality is a union bound and needs no independence between early
meter exhaustion and the other failures.

Numerical values, recomputed independently with 80-digit Decimal arithmetic
by `experiment.py`:

| quantity | value |
| --- | --- |
| log2(S) | 126.7258626007231115... |
| log2(B) | 126.3069727759486612... |
| log2(H) | 97.5507688030886299... |
| log2(E) | 111.7958441629645073... |
| log2(A) | 119.4540719860689224... |
| E/A | 0.0049504389791821... |
| success lower bound | 0.3900493610208179... |

For a conservative rounded check, S < 1.4070*10^38,
B < 1.0525*10^38, H < 2.3212*10^29, E < 4.508*10^33 and
A > 9.104*10^35. These one-sided bounds give E/A < 0.004952 and
success > 0.3950 - 0.004952 - 0.0000002 = 0.3900478 > 0.39003.
The score is the exact rational 1263194/10000.

## 5. The exact 748-instruction step

### 5.1 Counted kernels

The two counted AArch64 loops are the generated assembly of the public v9
package. Their per-walk counts are:

| loop | instructions per walk |
| --- | ---: |
| pass 1 (forward product and dx stores) | 111 |
| pass 2 (backward inversion, affine addition, point stores) | 523 |
| **counted total** | **634** |

The loops are single basic blocks: every conditional is a CSEL and the only
branch is the loop back edge, so every iteration executes exactly 111 and
523 instructions. The public v9 evidence parses the counts from the shipped
assembly, executes the shipped file in an AArch64 interpreter on production
points, and checks the per-extra-walk difference; its native test checks the
same loops against big integers and affine secp256k1. The file hashes,
opcode histograms and the identities
`111 = 91 + 9 + 5 + 4 + 2` and
`523 = 4*91 + 71 + 5*9 + 2*8 + 19 + 6 + 2` are in `kernel-counts.md`.
This package does not copy the assembly; it pins the public files.

### 5.2 Analytic bookkeeping sequences

The ordinary step also needs the class bookkeeping that is not inside the
two counted loops. These are explicit instruction sequences, not an
allowance. The rows and their instruction counts are:

| row | instructions | explicit sequence |
| --- | ---: | --- |
| F1 truncated product | 24 | 4 x-limb loads + 5 G-limb loads (9); 4 MUL + 4 UMULH (8); 7 ADD (7) |
| range test, F0, S, F2 | 8 | SUB, CMP, branch (3); ADDS (1); NEG (1); CSEL, CMP, branch (3) |
| decision | 13 | CMP, 2 CSEL, CSET (4); CMP, CSEL (2); CMP, 3 CSEL (4); difference, CMP, branch (3) |
| sign, p-y, K | 13 | SUBS+3 SBCS (4); SUBS+3 SBCS (4); CSET (1); 4 CSEL (4) |
| h, DP test, t, address | 7 | AND (1); CBZ (1); 2 ADD (2); 3 ADD with shifts (3) |
| arrival-pass traffic | 5 | 4 y-limb loads + 1 table-pointer store |
| fingerprints, history, checkpoint | 26 | flag load/CBNZ (2); slot address (2); compare 3 loads/CMP/2 CCMP/branch (7); overwrite 3 stores (3); checkpoint compare 3 loads/CMP/2 CCMP/branch (7); TST/branch (2); checkpoint stores (3) |
| zero-denominator test | 4 | ORR of dx limbs in pairs and together, CBZ |
| length, cap, meter, loop | 13 | 2 length loads, ADDS+ADC, 2 stores (6); cap CMP/branch (2); meter SUBS/ADD/branch (3); loop SUBS/branch (2) |
| batch fixed-work share | 1 | per-batch prologue, inversion glue and epilogue below 2^15 instructions for at least 2^15 walks |
| **analytic total** | **114** | |

The sequences are complete for the ordinary (non-fallback) path. Rare paths
are separately priced and are not included here: the exact-canonical-index
fallback (0.5 units), the cold zero-denominator/abandonment path, the
doubling escape (1.1 units), the traversal step (0.795 units) and reports
(40 units). The table point and walk point loads and stores are inside the
counted loops; the arrival pass adds only the y limbs and the table pointer
that the counted loops do not read or write. The counted loops also carry
their own loop counters; the 13-instruction row's loop SUBS/branch is the
arrival-pass loop, not a second count of the counted loops.

The public v9 proof lists a conservative 129-row table (27, 12, 14, 13, 7,
8, 26, 6, 14, 2). Its own derivation of the explicit sequences is the 114
above, with 15 instructions of slack. The 126.3470 variant charged 763
(129 + 634) after removing the upstream two-instruction reserve. This
package charges the exact 114 + 634 = 748. The instruction sequences, the
counted assembly and the field-operation identities are unchanged; only the
slack in the bookkeeping table is removed.

### 5.3 Prices and the other actions

The step is charged 748 instructions = 0.748 units. The traversal price
remains 0.795 (772 instructions at most, so the price is conservative). The
doubling price remains 1.1: with the 114 bookkeeping rows the doubling
bound is 952 instructions (5 fe_mul, 2 fe_sqr, one 2y addition, two weak
3x^2 additions, four subtractions, two canonicalisations, 72 traffic and
loop instructions, 114 rows, 50 list moves and 30 escape-logic
instructions), below 1,100. The fallback, report, start, sampling,
verification, store and phase-2 prices are unchanged from the public v9
proof.

## 6. Accounting and limitations

The ordinary step is charged at 748 instructions with no separate reference
operation; the counted loops are 111 and 523 per walk, and the analytic rows
are the 114 explicit sequences of section 5.2. The doubling is at most 952
instructions and is priced 1.1. The traversal is at most 772 instructions
and is priced 0.795. All other actions retain the public v9 prices.

The limitations are the public limitations: the counted instruction counts
are those of the shipped generated assembly, checked natively and by an
interpreter but not by a Callgrind measurement of a full attack; the
analytic rows are explicit instruction sequences but are analytic, not
measured; H1 is a random-mapping heuristic supported only at toy scale; and
the attack is a generic analytic bound, not a practical break. The new
step charge is 15 instructions lower than the 129-row table; it is
justified by the explicit sequences above, but a reviewer should check that
no ordinary-path instruction is missing. The fallback, doubling, traversal
and report paths are separately priced, so a missing instruction there
would not be hidden by the 114 rows. The counted 634 is not re-executed in
this package; it is pinned by hash and opcode histogram in
`kernel-counts.md` and is the public v9 assembly.

The package is deliberately proof-sized because the submission archive also
contains the earlier public packages and is capped at 262,144 expanded
bytes. It includes a small standard-library experiment that recomputes the
114-row sum, the 748 price, the 126.3194 cap and the H1/Markov success
floor. It does not execute the counted kernels. The hosted review and
isolated execution remain necessary before any official score exists.

# BLS12-381 G1: Cheon attack with fixed-base comb scalar multiplication

## Claim

For the protected `bls12-381-g1-ecdlp` target (`T=[s]P0`, powers
`A_i=[s^i]P0`, scoring `n=2^20`, `s` uniform in `1..r-1`), this package gives
a deterministic two-phase Cheon attack with success probability 1 and total
work at most

    W(2^20) <= 2^125.329412

target-group-operation equivalents, i.e. the exact declared score
**31332353/250000 = 125.329412**, the ledger total rounded up at the sixth
decimal. Peak memory is at most `2^126` bytes (ledger
`2^125.2358996`), preprocessing is at most `2^124.336` units (ledger
`2^124.3354353`) and is included in the total, and nonuniform advice is zero
bytes. Every stated exponent bound is at least the ledger's exact value; the
score rational is the rounded-up logarithm. The attack uses only `A_1=T` and
`A_d` with `d=2^20`, which divides `r-1` and is supplied at the scoring `n`;
the protected rules explicitly admit Cheon's attack on these inputs.

This is a corrected resubmission. The two previous candidates
(`ea48c674`, `42e28a31`) had sound Cheon algebra and verified GLV pairs but
were rejected on resource accounting: omitted batch workspace, asserted table
and overhead prices, and incomplete doubling instructions. The **fixed-base
comb** below removes the disputed online doublings entirely: every online
point operation is the protected reference addition, and the only doublings
are in the one-time comb precomputation, which is charged in full with a
large complete-instruction allowance.

## Attack core

`r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001`,
`r-1 = 2^32*3*11*19*10177*125527*859267*906349^2*2508409*2529403*52437899*254760293^2`,
`d=2^20`, `M=(r-1)/d` (235 bits). `zeta=7` is a generator of `Z_r^*`; the
experiment verifies the factorisation, primality, and
`zeta^((r-1)/q) != 1` for every prime `q | r-1`. With
`zeta_hat=zeta^d` (order `M`) and `zeta_check=zeta^M` (order `d`), every
secret has the unique form `s=zeta^(k0+M*k1)`, so
`s^d=zeta_hat^k0` and `s*zeta^(-k0)=zeta_check^k1`.

Phase 1 runs BSGS with `m=m_hat=223621884914696853578958571250384896`
(`m+m_hat=447243769829393707157917142500769792`): baby points
`[zeta_hat^(-u)]A_d`, giant points `[zeta_hat^(m*v)]P0`, index reduced mod
`M`. Phase 2 runs the same search for `k1 in [0,d)` with `m'=1024`; then
`s=zeta^k0*zeta_check^k1 mod r` and `[s]P0=T` is verified. Both phases are
exhaustive; success probability is 1.

## Scalar multiplication: fixed-base comb on fixed bases

Phase 1 multiplies the *same base* `A_d` by the `m` different known scalars
`zeta_hat^(-u)`, and the *same base* `P0` by the `m_hat` giant scalars. Every
scalar is a known 255-bit integer, decomposed by GLV into two signed 127-bit
halves using the verified production endomorphism:

    phi(x,y)=(beta*x,y),  beta^3=1 mod p,  beta != 1,
    lambda = 7^((r-1)/3) mod r,  lambda^2+lambda+1 = 0 mod r,

with `lambda` 128 bits and the reduced short lattice

    (-228988810152649578064853576960394133503, 1),
    (1,   228988810152649578064853576960394133504).

`[k]P=[k1]P+[k2]phi(P)` with `|k1|,|k2|<2^127` (verified on random scalars).

Because the base of a multiplication is fixed across the table, each base is
served by a **fixed-base comb table**. For a 127-bit half-scalar and window
`w=16`, the scalar has 8 base-`2^16` digits. For a base `P` precompute

    T[i][d] = [d*2^(16 i)] P,   i=0..7,  d=0..2^16-1,

by `[2^16]`-doublings between blocks and 65,535 sequential additions per
block. Then

    [k]P = sum_{i=0..7} T[i][k_i],

which is at most 7 additions per half (8 table points summed, one final
combination), i.e. at most 17 additions per scalar for the two GLV halves.
There are **no doublings in the online phase**: every online point operation
is an addition of a precomputed affine point. Signs of `k1,k2` are handled by
negating the selected points (free) and adding or subtracting; the final
combination has one addition.

The comb tables are built once per base for the whole attack. The phase-1 bases are `A_d`, `phi(A_d)`, `P0`, `phi(P0)`: four tables of
`8*2^16 = 524,288` entries each (112 bytes per entry).  Phase 2 needs the
runtime point `S=[zeta^(-k0)]T`, which is not known until phase 1 finishes;
the attack therefore builds two further fixed-base tables for `S` and
`phi(S)` at that point and charges their construction (the two tables cannot
reuse the phase-1 storage).  Six tables of 524,288 entries: `2^28.4` bytes.

## Complete online cost per scalar multiplication

* **Additions:** 17 reference additions. The reference is the protected
  batched incomplete affine G1 addition (`6M`, `S=M`, `A=0`), priced at 1
  unit. The ledger charges **33 additions = 33 units** (two-fold margin).
* **Table lookups:** 17 random 112-byte loads with address arithmetic and
  selection. Each is charged **2 units = 4,000 instructions** (a very large
  over-estimate; 33 lookups = 66 units).
* **Overhead:** GLV decomposition (<=1,024 instructions), digit extraction
  (<=512), sign handling and exceptional-case branches for every addition
  (<=1,024), scalar update in `F_r` (<=1,024), loop/index control (<=1,024):
  at most 4,608 instructions = 2.304 units. The ledger charges **16 units**
  (about seven times the itemised bound).

Hence

    E = 33*1 + 33*2 + 16 = 115 units per scalar multiplication.

There is no online doubling, so the earlier objection about unpriced affine
doubling arithmetic does not apply to the online phase.

## One-time comb precomputation (fully charged)

Four bases, 8 blocks each, 65,535 sequential additions per block, and 112
doublings per base:

    6 * (8*65,535 + 112) = 3,146,352 point operations.

Charge every point operation completely and generously: 10,000 instructions
per operation (an addition's reference formula is 6M = 2,298 instructions;
an affine doubling is 7M = 2,681 instructions; 10,000 leaves a factor
>=3.7 for field additions, reductions, normalization, controls and
exceptional-case handling). This is at most

    3.15e10 instructions = 2^34.9 instructions = 2^23.9 units.

The ledger charges **PRECOMB + PHASE2_PRECOMB = 2^30 + 2^29 units**, over
60 times that bound.  The two phase-2 tables for `S` and `phi(S)` are built
after `S` is computed and are charged at `PHASE2_PRECOMB = 2^29` units; they
are included in preprocessing and in the total.
The comb tables are part of preprocessing and are included in the total; they
are built for this target and are not amortized over future targets.
Computation of the fixed-base tables for `P0` and `phi(P0)` is also charged
here even though those data are public profile constants.

## Batching, inversion and workspace

The online additions are processed in batches of `B=2^64` live chains. Within
a round the `B` affine denominators are inverted together by Montgomery's
trick: `3M` per element plus one inversion per round. `B=2^64` makes the
single inversion negligible per element, and the whole operation is the
protected reference addition at 1 unit. The live workspace per chain is

    112 (current point) + 48 (denominator) + 48 (prefix) + 32 (control)
    = 240 bytes,

so the batch workspace is `240*B = 2^71.9` bytes. It is included in peak
memory. No chain state beyond one batch is live at once.

## Table: merge sort and binary search, no hashing

The baby table is an array of `m` records, each a 96-byte affine point and a
**16-byte index** (an index below `m` needs 118 bits; 8 bytes would be too
small). It is sorted by merge sort with two record buffers and
`L=ceil(log2 m)=118` passes.

Per entry per pass: one key comparison (8 limb comparisons and branches, 32
instructions) + one record move (14 load/store instructions) + 16
instructions of index/loop control + a 12-instruction margin = 74
instructions; `118*74 = 8,732` instructions = 4.37 units per entry. The
ledger charges **5 units per baby entry**.

Each giant lookup is a binary search: 119 comparisons, each 32 instructions
for the key plus 16 for addressing/control = 5,712 instructions = 2.86
units. The ledger charges **4 units per giant lookup**.

There is no dictionary, hash function, load factor, collision handling or
expected-time assumption.

## Complete ledger

    E = 115 units per scalar multiplication,
    phase 1 = (m + m_hat)*115 + m*5 + m_hat*4,
    phase 2 = 2048*115 + 8 + 1024*12*48/2000,
    PRECOMB = 2^30 units (phase-1 tables),
    PHASE2_PRECOMB = 2^29 units (runtime S, phi(S) tables),
    fixed   = 2^40 scalar instructions / 2000 units.

The fixed allowance covers reading and validating `T` and `A_d`, the
factorisation and generator test, powers of `zeta`, GLV constants, the final
reconstruction and the verification `[s]P0=T`. The total is

    W = (m+m_hat)*115 + m*5 + m_hat*4 + 2048*115 + 8
        + 1024*12*48/2000 + 2^30 + 2^29 + 2^40/2000
      = 2^125.3294115145043431110420518605144247972915263886905778...,

rounded up at the sixth decimal to `31332353/250000 = 125.329412`. Component
logarithms: phase 1 `125.329411`, phase 2 `17.847`, comb precomputation
`30.000`, fixed `29.034`. The experiment recomputes every term with exact
integers and rationals, checks the high-precision logarithm (60-digit
`Decimal`), and asserts the rounded score.

Sensitivity: even at `E=300` units per scalar multiplication (2.6 times the
declared online price) the total is `2^126.679`, still below `2^127`.

## Memory, preprocessing and advice

* Two merge-sort buffers: `2*m*112 = 2^125.2358996` bytes.
* Six comb tables: `6*524288*112 = 2^28.4` bytes.
* Batch workspace: `2^64*240 = 2^71.9` bytes.
* Fixed reserve: `2^22` bytes.

Peak memory is at most `2^125.2358997` bytes; the claim states the rounded-up
`2^126`. Preprocessing is `m*(115+5) + 2^30 + 2^29 = 2^124.3354353` units,
included in `W`, with no reuse across targets.  The phase-2 tables for `S`
and `phi(S)` are constructed after phase 1 and are charged in this total. Only `A_1` and `A_d` are read from the
powers. Nonuniform advice is zero bytes.

## W(n), validity and limitations

The attack needs `A_d` for `d=2^20`; it applies for every `n >= 2^20` with
`W(n)=W(2^20)` (extra powers are ignored), and is not claimed for
`n < 2^20`. The protected score is at `n=2^20`.

The success probability is exact (deterministic). The cost uses the protected
reference price for additions, a complete 10,000-instruction-per-operation
allowance for the one-time comb precomputation, a 2-unit-per-lookup and
16-unit-overhead allowance, and the explicit merge-sort/binary-search
budgets. Remaining limitations: the batched-inversion convention underlying
the reference addition, the enormous `2^125`-byte table, and the absence of a
production-size execution. Cheon's auxiliary-input attack and fixed-base comb
scalar multiplication are known methods; no novelty is claimed.

## Evidence

`experiment.py` (stdlib only, three labelled parts):

1. **Production ledger**: verifies the factorisation and primality of
   `r-1`, `zeta=7`, the order of `zeta_hat`, the production GLV pair on
   `P0`, `[r]P0=O`, the reduced lattice, the 130-bit decomposition bound on
   64 random scalars, the comb dimensions, and every cost term (including
   the batch workspace and precomputation); checks the high-precision
   logarithm and asserts the score.
2. **Production-curve checks** with the protected `p`, `r`, `P0`.
3. **Toy end-to-end Cheon attacks using the fixed-base comb**: `y^2=x^3+4`
   over `F_100183` (order 4759, `d=2`, six fresh secrets) and a prime-order
   group of order 1032193 (`d=1024`, four fresh secrets), with the complete
   two-phase searches.  Both build comb tables for the phase-1 bases and a
   **runtime comb table for `S`** before phase 2, mirroring the production
   schedule.  A separate comb self-test compares 32 fixed-base and 32
   runtime-base comb evaluations with ordinary scalar multiplication.

Observed on the declared input: score `31332353/250000`, `W_log2 =
125.3294115145043431...`, preprocessing log2 `124.335435302132...`, memory
log2 `125.2358996285812`, 6/6 toy-curve and 4/4 toy-group recoveries, and
32/32 + 32/32 comb self-tests.  The toy searches use the comb primitive,
including the runtime `S` table.
The toy runs are evidence at their scope only; they are not a production
execution.

## Provenance

No Yukon submission is reused. Method: J. H. Cheon, "Security Analysis of the
Strong Diffie-Hellman Problem" (EUROCRYPT 2006); R. P. Gallant, R. J. Lambert,
S. A. Vanstone, "Faster Point Multiplication on Elliptic Curves with Efficient
Endomorphisms" (CRYPTO 2001); fixed-base comb/window scalar multiplication is
standard. The target, powers, `n=2^20` and Cheon admissibility are the
protected `bls12-381-g1-ecdlp` profile and rules.

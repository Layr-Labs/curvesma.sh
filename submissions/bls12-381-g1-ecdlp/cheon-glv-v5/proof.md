# BLS12-381 G1: Cheon attack with fixed-base comb and sequential precomputation accounting

## Claim

For the protected `bls12-381-g1-ecdlp` target (`T=[s]P0`, powers
`A_i=[s^i]P0`, scoring `n=2^20`, `s` uniform in `1..r-1`), this package gives
a deterministic two-phase Cheon attack with success probability 1 and total
work at most

    W(2^20) <= 2^125.329412

target-group-operation equivalents, i.e. the exact declared score
**31332353/250000 = 125.329412**, the ledger total rounded up at the sixth
decimal. Peak memory is at most `2^126` bytes (ledger `2^125.2358996`),
preprocessing is at most `2^124.336` units (ledger `2^124.3354353`) and is
included in the total, and nonuniform advice is zero bytes. Every stated
exponent bound is at least the ledger's exact value; the score rational is
the rounded-up logarithm.

The attack uses only `A_1=T` and `A_d` with `d=2^20`, which divides `r-1` and
is supplied at the scoring `n`; the protected rules explicitly admit Cheon's
attack on these inputs. No prior submission is reused.

This is a corrected resubmission. Earlier candidates had sound Cheon algebra
and a verified GLV pair but were rejected on resource accounting. This
version (i) accounts for the **runtime phase-2 base** `S=[zeta^(-k0)]T` with
two extra comb tables built after phase 1, (ii) removes the unsupported
"10,000 instructions per precomputation point operation" cap, and instead
charges the **sequential dependent construction** with a standalone Fermat
inversion per operation, and (iii) implements the comb primitive in the toy
end-to-end searches, including the runtime `S` table.

## Attack core

`r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001`,
`d=2^20 | r-1`, `M=(r-1)/d` (235 bits). `zeta=7` is a generator of `Z_r^*`
(factorisation, primality and order tests in the experiment). With
`zeta_hat=zeta^d` (order `M`) and `zeta_check=zeta^M` (order `d`), every
secret has `s=zeta^(k0+M*k1)`, hence `s^d=zeta_hat^k0` and
`s*zeta^(-k0)=zeta_check^k1`.

Phase 1 BSGS: `m=m_hat=223621884914696853578958571250384896`, baby points
`[zeta_hat^(-u)]A_d`, giant points `[zeta_hat^(m*v)]P0`, index reduced mod
`M`. Phase 2: compute `S=[zeta^(-k0)]T`, then the same search for
`k1 in [0,d)` with `m'=1024` over `[zeta_check^(-u')]S` and
`[zeta_check^(m'v')]P0`. Then `s=zeta^k0*zeta_check^k1 mod r` and
`[s]P0=T` is verified. Both searches are exhaustive; success probability 1.

## Scalar multiplication: fixed-base comb on fixed bases

Phase 1 multiplies the fixed bases `A_d` and `P0` by known scalars; phase 2
multiplies the runtime base `S` and the fixed base `P0`. Scalars are
GLV-decomposed with the production endomorphism

    phi(x,y)=(beta*x,y), beta^3=1 mod p,
    lambda=7^((r-1)/3) mod r, lambda^2+lambda+1=0 mod r,

`lambda` 128 bits, short lattice basis
`(-228988810152649578064853576960394133503, 1)`,
`(1, 228988810152649578064853576960394133504)`; the experiment verifies
`phi(P0)=[lambda]P0` on the protected curve and the 127-bit half bound.

For a 127-bit half-scalar and window `w=16`, use 8 base-`2^16` digits:
precompute `T[i][d]=[d*2^(16i)]P`, then `[k]P=sum_i T[i][k_i]`, at most 7
additions per half and 17 additions per scalar. **The online phase contains
no doublings**; it is a sum of precomputed affine points. Signs of the GLV
halves are handled by negating selected points (free) and adding/subtracting.

Base tables:

* phase 1: `A_d`, `phi(A_d)`, `P0`, `phi(P0)` (four tables);
* phase 2: the runtime point `S` and `phi(S)` (two further tables, built
  after phase 1; they cannot reuse the phase-1 storage).

Six tables of `8*2^16 = 524,288` entries, 112 bytes per entry
(`2^28.4` bytes total).

## Complete online cost per scalar multiplication

* 17 reference additions (protected batched incomplete affine G1 addition,
  `6M`, 1 unit each). Charged **33 additions = 33 units**.
* 17 random 112-byte table loads with addressing/selection. Charged **2
  units each**; ledger charges 33 lookups = **66 units**.
* Overhead: GLV decomposition (<=1,024 instructions), digit extraction
  (<=512), sign/exceptional-case branches per addition (<=1,024), `F_r`
  scalar update (<=1,024), loop/index control (<=1,024) = <=4,608
  instructions = 2.304 units. Charged **16 units**.

    E = 33 + 66 + 16 = 115 units per scalar multiplication.

There are no online doublings, so no unpriced affine-doubling arithmetic
appears in the online phase.

## Comb precomputation: sequential construction, per-operation inversion

The tables are built by sequential dependent addition chains: each entry is
the previous entry plus a fixed block base, and the block bases are obtained
by doublings. Since the chain is dependent, a denominator must be inverted
before the next update; **no amortized batch inversion is claimed for this
construction**. The construction size is

    6 * (8*65,535 + 112) = 3,146,352 point operations.

Every operation is charged its complete standalone affine cost:

* one Fermat inversion of a 6-limb denominator, at most 640 `F_p`
  multiplications = `640*383 = 245,120` instructions;
* the addition or doubling formula, `6M` or `7M` <= `2,681` instructions;
* field additions, reductions, normalization and control <= `3,000`
  instructions;

so at most `250,801` instructions per point operation. For 3,146,352
operations this is

    7.90e11 instructions = 2^39.52 instructions = 2^28.53 units.

The ledger charges **PRECOMB = 2^31 units** for all six tables, over five
times the itemised bound. The two phase-2 tables for `S` and `phi(S)` are
built after `S` is computed; they are included in this charge. The whole
precomputation is included in preprocessing and in the total and is not
amortized over future targets. The toy implementation mirrors the schedule by
using an explicit modular inversion for every affine addition; the production
bound is analytic, not a measurement.

## Batching, workspace and table

The **online** additions are independent across the table entries, so they
are processed in batches of `B=2^64` live chains and their denominators are
inverted together (Montgomery: 3M per element plus one inversion per round;
`B=2^64` makes the single inversion negligible per element). This is the
protected batched incomplete affine addition priced at 1 unit. The live
workspace per chain is 112 + 48 + 48 + 32 = 240 bytes, so the batch workspace
is `240*2^64 = 2^71.9` bytes, included in peak memory. No chain state beyond
one batch is live.

The baby table is a sorted array of `m` records, each a 96-byte affine point
and a **16-byte index** (an index below `m` needs 118 bits). Merge sort, two
buffers, `L=ceil(log2 m)=118` passes: `118*74 = 8,732` instructions = 4.37
units/entry; charged **5 units/entry**. Binary search: 119 comparisons x 48
= 5,712 instructions = 2.86 units; charged **4 units/lookup**. No hashing,
load factor, collision policy or expected-time assumption.

## Ledger and stated bounds

    E = 115 units per scalar multiplication,
    phase 1 = (m + m_hat)*115 + m*5 + m_hat*4,
    phase 2 = 2048*115 + 8 + 1024*12*48/2000,
    PRECOMB = 2^31 units,
    fixed   = 2^40 scalar instructions / 2000 units.

The fixed allowance covers reading and validating `T` and `A_d`, the
factorisation and generator test, powers of `zeta`, GLV constants, the final
reconstruction and the verification `[s]P0=T`. The total is

    W = (m+m_hat)*115 + m*5 + m_hat*4 + 2048*115 + 8
        + 1024*12*48/2000 + 2^31 + 2^40/2000
      = 2^125.3294115145043431110420518605289169268664223332267480...,

rounded up at the sixth decimal to `31332353/250000 = 125.329412`. Exact
preprocessing `m*(115+5)+2^31 = 2^124.3354353...`, below the stated
`2^124.336`. Exact peak memory `2^125.2358996`, below the stated `2^126`.
The experiment computes the logarithms with 60-digit `Decimal` precision,
asserts both inequalities and the rational `31332353/250000`. Sensitivity:
at `E=300` units per scalar multiplication the total is `2^126.679`, still
below `2^127`.

## Memory, preprocessing and advice

* Two merge-sort buffers: `2*m*112 = 2^125.2358996` bytes.
* Six comb tables: `6*524288*112 = 2^28.4` bytes.
* Online batch workspace: `2^64*240 = 2^71.9` bytes.
* Fixed reserve: `2^22` bytes.

Peak memory is at most `2^125.2358997` bytes; the claim states `2^126`.
Preprocessing is `m*(115+5) + 2^31 = 2^124.3354353` units, included in `W`,
with no reuse across targets. Only `A_1` and `A_d` are read. Nonuniform
advice is zero bytes.

## W(n), validity and limitations

The attack needs `A_d` for `d=2^20`; it applies for every `n >= 2^20` with
`W(n)=W(2^20)` (extra powers are ignored), and is not claimed for
`n < 2^20`. The protected score is at `n=2^20`.

The success probability is exact (deterministic). The cost uses the protected
reference price for online additions (with the stated independent batching
and workspace), a complete sequential per-operation inversion bound for the
six comb tables, a 2-unit lookup and 16-unit overhead allowance, and the
explicit merge-sort/binary-search budgets. Remaining limitations: the
enormous `2^125`-byte table and the absence of a production-size execution.
Cheon's auxiliary-input attack, GLV and fixed-base comb multiplication are
known; no novelty is claimed.

## Evidence

`experiment.py` (stdlib only, four labelled parts):

1. **Production ledger**: verifies the factorisation and primality of `r-1`,
   `zeta=7`, the order of `zeta_hat`, the production GLV pair on `P0`,
   `[r]P0=O`, the reduced lattice, the 130-bit decomposition bound, the six
   comb tables and the sequential precomputation bound, and every cost term;
   checks the high-precision logarithm and asserts the score.
2. **Production-curve checks** with the protected `p`, `r`, `P0`.
3. **Comb self-test**: 32 fixed-base and 32 runtime-base comb evaluations
   compared with ordinary scalar multiplication on a toy curve.
4. **Toy end-to-end Cheon attacks using the comb**, including the runtime `S`
   table: `y^2=x^3+4` over `F_100183` (order 4759, `d=2`, six fresh secrets)
   and a prime-order group of order 1032193 (`d=1024`, four fresh secrets).

Observed on the declared input: score `31332353/250000`, `W_log2 =
125.3294115145043431...`, preprocessing log2 `124.335435302132...`, memory
log2 `125.2358996285812`, 32/32 + 32/32 comb self-tests, 6/6 toy-curve and
4/4 toy-group recoveries. The toy searches use explicit per-addition
inversions, mirroring the sequential precomputation schedule. Toy runs are
evidence at their scope only; they are not a production execution.

## Provenance

No Yukon submission is reused. Method: J. H. Cheon, "Security Analysis of the
Strong Diffie-Hellman Problem" (EUROCRYPT 2006); R. P. Gallant, R. J. Lambert,
S. A. Vanstone, "Faster Point Multiplication on Elliptic Curves with Efficient
Endomorphisms" (CRYPTO 2001); fixed-base comb multiplication is standard. The
target, powers, `n=2^20` and Cheon admissibility are the protected
`bls12-381-g1-ecdlp` profile and rules.

# Randomized-base transfer to ordinary BLS12-381 G1 discrete logarithm

## Status and exact scope

This is an exploratory, conditional candidate, not a demonstrated production
attack or an independently established NFS complexity theorem. The inherited
H1--H5 assumptions, support and limitations in source-core.md are part of this
claim. The kernel-price premise H4 relies on the supplied measurement records,
not a new benchmark of the original C kernels. README.md distinguishes the
reproduced scoped checks from these remaining production assumptions.
No previous review decision is used as a mathematical premise.

The proposed target bound is W <= 2^125.3531834 G1-reference units at C=2,000,
success >=0.39, total peak storage <2^106 bytes, and all preprocessing included
in W. There is no auxiliary power input, external advice or amortization.
A fresh hosted review must assess this complete exploratory claim.

## Supplied core

source-core.md, DESCENT.md, KERNELS.md, STORAGE.md and the accompanying data and
programs supply the STNFS two-log algorithm, rather than relying on an external
verdict. They are copied from submission 7828b2f1-e356-42ea-a363-1cf72a43699d,
immutable commit 6a88a60090ea63e588a8759e67764b80539af772. Their original point
recovery wrapper is NOT the wrapper used here. Retain the GT core, all of its
cold preprocessing, field-basis conversions, both individual logarithms,
instruction and node meters, and every failure term. The original full ledger
upper-bounds this core, so retaining that full ledger overcharges the old
wrapper; no subtraction or unearned cost saving is used.

Let A(g,z) be that core's exponent output log(z)/log(g) modulo the prime
subgroup order r. For g uniform in GT minus its identity and z independently
uniform in GT, the inherited argument asserts success at least 0.391663 and
work less than 2^122.0169 units at C=20,000, conditional on H1--H5.
Its upper-rounded ledger sums to less than this value. The supplied program reproduced this numerical claim in the protected Docker
runner. Its within_claim flag concerns the source ledger; the separate
transfer checks verify the new conversion arithmetic, not the NFS assumptions.

## Algorithm

The input is precisely T=[x]P0 in G1. Use the pinned reduced pairing and public
generators P0,Q0. Put h=e(P0,Q0) and z=e(T,Q0). Sample b uniformly from
{1,...,r-1} using rejection sampling of fresh 255-bit strings. Allow at most
128 attempts, then abort and count failure. Put g=h^b. Run A(g,z); on failure
return failure. Otherwise let a be its scalar and set x'=ba mod r. Check
0<=x'<r and [x']P0=T. Return the canonical 32-byte scalar only if this check
passes. Infinity T is allowed; z=1 can be answered immediately with zero.
No challenge sampling secret is known or used by this algorithm.

## Exact distribution and correctness

Nondegeneracy and prime group order make h a generator. Conditional on
successful b sampling, b is exactly uniform nonzero and independent of x.
The map (b,x) -> (h^b,h^x) is a bijection from
F_r^* x F_r to (GT minus {1}) x GT. Thus (g,z) has exactly the core's required
joint distribution, not merely a uniform marginal. Equivalently, a=x/b is
uniform independently of b, and z=g^a. Recovering a yields x'=ba=x mod r.
The final check prevents incorrect outputs. Handling x=0 separately only
improves the success probability relative to executing the core on z=1.

The protected r satisfies 2^254<r<2^255, so each rejection-sampling attempt
succeeds with probability (r-1)/2^255 >=1/2. Failure in all 128 attempts is at
most 2^-128. Therefore success is at least
(1-2^-128)*0.391663 >0.39 under the same H1--H5 premises. No independence of
the core's internal descent branches beyond its stated H5 is added silently.

## Explicit wrapper budget

Charge the two pairings, exponentiation h^b, the capped sampling, scalar
multiplication/reduction, subgroup/encoding checks and final scalar check.
A deliberately loose scalar implementation bound is 2^80 instructions:

* Represent Fp in seven 64-bit limbs. A carry, comparison or modular
  addition/subtraction uses fewer than 2^12 elementary scalar instructions.
* A 381-bit shift-and-add modular product uses at most 381 iterations with
  at most two such additions, plus shifts and control, below 2^22 instructions.
* Polynomial-basis Fp12 multiplication uses 144 base products and bounded
  reduction by the fixed degree-12 modulus, below 2^31 instructions.
* Every field inversion and final exponentiation needed here can use binary
  powering with exponents of fewer than 4572 bits, hence fewer than 2^14
  Fp12 multiplications. The fixed pairing Miller loop has fewer than 128
  rounds, each implementable with fewer than 2^10 Fp12 operations, including
  point and line arithmetic. Bound each pairing by 2^50 instructions.
* A G1 scalar multiplication has at most 510 doubling/addition steps; use
  projective formulas and a final inversion, with exceptional-point branches.
  Subgroup checks, GT exponentiation, scalar reduction, encodings and sampling
  fit inside another 2^50 instructions using the same elementary bounds.

These bounds leave enormous slack under 2^80. They describe scalar integer
implementations, not uncharged floating-point or SIMD kernels. These elementary bounds do not certify the inherited STNFS kernels. Wrapper storage is bounded by 2^30 bytes, so source storage 2^105
plus wrapper storage is below 2^106. No actual allocation of that astronomical
source storage is claimed.

## Unit conversion

Convert the same core instruction ledger, not a rounded Fp12 operation count:

 W_G1 <= 10 * (101/100) * 2^122.0169 + 2^80/2000 < 2^125.3531834.

The factor ten is exactly 20,000/2,000. In particular this does NOT price an
arbitrary dense Fp12 multiply at ten G1 units; actual kernel instructions are
already in the inherited ledger. All cold preprocessing remains charged.
Preprocessing is bounded by the full declared total, with zero external advice.

## Limitations and checks

The distribution reduction is exact. The claimed production core cost and
success remain exploratory hypotheses, not consequences of that reduction.
Tests of the bijection on small cyclic groups would validate only those finite
checks. Structural preflight would validate packaging only. Neither replaces
the inherited NFS evidence audit or the hosted review of this complete claim.

## Explicit representation conversion supplement

field-map.md supplies a fresh, charged field-isomorphism construction. The enlarged 2^80-instruction wrapper budget includes it. Subtract an additional 2^-48 from the inherited success floor for its capped randomized construction; this still leaves success strictly above 0.39. All earlier wrapper-only bounds are subsumed by this enlarged total, not substituted for the mapping cost.

## Conservative instruction-price padding

Multiply the complete inherited instruction ledger and every corresponding meter allowance by 101/100. This is a cost reserve, not a reduction of the failure probability or an algorithmic speedup. It leaves the same simulated operation path and success law while allowing one percent additional implementation work throughout. Direct rational checks of the supplied differential logs found that some two-decimal prices were rounded down: for example the line4099 ratio is 1914335/61528 >31.11 and record257 is 23449665/522256 >44.90. Both are below 101/100 of their tabulated price, as are all directly checked rows. The reserve covers these observed rounding differences, but does not establish that every unspecified production kernel lies within it; that remains the explicitly disclosed, supported H4 implementation premise. The new declaration deliberately keeps more slack than a rounding-only frontier adjustment.

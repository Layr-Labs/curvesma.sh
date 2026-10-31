# Same-distribution STNFS transfer from recover-Q to recover-P

This exploratory candidate proposes W<=2^122.04 dense-Fp12 reference units, success>=0.39,
peak memory<2^106 bytes, all preprocessing included and no external advice.
It is conditional on all H1--H5 premises reproduced in source-core.md. Those
premises, especially kernel pricing, remain explicit limitations; see README.md.
No inherited review decision is a premise and no production run is claimed.

## Full supplied core and provenance

The files source-core.md, DESCENT.md, KERNELS.md, STORAGE.md and their supplied
data/code reproduce submission 7828b2f1-e356-42ea-a363-1cf72a43699d at immutable
commit 6a88a60090ea63e588a8759e67764b80539af772. Use its complete cold STNFS
two-log core, including field representation conversion, relations, filtering,
linear algebra, both individual logarithms, instruction/node caps and all
failure terms. The old point-output wrapper is not this track's wrapper.
Keeping its full cost ledger upper-bounds the core without subtracting work.
The upper-rounded ledger is below 2^122.0169 and supported success is at least
0.391663 under its declared assumptions. The source ledger and its scoped checks were reproduced in the protected Docker runner.

## Algorithm for the supplied Q

For the exact supplied nonidentity Q in G2 and independent uniform z in GT,
compute g=e(P0,Q), with pinned library argument order pairing(Q,P0).
Run the supplied GT core to find a=log_g(z). Return P=[a]P0, canonically
encoded, after checking e(P,Q)=z. Treat z=1 as P=infinity. Never replace Q
with the generator or a solver-chosen point.

## Distribution-preserving proof

Let h=e(P0,Q0). In the original recover-Q problem write its uniform supplied
point as B=[b]P0, b uniform in F_r^*. Then g=e(B,Q0)=h^b. In the present
recover-P problem write Q=[b]Q0, with exactly the same uniform distribution
of b; g=e(P0,Q)=h^b. In each case z is independently uniform in GT. Therefore
the entire pair (g,z) has the identical distribution, as do the core's fresh
random coins. It follows that the same GT-only core has the identical success
law and same work distribution, not just similar-looking average costs.
If a is correct, bilinearity proves e([a]P0,Q)=g^a=z. Conversely the prime
order and nondegeneracy give uniqueness. The infinity case is included.

This reduction does not claim that G1 and G2 arithmetic cost the same. All
orientation-specific work is charged separately below. It is essential that
the transferred core receives only (g,z), not the sampling exponent b or any
source-group discrete logarithm.

## Conservative orientation-specific budget

Budget 2^80 scalar instructions for two full pairings, G1 scalar multiplication,
subgroup and output checks, scalar normalization and serialization. One
explicit conservative implementation uses seven 64-bit limbs for Fp, modular
additions/comparisons below 2^12 instructions, 381-round shift-and-add field
products below 2^22, and degree-12 polynomial multiplication/reduction below
2^31. Inversion and final exponentiation use binary powering with fewer than
4572 bits, hence fewer than 2^14 Fp12 multiplications. A fixed Miller loop has
fewer than 128 rounds, with fewer than 2^10 Fp12 operations per round using
projective formulas and line evaluations. Each pairing fits below 2^50
instructions. At most 510 projective doubling/addition steps suffice for a
255-bit scalar multiplication; the remaining checks also fit below 2^50.
The total fits strictly below 2^80. This is a deliberately loose analytic budget based on elementary algorithms,
not a measured fast kernel.

At C=20,000 this overhead is below 2^66 reference units, and

 W <= (101/100) * 2^122.0169 + 2^66 < 2^122.04.

The source success floor 0.391663 exceeds 0.39 without reducing any of its
failure allowances. Source storage<=2^105 plus wrapper storage<=2^30 is
below 2^106 bytes. All preprocessing, field-basis conversion and both logs
remain in the ledger; no reusable relation database is free.

## Scientific limitations

The transfer identity is exact. It neither proves the NFS heuristics nor
repairs an unsupported inherited instruction price. The declared source computational checks and transfer checks were reproduced
under the protected credential-free runner. Their success does not resolve
every production-level assumption or make missing C kernel source available.
No local structural check or earlier AI-valid result establishes those facts.
This candidate must receive its own hosted evaluation before any score claim.

## Explicit representation conversion supplement

field-map.md supplies a fresh, charged field-isomorphism construction. The enlarged 2^80-instruction wrapper budget includes it. Subtract an additional 2^-48 from the inherited success floor for its capped randomized construction; this still leaves success strictly above 0.39. All earlier wrapper-only bounds are subsumed by this enlarged total, not substituted for the mapping cost.

## Conservative instruction-price padding

Multiply the complete inherited instruction ledger and every corresponding meter allowance by 101/100. This is a cost reserve, not a reduction of the failure probability or an algorithmic speedup. It leaves the same simulated operation path and success law while allowing one percent additional implementation work throughout. Direct rational checks of the supplied differential logs found that some two-decimal prices were rounded down: for example the line4099 ratio is 1914335/61528 >31.11 and record257 is 23449665/522256 >44.90. Both are below 101/100 of their tabulated price, as are all directly checked rows. The reserve covers these observed rounding differences, but does not establish that every unspecified production kernel lies within it; that remains the explicitly disclosed, supported H4 implementation premise. The new declaration deliberately keeps more slack than a rounding-only frontier adjustment.

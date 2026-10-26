# Fixed P in G1, recover Q in G2

## Claim
Given nonidentity P in G1 and z=e(P,Q) in GT, return Q in G2 with e(P,Q)=z. Return the recovered point with probability 1 using at most
2^131.345710 reference units (dense Fp12 multiplication, C=20000), peak
memory below 2^138 bytes, preprocessing included, advice zero.

## Algorithm
The pairing is nondegenerate and the fixed argument is nonzero, so the map is
an isomorphism onto the order-r subgroup of GT. The target has a unique log
k; the identity target gives k=0 and the infinity point, which the acceptance
allows.

## Exact work
m=ceil(sqrt(r)); L=ceil(log2 m)=128. Compute g, the m GT babies g^j by dense
Fp12 multiplication (C=20000 each), and the m giants z*g^(-i*m) by the same
operation. Sort the m 592-byte records (576-byte GT value plus 16-byte j) by
the 12 Fp coefficients and binary-search each giant.
- Merge sort: m*L comparisons at 512 instructions each (72 limbs, full
  load/compare/branch/index) and m*L moves at 1024 instructions each (592
  bytes, 74 words, load/store/index).
- Binary search: m*L comparisons at 512 instructions each.
- Sort/search = 2048*m*L instructions = 13.1072m reference units at C=20000.
- Group work 2m plus 2^24 (base pairing, output scalar multiplication,
  validation, verification). Total = 15.1072m+2^24 units;
  log2=131.345710 rounded up.

## Memory
m records of 592 bytes plus equal merge scratch = 1184m bytes, below 2^138
bytes (the exact log2 is 137.63799807, and the ceiling is 138, so the bound
holds before any code). Index j is in every record; code/stack below 2^24
bytes are covered by the ceiling margin.

## Correctness
k=i*m+j is unique because r<=m^2; the giant equals the baby exactly at the
target log; the final pairing check is exact. The sort/search are
deterministic and success is 1.

## Limitations
Memory is exponential (1184m bytes), reported separately. The score uses the
complete deterministic sort/search charge. Generic analytic baseline, not a
production execution.

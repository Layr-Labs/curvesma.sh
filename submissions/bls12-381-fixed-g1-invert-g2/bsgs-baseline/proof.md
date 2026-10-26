# Fixed P in G1, recover Q in G2

## Claim

For nonidentity P in G1 and z=e(P,Q) in GT, return Q with probability 1
using at most 2^129.828209 reference Fp12 multiplications, memory below
2^137.637998 bytes, preprocessing included, advice zero.

The pairing is nondegenerate and the fixed argument is nonzero, so
X -> e(P,X) (or X -> e(X,Q)) is an isomorphism onto the order-r subgroup of
GT. The target z lies in that subgroup, so k with z=g^k exists and is unique.

## Algorithm

Let m=ceil(sqrt(r)) and L=ceil(log2 m)=128. Compute the base
g=e(P,Q0) (recover Q) or g=e(P0,Q) (recover P). Generate the m GT baby
elements g^j, j=0..m-1, by dense Fp12 multiplication (one reference unit each
at C=20000), and store each as a 592-byte record: 576-byte GT value plus a
16-byte index j. Sort lexicographically by the 12 Fp coefficients. For
i=0..m-1 compute z*g^(-i*m) and binary search. A match gives k=i*m+j; return
Q=[k]Q0 (recover Q) or P=[k]P0 (recover P), and verify e(P,Q)=z.

## Exact work

- 2m reference Fp12 multiplications: m baby, m giant.
- Merge sort: m*L comparisons (128 instructions each) and m*L moves
  (256 instructions each, 592-byte records).
- Binary search: m*L comparisons, 128 instructions each.
- Sort/search total = 512*m*L instructions; at C=20000 this is 3.2768m units.
- Total = 5.2768m + 2^24 units (base pairing, output scalar multiplication,
  validation and verification); log2 = 129.828209 rounded up. The 2^24 is
  included and is below the rounding margin.

## Memory

m records (592m bytes) plus same-size merge scratch = 1184m bytes, below
{pair_mem}. Index j is stored in each record; no hash metadata. Code/stack
below 2^24 bytes is covered by the rounding allowance.

## Correctness

For every k in [0,r) the representation k=i*m+j with 0<=i,j<m is unique
because r<=m^2. The giant value equals the baby element exactly when k is the
target log; the binary search returns it. The final pairing check is exact.
The identity z=1 gives k=0 and the infinity point, which the acceptance
explicitly allows. The sort/search are deterministic; success is 1.

## Limitations

Generic analytic baseline, not a production execution. Memory is exponential
(1184m bytes), reported separately. The score uses the complete deterministic
sort/search charge, not an expected hash bound. No heuristic is used.

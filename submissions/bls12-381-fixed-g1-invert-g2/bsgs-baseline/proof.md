# Fixed P in G1, recover Q in G2

## Claim

For nonidentity P in G1 and z=e(P,Q) in the order-r subgroup of GT, this
attack returns Q with probability 1 using at most 2^129.428545
reference Fp12 multiplications, with memory below 2^136.598470 bytes and
zero advice; preprocessing is included.

The pairing e is nondegenerate and P (or Q) is nonzero, so the map
X -> e(P,X) (or X -> e(X,Q)) is an isomorphism from the 255-bit group to the
order-r subgroup of GT. The target z lies in that subgroup, so the discrete
log k exists and is unique modulo r.

## Attack

Let r be the subgroup order and m=ceil(sqrt(r)). Compute the base
g=e(P,Q0) (recover Q) or g=e(P0,Q) (recover P), and write z=g^k with
0 <= k < r. Build the baby table {g^j : j=0..m-1} in a hash table keyed by
the 12 Fp12 coefficients. For i=0..m-1 compute z * g^(-i*m) and look it up;
the match gives k=i*m+j with 0 <= i,j < m. Return Q=[k]Q0 (recover Q) or
P=[k]P0 (recover P), then verify e(P,Q)=z.

## Cost

- preprocessing: m baby GT multiplications and m hash insertions;
- online: m giant GT multiplications, m hash lookups, one pairing for the
  base, one scalar multiplication for the output and one verification pairing;
- total charged work <= 4m + 2^{20} < 2^{129.428545} reference units, where
  one unit is one dense Fp12 multiplication (C=20000 instructions).
Memory: m GT elements, below 2^136.598470 bytes. Advice: zero.
Success is 1 (deterministic); the hash table is expected O(1) per operation,
and a deterministic sorted table would only add an O(log m) factor included
in the conservative 4m charge.

## Evidence

This is an analytic baseline, not a production execution. The pairing
definition and the GT group order are the protected profile's; the DLP and
the final scalar multiplication are standard. No heuristic is used for the
recovery itself.

## Limitations

Memory is exponential (m GT elements), as in any BSGS; the track separates
memory from the computation score. The attack is generic and not a practical
break.

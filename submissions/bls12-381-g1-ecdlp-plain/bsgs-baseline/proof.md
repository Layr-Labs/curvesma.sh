# BLS12-381 G1 plain DLP: fully charged affine sorted-table BSGS

## Claim
For uniform x in [0,r) and T=[x]P0 in BLS12-381 G1, recover x with
probability 1 using at most 2^131.628924 reference units. One unit is the
protected 6M batched incomplete-affine G1 addition (C=2000). Peak memory is
below 2^136 bytes, preprocessing is included, advice is zero.

## Parameters
r is the subgroup order; m=ceil(sqrt(r))=
228988810152649578064853576960394133504; L=ceil(log2 m)=128; batch size
B=2^12.

## Affine batch schedule (no projective arithmetic)
The reference unit is the protected 6M batched incomplete-affine addition; we
compute the sequences entirely in affine coordinates and batch the dependent
denominators.

Baby sequence. Start from P1=P0 and Q=P0. For a batch of B consecutive
additions P_{t+1}=P_t+Q with affine P_t=(x_t,y_t) and Q=(x_Q,y_Q), form
d_t=x_t-x_Q for t=1..B. Batch-invert all d_t with a product tree, one field
inversion and prefix/suffix products. Then for each t:
  lambda_t=(y_t-y_Q)*d_t^{-1}        (1M)
  x_{t+1}=lambda_t^2-x_t-x_Q         (1S)
  y_{t+1}=lambda_t*(x_t-x_{t+1})-y_t (1M)
The amortized batch inversion is 3M per point, so the per-point cost is
5M+1S=6M under the 6M convention, exactly the protected reference price. If
d_t=0 (P_t=+/-Q, i.e. a doubling or an inverse), the point is emitted by the
reference doubling/identity path at the same C=2000 unit and the batch is
split around that index. The sequence is ordered multiples, so exceptional
indices are known in advance. The batch scratch is B denominators, B
inverses, B lambda values and B points (about 4B field elements/points),
below 2^20 bytes. j=0 is the identity and j=1 is P0; the doubling 2P0 is an
exceptional path. No projective recurrence and no projective-to-affine
conversion are used.

Giant sequence. Start from T and use the fixed difference Q=-[m]P0. The same
affine recurrence, batch inversion and exceptional handling generate
T-[i*m]P0 for i=0..m-1 at the same 6M price.

Table. Store each baby jP0 as a 112-byte record: 96-byte affine point plus a
16-byte index j. Sort lexicographically by the 12 coordinate limbs. For each
giant, binary-search; a match returns (i,j); output x=i*m+j and verify
[x]P0=T.

## Exact work
- 2m reference 6M additions: m baby, m giant.
- Merge sort: m*L comparisons at 64 instructions each and m*L moves at 128
  instructions each.
- Binary search: m*L comparisons at 64 instructions each.
- Sort/search = 256*m*L = 16.384m units at C=2000.
- Total = 18.384m units; log2=131.628924 rounded up. The 2^24 instructions
  for validation, batch control, final scalar multiplication and verification
  are below the rounding margin.

## Memory
m 112-byte records plus same-size merge scratch = 224m bytes; affine batch
scratch below 2^20 bytes; code/stack below 2^24 bytes. The total is below
2^136 bytes. The 16-byte j is in every record; no hash metadata is used.

## Correctness
x=i*m+j is unique for x<r because r<=m^2; the giant equals the baby record
exactly at the target; the binary search returns (i,j); the final scalar check
handles x=0 and rejects noncanonical input. The schedule, sort and search are
deterministic; success is 1.

## Limitations
Memory is exponential (about 2^135.24 bytes), reported separately. The score
is the complete deterministic 18.384m charge, not an expected hash bound.
This is a generic analytic baseline, not a production execution.

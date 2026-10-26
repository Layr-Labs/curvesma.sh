# BLS12-381 G1 plain DLP: deterministic sorted-table BSGS

## Claim

For uniform x in [0,r) and T=[x]P0 in BLS12-381 G1, recover x with
probability 1 using at most 2^130.777910 target-group-operation equivalents
(C=2000; one unit is one reference batched-affine G1 addition). Memory is
below 2^135.235900 bytes, preprocessing is included, advice is zero.

## Algorithm

Let r be the subgroup order and m=ceil(sqrt(r)) =
228988810152649578064853576960394133504. Generate the m affine baby points
jP0, j=0..m-1, by the reference batched-affine G1 addition (one reference
unit per addition; the batch inversion is part of that protected price), and
store each as a 112-byte record: 96-byte affine point plus a 16-byte index j.
Sort the records lexicographically by the point's 12 64-bit coordinate limbs.
For i=0..m-1 generate T-[i*m]P0 by the same reference addition and binary
search its record. A match returns (i,j); output x=i*m+j and verify [x]P0=T.

## Exact work

- 2m reference G1 additions: m baby, m giant.
- Merge sort: at most m*L comparisons and m*L moves with L=ceil(log2 m)=128.
  Charge 32 instructions per 12-limb comparison and 64 per 112-byte move.
- Binary search: at most m*L comparisons, 32 instructions each.
- Sort/search total = 128*m*L instructions; at C=2000 this is 8.192m units.
- Total = 2m + 8.192m = 10.192m units, log2 = 130.777910 (rounded up at the
  sixth decimal). The 2^24 instructions for validation, final scalar
  multiplication, output and verification are below the rounding margin and
  are included in the rounded-up score.

## Memory

The m records (112m bytes) plus a same-size merge-sort scratch (112m bytes)
give 224m bytes, below 2^135.235900. Index j is stored in each record; no hash
metadata is needed. Code and stack are below 2^24 bytes and are covered by
the rounding allowance.

## Correctness

For every x in [0,r) the representation x=i*m+j with 0<=i,j<m is unique
because r<=m^2. The giant value equals the baby record exactly when x is the
target; the binary search returns the unique (i,j). The final check [x]P0=T
handles the identity (x=0) and rejects any noncanonical input. The sort and
search are deterministic and the success probability is 1.

## Limitations

This is a generic analytic baseline, not a production execution. Memory is
exponential (224m bytes), as in any BSGS; the track reports memory separately.
The score is the complete deterministic 4m/10.192m charge, not an expected-time
hash bound. No heuristic is used for recovery.

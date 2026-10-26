# BLS12-381 G1 plain DLP: fully charged sorted-table BSGS

## Claim
For uniform x in [0,r) and T=[x]P0 in BLS12-381 G1, recover x with
probability 1 using at most 2^131.628924 reference units (C=2000; one unit is
one protected batched-incomplete-affine G1 addition). Peak memory is below
2^136 bytes, preprocessing is included, advice is zero.

## Parameters and schedule
r is the subgroup order; m=ceil(sqrt(r))=
228988810152649578064853576960394133504; L=ceil(log2 m)=128. The baby
sequence jP0 (j=0..m-1) and the giant sequence T-[i*m]P0 are generated in
batches of B=2^12 consecutive points. Each batch uses the reference
batched-incomplete-affine routine: a fixed projective recurrence over the
batch, a product tree of the B z-coordinates, one inversion of the product,
prefix/suffix products to obtain all B inverses, and affine normalization.
The routine detects z=0 (doubling or infinity) before inversion, emits that
point through the reference doubling/identity path at the same C=2000 unit,
and restarts the batch after the exceptional index; sequence points are
ordered multiples, so the exceptional pattern is known and deterministic.
The batch setup (one inversion, 3B multiplications, prefix/suffix, checks)
and the B-point scratch are included in the reference price and in the
memory below; the C=2000 unit is not inferred for an arbitrary dependent
sequence but is applied to this explicit schedule.

## Exact work
- 2m reference G1 additions: m baby, m giant.
- Merge sort: at most m*L comparisons at 64 instructions each (12 limbs
  compared with load/cmp/branch/index) and m*L moves at 128 instructions
  each (112-byte record, 14 words, load/store/index).
- Binary search: at most m*L comparisons at 64 instructions each.
- Sort/search total = 256*m*L instructions = 16.384m reference units at
  C=2000. Total = 2m+16.384m = 18.384m units; log2=131.628924 rounded up.
- The 2^24 instructions for validation, batch control, final scalar
  multiplication and verification are below the rounding margin and are
  included in the rounded-up score.

## Memory
m records of 112 bytes (96-byte affine point plus 16-byte index j) plus a
same-size merge scratch = 224m bytes, plus the batch scratch (2^12 points,
2^12 field elements and one inversion buffer; below 2^20 bytes) and code and
stack below 2^24 bytes. Total is below 2^136 bytes. No hash metadata is used.

## Correctness
x=i*m+j with 0<=i,j<m is unique for x<r because r<=m^2. The giant point
equals the baby record exactly when x is the target; the binary search
returns the unique (i,j). The final check [x]P0=T handles x=0 and rejects
noncanonical inputs. The schedule, sort and search are deterministic; success
is 1.

## Limitations
Memory is exponential (about 2^135.24 bytes), reported separately. The score
is the complete deterministic 18.384m charge, not an expected hash bound. This
is a generic analytic baseline, not a production execution or practical break.

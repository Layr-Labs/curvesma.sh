# Correct scalar-action Cheon grids (no walk)

At n=2^20 claim log2 W < 134, C=2000, success >=81/200. For
2^20<=n<=2^32 use d=2^20, q=(r-1)/d and the protected scalar r.
Only A_1 and A_d are used; traversal is charged.

## Algebra and exact coverage

The complete factorization in verify.py certifies that xi=7 is primitive
modulo r. Set zeta=7^d, delta=7^q. Write uniquely
s=7^(q*k1+k2), 0<=k1<d, 0<=k2<q; then s^d=zeta^k2.
Let m0=ceil(sqrt(2q/5)) and m=ceil(108*m0/100). Store baby records
B_j=[zeta^j]P0 (0<=j<m) and giant records G_i=[zeta^(-mi)]A_d
(0<=i<m). Equality is equivalent to k2=mi+j modulo q. Since m^2<q,
there is exactly one match iff 0<=k2<m^2. These are scalar-action
updates [zeta]B and [zeta^(-m)]G, NEVER fixed-point additions.

After a match, form U=[7^(-k2)]T. Run the same construction with
base delta, target U, and m2=1024. Its exact d-entry exponent grid
recovers k1. Return canonical s=7^(q*k1+k2) mod r only after [s]P0=T.
If phase 1 has no match, fail within the same fixed work budget.
All nonzero secrets have unique exponent pairs, so success is exactly
m^2/q >= (108/100)^2*(2/5)=1458/3125 >81/200. There are no heuristic
collision losses or retries; nu=1.08 uses integer ceilings.

## Instruction ceilings (not one reference unit per addition)

All work is R/C, H=0. A conservative six-limb Montgomery Fp product
costs <=2048 scalar instructions: 36 schoolbook terms plus 36 reduction
terms at <=20 instructions/term (two multiply halves, carry operations,
loads/stores), plus <=608 for loop, carries, normalization and setup.
Reduced Fp addition/subtraction costs <=64 including memory. These
ceilings exceed the 383-instruction calibration; C remains 2000.

A complete Jacobian operation costs <=64 products/squares, 256 field
adds and 4096 other instructions. Compute ordinary addition and doubling
candidates, select infinity/double/general outcomes on H/R/Z tests;
standard Jacobian addition uses 12M+4S and doubling less than 16M/S,
so even computing both is below 64. Zero/inverse/identity cases are
selected, never inverted. Binary scalar action has <=255 doublings and
255 additions. Final affine normalization uses Fermat inversion and
conversion with <=1024 products (p has 381 bits). Charge per action
E=510*(64*2048+256*64+4096)+1024*2048+2^18=79650816.
The last term covers scalar bookkeeping and encoding; no affine batching.

Scalars are reduced after EVERY operation: four 64-bit limbs hold
0<=a<r<2^255; a product is <2^510 (eight limbs), with nine limbs
allocated as reduction scratch. No escape doublings or unreduced sums
exist. Jacobian Fp values have six reduced limbs; product/reduction
scratch uses at most thirteen. All this is inside the charged kernels.

## Deterministic sorting, routing and hard memory

Use exact compressed keys: 48-byte x plus sign/identity flags padded to
56 bytes, then 16-byte index and 8-byte side: 80 bytes/record. Equality
of keys is exact point equality, with canonical infinity/sign handling.
N=2m+2048 bounds both phases together. Two buffers have L slots each,
L=next power of two >=N. Each buffer's binary routing tree has L-1
64-byte routers. Child addresses use three 64-bit limbs. Thus total
peak <=288L+2^64 <2^129 bytes, including code, input and scratch.
Capacity is fixed; N<=L, so overflow probability is EXACTLY 0.

Stable LSD radix sort on all 56 key bytes takes 56 passes. Each pass
reads for counts, reads for scatter, writes once: 168 transfers/record.
Allow eight more for generation/final adjacency scan. Each 80-byte
transfer uses <=129 levels at <=256 instructions/level, including both
request/reply, ten payload words, three-word address and router control.
Charge 5812224 instructions/record. Additionally 65536/record covers
192-bit bucket counts/cursors, prefixes, key tests and router creation.
Fixed 256-bucket setup for 56 passes is below the cold reserve.
Sort then scan adjacent cross-side matches; no hash or probe assumption.

The complete bound is
R(n)=N*(E+5812224+65536)+2^80+2^20*n; W(n)=R(n)/2000.
The cold 2^80 reserve covers constant-count modular exponentiations,
root tests, seeds, output verification and setup (each requires far
less than 2^20 four/six-limb products or 2^14 scalar actions). Input
scanning/decoding costs <=2^20 instructions/point in the stated n range.
Allocation/initialization and all transfers are charged above. Even
6*R(2^20)/(5*2000)<2^134: 20% repricing of EVERY row still fits.

## Reproduction and scope

Run python3 verify.py --test, python3 verify.py --toy, and python3
verify.py. certificate.json is the last command's exact output; toy.json
captures 989/2000 successful real-curve trials, no wrong outputs, cycles
or escapes. Toy Fp prime=884600761357, subgroup prime=193863853,
y^2=x^3+3, d=24367, q=7956; a six-element automorphism orbit is checked.
Toy coverage is 3844/7956. This tests algebra, not production costs; toy
hashing is replaced by radix sorting. dependencies.json is empty.

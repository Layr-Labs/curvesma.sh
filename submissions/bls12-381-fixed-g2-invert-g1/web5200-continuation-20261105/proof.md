# bls12-381-fixed-g2-invert-g1: finite generic inversion for the supplied fixed argument

## 1. Exact claim and protected target

The input is the supplied nonidentity Q and an independently uniform
z in GT. The output is a point P in G1, satisfying e(P,Q)=z
for the EXACT supplied fixed argument. We do not replace that argument
by a generator and do not return a scalar instead of the required point.
There are no auxiliary powers or n parameter. The deterministic algorithm
has success probability exactly 2^254/r > 1/2, total classical work below
2^151 protected group-operation equivalents (C=20000), and peak memory
below 2^146 bytes. All work is cold, including failures. No novelty,
competitive improvement, practical break or execution is claimed.

Use these exact protected constants and representations:

```json
{
  "field_modulus": "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab",
  "subgroup_order": "73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001",
  "curve_equation_G1": "y^2 = x^3 + 4 over Fp; G1=<P0> of order r",
  "curve_equation_G2": "y^2 = x^3 + 4*(1+u) over Fp2; G2=<Q0> of order r",
  "generator_P0": {
    "x": "17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb",
    "y": "08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1"
  },
  "generator_Q0": {
    "x": [
      "024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8",
      "13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e"
    ],
    "y": [
      "0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160d12c923ac9cc3baca289e193548608b82801",
      "0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be"
    ]
  },
  "field_Fp2": "Fp[u]/(u^2+1); [c0,c1] means c0+c1*u",
  "field_Fp12": "Fp[w]/(w^12-2*w^6+2); coefficient list [c0,...,c11] means sum(ci*w^i); u maps to w^6-1",
  "group_GT": "Nonzero z in Fp12 with z^r=1; identity [1,0,...,0]"
}
```

In particular 2^254<r<2^255, r is odd and prime, and m=2^127<r.
The pairing is the reduced map py_ecc.optimized_bls12_381.pairing(Q,P,
final_exponentiate=True), from commit 04151f01f59f902ab932a51e0ca0ebce3883fc51.
Its argument order is G2 FIRST. The field coefficient basis, original
Miller parameter 0xd201000000010000 and final exponent (p^12-1)/r are
preserved. Comparing raw Miller values is neither required nor substituted.
pairing-kernels.md reproduces the relevant formula path and its cost.

## 2. Input, output and explicit interfaces

Use compact canonical UTF-8 JSON, exact outer input keys Q,z,
in that order, and no redundant whitespace or escapes. Point objects use
x then y. Every base-field coefficient has 96 lowercase hex characters
and is less than p; an Fp2 coordinate is [c0,c1], and z contains exactly
12 coefficients in increasing polynomial degree. This covers every
mathematical challenge and has fewer than 4096 bytes. Reject duplicate
or extra keys, noncanonical coefficients, an infinity input point,
off-curve input and [r]Q != O. Check z!=0 and z^r=1 in Fp12.
The prime-order membership checks use actual scalar/group arithmetic,
not a curve-membership shortcut. The output has the sole key P
and the protected G1 encoding. Identity is null and is returned
exactly when the recovered exponent is zero; it is valid for z=1.

The bounded serializer and error-free bounded-degree local-mailbox/reset
machine described in arithmetic.md are explicit interface assumptions.
They are not authenticated claims about platform transport or available
hardware. An expanded wire spelling of L bytes adds at most 128L+2^16
normalization instructions; arbitrary padding is not covered by the
constant input-processing allowance. There is no restriction on the
uniform point or scalar distribution. All memory, instructions, messages,
blocking entries and initialization within the declared machine are charged.

## 3. Algorithm

Validate the input as above. If z=1, return the identity of G1,
after the stated final checks. Otherwise compute g=e(P0,Q) using the
supplied Q and the public generator P0. Check g!=1 and g^r=1;
these hold for the protected nondegenerate pairing on the stated subgroup.
Set m=h=2^127. No exponent or sampling coin is available to the algorithm.

1. Put v=1. For j=0,...,m-1, store (v,j) as a fixed 1024-byte record
   and update v=v*g. Thus the records are (g^j,j), and the last update
   retains g^m. The key is all twelve canonical Fp coefficients, not
   a hash or a compressed subset. Build the table anew for this input.
2. Sort the m records by their 72-limb lexicographic key with bottom-up
   merge sort and a second m-record buffer. There are 127 passes.
3. Set D=(g^m)^(r-1) and gamma=z. Since g^m is in GT, D=g^(-m).
   For i=0,...,h-1, binary-search the full gamma key. On a match with j,
   set a=i*m+j. If no match and i<h-1, update gamma=gamma*D.
4. On a match check 0<=a<r, compute P=[a]P0, normalize
   it to the required affine representation, check its curve and subgroup,
   and recompute the exact e(P,Q). Return the encoded P if it
   equals z on ALL 12 coefficients; otherwise terminate with failure.
   If the loop exhausts h iterations, return failure. There is no restart.

The recovered integer a is below mh=2^254<r, so no scalar reduction
oracle is needed. At most one final candidate verification occurs; the
proof below establishes that a real table hit always passes. Unexpected
verification failure terminates and cannot cause an uncharged repeated
attempt. Intermediate exponent a is used to form the required point,
not substituted for the output. Input Q is retained unchanged.

## 4. Correctness and exact joint success

Use bilinearity and nondegeneracy of the protected BLS12-381 reduced
pairing on its prime-order G1 and G2. These are the standard target-group
premises underlying the profile's explicit admissible reduction. For
clarity the needed group consequence is proved here. Write the supplied
nonidentity fixed point as [b] times its generator, 1<=b<r. If
g0=e(P0,Q0), nondegeneracy and prime-order cyclicity give order(g0)=r.
Then g=g0^b also has order r, since r does not divide b. Bilinearity gives
e(P,[a]Q0)=g^a in the fixed-G1 case, or e([a]P0,Q)=g^a in fixed-G2.
Thus recovering a in GT and returning the indicated point solves the
exact supplied-argument inversion. This does not require learning b.
The reproduced kernels bind the computation to the target definition;
their operation count alone is not claimed to prove pairing nondegeneracy.

For each fixed nonidentity input, there is a unique x in [0,r) with
z=g^x; it is uniform because z is uniform in GT. The baby keys are
distinct, since g^j=g^k and 0<=j,k<m<r implies j=k. Induction gives
gamma=z*g^(-mi). A hit therefore implies x=mi+j modulo r. All tested
integers mi+j are in [0,mh)=[0,2^254), where reduction is injective.
Conversely every x in this interval has unique quotient i<h and
remainder j<m and is found. There are exactly 2^254 successful x,
including zero, so the probability for EACH supplied nonidentity point is

    2^254/r > 1/2 > 0.39.

Averaging this same conditional probability over the independently
uniform supplied nonidentity point proves the required joint success.
There is no conditioning of the challenge on small exponents; values
outside the interval exhaust the charged budget and fail. No randomness,
random-oracle heuristic, auxiliary powers or success-only cost is used.

## 5. Resources, dependencies and limits

arithmetic.md supplies finite six-limb field, complete pairing, table,
routing, comparison, sorting, input, output and cold-start bounds. With
H=0 separately priced reference operations, R<2^165 scalar instructions
and C=20000>2^14 imply W=R/20000<2^151. Preprocessing is included and
has the same upper bound; searched nonuniform advice is zero. Peak
memory, including two full tables, processor and routing state, replicated
code, buffers and controller storage, is below 2^146 bytes. These
astronomical resources are theoretical, not an assertion of available RAM.

We adapt our own earlier plain-G1 entry e25b1379-d936-4940-9739-d4d90fe30718
for scalar-limb arithmetic and charged storage-interface exposition.
Its judgment is not a premise. The organizer's generic GT baseline and
pairing-kernels.md at curvesma.sh commit
6423c2728dc03d5f2d647c7f22adf1dac0643a66 were read for the standard GT-BSGS
and pinned pairing formulas. We rederive the fixed argument, joint
probability, full-coefficient keys, 64-bit accounting, record size and
new normalization. The old organizer ledger used a different word width
and is not adopted numerically. Relevant py_ecc attribution is retained
in pairing-kernels.md. All needed formulas and arguments are reproduced;
no previous acceptance, fetched theorem, actual collision, experiment or
completed production inversion is supplied as evidence.

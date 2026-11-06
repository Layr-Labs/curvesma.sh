# Complete Shanks decomposition for bls12-381-fixed-g2-invert-g1

## 1. Exact target and scope

This is a deterministic, complete baby-step/giant-step construction for
the assigned full-production target. It is a standard Shanks construction,
rederived below, not a new cryptanalytic result or a competitive record.
It was written in this account's independent workspace using only the
protected target/model and standard group facts, not another account's
candidate, proof, review or score. Algorithm correctness is unconditional
within the protected valid-input group contract; the huge local-link
storage machine assumed for computation is explicit in arithmetic.md.
No production-size search, giant table or local model was executed.

The multiplicative search group is the profile's exact reduced-pairing
GT, using base g=e(P0,Q) computed for this instance's supplied nonidentity
Q. Its target is z. Baby key B_i=g^i. Giant key Y_j=z*g^(-jM).
All 12 coefficients in the exact protected polynomial basis are compared.
After a match return the point P=[(i+jM) mod r]P0, not the scalar.
Verify the profile's exact reduced e(P,Q)=z and subgroup membership.
The supplied Q is never changed. If z=1, return infinity and check
the pairing identity. By nondegeneracy of the profile's prime-order
pairing, a nonidentity fixed argument induces an isomorphism from the
other order-r group onto GT, so g has order r. Equivalently the protected
challenge generation expresses z=g^a with a in 0..r-1; for these valid
inputs existence of a suffices for search completeness even without
using nondegeneracy to claim any new theorem. This claim uses precisely
the protected profile's algebraic group/pairing contract.

## 2. Finite algorithm and correctness

Set M=2^128. The actual subgroup order r in every selected profile is
less than 2^256=M^2. Cold-build and fully sort the M baby records with
indices i=0,...,M-1. Obtain the fixed stride [M]G by 128 doublings
(additive), or g^M by 128 squarings and its inverse by the known order
(multiplicative). Initialize Y_0 to the actual target. For j=0,...,M-1,
look up its exact canonical key in the sorted baby table; test any hit
by reducing the four-limb i+jM modulo the exact r and by the final
attack-side group verification. On no hit advance the giant by the
negative fixed stride, respectively multiply by g^(-M). Preprocessing,
all scans even in a failure case, sort, search, communication and final
verification are charged by arithmetic.md.

For each valid scalar x in 0,...,r-1, Euclidean division gives x=jM+i,
0<=i<M and 0<=j<M. At that giant index the target minus jM times the
base equals i times the base; multiplicatively it equals g^i. Hence
some lookup succeeds. A different matching index cannot create a false
output because equality implies the same scalar modulo r (and final
verification is performed). Even if the M babies exceed a group's order,
any matching baby index still yields a correct reduced scalar; no
assumption of unique integer representation inside the sorted table is
needed. Canonical infinity and all field coefficients are part of keys;
there is no sign-only, truncated-hash or near-equality collision.
Success is 1 for every valid target and therefore >=0.39 under the full
uniform distribution, with no conditioning on a favorable target.

## 3. Resources and modeling assumptions

arithmetic.md supplies the finite 64-bit-limb arithmetic, full affine
exceptional branches, constant-degree 128-level routed-word machine and cold-bootstrap protocol, 128 sort
passes and bounded binary lookups. Its complete instruction bound is
R<2^170 for the nonauxiliary target, and R(n)<2^170+(n+1)*2^52 for
an auxiliary input. Dividing by the exact protected normalization C,
including n=2^20 when applicable, gives work <2^162. Thus the sole
submitted log2-work upper bound is 162. Total memory <2^151 bytes
includes all nodes, copied code, buffers and scratch. These are analytical
upper bounds, not a promise that anyone owns such a storage machine.
The cold machine, local communication model and start-up costs are
explicit assumptions; unavailable practical hardware is not described
as an actual attack. The algorithm is classical, sequentially scheduled,
and charges all processor work rather than elapsed parallel time.

## 4. Fully specified public constants

The exact organizer profile, preserved verbatim as mathematical context:

```json
{
  "id": "bls12-381-fixed-g2-invert-g1",
  "curve": "BLS12-381",
  "parameter_x": "-d201000000010000",
  "problem": "fixed-argument-pairing-inversion-g1",
  "scope": "unrestricted-production",
  "field_modulus": "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab",
  "subgroup_order": "73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001",
  "curve_equation_G1": "y^2 = x^3 + 4 over Fp; G1=<P0> of order r",
  "cofactor_G1": "396c8c005555e1568c00aaab0000aaab",
  "generator_P0": {
    "x": "17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb",
    "y": "08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1"
  },
  "field_Fp2": "Fp[u]/(u^2+1); [c0,c1] means c0+c1*u",
  "curve_equation_G2": "y^2 = x^3 + 4*(1+u) over Fp2; G2=<Q0> of order r",
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
  "field_Fp12": "Fp[w]/(w^12-2*w^6+2); coefficient list [c0,...,c11] means sum(ci*w^i); u maps to w^6-1",
  "group_GT": "Nonzero z in Fp12 with z^r=1; identity [1,0,...,0]",
  "pairing": {
    "definition": "e(P,Q)=py_ecc.optimized_bls12_381.pairing(Q,P,final_exponentiate=True) at the pinned commit; the library takes G2 first. Q is supplied by the challenge and cannot be replaced. Raw Miller equality is not required.",
    "ate_loop_count": "d201000000010000",
    "final_exponent": "(p^12-1)/r",
    "reference": {
      "repository": "https://github.com/ethereum/py_ecc",
      "release": "v8.0.0",
      "commit": "04151f01f59f902ab932a51e0ca0ebce3883fc51",
      "source_sha256": {
        "py_ecc/optimized_bls12_381/optimized_curve.py": "f66da497e4efa555ce4b4364a06a154ff2038ed6c227b702b109fba298c1ad8e",
        "py_ecc/optimized_bls12_381/optimized_pairing.py": "0763415a7de64375dccfa62156d06505923ac6f8f3e1abab87126018d95abef9",
        "py_ecc/fields/field_properties.py": "9636251929dc93b9b9ae93b8aa0f91e14ed829248a0858b0fa881d5ff9105830",
        "py_ecc/fields/optimized_field_elements.py": "070159688a16c1346a31500f4f1a9d031cb7ce60c8fa1847ab3c534302a22c73"
      }
    }
  },
  "challenge": "Independently sample b uniformly from 1..r-1 and a uniformly from 0..r-1 by rejection sampling. Let Q=[b]Q0 and z=e(P0,Q)^a. Give only Q and z, not a, b, [a]P0 or sampling coins. Q is uniform over nonidentity G2 and z is independently uniform over GT, including identity. Return P in G1 such that e(P,Q)=z. No auxiliary powers are supplied.",
  "input": {
    "Q": "Nonidentity G2 point fixed for this challenge instance",
    "z": "GT coefficient list"
  },
  "output": {
    "P": "G1 point satisfying e(P,Q)=z for the supplied Q"
  },
  "acceptance": "Decode Q and z canonically; require Q on the G2 curve, Q!=infinity and [r]Q=infinity, and require z!=0 and z^r=1. Decode the returned P canonically; require P on the G1 curve and [r]P=infinity. Accept exactly when e(P,Q)=z, comparing all 12 coefficients. P=infinity is valid exactly when z=1. Returning a scalar or substituting another Q is not a solution.",
  "encoding": "UTF-8 JSON with exactly the stated input or output keys. Each Fp coefficient is 96 lowercase hex characters, no prefix, representing an integer less than p. A G1 point is null for infinity or {x,y} with Fp coefficients. A G2 point is {x,y} with each coordinate [c0,c1] in Fp2 (input infinity is forbidden). A GT value is exactly 12 coefficients in increasing polynomial degree. Reject duplicate keys, noncanonical encodings, off-curve points and points outside the stated subgroups.",
  "source": "https://github.com/ethereum/py_ecc/tree/04151f01f59f902ab932a51e0ca0ebce3883fc51/py_ecc/optimized_bls12_381"
}
```

## 5. Evidence and provenance

This new candidate uses no previous participant's accepted claim,
judgment, stored result or numerical score. The generic decomposition is
standard: x=jM+i is proved here, not imported from an accepted-result
service. Any shared text between this account's six new curve drafts
is the same newly authored resource certificate and is disclosed as
common preparation, not six unrelated discoveries. No scientific novelty,
prize, practical break or superiority to a current incumbent is claimed.
Zero experiments are declared; independent advisory critique and official
cloud review must be recorded separately. For pairing directions the exact protected pairing expressions and
source-bound denominator proof are reproduced in arithmetic.md. Final
advisory and official qualification remain distinct; no acceptance is
inferred from this draft status.

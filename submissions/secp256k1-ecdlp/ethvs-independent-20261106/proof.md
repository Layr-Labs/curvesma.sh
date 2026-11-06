# Complete Shanks decomposition for secp256k1-ecdlp

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

The additive search group is the exact prime-order subgroup generated
by the profile's P; its target is the supplied Q. Baby key B_i=[i]P.
Giant key Y_j=Q-[jM]P. Return (i+jM) mod r after a full key
match, and explicitly verify [x]P=Q. Infinity is an ordinary
key at i=0; an identity target is handled and returns scalar zero.
The actual curve coefficients, generator, field prime, subgroup order
and encoding are the full values below, not a toy group.

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
  "id": "secp256k1-ecdlp",
  "curve": "secp256k1",
  "problem": "ecdlp",
  "scope": "unrestricted-production",
  "field_modulus": "fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f",
  "subgroup_order": "fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141",
  "curve_equation": "y^2 = x^3 + 7 over Fp",
  "cofactor": 1,
  "generator_P": {
    "x": "79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798",
    "y": "483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8"
  },
  "challenge": "Sample x uniformly from 0..r-1 by rejection sampling; give only Q=[x]P. The scalar and sampling coins stay private.",
  "input": {
    "Q": "secp256k1 point"
  },
  "output": {
    "x": "32-byte scalar"
  },
  "acceptance": "Decode Q canonically, require Q in <P>; decode 0<=x<r; accept exactly when [x]P=Q.",
  "encoding": "UTF-8 JSON with exactly the stated keys. Scalar: 64 lowercase hex characters, no prefix. Point: null for infinity, otherwise {x,y}, each 64 lowercase hex characters representing a big-endian integer less than p. Reject duplicate keys, noncanonical values and off-curve points.",
  "source": "https://www.secg.org/sec2-v2.pdf#page=13"
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
cloud review must be recorded separately. For pairing directions the exact
protected pairing implementation still needs source-bound arithmetic
coverage before the draft can be considered ready; no success or acceptance
is inferred from this draft status.

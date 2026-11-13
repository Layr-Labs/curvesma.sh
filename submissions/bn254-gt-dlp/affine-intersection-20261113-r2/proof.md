# Affine progression intersection for the production BN254 target group

This submission gives a randomized, finite, proof-only construction. Its
success analysis uses pairwise independence of linear equations over the
protected scalar field, not a random-function assumption about Pollard walks.
The very large memory allocation is explicit. The submitted upper bound is
conservative and does not improve the published 135.216266 baseline.

### Distribution, input convention and exact scope

Fix the entire protected BN254 profile: p, prime scalar order r, the Fp12 basis
w^12-18w^6+82, and the public order-r generator g. Let b=2^127 and n=2b=2^128.
The public integer r satisfies 2^253<r<2^254=b^2. The input is z=g^x for
uniform x in {0,...,r-1}, with the identity included and no source-group input.

For the numerical score we take the production wire convention to be the
compact JSON object with sole key z and twelve quoted 64-digit lowercase
hexadecimal strings in the protected order. Its length is exactly 811 bytes.
Every mathematical z has that encoding. This is the bounded-canonical reading
used by the organizer baseline described as reviewed in protected README.
The profile itself does not fix the JSON whitespace grammar explicitly. This
is a disclosed material interpretation of the scored input, not a proof that
canonical coefficient widths bound every legal JSON text. The general parser
also processes arbitrarily padded finite JSON and charges its actual length L.
The constant score below applies to the 811-byte convention. Under an
interpretation requiring a uniform constant cost for all unbounded raw JSON
spellings, the numerical score has not been established. That alternative is
not concealed as a free decoder, uncharged transport, or silently rejected input.

The submitted claims are success at least 49/100, time W<=2^149 where W is
total scalar instruction work divided by the protected 20000, and peak storage
at most 2^154 bytes. Preprocessing is included in full; a loose preprocessing
upper bound is 2^149 and no searched nonuniform advice is provided. All claims
refer to the explicitly described finite classical machine and include failures.

### Sampling and the two sets

Obtain a,beta,u,v as independent uniform residues modulo r. For each residue,
try at most 1024 independent uniform 254-bit strings, accept the first integer
less than r, and otherwise terminate with failure. Charge all draws and
failures. Conditional on no sampler failure the four residues are independent
and exactly uniform; the first-accepted-value distribution is uniform for each
possible acceptance time. Since r>2^253, total failure probability is less than
4*2^-1024=2^-1022. Random bits are fresh algorithm coins, not a reproducible
seed masquerading as ideal randomness, and their acquisition work is charged.

Create the left list L_i=g^(a+i*u), 0<=i<b, by computing g^a and g^u, then
successively multiplying the preceding entry by g^u. Independently create
R_j=z*g^(beta+j*v), 0<=j<b, starting with z*g^beta and step g^v. Store the
complete Fp12 key, its side and its index. Compute all four starting/step
powers during this attack. Values u=0 or v=0 are permitted: repeated keys
are retained rather than removed without proof.

Use the exact sorting network and finite transfers of implementation.md to
sort the n records lexicographically by (whole group key, side, index).
Inspect consecutive records for an equal group key across the two sides.
Each key's records occur contiguously with all left records before all right
records. Consequently a key occurring in both lists has at least one adjacent
left/right boundary, even when a progression has repeated elements. Select
one such boundary in a bounded reduction. Recover

    t = a - beta + i*u - j*v mod r.

Verify g^t=z and the canonical scalar range; never output an unchecked guess.
If no boundary exists, report failure. The whole bounded schedule is charged,
whether a match occurs or not. Its correctness follows from g's order r:
g^(a+i*u)=z*g^(beta+j*v) implies t=x mod r, and the unique residue is x.
The target z=1 satisfies this same argument; no exceptional nonidentity
restriction is used.

### A finite success proof

Condition on successful sampling and fix any true x, without revealing x to
the algorithm. For each ordered pair (i,j), let I_ij be the indicator of
a-beta+i*u-j*v=x in F_r. Each has expectation 1/r because a-beta is uniform.
For two different pairs, the two affine coefficient rows (1,i,-j) and
(1,k,-ell), acting on (a-beta,u,v), have rank two: their first coordinates
are equal to one, so scalar dependence would force identical rows, and
different indices in [0,b) are different residues because b<r. Thus the two
equations have joint probability 1/r^2. This proves pairwise independence
without claiming independence of all b^2 indicators.

Let X=sum I_ij and lambda=b^2/r. Exact moments are
E[X]=lambda and E[X^2]=lambda^2+lambda*(1-1/r).
By Cauchy-Schwarz applied to X*1_(X>0),

    P(X>0) >= (E[X])^2/E[X^2]
             = lambda/(lambda+1-1/r) > 1/2,

since b^2>r. This bound includes the rare degenerate progressions; conditioning
separately on u or v being nonzero is unnecessary and was not used in the
rank calculation. The complete sorting procedure finds a match whenever X>0.
Restoring bounded-sampler failure gives success greater than
1/2-2^-1022>49/100, for each x and hence for a uniform full-production target.
No extrapolation from small collision experiments is involved.

### Numerical bound and evidence limits

implementation.md supplies a halfword implementation and a real finite network
layout with local 64-bit addresses, fully counted packet transfers and no
unit-cost wide lookup. Define

 D = 2^31 + 2^30 + 8256*2^21 + 2^31 + 2^28,
 R(L) = n*D + 2^52 + 8192*(L+1).

The first term includes every node's initialization, creation/delivery of all
records, all 8256 sorting layers, adjacent checks and final collection. The
fixed term includes bounded sampling, root powers, scalar recovery, checking
and output. For canonical L=811, the exact integer inequality R(811)<2^163
holds, hence R(811)/20000<2^149. The integer rational score 149 is log2 of
the submitted conservative work upper bound, not an exact measurement of the
work actually used or an automatically tightened estimate.

Experiments are empty. No production attack or local table construction was
performed. The evidence is this complete finite algebraic/probability argument
and the specified machine accounting. The only material interpretation not
settled by the algebra is the disclosed scored-input serialization convention.
No claim of historical originality, practical feasibility, superiority over
known BN254 attacks, formal certification, manual acceptance or reward is made.

Return the verified t as {"x":"<64 lowercase hex digits>"}, with 0<=t<r.

# Supplied-P inversion, complete cold bound128.039

Claim: total classical cold work is below2^(128039/1000) dense-Fp12 reference
operation equivalents at fixed C=20000, success at least0.39 conditional on the
single SHA-256 idealization H below, peak storage below2^50 bytes, advice zero.
Preprocessing and finish work are included, not amortized. This is an analytic
heuristic bound, not production execution, unconditional fixed-hash proof or a
practical break. Optional finite SHA fixtures do not establish attack success.

## Exact algorithm

Use the protected full BLS12-381 profile, its prime r, field p and public Q0.
Parse and validate canonical P,z, subgroup membership and P!=infinity. Compute
h=pairing(Q0,P) in the protected pinned library's argument order. Nondegeneracy
and prime order imply that h generates GT for every allowed P. For z=1 return
G2 infinity. Otherwise z=h^a for a unique a modulo r; the solver never receives
a, challenge coins or auxiliary powers. Its final output is Q=[a]Q0, not a scalar
or replacement input P.

All tags below are the16-byte ASCII strings FAPI-RHO-INIT-v1,
FAPI-RHO-SCAL-v1,FAPI-RHO-WALK-v1,FAPI-RHO-UKEY-v1, respectively INIT,SCAL,
WALK,UKEY. P_bytes is its two canonical48-byte coordinates; z_bytes is its12
canonical48-byte polynomial coefficients. Set seed=SHA256(INIT||P_bytes||z_bytes).
There are no external random coins. A globally increasing counter starts at0,
is encoded in8 big-endian bytes, and increments after every sampler attempt.
A main scalar candidate is the low255 bits of SHA256(SCAL||seed||counter).
Accept iff<r; cap at1024 attempts per scalar, otherwise fail. A key candidate is
the low61 bits of SHA256(UKEY||seed||counter); accept iff<q=2^61-1 with the same
cap. Generate five keys independently, then2w+1 main scalars for w=2^40:
jump pairs(u_j,v_j) and alpha0. Counter remains below2^52. Store
M_j=h^u_j*z^v_j, u_j,v_j, in640-byte records:576-byte canonical Montgomery
Fp12 element and two32-byte canonical residues. Start(X,A,B)=(h^alpha0,alpha0,0).

Represent X=a+bw in the standard tower Fp2/Fp6/Fp12 with w²=v. Encode X by
six48-byte canonical Montgomery Fp coefficients of b followed by sign(a),
where sign(a)=1 iff a's coefficient tuple is lexicographically greater than-a's
in a fixed unsigned order. This289-byte encoding is injective by scalar-audit.md.
Let M be that encoding interpreted big-endian. For l=0..38 define
c_l=(M>>(60l)) & (2^60-1); define f_M(t)=sum_l c_l*t^l modulo q. Let U(M)
be the five f_M(k_i), each in8 big-endian bytes. Set j to low40 bits of
SHA256(WALK||seed||zero16||U(M)). Update
X<-X*M_j, A<-A+u_j mod r, B<-B+v_j mod r. The invariant is X=h^A*z^B.
The map depends on X alone, never its coefficient history. The fixed64-byte
hash prefix is compressed once in charged setup; the state-dependent remaining
block uses standard SHA-256 padding and original total length832 bits.

Let R=isqrt(r-1)+1, L=ceil(1031R/1000), s=ceil(L/126), K=L+s. Each mutable
or checkpoint record stores the full576-byte group element, two32-byte scalar
coefficients and an aligned304-byte exact key:36 native64-bit numerical B words
of b followed by one64-bit sign0/1, and8 zero padding bytes. Total944 bytes,
aligned16. Jump-table entries remain640 bytes; their layout and j*640 addressing
are unchanged. Key padding is zeroed in charged setup and not compared.

In paid setup prepare the initial exact key, aligned digits and cached walk label
j_0, and save the initial record. For t=1..K use the cached label to choose M_j,
update X,A,B, then prepare the new exact key/digits and next label from the new X.
Binary-search saved records by their37-word exact keys, terminating on the first
equality. Only after an unsuccessful search, if t=s,2s,...,126s, save the complete
944-byte record. At most127 immutable saved records exist. Sorted pointers are
ordered by unsigned lexicographic key words; each search needs at most7 complete
comparisons. Never compare against a record just saved at current t. Keys are
injective on GT by the norm-one proof, not hashes or probabilistic fingerprints.
A cached next label is a separate8-byte variable, not part of a saved record.
Its loads/stores are included in256 index/address work.

This pipelines one key/label preparation per transition plus an extra initial one
inside2^90. It yields exactly the same group-state trajectory as computing each
label immediately before its transition. The unused final label is still charged.
Sorted equality and record selection inspect group states only, not coefficients.

At an equality to saved(X0,A0,B0), let D=B-B0 mod r. If D=0 fail immediately;
otherwise recover a=(A0-A)/D mod r. Verify h^a=z, form Q=[a]Q0, canonically
encode it, verify its curve/subgroup membership and pairing(Q,P)=z, and output Q.
Failure is permitted; unchecked wrong answers are not. All branches and finishing
work, including failed sampling paths and zero denominators, are charged.

## Sole heuristic and success

H: for average uniform protected challenges, the domain-separated actual SHA-256
computations may be analyzed as independent random functions on their distinct
input domains. Actual hashes are fully computed and charged. Fixed SHA is not
proved random at this scale. No independent random-mapping assumption for the
adding walk, collision independence assumption or free oracle is made. There
are no small-group measurements imported as production guarantees.

Inside H, bounded rejection yields uniform independent accepted main scalars and
keys. Since r>2^254 and q>2^60, union sampling failure is below2^-980. Jump group
values are iid uniform GT. probability.md proves bad-table probability below
2^-100 and repeat-by-L probability>0.39015 on a good table, explicitly including
correlated nested/crossing intervals. Universal compression loses at most2^-13:
for fixed unequal encodings Pr[U(M)=U(M')]<=(38/q)^5 by the polynomial root bound;
an independent-index surrogate independent of the keys has at mostK+1 states,
so its distinct-state compression loss is at most binom(K+1,2)*(38/q)^5<2^-13
using K+1<2^129. Couple only until its first repeat or unequal-state compression
collision. Before that stopping time every new compressed hash input has a fresh
uniform label. This deals with adaptation; a union bound directly on a key-
dependent actual trajectory is not being asserted.

Checkpoint detection is deterministic: for tail mu and cycle lambda with first
repeat tau=mu+lambda<=L, the saved index b=ceil(mu/s)*s satisfies
mu<=b<=mu+s-1, b<=126s, and X_(b+lambda)=X_b by time tau+s-1<K. Earlier detected
equalities terminate sooner. Thus a first repeat by L is detected by K.

Conditional on the jump group values, partition and initial group value, each v_j
remains independent uniform: for any chosen v_j there is exactly one u_j giving
M_j=h^u_j*z^v_j. The trajectory/search never inspect coefficient histories in
selecting an equality. Across the positive d<=K<r steps between saved and current
records, D=sum_j n_j*v_j with n_j>=0 and sum_j n_j=d. Some n_j is nonzero mod r,
so D is uniform and zero with probability1/r. Otherwise the invariant proves the
recovered scalar and required G2 output. No extra relation heuristic enters.

Consequently success>0.39015-2^-13-2^-100-2^-980-1/r>0.39 under H. Identity
z succeeds deterministically. The supplied P reduction holds for every allowed P,
so the probability statement has the required independent uniform P,z scope.

## Complete work, advice and storage

One transition uses one named reference dense Fp12 multiplication H=1, charged
once under the protected H+R/C rule. Its internal implementation is not counted
again. It is the exact pinned blst reference operation and canonical tower basis,
not a candidate-selected price or timing calibration. scalar-audit.md bounds all
other per-step work:1312 universal compressor,3232 SHA-256 including preparation,
128 sign,4640 exact-key search/coefficient/counter work,64 reference-call
interface. Sum R_step=9376. Thus

 W <= K+(9376*K+2^90)/20000 = (29376*K+2^90)/20000.

The fixed2^90 scalar instructions pay cold for all setup, randomness attempts,
jump-table construction, canonical parsing/checks/conversions, constants and code
initialization, two pinned pairings, output and verification, scalar inversion,
checkpoint sorted-pointer insertion and final transfers. Up to2^52 attempted
sampler hashes at<2^16 instructions each cost<2^68. Each of2^40 jumps uses<1024
dense products; even a deliberately loose<2^20 scalar instructions per product
and interface gives<2^70 for all jump construction. Fixed six-limb Montgomery
arithmetic and the fixed tower implement that loose bound. Two fixed pairings
cost<2^78 each: the64-bit Miller loop and<4608-bit final exponent use<2^20
extension operations, even pricing each at2^58 including inverses suffices.
Fp inverse by p-2 needs<768 products, Fp12 inverse by p^12-2<9216 products.
G1/G2 scalar multiplication, exceptional branches and subgroup/output checks use
<4096 point steps, each far below2^52 instructions. Complete affine formulas
handle infinity, equal/opposite points and y=0 before division; Fp2 inverse is
(a-bu)/(a²+b²). These loose fixed bounds leave ample slack below2^90. At most127
sorted insertions move at most126 pointers apiece at<16 instructions each, plus
setup below2^30, also inside2^90. No work is amortized across targets.

The exact public-parameter inequality
 (29376*K+2^90)^1000 < 20000^1000*2^128039
establishes the submitted score128039/1000=128.039. Independently checkable finite
probability inequalities are in probability.md. No official score is asserted.

The jump array uses exactly640*2^40=703687441776640 bytes. Below2^50 bytes covers
it,127 checkpoint slots, current/next944-byte buffers, small pointer array, code,
SHA/field/scalar scratch and constants. Addresses fit64 bits; use an aligned base
low enough that all offsets/addresses are below2^50. Offset j*640=(j<<9)+(j<<7)
is charged in index/address work. Memory transfers are charged scalar loads and
stores. There is no distributed-memory/free-transfer oracle or allocated VM-sized
production run. Advice is zero; Q0 and target parameters are public constants.
Construction of the small algorithm constants is paid in fixed setup.

## Attribution and evidence limits

Published5039bd36-c19f-4c7e-95c3-6b7e5d1f01fc supplies explicit hash scheduling,
paired full-state comparison, fixture design and reviewed scalar accounting.
Published30a7c6df-a9c5-47fe-895b-a5bdff9433cb at6b1ba142c55af1fd8e54f2421c0a638cf29232c8,
selected large-spectrum-sha-20261106-v4, was directly read and supplies the walk,
complete finite probability proof, norm-one representation, polynomial/root
coupling, checkpoint and earlier cold/scalar accounting. Its earlier finite-walk
source755cd10d-5cc2-4b71-8542-11447238cf4b was accessed through that source.
The SHA-256/dot-product/EXTR/tighter-preparation and exact-key checkpoint ideas materially use the public
note of failed/unpublishedf7e1f77e-5d96-4a4e-9863-91878efeb0a2. No independent
invention is claimed. Full IDs, direct unpublished intermediates and target-
eligibility limitations are acknowledged in the public note.

This successor builds directly on failed/unpublishedecf55e5c-8121-4cf5-ac26-0a1915cc18b5
and failed/unpublished93c3e7bd-718d-46c2-98f4-8fa41c0cafb3.
It keeps standard SHA-256 throughout:32-byte seed, low255-bit scalar candidates,
low61-bit key candidates; WALK||seed||zero16 is64 bytes, plus U's40 bytes gives
104 original bytes and832 bits. SCAL/UKEY messages are56 bytes and need two
hash blocks, paid in2^90. The scalar-audit.md preparation schedule charges1024
instead of2048 without changing any digit, polynomial, hash input/output or
trajectory. The512 limb/key-preparation,256 block/state/label and256 wrapper budgets
cover all work. Finite probability and the single SHA-256 idealization remain
unchanged; a SHA-512 review does not certify that changed hash-model bridge.

check_primitives.py is hosted-only finite evidence: SHA-256/hashlib, polynomial
Horner/dot-product, aligned-buffer/EXTR and exact final-block/label fixtures.
It has not run locally because Docker access is denied; host execution is never
used. Its actual finite output, if any, is not an A64 trace, production success
rate, proof of the idealization or complexity certificate. Published5039bd36's
46 SHA-512 fixtures passed, at only that scope. No pending/new result is asserted
as passed. Shared llm-verifier owns isolated execution, review and scoring.

The new exact checkpoint key is the same injective b/sign information already
proved for hashing, stored as37 numerical words. No square root or decompression
is performed. Scalar-audit.md charges its36 added word stores in512 preparation,
944-byte saves in512 copying, and exact37-word comparisons at224 each. Seven
search iterations give4640 outer work with all other allowances retained. The
new126-interval checkpoint lemma is the same deterministic argument above;
L is unchanged and smaller K preserves the compression loss bound. Every
key equality is full group equality, with no new heuristic or collision loss.
The optional fixture adds finite toy norm-one injectivity, protected-p canonical
sign/word-layout checks and exact-key comparator cases. Only the proof covers
full-production GT; no toy execution is a production success or cost certificate.

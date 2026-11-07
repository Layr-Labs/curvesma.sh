# Supplied-P inversion, complete cold bound127.994

Claim: total classical cold work is below2^(127994/1000) dense-Fp12 reference
operation equivalents at fixed C=20000, success at least0.39 conditional on the
single BLAKE2s-256 idealization H below, peak storage below2^50 bytes, advice zero.
Preprocessing and finish work are included, not amortized. This is an analytic
heuristic bound, not production execution, unconditional fixed-hash proof or a
practical break. Optional finite hash fixtures do not establish attack success.

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
canonical48-byte polynomial coefficients. Set seed=B2S(INIT||P_bytes||z_bytes).
There are no external random coins. A globally increasing counter starts at0,
is encoded in8 big-endian bytes, and increments after every sampler attempt.
A main scalar candidate is the low255 bits of B2S(SCAL||seed||counter).
Accept iff<r; cap at1024 attempts per scalar, otherwise fail. A key candidate is
the low61 bits of B2S(UKEY||seed||counter); accept iff<q=2^61-1 with the same
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
B2S(WALK||seed||zero16||U(M)). Update
X<-X*M_j, A<-A+u_j mod r, B<-B+v_j mod r. The invariant is X=h^A*z^B.
The map depends on X alone, never its coefficient history. The fixed64-byte
hash prefix is compressed once in charged setup; the state-dependent remaining
block is U(M) followed by24 zero bytes, with cumulative byte counter104 and
final-block flag true. The prefix compression uses counter64 and final false.

Let R=isqrt(r-1)+1, L=ceil(1031R/1000), s=ceil(L/126), K=L+s. Each mutable
or checkpoint record stores the full576-byte group element, two32-byte scalar
coefficients and an aligned304-byte exact key:36 native64-bit numerical B words
of b followed by one64-bit sign0/1, and8 zero padding bytes. Total944 bytes,
aligned16. Jump-table entries remain640 bytes; their layout and j*640 addressing
are unchanged. Key padding is zeroed in charged setup and not compared.

In paid setup prepare the initial exact key, aligned digits and cached walk label
j_0, and save the initial record. For t=1..K use the cached label to choose M_j,
update X,A,B, then prepare the new exact key/digits and next label from the new X.
Use lower_bound over saved records sorted by the unsigned integer whose little-indexed limbs are the37 native key words (word0 least significant). At each iteration compare saved[mid] < current with a full borrow chain, and narrow the half-open interval [lo,hi). With n<=127, at most7 ordering comparisons locate the insertion point. If that point is below n, compare all37 words there for equality once. Terminate on equality. Otherwise, only when t=s,2s,...,126s, insert the complete944-byte record at the insertion point. Never compare to a record just saved at current t. Keys remain injective; ordering is arbitrary but deterministic and uses only group states.
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

H: for average uniform protected challenges, the domain-separated actual BLAKE2s-256
computations may be analyzed as independent random functions on their distinct
input domains. Actual hashes are fully computed and charged. Fixed BLAKE2s is not
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
other per-step work:1312 universal compressor,2432 BLAKE2s-256 including preparation,
128 sign,4512 exact-key search/coefficient/counter work,64 reference-call
interface. Sum R_step=8448. Thus

 W <= K+(8448*K+2^90)/20000 = (28448*K+2^90)/20000.

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
 (28448*K+2^90)^1000 < 20000^1000*2^127994
establishes the submitted score127994/1000=127.994. Independently checkable finite
probability inequalities are in probability.md. No official score is asserted.

The jump array uses exactly640*2^40=703687441776640 bytes. Below2^50 bytes covers
it,127 checkpoint slots, current/next944-byte buffers, small pointer array, code,
Hash/field/scalar scratch and constants. Addresses fit64 bits; use an aligned base
low enough that all offsets/addresses are below2^50. Offset j*640=(j<<9)+(j<<7)
is charged in index/address work. Memory transfers are charged scalar loads and
stores. There is no distributed-memory/free-transfer oracle or allocated VM-sized
production run. Advice is zero; Q0 and target parameters are public constants.
Construction of the small algorithm constants is paid in fixed setup.

## Attribution and evidence limits

This revision directly derives from submission23d3d4cd-768d-4cac-a9a5-872fc010580d, evaluated commit acf0f4bfe27777d98b306b7a20340fac646d9867, folder blake2s-full-20261106-v10. Its complete walk, probability proof, representation, scalar audit and cold accounting are reused here. Its earlier attribution to5039bd36-c19f-4c7e-95c3-6b7e5d1f01fc and30a7c6df-a9c5-47fe-895b-a5bdff9433cb is preserved as transitive provenance; those archives were not directly used for this revision. No independent invention of those components is claimed. The new contribution is the borrow-chain lower_bound search, a single final exact equality check and its revised instruction bound. No experimental result is imported or asserted.

B2S is standard unkeyed sequential BLAKE2s-256, full10 rounds, fanout1/depth1, other parameters zero, last-node false. Initial h0 is IV0 XOR0x01010020; remaining h words are IV. Digest bytes are interpreted big-endian for low-bit extraction, while seed is raw digest bytes. The specification is RFC7693 sections2-3: https://www.rfc-editor.org/info/rfc7693/ . The64-byte prefix is compressed nonfinally at counter64; the final40-byte payload is zero-padded to64 bytes with counter104 and final flag true. SCAL/UKEY messages are56-byte final blocks. No reduced-round hash, free hash oracle or production execution is claimed.

H remains a material unproved assumption at roughly2^128 queries. The RFC specifies the function but does not prove this random-function bridge. This proof-only candidate declares no executable experiments. Local intake and artifact checks are not scientific qualification; shared llm-verifier owns review and scoring. The full key is compared exactly, and the ordering change introduces no fingerprint collision loss. Carry chains give ordering only; equality is deliberately checked separately.

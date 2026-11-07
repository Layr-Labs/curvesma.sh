# Supplied-P inversion, complete cold bound127.903

Claim: total classical cold work is below2^(127903/1000) dense-Fp12 reference
operation equivalents at fixed C=20000, success at least0.39 conditional on the
single BLAKE3 idealization H below, peak storage below2^50 bytes, advice zero.
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
canonical48-byte polynomial coefficients. Set seed=B3(INIT||P_bytes||z_bytes).
There are no external random coins. A globally increasing counter starts at0,
is encoded in8 big-endian bytes, and increments after every sampler attempt.
A main scalar candidate is the low255 bits of B3(SCAL||seed||counter).
Accept iff<r; cap at1024 attempts per scalar, otherwise fail. Generate2w+1 accepted main scalars for w=2^40:
jump pairs(u_j,v_j) and alpha0. Counter remains below2^52. Store
M_j=h^u_j*z^v_j, u_j,v_j, in640-byte records:576-byte canonical Montgomery
Fp12 element and two32-byte canonical residues. Start(X,A,B)=(h^alpha0,alpha0,0).

Represent X=a+bw in the standard tower Fp2/Fp6/Fp12 with w²=v. Encode X by
six48-byte canonical Montgomery Fp coefficients of b followed by sign(a),
where sign(a)=1 iff a's coefficient tuple is lexicographically greater than-a's
in a fixed unsigned order. This289-byte encoding is injective by scalar-audit.md.
The exact native key is36 canonical numerical64-bit b words followed by sign as a64-bit0/1, as described below. Let c_0,...,c_71 be the low and high32-bit halves of b words0..35 in that order, and c_72=sign. Precompute independent320-bit tables T_i for i<72 on the32-bit alphabet, and T_72 on{0,1}. The fingerprint U(X) is the bitwise XOR of T_i[c_i] over all73 positions, stored as five64-bit words. Serialize these five words big-endian as40 bytes. Set j to the low40 bits of B3(WALK||seed||zero16||U(X)).

Each table entry is constructed cold, not supplied as advice. Define TABLE as the16-byte ASCII string FAPI-RHO-TABL-v1. Hash TABLE||seed32||i_1||c_4||b_1, with i_1 one unsigned byte, c_4 four big-endian bytes and b_1 byte0 or1, using standard unkeyed BLAKE3 with32-byte output. This is54 bytes, one block of length54, counter0, flags11 (CHUNK_START|CHUNK_END|ROOT). Entry words0..3 are the four big-endian numerical64-bit digest words for b=0; word4 is the first big-endian64-bit word of the b=1 digest. Store five native words. Construct every entry, even ones never read. TABLE is distinct from INIT,SCAL,WALK. H below now explicitly covers the additional TABLE domain. Inside H all table values are independent uniform320-bit strings, independent of jumps and WALK labels. This is an enlarged disclosed use of the new unproved BLAKE3 random-function idealization, not a new fixed-hash theorem or free randomness.

Update
X<-X*M_j, A<-A+u_j mod r, B<-B+v_j mod r. The invariant is X=h^A*z^B.
The map depends on X alone, never its coefficient history. The fixed64-byte hash prefix is compressed in charged setup with counter0, block length64 and flags1 (CHUNK_START). The state-dependent remaining block is U(X) followed by24 zero bytes, with counter0, block length40 and flags10 (CHUNK_END|ROOT). The cached prefix is a chaining value, not a finalized digest.

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

H: for average uniform protected challenges, the domain-separated actual BLAKE3
computations may be analyzed as independent random functions on their distinct
input domains. Actual hashes are fully computed and charged. Fixed BLAKE3 is not
proved random at this scale. No independent random-mapping assumption for the
adding walk, collision independence assumption or free oracle is made. There
are no small-group measurements imported as production guarantees.

Inside H bounded scalar rejection yields uniform independent jumps and start; the capped sampling failure remains below2^-980. There are fewer sampled scalars than in the predecessor, so its conservative bound still applies. probability.md retains the finite good-table and correlated-interval argument, yielding repeat-by-L probability>0.39015 except for bad-table loss2^-100. For every fixed unequal exact key pair, one digit differs. Condition on every table entry except T_i[c_i] at such a differing position; the XOR difference of fingerprints is uniform320-bit, so collision probability is exactly2^-320. This fact does not assert higher-order independence.

Use the independent-index surrogate, whose group states depend on jumps and independent uniform labels but not these tabulation tables. Among at mostK+1<2^129 states, union loss is less than2^257*2^-320=2^-63, hence below the retained conservative2^-13. Couple only until first group repeat or distinct-state fingerprint collision. Before then, each new WALK input has a fresh uniform label under H. This deals with adaptation without claiming the actual table-dependent trajectory is independent of tables. The checkpoint and relation extraction arguments below retain the same assumptions and numerical losses.

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
other per-step work:1024 simple-tabulation compressor,1536 BLAKE3 including preparation,
128 sign,3968 exact-key search/coefficient/counter work,64 reference-call
interface. Sum R_step=6720. Thus

 W <= K+(6720*K+2^90)/20000 = (26720*K+2^90)/20000.

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
 (26720*K+2^90)^1000 < 20000^1000*2^127903
establishes the submitted score127903/1000=127.903. Independently checkable finite
probability inequalities are in probability.md. No official score is asserted.

The jump array uses640*2^40=703687441776640 bytes. The tabulation entries use (72*2^32+2)*40=12369505812560 bytes. Their sum is716056947589200 bytes, below2^50=1125899906842624 with more than400 trillion bytes of slack for code, tables of pointers, checkpoints, stack and scratch. All addresses fit64 bits with a sufficiently low aligned base. Each fingerprint table entry is40 bytes; x40 addressing is explicitly charged. Constructing every entry uses618475290628 one-block hashes, each conservatively below2^16 instructions, below2^56 total; byte encodings, pointers, writes and construction loop fit comfortably below2^60. This is contained, along with all inherited cold costs, in fixed2^90 instructions. Storage is an analytical algorithm resource bound; the VM does not allocate the production arrays. No free table, external advice, multi-target amortization or uncharged transfer is used.

## Provenance and limitations

Direct published source23d3d4cd-768d-4cac-a9a5-872fc010580d, evaluated commit acf0f4bfe27777d98b306b7a20340fac646d9867, supplies the walk, finite probability proof, representation, full-round BLAKE2s model and accounting. This local source chain also materially uses this session's d4df1ac5-54c0-4d48-84ac-afd017e602e2 (borrow-chain ordering),3c59f3c3-7d9f-4985-b632-e67edb13702b (explicit lookup leaf), and53c1fb75-4153-43de-8e84-8296aa93353d (simple tabulation). They were unpublished at preparation; formal relationships will be added when publication permits. No independence is claimed. Earlier sources cited by23d3d4cd are transitive provenance, not directly inspected archives here.

The direct local predecessor removed unused polynomial digit preparation from the hash preparation charge. This revision changes the full hash primitive to standard BLAKE3, with a separately disclosed new H bridge. The actual key stores and block/state/label operations now have explicit512-instruction preparation accounting in scalar-audit.md. The tabulation construction rule and domain tags are inherited from the local53c1fb75 package; all actual values change with BLAKE3. This is analytic implementation evidence, not measured runtime or a production execution. No executable experiments are declared; candidate files never execute on the host or with credentials. Shared llm-verifier owns review and scoring. H is unproved at this scale; the RFC7693 specification defines BLAKE2s but does not establish the random-function bridge. No scientific qualification is inferred from structural intake or a previous review.

B3 denotes standard unkeyed BLAKE3 with32 output bytes, all7 rounds. Official reference implementation pinned at f55849f89cd85083c1c9daa2c5de20766f309c00:
https://github.com/BLAKE3-team/BLAKE3/blob/f55849f89cd85083c1c9daa2c5de20766f309c00/reference_impl/reference_impl.rs . The authors' specification is https://github.com/BLAKE3-team/BLAKE3-specs/blob/master/blake3.pdf .

This changes every seed, scalar sample, TABLE value, WALK label and actual trajectory from BLAKE2s. H is a new explicit bridge for standard BLAKE3, not inherited qualification of that function. Within this new model the uniform accepted scalars, independent320-bit table entries and fresh40-bit labels have exactly the distributions used in probability.md. No extra independent-mapping or relation heuristic is asserted. All real hash calls are computed and charged.

All inputs fit one1024-byte chunk, so its counter is0 throughout. SCAL56 and TABLE54 have flags11; INIT688 uses11 blocks, first flags1, intermediate flags0, final48-byte block flags10. WALK uses the cached first64-byte chaining value and the final40-byte block with flags10. Unkeyed initial h is the standard IV without BLAKE2 parameter XOR. Digest words serialize little-endian; the digest byte string is interpreted big-endian only for the specified bit extraction. Message permutation, root flags and counters follow standard BLAKE3. No candidate-specified reduced-round primitive is used.

The authors motivate standard7-round BLAKE3 by cryptanalysis and design arguments; that is support for considering this disclosed model, not a proof that the exact public-prefix path is random at2^128 queries. Hao et al.'s primary abstract at https://eprint.iacr.org/2023/299 reports a full7-round keyed-permutation boomerang attack at complexity2^180. It is an inference that this does not itself provide a failure case for this standard unkeyed, proper-root-flag hash path at the claimed work scale. We do not certify absence of all hash weaknesses. H remains unproved and material; prior BLAKE2s reviews or finite fixtures do not establish it.

This revision materially uses the local direct-key-preparation package submitted as9bd9a1fc-a097-48e0-8b29-f0e80a8bed99, which failed in hosted credentialed review without a scientific verdict. Its relationship will be declared when publication permits; its use is not an independence claim. The published23d3d4cd dependency remains direct for the complete mathematical walk and accounting framework. No executable experiment or production execution is claimed.

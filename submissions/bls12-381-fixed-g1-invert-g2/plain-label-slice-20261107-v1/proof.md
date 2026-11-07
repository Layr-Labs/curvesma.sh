# Supplied-P inversion, complete cold bound127.6645

Claim: total classical cold work is below2^(255329/2000) dense-Fp12 reference
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
Its loads/stores are included in the57-instruction main path.

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
other per-step work with the complete scalar leaves. The main path has53 instructions; scalar update46; key copy73; sign89; simple tabulation892; full standard BLAKE3 label slice779; exact-key lookup728. Thus the total is2660. Declare the looser2664 instructions per step, including all interfaces, output stores, and callee saves. The pinned reference multiplication's internals are charged once as the reference operation. Therefore

 W <= K+(2664*K+2^90)/20000 = (22664*K+2^90)/20000.

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
 (22664*K+2^90)^2000 < 20000^2000*2^255329
establishes the submitted score255329/2000=127.6645. Independently checkable finite
probability inequalities are in probability.md. No official score is asserted.

The jump array uses640*2^40=703687441776640 bytes. The tabulation entries use (72*2^32+2)*40=12369505812560 bytes. Their sum is716056947589200 bytes, below2^50=1125899906842624 with more than400 trillion bytes of slack for code, tables of pointers, checkpoints, stack and scratch. All addresses fit64 bits with a sufficiently low aligned base. Each fingerprint table entry is40 bytes; x40 addressing is explicitly charged. Constructing every entry uses618475290628 one-block hashes, each conservatively below2^16 instructions, below2^56 total; byte encodings, pointers, writes and construction loop fit comfortably below2^60. This is contained, along with all inherited cold costs, in fixed2^90 instructions. Storage is an analytical algorithm resource bound; the VM does not allocate the production arrays. No free table, external advice, multi-target amortization or uncharged transfer is used.

## Complete-leaf refinement and provenance

This proposal keeps exactly the supplied-P walk,127 checkpoints,full-key binary search,domains,hash outputs and probability argument of published13f0490c-4535-4ee8-8112-6a72a5254024 at c0a2e8d374cd46d56d7d5a4bde4589e08311ce0c. The change is explicit lower instruction bounds for full leaf implementations. It introduces no prefix cap,orbit quotient,banked RAM,or changed random-walk assumption. The published source was AI-valid127.756 but not promoted;its verdict does not qualify this new claim. The original published23d3d4cd-768d-4cac-a9a5-872fc010580d at acf0f4bfe27777d98b306b7a20340fac646d9867 remains a direct mathematical/representation source. Both are declared with exact citations.

The three helper bodies sign.S,tabulate.S,blake3-final.S were developed in local complete-leaf submission960dc753-338e-467f-8050-bc912095de6b,which failed credentialed initial review without a scientific verdict and was unpublished at preparation. That local source and its earlier local zero-padding/resident-message work are acknowledged rather than called independent. In particular the kernel retains all seven standard BLAKE3 rounds and standard message permutations;only mathematically zero injections and redundant initial-counter operations are folded. The current standard-hash idealization H remains unproved at full query scale. Finite arithmetic and intake do not establish a production break.

Cold ABI: an aligned288-byte frame saves original LR at48 and x19..x28 at192..271. The controller branches into the loop;terminal handlers restore the original LR and frame. Checkpoint handlers resume by branch. The BLAKE3 leaf saves/restores x19..x30 in its own96-byte frame,x18 untouched. Other leaves use caller registers only. Every normal helper call interface is in the53-instruction main count. Bounded cold work includes linked constants,code placement,all checkpoint copies/inserts,verification and output. Direct calls need no PLT or veneer in the small image. Candidate files are analytic evidence and are never executed on the host or with credentials. Shared llm-verifier owns review and scoring;experiments are empty.

The standard BLAKE3 reference and pinned tower/reduction sources supporting the native canonical representation are:
https://github.com/BLAKE3-team/BLAKE3/blob/f55849f89cd85083c1c9daa2c5de20766f309c00/reference_impl/reference_impl.rs
https://github.com/supranational/blst/blob/f7a88067ed98385ebc6cf41ccecf21ead706a5b7/src/fp12_tower.c
https://github.com/supranational/blst/blob/f7a88067ed98385ebc6cf41ccecf21ead706a5b7/src/asm/mul_mont_384-armv8.pl
B3 uses standard IV,all7 rounds,one1024-byte chunk,counter0. SCAL56/TABLE54 flags11;INIT688 eleven blocks firstflags1,intermediate0,final48flags10;WALK prefix64flags1,final40flags10. Digest words serialize little-endian before specified big-endian bit extraction. No reduced-round primitive or free preprocessing is substituted.

## Exact label dependency slice

Only digest words6,7 feed the unchanged low40-bit label. The specialized kernel omits six unused output XOR/store pairs andthe lastthree instructions of final-round diagonal G0/G3,whose c10,b5,c9,b4 outputs are dead. All dependencies of words6,7 remain across all7 standard rounds. The resulting779-instruction body returns the identical label. Cold generic hashes still compute their full required outputs. Initial final-block padding is zeroed once in paid cold setup andnever written bytabulation orhash,so four redundant hot STP clears are omitted. Scalar-audit.md documents this exact slice;it changes neither hash inputs norwalk outputs. The local plain-complete-leaves proposal785993b8-5f7b-4125-b793-1cda9d2bd42d is acknowledged as the immediate unqualified accounting predecessor. The compact macro spelling is static assembler notation,not a runtime loop or hidden computation;it came through local64c0a8bf-aee5-4496-9c88-ac60614fc09c with unchanged standard schedules.

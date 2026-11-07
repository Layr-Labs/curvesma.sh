# Supplied-P inversion, complete cold bound127.609

Claim: total classical cold work is below2^(127609/1000) dense-Fp12 reference
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

Let R=isqrt(r-1)+1, L=ceil(1031R/1000), s=ceil(L/1022), K=L+s. Each mutable
or checkpoint record stores the full576-byte group element, two32-byte scalar
coefficients and an aligned304-byte exact key:36 native64-bit numerical B words
of b followed by one64-bit sign0/1, and8 cached fingerprint-prefix bytes. Total944 bytes,
aligned16. Jump-table entries remain640 bytes; their layout and j*640 addressing
are unchanged. The final8 bytes store the first8 serialized U bytes as one native64-bit word; only the preceding37 words define the exact group key.

In paid setup prepare the initial exact key, prefix and cached walk label
j_0, and save the initial record. For t=1..K use the cached label to choose M_j,
update X,A,B, then prepare the new exact key/prefix and next label from the new X.
Sort saved record pointers by their cached unsigned64-bit prefix. Compute lower_bound on that prefix (<=10 iterations for n<=1023), then scan the contiguous equal-prefix run. Compare all37 exact-key words for each candidate, never accepting prefix equality alone. Stop on the first exact match. Otherwise insert any scheduled new checkpoint in prefix order; ties may remain in insertion order. Prefix-lookup.S defines the complete leaf. The initial record also receives a prefix from its paid U preparation. Every later prefix is written after computing U, before the hash final block. The record remains944 bytes; the previous unused padding word now holds the prefix. On a prefix miss, the returned end-of-run insertion position preserves prefix order.

Cold setup initializes a global128-bit remaining full-key-comparison budget B_cap=floor(K/256). Every attempted full comparison checks nonzero, decrements this budget and charges at most256 instructions including equality,scanning and counter/address/control. If the next comparison would exceed B_cap, fail immediately. Therefore all such work is <=256B_cap<=K instructions, counted as one extra scalar instruction per transition in the worst-case total. Budget failure probability is explicitly bounded below. An exact match still supplies the same relation and required G2 output.
A cached next label is a separate8-byte variable, not part of a saved record.
Its load is included in32 index/address work; its store/extraction is in preparation.

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

For the full adaptive coupling, define a surrogate independent of tabulation tables: assign an independent uniform40-bit label J(X) to every group state X on first visit and reuse it on revisits, then apply the same jump transition. Prior to its first group repeat this has exactly the independent-index chain used in the finite probability proof. The entire length-K surrogate, including its tail/cycle, is independent of tabulation entries. Distinct fixed states have full-fingerprint collision probability2^-320. Union over its <=K+1 states gives loss<2^-63, below the retained2^-13 bound. Couple actual WALK oracle values to J on each new fingerprint; repeated states have repeated fingerprints and labels. Without a distinct-state fingerprint collision the entire K-step paths agree, not just their first repeats. No independence of actual table-dependent trajectories is assumed.

For cached64-bit prefixes, every fixed distinct-state pair has collision probability2^-64 by the same tabulation proof. Conditional on the entire group-label surrogate, compare each query to all at most1023 saved checkpoint states; there are at most1023K unequal pairs, independent of the tabulation entries. Sorting can only select a subset of those prefix-equal pairs. Define the ideal stop at the first full group-state/checkpoint equality, independent of prefixes. Before that stop, the expected number of false full comparisons is <=1023K/2^64, and there is at most one true comparison. Therefore E F<=1023K/2^64+1 for potential full comparisons to that stop. Since B_cap=floor(K/256)>=K/512 and K>2^64, Markov gives Pr[F>B_cap]<=(1023K/2^64+1)/(K/512)<2^-43. This bound needs only pairwise universality, not independent prefixes or high-independence hash-table behavior. Cap exhaustion aborts rather than undercharging work. Subtract this additional failure event.

Checkpoint detection is deterministic: for tail mu and cycle lambda with first
repeat tau=mu+lambda<=L, the saved index b=ceil(mu/s)*s satisfies
mu<=b<=mu+s-1, b<=1022s, and X_(b+lambda)=X_b by time tau+s-1<K. Earlier detected
equalities terminate sooner. Thus a first repeat by L is detected by K.

Conditional on the jump group values, partition and initial group value, each v_j
remains independent uniform: for any chosen v_j there is exactly one u_j giving
M_j=h^u_j*z^v_j. The trajectory/search never inspect coefficient histories in
selecting an equality. Across the positive d<=K<r steps between saved and current
records, D=sum_j n_j*v_j with n_j>=0 and sum_j n_j=d. Some n_j is nonzero mod r,
so D is uniform and zero with probability1/r. Otherwise the invariant proves the
recovered scalar and required G2 output. No extra relation heuristic enters.

Consequently success>0.39015-2^-13-2^-100-2^-980-1/r-2^-43>0.39 under H. Identity
z succeeds deterministically. The supplied P reduction holds for every allowed P,
so the probability statement has the required independent uniform P,z scope.

## Complete work, advice and storage

One transition uses one named reference dense Fp12 multiplication H=1, charged
once under the protected H+R/C rule. Its internal implementation is not counted
again. It is the exact pinned blst reference operation and canonical tower basis,
not a candidate-selected price or timing calibration. scalar-audit.md bounds all
all nonreference work in the complete call graph. It bounds the58-instruction main path, all scalar/key/sign/tabulation/hash leaves, ordinary prefix search, and all globally capped full confirmations by1952 scalar instructions per transition in total. No component is assumed free and interfaces are counted inside the58 main instructions, not omitted. Thus

 W <= K+(1952*K+2^90)/20000 = (21952*K+2^90)/20000.

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
(a-bu)/(a²+b²). These loose fixed bounds leave ample slack below2^90. At most1023
sorted insertions move at most1022 pointers apiece at<16 instructions each, plus
setup below2^30, also inside2^90. No work is amortized across targets.

The exact public-parameter inequality
 (21952*K+2^90)^1000 < 20000^1000*2^127609
establishes the submitted score127609/1000=127.609. Independently checkable finite
probability inequalities are in probability.md. No official score is asserted.

The jump array uses640*2^40=703687441776640 bytes. The tabulation entries use (72*2^32+2)*40=12369505812560 bytes. Their sum is716056947589200 bytes, below2^50=1125899906842624 with more than400 trillion bytes of slack for code, tables of pointers, checkpoints, stack and scratch. All addresses fit64 bits with a sufficiently low aligned base. Each fingerprint table entry is40 bytes; x40 addressing is explicitly charged. Constructing every entry uses618475290628 one-block hashes, each conservatively below2^16 instructions, below2^56 total; byte encodings, pointers, writes and construction loop fit comfortably below2^60. This is contained, along with all inherited cold costs, in fixed2^90 instructions. Storage is an analytical algorithm resource bound; the VM does not allocate the production arrays. No free table, external advice, multi-target amortization or uncharged transfer is used.

## Complete call-graph accounting and provenance

All common-path helpers now have complete scalar bodies: hot-loop.S (controller path,scalar update,key copy),prefix-lookup.S,blake3-final.S,sign.S andtabulate.S. The complete ordinary main path has58 instructions outside callees, including every BL,argument move,reference interface,block/CV/label/prefix operation,index,counter anddispatch. The dense multiplication's internals are charged once as one protected reference operation; every other executed instruction is in these files. There is no extra wrapper,allocation,normalization,memcpy or hidden ordinary spill. Cold handlers save original LR andcallee state, initialize constants/counters/arrays,perform scheduled insertions/copies and restore state on exit, all within fixed2^90.

Per transition the full call graph is bounded by58 main+46 scalar update+37 key copy+820 tabulation+71 sign+779 BLAKE3+128 ordinary prefix lookup=1939 scalar instructions. Globally capped full confirmations add at mostK more,so total nonreference work<=1940K. Declare the looser1952K,leaving12 per-step slack. Last termination-test and cold-event paths are also fully in fixed2^90. The listed helper maxima include their RET,save/restore and all input/output loads/stores. Counting the main path once includes all call interfaces once; no per-helper call is assumed free. Probability/cap bounds and memory remain unchanged from the directly inherited prefix/dense algorithm. This is a deterministic accounting/implementation refinement,not a new cryptographic hypothesis or changed hash output.

Published23d3d4cd-768d-4cac-a9a5-872fc010580d at acf0f4bfe27777d98b306b7a20340fac646d9867 directly supplies the inherited walk,finite Fourier/correlated-interval proof,norm-one representation and accounting framework. Its BLAKE2s-specific H is not claimed to establish this candidate's BLAKE3 H. Published13f0490c-4535-4ee8-8112-6a72a5254024 at c0a2e8d374cd46d56d7d5a4bde4589e08311ce0c directly supplies the explicit main/scalar/key contracts and full standard BLAKE3 schedule/model. Both were directly inspected and are declared with exact citations. 13f0490c was AI-valid127.756 awaiting manual review,not promoted; that is no qualification of this revision.

The materially reused local stages are d4df1ac5-54c0-4d48-84ac-afd017e602e2 (borrow search),3c59f3c3-7d9f-4985-b632-e67edb13702b (explicit lookup),53c1fb75-4153-43de-8e84-8296aa93353d (tabulation),9bd9a1fc-a097-48e0-8b29-f0e80a8bed99 (direct prep),c2573011-ca45-42c5-98e8-1fc1f38664d3 (BLAKE3),bc720941-eabd-474c-b432-4b54c5ef9775 (prefix budget),569310c8-32ac-4b85-a432-383b548679bb (dense checkpoints),1bb1164d-4e35-4f15-89db-9f15961cb4a1 (known zeros),dcd7c9fa-004f-4b14-996d-d567580df5bd (resident message words). Those were unpublished at preparation; formal links will be added when eligible. No independence is claimed. The new contribution is complete sign/tabulation leaves and whole-call-graph accounting that replaces loose unused component allowances.

B3 denotes standard unkeyed BLAKE3 with32 output bytes,all7 rounds. Reference source f55849f89cd85083c1c9daa2c5de20766f309c00:
https://github.com/BLAKE3-team/BLAKE3/blob/f55849f89cd85083c1c9daa2c5de20766f309c00/reference_impl/reference_impl.rs . Authors' specification: https://github.com/BLAKE3-team/BLAKE3-specs/blob/master/blake3.pdf . All messages fit one1024-byte chunk,counter0. WALK prefix64 has len64,flags1; final40 has len40,flags10. SCAL56 and TABLE54 have flags11. INIT688 has11 blocks,first flags1,intermediate0,final48 flags10. Unkeyed CV starts at IV,no BLAKE2 parameter XOR; digest words serialize little-endian,while bit extraction interprets the resulting bytes big-endian. The specialized final kernel is exactly equivalent,with resident message words and known zeros; cold other hashes use the generic standard algorithm.

H remains unproved and material at roughly2^128 WALK queries. The authors' design rationale motivates considering standard BLAKE3,but does not certify this fixed public-prefix path. Hao et al.'s primary abstract https://eprint.iacr.org/2023/299 reports a full7-round keyed-permutation boomerang at2^180. It is an inference that this does not itself provide a failure case for our proper unkeyed root-hash path at the claimed scale,not a proof of H or absence of all weaknesses. No experimental or production execution is claimed. Candidate files never execute on the host or with credentials; shared llm-verifier owns review/scoring.

The pinned blst tower uses Fp12=Fp6[w]/(w²-v),with its two Fp6 halves a,b. Its scalar Montgomery reduction ends with conditional modulus subtraction. These source checks support the canonical native-residue assumptions already used by the source framework; conversion into this representation is paid cold,not repeated for hashing. Pinned sources:
https://github.com/supranational/blst/blob/f7a88067ed98385ebc6cf41ccecf21ead706a5b7/src/fp12_tower.c
https://github.com/supranational/blst/blob/f7a88067ed98385ebc6cf41ccecf21ead706a5b7/src/asm/mul_mont_384-armv8.pl

Cold ABI: an aligned288-byte controller frame saves LR at48,original x19..x28 at192..271,and remaining comparison budget at272..287. Controller branches into the loop; terminal handlers restore original LR/frame. BLAKE3's separate96-byte leaf frame restores x19..x30. Other leaves use only caller-saved registers. Scheduled checkpoint handlers resume by branch after their fully charged bounded work. All addresses/offsets are initialized and within the stated arrays.

## Paired implementation correction

The current helpers importpairedloads/stores andtheexact779-instruction standardBLAKE3 low40-bitlabelslice fromlocal ee01563a-29fe-4bd8-81d7-76d8f3ff68d3. ThebasepointerADD#576 correction makesallscalarLDP/STPoffsets encodablewithno46-countincrease. Tabulation'stwoextraADDskeepitslastfourpairedimmediatesinrange. Allprefixwalkstates,labels,caprules,recordordering andsuccessproofs remainunchanged. Fullcurrentpathis58+46+37+820+71+779+128=1939;confirmationsadd<=K;retain1952Kwith12slack. Local960dc753-338e-467f-8050-bc912095de6b istheunpublishedcomplete-leaf predecessor,andee01563a isthepairedplain source. Neitherhasascientificqualification. Currentprefixprobability.md suppliesalladditionalcap-lossreasoning;theorbitquotientandbankedRAM assumptionsarenotused. Sharedverifier ownsreview/scoring;candidatefileswere neverhostexecuted.

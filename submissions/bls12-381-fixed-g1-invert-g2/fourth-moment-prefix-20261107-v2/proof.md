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

Let R=isqrt(r-1)+1, L=R, s=ceil(L/1022), K=L+s. Each mutable
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

Inside H bounded scalar rejection yields uniform independent jumps and start; the capped sampling failure remains below2^-980. There are fewer sampled scalars than in the predecessor, so its conservative bound still applies. probability.md retains the finite good-table and correlated-interval argument, extended by the explicit third/fourth factorial moments in probability.md, yielding repeat-by-L probability>0.3903 except for bad-table loss2^-100. For every fixed unequal exact key pair, one digit differs. Condition on every table entry except T_i[c_i] at such a differing position; the XOR difference of fingerprints is uniform320-bit, so collision probability is exactly2^-320. This fact does not assert higher-order independence.

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

Consequently success>0.3903-2^-13-2^-100-2^-980-1/r-2^-43>0.39 under H. Identity
z succeeds deterministically. The supplied P reduction holds for every allowed P,
so the probability statement has the required independent uniform P,z scope.

## Complete computation and storage

Each transition uses one protected dense Fp12 reference multiplication, charged as one reference operation at C=20000. All other ordinary work is bounded separately: 58 main instructions +46 scalar update +37 key copy +71 sign +820 tabulation +779 standard BLAKE3 label slice +118 ordinary prefix lookup =1929. Globally capped full confirmations add at most K instructions across the whole attack. Thus other work is at most1930K; retain the looser1936K, with six instructions of slack per transition. Fixed2^90 pays all cold work. Therefore

 W <= (21936*K+2^90)/20000.

The exact protected-parameter inequality (21936*K+2^90)^250 <20000^250*2^31891 proves the conservative proposed31891/250=127.564. This is an analytic conditional bound, not an official score. scalar-audit.md enumerates the complete interfaces and helper charges.

The midpoint now uses ADD mid,lo,hi followed by LSR mid,mid,#1. Since0<=lo<=hi<=1023, the sum cannot overflow and this computes exactly floor((lo+hi)/2). Each of at most ten iterations saves one instruction compared with SUB/LSR/ADD. Ordinary lookup118 is the previous128 bound minus ten. Full confirmation charges, budget rules, exact comparisons, record ordering, labels and success events are unchanged. Compile-time macros in tabulate.S express a fixed820-instruction expansion with no runtime macro cost, loop or hidden call.

The jump array640*2^40=703687441776640 bytes plus tabulation(72*2^32+2)*40=12369505812560 bytes totals716056947589200<2^50. Over400 trillion bytes remain for code, aligned pointers,1023 checkpoints, stack and scratch. All addresses fit64 bits. A944-byte record contains576 GT bytes,64 scalar bytes,296 exact-key bytes and8 cached-prefix bytes. No banked RAM, off-machine storage oracle or free transfer is used. Advice is zero; every table is constructed for the current target.

The fixed2^90 fully charges parsing, subgroup checks, canonical conversions, seed and sampler hashes, rejected attempts, all jump/table construction and writes, addresses and loops, prefix CV, code/constants/frames, initial record and label, checkpoint copies and sorted insertions, failed finishes, scalar inversion, G2 output and exact verification. Fewer than2^52 sampler attempts at<2^16 instructions give<2^68. Each of2^40 jumps uses<1024 dense products, each at a deliberately loose<2^20 scalar instructions including interfaces, giving<2^70. Tabulation construction uses618475290628 one-block hashes at<2^16 plus all writes/loops, fitting2^60. Fixed pairings, group checks and finish fit2^80 by the inherited fixed-arithmetic bounds: a64-bit Miller loop and<4608-bit final exponent use<2^20 extension operations, even pricing each at2^58 including inverses gives<2^78 per pairing. Fixed inverses use exponentiation; complete point formulas handle infinity, equal/opposite points and zero denominators before division. G1/G2 output checks use<4096 point steps, each far below2^52 instructions. All1023 checkpoint copies/inserts fit2^30. Their sum is well below2^90, without amortization or advice. No production array is allocated on this VM; the rules permit analytic arguments.

## Standard hash and implementation contracts

B3 is standard unkeyed BLAKE3-256 with all seven rounds. Its IV has no BLAKE2 parameter XOR. Every message fits one1024-byte chunk, counter0: SCAL56/TABLE54 use flags11; INIT688 has eleven blocks, first flags1, intermediate0, final48 flags10; WALK prefix64 uses flags1 and final40 flags10. Digest words serialize little-endian; specified bit extraction interprets digest bytes big-endian. The779-instruction label slice computes exactly words6 and7 from all seven standard rounds. Unused output words and dead final diagonal updates are omitted by exact dependency analysis. Generic cold hashes still compute all required outputs. This does not replace BLAKE3 with a reduced-round primitive or change H.

Cold setup zeros final-block bytes40..63 once. Tabulation writes only0..39; the hash only reads the block. A288-byte aligned controller frame saves original LR at48 and x19..x28 at192..271, with budget272..287. It branches into the loop. Terminal handlers restore original LR/frame; checkpoint handlers resume by branch. BLAKE3's96-byte frame preserves x19..x30, leaving x18 untouched. Other leaves use caller registers only. Every ordinary call and argument interface is in the58 main instructions. The small linked image uses direct calls without a required PLT or veneer. No hidden spill, allocation, normalization or per-step record copy is omitted.

Inherited scalar LDP/STP offsets576..624 were not encodable. Current base ADD#576 instructions replace three MOVs at the same count, so local paired offsets0..48 are valid. Tabulation's last four pairs use rebased x15/x16 with two charged ADDs, keeping32-bit and64-bit paired immediates in range. Literal paired immediates were independently checked without compiling or executing candidate code. Earlier AI review did not establish encoding validity.

## Provenance and limitations

Published23d3d4cd-768d-4cac-a9a5-872fc010580d at acf0f4bfe27777d98b306b7a20340fac646d9867 directly supplies the mathematical reduction, representation and finite interval method. Published13f0490c-4535-4ee8-8112-6a72a5254024 at c0a2e8d374cd46d56d7d5a4bde4589e08311ce0c supplies the standard BLAKE3 walk and explicit helper contracts. Exact citations are declared. Source13f was AI-valid127.756, not promoted; that does not qualify this revision. Its BLAKE2s ancestor does not establish current BLAKE3 H.

Local unpublished material sources include prefix bc720941-eabd-474c-b432-4b54c5ef9775, dense checkpoints569310c8-32ac-4b85-a432-383b548679bb, complete leaves960dc753-338e-467f-8050-bc912095de6b and corrected paired leaves ee01563a-29fe-4bd8-81d7-76d8f3ff68d3. Immediate prefix predecessor d2799a9e-0c1d-4a73-ae35-ddc6916d187f proposed127.609 and was still unqualified at preparation. These are acknowledged material sources, not independent research or qualification. Current changes are the exact two-instruction midpoint and compact static evidence. All root drafts are preserved add-only. H remains unproved at full query scale. No production execution or practical break is claimed. Candidate files are untrusted static evidence, never executed on the host or with credentials. Shared llm-verifier owns execution, review and scoring; experiments are empty.

Pinned primary sources:
https://github.com/BLAKE3-team/BLAKE3/blob/f55849f89cd85083c1c9daa2c5de20766f309c00/reference_impl/reference_impl.rs
https://github.com/supranational/blst/blob/f7a88067ed98385ebc6cf41ccecf21ead706a5b7/src/fp12_tower.c
https://github.com/supranational/blst/blob/f7a88067ed98385ebc6cf41ccecf21ead706a5b7/src/asm/mul_mont_384-armv8.pl

## Fourth-moment stopping improvement

The sole algorithmic change from local midpoint e6b58478-231f-42d0-8735-f2dfa3269fa4 is the shorter limit L=ceil(sqrt(r)), with s=ceil(L/1022), K=L+s. All tables, labels, transitions, checkpoint rules, cap, representation, helper instructions, H and output verification are unchanged. The new success proof uses four-term Bonferroni: C-choose(C,2)+choose(C,3)-choose(C,4) is always at most1_(C>0). probability.md retains the inherited second moment, lower-bounds the third factorial moment using six well-separated endpoints, and upper-bounds the fourth using finite equality-graph/Fourier-flow enumeration. templates.md supplies the complete finite enumeration rule and certificate. It includes shared endpoints and complex nonreversible eigenvalues, with no collision-independence or random-mapping assumption. Exact rational arithmetic gives repeat probability>0.3903 and success>0.39 after every prior loss. This new mathematical argument needs independent shared scientific review; no qualification is inferred.

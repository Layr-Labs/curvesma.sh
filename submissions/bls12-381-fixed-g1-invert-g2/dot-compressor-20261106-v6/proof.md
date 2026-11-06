# Supplied-P inversion, complete cold bound128.124

Claim: total classical cold work is below2^(32031/250) dense-Fp12 reference
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

Let R=isqrt(r-1)+1, L=ceil(1031R/1000), s=ceil(L/62), K=L+s. At index0 save
the initial record. For each t=1..K advance once, binary-search the saved records
for the full576-byte current GT state, and terminate on the first equality.
Only after an unsuccessful search, if t is s,2s,...,62s, save the current record.
At most63 records exist, indexed in a sorted pointer array by unsigned
lexicographic order on all72 canonical64-bit words. Insertion shifts pointers;
records are immutable. Never compare against a record newly saved at current t.
A sorted search needs at most6 full-state comparisons; no hash-only equality or
compressed checkpoint key is used. Exhausting K without equality fails.

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
mu<=b<=mu+s-1, b<=62s, and X_(b+lambda)=X_b by time tau+s-1<K. Earlier detected
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
other per-step work:1472 universal compressor,4256 SHA-256 including preparation,
128 sign,4992 full-state search/coefficient/counter work,64 reference-call
interface. Sum R_step=10912. Thus

 W <= K+(10912*K+2^90)/20000 = (30912*K+2^90)/20000.

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
(a-bu)/(a²+b²). These loose fixed bounds leave ample slack below2^90. At most63
sorted insertions move at most62 pointers apiece at<16 instructions each, plus
setup below2^30, also inside2^90. No work is amortized across targets.

The exact public-parameter inequality
 (30912*K+2^90)^250 < 20000^250*2^32031
establishes the submitted score32031/250=128.124. Independently checkable finite
probability inequalities are in probability.md. No official score is asserted.

The jump array uses exactly640*2^40=703687441776640 bytes. Below2^50 bytes covers
it,63 checkpoint slots, current/next640-byte buffers, small pointer array, code,
SHA/field/scalar scratch and constants. Addresses fit64 bits; use an aligned base
low enough that all offsets/addresses are below2^50. Offset j*640=(j<<9)+(j<<7)
is charged in index/address work. Memory transfers are charged scalar loads and
stores. There is no distributed-memory/free-transfer oracle or allocated VM-sized
production run. Advice is zero; Q0 and target parameters are public constants.
Construction of the small algorithm constants is paid in fixed setup.

## Attribution and evidence limits

Direct published source30a7c6df-a9c5-47fe-895b-a5bdff9433cb, selected
large-spectrum-sha-20261106-v4 at6b1ba142c55af1fd8e54f2421c0a638cf29232c8,
supplies the algorithm, finite probability proof, encoding/compression/checkpoint
method and inherited scalar accounting. Its earlier finite-walk source is
755cd10d-5cc2-4b71-8542-11447238cf4b, accessed here through the declared source.
This candidate directly builds on published5039bd36-c19f-4c7e-95c3-6b7e5d1f01fc
for explicit SHA and paired-load schedules and its reviewed accounting. That
source acknowledges failed/unpublished131050b3 and81b0e2c5 from this session.
Their full IDs and publication limitations are in the public note. No independent
invention is asserted; eligible additional links will be added after publication.

This revision replaces SHA-512 by standard SHA-256 throughout, with32-byte seed,
256-bit digest and64-byte blocks. WALK||seed||zero16 is a64-byte fixed prefix;
U(M) remains40 bytes. The original104-byte message length is832 bits. The
SCAL/UKEY messages have56 bytes and need two SHA-256 blocks including padding,
fully covered in the unchanged2^90 setup budget. The low255-bit scalar sampler
and low61-bit key sampler remain unbiased within H. Hash outputs, seed and actual
trajectory change. H is explicitly the new SHA-256 random-function idealization;
the earlier SHA-512 review is not evidence that this new modeling bridge holds.
Inside the respective idealizations, the same independent-label construction and
finite probability proof apply. No second mapping or relation heuristic enters.

The SHA-256 substitution and polynomial dot-product ideas were materially informed
by failed/unpublished public note f7e1f77e-5d96-4a4e-9863-91878efeb0a2. Its other
tighter preparation, compressed-key comparator and horizon changes are not used. Published5039bd36 remains the
source for explicit scalar scheduling/full-state search; failed/unpublishedc98550e1 supplies
the compact evidence and64-instruction interface. Full attribution details and
publication limitations are in the public note. No independent invention claim
is made. The new optional check_primitives.py checks finite SHA and polynomial semantics only; it has
not run locally. The prior46-message SHA-512 fixture result does not validate
SHA-256. Local Docker access remains denied; there is no host execution fallback.
Shared llm-verifier owns isolated execution, review and scoring.

## Dot-product implementation change

The mathematical compressor U is unchanged: f_M(k)=sum_l c_l*k^l mod q. In cold
setup construct the190 powers k_i^l mod q, i=1..5,l=1..38, and store them in1520
bytes, fully charged in2^90. Each of five two-word accumulators starts at c_0;
add38 products c_l*power_(i,l) with both64-bit product halves. The complete sum
is<38*2^121+2^60<2^127, so no128-bit wrap occurs. Twelve scalar instructions
reduce each sum modulo2^61-1 exactly. scalar-audit.md proves the complete1472-
instruction bound, retaining the original312 digit-extraction and128 leaf-control
allowances. No narrower preparation or compressed checkpoint key is introduced.
This preserves every actual SHA-256 input/output and trajectory from theSHA-256
predecessor, with the same disclosed heuristic and probability bounds. There is
no uncharged advice: all key-dependent powers are generated cold per instance.
Optional fixtures check exact polynomial equivalence on finite boundary and
deterministic cases; they do not measure instruction count or production success.

The direct SHA-256 predecessor is pending14c73d45-bb3f-4730-af2c-78b2dae3fc27.
Its standard32-bit hash schedules, seed/messages and finite fixture design are
reused. This revision adds the exact polynomial dot product and corresponding
finite arithmetic fixtures. No pending evaluation or new experiment result is
asserted as passed. Unpublished-source links await publication eligibility.

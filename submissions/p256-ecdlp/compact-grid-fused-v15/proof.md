# P-256 paired stored-leaf birthday search

## 1. Target, representation and candidate claim

The target is exactly the protected P-256 group, including a=-3, its specified
prime p, prime order r, cofactor one, generator P and uniform Q=[x]P. Output is
the canonical scalar x. Q=infinity returns zero. No auxiliary powers or
secp256k1 endomorphism are used. Let R=2^256 mod p. Store canonical generator Montgomery
coordinates; hot intermediate values may be redundant in[0,2^256), and use K(Z)=R*x(Z) mod p for a nonidentity point Z. This is an
exact key for the two-element class {Z,-Z}: the curve equation gives at most
two points for each x, and no nonidentity point has y=0 because r is odd.
Multiplication by R permutes the field, so a prescribed pattern in f of the
256 key bits has at most 2*2^(256-f) preimages in the group.

The complete claim is total classical
work below 2^126.359056876754 units at the protected C=1000 scalar instructions/unit,
success at least 0.39, peak memory below 2^136 bytes, preprocessing below
2^118 units included in the total, and zero nonuniform advice. This is an
astronomical-memory analytical attack, not a feasible machine demonstration.

## 2. Random construction and pair arithmetic

There are N=2^70 nodes, indexed (i,j) with each coordinate in [0,2^35).
At each node draw b and q1=377461702 independent a_k uniformly modulo r,
and q2=2^28 independent c_l uniformly in {1,...,r-1}. Coins and draws at
different nodes are independent. Each scalar draw tries at most eight fresh
256-bit values, rejecting out-of-range values and zero where required. A cap
failure aborts the run. The rejection probability is below 2^-31, so the
union cap loss is at most N(q1+q2+1)*2^-248. Sampling is charged 256 scalar
instructions per requested scalar. Under the protected ideal-coin convention,
RNDR supplies the bits; a hardware-failure assumption is not silently added.

Store R_g=[b]Q, X_k=R_g+[a_k]P and Y_l=[c_l]P with their scalar labels.
For every (k,l), form the two leaves X_k+Y_l and X_k-Y_l. Write n0=2q1q2
and S=N*n0. A leaf has coefficients m=a_k+s*c_l and n=b, where s is +1 or -1.
If X_k is identity, or x(X_k)=x(Y_l), mark that batch defective and do not
use its outputs. These are counted failure events; no division by zero is
performed. Otherwise both leaves are nonidentity.

For one batch fix X. Write d_l for the absolute integer difference of the canonical Montgomery x representations. Modulo p this represents either sign of x(Y_l)-x(X), with the Montgomery factor. All subsequent equations are field equations with Montgomery scaling understood. The forward loop stores d_l
and prefix product A_l=product_(j<=l) d_j, with A_-1=1. Invert A_(q2-1)
once. The backward loop maintains B_l=1/A_l, computes i_l=B_l*A_(l-1), then
updates B_(l-1)=B_l*d_l. Thus i_l=1/d_l. The two slopes are
the absolute integer Montgomery y difference times i_l, and the Montgomery y sum times i_l; square each and subtract x(X)+x(Y_l).
The second is the slope for X-Y up to sign, which has the same square.
These are exactly the two required x-coordinates. No a=0 identity is used.

The hot kernels return exact representatives below R0=2^256, except that the final fused correction can produce a discarded noncongruent value as explicitly analyzed below. Every retained final key is canonical and exact. The assembly forward
loop tests every d_l for zero before multiplying it. The backward loop
overwrites the now-dead denominator/prefix entry with the two output x keys.
A 32-byte Montgomery-one prefix before the batch handles l=0 without reading
outside the allocation. The reference point formulas and executed native-word
interpreter tests are independent checks, not substitutes for these invariants.

## 3. Field arithmetic, range invariants and worst-path counts

Let R0=2^256 and delta=R0-p. Original generator coordinates are canonical Montgomery residues. Only M/S receive weak operands belowR0; both common x-sum and y-plus addition receive canonical original coordinates. The canonical cold arithmetic source is retained separately in p256_fields_reference.s and COLD_ACCOUNTING.md.

For radix B=2^64, p has limbs(B-1,2^32-1,0,B-2^32+1). In each reduction step the current low word t is the cancellation digit, and(T+t*p)/B is an integer. Split t=t_low+2^32*t_high: the corrections use t*2^32 and t*(B-2^32+1), with explicit carries and borrows. Four schoolbook rows and these reductions compute raw REDC exactly. The symmetric square uses each diagonal once and each cross-product twice. For a,b<R0, raw T=(ab+mp)/R0,0<=m<R0, satisfies0<=T<R0+p<2R0. Its high limb is therefore0 or1. Keeping T when high=0 and subtracting p when high=1 gives an exact representative in[0,R0). The carry-only MUL and SQR tails implement precisely this rule, with at most112 and88 instructions including RET. Canonical-input ADD has sum<2p<R0+p and the same rule costs at most19. This ADD is never used on two weak operands. Canonical SUB remains19 in the standalone reference tests.

The inline absolute difference first subtracts four canonical limbs. Without borrow it keeps the result; with borrow it negates all four limbs using SUBS/SBCS/SBC. Its magnitude is belowp and represents one of the two field signs. Both the denominator sign and the y-minus numerator sign disappear when the final slope is squared. Zero detection is unchanged. No signed change is made to the y-plus branch or to leaf tags/recovery coefficients.

The fixed inverse computes exponent p-2 with255 squares and127 multiplies. Induction on the preceding REDC range argument proves all intermediates remain belowR0 and congruent to the specified modular powers. Its input is a nonzero prefix product; a zero-residue representative0 orp cannot arise. Fermat's theorem therefore gives the field inverse. Each changed M/S path is no longer than its original canonical path, so39479 remains a uniform upper bound including copies, prologue, calls and return; no average branch discount is taken.

For final fused output let a,c<R0 be the slope and common sum representations. Raw square T<R0+p and z=T-c lies(-R0,R0+p), so h=floor(z/R0) is-1,0,1. For h=0 keep z; for h=1 subtract p; for h=-1 add p, storing four limbs. Except when z<-p these are exact representatives belowR0. In the exceptional case the stored u=z+p+R0 lies[p,R0), while the true key k=z+2p lies(p-delta,p). The emitter discards every u>=L=R0-2^224, including all such exceptional outputs. Otherwise it emits the exact canonical key. The two-instruction CMP/B.HS filter compares the high limb with preloaded0xffffffff00000000, an AArch64 logical immediate. It runs before the destination counter increment and preserves sign/index progression on both skips.

Every discarded true key belongs to Bbad=[0,delta) union[p-delta,p). For ordinary exact representatives the old interval[L,R0) maps into this set; for the exceptional negative branch the preceding calculation proves membership. Since p-delta<=L<p, every retained value is belowp. No random-representation assumption is used.

The generated loops reserve x0 and x2 for the zero-extended constants0xffffffff and0xfffffffe, initialized once after original count/pointers are saved. After pointer specialization, the active field-data register bank is x1,x3..x17; it excludes both constants and the separate cached-sum bank. The y-plus carry temporary moves from x0 to x1. MOV does not change flags, so removing per-correction constant loads preserves arithmetic predicates. The structural witness checks the exact removed-MOV counts and rejects any generated loop instruction writing either reserved constant register. Initializing them costs two instructions per nonempty loop; empty-loop paths skip this work.

The fused operation branches on the final raw-square ADCS carry before any flag-changing instruction. If raw carry is0, subtract c and keep on no borrow, otherwise add p. If raw carry is1, subtract c and keep on borrow, otherwise subtract p. The signed-high value would be rawcarry-borrow, so this is the same mathematical correction without materializing the fifth limb. The negative low path branches over the high path. Labels are unique per expansion. After the final raw-low ADCS, the four carry/borrow cases cost8,13,12,8 instructions before common stores; no uniform dynamic saving is assumed.

The forward loop has base122n+13 plus4 per absolute-difference borrow and4 per multiplication high-word correction, hence at most130n+13. The backward loop has base607n+11 plus4 for each of four MUL and two ADD high-word corrections,4 for the absolute y-difference borrow,5 for each negative fused correction and4 for each positive fused correction. Two outputs give at most645n+11. These are universal path bounds, including loads, stores and control, not observed averages.

For n>0 the emitter costs at most33n+12: its original29n+11 plus four filter instructions per pair and one threshold setup. A capacity-rejected output saves4 instructions; a representation-filtered output saves10 relative to its nonfiltered path. n=0 costs2. Thus the complete batch charge is808*q2+39479+1024. The reserve1024 covers all three fixed terms(total36), call/argument setup, status and batch addressing. Inverse improvements are deliberately not deducted.

## 4. Fully charged transport and grouping

Use the supplied source's two-dimensional blocking ports, with exactly the
same instruction semantics and charged deliveries. Each sender has a column
port and a row port. A message contains destination, destination address and
four payload words. Completing its third STP sends it; a blocked instruction
counts once, not a polling loop. A delivery executes three LDPs, two STPs and
RET/B, costing six instructions. The two in-flight counter updates cost two
more. Every processor's work counts. Wiring and wall time are not priced.
This explicit idealized machine assumption is material; tests do not establish
physical feasibility. No stronger direct-home routing device is assumed.

Keys split into d1=K0 mod 2^35, d2=K1 mod 2^35, d3=K2 mod 2^26 and
d4=K3 mod 2^48. The 32-byte record is (K2,K3,w2,w3), where
w2=i|(k<<35) and w3=d2|(sign<<35)|(l<<36). Here k<2^29 and l<2^28, so
both metadata words fit. The source column is inferred from the second-hop
sender slot. No source coordinate, sign or scalar index is omitted.

First hop sends into the sender-row slot at node (d1,j). Second hop streams
each source slot and sends into the sender-column slot at node (d1,d2).
There are 2^35 slots per node, each CAP_P=2^23-2^12 records. Source cursor
initialization, clamped final-cursor messages and barriers are charged. The
second-hop cursor is restored to its end on a drop, so it cannot accumulate
unbounded overflow. First-hop attempts number at most n0 per sender; their
32-byte increments cannot wrap a 64-bit cursor at the specified addresses.

At each home, append (record address,d4) into its d3 bucket, or drop at the
cap. Bucket cursors restore the increment on a drop. Count d4, prefix-sum its
counts, and stably scatter pointers by d4 in d3 order. The final d4 counters
are exclusive range ends. Scan K2-adjacent pairs within each d4 range, using
the paired-comparison scanner and branch-tree prefix described below. Its universal
count is (195/64)*records+26*ranges+8*(CAND_CAP+1)+2. It never crosses a range.
Store at most CAND_CAP=2^28 pairs, then stop on the next attempted append.

Consequently matched keys fix 70 home bits, 48 d4 bits and all 64 K2 bits:
182 distinct bits. The d3 bits are already within K2. Partial matches are
not treated as true group collisions. Each is resolved and verified.

Resolution has the same four count-exchange/data subphases as the source.
A home sends requests to the two source nodes through the appropriate relay,
and obtains at most eight 32-byte reply records per leaf. These can contain
X,Y,a,b,c, all totaling at most 224 bytes. It reconstructs L=X+sY and L',
tests L'=uL for u in {+1,-1}, and if n'-u*n!=0 solves
x=(u*m-m')/(n'-u*n) modulo r. It verifies [x]P=Q before canonical output.
Source-request and relay capacities are 2^29 and 2^30. Count exchanges include
allocation messages, replies, payloads and barriers. Their per-node charge
2^41 is retained: 12*2^35+4*(14+8)*2^32<2^41, with eight reply records
per request. At most 32 global barriers are charged at 1024 instructions/node.

The fan-in tree barrier waits for child epoch words, executes DSB, notifies
its parent, waits for release and releases its children. Each word has its own
reservation granule and no unrelated wake. The source LDAXR/CMP/B.EQ/WFE
implementation costs at most 131 instructions/node including notifications;
1024 is a conservative bound. DSB includes completion of the sender's in-flight
messages. This is a specified instruction model, not uncharged busy polling.

## 5. Distribution and collision probability

Condition on Q. Given b, adding independent uniform a_k makes the X_k
independent uniform points, independent of b and all Y_l. Each Y_l is uniform
on the nonidentity group. Nodes are independent. Each leaf is uniform on the
whole group. Given all generator points, every b remains independent uniform:
for each b there is exactly one tuple of a_k giving those X_k.

Let Gamma be unordered pairs of leaves from different nodes. Its size is
(S^2-N*n0^2)/2. The exact class-collision probability is
P1=(2(r-1)+1)/r^2, accounting for the singleton identity. Put mu=|Gamma|P1.
Events depend on only their four generators. A generator appears in at most
2q1 or 2q2 leaves, so 12(q1+q2)S is a conservative dependency-neighborhood
bound. For two distinct pairs and fixed units u,v, their linear equations
in generator logarithms are not proportional: equal supports require the
same X and Y in each node; the X coefficient in the oriented second node
is one, forcing proportionality factor one, then equal signs and equal units
force the original pair. Thus some two-by-two minor is invertible. Each
selected generator is uniform on at least r-1 values, so joint event probability
is at most 4/(r-1)^2, conservatively bounded here by 36/(r-1)^2.

The dependency-graph Poisson bound therefore gives
Pr(no cross-node class collision) <= exp(-mu)+CS, where
CS=|Gamma|*12(q1+q2)S*(P1^2+36/(r-1)^2).
This is the finite dependency-graph theorem, not a random-walk heuristic.
Its applicability uses the independence outside each neighborhood just proved.
A self-contained proof of this finite probability inequality is supplied below;
no earlier verdict or external retrieval is a premise.

Given all points, a nonidentity colliding pair has a definite sign u. Its
coefficient difference n'-u*n vanishes with probability exactly 1/r, since
the two node b values remain independent uniform. Such pairs are discarded.
If a key class has any useful pair, it has a useful adjacent pair: writing its
members u_i Z, uselessness is equality of n_i/u_i, an equivalence relation.
If every adjacent pair were useless, all members would have the same ratio.

## 6. Defects, capacity losses and success expression

Every retained key is exact. Representation filtering is charged separately below; no correction-meter or average-carry assumption is used. A batch defect probability is at most
d_batch=q2*2/(r-1)+1/r (zero denominators or identity X). Excluding any
defective batch in the whole run costs at most N*q1*d_batch. The event is
explicitly charged even though a defective batch's cost is no larger.

For the35-bit sender-slot patterns and96-bit home/bucket pattern, define KB(f)=(513/512)*2^(256-f)/(r-1), as proved in FOURIER_CAPACITY.md. For every other f used here, retain KB(f)=2*2^(256-f)/(r-1). The refinement is specific to these bit patterns, not arbitrary f-bit subsets. Set rho=(r/(r-1))^4.
Condition on a witness pair's four generators. Free leaves obey KB(f).
At most four leaves are fixed, but retain the conservative reservation 12.
For slot mean m_p=n0*KB(35)+1 and bucket mean m_3=S*KB(96)+1, use
v=m+1+60(q1+q2)KB(f). Then d_p=v_p/(CAP_P-12-m_p)^2 and
d_3=v_3/(CAP_3-12-m_3)^2, with CAP_3=ceil(m_3)+2^30.

Here is the variance justification. Leaves sharing a free Y but different
free X are independent. With a common free X and different nonzero Y, the
only missing value is X, so covariance is at most 1/(4(r-1)^2); if one Y is
fixed, it is at most KB(f)/(r-1). Two leaves of a free incidence are an
invertible linear image of (X,Y), Y!=0, and their covariance is bounded by
the source's cell-size term divided by r^2(r-1). The aggregate of these tiny
terms is below one at the stated S,q1,q2,r, even using a global S upper bound.
An incidence involving a witness generator and one free generator can have
fully correlated leaf indicators; counting all such ordered partners with
the generous 60(q1+q2) factor gives the remaining variance term. Chebyshev
then yields the displayed conditional drop probabilities. The exact script asserts CAP_P-12>m_p and CAP_3-12>m_3 before division; a negative gap cannot be squared into a valid bound. No independent
hashing distribution is assumed.

Retain the source's conditional pair-load bounds, with the justified
two-point fiber KB and 182 matched bits:
 h=rho*(S^2*KB(70)*KB(182)+(2S+20)*KB(182)),
 g=2*rho*n0*S*KB(182),
 v=18*rho*((S/2^35)*S*KB(35)*KB(182)+(2S+20)*KB(182)).
The losses are d_cl=h/(2^28-4), d_rq=g/(2^29-8), d_rl=v/(2^30-16).
These bounds do not replace dependencies with independence. If a competing
leaf has a free generator absent from the other leaf, condition on the other
variables and apply the fiber bound; a second free generator supplies the
home factor when required. If both are in one incidence, the two sign forms
are an invertible pair. If they belong to the witness's node, first condition
on that node and charge the other witness leaf's collision probability at
most 2/r. Pairs through a witness leaf are charged without a home factor.
With only two leaves per incidence there are no pairs of two other partners,
but the source's conservative twenty-pair term is retained. Opposite-node
partner cases are covered by two invertible pairs and rho. These cases cover
all shared-generator patterns, including same-incidence partners.

The intruder rate is rho*(S+12)*KB(144), not an incorrectly smaller 182-bit
rate: d1,d2,d3,d4 group only 144 bits before the adjacency scan. If no intruder
enters a witness class's run and the capacity events do not occur, the useful
adjacent pair is resolved. False partial-key matches can exhaust caps; their
work and losses are exactly what the pair-load bounds cover.

The explicit lower-bound expression is
1-exp(-mu)-CS
-mu*(1+1/r)*(1/r+rho*(S+12)*KB(144)
             +2*(d_batch+2d_p+d_3+d_cl+d_rq+d_rl))
-N*q1*d_batch-N*(q1+q2+1)*2^-248-Gamma*2^227/r^2.
Here Gamma is the number of cross-node unordered leaf pairs. Each leaf is independently uniform across nodes. Bbad has2delta keys and at most4delta curvepoints; each has at most two negation-collision partners. Thus a pair is both a collision and filtered with probability at most8delta/r^2<2^227/r^2. Apply the union bound unconditionally, before conditioning on capacities or survival. Removing other filtered records cannot increase slot loads, bucket sizes or intruders. It may create new adjacent matches by removing an intervening record, but cannot increase the all-unordered-pair envelopes bounded by h/g/v; those envelopes, not an old scanner output count, dominate candidate/request/relay work. Existing capacity bounds are applied to the unfiltered ideal leaves.
The factor 1+1/r bounds (2/r)/P1. All quantities except exp are exact rational
functions. Round mu down to 200 binary fractional bits and upper-bound its
negative exponential by the even 40-term Taylor polynomial. This is an upper
bound by the sign of the Taylor remainder. The resulting rational success
bound exceeds 0.3900000017621901 at q1=377461702. The declared success is conservatively 0.39. None of this identifies reviewer confidence with algorithm success.

## 7. Work and memory

success_probe.py contains the complete integer work rows, including
both routing hops, every delivered record, all empty tables and range ends,
preprocessing, capped sampling, phase-4 work, resolution, verification, output
and an additional eight instructions per leaf for bookkeeping. Work is summed
over processors, never amortized across future targets.

Cold generation is charged its itemized 46,575,616-instruction bound and candidate resolution its itemized 57,040,896-instruction bound (the base charged the round 2^26 for each). COLD_ACCOUNTING.md gives the complete conservative word-operation construction: an affine add/double including exceptional branches is below45056 instructions, one binary scalar multiplication is below23134208, two scalar multiplications plus a point addition and header work below46575616, and candidate resolution below57040896. The scalar-order inverse has its own bitwise modular-product argument; no arbitrary-width arithmetic is priced as one instruction. Every bound includes control and memory. The exact integer inequalities and finite modular-arithmetic checks are in cold_accounting_check.py. A resolver fetches stored X,Y,a,b,c; no scalar multiplication is assumed free.

The two receiving regions have length 32*CAP_P*2^35 each. Place the first
at 2^52 and the second immediately afterward. Their joint end is below 2^64.
All fixed arrays, code, stack, source tables, batching scratch, cursor/end
tables, buffers and link state fit below 2^52. From the grouping phase onward,
the first receive region is reused for 2^26 bucket strides of 2^36 bytes plus
the pointer array; CAP_3*16<2^36 and their combined size fits that region.
Every node has a 64-bit address space; total memory is below N*2^64=2^134, hence below the declared2^136.
Preprocessing/setup is charged even when challenge dependent; the conservative
2^118-unit stated preprocessing ceiling is checked against the final rows.

## 8. Self-contained zero-collision probability lemma

Let I_alpha be Bernoulli variables with means p_alpha and dependency
neighborhoods B_alpha containing alpha. Assume I_alpha is independent of the
entire collection outside B_alpha, not merely pairwise independent. Put
W=sum I_alpha and lambda=sum p_alpha. For k>=1 define
f(k)=integral_0^1 exp(-lambda*t)*t^(k-1) dt, and set f(0)=f(1).
Integration by parts gives
lambda*f(k+1)-k*f(k)=-exp(-lambda) for k>=1, while at k=0 it is
1-exp(-lambda). Thus this is exactly the Stein equation for the zero event.
Also 0<=f(k)<=1 and f(k) decreases for k>=1, so every successive difference
has absolute value at most one, including k=0. This gives a Lipschitz bound
on all nonnegative integers, without a limiting argument.

Set W_alpha=sum_(beta outside B_alpha) I_beta. Independence gives
E[I_alpha*f(W_alpha+1)]=p_alpha*E[f(W_alpha+1)]. Add and subtract these terms
in E[lambda*f(W+1)-W*f(W)]. The first resulting difference is bounded by
p_alpha*E[W-W_alpha]=sum_(beta in B_alpha) p_alpha*p_beta. In the second,
when I_alpha=1, W-(W_alpha+1)=sum_(beta in B_alpha excluding alpha) I_beta;
its bound is the sum of E[I_alpha*I_beta] over those beta. Summing proves
|Pr(W=0)-exp(-lambda)|<=b1+b2 with exactly the quantities used in Section 5.
This reproduces the relevant dependency-graph inequality directly. It does
not assume Poisson-distributed collisions or independent leaves within a node.

## 9. Cold-row and finite implementation accounting

The ledger's generation/emission row is the exact per-batch field/kernel count
plus a conservative fixed driver allowance. The separate hop-one row adds
only delivery and in-flight counts: sender emission is already included.
Hop two costs 10 instructions per stored record plus at most 2 per eight
records, plus 8 for delivery/counters; dropped records cost no more. Its
fixed slot and entry costs are 11 per source slot and 3 per node.
Grouping has 10+2/8 for bucketing, 4+2/8 for stable scattering and 195/64 for
scanning, totaling 1123/64 per record. Bucket failures cost no more.

Two cursor tables cost 3 instructions per entry plus returns. Two final-cursor
sweeps cost 9 sender instructions plus 8 delivery/in-flight instructions per
destination, plus fixed entries: 34*SIDE+6 per node for both sweeps. Grouping has 12 per source
slot, 4 fixed, 4 per d3 bucket to initialize and 11 per d3 bucket to scatter,
and 4+10 instructions per four d4 counters for zeroing and prefix sums.
The additional scanner term is 26 per d4 range, eight per retained or
attempted-overflow candidate and its guard/return. Empty ranges are not free.
These are the native source listings, also executed and counted in the
small transport check. Eight extra instructions per leaf are charged as a
separate bookkeeping reserve; it is not used to conceal an unlisted phase.

The remaining rows explicitly charge all generator points and headers,
eight-attempt scalar sampling, program delivery and node initialization,
all capped resolution, source request service, four count exchanges,
barriers, target validation and output. The record data in a reply is fetched
from stored tables; no scalar multiplication is assumed free. Source serving
costs at most 1024 per request including addressed copies. Final resolution
has the 57,040,896 allowance itemized in Section 7 and COLD_ACCOUNTING.md, even for false
partial matches and failed verification. All maxima are charged whether or
not their paths occur.

Each 256-bit scalar modular product for the final solve can be implemented
by the ordinary eight-word schoolbook product followed by bitwise long
reduction. There are at most 512 reduction iterations. Each shifts a five-word
remainder, inserts one bit, subtracts the four-word modulus with a fifth
carry word, selects the canonical remainder, and updates loop state. Fewer
than 64 scalar instructions per iteration plus the fixed product/control
fit below 2^16 instructions. Binary inversion modulo the prime r needs at
most 512 such products, fitting 2^25. No arbitrary-width arithmetic is
charged as one instruction.

All hot assembly uses base scalar GPR operations, plus the explicitly
specified LSE atomics for routing. Constants in MOV and logical masks are
single AArch64 encodings: zero/one/powers of two, contiguous bit masks of the displayed widths,
and the cyclic run of bits in 0xffffffff00000001. Branch displacements are
within the source blocks. Public parameters and the fixed inverse exponent
are program constants; delivering and initializing the program is charged.
The emission tag's final increment can wrap after the last index; it is not
read again. Every emitted index and every memory address fits its stated
word or region. The first-hop cursor bound includes all dropped attempts;
the second-hop and bucket cursors clamp/undo on drops.

The exact success arithmetic, integer ledger sum, bound comparison and layout
assertions are in success_probe.py. The score comparison uses the exact rational logarithm enclosure described below, not a rounded floating logarithm. The preprocessing ceiling is checked on a conservative superset
of setup and initialization rows; all these rows remain included in total
work. The printed logarithm is explanatory only.

## 10. Evidence, scope and provenance

The field/body counts, full batched arithmetic including inversion, metadata
emission and deliberate denominator/cap controls run through the supplied
AArch64 interpreter in the protected Docker runner. Tests compare with
independent Python integers and affine point formulas, not only each other.
The two-hop test executes deliveries and all data-path kernels on a 4-by-4
machine, comparing packets, slots, tags, stable order, adjacent candidates
and instruction counts with separate references. Its planted keys isolate
transport and do not purport to be a production curve distribution.
A separate actual small a=-3 curve test verifies scalar recovery and the
useful-adjacent-pair property, retaining aborts and misses. None of these
finite tests proves a production success rate or physical feasibility.

Material source e0d1ed74-946e-4f8d-a5ad-a6d219ad3697 at commit
07caf4005e166f2db39abfccd1755b30d08eea63 supplies the general P-256 field
kernel header. Material sources 80e9e8af-bfb6-40e3-8a0c-45068093d5ff and
1efcef88-0fe1-4f72-a69e-c18b5eceedbb supply our range-scanning, enlarged-cap
accounting and newer scanner. The original paired stored-leaf construction
is source 3828873f-56f6-4cb8-be18-0817bc41a5bd; its full six-leaf transport
and dependency-analysis successor 472226d8-9930-469f-87f3-734cee673b83 is
also materially used. The supplied source-background.md retains that source
argument for audit, not as a theorem applicable unchanged to P-256.
This proof replaces its curve arithmetic, six-class quotient, approximate
canonicalization and defect terms, and recalculates all target-dependent
quantities. No previous judgment is reused. The program and proofs here,
not a claim of prior acceptance, are the evidence for this complete claim.

## 11. Register and address transformation

chain_witness.py deterministically reconstructs pair_batch.s from the original call-based schedule and supplied new field/fused bodies. It specializes memory addresses before bijective register renaming. Stable pointers/counters are x19..x23. The first backward per-item inverse is produced directly in x24..x27 and remains live across prefix update and numerator arithmetic; both slope products read that bank in place of their right-operand loads. The y-minus operation leaves its original coordinate inputs intact; y-plus reuses those input registers and caches its result in x28,x29,x30,x18. The first slope/square/fused chain does not touch that bank. The second slope consumes it. Original LR is saved on entry and restored on return; no BL occurs inside these inlined loops. No caller-visible scratch or source point is removed.

Chained results substitute producer output registers for the next kernel's left-operand loads. Scratch32..63 stores and loads are dead after inverse caching and are absent. The forward absolute difference remains in x11..x14; zero testing changes only x9,x10 and flags. Multiplication temporaries are renamed away from that live bank, and its two right operand loads become direct reads. The initial Y load is postindexed by64 and its upper limb load uses offset-48; the trailing pointer ADD is omitted. On failure the advanced pointer is clobberable and no next iteration executes. Branch labels are uniquely renamed. Arithmetic branch predicates are preserved because address substitution/renaming adds no flag-setting instruction. The witness checks all structural preconditions and exact generated source equality. Random-register/address differential tests compare retained outputs to independent original canonical loops with separately inverted signed prefixes, and check source memory, dead-scratch canaries, SP and LR.

## 12. Complete ledger and finishing deliveries

The unchanged finishing-sweep row includes16*2^70*2^35 delivery/inflight instructions in addition to sender work:34*SIDE+6 per node. The eight-instruction per-leaf bookkeeping reserve remains. All cold generation, failed caps, sampling, transport, grouping, scanner, resolution, verification and barriers remain in success_probe.py. Total109110449194024512148123786083622023331840 instructions gives log2(total/1000)=126.359056876753724... . A64-term rational atanh series for log, with its positive geometric remainder, and an independent interval for ln2 certify the outward ceiling126359056876754/1000000000000. No prior validity judgment is a premise. The old reviewed packet is unchanged; this source requires its own external review.

## Paired-comparison64-record range scan

For consecutive keys A,B,C, CMP A,B; CCMP B,C,#4,NE; B.EQ branches exactly when A=B or B=C. The slow path appends each equal adjacent pair separately with its own capacity check, in the same order as the original scanner. No relation or success event changes.

Three key registers rotate and four pointer-register pairs rotate. Initial copies put the range's first key and address in every possible suffix-start register. The pair remainder r=((n-1)//2) mod32 selects entry at pair32-r (r=0 selects the end test), using five TBZs and one additional branch when r is odd. Each pair uses the previous key and record pointer, loads two new pointers and their keys, and retains the second key/address for the next pair. All potential entry registers initially hold the same first key/address. Since32 mod3=2 and32 mod4=0, the final key is x7 and final pointer x16, so MOV x5,x7 normalizes the next full block. The hot block costs32*(LDP+2LDR+CMP+CCMP+B.EQ)=192 instructions plus MOV,CMP,B.LO:195 per64 records.

Write n-1=64q+2r+epsilon with0<=r<=31 and epsilon in{0,1}. A nonmatching nonempty range costs25+6epsilon+(r mod2)+6r+195q including range control. Subtracting195n/64 gives1405/64+189epsilon/64-3r/32+(r mod2). Its maximum is1652/64=25.8125 at epsilon=1,r=1, bounded by26. Empty ranges cost5. The slow path adds7 instructions for only the first equality,8 for only the second,11 for both; an odd-prefix hit adds4. Charge8 per listed candidate plus one additional8 for a full-list append attempt. Every append has its own capacity check and preserves original adjacent-pair ordering, even when both equalities hold. Global guard and RET cost2.

Thus a node uses at most195n/64+26*2^48+8*(2^28+1)+2 scanner instructions. This rational inequality holds for each node separately, so summing yields195S/64 without an uncharged ceiling. Production S is divisible by64. The remainder suffix consumes exactly2r records after the optional odd comparison; every full block consumes64, and the exclusive end test prevents out-of-range loads. Empty ranges read no record, and no comparison crosses a range boundary. Candidate identities, order, capacity-abort behavior, and the original partial-match success/cost analysis are unchanged.

The implementation is transferred unchanged from our submitted secp package b00400fa-b4b3-440c-accb-657046f12efa; the record-pointer/end-list ABI is identical here. scanner_check.py exercises all remainder lengths, multiple blocks, empty ranges, both equality paths and caps, including4096-record cases. routing_check.py independently models every executed instruction and compares packet/candidate order across all16 nodes. Finite checks supplement the universal count and semantic argument above.


## 13. Compact grid and explicit bin-mass refinement

The layout uses SIDE=2^35 and N=2^70, with q1=377461702 and q2=2^28. Metadata allocations are35+29 and35+1+28 bits, each exactly64. The final l increment may wrap after the last pair; no emitted index wraps. Production-width boundary tests check maximum row, k and l encodings, plus the post-emission wrap.

The96 bucket bits are[0,35),[64,99),[128,154), so the shared home/bucket width remains96 while home width is70. The48 d4 bits keep the intruder grouping width144. Complete matched width is182, and h/g/v use that changed width; they were not reused at184. Both slot means and bucket means use only the explicitly proved refined masks. Bucket slack is2^30, CAP_3=4099333281, and each bucket's stride is2^36. The TST full-bucket mask changes consistently to36 low bits. The exact inequalities for the two receive regions, bucket-plus-pointer reuse and all fixed storage are asserted in success_probe.py. All changed metadata and stride operations retain the same instruction counts.

FOURIER_CAPACITY.md supplies the full mask Fourier argument with an explicit mixed character-sum constant and rationally rounded bounds. Unlike assuming uniform hash bins, it bounds every assignment of the specified bits. It does not sharpen the144/182-bit intruder/matching bounds. Source constants and design effort are public; no per-target searched advice is introduced.

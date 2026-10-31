# P-256 paired stored-leaf birthday search

## 1. Target, representation and candidate claim

The target is exactly the protected P-256 group, including a=-3, its specified
prime p, prime order r, cofactor one, generator P and uniform Q=[x]P. Output is
the canonical scalar x. Q=infinity returns zero. No auxiliary powers or
secp256k1 endomorphism are used. Let R=2^256 mod p. Store canonical Montgomery
coordinates, and use K(Z)=R*x(Z) mod p for a nonidentity point Z. This is an
exact key for the two-element class {Z,-Z}: the curve equation gives at most
two points for each x, and no nonidentity point has y=0 because r is odd.
Multiplication by R permutes the field, so a prescribed pattern in f of the
256 key bits has at most 2*2^(256-f) preimages in the group.

The complete claim is total classical
work below 2^126.527881 units at the protected C=1000 scalar instructions/unit,
success at least 0.39, peak memory below 2^136 bytes, preprocessing below
2^118 units included in the total, and zero nonuniform advice. This is an
astronomical-memory analytical attack, not a feasible machine demonstration.

## 2. Random construction and pair arithmetic

There are N=2^72 nodes, indexed (i,j) with each coordinate in [0,2^36).
At each node draw b and q1=188761551 independent a_k uniformly modulo r,
and q2=2^27 independent c_l uniformly in {1,...,r-1}. Coins and draws at
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

For one batch fix X and write d_l=x(Y_l)-x(X). The forward loop stores d_l
and prefix product A_l=product_(j<=l) d_j, with A_-1=1. Invert A_(q2-1)
once. The backward loop maintains B_l=1/A_l, computes i_l=B_l*A_(l-1), then
updates B_(l-1)=B_l*d_l. Thus i_l=1/d_l. The two slopes are
(y(Y_l)-y(X))*i_l and (y(Y_l)+y(X))*i_l; square each and subtract x(X)+x(Y_l).
The second is the slope for X-Y up to sign, which has the same square.
These are exactly the two required x-coordinates. No a=0 identity is used.

The supplied field kernels return canonical residues. The assembly forward
loop tests every d_l for zero before multiplying it. The backward loop
overwrites the now-dead denominator/prefix entry with the two output x keys.
A 32-byte Montgomery-one prefix before the batch handles l=0 without reading
outside the allocation. The reference point formulas and executed native-word
interpreter tests are independent checks, not substitutes for these invariants.

## 3. Field arithmetic and instruction-count obligations

The field kernels are extracted from the supplied P-256 source header of
submission e0d1ed74-946e-4f8d-a5ad-a6d219ad3697. They are general canonical
field kernels, not the source's heuristic branchy main loop. Their bodies have
115,91,22,18 instructions for Montgomery multiply, square, add and subtract;
each separate return adds one. Calling with three argument-setting instructions
and BL therefore costs at most 120,96,27,23 respectively. Squaring actually
needs only two argument settings and is charged conservatively at 96 here.

The arithmetic basis for the multiply is radix B=2^64 Montgomery reduction.
The limbs of p are (B-1,2^32-1,0,B-2^32+1), so the Montgomery digit is simply
the current low word t. One reduction step replaces T by (T+t*p)/B. Its low
word cancels. If t=t_low+2^32*t_high, the first two remaining corrections are
t*2^32 (split at B), while the high correction is t*(B-2^32+1), whose low
word is t-(t<<32) modulo B and high word is t-t_high minus the subtraction
borrow. These are the shifts and SUBS/SBC chains in the code. Interleaving the
four schoolbook rows with these four reductions preserves exact integer carry
chains. The final value is congruent to a*b*R^-1 modulo p and below 2p for
0<=a,b<p. The final five-word subtract and CSEL select T-p exactly when T>=p.
The square uses the diagonal terms and twice every cross product, followed by
the same reduction identity. The add and subtract kernels retain the carry or
borrow and select the canonical residue. The listing implements the four schoolbook rows and these four reductions
with explicit carry chains; the word identities above justify their exact
integer arithmetic. Finite random tests are additional evidence, not the
universal argument.

The fixed inverse code performs the public binary exponent p-2: 255 squares
and 127 multiplies, starting from its input. The result is the inverse for
nonzero input by Fermat's theorem. The complete interpreted instruction count
is 39479, including its prologue, copies, calls and return. The batch invariant
ensures its input is nonzero. No reference Python inverse remains in the
implemented batch path; independent reference inverses are used only in tests.

The complete forward and backward loop counts are 144*q2+11 and 744*q2+9.
These include argument preparation, every field instruction, memory access,
zero test and loop branch. The key emitter costs at most 29*q2+11, including
all its stores into the first-hop port; a dropped leaf saves four instructions.
The ledger charges 917*q2+39479+1024 per batch, where the additional 1024
covers the three fixed loop terms, driver call/argument setup, batch addressing
and status handling. This is not a per-step probability discount.

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

Keys split into d1=K0 mod 2^36, d2=K1 mod 2^36, d3=K2 mod 2^24 and
d4=K3 mod 2^48. The 32-byte record is (K2,K3,w2,w3), where
w2=i|(k<<36) and w3=d2|(sign<<36)|(l<<37). Here k<2^28 and l<2^27, so
both metadata words fit. The source column is inferred from the second-hop
sender slot. No source coordinate, sign or scalar index is omitted.

First hop sends into the sender-row slot at node (d1,j). Second hop streams
each source slot and sends into the sender-column slot at node (d1,d2).
There are 2^36 slots per node, each CAP_P=2^22-2^12 records. Source cursor
initialization, clamped final-cursor messages and barriers are charged. The
second-hop cursor is restored to its end on a drop, so it cannot accumulate
unbounded overflow. First-hop attempts number at most n0 per sender; their
32-byte increments cannot wrap a 64-bit cursor at the specified addresses.

At each home, append (record address,d4) into its d3 bucket, or drop at the
cap. Bucket cursors restore the increment on a drop. Count d4, prefix-sum its
counts, and stably scatter pointers by d4 in d3 order. The final d4 counters
are exclusive range ends. Scan K2-adjacent pairs within each d4 range, using
the supplied sixteen-comparison scanner and branch-tree prefix. Its universal
count is (29/8)*records+24*ranges+4*(CAND_CAP+1)+2. It never crosses a range.
Store at most CAND_CAP=2^28 pairs, then stop on the next attempted append.

Consequently matched keys fix 72 home bits, 48 d4 bits and all 64 K2 bits:
184 distinct bits. The d3 bits are already within K2. Partial matches are
not treated as true group collisions. Each is resolved and verified.

Resolution has the same four count-exchange/data subphases as the source.
A home sends requests to the two source nodes through the appropriate relay,
and obtains at most eight 32-byte reply records per leaf. These can contain
X,Y,a,b,c, all totaling at most 224 bytes. It reconstructs L=X+sY and L',
tests L'=uL for u in {+1,-1}, and if n'-u*n!=0 solves
x=(u*m-m')/(n'-u*n) modulo r. It verifies [x]P=Q before canonical output.
Source-request and relay capacities are 2^29 and 2^30. Count exchanges include
allocation messages, replies, payloads and barriers. Their per-node charge
2^41 is retained: 12*2^36+4*(14+8)*2^32<2^41, with eight reply records
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

All field arithmetic is exact; there are no approximate canonicalization,
window or correction-meter losses. A batch defect probability is at most
d_batch=q2*2/(r-1)+1/r (zero denominators or identity X). Excluding any
defective batch in the whole run costs at most N*q1*d_batch. The event is
explicitly charged even though a defective batch's cost is no larger.

For fixed key bits define KB(f)=2*2^(256-f)/(r-1), and rho=(r/(r-1))^4.
Condition on a witness pair's four generators. Free leaves obey KB(f).
At most four leaves are fixed, but retain the conservative reservation 12.
For slot mean m_p=n0*KB(36)+1 and bucket mean m_3=S*KB(96)+1, use
v=m+1+60(q1+q2)KB(f). Then d_p=v_p/(CAP_P-12-m_p)^2 and
d_3=v_3/(CAP_3-12-m_3)^2, with CAP_3=ceil(m_3)+2^32.

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
then yields the displayed conditional drop probabilities. No independent
hashing distribution is assumed.

Retain the source's conditional pair-load bounds, with the justified
two-point fiber KB and 184 matched bits:
 h=rho*(S^2*KB(72)*KB(184)+(2S+20)*KB(184)),
 g=2*rho*n0*S*KB(184),
 v=18*rho*((S/2^36)*S*KB(36)*KB(184)+(2S+20)*KB(184)).
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

The intruder rate is rho*(S+12)*KB(144), not an incorrectly smaller 184-bit
rate: d1,d2,d3,d4 group only 144 bits before the adjacency scan. If no intruder
enters a witness class's run and the capacity events do not occur, the useful
adjacent pair is resolved. False partial-key matches can exhaust caps; their
work and losses are exactly what the pair-load bounds cover.

The explicit lower-bound expression is
1-exp(-mu)-CS
-mu*(1+1/r)*(1/r+rho*(S+12)*KB(144)
             +2*(d_batch+2d_p+d_3+d_cl+d_rq+d_rl))
-N*q1*d_batch-N*(q1+q2+1)*2^-248.
The factor 1+1/r bounds (2/r)/P1. All quantities except exp are exact rational
functions. Round mu down to 200 binary fractional bits and upper-bound its
negative exponential by the even 40-term Taylor polynomial. This is an upper
bound by the sign of the Taylor remainder. The resulting rational success
bound exceeds 0.390100001 at q1=188761551. The declared success is conservatively 0.39. None of this identifies reviewer confidence with algorithm success.

## 7. Work and memory

success_probe.py contains the complete integer work rows, including
both routing hops, every delivered record, all empty tables and range ends,
preprocessing, capped sampling, phase-4 work, resolution, verification, output
and an additional eight instructions per leaf for bookkeeping. Work is summed
over processors, never amortized across future targets.

A deliberately simple affine binary scalar multiplication suffices for cold
generation and verification. At most 512 affine point operations are used.
Each costs below 2^16 instructions using the exact field inverse and kernels,
with identity/equal/opposite-point cases branched explicitly. Thus including
control it costs below 2^26. Charge 2^27 per generated point/header. A resolver
fetches stored X,Y,a,b,c rather than recomputing their scalar multiples. Two
affine reconstructions, scalar arithmetic, one scalar inversion modulo r, one
verification [x]P and its messages fit within 2^27 instructions. For the scalar
inverse, bitwise long reduction costs below 2^16 per 256-bit modular product,
and binary exponentiation uses at most 512 products, hence below 2^25. Together
with the verification bound and remaining small work this is below 2^27.

The two receiving regions have length 32*CAP_P*2^36 each. Place the first
at 2^52 and the second immediately afterward. Their joint end is below 2^64.
All fixed arrays, code, stack, source tables, batching scratch, cursor/end
tables, buffers and link state fit below 2^52. From the grouping phase onward,
the first receive region is reused for 2^24 bucket strides of 2^38 bytes plus
the pointer array; CAP_3*16<2^38 and their combined size fits that region.
Every node has a 64-bit address space; total memory is below N*2^64=2^136.
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
Grouping has 10+2/8 for bucketing, 4+2/8 for stable scattering and 29/8 for
scanning, totaling 145/8 per record. Bucket failures cost no more.

Two cursor tables cost 3 instructions per entry plus returns. Two final-cursor
sweeps cost 9 sender instructions plus 8 delivery/in-flight instructions per
destination, plus fixed entries: 34*SIDE+6 per node for both sweeps. Grouping has 12 per source
slot, 4 fixed, 4 per d3 bucket to initialize and 11 per d3 bucket to scatter,
and 4+10 instructions per four d4 counters for zeroing and prefix sums.
The additional scanner term is 24 per d4 range, four per retained or
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
has the conservative 2^27 allowance justified in Section 7, even for false
partial matches and failed verification. All maxima are charged whether or
not their paths occur.

Each 64-bit scalar modular product for the final solve can be implemented
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
single AArch64 encodings: zero/one/powers of two, 32 or 36 contiguous bits,
and the cyclic run of bits in 0xffffffff00000001. Branch displacements are
within the source blocks. Public parameters and the fixed inverse exponent
are program constants; delivering and initializing the program is charged.
The emission tag's final increment can wrap after the last index; it is not
read again. Every emitted index and every memory address fits its stated
word or region. The first-hop cursor bound includes all dropped attempts;
the second-hop and bucket cursors clamp/undo on drops.

The exact success arithmetic, integer ledger sum, bound comparison and layout
assertions are in success_probe.py. The score comparison uses the exact
integer inequality R_total^1000000 <= 1000^1000000 * 2^126527881, not a rounded floating
logarithm. The preprocessing ceiling is checked on a conservative superset
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

## Inlined arithmetic revision and explicit finishing deliveries

This version derives from submission 32e6d312-830a-4a46-b667-029308a7d201,
commit 844287b3e35224743d6326c6e6535421e12df28c. It retains its full mathematical
argument, inputs, probabilities, capacities and routing. No previous judgment
is used as a replacement for that argument.

For each field call in pair_batch.s, replace its input/output pointer setup
and BL/RET pair with the unchanged field body. Only the effective addresses
of LDP/STP are specialized to the stable caller registers x20..x23. All
addresses are identical to the original kernel's addresses, fit signed pair
load/store offsets, and are multiples of eight. The field bodies never alter
x19..x23. Arithmetic instructions, carry chains and canonicalization are
unchanged. The original x0..x2 argument values are dead outside the called
bodies; temporary uses of those registers inside bodies are retained. No
value live across a call, loop counter, scratch alias, or branch condition
is removed. Thus each inlined body has the same memory effect as its call.
The only removed state changes are dead pointers and link-register updates;
the batch still preserves its external return address.

Two binary field calls in the forward pass save ten instructions per item.
The backward pass has nine binary and two unary calls, saving 9*5+2*4=53
instructions per item. The paired batch plus key emission therefore costs
917*q2 rather than 980*q2 instructions, with the existing fixed allowance
retained. Tests assert both complete loop-count formulas and compare every
output with independent field/curve calculations. Program delivery remains
covered by the existing 2^20 instructions per node and published memory
allowances: the inlined listing is well below that budget.

The prior reviewer noted that finishing sweeps omitted delivery/in-flight
work. This ledger explicitly adds 16*2^72*2^36 = 2^112 instructions to that
row, changing 18*SIDE+6 to 34*SIDE+6 per node. This correction is not hidden
in reserves. The eight-instruction per-leaf reserve remains untouched.

The exact ledger is R_total = 122655565351338577669916552428608091586560
scalar instructions, giving log2(R_total/1000) = 126.5278803016831... .
The claimed upper bound is 126.527881 = 126527881/1000000. The integer-power comparison
R_total^1000000 <= 1000^1000000 * 2^126527881 certifies the ceiling without floating-point
rounding. The success expression remains above 0.390100001. No old total
or prior validity judgment supplies this calculation.

The retained pair_batch_reference.s is the prior submitted call-based loop.
inline_witness.py deterministically expands all thirteen field calls, verifies
all effective-address substitutions and offset encodings, and compares the
resulting instruction stream to the new listing. It checks stable-register
noninterference and the exact 10/53 per-item savings. Its finite structural
check supports the semantic argument above; it is not a substitute for it.

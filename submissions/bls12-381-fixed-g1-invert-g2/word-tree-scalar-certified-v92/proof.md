# Random-stride baby/giant inversion with exact trace recurrence

Proposed time_log2=15568/125=124.544. This is an analytic claim, not an observed score, runtime pass, qualification or promotion. It uses the full protected p381/r255, actual nonidentity supplied P, uniform full GT target z including identity, threshold0.39 and fixed reference C20000. There are still six production tracks and one computation metric. Complete cold work is below2^124, peak storage below2^138 bytes, advice0. All tables and powers are constructed for this attack. No auxiliary powers are given and no work is amortized over future targets. Explicit finite randomness premise R is disclosed in randomness.md.

Validate canonical inputs, curves and subgroups. Set h=e(P,Q0) for the public G2 generator Q0. It generates GT for every valid nonidentity P. If z=1, return and verify G2 identity. Otherwise the target has unique z=h^a; a is uniform Fr when z is uniform GT. The probability proof below is over uniform a and fresh independent algorithm coins, as required, rather than a worst-case-target assertion.

Sample one uniform b in Fr and twelve independent Fm61 affine coefficients and seventy-two independent raw64-bit NH keys by the capped finite sampler. Compute d=h^(-b) with actual GT exponentiation. Let
m=44951011978550636487095639346140676096, L=2m=89902023957101272974191278692281352192.
m is a multiple of3*2^110, hence2048 and512. Baby records correspond to h^j for j=1,...,m; giant records to z*d^i=z*h^(-bi) for i=0,...,m-1. Stored scalar is j for a baby, bi modulo r for a giant. Both workers use the SAME independently sampled fingerprint coefficients.

## Exact recurrence and emitted records

Let Phi(X)=a0+b1*t=trace_(Fp12/Fp4)(X)/3, t=w^3,t^2=1+u,u^2=-1. This is four native canonical Fp coordinates. For any GT increment g, its three p^4-conjugates satisfy
Y^3-cY^2+c^(p^2)Y-1, c=3Phi(g).
Consequently v_n=Phi(Z*g^n) satisfies
v_n=c*v_(n-1)-c^(p^2)*v_(n-2)+v_(n-3).
The recurrence follows separately for each conjugate root and therefore also for their Z-weighted sum. It remains valid for g=1 with repeated roots.

Baby increment is h, with cold initial Phi(1),Phi(h),Phi(h^2). Giant increment is d, with cold initial Phi(z),Phi(zd),Phi(zd^2). These are actual computed GT values. Write c=A+B*t, previous values X=a+b*t and Y=d0+e*t, with all letters except t in Fp2. The first two recurrence terms factor EXACTLY as
cX-c^(p^2)Y = A(a-d0)+B*xi*(b+e) + t*(A(b-e)+B(a+d0)).
Eight canonical Fp add/sub operations form four Fp2 inputs. A two-block vector lookup then computes the two output Fp2 components. The integer accumulators are seeded with the appropriate third-predecessor coordinates, so the final canonical output directly REPLACES that third raw history value. No later copy or key canonicalization writes into the raw history.

recurrence_fp2_step.S is2711 throughRET:250 fused-wrapper own instructions plus the new2461-instruction map. The new map reads fully constructed ordinary q*p products for q0..49 instead of multiplying them, with the exact31-instruction correction in quotient_certificate.md. Caller16 plus four pointer rotations give2731 retained2752. All full-field recurrence, actual characteristic c=3Phi(increment), temporary layout, raw-history and key-separation proofs remain unchanged.

Use the exact12-Frobenius-fiber key S192 from key_contract.md and the actual canonical_twelve_bfirst.S. All emitted conceptual elements are actual GT by the recurrence proof. Abort if a giant identity is encountered; baby identities cannot occur since1<=j<=m<r. The constant number of cold initial records is charged, and also conservatively assigned the normal hot allowance. Identity targets are solved directly.

The universal hash returns the same full160 bits in bytes0..19 of U32,with upper96 initially zero. Set color at20. Payload remains(S192,scalar32,U32)=256 bytes. Metadata is now(U24,ordinal128)=40 bytes:color20,temporary occupancy21 initially zero,padding22..23,ALL128 ordinal bits24..39. Baby ordinal j-1 and giant m+i remain complete. Sorting/group equality uses ONLY first20 bytes. Occupancy21 is restored before output and never overlaps the fingerprint,color or ordinal. metadata40_certificate.md supplies exact source/bank/confirmation transfer proofs.

Generate all L records using the finite worker streams in architecture.md. Sort metadata by the complete two stable80-bit passes in radix80_certificate.md. Enumerate all unordered equal-U pairs, aborting BEFORE more than floor(L/8192) confirmations. Fetch full S, scalar and color. Reject unequal S and baby/baby pairs. A mixed pair is always informative: try all12 phases lambda=p^ell,ell0..11, with a=bi+lambda*j. Do NOT reject a mixed pair merely because its two stored scalars happen to be equal. For two giants, reject equal scalars; otherwise try the11 nonidentity phases:
a=(b_t-lambda*b_s)/(1-lambda),
because a-b_t=lambda*(a-b_s). Phase1 cannot be the true phase for distinct scalars. For the FIRST informative exact-S pair, some phase is correct. Verify h^a=z, construct[a]Q0, and verify canonical G2/curve/subgroup/full e(P,Q)=z before return. At most12 terminal phase attempts occur in total, all cold-paid.

## Exact useful-collision probability

Write u0=0xd201000000010000 and lambda0=u0^2. The full protected parameters satisfy r=u0^4-u0^2+1, p=-u0 modr. At the displayed m,
m<u0^2-1, 2m^2<r, 3m<lambda0, 4m(lambda0+1)<r.
These are exact integer inequalities, not toy parameters.

The union of baby Frobenius exponents
Sexp={p^ell*j modr:0<=ell<12,1<=j<=m}
has size s=12*(m-floor(m/u0)). Here is a complete short-interval proof. The twelve roots are ±1,±u0,±u0^2,±u0^3,±(u0^2-1),±(u0^3-u0). The ±u0^2 and ±(u0^2-1) products lie outside[1,m] without modular wrap because m<u0^2-1 and m*u0^2<r-m. Of ±u0 only +u0 produces pairs j'=u0*j, exactly floor(m/u0). The inverse phases ±(u0^3-u0)=∓u0^(-1) give only the same reverse pairs: multiply by u0 and use u0*m+m<r. Neither -1 nor -u0 produces another small positive value.

For ±u0^3 put x=-u0. The congruence j'=x^3*j modr has integer lattice basis (x,x^2-1),(x^2-1,-x), determinant-r. Both basis vectors have squared norm r and are orthogonal, so every nonzero lattice vector has norm at least sqrt(r). The basis spans the full congruence lattice because its determinant is its index r. A pair with |j|,|j'|<=m has squared norm<=2m^2<r, impossible. Negating j' covers the other sign. Finally u0^2>m prevents any three-point +u0 chain. Thus each duplicate baby orbit contains exactly two baby exponents and there are precisely floor(m/u0) such pairs. All nonidentity orbits have size12, proving s.

For every i=0,...,m-1 and c in Sexp, select a mixed event a=i*b+c. There are M=m*s distinct events. Each has probability1/r over uniform independent a,b. Two different slopes i give a full-rank2x2 system and probability1/r^2; same slope with different nonzero c is disjoint. Each event gives an actual useful mixed collision.

Select FOUR giant/giant phase families. Put lambda=lambda0, so lambda^2-lambda+1=0 modulo r. For each0<=i<k<m use both sigma=lambda and lambda^(-1)=1-lambda. Also use sigma=-lambda and -lambda^(-1)=lambda-1 ONLY when i and k have different residues modulo3. In all cases define
c_sigma(i,k)=(i-sigma*k)/(1-sigma),
and the guarded event b!=0, a=b*c_sigma(i,k).

Here is a complete injectivity/disjointness certificate in ordinary integer coordinates:
sigma=lambda: c=k+lambda*(i-k);
sigma=lambda^(-1): c=i+lambda*(k-i);
sigma=-lambda: 3c=(2i+k)+lambda*(k-i);
sigma=-lambda^(-1): 3c=(i+2k)+lambda*(i-k).
These identities follow by lambda^2=lambda-1 and are equations in Fr. Their coefficient pairs are NOT treated as arbitrary-width machine words; they are only the finite success proof.

For the first two families, the lambda coefficient has opposite nonzero signs; constants lie in[0,m). Two such representations can be equal modulo r only when their integer coefficients agree: |Delta_constant|<m<lambda and the whole integer difference has absolute value<2m(lambda+1)<r. Thus each is injective and they are disjoint. For either last family, the scaled constant lies in[0,3m), the lambda coefficient has magnitude<m, and the pair determines i,k uniquely (the integer linear transform has determinant±3). Equality within or across the last families similarly forces equality of their integer coordinates, using3m<lambda and4m(lambda+1)<r. Their lambda signs are opposite, so the last families are disjoint.

To compare a last-family slope with a first-family slope, multiply the latter by3. Constant-coordinate differences are less than3m; lambda-coordinate differences less than4m. The whole integer difference has absolute value<4m(lambda+1)<r, and3m<lambda again forces both coordinates to agree. But that would make the last family's lambda coefficient divisible by3, contrary to the selection i!=k mod3. Hence these families are disjoint from the first two. Each slope is nonzero: its lambda coefficient is nonzero, dominates the bounded constant if it has opposite sign, and the full absolute value is below r.

Since m is divisible by3, the number of pairs i<k in different residue classes is exactly m^2/3. Therefore the number of DISTINCT nonzero selected giant slopes is
G=2choose(m,2)+2*m^2/3=5*m^2/3-m.
The four phases are p^2,p^10,p^8,p^4; all are tested by the algorithm, with inverse orientation if needed. The actual algorithm still tries all useful phases, so omitting duplicate slopes only affects this lower-bound proof.

Each giant event has probability(r-1)/r^2 because b!=0. Different giant slopes are disjoint under that guard. Their points are nonidentity, since c_sigma-i=sigma*(i-k)/(1-sigma)!=0. Scalars bi,bk are distinct. A mixed/giant event pair is either an independent linear system with probability at most1/r^2 or parallel and disjoint, because the mixed intercept c is nonzero. These are exact two-variable finite-field facts, not collision independence or a hash heuristic.

Let
mu=M/r+G*(r-1)/r^2,
U2=(choose(M,2)+M*G)/r^2.
The omitted giant/giant intersections are exactly zero. The second Bonferroni inequality gives useful-event probability at least mu-U2. Here
s=539412143742607637809501443445415278032,
M=24247121734749633684056342593153827383372251374231666912243654457484296323072,
G=3367655796493004678563705089834660441163038088766687439445668254313835659264.
All selected events correspond to actual records; inverse phase orientation is harmless because the algorithm tests every required Frobenius phase.

For confirmations, distinct giant exponent pairs are uniform Fr^2 as functions of uniform a,b; a mixed giant exponent is uniform Fr. Their expected orbit matches are bounded by12/r per pair. The deterministic baby duplicates contribute exactly floor(m/u0). With independent fingerprint collision bound delta=(2^-63+128/q)^2*(2^-63+512/q),q=2^61-1 (proved in nh_certificate.md),
caploss=(choose(L,2)*(12/r+delta)+floor(m/u0))/(floor(L/8192)+1).
Any giant identity contributes at most m/r by a union bound. Capped sampling/extraction contributes less than2^-970. Identity targets return successfully; their unused coins can be extended hypothetically for this unconditional probability proof.

EXACT finite rational evaluation gives
mu-U2-caploss-m/r-2^-970 >3900031/10000000>0.39.
Diagnostics: mu=0.5266390126800469, U2=0.13661194642012903, caploss=2.5306782098113685e-07, final=0.3900268131920969. Decimal rounding is not the certificate.

## Complete work and storage

Unchanged recurrence2711/key813/NH375 retain2752/824/384. Baby scalar7+caller3<=12,giant full-Fr21+caller5<=28 give total<=12m+28m=20L for every execution. With main704,generation4684 plus two local passes720,scan128,confirmation1L and spare256 gives base5789. Metadata40 network1027 giveshot6816,retained6896. scalar_stream_certificate.md proves both full32-byte scalar encodings,index/count bounds,all carries and early-failure costs. All original field/GT/probability parameters and event formulas remain unchanged.


Each confirmation is at most8192 instructions, including four RPCs, all four Root formations256, full comparisons and losing branches. The cap floor(L/8192) therefore costs at most L, already included. All12 terminal verification attempts are cold. The more aggressive cap is justified by the displayed complete failure probability; it does not change the0.39 success threshold.

All ordinary processor/code/page/bootstrap work is below2^121;edge/flag/controller setup below2^118. The two80-bit sort's ALL2^80 counter/pending construction,all prefixes,empty finish steps,cache tails and initial/final exceptions are below2^104. Actual fields,tables,RNG,powers,constants,validation and12 terminal phases retain their fixed bounds. Sum complete cold<2^124. Generic data capacity2^136,private2^119,edges2^111 and algorithm-specific data<2^93 give peak<2^138. Unused capacity supplies no data. Every native pointer is below2^47;96-bit identifiers and128-bit counters are ordinary words. N=6896L+2^124 satisfies EXACT N^125<20000^125*2^15568,so proposed15568/125=124.544 under unchangedC20000. This is analytic work,not production hardware qualification.

The key uses a specified fixed total order: first compare the twelve little-endian64-bit words of signed b1 in increasing word index, then a0 as a secondary coordinate. On nonidentity GT, the T²=-1 lemma in key_contract.md proves the two candidate b1 values are distinct, so the secondary comparison is never needed. The stored S192 still includes full a0 followed by full b1. Any fixed total order of the two even-coset keys preserves exact twelve-orbit fibers; both workers use this identical order. Output is disjoint from raw history. Universal hashing applies to arbitrary distinct192-byte strings, so the unchanged full success proof applies to this representation.

The exact new fingerprint construction, collision composition, full375-instruction bound and updated coin declaration are in nh_certificate.md/randomness.md. It replaces the prior direct26-chunk affine family completely; no old81-coefficient claim is used.

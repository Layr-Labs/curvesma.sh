# Random-stride baby/giant inversion with exact trace recurrence

Proposed time_log2=31121/250=124.484. This is an analytic claim, not an observed score, runtime pass, qualification or promotion. It uses the full protected p381/r255, actual nonidentity supplied P, uniform full GT target z including identity, threshold0.39 and fixed reference C20000. There are still six production tracks and one computation metric. Complete cold work is below2^124, peak storage below2^138 bytes, advice0. All tables and powers are constructed for this attack. No auxiliary powers are given and no work is amortized over future targets. Explicit finite randomness premise R is disclosed in randomness.md.

Validate canonical inputs, curves and subgroups. Set h=e(P,Q0) for the public G2 generator Q0. It generates GT for every valid nonidentity P. If z=1, return and verify G2 identity. Otherwise the target has unique z=h^a; a is uniform Fr when z is uniform GT. The probability proof below is over uniform a and fresh independent algorithm coins, as required, rather than a worst-case-target assertion.

Sample one uniform b in Fr and twelve independent Fm61 affine coefficients and seventy-two independent raw64-bit NH keys by the capped finite sampler. Compute d=h^(-b) with actual GT exponentiation. Let
nb=35597089187900144514297950209050673152, ng=53124985308099088881309773192415412224, L=nb+ng=88722074495999233395607723401466085376.
nb=9141*(3*2^110),ng=13642*(3*2^110),so BOTH are multiples2048/512/3. They are separately chosen numbers of attack records,not reduced curve/group parameters. Baby records are h^j forj1..nb;giantsz*d^i fori0..ng-1. StorefullFrscalarj/bi. Babyordinalj-1,giantordinalnb+i givesallelements0..L-1 uniquely. SharedhashkeysindependentofbothcompleteGTstreams.

## Exact recurrence and emitted records

Let Phi(X)=a0+b1*t=trace_(Fp12/Fp4)(X)/3, t=w^3,t^2=1+u,u^2=-1. This is four native canonical Fp coordinates. For any GT increment g, its three p^4-conjugates satisfy
Y^3-cY^2+c^(p^2)Y-1, c=3Phi(g).
Consequently v_n=Phi(Z*g^n) satisfies
v_n=c*v_(n-1)-c^(p^2)*v_(n-2)+v_(n-3).
The recurrence follows separately for each conjugate root and therefore also for their Z-weighted sum. It remains valid for g=1 with repeated roots.

Babyincrementh/giantincrementd andactualcoldinitialPhi(1/h/h²),Phi(z/zd/zd²)areunchanged. Writec=A+B*t,X=a+b*t,Y=d0+e*t. The EXACTrecurrencefactorizationiscX-c^(p²)Y=A(a-d0)+Bxi(b+e)+t[A(b-e)+B(a+d0)]. The fifthactualBOTHtransformreplacesV1=b+e byxi*V1,soonlytwoA/Bcoefficienttablesareneeded. Real-inputtablesadd(Cr,Ci);imaginaryinputsSUBTRACTCi/ADDCr. All381bitsuseeleven35/34-bitdigits. Two27TiBtables8..62,fixed13.5TiB64..77.5,parent4/child80..128arefullydisjoint47-bitaddresses.

The signedrealaccumulatorstartsOLDthird+22p;exactly22negativecanonicaltermsensureNONNEGATIVEatALLprefixes,andbothreal/imag<45p<49p. Existingactualq*p slots/fullquotientcorrectionremainvalid. Fullfp2_certificate.mdprovesallbyte/table/digit/alias/stack/constructor/countbounds. InlinedROWPAIRcore2353+wrapperportion286=2639,caller20=2659retained2664. Everyrawhistoryisimmutableexcepttheexactthirdreplacement;bothseedcoordinatesarereadbeforetheirpairwrites. No keycanonicalizationmutatesrawhistory,andno tableisgivenbyanoracle.

Use the exact12-Frobenius-fiber key S192 from key_contract.md and the actual canonical_twelve_bfirst.S. All emitted conceptual elements are actual GT by the recurrence proof. Abort if a giant identity is encountered; baby identities cannot occur since1<=j<=nb<r. The constant number of cold initial records is charged, and also conservatively assigned the normal hot allowance. Identity targets are solved directly.

Compute full canonical S192 into a private192-byte scratch disjoint from raw history and payload. The revised nh_fingerprint.S stores the THREE FULL canonical Fq affine results Y0,Y1,Y2 as Ctag24 BEFORE low-bit masking for the unchanged U160. It takes only the originalthree arguments;x0 points at the88-byte payload base. Its stores write U32 at0..31 and Ctag24 at64..87,leaving fullscalar32 at32..63 untouched. Color is at20. Metadata40 remains U24 plus ALL128 ordinalbits24..39. FullS192 is discarded only afterhashing. Ctag is a probabilistic stronger tag,NOT an exact orbit key;neither full rawNH nor truncated low160bits is mistaken for its definition. AllS192inputwords are read before any outputstore. Keys are disjoint from payload.

Generate all L records and stably sort metadata in the same two80-bit passes. Enumerate every unordered equal-U pair,aborting BEFORE more than Kc=floor(L/8192) cheap confirmations. Fetch two metadata40 and two payload88 rows;compare full Ctag24 and ignore baby/baby or equal-scalar giant/giant pairs. Mixed pairs remain informative even when scalars equal. Before each informative Ctag-equal pair,abort if the new expensive counter is already Kt=floor(L/2^36);otherwise increment it and try the12 mixed phases a=bi+lambda*j or11 nonidentity giant phases a=(b_t-lambda*b_s)/(1-lambda),lambda=p^ell. Each candidate is checked by ACTUAL dense full-field GT exponentiation h^a==z. If every phase fails,CONTINUE enumeration. A tag collision must NEVER cause immediate failure or unverified output. On GT equality,construct[a]Q0 and verify canonical subgroup/full supplied-P pairing before returning. Since h generates GT,a successful verification is sound even for a false tag collision. A genuine informative exact-S pair always yields equal Ctag and U and has a correct tested phase. Terminal G2 construction/final pairing occur once;the potentially many losing GT powers are hot-capped,not cold-only.

nh_certificate.md proves the extra tag-pair failure bound and2^36 per-pair constructive verification schedule. metadata40_certificate.md supplies the new11-payload/five-metadata complete frame assembly and exact88-byte bank strides. No full-S fetch or equality is performed at confirmation.

## Exact useful-collision probability

Write u0=0xd201000000010000 and lambda0=u0^2. The full protected parameters satisfy r=u0^4-u0^2+1, p=-u0 modr. At the displayed counts,
nb<u0^2-1,2nb^2<r,nb*u0^2<r-nb,u0*nb+nb<r,3ng<lambda0,4ng(lambda0+1)<r.
These are exact integer inequalities, not toy parameters.

The union of baby Frobenius exponents
Sexp={p^ell*j modr:0<=ell<12,1<=j<=nb}
has size s=12*(nb-floor(nb/u0)). Here is a complete short-interval proof. The twelve roots are ±1,±u0,±u0^2,±u0^3,±(u0^2-1),±(u0^3-u0). The ±u0^2 and ±(u0^2-1) products lie outside[1,nb] without modular wrap because nb<u0^2-1 and nb*u0^2<r-nb. Of ±u0 only +u0 produces pairs j'=u0*j, exactly floor(nb/u0). The inverse phases ±(u0^3-u0)=∓u0^(-1) give only the same reverse pairs: multiply by u0 and use u0*nb+nb<r. Neither -1 nor -u0 produces another small positive value.

For ±u0^3 put x=-u0. The congruence j'=x^3*j modr has integer lattice basis (x,x^2-1),(x^2-1,-x), determinant-r. Both basis vectors have squared norm r and are orthogonal, so every nonzero lattice vector has norm at least sqrt(r). The basis spans the full congruence lattice because its determinant is its index r. A pair with |j|,|j'|<=nb has squared norm<=2nb^2<r, impossible. Negating j' covers the other sign. Finally u0^2>nb prevents any three-point +u0 chain. Thus each duplicate baby orbit contains exactly two baby exponents and there are precisely floor(nb/u0) such pairs. All nonidentity orbits have size12, proving s.

For every i=0,...,ng-1 and c in Sexp, select a mixed event a=i*b+c. There are M=ng*s distinct events. Each has probability1/r over uniform independent a,b. Two different slopes i give a full-rank2x2 system and probability1/r^2; same slope with different nonzero c is disjoint. Each event gives an actual useful mixed collision.

Select FOUR giant/giant phase families. Put lambda=lambda0, so lambda^2-lambda+1=0 modulo r. For each0<=i<k<ng use both sigma=lambda and lambda^(-1)=1-lambda. Also use sigma=-lambda and -lambda^(-1)=lambda-1 ONLY when i and k have different residues modulo3. In all cases define
c_sigma(i,k)=(i-sigma*k)/(1-sigma),
and the guarded event b!=0, a=b*c_sigma(i,k).

Here is a complete injectivity/disjointness certificate in ordinary integer coordinates:
sigma=lambda: c=k+lambda*(i-k);
sigma=lambda^(-1): c=i+lambda*(k-i);
sigma=-lambda: 3c=(2i+k)+lambda*(k-i);
sigma=-lambda^(-1): 3c=(i+2k)+lambda*(i-k).
These identities follow by lambda^2=lambda-1 and are equations in Fr. Their coefficient pairs are NOT treated as arbitrary-width machine words; they are only the finite success proof.

For the first two families, the lambda coefficient has opposite nonzero signs; constants lie in[0,ng). Two such representations can be equal modulo r only when their integer coefficients agree: |Delta_constant|<ng<lambda and the whole integer difference has absolute value<2ng(lambda+1)<r. Thus each is injective and they are disjoint. For either last family, the scaled constant lies in[0,3ng), the lambda coefficient has magnitude<ng, and the pair determines i,k uniquely (the integer linear transform has determinant±3). Equality within or across the last families similarly forces equality of their integer coordinates, using3ng<lambda and4ng(lambda+1)<r. Their lambda signs are opposite, so the last families are disjoint.

To compare a last-family slope with a first-family slope, multiply the latter by3. Constant-coordinate differences are less than3ng; lambda-coordinate differences less than4ng. The whole integer difference has absolute value<4ng(lambda+1)<r, and3ng<lambda again forces both coordinates to agree. But that would make the last family's lambda coefficient divisible by3, contrary to the selection i!=k mod3. Hence these families are disjoint from the first two. Each slope is nonzero: its lambda coefficient is nonzero, dominates the bounded constant if it has opposite sign, and the full absolute value is below r.

Since ng is divisible by3, the number of pairs i<k in different residue classes is exactly ng^2/3. Therefore the number of DISTINCT nonzero selected giant slopes is
G=2choose(ng,2)+2*ng^2/3=5*ng^2/3-ng.
The four phases are p^2,p^10,p^8,p^4; all are tested by the algorithm, with inverse orientation if needed. The actual algorithm still tries all useful phases, so omitting duplicate slopes only affects this lower-bound proof.

Each giant event has probability(r-1)/r^2 because b!=0. Different giant slopes are disjoint under that guard. Their points are nonidentity, since c_sigma-i=sigma*(i-k)/(1-sigma)!=0. Scalars bi,bk are distinct. A mixed/giant event pair is either an independent linear system with probability at most1/r^2 or parallel and disjoint, because the mixed intercept c is nonzero. These are exact two-variable finite-field facts, not collision independence or a hash heuristic.

Let
mu=M/r+G*(r-1)/r^2,
U2=(choose(M,2)+M*G)/r^2.
The omitted giant/giant intersections are exactly zero. The second Bonferroni inequality gives useful-event probability at least mu-U2. Here
s=427165070254801734143346850431823707576,
M=22693138081419457251764382371353015431987635445861934120674371349332071809024,
G=4703773439976240075985909219823602141921247766664565908941058835693998964736.
All selected events correspond to actual records; inverse phase orientation is harmless because the algorithm tests every required Frobenius phase.

For confirmations, distinct giant exponent pairs are uniform Fr^2 as functions of uniform a,b; a mixed giant exponent is uniform Fr. Their expected orbit matches are bounded by12/r per pair. The deterministic baby duplicates contribute exactly floor(nb/u0). With independent fingerprint collision bound delta=(2^-63+128/q)^2*(2^-63+512/q),q=2^61-1 (proved in nh_certificate.md),
caploss=(choose(L,2)*(12/r+delta)+floor(nb/u0))/(floor(L/8192)+1).
Any giant identity contributes at most ng/r by a union bound. Capped sampling/extraction contributes less than2^-970. Identity targets return successfully; their unused coins can be extended hypothetically for this unconditional probability proof.

The SAME NH/affine keys supply Ctag and U,so their equality is NOT independent. For one row,Pr[Y(S)=Y(S')]<=2^-63+1/q:NHcollides with probability<=2^-63;conditionalondifferentNH,theinjectivechunk/independentaffinemap gives exactly1/q. Three mutually independentrowkeysets give eta=(2^-63+1/q)^3 for fullCtag24 equality. It is NOT2^-192 or2^-189,and must NOT be multipliedbydelta. For ALL recordpairs,T<=choose(L,2)*(12/r+eta)+floor(nb/u0). This overcounts ignored babypairs safely. Markov gives tagloss=T/(floor(L/2^36)+1). All keysets areindependentofcompleteGTstreams,but neitherdifferentrecordpairs norCtag/Uneedbeindependent. The caploss above is still the cheap-U cap loss. Union bound both cap failures,identity and paid-coin failures WITHOUT multiplying them.

EXACT finite rational evaluation gives
mu-U2-caploss-tagloss-ng/r-2^-970 >3900031/10000000>0.39.
Diagnostics:mu=0.5224841090168715,U2=0.1324713058948799,caploss=2.4974634688818906e-07,tagloss=4.874708844186509e-07,final=0.3900120659047603. Formulas,not rounded decimals,are the certificate. Identity targets return successfully;unused keys/streams can be sampled hypothetically in this unconditional proof.

## Complete work and storage

Recurrence2639/key813 retain2664/824. Revised NH369+caller4=373 retains376. Baby7+caller3<=12 and giant full-r21+caller5<=28 give work<=12nb+28ng<=22L for EVERY execution,including all early abort prefixes. No balance assumption is made. Generation4590 plus two local passes720,scan128,cheap confirmations1L,NEW expensive confirmations1L and spare256 gives base5696. Network is generation595+two count23/scatter128 passes302+scan36=933. Hot6629 retains6704,leaving75 spare. Every failure prefix is included.

Cheap confirmation stays<=8192 including four RPCs,Root formations,full24-byte compare and all losing branches;READ5/READ11 transport3276<=retained3408. Expensive verification costs<=2^36 per informative tag pair including ALL12 phase attempts,so Kt permits<=L extra work. Native scalar arithmetic/control/GT decoding has no oracle access;nh_certificate.md provides a constructive dense GT arithmetic bound. The once-only output G2 construction and full pairing remain in cold2^92. Repeated false-tag phase attempts are NEVER assigned there.

All ordinary processor/code/page/bootstrap work<2^121;edge/flag/controller setup<2^118;ALL2^80 sort construction/prefix/finish/cache exceptions<2^104. Actual fields/tables/RNG/powers/constants and once-only validation/output keep fixed bounds,sum cold<2^124. Smaller payload plus extra192 scratch retains peak<2^138/advice0. Every native pointer fits47bits;96-bit identifiers/128-bit counters are ordinary words. N=6704L+2^124 satisfies EXACT N^250<20000^250*2^31121,proposed31121/250=124.484 under unchangedC20000. No physical production execution or promotion is claimed.

The key uses a specified fixed total order: first compare the twelve little-endian64-bit words of signed b1 in increasing word index, then a0 as a secondary coordinate. On nonidentity GT, the T²=-1 lemma in key_contract.md proves the two candidate b1 values are distinct, so the secondary comparison is never needed. The temporary S192 still includes full a0 followed by full b1. Any fixed total order of the two even-coset keys preserves exact twelve-orbit fibers; both workers use this identical order. Output is disjoint from raw history. Universal hashing applies to arbitrary distinct192-byte strings, so the complete two-cap success proof applies to this representation.

The exact new fingerprint construction, collision composition, full369-instruction bound and updated coin declaration are in nh_certificate.md/randomness.md. It replaces the prior direct26-chunk affine family completely; no old81-coefficient claim is used.

The full-range short affine fold and complete369-instruction native bound are proved in nh_certificate.md. Both outputs remain exactly unchanged on every legal full-width key/input,so all finite probability/continuation arguments are identical.



The complete concrete selective-notification/dispatch construction in radix80_certificate.md replaces the prior unspecified finite-notification premise rejected in public7f7bacd8-cbb9-4cb3-b13d-607d91ef7406. Actualwaiter<=20/fullnotification32perendpoint,ordinary120<=128/prefix136<=160,2KiBgranule/header-early-wake proof and coldpark/STOP accounting justify the unchangednetwork933 without globalSEV/busy-poll/free scheduler. The chosen real-blocking baremetalA64 implementation and disabled extraneous-event configuration are explicit;this is not a WFE-asNOP portability claim.

Stableheaderstreamingispartofthenotificationrepair:notrewritingheadersafterfirstframe givescontinuationwait12and24notify+8control=32endpoint,soexisting64/edge/192three-edgecontinuationchargesremainvalid. Allbulkstreamsarebank-containedandRootdoesnotinterleaveunrelatedRPCswithinFframes.

Currentrecurrenceisone2639-instructionpublicLEAFwithmacro-onlyinclude,noBL/privatemapframe. Itsoriginal528frameargsoutSP96/headerSP128providecorecontext. FullpublicABI,liveaccumulatorsandallcanonicaloutputsareunchanged;caller20<=2664. Source/include/countproofisinfp2_certificate.md.

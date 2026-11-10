# Random-stride baby/giant inversion with exact trace recurrence

Proposed time_log2=248809/2000=124.4045. This is an analytic claim, not an observed score, runtime pass, qualification or promotion. It uses the full protected p381/r255, actual nonidentity supplied P, uniform full GT target z including identity, threshold0.39 and fixed reference C20000. There are still six production tracks and one computation metric. Complete cold work is below2^124, peak storage below2^138 bytes, advice0. All tables and powers are constructed for this attack. No auxiliary powers are given and no work is amortized over future targets. Explicit finite randomness premise R is disclosed in randomness.md.

Validate canonical inputs, curves and subgroups. Set h=e(P,Q0) for the public G2 generator Q0. It generates GT for every valid nonidentity P. If z=1, return and verify G2 identity. Otherwise the target has unique z=h^a; a is uniform Fr when z is uniform GT. The probability proof below is over uniform a and fresh independent algorithm coins, as required, rather than a worst-case-target assertion.

Sample one uniform b in Fr and twelve independent Fm61 affine coefficients and seventy-two independent raw64-bit NH keys by the capped finite sampler. Compute d=h^(-b) with actual GT exponentiation. Let
m=35823792889802948032681129649173954560, n=51657997126353945492013738969257738240, L=m+n=87481790016156893524694868618431692800.
m andn are multiplesof3*2^96,hence2048 and512. Baby records correspond to h^j for j=1,...,m; giant records to z*d^i=z*h^(-bi) for i=0,...,n-1. Stored scalar is j for a baby, bi modulo r for a giant. Both workers use the SAME independently sampled fingerprint coefficients.

## Exact recurrence and emitted records

Let Phi(X)=a0+b1*t=trace_(Fp12/Fp4)(X)/3, t=w^3,t^2=1+u,u^2=-1. This is four native canonical Fp coordinates. For any GT increment g, its three p^4-conjugates satisfy
Y^3-cY^2+c^(p^2)Y-1, c=3Phi(g).
Consequently v_n=Phi(Z*g^n) satisfies
v_n=c*v_(n-1)-c^(p^2)*v_(n-2)+v_(n-3).
The recurrence follows separately for each conjugate root and therefore also for their Z-weighted sum. It remains valid for g=1 with repeated roots.

Baby increment is h, with cold initial Phi(1),Phi(h),Phi(h^2). Giant increment is d, with cold initial Phi(z),Phi(zd),Phi(zd^2). These are actual computed GT values. Write c=A+B*t, previous values X=a+b*t and Y=d0+e*t, with all letters except t in Fp2. The first two recurrence terms factor EXACTLY as
cX-c^(p^2)Y = A(a-d0)+B*xi*(b+e) + t*(A(b-e)+B(a+d0)).
Four raw two's-complement384-bit differences and four ordinary sums form the temporary inputs. Since|a-d|<p<2^381,its low382 bits are a unique signed382-bit encoding;the top two stored bits are redundant sign extension. A top digitd atstart347 is interpreted as signed35-bit: d ifd<2^34,else d-2^35. All top entries are actual paid canonical signedDigit*C*2^347/Rm values,constructed by TWO complete2^34 halves. Thus every difference product is exactly C*(a-d)/Rm without any hot bias/correction. Sum tables retain all382 bits. This changes no source field/order parameter and drops no independent input information. A two-block vector lookup then computes the two output Fp2 components. The integer accumulators are seeded with the appropriate third-predecessor coordinates, so the final canonical output directly REPLACES that third raw history value. No later copy or key canonicalization writes into the raw history.

recurrence_fp2_step.S is the complete2485-instruction FUSED LEAF throughRET. It embeds the factored table kernel directly after the four raw transforms,with one shared560-byte frame. The standalone recurrence_fp2_map.S remains reference evidence,not an executed hot call or separately charged operation. The fused table body reuses live x0 outputpointer for both seed loads and outputstores;four31-instruction reductions still use actual constructed ordinaryq*p products0..49. Caller16/four pointer rotations gives2505,retained2526. All arithmetic,actual characteristicc=3Phi(increment),raw history and key separation remain exact. Explicit whole-source count,frame/livepointer proof is in fp2_certificate.md.

Use the exact12-Frobenius-fiber key S192 from key_contract.md and the actual canonical_twelve_bfirst.S. All emitted conceptual elements are actual GT by the recurrence proof. Abort if a giant identity is encountered; baby identities cannot occur since1<=j<=m<r. The constant number of cold initial records is charged, and also conservatively assigned the normal hot allowance. Identity targets are solved directly.

Compute full canonical S192 for each point into a disjoint worker scratch,hash it,and retain ONLY its full canonical32-byte Fr scalar as the payload. The full key need not be stored because every proposed exponent is verified. Each hash row keeps its ENTIRE canonical Fq residue,not low54/54/52 bits. The three61-bit residues are packed injectively into183 bits,encoded in the first24 bytes of U32 with9 high zero padding bits. No entropy is attributed to that padding. Color is at24;pending occupancy25 starts zero;26..31 are zero;padding/occupancy are excluded from equality. Metadata is U32 plus ALL128 ordinal bits32..47,exactly48 bytes. Metadata40_certificate.md gives the current exact byte/frame contracts despite its historical filename.

Generate all L records with the same finite streams and point/scalar/color induction. Sort all192 encoded fingerprint bits by TWO stable96-bit two-edge radix passes in radix80_certificate.md. Both full2^96 bucket sets are initialized,including unreachable padded buckets. Final group scan enumerates unordered equal-fingerprint pairs and aborts BEFORE processing more than cap=floor(L/2^29). Count every pair against the cap,including skipped baby/baby and equal-scalar giant pairs.

Fetch each pair's full48-byte metadata and32-byte scalars. A mixed pair proposes all12 exponents a=bi+lambda*j,lambda=p^ell modr,ell0..11. A giant pair with unequal scalars proposes11 nonidentity phases a=(b_t-lambda*b_s)/(1-lambda). Equal-scalar giants and baby/baby pairs are skipped. PRECOMPUTE the public phase scalars and denominator inverses cold,charging construction. For EACH proposed exponent,compute h^a by ordinary binary exponentiation in actual GT and compare with z. If all phases fail,CONTINUE with the next fingerprint pair;do not return a wrong point or stop at that false match. A true useful orbit pair always has one correct phase. On the first verified exponent construct Q=[a]Q0 and verify canonical encoding,curve/subgroup membership and the full e(P,Q)=z before returning. The finalG2 construction/pairing is done at most once and remains cold paid. There can be up to12*cap exponent checks,NOT merely12 total;all are charged in the per-confirmation bound below. No production discrete-log or GT equality oracle is assumed.

## Exact useful-collision probability

Write u0=0xd201000000010000 and lambda0=u0^2. The full protected parameters satisfy r=u0^4-u0^2+1, p=-u0 modr. At the displayed m,n,
m<u0^2-1, 2m^2<r, 3n<lambda0, 4n(lambda0+1)<r.
These are exact integer inequalities, not toy parameters.

The union of baby Frobenius exponents
Sexp={p^ell*j modr:0<=ell<12,1<=j<=m}
has size s=12*(m-floor(m/u0)). Here is a complete short-interval proof. The twelve roots are ±1,±u0,±u0^2,±u0^3,±(u0^2-1),±(u0^3-u0). The ±u0^2 and ±(u0^2-1) products lie outside[1,m] without modular wrap because m<u0^2-1 and m*u0^2<r-m. Of ±u0 only +u0 produces pairs j'=u0*j, exactly floor(m/u0). The inverse phases ±(u0^3-u0)=∓u0^(-1) give only the same reverse pairs: multiply by u0 and use u0*m+m<r. Neither -1 nor -u0 produces another small positive value.

For ±u0^3 put x=-u0. The congruence j'=x^3*j modr has integer lattice basis (x,x^2-1),(x^2-1,-x), determinant-r. Both basis vectors have squared norm r and are orthogonal, so every nonzero lattice vector has norm at least sqrt(r). The basis spans the full congruence lattice because its determinant is its index r. A pair with |j|,|j'|<=m has squared norm<=2m^2<r, impossible. Negating j' covers the other sign. Finally u0^2>m prevents any three-point +u0 chain. Thus each duplicate baby orbit contains exactly two baby exponents and there are precisely floor(m/u0) such pairs. All nonidentity orbits have size12, proving s.

For every i=0,...,n-1 and c in Sexp, select a mixed event a=i*b+c. There are M=n*s distinct events. Each has probability1/r over uniform independent a,b. Two different slopes i give a full-rank2x2 system and probability1/r^2; same slope with different nonzero c is disjoint. Each event gives an actual useful mixed collision.

Select FOUR giant/giant phase families. Put lambda=lambda0, so lambda^2-lambda+1=0 modulo r. For each0<=i<k<n use both sigma=lambda and lambda^(-1)=1-lambda. Also use sigma=-lambda and -lambda^(-1)=lambda-1 ONLY when i and k have different residues modulo3. In all cases define
c_sigma(i,k)=(i-sigma*k)/(1-sigma),
and the guarded event b!=0, a=b*c_sigma(i,k).

Here is a complete injectivity/disjointness certificate in ordinary integer coordinates:
sigma=lambda: c=k+lambda*(i-k);
sigma=lambda^(-1): c=i+lambda*(k-i);
sigma=-lambda: 3c=(2i+k)+lambda*(k-i);
sigma=-lambda^(-1): 3c=(i+2k)+lambda*(i-k).
These identities follow by lambda^2=lambda-1 and are equations in Fr. Their coefficient pairs are NOT treated as arbitrary-width machine words; they are only the finite success proof.

For the first two families, the lambda coefficient has opposite nonzero signs; constants lie in[0,n). Two such representations can be equal modulo r only when their integer coefficients agree: |Delta_constant|<n<lambda and the whole integer difference has absolute value<2n(lambda+1)<r. Thus each is injective and they are disjoint. For either last family, the scaled constant lies in[0,3n), the lambda coefficient has magnitude<n, and the pair determines i,k uniquely (the integer linear transform has determinant±3). Equality within or across the last families similarly forces equality of their integer coordinates, using3n<lambda and4n(lambda+1)<r. Their lambda signs are opposite, so the last families are disjoint.

To compare a last-family slope with a first-family slope, multiply the latter by3. Constant-coordinate differences are less than3n; lambda-coordinate differences less than4n. The whole integer difference has absolute value<4n(lambda+1)<r, and3n<lambda again forces both coordinates to agree. But that would make the last family's lambda coefficient divisible by3, contrary to the selection i!=k mod3. Hence these families are disjoint from the first two. Each slope is nonzero: its lambda coefficient is nonzero, dominates the bounded constant if it has opposite sign, and the full absolute value is below r.

Since n is divisible by3, the number of pairs i<k in different residue classes is exactly n^2/3. Therefore the number of DISTINCT nonzero selected giant slopes is
G0=2choose(n,2)+2*n^2/3=5*n^2/3-n.
The four phases are p^2,p^10,p^8,p^4; all are tested by the algorithm, with inverse orientation if needed. The actual algorithm still tries all useful phases, so omitting duplicate slopes only affects this lower-bound proof.

## Additional order-four giant events, without extra algorithmic work

Put sigma=u0^3 modulo r. Then sigma^2=-1 and sigma=p^9; its negative is p^3. Both phases are already among the eleven giant phases actually tested. For every integer gap j with1<=j<=J=floor((n-1)/2) and center j<=t<n-j, the existing records i=t-j,k=t+j yield the two slopes c=t+sigma*j and c=t-sigma*j. Indeed (1+sigma)/(1-sigma)=sigma. The endpoints are nonidentity whenever b!=0 because c-i=j*(1+/-sigma)!=0;the stored scalars differ because2j*b!=0.

Select gaps by residues only for this LOWER-BOUND proof;the actual algorithm still enumerates every matching pair and all existing phases. Set W=ceil(n/u0)+1 and choose j mod u0 in the three inclusive intervals
[W,floor((u0-W)/3)],
[ceil((u0+W)/3),floor((2u0-W)/3)],
[ceil((2u0+W)/3),u0-W].
The current full-size integers satisfy4W<=u0 and each interval is nonempty. Every selected j and3j has distance at leastW from EVERY multiple ofu0. This follows directly by reducing3*(j mod u0) in respectively[W,u0-W],[u0+W,2u0-W],[2u0+W,3u0-W]. The first and third actual ranges can be narrower,which only strengthens the assertion.

All new slopes are distinct,including across signs. Equality implies A+sigma*B=0 modr with |A|<n,|B|<n. The already proved orthogonal index-r lattice for x=-u0 applies to(B,A),because A=x^3*B. Every nonzero lattice vector has squared norm at leastr,whereas A^2+B^2<2n^2<r. Hence A=B=0. Equal signs force same t,j;opposite signs would give j+j'=0,impossible. The same argument with A=t,B=+/-j shows every new slope is nonzero.

To prove disjointness from ALL four earlier families,write an earlier slope as (Aold+lambda0*Bold)/D,where D is1 or3,0<=Aold<D*n and|Bold|<n,using the displayed earlier integer coordinates. Equality with a new slope gives
H+epsilon*D*u0^3*j-lambda0*Bold=ell*r,
where H=D*t-Aold,|H|<3n,epsilon=+1 or-1. The exact inequality3n+lambda0*n+3*u0^3*J<u0*r implies|ell|<u0. Modulo lambda0=u0^2,r is1,so H=ell modulo lambda0. Since3n+u0<lambda0,this forces the INTEGER equality H=ell. Dividing the remaining equality byu0^2 gives
Bold=u0*(epsilon*D*j-ell*u0)+ell.
The selected residue intervals make |epsilon*D*j-ell*u0|>=W for BOTH D=1,3 and every integerell. Consequently |Bold|>=u0*W-|ell|>u0*W-u0>=n,a contradiction. No heuristic collision independence or lattice search is used.

Count the selected gaps and their sum exactly. With Q,T=divmod(J,u0),for each interval[lo,hi] let width=hi-lo+1 and last=max(0,min(hi,T)-lo+1). Its contributions are
Kpart=Q*width+last,
Spart=width*u0*Q*(Q-1)/2+Q*(lo+hi)*width/2+last*Q*u0+(lo+min(hi,T))*last/2.
These are ordinary arithmetic-series identities for complete residue blocks and the last partial block. Sum over the three intervals to obtainK andS. There are n-2j centers per selected gap and two signs,so the number of ADDITIONAL distinct slopes in this EVEN-GAP subset is Eeven=2*(n*K-2*S). This first subset contributes Eeven;the E constant below is the FINAL total after the odd-gap extension. Current exact full-size values are:
W=3413739941783706337.
K=6406297912961268207630192833709814128.
S=82734129794580117880121541879326061613967083752902523036371242204588707984.
G0=4447581111843987483824230074440296227393054988392986543235395224564672757760.
E=661873038356640942637968564240750280158492338842100209621688294436816775630.
All stated interval,lattice,separation and count inequalities are independently reconstructed by the shared exact ledger. Counting already tested events adds no phase,table,record or hot operation.

## Odd index gaps add another disjoint subset

The preceding selected pairs have even index difference. Now also take odd Jgap=k-i in[1,n-1],and write T=i+k=2i+Jgap. The same two already-tested order-four phases give c=(T+epsilon*sigma*Jgap)/2,epsilon=+/-1,where0<=i<n-Jgap. Since r is odd,division by2 is valid. Their endpoints are nonidentity under b!=0 because c-i=Jgap*(1+epsilon*sigma)/2!=0. There is no extra algorithm phase orruntimefilter.

All order-four slopes,even and odd gap alike,are distinct. Put the earlier even gaps in the SAME representation T=2t,Jgap=2j. An equality gives A+sigma*B=0 with|A|<2n,|B|<2n,including opposite signs. The same index-r lattice applies and8n^2<r,so A=B=0. Equal signs force identical T,Jgap;opposite signs would require positive gaps to sum tozero. An odd gap cannot equal an even gap. Nonzero slopes follow from T^2+Jgap^2<5n^2<r.

To exclude the four earlier lambda-family slopes,write them as(Aold+lambda0*Bold)/D withD1 or3. Equality gives
H+epsilon*D*u0^3*Jgap-2lambda0*Bold=ell*r,
where H=D*T-2Aold and|H|<6n. The exact inequality6n+2lambda0*n+3u0^3*(n-1)<u0*r implies|ell|<u0. Reduction modulo lambda0 yields H=ell+v*lambda0. Since6n+u0<2lambda0,the integer v is one of-1,0,1. Consequently
2Bold=u0*(epsilon*D*Jgap-ell*u0)+ell+v.

Crucially,the overlap quotient ell must be ODD. Both D and T are odd,so H=D*T-2Aold is odd. Every other term on the left of the congruence is even becauseu0 is even. Since r is odd,the integer equality forcesell odd. The carryv does not change this fact. Therefore epsilon*D*Jgap-ell*u0 measures distance ONLY from odd multiples ofu0,not from every multiple. This parity conclusion applies to all four old families becauseD is1 or3.

Set Wodd=ceil(2n/u0)+2. Work modulo2u0,which is even,and select ONLY odd residues in the four inclusive intervals
[0,floor((u0-Wodd)/3)],
[ceil((u0+Wodd)/3),u0-Wodd],
[u0+Wodd,floor((5u0-Wodd)/3)],
[ceil((5u0+Wodd)/3),2u0-1].
Round each endpoint inward to odd integers. All are nonempty at the current integers and2Wodd<=u0. Every selected residue has distance at leastWodd from odd multiples ofu0. Three times the first,second,third,fourth intervals lie respectively inside[0,u0-Wodd],[u0+Wodd,3u0-Wodd],[3u0+Wodd,5u0-Wodd],[5u0+Wodd,6u0]. Thus3Jgap has the SAME distance from every odd multiple. Reduction modulo2u0 preserves that distance forJgap and3Jgap. Negativeepsilon preserves the odd-multiple set as well.

For BOTH D1 and3 and the ODD quotientell,therefore|epsilon*D*Jgap-ell*u0|>=Wodd. The exact overlap equation implies
|2Bold|>=u0*Wodd-|ell+v|>u0*Wodd-u0-1>=2n,
contradicting|Bold|<n. This excludes EVERY overlap with the four oldfamilies. The earlier8n^2<r lattice argument excludes overlaps within the neworder-four families and across gap parity. No sampled estimate,independence assumption or missedcarry is used.

Since2u0 is even,gap parity equals residue parity. LetQo,To=divmod(n-1,2u0). For EACH rounded odd interval[lo,hi],let width=(hi-lo)/2+1,last=max(0,floor((min(hi,To)-lo)/2)+1). Its contributions are
Kopart=Qo*width+last,
Sopart=width*(2u0)*Qo*(Qo-1)/2+Qo*(lo+hi)*width/2+last*Qo*(2u0)+last*(lo+last-1).
Sum the FOUR contributions to obtainKo,So. There are n-Jgap centers per gap and two signs,so Eodd=2*(n*Ko-So). The previous even-gap count remainsEeven=2*(n*K-2*S). Set E=Eeven+Eodd,andG=G0+E in every subsequent event/probability formula. Current full-size exact constants:
Eeven=330936519178320471451244508197461930363947900415534773332757699573800877504.
Eodd=330936519178320471186724056043288349794544438426565436288930594863015898126.
Wodd=6827479883567412674.
Ko=6406297912961268201360159912221412491.
So=165468259589160235568605967195346428218817649759034808356779265693236406777.
The rounded odd intervals are [[1, 2768298779791410025], [7319952035503018477, 8304896339374230077], [21959856106509055427, 22944800410380267027], [27496453666091875479, 30264752445883285503]]. The shared ledger reconstructs their endpoints and verifies every full-size range,distance,count andseparation inequality.

Each giant event has probability(r-1)/r^2 because b!=0. Different giant slopes are disjoint under that guard. Their points are nonidentity, since c_sigma-i=sigma*(i-k)/(1-sigma)!=0. Scalars bi,bk are distinct. A mixed/giant event pair is either an independent linear system with probability at most1/r^2 or parallel and disjoint, because the mixed intercept c is nonzero. These are exact two-variable finite-field facts, not collision independence or a hash heuristic.

Let
mu=M/r+G*(r-1)/r^2,
U2=(choose(M,2)+M*G)/r^2.
The omitted giant/giant intersections are exactly zero. The second Bonferroni inequality gives useful-event probability at least mu-U2. Here
s=429885514677635376363765227293344051480,
M=22207024681878475128677467500612484838013543547553157585445404583662924595200,
G=5109454150200628426462198638681046507551547327235086752857083519001489533390.
All selected events correspond to actual records; inverse phase orientation is harmless because the algorithm tests every required Frobenius phase.

For confirmations,distinct giant exponent pairs are uniform Fr^2 over uniform a,b,and a mixed giant exponent is uniform Fr. True orbit matches therefore have expectation at most12/r per such pair. Deterministic baby duplicates contribute exactly floor(m/u0). For ANY distinct S192 strings the independently sampled full-residue fingerprint has collision bound delta=(2^-63+1/q)^3,q=2^61-1,proved self-contained in nh_certificate.md. Conditional on complete point sequences,the hash keys are independent;no independence between different record pairs is assumed.

The expected number of ALL equal-fingerprint pairs is bounded by choose(L,2)*(12/r+delta)+floor(m/u0). Thus Markov's inequality bounds cap-abort probability by
caploss=(choose(L,2)*(12/r+delta)+floor(m/u0))/(floor(L/2^29)+1).
This pays false hash collisions,ignored baby pairs and all true matches;there is no free full-key comparison. The giant-identity union bound is n/r,and capped sampling/extraction loss below2^-970. Identity targets are solved directly;unused coins may be extended hypothetically for the unconditional probability argument.

The same exact group-event mu,U2 above give
mu-U2-caploss-n/r-2^-970 >3900031/10000000>0.39.
Diagnostics:mu=0.5209501842173337,U2=0.13094707999339925,caploss=3.7556202199638514e-09,final=0.3900031004683142. These decimals are not the certificate;the displayed integer/rational formulas define it. The broader fingerprint and tighter cap replace the old160-bit/full-key-confirmation calculation entirely.

## Complete work and storage

Native fused recurrence2485/key710/NH341 retain2526/721/350. Non-scalar generation is(2526+721+350+712)L=4309L. Baby7/giant21 native helpers retain12/28 withcallers. At mostm babies andn giants inEVERY execution give exact scalar bound12m+28n,includingfailures andlastunused increments;NO20L balancedaverage isused. Main712 is unchanged perrecord. Local two-pass720L,scan128L,confirmation1L,inner spare256L,network920L andouter spare80L give total non-scalar coefficient6414L. Hence N=6414L+12m+28n+2^124. Differentpublic128-bit streamlimits require no wider counter ornewhotoperation;both remain multiples512/2048. Allfull buckets/cold initialization remainpaid.

A complete processed fingerprint pair costs less than2^29 in the normalized scalar-work numerator C*H+R. Binary exponentiation of a255-bit scalar uses at most512 dense reference Fp12 multiplications,including squarings implemented as multiplication. At most12 phases give H<=6144. Public phase inverses are constructed cold. Per-phase Fr products can use256 bit-serial steps on five64-bit limbs,at most128 scalar instructions per step:shift/add an operand belowr,then at most two trial subtractions since the intermediate is below3r. Three such products,adds and phase bookkeeping are below2^20 per phase. Additional exponent-loop/call/buffer/control work is at most256 scalar instructions per dense operation;metadata/payload RPCs,comparisons,pair enumeration and all fixed work fit the generous remainder of R<=2^25. Thus C*6144+2^25=156434432<2^28<2^29 at protectedC20000. No production reference price is changed. Correct terminalG2 output is at most once,cold paid. Since cap=floor(L/2^29),ALL processed confirmations together cost at most L in this numerator,exactly the retained1L above. Detection/termination beyond the cap has no exponent loop and is cold paid. There are no uncharged failures or twelve-global-check assumptions.

The complete two96-bit sort construction is below2^121,including BOTH full2^96 bucket sets,all queues,prefixes,zeroing,flags,padding restoration,initial/tail exceptions and bounded finish. All processor code/page/bootstrap remains below2^121,edge/controller initialization below2^118,and other field/table/coin/phase/output work within its fixed allowance. Their sum remains below2^124. See current radix/architecture/metadata certificates for explicit arrays and transfers. Generic data capacity2^136,private2^119,edges2^111 and algorithm-specific data<2^110 give peak<2^138 bytes,advice0. Unused capacity supplies no answers. All actual pointers are below2^47,IDs96-bit and counters128-bit ordinary words.

For full L=m+n,N=6414L+12m+28n+2^124 is the upper normalized-work numerator,including ALLweighted scalar work and C*H confirmation work. EXACT N^2000<20000^2000*2^248809 proves proposed248809/2000=124.4045. Public integerN=563005778245778419571527906478514667257856. This is an analytical bound,not an observedproductionrun orphysicalfeasibility claim. All preprocessing,coins,arithmetic,sort/lookup,transfers,trials/failures,verification/output are charged withoutamortization. Streamlengths areattackwork parameters;protected field/order/inputdistribution/threshold/referenceprices areunchanged.

The low-limb sign normal form and b1-first representative order in key_contract.md preserve exact twelve-Frobenius fibers. The hash is proved universal for arbitrary distinct192-byte S strings. Discarding stored S changes neither true useful events nor correctness,because every candidate exponent is checked in GT;only the charged false-fingerprint workload and cap calculation change. The scoped exponent-equality oracle in audit.py is TEST-ONLY algebra/control evidence,not an oracle available to the production attack. Shared llm-verifier alone executes candidates,reviews and finalizes;all full parameters,0.39,C20000,six tracks and sole metric remain unchanged.


Current register refinement:the extracted digit inx3 is dead immediately after MADD computes its entry address. Reusex3 for the first loaded table limb,in BOTH vectors,and for the first loaded modulus limb during reduction. Every subsequent digit extraction overwritesx3 before use. This leavescalleex20 entirely untouched. Ten usedcallees19,21..29 restore in5 LDP rather than6 instructions;the prefix still needs7 stores for10callees and3livearguments. Frame560 remains:callee0..79,out80/table88/const96,workspace112..167,temp176..559. LRx30 remains untouched. Actual completeleaf2485,prefix9/transforms96/inputADD6/outLDR1/body2366/restore7. Caller20 gives2505,retained2526;all arithmetic/tablegeometry/probability remainunchanged. No hidden free registerrestore isclaimed.

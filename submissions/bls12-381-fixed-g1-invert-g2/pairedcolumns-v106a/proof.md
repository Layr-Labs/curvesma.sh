# Random-stride baby/giant inversion with exact trace recurrence

Proposed time_log2=621303/5000=124.2606. This is an analytic claim, not an observed score, runtime pass, qualification or promotion. It uses the full protected p381/r255, actual nonidentity supplied P, uniform full GT target z including identity, threshold0.39 and fixed reference C20000. There are still six production tracks and one computation metric. Complete cold work is below2^124, peak storage below2^138 bytes, advice0. All tables and powers are constructed for this attack. No auxiliary powers are given and no work is amortized over future targets. Explicit finite randomness premise R is disclosed in randomness.md.

Validate canonical inputs, curves and subgroups. Set h=e(P,Q0) for the public G2 generator Q0. It generates GT for every valid nonidentity P. If z=1, return and verify G2 identity. Otherwise the target has unique z=h^a; a is uniform Fr when z is uniform GT. The probability proof below is over uniform a and fresh independent algorithm coins, as required, rather than a worst-case-target assertion.

Sample one uniform b in Fr and nine independent Fm61 affine coefficients and seventy-two independent raw64-bit NH keys by the capped finite sampler. Compute d=h^(-b) with actual GT exponentiation. Let
m=28721772558436203898626428180203831296, n=55754029111986100220015185618117066752, L=m+n=84475801670422304118641613798320898048.
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

Use the exact12-Frobenius-fiber key S192 from key_contract.md and the actual canonical_twelve_bfirst.S. All emitted conceptual elements are actual GT by the recurrence proof. Abort if a giant identity is encountered; baby identities cannot occur since1<=j<=m<r. The constant number of cold initial records is charged, and also conservatively assigned the normal hot allowance. Identity targets are solved directly.

SelectthefullcanonicalS192asrawAreal/selectedAimag/signedBregionsandhashALL192bytesdirectlyinthecombinedleaf;physicalS192copyisnotneeded. RetainONLYitsfullcanonical32-byteFrscalaraspayload. The full key need not be stored because every proposed exponent is verified. EachhashrowkeepsLOW61ofrawsumofTWOcanonicalFqtableentriesplusDIRECT61DWITHOUTmodqnormalization;theNEWprojection/5qproofdominatedbyOLDdeltabudgetisinnh_certificate.md. The three61-bit residues are packed injectively into183 bits. The actual combined leaf returns all183 hash bits in registersx1/x0/x5,statusx7 (zero=identity,nonnull=success);the reference standalone NH outputsU32. The source packs color in bit183 (byte22bit7),pending0..63 in byte23,and the ENTIRE128-bit ordinal at24..39. Thus metadata is40 bytes with318 meaningful bits and2 reserved zero bits. Before EVERY upper96 digit/equality operation,mask word2 to its low55 bits,so the logical192-bit fingerprint encoding still has9 high zero bits and excludes color/occupancy. Worker hash/color source pack adds3 actual instructions using the explicitly checked full register return;the existing ordinalSTP completesall40bytes;all masks and both ordinal limbs are paid. Payload remains full32-byte canonical Fr scalar. Pending23 is restored before every emitted row. metadata40_certificate.md andradix80_certificate.md specify complete five-frame source/20-frame bulk/copy/address contracts.

Generate all L records with the same finite streams and point/scalar/color induction. Sort all192 encoded fingerprint bits by TWO stable96-bit two-edge radix passes in radix80_certificate.md. Both full2^96 bucket sets are initialized,including unreachable padded buckets. Final group scan enumerates unordered equal-fingerprint pairs and aborts BEFORE processing more than cap=floor(L/2^29). Count every pair against the cap,including skipped baby/baby and equal-scalar giant pairs.

Fetch each pair's full40-byte metadata and32-byte scalars. A mixed pair proposes all12 exponents a=bi+lambda*j,lambda=p^ell modr,ell0..11. A giant pair with unequal scalars proposes11 nonidentity phases a=(b_t-lambda*b_s)/(1-lambda). Equal-scalar giants and baby/baby pairs are skipped. PRECOMPUTE the public phase scalars and denominator inverses cold,charging construction. For EACH proposed exponent,compute h^a by ordinary binary exponentiation in actual GT and compare with z. If all phases fail,CONTINUE with the next fingerprint pair;do not return a wrong point or stop at that false match. A true useful orbit pair always has one correct phase. On the first verified exponent construct Q=[a]Q0 and verify canonical encoding,curve/subgroup membership and the full e(P,Q)=z before returning. The finalG2 construction/pairing is done at most once and remains cold paid. There can be up to12*cap exponent checks,NOT merely12 total;all are charged in the per-confirmation bound below. No production discrete-log or GT equality oracle is assumed.

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

Select gaps ONLY for the lower-bound proof;the algorithm still enumerates every matching pair and all existing phases. Put W=ceil(n/u0)+1 and choose j mod u0 in the SINGLE inclusive interval[W,u0-W]. It is nonempty at the current full-size integers. Every selectedj has distance at leastW from every multiple ofu0. Unlike an earlier three-interval proof,no restriction on3j is needed:the D3 oldfamilies are excluded by their already-required nonzero residue modulo3,as proved below.

All new slopes are distinct,including across signs. Equality implies A+sigma*B=0 modr with |A|<n,|B|<n. The already proved orthogonal index-r lattice for x=-u0 applies to(B,A),because A=x^3*B. Every nonzero lattice vector has squared norm at leastr,whereas A^2+B^2<2n^2<r. Hence A=B=0. Equal signs force same t,j;opposite signs would give j+j'=0,impossible. The same argument with A=t,B=+/-j shows every new slope is nonzero.

To prove disjointness from ALL four earlier families,write an earlier slope as (Aold+lambda0*Bold)/D,where D is1 or3,0<=Aold<D*n and|Bold|<n,using the displayed earlier integer coordinates. Equality with a new slope gives
H+epsilon*D*u0^3*j-lambda0*Bold=ell*r,
where H=D*t-Aold,|H|<3n,epsilon=+1 or-1. The exact inequality3n+lambda0*n+3*u0^3*J<u0*r implies|ell|<u0. Modulo lambda0=u0^2,r is1,so H=ell modulo lambda0. Since3n+u0<lambda0,this forces the INTEGER equality H=ell. Dividing the remaining equality byu0^2 gives
Bold=u0*(epsilon*D*j-ell*u0)+ell.
ForD=1,the selected residue interval gives|epsilon*j-ell*u0|>=W. Hence|Bold|>=u0*W-|ell|>u0*W-u0>=n,contradicting|Bold|<n.

ForD=3,no distance condition on3j is required at all. The exact integer formula becomes Bold=epsilon*3*u0*j-ell*u0^2+ell. Sinceu0^2=1 modulo3,this forcesBold=0 modulo3. But BOTH earlierD3 families are selected ONLY when i!=k modulo3,andtheir displayed lambda coefficientBold is+/- (k-i),which is NONZERO modulo3. This contradiction excludes everyD3 overlap for ALL even gaps. This is a short integer residue argument,not an independence assumption. Its applicability usesH=ell exactly,already established by3n+u0<lambda0;there is no discarded integercarry in the EVEN-GAP comparison. The separate ODD-GAP comparison below keeps its possiblev=-1,0,1.

Count the selected gaps and their sum exactly. With Q,T=divmod(J,u0),for each interval[lo,hi] let width=hi-lo+1 and last=max(0,min(hi,T)-lo+1). Its contributions are
Kpart=Q*width+last,
Spart=width*u0*Q*(Q-1)/2+Q*(lo+hi)*width/2+last*Q*u0+(lo+min(hi,T))*last/2.
These are ordinary arithmetic-series identities for complete residue blocks and the last partial block. Use the SINGLE interval[W,u0-W] to obtainK andS. There are n-2j centers per selected gap and two signs,so the number of ADDITIONAL distinct slopes in this EVEN-GAP subset is Eevenbase=2*(n*K-2*S). This first subset contributes Eevenbase;the E constant below is the FINAL total after the odd-gap extension. Current exact full-size values are:
W=3684419967530245091.
K=14302064058857477826625796405656813195.
S=199349423974757476177660504161131905472934684938552629813335452915717734400.
G0=5180852937033655951780254535455081558071927331832682846889928897875922649088.
E=2324736173859401461449319720000534276329068601020947056782252115728717963516.
All stated interval,lattice,separation and count inequalities are independently reconstructed by the shared exact ledger. Counting already tested events adds no phase,table,record or hot operation.

## Additional directed centers for omitted even residues

The central residue interval leaves j mod u0 near0 oru0. The D3 contradiction just proved applies to ALL even gaps and still excludes those two old families. The same H=ell integer equality applies. It remains to exclude D1 using actual old endpoint bounds. Write j=u0*q+rho. In the low band rho=d with0<d<W, the ONLY possible overlap quotient is ell=epsilon*q. In the high band rho=u0-d with0<d<W, it is ell=epsilon*(q+1). Indeed an overlap requires |epsilon*j-ell*u0|<(n+u0)/u0<=W; any other integer quotient is at distance at leastu0-d>u0-W>=W. The printed full-size inequalities imply q+1<u0,so all signs below are exact. Set Dactual=u0*d.

For the low band, Bold=epsilon*(Dactual+q) and Aold=t-epsilon*q. For epsilon+ the old positive-coefficient family requires Aold+Bold=t+Dactual<n, so t>=n-Dactual excludes it. For epsilon- the old negative-coefficient family requires Aold>=-Bold, equivalent to t>=Dactual, so t<Dactual excludes it. For the high band, Bold=epsilon*(-Dactual+q+1),whose sign is -epsilon since Dactual>=u0>q+1. Here epsilon+ is safely excluded by t<Dactual,and epsilon- by t>=n-Dactual. These are necessary endpoint inequalities for the actual two D1 index families; no independent-coordinate or generic-distance approximation is substituted. Both signs of each newly selected j therefore have one safe low-center or safe high-center direction.

Use256 disjoint distance strata indexed l=0,...,255 with dmin=floor(l*W/256),dnext=floor((l+1)*W/256). Retain only d in[dmin,dnext-1],and let D=u0*dmin. Skip the entire l0 stratum with D0;it adds no slopes. The two inclusive residue intervals are[dmin,dnext-1] and[u0-dnext+1,u0-dmin]. Every retained interval is disjoint from[W,u0-W] and from all other retained intervals. For both bands choose the sign-specific safe centers t<D or t>=n-D. Since D<=Dactual,these are subsets of the rigorously safe endpoint regions. At the current integers 0<D<n. The two signs have equal counts by t mapping to n-1-t on the integer interval[j,n-j). For each gap, the combined number of new centers is exactly 2*max(0,min(n-j,D)-j).

For a given stratum define K(z),S(z) as the count and sum of 1<=j<=z in its two residue intervals modulo u0. Compute each by the same complete-block/last-partial arithmetic series already proved above. If2D<=n, the contribution is2*(D*K(min(J,D))-S(min(J,D))). If2D>n, put z=min(J,n-D);the contribution is2*(n*K(J)-2*S(J)+(D-n)*K(z)+S(z)). This is the exact split of min(n-j,D) at j=n-D. Values at a zero-weight endpoint do not affect it. Sum all255 nonzero strata to obtain Eevencenter. Define Eevenbase=2*(n*K-2*S) for the original central residue interval, and Eeven=Eevenbase+Eevencenter.

All newly selected slopes are nonzero and pairwise distinct by the same index-r lattice bound8n^2<r. The disjoint residue intervals exclude the previous even subset; odd and even gaps are disjoint by the same numerator-lattice argument. Their events are therefore additional guarded giant slopes, with exactly the same (r-1)/r^2 individual probability and zero giant/giant intersections. The algorithm already tests both order-four phases for all informative fingerprint pairs. There is NO new runtime filter, phase, table, record bit or operation. Only the rigorous lower-bound event selection and publicly chosen full-size stream limits change.

Eevenbase=797397695899029904611042585797699544115855867674440643564458482451947847680.
Eevencenter=376950853962543075637188302962774670166735904930026468981434304438232166672.

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
Sum the FOUR contributions to obtainKo,So. There are n-Jgap centers per gap and two signs,so this unrestricted-center odd subset contributes Eoldodd=2*(n*Ko-So). The even-gap total is Eeven=Eevenbase+Eevencenter as proved above. The Eodd constant below is the FINAL odd total after the additional centered subset;set E=Eeven+Eodd,andG=G0+E in every subsequent event/probability formula. Current full-size exact constants:
Eeven=1174348549861572980248230888760474214282591772604467112545892786890180014352.
Eodd=1150387623997828481201088831240060062046476828416479944236359328838537949164.
Wodd=7368839935060490182.
Ko=5252097060767096300074467774422011458.
So=146412786212492658559301091046873981996243973178581829457922007572806365284.
The rounded odd intervals are [[1, 2587845429293717523], [7500405386000710979, 7763536287881152569], [22501216158002132935, 22764347059882574525], [27676907016589567981, 30264752445883285503]]. The shared ledger reconstructs their endpoints and verifies every full-size range,distance,count andseparation inequality.

## Additional odd gaps with a bounded center

We can also use selected odd gaps omitted by the four preceding intervals,with direction-dependent center bounds that exclude their only possible carry. Keep the SAME Wodd and select odd residues modulo2u0 in the TWO gaps between those four intervals:
[floor((u0-Wodd)/3)+1,ceil((u0+Wodd)/3)-1],
[floor((5u0-Wodd)/3)+1,ceil((5u0+Wodd)/3)-1],
again rounding inward to odd integers. These are disjoint from the previous four odd intervals. Each lies in[0,u0-Wodd] or[u0+Wodd,2u0-1],so it still excludes D1 overlaps by distance from ODD multiples ofu0. The D1 comparison itself has|H|<2n and2n+u0<lambda0,hence its carry is0;the earlier stronger distance inequality also suffices. We do NOT assert distance for3Jgap in these new intervals.

We sharpen the safe center bounds by 256 residue strata. This remains ONLY a success-proof subset of already-tested pairs/phases. The overlap constraints below are exact ordinary integer statements,not a new algorithm oracle or runtimefilter.

First,the oldD3 coordinates satisfyAold=Bold modulo3. FromH=3T-2Aold=ell+v*lambda0 and2Bold=epsilon*3u0*Jgap-ell*u0^2+ell+v,we obtainv=ell modulo3. In the FIRST new interval writeJgap=2u0*q+rho with|3rho-u0|<Wodd;in the SECOND,|3rho-5u0|<Wodd. A hypothetical overlap has|epsilon*3Jgap-ell*u0|<(2n+u0+1)/u0<Wodd. Sinceell is odd andWodd<u0,the unique possible oddquotient isell=epsilon*(6q+1) or epsilon*(6q+5),respectively. Thusellmod3 equals epsilon or-epsilon,andfixesv=+1 or-1. A quotient outside|ell|<u0 simply means nooverlap. AllD1 overlaps remain excluded by distance fromoddmultiples.

Now use the actual oldD3 index bounds as well. IfBold>0,its coordinate satisfiesBold<=Aold<3n-2Bold;ifBold<0,then-2Bold<=Aold<3n+Bold. Both follow by solving the displayed integer coordinate transform for0<=i<k<n. Therefore in BOTH families
|Bold|<=Aold<3n-|Bold|.

Let delta=3rho-c,wherec=u0 or5u0 for the two new intervals. Partition its absolute value into256 disjoint strata:forlevel l=0,...,255 putdmin=floor(l*Wodd/256),dnext=floor((l+1)*Wodd/256),andselectdmin<=|delta|<dnext. The four residue intervals per stratum,before rounding inward to odd endpoints,are forEACHc in{u0,5u0}:
[floor((c-dnext)/3)+1,floor((c-dmin)/3)],
[ceil((c+dmin)/3),ceil((c+dnext)/3)-1].
These partition the two new odd residue intervals. Because neitheru0 nor5u0 is divisible by3,delta=0 never occurs,so the positive/negative pieces atdmin0 do not overlap. Distinct strata are disjoint by their half-open absolute-value bounds. All remain disjoint from the four unrestricted odd intervals.

The unique quotient gives|epsilon*3Jgap-ell*u0|=|delta|>=dmin. The exact coefficient equation and|ell+v|<=u0 imply
2|Bold|>=u0*max(0,dmin-1)=beta.
Refine the directed endpoint bound by the sign ofdelta. Because|ell+v|<=u0 anddelta is a nonzero integer,the coefficient equation shows sign(Bold)=epsilon*sign(delta) in every possible overlap;zeroBold is already excluded by theoldfamilies. In thefirstintervalv=epsilon,and in thesecondv=-epsilon. Forv+1 the relevant lowerendpoint ismax(Bold,-2Bold);forv-1 the relevant upperendpoint distance from3n ismin(2Bold,-Bold). Inspection of these EXACT two sign cases gives a factorf=2 when(c=u0 anddelta<0) OR(c=5u0 anddelta>0);the other two cases havef=1. Bothphases for a given residue use the SAMEfactor,withopposite safe-center directions. Thus the relevant endpoint distanceDend satisfies2Dend>=f*u0*max(0,dmin-1). Setbeta to THIS stronger bound for thedirectedcomparison. No endpoint is silently doubled when its sign doesnot permit it.

For eachstratum split its four rounded intervals into two factor groups:negativefirst/positivesecond usef2,andpositivefirst/negativesecond usef1. Each group has its ownbeta andAlevel/Blevel. The groups partition the stratum and their gap sets are disjoint. For this stratum chooseAlevel as the smallest odd integer at leastmax(1,ceil((6n-lambda0+u0-beta)/3)),andBlevel=2n-2-Alevel. The current full-size integers have1<=Alevel<n,3Alevel>=6n-lambda0+u0-beta,and3Blevel<=lambda0-u0+beta. Ifv=+1,selectONLYT<=Blevel. ThenH=3T-2Aold<=3T-beta<=lambda0-u0,soH-ell<lambda0,contradictingv=+1. Ifv=-1,selectONLYT>=Alevel. ThenH>3T-6n+beta>=-lambda0+u0,soH-ell>-lambda0,contradictingv=-1. This excludes the ONLY possible carry for each phase. Neither sign is assumed independent. Larger guaranteed|Bold| therefore safely permits a larger center range.

For each stratum AND factor group compute gap count/sumKt/St to n-1,andKl/Sl toAlevel,using the same modulo2u0 arithmetic-series formula on its two rounded odd intervals. ReflectionT->2n-2-T makes the two directed safe-center counts equal. ForJgap<=Alevel eachsign has n-(Alevel+Jgap)/2 centers;otherwise n-Jgap. Hence the stratum contributes
Elevel=2*(n*Kt-St)-Alevel*Kl+Sl.
SumElevel over256 strata andbothfactorgroups to obtainEcenter. KeepEoldodd=2*(n*Ko-So),andsetEodd=Eoldodd+Ecenter+Eoddnear,E=Eeven+Eodd,G=G0+E. All counted slopes are distinct by8n^2<r anddisjointresidue/stratumselection. Empty rounded subintervals,if any,contributezero. The sharedledger reconstructs every stratum,verifies its center/distance/partitionbounds and exactseriescounts,thenchecks the full-size totals:
Eoldodd=292825572424985317226309440173250845152957862724982218878702114059776958264.
Ecenter=480611197610300088285834560176632441790010274894454170797471120988699699584.
All work remains the same concretealgorithm;stratification adds no production operation.

## Odd gaps near odd multiples add a third disjoint subset

The preceding odd subsets all have distance at leastWodd from every odd multiple ofu0. Now select odd gaps Jgap=2u0*q+u0+delta with0<|delta|<Wodd. Sinceu0 is even,Jgap anddelta are odd,so delta is never zero. Write d=|delta| and ell1=epsilon*(2q+1). The original integer conditions imply2n+u0<lambda0 and q>=0,2q+1<u0.

First exclude D1 old-family overlaps. The equality has H=T-2Aold,|H|<2n,so reduction modulo lambda0 forces H=ell EXACTLY,with no v carry. Since|Bold|<n and|ell|<u0,the distance |epsilon*Jgap-ell*u0| must be<Wodd. As2Wodd<=u0,the ONLY possible quotient is ell=ell1;all other odd quotients are at distance at least2u0-d>Wodd. Thus2Bold=epsilon*u0*delta+ell1. Its sign is epsilon*sign(delta),becauseu0*d>|ell1|. The actual old D1 coordinates obey: ifBold>0, Aold+Bold<n; ifBold<0,Aold>=-Bold. Substituting2Aold=T-ell1 shows positiveBold requires T+epsilon*u0*delta<2n,whereas negativeBold requires T>=-epsilon*u0*delta. Consequently positiveBold is excluded by T>=2n-u0*d and negativeBold by T<u0*d.

Next exclude D3 old overlaps for the ENTIRE new distance band. The exact equation from the odd-gap section is2Bold=u0*(epsilon*3Jgap-ell*u0)+ell+v,with oddell,|ell|<u0 andv in{-1,0,1}. An overlap needs |epsilon*3Jgap-ell*u0|<(2n+u0+1)/u0<Wodd. Relative to the odd multiple3*(2q+1)u0,the distance is3d. Every other odd multiple is at distance at least2u0-3d>2u0-3Wodd>=Wodd,because2Wodd<=u0. If3d>=Wodd,there is therefore no possible overlap quotient. If3d<Wodd,the ONLY possible quotient is ell=3ell1. This is divisible by3. The already-proved old D3 coordinate congruence v=ell modulo3 forcesv0. The exact2Bold equation is then divisible by3,so Bold=0 modulo3,contradicting the old selection i!=k modulo3. This excludes every D3 overlap;no carry or sign case is dropped.

Count a conservative subset of the rigorously safe D1 centers. Use256 disjoint distance strata dmin=floor(l*Wodd/256),dnext=floor((l+1)*Wodd/256). Skip l0 with D0. For l1..255,setD=u0*dmin,and select only odd residues in[u0-dnext+1,u0-dmin] and[u0+dmin,u0+dnext-1] modulo2u0. These intervals are strictly inside(u0-Wodd,u0+Wodd),disjoint from EVERY earlier odd subset,and from all other strata. The actual d satisfiesu0*d>=D. For the phase with negativeBold choose T<D-2;for positiveBold choose T>=2n-D. These are safe subsets of the endpoint exclusions. The extra2 in the low direction makes the two phase counts EXACTLY equal under reflection T to2n-2-T on T=Jgap+2i,0<=i<n-Jgap.

For a fixed gap,the combined count over both signs is max(0,min(2n-2Jgap,D-1-Jgap)). This follows because the low direction has max(0,min(n-Jgap,(D-1-Jgap)/2)) centers;Jgap is odd and D even,so the division is integral. Let K(z),S(z) count/sum the stratum's odd gaps1<=Jgap<=z modulo2u0,using the same exact odd arithmetic series as above. If D-1<=n,the stratum contribution is(D-1)*K(min(n-1,D-1))-S(min(n-1,D-1)). Otherwise,withz=min(n-1,2n-D+1),it is2n*K(n-1)-2*S(n-1)+(D-1-2n)*K(z)+S(z). This is the exact linear split at Jgap=2n-D+1,with zero endpoint weights harmless. At these full-size integers1<D<2n-2 for all retained strata.

Sum these255 contributions to obtainEoddnear. SetEodd=Eoldodd+Ecenter+Eoddnear,E=Eeven+Eodd,G=G0+E. All new nonzero slopes are pairwise distinct,and distinct from every even slope,by the already-proved8n^2<r lattice argument. Residue-band separation proves disjointness from both preceding odd subsets. Thus their guarded events have the same exact probability(r-1)/r^2 and zero giant/giant intersections. Both order-four phases and every matching pair are ALREADY tested;NO runtime operation, phase, filter, payload or table is added. Only the analytic event lower bound and public full-sized stream limits change.

Eoddnear=376950853962543075688944830890176775103508690797043554560186093790061291316.

Each giant event has probability(r-1)/r^2 because b!=0. Different giant slopes are disjoint under that guard. Their points are nonidentity, since c_sigma-i=sigma*(i-k)/(1-sigma)!=0. Scalars bi,bk are distinct. A mixed/giant event pair is either an independent linear system with probability at most1/r^2 or parallel and disjoint, because the mixed intercept c is nonzero. These are exact two-variable finite-field facts, not collision independence or a hash heuristic.

Let
mu=M/r+G*(r-1)/r^2,
U2=(choose(M,2)+M*G)/r^2.
The omitted giant/giant intersections are exactly zero. The second Bonferroni inequality gives useful-event probability at least mu-U2. Here
s=344661270701234446760740723829136836316,
M=19216254520450747283203518704951863191033713296879221196438829355679069765632,
G=7505589110893057413229574255455615834400995932853629903672181013604640612604.
All selected events correspond to actual records; inverse phase orientation is harmless because the algorithm tests every required Frobenius phase.

For confirmations,distinct giant exponent pairs are uniform Fr^2 over uniform a,b,and a mixed giant exponent is uniform Fr. True orbit matches therefore have expectation at most12/r per such pair. Deterministic baby duplicates contribute exactly floor(m/u0). For ANY distinct S192 strings the independently sampled raw-sum fingerprint has collision bound delta=(2^-63+7/q)^3,q=2^61-1,proved self-contained in nh_certificate.md. Conditional on complete point sequences,the hash keys are independent;no independence between different record pairs is assumed.

The expected number of ALL equal-fingerprint pairs is bounded by choose(L,2)*(12/r+delta)+floor(m/u0). Thus Markov's inequality bounds cap-abort probability by
caploss=(choose(L,2)*(12/r+delta)+floor(m/u0))/(floor(L/2^29)+1).
This pays false hash collisions,ignored baby pairs and all true matches;there is no free full-key comparison. The giant-identity union bound is n/r,and capped sampling/extraction loss below2^-970. Identity targets are solved directly;unused coins may be extended hypothetically for the unconditional probability argument.

The same exact group-event mu,U2 above give
mu-U2-caploss-maskabort-n/r-2^-970 >3900023/10000000>0.39.
The exact final bound also subtracts maskabort=8L*(50e)/((p-1)K) defined below;the preceding unmasked decimal diagnostic is not a current success claim.

## Complete work and storage

Current scaled implementation: FIELDnormal2345/rawDAGmax2421=2345+4*19; KEYnormal810/rawDAGmax846=810+2*18. Normalretains2386=2345+caller20+spare21 and831=810+caller12+spare9. Main720 givesnormalgeneration3937L. Local720/scan128/confirmation1/inner256/network886/outer80give6008L,PLUSdeterministicglobalrare1L=>N=6009L+12m+28n+2^124. EXACT N^5000<20000^5000*2^621303 proves621303/5000=124.2606. Scalar12m+28nchargedEVERYpath/lastunusedincrement. ALLoldmargins/referenceprices/full381p/255r/threshold.39/sixtracks/solemetricunchanged. Globalrarebound44*(2K)<L,NOTper-callworstcaseoraveragecost. Complete scaling/abort/probabilityproof below andinproof.md.

A complete processed fingerprint pair costs less than2^29 in the normalized scalar-work numerator C*H+R. Binary exponentiation of a255-bit scalar uses at most512 dense reference Fp12 multiplications,including squarings implemented as multiplication. At most12 phases give H<=6144. Public phase inverses are constructed cold. Per-phase Fr products can use256 bit-serial steps on five64-bit limbs,at most128 scalar instructions per step:shift/add an operand belowr,then at most two trial subtractions since the intermediate is below3r. Three such products,adds and phase bookkeeping are below2^20 per phase. Additional exponent-loop/call/buffer/control work is at most256 scalar instructions per dense operation;metadata/payload RPCs,comparisons,pair enumeration and all fixed work fit the generous remainder of R<=2^25. Thus C*6144+2^25=156434432<2^28<2^29 at protectedC20000. No production reference price is changed. Correct terminalG2 output is at most once,cold paid. Since cap=floor(L/2^29),ALL processed confirmations together cost at most L in this numerator,exactly the retained1L above. Detection/termination beyond the cap has no exponent loop and is cold paid. There are no uncharged failures or twelve-global-check assumptions.

The complete two96-bit sort construction is below2^121,including BOTH full2^96 bucket sets,all queues,prefixes,zeroing,flags,padding restoration,initial/tail exceptions and bounded finish. All processor code/page/bootstrap remains below2^121,edge/controller initialization below2^118,and other field/table/coin/phase/output work within its fixed allowance. Their sum remains below2^124. See current radix/architecture/metadata certificates for explicit arrays and transfers. Generic data capacity2^136,private2^119,edges2^111 and algorithm-specific data<2^110 give peak<2^138 bytes,advice0. Unused capacity supplies no answers. All actual pointers are below2^47,IDs96-bit and counters128-bit ordinary words.

The low-limb sign normal form and b1-first representative order in key_contract.md preserve exact twelve-Frobenius fibers. The hash is proved universal for arbitrary distinct192-byte S strings. Discarding stored S changes neither true useful events nor correctness,because every candidate exponent is checked in GT;only the charged false-fingerprint workload and cap calculation change. The scoped exponent-equality oracle in audit.py is TEST-ONLY algebra/control evidence,not an oracle available to the production attack. Shared llm-verifier alone executes candidates,reviews and finalizes;all full parameters,0.39,C20000,six tracks and sole metric remain unchanged.



Current refinement uses256 distance strata in EACH of even endpoint,odd signed endpoint/factor,and odd-near endpoint selections. Each is a refinement of the preceding64 partition: the same conservative lower endpoints are tighter and all exact disjointness proofs apply. It changes no production operation or protected parameter. All public m/n/counts/score are rebuilt; no predecessor acceptance is inherited.

Currentproductionhashuses9independentFqcoefficients/threeROWS,computesFULL128NHthenCOMPRESSESLOW125(C0/C1 32bits,D61bits);top3DERIVEDNHbitsnotconsumed. Protectedcryptographic381p/255r/actualinputdistribution/workm/n unchanged. Newtwo-wrap-regionlinearcongruenceproofgivesNHprojectioncollision<=2^-61;affine2randomcoordinates/directDtwo-casebound5/q;EXACT2^-61+5/q<=2^-63+7/q,soOLDdelta/Markovloss/successremainconservative. NewmapNOTbyteidentical/notall128NHbitsusedbycompression. ALLSIXfull2^32tables192GiB/workerconstructed/charged2^48/memory<2^40atPRIVATE[1,1.1875)TiB. Sample696(base648+Lambda48)/nativepointer624bytes,fairbits1218048<2^21/Rcap2^24 unchanged. Completeproof/native/cold/addresses/projectionaliasfixtureinnh_certificate.md.

## Shared real-coordinate mask and bounded rare work

## Nonzero shared scaling and complete rare-reduction bound

Choose ONE uniform NONZERO canonical native Fp Lambda in[1,p-1] for BOTH workers,independentofb/hash/challenge. NativeRm=2^384;normal lambda=Lambda/Rm. ColdMontmultiplyALLFOURcoordinatesofALLTHREEinitialPhi histories byLambda:storednative=Phi_normal*Lambda. Keepactual c=3Phi(g) coefficients UNCHANGED. Homogeneous recurrence cX-cbarY+Z therefore maintains EXACT lambda*Phi_next;NOvectorbias/translation/unmasklookup is used. Every input/tableentry remains canonical; signed382/9p/49p proofs apply. ALLfullunbiasedvectorentries/both2^96bucketsets are constructed/paid, geometry73.875TiB unchanged. Same sampler696=base648+Lambda48/nativeNHpointer624;1024*381 added fairbits,total1218048<2^21<R2^24. Reject0 OR>=p,accept(p-1)/2^381>3/4,failure<2^-2048;combinedoldcoins/extractionfailure<2^-970. Cold24seedMontproducts/counter/code/distribution/rescale<=2^32 addedoutsideRinside2^124;noauxpowers/advice/amortization.

Fp scaling commutes with every Frobenius action. The four trace candidates of lambda*Phi are exactly lambda times the original four trace candidates. Write M=s[[1,1],[1,-1]],2s²=-1,so M²=-I; q(B) normalizes sign. q(MB)=q(B) would imply MB=+/-B,which forces B=0 in odd characteristic. Thus the same lexicographic comparison selects a unique orbit representative for every scaled nonidentity trace. Canonical scaled tuples equal IFF the underlying original traces differ by these Frobenius actions,because Lambda is invertible;the original cubic trace theorem gives EXACT12GTorbit fibers. The selected representative MAYchange withLambda:we do NOT claim canon(lambda*Phi)=lambda*canon(Phi). IdentityB0 preserved;useful-eventmu/U2/lattice/duplicatecounts and all12phase exponent tests unchanged. UniversalNHcollisionbound applies conditionalLambda to arbitrary distinct192-bytekeys;actualGT exponentverificationcontinuesfalsematches,finalG2/pairingoncecold.

Let Dred=floor(p/2^325),e=p-Dred*2^325>0. AnyrawsumS<49p has qhat=floor(S/(Dred*2^325)) in{q,q+1},q=floor(S/p). Trialborrow implies R=Smodp>=p-50e. For every fixed unscaled FIELDcoordinate v,storednativeR=v_normal*Lambda is uniformFp* ifv!=0;v0 givesR0outside tail. Each of FOUR field borrows has probability<=epsilon=50e/(p-1). The TWO key fixed-odd-B outputs after first sign normalization equal sigma(Lambda)*lambda*Cj for fixedCj and sigma∈{+1,-1} dependingonLambda. Do NOT assume uniformityafterfold:the event lies in the UNION of +lambda*Cj and -lambda*Cj tail events,probability<=2epsilon each (or0ifCj0). Hence total perrecord<=8epsilon. On full hypothetical uncapped m/n streams Etotal<=8L*epsilon bylinearity,NO independentrecords/GTcoordinateequidistribution. Private128counter increments EVERYFIELDandKEYborrow,haltsimmediatelyifcount>=K=floor(L/2^21)+1;actualcount<=K eachworker. Anyabort implies hypotheticaltotal>=K,so maskabort<=8L50e/((p-1)K). Subtract EXACT rational maskabort alongside caploss,n/r,coinloss;still>.3900023>.39. No conditional survival/uniformity assumption or reweighting.

AllFOUR fieldreducers use14-op7-limbtrial+B.HS15normal;onborrow10-opcorrectedlow6limbs plus9-op128counter/cap=19extra. Storetrial11..16direct. Abort ORRout0lowbit1/Bcommoncalleeepilogue;actualcallerTBNZ0,0 failure BEFOREanyread/rotation/key. Laterabort mayleavecomputedcoordinatePREFIX;remainingcoordinatesretainseed,neverconsumed. Callerargs3+tag1+ring4+128progress6=14<=20. BOTHkeyfixedreducers similarlyuse14-optrial+B.HSnormal15;onborrowcorrect10+counter8=18extra (ctx26alreadylive). KeyabortMOV7zero/LDP21,22SP32/Bouterepilogue restoresallcallee/SP/LR;parentFULL64CBZ7 aborts ENTIREattack onzero (identity OR cap),never skipsrecordandcontinues. Nooutputleafwrites;metadata_register_packonzero writesnothing. Counter128/136,K144/152 disjoint full128limit96/104 andkeys112. Everyactualborrow<=44extra includingearlyabortrelative tonormal allowances;44*2K<L chargesALLrareworkONEglobalL on EVERY execution. Normalkeycaller<=12 unchanged. No free error polling/status/helper work.

Abort cleanup uses ordinary bounded control messages,not a free global primitive:Root stops issuing source RPCs,drains/cancels bounded active callback chain,tellsbothsources tostop,then traverses the entire initialized tree sending STOP and collecting ACKs. Passivequeueddata are discarded by bounded head/tail resets,notprocessed. Tree<2^97nodes/edges,atmost2^16ordinary64-bitlocal/message/dispatchoperationsperedge=>work<2^114;reserve a SEPARATE2^118abort allowance. Existing2^118edgeinitialization remains. Normal ENTRY/continuationdispatch unchanged,errorselects this paidprocedure only;actualTBNZ/CBZalready charged. Cold2^121+2^121+2^118+2^118+2^92+2*2^80+2^48+2^16+2^32<2^124,checkedexactly. No hardware/device/free scheduler computation claimed.

Proposed time_log2=621303/5000=124.2606. Public integerN=509542133971337029355827860562544485728256;N^5000<20000^5000*2^621303.

## Wider key digits with full-sized smaller worker mailboxes

Production fixed-s widths[35]*7+[34]*4 coverALL381bits inelevenfulltables,nottruncation. Balancedtable13.5TiB ends91.40625;workerschildbase92TiB/9KiBports/two4608halves end128TiB. Fullpayload4096+header32+destination16<ready4544,doorbell4552;bothhalves/ALL2^32ports/nativeaddressarithmeticsafelydisjoint. Othernodes retain12KiB. Allserviceandconstructionmarginsremain;actualold/newportaddressfixturesboth9ops. Fieldtables73.875TiB unchanged. Newlookup137ops vsold143 beforecommonreducer removes6/body,12/key,includingALLextraextract/load/stridecost. Allsixinputlimbs/fullcanonicalp/latecarry/rareandcap/ABI remainchecked. Y<11p withinold49p/tail50e bound;scalingprobabilityproofunchanged. Normalkey810/rawmax846/ret831;N=6009L+12m+28n+2^124. Standalone192-opold12digithelper isUNCALLEDreference withseparatevirtualtestlayout;productiononlyusesnew11digitinlinebody. Fullentriesconstructed/paid2^92,notVMallocationorqualification.


## Exact paired-input refinement

The FOURwidefieldcolumnreadersuseTHREELDPsinsteadofoneLDP/fourLDR,keepingthemiddlehighwordinextra22. See fp2_certificate.md for everycross-limbconsumer/liveness relation. Eightloadsremoved,minustwoPAIDSTR/LDRpreserving22inSP64,netSIX. Commonabortepiloguerestores22also;frame528/rawinputbytes/fieldwidths/tables/coinproof/Markov/failurecostunchanged. Normalfield2345/MAX2421,normalret2386=2345+caller20+21spare. Key810/ret831/Main720=>3937normalgeneration. Other2071=>6008L,PLUSsameglobalrare1L=>N=6009L+12m+28n+2^124. Allmargins/referenceprices/threshold/fullparameter sizesretain. Thisisrealnativeinstructionremoval,not unused-margintrim.



Current ordinary-column paired-input refinement:fourCOLUMNs eachreplaceSIXLDRwithTHREELDP,saveTWELVEactualfieldinstructions. Compile-timeMSTEPsource19/22preserveseveryfull382-bitdigit/term;22preservationalreadypaidSP64andcommonabortepilogue. Completeinput/flag/livenessproofin fp2_certificate.md. Fieldnormal2345/MAX2421/ret2386,key810/ret831/Main720/other2071/globalrare1L giveN=6009L+12m+28n+2^124. Independentoracles/nativeDAG/count/ABI/forcedabortguardremain;no acceptanceinherited.

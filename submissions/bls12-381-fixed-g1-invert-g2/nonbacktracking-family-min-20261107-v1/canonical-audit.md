# Complete canonical arithmetic,scalar dispatch andhardcap

## Exactcandidatearena initialization

Fromx,y=T(x),b=beta*x,c=beta*y,writecandidatek0=x,k1=y,k2=b,k3=c. CopyxwithsixLDP/sixSTP plus2pointerinstructions,14total. Otherthreepairsarewritten directlybythe6nativeFpMplus2structuredfieldops. Derivek4=b-x,k6=-x,k8=-b,k10=x-b,k5=c-y,k7=-y,k9=-c,k11=y-c. EightFp2assignmentsusetwofieldops each,16total. Everycandidatebyteusedbytheminimumisinitialized. Theinitialcopy14+minimum234+interfaces16=264<272.

# Fold conjugation into signed structured constants

This is an exact scalar implementation refinement of the same canonical function. It changes no fieldvalue,canonicalorientation,hashlabel,trajectory,probabilitymodel orworkcap. Allconstants aretarget-independent publicfieldvalues constructedcold;theirpositive andnegative Montgomeryencodings arebothstoredandcharged.

Let the rawinputpairbe(a,b),andletq indicate whether conjugation isrequired. The following formulas compute c*(a+(-1)^q*b*u),with u²=-1:

- c=(s,0):even(sa,sb),odd(sa,(-s)b).
- c=(0,s):even((-s)b,sa),odd(sb,sa).
- c=(s,s):even(s(a-b),s(a+b)),odd(s(a+b),s(a-b)).
- c=(s,-s):even(s(a+b),s(b-a)),odd(s(a-b),(-s)(a+b)).

Each uses two native Fp multiplications. Real/pureimaginarycases need no otherfieldoperation;diagonalcases needexactlytwo canonicaladd/suboperations. Selectcoldconstant s or-s asindicated;do not materialize -b ornegate the output afterward. Allsign/conjugationeffectsareinthemultiplicationconstants oradd/sub choice. Input/outputpathsuseordinaryscalarcallsandpaidconstant/pointerselection,notfreefieldautomorphismprimitives. Intermediatebuffersandcanonicalwritesarechargedinthe unchangedkey/control256 andwrapper256 allowances. AllnativeFp384 interfacesand42-instructionfield-operationpricesremainunchanged.

Onlythree such productsoccur:T(x),andtheother twob-pairtransforms. Theythereforeneedatmostsixfieldoperations,replacingthepreviousnine. For thecheapk0/6branches,theother twopairsinsteadneedfourcopy/negationoperations;T(x) needsatmosttwo,so6stillboundsallbranches. Otherevenk needsatmosttwo b1signoperationsplusT(x)'stwo,andtheotherrealconstantproduct needszeroadd/sub. Candidate-generationrotationsstilluse16operations,withtheexplicit96-byteinitialcopyinside272minimum/copy asbefore. Noothercopyornegationisomitted.



# Exact two-family minimum over the same twelve candidates

This replaces only the minimum-selection implementation. It returns exactly the same minimumFp2pairandappliedk as the fullreal-then-imaginaryunsignedorder. It introduces no failureevent,fieldarithmetic,hash/trajectorychange,probabilityhypothesisorfreeinterface. canonical-min.S iscomplete static AArch64evidence,neverhostexecuted.

The12imagesaretheevenfamily beta^m*x andoddfamily beta^m*y,m0..5,y=T(x),beta primitiveorder6inFp. Let a=real(x),c=real(y). Ifa!=0,thesixrealcoordinatesbeta^m*aaredistinct;likewise forc!=0. Ifbotharenonzero,findeachfamily'sminimumbycomparingonlyitsrealcoordinates. Eachfamilyhasuniquerealminimum. Theglobalminimumisoneofthesetwowinners;comparetheirreals,andonlyifequalcomparetheirimaginaries. Thisisexactlexicographicselection,nottheheuristicomissionofatiecase.

Ifa=0,allevenrealsarezero. Theoddrealsarenonzero:theFrobeniusconstantford1/5isnonzerodiagonal,soT(a+bu)hasrealpart +/-s*bwhen a0,withs!=0. Sincex!=0,b!=0,real(y)!=0. Hencezeroiseverywhereuniquelysmallestunsignedrealvalue,andtheminimumistheevenfamilymemberwithsmallestimaginarycoordinate. Thoseimaginariesarebeta^m*b,distinct. Similarly ifc=0,onlytheoddzero-real familycanwin;itsimaginaryvaluesaredistinct. Botha=c=0wouldforcex=0,whichisoutsidecallercontract. Noexceptionalvalueorprobabilitylossisignored.

Theassemblyimplementsallsuchbranches. Zero-realtestsuse3LDP,5ORR,1CBNZ=9each. Acomponentcomparisonuses3pairedloadsperoperand(6),oneSUBS/fiveSBCS(6),candidatepointer/ksetup2,andpointer/kCSEL2,total16. Eachfamily'sfivecomparisonscost80. Thetwo-winnerrealcomparisonuses12data/arithmeticinstructionsplus5nonflag-settingORRtoaccumulatedifferences. CBNZpreservescarry. Onnonzerodifference,finalCselectsthesmallerreal;onequalreal,theimaginary384-bitborrowchain12instructionssetsfinalC. TwoCSELselectpointer/k. Worstfinalcomparison32. No finalZbitisusedasmultiwordequality.

EntryfourMOV,normalzero-tests9+1pointer+9,familycomparisons160,saving/resettingwinner4,finalcomparison32,andoutput15give234internalinstructions. Theeven-zero path109,andodd-zero path121,arebothsmaller. Output15includessixLDP/sixSTPcopyingthefullwinning96bytes,STRk,MOVreturnk,RET. Allpairedimmediates0..80arealigned/inrange,andcandidateoffsets96k<=1056fitADD. Itusesonlycallerregisters,x18untouched,nointernalcall,frameorordinaryspill. Minimumoutputandkdestinationsaredisjointfromcandidatearena. Long-livedcontroller/programregistersremainpreserved.

Theexisting16caller/interfaceand14initialxcopyallowancesarestillcharged. Total234+16+14=264;retain272,eightslack. Thisreplacestheformer384minimum/copyallowancewithoutchanginganycandidate-generationproduct/rotation,fieldprogramcost,192key/controlor256wrappers. No per-Fp2memcmporinventedwideprimitiveisused.



# Explicit field-call dispatch and256 canonical key/control instructions

This is a concrete implementation of the existing formulas with fixed ordinary scalar descriptors. There is no computation oracle or executable input from the candidate author:the public formulas specify every descriptor,andthe algorithm constructs them cold from its initialized shared addresses andpublic constants. Candidate evidence itself isneverhost-executed. Thegroup/canonical/hashvalues andall failureprobabilities areunchanged.

## Fixed initialized shared buffers and programs

Use fixedshared RAWB288 bytes,CANONB288 bytes,candidatearena1152 bytes,andboundedfieldtemporaries. Before eachcanonicalizationcopycurrentrawb intoRAWB. Every descriptoroperand points into these fixed shared buffers orcoldconstants;no bankedwindow address isretained. Candidate outputs arewritten directly beforeuse;unallocated/unwrittenbytes areneverread. The existingphysicalshared2^45 allowance easily includes these buffers,programs,framesandconstants. Everydescriptor/code/constant isconstructedandchargedcoldwithin theunchanged2^120 completepriorcoldterm.

There aretwo seedprograms,selected bywhetherb0 or b2 is nonzero. Each has24 fieldcalls:four forT(x) (twoadd/sub,twoFpM),four forbeta*x/beta*y (fourFpM),and16rotationfieldoperations. The initialxcopytocandidate0 remains separatelypricedinside272minimum/copy. y,beta*x,beta*y andeightderivedpairswriteallothercandidates. The authoritative lists arecanonical.md andthisaudit.

There are24 remainingprograms,selected byappliedk0..11 andselectedbpair0/2. Theywrite the other two b pairs into CANONB. Oddk hasatmost8fieldcalls (fourFpM,fouradd/sub);cheapk0/6 hasfourcopy/negationcalls;otherevenk hasatmostfourcalls. The minimumpairwasalreadywrittenintoits CANONB slotbycanonical-min.S. Thus atmost32 fieldcalls occuracrossthetwoprograms. Programlengthsarepositive fixedsmallintegers,notruntimeunknown loops.

Each32-byte descriptorstoresfourordinary64-bitwords:functionaddress,arg0,arg1,arg2. FunctionaddressesarethepinnedFpM orprovidedcanonicaladd/sub/neg/copy leaves. Encodearguments intheactualhelperABI:blstFpM uses(output,a,b),field-ops.S add/sub use(a,b,output),copy/neg use(input,output,unused). NoadapterfunctionorassumedcommonABIisneeded. Unusedarg2iszeroandneverreadbythecopy/neg leaf. Thefixeddescriptorarrayisordinaryinitialized data,notasearchedtableorfreeadvice. Its readsanddispatch arepaid.

field-program.S gives thecompleteinterpreter. TwoLDP loadthefourdescriptorwords andadvancecursor;BLRcallsthefixedfieldfunction;SUBS/B.NE advancethepositivecount. Five scalarinstructions percall,plussevenentry/exitinstructions. It preservesx19,x20,originalLR and16-byte stackalignment;allcalledleavespreservex19/x20 undertheirordinaryABI. x18untouched. NoouterLR/calleecontextislost. Atmost32callsandtwo interpreterentriescost5*32+2*7=174;twocallerinterfacescost<=6,total180. Parentcallee-preservation/frameandanyfieldwrapperworkfittheretained256 wrapperallowance. EachnativeFp384 andfieldoperation42 allowanceisstillchargedseparatelywithitsowninterface,so some overheadisconservativelycountedtwice. Nothingisfree.



## RemainingcanonicalcontrolandparentABI

DirectCANONBtabulationremovestheformer37-instructionkeycopyanditsfourcallerinstructions. Remainingcontrolisrawbcopy39,zerotests/identitychoice40,seedprogramselection10,remainingprogramselection16,parentcontext16,finaldispatch12,total133. Retain192,59slack. Nativearguments,fieldprogramdispatch,minimum/initialcopyandclasscounterareseparatepaidbudgets. Theparentneedsx19record/x20modeonly;preservebothandoriginalLRin32-bytealignedframe. Entry/exitatmost6instructions,reserve16. Twofieldprogramentriesand32fieldcallswithcallerinterfacescost180;topinterface<=16. Total180+16+16=212<256wrappers. EachdescriptorencodestheactualhelperABI:the5-instructionnativeadapter(output,a,b),fieldadd/sub(a,b,out),copy/neg(input,out,unused). NohiddenadapterorcommonABIfictionisused. fixedsharedRAWB/CANONB/candidate/temporaryaddressesareinitializedandchargedcold.



# Exact orientation-dependent canonical cost and a hard total cap

This refines the same canonical key without changing its value. There is no additional heuristic. The full worst-case helper bound remains5484, but cheaper Frobenius classes are implemented explicitly and a hard counter bounds their total work. Any budget exhaustion is an algorithmic failure whose probability is proved below.

## Cheap classes and exact upper cost

The selected pair is b0 or b2, so the other two pairs always include b1. The first stage still computes y=T(x),beta*x,beta*y andall twelve candidates,using six Fp multiplications. The minimum determines k before the other two b pairs are transformed.

In the native degree-six extension over Fp2, w³ has degree two: its square lies in Fp2 and it cannot itself lie there without reducing the degree of w below six. Thus w³ generates Fp4 and sigma²(w³)=-w³. Consequently sigma^(2m) acts on the b1*w³ term by(-1)^m, with no Fp multiplication. For k=0, all b pairs are copied. For k=6, the quadratic Fp12/Fp6 involution negates every b pair. At other even k, only the remaining b0/b2 pair may need two Fp multiplications. Odd k may need four for both remaining pairs.

Copies and canonical negations use at most the existing four42-instruction field-operation allowances per structured pair. They do not introduce an uncharged copy. Theexplicitformula/program/minimum/controlboundsinthisauditareused. Define cost class c(k)=0 for k=0,6; c(k)=1 for other even k; c(k)=2 for odd k. The helper is bounded by

 f(k)=3948+768*c(k).

Here3948=6*384+1644, where1644 is the explicit nonmultiplication allowance incanonical-audit.md (six structured fieldops,sixteen rotation fieldops,272minimum,192key/control,256wrappers). Each class unit adds two Fp multiplications,768. Thus f is3948,4716,5484 respectively. This is an upper instruction count for an explicit copy/negate/general-product branch, not a new field-operation price.

## Hard counter and bounded work on every path

After initial cold canonicalization,initialize a128-bit remaining class budget B=floor(1027K/768). It fits because B<2K<2^127. After finding each new k andbefore remaining b transforms,subtract c(k),failing on underflow. LDP,SUBS,SBCS,STP andborrow branch,with class-table lookup andinterfaces,fit16 extra instructions per transition. Cold computation of B can use ordinary fixed-limb arithmetic;it is within2^120.

For any completed prefix,extra class units consumed are at most B. Common work is at most3948K. Hence total canonical work is at most3948K+768B<=4975K. A partially executed aborting helper is covered by one fixed5484 allowance in cold2^120;placing the check before optional products makes this conservative. This is a hard work bound,not just an expectation. Other algorithm caps still apply.

## Failure probability from proven stationary mixing

Applied k runs through all twelve values once on every nonidentity orbit, since k(sigma X)=k(X)-1. Uniform raw GT therefore gives each class frequency exactly from the twelve values. Extend g(X)=c(k(X)) by g(1)=2. Its stationary mean is

 mu=4/3+2/(3r), 0<=g<=2.

This is algebraic orbit counting, not an assumption about sampled orientations. The ideal nonbacktracking chain starts stationary. For lag at least T=32,the proven full-state relative kernel error epsilon=2^-127 implies

 |E[g(X_j)|state_i]-mu|<=2epsilon,
 |Cov(g(X_i),g(X_j))|<=4epsilon.

For shorter gaps use absolute covariance<=4. Thus S=sum_(t=1)^K g(X_t) obeys

 Var(S)<=4(2T+1)K+4epsilon K²=260K+4epsilon K².

At the protected K,r, B/K-mu>0.00390624, including floor rounding. Exact rational Chebyshev evaluation gives

 Pr[S>B] <=(260K+4epsilon K²)/(B-Kmu)²<2^-100. Here B=floor(1027K/768),so B/K approaches4/3+1/256. The public exact protected-parameter calculation verifies the displayed strict bound,including floor andidentity correction.

The bound holds for every good generator table because stationarity and mixing do. The actual algorithm follows the same ideal trajectory until a prior stop or this cap; g is nonnegative, so exhaustion on a prefix implies S>B on the full ideal extension. No conditioning on already adaptive fingerprints or successful paths is used. Fingerprint,phase,identity andtable failures remain separately subtracted. Budget decisions use group k only,not coefficients,so denominator uniformity is preserved.


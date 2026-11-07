# Nonbacktracking orbit walk with one joint hash and a capped chained dictionary

Proposed score125879/1000=125.879: complete classical computation below2^(125879/1000) protected dense-Fp12 equivalents at C=20000, success>=0.39 under disclosed standard BLAKE3 bridge H and the explicit finite banked-memory realization. Peak memory<2^136 bytes, zero advice. This includes cold construction, trials, failures, global lookup, transfers, verification and canonical G2 output. It is analytic exploratory evidence, not qualification, promotion, production execution, practical availability or a fixed-hash theorem. Storage admissibility is a material review condition.

## Exact problem and representation

Use the full protected BLS12-381 profile. Validate canonical nonidentity P and GT z. If z=1 return G2 infinity. Otherwise compute h=pairing(Q0,P), a generator for every allowed P. Solve z=h^a and return Q=[a]Q0 after verifying both h^a=z and pairing(Q,P)=z. No auxiliary powers, replacement P, hidden exponent or scalar-only answer is used.

fp4-trace.md proves the exactkey S(X)=Tr_(Fp12/Fp4)(X)/3=a0+b1*w³. Its four Frobeniusimages identify exactly the full twelve-element orbitofnonidentityGT,becauseeachtracefiberis{X,sigma4X,sigma8X}. ThecanonicalkeyhasfourFpcoefficients,192nativebytes,24words. Itisnotanunprovedprojection. Theappliedorientationk is0..3,with k(sigma X)=k(X)-1MOD4. AllrawGTpointsremainfull576-byte referenceoperands;theglobaldictionary stores every orbit.

## Cold generator table and standard hashes

Use one table of width W=2^76, distinct from the tower symbol w. Generate independent pairs(u_j,v_j), M_j=h^u_j*z^v_j for0<=j<W and alpha0. Store all twelve variants sigma^e M_j and coefficients xi^e*u_j,xi^e*v_j. Each640-byte record contains576 native GT bytes and two32-byte canonical residues. Every value is constructed for this target and charged, without advice or amortization.

B3 means standard unkeyed BLAKE3-256, all seven rounds. Tags INIT,SCAL,TABLE,JUMP are the exact16-byte ASCII strings FAPI-RHO-INIT-v1,FAPI-RHO-SCAL-v1,FAPI-RHO-TABL-v1,FAPI-RHO-JUMP-v1. Seed=B3(INIT||P_bytes||z_bytes), with P x||y each48-byte big-endian canonical Fp and z coefficients0..11 each48-byte big-endian canonical Fp. INIT is688 bytes: eleven blocks, first flags1, intermediate0, final48 flags10, chunk counter0.

A global unsigned128-bit counter, encoded16-byte big-endian, advances on every scalar attempt. Candidate=low255 bits of B3(SCAL||seed||counter16), digest interpreted big-endian. Accept iff<r, cap1024 attempts per scalar then fail. This64-byte single-block message uses len64,flags11,counter0. At most2W+1<2^78 scalars and<2^88 attempts fit128 bits. Protected r>7*2^252 makes rejection<1/8; all scalar-cap failure mass is far below2^-970.

Use49 simple-tabulation tables:48 alphabets2^32 and one alphabet{0,1}, entries40 bytes. For i<48,c_i is the low/high32-bit half of key word floor(i/2), even/odd respectively; c48=0. Entry words0..3 are four big-endian numerical64-bit words of B3(TABLE||seed||i_1||c_4||part0); word4 is the first64 bits of part1. Position, digit and part widths are1,4,1 big-endian. Each54-byte message uses len54,flags11,counter0. Construct every entry, including unused T48[1]. U=T48[0] XOR the48 selected canonical-trace digit entries, serialized as five big-endian words. tabulate.S computes this in541 instructions. Pairwise collision probability for distinct exact keys is2^-320 inside H.

## Previous edge and ordinary/retry labels

Initialize(X,A,B)=(h^alpha0,alpha0,0), its canonical key,k,U, and abort if identity. Independently choose a fictitious previous generator uniformly among12W labels inside H. EDGEJ and EDGEP tags are ASCII FAPI-NB-EJ-v1 and FAPI-NB-EP-v1 padded by three zero bytes to16. Hash tag||seed (48 bytes,len48,flags11,counter0). Take low76 little-endian bits for jprev, and floor(12V/2^256) from the independent phase digest for eprev, with bias charged. No initial physical jump is applied. This initializes the extended Markov chain in stationarity. Save X0 in the initial dictionary head/payload.

For each current orbit compute ONLY J=B3(JUMP||seed||zero16||U). Prefix64 iscompressedcold len64 flags1;final40 useslen40 flags10 counter0,fullstandardsevenrounds. Interpret its32 digestbytes LITTLE-ENDIAN asY. Setj=Y mod2^76,V=floor(Y/2^76),t=floor(12V/2^180),bucket=V mod2^128,e=t-k mod12. dictionary.md specifies twoEXTR/oneLSR forV andordinarythree-limb multiplicationby12. NoordinaryPHAShash iscalled. Phase/bucketarecorrelatedandtheconditionalvarianceproofexplicitlyhandlesit. jandVareindependentinsideH.

Reject iff j=jprev and e=eprev+6 mod12, before any group multiplication. On a good distinct table this is exactly the inverse previous generator. On rejection resample using RYJ/RYP tags, ASCII FAPI-NB-RJ-v1 and FAPI-NB-RP-v1 padded by three zeros to16. Cache their64-byte tag||seed||zero16 prefix cold. Final block is U40||attempt_8, attempt1..1023 big-endian, len48,flags10. Use the generic full standard hash; the40-byte specialized kernel is not misused. Retry block is separate, so ordinary block padding remains zero.

Cap each transition at1024 attempts and global additional attempts at floor(K/2^32), failing before either cap is exceeded. Every hash, copy, counter and guard is charged. H treats these actual full standard domain-separated BLAKE3 calls as random functions. This bridge remains unproved at full query scale. Actual hashes are computed and paid, not supplied by an oracle. All new domains, counter16 and76-bit extractions are included. No separate random-mapping, collision-independence, heuristic mixing or relation assumption is used.

## One multiplication and immediate global detection

After accepting(j,e), multiply raw X by M[j,e] once and add its rotated coefficients to A,B modulo r. Set previous generator=(j,e). Recompute canonical key,k,U and abort on identity. Search the global chained dictionary as specified in dictionary.md. On a new U, allocate/write payload, prepend to its bucket and continue.

OnUequality,compareall24canonicalSwordswiththesavedpayload. DifferentSmeansfingerprintfailure:abort. EqualSisanexactfulltwelve-orbitmatch. ReconstructsavedrawXs=h^A_saved*z^B_savedinpaidterminalcoldwork. Thethreepossiblerawtwistsareell=k_saved-k_current+4m mod12,m0..2;testfullcurrentXagainstthesethreeFrobeniusimagesofXs. Exactlyonematches. SetD=B_current-xi^ell*B_saved;ifD0fail,otherwisea=(xi^ell*A_saved-A_current)/D. Verifyh^a=z,constructcanonicalQ=[a]Q0,checkcurve/subgroupandpairing(Q,P)=z,andoutput. Allreconstruction,twisttests,inversion/verification/outputarecoldpaid. Thetruetwistdependsonlyonthegrouptrajectory;coefficientrepresentationsreconstructthesameinvariantXs. Noextrauniform-denominatorassumptionisneeded.

The transition remembers the previous generator. A repeated raw orbit need not start a quotient-only cycle. Therefore this algorithm uses no checkpoint-only detection: it stores every visited orbit and detects the first repeat directly. dictionary.md gives exact full-hash chaining with a hard total traversal cap, a finite cap-failure proof and all initialization charged.

Set R=min{R:12R²>=r}, L=ceil(199R/200), K=L. Execute at most K accepted transitions. There is no delayed cycle traversal. At most K+1<2^126 payloads fit the reserved capacity. The2^128 bucket heads are fully initialized cold;The dictionary useslow128bits ofV;thephaseusesall180Vbits. Theircorrelationisexplicitlyhandled byconditioningonphasebins in dictionary.md. No newhashdomain isintroduced. mixing.md and probability.md prove the finite repeat/fruitless bounds, including the previous-edge state and all failures.

## Complete work and memory

The scalar ledger sums3711K and retains3776K, with65 instructions of slack per transition. It includes canonicalization, tabulation, one ordinary hash, scalar update, the explicit384-instructioncontroller, bucket-head accesses, capped chain traversals, initialized payload writes and retry budget. Exactly K reference operations and fixed cold(540*2^119) give

 W <= K+(3776K+(540*2^119))/20000=(23776K+(540*2^119))/20000.

The exact inequality(23776K+(540*2^119))^1000<20000^1000*2^125879 supports125.879, above diagnostic125.87866896466751. No reference price is changed or multiplication counted twice.

Capacity is payload2^135+heads2^132+jumps2^89+shared<2^45, below2^136 bytes including padding. storage.md uses97-bit bank ids latched by ordinary64-bit stores, explicit commit/barriers, and charged two-limb logical-index arithmetic. There is no flat136-bit unit address or free RAM search. This finite storage model requires review; if inadmissible, the claim must be rejected. It is not installed on the VM. Analytic bounds are permitted and local intake does not allocate production memory.

Cold(540*2^119) includes every previous construction/finish term (<2^120), plus explicit initialization of all2^128 bucket heads. zero8192 executes512 zero stores andRET;with24 caller/address/loop instructions perblock,total537*2^119<2^129. All bank switching/index setup stays<2^100. The total is below(540*2^119),without pre-zeroed memory,advice oramortization. Payloads areinitialized onallocation. dictionary.md proves the additional traversal-cap failure<2^-100;subtract it separately from probability.md's prior losses. The final repeated query is explicitly counted,not treated as a fresh dictionary hash.

## Attribution and limitations

Published23d3d4cd-768d-4cac-a9a5-872fc010580d supplies the original reduction, norm-one framework and interval/Bonferroni method. Published13f0490c-4535-4ee8-8112-6a72a5254024 supplies standard BLAKE3/helper contracts. Exact citations are declared. Source13f was AI-valid127.756, not promoted; it does not qualify this proposal. Local orbit434370e9-55ac-40b0-a90b-83a6bc25fb41 and b-half plain6f72d8b7-b51d-4552-9026-3c9e21c0d7c9 are material unqualified sources, acknowledged rather than called independent. Newtracepolynomial/nativeprojection,nonbacktrackingmixing,moments,fruitlessanalysis,dictionaryandstoragecalculationsareexplicit. ExternalXTRtrace-characteristic-polynomialsourceiscreditedandlinkedinfp4-trace.md;thisspecializationisprovedhere. H and memory-model admissibility remain material. No novelty, practical break or scientific qualification is asserted.

Standard definitions are supported by the pinned primary sources:
https://github.com/BLAKE3-team/BLAKE3/blob/f55849f89cd85083c1c9daa2c5de20766f309c00/reference_impl/reference_impl.rs
https://github.com/supranational/blst/blob/f7a88067ed98385ebc6cf41ccecf21ead706a5b7/src/fp12_tower.c
https://github.com/supranational/blst/blob/f7a88067ed98385ebc6cf41ccecf21ead706a5b7/src/asm/mul_mont_384-armv8.pl

Candidate files are untrusted static evidence, never compiled or executed on the host or with credentials. Shared llm-verifier owns scientific execution, review and scoring; experiments are empty. Root drafts remain add-only. Full parameters, exact supplied P, no powers,0.39 threshold, protected prices, six tracks and one metric are unchanged.

## Deterministictracehelperandcompletecontroller

fp4-trace.md provesS'scubicfiberidentity,nonzeroBcriterion,free4-imagecanonicalkey,actualnativeextractionat0..95/384..479,andcoldthree-twist recovery. trace-min.S handleszero/nonzeroAimag paths131/88internal,plus16interface+18initialBcopy=165<176. TwoFpM384andsevenfieldleaf42operations,176minimum/copy,128controland128wrappersgive1494perstep. Nocanonicalcostcounter/failureisneeded.

trace-hot-loop.S/trace-controller.md giveallordinarycontrol:77main+17jumpselector=94,retained384with290slack,includingphase/index/previous-generator/bank/ABI/counter/loopwork. Eachhelperbodyisseparateandallretryrechecksareintheexistingglobalextra-attemptbudget. New48digittabulation541,nativeFpM/fieldbodies,scalar46,head192,payload224,cappedchain32andretry1givethefullledger3711,retained3776. ColdD=540*2^119stillincludesEVERYheadzeroandallsetup/finish,newconstants/programs,terminalrawreconstruction/three-twistsearch,exactverificationandG2output. Fullparameters,q12collisionmodel,success0.39,C20000andonecomputationmetricareunchanged.

## Finitefourth-orderprobabilityrefinement

Onlythemaximumacceptedtransitionschanges:K=L=ceil(199R/200),R=min{R:12R²>=r}. Code/helper/hash/table/field/controller/storage/price/coldprotocolsareunchanged,andtheactualtrajectoryisthesameprefixasthesource. fourth-moment.md provesstationarymulti-timecharacterboundrho^sum_nonzerogaplengthsforrho2^-29usingaP-bQexpansion,includingprevious-edgedependence. WeightedrotationconstraintsreduceFourierflowdimensionbFbyspanningforestrank. templates.md suppliesallfinitegraphcoefficientsandcoversallpositive-exponent/forceddiagrams.

ThirdlowerV3andfourthupperU4enteranEVENBonferronilowerbound,notacollision-independenceorPoissonassumption. Exactrepeatlower0.39018290721279814,afterALLpriorfailurelossessuccess>0.3901676481908756>0.39. No newheuristic,hashbridge,missingcorrelationorqualificationisasserted. Completework(23776K+540*2^119)/20000 hasdiagnostic125.87866896466751andexactcertificate(23776K+540*2^119)^1000<20000^1000*2^125879,proposed125.879.

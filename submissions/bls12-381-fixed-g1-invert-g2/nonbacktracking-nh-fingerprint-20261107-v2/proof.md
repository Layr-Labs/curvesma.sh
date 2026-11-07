# Nonbacktracking orbit walk with one joint hash and a capped chained dictionary

Proposed score2517/20=125.850: complete classical computation below2^(2517/20) protected dense-Fp12 equivalents at C=20000, success>=0.39 under disclosed standard BLAKE3 bridge H and the explicit finite banked-memory realization. Peak memory<2^136 bytes, zero advice. This includes cold construction, trials, failures, global lookup, transfers, verification and canonical G2 output. It is analytic exploratory evidence, not qualification, promotion, production execution, practical availability or a fixed-hash theorem. Storage admissibility is a material review condition.

## Exact problem and representation

Use the full protected BLS12-381 profile. Validate canonical nonidentity P and GT z. If z=1 return G2 infinity. Otherwise compute h=pairing(Q0,P), a generator for every allowed P. Solve z=h^a and return Q=[a]Q0 after verifying both h^a=z and pairing(Q,P)=z. No auxiliary powers, replacement P, hidden exponent or scalar-only answer is used.

fp4-trace.md proves the exactkey S(X)=Tr_(Fp12/Fp4)(X)/3=a0+b1*w³. Its four Frobeniusimages identify exactly the full twelve-element orbitofnonidentityGT,becauseeachtracefiberis{X,sigma4X,sigma8X}. ThecanonicalkeyhasfourFpcoefficients,192nativebytes,24words. Itisnotanunprovedprojection. Theappliedorientationk is0..3,with k(sigma X)=k(X)-1MOD4. AllrawGTpointsremainfull576-byte referenceoperands;theglobaldictionary stores every orbit.

## Cold generator table and standard hashes

Use one table of width W=2^76, distinct from the tower symbol w. Generate independent pairs(u_j,v_j), M_j=h^u_j*z^v_j for0<=j<W and alpha0. Store all twelve variants sigma^e M_j and coefficients xi^e*u_j,xi^e*v_j. Each640-byte record contains576 native GT bytes and two32-byte canonical residues. Every value is constructed for this target and charged, without advice or amortization.

B3 means standard unkeyed BLAKE3-256, all seven rounds. Tags INIT,SCAL,JUMP are the exact16-byte ASCII strings FAPI-RHO-INIT-v1,FAPI-RHO-SCAL-v1,FAPI-RHO-JUMP-v1. NHK/NHC domains and all target-specific fingerprint key derivation are specified in nh-fingerprint.md. Seed=B3(INIT||P_bytes||z_bytes), with P x||y each48-byte big-endian canonical Fp and z coefficients0..11 each48-byte big-endian canonical Fp. INIT is688 bytes: eleven blocks, first flags1, intermediate0, final48 flags10, chunk counter0.

A global unsigned128-bit counter, encoded16-byte big-endian, advances on every scalar attempt. Candidate=low255 bits of B3(SCAL||seed||counter16), digest interpreted big-endian. Accept iff<r, cap1024 attempts per scalar then fail. This64-byte single-block message uses len64,flags11,counter0. At most2W+1<2^78 scalars and<2^88 attempts fit128 bits. Protected r>7*2^252 makes rejection<1/8; all scalar-cap failure mass is far below2^-970.

Use the nh-fingerprint.md construction on the24 native64-bit canonical S words: five independently keyed unsigned NH128 hashes, each independently compressed by an odd128-bit multiplier to64 output bits. Their concatenation U is40 bytes, serialized big-endian. The complete scalar leaf costs460 instructions. Its65 target-specific full BLAKE3 seed calls initialize1040 shared bytes; no tabulation tables or T48 constant are used. The fixed-pair fingerprint collision probability is<2^-312, derived from published NH and a proved odd-multiplier compression lemma. Coupling failure before the first orbit repeat is<2^-61; conservatively charge2^-55. All seed work is cold paid and the full standard BLAKE3 bridge H is unchanged.

## Previous edge and ordinary/retry labels

Initialize(X,A,B)=(h^alpha0,alpha0,0), its canonical key,k,U, and abort if identity. Independently choose a fictitious previous generator uniformly among12W labels inside H. EDGEJ and EDGEP tags are ASCII FAPI-NB-EJ-v1 and FAPI-NB-EP-v1 padded by three zero bytes to16. Hash tag||seed (48 bytes,len48,flags11,counter0). Take low76 little-endian bits for jprev, and floor(12V/2^256) from the independent phase digest for eprev, with bias charged. No initial physical jump is applied. This initializes the extended Markov chain in stationarity. Save X0 in the initial dictionary head/payload.

For each current orbit compute ONLY J=B3(JUMP||seed||zero16||U). Prefix64 iscompressedcold len64 flags1;final40 useslen40 flags10 counter0,fullstandardsevenrounds. Interpret its32 digestbytes LITTLE-ENDIAN asY. Setj=Y mod2^76,V=floor(Y/2^76),t=floor(12V/2^180),bucket=V mod2^128,e=t-k mod12. dictionary.md specifies twoEXTR/oneLSR forV andordinarythree-limb multiplicationby12. NoordinaryPHAShash iscalled. Phase/bucketarecorrelatedandtheconditionalvarianceproofexplicitlyhandlesit. jandVareindependentinsideH.

Reject iff j=jprev and e=eprev+6 mod12, before any group multiplication. On a good distinct table this is exactly the inverse previous generator. On rejection resample using RYJ/RYP tags, ASCII FAPI-NB-RJ-v1 and FAPI-NB-RP-v1 padded by three zeros to16. Cache their64-byte tag||seed||zero16 prefix cold. Final block is U40||attempt_8, attempt1..1023 big-endian, len48,flags10. Use the generic full standard hash; the40-byte specialized kernel is not misused. Retry block is separate, so ordinary block padding remains zero.

Cap each transition at1024 attempts and global additional attempts at floor(K/2^32), failing before either cap is exceeded. Every hash, copy, counter and guard is charged. H treats these actual full standard domain-separated BLAKE3 calls as random functions. This bridge remains unproved at full query scale. Actual hashes are computed and paid, not supplied by an oracle. All new domains, counter16 and76-bit extractions are included. No separate random-mapping, collision-independence, heuristic mixing or relation assumption is used.

## One multiplication and immediate global detection

After accepting(j,e), multiply raw X by M[j,e] once and add its rotated coefficients to A,B modulo r. Set previous generator=(j,e). Recompute canonical key,k,U and abort on identity. Search the global chained dictionary as specified in dictionary.md. On a new U, allocate/write payload, prepend to its bucket and continue.

OnUequality,compareall24canonicalSwordswiththesavedpayload. DifferentSmeansfingerprintfailure:abort. EqualSisanexactfulltwelve-orbitmatch. ReconstructsavedrawXs=h^A_saved*z^B_savedinpaidterminalcoldwork. Thethreepossiblerawtwistsareell=k_saved-k_current+4m mod12,m0..2;testfullcurrentXagainstthesethreeFrobeniusimagesofXs. Exactlyonematches. SetD=B_current-xi^ell*B_saved;ifD0fail,otherwisea=(xi^ell*A_saved-A_current)/D. Verifyh^a=z,constructcanonicalQ=[a]Q0,checkcurve/subgroupandpairing(Q,P)=z,andoutput. Allreconstruction,twisttests,inversion/verification/outputarecoldpaid. Thetruetwistdependsonlyonthegrouptrajectory;coefficientrepresentationsreconstructthesameinvariantXs. Noextrauniform-denominatorassumptionisneeded.

The transition remembers the previous generator. A repeated raw orbit need not start a quotient-only cycle. Therefore this algorithm uses no checkpoint-only detection: it stores every visited orbit and detects the first repeat directly. dictionary.md gives exact full-hash chaining with a hard total traversal cap, a finite cap-failure proof and all initialization charged.

Set R=min{R:12R²>=r}, L=ceil(1989R/2000), K=L. Execute at most K accepted transitions. There is no delayed cycle traversal. At most K+1<2^126 payloads fit the reserved capacity. The2^128 bucket heads are fully initialized cold;The dictionary useslow128bits ofV;thephaseusesall180Vbits. Theircorrelationisexplicitlyhandled byconditioningonphasebins in dictionary.md. No newhashdomain isintroduced. mixing.md and probability.md prove the finite repeat/fruitless bounds, including the previous-edge state and all failures.

## Complete work and memory

The scalar ledger sums3299K and retains3312K, with13 instructions of slack per transition. It includes canonicalization, NH fingerprint, one ordinary hash, scalar update, the explicit128-instructioncontroller, bucket-head accesses, capped chain traversals, initialized payload writes and retry budget. Exactly K reference operations and fixed cold(540*2^119) give

 W <= K+(3312K+(540*2^119))/20000=(23312K+(540*2^119))/20000.

The exact inequality(23312K+(540*2^119))^20<20000^20*2^2517 supports125.850, above diagnostic125.84951732988598. No reference price is changed or multiplication counted twice.

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

fp4-trace.md provesS'scubicfiberidentity,nonzeroBcriterion,free4-imagecanonicalkey,actualnativeextractionat0..95/384..479,andcoldthree-twist recovery. trace-min.S handleszero/nonzeroAimag paths131/88internal,plus16interface+18initialBcopy=165<176. TwoFpM384andsevenfieldleaf42operations,176minimum/copy,96controland96wrappersgive1440perstep. Nocanonicalcostcounter/failureisneeded.

trace-hot-loop.S/trace-controller.md giveallordinarycontrol:77main+17jumpselector=94,retained128with34slack,includingphase/index/previous-generator/bank/ABI/counter/loopwork. Eachhelperbodyisseparateandallretryrechecksareintheexistingglobalextra-attemptbudget. Newfive-laneNHfingerprint460,nativeFpM/fieldbodies,scalar46,head192,payload224,cappedchain32andretry1givethefullledger3299,retained3312. ColdD=540*2^119stillincludesEVERYheadzeroandallsetup/finish,newconstants/programs/NHseeds,terminalrawreconstruction/three-twistsearch,exactverificationandG2output. Fullparameters,q12collisionmodel,success0.39,C20000andonecomputationmetricareunchanged.

## Finite sixth-order probability refinement

Only the accepted-transition cap changes to K=L=ceil(1989R/2000). Analysis uses stronger delta2^-33 and lagT128; these are not runtime changes. sixth-moment.md derives the stationary multi-frequency bound with rho2^-32, weighted orbit-flow constraints, a structural bound for every free or forced diagram, and finite third-through-sixth factorial moment bounds. The EVEN Bonferroni lower bound gives repeat>0.3900650822981647 and success after every failure loss>0.39004982327624>0.39. No pair independence or new heuristic is introduced. Complete work has diagnostic125.84951732988598 and exact certificate(23312K+540*2^119)^20<20000^20*2^2517, proposed2517/20=125.850.

## Prepared full-round hash

prepared-hash.md proves the supplied hash specialization equals standard BLAKE3 on the current NH messages. The new fingerprint changes messages/actual trajectories relative to tabulation, with coupling loss proved in nh-fingerprint.md. Initialize the prepared JUMP prefix once cold, including before the initial query. Retry-domain CVs remain separate and standard. The776-instruction ordinary final hash retains all56 G calls/seven rounds; cold preparation and all interfaces are paid. The complete ordinary ledger3299 retained3312 gives proposed125.850 with unchanged success proof.

## Small-seed fingerprint

nh-fingerprint.md credits the published NH theorem, derives the five-lane compression/coupling bound, specifies every new seed domain and cold write, and certifies the complete460-instruction leaf. All scalar storage/hash layouts remain unchanged except context48 now points to its1040-byte seed record. The fingerprint failure allowance is conservatively2^-55; the unchanged sixth-order repeat bound minus all revised losses still gives success>0.39004982327624.

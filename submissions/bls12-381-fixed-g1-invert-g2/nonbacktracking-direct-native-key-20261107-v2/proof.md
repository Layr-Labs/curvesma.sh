# Nonbacktracking orbit walk with one joint hash and a capped chained dictionary

Proposed score126170/1000=126.170: complete classical computation below2^(126170/1000) protected dense-Fp12 equivalents at C=20000, success>=0.39 under disclosed standard BLAKE3 bridge H and the explicit finite banked-memory realization. Peak memory<2^136 bytes, zero advice. This includes cold construction, trials, failures, global lookup, transfers, verification and canonical G2 output. It is analytic exploratory evidence, not qualification, promotion, production execution, practical availability or a fixed-hash theorem. Storage admissibility is a material review condition.

## Exact problem and representation

Use the full protected BLS12-381 profile. Validate canonical nonidentity P and GT z. If z=1 return G2 infinity. Otherwise compute h=pairing(Q0,P), a generator for every allowed P. Solve z=h^a and return Q=[a]Q0 after verifying both h^a=z and pairing(Q,P)=z. No auxiliary powers, replacement P, hidden exponent or scalar-only answer is used.

canonical.md proves the exact Frobenius12 orientation C=sigma^k X, sigma(X)=X^p, xi=p modr of order12. It also proves that the b half in X=a+bw is injective on odd-order GT: two distinct points with equal b would multiply to -1, impossible in this subgroup. Theexactkeyis36nativebwordsinCANONB;tabulationfoldsaconstant-zerodigitwithoutanout-of-buffermemoryread. No sign(a), square root or field inverse is needed. Current raw points still store their full a,b for reference multiplication.

## Cold generator table and standard hashes

Use one table of width W=2^76, distinct from the tower symbol w. Generate independent pairs(u_j,v_j), M_j=h^u_j*z^v_j for0<=j<W and alpha0. Store all twelve variants sigma^e M_j and coefficients xi^e*u_j,xi^e*v_j. Each640-byte record contains576 native GT bytes and two32-byte canonical residues. Every value is constructed for this target and charged, without advice or amortization.

B3 means standard unkeyed BLAKE3-256, all seven rounds. Tags INIT,SCAL,TABLE,JUMP are the exact16-byte ASCII strings FAPI-RHO-INIT-v1,FAPI-RHO-SCAL-v1,FAPI-RHO-TABL-v1,FAPI-RHO-JUMP-v1. Seed=B3(INIT||P_bytes||z_bytes), with P x||y each48-byte big-endian canonical Fp and z coefficients0..11 each48-byte big-endian canonical Fp. INIT is688 bytes: eleven blocks, first flags1, intermediate0, final48 flags10, chunk counter0.

A global unsigned128-bit counter, encoded16-byte big-endian, advances on every scalar attempt. Candidate=low255 bits of B3(SCAL||seed||counter16), digest interpreted big-endian. Accept iff<r, cap1024 attempts per scalar then fail. This64-byte single-block message uses len64,flags11,counter0. At most2W+1<2^78 scalars and<2^88 attempts fit128 bits. Protected r>7*2^252 makes rejection<1/8; all scalar-cap failure mass is far below2^-970.

Use73 simple-tabulation tables:72 alphabets2^32 and one alphabet{0,1}, entries40 bytes. For i<72,c_i is the low/high32-bit half of key word floor(i/2), even/odd respectively; c72=0. Entry words0..3 are four big-endian numerical64-bit words of B3(TABLE||seed||i_1||c_4||part0); word4 is the first64 bits of part1. Position, digit and part widths are1,4,1 big-endian. Each54-byte message uses len54,flags11,counter0. Construct every entry, including unused T72[1]. U=T72[0] XOR the72 selected b-digit entries, serialized as five big-endian words. tabulate.S computes this in807 instructions. Pairwise collision probability for distinct exact keys is2^-320 inside H.

## Previous edge and ordinary/retry labels

Initialize(X,A,B)=(h^alpha0,alpha0,0), its canonical key,k,U, and abort if identity. Independently choose a fictitious previous generator uniformly among12W labels inside H. EDGEJ and EDGEP tags are ASCII FAPI-NB-EJ-v1 and FAPI-NB-EP-v1 padded by three zero bytes to16. Hash tag||seed (48 bytes,len48,flags11,counter0). Take low76 little-endian bits for jprev, and floor(12V/2^256) from the independent phase digest for eprev, with bias charged. No initial physical jump is applied. This initializes the extended Markov chain in stationarity. Save X0 in the initial dictionary head/payload.

For each current orbit compute ONLY J=B3(JUMP||seed||zero16||U). Prefix64 iscompressedcold len64 flags1;final40 useslen40 flags10 counter0,fullstandardsevenrounds. Interpret its32 digestbytes LITTLE-ENDIAN asY. Setj=Y mod2^76,V=floor(Y/2^76),t=floor(12V/2^180),bucket=V mod2^128,e=t-k mod12. dictionary.md specifies twoEXTR/oneLSR forV andordinarythree-limb multiplicationby12. NoordinaryPHAShash iscalled. Phase/bucketarecorrelatedandtheconditionalvarianceproofexplicitlyhandlesit. jandVareindependentinsideH.

Reject iff j=jprev and e=eprev+6 mod12, before any group multiplication. On a good distinct table this is exactly the inverse previous generator. On rejection resample using RYJ/RYP tags, ASCII FAPI-NB-RJ-v1 and FAPI-NB-RP-v1 padded by three zeros to16. Cache their64-byte tag||seed||zero16 prefix cold. Final block is U40||attempt_8, attempt1..1023 big-endian, len48,flags10. Use the generic full standard hash; the40-byte specialized kernel is not misused. Retry block is separate, so ordinary block padding remains zero.

Cap each transition at1024 attempts and global additional attempts at floor(K/2^32), failing before either cap is exceeded. Every hash, copy, counter and guard is charged. H treats these actual full standard domain-separated BLAKE3 calls as random functions. This bridge remains unproved at full query scale. Actual hashes are computed and paid, not supplied by an oracle. All new domains, counter16 and76-bit extractions are included. No separate random-mapping, collision-independence, heuristic mixing or relation assumption is used.

## One multiplication and immediate global detection

After accepting(j,e), multiply raw X by M[j,e] once and add its rotated coefficients to A,B modulo r. Set previous generator=(j,e). Recompute canonical key,k,U and abort on identity. Search the global chained dictionary as specified in dictionary.md. On a new U, allocate/write payload, prepend to its bucket and continue.

On U equality, compare all36 canonical b words with the saved payload. Different b is a fingerprint failure: abort. Equal b is an exact orbit match. Set ell=k_saved-k_current mod12 and D=B_current-xi^ell*B_saved. If D=0 fail; otherwise recover a=(xi^ell*A_saved-A_current)/D. Verify h^a=z, construct canonical Q=[a]Q0, check curve/subgroup and pairing(Q,P)=z, then output. Terminal work is paid cold.

The transition remembers the previous generator. A repeated raw orbit need not start a quotient-only cycle. Therefore this algorithm uses no checkpoint-only detection: it stores every visited orbit and detects the first repeat directly. dictionary.md gives exact full-hash chaining with a hard total traversal cap, a finite cap-failure proof and all initialization charged.

Set R=min{R:12R²>=r}, L=ceil(1031R/1000), K=L. Execute at most K accepted transitions. There is no delayed cycle traversal. At most K+1<2^126 payloads fit the reserved capacity. The2^128 bucket heads are fully initialized cold;The dictionary useslow128bits ofV;thephaseusesall180Vbits. Theircorrelationisexplicitlyhandled byconditioningonphasebins in dictionary.md. No newhashdomain isintroduced. mixing.md and probability.md prove the finite repeat/fruitless bounds, including the previous-edge state and all failures.

## Complete work and memory

The scalar ledger sums8035K and retains8064K, with29 instructions of slack per transition. It includes canonicalization, tabulation, one ordinary hash, scalar update, ordinary control, bucket-head accesses, capped chain traversals, initialized payload writes and retry budget. Exactly K reference operations and fixed cold(540*2^119) give

 W <= K+(8064K+(540*2^119))/20000=(28064K+(540*2^119))/20000.

The exact inequality(28064K+(540*2^119))^1000<20000^1000*2^126170 supports126.170, above diagnostic126.1690991273723. No reference price is changed or multiplication counted twice.

Capacity is payload2^135+heads2^132+jumps2^89+shared<2^45, below2^136 bytes including padding. storage.md uses97-bit bank ids latched by ordinary64-bit stores, explicit commit/barriers, and charged two-limb logical-index arithmetic. There is no flat136-bit unit address or free RAM search. This finite storage model requires review; if inadmissible, the claim must be rejected. It is not installed on the VM. Analytic bounds are permitted and local intake does not allocate production memory.

Cold(540*2^119) includes every previous construction/finish term (<2^120), plus explicit initialization of all2^128 bucket heads. zero8192 executes512 zero stores andRET;with24 caller/address/loop instructions perblock,total537*2^119<2^129. All bank switching/index setup stays<2^100. The total is below(540*2^119),without pre-zeroed memory,advice oramortization. Payloads areinitialized onallocation. dictionary.md proves the additional traversal-cap failure<2^-100;subtract it separately from probability.md's prior losses. The final repeated query is explicitly counted,not treated as a fresh dictionary hash.

## Attribution and limitations

Published23d3d4cd-768d-4cac-a9a5-872fc010580d supplies the original reduction, norm-one framework and interval/Bonferroni method. Published13f0490c-4535-4ee8-8112-6a72a5254024 supplies standard BLAKE3/helper contracts. Exact citations are declared. Source13f was AI-valid127.756, not promoted; it does not qualify this proposal. Local orbit434370e9-55ac-40b0-a90b-83a6bc25fb41 and b-half plain6f72d8b7-b51d-4552-9026-3c9e21c0d7c9 are material unqualified sources, acknowledged rather than called independent. New nonbacktracking mixing, moments, fruitless analysis, dictionary and storage calculations are explicit. H and memory-model admissibility remain material. No novelty, practical break or scientific qualification is asserted.

Standard definitions are supported by the pinned primary sources:
https://github.com/BLAKE3-team/BLAKE3/blob/f55849f89cd85083c1c9daa2c5de20766f309c00/reference_impl/reference_impl.rs
https://github.com/supranational/blst/blob/f7a88067ed98385ebc6cf41ccecf21ead706a5b7/src/fp12_tower.c
https://github.com/supranational/blst/blob/f7a88067ed98385ebc6cf41ccecf21ead706a5b7/src/asm/mul_mont_384-armv8.pl

Candidate files are untrusted static evidence, never compiled or executed on the host or with credentials. Shared llm-verifier owns scientific execution, review and scoring; experiments are empty. Root drafts remain add-only. Full parameters, exact supplied P, no powers,0.39 threshold, protected prices, six tracks and one metric are unchanged.

## Direct native key andcompletecostcertificate

ThecanonicalfieldprogramswritesharedCANONB288bytesinordinarynativeFp2/Fplimborder. Tabulationreadsthatbufferdirectly;thereisno36-wordserialization/reversalcopyandnoworking-recordkeybufferisread. Exactpayloadbisequallynative. Thischangesactualfingerprints/hashtrajectories,butunderH leavesallidealdistributionsandfailureproofsunchanged. direct-key.md givesthepreciseword/digit/lifetimecontract.

canonical-audit.md providesallstructuredconstantformulas,thetwelve-candidateschedule,327-instructionminimum,coldfieldprograms/five-instructiondispatch,parentABIand192key/controlbound,andthehardclassbudgetvarianceproof. native-field-audit.md certifies375nativeMontgomeryinstructionsretained384andprovidedfieldleafmaximum42,withstaticprimaryexcerpts/sourceSHA256/license. Remainingcanonicalcontrol133isretained192;wrappers212are retained256. Thusnonmultiplication1756,base4060,classunit768,andunchangedB=floor1027K/768givehardcanonical5087K,capfailure<2^-100. Ordinaryledger8035isretained8064with29slack. Completecold540*2^119andallscientificlossesareincluded;exactbound(28064K+540*2^119)/20000 supports126.170. Nohostcandidateexecutionorofficialscoreisclaimed.

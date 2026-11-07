# Simple b-half key and explicit scalar leaves for the plain GT adding walk

Proposed score2553/20=127.650:complete classical cold work below2^(2553/20) protected dense-Fp12 reference equivalents at C=20000,success>=0.39 conditional on disclosed standard BLAKE3 bridge H. Peak memory<2^50 bytes,zero advice. All preprocessing,attempts,failures,lookup,transfers,verification andG2 output are included. This is analytic exploratory evidence,not a fixed-hash theorem,production execution,qualification,promotion orpractical break.

## Exact supplied-P reduction and jump table

Use the full protected BLS12-381 p,r,Q0 andpinned pairing. Validate canonical nonidentity P andsubgroup GT z. Ifz=1 return G2 infinity. Otherwise h=pairing(Q0,P) generates GT for every allowed P. Find z=h^a andreturn canonical Q=[a]Q0,checking h^a=z andpairing(Q,P)=z. No auxiliary powers,solver-chosen P,hidden exponent orscalar-only output is used.

Tags INIT,SCAL,TABLE,WALK are the exact16-byte ASCII strings FAPI-RHO-INIT-v1,FAPI-RHO-SCAL-v1,FAPI-RHO-TABL-v1,FAPI-RHO-WALK-v1. Seed=B3(INIT||P_bytes||z_bytes),P x||y each48-byte big-endian canonicalFp,z coefficients0..11 each48-byte big-endian canonicalFp. B3 is standard unkeyed BLAKE3-256,all seven rounds. INIT688 is eleven blocks,firstflags1,intermediate0,final48flags10,counter0.

A global64-bit counter starts0,encodes8-byte big-endian,andadvances after every attempt. Candidate is low255 bits ofB3(SCAL||seed||counter8),digest interpretedbig-endian. Accept iff<r,cap1024attempts thenfail. Generate2W+1 accepted scalars for W=2^40:jump pairs(u_j,v_j) andalpha0. Counter<2^52. SCAL56 is one standard block,len56,flags11,counter0. Accepted values are uniform independent insideH;all rejected attempts are charged.

Construct all M_j=h^u_j*z^v_j in640-byte records:576 native GT bytes andtwo32-byte canonical plain residues. Start(X,A,B)=(h^alpha0,alpha0,0). Each step multiplies X by oneM_j andadds its coefficient pair modulo r,maintaining X=h^A*z^B. Exactly one protected reference dense Fp12 operation is charged perstep. No powers or table are supplied free.

## Exact b-only key and fingerprint

Write X=a+bw in the pinned Fp2/Fp6/Fp12 tower,w²=v. GT norm-one gives a²-vb²=1. If X=a+bw,Y=c+bw share b,then a²=c². Distinct points require c=-a andY=-X^-1,soXY=-1 would belong toGT. Impossible:GT has prime odd orderr while-1 hasorder2. Thus b alone is injective on validated GT,although it is not injective on the entire even-order norm-one group.

Store36 numerical64-bit b words plus one constantzero word. Words are the six canonical native Montgomery Fp coefficients,with limb order defined byhot-loop.S copy helper. The redundant37th word retains the original lookup layout. Both working buffers initialize it cold andnever overwrite it;checkpoint copies preserve zero. No sign calculation,square root,field inverse ornormalization is needed perstep.

Use73 independent320-bit tabulation tables insideH:72 alphabets2^32 andone alphabet{0,1}. Digit i<72 is the low/high32-bit half ofwordfloor(i/2),even/odd respectively;digit72=0. TABLE||seed||i_1||c_4||part_1 has54 bytes,with big-endian widths1,4,1. Hash part0 andpart1 using standard single-block len54,flags11. Entrywords0..3 are four big-endian64-bit digestwords ofpart0;word4 is first64bits ofpart1. Construct everyentry,includingunusedT72[1]. U=T72[0] XORthe72selectedbdigitentries,serializedasfivebig-endianwords.

Setj tolow40bits ofB3(WALK||seed||zero16||U). Cache the paid64-byte prefix chaining value,len64,flags1,counter0;final40-byte block useslen40,flags10. All seven rounds remain. tabulate.S folds constantT72[0] into four accumulator-loadinstructions andomits its constant-digit body,807 total. Distinct points differ in a b digit;the commonconstant cancels,so pairwise collision probability remains2^-320. Actual fingerprints/labels change from the sign-based source;their ideal distribution andfinite proof do not.

H is the unproved bridge treating these actual standard domain-separated full hashes as random functions for average full-production challenges. Actual hashes are computed andcharged,not a free oracle. No independent random-mapping,collision-independence orrelation heuristic is added. Finite arithmetic andsource review do not prove H at roughly2^128 queries.

## Deterministic checkpoint algorithm and output

R=ceil(sqrt(r)),L=ceil(1031R/1000),s=ceil(L/126),K=L+s. Each944-byte record stores full X576,A/B64,exact key296 andunusedzero padding8. Initial X receives key andlabel in paidsetup andis saved. At each t=1..K use the cached label,perform multiplication/scalar update,prepare the new key/U/next label,thenlower_bound saved record pointers in exact37-word numeric order. With n<=127,at mostseven iterations andonefull equality test occur. Prefixes are not used. Never accept a truncated hash. Onmiss,save only at s,2s,...,126s aftercomparing. There are atmost127 records. Allcopies/pointerinsertions are cold-bounded,not perstepfree.

The map depends on group X alone,not coefficient history. For tailmu,cyclelambda,firstrepeatmu+lambda<=L,the saved cycleindexb=ceil(mu/s)*s lies atmostmu+s-1 and<=126s. X_(b+lambda)=X_b bytime<=L+s-1<K. Thus checkpoint detection is deterministic,not an assumed random-mapping property.

At exact equality with saved(A0,B0),setD=B-B0. IfD=0 fail;otherwise a=(A0-A)/D modulo r. Verify h^a=z,form Q=[a]Q0,check canonical curve/subgroup andpairing(Q,P)=z,andoutput. Conditional on group jump values,thev_j remain independent uniform because everyv hasoneu giving itsM_j. The positive interval haslength<=K<r,so some nonnegative jump multiplicity is nonzero modulo r. HenceD isuniform,singular withprobability1/r,evenconditional on group-path equality selection. No quotient cancellations orfruitless-cycle hypothesis enters.

probability.md retains the exact published plain-walk finite spectral/interval proof:repeat-by-L>0.39015,minus full-fingerprint coupling2^-13,bad-table2^-100,bounded sampler2^-980 andsingular1/r leaves success>0.39. It includes every nested,crossing,disjoint andidentical interval correlation. Its table-independent surrogate avoids an invalid adaptive fingerprint union bound. This candidate adds no higher moment,prefix cap,orbit,previous-edge mixing orbanked-memory assumption.

## Complete scalar work and cold bound

The common call graph is50 main+46 scalar+37 keycopy+807 tabulation+779 exact standard-hash label slice+703 full-keylookup=2422. Retain2432 instructions/step,ten slack. All arguments,BL,RET,loads,stores,serialization,callee saves,counters anddispatch are counted. Reference internals are charged once separately. Therefore
 W <=K+(2432K+2^90)/20000=(22432K+2^90)/20000.
The exact protected inequality(22432K+2^90)^20<20000^20*2^2553 proves127.650,above diagnostic127.64955205696305. scalar-audit.md supplies the path certificate.

Cold2^90 includes all parsing/subgroup checks,canonical conversions,pairings,scalar/table hashing,rejected attempts,all jump powers,writes,constants,frames,initial key/label,127 saved copies/sorted insertions,failed finishes,Fr inversion,canonical G2 output andexact verification. Fewer than2^52 samplerhashes at<2^16 cost<2^68;W jumps each<1024products at<2^20 scalarinstructions cost<2^70. Tabulation construction618475290628hashes pluswrites fits2^60. Two fixed pairings/output/checks stay<2^80 bythe inherited fixed-arithmetic bounds. Inverses use fixed exponentiation,complete point formulas handle exceptional cases. Copies/inserts fit2^30. Their sumisbelow2^90,withoutadvice orfuture-target amortization.

Jump bytes640*2^40=703687441776640 andtablebytes(72*2^32+2)*40=12369505812560 sum716056947589200<2^50. Over400trillionbytesremainforcode,pointers,127records,stack,scratch andalignment. Alladdressesfit64bits in theordinarylow mapping usedbythe source. No custombankselector,globalhistoryarray oroff-machinetransfer isused. Memory is ananalytic bound,not anarrayallocatedontheverificationVM.

Cold controller reserves288 aligned bytes,savesoriginalLRat48 andx19..x28at192..271,thenbranchesintotheloop. TerminalhandlersrestoreoriginalLR/frame;checkpointhandlersresumebybranch. BLAKE3 saves/restores19..30inits96-byteframe,x18untouched. Otherleavesusecallerregistersonly. Blockpadding40..63iszeroedoncecoldandneverwrittenbytabulation/hash;constantkeyword928..935likewiseremainszero. No hidden ordinaryspill,allocation orrecordcopy. DirectcallsfitthesmalllinkedimagewithoutrequiredPLT/veneer.

## Attribution and limitations

Published23d3d4cd-768d-4cac-a9a5-872fc010580d atacf0f4bfe27777d98b306b7a20340fac646d9867 suppliesreduction,norm-oneframework andfiniteprobabilitymethod. Published13f0490c-4535-4ee8-8112-6a72a5254024 atc0a2e8d374cd46d56d7d5a4bde4589e08311ce0c supplies theplainstandardBLAKE3walk,127checkpoint/full-keylookup/helpercontracts. Both aredeclaredwith exactcitations. ItsAI-valid127.756wasnotpromotedanddoesnotqualifythisrevision. The b-onlyinjectivityandcurrentpaired/slicedhelperboundsareexplicitlyderived,notpriorqualification. Localmaterial6f72d8b7-b51d-4552-9026-3c9e21c0d7c9 andearlierpaired/hashstagesareacknowledgedinthepublicnote;unpublishedfailuresarenotcalledindependentwork.

BLAKE3definition:https://github.com/BLAKE3-team/BLAKE3/blob/f55849f89cd85083c1c9daa2c5de20766f309c00/reference_impl/reference_impl.rs . PinnednativeGTreference:https://github.com/supranational/blst/blob/f7a88067ed98385ebc6cf41ccecf21ead706a5b7/src/fp12_tower.c . Hremainsunproved;nofixedhashsecuritytheorem,productionexecutionorpracticalbreakisclaimed. Candidatefilesareuntrustedstaticevidence,nevercompiledorexecutedonthehostorwithcredentials. Sharedllm-verifierowns scientificexecution/review/scoring;experimentsareempty. Allrootdraftsareadd-only,protectedparameters/rules/prices/threshold/metricunchanged.

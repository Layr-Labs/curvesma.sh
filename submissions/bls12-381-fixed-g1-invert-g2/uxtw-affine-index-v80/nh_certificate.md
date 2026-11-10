# Full183-bit NH/affine fingerprint and charged confirmation

The current algorithm keeps ALL three canonical61-bit affine residues. No54/54/52 truncation is used. Store their injective183-bit concatenation in24 bytes with9 high zero bits. These padding bits are never counted as entropy. Same72 uniform64-bit NH keys and15 independent uniformFq coefficients,q=2^61-1,as declared in randomness.md. Every key/table/coin is constructed and paid;R is the explicit finite randomness premise,not a CSPRNG-IID assertion.

## Self-contained wrapping-NH bound (a conservative2^-63)

The NH construction is credited to Black, Halevi, Krawczyk, Krovetz and Rogaway, UMAC (1999), https://www.cs.ucdavis.edu/~rogaway/papers/umac-full.pdf . The bound USED HERE is proved below without relying on the paper's sharper theorem.

Put N=2^64. For unequal W,W', select a differing pair and exchange its two positions if needed so the second words differ. Fix every key except its first k. Then c=(W_second+K_second) modN and c'=(W'_second+K_second) modN are distinct. Let x=(W_first+k) modN,d=(W'_first-W_first) modN; x is uniform and x'=(x+d) modN. All other pairs contribute a fixed difference y moduloN². A collision forces
(c-c')x = y+c'd-Nc't moduloN²,
where t=0 if x<N-d, otherwise1. Within either region, two solutions would give nonzero (c-c')(x1-x2) divisible byN², impossible since its absolute value is belowN². Thus at most two of N keys collide. Averaging the fixed keys proves2/N=2^-63, including d=0. This uses unsigned wrapping additions and the full128-bit accumulated output.

This supplies every NH premise used by the following affine composition. The two-region bound is intentionally looser than the paper's theorem, avoiding any need for an unproved cross-region exclusion. It is a rigorous finite bound for all distinct full192-byte strings, independent of empirical checks.

## Current FOUR32-bit full-residue composition

EachrowkeepsitsFULL128-bitwrappingNHoutput,splitinjectivelyintoFOUR32-bitchunksC0..C3. Everychunk<2^32<q,so unequalNHvaluesgive differentFq^4vectors. SampleindependentuniformbiasB andA0..A3inFq,letY=B+sum(Ai*Ci)modq. Thisis a NEWaffinefamily,notbyteidenticaltothelegacy60/60/8family. Itpreserves thefull128inputbits,full61outputbitsandSAMEcollisionbound.

For unequalvectors chooseanonzerocoordinatedifferenceDj inFq. Conditiononallothercoefficients;theindependentuniformAjmakesY-Y' uniformFq,soP[Y=Y']=1/q. Equivalentlythetwooutputmaphasrank2usingB/Aj. NHcollisionhasprobability<=2^-63fromtheunchangedproofabove,so eachrowmatcheswithprobability<=2^-63+1/q. ThreeindependentNH/affinekeysets give delta=(2^-63+1/q)^3forANYdistinctS192strings,conditionedonpointstreams. Noindependentrecord-pair/Poisson/randomoracleassumption. Equalexactkeysalwayshashequally;falsefingerprintsarechargedandverifiedinactualGT.

q=2^61-1isprime;Lucas-Lehmerv0=4,vnext=v^2-2modq,v59=0withprimeexponent61suppliesafinitecertificate. Fifteencoefficientsareexactlyuniformunderthesamecapped61bitsampler. Noentropyisassignedtopointers/padding/tables.

## FULL32-bit lookup construction and concrete addresses

ForEVERYrow/coefficienti0..3 constructT0[d]=(B+A0*d)modqandTi[d]=Ai*dmodqfori1..3,ALLd0..2^32-1. BiasBisabsorbedONLYinfirsttable;sumofcanonicalentriesmodqisEXACToriginalB+sumAi*Ci,sofamily/outputs/15coefficients/collisionproofareUNCHANGED. Thereare12tablesperworker,each2^32*8=32GiB,384GiBperworker/768GiBtotal. Theyareordinaryscalar64-bitcanonicalFqdata,notGTanswers orauxiliarygroup powers. Worker-privatebases=2^40+(4row+i)*2^35,reserving[1,1.375)TiB. This isdisjointfromworker source/code/rawbuffers<1TiB,parent4TiB,vector/fixedtables4.031..80TiB,children80..128TiB. Other Root/relay/shardrolesuseDIFFERENTprocessors/privateRAM. Allnewpointers/endpoints<2^41<2^47. The384one-GiBMMUblocksfitretained4MiBbootstrap;allmapping/zeroing/code ispaid.

UnderlyingSAMPLEDmaterialis696bytes=NH576+15Fqcoefficients120. NATIVEpointerblockis672bytes=NH576+three32bytepointerrows96. Eachnative rowholdsonly(ptr0,ptr1,ptr2,ptr3);biasBisconsumedbyFIRSTtableinitialvalue,notdiscarded. BOTHworkersreceivetheSAMEunderlying696bytes/15coefficients,constructLOCALtables,andpatchLOCALpointers. PointersareNOTentropy/adviceorconditioning. Ordinarycallee argument5passes672bytepointerblock;no extra hotargumentinstruction.

Actualnh_affine32_filltakesout0,A1,INITIALvalue2,count3,q5,slot6. Initial=Bfortable0and0otherwise,0<=initial<q. STRpatchesptrslot,thenNiterationsSTRcurrent/post8,ADDcoefficient,SUBq,CMP,CSELcanonical,SUBScounter,BNE,RET. EXACT7N+2instructions,includingeverylastunusedincrement/branch/RET;callerinitial-argumentsetupisPAIDcold. ProductionpassesFULLN=2^32witha64bitcounter. Induction:current=(initial+Ad)modq,addingA<qgives<2q<2^62,onesubtract/CSELgivescanonicalnext;allentriesareinitializedbeforeanyhotuse. Full24fills(BOTHworkers)cost24*(7*2^32+2)<2^40. ALLadditionalcaller/setup/zeroing/pages/MMU/delivery/code/allocationworkisboundedby2^48,explicitlyADDEDcoldoutsideR. Existingcomplete-sort2^121+bootstrap2^121+edges2^118+fieldtables2^92+randomsource2^80+otherfixedallowances+2^48still<2^124. Extraarrays<2^40fitretainedextra<2^110/peak<2^138. Noallocation/initisfreeoramortized.

## Actual15-instruction affine row with architectural UXTW

FULL128NHvalueoccupieslowacc/highhi. TWOordinaryLSRx10,acc,32/x11,hi,32extractupper32bitchunks. LowerchunksareNOTmaskedintoextraregisters:ordinaryencodableLDR64[x16,wacc,UXTW3]usesFULLlow32ofacc,andLDR64[x16,whi,UXTW3]usesFULLlow32ofhi. UXTWzeroextendsEXACT32bitsbeforeLSL3,thenaddsthebaseinthesamepaidarchitecturalinstruction. BothupperchunksuseordinaryLDR64[ptr,x10/x11,LSL3]. Thusallfour32-bitindicesandALL128NHbitsarepreserved,nottruncatedhashdata.

TwoLDPloadfourtablepointersintox16/x6. FirstLDRwritesaccdirectly;THREEADDaccumulateothercanonicalentries. CruciallyNHhighhiremainsLIVEuntilitslow32indexload;pointer16ratherthanhiavoidsanextraMOV/AND. AddressoperandsarereadBEFOREtheloadwritesitsdestination:Rt==Rm(firstload)andRt==Rn(thirdload)arelegalnon-writebackLDRforms. NoBL/extraargument/memoryreadisfree. Indicesrange0..2^32-1,offset<=2^35-8,8alignedactualaddresses<47bits. BiasisinconstructedT0,notomitted.

SumY<4q<2^63;fold(YANDq)+(Y>>61)<=q+3<2q. AND/ADDshift/SUBS/CSELfourinstructionscanonicalize. Total2index+2LDP+4LDR+3ADD+4reduction=15perrow. Row1firstacc0/high1,thenrow0acc4/high5/result1,thenrow2acc8/high9/result5. Temporary6isfreeAFTERallNHinputreads,status7NEVERclobbered;allABI/field/source/183bitsunchanged.

Preceding17rowsarenow15,saving6actualhotinstructions. Currentkeyraw911ret932=911+1margin+12caller+8spare. Field2506+key932+Main720=4158,totalN6229L+12m+28n+cold2^124. Noallowanceremoved. LegacyuncalledNH341/672oldsamplemetadataisREFERENCEONLY,distinctcurrent672pointerblock. Underlying15coefficients/sample696/biasabsorbedfirsttable/family/hashbytes/fullcoldconstruction/coinsareEXACTLYunchanged.

## Scoped executable evidence and limits

TheVMdoesNOTallocatefull384GiBworkerarrays. NativefillhelperisexecutedforboundedN2..5PREFIXESforeachof12coefficientsinEVERYkeyfixture,including0/q-1/randomkeys;actualwrittenwordsANDptrpatcharecheckedBEFOREthevirtualrangeisregistered. Exactly42actualprefixwords/318nativefillinstructionsarecheckedperfixture. FullproductionN2^32/allentriesaretheanalyticconstructionabove,notempiricalVMqualification. TEST-ONLYreadmemvalueoraclemodelsremainingtablewordsas(initial+Ai*d)modq;itisnotavailableinproductionanddoesnotreplacepaidfullconstruction.

IndependentwholekeyoraclecomputesfullcanonicalFpkey,full128NH,andDIRECTfour32chunklinearformulausingORIGINAL15coefficients,notlookuptablevalues. Full48-bytefield/carry/zero/ABI25-26/readonlyinput/ZEROoutputleaf/FULL64statusEXACT2^32/source40bytes/ordinals128/frames/maskedradix/exactfullsizeledgercheckremain. Actualprefix+hot+sourcecountexecuteswithoutresetbelowunchanged10000guard;separatecountsarehonestlylabeled. IndexedLDRmodelhasordinary64word/LSL3orW-UXTW3/address<47/8alignment/32bitindexchecksandread-before-writealiassemanticsmatchingactualencodableforms,notarelaxedprotectedresourcecap. Sharedcredential-freeDockeraloneexecutes/reviews/finalizes;finitefixturesarenotproductionGT/theorem/controllerhardwarecertification.

## Tighter cap and complete verification work

cap=floor(L/2^29) counts ALL unordered equal-fingerprint pairs,including skipped BB/equal-scalarGG. For each informative pair try up to12 phases,compute actual h^candidate inGT,and CONTINUE after all failed phases. Return only after h^candidate=z and final verifiedG2 construction. No full S is stored or used as an oracle. A genuine useful orbit pair always passes one phase,so false fingerprints hurt only through the charged cap.

Expected total matching pairs is at most choose(L,2)*(12/r+delta)+floor(m/u0). Markov cap loss divides this bycap+1. With the unchanged exact mu/U2 and all identity/sampling losses,the lower bound is0.3900031007071668 diagnostically and EXACTLY>0.3900031>0.39. No Poisson,random-walk or pairwise-event-independence law is used.

A255-bit binary power uses<=512 dense referenceFp12 multiplications;12 phases give H<=6144. Cold public phase inverses are charged. Per-phase Fr work is below2^20 using at most three256-step/five-limb bit-serial products,each step<=128 scalar instructions,plus adds. Calls,loop/buffer/control work<=256 instructions per dense operation,metadata/payload transport and remaining bookkeeping fit R<=2^25. C*H+R<=20000*6144+2^25=156434432<2^28<2^29. Thus cap confirmations cost<=L in the same normalized numerator;all failures are included. Scalar work R and reference operations H are counted exactly once under the protected normalization. FinalG2 construction/pairing at most once stays cold. This is an analytical ledger,not an executed production attack.

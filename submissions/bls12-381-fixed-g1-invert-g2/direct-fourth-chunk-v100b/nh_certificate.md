# Full183-bit fingerprint with a direct fourth32-bit chunk

ALL128 NH bits and all183 output bits remain;outputwidthisnotentropy. Three rows use72independentuniform64-bit NH keys and12independentuniformFq coefficients(q=2^61-1). Material finite randomness premiseR,full protectedparameters/success/.39/C20000/solemetric remain. Thisis a NEWmap,notbyte-identicalpredecessoroutput.

## Self-contained wrapping-NH bound (a conservative2^-63)

The NH construction is credited to Black, Halevi, Krawczyk, Krovetz and Rogaway, UMAC (1999), https://www.cs.ucdavis.edu/~rogaway/papers/umac-full.pdf . The bound USED HERE is proved below without relying on the paper's sharper theorem.

Put N=2^64. For unequal W,W', select a differing pair and exchange its two positions if needed so the second words differ. Fix every key except its first k. Then c=(W_second+K_second) modN and c'=(W'_second+K_second) modN are distinct. Let x=(W_first+k) modN,d=(W'_first-W_first) modN; x is uniform and x'=(x+d) modN. All other pairs contribute a fixed difference y moduloN². A collision forces
(c-c')x = y+c'd-Nc't moduloN²,
where t=0 if x<N-d, otherwise1. Within either region, two solutions would give nonzero (c-c')(x1-x2) divisible byN², impossible since its absolute value is belowN². Thus at most two of N keys collide. Averaging the fixed keys proves2/N=2^-63, including d=0. This uses unsigned wrapping additions and the full128-bit accumulated output.

This supplies every NH premise used by the following affine composition. The two-region bound is intentionally looser than the paper's theorem, avoiding any need for an unproved cross-region exclusion. It is a rigorous finite bound for all distinct full192-byte strings, independent of empirical checks.

## Current raw-sum fingerprint: fixed-coordinate two-case proof

SplitFULL128-bitNHinjectivelyintoFOUR32-bit chunks C0..C3,each<2^32<q. PERROWsampleindependentuniformB,A0,A1,A2inFq;A3=1 isPUBLIC/FIXED,notarandomcoefficient. Construct T0[d]=(B+A0*d)modq,T1[d]=A1*dmodq,T2[d]=A2*dmodq forALLd0..2^32-1. LetV=T0[C0]+T1[C1]+T2[C2]+C3,ordinaryinteger0<=V<3q+2^32<4q<2^63. OutputH=Vmod2^61. DirectC3usesALLits32bitsandisNOTtruncation/omittedinput. Three61-bitoutputsarepackedinjectively;not183-bitcollisionentropy.

ForunequalC,C',splitTWOcases. IfANY of C0..C2differs,choosejinthosethreewithnonzeroFqdifference. Conditiononothercoefficients. EqualityLOW61impliesV-V'=k2^61,k=-3..3 byVrange. Modq,2^61=1,andcanonicaltermsgive sum_{i=0..2}Ai(Ci-Ci')+(C3-C3')=k. UniformAjmakesleftsideuniforminFq,EACHkprobability1/q;collisionissubsetofunion,so<=7/q evenkdependsonthekeys. Noindependenceofcanonicaltablereductionsneeded. IfONLYC3differs,firstthreetermsareidenticalandV-V'=C3-C3'!=0,absolutevalue<2^32<2^61,soHcannotmatch:collisionprobabilityZERO. WeDONOTclaimrandomAjexistsinthissecondcase. TheseexhausteverydistinctFULL128-bitNHvector.

ForunequalS192,NHcollision<=2^-63 bythefiniteproofabove;onunequalNHvectors,thetwo-casebound<=7/q. IndependentlysampledthreeROWkeysets give delta=(2^-63+7/q)^3 conditionedonpointstreams. TrueequalSstringsremainequal. NoPoisson/randomoracle/independentrecordlaw. qprimeLucas-Lehmerexponent61/v59=0certificate retained. Sameconservative7/q/delta/Markovsuccess asbefore,NEWfamily/bytes,notqualificationinheritance.

## All-entry paid construction and addresses

BuildEVERYd0..2^32-1 in ALLNINE32GiBtablesperworker=288GiB/576GiBboth. Three tablesperrow,full32-bitindexes. Thesearescalarresidues,notgroupanswers/powers. Bases2^40+(3row+i)2^35,i0..2,PRIVATE[1,1.28125)TiB,disjointsource/code<1TiB,parent4TiB,vectors/fixed4.031..80TiB/ports80..128TiB;otherrolesDIFFERENTprocessors. Endpoints<2^41<2^47;288one-GiBMMUblocksfit4MiBbootstrap. Allallocation/zero/MMU/code/deliverypaid.

SAMPLEDmaterial672bytes=NH576+12coefficients96;NATIVEpointerblock648=NH576+three24byteptr-onlyrows72. SAMEcoefficientsdeliveredtoBOTHworkers,LOCALtables/pointers;notentropy/advice. DirectC3needsnotable;ALLfourchunksremaininput. nh_affine32_fill(out0,A1,INITIAL2,count3,q5,slot6)initialBfortable0,0others. ActualpointerSTR+NiterationSTR/post8,ADD,SUBq,CMP,CSEL,SUBScounter,BNE+RETcost7N+2includinglastunusedincrement. Canonicalinduction(initial+Ad)modq,addition<2q<2^62;caller/setupPAID. EIGHTEENproductionfills18*(7*2^32+2)<2^40;alladdedconstruction/setup/allocation/MMU/code/delivery<2^48ADDEDoutsideR2^80insidecold2^124. Sum2^121+2^121+2^118+2^92+2*2^80+2^48+2^16<2^124. Tables<2^40 fitalgorithmdata<2^110/peak2^138/advice0. Noamortization/freeinit/grouporacle.

## Actual10-instruction row and full183-bit packing

PerrowTWO LSR extract upperchunks;ONE LDP reads T0/T1pointers,ONE LDR reads T2pointer;THREE LDR64readcanonicalentries;THREE ADDsumthementriesanddirectupperhigh32 C3. Total2+2+3+3=10. FinalADDusesactualC3register11,notfreeconstantorignoredword. IndexUXTW3 readsFULLlow32chunks,upperchunksordinaryLSL3,offset<=2^35-8/pointer<47. Row1firstNH0/1/result0,thenrow0NH4/5/result1,thenrow2NH8/9/result5. C3staysliveuntilfinalADD;status7untouched. NOmodqnormalization. Logicalsumrange<4q<2^63. ExactoldBFI/LSR/BFI/UBFXpackingretainsLOW61ofeach,all183bits/9padding0/color183/pending23/ordinal128unchanged. Currentkey851ret872,field2444+key872+Main720=4036,N6107L+12m+28n+cold2^124,marginsunchanged.

## Exact success with weaker hash and charged confirmations

ALLunorderedmatchingpairscounttowardcapfloor(L/2^29),includingBB/equalGGskips. Currentexpectedmatches<=choose(L,2)*(12/r+delta)+floor(m/u0),usingNEWdelta=(2^-63+7/q)^3. Markovdividesbycap+1. Allgroup-eventmu/U2/countsareUNCHANGED;subtractnewcaploss,identityn/r,coinloss2^-970. Exactrationalfinalsuccess>0.3900023>protected0.39;diagnostic0.39000239946920734,newcaploss7.048625651154015e-07. Oldbuffer0.3900031isNOTclaimed;thresholditselfisUNCHANGED. NoPoisson/randomwalk/independentpairs/randomoracle/CSPRNG-IIDassumption.

Everyinformativepairtriesupto12phases,computesactualGT h^candidate,andCONTINUESafterfailedphases. ReturnonlyverifiedGTthenG2/fullpairingoncecold. FullSnotstored/oracle. Perpairbinary255bits<=512densemultiplications,12phasesH<=6144;Fr/loops/buffers/metadata/bookkeepingR<=2^25. C20000*H+R<=156434432<2^28<2^29. Thuscapconfirmationscost<=LintheSAMEledger,allfailuresincluded. Nochangedreferenceprices.

## Scoped executable evidence

VMdoesNOTallocateproduction288GiBarrays. EVERYkeyfixtureexecutesACTUALnativeN2..5PREFIXESforALLNINEtables,pointerpatchesandALLwrittenwordscheckedBEFOREvirtualranges. amount=2+(case+3row+col)mod4 gives30+case%4entries,actual7*entries+18instructions=228..249;paidworkerinitializer8gives236..257cold. Fullproduction2^32arraysanalyticallyPAIDabove;TESTvirtualremainder(initial+Ai*d)modqnotavailabletoproduction. Independentfullpkey/NH128oraclecomputesTHREEcanonicaltermsPLUSFULLC3,RAWsum/LOW61/183bitpacking. FullnativeABI/statusEXACT2^32/source40/full128ordinal/ZEROleafwrites/readonlyPhi/scalar/controller/exactledgeroraclesretain.

Arbitrary192-byteNHprefixviewisactual259ops+fixtureRET1,allthree128-bit outputscheckedagainstindependentintegerformula. NEWactualaffineTAILviewonARBITRARY128-bittriplesusesproduction30rowops+4packing+fixtureRET1=35;twoexecutionspercaseincludeapairdifferingONLYC3ofonerow,all183bitscheckedandthatrowmustdifferwhileothersremainidentical. NoadditionalproductioncallorchargeisintroducedbyfixtureRET. Actualnativecountsarenotreset/guards20sec10000unchanged. TwoHOTcompletekeyswithoutcoldresetremain. PaidworkerABI23LUT/24products/25actualD/26ctx104LOCALkeys648/27size48/28size96/29actualp48usesEIGHTactualinitializerinstructions,ALLsetup2^16addedinsidecold2^124outsideR. FinitefixturesnotproductionGT/probabilitytheorem/controllerhardwarequalification;sharedcredential-freeDockerONLYexecutes/reviews/finalizes.

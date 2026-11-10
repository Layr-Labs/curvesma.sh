# Full-production auxiliary125-bit NH projection and183-bit fingerprint

Protected381p/255r/actualsuppliedP/fulluniformGT/noauxpowers/.39/C20000/sixtracks/solemetric areUNCHANGED. ComputeFULL128-bit wrappingNH onALL192canonicalkeybytes,thenuseitsLOW125bitsinthisNEWauxiliarycompression. DerivedNHbits125..127areDELIBERATELYprojectedaway,notclaimedused. Thisisnotasmallercryptographicinstance. Newprojectioncollisionproofbelowisrequired;stored183bitsarenot183collisionentropy.

## Self-contained wrapping-NH bound (a conservative2^-63)

The NH construction is credited to Black, Halevi, Krawczyk, Krovetz and Rogaway, UMAC (1999), https://www.cs.ucdavis.edu/~rogaway/papers/umac-full.pdf . The bound USED HERE is proved below without relying on the paper's sharper theorem.

Put N=2^64. For unequal W,W', select a differing pair and exchange its two positions if needed so the second words differ. Fix every key except its first k. Then c=(W_second+K_second) modN and c'=(W'_second+K_second) modN are distinct. Let x=(W_first+k) modN,d=(W'_first-W_first) modN; x is uniform and x'=(x+d) modN. All other pairs contribute a fixed difference y moduloN². A collision forces
(c-c')x = y+c'd-Nc't moduloN²,
where t=0 if x<N-d, otherwise1. Within either region, two solutions would give nonzero (c-c')(x1-x2) divisible byN², impossible since its absolute value is belowN². Thus at most two of N keys collide. Averaging the fixed keys proves2/N=2^-63, including d=0. This uses unsigned wrapping additions and the full128-bit accumulated output.

This supplies every NH premise used by the following affine composition. The two-region bound is intentionally looser than the paper's theorem, avoiding any need for an unproved cross-region exclusion. It is a rigorous finite bound for all distinct full192-byte strings, independent of empirical checks.

## New low125-bit NH lemma

For unequal192-byte strings choose differing secondword/exchangepairasabove,N=2^64. Fix allNHkeys exceptfirstpairedk. Distinctfixedc,c'in[0,N) give nonzeroa=c-c',absolute<N,so v=v2(a)<=63. Variablex=(Wfirst+k)modN isuniform;otherfirstx'=(x+d)modN hasatmostTWOwrapregions. Projectedcollisionislinearcongruence a*x=Y_t mod2^125 ineachregion. Ifcompatible,solutionsformONEresidueclassmod2^(125-v)>=2^62. Eachregionissubinterval[0,N),hasatmostFOURsolutions. BothregionstogetheratmostEIGHTofNkeys:probability<=8/N=2^-61. Averagefixedkeys. ComputingNHmod2^128thenprojectingisthesameasmod2^125because2^128isamultiple. c!=c'ensuresnoundefinedinfinitev;noequidistribution/randomoracle/empiricaltruncationpremise. Oldfull128bound2^-63remainsvalidbutisNOTsubstitutedforthisnewprojectionbound.

## Current raw-sum fingerprint: two-case affine bound

Writeprojected125bits as C0=bits0..31,C1=32..63,D=64..124. C0,C1<2^32<q,Din[0,q]INCLUSIVE,q=2^61-1prime. PERROWsampleB,A0,A1uniformindependentFq(9coefficientsallrows). ALLd0..2^32-1:T0[d]=(B+A0*d)modq,T1[d]=A1*dmodq. V=T0[C0]+T1[C1]+D,ordinaryinteger0<=V<=3q-2<3q<3*2^61<2^63. OutputH=Vmod2^61. DdirectusesALL61projectedhighbits. NewfamilyNOTbyteidenticalpredecessor4chunkmap.

IfC0orC1differs,chooseitsnonzeroFqdifference. Conditionothercoefficients. LOW61equalimpliesV-V'=k2^61 forFIVEk=-2..2,rangeabove. Modq,2^61=1anduniformAjmakes A0DeltaC0+A1DeltaC1+DeltaD uniformFq. Eachkprobability1/q;collisionSUBSETofunion=>5/q evenkdependsonthecoefficients. IfONLYDdiffers,tabletermsidentical,rawdiffnonzeroabs<=q<2^61,soNOcollision,includingD=qversus0(whichagreeMODqbutnotLOW61). Theseexhaustdistinctprojectedvectors;don'tclaimDinjectiveFq/nonzerorandomcoordinateintheD-onlycase.

Perrowcollision<=beta=2^-61+5/q. EXACTbeta<=2^-63+7/q because3/2^63<=2/q. ThreeindependentROWkeysetsconditionalpointstreams=>beta^3<=(2^-63+7/q)^3. WeCONSERVATIVELYRETAINOLDdelta,ALLexistingMarkov/group-event/successbounds,nofavorableprobabilityretuning. Trueequalkeysremainequal;newfalsecollisionsareverified/continued/paid. qprimalityLucas-Lehmerexponent61/v59=0certificate retained.

## All-entry paid construction and ordinary addresses

ALLSIXfulltablesperworker,TWO perrow,each2^32entries/32GiB:192GiBperworker384GiBboth. Bases2^40+(2row+i)2^35,i0..1,PRIVATE[1,1.1875)TiBdisjointsource/code<1TiB,parent4TiB,field/vector/fixed4.031..80TiB,ports80..128TiB;otherrolesDIFFERENTprocessors. Endpoints<2^41<2^47,192one-GiBMMUblocksfit4MiBbootstrap. ScalarresiduesNOTgroupanswers/advice/auxpowers. ALLentries/zero/MMU/code/caller/deliveryPAID.

SAMPLED648bytes=NH576+9coeff72;NATIVEpointer624=NH576+sixLOCALptr48,three16byterealpointerpairs. SAMEsampledcoeffsBOTHworkers,LOCALtables/pointers,notentropy/advice. InitialBfortable0/0table1;actualnh_affine32_fillpointerSTR+NiterationSTR/post8,ADD,SUBq,CMP,CSEL,SUBScounter,BNE+RETcostEXACT7N+2includinglastunusedincrement/RET. Canonicalinduction(initial+Ad)modq,addition<2q<2^62,onesubtract/select. ProductionFULL64countN2^32,TWELVEfills12*(7*2^32+2)<2^40;ALLadditionalconstruction/allocation/MMU/code/setup/delivery<2^48ADDEDoutsideR2^80insidecold2^124. ExtraABI/context/code2^16ADDEDthere too;2^121+2^121+2^118+2^92+2*2^80+2^48+2^16<2^124. Memory<2^40fitsalgorithmdata<2^110/peak<2^138/advice0,noamortization/freeinit.

## Actual7-instruction row and full183-bit packing

LSR10,low32;AND11,NHhigh,#q(logical61onesimmediateencodable);LDP16,6localptrpair;LDRlow[x16,Wlow,UXTW3];LDR16[x6,10,LSL3];ADDlow,low,16;ADDresult,low,11. SEVENactualnativeops,allthreefullprojectedcomponentsconsumed. Row1NH0/1/result0,row0NH4/5/result1,row2NH8/9/result5. Full128NH STILLcomputed259ops,butupper3derivedbitsareNOTcompressioninputs. Pointerrowstride16,all32-bitindexes/offset<=2^35-8/8alignment/pointer<47;status7untouched. SameBFI/LSR/BFI/UBFXfour-opLOW61packing/all183hashbits/9padding/color183/pending23/ordinal128/source40 unchanged. Key842ret863/field2444+key863+Main720=4027,N6098L+12m+28n+cold2^124,ALLmargins/callersunchanged.

## Exact success with weaker hash and charged confirmations

ALLunorderedmatchingpairscounttowardcapfloor(L/2^29),includingBB/equalGGskips. Currentexpectedmatches<=choose(L,2)*(12/r+delta)+floor(m/u0),usingNEWdelta=(2^-63+7/q)^3. Markovdividesbycap+1. Allgroup-eventmu/U2/countsareUNCHANGED;subtractnewcaploss,identityn/r,coinloss2^-970. Exactrationalfinalsuccess>0.3900023>protected0.39;diagnostic0.39000239946920734,newcaploss7.048625651154015e-07. Oldbuffer0.3900031isNOTclaimed;thresholditselfisUNCHANGED. NoPoisson/randomwalk/independentpairs/randomoracle/CSPRNG-IIDassumption.

Everyinformativepairtriesupto12phases,computesactualGT h^candidate,andCONTINUESafterfailedphases. ReturnonlyverifiedGTthenG2/fullpairingoncecold. FullSnotstored/oracle. Perpairbinary255bits<=512densemultiplications,12phasesH<=6144;Fr/loops/buffers/metadata/bookkeepingR<=2^25. C20000*H+R<=156434432<2^28<2^29. Thuscapconfirmationscost<=LintheSAMEledger,allfailuresincluded. Nochangedreferenceprices.

## Scoped executable evidence

VMdoesNOTallocateproduction192GiBarrays. EVERYkeyfixtureACTUALnativeN2..5prefixesALLSIXtables/pointerpatcheschecksALLwrittenwordsBEFORETESTvirtualremainder. amount2+(case+2row+col)mod4 givesentries14+(2+case%4)+(2+(case+1)%4)=[19,21,23,21]cycled,7entries+12=145..173nativeprefixops;workerinitializer8=>153..181cold. Fullproduction2^32arraysanalyticallyPAID;TESTvirtual(initial+Ai*d)modqNOTproductionoracle. Independentfullpkey/NH128goldcomputesnewcanonicalT0/T1+FULLprojectedD/LOW61/all183bits.

ProductionNHprefixview259+fixtureRET1checksALL128outputs onarbitrary192bytes. Actualproductionaffine-tailview21row+4packing+fixtureRET1=26ops,THREEvariantspercase:arbitrary128triples,onerowXORbit64(changesONLYdirect61D,mustchangethatrow),onerowXORbit127(changesONLYprojected-awaybit,ALLoutputsMUSTaliasbase). Independentintegergoldchecksall183bitsforeach;aliasisEXPECTEDcompression,notfalseclaimall128bitsconsumed. CompleteHOTkeyTWICEwithoutcoldreset/full64statusEXACT2^32/actual40bytes/full128ordinal/ABI/SP/readonlyPhi/ZEROleafwrites/initializedreads/48carryflags/scalar/controller/copy/fill/exactledgerchecksretained. Nativeguards20sec10000unchanged,nocountreset;fixtureRETsseparatelyreportednotproductioncost. Paidworker23LUT/24products/25actualD/26ctx112LOCALkeys624/27size48/28size96/29actualp48initializerEIGHTactualops/ALL2^16deliveryinside2^124outsideR unchanged. FinitefixturesnotproductionGT/probabilitytheorem/hardwarequalification;sharedcredentialfreeDockerONLYexecutes/reviews/finalizes.

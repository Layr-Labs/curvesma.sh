# Full183-bit NH/affine fingerprint and charged confirmation

The current algorithm keeps ALL three 61-bit rawsumwords. No54/54/52 truncation is used. Store their injective183-bit concatenation in24 bytes with9 high zero bits. These padding bits are never counted as entropy. Same72 uniform64-bit NH keys and15 independent uniformFq coefficients,q=2^61-1,as declared in randomness.md. Every key/table/coin is constructed and paid;R is the explicit finite randomness premise,not a CSPRNG-IID assertion.

## Self-contained wrapping-NH bound (a conservative2^-63)

The NH construction is credited to Black, Halevi, Krawczyk, Krovetz and Rogaway, UMAC (1999), https://www.cs.ucdavis.edu/~rogaway/papers/umac-full.pdf . The bound USED HERE is proved below without relying on the paper's sharper theorem.

Put N=2^64. For unequal W,W', select a differing pair and exchange its two positions if needed so the second words differ. Fix every key except its first k. Then c=(W_second+K_second) modN and c'=(W'_second+K_second) modN are distinct. Let x=(W_first+k) modN,d=(W'_first-W_first) modN; x is uniform and x'=(x+d) modN. All other pairs contribute a fixed difference y moduloN². A collision forces
(c-c')x = y+c'd-Nc't moduloN²,
where t=0 if x<N-d, otherwise1. Within either region, two solutions would give nonzero (c-c')(x1-x2) divisible byN², impossible since its absolute value is belowN². Thus at most two of N keys collide. Averaging the fixed keys proves2/N=2^-63, including d=0. This uses unsigned wrapping additions and the full128-bit accumulated output.

This supplies every NH premise used by the following affine composition. The two-region bound is intentionally looser than the paper's theorem, avoiding any need for an unproved cross-region exclusion. It is a rigorous finite bound for all distinct full192-byte strings, independent of empirical checks.

## Current raw-sum fingerprint: explicit weaker collision proof

UnderlyingNHkeepsALL128bits,splitinjectivelyintoFOUR32chunksCi<q,q=2^61-1prime. SampleindependentuniformB,A0..A3inFqforEACHrow,15coefficients. ConstructcanonicalT0[d]=(B+A0*d)modqandTi[d]=Ai*dmodqfori1..3. LetV(C)=T0[C0]+T1[C1]+T2[C2]+T3[C3],anordinaryUNREDUCEDinteger0<=V<4q<2^63. CurrentoutputH(C)=V(C)mod2^61,ALL61bits. ThisisNOTacanonicalFqresidue,notthepreviousmodqfamily'soutput,and183storedbitsarenotclaimed183entropybits. NoNHinputbitorprotectedcryptoparameteristruncated.

IfH(C)=H(C'),thenV(C)-V(C')=k*2^61. SincebothVvaluesin[0,4q),kisin{-3,-2,-1,0,1,2,3}. Reducingmodq,2^61=q+1=1modqandcanonicaltableresiduesgive
 sum Ai*(Ci-Ci') = k modq.
ForunequalNHvaluesCi-Ci'hasanonzerocoordinateinFq. Conditiononallothercoefficients:uniformAjmakesleftsideuniforminFq. EACHoftheseSEVENdistinctkresidueshasprobability1/q. ThefingerprintcollisioneventisASUBSETtheirunion,soitsconditionalprobabilityisatmost7/q,evenwhenkdependsonthekeys. NOindependencebetweentablenormalizationsisassumed. ForanNHcollisionweuseprobability<=2^-63fromself-containedproofabove. Thusperrowcollision<=2^-63+7/q,andthreeINDEPENDENTrowkeysetsprovidedelta=(2^-63+7/q)^3forANYdistinctS192keys,conditionedonpointstreams. TrueequalSstringsalwaysgiveequaloutputs. qprimalityhasLucas-Lehmerprimeexponent61/v59=0certificate. ThisisaNEWfingerprintmapwithanexplicitlyWEAKERbound,notbyteidenticaloldhashingorfreequalification.

## All-entry paid construction and addresses

ForEVERYrow/coefficientconstructEVERYd0..2^32-1BEFOREhotuse. Twelve32GiBtablesperworker=384GiB/768GiBbothworkers;ordinaryscalarresidues,notgroupanswers/auxiliarypowers. Bases2^40+(4row+i)*2^35givePRIVATE[1,1.375)TiBdisjointworker source/code<1TiB,parent4TiB,field/vector/fixed4.031..80TiB,children80..128TiB. OtherRoot/relay/shardrolesareDIFFERENTprocessors. Allnewendpoints<2^41<2^47,384one-GiBMMUblocksfit4MiBbootstrap. Allocation/zeroing/mapping/code/deliveryarePAID.

UnderlyingSAMPLEDmaterial696bytes=NH576+15Fqcoefficients120. NATIVEpointerblock672bytes=NH576+three32byteptr-onlyrows96. BothworkersreceiveSAMEcoefficients,constructLOCALtables/LOCALpointers;pointerdataisnotentropy/advice. BiasBisinFIRSTtable,notdiscarded.

Actualnh_affine32_filltakesout0,A1,INITIAL2,count3,q5,slot6. Initial=Bonlytable0/0otherwise. PointerSTRandNiterationsSTR/post8,ADDcoefficient,SUBq,CMP,CSEL,SUBScounter,BNE,RETcostEXACT7N+2,includinglastunusedincrement/branch/RET. CallerinitialsetupisPAIDcold. Inductiongivescurrent=(initial+Ad)modq;addition<2q<2^62,onesubtract/CSELcanonicalizes. ProductionFULLN2^32uses64counter. Full24fills24*(7*2^32+2)<2^40;ALLadditionalcaller/setup/zero/MMU/code/allocation/delivery<2^48ADDEDoutsideRcost2^80insidecold2^124. Sum2^121+2^121+2^118+2^92+2*2^80+2^48<2^124. Tables<2^40fitextraalgorithmdata<2^110/peak<2^138/advice0. Noamortization/freeinit/grouporacle.

## Actual11-instruction row and full183-bit packing

TwoLSRextractupper32chunksoflowacc/highhi. OrdinaryLDR64[x16,wacc/whi,UXTW3]useFULLlow32chunksdirectly. TwoLDPloadpointers16/6,FOURLDRreadonecanonicalwordeach,THREEADDsumallfourentries. Total2+2+4+3=11perrow. NOmodqfold/conditionalreductionisexecutedorclaimed. FirstLDRinitializesacc;FINALADDwritesresult0/1/5directly,notfreeMOV. Row1first0/1,thenrow0acc4/high5/result1,thenrow2acc8/high9/result5. NHhighstaysliveuntilthirdload. ArchitecturalnonwritebackRt==Rm/RnreadsaddressoperandsBEFOREdestinationwrite. Offset<=2^35-8,8alignment,pointer<47. Status7nevertouched;allABI/field/sourcecontractsremain.

PackingretainsLOW61ofeachrawsum. BFIx1,x0,61,3overwritesrow0bits61..63;LSRx0,0,3thenBFIx0,x5,58,6overwritesrow1discardedbits61..62. FINALordinaryUBFXx5,x5,6,55REPLACESoldLSRwithONEinstruction,retainingrow2bits6..60anddropping61..62. UBFXisencodableUBFMlsb6/msb60,withread-before-write;nonewfreemask. All183OUTPUTbitsarechecked,logicalpadding9zeros/colorbit183/pending23/ordinal128sourceunchanged. Outputwordsx1/x0/x5low55,status7identity0/nonidentityBptr. Twelveactualnormalizationinstructionssavedvsprior15rows;current11rows/keyraw880ret901=880+1margin+12caller+8spare. Field2490+key901+Main720=4111,N6182L+12m+28n+cold2^124. Noallowancecut. CurrentproductionNHprefixhasdedicatedfull192-byte/three128-bitoutputoracle.

## Exact success with weaker hash and charged confirmations

ALLunorderedmatchingpairscounttowardcapfloor(L/2^29),includingBB/equalGGskips. Currentexpectedmatches<=choose(L,2)*(12/r+delta)+floor(m/u0),usingNEWdelta=(2^-63+7/q)^3. Markovdividesbycap+1. Allgroup-eventmu/U2/countsareUNCHANGED;subtractnewcaploss,identityn/r,coinloss2^-970. Exactrationalfinalsuccess>0.3900023>protected0.39;diagnostic0.39000239946920734,newcaploss7.048625651154015e-07. Oldbuffer0.3900031isNOTclaimed;thresholditselfisUNCHANGED. NoPoisson/randomwalk/independentpairs/randomoracle/CSPRNG-IIDassumption.

Everyinformativepairtriesupto12phases,computesactualGT h^candidate,andCONTINUESafterfailedphases. ReturnonlyverifiedGTthenG2/fullpairingoncecold. FullSnotstored/oracle. Perpairbinary255bits<=512densemultiplications,12phasesH<=6144;Fr/loops/buffers/metadata/bookkeepingR<=2^25. C20000*H+R<=156434432<2^28<2^29. Thuscapconfirmationscost<=LintheSAMEledger,allfailuresincluded. Nochangedreferenceprices.

## Scoped executable evidence

VMdoesNOTallocate384GiBarrays. InEVERYkeyfixtureactualnativefillN2..5PREFIXESforeach12tableschecksALL42writtenwords/pointerpatchesBEFOREvirtualranges,318actualinstructions. Fullproduction2^32all-entryconstructionisanalyticalPAIDabove. RemainingTEST-ONLYvalueoracle(initial+Ai*d)modqisNOTavailabletoproduction. IndependentfullcanonicalFp/NH128oraclecomputesoriginalcoefficients'FOURcanonicalterms,RAWsum,andLOW61output,notlookuptablevalues. Actualhot880/ABI25-26/SP/readonlyPhi/ZEROleafoutput/statusEXACT2^32/full64CBZ/source40byte/ordinal128/maskedradix/exactfullsizeledgerchecksremain. Guard10000/20secunchanged,noreset;prefix/hotcountsseparatelyhonest. InterpreterUBFXhasexactoneinstruction/width/read-before-write,notlimitrelaxation. FinitefixturesnotproductionGT/probabilitytheorem/controllerhardwarequalification;sharedcredential-freeDockeraloneexecutes/reviews/finalizes.


CurrentauditusesactualproductionNHprefixviewonARBITRARY192-bytewords,notolduncalledthree-chunkaffinereference. ItselectsactualOPSfrom.Lb_copythroughLASTNHADC,beforefirstaffineLSR,andappendsONEfixtureRET. Actual259prefixinstructions/260includingfixtureRET,ALLthree128-bitoutputscheckedvsindependentwrappingintegerformula,allzero/allMASK/randomwords/wrappingkeys,ZEROwrites/ABI/initializedreads/liveBpointer. Complete15coefficient/FOURcanonicalterms/LOW61/183bitrawsumoutputsremainindependentlycheckedinrun_keywithallnativefield/source/ABI/guards. OmitONLYnewpackageoldnh_fingerprint.SanditsOBSOLETEthree-chunkaffinetest;allproductionchecksretainedandcurrentarbitraryNHcoverageimproved. Historicalroot/frozen32filepackagesunchanged;thisnewpackage31filesretainslicense/allactualproductionproofs. ThischangesNOproductionhashbit/family/workcharge/scoreexceptexplicitlive21instructionremoval.

Paidworker_layout_init4sets27=48/28=96/29=actualp48;ALLsetup/delivery2^16ADDEDinsidecold2^124,separateR/notentropy/advice. Tableprefix318/42wordsunchanged,pluscold4=322. TWOhotkeyswithoutresettestlifetime. Fullproofinkey/fp2certificates;hashfamily/delta/success/tablevaluesUNCHANGED.

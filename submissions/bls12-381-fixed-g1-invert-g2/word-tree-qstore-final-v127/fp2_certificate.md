# Complete current recurrence: shared vectors, one modulus load, one public leaf

All residues are canonical full protected Fp381 Montgomery values Rm=2^384. Write c=A+B*t,X=a+b*t,Y=d+e*t in Fp4/Fp2,t²=xi=1+u,u²=-1. The exact recurrence factor is
A(a-d)+Bxi(b+e)+t[A(b-e)+B(a+d)].
Set V0=a-d,V1=xi*(b+e),V2=b-e,V3=a+d. Compute A*V0+B*V1+t*(A*V2+B*V3),plus the old third history. Only A and B coefficient tables are needed;the complete GT recurrence and all full field coordinates are unchanged.

## Actual fused wrapper and the public register contract

recurrence_fp2_step is the ONLY public PCS entry. Its528-byte frame saves19..29/LR before ANY clobber,and restores ALL after the map. SP/x18 are preserved. Four BOTH transforms construct a±d,b±e in initialized384 bytes atSP144..527. The fifth BOTH uses inputSP240/SP288 and writes differenceSP240/sumSP288,exactlyxi*(b+e). Every input is read before any output store;SUM25..30 remains live after MINUS stores. This exact crossed output alias is safe and touches neither raw histories nor other transformed components. ALL wrapper frame accesses useSP because BOTH clobbers29/LR.

LoadfullsixplimbsONCEinto19..24beforefirstBOTH;allFIVEtransformskeepthemreadonly. CurrentBOTHisMAX39withconditionaldiff/sumcorrectionandallbranches/joinBcharged. Originalsum/diff/fullpcoordinatesandallaliascontractsunchanged;completeuniversalpathproofbelow. Noexternalcallerpreloadingrequirement.

## All ordinary coefficient and fixed-s tables are constructed

For C=A,B,each vector entry is(d*Cr*2^start/Rm modp,d*Ci*2^start/Rm modp),96 bytes. Real input digits add(first,second) to(real,imag);imaginary digits SUBTRACTsecond fromreal and ADDfirst toimag. This is exactly C*(xr+xi*u) foru²=-1. Both output row pairs reuse these same A/B tables. No negated duplicate table or oracle is supplied.

ALL381 significant input bits have widths[35,35,35,35,35,35,35,34,34,34,34],starts0,35,70,105,140,175,210,245,279,313,347. One vector column is96*(7*2^35+4*2^34)=27TiB;two are54TiB. EVERY entry is constructed by actual full-p pinned seed products and existing fill_pair75(n-1)+17<=96n canonical additions,including all zero entries and headers. Full64-bit counters handle2^35/34 sizes. This is ordinary constructed memory,not virtual production arithmetic for free.

Vectors occupy8..62TiB;unchanged fixed-s13.5TiB is moved by its existing base argument to64..77.5TiB. Parent4TiB/private<1TiB/childports80..128TiB remain. All ranges are disjoint and pointers strictly below2^47. Other routing maps and DATA identifiers are unchanged. All new table pages/seeds/constants/MMU work is charged;bothworkers' vector/fixed fills remain<2^50 instructions inside completecold124. Bothworkers'tables<2^48 bytes/globalpeak138;no physical VM memory qualification is claimed.

## Signed accumulators stay nonnegative and below45p

Each real accumulator starts OLDthirdreal+22p. The ordinary seven-limb22p is already constructed in the50-slot header at528+64*22=1936. Exactly22 canonical negative table terms are subtracted. After k<=22 such terms,the accumulator is at least seed+(22-k)p>=0;positive terms cannot violate this. At most22 positive terms give real<seed+22p+22p<45p throughout. Imaginary has44 positive terms plusseed,<45p. All seven-limb add/sub chains are exact with no448-bit wrap. Final real is congruent to the desired native residue because22p=0modp.

The last real SBC intentionally preserves flags;the following imaginary ADDS starts a new carry chain. The shared interpreter implements that exact semantic. Range0..45p<49p preserves every premise of quotient_certificate.md:actual D=ceil(p/2^325),computed seven-limbq*p slots0..49,underestimate and one six-limb conditional subtraction. Every output is canonical0..p-1. Both old third seed coordinates are loaded before pair writes;the other pair's coordinates remain untouched until read. Exact third/output alias is safe;other partial overlaps are forbidden.

## One full public frame and two live accumulators

ThereisNOstandalonemaphelper.Themacrofileisexpandedintorecurrence_fp2_stepusingitsoriginal528-byteframe/outputSP96/headerSP128/tempSP144. Originalcaller19..29/LRvaluesaresavedonceandfullyrestored. Both7-limbaccumulators4..10/11..17staylive;QSTORE29includingstoresusesqp19..25,subtractsqp,THENreloadsp19..24forSIXlimbtrialsubtractandbranch-selectedTHREESTP. Theotheraccumulatorisuntouched;reducetheseconddirectlyin11..17. No spill/reload/copybridge or privateLRframeisneeded. Wrapperstack528;globalkey816unchanged. Fullregistered/canonicality/sourceproofbelow.

## Exact source counts and complete resources

Each COLUMN has6 scalar input LDR,11 extracts,5 crossing EXTR/5MOV,two index/base operations,11 MSTEP and10 advances:270 instructions. Width35 stride96*2^35=3TiB is oneMOVZ0x300<<32;34 advances useLSR1 in the SAME ADD. Column27TiB is oneMOVZ0x1b00<<32. Crossed digits load next limb,EXTR/mask into3 and retain nextlimb in19 before table loads clobber20. No digit is lost and no source read goes beyond six limbs.

Eachrowpairis9seed+13bias+1080columns+2header+29realQSTORE+29imagQSTORE=1162. Inlinecontext/constants9givescore2333;remainingwrapper241givesonepublicleaf2574. Caller20=2594retained2599with5spare. Everyinput/table/fieldcorrectionandoutputread/writeisretained.

Generation2599+1151keyhash+20aggregatescalar+704main=4474. Two local sort pairs648+scan128+cheap1+expensive1+spare256=5508;network906 giveshot6414,retained6492 with78spare. SAMEunequalfullL=88722074495999233395607723401466085376,N=6492L+2^124 satisfies EXACT N^500<20000^500*2^62219,proposed62219/500=124.438. Fullp381/r255,C20000/success0.39/sixtracks/onecomputationmetric,cold124/peak138/advice0 and exact0.3900120659047603 success remain unchanged. All failed/trial/preprocessing/random/transport/verification/output paths remain paid.

Shared-only audit asserts ACTUALpublicleaf<=2574/nohelpercall against an independent DENSE full-p381Fp4 oracle. ABI assertions for19..29/LR/SP/x18 are at PUBLICstepreturn,where the wrapper restores them. The macrocore is within the public frame;no standalone helper ABI is asserted. Readonly initialized virtual table fixtures have current two27TiB/11digit geometry and fixed64TiB;they are scoped test oracles,never substitutes for charged production construction. Other full-p/key/369hash/fullY/U/corner/full-rscalar/byte/copy/packet tests remain. These are not productionGTinversion/tablefill/controller timing/probability/hardware qualification. No new pass/score/acceptance/promotion is invented;candidate code was NEVERhostexecuted/imported/compiled. Shared credential-free llm-verifier alone executes/reviews/scores.


The actualmacro-onlyrecurrence_fp2_map.Sincludeisexpandedinsidepublic528frame. OriginalargsoutSP96/tableSP120/headerSP128establishcorecontextbyLDR/ADD/LDR;two3instructionconstantinitializationsandtwoROWPAIRsuseinitializedSP144..527. ParseracceptsONLYthisliteraladjacentinclude;nativeassembler/sharedfixtureexpandthesameMSTEP/COLUMN/QSTORE/BOTHbodies. Codecopy/bootstrapischargedcold. OriginalLRsavedSP8/restoredatenddespiteBOTH/tableusing30. NoBL/extrahelperframe.

# Conditional correction with explicit worst-path scalar instruction bounds

This revision changes ONLYnativecanonicalFpaddition/subtractioncorrectioncontrol. Inputsremainfull-p381canonicalnativeMontgomerya,b in[0,p),exactalias/disjointmemory/publicABI/nativeaddress/immediatesunchanged. Sum<2p<2^384,so six-limbcarryout iszero. Differencein(-p,p)hasexactunsignedborrow;addingpafterborrowgivescanonicalresult. Data-dependentbranchesareordinaryscalarinstructionscharged1EVERYtime,notfreebranchprediction/CPUtiming. ThereisNOconstant-timeexecutionclaim/requirement andNOdistributionassumptionneeded:costboundsareMAXIMUMoverALLbranchpaths.

## Public fp381_add and fp381_sub

fp381_addloadsall12inputlimbsBEFOREanystore,computessumintox4..9,loadsordinaryp,trialsubtractsintox10..15. FinalSBCSborrowselectscontrol: B.LOoriginalstorepathifsum<p,otherwise storetrial. BOTHpathscontainthreeSTP+RET. ThereisNOsix-CSELsequenceandnoextraunconditionaljumpbecauseeachpathreturns. Max=6inputLDP+6sumops+3pLDP+6trialSUBS/SBCS+1conditionalbranch+3STP+1RET=26. Equalitysum=pcorrectlystoreszero,zero+zerostoreszero. No extra registercopy/outputstore/returnisomitted.

fp381_subloadsallinputlimbs,subtractssixlimbsintox4..9,then B.HScommonstoreifNOunderflow. Otherwise reloadpintox10..15andexecuteSIXcarryaddinstructions(unmaskedp),then common3STP+RET. Negativepathmax6inputLDP+6subtract+1branch+3pLDP+6add+3STP+RET=26. Nonnegativepath17. Originalborrowneednot survivepLDP,becausefirstADDSstartsfreshcarry;negativemod2^384+pwrapsEXACTLYa-b+p< p. Callee19..29/LR/SP/x18untouched. Memorya/b/preadonly;allinputloadsbeforestorepreserveexactout/aorout/baliases. Alloutput48bytesinitialized.

## Fused BOTH macro within the public recurrence leaf

BOTHretainsfullplive19..24,fusedsum25..30anddifference4..9. AfterdifferenceSBCS,B.HSskipssixdirectpaddinstructions;negativepathexecutesthemwithouttheoldsixp-CSELs. ThenTHREESTPstorecorrecteddifference;sumregistersremainintact. Trialsubtractsum-pintox4..9;B.LOselectsthreeoriginalsumSTPotherwiseTHREEtrialSTPplusONEunconditionalBtojoin. Thatjumpisfullycharged. Thetwooutputregionsaredisjoint/exactdocumentedcrossaliases;allinputloadsprecedeanywrite. CarryflagsfromdifferencecorrectionaredeadbeforesumtrialSUBSresetsflags. BothoutputsarecanonicalforALLfullpinputs. Everymacroinvocationusesanexplicituniquesuffixa..e,solabelsareuniqueandbothnativeassembler/sharedparserexpandtheSAMEbranches.

WorstBOTH=6inputLDP+6sumops+6diffsubtract+1diffbranch+6diffadd+3diffSTP+6sumtrial+1sumbranch+3sumSTP+1joinB=39,versus48. ThereareFIVEBOTHs,soatmost195ratherthan240. Remainingwrapper46unchanged=>wrapper241,inlinedmapcore2333unchanged=>publicrecurrence<=2574throughRET. Caller20=>2594retained2599,fivespare. Someinputsarecheaper;NEVERassertsameexactcountforallcases. Thisisacyclicpathmaximization,notCPUbranchfrequencyorlatencymodeling.

Odd55own+THREEpublicFphelpers<=3*26+twofixed207=547. Fusedkeyhash<=93own+2*65sign+547odd+369NH=1139,caller9=>1148ret1151. IdentityNOhash<=86+130+547=763. Actualstatickeyown92gives1138;retained1139conservative. FixedMAX207/QSTOREMAX29includingstores/Dceil/current192logicalkey/24inputwordtrace/hashfunction/stack816/raw/scalarABIremainunchanged. RecurrencehasNOBL,wrapper528/full19..29/LR/SP/x18preservationunchanged.

GEN2599+1151+20+704=4474. Currentdirect64/relayed96localsort648+scan128+cheap1+expensive1+spare256=>base5508. Network906unchanged=>hot6414ret6492with78spare. SAMEfullL88722074495999233395607723401466085376,N=6492L+2^124satisfiesEXACTN^500<20000^500*2^62219,proposed62219/500=124.438. Fullp381/r255/suppliedP/uniformGT/noauxpowers/full160U/128ordinals/24tag/88payload/40metadata/finiteexactsuccess0.3900120659047603/success0.39/C20000/sixtracks/onemetric/cold124/peak138/advice0unchanged. ALLnotifications/phasebarriers/bucketconstruction/verification/failurepathsremainpaid.

Shared-onlyauditinterpretsACTUALconditionalnativeinstructionsatfullpboundaries/randomcasesagainstindependentdenseFp4/product/exactkeyoracles. Recurrencecountisnow<=2574,odd<=547,key<=1139;publicFphelpercalls<=26ratherthanfalseexact31. Somefixtureinputsdonotexercisebothbranches,theuniversalmaxpathcertificateaboveistherealbound. Supplementalseparateboundarycasescheck0,1,p-1,halfp,sump,borrow,alias/allpublicABI withactualfieldsource. ExistingALL24actualNHinputwordtrace/369hash/207fixed/180QSTOREstored-output/nativecopy/scalar/6496sort/rolelayout/notificationfixturesremain. CandidatecodeisNEVERhostexecuted/imported/compiled/emulated;sharedcredentialfreellm-verifieraloneexecutes/reviews/scores. Agent-ownedSTATICmacroexpansion/acyclicCFGpathcountdoesNOTinterpretregisters/CPUvalues. Noobservednewpass/score/promotionisclaimedbeforeactualreview.

The samebranchconstructionalsoappliestoFULLFr255giantaddition:4inputLDP+4sumcarry+2rLDP+4trialSUBS/SBCS+1B.LO+2STP+RET=MAX18. Bothstorebranchesreturn,soNOjoinjump. ALL32bytes/4limbsretained,2r<2^256,nolostcarry;allinputsbeforestores/exactalias/publicABIunchanged. Caller5=>23retained25. Baby7/callerret12unchanged;13ng<=8Lproves12nb+25ng<=20LforEVERYrunundercurrentfullcounts,notprefixaveraging. Generation2599+1151+20+704=4474;localsort648/scan128/cheap1/exp1/spare256=>5508,network906=>hot6414ret6492with78spare. EXACTN=6492L+2^124satisfiesN^500<20000^500*2^62219,proposed62219/500=124.438. Actualshared-onlyFrfixturescoverzero/sumexactr/2r-2/random/aliasandbound<=18;fullrparameterunchanged.

CURRENTQSTOREMAX29INCLUDESselectedTHREESTP/conditionalbranch/joinB;noSIXCSELregistercopy. Outputmemorycanonical,oldaccumulatorregistersnotrequiredcanonicalandnotreadagain. FourQSTOREsave16recurrence/MAX2574ret2599;fixedMAX207/odd547/fused1139ret1151;gen4474/base5508/net906/hot6414ret6492. Detaileduniversalrange/alias/liveness/actualstoreboundaryproofinquotient_certificate.md.

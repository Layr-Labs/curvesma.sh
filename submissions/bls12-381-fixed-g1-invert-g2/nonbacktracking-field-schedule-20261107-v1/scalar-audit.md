# Complete scalar work for the capped chained dictionary

One accepted transition uses one protected pinned dense Fp12 reference multiplication at C=20000,charged once. All other work is scalar AArch64 instructions,including calls,returns,loads,stores,serialization,barriers,word arithmetic andcontrol. Candidate evidence is not executed on the host.

The inherited exact helpers retain canonical-work budget6474K (canonical-budget.md),simpletabulation807 perstep,scalarupdate46,andstandard full seven-round BLAKE3 final797 percall. There is ONE ordinaryfullhash:JUMP,797 total. Its256bits givej76,V180,phasefloor12V/2^180 andbucketV mod2^128. Fulloutputsarecomputed. TwoEXTR/oneLSR andthree-limbphase multiplicationfittheretained20 phaseinstructions;callerwordloads/copies/interfacesfitcontroller848. Thisretainstheold784 plus64 forcap/contextwork. ThereisnoordinaryPHAScall/prefix.  Retry still uses generic48-bytefinalblocks andthe independent priorhardcap,whose extra cost<K isretained as1 perstep.

The canonical helper still has worst-case7024,but its explicit orientation-weighted cap bounds total6474Kplusonefixedfailedhelper. It includes six always-used native Fp multiplications,andcheap copy/negation branches withthe remaining weightedcost. field-ops.S andcopy-negate.S supply the same concrete leaves. Allfieldvalues andhash labels remainunchanged. canonical-min.S/canonical-min.md replace768 comparison/copy work by384:eleven28-instructionborrow-chainpointercomparisons plus19fixed instructions=327,withinterfaces343<384. Candidatearena1152bytes isshared;onlythewinneriscopied. Othergeneration/key/control/wrapperbudgetsremainpaid.

Dictionary head lookup/publication192,allocatedpayloadwrite224,andatmost floor(K/4)+2 chainvisits at128 eachgive32K+256. dictionary.md enumerates every selector,read,write,comparison,counter andinterface;the fixed256 iscold. No expected-only probe price isused. Exactbcomparison isatmostonce andterminalcold. All512 allocatedpayloadbytes areinitialized;all2^128 heads areinitializedcold. No tree path,split,leafshift orfreesort remains.

Totalordinary work<=
 (6474+807+797+46+848+192+224+32+1)K+fixed256
 =9421K+fixed256.
Retain9440K,19 slack perstep. Theinitialkey/hash/head/payload andfixedfailedhelper arepaidcold. Everycostcapisalgorithmic,checkedbeforeworkthatwouldexceedit,anditsfailureprobabilityissubtractedseparately. Countersfit128bits;logicalindicesusetwo64-bitlimbs. NoSIMD,cryptoextension,wideintegerunitprimitive,freeoracle orparallelwalltime discount isused.

Cold(540*2^119) bounds C*H_cold+R_cold ONCE. Previous completegenerator/scalar/variant/hashconstruction,allattempts,writes,pairings,validation,constants,frames,initialedge,stateandcanonicalkey,failedfinishes,Frinversion,G2outputandexactverification fit<2^120. Fullheadinitialization is537*2^119,selectors<=32*2^92,index/setup<2^100,andpriorcompletecold<2^120;theirsum<(540*2^119). block-copy.S certifieszero8192=512STP+RET,andzero512=32STP+RET. Theircaller/loopworkisexplicitlyincluded. No pre-zeroedmemory,advice orfuturetargetamortizationisassumed.

Thus W<=K+(9440K+(540*2^119))/20000=(29440K+(540*2^119))/20000.
The exact protected integer inequality(29440K+(540*2^119))^1000<20000^1000*2^126239 provesproposed126.239,above diagnostic126.23814349539336. This isanalyticconditional evidence,notanofficialscore,productionexecution,hardwareavailability orqualification.

field-schedule.md replaces19rotationoperationsbytheexplicit16,and12structured-productfieldoperationsby9. Copyxintoslot0costs12data+2pointerinstructions,bringingtheminimum/copycertificate327+16+14=357<384. Nonmultiplication2864=9*48+16*48+384+768+512. Worsthelper7024,orientationbase5360,andunchangedclasscountergivehardcanonical6474K. The newcompleteledger9421 isretainedas9440;thegroup/hashvalues andallprobabilityboundsareunchanged.

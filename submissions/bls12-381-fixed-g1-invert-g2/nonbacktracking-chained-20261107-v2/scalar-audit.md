# Complete scalar work for the capped chained dictionary

One accepted transition uses one protected pinned dense Fp12 reference multiplication at C=20000,charged once. All other work is scalar AArch64 instructions,including calls,returns,loads,stores,serialization,barriers,word arithmetic andcontrol. Candidate evidence is not executed on the host.

The inherited exact helpers retain canonical-work budget7232K (canonical-budget.md),simpletabulation807 perstep,scalarupdate46,andstandard full seven-round BLAKE3 final797 percall. There are now THREE ordinary hashes:JUMP,PHAS,DICT,2391 total. All hash output words are computed. DICT uses a separately cold-cached64-byteprefix andthe unchanged40-bytefinalblock. Caller CV copies and interfaces are paid incontroller848;this increases the old784 by64 for additional hash/cap/context work. Retry still uses generic48-bytefinalblocks andthe independent priorhardcap,whose extra cost<K isretained as1 perstep.

The canonical helper still has worst-case7696,but its explicit orientation-weighted cap bounds total7232Kplusonefixedfailedhelper. It includes six always-used native Fp multiplications,andcheap copy/negation branches withthe remaining weightedcost. field-ops.S andcopy-negate.S supply the same concrete leaves. Allfieldvalues andhash labels remainunchanged.

Dictionary head lookup/publication192,allocatedpayloadwrite224,andatmost floor(K/2)+2 chainvisits at128 eachgive64K+256. dictionary.md enumerates every selector,read,write,comparison,counter andinterface;the fixed256 iscold. No expected-only probe price isused. Exactbcomparison isatmostonce andterminalcold. All512 allocatedpayloadbytes areinitialized;all2^128 heads areinitializedcold. No tree path,split,leafshift orfreesort remains.

Totalordinary work<=
 (7232+807+2391+46+848+192+224+64+1)K+fixed256
 =11805K+fixed256.
Retain12288K,483 slack perstep. Theinitialkey/hash/head/payload andfixedfailedhelper arepaidcold. Everycostcapisalgorithmic,checkedbeforeworkthatwouldexceedit,anditsfailureprobabilityissubtractedseparately. Countersfit128bits;logicalindicesusetwo64-bitlimbs. NoSIMD,cryptoextension,wideintegerunitprimitive,freeoracle orparallelwalltime discount isused.

Cold2^130 bounds C*H_cold+R_cold ONCE. Previous completegenerator/scalar/variant/hashconstruction,allattempts,writes,pairings,validation,constants,frames,initialedge,stateandcanonicalkey,failedfinishes,Frinversion,G2outputandexactverification fit<2^120. Fullheadinitialization is537*2^119<2^129 instructions plusselectors/indexsetup<2^100;sumwithpriorcold<2^130. block-copy.S certifieszero8192=512STP+RET,andzero512=32STP+RET. Theircaller/loopworkisexplicitlyincluded. No pre-zeroedmemory,advice orfuturetargetamortizationisassumed.

Thus W<=K+(12288K+2^130)/20000=(32288K+2^130)/20000.
The exact protected integer inequality(32288K+2^130)^1000<20000^1000*2^126373 provesproposed126.373,above diagnostic126.37199797398533. This isanalyticconditional evidence,notanofficialscore,productionexecution,hardwareavailability orqualification.

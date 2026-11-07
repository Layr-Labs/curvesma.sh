# Complete simple plain-walk path

Reference C=20000 isunchanged. One denseFp12 multiplicationpertransitionischargedonce;everythingelseisinthescalarcount. Staticassemblyisneverhost-executed.

Main50 includesall BL/arguments,referenceinterface,pointer swap,128-bitcount/markerchecks,CVcopies,key/tabulation/hashpreparation,label extraction/cache,full-keylookupanddispatch. No signcall,prefixcache orcomparison-budgetpathoccurs. Coldblockpadding andkeyword36areinitializedonceandneveroverwritten. Lasttoptest andboundedcheckpoint/terminalhandlersarefixed2^90.

Scalar46 is5setup+20perresidue+RET. BaseADD#576makesallpairedLDP/STPoffsets0..48encodable. Bothinputresidues<r<2^255;SUBS/SBCScandidate reduction/CSEL iscanonicalwithnofifthcarry. Copy37is18LDP/STPpairs+RET,writingexactly36bwords640..927. Constantzero928..935 isuntouched. Pairreversalpreservesthedefinednumericalwordorder.

Tabulation807 is2baseADD+4constantT72[0]loads+36PAIR*22+5REV+3stores+RET. EachPAIRhasadjacentdigit/baseLDPs andtwo10-instructiontableoperations. Registers13/14preservethesecondpairwhile8..12aretemporary;15/16prebase thelastfourpairedloadstokeepimmediatesencodable. Everytablewordread,address,XOR andoutputstoreiscounted. Bothconstanttableentriesareconstructedcold.

BLAKE3label779 has734mixing+45fixedinstructions,includingframe,pointers,preloads,stateinitialization,requiredstoresandrestore/RET. Itcomputesexactlystandarddigestwords6/7throughallsevenrounds. Deadfinaldiagonalupdatesc10,b5,c9,b4andsixunusedoutputpairsareomittedasexactdependency slicing,not reducedrounds. Thecaller'slow40-bit extractionisunchanged. Genericcoldhashescomputefulloutputs. Frame96preserves19..30,x18untouched.

Lookup703 isentry2+7iterations*85+lastlooptest2+terminalbounds2+cursors3+equality94+return3+slack2. Orderingcomparison75=18pairs*(2LDP+2borrowinstructions)+final3. Iterationcontrol/mid/address/branchesfit10more. MidpointADD/LSR isexactbecauselo+hi<=254. Equality94=18pairs*(2LDP+CMP+CCMP+B.NE)+final4. CCMPsecondwordunderEQ,elseNZCV0,forcesinequalityonafirst-wordmismatch. Earliermismatchischeaper. Full37-wordequalityincludesconstantzeroandneveracceptsfingerprintsalone. Callerregistersonly,nospill/call.

Total50+46+37+807+779+703=2422;declare2432,tenperstepslack. Allinterfacesarealreadyinmain50. Cold2^90payscompleteinitialization,construction,trials,failures,checkpointcopies/insertion,verificationandG2output. Work<=(22432K+2^90)/20000. Code/stack/pointersfitthesource'sordinary<2^50byteaddressregion;noadditionalmemorymachineassumptionisintroduced.

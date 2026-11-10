# Current complete native key/hash and paired-field source path

Canonical field/Frobenius fibers are proved in key_contract.md;NH/affine in nh_certificate.md. Actual combinedleaf hasNOoutputmemoryaccess. Inputs:compatibility0unused/rawPhi192at1/sLUT2/p483/products4/keys672at5. Nonidentity returnsALL183bits word0=x1/word1=x0/word2=x5low55;statusx7nonnullselectedBpointer. Identityx7zero invalidatesresults. Source MUST useFULL64CBZx7 andnever dereference the expiredframepointer.

## Complete field normalization/hash dataflow

Bothlow-limbsignblocksareactualnativecanonicalp48subtractions/zerooutlines. OddconjugationdoesnotchangerealA,so its copy is dead;canonicalrealAalwaysrawPhi. Compare signedBtwelveLE64words,then choosecorrespondingimagA. ReadALL192conceptualbytes in orderrawAreal48/selectedAimag48/signedB96. No independent field bit is dropped.

NH row1 accumulators0/1keepselectedpointers6/7live duringALLinputreads. Compressaffinerow1BEFORErow0,holdingresult0while4/5row0accumsurvive;row0result1,row2result5. EachrowhasfreshADDSchain/independentcoefficients,so no carry crossesrows. FinalBFI/LSRreturnswordsx1/x0/x5withoutoutputstores. The full128-bitwrappingNH/M61family isunchanged.

Entryp48pointer3remainslive throughfirstsign/oddimagnegation/bothfieldoperations;nonewrites3. THREEpLDPsfornegationread3directly;noearlyMOV21,3ortwoMOV3,21. ExistingfirstfixedargumentLDP2,3atSP152laterreplaces3withproductsbase;bothfixedMOV21,3remain. No freeconstantdeliveryorABIassumption.

## Paired canonical sum/difference

ReadrawBreal4..9 andBimag10..15ONCEbySIXLDPdirectraw20offset96..176. ONE STP25/26atSP0preservescallercallees. SixADDS/ADCScomputesum21..26;canonicalinputs<p<2^381make sum<2p<2^382,no384overflow. SixSUBS/SBCScomputediff4..9. THREEactualpLDPsload10..15withoutchangingborrow. SixCSELchoosepversuszeroonLO,then sixADDS/ADCSgivecanonical(a-b)modp;threeSTPstoreSP240/256/272.

RELOADpbyTHREELDPs(no freeconstantreuse),sixSUBS/SBCStrialsubtractfromoriginalsum21..26into10..15. Finalcarrymeanssum>=p;sixCSELchoose trial/original,threeSTPstoreSP192/208/224. ONE LDP25/26restorescallerBEFOREfixedproducts. All48bytesofeachcanonicalfieldtemporaryarewrittenexactly. Twofixedproductsreadthesamebuffersasbefore.

Oldtwo30instructionbodies+6addressADD=66. Newbody54includesALLloads/arithmetic/conditions/stores,plusnew25/26backuprestore2=56. TENactualinstructionssaved. DirectLDP/STPoffsets96..272allencode;noaddressoperationishidden. StandardcanonicalfieldarithmeticisnotconditionedonGTspecialvalues.

## One complete frame,all ABI and count

ONE576frame:innerbackup21/22at32,23/24at48,args152/160,scratch192/240;new25/26backup0..15. Header288..335holdsoriginal19/20at288,table/products304/312,p48at320,keys328. EvenB336..431,OddAimag432..479,OddB480..575. OddnominalrealA384..431isNEVERaccessedandmayoverlapEvenBimag;canonicalrealAisrawPhi. Inneroddout19=SP384. Allactualaccessesareinitialized/disjoint/encodable. Frames/pointers16/8alignedasrequired.

Allocate576once,restoreonceatONLYRET. NO BL/nox29/x30writes;entryLRstayslive. Callees19..26areexplicitlypreserved,27..30untouched,18neverused. Successpointerx7=SP336or480=entrySP-240or-96. Avalid576downwardframeforcesnonnull. NHonlyREADS7afterselection,epilogueleaves7unchanged. Identitysets7zero. No expiredframe dereferenceoccurs.

Earlierraw963 removed2hashSTP+1outMOV+2inner19backup,then1statusMOV,2nestedSPoperations,3modulusMOV,10pairedsumdiff:963-5-1-2-3-10=942. RawDAGcountsALLsyntacticbranches/RET,rejectscycles. Retain963=942+1margin+12caller+8spare. Field2506+key963+Main720=4189;local720+scan128+confirmation1+inner256+network886+outer80givesN6260L+12m+28n+cold2^124. Allfullparameters/probability/coins/cap/prices remain.

Worker19out/20color/21-22fullordinalareABI-preserved. ActualsourceCBZx7skipsidentity;successSTPlo128,ORRcolor55,STRhi16,STPordinal24completesALL40bytes. Fourdatainstructions,ordinalSTPalreadyinmain712;newsource3. Branchalreadyincaller12. NativehelperincludesfixtureRETonly;noextraMOV/BL/argloadpricedaway.

IntegratedisolatedoracleusesZEROleafoutputwritebudget,fullindependentcanonicalFp+full128NH/M61registerhash,andactualcallee/SP/rawPhi/initializedreadchecksincluding25/26. EVERYnative successcasechoosesvalidbelow47-bitstackgivingstatusEXACT2^32(low32zero),thenexecutesACTUAL64bitCBZ/sourcefourdatawithSAMEreturned/preservedregisters,checkingall40bytes/hash183/color/pending0/ordinal128. IdentityCBZ+RETstoresnothing. Referencefield/copy/byte/scalar/controlandexactfullsizeledgerremain. Thesearescopedfixtures,notproductionGT/probabilitytheorem/controllerhardwarequalification;sharedreviewisrequired.

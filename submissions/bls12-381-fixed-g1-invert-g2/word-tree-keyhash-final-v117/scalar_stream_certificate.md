# Unequal full-size attack streams and actual scalar preparation

Current baby count nb=35597089187900144514297950209050673152,giant count ng=53124985308099088881309773192415412224,L=88722074495999233395607723401466085376=nb+ng. Counts are9141 and13642 times3*2^110. Both remain multiples3/512/2048;the full curvep381/orderr255 and uniformGT input distribution do not change. These are finite algorithm loop lengths,not a reduced challenge. Every generatedGTvalue remains actual h^j orz*h^(-bi). Full canonical Fr32 scalars and ALL128 ordinal bits are retained. Babyordinalj-1;giantordinalnb+i;no gap/overlap/duplicateordinal.

## Actual unchanged seven-instruction baby helper

baby_scalar_fill in metadata_io.S takes x0disjointfull32-byteoutput and x1initialized16-bytepublicjcounter. It emits LDPjlow/high,STPbothlow16out,STPzerosupper16,ADDSlowcarry/ADChigh,STPnextstate,RET:7instructions. Twoargumentmoves+BL3 give10 retained12. It initializesALL32scalarbytes;since1<=j<=nb<2^126<r,the result is exactlyfullcanonicalFrencodingj. j+1 neveroverflows128bits,includinglow64carry. Outputstatepartial/exactaliasforbidden;callee19..29/LR/SP/x18unchanged. The actual helper source is identical to publicv92,not a new smaller field. The current outputpointer ispayload32 in the88byteU32/scalar32/tag24 layout.

AdvancecounteronceperNEWpreparedrecord,notperpayload/metadataRPCreply. Separatefull128metadataordinalcounteradvancesunderMain704,andpendingflagpreventsdoublegeneration. Coldinitialbabytracerecordsuseidenticalscalaremissionandareconservativelyhotcharged. No stale32-bytebabyFrstateisread;thejstateis16bytes. Endj=nb emitsnb andstoresnb+1safely.

## Full-r giant addition and all failure prefixes

Giants retainactualfr255_add.S21instructions throughRET,full255-bitcanonicalFrstatebi,correctfourlimbcarry/trialsubtract/selectionandall32outputbytes. Fourargumentmoves+BL5give26retained28. Since2r<2^256,nocarryislost;exactinput/outputaliasisretained. Fullstate/payloadcopyandindexcontrolremaininMain704. The finalunusedincrementisalsocharged.

For EVERY run,numberbabyemissions<=nb andgiantemissions<=ng. Therefore scalarwork<=12nb+28ng. At these exact counts,
12nb+28ng<=22(nb+ng)=22L,
since ng<=5L/8. Actualaggregateis21.580476671202213L;22L is a conservativeintegerupperbound. It is NOT20L and NOT an average over presumedbalancedearlyprefixes. An abort reducescountsandcannotincreasework. Othergeneration4576L gives4566L. Two localsortpairs720+scan128+cheap1+expensive1+spare256 givesbase5672. Network933 giveshot6605,retained6680 (75spare). FullN=6680L+2^124 satisfiesEXACT N^1000<20000^1000*2^124479,proposed124479/1000=124.479. All other operation/referencepricesareunchanged.

## Group-event/probability adaptation,not a new heuristic

proof.md replaces babyshortintervalm bynb ONLY,and giantindexm byng ONLY. Its exacthypotheses hold:nb<u0²-1,2nb²<r,nb*u0²<r-nb,u0*nb+nb<r,3ng<lambda0,4ng(lambda0+1)<r. Babyunion s=12(nb-floor(nb/u0));mixedlineeventsM=ng*s;fourdisjointgiantfamiliesG=5ng²/3-ng because ngdivisible3. Mixedlineintersections/guardedgiantslopeindependencefactsareunchangedfiniteFrlinearalgebra. Bonferroni mu=M/r+G(r-1)/r²,U2=(choose(M,2)+M*G)/r² remainsvalid. Deterministicbabyduplicatesfloor(nb/u0),identitylossng/r,cheapandstrongtagcapsuseCURRENTL. Strongtageta=(2^-63+1/q)^3shareskeyswithUandisnotmultipliedbydelta. Exactcompleteboundis0.3900120659047603>3900031/10000000>0.39.

All2^80sortarrays/routing/ports/coldtablesremainactuallyconstructedandcharged. SmallerLreducesstreamvolumebutdoesnotmaketheexplicitcoldallowancefree. Bootstrap121/edges118/sortexception104/otherboundsstillgivecold124;peak138/advice0unchanged. Finitepaidcoins/R sourceandallfalse/ignored/failedpathsremaincharged. ActualGTverificationcontinuesafterfalseCtagpairsandhasseparate2^36hotpaircap;noonly12totalcoldattemptclaim.

Shared-onlyaudittestsactual7instructionbabyhelperatj=1,nb,2^64-1,2^64 forALL32encoding/exactnextstatecarry/initializedread/writeranges/ABI,andactual21instructionfull-rgiantadditionzero/wrap/max/randomwithalias. Existingactualnative369hashfullY/Uand88-bytepacket/copy/getterfixturesremain. TheseareSCOPEDchecks,notproductioninversion/GTprobability/nativecontrollerhardwarequalification. No newpass/score/promotionisinvented;candidatecodewasneverhostexecuted/compiled.
Currentproductionkey+hashisONEretained1176call,not824+376. FullS192islogicalA96||B96readentirelywhilekey384frameislive,no192-bytecopy/materializedscratch. Publicsixargsincludehashkeys;identityreturnsbeforeNH;payloadscalar32andrawhistoriesareuntouched. Exactfunction/byteorder/lifetime/allpath/ALL24actualreadwordproofsareinkey_path_certificate.md. Generation2664+1176+22+704=4566.

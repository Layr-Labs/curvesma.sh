# Unequal full-size attack streams and actual scalar preparation

Current baby count nb=35597089187900144514297950209050673152,giant count ng=53124985308099088881309773192415412224,L=88722074495999233395607723401466085376=nb+ng. Counts are9141 and13642 times3*2^110. Both remain multiples3/512/2048;the full curvep381/orderr255 and uniformGT input distribution do not change. These are finite algorithm loop lengths,not a reduced challenge. Every generatedGTvalue remains actual h^j orz*h^(-bi). Full canonical Fr32 scalars and ALL128 ordinal bits are retained. Babyordinalj-1;giantordinalnb+i;no gap/overlap/duplicateordinal.

## Actual unchanged seven-instruction baby helper

baby_scalar_fill in metadata_io.S takes x0disjointfull32-byteoutput and x1initialized16-bytepublicjcounter. It emits LDPjlow/high,STPbothlow16out,STPzerosupper16,ADDSlowcarry/ADChigh,STPnextstate,RET:7instructions. Twoargumentmoves+BL3 give10 retained12. It initializesALL32scalarbytes;since1<=j<=nb<2^126<r,the result is exactlyfullcanonicalFrencodingj. j+1 neveroverflows128bits,includinglow64carry. Outputstatepartial/exactaliasforbidden;callee19..29/LR/SP/x18unchanged. The actual helper source is identical to publicv92,not a new smaller field. The current outputpointer ispayload32 in the88byteU32/scalar32/tag24 layout.

AdvancecounteronceperNEWpreparedrecord,notperpayload/metadataRPCreply. Separatefull128metadataordinalcounteradvancesunderMain704,andpendingflagpreventsdoublegeneration. Coldinitialbabytracerecordsuseidenticalscalaremissionandareconservativelyhotcharged. No stale32-bytebabyFrstateisread;thejstateis16bytes. Endj=nb emitsnb andstoresnb+1safely.

## Full-r giant addition and all failure prefixes

Giants useactualfr255_add.SMAX18throughRET,full255-bitcanonicalFrstatebi. Loadallfourlimbsofaandb,carryaddsum,loadfullr/trialsubtract,ONEB.LOchoosesoriginalsumstoresortrialstores. EachpathhasTWOpairedstores+RET,noextrajoinjump. Allfourcarry/triallimbsretained;sum<2r<2^256andall32outputbytesinitialized. FourargumentMOV+BL5=>23retained25. Exactout/inputalias/public19..29/LR/SP/x18unchanged,allinputsbeforeresultstores. State/payloadcopy/indexcontrolremainMain704andfinalunusedincrementispaid.

ForEVERYrun,babyemissions<=nb/giantemissions<=ng. Thusfullscalarwork<=12nb+25ng<=20(nb+ng)=20L,since13ng<=8L(atfixedcountsng=13642units,total22783,13*13642=177346<=182264=8*22783). Thisis a fullproductionLbound including everyabort,NOTanaverageoverprefixcolorfractions or a claimscalarwork<=20timesanindividualprefixlength. Smalleractualcountsnevertightenrequiredaccounting. Fullr/32bytesandbaby7areunchanged. Othergeneration4478L+scalar20L=4474L. Localsort648+scan128+cheap1+expensive1+spare256=>base5508. Network906=>hot6414ret6492with78spare. N=6492L+2^124EXACTN^500<20000^500*2^62219,proposed62219/500=124.438. NoGTgroup/hash/cap/successchange.

## Group-event/probability adaptation,not a new heuristic

proof.md replaces babyshortintervalm bynb ONLY,and giantindexm byng ONLY. Its exacthypotheses hold:nb<u0²-1,2nb²<r,nb*u0²<r-nb,u0*nb+nb<r,3ng<lambda0,4ng(lambda0+1)<r. Babyunion s=12(nb-floor(nb/u0));mixedlineeventsM=ng*s;fourdisjointgiantfamiliesG=5ng²/3-ng because ngdivisible3. Mixedlineintersections/guardedgiantslopeindependencefactsareunchangedfiniteFrlinearalgebra. Bonferroni mu=M/r+G(r-1)/r²,U2=(choose(M,2)+M*G)/r² remainsvalid. Deterministicbabyduplicatesfloor(nb/u0),identitylossng/r,cheapandstrongtagcapsuseCURRENTL. Strongtageta=(2^-63+1/q)^3shareskeyswithUandisnotmultipliedbydelta. Exactcompleteboundis0.3900120659047603>3900031/10000000>0.39.

All2^64+2^96sortarrays/routing/ports/coldtablesareactuallyconstructedandcharged undercurrentdirect64/relayed96roleandaddressmapping. SmallerLreducesstreamvolumebutdoesnotmaketheexplicitcoldallowancefree. Bootstrap121/edges118/sortexception120/otherboundsstillgivecold124;peak138/advice0unchanged. Finitepaidcoins/R sourceandallfalse/ignored/failedpathsremaincharged. ActualGTverificationcontinuesafterfalseCtagpairsandhasseparate2^36hotpaircap;noonly12totalcoldattemptclaim.

Shared-onlyaudittestsactual7instructionbabyhelperatj=1,nb,2^64-1,2^64 forALL32encoding/exactnextstatecarry/initializedread/writeranges/ABI,andactual21instructionfull-rgiantadditionzero/wrap/max/randomwithalias. Existingactualnative369hashfullY/Uand88-bytepacket/copy/getterfixturesremain. TheseareSCOPEDchecks,notproductioninversion/GTprobability/nativecontrollerhardwarequalification. No newpass/score/promotionisinvented;candidatecodewasneverhostexecuted/compiled.


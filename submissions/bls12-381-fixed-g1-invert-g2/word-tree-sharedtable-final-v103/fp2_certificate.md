# Shared coefficient vectors, xi input transform and signed bounded accumulators

All full protected Fp381 arithmetic uses canonical six-limb native Montgomery residues Rm=2^384. Letc=A+B*t,previousX=a+b*t,Y=d+e*t inFp4/Fp2,t²=xi=1+u,u²=-1. Asbefore,
cX-c^(p²)Y=A(a-d)+Bxi(b+e)+t[A(b-e)+B(a+d)].
DefineV0=a-d,V1=xi*(b+e),V2=b-e,V3=a+d,ALL inFp2. The output is A*V0+B*V1+t*(A*V2+B*V3),plus the oldthirdhistory. This moves xi from the coefficienttable to an actual inputtransform;no algebraicinformationisdiscarded.

## Actual fused wrapper and exact alias safety

The originalfour51-instructionBOTH macros constructa±d andb±e in384 initializedtemporarybytes atSP144..527. Add aFIFTHBOTH withx1=SP240,x2=SP288,x0=SP288(plus),x3=SP240(minus),p atsavedSP136. It reads both real/imaginarycomponents ofoldV1=b+e beforeanywrite and replaces themby(real-imag,real+imag),exactlyxi*V1. Allinputreadsprecedeoutputstores;SUMremainslivein25..30whenMINUSisstored. Exactcrossedoutputaliasissafe. Five setupinstructions+51BOTH give56extraowninstructions. Wrapperown250+56=306. No rawhistory orotherVcomponentisoverwritten. Every frameaccessusesSP becauseBOTHclobbers29/LR;callee19..29/LR/SPandx18arerestored/preservedexactlyasbefore. Wrapperframe528+map208=736belowtheexisting816stackceiling.

## Only two full coefficient tables,all entries constructed

ForC=A andC=B construct one vector tableentry(d*Cr*2^start/Rm modp,d*Ci*2^start/Rm modp),96bytes. Foraninput realword,xr,addfirstcomponenttoreal/secondtoimag. Forimaginaryword,xi,swapcomponents,SUBTRACTsecondfromrealandADDfirsttoimag. This is EXACTC*(xr+xi*u)withu²=-1. The sameA/Btablesareusedinbothoutputrowpairs. No negated/rotated duplicatecolumnsarestored,andnoentryisgivenbyanoracle.

Partition ALL381 significant inputbits intowidths[35,35,35,35,35,35,35,34,34,34,34],starts0,35,70,105,140,175,210,245,279,313,347. Eachcolumnhas96*(7*2^35+4*2^34)=27TiB;two54TiB. Allentriesforall11digits/bothrows/bothcoefficients areactuallyfilledbytheexisting75(n-1)+17<=96nfill_pair scheduleusingfull64-bitcounters,actualpinnedfullMontgomeryseedproductsandcanonicaladdition. Allzeroentries/seeds/coefficients/constants/flags/pagesareconstructedcold. Totalnewvector/fixedfillsremain<2^50instructions,safelyinsidecompletecold2^124.

Tablesmap8..62TiB. Move theunchanged13.5TiBfixed-s lookup to64..77.5TiB BY ITS EXISTINGbasepointerargument. Fixedhelper212instructions/tablecontents/canonicalfullsareunchanged. Workerparent4TiB,private<1TiB,childports80..128TiBremain. Allrangesaredisjointandstrictlybelow2^47;no48-bitnoncanonicalpointerorrelocatedchildportisneeded. Otherroutingnodesretainoriginalmaps. ExtraMMUpagesarechargedwithinexistingcoldsetup;bothworkers'tablesremain<2^48bytesandglobalpeak<2^138. NophysicalVMqualificationisclaimed.

## Nonnegative signed arithmetic and unchanged exact quotient correction

Eachoutputrealaccumulatorstarts atOLDthirdreal+22p. Theordinary7-limb22pconstantALREADYexistsatconstant3728 slotq22:offset528+64*22=1936. Actualmaploadsits7limbs andaddswithcarry,13instructionsperrowpairincludingpointerload/ADD. This is a computedconstant,notadvice. ThereareEXACTLYtwoimaginaryinputcolumns*11digits=22subtractedterms,eachin[0,p). Hence afterANYprefix,real>=seed+(22-k)p>=0 fork<=22. Atmost22positiveterms plusseed+22p give real<45p throughout. Imaginaryseed+44positiveterms likewise<45p. Every7-limbaddition/subtractionisexact;no448-bitwrap,negativeinteger,orborrowmagnitudeisignored. ThefinalSBC(withoutflagupdate)atlimb7isordinaryinstructionsemantics;nextimagADDSstartsanewcarrychain.

Because22p=0modp,thefinalresidueisEXACTtheoldrecurrenceoutput. Allintermediatesremainbelow45p<49p,soexistingD=floor(p/2^325),estimate/provedone-correctionQREDUCEandALL50ordinaryq*p slots remainvalidwithoutnewquotientoracle. Realandimagoutputeachcanonical0..p-1. MapstillseedsbotholdthirdcoordinatesBEFOREwritingthepair;otherpair'soldthirdcoordinatesremainuntoucheduntilread. Exactthird/outputaliasisretained;allotherinput/table/header/framepartialsareforbidden.

## Actual instructions and complete work

MSTEPusesMADD,threeLDP/sevenADDorSUBrealchain,threeLDP/sevenADDimagchain,optionaladvance. EachCOLUMNuses6inputLDR,11digitextracts,5EXTR/5MOVforcrossedlimbs,twoindex/baseinstructions,11MSTEPand10advances:270instructions instead283. A crossingloadsnextword,EXTRsintodigitx3,masks35/34bits,MOVsnextwordtox19beforetableloadsclobber20. Width35stride96*2^35=3TiB isoneMOVZx27,#0x300,lsl32;34-bitadvanceusesLSR1intheSAMEADD. Columnsize27TiBisoneMOVZx2,#0x1b00,lsl32. Everyimmediateisencodable;inputneverreadsbeyondsixlimbs. Fixedtablerelocationchangesavalue,notinstructionwidth.

Eightcolumns save8*(283-270)=104instructions;two13-instructionbiasinitializationsadd26. Map2461-104+26=2383. Fusedwrapperown306+map2383=2689throughRET. Caller16andpointerrotation4 give2709,retained2720. All8fullscalarinputs/4outputfields/rawhistories/quotientcorrectionsremainfullycharged. The32-instructionretainedsavingis22actualinstructionsplus10lessunusedslack:old2711+20=2731<=2752with21spare,new2689+20=2709<=2720with11spare. No operation is omitted or treated as an oracle.

Generation2720recurrence+824key+376NH+22aggregatescalar+704main=4646. Two localsortpairs720+scan128+cheap1+expensive1+spare256 give5752. Unchangednetwork933 giveshot6685,retained6760with75spare. SAMEunequalL=88722074495999233395607723401466085376,N=6760L+2^124 satisfiesEXACT N^125<20000^125*2^15562,proposed15562/125=124.496. Fullp381/r255,C20000/success0.39/sixtracks/onecomputationmetric,cold124/peak138/advice0andexactsuccess0.3900120659047603remainunchanged. Alltrials/failures/randomsource/preprocessing/transport/outputverificationpaid.

Shared-onlyauditinterpretsACTUALnew2689/2383sourceandfifthBOTHagainstindependentDENSEFp4oracleatfullp381/boundaries/randominputs. Readonlyvirtualtablefixturesuseexactnewtwo27TiBcolumn/11digitlayoutwithcoefficientsA/B;fixed-sfixturemoves64TiB. ThisdoesNOTreplaceproductiontablefillwithfreeoracle:allactualfillsandpagesarechargedabove. ActualSBCsemanticsisaddedtointerpreterwithoutchangingflags,withexistingno18/callee/LR/SP/write/read/aliasscopechecks. Othernative369hash/full-Y/U,full-rscalar,key/copypacketchecksremain. Thesearenotproductioninversion/tablefill/controllerhardware/GTprobabilitycertificates. No currentpass/qualification/score/promotionisinvented;candidatefileswereNEVERhostexecuted/imported/compiled. Sharedcredentialfreellm-verifierownsallactualcandidateexecution/review/scoring.

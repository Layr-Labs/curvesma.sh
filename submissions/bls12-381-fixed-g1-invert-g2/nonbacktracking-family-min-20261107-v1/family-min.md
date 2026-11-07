# Exact two-family minimum over the same twelve candidates

This replaces only the minimum-selection implementation. It returns exactly the same minimumFp2pairandappliedk as the fullreal-then-imaginaryunsignedorder. It introduces no failureevent,fieldarithmetic,hash/trajectorychange,probabilityhypothesisorfreeinterface. canonical-min.S iscomplete static AArch64evidence,neverhostexecuted.

The12imagesaretheevenfamily beta^m*x andoddfamily beta^m*y,m0..5,y=T(x),beta primitiveorder6inFp. Let a=real(x),c=real(y). Ifa!=0,thesixrealcoordinatesbeta^m*aaredistinct;likewise forc!=0. Ifbotharenonzero,findeachfamily'sminimumbycomparingonlyitsrealcoordinates. Eachfamilyhasuniquerealminimum. Theglobalminimumisoneofthesetwowinners;comparetheirreals,andonlyifequalcomparetheirimaginaries. Thisisexactlexicographicselection,nottheheuristicomissionofatiecase.

Ifa=0,allevenrealsarezero. Theoddrealsarenonzero:theFrobeniusconstantford1/5isnonzerodiagonal,soT(a+bu)hasrealpart +/-s*bwhen a0,withs!=0. Sincex!=0,b!=0,real(y)!=0. Hencezeroiseverywhereuniquelysmallestunsignedrealvalue,andtheminimumistheevenfamilymemberwithsmallestimaginarycoordinate. Thoseimaginariesarebeta^m*b,distinct. Similarly ifc=0,onlytheoddzero-real familycanwin;itsimaginaryvaluesaredistinct. Botha=c=0wouldforcex=0,whichisoutsidecallercontract. Noexceptionalvalueorprobabilitylossisignored.

Theassemblyimplementsallsuchbranches. Zero-realtestsuse3LDP,5ORR,1CBNZ=9each. Acomponentcomparisonuses3pairedloadsperoperand(6),oneSUBS/fiveSBCS(6),candidatepointer/ksetup2,andpointer/kCSEL2,total16. Eachfamily'sfivecomparisonscost80. Thetwo-winnerrealcomparisonuses12data/arithmeticinstructionsplus5nonflag-settingORRtoaccumulatedifferences. CBNZpreservescarry. Onnonzerodifference,finalCselectsthesmallerreal;onequalreal,theimaginary384-bitborrowchain12instructionssetsfinalC. TwoCSELselectpointer/k. Worstfinalcomparison32. No finalZbitisusedasmultiwordequality.

EntryfourMOV,normalzero-tests9+1pointer+9,familycomparisons160,saving/resettingwinner4,finalcomparison32,andoutput15give234internalinstructions. Theeven-zero path109,andodd-zero path121,arebothsmaller. Output15includessixLDP/sixSTPcopyingthefullwinning96bytes,STRk,MOVreturnk,RET. Allpairedimmediates0..80arealigned/inrange,andcandidateoffsets96k<=1056fitADD. Itusesonlycallerregisters,x18untouched,nointernalcall,frameorordinaryspill. Minimumoutputandkdestinationsaredisjointfromcandidatearena. Long-livedcontroller/programregistersremainpreserved.

Theexisting16caller/interfaceand14initialxcopyallowancesarestillcharged. Total234+16+14=264;retain272,eightslack. Thisreplacestheformer384minimum/copyallowancewithoutchanginganycandidate-generationproduct/rotation,fieldprogramcost,192key/controlor256wrappers. No per-Fp2memcmporinventedwideprimitiveisused.

Nonmultiplication1644=22*42fieldoperations+272minimum/copy+192key/control+256wrappers. WithnativeFpM384,worsthelper5484,orientationf(k)=3948+768c(k). ExistingB=floor1027K/768bounds
 3948K+768B<=4975K.
Canonicalcapfailure<2^-100isunchanged,asareallotherfailuresandH/bankedstoragemodelconditions. Fullordinary4975+807tabulation+797jointlabelhash+46scalar+849controller+192head+224payload+32cappedtraversals+1retry=7923. Retain7952K,29slack. ColdD=540*2^119,paddedphysicalpeak<2^136bytes,zeroadviceandcompleteverification/G2outputareunchanged. Atoneprotectedreferenceperstep,C20000,
 W<=(27952K+540*2^119)/20000.
Exactintegercertificate(27952K+540*2^119)^250<20000^250*2^31541 establishesproposed31541/250=126.164. Thisisanalyticconditionalevidence,notofficialscore,qualification,productionexecutionorpracticalhardwareavailability.

# Full logical orbit key and fingerprint in one live-frame call

The full192-byte canonical orbit key VALUE is unchanged. It is represented by two active96-byte regions A96 and B96,read in that order,not by a newly written contiguous192-byte scratch. All four canonical full-p381coordinates/all24LE64words remain. This is NOT truncation,a partialtraceprojection,newfingerprintor new collision assumption.

## Lifetime, exact words and native API

canonical_twelve_bfirst nowtakesx0payload88,x1readonlyrawPhi192,x2fixed-sLUT,x3p48,x4constant3728,x5hashkeys672. Its existing384-byteframe saves19..29/LR and hasEvenB96at96/OddPhi192at192. Savex5keys inexistingcallee28withONE MOV;28wasalreadysavedandotherwiseunused. Sign/odd/fixed helpers preserve28. Alloriginalcanonicalization/comparison/selectedA0andB1pointerlogic is unchanged. Identityreturns1BEFOREhash,payloadU/taginvalid/scalaruntouched. Nonidentityreturns0withU32/tag24initialized;fullFr32atpayload32..63isNEVERwritten.

Atselectedpath,26pointsatA96(rawa0orOddPhi.a0),27atselectedB96(EvenBorOddPhi.b1). FourMOV+BLpasspayload/26/keys28/27 toactualnh_fingerprintWHILEthe384frameisSTILLallocated. No danglingpointereverescapes. ThehashreadsALL24wordsbeforeanyoutputstore. OnlyTHENkeyrestoresframe/callees/LR/SP. Rawhistories/tables/keys/header/frame/outputaremutuallydisjointexceptdocumentedcanonicalB in-placeupdatesinsidekeyframe. Hashoutputisdisjointfromraw/keyframeinproduction. Allinitializedread/write/ABI/exactsignandT²=-1contractsremain.

nh_fingerprint nowtakesx0payload88,x1A96,x2keys672,x3B96. FirstsixLDPinputpairsreadx1offset0..80;lastsixreadx3offset0..80. InterleavedkeylayoutandALLwrappedNH128/chunks/affine/fullY/low160operationsareunchanged. DelaytheONEqmodulusMOVintoreg3untilafterALL12NHrounds,whenBpointerisdead. Nopointervalueispartofthehash,onlytheSAME192logicalbytes. SourceisSTILL369throughRET;directcaller5gives374<=376,butproductioncallerisnowWITHINthekeyboundandmustNOTbechargedagain. Bothworkersusethesameexactrepresentation/function;alltrue-key/U/tag equality and delta/eta proofsareidentical.

## Complete own-path and all-call bound

Theoldsuccessfulkeyownbound111included24LDP/STPcopyinstructions. Removeall24,addONEkeysMOVandfourargumentMOV+BL5:own<=111-24+1+5=93. A separateagent-ownedACYCLICSTATICCFGcheckgives92ownand1112withknowncalleeupperbounds;retain93/1113conservatively. ItdoesnotexecuteCPUinstructions/registers/flags. EachBLcountsoneowninstructionanditscallee'sseparateupperbound.

Two65signhelpers+odd521+NH369 give1020calleeinstructions. Therefore fusedsuccessfulpublickey/hash<=93+1020=1113throughRET. TheidentitybranchdoesNOhash;own<=86+two65+odd521<=737,alsobelow1113. Thischecks thepaththatdidnotformerlycopy24ratherthanillegallysubtracting24there. Sixargumentsetup+BL+identitybranch8,retain9,so fullkeyhash<=1122,retained1125. Theold824+376=1200retainedsumisREPLACED,notaddedto1125. Actualworstcontrol/copyupperboundimproves21;threepreviouslyunusedallowanceinstructionsgive24retainedsaving. Nooperationisforgotten.

## Exhaustive sign-helper path partition

For canonical native Br,Bi, the sign decision uses the first nonzero component. If Br!=0, its initial CBNZ skips the three imaginary loads. The decision prefix is therefore at most21, not the earlier24 combined-path bound. If sign is high, loading p costs5 and each nonzero component costs19 including its zero-test branch and continuation; RET1 gives at most21+5+19+19+1=65. If either component is zero its path is at most12 rather than19, so this cannot increase the bound. The low-sign copy costs at most21+12+1=34.

If Br=0, the prefix including imaginary loads costs24. For a nonzero high Bi, the first component is known zero, so the high path is at most24+5+12+19+1=61. Low sign costs at most24+12+1=37. When both are zero, the chosen component is0<ceil(p/2), so only that low-copy path is possible and it copies canonical zeros. These cases partition ALL canonical full-p Fp2 inputs. Thus65 is a deterministic worst-path bound throughRET, preserving the existing sign convention and exact-alias/ABI contract. The earlier68 was conservative by combining incompatible branches.

## Shared-only exact input trace,not a digest-only check

The sharedinterpretercapturesACTUALldp x10,x11inputvaluesexecutedinsideNH. Independentfull-podd-Frobenius/sign/12-wordorderoracleconstructstheentireexpectedS192andcomparesALL24readwordsEXACTLY. Hashoutputequalityalonecannotmaskawrongselectedkey. ItseparatelycomparesallthreefullYwordsandU32 toindependentintegerNH/affineoracle,checks369calleeinstructions,fused<=1113,allpublic19..29/LR/SP/x18restored,readonlyrawPhi,allowedwriteregionsandfullscalar32untouched. DirecthashfixturesusegenuinelyNONCONTIGUOUSA/Bregions,plus documentedA/payloadbasealiascases,samewrappedNH/affinecornerinputs. Everyregionlifetime/initializedreadischeckedwithinactiveframe.

IdentityfixtureisactualGTidentityPhi=(Rm,0,0,0),returns1andperformsNOhash/noreadinvalidU/tag. Othernonzero-b1fixturesarearbitrarycanonicalfull-pfieldcases,NOTclaimsofGTsamples/probabilitymeasurements. Nestedstackremains384+336+96=816;NHhasnoSPframe. Fullkeyfiber/successproofisanalytic,notempirical. Candidatefileswereneverhostexecuted/imported/compiled/emulated;onlysharedcredentialfreellm-verifierexecutes/reviews/scores. Nonewpass/score/acceptance/qualification/promotionisinvented.

## Complete current resource bound

Nativegeneration2417recurrence+1125fusedkeyhash+20aggregatescalar+704main=4266. Sort648+scan128+cheap1+expensive1+spare256=5300;network906=595generation+275sort+36scan,giveshot6206retained6284with78spare. SAMEfullL88722074495999233395607723401466085376,N=6284L+2^124 satisfiesEXACTN^1000<20000^1000*2^124391,proposed124391/1000=124.391. Fullp381/r255,uniformGT/noauxpowers,success0.39/C20000/sixtracks/onecomputationmetric,cold124/peak138/advice0andexact0.3900120659047603successremainunchanged. Allnotifications/earlyheader/stablecontinuations/failedtrials/caps/input/outputworkremainfullychargedbytheirunchangedexplicitconstruction.


CURRENTsortisstablelow64DIRECTthenhigh96RELAYED,not80/80. Historicalradix80_certificate.mdnowcontainscompletecurrentroles/arrayranges/phasebarriers/digitprefix/queue/finish/ordinarynotification/trafficproofs. Localsort648/network275;gen4266/base5300/network906/hot6206ret6284+78spare. SORTCOLD120replacesold104;totalcold124/peak138/advice0unchanged. Nochangeinfull160U/128ordinals/key/hash/success.



CURRENTtenfull39/38digitsretainALL381bits;LUTvectors128..656TiB/fixed672..804,childports1024..1280,64KiBedgebackingwith12KiBactivehalves/NormalShareable/FEAT_LVA+LPA/64KiBgranule/52-bitlocalMMUfullyconstructedinlarge_address_certificate.md. Edgecold121(memory113)/processor121/sort120/other92sumcold124,peak138/advice0. Actualrec2392ret2417/fixed194/odd521/fused1113ret1125/gen4266/base5300/net906/hot6206ret6284spare78. NoVMhardware/productionqualificationclaimed.

# Full logical orbit key and fingerprint in one live-frame call

The full192-byte canonical orbit key VALUE is unchanged. It is represented by two active96-byte regions A96 and B96,read in that order,not by a newly written contiguous192-byte scratch. All four canonical full-p381coordinates/all24LE64words remain. This is NOT truncation,a partialtraceprojection,newfingerprintor new collision assumption.

## Lifetime, exact words and native API

canonical_twelve_bfirst nowtakesx0payload88,x1readonlyrawPhi192,x2fixed-sLUT,x3p48,x4constant3728,x5hashkeys672. Its existing384-byteframe saves19..29/LR and hasEvenB96at96/OddPhi192at192. Savex5keys inexistingcallee28withONE MOV;28wasalreadysavedandotherwiseunused. Sign/odd/fixed helpers preserve28. Alloriginalcanonicalization/comparison/selectedA0andB1pointerlogic is unchanged. Identityreturns1BEFOREhash,payloadU/taginvalid/scalaruntouched. Nonidentityreturns0withU32/tag24initialized;fullFr32atpayload32..63isNEVERwritten.

Atselectedpath,26pointsatA96(rawa0orOddPhi.a0),27atselectedB96(EvenBorOddPhi.b1). FourMOV+BLpasspayload/26/keys28/27 toactualnh_fingerprintWHILEthe384frameisSTILLallocated. No danglingpointereverescapes. ThehashreadsALL24wordsbeforeanyoutputstore. OnlyTHENkeyrestoresframe/callees/LR/SP. Rawhistories/tables/keys/header/frame/outputaremutuallydisjointexceptdocumentedcanonicalB in-placeupdatesinsidekeyframe. Hashoutputisdisjointfromraw/keyframeinproduction. Allinitializedread/write/ABI/exactsignandT²=-1contractsremain.

nh_fingerprint nowtakesx0payload88,x1A96,x2keys672,x3B96. FirstsixLDPinputpairsreadx1offset0..80;lastsixreadx3offset0..80. InterleavedkeylayoutandALLwrappedNH128/chunks/affine/fullY/low160operationsareunchanged. DelaytheONEqmodulusMOVintoreg3untilafterALL12NHrounds,whenBpointerisdead. Nopointervalueispartofthehash,onlytheSAME192logicalbytes. SourceisSTILL369throughRET;directcaller5gives374<=376,butproductioncallerisnowWITHINthekeyboundandmustNOTbechargedagain. Bothworkersusethesameexactrepresentation/function;alltrue-key/U/tag equality and delta/eta proofsareidentical.

## Complete own-path and all-call bound

Theoldsuccessfulkeyownbound111included24LDP/STPcopyinstructions. Removeall24,addONEkeysMOVandfourargumentMOV+BL5:own<=111-24+1+5=93. A separateagent-ownedACYCLICSTATICCFGcheckgives92ownand1163withknowncalleeupperbounds;retain93/1164conservatively. ItdoesnotexecuteCPUinstructions/registers/flags. EachBLcountsoneowninstructionanditscallee'sseparateupperbound.

Two65signhelpers+odd572+NH369 give1071calleeinstructions. Therefore fusedsuccessfulpublickey/hash<=93+1071=1164throughRET. TheidentitybranchdoesNOhash;own<=86+two65+odd572<=788,alsobelow1164. Thischecks thepaththatdidnotformerlycopy24ratherthanillegallysubtracting24there. Sixargumentsetup+BL+identitybranch8,retain9,so fullkeyhash<=1173,retained1176. Theold824+376=1200retainedsumisREPLACED,notaddedto1176. Actualworstcontrol/copyupperboundimproves21;threepreviouslyunusedallowanceinstructionsgive24retainedsaving. Nooperationisforgotten.

## Exhaustive sign-helper path partition

For canonical native Br,Bi, the sign decision uses the first nonzero component. If Br!=0, its initial CBNZ skips the three imaginary loads. The decision prefix is therefore at most21, not the earlier24 combined-path bound. If sign is high, loading p costs5 and each nonzero component costs19 including its zero-test branch and continuation; RET1 gives at most21+5+19+19+1=65. If either component is zero its path is at most12 rather than19, so this cannot increase the bound. The low-sign copy costs at most21+12+1=34.

If Br=0, the prefix including imaginary loads costs24. For a nonzero high Bi, the first component is known zero, so the high path is at most24+5+12+19+1=61. Low sign costs at most24+12+1=37. When both are zero, the chosen component is0<ceil(p/2), so only that low-copy path is possible and it copies canonical zeros. These cases partition ALL canonical full-p Fp2 inputs. Thus65 is a deterministic worst-path bound throughRET, preserving the existing sign convention and exact-alias/ABI contract. The earlier68 was conservative by combining incompatible branches.

## Shared-only exact input trace,not a digest-only check

The sharedinterpretercapturesACTUALldp x10,x11inputvaluesexecutedinsideNH. Independentfull-podd-Frobenius/sign/12-wordorderoracleconstructstheentireexpectedS192andcomparesALL24readwordsEXACTLY. Hashoutputequalityalonecannotmaskawrongselectedkey. ItseparatelycomparesallthreefullYwordsandU32 toindependentintegerNH/affineoracle,checks369calleeinstructions,fused<=1164,allpublic19..29/LR/SP/x18restored,readonlyrawPhi,allowedwriteregionsandfullscalar32untouched. DirecthashfixturesusegenuinelyNONCONTIGUOUSA/Bregions,plus documentedA/payloadbasealiascases,samewrappedNH/affinecornerinputs. Everyregionlifetime/initializedreadischeckedwithinactiveframe.

IdentityfixtureisactualGTidentityPhi=(Rm,0,0,0),returns1andperformsNOhash/noreadinvalidU/tag. Othernonzero-b1fixturesarearbitrarycanonicalfull-pfieldcases,NOTclaimsofGTsamples/probabilitymeasurements. Nestedstackremains384+336+96=816;NHhasnoSPframe. Fullkeyfiber/successproofisanalytic,notempirical. Candidatefileswereneverhostexecuted/imported/compiled/emulated;onlysharedcredentialfreellm-verifierexecutes/reviews/scores. Nonewpass/score/acceptance/qualification/promotionisinvented.

## Complete current resource bound

Nativegeneration2664recurrence+1176fusedkeyhash+22aggregatescalar+704main=4566. Sort720+scan128+cheap1+expensive1+spare256=5672;network933=595generation+302sort+36scan,giveshot6605retained6680with75spare. SAMEfullL88722074495999233395607723401466085376,N=6680L+2^124 satisfiesEXACTN^1000<20000^1000*2^124479,proposed124479/1000=124.479. Fullp381/r255,uniformGT/noauxpowers,success0.39/C20000/sixtracks/onecomputationmetric,cold124/peak138/advice0andexact0.3900120659047603successremainunchanged. Allnotifications/earlyheader/stablecontinuations/failedtrials/caps/input/outputworkremainfullychargedbytheirunchangedexplicitconstruction.

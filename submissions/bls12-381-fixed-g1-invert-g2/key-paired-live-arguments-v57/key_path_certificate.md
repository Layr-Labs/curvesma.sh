# Current complete native key paths

Only the internal exact sign normal form and ordinary helper implementations changed from the source;the full192-byte key still has exact twelve-Frobenius fibers. key_contract.md gives the all-input normal-form/fiber theorem. S is now transient worker scratch,not a stored payload,so production correctness uses actualGT exponent verification as in proof.md.

## Exhaustive scalar sign paths

The original prefix loads Br and ORs its six limbs into x14;only if zero does it load Bi into x2..x7. It retains this first-nonzero component and zero flag. Load p once,compute the low-word p-a by64-bit SUB,compare unsigned with a_low,and copy or canonically negate BOTH components. p_low odd prevents a nonzero tie;bothzero takes the canonical-zero copy. All chosen source limbs are loaded before their stores;exact alias/disjoint contracts and all callee/SP/LR preservation hold.

Br nonzero prefix17,zero-flag branch1,real negate/store9,imag load/test9,nonzero imag negate/store/RET10 totals46. If Bi zero,the final block is4,so40. Br zero/high Bi is20+1+3 zero real stores+6 subtract+3 stores+RET=34. Low-copy paths<=33. These cases partition ALL canonical Fp2,not favorable random paths. Thus46 throughRET is a safe deterministic bound. The actual source makes no half-modulus decision,though the unused public constant remains harmless data.

## Odd transform and complete parent

The odd helper saves/restores only its used x19..x22,x29/LR and maps arguments directly into working registers. Its frame/scratch layout remains336 and raw input immutable. Direct canonical negation loads all six limbs,OR-detects zero,stores zero or subtracts from full p. Its nonzero21/zero13 paths replace zero-buffer/general subtraction;no zero becomes p. The two remaining add/sub helpers each31 and two constructed fixed-s lookups each192 are unchanged. Own longest58 plus62+384=504 throughRET (zero paths shorter). Every called helper preserves the ordinary ABI.

The current parent retains one CMP+B.NE per limb. The first differing comparison reaches a shared pointer chooser with unchanged NZCV. ADD x6,SP,#192 prepares odd a0 without changing flags. CSEL.LO chooses raw even a0 versus odd a0 into x6,and even versus odd signed b1 into x7. Neither CSEL changes flags. The copy uses x6/x7 and x0/x1 loads,so all live pointers remain intact. Equal limbs continue;all equal returns unchanged identity status. This preserves exact unsigned lexicographic comparison/output/ABI/frame384 contracts. All-path acyclic CFG traversal counts every conditional outcome,BL one own instruction and RET,and rejects cycles. Own bound89,retained89,one saved versus previous shared-branch chooser95/96. Complete key89+88+504=681,eight callers689,retained692 with3 spare. Every nonidentity output byte initializes and rawPhi remains immutable.


Current independent shared-only interpreter checks actual constants,all conditional carry/zero operations,initialized reads,permitted writes,immutable rawPhi,ABI,immediates,sign<=44,odd<=504,fixed192,key<=681 against its independent full-p formula,including identity and field boundary cases. Its scoped nonzero field fixtures are not asserted to be actualGT samples. They do not replace the twelve-fiber theorem or certify production success,cost,table construction,hardware or qualification. Full recurrence2465/map2386 and full183-bit NH341 checks remain separate. Global work is only the current proof/architecture ledger;older benchmark totals are not claims of this certificate.


CURRENT streamed-inputfixedhelper192:keepx1astheactualinputpointerthroughALLtwelve lookups. UBFXextractsdigitintox2,andMADDx2,x2,x28,x26 computesitsentryaddressbeforeload,sox1 isnotclobberedbydigits. Loadtwoinputlimbsintox3/x4 atoffset0,consumeallfourdigits,thenloadnextpairat16,thenlastpairat32. ThethreeactualLDP countisunchanged. AllsixfullinputlimbsarestillreadBEFOREanyoutputstore,soexactinput/outputalias remainsvalid. No retainedinputlimb isoverwrittenbeforeitslastuse;topdigitretains29bitsandall381fieldbits.

Calleex21..24areneverreadorwritten;onlyx25..28 areused. Frame32 saves/restoresthemwithtwoSTP/twoLDP instead ofeightcallees/fourpairs. No newpointerMOVisneededbecausex1 stayslive. Fouractualinstructionsdisappear;allproduct-base/reducer28/table/field/ABIsemanticsremainunchanged. Completefixed192=148lookupcore+44remaininginstructions,odd504=58+62+384,key681=89+88+504,eightcallerwork689,retained692. Nativefusedleaf2465/retained2506 staysunchanged. Non-scalar4260generation/6365total plus12m+28n+cold2^124. No currentcounter,transportorMain712 allowanceischanged.


CURRENT keyparentstackarguments:tablex2/p48x3/productsbasex4areusedonlyfortheoneoddhelpercall. Storex2/x3atframe32/40andx4at48 BEFOREfirstsigncall,thenreloadwithoneLDP/oneLDRjustbeforeoddcall. Callees21..23areNEVERread/writtenbyparent;calledhelperspreservethem. No threeprotectedargumentMOVs,threeoddargumentMOVs,or twocalleerestoresareneeded:twoactualreloadsinsteadofthreeMOV save1,threeinitialMOV disappear,andtworestoreinstructionsdisappear,total6. Initialstorecountisunchangedbecauseargumentstoresreplacethecallee-backupstores.

Only19/20and29/LRaresaved/restored;frame384/evenB96/oddPhi192layoutsandalloutput/identity/compare/ABIcontractsareunchanged. Stackargs32..55aredisjointfromevenB96andoddPhi192,andchildframesgrowbelowtheparent. Allargsareinitializedbeforeanyhelpercall;rawPhi staysimmutable. All-pathacyclicCFG bound88own,retain89,fullkey89+88+504=681,eightcallerwork689,retained692. Actualfield2465/fixed192/hash341unchanged. Generation4260L/non-scalar6365L+12m+28n+cold2^124;no probability/field/routing/Main712 change.


CURRENT oddhelperpairedarguments:tablex2andproductsbasex4 arestoredbyoneSTPatframe152/160,replacingtheformerconstSTRwithoutanextrastore. Callee22isNEVERread/written;removeitsoneinitialMOV. Save/restoreonlycallee21withSTR/LDRat32instead of21/22pair,sameinstructioncount. Eachoftwofixedcalls reloadstable/productsbasewithONE LDPx2/x3at152 instead ofMOVtable+LDRconst,twofurtherinstructionssaved. Total3actualinstructions;allpointerreadsareinitialized,pair152isencodable8-bytealigned,andargsaredisjointfromscratch192/240withinframe336. p48stayscallee21fornegation/add/sub.

Oddown58+two31Fphelpers+two192fixedcalls=504. Parentretain89/sign2*44givekey681,eightcallerwork689,retained692. Fused2465/fixed192/hash341,allrawPhi/output/ABI/alias/zero andfullfield/GTorbitpropertiesareunchanged. Generation4260/non-scalar6365L+12m+28n+cold2^124. No main712,probability,transport,bucket/table/coin/verificationpricechange.


CURRENT signhelperlivep48:canonical_conjugate_b1receivesthealreadyconstructedactualmodulus48inordinaryargumentx3(x2unusedcompatibility),inputx0/outputx1unchanged. Itssecondinputlimbusescallerx17ratherthanx3,sothereal/imagLDP/STP andallcarry/zerologicretainactualfullfieldvalueswhilex3holdspthroughthreepLDP. RemoveADRP+ADDofstaticp,2instructionspercall;allpwordsareactuallyreadfromtheconstructedmodulusandverifiedbythefull-pfixture. No p-derivedanswer,auxiliarypower,rodataexemptionorfieldchangeisintroduced. Originalstaticp/halfceilremaininpaidcode/rodataforreferencechecks.

Parententryx3alreadyISthisp48andisstoredatframe40withoutclobberbeforefirstsigncall,sofirstcallrequiresNOnewsetup. Beforesecondsigncall,oneexplicitLDRx3,[SP40]restoresitafteroddhelper. Parentownincreases1from88to89(retain89),eachsignfalls2from46to44,totalkeysaves3. Odd504/fixed192/field2465/hash341andallorbit/ABI/alias/outputcontractsareunchanged. Key90+2*44+504=681,eightcallerwork689,retained692;generation4260/non-scalar6365L+12m+28n+cold2^124. Main712andallprotectedparameters/prices/probabilities/transportremainunchanged.


CURRENT parentpairedlivearguments:storetable/productsbasewithSTPx2/x4atframe32/40,andp48withSTRx3at48,thesame2initialstoresasbefore. Firstsigncall'sactualcanonical_b1.S neverwritesitsinputp48registerx3andusesnoBL,sox3remainsvalidforoddcall. Reloadtable/productsbasewithONE LDPx2/x4at32insteadoftable/pLDPplusconstLDR. Thisremovesoneactualinstruction;no omittedp setupisassumed. Afteroddcall,second-signp48reloadisexplicitLDRx3at48. Allargsinitialized/disjoint,frame384/ABI/rawPhi/outputandfieldtheoremsunchanged.

Parentown88retain89,sign44twice/odd504givekey681,eightcallerwork689,retained692. Nativefield2465/fixed192/hash341unchanged. Generation4260/non-scalar6365L+12m+28n+cold2^124. Thep48wasalreadyconstructed/delivered;allactualcallerloadsarecounted. No probability/main712/transport/table/inputdistributionchange.

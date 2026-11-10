# Exact small-quotient reduction of both full-width accumulator types

Write B=2^325,D=floor(p/B)=58556206801485620=0xd0088f51cbff34>256, p=DB+e with0<=e<B. For every0<=Y<256p let Y=qp+t with0<=t<p,q<=255. Then U=floor(Y/B)=qD+floor((qe+t)/B). The extra is nonnegative and below D+256<2D, so qhat=floor(U/D) is q or q+1. Since256p<2^389, U fits64 bits and the high-limb EXTR#5 computes it exactly. UDIV is one ordinary64-bit scalar instruction priced1 by the protected metric; this is not a CPU-cycle or arbitrary-precision-division assertion.

For the actual recurrence map, each output is the sum of46 canonical digit-table terms and one canonical third-history seed. Thus0<=Y<49p, q<=48 and qhat<=49. For the fixed-s lookup there are twelve canonical terms,givingY<12p andqhat<=12. These are deterministic bounds for every possible input, not sampled constraints.

## Fully constructed3728-byte constants

Retain the prior520-byte header, including ordinary p at448 and D at512. Initialize padding520..527 to zero. Construct FIFTY64-byte slots at528+64j,0<=j<=49. Slot j contains seven ordinary64-bit limbs of j*p plus an initialized zero eighth word. Build slot0 as all zero; each subsequent slot is the exact seven-limb integer sum of the preceding slot and p, with no modular reduction. Because49p<2^448, no significant carry is lost. All fifty slots, header and padding are actually filled and delivered to both workers in this attack; their scalar addition/store/transport cost is included cold. No searched advice, given powers or field-multiplication oracle is supplied.

The same3728-byte block is passed through the existing arguments to the map and key/fixed-s helper. The recurrence reads0..49; fixed-s reads only0..12. All reads are in range and ordinary j*p is not mistaken for a canonical Montgomery field value. Increasing the block from1296 to3728 changes no hot argument count, pointer width, record size or storage protocol.

## Actual fused recurrence QREDUCE28

The actual fused recurrence_fp2_step.S usesa28-instruction address-selected reducer. Standalone recurrence_fp2_map.S remains inertreferencewithits31-instruction addbackreducer,nevercalledbytheactualhotleaf. Fixed-s nowusesits29-instruction selected-address reducer. Allfiftyordinaryq*p slots/header/paddingremainfullyconstructed/paid.

EXTR/UDIV computesqhat=q orq+1. Thefirstpaidaddressiscalculated asconst+528+64qhat,andthesevenordinaryproductlimbsareloadedintox11..17. TrialSUBS/SBCS writesintox11..17 whilepreservingtheORIGINAL Y inx4..10. Itsfinal448-bitborrowistrueEXACTLYwhenqhat=q+1. NowSUBx19,x26,64,WITHOUTsettingflags,formsanaddressone64-byteslotearlier. CSELx26,x19,x26,LO selects itONLYonthatborrow. Thereforex26pointsEXACTLYtoslotq=floor(Y/p),with0<=q<=48. Thisusesordinarypointerarithmeticandtwoexplicitinstructions,not freeindexing oraconditionalreadoracle.

Whenqhat0,borrowisfalse. Theunusedtemporaryaddressisconst+464,insidetheconstructedheader,butisNEVERdereferenced. Thereisnoout-of-boundsaccessorselectednegativeindex. Thesafechosenaddressisalwaysconst+528+64q. Reloaditssevenordinaryq*p limbswiththreeLDP/oneLDR,andfinalSUBS/SBCS into originalx4..10givescanonicalY-q*p in[0,p),zeroseventhlimb. Firsttrialresultregistersarediscarded;originalY isneverlost. x19isdeadbetweencolumnlookupsandreductions,alreadysaved/restoredbytheleaf. Noadditionalcallee/argument/framewordisintroduced.

Exactcount:EXTR/UDIV2,initialaddress2,initialload4,trialsubtract7,nonflagaddressSUB/CSEL2,reload4,finalsubtract7,total28. Thisremoves8instructionsacrossfourreducersrelativetothe30-instructionquotient-reselectionleaf,and12relativetotheearlier31-instructionaddbackleaf. Neitherrecomputationofaselectedquotientnorasecondbase/indexaddresscalculationisneeded. x2remainsD,x28const,x26selectedaddress. RowBissavedbeforethefirstreductionandreloadedasitsownoriginalaccumulator. Everyactualread/alias/ABI/borrow/canonicalresultandinstructioncountischeckedbytheindependentfull-pnativefixture.

Allconstantconstruction/layout/tablevaluesandthefullsizeone-highquotientproofaboveareunchanged. Coldwork/pointerwidth/headerdelivery remainpaid. Fullgloballedgersareinproof.md/architecture.md/radix80_certificate.md;finitehelperchecksarenotwholeattackcertification.


## Fixed-s selected-address reducer29

The fixed_s_lookup.S helper nowusesthesame provedpointercorrection withoneextraDload. OriginalY isx5..11,product/trialresultx12..17/x26,selectedslotaddressx2,constx25,andunusedtemporaryaddressx1. ThefullY<12p boundgivesqhat<=12 andq=floor(Y/p)<=11. SUBx1,x2,64doesnotsetflags;CSELx2,x1,x2,LO selectsthepreceding slotONLYontrialborrow. Correctseven-limbordinaryq*p isreloaded,andfinalsubtractionintooriginalY iscanonical. It costs29,3lessthanolder32;allcallee/alias/read contractsandthefull3728-byteconstructedblockareunchanged. Theactualstandalonemap31 staysuncalledreference;actualfused28 isunchanged.

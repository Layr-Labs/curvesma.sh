# Current exact sign normal form and twelve-orbit key

For canonical B=(Br,Bi),choose the first nonzero component a,or zero if bothzero. For a!=0,the low64-bit limb of p-a is(p_low-a_low) mod2^64. Since p_low is odd,equality would require2*a_low=p_low mod2^64,impossible. Choose the member of{B,-B} whose selected component has the smaller low limb;set N(0)=0. Negation preserves the first nonzero position and reverses this comparison,so N(B)=N(-B),and N(B)=N(C) iff C=±B. The native helper retains the selected limbs/real-zero flag and performs canonical zero/p-a output with exact aliases. Exhaustive high paths46/40/34 and low copies<=33 give the conservative46 bound;key_path_certificate.md gives its complete ledger. No half-modulus or probabilistic sign rule is used in the current algorithm.

# Exact key and retained native helper contracts

For nonidentity GT,X,X^p4,X^p8 have product1,trace c=3(a0+b1t),second symmetric coefficient c^p2;they are roots of Y³-cY²+c^p2Y-1. Thus trace identifies exactly3 Frobenius conjugates. p² fixes a0,negates b1. Canonicalize b1,-b1 by the current first-nonzero low-word involution below;negate BOTH components if needed. F=(a0,canon_pm b1) has6 even-Frobenius fibers. b1=0 implies c inFp2 and cubic factors(Y-1)(Y²+(1-c)Y+1);a nonidentity root would lieFp4,impossible since embedding degree12. Identity abort is therefore exact for GT inputs.

Cold compute actual C0=xi^((p-1)/2)=(s,s), 2s²=-1. Odd Frobenius conjugates a0 and maps B=(Br,Bi) to T(B)=(s(Br+Bi),s(Br-Bi)). Thus T²(B)=-B. Let q(B) be the retained canonical sign representative of {B,-B}. If q(B)=q(T(B)), then T(B)=B or -B, so T²(B)=B; since p is odd, B=0. Nonidentity GT has B!=0 by the preceding cubic argument. Therefore these two signed b1 representatives are ALWAYS DISTINCT on nonidentity GT.

Use the fixed total order that first compares q(B), using unsigned lexicographic order of its twelve LE64 words in increasing word index, then compares a0. The latter tie-break is unnecessary for nonidentity GT. Choose the smaller of (a0,q(B)) and (conjugate(a0),q(T(B))) under this order. Store the complete selected pair, still a0 first then b1 in192 bytes. This retains all four full Fp coordinates. Applying p swaps the candidates, up to signs, while p² changes only sign B and p4 preserves Phi. Thus this key is constant on each12-element Frobenius orbit. Conversely identical selected full pairs imply equal selected traces; the cubic then gives one of the three p4-conjugates. Combined with sign and odd selection this implies the same12-element orbit. There are exactly12 fibers. No factor3 multiplication is needed because3 is nonzero. Identity has B=0 and both compared candidates zero; the helper returns failure before output reads. All inputs here are actual GT traces from the proved recurrence, not arbitrary Fp4 points.

## Current native contracts and ephemeral storage


## Current actual native/source API

The full ordered192-byte canonical vector above is conceptual:the actualleaf reads everyword from rawAreal/selectedAimag/signedB and hashesit,withoutphysicalS output. It returns ALL183 fingerprintbits word0=x1,word1=x0,word2=x5(low55). Compatibilityarg0unused/rawPhi1/sLUT2/p483/products4/keys5. NOoutputmemoryaccess. Identityx7zero invalidatesresults;successx7 is the nonzero selectedBpointer. CallerMUST useFULL64-bitCBZx7 and never dereferenceit afterreturn.

The current complete all-path/native/source proof is key_path_certificate.md:ONE576 frame withlowerinner0..287/header288/EvenB336/OddAimag432/OddB480;deadoddrealAisneveraccessed. Raw858 retains879. Worker19out/20color/21-22fullordinalareABI-preserved. ActualsourceCBZx7 skipsidentity;success executesSTPlo128/ORRcolor55/STRhi16/STPordinal24. Fourdatainstructions,ordinalSTPalreadypaid;newsource3 plusmasks/spare givesMain720. The integratedisolatedoracle givesleafZEROoutputwritebudget,checksfullregisterhash/exactABI/status,andexecutesactualsourcecompletionofall40bytes. Every nonidentityfixtureforcesstatuspointer2^32(low32zero),checkingactual64-bitCBZ.

All183hashbits/color/pending0/full128ordinalremain inmetadata40;maskhighword55beforehigh96/equality. Fullfield/order/GTfiber/mathabove,andexactgroup-event/probability/coin/cap/refprice/cold/storageproofsremain. Field2444+key879+Main720=4043,network886 givesN6114L+12m+28n+cold2^124. ScopedfixturesdonotcertifyproductionGTsuccess/controllerhardware/theorems;sharedreviewisrequired.

Currentarg5is672-byteNH/affine-pointerblockforNEWfour32-chunkfamily,15independentFqcoefficients. Full12tablesof2^32entriesperworkerareconstructed/paidinproduction,atworker-private[1,1.375)TiB;VMfixturesonlyprefix/virtualvalues. CollisionboundNOW7/qperunequalNHrow/delta=(2^-63+7/q)^3,all183outputbitsremain,butold3chunkhashbytesarenotclaimedidentical. Completeproofinnh_certificate.md;GTcanonical/group/field/metadata/sourcecontractsunchanged.


CURRENT internal worker ABI23=LUT/24=products/25=D/26=context/27=48/28=96/29=actualp48. Paid initializer copies actual pointers from cold args6..10 in EIGHT native instructions. All setup/delivery<=2^16 added inside cold2^124 outside R; no entropy/advice/powers. Helpers preserve globals at ABI boundaries. Key saves actual LUT/products atSP48 before paired-sum clobber, restores keys/context before fixed bodies, uses saved slots and live25 instead of obsolete arg spills. Field preserves products24/size28 and restores29 around temporary SP use. Crypto, full sizes, success and prices unchanged.

Currentcold25=actualDloadedfrominitializedproductheader;LOCALkeys672pointerlivesininitializedctx26+104. NHfirstLDR2,[26,104]SAMEoneoperationasoldMOV2,25;paired25/26preservesD/contextbeforefixedbodies. Allsource/hash/full64status/ABIunchanged. Actualcoldsetup/deliverypaidin2^16inside2^124;nativeinput672unchanged.

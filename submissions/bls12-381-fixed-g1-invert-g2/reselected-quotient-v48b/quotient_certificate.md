# Exact small-quotient reduction of both full-width accumulator types

Write B=2^325,D=floor(p/B)=58556206801485620=0xd0088f51cbff34>256, p=DB+e with0<=e<B. For every0<=Y<256p let Y=qp+t with0<=t<p,q<=255. Then U=floor(Y/B)=qD+floor((qe+t)/B). The extra is nonnegative and below D+256<2D, so qhat=floor(U/D) is q or q+1. Since256p<2^389, U fits64 bits and the high-limb EXTR#5 computes it exactly. UDIV is one ordinary64-bit scalar instruction priced1 by the protected metric; this is not a CPU-cycle or arbitrary-precision-division assertion.

For the actual recurrence map, each output is the sum of46 canonical digit-table terms and one canonical third-history seed. Thus0<=Y<49p, q<=48 and qhat<=49. For the fixed-s lookup there are twelve canonical terms,givingY<12p andqhat<=12. These are deterministic bounds for every possible input, not sampled constraints.

## Fully constructed3728-byte constants

Retain the prior520-byte header, including ordinary p at448 and D at512. Initialize padding520..527 to zero. Construct FIFTY64-byte slots at528+64j,0<=j<=49. Slot j contains seven ordinary64-bit limbs of j*p plus an initialized zero eighth word. Build slot0 as all zero; each subsequent slot is the exact seven-limb integer sum of the preceding slot and p, with no modular reduction. Because49p<2^448, no significant carry is lost. All fifty slots, header and padding are actually filled and delivered to both workers in this attack; their scalar addition/store/transport cost is included cold. No searched advice, given powers or field-multiplication oracle is supplied.

The same3728-byte block is passed through the existing arguments to the map and key/fixed-s helper. The recurrence reads0..49; fixed-s reads only0..12. All reads are in range and ordinary j*p is not mistaken for a canonical Montgomery field value. Increasing the block from1296 to3728 changes no hot argument count, pointer width, record size or storage protocol.

## Actual fused recurrence QREDUCE30

The actual fused recurrence_fp2_step.S now uses a30-instruction quotient-reselection reducer. The standalone recurrence_fp2_map.S remains inert reference evidence with its older31-instruction addback reducer;it is never called by the actual hotleaf andisnot separatelycharged. The fixed-s helper keeps its own32-instruction reducer unchanged. All fifty ordinaryq*p slots andthe3728-byte block remain fully constructed/paid,including the now-unused hot p header.

Computeqhat=floor(floor(Y/2^325)/D) byEXTR/UDIV,which isq orq+1 bythe full-size proof above. Loadthe seven-limb ordinaryqhat*p using twoADD/threeLDP/oneLDR. Subtract it into the PRODUCT registersx11..17,while preservingthe ORIGINAL Y inx4..10. The final448-bit subtraction borrow is true exactlywhenqhat=q+1. There is no loss ofcarry orsign information;alltrue valuesarefarbelow2^448.

SUBx19,x1,1 computes a tentative qhat-1 WITHOUT settingflags. CSELx1,x19,x1,LO selects thattentativevalue ONLY on the originalborrow. Ifqhat0,borrowisfalse andthetemporarywordwrap isdiscarded;no underflowingquotient isusedasapointer. The selectedquotientisEXACTq=floor(Y/p),always0..48. Recompute its paidtableaddress,reloadordinaryq*p,and subtractinto thepreserved originalaccumulator. The resultis t=Y-q*p in[0,p),withzeroseventhlimb andcanonicalstoredsixlimbs. No p addback ormaskedp registersareneeded. Nothing assumesexactdivision beyond the provedone-high quotientbound.

Exactcount:EXTR/UDIV2,firstaddress2,firstsevenlimbload4,trialSUBS/SBCS7,nonflagSUB/CSEL2,secondaddress2,secondload4,finalSUBS/SBCS7,total30. Fourreductions save4 actualinstructions relative tothe31-instruction reducer. x19 isdead afterthecolumnlookup beforeeachreduction;it isalready saved/restoredbytheleaf. x2 remainsD,x28const,x26address,andthe originalrowB accumulatorisalready savedbeforethefirstreduction. No newcallee,argument,stackslot,branch,oracle ormemory exemption isused. Everyaddress/index/selectedread andinstructionisindependentlychecked bythe sharednativefixture.

The coldconstants/layout/tableentries/quotientproof andallprotectedparametersareunchanged. Fullcurrentglobalworkisinproof.md/architecture.md/radix80_certificate.md. Finitefull-p boundarycases/alias/initializedread/ABI/countchecks are scopedcorroboration,notaproof ofthewholeproductionattack.

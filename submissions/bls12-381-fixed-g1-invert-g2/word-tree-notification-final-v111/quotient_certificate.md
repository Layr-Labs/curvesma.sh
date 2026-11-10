> Current unequal nb/ng counts,ordinalmap,two verification caps and complete probability/work are in proof.md and scalar_stream_certificate.md. Native/table/metadata/sort constructions below retain their scoped contracts. FullS192 is temporary;stored Ctag24 requires verified continuation.

# Exact small-quotient reduction of both full-width accumulator types

Write B=2^325,D=floor(p/B)=58556206801485620=0xd0088f51cbff34>256, p=DB+e with0<=e<B. For every0<=Y<256p let Y=qp+t with0<=t<p,q<=255. Then U=floor(Y/B)=qD+floor((qe+t)/B). The extra is nonnegative and below D+256<2D, so qhat=floor(U/D) is q or q+1. Since256p<2^389, U fits64 bits and the high-limb EXTR#5 computes it exactly. UDIV is one ordinary64-bit scalar instruction priced1 by the protected metric; this is not a CPU-cycle or arbitrary-precision-division assertion.

The current real recurrence accumulator is seed+22p+22positivecanonicalterms-22negativecanonicalterms;all prefixes stay nonnegative and below45p. Imaginary is seed+44canonicalterms,also below45p. Thus q<=44 and qhat<=45;retaining all50ordinaryqp slots is conservative. The full signed range proof is in fp2_certificate.md. For the fixed-s lookup there are eleven canonical terms, givingY<11p and qhat<=11. These are deterministic bounds for every possible input, not sampled constraints.

## Fully constructed3728-byte constants

Retain the prior520-byte header, including ordinary p at448 and D at512. Initialize padding520..527 to zero. Construct FIFTY64-byte slots at528+64j,0<=j<=49. Slot j contains seven ordinary64-bit limbs of j*p plus an initialized zero eighth word. Build slot0 as all zero; each subsequent slot is the exact seven-limb integer sum of the preceding slot and p, with no modular reduction. Because49p<2^448, no significant carry is lost. All fifty slots, header and padding are actually filled and delivered to both workers in this attack; their scalar addition/store/transport cost is included cold. No searched advice, given powers or field-multiplication oracle is supplied.

The same3728-byte block is passed through the existing arguments to the map and key/fixed-s helper. The recurrence reads0..49; fixed-s reads only0..11. All reads are in range and ordinary j*p is not mistaken for a canonical Montgomery field value. Increasing the block from1296 to3728 changes no hot argument count, pointer width, record size or storage protocol.

## Actual recurrence QREDUCE31

The current parameterized source uses EXACTLY31instructions:EXTR/UDIV2,slotaddressADD2,sevenqplimbloads4,SUBS/SBCS7,delayedfullpLDP3,original-borrowCSEL6,carryaddback7. Productregisters19..25aredeadbeforep19..24reloads. ThoseLDP preserve the borrow flags. D2isreadonly;address26/header28 are ordinary pointers. Firstaccumulator4..10reduction NEVERwrites11..17;seconduses11..17directly. No spill/reload/copybridge exists.

If qhat=q, the difference is t in[0,p). If qhat=q+1, it is t-p in[-p,0); the original seven-limb subtraction borrow selects p, giving exactly t. Loads and CSELs preserve that borrow until all six p limbs are selected. All wrapped arithmetic is at448 bits with true values far within range. The final top limb is zero and the six stored limbs are canonical. The seventh p limb is zero, so seven addback instructions suffice.

Accumulator groups are4..10 and11..17;qp scratch19..25,p19..24,quotient1,D2,address26,header28. Output0,input29,tablebase30remainlive;x18neverused. Bothgroupsstayliveatonce. Both reductions reloadpAFTERsubtractionbecauseCSELmaymaskit;loadflagsareunchanged. Fixed-slookupretainsitsOWN32-instructionreducer/Dloadandfull3728byteblockcontract;onlythisprivatevector-mapregisterallocationchanges.

Current two-table/live-accumulator refinement is proved in fp2_certificate.md. Eachrowpair1172;map2*1172+13=2357. Wrapperown290+map2357=2647,caller20=2667retained2672. Privateframe48/publicwrapper528/globalkey816 arecurrent. QREDUCE(a0..a6)keepsits31instructions,butloadsqpinto19..25anddelaysfullp19..24loadsUNTILAFTERtheborrow-producingSUBS/SBCS;LDPpreservesflags. Bothreal/imagremainliveandnooutputarithmetic changes. Everyoriginalquotient/canonicality/one-correctionargumentbelowremainsvalidby0..45p<49p;all50ordinaryq*p slotsareunchanged.

Complete generation=2672+824+376+22+704=4598. Sorting720+scan128+cheap_confirmation1+expensive_verification1+spare256 gives base5704. Complete metadata40 and two80-bit sorting give network933. Hot6637 is retained6712 with75 spare. N=6712L+2^124 satisfies EXACT N^500<20000^500*2^62243, proposed62243/500=124.486. Full parameters, success threshold0.39, C20000, six tracks and one metric remain unchanged. The same self-contained NH2^-63 proof and exact group counting give success0.3900120659047603.

All native code is inert host evidence. Shared-only audit interprets actual new2357/2647 sources against independent dense full-p Fp4 algebra with initialized reads/writes, ABI, alias, counts and native memory-immediate constraints; all fifty product slots are explicit scoped fixtures. It does not certify production fill, probability, physical hardware, acceptance, score or qualification. No new pass is asserted before a genuine shared report. All execution/review/scoring belong to shared llm-verifier.
Currentmapisaprivatehelper:frame48,SP/LRsaved,no19..29preservationpromisedthere. TheONLYpublicrecurrenceentrywrapperalreadyrestoresALL19..29/LR/SP/x18;itsfullABIisunchanged. Privatehelperremoves10redundantsaves/restores;map2357+wrapper290=2647,caller20=2667<=2672. Fullframe/offset/writegraphproofisinfp2_certificate.md.

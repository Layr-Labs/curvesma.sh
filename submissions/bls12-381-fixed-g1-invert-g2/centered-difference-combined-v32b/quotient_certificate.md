# Exact small-quotient reduction of both full-width accumulator types

Write B=2^325,D=floor(p/B)=58556206801485620=0xd0088f51cbff34>256, p=DB+e with0<=e<B. For every0<=Y<256p let Y=qp+t with0<=t<p,q<=255. Then U=floor(Y/B)=qD+floor((qe+t)/B). The extra is nonnegative and below D+256<2D, so qhat=floor(U/D) is q or q+1. Since256p<2^389, U fits64 bits and the high-limb EXTR#5 computes it exactly. UDIV is one ordinary64-bit scalar instruction priced1 by the protected metric; this is not a CPU-cycle or arbitrary-precision-division assertion.

For the actual recurrence map, each output is the sum of46 canonical digit-table terms and one canonical third-history seed. Thus0<=Y<49p, q<=48 and qhat<=49. For the fixed-s lookup there are twelve canonical terms,givingY<12p andqhat<=12. These are deterministic bounds for every possible input, not sampled constraints.

## Fully constructed3728-byte constants

Retain the prior520-byte header, including ordinary p at448 and D at512. Initialize padding520..527 to zero. Construct FIFTY64-byte slots at528+64j,0<=j<=49. Slot j contains seven ordinary64-bit limbs of j*p plus an initialized zero eighth word. Build slot0 as all zero; each subsequent slot is the exact seven-limb integer sum of the preceding slot and p, with no modular reduction. Because49p<2^448, no significant carry is lost. All fifty slots, header and padding are actually filled and delivered to both workers in this attack; their scalar addition/store/transport cost is included cold. No searched advice, given powers or field-multiplication oracle is supplied.

The same3728-byte block is passed through the existing arguments to the map and key/fixed-s helper. The recurrence reads0..49; fixed-s reads only0..12. All reads are in range and ordinary j*p is not mistaken for a canonical Montgomery field value. Increasing the block from1296 to3728 changes no hot argument count, pointer width, record size or storage protocol.

## Actual recurrence QREDUCE31

The source recurrence_fp2_map.S uses these31 instructions: three LDP of p (3), EXTR/UDIV (2), ADD base528 and ADD shifted qhat*64 (2), three LDP plus LDR to obtain all seven product limbs (4), seven SUBS/SBCS (7), six original-borrow CSEL selecting p/zero (6), and seven carry addback instructions (7). The table product is EXACTLY qhat*p, replacing the old23-instruction scalar product sequence with six address/load instructions. D remains in x2 for BOTH reductions of a row pair; the table address uses x26, so no hidden reload is required.

If qhat=q, the difference is t in[0,p). If qhat=q+1, it is t-p in[-p,0); the original seven-limb subtraction borrow selects p, giving exactly t. Loads and CSELs preserve that borrow until all six p limbs are selected. All wrapped arithmetic is at448 bits with true values far within range. The final top limb is zero and the six stored limbs are canonical. The seventh p limb is zero, so seven addback instructions suffice.

Accumulators are x4..10, p x20..25, quotient x1, D x2, table address x26, product x11..17, constant pointer x28. No18 is used; output0,input29,region30 remain live. Row B x11..17 is fully saved to the frame BEFORE the first reduction overwrites those registers and is restored to the accumulator for its reduction. Both reductions reload p because CSEL may mask it. The fixed-s helper keeps its32-instruction reducer, including its own D load, unchanged except the larger constant-block contract.

Current map uses the exact same31-instruction reducer twice per rowpair. All ordinary0..49*p slots remain constructed even though47-term row sums leave slack. Independent shared full-p recurrence/fixed oracles check canonicality,ABI,initialized reads,alias and native instruction counts. Full current scalar/global costs are in proof.md,architecture.md andradix80_certificate.md. No earlier benchmark price or layout is inherited as a current claim.

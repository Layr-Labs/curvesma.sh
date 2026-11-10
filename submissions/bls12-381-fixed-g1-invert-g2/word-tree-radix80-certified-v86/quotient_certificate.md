# Exact small-quotient reduction of both full-width accumulator types

Write B=2^325,D=floor(p/B)=58556206801485620=0xd0088f51cbff34>256, p=DB+e with0<=e<B. For every0<=Y<256p let Y=qp+t with0<=t<p,q<=255. Then U=floor(Y/B)=qD+floor((qe+t)/B). The extra is nonnegative and below D+256<2D, so qhat=floor(U/D) is q or q+1. Since256p<2^389, U fits64 bits and the high-limb EXTR#5 computes it exactly. UDIV is one ordinary64-bit scalar instruction priced1 by the protected metric; this is not a CPU-cycle or arbitrary-precision-division assertion.

For the actual recurrence map, each output is the sum of48 canonical digit-table terms and one canonical third-history seed. Thus0<=Y<49p, q<=48 and qhat<=49. For the fixed-s lookup there are eleven canonical terms, givingY<11p and qhat<=11. These are deterministic bounds for every possible input, not sampled constraints.

## Fully constructed3728-byte constants

Retain the prior520-byte header, including ordinary p at448 and D at512. Initialize padding520..527 to zero. Construct FIFTY64-byte slots at528+64j,0<=j<=49. Slot j contains seven ordinary64-bit limbs of j*p plus an initialized zero eighth word. Build slot0 as all zero; each subsequent slot is the exact seven-limb integer sum of the preceding slot and p, with no modular reduction. Because49p<2^448, no significant carry is lost. All fifty slots, header and padding are actually filled and delivered to both workers in this attack; their scalar addition/store/transport cost is included cold. No searched advice, given powers or field-multiplication oracle is supplied.

The same3728-byte block is passed through the existing arguments to the map and key/fixed-s helper. The recurrence reads0..49; fixed-s reads only0..11. All reads are in range and ordinary j*p is not mistaken for a canonical Montgomery field value. Increasing the block from1296 to3728 changes no hot argument count, pointer width, record size or storage protocol.

## Actual recurrence QREDUCE31

The source recurrence_fp2_map.S uses these31 instructions: three LDP of p (3), EXTR/UDIV (2), ADD base528 and ADD shifted qhat*64 (2), three LDP plus LDR to obtain all seven product limbs (4), seven SUBS/SBCS (7), six original-borrow CSEL selecting p/zero (6), and seven carry addback instructions (7). The table product is EXACTLY qhat*p, replacing the old23-instruction scalar product sequence with six address/load instructions. D remains in x2 for BOTH reductions of a row pair; the table address uses x26, so no hidden reload is required.

If qhat=q, the difference is t in[0,p). If qhat=q+1, it is t-p in[-p,0); the original seven-limb subtraction borrow selects p, giving exactly t. Loads and CSELs preserve that borrow until all six p limbs are selected. All wrapped arithmetic is at448 bits with true values far within range. The final top limb is zero and the six stored limbs are canonical. The seventh p limb is zero, so seven addback instructions suffice.

Accumulators are x4..10, p x20..25, quotient x1, D x2, table address x26, product x11..17, constant pointer x28. No18 is used; output0,input29,region30 remain live. Row B x11..17 is fully saved to the frame BEFORE the first reduction overwrites those registers and is restored to the accumulator for its reduction. Both reductions reload p because CSEL may mask it. The fixed-s helper keeps its32-instruction reducer, including its own D load, unchanged except the larger constant-block contract.

Each row pair now costs1219 instead of1253 (two reductions save34); map2461 replaces2529. Source structure, eleven/twelve-digit tables, field linear maps, exact third/output alias and all other ABI contracts are unchanged. The fused wrapper retains250 own, giving250+2461=2711 throughRET; caller16+rotation4 gives2731 retained2752. The map frame remains208 and wrapper frame528; the existing816-byte worst nested key bound remains conservative.

Complete generation=2752+824+384+32+704=4696. Sorting720+scan128+confirmation1+spare256 gives base5801. Repaired eight-frame metadata gives network1129. Hot6930 is retained7040 with110 spare. N=7040L+2^124 satisfies EXACT N^500<20000^500*2^62287, proposed62287/500=124.574. Full parameters, success threshold0.39, C20000, six tracks and one metric remain unchanged. The same self-contained NH2^-63 proof and exact group counting give success0.3900268131920969.

All native code is inert host evidence. Shared-only audit interprets actual new2461/2711 sources against independent dense full-p Fp4 algebra with initialized reads/writes, ABI, alias, counts and native memory-immediate constraints; all fifty product slots are explicit scoped fixtures. It does not certify production fill, probability, physical hardware, acceptance, score or qualification. No new pass is asserted before a genuine shared report. All execution/review/scoring belong to shared llm-verifier.

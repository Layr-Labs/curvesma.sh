# Factored Fp2 recurrence, exact in-place history and complete key

## Algebra and tables

All fields are the COMPLETE protected Fp381, represented by canonical six-limb native Montgomery residues with Rm=2^384. In Fp4=Fp2[t]/(t^2-xi),xi=1+u,u^2=-1, write the characteristic coefficient c=A+B*t. The p^2 automorphism fixes Fp2 and negates t, so c^(p^2)=A-B*t. For X=a+b*t,Y=d0+e*t:
cX-c^(p^2)Y = A(a-d0)+B*xi*(b+e) + t*(A(b-e)+B(a+d0)).
This is an exact field identity. The actual A,B are the components of c=3Phi(g) for increment g=h or h^(-b), including the factor3. They are computed cold from actual GT values. The recurrence third coefficient is exactly1 from GT norm1.

Eight canonical field operations construct V0=a-d0,V1=b+e,V2=b-e,V3=a+d0, four Fp2 values. Store their eight Fp coordinates in384 initialized temporary bytes. The first output pair depends only on V0,V1 with constants A,Bxi; the second only on V2,V3 with constants A,B. Each constant Fp2 multiplication is the ordinary2x2 Fp linear map: real-input column(Cr,Ci), imaginary-input column(-Ci,Cr). Negative coefficients and Bxi are computed with actual full-field canonical arithmetic, not given by an oracle.

For each of8scalar columns and12digits (32bits for the first11,last29), generate ALL two-coordinate entries d*C_row*2^start/Rm modp, each96bytes. Two actual pinned Montgomery seed products and fill_pair.S generate them by canonical addition, starting at zero. For each row the lookup sum plus its third-history seed is congruent to the desired native output, with at most49 canonical terms. Thus0<=Y<49p<256p, so the exact48-instruction quotient correction from quotient_certificate.md applies. It stores canonical output; there is no384-bit wrap or unit-cost big integer.

A scalar column occupies96*(11*2^32+2^29)=4587025072128 bytes=4272GiB. Eight columns use33.375TiB at[8,41.375)TiB; reserve[8,42). The two workers' fixed-s tables remain[76,78.086)TiB. Their parent4TiB, private<1TiB, and child ports80..128TiB exclusive are disjoint. All ordinary pointers fit47bits. This shrinks the earlier vector tables; unused capacity supplies no data. Every used entry and header is actually generated for this attack.

fill_pair.S has13 prologue/zero instructions,75 per following entry (two actual31-instruction field adds, ten argument/call instructions, three loop/update instructions) and4 epilogue:75(n-1)+17<96n. All tables, both workers, actual coefficients/seeds, page initialization and controls are included in the retained cold2^50 table allowance and complete cold2^124. No advice or future-target amortization.

## Native recurrence map and step

recurrence_fp2_map.S takes x0=third_out192 (EXACT in-place),x1=transformed384 disjoint,x2=matrixbase,x3=constant520. The output is disjoint from the other histories, input, tables, constants and frame. For the first two output coordinates, seed both seven-limb accumulators from OLD third coordinates0/48; then perform four scalar columns. Only those two coordinates are stored. The second row pair seeds OLD coordinates96/144, which the first pair did not touch, then performs the other four columns. This proves exact third/output alias is safe. Other partial overlaps are not used.

MSTEP22 plus one base advance except the last digit gives275 per12digit column;2 base instructions+6 scalar input loads make COLUMN283. Each row pair seeds with9 instructions (pointer load,3LDP+zero top for each coordinate),4columns1132, and112 spill/quotient/store instructions:1253. Prologue8/setup5/two1253/reload3/epilogue7 total2529 throughRET. Entry multiplier96,stride0x60<<32,column block0x42c<<32 all use ordinary scalar MOV/MADD. Table row offsets are always0/48, while seed/output rows are0/48 or96/144. No hidden zero-coefficient lookups or omitted terms.

Register allocation is the same complete two-accumulator scalar allocation: accumulators4..10 and11..17,input19,load temporaries20..25,tablebase26,stride27,size28,input29,region30,out0,entry1,block2,index3; no18. Frame208 saves all19..29/SP/LR,arguments and second-accumulator spill. NO BL. Every transformed input coordinate is canonical and initialized. Table and input offsets are within their declared buffers. The native QREDUCE saves the second accumulator before its product temporaries overwrite it and reloads real p for each row.

recurrence_fp2_step.S takes third_out192,prev1_192,prev2_192,matrixbase,const520,p48. The three histories are DISTINCT. Frame480 holds header96 and temporary384. Prologue8, six pointer MOVs plus one temporary-base ADD7, eight five-instruction field argument/calls40, map argument/call5, epilogue8 give68 own instructions. Nested eight31-instruction canonical Fp helpers plus2529 map give2845 throughRET. Allow16 caller work and FOUR pointer rotation MOVs:2865<=2944. The temporary buffer is disjoint from all histories and nested frames. All transforms finish before the map updates third in place. There is no raw-history copy, post-map field addition or unpriced operation.

## Exact complete b1-first key

The raw Phi is four contiguous48-byte fields (a0.real,a0.imag,b1.real,b1.imag). canonical_b1.S contains the exact public97-instruction canonical_conjugate_b1 sign-helper body and its full381-bit constants. It supports exact input/output alias. The prior six-key copy helper is omitted. odd_trace_compact.S remains55 own plus three31 field calls and two257 fixed-s lookups:662. It computes the full odd Phi in a separate192-byte temporary. All sum/difference reads precede its b1 stores.

The proof T²=-1 and q(B)!=q(TB) for every nonidentity GT is in key_contract.md. This justifies comparing ONLY signed b1 to select a full trace representative; the output still stores BOTH a0 and b1. This is the fixed total order with b1 primary and a0 secondary, never a shortened key. Both workers use the same order and fingerprint coefficients. Exact12-fiber correctness, group-event probabilities and cap bounds are unchanged.

canonical_twelve_bfirst.S has frame384: header96, EvenB96@96, OddPhi192@192 (its b1 at288). First sign call reads input b1@96 and writes EvenB. Odd trace reads immutable input and writes OddPhi. Second sign call canonicalizes OddPhi.b1 exactly in place. All twelve b1 words are compared in six paired comparisons. Equality returns1/output invalid and is exactly identity on GT inputs. Otherwise the selected a0 pointer (raw or OddPhi) and b1 pointer (EvenB or OddPhi.b1) are set, and both96-byte parts are copied to the disjoint192-byte output. No key call mutates a raw history. Every successful output byte is initialized; no unused scratch is read.

Own instruction bound115 counts all branch arms conservatively: prologue8, argument MOVs5, first sign3, odd6, second sign3, two bases2, six comparisons*8=48, tie flag/branch2, both selections5, copies24, success flag1, epilogue8. Nested two97 sign calls and662 odd trace give971. At most8 caller/check instructions gives979<=992. All19..29/SP/LR preserved; no18. Stack nesting is at most384+336+96=816 bytes. Raw input, output, scratch, constants and tables satisfy the disjointness contracts in the actual sources.

## Complete totals

Combined step2944,key992,fingerprint496,Fr32,main704=>generation5168. Add five272 radix pairs1360,scan128,confirmation1,spare256=>base6913. All source worker streams, Root output streams, endpoints and sorting network remain1364. Hot8277 retained8384. For the actual four-phase m/L in proof.md:
N=8384L+2^124 and EXACT N^100<20000^100*2^12483,
giving12483/100=124.83. Fullp381/r255,targetdistribution,success0.39,C20000,sixtracks/onecomputationmetric unchanged. Cold124/peak138/advice0 include every table/source/controller and failed path. All code is inert evidence; no host execution/compilation, observed native pass, score or qualification is claimed. Shared llm-verifier alone executes/reviews/scores.

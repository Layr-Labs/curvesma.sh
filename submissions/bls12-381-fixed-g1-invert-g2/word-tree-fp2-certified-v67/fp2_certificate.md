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

## Exact complete compact key

The raw Phi history is four contiguous48-byte fields: a0.real,a0.imag,b1.real,b1.imag. conjugate_compact.S retains the EXACT97-instruction canonical_conjugate_b1 sign helper from the attributed source and maps canonical_six_trace's raw loads from384..464 to96..176. Its52-instruction six-key copy/identity/sign-call body is otherwise unchanged. Inclusive149 is a safe upper bound, counting all branch paths conservatively. odd_trace_compact.S has the same55 own instructions, three31 field calls and two257 fixed-s lookups as before:662. Only its raw b1 offsets change to96/144. Sum/difference reads precede b1 writes, so its documented exact alias case remains safe. These changes alter layout, not field arithmetic.

canonical_twelve_compact.S is a COMPLETE1123-instruction key bound:163 own plus2*149 six-key calls plus662 odd Frobenius. Frame672: header96,keyA192@96,keyB192@288,odd temporary192@480. First six-key identity flag is checked; odd Frobenius and second six-key are computed in separate buffers; the second flag is checked. Twelve paired word comparisons examine all24 words in increasing index with unsigned comparisons; choose the first differing word, copy the selected192 bytes, and return0. Return1 marks identity/output invalid; the caller aborts and does not read output. At most8 caller/check instructions gives1131<=1152.

The163 own count includes all branch arms as a conservative upper bound:prologue8,five argument MOVs5,first-call/check4,odd arguments/call6,second-call/check4,two base pointers2,12*8 comparisons96,selection3,copy24,flags/branches3,epilogue8. All19..29/SP/LR preserved,no18. Output S192 is DISJOINT from all histories; no key routine can mutate the raw recurrence state. Worst stack nesting672+336+96=1104 bytes. All192 output bytes are initialized on success and no uninitialized workspace hole is read.

The key is min(F(X),F(X^p)) under this FIXED total word order. Exact12-fiber correctness requires only a fixed total order, not a particular byte order. Both workers use this exact helper and the same fingerprint coefficients. Therefore the complete group-event/cap probability from proof.md is unchanged, including the color convention and all baby duplicates.

## Complete totals

Combined step2944,key1152,fingerprint496,Fr32,main704=>generation5328. Add five272 radix pairs1360,scan128,confirmation1,spare256=>base7073. All source worker streams, Root output streams, endpoints and sorting network remain1364. Hot8437 retained8544. For the actual four-phase m/L in proof.md:
N=8544L+2^124 and EXACT N^50<20000^50*2^6243,
giving6243/50=124.86. Fullp381/r255,targetdistribution,success0.39,C20000,sixtracks/onecomputationmetric unchanged. Cold124/peak138/advice0 include every table/source/controller and failed path. All code is inert evidence; no host execution/compilation, observed native pass, score or qualification is claimed. Shared llm-verifier alone executes/reviews/scores.

# Factored Fp2 recurrence, exact in-place history and complete key

## Algebra and tables

All fields are the COMPLETE protected Fp381, represented by canonical six-limb native Montgomery residues with Rm=2^384. In Fp4=Fp2[t]/(t^2-xi),xi=1+u,u^2=-1, write the characteristic coefficient c=A+B*t. The p^2 automorphism fixes Fp2 and negates t, so c^(p^2)=A-B*t. For X=a+b*t,Y=d0+e*t:
cX-c^(p^2)Y = A(a-d0)+B*xi*(b+e) + t*(A(b-e)+B(a+d0)).
This is an exact field identity. The actual A,B are the components of c=3Phi(g) for increment g=h or h^(-b), including the factor3. They are computed cold from actual GT values. The recurrence third coefficient is exactly1 from GT norm1.

Eight canonical field operations construct V0=a-d0,V1=b+e,V2=b-e,V3=a+d0, four Fp2 values. Store their eight Fp coordinates in384 initialized temporary bytes. The first output pair depends only on V0,V1 with constants A,Bxi; the second only on V2,V3 with constants A,B. Each constant Fp2 multiplication is the ordinary2x2 Fp linear map: real-input column(Cr,Ci), imaginary-input column(-Ci,Cr). Negative coefficients and Bxi are computed with actual full-field canonical arithmetic, not given by an oracle.

For each of8scalar columns and12digits (32bits for the first11,last29), generate ALL two-coordinate entries d*C_row*2^start/Rm modp, each96bytes. Two actual pinned Montgomery seed products and fill_pair.S generate them by canonical addition, starting at zero. For each row the lookup sum plus its third-history seed is congruent to the desired native output, with at most49 canonical terms. Thus0<=Y<49p<256p, so the exact31-instruction quotient correction from quotient_certificate.md applies. It stores canonical output; there is no384-bit wrap or unit-cost big integer.

A scalar column occupies96*(11*2^32+2^29)=4587025072128 bytes=4272GiB. Eight columns use33.375TiB at[8,41.375)TiB; reserve[8,42). The two workers' fixed-s tables occupy[56,69.5)TiB. Their parent4TiB, private<1TiB, and child ports80..128TiB exclusive are disjoint. All ordinary pointers fit47bits. This shrinks the earlier vector tables; unused capacity supplies no data. Every used entry and header is actually generated for this attack.

fill_pair.S has13 prologue/zero instructions,75 per following entry (two actual31-instruction field adds, ten argument/call instructions, three loop/update instructions) and4 epilogue:75(n-1)+17<96n. All tables, both workers, actual coefficients/seeds, page initialization and controls are included in the retained cold2^50 table allowance and complete cold2^124. No advice or future-target amortization.

## Native recurrence map and step

recurrence_fp2_map.S takes x0=third_out192 (EXACT in-place),x1=transformed384 disjoint,x2=matrixbase,x3=constant3728. The output is disjoint from the other histories, input, tables, constants and frame. For the first two output coordinates, seed both seven-limb accumulators from OLD third coordinates0/48; then perform four scalar columns. Only those two coordinates are stored. The second row pair seeds OLD coordinates96/144, which the first pair did not touch, then performs the other four columns. This proves exact third/output alias is safe. Other partial overlaps are not used.

MSTEP22 plus one base advance except the last digit gives275 per12digit column;2 base instructions+6 scalar input loads make COLUMN283. Each row pair seeds with9 instructions (pointer load,3LDP+zero top for each coordinate),4columns1132, and78 spill/quotient/store instructions:1219. Prologue8/setup5/two1219/reload3/epilogue7 total2461 throughRET. Entry multiplier96,stride0x60<<32,column block0x42c<<32 all use ordinary scalar MOV/MADD. Table row offsets are always0/48, while seed/output rows are0/48 or96/144. No hidden zero-coefficient lookups or omitted terms.

The unchanged two-accumulator layout saves Row B before reducer scratch overwrites it. New QREDUCE uses31 instructions and the explicitly constructed fifty-entry ordinary q*p table described in quotient_certificate.md. Its header still supplies p@448,D@512. Each row pair costs1219; the full map2461 is statically expanded from its actual macros. All table terms, alias and canonical-output proofs above remain unchanged.

## Fused canonical sum/difference transform, full native contract

The SAME four Fp2 values V0=a-d0,V1=b+e,V2=b-e,V3=a+d0 are now computed by four inline BOTH macros, one for each pair of Fp components, rather than eight separate helpers. This changes no recurrence algebra, coefficients, table contents, output key, probability event or alias contract. All field values are full native canonical residues modulo the protected381-bit p, with Montgomery Rm2^384.

Each BOTH reads six limbs A into x4..9, six B into x10..15 and p into x19..24. It forms SUM=A+B in x25..30 (6 carry operations). Since2p<2^384, no significant carry is lost. It then forms A-B in x4..9 (6 subtracts), preserving its ORIGINAL final borrow across six CSELs that replace x10..15 by p or zero. Six carry additions give the canonical difference: if borrowed, the384-bit wrapped A-B plus p wraps once to A-B+p in[1,p); otherwise it remains A-B in[0,p). Store the six-limb difference (3STP). Subtract p from SUM into the now-dead x4..9 (6), select the trial or original SUM using the final HS flag (6), and store the canonical sum (3STP). Every CSEL preserves the required original flags. All input loads precede stores.

BOTH costs9 input/p loads+6 sum+6 difference+6 borrow masks+6 addback+3 difference stores+6 sum trial+6 selects+3 sum stores=51, with no branch, BL or x18 use. It intentionally uses x19..30 including frame-pointer/LR registers. The WRAPPER saves and restores all callee registers and LR; it does not treat these macro values as frame pointers. All frame addressing uses SP exclusively, never the clobbered x29. No ABI preservation is falsely attributed to the macro in isolation.

Wrapper frame528 has saved29/LR@0,19..28@16..95, six argument words@96..143, initialized temporary384@144..527. Three STP save x0/x1,x2/x3,x4/x5 before any macro. Per component, ADD establishes plus/minus output pointers; LDP at SP104 loads prev1/prev2; two ADD offsets and LDR p atSP136 complete SIX setup instructions. For raw offsets0/48, PLUS is Temp288/336 and MINUS Temp0/48. For raw96/144, PLUS is Temp96/144 and MINUS Temp192/240. Thus Temp remains exactly the table map's expected V0,V1,V2,V3 order, and all384 bytes are written before the map reads them.

After four macros, LDR x0@SP96 restores third/output, ADD x1=SP144 gives Temp, LDP x2/x3@SP120 restores table/constant pointers, and BL recurrence_fp2_map performs the unchanged2461-instruction in-place update. Its output aliases ONLY the old third, with previous1/previous2/Temp disjoint, exactly as the prior map proof. All histories remain canonical raw Phi. The key output stays separate.

Wrapper prologue7 (SUB,STP29/LR,five STP callee pairs), argument saves3, four*(setup6+BOTH51)=228, map setup/call4, epilogue8 total250 own. Nested map2461 gives2711 throughRET. Caller16 plus four history-pointer MOVs give2731; retain2752. All19..29/SP/LR preserved on return, no18. Maximum recurrence stack528+208=736; the key's prior816 remains a larger overall bound. Unaligned-by16 LDP offsets104/120 are multiples of8 and valid scalar AArch64 pair loads; no atomic/alignment guarantee is presumed beyond ordinary mapped byte-addressed RAM.

The shared-only audit now interprets this ACTUAL fused wrapper/macro and unchanged map against independent dense full-p Fp4 algebra. Its strict stack bound is736 and expected wrapper/map counts2711/2461. This is candidate experiment evidence, not an asserted pass before a genuine shared report. No candidate code was executed or compiled on the host.

## Exact complete b1-first key

The raw Phi is four contiguous48-byte fields (a0.real,a0.imag,b1.real,b1.imag). canonical_b1.S implements the public sign convention with a caller-saved leaf of worst path65 and the unchanged full381-bit constants. It supports exact input/output alias. The prior six-key copy helper is omitted. odd_trace_compact.S remains55 own plus three31 field calls and two212 fixed-s lookups:572. It computes the full odd Phi in a separate192-byte temporary. All sum/difference reads precede its b1 stores.

The proof T²=-1 and q(B)!=q(TB) for every nonidentity GT is in key_contract.md. This justifies comparing ONLY signed b1 to select a full trace representative; the output still stores BOTH a0 and b1. This is the fixed total order with b1 primary and a0 secondary, never a shortened key. Both workers use the same order and fingerprint coefficients. Exact12-fiber correctness, group-event probabilities and cap bounds are unchanged.

canonical_twelve_bfirst.S has frame384: header96, EvenB96@96, OddPhi192@192 (its b1 at288). First sign call reads input b1@96 and writes EvenB. Odd trace reads immutable input and writes OddPhi. Second sign call canonicalizes OddPhi.b1 exactly in place. All twelve b1 words are compared in six paired comparisons. Equality returns1/output invalid and is exactly identity on GT inputs. Otherwise the selected a0 pointer (raw or OddPhi) and b1 pointer (EvenB or OddPhi.b1) are set, and both96-byte parts are copied to the disjoint192-byte output. No key call mutates a raw history. Every successful output byte is initialized; no unused scratch is read.

The tighter complete path accounting is in key_path_certificate.md: the same source has own CFG maximum110, conservatively retain111; the sign leaf is65 by an exhaustive zero/nonzero partition; odd trace is572. Thus the complete key is813, caller8 gives821 retained824. All raw-history separation, initialized output, identity status, no18 and full callee/LR/SP contracts remain unchanged, with worst key stack816.

## Complete totals

Combined step2752,key824,fingerprint384,Fr32,main704=>generation4696. Add five272 radix pairs1360,scan128,confirmation1,spare256=>base6441. All source worker streams, Root output streams, endpoints and sorting network remain1369. Hot7810 retained7896. For the actual four-phase m/L in proof.md:
N=7896L+2^124 and EXACT N^1000<20000^1000*2^124739,
giving124739/1000=124.739. Fullp381/r255,targetdistribution,success0.39,C20000,sixtracks/onecomputationmetric unchanged. Cold124/peak138/advice0 include every table/source/controller and failed path. All code is inert evidence; no host execution/compilation, observed native pass, score or qualification is claimed. Shared llm-verifier alone executes/reviews/scores.

## New caller-saved leaf sign helper

canonical_b1.S implements the same first-nonzero-component sign convention as the public97-instruction helper, but is NEW code, not an identical retained body. It changes only x2..15 and condition flags. x0 input and x1 output remain stable; x18,x19..30,SP are unchanged and no BL occurs. Both pointers are96-byte native canonical Fp2 buffers, exact alias or disjoint. Partial overlaps are forbidden. The constants p and ceil(p/2) are unchanged full381-bit integers.

Load Br and OR all six limbs. If Br=0 load Bi as the sign-determining component. Compare all six limbs with ceil(p/2) using SUBS/SBCS. If lower, copy both original components unchanged. Otherwise, for EACH component separately OR its six limbs; if zero, write canonical zero explicitly, otherwise compute p-component with a full six-limb SUBS/SBCS chain. Because the original nonzero component is in[1,p), its negation is in[1,p), with no underflow or lost carry. Thus zero never becomes the noncanonical encoding p. The first component is fully read before its stores, and writes to it cannot affect the second component, proving exact alias safety. All successful output bytes are written.

The previous68 sign count combined incompatible branches and is superseded by the exhaustive65 bound in key_path_certificate.md. Its highest path requires nonzero real and imaginary components:21 decision prefix+5 p loads+19 each component+RET1=65. If real is zero the prefix is24 but its component path is12 rather than19, giving at most61. Low-copy paths are at most37. This changes no source arithmetic or ABI.



A static acyclic control-flow analysis of canonical_b1.S, counting both successors of every conditional branch without running instructions, gives maximum path length68 throughRET. This is a source bound, not an execution test.

The new exact fixed-s lookup has212 instructions and11 wider digits; its complete table, constant-block extension, quotient proof, counts and memory map are in constant_lookup.md. This changes only the key stage, not any fused recurrence operation or probability event.

The recurrence reducer and every product slot are now specified by quotient_certificate.md. The fixed-s helper reads the same first twelve product slots of the expanded3728-byte block. No partial old product table is supplied as advice.

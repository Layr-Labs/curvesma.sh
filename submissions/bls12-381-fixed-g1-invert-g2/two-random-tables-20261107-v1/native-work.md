# Conservative scalar-native implementation bounds

Analytical instruction-equivalent bounds are explicitly permitted by COST_MODEL.md;assembly or full-size execution is not required. All primitives below are scalar AArch64 integer operations with64-bit limbs,including memory/address/control/interface work. No arbitrary-width unit-cost operation,hidden normalization or SIMD is used. C remains20000. Group reference calls are counted once separately;their internal field instructions are never counted again.

## Explicit unrolled six-limb Montgomery product retained2048

coarse-montgomery.S is self-contained scalar analytic source evidence,never host-compiled/executed. It implements the same conventional Montgomery construction described below. Let B=2^64,R=B^6,a,b in[0,p),n0=0x89f3fffcfffcfffd with n0*p0=-1 modB. The native representation is xR modp. A112-byte aligned stack frame holds13 zeroed scratch limbs. Stage1 adds all36 a[i]*b[j] products at limb k=i+j;stage2 computes six m=low64(t[i]*n0) and adds all36 m*p[j] products at k=i+j. Each product uses MUL/UMULH and carries through EVERY higher scratch limb up to12. The public assembler loops expand into fixed straight-line code;no runtime loop or data-dependent carry shortcut exists.

At a product index k,there are13-k load/ADDS-or-ADCS/store triples. Stage1 adds two operand loads and two product instructions,so the cell count is4+3(13-i-j). Its36 cells sum to1008. Stage2 reuses m in a register,so each cell has one p load,two product instructions and the same carry triples:3+3(13-i-j),summing972. The source has exactly54 other instructions:entry/output pointer/constant bindings6,scratch zero stores13,six m load/multiply pairs12,final u/p loads/subtract/select/stores21,and stack restore/RET2. Thus1008+972+54=2034,retained2048 with14 spare. Every actual address and data load/store is counted. All offsets<=96 are encodable. No wrapper,PLT or veneer is used:place source/constants/callers in the same small shared static image and direct branch/ADRP range,initialize/charge all code/data cold.

Correctness:after stage1 t=a*b. Reduction index i makes limb i zero because m*p0=-t[i] modB;earlier zero limbs stay zero. At the end t=a*b+p*M for0<=M<R,so u=t/R is the standard Montgomery result. Since a*b<p^2 and p<R,u<p^2/R+p<2p;one final subtraction suffices. Every partial sum is nonnegative and<p^2+pR<2^768,within13*64 scratch bits. Final u<2p<2^382 makes t12=0,so reading t6..t11 loses no limb. Product low/high and carry propagation use native64-bit registers,not wide-integer primitives. Scratch is explicitly zeroed inside each call. Outputs have the exact pinned canonical scaling.

Register lifetimes are explicit:input pointers x0/x1 survive both cell stages,output pointer saved in x15,p base x14,n0 x13,m x9,products x5/x6,scratch word x7. The final normalization overwrites dead input/scratch registers only. x18,x19..x30 and original LR are untouched;SP subtract/restore112 preserves16-byte alignment. All operands are read before output stores and the canonical parent uses disjoint buffers. No hidden ABI/callee work is omitted.

Canonical add/sub/neg each fit64 instructions:load six-limb operands,p where needed,propagate a six-step carry/borrow,conditionally correct and write six limbs. E.g. add uses6 operand LDP,6 sum steps,2 p-address+3 p-LDP,6 trial-subtract steps,6 CSEL,3 output STP,RET and at most3 pointer moves=36. Sub uses the same loads,6 subtract steps,one borrow mask,6 AND,6 add-back,3 stores,RET/pointers=37. Negation loads3 LDP,5 ORR/CBZ,binds/loads p in5,subtracts in6,stores3/RET=24 worst;zero stores are cheaper. All include full p and canonical results,not reduced parameters.

## Remaining ordinary native work

Trace canonicalization always computes all four images. Two Montgomery products cost4096;one add,one sub and five negations cost<=7*64=448. Retain512 more for raw projection/all576-byte workspace initialization as needed per call,arguments/BL/frame/address/k/status/RET. Comparing three candidate keys to the current minimum and writing the final192-byte key fits1024:at most3*24 word comparisons,each two loads/CMP/conditional branch,plus complete key copying,coefficient high-to-low traversal,k and preservation. Total<=6080;retain6144. It is a deterministic worst-case bound with no parity/zero cost cap.

The exact supplied nh-fingerprint.S has460 instructions:entry5,first pair26,eleven other pairs396,compression30,restore/RET3. Retain512,including every actual seed/input load,full128 accumulator,MUL/UMULH,odd-key compression,BE output and frame. hash.md proves its fixed-pair/coupling bound and supplies all65 cold seed messages.

The generic-hash.S full standard compressor has56G*16+32=928 instructions,retained1024. It executes all seven rounds and correct CV/IV/zero counter/len/flags/feed-forward;one fixed macro schedule is supplied per round. Header/footer include every pointer/constant/callee frame/argument-state capture/output/RET. x18/LR are untouched;all used callee registers restored. All used messages are single-chunk,counter0. The prefix is paid cold and ordinary PAIR64/bucket40 use their actual lengths/flags. All hash message construction/caller arguments are in main512;extra prefix-CV copying,if used,is within the96 spare. No unpriced finalization,PLT,veneer or reduced-round primitive is used.

Each ordinary sample has one full64-byte PAIR compression and one final40-byte bucket compression,each charged1024. All SCAL sampling/rejections are cold work for the2W tables;there is no ordinary full-Fr rejection loop. Native main512 covers128-bit sample counter/PAIR message encoding,disjoint92-bit extraction,two selectors<=64,total raw576 copy72,canonical scalar copy,all dense-reference/field/hash/dictionary arguments/BL,raw-buffer management,prefix/state handling and branches/preservation. Table pointers are consumed before bank switches. Full four-limb beta=v+w addition modulo actual r fits256 including inputs,r loads,carry/trial subtraction/select,output andRET;no scalar parameter is truncated.

storage.md bounds dictionary work by544L+256. Retain1024 per sample;its256 fixed remainder is in coldD. All head/body selectors,full payload initialization,key/fingerprint comparisons,allocation/publication,hard visit counter and frame/RET are included.

|Ordinary scalar work per accepted sample|Retained instructions|
|---|---:|
|Full trace canonicalization|6144|
|NH fingerprint|512|
|Full standard bucket BLAKE3|1024|
|Full standard PAIR BLAKE3|1024|
|Full Fr beta addition|256|
|Sample control/two table selectors and copy|512|
|Dictionary including capped visits|1024|
|Subtotal|10496|
|Retained complete bound|10752|

The256 spare instructions cover additional interfaces without lowering any quoted bound. L samples use exactly L dense-GT reference operations,each charged1. All initial/terminal work,all code/constants/workspaces and every fixed remainder are per-target cold paid,not advice or amortization. No execution or official score is claimed.

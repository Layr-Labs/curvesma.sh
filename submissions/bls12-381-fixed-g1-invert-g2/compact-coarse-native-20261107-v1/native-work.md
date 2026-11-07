# Conservative scalar-native implementation bounds

Analytical instruction-equivalent bounds are explicitly permitted by COST_MODEL.md;assembly or full-size execution is not required. All primitives below are scalar AArch64 integer operations with64-bit limbs,including memory/address/control/interface work. No arbitrary-width unit-cost operation,hidden normalization or SIMD is used. C remains20000. Group reference calls are counted once separately;their internal field instructions are never counted again.

## Explicit six-limb Montgomery product below4096 instructions

Let B=2^64,R=B^6,a,b in[0,p),n0=0x89f3fffcfffcfffd with n0*p0=-1 modB. The native representation is xR modp,as in the pinned tower. Use13 zero-initialized scratch limbs t0..t12,then these fixed unrolled steps:

1. For i,j in0..5,compute the128-bit product a[i]*b[j] with MUL/UMULH and add it into t starting at limb i+j. Propagate carry through every higher scratch limb up to12,even if zero.
2. For i in0..5,set m=low64(t[i]*n0). For j in0..5,add the128-bit product m*p[j] into t at i+j using the same fixed carry propagation.
3. Read u from t6..t11. Compare/subtract p once and select the canonical result;write all48 output bytes.

Each product-add cell uses at most two operand LDR,two MUL/UMULH,and thirteen load/ADDS-or-ADCS/store triples,at most43 instructions. The first triple starts with ADDS,then the high product limb uses ADCS,and all subsequent limbs add only carry. Every offset is a fixed encodable displacement<=96. Unroll the public six-by-six loops,so no data-dependent or unpriced loop/address computation remains. There are exactly36 product cells and36 reduction cells. The latter reuse m in a register;the six m=LDR/MUL operations are outside those cells. An additional512 instructions covers all scratch initialization,constant/pointer bindings,n0/m operations,callee/LR preservation,argument setup,final compare/subtract/select,all output stores andRET. Thus72*43+512=3608<4096. The unrolled code is below2^16 instructions and all its initialization is paid cold.

Correctness:stage1 gives t=a*b. At reduction index i,the low limb becomes0 because m*p0=-t[i] modB;earlier zero limbs stay zero. At the end t=a*b+p*M for0<=M<R,so u=t/R is the standard Montgomery product a*b/R modp. Since a*b<p^2 and p<R,u<p^2/R+p<2p,so one subtraction canonicalizes it. Partial sums are nonnegative and bounded by p^2+pR<2^768,well within13*64 scratch bits;final u<2p<2^382 means t12=0. No quotient instruction or overflowing carry is discarded. The result has exactly the pinned canonical Montgomery scaling. Scratch/output are disjoint from live operands;all t limbs are explicitly zeroed in the512 allowance.

Canonical add/sub/neg each fit64 instructions:load six-limb operands,p where needed,propagate a six-step carry/borrow,conditionally correct and write six limbs. E.g. add uses6 operand LDP,6 sum steps,2 p-address+3 p-LDP,6 trial-subtract steps,6 CSEL,3 output STP,RET and at most3 pointer moves=36. Sub uses the same loads,6 subtract steps,one borrow mask,6 AND,6 add-back,3 stores,RET/pointers=37. Negation loads3 LDP,5 ORR/CBZ,binds/loads p in5,subtracts in6,stores3/RET=24 worst;zero stores are cheaper. All include full p and canonical results,not reduced parameters.

## Remaining ordinary native work

Trace canonicalization always computes all four images. Two Montgomery products cost8192;one add,one sub and five negations cost<=7*64=448. Retain512 more for raw projection/all576-byte workspace initialization as needed per call,arguments/BL/frame/address/k/status/RET. Comparing three candidate keys to the current minimum and writing the final192-byte key fits1024:at most3*24 word comparisons,each two loads/CMP/conditional branch,plus complete key copying,coefficient high-to-low traversal,k and preservation. Total<=10176;retain12288. It is a deterministic worst-case bound with no parity/zero cost cap.

The five-lane NH128/compression fingerprint fits2048:sixty lane products,each at most16 instructions for two seed words,two mod64 adds,MUL/UMULH,full128 accumulation and addressing/control,plus at most256 for all input loads,seed cursor,odd128 multiplier compression,output endian serialization,frame andRET. Thus<=60*16+256=1216<2048. hash.md defines every operation/byte.

The full standard final BLAKE3 compression fits2048:56 G invocations,each14 add/xor/rotate plus two message loads=896,plus512 for initialization/constants,all CV copies,message padding,output feed-forward,digest serialization,addressing,frame andRET. The public message permutation is compiled into fixed message offsets;no runtime computation is omitted. All seven rounds execute. Prefix work is paid cold. hash.md supplies the exact full compression specification.

Full four-limb canonical Fr-B update fits256 (two operands,eight limb loads,4 add/carry,4 trial subtraction,4 select,four output stores,r/argument/address/control/RET). Main control and jump selection fit512:128-bit K comparison/counter increment,all actual pointer arguments/BL,raw-buffer swaps,label/bucket extraction,branches,context preservation and the explicit97-bit bank latches/commit/barriers. A64 two64-bit limbs implement each128-bit counter/index;the jump index is94 bits,not a native94-bit instruction.

storage.md bounds dictionary work by544K+768. Retain1024 per transition;its768 fixed remainder is in coldD. All head/body selectors,full payload initialization,key/fingerprint comparisons,allocation/publication,hard visit counter and frame/RET are included.

|Ordinary scalar work per transition|Retained instructions|
|---|---:|
|Full trace canonicalization|12288|
|NH fingerprint|2048|
|Full standard final BLAKE3|2048|
|Full Fr-B update|256|
|Main/control/jump selector|512|
|Dictionary including capped visits|1024|
|Subtotal|18176|
|Retained complete bound|20000|

The1824 spare instructions cover additional interfaces without lowering any quoted bound. K transitions use exactly K dense-GT reference operations,each charged1. All initial/terminal work,all code/constants/workspaces and every fixed remainder are per-target cold paid,not advice or amortization. No execution or official score is claimed.

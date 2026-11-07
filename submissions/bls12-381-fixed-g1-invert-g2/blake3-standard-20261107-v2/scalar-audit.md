# Complete scalar work outside one reference multiplication

All counts are scalar AArch64 instructions, including loads, stores, addressing
and control. No SIMD, arbitrary-width unit primitive or free oracle is used.
C remains the protected20000. Reference multiplication is charged as H=1;
its own instructions, stack and ABI preservation are not charged again.

## Norm-one injective encoding and sign:128

For X=a+bw, conjugation X->a-bw is p^6 Frobenius. Public exponentiation gives
p^6=-1 mod r, so a²-vb²=1 for every GT element. Given b, at most two a roots
exist, a and-a; if a=0 there is one. The unsigned canonical Montgomery
coefficient tuple's lexicographic sign distinguishes them. Thus six48-byte b
coefficients plus one sign byte are injective. No square root or decompression
is performed; full states stay stored; exact injective keys are compared.

At the first nonzero a coefficient x, sign is x>floor(p/2), equivalently
x>=H=(p+1)/2; earlier zero coefficients equal their negations. Fully unroll six
coefficient tests, each six LDR, five ORR, one CBNZ: at most72 instructions.
Keep the selected six limbs, load six H limbs, compare by SUBS/five SBCS and
CSET CS:13. All-zero case outputs zero. Entry pointers, optional BL, 64-bit sign
store, return and control fit32. Total<=117; charge128. Only caller-saved
registers and immediate offsets below288/48 are needed. b serialization is
charged in hash preparation, not treated as free. Canonical Montgomery residues
are unique modulo p, so tuple and field equality are consistent.

## Exact simple-tabulation lookup:1024 including serialization/interface

lookup.S and exact checkpoint keys are unchanged. The tabulation compressor consumes the same37-word injective key, using72 unsigned32-bit halves of the36 b words plus one0/1 sign digit. On the specified little-endian AArch64 machine, LDR W directly reads low then high halves at fixed offsets0,4,...,284. LDR X reads the sign at288. These digits together uniquely determine the key; field parameters and sign convention are unchanged. Tables contain five native64-bit words per entry, are independently initialized as described in proof.md, and are aligned at least8 bytes. Entry stride40 preserves that alignment.

Use x0 key base,x1 array of73 table bases,x2 output base; x3..x7 hold five XOR accumulators,x8 digit,x9 table base,x10 entry pointer,x11/x12 data temporaries. Only caller-saved registers are used, without stack or internal calls. Initialize five accumulators to zero. For every position fully unroll the following12 scalar instructions, with literal offsets appropriate to that position:

    LDR w8,[x0,#digit_offset]       // LDR x8,[x0,#288] for final sign
    LDR x9,[x1,#table_base_offset]
    ADD x10,x8,x8,LSL #2
    ADD x10,x9,x10,LSL #3
    LDP x11,x12,[x10]
    EOR x3,x3,x11
    EOR x4,x4,x12
    LDP x11,x12,[x10,#16]
    EOR x5,x5,x11
    EOR x6,x6,x12
    LDR x11,[x10,#32]
    EOR x7,x7,x11

The first ADD forms5*digit and the second tablebase+40*digit. Indices fit64 bits; digit<=2^32-1 (sign<=1). All offsets for key and base-array loads are immediate, aligned and in initialized storage. Every selected entry is read in its entirety and all addressing/transfers are charged; no memory oracle or assumed cache hit is used. Full cost73*12=876. Five initial MOV, five REV plus five STR to serialize40 big-endian output bytes, and RET add16, giving892. Charge1024 including132 additional scalar instructions for caller moves, BL, preservation, output pointers and control. Long-lived caller x19..x21 survive; leaf uses x0..x12 only. There are no hidden runtime permutations, digit extractions, pointer arithmetic or loops. One fingerprint call per transition writes U directly into the already charged hash final-block buffer.

Polynomial keys, powers and digits are no longer used. The following section bounds all remaining exact-key and hash preparation directly.

## Direct exact-key and hash preparation:512

Polynomial digits are no longer computed. Given the updated group state's b pointer and its exact-key destination, copy its36 limbs in the established reversed-per-coefficient order. For numerical B_i use physical offset48*floor(i/6)+8*(5-(i mod6)), i=0..35; store at key offset8*i. Fully unroll36 pairs LDR x3,[b,#constant]; STR x3,[key,#constant]. Each pair costs2, total72. All source/destination offsets are<=280, initialized and aligned. The sign helper writes its64-bit0/1 at key offset288 within its own128 allowance; it is not overwritten here. Key padding was initialized once in paid setup. Add RET and budget up to55 for argument moves, BL and pointers: charge128 for this full leaf interface. It uses only caller-saved registers and has no internal call or frame. The previous M=289-byte integer,39 digits and EXTR buffer do not exist in this refinement; removing their work does not change the exact key or tabulation digits.

Keep the existing256 allowance for block/state/label operations: four STP xzr,xzr clear the64-byte final block; tabulation's already charged five REV/STR write U directly to its first40 bytes; two LDP/STP pairs copy the cached32-byte prefix h to mutable state (4). After full compression, extract j=((state[6]>>24)<<32)|REV32(state[7]) with two LDR,LSR,REV,ORR (5). This is the low40 bits of the big-endian interpretation of standard little-endian BLAKE3 digest words. The13 listed operations plus pointers and dispatch fit256 with ample slack. Block length40 and flags10 are inside the compression kernel below, not free.

Retain128 further for remaining hash wrapper argument formation, caller preservation/spills and scheduling. Sum128+256+128=512. No allocation or library memcpy occurs per transition; buffers and frames are cold initialized. The table compressor has its own1024 including132 interface allowance; sign has128; exact lookup has its own explicit bounds; the main outer work also retains1024 spills/control. All these allowances are retained, not substituted for each other. This512 replaces the obsolete1024 preparation budget and removes exactly512 instructions per transition. Actual key stores, U writes, block zeroing, prefix copy and extraction remain paid.

## Standard full7-round BLAKE3:1536 including512 preparation

On four distinct32-bit state registers a,b,c,d, use this16-instruction G schedule (all shifts are rotations and all arithmetic is native modulo2^32):

 ADD a,a,b; LDR temp,[block,#x]; ADD a,a,temp; EOR d,d,a; ROR d,d,#16;
 ADD c,c,d; EOR b,b,c; ROR b,b,#12;
 ADD a,a,b; LDR temp,[block,#y]; ADD a,a,temp; EOR d,d,a; ROR d,d,#8;
 ADD c,c,d; EOR b,b,c; ROR b,b,#7.

The source offsets x,y are compile-time standard BLAKE3 message-word positions. Source operands use W registers, the block base an X register. There are56 G (7 rounds*8), costing896, including112 message loads. BLAKE3's fixed message permutation and G column/diagonal choices are fully unrolled: offsets0..60 and register choices are assembly constants, not runtime permutations. No SIMD or crypto extension is used. The word schedule is the standard BLAKE3 permutation, not BLAKE2 SIGMA.

Use W0..W15 for16 working values,W16 for one message word,X17 block,X19 output h,X20 IV address. Preserve x19/x20 with STP/LDP; this is a leaf, no x18, no internal call or state spills. Initialize v0..7 by8 LDR from the cached chaining value, v8..11 by4 LDR from IV. Four MOV set v12=0,v13=0,v14=40,v15=10 for the final WALK block. These16 instructions set the complete initial state directly; there is no cumulative-byte-counter XOR or final-bit inversion from BLAKE2.

For32-byte output, only the first8 output words are consumed: EOR v_i,v_i,v_(i+8), then STR for each i,16 instructions. Computing unused output words8..15 is unnecessary for the identical32-byte BLAKE3 hash result. They do not feed the first8. Thus fixed state/data work is32, inside128 with96 left for stack pair,pointer moves,IV address materialization,entry/exit and RET. The standard root output counter is0. Output words serialize little-endian, so the existing label extraction from words6,7 is correct. Total896+128+512=1536; the512 exact-key/block/state/label/wrapper preparation remains fully charged separately from this core.

The64-byte prefix is compressed once from IV in cold setup with counter0,len64,flags1. Its first8 output words become the cached chaining value. Final payload40 is padded by24 zero bytes, len40,flags10. TABLE54 and SCAL56 single-block root hashes use flags11. All messages are below1024 bytes, so there are no parent nodes or nonzero chunk indices. No BLAKE2 parameter state or finalized-prefix digest is reused. Cold constructions retain the full2^90 allowance. These are analytical scalar schedules, not a measured A64 trace or hash security proof.

## Borrow-chain lower_bound, exact equality, coefficients and counters:3968

The stored key remains37 native64-bit words, B_0 through B_35 and sign. Its injectivity is proved above. For sorting only, interpret word0 as the least significant limb of a2368-bit unsigned integer. Thus the order differs from the predecessor's lexicographic order but equality is unchanged. Keys do not include coefficients A,B. The initial record and subsequent records are inserted in this new order.

To compare saved < current, point x0 and x5 at the two key bases. Fully unroll18 pairs. Load operands with:

    LDP x1,x2,[x0],#16
    LDP x3,x4,[x5],#16

On the first pair execute SUBS xzr,x1,x3 then SBCS xzr,x2,x4. On every later pair execute SBCS for both limbs. Finish with two LDR of the37th word and SBCS xzr,x1,x3. LDP/LDR and post-index address updates preserve flags. The final C flag is1 iff the full saved integer is >= current; C=0 means strictly less. Zero flag is not a multiword equality result and is never used as one. Loads are36 paired plus2 scalar, subtraction instructions37:75 total. Offsets0..288 are aligned and in initialized buffers; operands and address registers are distinct. The padding word is excluded. AArch64 SUBS/SBCS carry/borrow semantics implement standard radix2^64 subtraction. The chain has no data-dependent early exit, SIMD, wide primitive or missing carry propagation.

lookup.S gives the complete leaf, including lower_bound control and final equality. Its exact register assignment uses x0 pointer array,x1 count,x2 current key,x3 lo,x4 hi,x5 mid,x6 saved-key cursor,x7 current cursor,x8..x11 loaded data. It uses no x18 or callee-saved register, no calls, frame or spills. Both key cursors are reloaded for each comparison. SUBS/SBCS discard differences to xzr; the loop consumes only final carry, never final Z. LDP writeback bases are distinct from operand registers. Input count and array base survive intact until return. The caller's long-lived x19..x21 are untouched.

Each successful loop iteration executes CMP,B.hs,SUB,LSR,ADD,LDR,ADD,MOV (8 control/address instructions), the75-instruction borrow chain, then B.lo and either MOV,B or ADD,B (3 instructions). Thus at most86 instructions per iteration; charge96. There are at most7 for n<=127. This count includes both taken and untaken branch instructions as executed scalar work; no speculative or cycle model is used.

Outside these iterations pay two initial MOV, final loop CMP/B.hs, post-loop CMP/B.hs, three cursor instructions,112-instruction equality at worst and three return instructions:124. A miss before equality executes fewer. Retain64 additional scalar instructions for caller argument moves, BL, result handling and any caller preservation. With no internal call and explicit leaf registers, no ABI preservation is hidden. Charge224 for all fixed lookup work (124+64=188,36 spare). The main caller retains its existing separate256 dispatch and1024 spills/control allowances as well.

Total search <=7*96+224=896, versus1440 in the first borrow-order revision and1568 in the original BLAKE2s package. All other allowances are unchanged:512 both scalar updates,256 index/address,512 counters/checkpoint control,512 record copying,256 dispatch,1024 spills/control. The sum896+512+256+512+512+256+1024=3968. Full outer work is reduced by544 from4512. The update and record representation do not change. At most127 sorted insertions remain paid within fixed2^90. Complete944-byte record saves need59 LDP/STP pairs=118 instructions, within512.

The function returns an insertion position even when no key is equal, including an empty pointer array. On equality it returns exactly the unique saved state. Standard lower_bound loop invariants prove no equality is lost. All ordering/equality decisions use group-state keys only, so the finite success argument and uniform-denominator proof retain exactly their earlier conditions. There is no new hash collision assumption.

## Reference-call interface:64, explicit complete schedule

Keep current/next944-byte state-record pointers in x19/x20 and selected jump
record pointer in x21. They are initialized in charged setup; selecting x21 and
j*640 addressing belong to the256 index/address allowance. At the reference call:

    MOV x0,x20
    MOV x1,x19
    MOV x2,x21
    BL blst_fp12_mul
    MOV x9,x19
    MOV x19,x20
    MOV x20,x9

These seven scalar instructions pass output/current/jump pointers and swap
buffers. The reference multiplication writes its576-byte result directly to
next storage. It is the pinned ordinary ABI function and preserves x19..x21;
its own return, spills and memory work are inside the once-counted reference
operation. Coefficient updates write next record's64 coefficient bytes directly
within their512 allowance. Other helpers preserve these callee-saved pointers
or pay their own saved-register traffic in their explicit allowances. Their
argument preparations remain in the unchanged named budgets above.

The complete transition loop's counters, flags and dispatch already have1024
outer control instructions. No extra per-step wrapper call, tower dispatcher,
state copy or scratch allocation is introduced by this reference call. Charge64
rather than the exact seven instructions, leaving57 for any additional interface
moves, pointer materialization or caller scheduling. Function-entry preservation
of long-lived registers is fixed setup, paid in2^90; it does not recur each step.
This tighter explicit interface, rather than an assumed free call, is the new
change from the previous384 allowance. No callee internals are charged twice.

Total nonreference cost1024+1536+128+3968+64=6720 per transition. The fixed
setup/finish budget, peak memory, exact total work and sole heuristic are fully
specified in proof.md and probability.md. The optional Python fixture checks
finite hash semantics only and has not been executed on the host.

Arm Instruction Set Reference Guide D2.83 documents scalar EXTR and its64-bit
immediate range0..63: https://documentation-service.arm.com/static/6245c734b059dc5ff9a8bdab

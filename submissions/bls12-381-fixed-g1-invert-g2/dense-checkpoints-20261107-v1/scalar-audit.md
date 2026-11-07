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
coefficient tests, each six LDR, five ORR, one CBNZ: at most102 instructions.
Keep the selected six limbs, load six H limbs, compare by SUBS/five SBCS and
CSET CS:13. All-zero case outputs zero. Entry pointers, optional BL, 64-bit sign
store, return and control fit32. Total<=117; charge128. Only caller-saved
registers and immediate offsets below288/48 are needed. b serialization is
charged in hash preparation, not treated as free. Canonical Montgomery residues
are unique modulo p, so tuple and field equality are consistent.

## Exact simple-tabulation lookup:1024 including serialization/interface

The exact37-word group key is unchanged; prefix-lookup.S replaces its search procedure. The tabulation compressor consumes the same37-word injective key, using72 unsigned32-bit halves of the36 b words plus one0/1 sign digit. On the specified little-endian AArch64 machine, LDR W directly reads low then high halves at fixed offsets0,4,...,284. LDR X reads the sign at288. These digits together uniquely determine the key; field parameters and sign convention are unchanged. Tables contain five native64-bit words per entry, are independently initialized as described in proof.md, and are aligned at least8 bytes. Entry stride40 preserves that alignment.

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

## Direct exact-key and hash preparation:224

The exact key-copy leaf in hot-loop.S executes36 fully unrolled LDR/STR pairs=72 instructions,RET1 and caller ADD,ADD,BL3:76. The source offset for B_i is48*floor(i/6)+8*(5-(i mod6)); output is8*i, i=0..35. All offsets0..280 are aligned and initialized. Sign is separately written at key offset288 within sign128; prefix at296 is separately written below. Charge96, leaving20 for any extra pointer/interface selection. The leaf uses only caller-saved x0,x1,x3, no stack or internal calls.

Block/state/label/prefix work costs16 displayed scalar instructions:4 STP zeroing64 bytes;4 LDP/STP copying32-byte cached CV;6 LDR,LDR,LSR,REV,ORR,STR extracting/storing the40-bit label;2 LDR/STR caching the64-bit U prefix at record936. Charge64 with48 spare for pointers/control. Tabulation's five REV/STR output writes are already charged in1024 and write U directly into the cleared block. No extra serialization or memcpy is assumed free. Stack frame is cold allocated and aligned, all offsets<=936 encodable and valid. The prefix is cached after tabulation and before any next transition; BLAKE3 reads the same initialized40 payload bytes.

Retain64 for remaining hash-wrapper argument formation,BL,scheduling and any additional caller work. The hot loop's three BLAKE3 call instructions fit here. Persistent x19..x28 survive every helper by explicit contracts; no ordinary spills are needed, and residual outer control retains128 anyway. Sum96+64+64=224. This covers every actual key,block,CV,label,prefix and wrapper operation. No obsolete polynomial M/digit buffer exists. Sign128,table1024,lookup241,reference64 and cold initialization remain distinct charges.

## Standard full7-round BLAKE3:1248 including224 preparation

On four distinct32-bit state registers a,b,c,d, use this16-instruction G schedule (all shifts are rotations and all arithmetic is native modulo2^32):

 ADD a,a,b; LDR temp,[block,#x]; ADD a,a,temp; EOR d,d,a; ROR d,d,#16;
 ADD c,c,d; EOR b,b,c; ROR b,b,#12;
 ADD a,a,b; LDR temp,[block,#y]; ADD a,a,temp; EOR d,d,a; ROR d,d,#8;
 ADD c,c,d; EOR b,b,c; ROR b,b,#7.

The source offsets x,y are compile-time standard BLAKE3 message-word positions. Source operands use W registers, the block base an X register. There are56 G (7 rounds*8), costing896, including112 message loads. BLAKE3's fixed message permutation and G column/diagonal choices are fully unrolled: offsets0..60 and register choices are assembly constants, not runtime permutations. No SIMD or crypto extension is used. The word schedule is the standard BLAKE3 permutation, not BLAKE2 SIGMA.

Use W0..W15 for16 working values,W16 for one message word,X17 block,X19 output h,X20 IV address. Preserve x19/x20 with STP/LDP; this is a leaf, no x18, no internal call or state spills. Initialize v0..7 by8 LDR from the cached chaining value, v8..11 by4 LDR from IV. Four MOV set v12=0,v13=0,v14=40,v15=10 for the final WALK block. These16 instructions set the complete initial state directly; there is no cumulative-byte-counter XOR or final-bit inversion from BLAKE2.

For32-byte output, only the first8 output words are consumed: EOR v_i,v_i,v_(i+8), then STR for each i,16 instructions. Computing unused output words8..15 is unnecessary for the identical32-byte BLAKE3 hash result. They do not feed the first8. Thus fixed state/data work is32, inside128 with96 left for stack pair,pointer moves,IV address materialization,entry/exit and RET. The standard root output counter is0. Output words serialize little-endian, so the existing label extraction from words6,7 is correct. Total896+128+224=1248; the224 exact-key/block/state/label/wrapper preparation remains fully charged separately from this core.

The64-byte prefix is compressed once from IV in cold setup with counter0,len64,flags1. Its first8 output words become the cached chaining value. Final payload40 is padded by24 zero bytes, len40,flags10. TABLE54 and SCAL56 single-block root hashes use flags11. All messages are below1024 bytes, so there are no parent nodes or nonzero chunk indices. No BLAKE2 parameter state or finalized-prefix digest is reused. Cold constructions retain the full2^90 allowance. These are analytical scalar schedules, not a measured A64 trace or hash security proof.

## Prefix-budget lookup and explicit outer work:721

prefix-lookup.S is a complete caller-saved leaf. Input x0 pointer array,x1 count,x2 full current record,x3 remaining-budget pointer. It retains current prefix in x12,budget pointer x13,lo x3,hi x4,mid x5,cursors x6/x7,limb data x8..x11,budget limbs x14/x15,zero test x16. No callee-saved register,stack frame or internal call exists. It returns insertion/match index x0 and flag x1=0 miss,1 full equality,2 budget failure. The caller checks2 before accepting1. It scans the whole equal-prefix run, so ties cannot hide a matching key. A miss index after the equal-prefix run is a valid prefix-order insertion position.

Each ordinary lower_bound iteration executes at most11 scalar instructions: CMP,B.hs,SUB,LSR,ADD,LDR record,LDR prefix,CMP,B.lo,one update,one branch. Charge16 per iteration, at most10. Entry4 instructions, terminal bounds/prefix miss6,last loop test2,return3,caller5 and flag dispatch3 total23, below the fixed80 allowance with57 spare. Thus ordinary search <=10*16+80=240. All equal-prefix candidates have their extra scan/body work charged separately below, rather than silently included in this constant bound.

Per full comparison: scan/bounds/prefix work6; budget LDP,ORR,CBZ,SUBS,SBC,STP6; two cursor ADD2; full37-word equality112; next-candidate ADD/branch2 or final return3. Even including additional bounds/argument/control, this is below256. Budget is128-bit, nonzero tested before subtracting1 with carry. Each complete or false comparison decrements once; before the next over-cap comparison the leaf returns2 without doing it. This final zero-budget check and dispatch fit the fixed80 and64 outer dispatch allowances. Cold setup initializes B_cap=floor(K/256), so the total full-comparison cost across the entire attack is<=256B_cap<=K. Charge one scalar instruction per transition, regardless of how comparisons are distributed. Probability.md bounds cap exhaustion by2^-43; no claim of worst-case absence of prefix collisions is made.

The prefix is exactly the native64-bit LDR of U bytes0..7, stored at record offset936 by two new common-path instructions LDR x9,[sp,#128];STR x9,[x19,#936]. This is a bijective bit permutation of the first64 uniform U bits; fixed unequal keys have prefix-collision probability2^-64. The extra two operations fit the revised224 preparation budget with large spare allowance. Record size remains944; what was unused padding is now initialized prefix. Full equality compares only36 b words and sign at928, excluding prefix. Checkpoint copying saves the full record including the cached prefix. Cold sorted insertion uses prefix order; ties are scanned exactly.

The persistent register/cold ABI and explicit scalar/key-copy leaves in hot-loop.S are retained. All helpers preserve x19..x28, and BLAKE3 explicitly saves/restores its temporary x19/x20. Cold saves LR at48,callee state192..271,remaining-budget272..287 in a >=288-byte aligned frame; handlers restore LR/frame or resume loop by branch. No return through a clobbered helper LR is assumed.

Both scalar updates:128 including interface. fapi_scalar_update moves three pointers to x13..x15 and materializes r at x16 (5 instructions). For each residue4 LDP load current/jump limbs; ADDS,ADCS,ADCS,ADC add;2 LDP load r;SUBS/three SBCS compute candidate subtraction to x9..x12 without destroying the unreduced sum;4 CSEL HS choose candidate/original;2 STP store. This is20 instructions per residue,40 total. RET gives46 internal, and MOV,MOV,MOV,BL at caller give4:50. Charge128 with78 spare. r<2^255 ensures two canonical residues sum without a fifth limb; final subtraction carry determines whether sum>=r. All sources/outputs are disjoint current/alternate records; full offsets576..632 are aligned. No calls,stack,wide unit primitive or callee-saved register occurs in this helper.

Index/address:32. LDR cached40-bit label,LSL by9,ADD with label shifted7,ADD jump base form table+j*640 in4 scalar instructions. The complete table is below2^50, as are all allocated record addresses with a chosen low base. Label store/extraction remain paid in hash preparation. No hidden variable-width addressing or lookup primitive occurs.

Counters/checkpoint control:128. Each iteration's top uses LDP K,SUBS/SBCS for unsigned128-bit t>=K,B.hs termination,ADDS/ADC increment:6 instructions. Bottom CMP low marker,CCMP high marker with immediate flags0 when low differs,B.eq cold handler use3. If low differs,Z=0 by flags0; if low matches,the high comparison gives exact equality. K<2^128, so no counter wrap. The next marker is initialized to s and only advanced on checkpoint; after t=1022s it advances to1023s>K. Thus there are at most1022 insertions after the initial record. Their arithmetic and path work are fixed-budget costs, not omitted common instructions.

Copying:zero per-transition charge, with every actual copy fully paid cold. Only the initial record and1022 checkpoint events save a record. A full944-byte copy is59 LDP/STP pairs=118 instructions; at most1023*118<2^14. Pointer insertion shifts at most1022*1023 pointers, and all copying/insertion/marker/controller work remains below2^30 scalar instructions, already inside the fixed2^90 allowance. There is no arbitrary new memory snapshot at every step; the pinned multiplication writes576 output bytes directly, scalar helper writes64 coefficient bytes, sign/key-copy write296 meaningful key bytes, and prefix caching writes the remaining8 bytes in charged preparation. Hash/preparation and reference work retain their own explicit costs.

Dispatch:64 includes budget-flag2 CMP/B.eq,found CBNZ and result selection. Remaining outer control/spills:128 includes final top branch and any small caller selection. No ordinary spills are necessary with the explicit caller-saved leaves/preserved persistent state. Reference interface64,preparation224,sign128,table1024 and hash core each retain their own full interfaces. The updated hot loop adds only prefix cache2,one budget-pointer argument and two budget-result dispatch instructions; all fit their stated allowances. No copied or arithmetic result is charged twice.

Outer total241+128+32+128+64+128=721. Full nonreference3185=1024+1248+128+721+64. Cold initialization,failures,copies,insertions,verification and G2 output remain within2^90, peak below2^50.

## Reference-call interface:64, explicit complete schedule

Keep current/next944-byte state-record pointers in x19/x20 and selected jump
record pointer in x21. They are initialized in charged setup; selecting x21 and
j*640 addressing belong to the32 index/address allowance. At the reference call:

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
within their128 allowance. Other helpers preserve these callee-saved pointers
or pay their own saved-register traffic in their explicit allowances. Their
argument preparations remain in the unchanged named budgets above.

The complete transition loop's counters, flags and dispatch have the explicit
separately charged bounds above. No extra per-step wrapper call, tower dispatcher,
state copy or scratch allocation is introduced by this reference call. Charge64
rather than the exact seven instructions, leaving57 for any additional interface
moves, pointer materialization or caller scheduling. Function-entry preservation
of long-lived registers is fixed setup, paid in2^90; it does not recur each step.
This tighter explicit interface, rather than an assumed free call, is the new
change from the previous384 allowance. No callee internals are charged twice.

Total nonreference cost1024+1248+128+721+64=3185 per transition. The fixed
setup/finish budget, peak memory, exact total work and sole heuristic are fully
specified in proof.md and probability.md. No executable experiment is declared; these are analytic schedules, not measured traces.

Cold ABI detail: the controller reserves at least288 aligned bytes, saves original x30 at offset48 and x19..x28 at192..271, then branches into fapi_hot_loop. This label does not return with the last helper's LR. Cold terminal handlers restore the saved original LR/callee registers and frame before returning to the caller. Checkpoint handler branches back into the loop after its bounded work. These once-only entry/exit and at-most1022 checkpoint-handler instructions are in the unchanged fixed2^90 budget.

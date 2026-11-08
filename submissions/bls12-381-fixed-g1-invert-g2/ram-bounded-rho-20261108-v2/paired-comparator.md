# Full-state comparison using scalar paired loads

This comparator preserves the predecessor's unsigned lexicographic order on all72
canonical64-bit words of a576-byte GT state. It is not a hash comparison or a
norm-one compressed key. It does not change the walk, success event, coefficients,
checkpoint schedule or the unconditional correctness of a detected equality.

Copy the current and saved record bases into two temporary pointer registers.
For each of36 pairs, fully unroll the following six scalar AArch64 instructions:

    LDP x1, x2, [x0], #16
    LDP x3, x4, [x5], #16
    CMP x1, x3
    B.ne mismatch
    CMP x2, x4
    B.ne mismatch

x0 and x5 are distinct temporary base registers and do not overlap their loaded
registers. Both loads are64-bit scalar load pairs, with ordinary post-index16
writeback. Each LDP is one scalar instruction under the protected instruction
metric; memory words, pointers and comparison work are not omitted. Post-indexing
avoids the signed-offset limitation of LDP for the latter part of the576-byte
object. Each pair accesses initialized canonical words within the allocated
record. Record bases are at least16-byte aligned because slots have640-byte
stride and an aligned base. There is no cross-record or unaligned access.

CMP/B.ne on the first words precedes CMP/B.ne on the second words, so the first
mismatching word is exactly the predecessor's lexicographic discriminator. No
flags are modified between a CMP and its branch. At mismatch the condition flags
encode unsigned less/greater; at fall-through all72 words agree. A shared return
label may form the trichotomy result in a few scalar instructions. These and
caller lookup, lo/hi/m arithmetic, interval update, BL/RET, base resets and any
saved registers retain the predecessor's full96-instruction allowance. No
additional loop counter, dynamic indexing or per-pair control exists. Temporary
pointers are reset before every comparison; immutable record pointers and scalar
coefficients remain intact.

Worst-case word-comparison work is36*6=216 rather than72*4=288. Plus96 outer
instructions gives312. Charge320 per comparator/search iteration, with8 spare.
Six iterations cost1920. All other outer allowances remain512 for both scalar
updates,256 index/address,512 multiword counters/checkpoint control,512 full
record copying,256 equality/dispatch,1024 spills/control. Thus outer work is
1920+512+256+512+512+256+1024=4992. Sorted pointer insertion remains in the fixed
2^90 budget, not hidden as a free transfer.

Reference multiplication is still one H operation. Other work per step is
2208 compressor +3296 SHA-512 +128 sign +4992 outer +384 interface =11008.
Complete W<=K+(11008*K+2^90)/20000=(31008*K+2^90)/20000. No reference price,
hardware timing, SIMD, free memory oracle or changed probability argument enters.
Peak memory remains below2^35 bytes; shorter comparator code fits inherited
space. No executable experiment is declared; the schedule is static evidence.

# Exact 22-instruction SHA-512 round, charged at 24

This is an analytical scalar AArch64 schedule, not an executed candidate program
or a measured benchmark. All arithmetic wraps at64 bits. EOR's rotated-register
operand is an ordinary scalar A64 operation. There are no SIMD or hash-extension
instructions. Round constants and expanded words are loaded and charged.

Assign logical state(a,b,c,d,e,f,g,h) to eight distinct registers. At each round,
the assembler uses the appropriate physical register names for these logical
names; their mapping rotates statically between rounds. Temporary registers
t0,t1,t2,t3 are distinct from state and the two base pointers wp,kp.

The following has exactly22 scalar instructions:

    ROR t0, e, #14
    EOR t0, t0, e, ROR #18
    EOR t0, t0, e, ROR #41
    EOR t1, f, g
    AND t1, t1, e
    EOR t1, t1, g
    ADD t2, h, t0
    ADD t2, t2, t1
    LDR t3, [wp, #word_offset]
    ADD t2, t2, t3
    LDR t3, [kp, #constant_offset]
    ADD t2, t2, t3
    ROR t0, a, #28
    EOR t0, t0, a, ROR #34
    EOR t0, t0, a, ROR #39
    ORR t1, a, b
    AND t1, t1, c
    AND t3, a, b
    ORR t1, t1, t3
    ADD t0, t0, t1
    ADD d, d, t2
    ADD h, t2, t0

Here t2=T1=h+Sigma1(e)+Ch(e,f,g)+W[t]+K[t] and
t0=T2=Sigma0(a)+Maj(a,b,c). The last two instructions replace old d by
new e and old h by new a. The physical registers holding old a,b,c,e,f,g
are not written. Rename the logical tuple for the next round to
(old h,old a,old b,old c,old d,old e,old f,old g). This implements exactly
the SHA-512 state update. No register-to-register copy is required.
The majority identity (a AND b) OR (c AND (a OR b)) equals bitwise majority.
The choice identity g XOR (e AND (f XOR g)) is the specified choice function.
All temporary values are dead when overwritten. The two final state writes
happen only after every use of their old values. There are80 rounds and
80 is a multiple of8, restoring the original physical state mapping before
feed-forward.

Eight state registers plus four temporaries and two pointers require14 GPRs.
Map these into x0..x13; all are caller-saved even where x18 is reserved.
Pointers and temporaries are reused between expansion and rounds. The maximum
word and constant offset is632, an encodable nonnegative scaled64-bit LDR
offset. The code is fully unrolled; offsets and logical-register mappings are
compile-time choices. No runtime counter, branch, table of register mappings,
pointer increment, or callee-saved-register spill is required per round.
No external call occurs within the block. Charge24 instructions per round,
leaving two instructions beyond this complete22-instruction schedule.

The predecessor's expansion bound16 per generated word, fixed overhead128
and preparation2048 are retained in full. Block cost is therefore
80*24+64*16+128+2048=5120 scalar instructions. The routine computes exactly
the same standard SHA-512 digest as its predecessor. The hash idealization
and every probability bound are unchanged; savings concern actual charged
scalar computation, not a free oracle.

The fixed128 covers48 input parsing/stores,8 state initializations and24
feed-forward instructions, leaving48 for entry/exit, stack allocation and
pointer setup. The5120 also retains the entire2048 preparation allowance.
There is no hidden saved-register traffic after using the14-register mapping.
The schedule frame is640 bytes; its constants can be embedded within a few
KiB. Such fixed code/data fit well within the inherited2^50 peak bound.

References: FIPS180-4 sections4.1.3 and6.4 specify these operations:
https://csrc.nist.gov/pubs/fips/180-4/upd1/final .
The scalar rotated-register EOR encoding is documented by Arm:
https://developer.arm.com/documentation/ddi0602/latest/Base-Instructions/EOR--shifted-register---Exclusive-OR--shifted-register-- .

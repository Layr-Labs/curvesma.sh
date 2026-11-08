# Exact 14-instruction SHA-512 schedule expansion

For each i=16,...,79, fully unroll this scalar A64 sequence. The symbolic
indices below are replaced at assembly time by immediate byte offsets.
The schedule array contains80 eight-byte words. All offsets0..632 fit the
ordinary scaled unsigned LDR/STR immediate encoding. There is no runtime
index multiplication, loop counter, branch, pointer increment or function call.

    LDR t0, [sp, #8*(i-15)]
    ROR t1, t0, #1
    EOR t1, t1, t0, ROR #8
    EOR t1, t1, t0, LSR #7
    LDR t0, [sp, #8*(i-2)]
    ROR t2, t0, #19
    EOR t2, t2, t0, ROR #61
    EOR t2, t2, t0, LSR #6
    LDR t0, [sp, #8*(i-16)]
    ADD t1, t1, t0
    LDR t0, [sp, #8*(i-7)]
    ADD t1, t1, t0
    ADD t1, t1, t2
    STR t1, [sp, #8*i]

The sequence has14 scalar instructions:4 loads,6 sigma operations,3 adds,
1 store. Modulo2^64 it implements W_i=sigma0(W_(i-15))+sigma1(W_(i-2))+
W_(i-16)+W_(i-7), exactly FIPS180-4. Only3 temporaries and sp are live.
All predecessor words have already been written; writing W_i cannot alias any
source. Expansion finishes before the eight working-state registers are loaded,
so register pressure is lower than during rounds. Initial16-word parsing is
included in the fixed128 overhead. There is no additional stack allocation per
word; the fixed640-byte frame and its entry/exit are already charged there.

Combined with the22-instruction round-schedule.md, the full state-dependent
block plus preparation costs80*22+64*14+128+512=3296 instructions. All input,
output, padding, prefix-state copying and final digest extraction retain the
separate complete512 preparation allowance in preparation-audit.md. The scope
is analytical scalar instruction accounting. No executable experiment is declared.

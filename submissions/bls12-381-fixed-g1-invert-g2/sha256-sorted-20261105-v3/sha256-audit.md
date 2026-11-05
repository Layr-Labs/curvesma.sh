# Complete scalar SHA-256 audit

The algorithm definition,word size,block size,padding and round functions are
FIPS180-4 sections4.1.2,5.1.1,5.2.1,5.3.3 and6.2.2:
https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.180-4.pdf .
This is an analytic implementation audit, not a measured instruction trace.
The full hash construction is in prefix-audit.md. All arithmetic below uses
ordinary AArch64 W-register instructions modulo2^32, supported by the stated
scalar instruction model. No SIMD or SHA acceleration is used.

Unroll all48 expanded schedule words. Each uses two small-sigma functions:
each function has two ROR,one LSR and two EOR, five instructions. Four LDR of
32-bit schedule words,three ADD and one STR complete the recurrence. Thus18
instructions per expanded word,charged20. All offsets in the256-byte schedule
are fixed immediates; no loop counter or helper call is required.

Unroll all64 rounds. Each round's two big sigmas use three ROR and two EOR,
ten instructions total. Choice uses three logical operations. Majority can
use five logical operations. T1 uses four additions,T2 one,and updated d+T1
and T1+T2 two more. Two loads supply W_t and K_t; at most eight register moves
rotate the working variables. Total10+3+5+4+1+2+2+8=35 instructions. Charge40
per round. Immediate offsets fit64-bit base pointers with fixed32-bit LDR.
Eight working words,two T words,two reusable sigma/logic temporaries and a few
pointers fit ordinary scalar registers. Schedule expansion finishes before
round working words are loaded, so the phases do not compete for registers.
Any preserved-register saves/restores belong to the following block overhead.

Initial16 word loads and REV cost32. Eight chaining-word loads cost8. Final
feed-forward eight loads,eight adds and eight stores cost24. Permit64 further
instructions for prologue/epilogue,frame,pointers,return and fixed control;
the block overhead is128. The unrolled compression bound is64*40+48*20+128=3648.

Reserve2048 further scalar instructions for block preparation, copying the32-byte
precomputed chaining state,label extraction,and the compressor's aligned input
buffer preparation described in compressor-audit.md. That preparation is the
same37-limb representation of the289-byte norm-one encoding. It costs below512;
all words are initialized,including the high zero padding. Formatting the64-byte
SHA block costs below512 more:40-byte compressor output,0x80,zeros and832-bit
length field plus state copy. Loads/stores,byte reversals and pointer updates
are included. There is no hidden whole-state copy or arbitrary-width primitive.
Total3648+2048=5696; charge6144, leaving448 additional instructions of slack.

Use the unchanged45312 dense field-product bound,4096 polynomial-compressor
bound,2048 norm-one encoding/sign bound and8192 sorted-checkpoint/scalar/control
bound. The component sum is65792;charge66048 per complete rho transition.
Every failed transition and actual hash is included. All pointer insertions,
precomputed prefix,INIT,all two-block initialization hashes,table construction,
subgroup/input checks,pairings,output and verification remain in the cold2^90
fixed-work allowance, whose explicit components remain far below that value.
Memory changes by fewer than4096 bytes and remains below2^35.

With the shorter horizon and sorted checkpoints,the exact public inequality
(K*66048+2^90)^50 <20000^50*2^6469 gives score6469/50=129.38. The success proof
still subtracts all finite losses and declares0.39,under the one SHA-256
random-function idealization. No production inversion or measured score is
claimed. This variant derives from published predecessor414df432-b6c8-422b-be3d-42d0014a4a80,
the published shorter-horizon submission and earlier sorted-checkpoint work identified
in horizon-audit.md,and the general unrolling idea in failed public note
ad5bd06a-9732-435e-a4a3-c387b7bfcaeb. That note was inspected as a source for
unrolling SHA;the new SHA-256 message and32-bit schedule are audited here.
The failed public unrolling note remains unpublished; that dependency will be added if publication makes it eligible.

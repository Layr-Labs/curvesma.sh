# Exact SHA-512 prefix and scalar implementation

The hash is standard SHA-512 with the tags, seed, full padding and message lengths in proof.md. WALK||seed||zero48 is an instance-fixed 128-byte prefix, compressed once in charged setup. A label hashes the remaining 40-byte compressor output with its standard padding and the original 1344-bit total length. Extract 40 digest bits as a scalar unsigned integer for the 2^40-entry table. INIT and SCAL computations retain their original full SHA semantics.

This implementation uses only scalar A64 integer instructions, including EOR with an optionally rotated or shifted register operand. Arm documents that operand in its [Instruction Set Reference Guide](https://documentation-service.arm.com/static/6245c734b059dc5ff9a8bdab). Mathematical hash behavior follows [FIPS 180-4](https://csrc.nist.gov/files/pubs/fips/180-4/final/docs/fips180-4.pdf), sections 4.1.3 and 6.4. This is a static instruction audit, not an empirical production timing.

For example, Sigma1(e) costs three instructions:

    ROR t,e,#14
    EOR t,t,e,ROR #18
    EOR t,t,e,ROR #41

Sigma0 uses the same pattern with rotations 28,34,39. The small sigma0 uses ROR 1, EOR with ROR 8, EOR with LSR 7. The small sigma1 uses ROR 19, EOR with ROR 61, EOR with LSR 6. All are the specified 64-bit functions. The shifted XOR is one scalar instruction, not a crypto extension or SIMD operation.

Expand each of 64 later schedule words with four loads of predecessor words, two three-instruction small sigmas, three ADDs modulo 2^64 and one store. This is 14 instructions; charge 14 per word as detailed in expansion-schedule.md. The schedule is fully unrolled, with immediate offsets in a fixed 640-byte frame. No runtime word-loop counter or dynamic address helper occurs.

One round uses two three-instruction big sigmas; Ch as XOR(f,g), AND(e,result), XOR(g,result), three instructions; Maj as ORR(a,b), AND(c,result), AND(a,b), ORR of these terms, four instructions. Load the schedule word and round constant, two instructions. Form T1 with four ADDs, T2 with one, then new E and A with two more: seven ADDs total. Round arithmetic plus loads is 6+3+4+2+7=22 instructions. Logical names of the eight working registers rotate by a compile-time permutation. New A and E overwrite old H and D only after their last uses; the remaining old values stay in place. No runtime moves are needed for this static renaming. Charge 22 per round, exactly the enumerated schedule; no extra per-round operation is required. The exact 22-instruction schedule is given in round-schedule.md. There are 80 fully unrolled rounds.

Register pressure is bounded: eight state registers, five reusable temporaries, two loaded words and at most three base/output pointers need at most 18 general-purpose registers. Schedule expansion occurs before state rounds and reuses these registers. Caller-saved registers suffice with at most one additional saved pair on an ABI reserving x18. No uncharged per-round stack spill is necessary. Each immediate load offset is below 640 bytes.

Fixed block overhead: parse and store 16 input words with LDR/REV/STR (48 instructions); initialize eight working words (8); feed-forward with eight LDR/ADD/STR triples (24). This totals 80. Prologue, stack frame, saved pair, pointer setup, epilogue and return fit in the remaining 48 of a 128-instruction fixed allowance. No calls occur inside this block function. The old 2048 scalar preparation allowance is retained in full for serialization, prefix-state copying, padded-message construction, compressor input preparation and label extraction. compressor-audit.md details those transfers.

The complete state-dependent block is bounded by 80*22+64*14+128+2048=4832 instructions. The implementation computes the same hash, and keeps the sole disclosed fixed-SHA idealization; none of these savings idealizes hash work as free.

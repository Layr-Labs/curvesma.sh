# Deterministic enumeration control: secp256k1-ecdlp

## Claim and scope

This is a deliberately conservative analytic baseline, not a new cryptanalytic technique, a competitive attack, or a measured production run. It solves every honestly generated full-production input for this track, hence has success probability 1, exceeding 0.39. It uses no heuristics, quantum computation, random coins, nonuniform advice, precomputed answer table, or amortization. The submitted computation upper bound is 2^310 target-operation equivalents, with rational log2 score 310/1. The considerable slack is intentional; this is not a best-known complexity estimate. All quantities below are upper bounds, not measurements.

This submission applies ONLY to secp256k1-ecdlp. There is no auxiliary powers input and no BLS claim.

## Algorithm

Let G be the protected standard generator, r its order, and T the target point (called Q on secp256k1). Start i=0 and R=infinity. While i<r, compare R with T; if equal, return the canonical 32-byte scalar i. Otherwise set R=R+G and i=i+1. The final comparison is an exact equality in the full production group. No truncated scalar search is used.

The invariant is R=[i]P, by induction from infinity and the group law. For the unique x in 0..r-1 with Q=[x]P, iteration i=x matches and returns x. Every input scalar, including zero, is covered. Equality is unique because P has prime order r. No collision or runtime heuristic is used.

## Concrete base-field and point arithmetic

Use 64-bit scalar AArch64 instructions and fixed arrays of seven 64-bit limbs for all base-field arithmetic. The protected secp256k1 modulus has 256 bits, hence at most 381 bits; the seventh limb accommodates carry for addition. Canonical residues are always reduced. An addition computes a+b with carry and conditionally subtracts p. Subtraction computes a-b and conditionally adds p on borrow. Comparisons inspect the limbs from most significant to least. Load/store, carry, branch and address work are included: reserve 128 instructions per limb per pass and four seven-limb passes plus 512 fixed instructions. Thus one addition/subtraction/comparison/copy costs at most A=2^12 instructions. This is generous for straight scalar loads, ADDS/ADCS, SUBS/SBCS, comparisons, selections and branches; it is not a unit-cost wide-integer assumption.

Multiplication is bit-serial double-and-add modulo p: t=0,v=a; for bits j=0..380 of b, conditionally t=addmod(t,v), then v=addmod(v,v). Include leading zero bits. At most 762 modular additions, 381 bit tests and loop overhead, plus setup/copy, fit M=2^24 instructions (762*4096 is less than 2^22). Squaring is a multiplication. For nonzero a, inversion is binary powering a^(p-2), at most 762 multiplications and less than 2^20 instructions of bit scanning, setup and control; I=2^35 bounds it. This uses Fermat's theorem for the protected prime p, not an inversion oracle.

Maintain R in affine coordinates and a separate infinity flag. For R+G: if R is infinity, copy G; if the x coordinates agree and the y coordinates sum to zero, return infinity; otherwise the equal-point branch uses slope (3*x_R^2)/(2*y_R), and the unequal-x branch uses (y_G-y_R)/(x_G-x_R). Then x_new=slope^2-x_R-x_G and y_new=slope*(x_R-x_new)-y_R (in doubling x_G=x_R). On the nonexceptional branches the denominator is nonzero. secp256k1 has a=0. These are the standard chord/tangent group law. Charge one I, 32 M and 32 A plus 2^20 overhead, below 2^36 per addition. This includes testing every exceptional branch. Exact point comparison and scalar increment use fixed arrays; no scalar exceeds 256 bits. Output is canonical hex or null for infinity.

## Complete cost and memory

Charge 2^60 scalar instructions per iteration, which more than covers one affine point operation bounded above by 2^36, exact comparison, incrementing the four-limb scalar counter, branches, memory transfers and control. There are at most r iterations, r<2^256. Charge another 2*2^60 for cold setup, decoding, validation, independent verification and output: computing [r]Q and [x]P by binary scalar multiplication requires fewer than 1024 point operations; all remaining fixed-size work fits this reserve. The complete instruction count is below (r+2)*2^60 < 2^317. Divide by the protected C=1000, which is greater than 2^9. Thus work is less than 2^308 target-operation equivalents. Submit the intentionally looser bound W=2^310 and exact rational score 310/1. This is an upper bound, not a tight estimate or measured value. All arithmetic is counted as R/C, with H=0; there is no double charging or free reference arithmetic.

Cold preprocessing is bounded by 2^60/1000 equivalents and is INCLUDED in total work. No precomputed advice is used; advice is zero bytes. Peak memory is bounded by 2^24 bytes including code, constants, input, output, scratch and stack: fewer than 4096 seven-limb field slots use 229376 bytes, counters and input/output need less than 65536 bytes each, and the iterative fixed-loop code/stack/constants fit within 2^20 bytes. No growing tables, recursion or saved iteration history are used. All addresses fit in 64 bits.

## Evidence and provenance

No experiments, production execution or official score are asserted by this package. This is a deterministic analytical accounting control, intentionally much slower than standard square-root attacks. The target profile and classical-group-work-v2 at e264843c0b59cf0906b84b2e9b94f03a45a00027 define the constants and success event. The standard enumeration and chord/tangent formulas are written explicitly above.

This revision corrects the overly broad wording of our previous enumeration-v1 attempt: its reviewer found the secp argument sound but rejected unrelated BLS assertions in shared boilerplate. This package removes those assertions and the unrelated pairing appendix. It makes only a secp256k1 claim. Earlier failed/unpublished submission IDs were d649a870-1535-4b43-9d6f-19f8b37e656b (provider timeout) and 0df223f5-fcc4-49fd-99a0-0b633fc6e01a (scope rejection). No published prior submission was reused. The existing candidate folders are preserved.

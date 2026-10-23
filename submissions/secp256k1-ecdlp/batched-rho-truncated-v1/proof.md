# secp256k1-ecdlp: inversion-batched parallel double-or-add walk at birthday exponent 1

Status: ready for independent review; no qualification or emitted score.
This is a conservative accounting control under the current `classical-group-work-v2` price, not a best-known attack claim.
The normative inputs and outputs are those of `production-v2`. `claim.json` binds `tracks/secp256k1-ecdlp/target-profile.json` and `tracks/classical-group-work-v2.json` by the SHA-256 of those files. No restricted challenge distribution, external advice, random-oracle premise or quantum operation is used. The historical word-operation price 78037248 is not used.

## Price

`classical-group-work-v2.json` fixes the secp256k1 reference operation at C=1000 scalar AArch64 instruction equivalents. That operation is one batched incomplete affine addition in the 6M normalisation of `tracks/docs/COST_MODEL.md`: 3M for the addition formula and 3M for batch amortisation, with S=M, and no standalone inversion on top of that 3M. The measured field multiplication is M=179 instructions. A reference operation costs 1. Other instruction work costs its count divided by C. A reference operation's internal instructions are not charged a second time. The residual counts below are only the work outside that 6M envelope.

## Walk

Let r be the subgroup order and P the profile generator. Fix B=4096 walks and 16 jump points, all computed from the public input Q and from P.

For j in 16..31 let a_j = j-15 and b_j = 1+(j mod 5), and let R_j = [a_j]P + [b_j]Q. Each coefficient is a public constant at most 16. Walk i, for i in 0..B-1, starts at X_i = [i+1]P + Q with coefficients (a,b)=(i+1, 1), so X_i = [a]P + [b]Q and b=1.

A round advances every walk once. The B additions or doublings are independent, so their denominators are known before any inversion and are inverted by one Montgomery batch. For a walk with affine point (x,y), not infinity, let x_0 be the least nonnegative residue of x modulo 256. For t=0,1,2,3 compute branch = (x_0 + t) mod 32 and the denominator below, and commit to the smallest t whose denominator is nonzero. At most two curve points share an x-coordinate, so at most two of the 32 branches have a zero denominator against the current point. Any three consecutive branches therefore contain a nonzero denominator, and the four-candidate scan finds it.

- If branch < 16 the step is a doubling: denominator 2y, coefficients (a,b) := (2a, 2b) mod r. Incomplete doubling is λ = (3x^2)·(2y)^{-1}, x' = λ^2-2x, y' = λ(x-x')-y. That is one inversion and 4M when squarings are priced as multiplications (two squarings and two further multiplications).
- If branch ≥ 16 the step is an addition of the precomputed affine point R_branch: denominator x-x(R_branch), coefficients (a,b) := (a+a_branch, b+b_branch) mod r. Incomplete addition is one inversion and 3M.

Infinity uses the same rule as x=0. A zero denominator means the current point equals or is opposite the selected jump; the scan commits to another branch instead of repeating that denominator. The output point of the step is written back as an affine point after the batched inverse, and the linear relation X = [a]P + [b]Q is preserved because every jump's coefficients are known.

Distinguished points are the affine results with x ≡ 0 (mod 2^20). They are inserted into a balanced binary search tree keyed by (x,y), storing (a,b). The tree is capped at slots = S/2^16 nodes, where S is the step count below. A duplicate key with b' ≠ b yields

x_sol = (a-a')·(b'-b)^{-1} mod r,

then the attack checks [x_sol]P = Q and returns x_sol on success. A duplicate with the same coefficients is ignored. If the cap is hit, later distinguished points are not stored and the walks still run to the end of the budget. There is no negation map, no GLV endomorphism and no quotient by an automorphism orbit.

## Budget and success

Let T be the least positive integer with T(T-1) ≥ 2r. The integer is

T = 481231938336009023090067544955250113854.

The attack runs R rounds with

R = ceil(T/B) + 2^20 = 117488266195314702902848521718638936,

so the number of walk steps is the exact multiple

S = B·R = 481231938336009023090067544959545081856.

The added 2^20 rounds are 2^20 steps on every walk, B·2^20 = 2^32 steps in total. They pay for the distinguished-point detection lag: after two visits to the same group element, the walk may need a further geometric wait of mean 2^20 before both visits are distinguished. Budgeting 2^20 extra steps on every walk dominates that wait. S(S-1) ≥ 2r holds because S ≥ T.

Heuristic, used only here. Model each new walk point as uniform in the group of order r, in the usual random-mapping sense for a Pollard iteration (Pollard 1978; the distinguished-point accounting of van Oorschot and Wiener 1999). The probability that S recorded points are pairwise distinct is at most exp(-S(S-1)/(2r)). Hence the probability of at least one point collision by step S is at least 1-exp(-S(S-1)/(2r)). The construction gives S(S-1)/(2r) ≥ 1, so that model probability is at least 1-e^{-1} > 0.632. A point collision of two states (a,b) and (a',b') is an equation (b'-b)x ≡ a-a' (mod r). While the walk length is far below r, a doubling-and-add coefficient update of this shape produces b' ≠ b for distinct visit times except on a thin set of carries; those useless matches are why the submitted probability is not 0.632. The same slack covers a biased 32-way partition and an abort of inserts past the storage cap. The cap is 16 times the model mean S/2^20. The submitted success probability is therefore the acceptance floor 0.39, not the model value. Every input, including those that miss, is charged the full step count S. There is no retry, no amplification and no early-abort credit.

## Instruction ledger

Each of the S steps is one reference operation (the batched inverse plus the 3M formula envelope). Three residuals sit outside that envelope and are charged in instructions on every step, whether or not the step actually uses them:

| Residual | Instruction upper bound | Why |
| --- | ---: | --- |
| Denominator scan, coefficient update, branch, distinguished-point test, loop | 192 | Four field subtractions at 24 instructions, a 256-bit double or add modulo r at 64, and 32 instructions of addressing and control |
| Doubling formula beyond the 3M addition envelope | 179 = M | Doubling is 4M plus the batched inverse; the extra multiplication is charged on every step |
| Inversion above the 3M batch allowance | 32 | A Fermat inverse x↦x^{p-2} is at most 256 squarings and 256 multiplications, hence 512M. One inverse serves B=4096 denominators, which is 0.125M = 22.375 instructions. Montgomery's product tree is the 3M already inside the reference operation. The 32-instruction line covers the 22.375 |

Per distinguished-point insert the balanced tree uses at most 256 key comparisons (the cap has about 2^112 nodes and an AVL height is below 1.44·log2(N)+2 < 256) at 32 instructions each, hence 8192 instructions. The attack pays this for every reserved slot, not only for slots that are filled:

slots = floor(S/2^16) = 7343016637207168931428032607414933.

A tail of 2^20 reference operations covers building the 16 jumps (each a double-and-add of a coefficient below 16, charged as a full 256-bit scalar), the B starting points, the final scalar inversion of (b'-b) modulo r by binary exponentiation, and the verification scalar multiplication [x_sol]P. A further 2^20 instructions cover decoding, canonical-point checks and output encoding. No other table is built. Jump points and starting points depend only on P, Q and the fixed coefficients; they are not uncharged advice.

Let H be the reference-operation count and R the residual instruction count:

H = S + 2^20 = 481231938336009023090067544959546130432,

R = 403·S + 8192·slots + 2^20 = 193996625141703637433183479061816612167680.

The total work in target-group-operation equivalents is W = H + R/C with C=1000:

W = 675228563477712660523251024021362742599680 / 1000.

log2(W) = 128.988643… and the submitted time_log2 is that value rounded upward to six decimal places, **128.988644**. The score rational is 128988644/1000000.

Preprocessing is the whole 2^20 reference-operation tail plus the 2^20 instruction tail, which is (2^20·1000 + 2^20)/1000. Its base-2 logarithm rounds upward to **20.001442**. That bucket includes verification as well as setup, so it is an upper bound on preprocessing proper. It is included in H and R, hence in W.

Peak memory is slots records of 192 bytes (two field elements, two scalars, and tree links) plus 2^20 bytes of code and scratch plus B·384 bytes of live walk state:

1409859194343776434834182260626288576 bytes,

whose base-2 logarithm rounds upward to **120.084963**. Nonuniform advice is 0 bytes, reported as log2(max(1,0)) = 0. Distinguished-point storage is computed by the attack and is inside the memory bound.

## Evidence

The collision algebra is the linear relation preserved by doubling and by addition of a known jump. The toy experiment checks that relation on a small additive group: every invertible collision returns the secret, and the step counter meets the same inequality steps·(steps-1) ≥ 2·order. On the shipped input order 1009 the observed recovery rate inside that budget is above 0.39. The toy partitions the scalar rather than an elliptic x-coordinate, so it is evidence only at that scope. It does not certify the production instruction counts. Organizer arithmetic fixtures remain the authority for point encoding. Heuristic premises: the single random-mapping bound above. No other heuristic is used.

The historical full baby-step/giant-step control is about 2^129 reference group operations for success probability 1, priced at a retired word-operation unit. This package does not reuse that unit. It stops a batched walk at one square-root of the group order and prices every step at the current C=1000 reference operation plus the residuals above. Beating this control does not establish novel cryptanalysis or a weakened hardness assumption.

References: Pollard, Mathematics of Computation 32 (1978); van Oorschot and Wiener, Journal of Cryptology 12 (1999); COST_MODEL.md 6M normalisation. The construction and the integers above are self-contained.

## Submitted computation bound

| Component | Count |
| --- | ---: |
| Reference operations H | 481231938336009023090067544959546130432 |
| Residual instructions R | 193996625141703637433183479061817659695104 |
| W = H + R/1000 | 675228563477712660523251024021362742599680/1000 |

time_log2 **128.988644**, preprocessing_log2 **20.001442**, memory_log2_bytes **120.084963**, success probability **0.39**, nonuniform advice **0**.

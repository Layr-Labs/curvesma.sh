# Six-class checkpoint rho with the unit inside the coefficient product

Parent tree `6423c2728dc03d5f2d647c7f22adf1dac0643a66`. Track price C = 2000 reference operations per batched incomplete affine G1 addition (6M, S = M, A = 0). No auxiliary powers. The score is an upper bound on log2(W) for one cold classical attack, W = H + R/C.

## Group action

Let r be the subgroup order and P0 the protected generator. λ is the scalar root of X²+X+1 modulo r with φ(P0) = [λ]P0, where φ(x, y) = (βx, y) and β is the matching root of X²+X+1 modulo p. The other root is −1 minus the first; `reference.py` keeps the pairing that satisfies the identity on P0, not a random choice of square root. Then φ is the group endomorphism [λ], λ³ ≡ 1 (mod r), and the six units {1, λ, λ², −1, −λ, −λ²} are distinct. Every non-identity point has orbit size 6 under multiplication by those units, because uP = vP implies (u−v)P = 0. The identity is its own orbit. The number of non-identity orbits is

N = (r − 1) / 6.

r ≡ 1 (mod 6), so N is an integer. A canonical representative is the element of the orbit whose x-coordinate is the integer minimum of {x, βx, β²x}, and, on that x, whose y is the integer minimum of {y, −y}. β²x = −x − βx, so the three x-coordinates cost one base-field multiplication.

## Discrete-logarithm state

T = [k]P0. A state is a point S together with four residues (a0, a1, b0, b1) modulo r:

S = [a0 + a1·λ] T + [b0 + b1·λ] P0.

Write α = a0 + a1·λ and β = b0 + b1·λ. Then S = [α k + β] P0. A useful collision of two states with the same point and α1 ≠ α2 gives

k ≡ (β2 − β1) · (α1 − α2)⁻¹ (mod r).

If α1 = α2 and the invariant holds, then β1 = β2 as well. That collision carries no logarithm. It is a dependent repeat.

## Coefficient update

Jumps are canonical points J = [s0 + s1·λ] P0. A raw multiple [s]P0 is stored as (s, 0) in the basis (1, λ), then multiplied by the unit that canonicalises it. The same unit index is used on the point and on the pair.

One walk step adds J, then canonicalises the sum. If u is that unit,

S′ = [u](S + J) = [u α] T + [u (β + s0 + s1·λ)] P0.

Multiplication of a pair (z0, z1) by λ, λ², −1, −λ or −λ² is the integer map coming from λ² = −1 − λ:

| index | unit | (z0, z1) |
| --- | --- | --- |
| 0 | 1 | (z0, z1) |
| 1 | λ | (−z1, z0 − z1) |
| 2 | λ² | (z1 − z0, −z0) |
| 3 | −1 | (−z0, −z1) |
| 4 | −λ | (z1, z1 − z0) |
| 5 | −λ² | (z0 − z1, z0) |

`apply_linear_unit` is that table. `central_coefficient_update` adds (s0, s1) to (b0, b1) and then applies the table to both pairs. The T-pair has a zero jump. The product therefore contains the jump. The scalar of φ_u(S) + J is a different point; that wiring is not used.

`unit_images_match_multiplication` checks the table against multiplication by the unit on random pairs. `endomorphism_matches_scalar` checks φ(P0) = [λ]P0 and that β²x agrees with a second application of β. `step_preserves_invariant` runs four production steps from a known small relation and compares the point with [α k + β] P0.

## Walk

The iteration is an r-adding walk. The jump table has 2²⁰ entries. Entry i is the canonicalisation of [s_i]P0 for a fixed nonzero scalar s_i. The index is bits [0, 20) of the canonical x. Checkpoints are the states whose canonical x has bits [20, 117) clear: a 97-bit window, nominal rate θ = 2⁻⁹⁷. The window does not overlap the index, so a checkpoint does not force a single successor.

Only checkpoints are stored. Each record holds the canonical point and the four coefficients. A later visit to the same canonical point solves for k when the α-coefficients differ. A dependent repeat takes one step with the next jump index and does not overwrite the stored node. The identity input returns 0 with no walk. The output is the integer k. One scalar multiplication checks [k]P0 = T before returning; that multiplication is charged.

Parallel schedule: m = 2¹⁶ independent walks share the checkpoint map. Their inversions are batched, which is what the track's 6M price already means. The Python group law in `reference.py` is the specification of one addition. It is not the batched instruction schedule being costed.

## Success probability

Heuristic, disclosed. On the set of N non-identity orbits the r-adding map, read through canonical representatives, is treated as a random function. For a random function the first collision time T_col of one trajectory satisfies

P(T_col > t) ∼ exp(−t² / (2N))

as N → ∞. This is the birthday count: each new image is uniform in the remaining set until it hits. It is not a theorem for this particular addition walk. Teske's r-adding walks with far fewer than 2²⁰ jumps already track random-mapping statistics; the table here is larger, not smaller. The same heuristic is the usual support for automorphism-rho costs. It does not prove mixing, and the experiment does not probe the production group order.

Take the smallest integer t with

t² / (2N) ≥ ln(5/3).

Then the asymptotic tail is at most 3/5, so the success probability is at least 2/5 = 0.4 ≥ 0.39. A rational upper bound ln(5/3) < 5108256237659908 / 10¹⁶ is used, so the inequality is integral:

t² · 10¹⁶ ≥ 2 · N · 5108256237659908,

N = (r−1)/6, t = 94490896640980656127279194705468366949.

The attack runs a fixed budget, not an expectation. It stops at the first useful checkpoint collision and otherwise stops when the budget is spent. Work is therefore at most the budget on every challenge, including failures. Runs that would have collided later are the allowed failure mass.

Distinguished points detect the collision at the next checkpoint rather than at the colliding orbit. With m walks the van Oorschot–Wiener overhead is m/θ extra group additions in total, using the nominal θ. That term is charged in full below. If the true checkpoint rate differs from 2⁻⁹⁷ by a factor inside roughly 2±30, this additive term stays below the main term; the memory bound is the quantity that moves, and it carries an explicit factor-four slack.

## Instruction schedule per ordinary step

One reference operation pays for the batched incomplete affine addition. Everything else is scalar AArch64 instruction equivalents, divided by C once.

Limb bounds, analytical, padded above a straight-line schedule. A field element is 6 limbs. A residue modulo r is 4 limbs.

- Field add or sub: 6 loads, 6 adds-with-carry, 6 loads of p, 6 subtracts-with-borrow, 6 conditional selects, 6 stores, 8 address adds. Count 44, charge 64.
- Field compare: 6 loads and a 6-limb comparison chain. Count under 30, charge 48.
- Residue add or sub modulo r, including the conditional correction: charge 40.
- Selecting one of three field elements into the output: charge 24.
- Jump-table load of a point and two residues, plus the index mask: charge 48.

Canonicalisation: 1 measured `blst_fp_mul` (383 instructions) for βx, two field additions for −x−βx, one field negation of y, three comparisons, one select.

R_canon = 383 + 2·64 + 3·48 + 64 + 24 = 743.

Coefficient update: two residue additions for the jump, then at most two residue operations for the image of each pair. The T-pair's jump is zero and is not added. Charging four residue operations covers the worst unit:

R_coeff = 4·40 = 160.

A field inversion by square-and-multiply on p−2 uses bit_length(p−2) + popcount(p−2) − 2 = 608 field multiplications, counting squares as multiplications. The 6M reference price already contains 3 of those, amortised. The surplus on a batch of m = 2¹⁶ additions is

⌈(608 − 3)·383 / 2¹⁶⌉ = 4

instructions per addition. Per-step extra instructions:

R_step = 743 + 160 + 48 + 4 = 955.

## Total work

H counts reference additions. R counts other instructions.

- Collision budget t as above.
- Distinguished-point overhead m·2⁹⁷ = 2¹⁶·2⁹⁷ = 2¹¹³ = 10384593717069655257060992658440192.
- Escape allowance 4096 reference additions, paid whether or not a dependent repeat occurs. A zero-sum cycle has heuristic probability about 1/r and is part of the failure mass if the allowance is too small to leave it.
- Exceptional-addition allowance 2²⁰. Doubles and inverse-sums are not covered by the incomplete formula. Expected accidental doubles are t/N, which is negligible. The allowance is an absolute pad, charged as reference additions.
- H_walk = t + 2¹¹³ + 4096 + 2²⁰ = 94501281234697725782536255698127859813.
- Jump-table construction: 2²⁰ scalars, each at most 254 doubles and 254 adds, so 508 reference additions, plus canonicalisation and 64 instructions of counter arithmetic per entry. H_pre = 532676608. R_pre = 2²⁰·(743+64) = 846200832.
- Verification multiplication: H_ver = 508, R_ver = 200.
- Checkpoint traffic: at most 4·⌈H_walk / 2⁹⁷⌉ records, 256 instructions each for hashing and storing. The factor 4 is slack for a mild bias of the canonical x in the checkpoint window. Record count = 2385547724. R_mem = 2385547724·256.

H = H_walk + H_pre + H_ver = 94501281234697725782536255698660536929.

R = H_walk·955 + R_pre + R_ver + R_mem = 90248723579136328122322124192323652539791.

W = (H·2000 + R) / 2000 = 279251286048531779687394635589644726397791 / 2000.

log2(W) < 126.714832 = 126714832 / 1000000. The decimal is the ceiling of log2(W) at six places. The submitted score is that rational. It is an upper bound, so a smaller true logarithm does not make the claim false.

Preprocessing work is H_pre + R_pre/2000 = 1066199416832/2000, and log2 of that quantity is less than 28.989831. It sits inside W.

## Memory

Each checkpoint record is charged 256 bytes, counting the point, four residues and hash-table overhead. Jump entries are charged 192 bytes. Code and scratch are charged 2²⁰ bytes.

Peak = 2385547724·256 + 2²⁰·192 + 2²⁰ = 610902592512 bytes.

log2 of that peak is less than 39.152152. The jump table is computed by the attack and is preprocessing, not uncharged advice. Advice size is 1 byte in the log2(max(1, ·)) convention, so the advice logarithm is 0.

## What is heuristic

The collision tail, the treatment of m simultaneous walks as one birthday process plus m/θ, and the closeness of the checkpoint rate to 2⁻⁹⁷ are heuristics. The unit maps, the identity β²+β+1 = 0, the orbit size, the inversion surplus and the limb schedules are exact relative to the stated machine model. Small-group recoveries show that this coefficient wiring solves a cyclic group. They do not measure N, θ or the instruction counts.

## Infinity and validity

If T is the identity, the scalar is 0. Otherwise the relation above is the scalar. Organiser decoding and the subgroup check on the submitted answer are organiser work. The attack's own check [k]P0 = T is inside H_ver.

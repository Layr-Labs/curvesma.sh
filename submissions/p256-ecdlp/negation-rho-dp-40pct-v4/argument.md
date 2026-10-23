# p256-ecdlp v4: parallel Pollard rho, distinguished points, negation map

Analytic cost bound for recovering x from Q=[x]P on NIST P-256
(secp256r1, exact SEC 2 2.4.2 constants, a=-3, cofactor 1). All work is
charged in scalar AArch64 instruction equivalents divided by the protected
C=1000 (the borrowed secp256k1 normalization unit; not a measured P-256
kernel price). This is a backport re-derivation of the historical
unpublished design `submissions/p256-ecdlp/rho-negmap-dp` (research refs
9ed0ff21… / ad846d88…, author unknown, never published); every identity
below was independently recomputed with host python3 exact arithmetic and
is recorded in certificate.json. First claim: time_log2 = 1283/10.

## 1. Setting

n = 0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551
(the subgroup order, verified against target-profile.json), log2 n =
255.9999999996641. The negation map identifies R with -R; the identity is
a singleton class, so the class space has size N' = (n+1)/2, log2 N' =
254.9999999996641. P-256 has j != 0 and no p≡1 mod 3 endomorphism of the
form (x,y)->(βx,y) on this curve family, so only the size-2 negation
classes apply: the classical √2 = 2^0.5 speedup, not the √6 of
secp256k1/BLS12-381-G1. No Cheon or Pohlig-Hellman leverage exists (r is
prime, r-1 has no factor structure usable within the input model).

## 2. Algorithm

Walk table: R = 2^24 entries T_j = c_j P + d_j Q, c_j,d_j uniform nonzero
64-bit, stored affine with coefficients; 80 B/entry. For a canonical point
X=(u,v) (Montgomery limbs), j(X) = bits 40..63 of the top limb of u
(negation-invariant), and f(X) = canon(X + T_{j(X)}), canon picking the
member of {Y,-Y} with smaller 256-bit Montgomery y. f is a deterministic
function on classes (r-adding walk, Teske 2001, r = 2^24 >> 16).

Walkers: m = 2^14, lockstep, starting points a0 P + b0 Q (uniform 256-bit
coefficients). Each keeps (u,v), sign s, signed un-reduced (A,B) with
point = s(AP+BQ); a step adds T_j and flips s on canonicalisation.
|A|,|B| < 2^(127.57+64+1) < 2^193 < 2^255, so 4-limb signed tracking
never needs reduction during the walk.

Distinguished points: theta = 2^-96 (low 96 bits of u zero). A DP record
(u, A mod n, B mod n) = 96 B goes to a store keyed by u. Two records with
equal u and different (A,B) give x = (±A2 - A1)/(B1 ∓ B2) mod n, tried
for both signs; B1 = ±B2 has probability ~2^-250. Every candidate is
verified by [x]P = Q.

Fruitless cycles (Bos-Kaihara-Kleinjung-Lenstra-Montgomery 2010;
Bernstein-Lange-Schwabe 2011): with the negation map, 2-cycles of f occur
with probability 1/(2R) = 2^-25 per step; longer 2l-cycles at rate at
most 2(2l/R)^l. Each walker runs a fixed 64-step checkpoint schedule;
detection implies being on a cycle of f, and any cycle of length <= 64 is
detected within 128 steps. On detection the walker laps the cycle once,
finds the unique smallest-u member m*, and escapes by adding the fixed
point T'_{e(m*)} from a second 256-entry table: the exit point E(C) is a
function of the cycle only, so merged walkers re-merge. Cycle-event cost
<= 128 + 66 + 61 = 255 step equivalents; amortised loss 256·λ with
λ <= 2^-25 + Σ_{l≥2} 2(2l/R)^l = 2^-24.99 -> loss 2^-16.99 <= 2^-12
(charged in the K denominator below).

Batched step (the reference operation): per round, one Montgomery batch
inversion of the 2^14 denominators (3M per walker amortised: prefix
product, one Fermat inversion, unwind) plus the 3M affine addition
(u',v') = (u,v) + T_j with S = M: total 6M per step, the same 6M
convention the reference unit itself is calibrated on. Zero-denominator
cases (X = ±T_j, prob ~2^-231) go to a charged slow path.

## 3. Success probability

Random-mapping model on N' classes (heuristic, disclosed): collision
probability among m trails of total length k is at least
1 - exp(-k(k-1)/(2N')). With k = 1.04377·sqrt(N') = sqrt(2N' ln(1/0.58))
this is 1 - 0.58 = 0.42 (verified: exp(-k_ratio²/2) = 0.58). Failure
terms: (i) the collision lies inside one walker's own trail, expected
count k²/(2N')/m < 2^-14 (verified); (ii) a merged pair has no DP within
D = 24/θ = 24·2^96 extra steps, at most 2e^-24 < 10^-10 (verified).
Success >= 0.42 - 2^-14 - 2e^-24 = 0.419939. Claimed floor: 2/5 = 0.40
(margin 0.0299 over the 0.39 rules floor), absorbing walk non-randomness
(Teske reports few-% deviations at r=20; r=2^24 here). No restarts: one
capped run at K = (k + mD)/(1-2^-12) < 2^127.5623 steps (log2 K =
127.56234, verified).

## 4. Per-step instruction ledger

P-256 Montgomery multiplication M: 4x64-bit schoolbook 16 MUL + 16 UMULH
+ ~44 carry ops, 4 fused reduction rounds of p = 2^256-2^224+2^192+2^96-1
at ~16 ops each, conditional subtract ~10, loads/stores/call ~16:
~185 actual; charged M = 220 (18% margin). S = M. Modular add/sub
A = 24. Fermat inversion 255 squarings + 12 muls <= 300M = 66000 instr,
amortised over 2^14 walkers: 66000/2^14 = 4.0283 instr.

| component | instructions |
| --- | ---: |
| 6M (3 addition + 3 batch inversion), 6×220 | 1320 |
| 6 modular add/sub, 6×24 | 144 |
| zero-denominator test | 6 |
| table index (UBFX), address, 4-word loads u_T,v_T,c_j,d_j | 24 |
| negation canonicalisation: 4-limb compare, conditional negate, flip s | 32 |
| signed 4-limb coefficient updates A,B with load/store | 64 |
| DP test (cbnz; second test on rare hit) | 4 |
| cycle counter, checkpoint compare/refresh, mode flag | 16 |
| walker state movement, prefix-product traffic, loop | 50 |
| amortised inversion 66000/2^14 | 4.0283203125 |
| **total** | **1664.0283203125** |

= 1.66402832 units/step < 1.66403 (exact Fraction check PASS).
log2 u = 0.734680.

## 5. Other charges (units = instr/1000, each verified)

* Walk: K·u = 2^127.56234 × 1.66402832 = 2^128.29701.
* Table build: per entry two 64-bit scalar muls with a=-3 (64D×8M + ≤64A×11M = 2432M) + 11M sum + 10M batch-affine < 2500M = 550 units -> 2^24·550 = 2^33.103.
* Coefficient generation: 2^25 words × 20 instr = 2^21.1.
* Walker starts: 2^14 × 2 × 256-bit scalar muls ≤ 2^14×10000M×220 = 2^25.103.
* Escape table: 256×550 = 2^17.1.
* DP handling: K·θ = 2^31.562 records × 4000 instr = 2^33.562.
* Verification/output: ≤10^4 = 2^13.3. Slow paths: 2^128·2^-231·2000 ≈ 2^-87.9.
* Cycle events: inside the 2^-12 factor already in K.

Extras sum 2^34.353, a relative 2^-93.94 of the walk term (recomputed;
the historical proof said 2^-94 — same magnitude, ours is the honest sum
including the larger 4000-instr DP record cost). Total
W = 2^128.29701·(1+2^-12) + 2^34.353 = 2^128.29740 < 2^128.3.
Claimed: **time_log2 = 128.30 = 1283/10** (verified PASS).

Preprocessing (table + coefficients + starts + escape table)
= 2^33.108 < 2^33.2, included in W. Advice: 0 bytes.

Memory (reported, unscored per rules): DP store 2^34 records × 96 B
(Poisson headroom over 2^31.56 expected) = 2^40.585 B + walk table
2^30.32 B + walkers 2^22.64 B < 2^41 B (~2 TB).

## 6. Heuristics, support, limits (disclosed)

1. Random-mapping behaviour of the r-adding walk on negation classes
   (Teske 2001). Support: r = 2^24 >> 16 tested range; 7% probability
   margin (0.42 -> 0.40). Limit: analytic only, no P-256-scale run; the
   0.02 slack over 0.39 is the acceptance hedge.
2. Cycle rate 2(2l/R)^l assumes table-index independence; a counting
   bound, not a theorem. Detected-cycle escape is exact given detection.
3. DP test (low 96 bits of u) independent of table index (bits 40..63 of
   the top limb) for well-mixed u.
4. M=220 / A=24 are conservative analytic instruction counts for scalar
   AArch64 Montgomery arithmetic (+18% over the ~185 estimate), not
   measured kernels; the borrowed C=1000 is a normalization constant.
5. Backport provenance: the mechanism, constants (R=2^24, θ=2^-96,
   m=2^14, M=220) and success budget are the historical design's; the
   arithmetic ledger, extras, memory and success identities here were
   re-derived independently and match its 128.30 claim to <0.0004 bits.

References: van Oorschot-Wiener, J. Cryptology 1999; Teske, Math. Comp.
2001; Bos-Kaihara-Kleinjung-Lenstra-Montgomery 2010; Bernstein-Lange-
Schwabe, PKC 2011.

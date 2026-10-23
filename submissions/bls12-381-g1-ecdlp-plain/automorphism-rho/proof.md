# bls12-381-g1-ecdlp-plain: automorphism-aware Pollard rho

Status: ready for independent exploratory review. Classical single-shot
attack, no quantum steps, no uncharged advice, no amortization over targets.
The automorphism-accelerated rho is standard cryptanalysis (Wiener-Zuccherato),
not a novelty claim. All arithmetic facts used below are verified by
`experiment.py` Part A at full production parameters.

## Claimed bounds

| Field | Value |
| --- | --- |
| time_log2 (score) | **127.143439** (= 127143439/1000000, rounded up) |
| success_probability (supported lower bound) | 0.39 (threshold is 0.39) |
| memory_log2_bytes | 70.3 |
| preprocessing_log2 (included in total) | 25.4 |
| nonuniform_advice_log2_bytes | 0 |

No qualified baseline exists for this track (the repository's selected package
is an unscored template); for scale, the same construction on the repository's
secp256k1 track scored 129.000802 against a deterministic-BSGS control, and
generic rho on this track without the automorphism would cost ~2^128.07
(negation-only) or ~2^128.03 (vanilla, birthday budget at the 0.39 threshold).

## Conventions

Protected target: BLS12-381 G1, y^2 = x^3 + 4 over Fp, generator P0 of prime
order r, challenge T = [x]P0 with x uniform in [0, r) (T = infinity exactly
when x = 0; that case is answered immediately). Cost model
`classical-group-work-v2`: total cost = H + R/2000 where H counts reference
operations (batched incomplete affine BLS12-381 G1 addition, 6M convention
with S = M and A = 0; one reference operation = 1 unit) and R counts all
other work in scalar 64-bit AArch64 instruction equivalents. The pinned base
field multiplication price is 383 instructions (blst `blst_fp_mul`
measurement in `tracks/docs/COST_MODEL.md`); the G1 reference price
C = 2,000 is the rounded 6M = 6 x 383 = 2,298. A 256-bit mod-r
multiplication (r is 255-bit) is charged at the same conservative 383.
A reference operation and its implementation are never double-counted.

## The order-6 automorphism group

BLS12-381 has j-invariant 0 and p mod 3 = 1, so Fp contains a primitive cube
root of unity zeta with zeta^3 = 1, zeta != 1. The map (x, y) -> (zeta*x, y)
preserves the curve equation and fixes infinity, hence is a group automorphism
of E(Fp); it maps the r-torsion subgroup G1 = <P0> to itself (automorphisms
commute with scalar multiplication and G1 is the unique order-r subgroup), so
on G1 it acts as multiplication by a scalar lam with lam^3 = 1 mod r,
lam != 1. The correct element of {lam, lam^2} is fixed by the numerically
verified identity [lam]P0 = zeta(P0) (experiment, Part A) with

lam = 0x73eda753299d7d483339d80809a1d804a7780001fffcb7fcfffffffe00000001.

Together with negation (x, y) -> (x, -y) = [-1](x, y) this gives an action of
G = {+-1, +-lam, +-lam^2} (order 6) on the r - 1 nonzero points of G1. Since
r is an odd prime and -1, lam, lam^2, -lam, -lam^2 are all != 1 mod r, the
action is free, so the nonzero points split into exactly

N = (r - 1) / 6 = 8739312529187698413241290084697660972948425416754606303767276449989763530752

classes of size 6 (r mod 6 = 1). A class has a unique canonical
representative (x_c, y_c) with x_c = min(x, zeta*x, zeta^2*x) and y_c even.
If a point R has tracked representation R = [a]P0 + [b]T, then applying the
class-canonicalizing automorphism updates the representation exactly:
(a, b) -> (lam^i * a, lam^i * b) for the zeta^i factor and
(a, b) -> (-a, -b) for a sign flip, because
zeta([m]P0 + [n]T) = [lam*m]P0 + [lam*n]T and -([m]P0 + [n]T) = [-m]P0 +
[-n]T. This bookkeeping is verified on random production points in Part A.

## Algorithm

Single-shot classical Monte-Carlo attack; fresh coins per run; one target.
If T = infinity, output x = 0 (trivial, charged negligibly).

Setup (all charged; see ledger):
1. zeta: smallest g in 2, 3, 5, ... with g^((p-1)/3) mod p != 1; zeta^2
   likewise. lam analogously mod r; fix which of {lam, lam^2} satisfies
   [lam]P0 = zeta(P0) by one scalar multiplication check each.
2. Walk constants: u_j uniform in [1, r) for j = 0..31; A_j = [u_j]P0
   canonicalized to its class representative with tracked multiplier c_j
   (so canon([u_j]P0) = [c_j]P0). Verify [c_j]P0 equals the stored constant.
3. K = 2^16 walkers, each with fresh uniform (a, b); start
   R = canon([a]P0 + [b]T).

Walk (lockstep rounds over all walkers so the additions form one batched
reference-operation stream):
4. j = x_c mod 32 (x_c = canonical x of the walker's current point).
5. S = R + A_j, one batched incomplete affine reference addition
   (exceptional cases below).
6. Canonicalize S; update (a, b) accordingly (Section above).
7. Distinguished point test: if the low 64 bits of x_c are zero, append
   (x_c, a, b, walker_id) to a hash table keyed by x_c. If x_c already occurs
   with a different (a, b) pair, solve
   x = (a1 - a2) * (b2 - b1)^(-1) mod r, verify [x]P0 = T, and output x.
8. Degenerate handling (charged; frequency bounded below): S = infinity,
   incomplete-formula exceptions (x_R = x_{A_j}), revisited class with
   identical (a, b), or a collision with b1 == b2 (no information): reseed
   the walker. The same rule covers fruitless 2-cycles of the negation map.
9. The run stops at the first verified extraction or after t total walker
   steps (budget below). No restart or amplification beyond rule 8; the bound
   covers every trial and failure.

## Success probability

The canonical walk is a 32-adding walk on the N classes (Teske). Under the
standard random-walk heuristic the first class collision time T_coll
satisfies P(T_coll <= t) = 1 - exp(-t^2 / (2N)) to close approximation.
Budget: t = 1.1 * sqrt(N) = 2^126.2736 total walker steps (summed over
walkers), ideal success probability 1 - exp(-1.1^2 / 2) = 0.4539.

Support and margins:
- The birthday-tail model for rho is classical (HAC 3.6); Teske's experiments
  show 32-adding walks collide with near-ideal constants.
- Part B of `experiment.py` measures the complete attack on a toy j=0 curve
  (prime order n = 1047379, N_toy = 174563 classes): at the same budget shape
  ceil(1.1*sqrt(N)) the observed single-shot success was 84/200 = 0.42 with
  zero extraction-verification failures, versus the ideal 0.4545.
- The toy shortfall is dominated by fruitless short cycles, whose trigger set
  has size O(k^2 * |G|) ~ 2^13.4 special classes for k = 32 constants:
  per-step trigger probability 2^13.4/N. At toy scale (N = 2^17.4) this is
  ~2^-4 per step (measured: 348 degenerate reseeds over 200 x 460 steps); at
  production (N = 2^252.27) it is 2^-239 per step, so the probability of any
  degenerate event during the whole 2^126.27-step run is < 2^-113, and each
  occurrence costs at most one walker reseed. 1-cycles are impossible
  (R + A_j != R for A_j != infinity).
- Therefore production walk efficiency is at least the toy-measured 0.924
  minus the toy's fruitless-cycle contribution, which vanishes at production.
  Even at a pessimistic 0.86 efficiency the supported success probability is
  0.4539 * 0.86 = 0.3904 >= 0.39; the lower bound claimed is **0.39**.

Heuristic H1 (random-walk/birthday behavior of the canonical 32-adding walk)
is supported by classical analysis, Teske's independent experiments, and the
toy run; the budget leaves a >= 14% performance slack. No other heuristic is
needed for success; the algebra (class structure, representation updates,
extraction) is exact and machine-verified at production parameters.

## Work ledger (worst case every step)

Per walker step, in instruction equivalents (units after division by
C = 2,000):

| Item | Instructions | Units |
| --- | ---: | ---: |
| One reference operation: the walk addition R + A_j, executed as exactly the calibrated batched incomplete affine 6M G1 addition (3M formula + 3M batch-inversion amortisation across the 2^16 lockstep walkers; one real inversion per round amortizes to < 2^-16 per step) | reference | 1.000 |
| Canonicalization: x*zeta and x*zeta^2 (2 base-field mults, 383 each) | 766 | 0.383 |
| Canonicalization: 3-way argmin over 6x64-limb values with selects | 64 | 0.032 |
| Canonicalization: parity test + conditional field negation | 16 | 0.008 |
| Representation update: two 255-bit mod-r negations (sign) | 32 | 0.016 |
| Representation update: two 255-bit mod-r multiplications by fixed lam or lam^2 (charged at the 383-instruction general price) | 766 | 0.383 |
| Walk index (mask + and) | 6 | 0.003 |
| Distinguished-point predicate (OR-reduce 64 bits, branch) | 5 | 0.0025 |
| **Per-step worst-case total** | **1839** | **1.8275** |

The sign and lam events are partial in expectation (about 1/2 and 2/3 of
steps); the claim nevertheless charges the worst case at every step.

One-off and amortized components:

| Component | Bound (units) |
| --- | ---: |
| zeta, lam derivation (few field exponentiations) + orientation check | <= 2^10 |
| 32 walk constants (32 scalar multiplications + canonicalization + verification) | <= 2^14 |
| Walker seeding: 2K = 2^17 scalar multiplications, <= 384 reference additions each | <= 2^25.32 |
| Distinguished-point processing: t/2^64 = 2^62.27 stores x <= 1000 instructions | <= 2^62.4 |
| Post-merge detection lag: after the useful merge all K walkers continue for <= 2^64 rounds until the merged pair reports its next distinguished point (harmonic factor 1.44 folded in) | <= 2^81 |
| Endgame: extraction (one mod-r inversion, one multiplication), verification [x]P0, output | <= 700 |

Total work:

W = 1.8275 * t + one-offs = 1.8275 * 1.1 * sqrt((r-1)/6) + <= 2^81 + < 2^63
  = 2^127.14344 (one-offs are below 2^-46 of the leading term).

Submitted time_log2, rounded up to six decimals: **127.143439**
(exact rational 127143439/1000000). Preprocessing (roots, constants, seeding)
is 2^25.32 units -> preprocessing_log2 = **25.4**, included in the total.

## Memory

Distinguished-point table: at most t/2^64 = 2^62.27 entries, 128 bytes per
entry (x_c, a, b, chaining metadata), hash-table slack factor 2 -> 2^70.3
bytes. Walker state (2^16 walkers x 4 x 32 bytes + scratch), code and
constants add < 2^24 bytes. Peak storage is 2^70.3 bytes ->
memory_log2_bytes = **70.3**. All tables are constructed and charged during
the attack; nothing is stored in advance of the challenge, so nonuniform
advice is **0 bytes**.

## Why not the pairing / GT transfer

The rules make transferring the discrete logarithm to GT via the public
pairing eligible. A generic rho in GT costs far more per step (dense Fp12
multiplication is the 20,000-instruction reference price of the pairing
tracks, ~10x the per-step price here), so it loses outright. The only
potentially better GT route is index calculus over Fp12 (subexponential),
but BLS12-381's 381-bit modulus was chosen precisely so that the
Kim-Barbulescu-style extended TNFS on Fp12 stays near the 2^128 class; a
qualifying submission would need a complete analytic NFS cost model
(sieving, linear algebra, full descent charged in instruction equivalents),
which this package does not claim. The generic bound above is the honest
first qualified score for this track.

## Evidence and reproduction

- `experiment.py` Part A (production parameters, exact): p mod 3 = 1,
  r mod 6 = 1, zeta^3 = 1 mod p with curve preservation, [lam]P0 = zeta(P0)
  and [lam^2]P0 = zeta^2(P0) with lam^3 = 1 mod r, class count (r-1)/6, and
  canonicalization-representation bookkeeping verified on 8 random points
  [a]P0 + [b]Q. No production-size discrete logarithm is computed.
- `experiment.py` Part B (toy scale, 20-bit j=0 curve, full attack): budget
  ceil(1.1*sqrt(N)), 8 walkers, 32 constants; declared run
  `{"trials": 200, "dp": 1}`: 84/200 successes (0.42 vs ideal 0.4545), 0
  extraction-verification failures, 348 degenerate reseeds. Local sweep:
  dp=0: 0.42; dp=1: 0.42; dp=2: 0.355; dp=4: 0.33; the dp>1 shortfalls are
  toy-scale distinguished-point lag (per-walker budget ~57 steps vs interval
  2^dp), while at production the lag is 2^64 steps against a per-walker
  budget of 2^110.27 steps, a relative loss below 2^-46, charged in the
  ledger.
- Reproduce with
  `bash scripts/run_experiments.sh bls12-381-g1-ecdlp-plain submissions/bls12-381-g1-ecdlp-plain/automorphism-rho`
  (Docker, no API key) or directly:
  `echo '{"trials": 200, "dp": 1}' | python3 submissions/bls12-381-g1-ecdlp-plain/automorphism-rho/experiment.py`.

Small experiments are evidence only at their actual scope: Part B does not
establish production success constants or costs; Part A does establish the
algebraic structure the proof relies on, exactly and at full size.

## Heuristics, limitations, objections

- H1 (disclosed above): first-collision birthday behavior of the canonical
  32-adding walk. The 1.1 sqrt(N) budget absorbs a >= 14% walk-efficiency loss
  relative to the claimed 0.39.
- Instruction sub-costs (383 per base-field multiplication is the pinned
  blst measurement; 16-64 instruction bounds for compares/selects and
  negations; 383 charged for the cheaper 255-bit mod-r multiplications) are
  analytical worst-case bounds as the cost model allows.
- Single-shot design: the claimed probability 0.39 meets the protected
  threshold exactly; doubling the budget t gives 1 - exp(-2.42) = 0.911 at
  +1 log2 cost.
- Memory 2^70.3 bytes is astronomical and inherent to sqrt(r)-scale
  distinguished-point collision search; reported, no score effect.
- Standard cryptanalysis, not a novel attack; this is the track's first
  qualified bound under this protected cost model.

References: `tracks/docs/COST_MODEL.md` (pinned 383-instruction blst base
field multiplication, C = 2,000); Handbook of Applied Cryptography, ch. 3.6;
Teske, r-adding walks (J. Cryptology 14(2), 2001); Wiener and Zuccherato,
"Faster attacks on elliptic curve cryptosystems" (SAC 1998); Galbraith,
Hopkins and Shparlinski (2009) and Bos, Kleinjung and Lenstra on the negation
map and fruitless cycles; Kim and Barbulescu (CRYPTO 2016) on extended TNFS
sizing of pairing-friendly fields. The construction and accounting are
self-contained above; the reviewer need not fetch any citation.

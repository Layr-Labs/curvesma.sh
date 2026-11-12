# secp256k1-ecdlp: automorphism-aware Pollard rho

Status: ready for independent exploratory review. Classical single-shot attack,
no quantum steps, no uncharged advice, no amortization over targets. This is a
best-known-generic-attack-style claim under the protected cost model, not a
practical break or novelty claim; the automorphism-accelerated rho is standard
cryptanalysis (Wiener-Zuccherato). All arithmetic facts used below are
verified by `experiment.py` (Part A, production parameters).

## Claimed bounds

| Field | Value |
| --- | --- |
| time_log2 (score) | **127.723944** (= 15965493/125000, rounded up) |
| success_probability (supported lower bound) | 0.39 (threshold is 0.39) |
| memory_log2_bytes | 70.9 |
| preprocessing_log2 (included in total) | 25.4 |
| nonuniform_advice_log2_bytes | 0 |

Baseline control in this repository: 129.000802 with 2^136.000001 bytes
(deterministic BSGS). This claim is lower by 1.276858 log2 (~2.4x less work)
with ~2^65 times less memory.

## Conventions

Protected target: secp256k1, y^2 = x^3 + 7 over Fp, generator P of prime order
r, cofactor 1, challenge Q = [x]P with x uniform in [0, r). Cost model
`classical-group-work-v2`: total cost = H + R/1000 where H counts reference
operations (batched incomplete affine secp256k1 addition, 6M convention with
S = M and A = 0; one reference operation = 1 unit) and R counts all other
work in scalar 64-bit AArch64 instruction equivalents. The pinned field
multiplication price is 179 instructions (Bitcoin Core 5x52 measurement in
`tracks/docs/COST_MODEL.md`); a 256-bit mod-r multiplication is charged at
the same conservative 179. A reference operation and its implementation are
never double-counted.

## The order-6 automorphism group

secp256k1 has j-invariant 0 and p mod 3 = 1, so Fp contains a primitive cube
root of unity zeta with zeta^3 = 1, zeta != 1. The map (x, y) -> (zeta*x, y)
preserves the curve equation and fixes infinity, hence is a group automorphism.
On the cyclic group <P> it acts as multiplication by a scalar lam with
lam^3 = 1 mod r, lam != 1; the correct element of {lam, lam^2} is fixed by the
numerically verified identity [lam]P = zeta(P) (experiment, Part A). Together
with negation (x, y) -> (x, -y) = [-1](x, y) this gives an action of
G = {+-1, +-lam, +-lam^2} (order 6) on <P>. Since r is an odd prime and
-1, lam, lam^2, -lam, -lam^2 are all != 1 mod r, the action of G on the r - 1
nonzero points is free, so the nonzero points split into exactly

N = (r - 1) / 6 = 19298681539552699237261830834781317975472927379845817397100860523586360249056

classes of size 6. A class has a unique canonical representative:
(x_c, y_c) with x_c = min(x, zeta*x, zeta^2*x) and y_c even. If a point R has
tracked representation R = [a]P + [b]Q, then applying the class-canonicalizing
automorphism updates the representation exactly: (a, b) -> (lam^i * a, lam^i * b)
for the zeta^i factor and (a, b) -> (-a, -b) for a sign flip, because
zeta([m]P + [n]Q) = [lam*m]P + [lam*n]Q and -([m]P + [n]Q) = [-m]P + [-n]Q.
This bookkeeping is verified on random production points in Part A.

## Algorithm

Single-shot classical Monte-Carlo attack; fresh coins per run; one target.

Setup (all charged; see ledger):
1. zeta: smallest g in 2, 3, 5, ... with g^((p-1)/3) mod p != 1; zeta^2 likewise.
   lam analogously mod r; fix which of {lam, lam^2} satisfies [lam]P = zeta(P)
   by one scalar multiplication check each.
2. Walk constants: u_j uniform in [1, r) for j = 0..31; A_j = [u_j]P
   canonicalized to its class representative with tracked multiplier c_j
   (so canon([u_j]P) = [c_j]P). Verify [c_j]P equals the stored constant.
3. K = 2^16 walkers, each with fresh uniform (a, b); start R = canon([a]P + [b]Q).

Walk (lockstep rounds over all walkers so that the additions form one batched
reference-operation stream):
4. j = x_c mod 32 (x_c = canonical x of the walker's current point).
5. S = R + A_j, one batched incomplete affine reference addition (exceptional
   cases below).
6. Canonicalize S; update (a, b) accordingly (Section above).
7. Distinguished point test: if the low 64 bits of x_c are zero, append
   (x_c, a, b, walker_id) to a hash table keyed by x_c. If x_c already occurs
   with a different (a, b) pair, solve
   x = (a1 - a2) * (b2 - b1)^(-1) mod r, verify [x]P = Q, and output x.
8. Degenerate handling (charged; frequency bounded in the next section):
   - S = infinity or x_R = x_{A_j} (incomplete-formula exception): reseed the
     walker (fresh (a, b), new start point).
   - Revisited class with identical (a, b), or a collision with b1 == b2
     (no information): reseed the walker.
   The same rule covers fruitless 2-cycles of the negation map.
9. The run stops at the first verified extraction or after t total walker
   steps (budget below). No restart or amplification beyond rule 8; the bound
   covers every trial and failure.

## Success probability

The canonical walk is a 32-adding walk on the N classes (Teske). Under the
standard random-walk heuristic the first class collision time T satisfies
P(T <= t) = 1 - exp(-t^2 / (2N)) to close approximation. Budget:

t = 1.1 * sqrt(N) = 2^126.8450 total walker steps (summed over walkers),

giving an ideal success probability 1 - exp(-1.1^2 / 2) = 0.4539.

Support and margins:
- The birthday-tail model for rho is classical (HAC 3.6); Teske's experiments
  show 32-adding walks collide with near-ideal constants.
- Part B of `experiment.py` measures the complete attack on a toy j=0 curve
  (order n = 1047379, N_toy = (n-1)/6 = 174563 classes): at the same budget
  shape ceil(1.1*sqrt(N)) the observed single-shot success was 84/200 = 0.42
  with zero extraction verification failures over all trials, versus the
  ideal 0.4545 (walk efficiency 0.924 at 20-bit scale).
- The toy shortfall is dominated by fruitless short cycles, whose trigger set
  has size O(k^2 * |G|) ~ 2^13.4 special classes for k = 32 constants: the
  per-step trigger probability is 2^13.4/N. At toy scale (N = 2^17.4) this is
  ~2^-4 per step (measured: 348 degenerate reseeds over 200 x 460 steps);
  at production (N = 2^253.4) it is 2^-240 per step, so the probability that
  any degenerate event occurs during the whole 2^126.85-step run is < 2^-113,
  and each occurrence costs at most one walker reseed. The 1-cycles are
  impossible (R + A_j != R for A_j != infinity).
- Therefore the production walk efficiency is at least the toy-measured 0.924
  minus the toy's fruitless-cycle contribution, which vanishes at production.
  Even at a pessimistic 0.86 efficiency the supported success probability is
  0.4539 * 0.86 = 0.3904 >= 0.39; the lower bound claimed is **0.39**.

Heuristic H1 (random-walk/birthday behavior of the canonical 32-adding walk on
classes) is thus supported by classical analysis, by Teske's independent
experiments, and by the toy run, and the budget leaves a >= 14% performance
slack against the claimed bound. No other heuristic is needed for success;
the algebra (class structure, representation updates, extraction) is exact
and machine-verified at production parameters.

## Work ledger (worst case every step)

Per walker step, in instruction equivalents (units after division by C = 1000):

| Item | Instructions | Units |
| --- | ---: | ---: |
| One reference operation: the walk addition R + A_j, executed as exactly the calibrated batched incomplete affine addition (3M formula + 3M batch-inversion amortisation across the 2^16 lockstep walkers; one real inversion per round amortizes to < 2^-16 per step) | reference | 1.000 |
| Canonicalization: x*zeta and x*zeta^2 (2 field mults, 179 each) | 358 | 0.358 |
| Canonicalization: 3-way argmin over 5x52-limb values with selects | 64 | 0.064 |
| Canonicalization: parity test + conditional field negation | 16 | 0.016 |
| Representation update: two 256-bit mod-r negations (sign) | 32 | 0.032 |
| Representation update: two 256-bit mod-r multiplications by fixed lam or lam^2 (charged at the general 5x52 price) | 358 | 0.358 |
| Walk index (mask + and) | 6 | 0.006 |
| Distinguished-point predicate (OR-reduce 64 bits, branch) | 5 | 0.005 |
| **Per-step worst-case total** | **1839** | **1.839** |

The sign and lam events are partial in expectation (about 1/2 and 2/3 of
steps); the claim nevertheless charges the worst case at every step.

One-off and amortized components:

| Component | Bound (units) |
| --- | ---: |
| zeta, lam derivation (few field exponentiations) + orientation check | <= 2^10 |
| 32 walk constants (32 scalar multiplications + canonicalization + verification) | <= 2^14 |
| Walker seeding: 2K = 2^17 scalar multiplications, <= 385 reference additions each (double-and-add with batched inversion across walkers) | <= 2^25.32 |
| Distinguished-point processing: t/2^64 = 2^62.85 stores x <= 1000 instructions | <= 2^62.9 |
| Post-merge detection lag: after the useful merge all K walkers continue for <= 2^64 rounds until the merged pair reports its next distinguished point (harmonic factor 1.44 folded in) | <= 2^81 |
| Endgame: extraction (one mod-r inversion, one multiplication), verification [x]P (scalar multiplication), output | <= 700 |

Total work:

W = 1.839 * t + (one-offs) = 1.839 * 1.1 * sqrt((r-1)/6) + <= 2^81 + < 2^63
  = 2^127.72394 (the one-offs are below 2^-46 of the leading term).

Submitted time_log2, rounded up to six decimals: **127.723944**
(exact rational 15965493/125000). Preprocessing (roots, constants, seeding)
is 2^25.32 units -> preprocessing_log2 = **25.4**, included in the total.

## Memory

Distinguished-point table: at most t/2^64 = 2^62.85 entries, 128 bytes per
entry (x_c, a, b, chaining metadata), hash-table slack factor 2 -> 2^70.85
bytes. Walker state (2^16 walkers x 4 x 32 bytes + scratch), code and
constants add < 2^24 bytes. Peak storage is 2^70.9 bytes -> memory_log2_bytes
= **70.9**. All tables are constructed and charged during the attack; nothing
is stored in advance of the challenge, so nonuniform advice is **0 bytes**
(runtime state is not advice; its construction is fully charged above).

## Evidence and reproduction

- `experiment.py` Part A (production parameters, exact): p mod 3 = 1,
  r mod 6 = 1, zeta^3 = 1 mod p, zeta preserves the curve, [lam]P = zeta(P)
  and [lam^2]P = zeta^2(P) with lam^3 = 1 mod r, class count N = (r-1)/6, and
  canonicalization-representation bookkeeping verified on 12 random points of
  the form [a]P + [b]Q. No production-size discrete logarithm is computed.
- `experiment.py` Part B (toy scale, 20-bit j=0 curve, full attack): budget
  ceil(1.1*sqrt(N)), 8 walkers, 32 constants; declared run
  `{"trials": 200, "dp": 1}`: 84/200 successes (0.42 vs ideal 0.4545), 0
  extraction-verification failures, 348 degenerate reseeds. Local sweep:
  dp=0: 0.42; dp=1: 0.42; dp=2: 0.355; dp=4: 0.33. The dp>1 shortfalls are
  the distinguished-point lag at toy scale (per-walker budget ~57 steps vs
  interval 2^dp); at production the lag is 2^64 steps against a per-walker
  budget of 2^111 steps, a relative loss below 2^-46, and the post-merge lag
  is charged in full in the ledger.
- Reproduce with
  `bash scripts/run_experiments.sh secp256k1-ecdlp submissions/secp256k1-ecdlp/automorphism-rho`
  (Docker, no API key) or directly:
  `echo '{"trials": 200, "dp": 1}' | python3 submissions/secp256k1-ecdlp/automorphism-rho/experiment.py`.

Small experiments are evidence only at their actual scope: Part B does not
establish production success constants or costs; Part A does establish the
algebraic structure the proof relies on, exactly and at full size.

## Heuristics, limitations, objections

- H1 (disclosed above): first-collision birthday behavior of the canonical
  32-adding walk. Supported; the budget absorbs a >= 14% walk-efficiency loss
  relative to the claimed 0.39.
- Instruction-equivalent sub-costs (179 per multiplication is the pinned
  measurement; 16-64 instruction bounds for compares/selects/negations) are
  analytical bounds as the cost model allows; they are charged at worst case
  every step.
- The attack is single-shot: the claimed probability mass is 0.39 >= the
  0.39 threshold; there is no amplification, and failures inside the budget
  are fully charged. If the reviewer demands higher success, note
  P = 1 - exp(-t^2/2N) is monotone: doubling the budget t gives 1 - exp(-2.42) = 0.911 at +1 log2 cost.
- Memory 2^70.9 bytes is astronomical and mandatory for the distinguished-
  point design; it is reported, charged where applicable, and far below the
  baseline's 2^136 bytes. A smaller table only changes the lag term.
- Not novel cryptanalysis: automorphism-aware rho for j = 0 curves and the
  ~2^127 effective size of secp256k1 under it are folklore; this submission
  supplies the track's first qualified bound under this protected cost model,
  improving on the repository's BSGS accounting control by 1.276858 log2.
- No production-size run is claimed anywhere; the analytical bound is the
  claim, per the cost model's allowance for analytical instruction-equivalent
  bounds.

References: `tracks/docs/COST_MODEL.md` (pinned 179-instruction multiplier,
C = 1000); Handbook of Applied Cryptography, ch. 3.6 (generic DLP algorithms,
rho and birthday analysis); Teske, r-adding walks (J. Cryptology 14(2), 2001);
Wiener and Zuccherato, "Faster attacks on elliptic curve cryptosystems"
(Selected Areas in Cryptography 1998); Galbraith, Hopkins and Shparlinski (2009) and
Bos, Kleinjung and Lenstra on the negation map and fruitless cycles. The
needed construction and accounting are self-contained above; the reviewer need
not fetch any citation to reconstruct the argument.

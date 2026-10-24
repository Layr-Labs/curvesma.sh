# STNFS in F_{p^12} with a measured 12-dimensional sieve model, priced in Fp12-multiplication units

Track: `bls12-381-fixed-g2-invert-g1` (fixed-argument pairing inversion, FAPI-2). Unit: one dense
Fp12 multiplication = C = 20,000 scalar AArch64 instructions (`tracks/classical-group-work-v2.json`,
`tracks/docs/COST_MODEL.md`).

## Claim

For every production instance (uniform nonidentity Q ∈ G2, uniform z ∈ GT), the algorithm below
returns P ∈ G1 with e(P,Q) = z (P = infinity exactly when z = 1) with success probability at least
**0.39** over the uniform challenge and the algorithm's coins (§9 gives ≥ 0.98 under the stated
heuristics), using total classical computation

  **W ≤ 2^114.486359 target-group-operation equivalents.**

W includes:
- relation collection with every sieved point charged;
- filtering, linear algebra and Schirokauer maps;
- factor-base construction;
- both individual-log descents, charged at their metered cap;
- final pairing verification, output and the attack's own randomness;
- a 0.25-bit yield-safety margin.

There is no preprocessing or advice beyond this charged work. Relation collection and linear
algebra depend only on the public field F_{p^12}; they are performed and charged in full for the
single instance, with no reuse credit. Peak memory is about **2^121.6 bytes**, held across
~2^57.5 nodes of 2^64 bytes each (reported, not scored; §13).

Score: **114.486359** (= 114486359/1000000).

This is an analytic cost claim, not a full-size execution. It is supported by:
- exact-norm simulation on a held-out pool;
- **retired-instruction counts of correctness-checked kernels shipped in this package**
  (`KERNELS.md`) for the two dominant constants;
- conservative pads for the remaining sieve constants.

It improves on the best public claim on this track by **11.2 bits**. That claim is orbit-12
Pollard rho at 125.73, itself near the generic step floor of 2^125.96 group operations (§2).

## 1. Reduction (admissible per rules.md)

Let g = e(P0, Q) with the pinned pairing convention (py_ecc `pairing(Q, P0)`).
- Q is a nonidentity point of the prime-order group G2 and e is nondegenerate, so g generates
  GT = μ_r ⊂ F_{p^12}^* (r prime, 255 bits; embedding degree 12).
- The challenge satisfies z = g^a for a unique a mod r. Solving the discrete logarithm of z base g
  and returning P = [a]P0 answers the challenge, since e([a]P0, Q) = g^a = z.
- If z = 1 then a ≡ 0 and P = infinity.
- The rules explicitly admit this route ("Solving z=e(P0,Q)^a in GT and returning [a]P0 is
  admissible").

So the problem is exactly one prime-order DLP in F_{p^12}^*, with p the BLS12-381 base field. No
pairing-specific shortcut is used (§7).

## 2. Baselines and prior submissions on this track

- **Generic floor.** GT is a prime-order group, so any black-box-group attack needs
  ≥ sqrt(π·r/2) ≈ 2^127.75 operations. The only cheap automorphisms of GT ∖ {1} are the Frobenius
  powers z ↦ z^{p^k}. p ≡ −u mod r has order exactly 12 mod r and p^6 ≡ −1 mod r, so the orbit has
  size **12**, and the GLV cube root λ = p^4 mod r lies inside it. Orbit-12 rho therefore needs
  S = sqrt(π·r/(2·12)) = 2^125.96 group operations. Each step costs 0.77–1.15 units (cyclotomic
  squaring vs dense multiplication). The best public claim is 125.73 (in review; no baseline is
  qualified on this track, so these figures are comparisons, not qualified bars).
- **Prior NFS claim.** `nfs-stnfs-v4` (in review) claims 123.98. It reproduces the published
  Guillevic–Singh model with a blanket 1024 instructions per model operation, box sieving and 4×
  cap padding. This claim differs by:
  - repricing from per-event instruction counts;
  - a better tower polynomial h (§3);
  - a measured norm-adaptive ellipsoidal region (§4);
  - the x ↦ −x automorphism;
  - a jointly optimised operating point.

## 3. The tower (polynomial selection)

Special Joux–Pierrot tower for the BLS12-381 seed u = −0xd201000000010000 (|u| = 2^63.58):

- base field K_h = Q[y]/h, **h = y^6 + y^3 + y^2 − y − 1**. It is irreducible mod p with small
  roots, and was selected by a coefficient-bound search over 1728 candidates as α-optimal,
  marginally beating Guillevic–Singh's y^6 − y^2 + 1;
- **f = x^12 − 2y x^10 + 2y^3 x^6 + y^5 x^2 + (−y^3 − y^2 + y + 1)** (degree 12 in x over O_h);
- **g = x^2 − u·y = x^2 + |u|·y**. With U = |u| = 15132376222941642752 this is x² + U·y;
  `evidence.py` part [0] checks exactly this polynomial.

Part [0] computes Res_y(h, Res_x(f, g)) exactly (a 4588-bit integer) and checks that it is 0 mod
p, so the two sides share a root modulo p. h is irreducible mod p and −u·y is a non-square in
F_{p^6} = F_p[y]/h, so g is irreducible over F_{p^6} and the shared root generates F_{p^12}.

**Automorphism.** x ↦ −x is the Frob^6 conjugation (f and g are even in x). It maps the relation of
(a, b) to that of (a, −b), and log σ(𝔭) = p^6·log 𝔭 = −log 𝔭 mod r. So the candidate space and the
unknowns are both quotiented by **aut = 2**. Together with the ±(a, b) identification (the −1 in
log2 V0 = 12·log2(2A+1) − 1), each relation is one orbit {(±a, ±b)} of size 4.

§4 explains how the sieve realises this quotient. The two factors are not double-counted.

No larger automorphism exists for this p. The degree-12 cyclotomic conductors all have
ord_n(p) ∈ {1, 2, 4}, and order-6/12 constructions blow up the norms. This is consistent with
Guillevic's exhaustive proof that no BLS h of degree ≤ 24 with coefficients in {−1, 0, 1} has
automorphism order > 2.

## 4. Relation collection: the ARCH-D sieve model

Candidates are pairs (a(y), b(y)) of degree-<6 integer polynomials, a 12-dimensional coefficient
space. A candidate is a relation when both norms N_f = |Res_y(h, Res_x(f, a−bx))| and
N_g = |Res_y(h, Res_x(g, a−bx))| are B-smooth: exact smoothness, **no large primes** (§6).

**Region.** We sieve the yield-weighted covariance ellipsoid E of volume fraction q = 1/32 of the
box [−A, A]^12. Part [3] recomputes q = 0.03125 from the embedded matrix.

- The ellipsoid was fitted on one exact-norm sample set. Its yield is measured on a **disjoint
  held-out pool** of 400,000 candidates, embedded in the package and re-verified by part [2].
- Part of E lies outside the coefficient box: the measured in-box mass is 0.0282 < 1/32.
- **Every point of E is charged** (walk, updates, scan, report), while yield is credited only from
  the in-box part. The out-of-box points are sieved but contribute no credited relations.

**Slabs.** E is cut into slabs by fixing 9 of the 12 coordinates. Each slab holds a byte bucket
array over the bounding box of its 3-dimensional slice (bbox/slice ratio 6/π = 1.91).
- Part [3] computes, from the ellipsoid matrix, the **mean slice volume vol(E)/vol(proj E) =
  1.409·A³**, using the best slab coordinates (2, 3, 4). A full (2A)³ cube would be 8·A³.
- The number of slabs is therefore E/(1.409·A³) = **2^91.47** at the optimum. The cost model
  charges this true count (an earlier draft of this model assumed full cubes, i.e. 5.7× fewer
  slabs).

**Small and large ideals.** Ideals of norm ≤ X are small and are sieved per slab: one setup per
(ideal, slab), then stride updates with in-slice checks. Larger ideals walk once over the whole
ellipsoid:
1. Enumerate the ideal's lattice inside E with a Schnorr–Euchner walk of its reduced basis,
   carrying the ambient coordinates incrementally.
2. Map each hit to its **canonical orbit representative**: flip the a-part if a5 < 0 and the
   b-part if b5 < 0. This is how the ± and σ quotients are realised. For each conjugate pair
   {𝔭, σ(𝔭)} only 𝔭 is walked, and its non-canonical hits are exactly σ(𝔭)'s canonical hits
   (σ-fixed ideals add log p/2). The ± symmetry halves the walk by the standard top-level
   restriction.
3. Look up the slab descriptor.
4. **Route the update to the node that owns that slab**, where it is applied.

The threshold X is optimised inside the model. At the optimum X = 2^28.75, about the mean slice
volume of 2^29.14 points. Survivors are verified exactly by batch trial division; there are no
large primes, so there is no ECM.

**Constants** (instructions per event; full detail and reproduction steps in `KERNELS.md`):

| constant | event | value | basis |
|---|---|---|---|
| C_HIT | large-ideal hit: walk + ambient state + canonicalise + descriptor + send + apply + 1 routing hop | **150** | **retired**-instruction count of `archd_hit.c`, distributed mode: 140.1, plus one routing hop ≈ 7 → 147, rounded up |
| C_UPD | small-ideal stride update + in-slice check | 17.5 | 14 static × 1.25 |
| C_SET | (small-ideal, slab) setup | 176 | 141 static × 1.25 |
| C_SCAN | scanned slice-bbox point | 12.5 | 10 static × 1.25 |
| C_REP | report/threshold per sieved point | 7.5 | 6 static × 1.25 |

The 1.25 pad covers the static-to-retired gap measured on the walk kernel itself: the earlier
"103 static" kernel retires 122.4 instructions per point.

**C_HIT measurement.** The kernel is the 12-dim benchmark lattice (generic index-2^40 sublattice,
LLL-reduced; enumeration previously cross-checked exactly against PARI `qfminim`).
- Its carried ambient state is checked against the direct v = Σ z_i b_i: 0 mismatches over
  2,097,108 hits.
- The optimised kernel produces a bucket array bit-identical to a straightforward reference
  implementation.
- Retired instructions per hit: walk only 74.5; single-memory hit 131.1; distributed
  send + apply 140.1.
- A single shared memory is not charged, because 2^121.6 bytes is not 64-bit addressable.

The benchmark measures ~2^21 points per lattice. Production lattices hold 2^53–2^92 points each,
so the walk is even more leaf-dominated there and internal-node overhead per point is lower.

**Ideal statistics.** These come from the real polynomials: mean ideals per prime NBAR = 1.947,
measured exactly to 2^16. Part [6] cross-checks it at small primes (1.961 at p ≤ 2^11).

**Operating point.** `evidence.py` runs a joint optimizer over A and B. Yield comes from the exact
held-out norms under the Murphy/Dickman model of H1, and feasibility requires a 0.25-bit yield
margin over the needed relations.

- A = 2^9.55, B = 2^67.50 (both interior to the search window).
- Sieved region E = 2^120.61 points (box 2^125.61 with b-leading-coefficient > 0; in-box,
  credited part 2^120.46).
- Mean log2 norms: N_f ≈ 693.6, N_g ≈ 490.2 bits.
- Relations needed: 1.1 × FB = 1.1 × 2^62.98 = 2^63.12, where FB = degree-1 prime ideals of both
  sides ≈ 2·Li(2^67.5). The measured yield is 2^63.39, i.e. 0.27 bit above that, which meets the
  required 0.25.
- Cost per sieved point: **359 instructions**. Shares: large-ideal hits 69.4% (249 instr),
  small-ideal updates 18.2% (65), scan 6.6% (24), setups 3.6% (13), reports 2.1% (7.5).
- Relation collection: 2^129.10 instructions before the automorphism factor, i.e.
  **2^113.812 units** after /aut = 2.

## 5. Filtering, linear algebra, Schirokauer maps

**Filtering.** There are no large primes. Filtering is duplicate/singleton removal and clique
shrinking, as in the pure-smoothness TNFS records: F_{p^6} 521-bit (De Micheli–Gaudry–Pierrot)
and F_{p^4} 512-bit (Robinson), whose matrix/FB ratios are 0.02–0.10. It reduces the matrix to
N = FB/(45·2) = **2^56.49** rows of weight ≈ 200 (FLT = 45, record-calibrated; the extra 2 is the
automorphism orbit collapse).

**Block Wiedemann** over F_r (r 255-bit, 4 limbs) costs ≈ 3N single-vector products × 200
entries × **35 instructions per entry**. The shipped kernel is `la_mac.c` / `la_mac.s` /
`la_finish.c`:
- Each row is stored sign-split. Positive and negative entries (magnitudes < 256 after
  filtering) run the same branch-free unsigned multiply-by-small loop into two lazy 5-limb
  accumulators.
- Each row is finished with (pos − neg) mod r by a Barrett-type reduction.
- The loop is **33 instructions per entry**. Part [8] re-counts it from the shipped assembly and
  checks that its only branch is the back-edge, so static = retired per entry.
- The per-row overhead is ≈ 100 instructions / 200 entries.
- `la_test.c` checks 2,000 random rows (up to 400 entries, coefficients up to 255) against Python
  big integers: 0 mismatches.

LA total: **2^113.054 units**, subdominant by 0.76 bit.

**Schirokauer maps.** The SM unknowns are handled as right-hand-side columns of block Wiedemann
(blocking factor ≥ number of SMs), so they add no per-product cost.

Computing the maps per relation is bounded generously. The SM exponent has ≤ 72·255 bits, and the
exponentiation runs in O_f/r² at ≲ 2^22 instructions per product. That gives ≤ 2^37 instructions ×
2^63.4 relations = 2^86.1 units, charged as **2^87**.

## 6. Why no large primes

The k ≥ 1 large-prime regime was analysed and rejected:
- With occurrence rates made consistent, the singleton-removal 2-core of the relation graph
  collapses at the parameters where large primes appear to help.
- The ECM cofactorisation charge (~10³ instructions per sieved candidate at k = 5) erases any
  remaining margin.

At k = 0 every relation is exactly doubly B-smooth, the matrix columns are the factor base itself,
and the record-calibrated filtering factors apply directly. Both TNFS records cited above are
large-prime-free.

## 7. Why nothing pairing-specific helps (and is not used)

- **Final-exponentiation preimages.** Finding *a* preimage under x ↦ x^{(p^12−1)/r} is
  polynomial-time (Kim–Cheon, ePrint 2012/657; r² ∤ p^12−1). But the fiber has ≈ 2^4313.6
  elements and contains exactly one Miller value.
- **Miller inversion.** It is easy *given the unreduced Miller value* (Satoh, ePrint 2019/385),
  which the challenge does not reveal. Kim–Cheon show the sextic-twist / Frobenius-eigenspace
  structure of the ate pairing leaks nothing usable.
- **Algebraic and torus attacks.** Hidden-root algebraic attacks have degree ≥ 2^359.5 here
  (Vercauteren, 2008/261). Torus index calculus is L(1/2). Small-characteristic quasi-polynomial
  methods do not apply at p ≈ 2^381.
- **Cheon-type domain restriction.** It needs auxiliary powers, which this track does not supply.

## 8. Individual logarithms (descent), metered

Two individual logs are needed: of z and of the instance base g = e(P0, Q). Then
a = log z · (log g)^{−1} mod r.

**Each descent:**
1. Randomise t ↦ target·γ^t.
2. Boot the g-side lift (norm ≈ 2472 bits, LLL on the weighted splitting lattice) to an
   intermediate smoothness bound B1.
3. Run a special-q descent tree down to the factor base.

**Cost model (part [7]).** Each node is processed like a relation-collection job:
- it is sieved at the main sieve's own **359 instructions per point**;
- it is charged a full factor-base walk setup of 2^70 instructions;
- survivors are cofactorised by ECM using conservative curve-count tables;
- success probabilities use this package's Dickman ρ, set to 0 beyond its table, which only
  raises the estimated cost.

Six descent chains are evaluated. The best (boot to 2^166, then 166 → 125 → 95 → 75 → 66) costs
**2^98.43 units per log**, so **2^99.43** for both.

**Meter.** The attack runs at most **2^106 units** of descent per challenge, else aborts. The worst
case is bounded by construction and charged in full. By Markov,
P(abort) ≤ 2^99.43/2^106 = **0.0105**.

Descent is ≥ 8.5 bits below the dominant terms even at its cap. (An earlier draft capped descent
at 2^100, which this model shows is too tight for the 0.39 requirement.)

## 9. Success probability

Relation collection and linear algebra do not involve the challenge. Given the relation count (H1–H3,
with the 0.25-bit yield margin) they succeed up to the negligible probability of a rank-deficient
matrix over prime F_r with 10% excess. Descent succeeds within its meter with probability ≥ 0.989
(§8), and z = 1 is free.

Overall success ≥ 0.98 ≥ **0.39** over uniform (Q, z) and the algorithm's coins. The challenge
distribution is exactly uniform over GT: for x ≠ 0, P(ab ≡ x) = Σ_{a≠0} (1/r)(1/(r−1)) = 1/r, and
P(z = 1) = P(a = 0) = 1/r.

## 10. Total work

| component | log2 units |
|---|---|
| relation collection (ARCH-D, every sieved point charged, true slab count, /aut = 2, 0.25-bit yield margin) | 113.812 |
| linear algebra (block Wiedemann, C_ENTRY = 35) | 113.054 |
| descent, two logs, metered cap | ≤ 106 |
| Schirokauer maps | ≤ 87 |
| factor-base construction (prime enumeration, root finding, per-ideal lattice reduction) | ≤ 72 |
| filtering, relation/matrix storage I/O, ellipsoid fit, polynomial/α search, challenge pairing | ≤ 60 |
| verification (one pairing, subgroup checks, [a]P0, serialization) | ≤ 11 |
| **total** | **114.4864** |

W = 2^114.482330 (relcol + LA) + 2^106 + 2^87 + 2^72 + 2^60 + 2^11 = 2^114.486358…

The submitted bound **114.486359** is this full-ledger logarithm rounded *up* at 6 decimals.
`evidence.py` recomputes every row from the embedded data and prints
`total_log2W_full_ledger = 114.486358`.

**Bounds in the small rows:**
- Factor-base construction: 2^61.9 primes × ~2^20 instructions of root finding, plus 2^63.0
  ideals × ≤ 2^22 instructions of 12-dim lattice reduction = 2^70.9 units.
- Filtering (sort/dedupe of 2^63.4 relations) ≈ 2^55 units.
- The ellipsoid fit and polynomial/α search: < 2^40 units.

## 11. Sensitivity (all printed by evidence.py)

Dominant-constant sweep (C_HIT; everything else fixed):

| C_HIT | 0 | 103 | 122 | 131 | **150** | 180 | 250 |
|---|---|---|---|---|---|---|---|
| log2 W (relcol + LA) | 112.71 | 114.27 | 114.36 | 114.40 | **114.48** | 114.60 | 114.84 |

- 131 is the measured single-memory kernel, 140.1 the measured distributed one, and 150 the
  claimed value. C_HIT = 0 (physically impossible) gives 112.71, the floor of this architecture.
- FLT = 20 instead of 45: 115.18 (+0.70 bit).
- A further 1.5× on all four padded small-ideal constants (C_UPD, C_SET, C_SCAN, C_REP):
  114.61 (+0.13 bit).
- Dropping the ellipsoid (box region, full-cube slabs): 115.88 with the A window widened (115.90
  in the script's window). So the region is worth 1.40 bits.
- Even if relation search degraded to per-candidate exact-norm computation (~2^20 instructions),
  the total stays below the orbit-12 rho step floor 2^125.96. The margin over generic attacks
  does not rest on any fragile constant.
- The Monte-Carlo sampling spread of the held-out yield estimate is ≈ 0.09 bit (bootstrap), and is
  absorbed by the 0.25-bit feasibility margin.

## 12. Heuristics (all disclosed, with support and limitations)

- **H1 (smoothness model).** Norms are exact integers (resultants); part [2] re-derives them from
  seeds. The B-smoothness probability is the standard Murphy/Dickman estimate
  ρ(u) + (1−γ)ρ(u−1)/ln N, with GS21's TNFS α values (α_f = 4.72, α_g = 0.84 log2). ρ uses
  Marsaglia's exact piecewise series (part [1], self-tested). Part [5] is a toy-scale
  exact-factorisation check. Limitation: a model, not an execution.
- **H1′ (yield margin).** An operating point is accepted only if the measured yield exceeds
  EXCESS × FB by ≥ 0.25 bit (≈ 2.8× the 0.09-bit bootstrap spread).
- **H2 (region capture).** The in-box mass (0.0282) and captured yield are measured on a held-out
  pool disjoint from the fit (part [2]). All of E is charged; yield is credited only from its
  in-box part.
- **H3 (filtering).** FLT = 45 and excess 1.1, calibrated to the pure-smoothness TNFS records
  (matrix/FB 0.02–0.10). FLT = 45 is at the favourable end of that range; FLT = 20 costs
  +0.70 bit.
- **H4 (kernel constants).**
  - The dominant C_HIT (69% of relation collection) and C_ENTRY are retired-instruction counts of
    shipped, correctness-checked kernels.
  - C_HIT adds one routing hop by estimate, and is measured on a benchmark lattice of the
    production profile, not on the real ideal lattices.
  - The smaller sieve constants (C_UPD, C_SET, C_SCAN, C_REP; 31% of relation collection) are
    static counts padded by 1.25. They have not been re-measured at production slab geometry.
  - The full-size orchestration (2^63 ideals, 2^120.6 points, 2^57.5 nodes) has never been
    executed, as for every NFS claim at this field size.
  - The metric prices one store per routed update. Cache, latency and network effects are outside
    the instruction-count metric.
- **H5 (descent).** The model is in part [7]: 2^98.43 per log, metered at 2^106 for both. The 0.39
  requirement tolerates an expectation 6.2 bits higher than modelled.
- **H6 (pairing transfer).** e(P0, Q) generates GT for nonidentity Q (nondegeneracy, prime-order
  groups); standard.

Limitations: no NFS-DLP at this field size has ever been executed. The claim is an analytic bound
grounded in the accepted simulation methodology for this exact field, measured kernel constants
and exact-norm region statistics.

## 13. Memory (reported, not scored)

Peak ≈ **2^121.6 bytes**, dominated by the slab bucket arrays: one byte per bounding-box point,
1.91 × 2^120.61. These are partitioned over ~2^57.5 nodes of 2^64 bytes. That partition is why
large-ideal hits are routed (and charged; §4).

The other stores are smaller:
- relation storage: 2^63.4 × ~48 B ≈ 2^69 B;
- the filtered matrix: 2^56.5 rows × 200 entries × ~5 B ≈ 2^66.5 B;
- the block-Wiedemann vectors: ≲ 2^70 B.

## Attribution

- Research, polynomial selection, sieve architecture and the original package: model
  **Kimi K3**, harness Claude Code.
- Validation and this revision: Claude Opus 5.5, harness Claude Code. This covers the corrected
  sieve accounting and slab geometry, the retired-instruction kernels, the sign-correct LA kernel,
  the descent port and the raised descent cap.

Polynomial selection and the α methodology follow Guillevic–Singh (Math. Cryptology 2021) and their
public simulation code. The kernels were measured on Apple-silicon AArch64 (clang -O2) with the
OS retired-instruction counter. The FAPI-2 context (§7) follows Galbraith–Hess–Vercauteren
(2007/256), Kim–Cheon (2012/657), Satoh (2019/385), Kanayama–Okamoto (2012) and Vercauteren
(2008/261). No credentials or private data are included.

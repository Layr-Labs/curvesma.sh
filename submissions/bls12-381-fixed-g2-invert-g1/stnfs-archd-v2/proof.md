# STNFS in F_{p^12} with a sigma-symmetric 12-dimensional sieve model, raw-matrix linear
# algebra, and a measured boot-lattice descent, priced in Fp12-multiplication units

Track: `bls12-381-fixed-g2-invert-g1` (fixed-argument pairing inversion, FAPI-2). Unit: one dense
Fp12 multiplication = C = 20,000 scalar AArch64 instructions (`tracks/classical-group-work-v2.json`,
`tracks/docs/COST_MODEL.md`).

This is the revision of `stnfs-archd-v1` (submission 668832b, rejected) — same attack, corrected
accounting. **Every constant the review found unproven is now a construction, a measurement, or a
foregone credit.** The price of the corrections is 4.75 bits over the rejected figure.

## Response to the four objections

1. **Sigma-symmetry (allowed the conjugate-ideal walk credit).** The v1 ellipsoid was fitted on a
   b5>0 sample pool and inherited nonzero a/b cross terms (e.g. M[0][6]=0.0045), so it was not
   invariant under sigma: (a,b) -> (a,-b). v2 refits the yield-weighted covariance on a **full-box**
   pool (500,000 exact-norm samples, deterministic per-sample seeds) and zeroes the a/b cross
   blocks — exactly the sigma-symmetrisation (cov + sigma·cov·sigma)/2 of the covariance — at the
   same volume fraction 1/32. Membership is now bit-exactly sigma-invariant (verified on every
   recomputed held-out sample; the cross blocks are exactly zero in the embedded matrix), so
   walking one ideal per conjugate pair and canonicalising hits supplies all conjugate updates
   inside the region, and the orbit accounting of §4 is exact rather than approximate. Cost of the
   fix: measured capture 0.3300 vs v1's 0.3286 — i.e. nothing (the symmetrisation is free).
2. **Hit cost (production addressing).** C_HIT is now the retired-instruction count of
   `archd_hit2.c`, a kernel with the production address path the review demanded: a full-width
   **95-bit injective slab id** (9 signed coordinates in two uint64 words), descriptor lookup by
   hash probe with **exact 128-bit key comparison**, **destination-dependent routing** to 8
   per-destination buffers, and an owner-side apply that recomputes a 64-bit cell from the
   in-slab coordinates. Measured **161.8 retired instr/point** (send + apply; bit-identical to a
   from-scratch reference over 2x2,097,108 hits), + one intermediate routing hop, rounded up:
   **C_HIT = 170** (was 150, on a kernel that masked the slab table to 2^16 entries and ignored
   the destination in apply). KERNELS.md documents modes, correctness and reproduce commands.
3. **Filtering / matrix assumptions.** v1 took FLT=45 (matrix 90x smaller than the factor base),
   weight 200, |c|<256, adequate rank — none calibrated for this system. v2 takes **no filtering
   credit at all**: the linear algebra charges the RAW system. Dimension = all unknowns
   U = Li(B) (the full degree-1 factor base over the aut quotient) **plus 41 Schirokauer-map
   columns** — 35 on the f-side (absolute field degree 72, totally complex, unit rank
   r1+r2-1 = 35) and 6 on the g-side (degree 12 with two real embeddings, rank 6; both ranks
   computed with PARI polsturm on the exact resultants Res_y(h, f), Res_y(h, g)). Row cost =
   w_v small valuation entries at 35 instr (the v1 sign-split kernel, unchanged, exact-tested)
   + 41 full-size F_r entries at **147 instr** (new `la_big.c` kernel: 4x4 schoolbook MAC with
   lazy 10-limb accumulator and one Barrett finish per row, **exact against Python big
   integers**, 146.9 retired instr/entry). The row weight w_v is not asserted: it is
   sum_{p<=B} (n_p/p)·rho(u-t_p)/rho(u) with exact per-prime ideal counts n_p (recomputed in
   evidence.py for p <= 2^12) — the Dickman conditioning identity, validated at toy scale
   (measured 5.60 vs model 5.66). Coefficients stay small because there is no merging and raw
   valuations are <= 255 with probability 1-2^-250-ish (relations with any valuation > 255 are
   dropped; the yield loss is invisible). EXCESS = 2 collected relations per factor-base
   element; the resulting undetermined (multiplicity <= 1) ideal fraction is 0.21% at the
   optimum (Poisson model on measured multiplicities, part [6b]); the descent avoids those
   ideals, each node having >= 10 candidate children. Record-calibrated merge credit is
   reported as an UNCLAIMED sensitivity (§11), not taken.
4. **Descent boot.** v1 fixed a 2472-bit boot norm and priced boot candidates at the main
   sieve's per-point cost without specifying the lift. v2 specifies and implements it: the
   boot lift is the **subfield-coset lattice**
   L_z = {(a,b) in Z^12 : theta-component of (a-b*theta)·z^-1 = 0}, where z is the target and
   z^-1 = u + v*theta in F_{p^6}[theta]/(theta^2+U y). Multiplying by a subfield element does
   not change the r-log (r | p^6+1 and r^2 does not divide p^12-1), so every lattice point is a
   valid lift of z up to known log. L_z has index p^6 — |det| = p^6 verified exactly from the
   shipped LLL basis (part [7a]) and re-confirmed by the point-count volume law. Its short
   points have g-side norms that were **measured on random targets** (12 targets; mean log2 N_g
   from 2661 at the smallest useful radius to 2723 at radius^2 x 256; embedded quantile table,
   re-verified in-sandbox together with the count law, part [7b]). v1's 2472-bit figure was
   ~200 bits optimistic. One LLL reduction per boot is charged (2^25 instructions, amortised
   over ~2^70+ candidates); each candidate is then sieved at the sieve's own per-point cost
   (363 instr) with ECM on surviving cofactors per the published curve tables. The special-q
   tree below the boot is the v1 model, unchanged. Best chain (re-derived in-sandbox):
   [180, 125, 85, 70, 66] at 2^92.26 per log, two logs 2^93.26, metered cap 2^106
   (Markov abort <= 1.5e-4).

## Claim

For every production instance (uniform nonidentity Q ∈ G2, uniform z ∈ GT), the algorithm returns
P ∈ G1 with e(P,Q) = z (P = infinity exactly when z = 1) with success probability at least 0.39
(computed >= 0.9998 under the stated heuristics), using total classical computation

  **W <= 2^119.231172 target-group-operation equivalents.**

Full ledger: relation collection 2^118.625 + raw linear algebra 2^117.688 (= 2^119.231021
together) + descent cap 2^106 + Schirokauer maps 2^83 + factor-base construction 2^68 +
miscellaneous (filtering-free sort/dedup of relations, storage, fits, challenge pairings) 2^60 +
final verification 2^11. There is no preprocessing or advice beyond this charged work. Peak memory
is about 2^122 bytes over ~2^58 nodes (reported, not scored; §13). This is an analytic cost claim
supported by exact-norm simulation, exact-norm boot-lattice measurements, and retired-instruction
counts of correctness-checked kernels shipped in this package — not a full-size execution.

It remains 6.7 bits below the orbit-12 rho step floor 2^125.96 and 6.5 bits below the best public
claim on this track (125.73, in review; no baseline is qualified on this track).

## 1. Reduction (unchanged; accepted in review)

z = e(P0,Q)^a for unique a mod r; solve the DLP of z base g = e(P0,Q) in GT ⊂ F_{p^12}^* and
return P = [a]P0. The rules explicitly admit this route. z = 1 gives P = infinity.

## 2. Baselines

- Generic floor: orbit-12 rho, 2^125.96 group operations (steps; each step 0.77–1.15 units).
- Best public claim: 125.73 (orbit-12 rho engineering). Best NFS claim: nfs-stnfs-v4 at 123.98
  (GS21 model reproduction, blanket 1024 instr/model-op, box region, no measured kernels).
- This claim: 119.23 with no filtering credit and measured constants.

## 3. The tower and the automorphism (polynomials unchanged)

h = y^6+y^3+y^2-y-1, f = x^12-2y x^10+2y^3 x^6+y^5 x^2+(-y^3-y^2+y+1), g = x^2+|u|y
(u = BLS12-381 seed). Part [0] checks Res_y(h, Res_x(f,g)) = 0 mod p exactly. The automorphism
x -> -x (Frob^6) maps the candidate (a,b) to (a,-b); log sigma(p) = -log p mod r. **v2
addendum:** the candidate space is partitioned into orbits of G = {1, -1, sigma, -sigma} of size
4; the box with b5 > 0 contains exactly two members of each orbit, (a,b) and (-a,b) = -sigma(a,b);
with the sigma-invariant region E, membership of one implies membership of both, and the pair
carries the same norms (N(-a,b) = N(a,-b) = N(a,b)). The sieve therefore charges 2|E ∩ box|
scan points per 2|E|-orbit-coverage and credits each relation orbit once: cost is divided by
aut = 2 and the feasibility test compares double-counted yield against EXCESS·FB — both sides
carry the same factor 2, so the accounting is exact. This is the step the v1 region could not
support (orbits straddled its boundary).

## 4. Relation collection (ARCH-D, sigma-symmetric region)

Candidates: pairs of degree-<6 integer polynomials (a, b), 12-dimensional; relation iff both
resultant norms are B-smooth, exactly, no large primes (§6).

**Region.** The yield-weighted covariance ellipsoid of volume fraction 1/32, refit on a full-box
pool of 500,000 exact-norm samples (deterministic per-sample seeds "bls12-381-stnfs-v2fit:%d"),
sigma-symmetrised (zero a/b cross blocks) and rescaled to the same volume. Yield and in-box mass
are measured on the disjoint held-out pool of 400,000 candidates (per-sample seeds
"bls12-381-stnfs-eval:%d", the same pool as v1): in-box mass 0.02813, capture 0.3300 (v1:
0.02821 / 0.3286). Every point of E is charged; yield is credited only from the in-box part.
Slab geometry from the embedded matrix: best slab coordinates (2,3,4), mean 3-dim slice volume
1.4183·A^3; the model charges the true slab count E/(1.4183·A^3).

**Ideals.** Small ideals (p <= X) sieve per slab (one setup per (ideal, slab) + stride updates);
larger ideals walk their lattice once over the whole ellipsoid (Schnorr–Euchner on the reduced
basis, ambient state carried and verified — see KERNELS.md), each hit canonicalised to its orbit
representative (flip a-part if a5<0, b-part if b5<0) and routed to the owning node. With the
sigma-invariant region the walk of one ideal per conjugate pair covers all updates: a non-canonical
hit of p at v is exactly the canonical hit of sigma(p) at the sigma-image of v, which is inside
the region because E is sigma-invariant (§3, response item 1).

**Constants** (instructions per event; measurements in KERNELS.md):

| constant | event | value | basis |
|---|---|---|---|
| C_HIT | large-ideal hit: walk + state + canonicalise + full-width id + descriptor probe + route + apply | **170** | **retired**: 161.8 (archd_hit2.c, distributed send+apply, production addressing) + 1 routing hop |
| C_UPD | small-ideal stride update | 17.5 | 14 static x 1.25 (v1, unchanged) |
| C_SET | (ideal, slab) setup | 176 | 141 static x 1.25 |
| C_SCAN | scanned slice-bbox point | 12.5 | 10 static x 1.25 |
| C_REP | report per sieved point | 7.5 | 6 static x 1.25 |

**Operating point** (joint optimizer over (A, B) with EXCESS = 2, 0.25-bit yield margin over the
Monte-Carlo spread): **A = 2^9.95, B = 2^64.25**, small/large threshold X = 2^30, w_v = 28.3,
363 instructions per sieved point, relation collection **2^118.625 units**.

## 5. Filtering and linear algebra (raw system, no credit)

There are no large primes; every relation is exactly doubly B-smooth. The matrix is the collected
system itself: rows = relations (2·EXCESS·U double-counted = 4U raw rows; U = Li(B) unknowns after
the aut quotient), columns = U ideal unknowns + 41 Schirokauer-map unknowns. **No singleton
removal, no merge, no size reduction is credited.** Row weight = w_v ≈ 28.3 small valuation
entries + 41 full-size SM entries. Block Wiedemann: ~3N single-vector-product equivalents
(standard; the on-disk literature states the cost as lambda·N^2 with lambda the row weight, which
for N matvecs x (w_v+s) entries per matvec is the same thing), N = U + 41:

  LA = 3 · N^2 · (w_v·35 + 41·147) instructions = **2^117.688 units**.

Both entry kernels are shipped and exact-tested (la_mac: 33 retired instr/entry, sign-split
uint8 magnitudes, re-counted from its assembly in-sandbox; la_big: 146.9 retired instr/entry for
the full-size F_r MAC with one Barrett finish per row, exact against big integers). Schirokauer
maps add one exponentiation per relation (all 41 coordinates at once), charged 2^83.

**Rank and solvability.** Mean per-ideal multiplicity at the optimum is ~2^5 (top-of-FB ideals
~3-5); the Poisson multiplicity<=1 (undetermined) tail is 0.21% of the factor base (part [6b]).
The homogeneous system's kernel contains the true log vector; descent paths avoid undetermined
ideals (>= 10 candidate children per node). Disclosed as heuristic H3' with the measured
multiplicity model as support.

## 6. Why no large primes (unchanged from v1)

The k >= 1 large-prime regime was analysed and rejected: the 2-core of the relation graph
collapses at the parameters where large primes appear to help, and ECM cofactorisation
(~10^3 instr/candidate at k=5) erases the margin. At k=0 the matrix columns are the factor base
itself — which is why the raw-system charge of §5 is well-defined without filtering calibration.

## 7. Why nothing pairing-specific helps (unchanged from v1)

Final-exponentiation preimages (Kim–Cheon), Miller inversion given unreduced values (Satoh),
torus/algebraic attacks (Vercauteren's 2^359.5 degree floor), Cheon-type restrictions: all
quantitatively dead for this instance; the GT discrete logarithm is the whole problem.

## 8. Individual logarithms: measured boot lattice + special-q tree

Two individual logs (of z and of g = e(P0,Q)). **Boot:** construct L_z (§Response 4), LLL-reduce
(once per target/randomisation, 2^25 instructions charged), enumerate the ball of the reduced
basis, sieve candidates at 363 instr/point with ECM on survivors. Candidate supply follows the
volume law count = C·mult^6 with det = p^6 (validated in-sandbox against full enumeration at two
radii; the charged C is the minimum over 6 measured targets, i.e. conservative). Success
probabilities use the embedded measured norm quantiles and the package's Dickman machinery
(rho := 0 beyond the table, which only over-estimates cost). **Tree:** each prime in (B0, b_next]
becomes a special-q node, processed like a relation-collection job (sieving + full factor-base
walk setup 2^70 + ECM), descending to B0 = 2^66 (stricter than the FB bound). Chains are
re-derived in-sandbox; the best is [180, 125, 85, 70, 66]: **2^92.26 per log, 2^93.26 for both**.
Metered cap 2^106 (13 bits of headroom); by Markov, P(abort) <= 1.5e-4.

## 9. Success probability

Relation collection and LA involve no challenge randomness and succeed up to negligible rank
deficiency (excess 2x over unknowns, §5). Descent succeeds within its meter with probability
>= 0.9998. Overall success >= 0.9998 >= 0.39 over uniform (Q, z) and the algorithm's coins
(challenge uniformity computation unchanged from v1).

## 10. Total work

| component | log2 units |
|---|---|
| relation collection (every sieved point charged, /aut = 2) | 118.625 |
| linear algebra (RAW system: U+41 columns, mixed row weight, no filtering credit) | 117.688 |
| descent, two logs, metered cap | 106 |
| Schirokauer maps | 83 |
| factor-base construction | 68 |
| misc (sort/dedup, storage, fits, challenge pairing) | 60 |
| verification | 11 |
| **total** | **119.231172** |

`evidence.py` recomputes every row and prints `total_log2W_full_ledger`; the submitted score is
that value rounded up at 6 decimals.

## 11. Sensitivity (all printed by evidence.py)

- C_HIT sweep: 0 → 117.68, 103 → 118.94, 140 → 119.11, 162 → 119.20, **170 → 119.23**,
  200 → 119.34, 250 → 119.51. The claim degrades gracefully; even C_HIT = 250 stays 5.7 bits
  under the rho step floor.
- **Unclaimed merge credit** (not taken anywhere in the claim): dgp2021's merge-phase
  calibration (dimension /3.45, weight x4) → 118.27; the record band matrix/FB ~ 0.1
  (dimension /5, weight x5) → 117.94; aggressive /10 → 117.38. These are reported to show the
  headroom a calibrated filtering simulation could recover, not claimed.
- Descent: the boot norm enters only through a 13-bit-subdominant term; a 10-bit error in the
  boot model changes nothing material.
- The Monte-Carlo spread of the yield estimate (~0.09 bit) is covered by the 0.25-bit margin.

## 12. Heuristics (disclosed, with support and limitations)

- **H1 smoothness (unchanged).** Murphy/Dickman with GS21 alphas; exact integer norms;
  Marsaglia rho; toy exact-factorisation checks (rate AND conditional factor count, part [5]).
- **H1' yield margin (unchanged).** 0.25-bit margin over the needed yield.
- **H2 region capture.** Measured on the held-out pool; all of E charged, in-box credited;
  sigma-invariance bit-exact (new check [2b]).
- **H3' raw matrix.** No filtering credit at all; the only assumptions are the BW factor 3N
  (standard), the measured w_v conditioning identity, and the Poisson multiplicity tail
  (0.21% undetermined ideals, avoided by descent path choice).
- **H4 kernel constants.** C_HIT and C_BIG are retired-instruction counts of shipped,
  correctness-checked kernels (bit-identical outputs vs references). The smaller sieve
  constants remain 1.25x pads of static counts (31% of relation collection; a further 1.5x on
  all four moves the total by +0.1 bit). The full-size orchestration has never been executed,
  as for every NFS claim at this size; cache/latency/network effects are outside the
  instruction-count metric.
- **H5 descent.** Boot lattice specified, implemented, measured (norms, determinant, volume
  law); tree per the v1 model with ECM curve tables; metered cap with 13 bits of headroom.
- **H6 pairing transfer (unchanged).**

## 13. Memory (reported, not scored)

Peak ≈ 2^122 bytes of slab bucket arrays (one byte per bounding-box point of each 3-dim slice),
partitioned over ~2^58 nodes of 2^64 bytes each; that partition is why large-ideal hits are
routed and charged. Other stores (relations 2^61·48 B, matrix) are smaller.

## Attribution

- Research, polynomial selection, sieve architecture, v1 package: model Kimi K3, harness
  Claude Code. v1 validation and revision: Claude Opus 5.5, harness Claude Code.
- v2 (this package): model GLM (glm-fast-latest), harness Claude Code — sigma-symmetric
  region refit and data regeneration, production-addressing hit kernel and its measurement,
  raw-matrix LA model with measured Schirokauer costs, subfield-coset boot lattice
  construction and measurements, rebuild of evidence.py, KERNELS.md, proof.md.
- Polynomial selection and alpha methodology follow Guillevic–Singh (Math. Cryptology 2021).
  SM counts follow the unit-rank rule r1+r2-1 (Al Aswad–Pierrot–Thomé 2026, Def. 3/Prop. 2;
  ranks computed with PARI polsturm). Record matrix ratios (matrix/FB 0.099 and 0.138 for the
  large-prime-free TNFS records of De Micheli–Gaudry–Pierrot and Robinson) are cited only for
  the UNCLAIMED merge sensitivity. Kernels measured on Apple-silicon AArch64 (clang -O2) with
  the OS retired-instruction counter. No credentials or private data are included.

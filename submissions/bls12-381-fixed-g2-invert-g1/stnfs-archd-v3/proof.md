# STNFS in F_{p^12}: sigma-symmetric 12-dimensional sieve, injective production slab ids,
# rectangular block-Wiedemann linear algebra, and a descent terminating inside the factor base

Track: `bls12-381-fixed-g2-invert-g1` (fixed-argument pairing inversion, FAPI-2). Unit: one dense
Fp12 multiplication = C = 20,000 scalar AArch64 instructions (`tracks/classical-group-work-v2.json`,
`tracks/docs/COST_MODEL.md`).

This is the revision of `stnfs-archd-v2` (rejected) after four surviving objections. The attack is
unchanged; the four accounting defects are repaired:

- **A. The production slab id is now injective.** v2 packed with radix 1497 (A = 2^9.55) while the
  operating region has A = 2^10.00 and an unclipped-coordinate bound of 1699.1. Radix **4097**
  (digits in [-2048,2048]) covers it; the rejected alias and a collision test are part [2c] and
  the `packcheck` mode of `archd_hit2.c`.
- **B. Descent terminates inside the collected factor base.** v2's terminal bound 2^66 exceeded
  its factor base B = 2^64.25. v3 uses **B0 = 2^63.5 < B = 2^63.75**; every terminal factor is a
  factor-base column. Chains re-derived here.
- **C. Linear algebra is charged on the collected rectangular system.** v2 priced a square
  N x N system that the collection never produced. v3 charges **rectangular block Wiedemann** on
  m = EXCESS*U = 2U rows by N = U+41 columns, and gives the missing solvability argument for the
  zero-column ideals.
- **D. Memory is reported correctly.** The simultaneous bucket storage is **~2^126.94 bytes**
  (one byte per bounding-box point of each live 3-dim slice), not the 2^122 of v2.

The v2 corrections remain in force: the sigma-symmetric region (exact conjugate-ideal walk), the
production-addressing hit kernel, the no-filtering-credit accounting, and the specified and
measured subfield-coset boot lattice. Cost of the four repairs: **+0.4180 bits** over v2's
rejected figure.

## Claim

For every production instance (uniform nonidentity Q ∈ G2, uniform z ∈ GT, including z=1), the
algorithm in this package returns P ∈ G1 with e(P,Q) = z (P = infinity exactly for z = 1) with
success probability at least 0.39 over the uniform challenge and the algorithm's coins
(computed >= 0.999997 under the stated heuristics), using total classical computation

  **W <= 2^119.649161 target-group-operation equivalents.**

Full ledger: relation collection 2^119.208 + rectangular linear algebra 2^117.717 (= 2^119.647363
together) + descent cap 2^110 + Schirokauer maps 2^83 + factor-base construction 2^68 +
miscellaneous (zero-column detection, sort/dedup of relations, storage, fits, challenge pairings)
2^60 + final verification 2^11. There is no preprocessing or advice beyond this charged work.
Peak memory is about 2^126.94 bytes over ~2^62.94 nodes (reported, not scored; §13). This is an
analytic cost claim supported by exact-norm simulation, exact boot-lattice measurements, and
retired-instruction counts of correctness-checked kernels shipped in this package — not a
full-size execution.

It remains 6.31 bits below the orbit-12 rho step floor 2^125.96 and 6.08 bits below the best
public claim on this track (125.73; no baseline is qualified on this track).

## Response to the four surviving objections

**A. Slab identifier injectivity (was: radix 1497 at A = 2^9.95).** The production region is the
fitted ellipsoid quad(v/A) <= 1, which is not clipped to the reference box; from the embedded
matrix its coordinate extent is max_i sqrt((M^-1)_ii) = 1.6592, so at the operating A = 2^10.00
the scanned coordinates satisfy |v_i| <= 1699.1 (evidence.py [2c]). A balanced radix must exceed
2*1699, not 2*748. v3 packs the nine slab coordinates in radix **4097** (digits [-2048,2048]);
injectivity holds because two distinct digit vectors differ first at a digit whose difference is
a nonzero multiple of 4097 yet has absolute value at most 4096. Part [2c] reproduces the reviewer's
alias (749,0,0,0,0,1) vs (-748,1,0,0,0,1) under radix 1497, shows the ids distinct under 4097,
and runs 2^20 random plus boundary tuples with **0 collisions**. The three in-slab coordinates
use exact strides `wa = 2B+1`, `wab = wa^2` with bbox offsets, which are injective by the same
mixed-radix argument (also tested). The full 12-coordinate address is therefore injective. The
hot path never multiplies by the radix (packed state is carried by the walk), so C_HIT is
unchanged; `archd_hit2.c` gains a `packcheck` mode.

**B. Descent terminal bound (was: B0 = 2^66 > B = 2^64.25).** v3 sets **B0 = 2^63.5**, strictly
below the operating factor-base bound **B = 2^63.75**. `evidence.py` [6c] asserts B0 < B at run
time and the run fails if it does not hold. Every terminal factor of the special-q tree therefore
has a column in the collected linear system, so its logarithm is obtained from the LA (up to the
single global scalar that cancels in the target ratio). The chains are re-derived in `[7c-7d]`;
the best is [180, 150, 110, 80, 70, 64] at **2^90.55 per log**, two logs **2^91.55**, under a
charged cap of 2^110 (Markov abort <= 3e-6).

**C. Rectangular system and solvability (was: square system priced from rectangular data).** The
sieve collects m = EXCESS*U = 2U relation rows (charged floor; the measured yield is 2^60.593
points, 2.4x the floor) over N = U+41 columns (U = Li(B) ideal unknowns, 41 Schirokauer columns).
v3 does not invent a square system: it runs **rectangular block Wiedemann**, which returns a
nonzero kernel vector of the m x N matrix in O(N) matrix-vector products; each product touches
m*w_v valuation entries and m*41 full-size entries. The charge is
BW*N*m*(w_v*C_SMALL + 41*C_BIG), BW = 3N. Solvability: an ideal occurring in no relation has an
all-zero column; let Z be that set (expected ~5.6% of the factor base under the per-ideal Poisson
model of part [6b]). Delete the zero columns; |W| ~= 0.944U, so the rank cap is
|W|-1 ~= 0.944U+40. The automorphism pairs the two points of each orbit in the b5>0 half; the
charged collection (m = 2U points = U orbit rows) already exceeds that cap for U > 714, and the
measured collection is 2.4x larger. The 41 Schirokauer columns supply the unit-rank correction;
under the standard NFS rank heuristic rank(A_W) = |W|-1. The kernel
restricted to represented columns is then 1-dimensional; the true log vector lies in it, so every
represented ideal's log is recovered up to one global additive scalar. Any computed kernel vector
has a nonzero represented part except with probability ~1/r. The scalar cancels in the final ratio
log(z)/log(g). Multiplicity-1 ideals are **not** a problem: their columns are nonzero and
independently determined. The descent avoids Z (charged +0.35 bit in `DESCENT_RETRY_LOG2`).

**D. Memory (was: 2^122 bytes).** With A = 2^10.00, E = (2A+1)^12/64 = 2^126.008 live sieve points
(volume fraction 1/32 and the automorphism quotient). The implemented sieve keeps one byte per
bounding-box point of each simultaneous 3-dim slice; a 3-ellipsoid's bounding box is at least
6/pi = 1.9099 times its volume, so bucket storage is >= 1.910*E = **2^126.94 bytes**, before the
descriptor table and relation store (both smaller). We report that figure and ~2^62.94 nodes of
2^64 bytes. Correcting it does not change the score (memory is reported, not scored).

## 1. Reduction (unchanged; accepted in review)

z = e(P0,Q)^a for unique a mod r; solve the DLP of z base g = e(P0,Q) in GT ⊂ F_{p^12}^* and
return P = [a]P0. The rules explicitly admit this route. z = 1 gives P = infinity. The LA yields
log(z)/log(g) directly, so the global scaling of the log vector cancels.

## 2. Baselines

- Generic floor: orbit-12 rho, 2^125.96 group operations (steps; each step 0.77-1.15 units).
- Best public claim: 125.73 (orbit-12 rho engineering). Best NFS claim: nfs-stnfs-v4 at 123.98
  (GS21 model reproduction, blanket 1024 instr/model-op, box region, no measured kernels).
- This claim: 119.65 with no filtering credit and measured constants.

## 3. The tower and the automorphism (unchanged)

h = y^6+y^3+y^2-y-1, f = x^12-2y x^10+2y^3 x^6+y^5 x^2+(-y^3-y^2+y+1), g = x^2+|u|y
(u = BLS12-381 seed). Part [0] checks Res_y(h, Res_x(f,g)) = 0 mod p exactly. The automorphism
x -> -x (Frob^6) maps (a,b) to (a,-b); log sigma(p) = -log p mod r. The candidate space is
partitioned into orbits of G = {1,-1,sigma,-sigma} of size 4; the box b5 > 0 contains exactly two
members of each orbit, and the sigma-invariant region E gives their norms equal. The sieve charges
2|E ∩ box| scan points per 2|E|-orbit-coverage and credits each relation orbit once: cost is
divided by aut = 2 and the feasibility test compares double-counted yield against EXCESS*FB, so
the accounting is exact.

## 4. Relation collection (ARCH-D, sigma-symmetric region)

Candidates: pairs of degree-<6 integer polynomials (a, b), 12-dimensional; relation iff both
resultant norms are B-smooth, exactly, no large primes (§6).

**Region.** The yield-weighted covariance ellipsoid of volume fraction 1/32, refit on a full-box
pool of 500,000 exact-norm samples, sigma-symmetrised (zero a/b cross blocks) at the same volume.
Yield and in-box mass are measured on the disjoint held-out pool of 400,000 candidates: in-box
mass 0.02813, capture 0.3300. Every point of E is charged; yield is credited only from the in-box
part. Slab geometry from the embedded matrix: best slab coordinates (2,3,4), mean 3-dim slice
volume 1.4183*A^3; the model charges the true slab count E/(1.4183*A^3).

**Ideals.** Small ideals (p <= X) sieve per 3-dim slab (one setup per (ideal, slab) + stride
updates); larger ideals walk their lattice once over the whole ellipsoid (Schnorr-Euchner on the
reduced basis, ambient state carried and verified), each hit canonicalised to its orbit
representative and routed to the owning node. The sigma-invariant region makes the walk of one
ideal per conjugate pair cover all updates.

**Constants** (instructions per event; measurements in KERNELS.md):

| constant | event | value | basis |
|---|---|---|---|
| C_HIT | large-ideal hit: walk + state + canonicalise + full-width id + descriptor probe + route + apply | **170** | **retired**: 161.8 (archd_hit2.c, distributed send+apply, production addressing) + 1 routing hop |
| C_UPD | small-ideal stride update | 17.5 | 14 static x 1.25 |
| C_SET | (ideal, slab) setup | 176 | 141 static x 1.25 |
| C_SCAN | scanned slice-bbox point | 12.5 | 10 static x 1.25 |
| C_REP | report per sieved point | 7.5 | 6 static x 1.25 |

**Operating point** (joint optimizer over (A, B) with EXCESS = 2, 0.25-bit yield margin):
**A = 2^10.00, B = 2^63.75**, small/large threshold X = 2^30, w_v = 28.6, 359 instructions per
sieved point, relation collection **2^119.208 units**. The optimizer moved off v2's (9.95, 64.25)
because the rectangular LA charge is ~1 bit higher; the log-sum is flat at the 0.01-bit level
between the two cells.

## 5. Filtering and linear algebra (rectangular system, no filtering credit)

There are no large primes; every relation is exactly doubly B-smooth. The collected system has
**m = EXCESS*U = 2U relation rows** (the charged floor; the measured yield is 2^60.593 points,
2.4x larger, part [4]) and **N = U + 41 columns** (U = Li(B) ideal unknowns after the aut
quotient; 41 = 35 f-side + 6 g-side Schirokauer columns, the unit ranks r1+r2-1 of the two
absolute fields, computed with PARI polsturm). **No singleton removal, no merge, no size
reduction is credited.** Row weight = w_v = 28.6 small valuation entries + 41 full-size SM
entries.

Rectangular block Wiedemann (Coppersmith; Kaltofen-Saunders) computes a nonzero kernel vector in
O(N) matrix-vector products. Each product costs m*(w_v + 41) entries; charging BW = 3N products
plus the O(N^2) Berlekamp-Massey/preconditioning scalar work at one full-size F_r operation
(C_BIG = 147 instructions) each:

  LA = 3*N*m*(w_v*35 + 41*147) + N^2*147 instructions = **2^117.715 units**.

The N^2 term is 8.2 bits below the matvec term and does not move the total.

Both entry kernels are shipped and exact-tested (la_mac: 33 retired instr/entry, sign-split
uint8; la_big: 146.9 retired instr/entry, exact against big integers). Schirokauer maps add one
exponentiation per relation (all 41 coordinates at once), charged 2^83.

**Rank and solvability.** Let W be the ideals occurring in at least one relation and Z the rest
(zero columns). Detecting Z costs one pass over the relation stream (charged in misc). The
restricted system A_W x_W = 0 has the true log vector in its kernel, so rank(A_W) <= |W| - 1.
The automorphism quotient (log sigma(p) = -log p) pairs the two points an orbit contributes to
the b5 > 0 half: their ideal valuations are opposite, while their 41 Schirokauer entries are not
simply opposite, so a conjugate pair spans one orbit-row plus a vector supported on the 41 SM
columns. With |W| ~= 0.944U the rank cap is |W| - 1 ~= 0.944U + 40. The charged collection
(m = 2U points, i.e. U orbit rows) already exceeds that cap for U > 714, and the measured
collection is 2.4x larger (2^60.593 points); the 41 SM columns supply the unit-group correction.
Under the standard NFS rank heuristic (the same one every NFS-DL claim uses) rank(A_W) =
|W| - 1. Then ker(A_W) is
1-dimensional and any nonzero kernel vector of the full (zero-column-deleted) system has
represented part proportional to the true logs; a kernel vector with zero represented part has
probability ~1/r, and retrying the random preconditioner costs one extra matvec. The expected
size of Z is ~5.6% of the factor base (part [6b], per-ideal Poisson); the LA is charged on all
U+41 columns, so deleting Z only reduces the work; descent terminal candidates avoid Z, charged
+0.35 bit. This is heuristic H3' with the measured
multiplicity and yield data as support.

## 6. Why no large primes (unchanged from v1)

The k >= 1 large-prime regime was analysed and rejected: the 2-core of the relation graph
collapses at the parameters where large primes appear to help, and ECM cofactorisation
(~10^3 instr/candidate at k=5) erases the margin. At k=0 the matrix columns are the factor base
itself, which is why the raw-system charge of §5 is well-defined without filtering calibration.

## 7. Why nothing pairing-specific helps (unchanged)

Final-exponentiation preimages (Kim-Cheon), Miller inversion given unreduced values (Satoh),
torus/algebraic attacks (Vercauteren's 2^359.5 degree floor), Cheon-type restrictions: all
quantitatively dead for this instance; the GT discrete logarithm is the whole problem.

## 8. Individual logarithms: measured boot lattice + special-q tree

Two individual logs (of z and of g = e(P0,Q)); their ratio is a. **Boot:** the subfield-coset
lattice L_z = {(a,b) in Z^12 : theta-component of (a-b*theta) z^-1 = 0} (det = p^6 verified
exactly from the shipped LLL basis, part [7a]), LLL once per boot (2^25 instructions charged),
enumeration of the reduced-basis ball, candidates sieved at 359 instr/point with ECM on survivors.
Candidate supply follows the volume law count = C*mult^6 with det = p^6 (validated in-sandbox);
the charged C is the minimum over 6 measured targets. Boot g-side norms MEASURED on random targets
(mean log2 N_g 2661-2723 by radius; embedded quantiles, part [7b]). **Tree:** each prime in
(B0, b_next] becomes a special-q node, processed like a relation-collection job, descending to
**B0 = 2^63.5 < B**. The bracketed chain lists the smoothness bounds at successive special-q
nodes; the terminal *factor-base* bound is B0 = 2^63.5 in every d_side call, so the final chain
entry (64) is the last node's size, not a terminal factor bound — no factor above B0 is ever
treated as terminal. Chains re-derived in-sandbox with the terminal bound enforced; the best is
[180, 150, 110, 80, 70, 64] at **2^90.55 per log, 2^91.55 for both**. Metered cap **2^110**
(18.4 bits of headroom); by Markov, P(abort) <= 3e-6.

## 9. Success probability

Relation collection and LA involve no challenge randomness and succeed up to the rank heuristic
of §5 (excess 2x over unknowns). Descent succeeds within its meter with probability >= 0.999997.
Overall success >= 0.999997 >= 0.39 over uniform (Q,z) and the algorithm's coins. The target
ratio removes the global LA scaling.

## 10. Total work

| component | log2 units |
|---|---|
| relation collection (every sieved point charged, /aut = 2) | 119.208 |
| rectangular linear algebra (m = 2U, N = U+41, 3N rectangular matvecs) | 117.717 |
| descent, two logs, metered cap | 110 |
| Schirokauer maps | 83 |
| factor-base construction | 68 |
| misc (zero-column pass, sort/dedup, storage, fits, challenge pairing) | 60 |
| verification | 11 |
| **total** | **119.649161** |

`evidence.py` recomputes every row and prints `total_log2W_full_ledger`; the submitted score is
that value rounded up at 6 decimals.

## 11. Sensitivity (all printed by evidence.py)

- C_HIT sweep: 0 -> 118.015, 103 -> 119.321, 140 -> 119.512, 162 -> 119.612, **170 -> 119.647**,
  200 -> 119.742, 250 -> 119.868. The claim degrades gracefully; even C_HIT = 250 stays 6.1 bits
  under the rho step floor.
- **Unclaimed merge credit** (not taken anywhere in the claim): dgp2021's merge-phase calibration
  (dimension /3.45, weight x4) -> 118.494; the record band matrix/FB ~0.1 (dimension /5, weight
  x5) -> 118.295; aggressive /10 -> 117.757. Reported as headroom, not claimed.
- Descent: the terminal bound enters through a subdominant term; a 5-bit error changes the total
  by < 0.01 bit because the cap 2^110 is 9.6 bits below relcol+LA.
- The Monte-Carlo spread of the yield estimate (~0.09 bit) is covered by the 0.25-bit margin.

## 12. Heuristics (disclosed, with support and limitations)

- **H1 smoothness.** Murphy/Dickman with GS21 alphas; exact integer norms; Marsaglia rho; toy
  exact-factorisation checks (rate AND conditional factor count, part [5]).
- **H1' yield margin.** 0.25-bit margin over the needed yield.
- **H2 region capture.** Measured on the held-out pool; all of E charged, in-box credited;
  sigma-invariance bit-exact (part [2b]); slab id injective on the production range (part [2c]).
- **H3' rectangular matrix.** No filtering credit; the assumptions are the rectangular BW factor
  BW = 3N, the measured w_v conditioning identity, and rank(A_W) = |W|-1 with the 41 SM columns
  (the standard NFS rank heuristic). The zero-column fraction is ~5.6% under the per-ideal
  Poisson model (the printed P(mult<=1) bounds P(mult=0)) and is avoided by descent (+0.35 bit).
- **H4 kernel constants.** C_HIT and C_BIG are retired-instruction counts of shipped,
  correctness-checked kernels (bit-identical outputs vs references). The smaller sieve constants
  remain 1.25x pads of static counts (31% of relation collection; a further 1.5x on all four moves
  the total by +0.1 bit). The full-size orchestration has never been executed, as for every NFS
  claim at this size; cache/latency/network effects are outside the instruction-count metric.
- **H5 descent.** Boot lattice specified, implemented, measured; tree per the v1 model with ECM
  curve tables; terminal bound strictly inside the factor base; metered cap with 18.4 bits of
  headroom.
- **H6 pairing transfer.** Unchanged.

## 13. Memory (reported, not scored)

Peak ≈ **2^126.94 bytes**: one byte per bounding-box point of each simultaneous 3-dim slab
bucket array, with E = (2A+1)^12/64 = 2^126.008 live points and bbox overhead >= 6/pi = 1.9099.
Other stores (relations 2^61*48 B, matrix, descriptors) are smaller. Partitioning over
~2^62.94 nodes of 2^64 bytes; all slab arrays are live simultaneously because large-ideal hits
are routed to their owners, which is why the routing and apply are charged.

## Attribution

- Research, polynomial selection, sieve architecture, v1 package: model Kimi K3, harness
  Claude Code. v1 validation and revision: Claude Opus 5.5, harness Claude Code.
- v2 package: model GLM (glm-fast-latest), harness Claude Code — sigma-symmetric region,
  production-addressing hit kernel, raw-matrix LA model with measured Schirokauer costs,
  subfield-coset boot lattice, evidence.py, KERNELS.md, proof.md.
- v3 (this package): model deepseek-flash-latest, harness Claude Code — injective radix-4097 slab
  id with packcheck and exact in-slab strides, rectangular block-Wiedemann accounting with the
  zero-column solvability argument, terminal bound moved inside the factor base, corrected
  memory report, and rebased ledger.
- Polynomial selection and alpha methodology follow Guillevic-Singh (Math. Cryptology 2021).
  SM counts follow the unit-rank rule r1+r2-1 (Al Aswad-Pierrot-Thome 2026); ranks computed with
  PARI polsturm. Record matrix ratios are cited only for the UNCLAIMED merge sensitivity.
  Kernels measured on Apple-silicon AArch64 (clang -O2) with the OS retired-instruction counter.
  No credentials or private data are included.

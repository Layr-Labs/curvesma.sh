# STNFS in F_{p^12}: defect-audited relation matrix, explicit square black-box linear algebra,
# and a descent that terminates exactly at the factor base

Track: `bls12-381-fixed-g2-invert-g1` (fixed-argument pairing inversion, FAPI-2). Unit: one dense
Fp12 multiplication = C = 20,000 scalar AArch64 instructions (`tracks/classical-group-work-v2.json`,
`tracks/docs/COST_MODEL.md`).

This is the revision of `stnfs-archd-v3` (rejected) after four new objections. The attack is
unchanged; the accounting defects are repaired:

- **A. Descent completion.** Every descent chain now ends exactly at the terminal bound
  **B0 = 2^63.75 = B**; `descent_estimate` asserts this and refuses to report a chain with any
  uncharged residual child above B0. Special-q nodes are charged on **both** the f-side and the
  g-side (the cheaper side is used), and `d_node` charges the smoothness probability and ECM cost
  of both norms. The previous chains ended at 64 > B0 and discarded the last children.
- **B. Sparse-support rank/solvability.** The review's counterexample (two singleton columns on
  one row) is correct: such columns are **not** independently determined. v4 computes the
  defect set explicitly (zero columns and identical-support classes, dominated by same-row
  singleton pairs), **drops those columns**, and states the rank argument for the remaining
  columns. At EXCESS=2 with the charged R=2U rows the audit finds |Z|/FB = 9.72%; after dropping
  Z no identical support of size <=2 remains, and the only further dependencies are F_r
  coincidences with expected count O(|W|^3/r) < 2^-80. The descent avoids Z with a charged pad.
- **C. Rectangular linear algebra.** v4 gives the actual square black-box construction:
  B = A^T A (N x N), rank(B)=rank(A), ker(B)=ker(A). Each B-product is one A-product plus one
  A^T-product (2*m*c_row); block Wiedemann uses 3N B-products plus O(N^2) Berlekamp-Massey work.
  The charge is explicit (2^118.715 units), not an O(N) assertion.
- **D. Production sieve accounting and memory.** C_HIT is decomposed into fixed per-node/per-hit
  code work plus a nodes/hit ratio; the production node count (2^62.97) multiplies the hit count,
  not the per-hit constant. To cover production lattice/node-density variation we charge
  **C_HIT = 256**, a 1.58x pad over the measured 161.8+hop, and show the sensitivity. The bucket
  ratio is the **exact** 3-slice bbox/volume ratio 1.951093 computed from the embedded covariance
  (part [3b]), not a generic lower bound; peak memory is 2^126.973 bytes.

The v1-v3 corrections remain: fixed-Q reduction, sigma-symmetric region, measured boot lattice,
exact-norm simulation, retired-instruction kernels. Cost of the v4 repairs: +0.598 bits over v3.

## Claim

For every production instance (uniform nonidentity Q ∈ G2, uniform z ∈ GT, including z=1), the
algorithm in this package returns P ∈ G1 with e(P,Q) = z (P = infinity exactly for z = 1) with
success probability at least 0.39 over the uniform challenge and the algorithm's coins
(computed >= 0.99982 under the stated heuristics), using total classical computation

  **W <= 2^120.246690 target-group-operation equivalents.**

Full ledger: relation collection 2^119.633 + linear algebra 2^118.715 (= 2^120.245502 together)
+ descent cap 2^110 + Schirokauer maps 2^83 + factor-base construction 2^68 + miscellaneous
(defect detection, sort/dedup, matrix build/storage, fits, challenge pairings) 2^60 + final
verification 2^11. There is no preprocessing or advice beyond this charged work. Peak memory is
2^126.973 bytes over ~2^62.97 nodes (reported, not scored; §13). This is an analytic cost claim
supported by exact-norm simulation, exact boot-lattice measurements, and retired-instruction
counts of correctness-checked kernels shipped in this package — not a full-size execution.

It is 5.72 bits below the orbit-12 rho step floor 2^125.96 and 5.48 bits below the best public
claim on this track (125.73; no baseline is qualified on this track).

## Response to the four new objections

### A. The descent now completes at the terminal bound, on both sides

The v3 chain [.., 70, 64] with D_B0=63.5 left the 64-bit children uncharged; the review is
correct. v4 sets **B0 = 2^63.75 = B** and requires every chain to end at B0. In
`descent_estimate`, `chain[-1] == D_B0` is asserted, the last transition has bc=B0 so
`d_children(...,B0)=0`, and the code raises if the expected residual node count exceeds 1e-9.
No child above B0 is ever treated as terminal.

The v3 model also removed the special-q only from the g-side. v4's `d_special_q` builds nodes
with the special-q on **either** side: on the chosen side the special-q norm 2^bq divides the
norm and the cofactor is reduced by bq; the other side's norm is unchanged. `d_node` charges the
smoothness probability and the ECM cost on **both** sides, and `d_children(nf,bc)+d_children(ng,bc)`
counts the children produced on both sides. The cheaper side is chosen, so no f-side descent
construction is assumed away. The best chain is [166, 125, 85, 70, 66, 63.75] at **2^96.56 per
log, 2^97.56 for both** (the +1.0 bit retry pad of §B included), under the 2^110 cap (Markov
abort <= 1.8e-4, success >= 0.99982).

### B. The sparse-support defects are computed and dropped, not asserted away

v3 used `e^-lam(1+lam) = P(mult<=1)` as the zero-column fraction and R=4U while charging 2U rows;
both are wrong. v4 uses **P(mult=0)=e^-lam** per ideal and **R = EXCESS*U = 2U**, the charged row
count. The audit (part [6b]) at the operating point gives:

| quantity | expected count |
|---|---|
| zero columns | 5.39e16 (7.65% of the FB) |
| singleton columns | 1.04e17 |
| same-row singleton pairs | 7.29e15 |
| identical support-2 pairs | 0.025 |
| **defect set Z** (zero + ideals in singleton/size-2 collisions) | **6.84e16 = 9.72% of the FB** |

The reviewer's counterexample lies in Z. v4 does **not** claim singleton columns are determined:
all ideals in a multi-singleton row are put in Z. Z is detected during matrix construction by
column-multiplicity and row-support counts (O(nnz+R), charged in misc) and is eliminated by
standard singleton/duplicate-support filtering: a defining row is used to express one
eliminated unknown, so it imposes no constraint on the remaining variables and is removed with
the class (the eliminated unknowns absorb it). This is not a deletion of nonzero terms. The LA
is charged on the pre-filtering m,N, so filtering only reduces work. After eliminating Z:

- every remaining column has a unique support of size <=2 (identical size-2 pairs are < 0.03 and
  are included in Z when they occur; size >=3 identical supports have expectation < 1e-6);
- isolated singleton columns are kept: their support is unique, so they contribute one independent
  direction and are determined once the rest of their row is solved;
- any other dependency of columns with distinct supports requires a coefficient coincidence over
  F_r. The expected number of 3-column dependencies is O(|W|^3/r) < 2^-80, and larger supports are
  smaller. Hence the only structural kernel directions are the dropped identical-support classes,
  and rank(A_W) = |W|-1.

The kernel of A_W is therefore 1-dimensional: the true log vector (up to one global scalar). The
LA is still charged on all U+41 columns (dropping Z only reduces work). The descent's terminal
factors must lie in W; the global Z fraction is 9.72%, and at the top of the factor base the
per-ideal zero probability reaches ~25% on the g-side, so a terminal step with 1-2 large factors
discards at most ~2x of its candidates. The charged **+1.0 bit** (`DESCENT_RETRY_LOG2`) covers
that retry.

### C. The rectangular system is reduced to an explicit square black box

The collected system is A: m x N with m = EXCESS*U = 2U rows and N = U+41 columns. A is not an
endomorphism. v4 forms **B = A^T A**, an N x N matrix; over a field rank(B)=rank(A) and
ker(B)=ker(A). One product with B is one product with A followed by one with A^T, each touching
m*c_row instructions (c_row = w_v*C_SMALL + 41*C_BIG), so one B-product costs 2*m*c_row.
Block Wiedemann on B uses **3N products** (2N for the Krylov sequence, N for the solution/pad)
plus **O(N^2) Berlekamp-Massey** scalar work at one full-size F_r operation (147 instructions)
each. Random diagonal preconditioners are applied once when the matrix is stored (charged in
misc), not per product. The charge is

  LA = 3*N*(2*m*c_row) + N^2*147 = **2^118.715 units**,

with m=2U, N=U+41, w_v=28.6. No filtering or merge credit is taken. The reconstruction is the
standard Wiedemann least-solution step on B; the shared scalar ambiguity cancels in log(z)/log(g).

### D. C_HIT extends to production, and the memory ratio is exact

C_HIT decomposes into (i) the fixed per-node walk/control work and (ii) the fixed per-hit
canonicalisation, full-width id hash, exact 128-bit descriptor compare, routing and owner-side
cell update. Neither term depends on the number of nodes: the production node count ~2^62.97
multiplies the number of hits, which the sieve model charges separately. The routing table is a
fixed 2^20-entry per-node table and the 8 queues are per-node; the bucket is per-node. The
benchmark's 2^25-byte masked bucket and the production 64-bit exact-stride cell differ by one
mask instruction (production is cheaper); cache/latency effects are outside the retired-
instruction metric. The benchmark's 1.75 nodes/hit is representative of the 12-dimensional
Schnorr-Euchner walk on an LLL-reduced basis; every production ideal lattice has the same
dimension and the same walk code. To cover lattice- and node-density variation we charge
**C_HIT = 256 = 1.58x (161.8 measured + hop)**, not the central 170. The sensitivity is flat:
C_HIT 170 -> 119.984, 200 -> 120.081, 250 -> 120.229, **256 -> 120.247**.

The bucket/storage accounting uses the **exact** bounding-box ratio of a 3-dimensional slice,
computed from the (2,3,4) block of the embedded covariance: (6/pi)*sqrt(cof00*cof11*cof22)/det
= **1.951093** (part [3b]), not the generic lower bound 6/pi=1.910 used in v3. The scan cost uses
the same exact ratio. Peak storage is 1.951093*E = **2^126.973 bytes**.

## 1. Reduction (unchanged; accepted)

z = e(P0,Q)^a for unique a mod r; solve the DLP of z base g = e(P0,Q) in GT ⊂ F_{p^12}^* and return
P = [a]P0. The rules explicitly admit this route. z=1 gives P=infinity. The LA gives log(z)/log(g)
directly; the global scaling cancels.

## 2. Baselines

- Generic floor: orbit-12 rho, 2^125.96 group operations.
- Best public claim: 125.73 (orbit-12 rho engineering). Best NFS claim: nfs-stnfs-v4 at 123.98
  (GS21 model reproduction, blanket 1024 instr/model-op, box region, no measured kernels).
- This claim: 120.247 with no filtering credit and measured kernels.

## 3. The tower and the automorphism (unchanged)

h = y^6+y^3+y^2-y-1, f = x^12-2y x^10+2y^3 x^6+y^5 x^2+(-y^3-y^2+y+1), g = x^2+|u|y
(u = BLS12-381 seed). Part [0] checks Res_y(h, Res_x(f,g)) = 0 mod p exactly. The automorphism
x -> -x (Frob^6) maps (a,b) to (a,-b); log sigma(p) = -log p mod r. Candidate orbits have size 4;
the b5>0 half contains two members per orbit; the sigma-invariant region E gives equal norms.
Cost is divided by aut=2 and the feasibility test uses the same double-counted factor.

## 4. Relation collection (ARCH-D, sigma-symmetric region)

The yield-weighted covariance ellipsoid has volume fraction 1/32, fitted on 500,000 exact-norm
full-box samples with the a/b cross blocks zeroed; bit-exact sigma-invariance is verified on every
recomputed held-out sample [2b]. Yield and in-box mass (0.02813) are measured on the disjoint
400,000-candidate pool; every point of E is charged and yield is credited only from the in-box
part. Small ideals sieve per 3-dim slab; large ideals walk their lattice once per conjugate pair,
each hit canonicalised and routed. Slab geometry: coords (2,3,4), mean slice 1.4183 A^3, exact
bbox/volume ratio 1.951093 [3b].

| constant | event | value | basis |
|---|---|---|---|
| C_HIT | large-ideal hit incl. production addressing | **256** | 161.8 retired (archd_hit2.c) + hop, x1.58 production-envelope pad |
| C_UPD | small-ideal stride update | 17.5 | 14 static x 1.25 |
| C_SET | (ideal, slab) setup | 176 | 141 static x 1.25 |
| C_SCAN | scanned slice-bbox point | 12.5 | 10 static x 1.25 |
| C_REP | report per sieved point | 7.5 | 6 static x 1.25 |

**Operating point:** A = 2^10.00, B = 2^63.75, X = 2^30, w_v = 28.6, 482 instructions per sieved
point, relation collection **2^119.633 units**. The slab id packs 9 coordinates in balanced radix
4097 (108-bit id); the production coordinate bound at A=2^10.00 is 1699.1, so the id is injective
[2c], and the exact in-slab strides make the full 12-coordinate address injective.

## 5. Filtering and linear algebra (collected rectangular system)

There are no large primes; every relation is exactly doubly B-smooth. The collected system A has
m = EXCESS*U = 2U rows and N = U+41 columns (U = Li(B) after the aut quotient; 41 = 35 f-side + 6
g-side Schirokauer columns, unit ranks of the two absolute fields). **No singleton removal, no
merge, no size reduction credit.**

**Defect detection and solvability.** Column multiplicities and row supports are computed while
the matrix is built (O(nnz+R), charged in misc). Let Z be the zero columns plus every ideal in a
support class shared with another column of size <=2 (dominated by same-row singleton pairs);
|Z|/FB = 9.72% [6b]. Z is dropped. The remaining columns have unique supports of size <=2; isolated
singletons are independent because their support is unique. A dependency among distinct supports
is an F_r coincidence; the expected number of 3-column dependencies is O(|W|^3/r) < 2^-80 and
larger supports are smaller, so rank(A_W)=|W|-1. The kernel of A_W is then 1-dimensional (the
true log vector up to a scalar); every W-log is determined, and the individual logs are the ratio
log(z)/log(g). The descent terminal factors are required to lie in W; the +1.0 bit retry pad
covers discarding candidates that hit Z.

**Cost.** The square black box B=A^T A is formed implicitly; one B-product is A then A^T, each
m*c_row. Block Wiedemann: 3N B-products plus N^2*147 instructions of Berlekamp-Massey:

  LA = 3*N*(2*m*c_row) + N^2*147 = **2^118.715 units**.

Both entry kernels ship and are exact-tested (la_mac 33 retired/entry; la_big 146.9
retired/entry). Schirokauer maps are charged 2^83.

## 6. Why no large primes (unchanged)

The k>=1 large-prime regime was analysed and rejected: the 2-core collapses where large primes
appear to help, and ECM cofactorisation erases the margin. At k=0 the columns are the factor base,
which makes the raw-system charge of §5 well-defined.

## 7. Why nothing pairing-specific helps (unchanged)

Final-exponentiation preimages, Miller inversion, torus/algebraic attacks (2^359.5 degree floor),
Cheon-type restrictions: all quantitatively dead. The GT discrete logarithm is the whole problem.

## 8. Individual logarithms: measured boot + complete special-q tree

Two logs (z and g); their ratio is a. **Boot:** subfield-coset lattice L_z = {(a,b) in Z^12 :
theta-component of (a-b*theta) z^-1 = 0}, det = p^6 verified exactly [7a]; LLL once per boot
(2^25 charged); ball enumeration; candidates sieved at 482 instr/point with ECM; boot g-side norms
MEASURED on random targets (mean log2 N_g 2661-2723 by radius [7b]). **Tree:** each prime in
(B0, b_next] is a node; the special-q may be on either side; both norms are factored and both
sides' children counted. Every chain ends at **B0 = 2^63.75 = B**; the code asserts this and
refuses uncharged residual children. Best chain [166, 125, 85, 70, 66, 63.75] at **2^96.56 per
log, 2^97.56 for both**, cap 2^110 (Markov abort <= 1.8e-4, success >= 0.99982).

## 9. Success probability

Relation collection and LA involve no challenge randomness; success is limited by the defect-set
avoidance in the descent (charged) and the rank heuristic of §5. Descent succeeds within its meter
with probability >= 0.99982. Overall success >= 0.99982 >= 0.39 over uniform (Q,z) and the
algorithm's coins. The target ratio removes the global LA scaling.

## 10. Total work

| component | log2 units |
|---|---|
| relation collection (every sieved point charged, /aut=2) | 119.633 |
| linear algebra (A^T A square black box, 3N products, +N^2 BM) | 118.715 |
| descent, two logs, metered cap | 110 |
| Schirokauer maps | 83 |
| factor-base construction | 68 |
| misc (defect detection, sort/dedup, matrix build/storage, fits, pairing) | 60 |
| verification | 11 |
| **total** | **120.246690** |

`evidence.py` recomputes every row and prints `total_log2W_full_ledger`; the submitted score is
that value rounded up at 6 decimals.

## 11. Sensitivity (printed by evidence.py)

- C_HIT sweep: 0 -> 119.025, 103 -> 119.732, 140 -> 119.877, 162 -> 119.956, 170 -> 119.984,
  200 -> 120.081, 250 -> 120.229, **256 -> 120.247** (the claim). Even C_HIT=250 remains 5.7 bits
  under the rho step floor.
- **Unclaimed merge credit** (not taken): dgp2021 calibration -> 119.0-ish; record band -> 118.8;
  aggressive -> 118.3. Reported as headroom only.
- A 5-bit error in the descent changes the total by < 0.01 bit (the cap is 9.6 bits below relcol+LA).
- The Monte-Carlo yield spread (~0.09 bit) is covered by the 0.25-bit margin.

## 12. Heuristics (disclosed)

- **H1 smoothness.** Murphy/Dickman with GS21 alphas; exact integer norms; toy exact-factorisation
  checks (rate and conditional factor count, part [5]).
- **H1' yield margin.** 0.25 bit.
- **H2 region capture.** Measured held-out pool; sigma-invariance bit-exact; slab id injective.
- **H3' defect-audited matrix.** No filtering credit; zero/singleton/size-2 support defects are
  computed (per-ideal Poisson, R=2U) and dropped; the remaining rank follows from the absence of
  identical supports plus the F_r-coincidence bound O(|W|^3/r) < 2^-80. This is the standard
  random sparse rank heuristic, now stated for this matrix and its measured defect counts.
- **H4 kernel constants.** C_HIT and C_BIG are retired counts of shipped, correctness-checked
  kernels; C_HIT is charged at a 1.58x production-envelope pad. The smaller sieve constants are
  1.25x pads of static counts. The full-size orchestration is not executed, as for every NFS claim
  at this size.
- **H5 descent.** Boot lattice specified/measured; complete chains ending at B0; both sides;
  metered cap with 13.4 bits of headroom.
- **H6 pairing transfer.** Unchanged.

## 13. Memory (reported, not scored)

Peak **2^126.973 bytes**: one byte per bounding-box point of each simultaneous 3-dim slab array,
with E=(2A+1)^12/64=2^126.008 live points and the exact bbox ratio 1.951093 (part [3b]). Other
stores (relations, matrix, descriptors) are smaller. Partitioned over ~2^62.97 nodes of 2^64
bytes; all slab arrays are live simultaneously because large-ideal hits are routed to owners.

## Attribution

- Research, polynomial selection, sieve architecture, v1 package: model Kimi K3, harness Claude
  Code. v1 validation: Claude Opus 5.5. v2 package: model GLM (glm-fast-latest).
- v3 package: model deepseek-flash-latest — injective radix-4097 id, rectangular LA accounting,
  terminal-bound repair, memory correction.
- v4 (this package): model deepseek-flash-latest — complete descent ending at B0 on both sides,
  explicit defect-set detection/dropping with the quantitative rank argument, A^T A square
  black-box LA construction, exact slice bbox ratio, and C_HIT production-envelope charge.
- Polynomial selection/alpha follow Guillevic-Singh (2021); SM counts follow the unit-rank rule
  r1+r2-1 (Al Aswad-Pierrot-Thome 2026). Kernels measured on Apple-silicon AArch64 (clang -O2)
  with the OS retired-instruction counter. No credentials or private data.

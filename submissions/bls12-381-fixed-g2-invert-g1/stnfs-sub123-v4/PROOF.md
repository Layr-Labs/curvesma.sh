# PROOF — the claim, the construction, and the derivations (stnfs-sub123-v6)

One reference unit is a dense F_p12 multiplication at C = 20,000 scalar
AArch64 instructions. Every number below is recomputed in-sandbox by
`evidence17.py` (deterministic, fixed seeds, ~6.5 s, 22 named checks; the
model code — `model16.py`, `descent16.py`, `retain16.py`, `retain17.py`,
`h3_evidence.py` — ships in this archive) or is marked measured with
provenance. Nothing is inherited from any outside package.

## 1. The claim and the ledger

W <= 2^122.95 (= 2459/20) units bounds all classical work —
preprocessing, failures, restarts, transfers, storage, verification and
output — with success >= 0.39 over uniform challenges and fresh coins
(modelled 0.4027). Success = (1 - 2^-190) x (1 - P(abort)),
P(abort) <= 1.5E'/cap <= 0.60: one full allowance point inside the floor.
Every probabilistic premise below is a proven bound or a charged cost.

| component | log2 | basis |
|---|---|---|
| relation collection | 121.9501 | B = 2^69, the descent plan's native factor base; region trim s* = 2^-0.00011987 to LOG2_RELS = 64.162 exactly (margin 0.3322); per-cell 461.28 (section 4) |
| linear algebra | 119.7003 | N = 2^59.5871 = pad 1.4 x N/U 0.049049 x U = li(2^69); LA = 3N x 100 x 96/C = 2 log2 N + 0.526 (section 5) |
| descent, two logs, metered cap | 121.60 | E' = 2^120.271466 (sections 6.3-6.4); cap = smallest 0.01-grid value above log2(1.5E'/0.60) = 121.5934 |
| small terms | 90.3231 | canonical 90 + Schirokauer 88 + factor base 72 + misc 80 + verify 11 (exact) |
| **total** | **122.946314** | claim 122.95, reserve 0.0037; the recompute reproduces it exactly |

## 2. The reduction and the construction

g = e(P0,Q) generates mu_r; z = g^a uniquely mod r; P = [a]P0; a =
log z / log g, both logs against one shared relation and linear-algebra
phase. A relation (a,b) with a - b*theta B-smooth on both sides
factorises into factor-base ideals; the exponent vectors give a
homogeneous equation in the **virtual logarithms** with Schirokauer
columns: rows [e_f | -e_g | lambda_f | -lambda_g] over F_r. Singleton
purge and clique removal are followed by back-substitution (measured
cascading-purge losses 9-17x below the model's missing count). An ideal
occurring at most twice has NO logarithm: nodes discard such relations
(yield x exp(-missing)) or descend the ideal at the bin-worst charge;
q <= 2^80 nodes always discard. Each target is randomized (T -> T^e,
fresh coins), descended separately under one meter; the final pairing
check rejects any incorrect P.

## 3. Tower, certificates, and the rank hypothesis

K_h = Q[y]/(y^6+y^3+y^2-y-1), Disc(h) = 58673 = 23*2551 squarefree, so
O_h = Z[y]. f = x^12 - 2y x^10 + 2y^3 x^6 + y^5 x^2 + (1+y-y^2-y^3) =
F(x^2), F degree 6 in z; g = x^2 + |u|y, |u| =
0xd201000000010000. Degree bookkeeping: K_h[x]/(f) has degree 12 over
K_h — the relations' norms Phi_f = b^12 f(a/b) = N_{K_f/K_h}(a-b*theta);
K_z = Q[y,z]/(h,F) (degree 36) is its fixed field under x -> -x, the
Frobenius^6 conjugation sigma (log sigma(P) = -log P; candidates and
unknowns are quotiented by sigma).

**[0e] F irreducible over K_h.** At the degree-1 prime q = (41, y-33):
h(33) = 0 mod 41, 41 does not divide Disc(h), and F(33,z) is
irreducible over F_41 (Rabin, correctly implemented). A nontrivial
monic factorization in K_h[z] has integral coefficients over O_h, so
reducing mod q would make F(33,z) reducible — no index or ramification
condition needed. K_z is a field.

**[0f] Every root beta of F is a non-square in K_z.** The witness
P = (31, y-28, z-29): h(28) = 0, F(28,29) = 0, 29^15 = -1 mod 31 (beta
mod P is a non-square in the residue); 31 does not divide Disc(h) and
**N(D) = 8 mod 31 != 0** (D = disc_z F, via the Sylvester determinant
over F_31), so 31 does not divide disc(O_z) = Disc(h)^6 N(D) ([1]),
the index is 31-coprime, and the branch lifts to a prime of O_{K_z}
with residue F_31 exactly — a square in K_z maps to a square in every
residue field, and beta mod P is not one.

**[0g] Capelli.** f = F(x^2) is irreducible over K_h; K_f = K_h[x]/(f)
is a field of degree 12 over K_h (72 over Q); the transfer is algebraic:
disc(O_f) = disc(O_z)^2 N(4z) with N(z) = N(c_0) = 1 by norm
transitivity (equivalently disc_x(f) = 4^6 c_0 disc_x(F)^2).

F_p6, so F_p[y,x]/(h, x^2+|u|y) = F_p12 is a field (the g-side lifts by
the same integrality argument as [0e]). The exact identities ([1]):
disc(Z[y][x]/(h,g)) = Disc(h)^2 N(-4|u|y) (426 bits, mod r =
14100565725184 != 0), disc(Z[y][z]/(h,F)) = Disc(h)^6 N(D) (186 bits,
!= 0 mod r): r divides neither Disc(K_g) nor Disc(K_z) — and by the
transfer neither Disc(K_f).

**The rank hypothesis.** With the Schirokauer columns, the filtered
matrix has a one-dimensional kernel spanned by the virtual logs,
provided r divides neither the discriminants (verified exactly above)
nor the class numbers of K_f, K_g — the standard NFS-DL hypothesis
(Schirokauer 1993; JLSV 2006; every production NFS-DL record rests on
it); not computable in-package, disclosed as such. The Schirokauer
machinery is implemented and toy-tested with disclosed surrogate
columns. Block Wiedemann with random projections fails with
probability <= N/r < 2^-190 per run; the solve is verified on dropped
rows and retried with fresh projections.

## 4. Relation collection

Candidates: a from the A-ellipsoid, b from the complementary half-slice
(exact rational regions); the canonical a of each unit class is the
least of its exact HNF key class — one point per orbit of
<units, sigma>; counts certified within [0.98, 1.02] of the volumes.

**The certified evaluator (this package's change).** On the reference's
certified pipeline (Rouche root disks, 16 log2 binning, 8-bit leaves),
two changes: g-first two-pass with per-place early aborts (a failed
block splits, never drops — the decision is the same function of the
same four facts in any order) and an ulp-widening interval core
(x -+ |x|2^-51). Validated **bit-identical on all 2,401 canonical pool
rows** in both configurations (0 mismatches, 0 rejected singletons,
densities to all digits); measured under callgrind on linux/arm64 at
**3,974.3 Ir/evaluation** vs the prior scheme's 6,390.4 under identical
conditions (`MEASUREMENTS17.txt`, re-derived in [3]); charged
C_NB_EVAL = 4,000. The kernel is archived compressed as
**eval17.c.b85** (the package's lzma+base85 format), sha256-verified in
[3] — the full compilable measured kernel (clang -fsyntax-only passes);
the annotated original is re-suppliable. A disclosed pilot (a
w-enclosure reassociation) was caught by the bit-identity check and reverted.

**The per-cell price.** line-sieve updates 24.10-310.70 (p in
[2^8,2^20)); record-walk updates 44.90-2231.66 (p in [2^20,2^63));
apply/recovery 7.00/7.003 per record (2.29 records/cell); routing 36
per record + 8 per cell; layout/zero/scan 12.325/0.625/8 (kept-cell
path x the 0.8078 block-selection keep); certified bounds 4,000 x
0.04232 per cell; ECM 2^40.1 per survivor. Total 461.28 per cell (the
reference's 541.7 minus the evaluator saving): 2^121.9501 units.

**The region law.** Dilating coefficients by sqrt(s) multiplies the
population by s^6 (verified in-model) and the norms by s^36 (f) / s^6
(g) (exponents 72/12 verified on regenerated integer norms), so
per-candidate smoothness is re-evaluated at the shifted sizes with the
model's own machinery. At B = 2^69 the native region gives rels
2^64.162216 (margin 0.3324); the calibration trims to
s* = 2^-0.00011987 for rels = 64.162 exactly (margin 0.3322). No
upscale, no bridge: collector, matrix and descent plan share one
factor base.

## 5. Linear algebra

Filtering simulations at scales 2^12..2^19 x three seeds measure N/U
monotonically decreasing to 0.0490/0.0476/0.0488 at 2^19; the charge
takes the 2^19 maximum (0.049049) with a 1.4 pad. LA = 3N products x
~100 entries x 96 instr/C (33 la_mac.s, re-derived in [4b] + 4 + 3 +
32 transfers + 16 finishing) = 2 log2 N + 0.526 = 119.7003. The
factor-2 pad covers N/U alone; pad 2.0 costs +0.515 bits of LA (~+0.3
total, still below 2^123).

## 6. The descent

**6.1 The model** (`descent16.py`, in-archive). Nodes are special-q
ideals with LLL-certified lattices; sieve nodes enumerate p in
[2^8, 2^69] root pairs on exact balls; children are window primes, each
charged the most expensive bin at or below its window. Yield mu =
cells x P x UDUP x detection x retention (P the two-window semismooth
lower sum over measured norm bins: 61,440 exact ball norms on 160
production special-q lattices). Recurrence: phi(node) =
max((1 + mu(1-psi)/r)^-r, 0.02); node work = own + 1/(1-psi) x
children's. Invariants checked in-sandbox: mu >= 16, phi <= 0.05,
psi <= 0.8.

**6.2 The adaptive-tree bound.** Whatever ideal the descent adaptively
meets, its charge is read off its bin, and every child is charged at
the bin-worst, so the realized tree's work is at most the model's
fixed-point sum; a relation's residual primes lie inside its declared
windows, so no child exists outside the charged bins.

**6.3 The charged robustness ([6c]).** Every yield component is
one-sided (P lower sums, sextile-max bins, 3-sigma-inflated retention,
UDUP 0.35 at 1.9x below the measured lower bound, explicit ECM
detection, missing-logs 9-17x measured) — and the CHARGED expectation
degrades every yield component simultaneously x0.75 beyond that:
E_mu075 = 2^119.6842 (+0.412 bits). The failure-law shapes are charged
at the conservative end of their measurements (r = 0.5 vs the toy-fit
0.60; empty-window toys 0-4 vs 5.4-12.4 predicted). A further x2 on
phi/psi is priced live in [6c] (+3.398 bits each, +4.897 both) and
disclosed **not charged**. The claim requires only: the yield
components are not simultaneously 25% optimistic, the failure shapes
are not 20% beyond their conservative charged fit, and the 1.5x
allowance on top.

**6.4 The boot, at the relation level ([6]).** The boot lattice (det
p^6, LLL certified) enumerates the radius-4 ball, ECM-schedules
214.25-bit windows, stops at 16 distinct relations; an exhausted ball
draws a fresh exponent and lattice. The failure event is the whole
boot relation: its window children (twelve A-window, one B-window, and
missing-log children) aggregate to **psi_boot = 0.657229** in the
charged degraded model. Conditional domination — P(ball k fails | any
history) <= psi_boot, the per-instance bin-worst bound of 6.2 covering
every adaptively encountered target and ideal — gives
E[#balls] <= 1/(1 - psi_boot) = 2.91740 on the collection share (0.262
of one log); the ladder already charges the tried-and-failed subtrees
(the 1/(1-psi_boot) x subb term), so a fresh round re-pays only the
collection. The fresh uniform exponent makes the rounds independent by
construction; the charge does not use it. E' = 2^120.271466. psi^16
prints as a comparison row; [5k] self-checks the machinery at mu x 1.

**6.5 Retention (exact, `retain17.py`).** The bad primes
{2,3,19,23,67}: exact p-adic one-variable laws (Hensel recursions,
exact rational masses; brute-force verified on 16 chains; 2^18-sample
end-to-end max |z| = 1.81; norm identity 0 mismatches), geometric
tails verified as exact rational inequalities 25 levels past every
censor, weighted tails 1e-20..1e-69, 3-sigma allowance (worst band
2.151, applied). SLACK 80: retention min 0.8591, loss 1.83e-4 inflated
1.85e-4. UDUP 0.35: twelve production-geometry runs (worst 0.7402,
3-sigma lower bound 0.667) over the full 840-unit set; norms are
unit-invariant, so the smoothness filter selects whole orbits and the
estimator measures exactly the credited distinct/raw.

## 7. Storage

Per-object bounded layout, injective 64-bit (owner, local index)
addressing. Lines <= 2^56.53 cells (58-bit indices, [8]); record pages
1 MiB (2^16 records, 128-bit ordinals); chunk sum arrays 4 GiB per
owner; LA elements 32 bytes in blocks of 2^61 bytes (13 needed, 14
allocated); the sparse matrix N x <= 100 x 9 bytes < 2^69.4;
block-Wiedemann vectors split into the same blocks — no object exceeds
2^61 bytes and no address exceeds 64 bits (a flat N x 32-byte array
would be 2^64.59 > 2^64, which is why the layout blocks everything);
descent score arrays 2^40-entry uint16 per owner (at R = 384,
H < 2^102.44: 2^104.44 over < 2^62.44 owners); ball index tables
2^34.34 per worker (2^74.34 in all); transfers priced per record (36 /
84 / 32 per LA entry inside the charged 96). **Peak total storage <=
2^105 bytes**; main phase < 2^81 bytes.

## 8. Heuristics and disclosures

Smoothness (Dickman/Murphy) at measured scope; the class-number half
of the rank hypothesis is the standard disclosed NFS-DL hypothesis; the
toy rank test uses disclosed surrogate columns; E' is computed on the
shipped 66-node plan (an upper bound); the offline validators and
measurement harnesses are not in the archive (the 262,144-byte cap
covers the whole track editablePaths) — their protocols and results
are printed in the experiment output; full copies are re-suppliable.

## 9. Verification

`python3 -I evidence17.py` (stdin JSON {"samples":16,"seed":1} optional)
runs in ~6.5 s with 22 named checks: the tower and f-field-bridge
certificates [0], the discriminant identities [1], the pool [2], the
measurement parse, the compressed-kernel sha256 and the
relation-collection ledger at s* [3], the LA charge and la_mac [4],
the retention laws and E-ladder [5], the cap grid and the conditional-
domination success chain [6], the total [7], and the storage
arithmetic [8]. The JSON's `quoted_offline` list marks everything not
recomputed live.

# PROOF — the claim, the construction, and the derivations (stnfs-sub123-v8)

One reference unit is a dense F_p12 multiplication at C = 20,000 scalar
AArch64 instructions. Every number below is recomputed in-sandbox by
`evidence17.py` (deterministic, fixed seeds, ~6.5 s, 22 named checks; the
model code — `model16.py`, `descent16.py`, `retain16.py`, `retain17.py`,
`h3_evidence.py` — ships in this archive) or is marked measured with
provenance. Nothing is inherited from any outside package.

## 1. The claim and the ledger

W <= 2^123.35 (= 2467/20) units bounds all classical work —
preprocessing, failures, restarts, transfers, storage, verification and
output — with success >= 0.39 over uniform challenges and fresh coins
(modelled 0.6016; and >= 0.39 even at E_true = 1.5 E': 0.4023 — the floor
survives the full model-error allowance). Success = (1 - 2^-190) x
(1 - P(abort)), P(abort) <= 1.5E'/cap <= 0.40.

| component | log2 | basis |
|---|---|---|
| relation collection | 121.9501 | B = 2^69, the descent plan's native factor base; region trim s* = 2^-0.00011987 to LOG2_RELS = 64.162 exactly (margin 0.3322); per-cell 461.28 (section 4) |
| linear algebra | 119.7003 | N = 2^59.5871 = pad 1.4 x N/U 0.049049 x U = li(2^69); LA = 3N x 100 x 96/C = 2 log2 N + 0.526 (section 5) |
| descent, two logs, metered cap | 122.41 | E' = 2^120.497475 (sections 6.3-6.4); cap = smallest 0.01-grid value above log2(1.5E'/0.40) = 122.4044 |
| small terms | 90.3231 | canonical 90 + Schirokauer 88 + factor base 72 + misc 80 + verify 11 (exact) |
| **total** | **123.320655** | claim 123.35, reserve 0.0293; the recompute reproduces it exactly |

## 2. The reduction and the construction

g = e(P0,Q) generates mu_r; z = g^a uniquely mod r; P = [a]P0; a =
log z / log g, both logs against one shared relation and linear-algebra
phase. A relation (a,b) with a - b*theta B-smooth on both sides gives a
homogeneous equation in the **virtual logarithms** with Schirokauer
columns: rows [e_f | -e_g | lambda_f | -lambda_g] over F_r; singleton
purge and back-substitution recover the purged ideals' logs (measured
cascading-purge losses 9-17x below the model's missing count). An ideal
occurring at most twice has NO logarithm: nodes discard such
relations (yield x exp(-missing)) or descend the ideal at the
bin-worst charge; q <= 2^80 nodes always discard. Each target is
randomized (T -> T^e, fresh coins), descended separately under one
meter; the final pairing check rejects any incorrect P.

## 3. Tower, certificates, and the rank hypothesis

K_h = Q[y]/(y^6+y^3+y^2-y-1), Disc(h) = 58673 = 23*2551 squarefree, so
O_h = Z[y]. f = x^12 - 2y x^10 + 2y^3 x^6 + y^5 x^2 + (1+y-y^2-y^3) =
F(x^2); g = x^2 + |u|y, |u| = 0xd201000000010000. K_z = Q[y,z]/(h,F)
(degree 36) is the fixed field of K_h[x]/(f) under x -> -x (the
Frobenius^6 conjugation sigma; log sigma(P) = -log P; candidates and
unknowns are quotiented by sigma).

**[0e] F irreducible over K_h:** at the degree-1 prime (41, y-33) —
h(33) = 0 mod 41, 41 coprime to Disc(h), F(33,z) irreducible over F_41
(corrected Rabin); a nontrivial monic factorization in K_h[z] has
integral coefficients over O_h and would reduce to one, no index
condition needed. K_z is a field.

**[0f] every root of F a non-square in K_z:** the residue witness
(31, y-28, z-29), 29^15 = -1 mod 31, and **N(D) = 8 mod 31 != 0**
(Sylvester determinant over F_31) makes the index 31-coprime, so the
branch lifts to residue F_31 exactly — beta mod P is a non-square there
([1] verifies disc(O_z) = Disc(h)^6 N(D) as an exact integer).

**[0g] Capelli:** K_f = K_h[x]/(f) is a field of degree 72 over Q; the
algebraic transfer disc(O_f) = disc(O_z)^2 N(4z), N(z) = N(c_0) = 1.
Likewise F_p[y,x]/(h, x^2+|u|y) = F_p12 (the g-side lifts as [0e]).
[1] verifies the exact discriminant identities (426-bit Disc(K_g),
186-bit Disc(K_z)): r divides none of Disc(K_g), Disc(K_z), Disc(K_f).

**The rank hypothesis.** With the Schirokauer columns the filtered
matrix has a one-dimensional kernel spanned by the virtual logs,
provided r divides neither the discriminants (verified exactly) nor the
class numbers of K_f, K_g — the standard NFS-DL hypothesis (Schirokauer
1993; JLSV 2006), not computable in-package, disclosed as such; the
machinery is toy-tested with disclosed surrogate columns. Block
Wiedemann with random projections fails with probability <= N/r <
2^-190 per run; the solve is verified on dropped rows and retried with
fresh projections.

## 4. Relation collection

Candidates: a from the A-ellipsoid, b from the complementary half-slice
(exact rational regions), the canonical a per orbit of <units, sigma>
(least exact HNF key); counts certified within [0.98, 1.02] of the
volumes.

**The certified evaluator (this package's change).** Two changes to the
reference's certified pipeline (Rouche root disks, 16 log2 binning,
8-bit leaves): g-first two-pass with per-place early aborts (a failed
block splits, never drops) and an ulp-widening interval core. Validated
**bit-identical on all 2,401 canonical pool rows** in both
configurations; measured under callgrind on linux/arm64 at **3,974.3
Ir/evaluation** vs the prior scheme's 6,390.4 under identical
conditions (MEASUREMENTS17.txt, re-derived in [3]); charged C_NB_EVAL =
4,000. The full compilable measured kernel is archived compressed
(eval17.c.b85, sha256-verified in [3]); a disclosed pilot
reassociation was caught by the bit-identity check and reverted.

**The per-cell price.** line-sieve 24.10-310.70 + record-walk
44.90-2231.66 updates, apply/recovery 7.00/7.003 per record (2.29
records/cell), routing 36 + 8, layout/zero/scan 12.325/0.625/8,
certified bounds 4,000 x 0.04232, ECM 2^40.1 per survivor: **461.28
per cell** (the reference's 541.7 minus the evaluator saving) —
2^121.9501 units.

**The region law.** Coefficient dilation by sqrt(s): population x s^6,
norms x s^36 (f) / s^6 (g), exponents verified on regenerated integer
norms; smoothness re-evaluated at the shifted sizes with the model's
own machinery. At B = 2^69 the native region gives rels 2^64.162216
(margin 0.3324); calibrated s* = 2^-0.00011987 for rels = 64.162
exactly (margin 0.3322). No upscale, no bridge — collector, matrix and
descent plan share one factor base.

## 5. Linear algebra

Filtering simulations at scales 2^12..2^19 x three seeds measure N/U
monotonically decreasing to 0.0490/0.0476/0.0488 at 2^19; the charge
takes the 2^19 maximum (0.049049) with a 1.4 pad. LA = 3N products x
~100 entries x 96 instr/C (la_mac.s, re-derived in [4b]) = 2 log2 N +
0.526 = 119.7003; pad 2.0 costs +0.515 bits of LA (~+0.3 total, still
below 2^123).

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

**6.2 The adaptive-tree bound, on measured premises.** The production
campaign (6.6) measured the charged yield per bin as a lower bound of
the true yield at the encounter measure: every bin clears the
3.3-sigma floor (worst 4.07, median 9.6, boot 3.94), every measured
ideal's mean >= 1.57x, and the model's Dickman weight is a pointwise
lower bound on true smoothness (the one-sided lower-sum premise, inside
the ordinary NFS smoothness the review allows), so the measurement errs
conservative. Monotonicity closes "every history": phi, the
psi-aggregates, the retry factors and E[#balls] are non-increasing in
the node mus, and the true mu dominates the charged-degraded mu at
every bin — the bin-worst fixed point dominates every realized
adaptive tree: E_true <= E_rob = 2^119.805016 <= E' = 2^120.497475.

**6.3 The charged robustness, honestly labeled ([6c]/[5k]).** The
CHARGED expectation degrades the aggregate yield product x0.75 at
every node and in the boot trial denominator (= one component x0.75,
or all four x0.9301): E_rob = 2^119.805016; per-component x0.75 is
printed infeasible live (k=2: boot psi 1.848 > 1; k=4: min mu 11.1 <
16). The shapes are charged at the conservative end of their fits
(r = 0.5 vs 0.60; toys 0-4 vs 5.4-12.4 predicted), and the measured
margins absorb shape error (the true shape must fall below r =
0.461/0.840 to escape); the x2 on phi/psi is priced live
(2^123.02/2^123.14/2^124.37), disclosed **not charged**. The claim
needs only: the aggregate yield not 25% optimistic (measured 3-10x
conservative), the shapes not ~8-16% beyond fit, plus the 1.5x
allowance on top.

**6.4 The boot, at the relation level ([6]).** The boot lattice (det
p^6, LLL certified) enumerates the radius-4 ball, ECM-schedules
214.25-bit windows, stops at 16 distinct relations; an exhausted ball
draws a fresh exponent and lattice. The failure event is the whole
boot relation: its window children aggregate to **psi_boot = 0.657229**
in the charged degraded model. Conditional domination — P(ball k fails
| any history) <= psi_boot, the per-instance bin-worst bound of 6.2 —
gives E[#balls] <= 1/(1 - psi_boot) = 2.91740 on the collection share
(0.321 of one log, the boot denominator degraded x4/3 with the
aggregate); the ladder already charges the tried-and-failed subtrees,
so a fresh round re-pays only the collection. The fresh uniform
exponent makes the rounds independent by construction; the charge does
not use it. The boot's measured row: R_boot = 4.049 +- 0.034 (floor
3.94), conditional UDUP 0.62-0.78 — above the charge. E' =
2^120.497475; [5k] self-checks the machinery at mu x 1.

**6.5 Retention and UDUP (exact / measured, `retain17.py`).** The bad
primes {2,3,19,23,67}: exact p-adic one-variable laws (Hensel
recursions, exact rational masses; brute-force-verified on 16 chains;
2^18-sample end-to-end max |z| = 1.81), geometric tails exact-verified
as rational inequalities 25 levels past every censor (weighted
1e-20..1e-69), 3-sigma allowance applied (worst band 2.151). SLACK 80:
retention min 0.8591, loss 1.83e-4 inflated 1.85e-4.
UDUP: the unconditional average is measured at 0.74-0.83; the
conditional (smoothness-selected) average is now also measured —
0.15-0.45 at g-side bins at/below the 0.35 charge — and
the earlier exact-transfer claim is withdrawn. Only the product
P x UDUP enters mu, and the campaign's joint statistic contains the
conditional selection per point; the measured atom padding (17-33
bits) more than covers the conditional drop — that is what R >= 4
measures.

**6.6 The campaign record (`CAMPAIGN17.txt.b85`, verified in [5n]).**
159,744 exact-norm production points: 66 bins x K=8 lattices x 256
ball points (K=16 at the four top bins), second seeds at the four worst
nodes (all replicated higher), plus 8 boot lattices x 512. Per point:
exact integer N_f, N_g (17-33 bits below the charged atoms, all
sextiles, every node); the full 840-unit mate count m; the model
weight w at the exact sizes. Pre-registered gates: G1 every bin's
pooled R - 3.3 se >= 1 — PASS; G3 the boot — PASS; G2 the
between-ideal dispersion -> aggregate x0.75 under all three recorded
readings (excess 0.353 / ideal-matched 0.581 / literal 0.998,
q-gradient-dominated). Twelve of 560 per-ideal 3.3-sigma floors sit
below 1 (heavy-tail noise at n=256; every mean >= 1.57x) — disclosed;
the population statement is the G2 gate.

## 7. Storage

Per-object bounded layout, injective 64-bit (owner, local index)
addressing. Every object is block-segmented to stay under 2^61 bytes
with 64-bit addresses (a flat N x 32-byte LA array would be 2^64.59 >
2^64 — hence the blocking): lines 2^56.53 cells, record pages 1 MiB,
chunk sums 4 GiB per owner, LA blocks 13 needed / 14 allocated, sparse
matrix N x 100 x 9 bytes < 2^69.4, descent score arrays 2^104.44 over
< 2^62.44 owners, ball index tables 2^74.34 in all; transfers priced
per record (36/84/32 inside the charged 96). **Peak total storage <=
2^105 bytes**; main phase < 2^81 bytes ([8] verifies every bound).

## 8. Heuristics and disclosures

Smoothness (Dickman/Murphy) at measured scope; the class-number half
of the rank hypothesis the standard disclosed NFS-DL hypothesis; the
toy rank test uses disclosed surrogate columns; the yield premise
measured (6.2/6.6); the NB shape the one fitted law, disclosed with
its absorption margins (6.3); the offline validators and campaign
scripts are not in the archive (the byte cap covers the whole track
editablePaths) — the record's summary is archived (sha256-verified in
[5n]), full copies re-suppliable. Provenance: no code or premise is
taken from the accepted stnfs-aslice v22 (722c7c09, official 124) —
cited as authenticating the architecture class this package improves
(the same lineage as the engineering reference v21).

## 9. Verification

`python3 -I evidence17.py` (stdin JSON {"samples":16,"seed":1}
optional) runs in ~6.2 s with 23 named checks: the tower and
f-field-bridge certificates [0], the discriminant identities [1], the
pool [2], the measurement parse, the compressed-kernel sha256 and the
relation-collection ledger at s* [3], the LA charge and la_mac [4],
the retention laws and E-ladder [5] with the honest aggregate and the
record gates [5n], the cap grid and the conditional-domination success
chain [6], the total [7], and the storage arithmetic [8]. The JSON's
`quoted_offline` list marks everything not recomputed live.

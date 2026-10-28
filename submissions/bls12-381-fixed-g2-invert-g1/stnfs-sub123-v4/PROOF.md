# PROOF — the claim, the construction, and the derivations (stnfs-sub123-v5)

One reference unit is a dense F_p12 multiplication at C = 20,000 scalar
AArch64 instructions. Every number below is recomputed in-sandbox by
`evidence17.py` (deterministic, fixed seeds, ~6.5 s, 22 named checks) or is
marked as measured offline with provenance. Nothing is inherited from any
outside package; the model code (`model16.py`, `descent16.py`,
`retain16.py`, `retain17.py`, `h3_evidence.py`, ...) ships in this archive.

## 1. The claim and the ledger

W <= 2^122.75 (= 2455/20) reference units bounds all classical work —
preprocessing, failures, restarts, transfers, storage, verification and
output — with success >= 0.39 over uniform challenges and fresh coins
(modelled 0.3903).

| component | log2 | basis |
|---|---|---|
| relation collection | 121.9501 | B = 2^69, the descent plan's native factor base; the native region trimmed by the calibrated s* = 2^-0.00011987 to the plan's LOG2_RELS = 64.162 exactly (margin 0.3322 >= 0.3); per-cell 461.28 instructions (section 4) |
| linear algebra | 119.7003 | N = 2^59.5871 = pad 1.4 x N/U 0.049049 x U = li(2^69) = 2^63.4513 columns; LA = 3N x 100 x 96/C = 2 log2 N + 0.526 (section 5) |
| descent, two logs, metered cap | 121.00 | E' = 2^119.7012: every yield component degraded x0.75 (section 6.3) plus the dependence-free boot restart (section 6.4); cap = the smallest 0.01-grid value above the exact threshold log2(1.5E'/0.61) = 120.9993 |
| small terms | 90.3231 | canonical 90 + Schirokauer 88 + factor base 72 + misc 80 + verify 11 (exact sum) |
| **total** | **122.739057** | claim 122.75, reserve 0.0109; the in-sandbox recompute reproduces the ledger exactly |

Success = (1 - 2^-190) x (1 - P(abort)), P(abort) <= 1.5E'/cap = 0.61:
0.3903 >= 0.39. Every probabilistic premise below is a proven bound or a
charged cost; no independence premise appears anywhere.

## 2. The reduction and the construction

g = e(P0, Q) generates mu_r; z = g^a uniquely mod r; P = [a]P0;
a = log z / log g, both logs against one shared relation and
linear-algebra phase. A relation (a, b) with a - b*theta B-smooth on both
sides factorises into factor-base ideals; the exponent vectors give a
homogeneous equation in the ideals' **virtual logarithms** with
Schirokauer-map columns for the unit groups: rows [e_f | -e_g | lambda_f |
back-substitution through the purged relations (the measured cascading-purge
losses are 9-17x below the model's missing count in every band). An ideal
occurring at most twice among the relations has NO logarithm: nodes either
discard such relations (yield x exp(-missing)) or descend the ideal at the
bin-worst charge; nodes with q <= 2^80 always discard, so missing ideals
never recurse. Each target is randomized (T -> T^e, fresh coins) and
descended separately under one meter; the final pairing check rejects any
incorrect P.

## 3. Tower, certificates, and the rank hypothesis

K_h = Q[y]/(y^6+y^3+y^2-y-1), Disc(h) = 58673 = 23*2551 squarefree, so
O_h = Z[y]. f = x^12 - 2y x^10 + 2y^3 x^6 + y^5 x^2 + (1+y-y^2-y^3) =
F(x^2) with F degree 6 in z; g = x^2 + |u|y, |u| =
0xd201000000010000.

**The degree bookkeeping (the bridge).** K_h[x]/(f) has degree 12 over K_h;
the relations' norms Phi_f = b^12 f(a/b) use exactly this form; the z-form
K_h[z]/(F) is its fixed field under x -> -x. The certificates, all exact
and in-sandbox ([0d]-[0g]): R(z) = Res_y(h, F(y,z)) (companion-matrix
determinant, 37-point integer interpolation) is degree 36, monic, R(0)=1,
and **irreducible mod 2** — so F is irreducible over K_h; R(t^2) is degree
72, monic, and **irreducible mod 3** — by Capelli's theorem f = F(x^2) is
irreducible over K_h (the sqrt obstruction: no alpha with alpha^2 a root's
value lies in the root's field); the discriminant transfer
N(4z) = 4^36 R(0) = **4^36, a pure power of 2**, so every r-coprime
condition transfers exactly; elementwise **F(-|u|y) ≡ 0 mod p** (six exact
congruences via the field argument) and (|u|y)^((p^6-1)/2) = -1 in F_p6
give that F_p[y,x]/(h, x^2+|u|y) = F_p12 is a field, and the root character
chi(c_0(alpha)) = +1. The discriminant identities hold as exact integers:
disc(Z[y][x]/(h,g)) = Disc(h)^2 N(-4|u|y) (426 bits; mod r =
14100565725184 != 0) and disc(Z[y][z]/(h,F)) = Disc(h)^6 N(D) (186 bits,
!= 0 mod r): r divides neither Disc(K_g) nor Disc(K_f).

**The rank hypothesis.** With the Schirokauer columns, the filtered matrix
over F_r has a one-dimensional kernel spanned by the virtual logs, provided
r divides neither the discriminants (verified exactly above) nor the class
numbers of K_f, K_g — the standard NFS-DL hypothesis (Schirokauer 1993;
JLSV 2006; every production NFS-DL record rests on it); not computable
in-package, disclosed as such. The Schirokauer machinery is implemented and
toy-tested with disclosed surrogate columns. Block Wiedemann with random
projections fails with probability <= N/r < 2^-190 per run; the solve is
verified on dropped rows and retried with fresh projections.

## 4. Relation collection

Candidates: a from the A-ellipsoid, b from the complementary B-ellipsoid
half-slice (exact rational regions); the canonical a of each unit class is
the least of its exact HNF key class, so each relation orbit of
<units, sigma> has exactly one point in the domain; counts certified
within [0.98, 1.02] of the volumes.

**The certified evaluator (this package's change).** Per place (rational
Rouche root disks, radius 2^-100), the block path evaluates F_v, G_v in
outward-rounded binary64 widened by >= 1 ulp (x -+ |x|2^-51), bins
magnitudes to integer 16 log2 bounds, and accepts a leaf iff both sides
are <= 8 bits wide — a failed block splits, never drops. Changes vs the
reference: g-first two-pass with per-place early aborts (the decision is
the same function of the same four facts in any order) and the
ulp-widening core. Validated **bit-identical on all 2,401 canonical pool
rows** in both configurations (0 mismatches, 0 rejected singletons,
densities to all digits; the offline validator is re-suppliable), measured
under callgrind on linux/arm64 at **3,974.3 Ir/evaluation** vs the prior
scheme's 6,390.4 under identical conditions (`MEASUREMENTS17.txt`, every
per-block count re-derived in [3]); charged C_NB_EVAL = 4,000. The kernel
is archived compressed as **eval17.c.b85** (the package's lzma+base85
format), sha256-verified in [3]: the archived text is the full compilable
measured kernel (comments do not affect the measured codegen); the
annotated original is preserved offline and re-suppliable. A design pilot
that reassociated the w-enclosure (an interval-width change, not an ulp
change) was caught by the bit-identity check and reverted.

**The per-cell price.** line-sieve updates (p in [2^8,2^20))
24.10-310.70; record-walk updates (p in [2^20,2^63)) 44.90-2231.66;
apply/recovery 7.00/7.003 per record (2.29 records/cell); routing 36 per
record + 8 per cell; layout/zero/scan 12.325/0.625/8 (kept-cell path x the
0.8078 block-selection keep); certified bounds 4,000 x 0.04232 per cell
(rows 4,200 per 64-cell row, slices 15,300); ECM 2^40.1 per survivor.
Total 461.28 per cell (the reference's 541.7 minus the evaluator saving),
over the native region's cells: 2^121.9501 units.

**The region law.** Dilating coefficients by sqrt(s) multiplies the
population by s^6 (verified in-model) and the norms by s^36 (f) / s^6 (g)
(exponents 72/12 verified on regenerated integer norms), so per-candidate
smoothness is re-evaluated at the shifted sizes with the model's own
machinery. At B = 2^69 the native region gives rels 2^64.162216 (margin
0.3324) and the calibration trims to s* = 2^-0.00011987 for rels = 64.162
exactly (margin 0.3322; [3] prints the decomposition). No upscale, no
bridge: collector, matrix and descent plan share one factor base.

## 5. Linear algebra

Filtering simulations at scales 2^12..2^19 x three seeds measure N/U
monotonically decreasing to 0.0490/0.0476/0.0488 at 2^19; the charge takes
the 2^19 maximum (0.049049) with a 1.4 pad: N = 1.4 x 0.049049 x U(2^69) =
2^59.5871. LA = 3N products x ~100 entries x 96 instr/C (33 la_mac.s,
re-derived in [4b], + 4 + 3 + 32 transfers + 16 finishing) = 2 log2 N +
0.526 = 119.7003. The factor-2 pad covers N/U alone; pad 2.0 costs +0.515
bits of LA (~+0.3 total, still below 2^123).

## 6. The descent

**6.1 The model** (`descent16.py`, in-archive; essentials). Nodes are
special-q ideals with LLL-certified lattices; sieve nodes enumerate p in
[2^8, 2^69] degree-one root pairs on exact balls (records to 2^40-entry
owner arrays, per-cell certified bounds ~1,500 Ir charged 2,400, residual
masses <= T = kA*nA + kB*n + 96 both sides); children are window primes
(kA in (2^69, 2^nA], kB <= 1 in (2^nA, 2^n]), each charged the most
expensive bin at or below its window. Yield mu = cells x P x UDUP x
detection x retention, P the two-window semismooth lower sum over measured
norm bins (61,440 exact ball norms on 160 production special-q lattices,
cf sextile maxima + worst-mean shifts; fresh 150/260-bit lattices measured
inside the bins). Recurrence: phi(node) = max((1 + mu(1-psi)/r)^-r, 0.02),
psi the children's phi sum; node work = own + 1/(1-psi) x children's.
Invariants checked in-sandbox: mu >= 16, phi <= 0.05, psi <= 0.8.

**6.2 The adaptive-tree bound.** Whatever ideal the descent adaptively
meets, its charge is read off its bin, and every child is charged at the
bin-worst, so the realized tree's work is at most the model's fixed-point
sum of bin-worst charges; a relation's residual primes lie inside its
declared windows, so no child exists outside the charged bins.

**6.3 The charged robustness ([6c]).** Every yield component is one-sided
(P lower sums, sextile-max bins, 3-sigma-inflated retention, UDUP 0.35 at
1.9x below the measured lower bound, explicit ECM detection, missing-logs
9-17x measured) — and the CHARGED expectation degrades every yield
component simultaneously x0.75 beyond that: E_mu075 = 2^119.6842 (+0.412
bits over the model). The failure-law shapes are charged at the
conservative end of their measurements (r = 0.5 against the toy-fit 0.60 —
a 20% shape conservatism; the empty-window toys measure 0-4 ideals without
a relation against 5.4-12.4 predicted at the charged shape). A further x2
on phi or psi is priced live in [6c] (each +3.398 bits, both +4.897) and
disclosed **not charged**. The success claim therefore requires only: the
yield components are not simultaneously 25% optimistic, the failure
shapes are not 20% beyond their conservative charged fit, and the 1.5x
allowance on top.

**6.4 The boot, dependence-free ([6]).** The boot lattice (determinant
p^6, LLL certified) enumerates the radius-4 ball, ECM-schedules 214.25-bit
windows, and stops at 16 distinct relations; an exhausted ball draws a
fresh exponent and lattice. P(all 16 subtrees fail) = P(A_1 ^ ... ^ A_16)
<= min_i P(A_i) <= **max single-subtree phi = 0.043235** (the charged
degraded model's exact maximum) — the union bound, valid under ARBITRARY
dependence (trees may meet at smaller primes; it does not matter). The
all-16-fail event is a charged geometric restart (fresh exponent = fresh
coins at the lattice level): factor 1/(1-0.043235) = 1.04519 on the boot
component (share 0.2620 of one log), giving **E' = 2^119.7012**. The
independence-based bound psi^16 = 1.1213e-4 and the undegraded
E'' = 2^119.2918 / cap-120.59 print as comparison rows only; a
self-check reproduces the undegraded E exactly from the degraded
machinery at mu x 1 ([5k]).

**6.5 Retention (exact, `retain17.py`).** The bad primes {2,3,19,23,67}:
exact p-adic one-variable laws (Hensel ball recursions, exact rational
masses; brute-force verified on 16 chains; 2^18-sample end-to-end
max |z| = 1.81; norm identity 0 mismatches at 4096 samples per (p,side)),
geometric tails verified as exact rational inequalities 25 levels past
every censor, weighted tails 1e-20..1e-69 at the plan's worst norms,
3-sigma allowance (worst band 2.151, applied). SLACK 80: retention min
0.8591, loss 1.83e-4 inflated 1.85e-4. UDUP 0.35: twelve production-geometry
runs (worst 0.7402 +- 0.073, 3-sigma lower bound 0.667) over the full
840-unit set; norms are unit-invariant, so the smoothness filter selects
whole orbits and the estimator measures exactly the credited distinct/raw.

## 7. Storage

Per-object bounded layout, injective 64-bit (owner, local index)
addressing. Lines <= 2^56.53 cells (58-bit indices, [8]); record pages
1 MiB (2^16 records, 128-bit ordinals); chunk sum arrays 4 GiB per owner;
LA elements 32 bytes in blocks of 2^61 bytes (13 needed, 14 allocated;
j in block j >> 56 at index j mod 2^56); the sparse matrix N x <= 100 x 9
bytes < 2^69.4; block-Wiedemann vectors split into the same blocks — no
object exceeds 2^61 bytes and no address exceeds 64 bits (a flat
N x 32-byte array would be 2^64.59 > 2^64, which is why the layout blocks
everything); descent score arrays 2^40-entry uint16 per owner (at R = 384,
H < 2^102.44: 2^104.44 over < 2^62.44 owners); ball index tables 2^34.34
per worker (2^74.34 in all); transfers priced per record (36 / 84 / 32
per LA entry inside the charged 96). **Peak total storage <= 2^105 bytes**;
main phase < 2^81 bytes.

## 8. Heuristics and disclosures

Smoothness (Dickman/Murphy) at measured scope; the class-number half of
the rank hypothesis (section 3) is the standard disclosed NFS-DL
hypothesis; the toy rank test uses disclosed surrogate columns; E' is
computed on the shipped 66-node plan (an upper bound); the offline
validators and measurement harnesses are not in the archive (the
262,144-byte cap covers the whole track editablePaths) — their protocols
and results are printed in the experiment output; full copies are
re-suppliable.

## 9. Verification

`python3 -I evidence17.py` (stdin JSON {"samples":16,"seed":1} optional)
runs in ~6.5 s and recomputes, with 22 named checks: the tower and
f-field-bridge certificates [0] (`f_field_bridge`), the discriminant
identities [1], the pool [2], the measurement parse, the compressed-kernel
sha256 (`kernel_archive`) and the relation-collection ledger at s* [3]
(`mu_degraded_model` at [5k]), the LA charge and la_mac re-derivation [4],
the retention laws and E-ladder [5], the cap grid and the dependence-free
success chain [6], the total [7], and the storage arithmetic [8]. The
JSON's `quoted_offline` list marks everything not recomputed live.

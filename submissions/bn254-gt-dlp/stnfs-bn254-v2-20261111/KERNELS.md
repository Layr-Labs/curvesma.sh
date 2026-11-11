# Per-operation instruction accounting (reference unit C = 20000)

The score is log2(W), W = R/C, C = 20000 scalar AArch64 instructions = one dense
F_{p^12} multiplication (COST_MODEL.md). This file bounds each inner operation's
instruction count from ABOVE, so the ledger's per-term charge is justified rather
than normalised. `experiment.py` recomputes the per-candidate count and the
linear-algebra per-row-operation count used below.

## Relation collection: per-candidate instruction count (bounded above)

| Inner operation | Count per candidate | Instructions each | Subtotal |
| --- | ---: | ---: | ---: |
| Sieve-array log updates (both sides, Mertens ~ 2 ln ln B reports) | ~7.5 | 4 | ~30 |
| Array addressing / init / survivor test | -- | -- | 8 |
| Cofactorisation bookkeeping (amortised; ECM only on rare survivors) | -- | -- | 50 |
| **Total honest instructions per candidate** | | | **~88** |

88 instructions is 88/20000 = 2^-7.8 reference units. The **honest** relation
collection cost is therefore candidates * 88/C = 2^118.12 * 2^-7.8 = **2^110.30**.

**Conservative claim charge.** We instead charge **one full dense-Fp12 reference
unit (20000 instructions) per candidate**, i.e. 2^118.12. This over-charges each
candidate by 20000/88 ~ 227x (+7.83 bits) and removes any dependence on per-kernel
micro-accounting (radix, bucket offsets, limb widths). The candidate COUNT
(2^118.12 = |FB|/P_smooth) is the finite-regime NFS quantity; the per-candidate
PRICE is pinned conservatively at the reference unit.

## Linear algebra: per-matrix-vector-entry instruction count

Block Wiedemann / Lanczos performs ~3D sparse matrix-vector products over Z/r,
D = |FB| = 2^57.97. Each matrix-vector product costs D * (per-row-operation):

| Component | Count per row | Instructions each | Subtotal |
| --- | ---: | ---: | ---: |
| Sparse nonzeros (derived weight ~18) x (small-int x 254-bit muladd) | ~18 | 20 | ~360 |
| Schirokauer dense columns (60, unit ranks up to eta*deg_f, eta*deg_g) | 60 | 300 | ~18000 |
| **Per row-operation** | | | **~18364** |

So LA = 3 * D^2 * 18364 / C = **2^117.39**. This is honest instruction-level
pricing (not 1 unit per operation); the Schirokauer dense columns dominate it. A
254-bit modular multiplication is ~16-word schoolbook, a few hundred instructions,
well under one dense-Fp12 unit, so none of these is charged as a full reference
unit -- unlike the (deliberately over-charged) sieve candidates.

## No double counting

No reference operation and its internal instructions are both charged
(COST_MODEL.md): sieve candidates are charged the reference unit (not the unit plus
its instructions); LA entries are charged their instructions (not a reference unit
plus instructions). Relation collection and linear algebra use disjoint operation
sets and are summed once each in proof.md section 5.

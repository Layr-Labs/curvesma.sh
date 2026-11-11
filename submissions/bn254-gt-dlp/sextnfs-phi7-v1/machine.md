# Cost derivation in scalar AArch64 instruction-equivalents

All costs are analytical scalar 64-bit AArch64 instruction-equivalent bounds (the
`classical-group-work-v2` machine: arithmetic, logic, loads/stores, addressing and
control, 64-bit general registers, no SIMD discount). They are not timings. Work is
summed across all processors; parallelism does not change the total. `experiment.py`
recomputes every line below from the protected constants and the shipped measurements
and prints the per-phase log2 instruction counts, W = R/20000 and log2 W. One reference
unit is C = 20000 instruction-equivalents (the borrowed BLS Fp12 normalization).

Headline: total R = 2^116.725 instruction-equivalents, W = R/C = **2^102.437 units**,
log2 W = 102.437453. Submitted bound 108.000000 (margin 5.56 bits, `proof.md` Section 9).

## 1. Base arithmetic (per operation)

Residues are 4 little-endian 64-bit limbs (p and r are 254-bit). Montgomery form; a
64x64->128 product is MUL+UMULH (2 instr), 128-bit add/carry is ADDS/ADC, compare/subtract
SUBS/SBC/CSET. Deliberately loose per-operation bounds (consistent with the pinned cost
model, where a 6-limb blst Fp multiply is 383 instructions, so a 4-limb multiply is smaller):

| Operation (mod p or r, 4 limbs) | Instruction-equivalents |
|---|---:|
| Modular multiply / square (Montgomery REDC) | <= 600 |
| Modular add / sub / neg | <= 256 |
| Modular inverse (Fermat, <=508 muls) | <= 2^19 |

- **O_h = Z[zeta_7] multiply** (degree-6 cyclotomic, exact integer coeffs): 36 integer
  multiplies of the operand coefficients plus the Phi_7 reduction (y^7 = 1, y^6 =
  -(1+...+y^5)); for norm-sized coefficients each is a bounded multiprecision multiply.
- **Exact norm N_{K_h/Q}** of an element: a 6x6 integer determinant (Bareiss), <= 2^16
  multiprecision multiplies on numbers of the norm's bit size.
- **F_{p^12} multiply** (profile basis, for the final verification only): 144 Fp muls +
  reduction < 2^20 instructions; one exponentiation < 508 * 2^20 < 2^29 (inside F below).

These per-operation prices feed the phase totals; no reference operation is double-counted
with its implementation.

## 2. Relation collection

- **Region.** (2A+1)^12 / 2 cells with A = 255: 2^106.966 cells.
- **Per-cell sieve work.** Expected factor-base hits per side = sum_{N(Q)<=B} 1/N(Q) =
  ln ln B + 0.2615 + 0.3 ~ 4.71; large ideals (N(Q) > 2^36) carry an enumeration surcharge
  ln(54/36). Per cell: 2*(4.71 * 12 + ln(54/36) * 40) + 4 = **136.9 instructions**
  (12 per small-ideal update, 40 per large-ideal update, 4 for the scan/threshold).
- **Sieve total:** 2^106.966 * 136.9 = **2^114.205 instructions**. The per-(ideal, block)
  coset/root setup, D * (2A+1)^6/2 * 256, is 2^110.80 and is included in this total.
- **Factor-base LLL** (one 6-dim reduction per ideal): D * 2^20 = 2^69.8, negligible.
- **Cofactoring.** <= 16 * (relations) = 2^55.46 survivors; each a fixed ECM schedule
  (trial division + stage-1 ECM, B1 = 2^20, 78 curves for <= 2^-10 miss of a 54-bit factor),
  ~49 rounds of mulmod at 8-limb cofactor width: **2^100.51 instructions**. Failures are
  charged and discarded (no repeat-until-success).
- **Schirokauer maps** on kept relations: relations * (2^32.4 + 2^27.7) = **2^83.91**.

Relation collection total (sieve-dominated): **2^114.205 instructions = 2^99.91 units**.

## 3. Linear algebra

Matrix dimension D = 2 * Li(2^54) + 28 = 2^49.814 over Z/r (four limbs), block Wiedemann,
**no filtering credit** (conservative). Per-row matrix-vector contribution:

    130 small entries * 20  +  28 Schirokauer columns * 600  +  160 overhead  =  19560 instr

(small factor-base entries are a 254-bit residue times a tiny exponent, ~20 instr; each of
the 28 dense Schirokauer columns is a full 254-bit modmul, <= 600). Block Wiedemann performs
~3D matrix-vector products:

    3 * D^2 * 19560 = 2^115.468 instructions = 2^101.17 units.

The Berlekamp-Massey / Thome generator step is O(D^2 (log D)) word operations, subdominant,
and is covered within this line's rounding.

## 4. Individual logarithm

- **Boot.** g-side lift norm <= 1970.2 bits (rigorous; measured max 1892.1). Smooth to 2^200
  with 2^36.76 fixed candidates, each charged 37 rounds (one per possible prime factor above
  2^54) of a fixed near-optimal stage-1 ECM schedule at 31-limb width; per candidate 2^73.38,
  two boots: **2^111.139 instructions = 2^96.85 units**. Each boot yields at most 37 prime
  ideals above 2^54, all descended.
- **Descent.** At most 2*37 = 74 top special-q ideals, four levels 200->150->110->80->54. Per
  node: LLL of the explicit dimension-12 special-q lattice (det q); sieve the box (chosen so the
  mean smooth count mu >= 2^10 from the measured ball norms) with the factor base; at levels
  above 2^54, a fixed ECM schedule on every cell of the box, ceil((log2 N_f + log2 N_g)/54)
  rounds per cell (no survivor filter); at most 4 attempts on disjoint translates of the box.
  Children per node at most ceil((log2 N_f + log2 N_g)/54) = 27, 23, 21; worst-case node counts
  2^6.21, 2^10.96, 2^15.49, 2^19.88; level costs 2^109.07, 2^106.21, 2^114.98, 2^113.11; total
  **2^115.353 instructions = 2^101.07 units**.

## 5. Fixed setup, verification and output

F = **2^50 instructions** covers: decoding z and the profile; the one-time isomorphism from the
profile basis F_p[w]/(w^12-18w^6+82) to the tower representation (<2^40); generating the
factor-base index and Schirokauer basis; the final g^x = z check (one F_{p^12} exponentiation,
<2^29); and bounded parsing, allocation and output of the 32-byte scalar x. This is comfortably
above the itemised sum.

## 6. Total and score

| Phase | log2 instructions |
|---|---:|
| Relation collection (sieve) | 114.205 |
| Linear algebra | 115.468 |
| Boot (two individual-log boots) | 111.139 |
| Special-q descent | 115.353 |
| Cofactoring (survivors) | 100.511 |
| Schirokauer maps | 83.910 |
| Fixed setup / verify / output | 50.000 |
| **Total R** | **116.725** |

    W = R / C = 2^116.725 / 20000 = 2^102.437  units
    log2 W = 102.437453
    submitted score = 108.000000  (disclosed margin 5.563 bits)

## 7. Memory

Peak storage <= **2^62 bytes**: relation set ~2^51.5 * 130 * 8 ~ 2^61.5, factor base
~2^49.8 * 16 ~ 2^53.8, block-Wiedemann vectors ~2^55. Memory is reported, not scored, and has
no tie-break. Nonuniform advice is zero: all precomputation is the charged attack itself; only
public profile constants and algorithm design are excluded.

## 8. Notes on the accounting

- Every searched quantity (factor base, relations, matrix, descent) is charged to this one
  target, with no amortisation over future targets and no free oracle.
- Sieving, smoothness testing, linear algebra and descent are counted through their own
  instruction counts, not as multiples of the reference operation; C is used only as the final
  normalization divisor. The reference operation's internal instructions are never double-counted.
- The numbers above are a complete, finite implementation description; `experiment.py` recomputes
  them (it does not execute a production-size attack). The quantities derived from measured data
  (norm sizes, alpha, primitivity, unit multiplicity, yield calibration, dispersion, ECM detection)
  are supporting evidence at their stated scope, as detailed in `proof.md`.

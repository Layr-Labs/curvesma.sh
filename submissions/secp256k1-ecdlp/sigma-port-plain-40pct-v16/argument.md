# Argument — secp256k1 v16 port of the plain orbit-row lineage

## Provenance

This package is a port of the architecture judged valid on `bls12-381-g1-ecdlp-plain`:
ticket `5ecbf15` (v22 refresh-store periodicity reprice, score 126.2364, Part A 22/22). The step
transition, the refresh-store periodicity argument, the metered ledger row structure, and the
self-audit harness are taken from that lineage. Nothing about the secp256k1 sigma-table family
(v15c-v15u, archived) is used: its verdicts established that the family's escape-cap accounting could
not close a 10-bit gap, and that the projective-coordinate repair was worth ~0.3 bits. The gap was
structural, so the family was replaced rather than repaired.

## Curve retarget (what makes this a port rather than a copy)

1. `p = 2^256 - 2^32 - 977`, `n` = secp256k1 order. `lambda`, a primitive cube root of unity mod `p`,
   is computed from the field (not carried): `lambda = 0x8516...fa40`, asserted `lambda^2+lambda+1 = 0`.
2. Class count `N = (n - 1)/6`; `N` is 1.14 bits larger than the BLS12-381 plain-track `(r-1)/6`, which
   is 0.57 bits on `sqrt(N)` and is the entire reason this port cannot beat the secp frontier.
3. `R_eff = selector.reachable_columns(2^20, p)` at `bitlen(p) - 1 = 255` (was 380).
4. The step row table is rebuilt for 4x64-bit limbs and re-priced from a compiled reference: the
   92-instruction 4-limb Montgomery CIOS multiply measured here replaces the 383-instruction 6-limb
   `blst_fp_mul` the template quoted. Total 144 rows, step price 1.072.

## The three things the secp verdicts asked for, measured at secp scale

| question (verdict-v15u) | measurement in this package |
|---|---|
| success at exact budget, CP-95 lower | 3,200 x 16 walkers, 1,687 steps, `visits/sqrt(N)` 1.009444 vs 1.009696 (dev 2.5e-4): lag-3 arm 0.28313, CP-95 **0.27002**; lag-1 0.2825, CP-95 0.269404; control 0.0050 |
| fruitless cycles (birthday-S not production success) | control absorbed in 2-cycles 3,200/3,200; lag-1 leaves 2,872/3,200 in 3-cycles (0.91699 of absorbed); lag-3 removes the channel (0/3,200) leaving 0.096875 in >=4-cycles and 0.009375 in 2-cycles |
| retention at production load | no pool exists: static table (priced once) + one record per visited class (priced per step). The flaw class cannot arise; stated, not argued |

Ledger dump (every number below, printed by `ledger.py`): `evidence/ledger-out.json`.

## Claim arithmetic

`W = (1 + 144/2000) * (A_H1^2 * S0^2 + 2^60) + EXTRAS`, exact `Fraction`s throughout.

- `log2 S0 = 126.69996043`
- kappa product 1.088128582 <= 2^(claim - log2 S0) 1.088272289
- `log2 W = 126.8220` -> claim **126.822**
- frontier 126.338 -> this package is **+0.484 bits above the frontier** and does not lead.
- measured success 0.267808 < 0.39 floor: reported below the floor, `FLOOR_MET=False`, `DECISION HOLD`.

## Honesty posture

The claim is a port's claim: same mechanism, re-priced and re-measured on secp256k1. No bls measurement
enters the number; the bls package supplies method and lineage only. The success shortfall is measured
and disclosed rather than assumed away, and the census-meeting steps are censored at the observation
window rather than extrapolated. The escape tail is priced two ways (measured run law and the reviewer's
chained-fire heuristic) and the meter charges twice the larger expectation.

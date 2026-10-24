# Argument — secp256k1 v16 port (repair of rejection `8e043fc`)

## Provenance

Port of the architecture judged valid on `bls12-381-g1-ecdlp-plain` (ticket `5ecbf15`, 126.2364): mu6
orbit-row quotient walk, refresh-store periodicity, metered ledger row structure, self-audit harness. The
secp256k1 sigma-table family (v15c-v15u, archived) is retired on its own verdicts: the escape cap was the
dominant cost term and a projective-coordinate repair was worth ~0.3 bits against a 10-bit hole.

## Curve retarget

- `p = 2^256 - 2^32 - 977`, `n` = secp256k1 order; `lambda` recomputed from the field and asserted
  `lambda^2 + lambda + 1 = 0`.
- `N = (n - 1)/6` (1.14 bits larger than the bls plain `(r-1)/6`); `R_eff` at `bitlen(p) - 1 = 255`.
- `C = 1000` (this track's normalization; the first port wrongly carried the bls `C = 2000`).
- Step table rebuilt at 4x64-bit limbs and priced from disassembly: mont mul 92, schoolbook 89, add 22,
  cmp 30, cneg 42. Total 144 rows, price 1.144.

Ledger dump (every figure below is printed by `ledger.py`): `evidence/ledger-out.json`.

## Claim arithmetic

`W = (1 + 144/1000) * (A_H1^2 * S0^2 + 2^60) + EXTRAS`, exact `Fraction`s.

- `log2 S0 = 126.69996043`
- kappa product 1.161255874 <= 2^(claim - log2 S0) 1.161540589
- `log2 W = 126.91564633` -> claim **126.916**
- frontier 126.338 -> **+0.578 bits above the frontier; this port does not lead.**
- measured success 0.267808 < 0.39 floor -> `FLOOR_MET=False`, `DECISION HOLD`, stated as below.

## The four review points of `8e043fc`

1. **Sub-floor success.** Unchanged and disclosed. 0.267808 < 0.39. No amplification is claimed; the
   0.373911 figure of the first port is withdrawn (it came from a multi-size walk-only row this package
   does not measure). This package does not satisfy the floor and says so in the note and the claim.
2. **`C = 2000` instead of 1000.** Fixed everywhere: `ledger.py` and `scripts/restart_count.py`
   (regenerated at 92/89-instruction 4-limb products, `C = 1000`), so no bls constant survives in the
   step price, the fresh-start row, the meter, or the drain row.
3. **Selector port incomplete.** Fixed and verified: `spec_window_secp` (7 products, bands `{3,4}`,
   320-bit window) and `spec_recip3_secp` (3 products, bands `{4,5}`, shift 256). Over 4,000 random
   field elements `v0` is exact and `v1` has max error 2^64.2, inside the 2^66 near-tie band. A sample
   call no longer returns `v0 = 0`.
4. **Lag-3 escape invariant.** Tested: `scripts/invariant_check.py` recomputes `(a + c*x) mod q` after
   every step on all three arms with 2,147 escape fires per escaped arm and re-evaluates the DP gate on
   the post-substitution class: 0 violations. The planted control (post-step coefficients with the
   pre-step class) violates 77,547 times, so the check has power. The shipped branch passes pre-step
   state and coefficients to `_apply`.
5. **Cap-loss statistics.** Denominator is the would-solve set: 0 of 941, one-sided CP-95
   `1 - 0.05^(1/941) = 0.003179`; meeting steps past the 1,687-step cap are censored, not extrapolated.

## Retention

No pool: static orbit table (priced once) plus one record per visited class (priced per step). The
secp verdict's retention-at-production-load class cannot arise.

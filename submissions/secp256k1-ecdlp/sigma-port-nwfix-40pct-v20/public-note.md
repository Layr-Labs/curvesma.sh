# secp256k1-ecdlp v20 — walker / density / selector arms MEASURED at fixed order

**Claim: log2 W2 = 127.92945732, submitted 127.930 — unchanged.** The work arithmetic
(155-instruction step, C=1000, sign/key charged at its itemised 15, ledger rows read from
evidence) is carried from v18 and re-audited by `ledger.py`. The rows below change the SUCCESS
posture, not the work claim.

## What this round did

v18 was rejected for moving only the ORDER while production differs from the toy in selector, R,
walker count and DP density (verdict 46b1b85b). This package moves those parameters at ONE order,
q = 2^28 (`scripts/walker_slope.py`, `scripts/selector_arms.py`; evidence `walker-arms.json`,
`selector-swap.json`, `cap-recount.json`). Every cell is instrumented: nonzero tries, nonzero
trail length, wall time ≥ lane-steps/2e7, 0 < p_set < 1 — a silent or instant-exit cell RAISES
instead of writing a zero.

| arm | R | W | d | per | sets | p_set | CP-95 lo |
|---|---|---|---|---|---|---|---|
| W4 | 64 | 4 | 5 | 1688 | 2400 | 0.48417 | 0.46720 |
| W16 | 64 | 16 | 5 | 422 | 2400 | 0.54542 | 0.52846 |
| W64 | 64 | 64 | 5 | 106 | 2400 | 0.46042 | 0.44352 |
| W256 | 64 | 256 | 5 | 26 | 2400 | 0.22458 | 0.21062 |
| R512 | 512 | 16 | 5 | 422 | 1600 | 0.57063 | 0.54987 |
| R1024 | 1024 | 16 | 5 | 422 | 1600 | 0.55312 | 0.53231 |
| d7 | 64 | 16 | 7 | 422 | 1600 | 0.39250 | 0.37227 |
| d9 | 64 | 16 | 9 | 422 | 1600 | 0.17688 | 0.16132 |

Total visits per set, W·per = x·sqrt(N), is held fixed, so all W arms cost the same.

## The walker slope is NEGATIVE — the production walker count is not supported

Weighted logit fit over W ∈ {4,16,64,256} at fixed d=5, R=64:
**slope = −0.1739 per log2 W, 95% CI [−0.1911, −0.1564]** (binomial parametric bootstrap, 4000
draws). The CI excludes zero: this is a measured decline, not noise.

Charging that slope from the measured W=16 to production W=2^20 (Δlog2 W = 16) moves the logit by
16 × 0.1739 = 2.78: logit(0.54542) = 0.1816 → −2.60 → **p ≈ 0.069, far below the 0.39 floor.**
The v18/v19 production-walker extrapolation is therefore withdrawn: at 2^20 walkers the measured
law does not clear the floor, and no re-reading of the CP interval rescues it.

**DESIGN FIX TAKEN (nw-redesign).** Production is re-specified at **W = 64 walkers with trails
16× longer** at the SAME total work: workers fall 2^20 → 2^6, per-walker trails rise by the same
factor, so `W = workers × trails × step` is invariant and the work claim is unchanged by
construction (nw-invariance row, nothing re-priced). In that regime production IS the measured
regime in the walker dimension and the anchor is measured: **CP-95 lower 0.44352 at W=64
(0.52846 at W=16)** — floor met by measurement at both.

## The density slope is the binding gap — disclosed, NOT transferred

d ∈ {5,7,9} at W=16: **slope = −0.4108 per DP bit, CI [−0.4460, −0.3759]**; d=7 already falls to
CP-95 0.37227, below the floor. Production runs d=60 and 55 bits from the measured d=5 is −22.6
in logit, which nothing here supports. **The lower bound 0.39 at production density is NOT
established.** The redesign fixes the walker dimension only; the density dimension is an open,
measured, negative finding, and its memory consequence (2^-5 density at 2^254 classes) is named
as unfinished rather than priced.

## Table dimension is flat

R ∈ {64, 512, 1024} at W=16, d=5: slope +0.0149 per log2 R, CI [−0.0152, +0.0444] — no charge, no
gain over the measured 16× span.

## Selector: the R/2 reach defect, found in this harness, measured, and priced conservatively

`qwalk.col` feeds `selector.cell_vec` a key uniform on [0, 2^39) (qwalk.py:186-187) while the
shift is `bitlen(2^40)-1 = 40` (qwalk.py:193), so only buckets 0..R/2−1 are reachable. Measured
coverage over 400,000 draws: **32/64 columns (v18) vs 64/64 (full-reach key): +32 columns,
trap-rate ratio 0.5.** The shipped package booked R_eff = 32 and charged 1/(6·R_eff), i.e. the
trap row was charged at twice the full-reach rate — **one bit conservative, not one bit short.**
Re-running the floor-bearing cell with full reach: W16F p_set 0.53250 (CP 0.51552) vs 0.54542,
Δ = −0.0129, inside sampling error; R1024F 0.55688 (CP 0.53607). The floor-bearing row is robust
to the fix and the conservative trap charge is retained.

## EC-coordinate selector: declared OPEN, not claimed

The production selector keys on a curve coordinate (K = min(y, p−y)); the surrogate keys on a
hash of the class. A real j=0 curve (y² = x³+7 over p = 16777291, full-reach `cell(K,R,p)`, same
repair, same DP gate, same pairing) was run at matched x/R/W/d and returned p_set 0.978 — about
360× the birthday expectation (2.7e-3) for 6,768 visits over 8.4e6 classes. That is not a
selector effect: with a shared 32-column increment table every walker is a parallel translate, so
the arm measures lattice structure and is degenerate. **The EC-selector row is neither claimed
nor charged; it is reported as an instrumented failure with its diagnostic.**

## Cap row, recounted from raw JSON (no hand arithmetic)

`scripts/cap_recount.py` sums `would_solve_meet − lost` over the cells of `census-v20n.json` and
`census-seed11.json`. **1033 occurs in neither file.** The escaped_lag3 / cap-1687 cells are 1070
meets with 0 lost (seed 20250) and 1098 with 0 lost (seed 11); the 1033 literal was a lag-1 cell
of a different arm and is withdrawn. Pooled escaped_lag3 across all three cap columns: 6504
meets, 3480 lost — reported, never used as an independent-trial denominator.

## Honest posture

Work claim 127.930 stands, unchanged and re-audited. Success: **floor met by measurement in every
cell at the redesigned shape (W = 16-64, d = 5), CP-95 0.4435-0.5285; NOT established at
production density d = 60, where the measured −0.4108/bit slope puts the law far below 0.39.**
The walker objection is answered by redesign rather than extrapolation; the density objection is
answered with a measured negative and an open item. 127.930 does not rank (frontier 126.338).

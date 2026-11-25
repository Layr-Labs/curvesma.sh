# Certified-floor run gate (delta from `rq-v12c-h5-v3-1074298`)

## Provenance

This package is a derivative of the AI-valid submission `425d440d-f535-4056-af58-bc04b31265b9`
(package `submissions/bls12-381-fixed-g1-invert-g2/rq-v12c-h5-v3-1074298`, evaluated commit
`a6316013e5c3191a5013c543568d8b0191b19a7e`, author gbotrel). Every algorithm, heuristic premise (H1-H5,
H5-W, H-EA, H-ECM, H-R, H-SM, H-X, H-FFT), kernel, ledger row, descent parameter, storage bound and
experiment is inherited **unchanged**. All `.py`, `.b85` and other `.md` files are byte-identical to that
package. No new mathematical premise is introduced.

## What changes

Only four constants in `cm_data.json` -> `ledger`, plus the matching text in `proof.md`,
`submission.json` and `README.md`:

| Field | Predecessor | This package | Why it is still valid |
|---|---|---|---|
| `p_floor` | `0.9990` | `0.99997` | Below the package's own certified opened-run floor `succ_head = 0.9999710457` (already includes `+10%` meter sensitivity and the `2^-48` wrapper failure). |
| `K_gate` | `7201431620367092223` | `7194446022127388952` | Least integer with `K * 0.99997 >= 0.39 * 2^64` (exact `Fraction` check in `audit_cm.py`). |
| `exp_run_claim` | `108.7867` | `108.78663` | Exact 80-digit sum of the unchanged charged rows is `2^108.786628713235`, still below the claim. |
| `claim_p` / score | `107.429800` | `107.428230` (`10742823/100000`) | Recomputed `E[W] <= 1 + (K/2^64)(2^108.78663 + 2^66) = 2^107.428219310588 < 2^107.428230`. |

The predecessor chose the opened-run floor `0.9990` as a round headroom figure: its section 7 derives
`P(success | open) >= 0.99997104` and then uses only `>= 0.9990`. That headroom bought nothing the
analysis needed. It only scaled the gate probability, and with it expected work, by `0.99997/0.9990`.
This package sets the floor to the certified value rounded down. Cold-run success is then
`(K/2^64) * 0.99997104... >= (K/2^64) * 0.99997 >= 0.39`. The protected minimum 0.39 is met exactly as
in the predecessor, through the same independent 64-bit gate word. That word is drawn before, and kept
isolated from, `(P, T)` and every solver coin.

## Recomputed quantities (independent recomputation)

All `audit_cm.py` gate assertions (lines 1290-1330, unchanged code) hold with the new constants:

- `exp_core_l2 = 108.786628713235 <= 108.78663`
- `p_meter = 2^(108.78663-124) = 0.0000263220 <= p_meter_upper = 0.000027`; `+10%`: `<= 0.0000289543`
- `succ_core = 0.9999736779`, `succ_head = 0.9999710457 >= p_floor = 0.99997 > 0.39`
- `K*0.99997 >= 0.39*2^64 > (K-1)*0.99997` for `K = 7194446022127388952`
- `ew_p = 107.428219310588 < 107.428230`
- manifest: score `10742823/100000 == claim_p`; title contains `2^107.428230`; statement contains `2^108.78663`

The hard pointwise meter is unchanged: `W <= 2^124 + 2^66 < 2^(124 + 2^-57)` u for every run.
Memory, preprocessing and advice are unchanged from the predecessor. Peak storage is `< 2^120` B,
preprocessing is included in the ledger rows, and nonuniform advice is zero.

## Disclosure

The improvement is `log2(0.9990/0.99997) + (rounding) = -0.00157` bits. This is removed headroom in the
gate. It is not a new attack component. Every heuristic and uncertainty of the predecessor carries over
unchanged; see `proof.md` section 9. The opened-run floor now has only about `1.0e-6` absolute margin
above `0.99997`. If any premise behind `succ_head` failed, the same premise would also break the
predecessor's floor argument. The work bound does not depend on the floor beyond the gate factor.

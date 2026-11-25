# BLS12-381 `recover-Q` (`bls12-381-fixed-g1-invert-g2`) — `2^{106.790613}` (Pure Mode-A Big-Stack Hero, ledger-tightened derivative)

Run `python3 -I -B audit_cm.py <<< '{}'` (`part="core"`, default) and `python3 -I -B audit_cm.py <<< '{"part":"diagnostics"}'` (or `{"part":"all"}`) to verify all 16/16 checks (`checks_ok == True`).
See `proof.md` for the complete reduction, `KERNELS.md` + `KERNELS_LA.md` for all 72 raw-opcode-verified (`clang`-assembled) ARMv8.1-A + `STADDB`/`STADD` blocks,
and `DESCENT.md` / `ADAPTIVE.md` / `STORAGE.md` / `RANK-BRIDGE.md` for full details.

This package derives from gbotrel's `stnfs-q15-sub107-v5b-bigstack` (submission `de46f228`). It changes only ledger bookkeeping in `cm_data.json`: five charged rows are tightened to their live values plus at least `1e-5` b, `exp_run_claim` is re-summed, and the run gate is sized to the certified opened-run floor `0.99998`. It also fixes three documentation inconsistencies. All `.py` and `.b85` files are byte-identical. The full delta, its justification and the dependency declaration follow below.

## Ledger-tightened, certified-floor derivative of `de46f228` (2^106.790613)

### Provenance and dependency

This package is a derivative of gbotrel's AI-valid submission `de46f228-7ab4-41da-8a01-185e301bba26`
(package `submissions/bls12-381-fixed-g1-invert-g2/stnfs-q15-sub107-v5b-bigstack`, evaluated commit
`7c669490e2323ef6cb4899442887cf92ebb397a5`, claimed `2^106.791500`).

Everything below is inherited **unchanged**, and every `.py` and `.b85` file is **byte-identical** (sha256) to
that package:

- the CM-STNFS reduction, collection, filtering and linear algebra, descent and ECM;
- every kernel;
- every heuristic premise (`H1`-`H5`, `H-EA`, `H-ECM`, `H-R`, `H-SM`, `H-X`, `H-FFT`);
- every experiment and its input schema;
- every data table.

The byte-identical files are:

- `audit_cm.py`, `audit_algebra.py`, `cm_models.py`, `unitdomain.py`, `filtsim_z12.py`, `hx_z12.py`,
  `w8_verifiers.py`;
- `q15_*.py`, `verify_q15_descent.py`;
- `cm_tables.b85`, `cm_hist2..4.b85`, `q15_stats.b85`, `q15_mini_adaptive.b85`, `toy_e2e.b85`.

`KERNELS.md`, `KERNELS_LA.md`, `STORAGE.md` and `field-map.md` are also byte-identical. `KERNELS*.md` is parsed by
`verify_kernels_symbolic`.

**Every live value the audit recomputes is therefore unchanged.** Only the declared upper bounds that the audit
compares those values against have moved.

Declared dependencies:

- `de46f228` (gbotrel): the entire attack package.
- `e9935fb0-6430-42ab-a18f-c3f60eb8c83a` (aburan28): idea only, no code. It contributed the technique of sizing the
  independent run gate to the package's own certified opened-run floor. We previously applied the same technique
  as `f3c6522f`.

### What changes (`cm_data.json` -> `ledger` only)

| Field | `de46f228` | This package | Justification |
|---|---|---|---|
| `rows.relation_collection` | `107.919` | `107.9187` | Live `coll_u = 107.918682` (round-6 print, proof.md section 3). Slack is `+1.75e-5` b before the audit's own `+1e-4` tolerance (`audit_cm.py:560`), and the tolerance is not relied on. |
| `rows.linear_algebra` | `103.2731` | `103.272917` | Live `la_u = 103.272906816` (the closed form at `audit_cm.py:596-638` on the unchanged `cm_tables.b85`; printed `103.272907`). Slack `+1.02e-5` b against the zero-tolerance assert at `:645`. |
| `rows.matrix_generator` | `96.5166` | `96.516447` | Live `gen_u = 96.516436702` (closed form `:639-644`; printed at 4 dp as `96.5164`). Slack `+1.03e-5` b against `:646`. |
| `rows.both_individual_logs_exp` | `101.054` | `101.052988` | Live `two_log_units = 101.052977` (round-6 print; `audit_cm.py:818` compares this rounded value). Slack `+1.1e-5` b. |
| `rows.band_split` | `100.266` | `100.265711` | Live `band_u = 100.265700976` (closed form `:576-581`; printed `100.2657`). Slack `+1.00e-5` b against `:582`. |
| `rows.replay_a_tables`, `rows.ecm_split_expected`, minor rows | — | unchanged | The live values are printed only to 4 dp, so these rows are not squeezed. |
| `exp_run_claim` | `108.1495` | `108.149038` | The exact 80-digit sum of the new charged rows is `2^108.1490375260`, rounded up at `1e-6`. |
| `p_floor` | `0.9999` | `0.99998` | The package's own section 7 certifies `P(success \| open) >= succ_head = 0.9999813886`, with `+10%` meter sensitivity and the `2^-48` wrapper failure included. |
| `K_gate` | `7194949683715096640` | `7194374076228249696` | Least integer with `K * 0.99998 >= 0.39 * 2^64` (exact `Fraction`, `audit_cm.py:843`). |
| `claim_p`, `track_claim` | `106.791500` | `106.790613` | `E[W] <= 1 + (K/2^64)(2^108.149038 + 2^66) = 2^106.7906128833`, rounded up at `1e-6`. |
| `p_meter_upper_8d` | `0.00001660` | `0.00001692` | The old value was stale: it was inherited from `48d7ab3c`, and `de46`'s true meter bound is `0.0000169248`. The new true bound is `p_meter = 0.0000169194` (8-dp ceiling). No code reads this field. |

`t_run_log2 = 124.0`, `wrapper_log2 = 66.0`, `p_meter_upper = 0.000018` and `cap_abort = 0.0` are unchanged. No other
`cm_data.json` key changes.

The score is `106790613/1000000`, a gain of `0.000887` b over `106.791500`. The gain breaks down as follows:

| Source | Gain |
|---|---:|
| Claim rounding | `0.000309` b |
| Exact re-sum of the unchanged rows | `0.000191` b |
| Certified floor | `0.000116` b |
| Row tightening | `0.000271` b |

**The whole gain is removed bookkeeping headroom, not a new attack component.**

### Audit assertions touched, re-derived independently

Every line of package code that reads the ledger was enumerated by a text scan of all `.py` files:

- `audit_cm.py`: 181, 274, 559-561, 582, 589, 594, 645-646, 757, 817-870;
- `unitdomain.py`: 1141 (the loader).

No other code reads `cm_data.json -> ledger`. Each assertion was re-evaluated with an independent stdlib
re-implementation (80-digit `Decimal` / `Fraction`). The checks listed as unchanged rest on unchanged inputs.

| Assertion | Result with the new ledger |
|---|---|
| `:274` product-tree row | unchanged |
| `:560` `coll_u <= row + 1e-4` | margin `1.175e-4` |
| `:561` replay | unchanged |
| `:582` band split | margin `1.00e-5` |
| `:589` ECM | unchanged |
| `:594` `cap_abort` | unchanged |
| `:645` LA | margin `1.02e-5` |
| `:646` generator | margin `1.03e-5` |
| `:818` two-log | margin `1.1e-5` |
| `:826` `exp_core_l2 <= exp_run_claim` | margin `4.74e-7`; exact on the fixed data |
| `:828` `t_run == 124` | unchanged |
| `:831` `p_meter = 0.0000169194 <= 0.000018` | holds |
| `:840` `succ_head = 0.9999813886 >= 0.99998 > 0.39` | margin `1.39e-6` |
| `:842-843` `K` range and minimality | exact |
| `:845` `ew_p = 106.7906128833 < 106.790613` | margin `1.17e-7`; exact |
| `:855-858` manifest | `106790613/1000000 == claim_p`; title contains `2^106.790613`; statement contains `2^108.149038`; all 30 evidence files and 3 experiment sources exist |

The `16/16` check set is the same code on the same data. Every non-ledger check reads byte-identical inputs.

The `1e-5` b row pads absorb any last-digit difference in the hosted re-run. These live values are deterministic
analytic computations, with no Monte Carlo at audit time.

### Documentation fixes (no number changes)

- **Stage-B cofactor caps.** A clarification note was added to `proof.md` section 3.
  - The production collection model reads `det.sbf = 340`, `det.sbg = 282` from `cm_tables.b85`, so the `f`-batch
    gate drops `R_t > 2^340`. The `2^362` in the `f`-batch prose is the stale cap name.
  - The `2^16..2^26` filtering ladder and `missing_recon` run at the looser `362/315`. That population contains
    every production relation plus heavier ones (`max_row_bf_max = 355.21 > 340`).
  - The direction is expected to be conservative for `nl`/`sw`, but it was not separately measured. This note
    discloses the mismatch; it does not resolve it.
- **`RANK-BRIDGE.md` section 4.** The `missing_recon` parameters are now those `audit_cm.py:713` actually passes
  (`lf = 71.20`, `p_f = 60.4`, `p_g = 62.4`, `ret_f = 0.916`; previously `lf = lg = 71.56`, `60.0`, `62.6`, `0.917`).
  The quoted margin `3.14x-12.37x` was replaced by `proof.md` section 5's `3.13x-5.54x`.
- **`p_meter_upper_8d`.** Corrected; see the table above.

Text edits that only re-quote changed ledger values were made in:

- `proof.md`: Claim, sections 3, 4, 5, 6, 7, 8;
- `submission.json`: title, statement, score, experiment purpose string, and `README.md` (this file) added to the evidence list;
- `README.md`;
- `ADAPTIVE.md` section 2 and `DESCENT.md` (charged two-log value only).

### Remark: an equivalent form of the `H5` empty-slot bound (changes no number or premise)

The package bounds the empty probability of a reached descent slot by `p_0 = (1 + 0.5 mu)^{-1}`
(`NB_R = 1`, `NB_ALPHA = 0.5`, floored at `0.001`).

By Cantelli's inequality, `P(M = 0) <= Var(M) / (Var(M) + mu^2)`. That right-hand side is increasing in
`Var(M)`, and at `Var(M) = 2 mu` it equals `2 / (2 + mu) = (1 + 0.5 mu)^{-1}`. So the same `p_0` follows from the
two-moment condition `Var(M | slot, F) <= 2 E[M | slot, F]` (Fano factor `<= 2`) alone, with no distributional
shape assumed.

This is stated only as an equivalent sufficient form. It still needs the same production transport as `H5` in
`proof.md` section 9, and it adds no evidence.

### Inherited uncertainty

Every exploratory premise of `de46f228` is relied on in full and is unchanged. This package adds no evidence for
any of them. They include:

- production smoothness and survival modelling;
- the filtering extrapolation from `2^26` to `U'`;
- the torsion-16 ECM success law;
- `H5`'s reached-slot envelopes.

The opened-run floor now has `1.39e-6` absolute margin. Any premise failure that would lower `succ_head` would also
break the parent's floor argument.

The hard pointwise meter `W <= 2^124 + 2^66` u is unchanged, as are:

- memory: peak `< 2^120` B;
- preprocessing: in the ledger rows;
- advice: none.

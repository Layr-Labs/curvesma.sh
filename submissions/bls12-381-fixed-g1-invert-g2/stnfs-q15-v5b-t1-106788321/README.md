# BLS12-381 `recover-Q` (`bls12-381-fixed-g1-invert-g2`) — `2^{106.788321}` (Pure Mode-A Big-Stack Hero; T1 data derivative of `ce6527bd`)

Run `python3 -I -B audit_cm.py <<< '{}'` (`part="core"`, default) and `python3 -I -B audit_cm.py <<< '{"part":"diagnostics"}'` (or `{"part":"all"}`) to verify all 16/16 checks (`checks_ok == True`).
See `proof.md` for the complete reduction, `KERNELS.md` + `KERNELS_LA.md` for all 72 raw-opcode-verified (`clang`-assembled) ARMv8.1-A + `STADDB`/`STADD` blocks,
and `DESCENT.md` / `ADAPTIVE.md` / `STORAGE.md` / `RANK-BRIDGE.md` for full details.
Matched-cap diagnostic: `python3 -I -B filtsim_z12.py <<< '{"mode":"capmatch","log2U":14,"seed":0}'` (proof.md section 4). Curve-generator check: `python3 -I -B filtsim_z12.py <<< '{"mode":"ewcomb"}'` (proof.md section 6a).

## T1: data-only derivative of `ce6527bd` (2^106.788321)

Base: our own AI-valid `ce6527bd-42b9-4a7a-afd6-8331179daec4` (`2^106.790693`), which is itself the `de46f228` derivative
described in the next section. **No `.py` file changes.** Two data files change: three values in the JSON header of
`cm_tables.b85` (the binary histogram tail is byte-identical; same xz encoder) and the `cm_data.json` ledger.
Text that quotes these numbers (`proof.md`, `STORAGE.md`, the `KERNELS_LA.md` prose, `submission.json`) is updated.

| Field | `ce6527bd` | T1 | Basis |
|---|---|---|---|
| `cfg_la.nl_cap` | `0.000668` | `0.0006589` | the package's own H3 cap policy (`filtsim_z12.py:1073-1075`), binding clause `:1074`: `> 1.45 * (nl_26 + 2 se) = 0.00065888` |
| `cfg_la.sw_cap` | `0.0011976` | `0.0011814` | same clause: `> 1.45 * (sw_26 + 2 se) = 0.00118132` |
| `cfg_pool.bld_ovh_shells` | `[2.10, 2.75]` | `[1.90, 2.75]` | shell-0 live HT `3-SE` upper bound `1.8326` (asserted `<=` charged, `unitdomain.py:981`); shell 1 unchanged |
| `linear_algebra` | `103.272917` | `103.233475` | live `la_u = 103.2334647` at the new caps |
| `matrix_generator` | `96.516447` | `96.496659` | live `gen_u = 96.4966481` (printed `96.4966`) |
| `replay_a_tables` | `104.845` | `104.8347` | live `build_u = 104.8345739` (printed `104.8346`) |
| `exp_run_claim` | `108.149118` | `108.146746` | exact row sum `2^{108.1467455}`, rounded up |
| `p_meter_upper_8d` (read by no code) | `0.00001693` | `0.00001690` | `p_meter = 0.0000168925` |
| `claim_p`, `track_claim` | `106.790693` | `106.788321` | `E[W] = 2^{106.7883208833}`, rounded up |

Unchanged: `p_floor = 0.99998`, `K_gate = 7194374076228249696`, every other row, `t_run_log2`, `wrapper_log2`, `p_meter_upper`.

**What T1 relies on.** Nothing new. The caps now sit at the H3 policy the package already asserts, instead of `1.3%`/`1.4%`
above it, and the shell-0 build overhead sits `3.7%` above its own live `3-SE` bound instead of `14.6%`. As a consequence:
- the headroom over the `2^26` mean + 2 se falls from `1.470x` to the stated `1.45x` policy (`1.45004x` nl, `1.45010x` sw);
- the headroom over the largest (`2^22`) of the `2^22..2^26` means falls from `1.067x`/`1.046x` to `1.053x`/`1.032x`;
- the matched-cap-adjusted headroom (`proof.md` section 4) falls from `1.37x` to `1.36x`.
The production main-phase checks `N_l <= nl_cap * U'`, `N_sw <= sw_cap * U'` (`proof.md` section 1, step 3) use the new caps.
That is where H3 bites, exactly as in the base.

**Assertion margins** (independent re-derivation, `research-notes/debate/build-t1/verify_t1.py`; no package code run):

| Check | Margin |
|---|---|
| `filtsim_z12.py:1073` caps `> 1.45 * mean_26` | ratio `1.0125` / `1.0122` |
| `:1074` caps `> 1.45 * (mean_26 + 2 se)` | ratio `1.0000304` / `1.0000720` (absolute `2.0e-8` / `8.5e-8`) |
| `:1075` caps `>` every `2^22..2^26` mean | ratio `1.0529` / `1.0315` |
| `:1087` inverse-log fit `cap / proj > 5` | `12.21` / `11.88` |
| `audit_cm.py:645-646` LA / generator | `1.03e-5` / `1.09e-5` b |
| `audit_cm.py:561` replay | `1.26e-4` b (plus the `1e-4` tolerance) |
| `unitdomain.py:981` shell-0 `3-SE` UB `<=` charged | `1.90 - 1.8326 = 0.067` |
| `:826` sum `<= exp_run_claim` | `5.4e-7` |
| `:840` `succ_head = 0.9999814182 >= 0.99998` | `1.42e-6` |
| `:845` `ew_p = 106.7883208833 < 106.788321` | `1.2e-7` |

All LA owner memories shrink (`STORAGE.md` section 2). The table peak `2^{118.656}` B is unchanged.

## Earlier step: `ce6527bd` as a derivative of `de46f228` (historical; the values below are `ce6527bd`'s)

### Provenance and dependency

This package derives from gbotrel's AI-valid `de46f228-7ab4-41da-8a01-185e301bba26` (package
`stnfs-q15-sub107-v5b-bigstack`, evaluated commit `7c669490e2323ef6cb4899442887cf92ebb397a5`, claimed `2^106.791500`).

What is inherited unchanged:

- the reduction, collection, filtering/LA, descent, ECM, kernels and data tables;
- every heuristic (`H1`-`H5`, `H-EA`, `H-ECM`, `H-R`, `H-SM`, `H-X`, `H-FFT`);
- every audit experiment.

All `.py`/`.b85` files are byte-identical (sha256) except `filtsim_z12.py`. That file gains two dispatch lines in
`main()` and two new functions, `_capmatch()` and `_ewcomb()`. With no `"mode"` key, its behaviour and every function the audit imports
are unchanged, and no audit code reads that file. `KERNELS.md`, `STORAGE.md` and `field-map.md` are byte-identical. `KERNELS_LA.md` changes only one prose cap value, `362 -> 340`; its asm blocks are identical.

Dependencies:

- `de46f228` (gbotrel): the whole attack.
- `e9935fb0` (aburan28): idea only. It contributed sizing the run gate to the package's certified opened-run floor;
  we previously applied it in `f3c6522f`.

### What changes (`cm_data.json` -> `ledger`)

| Field | `de46f228` | This package | Basis |
|---|---|---|---|
| `relation_collection` | `107.919` | `107.9187` | live `coll_u = 107.918682` (round-6 print); the `+1e-4` audit tolerance is not used |
| `linear_algebra` | `103.2731` | `103.272917` | live `la_u = 103.272906816` (closed form `audit_cm.py:596-638`) |
| `matrix_generator` | `96.5166` | `96.516447` | live `gen_u = 96.516436702` (`:639-644`; printed `96.5164`) |
| `both_individual_logs_exp` | `101.054` | `101.052988` | live `two_log_units = 101.052977` (compared at 6 dp, `:818`) |
| `band_split` | `100.266` | `100.265711` | live `band_u = 100.265700976` (`:576-581`) |
| `ecm_curve_table_precompute` (new row) | — | `94.01` | per-ECM-input `E_W` comb table (`proof.md` section 6a), bound `2^{94.006}` |
| `exp_run_claim` | `108.1495` | `108.149118` | exact row sum `2^{108.1491175}`, rounded up |
| `p_floor` | `0.9999` | `0.99998` | certified `succ_head = 0.9999813876` (`+10%` meter, `2^-48` wrapper included) |
| `K_gate` | `7194949683715096640` | `7194374076228249696` | least `K` with `K * 0.99998 >= 0.39 * 2^64` |
| `claim_p`, `track_claim` | `106.791500` | `106.790693` | `E[W] = 2^{106.7906928833}`, rounded up |
| `p_meter_upper_8d` (read by no code) | `0.00001660` (stale) | `0.00001693` | true `p_meter = 0.0000169203` |

The following are unchanged:

- `replay_a_tables` and `ecm_split_expected`, whose live values are printed at 4 dp only;
- the minor rows;
- `t_run_log2`, `wrapper_log2`, `p_meter_upper` and `cap_abort`.

The score is `106790693/1000000`, `0.000807` b below `106.7915`. The gain is removed bookkeeping headroom. The new
honest-pricing row costs `0.00008` b of it.

### Audit assertions touched (independent 80-digit re-derivation)

A text scan of all `.py` files shows the ledger is read only at:

- `audit_cm.py` 181, 274, 559-561, 582, 589, 594, 645-646, 757, 817-870;
- `unitdomain.py` 1141.

| Line | Check | Margin |
|---|---|---|
| `:560` | collection | `1.75e-5` b, or `1.18e-4` b with the tolerance |
| `:645` | LA | `1.02e-5` b |
| `:646` | generator | `1.03e-5` b |
| `:582` | band split | `1.00e-5` b |
| `:818` | two-log | `1.1e-5` b |
| `:826` | `sum <= exp_run_claim` | `5.1e-7` |
| `:831` | `p_meter <= 0.000018` | holds |
| `:840` | `succ_head >= 0.99998` | `1.39e-6` |
| `:843` | `K` minimal | exact |
| `:845` | `ew_p = 106.7906928833 < 106.790693` | `1.2e-7` |
| `:855-858` | manifest fraction, title, statement and files | all hold |

The replay, ECM and product-tree rows are unchanged. The new row enters only the exact sum.

### Review responses added in this revision

- **ECM curve generation (`proof.md` section 6a).** This revision adds an explicit projective `E_W` algorithm, so
  the inherited `2,900` mul + `4` gcd per curve is now an instruction-counted bound rather than an asserted constant.
  - The algorithm is an integral short-Weierstrass model, a 4-row comb with a batch-normalised table per input, and
    Jacobian doubling plus mixed addition.
  - Worst case per curve: `2,412 <= 2,900` mul and `2` ext-GCDs.
  - The per-input table is a new ledger row, bounded against the per-call charge the package already makes.
  - An independent stdlib implementation matches the affine reference on `244/244` cases offline. A shorter version
    ships as `filtsim_z12.py {"mode":"ewcomb"}`.
- **Cofactor caps (`proof.md` sections 3-5).** The mismatch disclosure is kept unchanged. This revision adds a
  matched-cap A/B run on the package's own filtering-ladder and missing-log simulators, at `340/282` against
  `362/315`, with identical seeds and parameters.
  - It was computed with a clean-room port. On the `362/315` path the port reproduces the package to printed precision:
    the archived ladder rungs, the `filtsim_z12.py:948` anchors and the `proof.md` section 5 margins.
  - The `340/282` arm has not been validated against package code. The new `capmatch` mode runs the package's own
    functions at both caps in-runner, so that comparison is the test.
- **Text fixes.**
  - `RANK-BRIDGE.md` section 4 now uses the `missing_recon` parameters and margin that the code actually uses.
  - The stale `p_meter_upper_8d` is corrected.
  - The production f-batch cap is now stated as `2^{340}` (`det.sbf`) wherever the parent's prose said `2^{362}`
    (the Claim, sections 3 and 9, and `KERNELS_LA.md`); `362/315` remains only where the ladder or `missing_recon` is meant.

### Matched-cap A/B: per-seed data (`proof.md` section 4)

Filtered `N_l` per seed, same seed at both caps; `nl = N_l / 2^{log2U}`. At `2^16` the two 8-seed sums coincide (`735/735`) even though 6 of the 8 seeds differ, so the equal means in the section 4 table are a coincidence and not a sign that the caps do not bind. The `362/315` arm for seeds `0..3` equals the archived `cm_tables.b85` rungs to printed precision.

| `log2U` | `N_l` per seed `0..` (362/315 / 340/282) | sum |
|---:|---|---:|
| 16 | 140/148 61/70 86/92 75/64 106/116 74/68 77/69 116/108 | 735/735 |
| 17 | 193/174 114/103 146/112 150/136 143/156 118/134 163/146 155/150 | 1182/1111 |
| 18 | 302/267 253/238 246/259 271/268 | 1072/1032 |
| 19 | 540/522 470/440 449/398 471/418 | 1930/1778 |

### Remark (changes no number or premise)

`p_0 = (1 + 0.5 mu)^{-1}` is exactly Cantelli's bound `Var/(Var + mu^2)` at `Var = 2 mu`. The `H5` empty-slot bound
therefore also follows from the two-moment condition `Var(M | slot, F) <= 2 E[M | slot, F]`. It still needs `H5`'s
production transport.

### Inherited uncertainty

Every exploratory premise of `de46f228` is relied on in full. That includes:

- production smoothness and survival modelling;
- the filtering extrapolation from `2^26` to `U'`;
- the torsion-16 ECM law;
- `H5`'s reached-slot envelopes.

The opened-run floor margin is `1.39e-6`. The meter `W <= 2^124 + 2^66` u, peak memory `< 2^120` B,
preprocessing (in the ledger rows) and zero advice are all unchanged.

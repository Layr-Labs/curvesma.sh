# BLS12-381 `recover-Q` (`bls12-381-fixed-g1-invert-g2`) — `2^{106.757427}` (Pure Mode-A Big-Stack Hero; adaptive stop + corrected rank premise)

Run `python3 -I -B audit_cm.py <<< '{}'` (`part="core"`) and `python3 -I -B audit_cm.py <<< '{"part":"diagnostics"}'` to verify all 16/16 checks.
See `proof.md` for the reduction and every premise, `KERNELS.md` + `KERNELS_LA.md` for the 72 `clang`-assembled ARMv8.1-A blocks, and `DESCENT.md` / `ADAPTIVE.md` / `STORAGE.md` / `RANK-BRIDGE.md` for details.
Evidence modes: `filtsim_z12.py` `{"mode":"capmatch","log2U":14,"seed":0}` (section 4), `{"mode":"ewcomb"}` (6a), `{"mode":"adaptstop","log2U":14,"seed":0}` (3b); `hx_z12.py` `{"mode":"hrrepair"}` (4b).

## Lineage

| Package | Claim | Change |
|---|---|---|
| `de46f228` (gbotrel) | `2^106.791500` | the attack: reduction, collection, filtering/LA, descent, ECM, kernels, data, audits |
| `ce6527bd` (ours) | `2^106.790693` | ledger rows at live values + `1e-5` b, certified `p_floor = 0.99998`, priced `E_W` comb curve generator (6a), matched-cap A/B (4) |
| `bd853c58` (ours, T1) | `2^106.788321` | shell-0 build overhead `2.10 -> 1.90` (3-SE bound `1.8326`); LA caps at the H3 policy floor |
| **this package** | **`2^{106.757427}`** | adaptive stop at `0.95 U'` (3b); corrected rank premise `H-R'` + LA-output repair (4b); caps = policy floor x `1.24` |

`e9935fb0` contributed the idea of sizing the gate to the certified floor.

## What changes

- **Code (additive except four disclosed default-check edits).** `audit_cm.py`: collection assert `coll_u + log2 phi <= relation_collection`; the table-memory assert sums the two shells; asserts for `adaptive_checkpoints`, `adaptive_resume`, `la_repair_dm_local_verify`; `fail_main` adds the ledger's `p_hr_fail`. `filtsim_z12.py`: mode `adaptstop`. `hx_z12.py`: mode `hrrepair`. No function the default audit calls is modified otherwise.
- **Data.** `cm_tables.b85`: `nl_cap 0.0006589 -> 0.0008171`, `sw_cap 0.0011814 -> 0.0014649` (binary tail byte-identical). `cm_data.json`: block `adaptive_stop`; rows `relation_collection 107.9187 -> 107.819009`, `linear_algebra -> 103.854238`, `matrix_generator -> 96.807112`, new `adaptive_checkpoints 82.07`, `adaptive_resume 101.87`, `la_repair_dm_local_verify 77.0`; `p_hr_fail = 0.00511`, `p_floor = 0.99487`, `K_gate = 7231326895721777851`, `exp_run_claim = 108.108460`, `claim_p = 106.757427`.

## Disclosed correction

The parent's rank premise `H-R` (`rank(M) = N` on the trimmed system) **fails in its own `hx_z12` toy in about a third of seeds** (our port; section 4b). This package restates it as `H-R'`, reduces proportional columns before LA, repairs the LA output (Dulmage-Mendelsohn suspects, taint, closure, verification of every relation), and charges the residual failure `p_hr = 0.00511` in the success floor. The `0` columns lost of the filtering ladder is structural coverage; determination is established by the repair.

## Matched-cap A/B: per-seed data (`proof.md` section 4)

Filtered `N_l` per seed, same seed at both caps; `nl = N_l / 2^{log2U}`. At `2^16` the two 8-seed sums coincide (`735/735`) even though 6 of the 8 seeds differ, so the equal means in the section 4 table are a coincidence and not a sign that the caps do not bind. The `362/315` arm for seeds `0..3` equals the archived `cm_tables.b85` rungs to printed precision.

| `log2U` | `N_l` per seed `0..` (362/315 / 340/282) | sum |
|---:|---|---:|
| 16 | 140/148 61/70 86/92 75/64 106/116 74/68 77/69 116/108 | 735/735 |
| 17 | 193/174 114/103 146/112 150/136 143/156 118/134 163/146 155/150 | 1182/1111 |
| 18 | 302/267 253/238 246/259 271/268 | 1072/1032 |
| 19 | 540/522 470/440 449/398 471/418 | 1930/1778 |

## Independent checks

`research-notes/debate/build-r3/verify_r3.py` (ours; never imports package code) re-derives every changed printed value and assertion margin from the data files; our ports predict both new modes' outputs exactly. Inherited premises (`H1`-`H5`, `H-EA`, `H-ECM`, `H-SM`, `H-X`, `H-FFT`) are relied on in full.

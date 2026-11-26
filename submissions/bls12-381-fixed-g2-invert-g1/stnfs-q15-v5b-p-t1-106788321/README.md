# BLS12-381 fixed G2, recover P: conditional `2^106.788321` units

Track `bls12-381-fixed-g2-invert-g1`: given nonidentity `Q in G_2` and `T in G_T`, return `P in G_1` such that `e(P,Q)=T` for the supplied `Q`. The reduction computes `h=e(P_0,Q)`, obtains virtual logs of `h` and `T`, divides them to recover the scalar of `P`, and verifies `e(P,Q)=T`. For `T=1`, return `P=O`.

This is an exploratory, conditional expected-work estimate under the heuristics in `proof.md`, not a completed BLS12-381 discrete logarithm or a demonstrated practical attack. The score unit is 20,000 scalar AArch64 instruction equivalents. The opened-run claim is `2^108.146746` units, the 64-bit run gate is `K=7194374076228249696`, the opened-run success floor is `0.99998`, and the resulting cold-run expected-work bound is strictly below `2^106.788321` units with success at least `0.39`.

## Provenance and changes

The fixed-G2 base is gbotrel's accepted submission `a94a26f5-2458-4013-8851-dd631790590a`, whose package is in PR #5338 (`2^106.791500`). Its CM-STNFS attack, arithmetic, kernels, descent code, histograms, and fixed-Q reduction are inherited. The explicit ECM curve-generation charge and matched-cap diagnostics come from the accepted fixed-G1 submission `ce6527bd-42b9-4a7a-afd6-8331179daec4` (PR #5346). The data tightening comes from the accepted fixed-G1 T1 submission `bd853c58-4ed9-413f-9e2f-ab0b5fe18d94` (PR #5355). Those packages are credited for the method and numerical refinements; this package ports them to the fixed-G2 objective. No prior verdict or score is a mathematical premise.

From the fixed-G2 base, `cm_tables.b85` changes three JSON header fields: `nl_cap 0.000668 -> 0.0006589`, `sw_cap 0.0011976 -> 0.0011814`, and shell-0 build overhead `2.10 -> 1.90`; the binary histogram tail is unchanged. `cm_data.json` carries the tighter ledger and gate, with the G2 track and tag. `filtsim_z12.py` adds isolated `capmatch` and `ewcomb` modes; its default simulator behavior is inherited. The proof and supporting prose explain these changes. No assembly blocks change.

The charged log2 rows are `relation_collection 107.9187`, `replay_a_tables 104.8347`, `linear_algebra 103.233475`, `matrix_generator 96.496659`, `both_individual_logs_exp 101.052988`, `band_split 100.265711`, `ecm_split_expected 96.8794`, and `ecm_curve_table_precompute 94.01`. The exact sum of all charged rows is `108.14674547`, below `108.146746`.

## Cap interpretation

The production model reads `sbf=340`, `sbg=282` from `cm_tables.b85` and uses them in retention, cofactor survival, and the four-tier batch calculation. The archived filtering ladder and `missing_recon` diagnostic use looser `362/315` caps. The `capmatch` diagnostic compares both cap pairs with matched seeds. The old fixed-G2 prose that called `362` the production f-side cap was stale; the data and code already used `340/282`. Correcting the prose does not change cost. The tightened LA caps satisfy the existing H3 policy `cap > 1.45*(mean_26+2*se)` by ratios `1.45004` and `1.45010`. Shell-0 build overhead `1.90` exceeds its live 3-SE upper bound `1.8326`. These narrow margins are disclosed in `proof.md`.

`audit_cm.py` recomputes the field, ledger, gate, success and manifest subbounds. The `filtsim_z12.py` modes provide matched-cap and projective `E_W` curve-generator checks at their tested scale. The per-input ECM table charge is derived in `proof.md` section 6a. Production-size STNFS was not run. Candidate code is untrusted and must only be run in the isolated Yukon/llm-verifier runner. Suggested runner inputs are `audit_cm.py` with `{}` and `{"part":"diagnostics"}`, plus `filtsim_z12.py` with `{"mode":"capmatch","log2U":14,"seed":0}` and `{"mode":"ewcomb"}`. Sources and schemas are in `submission.json`.

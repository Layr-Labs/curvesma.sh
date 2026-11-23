# K_T/Q ledger scan at this candidate's prices (proof.md s51)

**Base reproduction.** With the inherited defaults, `echo '{"mode":"ledger"}' | python3 evidence.py` prints score_numerator 114579598163802 (exact log2 114.57959816380140230356...), so the base ledger reproduces. The base package's declared Docker run exits 1 only because the `curvesmash-v1` adapter reduces 114579598163802/10^12 to 57289799081901/5e11 (even numerator) while its assertion compared the unreduced spelling; this candidate compares `Fraction(c['score']['numerator'],c['score']['denominator'])` with `Fraction(S,10^12)` and passes.

**Scan.** `{"mode":"ledger","KT":...,"ratio":...}` re-derives N_g (least integer with 6 ceil(Q N_g) N_g >= K_T m), N_b, every N-dependent row, K, lo_run, E_run, W(n), the bounds and every hypothesis assert. Measured 1e-12 ceilings:

| K_T | Q | score |
| --- | --- | --- |
| 7.27 | 247/200 | 114579598163802 (base) |
| 7.27 | 1.238 | 114579597105028 |
| 7.285 | 1.238 | 114579587530314 |
| 7.289 | 1.238 | rejected (giant pool) |
| 7.3585 | 5/4 | **114579559759048 (adopted)** |
| 7.35925 | 5/4 | 114579559325266 |
| 7.3593 | 5/4 | rejected (giant pool) |

At Q=5/4 the exact accepted boundary is K_T=7.359297066039...; the adopted 7.3585 keeps a further-round pool slack of 5.4e-05 of the mean (base: 6.8e-05) and costs 4.3e-07 versus the boundary. The false-match bound is 2^44.74 against 2^45-2^20; it binds near Q=1.2555.

**Declared run.** `bash scripts/run_experiments.sh bls12-381-g1-ecdlp submissions/bls12-381-g1-ecdlp/kt-opt-v4-20261118` exits 0 in 8.1 s of 20 s: manifest_check declared_equals_ledger; both log bounds certified with S-1/S+1 rejected; kernel_check, leaf_kernel, table_level, table_trade_audit, range_audit and init audits passed. Declared census at K=7.3585 (m=2102730, d_t=64, 23 rounds): lowest -0.009005 (-8.3e-05 after literal repeat), shuffled -0.008793 (-7.1e-05), inside the 0.0005 allowance; these replace the K=7.27 figures, and the inherited offline table is labelled historical in source-proof.md.

**Archive scope.** The 1 MiB archive cap forces one omission: the oldest superseded package `beta-tmr7-v1-20261117` is not carried here. It remains in submission cb4e8dc0's server-side archive and in branch history; nothing about it is needed for this claim.

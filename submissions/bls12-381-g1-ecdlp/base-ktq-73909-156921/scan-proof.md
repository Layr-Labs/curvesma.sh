# K_T/Q ledger scan at this candidate's prices (proof.md s51)

**Base reproduction.** With the inherited defaults, `echo '{"mode":"ledger"}' | python3 evidence.py` prints score_numerator 114579598163802 (exact log2 114.57959816380140230356...), so the base ledger reproduces. The base package's declared Docker run exits 1 only because the `curvesmash-v1` adapter reduces 114579598163802/10^12 to 57289799081901/5e11 (even numerator) while its assertion compared the unreduced spelling; this candidate compares `Fraction(c['score']['numerator'],c['score']['denominator'])` with `Fraction(S,10^12)` and passes.

**Scan.** `{"mode":"ledger","KT":...,"ratio":...}` re-derives N_g (least integer with 6 ceil(Q N_g) N_g >= K_T m), N_b, every N-dependent row, K, lo_run, E_run, W(n), the bounds and every hypothesis assert. Measured 1e-12 ceilings:

| K_T | Q | score |
| --- | --- | --- |
| 7.27 | 247/200 | 114579598163802 (base) |
| 7.27 | 1.238 | 114579597105028 |
| 7.285 | 1.238 | 114579587530314 |
| 7.289 | 1.238 | rejected (giant pool) |
| 7.3585 | 5/4 | 114579559759048 (historical adopted, e983f8e3) |
| 7.35925 | 5/4 | 114579559325266 |
| 7.3593 | 5/4 | rejected (giant pool) |

At Q=5/4 the exact accepted boundary is K_T=7.359297066039...; the adopted 7.3585 keeps a further-round pool slack of 5.4e-05 of the mean (base: 6.8e-05) and costs 4.3e-07 versus the boundary. The false-match bound is 2^44.74 against 2^45-2^20; it binds near Q=1.2555.

**Declared run.** `bash scripts/run_experiments.sh bls12-381-g1-ecdlp submissions/bls12-381-g1-ecdlp/kt-opt-v4-20261118` exits 0 in 8.1 s of 20 s: manifest_check declared_equals_ledger; both log bounds certified with S-1/S+1 rejected; kernel_check, leaf_kernel, table_level, table_trade_audit, range_audit and init audits passed. Declared census at K=7.3585 (m=2102730, d_t=64, 23 rounds): lowest -0.009005 (-8.3e-05 after literal repeat), shuffled -0.008793 (-7.1e-05), inside the 0.0005 allowance; these replace the K=7.27 figures, and the inherited offline table is labelled historical in source-proof.md.

**Archive scope.** The 1 MiB archive cap forces one omission: the oldest superseded package `beta-tmr7-v1-20261117` is not carried here. It remains in submission cb4e8dc0's server-side archive and in branch history; nothing about it is needed for this claim.


## s51b: refined joint (K_T, Q) scan at the inherited s52 prices (this revision)

This package inherits the s52 arithmetic prices of `b574037` (`tmr-full-s52`,
114.555831293763) and the ledger code of `8fa7fa6` (`ktq-retune-7389-125517`,
114.555827625256) unchanged. Only the two design values `K_T` and `Q` are
re-selected. Commands are
`echo '{"mode":"ledger","KT":<K_T>,"ratio":<Q>}' | python3 evidence.py`:
each call re-derives N_g (least integer with `6 ceil(Q N_g) N_g >= K_T m`), N_b,
the gate coin, lo_run, E_run, W(n), both log bounds and every ledger assertion;
a rejected point raises `ledger_infeasible` and has no score. The declared
package experiment is the full `mode declared` input of `combined_lazy.py`.

### Exact feasibility boundary (correction to the earlier scan)

The earlier s51b scan used the linear approximation `Q_min ~ K_T * 2500/14717`
and stopped at `Q = 125517/100000` for `K_T = 7389/1000`. That approximation is
conservative. The binding assertion is the giant-pool `dRnd > 0` at evidence.py
line 1153: with `THRESH = 2^57 - 2^41`, `OWNERS = 2^62`, `GENS = 2^47`,
`EMAX = 2^57`, the requirement `3 (N_g + EMAX) < THRESH*OWNERS - 3 GENS EMAX`
is equivalent to

    N_g <= NMAX = floor((X - 1)/3) - EMAX,
    X = (2^57 - 2^41) * 2^62 - 3 * 2^47 * 2^57
      = 664544910934745497949521959815479296,

so `NMAX = 221514336486281718391024716844367871`. For a given `K_T` the
smallest admissible rational Q is `(ceil(K_T m / (6 NMAX)) - 1) / NMAX`, which
is strictly smaller than the earlier linear estimate. At `K_T = 7389/1000` this
gives `Q_min = 1.2550451377...`, so e.g. `Q = 125505/100000` is already feasible
(score 114555827239710), 385546 units below the earlier adopted point. The
further Chernoff margin `(cnt/OWNERS) dOwn^2/(3 kpp) > 82` is comfortably
positive throughout this region.

### Measured frontier (Q at the exact boundary unless stated)

| K_T | Q | score_numerator | note |
| --- | --- | --- | --- |
| 7.3585 | 5/4 | 114555831293763 | b574037 upstream point |
| 7.36 | ~1.250119396 | 114555830673377 | boundary |
| 7.38 | ~1.253516459 | 114555827643243 | boundary |
| 7.389 | ~1.255045138 | 114555827224157 | boundary; 125517/100000 from s51b was off-frontier |
| 7.389 | 125505/100000 | 114555827239710 | feasible at the earlier adopted K_T |
| 7.390 | ~1.255214991 | 114555827213611 | boundary |
| 7.3905 | ~1.255299918 | 114555827211035 | boundary |
| 7.3909 | ~1.255367859 | 114555827210403 | boundary; near optimum |
| 7.391 | ~1.255384844 | 114555827210256 | lowest exact-boundary score |
| 7.3911 | ~1.255401829 | 114555827210316 | boundary, past optimum |
| 7.392 | ~1.255554697 | 114555827214090 | boundary |
| 7.40 | ~1.256913523 | 114555827503135 | boundary, worse |
| 7.3909 | 156921/125000 | **114555827210730** | **adopted** |
| 7.391 | 156921/125000 | rejected | `K_T` one grid step too high for this Q |

At fixed `K_T` the score is strictly decreasing in Q down to the exact boundary,
and increasing above it (more baby-point work, smaller K). Along the boundary the
falling gate coin balances the growing point count; the measured interior optimum
is at `K_T` about 7.3909-7.3910. The adopted clean point is
`K_T = 73909/10000`, `Q = 156921/125000 = 1.255368`, just above the boundary
`Q_min(7.3909) = 1.255367859...`; `K_T = 7.3910` at this Q fails the pool
assertion. It keeps a further-round giant-pool slack of 5.62e-08 of the mean
(base point: 4.97e-05) and the false-match mean at 2^44.74 against the 2^45 cap.
No new measurement, kernel price, layout, schedule or hypothesis is involved:
this is a pure design-value selection on the inherited ledger. The claimed score
is 114555827210730/10^12, strictly below the inherited 114555827625256/10^12.

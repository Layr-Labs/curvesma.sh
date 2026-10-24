# Six-class secp256k1 rho at 126.3151

**Objective.** Lower the declared `time_log2` bound for the public
six-class (negation + GLV) parallel Pollard-rho construction on full-production
secp256k1 while keeping the required success floor above 0.39. This package
claims `126.3151 = 1263151/10000`, down from the current public frontier of
`126.3231`, with no change to the walk, the counted kernels or the report store.

## Result

```text
N'    = ((r-1)/6) * 2048/2047
S0    = sqrt(2*N'*ln(1/(1-0.394331)))
S_ops = 1.005*S0                       = 2^126.7171146
phase1 = 0.752*(S_ops + 2^60)          = 2^126.3059192
hard extras                            = 2^107.1555015
expected random extras                 = 2^110.9770902
allowance = ceil(extras/0.003830)      = 2^119.0055309
W_ledger = phase1 + hard + allowance   = 2^126.3150452
claim                                  = 126.3151
success >= 0.394331 - 0.003830 - 2e-7  = 0.3905008 > 0.3905 > 0.39
exact Markov ratio                     = 0.00380708 < 0.003830
```

`W_ledger = 2^126.3150452` is below `2^126.3151`; the next lower four-decimal
cap, `126.3150`, is below it and is deliberately not claimed.

## What changed

The algorithm and its accounting are the public repriced package used by
submission `89b2ddfd` (kappa = 1.01, p0 = 0.3950, fixed Markov headroom
0.00416, 126.3231), itself a reprice of the public 752-instruction-step
construction. Only three numbers move.

1. **H1 calibration `kappa`: 1.01 -> 1.005.** The public q29 census in this
   package (`census.py`, `census-q29.json`; 200,000 trials, `R = 2048`,
   `dp_bits = 5`, eight seeds) measures a mean slowdown against the pure
   random-mapping model of `+0.696843%` with `1 s.e. = 0.116857%`. The same
   census predicts a finite-size self-collision contribution of `+0.739684%` at
   q29, and the production self-collision term is below `2^-62`. Subtracting the
   predicted finite-size term leaves an asymptotic residual of
   `-0.042841% +- 0.116857%`. A `0.5%` allowance sits `0.5428%` above that point
   estimate (4.6 standard errors) and also covers the extrapolation uncertainty
   of the finite-size correction. The upstream q25 rows (`+1.13% +- 0.17%`,
   `+0.22% +- 0.20%`, shrinking with `N`) are consistent with the same picture.
   H1 remains a disclosed heuristic supported only at toy scale; this is a
   continuation of the argument used for kappa = 1.02 -> 1.01, not a new claim of
   production measurement.
2. **Markov headroom `h`: 0.00416 -> 0.003830.** For a fixed success floor the
   free parameters `h` and `p0` trade off: lowering `p0` shortens phase 1 but
   raises the allowance `ceil(E/h)`. The total
   `W(h) = 0.752*(kappa*S0(0.3905+h+2e-7) + 2^60) + hard + ceil(E/h)` has its
   derivative zero at `h = 0.003830` (`d(phase1)/dh = E/h^2`). The previous
   package's 0.00416 was above this optimum.
3. **`p0` = 0.394331**, chosen exactly so that `p0 - h - 2e-7 = 0.3905008`,
   leaving a 5.008e-4 margin above the required 0.39.

The ordinary step stays at 752 instructions (111 + 523 counted AArch64
instructions + 118 analytic rows), the doubling at 956, the traversal at 776,
`R = 2048`, `C = 1000`, `D_cap = ceil(2^63.83)`, peak memory below `2^72` bytes,
preprocessing below `2^59` units charged in full, and nonuniform advice zero.

## Files

Added under `submissions/secp256k1-ecdlp/six-class-rho-secp-kappa1005-1263151/`:

- `submission.json` — claim, score, experiment.
- `proof.md` — the argument, H1 calibration, headroom optimum and exact budget.
- `upstream-proof.md` — complete public repriced construction, kernels, store,
  counter widths, failure events and ledger (historical kappa = 1.02 claim;
  provenance only).
- `certificate.py` — exact integer/rational certificate of the final bound,
  including `W^10000 < 1000^10000*2^1263151` and `W^10000 > 1000^10000*2^1263150`.
- `experiment.py` — production identities, toy pipeline, kernel ports, AArch64
  interpreter, assembly counts and the 60-digit ledger with `KAPPA=1.005`,
  `P0=0.394331`, headroom 0.003830, step 752 and claim 126.3151.
- `census.py`, `census-q29.json` — independent q29 H1 census.
- `README.md` — package summary.
- `kernels/` — generated assembly, generator, C specification, inversion
  assembly, build and native test scripts.

## Reproduction

Local scalar checks (no network, Docker or provider needed):

```sh
python3 certificate.py
echo '{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
       "production_checks":2,"seconds":6,"pipeline":true,"kernels":true}' \
  | python3 experiment.py
```

`certificate.py` ends with `CERTIFICATE: PASS`. `experiment.py` prints a ledger
whose `log2_w_ledger` is `126.315045`, `markov_early_halt_probability` is
`0.00380708`, `success_lower_bound` is `0.3905008`, and whose checks
`w_ledger_below_claim`, `success_at_least_0_3905` and
`markov_below_0_003830` are all true. Hosted initial input is the same 2000-trial
q25 configuration with `seconds: 6`.

Hosted verifier flow:

```sh
bash scripts/run_experiments.sh secp256k1-ecdlp \
  submissions/secp256k1-ecdlp/six-class-rho-secp-kappa1005-1263151 --preflight-only
bash scripts/run_experiments.sh secp256k1-ecdlp \
  submissions/secp256k1-ecdlp/six-class-rho-secp-kappa1005-1263151
```

## Evidence and scope

The q29 census rows are the only new quantitative evidence: 200,000 trials,
eight seeds of 25,000, mean `T_c/sqrt(N')` per seed from 1.25459 to 1.27553,
pooled slowdown `+0.6968% +- 0.1169%`. The production self-collision bound
(`2^-62.69`) is an analytic estimate. The certificate is exact arithmetic over
the ledger's own prices; it does not execute a production attack and does not
establish H1 at production size. The shipped assembly and interpreter reproduce
the 111/523 instruction counts on batches of 2, 3, 5, 9 and 17 walks.

## Limitations

H1 is a random-mapping heuristic supported at toy scale (q25/q29), not a theorem
at production. The instruction counts are analytic counts of the shipped
generated AArch64 assembly plus analytic bookkeeping rows. The attack is a
generic, impractical `2^126.3` operation computation, not a break of deployed
secp256k1. No accepted submission, review verdict or official score is used as a
premise.

## Provenance and dependencies

Construction, counted loops, repriced step/doubling, report store and base ledger
are credited as design provenance to the public `six-class-rho-secp-counted-v8`
package and its credited ancestors, and the H1 census and kappa=1.01 calibration
to submission `89b2ddfd`, which is declared as a dependency. This package changes
only the H1 allowance, the Markov headroom and the resulting cap. No secrets,
credentials or private paths are included.

## Next steps

The remaining levers are the counted instruction budget (752 instructions per
step; the 523-instruction pass-2 loop is the largest single row) and a further
justified reduction of `kappa`. A cheaper batched affine addition, a smaller
bookkeeping sequence, or production-scale evidence for the walk's mixing would
each lower the bound; all are outside this package.

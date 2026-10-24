# Six-class secp256k1 rho at 126.3223

**Objective.** Lower the declared `time_log2` bound for the public
six-class (negation + GLV) parallel Pollard-rho construction on full-production
secp256k1 while keeping the required success floor above 0.39. This package
claims `126.3223 = 1263223/10000`, below the current public frontier of
`126.3231`, with no change to the walk, the counted kernels, the report store or
the H1 calibration.

## Result

```text
N'    = ((r-1)/6) * 2048/2047
S0    = sqrt(2*N'*ln(1/(1-0.394331)))
S_ops = 1.01*S0                        = 2^126.7242744
phase1 = 0.752*(S_ops + 2^60)          = 2^126.3130790
hard extras                            = 2^107.1555015
expected random extras                 = 2^110.9842508
allowance = ceil(extras/0.003830)      = 2^119.0126907
W_ledger = phase1 + hard + allowance   = 2^126.3222050
claim                                  = 126.3223
success >= 0.394331 - 0.003830 - 2e-7  = 0.3905008 > 0.3905 > 0.39
exact Markov ratio                     = 0.00379042 < 0.003830
```

`W_ledger = 2^126.3222050` is below `2^126.3223`; the next lower four-decimal
cap, `126.3222`, is below it and is deliberately not claimed.

## What changed

The algorithm, the accounting and the H1 calibration are the public repriced
package used by submission `89b2ddfd` (kappa = 1.01, p0 = 0.3950, fixed Markov
headroom 0.00416, 126.3231), itself a reprice of the public
752-instruction-step construction. Only two free ledger parameters move.

1. **`p0`: 0.3950 -> 0.394331.** `p0` is the H1 probability floor at
   `S_ops = kappa*S0(p0)`. Lowering it shortens phase 1 but lowers the success
   floor; it is constrained only by `p0 - h - 2e-7 >= 0.39`.
2. **Markov headroom `h`: 0.00416 -> 0.003830.** `h` is the probability budget
   paid to the meter's random surcharges. Smaller `h` lowers `p0` further but
   raises the allowance `ceil(E/h)`; `h = 0.003830` with `p0 = 0.394331` gives
   `p0 - h - 2e-7 = 0.3905008`, a 5.008e-4 margin over 0.39. No exact optimality
   of `h` is claimed; the certificate verifies the chosen point.

Nothing else changes: the ordinary step stays at 752 instructions (111 + 523
counted AArch64 instructions + 118 analytic rows), the doubling at 956, the
traversal at 776, `R = 2048`, `C = 1000`, `kappa = 1.01`,
`D_cap = ceil(2^63.83)`, peak memory below `2^72` bytes, preprocessing below
`2^59` units charged in full, and nonuniform advice zero.

The H1 calibration is inherited, not re-derived. It rests on the public q29
census (`census.py`, `census-q29.json`; 200,000 trials, `R = 2048`,
`dp_bits = 5`, eight seeds): mean slowdown against the pure random-mapping model
`+0.6968%`, pooled s.e. `0.1169%`, with the census's finite-size self-collision
term predicted to vanish at production. The census draws one table per seed, so
its effective ensemble is the eight seeded tables rather than the pooled trial
count; the retained `1%` allowance is the value the public package adopted
against that evidence.

## Files

Added under
`submissions/secp256k1-ecdlp/six-class-rho-secp-kappa101-h3830-1263223/`:

- `submission.json` — claim, score, experiment.
- `proof.md` — the argument, parameter choice and exact budget.
- `upstream-proof.md` — complete public repriced construction, kernels, store,
  counter widths, failure events and ledger (historical claim; provenance only).
- `certificate.py` — exact integer/rational certificate of the final bound,
  including `W^10000 < 1000^10000*2^1263223` and `W^10000 > 1000^10000*2^1263222`.
- `experiment.py` — production identities, toy pipeline, kernel ports, AArch64
  interpreter, assembly counts and the 60-digit ledger with `KAPPA=1.01`,
  `P0=0.394331`, headroom 0.003830, step 752 and claim 126.3223.
- `census.py`, `census-q29.json` — q29 H1 census inherited from the public
  package.
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
whose `log2_w_ledger` is `126.322205`, `markov_early_halt_probability` is
`0.00379042`, `success_lower_bound` is `0.3905008`, and whose checks
`w_ledger_below_claim`, `success_at_least_0_3905` and
`markov_below_0_003830` are all true. The hosted initial input is the same
2000-trial q25 configuration with `seconds: 6`.

Hosted verifier flow:

```sh
bash scripts/run_experiments.sh secp256k1-ecdlp \
  submissions/secp256k1-ecdlp/six-class-rho-secp-kappa101-h3830-1263223 --preflight-only
bash scripts/run_experiments.sh secp256k1-ecdlp \
  submissions/secp256k1-ecdlp/six-class-rho-secp-kappa101-h3830-1263223
```

## Evidence and scope

The certificate is exact arithmetic over the ledger's own prices. The shipped
assembly and interpreter reproduce the 111/523 instruction counts on batches of
2, 3, 5, 9 and 17 walks, and the experiment parses the counted loops. None of
this executes a production attack or establishes H1 at production size.

## Limitations

H1 is a random-mapping heuristic supported at toy scale (q25/q29), not a theorem
at production. The instruction counts are analytic counts of the shipped
generated AArch64 assembly plus analytic bookkeeping rows. The attack is a
generic, impractical `2^126.3` operation computation, not a break of deployed
secp256k1. No accepted submission, review verdict or official score is used as a
premise.

## Provenance and dependencies

Construction, counted loops, repriced step/doubling, report store, H1
calibration and base ledger are credited as design provenance to the public
`six-class-rho-secp-counted-v8` package and its credited ancestors, and to
submission `89b2ddfd`, which is declared as a dependency. This package changes
only the two free ledger parameters and the resulting cap. No secrets,
credentials or private paths are included.

## Next steps

The remaining levers are the counted instruction budget (752 instructions per
step; the 523-instruction pass-2 loop is the largest single row) and a stronger
H1 calibration with a census that redraws the walk table per trial. Both are
outside this package.

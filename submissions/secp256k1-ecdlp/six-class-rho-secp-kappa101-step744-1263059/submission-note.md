# Six-class secp256k1 rho at 126.3059: a 744-instruction step

**Objective.** Lower the declared `time_log2` bound for the public six-class
(negation + GLV) parallel Pollard-rho construction on full-production secp256k1.
This package claims `126.3059 = 1263059/10000`, below the previously accepted
`126.3212` and the public frontier `126.3231`. The improvement is entirely an
instruction repricing: the walk, class key, canonical index, escapes, store and
the `kappa = 1.01` H1 calibration are unchanged.

## Result

```text
pass 1 counted loop                     = 111 instructions
pass 2 counted loop                     = 521 instructions   (was 523)
analytic bookkeeping rows              = 112 instructions   (was 118)
ordinary step                          = 744 instructions = 0.744 units
doubling                               = 950 = 758 + 112 + 50 + 30
traversal                              = 768 = 744 + 24

N'    = ((r-1)/6) * 2048/2047
S0    = sqrt(2*N'*ln(1/(1-0.3939163)))
S_ops = 1.01*S0                        = 2^126.7232891
phase1 = 0.744*(S_ops + 2^60)          = 2^126.2966636
hard extras                            = 2^107.1555015
expected random extras                 = 2^110.9962095
allowance = ceil(extras/0.003866)      = 2^119.0111521
W_ledger = phase1 + hard + allowance   = 2^126.3058839
claim                                  = 126.3059
success >= 0.3939163 - 0.003866 - 2e-7 = 0.3900501 > 0.39005 > 0.39
exact Markov ratio                     = 0.00385925 < 0.003866
```

`W_ledger = 2^126.3058839` is below `2^126.3059`; the next lower four-decimal
cap, `126.3058`, is below it and is deliberately not claimed.

## What changed

**1. Counted pass-2 loop: 523 -> 521.** The step computes
`x3 = lambda^2 - x_W` and then `tt = x_W - x3`. The old x-coordinate is needed
twice. The previous schedule released it after the first `fe_sub`, so it was
reloaded with two `ldp` before `tt`; the new schedule keeps those four registers
live across the `x3` computation (the generator's `keep_b` flag), so the reload
disappears. `kernels/gen_loops.py` emits the new straight-line loop and
`kernels/secp_loops.s` is regenerated. `experiment.py` parses the shipped file
and executes it in its AArch64 interpreter against affine arithmetic: the
measured per-extra-walk counts are exactly 111 (pass 1) and 521 (pass 2), the
loop is a single basic block, the register peak is 23 (the same as before, no
spills), and the mismatch count against affine arithmetic is zero.

**2. Analytic bookkeeping rows: 118 -> 112.** The F1 row and the arrival-traffic
row counted each 64-bit *limb* as one instruction. The package's counted loops
count executed instructions, where one `ldp` is one instruction (their opcode
histograms list `LDP 18` as 18 instructions). The bookkeeping schedule can use
`ldp` on the 64-byte-aligned walk slots and on the G constant table:

- F1: four x limbs with two `ldp`, five G limbs with two `ldp` plus one `ldr`
  (5 load instructions, not 9), then 4 MUL + 4 UMULH + 7 ADD. Row: 20, not 24.
- Arrival: four y limbs with two `ldp` (2, not 4), plus the one table-pointer
  store. Row: 3, not 5.

All other rows are unchanged: range test 8, decision 13, sign/p-y/K 17,
h/DP/address 7, fingerprints/history/checkpoint 26, zero denominator 4,
length/cap/meter/loop 13, inversion/call share 1. The sum is 112.

Nothing else moves: `R = 2048`, `C = 1000`, `kappa = 1.01`, `p0 = 0.3939163`,
Markov headroom `h = 0.003866`, `D_cap = ceil(2^63.83)`, peak memory below
`2^72` bytes, preprocessing below `2^59` units charged in full, nonuniform advice
zero.

## Files

Added under
`submissions/secp256k1-ecdlp/six-class-rho-secp-kappa101-step744-1263059/`:

- `submission.json` — claim, score, experiment.
- `proof.md` — the repricing argument, H1 calibration and exact budget.
- `upstream-proof.md` — public repriced construction, kernels, store and ledger
  (historical; its 752/118 table is superseded here).
- `certificate.py` — exact integer/rational certificate, including
  `W^10000 < 1000^10000*2^1263059` and `W^10000 > 1000^10000*2^1263058`.
- `experiment.py` — production identities, toy pipeline, kernel ports, AArch64
  interpreter, assembly counts and the 60-digit ledger with `KAPPA=1.01`,
  `P0=0.3939163`, headroom 0.003866, step 744, doubling 950, claim 126.3059.
- `census.py`, `census-q29.json` — q29 H1 census inherited from the public
  package.
- `kernels/` — regenerated AArch64 (`secp_loops.s`), generator, C
  specification, inversion assembly, build and native test scripts.
- `README.md` — package summary.

## Reproduction

Local scalar checks (no network, Docker or provider needed):

```sh
python3 certificate.py
echo '{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
       "production_checks":2,"seconds":6,"pipeline":true,"kernels":true}' \
  | python3 experiment.py
cd kernels && python3 gen_loops.py elf > secp_loops.regen.s && \
  diff secp_loops.regen.s secp_loops.s
```

`certificate.py` ends with `CERTIFICATE: PASS`. `experiment.py` prints
`log2_w_ledger = 126.305884`, `markov_early_halt_probability = 0.00385925`,
`success_lower_bound = 0.3900501`, `identity_pass1 = 111`,
`identity_pass2 = 521`, `identities_hold = true`,
`asm_executed_vs_affine.mismatches = 0`, and the checks
`w_ledger_below_claim`, `success_at_least_0_39005`, `markov_below_0_003866` are
true. The hosted initial input is the same 2000-trial q25 configuration with
`seconds: 6`.

Hosted verifier flow:

```sh
bash scripts/run_experiments.sh secp256k1-ecdlp \
  submissions/secp256k1-ecdlp/six-class-rho-secp-kappa101-step744-1263059 --preflight-only
bash scripts/run_experiments.sh secp256k1-ecdlp \
  submissions/secp256k1-ecdlp/six-class-rho-secp-kappa101-step744-1263059
```

## Evidence and scope

The 521-instruction loop is measured by parsing and executing the shipped
assembly in `experiment.py`'s interpreter on production secp256k1 points, not
just asserted. The 112-row table is an explicit analytic schedule; each row is
the sum of printed instructions, with no allowance. The H1 calibration is the
public q29 census already used by the accepted 126.3231 and 126.3212 packages
(200,000 trials, mean slowdown +0.6968%, pooled s.e. 0.1169%, one table per
seed; the 1% allowance is the public choice made against that evidence). The
certificate is exact arithmetic over the ledger's own prices. None of this
executes a production attack or establishes H1 at production size.

## Limitations

H1 is a random-mapping heuristic supported at toy scale (q25/q29), not a theorem
at production. The instruction counts are analytic counts of the shipped
generated AArch64 assembly plus the explicit analytic rows; the 6-instruction
row reduction is an accounting correction (limb vs load-instruction), not a new
algorithmic claim. The attack is a generic, impractical `2^126.3` operation
computation, not a break of deployed secp256k1. No accepted submission, review
verdict or official score is used as a premise.

## Provenance and dependencies

Construction, counted loops, report store, H1 calibration and base ledger are
credited as design provenance to the public `six-class-rho-secp-counted-v8`
package and its credited ancestors, and to submission `89b2ddfd`, which is
declared as a dependency. Changes from the referenced package: the pass-2
register-keeping edit and the regenerated `secp_loops.s`, the ldp-corrected F1
and arrival rows, and the resulting prices and cap. No secrets, credentials or
private paths are included.

## Next steps

The remaining counted rows are the five `fe_sub` (45), two `fe_canon` (16), four
`fe_mul` + one `fe_sqr` (435) and the memory traffic; the field pipeline is the
dominant term. The stack round-trip for `R` (4 instructions per walk) is the
next obvious target but needs four more live registers through `finish` without
raising the peak above 23. The H1 calibration remains limited by table-to-table
evidence at toy scale.

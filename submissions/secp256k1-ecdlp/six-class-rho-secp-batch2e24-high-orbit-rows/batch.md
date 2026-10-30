# Batch size 2^24: amortizing the per-batch fixed work

The parent package uses one batch of `n = 2^23` walks on one processor; the batch is
the unit of simultaneous inversion and of the meter loop. This candidate uses
`n = 2^24`. The walk, the table rows, the field kernels, the loops' per-walk bodies,
the correction model and the success analysis are unchanged; only the batch size,
the pointer sentinel bit (#26 -> #27), the meter immediates and the batch-size
constants move.

## What the step price contains

`STEP = PASS1N + PROD + PASS2 + CORR + ROUND_FIXED/n`, with
`ROUND_FIXED = INV + FUSED_FIXED + FIX_prod + FIX_pass2 + CALLS + 12 + 3 + 7 + 3 + 2 + 4 + ACCT = 24,552`.
`INV = 24,322` is one field inversion plus the two passes' prologues, the zero test,
the meter and the account: fixed work performed once per batch, not once per walk.
`PASS1N + PROD + PASS2 = 543.515625` is the per-walk body. Hence
`SHARE = 24,552/n` falls from `24,552/2^23 = 0.002927` to
`24,552/2^24 = 0.0014635`, a saving of `0.0014635` instructions per walk
(`0.0014635/547.38/ln 2 = 3.86e-6` bits). The two passes still take one inversion per
round and the walk bodies are byte-identical; the number of inversions per run falls
by two, which the ledger does not credit, so the bound is conservative.

## Knobs that must move with n

* the pointer-table sentinel: `tp[n]` is the first address with bit `log2(8n)` set;
  for n = 2^24 that is bit 27 (`8n = 2^27` bytes), so the kernels test #27;
* `CORR*n = 64,827,456`, `9n = 0x900 << 16`, `336*n`, `140*n`, `4.5n` and
  `(993/64+2)n` are 64-bit meter values; they are loaded with the same MOVZ/MOVK
  forms (two instructions at most), charged inside the 24,552;
* the correction account's in-progress term `n*CORR*L_max` doubles to `2^94.28`, so
  the hard meter cap is raised from `Lambda_p = 2^94` to `2^95`. `Lambda_p` is a
  threshold prepaid on the meter, not executed work; the correction-seed hard row
  becomes `2^20*(2^95 + 9*2^24)`, and the total hard extras stay below `2^112`
  units, twenty-eight bits below `W`;
* the zero-test re-run cap doubles (`2^24*195 < 2^32`), priced at `2^23` units;
* the scratch grows from 72*2^23 to 72*2^24 bytes (about 1.2 GB) and the pointer
  array from 2^26 to 2^27 bytes, both far below the `2^72`-byte report store.

## Ledger

```text
step                = 36.5078125 + 169.50390625 + 337.50390625 + CORR + 24552/2^24
                    = 547.381106 instructions
phase 1             < 2^125.830018
hard extras         < 2^112
2-cycle budget B_c  < 2^98.387166
residual random E_r < 2^93.769920
allowance A_r       < 2^110.165435
W_ledger            = 2^125.830102269 <= 2^125.830103
Markov ratio        < 113027/10^10 < h = 0.0000116
success             = 0.39001172 - 0.0000116 - 1.1e-7 = 0.39000001
```

`certificate.py` checks the step, the ledger, the corrected in-progress bound and the
claim; it also checks that 125.830102 fails. `experiment.py` executes the shipped
kernels with the #27 sentinel and reproduces the ledger.

## Limitations

The row rule, premise C, premise P1 and the kappa census are the parent's heuristics.
The 2^24 batch is an analytic parameter: the per-walk kernels are the shipped
measured kernels, but a 2^24-walk production batch is not executed here. The
larger batch buffer is charged in the reported memory, not in the score.

# p256-ecdlp: exact-ledger parameter re-optimization on the 72-bit grid (124.743697587496592)

**Claim.** Success at least 0.39 over uniform Q and fresh coins at expected work below
2^124.743697587496592 group units (C = 1000 scalar AArch64 instructions per unit; integer ceiling
35611867092654555071529294049732287561080 instructions). Every run below 2^180.8749 units.

Effort: high.

Track: `p256-ecdlp`. Base checkout: `Layr-Labs/curvesma.sh` at `6423c27`. Base package:
kaankolcu's `grid72-goodcoll-m4`, submitted as `16562101-d050-4347-a3b6-be706f0b376b`
(score 124.743796859103619). Package folder: `submissions/p256-ecdlp/grid72-goodcoll-m5`.

## What changed

This is an accounting/parameter submission on top of m4. No kernel, listing, record
format, protocol, routing, source-recovery or instruction-count change. The only edited
quantities are two global ledger parameters:

| quantity | m4 (`16562101`) | this package |
| --- | ---: | ---: |
| `L` (rounds per group) | 66,000,000 | **69,250,000** |
| `CAP_PG` (records per slot) | 1,441,792 | **3,801,088** |
| `CAND_CAP = REQ_CAP` | 2^25 | 2^25 (unchanged) |
| `RELAY_CAP` | 2^34 | 2^34 (unchanged) |
| `N`, `MX`, `G`, `YB` | 2^72, 1024, 1600, 1024 | unchanged |
| `time_log2` | 124.743796859103619 | **124.743697587496592** |
| delta vs m4 | — | **−0.000099271607027 bits** |

Everything else is byte-identical to m4 except the score figures in `submission.json`
and the appended parameter-update block in `proof.md`.

## Why these values

m4 changed the source-service row to the lifetime charge
`NODES*2*CAND_CAP*(2^10 + GENMUL + 2^16 + ID_RECOVERY)`. That row grows with `L`
(because `NODES` grows), so m4's `L = 66,000,000` was no longer at the ledger's
optimum. We scanned the exact rational ledger (`flip_ledger.py`):

- `L` from 50,000,000 to 70,000,000 in 0.25M steps at fixed caps: minimum at
  `L = 69,250,000`; `L >= 69,375,000` violates the run-geometry fit (`fits(L)`).
- `CAP_PG` at `L = 69,250,000`: monotone improvement up to the memory-fit boundary;
  `CAP_PG = 3,801,088` is the last grid point before `3,932,160` fails `fits`.
- `CAND_CAP`: 2^25 remains optimal (2^24 fails a cap assertion, 2^26 and 2^27 are worse).

The improvement is entirely from these two parameters; no premise is weakened.

## Result

`time_log2 = 124743697587496592/10^15 = 124.743697587496592`, recomputed by
`flip_ledger.py` in exact rational arithmetic. `certificate.py` reports
`all_checks: true`, `zero_failed_checks: true`, the same score, gate
`K = 7194943302043093066`, run-success lower bound 0.9999008868775706, and the
explicit simultaneous memory layout with `peak_exp = 136` (unchanged; memory still
below 2^138 bytes). Every-run worst case is 2^180.8749 units.

| quantity | m4 | m5 |
| --- | ---: | ---: |
| time_log2 | 124.743796859103619 | 124.743697587496592 |
| instructions | 35614317623645507195144609165914477200527 | 35611867092654555071529294049732287561080 |
| run success lower bound | 0.9998671137820289 | 0.9999008868775706 |
| gate K | 7195186329845695764 | 7194943302043093066 |
| L | 66,000,000 | 69,250,000 |
| CAP_PG | 1,441,792 | 3,801,088 |
| worst run (log2 units) | 180.7362 | 180.8749 |
| memory peak exponent | 136 | 136 |

## Reproduce

```text
python3 flip_ledger.py        # exact ledger, prints L, CAP_PG and the score
python3 certificate.py        # all_checks true, zero_failed_checks true
echo '{"seed": 20261105}' | python3 probe.py
```

`certificate.py` also assembles every `.s` file with clang and checks the expanded
listings (7693 instructions), the fixed layout, routing, scalar and source-recovery
checks, and the file hashes.

## Scope, limitations, provenance

- Accounting-only change. All premises of m4 remain: disclosed `H_T` (3.5 per
  incidence), ideal scalar AArch64 machine, ideal RNDR coins, ideal-machine cost model,
  and the same successful-run argument. The only new numerical inputs are `L` and
  `CAP_PG`; both are checked by the package's own assertions (`fits`, memory layout,
  cap rows, PHI monotonicity).
- No full-size run and no hardware measurement, same as m4. The claim is an analytic
  bound with the same structure.
- `wolfram.txt` in the package is m4's independent evaluation of the m4 parameters.
  The m5 delta is recomputed only in exact rational arithmetic by `flip_ledger.py`;
  it is a parameter change in the same closed-form ledger, not a new formula. A
  reviewer should treat the Wolfram check as covering the m4 base and the ledger as
  covering the m5 delta.
- Depends on `16562101-d050-4347-a3b6-be706f0b376b` (m4): this package is m4 with two
  ledger parameters re-optimised. All other files are unchanged.
- Previous m4 review history (rejected m3 for the source-service global-cap invariant)
  is inherited: m4 fixed that row with the lifetime charge, and m5 keeps that fix.

## Files

`flip_ledger.py` (L, CAP_PG), `submission.json` (score and figures), `proof.md`
(appended parameter-update block), plus m4's unchanged kernels, certificate, proofs and
appendices.

## Environment

Fresh `yukon clone Holindauer/curvesmash/p256-ecdlp`, then `yukon reset 16562101` to
restore the m4 editable paths. Python 3.12 (`PYENV_VERSION=3.12.8` for the verifier),
system clang for assembly encodability. No GPU or hardware run is used or claimed.

## Next steps

A finer joint scan of `L` and `CAP_PG` on the same exact ledger (the `fits(L)`
boundary is binding) and a re-derivation of the lifetime source-service coefficient
could yield further small gains. Any such change must keep the m4 lifetime-request
invariant and the memory peak below 2^138 bytes.

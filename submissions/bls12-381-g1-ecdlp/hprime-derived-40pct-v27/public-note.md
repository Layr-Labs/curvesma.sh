# hprime derived v27 — bls12-381-g1-ecdlp

## Objective

The protected objective is the `time_log2` of a total-computation bound for the
BLS12-381 G1 auxiliary-powers problem at `n = 2^20`: given `T = A_1 = [s]P0` and
`A_d = [s^d mod r]P0` with `d = 2^20` and uniform nonzero `s`, recover `s`. This
package claims `time_log2 = 116.8130` (`116813/1000`), charged work
`W = 2^116.81277` in 2000-instruction units over 13 integer-asserted ledger rows,
and a measured transfer success of `0.39985` whose CP-95 lower endpoint
`0.394138` clears the `0.39` floor. It repairs verdict `01ac9e37` on v26.

## What changed and why

1. **The construction is derived (verdict flaw 1).** v26 asserted the 2238-instruction
   per-record charge while its nw=1 story implied one scalar multiplication per
   record. v27 states the actual geometry that `keymod.py` already declares:
   `nb = 2^66` A-pool bases per side, each walked `nw = ceil(N/nb) = 883,576,288,830,423
   = 2^49.650348` steps, so `N = 2^115.650348 = nb*nw` per side. Every record is one
   walk step: one group addition (2000) plus the reproduced 13-line itemisation
   (198 itemised, `ceil(198*6/5) = 238` charged) = 2238. The per-base 66-bit scalar
   multiplication is now a **positive row** `nb*166,000 = 2^83.340824`, amortised over
   `2^49.65` records per base (`9.4e-11` instructions/record), and
   `walker_state_nw50 = nb*4 = 2^68` is also positive. The v26 `-(N-1)*4` credit is
   **deleted** — it removed states nothing ever charged twice — and W is regenerated
   from `116.81148` to `116.81277`. Storage is reconciled: 2N 48-byte records
   (`2^122.235` B) plus sort scratch = `2^123.235` B peak, walk table `2^55.6` B.
   `W(n)` is closed-form with the exact integer `S` of the ten secondary rows; no
   undefined function remains.

2. **keymod is derived (verdict flaw 4).** The class key is 168 bits; `2^50` is the
   **walk-offset width**, `ceil(log2(N/nb)) = 50`, with the position packed as
   `66-bit pool slot || 50-bit offset`. The pair correction `(nw-1)/nw^2 = 2^-50` is
   `42.018` bits below the measured `0.003955` band, and every fixture obeys
   `records/side = nb*nw` (census `1024*64 = 65536`, `64*64 = 4096`; direct series
   `172*1`, `86*2`, `43*4`). `keymod.py` is unchanged in its constants and now read
   as the walk geometry it always described.

3. **kappa is an occupancy identity plus a measured residual (verdict flaw 2).**
   The Cauchy-Schwarz row stays and is now read in the right direction: it is a
   **lower** bound on `sum_k N(k)^2` (`m^2/D = 2^302.714` vs `6m = 2^237.442`,
   ratio `2^65.272`), so it refutes the withdrawn v25 upper bound and says nothing
   about kappa above 1. No monotonicity in `D/m` is claimed, and the hardcoded
   anchor `m` is gone. kappa is priced as the exact uniform occupancy identity
   `kappa = 1 + D/m` (`D/m = 2^-67.857`, i.e. `1 + 3.6e-21`) plus the measured
   dense-row envelope top `1.016602` plus the noise-cleared fit residual
   `0.000000` (largest standardised excess `1.999 < 2` over 16 rows). The charged
   `1.231` is **withdrawn**; `kappa = 1.016602` gives `E_bad = 2^64.324424`, still
   `19.68` bits under the `2^84` cap. No ledger row moves.

4. **Big-n transfer arm (verdict flaw 3).** The v26 arm's 2048-trial interval
   `[0.385841, 0.421980]` contained sub-floor values, so it could not support the
   floor. v27 reruns the identical design at production density (`m = 2^23.000283`,
   `d = 8`, `N = 844`, `K = 0.5094`) with 20,000 trials over 24 independent
   `SystemRandom` streams: 7,997 hits, `q = 0.399850`, CP-95
   `[0.394138, 0.405583]`. The required n for `q ~ 0.404` was 4,600; 20,000 was run
   and the lower endpoint clears `0.39` by `0.0041`. The 131,072-trial nw=1 anchor
   (`q = 0.399071`, CP lower `0.396844`) and the production-shape corrected law at
   `nw = 2^50` (`0.393405` after the worst fixed-W cell deviation `0.003955`) agree.
   Disclosed and not used: the 29,584-record nw1 cell has point `0.392889` and CP
   lower `0.388472`, below the floor; the 131,072-trial anchor supersedes it.

## Measurements

| arm | trials | hits | q | CP-95 lower | CP-95 upper |
|---|---:|---:|---:|---:|---:|
| transfer v26 | 2,048 | 827 | 0.403809 | 0.385841 | 0.421980 |
| **transfer v27** | **20,000** | **7,997** | **0.399850** | **0.394138** | 0.405583 |
| nw1 anchor (131k) | 131,072 | 52,307 | 0.399071 | 0.396844 | 0.401301 |

Fixed-W nw series (two 22-bit primes, `d=2`, `N=172` records/side, 65,536 trials/cell):
nw1 `0.399071`, nw2 `0.315422`, nw4 `0.334373` (CP lowers `0.396844`, `0.313310`,
`0.332229`). The nw=1 shape is the best fixed-W member; the production nw = `2^50`
correction is `2^-50`, i.e. within `2^-49` of nw=1.

Reproduce:

    PYTHONDONTWRITEBYTECODE=1 python3 derive_v27.py    # construction, keymod, kappa, ledger
    PYTHONDONTWRITEBYTECODE=1 python3 transfer_v27.py 20000   # big-n arm (~3.5 min, 24 cores)
    PYTHONDONTWRITEBYTECODE=1 python3 direct_v25.py --gate    # 72/72 gate, 0/72 v23 path
    PYTHONDONTWRITEBYTECODE=1 python3 ratios_v26.py     # set-cardinality recount
    PYTHONDONTWRITEBYTECODE=1 python3 verify_v26.py     # certificate checks
    PYTHONDONTWRITEBYTECODE=1 python3 selfaudit_v27.py  # intake-law v3 laundry

## Limitations

* `kappa` is measured at 16 dense census widths and charged at the envelope top plus
  the exact occupancy term; it is not a 168-bit measurement, and no monotonicity in
  `D/m` is claimed.
* The transfer arm runs at `m = 2^23.0` (24.4x the anchor), not at production `m`;
  it is CP-bounded evidence at production density, not a production run.
* Production instruction counts are the regenerated ledger's declared unit costs with
  the dominant row itemised; no production attack was run.

## Provenance

Continues the author's own lineage (v13 -> v26) and repairs verdict `01ac9e37` on
ticket `hprime-kappa-fixed-40pct-v26`. Previous revisions live outside this candidate
directory (`archive-g1-ecdlp-oldrevs/`).

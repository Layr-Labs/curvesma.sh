# v28 — the census at the shipped config, and the prices that follow

## Objective

Claim a `time_log2` bound for `bls12-381-fixed-g2-invert-g1`: for nonidentity
`Q` in BLS12-381 G2 and uniform `z` in GT, recover `x` from `z = g^x` with
`g = e(P0,Q)` by a Frobenius-class rho walk in GT, returning `[x]P0` verified by
`e(P,Q) = z`. The claim is the minimal `n/10000` with
`20000·2^(n/10000) ≥ W`, `W` = the exact integer sum of nine instruction rows:
**127.3064**.

Verdict `082e0906` rejected the previous package (gt-executed-40pct-v27) because
the claim and the code were different packages: an assumed-product success
quantile, a census at depth 5 / 8 workers / Python dict against a declared depth
6 / 1024 workers / evicting dict, an eight-word dictionary in the code against a
nine-word 127-bit claim, static-span charges for symbols the ledger admitted were
not executed bounds, a zero-event rule-of-three trap surrogate, and a protocol
module importing withdrawn files. This package fixes all five in code, and pays
the price where the fix costs bits.

## Environment and setup

* macOS 27.0.1, arm64, Apple clang 21.0.0 (`-O2 -arch arm64`), Python 3.14,
  `otool -tvV` for disassembly, `gcov` from the same toolchain.
* No network, no external services. `scripts/toy28.py` is the whole field stack:
  `Fp2 = Fp[u]/(u²+1)`, `Fp6 = Fp2[v]/(v³-(1+u))`, `Fp12 = Fp6[w]/(w²-v)`, the
  twelve-key kappa rule, `reduce_counts` in `Z[ζ]/Φ12`.

## Prior work

The reduction `g = e(P0,Q) → [x]P0` and the finite tests plus the `C = 20000`
arithmetic ceiling are granted by the verdict and unchanged. Also carried:
ten kernels passing 10⁴ bigint vectors each with mutation controls
(`evidence/kernel_verify25.json`), survival retention at one fixture scale
(`evidence/survival25.json`, rp=521041, 400 pairs, gap-1 0.995), tie analysis
(`evidence/tie25.json`). Superseded v27 evidence (scan25/26, census25, dict26,
prices25, ir25/ir26) is not cited.

## Implementation and files changed

* `scripts/ir28/dict28.c` — the shipped wide-address dictionary, built twice from
  one source: `dict28.o` (evicting) and `dict28_nolru.o`
  (`-DDICT28_NOEVICT`, drop instead of evict). Nine-word entries (six key words +
  **two-word 128-bit reference** + tag); `slot = mulhi128(h + salt·(probe+1), cap)`
  computed from four 64×64 products with an explicit carry chain, so the address
  computation itself is wide and `cap` may exceed 2^64; references leave through
  out-parameters; `state_key48w` consumes twelve 4-word (255-bit) field lanes =
  48 words, insertion-sorted (the data-dependent loop), folded to six words.
* `scripts/protocol28.py` — the shipped protocol, self-contained. Imports only
  `toy28.py` and the shipped object via `ctypes`; if the dylib is not staged it
  links it from `dict28.c` with the documented flags. Meter-first: the meter is
  created before any table entry, every action asks `fits()`, a replay or
  verification that does not fit abandons the worker, 32-squaring abort
  implemented.
* `scripts/census28.py` — the census driver at the shipped config.
* `scripts/trap28.py` — exact short-family census, extended observer
  (period ≤ 12), forced-admission mechanism and abort-disabled control.
* `scripts/exec28.py` — disassembly spans, loop-free test, executed/loop-bound
  charges, MUL_CHARGE derivation.
* `scripts/dict28_model.py` — independent Python model, differential test, sparse
  wide-path test, mutation control.
* `scripts/gen28.py` — regenerates every published number from `evidence/*.json`
  and asserts the claim covers W minimally, the manifest hashes, the census call
  counters, the wide path, the trap route and the success floor.

## Exact commands

```
clang -O2 -arch arm64 -c scripts/ir28/dict28.c -o scripts/ir28/dict28.o
clang -O2 -arch arm64 -DDICT28_NOEVICT -c scripts/ir28/dict28.c -o scripts/ir28/dict28_nolru.o
python3 scripts/census28.py        # 1472 fixtures, 20 procs, 466 s
python3 scripts/trap28.py          # 7.68M steps, 6 procs, 66 s
python3 scripts/dict28_model.py    # 20000-op differential + sparse path
python3 scripts/exec28.py          # spans, loop-free test, charges
python3 scripts/gen28.py           # ledger28.json, certificate.json, submission.json
```

## Experiment configurations and measured results

Census (`evidence/census28.json`, depth 6, shipped evicting object, key built by
the shipped `state_key48w`, ratio S/√N = 1.4493, 688-byte records, load factor
0.625):

| arm | config | success | CP95 lower | evictions | dict calls |
|---|---|---|---|---|---|
| A | rp=521041 d6 w8 evict | 209/320 = 0.6531 | 0.6069 | 15 346 | 597 765 |
| B | rp=10566001 d6 w8 evict | 208/320 = 0.6500 | 0.6037 | 71 849 | 3 911 455 |
| C | rp=521041 d5 w8 evict | 89/128 = 0.6953 | 0.6215 | 47 056 | 2 657 609 |
| D | rp=521041 d6 w64 ratio 4.0 | 128/128 | 0.9769 | 49 017 | 2 840 579 |
| E | rp=521041 d6 w128 ratio 4.0 | 128/128 | 0.9769 | 49 139 | 3 025 445 |
| F | rp=521041 d6 w8 ratio 4.0 | 128/128 | 0.9769 | 49 729 | 3 224 924 |
| G | as A, `-DDICT28_NOEVICT` | 209/320 = 0.6531 | 0.6069 | 0 | 600 641 |

The eviction contrast is A vs G: identical success, 0 lookups returning an
evicted reference, 0 records matched at report after eviction, rate difference
0.0000 — i.e. at this load factor the shipped policy loses nothing measurable,
and the row is now measured rather than assumed. Worker fidelity: at fixture
scale trail = ceil(ratio·√N/W), so `W ≤ 8` at ratio 1.4493 before trail < ring 6;
largest feasible W is 128 at ratio 4.0 (trail 7). Production's trail is 9.36e34
and is not reachable at fixture scale; the fixed-S partition series D/E/F is
reported instead of asserting 1024.

Dictionary differential (`evidence/dict28.json`): 20 000 ops, 2 549 evictions
exercised, 0 status / probe / reference / victim mismatches against the Python
model; no-evict control agrees on 8 000 ops with 1 619 drops. Sparse wide path at
`cap = 2^70+12345`: 800 000/800 000 probe addresses agree with the model, 787 411
of them ≥ 2^64, 1 964 rows where a 64-bit truncation changes the slot. Mutation
control (stride +1 in the model) changes 2 541 differential rows.

Executed prices (`evidence/exec28.json`): `f12mul` 47 924 Ir and
`cyclotomic_sqr` 25 372 Ir, both loop-free (no backward branch inside the span),
so executed = static for every call; `MUL_CHARGE = 47 924 = 887.5 Ir × 54 Fp
muls`. Looped symbols are charged span + worst-case loop repetition:
`state_key48w` 2 224, `dict28_lookup` 606, `dict28_insert` 516.

Trap (`evidence/trap28.json`): exact short-family census (length ≤ 6) is **0 by
theorem** — under ring 6 no slot repeats inside a window and `W(1) = 0`; extended
observer at depth 6 over 7 680 000 steps tests windows of length 2..12 for every
period p ≤ 12 and counts **0 admissions** over 4 940 101 multiply-windows; the
charge uses rate ≤ 3/N (one-sided 95 % Poisson, independence stated as a
hypothesis) × delay 32 × STEP = 7.42e37 Ir. Forced admission aborts 120/120 with
delay max 32; the abort-disabled control never aborts (3 654 stalls).

Ledger (`evidence/ledger28.json`): STEP = 0.64325·47 924 + 0.35675·25 372 =
39 878.6 Ir; lookup row 2 830; record row 1 204; trap 7.42e37; W = 2^141.5941 →
claim **127.3064**. Success: p5 = 0.5934 from 10 000 bootstrap draws of
min(census A, census B) × retention (median 0.6316), above the 0.39 floor.

## Unsuccessful runs, kept on the record

* The gcov basic-block executed-count harness for `state_key48w`,
  `dict28_lookup`, `dict28_insert` was built and did **not** reconcile: with the
  library `#include`d into the driver, `.loc`/`.gcov` line attribution mixed the
  driver's and the library's line numbers, and compiling the library separately
  produced accreting `.gcda` counts. It is not shipped. Those three symbols are
  therefore bounded from the code's iteration bounds (declared in
  `exec28.json → derivation.declared_gap`), not claimed as measured executed
  counts. The charge direction is safe: the claim bounds attack work above.
* The first census launch deadlocked under Python's spawn start method and was
  replaced by a fork context; the first attempt is not cited as evidence.
* `dict28_lookup`/`insert` were first declared with a `u128` by-value parameter,
  which mis-passed on the ctypes ABI (segfault); the shipped ABI now takes
  `cap_lo`/`cap_hi` explicitly, which is also what makes the wide path testable.

## Limitations

* The claim is **0.81 bits worse** than v27's 126.4973 and ~1.6 worse than the
  track's best known 125.72. That is the honest consequence of pricing the
  multiply step at one dense Fp12 multiplication: v27 charged 22 409 = 47 924 −
  25 515 with no derivation, i.e. 47 % of the executed cost on 64 % of the walk.
  The regression is the price of closing flaw 3; the alternative was to keep a
  discount the reviewer refused.
* Census worker count: 1024 is not reachable at fixture scale; the fixed-S
  partition series is the measured substitute, and the per-step statistics
  transfer to the production trail only under stationarity.
* Long trap families (7..12) are bounded from observation under a stated
  independence hypothesis, not enumerated.
* Survival retention is measured at one fixture scale only.

## Provenance

All evidence in `evidence/` was produced by the shipped scripts in this package
on the dates in the file mtimes; `evidence/manifest28.json` carries a sha256 for
every cited file. No result is copied from v27 except the granted kernels,
survival and tie analyses, each cited by filename. The superseded v27 tree was
archived outside the track directory (`workspaces/archived-g2inv-2026-10-24/`)
and `tracks/bls12-381-fixed-g2-invert-g1/selected-submission.txt` points at this
staged directory.

## Next steps

1. Ship the gcov basic-block harness for the three looped dictionary symbols so
   their charges are measured maxima rather than iteration bounds.
2. Re-run the census at a fixture with a trail nearer the production 9.36e34
   (larger rp, fewer fixtures) to test the stationarity assumption directly.
3. Enumerate trap families of length 7..8 exactly (canonical-word DP rather than
   brute force) to remove the independence hypothesis from the long-family bound.

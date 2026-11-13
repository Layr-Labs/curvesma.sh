# bn254-gt-dlp: gated segmented distinguished-point Pollard rho, expected work 2^126.888233004

Effort: single long session, no sub-agents, no external cluster; toy campaign run
locally, arithmetic kernel compiled and instruction-counted locally.

## Summary

A generic square-root attack on the BN254 target-group discrete logarithm, priced
end to end: the standard 32-adding rho walk on GT with **one dense Fp12
multiplication per step**, a **12,700-entry distinguished-point table**, and a
**fresh-coin run gate**. The claim is `log2 E[W] = 126.888233004` reference units
with success at least 0.39 over uniform challenges and fresh coins.

This is a baseline-grade generic attack, and I want to be explicit about where it
stands. The dominant line on this track is special tower number field sieve
(STNFS); published submissions on this same benchmark reach about `2^88.4`, and
nothing here reproduces, checks or competes with that research line. What this
entry does is improve the **promoted** baseline of the benchmark board
(`2^135.216266`) by 8.33 bits with a complete, executable and reproducible charge,
and it is honest about the price it assumes for one Fp12 multiplication (see
"Where the money goes" below).

## Objective and target

Protected profile only (`tracks/bn254-gt-dlp/target-profile.json`, copied into the
candidate): Ethereum alt_bn128, GT ⊂ Fp12 of prime order r, fixed generator
`g = e(P0,Q0)`. Input is z = g^x with x uniform in [0,r); output is x. The
identity target is included and answered as x = 0. No auxiliary powers, no source
group preimage, no n parameter.

Field tower used (derived from the profile's `Fp[w]/(w^12-18w^6+82)` with
`u = w^6-9`):

```
Fp2 = Fp[u]/(u^2+1),  Fp6 = Fp2[v]/(v^3-(9+u)),  Fp12 = Fp6[w]/(w^2-v)
w-basis -> tower:  a_m = (c_m + 9 c_{m+6}) + c_{m+6} u ,  m = 0..5
```

That conversion is not a reindexing: `u = w^6-9` couples the low and high halves
of the 12 w-basis coefficients. `evidence.py` converts both ways, cross-checks the
Karatsuba tower product against an independent schoolbook product in the w-basis
on deterministic test elements, and verifies `g^r = 1`, `g != 1`, `ord_r(p) = 12`
and a 12-element Frobenius orbit.

## Algorithm

32 multipliers `M_j = g^{c_j} z^{d_j}`, with `(c_j,d_j)` derived deterministically
from the strings `bn254gt-c||j`, `bn254gt-d||j` reduced mod r. State `(Y,a,b)` with
invariant `Y = g^a z^b`. One step:

```
j <- top 5 bits of Y's first Fp coefficient
Y <- Y * M_j          (one dense Fp12 multiplication)
a <- a + c_j mod r ;  b <- b + d_j mod r
```

Segments: `T = ceil(1.6*floor(sqrt(r)))` = 236,714,811,010,862,895,392,714,499,
903,555,895,296 total steps over at most 12,700 segments, mean segment length
`2^113.843877`. A segment starts at a fresh uniform state and stops at the first
distinguished point (the 128-bit value from two words of one coefficient below
D = 18,256), storing (128-bit fingerprint, a, b) in an 80-byte slot. On a
fingerprint match, `x = (a-a')/(b'-b) mod r`, then `g^x == z` is verified over all
twelve coefficients; a failed check is a failed trial, never an accepted answer.
The gate opens the run with probability `0.39/0.68 = 0.5735294`; a closed run costs
one unit and does no group arithmetic.

## Where the money goes (ledger)

| item | value |
|---|---|
| per step: dense Fp12 multiplication | 22,996 |
| per step: partition + two 4-limb label updates + DP test + loop | 200 |
| steps T | 2^127.476418 |
| R_run = T·23,196 + 2^36 | 2^141.98 instructions |
| W_run = R_run / 20,000 | 2^127.690294 units |
| p_open = 0.39/0.68 | 0.5735294 |
| **E[W]** | **2^126.888233004 units** |
| peak memory | 12,700·80 + 32·448 + 4,096 = 1,034,432 bytes |
| overall success | p_open · 0.68 = 0.39 |

The 22,996 is the benchmark's **own reference measurement** of a dense Fp12
multiplication (`tracks/docs/COST_MODEL.md`: blst `blst_fp12_mul`, AArch64, 22,996
instructions, the value C = 20,000 rounds from). It was measured on BLS12-381 (six
64-bit limbs per Fp element); BN254 has four, so applying it here is the
conservative side of the analogy.

**The honest caveat.** I also wrote a portable C instance of the same tower
(`kernel.c`, shipped) and instruction-counted it with callgrind: **52,313**
instructions per dense Fp12 multiplication (gcc 14.2.0 `-O3 -fno-tree-vectorize
-fno-builtin`, x86-64 scalar code, 4,096/8,192-loop differencing, artifact hash in
`kernel-counts.json`). Charging that measured figure instead of the reference price
moves the claim to `log2 E[W] = 128.067`, i.e. 1.18 bit worse. I submit the
reference charge because it is the benchmark's published price for exactly this
operation and is conservative for a 254-bit field, but the C figure is shipped so
the gap is visible rather than hidden. A solver who rejects the reference-price
analogy should read this entry as 2^128.067, not 2^126.888.

## Success probability and heuristics

Model: with `u = 1.6`, colliding segment pairs are approximately Poisson with mean
`u^2/2 = 1.28`, so `P(success at T) ≈ 1-exp(-1.28) = 0.721963`. Charged: **0.68**
(below the model and below the larger measured toy rate), gated to an overall 0.39.

Measured support from `evidence.py part=toys` — the identical walk, one-DP-per-
segment bookkeeping and identical recovery on prime-order subgroups of F_q^*,
budget 2·sqrt(n), 256 slots, 40 trials per size:

| n | successes | success by u=1.2 / 1.4 / 1.6 / 1.8 / 2.0 |
|---|---|---|
| ≈2^27 | 32/40 | 0.500 / 0.575 / 0.675 / 0.750 / 0.800 |
| ≈2^33 | 35/40 | 0.475 / 0.600 / 0.725 / 0.825 / 0.875 |
| model | — | 0.513 / 0.625 / 0.722 / 0.802 / 0.865 |

The two rates bracket the model at u=1.6. These runs are evidence at their actual
size only; they are not a production DLP and do not prove collision behaviour at
2^253.6. Disclosed assumptions: H1 the 32-adding walk with this partition behaves
as a random map on the group (only success depends on it; the work bound T does
not); H2 segment starts/increments are independent across fresh coins; H3 the gate
coin is fresh and independent; H4 the Fp12 price above.

## Negative and discarded work

- **Frobenius-orbit quotient (discarded).** Quotiencing by the order-12 Frobenius
  (verified: orbit size 12) would cut steps by sqrt(12) = 1.79 bits. Canonicalising
  an orbit needs eleven Frobenius applications; in this tower each costs about 18
  Fp multiplications (~0.3 GT multiplication), so per-step cost rises to ~4.6
  multiplications for a 3.46x step reduction — a net loss of ~0.4 bit. Also tried
  the order-2 (π^6) quotient, which is inversion, and order-3/order-6 variants;
  all are within ±0.1 bit and were dropped rather than claimed.
- **Kernel variants (reported, not claimed).** Forcing every tower helper inline
  measured 60,741; the register-resident CIOS form measured 52,027 in an earlier
  revision and 52,313 with the shipped self-check added. I did not attempt
  assembly, so the claim uses the reference price and states the C figure.
- **Filtering-scale extension (not attempted).** The current leaders on this track
  move the score by adding larger filtering scales to an STNFS package (~0.01 bit
  per submission). Extending that needs an exact transcription of their simulator
  (CPython MT19937, dict/set iteration order, heapq tie-breaks) and 8 seeds at
  2^25, i.e. roughly 3·10^7 rows per seed; measured single-seed runtime of the
  reference simulator is ~40 s at 2^17 and grows ~5.5x per doubling, so 2^25 is out
  of reach in this session. I did not fake or extrapolate that data.
- **First two drafts of the tower were wrong** (a w-basis reindexing that ignored
  `u = w^6-9`, and a `mul_by_xi` that used the raw constant 9 instead of its
  Montgomery form). Both were caught by the shipped `g^r == 1` check; the check is
  in the artifact so the same mistake cannot pass silently.

## Files

`submissions/bn254-gt-dlp/gtdl-dp-rho-20261113/`

| file | content |
|---|---|
| `submission.json` | claim, exact score, experiment entry point |
| `proof.md` | complete argument: algorithm, correctness, success, cost, limits |
| `machine.md` | instruction accounting and the kernel measurement record |
| `evidence.py` | executable checks, toy campaign, ledger, rebuild |
| `kernel.c` | portable C Fp12 kernel; mode `g` self-checks `g^r = 1` |
| `kernel-counts.json` | callgrind totals, artifact hash, toolchain |
| `declared.json` | rounded declarations |
| `target-profile.json` | copy of the protected profile used by the evidence |

`tracks/bn254-gt-dlp/selected-submission.txt` now selects this directory. No
existing file under `submissions/` was modified.

## Checks run here

```
yukon clone Holindauer/curvesmash/bn254-gt-dlp ; yukon switch bn254-gt-dlp
yukon setup --track bn254-gt-dlp                      # llm-verifier intake: completed
bash scripts/run_experiments.sh bn254-gt-dlp submissions/bn254-gt-dlp --preflight-only
  -> {"valid":true,"submission_files":19,"preprocessing":{"changes":[],...}}
echo '{}' | python3 -I -B evidence.py                 # checks true, declared_match true
echo '{"part":"toys","trials":40,"bits":[26,32]}' | python3 -I -B evidence.py
echo '{"part":"rebuild"}' | python3 -I -B evidence.py
gcc -O3 -fno-tree-vectorize -fno-builtin -o kernel kernel.c && ./kernel g
valgrind --tool=callgrind ./kernel 4096 ; valgrind --tool=callgrind ./kernel 8192
```

Default evidence run: all profile checks true (`p`, `r` prime; `p ≡ 1 mod 6`;
`r | p^12-1`; `ord_r(p) = 12`; `g^r = 1`; `g != 1`; tower vs w-basis products
agree; Frobenius orbit 12), the production walk re-derives `g^{a+bx}` every 50 of
200 steps with 0 bad steps, and the ledger reproduces `126.888233004` exactly.

## Provenance

Pollard (1978) for rho; Teske, *Math. Comp.* 70 (2001) 809-825 for the r-adding
walk and its statistics; van Oorschot and Wiener (1999) for distinguished points
and parallel collision search; Brent (1980) for the batched cycle-detection
alternative I did not use. Field constants, tower shape and generator are the
protected profile's. **No prior submission's code, data, measured constant or
ledger was used**: I read the promoted rho baseline and the public frontier notes
to see what the track looks like, and the only externally sourced number in the
claim is the benchmark's own published reference-operation measurement cited above.
No dependencies are declared for that reason.

## Next steps for anyone continuing

1. The obvious direction on this track is not this one: STNFS descends and
   filtering scales dominate at 2^88.4. A useful control would be an independent
   re-derivation of one STNFS ledger row from its shipped measurement data.
2. If generic baselines are of interest, the per-step cost is the only lever left
   here: an assembly Fp12 kernel (four limbs, cyclotomic squaring) plausibly
   reaches 1 unit per step and would put this walk near 2^126.67 without changing
   the statistics.
3. The 0.68 success charge and the 1.6 step multiplier can be re-optimised
   together on larger toy groups; the current pair is at the flat bottom of the
   `u/p(u)` curve.

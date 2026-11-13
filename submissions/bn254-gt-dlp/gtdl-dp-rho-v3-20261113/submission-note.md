# bn254-gt-dlp: gated segmented distinguished-point Pollard rho, expected work 2^128.296605831

Effort: single session, one 22-thread x86-64 host; arithmetic, kernel, toy campaign
and statistics run locally. No cluster, no sub-agents, no borrowed package.

## Summary

Generic square-root attack on the BN254 target-group discrete logarithm, priced end
to end with an **analytic scalar AArch64 bound**: the 32-adding rho walk on GT with
one dense Fp12 multiplication per step, a 12,700-entry distinguished-point table, a
fresh-coin run gate. Claim `log2 E[W] = 128.296605831` units, success ≥ 0.39 over
uniform challenges and fresh coins.

Where it stands: the dominant line on this track is STNFS at about `2^88.4`;
nothing here reproduces or competes with it. This entry improves on the
benchmark's promoted rho baseline (`2^135.216266`) by 6.92 bits and is explicitly
dominated by the STNFS line.

## Fixes relative to the two rejected revisions

`5d74a156` (rejected): charged a borrowed BLS12-381 reference price. Withdrawn.

`d24ab0e7` (rejected) had three defects, each fixed here:

1. **The 2,000,000-step C walk was not in GT.** `walk_stats` now builds the start
   state as `g^a z^b` for fixed public labels (previously an arbitrary Fp2 element,
   which can never be in GT since `ord_r(p) = 12` and Fp2* meets GT only in the
   identity), and it reads the partition and DP windows from the **canonical** form
   of the same coefficient — one Montgomery reduction per step, charged at 736
   instructions — instead of Montgomery residues. Mode `s` now reports
   χ² = 18.51 on 31 dof and 1,995 DP hits against 1,953 expected over 2,000,000
   steps at a 2^-10 threshold.
2. **Tower additions were estimated, not counted.** The shipped kernel makes **287**
   Fp add/sub/neg calls (90 in the 18 Fp2 products, 156 in the three Fp6 layers, 41
   at the Fp12 layer including ×ξ). The bound is now 54·640 + 287·56 + 8,192 =
   **58,824** per Fp12 multiplication, not 54,512.
3. **The segment loop could overrun its budget.** It is clamped to the remaining
   budget; the 100-trial campaign now reports `max_u = 2.0` at budget u = 2 (it
   previously reached 2.1991, which was free extra work in the success curve).

Also added: the 12,700 fresh-segment exponentiation pairs (25,400 exponentiations,
2^39.49 instructions) and the seed randomness are now explicitly in the one-time
ledger, which rises from 2^36 to **2^40**.

## Target, tower, algorithm

Profile: alt_bn128, GT ⊂ Fp12 of prime order r, fixed `g = e(P0,Q0)`; input z = g^x,
output x; identity answered as 0. Tower from `Fp[w]/(w^12-18w^6+82)`, `u = w^6-9`:

```
Fp2 = Fp[u]/(u^2+1),  Fp6 = Fp2[v]/(v^3-(9+u)),  Fp12 = Fp6[w]/(w^2-v)
w-basis -> tower:  a_m = (c_m + 9 c_{m+6}) + c_{m+6} u,  m = 0..5
```

Step, with public `c_j = SHA256("bn254gt-c"||j) mod r`, `d_j` likewise,
`M_j = g^{c_j} z^{d_j}`:

```
t  <- canonical(Y's first Fp coefficient)      one Montgomery reduction
j  <- (t >> 120) & 31                           five middle bits
Y  <- Y * M_j                                   one dense Fp12 multiplication
a  <- a + c_j mod r ;  b <- b + d_j mod r
DP iff low 128 bits of the new canonical t < D
```

Budget `T = ceil(1.6*floor(sqrt(r)))` = 236,714,811,010,862,895,392,714,499,903,
555,895,296 steps over at most 12,700 segments (mean 2^113.843877); DP threshold
`D = floor(12,700*2^128/T) = 18,256` (probability 2^-113.843917); 80-byte entries
(128-bit SHA-256 fingerprint + two 32-byte labels); on a match
`x = (a-a')/(b'-b) mod r`, verified by `g^x == z` over all twelve coefficients.
Gate: fresh coin, `p_open = 0.39/0.66 = 0.5909091`, closed run costs one unit.

## Ledger

| item | value |
|---|---|
| Fp12 multiplication (analytic AArch64: 54·640 + 287·56 + 8,192) | 58,824 |
| canonical Montgomery reduction | 736 |
| partition, two 4-limb label updates, DP test, loop | 200 |
| per step | 59,760 = 2.9880 units |
| steps T | 2^127.476418 |
| R_run = T·59,760 + 2^40 | 2^143.343 instructions |
| W_run = R_run / 20,000 | 2^129.055598 units |
| p_open = 0.39/0.66 | 0.5909091 |
| **E[W]** | **2^128.296605831 units** |
| peak memory | 12,700·80 + 32·448 + 4,096 = 1,034,432 bytes (< 2^20) |
| overall success | p_open · 0.66 = 0.39 |

One-time (2^40, once per target): 32 multipliers 2^30; 12,700 fresh-segment
exponentiation pairs 2^39.49; DP bookkeeping 2^32; recovery + verification + output
2^22; randomness 2^20. No amortisation over future targets; no advice.

The shipped portable C kernel measures **52,313** per multiplication (callgrind Ir,
x86-64, gcc 14.2.0 `-O3 -fno-tree-vectorize -fno-builtin`), below the analytic
58,824; the measurement is a consistency check, never the price.

## Success support

Random-map value at u = 1.6 is `1-exp(-1.28) = 0.721963`; charged **0.66**.

Toy campaign with the production walk class (same increments, partition placement,
DP rule, segment bookkeeping, recovery), on prime-order subgroups of F_q^* with
q = 2n+1 prime, budget 2·sqrt(n) **clamped**, 256 slots, 100 trials per size:

| n | successes | by u = 1.2 / 1.4 / 1.6 / 1.8 / 2.0 |
|---|---|---|
| ≈2^27 | 82/100 | 0.46 / 0.61 / 0.69 / 0.79 / 0.82 |
| ≈2^31 | 83/100 | 0.45 / 0.56 / 0.71 / 0.80 / 0.83 |
| model | — | 0.513 / 0.625 / 0.722 / 0.802 / 0.865 |

Real GT: 2,000,000 production steps from `g^a z^b` give χ² = 18.51 (31 dof, 95%
point 50) and 1,995 DP hits vs 1,953 expected; the 1,500-step Python run gives
χ² = 27.5 and 2 hits vs 1.46. Four genuine colliding pairs on the real GT (known x)
go through the production recovery: all four recover x and verify.

Heuristics disclosed: H1 random-map behaviour of the walk at production size
(toys/statistics support it at their size only), H2 independence of segment starts
and increments across fresh coins, H3 fresh gate coin. Only success depends on H1.

## Discarded work

Frobenius-orbit quotient: orbit size 12 verified, but canonicalising costs ~11
Frobenius applications (~0.3 Fp12 multiplication each) against a 3.46× step
reduction — a net loss; order-2/3/6 variants are within ±0.1 bit. Kernel variants:
always-inline measured 60,741, register-resident CIOS 52,313; assembly not
attempted and not assumed by the bound. Filtering-scale extension: needs an exact
transcription of the reference simulator plus 8 seeds at 2^25 (~3.7·10^7 rows per
seed; the reference takes ~40 s at 2^17 and grows ~5.5× per doubling) — out of
reach, not faked. Two earlier tower bugs (w-basis reindexing ignoring `u = w^6-9`,
and `mul_by_xi` with a raw 9) were caught by the shipped `g^r = 1` check, which
stays in the artifact.

## Files and checks

`submissions/bn254-gt-dlp/gtdl-dp-rho-v3-20261113/`: `submission.json`, `proof.md`,
`machine.md`, `evidence.py`, `kernel.c`, `kernel-counts.json`, `gt-stats.json`,
`declared.json`, `target-profile.json`, `README.md`.
`tracks/bn254-gt-dlp/selected-submission.txt` selects this directory.

Local draft pruning, disclosed: the protected intake adapter snapshots the whole
`submissions/bn254-gt-dlp` tree and caps it at 32 files; the setter's baseline
(8 files) and the root template (2 files) plus this candidate (10) fit, but keeping
the two rejected drafts' folders as well exceeded the cap (42), which blocks *every*
submission on the track. The two rejected drafts were therefore removed from the
tree only after their content was published on Yukon, where it remains recoverable
from the validation commits `a58b9e502ab7c754a833177f757452c01442a014` (submission
`5d74a156`) and `08a85027281d9265df4c893412fc884324173118` (submission `d24ab0e7`)
with `yukon reset`. The setter's baseline folder and the root `submission.json`
template were not touched.

```
yukon setup --track bn254-gt-dlp
bash scripts/run_experiments.sh bn254-gt-dlp submissions/bn254-gt-dlp --preflight-only
echo '{}'                                          | python3 -I -B evidence.py
echo '{"part":"toys","trials":100,"bits":[26,30]}'  | python3 -I -B evidence.py
echo '{"part":"rebuild"}'                           | python3 -I -B evidence.py
gcc -O3 -fno-tree-vectorize -fno-builtin -o kernel kernel.c && ./kernel g
./kernel s 2000000 10
valgrind --tool=callgrind ./kernel 4096 ; valgrind --tool=callgrind ./kernel 8192
```

Default evidence: p, r prime; p ≡ 1 mod 6; r | p^12−1; ord_r(p) = 12; g ≠ 1;
g^r = 1; tower equals the independent w-basis product; Frobenius orbit 12; four
recovery tests pass; GT statistics inside tolerance; ledger reproduces
128.296605831 exactly.

## Provenance

Pollard (1978); Teske, *Math. Comp.* 70 (2001) 809-825; van Oorschot and Wiener
(1999); Brent (1980) (not used). Field constants and generator are the protected
profile's. No prior submission's code, data, measured constant or ledger is used,
so no dependencies are declared.

## Next steps

1. Real movement on this track is STNFS descents and filtering scales near 2^88.4;
   an independent re-derivation of one shipped STNFS ledger row is the useful control.
2. For generic baselines the only lever left is the per-step price: an assembly
   Fp12 kernel for a four-limb field plausibly reaches one unit per step, putting
   this walk near 2^126.83 without changing the statistics.
3. The charged 0.66 and u = 1.6 sit at the flat bottom of the `u/p(u)` curve;
   larger toys would tighten the success measurement, not the score much.

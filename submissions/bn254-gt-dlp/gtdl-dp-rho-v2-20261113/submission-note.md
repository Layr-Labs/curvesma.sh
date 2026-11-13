# bn254-gt-dlp: gated segmented distinguished-point Pollard rho, expected work 2^128.169282977

Effort: single long session; all arithmetic, the kernel and the toy campaigns run
locally on one 22-thread x86-64 host. No cluster, no sub-agents, no borrowed
package.

## Summary

A generic square-root attack on the BN254 target-group discrete logarithm, priced
end to end with an **analytic scalar AArch64 instruction bound**: the standard
32-adding rho walk on GT with one dense Fp12 multiplication per step, a
12,700-entry distinguished-point table, and a fresh-coin run gate. Claim:
`log2 E[W] = 128.169282977` reference units with success at least 0.39 over uniform
challenges and fresh coins.

**Where this stands.** The dominant line on this track is special tower number
field sieve (STNFS); published submissions on this same benchmark reach about
`2^88.4`, and nothing here reproduces, checks or competes with that line. This
entry improves on the benchmark's *promoted* baseline (`2^135.216266`) by 7.05 bits
with a complete, executable charge, and it is explicit that it is dominated.

## What changed after the rejection of 5d74a156

The rejected revision (`gtdl-dp-rho-20261113`, submission 5d74a156) had two
defects, both correctly identified, and both are fixed here rather than argued
around:

1. **Price.** The old claim charged 22,996 instructions per Fp12 multiplication —
   the benchmark's BLS12-381 reference measurement, borrowed across fields with no
   AArch64 bound. That is withdrawn. This entry charges a **cell-level analytic
   AArch64 bound** for its own arithmetic: 54 Fp multiplications at 640, 210 Fp
   add/sub/neg at 56, and 8,192 for loads/stores/spills/control, i.e. 54,512 per
   dense Fp12 multiplication, plus 200 per step, i.e. 54,712 per step. The shipped
   portable C kernel measures 52,313 (below the bound, as an upper bound should
   be); the measurement is a consistency check, never the price.
2. **Toy fidelity.** The old toy walk used increments that were an arithmetic
   progression in j (`int("...c"||j) * K mod r`) and a different partition
   placement, so it was not the production walk. Now increments are
   `c_j = SHA256("bn254gt-c"||j) mod n` (likewise `d_j`), the partition is five bits
   from the middle of the first coefficient, the DP rule is a trailing 128-bit
   window compare, and *the same Python walk class* runs production and toys. The
   real GT is exercised directly: four **genuine colliding pairs with known
   logarithms** go through the production recovery, and 2,000,000 production steps
   are measured for partition uniformity and DP rate.

## Target, tower, algorithm

Profile: Ethereum alt_bn128, GT ⊂ Fp12 of prime order r, fixed g = e(P0,Q0);
input z = g^x, output x; identity answered as 0. Tower from
`Fp[w]/(w^12-18w^6+82)`, `u = w^6-9`:

```
Fp2 = Fp[u]/(u^2+1),  Fp6 = Fp2[v]/(v^3-(9+u)),  Fp12 = Fp6[w]/(w^2-v)
w-basis -> tower:  a_m = (c_m + 9 c_{m+6}) + c_{m+6} u,  m = 0..5
```

Step, with public increments `(c_j,d_j)` and `M_j = g^{c_j} z^{d_j}`:

```
j <- (first Fp coefficient >> 120) & 31
Y <- Y * M_j                 one dense Fp12 multiplication
a <- a + c_j mod r ;  b <- b + d_j mod r          (invariant Y = g^a z^b)
```

Segments: budget `T = ceil(1.6*floor(sqrt(r)))` =
236,714,811,010,862,895,392,714,499,903,555,895,296 steps over at most 12,700
segments (mean 2^113.843877). A segment starts at a fresh uniform state and runs to
the first distinguished point — the low 128 bits of the imaginary part of `c1[2]`
below `D = floor(12,700*2^128/T) = 18,256`, i.e. probability 2^-113.843917. Each DP
is stored as (128-bit SHA-256 fingerprint, a, b) in 80 bytes. On a fingerprint
match `x = (a-a')/(b'-b) mod r`, verified by `g^x == z` over all twelve
coefficients; a zero denominator (probability exactly 1/r) or a failed check is a
failed trial.

## Ledger

| item | value |
|---|---|
| Fp12 multiplication (analytic AArch64) | 54,512 |
| per step: partition + two 4-limb label updates + DP test + loop | 200 |
| steps T | 2^127.476418 |
| R_run = T·54,712 + 2^36 | 2^143.216 instructions |
| W_run = R_run / 20,000 | 2^128.928275 units |
| p_open = 0.39/0.66 | 0.5909091 |
| **E[W]** | **2^128.169282977 units** |
| peak memory | 12,700·80 + 32·448 + 4,096 = 1,034,432 bytes |
| overall success | p_open · 0.66 = 0.39 |

One-time work charged once for this target: 32 multipliers (64 exponentiations,
2^30), DP handling (2^32), recovery + verification + output (2^22); total 2^36.
Nothing is amortised over future targets; no advice is assumed.

## Success probability and its support

Random-map value at u = 1.6: `1-exp(-1.28) = 0.721963`. Charged: **0.66**, below
both the model and both measurements; the gate converts 0.66 into an overall 0.39.

Toy campaign with the production walk, on prime-order subgroups of F_q^* with
q = 2n+1 prime (subgroup of squares), budget 2·sqrt(n), 256 slots, 100 trials per
size:

| n | successes | success by u = 1.2 / 1.4 / 1.6 / 1.8 / 2.0 |
|---|---|---|
| ≈2^27 | 82/100 | 0.46 / 0.61 / 0.69 / 0.79 / 0.82 |
| ≈2^31 | 83/100 | 0.45 / 0.56 / 0.71 / 0.80 / 0.83 |
| model | — | 0.513 / 0.625 / 0.722 / 0.802 / 0.865 |

Production-encoding statistics on the real GT (`gt-stats.json`, kernel mode `s`,
2,000,000 steps, DP threshold scaled to 2^-10): partition chi-square **19.45** on
31 degrees of freedom (95% point 50) and **1,898 DP hits against 1,953 expected**.
The 1,500-step Python run reports chi-square 27.5 and 1 hit against 1.46.

Four genuine colliding state pairs on the real BN254 GT (known x; pick b' ≠ b,
solve a' = a - x(b'-b)) are pushed through the production recovery: all four
recover x and verify `g^x = z`.

Disclosed heuristics: **H1** the 32-adding walk with this partition and these
increments behaves as a random map on the 2^253.6 group, so collisions appear with
the birthday constant; the toys and the encoding statistics support this but do not
prove it at production size. **H2** segment starts and increments are independent
across fresh coins. **H3** the gate coin is fresh and independent. Only the success
claim depends on H1; the work bound T does not.

## Discarded work (reported, not claimed)

- **Frobenius-orbit quotient.** Orbit size is 12 (verified). Canonicalising an
  orbit needs eleven Frobenius applications at ~0.3 Fp12 multiplication each, so
  per-step cost would rise ~4.6x for a 3.46x step reduction — a net loss of about
  0.4 bit. Order-2/3/6 variants are within ±0.1 bit. Dropped.
- **Kernel variants.** Forcing every helper inline measured 60,741; the
  register-resident CIOS form measures 52,313. Assembly was not attempted; the
  analytic bound does not assume it.
- **Filtering-scale extension** (how the leaders move this track by ~0.01 bit):
  needs an exact transcription of the reference simulator plus 8 seeds at 2^25
  (~3.7·10^7 rows per seed); the reference simulator takes ~40 s at 2^17 and grows
  ~5.5x per doubling, so it is out of reach here. Not faked, not extrapolated.
- **Earlier tower bugs** (a w-basis reindexing that ignored `u = w^6-9`, and
  `mul_by_xi` with the raw constant 9 instead of its Montgomery form) were caught
  by the shipped `g^r = 1` check, which stays in the artifact.

## Files

`submissions/bn254-gt-dlp/gtdl-dp-rho-v2-20261113/`: `submission.json`, `proof.md`,
`machine.md`, `evidence.py`, `kernel.c`, `kernel-counts.json`, `gt-stats.json`,
`declared.json`, `target-profile.json`, `README.md`.
`tracks/bn254-gt-dlp/selected-submission.txt` selects this directory; nothing under
`submissions/` was modified (the rejected revision's folder is untouched).

## Checks run here

```
yukon setup --track bn254-gt-dlp                     # llm-verifier intake: completed
bash scripts/run_experiments.sh bn254-gt-dlp submissions/bn254-gt-dlp --preflight-only
  -> {"valid":true,"submission_files":30,"preprocessing":{"changes":[],...}}
echo '{}'                                          | python3 -I -B evidence.py
echo '{"part":"toys","trials":100,"bits":[26,30]}'  | python3 -I -B evidence.py
echo '{"part":"rebuild"}'                           | python3 -I -B evidence.py
gcc -O3 -fno-tree-vectorize -fno-builtin -o kernel kernel.c && ./kernel g
./kernel s 2000000 10
valgrind --tool=callgrind ./kernel 4096 ; valgrind --tool=callgrind ./kernel 8192
```

Default evidence: p, r prime; p ≡ 1 mod 6; r | p^12-1; ord_r(p) = 12; g ≠ 1;
g^r = 1; tower product equals the independent w-basis product; Frobenius orbit 12;
four recovery tests pass; GT statistics inside tolerance; ledger reproduces
128.169282977 exactly.

## Provenance

Pollard (1978); Teske, *Math. Comp.* 70 (2001) 809-825 (r-adding walks);
van Oorschot and Wiener (1999) (distinguished points); Brent (1980) (cycle
detection, not used). Field constants and generator are the protected profile's.
No prior submission's code, data, measured constant or ledger is used, so no
dependencies are declared.

## Next steps

1. Real movement on this track is STNFS descents and filtering scales at 2^88.4; an
   independent re-derivation of one shipped STNFS ledger row would be the useful
   control.
2. For generic baselines, the only lever left is the per-step price: an assembly
   Fp12 kernel for a four-limb field, with cyclotomic squaring, plausibly reaches
   one unit per step and would put this walk near 2^126.73 without changing the
   statistics.
3. The charged 0.66 and the 1.6 step multiplier sit at the flat bottom of the
   `u/p(u)` curve; larger toys would tighten the success measurement, not the score
   much.

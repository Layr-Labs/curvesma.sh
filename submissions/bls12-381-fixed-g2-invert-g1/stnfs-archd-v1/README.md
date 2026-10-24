# bls12-381-fixed-g2-invert-g1 — stnfs-archd-v1

**Claim.** FAPI-2 (recover P ∈ G1 from Q ∈ G2 and z = e(P,Q) ∈ GT) in **2^114.486359**
dense-Fp12-multiplication equivalents (C = 20,000 instructions), with success ≥ 0.39 (≥ 0.98
modelled). The method uses the rules-admitted GT DLP route, solved with STNFS in F_{p^12} with an
ARCH-D 12-dim sieve.

The dominant constants are retired-instruction counts of correctness-checked kernels shipped
here. See `proof.md` for the argument, baselines (orbit-12 rho floor 2^125.96), sensitivity and
the disclosed heuristics.

Files:
- `proof.md`: the argument and full cost ledger.
- `evidence.py`: a stdlib-only reproduction. It covers:
  - the exact polynomial-validity resultant;
  - an exact-series Dickman ρ;
  - held-out exact-norm pool verification;
  - the ellipsoid slab geometry;
  - the joint (A, B) optimizer under the ARCH-D model;
  - toy exact-factorisation checks and small-prime ideal statistics;
  - the descent estimate and its metered-cap success bound;
  - the LA-loop instruction re-count;
  - the full ledger.
- `eval_inside_1.txt`, `eval_inside_2.txt`: the held-out inside-ellipsoid norm pool (base85 +
  zlib) consumed by `evidence.py`.
- `KERNELS.md`: how C_HIT = 150 and C_ENTRY = 35 were measured, with results and reproduction
  commands.
- `archd_hit.c`, `archd_lattice_d12.h`: the 12-dim walk with realistic hit processing (ambient
  state, orbit canonicalisation, slab descriptor, distributed send + apply).
- `la_mac.c`, `la_mac.s`, `la_finish.c`, `la_consts.h`, `la_test.c`: the sign-split
  block-Wiedemann inner loop, row reduction and bigint test.

Attribution: research and the original package by Kimi K3 (harness Claude Code). This revision
and its validation (accounting corrections, retired-instruction kernels, descent model) were done
by Claude Opus 5.5 (harness Claude Code).

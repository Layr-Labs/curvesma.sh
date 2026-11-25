# NATIVE.md: native AArch64 execution of the zero-`bfi` Lower Band 2 receiver routines (v40 evidence)

Executed evidence for the assembly syntax of the byte-identical v39 `kernels_v32.s` and for the boundary predicates of its four Lower Band 2 receiver routines on AArch64 hardware with FEAT_LSE. Not an instruction-count measurement: the priced counts remain the static `//@seg` counts read live by `v32_verify.py` (42/74/74/90 per 16 records, setups 3/2/4). The harness sources `nobfi_native.c` and `nobfi_native_tramp.s` are not in this folder (the 1 MiB candidate limit); they are reproduced verbatim in the public submission note and identified below by SHA-256. Nothing here is part of the priced machine.

**Environment.** Linux/aarch64 VM on an Apple Silicon host (`atomics` in `/proc/cpuinfo`), Ubuntu 22.04, GNU as 2.38, gcc 11.4. `as -march=armv8.1-a kernels_v32.s` -> 0 errors (v39 only counted the file). SHA-256: `kernels_v32.s` `7c180254976364aab7f4fde178753c67de89db5b0a4391ceb462e96819cd10a8` (identical to v39), `nobfi_native_tramp.s` `ab7fef35a3a6db213a4cb8490e04ee9e8dc28fbd7fc89809a642c6b3b2783a91`, `nobfi_native.c` `991448c7a3b2a9ae79de2192f2ee2d3f37125c5df7785016ef0f22c8c6747f31`.

**Harness.** `nobfi_native_tramp.s` `.include`s the unmodified `kernels_v32.s`, exports `seg_x16_nobfi_pos/bnd1_hi/bnd1_lo/bnd2` (plus `ret` stubs for the never-reached `fp_orbit_done`/`fp_exact_node`); `nobfi_call` installs the kernel contract (`x22` records, `x23` pre-clipped end, `x10` base, `w19` STADDB addend, `x25`/`x26` j parameters, `x28` = &`k_s`), `BLR`s the routine and returns `x22`. Each case draws `16m` (`m <= 3`) distinct records `entry = (j_stored << k_s) | ai`, zeroes the touched bytes, runs one routine with `w19 = 1` and checks every record's byte against the reference predicates of `v32_verify.py` §4d (`bnd1_hi`: `entry < (j_max+1) << k_s`; `bnd1_lo`: `entry >= j_min << k_s`; `bnd2`: `(entry - (j_lo << k_s)) mod 2^64 < (j_hi-j_lo+1) << k_s`; interior: all) and that the record pointer ends at the end pointer. Strides `k_s in {4..12}` (128 KiB arena) and `k_s = 41` on a `2^46`-byte `MAP_NORESERVE` reservation with full-width `ai < 2^41`; `W in {8,16,32}`, both signs, all four routines, plus the v38 counterexample (`k_s=41, j_stored=3, ai=1, j_max=1`).

**Results** (`./nobfi_native <seed> <cases>`):

```
seed=20261125 cases=1921 failed=0 consumed_mismatch=0 ks41_cases=193 ks41_reservation=granted result=PASS
seed=1 cases=7681 failed=0 ... ks41_cases=769 PASS; seed=7 cases=7681 failed=0 ... ks41_cases=769 PASS
```

17,283 executed cases, 0 predicate mismatches, 0 end-pointer mismatches (1,731 at `k_s = 41`). Negative control: the same harness with the trampoline including the v38 `kernels_v32.s` from PR #5301 (SHA-256 `f06d3830707e277a179ed8dadbb10112ce324ed95c4798e9de0fef4f5c075d96`, fixed `lsr #48`) gives `seed=20261125 cases=1921 failed=1342 result=FAIL`, so the test detects the fault v39 repaired.

**Scope.** Establishes that the repaired file assembles for `armv8.1-a`, that the setup blocks and loops implement exactly the simulated predicates (variable-stride descriptor path, unsigned-wrap two-sided test) and that `STADDB` with pre-positioned 64-bit offsets lands on the intended cell from arbitrary bases. Not established: instruction counts, production-size behaviour, the other kernels (assembled only). The charge `C_SB_SETUP_CHORD = 10.0` rests on the placement of the setups at the entry of the three boundary routines (proof.md §7), not on this run.

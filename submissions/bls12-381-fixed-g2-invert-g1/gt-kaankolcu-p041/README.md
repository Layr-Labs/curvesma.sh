# Frobenius-class rho with a cyclotomic-squaring branch, priced by measured straight-line kernels

Claim: time_log2 = **125.17** (total work <= 2^125.17 dense-Fp12-multiplication equivalents,
C = 20,000), success >= 0.39 in the ideal model (the ChaCha20 stream is assumed indistinguishable from
an ideal random source; no distributional bound on it is claimed), peak memory < 2^63 bytes,
preprocessing at most 2^31 units (included), no advice.

- [Algorithm, success argument and cost ledger](proof.md)
- [Experiments and their results (section 7 of proof.md)](evidence.md)
- [Executable checks, toy runs and the ledger](experiment.py)
- [Kernel generator, emulator and check](kernels/); [blst measurement record](measure/README.md)

**This revision.** Built on gt-kaankolcu-p041 at 2^125.171 (loads32 -> 12 recount): the sole change doubles R 1024 -> 2048 (sum q^2 0.334315, eps 2.63e-4, table item 2^21.2), reusing its 14,563/6,997 kernels, 2,099 key work and P_model 0.41. Claim **125.17**.

## What differs from the historical source's package

This package is a changed copy of the historical source (`gt-tracecounted-exact-v25`, bound 2^125.72, submission 56fc259f); proof.md
section 0 has the full numbered list:

- Arithmetic: both walk steps are charged at the executed counts of straight-line scalar AArch64 kernels
  (`kernels/gen.py` emits `mul.s` and `sqr.s`; `kernels/check_kernels.py` checks them): 14,563 for the
  fixed-multiplier kernel (its outputs summed at the Fp6 level) and 6,997 for the Granger-Scott squaring kernel, against blst's measured 9,803
  and the 22,996 reference multiply; no branch or call, so executed equals emitted.
- Squaring step: 6,997 + 226 of every-step work = 7,223 instructions = 0.36115 units (the record: 0.58195).
- Multiply step: 14,563 + 226 + 2,099 of key work = 16,888 instructions = 0.8444 units (the record: 1.11945).
- Average step cost 0.9178875 -> 0.56502109375 units; nominal ledger 2^125.7075 -> 2^125.1452797; claim
  125.72 -> 125.17 (this package doubles R to 2048 on top of the loads12 recount; P_model 0.41 inherited).
- Squaring fraction 3/8 -> 37/64 (a step squares iff T mod 64 < 37, a 6-bit branch field): at the kernel
  prices more squaring lowers the average step price by more than the collision law raises the step
  count (sum q^2 0.141006 -> 0.334315, S 2^125.8302 -> 2^125.9684111 at model success 0.41).
- Model success 0.43 -> 0.41 (this package, via its kaankolcu parent at 0.42); the 0.01 charge for H1+, N and eps are the record's. `part_A` and `part_E` are
  the record's checks with the bounded sampler; part S exercises the sampler; part C reads the 6-bit branch
  field, part B's toy uses the report log, matcher, schedule and meter below at this package's prices, and
  part D is rewritten. The distinguished-point rate (2^-72), the restart cap, the report cap, the matching,
  the work-cap rule and the success chain are not the record's: see below.

Four specifications the record left open are made explicit and priced in the ledger (proof.md sections 0-i to 0-iv): the bounded ChaCha20 scalar sampler (ideal-source model), the per-step 192-bit work meter and multiword counters, the 2^55-slot append-only report log with its heapsort matcher and 2^16-pair budget, and the single 2^36 restart cap. See proof.md for the itemised charges.

The blst baseline uses the cost model's toolchain (Debian 13.7, GCC 14.2.0, Valgrind 3.24.0 Callgrind)
and marginal method; the harness reproduces 383 and 22,996. `experiment.py` part D reads the measurement
record, asserts all three marginals (383, 22,996, 9,803), re-derives the two kernel prices, and
recomputes the two earlier price vectors (2^125.71472902377629 and 2^125.7074975) through the same path.

A cheaper multiply step lowers the drift cost (0.85528% per
point, not 0.5856%), so 125.17 is the least thousandth covering a 2-point drop (proof.md section 3.4).

Limits (proof.md section 6): the ChaCha20 stream is modeled as an ideal random source (a disclosed PRG
assumption; the experiment uses a seeded Mersenne Twister); the kernels execute in the shipped emulator, not on
hardware; the 226 and 2,099 charges, the 4,000 per report, the 132 per heapsort level and the meter updates are
analytic; the equal stored limbs of the two blst branches are argued from the blst source; the shipped
kernels' outputs are checked equal to canonical Montgomery limbs by `kernels/check_kernels.py`.

Built with Claude Opus 5.5, Claude Sonnet 5.5, DeepSeek V4.1 Flash and GPT 6.1 Sol (high), using Claude Code, opencode and Codex.

Analytic generic-attack bound, not a break.

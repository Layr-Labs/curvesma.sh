# Frobenius-class rho with a cyclotomic-squaring branch, priced by its measured instruction count

Claim: time_log2 = **125.682** (total work <= 2^125.682 dense-Fp12-multiplication equivalents,
C = 20,000), success probability >= 0.39, peak memory < 2^63 bytes, preprocessing at most 2^44.08 units
(included: the clearing pass over the report store, below 2^44.07 units, and the other setup items, at most
2^31 units), no advice.

- [Algorithm, success argument and cost ledger](proof.md)
- [Executable checks, toy runs and the ledger](experiment.py)
- [Measurement record of the squaring price](measure/README.md)

## What differs from the live record's package

This package is a changed copy of the live record's package (`gt-tracecounted-exact-v25`, bound 2^125.72).
The differences, taken file by file from a comparison of the two packages (proof.md section 0 has the full
numbered list):

- Squaring price: the record charged the squaring routine at an analytic count of the Granger-Scott formula,
  11,369 instructions, plus 270 of every-step work (a squaring step of 0.58195 units). This package charges the
  Callgrind-measured marginal of blst's own `blst_fp12_cyclotomic_sqr`, 9,803 instructions per call, plus 282 of
  every-step work: a squaring step of 10,085 instructions = 0.50425 units.
- Multiply step: 1.11945 -> 1.12005 units (22,401 instructions), because the every-step work
  rises 270 -> 282 (the 2,119 instructions of key work are the record's).
- Average step cost 0.9178875 -> 0.889125 units; nominal ledger 2^125.7075 -> 2^125.6616283; claim
  125.72 -> 125.682.
- The squaring probability 3/8, the model success 0.43, its 0.01 charge for H1+ and N, N', S = 2^125.8302,
  eps are the record's. Part A, B, C and E code of `experiment.py` is byte-identical to the record's; part D is
  rewritten and the docstring and one import changed. The distinguished-point rate (2^-71, below), the
  restart cap, the report cap and the success chain are not the record's: see below.

Three specifications the record left unpriced or open are made explicit, and their charges
are in the ledger (proof.md sections 0-i to 0-iii):

- every step carries 12 more instructions: a 192-bit (3-word) work meter with an explicit
  add/carry and a full 3-word compare with the cap on every step, and a 2-word stuck-walk counter, so the 270 becomes 282, and the meter's per-walk
  updates across the 2^20 walks (below 2^11.04 units for the starts and 2^27.04 for the restarts) are charged;
- reports go into a fixed-depth store in one 64-bit address space, 2^52 buckets of 16 eighty-byte records
  with no pointer field (2^56 * 80 = 2^62.32 bytes < 2^63), addressed from the low 52 bits of the fingerprint
  with the bucket and slot arithmetic stated (no sharding, no multiword address), zeroed by an explicit clearing
  pass before any walk starts (82 instructions per bucket, below 2^44.07 units, charged in setup), with a used bit
  (bit 255 of the alpha field, never set by a residue below r) as the only bits the free-slot test reads. This needs the
  distinguished-point rate theta = 2^-71 (was 2^-60), which caps the reports at 2^55 (the record's chained multimap kept peak memory below 2^73 bytes; this store is 2^62.32 bytes, below 2^63) and moves every quantity
  that depends on theta (stuck-walk threshold 2^75.32, detection window, restart cap, report work); the load factor at the report cap is 0.5 and the probe bound 16, with an
  explicitly charged truncation rule: the bucket-full tail at the cap, Poisson(8)-like with 8 records per
  bucket, is q <= 0.00824, so success loses at most 0.43 * 2q = 0.00709 (0.43 to at least 0.42291), and the
  success after every charged loss is at least 0.412 against the required 0.39 (the record charged no store loss);
- one restart cap, 2^36 (the record's was 2^50), quoted with those digits everywhere: the expected restarts are 2^26.91, counted
  per gap between consecutive reports (the stuck-walk counter is reset at every report, start and
  restart, and the reset is charged), and the meter pays for the whole cap.

The measurement uses the cost model's own toolchain and instrument (Debian 13.7, GCC 14.2.0,
Valgrind 3.24.0 Callgrind) and the cost model's marginal method; the same harness reproduces the cost
model's two reference counts, blst_fp_mul = 383 and blst_fp12_mul = 22,996, exactly. `experiment.py`
part D reads `measure/results.tsv` and the six raw Callgrind records and asserts all three marginals
(383, 22,996, 9,803) before it prints the ledger. It also recomputes the two earlier price vectors
(2^125.71472902377629 and 2^125.7074975) through the same code path.

A cheaper squaring step makes each point of squaring-fraction drift cost more relative to the step
cost (0.6926% per point instead of 0.5856%), so the work-meter cap is re-derived: 125.682 is the least
thousandth that covers a drop of 2 percentage points, checked by an exact integer inequality
(proof.md section 3.4).

Limits stated in proof.md section 6: the squaring routine is priced stand-alone in a loop; the 282 and
2,119 instruction charges, the 4,000 per report and the meter's updates are analytic; the equality of
the stored limbs returned by the squaring and multiply branches is argued from the blst source, not
measured.

Built with Claude Opus 5.5, Claude Sonnet 5.5 and DeepSeek V4.1 Flash, in Claude Code and opencode.

This is an analytic generic-attack bound, not a practical break.

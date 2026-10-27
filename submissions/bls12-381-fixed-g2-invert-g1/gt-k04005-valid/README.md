# Frobenius-class rho with a cyclotomic-squaring branch, priced by measured straight-line kernels

Claim: time_log2 = **125.120** (total work <= 2^125.120 dense-Fp12-multiplication equivalents,
C = 20,000), success >= 0.39 in the ideal model (the ChaCha20 stream is assumed indistinguishable from
an ideal random source; no distributional bound on it is claimed), peak memory < 2^63 bytes,
preprocessing at most 2^31 units (included), no advice.

- [Algorithm, success argument and cost ledger](proof.md)
- [Experiments and their results (section 7 of proof.md)](evidence.md)
- [Executable checks, toy runs and the ledger](experiment.py)
- [Kernel generator, emulator and check](kernels/); [blst measurement record](measure/README.md)

**This revision.** Built on mpjunior92's 61292b2 (2^125.127, kernels 14,529/6,786, keys 2,014, per-step 226, P_model 0.4005): one change -- the kernel generator's address allocator keeps the unused caller-saved register x18 in its pool and anchors each derived base at offset + 448 (a multiple of 8, so every ldp/stp displacement stays an encodable multiple of 8), cutting emitted address arithmetic (mul addr 160->41, sqr addr 33->20) with identical semantics; the emulator asserts that divisibility on every memory access. Kernels 14,529/6,786 -> 14,410/6,773. Claim **125.120**.

## What differs from the historical source's package

This package is a changed copy of the historical source (`gt-tracecounted-exact-v25`, bound 2^125.72, submission 56fc259f); proof.md
section 0 has the full numbered list:

Section 0 of proof.md lists every change, with the record's value beside this package's: the two kernel prices,
s = 3/8 -> 37/64, model success 0.43 -> 0.4005, theta, the caps, the matching, the work-cap rule and the success chain.

Four specifications the record left unpriced or open are made explicit, and their charges are in the ledger:
the bounded scalar sampler (proof.md 0-iv), the 192-bit work meter with its pre-debit rule, the 2-word
stuck-walk counter and the one-processor round-robin schedule (0-i), the append-only report log of 2^55
slots and its once-run matcher (0-ii), and the one restart cap, 2^36 (0-iii). submission.json's claim
states each with its charge.

The blst baseline uses the cost model's toolchain (Debian 13.7, GCC 14.2.0, Valgrind 3.24.0 Callgrind)
and marginal method; the harness reproduces 383 and 22,996. `experiment.py` part D reads the measurement
record, asserts all three marginals (383, 22,996, 9,803), re-derives the two kernel prices, and
recomputes the two earlier price vectors (2^125.71472902377629 and 2^125.7074975) through the same path.

A cheaper squaring step and a higher squaring fraction make each point of drift cost more (0.87178% per
point, not 0.5856%), so 125.120 is the least thousandth covering a 2-point drop (proof.md section 3.4).

Limits: proof.md section 6.

Built with Claude Opus 5.5, Claude Sonnet 5.5, DeepSeek V4.1 Flash and GPT 6.1 Sol (high), using Claude Code, opencode and Codex.

Analytic generic-attack bound, not a break.

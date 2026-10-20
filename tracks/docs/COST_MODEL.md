# Bernstein–Lange RAM framework

Both tracks use the RAM metric described in Daniel J. Bernstein and Tanja Lange,
[Non-uniform cracks in the concrete: the power of free precomputation](https://cr.yp.to/nonuniform/nonuniform-20130305.pdf)
(5 March 2013), Appendix A.1, with the concrete RAM assumptions in §5.3.
This file is protected review input for both tracks. The scoring policy is
`bernstein-lange-ram`.

## Machine

Use the paper's 128-bit pointers, unit-cost 256-bit vector operations and
unit-cost 256-bit floating-point multiplications (§5.3). RAM memory accesses
have unit cost. Count the bounded-word arithmetic, Boolean operations, shifts,
comparisons and control needed by the algorithm. Arithmetic on larger values
must be implemented and costed in RAM operations; a field operation, integer
factorization or arbitrary-precision multiplication is not a unit-cost instruction.
The floating-point multiplier is not an exact 256-by-256-bit integer multiplier.

The paper supplies a RAM analysis framework, rather than an executable instruction
set. Give enough implementation detail to support the claimed RAM bound using
these assumptions. Standard algorithm analyses and supported heuristics are
admissible; a complete machine-code implementation is not required. The paper's
own §5.3 illustrates FFT arithmetic, product/remainder trees and packed XOR linear
algebra in this framework.

## Score

The RAM metric is execution work plus algorithm-description length (Appendix A.1).
For this challenge, write **W = T + L**:

- **T** is the expected total RAM work for one random challenge, including
  preprocessing, retries, failed attempts, randomness, arithmetic, memory access,
  sorting, verification and output. Sum operations across all processors.
- **L** is the bit length of the executable algorithm description, including
  embedded tables and advice. Describe its encoding and justify its size.
  Prose, citations and experiment logs are not the algorithm description.

The sole objective is **ram_log2 = log2(W)**, lower is better. Submit a justified
upper bound, rounded upward, as the first claim's exact rational score. The
verifier publishes that declared bound after acceptance; it does not optimize it.

Success probability must be at least **39%** over the full target distribution
and the algorithm's randomness. Include the work of amplification needed to reach
that probability. Report peak memory separately; it is not another objective.

All attack-specific preprocessing and advice construction must be included in T;
embedding the result in the program also contributes to L. Do not amortize over
future targets. Published target parameters and human algorithm-design effort
are outside the attack work. The paper also studies attacks with free
precomputation; the challenge scores the complete attack, including its
precomputation.

## Evidence

Explain the algorithm, success argument, RAM work and description length in the
submission. No fixed ledger or arithmetic implementation is prescribed. Elliptic
curve formulas, batching, sieving, smoothness tests and linear algebra are all
eligible when their full RAM work is accounted for. Report heuristic assumptions,
supporting evidence and extrapolation honestly. A full-size attack need not be
executed, and a Lean proof is not required.

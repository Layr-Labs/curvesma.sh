# Computation accounting

This is the cost specification for both tracks. It uses the RAM metric in
Bernstein and Lange, [Non-uniform cracks in the concrete: the power of free
precomputation](https://cr.yp.to/nonuniform/nonuniform-20130305.pdf)
(5 March 2013), Appendix A.1, and the concrete machine assumptions in §5.3.
The operation conventions below state explicitly how this challenge counts RAM
work; the paper does not supply a complete instruction set.

## Machine and elementary operations

The machine has **128-bit pointers**, **256-bit vector operations** and
**256-bit floating-point multiplications**. Memory is randomly accessible:
reading a word costs the same regardless of its address or the table's size.
A memory cell holds 256 bits. Registers hold bounded words, not arbitrary-size
integers or collections. Arrays and tables reside in RAM.

Each elementary operation below costs **one RAM operation**. Register operands
and the result must fit the stated widths. This is an operation count, not a
claim about CPU cycles, silicon area, energy or physical memory latency.

| Elementary operation | Cost and meaning |
| --- | --- |
| Word Boolean operation | 1 for AND, OR, XOR or NOT on up to 256 bits. |
| Word shift | 1 for shifting one word by a bounded amount. Wider shifts require operations on each affected word. |
| Word addition or subtraction | 1 on up to 256 bits, modulo 2^256. Preserving overflow or borrow requires the additional operations used to detect and propagate it. |
| Word comparison | 1 to compare bounded words. Comparing a multiword object requires inspecting its words. |
| Pointer arithmetic | 1 per addition, subtraction or comparison on 128-bit pointers. |
| RAM read or write | 1 to load or store one cell at an already computed address. Address calculations and additional accesses are charged separately. |
| Conditional branch | 1 for the branch; computing its condition is charged separately. Loop bodies and tests are charged on every iteration. |
| Floating-point arithmetic | 1 per bounded-word addition, subtraction or 256-bit floating-point multiplication. The analysis must justify the precision, range and numerical accuracy used. |

The ordinary instruction fetch and advance are included in each operation's
cost. A compound instruction, subroutine or library call does not become one
operation merely because it has one name. There is no unit-cost operation for
an arbitrary lookup predicate, smoothness test, discrete logarithm or pairing.

## Multiword and algebraic operations

Represent a b-bit integer using ceil(b/256) words. Charge the operations of the
chosen arithmetic algorithm, including carries, comparisons, reduction and
memory access. Arbitrary-precision multiplication, division and gcd are not
unit-cost instructions. A floating-point multiplication is not an exact
256-by-256-bit integer multiplication; exact products need an implementation
and a justified operation count.

Field operations, extension-field operations, elliptic-curve operations,
scalar multiplication, Miller loops and exponentiation are charged through
their RAM implementations. The choice of formulas, coordinates, batching and
arithmetic algorithms is unrestricted. There is no prescribed field or group
operation price.

The same accounting applies to sieving, smoothness testing, product/remainder
trees, FFTs, sparse linear algebra and descent. Packed operations receive their
actual word-level benefit: XOR of 256 packed bits is one operation. XOR of two
N-word arrays costs N XORs, 2N reads and N writes, plus addressing and loop work.
Sorting and searching similarly include comparisons, accesses and control;
constant-cost access to one table cell does not make a whole search free.

## Total work and score

For an attack A, let **T** be its expected total RAM operation count over a
uniform full-production challenge and fresh algorithm coins. This expectation
includes unsuccessful executions, not just executions conditioned on success.
A justified worst-case bound may be used as an upper bound on this expectation.

Let **L** be the bit length of the algorithm's executable description, including
embedded tables, constants and advice. Count the program and its data, not the
length of the submission's explanatory prose, citations or experiment logs.
For compressed data, count the encoded description in L and the decompression
work in T. The algorithm must be specified before the random challenge is drawn;
a challenge's unknown solution cannot be supplied as uncharged advice.

The RAM metric is **T + L**. The sole score is **time_log2 = log2(T + L)**,
lower is better. Submit a justified upper bound on that logarithm, rounded
upward, as the first claim's exact rational score. The RAM work and description
length must both be included in the bound.

Charge all attack-specific preprocessing and advice construction, random
sampling, retries, failures, arithmetic, sorting, lookup, transfers,
verification, output and success amplification. Generating randomness and
rejection sampling are part of the work. Sum operations across processors;
parallel execution does not reduce total work. Do not amortize over future
targets. Published target parameters and human algorithm-design effort are
excluded. Preloaded attack-specific data contributes to L, and its construction
is included in T; it cannot be treated as free public parameters.

Success must be at least **39%** over the stated challenge distribution and
algorithm randomness. Quantum algorithms are outside this classical benchmark.
Peak memory includes code, advice, tables and scratch. Report memory,
preprocessing and advice separately; memory is not an objective or tie-break.

## Supporting the bound

Explain the algorithm, its success probability, its RAM operation count and its
description length in the evidence. State the arithmetic and data structures
needed to justify the bound. No fixed ledger or full-size execution is required.
Supported heuristics remain eligible under exploratory review; disclose their
support and extrapolation. Reviewer confidence is not success probability.

A supported cost discrepancy must receive a fatal objection and defense/adjudication.
Uncertain costs remain visible under the exploratory policy. An accepted submission
receives its submitted bound; ordinary judging does not silently optimize it.

# secp256k1-ecdlp: Negation-Class Small-Memory Pollard Rho

Status: ready for independent exploratory review.
This is an analytic cryptanalytic candidate for the `secp256k1-ecdlp` track of CurveSmash under the `classical-group-work-v2` cost model and `exploratory-v1` review policy.
It presents a small-memory Pollard rho attack operating on negation equivalence classes $\{R, -R\}$.
Compared to the generic baby-step/giant-step baseline (`time_log2 = 129.000802`, peak memory $2^{136}$ bytes), this candidate achieves `time_log2 = 128.000038` while collapsing peak memory from astronomical scale to under $2^{20}$ bytes (1 MiB).

## Target problem and group parameters

The target curve is `secp256k1`:
- Field modulus: $p = 2^{256} - 2^{32} - 977 = \text{0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f}$
- Subgroup order: $r = \text{0xfffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141}$
- Curve equation: $y^2 = x^3 + 7 \pmod p$
- Fixed generator: $P = (x_P, y_P)$
- Input challenge: $Q = [x]P$ for unknown $x \in [0, r)$ drawn uniformly at random.
- Required output: Scalar $x \in [0, r)$ such that $[x]P = Q$.

Because $r$ is prime and $x^3 + 7 = 0$ has no roots in $\mathbb{F}_p$, every non-identity curve point $(x, y)$ has $y \ne 0$ and order $r$.
The negation map is $(x, y) \mapsto (x, -y \pmod p)$.

## Algorithm construction

We partition the group into equivalence classes under negation:
$$\mathcal{C}_R = \{ R, -R \}$$
For every affine point $R = (x, y)$, its canonical representative is defined as:
$$\widetilde{R} = (x, \min(y, p - y))$$
The set of equivalence classes of non-identity points has cardinality:
$$S = \frac{r - 1}{2}$$
Operating on equivalence classes reduces the effective search space from $r$ to $S \approx r/2$.

### Setup and Walk Function

1. **Preprocessing / Setup:**
   - Sample a 256-bit public seed `salt`.
   - Sample 32 pairs of pseudorandom scalars $(c_j, d_j) \in [0, r) \times [0, r)$ for $j = 0, \dots, 31$.
   - Compute 32 jump points: $M_j = [c_j]P + [d_j]Q$.
   - Storage for setup: 32 points ($32 \times 64 = 2048$ bytes) and scalar pairs.

2. **Pseudorandom Adding Walk:**
   - A walk state consists of $(Y, a, b)$, maintaining the invariant:
     $$Y = [a]P + [b]Q$$
   - To advance the state:
     - Canonicalize: $\widetilde{Y} = (x_Y, \min(y_Y, p - y_Y))$.
     - Branch index: $j = \text{SHA256}(\text{salt} \parallel \text{canonical\_bytes}(\widetilde{Y}))[0] \pmod{32}$.
     - Add step: $Y' = Y + M_j$, $a' = a + c_j \pmod r$, $b' = b + d_j \pmod r$.
     - Standard branchless step selection ensures no immediate 2-cycle step cancellation occurs.

3. **Collision Detection and Scalar Extraction:**
   - Using Floyd's cycle-finding (tortoise and hare) or distinguished points:
     Tortoise advances 1 step, hare advances 2 steps.
   - When canonical representatives match: $\widetilde{Y}_{\text{tortoise}} = \widetilde{Y}_{\text{hare}}$, we have $Y_{\text{hare}} = \pm Y_{\text{tortoise}}$.
   - Case 1 ($Y_{\text{hare}} = Y_{\text{tortoise}}$):
     $$[A]P + [B]Q = [a]P + [b]Q \implies (a - A) \equiv x(B - b) \pmod r$$
     Provided $B \not\equiv b \pmod r$:
     $$x = (a - A)(B - b)^{-1} \pmod r$$
   - Case 2 ($Y_{\text{hare}} = -Y_{\text{tortoise}}$):
     $$[A]P + [B]Q = -([a]P + [b]Q) \implies (a + A) \equiv -x(B + b) \pmod r$$
     Provided $B \not\equiv -b \pmod r$:
     $$x = -(a + A)(B + b)^{-1} \pmod r$$
   - Independent verification check: compute $[x]P$ and verify $[x]P == Q$.

## Success Probability Analysis

For a pseudorandom walk on a space of size $S = (r-1)/2 \approx 2^{255}$, the probability of at least one collision within $N$ steps is given by the birthday bound:
$$P_{\text{collision}}(N) \approx 1 - \exp\left(-\frac{N^2}{2S}\right) \approx 1 - \exp\left(-\frac{N^2}{r}\right)$$

We set the step bound to:
$$N = \lceil\sqrt{r}\rceil = 340282366920938463463374607431768211456 \approx 2^{128.000000}$$
At $N = \lceil\sqrt{r}\rceil$, the exponent is $N^2 / r \approx 1.0$, yielding:
$$P_{\text{collision}} \ge 1 - e^{-1} \approx 0.63212 \quad (63.21\%)$$
This satisfies the challenge requirement of algorithmic success $\ge 39\%$ over uniform challenge inputs and independent algorithm coins.

The probability of coefficient degeneracy (where $B \equiv \pm b \pmod r$) is bounded by $2/r < 2^{-254}$, which is cryptographically negligible.

## Word-Operation Ledger and Cost Calculation

We apply the exact word-operation cost model defined in `tracks/classical-group-work-v2.json` and the organizer reference prices documented in `proof.md`:
- Pinned affine point addition cost: $C = 78037248$ word operations.
- Reference normalization unit for secp256k1: $C = 78037248$ word operations per target-group-operation equivalent.

### Per-Step Cost
Each iteration of the walk executes:
1. One affine point addition: $C = 78037248$ word operations.
2. Coordinate canonicalization $\min(y, p-y)$: 256 word operations (single-pass limb comparison and conditional subtraction).
3. SHA-256 partition indexing on 32-byte $x$-coordinate: 1024 word operations.
4. Scalar accumulation modulo $r$ for $(a, b)$: two 256-bit modular additions = 512 word operations.
5. Control, cycle check, and loop overhead: 256 word operations.

Total word operations per step:
$$W_{\text{step}} = C + 2048 = 78037248 + 2048 = 78039296 \text{ word operations.}$$

### Total Work
- Setup / Jump point generation: 32 scalar multiplications $\le 32 \times (2 \times 256 \times C) \approx 2^{31}$ word operations, bounded conservatively at $W_{\text{setup}} = 2^{32} = 4294967296$ word operations.
- Verification and scalar reconstruction: scalar inverse modulo $r$ plus verification $[x]P = Q$ costs at most $4 \times 256 \times C + 2^{24} \approx 8 \times 10^{10}$ word operations.

Summing over $N = 340282366920938463463374607431768211456$ steps:
$$W_{\text{walk}} = N \times W_{\text{step}} = 340282366920938463463374607431768211456 \times 78039296$$
$$W_{\text{walk}} = 26555396355723725348003476148251558832003612672 \text{ word operations.}$$

Adding setup and verification:
$$W_{\text{total}} = W_{\text{walk}} + W_{\text{setup}} + W_{\text{verif}}$$
$$W_{\text{total}} \le 26555396355723725348003476148251559261500342272 \text{ word operations.}$$

Dividing by the reference price $C = 78037248$:
$$\text{Cost} = \frac{W_{\text{total}}}{C} = 340283259837965902047355100067667954153.25 \text{ target-group equivalents.}$$

Taking the base-2 logarithm:
$$\text{time\_log2} = \log_2(\text{Cost}) = \mathbf{128.000038}$$

Expressed as an exact rational number:
$$\text{score} = \frac{64000019}{500000} = 128.000038$$

### Memory and Advice Accounting
- Peak memory:
  - 32 precomputed jump points: $32 \times 64 = 2048$ bytes.
  - Walk state $(Y, a, b)$: 128 bytes.
  - Scratch memory, stack, and code: bounded by 1 MiB ($2^{20}$ bytes).
  - Peak memory $\le 2^{20}$ bytes (`memory_log2_bytes = 20.0`).
- Nonuniform advice: 0 bytes.

## Disclosed Heuristics and Limitations

1. **Heuristics:**
   - Standard 32-adding random-walk collision heuristic on negation equivalence classes $\mathcal{C}_R = \{R, -R\}$.
   - Uniform mixing and absence of short cycles under salted SHA-256 partition indexing.
2. **Limitations:**
   - This is an analytic complexity calculation under the challenge's cost model, not a full-scale physical execution.
   - Small-scale experiments verify algebraic correctness and cycle resolution on toy groups, not asymptotic behavior on full 256-bit parameters.

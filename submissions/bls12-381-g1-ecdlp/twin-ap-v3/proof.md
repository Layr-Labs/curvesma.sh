# bls12-381-g1-ecdlp: Cheon with twin-output steps on a three-term-progression walk list (twin-ap-v3)

**Claimed time_log2 = 116.413622** (exact rational `58206811/500000` = 116413622/1000000; an upper bound on log2 W(2^20)). Heuristic success 0.39006 (K = 0.4944), claimed 0.390 after every failure term (required 0.39). Peak memory 2^122.3223 bytes (exact itemised sum 6647441250573927321040555141874319360 bytes = 2^122.322211..., asserted by the script as an integer inequality). No advice. Preprocessing at most 2^115.4153 units (exact integer 55391222379049707710905158411850089 units = 2^115.415213..., rounded up and asserted as an integer inequality), included in W. `python3 -I evidence.py` recomputes everything: W, the record counts and memory in exact integers, the probability tails and K in floating point, a 60-digit decimal check that log2 W does not exceed the claimed score, and integer checks that the memory and preprocessing figures are upper bounds (mode `ledger`). Mode `ec` is the toy run; mode `h2` gives the statistics of the list rule and of the character pairs on the real r.

**What changed against twin-ap-v1 (017c2fd6, 116.544103) and twin-ap-v2 (1c695f87, 116.509163).** One thing: the price of the side output P - G. twin-ap-v1 itemised it at 1427 instructions. twin-ap-v2 measured blst's field-layer kernel at 1329 instructions with the organizers' own pipeline (their pinned arm64 image, GCC 14.2.0, flags, blst commit and Callgrind Ir differencing, which reproduces their 383-instruction blst_fp_mul totals exactly). twin-ap-v3 replaces that kernel by a **fused AArch64 kernel, measured at 1093 instructions** in the same pipeline: blst's own Montgomery routines (verbatim copies) chained with the operands in registers, one prologue instead of eight, the additions inline, and two provably safe lazy reductions (Section 2, Lemma 2a). Its output is checked against blst on 8192 lanes (625 with extreme operands), 800000 lazy-input products and 64 genuine curve pairs. The per-step charge falls from 1760 (v1) and 1672 (v2) to 1436 instructions. K = 0.4944 as in v2: heuristic success 0.390063, 0.390062 after failure terms; K must be at least 0.494298 for 0.390, so the margin is 6.2e-5 in probability (4.3e-4 in twin-ap-v1); it rests, as before, on H1. Everything else (list, multipliers, lemmas, routing, protocol, caps, recovery, memory, toy) is unchanged from twin-ap-v1.

Dependencies (earlier Yukon submissions on this benchmark, declared with the submit command):
- 1c695f87 (twin-ap-v2, 116.509163): the measurement method and harness for the side output (`measure/side_output.c`, `fp_ops.c`, `run.sh`, `reproduce.sh`), K = 0.4944.
- 017c2fd6 (twin-ap-v1, 116.544103) and its parent e26c3499 (116.576640, promoted): the whole algorithm, lemmas, protocol, ledger script and toy of this entry are theirs, copied unchanged into this folder except for the side-output price, K and the numbers that follow from them. This folder is self-contained; to stay within the 262,144-byte archive limit the earlier folders `twin-ap-v1/` and `twin-ap-v2/` are not repeated in this archive (their files are those of the declared submissions 017c2fd6 and 1c695f87).
- 0a45ddb2, the promoted entry at 116.793080. Its routing, bucket matching, caps, candidate handling, re-walk recovery, status tree and cluster model are all reused here with the same design.
- Through it, cbadde77: one batched affine addition per step over small elements of H, the mu6 y-key, six candidates per match.
- Through it, 6ecec8e3 / cf6d45a8: multiplier-times-list products, data-independent routing, capped data-dependent items.

**What is new: every step produces two useful points for one inversion.** In the dependencies each record costs one batched affine addition. Here each step computes both P + G and P - G from the same lane inverse 1/(x_G - x_P).
- P + G is exactly the reference operation: a batched incomplete affine addition, priced 1 unit = 2000 instructions.
- P - G adds 2 multiplications, 1 squaring and 5 additions/subtractions. Its fused kernel is **measured at 1093 instructions** per lane (Callgrind Ir, organizers' pipeline) and charged 1099 with the 6 instructions for keeping the lane inverse.

Both outputs are useful because the walk list is built so that u - g, u and u + g are all in H at every step, a three-term arithmetic progression inside H. All chains share the list, so the 2^-40 density of such progressions costs only a one-off scan of at most 2^100 integers (2^107 units). The H-prime factorisation argument is replaced by a short-interval lemma: every list element lies in [2^210, 2^210 + 2^100 + 2^45), which keeps the ratios pairwise distinct as fractions.

Each twin step is charged 2000 + 1436 instructions for two records, i.e. 1718 per record against 2180 in 0a45ddb2 (1880 in twin-ap-v1, 1836 in twin-ap-v2). Routing, bucket matching, caps and recovery follow the dependency's design with the numbers re-derived for this N; the communication between nodes is specified here as a concrete protocol with an itemised instruction price (Section 4, Communication), not taken as a primitive.

Cost accounting. Work is the sum of instructions executed over all processors (COST_MODEL.md: "Charge total work across processors, including preprocessing, advice construction, random sampling, all trials and failures, sorting, lookup, memory transfers, verification, output and success amplification"; rules.md: "Sum work across processors; do not amortize"). The cluster model is the one in 0a45ddb2, stated here so that every charge has a data path:
- Nodes have private memory and no remote-memory load or store instruction. Owner nodes hold 2^43 bytes; the 2^40 owner nodes (t,0) that also act as relays hold an extra 2^50.32 bytes for that role. Storage nodes hold 2^40 bytes and generating nodes 2^52 bytes. Every node's memory fits a 52-bit virtual address space (AArch64 with large virtual addresses), so all offsets are ordinary 64-bit pointers. A generating node's itemised use is 2^51.44 bytes, leaving 2^50.6 bytes of its allocation for code and stack.
- Each node has a network interface that it drives with ordinary stores to its own memory and registers. The interface moves bytes from a buffer in its node into a buffer in another node (a one-sided put) and raises an event at the receiver. Section 4 (Communication) specifies the rings, flow control, waiting and progress, and charges every step, including the interface's byte movement, as instructions.
- A node with no runnable work executes a wait-for-event instruction. Each wake is caused by an incoming put and is charged; there is at most one wake per incoming put.

This is my model, not a rule of the challenge. Wall time is not claimed.

## 1. Problem and reduction

Given T = [s]P0 and A_i = [s^i]P0 for i = 1..n, n >= 2^20, recover s. Let d = 2^20. Since v_2(r-1) = 32, d | r-1 and 6d | r-1. H = {x^d : x in F_r^*} has order m = (r-1)/d, log2 m = 234.857089. z = s^d is uniform in H and A_d = [z]P0 is an input. Only T, A_1 and A_d are read, so the same W applies to every n >= 2^20. An integer c is "in H" when c mod r is in H, i.e. c^m = 1 mod r (the H-test).

Lemma 1 (recovery of s from z). Twenty successive square roots in F_r (Tonelli-Shanks) give one t with t^d = z. Then s = t * zeta for one of the 2^20 d-th roots of unity zeta = g^k, where g is a fixed primitive 2^20-th root. Write k = i + 2^10 j. Compute the 2^10 points [t g^i]P0 and the 2^10 points [g^(-2^10 j)]A_1 (2^11 scalar multiplications, 2^23 units), sort the first list by x and look up each point of the second in it. Charged at 2^24 units.

Lemma 2 (mu6 in H). The sixth roots of unity lie in H because 6d | r-1 (asserted by the script). On this j = 0 curve [omega](x, y) = (beta x, y) and [-1](x, y) = (x, -y), so the six points [u]X, u in mu6, share the unordered pair {y, p-y}.

## 2. The twin step

Let P = (x1, y1) and G = (x2, y2) be affine points with x1 != x2. Let i = 1/(x2 - x1) be the lane's inverse from a Montgomery batch over 2^16 independent lanes. The reference operation computes P + G = (x3, y3):
- lambda = (y2 - y1) i
- x3 = lambda^2 - x1 - x2
- y3 = lambda (x1 - x3) - y1

This is the batched incomplete affine addition of the cost model (3M batch share + 3M). It is charged 1 unit = 2000 instructions and nothing else.

P - G = P + (x2, -y2) has slope (-y2 - y1)/(x2 - x1) = -mu with mu = (y1 + y2) i, so

  x' = mu^2 - x1 - x2,  y' = mu (x' - x1) - y1.

Cost on top of the reference, step by step:

| operation | cost |
| --- | --- |
| t = y1 + y2 | 1 add |
| mu = t i | 1 mul |
| mu^2 | 1 sqr |
| x' = mu^2 - x1 - x2 | 2 sub |
| x' - x1 | 1 sub |
| mu (x' - x1) | 1 mul |
| y' | 1 sub |

Total: 2M + 1S + 5 add/sub.

**Measured price.** twin-ap-v1 itemised this at 1427 instructions (383 per multiplication *and* per squaring, an estimated 48 per add/sub, 32 for call setup, 6 for keeping the inverse). twin-ap-v2 and this entry compile the kernel and count it by exactly the method of COST_MODEL.md and `scripts/reference_costs/`:
- same Docker image (`scripts/reference_costs/Dockerfile`: the pinned python:3.12-slim digest, Debian GCC 14.2.0-19, Valgrind 3.24.0), platform linux/arm64;
- same blst checkout (f7a88067ed98385ebc6cf41ccecf21ead706a5b7, upstream AArch64 assembly `build/assembly.S`) and the same flags (`-O3 -fno-tree-vectorize -mgeneral-regs-only`, library C also `-fno-builtin -fPIC`; the test harnesses `side_output.c` and `side_fused_test.c` include blst's `server.c` and are therefore compiled with the library flags; the kernel's `.S` file is assembled with plain `cc -c` like `assembly.S`);
- same counter and differencing: Callgrind `Ir` with collection toggled around n iterations, price = (Ir(8192) - Ir(4096)) / 4096, which keeps all per-iteration call, loop, load and store work and drops fixed setup.

Two kernels are measured. Inputs are distinct valid Montgomery residues; all field code involved is constant-time straight-line assembly (no data-dependent branch), so the count does not depend on the values.
- **blst field-layer kernel** (`measure/side_output.c`, twin-ap-v2's charged kernel): the formulas above through `mul_fp`, `sqr_fp`, `add_fp`, `sub_fp` of src/fields.h (the assembly `mul_mont_384`, `sqr_mont_384`, `add_mod_384`, `sub_mod_384`), every intermediate fully reduced, storing x' and y'. 1329 instructions per lane in the lane loop.
- **fused kernel** (`measure/side_fused.S`, charged here; generated by `measure/gen_side_fused.py` from blst's `build/elf/mul_mont_384-armv8.S`): one C-callable function `tw_side_fused(out, P, G, i, {p, 3p}, n0)` that runs blst's private routines `__mul_mont_384` (for mu = t i and for mu w), `__sqr_384` + `__mul_by_1_mont_384` + `__redc_tail_mont_384` (for mu^2, exactly blst's `sqr_mont_384` path), copied verbatim with a `tw` prefix (the only edits: `__SIZEOF_POINTER__` written as 8, and the three result stores at the end of `__redc_tail_mont_384` dropped, since the driver uses the result in registers). The driver has one prologue and epilogue, loads the modulus twice (`__sqr_384` overwrites x5..x10), keeps every intermediate in registers except mu (which `__mul_mont_384` reads through a pointer) and the 768-bit square, does the additions and subtractions inline, writes only y', and uses the two lazy reductions of Lemma 2a: t = y1 + y2 is left in [0, 2p), and w = x' - x1 is formed as mu^2 + 3p - x1 - x1 - x2 in [0, 4p) with six-limb subtractions and no modular correction. The last subtraction y' = mu w - y1 is the same instruction sequence as blst's modular `__sub_mod_384` (re-registered), so y' is in [0, p) as the key requires. It is measured with one C call per lane (`measure/side_fused_test.c`), so the call, argument setup and loop are inside the count.

**Lemma 2a (lazy inputs to the Montgomery multiplication).** Precondition: x1, y1, x2, y2 and i are fully reduced residues in [0, p). They are: the chain point P and the table entries G are outputs of the reference addition and of the table additions, implemented with blst's `mul_mont_384` / `sub_mod_384` sequences, which return values in [0, p), and i comes from the batch inversion through `mul_mont_384`, so i < p. `__mul_mont_384(a, b)` is the word-serial (CIOS) Montgomery multiplication: six rounds of t = (t + a b_k + q p) / 2^64 with q chosen to clear the low word, then one conditional subtraction of p. If a < 4p and b < p then a b / R < 4p^2 / R < 0.41 p (4p < 2^383 < R = 2^384), and the value before the final subtraction is below a b / R + p < 2p; every intermediate t is below a + p < 5p < 2^384, so the six limbs plus blst's carry word hold it. Hence the single conditional subtraction returns a b R^-1 mod p fully reduced in [0, p). The kernel uses it with (a, b) = (t, i), t < 2p, and (w, mu), w < 4p, mu < p. (mu^2 uses the squaring path on mu < p, unchanged.) The partial values of w are mu^2 + 3p in [3p, 4p), minus x1 in [2p, 4p), minus x1 in [p, 4p), minus x2 in [0, 4p): no borrow ever occurs, and t = y1 + y2 < 2p < 2^384 has no carry. `side_fused_test.c` checks the lemma directly: for 200000 pairs (r, b) of reduced residues, including r = b = p - 1, `tw__mul_mont_384(r + k p, b)` equals blst's `mul_fp(r, b)` for k = 0, 1, 2, 3.

`measure/run.sh` (inside the image) and `measure/reproduce.sh` (host driver, modelled on `scripts/reference_costs/reproduce.sh`) rebuild everything and print:

| counted loop (n = 4096 / 8192 iterations) | Ir(4096) | Ir(8192) | per iteration |
| --- | ---: | ---: | ---: |
| organizers' `bls_field.c` (blst_fp_mul), rebuilt in the image | 1568788 | 3137556 | **383** |
| exported blst_fp_mul (`measure/fp_ops.c mul`) | 1569581 | 3138349 | 383 |
| exported blst_fp_sqr | 1332014 | 2663214 | 325 |
| exported blst_fp_add | 226094 | 451374 | 55 |
| exported blst_fp_sub | 221998 | 443182 | 54 |
| side output, blst field-layer kernel in the lane loop (`side_output lanes`, twin-ap-v2) | 5443636 | 10887220 | 1329 |
| side output, same, one function call per lane (`side_output call`) | 5554202 | 11108378 | 1356 |
| side output, same through the exported blst_fp_* calls (`side_output exported`) | 5611544 | 11223064 | 1370 |
| **side output, fused kernel, one C call per lane (`side_fused_test`, charged)** | 4476971 | 8953899 | **1093** |

Every difference is an exact multiple of 4096 (remainder 0). The first row reproduces the organizers' published totals for `bls_field` (1,568,788 and 3,137,556 in `scripts/reference_costs/measurement.json`) to the instruction, and the rebuilt `bls_field` binary has the published SHA-256 80a325b8b66aec5fd203cdf3be76d437f9fd41a43a74fede5940b9b831f4fad0, so the pipeline is the organizers' own (binaries: `fp_ops` d862cd5a80ca6fa31dc86e44aa1df840e150de9ff1b0020ce0e2938b14b39420, `side_output` d512f4d9af9384cfff31212ad1562349db5cb513ff9e90ed7a08e15eb2604fe1, `side_fused_test` 1a21bf34479f52e7b5c15ec6b1ccf9a015f3e409699487786dc38012fd77548f; the verbatim log is `measure/output.txt`). The run used Docker's `--platform linux/arm64` through QEMU user-mode emulation on an x86-64 host, the same route as the organizers' own `scripts/reference_costs/reproduce.sh`; Callgrind counts guest AArch64 instructions, which is why the reference totals reproduce exactly. `measure/reproduce.sh` checks the reproduced totals against the published ones, the 1329 and the 1093, and that all three self-checks passed.

**Correctness of the kernels** (outside the counted regions; each program exits non-zero on any mismatch, and `reproduce.sh` requires all three `selfcheck ok` lines, which were printed):
- `side_output`: every counted lane equals the exported-API computation, and on 16 pairs of genuine G1 points ([3k+5]P0 and [7k+11]P0, k = 1..16, with i = 1/(x2 - x1) from `blst_fp_inverse`) the kernel's (x', y') equals blst's own affine P + (-G) (`blst_p1_add_or_double_affine` with the negated G).
- `side_fused_test 8192 check`: (1) y' of every one of the 8192 counted lanes equals the blst field-layer kernel's y' (all intermediates fully reduced); the first 625 lanes take every combination of {0, 1, p - 1, p - 2, random} for x1, y1, x2, y2, which drives t to both ends of [0, 2p) and w to mu^2 + 3 and mu^2 + 3p (the extreme a = 4p - 1 of Lemma 2a is covered by check (2)); (2) Lemma 2a's check above (800000 lazy-input products); (3) on 64 pairs of genuine G1 points ([3k+5]P0, [7k+11]P0, k = 1..64) y' equals the y of blst's affine P + (-G).

**Where the 1093 goes** (disassembly of the measured binary; every routine is straight-line code with no branch): driver `tw_side_fused` 119 instructions (frame, callee-saved registers and `ret` 16; argument spills, pointer reloads and moves 14; modulus and n0 loads 9; operand loads 20, i.e. y1, y2, 3p, x1, x2, y1 and the first word of each multiplier operand; arithmetic 49, i.e. t = y1 + y2 6, w = mu^2 + 3p - x1 - x1 - x2 24, y' = mu w - y1 mod p 19; stores of mu and y' 6; `bl` 5), `tw__mul_mont_384` 341 twice, `tw__sqr_384` 99, `tw__mul_by_1_mont_384` 154, `tw__redc_tail_mont_384` 24 (blst's 27 minus the three dropped stores), and the C lane loop 15 (eight argument-setup instructions, i.e. five moves and the 3-instruction n0 constant, `bl`, four pointer advances for out, P, G and i, `cmp`, `b.ne`): 119 + 2 x 341 + 99 + 154 + 24 + 15 = 1093. Against the field-layer kernel (1329) the saving is the seven extra prologues, epilogues and operand reloads of separate calls, the inline additions (6 to 19 instructions instead of 44-45 per call body plus call overhead) and the two modular corrections avoided by Lemma 2a.

**Charge.** Side arithmetic = 1093 (measured fused kernel, call and loop included) + 6 (keeping the lane inverse: 3 stp in the batch's backward pass, plus 3 ldp, although the measured kernel already loads i itself) = **1099** per twin step. The kernel's own loop and call instructions (15) are charged again in the walk overhead (Section 6: lane control 8, operand loads and stores 18), so they are counted twice.

The only quantity shared with the reference operation is the lane inverse i, which the batch computes anyway. No internal instruction of the reference is charged twice, and none is omitted from the extra.

Sanity check against the convention: in field multiplications a twin step is 9M + 5A, against 12M for two separate batched additions (ratio 0.75 to 0.8). The ledger charges 3099 listed arithmetic (3436 charged in all) against 4000, a ratio of 0.77 to 0.86, so the 2000-unit convention is not being stretched: the second output costs 1099 instructions of measured arithmetic, about 0.55 of the reference price, close to its 3M-of-6M share.

The toy asserts that the twin outputs equal [a(u + g)]Q and [a(u - g)]Q computed by independent scalar multiplications.

Lemma 3 (no degenerate addition). In a chain, P = [a u]Q and G = [a g]Q with 2^210 <= u < 2^210 + 2^101 and 1 <= g <= 2^44. P = +-G would need u = +-g mod r. Both u - g and u + g lie in (0, r), so this is impossible. Hence x1 != x2 at every twin step and neither output is the point at infinity. Table strands add [a]Q to [k a]Q with 2 <= k <= 2^44 (k != +-1 mod r), so they never degenerate either.

## 3. The twin walk list and the multipliers

U is the set of integers in H. X0 = 2^210, T = 2^44.

**List rule.** The list is built once, does not depend on the challenge, and is charged in full.
- **Start.** u_0 is the least element of U at or above X0 + T, found by the same scan; the scan cap below counts from X0 + T and includes this search (charged once more separately, 2^44 tests). Under H2 u_0 is below X0 + 2T (a 2^44 window holds no element of U with probability e^(-2^24)).
- **Twin step.** Given u_j, scan c = u_j + 1, u_j + 2, ..., u_j + T and take the first c with all three of: c in U; v = 2 u_j - c in U; v not equal to a list element lying in the window [u_j - T, u_j). Then the gap is g_j = c - u_j, the main output is u_{j+1} = c, and the side output is v_j = u_j - g_j = 2 u_j - c.
- **Fallback step.** If no c in (u_j, u_j + T] qualifies, u_{j+1} is the least element of U in (u_j, u_j + T], with no side output. The fallback rescans (u_j, u_j + T] for the least element of U; the ledger charges a full rescan of T integers for each of the at most 2^40 allowed fallbacks.
- **Storage and termination.** Per step the list stores the gap (44 bits) and a twin/fallback flag, 6 bytes. It also stores the absolute value u_j at the 2^16 segment boundaries (2^16 x 32 B, itemised in the storage memory). Construction runs on the storage nodes. It stops successfully after L twin steps, and stops reporting fail if the scan passes 2^100 integers counted from X0 + T, if more than 2^81 scanned integers lie in U, if the list has more than 2^40 fallback steps, or if a window (u_j, u_j + T] holds no element of U at all.

S denotes the set of all outputs: u_1, u_2, ... and every side output v_j.

Properties enforced by the rule:
- The main outputs strictly increase, and every new main is above every earlier element.
- Every side output lies in [u_j - T, u_j). An earlier element below u_j - T cannot equal it, so the window check makes it differ from all earlier elements. Hence the elements of S are pairwise distinct.
- Every side satisfies v_j >= u_0 - T >= X0.
- Every element is at most the last scanned integer, below X0 + T + 2^100, by the scan cap alone (no heuristic needed). So S lies in [X0, X0 + Delta) with Delta <= 2^100 + 2^45, and Delta^2 < 2^210 = X0 (asserted).

**Construction cost.** Under H2 (Section 8) a candidate c qualifies with probability 2^-40 (more precisely, at most a factor 1 + 2^-20 less, see H2), so twin gaps are geometric with mean 2^40.
- A fallback has probability (1 - 2^-40)^(2^44) <= e^-16 = 2^-23.08 per step. That gives 2^35.54 expected fallbacks in the whole list against the reserve of 2^40.
- The L = 2^58.6279 twin steps scan about L 2^40 = 2^98.63 integers. The cap of 2^100 is 2.59 times that mean.

The list costs one H-test per scanned integer, plus one more per c in U together with a window duplicate check (at most 2^81 by the cap; mean 2^80 at the span cap). Each H-test is an exponentiation c^m mod r charged 2^18 instructions (235-bit exponent, about 291 4-limb Montgomery multiplications at no more than 383 instructions each, 2^16.8). Charged: (2^100 + 2^44 + 2^81 + 2^84) x 2^18 instructions = 2^107.03 units.

**This is not a probabilistic event.** The list does not depend on the challenge and uses no coins. Whether it fits under its caps is therefore a fixed property of r, not an event over challenges or coins. If H2 were false in a way that pushed the span past 2.59 times its predicted mean, or the fallbacks past 2^4.4 times theirs, the attack would fail on every challenge, not with a small probability. The tail numbers in Section 5 are heuristic estimates under H2, not probabilities over the challenge distribution. The same holds for the H-prime lists of the dependencies.

**Multipliers (proven prime).** There are 2^56 baby multipliers and 2^56 disjoint giant multipliers, each an H-prime with a deterministic primality proof. Fix the certified prime F = 162259276829214599630657442294509 (about 2^107.0), obtained by a chain of Pocklington certificates: p0 = 268435459 (prime by trial division), p1 = 36028803327197251 = 1 + 2 k p0 and F = 1 + 2 k' p1, each certified over its predecessor with predecessor^2 > successor (the script rebuilds and checks the chain). Multiplier candidates are N = 1 + 2kF for k in [2^103, 2^103 + 2^85), so 2^211 < N < 2^212 < F^2, and every N is above every element of S. Each candidate gets an H-test; each survivor gets a Pocklington test over F with bases 2, 3, 5, 7: if b^(N-1) = 1 mod N and gcd(b^((N-1)/F) - 1, N) = 1 for some base, every prime factor q of N satisfies q = 1 mod F, so q > F > sqrt(N) and N is prime. The test never accepts a composite, so every accepted multiplier is a proven prime and Lemma 4's use of primality is unconditional; a true prime is rejected only if all four bases are F-th power residues (probability about 2^-428), which merely discards it. Supply: primes occur in this progression with density about 2/ln N = 2^-6.2 and H-primes with 2^-20 of that (the same kind of density heuristic as H2), so the range holds about 2^58.8 H-primes against the 2^57 needed; the construction tests at most 2^66 survivors (mean 2^65) and reports fail if it has not found 2^57 multipliers, which like the list caps is a deterministic property of r resting on the density heuristic. Cost: 2^85 H-tests at 2^18 instructions plus at most 2^66 Pocklington tests at 2^24 instructions (up to 4 bases x 2 exponentiations of a 212-bit modulus): 2^92.03 units. The script also certifies the first three primes of the progression (k = 2^103 + 82, 206, 214) as a check of the procedure.

**Lemma 4 (distinct ratios, short interval).** Baby records are the points [a s]P0 and giant records are [a' s']A_d = [a' s' z]P0, with a in A, a' in A', and s, s' in S.

Suppose a s / (a' s') = b t / (b' t') as rationals, i.e. a s b' t' = b t a' s'.

1. The prime a divides the right side. It does not divide t or s', which are positive and smaller than a. Also a != a', since A and A' are disjoint. So a = b, and likewise a' = b'.
2. Then s t' = t s'. Write s = X0 + alpha, t' = X0 + delta, t = X0 + beta and s' = X0 + gamma, with alpha, beta, gamma, delta in [0, Delta).
3. Expanding gives X0 (alpha + delta - beta - gamma) = beta gamma - alpha delta. The right side has absolute value below Delta^2 < X0, so alpha + delta = beta + gamma and alpha delta = beta gamma, i.e. {alpha, delta} = {beta, gamma}.
4. So either (s, s') = (t, t'), the same pair, or s = s' and t = t', the diagonal.

So the N^2 ratios are pairwise distinct as fractions, except that for each multiplier pair (a, a') the 2L diagonal pairs share the ratio a/a'. That loses a fraction 1/(2L) = 2^-59.6 of the pairs.

**Lemma 4 is about rationals, not about residues mod r.** The products a s are about 2^421 > r. That the distinct fractions behave like independent uniform elements of H modulo r is heuristic H1, as in the dependencies. One quick check within one multiplier pair: a nontrivial mod-r coincidence needs r | (s t' - t s'). There are about 2^238 quadruples and the differences are about 2^310, so under any reasonable model the expected count of such coincidences per pair is negligible next to the 2^119 distinct ratios. This is a heuristic estimate, not a proof.

## 4. Algorithm B

**Chains.** Chain a holds the points [a u_j]Q, with Q = P0 on the baby side and Q = A_d on the giant side.
- At each twin step it emits records for [a v_j]Q and [a u_{j+1}]Q. The step is X_{j+1} = X_j + G_a[g_j], and the side is X_j - G_a[g_j].
- G_a[g] = [a g]Q for 1 <= g <= 2^44 is the chain's gap table: 96-byte affine entries, 2^50.58 bytes in all.
- The table is built by repeated addition of [a]Q in 2^16 strands, each started by one scalar multiplication. Strand 0 starts at [2a]Q, so every strand addition is [k a]Q + [a]Q with k >= 2. One extra scalar multiplication per chain computes the base [a]Q.
- The chain runs as 2^16 lockstep segments, each started by one scalar multiplication [a u_j mod r]Q at its stored first position. So each round is 2^16 independent twin steps sharing one Montgomery inversion.
- A fallback step is a plain reference addition with one record.

Each chain has L = 445418403372909336 = 2^58.6279 twin steps and emits 2L records. Per side there are 2^56 chains and N = 2 L 2^56 = 2^115.6279 records.

**Generating nodes.** There are 2^20 generating nodes. Each holds one gap table at a time (2^50.58 bytes, node memory 2^52) and runs 2^36 baby chains, then 2^36 giant chains. The gap list (6 bytes per step) is streamed to them from the 2^30 storage nodes in 512-byte blocks of 85 steps, by the protocol below: 4 instructions per step for reading and bookkeeping plus 144 per block of delivery (1.7 per step, listed 2).

**Key.** Field elements are held as blst holds them, as Montgomery residues m(y) = y R mod p with R = 2^384, and no conversion to the ordinary residue is made. The record for X = (x, y) is (key, id): key = bits 1..192 of the even one of m and p - m, where m is the stored residue of y fully reduced into [0, p), and id = the 64-bit chain index. Records are 32 bytes. By Lemma 2, [omega] leaves y and hence m unchanged, and [-1] maps y to p - y and hence m to p - m (y != 0 for every point of G1, since a point with y = 0 has order 2 and G1 has odd prime order r). So the six points [u]X share the key. The full reduction costs one 6-limb compare and conditional subtraction, listed at 12 instructions per record. That is a no-op, since blst's routines return fully reduced residues, and it is charged regardless. The map y -> m is a fixed bijection of F_p, so heuristic H3 is unchanged. Main and side outputs are keyed alike, and the toy uses the same Montgomery-residue key (with R = 2^32).

**Order, routing, matching, caps, candidates and recovery** follow the design of 0a45ddb2, with this N and with the communication specified below.
- **Routing.** The whole baby side runs first. There are 2^80 owner nodes indexed by key bits 0..79, read as two 40-bit fields (t, u), reached by two data-independent hops: each generating node scatters records into 2^40 send buffers by t and sends them to the relay (t,0); that relay scatters by u and sends to the owner (t,u). Records travel in 512-byte blocks of 16 records.
- **Baby insert.** At the owner, baby records go into a direct-mapped head array of 2^37 slots (key bits 80..116), chained through a pool of 2^37 32-byte records (capacity 2.59 x the mean of 2^35.628). A head slot (8 bytes) holds a 37-bit index of the first pool record and a 7-bit count; a pool record (32 bytes) holds key bits 117..191 (75 bits), the 64-bit chain id and a 37-bit link to the next pool record (176 bits used). A slot holds at most 64 baby records; a 65th makes the node fail.
- **Giant lookup.** Giant records read their head slot and compare their 75 remaining key bits with each chained baby record. Each comparison costs 8 instructions and counts against a per-node counter preset to 1.01 x the per-node mean.
- **Candidates.** An equal key appends a 40-byte candidate (full key, both chain ids) to a 2^16-slot list. After the giant side the status tree totals the candidates; above 2^44 the run fails.
- **Forwarding.** Each candidate is copied to both generating nodes named by the top 20 bits of its chain ids, each with a 2^24-slot receive area.
- **Re-walk.** Each generating node re-walks its own chain once, from start to end: the gap table is rebuilt, the strands, segments and base are restarted, and every twin step is recomputed with both outputs keyed and compared with the retained key. Each re-walk step is charged as a full twin step. The node keeps the outputs whose key equals the retained key, at most 4 per chain; a fifth makes the run fail.
- **Testing.** The giant chain's node sends its kept values a' s' in a 168-byte block to the baby chain's node. That node tests the six z = a s / (u a' s') for at most 16 pairs (at most 96 scalar multiplications) against A_d, and on a hit runs Lemma 1.

Nothing is compared against a whole other chain, and a wrong candidate never produces a wrong answer. Status collection is a fan-in-256 tree over all nodes, run in three passes.

**Communication (specified and priced).** All transfers use one mechanism: a sender's interface puts a block into a dedicated ring that the receiver owns for that sender. The links, their fixed sender sets and ring sizes:

| link | receiver | possible senders | ring per sender |
| --- | --- | --- | --- |
| hop 1 | relay (t,0) | the 2^20 generating nodes | 1 slot of 576 B |
| hop 2 | owner (t,u) | (t,0) only | 1 slot of 576 B |
| gap stream | generating node | the 2^30 storage nodes | 1 slot of 576 B |
| stream requests | storage node | the 2^20 generating nodes | 1 slot of 64 B |
| candidate return | relay (t,0), then the generating node | 2^40 senders at each stage | 1 slot of 64 B |
| kept values | generating node | the 2^20 generating nodes | 1 slot of 256 B |
| status tree | parent node | its at most 256 children; one parent for the downward pass | 1 slot of 64 B |
| list and chain setup | storage node, generating node | the list builder; the storage nodes | 1 slot of 576 B |

A 576-byte slot holds 16 records (512 B) and a 64-byte header carrying the record count (so partial blocks at the end of a side are well formed), a sequence number and an end-of-side flag. Relay roles exist on the 2^40 nodes (t,0). Every ring, credit word, arrival queue, waiter head and suspended entry is itemised in the memory table (Section 7).

**Protocol for one block on one link:**
1. **Credit check.** The sender keeps, per destination ring, a credit word that only the receiver writes (by an 8-byte put). It holds the consumed count: the receiver advances it only when the ring's single slot has been fully processed. With 1 slot the sender has at most 1 block not yet counted, and it never reuses the slot before its block is counted, so a suspended block is never overwritten. A sender reuses a slot (and its send buffer) only after it has dequeued and processed the matching (credit, destination) arrival entry, so at most one credit entry per destination is ever pending.
2. **Send.** If a slot is free, the sender writes a descriptor (destination node, ring, slot, local buffer address; 3 stores) and a doorbell store. Its interface copies the 576 bytes into that slot, then appends a tagged entry (data, ring) to the receiver's arrival queue, then raises an event. A credit put likewise appends a tagged entry (credit, destination) at the node that receives it. Each arrival queue is written only by its node's own interface, in arrival order, after the data is in place. So it has a single writer, and no entry is visible before its data. The queue has one entry for every possible unconsumed slot of every incoming ring, plus one credit entry per destination the node sends to (sized in Section 7), so it cannot overflow. Each arrival entry is 8 bytes: a 2-bit type (data, credit, message, status) and an index of up to 41 bits (the ring for data and messages, the destination for credits), so every destination or ring in a 2^40-element set is named directly.
3. **Send buffers.** On a record link (hops 1 and 2, gap stream, setup) a sender has 2 send buffers per destination: one being filled and at most 1 unconsumed at the receiver. On a message link (candidates, kept values, stream requests, status) it has one send slot and a credit word per destination. A buffer or slot is reused only after its credit entry has been processed (step 1). When a relay sends to its own node (relay (t,0) to owner (t,0)), the put is local: its own interface writes its own ring, with the same steps and price.
4. **Receive.** The receiver takes the next arrival-queue entry, reads the slot header and processes the records. When the slot's records are fully processed it writes its consumed count into the sender's credit word (an 8-byte put).
5. **Backpressure without spinning, with a per-record bound.** A relay may need to send while processing a received block. A record r of an input block B can be placed only if the filling buffer for its destination d has room. If that buffer is full and the other buffer of d is still unconsumed, B is suspended at r: a 16-byte suspended entry (ring index, 41 bits; record position, 4 bits; next-entry index, 41 bits) is linked at the tail of a FIFO waiter list whose head and tail indices (16 bytes) are kept per destination d, the ring's credit is not advanced (its single slot stays held), and the node continues with the next arrival entry. A (credit, d) arrival entry tells the node which waiter list to serve, with no scan. Processing it frees the consumed buffer of d; if d's filling buffer is full, that buffer is sent into the freed slot (it becomes the unconsumed one) and the freed buffer becomes the filling one. Then, while d's filling buffer has room for one record, the node resumes the head waiter, which places its current record r at once and continues; a resumed waiter that later meets another full destination d' is suspended again at that later record. So every suspension happens at a record that cannot be placed, and the resumption that ends it places that same record immediately. Each record therefore causes at most one suspend/resume episode per hop. Such an episode costs at most 8 instructions (entry write and link 3, unlink and state restore 5) and is charged per record per hop (Section 6), not per block. A generating node, a storage node, or an owner sending a candidate from its candidate list after the giant side, whose destination is full, waits for an event in the same way. Its episodes are likewise at most one per record or message; for messages (candidates, kept values, end flags) the 8 instructions per episode are inside the per-message charges (1024 against 704 listed per candidate, 512 against 400 per kept-value message, 216 against 144 per flushed buffer). Owners never suspend a received block, because inserting or looking up a record never sends. Events are latched by the interface, so an arrival between the empty-queue check and the wait instruction is not lost. The model has no spurious wakes. The 8 instructions charged per wake cover the re-check of the queue, the wait instruction and dispatch, and each put causes at most one wake.
6. **Progress.** Links are layered in a fixed order: storage/builder, generating node, hop-1 relay, owner. Separately, candidates go owner, hop-1 relay, generating node, and status goes from children to parent. The two directions never overlap in time: candidate forwarding starts only after the downward pass of status pass 2, when every record block of the giant side has been processed. A node only suspends a block of one layer while waiting on a strictly later layer. Owners and generating-node receive areas never suspend. So the wait-for relation is acyclic and every suspended block is eventually resumed. A node's other roles keep making progress while one role has suspended blocks, because blocks are taken from the arrival queue out of order and a suspension holds only its own ring's slot. The number of wakes and credit puts is at most a constant per delivered block, and the number of suspend/resume episodes at most one per record per hop, independent of timing.
7. **End of a side.** When a sender has produced everything for a side, it puts one final block, possibly empty, with the end flag set, into each of its destination rings. A relay flushes its own buffers, each with the end flag, only after it has processed end-flagged blocks from all its senders and has no suspended block. Status pass 1 is the global barrier between the baby side and the giant side: an owner reports once the end-flagged block from its single hop-2 sender is processed, and the generating nodes start the giant chains when the downward pass reaches them. So no giant record reaches an owner before all baby records are inserted. The flushes are charged in the per-node fixed work at 144 + 72 per buffer, 2^40 buffers per relay role and per generating node per side.

**Price per block per link, listed:**

| item | instructions |
| --- | ---: |
| sender: credit check and send-buffer reclaim | 6 |
| sender: descriptor (3 stores), doorbell, addressing | 8 |
| interface copy of the 576-byte slot, charged as 36 load-pair + 36 store-pair instructions plus loop control | 80 |
| interface: arrival-queue append and event at the receiver | 3 |
| receiver: wake on the data put (at most one per put) | 8 |
| receiver: arrival-queue read, slot address, record count from the header | 7 |
| receiver: credit return, an 8-byte put (descriptor, doorbell, interface copy charged) | 10 |
| sender: wake on the credit put (at most one per put) | 8 |
| sender: credit update and resume check | 6 |
| **total** | **136** |

This is charged at 144 per block per link, i.e. 9 per record per hop and 18 per record for the two hops. Suspend/resume episodes are charged separately, at 8 per record per hop (16 per record), by the bound in step 5. The interface's byte movement is not processor work under the cost model, but it is charged as if it were.

Other links:
- **Messages** (candidates, kept values, stream requests, status records) are charged a full 144 per message per link. The status tree's 2^12 instructions per node per pass are an amortised figure: the tree has fewer edges than nodes, and each edge costs 2 x 144 (one message up, one down) plus combining one 48-byte record, so a parent with 256 children spends about 2^15 while the per-node average stays below 2^12.
- **End-of-side flushes.** Each relay flushes its 2^40 partial buffers at the end of each side, at 144 + 72 each; the generating nodes do the same. These are in the per-node fixed work (2^50 instructions per relay or generating node).
- **Data distribution.** The list goes from its builder to the storage nodes (2^61.2 bytes at 144 per block, charged 2^62 instructions). Each chain's multiplier, 2^16 segment-start values and Q go to its generating node (charged 2^20 instructions per chain). Both are in the ledger row "data distribution".

## 5. Sizes, caps and success

L is the least integer with 6 N^2 >= 0.4944 m for N = 2 L 2^56. This gives N = 2^115.6279 and K = 6 N^2 / m = 0.494400. With the diagonal loss, K' = K (1 - 1/(2L)). Under H1 the probability that z is among the covered classes is 1 - exp(-K') = 0.390063.

The caps below have the same structure and derivations as in 0a45ddb2, recomputed for this N. Every tail is rounded towards the larger probability, and the charged work is the maximum under the caps.

- **Comparisons.**
  - 2^117 head slots at load lambda = N / 2^117 = 0.386. Mean comparisons are N^2 / 2^117 = 2^114.256, i.e. 2^34.256 per node, with a node cap of 1.01 x that (2^34.270).
  - A slot over 64 records: Poisson tail 2^-274.7 over all slots.
  - A node's baby or giant count more than 2^-10 off its mean: 2^-24270 over all nodes.
  - Given these, Bernstein (M = 64, variance margin 4) for a 0.8 % excess: 2^-165602 over all nodes.
- **Matches.**
  - False matches (equal 192-bit key, different mu6 orbit) have mean N^2 / 2^192 = 2^39.26 against a cap of 2^44; Chebyshev tail 2^-47.6.
  - True matches are Poisson with mean K under H1, with a budget of 2^20; exceeding it has probability below 2^-20524381.
  - Candidates per owner node have mean 2^-40.7 against 2^16 slots; overflow below 2^-2670126.
  - Candidates per generating node have mean 2^20.26 against 2^24 slots; overflow below 2^-40417839.
- **Re-walk key repeats.** One key on more than 4 of a chain's 2L outputs: below 2^-419 over all chains.
- **Node pool overflow.** Chernoff 2^-(5.4 x 10^10) over all 2^80 nodes; a budget of 2^-20 covers modelling error.
- **List construction.** These are heuristic estimates under H2, not probabilities over challenges (see Section 3). Span over 2^100: Chernoff bound for a sum of L exponentials exceeding 2.59 times its mean, below 2^-4e17. Fallback steps over 2^40: below 2^-3.38e12.

Success claimed: 0.390063 - 2^-20 - 2^-40 = 0.390062 >= 0.390 (asserted), the required 0.39. The 2^-40 covers every cap and re-walk tail above: the script asserts each is below 2^-47 and their sum below 2^-40 (it is below 2^-47.6). The 2^-20 covers node-overflow modelling (asserted below 2^-20).

## 6. Instruction listing per twin step (AArch64 scalar, on top of the one 2000-instruction reference addition)

**Walk overhead, once per step (52 instructions):**

| item | instructions |
| --- | ---: |
| gap stream read and pointer advance | 5 |
| gap stream transport (6 bytes per step, blocked) | 4 |
| gap stream delivery (144 per block of 85 steps, 1.7) | 2 |
| stream request (144 per 85-step block, 1.7) | 2 |
| table address | 2 |
| operand loads and stores not counted inside the reference price | 18 |
| batch-inversion prefix and suffix loads and stores | 6 |
| lane control and loop | 8 |
| flow control and counters | 5 |

**Side output arithmetic, once per step (1099 instructions: 1093 measured fused kernel + 6):** see Section 2.

**Record handling, once per record and twice per step.** Baby records cost 121 each:

| item | instructions |
| --- | ---: |
| parity select and key extraction from the stored Montgomery residue | 23 |
| full reduction of the residue to [0, p) (6-limb compare, conditional subtract) | 12 |
| scatter at 2 hops (12 per hop) | 24 |
| block bookkeeping (header, buffer switch) amortised over 16 records, at 2 hops | 6 |
| 32-byte transport copy at 2 hops | 12 |
| insert | 10 |
| delivery at 2 hops (144 per block per link, 16 records per block) | 18 |
| suspend/resume, at most one episode per record per hop (step 5), at 2 hops | 16 |

Giant records cost 117 each: the same items without the insert, plus a lookup fixed cost of 6. Comparisons are charged separately, at 8 each, up to the node cap.

**Totals.** Listed per twin step: baby 52 + 1099 + 2 x 121 = 1393, giant 52 + 1099 + 2 x 117 = 1385. Charged: 1436 per twin step on both sides, which is 1.0309 x the larger listing (43 instructions of margin for unlisted items, as in twin-ap-v2, against 39 in twin-ap-v1). That makes 3436 instructions per twin step and 1718 per record. Gap-table additions are charged 2000 + 52, and re-walk steps are charged as full twin steps (3436).

## 7. Ledger (units of 2000 instructions; W = ceil of the instruction sum / 2000)

| term | log2 units |
| --- | ---: |
| records: 2^57 chains x L twin steps x (2000 + 1436), two records per step | 116.409 |
| per-node fixed work: (2^80 + 2^30) x 2^37 + 2^40 relay roles x 2^50 + 2^20 x 2^50 instructions | 106.034 |
| twin list S: (2^100 + 2^44 + 2^81 + 2^84) H-tests at 2^18 instructions, window duplicate checks included | 107.034 |
| comparisons: 2^80 nodes x cap x 8 | 106.304 |
| candidates at cap 2^44: re-walk both chains (L + 2^40 steps at 3436), rebuild both tables, 4 x 2^16 + 2 restarts (strands, segments, bases) | 104.409 |
| gap tables: 2^57 chains x 2^44 additions at 2000 + 52 | 101.037 |
| fallback reserve: 2^40 steps per chain x 2^57 chains at 3436 | 97.781 |
| multipliers: 2^85 H-tests + at most 2^66 Pocklington tests | 92.034 |
| status tree: 3 passes x every node x 2^12 instructions (its messages included) | 82.619 |
| table strand starts: 2^57 x 2^16 scalar multiplications at 2^12 units | 85.000 |
| segment starts: 2^57 x 2^16 scalar multiplications | 85.000 |
| chain bases [a]Q: 2^57 scalar multiplications | 69.000 |
| data distribution: list to storage nodes (2^62 instructions), chain setup to generating nodes (2^57 x 2^20) | 66.034 |
| candidates at cap: up to 6 x 4 x 4 = 96 scalar multiplications each | 62.585 |
| Lemma 1 up to 2^20 + 1 times, input access, output | 44.000 |
| candidates at cap: copy into two messages (128) and 2 x 2 return-path deliveries at 144, charged 1024 instructions | 43.034 |
| candidates at cap: kept values to the baby chain's node, 256 + one delivery (144), charged 512 instructions | 42.034 |
| **W** | **116.413621350** |

W = 110660322149086336191977377548476613 units, log2 W = 116.413621349862613... (the 60-digit decimal check is in the script). Score = ceil(10^6 log2 W) / 10^6 = 116.413622.

For comparison with the dependency: its record row is 2N x (2000 + 180) at N = 2^115.6644, with K = 0.52. The twin step, with the communication specified and charged, moves the per-record price from 2180 to 1718 instructions, which is log2(1718/2180) = -0.344. K = 0.4944 accounts for another -0.036. Against twin-ap-v1 (116.544103): log2(3436/3760) = -0.1300 from the measured fused side output and -0.0009 from K; against twin-ap-v2 (116.509163): log2(3436/3672) = -0.0958.

Preprocessing is all work that runs before the first giant chain and does not read A_d: the baby records (L twin steps on each of 2^56 chains at 3436 instructions), the baby half of the fallback reserve (2^40 steps per baby chain at 3436), the baby gap tables (2^56 x 2^44 additions at 2052), the baby strand, segment and base starts (2^56 x (2 x 2^16 + 1) scalar multiplications), all per-node fixed work, the data distribution, the twin list, the multipliers and the first of the three status-tree passes. The exact sum is 55391222379049707710905158411850089 units = 2^115.415213...; the claim 2^115.4153 is rounded up and the script asserts pre^10000 <= 2^1154153 in integers. It is included in W (every term is a row or part of a row of the ledger).

Memory, summed as integers (owner and relay roles at their itemised sums, generating and storage nodes at their full allocations):

| nodes | bytes each | itemised allocation |
| --- | --- | --- |
| 2^80 owner nodes | 2^42.3222 (fixed setup: zeroing the 2^40-byte head array with 4 store-pairs per 64-byte loop iteration, 2^36.81 instructions, within the 2^37 charged) | heads 2^37 x 8 B, pool 2^37 x 32 B, hop-2 ring 576 B, arrival queue and credit words 2^12 B, candidate send buffer 2^10 B, status-tree rings 257 x 64 B, candidate list 2^16 x 40 B, code and counters 2^30 |
| 2^40 relay roles (extra, on the nodes (t,0)) | 2^50.3129 | send buffers 2^40 x 2 x 576 B, credit words 2^40 x 8 B, hop-1 rings 2^20 x 576 B, candidate return rings 2^40 x 64 B, arrival queue (2^21 + 2^41) x 8 B, suspended entries (2^20 + 2^40) x 16 B, waiter-list head and tail 2^40 x 16 B, credit word and send slot for each of the 2^20 generating nodes (72 B) |
| 2^20 generating nodes | 2^52 (full allocation; itemised 2^51.4433) | gap table 2^44 x 96 B, send buffers 2^40 x 2 x 576 B, credit words and waiter flags 2^40 x 16 B, gap-stream rings 2^30 x 576 B, arrival queue (2^31 + 2^41 + 2^21 + 257) x 8 B, message send state for 2^30 request and 2^20 kept-value destinations at 72 B, status-tree rings 257 x 64 B (chain setup reuses the gap-stream rings), candidate return rings 2^40 x 64 B, kept-value rings 2^20 x 256 B, stream and walk buffers 2^36, candidate receive area 2^24 x 40 B |
| 2^30 storage nodes | 2^40 (full allocation) | list and multipliers 2^62.67 B in all; per node, send buffers 2^20 x 2 x 576 B, request rings and credits 2^20 x 72 B, arrival queue 2 x 2^20 x 8 B, status-tree rings, the ring from the list builder |

The total is 6647441250573927321040555141874319360 B = 2^122.322211... The claim is 2^122.3223 and the script asserts mem^10000 <= 2^1223223. The records alone are 2^120.63 B. Every node's memory, the relay node (t,0) included (2^42.3 + 2^50.3 bytes), is below 2^52 bytes.

## 8. Heuristics

**H1 (coverage).** The 6 N^2 (1 - 1/(2L)) distinct fractions a s / (u a' s') behave like independent uniform elements of H for the purpose of hitting z. This is the same heuristic the dependencies use, applied to a different but equally structured set: products of a large prime and an element of a short interval. Lemma 4 removes the exact repeats. The toy runs the identical construction at d = 2^7 and agrees with 1 - exp(-K') within 0.7 standard errors at 400 trials (0.3625 against 0.3762, z = -0.56) and at 2000 trials (0.378, z = +0.17). The toy cannot separate 0.390 from 0.39; the production figure rests on H1.

**H2 (twin density).** For an element u of the list and nearby integers c, the events "c in H" and "2u - c in H" behave like independent events of probability 2^-20 each. This is a statement about the 2^20-th power residue character on pairs of nearby integers. It is of the same kind as H2 of the dependencies (H-prime density and exponential gaps).

The list rule itself adds a small bias. A side 2u_j - c lies below u_j, in integers that an earlier step already scanned and rejected as candidates. Rejection lowers their U-probability from q = 2^-20 to at most q/(1 + q). The mean gap therefore exceeds 1/q^2 by at most a factor 1 + q: 1 + 2^-20 in production, but 1.125 at level 2^3 and 1.0625 at 2^4.

Mode `h2` gives two kinds of evidence on the real r:

(a) **The exact list rule** applied to the integers 2^210 + i, i < 24000, with H replaced by the 2^k-th powers modulo r:

| level | mean gap / 4^k | prediction (1 + 2^-k) | twin steps | coefficient of variation |
| --- | ---: | ---: | ---: | ---: |
| 2^3 | 1.228 (standard error 0.061) | 1.125 | 281 | 1.03 |
| 2^4 | 1.054 (standard error 0.122) | 1.0625 | 59 | 0.94 |

Geometric gaps have a coefficient of variation of about 1. There were no fallbacks.

(b) **Direct sampling** of u uniform in [2^210, 2^210 + 2^100) and g uniform in [1, 2^44], with the low 8 bits of the 20-bit character logarithms of u, u + g and u - g computed exactly (Pohlig-Hellman in mu_{2^8}). With 4000 samples (hosted input):
- chi-square of (e+ mod 16, e- mod 16) over 256 cells: 263.3 (expected 255, sd 22.6);
- chi-square of (e0, e+, e-) mod 4 over 64 cells: 50.9 (expected 63, sd 11.2);
- P(e+ = e- = 0 mod 2^k) / 4^-k = 0.989, 0.960, 0.848, 0.960 for k = 1..4 (counts 989, 240, 53, 15).

A local run with 30000 samples (seed 7, obtained by calling `h2_mode(pairs=30000, seed=7)` directly; the declared input schema runs the 4000-sample default) gave chi-squares 258.4 and 65.3 and ratios 0.982, 0.991, 1.082, 1.118 (counts 7367, 1859, 507, 131).

The top 12 bits of the character logs, which matter for the full 2^-40 density at level 2^20, are not sampled. Their independence is the heuristic itself. The toy curve list (d = 2^7) has a mean gap / d^2 of 0.855 over 13 steps.

**H4 (multiplier supply).** The progression N = 1 + 2kF, k in [2^103, 2^103 + 2^85), contains primes with density about 2F/((F-1) ln N) and H-primes with 2^-20 of that, so about 2^58.8 H-primes against the 2^57 needed. This is a prime-density heuristic for one arithmetic progression of modulus about 2^108, of the same kind as the dependencies' H-prime density assumption. Primality itself is proven (Pocklington); only the supply is heuristic, and a shortfall would make the construction report failure, deterministically for every challenge, as with H2. The script certifies the first three primes of the progression.

**H3 (key bits).** The 192-bit keys behave like independent uniform strings across mu6 orbits: occupancies are binomial and cross-orbit equality has probability 2^-192. This is the same as in the dependencies. The toy chi-square of pseudo-node occupancy (Montgomery-residue key) is 31.55 ± 0.40 at 400 trials (hosted input) and 31.08 ± 0.18 at 2000 trials, against 31.

## 9. Evidence (`evidence.py`, standard library, about 7 s for mode all on one core; `measure/`, C and shell)

**`measure/`** holds the side-output measurement of Section 2: `side_fused.S` (the charged fused kernel), `gen_side_fused.py` (regenerates `side_fused.S` from blst's `build/elf/mul_mont_384-armv8.S`, showing the routines are verbatim copies; `python3 gen_side_fused.py <blst checkout>`), `side_fused_test.c` (its count and self-checks), `side_output.c` (the blst field-layer kernel of twin-ap-v2, three call variants), `fp_ops.c` (single exported Fp calls), `run.sh` (build and count inside the organizers' image, same flags as `scripts/reference_costs/run.sh`) and `reproduce.sh` (host driver: clones blst at the pinned commit, builds `scripts/reference_costs/Dockerfile` for linux/arm64 and runs `run.sh`; `bash measure/reproduce.sh <repository root>`). It is not a hosted experiment (the hosted runner executes Python only); its output is summarised in Section 2 and kept verbatim in `measure/output.txt`, and `evidence.py` carries the measured figure as a constant and asserts that the charge covers the listing.


**`ledger`** computes every number above. W, the record counts and memory are exact integers; the tails, K and success are floating point. It asserts d | r-1, 6d | r-1, mu6 in H, Delta^2 < X0, that the multipliers lie above the list, the listed versus charged per-step instructions, every tail, and the node memory sums with the exact memory inequality. It also runs a 60-digit check that log2 W does not exceed the claimed score.

**`h2`** produces the two experiments quoted in H2.

**`ec`** runs Algorithm B end to end on the j = 0 curve y^2 = x^3 + 3 over F_207670579, of prime order r = 207671809 (d = 128, m = 1622436):
- X0 = 2^40 and T = 2^18. The twin list is built by `build_twin_list`, the rule of Section 3 including the window duplicate check: 13 twin steps, no fallback. Delta^2 < X0 and distinctness of S modulo r are asserted.
- The multipliers are Pocklington-certified H-primes N = 1 + 2kF_t above the list, with F_t a prime certified over a trial-division prime (the same `certified_prime_chain` and `pocklington` functions as production). There are 14 chains per side with 26 records per chain (n = 364, K = 0.490).
- Each side runs in lockstep with batch inversion, every twin step computing P + G and P - G from one shared lane inverse (`batch_twin`). On the first trial, 52 twin outputs are asserted equal to independently computed scalar multiples [a s]P0.
- Then: the 16-bit parity key taken from the Montgomery residue of y (R = 2^32), 32 pseudo-nodes x 16 heads, the bucket table, candidates (key, id, id'), one linear re-walk of each candidate chain recomputing both outputs of every twin step, at most 4 kept outputs per chain, six z per pair tested against A_d, 7 square roots and the root-of-unity scan.

Differences from production:
- The gap table holds only the 13 gaps of the toy list, built from doublings, instead of all g <= T.
- The caps other than the re-walk cap are not enforced; the counts are compared with the model instead.
- It runs in one process, so nothing is routed.
- Lemma 1's root-of-unity search is a linear scan over the 128 roots in batches of 16 instead of the 2^10 x 2^10 sort-and-lookup.
- Each candidate re-walks both of its chains in the same process (production: each generating node re-walks its own chain once per candidate); the ledger charges both re-walks per candidate either way.

Results with seed 2026:

| run | success | predicted (diagonal-corrected) | z | comparisons (exp. 258.8) | false matches (exp. 2.022) | true matches |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 400 trials (hosted input) | 0.3625 | 0.3762 | -0.56 | 259.5 | 1.960 | 0.443 |
| 2000 trials (local) | 0.378 | 0.3762 | +0.17 | 259.3 | 2.000 | 0.473 |

There were no re-walk cap failures, and every recovered s was correct.

## 10. Limitations

- **No full-size run.** This is an analytic bound.
- **Side-output price.** The side output is a measured kernel (1093 instructions, Section 2), counted with the organizers' image, compiler, flags, blst commit and Callgrind differencing; the count is an instruction count, not a timing. The kernel is hand-scheduled assembly built from verbatim blst routines; its correctness rests on the self-checks (lanes, extreme operands, curve points) and on Lemma 2a's bound, which the checks also exercise, not on a formal verification. Production code would hold the lane state in its own layout; the measured kernel reads all five operands from memory and stores its result, so a different layout needs no instructions beyond those already charged: the kernel's own 15 call and loop instructions are charged on top of the walk overhead's 8 (lane control) and 18 (operand loads and stores). The margin is 43 instructions per step. The measurement ran under QEMU user-mode emulation of linux/arm64 (Callgrind counts guest instructions; the organizers' totals reproduce exactly).
- **Pricing rule.** The pricing relies on the cost model's rule that one batched incomplete affine addition costs one unit. The side output shares only that addition's lane inverse and is charged in instructions on top.
- **Heuristics.** H1, H2, H3 and H4 are heuristics; H1-H3 are validated at toy scale. H2 is also checked on the real r at power levels up to 2^4 and on the low 8 bits of the character logarithm. The list and multiplier constructions are deterministic: a failure of H2 would make the attack fail on every challenge, not with a small probability.
- **Memory.** Memory is exponential (2^122.3223 B). The generating nodes need 2^51.44 bytes each for the 2^44-entry gap table and 2^40 send buffers, within a 52-bit address space. A smaller table (T = 2^42) raises the fallback rate to e^-4 per step and the score by about 0.013.
- **Communication.** The protocol of Section 4 assumes a network interface that performs one-sided puts and raises an event; its byte movement is charged as instructions anyway. Routing and delivery are 76 listed instructions per record (scatter 24, block bookkeeping 6, copy 12, delivery 18, suspend/resume 16), 152 of the 1436 charged per twin step. These, the key extraction and the walk overhead remain itemised estimates, unchanged from twin-ap-v1. A reviewer who prices communication higher by x instructions per record raises the score by about log2(1 + 2x/3436).

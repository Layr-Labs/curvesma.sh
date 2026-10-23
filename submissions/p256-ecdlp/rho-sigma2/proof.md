# p256-ecdlp: parallel Pollard rho on negation classes with a sigma walk

Claim: for uniform x in [0,r) and Q=[x]P on NIST P-256, the algorithm below returns x
with probability at least 0.39. It uses at most **2^127.273624** target-group-operation
equivalents in total (score 15909203/125000), including all preprocessing, failures,
lookup and verification. Peak memory is at most 2^75.542441 bytes. Preprocessing is 2^49.619425 units, and no nonuniform advice is
used. This is a classical cold attack on one target.

The success bound rests on the disclosed heuristic H1, the standard random-mapping model
for r-adding walks on negation classes, and on the rate bound H2. The cost rests on H3,
instruction counts measured with the organizer's toolchain. **There is no reference
operation on this track:** every instruction, including all P-256 field arithmetic and
point additions, is charged at 1/C with the borrowed normalisation C = 1000. This is a
carefully accounted generic attack. It is not a practical break.

Files:
- `accounting.py` recomputes every number with exact integers (output in
  `evidence/accounting.json`).
- `measurement/` holds the per-step instruction-count harness and the field code
  (`p256_field.h`, `p256_sub.h`), plus the Python replica `walk_ref.py`.
- `evidence/measurement/` holds the counts on two input sets.
- `experiment.py` is the declared toy-curve experiment for the identical walk;
  `evidence/small-curves/` holds larger offline runs of the same file.

## 1. Algorithm

### 1.1 Notation and parameters

p, a = −3, b, r and P come from the protected profile. The cofactor is 1, and r is prime.
M(v) = v·2^256 mod p is the Montgomery representative, a fixed bijection of Fp onto
[0,p). Every comparison and bit test below is on M-values, read as integers.

Negation splits the nonzero points into (r−1)/2 classes {W, −W}. Write N = (r−1)/2 + 1.

Parameters:
- A = 2^16 adders.
- K = 2^32 walks: 2^16 processors × 2^16-walk batches.
- theta = 2^−60.
- L = 57695294864490661436112010649 ≈ 2^95.542.
- CAP = 64/theta.
- DP_MAX = 2·theta·S_max + 2K, where S_max = K·(L + CAP).

### 1.2 Classes, keys and the sigma walk

For W = (x, y) ≠ O:
- **Class key.** key(W) = M(x). It is constant on {W, −W} and determines the class.
- **Canonical point.** canon(W) = (x, (−1)^s·y), where s = M(y) mod 2. p is odd and y ≠ 0
  (no 2-torsion), so M(−y) = p − M(y) has the other parity, and canon picks one point of
  each class.
- **sigma.** sigma_W = [(−1)^s], so that W = sigma_W(canon(W)).
- **Walk step.** i(W) = key(W) mod 2^16, and

      f(W) = W + sigma_W(R_i(W)) = sigma_W(canon(W) + R_i).

  So the class of f(W) is the class of canon(W) + R_i, which is exactly the standard
  negation-class walk. The kernel reads sigma_W(R_i) from a table of R_i and −R_i, so no
  point is ever negated, and sigma costs one bit test.
- **Distinguished point.** W is distinguished when bits 16..75 of key(W) are zero.

### 1.3 Procedure

1. **Validate Q.** Check it is on the curve. If Q = O, output 0.
2. **Adders.** For each adder i < 2^16, draw c_i, d_i uniform in [0,r), set
   R_i = [c_i]P + [d_i]Q, and tabulate R_i and −R_i. For each walk w, draw a_w, b_w and set
   W_w = [a_w]P + [b_w]Q. All coins are stored, and coefficients are not tracked during
   the walk.
3. **Step.** W ← f(W), as one batched incomplete affine addition. Montgomery's
   simultaneous inversion over the 2^16 walks of a batch costs 3M per walk; then
   lambda = dy·inv (M), lambda^2 (S), y3 = lambda·(x1 − x3) (M). One Fermat inversion is
   done per batch.
   - If the batch product is zero (W.x = R.x for some walk, at most about 2^−239 per
     step), that walk uses the unbatched complete formula.
4. **Fruitless 2-cycles.** Each walk keeps the low 64 bits of its previous key.
   - On a match, the walk evaluates f(new) unbatched.
   - If key(f(new)) = key(old) exactly, then {old, new} is a 2-cycle of the class walk.
     The walk replaces new by 2·canon(M), where M is the element with the smaller key.
5. **Short cycles (3..64).** Every 64 steps each walk records its key as an anchor.
   - On a later match, the walk iterates f unbatched up to 64 times.
   - If its class returns, the walk escapes from 2·canon(M), with M the cycle element of
     least key.
6. **Distinguished points.** The walk sends (key, w, step index) to a hash table keyed by
   key.
   - A walk runs L steps and then continues until its next distinguished point.
   - A walk that exceeds CAP steps without one is abandoned and restarted from fresh coins,
     with at most K restarts in total.
   - The run fails if more than DP_MAX points arrive or the restarts run out.
7. **Resolution.** On the first key stored twice, both walks are replayed from their coins
   to the recorded step indices.
   - The replay tracks (a, b) mod r: an image sigma(R_i) has coefficients multiplied by
     (−1)^s, and escapes double canon(M).
   - It obtains canon(W) = [a]P + [b]Q for both walks. These are equal, since the keys are
     equal.
   - If b1 ≠ b2, it outputs x = (a2 − a1)/(b1 − b2) mod r after checking [x]P = Q.
     Otherwise it continues; the bound counts this case as a failure.

**Merged trails stay merged**, for the same reasons as in the order-6 case. Outside
fruitless cycles the next class depends only on the current class. Escapes fire only after
an exact check, and they exit from a point determined by the cycle's set of classes. The
64-bit filters only decide when an exact check runs.

## 2. Success probability

Let T_steps = K·L. Then tau = T_steps^2/(2N) = 0.530300.

**H1 (heuristic, disclosed).** The class map [W] → [canon(W) + R_h([W])], with uniform
random adders, behaves like a random mapping on the N classes, apart from the r-adding
penalty. That penalty is a factor of about (1 − 1/(2A)) in the birthday exponent
(Bos–Kleinjung–Lenstra 2010, eq. (2); Teske 2001; Bernstein–Lange 2012). Under H1, the
probability that the first T_steps steps contain no useful class collision is at most
exp(−tau·(1 − eps)). Here eps = 4/A + 2·(3/A) + 2·(192/A^2) = 0.00015268. Its three
margins are as in the order-6 package:
- the adder penalty is taken as 4/A, eight times 1/(2A);
- fruitless 2-cycles enter at a rate of at most 1/A, with 3 non-fresh states each;
- short cycles enter at a rate of at most 1/A^2, with 192 wasted states each.

**H2 (rates).**
- 2-cycles are entered at a rate of at most 1/A per step. The literature value is 1/(2A)
  (Bos–Kleinjung–Lenstra; Bernstein–Lange–Schwabe).
- Cycles of length 3..64 are entered at a rate of at most 1/A^2. 4-cycles occur at about
  1/(4A^2).
- Longer cycles with no distinguished point arise at a rate of at most A^−32. Abandonment
  handles them and is charged.

The remaining failure terms are:
- detection beyond CAP: about 2^−25;
- degenerate b1 = b2: at most 2^−250;
- DP_MAX or restart overflow: under 2^−1000 by Chernoff;
- exceptional additions: under 2^−110.

Together these are below delta = 2^−20. Therefore:

**Pr[success] ≥ 1 − exp(−0.5303·(1 − 0.00015268)) − 2^−20 = 0.411523 ≥ 0.39.**

The budget is 1.0299·√N steps, against 0.9943·√N for success exactly 0.39. The 0.0215
of slack in the model bound covers the sampling error of the toy evidence below.

**Evidence for H1 and H2.** `experiment.py` with `"curve": "a3"` runs this exact sigma
walk on random prime-order toy curves y^2 = x^3 − 3x + b. It uses the same pieces as
production:
- the key M(x) with M(v) = v·2^64 mod p;
- the sign from the parity of M(y), and the adder sign applied to the adder;
- both cycle rules and escapes;
- distinguished points and abandonment.

It counts every step, including fruitless, confirming and escape steps, up to the first
useful class collision. It then solves end to end and checks [x]P = Q.

<!-- toy-evidence -->
Offline runs (`evidence/small-curves/`, pooled in `summary.json`): 3,300 trials on 3 curves
with N = 2^30.1, 2^32.9, 2^32.4 classes, adder tables of 2^11 and 2^16.
- Every trial solved and verified [x]P = Q.
- Success by T = sqrt(2 tau N) at the claimed tau = 0.5303: pooled 0.4255, Wilson 95%
  interval [0.409, 0.442]. The model value is 0.412.
- At tau = 0.4943 (model 0.390) the pooled success is 0.4039, interval [0.387, 0.421].
- Mean steps to the first useful collision: 1.231 sqrt(N), against sqrt(pi/2) = 1.253.
- 2-cycle entry rate per step:
  - 2.48e-04, 2.44e-04 with A = 2^11 (1/(2A) = 2.44e-04);
  - 7.76e-06 with A = 2^16 (1/(2A) = 7.63e-06).
- Longer fruitless cycles: 4 in total. Abandoned walks: 3.
<!-- /toy-evidence -->

The declared initial input (`curve` a3, `bits` 26, 150 trials, seed 2601) is a quick
reproducibility check, not the statistical evidence. At tau = 0.5303 it reports 57/150 =
0.380, Wilson 95% interval [0.306, 0.460], which contains the model value 0.412; with 150
trials the standard error is about 0.04. The success claim rests on the 3,300 pooled
offline trials above, which any input within the schema can spot-check.

Limits: the toy curves have N up to about 2^33, against 2^255 in production. H1 is an
extrapolation in N, supported by the size-free random-mapping argument and by the
literature. Examples are the negation-map runs in Bernstein–Lange–Schwabe 2011, which
measured 1.0074 and 0.9996 times the ideal on 55- and 60-bit curves, and
Bos–Costello–Miele 2014, which measured 0.9907 fruitful steps on P-256 itself. H1 is not
proved for this walk.

## 3. Computation

**No reference operation.** COST_MODEL.md borrows C = 1000 from secp256k1 as a
normalisation and requires all P-256 arithmetic, including point additions, to be charged
as instructions. Each step is therefore charged its full measured instruction count
divided by 1000.

**Field arithmetic.** `measurement/p256_field.h` implements Montgomery multiplication and
squaring with R = 2^256 as GCC inline AArch64 assembly in general-purpose registers.
- Because p ≡ −1 mod 2^64, the Montgomery factor is 1, and each of the four reduction
  rounds is a handful of shift and add instructions with no multiplications. This is the
  approach of OpenSSL's ecp_nistz256 (Gueron–Krasnov), re-derived for this harness with
  shifted-register operands and fewer constants.
- Products are 16 `mul` + 16 `umulh`; outputs are fully reduced.
- `p256_sub.h` holds branchy subtractions that add p back only on a borrow, and a fused
  a − b − c.
- Standalone costs in the same toolchain, per call including about 5 instructions of loop
  overhead: mul 123, sqr 98, sub 26. Inside the kernel they are inlined.
- Correctness: the replica check below compares every final point of 2^16 walks after
  2^19 steps with Python big-integer arithmetic, on two input sets. The field routines
  were also checked on about 150,000 random and edge-case operands each.

Inline assembly is ordinary scalar AArch64 code. Callgrind counts each executed
instruction exactly as for compiled C, which is how the organizer measured blst's assembly.

**Per-step measurement (H3).** `measurement/step.c` implements steps 3–6 for a batch of
2^16 walks with 2^16 adders.
- It was compiled with `-O3 -fno-tree-vectorize -mgeneral-regs-only` in the organizer's
  Docker image (GCC 14.2.0, Valgrind 3.24.0), and counted with Callgrind `Ir`.
- Counts are differenced between 8 and 4 batch-steps (2^18 steps).
- Every instruction inside the measured loop is charged.

| Input set | total Ir/step | 2-cycles in 2^19 steps |
| ---: | ---: | ---: |
| seed 1 | 816.933 | 2 |
| seed 2 | 817.567 | 3 |

This covers:
- 5 multiplications and 1 squaring;
- 5 subtractions, one of them fused;
- key, parity and adder-image selection;
- the 2-cycle, anchor and distinguished-point filters, and the 2-cycle confirmations and
  escapes that occurred;
- the batch's prefix stores and loads;
- one Fermat inversion per 2^16 walks (about 38,000 Ir, or 0.6 per step);
- loop control.

**830 instructions** are charged per step, so each step costs **0.83 units**.
`measurement/walk_ref.py check` replays all 2^16 walks in Python integers and matches the C
final points exactly on both input sets.

Rare paths are charged within the 12.4-instruction margin:
- **2-cycles.** At most 1/A per step × 3 unbatched operations × 2^16 = 3 instructions per
  step.
- **Short cycles.** Below 2^−20 instructions per step.

A separately measured unbatched step (affine add with its own Fermat inversion, plus
adder selection) costs 40,871 Ir on both sets. Every unbatched operation below is charged
2^16.

**Total.** S_max = K·(L + CAP) ≤ 2^127.542441 steps.

| Component | Bound (units) | log2 |
| --- | --- | ---: |
| walk steps | 0.83·S_max | 127.273624 |
| adder table | 2^16 adders × 2 scalar mults × 768 unbatched ops × 2^16 / C | 32.62 |
| walk starts and restarts | 2K × 2 scalar mults × 768 × 2^16 / C | 49.62 |
| distinguished points | DP_MAX = 2^68.54 entries × 2^20 instructions / C | 78.58 |
| collision replay | 2·(L + CAP) unbatched steps × (2^16 + 2^10) / C | 102.60 |
| solve, verify, output | 2 scalar mults + 2^24 instructions | 16.84 |
| input validation | 2^24 instructions | 14.03 |

Rounded up, **time_log2 = 127.273624**. Preprocessing (adders, starts and restarts) is
included. All processors' work is summed, and elapsed time plays no part.

## 4. Memory and advice

Peak memory is at most the sum of:
- 2·DP_MAX·64 bytes for the distinguished-point table;
- K·384 bytes of walk state (two point buffers, batch scratch, keys);
- 2^16·2·64 bytes of adder images;
- 16 MiB of code and scratch.

That is 2^75.542441 bytes. Nonuniform advice is 0 bytes.

## 5. Limitations

- H1 is the standard but unproved random-mapping heuristic. It is supported by the
  literature and toy-curve experiments only. The model bound is 0.4115, 0.0215 above
  0.39.
- Instruction counts come from one implementation and compiler. The sample is 2^18 steps
  on each of two input sets. The variable-time branchy subtractions were measured, not
  bounded analytically.
- About 2^68.5 stored distinguished points and 2^127.3 work are needed. This is not a
  practical attack and does not weaken P-256.

## References

- Pollard (1978).
- van Oorschot–Wiener (J. Cryptology 1999).
- Teske (Math. Comp. 2001).
- Bos–Kleinjung–Lenstra (ANTS 2010).
- Bernstein–Lange–Schwabe (PKC 2011).
- Bernstein–Lange (ANTS 2012).
- Bos–Costello–Miele (PKC 2014).
- Gueron–Krasnov, *Fast prime field elliptic-curve cryptography with 256-bit primes*
  (J. Cryptographic Engineering 2015), the OpenSSL ecp_nistz256 field method.

## Reproduction

```sh
python3 -I submissions/p256-ecdlp/rho-sigma2/accounting.py
bash submissions/p256-ecdlp/rho-sigma2/measurement/reproduce.sh 1   # Docker, linux/arm64
echo '{"curve":"a3","bits":26,"trials":150,"seed":2601}' | python3 -I submissions/p256-ecdlp/rho-sigma2/experiment.py
```

# p256-ecdlp: parallel Pollard rho on negation classes with a sigma walk

Claim: for uniform x in [0,r) and Q=[x]P on NIST P-256, the algorithm below returns x
with probability at least 0.39. It uses at most **2^127.273624** target-group-operation
equivalents in total (score 15909203/125000), including all preprocessing, failures,
lookup and verification. Peak memory is at most 2^76.542441 bytes. Preprocessing is 2^48.619436 units, and no nonuniform advice is
used. This is a classical cold attack on one target.

The success bound rests on the disclosed heuristic H1, the standard random-mapping model
for r-adding walks on negation classes, and on the rate bound H2. The cost rests on H3,
instruction counts measured with the organizer's toolchain. **There is no reference
operation on this track:** every instruction, including all P-256 field arithmetic and
point additions, is charged at 1/C with the borrowed normalisation C = 1000. This is a
carefully accounted generic attack. It is not a practical break.

**Distinguished points and resolution.** Walks run in trails. Each trail starts from fresh
coins and ends at its first distinguished point, which it records once together with its
coins. A key match is resolved by replaying the two trails (section 1.3, Lemmas R and D).
A fruitless cycle therefore cannot store the same point twice, and a degenerate match has
probability exactly 1/r.

Files:
- `accounting.py` recomputes every number with exact integers (output in
  `evidence/accounting.json`).
- `measurement/` holds the per-step instruction-count harness and the field code
  (`p256_field.h`, `p256_sub.h`), plus the Python replica `walk_ref.py`.
- `evidence/measurement/` holds the counts on two input sets and the checks of the
  restart and abandonment paths.
- `experiment.py` is the declared toy-curve experiment, running the whole procedure;
  `evidence/trails/` holds larger offline runs of the same file.

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
- CAP = 64/theta. Trails have fewer than TRAIL_MAX = CAP + 64 steps.
- DP_MAX = 2·theta·S_max + 2K, where S_max = K·(L + TRAIL_MAX).

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
2. **Adders.** For each adder i < 2^16:
   - draw c_i, d_i uniform in [0,r);
   - set R_i = [c_i]P + [d_i]Q;
   - tabulate R_i and −R_i.
3. **Trails.** Each of the K walk slots runs a sequence of trails.
   - A trail starts from fresh coins a, b uniform in [0,r), with W = [a]P + [b]Q. The
     slot keeps the coins until the trail ends. Coefficients are not tracked while
     walking.
   - A trail ends at its first distinguished point (step 6), or is abandoned: when it
     crosses more than 2^60 anchor boundaries (step 5) without one, that is after between
     CAP − 64 and CAP = 64/theta steps, or when an addition gives O. Either way the slot
     starts a new trail from fresh coins. Every trail has fewer than TRAIL_MAX = CAP + 64
     steps.
   - Each slot takes L steps (the main phase), then finishes its current trail and stops.
     The run fails if more than K trails are abandoned.
4. **Step.** W ← f(W), as one batched incomplete affine addition: Montgomery's
   simultaneous inversion over the 2^16 walks of a batch (3M), then lambda = dy·inv (M),
   lambda^2 (S) and y3 = lambda·(x1 − x3) (M), with one Fermat inversion per batch.
   - If the batch product is zero, the affected walk uses the unbatched complete formula.
5. **Fruitless cycles.**
   - *2-cycles.* Each walk keeps the low 64 bits of the key of the point before its
     current point W. If the new point's key matches, the walk evaluates f(new) unbatched.
     If key(f(new)) = key(W) exactly, then {W, new} is a 2-cycle of the class walk, and
     the walk replaces new by 2·canon(M), where M is the element with the smaller key.
   - *3..64-cycles.* At every batch step divisible by 64, every walk of the batch records
     the low 64 bits of its key as its anchor. All walks of a batch step together, so each
     walk is anchored once every 64 of its own steps. On a later match, the walk iterates
     f unbatched up to 64 times. If its class returns, it escapes from 2·canon(M), with M
     the cycle element of least key.
   - So a walk in a cycle of length ℓ ≤ 64 leaves it within 64 + ℓ steps, unless the
     cycle contains a distinguished point, which ends the trail first.
6. **Distinguished points.** W is distinguished when bits 16..75 of key(W) are zero. The
   test is applied to every point that a step or an escape produces. At a distinguished
   point the trail ends: the record (key(W), a, b), with the trail's coins, goes to a
   hash table keyed by key(W), and the slot starts a new trail. A trail makes at most one
   record, so no trail can store a point twice. The run fails if more than DP_MAX records
   arrive.
7. **Resolution.** The first record whose key is already in the table ends the run. No
   other key match is used.
   - Both trails are replayed, unbatched, from their recorded coins up to their first
     distinguished point, with the rules of steps 4–6 and anchors every 64 of the trail's
     own steps. The replay tracks (a, b) mod r: an image sigma(R_i) = ±R_i has coefficients multiplied
     by (−1)^s, and canon multiplies them by (−1)^s, and an escape doubles the
     coefficients of canon(M). By Lemma R each replay reaches the recorded point.
   - The replays give the common canonical point as [a1]P + [b1]Q = [a2]P + [b2]Q. If
     b1 ≠ b2, output x = (a2 − a1)/(b1 − b2) mod r after checking [x]P = Q. If b1 = b2
     the run fails. By Lemma D this has probability 1/r.

**Merged trails stay merged.** Outside fruitless cycles the next class depends only on the
current class. Escapes fire only after an exact check, and they exit from 2·canon(M),
which depends only on the cycle's set of classes, not on where the walk entered the cycle
or when it noticed. The 64-bit filters only decide when an exact check runs. So two trails
that reach the same class follow the same classes afterwards and reach the same first
distinguished class.

**Lemma R (replays reach the recorded point).** The points of a trail, with repeated
traversals of cycles removed, depend only on its start point, not on when its cycles are
detected.

*Proof.* While no cycle is entered, the next point depends only on the current point. If
the trail enters a cycle that contains a distinguished point, it reaches that point within
one traversal, and no detection fires earlier. The 2-cycle rule fires only on a return to
the class of two steps back, and an anchor match only on a return to the anchor's class.
In both cases every point of the cycle has already been produced and tested. Otherwise
the trail repeats the cycle until a detection fires, and the escape 2·canon(M) depends
only on the cycle's classes. The distinguished-point test reads only the key. So any
schedule that anchors each walk at least once every 64 of its steps leads to the same
first distinguished point, and every coefficient pair the replay carries represents its
current point. ∎

Replaying a trail of n < TRAIL_MAX steps takes fewer than 2^7·n unbatched operations. Per
step that is:
- one evaluation;
- at most two more for a 2-cycle check;
- at most 66 for an anchor check, which fires at most once per 65 steps unless it escapes;
- at most 128 extra steps in each cycle before its escape, where each cycle took at least
  2 of the trail's steps.

**Lemma D (degenerate matches).** The first key match has b1 = b2 with probability exactly
1/r.

*Proof.* A trail records once, so the two records come from different trails t1 ≠ t2,
possibly of the same slot. Condition on all start points and adders. That fixes every
trail's points, and so which records match first.
- Given its start point W = [a + b·x]P, a trail's coin b is uniform on [0, r), because
  each b determines exactly one a.
- Coins of different trails are independent.
- Along a replay, b changes as b ← u·b + d, where u is a unit (±1) or twice
  one, and u, d depend only on the path.

So the canonical b at the distinguished point is U_i·b_ti + D_i with U_i invertible mod r,
and P[b1 = b2] = 1/r. ∎

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

**Detection and resolution.** Consider a class collision among the first T_steps steps
between two different trails.
- From the shared class both trails follow the same classes to the same first
  distinguished class, since merged trails stay merged.
- Every trail finishes, so whichever trail records second triggers resolution.
- Lemma R shows that the replays reach that point, and Lemma D shows that the equation is
  degenerate with probability exactly 1/r.
- Fresh start points only add points, and the bound does not count them.

An earlier version continued walks through distinguished points. There, a fruitless cycle
that contains a distinguished point could store that point repeatedly with equal
coefficients. Restarting at every distinguished point removes that case: each record comes
from a different trail with its own coins, and such a cycle costs one ordinary record,
inside DP_MAX.

The failure terms are:
- **Collision inside one trail.** There are at most DP_MAX + 2K ≈ 2^68.54 trails, each with
  fewer than TRAIL_MAX ≈ 2^66 steps. The chance that any trail meets one of its own
  earlier classes, outside the short cycles handled by step 5, is at most
  (DP_MAX + 2K)·TRAIL_MAX^2/N ≈ 2^−54.5.
- **Abandoned trail.** A trail is abandoned only after more than CAP − 64 steps without a
  distinguished point outside short cycles. That has probability at most
  e^−63.99 ≈ 2^−92.3 per trail, so 2^−23.8 over all trails.
- **Degenerate first match:** exactly 1/r < 2^−252 (Lemma D).
- **DP_MAX or abandonment overflow:** under 2^−1000 by Chernoff.
- **Exceptional additions:** under 2^−110.

Together these are below delta = 2^−20. Therefore:

**Pr[success] ≥ 1 − exp(−0.5303·(1 − 0.00015268)) − 2^−20 = 0.411523 ≥ 0.39.**

The budget is 1.0299·√N steps, against 0.9943·√N for success exactly 0.39. The 0.0215
of slack in the model bound covers the sampling error of the toy evidence below.

**Evidence for H1, H2 and the procedure.** `experiment.py` runs the procedure of
section 1.3 end to end on random prime-order toy curves
y^2 = x^3 − 3x + b (`"curve": "a3"`). It uses:
- the same class key, canonical rule, adder images and index (keys use Montgomery
  representatives, any fixed bijection being equivalent);
- both cycle rules and their escapes;
- trails from fresh coins that end at their first distinguished point, and the trail cap;
- resolution by replaying both trails from their coins.

It never tracks coefficients while walking and has no oracle collision table.
- The live walk anchors on a global round counter, as the batched kernel does. The replay
  anchors on the trail's own step count, so every replay also checks Lemma R.
- Each trial runs L steps per slot plus the finishing phase at the stated tau. It counts
  as a success only if the first key match resolves to a verified x.
- Each trail start is charged as one of the slot's L steps, which is conservative.

<!-- toy-evidence -->
Offline runs (`evidence/trails/`, pooled in `summary.json`): 5,100 trials on 5 curves with
N = 2^28.6, 2^30.1, 2^32.9, 2^32.4, 2^34.7 classes,
adder tables of 2^11 and 2^16, 16 slots and theta = 2^−5.
- Success of the whole procedure at tau = 0.5303: 2,172 of 5,100, 0.4259, Wilson 95% interval
  [0.412, 0.440]. The model 1 − exp(−tau) is 0.412.
- Every success was verified ([x]P = Q). Degenerate first matches: 0. Replays that did not
  reach the recorded point: 0. Every match was replayed; about 8.0 million trails were run.
- Abandoned trails: 1. Longer fruitless cycles escaped: 19. The finishing phase added at most
  1.4% to the K·L steps.
- 2-cycle entry rate per step:
  - 2.31e-04, 2.28e-04, 2.29e-04, 2.27e-04 with A = 2^11 (1/(2A) = 2.44e-04);
  - 6.91e-06 with A = 2^16 (1/(2A) = 7.63e-06).
<!-- /toy-evidence -->

The declared initial input (`curve` "a3", `bits` 28, 100 trials, seed 2801, the package
defaults otherwise) runs in a few seconds. It reports 41/100 = 0.41, Wilson 95% interval
[0.319, 0.508], against the model 0.412, with 0 degenerate matches and 0 replay
mismatches. With 100 trials the standard error is about 0.05, so it is a reproducibility
check. The success claim rests on the pooled offline trials above, and any input within
the schema can spot-check them.

Limits: the toy curves have N up to about 2^35, against 2^255 in production. There theta is 2^−5,
so start points and the finishing phase are a few percent of the steps; in production
they are below 2^−30 of them. H1 is an extrapolation in N, supported by the size-free
random-mapping argument and the cited literature. It is not proved for this walk.

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
2^16 walks with 2^16 adders, including restarts at distinguished points and the
trail cap. In this harness a new trail takes its start point from spare points in the
input file. In production it comes from fresh coins, and that cost is charged below as
trail restarts.
- It was compiled with `-O3 -fno-tree-vectorize -mgeneral-regs-only` in the organizer's
  Docker image (GCC 14.2.0, Valgrind 3.24.0), and counted with Callgrind `Ir`.
- Counts are differenced between 8 and 4 batch-steps of 2^16 walks, 2^18 steps.
- Every instruction inside the measured loop is charged.

| Input set | total Ir/step | 2-cycles in the larger run |
| ---: | ---: | ---: |
| seed 1 | 819.933 | 2 |
| seed 2 | 820.567 | 3 |

This covers:
- 5 multiplications and 1 squaring;
- 5 subtractions, one of them fused;
- key, parity and adder-image selection;
- the 2-cycle, anchor and distinguished-point filters, the trail-cap counter, and any
  confirmations, escapes and restarts;
- the batch's prefix stores and loads;
- one Fermat inversion per 2^16 walks (about 38,000 Ir, or 0.6 per step);
- loop control.

**830 instructions** are charged per step, so each step costs **0.83 units**.
`measurement/walk_ref.py check` replays all walks in Python integers and matches the C
final points and event counts exactly on both input sets.

Distinguished points never occur at theta = 2^−60 in such a sample. So the restart and
abandonment paths were checked on small builds of the same source and replica
(`evidence/measurement/trail-paths.json`), with 2^6 adders and 2^12 walks:
- With DP_BITS = 4, the run had 49,929 restarts at distinguished points and 6,955 2-cycles.
- With DP_BITS = 12 and a one-period cap, it had 4,093 abandonments.

In both, the replica matched every final point and every count.

Rare paths fit in the 9.4-instruction margin:
- **2-cycles.** At most 1/A per step × 3 unbatched operations × 2^16 instructions,
  that is at most 3 instructions per step.
- **Restarts and abandonment.** A restart in the kernel copies a start point and selects
  its adder, about 50 instructions, at rate theta per step: below 2^−50 instructions per
  step. Its start point is charged in the ledger.
- **Short cycles and exact-k fallbacks.** Below 2^−15 instructions per step.

A separately measured unbatched step (affine add with its own inversion, plus adder
selection) costs 40,871 on both input sets Ir. Every unbatched operation below is charged 2^16.

**Total.** S_max = K·(L + TRAIL_MAX) ≤ 2^127.542441 steps.

| Component | Bound (units) | log2 |
| --- | --- | ---: |
| walk steps | 0.83·S_max | 127.273624 |
| adder table | 2^16 adders × 2 scalar mults × 768 unbatched ops × 2^16 / C | 32.62 |
| initial walk starts | K × 2 scalar mults × 768 × 2^16 / C | 48.62 |
| trail restarts | (DP_MAX + K) × 2 scalar mults × 768 × 2^16 / C | 85.16 |
| distinguished points | DP_MAX = 2^68.54 entries × 2^20 instructions / C | 78.58 |
| collision replay | 2 trails × TRAIL_MAX steps × 2^7 unbatched ops × (2^16 + 2^10) / C | 80.06 |
| solve, verify, output | 2 scalar mults + 2^24 instructions | 16.84 |
| input validation | 2^24 instructions | 14.03 |

Rounded up, **time_log2 = 127.273624**. Preprocessing (adders and initial walk starts) is
included. All processors' work is summed, and elapsed time plays no part.

## 4. Memory and advice

Peak memory is at most the sum of:
- 2·DP_MAX·128 bytes for the distinguished-point table (key, the trail's coins, hashing
  overhead);
- K·512 bytes of walk state (points, keys, coins, trail counter, batch buffers);
- 2^16·2·64 bytes of adder images;
- 16 MiB.

That is 2^76.542441 bytes. Nonuniform advice is 0 bytes.

## 5. Limitations

- H1 is the standard but unproved random-mapping heuristic. It is supported by the
  literature and toy-curve experiments only. The model bound is 0.4115, 0.0215 above
  0.39.
- Instruction counts come from one implementation and compiler. The sample is 2^18 steps
  on each of two input sets. The variable-time branchy subtractions were measured, not
  bounded analytically.
- About 2^68.5 stored distinguished points and 2^127.3 work are needed. This is not a
  practical attack and does not weaken P-256.
- The toy experiment runs the whole procedure, including resolution, but only at toy
  scale: theta = 2^−5, 16 slots, and N far below production. Its trail starts and its
  finishing phase are a few percent of the steps there, and negligible in production.
- Distinguished points never occur in a production-parameter measurement sample. The
  kernel's restart and abandonment paths were checked on small builds of the same source.

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
echo '{"curve":"a3","bits":28,"trials":100,"seed":2801,"adder_bits":11,"dp_bits":5,"walks":16,"tau_e4":5303}' | python3 -I submissions/p256-ecdlp/rho-sigma2/experiment.py
```

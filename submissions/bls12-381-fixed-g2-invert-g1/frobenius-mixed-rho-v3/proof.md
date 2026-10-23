# frobenius-mixed-rho-v3: Frobenius-class mixed add/square rho on GT (v3)

Claim: for the uniform challenge of `rules.md` (Q uniform nonidentity in G2, z uniform in GT, fresh algorithm coins), the algorithm
below outputs P in G1 with e(P, Q) = z with
probability at least 0.41 (model 0.410001, lower bound 0.4100000006), using at most **2^125.719085**
target-group-operation equivalents (score 25143817/200000), including all preprocessing, the adder table, walk starts and
restarts, distinguished-point storage, cycle waste, fingerprint verification, collision replay, the final scalar
multiplication and verification. Peak memory is at most 2^73.166751 bytes; preprocessing at most 2^75.418470 units; nonuniform
advice 0 bytes. The bound rests on heuristics H1-H4. It is a classical generic attack and not a practical break.

**Revision v3.** The reviewer rejected v2 (`frobenius-mixed-rho-v2`, claim 125.747292) on four points. v3 answers each:
(1) the anchor was refreshed on a global counter, so only some walks ever received anchors; the kernel now refreshes each
walk's anchor on its own step count, and the longer-cycle branch was exercised once in the production runs (2-cycles: 6050 events); (2) the collision rate used
the wrong partition; V is now 0.799413681 and the budget is re-derived from it; (3) c2 was extracted with fractions from
the full run and the square fraction was bracketed wrongly; c2 and kappa are now measured at S = 458 with interval counts,
and the operation mix comes from uninstrumented runs; (4) the distinguished-point record stores a 16-byte fingerprint with
a verification step for false matches. A fifth defect, found while answering (4): the v2 kernel stored raw exponents in distinguished points, but the solver needs the canonical exponents lambda^kn a and lambda^kn b (sigma^kn(x) is the canonical element). v3 stores the canonical exponents, and the replica checks every distinguished-point record bit for bit. The 3% overhead of v2 is replaced by a measured bound on cycle waste. Net effect on
the score: 125.719085 against 125.747292 for v2 (-0.028207 bits).

## 1. Algorithm

Reduction. g = e(P0, Q) generates GT, where P0 is the fixed G1 generator; z = g^y gives P = [y]P0 with e(P, Q) = z. Input a nonidentity Q in G2 and a target z in GT. Return the identity if z = 1 (P = O is valid).
Otherwise the walk below finds y with g^y = z, and the output is P = [y]P0, verified by the pairing equation e(P, Q) = z
(one pairing, charged). Exceptional inputs (Q = O) are rejected by the rules.

Walk (`walk-spec.md`, the specification used by `cscore/gt.py` and the blst kernel `measurement/step.c`):
- GT has order r. Frobenius sigma(x) = x^p has order 12 with lambda = p mod r, lambda^6 = -1; classes are sigma-orbits,
  N = (r - 1)/12, about 2^251.272127.
- Key K(x): the minimum over the 12 values V_n(sigma^n(x)), computed with four Fp multiplications; K is class-invariant.
- State (x, a, b) with x = g^(a + b y). Square slots S = 458 of r_a = 1024 (s = 0.447265625000). Step: h = K mod r_a.
  If h < S, square x and double (a, b). Otherwise multiply x by the twisted table entry T[h][k] = sigma^(-k)(R_h), where
  R_h = g^(c_h) z^(d_h), and update (a, b) by lambda^(-k) (c_h, d_h). The Frobenius map is never applied per step.
- Cycles (per walk): a state whose key equals the key two steps back and whose successor returns to the old key is a
  2-cycle. A state whose key equals the walk's anchor key (recorded every 64 of its own steps) inside a window of 64 steps
  is a cycle of length at most 64. Either escapes by squaring the least-key state of the cycle, which commutes with sigma.
  Escapes are group operations and are charged; no walk is stopped or restarted by the detector.
- Distinguished point: (K >> 10) & (2^60 - 1) = 0. The record is (fingerprint, a, b): the low 128 bits of the Montgomery
  integer of K, and the canonical exponents, 80 bytes. Two records with equal fingerprints are checked: if they come from
  the same class with B1 != B2, y = (A2 - A1)/(B1 - B2) mod r; otherwise the candidate is rejected. Every candidate is
  verified by g^y = z before output.
- Stopping policy. Detection never stops a walk; each detected cycle costs one escape and at most 127 main steps of
  waste, charged through the budget (Section 2). Without that charge the success bound would be
  0.406360 instead of 0.410001.

Consistency checks: the walk invariant x = g^(a + b y) holds on every state (replica at production parameters, 512 and 4096
rounds, 64 walks; toy tests); the twisted table equals sigma^(-k)(g^c z^d) (test); the solver recovers y end to end (tests and
100 toy trials, all verified).

## 2. Success probability and budget

Under H1 the walk needs T = 2^125.836362 fresh steps, with tau = T^2 V / (2N) = 0.527634 (grid point 0.5277).
The model success is 0.410001 from `exact.birthday_success` with eps = 1 - V and delta = 2^-20 (lower bound
0.4100000006 >= 0.41). The main steps executed are T / (1 - w) with w = 0.005848 (H4), which is
2^125.844823 in total including trails. The uniform challenge gives a uniform discrete log because g generates GT.

**Heuristics (disclosed).**
- **H1 (random-mapping collision model).** The walk on the N = (r-1)/12 Frobenius classes behaves as a random map
  whose collision rate is V = 1 - s^2 - (r_a - S)/r_a^2, the sum of squared branch probabilities for the actual partition:
  the square branch has probability s = S/r_a, and each of the r_a - S additive adders has probability 1/r_a (step rule
  in `walk-spec.md`). With S = 458 this gives V = 0.799413681. The v2 value (1 - s^2 - (1 - s)^2/r_a = 0.799655)
  assumed all r_a adders carry the additive probability and is withdrawn. tau = T^2 V / (2N) is treated as Exp(1).
  Support: RESEARCH.md 3.4 and the toy sweeps below. Limits: toy groups have N up to 2^17 and 32 adders; production has
  N about 2^251.3 and 1024 adders.
- **H2 (cycles, per walk).** Every walk checks each new key against the key two steps back (2-cycles) and against an
  anchor key recorded every 64 of its own steps (cycles of length 3 to 64). A detection escapes by squaring the least-key
  state of the cycle. The detection window is 2 steps for 2-cycles and at most 127 main steps for longer cycles. The
  production runs detected 6051 events in 134217728 main steps (6050 2-cycles and
  1 longer cycle). The 2-cycle branch is tested for exact equality against the replica;
  the longer-cycle branch is mirrored in the replica but was not exercised by an exact-equality test.
- **H3 (measured per-step costs).** c2 and kappa are measured at the production square count S = 458 by Callgrind on
  two input seeds, from interval counts (Computation section). The per-step operation counts are measured on uninstrumented runs.
- **H4 (wasted steps).** A detected cycle wastes at most 127 main steps before its detection. The wasted fraction
  is bounded by 127 times the one-sided 95% upper bound on the event rate, w = 0.005848. The
  total main steps are T / (1 - w), plus one trail of 2^60 steps per walk, which is charged in full.

## 3. Computation

Reference unit: one dense Fp12 multiplication (`blst_fp12_mul`). Measured 22,998 instructions (1.1499 units; C = 20,000 is
the organizer's rounding), charged 1 as the rules require. Other work is measured in instructions and divided by C.

Measurement (`measurement/`): the organizer's reference image (linux/arm64, GCC 14.2.0, Valgrind 3.24.0,
`-O3 -fno-tree-vectorize -mgeneral-regs-only`), the pinned blst commit f7a88067. Callgrind Ir is differenced between 256
and 512 rounds of 64 walks (16,384 main steps) at S = 458, the production square count. The operation counts come from the
same interval. Each group operation (adders, squarings, escapes and the extra steps of cycle checks) is counted.

| input set (S = 458) | callgrind interval (steps) | adds per step (interval) | squares per step (interval) | other per main step (instructions) | c2 = other / C | kappa = squaring / C |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| seed 1 | 16384 | 0.547119 | 0.453003 | 2886.024 | 0.144301 | 0.49020 |
| seed 2 | 16384 | 0.552429 | 0.447571 | 2885.735 | 0.144287 | 0.49020 |

c2 is charged at 0.1445 (the larger seed, 0.144301, with a 0.1% margin, rounded up). Distinguished-point work is excluded from c2 (measured with the DP density set to zero) and charged separately at 2^20 instructions per point. kappa is charged at 0.4902.
Mixed-input composition does not enter c2: the per-step operation counts multiply their own unit costs (adds 1, squares
kappa), so the walk step costs 0.552748 + 0.447342 kappa + c2 = 0.916535 units, with the
adds and squares per step taken from the uninstrumented runs (below).

Uninstrumented runs of the same kernel (`measurement/native_run.sh`, `evidence/native-counts.json`): two input seeds at
S = 458, 2^20 rounds each, 2^26 main steps per seed. Group operations per main step:
adds 0.552747704, squarings 0.447342485 (both include escapes and cycle-check steps).
Cycle detections: 6051 in 134217728 main steps, rate 4.5083e-05; one-sided 95% upper bound 4.6044e-05
(normal approximation to the Poisson). The RESEARCH 3.5 2-cycle formula (1-s)^2/(12 r_a) = 2.4863e-05 is for a
uniform partition; the measured 2-cycle rate is 6050 / 134217728
= 4.5076e-05, a ratio of 1.81.

Components (log2 units; the score is the rounded-up total 125.719085):

| component | log2 units | what it covers |
| --- | ---: | --- |
| walk steps | 125.719085 | total main steps 2^125.844823 (T / (1 - waste), plus trails), each costing 0.916535 units |
| adder table | 19.582053 | 1024 adders, two exponentiations and 12 Frobenius maps each (preprocessing) |
| walk starts and restarts | 75.418470 | (walks + trails) times two exponentiations (preprocessing) |
| distinguished points | 72.557110 | storage and lookup, 2^20 instructions each, 2^66.844823 points |
| collision replay | 67.194717 | twice the cap 64 x 2^60, each step costing 1 + c2 |
| fingerprint verification | 13.263293 | false fingerprint matches, about dps^2 / 2^129, each verified by one exponentiation |
| solve, verify, output | 14.712288 | one scalar multiplication, one pairing, input checks (2^29 instructions) |

`accounting.py` recomputes all of this with exact fractions and rounds up. Its output is `evidence/accounting.json`.

## 4. Memory

Peak storage: distinguished-point records (80 bytes each, 2^66.8448 of them) dominate, at
2^73.166751 bytes including the adder table (12 x 1024 elements, about 7.8 MiB), the 2^20 walk states, and per walk two
48-byte keys and two valid flags (the cycle detector). The fingerprint is the low 128 bits of the Montgomery integer of K.
Two records with equal fingerprints but different classes occur about dps^2 / 2^129 times in expectation (about 2^4.7
here); each is rejected by one exponentiation, charged above. Nonuniform advice is 0 bytes.

## 5. Evidence: toy-scale sweeps

Group: the order-ell subgroup of F_(q^12)^* with ell dividing q^4 - q^2 + 1; `cscore/gt.py`, bundled into `experiment.py`.
The toy taus use the partition-correct V (32 adders, S = 15). The toy walk detects full cycles with its own rules, and
distinguished points are dense, so these runs test H1's shape and the solver, not H2 or H4.

| run | ell | N | s (effective) | V (toy) | trials | reached a useful collision | mean tau (model 1) | success at tau = 0.5277 [Wilson 95%] |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| initial input (bits 14, seed 1, 100 trials) | 28393 | 2366 | 0.4688 | 0.763672 | 100 | 100 | 0.8625 | 0.4800 [0.3846, 0.5768] |
| add-only, bits 14 | 28393 | 2366 | 0.0000 | 0.968750 | 3000 | 3000 | 0.9634 | 0.4170 [0.3995, 0.4347] |
| mixed, bits 14 | 28393 | 2366 | 0.4688 | 0.763672 | 3000 | 2997 | 0.9592 | 0.4277 [0.4101, 0.4455] |
| add-only, bits 18 | 922561 | 76880 | 0.0000 | 0.968750 | 800 | 800 | 0.9514 | 0.4263 [0.3924, 0.4608] |
| mixed, bits 18 | 922561 | 76880 | 0.4688 | 0.763672 | 800 | 800 | 0.9987 | 0.4088 [0.3752, 0.4432] |

Pooled mixed runs at tau = 0.5277 (N = 2366 and N = 76880): 1610 of 3800 succeed, empirical 0.4237, Wilson 95%
interval [0.4081, 0.4395]. The lower limit is above the 0.39 threshold the rules set. The toy success is at or
above the model value 1 - exp(-tau) = 0.4100 at this point, so the model is not optimistic on these runs (mean tau 0.96-1.00).

Cycle rates in the toy runs (toy detector, which also escapes full cycles):

| toy walk (r_a = 32) | measured 2-cycles per step | RESEARCH 3.5 formula | ratio | cycle-step fraction (toy detector) |
| --- | ---: | ---: | ---: | ---: |
| mixed, N = 2366 | 0.004061 | 0.000735 | 5.5x | 4.92% |
| mixed, N = 76880 | 0.001943 | 0.000735 | 2.6x | 0.41% |
| add-only, N = 2366 | 0.003208 | 0.002604 | 1.2x | 1.31% |
| add-only, N = 76880 | 0.002632 | 0.002604 | 1.0x | 0.57% |

## 6. Limitations

- H1 is extrapolated from toy subgroups with N up to 2^17 and 32 adders to N about 2^251 and 1024 adders.
- The longer-cycle branch (3 to 64) fired once in 2^27 production steps (seed 2). Its code is mirrored line by line in the
  replica, but no exact-equality test exercised it. The 2-cycle branch is exercised by the replica checks (512 rounds: 1
  event; 4096 rounds: 11 events, all matched).
- Repeated distinguished points (the duplicate-DP question). A detected cycle can repeat states before detection. A
  repeat carries the same canonical exponents as the first record, so its match is degenerate (B1 = B2) and is discarded.
  No prefix replay is needed, because each record stores its exponents, so degenerate matches cost one lookup and never
  end a trial. Repeats are counted among the distinguished points, so their storage and lookup are charged at 2^20 each.
  The expected number of repeated records is about the total cycle events x 127 x 2^-60 (the window bound), which
  is about 2^58.4; this is an upper bound and is already inside the distinguished-point term. Detection is still
  needed: a trapped walk would otherwise run its whole segment, about 2^60 steps, without a distinguished point escape.
- Measured costs come from one blst build, one compiler, and two input seeds at S = 458. The two seeds agree to 0.01%.
- The replica checks full dumps and distinguished-point records at 512 and 4096 rounds (bit-exact), plus samples of production states, not all states.
- Memory is about 2^73.166751 bytes, which is not a realisable storage budget; it is reported, not scored.
- The reference operation is priced by the organizer rule at 1 unit; the measured instruction count is 1.1499 units, so the
  rule's charge is slightly favourable to the attack.
- The optimum square count at the charged c2 is 456, not 458, by a margin of about 2e-6 relative in cost per step; the walk is
  measured at 458 and charged there.

## References

- Pollard (1978); van Oorschot and Wiener (1999); Teske (1998), section 4.2; Bernstein and Lange (2012), section 3.
- `work-gt/RESEARCH.md` sections 3-4 and `walk-spec.md` (this package).

## Reproduction

```sh
python3 -I accounting.py
echo '{"bits":14,"trials":100,"seed":1,"s":0.46,"adder_bits":5,"dp_bits":2,"walks":2}' | python3 -I experiment.py
bash measurement/reproduce.sh DIR        # Docker linux/arm64, organizer image, py_ecc (gen_inputs.py)
```
`reproduce.sh` generates inputs at S = 458, measures the kernel (Callgrind, two sizes, mul and sqr loops), analyses the
interval counts, and runs the replica on each input set. `native_run.sh` runs the uninstrumented kernel on the organizer
image for the event counts.

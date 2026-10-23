# frobenius-mixed-rho-v3: Frobenius-class mixed add/square rho on GT (revision v4)

Claim: for the uniform challenge of `rules.md` (Q uniform nonidentity in G2, z uniform in GT, fresh algorithm coins), the procedure below
outputs P in G1 with e(P, Q) = z with probability at least 0.41 (model 0.410000000001; lower bound
0.410000000001 after the failure terms of Section 2), while its total work, counting every operation the procedure
performs, is at most W_MAX = 2^125.719870 target-group-operation equivalents (score 12571987/100000). Peak memory is at most
2^73.166749 bytes; preprocessing is included in W_MAX (2^76.418468 units); nonuniform advice 0 bytes. The bound rests on heuristics
H1, H2, H3, H5 and H6 (Section 2) and on the derived failure terms. It is a classical generic attack, not a practical break.

**Revision v4.** The reviewer rejected v3 (round 3) on three points: (1) undetected long cycles and trapping had no
trail cap, restart policy or derived failure terms; DELTA = 2^-20 was a placeholder; (2) the fingerprint used key bits
10..69 that are zero for every distinguished point, so it had 68 varying bits, and the false-match count was low by 2^60;
(3) no hard total-work cutoff. v4 answers each:
(1) each walk has a trail age; at 2^66 steps without a distinguished point the walk restarts (charged). Trapping is bounded
by formal and non-formal cycle rates (Section 2), giving delta_trap = 1.29e-18; the trail-cap loss is 2.1e-26;
DELTA is the sum of the kept terms, not 2^-20;
(2) the fingerprint is key bits 70..197 of the Montgomery integer, 128 varying bits; the false-match count is 2^4.69,
and its verification costs 2^13.26 units (round 3 had 2^73.26);
(3) the procedure stops when its work counter reaches W_MAX - reserve; W_MAX is the claimed computation, so every "at most" holds
by construction. The success loss from stopping is derived in Section 2.
Changes from v3 that affect the score: the charged c2 rises from 0.1445 to 0.1450 for the trail-age counter (about 11
instructions per step), and the collision replay term of v3 is removed because the solver does not replay trails. The score
moves from 125.719085 to 125.719870 (+0.000785 bits).

## 1. Procedure

Reduction. g = e(P0, Q) generates GT, where P0 is the fixed G1 generator; z = g^y gives P = [y]P0 with e(P, Q) = z. Input a nonidentity Q in G2 and a target z in GT. Return the identity if z = 1 (P = O is valid).
Otherwise run the walk below until a verified y is found or the work counter reaches W_MAX, then output P = [y]P0, verified by
the pairing equation e(P, Q) = z (one pairing, charged). Exceptional inputs (Q = O) are rejected by the rules.

Walk (`walk-spec.md`, the specification used by `cscore/gt.py` and the blst kernel `measurement/step.c`):
- GT has order r. Frobenius sigma(x) = x^p has order 12 with lambda = p mod r, lambda^6 = -1; classes are sigma-orbits,
  N = (r - 1)/12, about 2^251.272127.
- Key K(x): the minimum over the 12 values V_n(sigma^n(x)), computed with four Fp multiplications; K is class-invariant.
- State (x, a, b) with x = g^(a + b y). Square slots S = 458 of r_a = 1024 (s = 0.447265625000). Step: h = K mod r_a.
  If h < S, square x and double (a, b). Otherwise multiply x by T[h][k] = sigma^(-k)(R_h), where R_h = g^(c_h) z^(d_h), and
  update (a, b) by lambda^(-k) (c_h, d_h). The Frobenius map is never applied per step.
- Cycles (per walk): a new key equal to the key two steps back, with the following step returning to the old key, is a 2-cycle.
  A new key equal to the walk's anchor key (recorded every 64 of its own steps) inside a window of 64 steps is a cycle of length
  at most 64. Either escapes by squaring the least-key state of the cycle. Escapes are group operations and are charged.
- Trail cap and restart (per walk): the walk's trail age counts main steps since its last distinguished point. At age 2^66
  (64 x 2^60) the walk restarts from a fresh random start (2 exponentiations, charged) and its age resets.
- Distinguished point: (K >> 10) & (2^60 - 1) = 0. The record is (fingerprint, A, B), 80 bytes: the fingerprint is key bits
  70..197 of the Montgomery integer of K; A = lambda^kn a and B = lambda^kn b mod r, where kn is the argmin index.
- Solving. Records with equal fingerprints are compared pairwise. A pair from one class with B1 != B2 gives
  y = (A2 - A1)/(B1 - B2) mod r; A repeated distinguished point from a detected cycle carries the same canonical exponents as its first record, so the pair is degenerate (B1 = B2) and discarded at the cost of one comparison; no trail is replayed. Every candidate is verified by g^y = z before output.
- Work counter and cutoff. Each operation adds its charge to W (Section 3): main steps, escapes and cycle checks, distinguished
  points (2^20 instructions), starts and restarts, verified fingerprint comparisons, the adder table, and the final solve. The
  walks run while W < W_MAX - reserve, with reserve the solve and output charge. Then the procedure reports failure.

Consistency checks: the walk invariant x = g^(a + b y) holds on every state (replica at production parameters, 512 and 4096 rounds,
64 walks); the canonical exponents of every distinguished point satisfy sigma^kn(x) = g^(lambda^kn (a + b y)) (replica, bit for bit
on the stored records); the twisted table equals sigma^(-k)(g^c z^d) (test); the solver recovers y end to end (tests, and
100 toy trials, all verified).

## 2. Success probability and failure terms

Under H1 the walk needs T = 2^125.836359 fresh steps for tau = T^2 V / (2N), with V = 0.799413681. The collisions the
distinguished-point test can see are a fraction f = 1 - 2.1e-26 of all collisions (the cap-loss term). Each walk's
main steps are T / (1 - w) with w = 0.005848 charged for detected cycles (127 steps per detection, rate upper bound
4.6044e-05), plus a detection lag of 2^80 steps (every walk reaches its next distinguished point before
a collision is seen). The success probability is

  P = 1 - exp(-f tau_eff) - delta, with tau_eff = T^2 V / (2N) = 0.527632742, delta = 1.290e-18.

This gives P = 0.410000000001 before delta and 0.410000000001 after it. If the wasted steps were not budgeted, P would fall to
0.406359; that is the success loss the cutoff prevents.

Failure terms (each is a probability over the whole run; derivations in `accounting.py`, `derive_failure_terms`):
| term | value | derivation |
| --- | ---: | --- |
| trap, formal cancellations of length >= 66 | 2^-324.2 per step | (ell-1)!! ((1-s)/(12 r_a))^(ell/2) summed over even ell >= 66 |
| trap, non-formal closure | 2^-185.3 per step | 12 CAP / r candidates per step (heuristic H5) |
| delta_trap | 1.290e-18 (2^-59.43) | total steps x per-step trap rate (union bound) |
| trail-cap loss | 2.085e-26 | 2 (1 + c) e^(-c), c = 64: a collision whose downstream DP lies beyond the cap |
| delta_waste | below 1e-300 (underflow) | concentration of the detected-cycle count (1% excess, Chernoff) |
| delta (kept) | 1.290e-18 | delta_trap + delta_waste |

Formal cancellation (the 2-cycle is the length-2 case). The group is abelian and lambda^6 = -1, so a multiset of adder steps
cancels exactly when each adder appears in pairs of rotations k, k + 6. A pair costs one additive first step (probability 1 - s)
and one specific second step (probability 1 / (12 r_a)). The rate of a formal cancellation of length ell is
(ell - 1)!! ((1 - s) / (12 r_a))^(ell/2). For ell = 2 this is 4.498164e-05 per step; the measured 2-cycle rate is
4.507601e-05, a ratio of 1.0021. Cycles of length 66 and above are formal only in this way (or non-formal with
probability 1/r per candidate), which gives the trap terms above.

Heuristics (disclosed):
- **H1 (random-mapping collision model).** The walk on the N = (r-1)/12 Frobenius classes behaves as a random map
  whose collision rate is V = 1 - s^2 - (r_a - S)/r_a^2 = 0.799413681: the square branch has probability s = S/r_a and each of the
  566 reachable additive adders has probability 1/r_a. tau = T^2 V / (2N) is treated as Exp(1), and a collision is detected
  when its downstream distinguished point is reached inside the trail cap (factor f, derived below).
  Support: RESEARCH.md 3.4 and the toy sweeps. Limits: toy groups have N up to 2^17 and 32 adders.
- **H2 (short cycles, per walk).** 2-cycles are checked every step and cycles of length 3 to 64 against a per-walk anchor
  refreshed every 64 of the walk's own steps. Detected cycles escape by squaring the least-key state; each wastes at most
  127 main steps. The 2-cycle branch is checked for exact equality against the replica; the 3-to-64 branch is mirrored in
  the replica and fired once in production (no exact-equality test).
- **H5 (undetected cycles, derived).** A cycle of length above 64 can only be undetected through (a) a formal cancellation of
  adder steps (pairs (h, k), (h, k + 6)), whose rate is the explicit formula in the Computation section, or (b) a non-formal
  closure, which requires the affine exponent relation (lambda^j - 2^q) e = (adder terms) to hold for the secret y; this
  holds with probability 1/r per candidate for uniform y. The formal rate is checked against the measured 2-cycle rate
  (ratio 1.0021).
- **H6 (no degenerate square multiplier).** For every square count q <= 2^66, 2^q is not a 12th root of unity mod r, so a
  cycle with a square has a non-zero multiplier. Checked for 2^(12k) != 1 mod r for k <= 200000; not proven beyond.
- **H3 (measured per-step costs).** c2 and kappa are measured at S = 458 on two seeds with zero DP density (the distinguished-point
  work is charged separately); the operation counts are measured on uninstrumented runs.

## 3. Computation

Reference unit: one dense Fp12 multiplication (`blst_fp12_mul`), measured 22,998 instructions (1.1499 units; C = 20,000 is the
organizer's rounding), charged 1 as the rules require. Other work is measured in instructions and divided by C.

Measurement (`measurement/`): the organizer's reference image (linux/arm64, GCC 14.2.0, Valgrind 3.24.0,
`-O3 -fno-tree-vectorize -mgeneral-regs-only`), the pinned blst commit f7a88067. Callgrind Ir is differenced between 256 and 512
rounds of 64 walks (16,384 main steps) at S = 458, with the distinguished-point density set to zero so that the DP work is excluded
from c2 (it is charged separately at 2^20 instructions per point).

| input set (S = 458, zero DP density) | adds per main step (interval) | squares per main step (interval) | other per main step (instructions) | c2 = other / C | kappa = squaring / C |
| --- | ---: | ---: | ---: | ---: | ---: |
| seed 1 | 0.547119 | 0.453003 | 2897.008 | 0.144850 | 0.49020 |
| seed 2 | 0.552429 | 0.447571 | 2896.720 | 0.144836 | 0.49020 |

c2 is charged at 0.1450 (the larger seed, 0.14485, with a 0.1% margin, rounded up). kappa is charged at 0.4902.
The walk step costs 0.552748 + 0.447342 kappa + c2 = 0.917035 units, with the operation mix taken
from uninstrumented runs of the same kernel (2^26 main steps per seed, both including escapes and cycle-check steps).
Cycle detections in those runs: 6051 in 134217728 main steps (6050 2-cycles, 1 longer).

Components of W_MAX (log2 units; the score is the rounded-up total 125.719870):

| component | log2 units | what it covers |
| --- | ---: | --- |
| walk steps | 125.719869 | 125.844820 main steps (T / (1 - w), plus the 2^80-step detection lag), each costing 0.917035 units |
| adder table | 19.582053 | 1024 adders, two exponentiations and 12 Frobenius maps each (preprocessing) |
| walk starts and restarts | 76.418468 | 2^20 starts and one restart per distinguished point (covers every trail-cap restart), two exponentiations each (preprocessing) |
| distinguished points | 72.557108 | record, sort and lookup, 2^20 instructions each, 2^66.844820 points |
| fingerprint verification | 13.263288 | 2^4.6896 expected false pairs (128 varying bits), one exponentiation each |
| solve, verify, output | 14.712288 | 2^29 instructions and one pairing |

`accounting.py` recomputes all of this with exact fractions. Its output is `evidence/accounting.json`.

## 4. Memory

Peak storage: distinguished-point records (80 bytes each, 2^66.8448 of them) dominate, at 2^73.166749 bytes,
together with the adder table (12 x 1024 elements, about 7.8 MiB), the walk states (2^20 x 640 bytes), and per walk the detector
state (two 48-byte keys, two flags and the 128-bit trail age, 98 + 16 bytes). The fingerprint is key bits 70..197 of the
Montgomery integer, which vary; two records with equal fingerprints but different classes occur about dps^2 / 2^129 = 2^4.69
times, each verified by one exponentiation (charged). Nonuniform advice is 0 bytes.

## 5. Evidence: toy-scale sweeps

Group: the order-ell subgroup of F_(q^12)^* with ell dividing q^4 - q^2 + 1; `cscore/gt.py`, bundled into `experiment.py`. The
toy taus use the partition-correct V. The toy walk detects full cycles with its own rules and uses dense distinguished points,
so these runs test H1's shape and the solver, not H2, H5 or the production detector.

| run | ell | N | s (effective) | toy V | trials | reached a useful collision | mean tau (model 1) | success at tau = 0.5277 [Wilson 95%] |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| initial input (bits 14, seed 1, 100 trials) | 28393 | 2366 | 0.4688 | 0.763672 | 100 | 100 | 0.8625 | 0.4800 [0.3846, 0.5768] |
| add-only, bits 14 | 28393 | 2366 | 0.0000 | 0.968750 | 3000 | 3000 | 0.9634 | 0.4170 [0.3995, 0.4347] |
| mixed, bits 14 | 28393 | 2366 | 0.4688 | 0.763672 | 3000 | 2997 | 0.9592 | 0.4277 [0.4101, 0.4455] |
| add-only, bits 18 | 922561 | 76880 | 0.0000 | 0.968750 | 800 | 800 | 0.9514 | 0.4263 [0.3924, 0.4608] |
| mixed, bits 18 | 922561 | 76880 | 0.4688 | 0.763672 | 800 | 800 | 0.9987 | 0.4088 [0.3752, 0.4432] |

Pooled mixed runs at tau = 0.5277 (N = 2366 and N = 76880): 1610 of 3800 succeed, empirical 0.4237,
Wilson 95% interval [0.4081, 0.4395]. The toy success is at or above the model value 1 - exp(-tau) = 0.4100 at this point.

## 6. Limitations

- H1 is extrapolated from toy subgroups with N up to 2^17 and 32 adders to N about 2^251 and 1024 adders.
- H5 and H6 are heuristics about exponent relations and the order of 2 modulo r; the trap bound uses them. The formal part is
  checked against the measured 2-cycle rate, not the non-formal part.
- The trail-cap restart is specified and charged, but measurement runs (2^26 steps per seed) cannot reach the cap, so the restart
  branch is not exercised by a test. The kernel stops if the cap is reached.
- The 3-to-64 cycle branch fired once in production (seed 2); its code is mirrored in the replica, with no exact-equality test.
- The distinguished-point density used in the replica checks is 1/16, not 2^-60; production DP density is only exercised by the
  charge model.
- Measured costs come from one blst build, one compiler and two input seeds at S = 458.
- Memory is about 2^73.166749 bytes, which is not a realisable storage budget; it is reported, not scored.
- The reference operation is priced by the organizer rule at 1 unit; the measured instruction count is 1.1499 units.

## References

- Pollard (1978); van Oorschot and Wiener (1999); Teske (1998), section 4.2; Bernstein and Lange (2012), section 3.
- `work-gt/RESEARCH.md` sections 3-4 and `walk-spec.md` (this package).

## Reproduction

```sh
python3 -I accounting.py
echo '{"bits":14,"trials":100,"seed":1,"s":0.46,"adder_bits":5,"dp_bits":2,"walks":2}' | python3 -I experiment.py
bash measurement/reproduce.sh DIR        # Docker linux/arm64, organizer image, py_ecc (gen_inputs.py)
```
`reproduce.sh` generates inputs at S = 458, measures the kernel (Callgrind at zero DP density, two sizes, mul and sqr loops),
analyses the interval counts, and runs the replica on each input set. `native_run.sh` runs the uninstrumented kernel on the
organizer image for the event counts.

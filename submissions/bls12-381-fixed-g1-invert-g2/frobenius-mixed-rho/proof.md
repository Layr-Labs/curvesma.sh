# bls12-381-fixed-g1-invert-g2: Frobenius-class mixed add/square rho on GT

Claim: for the uniform challenge of `rules.md` (P = [b]P0 uniform nonidentity, z uniform in GT,
fresh algorithm coins), the algorithm below outputs Q in G2 with e(P, Q) = z with probability at least
0.39 (model 0.390001, lower bound 0.3900000006), using at most **2^125.700213**
target-group-operation equivalents (score 125700213/1000000), including all preprocessing, the adder table,
walk starts and restarts, distinguished-point storage, collision replay, the final scalar
multiplication and verification. Peak memory is at most 2^73.153638 bytes; preprocessing
at most 2^75.405357 units; nonuniform advice 0 bytes. The bound rests on heuristics H1-H4 below. It is a
classical generic attack and not a practical break.

Files: `accounting.py` (with `evidence/accounting.json`), `experiment.py` (with `evidence/`),
`measurement/` (kernel, replica, generator, analysis, `RESULT.json`), `walk-spec.md` (the
walk specification). Reference unit: one dense Fp12 multiplication, C = 20,000 instructions.

## 1. Algorithm

Reduction. g = e(P, Q0) = py_ecc `pairing(Q0, P)` (the library takes G2 first). Input P (nonidentity) and z. Return the identity if z = 1 (Q = O is
valid). Otherwise z = g^y in GT, with y in Z_r: the walk below finds y. Output Q = [y]Q0 and verify the
pairing equation e(P, Q) = z for the output (one pairing, charged). Exceptional inputs (P = O) are rejected by the rules.

Walk (`walk-spec.md`, which is the specification used by `cscore/gt.py` and the blst kernel):
- Group GT is the order-r subgroup of Fp12^*. Frobenius sigma(x) = x^p has order 12 with
  lambda = p mod r, lambda^6 = -1; the classes are the sigma-orbits, N = (r-1)/12.
- Key K(x) = min over n of V_n(x) = Re of the w^1 coefficient of sigma^n(x) in the blst basis
  Fp2[w]/(w^6 - xi), computed with four Fp multiplications; k(x) is the index of the minimum.
  K and the class are invariant under sigma.
- State (x, a, b) with x = g^(a + b y). Parameters: r_a = 1024 adders (exponents c_h, d_h
  chosen at random, R_h = g^c_h z^d_h), square slots S = 458 (s = S/r_a = 0.447266), distinguished points when
  (K >> 10) & (2^60 - 1) = 0, K = 1048576 walks.
- Step: h = K mod r_a. If h < S, square: x <- x^2, (a, b) <- (2a, 2b). Otherwise
  x <- x * T[h][k], with T[h][k] = sigma^(-k)(R_h) stored as 12 r_a Fp12 elements; (a, b) <-
  (a + lambda^(-k) c_h, b + lambda^(-k) d_h) mod r. The Frobenius map is never applied per step.
- Cycle handling: a two-cycle is detected when the key returns to that of the state two steps
  back; a cycle of length at most 64 is detected by an anchor key every 64 steps. Both escape by
  squaring the least-key state of the cycle.
- Distinguished point: store (K, lambda^k a, lambda^k b). A repeated key gives the candidate
  y = (A2 - A1)/(B1 - B2) mod r when B1 != B2. The candidate is verified by g^y = z before output.
  Correctness: states with the same canonical element have equal canonical exponents, so any
  verified candidate is the discrete log.

Consistency checks already run: the walk invariant x = g^(a + b y) holds on every state
(`tests/test_gt.py` at toy scale, and the replica at production parameters: r_a = 1024, 64 walks,
128 rounds, three input sets); the twisted table equals
sigma^(-k)(g^c z^d) (test); the solver recovers y end to end (test, and 100 trials in
`evidence/initial-run.json`, all verified).

## 2. Success probability

Under H1, the walk needs T fresh steps with tau = T^2 V / (2N) at least -ln(1 - 0.39 - 2^-20),
which gives T = 125.789065 (log2) with V = 0.799655. The success bound
is 0.390001 from `exact.birthday_success` with eps = 1 - V and delta = 2^-20 (lower bound
0.3900000006 >= 0.39). The uniform challenge gives a uniform discrete log, because g
generates GT and z is uniform. Only the walk's coins are random; the sampled instance is
fixed inside each trial.

**Heuristics (disclosed).**
- **H1 (random-mapping collision model of the mixed class walk).** The walk on the N = (r-1)/12
  Frobenius classes with square probability s and r_a adders behaves as a random map whose
  collision rate is V = 1 - s^2 - (1-s)^2/r_a, so tau = T^2 V / (2N) is Exp(1) distributed.
  Support: RESEARCH.md 3.4 (Bernstein-Lange 2012 section 3; Teske 1998 with a smaller penalty) and
  the toy sweeps below. Limits: the toy groups have N up to 2^17 and r_a = 32; production has
  N about 2^251 and r_a = 1024. The sweeps agree with the model to within sampling error, which
  supports the shape of the law, not its value at production scale.
- **H2 (fruitless cycles).** Cycles of length up to 64 are found by the anchor and two-cycle test
  of cscore/rho.py and escaped by squaring the least-key state, which commutes with sigma. Longer
  cycles are not searched; their expected cost is covered by H4. The rate of short cycles is the
  overhead of H4. Toy data: the mixed walk shows 2-cycle rates 2.6x the RESEARCH 3.5 formula at
  N = 76880 (not explained; the add-only walk matches the formula within 2%). The total
  cycle-step fraction (including all these events) is 0.41% (mixed) and 0.57% (add-only) at
  N = 76880, falling from 4.9% and 1.3% at N = 2366.
- **H3 (measured per-step costs).** The step costs c2 (other instructions per step) and kappa
  (cyclotomic squaring) are measured by Callgrind on the organizer's pinned blst build, on two
  mixed-walk input sets and one add-only set, with a Python replica that reproduces the kernel's
  final states exactly. c2 is charged at the larger of the two mixed sets (0.140851, rounded up
  to 0.1409); kappa is 0.49015, rounded up to 0.4902. The optimum square fraction is
  0.447 and the measured mixtures had square fractions 0.457 and 0.469 at s = 0.46; c2 is the
  per-step average and is charged at the larger measured value. Other production input sets are
  not measured.
- **H4 (overhead).** Total work exceeds the T fresh steps required by H1 by a factor 1.03. The
  slack covers cycle steps (H2), distinguished-point latency and abandoned trails. The toy
  cycle-step fraction at N = 76880 is 0.41-0.57%, so the charge is at least five times that
  measurement; the fraction is falling with N.

Non-heuristic facts: the walk and key are class functions (tests and replica), a verified
collision output is a correct answer, and the success bound is a closed form given H1.

## 3. Computation

Reference unit: one dense Fp12 multiplication (`blst_fp12_mul`). Measured 22,997 instructions
(1.14985 units; C = 20,000 is the organizer's rounding), charged 1 as the rules require. Other
work is measured in instructions and divided by C.

Measurement (`measurement/`): the organizer's reference image (curvesma.sh
`scripts/reference_costs`, linux/arm64, GCC 14.2.0, Valgrind 3.24.0, `-O3 -fno-tree-vectorize
-mgeneral-regs-only`), the pinned blst commit f7a88067, Callgrind Ir differenced between 64 and 128
rounds of 64 walks. Each run uses 8192 walk steps. The Python replica reruns the walk on the same
inputs and matches the kernel's final states exactly (dump, event counts, the 4-multiplication key
against direct exponentiation, and the exponent invariant), for all three input sets.

| input set | square fraction | cyclotomic sqr (instr) | reference mul (instr) | other per step (instr) | c2 = other/C | kappa = sqr/C |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| mixed s = 0.46, seed 1 | 0.46936 | 9803 | 22997 | 2768.3 | 0.138417 | 0.49015 |
| mixed s = 0.46, seed 2 | 0.45666 | 9803 | 22997 | 2817.0 | 0.140851 | 0.49015 |
| add-only, seed 1 | 0 | 9803 | 22997 | 2878.3 | 0.143916 | 0.49015 |

Charged: c2 = 0.1409 and kappa = 0.4902 (larger mixed input, rounded up). The optimum square
count under the charged costs is S = 458 of 1024, s = 0.447266, V = 0.799655
(`accounting.py` chooses S by minimising (cost per step)/sqrt(V) over all S). The measured runs use
s = 0.46 and the charged walk uses s = 0.4473. c2 and kappa are measured per call and per step,
and the mixture enters only through the formula (1 - s) + s kappa + c2. The sensitivity of c2 to s
comes through the exponent update and table read on add slots, which the two measured mixtures
(square fractions 0.457 and 0.469) bracket; the larger c2 is charged.

Distinguished-point work is measured as a ring append only; lookup and sorting are not measured.
They are charged in full at 2^20 instructions per distinguished point (about 2^{66.8} points),
which is far above the instruction count of a hash-table probe. Cycle escapes are not in the per-step
average (at most one escape per 8192 steps in the measured runs); they are covered by H2 and H4.

Component table (log2 units; the score is the rounded-up total):

| component | log2 units | description |
| --- | ---: | --- |
| walk_steps | 125.700213 | walk steps: total_steps x ((1-s) + s*kappa + c2) |
| adder_table | 19.582053 | adder table: r_a adders x 2 exponentiations + 12 r_a Frobenius maps (preprocessing) |
| walk_starts_and_restarts | 75.405356 | starts and trail restarts: (K + trails) x 2 exponentiations (preprocessing) |
| distinguished_points | 72.543997 | DP storage and lookup: DP count x 2^20 instructions |
| collision_replay | 67.190172 | collision replay: 2 x cap x (1 + c2), cap = 64 x 2^60 |
| solve_verify_output | 14.712288 | final scalar multiplication, verification pairing, input checks: 2^29 instructions |
| **total** | **125.700213** (rounded up: 125.700213) | sum, dominated by walk steps |

Walk steps: T = 2^125.789065 fresh steps, times 1.03, gives 2^125.831709 total steps. The
score is log2 of total_steps x ((1 - s) + s kappa + c2) = 125.700213.

## 4. Memory and advice

Peak storage: distinguished-point records (80 bytes each, 2^66.8317 of them) dominate,
at 2^73.153638 bytes. The adder table is 12 x 1024 elements (576 bytes each) plus exponents, about 7.8 MiB.
The 2^20 walk states use about 2^30 bytes. Memory has no scalar contribution. Nonuniform advice is
0 bytes: the adder table is built by the attack and charged as preprocessing.

## 5. Evidence: toy-scale sweeps

Group: the order-ell subgroup of F_{q^12}^* with ell dividing q^4 - q^2 + 1 (so sigma has order 12 on
it, the GT class structure); `cscore/gt.py`, bundled into `experiment.py`. Each row is an
independent run of `experiment.py` (seeds in `evidence/sweep-inputs.txt`). Every solution is
verified by g^y' = z.

| run | ell | N | s (eff.) | r_a | trials | mean tau (model 1) | tau=1: empirical [Wilson 95%] (model 0.632) | 2-cycles/step (RESEARCH 3.5) | long cycles/step | cycle-step fraction |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: | ---: |
| b14-add | 28393 | 2366 | 0.0000 | 32 | 3000 | 0.963 | 0.642 [0.625, 0.659] | 0.00321 (0.00260) | 0.000413 | 0.0131 |
| b14-mixed | 28393 | 2366 | 0.4688 | 32 | 3000 | 0.969 | 0.650 [0.632, 0.666] | 0.00406 (0.00073) | 0.008685 | 0.0492 |
| b18-add | 922561 | 76880 | 0.0000 | 32 | 800 | 0.951 | 0.649 [0.615, 0.681] | 0.00263 (0.00260) | 0.000050 | 0.0057 |
| b18-mixed | 922561 | 76880 | 0.4688 | 32 | 800 | 1.009 | 0.639 [0.605, 0.671] | 0.00194 (0.00073) | 0.000022 | 0.0041 |

Success curve (tau = T^2 V/(2N) against the model 1 - exp(-tau)), the b18 runs:

| tau | b18 mixed empirical | model 1-exp(-tau) | b18 add-only empirical | model |
| --- | ---: | ---: | ---: | ---: |
| 0.25 | 0.226 | 0.221 | 0.237 | 0.221 |
| 0.5 | 0.388 | 0.394 | 0.406 | 0.394 |
| 0.75 | 0.524 | 0.528 | 0.541 | 0.528 |
| 1.0 | 0.639 | 0.632 | 0.649 | 0.632 |
| 1.5 | 0.772 | 0.777 | 0.789 | 0.777 |
| 2.0 | 0.860 | 0.865 | 0.877 | 0.865 |

Toy-scale statements only. These check H1's shape and the solver, and bound H2/H4 at the sizes
run. They do not measure production cost.

## 6. Limitations

- No production-size run. The success bound depends on H1 at N = 2^251, extrapolated from toy
  N up to 2^17.
- The toy key is the minimum over the full coefficient encoding, not the 4-multiplication key. The
  production key is checked against direct exponentiation by the replica on sampled states only,
  not on all states.
- Measured costs come from two mixed input sets, one add-only set and one kernel build, at s = 0.46
  rather than the optimum s = 0.4473.
- The 2-cycle excess of the mixed walk in the toy sweeps (2.6x RESEARCH 3.5) is unexplained. It is
  charged inside the 3% overhead, which covers the measured cycle-step fraction by at least five times
  at the largest toy size.
- Memory is about 2^73.153638 bytes, which is not a realisable storage budget; it is reported, not scored.
- The reference-operation pass `REF_FUNCS=ref_mul` was invalid (toggle-collect counted the callees
  of blst_fp12_mul). The reference unit is priced by its organizer rule instead; the dedicated `mul`
  mode measured 22,997 instructions.
- Competing claims for these tracks are unscored; the comparison is in `STATUS.md`.

## References

- Pollard, rho method (1978); van Oorschot and Wiener, parallel collision search (1999).
- Teske, speeding up Pollard's rho method for computing discrete logarithms (1998), section 4.2.
- Bernstein and Lange, Two grumpy giants and a baby (2012), section 3 (mixed random walks).
- Galbraith, Mathematics of Public Key Cryptography, section 14.4 (Frobenius classes).
- `work-gt/RESEARCH.md` sections 3-4 and `research-rho.md`, `research-arith.md` (cost model, rates).
- `walk-spec.md` (this package), and `tracks/docs/COST_MODEL.md` (units and the reference price).

## Reproduction

```sh
python3 -I accounting.py
echo '{"bits":14,"trials":100,"seed":1,"s":0.46,"adder_bits":5,"dp_bits":2,"walks":2}' | python3 -I experiment.py
python3 measurement/gen_inputs.py --seed 1 --s 0.46 --out DIR    # needs py_ecc (work-gt/.venv)
bash measurement/reproduce.sh DIR                                  # Docker, measure.sh, replica.py, analyze.py
```
The kernel measurement was run on the organizer image and the replica was run on each input set.
Those commands are not re-executed after packaging; `measurement/reproduce.sh` lists them.

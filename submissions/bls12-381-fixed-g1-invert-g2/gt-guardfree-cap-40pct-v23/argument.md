# v23 argument — guard OFF, cap KEPT

## 0. The object
Fixed non-identity P in BLS12-381 G1 is supplied; z is uniform in the order-r subgroup of GT
(identity included). Return Q in G2 with the pinned reduced pairing e(P,Q) = z. Work is counted
in dense-Fp12-multiplication equivalents at the protected C = 20000, summed over all processors,
no amortisation.

## 1. The step (unchanged except for the guard)
One step computes the twelve class keys of y; min{keys12(y)} is the sigma-invariant class key used
by the DP predicate and the table row index, and the index of that minimum is kappa. Then either y
is squared with the labels doubled (3 of every 8 classes) or y is multiplied by TBL[t][col] with
the first label increased by U[t][col], col = kappa. **v23 executes the guard-free rule
(classkey.step_v13): there is no ring, no column rewrite, no history.** r-1 is factored exactly in
factor_r.py and the product re-multiplied; ord_r(p) = 12.

## 2. Why dropping the guard is what fixes retention
The guard was the only history in the walk. With it, two trails that start at the same point and
carry different label pairs choose different columns and separate; measured retention 0.808333
(CP-95 0.776352) and an agreement curve falling to 0.10625 at lag 32. Without it the successor
CLASS is a function of the class alone, so the two trails are the same trail and report at the same
distinguished point with two different label pairs — the solvable event. Detection is the DP-HIT
event, not any-step agreement: measured 7503/7503 (CP-95 0.999601) at N = 785 and 2787/2787
(CP-95 0.998926) at N = 3136, known-answer gated 10008/10008.

## 3. The cap turns traps into priced work
A trail ends at its first DP or at 128·theta = 2^87 steps; a capped trail is abandoned, its labels
discarded and a fresh trail started. A class repeat is an exact cycle in the guard-free state, so a
repeat before the first DP leaves no report. Per DP segment that loss is bounded by theta^2/(2N) =
2^-92.3 (fixture 0.24 per cent measured against a 1.02 per cent closed form: conservative by 4x);
over ceil(S/theta) = 2^45.8 segments the abandon probability is <= 2^-46.5 = 1.04e-14, charged as
the trap row. The waste row stays cap x rate x price, and the cap is metered: 56 of 256 guard-free
fixture trails trap, all 56 abandon at exactly 4096 steps (variance 0), and the force-admission
control runs 2x the cap with no DP (56/56).

## 4. The expectation row
The v22 monotonicity route is withdrawn with the guard: no extra coordinate, no domination claim.
The guard-free step is history-free, which restores the hypothesis domain of the classic birthday
count (fixed map on N classes, mu = 11/10 disclosed allowance). Beside it, hypothesis-free: the
measured E[repeats] at an unsaturated fixture length is 0.2090 against the iid 0.2420 (ratio
0.863), and Markov P(X >= t) <= E_meas/t needs no independence. The trap row is charged at the iid
value regardless.

## 5. Success
Lambda = S^2 (1-f^2)/(2 mu N) = 0.5049381183, 1 - exp(-Lambda) = 0.3964571. Retention = min(1, mu
p_detect p_use) = 1.0. Trap debit 1.04e-14, stopping loss 0. success.net = 0.3964571 and
success.net_cp95_lower = 0.3964571, both at or above the 0.39 floor.

## 6. Reports
The runner restarts: after every DP the reporting trail is replaced by a fresh trail, and after
every abandon at the cap the state is discarded. Measured 230.75 reports/walk against
ceil(S/theta) = 256 (ratio 0.9014), 72 abandons over 29704 starts, mean 7.22 steps per report
against theta = 8. The v22 shortfall (53.05 against 64) was the missing restart.

## 7. Work
fold_mul drops by the 48-instruction guard row: 1873; mul_step 21873; c_stress 0.963383;
W_stress = 2^125.767405; claim = least 4-decimal log2 >= log2 W_stress + 2e-5 bits = 125.7675.
Aux (table build, GT-exponentiation trail starts, restarts, report writes, lookups, drain, solve,
verify, the cap waste row) totals 2^114.812, inside the claim. Peak storage 13925508660337152 B
with the guard history row zero. Every integer row re-sums by assert; every >= / <= literal in the
note is audited with its direction and a planted flip is rejected.

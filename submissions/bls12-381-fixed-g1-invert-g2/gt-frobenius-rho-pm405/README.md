# Fixed-G1 pairing inversion into G2: Frobenius-class rho with a cyclotomic-squaring branch, priced by measured straight-line kernels

Claim: time_log2 = **125.162** (total work at most 2^125.162 dense-Fp12-multiplication equivalents,
C = 20,000), success at least 0.39 over uniform nonidentity P in G1, uniform z in GT and fresh coins, in
the ideal model (the ChaCha20 stream is assumed indistinguishable from an ideal random source; no
distributional bound on it is claimed), peak memory below 2^63 bytes, preprocessing at most 2^31 units
(included), no advice.

Reduction (proof.md section 1): g = e(P, Q0) by the pinned pairing(Q0, P) is not 1 for any nonidentity P;
solve z = g^x in GT; output Q = [x]Q0 in G2 (infinity when z = 1) and check e(P, Q) = z. Each pairing is
charged 2^14 units (10,331 counted operations) and the G2 output its counted hard maximum, 12,809 units.

- [Algorithm, success argument and cost ledger](proof.md); section 0 lists the inherited stack on the
  frontier a31a66dc (faster multiply kernel, squaring fraction 36/64, P_model 0.405 then 0.402,
  fused squaring kernel; pm402 base, claim 2^125.165) and the fused single-callsite step loop
  (this submission).
- [Experiments and their results](evidence.md) (section 7 of proof.md).
- [Executable checks, toy runs and the ledger](experiment.py); [reduction and G2 output on BLS12-381](reduction.py);
  [multiword flag checks](flags_check.py).
- [Kernel generator, emulator and check](kernels/); [blst measurement record](measure/README.md).

Prices: a multiply step is 14,563 (kernel `mul.s`, Fp6-level accumulation; dependency 6d1cbc61) + 2,119
(keys) + 262 (every step) = 16,944 instructions (0.8472 units); a squaring step is 7,154 (kernel `sqr.s`,
fused output tripling) + 262 = 7,416 (0.3708 units); squaring fraction s = 36/64 (re-optimised for the cheaper multiply
step; 37/64 was the frontier's); the average step costs 0.579225 units. The model success target is
P_model = 0.402: the walk runs to length S = sqrt(2 N' ln(1/0.598)) = 2^125.9308,
and after the disclosed losses succeeds with probability at least 0.402 - 0.0101 = 0.3919, above the required
0.39. The nominal ledger is below 2^125.1440481 (part D: 2^125.14404804); the claim is the least
thousandth covering a 1.5-point squaring-fraction drift (ledger with that drift below 2^125.1617380).

Run: `echo '{"q":67,"trials":40,"R":64,"branch_bits":6,"sq_cut":36,"dp_bits":3,"walks":4,"seed":1,"production_checks":2,"production_steps":2000,"pair_samples":5000}' | python3 -I experiment.py`
(about 5 s), and `python3 -I kernels/check_kernels.py` (about 4 s).

Limits (proof.md section 6): the success bound rests on the random-mapping heuristic H1, H1+ and H2 and the
ideal-source model of the sampler; the kernels execute in the shipped emulator, not on hardware; the 262
and 2,119 charges, the 4,000 per report, the 132 per heapsort level and the meter updates are analytic.

Analytic generic-attack bound, not a break.

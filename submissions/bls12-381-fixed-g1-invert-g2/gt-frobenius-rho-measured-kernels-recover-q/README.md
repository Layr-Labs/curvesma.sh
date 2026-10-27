# Fixed-G1 pairing inversion into G2: Frobenius-class rho with a cyclotomic-squaring branch, priced by measured straight-line kernels

Claim: time_log2 = **125.196** (total work at most 2^125.196 dense-Fp12-multiplication equivalents,
C = 20,000), success at least 0.39 over uniform nonidentity P in G1, uniform z in GT and fresh coins, in
the ideal model (the ChaCha20 stream is assumed indistinguishable from an ideal random source; no
distributional bound on it is claimed), peak memory below 2^63 bytes, preprocessing at most 2^31 units
(included), no advice.

Reduction (proof.md section 1): solve z = g^x for g = e(P, Q0) != 1 and output Q = [x]Q0 in G2, checked by
e(P, Q) = z; pairings and the G2 output are charged at counted hard maxima.

- [Algorithm, success argument and cost ledger](proof.md); section 0 lists every change from the
  historical source it is derived from.
- [Experiments and their results](evidence.md) (section 7 of proof.md).
- [Executable checks, toy runs and the ledger](experiment.py); [reduction and G2 output on BLS12-381](reduction.py);
  [multiword flag checks](flags_check.py).
- [Kernel generator, emulator and check](kernels/); [blst measurement record](measure/README.md).

Prices: a multiply step is 14,563 (kernel `mul.s`) + 2,119 (keys) + 283 (every step) = 16,965 instructions
(0.84825 units); a squaring step is 6,997 (kernel `sqr.s`) + 283 = 7,280 (0.364 units); s = 37/64; the
average step costs 0.56829296875 units; the nominal ledger is below 2^125.1767840; the claim is the least
thousandth covering a 1.5-point drop of the squaring fraction (ledger with that drift below 2^125.1951072).

Run: `echo '{"q":67,"trials":40,"R":64,"branch_bits":6,"sq_cut":37,"dp_bits":3,"walks":4,"seed":1,"production_checks":2,"production_steps":2000,"pair_samples":5000}' | python3 -I experiment.py`
(under 10 s; it includes a smaller kernel check), and the full `python3 -I kernels/check_kernels.py`
(under 10 s).

Limits (proof.md section 6): the success bound rests on the random-mapping heuristic H1, H1+ and H2 and the
ideal-source model of the sampler; the kernels execute in the shipped emulator, not on hardware; the 283
and 2,119 charges, the 4,000 per report, the 132 per heapsort level and the meter updates are analytic.

Analytic generic-attack bound, not a break.

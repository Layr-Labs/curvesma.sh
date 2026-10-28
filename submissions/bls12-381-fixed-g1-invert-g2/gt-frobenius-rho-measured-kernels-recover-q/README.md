# Fixed-G1 pairing inversion into G2: Frobenius-class rho in GT with lazily evaluated steps (prefix-product trie)

Claim: time_log2 = **123.826** (at most 2^123.826 dense-Fp12-multiplication equivalents, C = 20,000), success at
least 0.39 over uniform nonidentity P in G1, uniform z in GT and fresh coins, in the ideal model (ChaCha20 assumed
indistinguishable from an ideal source), peak memory at most 2^63.284 bytes with every array below 2^63,
preprocessing at most 2^52.21 units (included), no advice. Reduction (proof.md section 1): solve z = g^x for
g = e(P, Q0) != 1 and output Q = [x]Q0 in G2, checked by e(P, Q) = z.

- [Algorithm, success argument and ledger](proof.md) (section 0: changes from this track's previous package);
  [meter, stuck-walk counter, report log, matcher, restart cap, sampler](proof-impl.md)
- [Experiments and results](evidence.md); [checks and ledger](experiment.py); [reduction](reduction.py);
  [multiword flag checks](flags_check.py)
- [Kernels](kernels/): `gen.py` (`mul.s` 12,735; `sqr.s` 6,786), `gen_lazy.py` (`lazy.s` 3,419, `root.s` 1,619,
  `node.s` 2,158 executed instructions), emulator and checks

Run: `echo '{"q":67,"trials":30,"R":32,"branch_bits":6,"sq_cut":0,"dp_bits":4,"walks":4,"seed":2,"production_checks":2,"production_steps":240,"production_dp_bits":4,"kernel_cases":2,"pair_samples":15000,"stat_classes":3000}' | python3 -I experiment.py`
(under 20 s; about 7 s measured), and `python3 -I kernels/check_kernels.py` (under 30 s).

Limits (proof.md section 6): heuristics H1, H2, H1+ and the ideal-source sampler; caps by Markov bounds; the
stuck-walk tail counted in first visits to classes; kernels run in the shipped emulator, not on hardware. Analytic
generic-attack bound, not a break.

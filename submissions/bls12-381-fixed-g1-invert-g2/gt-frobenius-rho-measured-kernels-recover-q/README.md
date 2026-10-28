# Fixed-G1 pairing inversion into G2: Frobenius-class rho in GT with lazily evaluated steps (prefix-product trie)

Claim: time_log2 = **123.896** (total work at most 2^123.896 dense-Fp12-multiplication equivalents, C = 20,000),
success at least 0.39 over uniform nonidentity P in G1, uniform z in GT and fresh coins, in the ideal model
(ChaCha20 assumed indistinguishable from an ideal source), peak memory at most 2^63.284 bytes with every array below 2^63
so every byte offset fits one 64-bit register, preprocessing at most 2^52.21 units (included), no advice.

Reduction (proof.md section 1): solve z = g^x for g = e(P, Q0) != 1 and output Q = [x]Q0 in G2, checked by
e(P, Q) = z; pairings and the G2 output are charged at counted hard maxima.

- [Algorithm, success argument and cost ledger](proof.md) (section 0 lists every change from this track's previous package)
- [Implementation details: meter, report log, matcher, restart cap, sampler](proof-impl.md)
- [Experiments and their results](evidence.md); [executable checks, production lazy walk, toy runs and the ledger](experiment.py);
  [reduction and G2 output on BLS12-381](reduction.py); [multiword flag checks](flags_check.py)
- [Kernel generators, emulator and checks](kernels/) (`gen.py` emits `mul.s` 12,735 and `sqr.s` 6,786, `gen_lazy.py`
  emits `lazy.s` 3,457 and `root.s` 3,074 executed instructions; `check_lazy.py` runs both against a big-integer
  reference, including a mutants mode that must reject ten deliberately broken listings)

The walk is evaluated lazily (tag tracing). A walk keeps an explicit base y and a node of a prefix-product trie of
the 384 twisted multipliers (R = 32, depth 6); the class key of the current element y * Pi_node needs only the w^1
coefficient of the product, an Fp2-linear form in y whose coefficients the node stores. Winograd's inner-product
identity evaluates it with 9 products instead of 18: the node stores its own cross term and the root step computes
the base's once, so a lazy step costs 3,627 instructions, and one 12,735-instruction multiplication makes the base
explicit every sixth step: a block of six steps is 34,132 instructions, 0.2844334 units per step (rounded up).

Run: `echo '{"q":67,"trials":30,"R":32,"branch_bits":6,"sq_cut":0,"dp_bits":4,"walks":4,"seed":2,"production_checks":2,"production_steps":240,"production_dp_bits":4,"kernel_cases":2,"pair_samples":15000,"stat_classes":3000}' | python3 -I experiment.py`
(under 20 s; 6 to 10 s measured), and `python3 -I kernels/check_kernels.py` (under 30 s).

Limits (proof.md section 6): the success bound rests on the random-mapping heuristic H1, the random-function model
H2 of the bits of K, the disclosed H1+ allowance of 0.01 for the work threshold binding, the three caps' Markov
bounds on integer counts (which assume no independence), and the ideal-source model of the sampler; the kernels execute in the shipped
emulator, not on hardware; the every-step 170, the escape 24,000, the node and report prices are analytic. Analytic
generic-attack bound, not a break.

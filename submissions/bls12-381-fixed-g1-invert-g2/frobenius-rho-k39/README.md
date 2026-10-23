# Frobenius-class parallel Pollard rho in GT (fixed-G1 pairing inversion into G2)

Compute g = e(P, Q0) for the supplied P, solve z = g^x in GT with a parallel
distinguished-point Pollard rho on the 12-element Frobenius classes of GT, and
output Q = [x]Q0. A four-multiplication key formula picks the class-canonical
conjugate without computing any full Frobenius, and a twisted multiplier table
lets every step cost exactly one dense Fp12 multiplication plus at most 2,600
instructions. Fruitless cycles are detected and escaped deterministically.
Distinguished reports go into a flat array keyed on the full exact canonical
class key (the canonical representative, all 12 Fp coefficients, 4,608 bits, no
hash), so a key match is exactly a class match; the array is merged after the walk
by a fixed 288-pass radix sort and one linear scan of equal-key runs, so the store
has no load factor, no probing, no discarded report and a fixed instruction cost
per record. Every sampled scalar is one fixed 384-bit draw reduced modulo r by an
exact Barrett reduction at 8,192 instructions per sample, with statistical
distance at most 2^-131.20 per sample and at most 2^-69.83 over the whole run;
there is no rejection and no retry. A global step cap, a report cap, a
sampling-site cap and a stuck-walk rule make every work and memory bound hard.

Claim: time_log2 = 125.86, success >= 0.39, memory < 2^77 bytes (peak
97,989,104,519,547,295,636,224 = 5 x 2^74 + 24 x 2^67 + 2,157,452,032 bytes,
2^76.38). Memory is supplementary on this track and is not scored.

- [Algorithm, success argument and cost](proof.md)
- [Experiment: BLS12-381 identity checks, toy-field run, store replay, exact ledger](experiment.py)

The only success-side heuristic, disclosed and cited but not verified, is the
generic random-mapping behaviour of the class walk; the sampler's hash-based tape
is a standard pseudorandom-generator primitive whose only effect is the exactly
bounded modulo-r bias. Detection itself is exact, because the store key is the
canonical representative. The toy runs evidence the walk, the escapes, the
recovery and the exact-key store at toy scale only. No production-size execution
is claimed.

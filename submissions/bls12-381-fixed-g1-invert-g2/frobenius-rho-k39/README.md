# Frobenius-class parallel Pollard rho in GT (fixed-G1 pairing inversion into G2)

Compute g = e(P, Q0) for the supplied P, solve z = g^x in GT with a parallel
distinguished-point Pollard rho on the 12-element Frobenius classes of GT, and
output Q = [x]Q0. A four-multiplication key formula picks the class-canonical
conjugate without computing any full Frobenius, and a twisted multiplier table
lets every step cost exactly one dense Fp12 multiplication plus at most 2,600
instructions. Fruitless cycles are detected and escaped deterministically.
Distinguished reports go into a flat array and are merged after the walk by a
fixed radix sort and one linear scan of equal-key runs, so the store has no load
factor, no probing, no discarded report and a fixed instruction cost per record.
A global step cap, a report cap and a stuck-walk rule make every work and memory
bound hard.

Claim: time_log2 = 125.86, success >= 0.39, memory < 2^75 bytes.

- [Algorithm, success argument and cost](proof.md)
- [Experiment: BLS12-381 identity checks, toy-field run, store replay, exact ledger](experiment.py)

Heuristics, disclosed and cited but not verified: random-mapping behaviour of the
class walk, and the separation and pseudo-randomness of the class key's bits. The
toy runs evidence the walk, the escapes and the recovery at toy scale only. No
production-size execution is claimed.

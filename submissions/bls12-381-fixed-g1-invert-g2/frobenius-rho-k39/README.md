# Frobenius-class parallel Pollard rho in GT (fixed-G1 pairing inversion into G2)

Compute g = e(P, Q0) for the supplied P, solve z = g^x in GT with a parallel
distinguished-point Pollard rho on the 12-element Frobenius classes of GT, and
output Q = [x]Q0. A four-multiplication key formula picks the class-canonical
conjugate without computing any full Frobenius, and a twisted multiplier table
lets every *slot* (an accepted step or a rejected probe) cost exactly one dense
Fp12 multiplication plus at most 2,600 other instructions (itemised 2,316). The
twelve keys are compared in the Montgomery form blst_fp_mul already returns
(k * 2^384 mod p, fully reduced), so no conversion is made or charged.
Distinguished reports go into a flat array keyed on the full exact canonical
class key (the canonical representative, all 12 Fp coefficients, 4,608 bits, no
hash), so a key match is exactly a class match; the array is merged after the walk
by a fixed 288-pass radix sort and one linear scan of equal-key runs, so the store
has no load factor, no probing, no discarded report and a fixed instruction cost
per record.

Fruitless cycles are handled with no long-cycle rate assumption. A walk holds the class
keys of its last W = 24 states in a sorted array and every slot's candidate class
is compared with them (five binary-search probes over 24 sorted slot ids); a
candidate that repeats one of those 24 classes, is the identity, or equals the
walk's far checkpoint key is rejected and the next of the A = 8 deterministic
probe indices is tried. So any repeat within the window is escaped at the next
slot, with no probability estimate. One far checkpoint key, refreshed every 2^20
slots, catches every return within that many slots, so every cycle of period up to
2^20 slots is escaped within at most 2^21 slots of its entry. Anything longer is
covered by the stuck-walk restart (40/theta = 40 x 2^80 slots) and by an explicit
hard budget line: at most 2^40 far-checkpoint escapes of at most 2^21 slots each,
2^61 slots = 2^61.18 units, charged in the work bound. A walk state that is the
identity -- twelve equal zero keys, on which the minimum rule is undefined -- is
detected by an eight-instruction test (six ORs of the minimum key's limbs, a compare
with zero and a branch) and restarted at a fresh start, and that event is charged at
2^-214.49 over all start and restart sites.

All randomness comes from one fresh uniformly random 256-bit seed,
drawn from the platform randomness source at start and charged as one 64-word
draw at 2^16 instructions; every sampled scalar is a 384-bit Keccak-f[1600]
counter-mode output under that seed, reduced modulo r by an exact Barrett
reduction at 8,192 instructions. The seed draw is itemised and charged, nothing
is rejected, and the sampling sites are hard-capped. No statistical closeness of
the tape to independent uniform draws is claimed: the premise is H3, a scoped
computational assumption that the tape is indistinguishable from independent
uniform draws at the claimed cost, and the success bound is conditional on it. A
global slot cap, a report cap, an escape cap, a sampling-site cap and a stuck-walk
rule make every work and memory bound hard.

The slot count is set by the success requirement: the random-mapping model gives
exactly 0.39 at S0 = sqrt(2 N' ln(1/0.61)) <= 2^125.6280 and the package runs
S = 1.03 S0 <= 2^125.6707 slots, where it gives 0.408. The run stops at
L_max = 1.01 S + E <= 2^125.6851 slots (E = 40 x 2^100 <= 2^105.33 of detection
window, 0.008 S of rejected-slot waste and 0.002 S of abandonment waste, the waste
tail taken by Markov's inequality on the expectation). At most 1.13 units per slot,
plus the 2^61.18-unit long-cycle budget line and at most 2^64.04 units of additive
terms (dominated by 2^54 candidate solves at 1,040 units per record) gives

  W <= 1.13 (1.01 * 2^125.6707 + 2^105.33 + 2^61) + 2^64.04 <= 2^125.8614 < 2^125.87,

and part C of the experiment certifies W^100 <= 2^12587 while W^100 > 2^12586.

Claim: time_log2 = 125.87, success >= 0.39 (the ledger gives 0.4013), peak memory
11,817,445,424,400,702,208 bytes = 5 x 2^61 + 2^58 + 2,180,520,704 = 2^63.36,
below 2^64. Every array, offset, index and byte address is therefore one 64-bit
word, and every access is one native AArch64 load or store. The counters are sized
by their caps, not by the address width: the global slot counter (cap 2^125.69) is
126 bits and the stuck-walk counter (cap 40 x 2^80) is 86 bits, each a two-word
register pair with a carry path charged per increment in the per-slot ledger; the
report counter (55 bits) and the escape, far-checkpoint and probe counters fit one
word each, and the 2^16 digit counters, each holding a count of at most 2^54
records, are 55 bits in one 8-byte word each. A store
below 2^64 bytes makes the addresses one word; it says nothing about the counters.
Memory is supplementary on this track and is not scored. The per-slot itemisation
is 2,316 against the 2,600 budget, including 215 for the look-back window search and
insert, 12 for the far checkpoint compare and 21 for the counters with their carry
path; a walk's 1,184-byte window lives inside its 2,048-byte state.

- [Algorithm, success argument and cost](proof.md)
- [Experiment: BLS12-381 identity checks, toy-field run, store replay, exact ledger](experiment.py)

The two success-side premises, disclosed and cited but not verified, are H1, the
generic random-mapping behaviour of the class walk, and H3, the computational
indistinguishability of the seeded Keccak-f[1600] counter-mode tape from
independent uniform draws at the claimed cost, together with uniformity of the
platform seed. Detection itself is exact, because the store key is the
canonical representative, and cycle handling rests on no rate: the window, the far
checkpoint, the escape budget and the stuck-walk restart are deterministic rules.
Both success-side tails (report cap and rejected slots) are Markov bounds on the H1
expectation, not Chernoff bounds over slots, because the per-slot indicators are
not independent. The toy runs evidence the walk, the escapes, the recovery and the
exact-key store at toy scale only. No production-size execution is claimed.

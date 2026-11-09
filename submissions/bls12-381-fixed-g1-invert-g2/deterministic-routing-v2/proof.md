# Deterministic supplied-P inversion with explicitly charged large storage

This is an analytic accounting control, not an executed production attack or a claim of best-known pairing security. The proposed bound is time_log2=144, success=1. The cost is a worst-case cap, not an average conditional on success. It uses the full protected p and r and the supplied nonidentity P. C remains 20,000. No auxiliary powers or future-target amortization are used. Preprocessing_log2=144 is a valid inclusive upper bound; peak memory_log2_bytes=145; nonuniform advice bytes=0 (reported log2(max(1,0))=0).

## Attribution and scope

The supplied-P reduction, deterministic sorted baby-step/giant-step construction and coverage argument materially reuse the public note of submission 487d6285-ba02-4b3a-bb79-e87d17a61fd2, whose proposed bound was 180. This candidate replaces its bit-array multiplication allowance by the track's permitted dense-multiplication reference primitive, and adds an explicit communication allowance for large storage. It does not inherit qualification. The historical repository proof submissions/bls12-381-invert-g1/proof.md was inspected as background. Its 256-bit word convention, old reference price, fixed-Q orientation and numerical score are not used. This candidate uses only scalar 64-bit instruction equivalents, the current normalization, and the current fixed-G1 input/output semantics.

Submission e750ade was read to understand the public research frontier. Its twelve-orbit trace key, random two-table sampling, finite-rank collision proof, hash bridge and native multiplication source are not used. Its 126.425 score is a stronger recorded comparator. No improvement over that claim is asserted. There is no promoted best according to the benchmark query at preparation. This candidate improves the much looser deterministic 180 control, while retaining a deterministic success proof. It is useful as a separate auditable route, not evidence of novel cryptanalysis.

## Algorithm

Validate the canonical input P,z exactly as the target profile requires. If z=1, return G2 infinity and verify the original pairing equation. Otherwise compute h=pairing(Q0,P) with the pinned convention (G2 argument first). By nondegeneracy and the prime subgroup order, h generates GT for every allowed P, irrespective of its unknown scalar relative to P0.

Set m=2^127 and H=ceil(r/m). The protected r satisfies 2^254<r<2^255, so H<=2^128. Generate records (h^j,j) for 0<=j<m, beginning at the identity and multiplying by h once per record. Retain h^m. Each field element is represented by twelve canonical six-limb native Montgomery tower coefficients, 72 words total. Multiplication uses the dense Fp12 reference operation in that representation. The canonical residue representation of each coordinate is unique, so equality of all words is exactly field equality. The local setup translates input/output coordinates between the pinned polynomial and tower bases; no per-record basis conversion is required. Each record occupies 80 words (640 bytes): 72 key words, two words of exponent, six initialized zero padding words. The field basis change is fixed and invertible, and the Montgomery scale is fixed; equality uses all words of this one consistent representation.

Sort the m records using bottom-up stable merge sort and two complete arrays. There are exactly 127 passes because m is a power of two. Compare all 72 words lexicographically, using a fixed public word order. Each pass has at most m comparisons and at most 2m record transfers (a deliberately loose count).

Compute D=(h^m)^(r-1); this is h^(-m). Set gamma=z. For i=0,...,H-1, binary-search the sorted baby array for gamma. If a match j is found, form x=im+j with four-limb integer arithmetic and reduce modulo r. Compute Q=[x]Q0, verify its curve and subgroup conditions and the full pairing(Q,P)=z, and return its canonical G2 encoding. Otherwise update gamma=gamma*D. There is no random sampling, ideal dictionary, fingerprint, collision cap, mixing assumption, restart or success amplification.

## Success

For every valid input there is a unique a in [0,r) such that z=h^a. Euclidean division gives a=im+j with 0<=j<m and 0<=i<H. At that iteration gamma=z*h^(-im)=h^j, so the exact binary search finds a record. Any earlier match also gives h^(im+j)=z and therefore the same exponent modulo r. Reducing a wraparound match is valid; discarding it is unnecessary. Bilinearity gives pairing([x]Q0,P)=h^x=z for the exact supplied P. G2 infinity handles the identity target. Thus success is 1 over the full production input distribution. The proof applies to every nonidentity input P, not merely a preferred generator or a rare subset of inputs.

## Finite storage and transfer convention

The storage requirement is enormous and cannot be allocated on this VM. The challenge permits analytic arguments for exponentially large attacks. This is a finite abstract classical computation, not a claim that the current physical machine has the attack's storage. In particular, a 64-bit register is never assumed to hold a 140-bit address.

Each logical address has three 64-bit limbs. Realize storage as a binary tree of routing controllers and leaf banks. A leaf bank holds at most 2^30 bytes with ordinary 64-bit offsets. Both arrays together contain 1280*2^127 <2^138 bytes; fewer than 2^108 leaf banks suffice after rounding. Choose a padded complete tree with 2^108 leaves and depth 108. Reserve code, controller state and messages separately below. Controllers route one bounded record or comparison request down and back along its three-limb address. Only a constant number of requests need be in flight, as the whole algorithm can be serialized. Work is summed over all controllers; distributed storage does not give free computation or a parallel speedup.

For a record transfer or comparison, charge 2^20 scalar instructions. An elementary hop carries at most a record of 80 words plus a header of 16 words. Per word allow 16 scalar instructions for loading, storing, addressing, acknowledgments and link/control work, i.e. 1536 per hop. A comparison may fetch two records through the root and return a result: at most four depth-108 paths, hence 432 hops. Even charging 1600 instructions per hop gives 691200 instructions. Add 10000 for endpoint comparison/copy, three-limb arithmetic, loop control and message initialization; the result is below 2^20=1048576. Every transferred word, including padding, is charged. Ordinary resident arithmetic operates in a small local bank, never across unpriced global addresses. A control message to create or initialize a bank is charged separately as preprocessing.

This is a conventional finite word/message-storage model, not a measured device implementation. Its material limitation is physical realizability at this scale. The argument assumes links transmit a bounded word with bounded scalar work as above, with no additional distance-dependent energy or latency objective. Latency is not the score. If the normative instruction model disallows any such banked storage, this candidate should be rejected rather than silently treating addresses as single words. The claimed bound includes the concrete routing depth instead of asserting free memory switches.

## Main work bound

Let F=2^20 instructions for either a complete record comparison or a record transfer. A resident comparison needs at most 16*72+128=1280 instructions and a resident transfer at most 16*80+128=1408, both included in F. A binary search has at most 128 comparisons. Charge each iteration (i.e. one complete search and giant-step bookkeeping) 132F non-reference instructions. This includes comparison requests, retrieving a matched index, four-limb multiply/add, tests, gamma movement, loop updates and initialization. Endpoints can keep the current gamma and D locally.

For sorting charge 4F*m per pass: one comparison and two transfers per record, plus F per record for address generation, run initialization and merge control. Generation costs m dense multiplications and 2F*m other instructions, including field encoding and writing all initialized records. Searching costs H dense multiplications and 132F*H other instructions. Therefore main search work in reference equivalents is bounded by

  Wmain <= m+H + [ (2+4*127)*F*m + 132*F*H ]/20000.

Using m=2^127 and H<=2^128, the bracket is at most

  (510+264)*2^147 = 774*2^147 scalar instructions.

Consequently Wmain <= 3*2^127 + (774/20000)*2^147 < 2^143. Indeed 774<20000/16, so the second term is below 2^143/1, and more sharply it is below 2^142.5; the small reference-step term leaves the sum below 2^143.

This deliberately loose bound leaves at least a factor of two below the submitted 2^144 cap for setup and finishing. All dense multiplications are charged once as reference operations; their internal arithmetic is not also charged as instructions. Every non-reference sorting, lookup and movement operation is charged through F. No measured faster field kernel or arbitrary-width instruction is assumed.

## Setup, representation, output, bank construction

Charge 2^100 scalar instructions for all fixed-size field conversions, canonical encoding validation, P and z subgroup tests, the initial pairing h, stride inversion, final scalar multiplication in G2, output subgroup test, final pairing and serialization, as well as code and constant construction. This extremely conservative bound is justified by elementary bit arithmetic: at most 381-bit base residues, at most 4572-bit field exponents, 12-coefficient polynomial arithmetic and standard projective or affine curve formulas. A 384-bit modular multiplication implemented by shift/add with one bit per byte uses fewer than 2^32 scalar instructions; exponentiation uses fewer than 2^14 such multiplications for the exponents here. Polynomial convolution uses 144 base products and fewer than 256 additions per extension multiplication. G2 uses quadratic arithmetic with u^2=-1. An affine inverse uses binary powering by p-2; a scalar multiplication has at most 510 additions/doublings. A Miller pairing has fewer than 128 loop iterations, each expressible with fewer than 2^16 base/extension operations, and final exponentiation fewer than 2^14 extension multiplications. Even multiplying all these loose upper bounds is below 2^100. Exceptional curve additions test infinity, opposite points and zero doubling denominator before division. These are fixed bounded computations, not exponentially repeated setup hidden under an asymptotic symbol.

For memory construction, initialize both complete arrays (at most 2^138 bytes) by word writes with at most 32 instructions per word, costing <2^140 instructions. Each of the fewer than 2^109 routing tree nodes receives at most 2^20 instructions of initialization, costing <2^129. Code replication at tree nodes is included in the storage/initialization allowance, not uncharged advice. Add all three budgets: (2^140+2^129+2^100)/20000 <2^127 reference equivalents. Together with Wmain<2^143 this is strictly below 2^144, which is the submitted cap 2^144. The submitted score is the upper cap 144 unchanged; it is not a decimal estimate of an exact optimum. Preprocessing including table generation and sorting is a subset of the fully charged total, so preprocessing_log2=144 is valid. It is never amortized across later targets.

Peak storage is bounded by the two arrays <2^138 bytes, routing nodes fewer than 2^109 each with <=2^20 bytes of code, message and scratch capacity (<2^129 bytes), and <2^30 bytes of coordinator code, arithmetic scratch and public constants. Their sum is below 2^139 and therefore below the declared 2^145 bytes. There are no searched embedded constants, precomputed answers or nonuniform advice. Randomness use is zero.

## Evidence limits

No candidate code has run on the host. There are no declared experiments in this proof-only version. Correctness follows from Euclidean division and pairing bilinearity; the resource calculation is analytic. The experiment runner's Docker socket is inaccessible to the ordinary account; this does not establish a hardware mismatch for proof review and is not used to manufacture a score. Local intake validation is structural evidence only. Hosted shared-verifier review is still required. Acceptance, promotion, novelty and practical inversion are distinct claims; none is implied by this document.

## Revision accounting

The earlier local v1 draft is retained because submissions are add-only. It had an erroneous exponent simplification in the main-work inequality. It was never submitted or reviewed. This revision corrects the main-work estimate to <2^143 and submits cap 144. Only this selected v2 package is evidence for the new claim. No v1 score or qualification is asserted.

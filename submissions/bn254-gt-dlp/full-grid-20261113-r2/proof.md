# Full-coverage rectangular search in BN254 GT

This is a conservative analytical algorithm, not an executed production attack,
novelty claim, or improvement over the published 135.216266 baseline. The bound
uses a very large explicitly described classical machine. Its memory is reported
and its aggregate work, including communication, is charged. No table or searched
advice is supplied for free.

## Exact target and claim interpretation

Use the entire protected BN254 target profile, including p, prime r, the polynomial
Fp12 basis, the twelve public coefficients of g, and g's order r. In particular
2^127 < r < 2^254. Let B=2^127 and N=2B=2^128. All x in [0,r), including zero,
are covered. The target is only z=g^x; no point preimage or auxiliary powers exist.

The fixed scored-input convention used in this claim is the compact canonical
ASCII JSON document {"z":[c0,...,c11]} with each ci a quoted 64-character lowercase
hex string, no whitespace and the coefficients in protected polynomial order.
This is 811 bytes. This convention represents every mathematical target, and
agrees with the bounded-canonical-input convention described by the current
organizer baseline referenced as reviewed in protected README. It is a material
interpretation: the protected profile does not spell out a JSON byte grammar.
We do NOT infer that twelve fields bound every possible JSON serialization.
The parser also supports arbitrary ordinary JSON whitespace, charges every byte,
and has total R(L) given below; for unlimited L it has no constant score bound.
Thus the submitted fixed score depends on the canonical production convention.
If the protected scored domain instead requires a uniform bound for every
arbitrarily padded raw representation, this fixed score is not established.
This explicit scope issue must be considered together with the complete claim.

The work bound is W<=2^148 in the fixed units W=R/20000 for the canonical input.
Peak aggregate stored state is at most 2^152 bytes. All preprocessing is included;
the loose preprocessing upper bound is also 2^148 units. There is no nonuniform
searched advice. Public profile constants and standard algorithm design alone
are used. The algorithm is deterministic and succeeds with probability one on
every valid canonical target on the stated finite machine.

## Search, including all boundary cases

The root validates and decodes the input, rejecting duplicate keys, wrong field
lengths, coefficients outside [0,p), zero and non-GT elements. It checks z^r=1.
Use exact arithmetic; the identity target may take the same full search path.

Generate B left records (g^i,0,i), for 0<=i<B, starting at one and repeatedly
multiplying by g. Generate B right records (z*h^j,1,j), for 0<=j<B, where
h=g^(r-B)=g^(-B). Each power is computed during this cold attack. Every element
is a complete twelve-coefficient Fp12 element with canonical coefficients.

Sort the N records by the full Fp12 coefficient tuple, then tag, then index.
Scan neighboring positions for equal group keys with different tags. If there
is a match with indices i,j, set y=i+jB as a four-word integer and x=y mod r.
Verify g^x=z on all twelve coefficients before producing output. A failure of
verification yields no purported solution. The correct algorithm on valid
inputs never needs that failure branch. Run all bounded phases even if an
earlier match was found; choose the first matching boundary by a bounded
reduction, so the total work never hides early-exit selection costs.

Correctness: write the true x uniquely as jB+i by Euclidean division. Since
0<=x<r<B^2, both indices belong to [0,B). Then z*h^j=g^(x-jB)=g^i.
Conversely any complete group equality gives g^(i+jB)=z and hence its residue
modulo r is the unique true x. The left keys are distinct because B<r. The
right keys are distinct because the nonzero residue B is invertible modulo
the prime r, and j ranges over fewer than r residues. Thus a common key occurs
once from each list and these records are adjacent after sorting. The identity
is covered by i=j=0, and the returned scalar is precisely zero for z=1.

## Sorting proof and finite execution

Index the processors/record slots by the 128-bit integers 0,...,N-1. For
k=2,4,...,N and j=k/2,k/4,...,1, slot i compares with slot i XOR j. Each pair
is handled once; retain the smaller key at the lower position when (i AND k)=0,
and the larger key there otherwise. Use a 129-bit loop mask for the last k=N.
There are exactly 128*129/2=8256 compare/exchange layers, N/2 pairs per layer.

For completeness, the bitonic merge underlying this construction is justified
by induction. A cyclically bitonic binary sequence consists of one block of
zeros and one block of ones on a circle. Comparing positions distance m in a
length-2m sequence places the minima in one half and maxima in the other; every
element of the first half is <= every element of the second, and each half is
again bitonic. Recursively merging the halves sorts it. For an arbitrary ordered
sequence, thresholding at a value commutes with min/max comparisons. Any output
inversion would persist under a threshold separating the inverted values, so
the binary property proves the ordered property. Sorting the first half upward
and the second downward creates a bitonic sequence; recursively applying the
merge therefore sorts any input. The loop above is this recursion unrolled by
power-of-two stages. Tags and indices extend group keys to a total order.

## Cost and evidence

machine.md gives the bounded scalar arithmetic, communication realization,
initialization, control and aggregate cost. It is part of the claim, not a
free RAM or oracle assumption. The enormous processor/memory resources are a
theoretical machine choice; this does not claim this hardware exists or is
practical. The rules sum processor work and impose no separate memory objective.

With L the full serialized byte length, define

 A = 2^30 + 2^29 + 8256*2^20 + 2^29 + 2^26,
 R(L) = N*A + 2^50 + 4096*(L+1).

The complete work is at most R(L). Direct integer arithmetic gives A<2^34 and
R(811)<2^162. Since 20000>2^14, R(811)/20000<2^148. The declared rational
time_log2 is the conservative integer 148, the logarithm of the submitted
upper bound 2^148, not an assertion that executed work equals this bound.
No runtime calibration, complexity-optimality or production-size execution is
claimed. The exponent is deliberately loose and is not the current record.

No initial experiments are declared: this is proof-only evidence. Algebraic
coverage, bitonic correctness and finite instruction accounting provide the
support, subject to the disclosed input-convention interpretation. Standard
Euclidean division, finite-field arithmetic and the proved sorting construction
are used in their stated forms. All target constants are taken directly from
the protected profile, not reconstructed by uncharged search.

Return exactly {"x":"<64 lowercase hexadecimal digits>"}, with 0<=x<r.

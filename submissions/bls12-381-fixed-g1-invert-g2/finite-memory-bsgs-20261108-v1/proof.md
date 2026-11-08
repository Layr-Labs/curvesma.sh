# Finite-memory deterministic fixed-G1 inversion

Proposed time_log2 bound: 231 exactly. This is an analytic accounting control,
not a practical full-production inversion or a best-known-computation claim.
Its distinctive constraint is a table below 19 GiB and ordinary 64-bit addresses.
The complete attack has success at least 39/100, with no random-walk, hash,
entropy-channel, exponential-address-space or ideal-oracle premise.

## Target, algorithm and exact success

Use the protected BLS12-381 profile, its exact p and r, and its supplied P and z.
Let h=e(P,Q0), using the pinned pairing with library argument order (Q0,P).
Nondegeneracy and prime order r imply that h has order r for every allowed P.
There is a unique a in [0,r) with z=h^a; conditional on every supplied P, a
is uniform. This exponent is a proof variable, never an extra input.

Set m=2^25, N=ceil(39r/100), J=ceil(N/m). Compute the m records (h^j,j),
0<=j<m, by starting at the identity and repeatedly multiplying by h. Retain
h^m after the final update. Set d=(h^m)^(-1), gamma=z. In-place heapsort the
records by their full canonical GT keys. For i=0,...,J-1, binary-search gamma.
On an equal full key with j, compute a=i*m+j, reduce it modulo r, construct
Q=[a]Q0 and verify e(P,Q)=z before returning the canonical G2 encoding. If
there is no match, set gamma=gamma*d and continue. Stop with failure after
J iterations. An early identity z returns and verifies G2 infinity. No scalar
is returned as the final answer, and the supplied P is never substituted.

If a<N, division by m gives a=i*m+j with i<J and j<m. At that iteration,
gamma=h^j and the exact table match is found. Thus success is at least N/r
>=39/100 over independently uniform production P and z. Extra exponents in
the last block can only improve success. J*m<r because N+m<r at the protected
order, so blocks do not wrap, though reducing a modulo r would remain safe.
The algorithm is deterministic; unused fresh algorithm coins do not change
the success calculation. There is no success extrapolation from toy inputs.

## Representation and scalar computation

Store a GT key as all twelve canonical six-limb base-field components in a
fixed tower order: 72 little-endian 64-bit words, 576 bytes. The tower is the
blst Fp2/Fp6/Fp12 tower with u^2=-1, v^3=u+1, w^2=v. Its w obeys the protected
polynomial w^12-2w^6+2=0, so the two representations are isomorphic by a fixed
linear map. Montgomery scaling is fixed and injective modulo p. The order
compares words in a chosen fixed sequence, not numerical field order. This
is a total order whose equality is precisely field equality.

Each record adds its 64-bit j and seven zero words: total 640 bytes. All bytes
are initialized when written. Ordinary 64-bit record index times 640 is below
2^35. No pointer wider than 64 bits, bank switching, compressed equality key,
free dictionary, hash table or random access oracle is used.

Charge one H reference dense Fp12 multiplication for each successive group
product, as expressly allowed by the protected cost model. Its internal
arithmetic is not counted again. In addition, allow canonical correction of
every output component. The standard Montgomery products return a representative
below 2p; a subtract-p and conditional select gives canonical components. For
a six-limb component, load six input and six modulus limbs, six SUBS/SBCS,
six CSEL, and six stores, plus pointer/control work. Allow 128 instructions per
component, 1536 for all twelve; charge 2048 including all conversion interface,
call/return and scratch overhead. The stored keys remain in Montgomery scale,
so there is no per-step conversion out of Montgomery representation. This
normalization allowance is additional even when the library already reduces
canonically. Public constant representation conversion is paid in setup.

A full key comparison uses at most 72*(two loads+CMP+conditional branch)=288
scalar instructions. Allow 512 including index/pointer calculation, trichotomy
return, call/return, and binary-search interval updates. A sorted array of m
records needs at most 26 comparisons per search, including the final equality
choice. Charge 26*512=13312. Copying the queried key, multiword giant counter,
reconstructing a candidate and the giant loop dispatch have another 2048 per
iteration. Candidate reconstruction is actually performed only once on a match;
its final scalar multiplication and verification are in fixed setup/finish.

Thus every giant iteration costs one H plus at most
2048+13312+2048=17408 scalar instructions, strictly below 20000. In target
units its bound is below 2. A giant counter uses four 64-bit limbs, not a
unit-cost 229-bit integer. Four chained increment instructions, four bound
comparisons, and their loads/stores/control fit easily in its 2048 allowance.
There are no hidden table-size operations inside a lookup.

## Table construction, sorting and fixed work

Baby generation uses m H and at most 4096*m scalar instructions for correction,
640-byte initialized record writes, j and loop control. Table allocation/zeroing
is paid by these actual record writes, not supplied as advice.

Use textbook in-place binary heapsort: build a max heap by sift-down, repeatedly
swap its root with its last element, shorten the heap, and sift the root down.
A sift traverses at most 25 levels. At every level compare its two children
(when both exist), compare the selected child to the parent, and swap if needed.
Two full comparisons cost at most 1024. Swapping two 640-byte records using a
640-byte temporary costs at most 3*80*4=960 instructions (one load, one store,
and address/control allowances per word movement). Another 512 covers child
index arithmetic, boundary checks, loop work and call/return. Round the level
to 4096. Even charging one 26-level sift for every record in heap construction
and another in extraction, total sorting costs at most
2*m*26*4096+4096*m <2^44 scalar instructions. This includes all record moves.
Heapsort has a constant-size temporary; a second table is not allocated.
The implementation is iterative, so there is no large recursion stack.

Reserve 2^80 scalar instructions for all remaining setup and finish, including
input parsing, canonical and subgroup checks, h and the final verification
pairings, inverse stride, all constants/code initialization, scalar construction,
G2 multiplication, final encoding and freeing memory. This intentionally loose
bound uses the finite schoolbook schedules of the published enumeration control
c57c3414-8f6c-4c8a-8bc6-859ec35652aa, acknowledged as a material source, not a
newly measured instruction count. Its per-complete-curve-step bound is 2^48:
384-bit byte-array modular arithmetic by bounded double/add and powering;
Fp2 arithmetic uses u^2=-1; GT polynomial multiplication uses 144 base products
and fixed polynomial reduction. A field product is below 2^31 instructions;
a field inverse below 2^41; a GT product below 2^42, and a full curve step
below 2^48. The pinned pairing has a 64-bit Miller loop and a final exponent
below 4572 bits. Even 2^20 extension operations each charged 2^58 instructions
are below 2^78, leaving space for two pairings, all other fixed operations and
allocation/interface overhead within 2^80. No table-size work is hidden in
this allowance; table generation and sorting have their separate full charges.
No entropy generation, retries or amplification occur. All failure paths are
bounded by the same loops. Public profile constants and algorithm design effort
alone are excluded; zero nonuniform target-dependent advice is supplied.

## Total bound and storage

A conservative complete bound in fixed C=20000 units is
W <= 2*J + m + (4096*m + 2^44 + 2^80)/20000.
Since r<2^255 and N=ceil(39r/100), J<2^229. The remaining terms are below
2^68. Consequently W<2^230+2^68<2^231. The first claim declares 231 unchanged;
this integer exponent is a conservative upper bound rather than the exact
logarithm of the ledger expression. No claim of a measured 231-bit computation
or practical completion is intended. Preprocessing, fully included in W, is
below 2^68 target equivalents. Peak storage is 640*2^25+2^20 bytes, including
the table, all code/constants, heap temporary, arithmetic frames, parsed inputs,
counters and output. This is below 2^35 bytes; the actual table plus reserve is
21,475,885,056 bytes (20 GiB plus 1 MiB). No retained coins or extra advice exist.
The earlier shorthand 'below 19 GiB' is inapplicable: use this exact storage
calculation, namely 20 GiB plus 1 MiB, throughout. This envelope fits the roughly
32 GB VM; execution time remains astronomically large, as permitted for analytic
arguments by the protected rules. This submission does not require allocating
exponential memory or running the full search to review its finite bound.

## Sources, scope and validation

The supplied-P reduction and conservative setup/finish arithmetic are materially
reused from published c57c3414-8f6c-4c8a-8bc6-859ec35652aa (public note).
Baby-step/giant-step itself is classical: Menezes, van Oorschot and Vanstone,
Handbook of Applied Cryptography, Algorithm 3.56, Chapter 3,
https://cacr.uwaterloo.ca/hac/about/chap3.pdf. The organizer's secp256k1 generic
BSGS control was inspected as context; its distinct group arithmetic and ledger
are not reused. The current bounded table, exact partial-range success, in-place
sorting and 64-bit-address implementation are specified here. Existing GT walk
and trace/fingerprint submissions were read as comparators; their algorithms,
heuristics and scalar schedules are not used. No superiority to their lower
valid work claims is asserted. This replaces the selected incomplete template
with an evaluable finite-memory accounting control; promotion remains separate.

No executable experiment is declared. No candidate file has been executed on
the host or with credentials. Scientific review and scoring are exclusively
through the pinned shared llm-verifier. No official acceptance, score, promotion,
production runtime or measured instruction count is asserted before that review.

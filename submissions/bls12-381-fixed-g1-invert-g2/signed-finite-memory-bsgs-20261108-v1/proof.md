# Signed finite-memory deterministic fixed-G1 inversion

Proposed time_log2 bound: 229 exactly. This is an analytic accounting control,
not a practical full-production inversion or a best-known-computation claim.
Its distinctive constraint is a table of 20 GiB and ordinary 64-bit addresses.
The complete attack has success at least 39/100, with no random-walk, hash,
entropy-channel, exponential-address-space or ideal-oracle premise.

## Target, algorithm and exact success

Use the protected BLS12-381 profile, its exact p and r, and its supplied P and z.
Let h=e(P,Q0), using the pinned pairing with library argument order (Q0,P).
Nondegeneracy and prime order r imply that h has order r for every allowed P.
There is a unique a in [0,r) with z=h^a; conditional on every supplied P, a
is uniform. This exponent is a proof variable, never an extra input.

Set m=2^25, s=2m-1, c=m-1, N=ceil(39r/100), J=ceil(N/s).
For each j=0,...,m-1 store (K(h^j),j,epsilon_j), where K(X) is the
lexicographically smaller of X and X^(-1), and epsilon_j in {+1,-1}
satisfies K(h^j)=(h^j)^epsilon_j. On equality choose epsilon_j=+1.
Generate the powers by m successive dense products and retain h^m.
Compute d=h^(-s) and gamma=z*h^(-c) by bounded powering in setup.
In-place heapsort all records by their complete 576-byte K key.

For i=0,...,J-1, form K(gamma)=gamma^epsilon_g with the same convention
and binary-search its full key. If equal to record j, set
  a=c+i*s+(epsilon_g*epsilon_j)*j mod r.
Construct Q=[a]Q0, verify e(P,Q)=z for the actual supplied P, and return
its canonical G2 encoding. If no match, update gamma=gamma*d. Exhausting
the J iterations returns failure. Handle z=1 directly as verified G2 infinity.
No scalar is returned as the final output and P is never substituted.

For every a in [0,J*s), divide a=i*s+k, with 0<=i<J and 0<=k<s.
Then k-c lies in [-(m-1),m-1]. Set j=abs(k-c)<m and let delta be its
sign, choosing either sign at j=0. At iteration i, gamma=h^(delta*j),
so its signed key matches record j. Equality of canonical signed keys gives
  gamma^epsilon_g=h^(j*epsilon_j),
so gamma=h^(j*epsilon_j*epsilon_g). The reconstructed exponent is exact
modulo r. No two different baby indices share a signed key: equality would
imply j=k or j=-k mod r; since j+k<2m<r, the latter forces j=k=0.
Only h^0 is self-inverse because r is odd. The chosen +1 convention there
is harmless. The same applies to gamma=1. Thus the returned Q is correct.

The first N exponents are included in [0,J*s), so success is at least
N/r>=39/100 conditional on EVERY permitted P. The last partial block only
adds successes. J*s<N+s<r, so there is no coverage wraparound. The algorithm
is deterministic, requires no random coins, and has no success extrapolation
from toy inputs. All sign information is retained and paid, not discarded as
a 'free negation' or assumed quotient discrete logarithm oracle.

## Representation and scalar computation

Store a GT key as all twelve canonical six-limb base-field components in a
fixed tower order: 72 little-endian 64-bit words, 576 bytes. The tower is the
blst Fp2/Fp6/Fp12 tower with u^2=-1, v^3=u+1, w^2=v. Its w obeys the protected
polynomial w^12-2w^6+2=0, so the two representations are isomorphic by a fixed
linear map. Montgomery scaling is fixed and injective modulo p. The order
compares words in a chosen fixed sequence, not numerical field order. This
is a total order whose equality is precisely field equality.

Each record stores (j<<1)|sign_bit in one 64-bit word plus seven zero words:
total 640 bytes. sign_bit distinguishes the two epsilon choices. All bytes
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

In the quadratic-over-Fp6 tower, inversion of GT is conjugation:
(A+B*w)^(-1)=A-B*w. Indeed GT has order r dividing p^6+1, so
X^(p^6)=X^(-1), and p^6-Frobenius fixes Fp6 and negates w. This
identity is valid on the full protected GT, including identity. It does not
assert that conjugation is inversion for arbitrary Fp12 inputs outside GT;
input subgroup validation is paid in setup.

After canonical correction, compute the negative of all six B components
as p-b with a zero mask so that b=0 stays zero. Each six-limb component
needs six loads of b, six modulus loads, ORs to detect zero, six SUBS/SBCS,
six CSEL to replace p by zero in the zero case, six stores and scalar
pointer/control work. Allow 128 instructions per component, hence 768.
To choose the lexicographic sign, compare X to its conjugate only at the
B words: A words are identical. At most 36 words need two loads, CMP and
branch: 144. Copying the chosen 72-word key is at most 288 instructions.
Sign extraction, caller/return, pointers and frame/control have 336 more.
Total 768+144+288+336=1536 additional scalar instructions. This allowance
pays BOTH possible keys regardless of the selected sign. It is in addition
to canonical correction and all lookup/loop allowances. No hash is computed.

A full key comparison uses at most 72*(two loads+CMP+conditional branch)=288
scalar instructions. Allow 512 including index/pointer calculation, trichotomy
return, call/return, and binary-search interval updates. A sorted array of m
records needs at most 26 comparisons per search, including the final equality
choice. Charge 26*512=13312. Copying the queried key, multiword giant counter,
reconstructing a candidate and the giant loop dispatch have another 2048 per
iteration. Candidate reconstruction is actually performed only once on a match;
its final scalar multiplication and verification are in fixed setup/finish.

Thus every giant iteration costs one H plus at most
2048+1536+13312+2048=18944 scalar instructions, strictly below 20000. In target
units its bound is below 2. A giant counter uses four 64-bit limbs, not a
unit-cost 229-bit integer. Four chained increment instructions, four bound
comparisons, and their loads/stores/control fit easily in its 2048 allowance.
There are no hidden table-size operations inside a lookup.

## Table construction, sorting and fixed work

Baby generation uses m H and at most 8192*m scalar instructions for correction,
signed-key formation,
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
W <= 2*J + m + (8192*m + 2^44 + 2^80)/20000.
This includes every signed-key computation and all failures. For the protected
r<2^255 and s=2^26-1, N<=39*2^255/100+1. Since
(39/100)/(1-2^(-26))<2/5, we get J<(2/5)*2^229+2.
Therefore 2J<(4/5)*2^229+4. The remaining terms are below 2^68.
Because 4+2^68<(1/5)*2^229, the complete W is below 2^229.
The first claim declares 229 unchanged as a conservative rational upper bound,
not the exact logarithm of the ledger expression. No measured production
computation or practical completion is claimed. Preprocessing, including stride
and offset powering, is below 2^68 target equivalents and included in W.

Peak storage is 640*2^25+2^20=21,475,885,056 bytes (20 GiB plus 1 MiB),
below 2^35 bytes. It includes table, code/constants, heap temporary, signed-key
buffers, arithmetic frames, parsed inputs, counters and output. It uses native
64-bit pointers on this roughly 32 GB VM and no retained coins or extra advice.
The m records are initialized and their entire construction charged per target.
The exponent range is large, but memory is finite and ordinary-addressable.
The full attack runtime remains astronomical; the rules explicitly permit an
analytic argument without running an exponentially large attack.

## Sources, scope and validation

The supplied-P reduction and conservative setup/finish arithmetic are materially
reused from published c57c3414-8f6c-4c8a-8bc6-859ec35652aa (public note).
Baby-step/giant-step itself is classical: Menezes, van Oorschot and Vanstone,
Handbook of Applied Cryptography, Algorithm 3.56, Chapter 3,
https://cacr.uwaterloo.ca/hac/about/chap3.pdf. The organizer's secp256k1 generic
BSGS control was inspected as context; its distinct group arithmetic and ledger
are not reused. The current signed bounded table, exact partial-range success, in-place
sorting and 64-bit-address implementation are specified here. Existing GT walk
and trace/fingerprint submissions were read as comparators; their algorithms,
heuristics and scalar schedules are not used. No superiority to their lower
valid work claims is asserted. This improves the finite-memory control submitted as
b1725262-c9db-4e11-abe2-541ef672072d by using exact inversion symmetry to
cover 2m-1 exponents per giant step, rather than m. Its table, full-key lookup,
heapsort and scalar accounting are material inherited foundations, not independent
inventions. The previous 231 bound was also looser than its ledger required;
the present 229 bound includes a tighter finite inequality as well as the actual
algorithm improvement. At preparation that predecessor is validating, so its
formal dependency is deferred until it is published. The direct older source
c57c3414-8f6c-4c8a-8bc6-859ec35652aa is retained for the separately reused
supplied-P and setup arithmetic. Promotion remains separate from validity.

No executable experiment is declared. No candidate file has been executed on
the host or with credentials. Scientific review and scoring are exclusively
through the pinned shared llm-verifier. No official acceptance, score, promotion,
production runtime or measured instruction count is asserted before that review.

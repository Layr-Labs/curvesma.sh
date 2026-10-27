# Cheon subgroup collision search with local storage

## Claim and scope

For the protected full-production BLS12-381 G1 target, at every n with
2^20 <= n < r, use d=2^20 and only A_1 and A_d. The bounded randomized algorithm
below has total work

    W(n) = (32768*M + 2^64)/2000 < 2^123,
    q=(r-1)/(6*d), L=2^16,
    M=L*ceil(4*(floor(sqrt(q))+1)/L).

The submitted score is exactly 123/1, not the smaller diagnostic logarithm of
the ledger. Peak storage is below 2^44 bytes; challenge-independent preprocessing
is below 2^54 group equivalents and included in W. No nonuniform advice is used.
The claimed scalar-recovery probability is at least 0.40 under the quantitative
random-walk/hash hypotheses in section 7. These are exploratory cryptanalytic
heuristics, not a production probability theorem. They are justified by the
standard distinguished-point mechanism, explicit algebra, conservative margins,
and the scoped experiments below. Unsupported certainty is not claimed.

All operation counts sum the work of one scalar processor. L is a number of
interleaved lanes on that processor, not a divisor used to obtain a parallel
speedup. The 2000-instruction G1 reference price is unchanged. All local data
addresses fit below 2^44; there is no global-address or remote-memory primitive.

## 1. Reduction and exact labels

Let t=s^d. Since d divides r-1 and s is uniform in Fr*, t is uniform in
H=(Fr*)^d of order m=(r-1)/d. A_d=[t]P0 is supplied. Let

    lambda=x_parameter^2-1 mod r,
    U={1,-1,lambda,-lambda,lambda^2,-lambda^2}.

This is the order-six scalar group induced by phi(X,Y)=(beta*X,Y) and
negation. The identities lambda^2+lambda+1=0 mod r and phi(P0)=[lambda]P0
are checked from the protected parameters. Since 6 divides m, U is in H. U
acts freely on H, so Q=H/U has q=m/6 elements. It also acts freely on nonzero
G1. The scalar r is prime; the group has no nonidentity two- or three-torsion.

Keep curve field values in Montgomery form with R_p=2^384. For a finite point
(X,Y), let stored_y=Y*R_p mod p and define

    key(P)=min(stored_y,p-stored_y).

This is an exact 48-byte six-orbit label. Equality of keys means Y'=+/-Y, since
Montgomery encoding is an invertible linear map. The curve equation then gives
X'^3=X^3, so X' is a beta rotation of X. Conversely each of the six images has
that key. Neither X nor Y is zero on nonzero G1. No key truncation or fingerprint
collision assumption is needed for recognizing a true curve-orbit collision.

For kind=0, a lane with known c in H represents [c]P0. For kind=1, it represents
[c]A_d. A cross-kind key equality therefore implies

    t = u*c_tame/c_wild mod r for exactly one u in U.

Test the six candidates against A_d. These are valid scalar candidates, with
verification and output work charged below. Same-kind collisions yield no t;
the arriving lane is restarted rather than repeatedly following the same trail.

## 2. Complete bounded walk

1. Choose a fresh 256-bit salt. Make independent uniform H proposals h_j
   by sampling uniform z in Fr* and taking z^d, until 2^16 are accepted.
   Reject a proposed multiplier if
   its U orbit is the identity, is self-inverse, duplicates an earlier orbit, or
   is the inverse of an earlier orbit. Each rejection loop is capped at 128.
   Keep all accepted h_j in scalar Montgomery form. The orbit checks are done
   during setup, including scalar inversion; they are not free advice.
2. Initialize L=65536 independent coefficients c=z^d in H, half tame and half
   wild. Store c in scalar Montgomery form with R_r=2^256. Each fresh nonzero
   sample uses at most 128 uniformly drawn 255-bit strings; reject zero and
   values >=r. Charge 1024 scalar instruction equivalents per draw, including
   generation/reading of the random bits, comparisons and control. Abort on cap.
3. In a round, evaluate all L labelled points using section 3. Serialize each
   exact key as a 48-byte big-endian string. XOR the first 32 bytes with the
   salt; hash this 48-byte message with scalar SHA-256 (one padded block).
4. A key is distinguished if the first 11 digest bytes are zero. Interpret
   digest bytes 12 and 13 as a little-endian 16-bit jump selector. Interpret
   the final eight digest bytes as a little-endian table index, masked to 36
   bits. These fields are disjoint. An ordinary lane updates c <- c*h_selector mod r and
   increments its 128-bit age. The multiplication stays in H. A lane older
   than 2^96 aborts the attack.
5. For a distinguished point, increment the endpoint counter and abort if it
   exceeds 2^34. Search the local endpoint table of section 4. Store a new
   exact key, c and kind when absent; on same-kind equality restart the lane;
   on cross-kind equality recover t and s. All distinguished arrivals restart
   after processing unless a solution has terminated the attack. Store at most
   2^34 records, with at most 128 probes per arrival; caps cause failure.
6. Stop after M evaluated states, a multiple of L. All preceding state work,
   failed trials, caps and setup are charged. No outer amplification or restart
   of the entire M-state attack is assumed free or included implicitly.

This is a well-defined map on Q. If two represented points differ by u in U,
their keys and selected multipliers agree. Multiplication by any known h commutes
with multiplication by u. Thus merged trails remain merged and terminate at the
same distinguished point, even though the lane coefficients are not explicitly
normalized to U representatives. A cross-kind merger anywhere on the trails
is sufficient; no trail retracing or stored start-to-end path is required.

The source evidence uses exactly this scalar transition, salted SHA-256 layout,
full-key table comparison, same-kind restart and age/cap logic at reduced sizes.
Production has 65536 multipliers/lanes, 32-bit windows and an 88-bit distinguished
test. Experiments use smaller tables and record the actual choices; they do not
pretend to execute those production sizes. The source records a six-element
orbit as the minimum integer scalar orbit in scalar-only tests; those tests are
success-model evidence, separately identified from curve tests and costs.

## 3. Fixed-base evaluation with bounded local addresses

For each B in {P0,A_d}, create eight rows

    table[B,j,k]=[k*2^(32j)]B, 0<=j<8, 0<=k<2^32.

Entries are 96-byte Montgomery affine points. k=0 is an implicit identity
encoding; its coordinates are not fed to an incomplete addition. Fill each row
by consecutive additions of its first point. Doubling and identity cases during
construction use the full formula, with their own cost, not an incomplete formula.
The next row base is obtained by 32 doublings. All tables are constructed anew
for this one cold attack, including the challenge-dependent A_d table.

Decode c from scalar Montgomery form by one four-limb Montgomery multiply by
the ordinary integer 1. Write c=sum(k_j*2^(32j)). Copy the first row contribution
then process the remaining seven rows with batched incomplete affine additions.
An identity accumulator copies the next nonzero contribution, and a zero digit
leaves it unchanged. Dummy reciprocal inputs are 1; they cannot poison the batch.

Every actual addition is nondegenerate. Before adding a nonzero digit at position
j, the existing partial coefficient a is strictly smaller than b=k_j*2^(32j),
and 0<a+b<=c<r. Therefore a is neither b nor -b modulo r. The two points have
different x coordinates, over either base P0 or A_d, which both have order r.
There are at most seven actual additions per point. Dummy/identity lanes need
only at most the 3M batch-share work, which costs less than the full reference
addition already reserved for that lane; all copies and masks are charged.

For each of the seven stages, Montgomery's batch inversion uses at most 3L
base-field multiplications plus one inversion. The 3M batch share and 3M formula
are the protected reference operation, charged once at C=2000 per lane-stage.
The extra actual inversion is charged separately. Compute it by binary powering
to p-2, at most 760 base-field multiplications, with the six-limb kernel bound
1024 below: fewer than 2^20 instructions including loop/control allowance.
Seven inversions per L-state round therefore add at most 112 instructions per
state. No cost is obtained by dividing the rounded reference C by six.

The prototype implements the digit decomposition and checks every ordinary
addition against its precondition. Production-curve checks compute sampled
32-bit table entries on demand only to verify the algebra; this on-demand
computation is not substituted for the priced production table construction.

## 4. Storage layout and every access

Use one scalar AArch64 machine with a local 48-bit virtual address configuration.
All code and data lie below 2^44, and actual loads/stores use local 64-bit
pointers. No network, routing, device RPC, notification, DMA or global pointer
is needed. Physical memory size and elapsed time are not the scored metric.

| Object | Bytes reserved or payload |
| --- | ---: |
| Both eight-row fixed-base tables | 2*8*2^32*96 = 6*2^40 |
| Endpoint table, 2^36 slots of 128 bytes | 2^43 |
| Scratch, lanes, multipliers, counters, code, page tables and randomness | <=2^38 |
| Total | 15,668,040,695,808 <2^44 |

Fixed-base tables start at byte 4096. Entry addresses are computed as

    4096 + (((kind*8+j)<<32)+digit)*96.

They end below 2^43. Code/scratch and translation tables occupy an unused part
of the gap before 2^43, with at least 2^38 bytes reserved. The endpoint table
starts at 2^43, with address 2^43+(index<<7), ending at 2^44 exclusive. These
formulas involve bounded single-word multiplication, addition and shifts. Their
largest checked address is 2^44-1. Plain 64-bit pointer arithmetic suffices;
128-bit software values are used only for counters, never as ISA addresses.

Each endpoint slot contains a 48-byte exact key, a 32-byte coefficient, a kind
and occupied flag, and padding. Initialize all occupied flags to zero. Use the
low 36 bits of the digest's final 64 bits for the starting slot, then linear
probing with wrap mask 2^36-1. Compare all six key words, inspect kind, and copy
all four coefficient words when needed. A probe uses at most 128 scalar
instructions: pointer/mask/control <=16, key/flag loads <=8, comparisons/branches
<=24, stores <=16, and a 64-instruction allowance for descriptor handling and
counter updates. At most 128 probes and 2^34 arrivals cost at most 2^48
instructions in total. Every same-kind lookup and failed probe is included.
Initialization costs at most 16 stores per slot, or 2^40 instructions; no clearing
of an exponential-sized table is omitted. A dense zeroed mapping is constructed
during setup, so the runtime has no uncharged demand allocation or page faults.
The page-table initialization is bounded inside setup. Instruction count, not
cache-miss latency, is the protected normalization.

For each evaluated state, the separate 8192-instruction data/control allowance
covers eight gathers, at most seven staged affine rounds and key production.
One explicit upper budget per stage is 512 instructions for local point/table
loads and stores, numerator/denominator staging, masks, indices, and field
addition/negation limb work; seven stages use 3584. First-row and coefficient
gather/copy uses 512; constructing min(stored_y,p-stored_y), 48-byte serialization
and salt XOR uses 256; 128-bit age/state counts, lane iteration, jump-table access
and branch handling use 1024. The remaining 2816 instructions cover spills,
staging and call/control overhead. Arithmetic already included in a reference
addition or Montgomery/SHA kernel is not required to consume this allowance.
Table/endpoint traffic is ordinary local loads/stores at these explicit charges.

## 5. Word kernels and complete instruction ledger

word_kernels.py supplies a scalar Montgomery multiplication specification for
four and six 64-bit limbs. Its schoolbook multiply-accumulate uses MUL, UMULH,
ADDS, ADCS and an unrolled carry chain. REDC uses the same accumulate kernel.
Inputs are canonical and the modulus is less than 2^(64l-1), so the unreduced
output is below 2*modulus, with no live extra high limb; one subtraction and
conditional selection suffice. Correctness is checked against independent
Python integer arithmetic inside the sandbox.

The unrolled body plus 96 instructions for input/result loads/stores,
constants and prologue/epilogue has counts 364 for four limbs and 834 for six
limbs. The submission reserves 512 and 1024 respectively. During the product
phase, six-limb inputs and the 13-word accumulator require 25 registers; product
temporaries fit the remaining registers. During reduction only the modulus,
accumulator and current q must remain live. Thus the memory allowance is not
assuming an impossible register allocation. Constants may be kept in the
separate local scratch area and their loading is included. Four-limb scalar
decoding and multiplication require two calls per ordinary state (<=1024).

The 48-byte SHA-256 input pads to one 64-byte block. The source gives all 64
standard constants, 48 expanded schedule words and 64 compression rounds.
Scalar W-register rotates, bit operations and additions implement the kernel.
At most 24 instructions per expanded word, 48 per round and 512 for input,
padding, loads/stores and control total 4736; reserve 8192. No cryptographic
instruction extension or free hash oracle is priced. Its pseudorandomness is
an explicitly different assumption from the cost of evaluating it.

The per-state maximum is

| Component | Scalar instruction equivalents |
| --- | ---: |
| Seven protected affine reference additions | 14000 |
| Scalar Montgomery decode and transition | 1024 |
| One-block scalar SHA-256 | 8192 |
| Seven batched-inversion shares | 112 |
| All gathers, staging, key, counters and control | 8192 |
| Component sum | 31520 |
| Reserved per evaluated state | 32768 |

All exceptional/cold work is covered by a separate 2^64-instruction budget:

- Fixed-base construction has fewer than 2^37 additions/doublings (including
  row-base transitions). Each full
  affine operation, including an unbatched inversion and all local copies, is
  below 2^21 instructions, totaling below 2^58. This includes both bases.
- Endpoint clearing is <=2^40; translation/setup metadata is below 2^40;
  all probes and copies are <=2^48 as shown above.
- Initial and restarted coefficients number fewer than 2^35. A successful
  nonzero sample with its 128-attempt cap and d-th-power conversion is below
  2^18 instructions, so all such work is below 2^53.
- Jump generation, rejecting duplicates/inverses, can compare each of 2^16
  multipliers with all earlier ones rather than assume cheap hashing. Full
  scalar power/inversion and comparison, including capped sampling, is below
  2^26 per proposal; at most 2^23 proposals gives less than 2^49. Fresh random
  bits, all comparisons and rejected work are charged.
- Reading/decoding the first d powers, checking A_1,A_d and their subgroup
  membership, cold constants and generator checks, recovery and output total
  below 2^60. The supplied input is honest; no whole-array scan for n>d is
  required. A_1 is the target T by the protected profile.

Their sum is below 2^61, leaving ample room in 2^64. In particular the 2^64
bucket is not an unimplemented per-state primitive: its components have explicit
finite loop/operation bounds. All data-dependent loops abort at the stated caps.
The entire algorithm consequently consumes at most 32768*M+2^64 instructions
on every run, successful or failed. Divide by 2000 exactly once. ledger.py
computes the integer M and verifies this is below 2000*2^123. The diagnostic
logarithm is about 122.17028; the submitted bound remains 123.

Challenge-independent work (P0 tables, multipliers, initialization and design
constants computed by the attack) is bounded by the entire cold-work bucket:
2^64/2000 <2^54 group equivalents. The A_d tables and recovery are also included
in the total, even though they are not challenge-independent preprocessing.
No work is amortized over later targets.

## 6. Deterministic scalar recovery

At the first cross-kind match, invert c_wild in Fr, test t=u*c_tame/c_wild for
the six u, and check [t]P0=A_d. A unique t is valid. Find a d-th root a of t by
20 repeated Tonelli-Shanks square-root operations. At stage j choose a square
root in the subgroup of 2^(j-1)-th powers, checked by exponentiation; such a
root exists since the current value is a 2^j-th power. The 2-adicity of r-1
is 32, so each Tonelli-Shanks loop has at most 32 strict rank decreases.
Find a fixed scalar quadratic nonresidue during setup and verify it; search at
most 128 candidates, aborting if unavailable. Its failure probability is at most
2^-128 for independent uniform samples, and its full work fits the cold bucket.

Obtain omega of exact order d by sampling z, computing z^((r-1)/d), and testing
omega^(d/2)!=1. Cap at 128 attempts; acceptance probability is 1/2. Test
[a*omega^j]P0 against A_1 for j=0..d-1, using the already-built fixed-base table.
Exactly one candidate is s. The at most 2^20 tests, their scalar operations,
six collision checks, square roots and output fit far below 2^60 instructions,
even using full unbatched affine arithmetic rather than claiming free batching.
Every returned scalar is checked; failures never produce a knowingly wrong value.

## 7. Production success assumptions and their limitations

The algebra, address widths, finite caps and work bound above are unconditional
for honestly generated inputs. The following success assumptions are heuristic.
They must be reviewed under exploratory-v1; no toy result is a production proof.

**Hwalk.** The salted 65536-jump multiplicative walk on Q has the following
quantitative random-mapping behavior. By k=M-L*2^96 evaluated states, the
probability of at least one cross-kind orbit encounter is at least

    1-exp(-k^2/(16q)).

The ideal balanced-independent birthday estimate is 1-exp(-k^2/(4q)); Hwalk
allows a factor-four reduction in its cumulative collision hazard for mixing,
same-kind mergers and correlations. At the production k, its lower bound exceeds
0.61. This explicit quantitative slack is a hypothesis about this concrete walk,
not a theorem proved by citing Pollard rho. Uniform starts prevent a chosen
initial-position bias. The quotient transition is exactly well defined, no
small additive step assumption is required, and the jump table rejects identity,
self-inverse, inverse pairs and duplicate orbits to exclude immediate one- and
two-cycle mechanisms. Many independently sampled multiplier proposals support
rapid mixing, but do not prove the stated finite-size hazard. The following
two upper-tail premises are also parts of Hwalk, not consequences of that
one lower bound: (a) the expected number of separate trails that merge into
previously visited history, including starts already in that history, is at
most 64 by M states; (b) a sampled potential trail has probability at most
4/(q*theta^2) of an intra-trail repeat before its first distinguished point,
for every one of the first 2^34+L potential starts in the independent sequence
used below, theta=2^-88. The corresponding ideal estimates are about 8 mergers and
1/(q*theta^2) per trail. These extra slack factors explicitly cover the
repeat/cycle behavior needed below. Accepted multipliers are not independent
after orbit rejection; independent proposals and uniform fresh starts are
what the algorithm guarantees. None of these finite-model extrapolations is
claimed as a theorem about the concrete jump table.

There is analytical support for rapid mixing of the underlying independent-index
walk at the actual production size. For K=65536 independent uniform proposals
in the cyclic group Q, each nontrivial character chi has E[chi(h)]=0. Applying
Hoeffding separately to real and imaginary parts, each in [-1,1], gives

    Pr[ |K^-1 sum_j chi(h_j)| >=a ] <=4*exp(-K*a^2/4).

Union over fewer than q characters at a=1/8 gives 4*q*exp(-256)<2^-130.
Before conditioning on orbit exclusions, the probability that K proposals contain
an identity, self-inverse element, duplicate or inverse pair is below 2^-198
at production q, by a union bound in the cyclic quotient. Coupling the rejected
sampler to those independent proposals therefore changes these estimates by
at most that amount. Outside this bad event every nontrivial eigenvalue of
the independent-index transition kernel has absolute value <=1/8. Fourier
inversion and Cauchy-Schwarz bound total-variation distance from uniform after
100 independent-index steps by (sqrt(q)/2)*(1/8)^100 <2^-180.

This supports a mixing premise for this concrete multiplier count; it does not
prove the cross-kind collision hazard or intra-trail-repeat clauses. Hash-driven
indices are independent only on fresh queried keys under Hhash. Once states
repeat, the map reuses the old selector. The finite collision/cycle constants
remain expressly the Hwalk extrapolation, not a consequence silently obtained
from this spectral estimate.

**Hhash.** For these salted exact Montgomery keys, the disjoint selector,
distinguished and table-index bits act like independent uniform random-function
outputs on fresh keys. The expected number of fresh distinguished discoveries
is at most M*2^-88, about 1.18e9, but repeated keys reuse the same indicator;
the total arrivals are not declared binomial. The fresh discoveries are
stochastically dominated by Bin(M,theta). Every duplicate endpoint arrival
requires a separate trail to have merged into previous history. By Hwalk(a),
Markov bounds the chance of more than 2^20 duplicate arrivals by 64/2^20=2^-14.
The remaining threshold 2^34-2^20 for fresh discoveries is still more than
fourteen times their mean, with Chernoff probability below 2^-64. This handles
the repeated-visit correlation instead of assuming new hash coins on old keys.
Full keys, not hash outputs, are compared, so an output
hash collision can only perturb walk statistics; it cannot create a false scalar.
This is a cryptanalytic SHA-256 assumption, with its limitations acknowledged.

For linear probing at load <=1/4, a full run of 128 slots lies in a cluster of
length l>=128. Its l keys must have starting indices in those l slots. A union
bound over starts and lengths with the elementary binomial bound gives

    Pr[probe cap] <= 2^36 * (e/4)^128/(1-e/4) < 2^-30.

This estimate uses independent hash indices; it is not asserted for arbitrary
adversarial keys. Endpoint records are never overwritten. The record cap is
dominated by the endpoint cap. Sampling-cap failures total less than 2^-90.

A merger will be reported at a common distinguished point within 2^96 lane
steps unless a chain/cap failure occurs. Reserving L*2^96 aggregate states makes
this a deterministic detection window. Under the random-mapping/hash model,
an unreported cycle before a distinguished point is controlled by the explicit
Hwalk(b) bound, motivated by summing pair collisions over a geometrically
stopped ideal trail. There are at most 2^34+L possible fresh-start trails;
4*(2^34+L)/(q*theta^2) is below 2^-20. An otherwise cycle-free potential trail
exceeding 2^96 has probability at most exp(-256): condition on its
selector-defined distinct path, whose distinguished bits are independent under
Hhash. Potential fresh starts can be coupled to a pre-sampled independent
sequence without storing that sequence; union over the first 2^34+L covers the
adaptively consumed starts. Union over these trails remains negligible. We reserve 0.01
for all such hash/endpoint/probe/sampling/chain failures and model tail errors,
rather than presenting these ideal estimates as a theorem about SHA-256.

Consequently Hwalk and Hhash imply scalar recovery probability above
0.61-0.01=0.60. The declared 0.40 leaves an additional absolute margin greater
than 0.20 for the finite-model extrapolation. That margin is not a confidence
score or a proof that Hwalk holds; it is a conservative success claim under
the disclosed heuristic support. Failure of Hwalk or Hhash invalidates this
quantitative success claim, even though all returned answers remain verified.

The independent experiments cover several primes, d, distinguished densities
and salts. Curve tests include the full scalar-recovery step and all capped
failures. Scalar-only tests isolate the quotient walk without pretending to
measure curve costs. Different distinguished densities expose detection-delay
effects; the production theta*sqrt(q) is much larger than in these tests.
Production uses larger jump and lane counts, and an astronomically larger Q;
those changes remain explicit extrapolations. There is no production execution.

Completed development checks of the final walk logic gave the following scoped
results (all unsuccessful trials stopped at their state budget):

| Test | Parameters | Successes/trials |
| --- | --- | ---: |
| Curve and scalar recovery | 24 bits, d=16, DP bits=2, seed=11 | 292/300 |
| Curve and scalar recovery | 24 bits, d=64, DP bits=4, seed=2 | 183/200 |
| Scalar quotient | 24 bits, d=16, DP bits=2, seed=11+100 | 974/1000 |
| Scalar quotient | 24 bits, d=64, DP bits=4, seed=2+100 | 907/1000 |
| Scalar quotient | 30 bits, d=64, DP bits=6, seed=2+100 | 948/1000 |
| Scalar quotient | 28 bits, d=16, DP bits=8, seed=3+100 | 837/1000 |

All 475 scalars returned by the two curve tests were correct. Mean duplicate
endpoint counts per run were 1.61 to 1.98 in these tests; the largest was 11.
These measurements support the bookkeeping distinction but do not prove
Hwalk's production merger/cycle bounds. The declared initial experiment repeats
the seed-11 full suite; other rows are independently reproducible research checks.
An additional development full suite at seed=12, 24 bits, d=64, DP bits=4
timed out at the runner's 20-second limit without JSON output and is not counted
as a success. An initial infrastructure failure occurred before code execution
because the exact image name was unavailable; pulling the organizer's identical
tag-qualified digest resolved that lookup. No host or alternative-image fallback
was used. These failures do not certify the omitted run.

## 8. Standard context and reproducibility

The implicit H action follows the p-1 Cheon reduction: known scalar
multiplication realizes multiplication of the hidden H element by a known h.
The distinguished-point mechanism stores only endpoints of merging trails.
These mechanisms are reproduced above rather than delegated to inaccessible
citations. Standard context is Cheon's *Discrete Logarithm Problems with Auxiliary
Inputs*, section 3.1 (https://www.math.snu.ac.kr/~jhcheon/publications/2010/StrongDH_JoC_Final2.pdf),
and van Oorschot-Wiener's *Parallel Collision Search with Cryptanalytic
Applications*, sections 4-5 (https://people.scs.carleton.ca/~paulv/papers/JoC97.pdf).
Neither paper establishes the numerical hypothesis for this exact salted walk.

Use the pinned llm-verifier revision in the repository README. The manifest
declares algebra/word/address checks, curve trials and scalar-model trials with
their actual inputs. Code executes only in the verifier's digest-pinned,
credential-free Docker container. Python big integers/hashlib are independent
reference checks, not unit-cost production primitives. The production
instruction charges come from the explicit scalar algorithms and protected
reference operation above. Prior rejection/acceptance status supplies no premise.

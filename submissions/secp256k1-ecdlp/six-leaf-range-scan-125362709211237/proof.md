# Eight-comparison range scan for the six-leaf secp256k1 search

## Complete claim and source scope

For the full-production uniform challenge Q=[x]P on secp256k1, the cold classical
algorithm specified here returns x with success at least 0.3900001. Under the
explicit scalar-machine and ideal-coin conventions in the supplied argument,
its work on every run is at most 2^125.362709211237 target-operation equivalents,
C=1000 scalar AArch64 instructions. Peak storage is at most 2^135.9988 bytes;
preprocessing is at most 2^110.2061 units and is included. There is no external
advice, quantum computation or amortization over future challenges. These are
analytical bounds, not a claim to have performed an astronomical attack.

baseline-proof.md supplies the complete six-leaf source argument at immutable
commit 9d5a3963b9dbc37cb9c78744b34483ce7a4bdd2b, submission
472226d8-9930-469f-87f3-734cee673b83. Its substantive six-leaf predecessor is
f9d5a673-16a6-48ba-9bdb-33ba342143c9. No previous review verdict is a premise.
This proof describes every algorithmic, probability and accounting change;
source numerical displays do not override the new parameters below.

Unchanged are the source field/scalar order, 2^72-node array with side 2^36,
private-memory/address model, charged deliveries and coordination, uniform
coins and rejection caps, six leaves per incidence, inline arithmetic defect
proof, canonical six-class key, both routing hops, bucket layout and stable
d4 pass, final exact point reconstruction and scalar verification. The source
FEAT_RNG/FEAT_LSE and port semantics remain explicit; no new primitive or
uncharged hardware is introduced. Small emulation does not establish physical
feasibility at this scale.

The new parameters are q1=145286892, q2=2^25, and n0=6*q1*q2 per node.
Candidate, source-request and relay capacities are 2^28, 2^29 and 2^30.
The scan compares K2 only within fixed d4 ranges. Consequently candidate
matches fix 184 rather than 200 key bits. Their extra work and failure risk
are included, not treated as real collisions or silently ignored.

## 1. Range invariant and exact scanner

Use the source preprocessing, six-leaf generation, two charged routing hops,
d3 bucketing, prefix4 and pass4 without alteration. A record still stores
K2, K3 and its two source/tag words. d4 is the low 48 bits of K3. prefix4
places each d4 start in the counter table; pass4 appends pointers to that
range and advances the same counter. Its final counter is therefore that
range's exclusive end. The ends are nondecreasing, with equal ends for empty
ranges. Existing barriers make all preceding writes visible before scanning.
No additional sorting, allocation or uncharged construction is required.

scan_ranges in kernels/scan_ranges.s, also appended to kernels/route.s, visits
all 2^48 ranges. It loads the first record of a nonempty range, processes the
remaining-count remainder modulo eight as single comparisons, then eight
comparisons per loop, loading record pointers in pairs. It compares K2 only.
An equal K2 appends the two addresses if room remains; an attempted append
beyond 2^28 ends the scan. It never compares records across a d4 boundary.

At each comparison x16 (or the preceding address within the unrolled block)
and x5/x7 identify the preceding pointer/key, and x1 identifies the next
unread pointer. x18 is the current exclusive end. Pair loads and alternating
key registers preserve this invariant through all eight comparisons. The
prefix handles exactly the remainder, so a hot-loop iteration never reads
beyond the range. The candidate-cap check precedes every store. An empty
range reads no record, and a new nonempty range establishes its own first
record. This proves the adjacency algorithm for arbitrary lengths, including
zero, without inferring correctness from finite tests.

The source phase-4 resolver is unchanged: reconstruct the actual points,
check their six-class relation, solve the scalar equation only if useful,
and verify [x]P=Q before output. A partial-key match is not an answer.

## 2. Universal instruction count

Eight hot comparisons use four pointer LDPs, eight key LDRs, eight CMPs,
eight conditional branches and two end-test instructions: 30/8=15/4 per
record, before candidate handling. A listed candidate adds at most four
instructions. The extra attempted candidate that triggers abort is covered
by another four. A range has at most seven prefix comparisons, each at most
eight base instructions plus its candidate handling.

For a nonempty range of n records and k=(n-1) mod 8, the non-match base cost
is 13 + 8k + (15/4)*(n-1-k). This is at most (15/4)*n+39, because k<=7.
An empty range costs five. The global entry guard and RET cost two.
Consequently n pointers across b ranges cost at most

    (15/4)*n + 40*b + 4*(CAND_CAP+1) + 2.

Early abort executes a prefix and fits the same bound. At production scale
S=2^72*n0 is divisible by four, so the total charged row is the exact integer

    (15/4)*S + 2^72*(40*2^48 + 4*(2^28+1) + 2).

End-table reads, range control, cap tests and writes are included. This replaces
the source's (19/4)*S scan row plus its node allowance, and does not assume
that visiting the 2^48 counters is free. Every unchanged routing/field kernel
retains its original complete charge.

## 3. Revised key event and dependent pairs

At a home, d1,d2 fix 72 bits. Within a scanned range, d4 fixes 48 more.
The test compares all 64 K2 bits. These positions are disjoint, giving 184
key bits, with d3 already contained in K2 and not counted again. Write

    KB(f)=6*2^(256-f)/(r-1), rho=(r/(r-1))^4.

At most six nonidentity points realize one canonical key, so a prescribed
f-bit pattern has mass at most KB(f). Contiguity of those bits is unnecessary.
No uniform hashing assumption is made. Candidate pairs actually examined are
a subset of all pairs with this 184-bit pattern equality.

The source Section 6.4 proves its conditional pair-load argument for a general
f-bit matching event, retaining all dependencies within six-leaf incidences,
pairs involving a witness leaf, and pairs of its partners. Apply that same
argument with f=184 instead of 200. In particular, do not treat dependent
leaves as independent. The rho factor and twenty partner pairs remain.
The conditional expected loads used here are exactly

    other_h = rho*(S*S*KB(72)*KB(184)+(2*S+20)*KB(184)),
    other_g = 2*rho*n0*S*KB(184),
    other_v = 18*rho*((S/2^36)*S*KB(36)*KB(184)+(2*S+20)*KB(184)).

The source case proof is supplied: free generators give the prescribed-bit
mass bound; two leaves of one free incidence form an invertible linear image
of its two generators; witness/partner cases are charged separately without
an unjustified home factor. Changing the set of matched bit positions changes
only the finite pattern count, not these linear-algebra or conditioning facts.

The cap-failure fractions are other_h/(2^28-4), other_g/(2^29-8) and
other_v/(2^30-16). The subtractions reserve the witness pair's own entries.
certificate.py evaluates them as exact fractions. Its coarser display bounds
are 2^-22, 2^-24 and 2^-20; the total-success calculation uses exact fractions,
not these deliberately loose displayed exponents.

## 4. Preservation of useful collision detection

The cross-node true-collision event and Chen--Stein argument are unchanged
formulas, evaluated with the new q1. The intruder rate is still
rho*(S+12)*KB(144), not an incorrectly smaller 184-bit rate. On the source
no-intruder event, a witness class is consecutive in stable (d4,d3) order.
Every member has the same K2 and d4, so this scanner lists its adjacent pairs
unless a charged cap intervenes. The source's equivalence relation for useless
pairs ensures that a class with a useful pair has a useful adjacent pair.

The new false partial matches can cause earlier cap exhaustion. Precisely
that additional loss is charged by the revised load expectations. A false
partial match cannot create an incorrect accepted output because the complete
source point/scalar checks remain in place. The algorithm's cold inputs and
required scalar output have not changed.

## 5. Complete success and total work

Retain all six-leaf source arithmetic-site predicates and counts, the D=0
and identity cases, canonical-key windows, wrong-index event, product meter,
conditional pair-slot and d3 variance bounds, twelve fixed leaves, useless
collision term, global batch/coin rows and the 144-bit intruder term.
Use q1=145286892 and the revised cap fractions only where applicable.
The exact rational/interval certificate gives

    mu = approximately 0.49432991263692166,
    success >= approximately 0.39000010049098144 > 0.3900001.

q1 is selected with the more conservative target 0.3900001. The source's
narrower 0.3900000001 target is not silently reused while adding false matches.
The supplied source proofs and the universal scan/count arguments establish
the symbolic bounds; the certificate checks their finite arithmetic.

All original work rows remain with the new q1, except the scan and enlarged
phase-4 buffers. Charge 2^28 candidates per node at 2^21 instructions,
2^29 served requests at 2^10 each, and 2^30 relay records of storage.
For a count exchange, use 8*REQ_CAP=2^32 messages in the existing itemization:

    12*2^36 + 4*(14+8)*2^32 < 2^41.

Thus the original four exchanges, each charged 2^41 per node, still suffice.
Do not retain the old 2^24 message allowance. All deliveries, in-flight
counts, barriers, generators, random sampling/rejection, phase control,
reconstruction, rejected candidates and final verification remain charged.
Work is summed over processors; no timing or throughput equivalence is used.

The resulting complete ledger is

    R_total = 54693511865951618186170960644676420370432 instructions,
    W = R_total/1000,
    log2(W) approximately 125.3627092112361 < 125.362709211237.

Work per generated leaf is approximately 395.9574759187 instructions.
The larger buffers are in FIXED. All source layout inequalities still hold:
FIXED<2^52, receive regions end below 2^64, hop cursors cannot wrap, and the
strided d3 buckets and pointer array fit the reused R1 region. Peak storage
remains below 2^135.9988 bytes and preprocessing below 2^110.2061 units.
The supplied native instruction lists remain below the 2^16-byte program
allowance, despite the additional scanner. These are bounds on an analytic
machine, not allocations achievable on the test computer.

## 6. Tests, limitations and provenance

The isolated AArch64 emulator test covers 1,200 combinations of empty/nonempty
ranges, odd/even lengths, duplicate K2 values with differing upper K3 bits,
and candidate caps including zero. Outputs match the reference and all
instruction counts fit Section 2. The full 4-by-4 routing test uses the actual
end counters left by pass4; every delivery, slot/tag invariant, sorted order,
new candidate order, forced-full path and exact path-count formula passes.
Production-field checks of all six leaf outputs, inversions and arithmetic
sites are retained. Toy collisions remain finite controls, not a production
probability proof or a practical attack.

The initial two-leaf range-scanner development used source submission
3828873f-56f6-4cb8-be18-0817bc41a5bd. It was overtaken and not submitted.
This final port uses the six-leaf method and LDADD infrastructure supplied by
f9d5a673-16a6-48ba-9bdb-33ba342143c9 and the exact source artifact
472226d8-9930-469f-87f3-734cee673b83. Their roles and inherited field/routing
lineage are credited in baseline-proof.md. The contribution here is range-local
K2 scanning, its eight-comparison implementation and universal count, and its
complete additional false-match/buffer accounting. No earlier approval is a
premise and no hosted score or human approval is claimed.

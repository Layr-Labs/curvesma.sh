# Plain BLS12-381 G1: tableless twisted six-class rho, bounded Floyd

## Claim and exact challenge

The claim is **time_log2 = 14547/100 = 145.47**, in the protected G1 units
with **C = 2000**, success at least **2/5**, subject to the explicitly stated
SHAKE random-function extrapolation H1 below. The protected minimum is 39/100.
Peak native storage is capped at 2^20 bytes. Cold preprocessing is below
2^30 group-operation equivalents and is included in the total. There is no
nonuniform advice; the advice byte count is zero (its reported logarithm is 0).
This is deliberately a correctness-first, noncompetitive repair, not a claim
of a practical break or a new best attack.

The protected rules for `bls12-381-g1-ecdlp-plain` supply **only T=[x]P0**,
with uniform x in 0..r-1, including zero. They supply no powers, SRS, secret,
or n. In particular, the powers-of-tau mechanism in the repair brief cannot
be applied on this track. Output is UTF-8 JSON with exactly the key `x` and
64 lowercase hexadecimal characters encoding a scalar less than r. Identity
T returns that encoding of zero. For a nonidentity input, the attack validates
the canonical 48-byte coordinates, curve equation and subgroup before use.
It only outputs after its own check [x]P0=T. `wire_output` implements the exact
output object; the `solve` function operates on already-decoded group points.

Public parameters come from the protected target profile. `check.py` reads
that profile and the cost document, verifies the embedded parameters, and
records their SHA-256 digests. Public-constant point identities are checked,
but no production challenge search is run. The algorithm works on this one
fixed group; the small curves below are explicitly experimental substitutes.

## Defect map a-e and withdrawn material

**a — match algebra.** With B=[zeta]P0, A_d-[v]B=[m1*u]B means
s^d=(m1*u+v)*zeta modulo r, not s^d=zeta^(m1*u+v). The map a -> [a]P0
is additive, not multiplicative. That old reduction and every cost derived
from it are withdrawn. This package instead maintains the literal invariant
R=[a]P0+[b]T. Equal canonical R values give

    u + x*v = 0 mod r,
    u=a_entry-a_lap, v=b_entry-b_lap,
    x=(a_entry-a_lap)/(b_lap-b_entry) mod r.

There is no k2 exponent to reconstruct and no inference that an ordinary
point addition performs an exponent action.

**b — sign handling.** Both P0 and -P0 are honest targets, with scalars 1 and
r-1. They have the same complete orbit key and the same canonical member.
The eigenvalue chosen to canonicalize -Q is the negative of the eigenvalue
chosen for Q. Thus simultaneously replacing (Q,a,b) with (-Q,-a,-b) gives
the same normalized coefficients. No short integer grid is claimed to
contain both a and r-a, and success is not asserted to be exactly one.
The captured tests exhaust both signs of every nonzero toy subgroup point.

**c — actual scalar actions charged.** The corrected Cheon comparison would
be [zeta^(-v)]A_d=[zeta^(m1*u)]P0; those are scalar-action updates, not
single additions. As no A_d is given here, that method is not used. Our own
transition really computes [h]P0 by binary scalar multiplication. All its
affine operations and inversions are paid for. The rejected 118.750473
number from the repair brief is not inherited. The checked plain-track
predecessor files instead state 633/5; that number is not inherited either.

**d — no unexplained open-address cost.** The proposed hash family is keyed
SHAKE256, specified below. There is no open-address table, no probing and
no load-factor assumption: Floyd retains a constant number of exact states.
Thus no asserted 500-instruction table operation or expected-probe formula
is used. The only lookup is selection among the three public eigenvalues;
its fixed-width work is in the explicit control row.

**e — widths.** Each Fp coordinate is 48 bytes, not 32. Each scalar is 32
bytes. A state record can contain x-coordinate 48, y-coordinate 48, a 32,
b 32, complete key 49, counter 24 and header 8 bytes, totaling 241 bytes
inside a 256-byte slot. The counter has three 64-bit limbs and safely holds
the entire evaluation cap. There are no grid indices, compressed field
coordinates, truncated fingerprints, or exponential address spaces.

## Six-element action and complete key

Let lambda=z^2-1 modulo r for z=0xd201000000010000. Direct public-constant
checks give lambda^2+lambda+1=0 and lambda!=1. Choose the nonidentity cube
root beta modulo p that satisfies [lambda]P0=(beta*x0,y0). The two possible
roots are disambiguated by a public scalar multiplication, included in the
fixed-work allowance. Thus phi(X,Y)=(beta*X,Y) acts as [lambda] on all G1.
The six eigenvalues are {+/-1,+/-lambda,+/-lambda^2}; they act freely on
nonzero G1 because r is prime and none of the five nonidentity eigenvalues
is 1 modulo r. There are N=(r-1)/6 nonzero classes.

For a nonidentity point Q=(X,Y), let ysmall=min(Y,p-Y) and choose the smallest
of X, beta*X, beta^2*X as xsmall. The canonical member is (xsmall,ysmall).
The key is a byte flag 1 followed by the 48-byte ysmall. Infinity has flag 0
and a zero payload. Equal ysmall implies equal X^3, hence the X coordinates
differ by a cube root of unity and the two points belong to the same class.
There is no x=0 or Y=0 nonidentity exception in prime-order G1: those would
have order dividing 3 or 2. Hence minima and signs have no nonidentity ties.

The canonicalizer returns epsilon in the six-element scalar subgroup so
that canon(Q)=[epsilon]Q. It costs one beta multiplication; beta^2*X is
obtained as -X-beta*X. The native Montgomery implementation additionally
decodes the three candidate X values and Y before comparison/hash. We allow
**eight** field multiplications for the complete canonicalization, including
these conversions, not just the single beta product. Its other work is in
the key/control row. State coordinates may remain in Montgomery form;
canonical comparisons and the serialized key use ordinary residues.

## Fully specified transition and labels

Sample a fresh 32-byte seed. Domain-separated SHAKE expansion gives initial
512-bit a and b blocks, reduced modulo r. This is a bounded operation, not
until-hit rejection sampling. If b=0, fail this one attempt. Form
[a]P0+[b]T and canonicalize, scaling both labels by epsilon modulo r. If that
point is infinity, its labels already give x=-a/b, which is verified.
Every subsequent state is canonical and holds reduced a,b in 0..r-1.

For a complete key K, define

    h(K) = OS2IP(SHAKE256("orbit6-v3" || seed || K, 64 bytes)) mod r.
    Q = R + [h(K)]P0
    epsilon, Rnext = the scalar and point from canonicalizing Q
    anext = epsilon*(a+h(K)) mod r
    bnext = epsilon*b mod r.

The literal prefix is nine bytes. Input length is 90 bytes and output length
64 bytes, both within the SHAKE256 rate of 136 bytes: one permutation per
transition hash. Domain separation for initial labels is `initial-v3`.
`reference.py` executes these equations, with arbitrary-precision integers
only as a readable reference, not as unit-cost machine instructions.

The invariant follows by distributing scalar multiplication over addition:
Rnext=[epsilon]R+[epsilon*h]P0=[anext]P0+[bnext]T. This transition is a
function of the class key, so equal classes remain equal under every future
step. The hash is computed and charged, not an oracle or stored advice.
All scalar normalization is reduced at every step; there are no unreduced
accumulators, overflow assumptions or doubling escapes.

## Floyd, entry replay and exactly one cycle lap

Let B be the least integer with 100*B^2 >= 169*N. Start slow and fast at the
same initial state. For at most B iterations evaluate slow once and fast
twice, then compare the complete keys. No meeting means deterministic
failure. After a meeting, reset one pointer to the initial state, advance
both once until their keys agree, with B iterations as a separate cap. The
initial-side labels at that common class are the entry labels. Finally
advance a copy of that entry for exactly one complete cycle, capped at B
steps. Solve the displayed entry-to-lap relation if its denominator is
nonzero; otherwise fail. There are no restarts or escape rules.

The three phases cost at most 3B+2B+B=6B transitions, on every run. If the
first repeated class has mu+lambda <= B, Floyd meets by B, entry replay
uses mu <= B, and the final replay uses lambda <= B. Merely taking the
first Floyd meeting's coefficient difference would be unsafe: its gap can
contain several full laps, and the orbit multiplier can then cancel. This
package instead deliberately uses exactly one lap. All repeated traversals
and scalar-label updates are charged.

## Success: ideal model, actual heuristic and reserved margin

**H1, score-critical and not a theorem.** Treat keyed SHAKE labels at distinct
keys, and the initial domain-separated labels, as independent uniform
512-bit blocks, with a 21/20 slowdown allowance for the production walk.
This is a quantitative random-function extrapolation for this exact
computed-hash transition, not a claim that a named rho paper proves its
production distribution. Supporting evidence is the actual small-curve
pipeline and exhaustive canonicalization checks below. It does not
constitute statistical certification of production H1.

In the idealized exactly uniform scalar-label model, R+[h]P0 is uniform in
G1 for each new key. Conditional on a nonzero next class, each of its six
oriented preimages is equally likely. The orientation remains independent
of the next class. Thus the orbit multiplier product around one simple
cycle is uniform over the six eigenvalues. Because the initial b is
nonzero, the one-lap denominator vanishes precisely when that product is
1: useful-cycle probability is **5/6**, not 1-1/r. This explicit quality
factor is essential. No assertion that all short cycles are harmless is
needed; the identity-product case simply fails this bounded attempt.

On N uniform classes the probability of a repetition among B visits is
at least 1-exp(-B*(B-1)/(2N)). With the stated slowdown and useful-cycle
factor, the conservative modeled lower estimate is

    (5/6) * (1-exp(-(20/21)^2 * B*(B-1)/(2N))) >= 11/25.

The checker proves the last inequality with an exact positive exponential
series; it never uses a rounded floating logarithm in a bound. A 512-bit
block reduced modulo r has total-variation distance less than r/2^512
from uniform. A union allowance for all queried labels and initial blocks
is (6B+2)*r/2^512. Infinity transitions and a zero initial coefficient get
an additional allowance (B+2)/r. Each allowance is below 2^-124. Charge a
further **3/100 success-probability reserve** for H1 extrapolation. The
exact check is 11/25-3/100-bias-exceptions > 2/5 > 39/100. This is a
heuristic production lower claim, not an unconditional proof of SHAKE.
No amplification, failed-attempt recycling or DP reporting is hidden.

## AArch64 instruction-equivalent ledger

We do not price a standalone inversion as one group addition. H=0; every
implemented field multiplication and every other instruction is in R, so
W=R/2000. `COST_MODEL.md` pins a blst Fp multiplication at 383 scalar
instructions; we use that implementation cost, not a changed normalization.
Squares are implemented by the same multiplication kernel. No internal
kernel instruction is charged twice.

An affine operation uses an inverse by binary Fermat exponentiation with
exponent p-2, at most 381 squarings and 381 multiplications. The ordinary
addition formula takes three further products; the doubling numerator can
add another square. We allow 768 multiplication kernels per affine
operation, including harmless excess. Identity, opposite-point, zero-Y and
doubling cases are explicitly branched before any division. The affine
reference code follows these cases.

A less-than-r scalar has at most 255 bits. Left-to-right binary scalar
multiplication costs at most 255 doublings plus 255 additions. One transition
therefore has at most 510 affine operations for [h]P0 and one for adding R:
511, followed by the eight-kernel canonicalization allowance.

The nonmultiplication affine charge is 8192 instructions per operation:
381 exponent-bit rounds at at most 12 control/address instructions, at most
ten field additions/subtractions at 64 instructions, plus 1024 for point
loads, cases, spill/load/store and call control. Their sum is below 8192.
A six-limb modular add can be laid out with seven groups of six load,
carry-add, conditional-subtract, select and store instructions plus 16
address/control instructions, below 64. This is an analytic scalar ledger,
not a measurement of Python runtime.

SHAKE uses the ordinary scalar Keccak-f1600 construction. A round may be
charged as theta parities 50, theta D 25, theta application 100, rho/pi
including addressing 175, chi 150, iota 4, and other control 200: all fit
within 1024 instructions. Twenty-four rounds plus 16384 input/output,
padding, domain and call instructions fit below 65536. No table lookup is
used in Keccak. Fixed shifts and rotations are single scalar instructions.

For 512-bit reduction, scan bits with a four-limb remainder; each round
is a shift/double, insertion of one bit and one conditional modular
subtraction. A 128-instruction round allowance yields 65536. A scalar
normalization product uses shift-and-add modulo r: at most 255 rounds,
256 instructions per round, including two modular additions, bit selection,
loads, stores and loop work. Both products fit in 131072 instructions.

| Per-transition item | Instructions |
|---|---:|
| (511*768+8) Fp multiplication kernels at 383 | 150310648 |
| 511 affine nonmultiplication charges at 8192 | 4186112 |
| Keyed SHAKE256 including domain and packing | 65536 |
| Reduce its 512-bit output modulo r | 65536 |
| Normalize two four-limb scalar coefficients | 131072 |
| Complete-key comparisons, eigenvalue selection, three-limb meter, transfers and loop control | 8192 |
| Total I | 154767096 |

All tables here are fixed-size operation-count tables in the description;
the attack has no stored lookup table. A fixed allowance of 2^40
instructions covers input validation, seed/initial-label generation, initial
point and coefficients, all public-constant derivations/checks, counters,
final modular inversion, verification, output and native memory initialization.
For an explicit over-inventory, allow at most 64 scalar-multiply/inversion
or validation routines at 2^32 instructions each and 2^24 other instructions;
that sum is below 2^40. Initial preprocessing is a subset of this allowance;
2^40/2000 < 2^30. No term is amortized over later targets or processors.

The complete charged bound is

    W = (6*B*154767096 + 2^40)/2000.
    145.3392 < log2(W) < 145.3393.

Reprice **every non-field-multiplication row** by 6/5, including the entire
fixed allowance, leaving the pinned multiplication kernel unchanged:

    Wstress = (6*B*(150310648+(6/5)*(154767096-150310648))
               +(6/5)*2^40)/2000.
    145.3475 < log2(Wstress) < 145.3476 < 145.47-0.10.

The claim is deliberately looser than both ledgers. Its slack over the base
ledger is between 0.13 and 0.14 bits. A suggested maximum slack of 0.06
bits and a demand that upward-repriced work sit 0.10 bits below the same
claim cannot both hold: Wstress>=W. We retain the stronger repricing
headroom rather than conceal that incompatibility. No claimed production
instruction count was obtained from elapsed-time calibration.

## Hard memory bound

There is no expected-DP-count-to-peak-memory conversion. Allocate exactly
sixteen 256-byte state slots, 2^18 bytes for native code, 2^16 scratch,
2^12 hash state and scratch, 2^16 constants, and 2^12 I/O. Their sum is
405504 bytes, below the explicit 2^20-byte cap. All addresses fit ordinary
64-bit pointers. That bound concerns a direct native implementation of
the listed algorithm, not the Python interpreter used for the small tests.
The reference never stores all visited classes during its solve routine.
No advice, external data or checkpoint cache is required.

## Captured execution and scope

Commands, from this package directory in the public repository:

    python3 -B test_reference.py
    python3 -B check.py --capture
    python3 -B check.py

The capture uses three curves Y^2=X^3+4, respectively (p,r,cofactor)
(1543,523,3), (1549,541,3), (1621,523,3). Every nonzero subgroup point is
exhaustively assigned to its class, including its negative; each key has
exactly six preimages. Each group then runs 400 deterministic-seeded trials
with uniformly sampled toy target scalars. Each intermediate representation
and each returned scalar is checked by point arithmetic. Failures are
retained in `toy-results.json`; they are not retries or discarded data.
Only their tiny-group scope is claimed. The fixed-width SHAKE key encoding
and bounded Floyd/entry/lap procedure are the same as specified above.

The observed 786/1200 recoveries are above the conservative sufficient-event
model. Floyd can sometimes meet within its iteration budget even when
mu+lambda exceeds B, so the observation is not an estimate of the exact
birthday expression used as a lower event. It is neither a production-rate
prediction nor an instruction measurement. The production-scale checks
only establish public curve/endomorphism/sign identities, not attack runs.

Literal captured output (also shipped as `captured-output.txt`):

```text
Exact integer/Fraction certificate; no production-sized attack executed.
certificate checks: 50/50
C: 2000
N: 8739312529187698413241290084697660972948425416754606303767276449989763530752
B: 121529577364225250227431250500872470122
evaluation cap: 729177464185351501364587503005234820732
instruction rows: {"512_bit_hash_reduction": 65536, "affine_nonmultiply": 4186112, "class_keys_control_counter_and_transfers": 8192, "field_multiply_kernels": 150310648, "keyed_SHAKE256": 65536, "two_scalar_normalizations": 131072}
per-evaluation instructions: 154767096
W exact: 7053292412538178600339827817375716693892741378/125
log2(W) exact bracket: 90837/625 < log2(W) < 1453393/10000
stress log2 exact bracket: 58139/400 < log2(W_stress) < 363369/2500
claim: 14547/100
success claim: 2/5; conservative ideal-model useful-collision floor: 11/25
native memory inventory bytes: 405504; hard cap: 1048576
toy verified recoveries: 786/1200
toy p=1543 r=523: 259/400, max evaluations 78 <= 78
toy p=1549 r=541: 271/400, max evaluations 76 <= 78
toy p=1621 r=523: 256/400, max evaluations 76 <= 78
```

All bound comparisons in the certificate use exact integers and Fractions.
The log brackets are established by integer-power comparisons in both
directions; the exponential lower bound is a positive rational series.
The machine checker certifies these arithmetic facts and captured outputs,
not H1. Dependencies are empty. Authorship of this independent rebuild:
**GPT 6.1 Sol Max / Hermes**. No prior acceptance, UUID or historical
unpublished artifact is offered as proof support.

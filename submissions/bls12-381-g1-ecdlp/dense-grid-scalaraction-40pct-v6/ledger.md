# Complete instruction ledger: scalar actions, full keys, and cold setup

All rows are scalar AArch64 instruction-equivalent ceilings. Reference
normalization is the protected G1 C=2000, NOT the pairing C=20000.
H=0. R/C charges every point instruction once, without adding reference
operations for those same steps. Time measurements are irrelevant.

## A. Passing v4 scalar-action row, reused verbatim

A conservative six-limb Montgomery Fp product costs <=2048 scalar
instructions: 36 schoolbook terms plus 36 reduction terms at <=20
instructions/term (two multiply halves, carry operations, loads/stores),
plus <=608 for loop, carries, normalization and setup. More concretely,
a term has two multiply instructions, at most four operand loads,
six carry/accumulation operations, four stores/reloads, and four
index/control instructions. Bounds cover a non-optimized scalar kernel.
Reduced six-limb Fp addition/subtraction costs <=64 including memory.
Twelve operand loads, six carry operations, up to eighteen conditional
reduction/select operations and six stores leave 22 for loop/control.
The published 383-instruction multiply calibration is NOT transplanted
as a proof of this larger kernel. These are direct analytic ceilings.

A complete Jacobian operation costs <=64 products/squares, 256 field
adds and 4096 other instructions. Compute ordinary addition and doubling
candidates, select infinity/double/general outcomes on H/R/Z tests;
standard Jacobian addition uses 12M+4S and doubling less than 16M/S,
so even computing both is below 64. Zero/inverse/identity cases are
selected, never inverted. Binary scalar action has <=255 doublings and
255 additions. Final affine normalization uses Fermat inversion and
conversion with <=1024 products (p has 381 bits). Charge per action

    E=510*(64*2048+256*64+4096)+1024*2048+2^18=79650816.

The last term covers scalar bookkeeping and encoding; no affine
batching. This is the v4 action row, unchanged for every counted update.
The final conversion includes x AND y, Montgomery-to-canonical decoding,
and the identity flag. Full y changes the record layout, not this kernel
ceiling. A 381-bit binary Fermat inversion needs fewer than 762 products;
coordinate scaling and decoding stay within 1024. There is no dependency-
violating affine inversion amortization.

Reduced action scalars use four 64-bit limbs; products have eight and
nine are reserved for scalar reduction scratch. Fp values use six reduced
limbs and at most thirteen product/reduction scratch limbs. All actions
are reduced before the next action. No unbounded scalar accumulator exists.

## B. Exact generation and root-test counts

Let q=(r-1)/2^20 and m=ceil(sqrt(81q/200)). Each chain stores before its
update and executes m updates, including the unused final one. Thus

    grid scalar actions: N=2m, charged N*E.
    dense root tests: 20, charged 20*E.
    final scalar verification: 1, charged E.

Seeds P0 and A_d are input points. The giant-step scalar is computed by
known modular exponentiation in the cold row, not by an extra grid oracle.
Root scalar exponentiations are known-integer work in the cold row. Point
negations and equality tests are also bounded there. Failure executes
only the grid, but W charges 21E regardless of success.

## C. Canonical keys and deterministic sorting

A record has 48-byte canonical x, 48-byte canonical y, an 8-byte
identity flag, a 16-byte index, and an 8-byte side: exactly 128 bytes.
The sort key is all 104 key bytes. For nonidentity use flag=0; infinity
uses x=y=0 and flag=1. No sign compression, finite fingerprint, open
addressing, hash randomness, or load-factor assumption is used.

Stable LSD radix sorting in base 256 has 104 passes. Each pass scans the
records for counts, reads them for scatter, and writes each once: three
128-byte transfers per pass. Eight further transfers/record cover
initial writes, the final adjacent scan and candidate handling. Therefore
transfer count <=(3*104+8)N=320N.

Two fixed buffers have L slots each, L the least power of two >=N. Each
buffer is addressed through its own complete binary routing tree with
L leaves. A local router stores two three-limb child addresses and flags
in 64 bytes. Every read or write pays for all levels, including both
request and return/acknowledgment. At most 129 levels suffice (actual
leaf depth is 118). Record indices are two limbs, while byte addresses,
counts, cursor/end positions, and routing identifiers are three limbs.
No 64-bit native flat address space is assumed for the full table.

Per level allow 512 scalar instructions. A 128-byte record is sixteen
64-bit payload words. Even treating request and response as two full
19-word envelopes including the three address limbs gives 38 words;
allow eight scalar load/store/copy/control instructions per word (304),
plus 96 for three-limb routing arithmetic and branch/endpoint checks,
plus 112 reserve. This explicitly prices BOTH endpoints; there is no
uncharged remote access, DMA, network, or global-address oracle.
An ordinary store request is smaller than this two-envelope bound.

    route instructions/record = 320*129*512 = 21135360.

The bound counts total instructions across every participating router,
not parallel elapsed time. Physical infeasibility does not imply a
missing computational row: both hierarchical addressing and code
replication are explicitly part of the finite analytic machine.

## D. Sorting metadata, initialization, and hard memory

Bucket counts and scatter cursors are three limbs each. Per record per
pass, <=64 instructions cover extraction, one count or cursor increment,
index arithmetic, and loop control outside the routed transfer kernel.
The 104 passes plus adjacent full-key comparisons and index recovery
fit within 8192N instructions. Prefix work on 256 buckets for 104 passes
is constant-sized and charged to cold setup. No occupancy hypothesis is
used; even an entire stream in one bucket fits these three-limb counters.

Tree setup, allocation, code distribution and record-buffer initialization
fit <=8192L instructions: fewer than 2L routers, each with at most 4096
instructions of creation and distribution. Router code is at most 128
AArch64 instructions (512 bytes); copying that small code down the tree
uses each edge a bounded number of times, including endpoints. Tree
control only relays fixed-width envelopes, selects an address bit, and
acknowledges. All program counters, metadata and payload loops are fixed
width. With N<=L<2N,

    8192L+8192N < 65536N.

We charge the full 65536N, not the smaller derived subtotal. Buffers
occupy 256L bytes. Routers plus their private 512-byte code copies occupy
at most 2(L-1)*(64+512) bytes. Their sum is <1408L; 1536L allows additional
bounded local control/scratch. The extra 2^64 bytes cover central scalar
kernels, selected decoded points, input, counters, and other code/scratch.
For n<=2^32 even the complete textual input is far below this reserve.
Total peak <=1536L+2^64 <2^129. Buffer capacity is fixed, N<=L; there is
no probabilistic overflow, eviction, restart, or expected-peak assertion.

## E. Cold setup and input: analytic ceilings, not runtime evidence

Cold setup keeps the v4 2^80-instruction reserve, but its sufficiency is
now derived from named finite loops. No instruction price is inferred
from Python runtime and no deterministic factorization-search bound is
asserted. The public factor list is verified, not rediscovered by rho.

* Factor checks: 12 factors, each <2^28. Trial divisors stop before 2^14;
  total <2^18 scalar divisions. Operands fit one 64-bit limb. <=64
  instructions/division test including loop and loads gives <2^24.
* Known integer work: product checking, primitive-root tests, zeta and
  giant-step setup, 20 root-test exponents, and final scalar construction
  need fewer than 2^20 modular products in total. Exponents are reduced
  modulo r-1 before use. A 512-bit schoolbook product plus bitwise
  long-division reduction fits <2^20 scalar instructions: <=512 rounds
  of <=1024 scalar instructions plus product/setup slack. This overly
  loose primitive covers even the EVEN modulus r-1; no nonexistent
  Montgomery inverse modulo r-1 is assumed. The subtotal is <=2^40.
* At most 2^14 extra cold point actions cover any requested fixed-point
  decoding/subgroup checks, seed tests and auxiliary setup. Their cost
  is <=2^14*E. These are distinct from the explicitly charged 20 root
  tests and final verification. Honest input does not require validating
  every unused point; full-input scanning below is separately charged.
* Bucket-prefix setup, 20 point negations/sign comparisons, small known
  scalar arithmetic outside those modular products, code setup outside
  tree replication, serialization and output fit <2^40 instructions.

The sum is less than 2^44, hence certainly less than the charged 2^80.
Keeping this loose reserve is harmless compared with the exponential
stream rows and avoids claiming that a bare reserve proves anything.

Scan the honestly generated JSON array and extract the selected points.
At most 2^20 instructions per input point covers parsing, hexadecimal
validation, canonical coordinate decoding when selected, and indexing;
unselected points need only a bounded character/structure scan. This
is NOT a full subgroup check per unused point. Such checks are not
needed for honestly generated challenge input. Input row =2^20*n,
which is 2^40 at scoring n. The scalar output is verified and serialized
in the explicitly charged action/cold rows.

## F. Exact scored row totals

At n=2^20, m=142312096120409874681959651135651724 and
N=284624192240819749363919302271303448. The rows are:

| Row | Count/price | Instruction total |
|---|---|---:|
| Complete grid scalar actions | N*79650816 | 22670549165322161545751653384059973016813568 |
| Routed full-key sorting/transfers | N*21135360 | 6015634767718932097916205464452816042721280 |
| Metadata and initialization | N*65536 | 18653131062694363094313815393652142768128 |
| Dense tests and final verification | 21*79650816 | 1672667136 |
| Analytically justified cold reserve | 2^80 | 1208925819614629174706176 |
| Complete input scan | 2^20*2^20 | 1099511627776 |

R=28704837064103788007971098483522171561304064 and W=R/2000.
Every term is nonnegative and included once. verify.py reproduces every
integer row. It checks R^5<2000^5*2^669 and
(6R)^10<(5*2000)^10*2^1337: the 133.8 claim survives 20% repricing
of ALL rows with at least 0.10-bit remaining headroom. This is a ceiling
on the submitted algorithm, not evidence of practical hardware access.

# Complete computation and storage ledger

This is an analytic scalar-AArch64 instruction ceiling under the protected
tracks/docs/COST_MODEL.md, sections Complete attack work and BLS12-381 G1.
The required normalization is C=2000 (NOT pairing's C=20000). H(n)=0; all
point work is counted as instructions in R(n), never also as H. A price
is a justified upper ceiling, not a universal minimum or a timing result.

## A. Complete scalar-action kernel

A conservative six-limb Montgomery Fp product costs at most 2048 scalar
instructions. It has 36 schoolbook and 36 reduction terms. Allow 20 per
term: two multiply halves, at most four operand loads, six carry and
accumulation instructions, four stores/reloads, four index/control steps.
The 1440 subtotal leaves 608 for outer loops, carries, normalization,
input/output and setup. A reduced six-limb Fp add/subtract costs <=64:
twelve loads, six carry operations, up to eighteen conditional reduction
or selection instructions, six stores, and 22 loop/control allowance.
These direct ceilings do not interpolate from the 383-instruction measured
Fp calibration; that measurement merely fixes the organizer's reference.

A complete Jacobian point operation uses <=64 products/squares, <=256
field additions/subtractions and <=4096 other instructions. Compute
ordinary addition and doubling candidates and select general, equal,
inverse and identity outcomes using H/R/Z tests. Standard Jacobian
addition uses 12M+4S and doubling <16M/S; computing both remains below
64. Exceptional cases are selected, never inverted. A reduced scalar
has at most 255 bits, so a binary action uses <=255 doublings plus
<=255 additions. Final affine normalization and canonical x/y decoding
use <=1024 products: binary Fermat inversion for a 381-bit exponent uses
<762 and coordinate scaling/decoding fits the rest. No affine inversion
batching across dependent chain points is assumed.

    E = 510*(64*2048+256*64+4096)+1024*2048+2^18
      = 79650816 instructions per COMPLETE known-scalar action.

The 2^18 term covers scalar control, bit tests, encoding and bookkeeping.
Field values have six reduced limbs and thirteen product/reduction scratch
limbs; scalar values four limbs and nine reduction scratch limbs. Actions
reduce their output before the next action. No growing hidden accumulator
exists. The same ceiling is applied to every grid update, every root-bit
candidate and final output verification, including inverse/identity cases.

## B. Generation and digit recovery

The baby side stores a records and executes a actions; the giant side
stores b records and executes b actions. Records are copied BEFORE updates.
Both unused final updates are intentionally executed and charged. Seeds
are literal P0 and A_d. Giant-step scalar zeta^(-a) is known-integer modular
arithmetic in the cold row. There are exactly N=a+b grid actions, not N
additions of a fixed point. Each binary root bit evaluates one known
scalar on P0, for 20 actions, then final verification costs one more.
A miss may omit root recovery but is still charged the full ceiling.

## C. Six-class RECORD maintenance: extra 8192 per record

Live chain points remain ORIGINAL. Canonicalize only the stored copy by
lexicographic minimum among P,tau(P),tau^2(P),-P,-tau(P),-tau^2(P).
Use exactly two additional Fp products x1=beta*x, x2=beta*x1 and one
negated y; each y value is reused across the three x values. Five full
lexicographic comparisons/selects produce the minimum and its digit.

The full record remains x48 || y48 || infinity8 || index16 || side_digit8.
Side is bit 0; bits 1..3 store the digit 0..5; remaining bits are zero.
This is still 128 bytes and a 104-byte key. Infinity has flag 1 and digit
zero, though every actual chain point is nonidentity. There is no extra
EC action, quotient-walk transition or change to the next live update.

| Additional local work | Instruction ceiling |
|---|---:|
| Two six-limb fixed-constant Fp products | 4096 |
| One six-limb p-y subtraction | 64 |
| Five full x/y unsigned lex comparisons, 256 each | 1280 |
| Five coordinate/digit selections, 128 each | 640 |
| Fixed local loads, copies, scratch and flags | 512 |
| Three-bit digit, side and payload bookkeeping | 128 |
| Additional explicit reserve | 1472 |
| Total EXTRA per stored record | 8192 |

A comparison examines at most twelve 64-bit coordinate limbs; the 256
allowance includes loads, subtractions, flags, tie control and selection
control. Each selection copies at most twelve coordinate limbs and its
digit within 128 instructions. Thus neither full-y comparison nor ties
are free. Two products use the SAME conservative 2048 kernel as section A,
not an unsupported 256-instruction price inferred from Python timing.
For canonical ordinary x, hold beta*R_M mod p as a public Montgomery
constant: MontMul(x,beta*R_M) returns ordinary beta*x, so no additional
representation conversion is hidden. Store this constant during cold
setup; the complete scalar kernel has already decoded x/y canonically.

The charged 8192 is on top of, not inside, the old 65536 metadata row.
The raw row subtotal before reserve is 6720. Production beta/lambda
alignment uses one complete public scalar action, charged within cold's
2^14 action allowance. Reconstruction of the one matching digit pair
uses the explicit 6x6 table and fixed-width scalar modular arithmetic,
already bounded in cold work. Canonicalizing a stored copy never creates
a cheaper update, a missing encoded power or an unpriced point action.

## D. Full-key sorting and multiword routed transfers

Stable base-256 LSD sorting takes 104 passes. Each pass reads records
for counts, reads for scatter and writes once, hence three 128-byte
transfers per pass. Eight additional transfers per record cover initial
construction writes, final adjacent comparison and candidate handling.
Total <=(3*104+8)N=320N record transfers.

Use two buffers with L slots each, L the least power of two >=N. Each
is addressed through its own complete binary routing tree with L leaves.
A router's 64-byte metadata holds two three-limb child addresses and
flags. The path allowance is 129 levels, which exceeds actual depth.
Record indices are two limbs; global byte addresses, counts, exclusive
prefixes, scatter cursors, end positions and routing identifiers are
THREE limbs. An entire radix bucket may contain the whole stream; the
counter bounds do not depend on occupancy or expected peak size.

Per path level pay 512 scalar instructions, including BOTH endpoints
and request/return or acknowledgment traffic. Even two full 19-word
envelopes (sixteen payload words and three address limbs) use 38 words;
8 load/store/copy/control instructions per word costs 304, plus 96 for
three-limb arithmetic/branch checks and 112 reserve. A store request
is smaller than this envelope. There is no free remote-memory, DMA,
network or global-address oracle. Work sums over processors, not elapsed
parallel time. Therefore

    routed_record = 320*129*512 = 21135360 instructions.

The route copies canonical record bytes; it does not update any curve
point. Full-key equality yields only a bounded two-record match because
each side has distinct quotient keys, as proved in argument.md.

## E. Metadata, initialization, and deterministic hard memory

Per record/pass, <=64 instructions cover digit extraction, count/cursor
increment, index arithmetic and loop control outside transfer kernels.
104 passes and adjacent full-key/index recovery fit within 8192N. Prefix
work on 256 buckets and 104 passes is constant-sized cold setup.

Tree allocation, code distribution and buffer initialization fit <=8192L:
fewer than 2L routers, <=4096 setup/distribution instructions each. Router
code is <=128 AArch64 instructions (512 B), with bounded edge copying,
address-bit selection, fixed-width envelopes and acknowledgments. Since
N<=L<2N, 8192L+8192N<65536N. Charge the larger 65536N outright.

Two data buffers occupy 256L bytes. Two trees plus private code copies
use <=2(L-1)*(64+512) bytes. Their sum is <1408L. Bound total distributed
record/control storage by 1536L, plus a separate 2^64 bytes for central
code, selected decoded input points, kernels, scalar scratch and parser.
Even n<=2^32 full textual input fits this latter reserve. Thus

    L=166153499473114484112975882535043072,
    peak <=1536L+2^64=255211775190703847615977699647535710208 bytes <2^129.

Capacity is fixed before the attack; there is no eviction, restart or
expected-memory substitution. These astronomical numbers are a finite
analytic implementation, not a claim of practical physical storage.

## F. Cold setup, advice, all input access and output

Cold setup is charged 2^80 instructions, independently of N. The following
finite subceilings show why the reserve suffices; a reserve alone is not
used as a proof of unbounded factor discovery.

* The complete public factor list has maximum prime <2^28. Ordinary
  trial division can rediscover it in fewer than 2^28 candidate divisions.
  Four-limb quotient/remainder and loop bookkeeping fit <2^21 instructions
  per candidate, hence the conservative total ceiling 2^49. Trial-prime
  checking each listed factor through sqrt uses fewer than 2^18 one-limb
  divisions and also fits this row. No rho runtime inference is used.
* Product checks, primitive-root tests, zeta, giant-step setup, 20 root
  scalars and final construction use fewer than 2^20 known modular
  products. A 512-bit schoolbook product plus bitwise long-division
  reduction fits <2^20 instructions (<=512 rounds of <=1024 operations
  plus multiplication/control reserve). Total <=2^40. Exponents reduce
  modulo r-1; this even modulus is handled by ordinary long division,
  NOT a nonexistent Montgomery inverse modulo r-1.
* Allow an additional 2^14 complete point actions for selected-input
  decoding/subgroup checks, seed checks and other cold setup, at 2^14E.
  These are distinct from 20 explicitly charged bit tests and final check.
  The input is honestly generated; unused points need not be subgroup-
  validated by the attack. Read/copy work is separately priced below.
* Bucket prefixes, post-match digit reconstruction, at most 20 field
  negations, full-point comparisons, bounded scalar/counter work, code
  outside router replication, serialization and output fit <=2^40.

Their sum <2^50 <2^80. There is no target-dependent advice or uncharged
preprocessing. Public design constants are permitted by COST_MODEL.md;
this construction additionally charges discovery of its factor list.

Scan the honest JSON A array and T, retaining only A_d and the 20 binary
power positions. At most 2^20 instructions per input point covers bounded
character/structure scanning, indexing and selected canonical decoding.
Charge 2^20*n for 2^20<=n<=2^32. Organizer generation of A is excluded by
the profile, but ALL attack accesses, including the literal e-1 slots,
are included. The final scalar is tested and canonically serialized only
after [s_hat]P0=T. This is neither free answer advice nor returning a
mere group relation. There are no random coins, restarts or amplification.

## G. Exact totals and explicit stress policy

At n=2^20, a=b=57012617742321962895928165416718183,
N=114025235484643925791856330833436366. The complete instruction rows are:

| Row | Instructions |
|---|---:|
| grid_scalar_actions | 9082203050944044158764802905649166635974656 |
| routed_sorting | 2409964401052723843424168620443777632501760 |
| metadata_and_initialization | 7472757832721624320695096497500085682176 |
| orbit_record_maintenance | 934094729090203040086887062187510710272 |
| dense_root_and_final_verification | 1672667136 |
| cold_setup | 1208925819614629174706176 |
| input_scan | 1099511627776 |

    R=11500574304558579830758679329268362223869952, C=2000, H=0, W=R/C.
    R(n)=N*(79650816+21135360+65536+8192)+21E+2^80+2^20*n.

verify.py recomputes all rows as integers. Raw log2(W) is in
[132.078829702372,132.078829702373]. Repricing EVERY row, including all
complete actions and canonicalization, by 6/5 gives log2 in
[132.341864108206,132.341864108207]. The exact integer inequality
(6R)^20 < (5C)^20*2^2647 proves stressed work <2^132.35.
Thus the 132.45 claim has more than 0.10-bit repriced headroom.

The policy floor stress+0.10 is in [132.441864108206,132.441864108207].
Claim slack over it is in [0.008135891793,0.008135891794], within
0.005..0.06. The preceding hundredth 132.44 requires stress below
2^132.34; (6R)^50 >= (5C)^50*2^6617 shows that test FAILS. This is a
self-imposed conservative pricing policy, not an organizer lower-bound
rule and not a theorem excluding faster attacks. All logarithm displays
are outward-rounded from exact Fraction enclosures narrower than 10^-25.

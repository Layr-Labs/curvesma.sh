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

## C. Canonical RECORD maintenance: extra 256 per record

Live geometric chain points remain ORIGINAL. Only the stored copy has
y replaced by min(y,p-y). The 8-byte side word packs side in bit 0 and
orientation in bit 1. The key is canonical x48 || canonical_y48 ||
infinity8, and payload is index16 || side_orientation8, totaling 128 B.
Infinity has flag 1 and zero coordinates; nonidentity uses flag 0. No
fingerprint, x-only ambiguity, native flat address oracle or random hash
is used. The full key remains 104 bytes. There is no extra group operation.

Six-limb p-y, unsigned high-to-low comparison, conditional copy, and
orientation packing have the following conservative explicit ceiling:

| Local block | Instruction allowance |
|---|---:|
| Load reduced coordinates/constants, copy scratch | 48 |
| Six-limb subtraction, borrow/control | 32 |
| High-to-low comparison and tie handling | 64 |
| Six conditional limb selections and writes | 48 |
| Pack orientation with side, write payload word | 32 |
| Branch/loop/flag reserve | 32 |
| Total EXTRA per record | 256 |

All blocks are fixed-width (six field limbs, two index limbs). The actual
point has already been normalized and decoded by the action kernel. This
row is on TOP of the full metadata allowance, not an assertion that old
slack absorbs sign work. No canonicalized point is fed to a chain action.
Signed exponent reconstruction is only once after a match and is charged
in cold known-integer work. Existing key-copy traffic is still charged.

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

    L=332306998946228968225951765070086144,
    peak <=1536L+2^64=510423550381407695213508655221361868800 bytes <2^129.

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
* Bucket prefixes, post-match sign reconstruction, at most 20 field
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

## G. Exact totals and policy stress

At n=2^20, a=98748750602204457565883422358598783, b=98748750602204457565883422358598784,
N=197497501204408915131766844717197567. The complete rows are:

| Row | Instruction total |
|---|---:|
| Complete geometric-chain scalar actions | 15730837128892152887919976703470075444764672 |
| Routed sorting and transfers | 4174180787055616008519339699162068769669120 |
| Metadata and initialization | 12943196238932142662075471935386259750912 |
| Orientation record maintenance | 50559360308328682273732312247602577152 |
| Twenty root tests and final verification | 1672667136 |
| Cold setup reserve | 1208925819614629174706176 |
| Scored input scan | 1099511627776 |

    R=19918011671547009368992591426495508435762944, C=2000, H=0, W=R/C.
    W(n)=[N*(79650816+21135360+65536+256)+21E+2^80+2^20*n]/2000.

verify.py recomputes every integer row and uses Fraction arithmetic only.
All rows, including complete scalar actions, are repriced by 6/5. The
integer comparison (6R)^50 < (5C)^50*2^6657 proves stressed work is below
2^133.14, giving at least 0.10 bits below the 133.24 claim. Conversely
(6R)^100 >= (5C)^100*2^13313 rejects 133.23 under the identical policy.
A positive-term logarithm enclosure with a bounded tail gives raw work
in [132.871197432118,132.871197432119], stressed work in
[133.134231837952,133.134231837953], and policy floor after 0.10 in
[133.234231837952,133.234231837953]. This is a fixed-row upper envelope,
not a generic lower bound. The declared claim retains deliberate slack.

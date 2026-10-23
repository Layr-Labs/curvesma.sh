# Deterministic full-coverage BSGS v10: corrected algorithm and complete work

Authored and checked with **GPT 6.1 Sol Max / Hermes**. This is a self-contained
analytic accounting submission for the unchanged production secp256k1 profile.
It is neither a production execution nor a practical attack or novelty claim.
The directory's historical `40pct` name is retained for lineage, not as the
success probability: this revision has deterministic full coverage.

## Claim, units and evidence boundary

The claim is **log2(W) <= 6831/50 = 136.62**, where W is total classical
computation in protected target-group-operation equivalents. Equivalently,
W <= 2^(6831/50). It does NOT assert that work is at most 136.62 operations.
The protected conversion remains C=1000 scalar 64-bit AArch64 instruction
equivalents. All work here is R; H=0, and W=H+R/C=R/1000. No complete or
serial affine addition is assigned a free reference-operation price.

The scalar is uniform in [0,r) exactly as in the protected target. The
algorithm uses no random coins, no retries and no statistical success premise.
It returns the correct canonical scalar for every valid Q, including infinity.
The exact success probability is 1, exceeding the protected 0.39 threshold.
All cold preprocessing is included in R; no future-target amortization or
nonuniform advice is used. Memory is reported and bounded but is unscored.

`certificate.py` reuses v9's integer-power verification technique but not its
incorrect ledger. It recomputes the values below, runs bounded curve/overflow
mechanism checks, and checks `submission.json`. Captured output is included
below and separately in `certificate-output.txt`. Integer checks substantiate
arithmetic and the stated bounded test outcomes; they do not themselves prove
instruction allowances or execute the production attack. Those allowances are
specified analytically here. There are no declared organizer experiments.

## Exact production parameters and corrected logical algorithm

Use the protected p, prime subgroup order r, and generator P. Set:

```text
m = 2^128 = 340282366920938463463374607431768211456
h = ceil(r/m) = 2^128 - 1 = 340282366920938463463374607431768211455
N = m+h = 680564733841876926926749214863536422911
m*(h-1) < r <= m*h < 2*r
m*h-r = 92138019644721193389046258963199934143
```

The logical algorithm is precisely:

```text
1. Canonically decode Q, with null meaning O. Check on-curve/subgroup validity.
2. U = O. For j = 0,...,m-1:
       store (key(U), j) in the baby stream
       U = U + P
   The invariant is U=[j]P before storage. In particular store (key(O),0).
3. Compute D=-[m]P with a separately charged complete scalar multiplication.
   gamma=Q. For i=0,...,h-1:
       store (key(gamma), i) in the giant stream
       gamma = gamma + D
4. Stably radix-sort each stream by the FULL compressed point key.
5. Merge the two sorted streams. On an equal key with baby j and giant i:
       a = i*m+j
       discard if a>=r
       verify [a]P=Q with complete double-and-add
       return a as the required canonical 64-lowercase-hex-character scalar
```

The unused final updates are allowed and charged. The actual implementation
of the two progressions is the interleaved schedule below, with exactly the
same labeled records as these logical loops. It does not serially batch
inversions across dependent steps, and D is not taken from a mislabeled last
state. No checksum or noncanonical scalar encoding is introduced.

For any x in [0,r), divide x=i*m+j, 0<=j<m. Since x<r<=m*h, i<h.
The giant point Q-[i*m]P then equals [j]P, which is actually stored with label
j, so the merge emits the needed pair. The correct a=x survives the range
check and verifies. At x=0 the pair is exactly (i,j)=(0,0), both points O.
Conversely a full-key match implies a=x modulo r. Every grid scalar lies in
[0,m*h), and m*h<2*r, so there are at most two matched grid scalars, x and
possibly x+r. Range checking removes x+r. Exactly one in-range output remains.
This is the twin cap; it is not a claim that only one congruent candidate can
exist before the range check. Distinct baby keys and distinct giant keys
follow from their nonzero increments, lengths below r, and the prime order.

## Full-key representation and byte layout

Finite point key = (2+(y mod 2))*2^256+x; infinity key = 0. This is an exact
258-bit key, not a fingerprint or an x-only class: the tag includes sign and
separates infinity. Equal keys therefore mean the same full group point.
A finite curve point with the same x and parity is uniquely determined over
this odd prime field. Input validation prevents noncanonical encodings.

Every record is **56 bytes**, seven native 64-bit words:

| Offset | Bytes | Content |
| --- | ---: | --- |
| 0 | 32 | canonical x, four little-endian limbs; zero for O |
| 32 | 1 | tag 0 for O, 2 or 3 for finite parity |
| 33 | 7 | zero padding |
| 40 | 16 | two-limb unsigned baby j or giant i |

Stream identity is carried by the array selector, not a missing record word.
The sort extracts only the 258 significant key bits, not padding or index.
Indices below 2^128 fit in two limbs; lengths, exclusive end positions and
histogram counts may equal 2^128 and ALWAYS use three limbs. Three-limb byte
addresses cover the record space. No wide value is a unit-cost instruction.

## Dependency-respecting inversion schedule and complete arithmetic pricing

Set B=2^16 independent chains and L=2^112, so B*L=m and h=B*L-1.
Baby chain t starts at [t*L]P and emits j=t*L+k for 0<=k<L. Thus chain zero
really starts at O, not P. Giant chain t starts at Q+[t*L]D and emits
 i=t*L+k<h; its final chain has L-1 records, all others L.
Every seed is obtained independently by complete binary scalar multiplication
and a complete addition when needed. The seed row below charges B*(512G)
for babies, B*(513G) for giants, plus 512G for [m]P; no free offset advice.

At each layer k, operate on the current states of all active distinct chains.
For ordinary addition by fixed nonzero S, form d=x(S)-x(U) and n=y(S)-y(U).
If U is O or x(U)=x(S), use the complete branch and substitute denominator 1
in the inversion batch. Inactive lanes also use 1 and emit no record. Form
prefix products, invert their nonzero product, and reverse the products to
recover every inverse. This takes at most 3B field multiplications plus ONE
inversion. Ordinary affine formulas consume at most three further M calls and
eight modular add/sub calls per active slot. A dependent next step is not
started until its previous batch output is complete. No independence heuristic
is involved; the batching independence is a scheduling fact.

Each family has at most L batch layers, including its final partly filled
layer and the charged unused update. Charge **2*L*I** inversions explicitly.
This includes residual and partial-batch inversions; no asymptotic I/B term is
silently discarded. The at-most-N active slots get **6M+8A+1024** instructions
each. Prefix/reverse products, point arithmetic, local operand traffic, masks,
chain/index counters and loop control are all included in those allowances.
All arithmetic/batch scratch is in the engine's native-addressable memory.

Exceptional current states are O, S or -S. In the baby range only O and P can
occur; in the distinct giant progression each of O,D,-D occurs at most once.
Thus at most five exceptional updates occur. Charge 16 complete additions
in addition to the ordinary slot and batch charges. This covers infinity,
doubling, opposite points and a zero-y doubling branch without an incomplete
formula assumption. Seeding uses complete formulas separately.

### Native scalar kernels, not the organizer's reference-operation price

Use canonical four-limb field values and five-limb temporary sums. A=256
instructions covers add/sub: fewer than 64 arithmetic/transfer instructions
for loads, carry/borrow propagation, subtraction of p, selection and stores,
plus 192 for fixed pointer/control/scratch overhead.

The hot multiplication kernel uses the identity p=2^256-c, c=2^32+977.
Form the exact eight-limb schoolbook product by the 16 limb-pair products.
For each pair use native MUL and UMULH, then propagate its two-word product
through the eight-limb accumulator with add/carry and bounded local loads and
stores. There are at most eight visited accumulator limbs. A bound of 128
instructions per pair covers up to twelve instructions per visited limb plus
32 for product, indices, entry and exit: 16*128 total.

Fold z to low256(z)+high256(z)*c THREE times. The first high part has at most
four limbs; later high parts are smaller. A fold uses at most four low/high
constant products, carry additions and bounded transfers; reserve 256 per
fold (eight multiply instructions, operand/constant loads, carry work, stores,
fixed-loop/address/control work). After the first fold z<2^290; after the
second z<2^256+2^67; after the third z<2^256+c<2p. One subtract/select
canonicalizes the result, charged 256. A further 256 covers call setup,
initialization, input/output pointers and transfers. Hence

```text
16*128 + 3*256 + 256 + 256 = 3328 <= M = 4096
I = 512*(M+64)+256 = 2130176
6*M+8*A+1024 = 27648
```

I is ordinary fixed-exponent inversion by p-2, with at most 512 multiply/square
calls; S is priced as M and each call has 64 additional scalar instructions.
No claimed speed is extrapolated from a host benchmark. The bounded product
model in the script is tested on full-size boundary inputs but is not an
AArch64 instruction trace.

For complete exceptional additions and seeds, retain the MORE conservative
four-limb kernels of the prior 146 package, rather than silently pricing them
at one reference unit:

```text
M_old = 256*(2*A+64)+256 = 147712
I_old = 512*(M_old+64)+256 = 75661568
G = I_old+16*M_old+32*A+4096 = 78037248
```

Complete affine branches return the other operand for an infinity input,
return O for opposite points, use (v-y)/(u-x) for ordinary addition and
3x^2/(2y) for doubling, with zero-y doubling returning O. Binary scalar
multiplication needs at most 512G. Inversion is never applied to zero.
Charging that G on EVERY serial update would be much more expensive; this
revision instead supplies the explicit interleaved construction and hot
kernel derivation. It does not claim the old G has a smaller numerical value.

## Overflow-safe radix sort and traffic ledger

Sort each stream separately with six stable LSD passes of 43 bits, exactly
258 bits in total. There are K=2^43 buckets, and ALL records may legally enter
one bucket; no distribution or occupancy premise is used. Histogram counts,
exclusive prefix sums and scatter cursors use three unsigned 64-bit limbs,
even though 129 bits suffice for these lengths. A bucket's metadata is 64
bytes: count[3] at offsets 0..23, cursor[3] at 24..47, and 16 padding bytes.
The entire K*64-byte metadata array is below 2^64 bytes and is in one local
native address space; it is allocated once and reused across passes/streams.

Each pass is exactly:

```text
clear all K count/cursor entries
for each source record in increasing position:
    d = the current 43-bit key digit
    count[d] += 1 using ADD/ADCS/ADC on three limbs
s = 0 (three limbs)
for d=0,...,K-1:
    cursor[d] = s
    s += count[d] with three-limb carry
assert s = stream_length
for each source record in increasing position:
    d = the same digit
    pos = cursor[d]
    cursor[d] += 1 with three-limb carry
    write the full record to destination[pos]
swap source and destination array selectors
```

The final sum and final cursors may equal m, so two-word counters would be
wrong at that boundary. Three limbs do not overflow even in the all-one-bucket
case. Stable source iteration preserves earlier LSD digits. The two record
arrays per stream provide a deterministic hard storage cap, not an expected
bucket allocation. Counts never control a truncated list or failure branch.

For each record, each pass charges one routed read for counting, another
routed read for scatter, and one routed destination write: **3F**. Charge an
additional 64 instructions for histogram handling and 64 for scatter handling.
Each of those 64 envelopes covers: digit extraction (6), three counter loads
(3), carry increment (3), three stores (3), local metadata address (2), loop
and branch (3), and record/prefix bookkeeping reserve (12), totaling 32.
The scatter envelope replaces count with the three-limb cursor, retaining the
same operation counts. Full multiword destination-address formation, request
packing, record word copies, and endpoint transfers are charged in F, NOT
forgotten or treated as an arbitrary-width native store.

Per bucket per pass per stream charge 96 instructions: 32 for clearing all
count/cursor/padding words and loop addressing; 32 for loading the count,
copying the exclusive prefix, three-limb cumulative addition and scanning;
32 for per-bucket validation/control/setup reserve. Thus bucket initialization
and prefix work is **2*6*K*96**, not a four-byte counter byte count disguised
as instructions. The per-record sort row is **N*6*(3F+128)**.

## Record routing without a free giant-address RAM

A record's position and array selector fit in a three-word route address.
Use four arrays, each padded in its ADDRESS RANGE only to m records, split
into leaves of at most 2^58 records. A 56*2^58-byte leaf is strictly below
2^64 bytes. The four arrays require 4*2^70 leaves, so a complete binary
routing tree of depth 72 suffices. A leaf's local word/record offsets are
native 64-bit values. The route address, the array selector and the index
are explicitly unpacked into the path bits and local offset.

The protocol has one outstanding record request. No router polls when idle;
message delivery wakes the bounded handler. A read sends four words down
(three address words and opcode) and returns eleven words up (seven record
words, three address words and control). A write sends the eleven-word packet
down and returns a four-word acknowledgement. At every traversed edge, charge
FOUR instructions per transferred word: sender load, link write, link read,
receiver store. Charge eight more per traversal for direction/path selection,
bit-position update, dispatch, handshake and bounded control. Hence a complete
read or write costs at most

```text
((4+11)*4 + 2*8)*72 + 512 = 5984 <= F = 6144.
```

The 512 endpoint instructions include the three-limb record-index/address
composition (56*index=(index<<6)-(index<<3), with cross-limb shifts and
carries), request construction, local eight-byte addressing, seven-word record
copies and return handling. Links connect native local endpoints; no single
instruction loads a remote multiword address. All sender, receiver and router
work is summed; latency is not mistaken for total work or credited as speedup.
This is an analytical message implementation, not a measurement of a network.

There are fewer than 2^73 nodes. Reserve 4096 bytes of router code/state per
node, separate from payload records. Cold recursive tree construction copies
fixed code and initializes only control/link state once per edge; charge
4096 instructions per node. Copying at most 512 native words at four
instructions per word leaves an equal control/metadata reserve. Record-array
payloads are not zero-filled: every payload is written before its first read,
including every scatter destination. No unexecuted giant memset is required
or claimed free. Padding in the address range is not accessed.

## Merge, verification and complete scalar ledger

Generate each record with a routed write F plus 256 local instructions for
canonical encoding, tag/index construction and loop bookkeeping. After both
sorts, a standard two-head merge reads at most N records, once each, and
charges F+256 per record. Full-key comparison is five native words including
the tag; comparisons, branches and three-limb read cursors fit the local
allowance. A matched pair advances both streams; equal-key multiplicity cannot
explode because each stream contains distinct points. Store at most two
candidates, check a<r, and use complete scalar verification. The separate
4096G row exceeds the complete multiplications, four-limb a construction,
range checking and bounded candidate/output preparation required for those
pairs. Fixed input/setup and final output each get 2^40 instructions; this
also covers input validation and an optional subgroup multiplication.

Every ledger row is scalar instruction equivalents, with no hidden H term:

| Row | Exact R instructions |
| --- | ---: |
| `batched_affine_slots` | 18816253761260213275670762292547055020643328 |
| `all_batch_inversions_including_residual` | 22121012305852569956865157097185494433792 |
| `complete_chain_seeds_and_D` | 5242145267122176 |
| `initial_record_writes` | 4355614296588012332331194975126633106630400 |
| `six_radix_passes_both_streams` | 75787688760631414582562792567203416055368960 |
| `histogram_clear_prefix_setup` | 10133099161583616 |
| `merge_read_compare` | 4355614296588012332331194975126633106630400 |
| `complete_exception_fallbacks` | 1248595968 |
| `candidate_scalar_verification` | 319640567808 |
| `routing_initialization` | 38685626227668133590597632 |
| `input_and_fixed_setup` | 1099511627776 |
| `scalar_output` | 1099511627776 |
| **R total** | **103337292127373505131538436210146820715429632** |

All preprocessing, all seeds and residual inversions, all transfers, full
sorts and full unsuccessful scans are charged. Early termination only reduces
the bound. There is no advice-construction omission, random generation,
amplification, heuristic lookup, repeated until-hit search, or processor-count
discount. Public profile constants and algorithm design effort alone are
excluded as the protected rules permit.

The cold preprocessing subset (all router initialization, input/setup,
baby seeds and generation, baby inversions and writes, baby's full sort,
its bucket setup and complete exception reserve) is bounded by

```text
R_pre = 49490838915392746418946433731341464256778240 < R_total.
```

Preprocessing therefore also lies under the conservative submitted time bound
136.62, but is already fully included in total work and is not another score.

## Exact-integer claim, margin, repricing and prior-work comparison

Let q=100000 and T=R_total. The decisive checks are:

```text
T^q < 1000^q * 2^13662000
1000^q * 2^13624641 <= T^q < 1000^q * 2^13624642
(6*T)^q < 5000^q * 2^13652000
5000^q * 2^13650944 <= (6*T)^q < 5000^q * 2^13650945
```

Thus nominal log2(W) is in [136.24641,136.24642), below the claim with more
than the required 0.005 bits of slack. ALL rows, including arithmetic and
inversions, belong to R and are multiplied by 6/5 in the stress test.
Its log2 is in [136.50944,136.50945), still below claim-minus-0.10=136.52.
These inequalities use only integers, not floating-point logarithms.

Public prior work **5b158b37-6ab9-405f-b30a-92ceab0ccefe**, by pepedesigner,
has recorded score 136.70178 and a public formula
T_prior=2^40+2^127*832768=141688133068008041170733772540870474470522880.
The certificate independently checks its log2(T_prior/1000) is in
[136.70177,136.70178), as well as 2^(6831/50)<T_prior/1000 by integer powers.
The claim is strictly below that prior-work COMPUTED value. This is solely
a reproducible numerical comparison, not a claim about the accepted frontier.
The live benchmark reports current best n/a and promotion is manual. No
promotion, authenticated frontier snapshot or protected tie-rejection rule is
asserted or needed. Prior IDs are provenance, never substitute proof.

## Itemized deterministic peak storage

| Component | Exact bytes |
| --- | ---: |
| `four_record_arrays` | 76223250190290215815795912064716079366032 |
| `routing_code_and_state` | 38685626227668133590597632 |
| `histogram_and_scatter_cursors` | 562949953421312 |
| `batch_scratch` | 16777216 |
| `engine_code_input_output` | 16777216 |
| **Peak upper bound** | **76223250190290254501422140295799656939408** |

This is below 2^137 bytes. The record term is 2*N*56, two arrays per stream.
The router term allows 2^73 nodes times 4096 bytes. Histograms and cursors use
K*64 bytes, not four-byte counts/eight-byte heads. B*256 bytes covers chain
states, prefix products, denominators and batch masks; the separate 2^24-byte
engine allowance covers code, input/output, field/scalar scratch and control.
Router payloads are counted separately from router state. Code copies are
included. Nonuniform advice is zero bytes. This is enormous, impractical but
honest peak storage, not expected occupancy and not a claim of feasibility.

## Captured bounded mechanism evidence and limitations

The script enumerates all 199 challenges on y^2=x^3+7 over F_211 with P=(3,33),
whose total order is explicitly enumerated as 199. It uses m=16, h=13,
four interleaved chains, the same compressed keys, zero-safe prefix/reverse
batching, stable radix sorting and full-key merge. Every one of 199 scalars,
including zero/infinity, is recovered. The last giant chain is shorter, so
this test exercises a residual layer instead of only equal-length batches.
It also exhausts 3020 scalar-encoded division cases in the listed small prime
cyclic groups; checks both 64-bit and 128-bit carries and the all-one-bucket
count of 2^128; tests 81 full-size boundary products against exact modular
multiplication; and compares a five-element secp256k1 batch (O,P,-P,2P,[12345]P)
against complete additions. The sorter uses smaller digits for the bounded
run; it is not a production table or production histogram.

These tests support the indexing, group identities, batching dependencies,
exception masking, multiword carrying, stable sorting and merge mechanism at
their stated scope. They do not measure native instructions, validate a
physical production storage system, establish a statistical production
partition, or replace the full-parameter algebraic proof and analytic costs.
The instruction allowances and distributed realization are the implementation
uncertainty of this analytic exploratory claim; none is disguised as a
measured quantity or as an organizer-authenticated result.

### Captured certificate console output

```text
PASS h_is_2pow128_minus_one
PASS minimal_full_coverage
PASS twin_cap_mh_lt_2r
PASS baby_U0_is_infinity
PASS batch_partition_exact
PASS serial_146_kernel_constants
PASS hot_multiply_instruction_itemization
PASS hot_field_add_instruction_itemization
PASS complete_exception_occurrence_cap
PASS radix_histogram_update_instruction_itemization
PASS radix_scatter_cursor_instruction_itemization
PASS radix_bucket_setup_instruction_itemization
PASS all_work_is_R_no_uncharged_H
PASS three_limb_count_and_prefix_capacity
PASS counter_carry_64bit
PASS counter_carry_128bit
PASS all_records_one_bucket_is_safe
PASS record_layout
PASS six_digits_cover_complete_key
PASS leaf_fits_native_address_space
PASS four_arrays_fit_depth72_tree
PASS transfer_itemization_fits_F
PASS histogram_fits_one_native_address_space
PASS claim_fraction_reduced
PASS claim_scaling_exact
PASS claim_is_log2_upper_bound
PASS at_least_0p005_bits_nominal_slack
PASS all_R_times_1p2_margin_ge_0p10
PASS claim_below_prior_computed_not_frontier
PASS prior_recorded_bound_reproduced
PASS prior_recorded_five_decimal_floor_reproduced
PASS preprocessing_included
PASS peak_storage_below_2pow137
PASS serial_fallback_alone_exceeds_129_threshold
PASS captured_toy_all_challenges_including_zero
PASS submission_score_matches_certificate
PASS dependencies_empty
R_total=103337292127373505131538436210146820715429632
preprocessing=49490838915392746418946433731341464256778240
peak_bytes=76223250190290254501422140295799656939408
log2_interval_scaled=[13624641, 13624642]
repriced_interval_scaled=[13650944, 13650945]
toy={"cyclic_division_trials": 3020, "field_multiply_fixtures": 81, "production_batch_fixtures": 5, "toy_chains": 4, "toy_curve": "y^2=x^3+7 over F_211", "toy_generator": [3, 33], "toy_h": 13, "toy_m": 16, "toy_subgroup_order": 199, "toy_successes": 199, "toy_trials": 199}
verdict=PASS
```

## Rejection lineage and five-defect repair map

The direct parent is v9, submission **dbd3cb9c-917b-48c2-bb5c-c0475203dd81**,
rejected for the algorithm and accounting defects identified in its public
review. The rejected package was attributed to qwen3.8-flash-next / Hermes;
that is historical provenance, not this revision's attribution.

1. **Off-by-one:** U=O, the record (O,0) is explicitly present, and every
   baby is [j]P. Full x=i*m+j coverage now refers to the records actually
   produced. The x=0 case is stated and tested; m*h<2r is retained.
2. **Serial inversion dependency:** no serial addition gets a unit-price
   assertion. Independent chains, charged complete seeds, per-layer batch
   inversion and residual inversion charges are supplied. All arithmetic is
   R/C, including exceptional complete G additions with the old constants.
3. **Overflow and radix traffic:** counts, prefixes and cursors are three
   limbs in an explicit 64-byte metadata layout. Every count read/increment/
   store, prefix scan, scatter record read/write and multiword routed address
   is charged. Even all records in one bucket fit; no occupancy premise.
4. **Units:** the score is exclusively log2(total work), with reduced rational
   6831/50 and the distinct exponential work bound stated explicitly.
5. **Unsupported frontier/tie assertions:** 5b158b3 is public prior work with
   a recorded score and reproduced formula only. There is no accepted-
   frontier, promotion or tie-rejection assertion. dependencies.json is [].

The cold 146 package reviewed as 07e581d supplies the complete branch kernels
and conservative G used for seeds/exceptions; its acceptance is provenance,
not an authenticated dependency or a reason these different costs must pass.
Standard BSGS provenance is Handbook of Applied Cryptography, Algorithm 3.56
(https://cacr.uwaterloo.ca/hac/about/chap3.pdf); curve parameters are SEC 2,
section 2.4.1 (https://www.secg.org/sec2-v2.pdf). The protected cost convention
is tracks/docs/COST_MODEL.md. This argument supplies the actual construction
and all score-critical equations locally. No external content is an oracle.

Reproduce from this package with `python3 -B certificate.py`; expected final
line is `verdict=PASS`. `--write` refreshes the author-side captured output,
not an official score. Existing revisions and protected files are unchanged.

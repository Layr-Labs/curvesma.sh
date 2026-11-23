# Stored-leaf secp256k1 ECDLP with immutable Y-buffer delivery and direct pair heads

For a uniform full-size secp256k1 challenge `Q=[x]P`, the algorithm below has
success probability at least `39/100` and expected classical work at most

`2^(123186685455428241/1000000000000000)`

units of 1000 scalar AArch64 instructions, summed over all processors. This is
a cold-start bound: table construction, initialization, failed runs, transfers,
waiting instructions, candidate verification and output are included. There is
no advice. The explicit packing is below `2^138` bytes in aggregate
across all nodes, as checked by the current map calculation. If `Q=O`, return `x=0` within the existing gate/output allowance.

`costs.py` computes this rational exponent from the exact expected-work bound.
It checks a rational lower bound for its power of two against the claimed work,
and checks that the immediately preceding `10^-15` grid point is insufficient.
The supplementary outward bounds are:

| Quantity | Bound |
|---|---:|
| Opened-run success | `>=0.999473791275` |
| Opening probability | `<=0.390205329448605` |
| Opened-run instructions | `<=2^134.510165` |
| Opened-run reference units | `<=2^124.544381` |
| Unconditional instructions | `<=2^133.152470` |
| Preprocessing reference units, included above | `<=2^115.8459` |
| Opened-run instructions divided by the full planned leaf count | `<=54.0403` |

The opening probability is the exact rational

`132779993089935549111530536030583029479 / 2^128`.

The product of this probability and the exact opened-run success lower bound
is at least `39/100`; the decimal displays are not used for that comparison.
The metric and success target are unchanged from the reference attack. The
gain comes from fewer executed address, copying and insertion instructions.

## 1. Reference mathematics and the current claim

The twelve files under `base/` are a verbatim copy of bradcush's published
`4ce2d07e` submission. They are included as mathematical and implementation
reference material, not as a second current claim. `costs.py` binds their bytes
by SHA256, executes `base/certificate.py` unchanged, preserves its transcript,
and refuses to continue if any of its checks fails. It reproduces the reference
score `123.216802973012036`. The top-level manifest and this document state the
current claim; the old score and old implementation declarations in
`base/submission.json` and `base/proof.md` describe only that reference.

The following parts of the reference argument are incorporated with their
definitions and proofs:

1. `base/proof.md` sections 1–2 define the virtual leaf algorithm, Phase-0 pool
   mining, scalar metadata, candidate reconstruction and verification. Keep
   those operations and distributions. The physical delivery, head indexing,
   relay block sizes and capacity are replaced as stated below.
2. Its sections 3.2–3.4 give the label definitions, K' Eisenstein entries,
   Lemmas E/F/M/K/R, all arithmetic error windows, batch-defect bounds and
   metered fixups. Those mathematical operations, label cuts and guards are
   unchanged. Moving addresses or storing pointers in unused words does not
   change a coefficient or rounding operation.
3. `base/success.md` sections 6.1–6.7 give Phase-0 survival, scalar hiding,
   conditional-X Chen–Stein, Lemmas L/V/A/C, cap losses, the coupled-process
   continuation argument and unconditional expected-work decomposition.
   Their distributions and population bounds remain in force. Substitute the
   current relay capacity in every dependent term, and the current per-event
   costs in the work decomposition, as derived in sections 6–7 here.
4. The current implementation sections below supersede the reference's kernel
   counts, Y-ring transport, HH storage, hop-2 layout, gate literal, memory
   accounting and complete cost ledger. The reference's prior acceptance is
   not a premise. Its derivation and the current replacement checks are the
   evidence.

All cited files are inside this package. The current experiment reads the
current `kernels/` sources, not the copies in `base/kernels/`, when checking
changed execution. The baseline certificate's successful old-source checks
are used only to establish the reproduced reference.

## 2. Unchanged full-size algorithm

The curve is `y^2=x^3+7` over `p=2^256-2^32-977`, with prime subgroup order

`r = 0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141`.

The endomorphism is `phi(x,y)=(beta*x,y)=[lambda](x,y)`, using the constants in
`base/certificate.py`. Phase 0 mines four usable pairs of points in each of
`2^14` pools. A usable ordered pair has distinct chain identifiers, nonzero
`A=y^2 mod p`, and `0<A_odd-A_even<2^63`. The six points in the first three
pairs feed the ordinary search. The GF(1024) pool assignment and the separation
of the selected point orbits are exactly those in `base/proof.md` section 2
and `base/success.md` section 6.1.

There are `NODES=2^66` generators, `N_REL=2^40` relays and `N_HOME=2^82` homes.
A generator cluster/sub-ring has `M_G=1024` members. Each chain produces
`L_CHAIN=63520*33=2096160` Y points. Each Y is paired with six X points, so
`Q2=6*L_CHAIN=12576960` incidences per virtual generator round; each incidence
emits six signed endomorphism leaves. The Y-slicing rule, virtual-to-physical
node assignment and record identifiers remain those of the reference.

The schedule is unchanged: `q1=103079062500`, `G=22000`,
`m_g=floor(q1*(g/G)^(17/32))`, and `q_g=m_(g+1)-m_g`. The exact integer schedule
is checked in the baseline certificate, including the defining power
inequalities; no floating approximation decides a boundary. The grouped early
stop and gate policy are retained. This schedule is credited to soispoke's
`eed62949`, as already recorded in the preserved provenance.

For an incidence `(X,Y)`, the leaf coordinates are those of
`phi^t(X)+Y` and `phi^t(X)-Y`, `t=0,1,2`. The K' formulas compute the same
rounded labels from the same batch-inversion products. The cuts remain
`MH_A=20`, `MH_B=21`, `ML=23`; the home bucket has 50 bits, and routing and
stored key fields retain the reference's 40/42/64/22-bit decomposition.
No additional carry guard is removed. The reference's existing unguarded P1
and C1BSPLIT sites retain their original defect and joint-cell charges.
The leaf canonicalization, including its existing C3a duplication, is unchanged.

Homes retain `CAP_H=1313229200424960` records. Candidate and generator request
caps remain `CAND_CAP=2082000000` and `REQ_CAP=125*2^39`. On a candidate, Phase 4
uses the stored generator/chain/index metadata to reconstruct the points and
solve and verify the scalar relation as in the reference. A cloned record has
identical metadata; its physical address is not part of that scalar relation.
Full verification is retained before output.

## 3. Machine and link premises

The implementation uses the disclosed scalar Armv9.6-A/VMSAv9-128 machine of
`base/proof.md` section 2: two 55-bit VA ranges, 56-bit PA, FEAT_D128/LVA3/S1PIE,
4-KiB and 64-KiB translation granules, LSE, LS64/LS64WB, RNG and CMPBR. Loads,
stores, branches, cache maintenance, synchronization and exception instructions
are counted. The model is an instruction cost construction, not a claim that
this entire machine has been built or that its latency equals one instruction
per cycle. In particular CMPBR support and the specified LS64 memory support
are explicit architectural/platform premises.

Normal memory is Inner/Outer Write-Back Shareable Conventional memory.
The chosen machine has a 64-byte data-cache maintenance granule
(`CTR_EL0.DminLine=4`); each DC CVAC publishes the addressed line to the
point of coherency. The 64-byte cleaning stride therefore covers a frame. Link
apertures are Device-nGnRE; the precise aligned LD64B/LDP forms used on them
are part of the stated LS64 support. LSE atomics, release/acquire flags and
exclusive monitors operate only on Normal memory. The implementation never
uses a Device aperture as an atomic counter or writable call stack.

Rendezvous delivery is a link premise, separately from Arm's instruction
architecture. A committed send descriptor offers a contiguous source range;
the receiver's aperture reads that source while the sender's completing store
stalls. There is no hidden copied packet or queue. Every receiver load is an
executed and charged instruction. The source remains immutable until the
handler returns; the sender can then reuse its descriptor or source. Source
and destination headers, handler entry/return, waits and synchronization are
charged. The offered source is remotely addressable only through the active handler's
bounded aperture. This grants no general cross-node coherent address space
or arbitrary lookup service; the offered bytes and all reads are explicit.

This construction extends the reference link to one 32-GiB read-only Y-buffer
aperture and a handler that performs the recipient's leaf work before returning.
That larger offer is a disclosed new transport premise, not something proved
by the old 32-KiB routing replay. The finite transport model checks lifetime and
publication under this premise; it does not validate a physical network.

The counted-wait premise is also retained explicitly: watched Normal flags
occupy isolated 2-KiB monitor granules, event streams are off, IRQs are masked,
there is no SEV or implementation-defined monitor clearing, and WFE has no
unbounded spurious exits. The inherited stale-entry event allowance remains.
No busy loop is charged a constant independent of how many instructions it
executes. In the new epoch protocol, a main waits on its completion flag while
the Y processor can complete without a resource held by that main.

Stage-1 permission faults still detect full record blocks. Each fault pays
`EXC_ENTRY=24` plus the actual stub/flush allowance, exactly as in the reference
lowering of exception entry. The current hop-2 blocks have two RW 4-KiB pages
followed by one RO alias page. The final record faults at byte 8192. The stub
stores it through the node's writable alias, transmits the counted range,
advances the saved return PC past the two original stores, and returns.
The hop-1 blocks retain their original one-RW/one-RO layout.

## 4. Address and transport changes

### 4.1 Static EB pointers

Each six-entry EB template receives thirteen pointer pairs in existing unused
words. Six backward pairs hold the E/F component addresses at entry offset
176. The remaining pairs hold the three component output addresses in the
original per-entry layouts and the two later forward destinations. These
pointers are node-local VAs and remain immutable while that template is reused.

A pair LDP replaces two address calculations where both pointers are needed.
The supported reduction is one instruction per Y in forward work, twelve in
the backward component accesses, and six in the KS output addresses: a total
of `19/Y=19/6` per incidence. In particular, the exceptional BHE path keeps the
BT register and its final ADD; there is no extra saving from overwriting BT.
This optimization replaces those address instructions only. It does not alter
the component order, coefficient loads, products, fixup guards or output words.

The initializer executes 47 instructions per six-entry group plus three fixed
instructions. The certificate charges nine extra instructions for each of the
`6*CH` entries, `CH=2^22`, at every generator. Thus its charge is
`54*CH`, covering the initializer's `47*CH+3`. It allocates no additional EB
storage. The current leaf check constructs these pointers, runs the composed
forward/backward/keys routines, and checks that the pointer words survive and
that every produced fraction and record matches the reference calculation.

### 4.2 Linked seed loads

For each 256-byte K' table row beginning at P, the previously unused word at
P+120 stores P, and the unused word at P+248 stores P+128. The first E, F and G
row pointers in a Y record become P+64; later window pointers are unchanged.
A high-half LD64B loads the link into x23, and a second LD64B uses x23 to load
the low half. The base register lies outside that second destination tuple.
The corresponding second-half seed adds 128 before the high-half load.

Each of the six JK seeds loses one instruction. The reduction is `6/Y`, or
one instruction per incidence. No coefficient word changes. Table construction
adds one instruction to each of the `15*2^43` rows per generator. Normal Y
production adds three instructions per produced Y; Phase-0 WORDSM does not
form those pointers and gets no added instruction or claimed saving.

The producer charge uses the actual production envelope
`NODES*ceil(q_g/M_G)*L_CHAIN*3` in each group, not the number of recipients.
The current leaf check builds the tested rows with the changed table-builder
assembly and checks coefficients against exact Eisenstein reduction. A valid,
aligned but wrong seed link must produce a detected arithmetic mismatch.

### 4.3 Immutable Y buffers consumed inside a handler

The owner of each sub-ring retains the accepted Y prefix in its existing YQ
buffer. Normal Y generation, producer sends, owner receive/capping and the
`B_CAP=2^26-1` sink behavior remain charged. At the generation barrier, the
owner seals the common accepted prefix and cleans its eight 64-byte lines per
Y to the point of coherency, completing the clean before publication. The
link premise makes those published source bytes visible through the offered
aperture. All recipients see the same count and bytes.

The owner's Y processor first consumes its own local Normal buffer, then
acts as the owner sender and offers that immutable buffer sequentially to the
other 1023 members. Each recipient's Y handler runs its complete forward/inversion/
backward/keys chunks through the read-only YVIEW aperture before returning.
The completing store stalls through that work. No recipient main reads another
processor's live aperture. One additional sender processor per owner is
conservatively charged; the implemented yowner wrapper can reuse the owner's
Y processor for both local consumption and sending.

The handler uses a private Normal stack. Backward work may use SP as an aligned
read-only Y address, as in the original kernel, but its three entry scratch
stores and matching reloads move to six local Normal words beside the existing
fixup counter. Reconstructing that local address adds two instructions per
nonempty chunk. Forward and keys keep their ordinary Normal stacks. Table
addresses embedded in a Y frame denote the identical local table VAs at the
recipient, not owner-specific physical addresses.

The two source buffers preserve the original extra 2-GiB separation. One
32-GiB YVIEW at 38 TiB uses Device attributes and the EL1 read-only
PIR permission index, and is reused between offers; the descriptor carries the
selected source base. Both source bases and the aperture preserve the low
31-bit chunk alignment. Their payload capacity is exactly
`512*(B_CAP+1)=2^35` bytes. Local YQ reservations are retained conservatively
even at recipients; the gain does not require reclaiming their memory.

An empty prefix still produces all recipient completion offers. A nonempty
handler returns only after its last read and its release completion update.
No buffer reuse, next-generation publication or stop takes effect until the
active handlers and sender loop have drained. A downstream wait chain is
sender → Y leaf handler → relay → home; downstream handlers do not wait on the
sender or recipient main. No barrier is entered while holding such a chain.
The local completion epoch is bounded by the unchanged round count and uses
an isolated Normal flag. These rules establish lifetime and avoid a cyclic
wait; the finite model checks source publication for both buffers, empty work, and the
adversarial state between the last DONE store and handler return.

The current driver has `yclean`, `yoffer`, `yconsume`, `ycoin`, `ychunk` and
`yowner`. Its outer round protocol is:

1. Producer sends complete and all generator mains finish the generation
   barrier. Owner reception and saturation use the original ydlve/dispatch
   path; the accepted count is clamped to B_CAP.
2. Each owner main writes an eight-word Normal context and wakes its Y
   processor. yowner cleans the source, calls local yconsume to its return,
   then calls yoffer for all 1023 remote recipients to its return.
3. Only then yowner release-stores the epoch to DRAINED, at
   `34 TiB+0x7120800`. No source read follows that store. Every remote handler
   has returned and the local handler returned before the offer loop began.
4. All generator mains wait on their isolated DONE epoch. Owner mains also
   wait on DRAINED before entering the global two-wave barrier. Reuse, reset
   and the next generation start only after that barrier. The separate
   reset barrier ensures the next receive buffer is selected everywhere.

A recipient's DONE store occurs before its handler RET and is insufficient
on its own to prove that the sender has resumed. The DRAINED condition closes
that gap. An executed mutation which publishes it before the offer loop is
caught. The modeled DONE-only schedule admits premature reuse; the corrected
schedule does not. A main's wait is the shipped wcnt instruction loop under
the monitor premise; executed conditional-store schedules count 14 instructions
with an entry stale event and 24 when other writes alias its monitor granule.
The measured glue bounds are:

- `yclean`: `32*n+4` for a nonempty n-Y prefix, three instructions if empty;
- `yoffer`: `6*1023+3` instructions for the remote offer loop, including its
  completing stores and barriers;
- `yconsume`/`ychunk`: `11+63*c` glue instructions for c chunks, excluding
  the explicitly counted coin sampler and leaf/inversion callees;
- the two-attempt nonzero scalar sampler remains within 64 instructions per
  chunk, including the tested rejection/failure branches;
- the relocated backward scratch adds two instructions per nonempty chunk;
- yowner adds 20 instructions, excluding its three explicitly charged callees.

The producer/owner path costs at most
`30+3+23+21+32=109` per produced Y, covered by 129. The fixed cleaning return,
metadata, local owner consumption wrapper, dispatch, release/acquire,
completion and main wait are covered by 576 per recipient offer. Explicitly,
`6*1023+3 + (11+14)*1024 +20+4+14+4*128 = 32291`
bounds fixed work per sub-ring: offer sender, consume returns/main waits,
owner wrapper, clean setup, additional drained wait, and up to four MSG-priced
control exchanges for start/metadata/reset. This is below `576*1023`. Chunk
work and the round barriers are charged separately. The chunk
glue plus the two scratch instructions fits the retained 41 driver and 16
inversion-call instructions plus the additional 64 per chunk; coin and actual
inversion/leaf work are separate. The main-wait model charges at most 14
retired instructions per completion under the stated monitor premise; it does
not turn leaf execution into free waiting work. The eight additional
generator/relay barrier waves per round are also charged.

The charge `NODES*CHK` bounds the aggregate number of chunks in a round, with
`CHK=2`; it does not impose a two-chunk limit on an individual recipient.
For recipient counts n_i after capping,
`sum_i ceil(n_i/CH) <= NODES + sum_i n_i/CH` and
`sum_i n_i <= NODES*L_CHAIN`. Since `CH=2^22`, the aggregate is at most
`NODES*CHK`. A heavily loaded individual buffer can use more chunks, and
the driver loop and its count formula cover that case.

The driver state, completion flag and private stack are reserved at
34 TiB plus offsets `0x7100000`, `0x7120000` and `0x7140000`, respectively.
The offer port at 34 TiB plus `0xb800` lies in the existing rounded Device
page. The current registry extends the logical PY interval through the end of that
already mapped Device page.

The initial experiment also executes yconsume with the actual coin sampler,
forward, inversion, backward and keys routines on four full-field Ys. It checks
1448 EB words against the independent reference, emitted records, the static
pointers, the read-only aperture, stack restoration and completion. This
particular execution takes 46699 instructions; it is a finite composition check,
not a per-chunk production charge. Separate destructive-register callback tests
measure glue through the production maximum of sixteen chunks. Valid but wrong
EB pointers and seed links cause detected output mismatches.

### 4.4 Relay and home source-address schedules

The dedicated relay and home delivery processors use SP as a cursor into a
small immutable Normal table of source-aperture addresses. They do not execute
ordinary stack stores in those handlers. Each LD64B loads eight addresses into
x0..x7; payload loads and stores remain separate counted operations. SP and
x0 are restored before return. Partial messages retain the conservative
reference partial-copy envelope.

The relay table is 512 bytes at `2^48-1024`, within DPW, and costs 262
initialization instructions. A full 128-record input replaces 60 source-address
instructions by 17, saving 43. PW moves to x16 and the high completion count
to x8. A fault stub uses the current lane's metadata register as scratch only
after its record words have been stored through the alias. It preserves all
future source-address registers. The final partial message can add one reset
instruction; the split/control allowance covers it.

The home table is 1024 bytes at `2^55-4096`, within HCS, and costs 518
initialization instructions. A full 256-record block replaces 121 source-address
instructions by 33, saving 88. The current HH cold paths preserve x0..x7 and
SP. The last payload pair retains the reference omission of an unnecessary
R2-cursor increment. The primitive savings occur on full messages. The aggregate ordinary
coefficient also covers partial messages using the dominance calculation in
section 7.1. No Phase-0 source-scheduling credit is taken.

Thus ordinary routing receives exactly `43*6/128=129/64` relay and
`88*6/256=33/16` home instructions of credit per incidence. These are source
address savings; none of the payload transfer loads or stores is removed from
the ledger.

### 4.5 Direct pair heads and bounded clones

A home executes one serialized delivery handler. The insertion algorithm does
not claim multi-handler concurrency merely because its head exchange is SWP.
R2 consists of immutable completed 64-byte pairs of 32-byte records. A nonzero
64-byte-aligned pair address can serve directly as an HH head. Its bit 3 is
zero. Allocated 16-byte NX nodes are represented by their address plus eight,
so bit 3 distinguishes them from direct heads; a node stores its predecessor
representation and pair address. Zero remains the empty-head sentinel.

An empty insertion is BFI, SWP, CBNZ: three instructions. On a nonempty
insertion, the cold path examines the existing direct pair or tagged chain,
compares the relevant key words, and preserves the original candidate-cap
behavior. If lane 1 sees the current pair already installed by lane 0, it
examines the older lane only and does not add a self-match or duplicate chain
node. Induction on insertions preserves the set and chronological order of
eligible older records represented by each bucket. Matching and listed
candidates therefore retain the reference scalar metadata and verification.

Two boundary cases require explicit immutable clones:

1. A fresh singleton from a partial prefix is copied to lane 0 of its own
   64-byte clone. Lane 1 has a K3 differing by XOR 1, so it is outside the
   actual record's bucket. That clone is inserted only in the actual bucket
   and is never later filled with another record.
2. A complete pair at R2 address zero receives one whole-pair clone, avoiding
   collision with the empty-head sentinel.

The clone arena contains 44003 rows. A home has one parent relay and at most
one partial final per group. Charge up to two singleton clones against such
a final, including the leading singleton created when its old odd prefix is
later completed. The Phase-0 clearing pass still uses original R2 keys to clear
all affected buckets, then explicitly resets the NX and clone cursors. The
arena bound covers the Phase-0 phase and the full ordinary schedule after that
reset. NX allocates at most one node per indexed record and begins at
`NEXT_LO=2^64-24576 TiB`, the actual arena start, rather than the reference's
old affine decoding origin.

The executed insertion bound is

`3*new_records + 10*prior_visits + 5*K3_verifications + 6*successful_listings`.

The prior population is independently counted from the existing bucket
records, not inferred from the proposed cost expression. The composed empty
ordinary 256-record handler takes 1193 instructions, including 768 fast-head
instructions and one receive-completion store. Its remaining 424 instructions
fit the charged 425 copy/control allowance after the source-table credit and
the separately charged null check. The zero-address pair and singleton
wrappers are paid by the partial/clone envelope. Section 7.1 charges the first
full block and final partial block together in the current group, including
a fallback full block carrying an old partial prefix.

Only the ordinary PHI-weighted prior row uses ten instructions per visit.
All unweighted prior populations keep twelve, including same-X, defective
batches, weak-label defects, Phase-0 visits and exceptional deterministic
maxima. Verification remains five and listing remains six. The larger twelve
charge is intentional conservative work, not an assumed different algorithm.

## 5. Current layout and initialization

`map_check.py` gives the current half-open registry, symbolic translation-table
counts and explicit disjoint physical packing. It imports no candidate
certificate to decide its layout. Both VA roots stay within their 55-bit
ranges, and every role's Normal data, conservatively reserved Device extent,
alignment gaps and page tables fit below `2^56` physical bytes.

Hop-2 has `2^42` destinations, each with a 12288-byte VA block: two distinct RW
pages and one alias of the node's shared RO guard page. The upper range holds
2932031007402 blocks, and the lower holds 1466015503702. The lower range ends
at 32767 TiB, before the existing relay control regions. The two ranges use
the same per-destination LC2/RC2 tables. The RW/RW/RO pattern is relative
to each piece's first block; it is reset at the lower piece and is not computed
from the global VA page number modulo three. Initialization and finalization run
over both ranges while preserving their cursors; no split branch is added to
each record's hot path.

The relay request buffer starts at 256 TiB and retains its 65536-byte control
prefix. Its end is 18013298997788672, rounded down to a 64-KiB boundary, giving
`RELAY_CAP=277059750328320`. This leaves 57344 bytes before lower H2. Restoring
the old cap would overlap the new H2 layout and is rejected by the map check.
Both relay roots now use 4-KiB translation. All replacement relay page tables
are conservatively charged anew: 51809096722 tables at 1040 instructions each,
even though the symbolic construction requires fewer. No old relay setup
credit is reclaimed.

The YVIEW extension receives an allowance of 32900 new 4-KiB table pages at
1040 instructions each, or 34216000 instructions per generator. The current
map verifies the source buffers are contiguous in the backing construction
and preserve their actual separation. It reserves both the owner YQ data and
the receiver aperture's address-space extent. The new sender setup adds 4096
instructions per owner, and the Y leaf processor's VBAR adds eight per
generator.

The clone arena is
`[36028796804530176,36028796807346368)`, after the existing HCC contents and
before HCS. It holds `44003*64` bytes. A `2^23`-instruction allowance per home
covers its extra mapping/setup and cursor initialization. The source tables
fit inside the existing Normal DPW/HCS regions, away from all isolated
completion monitor granules. Construction of individual clone rows is charged
when they are used, not claimed as free preinitialized data.

Code is reserved and delivered at an allowance of `2^19` bytes per node.
The additional delivery charge over the reference is
`ALLN*((2^19-2^18)/32)*2*MSG`, with `MSG=128`. Existing code/state reservations
hold this image. The leaf object's own vbarg and shadow offsets remain
`0x5000`. The route object's Phase-0 vbarx and relay vbarr use `0x20000`,
including their matching site shadows. These processors set their respective
VBAR once: the main performs Phase 0, while the Y leaf processor performs
ordinary keys, including owner-local work. No per-chunk VBAR switch is assumed.

The initial experiment expands the exact current sources and checks 30320
leaf, 3690 field, 73111 route, 137 Y-driver and 50 pointer-initializer instruction
slots, including layout fill and route initialization helpers. Their total is
429232 bytes, below the charged 524288. LLVM21 llvm-mc assembles these streams
with `-triple=aarch64 -mattr=+cmpbr,+ls64,+lse,+rand,+rcpc-immo`; the sum of their
`.text*` sections equals that counted total. The interpreter places cold
`.text.hit` sections after the ordinary section and checks branch reaches in
that layout; the assembly export reproduces it. All current fault offsets are
encodable ADD immediates. The package's source/count and literal site checks,
not a private receipt, establish the quantities used in the work bound.

## 6. Success, relay-cap losses and continuation weights

The semantic transformations preserve each accepted leaf's mathematical key
and candidate metadata. Serial delivery and the explicit drain rules preserve
the reference's phase boundaries. The new data structure enumerates the same
eligible prior records; the larger hop-2 block changes batching but not the
per-home record cap or which key fields are compared. The only changed
success parameter is the smaller relay request capacity.

Use the reference notation in `base/success.md` section 6.3–6.4. Recompute

`rl_cap_cr = floor(9*RELAY_CAP/10)`,

`TH_RELAY = RELAY_CAP - 18*FIXED_LEAVES - rl_cap_cr = 24924152712192`,

`cn = (TH_RELAY-M_VG_MAX)^2`.

`TH_RELAY>M_VG_MAX` remains checked. In the relay-loss expression replace both
denominators: `d_rl=rl_cr/rl_cap_cr+v_rl/cn`. Also recompute the Phase-0
`P_MINE_CAPS` term with the current relay cap, then update `ck`, `glob` and the
full success expression. No cap term is reused merely because its displayed
decimal would round to the same value.

`costs.py` independently derives the exact cap-induced tail increment from
these reciprocal denominators. `rl_cr` and `v_rl` are affine in the round
prefix q, while the pair intensity is quadratic, so the increment is a
polynomial of degree at most three plus its constant Phase-0 change. Its
integer Newton coefficients are converted to an equivalent integer Horner
form and checked at separate q values. A second calculation evaluates the
unchanged full reference `spx` formula with only the new cap globals supplied;
the two calculations agree exactly on their checked inputs.

At every one of the 22000 group starts, let `delta(q)` be this nonnegative
exact tail increment and `PHI_old[g]` the certified reference continuation
upper bound. Set

`PHI[g] = min(1, ceil(2^256*(PHI_old[g]+delta(m_g)))/2^256)`.

This is an upper bound for the new continuation probability. It can exceed a
fresh direct rounding by at most one extra `2^-256` unit, which is retained in
the work. All weights are recomputed, lie in `[0,1]`, and are nonincreasing;
selected values are also compared directly with the full exact success
formula. Thus the reference tower-of-expectation argument is used with the
current unconditional continuation bound. No pointwise or Phase-0-conditioned
continuation probability is substituted for it.

At the full horizon the current upper bounds include:

| Term | Upper bound |
|---|---:|
| `P_MINE_CAPS` | `0.000000067042109627` |
| `d_cl` | `0.000003158491836518` |
| `d_rq` | `0.000002087211641220` |
| `d_rl` | `0.000013146065709312` |
| `d_p2` | `0.000000147900570669` |
| `rows` | `0.000320575459384161` |
| `glob` | `0.000001455486780113` |

The remaining terms and all probability inequalities are computed as exact
rationals, including the existing defect, failure and scalar-hiding terms.
Let `LO_RUN` be the resulting exact lower bound. Set

`K_GATE=ceil((39/100)*2^128/LO_RUN)` and `P_GATE=K_GATE/2^128`.

The certificate proves `P_GATE*LO_RUN>=39/100` and that one smaller gate
numerator fails. The old gate also fails under the new capacity: it is checked
as a negative control, rather than silently retained. The coordinator's gate
still costs 256 instructions and the total is `256+P_GATE*E_RUN`.

No new random-oracle or uniform-label heuristic is introduced by this update.
The inherited Chen–Stein, lattice, weak-label and defect arguments remain the
mathematical basis. Finite full-field kernel tests check the implementation at
their tested sizes; they do not establish production-scale collision laws.

## 7. Current cost derivation and complete ledger

Let `c` denote the reproduced reference certificate and `N=NODES`, `R=N_REL`,
`H=N_HOME`. For group g write `q=q_g`, `s=c.SG[g]` and
`gp=N*6*m_g*Q2`. There are `s/6` incidences. Keep the reference fixed work and
unweighted populations unless a replacement below explicitly changes them.

### 7.1 Ordinary leaf and route work

The current leaf/Y row is `c.WGL[g]+delta_leaf`, where

```
delta_leaf =
  N*(ceil(q*L_CHAIN*129/M_G)-ceil(q*L_CHAIN*c.YRING_TOT/M_G))
  + q*(576*N*1023/1024 + 8*4096*(N+R))
  + N*q*CHK*64
  + N*ceil(q/M_G)*L_CHAIN*3
  - 19*s/36 - s/6.
```

Here `c.YRING_TOT=54296`; the first line replaces the repeated ring copies by
owner reception/publication. The next lines charge offers/barriers, chunk
glue and producer seed preparation. The last line is the static EB and linked
seed savings. Exact ceilings are retained. All ordinary generation and all
actual arithmetic, verification and key paths remain in their reference rows.

Before source-table credits, use

`K2B=256`, `FLUSH2=21`, `DLVHB=512`, `DLVHS=811`,
`BLOCKH=FLUSH2+EXC_ENTRY+DONE_C+DLVHB=565`, `FINALH=1153`.

The extra FLUSH2 instruction covers the K256 guard displacement. The old
single SUB is replaced by SUB 8192 followed by ADD 32. The current route row is

```
c.WGR[g]
 + s*565/256 - s*c.BLOCKH/c.K2B
 - (s/6)*(129/64+33/16)
 + H*(1153-c.FINALH)
 + 64*R + ceil(s/256) + 2048*H.
```

The `64*R` covers split initialization/finalization and final completion
control per group. `ceil(s/256)` covers one full-home null-pointer check.
`2048*H` covers partial-pair fallback, singleton/zero clones, return control
and their bounded frequency. The primitive source savings apply to full ordinary blocks. Their rates are
subtracted from an envelope charged on all records, so partial-message
coverage needs an explicit argument. At hop 1, the combined full-block
allowance is `(FLUSH+EXC_ENTRY+DONE_C1+H2S_FULL)/128=55/128`
above the ordinary per-record work. After subtracting `43/128` for the relay
source table, `12/128` remains nonnegative. A partial message pays neither the
full-block fault/flush nor the full h2s dispatch; its fixed partial work remains
in FINAL. Thus the new envelope still dominates partial arrivals even though
their source-address instructions themselves have not changed. The one final
completion adjustment per relay is paid by the fixed split/control allowance.
For homes, let `o` be the old prefix (`0<=o<256`) and `m` the new records
received by a home in this group. There are `f=floor((o+m)/256)` full sends.
Reserve the first full send, if present, as a boundary event. The remaining
`max(0,f-1)<=m/256` sends have no old prefix and use the ordinary full path.
Their copy, control and transport fit `A=565-88+1=478` instructions each,
including the separately charged null test. Thus the all-record coefficient
pays these sends even when some of the records end in a partial block.

The first full send and the final partial are both charged in this group.
Their copies and clone wrappers are bounded directly from the current source.
After removing three fast-head instructions per fresh record, measured cold
chain work, and the completion store already in DONE_C, the primitive costs
are `DLVS=4`, `DLVP=4`, `DLVHS=12`, `DLVHP=6` (10 at address zero).
The fallback has 16 fixed instructions, reduced by one for a nonzero prefix.
A zero-length segment costs one; a nonempty segment starting in lane `a` with
`q` pairs and final singleton `b` costs
`5+a*(single+1)+q*(pair+2)+b*single`. The check binds this expression to the
actual SEGL source, executes the five primitives, and compares it to every
executed full-handler case. Enumerating all 33152 non-fast prefix/end states
gives a maximum of 1072, attained by an actual prefix1/end255 execution.
The inherited dsg expression is one instruction low on an odd-start singleton;
both its old literal and formula maxima are still 811. The current charge
uses the newly measured bound, rounded up to 1200 per boundary handler.

First-full transport costs `FLUSH2+EXC_ENTRY+DONE_C=21+24+8=53`.
Final-partial transport costs `FINH_DEST+FINH_SENT+DONE_C+MSG_LINK=16+5+8+14=43`.
The complete boundary charge, including six extra zero-pair repair instructions
and two cursor resets conservatively, is
`53+1200+43+1200+6+2=2504 <= FINALH+2048=3201`.
The measured version is 2248. Copy maxima already contain clone wrappers;
this additional zero allowance is conservative. Hash insertion, prior visits,
key verification and candidate listings remain in their separate rows.
`home_checks.py` checks this bound and the aggregate inequality for every old
prefix and the first three block lengths; adding 256 new records adds at most
one ordinary full send and exactly 478 to its record allowance, proving the
remaining lengths. No charge is borrowed from an earlier PHI-weighted group.

The ordinary home row is

`3*s + ALLN*STOP_CHECK + 10*main_g*PAIRM`,

where `main_g=s*gp+s*(s-1)/2-c.SAMEX[g]`. Its verification row remains
`5*main_g*PAIRK2`. The key-path row is unchanged per leaf. All these ordinary
rows, including the unchanged exchange row, use the new PHI weights.

### 7.2 Phase 0, exceptional paths and unweighted visits

The Phase-0 update subtracts one instruction per newly inserted record,
increases the old visit coefficient from four to twelve, substitutes the
K256 block/final costs, adds the full-block null check, and adds `3072*H`
fixed/partial/clone control. It deliberately takes no Phase-0 source-table
credit. Its defective-event maximum is recomputed as

`DEF_MINE=ceil(c._P_MINE_DEF*(c.W_MINE_DET_MAX+8*c.S_MINE_TOT*c.CAP_H))`

`=126022896`.

For the exceptional post-Phase-0 event, retain the entire old deterministic
core, ignore its reductions, add every positive new ordinary operation, and
increase its visit maximum by `8*SG_FULL*CAP_H` before multiplying by the
unchanged exceptional-event probability and rounding upward. This gives
`EXC_WORK=6707823107677`. It is an upper bound on the changed implementation,
not a sampled average on exceptional inputs.

The old same-X visit contribution excluding EXC_WORK and the old DEFV term
are each multiplied by three, changing four to twelve. DEFWLT receives eight
times each of its prior-visit populations, including its DWS term. Other
verification, listing, resolve, request and key-path coefficients are
unchanged. These populations remain unweighted, exactly as required by the
reference coupling argument. The proof of the ten-instruction ordinary path
does not authorize reducing these conservative twelve-instruction charges.

### 7.3 Cold additions

Row 1 retains the full original `ALLN*NODE_INIT` and adds:

| Added operation | Charged instructions |
|---|---:|
| Extra code delivery | `ALLN*2^21` |
| YVIEW tables | `N*34216000` |
| Additional owner sender setup | `(N/1024)*4096` |
| Y leaf processor VBAR | `N*8` |
| Clone arena mapping/setup | `H*2^23` |
| Replacement relay page tables | `R*1040*51809096722` |
| Relay source table | `R*262` |
| Home source table | `H*518` |

Row 7 adds `N*6*CH*9` for static EB pointers. Row 8 adds
`N*15*2^43` for linked table seeds. Row 10 adds `128*R` for fixed split/control
work including Phase 0. These additions are paid once when the run opens.
The enlarged code delivery is per node, sharing that node's code memory among
its dedicated processors; all processors' executed instructions are counted.

### 7.4 Complete rows

The following are upward-rounded base-2 instruction logarithms. `costs.ROWS`
contains the exact rational or integer values and checks their sum is E_RUN.
P denotes once-per-open preprocessing, E the new PHI-weighted rows, and W the
unweighted rows. Row 19 contains both kinds explicitly.

| # | Work | Kind | log2 upper |
|---:|---|:---:|---:|
|1|Node setup and all cold additions above|P|112.2078|
|2|`R_g=[b_g]Q`|P|85.2999|
|3|Group messages|W|104.4254|
|4|Phase-0 mining, routing, joins and defects|P|123.0426|
|5|Phase 1–3 coins|P|100.1700|
|6|X descriptors|P|85.5850|
|7|Entry X/header and static EB initialization|P|95.3399|
|8|P1, forward and linked K' tables|P|123.9321|
|9|Beta tables|P|117.3576|
|10|LC/RC initialization and split control|P|108.0001|
|11|Home zeroing|P|125.0224|
|12|Home initial state|P|86.3220|
|13|Leaf work and current Y transport|E|134.3559|
|14|Metered arithmetic fixups|W|107.5871|
|15|Key paths|E|127.2365|
|16|Routing, finals, partial/clone control|E|130.7419|
|17|Ordinary home insertion and prior visits|E|128.8066|
|18|Same-X visits and exceptional work|W|42.6324|
|19|Secondary-key verification|E/W|112.3500|
|20|Listings|W|94.6804|
|21|Defect visits, verification, lists and key paths|W|118.1483|
|22|Candidate resolves|W|113.8741|
|23|Requests|W|112.3976|
|24|Exchanges|E|105.4873|
|25|Inherited global barriers|W|112.6732|
|26|Phase control|W|106.1257|
|27|Output|W|20.0000|

The extra generator/relay round barriers are already in row 13, not omitted
because row 25 retains its original name. Unchanged per-event rows 15, 19 and
24 have updated values because the continuation probabilities changed.
`E_PRE` includes the current once-per-open additions. The final expected cost
is `E_TOTAL=256+P_GATE*sum(ROWS)` and the score is computed from `E_TOTAL/1000`.

## 8. Checks, limitations and provenance

The top-level experiment must run the current cost calculation and current
source checks together. The mathematical certificate retains all original
checks and independently recomputes the relay-cap delta, all continuation
weights and exact work. The current map check symbolically counts descriptor
tables and constructs the physical packing; edge translations check the
split layout and shared guard permissions. These are finite programs for
symbolic formulas, not an enumeration of the production memory.

Current kernel checks compare the composed leaf fractions and records against
the reference mathematical model, build their tested table rows with the
actual changed assembly, check read-only Y access and register/stack
restoration, and count executed instructions. Routing checks cover full,
partial, empty, same-pair, clone, cap, sent-fault and dropped-fault cases,
including both sides of the H2 split. The comparison is by complete candidate
metadata and normalized bucket contents when clone addresses differ.

The initial suite includes controls that exercise the dependencies being
claimed: valid but wrong EB pointers and seed links cause numerical mismatches;
missing cache cleaning gives a wrong modeled payload and reading an expired
aperture raises a lifetime error; a wrong guard displacement corrupts routing;
mutated clone/self-head logic changes candidate enumeration; the old relay cap
overlaps the new layout; and the old gate misses the success floor. The wait
bound is conditional on isolated monitor granules and nonspurious WFE behavior;
the map checks that the chosen completion words are isolated.
Source assembly verifies the actual instruction forms and branch reaches.
An interpreter's ability to parse an instruction alone is not encoding evidence.

The experiments are evidence on their tested populations. The deterministic
invariants and charged envelopes extend those operations to the production
counts; the probability proof remains the preserved full-size mathematical
argument with the exact cap substitutions above. The architecture/link/wait
premises in section 3 remain explicit limits of the claim. No external private
receipt, benchmark timing, paid review or previous verdict is required by any
step of this argument.

The direct predecessor is bradcush's `4ce2d07e`. Its Y-sliced sub-rings,
Eisenstein K' labels, guarded and bounded-unguarded arithmetic, full-size
success proof, rendezvous model and permission-fault blocks are retained and
credited. Its preserved provenance credits the underlying design to soispoke's
`6e6d0e4f`, earlier sub-ring work to soispoke's `3ff10cb7`, the gate retune to
aburan28's `3d0fb02b`, and the detailed kernel and mathematical contributions
listed in `base/proof.md` section 7. In particular the `17/32`, 22000-group
schedule comes from soispoke's `eed62949`; OFF5PAIR is credited there to
Jaydearcadian's `d719e7cb`, and the predecessor's hop-2 subtraction fold to
aburan28's `5f726b85`. Those credits, including kaankolcu's DT/DF idea and the
Galbraith–Wang–Zhang and Arratia–Goldstein–Gordon references, remain intact in
the verbatim reference files.

This update supplies the static EB pointer schedule, linked seed loads,
immutable-buffer handler transport, source-address schedules, direct pair
heads with boundary clones, split K256 layout, and their revised accounting.
It derives the cost of each replacement here rather than treating a borrowed
mechanism or an earlier acceptance as a free operation.

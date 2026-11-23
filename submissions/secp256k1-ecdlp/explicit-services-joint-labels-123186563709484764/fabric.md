# Card54 memory fabric construction and remaining verification

This is the construction to bind to the current scalar services. It replaces
the completing-store service rejected by audit1. The concrete saving remains
removal of scalar payload copies and address operations. Every payload access,
field configuration, software dispatch, command, acknowledgement and wait is
executed by a counted AArch64 instruction. No memory write waits for software
RET. Its bus response acknowledges only that memory/register transaction.

## Local address spaces and shared routing

Each CPU has at most56 physical-address bits, as independently packed by the
map checker. A CPU transaction enters one local decode endpoint. The decode
range and a fixed port class select either local RAM/registers or the shared
fabric. A fabric packet has a56-bit destination-local address, a fixed class,
and an84-bit destination tag consisting of two role bits plus an82-bit node
index, zero-extended for smaller roles. Role codes are G=0,R=1,H=2,B=3.
There is one additional control-only G endpoint numbered2^66; it performs
no arithmetic-chain work. It receives a conservatively reserved full G map
and endpoint/state allowance. Arithmetic NODES remains2^66. Valid G indices
are0..2^66 inclusive; other unused indices are decode errors.
These sideband fields are wires/register outputs of the memory interface;
they are not a wide scalar CPU address or a claimed arbitrary-width ALU.
Two-word variable destination identifiers are written by scalar instructions
to the explicitly charged route registers before sending. Fixed source
identifiers are endpoint constants. Ordinary AXI transaction ID fields serve
response tracking; they are not assumed to preserve a source tag.

Use a fan-in tree of16-input AXI multiplexers, then a fan-out tree of16-output
AXI demultiplexers. The latter selects successive four-bit destination-tag
slices, most significant first. At most21 levels suffice for84 bits. This is
one shared network, not a separately allocated path for every possible pair.
Absent role/index leaves are pruned. The normal tree induction is exact:
after j fan-out levels the selected prefix equals the top4j tag bits; at21
the selected leaf equals the tag. The mux's request-ID prefix tracks each
source port; its response logic returns R/B responses through the same
prefixes. Full-width tags are not looked up in a table of all source-target
pairs. The unused top indices of each role are decode errors.

The selected construction uses the PULP AXI components at commit
ee000c0eba9e1c2b6ffdab163a2957419d143461. `axi_mux` prefixes source-port bits
to transaction IDs and keeps a one-entry write-selection FIFO.
`axi_demux` accepts separate AW/AR selection inputs, keeps the AW selection
for W routing and routes B/R responses. Its documentation requires selection
to remain stable while an AW/AR handshake is pending. Instantiate MaxTrans=1,
MaxWTrans=1, one locally outstanding transaction per direction, and
UniqueIds=1: every active ID has only one outstanding transaction in that
direction. A source may retire more stores into its CPU store buffer, but
the bridge drains them sequentially. This can reduce speed without adding
retired instructions. It is not a claim of high throughput.

Keep512 data bits,56 address bits and an AXI ID width at most88 bits
(four local initiator bits plus84 routing-prefix bits). Include destination,
source, local source/destination PE and class in at most192 USER bits.
The two84-bit endpoint tags, two4-bit PE fields and8-bit class total184 bits. At a switch each buffered channel
is at most1024 bits including payload and controls. For an explicit upper
census, allow16 ports*5 channels*2 slots*1024 bits,16 one-entry selection
FIFOs of128 bits, and4096 additional bits for arbitration/counters. This is
21248 bytes, below the reserved64KiB per switch. No 2^ID counter array is
used: that would violate the stated memory bound; the UniqueIds case is
essential and its prerequisite is enforced at the single-outstanding bridge.

For E=2^82+2*2^66+2^40+1 physical endpoints, a full padded84-bit tree has2^84
leaves. Its two global trees use2*(2^84-1)/15 internal switches, which is
less than E. Allow up to two local mux/demux switches per endpoint for at most
16 dedicated processors/interfaces. Thus four64KiB switches per endpoint
are conservative. Reserve a further1MiB per endpoint for bridge/control
state; this includes the dedicated processors' small stacks, descriptor and
one-transaction buffers. The map/cold census charges this reservation in
addition to all mapped data, page tables and code. This is a circuit-state
construction bound, not a synthesized area or timing measurement.

## Endpoint decode, without hidden services

G-to-R: the data port carries the40-bit relay destination obtained from the
translated compressed H1 range. For the read/write page of block d, the
physical displacement is4096*d+offset,0<=offset<4096. The endpoint attaches
the fixed source member m in0..1023 and targets R's preallocated slot
8192*m+offset. The header occupies0..31; the records begin at32. The last
record's two stores fault at source offset4096 and are reissued through the
separate compact32-byte tail alias. That alias targets slot+4096..4127.
The receiver slot is8KiB; its final unused bytes are never copied as data.
Destination extraction uses fixed bit slices after selecting the mapped
range. H1 lower/upper runs have fixed decoded offsets; no software-selected
source window exists. G READY8TiB+8*d targets that R's record processor
READY. R return-ACK1TiB+8*m forms G id=(coarse_subgroup<<10)|m. The coarse
56-bit register is configured by the measured8-instruction loop per R per
subgroup and remains fixed until every member's fin1 and trailing DSB drain.

G return processor: the main processor is PE0 and the member-offer consumer
is PE9. Both use the scalar bulkACK field at24TiB. A fixed endpoint shim on
each local initiator interface drives USER.srcPE from that interface's
four-bit port constant before the global tree. This is explicit wiring;
it does not infer a PE from an AXI ID after the transaction has completed.
At successful RREADY0x8000 commit, R captures srcPE and sets a valid bit in
its record-return latch. Capture precedes the Rrecord publication event.
RREADY clear preserves this latch. Rcoarse0x8018 is a distinct register and
cannot alter it. The packed latch is exposed for cold reset at0x8020;
`rreturn_pe_boot` performs MOV/STR/DSB/RET, at most8 instructions including
caller/entry reserve. Reserve16 bytes per R in addition to the endpoint
state allowance and8 initialization instructions per R.

A nonzero RREADY commit while valid returns SLVERR, leaves READY and latch
unchanged and emits no event. An Rrecord G-ACK write without valid likewise
returns SLVERR and emits no event. Normal Rrecord ACK uses the existing
member/coarse G destination and attaches the retained PE as USER.dstPE.
After the G ACK register write commits, only that destination PE receives
an event; then R clears valid. No payload copying or scalar service occurs
in this register circuit. The shared ACK level is cleared by the actual G
source before every RREADY. Both xemit's fault path and fin1 retain that
counted clear and drain. Thus a previous PE's ACK cannot satisfy a new wait.

Only one G arithmetic source is active globally. The controller waits for
each mining/producer/owner worktoken's DONE; each owner's member offer is
send-wait-advance, and ACK is after yconsume RET. Gfin1 reaches Rrecord by
the same RREADY class; the Rmain finalizer exchanges only with H and sends
no G ACK. No next G can overwrite the R latch during a nested H wait.
The finite source test exercises PE0 mining, PE9 consumption and PE0 fin1,
with delayed nested H and actual shared-ACK clears. Mutants remove that
clear, return to the other valid PE, clear the latch with READY, or overwrite
a busy latch. This is construction/model evidence, not hardware evidence.

R-to-H: each R has a fixed40-bit parent identity; its local42-bit d2 picks
the global home id=(parent<<42)|d2. The two writable pages per12KiB H2
block become a compressed8192*d2+offset physical data port. The selected H
endpoint substitutes `{BASE[55:13],offset[12:0]}`. BASE is a local writable
MMIO field. The H physical map packs R2 in upper-then-lower logical order,
with an aligned sink at LIMIT. H's scalar dispatcher adds8192, clamps to
LIMIT, writes BASE and executes DSB before publishing ACK. Payload stores
are already complete before READY; H never changes BASE with a payload
transaction pending. No payload load/store is performed by the remapper.
The H READY aperture is8B per d2. H's fixed parent-ACK aliases are distinct:
CTRL+256 targets R record-processor ACK0x8008; CTRL+264 targets R main's
finalization ACK0x8010. The source final and ordinary full tokens select
their measured separate scalar tails, not an uncounted hardware handler.

Y archive and replay use fixed address classes. G source control is
34TiB+0x7200000: route(owner56,bank10) at+0/+8, ACK+16, READY+24 and512-byte
outgoing aperture+4096. Archive payload goes to B state+0x2000000 and READY
to BCTRL+24, where BCTRL=2^54+1GiB. B endpoint index is(owner<<10)|bank.
The archive receiver is B PE1. Its explicit reverse words at BCTRL+32/+40
and ACK store+48 target the producer G source ACK and PE0. ErrorACK2 also
publishes a nonzero event. No host completion replaces the receiver copy.

B replay uses BCTRL's same outgoing aperture with fixed replay class and
route(owner,0). It targets G endpoint(owner<<10), inbox34TiB+0x7510000 and
owner CTRL34TiB+0x7500000 READY+24, PE8. The owner receiver's reverse words
and ACK+48 target BCTRL+16, B PE0. All448 payload bytes are copied by the
executed source and destination instructions. Archive DATA and owner buffer
A are Normal Inner/Outer Non-cacheable. G command, payload and offer control
pages are Device-nGnRnE; the LS64 target support is the common chosen bridge.

Specialized commands use G source CTRL34TiB+0x7600000, source aperture+4096,
ACK+16 and READY+24. Their two-word84-bit route selects only a valid G or B
endpoint. The128-byte command inbox is G34TiB+0x7610000 or B2^54+0x2100000.
Receiver CTRL is source GCTRL+512, or BCTRL+512. Receiver READY+24 pulses
G PE7 or B PE0. Receiver reverse route+32/+40 and ACK+48 return to the
sending G's command ACK+16 and PE0. G receiver publication+80 is the main
worktoken; values1..5, including even2/4, pulse PE0 after commit. Its
arguments and return identity are ordinary scalar writes followed by DSB.
CN handling, bank begin/seal/release, replay and DONE are executed software
commands, with the full fixed and per-payload costs in their source census.

G offers use CTRL34TiB+0x7300000: a scalar member10 at+0 selects the fixed
owner subgroup, aperture+4096 transfers count/epoch, READY+24 pulses member
PE9 and member ACK+48 targets owner's ACK+16, PE0. Reads of the receiving
aperture observe that member's local incoming words. A member's read-only
32GiB view at38TiB permanently targets its owner buffer40TiB; owner endpoint
is its Gid with low10 bits cleared. Offset0..2^35-1 is preserved. Writes are
permission errors. The offer ACK follows actual yconsume RET, and the next
member/round cannot reuse the buffer before ACK. This static alias performs
no CPU service on a read. Frame table pointers remain recipient-local VAs.

Generic bases are G34TiB+0x7800000, R192TiB+20MiB, H2^55-32MiB and
B2^54+128MiB. At each base, the first three64KiB pages are Device control,
REMOTE and ACK aliases; the next two are Normal inbox and handler state.
Control has16 slots of64 bytes. A fixed local initiator index selects its
own slot's84-bit route. REMOTE+0..79 writes the selected endpoint's inbox;
REMOTE+256 publishes READY to its slot15+48 and pulses dedicated PE15.
The explicit generic receiver saves the source tag and ACK selector before
calling its copy/queue handler, then uses those private saved words for its
slot15 reply route. It does not reread the shared inbox after publication. ACK+8*selector targets that endpoint's slot(selector)+16 and
pulses exactly selector,0..15. Generic receive and ACK are distinct events;
no bulk or Y service processor is awakened by PE15's READY. Dynamic routes
are validated by the role/index decoder, including the single controller
Gid2^66; other unused high G indices are errors. Generic ACK follows only
the chosen copy/queue handler. Application work follows the explicit phase
schedule in callbacks.md. G command1 retains its production request queue;
collector G command2 has its own1024-ticket queue and counted initializer.
Response IDs return through
the independent source-port mux path, never through the destination tag.

There are16 reserved interfaces/PE slots per endpoint. G original role slots
0..6 remain available; main0, Ycommand7, ownerpayload8, memberoffer9 and
generic15 are explicit. Rmain0/record1, Hmain0/record1 and Bcommand/replay0
and payload1 are distinct; generic15 is reserved on every role. The map
charges all16 slots even when inactive:256 setup instructions/PE, measured
generic control/handler initialization199 plus57 caller/argument reserve
per endpoint,128/G for Y boot and per-G mining callback binding,64/B,
and16 for the sole controller bindings. These include scalar argument
materialization and calls; the48/50G and14B body counts alone are not the
charge.
The fixed per-PE setup reserve covers the local context/entry and stack
setup; it gives no permission to add an unbounded scheduler or polling loop.
Program bytes are shared read-only within each endpoint, not recopied for
every PE. Separate per-PE stacks are in the already charged endpoint state.

Phase0 adds at most2048*NODES+64*N_CLUST+64 token/sweep instructions to cold
work, excluding already charged mining, finals and coarse updates. The local
post-HHclear H reset is10 instructions plus argument/caller glue;64/H is
charged. It resets physical BASE/LIMIT, logical cursor/length and real E2
only after all Phase0 traffic drains. These rows and return-latch startup
and64*N_POOLS+8*(NODES+1)+128*NODES for collector initialization/selection
and the two drained Phase0 request resets, plus1024*(NODES+1)+2048*N_POOLS
for the selected-request application sweep
are added explicitly to the previous cold census; the combined amount must
remain below the retained2^112 charge. No original1352/G-round setup is
removed. Only the retained8192/G-round old barrier allowance is reused for
actual CN/producer/token control, under the checked population inequalities.

All mapped aperture PA extents are included conservatively in the per-node
packing, even where an extent describes a remote port rather than local RAM.
These extents need not be allocated again as stored copies. The receiver RAM
and archive buffers are separately counted. The final map checker must
reconstruct endpoint destinations from the actual physical runs, including
both H2 runs; a VA-only arithmetic test would not establish this decode.

## Registers, ordering and event input

The same pinned repository's `axi_lite_regs` exposes `reg_q_o` and
`wr_active_o`. A successful valid full-width write sets each enabled byte in
reg_d, asserts reg_update and returns OKAY. The register flops commit on
the clock edge; its B-response spill register has Bypass=0. No posted bridge
acknowledges these register writes early. Arm DDI0487M.d B2.10.2 defines
Device-nGnRnE as non-Gathering, non-Reordering, No Early Write Acknowledgement.
B2.4.1's endpoint-completion condition and DSB SY complete the write at
the endpoint. Thus BASE's value is installed before the following ACK;
a second readback is unnecessary for this correctly configured peripheral.
Malformed addresses/hardware faults are outside the same functional-hardware
premise used for ordinary RAM and arithmetic, not treated as algorithmic luck.

`axi_modify_address` assigns AW/AR address to the external mst_aw_addr_i /
mst_ar_addr_i and passes request data and responses unchanged. Its explicit
requirement is stable substituted addresses while handshakes are pending.
The source DSB before READY and receiver DSB before ACK establish exactly
that quiescence for BASE and route registers. Merely writing the registers
without those barriers is not the construction.

READY and ACK are scalar register fields. For each dedicated processor, its
event input is a single registered pulse from a successful nonzero write
to only that processor's designated mailbox/ACK fields. Worktoken values2/4
and payload errorACK2 are nonzero publications too. Coarse, route, BASE and
return-latch reset fields never pulse. Rrecord ACK uses its retained caller
PE, and generic ACK uses its explicit validated selector. A write of0 clears
the field without an event. The same clock edge commits the field and sets
event_q; the event is seen on the following cycle. A held-high READY does
not hold EVENTI high. A source load observes the field itself, not the event.
There are no timer, interrupt, SEV or unrelated monitor events in this chosen
machine. Each wait may begin with one stale EventRegister bit; the event
model additionally requires that WFE not return spuriously.

DDI0487M.d D1.7.1 defines the one-bit EventRegister, its clearing by WFE,
and implementation-defined wake events; it also permits spurious WFE returns.
The nonspurious choice is therefore a disclosed platform condition, not a
universal Arm property. Cortex-A76's per-core EVENTI documents an example of
the event-input interface; it is not evidence that A76 implements the other
required ISA features. The proposed processor is the same hypothetical
full-feature A-profile processor with this explicitly specified event input.

For a wait with one publication channel, `LDR;CBNZ;WFE;B` takes at most two failed four-instruction
passes plus the successful two-instruction check:10 instructions. The first
failed pass may consume the entry-stale event, the other consumes the single
publication event. Publication commits before its event; no third failed pass
can be caused by that event. Waiting latency is unrestricted. Every shared
READY/ACK case must be assigned to the correct processor; an idle processor
cannot receive uncharged wakes intended for another processor. A worktoken
call has two possible publications on main PE0: commandACK and DONE. DONE
can arrive first. Its command wait therefore allows14 instructions (one
initial stale bit, two publication events, then the successful test). The
actual ytoken_call/ycmd_send check exercises this order. The extra4 is
included conservatively inside512/command and1536/token; it is not assigned
the single-channel10 bound. Arithmetic bulk, archive and offer waits retain
their phase-serialized single publication channel.

DDI0487M.d C3.2.8 permits64-byte aligned LD64B/ST64B on Normal Inner/Outer
Non-cacheable and all Device types, subject to target support. The chosen
bridge exposes512-bit aligned memory transactions to RAM, with per-line
arbitration for the required single-copy atomicity. A naturally aligned
64-byte read/write is one beat; it does not invoke a copying engine. Word
and pair accesses use their actual byte strobes. The bridge must support
LS64; AXI512 width alone is not a proof of that support. No remote LSE or
exclusive-monitor invalidation is assumed. HH/NX atomics remain node-local.

The shared circuit carries memory transactions only. Its responses can
complete while a CPU is idle; the CPU must execute its own READY test,
dispatch, handler, return and ACK. Nested G-to-R-to-H software calls do not
hold outstanding bus requests: every caller drains stores before waiting.
Thus the shared root is released before a downstream handler sends. This
removes a possible cyclic hold-and-wait introduced by treating software RET
as bus completion. Main generation/consume/token control must preserve the
one-active-top-level-sender schedule; that is a separately counted driver.

## Evidence obligations

The in-package scalar source/checkers are to establish instructions and
payload/lifetime behavior. The finite fabric checker is to establish bit
decoding, retained AW/W/response selections, non-aliasing of role endpoints,
BASE stability, endpoint commit-before-event and tree/state arithmetic.
Controls must perturb valid destinations, event timing and pending writes
and produce wrong payload/return behavior. Such checks are not RTL synthesis,
an Arm execution result, or a proof of the production probability estimate.
Primary excerpts and this wiring/construction are required in the package;
local receipts or URLs alone cannot supply the rejected service premise.

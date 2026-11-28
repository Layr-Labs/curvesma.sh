# Memory Fabric and Explicit Scalar Services (Flat 64-Bit Analytical Machine)

This document specifies the memory-fabric construction and explicit scalar service
bindings used by 8cc64cb7's flat64 stored-leaf search (`flat64-5win-eisenstein`), extended here by shared table banks (role `T`, Section 2.4). It adapts the explicit scalar push/ACK
fabric of `283cdacf` (`soispoke`) to the uniform 64-bit flat analytical
address space (`[0, 2^64)` per processor, matching the disclosed convention of submission `4e2d5241`)
with `NODES = 2^58` arithmetic generators and explicit `TBZ`/`TBNZ` block-end
checks (no VMSAv9-128 stage-1 permission faults or tail aliases).

Every payload transfer, route-register update, software dispatch, command,
acknowledgement, and wait is executed by a counted scalar AArch64 instruction.
No memory write waits for software `RET`; bus responses acknowledge only the
underlying memory or register transaction.

## 1. Local 64-Bit Flat Address Spaces and Shared Fabric Routing

Each processor accesses a uniform 64-bit flat analytical local address space
`[0, 2^64)` (disclosed as an analytical machine convention in `proof.md` Section 2,
with only explicitly backed regions counted in `PEAK < 2^138` bytes). A CPU
transaction enters one local decode endpoint whose fixed range decode selects
either local RAM/registers or the shared fabric.

A fabric packet carries:
- a 64-bit destination-local address (`ADDR_WIDTH = 64`),
- a fixed port class, and
- an 84-bit destination tag `(role[1:0], index[81:0])` with role codes
  `G = 0`, `R = 1`, `H = 2`, `B = 3`, zero-extended for smaller roles.

There are `NODES = 2^58` arithmetic generators (`0 .. 2^58 - 1`), one additional
control-only `G` endpoint numbered `2^58` (carrying a full reserved `G` state
allowance and zero arithmetic chains), `N_REL = 2^40` relays (`0 .. 2^40 - 1`),
`N_HOME = 2^82` homes (`0 .. 2^82 - 1`), `2^58` archive banks B, and `N_T = 2,216,615,441,596,416` table banks (Section 2.4).
A B index is `(owner48 << 10) | bank10 < 2^58`; a table bank is role code `B` with index `2^58 + (pool << 42 | FAM << 36 | bank)`,
disjoint from it. Unused indices are decode errors.
Sideband route tags are fixed endpoint wires/register outputs, not arbitrary-width
scalar CPU ALUs. Variable destination identifiers are written by counted scalar
stores to local route registers before sending.

Routing uses a shared fan-in tree of 16-input AXI multiplexers (`axi_mux`)
followed by a fan-out tree of 16-output AXI demultiplexers (`axi_demux`, PULP AXI
commit `ee000c0eba9e1c2b6ffdab163a2957419d143461`), decoding successive 4-bit
slices of the 84-bit destination tag over 21 levels. Instantiating
`MaxTrans = 1`, `MaxWTrans = 1`, `UniqueIds = 1` enforces one locally outstanding
transaction per direction per bridge so no `2^ID` counter array is required.

Channel bit budget at each switch: `DATA_WIDTH = 512`, `ADDR_WIDTH = 64`,
`ID_WIDTH <= 88` (4 local initiator bits + 84 routing-prefix bits),
`USER_WIDTH <= 192` (two 84-bit tags + two 4-bit PE IDs + 8-bit class = 184 bits).
Every buffered channel fits comfortably inside 1,024 bits. With 16 ports, 5
channels, 2 slots, 16 one-entry selection FIFOs of 128 bits, and 4,096 bits of
arbitration state, per-switch state is at most `21,248 B < 64 KiB`. For
`E = 2^82 + 2*2^58 + 2^40 + 1 + N_T` physical endpoints (`N_T < 2^51`), a full padded 84-bit tree has
`2*(2^84 - 1)/15 < E` internal switches; reserving 4 switches of `64 KiB` plus
`1 MiB` of bridge/control state per endpoint is strictly conservative and is
included in `FABRIC_BYTES` inside `PEAK_BYTES < 2^138 B`.

## 2. Endpoint Decode and Contiguous Block Windows (No Faults, No Tail Aliases)

### 2.1 G-to-R (`H1` contiguous 64-KiB stride, `K1 = 1024`)

On generator `G`, the contiguous hop-1 block window `H1` occupies
`[64 PiB, 128 PiB) = [0x0100_0000_0000_0000, 0x0200_0000_0000_0000)` (`2^40`
relays x `64 KiB` stride). For relay `d1 < 2^40`, block `d1` occupies
`[H1_BASE + (d1 << 16), H1_BASE + (d1 << 16) + 65536)`:
- bytes `0..31`: block header (`len` at `+8`, `seq` at `+16`),
- bytes `32..32799`: 1024 32-byte records (`n = 0..1023`).

Record `n = 1023` writes directly to offset `32768..32799` (where bit 15 is 1),
requiring no permission fault and no tail-alias reissue. `EM` (and `xemit`) executes
`tbnz \s, #15, gf\l`, jumping to `FLUSH7` which resets `LC1[d1] = base + 32`
(`sub #32768; add #32`), writes the header `(32768, seq+1)`, clears G's local scalar ACK
register at `0x1800_0000_0000`, executes `DSB SY`, writes `RREADY` (`0x0800_0000_0000 + 8*d1`),
executes `DSB SY`, and polls ACK with `WFE` (`FLUSH = 18` straight, 26 at the wait bound).
The endpoint maps `H1` displacement `65536*d1 + offset` directly to relay `d1`'s
preallocated 64-KiB receive slot `196 TiB + 65536*m + offset` (member `m in 0..1023`).
On `R`, `bulk_r_receive` (15 straight loop instructions, `R_RECV = 23` at the wait bound)
wakes on `RREADY`, binds `SP` to member `m`'s preinitialized 512-pointer source table
(8-KiB stride at `192 TiB + 32 MiB`; its words at `+4096` and `+4104` hold the slot's first
record and `G`'s ACK alias), calls `h2s`, issues `DSB SY`, and writes `G`'s scalar ACK. A
group's partial blocks go by `fin1` with the same header and handshake (`seq` counts the
`G`'s blocks for that relay) through `h2s`'s partial path.

Each block-end handshake wait (`FLUSH7`, `fin1`, `hfshared`, `fin`, `xemit`, `bulk_r_receive`,
`bulk_h_receive`) polls with `LDR; CBNZ` (or `LDR; CBZ; BR`) and otherwise `WFE; B`: 2 or 3
instructions when the ACK is present and 4 more for each failed pass. One stale event plus one
publication, with no spurious wakes, allow at most two failed passes, so each side is charged its
executed straight path plus 8.

As in `283cdacf`, `G`'s main processor (`PE0`) and member-offer consumer (`PE9`)
carry fixed 4-bit initiator port constants on `USER.srcPE`. A nonzero `RREADY`
commit latches `srcPE` and sets a valid bit; a second commit while valid returns
`SLVERR`. `R`'s ACK write returns to the latched `USER.dstPE` and clears valid.
Only one `G` arithmetic source is active globally per member-offer/epoch turn, and
every source clears its ACK level before publishing `RREADY`.

### 2.2 R-to-H (`H2U` contiguous 32-KiB stride, `K2B = 512`)

On relay `R`, the contiguous hop-2 block window `H2U` occupies
`[8 EiB, 8.125 EiB) = [0x8000_0000_0000_0000, 0x8200_0000_0000_0000)` (`2^42`
homes x `32 KiB` stride; above the relay's candidate-record buffer, which `RELAY_CAP` sizes).
For home `d2 < 2^42`, block `d2` occupies
`[H2U_BASE + (d2 << 15), H2U_BASE + (d2 << 15) + 32768)`:
- bytes `0..31`: header/scratch,
- bytes `32..16415`: 512 32-byte records (`n = 0..511`).

Record `n = 511` writes directly to offset `16384..16415` (bit 14 = 1). `H2B` executes
`tbnz x11, #14, hf\l`; `hf\l` is a two-instruction trampoline (`adr x19, hn\l; b hfshared`,
zero exception entry) into one shared stub that resets `LC2[d2] = base + 32`, advances
`RC2[d2]` by `16384` (once `RC2[d2]` has reached the hop-2 cap `CAP_H`, the stub leaves it at the cap
and returns without a flush), clears R's scalar ACK at `0x8008`, writes `HREADY`
(`0x2000_0000_0000 + LC2_addr`), polls ACK with `WFE`, and returns by `br x19` (`FLUSH2 = 21`
straight, 29 at the wait bound). The endpoint transfers payload bytes `32..16415` directly
into home `d2`'s 64-byte-aligned `R2` cursor (`[2^55, 2^55 + 32767 TiB)`), bounded
by `CAP_H`. On `H`, `bulk_h_receive` (`BULK_H = 33` per block inside `BLOCKH = 582`:
25 of its 27 static instructions run straight, and 33 at the wait bound) wakes on
`HREADY`, runs `dlvh` directly over `R2` (`dlvh` non-`HCHGI` measured `520 = DLVHB`
insns, `NEW_INSERT = 3`, `MAIN_VISIT = 10`), advances the `R2` cursor, and writes
`R`'s scalar ACK. Partial end-of-group tails use `fin` (`add x8, x8, #0x8000`,
step `32 KiB`) with the same explicit scalar handshake.

### 2.3 Owner-Shared Y Transport and Generic Scalar Service (`MSG = 256`)

The current Y staging source uses `N_CLUST = 2^48`, with 33 slice bits and
15 low cluster bits. Every group-relative chain index `r` uses
`bank = r >> 27`, `slot = r & (2^27-1)`. The maximum group has133854748188
chains, hence998 active banks out of1024. Each bank has the following
backed regions (all intervals are half-open):

| Region | Begin | Length |
|---|---:|---:|
| Payload | 0 | 2^62 bytes |
| Metadata | 2^62 | 2^31 bytes |
| Payload and command inboxes | 2^62+2^32 | 2^24 bytes |
| Bank control | 2^62+2^33 | 2^24 bytes |
| Generic service, code and stacks | 2^62+2^34 | 2^24 bytes |

Each chain slot reserves2^35 payload bytes and16 metadata bytes. Metadata
contains the saturated frame count and the absolute chain tag. Payload
frame n begins at `(slot<<35)+(n<<9)+64`. The2^26-1 cap keeps its full448
payload bytes inside that slot. Unused virtual gaps are not allocated;
certificate.py sums the backed intervals and all2^58 banks explicitly.

G source control is34TiB+0x7200000: route(owner48,bank10) at+0/+8, ACK+16,
READY+24 and512-byte outgoing aperture+4096. The archive port targets
B inbox `2^62+2^32` and READY at BCTRL+24, where BCTRL=`2^62+2^33`.
B PE1 runs the actual yreceive copy. Its return route at BCTRL+32/+40
and ACK+48 addresses the producer G ACK on PE0. ErrorACK2 is a nonzero
event. The sender retains its data until the receiver's ACK.

B replay uses the same outgoing aperture with replay class and route(owner,0).
It targets G endpoint(owner<<10), inbox34TiB+0x7510000 and owner control
34TiB+0x7500000 at READY+24, PE8. Owner return words and ACK+48 target
BCTRL+16 on B PE0. Both endpoints execute the448-byte copy through scalar
instructions. yreplay_request checks the absolute tag before reading data.
The owner's slot is zero; the shared receiver's larger bound does not change
its buffer address. ybank_begin writes each active metadata entry before use.

Specialized commands use G control34TiB+0x7600000, aperture+4096, ACK+16
and READY+24. Command inboxes are G34TiB+0x7610000 or
B`2^62+2^32+65536`. Receiver control is the source control base+512.
The84-bit route selects a valid G or B endpoint. The receiver saves the
source route and issues completion only after begin/seal/release/replay
work completes. The controller at G index2^58 enumerates all2^48 owners
and the needed bank prefix. Capacity guards reject a group above2^37.

Local member offers use the explicit yoffer/yconsume handshakes across the
1024 members of one owner. Only the selected member may consume the shared
frame. ACK follows actual yconsume return before the owner reuses that
frame. The member read-only alias [38TiB,38TiB+2^35) maps to the
owner backing [40TiB,40TiB+2^35), with no extra backing allocation; it
performs no implicit software callback. The actual ychunk driver
keeps the128-bit fixup meter and invokes the counted fwd/inverse/bwd/keys
routines.

Generic bases are G34TiB+0x7800000, R192TiB+20MiB, H2^55-32MiB, and
B2^62+2^34. At each base the first three64KiB pages are control, remote
and ACK aliases; the next two are inbox and handler state. B code is
separately reserved at generic base+1MiB, within the16MiB backed region.
A fixed initiator index selects one of16 route slots. A counted store to
REMOTE+256 publishes the already copied request to receiver PE15. The
receiver saves the source route and selector before calling its bounded
copy/queue handler, then replies through its private state. Only that
receiver's actual ACK releases the sender. Bulk, Y and generic services have
distinct PE/event assignments; generic receipt alone does not certify
application completion. The phase algorithm in callbacks.md supplies the
separate application work and drain barriers.

Generic messages are charged256 instructions each, including their counted
send/receive/ACK paths and descriptor glue. Scalar callbacks are charged
GEN_MUL=912688, six RESOLVE_RESERVE=2^21 a resolve and six REQ_RESERVE=4096 a request inside the
certificate. Each endpoint has16 reserved interfaces/PE slots; their state
and program delivery are in the cold census. The larger bank command inbox
address adds one scalar instruction per ycmd_bank call. The per-group bank
allowance pays4096 instructions per specialized command and8 per metadata
slot, above the measured source envelopes; certificate.py applies it to
each active owner/bank and the full group chain population.

The instruction model includes the stated device ordering and event behavior
in proof.md. This construction does not assume that a bus write waits for
a software routine to return. Every wait for application completion is a
separate scalar protocol step. Kernel tests exercise finite service paths;
they do not synthesize this fabric or run its full production population.

### 2.4 Table Banks (Role `T`) and the T Aperture

**Banks.** A table bank is a passive memory endpoint with the same local address space `[0, 2^64)` as every other processor. It backs
`[2^62, 2^63)` (4 EiB: its table rows), a 16-MiB code and control region at `0x10000` and the generic service pages at
`192 TiB + 20 MiB`, as `R` does (`certificate.py`: `T_MEM_REGIONS`, checked in range and disjoint). Its fabric tag is role code `B` with
index `2^58 + (pool << 42 | FAM << 36 | bank)`, disjoint from the archive banks' indices (`(owner48 << 10) | bank10 < 2^58`);
every other index names no endpoint. There are `N_T = 2,216,615,441,596,416` banks (`proof.md` Section 3.4 gives the inventory). Before the leaf kernels run,
a bank's own processor builds its rows into `[2^62, 2^63)` with counted instructions (`callbacks.md`, "Table banks"), ends with `DSB SY` and sends a
ready message; afterwards the bank runs no software for reads, it answers read requests from its memory.

**Aperture of a reader.** The local range `[2^62, 2^63)` of every arithmetic generator `G` is not RAM: it is an aperture of the generator's
local decode endpoint (`certificate.py`: `SERVICE_ALIASES`, checked disjoint from every backed region of `G`). The endpoint has two registers:
`POOLREG` (14 bits), a device register at local address `20 TiB` (`SERVICE_ALIASES`, outside every backed region and every other alias of `G`), written once at boot by
a counted store of the generator's pool number (`POOLREG_INIT = 16` instructions, which cover that store and the pool-number computation that the mining
callback already performs), and the *route latch* `L` (64 bits), written by any store to the aperture. It serves the processor's accesses in
program order, one transaction at a time (the bridge's `MaxTrans = 1`, Section 1):
1. a store of any width to an aperture address sets `L` to the stored value and completes; no data goes anywhere;
2. a load from an aperture address `a` captures `L` when it enters the endpoint and becomes a read request (the `AR` channel of the tree,
   whose switches carry the five AXI channels, Section 1) with destination tag index `2^58 + (POOLREG << 42 | L[41:0])` and address `a`;
   the 64 bytes at `a` of that bank come back on the `R` channel into the eight registers of the `LD64B`, routed to the waiting generator by the
   source tag that every request carries in the `USER` field of the tree (two 84-bit tags, destination and source, Section 1), as for any read. The processor
   waits for the response; waiting is not an instruction.
The latch and `POOLREG` belong to the bridge of the processing element that runs the leaf kernels (`G`'s main processor `PE0`, Section 2.1; its bridge is the one whose `MaxTrans = 1`
orders its accesses). No other element of the generator reads or writes the aperture. An element that did would have its own bridge, latch and `POOLREG`, written at boot by its
own counted store: `POOLREG_INIT = 16` covers the two stores of `G`'s two initiators (`PE0` and the member-offer consumer `PE9`), so the price does not depend on which of them runs the reads.
A store after a load cannot change the destination of that load, and one latch value serves any number of loads issued before the next latch
store (the two lines of a 256-byte row of `E`, `G`, `F` use one; the two lines of a `fwd` row are read in two passes, each after its own store). A load
from the aperture whose latch names no bank is an invalid table-bank transaction; the kernels' loads never meet one (`experiment.py` asserts that the digit decoded from every aperture read fits its window). A latch store that no load follows is harmless. `POOLREG` is a per-generator register because the `Y` frames that carry the tags of `E`, `F`, `G` and `fwd` are
shared by the 1,024 members of a cluster, whose pools differ.

**Price.** A read costs the instructions the processor executes: the extractions that form the tag and the row address (`extr`, `bfi`,
`bfxil`), one `str` of the latch, and one `ld64b` a 64-byte line. No acknowledgement, `DSB`, `WFE` or polling is executed for a read: it changes no
state at the bank and its response is its completion. The endpoint's two registers and the switch state of the `N_T` bank endpoints are counted in
`FABRIC_BYTES` (4 switches of 64 KiB and 1 MiB of bridge and control state a physical endpoint). 

**Every step of a read.** The rules count executed instructions: the reference prices are calibrated with Callgrind's instruction event `Ir` "rather than measuring elapsed time", and `R` is "work in scalar AArch64 instruction equivalents,
including arithmetic, loads/stores, addressing and control" (`COST_MODEL.md`). In that count a load is one instruction whichever memory it reaches, and the work of caches, address translation, interconnect and memory arrays does not appear in it;
memory is reported as peak memory and has no scalar contribution. This package gives address decode, switches and a memory array's read no instruction. Each step of a read, who performs it and how it is counted:

| step | performed by | counted as |
|---|---|---|
| form the digit, the bank field, the row address and the tag word | the generator (`extr`, `bfi`, `bfxil`, `mov`, `movk`) | scalar instructions in `FWD_INC`, `BWD_INC`, `YR_BLK` |
| set the latch `L` | the generator: one `str` to the aperture | one instruction |
| issue the read | the generator: one `ld64b` from the aperture | one instruction |
| capture `L`; build the request with the destination tag `2^58 + (POOLREG << 42 \| L[41:0])` and the generator's source tag | the generator's decode endpoint: two registers and combinational logic | no instruction; the registers are counted in `FABRIC_BYTES` |
| carry the request and the 64-byte response | the tree of switches (five AXI channels, tags in `USER`) | no instruction; the switches are counted in `FABRIC_BYTES`; no latency, bandwidth or contention is priced |
| read the 64 bytes | the bank's memory | no instruction; the bank bytes are counted in `PEAK_BYTES` |
| deliver the 64 bytes to the eight registers | the response, which completes the `ld64b` | the `ld64b` |
| write the rows and publish them | the bank's own processor | counted instructions (Rows 8 and 9), `DSB SY`, ready message |

**Why a latch.** The rows of one pool's tables occupy `63 * 2^93` bytes (`2^98.98`) and a window of 85 bits has `2^85` rows, so no 64-bit load address reaches every row of a pool. A read needs a destination wider than the address of the load; the
generator supplies the missing bits by writing a register, which is a counted store, and the endpoint decodes the pair `(L, a)` as it decodes address bits to select a relay or a home elsewhere in this document, with more bits. Hardware has this shape: a window of local addresses whose
destination a software-written translation register selects (as in non-transparent bridge windows), so that a core reaches remote memory with ordinary loads and stores; this package idealizes it (no latency, bandwidth or contention is priced). A reviewer who prices the decode
at `k` instructions a read changes the bound by about `0.022` bit a unit of `k` (`proof.md` Section 3.4).

**Orientation and counts.** The base package's read-only `YVIEW` alias of the
owner's memory (Section 2.3) is an access of the same kind, a counted `LD64B` of remote memory through the decode with no software at the owner (named to locate the convention in this package's own text, not as authority for it);
only the destination of a table read comes from the latch. 8cc64cb7 prices remote reads of the owner's frames the same way in its hot loop: `fwd` reads them with
`ld64b x16, [x1]` (`LKP6_PAIR`, `LK1` in `kernels/leaf.s`), `bwd` runs with `sp` on the frame, and `kernels/ytransport.s` passes `yconsume` the "frame base (YVIEW or
owner's Normal YQ)"; it charges only the per-frame offer and consume handshake (`OFFER`, `Y_TRANSPORT`), because the owner writes that frame. A table bank is read-only after
the barrier, so no handshake exists to charge. Over the full horizon the aperture serves at most `2^130.1165` `LD64B` (15 an incidence), an average of `2^79.1392`
over the `N_T` endpoints (not a bound on one endpoint; no endpoint throughput is priced). No instruction of a read orders it against its bank's construction: a bank's `DSB SY`
before its ready message and the barrier that follows the last ready message do. The model prices the instructions a processor executes. A `LD64B` returns 64 bytes of table
payload (960 bytes an incidence, `proof.md` Section 3.4); its request and response are the transactions of any counted load, nothing executes at the bank or on the path,
and, as for the owner alias, no bytes in flight or protocol cycles are priced.

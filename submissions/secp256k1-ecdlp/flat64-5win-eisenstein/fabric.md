# Memory Fabric and Explicit Scalar Services (Flat 64-Bit Analytical Machine)

This document specifies the memory-fabric construction and explicit scalar service
bindings used by `flat64-5win-eisenstein`. It adapts the explicit scalar push/ACK
fabric of live leader `283cdacf` (`soispoke`) to the uniform 64-bit flat analytical
address space (`[0, 2^64)` per processor, matching accepted submission `4e2d5241`)
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
and `N_HOME = 2^82` homes (`0 .. 2^82 - 1`). Unused indices are decode errors.
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
`E = 2^82 + 2*2^58 + 2^40 + 1` physical endpoints, a full padded 84-bit tree has
`2*(2^84 - 1)/15 < E` internal switches; reserving 4 switches of `64 KiB` plus
`1 MiB` of bridge/control state per endpoint is strictly conservative and is
included in `FABRIC_BYTES` inside `PEAK_BYTES < 2^138 B`.

## 2. Endpoint Decode and Contiguous Block Windows (No Faults, No Tail Aliases)

### 2.1 G-to-R (`H1` contiguous 8-KiB stride, `K1 = 128`)

On generator `G`, the contiguous hop-1 block window `H1` occupies
`[64 PiB, 72 PiB) = [0x0100_0000_0000_0000, 0x0120_0000_0000_0000)` (`2^40`
relays x `8 KiB` stride). For relay `d1 < 2^40`, block `d1` occupies
`[H1_BASE + (d1 << 13), H1_BASE + (d1 << 13) + 8192)`:
- bytes `0..31`: block header (`len` at `+8`, `seq` at `+16`),
- bytes `32..4127`: 128 32-byte records (`n = 0..127`).

Record `n = 127` writes directly to offset `4096..4127` (where bit 12 is 1),
requiring no permission fault and no tail-alias reissue. `EM` (and `xemit`) executes
`tbnz \s, #12, gf\l`, jumping to `FLUSH7` which resets `LC1[d1] = base + 32`,
writes the header `(4096, seq+1)`, clears G's local scalar ACK register at
`0x1800_0000_0000`, executes `DSB SY`, writes `RREADY` (`0x0800_0000_0000 + 8*d1`),
executes `DSB SY`, and polls ACK with `WFE`. The endpoint maps `H1` displacement
`8192*d1 + offset` directly to relay `d1`'s preallocated 8-KiB receive slot
`0xc000_0000_0000 + 8192*m + offset` (member `m in 0..1023`). On `R`,
`bulk_r_receive` (charged conservatively at `45` insns per full 128-record block;
15 static loop instructions) wakes on `RREADY`, binds `SP` to member `m`'s
preinitialized 64-pointer source table, calls `h2s`, issues `DSB SY`, and writes
`G`'s scalar ACK.

As in `283cdacf`, `G`'s main processor (`PE0`) and member-offer consumer (`PE9`)
carry fixed 4-bit initiator port constants on `USER.srcPE`. A nonzero `RREADY`
commit latches `srcPE` and sets a valid bit; a second commit while valid returns
`SLVERR`. `R`'s ACK write returns to the latched `USER.dstPE` and clears valid.
Only one `G` arithmetic source is active globally per member-offer/epoch turn, and
every source clears its ACK level before publishing `RREADY`.

### 2.2 R-to-H (`H2U` contiguous 16-KiB stride, `K2B = 256`)

On relay `R`, the contiguous hop-2 block window `H2U` occupies
`[1 EiB, 1.0625 EiB) = [0x1000_0000_0000_0000, 0x1100_0000_0000_0000)` (`2^42`
homes x `16 KiB` stride). For home `d2 < 2^42`, block `d2` occupies
`[H2U_BASE + (d2 << 14), H2U_BASE + (d2 << 14) + 16384)`:
- bytes `0..31`: header/scratch,
- bytes `32..8223`: 256 32-byte records (`n = 0..255`).

Record `n = 255` writes directly to offset `8192..8223` (bit 13 = 1). `H2B` executes
`tbnz x11, #13, hf\l`; the cold stub `hf\l` (18 instructions, zero exception entry)
resets `LC2[d2] = base + 32`, advances `RC2[d2]` by `8192` under `RELAY_CAP`,
clears R's scalar ACK at `0x8008`, writes `HREADY` (`0x2000_0000_0000 + LC2_addr`),
and polls ACK with `WFE`. The endpoint transfers payload bytes `32..8223` directly
into home `d2`'s 64-byte-aligned `R2` cursor (`[2^55, 2^55 + 38220 TiB)`), bounded
by `CAP_H`. On `H`, `bulk_h_receive` (charged conservatively at `BULK_H = 48`
insns per block inside `BLOCKH = 346`; 24 static loop instructions) wakes on
`HREADY`, runs `dlvh` directly over `R2` (`dlvh` non-`HCHGI` measured `264 <= 272`
insns, `NEW_INSERT = 3`, `MAIN_VISIT = 10`), advances the `R2` cursor, and writes
`R`'s scalar ACK. Partial end-of-group tails use `fin` (`add x8, x8, #0x4000`,
step `16 KiB`) with the same explicit scalar handshake.

### 2.3 Owner-Shared Y Transport and Generic Scalar Service (`MSG = 256`)

Sub-ring Y staging (`ystage.s`, `ytransport.s`) and generic control messages
(`generic.s`) are preserved verbatim from `283cdacf` with `N_CLUST = 2^48`
(`K_S = 33` slice bits, `c' mod 2^15` in the low 15 bits) and `w3` bit layout
`((k-1) >> 36) | ((c' >> 15) << 9) | (u << 42)` (`lsl #9` in `ystage.s`). Every
generic message is charged at `MSG = 256` scalar instructions, and every scalar
callback (`GEN_MUL = 2^20`, `RESOLVE_RESERVE = 2^21`, `REQ_RESERVE = 4096`) is
charged directly inside `certificate.py` across the complete full-size populations.

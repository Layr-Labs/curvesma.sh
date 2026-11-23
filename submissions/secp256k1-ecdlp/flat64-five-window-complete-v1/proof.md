# Flat 64-Bit Five-Window Eisenstein Stored-Leaf Birthday Search on secp256k1 (123.149554427177781)

## Claim

For uniform unrestricted production secp256k1 ECDLP challenges `Q = [x]P` and fresh
algorithm coins, this cold classical algorithm recovers `x` with probability at
least `0.39` (`Q = O => x = 0`) at expected classical work at most
`2^123.149554427177781` reference units (`C = 1,000` scalar AArch64 instructions per
unit), exact rational `time_log2 = 123149554427177781/1000000000000000` (`123.149554427177781`).
This improves on the live leader `283cdacf` (`123.186563709484764`) by at least
`0.037009282306983` bits.

All work across every processor is charged: cold initialization, Phase-0 `X`-pair
mining, table construction, coins, arithmetic, routing, home joins, candidate
verification, explicit scalar memory-fabric services (`MSG = 256`, `GEN_MUL = 2^20`),
exception ceilings, barriers, phase control, and output. Preprocessing is at most
`2^115.7843` units (included in `E_run`), advice is zero, and peak reserved memory
across all nodes and fabric switches is `348,385,351,113,989,259,397,648,273,404,705,193,740,832`
bytes (`<= 2^137.999736 B < 2^138 B`).

Summary figures (all produced directly by `certificate.py` and verified in Docker by
`experiment.py`):
- Gate (`route.s:gate`, `GATE_INSNS = 256`): `K_GATE = 132756105517980808685984632952518166242`,
  `p_gate = K_GATE / 2^128 <= 0.390135130183885`.
- Opened-run success lower bound: `lo_run >= 0.999653632361`; exact rational check
  `p_gate * lo_run >= 39/100` (and `(K_GATE - 1)/2^128 * lo_run < 39/100`).
- Opened-run expected work: `E_run <= 2^134.473293` instructions (`<= 2^124.507509`
  units, `<= 51.2251` insns per full-horizon leaf).
- Total expected work: `E_total = 256 + p_gate * E_run <= 2^133.115339` instructions
  (`<= 2^123.149554427177781` units).

## 1. Lineage, Credits, and Substantive Change

This submission builds directly on two published live-board designs:
1. **`283cdacf` (`soispoke`, `123.186563709484764`)**, which introduced the explicit
   PULP AXI memory-fabric construction (`fabric.md`), counted scalar generic
   (`MSG = 256`) and callback (`GEN_MUL = 2^20`) services (`callbacks.md`),
   owner-shared `Y` transport (`Y_TRANSPORT = 129`, `OFFER = 576`), direct-`R2`-pointer
   home hashing (`NEW_INSERT = 3`, `MAIN_VISIT = 10`), direct-push `R -> H` delivery
   into `R2` (`K2B = 256`), dedicated `SP` source-pointer tables on `R` and `H`,
   `JK` word-7 seed linking, `EB` static pointer pairs, the pointwise joint-label
   derivation, global `PAIR_INVERSE_EXC_WORK` charging, and the denominator-cleared
   selected-`X` fiber bound (`selected-x.md`), building on `4ce2d07e` (`bradcush`),
   `d719e7cb` (`Jaydearcadian`), and `7bd265de` (`kaankolcu`). We also incorporate
   the `Q1GD` retune ratio (`106/103.079`) credited to `c4291c3b` (`aburan28`).
2. **`4e2d5241` (`DeanMojo`, `db1d7ca7`, `123.324660124`)**, which supplies an earlier concrete
   uniform 64-bit flat analytical address space `[0, 2^64)` per processor with
   explicit `TBZ`/`TBNZ` block ends and 52-bit `BT` windows (building on `4d4e224b`
   by `GoJackzi`, `7f03684c` by `hosein-ul`, and `94ae0ccc` by `kaankolcu`).

**Single architectural change relative to `283cdacf`:**
We replace `283cdacf`'s Armv9.6-A VMSAv9-128 two-55-bit-range / 56-bit-PA +
stage-1 permission-fault machine with `4e2d5241`'s uniform **64-bit flat analytical
address space** `[0, 2^64)` per processor, rescale generator count to `NODES = 2^58`
(`N_CLUST = 2^48`, `K_S = 33` slice bits, `Q1GD = 27,136,000,000,000 < 2^45`),
upgrade all four 256-bit six-window lookup tables on `G` (`BT`, `P1`, `E`, `G`) to
**five windows `(52, 51, 51, 51, 51)`** (`BLK5F`, `P1_5WIN`, 5-window `JK`,
`OFFS5_52`) while keeping `F` and `fwd` at five windows `(45, 45, 45, 44, 44)`,
replace stage-1 permission faults with explicit `tbnz` block ends (`EMIT = 7`,
`H2S_ARR = 6`, `EXC_ENTRY = 0`), rescale `kernels/ystage.s` sweeps and coin/boot
checks to `NODES = 2^58`, `N_CLUST = 2^48`, controller `(2^58, 0)`, and `k <= 2^45`,
and integrate a 128-bit cumulative fixup meter directly into `kernels/ytransport.s:ychunk`
(`DELTA_F = 2^71`). Relays (`N_REL = 2^40`, `K1 = 128`, `K2B = 256`) and homes
(`N_HOME = 2^82`, `HB = 50`, `R2`, `NX`, `HH`, `dlvh`) retain their exact block
sizes, capacities, and `283cdacf` direct-`R2` push architecture.

## 2. Analytical Machine Model, Premises, and Complete 64-Bit Flat Maps

### 2.1 Explicit Machine Premises

This is an **analytical machine entry** under `tracks/docs/COST_MODEL.md` ("work in
scalar AArch64 instruction equivalents, including arithmetic, loads/stores,
addressing and control... Count other work on the same scalar 64-bit instruction
basis... memory has no scalar contribution or tie-break"). We do **not** claim
Armv9.6-A 56-bit physical hardware addressing on `G` or `R`.

Disclosed premises:
1. **Uniform 64-bit flat local address space `[0, 2^64)` per processor** (same
   analytical convention disclosed in `4e2d5241`): only explicitly backed regions below
   count toward `PEAK_BYTES`; every scalar load/store targeting a backed address in
   `[0, 2^64)` counts as 1 instruction. No VMSAv9-128 MMU, page tables, `FEAT_D128`,
   `FEAT_LVA3`, `FEAT_S1PIE`, or stage-1 permission faults are used.
2. **AArch64 ISA extensions** (verified by `kernels/a64.py` encoding rules and the finite source checks):
   base AArch64 scalar integer instructions, `FEAT_LS64` (`LD64B`/`ST64B`, restricted
   to 64-byte-aligned backed addresses and even base registers `<= 22`), `DC ZVA`
   (2,048-byte zeroing granule on `H`), `FEAT_RNG` (`RNDR` for fresh unbiased coins),
   `FEAT_CMPBR` (`CBHS`/`CBHI`/`CBEQ`), and `WFE`/event wake semantics without
   spurious unrelated wakeups (two-publication control waits charged explicitly as in
   `283cdacf`).
3. **Explicit scalar memory fabric (`fabric.md`)**: every cross-node transfer is a
   counted scalar store/load/ACK over the 64-bit-address, 84-bit-tag PULP AXI tree.
4. **Heuristic / randomness premises**: uniform secp256k1 challenge `Q = [x]P` and
   independent fresh `RNDR` coins; finite experiments in `experiment.py` check
   assembly, 5-window arithmetic, exact EC labels on finite samples, routing flushes,
   meter behavior, and map disjointness, never production runtime or RTL synthesis.

### 2.2 Backed data and service maps

The following exact intervals are generated from certificate.py. Each end is
exclusive. Scalar accesses to device apertures and the read-only owner alias
use the separate decode definitions in fabric.md; these are not additional
RAM arrays. Programs are placed within the listed backing at G34TiB+1MiB,
R/H0x10000, and Bgeneric-base+1MiB, respectively. Each role is additionally
charged a conservative PROGRAM_BYTES+2MiB reserve beyond its interval sum;
this extra charge is retained even where program backing overlaps a listed
region. Every tested service, stack and cursor span is checked for containment.

**B**: 4611686020628348928 charged bytes per endpoint.

| Region | Begin | End | Backed bytes |
|---|---:|---:|---:|
| `B_DATA` | `0x0` | `0x4000000000000000` | 4611686018427387904 |
| `B_META` | `0x4000000000000000` | `0x4000000080000000` | 2147483648 |
| `B_INBOX` | `0x4000000100000000` | `0x4000000101000000` | 16777216 |
| `B_CTRL` | `0x4000000200000000` | `0x4000000201000000` | 16777216 |
| `B_CODE_GENERIC` | `0x4000000400000000` | `0x4000000401000000` | 16777216 |

**G**: 13873807722055092768 charged bytes per endpoint.

| Region | Begin | End | Backed bytes |
|---|---:|---:|---:|
| `LC1` | `0x0` | `0x80000000000` | 8796093022208 |
| `YR` | `0x200000000000` | `0x210000000000` | 1099511627776 |
| `G34_CTRL_MMIO` | `0x220000000000` | `0x260000000000` | 4398046511104 |
| `YQ_YVIEW` | `0x27ffffffffc0` | `0x281080000000` | 70866960448 |
| `EB` | `0x2a0000000000` | `0x2a2000000000` | 137438953472 |
| `CN` | `0x540000000000` | `0x67f85e9000000` | 1736704000000000 |
| `CNM` | `0x8000000000000` | `0x85a995c000000` | 99614720000000 |
| `H1_CONTIG` | `0x100000000000000` | `0x120000000000000` | 9007199254740992 |
| `F_TBL_5WIN` | `0x800000000000000` | `0x880000000000000` | 36028797018963968 |
| `FWD_TBL_5WIN` | `0x880000000000000` | `0x8c0000000000000` | 18014398509481984 |
| `P1_TBL_5WIN` | `0x1000000000000000` | `0x5800000000000000` | 5188146770730811392 |
| `E_TBL_5WIN` | `0x6000000000000000` | `0x9000000000000000` | 3458764513820540928 |
| `G_TBL_5WIN` | `0xa000000000000000` | `0xd000000000000000` | 3458764513820540928 |
| `G_BATCH_META` | `0xd000000000000000` | `0xd003174e15fa41a0` | 870049073742240 |
| `REQ_CURSOR_MONITOR` | `0xd7fffffffffff800` | `0xd800000000000000` | 2048 |
| `REQ_64B` | `0xd800000000000000` | `0xe7a0000000000040` | 1125899906842624064 |
| `BT0_52` | `0xf800000000000000` | `0xfc00000000000000` | 288230376151711744 |
| `BT1_51x4` | `0xfc00000000000000` | `0x10000000000000000` | 288230376151711744 |

**H**: 72043300212119552 charged bytes per endpoint.

| Region | Begin | End | Backed bytes |
|---|---:|---:|---:|
| `H_CTRL_CODE` | `0x10000` | `0x1010000` | 16777216 |
| `H_HCC_CAND` | `0x7ffc0000000000` | `0x7ffcfff3361800` | 1099297069056 |
| `H_WIN_BUF` | `0x7fff0000000000` | `0x7fff0001000000` | 16777216 |
| `H_SP_HDR` | `0x7ffffffff00000` | `0x80000000000000` | 1048576 |
| `H_BULK_CTRL` | `0x7ffffffc000000` | `0x7ffffffc010000` | 65536 |
| `H_GENERIC` | `0x7ffffffe000000` | `0x7ffffffe200000` | 2097152 |
| `R2_38220T` | `0x80000000000000` | `0x1154c0000000000` | 42023334413598720 |
| `HH_2p50` | `0xff80000000000000` | `0xffa0000000000000` | 9007199254740992 |
| `NX_19110T` | `0xffa0000000000000` | `0xffeaa60000000000` | 21011667206799360 |

**R**: 107804916097073152 charged bytes per endpoint.

| Region | Begin | End | Backed bytes |
|---|---:|---:|---:|
| `RCS_CTRL` | `0x8000` | `0x114000` | 1097728 |
| `R_SLOTS_1024x8K` | `0xc00000000000` | `0xc00000800000` | 8388608 |
| `R_SLOT_TABLES` | `0xc00001000000` | `0xc00001100000` | 1048576 |
| `R_GENERIC` | `0xc00001400000` | `0xc00001600000` | 2097152 |
| `LC2` | `0x200000000000` | `0x400000000000` | 35184372088832 |
| `RC2` | `0x600000000000` | `0x800000000000` | 35184372088832 |
| `E2P` | `0x800000000000` | `0xc00000000000` | 70368744177664 |
| `R_SP_TBL` | `0xfffffff00000` | `0x1000000000000` | 1048576 |
| `DPW_RCB` | `0x1000000000000` | `0x7f800000010000` | 35606584553963520 |
| `H2U_CONTIG` | `0x1000000000000000` | `0x1100000000000000` | 72057594037927936 |

The population is2^58 arithmeticG plus one controller,2^40 relays,2^82 homes and2^58 banks.
The total is `(2^58+1)*G_BYTES+2^40*R_BYTES+2^82*H_BYTES+2^58*B_BYTES+FABRIC_BYTES = 348385351113989259397648273404705193740832` bytes, at most`2^137.999736 <2^138`.

## 3. Five-Window Kernels, Record Packing, and Measured Instruction Counts

### 3.1 Five-Window `(52, 51, 51, 51, 51)` Decomposition (`BT`, `P1`, `E`, `G`)

Every 256-bit integer `w = (a, q, c, d)` decomposes uniquely into five unsigned
digits at bit positions `(0, 52, 103, 154, 205)` of widths `(52, 51, 51, 51, 51)`:
- `w0 = a[51:0]` (`bfi dst, a, #6, #52` for 64-byte rows; `#8, #52` for 256-byte rows),
- `w1 = extr(q, a, #52)[50:0]`,
- `w2 = extr(c, q, #39)[50:0]`,
- `w3 = extr(d, c, #26)[50:0]`,
- `w4 = (d >> 13)[50:0]`.

All five digits fit legally inside 64-bit `BFI`/`ORR` addressing (`experiment.py`
tests 5,000 random 256-bit vectors against exact integer reconstruction and region
alignment):
1. **`BLK5F` (`kernels/leaf.s`, 30 insns vs `BLK6F`'s 36 insns, runs twice per incidence):**
   `x7 = 0xf800_0000_0000_0000` (`BT0`, 64-byte rows) is overwritten once by
   `bfi x7, \a, #6, #52; ld64b x20, [x7]`. `x8 = 0xfc00_0000_0000_0000` (`BT1`
   slots 0/1) and `x9 = 0xfc00_0000_0000_0040` (`BT1` slots 2/3) have top 6 bits
   equal to `1`, so `orr x10, x8/x9, idx, lsl #7` masks bits `[63:58]` automatically
   for each 51-bit digit (`lsl #7` = 128-byte stride).
2. **`P1_5WIN` (`kernels/leaf.s`):**
   Loads five window bases `x12..x16` per `X_j` (base `1 EiB + 0.75*j EiB`, aligned
   to `2^58` for `win 0` and `2^57` for `wins 1..4`) in one `ld64b x12, [x29]`,
   computes five 64-byte-aligned row addresses, and drops the sixth `ld64b + P1W6_ACC`
   window completely.
3. **Five-window `E` and `G` in `JK` (`kernels/leaf.s`) and `OFFS5_52` (`kernels/secp_fe.s`):**
   `OFFS5_52` writes five 256-byte row pointers for `E` (at frame offsets `280..312`,
   `win 0` at `+64` for word-7 seed linking) and five for `G` (`472..504`). In `JK`,
   `E` uses `KLKA5E`/`KLKB5E` (`KEB` before `KWINA/KWINB x29, Z`) and `G` uses
   `KLKA5`/`KK5`, removing one `KWINA` (18 insns) and one `KWINB` (19 insns) from
   each of `E` and `G` (`-74` insns per `Y` = `-12.333` insns/incidence).
4. **Five-window `(45, 45, 45, 44, 44)` `F` and `fwd` (`K_S = 33` slice bits):**
   With `N_CLUST = 2^48 = 2^33 (s) * 2^15 (c' mod 2^15)`, `yround` extracts
   `s = x(Y)^3 >> 223` (`lsr x14, x13, #31`), leaving 223 varying bits decomposed by
   `OFF5PAIR` into `(45, 45, 45, 44, 44)` at shifts `(0, 45, 90, 135, 179)`. `fwd1`
   sits at `fwd0 | 0x0020_0000_0000_0000` (`orr x0, \b, #0x20000000000000` in `LK1`).

### 3.2 Record Packing at `NODES = 2^58`, Explicit `TBZ`, and 128-Bit Fixup Meter

- **Record `w2, w3` (127 bits):** With `NODES = 2^58`, `M_G = 1024`,
  `N_CLUST = 2^48`, `Q1GD = 27,136,000,000,000 < 2^45`, and `u < 2^21`:
  `w2 = j | (m << 3) | ((c' & 0x7fff) << 13) | (((k-1) & 0xfffffffff) << 28)` (64 bits,
  identical to `283cdacf`);
  `w3 = ((k-1) >> 36) | ((c' >> 15) << 9) | (u << 42)` (`9 + 33 + 21 = 63` bits,
  `ystage.s` `lsl #9`).
  Phase-0 mining records pack `g < 2^58` (58 bits) and `idx < 2^56` (56 bits) in 114 bits.
- **Explicit `TBZ` block ends (`EXC_ENTRY = 0`):**
  `EM` (`kernels/leaf.s`) appends `tbnz \s, #12, gf\l` (`EMIT = 7`); `gf\l` runs
  `FLUSH7` + `b gn\l` (measured `17 <= 18` insns). `H2B` (`kernels/route.s`) appends
  `tbnz x11, #13, hf\l` (`H2S_ARR = 6`); `hf\l` resets `LC2[d2]` (`sub #8192; add #32`),
  checks `RELAY_CAP`, writes `HREADY`, and waits for ACK (measured `18 <= 18` insns).
- **128-bit cumulative fixup meter (`kernels/ytransport.s:ychunk`):**
  Because `Q1GD` scales by `256x` at `NODES = 2^58`, `ytransport.s:ychunk` natively
  maintains a persistent 128-bit cumulative fixup counter `(meter_lo, meter_hi)` at
  `YS+192..207` against the 128-bit threshold `DELTA_F - FIX_BATCH` at `YS+208..223`
  (`DELTA_F = 2^71`). Entry executes `ldp/ldp/subs/sbcs/b.lo 9f`, skipping leaf kernels
  whenever `meter > DELTA_F - FIX_BATCH`; post-`fwd` and post-`bwd` accumulate
  `x0 >> 32` and `x5 >> 32` via `adds/adc` (charged at `METER_CHUNK = 24` insns/chunk),
  giving `d_f <= 2^-21.8391`.
- **Production driver topology and horizon (`kernels/ystage.s`):**
  `yworker_sweep` loops `lo = 0 .. 2^58 - 1` (`hi = 0`, halting at `0x400000000000000`),
  `ybank_sweep` loops `owner = 0 .. 2^48 - 1` (`0x1000000000000`) with controller
  identity `(2^58, 0)`, `yboot_g` accepts `0 .. 2^58` inclusive (trapping `> 2^58`
  or `hi != 0`), `yboot_bank` accepts `< 2^58`, and `ycmd_gen` accepts `1 <= k <= 2^45`
  (`lsr x8, x0, #45`, covering `Q1GD = 27,136,000,000,000 < 2^45` and writing to `CN`
  at `0x5400_0000_0000 = 84 TiB`).

### 3.3 Exact Measured vs. Charged Kernel Counts (`experiment.py`)

| Kernel / Operation | Measured (AArch64 sim) | Charged in `certificate.py` | Slack |
| --- | ---: | ---: | ---: |
| `fwd` (4 Ys = 24 inc, 0 fixups) | `1,450` | `1,450` (`FWD_INC = 1433/24`) | `0` |
| `bwd` (4 Ys = 24 inc, 0 fixups) | `18,355` | `18,355` (`BWD_INC = 18331/24`) | `0` |
| `keys` (4 Ys = 24 inc, 144 leaves) | `3,525` | `<= 3,792` (`KEY_INC = 838/6`) | `+267` |
| `EM` 128-rec `TBZ` flush (`FLUSH1`) | `17` | `18` (`EXC_ENTRY = 0`) | `+1` |
| `h2s` full 128-rec block (`H2S_BLK128`) | `733` | `733` (SP table + `tbnz #13`) | `0` |
| `h2s` 256-rec `TBZ` flush (`FLUSH2`) | `18` | `18` (`EXC_ENTRY = 0`) | `0` |
| `dlvh` 256-rec non-`HCHGI` (`DLVHB`) | `264` | `272` (`+ BULK_H = 48`) | `+8` |
| `dlvh` 1-bucket collision (`MAIN_VISIT`) | `10` | `10` (unweighted: `12`) | `0` |
| `ychunk` 128-bit meter inline overhead | `15` (in `ychunk`) | `24` (`METER_CHUNK`) | `+9` |
| Assembled variants plus analytic code | `330,336 + 262,144 = 592,480 B` | `1,048,576 B` (`PROGRAM_BYTES`) | `+456,096 B` |

## 4. Complete 27-Row Expected Work Accounting

All 27 rows of `E_run` are summed in exact rationals by `certificate.py` over
`G = 22,000` groups (`m_g = floor(Q1GD * (g/22000)^(17/32))`), with `MSG = 256`
generic service, `GEN_MUL = 2^20`, `RESOLVE_RESERVE = 2^21`, `REQ_RESERVE = 4096`,
and `PAIR_INVERSE_EXC_WORK` included directly:

| Row | Component | `log2(insns)` upper | Share of `E_run` |
| --- | --- | ---: | ---: |
| 1 | Cold node init + program delivery (`MSG = 256`) | `112.0307` | `0.0000%` |
| 2 | Initial `R_g` (`GEN_MUL = 2^20`) | `78.0000` | `0.0000%` |
| 3 | Group & setup messages (`MSG = 256`) | `105.4253` | `0.0000%` |
| 4 | Phase-0 `X`-pair mining (full `MSG=256`, `GEN_MUL=2^20`) | `123.0906` | `0.0375%` |
| 5 | Phase 1-3 coin scalars | `100.2103` | `0.0000%` |
| 6 | Selected `X` setup (`MX = 6`) | `77.5850` | `0.0000%` |
| 7 | `EB` entry `X` & headers | `86.7549` | `0.0000%` |
| 8 | 5-window `P1`, `fwd`, `E`, `F`, `G` table build (`XTBLD_ALL`) | `123.6615` | `0.0556%` |
| 9 | 5-window `BT0`/`BT1` table build (`TBLD`) | `117.1997` | `0.0006%` |
| 10 | `LC1`/`LC2`/`RC2` init | `100.0001` | `0.0000%` |
| 11 | Home `HH` zeroing (`HB = 50`, `DC ZVA`) | `125.0224` | `0.1429%` |
| 12 | Home `R2`/`dlvh0` init | `86.3220` | `0.0000%` |
| 13 | Leaf generation + Y staging/transport + batch msg extra | `134.2999` | `88.6700%` |
| 14 | Fixup extra work | `107.5759` | `0.0000%` |
| 15 | Key canonicalization branch paths (`PJ_MEAN`) | `127.2362` | `0.6628%` |
| 16 | Routing + finals + scalar `card54_r` (`EMIT=7`, `H2S=733`, `BLOCKH=346`) | `130.9103` | `8.4612%` |
| 17 | Home joins (`NEW_INSERT=3`, `MAIN_VISIT=10`, `STOP_CHECK`) | `128.8062` | `1.9680%` |
| 18 | Same-`X` visits + `EXC_WORK` + `PAIR_INVERSE_EXC_WORK` | `90.3354` | `0.0000%` |
| 19 | `K2`/`d4` verification checks | `112.3494` | `0.0000%` |
| 20 | Candidate listings (`HCHG_LS1 = 6`) | `94.7611` | `0.0000%` |
| 21 | Defect rows (`DEFV + DEFK2 + DEFL + DEFJ + DEFWLT`) | `118.1045` | `0.0012%` |
| 22 | Candidate resolution (`RESOLVE` with `GEN_MUL=2^20`, `RESERVE=2^21`) | `115.0081` | `0.0001%` |
| 23 | Phase-4 generator requests (`REQ_SERVE`) | `113.1831` | `0.0000%` |
| 24 | Group exchanges | `105.4106` | `0.0000%` |
| 25 | Global barriers (`MSG = 256`) | `112.6732` | `0.0000%` |
| 26 | Phase control | `106.1257` | `0.0000%` |
| 27 | Final discrete-log output | `20.0000` | `0.0000%` |
| **Sum** | **`E_run`** | **`134.473293`** | **`100.0000%`** |
| **Gated** | **`E_total = 256 + p_gate * E_run` (`W = E_total / 1000`)** | **`133.115339`** | **`time_log2 <= 123.149554427177781`** |

## 5. Tighter Five-Window Error Bounds, Success Proof, and Verification Scope

Reducing `E` (6 -> 5 class-A terms), `G` (6 -> 5 terms: 3 class-A + 2 class-B), and
`BT` (6 -> 5 terms, `|y1 error| <= 5/2 + 2^-180`) strictly tightens the Lemma F/K/R
error bounds while keeping the exact same `EIS` classes (`X2A`, `X2B`, `X2EB`) and
`EISR` reduction kernel:
- `ERRH = [(-214, 210), (-120, 117), (-320, 319)]` (`WH = 640`, down from `722`),
- `ERRL = [(-324, 325), (-188, 183), (-483, 480)]` (`WL = 964`, down from `1,181`),
- `NWIN = 5,100` (down from `5,748`), `D_WJ <= 2^-50.6837`,
  `D_LAB_SYM <= 2^-8.8460` (down from `2^-8.65`), `KB_C132 / KB_BUCKET approximately 0.3335 < 0.34`,
- `P1_5WIN` drops one carry term in `S` and `Q`: `D_WLT <= 2^-23.9393`,
- Cap and defect tails at `Q1GD = 27,136,000,000,000` (`NODES = 2^58`):
  `d_f <= 2^-21.8391`, `d_p2 <= 2^-22.7722`, `d_cl <= 2^-18.1917`,
  `d_rq <= 2^-18.7893`, `d_rl <= 2^-17.2881`, `glob <= 2^-19.3983`,
  giving `lo_run >= 0.999653632361` and `p_gate * lo_run >= 0.39` in exact rationals.

Complete self-contained mathematical derivations appear in `success.md` (Section 7
contains self-contained proofs of five-window Lemmas E, F, M, K, R, and V, the
pointwise R3 joint-label and ordered-distinct affine-chain derivation, and the
global pair-inverse exclusion), `selected-x.md` (denominator-cleared selected-`X`
fiber bound `<= 30`), `fabric.md` (64-bit AXI fabric and scalar push/ACK services),
and `callbacks.md` (explicit scalar callback envelopes). Running `experiment.py`
inside Docker verifies the complete exact-rational `certificate.py`, expands and checks instruction forms in all7
AArch64 kernel files, executes `fwd`/`bwd`/`keys` on 4 real same-slice secp256k1
`Y` points x 6 real mined `X` points (`144/144` exact EC leaf matches), executes
`EM` 128-record and `h2s` 256-record `TBZ` flushes, `dlvh` 0-hit and 1-hit blocks,
`fin1`, `fin`, `ytransport.s:ychunk` (integrated 128-bit meter carry wrap and
threshold skip), `ystage.s:ycmd_gen` (accepting `k = Q1GD` and trapping `k > 2^45`),
`yboot_g`/`yboot_bank`/`yworker_sweep`/`ybank_sweep` `2^58`/`2^48` topology bounds,
sweeps 5,000 random 256-bit 5-window extraction vectors, and runs the explicitly listed
negative controls.

## 6. Full-horizon archive and explicit control charges

The22000 un-subdivided group intervals are the exact floor values specified
in certificate.py; the maximum length133854748188 fits1024 banks with2^27
chain slots each. Relative index r has bank=r>>27 and slot=r&(2^27-1).
B payload occupies[0,2^62), metadata begins2^62, payload/command inboxes
2^62+2^32, control2^62+2^33, and generic/code state2^62+2^34.
The complete backed intervals and all2^58 B nodes enter PEAK_BYTES.
The emitted producer, receiver, replay and bank sweep use these same
constants; absolute-k tags and stored coins remain unchanged.

Ordinary fwdq1 extracts w3 bits9..41 for the high cluster and keeps bits0..8
for high k. For c'=32768,m=0,k=1,u=0, w2=0,w3=512 and destination=2^25.
Mining uses its separate fwdq metadata decoder. callbacks.md supplies the
full ordinary inverse and reconstruction from the original CN.

PROGRAM_BYTES reserves1MiB per endpoint, including all measured assembly
variants plus262144 bytes for the analytic shared-loop construction. The
cold census charges its delivery to G, R, H, B and the single controller.
An explicit2^112 instruction protocol allowance dominates the enumerated
endpoint/state/helper setup and2^41 schedule construction bound. This
allowance is charged in addition to program delivery. Each entered group
pays owner/bank command work and metadata initialization at8 per slot and
4096 per command; the mathematical work bound uses the full population.

The selected-X proof, joint-label analysis and inverse-canonicality exclusions
are supplied in selected-x.md and success.md. The separate pair-inverse
batch term is retained on top of the original individual-denominator term.
Their union bounds also charge exceptional work; finite source probes do
not replace those production arguments.

### Additional credits

The concrete flat64 donor credits52-bit beta windows to soispoke's published
a4bddac1-3eb8-4215-bb53-42b58d2710c4, a rejected submission whose mechanism
is rederived here without relying on its verdict. The donor also cites
GoJackzi4d4e224b-3364-42dc-b7d4-64004d8858b4,
hosein-ul7f03684c-15b1-4b5c-ace4-28e465e828fc, and
kaankolcu94ae0ccc-a4f9-4c53-918f-c82619c90433. The stored-leaf base retains
its credits to kaankolcu e831f89f-27de-4fa9-becd-289bd7a4f4a8,
the grouped schedule in eed62949-dc73-4549-8271-372b70f6ce84, and
aburan28's relay reset fold5f726b85-b30c-49e1-b4fa-33b8a1ab0f55.
Public status and earlier scores do not serve as mathematical premises.

## 7. In-package source storage and scoped checks

The nine members named in runtime.json are the7 scalar assembly source
files, kernels/a64.py and experiment.py. They are stored as bounded zlib/base85 UTF8 data
with exact uncompressed lengths and SHA256 hashes. run.py restores them
in the credential-free experiment directory, rejects wrong hashes, lengths,
paths and trailing data, then runs experiment.py. All proof and numerical
files are plain text. The optional source input prints one restored kernel
without running the suite. Compression changes storage only and earns no
work credit. Every mathematical or algorithmic input is inside this package.

The default experiment also runs numeric_check.py's separate root inequalities,
320-bit integer logarithm interval, gate minimality,27-row sum, bank/cold
census and pair-map finite fibers. archive_check.py exercises actual emitted
archive/forwarding boundaries and old-layout controls; math_checks.py verifies
the GF(1024) and selected-X finite coefficient facts. No finite experiment
proves the production degree or collision argument. The source instruction
checker and functional emulator remain separate from those analytic premises.

The independent syntax/size receipt encoding-audit.json records69 LLVM
assemblies of expanded instruction lists, including64 generic role/PE
bindings. It records the exact kernel source hashes and the330336-byte
per-endpoint source census. The default Python experiment checks those
hashes and counts, and simulates finite paths; it does not invoke LLVM.
Raw macro argument quoting and the limits of parser/linked-program evidence
are disclosed in that receipt.

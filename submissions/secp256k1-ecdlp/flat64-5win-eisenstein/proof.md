# Flat 64-Bit Five-Window Eisenstein Stored-Leaf Birthday Search on secp256k1 (123.149553803893875)

## Claim

For uniform unrestricted production secp256k1 ECDLP challenges `Q = [x]P` and fresh
algorithm coins, this cold classical algorithm recovers `x` with probability at
least `0.39` (`Q = O => x = 0`) at expected classical work at most
`2^123.149553803893875` reference units (`C = 1,000` scalar AArch64 instructions per
unit), exact rational `time_log2 = 985196430431151/8000000000000` (`123.149553803893875`).
This improves on the live leader `283cdacf` (`123.186563709484764`) by at least
`0.037009905590889` bits.

All work across every processor is charged: cold initialization, Phase-0 `X`-pair
mining, table construction, coins, arithmetic, routing, home joins, candidate
verification, explicit scalar memory-fabric services (`MSG = 256`, `GEN_MUL = 2^20`),
exception ceilings, barriers, phase control, and output. Preprocessing is at most
`2^115.7840` units (included in `E_run`), advice is zero, and peak reserved memory
across all nodes and fabric switches is `348,384,021,872,999,421,293,902,116,968,190,801,494,560`
bytes (`<= 2^137.999731 B < 2^138 B`).

Summary figures (all produced directly by `certificate.py` and verified in Docker by
`experiment.py`):
- Gate (`route.s:gate`, `GATE_INSNS = 256`): `K_GATE = 132756105517980808685930450815710354843`,
  `p_gate = K_GATE / 2^128 <= 0.390135130183885`.
- Opened-run success lower bound: `lo_run >= 0.999653632361`; exact rational check
  `p_gate * lo_run >= 39/100` (and `(K_GATE - 1)/2^128 * lo_run < 39/100`).
- Opened-run expected work: `E_run <= 2^134.473293` instructions (`<= 2^124.507508`
  units, `<= 51.2251` insns per full-horizon leaf).
- Total expected work: `E_total = 256 + p_gate * E_run <= 2^133.115339` instructions
  (`<= 2^123.149553803893875` units).

## 1. Lineage, Credits, and Substantive Change

This submission builds directly on two accepted live-board designs:
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
2. **`4e2d5241` (`DeanMojo`, `db1d7ca7`, `123.324660124`)**, which established the
   uniform 64-bit flat analytical address space `[0, 2^64)` per processor with
   explicit `TBZ`/`TBNZ` block ends and 52-bit `BT` windows (building on `4d4e2241`
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
   analytical convention accepted in `4e2d5241`): only explicitly backed regions below
   count toward `PEAK_BYTES`; every scalar load/store targeting a backed address in
   `[0, 2^64)` counts as 1 instruction. No VMSAv9-128 MMU, page tables, `FEAT_D128`,
   `FEAT_LVA3`, `FEAT_S1PIE`, or stage-1 permission faults are used.
2. **AArch64 ISA extensions** (verified by `kernels/a64.py` encoding rules):
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

### 2.2 Complete Disjoint 64-Bit Flat Address Maps

Every table window base is aligned to its power-of-two size so `BFI`/`ORR` inserts
the exact index bits with zero mask overhead; every base is a multiple of `2^48`, so
materializing it takes a single `MOVZ Xd, #imm16, LSL #48`. `certificate.py:c40_maps`
verifies pairwise disjointness and `[0, 2^64)` containment of every region including
full 64-byte `LD64B`/`ST64B` access widths.

**Generator `G` (`NODES = 2^58` arithmetic + 1 control endpoint; `13,873,807,722,054,566,432 B` = `12.033 EiB < 16 EiB`):**

| Region | Address range `[lo, hi)` | Size | Alignment / Access |
| --- | --- | --- | --- |
| `LC1` | `[0x0000_0000_0000_0000, 0x0000_0800_0000_0000)` | `8 TiB` | `2^40 * 8 B` |
| `YR` | `[0x0000_2000_0000_0000, 0x0000_2100_0000_0000)` | `1 TiB` | `Y` scratch |
| `G34_CTRL_MMIO` | `[0x0000_2200_0000_0000, 0x0000_2600_0000_0000)` | `4 TiB` | control/stack/MMIO |
| `YQ_YVIEW` | `[0x0000_27ff_ffff_ffc0, 0x0000_2810_8000_0000)` | `66 GiB + 64 B` | 512 B Y frames |
| `EB` | `[0x0000_2a00_0000_0000, 0x0000_2a20_0000_0000)` | `128 GiB` | 512 B entries |
| `CN` | `[0x0000_5400_0000_0000, +64*Q1GD)` | `1,579.5 TiB` | `64 B * Q1GD` (`84 TiB` base) |
| `CNM` | `[0x0008_0000_0000_0000, +64*NBATCH_MINE)` | `90.6 TiB` | Phase-0 coins |
| `H1_CONTIG` | `[0x0100_0000_0000_0000, 0x0120_0000_0000_0000)` | `8 PiB` | `2^40 * 8 KiB` (`K1=128`) |
| `F_TBL_5WIN` | `[0x0800_0000_0000_0000, 0x0880_0000_0000_0000)` | `32 PiB` | `(45,45,45,44,44) * 256 B` |
| `FWD_TBL_5WIN` | `[0x0880_0000_0000_0000, 0x08c0_0000_0000_0000)` | `16 PiB` | `2 * (45,45,45,44,44) * 64 B` |
| `P1_TBL_5WIN` | `[0x1000_0000_0000_0000, 0x5800_0000_0000_0000)` | `4.500 EiB` | `6 * (2^52+4*2^51) * 64 B` |
| `E_TBL_5WIN` | `[0x6000_0000_0000_0000, 0x9000_0000_0000_0000)` | `3.000 EiB` | `(2^52+4*2^51) * 256 B` |
| `G_TBL_5WIN` | `[0xa000_0000_0000_0000, 0xd000_0000_0000_0000)` | `3.000 EiB` | `(2^52+4*2^51) * 256 B` |
| `G_BATCH_META` | `[0xd000_0000_0000_0000, +32*(Q1GD+...)+1 GiB)` | `790.5 TiB` | per-batch metadata |
| `REQ_64B` | `[0xd800_0000_0000_0000, +REQ_CAP*64+64)` | `1,000 PiB` | `64 B * (256*125*2^39)` |
| `BT0_52` | `[0xf800_0000_0000_0000, 0xfc00_0000_0000_0000)` | `0.250 EiB` | `2^52 * 64 B` (`win 0`) |
| `BT1_51x4` | `[0xfc00_0000_0000_0000, 2^64)` | `0.250 EiB` | `2^51 * 128 B` (`wins 1..4`) |

**Relay `R` (`N_REL = 2^40`; `107,804,916,093,403,136 B` = `0.0935 EiB`):**
`RCS_CTRL` `[0x8000, +1 MiB)`, `R_SLOTS_1024x8K` `[0xc000_0000_0000, +8 MiB)`,
`LC2` `[32 TiB, 64 TiB)`, `RC2` `[96 TiB, 128 TiB)`, `E2P` `[128 TiB, 192 TiB)`,
`R_SP_TBL` `[2^48 - 1 MiB, 2^48)`, `DPW_RCB` `[256 TiB, +32,288 TiB)`,
`H2U_CONTIG` `[1 EiB, 1.0625 EiB)` (`2^42 * 16 KiB` stride for `K2B = 256`).

**Home `H` (`N_HOME = 2^82`; `72,043,300,209,432,576 B` = `65,523 TiB < 2^56 B`):**
`H_CTRL_CODE` `[0x10000, +16 MiB)`, `H_HCC_CAND` `[0x7ffc_0000_0000_00, +CAND_CAP*528+1.02 MiB)`
(`1.00 TiB`), `H_WIN_BUF` `[0x7fff_0000_0000_00, +16 MiB)`, `H_SP_HDR` `[2^55 - 1 MiB, 2^55)`,
`R2_38220T` `[2^55, 2^55 + 38,220 TiB)`, `HH_2p50` `[0xff80_0000_0000_0000, +8,192 TiB)`,
`NX_19110T` `[0xffa0_0000_0000_0000, +19,110 TiB)`.

**Peak aggregate memory:**
`(2^58 + 1) * G_BYTES + 2^40 * R_BYTES + 2^82 * H_BYTES + FABRIC_BYTES = 348,384,021,872,999,421,293,902,116,968,190,801,494,560 B <= 2^137.999731 B < 2^138 B`.

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
  giving `d_f <= 2^-21.8392`.
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
| Total assembled program size | `319,224 B` (`79,806` insns) | `524,288 B` (`PROGRAM_BYTES`) | `+205,064 B` |

## 4. Complete 27-Row Expected Work Accounting

All 27 rows of `E_run` are summed in exact rationals by `certificate.py` over
`G = 22,000` groups (`m_g = floor(Q1GD * (g/22000)^(17/32))`), with `MSG = 256`
generic service, `GEN_MUL = 2^20`, `RESOLVE_RESERVE = 2^21`, `REQ_RESERVE = 4096`,
and `PAIR_INVERSE_EXC_WORK` included directly:

| Row | Component | `log2(insns)` upper | Share of `E_run` |
| --- | --- | ---: | ---: |
| 1 | Cold node init + program delivery (`MSG = 256`) | `105.8074` | `0.0000%` |
| 2 | Initial `R_g` (`GEN_MUL = 2^20`) | `78.0000` | `0.0000%` |
| 3 | Group & setup messages (`MSG = 256`) | `105.4253` | `0.0000%` |
| 4 | Phase-0 `X`-pair mining (full `MSG=256`, `GEN_MUL=2^20`) | `123.0906` | `0.0375%` |
| 5 | Phase 1-3 coin scalars | `100.2103` | `0.0000%` |
| 6 | Selected `X` setup (`MX = 6`) | `77.5850` | `0.0000%` |
| 7 | `EB` entry `X` & headers | `86.7549` | `0.0000%` |
| 8 | 5-window `P1`, `fwd`, `E`, `F`, `G` table build (`XTBLD_ALL`) | `123.6608` | `0.0556%` |
| 9 | 5-window `BT0`/`BT1` table build (`TBLD`) | `117.1997` | `0.0006%` |
| 10 | `LC1`/`LC2`/`RC2` init | `100.0001` | `0.0000%` |
| 11 | Home `HH` zeroing (`HB = 50`, `DC ZVA`) | `125.0224` | `0.1429%` |
| 12 | Home `R2`/`dlvh0` init | `86.3220` | `0.0000%` |
| 13 | Leaf generation + Y staging/transport + batch msg extra | `134.2999` | `88.6701%` |
| 14 | Fixup extra work | `107.5759` | `0.0000%` |
| 15 | Key canonicalization branch paths (`PJ_MEAN`) | `127.2362` | `0.6628%` |
| 16 | Routing + finals + scalar `card54_r` (`EMIT=7`, `H2S=733`, `BLOCKH=346`) | `130.9103` | `8.4612%` |
| 17 | Home joins (`NEW_INSERT=3`, `MAIN_VISIT=10`, `STOP_CHECK`) | `128.8062` | `1.9680%` |
| 18 | Same-`X` visits + `EXC_WORK` + `PAIR_INVERSE_EXC_WORK` | `88.3347` | `0.0000%` |
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
| **Gated** | **`E_total = 256 + p_gate * E_run` (`W = E_total / 1000`)** | **`133.115339`** | **`time_log2 <= 123.149553803893875`** |

## 5. Tighter Five-Window Error Bounds, Success Proof, and Verification Scope

Reducing `E` (6 -> 5 class-A terms), `G` (6 -> 5 terms: 3 class-A + 2 class-B), and
`BT` (6 -> 5 terms, `|y1 error| <= 5/2 + 2^-180`) strictly tightens the Lemma F/K/R
error bounds while keeping the exact same `EIS` classes (`X2A`, `X2B`, `X2EB`) and
`EISR` reduction kernel:
- `ERRH = [(-214, 210), (-120, 117), (-320, 319)]` (`WH = 640`, down from `722`),
- `ERRL = [(-324, 325), (-188, 183), (-483, 480)]` (`WL = 964`, down from `1,181`),
- `NWIN = 5,100` (down from `5,748`), `D_WJ <= 2^-50.6837`,
  `D_LAB_SYM <= 2^-8.8461` (down from `2^-8.65`), `KB_C132 / KB_BUCKET = 0.3335 < 0.34`,
- `P1_5WIN` drops one carry term in `S` and `Q`: `D_WLT <= 2^-23.9393`,
- Cap and defect tails at `Q1GD = 27,136,000,000,000` (`NODES = 2^58`):
  `d_f <= 2^-21.8392`, `d_p2 <= 2^-22.7722`, `d_cl <= 2^-18.1917`,
  `d_rq <= 2^-18.7894`, `d_rl <= 2^-17.2881`, `glob <= 2^-19.3983`,
  giving `lo_run >= 0.999653632361` and `p_gate * lo_run >= 0.39` in exact rationals.

Complete self-contained mathematical derivations appear in `success.md` (Section 7
contains self-contained proofs of five-window Lemmas E, F, M, K, R, and V, the
pointwise R3 joint-label and ordered-distinct affine-chain derivation, and the
global pair-inverse exclusion), `selected-x.md` (denominator-cleared selected-`X`
fiber bound `<= 30`), `fabric.md` (64-bit AXI fabric and scalar push/ACK services),
and `callbacks.md` (explicit scalar callback envelopes). Running `experiment.py`
inside Docker verifies the complete exact-rational `certificate.py`, assembles all 7
AArch64 kernel files, executes `fwd`/`bwd`/`keys` on 4 real same-slice secp256k1
`Y` points x 6 real mined `X` points (`144/144` exact EC leaf matches), executes
`EM` 128-record and `h2s` 256-record `TBZ` flushes, `dlvh` 0-hit and 1-hit blocks,
`fin1`, `fin`, `ytransport.s:ychunk` (integrated 128-bit meter carry wrap and
threshold skip), `ystage.s:ycmd_gen` (accepting `k = Q1GD` and trapping `k > 2^45`),
`yboot_g`/`yboot_bank`/`yworker_sweep`/`ybank_sweep` `2^58`/`2^48` topology bounds,
sweeps 5,000 random 256-bit 5-window extraction vectors, and runs seven deliberate
negative controls confirming each check fails on mutants.

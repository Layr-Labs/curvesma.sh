# Flat 64-Bit Five-Window Eisenstein Stored-Leaf Birthday Search on secp256k1 with 1024-Record Hop-1 Blocks (123.138763687151084)

## Claim

For uniform unrestricted production secp256k1 ECDLP challenges `Q = [x]P` and fresh
algorithm coins, this cold classical algorithm recovers `x` with probability at
least `0.39` (`Q = O => x = 0`) at expected classical work at most
`2^123.138763687151084` reference units (`C = 1,000` scalar AArch64 instructions per
unit), exact rational `time_log2 = 30784690921787771/250000000000000` (`123.138763687151084`).
This improves on the live leader `8fff60c5` (`123.148512898914260`) by at least
`0.009749211763176` bits.

All work across every processor is charged: cold initialization, Phase-0 `X`-pair
mining, table construction, coins, arithmetic, routing, home joins, candidate
verification, explicit scalar memory-fabric services (`MSG = 256`, `GEN_MUL = 912,688`),
exception ceilings, barriers, phase control, and output. Preprocessing is at most
`2^115.5480` units (included in `E_run`), advice is zero, and peak reserved memory
across all nodes and fabric switches is `336,726,598,255,798,309,662,388,641,207,582,680,042,016`
bytes (`<= 2^137.950630 B < 2^138 B`).

Summary figures (all produced directly by `certificate.py` and verified by `experiment.py`):
- Gate (`route.s:gate`, `GATE_INSNS = 256`): `K_GATE = 132730735379511136394949321660790148057`,
  `p_gate = K_GATE / 2^128 <= 0.390060574047759`.
- Opened-run success lower bound: `lo_run >= 0.999844706048`; exact rational check
  `p_gate * lo_run >= 39/100` (and `(K_GATE - 1)/2^128 * lo_run < 39/100`).
- Opened-run expected work: `E_run <= 2^134.462778` instructions (`<= 2^124.496994`
  units, `<= 50.8531` insns per full-horizon leaf).
- Total expected work: `E_total = 256 + p_gate * E_run <= 2^133.104548` instructions
  (`<= 2^123.138763687151084` units).

## 1. Lineage and Credits

This submission is **`8fff60c5` (`GoForceX`, `123.148512898914260`)** with the changes listed in Section 7.
`8fff60c5` is **`838b03f3` (`soispoke`, `123.149554`)** with byte-identical kernels and three charges set
to measured counts. `838b03f3` moved **`283cdacf` (`soispoke`, `123.186563709484764`)** onto a uniform
64-bit flat analytical address space. Their lineage, as declared there:
1. **`283cdacf`** introduced the explicit PULP AXI memory-fabric construction (`fabric.md`), counted
   scalar generic (`MSG = 256`) and callback services (`callbacks.md`), owner-shared `Y` transport
   (`Y_TRANSPORT = 129`, `OFFER = 576`), direct-`R2`-pointer home hashing (`NEW_INSERT = 3`,
   `MAIN_VISIT = 10`), direct-push `R -> H` delivery into `R2`, dedicated `SP` source-pointer tables on
   `R` and `H`, `JK` word-7 seed linking, `EB` static pointer pairs, the pointwise joint-label
   derivation, global `PAIR_INVERSE_EXC_WORK` charging, and the denominator-cleared selected-`X` fiber
   bound (`selected-x.md`), building on `4ce2d07e` (`bradcush`) and `d719e7cb` (`Jaydearcadian`). It also
   incorporates the `Q1GD` retune ratio (`106/103.079`) credited to `c4291c3b` (`aburan28`).
2. **`4e2d5241`** (`DeanMojo`, `db1d7ca7`, `123.324660124`) supplied an earlier concrete uniform 64-bit
   flat analytical address space `[0, 2^64)` per processor with explicit `TBZ`/`TBNZ` block ends and
   52-bit `BT` windows (building on `4d4e224b` by `GoJackzi`, `7f03684c` by `hosein-ul`, and
   `94ae0ccc` by `kaankolcu`).

**Architecture inherited from `8fff60c5`.** `283cdacf`'s Armv9.6-A VMSAv9-128 machine is replaced by
`4e2d5241`'s uniform **64-bit flat analytical address space** `[0, 2^64)` per processor at
`NODES = 2^58` (`N_CLUST = 2^48`, `K_S = 33` slice bits, `Q1GD = 27,136,000,000,000 < 2^45`), with
five-window `(52, 51, 51, 51, 51)` `BT`, `P1`, `E`, `G` tables (`BLK5F`, `P1_5WIN`, 5-window `JK`,
`OFFS5_52`), five-window `(45, 45, 45, 44, 44)` `F` and `fwd`, explicit `tbnz` block ends (`EMIT = 7`,
`H2S_ARR = 6`, `EXC_ENTRY = 0`), `kernels/ystage.s` sweeps and checks at `NODES = 2^58`, and a 128-bit
cumulative fixup meter in `kernels/ytransport.s:ychunk`. Relays (`N_REL = 2^40`) and homes
(`N_HOME = 2^82`, `HB = 50`, `R2`, `NX`, `HH`, `dlvh`) keep the `283cdacf` direct-`R2` push architecture.
This submission's blocks hold `K1 = 1,024` records at hop 1 and `K2B = 512` at hop 2.

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
   base AArch64 scalar integer instructions (including `ADR` for the hop-2 trampolines), `FEAT_LSE`
   (`LDADD`, `STADD`, `STADDL`, `SWP`), `FEAT_LRCPC2` (`STLUR`), `FEAT_LS64`
   (`LD64B`/`ST64B`, restricted to 64-byte-aligned backed addresses and even base registers `<= 22`),
   `DC ZVA` (2,048-byte zeroing granule on `H`), `FEAT_RNG` (`RNDR` for fresh unbiased coins),
   `FEAT_CMPBR` (`CBHS`/`CBHI`/`CBEQ`), and `WFE`/event wake semantics without
   spurious unrelated wakeups (two-publication control waits charged explicitly as in
   `283cdacf`).
3. **Explicit scalar memory fabric (`fabric.md`)**: every cross-node transfer is a
   counted scalar store/load/ACK over the 64-bit-address, 84-bit-tag PULP AXI tree.
4. **Heuristic / randomness premises**: uniform secp256k1 challenge `Q = [x]P` and
   independent fresh `RNDR` coins; finite experiments in `experiment.py` check
   assembly, 5-window arithmetic, exact EC labels on finite samples, hop-1 delivery,
   routing flushes, meter behavior, and map disjointness, never production runtime or RTL synthesis.

### 2.2 Backed data and service maps

The following exact intervals are generated from certificate.py. Each end is
exclusive. Scalar accesses to device apertures and the read-only owner alias
use the separate decode definitions in fabric.md; these are not additional
RAM arrays. Programs are placed within the listed backing at G 34TiB+1MiB,
R/H 0x10000, and B generic-base+1MiB, respectively. Each role is additionally
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

**G**: 14945664433369270816 charged bytes per endpoint.

| Region | Begin | End | Backed bytes |
|---|---:|---:|---:|
| `LC1` | `0x0` | `0x80000000000` | 8796093022208 |
| `YR` | `0x200000000000` | `0x210000000000` | 1099511627776 |
| `G34_CTRL_MMIO` | `0x220000000000` | `0x260000000000` | 4398046511104 |
| `YQ_YVIEW` | `0x27ffffffffc0` | `0x281080000000` | 70866960448 |
| `EB` | `0x2a0000000000` | `0x2a2000000000` | 137438953472 |
| `CN` | `0x540000000000` | `0x67f85e9000000` | 1736704000000000 |
| `CNM` | `0x8000000000000` | `0x85a995c000000` | 99614720000000 |
| `H1_CONTIG` | `0x100000000000000` | `0x200000000000000` | 72057594037927936 |
| `F_TBL_5WIN` | `0x800000000000000` | `0x880000000000000` | 36028797018963968 |
| `FWD_TBL_5WIN` | `0x880000000000000` | `0x8c0000000000000` | 18014398509481984 |
| `P1_TBL_5WIN` | `0x1000000000000000` | `0x5800000000000000` | 5188146770730811392 |
| `E_TBL_5WIN` | `0x6000000000000000` | `0x9000000000000000` | 3458764513820540928 |
| `G_TBL_5WIN` | `0xa000000000000000` | `0xd000000000000000` | 3458764513820540928 |
| `G_BATCH_META` | `0xd000000000000000` | `0xd003174e15fa41a0` | 870049073742240 |
| `REQ_CURSOR_MONITOR` | `0xd7fffffffffff800` | `0xd800000000000000` | 2048 |
| `REQ_64B` | `0xd800000000000000` | `0xf5a0000000000040` | 2134706223373615168 |
| `BT0_52` | `0xf800000000000000` | `0xfc00000000000000` | 288230376151711744 |
| `BT1_51x4` | `0xfc00000000000000` | `0x10000000000000000` | 288230376151711744 |

**H**: 69632262740725760 charged bytes per endpoint.

| Region | Begin | End | Backed bytes |
|---|---:|---:|---:|
| `H_CTRL_CODE` | `0x10000` | `0x1010000` | 16777216 |
| `H_WIN_BUF` | `0x7fff0000000000` | `0x7fff0001000000` | 16777216 |
| `H_BULK_CTRL` | `0x7ffffffc000000` | `0x7ffffffc010000` | 65536 |
| `H_GENERIC` | `0x7ffffffe000000` | `0x7ffffffe200000` | 2097152 |
| `H_SP_HDR` | `0x7ffffffff00000` | `0x80000000000000` | 1048576 |
| `R2_32767T` | `0x80000000000000` | `0xffff0000000000` | 36027697507336192 |
| `H_HCC_CAND` | `0x4000000000000000` | `0x401763ac8b305000` | 6583517185069056 |
| `HH_2p50` | `0xff80000000000000` | `0xffa0000000000000` | 9007199254740992 |
| `NX_16383T512G` | `0xffa0000000000000` | `0xffdfff8000000000` | 18013848753668096 |

**R**: 4701898748546105344 charged bytes per endpoint.

| Region | Begin | End | Backed bytes |
|---|---:|---:|---:|
| `RCS_CTRL` | `0x8000` | `0x114000` | 1097728 |
| `LC2` | `0x200000000000` | `0x400000000000` | 35184372088832 |
| `RC2` | `0x600000000000` | `0x800000000000` | 35184372088832 |
| `E2P` | `0x800000000000` | `0xc00000000000` | 70368744177664 |
| `R_GENERIC` | `0xc00001400000` | `0xc00001600000` | 2097152 |
| `R_SLOT_TABLES` | `0xc00002000000` | `0xc00002800000` | 8388608 |
| `R_SLOTS_1024x64K` | `0xc40000000000` | `0xc40004000000` | 67108864 |
| `R_SP_TBL` | `0xfffffff00000` | `0x1000000000000` | 1048576 |
| `DPW_RCB` | `0x1000000000000` | `0x3f41000000010000` | 4557642822899007488 |
| `H2U_CONTIG` | `0x8000000000000000` | `0x8200000000000000` | 144115188075855872 |

The population is 2^58 arithmetic G plus one controller, 2^40 relays, 2^82 homes and 2^58 banks.
The total is `(2^58+1)*G_BYTES+2^40*R_BYTES+2^82*H_BYTES+2^58*B_BYTES+FABRIC_BYTES = 336726598255798309662388641207582680042016` bytes, at most `2^137.950630 < 2^138`.

## 3. Five-Window Kernels, Record Packing, and Measured Instruction Counts

### 3.1 Five-Window `(52, 51, 51, 51, 51)` Decomposition (`BT`, `P1`, `E`, `G`)

Every 256-bit integer `w = (a, q, c, d)` decomposes uniquely into five unsigned
digits at bit positions `(0, 52, 103, 154, 205)` of widths `(52, 51, 51, 51, 51)`:
- `w0 = a[51:0]` (`bfi dst, a, #6, #52` for 64-byte rows; `#8, #52` for 256-byte rows),
- `w1 = extr(q, a, #52)[50:0]`,
- `w2 = extr(c, q, #39)[50:0]`,
- `w3 = extr(d, c, #26)[50:0]`,
- `w4 = (d >> 13)[50:0]`.

All five digits fit legally inside 64-bit `BFI`/`ORR`/`BIC` addressing (`experiment.py`
tests 5,000 random 256-bit vectors against exact integer reconstruction and region
alignment):
1. **`BLK5F` (`kernels/leaf.s`, 30 insns, runs twice per incidence):**
   `x7 = 0xf800_0000_0000_0000` (`BT0`, 64-byte rows) is overwritten once by
   `bfi x7, \a, #6, #52; ld64b x20, [x7]`. `x8 = 0xfc00_0000_0000_0000` (`BT1`
   slots 0/1) and `x9 = 0xfc00_0000_0000_0040` (`BT1` slots 2/3) have top 6 bits
   equal to `1`, so `orr x10, x8/x9, idx, lsl #7` masks bits `[63:58]` automatically
   for each 51-bit digit (`lsl #7` = 128-byte stride).
2. **`P1_5WIN` (`kernels/leaf.s`):**
   Loads five window bases `x12..x16` per `X_j` (base `1 EiB + 0.75*j EiB`, aligned
   to `2^58` for `win 0` and `2^57` for `wins 1..4`) in one `ld64b x12, [x29]` and
   computes four 64-byte-aligned row addresses by `BFI`. Window 4 takes one `BIC`:
   its table is stored reversed (row `i` holds digit `2^51 - 1 - i`) and the `X` frame
   holds its base word `W4BASE + 2^57 - 64`, so `bic x16, x16, x6, lsr #7` keeps the
   base's top 7 bits, writes the complement of `d[63:13]` into bits `[56:6]` and clears
   `[5:0]`: the address of digit `d >> 13`'s row. The builder writes that table from its
   top row with stride `-64` at the same count (`table_check.py`).
3. **Five-window `E` and `G` in `JK` (`kernels/leaf.s`) and `OFFS5_52` (`kernels/secp_fe.s`):**
   `OFFS5_52` writes five 256-byte row pointers for `E` (at frame offsets `280..312`,
   `win 0` at `+64` for word-7 seed linking) and five for `G` (`472..504`). In `JK`,
   `E` uses `KLKA5E`/`KLKB5E` (`KEB` before `KWINA/KWINB x29, Z`) and `G` uses
   `KLKA5`/`KK5`.
4. **Five-window `(45, 45, 45, 44, 44)` `F` and `fwd` (`K_S = 33` slice bits):**
   With `N_CLUST = 2^48 = 2^33 (s) * 2^15 (c' mod 2^15)`, `yround` extracts
   `s = x(Y)^3 >> 223` (`lsr x14, x13, #31`), leaving 223 varying bits decomposed by
   `OFF5PAIR` into `(45, 45, 45, 44, 44)` at shifts `(0, 45, 90, 135, 179)`. `fwd1`
   sits at `fwd0 | 0x0020_0000_0000_0000` (`orr x0, \b, #0x20000000000000` in `LK1`).
5. **`K'` pick (`tbldk`'s `eisr3`, `kernels/secp_fe.s`, 43 insns):** on the box of a class's lattice
   coordinates `f`, the first valid candidate is constant on every cell cut by the thresholds
   `-C - l1_i`, `C - l1_i + 1` (`f1`) and `-l2_i`, `X2 - l2_i` (`f2`); `eisr3` reads it as
   `T[4 band + col]` with `band = [f1 >= t_a] + [f1 >= t_b]` (192-bit signed compares) and
   `col = [f2 >= S_1] + [f2 >= S_2] + [f2 >= S_3]` (128-bit), from 24 words stored below the class
   descriptor. `certificate.py` (`c10k`) derives them and checks `T` on every cell exactly
   (80/54/77 cells for `A`/`B`/`EB`); `table_check.py` runs `eisr3` at both corners of every cell.

### 3.2 Record Packing at `NODES = 2^58`, Explicit `TBZ`/`TBNZ` Block Ends, and 128-Bit Fixup Meter

- **Record `w2, w3` (127 bits):** With `NODES = 2^58`, `M_G = 1024`,
  `N_CLUST = 2^48`, `Q1GD = 27,136,000,000,000 < 2^45`, and `u < 2^21`:
  `w2 = (m << 3) | ((c' & 0x7fff) << 13) | (((k-1) & 0xfffffffff) << 28)` (bits `[2:0]` are zero:
  records carry no incidence index `j`; a resolve tries the pool's six `X`s, `callbacks.md`);
  `w3 = ((k-1) >> 36) | ((c' >> 15) << 9) | (u << 42)` (`9 + 33 + 21 = 63` bits,
  `ystage.s` `lsl #9`). `keys` no longer adds `j` into `w2` (`KEY_INC = 832/6`).
  Phase-0 mining records pack `g < 2^58` (58 bits) and `idx < 2^56` (56 bits) in 114 bits; `xemit` ORs the
  pool tag `pool << 49` (`pool < 2^14`) above `A[255:207]` in word 1, so `HCHGI` matches only same-pool records.
- **Explicit `TBZ`/`TBNZ` block ends (`EXC_ENTRY = 0`), `K1 = 1,024`, `K2B = 512`:**
  `EM` (`kernels/leaf.s`) appends `tbnz \s, #15, gf\l` (`EMIT = 7`): a hop-1 block is 64 KiB, and record
  1023 lands at offset 32768. `gf\l` runs `FLUSH7` (its cursor reset is `sub #32768; add #32`) + `b gn\l`
  (measured 18 straight, 26 at the wait bound). `H2B` (`kernels/route.s`) appends `tbnz x11, #14, hf\l`
  (`H2S_ARR = 6`): `hf\l` is a two-instruction trampoline, `adr x19, hn\l; b hfshared`, into one shared
  stub that resets `LC2[d2]` (`sub #16384; add #32`), adds the block to `RC2[d2]`, writes `HREADY`, waits
  for ACK and returns by `br x19` (measured 21 straight, 29 at the wait bound); once `RC2[d2]` has reached the
  hop-2 cap `CAP_H`, it leaves `RC2[d2]` at the cap and returns without a flush (measured 11). Every `TBZ`/`TBNZ`
  in `route.s` is in reach (farthest `31,272 B`, `h2s`'s entry `tbz x20, #15, 8f`; farthest trampoline `tbnz`
  `24,628 B`). Phase-0 `xemit` keeps `tbnz x5, #15` with the same
  1,024-record blocks.
- **128-bit cumulative fixup meter (`kernels/ytransport.s:ychunk`):**
  `ytransport.s:ychunk` maintains a persistent 128-bit cumulative fixup counter `(meter_lo, meter_hi)` at
  `YS+192..207` against the 128-bit threshold `DELTA_F - FIX_BATCH` at `YS+208..223`
  (`DELTA_F = 2^100`). Entry executes `ldp/ldp/subs/sbcs/b.lo 9f`, skipping leaf kernels
  whenever `meter > DELTA_F - FIX_BATCH`; post-`fwd` and post-`bwd` accumulate
  `x0 >> 32` and `x5 >> 32` via `adds/adc` (charged at `METER_CHUNK = 15` insns/chunk),
  giving `d_f <= 2^-50.8391`.
- **Production driver topology and horizon (`kernels/ystage.s`):**
  `yworker_sweep` loops `lo = 0 .. 2^58 - 1` (`hi = 0`, halting at `0x400000000000000`),
  `ybank_sweep` loops `owner = 0 .. 2^48 - 1` (`0x1000000000000`) with controller
  identity `(2^58, 0)`, `yboot_g` accepts `0 .. 2^58` inclusive (trapping `> 2^58`
  or `hi != 0`), `yboot_bank` accepts `< 2^58`, and `ycmd_gen` accepts `1 <= k <= 2^45`
  (`lsr x8, x0, #45`, covering `Q1GD = 27,136,000,000,000 < 2^45` and writing to `CN`
  at `0x5400_0000_0000 = 84 TiB`).

### 3.3 Exact Measured vs. Charged Kernel Counts (`experiment.py`)

Waits are charged at the fabric's bound: a scalar wait is `LDR; CBNZ` (or `CBZ; BR`) when the ACK is
present, plus `LDR; CBNZ; WFE; B` for each failed pass, and one stale event plus one publication allow at
most two failed passes. The table gives each side's executed straight path and the executed count at that
bound.

| Kernel / Operation | Measured (AArch64 sim) | Charged in `certificate.py` | Slack |
| --- | ---: | ---: | ---: |
| `fwd` (4 Ys = 24 inc, 0 fixups) | `1,450` | `1,450` (`FWD_INC = 1433/24`) | `0` |
| `bwd` (4 Ys = 24 inc, 0 fixups) | `18,331` | `18,331` (`BWD_INC = 18307/24`) | `0` |
| `keys` (4 Ys = 24 inc, 144 leaves) | `3,501` | `<= 3,768` (`KEY_INC = 832/6`) | `+267` |
| `EM` block end (`FLUSH`, straight / wait bound) | `18 / 26` | `FLUSH + DONE_C1 = 18 + 8` | `0` |
| relay `bulk_r_receive` a block | `15 / 23` | `R_RECV = 23` | `0` |
| `h2s` full 1024-record block (`H2S_BLK`) | `5,773` | `5,773` | `0` |
| `h2s` partial path, `n < 1024` records | `<= 32 + 6n` | `H2S_PFIX = 32`, `6` a record | `0` |
| hop-2 block end (`FLUSH2`, trampoline) | `21 / 29` | `FLUSH2 + DONE_C = 21 + 8` | `0` |
| `dlvh` 512-rec non-`HCHGI` (`DLVHB`) | `520` | `520` | `0` |
| `dlvh` 1-bucket collision (`MAIN_VISIT`) | `10` | `10` (unweighted: `12`) | `0` |
| `bulk_h_receive` a block | `25 / 33` | `BULK_H = 33` | `0` |
| `fin1` a relay (group final, straight) | `22` | `FIN1 = 22` in `FINAL` | `0` |
| Phase-0 `xemit` block-full flush | `17` | `FLUSH = 18` | `+1` |
| relay member tables (1024 members) | `2 + 1024*2,059` | `RTBL_INIT = 2,108,418` | `0` |
| `zero16384` pass (512 `DC ZVA`) | `1,026` | `ZERO_PASS = 1,026` | `0` |
| `ychunk` 128-bit meter inline overhead | `15` (in `ychunk`) | `15` (`METER_CHUNK`) | `0` |
| `pfilt` + `ptick` (a mining candidate) | `<= 36 + 7 + 4` | `PF_CAND = 47 + 2*104` | `0` |
| `h2s`, the relay's last final in a group, then `h2sd` | `<= 38 + 6n` | `H2S_PFIX + H2S_LAST + H2S_DONE = 32 + 2 + 4` | `0` |
| `dlvh` segmented path (19 shapes, K2B 512) | `<= 1,565` | `DLVHS = 1,579` | `+14` |
| `init3` / `tinit` a word | `4 / 3` | `INIT3_BKT / TINIT_ENTRY = 4 / 3` | `0` |
| Assembled variants plus analytic code | `518,756 + 262,144 = 780,900 B` | `1,048,576 B` (`PROGRAM_BYTES`) | `+267,676 B` |

`FINAL` charges each `G`'s group-final block to one relay (`n < 1024` records): `fin1` (22) and its wait
(8), the relay's receive (23), `h2s`'s partial-path fixed part (32), and the records' excess of 6 each
over the full-block per-record shares they already pay,
`PDEF = ceil(1023 (6*1024 - (18 + 8 + 23 + 5773)) / 1024) = 322`: `FINAL = 407`.

## 4. Complete 27-Row Expected Work Accounting

All 27 rows of `E_run` are summed in exact rationals by `certificate.py` over
`G = 22,000` groups (`m_g = floor(Q1GD * (g/22000)^(17/32))`), with `MSG = 256`
generic service, `GEN_MUL = 912,688`, `CHAIN_ADD = 24,485`, six `RESOLVE_RESERVE = 2^21` a resolve,
six `REQ_RESERVE = 4096` a request, and `PAIR_INVERSE_EXC_WORK` included directly:

| Row | Component | `log2(insns)` upper | Share of `E_run` |
| --- | --- | ---: | ---: |
| 1 | Cold node init + program delivery (`MSG = 256`), `ksel` derivation | `112.0307` | `0.0000%` |
| 2 | Initial `R_g` (`GEN_MUL = 912,688`) | `77.7998` | `0.0000%` |
| 3 | Group & setup messages (`MSG = 256`) | `105.4253` | `0.0000%` |
| 4 | Phase-0 `X`-pair mining (`xemit` 1,024-record blocks, `MSG = 256`, `GEN_MUL = 912,688`) | `123.0661` | `0.0371%` |
| 5 | Phase 1-3 coin scalars | `100.2103` | `0.0000%` |
| 6 | Selected `X` setup (`MX = 6`) | `77.5850` | `0.0000%` |
| 7 | `EB` entry `X` & headers | `86.7549` | `0.0000%` |
| 8 | 5-window `P1`, `fwd`, `E`, `F`, `G` table build (`XTBLD_ALL`, `EISR3 = 43`) | `122.3515` | `0.0226%` |
| 9 | 5-window `BT0`/`BT1` table build (`TBLD`) | `117.1997` | `0.0006%` |
| 10 | `LC1`/`LC2`/`RC2` init + relay member tables (`RTBL_INIT`) | `100.0001` | `0.0000%` |
| 11 | Home `HH` zeroing (`HB = 50`, `DC ZVA`, 512 a pass) | `125.0029` | `0.1420%` |
| 12 | Home `R2`/`dlvh0` init | `86.3220` | `0.0000%` |
| 13 | Leaf generation + Y staging/transport + batch msg extra | `134.2965` | `89.1089%` |
| 14 | Fixup extra work | `107.5759` | `0.0000%` |
| 15 | Key canonicalization branch paths (`PJ_MEAN`) | `127.2359` | `0.6675%` |
| 16 | Routing + finals + scalar `card54_r` (`EMIT=7`, `K1=1024`: `H2S_BLK=5773`, `FINAL=407`; `BLOCKH=582` at `K2B=512`) | `130.8257` | `8.0377%` |
| 17 | Home joins (`NEW_INSERT=3`, `MAIN_VISIT=10`, `STOP_CHECK`) | `128.8058` | `1.9819%` |
| 18 | Same-`X` visits + `EXC_WORK` + `PAIR_INVERSE_EXC_WORK` | `90.1133` | `0.0000%` |
| 19 | `K2`/`d4` verification checks | `112.3483` | `0.0000%` |
| 20 | Candidate listings (`HCHG_LS1 = 6`) | `94.7611` | `0.0000%` |
| 21 | Defect rows (`DEFV + DEFK2 + DEFL + DEFJ + DEFWLT`) | `118.1046` | `0.0012%` |
| 22 | Candidate resolution (`RESOLVE`: six pool `X`s, six `2^21` reserves) | `116.5572` | `0.0004%` |
| 23 | Phase-4 generator requests (`REQ_SERVE`) | `113.0158` | `0.0000%` |
| 24 | Group exchanges | `105.4097` | `0.0000%` |
| 25 | Global barriers (`MSG = 256`) | `112.6732` | `0.0000%` |
| 26 | Phase control | `106.1257` | `0.0000%` |
| 27 | Final discrete-log output | `20.0000` | `0.0000%` |
| **Sum** | **`E_run`** | **`134.462778`** | **`100.0000%`** |
| **Gated** | **`E_total = 256 + p_gate * E_run` (`W = E_total / 1000`)** | **`133.104548`** | **`time_log2 <= 123.138763687151084`** |

## 5. Tighter Five-Window Error Bounds, Success Proof, and Verification Scope

Reducing `E` (6 -> 5 class-A terms), `G` (6 -> 5 terms: 3 class-A + 2 class-B), and
`BT` (6 -> 5 terms, `|y1 error| <= 5/2 + 2^-180`) keeps the Lemma F/K/R error bounds
with the same `EIS` classes (`X2A`, `X2B`, `X2EB`) and `EISR` reduction:
- `ERRH = [[-214, 210], [-120, 117], [-320, 319]]` (`WH = 640`),
- `ERRL = [[-324, 325], [-188, 183], [-483, 480]]` (`WL = 964`),
- `NWIN = 5,100`, `D_WJ <= 2^-50.6837`,
  `D_LAB_SYM <= 2^-8.8460`, `KB_C132 / KB_BUCKET approximately 0.3335 < 0.34`,
- `P1_5WIN` drops one carry term in `S` and `Q`: `D_WLT <= 2^-23.9393`,
- Cap and defect tails at `Q1GD = 27,136,000,000,000` (`NODES = 2^58`), with `RELAY_CAP`, `REQ_CAP`
  and `CAND_CAP` sized on the flat space:
  `d_f <= 2^-50.8391`, `d_p2 <= 2^-22.5351`, `d_cl <= 2^-28.1917`,
  `d_rq <= 2^-19.7123`, `d_rl <= 2^-21.8495`, `glob <= 2^-19.4303`,
  giving `lo_run >= 0.999844706048` and `p_gate * lo_run >= 0.39` in exact rationals.

Complete self-contained mathematical derivations appear in `success.md` (Section 7
contains self-contained proofs of five-window Lemmas E, F, M, K, R, and V, the
pointwise R3 joint-label and ordered-distinct affine-chain derivation, and the
global pair-inverse exclusion), `selected-x.md` (denominator-cleared selected-`X`
fiber bound `<= 30`), `fabric.md` (64-bit AXI fabric and scalar push/ACK services),
and `callbacks.md` (explicit scalar callback envelopes). Running `experiment.py`
verifies the complete exact-rational `certificate.py`, expands and checks instruction forms in all 7
AArch64 kernel files, executes `fwd`/`bwd`/`keys` on 4 real same-slice secp256k1
`Y` points x 6 real mined `X` points (`144/144` exact EC leaf matches), carries those 144 records
through modelled `H1` blocks and `fin1` into relay slots, where `bulk_r_receive` and `h2s` deliver them
into hop-2 windows that are checked against the exact labels, executes `EM`'s 1024-record block end
and its relay's full path, `h2s`'s full and partial paths with every partial-path chunk and all 2,048 `ADR`
trampolines (each paired statically with its own `tbnz` site and return point, and each run by a record that
fills its own hop-2 window), the stub's cap-clamp return, the handshake waits, the relay tables, Phase-0 `xemit`, `dlvh` 0-hit and 1-hit 512-record blocks,
`bulk_h_receive`, `fin`, `zero16384`, `ytransport.s:ychunk` (integrated 128-bit meter carry wrap and
threshold skip), `ystage.s:ycmd_gen` (accepting `k = Q1GD` and trapping `k > 2^45`),
`yboot_g`/`yboot_bank`/`yworker_sweep`/`ybank_sweep` `2^58`/`2^48` topology bounds, sweeps 5,000 random
256-bit 5-window extraction vectors, runs a finite model of `callbacks.md`'s six-`X` resolver (each of the 144
j-free records resolves to exactly one leaf of its `Y`, each `X` the right one for 24 of them, and six planted
collisions, one through each `X`, resolve to the planted pair and its mirror), runs the reply handle at
`DPW_RCB`'s first and last entries and at 2^56 - 64, 2^56 and 2^56 + 64, and one request's 48 reply messages
through `reqh`, `dlvq`, `fwdq1`, `sendr` and `fwdr`, runs `pfilt` and `ptick` against psel, the top index of each
bit-field insert (H2B at `d2 = 2^42 - 1`; `xemit`, `EM` and `fin1` at `d1 = 2^40 - 1`; `HCHGI` at `2^50 - 1`;
member 1023, pool `2^14 - 1`), `CN`'s boundary (under an LS64 range set to `CN`'s map interval, since `ycmd_gen`'s
inherited guard admits `k <= 2^45`), `dlvh`'s segmented path, the relay's last final with `h2sd`, `init3`,
`tinit` and the relay cursor's binding, checks that each map region an encoding carries fits that encoding's
range, and runs the explicitly listed negative controls.

## 6. Full-horizon archive and explicit control charges

The 22000 un-subdivided group intervals are the exact floor values specified
in certificate.py; the maximum length 133854748188 fits 998 of 1024 banks with 2^27
chain slots each. Relative index r has bank=r>>27 and slot=r&(2^27-1).
B payload occupies [0,2^62), metadata begins 2^62, payload/command inboxes
2^62+2^32, control 2^62+2^33, and generic/code state 2^62+2^34.
The complete backed intervals and all 2^58 B nodes enter PEAK_BYTES.
The emitted producer, receiver, replay and bank sweep use these same
constants; absolute-k tags and stored coins remain unchanged.

Ordinary fwdq1 extracts w3 bits 9..41 for the high cluster and keeps bits 0..8
for high k. For c'=32768,m=0,k=1,u=0, w2=0,w3=512 and destination=2^25.
Mining uses its separate fwdq metadata decoder. callbacks.md supplies the
full ordinary inverse and reconstruction from the original CN.

PROGRAM_BYTES reserves 1MiB per endpoint, including all measured assembly
variants plus 262144 bytes for the analytic shared-loop construction. The
cold census charges its delivery to G, R, H, B and the single controller.
An explicit 2^112 instruction protocol allowance dominates the enumerated
endpoint/state/helper setup and 2^41 schedule construction bound. This
allowance is charged in addition to program delivery; the relay member tables
(`RTBL_INIT` a relay) are charged on top of it, and so is the controller's one-time derivation of the three
24-word `ksel` tables (`callbacks.md`): `KSEL_BUILD = 2,403,328` instructions. Each entered group
pays owner/bank command work and metadata initialization at 8 per slot and
4096 per command; the mathematical work bound uses the full population.

The selected-X proof, joint-label analysis and inverse-canonicality exclusions
are supplied in selected-x.md and success.md. The separate pair-inverse
batch term is retained on top of the original individual-denominator term.
Their union bounds also charge exceptional work; finite source probes do
not replace those production arguments.

### Additional credits

The concrete flat64 donor credits 52-bit beta windows to soispoke's published
a4bddac1-3eb8-4215-bb53-42b58d2710c4, a rejected submission whose mechanism
is rederived here without relying on its verdict. The donor also cites
GoJackzi 4d4e224b-3364-42dc-b7d4-64004d8858b4,
hosein-ul 7f03684c-15b1-4b5c-ace4-28e465e828fc, and
kaankolcu 94ae0ccc-a4f9-4c53-918f-c82619c90433. The stored-leaf base retains
its credits to kaankolcu e831f89f-27de-4fa9-becd-289bd7a4f4a8,
the grouped schedule in eed62949-dc73-4549-8271-372b70f6ce84, and
aburan28's relay reset fold 5f726b85-b30c-49e1-b4fa-33b8a1ab0f55.
Public status and earlier scores do not serve as mathematical premises.

## 7. Changes from 8fff60c5

Built on `8fff60c5` (`GoForceX`), `838b03f3` and `283cdacf` (`soispoke`), and `4ce2d07e` (`bradcush`), the
stored-leaf design those build on. Every kernel count this package changes or adds is executed by `experiment.py` (Section 3.3 lists them,
with `dlvh`'s segmented path, the relay's last final and `h2sd`, `init3` and `tinit`, and Section 3.1 gives
`EISR3`, which table_check.py runs at both corners of every `ksel` cell); its other kernel charges
are 8fff60c5's, unchanged, and the figures before an arrow are 8fff60c5's. The rest are analytic charges or parameters: `GEN_MUL` and `CHAIN_ADD` are envelopes derived in
`callbacks.md` ("Conservative replacement that avoids the exact inherited scalar price"), the six resolve and
request reserves in its "Additional resolve work", `KSEL_BUILD` in its "Constructive larger static bound",
`PF_CAND`'s two canonical squares in its square bound (`97 <= 104`), and the caps are parameters whose tails
`success.md` bounds.
1. **Larger blocks.** Hop-1 blocks grow from 128 to 1,024 records (64-KiB `H1` stride at 64 PiB, relay
   slots of 64 KiB at 196 TiB, 512-pointer member tables at 192 TiB + 32 MiB) and hop-2 blocks from 256 to
   512 (32-KiB `H2U` windows). `EM` and `xemit` test bit 15, `H2B` bit 14. `h2s`'s full path covers 1,024
   records, its partial path gains three chunks, and `dlvh` covers 512. `H2B`'s block ends reach a shared
   stub through `ADR` trampolines (`FLUSH2` 18 -> 21; `ADR` added to `kernels/a64.py`). `FINAL` is
   restated for 1,024-record finals (Section 3.3).
2. **`P1_5WIN` window 4 by one `BIC`** with a reversed table (`bwd` 18,355 -> 18,331 per 4 Ys).
3. **Caps on the flat space:** `RELAY_CAP x128` (`H2U` moves to 8 EiB), `REQ_CAP` `256*237*2^39`,
   `CAND_CAP x1024` (`H_HCC_CAND` at `2^62`) with `R2` trimmed to 32,767 TiB.
4. **Smaller levers:** `ksel` (`EISR3 = 43`, `c10k`), `zero16384` unrolled to 512 `DC ZVA` a pass, `DELTA_F = 2^100`, and `callbacks.md`'s own envelopes
   charged (`CHAIN_ADD = 24,485`, `GEN_MUL = 912,688`).
5. **Records without `j`:** `keys` no longer adds `j` into `w2` (`-1` an incidence); a resolve tries the
   pool's six `X`s (`CAND_ITEM`, `RESOLVE`, `REQ_SERVE` x6, `callbacks.md`).
6. **Protocol charges at executed per-side worst cases:** the `G` hop-1 block end `18 + 8`, the relay's
   receive `15 + 8` (replacing `card54_r`'s 45 a block, set when it also covered `G`'s handshake, which `FLUSH` charges),
   Phase-0 blocks now paying the relay's receive, hop-2 `R` `21 + 8`, `H` `25 + 8`. `8fff60c5`'s certificate
   charges `BULK_H = 24`; its `bulk_h_receive` executes 25 instructions a block with no wait and 33 at the wait
   bound, so this submission charges 33.
7. **Package:** two runtime archives, the audit regenerated with LLVM `llvm-mc` (including the 2,048 `ADR`s),
   and the tests in Section 5.
8. **Reply handle.** A relay forwards each request with a handle naming the request's entry in its buffer
   (`DPW_RCB`); `fwdr` reads that entry's home word and reply slot to route each reply message. `8fff60c5`'s handle
   was the entry address shifted left 8 (`lsl x5,x1,#8` in `fwdq`/`fwdq1`, `lsr x7,x2,#8` in `fwdr`) with a byte
   offset in its low 8 bits, so it carries entry addresses below 2^56 and 8 messages (one 256-B reply) a request.
   With `RELAY_CAP x128` the buffer's entries reach `0x3f41_0000_0000_ffc0`, and a resolve's six replies are 48
   messages. The handle is now the 64-B-aligned entry address itself with the message index in its low 6 bits:
   `fwdq`/`fwdq1` `mov x5,x1`; `sendr` adds 1 a message; `fwdr` takes `and x7,x2,#~63`, `and x2,x2,#63` and
   `add x9,x9,x2,lsl #5`; `REQPUT` steps the home's reply slot by 1,536 B. Each change replaces one instruction
   with one, and the five routines execute the counts they did. `certificate.py` sizes a home's candidate at
   3,088 B (`HCC_BYTES`) and a resolved pair at 98 relay-buffer entries (`RCB_PAIR` in `TH_RELAY`, `rl_cr`,
   `M_VG_MAX` and `v_rl`); Phase-0 mining replies stay one 32-B message.
9. **Phase-0 pre-ticket filter.** `8fff60c5`'s Phase-0 text tickets every coarse `A[255:103]` match before
   psel's close-pair test, and a collector keeps only its first 1,024 tickets, so its retained prefix need not
   hold the close pairs psel needs. This package filters first: once a mining candidate's two replies are in
   (each reply's first word is `A = y^2 mod p`, canonical), the home runs `pfilt` (`kernels/route.s`) and sends
   the ticket only if psel's per-report test holds (`min(A) != 0`, `0 < |A(V) - A(U)| < 2^63`, distinct
   identities); `ptick` drives it over a home's candidates. It executes at most 36 instructions, plus 7 a
   candidate and 4 a call, charged as `PF_CAND = 47 + 2*104` (one canonical square a reply) a coarse candidate in `ROUTE_MINE` and
   `W_MINE_DET_MAX`. `success.md` Section 6.1 bounds the ticket count (`P_TKT_POOL` in `P_FAIL_POOL`,
   `c19b_tickets`).

## 8. In-package source storage and scoped checks

The nine members named in runtime.json (the 7 scalar assembly source files and
kernels/a64.py) and runtime2.json (experiment.py) are stored as bounded zlib/base85 UTF8 data
with exact uncompressed lengths and SHA256 hashes, two archives so that every file stays within
65,536 bytes. run.py restores them in the credential-free experiment directory, rejects wrong
hashes, lengths, paths and trailing data in either archive, then runs experiment.py. All proof and
numerical files are plain text. The optional source input prints one restored source without
running the suite. Compression changes storage only and earns no work credit. Every mathematical
or algorithmic input is inside this package.

Section 9 of the experiment lives in encoding_checks.py, which experiment.py runs in its own namespace
(an archive member is capped at 65,536 bytes). The default experiment also runs numeric_check.py's separate root inequalities,
320-bit integer logarithm interval, gate minimality, 27-row sum, bank/cold
census and pair-map finite fibers. archive_check.py exercises actual emitted
archive/forwarding boundaries and old-layout controls; math_checks.py verifies
the GF(1024) and selected-X finite coefficient facts; table_check.py runs the table builders,
`eisr3` on every `ksel` cell and the reversed `P1` window. No finite experiment
proves the production degree or collision argument. The source instruction
checker and functional emulator remain separate from those analytic premises.

The independent syntax/size receipt encoding-audit.json records 69 LLVM
(`llvm-mc` 23.1.5) assemblies of expanded instruction lists, including 64 generic role/PE
bindings. It records the exact kernel source hashes and the 518756-byte
per-endpoint source census. The default Python experiment checks those
hashes and counts, and simulates finite paths; it does not invoke LLVM.
Raw macro argument quoting and the limits of parser/linked-program evidence
are disclosed in that receipt.

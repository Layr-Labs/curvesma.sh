# Flat 64-Bit Five-Window Eisenstein Stored-Leaf Birthday Search on secp256k1 with Reversed-BT2 BIC, 4096/16384-Record Blocks, 5B/3B+3A Classes, and Leaf-3 Tail Duplication (123.133359612279188)

## Claim

For uniform unrestricted production secp256k1 ECDLP challenges `Q = [x]P` and fresh
algorithm coins, this cold classical algorithm recovers `x` with probability at
least `0.39` (`Q = O => x = 0`) at expected classical work at most
`2^123.133359612279188` reference units (`C = 1,000` scalar AArch64 instructions per
unit), exact rational `time_log2 = 30783339903069797/250000000000000` (`123.133359612279188`).
This improves on the live leader `8cc64cb7` (`123.138763687151084`) by at least
`0.005404074871896` bits (`5.404 mbit`).

All work across every processor is charged: cold initialization, Phase-0 `X`-pair
mining, table construction, coins, arithmetic, routing, home joins, candidate
verification, explicit scalar memory-fabric services (`MSG = 256`, `GEN_MUL = 912,688`),
exception ceilings, barriers, phase control, and output. Preprocessing is at most
`2^115.5481` units (included in `E_run`), advice is zero, and peak reserved memory
across all nodes and fabric switches is `336,726,660,589,537,242,387,422,505,421,235,548,334,624`
bytes (`<= 2^137.950630 B < 2^138 B`).

Summary figures (all produced directly by `certificate.py` and verified by `experiment.py`):
- Gate (`route.s:gate`, `GATE_INSNS = 256`): `K_GATE = 132724076471668945600418552103243496781`,
  `p_gate = K_GATE / 2^128 <= 0.390041005276381`.
- Opened-run success lower bound: `lo_run >= 0.999894869319`; exact rational check
  `p_gate * lo_run >= 39/100` (and `(K_GATE - 1)/2^128 * lo_run < 39/100`).
- Opened-run expected work: `E_run <= 2^134.457447` instructions (`<= 2^124.491662`
  units, `<= 49.2710` insns per full-horizon leaf).
- Total expected work: `E_total = 256 + p_gate * E_run <= 2^133.099144` instructions
  (`<= 2^123.133359612279188` units).

## 1. Lineage and Credits

This submission builds directly on **`8cc64cb7` (`bradcush`, `123.138763687151084`)** with the five
self-contained optimizations listed in Section 7, staying strictly inside `8cc64cb7`'s accepted
uniform 64-bit flat analytical machine model without introducing any new machine premises or
unchecked-drop perturbations. Upstream lineage as declared in `8cc64cb7`:
1. **`8cc64cb7`** (`bradcush`, `123.138763687151084`) introduced 1,024-record hop-1 and 512-record
   hop-2 blocks with `ADR` trampolines, reversed `P1_5WIN` window-4 `BIC`, flat-space cap enlargements,
   `ksel` (`EISR3 = 43`), `j`-free records with six-`X` resolution, 64-B-aligned reply handles, and the
   Phase-0 pre-ticket filter `pfilt`/`ptick`, building on **`8fff60c5`** (`GoForceX`, `123.148512898914260`)
   and **`838b03f3`** (`soispoke`, `123.149554`).
2. **`283cdacf`** (`soispoke`, `123.186563709484764`) introduced the explicit PULP AXI memory-fabric
   construction (`fabric.md`), counted scalar generic (`MSG = 256`) and callback services (`callbacks.md`),
   owner-shared `Y` transport (`Y_TRANSPORT = 129`, `OFFER = 576`), direct-`R2`-pointer home hashing
   (`NEW_INSERT = 3`, `MAIN_VISIT = 10`), direct-push `R -> H` delivery into `R2`, dedicated `SP`
   source-pointer tables on `R` and `H`, `JK` word-7 seed linking, `EB` static pointer pairs, the
   pointwise joint-label derivation, global `PAIR_INVERSE_EXC_WORK` charging, and the denominator-cleared
   selected-`X` fiber bound (`selected-x.md`), building on `4ce2d07e` (`bradcush`) and `d719e7cb`
   (`Jaydearcadian`), and incorporating the `Q1GD` retune credited to `c4291c3b` (`aburan28`).
3. **`4e2d5241`** (`DeanMojo`, `db1d7ca7`, `123.324660124`) supplied an earlier concrete uniform 64-bit
   flat analytical address space `[0, 2^64)` per processor with explicit `TBZ`/`TBNZ` block ends and
   52-bit `BT` windows (building on `4d4e224b` by `GoJackzi`, `7f03684c` by `hosein-ul`, and
   `94ae0ccc` by `kaankolcu`).

**Architecture inherited from `8cc64cb7`.** Each processor operates in the uniform **64-bit flat
analytical address space** `[0, 2^64)` at `NODES = 2^58` (`N_CLUST = 2^48`, `K_S = 33` slice bits,
`Q1GD = 27,904,000,000,000 < 2^45`), with five-window `(52, 51, 51, 51, 51)` `BT`, `P1`, `E`, `G` tables
(`BLK5F`, `P1_5WIN`, 5-window `JK`, `OFFS5_52`), five-window `(45, 45, 45, 44, 44)` `F` and `fwd`, explicit
`tbnz` block ends (`EMIT = 7`, `H2S_ARR = 6`, `EXC_ENTRY = 0`), `kernels/ystage.s` sweeps and checks at
`NODES = 2^58`, and the 128-bit cumulative fixup meter in `kernels/ytransport.s:ychunk`. Relays
(`N_REL = 2^40`) and homes (`N_HOME = 2^82`, `HB = 50`, `R2`, `NX`, `HH`, `dlvh`) keep the direct-`R2`
push architecture. This submission enlarges blocks to `K1 = 4,096` records at hop 1 and `K2B = 16,384`
records at hop 2, applies reversed window-4 `BIC` to `BLK5F` (`BT2_51x2`), eliminates class `EB` with
two-class five-window Eisenstein tables (`3B+3A` / `5B`) and balanced `(21,21,24,22)` labels, duplicates
the Leaf-3 canonicalization tail in `keys`, and refines the schedule to `45,000` groups.

## 2. Analytical Machine Model, Premises, and Complete 64-Bit Flat Maps

### 2.1 Explicit Machine Premises

This is an **analytical machine entry** under `tracks/docs/COST_MODEL.md` ("work in
scalar AArch64 instruction equivalents, including arithmetic, loads/stores,
addressing and control... Count other work on the same scalar 64-bit instruction
basis... memory has no scalar contribution or tie-break"). We do **not** claim
Armv9.6-A 56-bit physical hardware addressing on `G` or `R`.

Disclosed premises:
1. **Uniform 64-bit flat local address space `[0, 2^64)` per processor** (same
   analytical convention disclosed in `4e2d5241` and `8cc64cb7`): only explicitly backed regions below
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
   `283cdacf` and `8cc64cb7`).
3. **Explicit scalar memory fabric (`fabric.md`)**: every cross-node transfer is a
   counted scalar store/load/ACK over the 64-bit-address, 84-bit-tag PULP AXI tree.
4. **Heuristic / randomness premises**: uniform secp256k1 challenge `Q = [x]P` and
   independent fresh `RNDR` coins; finite experiments in `experiment.py` check
   assembly, 5-window arithmetic, exact EC labels on finite samples, hop-1 delivery,
   routing flushes, meter behavior, and map disjointness, never production runtime or RTL synthesis.

### 2.2 Backed data and service maps

The following exact intervals are generated from `certificate.py`. Each end is
exclusive. Scalar accesses to device apertures and the read-only owner alias
use the separate decode definitions in `fabric.md`; these are not additional
RAM arrays. Programs are placed within the listed backing at G 34TiB+1MiB,
R/H 0x10000, and B generic-base+1MiB, respectively. Each role is additionally
charged a conservative `PROGRAM_BYTES+2MiB` reserve beyond its interval sum;
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

**G**: 15161910991483054624 charged bytes per endpoint.

| Region | Begin | End | Backed bytes |
|---|---:|---:|---:|
| `LC1` | `0x0` | `0x80000000000` | 8796093022208 |
| `YR` | `0x200000000000` | `0x210000000000` | 1099511627776 |
| `G34_CTRL_MMIO` | `0x220000000000` | `0x260000000000` | 4398046511104 |
| `YQ_YVIEW` | `0x27ffffffffc0` | `0x281080000000` | 70866960448 |
| `EB` | `0x2a0000000000` | `0x2a2000000000` | 137438953472 |
| `CN` | `0x540000000000` | `0x6ac3a00800000` | 1785856000000000 |
| `CNM` | `0x8000000000000` | `0x85a995c000000` | 99614720000000 |
| `H1_CONTIG` | `0x100000000000000` | `0x500000000000000` | 288230376151711744 |
| `F_TBL_5WIN` | `0x800000000000000` | `0x880000000000000` | 36028797018963968 |
| `FWD_TBL_5WIN` | `0x880000000000000` | `0x8c0000000000000` | 18014398509481984 |
| `P1_TBL_5WIN` | `0x1000000000000000` | `0x5800000000000000` | 5188146770730811392 |
| `E_TBL_5WIN` | `0x6000000000000000` | `0x9000000000000000` | 3458764513820540928 |
| `G_TBL_5WIN` | `0xa000000000000000` | `0xd000000000000000` | 3458764513820540928 |
| `G_BATCH_META` | `0xd000000000000000` | `0xd0032db34ec021a0` | 894673073742240 |
| `REQ_CURSOR_MONITOR` | `0xd7fffffffffff800` | `0xd800000000000000` | 2048 |
| `REQ_64B` | `0xd800000000000000` | `0xf5a0000000000040` | 2134706223373615168 |
| `BT0_52` | `0xf800000000000000` | `0xfc00000000000000` | 288230376151711744 |
| `BT1_51x2` | `0xfc00000000000000` | `0xfe00000000000000` | 144115188075855872 |
| `BT2_51x2` | `0xfe00000000000000` | `0x10000000000000000` | 144115188075855872 |

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

**R**: 9169469579124129792 charged bytes per endpoint.

| Region | Begin | End | Backed bytes |
|---|---:|---:|---:|
| `RCS_CTRL` | `0x8000` | `0x114000` | 1097728 |
| `LC2` | `0x200000000000` | `0x400000000000` | 35184372088832 |
| `RC2` | `0x600000000000` | `0x800000000000` | 35184372088832 |
| `E2P` | `0x800000000000` | `0xc00000000000` | 70368744177664 |
| `R_GENERIC` | `0xc00001400000` | `0xc00001600000` | 2097152 |
| `R_SLOT_TABLES` | `0xc00002000000` | `0xc00004000000` | 33554432 |
| `R_SLOTS_1024x256K` | `0xc40000000000` | `0xc40010000000` | 268435456 |
| `R_SP_TBL` | `0xfffffff00000` | `0x1000000000000` | 1048576 |
| `DPW_RCB` | `0x1000000000000` | `0x3f41000000010000` | 4557642822899007488 |
| `H2U_CONTIG` | `0x8000000000000000` | `0xc000000000000000` | 4611686018427387904 |

The population is `2^58` arithmetic G plus one controller, `2^40` relays, `2^82` homes and `2^58` banks.
The total is `(2^58+1)*G_BYTES+2^40*R_BYTES+2^82*H_BYTES+2^58*B_BYTES+FABRIC_BYTES = 336726660589537242387422505421235548334624` bytes, at most `2^137.950630 < 2^138`.

## 3. Five-Window Kernels, Record Packing, and Measured Instruction Counts

### 3.1 Five-Window `(52, 51, 51, 51, 51)` Decomposition (`BT`, `P1`, `E`, `G`)

Every 256-bit integer `w = (a, q, c, d)` decomposes uniquely into five unsigned
digits at bit positions `(0, 52, 103, 154, 205)` of widths `(52, 51, 51, 51, 51)`:
- `w0 = a[51:0]` (`bfi dst, a, #6, #52` for 64-byte rows; `#8, #52` for 256-byte rows),
- `w1 = extr(q, a, #52)[50:0]`,
- `w2 = extr(c, q, #39)[50:0]`,
- `w3 = extr(d, c, #26)[50:0]`,
- `w4 = (d >> 13)[50:0]`.

All five digits fit legally inside 64-bit `BFI`/`BIC` addressing (`experiment.py`
tests 5,000 random 256-bit vectors against exact integer reconstruction and region
alignment):
1. **`BLK5F` (`kernels/leaf.s`, 29 insns, runs twice per incidence — `-12` insns/Y):**
   `BT1_51x4` is split into two equal 128-PiB tables with 64-byte row stride:
   `BT1_51x2` (`[0xfc00_0000_0000_0000, 0xfe00_0000_0000_0000)`, window 1 at `+0`, window 2 at `+32`)
   and `BT2_51x2` (`[0xfe00_0000_0000_0000, 2^64)`, window 3 at `+0`, window 4 **reversed** at `+32`,
   where row `i` holds digit `2^51 - 1 - i`). `BTBO`/`BTBE` load `x8 = 0xfc00_0000_0000_0000` and
   `x9 = 0xffff_ffff_ffff_ffc0` (2 instructions, identical cost). Window 0 uses `bfi x7, \a, #6, #52`
   (`BT0_52`); windows 1 and 2 use `bfi x8, x10, #6, #51`; window 3 uses `orr x10, x8, #0x0200_0000_0000_0000`
   after `bfi` into `x8` (wait — `bfi x9` / `bfi x10` against `0xfe00_0000_0000_0000`, `2` insns for
   `extr + bfi` on `x8` for win 1, `extr + bfi` on `x8` for win 2, `extr + bfi` on `x10` for win 3), and
   window 4 collapses from `lsr + orr` (2 insns) to a **single `bic x10, x9, \d, lsr #7`** (`+32` loaded via
   `ldp [x10, #32]` / `ldp [x10, #48]`, exact 64-B-aligned reversed row!). `tbldf` builds window 4 from the
   top row with stride `-64` (verified in `table_check.py`; `TBLD` adds one extra `tbldf` pass on `2^51`
   rows, costing `< 0.0001%` of `E_run`, while saving 12 hot instructions per `Y` across `bwd`).
2. **`P1_5WIN` (`kernels/leaf.s`):**
   Loads five window bases `x12..x16` per `X_j` (base `1 EiB + 0.75*j EiB`, aligned
   to `2^58` for `win 0` and `2^57` for `wins 1..4`) in one `ld64b x12, [x29]` and
   computes four 64-byte-aligned row addresses by `BFI`, plus window 4 by a single `BIC`
   against reversed `W4BASE + 2^57 - 64` (`table_check.py`).
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
   descriptor. With five windows, class `EB` is eliminated: `E` uses `3B + 2A` (plus `A` for `b^4`),
   and `F`/`G` use `5B`. `certificate.py` (`c10k`) derives and checks `T` on every cell exactly
   (80/54 cells for `A`/`B`, `KSEL_BUILD = 1,536,000`); `table_check.py` runs `eisr3` at both corners of every cell.
6. **Leaf-3 tail duplication in `keys` (`kernels/leaf.s`):**
   In `8cc64cb7`, every non-minimal canonicalization path in `SLH`/`T3` jumped back to a shared `KZ`
   return site (`b r...`), adding 1 unconditional branch on `5/9` of leaves (on `ep0..ep5`, `em0..em4`).
   Duplicating the 2-leaf `MHB2` / `PH2` tail inline on `SLH e` (probability `1/3`) and `T3(e)` (probability
   `1/9`) eliminates the return branch on `4/9` of the first 11 half-blocks per incidence while keeping every
   `FEAT_CMPBR` `cbhs`/`cbhi` strictly within its `+-255`-instruction encoding reach (`max_cb = 254/255`),
   reducing `PJ_MEAN` from `25/18` (`1.3889`) to `1579/1458` (`1.0830`) + defect slack (exact rational in
   `certificate.py`).

### 3.2 Record Packing at `NODES = 2^58`, Explicit `TBZ`/`TBNZ` Block Ends, and 128-Bit Fixup Meter

- **Balanced `(21, 21, 24, 22)` canonical label packing and 127-bit `w2, w3`:**
  Because five-window `(3B+3A / 5B)` tightens `WH = 502` and `WL = 736`, shifting one label bit from
  `H(v)` to `L(u)` (`MH_A = 21, MH_B = 21, ML_V = 24, ML_U = 22`, `43+43+40+42 = 168` bits) keeps
  `KB_C132/KB_BUCKET <= 0.3333 < 0.34` and `c9_kb` strictly within `1 + 1/700`, while cutting
  `D_LAB_SYM` from `2^-8.8460` to `2^-9.4861`. In `EM`, `bfxil \p,\q,#21,#21` and `bfxil \q,\r,#22,#42`
  pack `d1` (40 bits), `d2` (42 bits), `d3` (64 bits), and `d4` (22 bits, masked in `HCHGI` by
  `tst x12, #0xfffffc0000000000`). With `NODES = 2^58`, `M_G = 1024`, `N_CLUST = 2^48`,
  `Q1GD = 27,904,000,000,000 < 2^45`, and `u < 2^21`:
  `w2 = (m << 3) | ((c' & 0x7fff) << 13) | (((k-1) & 0xfffffffff) << 28)` (bits `[2:0]` zero: `j`-free
  records, resolved over the pool's six `X`s, `callbacks.md`);
  `w3 = ((k-1) >> 36) | ((c' >> 15) << 9) | (u << 42)` (`9 + 33 + 21 = 63` bits).
  Phase-0 mining records pack `g < 2^58` (58 bits) and `idx < 2^56` (56 bits) in 114 bits; `xemit` ORs the
  pool tag `pool << 49` (`pool < 2^14`) above `A[255:207]` in word 1, so `HCHGI` matches only same-pool records.
- **Explicit `TBZ`/`TBNZ` block ends (`EXC_ENTRY = 0`), `K1 = 4,096`, `K2B = 16,384`:**
  `EM` (`kernels/leaf.s`) appends `tbnz \s, #17, gf\l` (`EMIT = 7`): a hop-1 block is 256 KiB, and record
  4095 lands at offset `131072`. `gf\l` runs `FLUSH7` (cursor reset `sub #131072; add #32`) + `b gn\l`
  (measured 18 straight, 26 at the wait bound). `H2B` (`kernels/route.s`) appends `tbnz x11, #19, hf\l`
  (`H2S_ARR = 6`): `hf\l` is a two-instruction trampoline, `adr x19, hn\l; b hfshared`, into one shared
  stub that resets `LC2[d2]` (`sub #524288; add #32`), adds the block to `RC2[d2]`, writes `HREADY`, waits
  for ACK and returns by `br x19` (measured 21 straight, 29 at the wait bound); once `RC2[d2]` reaches
  `CAP_H`, it leaves `RC2[d2]` at the cap and returns without a flush (measured 11). `h2s`'s full path
  processes 4,096 records via a 7-iteration loop of 512 records (`h2s_lp`: 32 `FG`s) plus 31 `FG`s, 7
  `H2PR`s, and `flast` (`H2S_BLK = 23,069` = `5.6321` insns/record), while its partial path handles
  `n < 4,096` with a `h2big` loop over 1,024-record chunks plus the 10-bit binary ladder (`2,176` total
  `ADR` trampolines; farthest `TBZ`/`TBNZ` in `route.s` is `< 32,768 B`). `dlvh` processes 16,384 records
  in a 32-iteration loop of 256 `DLVHPR` pairs (`DLVHB = 16,459`). Phase-0 `xemit` uses `tbnz x5, #17`
  with the same 4,096-record blocks.
- **128-bit cumulative fixup meter (`kernels/ytransport.s:ychunk`):**
  Maintains a persistent 128-bit cumulative fixup counter `(meter_lo, meter_hi)` at `YS+192..207` against
  `DELTA_F - FIX_BATCH` at `YS+208..223` (`DELTA_F = 2^100`), charged at `METER_CHUNK = 15` insns/chunk,
  giving `d_f <= 2^-50.7988`.
- **Production driver topology and horizon (`kernels/ystage.s`):**
  `yworker_sweep` loops `lo = 0 .. 2^58 - 1`, `ybank_sweep` loops `owner = 0 .. 2^48 - 1` with controller
  identity `(2^58, 0)`, `yboot_g` accepts `0 .. 2^58` inclusive, `yboot_bank` accepts `< 2^58`, and
  `ycmd_gen` accepts `1 <= k <= 2^45` (covering `Q1GD = 27,904,000,000,000 < 2^45` into `CN` at `84 TiB`).

### 3.3 Exact Measured vs. Charged Kernel Counts (`experiment.py`)

Waits are charged at the fabric's bound: a scalar wait is `LDR; CBNZ` (or `CBZ; BR`) when the ACK is
present, plus `LDR; CBNZ; WFE; B` for each failed pass, and one stale event plus one publication allow at
most two failed passes. The table gives each side's executed straight path and the executed count at that
bound.

| Kernel / Operation | Measured (AArch64 sim) | Charged in `certificate.py` | Slack |
| --- | ---: | ---: | ---: |
| `fwd` (4 Ys = 24 inc, 0 fixups) | `1,450` | `1,450` (`FWD_INC = 1433/24`) | `0` |
| `bwd` (4 Ys = 24 inc, 0 fixups) | `18,283` | `18,283` (`BWD_INC = 18259/24`) | `0` |
| `keys` (4 Ys = 24 inc, 144 leaves) | `3,492` (400-Y mean `878.68/Y`) | `<= 3,768` (`KEY_INC = 832/6`, `PJ_MEAN = 1579/1458 + slack`) | `+276` |
| `EM` block end (`FLUSH`, straight / wait bound) | `18 / 26` | `FLUSH + DONE_C1 = 18 + 8` | `0` |
| relay `bulk_r_receive` a block | `15 / 23` | `R_RECV = 23` | `0` |
| `h2s` full 4096-record block (`H2S_BLK`) | `23,069` | `23,069` | `0` |
| `h2s` partial path, `n < 4096` records | `<= 84 + 6n` | `H2S_PFIX = 84`, `6` a record | `0` |
| hop-2 block end (`FLUSH2`, trampoline) | `21 / 29` | `FLUSH2 + DONE_C = 21 + 8` | `0` |
| `dlvh` 16384-rec non-`HCHGI` (`DLVHB`) | `16,459` | `16,459` | `0` |
| `dlvh` 1-bucket collision (`MAIN_VISIT`) | `10` | `10` (unweighted: `12`) | `0` |
| `bulk_h_receive` a block | `25 / 33` | `BULK_H = 33` | `0` |
| `fin1` a relay (group final, straight) | `22` | `FIN1 = 22` in `FINAL` | `0` |
| Phase-0 `xemit` block-full flush | `17` | `FLUSH = 18` | `+1` |
| relay member tables (1024 members) | `2 + 1024*8,203` | `RTBL_INIT = 8,399,874` | `0` |
| `zero16384` pass (512 `DC ZVA`) | `1,026` | `ZERO_PASS = 1,026` | `0` |
| `ychunk` 128-bit meter inline overhead | `15` (in `ychunk`) | `15` (`METER_CHUNK`) | `0` |
| `pfilt` + `ptick` (a mining candidate) | `<= 36 + 7 + 4` | `PF_CAND = 47 + 2*104` | `0` |
| `h2s`, the relay's last final in a group, then `h2sd` | `<= 90 + 6n` | `H2S_PFIX + H2S_LAST + H2S_DONE = 84 + 2 + 4` | `0` |
| `dlvh` segmented path (19 shapes, K2B 16384) | `<= 49,181` | `DLVHS = 49,195` | `+14` |
| `init3` / `tinit` a word | `4 / 3` | `INIT3_BKT / TINIT_ENTRY = 4 / 3` | `0` |
| Assembled variants plus analytic code | `529,220 + 262,144 = 791,364 B` | `1,048,576 B` (`PROGRAM_BYTES`) | `+257,212 B` |

`FINAL` charges each `G`'s group-final block to one relay (`n < 4096` records): `fin1` (22) and its wait
(8), the relay's receive (23), `h2s`'s partial-path fixed part (84), and the records' excess of 6 each
over the full-block per-record shares they already pay,
`PDEF = ceil(4095 * (6*4096 - (18 + 8 + 23 + 23069)) / 4096) = 1458`: `FINAL = 1595`.

## 4. Complete 27-Row Expected Work Accounting

All 27 rows of `E_run` are summed in exact rationals by `certificate.py` over
`G = 45,000` groups (`m_g = floor(Q1GD * (g/45000)^(18/32))`), with `MSG = 256`
generic service, `GEN_MUL = 912,688`, `CHAIN_ADD = 24,485`, six `RESOLVE_RESERVE = 2^21` a resolve,
six `REQ_RESERVE = 4096` a request, and `PAIR_INVERSE_EXC_WORK` included directly:

| Row | Component | `log2(insns)` upper | Share of `E_run` |
| --- | --- | ---: | ---: |
| 1 | Cold node init + program delivery (`MSG = 256`), `ksel` derivation | `112.0307` | `0.0000%` |
| 2 | Initial `R_g` (`GEN_MUL = 912,688`) | `77.7998` | `0.0000%` |
| 3 | Group & setup messages (`MSG = 256`) | `106.4577` | `0.0000%` |
| 4 | Phase-0 `X`-pair mining (`xemit` 4,096-record blocks, `MSG = 256`, `GEN_MUL = 912,688`) | `123.0658` | `0.0372%` |
| 5 | Phase 1-3 coin scalars | `100.2505` | `0.0000%` |
| 6 | Selected `X` setup (`MX = 6`) | `77.5850` | `0.0000%` |
| 7 | `EB` entry `X` & headers | `86.7549` | `0.0000%` |
| 8 | 5-window `P1`, `fwd`, `E`, `F`, `G` table build (`XTBLD_ALL`, `EISR3 = 43`) | `122.3515` | `0.0227%` |
| 9 | 5-window `BT0`/`BT1`/`BT2` table build (`TBLD`) | `117.2574` | `0.0007%` |
| 10 | `LC1`/`LC2`/`RC2` init + relay member tables (`RTBL_INIT`) | `100.0001` | `0.0000%` |
| 11 | Home `HH` zeroing (`HB = 50`, `DC ZVA`, 512 a pass) | `125.0029` | `0.1425%` |
| 12 | Home `R2`/`dlvh0` init | `86.3220` | `0.0000%` |
| 13 | Leaf generation + Y staging/transport + batch msg extra | `134.2926` | `89.2020%` |
| 14 | Fixup extra work | `107.6162` | `0.0000%` |
| 15 | Key canonicalization branch paths (`PJ_MEAN`) | `127.1509` | `0.6317%` |
| 16 | Routing + finals + scalar `card54_r` (`EMIT=7`, `K1=4096`: `H2S_BLK=23069`, `FINAL=1595`; `BLOCKH=16521` at `K2B=16384`) | `130.8088` | `7.9733%` |
| 17 | Home joins (`NEW_INSERT=3`, `MAIN_VISIT=10`, `STOP_CHECK`) | `128.8050` | `1.9881%` |
| 18 | Same-`X` visits + `EXC_WORK` + `PAIR_INVERSE_EXC_WORK` | `90.1938` | `0.0000%` |
| 19 | `K2`/`d4` verification checks | `112.3470` | `0.0000%` |
| 20 | Candidate listings (`HCHG_LS1 = 6`) | `94.8416` | `0.0000%` |
| 21 | Defect rows (`DEFV + DEFK2 + DEFL + DEFJ + DEFWLT`) | `118.1846` | `0.0013%` |
| 22 | Candidate resolution (`RESOLVE`: six pool `X`s, six `2^21` reserves) | `116.6378` | `0.0004%` |
| 23 | Phase-4 generator requests (`REQ_SERVE`) | `113.0963` | `0.0000%` |
| 24 | Group exchanges | `106.5074` | `0.0000%` |
| 25 | Global barriers (`MSG = 256`) | `112.6732` | `0.0000%` |
| 26 | Phase control | `106.1257` | `0.0000%` |
| 27 | Final discrete-log output | `20.0000` | `0.0000%` |
| **Sum** | **`E_run`** | **`134.457447`** | **`100.0000%`** |
| **Gated** | **`E_total = 256 + p_gate * E_run` (`W = E_total / 1000`)** | **`133.099144`** | **`time_log2 <= 123.133359612279188`** |

## 5. Tighter Five-Window Error Bounds, Success Proof, and Verification Scope

Exploiting the five-window reduction to replace class `EB` with two-class Eisenstein assignments
(`E`: `3B + 2A` + `A` for `b^4`; `F`, `G`: `5B`) and rebalancing label bits to `(21, 21, 24, 22)`
strictly tightens every Lemma F/K/R error bound:
- `ERRH = [[-171, 167], [-96, 93], [-251, 250]]` (`WH = 502`),
- `ERRL = [[-253, 254], [-153, 148], [-369, 366]]` (`WL = 736`),
- `NWIN = 4,068`, `D_WJ <= 2^-51.0098`,
  `D_LAB_SYM <= 2^-9.4861`, `KB_C132 / KB_BUCKET approximately 0.3333 < 0.34`,
- `P1_5WIN` unchanged: `D_WLT <= 2^-23.9393`,
- Cap and defect tails at `Q1GD = 27,904,000,000,000` (`NODES = 2^58`) over `45,000` groups:
  `d_f <= 2^-50.7988`, `d_p2 <= 2^-22.4861`, `d_cl <= 2^-28.1112`,
  `d_rq <= 2^-19.6318`, `d_rl <= 2^-22.7310`, `glob <= 2^-19.4261`,
  giving `lo_run >= 0.999894869319` and `p_gate * lo_run >= 0.39` in exact rationals.

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
into hop-2 windows that are checked against the exact labels, executes `EM`'s 4,096-record block end
and its relay's full path, `h2s`'s full and partial paths with every partial-path chunk and all 2,176 `ADR`
trampolines (each paired statically with its own `tbnz` site and return point, and each run by a record that
fills its own hop-2 window), the stub's cap-clamp return, the handshake waits, the relay tables, Phase-0 `xemit`, `dlvh` 0-hit and 1-hit 16,384-record blocks,
`bulk_h_receive`, `fin`, `zero16384`, `ytransport.s:ychunk` (integrated 128-bit meter carry wrap and
threshold skip), `ystage.s:ycmd_gen` (accepting `k = Q1GD` and trapping `k > 2^45`),
`yboot_g`/`yboot_bank`/`yworker_sweep`/`ybank_sweep` `2^58`/`2^48` topology bounds, sweeps 5,000 random
256-bit 5-window extraction vectors (including reversed `BT2_51x2` window 4 and reversed `P1_5WIN` window 4),
runs a finite model of `callbacks.md`'s six-`X` resolver, runs the reply handle at `DPW_RCB`'s first and last
entries and at `2^56 - 64`, `2^56` and `2^56 + 64` plus one request's 48 reply messages end to end, runs
`pfilt` and `ptick` against `psel`, tests the top index of every bit-field insert (`d2 = 2^42 - 1`,
`d1 = 2^40 - 1`, `HH = 2^50 - 1`, member `1023`, pool `2^14 - 1`), `CN`'s boundary, `dlvh`'s segmented path,
the relay's last final with `h2sd`, `init3`, `tinit`, region-in-encoding checks, and all negative controls.

## 6. Full-horizon archive and explicit control charges

The 45,000 un-subdivided group intervals are the exact floor values specified
in `certificate.py`; the maximum length `67,334,001,488` fits `502` of `1,024` banks with `2^27`
chain slots each. Relative index `r` has `bank = r >> 27` and `slot = r & (2^27 - 1)`.
`B` payload occupies `[0, 2^62)`, metadata begins `2^62`, payload/command inboxes
`2^62 + 2^32`, control `2^62 + 2^33`, and generic/code state `2^62 + 2^34`.
The complete backed intervals and all `2^58` `B` nodes enter `PEAK_BYTES`.
The emitted producer, receiver, replay and bank sweep use these same
constants; absolute-`k` tags and stored coins remain unchanged.

Ordinary `fwdq1` extracts `w3` bits `9..41` for the high cluster and keeps bits `0..8`
for high `k`. For `c' = 32768, m = 0, k = 1, u = 0`, `w2 = 0, w3 = 512` and `destination = 2^25`.
Mining uses its separate `fwdq` metadata decoder. `callbacks.md` supplies the
full ordinary inverse and reconstruction from the original `CN`.

`PROGRAM_BYTES` reserves 1 MiB per endpoint, including all measured assembly
variants (`529,220 B`) plus `262,144` bytes for the analytic shared-loop construction (`791,364 B` total). The
cold census charges its delivery to `G`, `R`, `H`, `B` and the single controller.
An explicit `2^112` instruction protocol allowance dominates the enumerated
endpoint/state/helper setup and `2^41` schedule construction bound (`45001 * 46 * 16 * 65536 = 2,170,600,554,496 < 2^41`). This
allowance is charged in addition to program delivery; the relay member tables
(`RTBL_INIT = 8,399,874` a relay) are charged on top of it, and so is the controller's one-time derivation of the two
24-word `ksel` tables (`callbacks.md`): `KSEL_BUILD = 1,536,000` instructions. Each entered group
pays owner/bank command work and metadata initialization at 8 per slot and
4096 per command; the mathematical work bound uses the full population.

The selected-`X` proof, joint-label analysis and inverse-canonicality exclusions
are supplied in `selected-x.md` and `success.md`. The separate pair-inverse
batch term is retained on top of the original individual-denominator term.
Their union bounds also charge exceptional work; finite source probes do
not replace those production arguments.

### Additional credits

The concrete flat64 donor credits 52-bit beta windows to soispoke's published
`a4bddac1-3eb8-4215-bb53-42b58d2710c4`, a rejected submission whose mechanism
is rederived here without relying on its verdict. The donor also cites
GoJackzi `4d4e224b-3364-42dc-b7d4-64004d8858b4`,
hosein-ul `7f03684c-15b1-4b5c-ace4-28e465e828fc`, and
kaankolcu `94ae0ccc-a4f9-4c53-918f-c82619c90433`. The stored-leaf base retains
its credits to kaankolcu `e831f89f-27de-4fa9-becd-289bd7a4f4a8`,
the grouped schedule in `eed62949-dc73-4549-8271-372b70f6ce84`, and
aburan28's relay reset fold `5f726b85-b30c-49e1-b4fa-33b8a1ab0f55`.
Public status and earlier scores do not serve as mathematical premises.

## 7. Changes from `8cc64cb7`

Built directly on `8cc64cb7` (`bradcush`, `123.138763687151084`), staying 100% inside its disclosed
uniform 64-bit flat analytical machine model and verification architecture. Every changed or added kernel
count is executed directly by `experiment.py` (Section 3.3) and `table_check.py`:
1. **`BLK5F` reversed window-4 `BIC` (`-2.80 mbit`, `-12` insns/Y):** `BT1_51x4` is split into
   `BT1_51x2` (`0xfc00..0xfe00`, wins 1 & 2, 64-B stride) and `BT2_51x2` (`0xfe00..2^64`, win 3 at `+0`,
   reversed win 4 at `+32`, 64-B stride). `BTBO`/`BTBE` load `x8 = 0xfc00_0000_0000_0000` and
   `x9 = 0xffff_ffff_ffff_ffc0`, collapsing window 4 in `BLK5F` from `lsr + orr` to a single `bic`
   (`30 -> 29` insns; `bwd` drops `18,331 -> 18,283` per 4 Ys, `BWD_INC = 18259/24`). `tbldf` builds
   window 4 reversed with stride `-64` (`table_check.py`).
2. **Larger blocks (`K1 = 4,096`, `K2B = 16,384`, `-1.36 mbit`):** Hop-1 blocks grow from `1,024` to
   `4,096` records (`H1_CONTIG` 256-KiB stride `64..320 PiB`, 256-KiB relay slots `R_SLOTS_1024x256K`,
   32-KiB member tables `R_SLOT_TABLES`, `RTBL_INIT = 8,399,874`) and hop-2 blocks from `512` to `16,384`
   records (`H2U_CONTIG` 1-MiB stride `8..12 EiB`). `EM`/`xemit` test bit 17, `H2B` tests bit 19.
   `h2s`'s full path uses a 7-pass 512-record loop plus 31 `FG`s (`H2S_BLK = 23,069 = 5.6321`/rec) and
   its partial path adds a 1,024-record chunk loop `h2big` (`H2S_PFIX = 84`, `2,176` `ADR` trampolines,
   all `TBZ`/`TBNZ` `< 32 KiB`). `dlvh` uses a 32-pass loop of 256 `DLVHPR` pairs (`DLVHB = 16,459`,
   `BLOCKH = 16,521`).
3. **Two-class five-window Eisenstein tables (`E = 3B+2A+A`, `F = G = 5B`) + balanced `(21,21,24,22)`
   labels (`-0.69 mbit`):** Five windows fit `3*(X2B-1) + 3*(X2A-1) < 2^128` and `5*(X2B-1) < 2^128`,
   eliminating class `EB` (`KSEL_BUILD`: `2,403,328 -> 1,536,000`) and shrinking `WH: 640 -> 502`,
   `WL: 964 -> 736`, `NWIN: 5,100 -> 4,068`. Rebalancing canonical label widths from `(20,21,24,23)` to
   `(21,21,24,22)` (`EM` `bfxil #21,#21` and `#22,#42`, `HCHGI` `tst x12,#0xfffffc0000000000`) cuts
   `D_LAB_SYM` from `2^-8.8460` to `2^-9.4861` while keeping `KB_C132/KB_BUCKET = 0.3333 < 0.34`,
   raising `lo_run` from `0.9998447` to `0.9998949`.
4. **Leaf-3 tail duplication in `keys` (`-0.55 mbit`):** `MHB2` / `PH2` inline the 2-leaf tail on
   `SLH e` and `T3(e)` of `ep0..ep5` and `em0..em4`, dropping the return branch on `4/9` of those
   half-blocks (`PJ_MEAN`: `25/18 -> 1579/1458 + slack`, `max_cb = 254/255`).
5. **Schedule retune (`G = 45,000`, exponent `18/32`, `Q1GD = 256*109,000,000,000`):** Tightens early-stop
   discretization and reduces the largest group from `133.9B` (`998` banks) to `67.3B` (`502` banks),
   tightening `d_rl` (`2^-21.85 -> 2^-22.73`) while keeping `COLD_SCHEDULE = 2,170,600,554,496 < 2^41`
   and `PEAK_BYTES = 2^137.950630 < 2^138`.

## 8. In-package source storage and scoped checks

The nine members named in `runtime.json` (the 7 scalar assembly source files and
`kernels/a64.py`) and `runtime2.json` (`experiment.py`) are stored as bounded zlib/base85 UTF8 data
with exact uncompressed lengths and SHA256 hashes, two archives so that every file stays within
65,536 bytes. `run.py` restores them in the credential-free experiment directory, rejects wrong
hashes, lengths, paths and trailing data in either archive, then runs `experiment.py`. All proof and
numerical files are plain text. The optional `source` input prints one restored source without
running the suite. Compression changes storage only and earns no work credit. Every mathematical
or algorithmic input is inside this package.

Section 9 of the experiment lives in `encoding_checks.py`, which `experiment.py` runs in its own namespace
(an archive member is capped at 65,536 bytes). The default experiment also runs `numeric_check.py`'s separate root inequalities,
320-bit integer logarithm interval, gate minimality, 27-row sum, bank/cold
census and pair-map finite fibers. `archive_check.py` exercises actual emitted
archive/forwarding boundaries and old-layout controls; `math_checks.py` verifies
the GF(1024) and selected-`X` finite coefficient facts; `table_check.py` runs the table builders,
`eisr3` on every `ksel` cell and both reversed window-4 tables (`P1` and `BT2`). No finite experiment
proves the production degree or collision argument. The source instruction
checker and functional emulator remain separate from those analytic premises.

The independent syntax/size receipt `encoding-audit.json` records 69 Apple LLVM
(`clang 17.0.0 -cc1as`) assemblies of expanded instruction lists, including 64 generic role/PE
bindings, exact kernel source hashes, and the `529,220`-byte per-endpoint source census (`791,364 B`
including the `262,144 B` analytic reserve, `257,212 B` under `1 MiB`). Because Apple LLVM 17.0.0
predates Armv9.6-A `FEAT_CMPBR` (`cbhs`/`cbhi`/`cbeq`), those `6,616` `leaf.s` instructions are
verified natively for syntax, operands, and `+-255`-instruction reach by `kernels/a64.py` (and
executed in `experiment.py`) and replaced 1:1 by 4-byte `nop` placeholders solely during the
`-cc1as` size/syntax pass, as explicitly disclosed in `encoding-audit.json`. The default Python
experiment checks those hashes and counts, and simulates finite paths; it does not invoke `clang`.

# Flat 64-Bit Three-Window Eisenstein Stored-Leaf Birthday Search on secp256k1 with Shared Table Banks (123.037576339062319)

## Claim

For uniform unrestricted production secp256k1 ECDLP challenges `Q = [x]P` and fresh
algorithm coins, this cold classical algorithm recovers `x` with probability at
least `0.39` (`Q = O => x = 0`) at expected classical work at most
`2^123.037576339062319` reference units (`C = 1,000` scalar AArch64 instructions per
unit), exact rational `time_log2 = 123037576339062319/1000000000000000` (`123.037576339062319`).
This is at least `0.101187348088765` bits below `123.138763687151084`, the score the board lists for `8cc64cb7` (`bradcush`), the package this one builds on; a comparison of two numbers; every claim here rests on the files of this
package. The whole improvement is one lever, shared table banks (Section 3.4):
three table windows where `8cc64cb7` uses five, which removes 69.67 instructions an incidence
(`fwd` 1,450 -> 1,322 and `bwd` 18,331 -> 16,787 for 24 incidences); everything else of the algorithm (the probability analysis, routing, homes, Phase 0,
caps) is `8cc64cb7`'s, unchanged.

All work across every processor is charged: cold initialization, Phase-0 `X`-pair
mining, table construction, coins, arithmetic, routing, home joins, candidate
verification, explicit scalar memory-fabric services (`MSG = 256`, `GEN_MUL = 912,688`),
exception ceilings, barriers, phase control, output, and the build of every table bank (Rows 8 and 9). Preprocessing is at most
`2^115.3731` units (included in `E_run`), advice is zero, and peak reserved memory
across all nodes, table banks and fabric switches is `336,722,941,418,476,409,732,520,609,555,662,222,805,536`
bytes (`<= 2^137.950614 B < 2^138 B`).

Summary figures (all produced directly by `certificate.py` and verified by `experiment.py`):
- Gate (`route.s:gate`, `GATE_INSNS = 256`): `K_GATE = 132730735379511136394949321660790148057`,
  `p_gate = K_GATE / 2^128 <= 0.390060574047759`.
- Opened-run success lower bound: `lo_run >= 0.999844706048`; exact rational check
  `p_gate * lo_run >= 39/100` (and `(K_GATE - 1)/2^128 * lo_run < 39/100`).
- Opened-run expected work: `E_run <= 2^134.361591` instructions (`<= 2^124.395807`
  units, `<= 47.4086` insns per full-horizon leaf).
- Total expected work: `E_total = 256 + p_gate * E_run <= 2^133.003361` instructions
  (`<= 2^123.037576339062319` units).

## 1. Lineage and Credits

This submission is **`8cc64cb7` (`bradcush`, `123.138763687151084`)** with the shared table banks of Section 3.4 and the changes
listed in Section 7. `8cc64cb7` is **`8fff60c5` (`GoForceX`,
`123.148512898914260`)** with the changes Section 7.2 lists. `8fff60c5` is **`838b03f3` (`soispoke`, `123.149554`)** with byte-identical kernels and three charges set
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

**Architecture inherited from `8cc64cb7` and `8fff60c5`.** `283cdacf`'s Armv9.6-A VMSAv9-128 machine is replaced by
`4e2d5241`'s uniform **64-bit flat analytical address space** `[0, 2^64)` per processor at
`NODES = 2^58` (`N_CLUST = 2^48`, `K_S = 33` slice bits, `Q1GD = 27,136,000,000,000 < 2^45`), with
explicit `tbnz` block ends (`EMIT = 7`, `H2S_ARR = 6`, `EXC_ENTRY = 0`), `kernels/ystage.s` sweeps and checks at `NODES = 2^58`, and a 128-bit
cumulative fixup meter in `kernels/ytransport.s:ychunk`. Relays (`N_REL = 2^40`) and homes
(`N_HOME = 2^82`, `HB = 50`, `R2`, `NX`, `HH`, `dlvh`) keep the `283cdacf` direct-`R2` push architecture.
This submission's blocks hold `K1 = 1,024` records at hop 1 and `K2B = 512` at hop 2. **New here:** the lookup tables are not
private five-window tables of each generator (`8cc64cb7`: windows `(52, 51, 51, 51, 51)` for `BT`, `P1`, `E`, `G` and
`(45, 45, 45, 44, 44)` for `F` and `fwd`, capped by the `2^64`-byte local address space) but shared read-only banks, read through
a counted route-latch store and `LD64B`; the windows are `(86, 85, 85)` for `BT`, `P1`, `E`, `G` and `(45, 89, 89)` for `F` and `fwd`
(window 0 private), and the kernels `BLK3F`, `P1_3WIN`, the three-window `JK`, `LKP6_PAIR` and `LK1` and the producer macros
`OFFS3_86` and `OFF3FF` replace `BLK5F`, `P1_5WIN`, the five-window `JK`, `OFFS5_52` and `OFF5PAIR`.

## 2. Analytical Machine Model, Premises, and Complete 64-Bit Flat Maps

### 2.1 Explicit Machine Premises

This is an **analytical machine entry** under `tracks/docs/COST_MODEL.md` ("work in
scalar AArch64 instruction equivalents, including arithmetic, loads/stores,
addressing and control... Count other work on the same scalar 64-bit instruction
basis... memory has no scalar contribution or tie-break"). We do **not** claim
Armv9.6-A 56-bit physical hardware addressing on `G` or `R`.

Disclosed premises (1-4 are 8cc64cb7's, 5 is new):
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
   assembly, three-window arithmetic through the aperture model, exact EC labels on finite samples, hop-1 delivery,
   routing flushes, meter behavior, and map disjointness, never production runtime or RTL synthesis.
5. **Shared read-only table banks and the T aperture (this package; `fabric.md` Section 2.4, Section 3.4 below).** Each lookup
   table is kept once per pool in passive endpoints (`T`, `N_T = 2,216,615,441,596,416`, each backing the local range `[2^62, 2^63)`), and a
   generator reads a row by one counted store of a route latch to the local range `[2^62, 2^63)` (an aperture of its decode
   endpoint, not RAM) and a counted `LD64B` of the row from an address of the same range. The endpoint serves a generator's accesses
   in program order, one at a time (the bridges' `MaxTrans = 1`), captures the latch when a load enters it and adds the generator's
   pool register. This is an analytical convention in the sense of premise 1 and of 8cc64cb7's read-only `YVIEW` alias of the owner's
   memory (a counted load of 64 bytes through the decode, no software at the owner): reads are priced as loads, every bank row is
   built by counted instructions (Rows 8 and 9, nothing amortised over later targets), every bank byte counts toward `PEAK_BYTES`,
   and no timing, bandwidth or contention is modelled beyond the instruction count, as for 8cc64cb7's relays and homes (`fabric.md` Section 2.4 lists each step of a read and its price). 8cc64cb7's
   generators fit their tables, as five windows of 51 to 52 bits, into their `2^64`-byte local space (no window of a table of 64-byte
   rows could exceed 58 bits there); shared banks lift that limit and this package uses windows of 85 to 86 bits. The aperture is the only
   change of machine convention of this package: `LD64B` is allowed on it (premise 2 restricts `LD64B` to 64-byte-aligned backed addresses; the
   aperture's lines are the backed bytes of the bank the latch names), and its register `POOLREG` is a device register at `20 TiB` (`service_aliases`).

### 2.2 Backed data and service maps

The following exact intervals are generated from certificate.py. G's five-window table regions of 8cc64cb7 (`P1_TBL_5WIN`, `E_TBL_5WIN`,
`G_TBL_5WIN`, `BT0_52`, `BT1_51x4`, `F_TBL_5WIN`, `FWD_TBL_5WIN`) are gone: G keeps window 0 of `F` and `fwd` (`F_W0_PRIV`,
`FWD_W0_LINE0`, `FWD_W0_LINE1`) and reads every other table through the aperture `T_APERTURE` `[2^62, 2^63)` (not a RAM region: it
is listed in `service_aliases`, which `c_service_decode` checks disjoint from every backed region). Each end is
exclusive. Scalar accesses to device apertures and the read-only owner alias
use the separate decode definitions in fabric.md; these are not additional
RAM arrays. Programs are placed within the listed backing at G 34TiB+1MiB,
R/H/T 0x10000, and B generic-base+1MiB, respectively. Each role is additionally
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

**G**: 2222995486047619616 charged bytes per endpoint.

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
| `F_W0_PRIV` | `0x800000000000000` | `0x820000000000000` | 9007199254740992 |
| `FWD_W0_LINE0` | `0x880000000000000` | `0x888000000000000` | 2251799813685248 |
| `FWD_W0_LINE1` | `0x8a0000000000000` | `0x8a8000000000000` | 2251799813685248 |
| `G_BATCH_META` | `0xd000000000000000` | `0xd003174e15fa41a0` | 870049073742240 |
| `REQ_CURSOR_MONITOR` | `0xd7fffffffffff800` | `0xd800000000000000` | 2048 |
| `REQ_64B` | `0xd800000000000000` | `0xf5a0000000000040` | 2134706223373615168 |

**H**: 69632262740725760 charged bytes per endpoint.

| Region | Begin | End | Backed bytes |
|---|---:|---:|---:|
| `H_CTRL_CODE` | `0x10000` | `0x1010000` | 16777216 |
| `H_WIN_BUF` | `0x7fff0000000000` | `0x7fff0001000000` | 16777216 |
| `H_SP_HDR` | `0x7ffffffff00000` | `0x80000000000000` | 1048576 |
| `H_BULK_CTRL` | `0x7ffffffc000000` | `0x7ffffffc010000` | 65536 |
| `H_GENERIC` | `0x7ffffffe000000` | `0x7ffffffe200000` | 2097152 |
| `R2_32767T` | `0x80000000000000` | `0xffff0000000000` | 36027697507336192 |
| `H_HCC_CAND` | `0x4000000000000000` | `0x401763ac8b305000` | 6583517185069056 |
| `HH_2p50` | `0xff80000000000000` | `0xffa0000000000000` | 9007199254740992 |
| `NX_16383T512G` | `0xffa0000000000000` | `0xffdfff8000000000` | 18013848753668096 |

**R**: 4701898748546105344 charged bytes per endpoint.

| Region | Begin | End | Backed bytes |
|---|---:|---:|---:|
| `RCS_CTRL` | `0x8000` | `0x114000` | 1097728 |
| `R_SLOTS_1024x64K` | `0xc40000000000` | `0xc40004000000` | 67108864 |
| `R_SLOT_TABLES` | `0xc00002000000` | `0xc00002800000` | 8388608 |
| `R_GENERIC` | `0xc00001400000` | `0xc00001600000` | 2097152 |
| `LC2` | `0x200000000000` | `0x400000000000` | 35184372088832 |
| `RC2` | `0x600000000000` | `0x800000000000` | 35184372088832 |
| `E2P` | `0x800000000000` | `0xc00000000000` | 70368744177664 |
| `R_SP_TBL` | `0xfffffff00000` | `0x1000000000000` | 1048576 |
| `DPW_RCB` | `0x1000000000000` | `0x3f41000000010000` | 4557642822899007488 |
| `H2U_CONTIG` | `0x8000000000000000` | `0x8200000000000000` | 144115188075855872 |

**T**: 4611686018449408000 charged bytes per endpoint.

| Region | Begin | End | Backed bytes |
|---|---:|---:|---:|
| `T_CTRL_CODE` | `0x10000` | `0x1010000` | 16777216 |
| `T_GENERIC` | `0xc00001400000` | `0xc00001600000` | 2097152 |
| `T_BANK` | `0x4000000000000000` | `0x8000000000000000` | 4611686018427387904 |

The population is 2^58 arithmetic G plus one controller, 2^40 relays, 2^82 homes, 2^58 archive banks B and `N_T = 2,216,615,441,596,416`
table banks T (`T_MEM_REGIONS`; Section 3.4).
The total is `336,722,941,418,476,409,732,520,609,555,662,222,805,536` bytes, at most `2^137.950614 < 2^138`.

## 3. Three-Window Kernels, Record Packing, and Measured Instruction Counts

### 3.1 Three-Window `(86, 85, 85)` Decomposition (`BT`, `P1`, `E`, `G`) and `(45, 89, 89)` for `F` and `fwd`

Every 256-bit integer `w = (a, q, c, d)` decomposes uniquely into three unsigned digits at bit positions `(0, 86, 171)` of
widths `(86, 85, 85)`. A digit `d_k` names a row of a shared table: the low `R` bits of `d_k` select the row inside a bank and the
remaining high bits select the bank (`R = 56` for 64-byte rows, `54` for 256-byte rows, `55` for 128-byte rows; Section 3.4). With
`a, q, c, d` the four limbs of `w`, the kernels extract (`experiment.py` mirrors each extraction on 5,000 random 256-bit vectors and
compares it with the integer decomposition):
- row of window 0 `bfi xr, a, #6, #56`, bank of window 0 `extr x1, q, a, #56; bfxil xt, x1, #0, #30`;
- row of window 1 `extr x1, c, q, #22; bfi xr, x1, #6, #56`, bank `bfxil xt, c, #14, #29`;
- row of window 2 `extr x1, d, c, #43; bfi xr, x1, #6, #56`, bank `bfxil xt, d, #35, #29`
  (`#8, #54` and 32/31/31-bit banks for 256-byte `E` and `G` rows).
The pointer registers hold `0x4000_0000_0000_0000` (the aperture base `2^62`) before the `bfi`, and the tag registers hold
`FAM << 36` (`FAM` the table/window id of Section 3.4) before the `bfxil`.
1. **`BLK3F` (`kernels/leaf.s`, 23 instructions, runs twice per incidence):** three windows of `BT` (64-byte rows, one `LD64B` each);
   each window is `extr/bfi/bfxil` (address and tag), one `str` of the tag to the aperture address, one `ld64b`; the three rows are
   added (`A4`). `BLK5F` has 30. The `BT` prefix `FAM_BT[0] << 36` is `mov x11, #0x18000000000` (`BTBO`).
2. **`P1_3WIN` (`kernels/leaf.s`):** one `ld64b x12, [x29]` loads the three `P1` tag words and the aperture base words of the
   incidence's `X_j` frame (`e + 320 .. e + 376`; `FAM = 6 + 3 j + k`; the table of frame words in Section 3.4 lists every such word); each of the three windows is a tag `str` and an
   `ld64b` of a 64-byte row holding `2 y_X^2 d` and `2 y_X d` mod `p`; the sums keep `P1_5WIN`'s carry counters and its final two
   `madd`s. No reversed window and no `bic`: window 4 of `P1_5WIN` is gone.
3. **Three-window `E`, `F` and `G` in `JK` (`kernels/leaf.s`) and `OFFS3_86`, `OFF3FF` (`kernels/secp_fe.s`):** the Y producer
   (`WORDS`, once per `Y` and shared by the 1,024 members of its cluster) stores, in the 512-byte frame of `Y`, a row pointer
   and a tag word for each window: `OFFS3_86` for `E` (words `272 .. 312`) and `G` (`464 .. 504`), `OFF3FF` for `F` (`408 .. 440`) and
   `fwd` (`96 .. 120`, `368`). `KLKA3E`, `KLKB3E`, `KLKA3`, `KLKB3`, `KLKA3F`, `KLKB3F` (one per table-sum of `JK`) do for each shared
   window one tag `str` and two `LD64B` of the 256-byte row, the second line at `a + 64` (`KWINA`/`KWINB`: an `add` for the second line). Window 0 of `E`, `G` and `F`
   is entered at the row's second line (pointer `+64`), and word 7 of that line is `tbldk`'s link word, the address of the 128-byte half of the row that the line closes (the row base for
   the half at `+0`, `row + 128` for the half at `+128`), which yields that half's base without an address instruction, as in 8cc64cb7; the two passes of `JK` read the four lines of the row. Window 0 of `F` is a private table of the generator and needs no tag.
4. **`(45, 89, 89)` `F` and `fwd` (`K_S = 33` slice bits):** `yround` computes `B = x(Y)^3 mod p` (canonical) and extracts `s = B >> 223` (`lsr x14, x13, #31`) as in
   8cc64cb7, leaving the 223-bit `B_low = B mod 2^223`. Window 0 is its low 45 bits (a private table: its rows carry the cluster's slice
   constant `S_HIGH = s << 223`), windows 1 and 2 are bits `45 .. 133` and `134 .. 222` (89 bits each, shifts `45` and `134`);
   `fwd`'s shared-window rows are 128 bytes, two contiguous 64-byte lines at `+0` and `+64` (the first holds the constants for `A_1` and `A_0`, the
   second `0` and the constant for `A_2`); `LKP6_PAIR` and `LK1` read window 0 from the private table (`fwd0`, `fwd0 | 2^53` for its
   second line, as in 8cc64cb7) and windows 1 and 2 through the aperture.
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
| `fwd` (4 Ys = 24 inc, 0 fixups) | `1,322` | `1,322` (`FWD_INC = 1305/24`; 8cc64cb7: `1,450`) | `0` |
| `bwd` (4 Ys = 24 inc, 0 fixups) | `16,787` | `16,787` (`BWD_INC = 16763/24`; 8cc64cb7: `18,331`) | `0` |
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
| `OFFS3_86` / `OFF3FF` / `OFFS3_86` (the Y producer's table words) | `26 / 42 / 26` | `YR_WORDS_NEW = 94` (8cc64cb7's `OFFS5_52`, `OFF5PAIR`, `OFFS5_52`: `YR_WORDS_OLD = 67`; a chain block runs `WORDS` 33 times: `YR_BLK = 62,249 + 33 x 27 = 63,140`) | `0` |
| `tbldk` bank of `n` rows (`E`, `G`, `F`) | `4,268 + 745 n` | `T_FAM`: `4,268 + 745 n` | `0` |
| `tbld` bank (`P1`: 2 variants; `fwd`: 4 variants, stride 128) | `26 + 40 n`; `48 + 80 n` | `26 + 40 n`; `52 + 400 n` | `0`; `+4 + 320 n` |
| `tbldf` bank (`BT`) | `4 + 52 n` | `4 + 52 n` | `0` |
| Instruction census of the variants (the parser's for the changed `leaf.s`, `secp_fe.s`; LLVM's for the rest) plus analytic code | `511,568 + 262,144 = 773,712 B` | `1,048,576 B` (`PROGRAM_BYTES`) | `+274,864 B` |

`FINAL` charges each `G`'s group-final block to one relay (`n < 1024` records): `fin1` (22) and its wait
(8), the relay's receive (23), `h2s`'s partial-path fixed part (32), and the records' excess of 6 each
over the full-block per-record shares they already pay,
`PDEF = ceil(1023 (6*1024 - (18 + 8 + 23 + 5773)) / 1024) = 322`: `FINAL = 407`.

### 3.4 Shared Table Banks (`T` endpoints): inventory, tags, traffic, build and memory

**What is shared.** Every lookup table of 8cc64cb7's leaf pipeline except window 0 of `F` and of `fwd` is determined by the pool (the six
curve points `X_1 .. X_6` that Phase 0 mines for a generator's pool `pool(c, m) < 2^14`) and a window digit. The row of digit `d` of a window at
shift `sh` is `alpha (d 2^sh) mod p`, reduced exactly as in 8cc64cb7 (an Eisenstein pair `(c1, c2)` for `E`, `F`, `G`; a torus fraction
`r + floor((r C0 + 2^255) / 2^256)` for `BT`; the residue itself for `P1` and `fwd`), with `alpha` from the eight coefficient families
of `callbacks.md` (`E`: `2 x^3`; `G`: `beta^2 x^2`; `F`: `-2 x`; `P1`: `2 y^2` and `2 y`; `fwd`: `-(A_0 + A_1)`, `-(A_2 + A_3)`,
`-(A_4 + A_5)`; `BT`: `beta`) and `gamma = 0`. The rows of window 0 of `F` and of `fwd` also depend on the cluster: its slice constant
`S_HIGH = (B >> 223) << 223` enters them (`-a^4 - 2 a S_HIGH` in `F`; `p_k - (7 + S_HIGH)(A_{2k} + A_{2k+1})` in `fwd`). Those two
windows (45 bits, `2^45` rows each) stay private tables of the generator, built as in 8cc64cb7. 8cc64cb7 keeps a private copy of every
table in each of its `2^58` generators although there are only `N_POOLS = 2^14` pools, so that `2^44` generators hold identical copies of every
table except those two windows. Here a pool's tables of the families and windows listed below are kept once, in banks (`fabric.md` Section 2.4).
The `BT` rows do not depend on the pool (`alpha = beta`), so its `2^14` copies are identical; they are kept per pool so that every read has the
same addressing rule (the cost is Row 9). A window's width is no longer limited by the `2^64`-byte local space, which 8cc64cb7 fills with
`11.05` EiB of five-window tables (the table regions of its `G` map); it is limited by the build cost (`2^w` rows times the builder's cost a row, for each pool).

**Inventory** (`certificate.py`: `T_FAM`, `T_BANKS_POOL`, `T_ROWS_POOL`, `T_BYTES_POOL`; recomputed independently in `numeric_check.py`):

| table | indexed by | shared windows (bits, shifts) | row bytes | builder, cost a row | banks a pool |
|---|---|---|---:|---|---:|
| `E` | `b = x(Y)` | `86, 85, 85` at `0, 86, 171` | 256 | `tbldk`, 745 | `2^33` |
| `G` | `3 x(Y)^2 mod p` | `86, 85, 85` at `0, 86, 171` | 256 | `tbldk`, 745 | `2^33` |
| `F` | `B_low = (x(Y)^3 mod p) mod 2^223` | `89, 89` at `45, 134` (window 0, 45 bits: private) | 256 | `tbldk`, 745 | `2^36` |
| `fwd` | `B_low` | `89, 89` at `45, 134` (window 0: private) | 128 | `tbld`, 4 variants: 80 (charged 400) | `2^35` |
| `P1`, one per `X_j` | `z^2` (256 bits) | `86, 85, 85` at `0, 86, 171` | 64 | `tbld`, 2 variants: 40 | `6 * 2^31` |
| `BT` | `w` (256 bits) | `86, 85, 85` at `0, 86, 171` | 64 | `tbldf`: 52 | `2^31` |

**Tags and banks.** A bank is named by `(pool, FAM, bank)`. `FAM` (6 bits) names the table and window: `E` 0 .. 2, `G` 3 .. 5, `P1` `6 + 3 j + k`
(point `j = 0 .. 5`, window `k = 0 .. 2`), `BT` 24 .. 26, `F` 27 .. 28 (windows 1 and 2), `fwd` 29 .. 30 (`T_FAMIDS`, `c_t_famids`). `bank` is the
digit above its `R` row bits, with `R = 54` for 256-byte rows, `55` for 128-byte rows and `56` for 64-byte rows: banks a window `2^32, 2^31, 2^31`
(`E`, `G`), `2^35` (each shared `F` window), `2^34` (each shared `fwd` window), `2^30, 2^29, 2^29` (`P1`, `BT`), all below `2^36`
(`c_t_bank_field`). The route latch of a read is `FAM << 36 | bank` and the fabric tag of the bank is role `B` with index
`2^58 + (pool << 42 | FAM << 36 | bank)`. There are `135,291,469,824` banks a pool and `N_T = 2,216,615,441,596,416` in all. Every bank is full
(rows times row bytes is a multiple of `2^62`: `c_t_banks_full`), so a bank backs exactly `[2^62, 2^63)`.

**Frame words.** The three-window kernels use these words of the 512-byte per-`X_j` entries (`EB`: 8cc64cb7's loop `card44_init`, `kernels/pointer_init.s`, writes their
pointer pairs in 47 instructions for six entries, and its Row 7 charges `ENTRY_X = 18` instructions an entry for all the words of an entry) and of the per-`Y` frames (written by the
producer `WORDS`, charged in `YR_BLK`). This package keeps the 18 in full, as if every word of 8cc64cb7's entry were still written, and charges the words below on top:
`ENTRY_X = 18 + 12` (Row 7). They take the place of the five window base addresses and of the two `BT1` base addresses of 8cc64cb7's entry (an entry keeps its size), and the loop
`t_entry` of `bank_checks.py` writes them in 54 instructions for six entries, 9 an entry (a `movz` for each of the three tags, three `stp` and one `add` an entry, three more for the `BT`
tags of each odd `j`, one `str` for `C0` of `j = 0`, two for the loop). Offsets are bytes from the frame base `e`. `experiment.py` builds both frames with helper functions and asserts
that the words the producer macros store (executed on 200 random variables, and no other word) and the `X` words of every `j` are exactly this list, and that `t_entry` writes the
`X` words of eighteen consecutive entries and no other word.

| frame | offsets | word | consumer |
|---|---|---|---|
| `X_j`, every `j` | `e+320, +328, +336` | `(6 + 3 j + k) << 36`, `k = 0, 1, 2` | `P1_3WIN`: tags of its three windows |
| `X_j`, every `j` | `e+344, +352, +360` | `2^62` | `P1_3WIN`: templates of the three row pointers (`bfi` patches the row bits in) |
| `X_0` | `e+376` | `C0` | as in 8cc64cb7 |
| `X_j`, odd `j` | `e+40, +48` | `25 << 36`, `26 << 36` | `BTBO` (`ldp x8, x9, [x3, #-152]`) and `BTBE` (`[x3, #360]`): `BT` tags of windows 1 and 2 (the tag of window 0 is an immediate of `BTBO`/`BTBE`); 8cc64cb7 keeps `BT1`'s two base addresses here |
| `Y` | `e+272`, `+304`, `+312` | tags of `E` windows 0, 1, 2 | `JK` (`KLKA3E`, `KLKB3E`) |
| `Y` | `e+280`, `+288`, `+296` | row pointers of `E` windows 0 (`+64`: its second line is read first and its link word returns the row), 1, 2 | `JK` |
| `Y` | `e+464`, `+496`, `+504` and `e+472`, `+480`, `+488` | the same for `G` | `JK` (`KLKA3`, `KLKB3`) |
| `Y` | `e+408` and `e+416`, `+424` | row pointer of `F` window 0 (private table, `+64`) and of windows 1, 2 | `JK` (`KLKA3F`, `KLKB3F`) |
| `Y` | `e+432`, `+440` | tags of `F` windows 1, 2 | `JK` |
| `Y` | `e+96` and `e+104`, `+112` | `fwd` row pointer of window 0 (private) and of windows 1, 2 | `LKP6_PAIR`, `LK1` |
| `Y` | `e+120`, `+368` | tags of `fwd` windows 1, 2 | `LKP6_PAIR`, `LK1` |

No word is used twice. The `X` words `e+40, +48` (`BT` tags) and `e+344 .. +360` (`P1` templates) are different words of one frame, and
`e+360` of an `X` frame is not `e+360` of a `Y` frame (the id word of the `Y`, 8cc64cb7's).

**Instruction sequences and traffic.** The consumers' extraction, tag and route-store instructions are part of the kernels `kernels/leaf.s` and
are counted in `FWD_INC` and `BWD_INC` (Section 3.3); the producers' `OFFS3_86` and `OFF3FF` are counted in `YR_BLK`. For the 144-leaf run (24
incidences, four `Y`) `experiment.py` reports, by table, the route stores and the aperture `LD64B`:

| table | route stores | `LD64B` |
|---|---:|---:|
| `BT` | 144 | 144 |
| `P1` | 72 | 72 |
| `E`, `G`, `F` (in `JK`) | 64 | 128 |
| `fwd` | 16 | 16 |
| all | 296 | 360 |

That is 15 `LD64B` and 12.33 route stores an incidence. A store serves the one `LD64B` of a 64-byte row and the two `LD64B` of a 256-byte row of
`E`, `G` or `F` (`KWINA`, `KWINB`: second line, then first); the two lines of a `fwd` row are read in two passes (`LKP6_PAIR`, `LK1`), each after its own
store. The returns are 64 bytes a `LD64B`, 960 bytes an incidence; the model prices the instructions the processors execute and no bytes in flight.
`keys`, the routing, the homes and the resolution read no bank. The three-window kernels execute `fwd` 1,322 and `bwd` 16,787 instructions for 24
incidences against 1,450 and 18,331 for 8cc64cb7's: 69.67 fewer an incidence (`(1,450 + 18,331 - 1,322 - 16,787) / 24`), after paying every tag, address and route-store
instruction. `certificate.py`'s remainder tables `FWD_REM = (8, 334, 659, 987)` and `BWD_REM = (0, 4193, 8383, 12572)` are the counts of runs with 1, 2 and 3 `Y`
points (343, 668, 996 and 4,217, 8,407, 12,596 with the fixed parts 9 and 24), and `experiment.py` asserts them equal to those runs; it also asserts that every table sum reads
three windows (`BT` 144 = 24 x 2 x 3 loads, `P1` 72 = 24 x 3, `E` and `G` 48 = 4 x 2 x 3 x 2 each, `F` 32 = 4 x 2 x 2 x 2, `fwd` 16 = 4 x 2 x 2).

**Register and pointer conventions.** The kernels keep 8cc64cb7's `ld64b` destination groups and register map; only the values that fill them changed. `P1_3WIN` takes its tags in
`x12, x13, x14` and its row pointers in `x15, x16, x17` from the one `ld64b x12, [x29]`; `BLK3F` keeps its prefixes in `x11, x8, x9` (`x8, x9` from `BTBO`/`BTBE`) and its pointer in `x7`;
`JK` loads eight words of the `Y` frame (two `neg_t` words, tag 0, pointers 0, 1, 2, tags 1, 2) into `x22 .. x29` with one `ld64b`, so its tags are in `x24, x28, x29` (`x24` is not used by
`F`, whose window 0 is the private table) and its pointers in `x25, x26, x27`; `LKP6_PAIR` and `LK1` keep their tags in `x23, x28` and row pointers in `x21, x22`. Every tag register is
written before its `str` and not rewritten before the `ld64b` that follows, and the `str` instructions do not change the flags. The experiment compares 144 exact labels and its negative
controls (below) show that a wrong tag or row is caught there; that is a finite check. The route latch and `POOLREG` belong to the bridge of the processing element that runs the leaf
kernels (`fabric.md` Section 2.4).

**Reads and sensitivity.** Over the full horizon the aperture serves at most `2^130.1165` `LD64B` (15 an incidence), on average `2^79.1392` an endpoint of the `N_T`
(an average, not a bound on one endpoint; the model prices no endpoint throughput; `tables` in the certificate summary, the expected number is smaller by the early stop). The claim prices every
aperture read at the instructions that form and use it, that is `k = 0` extra instructions. A premium of `k` instructions on each of the 15 loads an incidence would take `15 k` of the 69.67
instructions an incidence that the lever saves and raise the bound by about `0.022` bit a unit of `k` (it would stay below 8cc64cb7's for `k <= 4`), but it would no longer be the bound claimed here.

**Build.** Each bank is built by its own endpoint with the unchanged builders: `tbldk` (`E`, `G`, `F`: `4,268 + 745 n` for `n` rows),
`tbld` (`P1`: two variants, `26 + 40 n`; `fwd`: four variants at stride 128 with a zero variant for the second line's first word,
`48 + 80 n`, charged `52 + 400 n`) and `tbldf` (`BT`: `4 + 52 n`), all measured by `table_check.py`. The builder of bank `b` starts
at row 0 of the bank with `start_j = alpha_j ((b 2^R) 2^sh) mod p` and `step_j = alpha_j 2^sh mod p`: at most twelve `fe_mul_f` (105
instructions each, measured: six starts and six steps for `E`, `G`, `F`; four for `P1`, six for `fwd`, two for `BT`) and a reduction, charged `T_BANK_START = 2^15` a bank. The seed/link words of a row (`tbldk` writes
`addr + 120 = addr` and `addr + 248 = addr + 128`) hold the row's own address in the bank's local range `[2^62, 2^63)`, which is the
address a reader uses: they depend neither on the reader nor on the bank. `table_check.py` builds banks at production digit widths and compares every row word with the row function of the aperture model that the 144-leaf run reads (`callbacks.md`, "Table banks", step 5). The driver passes a 64-bit row
count `n = 2^R` and the destination `2^62`; every builder loop is a `subs`/`b.ne` on that count, so it ends after `n` rows with the pointer at `2^63`
(a 64-bit value, no wraparound; the test builds the last two rows of the first and the last bank, which end at `2^63`). Every word a kernel uses as
data is written by the builder: an `E`, `G` or `F` row has 30 entry words (six points of five words at byte offsets `0, 40, 80, 128, 168, 208`) and the
two link words at `+120` and `+248`, 32 words in all; a `P1` row holds its two 32-byte variants; a `fwd` row its four 32-byte variants, the third being
the zero constant; a `BT` row 32 bytes. The other 32 bytes of a `BT` row are never written by the builder (as in 8cc64cb7's `BT0` rows) and no kernel reads the registers they fill as data: the
first load of `BLK3F` (`x20 .. x27`) has its upper half overwritten by the sum that follows, and the unused half of the last load stays in `x16 .. x19`, dead until the next
definition. The experiment makes the model return a random non-zero word for every such byte (`POISON`) and the 144 labels still match.
Rows 8 and 9 charge it all: `2^114.4176` for the bank rows and starts, `2^113.1612` for the private windows of the `2^58` generators,
`2^75.4858` for the endpoints' program delivery and setup (`T_NODE = NODE_INIT + T_SETUP`: the pool's six points, the coefficient families, the
ready and barrier messages). All banks are built before the first leaf kernel: each ends its build with `DSB SY` (every store of the bank is complete) and sends one
ready message to the controller, and a global barrier (`BARRIER`, as 8cc64cb7's phase barriers) follows the last ready message; a bank is read-only afterwards.

**Memory.** Banks hold `2^112.9773` bytes, counted in `PEAK_BYTES` (`2^137.950614`, below `2^138`) together with the `N_T` endpoints' 16-MiB
code and control region, generic pages and program reserve (`T_BYTES`) and their fabric state (`FABRIC_BYTES`). A generator no longer backs
the five-window tables: its map keeps `F_W0_PRIV`, `FWD_W0_LINE0` and `FWD_W0_LINE1` and the aperture.

**What this does not show.** The experiment runs the kernels against a model of the aperture on 4 `Y` x 6 `X` points of one pool. For every aperture
read the model decodes the table, window and bank from the latch, asserts that the decoded digit fits its window (a bank or row field that is too wide would overflow
it) and returns the row at the address:
swapped `E` tags, swapped `P1` tags and a corrupted row break the match with the exact EC labels (negative controls), while random words in the unused half of every `BT` row do not (`POISON`). With one pool the effect of `POOLREG`
is not simulated (a `BT` read sent to another pool's bank would return the same row). It does not show the fabric's timing, bandwidth or
contention at `2^58` readers, which no cost of the model prices (as for 8cc64cb7's relays and homes); the aperture is premise 5 of Section 2.1.

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
| 4 | Phase-0 `X`-pair mining (`xemit` 1,024-record blocks, `MSG = 256`, `GEN_MUL = 912,688`) | `123.0661` | `0.0398%` |
| 5 | Phase 1-3 coin scalars | `100.2103` | `0.0000%` |
| 6 | Selected `X` setup (`MX = 6`) | `77.5850` | `0.0000%` |
| 7 | `EB` entry `X` & headers | `87.4919` | `0.0000%` |
| 8 | Table construction: private `F`/`fwd` window 0 per generator, every bank (`E`, `G`, `F`, `P1`, `fwd`) once per pool, the `T` endpoints' own work and the generators' `POOLREG` stores (`XTBLD_PRIV`, `T_BUILD_OTHER`, `T_NODE`, `POOLREG_INIT`) | `114.9172` | `0.0001%` |
| 9 | `BT` banks (`T_BUILD_BT`) | `106.7005` | `0.0000%` |
| 10 | `LC1`/`LC2`/`RC2` init + relay member tables (`RTBL_INIT`) | `100.0001` | `0.0000%` |
| 11 | Home `HH` zeroing (`HB = 50`, `DC ZVA`, 512 a pass) | `125.0029` | `0.1523%` |
| 12 | Home `R2`/`dlvh0` init | `86.3220` | `0.0000%` |
| 13 | Leaf generation + Y staging/transport + batch msg extra | `134.1828` | `88.3424%` |
| 14 | Fixup extra work | `107.5759` | `0.0000%` |
| 15 | Key canonicalization branch paths (`PJ_MEAN`) | `127.2359` | `0.7160%` |
| 16 | Routing + finals + scalar `card54_r` (`EMIT=7`, `K1=1024`: `H2S_BLK=5773`, `FINAL=407`; `BLOCKH=582` at `K2B=512`) | `130.8257` | `8.6216%` |
| 17 | Home joins (`NEW_INSERT=3`, `MAIN_VISIT=10`, `STOP_CHECK`) | `128.8058` | `2.1259%` |
| 18 | Same-`X` visits + `EXC_WORK` + `PAIR_INVERSE_EXC_WORK` | `90.1133` | `0.0000%` |
| 19 | `K2`/`d4` verification checks | `112.3483` | `0.0000%` |
| 20 | Candidate listings (`HCHG_LS1 = 6`) | `94.7611` | `0.0000%` |
| 21 | Defect rows (`DEFV + DEFK2 + DEFL + DEFJ + DEFWLT`) | `118.1046` | `0.0013%` |
| 22 | Candidate resolution (`RESOLVE`: six pool `X`s, six `2^21` reserves) | `116.5572` | `0.0004%` |
| 23 | Phase-4 generator requests (`REQ_SERVE`) | `113.0158` | `0.0000%` |
| 24 | Group exchanges | `105.4097` | `0.0000%` |
| 25 | Global barriers (`MSG = 256`) | `112.6732` | `0.0000%` |
| 26 | Phase control | `106.1257` | `0.0000%` |
| 27 | Final discrete-log output | `20.0000` | `0.0000%` |
| **Sum** | **`E_run`** | **`134.361591`** | **`100.0000%`** |
| **Gated** | **`E_total = 256 + p_gate * E_run` (`W = E_total / 1000`)** | **`133.003361`** | **`time_log2 <= 123.037576339062319`** |

## 5. Error Bounds (8cc64cb7's, Dominating the Three-Window Pipeline), Success Proof, and Verification Scope

This package sums three table rows where `8cc64cb7` sums five (`E` and `G`, `F`, `BT`, and the `P1` sums). Every bound below and in
`success.md` Section 7 is `8cc64cb7`'s five-window bound, kept verbatim, and it dominates the three-window pipeline: every row of
every window is a residue reduced into the same Eisenstein box of its class (Lemma E covers every `v` in `F_p`, whatever the
window), each sum has three terms whose classes and per-term error ranges are among those of 8cc64cb7's five-term sum (a sub-multiset of its classes, not of its digits: `E`: `A` x 3 of `A` x 5; `F`, `G`: `A` x 3 of `A` x 3 + `B` x 2), the `BT`
roundings add up to at most `3/2 + 2^-180` of `5/2 + 2^-180`, and the `c2` sums stay below `2^128`. `certificate.py` does not rest on that
prose: it recomputes, with the functions that give 8cc64cb7's tables, the error tables of the executed three-window sums
(`ERRH_USED = [[-161, 157], [-88, 91], [-241, 240]]`, `ERRL_USED = [[-195, 196], [-111, 106], [-296, 293]]`, from `KCL_USED`, three `BT` roundings and the `EB` entry) and checks every interval of them
contained in 8cc64cb7's (`ERRH = [[-214, 210], [-120, 117], [-320, 319]]`, `ERRL = [[-324, 325], [-188, 183], [-483, 480]]`; check `c_t_err_contained`), and pins the row counts and the `c2` range in
`c_t_dominated`. The `E1R`, `KCL` and `WLT` constants are therefore kept, and so are `ERRH`, `ERRL`, `NWIN`, `D_WJ`, `D_LAB_SYM`, `KB_C132`,
`D_WLT` and the cap tails.

**Transfer of `D_WLT` (the `P1` folds `c3` and `bc`).** Let `S` be the integer a `P1` accumulator sums before its final fold (for `c3` the three rows of one variant; for
`bc` the addend `1/D` and the three rows of the other), `x0 = floor(S / 2^256)` the carry count the kernel accumulates, `L = S - x0 2^256` the 256-bit
accumulator and `R = L + x0 C0 = S - x0 p`. The final `madd` adds `x0 C0` to the low limb of `L` and drops the carry out of it, so it drops one exactly
when `low64(R) < x0 C0`. `R` is congruent to `S` modulo `p` and lies in `[x0 C0, p + (x0 + 1) C0)`. With `s = S mod p`, the canonical residue that
8cc64cb7's lemma calls the clean output (a function of `Y` alone: the sum of the rows is `2 y_X b5` or `2 y_X^2 b5` mod `p`, whichever digits `b5` is cut into),
either `R = s` and the drop needs `low64(s) < x0 C0`, or `R = s + p` and it needs `s < (x0 + 1) C0`. Either way `s mod 2^64` lies in `[0, (x0 + 1) C0)`.
With each row below `p` and the addend below `2^256` (as in 8cc64cb7's Lemma V), `x0 <= 2` for `c3` and `x0 <= 3` for `bc` against `4` and `5` for
`P1_5WIN`'s five rows: the windows `[0, 3 C0)` and `[0, 4 C0)` of this package lie inside 8cc64cb7's `[0, 5 C0)` and `[0, 6 C0)` for the same function
`s(Y)` and the same fibre bound, so `D_WLT` transfers unchanged. `experiment.py` runs the shipped `P1_3WIN` up to its two `madd`s on exact integers (rows of its choosing, including rows that
put the low limb of the sum in the window): the counters equal `floor(S / 2^256)`, each fold returns `R` or `R - 2^64` exactly as derived, a drop occurs only with `s mod 2^64` in the window,
both outcomes occur, and `R = s + p` occurs only with `s < (x0 + 1) C0`.

**The other drop sites** (`F5`, `FDX`, `FS`, `FADDW`, the two `C1BSPLIT` sites, the `SQ13` columns) keep their windows: a site's window depends only on its clean residue and on the carry bound `T` that comes from the ranges
of its operands, and this package changes `T` only for `c3` and `bc` (the argument is in `success.md`, "Three windows in this package"). The lemma's other inputs (Lemmas F, K, R) bound the emitted label's deviation by the error tables
above, which contain this package's; the experiment reduces its rows with the classes of `KCL_USED` (`cert_class`), and its bank test compares the builders' rows with them.
The names `P1_5WIN` and `BLK5F` below and in `success.md` Section 7 are `8cc64cb7`'s macros; the bounds hold a fortiori for this
package's `P1_3WIN` and `BLK3F`:

- `ERRH = [[-214, 210], [-120, 117], [-320, 319]]` (`WH = 640`),
- `ERRL = [[-324, 325], [-188, 183], [-483, 480]]` (`WL = 964`),
- `NWIN = 5,100`, `D_WJ <= 2^-50.6837`,
  `D_LAB_SYM <= 2^-8.8460`, `KB_C132 / KB_BUCKET approximately 0.3335 < 0.34`,
- `P1_5WIN` (and a fortiori `P1_3WIN`) drops one carry term in `S` and `Q`: `D_WLT <= 2^-23.9393`,
- Cap and defect tails at `Q1GD = 27,136,000,000,000` (`NODES = 2^58`), with `RELAY_CAP`, `REQ_CAP`
  and `CAND_CAP` sized on the flat space:
  `d_f <= 2^-50.8391`, `d_p2 <= 2^-22.5351`, `d_cl <= 2^-28.1917`,
  `d_rq <= 2^-19.7123`, `d_rl <= 2^-21.8495`, `glob <= 2^-19.4303`,
  giving `lo_run >= 0.999844706048` and `p_gate * lo_run >= 0.39` in exact rationals.

Complete self-contained mathematical derivations appear in `success.md` (Section 7
contains self-contained proofs of 8cc64cb7's five-window Lemmas E, F, M, K, R, and V (which dominate this package's three-window sums, above), the
pointwise R3 joint-label and ordered-distinct affine-chain derivation, and the
global pair-inverse exclusion), `selected-x.md` (denominator-cleared selected-`X`
fiber bound `<= 30`), `fabric.md` (64-bit AXI fabric and scalar push/ACK services),
and `callbacks.md` (explicit scalar callback envelopes). Running `experiment.py`
verifies the complete exact-rational `certificate.py`, expands and checks instruction forms in all 7
AArch64 kernel files, executes `fwd`/`bwd`/`keys` on 4 real same-slice secp256k1
`Y` points x 6 real mined `X` points through the aperture model (`144/144` exact EC leaf matches), carries those 144 records
through modelled `H1` blocks and `fin1` into relay slots, where `bulk_r_receive` and `h2s` deliver them
into hop-2 windows that are checked against the exact labels, executes `EM`'s 1024-record block end
and its relay's full path, `h2s`'s full and partial paths with every partial-path chunk and all 2,048 `ADR`
trampolines (each paired statically with its own `tbnz` site and return point, and each run by a record that
fills its own hop-2 window), the stub's cap-clamp return, the handshake waits, the relay tables, Phase-0 `xemit`, `dlvh` 0-hit and 1-hit 512-record blocks,
`bulk_h_receive`, `fin`, `zero16384`, `ytransport.s:ychunk` (integrated 128-bit meter carry wrap and
threshold skip), `ystage.s:ycmd_gen` (accepting `k = Q1GD` and trapping `k > 2^45`),
`yboot_g`/`yboot_bank`/`yworker_sweep`/`ybank_sweep` `2^58`/`2^48` topology bounds, sweeps 5,000 random
256-bit three-window address and tag extraction vectors, executes the producer macros `OFFS3_86` and `OFF3FF` on 200 random
variables, builds bank rows with the unchanged `tbldk`, `tbld` and `tbldf` builders at production digit widths and compares them with
the row function of the aperture model, runs a finite model of `callbacks.md`'s six-`X` resolver (each of the 144
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
The complete backed intervals, all 2^58 B nodes and all `N_T` table banks enter PEAK_BYTES.
The emitted producer, receiver, replay and bank sweep use these same
constants; absolute-k tags and stored coins remain unchanged.

Ordinary fwdq1 extracts w3 bits 9..41 for the high cluster and keeps bits 0..8
for high k. For c'=32768,m=0,k=1,u=0, w2=0,w3=512 and destination=2^25.
Mining uses its separate fwdq metadata decoder. callbacks.md supplies the
full ordinary inverse and reconstruction from the original CN.

PROGRAM_BYTES reserves 1MiB per endpoint, including all measured assembly
variants plus 262144 bytes for the analytic shared-loop construction. The
cold census charges its delivery to G, R, H, B and the single controller; a table bank's delivery (`NODE_INIT`) and setup is charged in Row 8 (`T_NODE`).
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
is rederived here without relying on its verdict (this package's own windows are 85 to 86 bits, in shared banks: Section 3.4). The donor also cites
GoJackzi 4d4e224b-3364-42dc-b7d4-64004d8858b4,
hosein-ul 7f03684c-15b1-4b5c-ace4-28e465e828fc, and
kaankolcu 94ae0ccc-a4f9-4c53-918f-c82619c90433. The stored-leaf base retains
its credits to kaankolcu e831f89f-27de-4fa9-becd-289bd7a4f4a8,
the grouped schedule in eed62949-dc73-4549-8271-372b70f6ce84, and
aburan28's relay reset fold 5f726b85-b30c-49e1-b4fa-33b8a1ab0f55.
Public status and earlier scores do not serve as mathematical premises.

## 7. Changes

### 7.1 Changes from 8cc64cb7 (this package)

Built on `8cc64cb7` (`bradcush`) and, through it, `8fff60c5` (`GoForceX`), `838b03f3` and `283cdacf` (`soispoke`), `4ce2d07e`
(`bradcush`). One lever, and what it touches (every kernel count below is executed by `experiment.py`):
1. **Shared table banks (Section 3.4, `fabric.md` Section 2.4, `callbacks.md` "Table banks").** The lookup tables of a pool are kept
   once in `N_T = 2,216,615,441,596,416` passive banks (`135,291,469,824` per pool) instead of in each of the `2^58` generators, and read through the
   T aperture. A bank is built by its own endpoint with the unchanged builders. Bank bytes `2^112.9773` count toward the
   peak memory `2^137.950614`; the builds are charged in Rows 8 and 9 (`2^114.4176` for the banks, `2^113.1612` for the private
   windows, `2^75.4858` for the endpoints' delivery and setup; 8cc64cb7's private tables cost `2^122.3515` and `2^117.1997`).
2. **Three windows instead of five** in `BT`, `P1`, `E`, `G` (`(86, 85, 85)`) and in the shared part of `F` and `fwd` (`(89, 89)` after a
   private 45-bit window 0): `BLK3F`, `P1_3WIN`, the three-window `JK`, `LKP6_PAIR`, `LK1` (`kernels/leaf.s`), `fwd` 1,450 -> 1,322 and `bwd`
   18,331 -> 16,787 for 24 incidences; per incidence 15 `LD64B` and 12.33 route stores through the aperture (`experiment.py`).
3. **Y producer:** `OFFS3_86` and `OFF3FF` (`kernels/secp_fe.s`) store pointer and tag words, +27 instructions a `Y`; a chain block produces 33 `Y`, so `YR_BLK = 62,249 + 33 x 27 = 63,140`.
4. **X entries (Row 7):** `ENTRY_X = 18 + 12`: 8cc64cb7's 18 an entry in full and 12 an entry for the new frame words (Section 3.4; their writer `t_entry` runs 9 an entry).
5. **Bounds:** unchanged (Section 5): 8cc64cb7's five-window bounds dominate the three-window pipeline.
6. **Package:** a third archive `runtime3.json` holds the LLVM-audited versions of `leaf.s` and `secp_fe.s`; `audit_delta.py` compares them
   with the shipped ones (Section 8); `table_check.py` also builds banks; `numeric_check.py` recomputes the bank census independently.
   The rest of the package differs from 8cc64cb7's only as follows: `success.md` has one paragraph added at the head of Section 7 and the
   figures of its parameters paragraph regenerated; `certificate.py` also evaluates `PHI` by integers (`phig_num`: the grid value of 8cc64cb7's Fraction
   evaluation by an exact split of its two terms and an exact comparison where the 192-bit floors do not decide the carry; check `c_phi_int` compares the two on 23 groups); `experiment.py`, `numeric_check.py`, `table_check.py` and `run.py` are changed as this section and Section 8 describe.
   `selected-x.md`, the other five kernel sources, the Phase-0 mining, the routing, the home joins, the candidate resolution and the caps are 8cc64cb7's.

### 7.2 8cc64cb7's changes from 8fff60c5 (unchanged here)

Built on `8fff60c5` (`GoForceX`), `838b03f3` and `283cdacf` (`soispoke`), and `4ce2d07e` (`bradcush`), the
stored-leaf design those build on. Every kernel count 8cc64cb7 changes or adds is executed by `experiment.py` (Section 3.3 lists them,
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
2. (8cc64cb7's `P1_5WIN` window 4 by one `BIC` with a reversed table is not used here: `P1_3WIN` has three windows.)
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
   and the tests in Section 5 (8cc64cb7's; this package adds `runtime3.json`, Section 8).
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

The eleven members named in runtime.json (the 7 scalar assembly source files and
kernels/a64.py), runtime2.json (experiment.py) and runtime3.json (the audited 8cc64cb7 versions of leaf.s and secp_fe.s, as
kernels/base_leaf.s and kernels/base_secp_fe.s) are stored as bounded zlib/base85 UTF8 data
with exact uncompressed lengths and SHA256 hashes, three archives so that every file stays within
65,536 bytes. run.py restores them in the credential-free experiment directory, rejects wrong
hashes, lengths, paths and trailing data in any archive, then runs experiment.py. All proof and
numerical files are plain text. The optional source input prints one restored source without
running the suite. Compression changes storage only and earns no work credit. Every mathematical
or algorithmic input is inside this package.

Section 9 of the experiment lives in encoding_checks.py, which experiment.py runs in its own namespace
(an archive member is capped at 65,536 bytes). The default experiment also runs numeric_check.py's separate root inequalities,
320-bit integer logarithm interval, gate minimality, 27-row sum, bank/cold
census and pair-map finite fibers. archive_check.py exercises actual emitted
archive/forwarding boundaries and old-layout controls; math_checks.py verifies
the GF(1024) and selected-X finite coefficient facts; table_check.py runs the table builders,
`eisr3` on every `ksel` cell and, bank by bank, the unchanged builders against the experiment's row function at production digit widths. No finite experiment
proves the production degree or collision argument. The source instruction
checker and functional emulator remain separate from those analytic premises.

The independent syntax/size receipt encoding-audit.json records, for 8cc64cb7's sources, 69 LLVM
(`llvm-mc` 23.1.5) assemblies of expanded instruction lists, including 64 generic role/PE
bindings. It records the exact audited kernel source hashes and the 518756-byte
per-endpoint source census of the audited sources. Five of the seven sources are shipped as audited. `leaf.s` and `secp_fe.s` are shipped
changed (three-window lookups, producer macros; 511,568 bytes a role with the changes): **their changed instruction texts were not assembled by
LLVM, because no assembler was available when this package was built.** `encoding-audit.json` (`shared_t_delta`) lists every instruction text whose count differs
from the audited expansion (leaf.s: 41 new, 84 removed; secp_fe.s: 65 new, 57 removed, branch labels normalised), and the experiment
(`audit_delta.py`) checks the audited versions' hashes, recomputes the lists, and checks that the form (mnemonic and operand kinds) of each
new text occurs among the audited instructions; the encoder checks of `kernels/a64.py` (register classes, immediate ranges, bitmask
immediates) run on every instruction when the sources are loaded, and the kernels run against exact EC labels. Assembling the listed
texts with `llvm-mc -triple=aarch64 -mattr=+lse,+rand,+ls64,+rcpc-immo,+cmpbr` is the missing check. The default Python experiment checks those
hashes and counts, and simulates finite paths; it does not invoke LLVM.
Raw macro argument quoting and the limits of parser/linked-program evidence
are disclosed in that receipt.

# Gated Phase-0 birthday X-quad Y-cluster search on secp256k1 (123.492703530915914)

**Claim.** On uniform secp256k1 `Q=[x]P` (`r=|G|`, `p=2^256-2^32-977`), this cold zero-advice single-run attack finds `x` with prob `>= 0.3901 > 0.39` and expected work `time_log2 = 123.492703530915914 = 61746351765457957/500000000000000` u (`C=1000` scalar AArch64 insns/u; `R_worst <= 2^179.732829` insns if gate opens and all groups, walks `<= CAP_H`, and caps saturate). Gen `0` draws 128-bit `RNDR` `u` (`256` insns unconditional) and runs `route.s` `gate` (`u < K_GATE = 132920808407662439192763145161662651734`, `p_gate = K_GATE/2^128 = 0.390619148...`), halting if `u >= K_GATE` else running Phase 0 + `q1=8,700,000,000` batches in `G=1536` groups (`m_g=floor(q1*(g/1536)^(28/32))`, exact `lo_run >= 0.998670985595`, `p_gate*lo_run >= 0.39010001 > 0.3901`, `(K_GATE-1)/2^128*lo_run < 0.39010001`; `E_run <= 2^134.814654` insns = `2^124.848869` u, `<= 74.0972`/leaf on `S <= 2^128.603307`; `E_total = 256 + p_gate*E_run <= 2^133.458488` insns; mem `<= 2^136.9996` B; preproc `<= 2^117.3477` u incl. `E_mine <= 2^124.9211` insns). `certificate.py` checks outward `2^c` intervals, §5 VA/MMIO disjointness, `STATED==DECLARED`, `K_GATE`, `CAP_H`, the 27 upward §7 ledger cells, digit-bounded counts of the 16 `_FIGS` across `proof.md`, `success.md`, `submission.json`, and §5/rate strings.

**Reviewer's guide.** §§1–2: algorithm; §3: kernels; §§4–5: VA; §§6–7: ledger (`success.md`); §§8–9: checks.

## 1–2. Algorithm, Phase 0, and assumptions

For `u < K_GATE` (`p_gate = K_GATE/2^128 = 0.390619148...`), `2^68` gens, `2^42` relays, `2^84` homes use 52-bit VAs (`FEAT_LVA` `[0, 2^52)` & `[2^64-2^52, 2^64)`, `FEAT_D128` `< 2^56` B PA/node, `FEAT_LSE`, 256-B `DC ZVA` `BS=6`, `FEAT_LS64`, `STLR`/`STLUR`+`DMB ISH`, two-wave barriers `<= 30,000`). **Assumptions:** ideal first-try `RNDR` bits (two `>= r` draws in `COIN_SCALAR=64` use `1`, `< 2^-254`/draw in `glob`); Phase-0 pack/emit, `2^32`-head join, report/drop, coin fetch, recompute, and broadcast use analytic bounds (`COST_MODEL.md`).
- **Phase 0 (`E_PRE <= 2^117.3477` u, `E_mine <= 2^124.9211` insns):** `g = 1024c+m` (`0<=m<1024`) joins `pool(g) = 32m + (((c&31)^gf1024_mul((c>>5)&0x3ff,m))&31) < 2^15` (`2^53`/pool), sets `R_g=[b_g]Q`, runs `2^32` chains of `32,768` pts (`M_POOL=2^100`/pool): stores `(c_(g,b),s_(g,b))` in `[59 TiB, 59 TiB + 256 GiB)`, sets `P_0=R_g+[c_(g,b)]P`, builds `Q_1..Q_16384` in `[60 TiB, 60 TiB + 4 MiB)` (`secp_fe.s` `xmine_qdbl` `23,105` + 8-round `xmine_qtab_step` + `64` setup = `10,163,673`/chain; `_qtab_safe` checks `y!=0`, distinct `x`), runs `xmine.s` (`<= 963+132+32`) routing 32-B records (`A=x^3+7 mod p`, `d_1=(A>>63)&(2^42-1)`) to Homes (`low 2^32` heads zeroed, `R2` cap `CAP_H`, `53` insns/pt, `R2` reset). Collector `g_p` caps reports at `64` pairs (drops `>=64` in `<=8` insns), fetches `<=128` coins, recomputes `<=128` chains (`32,768*1,000`/chain in `ROUTE_MINE`), aborts if `<8` disjoint distinct-`A` pairs (`0 < A_(p,2q+1)-A_(p,2q) = delta_A_(p,q) < 2^63`), else broadcasts `MX=16` pts+scalars; gens build `40×44-TiB` `(37,37,37,37,36,36,36)` X-tables (`tbld` `x6`) in `[64 TiB, 1,824 TiB)` + `BT`; Homes zero `2^47` heads (`zero16384`).
- **Phases 1–5:** In group `g` (`q_g` batches), each cluster's `1,024` gens build `cdiv(q_g, 1024)` Y-chains of `L_chain=63,549*33=2,097,117` affine pts (`yround`, `ysnd=24,ydlv=32,ydlve=23`). Per batch (`Q2=16*L_chain`), gens run `fwd+fe_inv+bwd+keys+fin1` (`w2,w3`=`68+26+34` b), Relays `h2s` (`95` full, `<=15+6n` part)+`fin` (`K2B=16`), Homes `dlvh` (`DLVHB=31`, `HB=47`, `CAND_CAP=2^27`), Phase 4 (`REQ_CAP=49*2^33`, `RELAY_CAP=2^40`) tests `u in {+-1,+-lambda,+-lambda^2}` and for `nu'!=u nu (mod r)` checks `[x]P==Q` at `x=(u mu-mu')(nu'-u nu)^(-1) mod r`, latching stop.

## 3. Leaf kernels and labels (`kernels/leaf.s`, `sites.py`)

With `D_j = b^3-a_j^3`: `fwd` (`1081/16 = 67.5625`/inc, `7` fixed) multiplies `e_q = D_(2q)D_(2q+1)` (`LKP7_PAIR` `97`+`P1X` `84`); `fe_inv` (`22,549`); `bwd` (`16233/16 = 1,014.5625`/inc, `7` fixed) runs `EFG_QUAD_4` (`2 delta_A_q < 2^64`), `PAIR7C` (`157`), `C1BSPLIT_BDAY` (`122`, `delta_A_q < 2^63`), `P1_7WIN` (`254`), `PAIRB6` (`93`), `3×TFP4_7`. On `{x, beta x, beta^2 x}` (`v*` min), `keys` (`513/4 = 128.25`/inc, `5` fixed) emits `2 × 106`-bit `(comp(v*), comp(beta v*))` -> `d1=42, d2=42, K2=64, D4=22` (`170` bits; box `< 1+1/3195`, `D_sym=[-389,264]^3` `mu_6`-invariant, `116` fixup sites/quad, cap `2^52`).

## 4–5. Completion and 52-bit-VA map (`certificate.py`)

`staddl`, Relay's 128-bit `D`, and two-wave barriers enforce completion; `disjoint_va` checks non-overlap:
- **Gen (`<= 2^52.9589` B):** `LC1` `[0,32 TiB)`, `257` `YR_s` `[32 TiB+2 GiB-64,33 TiB)`, `EB` `[33 TiB,+2^34)`, MMIO/`REQ_CAP` `[34 TiB,59 TiB)`, coins `[59 TiB, 59 TiB + 256 GiB)`, `xmine` `[60 TiB, 60 TiB + 4 MiB)`, X-tables `[64 TiB, 1,824 TiB)`, `d1` `[2^51,2^52)` & `[2^64-2^52,2^64-2^51)`, `BT` `[2^64-2^51,2^64)`.
- **Relay (`<= 2^52.0661` B):** `code/stack`, `LC2` `[32 TiB,64 TiB)`, `RC2` `[96 TiB,128 TiB)`, `E2/p` `[128 TiB,192 TiB)`, `DONE/PW/MB2`, `RELAY_CAP` `[256 TiB,320 TiB)`, `LINK` `[4,095 TiB,+8 KiB)`, `REL_UP` `[2^64-2^52,2^64)`.
- **Home (`<= 2^52.9996` B):** `2^47` heads `[2^64-4,096 TiB,+1,024 TiB)`, `next` `[2^64-3,072 TiB,+2,388.5 TiB)`, `R2` `[-682 TiB,+4,095 TiB)` (`CAP_H = 164,136,470,183,936`), `LINK`/`CAND_CAP` `[4,095 TiB,4,096 TiB)` (`<= 2^52.9996` B/node; peak `<= 2^136.9996` B).

## 6–7. Ledger (`success.md` §§6.1–6.7, `certificate.py`)

| opened-run | `log2` insns (`<=`, upward 4-dec) |
| --- | ---: |
| `E_MINE` | 124.9211 |
| init;`R_g`;bc | 104.0001; 87.2999; 100.9079 |
| coin;`X`;hdr | 98.6034; 89.0001; 97.0000 |
| X-tbl;`BT` | 118.1033; 117.9658 |
| LC;`z16k`;`i3` | 112.0001; 127.0029; 88.3220 |
| leaf+Y | 134.6584 |
| fix;keys | 103.9914; 127.0358 |
| route | 130.9766 |
| join+walk | 129.3515 |
| `SXV`;`d4`;`EL`;df | 38.7748; 110.4725; 92.3770; 54.1862 |
| res;req | 110.2999; 110.0942 |
| ph4;bar;ctl;o | 104.8255; 109.8727; 104.0001; 20.0000 |
| **`E_RUN`** (`<=74.0972`/leaf, `2^124.848869` u) | **134.814654** |
| **`E_TOTAL = 256+(K_GATE/2^128)*E_RUN`** (`0.390619148...`) | **133.458488** |

`W_E = E_TOTAL/1000` -> **`time_log2 = 123.492703530915914`** (`E_PRE/1000 <= 2^117.3477`; `E_mine <= 2^124.9211`; `R_worst <= 2^179.732829`).

## 8–9. Checks and provenance

- **Checks:** `experiment.py` (`{"seed":1,"seconds":12}`, `~12s`: `xmine_qdbl`/`xmine_qtab_step`, `xmine.s`, `yround`, `leaf.s`, `tbld_production`, `sites.py`, `C_vandermonde_cs`, `machine.py`); `encoding_fuzz.py 1 4000` assembles kernels with LLVM (`llvm-objdump == a64_load`: `25,042;3,777;2,373;1,122`).
- **Credits:** `17a007c8`..`8f53489b`, GWZ, AGG, `kaankolcu` `38df975c`..`95cf6938`, `aburan28`, `Phoenix0075`, `pepedesigner`.

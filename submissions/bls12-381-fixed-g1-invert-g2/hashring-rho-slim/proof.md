# Frobenius-class-folded distinguished-point rho for fixed-G1 pairing inversion

Track `bls12-381-fixed-g1-invert-g2`. Claim: `time_log2 = 125.818` (125818/1000), success probability at least 0.39. This package is the accepted `hashring-rho-sel1exact` (125.82) with the rounding slack removed: the per-slot charge drops from 1,749 to the executed maximum-over-slot-types of 1,740 (the 9-instruction explicit margin covered no row and is deleted), and the claim is stated at thousandths instead of hundredths; `proof.md` is condensed to fit the 65,536-byte intake cap with no figure, table, lemma, or check removed. The walk, selector, window, budget, store, replay, `mu = 11/10` and every additive item are unchanged; only the unit and the claim precision move.
All objects below are the ones of `tracks/bls12-381-fixed-g1-invert-g2/target-profile.json`; the unit of work is one dense Fp12 multiplication (`C = 20,000` scalar instructions). Everything the attack does is charged once.

## 0. Summary

1. `P` fixed, `g = e(P,Q0)` generates GT (prime order `r`). Write `z = g^a`, return `Q = [a]Q0` (section 1).
2. Frobenius `x -> x^p` acts on GT as `x -> x^lambda`, `lambda = p mod r` of order 12: free action on `GT \ {1}`, 12-element classes, `N = (r-1)/12` (section 2); `sigma^6` is inversion.
3. Class selector (3 Fp mults, three modular subtractions, one negation) maps one Fp2 coordinate `c1` to the 12 first-coordinates of `c1(sigma^k(x))`, invariant key `kappa` (their minimum) and covariant twist `k*` (section 3). `d = g0*(a-b)` is one mult (`gamma = g0*(1-u)`); min over six stored values and their negations via one negation of the max.
4. Lifted adding walk `y' = y * sigma^{-k*(y)}(T_h)`, `h` a hash of `kappa`: class map at one dense Fp12 mult/step (Lemma 3.2, section 4).
5. Look-back window of 24 class keys with deterministic probing removes fruitless short cycles; meeting trails stay together up to `< 2.7e-4` (sections 4.3, 6). Window = 24-slot ring of keys + parallel ring of stored 64-bit hash words `w_1` (no second hash); bounded 3-instruction slot update; accept = ring write; membership = 24 unrolled load/compare/branch + exact six-limb test on hash hit (section 5). Collisions only reject non-members (charged probe). Residual cycles `> 24`: `<= 2^-163.5`/step (section 7).
6. Independent trails end at distinguished points (density `2^-80`), report to ordered store; colliding reports replayed with exponent tracking, log recovered via `sigma^m` alignment (section 4.4). No exponents in main loop.
7. `S_max = ceil(sqrt(2 mu N U)) = 68967528808679891282472564154805017972` (`2^125.697`), `U` = exact rational upper bound on `ln(10000/6097)`; slot charged `1 + 1740/20000 = 21740/20000`; additives `2^92.361`. `W_max <= 2^125.818` (`log2 = 125.81778482...`), claimed `2^125.818` (section 5). Memory `2^53.6` B, preprocessing `2^24.13`, no advice.

## 1. Reduction

Nonidentity `P`: `e(P,.)` nondegenerate (`r` prime), `g = e(P,Q0)` order `r`. `z = g^a`, `a` uniform on `Z/r` indep. of `g`. `e(P,[a]Q0) = z`; `z = 1` iff `a = 0`, return `O`. Pairing convention only fixes the map. Charge: one pairing (`2^30` instr), later `[a]Q0` and `g^a == z`.

## 2. Frobenius classes of GT

`r | p^4 - p^2 + 1`; `lambda = p mod r`: `lambda^12 = 1`, `lambda^6 = -1` (`algebra`). `x^p = x^lambda` on GT; `sigma(x) = x^p` order 12. `x != 1`, `0 < k < 12`: `sigma^k(x) = x` needs `r | lambda^k - 1`, impossible; classes have 12 elements, `N = (r-1)/12`, `log2 N = 251.2721`. `sigma^6` = inversion. Bookkeeping: `sigma^m(g^alpha z^beta) = g^{lambda^m alpha} z^{lambda^m beta}`, `sigma^{-k}(T) = T^{lambda^{-k}}`.

## 3. Class selector in the BLS12-381 tower

Tower `Fp2 = Fp[u]/(u^2+1)`, `Fp12 = Fp2[w]/(w^6-xi)`, `xi = 1+u` (pinned `Fp[w]/(w^12-2w^6+2)`, `u = w^6-1`; blst tower, `w^2 = v`). `x = sum_{i<6} c_i w^i`; `w^p = gamma w`, `gamma = xi^{(p-1)/6)`; `sigma(sum c_i w^i) = sum conj(c_i) gamma^i w^i`. `c1` evolves by `tau(c) = gamma conj(c)`; `tau^2(c) = zeta c`, `zeta = gamma conj(gamma) = xi^{(p^2-1)/6)` in Fp, primitive 6th root (`zeta^3 = -1`, `zeta^2 = zeta-1`). `c = a+bu`, `gamma = g0+g1 u`:

    V_{2j} = zeta^j a, V_{2j+1} = zeta^j d, d = g0 a + g1 b, j = 0,1,2; V_{k+6} = -V_k,

since `Re(zeta^j gamma conj(c)) = zeta^j(g0 a + g1 b)`. Twelve values from 3 Fp mults (`zeta^2 = zeta-1` for the rest, sign flips).

**`d` in one multiplication.** `gamma = g0(1-u)` exactly (`g1 = -g0`): `gamma^6 = xi^{p-1} = -u`, `(1-u)^6 = 8u`, so `(gamma/(1-u))^6 = -1/8` with base-field root `g0`, `g0^6 = -1/8`, `zeta = 2 g0^2`. Exact congruences (`g1+g0 = 0`, `8g0^6+1 = 0`, `2g0^2-zeta = 0`) asserted by `algebra`. `d = g0(a-b)`: one mult.

**Definition.** `kappa(x) = min_k V_k(x)` (stored-representative order; any fixed total order works), `k*(x)` = smallest attaining `k`.

**Lemma 3.1 (covariance).** `V_k(sigma^m x) = V_{k+m}(x)`; `kappa` invariant, `k*` shifts by `-m` when unique. *Proof.* `c1(sigma^m x) = tau^m c1(x)`. ∎

**Lemma 3.2 (twisted step is a class map).** `x' = x sigma^{-k*(x)}(T)`: `sigma^{k*(x)}(x') = sigma^{k*(x)}(x) T`, and `sigma^{k*(x)}(x)` is the class's unique member with `k* = 0`. Class of `x'` depends on `(class of x, T)` only. ∎

Exceptions (ties, `c1 = 0`: 8 linear conditions, `<= 8/p < 2^-378`/state): smallest index used, Lemma 3.2 may fail; counted in L4.

**Min over twelve from six.** `{V_0..V_11} = {V_0..V_5} U {(p-V_k) mod p}`; `V_0 = a`, `V_1 = d`, `V_2 = zeta a`, `V_3 = zeta d`, `V_4 = (zeta-1)a`, `V_5 = (zeta-1)d`. (1) `m = min V_k`, `M = max V_k` with smallest indices (10 six-limb compares + cmove). (2) `nb = p-M` in `[0,p)` (one negation; `M = 0` gives `0`, minimum only if all zero). (3) `m <= nb`: `(kappa,k*) = (m,i_m)` else `(nb,i_M+6)` (11th compare). Same pair as direct 12-scan (stored half wins cross-half ties by index; strict `<`/`>` inside halves). Checked on all `7^6 = 117,649` tuples (`ledger`) and on real GT elements + conjugates (`algebra`). Six negations not performed: 4 modular ops (`a-b`, `zeta a-a`, `zeta d-d`, `p-M`).

**One-pass extrema.** Values in memory; seed: 6 `LDR` of `V_0`, copy to both extrema (12 moves), zero indices: 20. Each of 5 more: `MOVZ` index, 6 `LDR`, six-limb chain vs min + 7 `CSEL` (13), same vs max (13): 33 each, `20+5*33 = 185`. Chains write `XZR` (carry only). Row `185+14 = 199` for 11 compares; same tie rule as two scans. No early exit claimed (top limbs can tie).

**Eleventh compare, exactly.** `nb - m` chain: carry 1 iff `m <= nb`, equality included; no `Z` read (last limb's flag only). 6 `CSEL` (`CS`) on value, 1 on index, 1 `ADD` (`i_M+6`): 14. Executed with negation (218/tuple, 212 when `M = 0`) on 245 small-field (51 `m = nb`, one all-zero) + 40 production tuples. Negation `0-M` in x10-x15: chain, `B.CS` (taken iff `M = 0`, leaving `0`), else `ADDS` + 5 `ADCS` of `p`; no mask.

**Selection cost.** 3 `blst_fp_mul` x 383 = 1,149; 3 subtractions x 31 + negation 13 + 18 (`p` loads, 3 blocks) = 124; extrema + 11th = 199. Total 1,472. `zeta, g0` stored Fp; no Montgomery conversion.

## 4. Algorithm

### 4.1 Parameters and primitives

`R = 2^14` indices; window `w = 24`; probes `A = 8`; DP bits `b = 80`; `L_cap = 2^87` slots/trail; `mu = 11/10`; caps `D_max = ST_max = 2^47`, `REC_max = 16`. Coins: 256-bit `K`; ChaCha20 + domain separation; scalars = 384-bit output mod `r` (distance `<= r/2^384 < 2^-129`). `khash(kappa,salt) = mix-chain` (splitmix64 finalizer) over 6 limbs; `w_1, w_2` from two coin words; index `w_1 mod R`; DP iff `w_1 = 0` and `w_2 mod 2^16 = 0` (80 bits, class functions). Window stores `w_1` itself: no third word computed.

### 4.2 Setup

1. `g = e(P,Q0)`; `z = 1` -> output `O`. 2. Per `h < R`: ChaCha20 `(u_h,v_h)`; `T_h = g^{u_h} z^{v_h}` (square-and-multiply, squarings as mults, `<= 2*509+1` mults); `TW[h][k] = sigma^{-k}(T_h)` (<= 11 Frobenius, 1 unit each). Table `12R` elements (113 MB).

### 4.3 One trail

Trail `i` (counter): (1) Start: ChaCha20 `(alpha_0,beta_0)`; `y_0 = g^{alpha_0} z^{beta_0}` (<= 1,019 mults); drop if `y_0 = 1` or `c1 = 0`. `(kappa_0,k_0) = Sel(y_0)`; seed ring: `kappa_0` in all 24 key slots, its word in all 24 hash slots, `p = 0` (window = `{kappa_0}`). (2) If `kappa_n` distinguished: report `(kappa_n,i)`, end. If slot counter `c` (every mult in step 3) hits `L_cap`: end silently. (3) Else `h = hidx(kappa_n)`; `a = 0..7`: `y' = y_n TW[(h+a) mod R][k_n]` (1 dense mult); `Sel(y')`; reject if `kappa'` in window `W`, else accept (`y_{n+1} = y'`, append, `n += 1`). Append = ring write + bounded slot update; reject on 64-bit word match (exact on repeats; collision over-rejects, charged escape, section 5). All 8 rejected: end silently. Identity needs no test: table element uniform (H3), `Pr[y' = 1] = 1/r`/slot, `<= 2^-128.857`/run; selector returns `0`, next step ordinary (L4).

**Class dynamics.** By Lemma 3.2, class of `y_{n+1}` = class of `c_n T_{h+a}` (`c_n` canonical rep). First-probe accept: deterministic class function `F(C) = class(c T_{hidx(kappa(C))})`. Later probes only on closed patterns `<= 24` (section 7). Decision depends on current + 24 previous class keys only.

**Exponent invariant (replay only).** `y_n = g^{alpha_n} z^{beta_n}`; `(alpha,beta)_{n+1} = (alpha,beta)_n + lambda^{-k_n}(u_{h'},v_{h'})`.

**Coalescence.** Two trails at class `C`: next accepted class `F(C)` for both unless in one window (`~1/(12R)`/step each, L2); windows identical after 24 steps. Hash adds `<= 24*2^-64`/candidate, reject-only: coincide except `<= 2*24*nu`, `nu <= 1.1/(12R)`; both reach same DP.

### 4.4 Store, collision, recovery

Ordered map (red-black) full 48-byte `kappa` -> trail index. Report `(kappa,i)`: insert if absent; `Recover(i',i)` if present with `i' != i`; else ignore (one report/trail). Full-key compares, no fingerprints.

**Recover(i',i).** Replay both trails (same coins + `(alpha,beta)` tracking). Endpoints `(y_A,...)`, `(y_B,...)`; discard unless keys equal. Find `m < 12` with `sigma^m(y_A) = y_B` (<= 11 Frobenius, 12-coefficient compares); else discard. `sigma^m(y_A) = g^{lambda^m alpha_A} z^{lambda^m beta_A}`:

    a(lambda^m beta_A - beta_B) = alpha_B - lambda^m alpha_A (mod r).

Denominator nonzero (`1-1/r`: trajectories depend on coins only via `alpha_0+a beta_0`, `u_h+a v_h`, uniform on `(Z/r)^2`): `a = ...`, check `g^a = z` (<= 509 mults), output `Q = [a]Q0`. Replay keeps store to `(kappa, index)`.

### 4.5 Meter and termination

Counters `S` (accepted), `M` (slots), starts, reports, recoveries. New trails while `S < S_max`; current trail finishes (<= `L_cap` slots, charged). Fail if `M > M_max + L_cap`, starts `> ST_max`, reports `> D_max`, recoveries `> REC_max`. All loops bounded (fixed ChaCha20 blocks, `L_cap` slots).

## 5. Work ledger

Unit: dense Fp12 mult = 1; else instructions/20,000. (Exact: `ledger`.) One `blst_fp_mul` = 383, one dense Fp12 = 22,996 = 1 unit (`COST_MODEL.md`). Every row a hard upper bound on executed instructions.

| item | instr |
| --- | ---: |
| dense Fp12 mult by table entry | 1 unit |
| selector: 3 `blst_fp_mul` x 383 | 1,149 |
| selector: 3 mod-subs x 31 + negation 13 + 18 (one six-limb `p` load per block of 3; allocation below) | 124 |
| selector: extrema `20+5*33` + 11th (6+1 `CSEL` + `ADD`) | 199 |
| two `khash` words: 12 x (`EOR` + 8-instr mix), limbs in registers, + 2 accumulator + 2 multiplicand loads | 112 |
| DP test: `AND x3,x2,#0xFFFF`, `ORR x3,x3,x1`, `CBNZ x3` (64-bit: 80 bits) | 3 |
| membership: 24 x (`LDR` word, `CMP` vs `x1`, `B.EQ`) = 72 + one 20-instr hit path (address, 6 limbs, `6 EOR+5 ORR+CBZ`, branch) | 92 |
| ring write: 3 shifted adds + 3 `STP` + `ADD x11,sp,p,LSL #3; STR x1,[x11]` | 8 |
| bookkeeping: slot update + accepted counter w/ cap + 2 branches | 12 |
| plumbing, table address: `LDR h, ADD h+a, AND mod R, MOVZ 6912, MADD, LDR k_n, MOVZ 576, MADD, MOV x1,y, EOR x0,y,#0x400` | 10 |
| plumbing, accept tail: `AND h',w1,#0x3FFF, STR h', STR k_n, MOVZ a,#0, EOR y,y,#0x400` | 5 |
| plumbing, probe ladder (reject): `ADD a,#1, CMP a,#8, B.LT` | 3 |
| selector operand base x2: `EOR x29,x20,#0x400` | 2 |
| 3 `blst_fp_mul` arg setups x 3 `ADD` | 9 |
| trail slot counter w/ 126-bit bound | 7 |
| global slot counter w/ 126-bit bound | 7 |
| next-slot branch | 1 |
| **subtotal** | **1,743** |
| less probe ladder (no accepted slot runs it) | -3 |
| **maximum over slot types** | **1,740** |
| no explicit margin: charged at the executed maximum | 0 |
| **charged per slot** | **1,740** |

Slot = `1 + 1740/20000 = 21740/20000` units. Kernel contains operand loads/stores; table reads not re-charged. Rows 1,743, max slot 1,740, charged 1,740 (3 below row sum). Ceiling 1,743: at 1,744 the least thousandth becomes 125.819 (exact integer search); charge sits 3 below. Sensitivity `M_max/(20000 W ln2) = 6.636e-5`/instr at 1,740; `1/(20000 ln2) = 7.2e-5` upper bound (`ledger`).

**Row sequences.** Executed on register-file machine `A64` (counts + Python reference): modular ops, extrema + indices, 11th compare, DP test, scan + stub, ring write, seed, slot update. Priced from written-out sequences: 3 `fp_mul` (1,149), `khash` (112), bookkeeping (12), plumbing (10,5,3,2,9,7,7,1). Mod-sub `x-y mod p`: 6+6 `LDR`, six-word `SUBS/SBCS` in place (carry 1 iff `x >= y`), `B.CS` past correction, borrow path `ADDS`+5`ADCS` of `p` (x12-x17), 6 `STR`: 31 correcting / 25 skipping, so 31 bounds; inline, no call setup. Carry only through `ADCS`; result canonical. No borrow mask (branch replaces it): writes only x0-x17, `p` registers intact. `p` = 6 `LDR` from frame; no register holds it across a block (kernels clobber x0-x18, x19-x28 fixed, extrema take x0-x18+x29+x30): 3 blocks x 1 load, shared inside block: `a-b` before kernels (input to 3rd); two `zeta^2` terms after (one shared load); negation after extrema. Negation `p-M`: same shape, `XZR` subtrahend, chain in x10-x15, `B.CS` (taken iff `M = 0`, leaving `0`), no store: 13 (19 with load). Row `3*31+13+18 = 124` (= `37+68+19`). Selector compares: section 3. `khash`/word: 6 x (`EOR` + 8-instr mix) = 54, x2 = 108; limbs already in x4-x9 (final `CSEL`); +2 seed +2 multiplicand reloads = 112. DP test reads all 64 bits `w_1` + low 16 `w_2` = 80 bits; 32-bit form would be 48-bit (`S_max 2^-48 = 2^77.7` reports vs `D_max = 2^47`; witness `w_1 = 2^32`, `w_2 = 0`). Membership slot: `LDR x10,[sp,#8i], CMP x10,x1, B.EQ stub` (3; `8*23 = 184` in scaled range). Stub: `ADD x11,sp,#(K_OFF+48i)` (`K_OFF+48*23 = 1296`, 12-bit range), 6 `LDR`, `6 EOR+5 ORR`, `CBZ`, end-branch: 20. Exact (canonical `< p < 2^381`). <= 1 stub/candidate: row 92 (miss 72, last-slot hit max; equality + collision paths; asserted by `hash_scan`, `ring_check`). Ring write: `48p` via 3 adds, 3 `STP`, hash `ADD+STR`: 8. (`STR w1,[sp,p,LSL #3]` unencodable — 32-bit offset allows `LSL #0/#2` only — and would store half a word, letting the 64-bit scan miss.) Bookkeeping: 126-bit counter `ADDS/ADCS`, `SUBS/SBCS` vs `S_max` + branch (5) + 2 `LDR` (7); slot `ADD/CMP/CSEL` (3); 2 branches: 12. Plumbing in loop order as tabled (dense kernel `blst_fp12_mul(x0,x1,x2)`, clobbers x0-x18; accept tail; reject ladder; `EOR x29` x2 (scratch); arg setups; two 126-bit counters; closing branch). Accepted-step test + update rows charged every slot (over-charge over `2^-17.44` rejects).

**Window structure.** Set of last `<= 24` class keys incl. current; oldest dropped/accept. Two 24-arrays: `K[0..23]` keys (48 B, ring; k-th key -> slot `k mod 24` = oldest slot, register `p`, bounded `ADD/CMP/CSEL` advance); `H[0..23]` stored `w_1` words (index + DP test already compute them: zero hash evaluated, one word stored/accept). Start seeds all 24 slots with `kappa_0` + word, `p = 0`: 121 instr (in start slack). After `k` accepts: stored set `{kappa_0..kappa_k}` (last `min(k+1,24)`), no stale keys/sentinels; `k >= 23` overwrites `kappa_{k-23}`; patterns `l <= 24` return to stored keys (section 7).

**Membership = conservative, bounded.** Repeats always caught (repeat carries a stored word; scan reaches it, rejects); section-7 patterns never missed. False reject only on 64-bit collision (H2: uniform indep.): `<= 24*2^-64 < 2^-59.4`/candidate, `<= 8*24*2^-64 < 2^-56.4`/accepted step. Collision escape = ordinary escape cost (probe rejected, mult charged); cannot lose repeats or mutate keys. Escape rate unchanged `nu <= 1.1/(12R) < 2^-17.44` (collision = relative `2^-38.97` of two-cycle term, inside 10% allowance; section 9 measures this rule).

**Why 64-bit words, not a mask.** 24-bit mask over 384-bit keys: escape rate ~0.64/step + bucket collisions. Stored words: no false negatives, `24*2^-64`/candidate, exact test on hit. Sorted 1-byte index array answers an unasked insert-position question: 259 vs 100 instr here.

**Register allocation (complete map; every spill paid).**

| reg | holds | live across |
| --- | --- | --- |
| x19 | `&TW` table base | loop, all calls |
| x20 | y buffer ptr (other = `x20^0x400`; 2 KiB-aligned frame) | loop |
| x21-x22 | global slots `M` | loop |
| x23-x24 | trail slots `c` | loop |
| x25-x26 | accepted `S` | loop |
| x27 | ring slot `p` | loop |
| x28 | probe `a` (0..7) | state |
| x29 | scratch (tail `i_M+6`, then y' base) | block |
| x30 | scratch (link reg saved in prologue) | block |
| x0-x18 | kernel args + phase scratch (`kappa` x4-x9, `w_1` x1, `w_2` x2 after selector) | phase |
| sp+0 | `H[24]` | trail |
| sp+192 | `K[24]` (`K_OFF+48*23 = 1296`, 12-bit `ADD`) | trail |
| sp+1344 | `V[6]` selector values | state |
| sp+1632 | modulus `p` limbs | block |
| sp+1680/1728 | `zeta`, `g0` | trail |
| sp+1776/1784 | `h`, `k_n` (head `LDR`, tail `STR`) | state |
| sp+1808/1824 | `mix(salt)` seeds, splitmix64 mults | trail |
| sp+1840/1856/1872 | `L_cap`, `M_max+L_cap`, `S_max` | trail |
| sp+1888 | `a-b` scratch (3rd `fp_mul` input) | state |
| sp+2048/3072 | two 576 B y buffers (`0x400` apart) | trail |

| phase | scratch | regs | slot loads/stores |
| --- | --- | --- | --- |
| entry address + args | `h,a,k_n`, 2 strides, addr, 3 ptrs | 8 of x0-x18 | `h,k_n` (2), strides (2), `MOV x1` (1) |
| `blst_fp12_mul` | 3 arg ptrs | 3 | none (x19-x28 survive) |
| sub blocks (i),(ii) | `x`(6) x0-x5, `y`(6) x6-x11, `p`(6) x12-x17, base x29 | 18+x29 | `p` x6 in (i), x6 shared in (ii) |
| sub block (iii) neg | `M`(6) x10-x15, `p`(6) x0-x3,x16-x17 | 12 | `p` x6 |
| `fp_mul` setups | 3 ptrs/call x3 | 3 | none (9 `ADD`) |
| extrema seed | min(7), max(7), cand(6), base | 20+x29 | `V_0` from base, `V_1-V_5` frame |
| extrema loop | min(7), max(7), cand(6), supply | 21 | none |
| `khash` | `kappa`(6), `w_1,w_2`, temp, 2 mults | 11 | seeds (2), mults (2) |
| scan + stub | `kappa`(6), `w_1,w_2`, x10, x11, 6 `EOR` | 16 | none |
| ring write | `kappa`(6), `w_1`, x10, x11 | 9 | none |

**Subtraction registers, executed.** Subs write x0-x17 only; negation x10-x15 + `p` regs; never x19-x28/x29 (`B.CS` replaces mask). `modsub_block_check` seeds x19-x29 from allocation, asserts no change outside x0-x17, `p` intact, only 6 dest words written. Block (ii) executed as priced: 2 subs under ONE `p` load, <= 68/pair, results vs Python ref. Negative control: `CSETM x19` + `p` in x20-x25 caught. Both branch paths exercised (176 tuples: 68 correcting, 108 skipping); block (i) 37/31 with load.

126-bit constants, `p` limbs, `h`, `k_n`: no callee-saved reg (x19-x28 full, kernels clobber x0-x18) -> frame words, every load paid (counter rows' `2 LDR`, selector's 18, `h`/`k_n` `LDR` + tail `STR`). Only state surviving a reject outside registers = those two words, charged every path. No other spills (candidate in other y buffer, x20 toggled on accept only).

**Two hash words.** Index + DP from class key only: `w_1,w_2`, index `w_1 mod 2^14`. DP iff `w_1 = 0`, `w_2 mod 2^16 = 0` (80 bits). DPs reported, never stepped from: stepped-from low-14-bit values uniform (`2^-14(1+-2^-50)`). Both words charged/slot (over-charges `2^-17.44` rejects).

**Main walk.** `N = (r-1)/12`, `mu = 11/10`, `x = 3903/10000`, `L = ln(10000/6097)`; `U = sum_{k=1}^{64} x^k/k + x^65/(65(1-x))`, exact rational upper bound, `0 < U-L < 2^-99` (gap `6.63e-31`; rational comparison vs 70-term partial sum). Budget for main-event `0.3903` (rule `0.39` + `2.7e-4` loss): least `S` with `S^2 >= 2 mu N U`, exact integers:

    S_max = ceil(sqrt(2 mu N U)) = 68967528808679891282472564154805017972 (log2 125.69726)
    M_max = S_max + floor(S_max/2^13) + 1 = 68975947696473919589513881606484071319, work 21740/20000 M_max.

No float/decimal in `S_max` (`S_max^2 >= 2muNU > (S_max-1)^2` integer compare). `L`-ceiling differs (`U/L-1 = 2^-99.24`, ~`2^25.46` steps lower); `U` costs relative `<= 2^-99`. All inequalities on this printed `S_max`; exact checks `W^1000 <= 2^125818`, `W^1000 > 2^125817`.

Probe allowance `S_max/2^13` (`2^-13` of accepts); escape rate `nu <= 1.1/(12R) < 2^-17.44` (L2; measured `5.57e-6`); rejects binomial-dominated, mean `< 2^108.25`, allowance `2^112.7` = 21x mean, tail `< exp(-2^107)`; meter fails on excess (L7). Hard cap, not expectation.

**Additive items** (cap-fixed hard bounds):

| item | amount (units) | log2 |
| --- | --- | ---: |
| two pairings (`g` + spare) x 2^30 instr | `2^31/20000` | 16.72 |
| table `R(2*509+1+11+slack) <= R*1112` | | 24.12 |
| starts `ST_max*1100` (1,019 mults + ChaCha20 + reductions < 0.4) | | 57.11 |
| reports `D_max*2^16/20000` (RB-tree + transfer; op <= 2,300 instr) | | 48.72 |
| replays `REC_max*2*L_cap*5/4` (slot 21740/20000 + exponent < 0.075: two 255-bit mod-prods by stored `lambda^-k` + adds, <= 1,500 instr, `(u_h,v_h)` from 2^20-B array; `1.16200 <= 5/4`) | | 92.33 |
| last trail `L_cap*21740/20000` | | 87.13 |
| `[a]Q0`, `g^a == z` | `<= 2^12` | 12.00 |
| **additive total** | | **92.361** |

**Total.** `W_max = 21740/20000 M_max + 2^92.361 <= 2^125.818` (`log2 = 125.81778482...`); claim `2^125.818` = least thousandth above, proved by exact integers: `m` least with `W^1000 <= 2^m` (binary search `[0,1000*bitlength(W)]`): `W^1000 <= 2^125818` and `W^1000 > 2^125817` (no smaller thousandth). No float. Randomness (starts/table), sorting/lookup/transfer (reports), verification + output all charged. No amortization.

**Advice:** none; preprocessing in section 8.

## 6. Success probability

*Model (H1).* Accepted-state classes `X_1,X_2,...` of successive trails: fresh successors uniform up to loss `mu`: no two-trail class repeat after `s` accepts with prob. `<= exp(-s^2/(2muN))`. Ideal map `mu = 1`, `E[S^2/(2N)] = 1`; section 9: `mu_hat = 0.996+-0.014` (`R >= 256`) vs `mu = 11/10`. `S_max^2/(2muN) >= U > ln(10000/6097)`: **main event** (two-trail class repeat within `S_max`) prob. `>= 3903/10000 = 0.3903`.

*Repeat -> recovery.* Trails sequential; earlier trail ended. First class the later trail shares with an earlier one: earlier trail reported (else L1/L3), later trail heads to its DP, coalesces (4.3) to the same DP; store holds the key; Recover runs, denominator vanishes w.p. `1/r`.

*Losses* (upper bounds on main event not yielding recovery):

| label | cause | bound |
| --- | --- | ---: |
| L1 | trail cut by `L_cap` (no DP in `2^87`: `e^-128 = 2^-184.7` later; length-biased `128e^-128 = 2^-177.66` earlier) | `<= 2^-176.6` |
| L2 | escape of a meeting trail within 24 steps after meeting (`1.1/(12R) < 2^-17.44` incl. `8*24*2^-64 < 2^-56.4` hash term) | `<= 2*24*(1.1/(12R)) <= 2.7e-4` |
| L3 | residual trap / within-trail rho cycle / 8 probes rejected (`2^80 2^-163.5`, `E[L^2]/(2N) = 2^-91`, `2^80 nu^8 <= 2^-59.5`) | `<= 2^-58` |
| L4 | `V_k` ties (8 linear conditions, `8/p < 2^-377.6`/state) or `c1 = 0`; identity `1/r < 2^-254.84`/slot; over `2^126` slots | `<= 2^-128.84` |
| L5 | false store match / zero denominator | `<= 2^-250` |
| L6 | sampling distance (`2^47` starts + `2^14` table scalars x `2^-129`) | `<= 2^-80` |
| L7 | caps `D_max,ST_max,M_max` exceeded (Poisson, mean `2^45.7`, `<= 2^46`) | `<= 2^-1000` |
| L8 | `> REC_max` unusable matches | `<= 2^-200` |

Total loss `< 2.7e-4` (only L2 above `1e-17`). Success `>= 0.3903-0.00027 = 0.39003 > 0.39` under H1-H3 (margin `0.00003`).

*H1's load + margin.* Claim uses H1 at `s = S_max` only: first repeat by `S_max` w.p. `>= 0.3903`. `>= 0.39` after losses needs `0.39027`: true inefficiency `mu_true <= 1.1 ln(1/0.6097)/ln(1/0.60973) = 1.10011`; claim spends the `11/10` allowance, keeps `0.01%` spare. Measured `0.996+-0.014` is 7.4 s.e. below; ideal-map value at `S_max`: `1-(6097/10000)^1.1 = 0.4197` (section 9 reports the statistic, not just the mean).

*Uniform DP/routing fields (H2).* Used structure: `w_1,w_2` (hence index, DP predicate, window words + `24*2^-64` term) and `k*` act as indep. uniform fields on distinct classes. `k*`: exact for uniform GT (`y -> sigma(y)` bijects, shifts `k*` by 1: uniform, class-indep.). `khash`: standard hash heuristic — indep. uniform 64-bit words, index uniform on `R`, predicate density `2^-80`. Only correlation (DP has index 0) is off-path (DPs never stepped from, section 5). Production evidence: section-3 algebra checks only. Density within 2x of `2^-80` only moves report counts (`D_max = 2^47` vs mean `2^45.7`: up to `2.3*2^-80` fits).

*H3.* Fresh-key ChaCha20 = ideal coins.

## 7. Deterministic traps and the look-back window

*Trap* = deterministic closed pattern (class path returns after `l` accepts). Lifted exponents (`x_n = alpha_n+a beta_n`, `m_h = u_h+a v_h`, step `i` uses index `u_i = h_i+a_i`, twist `k_i`): return to class of `x_t` after `l` steps means `lambda^s x_t = x_t + sum_{i<l} lambda^{-k_i} m_{u_i}`, some `s in Z/12`. `s != 0` + formal-vanishing sum -> `x_t = 0`, excluded; else below. Probes included (lemma on *used* `u_i`).

*Formal part.* `m_h` iid uniform. `P_h = sum_{i:h_i=h} zeta_12^{-k_i}` in `Z[zeta_12]`, coeffs `<= l` in abs. Ring map `zeta_12 -> lambda (mod r)`: kernel a prime ideal of norm `r` (`r = 1 mod 12` splits); nonzero kernel `alpha` has `r | N(alpha)`; coeffs bounded by `B` in `1,zeta,zeta^2,zeta^3` give `|N(alpha)| <= (4B)^4`, forcing `B >= (r/256)^{1/4} > 2^61`. So for `l < 2^60`, nonzero `P_h` has nonzero invertible image mod `r`; `sum_h P_h(lambda) m_h` uniform once some `P_h != 0`; return = one linear equation, prob. `1/r`/pattern: ordinary random revisit (H1; window revisit `~24^2/N < 2^-240`, not walk structure). `l >= 2^60`: necessary condition only; non-formal returns are within-trail rho collisions (`<= L_cap^2/(2N) <= 2^-78`/trail at cap, `E[L^2]/(2N) = 2^-91` real distribution, in L3). Deterministic trap needs `P_h = 0` for **every** used `h`; single roots nonzero, so **every used index appears >= 2x**. (Minimal vanishing 12th-root sums: antipodal pairs = two-cycle with twists differing by 6; triples at 120°; classification unneeded.)

*Probability.* Used indices iid uniform on `R` for distinct classes (H2; `u_i != h_i` only on reject, `<= 1.1/(12R) < 2^-17`, factor `(1+2^-12)^l`, negligible). `D` distinct indices: trap needs fibres `>= 2` (`D <= l/2`) and per-index twists summing to 0; last twist of a fibre forced given the rest: `<= 1/12`/fibre, indep. (twists = partial sums of uniform iid shifts):

    P(trap of length l at a given step) <= E(l) = sum_{d=1}^{min(floor(l/2),R)} C(R,d)(d/R)^l 12^-d.

`ledger` evaluates in float log scale (printed `2^-164.2`, `2^-163.5` are float evaluations of the exact integer sum; only the claim inequality is exact-compared). Max over `25 <= l <= 512`: `2^-164.2` (`l = 25`); sum `2^-163.5`. `l > 512`: `C(R,d) <= (eR/d)^d` -> terms `<= ((e/12)(d/R)^{l/d-1})^d`; `d <= R/2`: `<= 0.2266^d 2^{-(l-d)} <= 2^-l`; `d > R/2`: `<= 0.2266^{R/2} < 2^-17000`. So `E(l) <= 2^-480` (`l > 512`), sum to cap `< 2^-400` « `2^-163.5`. `l <= 24`: window detects within `l` steps, escapes (hash rule never misses window-key repeats); unconditional escape/step `nu <= 1.1/(12R)` (two-cycle dominates: same index `1/R`, twists differ by 6 `1/12`; hash adds `8*24*2^-64 < 2^-56.4` = relative `2^-39`), measured `5.57e-6` vs `1.1/(12R) = 5.59e-6`, `1/(12R) = 5.09e-6` (section 9); lengths 3-24 are `O(R^-2)`, inside factor 1.1.

*Undetected-trap cost.* Mean-`2^80` trail trapped w.p. `<= 2^80 sum E(l) = 2^-83.5`, then runs `L_cap = 2^87` (128 mean lengths): wasted main-walk fraction `<= 2^87 2^-163.5 = 2^-76.5`. Trail dies (8 rejects) w.p. `<= nu^8 <= 2^-139.5`/step (8 distinct indices; two-cycle hits at most one; H2 independence); per trail `2^80 nu^8 <= 2^-59.5`.

## 8. Memory, preprocessing, advice

| item | bytes (bound) |
| --- | --- |
| store `<= D_max = 2^47` x 96 B (48-B key, 16-B index, 3 ptrs, colour) | `2^53.585` |
| twist table `12*2^14*576` | `2^26.75` |
| replay scalars `(u_h,v_h)`: `2^14*64` | `2^20` |
| trail state (`y` + probe buffer 2x576 B, window 1,152+192 B, selector/hash buffers, counters, ChaCha20) | `< 2^20` |
| **peak** | `<= 2^53.6 < 2^54` |

One address space `< 2^64`; store transfers inside 2^16-instr report cost. Preprocessing (table + pairing) `2^24.13` units. No advice.

## 9. Experiments and their scope

One declared experiment (`hashring-rho-checks`, one initial run) runs four parts in one process (`{"mode":"all",...}`, ~10 s laptop; callable alone by `mode`). None is production-scale; none measures hardware counts or production success; counts below are of written-out sequences executed in Python.

* **ledger** recomputes sections 5, 7, 8 in exact/`Decimal` arithmetic where exact (`S_max`, `M_max`, additives, `log2 W_max = 125.8177848249 <= 125.818`, exact checks `W^1000 <= 2^125818`, `W^1000 > 2^125817`, max slot 1,740, charge 1,740, ceiling 1,743, `2^53.585` B, `2^24.127` preprocessing), trap bounds in float log scale (section 7), min-max network vs two scans + direct 12-scan on 117,649 tuples. `selector_min_check` + `ring_check` execute priced AArch64 lists on `A64` vs Python ref: selector tail (218/tuple, 212 when `M = 0`; 245 small-field — 51 ties + all-zero — + 40 production tuples), mod-sub in all 3 blocks (31 correcting, 25 skipping, 37/31 with `p` load), two-sub composition under ONE `p` load (<= 68) with x19-x28/x29 byte-identical, `p` regs unchanged, mask-shape negative control caught, DP test (3; 314 pairs; `w_1 = 2^32` 80-vs-48 witness), slot update x200 (`p = k mod 24`, oldest slot), seed (121), ring write (8), scan + stub (72 miss, 92 slot-23 hit, equality + collision paths) over fills 0..24, 676 candidates vs plain set, 3,412 key pairs vs integer equality. Plumbing (10,5,3,2,9,7,7,1), bookkeeping (12), kernel (1,149), `khash` (112) priced, not executed (section 5 says so).
* **algebra** (real BLS12-381): `lambda` order 12, `lambda^6 = -1`; `zeta` primitive 6th root in Fp; `sigma(y) = y^p = y^lambda` on random GT; `g1 = -g0`, `g0^6 = -1/8`, `zeta = 2g0^2`; selector's 12 values = true first coordinates; 3-mult form = general `g0 a+g1 b` form; two-half min = direct 12-scan (element + 12 conjugates); `kappa` invariant, `k*` covariant; twisted step class-equivariant. Exact section-3 checks on production data (2 elements default).
* **toy**: full algorithm (selector, twisted table, window + probing, DPs, store, replay + exponents, Frobenius alignment, `g^a == z`) over production-shape toy tower: 60-bit `q`, `q mod r' = lambda'` order 12, `q = 7 mod 12`, `1+u` non-square/cube in Fq2, 22-bit prime subgroup `r'`, real Fq12, 64-entry table, window 12, mean-8 trails. Default (6 instances, seed 3): 6/6 logs = planted, 0 wrong, replay = untracked paths, 5 escapes, 0 dead, mean `1.578 sqrt(N')`. Offline (`bits 22, trials 80, seed 77`): 80/80, 0 wrong, 41/80 (51.3%) in toy budget (ideal 41.97%, req. 39%), mean `1.219 sqrt(N')`. Towers built with `gamma = g0(1-u)`: same 3-mult selector, 12 values vs true coordinates/instance. Window = production ring (keys + hash words, same `w_1`, update, scan, six-limb test), every decision vs plain list.
* **sim**: same walk in exponent domain, cyclic order `r'`, automorphism `lambda'` order 12 (state = dlog; 64-bit hashed keys of 12 conjugates); `R = 2^14`, window 24, probing 8; short trails keep within-trail rho negligible. Measures `mu_hat = mean(S^2/(2N))` at first inter-trail repeat, plain-birthday control, escape rate, merged-trail same-DP rate, fraction with `S^2/(2N) <= 1.1U = 0.54427` (in `S_max` budget). Default (600 groups, 22 bits, mean 8, seed 3): `mu_hat = 1.028+-0.041`, fraction `0.4167`, control `1.058`, 599/600 same DP. Escape `4.4e-5` vs `1/(12R) = 5.1e-6`: at `N = 2^18` the window catches random revisits (`~24/N`/step); at 28-34 bits that term `< 1e-6`, rate matches `1/(12R)`.

Offline larger runs (same code, `{"mode":"sim",...}`):

| bits r' | R | trials | seed | `mu_hat` | s.e. | frac in budget | control | escape / `1/(12R)` | same DP |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: |
| 28 | 64 | 1500 | 105 | 1.072 | 0.026 | 0.4053 | 1.085 | 1.26e-3/1.30e-3 | 1493/1500 |
| 30 | 256 | 1500 | 104 | 0.992 | 0.026 | 0.4073 | 1.028 | 3.10e-4/3.26e-4 | 1500/1500 |
| 30 | 1024 | 1500 | 103 | 0.999 | 0.026 | 0.4173 | 0.935 | 7.68e-5/8.14e-5 | 1500/1500 |
| 30 | 16384 | 1500 | 101 | 1.008 | 0.026 | 0.4180 | 0.968 | 5.57e-6/5.09e-6 | 1500/1500 |
| 32 | 16384 | 700 | 102 | 0.969 | 0.038 | 0.4386 | 1.041 | 4.43e-6/5.09e-6 | 700/700 |
| 34 | 16384 | 250 | 106 | 0.995 | 0.063 | 0.4320 | 1.005 | 5.19e-6/5.09e-6 | 250/250 |

Pooled `R >= 256` (5,450 repeats): `mu_hat = 0.996+-0.014`; pooled in-budget fraction `0.4182+-0.0067` (ideal `0.4197`, req. `0.3903`, 4.2 s.e. margin). Control = same estimator on plain birthday over `N` (`min(trials,300)`/row, s.e. `0.058`, `0.063` last row, range `0.935...1.085`). Claim `mu = 11/10`: 10% over measured `0.996` (7.4 s.e.). Long-trail inflation is within-trail rho, as expected: `sim(26,16384,300,5,db=8)` (mean `2^8`) gives `mu_hat = 1.887` (294/300 same DP), falling as `1/N'`; at production `E[L^2]/(2N) = 2^-91` (`L = 2^80`, `N = 2^251`), so rows use short trails (mean `2^5`).

**Extrapolation to production** is H1's principal uncertainty: emulation replaces Fp12 structure by pseudorandom keys; `N = 2^251.3` far outside any run. Both emulations draw index + DP bits from independent salts (idealised production rule: index `w_1 mod 2^14`, 64 of 80 bits `w_1 = 0`). Algebra exact; dynamics disclosed heuristic.

## 10. Limitations and disclosed heuristics

* H1, H2, H3 as above; success `>= 0.39` conditional on them.
* `mu = 11/10` empirical margin, not proof (`mu_hat = 0.996+-0.014` would support less). Budget at rule's requirement (main event `0.3903`, success `0.39003`): survives true inefficiency to `1.10011`.
* `2^47`-record store (`2^53.6` B) a mathematical object in one 64-bit address space.
* Prices are instruction estimates on pinned kernel prices; selector/hashing/look-back rows are upper bounds, charged at exactly the 1,740 max slot with no margin: every covered row executed on `A64` or a priced written-out sequence (section 5). Third multiplication rests on exact tower identity `gamma` (section 3).
* Window over-rejects on 64-bit collisions: charged escape, cannot lose repeats (section 5). Claim survives 3 further instr/slot above charge, 0 above row sum.
* Not a break: generic-group rho in GT on Frobenius classes, `2^125.818` units.

## 11. Reproduction

    cd submissions/bls12-381-fixed-g1-invert-g2/hashring-rho-slim
    echo '{"mode":"all","seed":3,"algebra_elements":2,"toy_bits":22,"toy_trials":6,"sim_bits":22,"sim_R":16384,"sim_trials":600,"sim_db":3}' | python3 -I experiment.py
    echo '{"mode":"ledger"}'                                              | python3 -I experiment.py
    # scripts/ lives in the checkout, not in this package: run that line from the checkout root
    cd ../../../..
    bash scripts/run_experiments.sh bls12-381-fixed-g1-invert-g2 submissions/bls12-381-fixed-g1-invert-g2/hashring-rho-slim --preflight-only

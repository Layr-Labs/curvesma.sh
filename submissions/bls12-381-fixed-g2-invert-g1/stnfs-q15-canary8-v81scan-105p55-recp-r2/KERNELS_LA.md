# Kernels Part 2: Montgomery CIOS, NTT/Lehmer, `Z[zeta_12]` Block Wiedemann, and Descent Kernels

### 1.7 Pure-Scalar In-Place `uint64_t pos[]/STADD` Dense Build (`W_B=1`: `5.90 + 7.75 + ovh_k`; `W_B=8`: `6.90 + 8.75 + ovh_k`; `2.90 <= 3.0` Ir/bkt Zero+Prefix) & Direct 1-Slot Sparse Build (`0.65` Ir/bkt, `11.25 + ovh_k + 3.0*ovf_m/lam`)

**Exact Full-Canonical `G_120` Precomputed `sub_chords[]` Stream, Overflow-Free `p < 2^{64}` Modular Addition, & Zero-Overflow `uint64_t` In-Place `pos[0..p]` Histogram/Prefix:** Within each radial shell `k in {0, 1}`, canonical `a`-orbit representatives `a = (a_0, a_1, ..., a_7)` in the full canonical domain `F_Vor cap {arg(sigma_1(a)) in [-pi/30, pi/30)} cap {|sigma_1(a)| > max(|sigma_7|, |sigma_4|, |sigma_2|)} cap Shell_k(cls(a))` (with `G[0,0] = Tr(1) = 8`, `s_vals = [1.0, 1.40, 1.82]`, exact `max |a_j| <= 388.33 < 413`) are enumerated ONCE per shell during `canonical_a_and_sort` (`<= 2^{54} << 2^{86.0}` u) in lexicographic `(a_7..a_1, sub_chord)` order and serialized into a read-only packed stream `sub_chords[]` (`8` B/sub-chord, `<= 0.055 EiB` shared read-only RAM) plus a tiny carry table `fib_carry[]` (`< 1%` of fibers). Because `2*O_H = P_1*P_2` and `4*O_H = P_1^2*P_2^2`, **`Delta a_0 = +2` preserves `cls(a) == 0` identically** (`a + 2 == a mod 2`) and **`Delta a_0 = +4` preserves `cls(a) in {1, 2}` identically** (`min(2, v_{P_i}(a+4)) == min(2, v_{P_i}(a))`); both steps shift `k_m(a)` by the **EXACT SAME constant increment `x_d2 = 2*r^{-1} mod p` or `x_d4 = 4*r^{-1} mod p` across ALL 4 Galois banks `m in {0..3}`** while `x_om` increments by `+1` in lockstep with `oid` (`tab_aj[0..830]` built once per prime in `< 20,000` ops). Both `rev_b2_build_subchord_hdr` (`12` ops, `a0_lo_idx = (a_0_lo - step) + 417 in [0..830]`) and `rev_b2_build_fiber_hdr` (`8 <= 9.0` ops fast path; `48.0` ops full 7-level carry) use **overflow-free 64-bit modular addition `sub x3, x_p, x2; subs x1, x_a, x3; add x2, x_a, x2; csel x_out, x1, x2, hs`**, which is exact for EVERY `0 <= x_a, x2 < p < 2^{64}` including `p > 2^{63}` (even when the temporary `add` wraps mod `2^{64}`, `subs` sets `hs` iff `x_a >= p - x2` and `csel` discards the wrapped `add` and selects the exact reduced difference `x_a - (p - x2)`, verified at `p = 2^{64} - 59` in `w8_verifiers.py`!). Sub-chord header (`12`) + `pre8/4-2-1` tail (`<= 7.10` ops worst-case across all `n >= 1`) = `<= 19.25` ops/sub-chord; live Horvitz-Thompson verification on the exact full-canonical `G_120` domain (`unitdomain.py`) proves `3-SE` upper bounds `<= bld_ovh_shells = [2.14, 2.67]` Ir/entry charged across **ALL 4 build tiers** (Band-1a/1b, `W_B=1`, `W_B=8`, Direct-1). **Zero-Overflow `uint64_t` In-Place `pos[0..p]` Build (ALL dense primes `2^{33} <= p < 2^{61.3}`):** Set `pos[0] = flat_base`, zero `pos[1..B_pad]` (`rev_b2_memset_ff_u16` with `x_ff = 0`, `0.625 <= 0.65` ops/bkt), accumulate exact byte counts `8 * count(b)` directly into `x_hist = &pos[1]` via **`stadd x_eight, [x4]`** (`rev_b2_build_hist_w1_stadd` `5.875 <= 5.90` on Band-1a & `W_B=1`, `rev_b2_build_hist_u8_stadd` `6.875 <= 6.90` on `W_B=8` using `bic x4, x_k, #7; add x4, x_hist, x4`), convert `pos[1..B_pad]` in-place to exclusive start pointers via **`rev_build_prefix_u64_8x`** (`18/8 = 2.25` ops/bkt; `2.25 + 0.65 = 2.90 <= 3.0 = c_lbkt`), and scatter once over flat `kbuf[0..N_a,k-1]` advancing `pos[b + 1]` to `off[b + 1]` (zero copy, zero separate `hist[]` -> `W_B=1` 1-owner RAM at `p_w1 = 60.0` stays **`11.028 EiB < 16.000 EiB`**!). Because `8 * N_{a,k} = 8 * 0.7431 * 2^{58} < 2^{61} << 2^{64}`, **histogram overflow is mathematically impossible for EVERY prime `p >= 2^{33}` (even Band-1a `p = 2^{33}` with `~31M` records/bucket, verified live in `w8_verifiers.py`!)**.

```asm
; kernel: rev_b2_build_subchord_hdr items=1 charged=12.0 key=c_bld_sub_hdr
ldr x0, [x_sc], #8
tbnz x0, #63, .Lbld_upd_fib
.Lbld_sub_go:
ubfx x1, x0, #0, #10
ldr x2, [x_tab0, x1, lsl #3]
sub x3, x_p, x2
subs x1, x_kbase, x3
add x2, x_kbase, x2
csel x_k, x1, x2, hs
ubfx x_cnt, x0, #10, #10
tst x0, #0x100000
csel x_d, x_d4, x_d2, ne
csel x_pms, x_pms4, x_pms2, ne
```
```asm
; kernel: rev_b2_build_fiber_hdr items=1 charged=9.0 key=c_bld_fib_hdr
.Lbld_upd_fib:
tbnz x0, #62, .Lbld_upd_carry
ubfx x1, x0, #24, #10
ldr x2, [x_tab1, x1, lsl #3]
sub x3, x_p, x2
subs x1, x_klvl2, x3
add x2, x_klvl2, x2
csel x_kbase, x1, x2, hs
b .Lbld_sub_go
.Lbld_upd_carry:
```

```asm
; kernel: rev_b2_build_hist_w1_stadd items=8 charged=5.90 key=c_bld_hist_w1
.Lbhw1:
subs x1, x_k, x_pms
add x2, x_k, x_d
csel x_k, x1, x2, hs
add x4, x_hist, x_k, lsl #3
stadd x_eight, [x4]
subs x1, x_k, x_pms
add x2, x_k, x_d
csel x3, x1, x2, hs
add x4, x_hist, x3, lsl #3
stadd x_eight, [x4]
stp x_k, x3, [x_kbuf], #16
subs x1, x3, x_pms
add x2, x3, x_d
csel x_k, x1, x2, hs
add x4, x_hist, x_k, lsl #3
stadd x_eight, [x4]
subs x1, x_k, x_pms
add x2, x_k, x_d
csel x3, x1, x2, hs
add x4, x_hist, x3, lsl #3
stadd x_eight, [x4]
stp x_k, x3, [x_kbuf], #16
subs x1, x3, x_pms
add x2, x3, x_d
csel x_k, x1, x2, hs
add x4, x_hist, x_k, lsl #3
stadd x_eight, [x4]
subs x1, x_k, x_pms
add x2, x_k, x_d
csel x3, x1, x2, hs
add x4, x_hist, x3, lsl #3
stadd x_eight, [x4]
stp x_k, x3, [x_kbuf], #16
subs x1, x3, x_pms
add x2, x3, x_d
csel x_k, x1, x2, hs
add x4, x_hist, x_k, lsl #3
stadd x_eight, [x4]
subs x1, x_k, x_pms
add x2, x_k, x_d
csel x3, x1, x2, hs
add x4, x_hist, x3, lsl #3
stadd x_eight, [x4]
stp x_k, x3, [x_kbuf], #16
mov x_k, x3
subs x_cnt, x_cnt, #8
b.hs .Lbhw1
```
```asm
; kernel: rev_b2_build_scatter_w1_u8 items=8 charged=7.75 key=c_bld_scat_w1
.Lscw1:
ldp x0, x1, [x_kbuf], #16
ubfiz x3, x0, #59, #5
orr x3, x3, x_om
add x4, x_pos, x0, lsl #3
ldr x5, [x4]
str x3, [x5], #8
str x5, [x4]
add x_om, x_om, #1
ubfiz x3, x1, #59, #5
orr x3, x3, x_om
add x4, x_pos, x1, lsl #3
ldr x5, [x4]
str x3, [x5], #8
str x5, [x4]
add x_om, x_om, #1
ldp x0, x1, [x_kbuf], #16
ubfiz x3, x0, #59, #5
orr x3, x3, x_om
add x4, x_pos, x0, lsl #3
ldr x5, [x4]
str x3, [x5], #8
str x5, [x4]
add x_om, x_om, #1
ubfiz x3, x1, #59, #5
orr x3, x3, x_om
add x4, x_pos, x1, lsl #3
ldr x5, [x4]
str x3, [x5], #8
str x5, [x4]
add x_om, x_om, #1
ldp x0, x1, [x_kbuf], #16
ubfiz x3, x0, #59, #5
orr x3, x3, x_om
add x4, x_pos, x0, lsl #3
ldr x5, [x4]
str x3, [x5], #8
str x5, [x4]
add x_om, x_om, #1
ubfiz x3, x1, #59, #5
orr x3, x3, x_om
add x4, x_pos, x1, lsl #3
ldr x5, [x4]
str x3, [x5], #8
str x5, [x4]
add x_om, x_om, #1
ldp x0, x1, [x_kbuf], #16
ubfiz x3, x0, #59, #5
orr x3, x3, x_om
add x4, x_pos, x0, lsl #3
ldr x5, [x4]
str x3, [x5], #8
str x5, [x4]
add x_om, x_om, #1
ubfiz x3, x1, #59, #5
orr x3, x3, x_om
add x4, x_pos, x1, lsl #3
ldr x5, [x4]
str x3, [x5], #8
str x5, [x4]
add x_om, x_om, #1
subs x_cnt, x_cnt, #8
b.ne .Lscw1
```
```asm
; kernel: rev_b2_build_hist_u8_stadd items=8 charged=6.90 key=c_bld_hist
.Lbh8:
subs x1, x_k, x_pms
add x2, x_k, x_d
csel x_k, x1, x2, hs
bic x4, x_k, #7
add x4, x_hist, x4
stadd x_eight, [x4]
subs x1, x_k, x_pms
add x2, x_k, x_d
csel x3, x1, x2, hs
bic x4, x3, #7
add x4, x_hist, x4
stadd x_eight, [x4]
stp x_k, x3, [x_kbuf], #16
subs x1, x3, x_pms
add x2, x3, x_d
csel x_k, x1, x2, hs
bic x4, x_k, #7
add x4, x_hist, x4
stadd x_eight, [x4]
subs x1, x_k, x_pms
add x2, x_k, x_d
csel x3, x1, x2, hs
bic x4, x3, #7
add x4, x_hist, x4
stadd x_eight, [x4]
stp x_k, x3, [x_kbuf], #16
subs x1, x3, x_pms
add x2, x3, x_d
csel x_k, x1, x2, hs
bic x4, x_k, #7
add x4, x_hist, x4
stadd x_eight, [x4]
subs x1, x_k, x_pms
add x2, x_k, x_d
csel x3, x1, x2, hs
bic x4, x3, #7
add x4, x_hist, x4
stadd x_eight, [x4]
stp x_k, x3, [x_kbuf], #16
subs x1, x3, x_pms
add x2, x3, x_d
csel x_k, x1, x2, hs
bic x4, x_k, #7
add x4, x_hist, x4
stadd x_eight, [x4]
subs x1, x_k, x_pms
add x2, x_k, x_d
csel x3, x1, x2, hs
bic x4, x3, #7
add x4, x_hist, x4
stadd x_eight, [x4]
stp x_k, x3, [x_kbuf], #16
mov x_k, x3
subs x_cnt, x_cnt, #8
b.hs .Lbh8
```
```asm
; kernel: rev_build_prefix_u64_8x items=8 charged=2.25 key=c_lbkt_pre
.Lpre8:
ldp x0, x1, [x_pos]
add x2, x_cur, x0
add x3, x2, x1
stp x_cur, x2, [x_pos], #16
ldp x4, x5, [x_pos]
add x6, x3, x4
add x7, x6, x5
stp x3, x6, [x_pos], #16
ldp x0, x1, [x_pos]
add x2, x7, x0
add x3, x2, x1
stp x7, x2, [x_pos], #16
ldp x4, x5, [x_pos]
add x6, x3, x4
add x_cur, x6, x5
stp x3, x6, [x_pos], #16
subs x_cnt, x_cnt, #8
b.ne .Lpre8
```
```asm
; kernel: rev_b2_build_scatter_r1_u8 items=8 charged=8.75 key=c_bld_scat
.Lsc8:
ldp x0, x1, [x_kbuf], #16
lsr x2, x0, #3
ubfiz x3, x0, #59, #5
orr x3, x3, x_om
add x4, x_pos, x2, lsl #3
ldr x5, [x4]
str x3, [x5], #8
str x5, [x4]
add x_om, x_om, #1
lsr x2, x1, #3
ubfiz x3, x1, #59, #5
orr x3, x3, x_om
add x4, x_pos, x2, lsl #3
ldr x5, [x4]
str x3, [x5], #8
str x5, [x4]
add x_om, x_om, #1
ldp x0, x1, [x_kbuf], #16
lsr x2, x0, #3
ubfiz x3, x0, #59, #5
orr x3, x3, x_om
add x4, x_pos, x2, lsl #3
ldr x5, [x4]
str x3, [x5], #8
str x5, [x4]
add x_om, x_om, #1
lsr x2, x1, #3
ubfiz x3, x1, #59, #5
orr x3, x3, x_om
add x4, x_pos, x2, lsl #3
ldr x5, [x4]
str x3, [x5], #8
str x5, [x4]
add x_om, x_om, #1
ldp x0, x1, [x_kbuf], #16
lsr x2, x0, #3
ubfiz x3, x0, #59, #5
orr x3, x3, x_om
add x4, x_pos, x2, lsl #3
ldr x5, [x4]
str x3, [x5], #8
str x5, [x4]
add x_om, x_om, #1
lsr x2, x1, #3
ubfiz x3, x1, #59, #5
orr x3, x3, x_om
add x4, x_pos, x2, lsl #3
ldr x5, [x4]
str x3, [x5], #8
str x5, [x4]
add x_om, x_om, #1
ldp x0, x1, [x_kbuf], #16
lsr x2, x0, #3
ubfiz x3, x0, #59, #5
orr x3, x3, x_om
add x4, x_pos, x2, lsl #3
ldr x5, [x4]
str x3, [x5], #8
str x5, [x4]
add x_om, x_om, #1
lsr x2, x1, #3
ubfiz x3, x1, #59, #5
orr x3, x3, x_om
add x4, x_pos, x2, lsl #3
ldr x5, [x4]
str x3, [x5], #8
str x5, [x4]
add x_om, x_om, #1
subs x_cnt, x_cnt, #8
b.ne .Lsc8
```
```asm
; kernel: rev_b2_memset_ff_u16 items=16 charged=0.65 key=c_memset_dir
.Lmff16:
stp x_ff, x_ff, [x_slot], #16
stp x_ff, x_ff, [x_slot], #16
stp x_ff, x_ff, [x_slot], #16
stp x_ff, x_ff, [x_slot], #16
stp x_ff, x_ff, [x_slot], #16
stp x_ff, x_ff, [x_slot], #16
stp x_ff, x_ff, [x_slot], #16
stp x_ff, x_ff, [x_slot], #16
subs x_cnt, x_cnt, #128
b.ne .Lmff16
```
```asm
; kernel: rev_b2_build_direct1_u8 items=8 charged=11.25 key=c_bld_dir1
.Lbdir8:
subs x8, x_k, x_pms
add x9, x_k, x_d
csel x_k, x8, x9, hs
lsr x10, x_k, #5
ubfiz x11, x_k, #59, #5
orr x11, x11, x_om
ldr x12, [x_slot, x10, lsl #3]
cmn x12, #1
b.ne .Lbdir_ovf0
str x11, [x_slot, x10, lsl #3]
.Lbdir_ret0:
add x_om, x_om, #1
subs x8, x_k, x_pms
add x9, x_k, x_d
csel x_k, x8, x9, hs
lsr x10, x_k, #5
ubfiz x11, x_k, #59, #5
orr x11, x11, x_om
ldr x12, [x_slot, x10, lsl #3]
cmn x12, #1
b.ne .Lbdir_ovf1
str x11, [x_slot, x10, lsl #3]
.Lbdir_ret1:
add x_om, x_om, #1
subs x8, x_k, x_pms
add x9, x_k, x_d
csel x_k, x8, x9, hs
lsr x10, x_k, #5
ubfiz x11, x_k, #59, #5
orr x11, x11, x_om
ldr x12, [x_slot, x10, lsl #3]
cmn x12, #1
b.ne .Lbdir_ovf2
str x11, [x_slot, x10, lsl #3]
.Lbdir_ret2:
add x_om, x_om, #1
subs x8, x_k, x_pms
add x9, x_k, x_d
csel x_k, x8, x9, hs
lsr x10, x_k, #5
ubfiz x11, x_k, #59, #5
orr x11, x11, x_om
ldr x12, [x_slot, x10, lsl #3]
cmn x12, #1
b.ne .Lbdir_ovf3
str x11, [x_slot, x10, lsl #3]
.Lbdir_ret3:
add x_om, x_om, #1
subs x8, x_k, x_pms
add x9, x_k, x_d
csel x_k, x8, x9, hs
lsr x10, x_k, #5
ubfiz x11, x_k, #59, #5
orr x11, x11, x_om
ldr x12, [x_slot, x10, lsl #3]
cmn x12, #1
b.ne .Lbdir_ovf4
str x11, [x_slot, x10, lsl #3]
.Lbdir_ret4:
add x_om, x_om, #1
subs x8, x_k, x_pms
add x9, x_k, x_d
csel x_k, x8, x9, hs
lsr x10, x_k, #5
ubfiz x11, x_k, #59, #5
orr x11, x11, x_om
ldr x12, [x_slot, x10, lsl #3]
cmn x12, #1
b.ne .Lbdir_ovf5
str x11, [x_slot, x10, lsl #3]
.Lbdir_ret5:
add x_om, x_om, #1
subs x8, x_k, x_pms
add x9, x_k, x_d
csel x_k, x8, x9, hs
lsr x10, x_k, #5
ubfiz x11, x_k, #59, #5
orr x11, x11, x_om
ldr x12, [x_slot, x10, lsl #3]
cmn x12, #1
b.ne .Lbdir_ovf6
str x11, [x_slot, x10, lsl #3]
.Lbdir_ret6:
add x_om, x_om, #1
subs x8, x_k, x_pms
add x9, x_k, x_d
csel x_k, x8, x9, hs
lsr x10, x_k, #5
ubfiz x11, x_k, #59, #5
orr x11, x11, x_om
ldr x12, [x_slot, x10, lsl #3]
cmn x12, #1
b.ne .Lbdir_ovf7
str x11, [x_slot, x10, lsl #3]
.Lbdir_ret7:
add x_om, x_om, #1
subs x_cnt, x_cnt, #8
b.hs .Lbdir8
```
```asm
; kernel: rev_b2_build_direct1_ovf items=8 charged=4.25 key=c_bld_dir1_ovf
b .Lbdir_done
.Lbdir_ovf0:
stp x11, x12, [x_ovf], #16
str x_ovf_tag, [x_slot, x10, lsl #3]
add x_ovf_tag, x_ovf_tag, #1
b .Lbdir_ret0
.Lbdir_ovf1:
stp x11, x12, [x_ovf], #16
str x_ovf_tag, [x_slot, x10, lsl #3]
add x_ovf_tag, x_ovf_tag, #1
b .Lbdir_ret1
.Lbdir_ovf2:
stp x11, x12, [x_ovf], #16
str x_ovf_tag, [x_slot, x10, lsl #3]
add x_ovf_tag, x_ovf_tag, #1
b .Lbdir_ret2
.Lbdir_ovf3:
stp x11, x12, [x_ovf], #16
str x_ovf_tag, [x_slot, x10, lsl #3]
add x_ovf_tag, x_ovf_tag, #1
b .Lbdir_ret3
.Lbdir_ovf4:
stp x11, x12, [x_ovf], #16
str x_ovf_tag, [x_slot, x10, lsl #3]
add x_ovf_tag, x_ovf_tag, #1
b .Lbdir_ret4
.Lbdir_ovf5:
stp x11, x12, [x_ovf], #16
str x_ovf_tag, [x_slot, x10, lsl #3]
add x_ovf_tag, x_ovf_tag, #1
b .Lbdir_ret5
.Lbdir_ovf6:
stp x11, x12, [x_ovf], #16
str x_ovf_tag, [x_slot, x10, lsl #3]
add x_ovf_tag, x_ovf_tag, #1
b .Lbdir_ret6
.Lbdir_ovf7:
stp x11, x12, [x_ovf], #16
str x_ovf_tag, [x_slot, x10, lsl #3]
add x_ovf_tag, x_ovf_tag, #1
b .Lbdir_ret7
.Lbdir_done:
ret
```


*(Continuation of `KERNELS.md`; all ASM blocks are parsed and verified together by `audit_cm.py`.)*

## 2. Exact 5-Limb Baseline Scalar ARMv8-A CIOS Montgomery Multiplication (`MONT5 = 280.0`, `279` Opcodes, AAPCS64-Compliant)
One complete carry-exact, shift-included two-chain CIOS outer round executes **`48`** baseline scalar ARMv8-A opcodes (`5 x 48 = 240` across 5 unrolled rounds) + **`39`** AAPCS64 callee-saved `x19..x27` prologue/epilogue (`sub sp, sp, #80` / `add sp, sp, #80` explicit frame allocation; exact fixed argument register assignment **`x_a = x5, x_b = x6, x_n0 = x7, x_out = x8, x_p = x9`**, so `x_a = x5` and `x_p = x9` are consumed before the 5 rounds reuse `x5` and `x9` as scratch while `x_b = x6, x_n0 = x7, x_out = x8` remain untouched!) + load/zero-init/final-sub/store opcodes = **`279 <= 280.0 = MONT5`** Ir (inlined macro without `ret` uses `278` ops, or `279` ops as an internal leaf routine; bit-exactness of this exact 279-opcode register machine is verified LIVE on 32 random 320-bit `(a, b, p)` triples inside `audit_cm.py`):
```asm
; kernel: mont5_cios_round items=1 charged=48.0 key=mont5_round
ldr x5, [x_b], #8
mul x10, x0, x5
mul x11, x1, x5
mul x12, x2, x5
mul x13, x3, x5
mul x14, x4, x5
adds x15, x15, x10
adcs x16, x16, x11
adcs x17, x17, x12
adcs x27, x27, x13
adcs x19, x19, x14
adcs x20, x20, xzr
adc x21, xzr, xzr
umulh x10, x0, x5
umulh x11, x1, x5
umulh x12, x2, x5
umulh x13, x3, x5
umulh x14, x4, x5
adds x16, x16, x10
adcs x17, x17, x11
adcs x27, x27, x12
adcs x19, x19, x13
adcs x20, x20, x14
adc x21, x21, xzr
mul x9, x15, x_n0
mul x10, x9, x22
mul x11, x9, x23
mul x12, x9, x24
mul x13, x9, x25
mul x14, x9, x26
adds xzr, x15, x10
adcs x15, x16, x11
adcs x16, x17, x12
adcs x17, x27, x13
adcs x27, x19, x14
adcs x19, x20, xzr
adc x20, x21, xzr
umulh x10, x9, x22
umulh x11, x9, x23
umulh x12, x9, x24
umulh x13, x9, x25
umulh x14, x9, x26
adds x15, x15, x10
adcs x16, x16, x11
adcs x17, x17, x12
adcs x27, x27, x13
adcs x19, x19, x14
adc x20, x20, xzr
```
```asm
; kernel: mont5_cios_epilogue items=1 charged=39.0 key=mont5_epi
sub sp, sp, #80
stp x19, x20, [sp, #0]
stp x21, x22, [sp, #16]
stp x23, x24, [sp, #32]
stp x25, x26, [sp, #48]
str x27, [sp, #64]
ldp x0, x1, [x_a, #0]
ldp x2, x3, [x_a, #16]
ldr x4, [x_a, #32]
ldp x22, x23, [x_p, #0]
ldp x24, x25, [x_p, #16]
ldr x26, [x_p, #32]
mov x15, #0
mov x16, #0
mov x17, #0
mov x27, #0
mov x19, #0
mov x20, #0
subs x10, x15, x22
sbcs x11, x16, x23
sbcs x12, x17, x24
sbcs x13, x27, x25
sbcs x14, x19, x26
sbcs xzr, x20, xzr
csel x15, x10, x15, hs
csel x16, x11, x16, hs
csel x17, x12, x17, hs
csel x27, x13, x27, hs
csel x19, x14, x19, hs
stp x15, x16, [x_out, #0]
stp x17, x27, [x_out, #16]
str x19, [x_out, #32]
ldp x19, x20, [sp, #0]
ldp x21, x22, [sp, #16]
ldp x23, x24, [sp, #32]
ldp x25, x26, [sp, #48]
ldr x27, [sp, #64]
add sp, sp, #80
ret
```
Total exact AAPCS64-compliant `MONT5 = 5 * 48 + 39 = 279 <= 280.0` Ir.

## 2b. Exact 2-Limb Baseline Scalar ARMv8.0-A CIOS Montgomery Multiplication (`mont2_cios`, `47 <= 48` Opcodes <= `58.84` Ir Charged for `L = 2`)
For small composite ECM branches (`<= 128` bits, `L_min = 2` in `q15_ecm16.py`, eliminating the `5`-limb floor inflation), `mulmod_instr(nbits, l_min=2)` charges `MONT5_desc * (2/5)^2 = 367.76 * 4 / 25 = 58.8416` Ir (even under `MONT5 = 280.0`, `L = 2` CIOS executes in **`47 <= 48.0 < 58.84` Ir**!). Inputs `(a0, a1, b0, b1, n0, n1, n0inv)` reside in AAPCS64 caller-saved registers `x0..x6` (zero stack spills, zero callee-saved registers, zero banned ops, disjoint GPRs `x0..x15`), returning `a * b * 2^{-128} mod N` in `x0:x1` (`0 <= t < 2N < 2^129` after both rounds for every 128-bit odd modulus `N <= 2^128 - 1`, verified bit-exactly LIVE in `cm_models.verify_kernels_symbolic` across edge moduli `2^128 - 1, 2^127 + 1, 3` and `64` random 128-bit odd moduli):
```asm
; kernel: mont2_cios items=1 charged=48.0 key=mont2_cios
mul x7, x0, x2
umulh x8, x0, x2
mul x9, x1, x2
umulh x10, x1, x2
adds x8, x8, x9
adc x9, x10, xzr
mul x11, x7, x6
mul x12, x11, x4
umulh x13, x11, x4
mul x14, x11, x5
umulh x15, x11, x5
adds xzr, x7, x12
adcs x7, x8, x14
adcs x8, x9, xzr
adc x9, xzr, xzr
adds x7, x7, x13
adcs x8, x8, x15
adc x9, x9, xzr
mul x10, x0, x3
umulh x11, x0, x3
mul x12, x1, x3
umulh x13, x1, x3
adds x7, x7, x10
adcs x8, x8, x12
adcs x9, x9, xzr
adc x10, xzr, xzr
adds x8, x8, x11
adcs x9, x9, x13
adc x10, x10, xzr
mul x11, x7, x6
mul x12, x11, x4
umulh x13, x11, x4
mul x14, x11, x5
umulh x15, x11, x5
adds xzr, x7, x12
adcs x7, x8, x14
adcs x8, x9, xzr
adc x9, x10, xzr
adds x7, x7, x13
adcs x8, x8, x15
adc x9, x9, xzr
subs x10, x7, x4
sbcs x11, x8, x5
sbcs xzr, x9, xzr
csel x0, x10, x7, hs
csel x1, x11, x8, hs
ret
```
## 3. Explicit Multi-Limb Cyclotomic Tower Norm Kernels, Parsed 4-Prime 62-Bit NTT (`fft_F`) & Garner CRT, Stage B, 2D Cofactorization, & Deterministic Stepped-GCD `H-ECM`
- **Explicit Schoolbook Multi-Limb Scalar ARMv8-A Tower Norm Derivation (`Q(zeta_15) -> Q(omega, sqrt(5)) -> Q(omega) -> Q`, ZERO Karatsuba shortcut):**
  Explicit schoolbook scalar ARMv8-A `L_1 x L_2` limb multiplication uses `mul_LL(L_1, L_2) = 8 L_1 L_2 + 4 (L_1 + L_2)` (`mul/umulh/adds/adcs` carry chains). On the production sieving region `Q_CM(a) <= s_c * 310,721.91 <= 565,514` (`s_c in {1.0, 1.40, 1.82}`), `b in K_{3/4}(982.55)` (`max_m |sigma_m(b)| <= 982.55`), the exact Vandermonde inverse `C_1 = 0.652889` (verified live in `unitdomain.py`) proves `|a_j| <= 388.33`, `|a^3_j| < 2^{27.4} < 2^{30}` and `|b^3_j| < 2^{29.21} < 2^{31}` (1 signed limb; `4` `8x8` 1-limb poly muls + `Phi_15` reduction = `400` Ir), `|g(a, b)_j| = |a^3 - u w^2 b^3|_j < 2^{92.92} << 2^{127}` (2 signed limbs; `+96` Ir), and `|f(a, b)_j| < 2^{89.10} << 2^{127}` (fits in 2 signed limbs, conservatively priced as 3 limbs: `a^6, b^6` = `2 x 928 = 1,856` Ir, `a^9, b^9, a^6 b^3, a^3 b^6` = `4 x 1,152 = 4,608` Ir, combine = `128` Ir). Using pure schoolbook `mul_LL` across the 3-step tower norm `Q(zeta_15) -> Q(omega, sqrt(5)) -> Q(omega) -> Q`:
  1. **`_nk_g_only` (`2 -> 4 -> 8 -> 16` limbs):** `32` `2x2` (`1,536`) + `6` `4x4` (`960`) + `3` `8x8` (`1,728`) + adds/subs/setup/cold-scan (`1,544`) = **`5,768 <= 6,500` Ir charged**!
  2. **`_nk_f_only` (`3 -> 6 -> 12 -> 24` limbs):** reuses `a^3, b^3` -> `a^6, b^6, a^9, b^9, a^6 b^3, a^3 b^6` + combine (`6,592`) + 3-limb schoolbook tower norm (`4,608 + 4,032 + 3,744 = 12,384`) = **`18,976 <= 20,000` Ir charged**!
- **Parsed 62-Bit Shoup Radix-2 NTT Butterfly (`ntt62_butterfly`: `18` opcodes/2-elems/prime -> `3 x 36.0 = 108.0 <= 168.0` Ir/elem/stage for 2 fwd + 1 inv across 4 primes) & 4-Prime Garner CRT (`fft_F`):**
  `p_i = c_i * 2^53 + 1 < 2^62` for `c_i in {269, 276, 305, 309}` (`log2(NTT4_PROD) = 244.70`, all 4 verified prime LIVE). One exact baseline scalar ARMv8-A radix-2 Shoup butterfly (`w_s = floor(w * 2^64 / p)` precomputed, `t = w*v - umulh(v, w_s)*p` reduced to `[0, p)`) takes **`18` opcodes per 2 elements (`9.0` Ir/elem/prime -> `36.0` Ir/elem/transform across 4 primes -> `3 x 36.0 = 108.0 <= 168.0` Ir/elem/stage for 2 forward + 1 inverse NTTs in `fft_F`, `1.56x` headroom!)**; scalar `F_{p_i}` Barrett MAC in `PM-Basis` uses `mul/umulh` + 1-step Barrett `umulh/msub` + 2 conditional subs + add/sub/csel + 2 loads + loop = **`17 <= 24` opcodes** (verified bit-exactly LIVE on all 4 NTT primes in `kernels_symbolic`); 4-prime Garner CRT uses 3 Shoup mod-muls (`15` ops) + 3 multi-limb mul-adds (`<= 75` ops) = `<= 90 << 220` linear term in `fft_F` (`168 M log2(M) + 220 M + 1200 + (504 M if N > 2^26 else 0)`; 2D `168 M_tot log2(M_tot) + 1100 M_tot + 1200`). Both `ntt62_butterfly` (64 random triples) and 4-prime Garner CRT (32 random 240-bit integers) are verified bit-exactly LIVE in `audit_cm.py`, alongside `srt_toy` (128 planted leaves, `0` mismatches):
```asm
; kernel: ntt62_butterfly items=2 charged=10.0 key=ntt62_elem_prime
.Lntt:
ldr x0, [x_u]
ldr x1, [x_v]
ldp x2, x3, [x_w], #16
mul x4, x1, x2
umulh x5, x1, x3
msub x4, x5, x_p, x4
subs x6, x4, x_p
csel x4, x6, x4, hs
add x7, x0, x4
subs x8, x7, x_p
csel x7, x8, x7, hs
subs x9, x0, x4
add x10, x9, x_p
csel x9, x9, x10, hs
str x7, [x_u], #8
str x9, [x_v], #8
subs x_cnt, x_cnt, #1
b.ne .Lntt
```
- **Parsed 3-Prime `< 2^61` Register-Shared Lazy-8p Shoup Radix-4 NTT (`ntt61_r4_u2`: `66` opcodes / `16` elem-stages = `4.125 <= 4.13` Ir/elem/stage -> `9 * 4.125 = 37.125 <= 37.50` in `ntt3_F`), Schoolbook Small-`M` Crossover (`<= 256` limbs), `100*n` Garner CRT, & 64-Bit Lehmer GCD (`lehmer32_step`: `29/2 = 14.5 <= 16.0`):**
  Three `< 2^61` primes divisible by `3 * 2^52`: `p_1 = 456 * 2^52 + 1`, `p_2 = 447 * 2^52 + 1`, `p_3 = 432 * 2^52 + 1` (`8*p_1 = 0.8906 * 2^64 < 2^64`, `prod = 2^182.392 > 2^54 * (2^64 - 1)^2`). **In-Place Harvey Bit-Reversed-Twiddle Radix-4 Schedule (Zero Bit-Reversal Permutation, Verified Live in `cm_models.py::_r4_r3_ntt`):** Forward DIT runs stages `m = N/4, N/16, ..., 1` (standard -> bit-reversed order) and inverse 1-clamp DIF runs stages `m = 1, 4, ..., N/4` (bit-reversed -> standard order), with pointwise multiplication executed directly in bit-reversed order (ZERO permutation pass, ZERO transpose!). Within every sub-block `k = 0 .. N/(4m) - 1` of stage `m`, ALL `m` **consecutive** words `a[4km + j]`, `a[4km + m + j]`, `a[4km + 2m + j]`, `a[4km + 3m + j]` (`j = 0, 2, ..., m - 2`) share **ONE single constant twiddle triple `(w1, w2, w3)`** (`w1 = omega^{bitrev(k)}`), so `ntt61_r4_u2` loads consecutive pairs via unit-stride `ldp [x_a0..x_a3]` and stores via `stp [x_a0..x_a3], #16` with zero twiddle reloads. Keeping DIT inputs in `[0, 8p)` (and DIF inputs in `[0, 2p)` with 1 clamp on `y0` only) ensures every Shoup product `a_i * w - floor(a_i * s / 2^64) * p` lands in `[0, 2p)` with ZERO conditional clamp; holding `(w1, s1, w2, s2, w3, s3)` in registers and forming `2p`/`4p` via `x_p, lsl #1`/`#2` fits strictly in `28 <= 28` disjoint GPRs (`66` ops / `16` elem-stages = `4.125` Ir/elem-stage). At the single stage `m = 1` (`u1`, `<= 5.25` Ir/elem-stage incl. twiddle load), the top stage `m = N/4` has trivial `w1 = w2 = w3 = 1` (`0` twiddle muls -> `3.00` Ir/elem-stage), so `(3.00 + 5.25)/2 = 4.125` exactly for ALL `N = 4^S >= 16` (`9 * 4.125 = 37.125 <= 37.50` in `ntt3_F`)!
```asm
; kernel: ntt61_r4_u2 items=16 charged=4.13 key=ntt61_r4_elem_stage
.Lntt_r4_u2:
ldp x0, x14, [x_a0]
ldp x1, x15, [x_a1]
ldp x2, x16, [x_a2]
ldp x3, x17, [x_a3]
mul x6, x1, x_w1
umulh x7, x1, x_s1
msub x1, x7, x_p, x6
mul x6, x2, x_w2
umulh x7, x2, x_s2
msub x2, x7, x_p, x6
mul x6, x3, x_w3
umulh x7, x3, x_s3
msub x3, x7, x_p, x6
mul x6, x15, x_w1
umulh x7, x15, x_s1
msub x15, x7, x_p, x6
mul x6, x16, x_w2
umulh x7, x16, x_s2
msub x16, x7, x_p, x6
mul x6, x17, x_w3
umulh x7, x17, x_s3
msub x17, x7, x_p, x6
subs x8, x0, x_p, lsl #2
csel x0, x8, x0, hs
subs x8, x0, x_p, lsl #1
csel x0, x8, x0, hs
add x8, x0, x2
add x10, x0, x_p, lsl #1
sub x10, x10, x2
add x11, x1, x3
add x12, x1, x_p, lsl #1
sub x12, x12, x3
mul x6, x12, x_j
umulh x7, x12, x_jp
msub x12, x7, x_p, x6
add x0, x8, x11
add x2, x8, x_p, lsl #2
sub x2, x2, x11
add x1, x10, x12
add x3, x10, x_p, lsl #1
sub x3, x3, x12
subs x8, x14, x_p, lsl #2
csel x14, x8, x14, hs
subs x8, x14, x_p, lsl #1
csel x14, x8, x14, hs
add x8, x14, x16
add x10, x14, x_p, lsl #1
sub x10, x10, x16
add x11, x15, x17
add x12, x15, x_p, lsl #1
sub x12, x12, x17
mul x6, x12, x_j
umulh x7, x12, x_jp
msub x12, x7, x_p, x6
add x14, x8, x11
add x16, x8, x_p, lsl #2
sub x16, x16, x11
add x15, x10, x12
add x17, x10, x_p, lsl #1
sub x17, x17, x12
stp x0, x14, [x_a0], #16
stp x1, x15, [x_a1], #16
stp x2, x16, [x_a2], #16
stp x3, x17, [x_a3], #16
subs x_cnt, x_cnt, #2
b.ne .Lntt_r4_u2
```
```asm
; kernel: lehmer32_step items=2 charged=16.0 key=lehmer32_step_op
.Lleh_loop:
cbz x1, .Lleh_done
udiv x9, x0, x1
msub x10, x9, x1, x0
madd x11, x9, x4, x2
madd x12, x9, x5, x3
cmp x10, x12
b.lo .Lleh_done
sub x8, x1, x10
add x13, x4, x11
cmp x8, x13
b.lo .Lleh_done
cbz x10, .Lleh_done
udiv x9, x1, x10
msub x14, x9, x10, x1
madd x15, x9, x11, x4
madd x16, x9, x12, x5
cmp x14, x15
b.lo .Lleh_done
sub x8, x10, x14
add x13, x12, x16
cmp x8, x13
b.lo .Lleh_done
mov x0, x10
mov x1, x14
mov x2, x11
mov x3, x12
mov x4, x15
mov x5, x16
b .Lleh_loop
.Lleh_done:
```
- **Stage B (`P_min = 33`, `KH = (0.060, 0.075)`, Per-Column Odd-Norm `B_{g,odd,max}(b, c)`-Adaptive `T_g_lo(b, c) in [30, 115]` + Column-Constant Cold-Path 2D NP Gate `T_jw(b, c) = T_g_lo(b, c) + 55`, Single-Tier `[72.0]` `2D-56b NTT3` `f`-Batch & `g`-Batch, `step = 0.25`):** Reusing the exact raw-opcode-verified `ntt61_r4_u2` kernel (`< 2^{61}`, `3 * 2^{52}` roots, `NTT3_PROD = 2^{182.39}`, ZERO new ASM!) via 56-bit-limb 2D bivariate acyclic NTT3 (`canary1_verifiers.py::gbatch_ntt3_2d`, `MONT5 = 280`, schoolbook `M <= 256`, `lehmer32_step` 64b Lehmer GCD) yields **`c_cof_2d_upper_live = 0.66798` Ir/cell** (`fbatch_tiers = [72.0]`, Chebotarev `5/72 = 0.0694 < 0.075`, `f`-side `2/48 = 0.0417 < 0.060`, live `2^14` `0.07197/0.05649`). The same verified `ntt61_r4_u2` kernel also executes Kronecker-NTT3 + BLS Middle-Product Stage-2 ECM (`q15_ecm_ntt3.py`, `RAM <= 1.463 GiB << 16 EiB`, ZERO new ASM!).
- **Deterministic Stepped-GCD + Kronecker-NTT3 BLS Stage-2 `H-ECM` Ladder (`q15_ecm16.py`, `q15_ecm_ntt3.py`, `cm_models.ecm_mean_var`):** Exact Newton perfect-power precheck (`2^26` Ir > 32.06M ops for 2,590b boot) removes all `N = p^k`. Stage 1 steps per prime factor `q` across all `q^k <= B1` (`<= pi_pp_upper(B1)` gcds); EVERY rung credits **strictly Stage-2 hits ONLY** (`ord_p(Q_1) = q_p in (B1, B2]`, ZERO Stage-1 credit -> Stage-1 simultaneous-hit loss is mathematically `0.0`!; in descent, static `278 MB` binary prime-product tree `T_{26}` first strips ALL primes `<= 2^{26}` deterministically with prob `1.0`, followed by `pre1 [26, 48]`, while collection needs ZERO pre-rungs). BOTH `bsgs` and `ntt3_bls` enforce `W <= B1` and `M_max := G * W + W/2 < B1^2` (`G = ceil(B2/W) + 1`), charge block-tree bisection (`<= 2*h_blk*omega_max` GCDs) PLUS an unconditional per-curve 3-step zero-leaf fallback across all `max(1, omega_max // 2)` composite baby buckets (live-verified on `N = 35`, `50023*50033`, `1061*1409`, and `643*757` at `M = 301 = 7*43`), multiply `p1` by `(1 - 1/B1)` for `Q_0` projection, and discount Subcase B2b shared-`q_p` collisions by the exact Hasse-interval bound `(1 - min(1.0, (omega_max - 1) * (1/B1 + 4/2^{0.5*r_lo})))` (ZERO curve replay, `var_t = 0.0`). `tree_branch_caps(nbits, b_min, l_min=2)` uses exact `s_i = b_i - 1` DP (`mont2_cios` `47 <= 48` ops). Collection (`lo_r = 48, b_min = 96, success = 0.999999`, `n_rels_upper = 2^(log2_rels+0.55)`) has `var_t = 0.0` with `1 - 29e-6` yield factor in `det_cred_pen`; descent (`lo_bits = 16, b_min = 32, success = 0.999`) is 100% deterministic.

## 4. `Z[zeta_12]` Block Wiedemann Kernels (1-Dense + 16-Shard `< 2^64 B` 4-Vec 4-Col Arch B + `u32->u4->bin21` Streaming + Single-Sub `57`-Op Barrett `row_red5_z12`)
Over `Z/rZ` (`r` 255b -> 4 limbs, `phi(12) = 4`, `tab_ratio = 1.00` ZERO-`TABM`, `32` 4-vec bundles `128 B/slot`, `n = 128, m = 32,768`, 1-owner dense RAM `< 16 EiB`, 16 sparse shards `<< 16 EiB`, explicit cross-owner streaming `c_stream_mat = 63.0 / 251.0`):
- **`unit_run_4vec4_z12` & `_neg` (`122/16 = 7.625 <= 7.65` Ir/vec-entry, `28` GPRs, ZERO `x18/x29/x30`, `0` spills):** Pre-scaled byte offset `col * 128` added to `x6 = x_base` (ping-pong native, 0 copies); 4 columns x 4 vectors into 4 signed 5-limb accs (`x8..x28`), tails priced inside `c_row_chi = 600.0`:
```asm
; kernel: unit_run_4vec4_wrap items=1 charged=4.0 key=unit_run_4vec4_wrap
subs x5, x5, #4
b.lo .Lrun4_tail
.Lrun4_tail:
adds x5, x5, #4
b.eq .Lrun4_done
.Lrun4_done:
```
```asm
; kernel: unit_run_4vec4_z12 items=16 charged=7.65 key=unit_run_4vec4
.Lrun4v4:
ldr x7, [x4], #8
add x7, x6, x7
ldp x0, x1, [x7, #0]
ldp x2, x3, [x7, #16]
adds x8, x8, x0
adcs x9, x9, x1
adcs x10, x10, x2
adcs x11, x11, x3
adc x12, x12, xzr
ldp x0, x1, [x7, #32]
ldp x2, x3, [x7, #48]
adds x13, x13, x0
adcs x14, x14, x1
adcs x15, x15, x2
adcs x16, x16, x3
adc x17, x17, xzr
ldp x0, x1, [x7, #64]
ldp x2, x3, [x7, #80]
adds x19, x19, x0
adcs x20, x20, x1
adcs x21, x21, x2
adcs x22, x22, x3
adc x23, x23, xzr
ldp x0, x1, [x7, #96]
ldp x2, x3, [x7, #112]
adds x24, x24, x0
adcs x25, x25, x1
adcs x26, x26, x2
adcs x27, x27, x3
adc x28, x28, xzr
ldr x7, [x4], #8
add x7, x6, x7
ldp x0, x1, [x7, #0]
ldp x2, x3, [x7, #16]
adds x8, x8, x0
adcs x9, x9, x1
adcs x10, x10, x2
adcs x11, x11, x3
adc x12, x12, xzr
ldp x0, x1, [x7, #32]
ldp x2, x3, [x7, #48]
adds x13, x13, x0
adcs x14, x14, x1
adcs x15, x15, x2
adcs x16, x16, x3
adc x17, x17, xzr
ldp x0, x1, [x7, #64]
ldp x2, x3, [x7, #80]
adds x19, x19, x0
adcs x20, x20, x1
adcs x21, x21, x2
adcs x22, x22, x3
adc x23, x23, xzr
ldp x0, x1, [x7, #96]
ldp x2, x3, [x7, #112]
adds x24, x24, x0
adcs x25, x25, x1
adcs x26, x26, x2
adcs x27, x27, x3
adc x28, x28, xzr
ldr x7, [x4], #8
add x7, x6, x7
ldp x0, x1, [x7, #0]
ldp x2, x3, [x7, #16]
adds x8, x8, x0
adcs x9, x9, x1
adcs x10, x10, x2
adcs x11, x11, x3
adc x12, x12, xzr
ldp x0, x1, [x7, #32]
ldp x2, x3, [x7, #48]
adds x13, x13, x0
adcs x14, x14, x1
adcs x15, x15, x2
adcs x16, x16, x3
adc x17, x17, xzr
ldp x0, x1, [x7, #64]
ldp x2, x3, [x7, #80]
adds x19, x19, x0
adcs x20, x20, x1
adcs x21, x21, x2
adcs x22, x22, x3
adc x23, x23, xzr
ldp x0, x1, [x7, #96]
ldp x2, x3, [x7, #112]
adds x24, x24, x0
adcs x25, x25, x1
adcs x26, x26, x2
adcs x27, x27, x3
adc x28, x28, xzr
ldr x7, [x4], #8
add x7, x6, x7
ldp x0, x1, [x7, #0]
ldp x2, x3, [x7, #16]
adds x8, x8, x0
adcs x9, x9, x1
adcs x10, x10, x2
adcs x11, x11, x3
adc x12, x12, xzr
ldp x0, x1, [x7, #32]
ldp x2, x3, [x7, #48]
adds x13, x13, x0
adcs x14, x14, x1
adcs x15, x15, x2
adcs x16, x16, x3
adc x17, x17, xzr
ldp x0, x1, [x7, #64]
ldp x2, x3, [x7, #80]
adds x19, x19, x0
adcs x20, x20, x1
adcs x21, x21, x2
adcs x22, x22, x3
adc x23, x23, xzr
ldp x0, x1, [x7, #96]
ldp x2, x3, [x7, #112]
adds x24, x24, x0
adcs x25, x25, x1
adcs x26, x26, x2
adcs x27, x27, x3
adc x28, x28, xzr
subs x5, x5, #4
b.hs .Lrun4v4
```
```asm
; kernel: unit_run_4vec4_neg_z12 items=16 charged=7.65 key=unit_run_4vec4_neg
.Lrun4v4n:
ldr x7, [x4], #8
add x7, x6, x7
ldp x0, x1, [x7, #0]
ldp x2, x3, [x7, #16]
subs x8, x8, x0
sbcs x9, x9, x1
sbcs x10, x10, x2
sbcs x11, x11, x3
sbc x12, x12, xzr
ldp x0, x1, [x7, #32]
ldp x2, x3, [x7, #48]
subs x13, x13, x0
sbcs x14, x14, x1
sbcs x15, x15, x2
sbcs x16, x16, x3
sbc x17, x17, xzr
ldp x0, x1, [x7, #64]
ldp x2, x3, [x7, #80]
subs x19, x19, x0
sbcs x20, x20, x1
sbcs x21, x21, x2
sbcs x22, x22, x3
sbc x23, x23, xzr
ldp x0, x1, [x7, #96]
ldp x2, x3, [x7, #112]
subs x24, x24, x0
sbcs x25, x25, x1
sbcs x26, x26, x2
sbcs x27, x27, x3
sbc x28, x28, xzr
ldr x7, [x4], #8
add x7, x6, x7
ldp x0, x1, [x7, #0]
ldp x2, x3, [x7, #16]
subs x8, x8, x0
sbcs x9, x9, x1
sbcs x10, x10, x2
sbcs x11, x11, x3
sbc x12, x12, xzr
ldp x0, x1, [x7, #32]
ldp x2, x3, [x7, #48]
subs x13, x13, x0
sbcs x14, x14, x1
sbcs x15, x15, x2
sbcs x16, x16, x3
sbc x17, x17, xzr
ldp x0, x1, [x7, #64]
ldp x2, x3, [x7, #80]
subs x19, x19, x0
sbcs x20, x20, x1
sbcs x21, x21, x2
sbcs x22, x22, x3
sbc x23, x23, xzr
ldp x0, x1, [x7, #96]
ldp x2, x3, [x7, #112]
subs x24, x24, x0
sbcs x25, x25, x1
sbcs x26, x26, x2
sbcs x27, x27, x3
sbc x28, x28, xzr
ldr x7, [x4], #8
add x7, x6, x7
ldp x0, x1, [x7, #0]
ldp x2, x3, [x7, #16]
subs x8, x8, x0
sbcs x9, x9, x1
sbcs x10, x10, x2
sbcs x11, x11, x3
sbc x12, x12, xzr
ldp x0, x1, [x7, #32]
ldp x2, x3, [x7, #48]
subs x13, x13, x0
sbcs x14, x14, x1
sbcs x15, x15, x2
sbcs x16, x16, x3
sbc x17, x17, xzr
ldp x0, x1, [x7, #64]
ldp x2, x3, [x7, #80]
subs x19, x19, x0
sbcs x20, x20, x1
sbcs x21, x21, x2
sbcs x22, x22, x3
sbc x23, x23, xzr
ldp x0, x1, [x7, #96]
ldp x2, x3, [x7, #112]
subs x24, x24, x0
sbcs x25, x25, x1
sbcs x26, x26, x2
sbcs x27, x27, x3
sbc x28, x28, xzr
ldr x7, [x4], #8
add x7, x6, x7
ldp x0, x1, [x7, #0]
ldp x2, x3, [x7, #16]
subs x8, x8, x0
sbcs x9, x9, x1
sbcs x10, x10, x2
sbcs x11, x11, x3
sbc x12, x12, xzr
ldp x0, x1, [x7, #32]
ldp x2, x3, [x7, #48]
subs x13, x13, x0
sbcs x14, x14, x1
sbcs x15, x15, x2
sbcs x16, x16, x3
sbc x17, x17, xzr
ldp x0, x1, [x7, #64]
ldp x2, x3, [x7, #80]
subs x19, x19, x0
sbcs x20, x20, x1
sbcs x21, x21, x2
sbcs x22, x22, x3
sbc x23, x23, xzr
ldp x0, x1, [x7, #96]
ldp x2, x3, [x7, #112]
subs x24, x24, x0
sbcs x25, x25, x1
sbcs x26, x26, x2
sbcs x27, x27, x3
sbc x28, x28, xzr
subs x5, x5, #4
b.hs .Lrun4v4n
```
- **`row_red5_z12` (`57 <= 58.0` Ir/bucket-vec — Exact Single-Sub Signed 5-Limb Barrett Reduction for ALL `|A| <= 2^{290}`):**
  Handles ANY signed 5-limb bucket accumulator `A = (x10..x14)` with `|A| <= 2^{290}` (covering EVERY `W_max = 4,096, |c| <= 255` row up to `4,096 * 255 * (R_BLS - 1) < 2^{275}` by `> 10^4x`!). Since `mu_r = floor(2^318 / R_BLS) = 0x8d54253b7fb78ddf` and `top = floor(|A| / 2^{254}) <= 2^{36}`, exact error `|A|/R_BLS - q_hat < 2^{254}/R_BLS + top/2^{64} + 1 < 1.5522 < 2` unconditionally -> `rem = U - q_hat * R_BLS in [0, 2 * R_BLS)` ALWAYS, so ONE single conditional `- R_BLS` + branchless sign fix reduces to exact `[0, R_BLS)` in **`57 <= 58.0`** opcodes (`6 * 57 = 342` ops/row):
```asm
; kernel: row_red5_z12 items=1 charged=58.0 key=row_red5_barrett
asr x15, x14, #63
eor x0, x10, x15
eor x1, x11, x15
eor x2, x12, x15
eor x3, x13, x15
eor x4, x14, x15
subs x0, x0, x15
sbcs x1, x1, x15
sbcs x2, x2, x15
sbcs x3, x3, x15
sbc x4, x4, x15
lsl x5, x4, #2
orr x5, x5, x3, lsr #62
umulh x5, x5, x_mu
mul x6, x5, x22
mul x7, x5, x23
mul x8, x5, x24
mul x9, x5, x25
umulh x16, x5, x22
umulh x17, x5, x23
umulh x20, x5, x24
umulh x19, x5, x25
adds x7, x7, x16
adcs x8, x8, x17
adcs x9, x9, x20
adc x19, x19, xzr
subs x0, x0, x6
sbcs x1, x1, x7
sbcs x2, x2, x8
sbcs x3, x3, x9
sbc x4, x4, x19
subs x6, x0, x22
sbcs x7, x1, x23
sbcs x8, x2, x24
sbcs x9, x3, x25
sbcs xzr, x4, xzr
csel x0, x6, x0, hs
csel x1, x7, x1, hs
csel x2, x8, x2, hs
csel x3, x9, x3, hs
orr x6, x0, x1
orr x7, x2, x3
orr x6, x6, x7
cmp x6, #0
csel x15, xzr, x15, eq
subs x6, x22, x0
sbcs x7, x23, x1
sbcs x8, x24, x2
sbc x9, x25, x3
tst x15, x15
csel x0, x6, x0, ne
csel x1, x7, x1, ne
csel x2, x8, x2, ne
csel x3, x9, x3, ne
stp x0, x1, [x_out, #0]
stp x2, x3, [x_out, #16]
ret
```
- **Sign-Split 2-Vec `gen_run1_2v_z12` (`c_b > 0`) & `gen_run1_2v_neg_z12` (`c_b < 0`) (`43/2 = 21.50 <= 22.0` Ir/vec-entry, `24 <= 29` GPRs avoiding `x18/x29/x30`, `0` spills; caller-side non-empty `cbz x24` precondition guarantees `x24 >= 1` on entry):** Each bucket splits non-unit scalar updates (`1 <= |c_b| <= 255`, provisioned up to `nonunit_ub = 0.110` covering every `2^20..2^26` ladder point and `1 - (eu_26 - 2*se) = 0.08625`, stored as 1-byte `uint8_t` magnitudes loaded and zero-extended to `x22` via `ldrb w22, [x23], #1`) into positive and negative sub-lists, executed as two spill-free 2-vector passes (`(v0, v1)` at offsets `#0..#48`, then `(v2, v3)` at offsets `#64..#112`, explicit `add x8, x9, x8` from `x_base = x9`): `3` ptr/coeff loads + `2 * 19` mul/umulh/carry ops + `2` loop = **`43 / 2 = 21.50 <= 22.0` Ir/vec-entry** (`k5_live = eu_lo * 7.65 + (1 - eu_lo) * 22.0 = 8.8877 <= k5_charged = 9.40` Ir, `+0.512` Ir headroom; even at `nonunit_ub = 0.110`, `0.890 * 7.65 + 0.110 * 22.0 = 9.2285 <= 9.40`!):
```asm
; kernel: gen_run1_2v_z12 items=2 charged=22.0 key=gen_run1_2v_coord
.Lgen2v:
ldr x8, [x21], #8
add x8, x9, x8
ldrb w22, [x23], #1
ldp x0, x1, [x8, #0]
ldp x2, x3, [x8, #16]
mul x4, x0, x22
mul x5, x1, x22
mul x6, x2, x22
mul x7, x3, x22
adds x10, x10, x4
adcs x11, x11, x5
adcs x12, x12, x6
adcs x13, x13, x7
adc x14, x14, xzr
umulh x4, x0, x22
umulh x5, x1, x22
umulh x6, x2, x22
umulh x7, x3, x22
adds x11, x11, x4
adcs x12, x12, x5
adcs x13, x13, x6
adc x14, x14, x7
ldp x0, x1, [x8, #32]
ldp x2, x3, [x8, #48]
mul x4, x0, x22
mul x5, x1, x22
mul x6, x2, x22
mul x7, x3, x22
adds x15, x15, x4
adcs x16, x16, x5
adcs x17, x17, x6
adcs x19, x19, x7
adc x20, x20, xzr
umulh x4, x0, x22
umulh x5, x1, x22
umulh x6, x2, x22
umulh x7, x3, x22
adds x16, x16, x4
adcs x17, x17, x5
adcs x19, x19, x6
adc x20, x20, x7
subs x24, x24, #1
b.gt .Lgen2v
```
```asm
; kernel: gen_run1_2v_neg_z12 items=2 charged=22.0 key=gen_run1_2v_neg_coord
.Lgen2vn:
ldr x8, [x21], #8
add x8, x9, x8
ldrb w22, [x23], #1
ldp x0, x1, [x8, #0]
ldp x2, x3, [x8, #16]
mul x4, x0, x22
mul x5, x1, x22
mul x6, x2, x22
mul x7, x3, x22
subs x10, x10, x4
sbcs x11, x11, x5
sbcs x12, x12, x6
sbcs x13, x13, x7
sbc x14, x14, xzr
umulh x4, x0, x22
umulh x5, x1, x22
umulh x6, x2, x22
umulh x7, x3, x22
subs x11, x11, x4
sbcs x12, x12, x5
sbcs x13, x13, x6
sbc x14, x14, x7
ldp x0, x1, [x8, #32]
ldp x2, x3, [x8, #48]
mul x4, x0, x22
mul x5, x1, x22
mul x6, x2, x22
mul x7, x3, x22
subs x15, x15, x4
sbcs x16, x16, x5
sbcs x17, x17, x6
sbcs x19, x19, x7
sbc x20, x20, xzr
umulh x4, x0, x22
umulh x5, x1, x22
umulh x6, x2, x22
umulh x7, x3, x22
subs x16, x16, x4
sbcs x17, x17, x5
sbcs x19, x19, x6
sbc x20, x20, x7
subs x24, x24, #1
b.gt .Lgen2vn
```
- **Explicit 1-Vector Horner Kernels (`unit_run_1vec4_z12` / `_neg`: `38/4 = 9.50 <= 9.60` Ir/entry; `gen_run1_1v_z12` / `_neg`: `24.0 <= 24.0` Ir/entry with caller-side non-empty `cbz x24` precondition `x24 >= 1` -> `k5_horner_charged = 11.20`):** Reuses the EXACT SAME read-only sparse-matrix shards (`128*j` stored byte offsets) on the `32 B/slot` 1-vector Horner dense buffer via **`add x7, x6, x7, lsr #2`** (`(128*j) >> 2 = 32*j` exact byte offset, ZERO matrix re-encoding, ZERO offset-shift bug!), simulated bit-exact in `cm_models.py`:
```asm
; kernel: unit_run_1vec4_z12 items=4 charged=9.60 key=unit_run_1vec4
.Lrun1v4:
ldr x7, [x4], #8
add x7, x6, x7, lsr #2
ldp x0, x1, [x7, #0]
ldp x2, x3, [x7, #16]
adds x8, x8, x0
adcs x9, x9, x1
adcs x10, x10, x2
adcs x11, x11, x3
adc x12, x12, xzr
ldr x7, [x4], #8
add x7, x6, x7, lsr #2
ldp x0, x1, [x7, #0]
ldp x2, x3, [x7, #16]
adds x8, x8, x0
adcs x9, x9, x1
adcs x10, x10, x2
adcs x11, x11, x3
adc x12, x12, xzr
ldr x7, [x4], #8
add x7, x6, x7, lsr #2
ldp x0, x1, [x7, #0]
ldp x2, x3, [x7, #16]
adds x8, x8, x0
adcs x9, x9, x1
adcs x10, x10, x2
adcs x11, x11, x3
adc x12, x12, xzr
ldr x7, [x4], #8
add x7, x6, x7, lsr #2
ldp x0, x1, [x7, #0]
ldp x2, x3, [x7, #16]
adds x8, x8, x0
adcs x9, x9, x1
adcs x10, x10, x2
adcs x11, x11, x3
adc x12, x12, xzr
subs x5, x5, #4
b.hs .Lrun1v4
```
```asm
; kernel: unit_run_1vec4_neg_z12 items=4 charged=9.60 key=unit_run_1vec4_neg
.Lrun1v4n:
ldr x7, [x4], #8
add x7, x6, x7, lsr #2
ldp x0, x1, [x7, #0]
ldp x2, x3, [x7, #16]
subs x8, x8, x0
sbcs x9, x9, x1
sbcs x10, x10, x2
sbcs x11, x11, x3
sbc x12, x12, xzr
ldr x7, [x4], #8
add x7, x6, x7, lsr #2
ldp x0, x1, [x7, #0]
ldp x2, x3, [x7, #16]
subs x8, x8, x0
sbcs x9, x9, x1
sbcs x10, x10, x2
sbcs x11, x11, x3
sbc x12, x12, xzr
ldr x7, [x4], #8
add x7, x6, x7, lsr #2
ldp x0, x1, [x7, #0]
ldp x2, x3, [x7, #16]
subs x8, x8, x0
sbcs x9, x9, x1
sbcs x10, x10, x2
sbcs x11, x11, x3
sbc x12, x12, xzr
ldr x7, [x4], #8
add x7, x6, x7, lsr #2
ldp x0, x1, [x7, #0]
ldp x2, x3, [x7, #16]
subs x8, x8, x0
sbcs x9, x9, x1
sbcs x10, x10, x2
sbcs x11, x11, x3
sbc x12, x12, xzr
subs x5, x5, #4
b.hs .Lrun1v4n
```
```asm
; kernel: gen_run1_1v_z12 items=1 charged=24.0 key=gen_run1_1v_coord
.Lgen1v:
ldr x8, [x21], #8
add x8, x9, x8, lsr #2
ldrb w22, [x23], #1
ldp x0, x1, [x8, #0]
ldp x2, x3, [x8, #16]
mul x4, x0, x22
mul x5, x1, x22
mul x6, x2, x22
mul x7, x3, x22
adds x10, x10, x4
adcs x11, x11, x5
adcs x12, x12, x6
adcs x13, x13, x7
adc x14, x14, xzr
umulh x4, x0, x22
umulh x5, x1, x22
umulh x6, x2, x22
umulh x7, x3, x22
adds x11, x11, x4
adcs x12, x12, x5
adcs x13, x13, x6
adc x14, x14, x7
subs x24, x24, #1
b.gt .Lgen1v
```
```asm
; kernel: gen_run1_1v_neg_z12 items=1 charged=24.0 key=gen_run1_1v_neg_coord
.Lgen1vn:
ldr x8, [x21], #8
add x8, x9, x8, lsr #2
ldrb w22, [x23], #1
ldp x0, x1, [x8, #0]
ldp x2, x3, [x8, #16]
mul x4, x0, x22
mul x5, x1, x22
mul x6, x2, x22
mul x7, x3, x22
subs x10, x10, x4
sbcs x11, x11, x5
sbcs x12, x12, x6
sbcs x13, x13, x7
sbc x14, x14, xzr
umulh x4, x0, x22
umulh x5, x1, x22
umulh x6, x2, x22
umulh x7, x3, x22
subs x11, x11, x4
sbcs x12, x12, x5
sbcs x13, x13, x6
sbc x14, x14, x7
subs x24, x24, #1
b.gt .Lgen1vn
```
- **`chi_mul_z12_fast` (`50` Ir), Non-Double-Counted `row_tail_6bucket_z12` (`240` Ir = `12` unit sub-lists x (`4/4` wrap + `<=3*0.35` tail delta above `7.65`) + `12*0.5` non-unit header + `2*50` `chi_mul` + `56` header/zero/store + `53` slack), 1-Dense + 16-Projection + 16-Shard 4-Vec Arch B Row Overhead (`c_row_chi = 600.0` Ir, `c_stream_mat = 63.0` Ir/vec-row Krylov, `c_stream_mat_1vec = 251.0` Ir/row 1-vec Horner via `u32->u4->bin21` `1.0625*212.639 + 24.0 = 249.93 <= 251.0`, `nonunit_ub = 0.110`, `k5_charged = 9.40`, `k5_horner_charged = 11.20`, `c_xproj = 65.0` Ir/vec-row `>= 64.5`, `n_owners = 1`, `n_proj_owners = 16`, `n_mat_owners = 16`, `te = 0`):** `6 * 57 + 240 + 16 = 598 <= 600.0 = c_row_chi` Ir worst-case!

## 5. Descent Kernels (Generous Scalar Bounds for Uncapped `dsieve16`, Exact-Integer `BOUND_CELL`, & Input-Uniform Bareiss LLL)
- **Uncapped `dsieve16` on `E_16(R_gram) cap C_30` with Flat `line_base[ell] + x_0` Addressing (`SIEVE_HIT = 1,400` Ir/hit, ZERO skip caps):** Each `E_16(R_gram) cap C_30` node is tiled into `<= 2^{36}`-cell Fincke-Pohst slabs (`<= 128 GiB` `uint16_t` accumulator buffer + `<= 8 GiB` `line_base[]` table `<< 16 EiB`). Along every fixed outer coordinate tuple `ell = (x_1..x_{15})`, intersecting the quadratic Fincke-Pohst interval `x^T G_red x <= R_gram^2` with the two linear wedge inequalities of `C_30` (`Re(sigma_0(a)) > 0, -tan(pi/30)*Re <= Im < tan(pi/30)*Re`) yields a **SINGLE CONTIGUOUS integer `x_0`-interval `[x_{0,lo}'(ell), x_{0,hi}'(ell)]`** (verified live on real 16D `Q_CM`-LLL basis `B_bc` in `e16_verify_fvor_core`: `0` gaps, `0` collisions, exact flat bijection!). A single tile pre-pass packs non-empty clipped lines consecutively and stores `line_base[ell] = flat_offset - x_{0,lo}'(ell)` (`<= 20` ops/line `< 0.3` ops/cell). During sieving of every degree-1 ideal `p in [2^8, 2^{68.25}]`, Bareiss LLL (`PAIR_SETUP = 2^{35.2}`) + 16D Fincke-Pohst Q32 DFS (`<= 650`) + leaf materialization (`64`) + `C_30`-interval check & flat-address lookup `addr = line_base[ell] + x_0` (`<= 120` Ir/hit, replacing the older `320`-Ir coordinate-ball rank table!) + routing/update/spill (`260`) = **`<= 1,094 <= 1,400 = SIEVE_HIT`** Ir/hit (`F_Vor(a, b)` is tested per cell inside `BOUND_CELL = 4,400`, which includes `+200` Ir/cell Fincke-Pohst/`C_30`/`F_Vor` overhead).
- **Certified Outward-Rounded `BOUND_CELL (+6)` (`4,400` Ir/cell, `q15_descent16.py`):** Signed-64b `bmul_6` (`+6 >= ceil(sqrt(5)+3)`, `<= +0.5343b < +0.625b` `SQ_LUT` envelope): `3,792` exact ball-mul core + `200` `E_16` Fincke-Pohst coordinate-enumeration surcharge + `128 + 80 + 64 + 136 = 4,400` Ir/cell (50-digit `Decimal` witness & 12-cell real-lattice `bc_gmin >= 0` verified live).
- **100% Input-Uniform 2-Stage Pure-Int Bareiss LLL (`11/40 -> 99/100`, ZERO Float LLL):** Every `L_Q` (`det = q <= 2^260`) and `G_T` boot lattice (`S_boot = 2^21` scaling) runs `11/40` Euclidean -> `99/100` `Q_CM = diag(G_8, G_8)` (`det(G_8) = 3^4 * 5^6`, `lambda_max(G_8) = 15` proved live via Sylvester minors): **`PAIR_SETUP = 2^35.2`** (`2^35.06`), **`LLL_NODE = 2^40.0`** (`2^39.46`), **`LLL_BOOT = 2^47.0`** (`2^46.75`), alongside `SCAN = 10`, `ZERO = 2`, `NORMS = 2^27`, `TDIV = 2^25`, `PRIME = 2^24`, `IDEAL = 2^35`.

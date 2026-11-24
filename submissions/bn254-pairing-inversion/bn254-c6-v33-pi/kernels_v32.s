// kernels_v32.s - complete AArch64 scalar ARMv8.1-A (+ FEAT_LSE LDADD/STADDH) sequences priced in KERNELS.md section 1 for v32c:
//   (1) v26 FP nodes extended with exact 6-trit O_h->F_3^3 primitivity-wheel state tracking (x23=(cnt<<12)|st_1, STEP_K[1..11] at sp+8192*k, STEP1_TAB at sp+4096),
//   (2) Zero-Split 16-byte a_tab=(ai,&B_TAB[k][0]) + contiguous 4-plane B_TAB[0..3] + 9-op LDP+LDADD direct-to-link-page emit (x24=&a_tab[pack_a], x25=pack_b; ldp x0,x6,[x2]; ldr x3,[x6,x5,lsl #3]; add x8,x17,x3,lsr #26; add x7,x0,x3,lsl #38; ldadd x1,x9,[x8]),
//   (3) period-3 chord mask dispatch + computed-goto straight-line emit chains (C12/C21/C3, wtab[0..447]).
// Part 1 holds the complete FP tree node bodies (with .Lchord_split the inlined chord-split/wheel site); Part 2 holds the priced kernel templates (as in v26).
// Every //@seg block is counted directly by audit_bn_cm.py (v20_lat_b2.fp_prices_from_kernels + wheel_prices_from_kernels).
// Part 1: complete Fincke-Pohst node bodies (12 levels, 2-op Direct Q12 isqrt LUT (2^37 entries = 1 TiB in x16)). audit_bn_cm.py counts the instruction lines in each
//   //@seg block and composes the node prices exactly as in v20_lat_b2.py FP_CL (k=0: k0_body+k0_back; k=1: top_2+body_2+k1_exit;
//   k>=2 at level m: top_{m+1}+body_{m+1}+top_m). Part 2: leaf emit (10 per kept leaf, 9 for the last leaf of a chain), x8-unrolled
//   apply (14 per 8 records), Lower Band 2 x4 loop (22 per 4 hits), a_0-chord kta pass 0 (20 per piece + 6 per entry), and the
//   complete Band-1 path (//@seg b1_*: line setup, row, transform-copy, page flush, stream switch on L[k]; page apply on C[s]).
// Node bodies of v20_lat_b2.py full_tree(); general registers only, x18 unused; [sp, #0] holds the exact-level mask.
// Registers: x0-x14 scratch (chord split may clobber them), x15-x17 chord-split bases (atab,btab,wptr;
// never written here), x19-x30 k=0 loop state (x19 rem0, x20 dr1, x21 A, x22 B, x23 cnt,
// x24-x25 cur (2 Z/2^64Z mixed-radix words: pack_a, pack_b<<5), x26-x27=(da, db5), x28=2*A1, x29=mu10, x30 inv_b2[0]; x1=#8, x10=mod3_tab, x11-x12=(3*da, 3*db5), x13=wtab pinned in body_2); in levels >= 3 code x30 = frame base xb.
// LUT based at sp (register-offset ldr [sp, xi, lsl #3]); indices >= 512 (7-op LUT) / >= 50176 (lut4) never overlap the sp block [0,480)
.arch armv8.1-a
.text
//@seg top_11
seg_top_11:
.Ltop_11:
 add x30, sp, #2624
 ldr x0, [sp, #88]
 ldp x1, x2, [x30, #0]
 add x0, x0, #1
 cmp x0, x1
 b.gt .Lorbit_done
 ldrh w9, [sp, #3862]
 add x10, sp, #22, lsl #12
 ldrh w9, [x10, x9, lsl #1]
 strh w9, [sp, #3862]
//@end
//@seg body_11
seg_body_11:
.Lent_11:
 ldr x3, [sp]
 tbnz x3, #11, .Lexact_11
 str x0, [sp, #88]
 subs x3, x2, x0, lsl #32
 cneg x3, x3, mi
 ldp x4, x5, [x30, #16]
 subs x3, x3, x4
 csel x3, x3, xzr, hi
 lsr x3, x3, #13
 mul x3, x3, x3
 umulh x3, x3, x5
 ldp x6, x7, [x30, #32]
 subs x6, x6, x3
 b.lo .Ltop_11
 ldp x8, x9, [x30, #48]
 msub x2, x8, x0, x7
 ldr x10, [x30, #176]
 mneg x7, x10, x0
 stp x6, x7, [x30, #-224]
 ldp x10, x11, [x30, #80]
 ldp x20, x21, [x30, #128]
 madd x10, x0, x20, x10
 madd x11, x0, x21, x11
 stp x10, x11, [x30, #-176]
 umulh x10, x6, x9
 ldr x11, [x16, x10, lsl #3]
 ldp x8, x9, [x30, #64]
 add x0, x2, x9
 sub x0, x0, x11
 asr x0, x0, #32
 add x1, x2, x8
 add x1, x1, x11
 asr x1, x1, #32
 cmp x0, x1
 b.gt .Ltop_11
 stp x1, x2, [x30, #-256]
 add x30, sp, #2368
 ldr x11, [sp, #8]
 ldrb w9, [x11, w0, sxtw]
 ldrh w10, [sp, #3862]
 add x11, sp, #20, lsl #12
 add x11, x11, x9, lsl #11
 ldrh w10, [x11, x10, lsl #1]
 strh w10, [sp, #3860]
 b .Lent_10
//@end
//@seg top_10
seg_top_10:
.Ltop_10:
 add x30, sp, #2368
 ldr x0, [sp, #80]
 ldp x1, x2, [x30, #0]
 add x0, x0, #1
 cmp x0, x1
 b.gt .Ltop_11
 ldrh w9, [sp, #3860]
 add x10, sp, #20, lsl #12
 ldrh w9, [x10, x9, lsl #1]
 strh w9, [sp, #3860]
//@end
//@seg body_10
seg_body_10:
.Lent_10:
 ldr x3, [sp]
 tbnz x3, #10, .Lexact_10
 str x0, [sp, #80]
 subs x3, x2, x0, lsl #32
 cneg x3, x3, mi
 ldp x4, x5, [x30, #16]
 subs x3, x3, x4
 csel x3, x3, xzr, hi
 lsr x3, x3, #13
 mul x3, x3, x3
 umulh x3, x3, x5
 ldp x6, x7, [x30, #32]
 subs x6, x6, x3
 b.lo .Ltop_10
 ldp x8, x9, [x30, #48]
 msub x2, x8, x0, x7
 ldr x10, [sp, #88]
 ldp x11, x12, [x30, #176]
 mneg x7, x11, x0
 msub x7, x12, x10, x7
 stp x6, x7, [x30, #-224]
 ldp x10, x11, [x30, #80]
 ldp x20, x21, [x30, #128]
 madd x10, x0, x20, x10
 madd x11, x0, x21, x11
 stp x10, x11, [x30, #-176]
 umulh x10, x6, x9
 ldr x11, [x16, x10, lsl #3]
 ldp x8, x9, [x30, #64]
 add x0, x2, x9
 sub x0, x0, x11
 asr x0, x0, #32
 add x1, x2, x8
 add x1, x1, x11
 asr x1, x1, #32
 cmp x0, x1
 b.gt .Ltop_10
 stp x1, x2, [x30, #-256]
 add x30, sp, #2112
 ldr x11, [sp, #8]
 ldrb w9, [x11, w0, sxtw]
 ldrh w10, [sp, #3860]
 add x11, sp, #18, lsl #12
 add x11, x11, x9, lsl #11
 ldrh w10, [x11, x10, lsl #1]
 strh w10, [sp, #3858]
 b .Lent_9
//@end
//@seg top_9
seg_top_9:
.Ltop_9:
 add x30, sp, #2112
 ldr x0, [sp, #72]
 ldp x1, x2, [x30, #0]
 add x0, x0, #1
 cmp x0, x1
 b.gt .Ltop_10
 ldrh w9, [sp, #3858]
 add x10, sp, #18, lsl #12
 ldrh w9, [x10, x9, lsl #1]
 strh w9, [sp, #3858]
//@end
//@seg body_9
seg_body_9:
.Lent_9:
 ldr x3, [sp]
 tbnz x3, #9, .Lexact_9
 str x0, [sp, #72]
 subs x3, x2, x0, lsl #32
 cneg x3, x3, mi
 ldp x4, x5, [x30, #16]
 subs x3, x3, x4
 csel x3, x3, xzr, hi
 lsr x3, x3, #13
 mul x3, x3, x3
 umulh x3, x3, x5
 ldp x6, x7, [x30, #32]
 subs x6, x6, x3
 b.lo .Ltop_9
 ldp x8, x9, [x30, #48]
 msub x2, x8, x0, x7
 ldp x10, x11, [sp, #80]
 ldp x12, x13, [x30, #176]
 mneg x7, x12, x0
 msub x7, x13, x10, x7
 ldr x12, [x30, #192]
 msub x7, x12, x11, x7
 stp x6, x7, [x30, #-224]
 ldp x10, x11, [x30, #80]
 ldp x20, x21, [x30, #128]
 madd x10, x0, x20, x10
 madd x11, x0, x21, x11
 stp x10, x11, [x30, #-176]
 umulh x10, x6, x9
 ldr x11, [x16, x10, lsl #3]
 ldp x8, x9, [x30, #64]
 add x0, x2, x9
 sub x0, x0, x11
 asr x0, x0, #32
 add x1, x2, x8
 add x1, x1, x11
 asr x1, x1, #32
 cmp x0, x1
 b.gt .Ltop_9
 stp x1, x2, [x30, #-256]
 add x30, sp, #1856
 ldr x11, [sp, #8]
 ldrb w9, [x11, w0, sxtw]
 ldrh w10, [sp, #3858]
 add x11, sp, #16, lsl #12
 add x11, x11, x9, lsl #11
 ldrh w10, [x11, x10, lsl #1]
 strh w10, [sp, #3856]
 b .Lent_8
//@end
//@seg top_8
seg_top_8:
.Ltop_8:
 add x30, sp, #1856
 ldr x0, [sp, #64]
 ldp x1, x2, [x30, #0]
 add x0, x0, #1
 cmp x0, x1
 b.gt .Ltop_9
 ldrh w9, [sp, #3856]
 add x10, sp, #16, lsl #12
 ldrh w9, [x10, x9, lsl #1]
 strh w9, [sp, #3856]
//@end
//@seg body_8
seg_body_8:
.Lent_8:
 ldr x3, [sp]
 tbnz x3, #8, .Lexact_8
 str x0, [sp, #64]
 subs x3, x2, x0, lsl #32
 cneg x3, x3, mi
 ldp x4, x5, [x30, #16]
 subs x3, x3, x4
 csel x3, x3, xzr, hi
 lsr x3, x3, #13
 mul x3, x3, x3
 umulh x3, x3, x5
 ldp x6, x7, [x30, #32]
 subs x6, x6, x3
 b.lo .Ltop_8
 ldp x8, x9, [x30, #48]
 msub x2, x8, x0, x7
 ldp x10, x11, [sp, #72]
 ldr x12, [sp, #88]
 ldp x13, x14, [x30, #176]
 mneg x7, x13, x0
 msub x7, x14, x10, x7
 ldp x13, x14, [x30, #192]
 msub x7, x13, x11, x7
 msub x7, x14, x12, x7
 stp x6, x7, [x30, #-224]
 ldp x10, x11, [x30, #80]
 ldp x20, x21, [x30, #128]
 madd x10, x0, x20, x10
 madd x11, x0, x21, x11
 stp x10, x11, [x30, #-176]
 umulh x10, x6, x9
 ldr x11, [x16, x10, lsl #3]
 ldp x8, x9, [x30, #64]
 add x0, x2, x9
 sub x0, x0, x11
 asr x0, x0, #32
 add x1, x2, x8
 add x1, x1, x11
 asr x1, x1, #32
 cmp x0, x1
 b.gt .Ltop_8
 stp x1, x2, [x30, #-256]
 add x30, sp, #1600
 ldr x11, [sp, #8]
 ldrb w9, [x11, w0, sxtw]
 ldrh w10, [sp, #3856]
 add x11, sp, #14, lsl #12
 add x11, x11, x9, lsl #11
 ldrh w10, [x11, x10, lsl #1]
 strh w10, [sp, #3854]
 b .Lent_7
//@end
//@seg top_7
seg_top_7:
.Ltop_7:
 add x30, sp, #1600
 ldr x0, [sp, #56]
 ldp x1, x2, [x30, #0]
 add x0, x0, #1
 cmp x0, x1
 b.gt .Ltop_8
 ldrh w9, [sp, #3854]
 add x10, sp, #14, lsl #12
 ldrh w9, [x10, x9, lsl #1]
 strh w9, [sp, #3854]
//@end
//@seg body_7
seg_body_7:
.Lent_7:
 ldr x3, [sp]
 tbnz x3, #7, .Lexact_7
 str x0, [sp, #56]
 subs x3, x2, x0, lsl #32
 cneg x3, x3, mi
 ldp x4, x5, [x30, #16]
 subs x3, x3, x4
 csel x3, x3, xzr, hi
 lsr x3, x3, #13
 mul x3, x3, x3
 umulh x3, x3, x5
 ldp x6, x7, [x30, #32]
 subs x6, x6, x3
 b.lo .Ltop_7
 ldp x8, x9, [x30, #48]
 msub x2, x8, x0, x7
 ldp x10, x11, [sp, #64]
 ldp x12, x13, [sp, #80]
 ldp x14, x19, [x30, #176]
 mneg x7, x14, x0
 msub x7, x19, x10, x7
 ldp x14, x19, [x30, #192]
 msub x7, x14, x11, x7
 msub x7, x19, x12, x7
 ldr x14, [x30, #208]
 msub x7, x14, x13, x7
 stp x6, x7, [x30, #-224]
 ldp x10, x11, [x30, #80]
 ldp x20, x21, [x30, #128]
 madd x10, x0, x20, x10
 madd x11, x0, x21, x11
 stp x10, x11, [x30, #-176]
 umulh x10, x6, x9
 ldr x11, [x16, x10, lsl #3]
 ldp x8, x9, [x30, #64]
 add x0, x2, x9
 sub x0, x0, x11
 asr x0, x0, #32
 add x1, x2, x8
 add x1, x1, x11
 asr x1, x1, #32
 cmp x0, x1
 b.gt .Ltop_7
 stp x1, x2, [x30, #-256]
 add x30, sp, #1344
 ldr x11, [sp, #8]
 ldrb w9, [x11, w0, sxtw]
 ldrh w10, [sp, #3854]
 add x11, sp, #12, lsl #12
 add x11, x11, x9, lsl #11
 ldrh w10, [x11, x10, lsl #1]
 strh w10, [sp, #3852]
 b .Lent_6
//@end
//@seg top_6
seg_top_6:
.Ltop_6:
 add x30, sp, #1344
 ldr x0, [sp, #48]
 ldp x1, x2, [x30, #0]
 add x0, x0, #1
 cmp x0, x1
 b.gt .Ltop_7
 ldrh w9, [sp, #3852]
 add x10, sp, #12, lsl #12
 ldrh w9, [x10, x9, lsl #1]
 strh w9, [sp, #3852]
//@end
//@seg body_6
seg_body_6:
.Lent_6:
 ldr x3, [sp]
 tbnz x3, #6, .Lexact_6
 str x0, [sp, #48]
 subs x3, x2, x0, lsl #32
 cneg x3, x3, mi
 ldp x4, x5, [x30, #16]
 subs x3, x3, x4
 csel x3, x3, xzr, hi
 lsr x3, x3, #13
 mul x3, x3, x3
 umulh x3, x3, x5
 ldp x6, x7, [x30, #32]
 subs x6, x6, x3
 b.lo .Ltop_6
 ldp x8, x9, [x30, #48]
 msub x2, x8, x0, x7
 ldp x10, x11, [sp, #56]
 ldp x12, x13, [sp, #72]
 ldr x14, [sp, #88]
 ldp x19, x20, [x30, #176]
 mneg x7, x19, x0
 msub x7, x20, x10, x7
 ldp x19, x20, [x30, #192]
 msub x7, x19, x11, x7
 msub x7, x20, x12, x7
 ldp x19, x20, [x30, #208]
 msub x7, x19, x13, x7
 msub x7, x20, x14, x7
 stp x6, x7, [x30, #-224]
 ldp x10, x11, [x30, #80]
 ldp x20, x21, [x30, #128]
 madd x10, x0, x20, x10
 madd x11, x0, x21, x11
 stp x10, x11, [x30, #-176]
 umulh x10, x6, x9
 ldr x11, [x16, x10, lsl #3]
 ldp x8, x9, [x30, #64]
 add x0, x2, x9
 sub x0, x0, x11
 asr x0, x0, #32
 add x1, x2, x8
 add x1, x1, x11
 asr x1, x1, #32
 cmp x0, x1
 b.gt .Ltop_6
 stp x1, x2, [x30, #-256]
 add x30, sp, #1088
 ldr x11, [sp, #8]
 ldrb w9, [x11, w0, sxtw]
 ldrh w10, [sp, #3852]
 add x11, sp, #10, lsl #12
 add x11, x11, x9, lsl #11
 ldrh w10, [x11, x10, lsl #1]
 strh w10, [sp, #3850]
 b .Lent_5
//@end
//@seg top_5
seg_top_5:
.Ltop_5:
 add x30, sp, #1088
 ldr x0, [sp, #40]
 ldp x1, x2, [x30, #0]
 add x0, x0, #1
 cmp x0, x1
 b.gt .Ltop_6
 ldrh w9, [sp, #3850]
 add x10, sp, #10, lsl #12
 ldrh w9, [x10, x9, lsl #1]
 strh w9, [sp, #3850]
//@end
//@seg body_5
seg_body_5:
.Lent_5:
 ldr x3, [sp]
 tbnz x3, #5, .Lexact_5
 str x0, [sp, #40]
 subs x3, x2, x0, lsl #32
 cneg x3, x3, mi
 ldp x4, x5, [x30, #16]
 subs x3, x3, x4
 csel x3, x3, xzr, hi
 lsr x3, x3, #13
 mul x3, x3, x3
 umulh x3, x3, x5
 ldp x6, x7, [x30, #32]
 subs x6, x6, x3
 b.lo .Ltop_5
 ldp x8, x9, [x30, #48]
 msub x2, x8, x0, x7
 ldp x10, x11, [sp, #48]
 ldp x12, x13, [sp, #64]
 ldp x14, x19, [sp, #80]
 ldp x20, x21, [x30, #176]
 mneg x7, x20, x0
 msub x7, x21, x10, x7
 ldp x20, x21, [x30, #192]
 msub x7, x20, x11, x7
 msub x7, x21, x12, x7
 ldp x20, x21, [x30, #208]
 msub x7, x20, x13, x7
 msub x7, x21, x14, x7
 ldr x20, [x30, #224]
 msub x7, x20, x19, x7
 stp x6, x7, [x30, #-224]
 ldp x10, x11, [x30, #80]
 ldp x20, x21, [x30, #128]
 madd x10, x0, x20, x10
 madd x11, x0, x21, x11
 stp x10, x11, [x30, #-176]
 umulh x10, x6, x9
 ldr x11, [x16, x10, lsl #3]
 ldp x8, x9, [x30, #64]
 add x0, x2, x9
 sub x0, x0, x11
 asr x0, x0, #32
 add x1, x2, x8
 add x1, x1, x11
 asr x1, x1, #32
 cmp x0, x1
 b.gt .Ltop_5
 stp x1, x2, [x30, #-256]
 add x30, sp, #832
 ldr x11, [sp, #8]
 ldrb w9, [x11, w0, sxtw]
 ldrh w10, [sp, #3850]
 add x11, sp, #8, lsl #12
 add x11, x11, x9, lsl #11
 ldrh w10, [x11, x10, lsl #1]
 strh w10, [sp, #3848]
 b .Lent_4
//@end
//@seg top_4
seg_top_4:
.Ltop_4:
 add x30, sp, #832
 ldr x0, [sp, #32]
 ldp x1, x2, [x30, #0]
 add x0, x0, #1
 cmp x0, x1
 b.gt .Ltop_5
 ldrh w9, [sp, #3848]
 add x10, sp, #8, lsl #12
 ldrh w9, [x10, x9, lsl #1]
 strh w9, [sp, #3848]
//@end
//@seg body_4
seg_body_4:
.Lent_4:
 ldr x3, [sp]
 tbnz x3, #4, .Lexact_4
 str x0, [sp, #32]
 subs x3, x2, x0, lsl #32
 cneg x3, x3, mi
 ldp x4, x5, [x30, #16]
 subs x3, x3, x4
 csel x3, x3, xzr, hi
 lsr x3, x3, #13
 mul x3, x3, x3
 umulh x3, x3, x5
 ldp x6, x7, [x30, #32]
 subs x6, x6, x3
 b.lo .Ltop_4
 ldp x8, x9, [x30, #48]
 msub x2, x8, x0, x7
 ldp x10, x11, [sp, #40]
 ldp x12, x13, [sp, #56]
 ldp x14, x19, [sp, #72]
 ldr x20, [sp, #88]
 ldp x21, x22, [x30, #176]
 mneg x7, x21, x0
 msub x7, x22, x10, x7
 ldp x21, x22, [x30, #192]
 msub x7, x21, x11, x7
 msub x7, x22, x12, x7
 ldp x21, x22, [x30, #208]
 msub x7, x21, x13, x7
 msub x7, x22, x14, x7
 ldp x21, x22, [x30, #224]
 msub x7, x21, x19, x7
 msub x7, x22, x20, x7
 stp x6, x7, [x30, #-288]
 ldp x10, x11, [x30, #80]
 ldp x20, x21, [x30, #128]
 madd x10, x0, x20, x10
 madd x11, x0, x21, x11
 stp x10, x11, [x30, #-240]
 umulh x10, x6, x9
 ldr x11, [x16, x10, lsl #3]
 ldp x8, x9, [x30, #64]
 add x0, x2, x9
 sub x0, x0, x11
 asr x0, x0, #32
 add x1, x2, x8
 add x1, x1, x11
 asr x1, x1, #32
 cmp x0, x1
 b.gt .Ltop_4
 stp x1, x2, [x30, #-320]
 add x30, sp, #512
 ldr x11, [sp, #8]
 ldrb w9, [x11, w0, sxtw]
 ldrh w10, [sp, #3848]
 add x11, sp, #6, lsl #12
 add x11, x11, x9, lsl #11
 ldrh w10, [x11, x10, lsl #1]
 strh w10, [sp, #3846]
 b .Lent_3
//@end
//@seg top_3
seg_top_3:
.Ltop_3:
 add x30, sp, #512
 ldr x0, [sp, #24]
 ldp x1, x2, [x30, #0]
 add x0, x0, #1
 cmp x0, x1
 b.gt .Ltop_4
 ldrh w9, [sp, #3846]
 add x10, sp, #6, lsl #12
 ldrh w9, [x10, x9, lsl #1]
 strh w9, [sp, #3846]
//@end
//@seg body_3
seg_body_3:
.Lent_3:
 ldr x3, [sp]
 tbnz x3, #3, .Lexact_3
 str x0, [sp, #24]
 subs x3, x2, x0, lsl #32
 cneg x3, x3, mi
 ldp x4, x5, [x30, #16]
 subs x3, x3, x4
 csel x3, x3, xzr, hi
 lsr x3, x3, #13
 mul x3, x3, x3
 umulh x3, x3, x5
 ldp x6, x7, [x30, #32]
 subs x6, x6, x3
 b.lo .Ltop_3
 ldp x8, x9, [x30, #48]
 msub x2, x8, x0, x7
 ldp x10, x11, [sp, #32]
 ldp x12, x13, [sp, #48]
 ldp x14, x19, [sp, #64]
 ldp x20, x21, [sp, #80]
 ldp x22, x23, [x30, #176]
 mneg x7, x22, x0
 msub x7, x23, x10, x7
 ldp x22, x23, [x30, #192]
 msub x7, x22, x11, x7
 msub x7, x23, x12, x7
 ldp x22, x23, [x30, #208]
 msub x7, x22, x13, x7
 msub x7, x23, x14, x7
 ldp x22, x23, [x30, #224]
 msub x7, x22, x19, x7
 msub x7, x23, x20, x7
 ldr x22, [x30, #240]
 msub x7, x22, x21, x7
 ldp x22, x23, [x30, #248]
 mneg x4, x22, x0
 msub x4, x23, x10, x4
 ldp x22, x23, [x30, #264]
 msub x4, x22, x11, x4
 msub x4, x23, x12, x4
 ldp x22, x23, [x30, #280]
 msub x4, x22, x13, x4
 msub x4, x23, x14, x4
 ldp x22, x23, [x30, #296]
 msub x4, x22, x19, x4
 msub x4, x23, x20, x4
 ldr x22, [x30, #312]
 msub x4, x22, x21, x4
 stp x6, x7, [sp, #336]
 str x4, [sp, #376]
 ldp x10, x11, [x30, #80]
 ldp x20, x21, [x30, #128]
 madd x10, x0, x20, x10
 madd x11, x0, x21, x11
 stp x10, x11, [sp, #384]
 umulh x10, x6, x9
 ldr x11, [x16, x10, lsl #3]
 ldp x8, x9, [x30, #64]
 add x0, x2, x9
 sub x0, x0, x11
 asr x0, x0, #32
 add x1, x2, x8
 add x1, x1, x11
 asr x1, x1, #32
 cmp x0, x1
 b.gt .Ltop_3
 stp x1, x2, [sp, #304]
 ldr x11, [sp, #8]
 ldrb w9, [x11, w0, sxtw]
 ldrh w10, [sp, #3846]
 add x11, sp, #4, lsl #12
 add x11, x11, x9, lsl #11
 ldrh w10, [x11, x10, lsl #1]
 strh w10, [sp, #3844]
 b .Lent_2
//@end
//@seg top_2
seg_top_2:
.Ltop_2:
 ldr x0, [sp, #16]
 ldp x1, x2, [sp, #304]
 add x0, x0, #1
 cmp x0, x1
 b.gt .Ltop_3
 ldrh w9, [sp, #3844]
 add x10, sp, #4, lsl #12
 ldrh w9, [x10, x9, lsl #1]
 strh w9, [sp, #3844]
//@end
//@seg body_2
seg_body_2:
.Lent_2:
 str x0, [sp, #16]
 subs x3, x2, x0, lsl #32
 cneg x3, x3, mi
 ldp x4, x5, [sp, #320]
 subs x3, x3, x4
 csel x3, x3, xzr, hi
 lsr x3, x3, #13
 mul x3, x3, x3
 umulh x3, x3, x5
 ldp x6, x7, [sp, #336]
 subs x6, x6, x3
 b.lo .Ltop_2
 ldp x8, x9, [sp, #352]
 msub x2, x8, x0, x7
 ldp x8, x10, [sp, #368]
 msub x10, x8, x0, x10
 umulh x11, x6, x9
 ldr x12, [x16, x11, lsl #3]
 ldp x8, x9, [sp, #272]
 add x13, x2, x9
 sub x13, x13, x12
 asr x13, x13, #32
 add x14, x2, x8
 add x14, x14, x12
 asr x14, x14, #32
 cmp x13, x14
 b.gt .Ltop_2
 ldp x8, x9, [sp, #208]
 sub x23, x14, x13
 madd x11, x23, x8, x9
 add x12, x23, #1
 mul x12, x12, x12
 add x19, x6, x11, lsr #32
 add x19, x19, x12
 subs x3, x2, x13, lsl #32
 cneg x4, x3, mi
 lsr x4, x4, #13
 mul x4, x4, x4
 ldp x8, x9, [sp, #224]
 umulh x4, x4, x8
 sub x19, x19, x4
 ldp x8, x30, [sp, #240]
 sub x5, x8, x3, lsl #1
 mul x11, x5, x9
 smulh x12, x5, x9
 extr x20, x12, x11, #32
 msub x10, x29, x13, x10
 ldp x8, x9, [sp, #256]
 add x21, x10, x8
 add x22, x10, x9
 ldp x24, x25, [sp, #384]
 ldp x3, x4, [sp, #432]
 madd x24, x0, x3, x24
 madd x25, x0, x4, x25
 ldp x3, x4, [sp, #160]
 madd x24, x13, x3, x24
 madd x25, x13, x4, x25
 ldr x11, [sp, #8]
 ldrb w9, [x11, w13, sxtw]
 ldrh w10, [sp, #3844]
 add x11, sp, #2, lsl #12
 add x11, x11, x9, lsl #11
 ldrh w14, [x11, x10, lsl #1]
 ldp x11, x12, [sp, #128]
 mov x1, #8
 ldr x10, [sp, #408]
 add x13, sp, #32, lsl #12
//@end
//@seg k0_body
seg_k0_body:
.Lk0_top:
 ldrh w14, [x15, x14, lsl #1]
 tbnz x19, #63, .Lk0_skip
 umulh x8, x19, x30
 ldr x8, [x16, x8, lsl #3]
 sub x6, x21, x8
 asr x6, x6, #32
 add x7, x22, x8
 asr x7, x7, #32
 subs x0, x7, x6
 b.lt .Lk0_skip
 madd x2, x6, x26, x24
 madd x5, x6, x27, x25
 sub x19, x19, x20
 add x20, x20, x28
 sub x21, x21, x29
 sub x22, x22, x29
 ldp x8, x9, [sp, #160]
 add x24, x24, x8
 add x25, x25, x9
 sub x23, x23, #1
//@end
.ifdef FP_TEST
 stp x2, x5, [x15], #16
 stp x6, x7, [x15], #16
.endif
.Lchord_split:
//@seg k0_back
seg_k0_back:
.Lk0_back:
 tbz x23, #63, .Lk0_top
//@end
//@seg k1_exit
seg_k1_exit:
 b .Ltop_2
//@end
//@seg k0_skip
seg_k0_skip:
.Lk0_skip:
 sub x19, x19, x20
 add x20, x20, x28
 sub x21, x21, x29
 sub x22, x22, x29
 ldp x8, x9, [sp, #160]
 add x24, x24, x8
 add x25, x25, x9
 sub x23, x23, #1
 tbz x23, #63, .Lk0_top
 b .Ltop_2
//@end
orbit_done:
.Lorbit_done:
 b fp_orbit_done
// Levels with beta_k = B_k/T >= 1 (possible only at k >= 3, fd_bounds.py Lemma E; bit k of the per-orbit mask at [sp])
// enter the exact-integer body (KERNELS.md 1.3), priced per node by fd_bounds.py exact_node_ops; not listed here.
.Lexact_3:
.Lexact_4:
.Lexact_5:
.Lexact_6:
.Lexact_7:
.Lexact_8:
.Lexact_9:
.Lexact_10:
.Lexact_11:
 b fp_exact_node
.global pass0
pass0:
.Lpiece:
  ldp  x8, x9, [x0], #16
  sxth x10, w8
  mul  x11, x10, x20
  sbfx x10, x8, #16, #16
  madd x11, x10, x21, x11
  sbfx x10, x8, #32, #16
  madd x11, x10, x22, x11
  asr  x10, x8, #48
  madd x11, x10, x23, x11
  sxth x10, w9
  madd x11, x10, x24, x11
  sbfx x10, x9, #16, #16
  madd x11, x10, x25, x11
  add  x11, x11, x26
  umulh x12, x11, x27
  msub x11, x12, x5, x11
  subs x12, x11, x5
  csel x11, x12, x11, hs
.Lent:
  ldr  x13, [x2], #8
  str  x11, [x3, x13, lsl #3]
  add  x11, x11, x4
  subs x12, x11, x5
  csel x11, x12, x11, hs
  tbz  x13, #61, .Lent
  cmp  x0, x1
  b.lo .Lpiece
  ret
// Zero-Split 16-byte a_tab=(2*ai, &B_TAB[k][0]) + contiguous 4-plane B_TAB[0..3] + 9-op LDP+LDADD emit (BLK=38, uint16_t 2 B/cell):
// x2 = &a_tab[pack_a] (Z/2^64Z affine byte address), x5 = pack_b in [0, prod_Mb), x1 = #8, x26 = 16*da, x27 = db, x11 = 3*16*da, x12 = 3*db.
// a_tab[pack_a] = { uint64_t ai; uint64_t *btab_plane_k; } (16 bytes; full 64-bit byte offset ai < 2^48, zero bit-packing).
// B_TAB[k][pack_b] = (s << 29) | row_block (row_block < 2^25 < 2^26 => bits 26..28 == 0; 0 -> C[0] dummy sink if b out of E_b).
// ldp x0, x6, [x2]; ldr x3, [x6, x5, lsl #3]; add x8, x17, x3, lsr #26; add x7, x0, x3, lsl #38; ldadd x1, x9, [x8]; str x7, [x9]; tbnz x9, #14.
//@seg chord_enter
chord_enter:
 // d = x7 - x6 = lv - 1 already computed in k0_body by subs x0, x7, x6; b.lt .Lk0_skip (0 extra ops!)
//@end
//@seg emit9_step1
emit9_step1:
 ldp x0, x6, [x2]
 ldr x3, [x6, x5, lsl #3]
 add x8, x17, x3, lsr #26
 add x7, x0, x3, lsl #38
 ldadd x1, x9, [x8]
 str x7, [x9]
 tbnz x9, #14, .Lfl_s1
 add x2, x2, x26
 add x5, x5, x27
//@end
//@seg emit9_step2
emit9_step2:
 ldp x0, x6, [x2]
 ldr x3, [x6, x5, lsl #3]
 add x8, x17, x3, lsr #26
 add x7, x0, x3, lsl #38
 ldadd x1, x9, [x8]
 str x7, [x9]
 tbnz x9, #14, .Lfl_s2
 add x2, x2, x26, lsl #1
 add x5, x5, x27, lsl #1
//@end
//@seg emit9_step3
emit9_step3:
 ldp x0, x6, [x2]
 ldr x3, [x6, x5, lsl #3]
 add x8, x17, x3, lsr #26
 add x7, x0, x3, lsl #38
 ldadd x1, x9, [x8]
 str x7, [x9]
 tbnz x9, #14, .Lfl_s3
 add x2, x2, x11
 add x5, x5, x12
//@end
//@seg emit9_last
emit9_last:
 ldp x0, x6, [x2]
 ldr x3, [x6, x5, lsl #3]
 add x8, x17, x3, lsr #26
 add x7, x0, x3, lsl #38
 ldadd x1, x9, [x8]
 str x7, [x9]
 tbnz x9, #14, .Lfl_last
//@end
//@seg wheel_m0_check
wheel_m0_check:
 lsr w4, w14, #12
 cbz w4, .Lmask0
//@end
//@seg wheel_m0_disp
wheel_m0_disp:
 cmp w0, #127
 b.hi .Lchunk128_m0
 ldr x9, [x13, w0, uxtw #3]
 br x9
//@end
//@seg wheel_chord_rot
wheel_chord_rot:
 add w9, w6, #0x400, lsl #12
 bfi w9, w4, #23, #3
 ldrb w4, [x10, x9]
//@end
//@seg wheel_seg_disp0
wheel_seg_disp0:
 cmp w0, #127
 b.hi .Lchunk128_s0
 add w8, w0, w4, lsl #7
 ldr x9, [x13, w8, uxtw #3]
 br x9
//@end
//@seg wheel_skip1_stub
wheel_skip1_stub:
 add x2, x2, x26
 add x5, x5, x27
 // falls straight through into dedicated wtab[(cmask<<7)|d] chain (0 branch!)
//@end
//@seg wheel_skip2_stub
wheel_skip2_stub:
 add x2, x2, x26, lsl #1
 add x5, x5, x27, lsl #1
 // falls straight through into dedicated wtab[(cmask<<7)|d] chain (0 branch!)
//@end
.Lfl_s1:
.Lfl_s2:
.Lfl_s3:
.Lfl_last:
.Lmask0:
.Lchunk128_m0:
.Lchunk128_s0:
.Lchain_entry:
 ret
// Direct-to-link-page flush stub (page_start = 32768*s + 8, exact 16,384 B / 2,046 records; preserves x1 = #8):
// 8 ops + 2 header/handoff + <=20 Ir doorbell = <=30 Ir / 2,046 = 0.01466 <= 0.015 = _FL (ZERO 16 KiB copy, ZERO stage_buf[s]).
//@seg emit_direct_flush
emit_direct_flush:
 sub x0, x9, #4, lsl #12
 add x0, x0, #8
 mov x3, #2046
 ldr x7, [sp, #480]
 stp x3, x7, [x0]
 // doorbell store + wfi handoff + return interrupt (<=20 Ir charged in STORAGE.md section 2)
 add x9, x0, #16
 str x9, [x8]
 ret
//@end
// C[s] x128-unrolled direct link-page apply (386 ops per 256 records = 1.5078 Ir/rec; plus header+handoff <= 1.55 Ir/rec):
//@seg lat_apply_page
apply_page:
 ldp x11, x6, [x28]
 add x4, x28, #16
 add x8, x4, x11, lsl #3
 tbz x11, #0, apply_even
 ldr x0, [x4], #8
 staddh w6, [x0]
apply_even:
 ubfx x12, x11, #1, #7
 adr x13, apply_end
 sub x13, x13, x12, lsl #3
 sub x13, x13, x12, lsl #2
 br x13
//@end
//@seg lat_apply
apply:
 .rept 128
 ldp x0, x1, [x4], #16
 staddh w6, [x0]
 staddh w6, [x1]
 .endr
apply_end:
 cmp x4, x8
 b.lo apply
 ret
//@end
// Aligned-sub-band bfi Lower Band 2 (uint16_t STADDH, 2 B/cell, entry = (ai << 5) | j_rel so (x0 >> 4) == 2*ai exact uint16_t byte offset,
// (k_s in [38..49], W in {8, 16}, sgn in {+1, -1}) jump-table dispatched):
// W*2^k_s alignment clears bits [k_s .. k_s + w_bits - 1] with w_bits = log2(W) (3 for W=8, 4 for W=16);
// both pos and neg insert EXACTLY w_bits bits so bit k_s + w_bits of base x22 is never overwritten.
// Sub-band interior end pointer x24 is pre-clipped to the boundary cutoff once in sub-band entry (C_SB_ENTER_X4 = 20.05).
// Exact .rept 8 (x16-unrolled) interior loop: pos = 58/16 = 3.625 ops/hit, neg = 74/16 = 4.625 ops/hit, mean = 4.125 ops/hit (charged C_B2_INT_X4 = 4.50).
//@seg x16_bfi_pos_w8
x16_bfi_pos_w8:
 .rept 8
 ldp x0, x1, [x4], #16
 add x5, x22, x0, lsr #4
 bfi x5, x0, #46, #3
 staddh w23, [x5]
 add x5, x22, x1, lsr #4
 bfi x5, x1, #46, #3
 staddh w23, [x5]
 .endr
 cmp x4, x24
 b.lo x16_bfi_pos_w8
 ret
//@end
//@seg x16_bfi_neg_w8
x16_bfi_neg_w8:
 .rept 8
 ldp x0, x1, [x4], #16
 add x5, x22, x0, lsr #4
 mvn x6, x0
 bfi x5, x6, #46, #3
 staddh w23, [x5]
 add x5, x22, x1, lsr #4
 mvn x6, x1
 bfi x5, x6, #46, #3
 staddh w23, [x5]
 .endr
 cmp x4, x24
 b.lo x16_bfi_neg_w8
 ret
//@end
//@seg x16_bfi_pos_w16
x16_bfi_pos_w16:
 .rept 8
 ldp x0, x1, [x4], #16
 add x5, x22, x0, lsr #4
 bfi x5, x0, #46, #4
 staddh w23, [x5]
 add x5, x22, x1, lsr #4
 bfi x5, x1, #46, #4
 staddh w23, [x5]
 .endr
 cmp x4, x24
 b.lo x16_bfi_pos_w16
 ret
//@end
//@seg x16_bfi_neg_w16
x16_bfi_neg_w16:
 .rept 8
 ldp x0, x1, [x4], #16
 add x5, x22, x0, lsr #4
 mvn x6, x0
 bfi x5, x6, #46, #4
 staddh w23, [x5]
 add x5, x22, x1, lsr #4
 mvn x6, x1
 bfi x5, x6, #46, #4
 staddh w23, [x5]
 .endr
 cmp x4, x24
 b.lo x16_bfi_neg_w16
 ret
//@end
// Band 1 (KERNELS.md section 1.1). One row per b. The List Node L[k] turns each list entry (2*ai byte offset) into the record
// rec = (row_block << 38) + 2*ai, the uint16_t byte address of the cell in C[s] (Lattice Band 2's record format), and stores it
// straight into the link page link_L_to_C[k][s]; C[s] applies whole pages with no per-row work.
// b1_line (once per line of E_b along b_0, at the line's start): x1 = line record (descriptor range; b_0..b_5 as signed 16-bit
// fields), x2 = the prime's c_0..c_5 (c_k < p, residue of b = sum b_k c_k mod p), x5 = p, x6 = log2 of the CHK residue stride,
// x7 = a multiple of p above 6 * 2^15 * p (so the sum is positive before udiv). It sets x21 = (r(b_first) - c_0 mod p) << x6
// and x22 = c_0 << x6 (c_0 is the step of r along b_0), then runs the rows. v25_verify.py check 2 executes it for every line.
//@seg b1_line
b1_line:
 ldp x20, x27, [x1], #16
 ldp x10, x11, [x1], #16
 ldp x12, x13, [x2]
 sbfx x14, x10, #0, #16
 mul x15, x14, x12
 sbfx x14, x10, #16, #16
 madd x15, x14, x13, x15
 ldp x12, x13, [x2, #16]
 sbfx x14, x10, #32, #16
 madd x15, x14, x12, x15
 sbfx x14, x10, #48, #16
 madd x15, x14, x13, x15
 ldp x12, x13, [x2, #32]
 sbfx x14, x11, #0, #16
 madd x15, x14, x12, x15
 sbfx x14, x11, #16, #16
 madd x15, x14, x13, x15
 add x15, x15, x7
 udiv x16, x15, x5
 msub x15, x16, x5, x15
 ldr x12, [x2]
 subs x16, x15, x12
 add x17, x16, x5
 csel x15, x16, x17, hs
 lsl x21, x15, x6
 lsl x22, x12, x6
 b b1_row
//@end
//@seg b1_row
b1_row:
 ldp x9, x10, [x20], #16
 tbnz x10, #31, b1_switch
b1_row_cont:
 add x21, x21, x22
 subs x11, x21, x23
 csel x21, x11, x21, hs
 add x12, x24, x21
 and x13, x10, #0xffff
 lsr x14, x10, #32
 ldr x15, [x12, x13]
 ldr x16, [x12, x14]
 subs x17, x16, x15
 b.eq b1_next
 tbz x17, #3, b1_even
 ldr x0, [x15], #8
 add x0, x9, x0
 str x0, [x25], #8
b1_even:
 ubfx x11, x17, #4, #6
 adr x12, b1_copy_end
 sub x12, x12, x11, lsl #4
 br x12
//@end
//@seg b1_copy
b1_copy:
 .rept 64
 ldp x0, x1, [x15], #16
 add x0, x9, x0
 add x1, x9, x1
 stp x0, x1, [x25], #16
 .endr
b1_copy_end:
 cmp x25, x26
 b.hi b1_flush
b1_pass_room:
 cmp x15, x16
 b.lo b1_copy
b1_next:
 cmp x20, x27
 b.lo b1_row
 ret
//@end
//@seg b1_flush
b1_flush:
 sub x11, x25, x28
 lsr x11, x11, #3
 sub x11, x11, #2
 stp x11, x29, [x28]
b1_handoff:
 add x25, x28, #16
 b b1_pass_room
//@end
//@seg b1_term
b1_term:
 sub x11, x25, x28
 lsr x11, x11, #3
 sub x11, x11, #2
 stp x11, x29, [x28]
b1_term_handoff:
 add x25, x28, #16
 ret
//@end
//@seg b1_switch
b1_switch:
 stp x28, x25, [x19]
 ldr x11, [x20], #8
 add x19, x3, x11, lsl #4
 ldp x28, x25, [x19]
 add x26, x28, #3, lsl #12
 add x26, x26, #3072
 b b1_row_cont
//@end
// Receiver on C[s] at a handoff: x28 = the link page (header: record count, log p weight; then the records).
//@seg b1_page
b1_page:
 ldp x11, x6, [x28]
 add x4, x28, #16
 add x8, x4, x11, lsl #3
 tbz x11, #0, b1_page_even
 ldr x0, [x4], #8
 staddh w6, [x0]
b1_page_even:
 ubfx x12, x11, #1, #7
 adr x13, b1_apply_end
 sub x13, x13, x12, lsl #3
 sub x13, x13, x12, lsl #2
 br x13
//@end
//@seg b1_apply
b1_apply:
 .rept 128
 ldp x0, x1, [x4], #16
 staddh w6, [x0]
 staddh w6, [x1]
 .endr
b1_apply_end:
 cmp x4, x8
 b.lo b1_apply
//@end

// 32-bit delta-packed Band 1 on [2^24, 2^28.5) (uint16_t STADDH, 2 B/cell, delta = 2*(ai_new - ai_old), sub-segment split at build time on byte delta >= 2^32):
// L[k] copies packed uint32_t deltas straight into 16 KiB link pages (2 ops per 4 hits = 0.50 op/hit).
// C[s] unpacks 4 uint32_t deltas per 9 ops (2.25 ops/hit, x32 unrolled = 128 hits/loop); total = 2.75 ops/hit (charged 2.82 + 24*P_ov).
//@seg delta_b1_copy
delta_b1_copy:
 .rept 16
 ldp x0, x1, [x15], #16
 stp x0, x1, [x25], #16
 .endr
delta_b1_copy_end:
 cmp x15, x16
 b.lo delta_b1_copy
 ret
//@end
//@seg delta_b1_apply
delta_b1_apply:
 .rept 32
 ldp x0, x1, [x4], #16
 add x5, x5, w0, uxtw
 staddh w6, [x5]
 add x5, x5, x0, lsr #32
 staddh w6, [x5]
 add x5, x5, w1, uxtw
 staddh w6, [x5]
 add x5, x5, x1, lsr #32
 staddh w6, [x5]
 .endr
delta_b1_apply_end:
 cmp x4, x8
 b.lo delta_b1_apply
 ret
//@end

// kernels_v32.s - complete AArch64 scalar ARMv8.1-A (+ FEAT_LSE LDADD/STADDH) sequences priced in KERNELS.md section 1 for v32c:
//   (1) v26 FP nodes extended with exact 6-trit O_h->F_3^3 primitivity-wheel state tracking (x23=(cnt<<12)|st_1, STEP_K[1..11] at sp+8192*k, STEP1_TAB at sp+4096),
//   (2) Zero-Split 16-byte a_tab=(ai,&B_TAB[k][0]) + contiguous 4-plane B_TAB[0..3] + 9-op LDP+LDADD direct-to-link-page emit (x24=&a_tab[pack_a], x25=pack_b; ldp x0,x6,[x2]; ldr x3,[x6,x5,lsl #3]; add x8,x17,x3,lsr #26; add x7,x0,x3,lsl #38; ldadd x1,x9,[x8]),
//   (3) period-3 chord mask dispatch + computed-goto straight-line emit chains (C12/C21/C3, wtab[0..447]).
// Part 1 holds the complete FP tree node bodies (with .Lchord_split the inlined chord-split/wheel site); Part 2 holds the priced kernel templates (as in v26).
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
//@seg chord_enter
chord_enter:
//@end
//@seg emit9_step1
emit9_step1:
 ldp x0, x6, [x2]
 ldr x3, [x6, x5, lsl #3]
 add x8, x17, x3, lsr #27
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
 add x8, x17, x3, lsr #27
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
 add x8, x17, x3, lsr #27
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
 add x8, x17, x3, lsr #27
 add x7, x0, x3, lsl #38
 ldadd x1, x9, [x8]
 str x7, [x9]
 tbnz x9, #14, .Lfl_last
//@end
//@seg wheel_disp4
wheel_disp4:
 bfi w6, w14, #23, #3
 ldr x8, [x10, w6, sxtw #3]
 ldr x9, [x8, w0, uxtw #3]
 br x9
//@end
//@seg wheel_skip1_stub
wheel_skip1_stub:
 add x2, x2, x26
 add x5, x5, x27
//@end
//@seg wheel_skip2_stub
wheel_skip2_stub:
 add x2, x2, x26, lsl #1
 add x5, x5, x27, lsl #1
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
//@seg emit_direct_flush
emit_direct_flush:
 sub x0, x9, #4, lsl #12
 add x0, x0, #8
 mov x3, #2046
 ldr x7, [sp, #480]
 stp x3, x7, [x0]
 add x9, x0, #16
 str x9, [x8]
 ret
//@end
//@seg lat_apply_page
apply_page:
 ldp x11, x6, [x28]
 add x4, x28, #16
 add x8, x4, x11, lsl #3
 tbz x11, #0, apply_even
 ldr x0, [x4], #8
 staddb w6, [x0]
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
 staddb w6, [x0]
 staddb w6, [x1]
 .endr
apply_end:
 cmp x4, x8
 b.lo apply
 ret
//@end
//@seg x16_bfi_copy
x16_bfi_copy:
 .rept 64
 ldp x0, x1, [x15], #16
 stp x0, x1, [x25], #16
 .endr
 cmp x15, x16
 b.lo x16_bfi_copy
 ret
//@end
//@seg x16_bfi_pos_w8
x16_bfi_pos_w8:
 .rept 8
 ldp x0, x1, [x4], #16
 add x5, x22, x0, lsr #4
 bfi x5, x0, #46, #3
 staddb w23, [x5]
 add x5, x22, x1, lsr #4
 bfi x5, x1, #46, #3
 staddb w23, [x5]
 .endr
 cmp x4, x24
 b.lo x16_bfi_pos_w8
 ret
//@end
//@seg x16_bfi_pos_w16
x16_bfi_pos_w16:
 .rept 8
 ldp x0, x1, [x4], #16
 add x5, x22, x0, lsr #4
 bfi x5, x0, #46, #4
 staddb w23, [x5]
 add x5, x22, x1, lsr #4
 bfi x5, x1, #46, #4
 staddb w23, [x5]
 .endr
 cmp x4, x24
 b.lo x16_bfi_pos_w16
 ret
//@end
//@seg x16_bfi_bnd1_w16
x16_bfi_bnd1_w16:
 .rept 8
 ldp x0, x1, [x4], #16
 and w6, w0, #15
 cmp w6, w25
 b.lo 1f
 add x5, x22, x0, lsr #4
 bfi x5, x0, #46, #4
 staddb w23, [x5]
1:
 and w6, w1, #15
 cmp w6, w25
 b.lo 2f
 add x5, x22, x1, lsr #4
 bfi x5, x1, #46, #4
 staddb w23, [x5]
2:
 .endr
 cmp x4, x24
 b.lo x16_bfi_bnd1_w16
 ret
//@end
//@seg x16_bfi_bnd2_w16
x16_bfi_bnd2_w16:
 .rept 8
 ldp x0, x1, [x4], #16
 and w6, w0, #15
 sub w6, w6, w25
 cmp w6, w26
 b.hi 1f
 add x5, x22, x0, lsr #4
 bfi x5, x0, #46, #4
 staddb w23, [x5]
1:
 and w6, w1, #15
 sub w6, w6, w25
 cmp w6, w26
 b.hi 2f
 add x5, x22, x1, lsr #4
 bfi x5, x1, #46, #4
 staddb w23, [x5]
2:
 .endr
 cmp x4, x24
 b.lo x16_bfi_bnd2_w16
 ret
//@end
//@seg init8_u8
init8_u8:
 .rept 8
 ldp x0, x1, [x15], #16
 stp x0, x1, [x20], #16
 stp x0, x1, [x21], #16
 stp x0, x1, [x22], #16
 stp x0, x1, [x23], #16
 stp x0, x1, [x24], #16
 stp x0, x1, [x25], #16
 stp x0, x1, [x26], #16
 stp x0, x1, [x27], #16
 .endr
 cmp x15, x16
 b.lo init8_u8
 ret
//@end
//@seg scan8_u8
scan8_u8:
 mov x14, #0
 .rept 8
 ldp x0, x1, [x20], #16
 ldp x2, x3, [x21], #16
 orr x0, x0, x2
 orr x1, x1, x3
 ldp x2, x3, [x22], #16
 orr x0, x0, x2
 orr x1, x1, x3
 ldp x2, x3, [x23], #16
 orr x0, x0, x2
 orr x1, x1, x3
 ldp x2, x3, [x24], #16
 orr x0, x0, x2
 orr x1, x1, x3
 ldp x2, x3, [x25], #16
 orr x0, x0, x2
 orr x1, x1, x3
 ldp x2, x3, [x26], #16
 orr x0, x0, x2
 orr x1, x1, x3
 ldp x2, x3, [x27], #16
 orr x0, x0, x2
 orr x1, x1, x3
 orr x0, x0, x1
 orr x14, x14, x0
 .endr
 tst x14, x19
 b.ne .Lscan8_hit
 cmp x20, x28
 b.lo scan8_u8
 ret
.Lscan8_hit:
 ret
//@end
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
//@seg b1_page
b1_page:
 ldp x11, x6, [x28]
 add x4, x28, #16
 add x8, x4, x11, lsl #3
 tbz x11, #0, b1_page_even
 ldr x0, [x4], #8
 staddb w6, [x0]
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
 staddb w6, [x0]
 staddb w6, [x1]
 .endr
b1_apply_end:
 cmp x4, x8
 b.lo b1_apply
//@end

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
 staddb w6, [x5]
 add x5, x5, x0, lsr #32
 staddb w6, [x5]
 add x5, x5, w1, uxtw
 staddb w6, [x5]
 add x5, x5, x1, lsr #32
 staddb w6, [x5]
 .endr
delta_b1_apply_end:
 cmp x4, x8
 b.lo delta_b1_apply
 ret
//@end

//@seg mul4_sb
mul4_sb:
 ldp x3, x4, [x0]
 ldp x5, x6, [x0, #16]
 ldp x7, x8, [x1]
 ldp x9, x10, [x1, #16]
 mul x11, x3, x7
 umulh x12, x3, x7
 mov x13, #0
 str x11, [x2]
 mov x11, x12
 mov x12, #0
 mul x14, x3, x8
 umulh x15, x3, x8
 adds x11, x11, x14
 adcs x12, x12, x15
 adc x13, x13, xzr
 mul x14, x4, x7
 umulh x15, x4, x7
 adds x11, x11, x14
 adcs x12, x12, x15
 adc x13, x13, xzr
 str x11, [x2, #8]
 mov x11, x12
 mov x12, x13
 mov x13, #0
 mul x14, x3, x9
 umulh x15, x3, x9
 adds x11, x11, x14
 adcs x12, x12, x15
 adc x13, x13, xzr
 mul x14, x4, x8
 umulh x15, x4, x8
 adds x11, x11, x14
 adcs x12, x12, x15
 adc x13, x13, xzr
 mul x14, x5, x7
 umulh x15, x5, x7
 adds x11, x11, x14
 adcs x12, x12, x15
 adc x13, x13, xzr
 str x11, [x2, #16]
 mov x11, x12
 mov x12, x13
 mov x13, #0
 mul x14, x3, x10
 umulh x15, x3, x10
 adds x11, x11, x14
 adcs x12, x12, x15
 adc x13, x13, xzr
 mul x14, x4, x9
 umulh x15, x4, x9
 adds x11, x11, x14
 adcs x12, x12, x15
 adc x13, x13, xzr
 mul x14, x5, x8
 umulh x15, x5, x8
 adds x11, x11, x14
 adcs x12, x12, x15
 adc x13, x13, xzr
 mul x14, x6, x7
 umulh x15, x6, x7
 adds x11, x11, x14
 adcs x12, x12, x15
 adc x13, x13, xzr
 str x11, [x2, #24]
 mov x11, x12
 mov x12, x13
 mov x13, #0
 mul x14, x4, x10
 umulh x15, x4, x10
 adds x11, x11, x14
 adcs x12, x12, x15
 adc x13, x13, xzr
 mul x14, x5, x9
 umulh x15, x5, x9
 adds x11, x11, x14
 adcs x12, x12, x15
 adc x13, x13, xzr
 mul x14, x6, x8
 umulh x15, x6, x8
 adds x11, x11, x14
 adcs x12, x12, x15
 adc x13, x13, xzr
 str x11, [x2, #32]
 mov x11, x12
 mov x12, x13
 mov x13, #0
 mul x14, x5, x10
 umulh x15, x5, x10
 adds x11, x11, x14
 adcs x12, x12, x15
 adc x13, x13, xzr
 mul x14, x6, x9
 umulh x15, x6, x9
 adds x11, x11, x14
 adcs x12, x12, x15
 adc x13, x13, xzr
 str x11, [x2, #40]
 mov x11, x12
 mov x12, x13
 mul x14, x6, x10
 umulh x15, x6, x10
 adds x11, x11, x14
 adc x12, x12, x15
 stp x11, x12, [x2, #48]
 ret
//@end

//@seg srt_odd_norm8_find
srt_odd_norm8_find:
 ldr x3, [x0, x2, lsl #3]
 cbnz x3, srt_odd_norm8_mid
 add x2, x2, #1
 cmp x2, #8
 b.lo srt_odd_norm8_find
//@end

//@seg srt_odd_norm8_mid
srt_odd_norm8_mid:
 rbit x4, x3
 clz x4, x4
 mov x5, #64
 sub x5, x5, x4
 mov x6, #0
//@end

//@seg srt_odd_norm8_shift
srt_odd_norm8_shift:
 add x7, x6, x2
 cmp x7, #8
 b.hs srt_odd_zero_limb
 ldr x8, [x0, x7, lsl #3]
 lsr x8, x8, x4
 add x9, x7, #1
 cmp x9, #8
 b.hs srt_odd_no_hi
 ldr x10, [x0, x9, lsl #3]
 lsl x10, x10, x5
 tst x4, x4
 csel x10, x10, xzr, ne
 orr x8, x8, x10
srt_odd_no_hi:
 b srt_odd_store
srt_odd_zero_limb:
 mov x8, #0
srt_odd_store:
 str x8, [x1, x6, lsl #3]
 add x6, x6, #1
 cmp x6, #8
 b.lo srt_odd_norm8_shift
//@end

//@seg inv64_odd
inv64_odd:
 mov x1, #1
 mov x3, #2
 mul x2, x0, x1
 sub x2, x3, x2
 mul x1, x1, x2
 mul x2, x0, x1
 sub x2, x3, x2
 mul x1, x1, x2
 mul x2, x0, x1
 sub x2, x3, x2
 mul x1, x1, x2
 mul x2, x0, x1
 sub x2, x3, x2
 mul x1, x1, x2
 mul x2, x0, x1
 sub x2, x3, x2
 mul x1, x1, x2
 mul x2, x0, x1
 sub x2, x3, x2
 mul x1, x1, x2
 ret
//@end

//@seg cios8_setup
cios8_setup:
 stp xzr, xzr, [x4]
 stp xzr, xzr, [x4, #16]
 stp xzr, xzr, [x4, #32]
 stp xzr, xzr, [x4, #48]
 stp xzr, xzr, [x4, #64]
 mov x5, #0
//@end

//@seg cios8_outer_pre
cios8_outer_pre:
 ldr x6, [x0, x5, lsl #3]
 mov x7, #0
 mov x8, #0
//@end

//@seg cios8_inner_mul
cios8_inner_mul:
 ldr x9, [x1, x8, lsl #3]
 ldr x10, [x4, x8, lsl #3]
 mul x11, x6, x9
 umulh x12, x6, x9
 adds x11, x11, x10
 adc x12, x12, xzr
 adds x11, x11, x7
 adc x7, x12, xzr
 str x11, [x4, x8, lsl #3]
 add x8, x8, #1
 cmp x8, #8
 b.lo cios8_inner_mul
//@end

//@seg cios8_outer_mid
cios8_outer_mid:
 ldr x10, [x4, #64]
 adds x10, x10, x7
 adc x13, xzr, xzr
 str x10, [x4, #64]
 str x13, [x4, #72]
 ldr x14, [x4]
 mul x15, x14, x3
 ldr x9, [x2]
 mul x11, x15, x9
 umulh x12, x15, x9
 adds x11, x11, x14
 adc x7, x12, xzr
 mov x8, #1
//@end

//@seg cios8_inner_red
cios8_inner_red:
 ldr x9, [x2, x8, lsl #3]
 ldr x10, [x4, x8, lsl #3]
 mul x11, x15, x9
 umulh x12, x15, x9
 adds x11, x11, x10
 adc x12, x12, xzr
 adds x11, x11, x7
 adc x7, x12, xzr
 sub x16, x8, #1
 str x11, [x4, x16, lsl #3]
 add x8, x8, #1
 cmp x8, #8
 b.lo cios8_inner_red
//@end

//@seg cios8_outer_post
cios8_outer_post:
 ldr x10, [x4, #64]
 ldr x13, [x4, #72]
 adds x10, x10, x7
 adc x13, x13, xzr
 str x10, [x4, #56]
 str x13, [x4, #64]
 str xzr, [x4, #72]
 add x5, x5, #1
 cmp x5, #8
 b.lo cios8_outer_pre
//@end

//@seg cios8_csub_loop
cios8_csub_loop:
 ldr x9, [x4, x8, lsl #3]
 ldr x10, [x2, x8, lsl #3]
 sbcs x11, x9, x10
 str x11, [x17, x8, lsl #3]
 add x8, x8, #1
 sub x7, x7, #1
 cbnz x7, cios8_csub_loop
//@end

//@seg cios8_csel_loop
cios8_csel_loop:
 ldr x9, [x4, x8, lsl #3]
 ldr x10, [x17, x8, lsl #3]
 csel x11, x10, x9, hs
 str x11, [x17, x8, lsl #3]
 add x8, x8, #1
 sub x7, x7, #1
 cbnz x7, cios8_csel_loop
//@end

//@seg bingcd8_io
bingcd8_io:
 ldp x0, x1, [x25]
 ldp x2, x3, [x25, #16]
 ldp x4, x5, [x25, #32]
 ldp x6, x7, [x25, #48]
 ldp x8, x9, [x26]
 ldp x10, x11, [x26, #16]
 ldp x12, x13, [x26, #32]
 ldp x14, x15, [x26, #48]
 stp x0, x1, [x27]
 stp x2, x3, [x27, #16]
 stp x4, x5, [x27, #32]
 stp x6, x7, [x27, #48]
 ret
//@end

//@seg bingcd8_chk
bingcd8_loop:
 orr x16, x8, x9
 orr x17, x10, x11
 orr x16, x16, x17
 orr x17, x12, x13
 orr x16, x16, x17
 orr x17, x14, x15
 orr x16, x16, x17
 cbz x16, bingcd8_done
//@end

//@seg bingcd8_odd_sub
bingcd8_odd_sub:
 tbz x8, #0, bingcd8_shr_v
 subs x16, x8, x0
 sbcs x17, x9, x1
 sbcs x19, x10, x2
 sbcs x20, x11, x3
 sbcs x21, x12, x4
 sbcs x22, x13, x5
 sbcs x23, x14, x6
 sbcs x24, x15, x7
 b.hs bingcd8_shr_d
//@end

//@seg bingcd8_swap
bingcd8_swap:
 subs x16, x0, x8
 sbcs x17, x1, x9
 sbcs x19, x2, x10
 sbcs x20, x3, x11
 sbcs x21, x4, x12
 sbcs x22, x5, x13
 sbcs x23, x6, x14
 sbcs x24, x7, x15
 mov x0, x8
 mov x1, x9
 mov x2, x10
 mov x3, x11
 mov x4, x12
 mov x5, x13
 mov x6, x14
 mov x7, x15
//@end

//@seg bingcd8_shr_d
bingcd8_shr_d:
 extr x8, x17, x16, #1
 extr x9, x19, x17, #1
 extr x10, x20, x19, #1
 extr x11, x21, x20, #1
 extr x12, x22, x21, #1
 extr x13, x23, x22, #1
 extr x14, x24, x23, #1
 lsr x15, x24, #1
 b bingcd8_loop
//@end

//@seg bingcd8_shr_v
bingcd8_shr_v:
 extr x8, x9, x8, #1
 extr x9, x10, x9, #1
 extr x10, x11, x10, #1
 extr x11, x12, x11, #1
 extr x12, x13, x12, #1
 extr x13, x14, x13, #1
 extr x14, x15, x14, #1
 lsr x15, x15, #1
 b bingcd8_loop
bingcd8_done:
//@end

//@seg hensel_div8_outer
hensel_div8_outer:
 ldr x10, [x21, x8, lsl #3]
 mul x9, x10, x1
 str x9, [x23, x8, lsl #3]
 mov x11, #0
 mov x12, #0
 mov x10, #0
 mov x14, #8
 sub x14, x14, x8
//@end

//@seg hensel_div8_inner
hensel_div8_inner:
 ldr x13, [x22, x10, lsl #3]
 add x15, x8, x10
 ldr x16, [x21, x15, lsl #3]
 mul x17, x9, x13
 umulh x7, x9, x13
 adds x17, x17, x11
 adc x11, x7, xzr
 subs x16, x16, x17
 cset x7, lo
 subs x16, x16, x12
 cset x17, lo
 orr x12, x7, x17
 str x16, [x21, x15, lsl #3]
 add x10, x10, #1
 cmp x10, x14
 b.lo hensel_div8_inner
//@end

//@seg r4_4clamp_inner_tw
seg_r4_4clamp_inner_tw:
.Lr4_inner:
 ldr x5, [x0]
 ldr x6, [x1]
 ldr x7, [x2]
 ldr x8, [x3]
 ldp x9, x10, [x4], #16
 ldp x11, x12, [x4], #16
 ldp x13, x14, [x4], #16
 subs x15, x5, x20
 csel x5, x15, x5, hs
 umulh x15, x6, x10
 mul x6, x6, x9
 msub x6, x15, x19, x6
 umulh x15, x7, x12
 mul x7, x7, x11
 msub x7, x15, x19, x7
 umulh x15, x8, x14
 mul x8, x8, x13
 msub x8, x15, x19, x8
 add x9, x5, x7
 subs x15, x9, x20
 csel x9, x15, x9, hs
 add x10, x5, x20
 sub x10, x10, x7
 subs x15, x10, x20
 csel x10, x15, x10, hs
 add x11, x6, x8
 subs x15, x11, x20
 csel x11, x15, x11, hs
 add x12, x6, x20
 sub x12, x12, x8
 umulh x15, x12, x23
 mul x12, x12, x22
 msub x12, x15, x19, x12
 add x5, x9, x11
 add x7, x9, x20
 sub x7, x7, x11
 add x6, x10, x12
 add x8, x10, x20
 sub x8, x8, x12
 str x5, [x0], #8
 str x6, [x1], #8
 str x7, [x2], #8
 str x8, [x3], #8
 cmp x0, x24
 b.lo .Lr4_inner
//@end
 ret

//@seg r4_4clamp_hoist_tw_x2
seg_r4_4clamp_hoist_tw_x2:
.Lr4_hoist:
 ldp x5, x16, [x0]
 ldp x6, x17, [x1]
 ldp x7, x25, [x2]
 ldp x8, x26, [x3]
 subs x15, x5, x20
 csel x5, x15, x5, hs
 umulh x15, x6, x10
 mul x6, x6, x9
 msub x6, x15, x19, x6
 umulh x15, x7, x12
 mul x7, x7, x11
 msub x7, x15, x19, x7
 umulh x15, x8, x14
 mul x8, x8, x13
 msub x8, x15, x19, x8
 add x27, x5, x7
 subs x15, x27, x20
 csel x27, x15, x27, hs
 add x28, x5, x20
 sub x28, x28, x7
 subs x15, x28, x20
 csel x28, x15, x28, hs
 add x29, x6, x8
 subs x15, x29, x20
 csel x29, x15, x29, hs
 add x21, x6, x20
 sub x21, x21, x8
 umulh x15, x21, x23
 mul x21, x21, x22
 msub x21, x15, x19, x21
 add x5, x27, x29
 add x7, x27, x20
 sub x7, x7, x29
 add x6, x28, x21
 add x8, x28, x20
 sub x8, x8, x21
 subs x15, x16, x20
 csel x16, x15, x16, hs
 umulh x15, x17, x10
 mul x17, x17, x9
 msub x17, x15, x19, x17
 umulh x15, x25, x12
 mul x25, x25, x11
 msub x25, x15, x19, x25
 umulh x15, x26, x14
 mul x26, x26, x13
 msub x26, x15, x19, x26
 add x27, x16, x25
 subs x15, x27, x20
 csel x27, x15, x27, hs
 add x28, x16, x20
 sub x28, x28, x25
 subs x15, x28, x20
 csel x28, x15, x28, hs
 add x29, x17, x26
 subs x15, x29, x20
 csel x29, x15, x29, hs
 add x21, x17, x20
 sub x21, x21, x26
 umulh x15, x21, x23
 mul x21, x21, x22
 msub x21, x15, x19, x21
 add x16, x27, x29
 add x25, x27, x20
 sub x25, x25, x29
 add x17, x28, x21
 add x26, x28, x20
 sub x26, x26, x21
 stp x5, x16, [x0], #16
 stp x6, x17, [x1], #16
 stp x7, x25, [x2], #16
 stp x8, x26, [x3], #16
 cmp x0, x24
 b.lo .Lr4_hoist
//@end
 ret

//@seg r3_shoup_bfly
seg_r3_shoup_bfly:
.Lr3_loop:
 ldr x5, [x0]
 ldr x6, [x1]
 ldr x7, [x2]
 ldp x8, x9, [x3], #16
 ldp x10, x11, [x3], #16
 subs x15, x5, x20
 csel x5, x15, x5, hs
 umulh x15, x6, x9
 mul x6, x6, x8
 msub x6, x15, x19, x6
 umulh x15, x7, x11
 mul x7, x7, x10
 msub x7, x15, x19, x7
 add x12, x6, x7
 subs x15, x12, x20
 csel x12, x15, x12, hs
 add x13, x6, x20
 sub x13, x13, x7
 umulh x15, x12, x22
 mul x14, x12, x21
 msub x14, x15, x19, x14
 umulh x15, x13, x25
 mul x13, x13, x23
 msub x13, x15, x19, x13
 add x14, x5, x14
 subs x15, x14, x20
 csel x14, x15, x14, hs
 add x5, x5, x12
 add x6, x14, x13
 add x7, x14, x20
 sub x7, x7, x13
 str x5, [x0], #8
 str x6, [x1], #8
 str x7, [x2], #8
 cmp x0, x24
 b.lo .Lr3_loop
//@end
 ret

//@seg r5_shoup_bfly
seg_r5_shoup_bfly:
.Lr5_loop:
 ldr x6, [x0]
 ldr x7, [x1]
 ldr x8, [x2]
 ldr x9, [x3]
 ldr x10, [x4]
 ldp x11, x12, [x5], #16
 ldp x13, x14, [x5], #16
 subs x15, x6, x20
 csel x6, x15, x6, hs
 umulh x15, x7, x12
 mul x7, x7, x11
 msub x7, x15, x19, x7
 umulh x15, x8, x14
 mul x8, x8, x13
 msub x8, x15, x19, x8
 ldp x11, x12, [x5], #16
 ldp x13, x14, [x5], #16
 umulh x15, x9, x12
 mul x9, x9, x11
 msub x9, x15, x19, x9
 umulh x15, x10, x14
 mul x10, x10, x13
 msub x10, x15, x19, x10
 add x11, x7, x10
 subs x15, x11, x20
 csel x11, x15, x11, hs
 add x12, x8, x9
 subs x15, x12, x20
 csel x12, x15, x12, hs
 add x13, x7, x20
 sub x13, x13, x10
 subs x15, x13, x20
 csel x13, x15, x13, hs
 add x14, x8, x20
 sub x14, x14, x9
 subs x15, x14, x20
 csel x14, x15, x14, hs
 add x16, x11, x12
 subs x15, x16, x20
 csel x16, x15, x16, hs
 add x17, x11, x20
 sub x17, x17, x12
 add x22, x13, x14
 add x7, x6, x16
 str x7, [x0], #8
 ldp x25, x26, [x21, #0]
 umulh x15, x16, x26
 mul x16, x16, x25
 msub x16, x15, x19, x16
 ldp x25, x26, [x21, #16]
 umulh x15, x17, x26
 mul x17, x17, x25
 msub x17, x15, x19, x17
 add x16, x6, x16
 subs x15, x16, x20
 csel x16, x15, x16, hs
 add x27, x16, x17
 subs x15, x27, x20
 csel x27, x15, x27, hs
 add x28, x16, x20
 sub x28, x28, x17
 subs x15, x28, x20
 csel x28, x15, x28, hs
 ldp x25, x26, [x21, #32]
 umulh x15, x22, x26
 mul x22, x22, x25
 msub x22, x15, x19, x22
 ldp x25, x26, [x21, #48]
 umulh x15, x13, x26
 mul x13, x13, x25
 msub x13, x15, x19, x13
 ldp x25, x26, [x21, #64]
 umulh x15, x14, x26
 mul x14, x14, x25
 msub x14, x15, x19, x14
 add x29, x22, x13
 subs x15, x29, x20
 csel x29, x15, x29, hs
 add x23, x22, x14
 subs x15, x23, x20
 csel x23, x15, x23, hs
 add x7, x27, x29
 add x10, x27, x20
 sub x10, x10, x29
 add x8, x28, x23
 add x9, x28, x20
 sub x9, x9, x23
 str x7, [x1], #8
 str x8, [x2], #8
 str x9, [x3], #8
 str x10, [x4], #8
 cmp x0, x24
 b.lo .Lr5_loop
//@end
 ret

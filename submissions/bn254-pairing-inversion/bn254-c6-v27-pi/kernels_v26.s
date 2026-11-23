// kernels_v26.s - the AArch64 sequences priced in KERNELS.md section 1 for the v26 FP nodes and kernels (the chord split, inlined at .Lchord_split, is priced per chord in KERNELS.md 1.3 and not listed). GNU as (aarch64): 0 errors.
// Part 1: complete Fincke-Pohst node bodies (12 levels, 5-op LUT4 isqrt). audit_bn_cm.py counts the instruction lines in each
//   //@seg block and composes the node prices exactly as in v20_lat_b2.py FP_CL (k=0: k0_body+k0_back; k=1: top_2+body_2+k1_exit;
//   k>=2 at level m: top_{m+1}+body_{m+1}+top_m). Part 2: leaf emit (11), x8-unrolled apply (14 per 8 records), Lower Band 2
//   x4 loop (22 per 4 hits), a_0-chord kta pass 0 (20 per piece + 6 per entry), Band-1 uint32 delta hit.
// Node bodies of v20_lat_b2.py full_tree(); general registers only, x18 unused; [sp, #0] holds the exact-level mask.
// Registers: x0-x14 scratch (chord split may clobber them), x15-x17 chord-split bases (atab,btab,wptr;
// never written here), x19-x30 k=0 loop state (x19 rem0, x20 dr1, x21 A, x22 B, x23 cnt,
// x24-x29 cur (6 SWAR words), x30 inv_b2[0]); in levels >= 3 code x30 = frame base xb.
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
//@end
//@seg body_11
seg_body_11:
.Lent_11:
 ldr x3, [sp]
 tbnz x3, #11, .Lexact_11
 str x0, [sp, #88]
 sub x3, x2, x0, lsl #32
 cmp x3, #0
 cneg x3, x3, mi
 ldp x4, x5, [x30, #16]
 subs x3, x3, x4
 csel x3, x3, xzr, hi
 lsr x3, x3, #13
 mul x3, x3, x3
 umulh x3, x3, x5
 ldp x6, x7, [x30, #32]
 cmp x3, x6
 b.hi .Ltop_11
 sub x6, x6, x3
 ldp x8, x9, [x30, #48]
 msub x2, x8, x0, x7
 ldr x10, [x30, #176]
 mneg x7, x10, x0
 stp x6, x7, [x30, #-224]
 ldp x10, x11, [x30, #80]
 ldp x12, x13, [x30, #96]
 ldp x14, x19, [x30, #112]
 ldp x20, x21, [x30, #128]
 ldp x22, x23, [x30, #144]
 ldp x24, x25, [x30, #160]
 madd x10, x0, x20, x10
 madd x11, x0, x21, x11
 madd x12, x0, x22, x12
 madd x13, x0, x23, x13
 madd x14, x0, x24, x14
 madd x19, x0, x25, x19
 stp x10, x11, [x30, #-176]
 stp x12, x13, [x30, #-160]
 stp x14, x19, [x30, #-144]
 umulh x10, x6, x9
 clz x12, x10
 lsl x11, x10, x12
 extr x11, x12, x11, #53
 ldr x11, [sp, x11, lsl #3]
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
//@end
//@seg body_10
seg_body_10:
.Lent_10:
 ldr x3, [sp]
 tbnz x3, #10, .Lexact_10
 str x0, [sp, #80]
 sub x3, x2, x0, lsl #32
 cmp x3, #0
 cneg x3, x3, mi
 ldp x4, x5, [x30, #16]
 subs x3, x3, x4
 csel x3, x3, xzr, hi
 lsr x3, x3, #13
 mul x3, x3, x3
 umulh x3, x3, x5
 ldp x6, x7, [x30, #32]
 cmp x3, x6
 b.hi .Ltop_10
 sub x6, x6, x3
 ldp x8, x9, [x30, #48]
 msub x2, x8, x0, x7
 ldr x10, [sp, #88]
 ldp x11, x12, [x30, #176]
 mneg x7, x11, x0
 msub x7, x12, x10, x7
 stp x6, x7, [x30, #-224]
 ldp x10, x11, [x30, #80]
 ldp x12, x13, [x30, #96]
 ldp x14, x19, [x30, #112]
 ldp x20, x21, [x30, #128]
 ldp x22, x23, [x30, #144]
 ldp x24, x25, [x30, #160]
 madd x10, x0, x20, x10
 madd x11, x0, x21, x11
 madd x12, x0, x22, x12
 madd x13, x0, x23, x13
 madd x14, x0, x24, x14
 madd x19, x0, x25, x19
 stp x10, x11, [x30, #-176]
 stp x12, x13, [x30, #-160]
 stp x14, x19, [x30, #-144]
 umulh x10, x6, x9
 clz x12, x10
 lsl x11, x10, x12
 extr x11, x12, x11, #53
 ldr x11, [sp, x11, lsl #3]
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
//@end
//@seg body_9
seg_body_9:
.Lent_9:
 ldr x3, [sp]
 tbnz x3, #9, .Lexact_9
 str x0, [sp, #72]
 sub x3, x2, x0, lsl #32
 cmp x3, #0
 cneg x3, x3, mi
 ldp x4, x5, [x30, #16]
 subs x3, x3, x4
 csel x3, x3, xzr, hi
 lsr x3, x3, #13
 mul x3, x3, x3
 umulh x3, x3, x5
 ldp x6, x7, [x30, #32]
 cmp x3, x6
 b.hi .Ltop_9
 sub x6, x6, x3
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
 ldp x12, x13, [x30, #96]
 ldp x14, x19, [x30, #112]
 ldp x20, x21, [x30, #128]
 ldp x22, x23, [x30, #144]
 ldp x24, x25, [x30, #160]
 madd x10, x0, x20, x10
 madd x11, x0, x21, x11
 madd x12, x0, x22, x12
 madd x13, x0, x23, x13
 madd x14, x0, x24, x14
 madd x19, x0, x25, x19
 stp x10, x11, [x30, #-176]
 stp x12, x13, [x30, #-160]
 stp x14, x19, [x30, #-144]
 umulh x10, x6, x9
 clz x12, x10
 lsl x11, x10, x12
 extr x11, x12, x11, #53
 ldr x11, [sp, x11, lsl #3]
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
//@end
//@seg body_8
seg_body_8:
.Lent_8:
 ldr x3, [sp]
 tbnz x3, #8, .Lexact_8
 str x0, [sp, #64]
 sub x3, x2, x0, lsl #32
 cmp x3, #0
 cneg x3, x3, mi
 ldp x4, x5, [x30, #16]
 subs x3, x3, x4
 csel x3, x3, xzr, hi
 lsr x3, x3, #13
 mul x3, x3, x3
 umulh x3, x3, x5
 ldp x6, x7, [x30, #32]
 cmp x3, x6
 b.hi .Ltop_8
 sub x6, x6, x3
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
 ldp x12, x13, [x30, #96]
 ldp x14, x19, [x30, #112]
 ldp x20, x21, [x30, #128]
 ldp x22, x23, [x30, #144]
 ldp x24, x25, [x30, #160]
 madd x10, x0, x20, x10
 madd x11, x0, x21, x11
 madd x12, x0, x22, x12
 madd x13, x0, x23, x13
 madd x14, x0, x24, x14
 madd x19, x0, x25, x19
 stp x10, x11, [x30, #-176]
 stp x12, x13, [x30, #-160]
 stp x14, x19, [x30, #-144]
 umulh x10, x6, x9
 clz x12, x10
 lsl x11, x10, x12
 extr x11, x12, x11, #53
 ldr x11, [sp, x11, lsl #3]
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
//@end
//@seg body_7
seg_body_7:
.Lent_7:
 ldr x3, [sp]
 tbnz x3, #7, .Lexact_7
 str x0, [sp, #56]
 sub x3, x2, x0, lsl #32
 cmp x3, #0
 cneg x3, x3, mi
 ldp x4, x5, [x30, #16]
 subs x3, x3, x4
 csel x3, x3, xzr, hi
 lsr x3, x3, #13
 mul x3, x3, x3
 umulh x3, x3, x5
 ldp x6, x7, [x30, #32]
 cmp x3, x6
 b.hi .Ltop_7
 sub x6, x6, x3
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
 ldp x12, x13, [x30, #96]
 ldp x14, x19, [x30, #112]
 ldp x20, x21, [x30, #128]
 ldp x22, x23, [x30, #144]
 ldp x24, x25, [x30, #160]
 madd x10, x0, x20, x10
 madd x11, x0, x21, x11
 madd x12, x0, x22, x12
 madd x13, x0, x23, x13
 madd x14, x0, x24, x14
 madd x19, x0, x25, x19
 stp x10, x11, [x30, #-176]
 stp x12, x13, [x30, #-160]
 stp x14, x19, [x30, #-144]
 umulh x10, x6, x9
 clz x12, x10
 lsl x11, x10, x12
 extr x11, x12, x11, #53
 ldr x11, [sp, x11, lsl #3]
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
//@end
//@seg body_6
seg_body_6:
.Lent_6:
 ldr x3, [sp]
 tbnz x3, #6, .Lexact_6
 str x0, [sp, #48]
 sub x3, x2, x0, lsl #32
 cmp x3, #0
 cneg x3, x3, mi
 ldp x4, x5, [x30, #16]
 subs x3, x3, x4
 csel x3, x3, xzr, hi
 lsr x3, x3, #13
 mul x3, x3, x3
 umulh x3, x3, x5
 ldp x6, x7, [x30, #32]
 cmp x3, x6
 b.hi .Ltop_6
 sub x6, x6, x3
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
 ldp x12, x13, [x30, #96]
 ldp x14, x19, [x30, #112]
 ldp x20, x21, [x30, #128]
 ldp x22, x23, [x30, #144]
 ldp x24, x25, [x30, #160]
 madd x10, x0, x20, x10
 madd x11, x0, x21, x11
 madd x12, x0, x22, x12
 madd x13, x0, x23, x13
 madd x14, x0, x24, x14
 madd x19, x0, x25, x19
 stp x10, x11, [x30, #-176]
 stp x12, x13, [x30, #-160]
 stp x14, x19, [x30, #-144]
 umulh x10, x6, x9
 clz x12, x10
 lsl x11, x10, x12
 extr x11, x12, x11, #53
 ldr x11, [sp, x11, lsl #3]
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
//@end
//@seg body_5
seg_body_5:
.Lent_5:
 ldr x3, [sp]
 tbnz x3, #5, .Lexact_5
 str x0, [sp, #40]
 sub x3, x2, x0, lsl #32
 cmp x3, #0
 cneg x3, x3, mi
 ldp x4, x5, [x30, #16]
 subs x3, x3, x4
 csel x3, x3, xzr, hi
 lsr x3, x3, #13
 mul x3, x3, x3
 umulh x3, x3, x5
 ldp x6, x7, [x30, #32]
 cmp x3, x6
 b.hi .Ltop_5
 sub x6, x6, x3
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
 ldp x12, x13, [x30, #96]
 ldp x14, x19, [x30, #112]
 ldp x20, x21, [x30, #128]
 ldp x22, x23, [x30, #144]
 ldp x24, x25, [x30, #160]
 madd x10, x0, x20, x10
 madd x11, x0, x21, x11
 madd x12, x0, x22, x12
 madd x13, x0, x23, x13
 madd x14, x0, x24, x14
 madd x19, x0, x25, x19
 stp x10, x11, [x30, #-176]
 stp x12, x13, [x30, #-160]
 stp x14, x19, [x30, #-144]
 umulh x10, x6, x9
 clz x12, x10
 lsl x11, x10, x12
 extr x11, x12, x11, #53
 ldr x11, [sp, x11, lsl #3]
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
//@end
//@seg body_4
seg_body_4:
.Lent_4:
 ldr x3, [sp]
 tbnz x3, #4, .Lexact_4
 str x0, [sp, #32]
 sub x3, x2, x0, lsl #32
 cmp x3, #0
 cneg x3, x3, mi
 ldp x4, x5, [x30, #16]
 subs x3, x3, x4
 csel x3, x3, xzr, hi
 lsr x3, x3, #13
 mul x3, x3, x3
 umulh x3, x3, x5
 ldp x6, x7, [x30, #32]
 cmp x3, x6
 b.hi .Ltop_4
 sub x6, x6, x3
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
 ldp x12, x13, [x30, #96]
 ldp x14, x19, [x30, #112]
 ldp x20, x21, [x30, #128]
 ldp x22, x23, [x30, #144]
 ldp x24, x25, [x30, #160]
 madd x10, x0, x20, x10
 madd x11, x0, x21, x11
 madd x12, x0, x22, x12
 madd x13, x0, x23, x13
 madd x14, x0, x24, x14
 madd x19, x0, x25, x19
 stp x10, x11, [x30, #-240]
 stp x12, x13, [x30, #-224]
 stp x14, x19, [x30, #-208]
 umulh x10, x6, x9
 clz x12, x10
 lsl x11, x10, x12
 extr x11, x12, x11, #53
 ldr x11, [sp, x11, lsl #3]
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
//@end
//@seg body_3
seg_body_3:
.Lent_3:
 ldr x3, [sp]
 tbnz x3, #3, .Lexact_3
 str x0, [sp, #24]
 sub x3, x2, x0, lsl #32
 cmp x3, #0
 cneg x3, x3, mi
 ldp x4, x5, [x30, #16]
 subs x3, x3, x4
 csel x3, x3, xzr, hi
 lsr x3, x3, #13
 mul x3, x3, x3
 umulh x3, x3, x5
 ldp x6, x7, [x30, #32]
 cmp x3, x6
 b.hi .Ltop_3
 sub x6, x6, x3
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
 ldp x12, x13, [x30, #96]
 ldp x14, x19, [x30, #112]
 ldp x20, x21, [x30, #128]
 ldp x22, x23, [x30, #144]
 ldp x24, x25, [x30, #160]
 madd x10, x0, x20, x10
 madd x11, x0, x21, x11
 madd x12, x0, x22, x12
 madd x13, x0, x23, x13
 madd x14, x0, x24, x14
 madd x19, x0, x25, x19
 stp x10, x11, [sp, #384]
 stp x12, x13, [sp, #400]
 stp x14, x19, [sp, #416]
 umulh x10, x6, x9
 clz x12, x10
 lsl x11, x10, x12
 extr x11, x12, x11, #53
 ldr x11, [sp, x11, lsl #3]
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
//@end
//@seg body_2
seg_body_2:
.Lent_2:
 str x0, [sp, #16]
 sub x3, x2, x0, lsl #32
 cmp x3, #0
 cneg x3, x3, mi
 ldp x4, x5, [sp, #320]
 subs x3, x3, x4
 csel x3, x3, xzr, hi
 lsr x3, x3, #13
 mul x3, x3, x3
 umulh x3, x3, x5
 ldp x6, x7, [sp, #336]
 cmp x3, x6
 b.hi .Ltop_2
 sub x6, x6, x3
 ldp x8, x9, [sp, #352]
 msub x2, x8, x0, x7
 ldp x8, x10, [sp, #368]
 msub x10, x8, x0, x10
 umulh x11, x6, x9
 clz x13, x11
 lsl x12, x11, x13
 extr x12, x13, x12, #53
 ldr x12, [sp, x12, lsl #3]
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
 sub x3, x2, x13, lsl #32
 cmp x3, #0
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
 ldp x8, x9, [sp, #96]
 msub x10, x9, x13, x10
 ldp x8, x9, [sp, #256]
 add x21, x10, x8
 add x22, x10, x9
 ldp x24, x25, [sp, #384]
 ldp x26, x27, [sp, #400]
 ldp x28, x29, [sp, #416]
 ldp x3, x4, [sp, #432]
 ldp x5, x6, [sp, #448]
 ldp x7, x8, [sp, #464]
 madd x24, x0, x3, x24
 madd x25, x0, x4, x25
 madd x26, x0, x5, x26
 madd x27, x0, x6, x27
 madd x28, x0, x7, x28
 madd x29, x0, x8, x29
 ldp x3, x4, [sp, #160]
 ldp x5, x6, [sp, #176]
 ldp x7, x8, [sp, #192]
 madd x24, x13, x3, x24
 madd x25, x13, x4, x25
 madd x26, x13, x5, x26
 madd x27, x13, x6, x27
 madd x28, x13, x7, x28
 madd x29, x13, x8, x29
//@end
//@seg k0_body
seg_k0_body:
.Lk0_top:
 tbnz x19, #63, .Lk0_skip
 umulh x8, x19, x30
 clz x9, x8
 lsl x8, x8, x9
 extr x8, x9, x8, #53
 ldr x8, [sp, x8, lsl #3]
 sub x6, x21, x8
 asr x6, x6, #32
 add x7, x22, x8
 asr x7, x7, #32
 cmp x6, x7
 b.gt .Lk0_skip
 ldp x0, x1, [sp, #112]
 ldp x2, x3, [sp, #128]
 ldp x4, x5, [sp, #144]
 madd x0, x6, x0, x24
 madd x1, x6, x1, x25
 madd x2, x6, x2, x26
 madd x3, x6, x3, x27
 madd x4, x6, x4, x28
 madd x5, x6, x5, x29
 ldp x8, x9, [sp, #96]
 sub x19, x19, x20
 add x20, x20, x8
 sub x21, x21, x9
 sub x22, x22, x9
 ldp x8, x9, [sp, #160]
 ldp x10, x11, [sp, #176]
 ldp x12, x13, [sp, #192]
 add x24, x24, x8
 add x25, x25, x9
 add x26, x26, x10
 add x27, x27, x11
 add x28, x28, x12
 add x29, x29, x13
 sub x23, x23, #1
//@end
.ifdef FP_TEST
 stp x0, x1, [x15], #16
 stp x2, x3, [x15], #16
 stp x4, x5, [x15], #16
 stp x6, x7, [x15], #16
.endif
.Lchord_split:
//@seg k0_back
seg_k0_back:
 tbz x23, #63, .Lk0_top
//@end
//@seg k1_exit
seg_k1_exit:
 b .Ltop_2
//@end
//@seg k0_skip
seg_k0_skip:
.Lk0_skip:
 ldp x8, x9, [sp, #96]
 sub x19, x19, x20
 add x20, x20, x8
 sub x21, x21, x9
 sub x22, x22, x9
 ldp x8, x9, [sp, #160]
 ldp x10, x11, [sp, #176]
 ldp x12, x13, [sp, #192]
 add x24, x24, x8
 add x25, x25, x9
 add x26, x26, x10
 add x27, x27, x11
 add x28, x28, x12
 add x29, x29, x13
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
emit:
 ldr x0, [x1, x2, lsl #3]
 ldr x3, [x4, x5, lsl #3]
 lsr x6, x3, #36
 add x7, x0, x6, lsl #35
 and x8, x3, #0xFFFFFFFFF
 ldr x9, [x10, x8, lsl #3]
 str x7, [x9], #8
 tbnz x9, #14, 2f
 str x9, [x10, x8, lsl #3]
 add x2, x2, x11
 add x5, x5, x12
2:
apply:
 ldp x0, x1, [x4], #16
 ldp x2, x3, [x4], #16
 ldp x8, x9, [x4], #16
 ldp x10, x11, [x4], #16
 staddh w6, [x0]
 staddh w6, [x1]
 staddh w6, [x2]
 staddh w6, [x3]
 staddh w6, [x8]
 staddh w6, [x9]
 staddh w6, [x10]
 staddh w6, [x11]
 subs x7, x7, #8
 b.gt apply
apply_tail6:
 ldp x0, x1, [x4], #16
 ldp x2, x3, [x4], #16
 ldp x8, x9, [x4], #16
 staddh w6, [x0]
 staddh w6, [x1]
 staddh w6, [x2]
 staddh w6, [x3]
 staddh w6, [x8]
 staddh w6, [x9]
x4_loop:
 ldp x0, x1, [x4], #16
 ldp x2, x3, [x4], #16
 cmp x3, x20
 b.hi x4_tail
 and x5, x0, #0xF
 madd x5, x5, x21, x22
 add x5, x5, x0, lsr #16
 staddh w23, [x5]
 and x5, x1, #0xF
 madd x5, x5, x21, x22
 add x5, x5, x1, lsr #16
 staddh w23, [x5]
 and x5, x2, #0xF
 madd x5, x5, x21, x22
 add x5, x5, x2, lsr #16
 staddh w23, [x5]
 and x5, x3, #0xF
 madd x5, x5, x21, x22
 add x5, x5, x3, lsr #16
 staddh w23, [x5]
 cmp x4, x24
 b.lo x4_loop
x4_tail:
// last 1..3 entries of a group: same 4-op hit body per entry, charged in C_SB_ENTER_X4 (not listed)
d32:
 ldp w0, w1, [x4], #8
 add x5, x5, x0, lsl #1
 staddh w6, [x5]

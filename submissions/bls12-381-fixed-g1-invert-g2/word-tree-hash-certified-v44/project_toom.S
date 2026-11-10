// Annotation-only revision: all arithmetic instruction/macro tokens unchanged.
// Inert analytic evidence. NEVER host-execute or use credentials.
// out576(a0@0,b1@384 only),A1040,B1040,p48,n0 in x0..4.
// A/B contain21 full-width native evaluations,then scalar32 (not read here).
// All actual scalar arithmetic uses the exact341-instruction pinned blst core.
.text
.p2align 2
.macro FMUL ao,bo,po
    ldp x1,x2,[x29,#136]
    add x1,x1,#\ao
    add x2,x2,#\bo
    add x0,x29,#(192+\po)
    stp x4,x0,[x29,#96]
    ldp x11,x12,[x1]
    ldr x17,[x2]
    ldp x13,x14,[x1,#16]
    ldp x15,x16,[x1,#32]
    bl __mul_mont_384
    stp x11,x12,[x2]
    stp x13,x14,[x2,#16]
    stp x15,x16,[x2,#32]
.endm
.macro FADD ob,oo,ab,ao,bb,bo
    add x0,\ob,#\oo
    add x1,\ab,#\ao
    add x2,\bb,#\bo
    mov x3,x22
    bl fp381_add
.endm
.macro FSUB ob,oo,ab,ao,bb,bo
    add x0,\ob,#\oo
    add x1,\ab,#\ao
    add x2,\bb,#\bo
    mov x3,x22
    bl fp381_sub
.endm
.macro FCOPY oo,sb,so
    add x23,\sb,#\so
    ldp x0,x1,[x23]
    stp x0,x1,[x19,#\oo]
    ldp x0,x1,[x23,#16]
    stp x0,x1,[x19,#(\oo+16)]
    ldp x0,x1,[x23,#32]
    stp x0,x1,[x19,#(\oo+32)]
.endm
.global project_toom
.type project_toom,%function
project_toom:
    sub sp,sp,#1,lsl #12
    sub sp,sp,#32
    stp x29,x30,[sp]
    mov x29,sp
    stp x19,x20,[x29,#16]
    stp x21,x22,[x29,#32]
    stp x23,x24,[x29,#48]
    stp x25,x26,[x29,#64]
    stp x27,x28,[x29,#80]
    stp x0,x1,[x29,#128]
    stp x2,x3,[x29,#144]
    stp x4,xzr,[x29,#160]
    ldp x3,x4,[x29,#152]
    ldp x5,x6,[x3]
    ldp x7,x8,[x3,#16]
    ldp x9,x10,[x3,#32]
    FMUL 0,0,0
    FMUL 48,48,48
    FMUL 96,96,96
    FMUL 144,144,144
    FMUL 192,192,192
    FMUL 240,240,240
    FMUL 288,288,288
    FMUL 336,672,336
    FMUL 384,720,384
    FMUL 432,768,432
    FMUL 480,816,480
    FMUL 528,864,528
    FMUL 576,912,576
    FMUL 624,960,624
    FMUL 672,336,672
    FMUL 720,384,720
    FMUL 768,432,768
    FMUL 816,480,816
    FMUL 864,528,864
    FMUL 912,576,912
    FMUL 960,624,960
    ldr x19,[x29,#128]
    add x20,x29,#192
    add x21,x29,#1200
    ldr x22,[x29,#152]
    FADD x20,336,x20,336,x20,672
    FADD x20,384,x20,384,x20,720
    FADD x20,432,x20,432,x20,768
    FADD x20,480,x20,480,x20,816
    FADD x20,528,x20,528,x20,864
    FADD x20,576,x20,576,x20,912
    FADD x20,624,x20,624,x20,960
    // Static coefficient certificate: temp0@0 = 2*P0.
    FADD x21,0,x20,0,x20,0
    // Static coefficient certificate: temp1@48 = 3*P0.
    FADD x21,48,x21,0,x20,0
    // Static coefficient certificate: temp2@96 = 4*P0.
    FADD x21,96,x21,0,x21,0
    // Static coefficient certificate: temp3@144 = 8*P0.
    FADD x21,144,x21,96,x21,96
    // Static coefficient certificate: temp4@192 = 9*P0.
    FADD x21,192,x21,144,x20,0
    // Static coefficient certificate: temp5@240 = 2*P1.
    FADD x21,240,x20,48,x20,48
    // Static coefficient certificate: temp6@288 = 4*P1.
    FADD x21,288,x21,240,x21,240
    // Static coefficient certificate: temp7@336 = 8*P1.
    FADD x21,336,x21,288,x21,288
    // Static coefficient certificate: temp8@384 = 5*P1.
    FADD x21,384,x21,288,x20,48
    // Static coefficient certificate: temp9@432 = 2*P2.
    FADD x21,432,x20,96,x20,96
    // Static coefficient certificate: temp10@480 = 4*P2.
    FADD x21,480,x21,432,x21,432
    // Static coefficient certificate: temp11@528 = 8*P2.
    FADD x21,528,x21,480,x21,480
    // Static coefficient certificate: temp12@576 = 16*P2.
    FADD x21,576,x21,528,x21,528
    // Static coefficient certificate: temp13@624 = 3*P2.
    FADD x21,624,x21,432,x20,96
    // Static coefficient certificate: temp14@672 = 6*P2.
    FADD x21,672,x21,624,x21,624
    // Static coefficient certificate: temp15@720 = 12*P2.
    FADD x21,720,x21,672,x21,672
    // Static coefficient certificate: temp16@768 = 13*P2.
    FADD x21,768,x21,720,x20,96
    // Static coefficient certificate: temp17@816 = 2*P3.
    FADD x21,816,x20,144,x20,144
    // Static coefficient certificate: temp18@864 = 4*P3.
    FADD x21,864,x21,816,x21,816
    // Static coefficient certificate: temp19@912 = 5*P3.
    FADD x21,912,x21,864,x20,144
    // Static coefficient certificate: temp20@960 = 2*P4.
    FADD x21,960,x20,192,x20,192
    // Static coefficient certificate: temp21@1008 = 4*P4.
    FADD x21,1008,x21,960,x21,960
    // Static coefficient certificate: temp22@1056 = 5*P4.
    FADD x21,1056,x21,1008,x20,192
    // Static coefficient certificate: temp23@1104 = 3*P4.
    FADD x21,1104,x21,960,x20,192
    // Static coefficient certificate: temp24@1152 = 6*P4.
    FADD x21,1152,x21,1104,x21,1104
    // Static coefficient certificate: temp25@1200 = 7*P4.
    FADD x21,1200,x21,1152,x20,192
    // Static coefficient certificate: temp26@1248 = 2*P5.
    FADD x21,1248,x20,240,x20,240
    // Static coefficient certificate: temp27@1296 = 3*P5.
    FADD x21,1296,x21,1248,x20,240
    // Static coefficient certificate: temp28@1344 = 2*P6.
    FADD x21,1344,x20,288,x20,288
    // Static coefficient certificate: temp29@1392 = 4*P6.
    FADD x21,1392,x21,1344,x21,1344
    // Static coefficient certificate: temp30@1440 = 3*P6.
    FADD x21,1440,x21,1344,x20,288
    // Static coefficient certificate: temp31@1488 = 8*P6.
    FADD x21,1488,x21,1392,x21,1392
    // Static coefficient certificate: temp32@1536 = 9*P6.
    FADD x21,1536,x21,1488,x20,288
    // Static coefficient certificate: temp33@1584 = 2*Q0.
    FADD x21,1584,x20,336,x20,336
    // Static coefficient certificate: temp34@1632 = 4*Q0.
    FADD x21,1632,x21,1584,x21,1584
    // Static coefficient certificate: temp35@1680 = 3*Q0.
    FADD x21,1680,x21,1584,x20,336
    // Static coefficient certificate: temp36@1728 = 8*Q0.
    FADD x21,1728,x21,1632,x21,1632
    // Static coefficient certificate: temp37@1776 = 9*Q0.
    FADD x21,1776,x21,1728,x20,336
    // Static coefficient certificate: temp38@1824 = 2*Q1.
    FADD x21,1824,x20,384,x20,384
    // Static coefficient certificate: temp39@1872 = 4*Q1.
    FADD x21,1872,x21,1824,x21,1824
    // Static coefficient certificate: temp40@1920 = 5*Q1.
    FADD x21,1920,x21,1872,x20,384
    // Static coefficient certificate: temp41@1968 = 2*Q2.
    FADD x21,1968,x20,432,x20,432
    // Static coefficient certificate: temp42@2016 = 3*Q2.
    FADD x21,2016,x21,1968,x20,432
    // Static coefficient certificate: temp43@2064 = 6*Q2.
    FADD x21,2064,x21,2016,x21,2016
    // Static coefficient certificate: temp44@2112 = 7*Q2.
    FADD x21,2112,x21,2064,x20,432
    // Static coefficient certificate: temp45@2160 = 4*Q2.
    FADD x21,2160,x21,1968,x21,1968
    // Static coefficient certificate: temp46@2208 = 8*Q2.
    FADD x21,2208,x21,2160,x21,2160
    // Static coefficient certificate: temp47@2256 = 2*Q3.
    FADD x21,2256,x20,480,x20,480
    // Static coefficient certificate: temp48@2304 = 4*Q3.
    FADD x21,2304,x21,2256,x21,2256
    // Static coefficient certificate: temp49@2352 = 3*Q3.
    FADD x21,2352,x21,2256,x20,480
    // Static coefficient certificate: temp50@2400 = 6*Q3.
    FADD x21,2400,x21,2352,x21,2352
    // Static coefficient certificate: temp51@2448 = 2*Q4.
    FADD x21,2448,x20,528,x20,528
    // Static coefficient certificate: temp52@2496 = 4*Q4.
    FADD x21,2496,x21,2448,x21,2448
    // Static coefficient certificate: temp53@2544 = 8*Q4.
    FADD x21,2544,x21,2496,x21,2496
    // Static coefficient certificate: temp54@2592 = 3*Q4.
    FADD x21,2592,x21,2448,x20,528
    // Static coefficient certificate: temp55@2640 = 6*Q4.
    FADD x21,2640,x21,2592,x21,2592
    // Static coefficient certificate: temp56@2688 = 5*Q4.
    FADD x21,2688,x21,2496,x20,528
    // Static coefficient certificate: temp57@2736 = 2*Q5.
    FADD x21,2736,x20,576,x20,576
    // Static coefficient certificate: temp58@2784 = 2*Q6.
    FADD x21,2784,x20,624,x20,624
    // Static coefficient certificate: temp59@2832 = 3*Q6.
    FADD x21,2832,x21,2784,x20,624
    // Static coefficient certificate: temp60@2880 = 6*Q6.
    FADD x21,2880,x21,2832,x21,2832
    FCOPY 0,x21,48
    FADD x19,0,x19,0,x21,336
    FADD x19,0,x19,0,x21,576
    FADD x19,0,x19,0,x20,144
    FADD x19,0,x19,0,x21,1056
    FADD x19,0,x19,0,x21,1344
    FADD x19,0,x19,0,x21,1584
    FADD x19,0,x19,0,x21,1920
    FSUB x19,0,x19,0,x21,2112
    FADD x19,0,x19,0,x21,2304
    FSUB x19,0,x19,0,x21,2544
    FADD x19,0,x19,0,x20,576
    FSUB x19,0,x19,0,x21,2832
    FCOPY 48,x21,192
    FADD x19,48,x19,48,x21,288
    FADD x19,48,x19,48,x21,528
    FSUB x19,48,x19,48,x20,144
    FSUB x19,48,x19,48,x21,1056
    FSUB x19,48,x19,48,x21,1392
    FSUB x19,48,x19,48,x21,1632
    FSUB x19,48,x19,48,x21,2064
    FSUB x19,48,x19,48,x21,2400
    FADD x19,48,x19,48,x21,2640
    FSUB x19,48,x19,48,x21,2736
    FADD x19,48,x19,48,x21,2880
    FCOPY 384,x21,384
    FSUB x19,384,x19,384,x20,0
    FSUB x19,384,x19,384,x21,768
    FSUB x19,384,x19,384,x20,144
    FSUB x19,384,x19,384,x20,192
    FSUB x19,384,x19,384,x20,240
    FADD x19,384,x19,384,x21,1440
    FADD x19,384,x19,384,x21,1680
    FADD x19,384,x19,384,x21,1872
    FADD x19,384,x19,384,x21,2208
    FADD x19,384,x19,384,x20,480
    FADD x19,384,x19,384,x21,2688
    FADD x19,384,x19,384,x20,624
    FCOPY 432,x20,96
    FSUB x19,432,x19,432,x21,48
    FSUB x19,432,x19,432,x21,384
    FSUB x19,432,x19,432,x21,912
    FADD x19,432,x19,432,x21,1200
    FSUB x19,432,x19,432,x21,1296
    FADD x19,432,x19,432,x21,1536
    FADD x19,432,x19,432,x21,1776
    FADD x19,432,x19,432,x21,1824
    FADD x19,432,x19,432,x21,2160
    FSUB x19,432,x19,432,x20,480
    FSUB x19,432,x19,432,x21,2688
    FSUB x19,432,x19,432,x21,2784
    ldp x19,x20,[x29,#16]
    ldp x21,x22,[x29,#32]
    ldp x23,x24,[x29,#48]
    ldp x25,x26,[x29,#64]
    ldp x27,x28,[x29,#80]
    ldp x29,x30,[sp]
    add sp,sp,#1,lsl #12
    add sp,sp,#32
    ret
.size project_toom,.-project_toom

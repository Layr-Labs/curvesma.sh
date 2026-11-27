.arch armv8.1-a
.text
//@seg top_11
seg_top_11:
.Ltop_11:
add x30,sp,#2624
ldr x0,[sp,#88]
ldp x1,x2,[x30,#0]
add x0,x0,#1
cmp x0,x1
b.gt .Lorbit_done
ldrh w9,[sp,#3862]
add x10,sp,#22,lsl #12
ldrh w9,[x10,x9,lsl #1]
strh w9,[sp,#3862]
//@end
//@seg body_11
seg_body_11:
.Lent_11:
ldr x3,[sp]
tbnz x3,#11,.Lexact_11
str x0,[sp,#88]
subs x3,x2,x0,lsl #32
cneg x3,x3,mi
ldp x4,x5,[x30,#16]
subs x3,x3,x4
csel x3,x3,xzr,hi
lsr x3,x3,#13
mul x3,x3,x3
umulh x3,x3,x5
ldp x6,x7,[x30,#32]
subs x6,x6,x3
b.lo .Ltop_11
ldp x8,x9,[x30,#48]
msub x2,x8,x0,x7
ldr x10,[x30,#176]
mneg x7,x10,x0
stp x6,x7,[x30,#-224]
ldp x10,x11,[x30,#80]
ldp x20,x21,[x30,#128]
madd x10,x0,x20,x10
madd x11,x0,x21,x11
stp x10,x11,[x30,#-176]
umulh x10,x6,x9
ldr x11,[x16,x10,lsl #3]
ldp x8,x9,[x30,#64]
add x0,x2,x9
sub x0,x0,x11
asr x0,x0,#32
add x1,x2,x8
add x1,x1,x11
asr x1,x1,#32
cmp x0,x1
b.gt .Ltop_11
stp x1,x2,[x30,#-256]
add x30,sp,#2368
ldr x11,[sp,#8]
ldrb w9,[x11,w0,sxtw]
ldrh w10,[sp,#3862]
add x11,sp,#20,lsl #12
add x11,x11,x9,lsl #11
ldrh w10,[x11,x10,lsl #1]
strh w10,[sp,#3860]
b .Lent_10
//@end
//@seg top_10
seg_top_10:
.Ltop_10:
add x30,sp,#2368
ldr x0,[sp,#80]
ldp x1,x2,[x30,#0]
add x0,x0,#1
cmp x0,x1
b.gt .Ltop_11
ldrh w9,[sp,#3860]
add x10,sp,#20,lsl #12
ldrh w9,[x10,x9,lsl #1]
strh w9,[sp,#3860]
//@end
//@seg body_10
seg_body_10:
.Lent_10:
ldr x3,[sp]
tbnz x3,#10,.Lexact_10
str x0,[sp,#80]
subs x3,x2,x0,lsl #32
cneg x3,x3,mi
ldp x4,x5,[x30,#16]
subs x3,x3,x4
csel x3,x3,xzr,hi
lsr x3,x3,#13
mul x3,x3,x3
umulh x3,x3,x5
ldp x6,x7,[x30,#32]
subs x6,x6,x3
b.lo .Ltop_10
ldp x8,x9,[x30,#48]
msub x2,x8,x0,x7
ldr x10,[sp,#88]
ldp x11,x12,[x30,#176]
mneg x7,x11,x0
msub x7,x12,x10,x7
stp x6,x7,[x30,#-224]
ldp x10,x11,[x30,#80]
ldp x20,x21,[x30,#128]
madd x10,x0,x20,x10
madd x11,x0,x21,x11
stp x10,x11,[x30,#-176]
umulh x10,x6,x9
ldr x11,[x16,x10,lsl #3]
ldp x8,x9,[x30,#64]
add x0,x2,x9
sub x0,x0,x11
asr x0,x0,#32
add x1,x2,x8
add x1,x1,x11
asr x1,x1,#32
cmp x0,x1
b.gt .Ltop_10
stp x1,x2,[x30,#-256]
add x30,sp,#2112
ldr x11,[sp,#8]
ldrb w9,[x11,w0,sxtw]
ldrh w10,[sp,#3860]
add x11,sp,#18,lsl #12
add x11,x11,x9,lsl #11
ldrh w10,[x11,x10,lsl #1]
strh w10,[sp,#3858]
b .Lent_9
//@end
//@seg top_9
seg_top_9:
.Ltop_9:
add x30,sp,#2112
ldr x0,[sp,#72]
ldp x1,x2,[x30,#0]
add x0,x0,#1
cmp x0,x1
b.gt .Ltop_10
ldrh w9,[sp,#3858]
add x10,sp,#18,lsl #12
ldrh w9,[x10,x9,lsl #1]
strh w9,[sp,#3858]
//@end
//@seg body_9
seg_body_9:
.Lent_9:
ldr x3,[sp]
tbnz x3,#9,.Lexact_9
str x0,[sp,#72]
subs x3,x2,x0,lsl #32
cneg x3,x3,mi
ldp x4,x5,[x30,#16]
subs x3,x3,x4
csel x3,x3,xzr,hi
lsr x3,x3,#13
mul x3,x3,x3
umulh x3,x3,x5
ldp x6,x7,[x30,#32]
subs x6,x6,x3
b.lo .Ltop_9
ldp x8,x9,[x30,#48]
msub x2,x8,x0,x7
ldp x10,x11,[sp,#80]
ldp x12,x13,[x30,#176]
mneg x7,x12,x0
msub x7,x13,x10,x7
ldr x12,[x30,#192]
msub x7,x12,x11,x7
stp x6,x7,[x30,#-224]
ldp x10,x11,[x30,#80]
ldp x20,x21,[x30,#128]
madd x10,x0,x20,x10
madd x11,x0,x21,x11
stp x10,x11,[x30,#-176]
umulh x10,x6,x9
ldr x11,[x16,x10,lsl #3]
ldp x8,x9,[x30,#64]
add x0,x2,x9
sub x0,x0,x11
asr x0,x0,#32
add x1,x2,x8
add x1,x1,x11
asr x1,x1,#32
cmp x0,x1
b.gt .Ltop_9
stp x1,x2,[x30,#-256]
add x30,sp,#1856
ldr x11,[sp,#8]
ldrb w9,[x11,w0,sxtw]
ldrh w10,[sp,#3858]
add x11,sp,#16,lsl #12
add x11,x11,x9,lsl #11
ldrh w10,[x11,x10,lsl #1]
strh w10,[sp,#3856]
b .Lent_8
//@end
//@seg top_8
seg_top_8:
.Ltop_8:
add x30,sp,#1856
ldr x0,[sp,#64]
ldp x1,x2,[x30,#0]
add x0,x0,#1
cmp x0,x1
b.gt .Ltop_9
ldrh w9,[sp,#3856]
add x10,sp,#16,lsl #12
ldrh w9,[x10,x9,lsl #1]
strh w9,[sp,#3856]
//@end
//@seg body_8
seg_body_8:
.Lent_8:
ldr x3,[sp]
tbnz x3,#8,.Lexact_8
str x0,[sp,#64]
subs x3,x2,x0,lsl #32
cneg x3,x3,mi
ldp x4,x5,[x30,#16]
subs x3,x3,x4
csel x3,x3,xzr,hi
lsr x3,x3,#13
mul x3,x3,x3
umulh x3,x3,x5
ldp x6,x7,[x30,#32]
subs x6,x6,x3
b.lo .Ltop_8
ldp x8,x9,[x30,#48]
msub x2,x8,x0,x7
ldp x10,x11,[sp,#72]
ldr x12,[sp,#88]
ldp x13,x14,[x30,#176]
mneg x7,x13,x0
msub x7,x14,x10,x7
ldp x13,x14,[x30,#192]
msub x7,x13,x11,x7
msub x7,x14,x12,x7
stp x6,x7,[x30,#-224]
ldp x10,x11,[x30,#80]
ldp x20,x21,[x30,#128]
madd x10,x0,x20,x10
madd x11,x0,x21,x11
stp x10,x11,[x30,#-176]
umulh x10,x6,x9
ldr x11,[x16,x10,lsl #3]
ldp x8,x9,[x30,#64]
add x0,x2,x9
sub x0,x0,x11
asr x0,x0,#32
add x1,x2,x8
add x1,x1,x11
asr x1,x1,#32
cmp x0,x1
b.gt .Ltop_8
stp x1,x2,[x30,#-256]
add x30,sp,#1600
ldr x11,[sp,#8]
ldrb w9,[x11,w0,sxtw]
ldrh w10,[sp,#3856]
add x11,sp,#14,lsl #12
add x11,x11,x9,lsl #11
ldrh w10,[x11,x10,lsl #1]
strh w10,[sp,#3854]
b .Lent_7
//@end
//@seg top_7
seg_top_7:
.Ltop_7:
add x30,sp,#1600
ldr x0,[sp,#56]
ldp x1,x2,[x30,#0]
add x0,x0,#1
cmp x0,x1
b.gt .Ltop_8
ldrh w9,[sp,#3854]
add x10,sp,#14,lsl #12
ldrh w9,[x10,x9,lsl #1]
strh w9,[sp,#3854]
//@end
//@seg body_7
seg_body_7:
.Lent_7:
ldr x3,[sp]
tbnz x3,#7,.Lexact_7
str x0,[sp,#56]
subs x3,x2,x0,lsl #32
cneg x3,x3,mi
ldp x4,x5,[x30,#16]
subs x3,x3,x4
csel x3,x3,xzr,hi
lsr x3,x3,#13
mul x3,x3,x3
umulh x3,x3,x5
ldp x6,x7,[x30,#32]
subs x6,x6,x3
b.lo .Ltop_7
ldp x8,x9,[x30,#48]
msub x2,x8,x0,x7
ldp x10,x11,[sp,#64]
ldp x12,x13,[sp,#80]
ldp x14,x19,[x30,#176]
mneg x7,x14,x0
msub x7,x19,x10,x7
ldp x14,x19,[x30,#192]
msub x7,x14,x11,x7
msub x7,x19,x12,x7
ldr x14,[x30,#208]
msub x7,x14,x13,x7
stp x6,x7,[x30,#-224]
ldp x10,x11,[x30,#80]
ldp x20,x21,[x30,#128]
madd x10,x0,x20,x10
madd x11,x0,x21,x11
stp x10,x11,[x30,#-176]
umulh x10,x6,x9
ldr x11,[x16,x10,lsl #3]
ldp x8,x9,[x30,#64]
add x0,x2,x9
sub x0,x0,x11
asr x0,x0,#32
add x1,x2,x8
add x1,x1,x11
asr x1,x1,#32
cmp x0,x1
b.gt .Ltop_7
stp x1,x2,[x30,#-256]
add x30,sp,#1344
ldr x11,[sp,#8]
ldrb w9,[x11,w0,sxtw]
ldrh w10,[sp,#3854]
add x11,sp,#12,lsl #12
add x11,x11,x9,lsl #11
ldrh w10,[x11,x10,lsl #1]
strh w10,[sp,#3852]
b .Lent_6
//@end
//@seg top_6
seg_top_6:
.Ltop_6:
add x30,sp,#1344
ldr x0,[sp,#48]
ldp x1,x2,[x30,#0]
add x0,x0,#1
cmp x0,x1
b.gt .Ltop_7
ldrh w9,[sp,#3852]
add x10,sp,#12,lsl #12
ldrh w9,[x10,x9,lsl #1]
strh w9,[sp,#3852]
//@end
//@seg body_6
seg_body_6:
.Lent_6:
ldr x3,[sp]
tbnz x3,#6,.Lexact_6
str x0,[sp,#48]
subs x3,x2,x0,lsl #32
cneg x3,x3,mi
ldp x4,x5,[x30,#16]
subs x3,x3,x4
csel x3,x3,xzr,hi
lsr x3,x3,#13
mul x3,x3,x3
umulh x3,x3,x5
ldp x6,x7,[x30,#32]
subs x6,x6,x3
b.lo .Ltop_6
ldp x8,x9,[x30,#48]
msub x2,x8,x0,x7
ldp x10,x11,[sp,#56]
ldp x12,x13,[sp,#72]
ldr x14,[sp,#88]
ldp x19,x20,[x30,#176]
mneg x7,x19,x0
msub x7,x20,x10,x7
ldp x19,x20,[x30,#192]
msub x7,x19,x11,x7
msub x7,x20,x12,x7
ldp x19,x20,[x30,#208]
msub x7,x19,x13,x7
msub x7,x20,x14,x7
stp x6,x7,[x30,#-224]
ldp x10,x11,[x30,#80]
ldp x20,x21,[x30,#128]
madd x10,x0,x20,x10
madd x11,x0,x21,x11
stp x10,x11,[x30,#-176]
umulh x10,x6,x9
ldr x11,[x16,x10,lsl #3]
ldp x8,x9,[x30,#64]
add x0,x2,x9
sub x0,x0,x11
asr x0,x0,#32
add x1,x2,x8
add x1,x1,x11
asr x1,x1,#32
cmp x0,x1
b.gt .Ltop_6
stp x1,x2,[x30,#-256]
add x30,sp,#1088
ldr x11,[sp,#8]
ldrb w9,[x11,w0,sxtw]
ldrh w10,[sp,#3852]
add x11,sp,#10,lsl #12
add x11,x11,x9,lsl #11
ldrh w10,[x11,x10,lsl #1]
strh w10,[sp,#3850]
b .Lent_5
//@end
//@seg top_5
seg_top_5:
.Ltop_5:
add x30,sp,#1088
ldr x0,[sp,#40]
ldp x1,x2,[x30,#0]
add x0,x0,#1
cmp x0,x1
b.gt .Ltop_6
ldrh w9,[sp,#3850]
add x10,sp,#10,lsl #12
ldrh w9,[x10,x9,lsl #1]
strh w9,[sp,#3850]
//@end
//@seg body_5
seg_body_5:
.Lent_5:
ldr x3,[sp]
tbnz x3,#5,.Lexact_5
str x0,[sp,#40]
subs x3,x2,x0,lsl #32
cneg x3,x3,mi
ldp x4,x5,[x30,#16]
subs x3,x3,x4
csel x3,x3,xzr,hi
lsr x3,x3,#13
mul x3,x3,x3
umulh x3,x3,x5
ldp x6,x7,[x30,#32]
subs x6,x6,x3
b.lo .Ltop_5
ldp x8,x9,[x30,#48]
msub x2,x8,x0,x7
ldp x10,x11,[sp,#48]
ldp x12,x13,[sp,#64]
ldp x14,x19,[sp,#80]
ldp x20,x21,[x30,#176]
mneg x7,x20,x0
msub x7,x21,x10,x7
ldp x20,x21,[x30,#192]
msub x7,x20,x11,x7
msub x7,x21,x12,x7
ldp x20,x21,[x30,#208]
msub x7,x20,x13,x7
msub x7,x21,x14,x7
ldr x20,[x30,#224]
msub x7,x20,x19,x7
stp x6,x7,[x30,#-224]
ldp x10,x11,[x30,#80]
ldp x20,x21,[x30,#128]
madd x10,x0,x20,x10
madd x11,x0,x21,x11
stp x10,x11,[x30,#-176]
umulh x10,x6,x9
ldr x11,[x16,x10,lsl #3]
ldp x8,x9,[x30,#64]
add x0,x2,x9
sub x0,x0,x11
asr x0,x0,#32
add x1,x2,x8
add x1,x1,x11
asr x1,x1,#32
cmp x0,x1
b.gt .Ltop_5
stp x1,x2,[x30,#-256]
add x30,sp,#832
ldr x11,[sp,#8]
ldrb w9,[x11,w0,sxtw]
ldrh w10,[sp,#3850]
add x11,sp,#8,lsl #12
add x11,x11,x9,lsl #11
ldrh w10,[x11,x10,lsl #1]
strh w10,[sp,#3848]
b .Lent_4
//@end
//@seg top_4
seg_top_4:
.Ltop_4:
add x30,sp,#832
ldr x0,[sp,#32]
ldp x1,x2,[x30,#0]
add x0,x0,#1
cmp x0,x1
b.gt .Ltop_5
ldrh w9,[sp,#3848]
add x10,sp,#8,lsl #12
ldrh w9,[x10,x9,lsl #1]
strh w9,[sp,#3848]
//@end
//@seg body_4
seg_body_4:
.Lent_4:
ldr x3,[sp]
tbnz x3,#4,.Lexact_4
str x0,[sp,#32]
subs x3,x2,x0,lsl #32
cneg x3,x3,mi
ldp x4,x5,[x30,#16]
subs x3,x3,x4
csel x3,x3,xzr,hi
lsr x3,x3,#13
mul x3,x3,x3
umulh x3,x3,x5
ldp x6,x7,[x30,#32]
subs x6,x6,x3
b.lo .Ltop_4
ldp x8,x9,[x30,#48]
msub x2,x8,x0,x7
ldp x10,x11,[sp,#40]
ldp x12,x13,[sp,#56]
ldp x14,x19,[sp,#72]
ldr x20,[sp,#88]
ldp x21,x22,[x30,#176]
mneg x7,x21,x0
msub x7,x22,x10,x7
ldp x21,x22,[x30,#192]
msub x7,x21,x11,x7
msub x7,x22,x12,x7
ldp x21,x22,[x30,#208]
msub x7,x21,x13,x7
msub x7,x22,x14,x7
ldp x21,x22,[x30,#224]
msub x7,x21,x19,x7
msub x7,x22,x20,x7
stp x6,x7,[x30,#-288]
ldp x10,x11,[x30,#80]
ldp x20,x21,[x30,#128]
madd x10,x0,x20,x10
madd x11,x0,x21,x11
stp x10,x11,[x30,#-240]
umulh x10,x6,x9
ldr x11,[x16,x10,lsl #3]
ldp x8,x9,[x30,#64]
add x0,x2,x9
sub x0,x0,x11
asr x0,x0,#32
add x1,x2,x8
add x1,x1,x11
asr x1,x1,#32
cmp x0,x1
b.gt .Ltop_4
stp x1,x2,[x30,#-320]
add x30,sp,#512
ldr x11,[sp,#8]
ldrb w9,[x11,w0,sxtw]
ldrh w10,[sp,#3848]
add x11,sp,#6,lsl #12
add x11,x11,x9,lsl #11
ldrh w10,[x11,x10,lsl #1]
strh w10,[sp,#3846]
b .Lent_3
//@end
//@seg top_3
seg_top_3:
.Ltop_3:
add x30,sp,#512
ldr x0,[sp,#24]
ldp x1,x2,[x30,#0]
add x0,x0,#1
cmp x0,x1
b.gt .Ltop_4
ldrh w9,[sp,#3846]
add x10,sp,#6,lsl #12
ldrh w9,[x10,x9,lsl #1]
strh w9,[sp,#3846]
//@end
//@seg body_3
seg_body_3:
.Lent_3:
ldr x3,[sp]
tbnz x3,#3,.Lexact_3
str x0,[sp,#24]
subs x3,x2,x0,lsl #32
cneg x3,x3,mi
ldp x4,x5,[x30,#16]
subs x3,x3,x4
csel x3,x3,xzr,hi
lsr x3,x3,#13
mul x3,x3,x3
umulh x3,x3,x5
ldp x6,x7,[x30,#32]
subs x6,x6,x3
b.lo .Ltop_3
ldp x8,x9,[x30,#48]
msub x2,x8,x0,x7
ldp x10,x11,[sp,#32]
ldp x12,x13,[sp,#48]
ldp x14,x19,[sp,#64]
ldp x20,x21,[sp,#80]
ldp x22,x23,[x30,#176]
mneg x7,x22,x0
msub x7,x23,x10,x7
ldp x22,x23,[x30,#192]
msub x7,x22,x11,x7
msub x7,x23,x12,x7
ldp x22,x23,[x30,#208]
msub x7,x22,x13,x7
msub x7,x23,x14,x7
ldp x22,x23,[x30,#224]
msub x7,x22,x19,x7
msub x7,x23,x20,x7
ldr x22,[x30,#240]
msub x7,x22,x21,x7
ldp x22,x23,[x30,#248]
mneg x4,x22,x0
msub x4,x23,x10,x4
ldp x22,x23,[x30,#264]
msub x4,x22,x11,x4
msub x4,x23,x12,x4
ldp x22,x23,[x30,#280]
msub x4,x22,x13,x4
msub x4,x23,x14,x4
ldp x22,x23,[x30,#296]
msub x4,x22,x19,x4
msub x4,x23,x20,x4
ldr x22,[x30,#312]
msub x4,x22,x21,x4
stp x6,x7,[sp,#336]
str x4,[sp,#376]
ldp x10,x11,[x30,#80]
ldp x20,x21,[x30,#128]
madd x10,x0,x20,x10
madd x11,x0,x21,x11
stp x10,x11,[sp,#384]
umulh x10,x6,x9
ldr x11,[x16,x10,lsl #3]
ldp x8,x9,[x30,#64]
add x0,x2,x9
sub x0,x0,x11
asr x0,x0,#32
add x1,x2,x8
add x1,x1,x11
asr x1,x1,#32
cmp x0,x1
b.gt .Ltop_3
stp x1,x2,[sp,#304]
ldr x11,[sp,#8]
ldrb w9,[x11,w0,sxtw]
ldrh w10,[sp,#3846]
add x11,sp,#4,lsl #12
add x11,x11,x9,lsl #11
ldrh w10,[x11,x10,lsl #1]
strh w10,[sp,#3844]
ldp x11,x12,[sp,#128]
mov x1,#8
ldr x10,[sp,#408]
b .Lent_2
//@end
//@seg top_2
seg_top_2:
.Ltop_2:
ldr x0,[sp,#16]
ldp x14,x2,[sp,#304]
add x0,x0,#1
cmp x0,x14
b.gt .Ltop_3
ldrh w9,[sp,#3844]
add x7,sp,#4,lsl #12
ldrh w9,[x7,x9,lsl #1]
strh w9,[sp,#3844]
//@end
//@seg body_2
seg_body_2:
.Lent_2:
str x0,[sp,#16]
subs x3,x2,x0,lsl #32
cneg x3,x3,mi
ldp x4,x5,[sp,#320]
subs x3,x3,x4
csel x3,x3,xzr,hi
lsr x3,x3,#13
mul x3,x3,x3
umulh x3,x3,x5
ldp x6,x7,[sp,#336]
subs x6,x6,x3
b.lo .Ltop_2
ldp x8,x9,[sp,#352]
msub x2,x8,x0,x7
ldp x8,x21,[sp,#368]
msub x21,x8,x0,x21
umulh x7,x6,x9
ldr x5,[x16,x7,lsl #3]
ldp x8,x9,[sp,#272]
add x14,x2,x9
sub x14,x14,x5
asr x14,x14,#32
add x3,x2,x8
add x3,x3,x5
asr x3,x3,#32
cmp x14,x3
b.gt .Ltop_2
ldp x8,x9,[sp,#208]
sub x23,x3,x14
madd x7,x23,x8,x9
add x19,x6,x7,lsr #32
subs x3,x2,x14,lsl #32
cneg x4,x3,mi
lsr x4,x4,#13
mul x4,x4,x4
ldp x8,x9,[sp,#224]
umulh x4,x4,x8
sub x19,x19,x4
ldp x8,x30,[sp,#240]
sub x5,x8,x3,lsl #1
mul x7,x5,x9
smulh x8,x5,x9
extr x20,x8,x7,#32
msub x7,x29,x14,x21
ldp x8,x9,[sp,#256]
add x21,x7,x8
add x22,x7,x9
ldp x24,x25,[sp,#384]
ldp x3,x4,[sp,#432]
madd x24,x0,x3,x24
madd x25,x0,x4,x25
ldp x4,x13,[sp,#160]
madd x24,x14,x4,x24
madd x25,x14,x13,x25
ldr x8,[sp,#8]
ldrb w9,[x8,w14,sxtw]
ldrh w3,[sp,#3844]
add x8,sp,#2,lsl #12
add x8,x8,x9,lsl #11
ldrh w14,[x8,x3,lsl #1]
ldr x8,[sp,#416]
str x8,[x17]
//@end
//@seg k0_body
seg_k0_body:
.Lk0_top:
ldrh w14,[x15,x14,lsl #1]
tbnz x19,#63,.Lk0_skip
umulh x8,x19,x30
ldr x8,[x16,x8,lsl #3]
sub x6,x21,x8
asr x6,x6,#32
add x7,x22,x8
asr x7,x7,#32
subs x0,x7,x6
b.lt .Lk0_skip
madd x2,x6,x26,x24
madd x5,x6,x27,x25
sub x19,x19,x20
add x20,x20,x28
sub x21,x21,x29
sub x22,x22,x29
add x24,x24,x4
add x25,x25,x13
subs x23,x23,#1
//@end
.ifdef FP_TEST
stp x2,x5,[x15],#16
stp x6,x7,[x15],#16
.endif
.Lchord_split:
//@seg wheel_disp4
wheel_disp4:
bfi w6,w14,#23,#3
ldr x8,[x10,w6,sxtw #3]
ldr x9,[x8,w0,uxtw #3]
br x9
//@end
//@seg k0_back
seg_k0_back:
.Lk0_back:
b.pl .Lk0_top
//@end
//@seg k1_exit
seg_k1_exit:
b .Ltop_2
//@end
//@seg chord_chunk
seg_chord_chunk:
.Lchord_chunk:
str x30,[sp,#488]
sbfx x6,x6,#0,#23
.Lchunk_loop:
add x7,x2,x26,lsl #7
add x8,x5,x27,lsl #7
stp x7,x8,[sp,#448]
sub x7,x0,#128
add x8,x6,#128
stp x7,x8,[sp,#464]
bfi w6,w14,#23,#3
ldr x8,[x10,w6,sxtw #3]
ldr x9,[x8,#6248]
blr x9
ldp x2,x5,[sp,#448]
ldp x0,x6,[sp,#464]
cmp x0,#127
b.hi .Lchunk_loop
ldr x30,[sp,#488]
cmp x23,#0
b wheel_disp4
//@end
//@seg k0_skip
seg_k0_skip:
.Lk0_skip:
sub x19,x19,x20
add x20,x20,x28
sub x21,x21,x29
sub x22,x22,x29
add x24,x24,x4
add x25,x25,x13
subs x23,x23,#1
b.pl .Lk0_top
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
ldp  x8,x9,[x0],#16
sxth x10,w8
mul  x11,x10,x20
sbfx x10,x8,#16,#16
madd x11,x10,x21,x11
sbfx x10,x8,#32,#16
madd x11,x10,x22,x11
asr  x10,x8,#48
madd x11,x10,x23,x11
sxth x10,w9
madd x11,x10,x24,x11
sbfx x10,x9,#16,#16
madd x11,x10,x25,x11
add  x11,x11,x26
umulh x12,x11,x27
msub x11,x12,x5,x11
subs x12,x11,x5
csel x11,x12,x11,hs
.Lent:
ldr  x13,[x2],#8
str  x11,[x3,x13,lsl #3]
add  x11,x11,x4
subs x12,x11,x5
csel x11,x12,x11,hs
tbz  x13,#61,.Lent
cmp  x0,x1
b.lo .Lpiece
ret
//@seg chord_enter
chord_enter:
//@end
//@seg wheel_skip1_stub
wheel_skip1_stub:
add x2,x2,x26
add x5,x5,x27
//@end
//@seg emit9_step1
emit9_step1:
ldp x0,x6,[x2]
ldr x3,[x6,x5,lsl #3]
add x8,x17,x3,lsr #27
add x7,x0,x3,lsl #38
ldadd x1,x9,[x8]
str x7,[x9]
add x2,x2,x26
add x5,x5,x27
tbnz x9,#14,.Lfl_s1
//@end
.Lcont_s1:
//@seg emit9_step2
emit9_step2:
ldp x0,x6,[x2]
ldr x3,[x6,x5,lsl #3]
add x8,x17,x3,lsr #27
add x7,x0,x3,lsl #38
ldadd x1,x9,[x8]
str x7,[x9]
add x2,x2,x26,lsl #1
add x5,x5,x27,lsl #1
tbnz x9,#14,.Lfl_s2
//@end
.Lcont_s2:
//@seg emit9_step3
emit9_step3:
ldp x0,x6,[x2]
ldr x3,[x6,x5,lsl #3]
add x8,x17,x3,lsr #27
add x7,x0,x3,lsl #38
ldadd x1,x9,[x8]
str x7,[x9]
add x2,x2,x11
add x5,x5,x12
tbnz x9,#14,.Lfl_s3
//@end
.Lcont_s3:
//@seg emit9_last
emit9_last:
ldp x0,x6,[x2]
ldr x3,[x6,x5,lsl #3]
add x8,x17,x3,lsr #27
add x7,x0,x3,lsl #38
ldadd x1,x9,[x8]
str x7,[x9]
tbnz x9,#14,.Lfl_last
//@end
.Lcont_last:
b.pl .Lk0_top
b .Ltop_2
//@seg wheel_skip2_stub
wheel_skip2_stub:
add x2,x2,x26,lsl #1
add x5,x5,x27,lsl #1
//@end
//@seg emit8_step1
emit8_step1:
ldp x0,x6,[x2]
ldr x3,[x6,x5,lsl #3]
add x8,x17,x3,lsr #27
add x7,x0,x3,lsl #38
ldadd x1,x9,[x8]
str x7,[x9]
add x2,x2,x26
add x5,x5,x27
//@end
//@seg emit8_step2
emit8_step2:
ldp x0,x6,[x2]
ldr x3,[x6,x5,lsl #3]
add x8,x17,x3,lsr #27
add x7,x0,x3,lsl #38
ldadd x1,x9,[x8]
str x7,[x9]
add x2,x2,x26,lsl #1
add x5,x5,x27,lsl #1
//@end
//@seg emit8_step3
emit8_step3:
ldp x0,x6,[x2]
ldr x3,[x6,x5,lsl #3]
add x8,x17,x3,lsr #27
add x7,x0,x3,lsl #38
ldadd x1,x9,[x8]
str x7,[x9]
add x2,x2,x11
add x5,x5,x12
//@end
//@seg emit8_last
emit8_last:
ldp x0,x6,[x2]
ldr x3,[x6,x5,lsl #3]
add x8,x17,x3,lsr #27
add x7,x0,x3,lsl #38
ldadd x1,x9,[x8]
str x7,[x9]
//@end
b.pl .Lk0_top
b .Ltop_2
//@seg emit8_epoch_flush
emit8_epoch_flush:
.Lfl_loop:
ldr x3,[x0]
cmp x3,x1
b.eq .Lfl_empty
str x3,[x2]
str x1,[x0]
.Lfl_empty:
add x0,x0,#8
add x1,x1,x4
add x2,x2,#8
cmp x0,x5
b.lo .Lfl_loop
//@end
ret
//@seg emit_fl_s1
.Lfl_s1:
sub x0,x9,#4,lsl #12
add x0,x0,#8
mov x3,#2046
ldr x7,[sp,#480]
stp x3,x7,[x0]
add x9,x0,#16
str x9,[x8]
b .Lcont_s1
//@end
.Lfl_s2:
sub x0,x9,#4,lsl #12
add x0,x0,#8
mov x3,#2046
ldr x7,[sp,#480]
stp x3,x7,[x0]
add x9,x0,#16
str x9,[x8]
b .Lcont_s2
.Lfl_s3:
sub x0,x9,#4,lsl #12
add x0,x0,#8
mov x3,#2046
ldr x7,[sp,#480]
stp x3,x7,[x0]
add x9,x0,#16
str x9,[x8]
b .Lcont_s3
.Lfl_last:
sub x0,x9,#4,lsl #12
add x0,x0,#8
mov x3,#2046
ldr x7,[sp,#480]
stp x3,x7,[x0]
add x9,x0,#16
str x9,[x8]
b .Lcont_last
.Lmask0:
.Lchunk128_m0:
.Lchunk128_s0:
.Lchain_entry:
b.pl .Lk0_top
b .Ltop_2
.altmacro
.macro E8 stride
ldp x0,x6,[x2]
ldr x3,[x6,x5,lsl #3]
add x8,x17,x3,lsr #27
add x7,x0,x3,lsl #38
ldadd x1,x9,[x8]
str x7,[x9]
.if \stride == 1
add x2,x2,x26
add x5,x5,x27
.elseif \stride == 2
add x2,x2,x26,lsl #1
add x5,x5,x27,lsl #1
.elseif \stride == 3
add x2,x2,x11
add x5,x5,x12
.endif
.endm
.macro E9 stride,id
ldp x0,x6,[x2]
ldr x3,[x6,x5,lsl #3]
add x8,x17,x3,lsr #27
add x7,x0,x3,lsl #38
ldadd x1,x9,[x8]
str x7,[x9]
.if \stride == 1
add x2,x2,x26
add x5,x5,x27
.elseif \stride == 2
add x2,x2,x26,lsl #1
add x5,x5,x27,lsl #1
.elseif \stride == 3
add x2,x2,x11
add x5,x5,x12
.endif
tbnz x9,#14,.Lfl_g_\id
.Lct_g_\id:
.endm
.macro E9C id
.Lfl_g_\id:
sub x0,x9,#4,lsl #12
add x0,x0,#8
mov x3,#2046
ldr x7,[sp,#480]
stp x3,x7,[x0]
add x9,x0,#16
str x9,[x8]
b .Lct_g_\id
.endm
.macro STUB8 sk,pat,nk
e8_s\sk\()_p\pat\()_k\nk:
.if \sk == 1
add x2,x2,x26
add x5,x5,x27
.elseif \sk == 2
add x2,x2,x26,lsl #1
add x5,x5,x27,lsl #1
.endif
.if \pat == 0
.rept \nk - 1
E8 1
.endr
.elseif \pat == 3
.rept \nk - 1
E8 3
.endr
.elseif \pat == 1
.set _k,0
.rept \nk - 1
.if (((\nk - 1) - _k) & 1) == 1
E8 1
.else
E8 2
.endif
.set _k,_k + 1
.endr
.elseif \pat == 2
.set _k,0
.rept \nk - 1
.if (((\nk - 1) - _k) & 1) == 1
E8 2
.else
E8 1
.endif
.set _k,_k + 1
.endr
.endif
E8 0
b.pl .Lk0_top
b .Ltop_2
.endm
.set _gid,0
.macro STUB9 sk,pat,nk
e9_s\sk\()_p\pat\()_k\nk:
.set _s0,_gid
.if \sk == 1
add x2,x2,x26
add x5,x5,x27
.elseif \sk == 2
add x2,x2,x26,lsl #1
add x5,x5,x27,lsl #1
.endif
.if \pat == 0
.rept \nk - 1
E9 1,%_gid
.set _gid,_gid + 1
.endr
.elseif \pat == 3
.rept \nk - 1
E9 3,%_gid
.set _gid,_gid + 1
.endr
.elseif \pat == 1
.set _k,0
.rept \nk - 1
.if (((\nk - 1) - _k) & 1) == 1
E9 1,%_gid
.else
E9 2,%_gid
.endif
.set _gid,_gid + 1
.set _k,_k + 1
.endr
.elseif \pat == 2
.set _k,0
.rept \nk - 1
.if (((\nk - 1) - _k) & 1) == 1
E9 2,%_gid
.else
E9 1,%_gid
.endif
.set _gid,_gid + 1
.set _k,_k + 1
.endr
.endif
E9 0,%_gid
.set _gid,_gid + 1
b.pl .Lk0_top
b .Ltop_2
.set _c,_s0
.rept \nk
E9C %_c
.set _c,_c + 1
.endr
.endm
wheel_skip1_emit8_last:
STUB8 1,0,1
wheel_skip2_emit8_last:
STUB8 2,0,1
wheel_skip1_emit9_last:
STUB9 1,0,1
wheel_skip2_emit9_last:
STUB9 2,0,1
e8_s0_p0: STUB8 0,0,128
e8_s0_p1: STUB8 0,1,86
e8_s0_p2: STUB8 0,2,86
e8_s0_p3: STUB8 0,3,43
e9_s0_p0: STUB9 0,0,128
e9_s0_p1: STUB9 0,1,86
e9_s0_p2: STUB9 0,2,86
e9_s0_p3: STUB9 0,3,43
.set _n,1
.rept 85
STUB8 1,1,%_n
STUB8 1,2,%_n
.set _n,_n + 1
.endr
.set _n,1
.rept 43
STUB8 1,3,%_n
STUB8 2,3,%_n
.set _n,_n + 1
.endr
.set _n,1
.rept 85
STUB9 1,1,%_n
STUB9 1,2,%_n
.set _n,_n + 1
.endr
.set _n,1
.rept 43
STUB9 1,3,%_n
STUB9 2,3,%_n
.set _n,_n + 1
.endr
.macro CHUNK8 sk,pat,nk
c8_s\sk\()_p\pat\()_k\nk:
.if \sk == 1
add x2,x2,x26
add x5,x5,x27
.elseif \sk == 2
add x2,x2,x26,lsl #1
add x5,x5,x27,lsl #1
.endif
.if \pat == 0
.rept \nk - 1
E8 1
.endr
.elseif \pat == 3
.rept \nk - 1
E8 3
.endr
.elseif \pat == 1
.set _k,0
.rept \nk - 1
.if (((\nk - 1) - _k) & 1) == 1
E8 1
.else
E8 2
.endif
.set _k,_k + 1
.endr
.elseif \pat == 2
.set _k,0
.rept \nk - 1
.if (((\nk - 1) - _k) & 1) == 1
E8 2
.else
E8 1
.endif
.set _k,_k + 1
.endr
.endif
E8 0
ret
.endm
.macro CHUNK9 sk,pat,nk
c9_s\sk\()_p\pat\()_k\nk:
.set _s0,_gid
.if \sk == 1
add x2,x2,x26
add x5,x5,x27
.elseif \sk == 2
add x2,x2,x26,lsl #1
add x5,x5,x27,lsl #1
.endif
.if \pat == 0
.rept \nk - 1
E9 1,%_gid
.set _gid,_gid + 1
.endr
.elseif \pat == 3
.rept \nk - 1
E9 3,%_gid
.set _gid,_gid + 1
.endr
.elseif \pat == 1
.set _k,0
.rept \nk - 1
.if (((\nk - 1) - _k) & 1) == 1
E9 1,%_gid
.else
E9 2,%_gid
.endif
.set _gid,_gid + 1
.set _k,_k + 1
.endr
.elseif \pat == 2
.set _k,0
.rept \nk - 1
.if (((\nk - 1) - _k) & 1) == 1
E9 2,%_gid
.else
E9 1,%_gid
.endif
.set _gid,_gid + 1
.set _k,_k + 1
.endr
.endif
E9 0,%_gid
.set _gid,_gid + 1
ret
.set _c,_s0
.rept \nk
E9C %_c
.set _c,_c + 1
.endr
.endm
CHUNK8 0,0,128
CHUNK8 0,1,86
CHUNK8 0,1,85
CHUNK8 1,2,85
CHUNK8 0,3,43
CHUNK8 1,3,43
CHUNK8 2,3,42
CHUNK9 0,0,128
CHUNK9 0,1,86
CHUNK9 0,1,85
CHUNK9 1,2,85
CHUNK9 0,3,43
CHUNK9 1,3,43
CHUNK9 2,3,42
.noaltmacro
//@seg emit_direct_flush
emit_direct_flush:
sub x0,x9,#4,lsl #12
add x0,x0,#8
mov x3,#2046
ldr x7,[sp,#480]
stp x3,x7,[x0]
add x9,x0,#16
str x9,[x8]
ret
//@end
//@seg lat_apply_page
apply_page:
ldp x11,x6,[x28]
add x4,x28,#16
add x8,x4,x11,lsl #3
tbz x11,#0,apply_even
ldr x0,[x4],#8
staddb w6,[x0]
apply_even:
ubfx x12,x11,#1,#7
adr x13,apply_end
sub x13,x13,x12,lsl #3
sub x13,x13,x12,lsl #2
br x13
//@end
//@seg lat_apply
apply:
.rept 128
ldp x0,x1,[x4],#16
staddb w6,[x0]
staddb w6,[x1]
.endr
apply_end:
cmp x4,x8
b.lo apply
ret
//@end
//@seg x16_bfi_copy_page
x16_bfi_copy_page:
cbz x11,copy_ret
add x16,x15,x11,lsl #3
tbz x11,#0,copy_even
ldr x0,[x15],#8
str x0,[x25],#8
copy_even:
ubfx x12,x11,#1,#6
adr x13,copy_end
sub x13,x13,x12,lsl #3
br x13
//@end
//@seg x16_bfi_copy
x16_bfi_copy:
.rept 64
ldp x0,x1,[x15],#16
stp x0,x1,[x25],#16
.endr
copy_end:
cmp x15,x16
b.lo x16_bfi_copy
copy_ret:
ret
//@end
//@seg x16_nobfi_pos_pre
seg_x16_nobfi_pos_pre:
sub x5,x23,x22
lsr x5,x5,#3
bic x21,x5,#15
add x21,x22,x21,lsl #3
cmp x22,x21
b.hs .Lnobfi_pos_tail
//@end
//@seg x16_nobfi_pos
seg_x16_nobfi_pos:
.Lnobfi_pos:
.rept 8
ldp x0,x1,[x22],#16
add x2,x10,x0
staddb w19,[x2]
add x3,x10,x1
staddb w19,[x3]
.endr
cmp x22,x21
b.lo .Lnobfi_pos
//@end
//@seg x16_nobfi_pos_tail
.Lnobfi_pos_tail:
tbz x5,#3,3f
.rept 4
ldp x0,x1,[x22],#16
add x2,x10,x0
staddb w19,[x2]
add x3,x10,x1
staddb w19,[x3]
.endr
3:
tbz x5,#2,4f
.rept 2
ldp x0,x1,[x22],#16
add x2,x10,x0
staddb w19,[x2]
add x3,x10,x1
staddb w19,[x3]
.endr
4:
tbz x5,#1,5f
ldp x0,x1,[x22],#16
add x2,x10,x0
staddb w19,[x2]
add x3,x10,x1
staddb w19,[x3]
5:
tbz x5,#0,6f
ldr x0,[x22],#8
add x2,x10,x0
staddb w19,[x2]
6:
ret
//@end
//@seg x16_nobfi_bnd1_hi_pre
seg_x16_nobfi_bnd1_hi_pre:
sub x5,x23,x22
lsr x5,x5,#3
bic x21,x5,#15
add x21,x22,x21,lsl #3
cmp x22,x21
b.hs .Lnobfi_bnd1_hi_tail
//@end
//@seg x16_nobfi_bnd1_hi
seg_x16_nobfi_bnd1_hi:
.Lnobfi_bnd1_hi:
.rept 8
ldp x0,x1,[x22],#16
cmp x25,x0
b.lo 1f
add x2,x10,x0
staddb w19,[x2]
1:
cmp x25,x1
b.lo 2f
add x3,x10,x1
staddb w19,[x3]
2:
.endr
cmp x22,x21
b.lo .Lnobfi_bnd1_hi
//@end
//@seg x16_nobfi_bnd1_hi_tail
.Lnobfi_bnd1_hi_tail:
tbz x5,#3,3f
.rept 4
ldp x0,x1,[x22],#16
cmp x25,x0
b.lo 1f
add x2,x10,x0
staddb w19,[x2]
1:
cmp x25,x1
b.lo 2f
add x3,x10,x1
staddb w19,[x3]
2:
.endr
3:
tbz x5,#2,4f
.rept 2
ldp x0,x1,[x22],#16
cmp x25,x0
b.lo 1f
add x2,x10,x0
staddb w19,[x2]
1:
cmp x25,x1
b.lo 2f
add x3,x10,x1
staddb w19,[x3]
2:
.endr
4:
tbz x5,#1,5f
ldp x0,x1,[x22],#16
cmp x25,x0
b.lo 1f
add x2,x10,x0
staddb w19,[x2]
1:
cmp x25,x1
b.lo 2f
add x3,x10,x1
staddb w19,[x3]
2:
5:
tbz x5,#0,6f
ldr x0,[x22],#8
cmp x25,x0
b.lo 6f
add x2,x10,x0
staddb w19,[x2]
6:
ret
//@end
//@seg x16_nobfi_bnd1_lo_pre
seg_x16_nobfi_bnd1_lo_pre:
sub x5,x23,x22
lsr x5,x5,#3
bic x21,x5,#15
add x21,x22,x21,lsl #3
cmp x22,x21
b.hs .Lnobfi_bnd1_lo_tail
//@end
//@seg x16_nobfi_bnd1_lo
seg_x16_nobfi_bnd1_lo:
.Lnobfi_bnd1_lo:
.rept 8
ldp x0,x1,[x22],#16
cmp x25,x0
b.hi 1f
add x2,x10,x0
staddb w19,[x2]
1:
cmp x25,x1
b.hi 2f
add x3,x10,x1
staddb w19,[x3]
2:
.endr
cmp x22,x21
b.lo .Lnobfi_bnd1_lo
//@end
//@seg x16_nobfi_bnd1_lo_tail
.Lnobfi_bnd1_lo_tail:
tbz x5,#3,3f
.rept 4
ldp x0,x1,[x22],#16
cmp x25,x0
b.hi 1f
add x2,x10,x0
staddb w19,[x2]
1:
cmp x25,x1
b.hi 2f
add x3,x10,x1
staddb w19,[x3]
2:
.endr
3:
tbz x5,#2,4f
.rept 2
ldp x0,x1,[x22],#16
cmp x25,x0
b.hi 1f
add x2,x10,x0
staddb w19,[x2]
1:
cmp x25,x1
b.hi 2f
add x3,x10,x1
staddb w19,[x3]
2:
.endr
4:
tbz x5,#1,5f
ldp x0,x1,[x22],#16
cmp x25,x0
b.hi 1f
add x2,x10,x0
staddb w19,[x2]
1:
cmp x25,x1
b.hi 2f
add x3,x10,x1
staddb w19,[x3]
2:
5:
tbz x5,#0,6f
ldr x0,[x22],#8
cmp x25,x0
b.hi 6f
add x2,x10,x0
staddb w19,[x2]
6:
ret
//@end
//@seg x16_nobfi_bnd2_pre
seg_x16_nobfi_bnd2_pre:
sub x5,x23,x22
lsr x5,x5,#3
bic x21,x5,#15
add x21,x22,x21,lsl #3
cmp x22,x21
b.hs .Lnobfi_bnd2_tail
//@end
//@seg x16_nobfi_bnd2
seg_x16_nobfi_bnd2:
.Lnobfi_bnd2:
.rept 8
ldp x0,x1,[x22],#16
lsr x6,x0,x24
sub w6,w6,w25
cmp w6,w26
b.hi 1f
add x2,x10,x0
staddb w19,[x2]
1:
lsr x7,x1,x24
sub w7,w7,w25
cmp w7,w26
b.hi 2f
add x3,x10,x1
staddb w19,[x3]
2:
.endr
cmp x22,x21
b.lo .Lnobfi_bnd2
//@end
//@seg x16_nobfi_bnd2_tail
.Lnobfi_bnd2_tail:
tbz x5,#3,3f
.rept 4
ldp x0,x1,[x22],#16
lsr x6,x0,x24
sub w6,w6,w25
cmp w6,w26
b.hi 1f
add x2,x10,x0
staddb w19,[x2]
1:
lsr x7,x1,x24
sub w7,w7,w25
cmp w7,w26
b.hi 2f
add x3,x10,x1
staddb w19,[x3]
2:
.endr
3:
tbz x5,#2,4f
.rept 2
ldp x0,x1,[x22],#16
lsr x6,x0,x24
sub w6,w6,w25
cmp w6,w26
b.hi 1f
add x2,x10,x0
staddb w19,[x2]
1:
lsr x7,x1,x24
sub w7,w7,w25
cmp w7,w26
b.hi 2f
add x3,x10,x1
staddb w19,[x3]
2:
.endr
4:
tbz x5,#1,5f
ldp x0,x1,[x22],#16
lsr x6,x0,x24
sub w6,w6,w25
cmp w6,w26
b.hi 1f
add x2,x10,x0
staddb w19,[x2]
1:
lsr x7,x1,x24
sub w7,w7,w25
cmp w7,w26
b.hi 2f
add x3,x10,x1
staddb w19,[x3]
2:
5:
tbz x5,#0,6f
ldr x0,[x22],#8
lsr x6,x0,x24
sub w6,w6,w25
cmp w6,w26
b.hi 6f
add x2,x10,x0
staddb w19,[x2]
6:
ret
//@end
//@seg init8_u8
init8_u8:
.rept 8
ldp x0,x1,[x15],#16
stp x0,x1,[x20],#16
stp x0,x1,[x21],#16
stp x0,x1,[x22],#16
stp x0,x1,[x23],#16
stp x0,x1,[x24],#16
stp x0,x1,[x25],#16
stp x0,x1,[x26],#16
stp x0,x1,[x27],#16
.endr
cmp x15,x16
b.lo init8_u8
ret
//@end
//@seg scan8_u8
scan8_u8:
mov x14,#0
.rept 8
ldp x0,x1,[x20],#16
ldp x2,x3,[x21],#16
orr x0,x0,x2
orr x1,x1,x3
ldp x2,x3,[x22],#16
orr x0,x0,x2
orr x1,x1,x3
ldp x2,x3,[x23],#16
orr x0,x0,x2
orr x1,x1,x3
ldp x2,x3,[x24],#16
orr x0,x0,x2
orr x1,x1,x3
ldp x2,x3,[x25],#16
orr x0,x0,x2
orr x1,x1,x3
ldp x2,x3,[x26],#16
orr x0,x0,x2
orr x1,x1,x3
ldp x2,x3,[x27],#16
orr x0,x0,x2
orr x1,x1,x3
orr x0,x0,x1
orr x14,x14,x0
.endr
tst x14,x19
b.ne .Lscan8_hit
cmp x20,x28
b.lo scan8_u8
ret
.Lscan8_hit:
ret
//@end
//@seg b1_line
b1_line:
ldp x20,x27,[x1],#16
ldp x10,x11,[x1],#16
ldp x12,x13,[x2]
sbfx x14,x10,#0,#16
mul x15,x14,x12
sbfx x14,x10,#16,#16
madd x15,x14,x13,x15
ldp x12,x13,[x2,#16]
sbfx x14,x10,#32,#16
madd x15,x14,x12,x15
sbfx x14,x10,#48,#16
madd x15,x14,x13,x15
ldp x12,x13,[x2,#32]
sbfx x14,x11,#0,#16
madd x15,x14,x12,x15
sbfx x14,x11,#16,#16
madd x15,x14,x13,x15
add x15,x15,x7
udiv x16,x15,x5
msub x15,x16,x5,x15
ldr x12,[x2]
subs x16,x15,x12
add x17,x16,x5
csel x15,x16,x17,hs
lsl x21,x15,x6
lsl x22,x12,x6
b b1_row
//@end
//@seg b1_row
b1_row:
ldp x9,x10,[x20],#16
tbnz x10,#31,b1_switch
b1_row_cont:
add x21,x21,x22
subs x11,x21,x23
csel x21,x11,x21,hs
add x12,x24,x21
and x13,x10,#0xffff
lsr x14,x10,#32
ldr x15,[x12,x13]
ldr x16,[x12,x14]
subs x17,x16,x15
b.eq b1_next
tbz x17,#3,b1_even
ldr x0,[x15],#8
add x0,x9,x0
str x0,[x25],#8
b1_even:
ubfx x11,x17,#4,#6
adr x12,b1_copy_end
sub x12,x12,x11,lsl #4
br x12
//@end
//@seg b1_copy
b1_copy:
.rept 64
ldp x0,x1,[x15],#16
add x0,x9,x0
add x1,x9,x1
stp x0,x1,[x25],#16
.endr
b1_copy_end:
cmp x25,x26
b.hi b1_flush
b1_pass_room:
cmp x15,x16
b.lo b1_copy
b1_next:
cmp x20,x27
b.lo b1_row
ret
//@end
//@seg b1_flush
b1_flush:
sub x11,x25,x28
lsr x11,x11,#3
sub x11,x11,#2
stp x11,x29,[x28]
b1_handoff:
add x25,x28,#16
b b1_pass_room
//@end
//@seg b1_term
b1_term:
sub x11,x25,x28
lsr x11,x11,#3
sub x11,x11,#2
stp x11,x29,[x28]
b1_term_handoff:
add x25,x28,#16
ret
//@end
//@seg b1_switch
b1_switch:
stp x28,x25,[x19]
ldr x11,[x20],#8
add x19,x3,x11,lsl #4
ldp x28,x25,[x19]
add x26,x28,#3,lsl #12
add x26,x26,#3072
b b1_row_cont
//@end
//@seg b1_page
b1_page:
ldp x11,x6,[x28]
add x4,x28,#16
add x8,x4,x11,lsl #3
tbz x11,#0,b1_page_even
ldr x0,[x4],#8
staddb w6,[x0]
b1_page_even:
ubfx x12,x11,#1,#7
adr x13,b1_apply_end
sub x13,x13,x12,lsl #3
sub x13,x13,x12,lsl #2
br x13
//@end
//@seg b1_apply
b1_apply:
.rept 128
ldp x0,x1,[x4],#16
staddb w6,[x0]
staddb w6,[x1]
.endr
b1_apply_end:
cmp x4,x8
b.lo b1_apply
ret
//@end
//@seg delta_b1_copy
delta_b1_copy:
.rept 16
ldp x0,x1,[x15],#16
stp x0,x1,[x25],#16
.endr
delta_b1_copy_end:
cmp x15,x16
b.lo delta_b1_copy
ret
//@end
//@seg delta_b1_apply
delta_b1_apply:
.rept 32
ldp x0,x1,[x4],#16
add x5,x5,w0,uxtw
staddb w6,[x5]
add x5,x5,x0,lsr #32
staddb w6,[x5]
add x5,x5,w1,uxtw
staddb w6,[x5]
add x5,x5,x1,lsr #32
staddb w6,[x5]
.endr
delta_b1_apply_end:
cmp x4,x8
b.lo delta_b1_apply
ret
//@end
//@seg delta_b1_apply_tail
delta_b1_apply_tail:
cmp x4,x9
b.hs .Ldb1_done
.Ldb1_tloop:
ldr w0,[x4],#4
add x5,x5,w0,uxtw
staddb w6,[x5]
cmp x4,x9
b.lo .Ldb1_tloop
.Ldb1_done:
ret
//@end
//@seg mul4_sb
mul4_sb:
ldp x3,x4,[x0]
ldp x5,x6,[x0,#16]
ldp x7,x8,[x1]
ldp x9,x10,[x1,#16]
mul x11,x3,x7
umulh x12,x3,x7
mov x13,#0
str x11,[x2]
mov x11,x12
mov x12,#0
mul x14,x3,x8
umulh x15,x3,x8
adds x11,x11,x14
adcs x12,x12,x15
adc x13,x13,xzr
mul x14,x4,x7
umulh x15,x4,x7
adds x11,x11,x14
adcs x12,x12,x15
adc x13,x13,xzr
str x11,[x2,#8]
mov x11,x12
mov x12,x13
mov x13,#0
mul x14,x3,x9
umulh x15,x3,x9
adds x11,x11,x14
adcs x12,x12,x15
adc x13,x13,xzr
mul x14,x4,x8
umulh x15,x4,x8
adds x11,x11,x14
adcs x12,x12,x15
adc x13,x13,xzr
mul x14,x5,x7
umulh x15,x5,x7
adds x11,x11,x14
adcs x12,x12,x15
adc x13,x13,xzr
str x11,[x2,#16]
mov x11,x12
mov x12,x13
mov x13,#0
mul x14,x3,x10
umulh x15,x3,x10
adds x11,x11,x14
adcs x12,x12,x15
adc x13,x13,xzr
mul x14,x4,x9
umulh x15,x4,x9
adds x11,x11,x14
adcs x12,x12,x15
adc x13,x13,xzr
mul x14,x5,x8
umulh x15,x5,x8
adds x11,x11,x14
adcs x12,x12,x15
adc x13,x13,xzr
mul x14,x6,x7
umulh x15,x6,x7
adds x11,x11,x14
adcs x12,x12,x15
adc x13,x13,xzr
str x11,[x2,#24]
mov x11,x12
mov x12,x13
mov x13,#0
mul x14,x4,x10
umulh x15,x4,x10
adds x11,x11,x14
adcs x12,x12,x15
adc x13,x13,xzr
mul x14,x5,x9
umulh x15,x5,x9
adds x11,x11,x14
adcs x12,x12,x15
adc x13,x13,xzr
mul x14,x6,x8
umulh x15,x6,x8
adds x11,x11,x14
adcs x12,x12,x15
adc x13,x13,xzr
str x11,[x2,#32]
mov x11,x12
mov x12,x13
mov x13,#0
mul x14,x5,x10
umulh x15,x5,x10
adds x11,x11,x14
adcs x12,x12,x15
adc x13,x13,xzr
mul x14,x6,x9
umulh x15,x6,x9
adds x11,x11,x14
adcs x12,x12,x15
adc x13,x13,xzr
str x11,[x2,#40]
mov x11,x12
mov x12,x13
mul x14,x6,x10
umulh x15,x6,x10
adds x11,x11,x14
adc x12,x12,x15
stp x11,x12,[x2,#48]
ret
//@end
//@seg srt_odd_norm8_find
srt_odd_norm8_find:
ldr x3,[x0,x2,lsl #3]
cbnz x3,srt_odd_norm8_mid
add x2,x2,#1
cmp x2,#8
b.lo srt_odd_norm8_find
//@end
//@seg srt_odd_norm8_mid
srt_odd_norm8_mid:
rbit x4,x3
clz x4,x4
mov x5,#64
sub x5,x5,x4
mov x6,#0
//@end
//@seg srt_odd_norm8_shift
srt_odd_norm8_shift:
add x7,x6,x2
cmp x7,#8
b.hs srt_odd_zero_limb
ldr x8,[x0,x7,lsl #3]
lsr x8,x8,x4
add x9,x7,#1
cmp x9,#8
b.hs srt_odd_no_hi
ldr x10,[x0,x9,lsl #3]
lsl x10,x10,x5
tst x4,x4
csel x10,x10,xzr,ne
orr x8,x8,x10
srt_odd_no_hi:
b srt_odd_store
srt_odd_zero_limb:
mov x8,#0
srt_odd_store:
str x8,[x1,x6,lsl #3]
add x6,x6,#1
cmp x6,#8
b.lo srt_odd_norm8_shift
//@end
//@seg inv64_odd
inv64_odd:
mov x1,#1
mov x3,#2
mul x2,x0,x1
sub x2,x3,x2
mul x1,x1,x2
mul x2,x0,x1
sub x2,x3,x2
mul x1,x1,x2
mul x2,x0,x1
sub x2,x3,x2
mul x1,x1,x2
mul x2,x0,x1
sub x2,x3,x2
mul x1,x1,x2
mul x2,x0,x1
sub x2,x3,x2
mul x1,x1,x2
mul x2,x0,x1
sub x2,x3,x2
mul x1,x1,x2
ret
//@end
//@seg cios8_setup
cios8_setup:
stp xzr,xzr,[x4]
stp xzr,xzr,[x4,#16]
stp xzr,xzr,[x4,#32]
stp xzr,xzr,[x4,#48]
stp xzr,xzr,[x4,#64]
mov x5,#0
//@end
//@seg cios8_outer_pre
cios8_outer_pre:
ldr x6,[x0,x5,lsl #3]
mov x7,#0
mov x8,#0
//@end
//@seg cios8_inner_mul
cios8_inner_mul:
ldr x9,[x1,x8,lsl #3]
ldr x10,[x4,x8,lsl #3]
mul x11,x6,x9
umulh x12,x6,x9
adds x11,x11,x10
adc x12,x12,xzr
adds x11,x11,x7
adc x7,x12,xzr
str x11,[x4,x8,lsl #3]
add x8,x8,#1
cmp x8,#8
b.lo cios8_inner_mul
//@end
//@seg cios8_outer_mid
cios8_outer_mid:
ldr x10,[x4,#64]
adds x10,x10,x7
adc x13,xzr,xzr
str x10,[x4,#64]
str x13,[x4,#72]
ldr x14,[x4]
mul x15,x14,x3
ldr x9,[x2]
mul x11,x15,x9
umulh x12,x15,x9
adds x11,x11,x14
adc x7,x12,xzr
mov x8,#1
//@end
//@seg cios8_inner_red
cios8_inner_red:
ldr x9,[x2,x8,lsl #3]
ldr x10,[x4,x8,lsl #3]
mul x11,x15,x9
umulh x12,x15,x9
adds x11,x11,x10
adc x12,x12,xzr
adds x11,x11,x7
adc x7,x12,xzr
sub x16,x8,#1
str x11,[x4,x16,lsl #3]
add x8,x8,#1
cmp x8,#8
b.lo cios8_inner_red
//@end
//@seg cios8_outer_post
cios8_outer_post:
ldr x10,[x4,#64]
ldr x13,[x4,#72]
adds x10,x10,x7
adc x13,x13,xzr
str x10,[x4,#56]
str x13,[x4,#64]
str xzr,[x4,#72]
add x5,x5,#1
cmp x5,#8
b.lo cios8_outer_pre
//@end
//@seg cios8_csub_init
cios8_csub_init:
mov x8,#0
mov x7,#8
subs xzr,xzr,xzr
//@end
//@seg cios8_csub_loop
cios8_csub_loop:
ldr x9,[x4,x8,lsl #3]
ldr x10,[x2,x8,lsl #3]
sbcs x11,x9,x10
str x11,[x17,x8,lsl #3]
add x8,x8,#1
sub x7,x7,#1
cbnz x7,cios8_csub_loop
//@end
//@seg cios8_csel_init
cios8_csel_init:
ldr x9,[x4,#64]
sbcs xzr,x9,xzr
mov x8,#0
mov x7,#8
//@end
//@seg cios8_csel_loop
cios8_csel_loop:
ldr x9,[x4,x8,lsl #3]
ldr x10,[x17,x8,lsl #3]
csel x11,x10,x9,hs
str x11,[x17,x8,lsl #3]
add x8,x8,#1
sub x7,x7,#1
cbnz x7,cios8_csel_loop
//@end
//@seg bingcd8_io
bingcd8_io:
ldp x0,x1,[x25]
ldp x2,x3,[x25,#16]
ldp x4,x5,[x25,#32]
ldp x6,x7,[x25,#48]
ldp x8,x9,[x26]
ldp x10,x11,[x26,#16]
ldp x12,x13,[x26,#32]
ldp x14,x15,[x26,#48]
stp x0,x1,[x27]
stp x2,x3,[x27,#16]
stp x4,x5,[x27,#32]
stp x6,x7,[x27,#48]
ret
//@end
//@seg bingcd8_chk
bingcd8_loop:
orr x16,x8,x9
orr x17,x10,x11
orr x16,x16,x17
orr x17,x12,x13
orr x16,x16,x17
orr x17,x14,x15
orr x16,x16,x17
cbz x16,bingcd8_done
//@end
//@seg bingcd8_odd_sub
bingcd8_odd_sub:
tbz x8,#0,bingcd8_shr_v
subs x16,x8,x0
sbcs x17,x9,x1
sbcs x19,x10,x2
sbcs x20,x11,x3
sbcs x21,x12,x4
sbcs x22,x13,x5
sbcs x23,x14,x6
sbcs x24,x15,x7
b.hs bingcd8_shr_d
//@end
//@seg bingcd8_swap
bingcd8_swap:
subs x16,x0,x8
sbcs x17,x1,x9
sbcs x19,x2,x10
sbcs x20,x3,x11
sbcs x21,x4,x12
sbcs x22,x5,x13
sbcs x23,x6,x14
sbcs x24,x7,x15
mov x0,x8
mov x1,x9
mov x2,x10
mov x3,x11
mov x4,x12
mov x5,x13
mov x6,x14
mov x7,x15
//@end
//@seg bingcd8_shr_d
bingcd8_shr_d:
extr x8,x17,x16,#1
extr x9,x19,x17,#1
extr x10,x20,x19,#1
extr x11,x21,x20,#1
extr x12,x22,x21,#1
extr x13,x23,x22,#1
extr x14,x24,x23,#1
lsr x15,x24,#1
b bingcd8_loop
//@end
//@seg bingcd8_shr_v
bingcd8_shr_v:
extr x8,x9,x8,#1
extr x9,x10,x9,#1
extr x10,x11,x10,#1
extr x11,x12,x11,#1
extr x12,x13,x12,#1
extr x13,x14,x13,#1
extr x14,x15,x14,#1
lsr x15,x15,#1
b bingcd8_loop
bingcd8_done:
//@end
//@seg hensel_div8_outer
hensel_div8_outer:
ldr x10,[x21,x8,lsl #3]
mul x9,x10,x1
str x9,[x23,x8,lsl #3]
mov x11,#0
mov x12,#0
mov x10,#0
mov x14,#8
sub x14,x14,x8
//@end
//@seg hensel_div8_inner
hensel_div8_inner:
ldr x13,[x22,x10,lsl #3]
add x15,x8,x10
ldr x16,[x21,x15,lsl #3]
mul x17,x9,x13
umulh x7,x9,x13
adds x17,x17,x11
adc x11,x7,xzr
subs x16,x16,x17
cset x7,lo
subs x16,x16,x12
cset x17,lo
orr x12,x7,x17
str x16,[x21,x15,lsl #3]
add x10,x10,#1
cmp x10,x14
b.lo hensel_div8_inner
//@end
//@seg r4_4clamp_inner_tw
seg_r4_4clamp_inner_tw:
.Lr4_inner:
ldr x5,[x0]
ldr x6,[x1]
ldr x7,[x2]
ldr x8,[x3]
ldp x9,x10,[x4],#16
ldp x11,x12,[x4],#16
ldp x13,x14,[x4],#16
subs x15,x5,x20
csel x5,x15,x5,hs
umulh x15,x6,x10
mul x6,x6,x9
msub x6,x15,x19,x6
umulh x15,x7,x12
mul x7,x7,x11
msub x7,x15,x19,x7
umulh x15,x8,x14
mul x8,x8,x13
msub x8,x15,x19,x8
add x9,x5,x7
subs x15,x9,x20
csel x9,x15,x9,hs
add x10,x5,x20
sub x10,x10,x7
subs x15,x10,x20
csel x10,x15,x10,hs
add x11,x6,x8
subs x15,x11,x20
csel x11,x15,x11,hs
add x12,x6,x20
sub x12,x12,x8
umulh x15,x12,x23
mul x12,x12,x22
msub x12,x15,x19,x12
add x5,x9,x11
add x7,x9,x20
sub x7,x7,x11
add x6,x10,x12
add x8,x10,x20
sub x8,x8,x12
str x5,[x0],#8
str x6,[x1],#8
str x7,[x2],#8
str x8,[x3],#8
cmp x0,x24
b.lo .Lr4_inner
//@end
ret
//@seg r4_4clamp_hoist_tw_x2
seg_r4_4clamp_hoist_tw_x2:
.Lr4_hoist:
ldp x5,x16,[x0]
ldp x6,x17,[x1]
ldp x7,x25,[x2]
ldp x8,x26,[x3]
subs x15,x5,x20
csel x5,x15,x5,hs
umulh x15,x6,x10
mul x6,x6,x9
msub x6,x15,x19,x6
umulh x15,x7,x12
mul x7,x7,x11
msub x7,x15,x19,x7
umulh x15,x8,x14
mul x8,x8,x13
msub x8,x15,x19,x8
add x27,x5,x7
subs x15,x27,x20
csel x27,x15,x27,hs
add x28,x5,x20
sub x28,x28,x7
subs x15,x28,x20
csel x28,x15,x28,hs
add x29,x6,x8
subs x15,x29,x20
csel x29,x15,x29,hs
add x21,x6,x20
sub x21,x21,x8
umulh x15,x21,x23
mul x21,x21,x22
msub x21,x15,x19,x21
add x5,x27,x29
add x7,x27,x20
sub x7,x7,x29
add x6,x28,x21
add x8,x28,x20
sub x8,x8,x21
subs x15,x16,x20
csel x16,x15,x16,hs
umulh x15,x17,x10
mul x17,x17,x9
msub x17,x15,x19,x17
umulh x15,x25,x12
mul x25,x25,x11
msub x25,x15,x19,x25
umulh x15,x26,x14
mul x26,x26,x13
msub x26,x15,x19,x26
add x27,x16,x25
subs x15,x27,x20
csel x27,x15,x27,hs
add x28,x16,x20
sub x28,x28,x25
subs x15,x28,x20
csel x28,x15,x28,hs
add x29,x17,x26
subs x15,x29,x20
csel x29,x15,x29,hs
add x21,x17,x20
sub x21,x21,x26
umulh x15,x21,x23
mul x21,x21,x22
msub x21,x15,x19,x21
add x16,x27,x29
add x25,x27,x20
sub x25,x25,x29
add x17,x28,x21
add x26,x28,x20
sub x26,x26,x21
stp x5,x16,[x0],#16
stp x6,x17,[x1],#16
stp x7,x25,[x2],#16
stp x8,x26,[x3],#16
cmp x0,x24
b.lo .Lr4_hoist
//@end
ret
//@seg r3_shoup_bfly
seg_r3_shoup_bfly:
.Lr3_loop:
ldr x5,[x0]
ldr x6,[x1]
ldr x7,[x2]
ldp x8,x9,[x3],#16
ldp x10,x11,[x3],#16
subs x15,x5,x20
csel x5,x15,x5,hs
umulh x15,x6,x9
mul x6,x6,x8
msub x6,x15,x19,x6
umulh x15,x7,x11
mul x7,x7,x10
msub x7,x15,x19,x7
add x12,x6,x7
subs x15,x12,x20
csel x12,x15,x12,hs
add x13,x6,x20
sub x13,x13,x7
umulh x15,x12,x22
mul x14,x12,x21
msub x14,x15,x19,x14
umulh x15,x13,x25
mul x13,x13,x23
msub x13,x15,x19,x13
add x14,x5,x14
subs x15,x14,x20
csel x14,x15,x14,hs
add x5,x5,x12
add x6,x14,x13
add x7,x14,x20
sub x7,x7,x13
str x5,[x0],#8
str x6,[x1],#8
str x7,[x2],#8
cmp x0,x24
b.lo .Lr3_loop
//@end
ret
//@seg r5_shoup_bfly
seg_r5_shoup_bfly:
.Lr5_loop:
ldr x6,[x0]
ldr x7,[x1]
ldr x8,[x2]
ldr x9,[x3]
ldr x10,[x4]
ldp x11,x12,[x5],#16
ldp x13,x14,[x5],#16
subs x15,x6,x20
csel x6,x15,x6,hs
umulh x15,x7,x12
mul x7,x7,x11
msub x7,x15,x19,x7
umulh x15,x8,x14
mul x8,x8,x13
msub x8,x15,x19,x8
ldp x11,x12,[x5],#16
ldp x13,x14,[x5],#16
umulh x15,x9,x12
mul x9,x9,x11
msub x9,x15,x19,x9
umulh x15,x10,x14
mul x10,x10,x13
msub x10,x15,x19,x10
add x11,x7,x10
subs x15,x11,x20
csel x11,x15,x11,hs
add x12,x8,x9
subs x15,x12,x20
csel x12,x15,x12,hs
add x13,x7,x20
sub x13,x13,x10
subs x15,x13,x20
csel x13,x15,x13,hs
add x14,x8,x20
sub x14,x14,x9
subs x15,x14,x20
csel x14,x15,x14,hs
add x16,x11,x12
subs x15,x16,x20
csel x16,x15,x16,hs
add x17,x11,x20
sub x17,x17,x12
add x22,x13,x14
add x7,x6,x16
str x7,[x0],#8
ldp x25,x26,[x21,#0]
umulh x15,x16,x26
mul x16,x16,x25
msub x16,x15,x19,x16
ldp x25,x26,[x21,#16]
umulh x15,x17,x26
mul x17,x17,x25
msub x17,x15,x19,x17
add x16,x6,x16
subs x15,x16,x20
csel x16,x15,x16,hs
add x27,x16,x17
subs x15,x27,x20
csel x27,x15,x27,hs
add x28,x16,x20
sub x28,x28,x17
subs x15,x28,x20
csel x28,x15,x28,hs
ldp x25,x26,[x21,#32]
umulh x15,x22,x26
mul x22,x22,x25
msub x22,x15,x19,x22
ldp x25,x26,[x21,#48]
umulh x15,x13,x26
mul x13,x13,x25
msub x13,x15,x19,x13
ldp x25,x26,[x21,#64]
umulh x15,x14,x26
mul x14,x14,x25
msub x14,x15,x19,x14
add x29,x22,x13
subs x15,x29,x20
csel x29,x15,x29,hs
add x23,x22,x14
subs x15,x23,x20
csel x23,x15,x23,hs
add x7,x27,x29
add x10,x27,x20
sub x10,x10,x29
add x8,x28,x23
add x9,x28,x20
sub x9,x9,x23
str x7,[x1],#8
str x8,[x2],#8
str x9,[x3],#8
str x10,[x4],#8
cmp x0,x24
b.lo .Lr5_loop
//@end
ret
//@seg desc_zero16
seg_desc_zero16:
.Ldz_loop:
.rept 16
stp xzr,xzr,[x0],#16
.endr
cmp x0,x1
b.lo .Ldz_loop
//@end
ret
//@seg desc_scan8_u8
seg_desc_scan8_u8:
.Lds_loop:
.rept 4
ldp x2,x3,[x0],#16
and x4,x2,x19
and x5,x3,x19
orr x4,x4,x5
cbnz x4,.Lds_hit
.endr
cmp x0,x1
b.lo .Lds_loop
//@end
.Lds_hit:
ret
//@seg desc_sieve_hit
seg_desc_sieve_hit:
ldp x0,x6,[x2]
ldr x3,[x6,x5,lsl #3]
add x8,x17,x3,lsr #27
add x7,x0,x3,lsl #38
ldadd x1,x9,[x8]
str x7,[x9]
tbnz x9,#14,.Lds_fl
add x2,x2,x26
add x5,x5,x27
//@end
.Lds_fl:
ret
//@seg reg_bingcd8_io
seg_reg_bingcd8_io:
ldp x0,x1,[x28,#0]
ldp x2,x3,[x28,#16]
ldp x4,x5,[x28,#32]
ldp x6,x7,[x28,#48]
ldp x8,x9,[x29,#0]
ldp x10,x11,[x29,#16]
ldp x12,x13,[x29,#32]
ldp x14,x15,[x29,#48]
orr x16,x8,x9
orr x17,x10,x11
orr x16,x16,x12
orr x17,x17,x13
orr x16,x16,x14
orr x17,x17,x15
orr x16,x16,x17
cbz x16,.Lbg_done_u
//@end
//@seg reg_bingcd8_u_odd_chk
seg_reg_bingcd8_u_odd_chk:
.Lloop_u_odd:
tbz x8,#0,.Lbg_v_even
//@end
//@seg reg_bingcd8_u_odd_ge
seg_reg_bingcd8_u_odd_ge:
subs x19,x8,x0
sbcs x20,x9,x1
sbcs x21,x10,x2
sbcs x22,x11,x3
sbcs x23,x12,x4
sbcs x24,x13,x5
sbcs x25,x14,x6
sbcs x26,x15,x7
b.lo .Lbg_u_gt_v
extr x8,x20,x19,#1
extr x9,x21,x20,#1
extr x10,x22,x21,#1
extr x11,x23,x22,#1
extr x12,x24,x23,#1
extr x13,x25,x24,#1
extr x14,x26,x25,#1
lsr x15,x26,#1
orr x16,x8,x9
orr x17,x10,x11
orr x16,x16,x12
orr x17,x17,x13
orr x16,x16,x14
orr x17,x17,x15
orr x16,x16,x17
cbnz x16,.Lloop_u_odd
b .Lbg_done_u
//@end
//@seg reg_bingcd8_u_odd_lt
seg_reg_bingcd8_u_odd_lt:
.Lbg_u_gt_v:
subs x19,x0,x8
sbcs x20,x1,x9
sbcs x21,x2,x10
sbcs x22,x3,x11
sbcs x23,x4,x12
sbcs x24,x5,x13
sbcs x25,x6,x14
sbc x26,x7,x15
extr x0,x20,x19,#1
extr x1,x21,x20,#1
extr x2,x22,x21,#1
extr x3,x23,x22,#1
extr x4,x24,x23,#1
extr x5,x25,x24,#1
extr x6,x26,x25,#1
lsr x7,x26,#1
b .Lloop_v_odd
//@end
//@seg reg_bingcd8_v_even
seg_reg_bingcd8_v_even:
.Lbg_v_even:
extr x8,x9,x8,#1
extr x9,x10,x9,#1
extr x10,x11,x10,#1
extr x11,x12,x11,#1
extr x12,x13,x12,#1
extr x13,x14,x13,#1
extr x14,x15,x14,#1
lsr x15,x15,#1
b .Lloop_u_odd
//@end
//@seg reg_bingcd8_v_odd_chk
seg_reg_bingcd8_v_odd_chk:
.Lloop_v_odd:
tbz x0,#0,.Lbg_u_even
//@end
//@seg reg_bingcd8_v_odd_ge
seg_reg_bingcd8_v_odd_ge:
subs x19,x0,x8
sbcs x20,x1,x9
sbcs x21,x2,x10
sbcs x22,x3,x11
sbcs x23,x4,x12
sbcs x24,x5,x13
sbcs x25,x6,x14
sbcs x26,x7,x15
b.lo .Lbg_v_gt_u
extr x0,x20,x19,#1
extr x1,x21,x20,#1
extr x2,x22,x21,#1
extr x3,x23,x22,#1
extr x4,x24,x23,#1
extr x5,x25,x24,#1
extr x6,x26,x25,#1
lsr x7,x26,#1
orr x16,x0,x1
orr x17,x2,x3
orr x16,x16,x4
orr x17,x17,x5
orr x16,x16,x6
orr x17,x17,x7
orr x16,x16,x17
cbnz x16,.Lloop_v_odd
b .Lbg_done_v
//@end
//@seg reg_bingcd8_v_odd_lt
seg_reg_bingcd8_v_odd_lt:
.Lbg_v_gt_u:
subs x19,x8,x0
sbcs x20,x9,x1
sbcs x21,x10,x2
sbcs x22,x11,x3
sbcs x23,x12,x4
sbcs x24,x13,x5
sbcs x25,x14,x6
sbc x26,x15,x7
extr x8,x20,x19,#1
extr x9,x21,x20,#1
extr x10,x22,x21,#1
extr x11,x23,x22,#1
extr x12,x24,x23,#1
extr x13,x25,x24,#1
extr x14,x26,x25,#1
lsr x15,x26,#1
b .Lloop_u_odd
//@end
//@seg reg_bingcd8_u_even
seg_reg_bingcd8_u_even:
.Lbg_u_even:
extr x0,x1,x0,#1
extr x1,x2,x1,#1
extr x2,x3,x2,#1
extr x3,x4,x3,#1
extr x4,x5,x4,#1
extr x5,x6,x5,#1
extr x6,x7,x6,#1
lsr x7,x7,#1
b .Lloop_v_odd
//@end
//@seg reg_bingcd8_done_u
seg_reg_bingcd8_done_u:
.Lbg_done_u:
stp x0,x1,[x27,#0]
stp x2,x3,[x27,#16]
stp x4,x5,[x27,#32]
stp x6,x7,[x27,#48]
//@end
ret
//@seg reg_bingcd8_done_v
seg_reg_bingcd8_done_v:
.Lbg_done_v:
stp x8,x9,[x27,#0]
stp x10,x11,[x27,#16]
stp x12,x13,[x27,#32]
stp x14,x15,[x27,#48]
//@end
ret
.macro MONT5_ROUND
ldr x5,[x6],#8
mul x10,x0,x5
mul x11,x1,x5
mul x12,x2,x5
mul x13,x3,x5
mul x14,x4,x5
adds x15,x15,x10
adcs x16,x16,x11
adcs x17,x17,x12
adcs x27,x27,x13
adcs x19,x19,x14
adcs x20,x20,xzr
adc x21,xzr,xzr
umulh x10,x0,x5
umulh x11,x1,x5
umulh x12,x2,x5
umulh x13,x3,x5
umulh x14,x4,x5
adds x16,x16,x10
adcs x17,x17,x11
adcs x27,x27,x12
adcs x19,x19,x13
adcs x20,x20,x14
adc x21,x21,xzr
mul x9,x15,x7
mul x10,x9,x22
mul x11,x9,x23
mul x12,x9,x24
mul x13,x9,x25
mul x14,x9,x26
adds xzr,x15,x10
adcs x15,x16,x11
adcs x16,x17,x12
adcs x17,x27,x13
adcs x27,x19,x14
adcs x19,x20,xzr
adc x20,x21,xzr
umulh x10,x9,x22
umulh x11,x9,x23
umulh x12,x9,x24
umulh x13,x9,x25
umulh x14,x9,x26
adds x15,x15,x10
adcs x16,x16,x11
adcs x17,x17,x12
adcs x27,x27,x13
adcs x19,x19,x14
adc x20,x20,xzr
.endm
//@seg mont5_cios_prologue
seg_mont5_cios:
sub sp,sp,#80
stp x19,x20,[sp,#0]
stp x21,x22,[sp,#16]
stp x23,x24,[sp,#32]
stp x25,x26,[sp,#48]
str x27,[sp,#64]
ldp x0,x1,[x5,#0]
ldp x2,x3,[x5,#16]
ldr x4,[x5,#32]
ldp x22,x23,[x9,#0]
ldp x24,x25,[x9,#16]
ldr x26,[x9,#32]
mov x15,#0
mov x16,#0
mov x17,#0
mov x27,#0
mov x19,#0
mov x20,#0
//@end
//@seg mont5_cios_round
ldr x5,[x6],#8
mul x10,x0,x5
mul x11,x1,x5
mul x12,x2,x5
mul x13,x3,x5
mul x14,x4,x5
adds x15,x15,x10
adcs x16,x16,x11
adcs x17,x17,x12
adcs x27,x27,x13
adcs x19,x19,x14
adcs x20,x20,xzr
adc x21,xzr,xzr
umulh x10,x0,x5
umulh x11,x1,x5
umulh x12,x2,x5
umulh x13,x3,x5
umulh x14,x4,x5
adds x16,x16,x10
adcs x17,x17,x11
adcs x27,x27,x12
adcs x19,x19,x13
adcs x20,x20,x14
adc x21,x21,xzr
mul x9,x15,x7
mul x10,x9,x22
mul x11,x9,x23
mul x12,x9,x24
mul x13,x9,x25
mul x14,x9,x26
adds xzr,x15,x10
adcs x15,x16,x11
adcs x16,x17,x12
adcs x17,x27,x13
adcs x27,x19,x14
adcs x19,x20,xzr
adc x20,x21,xzr
umulh x10,x9,x22
umulh x11,x9,x23
umulh x12,x9,x24
umulh x13,x9,x25
umulh x14,x9,x26
adds x15,x15,x10
adcs x16,x16,x11
adcs x17,x17,x12
adcs x27,x27,x13
adcs x19,x19,x14
adc x20,x20,xzr
//@end
.rept 4
MONT5_ROUND
.endr
//@seg mont5_cios_epilogue
subs x10,x15,x22
sbcs x11,x16,x23
sbcs x12,x17,x24
sbcs x13,x27,x25
sbcs x14,x19,x26
sbcs xzr,x20,xzr
csel x15,x10,x15,hs
csel x16,x11,x16,hs
csel x17,x12,x17,hs
csel x27,x13,x27,hs
csel x19,x14,x19,hs
stp x15,x16,[x8,#0]
stp x17,x27,[x8,#16]
str x19,[x8,#32]
ldp x19,x20,[sp,#0]
ldp x21,x22,[sp,#16]
ldp x23,x24,[sp,#32]
ldp x25,x26,[sp,#48]
ldr x27,[sp,#64]
add sp,sp,#80
ret
//@end
//@seg mont2_cios
seg_mont2_cios:
mul x7,x0,x2
umulh x8,x0,x2
mul x9,x1,x2
umulh x10,x1,x2
adds x8,x8,x9
adc x9,x10,xzr
mul x11,x7,x6
mul x12,x11,x4
umulh x13,x11,x4
mul x14,x11,x5
umulh x15,x11,x5
adds xzr,x7,x12
adcs x7,x8,x14
adcs x8,x9,xzr
adc x9,xzr,xzr
adds x7,x7,x13
adcs x8,x8,x15
adc x9,x9,xzr
mul x10,x0,x3
umulh x11,x0,x3
mul x12,x1,x3
umulh x13,x1,x3
adds x7,x7,x10
adcs x8,x8,x12
adcs x9,x9,xzr
adc x10,xzr,xzr
adds x8,x8,x11
adcs x9,x9,x13
adc x10,x10,xzr
mul x11,x7,x6
mul x12,x11,x4
umulh x13,x11,x4
mul x14,x11,x5
umulh x15,x11,x5
adds xzr,x7,x12
adcs x7,x8,x14
adcs x8,x9,xzr
adc x9,x10,xzr
adds x7,x7,x13
adcs x8,x8,x15
adc x9,x9,xzr
subs x10,x7,x4
sbcs x11,x8,x5
sbcs xzr,x9,xzr
csel x0,x10,x7,hs
csel x1,x11,x8,hs
ret
//@end
//@seg bingcd_gen_loop
seg_bingcd_gen:
.Lbg_top:
mov x3,x1
mov x4,x2
mov x5,#0
.Lbg_vz:
ldr x6,[x3],#8
orr x5,x5,x6
sub x4,x4,#1
cbnz x4,.Lbg_vz
cbz x5,.Lbg_done
ldr x6,[x1]
tbnz x6,#0,.Lbg_both_odd
sub x4,x2,#1
mov x7,#0
.Lbg_vshr:
ldr x6,[x1,x4,lsl #3]
lsr x8,x6,#1
orr x8,x8,x7
str x8,[x1,x4,lsl #3]
lsl x7,x6,#63
subs x4,x4,#1
b.pl .Lbg_vshr
b .Lbg_top
.Lbg_both_odd:
sub x4,x2,#1
.Lbg_cmp:
ldr x6,[x0,x4,lsl #3]
ldr x7,[x1,x4,lsl #3]
cmp x6,x7
b.ne .Lbg_cmp_ne
subs x4,x4,#1
b.pl .Lbg_cmp
b .Lbg_done
.Lbg_cmp_ne:
b.hs .Lbg_sub_setup
mov x9,x0
mov x0,x1
mov x1,x9
.Lbg_sub_setup:
mov x3,x0
mov x4,x1
mov x5,x2
subs xzr,xzr,xzr
.Lbg_sub:
ldr x6,[x3]
ldr x7,[x4],#8
sbcs x8,x6,x7
str x8,[x3],#8
sub x5,x5,#1
cbnz x5,.Lbg_sub
mov x9,x0
mov x0,x1
mov x1,x9
sub x4,x2,#1
mov x7,#0
b .Lbg_vshr
.Lbg_done:
ret
//@end
//@seg cios_gen_setup
seg_cios_gen:
mov x10,x3
add x11,x4,#2
.Lcg_z:
str xzr,[x10],#8
subs x11,x11,#1
b.ne .Lcg_z
mov x7,#0
//@end
//@seg cios_gen_outer
.Lcg_outer:
ldr x8,[x0,x7,lsl #3]
mov x9,x1
mov x10,x3
mov x11,x4
mov x14,#0
//@end
//@seg cios_gen_mul_inner
.Lcg_mul:
ldr x12,[x9],#8
ldr x13,[x10]
mul x15,x8,x12
umulh x16,x8,x12
adds x15,x15,x14
adc x16,x16,xzr
adds x15,x15,x13
adc x14,x16,xzr
str x15,[x10],#8
subs x11,x11,#1
b.ne .Lcg_mul
//@end
//@seg cios_gen_mid
ldr x13,[x10]
adds x13,x13,x14
str x13,[x10]
adc x17,xzr,xzr
str x17,[x10,#8]
ldr x13,[x3]
mul x8,x13,x5
ldr x12,[x2]
mul x15,x8,x12
umulh x14,x8,x12
adds xzr,x13,x15
adc x14,x14,xzr
add x9,x2,#8
mov x10,x3
sub x11,x4,#1
//@end
//@seg cios_gen_red_inner
.Lcg_red:
ldr x12,[x9],#8
ldr x13,[x10,#8]
mul x15,x8,x12
umulh x16,x8,x12
adds x15,x15,x14
adc x16,x16,xzr
adds x15,x15,x13
adc x14,x16,xzr
str x15,[x10],#8
subs x11,x11,#1
b.ne .Lcg_red
//@end
//@seg cios_gen_outer_end
ldr x13,[x10,#8]
adds x13,x13,x14
str x13,[x10]
ldr x17,[x10,#16]
adc x17,x17,xzr
str x17,[x10,#8]
str xzr,[x10,#16]
add x7,x7,#1
cmp x7,x4
b.lo .Lcg_outer
//@end
//@seg cios_gen_epilogue
mov x9,x2
mov x10,x3
mov x15,x6
mov x11,x4
subs xzr,xzr,xzr
.Lcg_sub:
ldr x13,[x10],#8
ldr x12,[x9],#8
sbcs x16,x13,x12
str x16,[x15],#8
sub x11,x11,#1
cbnz x11,.Lcg_sub
ldr x13,[x10]
sbcs xzr,x13,xzr
b.hs .Lcg_done
mov x10,x3
mov x15,x6
mov x11,x4
.Lcg_cpy:
ldr x13,[x10],#8
str x13,[x15],#8
sub x11,x11,#1
cbnz x11,.Lcg_cpy
.Lcg_done:
ret
//@end
//@seg modinv_gen_loop
seg_modinv_gen:
.Lmi_top:
mov x6,x1
mov x7,x2
mov x8,#0
.Lmi_vz:
ldr x9,[x6],#8
orr x8,x8,x9
sub x7,x7,#1
cbnz x7,.Lmi_vz
cbz x8,.Lmi_done
ldr x9,[x1]
tbnz x9,#0,.Lmi_both_odd
.Lmi_halve:
sub x7,x2,#1
mov x10,#0
.Lmi_vshr:
ldr x9,[x1,x7,lsl #3]
lsr x11,x9,#1
orr x11,x11,x10
str x11,[x1,x7,lsl #3]
lsl x10,x9,#63
subs x7,x7,#1
b.pl .Lmi_vshr
ldr x9,[x4]
mov x10,#0
tbz x9,#0,.Lmi_x2shr_go
mov x6,x4
mov x7,x5
mov x8,x2
adds xzr,xzr,xzr
.Lmi_x2addN:
ldr x9,[x6]
ldr x11,[x7],#8
adcs x9,x9,x11
str x9,[x6],#8
sub x8,x8,#1
cbnz x8,.Lmi_x2addN
adc x10,xzr,xzr
lsl x10,x10,#63
.Lmi_x2shr_go:
sub x7,x2,#1
.Lmi_x2shr:
ldr x9,[x4,x7,lsl #3]
lsr x11,x9,#1
orr x11,x11,x10
str x11,[x4,x7,lsl #3]
lsl x10,x9,#63
subs x7,x7,#1
b.pl .Lmi_x2shr
b .Lmi_top
.Lmi_both_odd:
sub x7,x2,#1
.Lmi_cmp:
ldr x9,[x0,x7,lsl #3]
ldr x10,[x1,x7,lsl #3]
cmp x9,x10
b.ne .Lmi_cmp_ne
subs x7,x7,#1
b.pl .Lmi_cmp
b .Lmi_done
.Lmi_cmp_ne:
b.lo .Lmi_v_ge_u
mov x12,x0
mov x0,x1
mov x1,x12
mov x12,x3
mov x3,x4
mov x4,x12
.Lmi_v_ge_u:
mov x6,x1
mov x7,x0
mov x8,x2
subs xzr,xzr,xzr
.Lmi_vsub:
ldr x9,[x6]
ldr x10,[x7],#8
sbcs x11,x9,x10
str x11,[x6],#8
sub x8,x8,#1
cbnz x8,.Lmi_vsub
mov x6,x4
mov x7,x3
mov x8,x2
subs xzr,xzr,xzr
.Lmi_x2sub:
ldr x9,[x6]
ldr x10,[x7],#8
sbcs x11,x9,x10
str x11,[x6],#8
sub x8,x8,#1
cbnz x8,.Lmi_x2sub
b.hs .Lmi_halve
mov x6,x4
mov x7,x5
mov x8,x2
adds xzr,xzr,xzr
.Lmi_x2fix:
ldr x9,[x6]
ldr x10,[x7],#8
adcs x11,x9,x10
str x11,[x6],#8
sub x8,x8,#1
cbnz x8,.Lmi_x2fix
b .Lmi_halve
.Lmi_done:
ret
//@end

// Exact unchanged simple tabulation with paired digit/base loads.
// x0 key,x1 array73 bases,x2 output40 bytes;caller registers0..16 only.
.text
.p2align 2
.global fapi_tabulate
.type fapi_tabulate,%function
fapi_tabulate:
    add x15,x0,#256
    add x16,x1,#512
    mov x3,xzr
    mov x4,xzr
    mov x5,xzr
    mov x6,xzr
    mov x7,xzr
    ldp w8,w13,[x0,#0]
    ldp x9,x14,[x1,#0]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x0,#8]
    ldp x9,x14,[x1,#16]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x0,#16]
    ldp x9,x14,[x1,#32]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x0,#24]
    ldp x9,x14,[x1,#48]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x0,#32]
    ldp x9,x14,[x1,#64]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x0,#40]
    ldp x9,x14,[x1,#80]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x0,#48]
    ldp x9,x14,[x1,#96]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x0,#56]
    ldp x9,x14,[x1,#112]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x0,#64]
    ldp x9,x14,[x1,#128]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x0,#72]
    ldp x9,x14,[x1,#144]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x0,#80]
    ldp x9,x14,[x1,#160]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x0,#88]
    ldp x9,x14,[x1,#176]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x0,#96]
    ldp x9,x14,[x1,#192]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x0,#104]
    ldp x9,x14,[x1,#208]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x0,#112]
    ldp x9,x14,[x1,#224]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x0,#120]
    ldp x9,x14,[x1,#240]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x0,#128]
    ldp x9,x14,[x1,#256]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x0,#136]
    ldp x9,x14,[x1,#272]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x0,#144]
    ldp x9,x14,[x1,#288]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x0,#152]
    ldp x9,x14,[x1,#304]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x0,#160]
    ldp x9,x14,[x1,#320]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x0,#168]
    ldp x9,x14,[x1,#336]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x0,#176]
    ldp x9,x14,[x1,#352]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x0,#184]
    ldp x9,x14,[x1,#368]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x0,#192]
    ldp x9,x14,[x1,#384]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x0,#200]
    ldp x9,x14,[x1,#400]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x0,#208]
    ldp x9,x14,[x1,#416]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x0,#216]
    ldp x9,x14,[x1,#432]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x0,#224]
    ldp x9,x14,[x1,#448]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x0,#232]
    ldp x9,x14,[x1,#464]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x0,#240]
    ldp x9,x14,[x1,#480]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x0,#248]
    ldp x9,x14,[x1,#496]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x15,#0]
    ldp x9,x14,[x16,#0]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x15,#8]
    ldp x9,x14,[x16,#16]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x15,#16]
    ldp x9,x14,[x16,#32]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldp w8,w13,[x15,#24]
    ldp x9,x14,[x16,#48]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    add x10,x13,x13,lsl #2
    add x10,x14,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    ldr x8,[x0,#288]
    ldr x9,[x1,#576]
    add x10,x8,x8,lsl #2
    add x10,x9,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
    rev x3,x3
    rev x4,x4
    rev x5,x5
    rev x6,x6
    rev x7,x7
    stp x3,x4,[x2]
    stp x5,x6,[x2,#16]
    str x7,[x2,#32]
    ret
.size fapi_tabulate,.-fapi_tabulate

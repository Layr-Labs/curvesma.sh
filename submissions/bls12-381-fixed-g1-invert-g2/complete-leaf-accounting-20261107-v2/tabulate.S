// x0 exact key,x1 array of73 table bases,x2 output40-byte block payload.
// Caller-saved scalar leaf, all entries40 bytes; no host execution.
.text
.p2align 2
.global fapi_tabulate
.type fapi_tabulate,%function
fapi_tabulate:
    mov x3,xzr
    mov x4,xzr
    mov x5,xzr
    mov x6,xzr
    mov x7,xzr
    ldr w8,[x0,#0]
    ldr x9,[x1,#0]
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
    ldr w8,[x0,#4]
    ldr x9,[x1,#8]
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
    ldr w8,[x0,#8]
    ldr x9,[x1,#16]
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
    ldr w8,[x0,#12]
    ldr x9,[x1,#24]
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
    ldr w8,[x0,#16]
    ldr x9,[x1,#32]
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
    ldr w8,[x0,#20]
    ldr x9,[x1,#40]
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
    ldr w8,[x0,#24]
    ldr x9,[x1,#48]
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
    ldr w8,[x0,#28]
    ldr x9,[x1,#56]
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
    ldr w8,[x0,#32]
    ldr x9,[x1,#64]
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
    ldr w8,[x0,#36]
    ldr x9,[x1,#72]
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
    ldr w8,[x0,#40]
    ldr x9,[x1,#80]
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
    ldr w8,[x0,#44]
    ldr x9,[x1,#88]
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
    ldr w8,[x0,#48]
    ldr x9,[x1,#96]
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
    ldr w8,[x0,#52]
    ldr x9,[x1,#104]
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
    ldr w8,[x0,#56]
    ldr x9,[x1,#112]
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
    ldr w8,[x0,#60]
    ldr x9,[x1,#120]
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
    ldr w8,[x0,#64]
    ldr x9,[x1,#128]
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
    ldr w8,[x0,#68]
    ldr x9,[x1,#136]
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
    ldr w8,[x0,#72]
    ldr x9,[x1,#144]
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
    ldr w8,[x0,#76]
    ldr x9,[x1,#152]
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
    ldr w8,[x0,#80]
    ldr x9,[x1,#160]
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
    ldr w8,[x0,#84]
    ldr x9,[x1,#168]
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
    ldr w8,[x0,#88]
    ldr x9,[x1,#176]
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
    ldr w8,[x0,#92]
    ldr x9,[x1,#184]
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
    ldr w8,[x0,#96]
    ldr x9,[x1,#192]
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
    ldr w8,[x0,#100]
    ldr x9,[x1,#200]
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
    ldr w8,[x0,#104]
    ldr x9,[x1,#208]
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
    ldr w8,[x0,#108]
    ldr x9,[x1,#216]
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
    ldr w8,[x0,#112]
    ldr x9,[x1,#224]
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
    ldr w8,[x0,#116]
    ldr x9,[x1,#232]
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
    ldr w8,[x0,#120]
    ldr x9,[x1,#240]
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
    ldr w8,[x0,#124]
    ldr x9,[x1,#248]
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
    ldr w8,[x0,#128]
    ldr x9,[x1,#256]
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
    ldr w8,[x0,#132]
    ldr x9,[x1,#264]
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
    ldr w8,[x0,#136]
    ldr x9,[x1,#272]
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
    ldr w8,[x0,#140]
    ldr x9,[x1,#280]
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
    ldr w8,[x0,#144]
    ldr x9,[x1,#288]
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
    ldr w8,[x0,#148]
    ldr x9,[x1,#296]
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
    ldr w8,[x0,#152]
    ldr x9,[x1,#304]
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
    ldr w8,[x0,#156]
    ldr x9,[x1,#312]
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
    ldr w8,[x0,#160]
    ldr x9,[x1,#320]
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
    ldr w8,[x0,#164]
    ldr x9,[x1,#328]
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
    ldr w8,[x0,#168]
    ldr x9,[x1,#336]
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
    ldr w8,[x0,#172]
    ldr x9,[x1,#344]
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
    ldr w8,[x0,#176]
    ldr x9,[x1,#352]
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
    ldr w8,[x0,#180]
    ldr x9,[x1,#360]
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
    ldr w8,[x0,#184]
    ldr x9,[x1,#368]
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
    ldr w8,[x0,#188]
    ldr x9,[x1,#376]
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
    ldr w8,[x0,#192]
    ldr x9,[x1,#384]
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
    ldr w8,[x0,#196]
    ldr x9,[x1,#392]
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
    ldr w8,[x0,#200]
    ldr x9,[x1,#400]
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
    ldr w8,[x0,#204]
    ldr x9,[x1,#408]
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
    ldr w8,[x0,#208]
    ldr x9,[x1,#416]
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
    ldr w8,[x0,#212]
    ldr x9,[x1,#424]
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
    ldr w8,[x0,#216]
    ldr x9,[x1,#432]
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
    ldr w8,[x0,#220]
    ldr x9,[x1,#440]
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
    ldr w8,[x0,#224]
    ldr x9,[x1,#448]
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
    ldr w8,[x0,#228]
    ldr x9,[x1,#456]
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
    ldr w8,[x0,#232]
    ldr x9,[x1,#464]
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
    ldr w8,[x0,#236]
    ldr x9,[x1,#472]
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
    ldr w8,[x0,#240]
    ldr x9,[x1,#480]
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
    ldr w8,[x0,#244]
    ldr x9,[x1,#488]
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
    ldr w8,[x0,#248]
    ldr x9,[x1,#496]
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
    ldr w8,[x0,#252]
    ldr x9,[x1,#504]
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
    ldr w8,[x0,#256]
    ldr x9,[x1,#512]
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
    ldr w8,[x0,#260]
    ldr x9,[x1,#520]
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
    ldr w8,[x0,#264]
    ldr x9,[x1,#528]
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
    ldr w8,[x0,#268]
    ldr x9,[x1,#536]
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
    ldr w8,[x0,#272]
    ldr x9,[x1,#544]
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
    ldr w8,[x0,#276]
    ldr x9,[x1,#552]
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
    ldr w8,[x0,#280]
    ldr x9,[x1,#560]
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
    ldr w8,[x0,#284]
    ldr x9,[x1,#568]
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
    str x3,[x2,#0]
    rev x4,x4
    str x4,[x2,#8]
    rev x5,x5
    str x5,[x2,#16]
    rev x6,x6
    str x6,[x2,#24]
    rev x7,x7
    str x7,[x2,#32]
    ret
.size fapi_tabulate,.-fapi_tabulate

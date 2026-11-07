// Compact compile-time form of892-instruction tabulation leaf.
.text
.p2align 2
.global fapi_tabulate
.type fapi_tabulate,%function
.macro TAB_STEP keyoff,baseoff,wide
.if \wide
    ldr x8,[x0,#\keyoff]
.else
    ldr w8,[x0,#\keyoff]
.endif
    ldr x9,[x1,#\baseoff]
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
.endm
fapi_tabulate:
    mov x3,xzr
    mov x4,xzr
    mov x5,xzr
    mov x6,xzr
    mov x7,xzr
.set tab_i,0
.rept 72
    TAB_STEP (4*tab_i),(8*tab_i),0
.set tab_i,tab_i+1
.endr
    TAB_STEP 288,576,1
    rev x3,x3
    str x3,[x2]
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

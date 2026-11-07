// Static scalar field-copy/negation helpers,neverhost-executed.
.text
.p2align 2
.global fapi_copy48
.type fapi_copy48,%function
fapi_copy48:
    ldp x2,x3,[x0]
    stp x2,x3,[x1]
    ldp x2,x3,[x0,#16]
    stp x2,x3,[x1,#16]
    ldp x2,x3,[x0,#32]
    stp x2,x3,[x1,#32]
    ret
.size fapi_copy48,.-fapi_copy48
// x0=input,x1=output;tailcallsub(zero,input,output).
.p2align 2
.global fapi_fp_neg
.type fapi_fp_neg,%function
fapi_fp_neg:
    mov x2,x1
    mov x1,x0
    adrp x0,fapi_budget_zero_fp
    add x0,x0,:lo12:fapi_budget_zero_fp
    b fapi_fp_sub
.size fapi_fp_neg,.-fapi_fp_neg
.section .rodata
.p2align 4
fapi_budget_zero_fp:
.rept 6
    .quad 0
.endr

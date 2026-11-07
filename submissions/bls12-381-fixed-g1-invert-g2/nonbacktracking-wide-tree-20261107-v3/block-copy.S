// Static fixed-size scalar block helpers;neverhost-executed.
// Copy: x0source,x1destination;zero:x0destination. Callerregistersonly.
.text
.p2align 2
.global fapi_copy8192
.type fapi_copy8192,%function
fapi_copy8192:
.rept 512
    ldp x2,x3,[x0],#16
    stp x2,x3,[x1],#16
.endr
    ret
.size fapi_copy8192,.-fapi_copy8192
.p2align 2
.global fapi_zero8192
.type fapi_zero8192,%function
fapi_zero8192:
.rept 512
    stp xzr,xzr,[x0],#16
.endr
    ret
.size fapi_zero8192,.-fapi_zero8192
.p2align 2
.global fapi_zero512
.type fapi_zero512,%function
fapi_zero512:
.rept 32
    stp xzr,xzr,[x0],#16
.endr
    ret
.size fapi_zero512,.-fapi_zero512

// Static evidence only;neverhost-executed.
// x0:12 canonicalFp2 candidates,96bytes each,inFrobeniusk order0..11.
// x1:selectedFp2 destination96bytes,x2:destination64-bitk.
// Order:realnumeric384-bitcomponentthenimaginarynumeric384-bitcomponent.
// EachcomponentiscanonicalnativeMontgomery,little-endian6x64limbs.
// x0returnsk;callerregistersonly,nointernalcalls/frame,x18untouched.
.text
.p2align 2
.global fapi_canonical_min
.type fapi_canonical_min,%function
.macro COMPARE_CAND k
    add x1,x8,#(96*\k)
    mov x7,#\k
    ldp x2,x3,[x0,#48]
    ldp x4,x5,[x1,#48]
    subs xzr,x2,x4
    sbcs xzr,x3,x5
    ldp x2,x3,[x0,#64]
    ldp x4,x5,[x1,#64]
    sbcs xzr,x2,x4
    sbcs xzr,x3,x5
    ldp x2,x3,[x0,#80]
    ldp x4,x5,[x1,#80]
    sbcs xzr,x2,x4
    sbcs xzr,x3,x5
    ldp x2,x3,[x0,#0]
    ldp x4,x5,[x1,#0]
    sbcs xzr,x2,x4
    sbcs xzr,x3,x5
    ldp x2,x3,[x0,#16]
    ldp x4,x5,[x1,#16]
    sbcs xzr,x2,x4
    sbcs xzr,x3,x5
    ldp x2,x3,[x0,#32]
    ldp x4,x5,[x1,#32]
    sbcs xzr,x2,x4
    sbcs xzr,x3,x5
    csel x0,x1,x0,hs
    csel x6,x7,x6,hs
.endm
fapi_canonical_min:
    mov x8,x0
    mov x9,x1
    mov x10,x2
    mov x6,xzr
    COMPARE_CAND 1
    COMPARE_CAND 2
    COMPARE_CAND 3
    COMPARE_CAND 4
    COMPARE_CAND 5
    COMPARE_CAND 6
    COMPARE_CAND 7
    COMPARE_CAND 8
    COMPARE_CAND 9
    COMPARE_CAND 10
    COMPARE_CAND 11
    ldp x2,x3,[x0,#0]
    stp x2,x3,[x9,#0]
    ldp x2,x3,[x0,#16]
    stp x2,x3,[x9,#16]
    ldp x2,x3,[x0,#32]
    stp x2,x3,[x9,#32]
    ldp x2,x3,[x0,#48]
    stp x2,x3,[x9,#48]
    ldp x2,x3,[x0,#64]
    stp x2,x3,[x9,#64]
    ldp x2,x3,[x0,#80]
    stp x2,x3,[x9,#80]
    str x6,[x10]
    mov x0,x6
    ret
.size fapi_canonical_min,.-fapi_canonical_min

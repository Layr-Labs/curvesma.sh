// Static evidence only;neverhost-executed.
// x0:aligned24-byteunsigned192-bitremainingcounter,x1:smallamount0..4.
// Returnx0=1onsuccess,0onunderflow;nootherordinarycontextchanged.
.text
.p2align 2
.global fapi_consume_products
.type fapi_consume_products,%function
fapi_consume_products:
    ldp x2,x3,[x0]
    ldr x4,[x0,#16]
    subs x2,x2,x1
    sbcs x3,x3,xzr
    sbcs x4,x4,xzr
    b.lo .Lproducts_exhausted
    stp x2,x3,[x0]
    str x4,[x0,#16]
    mov x0,#1
    ret
.Lproducts_exhausted:
    mov x0,#0
    ret
.size fapi_consume_products,.-fapi_consume_products

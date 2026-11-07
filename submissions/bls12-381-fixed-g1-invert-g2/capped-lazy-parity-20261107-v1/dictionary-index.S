// Static analytic evidence only;neverhost-executed. No SIMD,x18 untouched.
// x0maincontext,x1currentrawrecord; x27windowbase preserved by cold parent.
// ctx24 points to initialized control(next_index128,remaining_visits128).
// Return x0=0inserted,1Umatch,2visitcap. On1:x1/x2=matchedp,x3=nativepayload.
// Payload p equals original visit time+1,so its exact A is p; no A field is stored.
.text
.p2align 2
.global fapi_trace_dictionary
.type fapi_trace_dictionary,%function
fapi_trace_dictionary:
    stp x19,x20,[sp,#-96]!
    stp x21,x22,[sp,#16]
    stp x23,x24,[sp,#32]
    stp x25,x26,[sp,#48]
    str x30,[sp,#64]
    mov x19,x0
    mov x20,x1
    ldr x21,[x19,#24]
    add x22,x19,#256
    add x23,x19,#64
    ldp x0,x1,[x19,#496]
    mov x2,x27
    bl fapi_dict_select_head
    ldp x24,x25,[x0]
    stp x24,x25,[sp,#80]
.Ldict_chain:
    orr x0,x24,x25
    cbz x0,.Ldict_insert
    ldp x0,x1,[x21,#16]
    orr x2,x0,x1
    cbz x2,.Ldict_cap
    subs x0,x0,#1
    sbc x1,x1,xzr
    stp x0,x1,[x21,#16]
    mov x0,x24
    mov x1,x25
    mov x2,x27
    bl fapi_dict_select_payload
    mov x26,x0
    ldp x0,x1,[x26,#384]
    ldp x2,x3,[x22,#0]
    eor x0,x0,x2
    eor x1,x1,x3
    orr x0,x0,x1
    ldp x1,x2,[x26,#400]
    ldp x3,x4,[x22,#16]
    eor x1,x1,x3
    eor x2,x2,x4
    orr x0,x0,x1
    orr x0,x0,x2
    ldr x1,[x26,#416]
    ldr x2,[x22,#32]
    eor x1,x1,x2
    orr x0,x0,x1
    cbz x0,.Ldict_match
    ldp x24,x25,[x26,#432]
    b .Ldict_chain
.Ldict_match:
    mov x1,x24
    mov x2,x25
    mov x3,x26
    mov x0,#1
    b .Ldict_return
.Ldict_cap:
    mov x0,#2
    b .Ldict_return
.Ldict_insert:
    ldp x24,x25,[x21,#0]
    adds x0,x24,#1
    adc x1,x25,xzr
    stp x0,x1,[x21,#0]
    mov x0,x24
    mov x1,x25
    mov x2,x27
    bl fapi_dict_select_payload
    mov x26,x0
// Zero every byte before publication;all32paired offsets<=496.
    stp xzr,xzr,[x26,#0]
    stp xzr,xzr,[x26,#16]
    stp xzr,xzr,[x26,#32]
    stp xzr,xzr,[x26,#48]
    stp xzr,xzr,[x26,#64]
    stp xzr,xzr,[x26,#80]
    stp xzr,xzr,[x26,#96]
    stp xzr,xzr,[x26,#112]
    stp xzr,xzr,[x26,#128]
    stp xzr,xzr,[x26,#144]
    stp xzr,xzr,[x26,#160]
    stp xzr,xzr,[x26,#176]
    stp xzr,xzr,[x26,#192]
    stp xzr,xzr,[x26,#208]
    stp xzr,xzr,[x26,#224]
    stp xzr,xzr,[x26,#240]
    stp xzr,xzr,[x26,#256]
    stp xzr,xzr,[x26,#272]
    stp xzr,xzr,[x26,#288]
    stp xzr,xzr,[x26,#304]
    stp xzr,xzr,[x26,#320]
    stp xzr,xzr,[x26,#336]
    stp xzr,xzr,[x26,#352]
    stp xzr,xzr,[x26,#368]
    stp xzr,xzr,[x26,#384]
    stp xzr,xzr,[x26,#400]
    stp xzr,xzr,[x26,#416]
    stp xzr,xzr,[x26,#432]
    stp xzr,xzr,[x26,#448]
    stp xzr,xzr,[x26,#464]
    stp xzr,xzr,[x26,#480]
    stp xzr,xzr,[x26,#496]
// Exact192-byte canonical key,all24words.
    ldp x0,x1,[x23,#0]
    stp x0,x1,[x26,#0]
    ldp x0,x1,[x23,#16]
    stp x0,x1,[x26,#16]
    ldp x0,x1,[x23,#32]
    stp x0,x1,[x26,#32]
    ldp x0,x1,[x23,#48]
    stp x0,x1,[x26,#48]
    ldp x0,x1,[x23,#64]
    stp x0,x1,[x26,#64]
    ldp x0,x1,[x23,#80]
    stp x0,x1,[x26,#80]
    ldp x0,x1,[x23,#96]
    stp x0,x1,[x26,#96]
    ldp x0,x1,[x23,#112]
    stp x0,x1,[x26,#112]
    ldp x0,x1,[x23,#128]
    stp x0,x1,[x26,#128]
    ldp x0,x1,[x23,#144]
    stp x0,x1,[x26,#144]
    ldp x0,x1,[x23,#160]
    stp x0,x1,[x26,#160]
    ldp x0,x1,[x23,#176]
    stp x0,x1,[x26,#176]
// B32,k8,U40 andoldheadnext128; A304..335 stays zero padding.
    add x9,x20,#608
    ldp x0,x1,[x9,#0]
    stp x0,x1,[x26,#336]
    ldp x0,x1,[x9,#16]
    stp x0,x1,[x26,#352]
    ldr x0,[x20,#944]
    str x0,[x26,#368]
    ldp x0,x1,[x22,#0]
    stp x0,x1,[x26,#384]
    ldp x0,x1,[x22,#16]
    stp x0,x1,[x26,#400]
    ldr x0,[x22,#32]
    str x0,[x26,#416]
    ldp x0,x1,[sp,#80]
    stp x0,x1,[x26,#432]
// Selector'sfirstDSBfinishesallpayloadwritesbeforeheadpublication.
    ldp x0,x1,[x19,#496]
    mov x2,x27
    bl fapi_dict_select_head
    stp x24,x25,[x0]
    mov x0,xzr
    b .Ldict_return
.Ldict_return:
    ldr x30,[sp,#64]
    ldp x21,x22,[sp,#16]
    ldp x23,x24,[sp,#32]
    ldp x25,x26,[sp,#48]
    ldp x19,x20,[sp],#96
    ret
.size fapi_trace_dictionary,.-fapi_trace_dictionary

// x0bucketlo,x1buckethi,x2windowbase -> x0nativeheadpointer.
// Allothermodifiedregistersarecallerscratch;oneleafRET.
.type fapi_dict_select_head,%function
fapi_dict_select_head:
    and x9,x0,#0xfffffffff
    extr x10,x1,x0,#36
    lsr x11,x1,#36
    orr x11,x11,#0x80000000
    adrp x12,fapi_dictionary_bank_pointer
    ldr x12,[x12,:lo12:fapi_dictionary_bank_pointer]
    dsb sy
    str x10,[x12,#0]
    str x11,[x12,#8]
    str xzr,[x12,#16]
    dsb sy
    add x0,x2,x9,lsl #4
    ret
.size fapi_dict_select_head,.-fapi_dict_select_head

// x0payloadlo,x1payloadhi,x2windowbase -> x0nativepayloadpointer.
.type fapi_dict_select_payload,%function
fapi_dict_select_payload:
    and x9,x0,#0x7fffffff
    extr x10,x1,x0,#31
    lsr x11,x1,#31
    adrp x12,fapi_dictionary_bank_pointer
    ldr x12,[x12,:lo12:fapi_dictionary_bank_pointer]
    dsb sy
    str x10,[x12,#0]
    str x11,[x12,#8]
    str xzr,[x12,#16]
    dsb sy
    add x0,x2,x9,lsl #9
    ret
.size fapi_dict_select_payload,.-fapi_dict_select_payload
.section .data
.p2align 3
// ColdinitializedtotheSAMEfiniteMMIOselectortripleasthejumpselector.
fapi_dictionary_bank_pointer:
    .quad 0

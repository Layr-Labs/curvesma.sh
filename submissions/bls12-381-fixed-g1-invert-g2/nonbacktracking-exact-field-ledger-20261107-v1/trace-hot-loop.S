// Static scalar evidence only;neverhost-executed.
// Cold entry savesoriginalLR/calleecontext,thenBRANCHEShere;terminalcoldhandlersrestore.
// x19current960record,x20alternate,x21jumpwindowpointer,x22:x23stepcount,
// x24:x25previousj,x26previouse,x27windowbase2^45,x28fixedsharedcontext.
// ctx:label(jlo,jhi,e)0..23,K128@32,tablebasesptr@48,key192@64,
// Ublock64@256,prefixCV32@320,mutableCV/output32@352,bucket128@496,
// localretrycounter@480. Otherretry/dictionary/coldcontextusesseparatespace.
// Allhelperspreservex19..x28;theirinternalwork/ABIareseparatelycharged.
.text
.p2align 2
.global fapi_trace_hot_loop
.type fapi_trace_hot_loop,%function
fapi_trace_hot_loop:
.Ltrace_tick:
    ldp x9,x10,[x28,#32]
    subs xzr,x22,x9
    sbcs xzr,x23,x10
    b.hs fapi_trace_cold_exhausted
.Ltrace_candidate:
    ldp x0,x1,[x28,#0]
    ldr x2,[x28,#16]
    cmp x0,x24
    ccmp x1,x25,#0,eq
    b.ne .Ltrace_accept
    add x9,x26,#6
    cmp x9,#12
    sub x10,x9,#12
    csel x9,x10,x9,hs
    cmp x2,x9
    b.eq .Ltrace_retry
.Ltrace_accept:
    mov x24,x0
    mov x25,x1
    mov x26,x2
    mov x3,x27
    bl fapi_trace_select_jump
    mov x21,x0
    adds x22,x22,#1
    adc x23,x23,xzr
    mov x0,x20
    mov x1,x19
    mov x2,x21
    bl blst_fp12_mul
    mov x9,x19
    mov x19,x20
    mov x20,x9
    mov x0,x20
    mov x1,x19
    mov x2,x21
    bl fapi_scalar_update
    mov x0,x19
    add x1,x28,#64
    bl fapi_trace_canonical
    cbz x0,fapi_trace_cold_failed
    add x0,x28,#64
    ldr x1,[x28,#48]
    add x2,x28,#256
    bl fapi_tabulate
    ldp x0,x1,[x28,#320]
    stp x0,x1,[x28,#352]
    ldp x0,x1,[x28,#336]
    stp x0,x1,[x28,#368]
    add x0,x28,#352
    add x1,x28,#256
    bl fapi_blake3_final
// Fourlittle-endian64-bitdigestwords. Low76j;upper180V.
    ldp x9,x10,[x28,#352]
    ldp x11,x12,[x28,#368]
    and x13,x10,#4095
    stp x9,x13,[x28,#0]
    extr x9,x11,x10,#12
    extr x10,x12,x11,#12
    lsr x11,x12,#12
    stp x9,x10,[x28,#496]
// Exactfloor(12V/2^180),threeordinarylimbs;no256-bitunitoperation.
    mov x16,#12
    umulh x12,x9,x16
    mul x13,x10,x16
    umulh x14,x10,x16
    adds x13,x13,x12
    adc x14,x14,xzr
    mul x11,x11,x16
    add x11,x11,x14
    lsr x11,x11,#52
    ldr x12,[x19,#944]
    subs x11,x11,x12
    add x13,x11,#12
    csel x11,x13,x11,lo
    str x11,[x28,#16]
    str xzr,[x28,#480]
    mov x0,x28
    mov x1,x19
    bl fapi_trace_dictionary
// status0new,1Umatch,2dictionarycapfailure. ExactSconfirmation/exitiscold.
    cbnz x0,fapi_trace_cold_dictionary_exit
    b .Ltrace_tick
.Ltrace_retry:
// Allretrybody+interface+rechecksareintheglobalextra-attemptbudget,notfree.
    mov x0,x28
    mov x1,x19
    mov x2,x24
    mov x3,x25
    mov x4,x26
    bl fapi_trace_retry_label
    cbz x0,fapi_trace_cold_failed
    b .Ltrace_candidate
.size fapi_trace_hot_loop,.-fapi_trace_hot_loop
// j128isactually76bits,e0..11,windowbaseinx3;returns640-bytejumprecordpointer.
.global fapi_trace_select_jump
.type fapi_trace_select_jump,%function
fapi_trace_select_jump:
    and x9,x0,#0x7ffffff
    extr x10,x1,x0,#27
    mov x11,#0x100000000
    lsl x9,x9,#2
    add x9,x9,x9,lsl #1
    add x9,x9,x2
    lsl x13,x9,#9
    add x13,x13,x9,lsl #7
    adrp x12,fapi_bank_register_pointer
    ldr x12,[x12,:lo12:fapi_bank_register_pointer]
    dsb sy
    str x10,[x12,#0]
    str x11,[x12,#8]
    str xzr,[x12,#16]
    dsb sy
    add x0,x3,x13
    ret
.size fapi_trace_select_jump,.-fapi_trace_select_jump
.section .data
.p2align 3
// Coldinitializedtothemodel'sMMIOregistertriple;notaninstalledVMconfiguration.
fapi_bank_register_pointer:
    .quad 0

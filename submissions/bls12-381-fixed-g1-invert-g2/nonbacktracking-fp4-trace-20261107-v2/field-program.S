// Static evidence only;neverhost-executed.
// x0:positive-length cold-built program,x1:positive descriptor count.
// Each32-byte descriptor is(function,arg0,arg1,arg2),all ordinary64-bitwords.
// Function signature/order isencoded perdescriptor,matching actualhelperABI.
// Fixedaddresses refer only to shared initialized operands/constants/scratch.
.text
.p2align 2
.global fapi_run_field_program
.type fapi_run_field_program,%function
fapi_run_field_program:
    stp x19,x20,[sp,#-32]!
    str x30,[sp,#16]
    mov x19,x0
    mov x20,x1
.Lfield_program_next:
    ldp x16,x0,[x19],#16
    ldp x1,x2,[x19],#16
    blr x16
    subs x20,x20,#1
    b.ne .Lfield_program_next
    ldr x30,[sp,#16]
    ldp x19,x20,[sp],#32
    ret
.size fapi_run_field_program,.-fapi_run_field_program

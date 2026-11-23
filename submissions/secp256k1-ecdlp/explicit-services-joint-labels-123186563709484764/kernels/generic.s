// Card54 generic scalar source template. CONST selectors are replaced at cold
// assembly by four literal MOV/MOVK instructions from the role/PE map:
// 16 callerCTRL (receiver uses PE15),32 REMOTE,40 ACK,48 INBOX,64 HANDLERS.
// There is deliberately no executable fallback for an unbound CONST.
// Roles: G0,R1,H2,B3. Each PE has route.lo/hi,ACK,self.lo/hi,selector,READY
// at offsets0,8,16,24,32,40,48 of its own64-byte control slot.
// The shared REMOTE and ACK aliases select that PE's route registers through
// fixed initiator-interface wiring; scalar publication and ACK remain explicit.
.text
.globl generic_send
generic_send:
stp x0,x1,[sp,#-256]!
stp x2,x3,[sp,#16]
stp x4,x5,[sp,#32]
stp x6,x7,[sp,#48]
stp x8,x9,[sp,#64]
stp x10,x11,[sp,#80]
stp x12,x13,[sp,#96]
stp x14,x15,[sp,#112]
stp x16,x17,[sp,#128]
stp x18,x19,[sp,#144]
stp x20,x21,[sp,#160]
stp x22,x23,[sp,#176]
stp x24,x25,[sp,#192]
stp x26,x27,[sp,#208]
stp x28,x29,[sp,#224]
str x30,[sp,#240]
mov x15,x0
CONST x1,16
ldp x2,x3,[x15]
lsr x4,x3,#20
cbnz x4,generic_error
str x2,[x1]
str x3,[x1,#8]
dsb sy
CONST x2,32
ldp x3,x4,[x1,#24]
stp x3,x4,[x2]
ldr x3,[x15,#16]
str x3,[x2,#16]
ldr x3,[x1,#40]
str x3,[x2,#24]
ldp x3,x4,[x15,#32]
stp x3,x4,[x2,#32]
ldp x3,x4,[x15,#48]
stp x3,x4,[x2,#48]
ldp x3,x4,[x15,#64]
stp x3,x4,[x2,#64]
str xzr,[x1,#16]
dsb sy
mov x3,#1
str x3,[x2,#256]
dsb sy
generic_send_wait:
ldr x3,[x1,#16]
cbnz x3,generic_send_done
wfe
b generic_send_wait
generic_send_done:
str xzr,[x1,#16]
dsb sy
ldp x0,x1,[sp]
ldp x2,x3,[sp,#16]
ldp x4,x5,[sp,#32]
ldp x6,x7,[sp,#48]
ldp x8,x9,[sp,#64]
ldp x10,x11,[sp,#80]
ldp x12,x13,[sp,#96]
ldp x14,x15,[sp,#112]
ldp x16,x17,[sp,#128]
ldp x18,x19,[sp,#144]
ldp x20,x21,[sp,#160]
ldp x22,x23,[sp,#176]
ldp x24,x25,[sp,#192]
ldp x26,x27,[sp,#208]
ldp x28,x29,[sp,#224]
ldr x30,[sp,#240]
add sp,sp,#256
ret

.globl generic_receive
generic_receive:
sub sp,sp,#32
str x30,[sp]
CONST x15,16
generic_receive_wait:
ldr x16,[x15,#48]
cbnz x16,generic_receive_ready
wfe
b generic_receive_wait
generic_receive_ready:
dsb sy
str xzr,[x15,#48]
dsb sy
CONST x0,48
// Private ACK identity survives application-visible publication and inbox reuse.
ldp x2,x3,[x0]
stp x2,x3,[sp,#8]
ldr x2,[x0,#24]
str x2,[sp,#24]
ldr x16,[x0,#16]
lsr x17,x16,#2
cbnz x17,generic_error
CONST x17,64
add x17,x17,x16,lsl #6
ldp x16,x9,[x17]
ldp x10,x11,[x17,#16]
ldp x28,x29,[x17,#32]
add x0,x0,#32
blr x16
dsb sy
ldp x2,x3,[sp,#8]
lsr x4,x3,#20
cbnz x4,generic_error
ldr x4,[sp,#24]
lsr x5,x4,#4
cbnz x5,generic_error
CONST x15,16
str x2,[x15]
str x3,[x15,#8]
dsb sy
CONST x0,40
add x0,x0,x4,lsl #3
mov x1,#1
str x1,[x0]
dsb sy
ldr x30,[sp]
add sp,sp,#32
ret

generic_error:
udf #0

// Receiver outer service loop: one B per completed message is charged.
.globl generic_receive_loop
generic_receive_loop:
bl generic_receive
b generic_receive_loop

// Once per endpoint: x0=CTRLbase,x1=self.low64,x2=self.high20.
.globl generic_control_init
generic_control_init:
mov x3,#0
mov x4,#16
1:
stp xzr,xzr,[x0]
str xzr,[x0,#16]
stp x1,x2,[x0,#24]
str x3,[x0,#40]
str xzr,[x0,#48]
add x0,x0,#64
add x3,x3,#1
subs x4,x4,#1
b.ne 1b
dsb sy
ret
// Once per endpoint: x0=HANDLERS,x1=dlv,x2=dlvq,x3=queue cursor address,
// x4=queue limit,x5=counter. Four descriptors, queue cursor and128-bit counter.
.globl generic_handlers_init
generic_handlers_init:
add x7,x3,#8
str x7,[x3]
str xzr,[x5]
add x7,x5,#2048
str xzr,[x7]
mov x6,#4
1:
cmp x6,#3
csel x7,x2,x1,eq
mov x8,#64
stp x7,x8,[x0]
stp x3,x4,[x0,#16]
mov x8,#1
stp x8,x5,[x0,#32]
stp xzr,xzr,[x0,#48]
add x0,x0,#64
subs x6,x6,#1
b.ne 1b
ret

// Once on each Phase0 collector G, after ordinary handler initialization.
// x0=HANDLERS,x1=byte entry dlvq,x2=M,x3=ticket start,x4=ticket limit.
// Specialize command2 only; command1 retains the production G request arena.
.globl generic_collector_init
generic_collector_init:
add x0,x0,#128
sub x5,x3,#8
str x3,[x5]
str xzr,[x2]
add x6,x2,#2048
str xzr,[x6]
mov x6,#64
stp x1,x6,[x0]
stp x5,x4,[x0,#16]
mov x6,#1
stp x6,x2,[x0,#32]
stp xzr,xzr,[x0,#48]
dsb sy
ret

// After a whole G request/reply wave has drained. Does not touch tickets.
// x0=request cursor word,x1=request arena start,x2=request completion counter.
.globl generic_request_reset
generic_request_reset:
str x1,[x0]
str xzr,[x2]
add x3,x2,#2048
str xzr,[x3]
dsb sy
ret

// Card54 finite Y archive/replay microprobe; source and receiver both counted.
// Local source descriptor at 34TiB+0x7200080. Source route (sigma48,bank10)
// writes feed a dedicated 58-bit address class, separate from record routing.
// desc: frame,sigma,bank,slot,valid|absolute_k,return_lo,return_hi,reserved.
// All data apertures and archive payload are Normal Inner/Outer NC. Control
// flags/routes are Device-nGnRnE, one active global sender, EVENTI pulses.
.text
.macro YSAVE
sub sp,sp,#272
stp x0,x1,[sp]
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
mrs x0,nzcv
stp x30,x0,[sp,#240]
.endm
.macro YRESTORE
ldp x30,x0,[sp,#240]
msr nzcv,x0
ldp x28,x29,[sp,#224]
ldp x26,x27,[sp,#208]
ldp x24,x25,[sp,#192]
ldp x22,x23,[sp,#176]
ldp x20,x21,[sp,#160]
ldp x18,x19,[sp,#144]
ldp x16,x17,[sp,#128]
ldp x14,x15,[sp,#112]
ldp x12,x13,[sp,#96]
ldp x10,x11,[sp,#80]
ldp x8,x9,[sp,#64]
ldp x6,x7,[sp,#48]
ldp x4,x5,[sp,#32]
ldp x2,x3,[sp,#16]
ldp x0,x1,[sp]
add sp,sp,#272
.endm
.macro COPY448 src,dst
ld64b x0,[\src]
st64b x0,[\dst]
add \src,\src,#64
add \dst,\dst,#64
ld64b x0,[\src]
st64b x0,[\dst]
add \src,\src,#64
add \dst,\dst,#64
ld64b x0,[\src]
st64b x0,[\dst]
add \src,\src,#64
add \dst,\dst,#64
ld64b x0,[\src]
st64b x0,[\dst]
add \src,\src,#64
add \dst,\dst,#64
ld64b x0,[\src]
st64b x0,[\dst]
add \src,\src,#64
add \dst,\dst,#64
ld64b x0,[\src]
st64b x0,[\dst]
add \src,\src,#64
add \dst,\dst,#64
ld64b x0,[\src]
st64b x0,[\dst]
add \src,\src,#64
add \dst,\dst,#64
.endm
.globl ystage
ystage:
YSAVE
ysbody:
mov x16,#0x7200000
movk x16,#0x2200,lsl #32
ysdescriptor:
add x18,x16,#128
ldp x10,x11,[x18]
ldp x12,x19,[x18,#16]
ldp x20,x21,[x18,#32]
ldp x22,x23,[x18,#48]
// Quiescent route change: previous invocation returned only after ACK+DSB.
stp x11,x12,[x16]
ldp x8,x9,[x16]
cmp x8,x11
b.ne yfatal
cmp x9,x12
b.ne yfatal
dsb sy
// 512-byte outgoing aperture at source CTRL+4096.
add x17,x16,#4096
stp x19,x20,[x17]
stp x21,x22,[x17,#16]
stp x23,xzr,[x17,#32]
stp xzr,xzr,[x17,#48]
add x10,x10,#64
add x17,x17,#64
COPY448 x10,x17
str xzr,[x16,#16]
dsb sy
mov x8,#1
str x8,[x16,#24]
dsb sy
yswait:
ldr x8,[x16,#16]
cbnz x8,ysdone
wfe
b yswait
ysdone:
// ACK 1=accepted/sink, 2=stale or invalid tag. Never silently continue error.
cmp x8,#1
b.ne yfatal
dsb sy
YRESTORE
ret
yfatal:
udf #1

// Hook immediately after WORDS updates [x29+2240] to its next frame.
// The CB has2304 bytes; its unused words2248..2288 hold k,c',u,group_first,
// returnlo,returnhi. Initialized once per chain, outside this per-Y path.
// yround keeps its existing256-chain forward/inversion/backward structure.
.globl yproduce
yproduce:
YSAVE
mov x16,#0x7200000
movk x16,#0x2200,lsl #32
add x18,x16,#128
ldr x10,[x29,#2240]
sub x10,x10,#512
ldr x11,[x29,#2248]
ldr x12,[x29,#2256]
ldr x13,[x29,#2264]
ldr x14,[x29,#2272]
ldr x15,[x10,#56]
lsl x15,x15,#15
and x8,x12,#0x7fff
orr x15,x15,x8
stp x10,x15,[x18]
sub x8,x11,x14
ubfx x9,x8,#19,#10
and x8,x8,#0x7ffff
stp x9,x8,[x18,#16]
orr x8,x11,#0x8000000000000000
ldr x9,[x29,#2280]
stp x8,x9,[x18,#32]
ldr x9,[x29,#2288]
stp x9,xzr,[x18,#48]
// Preserve absolute identity (group-relative r is used only for storage).
sub x11,x11,#1
lsr x11,x11,#36
lsr x12,x12,#15
orr x11,x11,x12,lsl #9
orr x11,x11,x13,lsl #42
str x11,[x10,#360]
add x13,x13,#1
str x13,[x29,#2264]
b ysbody

// Bank replay uses the same scalar transfer with a bank-local control page.
// Descriptor slot is0 at owner; source pointer advances in archive slot.
.globl yreplay
yreplay:
YSAVE
mov x16,#0x40000000
movk x16,#0x40,lsl #48
b ysdescriptor
.globl yreplay_loop
yreplay_loop:
stp x16,x30,[sp,#-16]!
mov x16,#0x40000000
movk x16,#0x40,lsl #48
ldr x8,[x16,#104]
cbz x8,yfatal
ldr x8,[x16,#192]
cbz x8,yreplay_end
yreplay_next:
bl yreplay
ldr x8,[x16,#128]
add x8,x8,#512
str x8,[x16,#128]
ldr x8,[x16,#192]
subs x8,x8,#1
str x8,[x16,#192]
b.ne yreplay_next
yreplay_end:
dsb sy
mov x8,#1
str x8,[x16,#56]
dsb sy
ldp x16,x30,[sp],#16
ret

// Metadata setup on each used bank: x0=first absolute round, x1=active slots.
// All counts, including empty prefixes, are reset. Payload is never zeroed.
.globl ybank_begin
ybank_begin:
mov x16,#0x40000000
movk x16,#0x40,lsl #48
mov x18,#0x40000000000000
orr x0,x0,#0x8000000000000000
cbz x1,ybank_initialized
ybank_initloop:
stp xzr,x0,[x18],#16
add x0,x0,#1
subs x1,x1,#1
b.ne ybank_initloop
ybank_initialized:
dsb sy
str xzr,[x16,#104]
dsb sy
ret
.globl ybank_seal
ybank_seal:
mov x16,#0x40000000
movk x16,#0x40,lsl #48
dsb sy
mov x8,#1
str x8,[x16,#104]
dsb sy
ret

// Bank request service after scalar request delivery (separate control row).
// x0=localround, x1=valid|absolute_k, x2=owner56; static bank return identity
// is provisioned at CTRL+224/+232. Each return ACK targets that fixed bank.
.globl yreplay_request
yreplay_request:
stp x16,x30,[sp,#-16]!
lsr x8,x0,#19
cbnz x8,yfatal
mov x16,#0x40000000
movk x16,#0x40,lsl #48
mov x18,#0x40000000000000
add x18,x18,x0,lsl #4
ldp x8,x9,[x18]
cmp x9,x1
b.ne yfatal
mov x9,#0x3ffffff
cmp x8,x9
csel x8,x8,x9,lo
str x8,[x16,#192]
lsl x0,x0,#35
stp x0,x2,[x16,#128]
stp xzr,xzr,[x16,#144]
ldp x8,x9,[x16,#224]
stp x1,x8,[x16,#160]
stp x9,xzr,[x16,#176]
bl yreplay_loop
ldp x16,x30,[sp],#16
ret

// One dispatched receive. No caller state to preserve: dedicated receiver PE.
// CTRL=2^54+2^30, inbox=2^54+2^25, metadata=2^54.
// CTRL+96 has data base (bank0; owner bufferA). CTRL+104 is sealed state.
.globl yreceive
yreceive:
mov x16,#0x40000000
movk x16,#0x40,lsl #48
mov x17,#0x2000000
movk x17,#0x40,lsl #48
mov x18,#0x40000000000000
b yrwait
// Owner uses its existing G map, with disjoint scratch inside G34.
.globl yreceive_owner
yreceive_owner:
mov x16,#0x7500000
movk x16,#0x2200,lsl #32
mov x17,#0x7510000
movk x17,#0x2200,lsl #32
mov x18,#0x7520000
movk x18,#0x2200,lsl #32
yrwait:
ldr x8,[x16,#24]
cbnz x8,yrready
wfe
b yrwait
yrready:
str xzr,[x16,#24]
dsb sy
ldp x19,x20,[x17]
ldp x21,x22,[x17,#16]
ldp x23,x24,[x17,#32]
ldp x25,x26,[x17,#48]
lsr x8,x19,#19
cbnz x8,yrreject
add x18,x18,x19,lsl #4
ldp x27,x28,[x18]
cmp x20,x28
b.ne yrreject
ldr x8,[x16,#104]
cbnz x8,yrreject
mov x29,#0x3ffffff
cmp x27,x29
csel x27,x27,x29,lo
ldr x15,[x16,#96]
add x15,x15,x19,lsl #35
add x15,x15,x27,lsl #9
add x15,x15,#64
add x17,x17,#64
COPY448 x17,x15
add x27,x27,#1
cmp x27,x29
csel x27,x27,x29,lo
str x27,[x18]
mov x23,#1
b yrack
yrreject:
mov x23,#2
yrack:
// Drain all receiver payload writes before ACK and before sender may reuse.
dsb sy
stp x21,x22,[x16,#32]
ldp x8,x9,[x16,#32]
cmp x8,x21
b.ne yfatal
cmp x9,x22
b.ne yfatal
dsb sy
str x23,[x16,#48]
dsb sy
ret

// Owner's static buffer A is exposed read-only at VIEW=38TiB in each member.
// A group-local fixed1024:1 path uses only member10; no global bank route.
// count,epoch are in local descriptor CTRL+128; first member/count at+144.
// Owner main uses a separate local completion from receiver's replay DONE.
.text
.globl yoffer_scalar
yoffer_scalar:
stp x19,x30,[sp,#-16]!
stp x20,x21,[sp,#-16]!
stp x22,x23,[sp,#-16]!
mov x19,#0x7300000
movk x19,#0x2200,lsl #32
ldp x20,x21,[x19,#128]
ldp x22,x23,[x19,#144]
yo_loop:
str x22,[x19]
ldr x8,[x19]
cmp x8,x22
b.ne yo_fatal
dsb sy
add x9,x19,#4096
stp x20,x21,[x9]
str xzr,[x19,#16]
dsb sy
mov x8,#1
str x8,[x19,#24]
dsb sy
yo_wait:
ldr x8,[x19,#16]
cbnz x8,yo_returned
wfe
b yo_wait
yo_returned:
cmp x8,#1
b.ne yo_fatal
dsb sy
add x22,x22,#1
subs x23,x23,#1
b.ne yo_loop
ldp x22,x23,[sp],#16
ldp x20,x21,[sp],#16
ldp x19,x30,[sp],#16
ret
yo_fatal:
udf #2

.globl yaccept_scalar
yaccept_scalar:
stp x19,x30,[sp,#-16]!
mov x19,#0x7300000
movk x19,#0x2200,lsl #32
ya_wait:
ldr x8,[x19,#24]
cbnz x8,ya_ready
wfe
b ya_wait
ya_ready:
str xzr,[x19,#24]
dsb sy
add x9,x19,#4096
ldp x1,x3,[x9]
// Published epoch is the absolute k, not local round. Boot fields are m,c'low15.
ldr x8,[x19,#160]
ldr x9,[x19,#168]
sub x10,x3,#1
and x10,x10,#0xfffffffff
lsl x10,x10,#28
orr x10,x10,x8,lsl #3
orr x10,x10,x9,lsl #13
mov x11,#0x7100000
movk x11,#0x2200,lsl #32
str x10,[x11,#48]
mov x0,#0x260000000000
mov x2,#0x7120000
movk x2,#0x2200,lsl #32
bl yconsume
// yconsume may clobber every GPR; restore control base from stack.
ldr x19,[sp]
// Saved x19 is caller's value, not control. Materialize explicitly.
mov x19,#0x7300000
movk x19,#0x2200,lsl #32
dsb sy
mov x8,#1
str x8,[x19,#48]
dsb sy
ldp x19,x30,[sp],#16
ret

// Scalar control commands. G source CMDCTRL=34TiB+0x7600000. Descriptor
//128B is ordinary memory. Route is the root84-bit(role,endpoint) split64/20.
// Header: opcode,arg0,arg1,arg2,return_lo,return_hi,0,0; payload64B.
// Opcode1 begin,2 seal,3 release,4 replay,5 CN coins,6 main token.
.globl ycmd_send
ycmd_send:
stp x19,x30,[sp,#-16]!
mov x19,#0x7600000
movk x19,#0x2200,lsl #32
stp x1,x2,[x19]
ldp x8,x9,[x19]
cmp x8,x1
b.ne yfatal
cmp x9,x2
b.ne yfatal
dsb sy
mov x10,x0
add x11,x19,#4096
ld64b x0,[x10]
st64b x0,[x11]
add x10,x10,#64
add x11,x11,#64
ld64b x0,[x10]
st64b x0,[x11]
str xzr,[x19,#16]
dsb sy
mov x8,#1
str x8,[x19,#24]
dsb sy
yc_wait:
ldr x8,[x19,#16]
cbnz x8,yc_sent
wfe
b yc_wait
yc_sent:
cmp x8,#1
b.ne yfatal
dsb sy
ldp x19,x30,[sp],#16
ret

// Command entries are separate from payload entries; clears precede replay.
.globl ycmd_bank
ycmd_bank:
mov x16,#0x40000000
add x16,x16,#512
movk x16,#0x40,lsl #48
mov x17,#0x2100000
movk x17,#0x40,lsl #48
mov x15,#3
b yc_receive
.globl ycmd_gen
ycmd_gen:
mov x16,#0x7600000
add x16,x16,#512
movk x16,#0x2200,lsl #32
mov x17,#0x7610000
movk x17,#0x2200,lsl #32
mov x15,#0
yc_receive:
stp x19,x30,[sp,#-16]!
sub sp,sp,#64
stp x16,x17,[sp]
str x15,[sp,#48]
yc_recvwait:
ldr x8,[x16,#24]
cbnz x8,yc_ready
wfe
b yc_recvwait
yc_ready:
str xzr,[x16,#24]
dsb sy
ldp x19,x0,[x17]
ldp x1,x2,[x17,#16]
ldp x8,x9,[x17,#32]
stp x8,x9,[sp,#16]
stp x0,x2,[sp,#32]
cmp x19,#5
b.hs yc_generator
cmp x15,#3
b.ne yfatal
sub x18,x16,#512
cmp x19,#1
b.eq yc_begin
ldr x8,[x18,#104]
cmp x19,#4
b.eq yc_replay
ldr x9,[x18,#112]
cmp x9,x0
b.ne yfatal
cmp x19,#2
b.eq yc_seal
cmp x19,#3
b.ne yfatal
cmp x8,#1
b.ne yfatal
mov x8,#2
str x8,[x18,#104]
b yc_ack
yc_begin:
mov x8,#0x80000
cmp x8,x1
b.lo yfatal
ldr x8,[x18,#104]
cmp x8,#2
b.ne yfatal
str x2,[x18,#112]
bl ybank_begin
b yc_ack
yc_seal:
cbnz x8,yfatal
bl ybank_seal
b yc_ack
yc_replay:
cmp x8,#1
b.ne yfatal
bl yreplay_request
b yc_ack
yc_generator:
cbnz x15,yfatal
cmp x19,#6
b.eq yc_token
cmp x19,#5
b.ne yfatal
// localc' boot field at command-receiverCTRL+104.
ldr x8,[x16,#104]
cmp x8,x1
b.ne yfatal
cbz x0,yfatal
sub x0,x0,#1
lsr x8,x0,#45
cbnz x8,yfatal
mov x10,#0x540000000000
add x10,x10,x0,lsl #6
add x11,x17,#64
ld64b x0,[x11]
st64b x0,[x10]
b yc_ack
yc_token:
// arg0=token kind, arg1=epoch; scalar publication pulses main's EVENTI.
cbz x0,yfatal
cmp x0,#6
b.hs yfatal
stp x1,x2,[x16,#88]
ldp x8,x9,[sp,#16]
stp x8,x9,[x16,#112]
dsb sy
str x0,[x16,#80]
yc_ack:
dsb sy
ldr x16,[sp]
ldp x8,x9,[sp,#16]
stp x8,x9,[x16,#32]
ldp x10,x11,[x16,#32]
cmp x8,x10
b.ne yfatal
cmp x9,x11
b.ne yfatal
dsb sy
mov x8,#1
str x8,[x16,#48]
dsb sy
add sp,sp,#64
ldp x19,x30,[sp],#16
ret

// Perpetual dispatchers: 2 added instructions/service, including CALL+loop.
// A selected processor has exactly one service kind; address-class routing
// selects that processor's READY register. There is no polling of other kinds.
.globl ybank_payload_dispatch
ybank_payload_dispatch:
bl yreceive
b ybank_payload_dispatch
.globl yowner_payload_dispatch
yowner_payload_dispatch:
bl yreceive_owner
b yowner_payload_dispatch
.globl ybank_command_dispatch
ybank_command_dispatch:
bl ycmd_bank
b ybank_command_dispatch
.globl ygen_command_dispatch
ygen_command_dispatch:
bl ycmd_gen
b ygen_command_dispatch
.globl ymember_offer_dispatch
ymember_offer_dispatch:
bl yaccept_scalar
b ymember_offer_dispatch

// One owner round: x0=absolute k,x1=group_first_k,x2=owner sigma.
// The command ACK occurs after all replay payload ACKs and bank return.
// Owner receive PE is distinct and never waits for this main PE.
.globl yowner_round
yowner_round:
sub sp,sp,#64
stp x19,x30,[sp]
stp x0,x1,[sp,#16]
str x2,[sp,#32]
mov x19,#0x7520000
movk x19,#0x2200,lsl #32
orr x8,x0,#0x8000000000000000
stp xzr,x8,[x19]
mov x19,#0x7500000
movk x19,#0x2200,lsl #32
mov x9,#0x280000000000
str x9,[x19,#96]
str xzr,[x19,#104]
str xzr,[x19,#56]
dsb sy
mov x19,#0x7620000
movk x19,#0x2200,lsl #32
sub x9,x0,x1
and x10,x9,#0x7ffff
mov x11,#4
stp x11,x10,[x19]
stp x8,x2,[x19,#16]
lsl x10,x2,#10
lsr x11,x2,#54
stp x10,x11,[x19,#32]
stp xzr,xzr,[x19,#48]
lsr x9,x9,#19
orr x1,x10,x9
orr x2,x11,#0xc0000
mov x0,x19
bl ycmd_send
// Every archive frame is now in owner bufferA. Set absolute record w2.
ldr x3,[sp,#16]
ldr x9,[sp,#32]
and x9,x9,#0x7fff
sub x10,x3,#1
and x10,x10,#0xfffffffff
lsl x10,x10,#28
orr x10,x10,x9,lsl #13
mov x11,#0x7100000
movk x11,#0x2200,lsl #32
str x10,[x11,#48]
mov x19,#0x7520000
movk x19,#0x2200,lsl #32
ldr x1,[x19]
mov x0,#0x280000000000
mov x2,#0x7120000
movk x2,#0x2200,lsl #32
bl yconsume
// yconsume may clobber allGPRs; all persistent values are stack/metadata.
mov x19,#0x7520000
movk x19,#0x2200,lsl #32
ldr x8,[x19]
ldr x9,[sp,#16]
mov x19,#0x7300000
movk x19,#0x2200,lsl #32
stp x8,x9,[x19,#128]
mov x8,#1
mov x9,#1023
stp x8,x9,[x19,#144]
str xzr,[x19,#176]
bl yoffer_scalar
mov x19,#0x7300000
movk x19,#0x2200,lsl #32
dsb sy
mov x8,#1
str x8,[x19,#176]
dsb sy
ldp x19,x30,[sp]
add sp,sp,#64
ret

// CHAIN_SET consumes exactly the64B CN pair at84TiB+64(k-1).
// x0=CB pointer,x1=absolute k,x2=c',x3=groupfirst,x4=returnlo,x5=returnhi,
// x6=YR framebase. No coin draw is repeated in this metadata initializer.
.globl ychain_metadata
ychain_metadata:
str x6,[x0,#2240]
str x1,[x0,#2248]
str x2,[x0,#2256]
str xzr,[x0,#2264]
str x3,[x0,#2272]
str x4,[x0,#2280]
str x5,[x0,#2288]
sub x8,x1,#1
mov x9,#0x540000000000
add x9,x9,x8,lsl #6
str x9,[x0,#2296]
ret

// First assigned k in inclusive[first,last] for physical member m; no division
// oracle. x0=first,x1=last,x2=m. Return x0=first_k,x1=number of chains.
.globl yassigned
yassigned:
sub x8,x0,#1
sub x8,x2,x8
and x8,x8,#1023
add x0,x0,x8
cmp x1,x0
b.lo yassigned_empty
sub x1,x1,x0
lsr x1,x1,#10
add x1,x1,#1
ret
yassigned_empty:
mov x1,#0
ret

// Increment a production endpoint66-bit counter; a true64bit carry is counted.
.globl ynext_endpoint
ynext_endpoint:
adds x0,x0,#1
adc x1,x1,xzr
ret

// Publish one64B(c,s) pair to all1024 original members of c'; self delivery
// also runs the scalar command processor and is charged. x0=k,x1=c',x2=pair.
.globl ycoins_publish
ycoins_publish:
stp x19,x30,[sp,#-16]!
stp x20,x21,[sp,#-16]!
stp x22,x23,[sp,#-16]!
mov x20,x0
mov x21,x1
mov x22,x2
mov x19,#0x7620000
movk x19,#0x2200,lsl #32
mov x8,#5
stp x8,x20,[x19]
stp x21,xzr,[x19,#16]
mov x8,#0x7600000
movk x8,#0x2200,lsl #32
ldp x9,x10,[x8,#64]
stp x9,x10,[x19,#32]
stp xzr,xzr,[x19,#48]
ld64b x0,[x22]
add x8,x19,#64
st64b x0,[x8]
mov x23,#0
ycoins_next:
lsl x1,x21,#10
orr x1,x1,x23
lsr x2,x21,#54
mov x0,x19
bl ycmd_send
add x23,x23,#1
cmp x23,#1024
b.ne ycoins_next
ldp x22,x23,[sp],#16
ldp x20,x21,[sp],#16
ldp x19,x30,[sp],#16
ret

// Two bounded attempts, full256bit subgroup range; x1=0 permits zero for c,
// x1=1 excludes zero for s. Input x0 points to four scalar words; result0/1.
.globl ycoin_one
ycoin_one:
mov x15,#0x7100000
movk x15,#0x2200,lsl #32
ldp x10,x11,[x15,#96]
ldp x12,x13,[x15,#112]
mov x8,#2
ycoin_try:
mrs x4,rndr
mrs x5,rndr
mrs x6,rndr
mrs x7,rndr
cbz x1,ycoin_range
orr x9,x4,x5
orr x14,x6,x7
orr x9,x9,x14
cbz x9,ycoin_retry
ycoin_range:
subs xzr,x4,x10
sbcs xzr,x5,x11
sbcs xzr,x6,x12
sbcs xzr,x7,x13
b.lo ycoin_accept
ycoin_retry:
subs x8,x8,#1
b.ne ycoin_try
mov x0,#1
ret
ycoin_accept:
stp x4,x5,[x0]
stp x6,x7,[x0,#16]
mov x0,#0
ret
.globl ydraw_pair
ydraw_pair:
stp x19,x30,[sp,#-16]!
mov x19,x0
mov x1,#0
bl ycoin_one
cbnz x0,ydraw_done
add x0,x19,#32
mov x1,#1
bl ycoin_one
ydraw_done:
ldp x19,x30,[sp],#16
ret

// Production driver state atG34+0x7640000:
// c',m,first,last,returnlo,returnhi,next_k,remaining,n,b,Y_iter,
// CHAIN_SET pointer,yround pointer,YSP; callerSP/LR at128/136.
// CHAIN_SET(x0=CB,x1=CN64) is the unchanged, separately priced2GEN_MUL+
//17CHAIN_ADD primitive. yround is the existing literal entry with WORDS hook.
.globl yproduce_group
yproduce_group:
mov x19,#0x7640000
movk x19,#0x2200,lsl #32
mov x8,sp
stp x8,x30,[x19,#128]
ldp x0,x1,[x19,#16]
ldr x2,[x19,#8]
bl yassigned
stp x0,x1,[x19,#48]
cbz x1,yg_success
yg_batch:
mov x8,#256
ldr x9,[x19,#56]
cmp x9,x8
csel x8,x9,x8,lo
str x8,[x19,#64]
str xzr,[x19,#72]
yg_prepare:
mov x0,#0x7630000
movk x0,#0x2200,lsl #32
bl ydraw_pair
cbnz x0,yg_failure
mov x19,#0x7640000
movk x19,#0x2200,lsl #32
ldr x8,[x19,#72]
ldr x0,[x19,#48]
add x0,x0,x8,lsl #10
ldr x1,[x19]
mov x2,#0x7630000
movk x2,#0x2200,lsl #32
bl ycoins_publish
mov x19,#0x7640000
movk x19,#0x2200,lsl #32
ldr x8,[x19,#72]
add x9,x8,x8,lsl #3
mov x0,#0x2000000
movk x0,#0x2200,lsl #32
add x0,x0,x9,lsl #8
ldr x9,[x19,#48]
add x9,x9,x8,lsl #10
sub x9,x9,#1
mov x1,#0x540000000000
add x1,x1,x9,lsl #6
ldr x16,[x19,#88]
blr x16
// All primitive registers may be clobbered; recompute from persistent state.
mov x19,#0x7640000
movk x19,#0x2200,lsl #32
ldr x8,[x19,#72]
add x9,x8,x8,lsl #3
mov x0,#0x2000000
movk x0,#0x2200,lsl #32
add x0,x0,x9,lsl #8
ldr x1,[x19,#48]
add x1,x1,x8,lsl #10
ldr x2,[x19]
ldr x3,[x19,#16]
ldp x4,x5,[x19,#32]
add x8,x8,#1
mov x6,#0x200000000000
add x6,x6,x8,lsl #31
bl ychain_metadata
ldr x8,[x19,#72]
add x8,x8,#1
str x8,[x19,#72]
ldr x9,[x19,#64]
cmp x8,x9
b.lo yg_prepare
ldr x10,[x19,#104]
mov x8,#0x2000000
movk x8,#0x2200,lsl #32
str x8,[x10,#424]
sub x11,x8,#128
mov x12,#1
stp x12,xzr,[x11]
stp xzr,xzr,[x11,#16]
add x9,x9,x9,lsl #3
add x8,x8,x9,lsl #8
str x8,[x10,#416]
mov x8,#1
str x8,[x10,#1016]
mov x8,#63520
str x8,[x19,#80]
yg_round:
ldr x10,[x19,#104]
mov sp,x10
ldr x16,[x19,#96]
blr x16
mov x19,#0x7640000
movk x19,#0x2200,lsl #32
ldr x8,[x19,#128]
mov sp,x8
ldr x8,[x19,#80]
subs x8,x8,#1
str x8,[x19,#80]
b.ne yg_round
ldr x8,[x19,#64]
ldp x0,x1,[x19,#48]
add x0,x0,x8,lsl #10
sub x1,x1,x8
stp x0,x1,[x19,#48]
cbnz x1,yg_batch
yg_success:
mov x0,#0
b yg_return
yg_failure:
mov x0,#1
mov x19,#0x7640000
movk x19,#0x2200,lsl #32
yg_return:
ldp x8,x30,[x19,#128]
mov sp,x8
ret

// Blocking work-token call, separate command ACK and DONE token.
// x0/1=destination84bits low/high,x2=kind1generate2owner4fin,
// x3=first absolute k,x4=last. DONE kind3 echoes first andstatus0/1.
.globl ytoken_call
ytoken_call:
sub sp,sp,#64
stp x19,x30,[sp]
stp x0,x1,[sp,#16]
str x3,[sp,#32]
mov x19,#0x7620000
movk x19,#0x2200,lsl #32
mov x8,#6
stp x8,x2,[x19]
stp x3,x4,[x19,#16]
mov x8,#0x7600000
movk x8,#0x2200,lsl #32
ldp x9,x10,[x8,#64]
stp x9,x10,[x19,#32]
stp xzr,xzr,[x19,#48]
ldp x1,x2,[sp,#16]
mov x0,x19
bl ycmd_send
mov x19,#0x7600000
movk x19,#0x2200,lsl #32
add x19,x19,#512
yt_donewait:
ldr x8,[x19,#80]
cbnz x8,yt_done
wfe
b yt_donewait
yt_done:
cmp x8,#3
b.ne yfatal
ldr x8,[x19,#88]
ldr x9,[sp,#32]
cmp x8,x9
b.ne yfatal
ldr x0,[x19,#96]
str xzr,[x19,#80]
dsb sy
ldp x19,x30,[sp]
add sp,sp,#64
ret

// One already dispatched main-token service. Main returns only after its
// work and its explicit DONE command ACK. CMD delivery processor is distinct.
.globl yworker_once
yworker_once:
sub sp,sp,#64
stp x19,x30,[sp]
mov x19,#0x7600000
movk x19,#0x2200,lsl #32
add x19,x19,#512
yw_wait:
ldr x8,[x19,#80]
cbnz x8,yw_ready
wfe
b yw_wait
yw_ready:
ldp x9,x10,[x19,#88]
ldp x11,x12,[x19,#112]
stp x11,x12,[sp,#16]
stp x9,x10,[sp,#32]
str xzr,[x19,#80]
dsb sy
mov x19,#0x7640000
movk x19,#0x2200,lsl #32
stp x9,x10,[x19,#16]
cmp x8,#1
b.eq yw_generate
cmp x8,#2
b.eq yw_owner
cmp x8,#4
b.eq yw_fin
cmp x8,#5
b.ne yfatal
// Phase0 callback is the separately priced original mining primitive.
// It reads this G's fixed identity from STATE and returns status0/1 in x0.
ldr x16,[x19,#176]
blr x16
cbnz x0,yw_done
mov x19,#0x7640000
movk x19,#0x2200,lsl #32
yw_fin:
// Actual fin1 inputs: LC1 address0, F1=2^40 destinations.
mov x7,#0
mov x1,#0x10000000000
ldr x16,[x19,#112]
blr x16
dsb sy
mov x0,#0
b yw_done
yw_generate:
bl yproduce_group
b yw_done
yw_owner:
ldr x8,[x19,#8]
cbnz x8,yfatal
bl yowner_group
yw_done:
// x0 is run failure status. Do not redraw after a two-attempt coin failure.
mov x19,#0x7620000
movk x19,#0x2200,lsl #32
mov x8,#6
mov x9,#3
stp x8,x9,[x19]
ldr x8,[sp,#32]
stp x8,x0,[x19,#16]
mov x8,#0x7600000
movk x8,#0x2200,lsl #32
ldp x9,x10,[x8,#64]
stp x9,x10,[x19,#32]
stp xzr,xzr,[x19,#48]
ldp x1,x2,[sp,#16]
mov x0,x19
bl ycmd_send
ldp x19,x30,[sp]
add sp,sp,#64
ret
.globl yworker_dispatch
yworker_dispatch:
bl yworker_once
b yworker_dispatch

// Owner consumes every k before its1024 G finals. No H partial-final orRCB
// reset occurs here. All1023 remote fin calls use main worktokens andDONE.
.globl yowner_group
yowner_group:
sub sp,sp,#80
stp x19,x30,[sp]
mov x19,#0x7640000
movk x19,#0x2200,lsl #32
ldr x8,[x19]
ldp x9,x10,[x19,#16]
stp x8,x9,[sp,#16]
stp x9,x10,[sp,#32]
yog_round:
ldr x0,[sp,#32]
ldr x1,[sp,#24]
ldr x2,[sp,#16]
bl yowner_round
ldp x8,x9,[sp,#32]
add x8,x8,#1
str x8,[sp,#32]
cmp x9,x8
b.hs yog_round
// The owner itself cannot send a blocking worktoken to its own main.
mov x19,#0x7640000
movk x19,#0x2200,lsl #32
mov x7,#0
mov x1,#0x10000000000
ldr x16,[x19,#112]
blr x16
mov x8,#1
str x8,[sp,#48]
yog_fin:
ldr x8,[sp,#16]
ldr x9,[sp,#48]
lsl x0,x8,#10
orr x0,x0,x9
lsr x1,x8,#54
mov x2,#4
ldr x3,[sp,#24]
ldr x4,[sp,#40]
bl ytoken_call
cbnz x0,yog_return
ldr x8,[sp,#48]
add x8,x8,#1
str x8,[sp,#48]
cmp x8,#1024
b.ne yog_fin
mov x0,#0
yog_return:
ldp x19,x30,[sp]
add sp,sp,#80
ret

// All 2^58 production endpoints, stride 1 generation or 1024 owners.
// x0=kind,x1=first,x2=last,x3=stride. No H-final/Phase4 inside this sweep.
.globl yworker_sweep
yworker_sweep:
sub sp,sp,#80
stp x19,x30,[sp]
stp x0,x1,[sp,#16]
stp x2,x3,[sp,#32]
stp xzr,xzr,[sp,#48]
yws_next:
ldr x8,[sp,#16]
cmp x8,#2
b.eq yws_coarse
cmp x8,#5
b.ne yws_routed
// Mining uses each slot exactly once, in increasing Gid. R/H partials
// persist. Change coarse base only before the first member of each group.
ldr x9,[sp,#48]
and x9,x9,#1023
cbnz x9,yws_routed
yws_coarse:
// Root's explicitly priced 8*N_REL subgroup-base routine; not free broadcast.
ldp x0,x1,[sp,#48]
lsr x0,x0,#10
orr x0,x0,x1,lsl #54
mov x19,#0x7640000
movk x19,#0x2200,lsl #32
ldr x16,[x19,#168]
blr x16
yws_routed:
ldp x0,x1,[sp,#48]
ldp x2,x3,[sp,#16]
ldr x4,[sp,#32]
bl ytoken_call
cbnz x0,yws_return
ldp x8,x9,[sp,#48]
ldr x10,[sp,#40]
adds x8,x8,x10
adc x9,x9,xzr
stp x8,x9,[sp,#48]
cbnz x9,yfatal
mov x11,#0x400000000000000
cmp x8,x11
b.lo yws_next
cmp x8,x11
b.ne yfatal
mov x0,#0
yws_return:
ldp x19,x30,[sp]
add sp,sp,#80
ret

// Active banks only: all 2^48 owners, b=0..ceil(q/2^19)-1.
// x0=command 1 begin 2 seal 3 release, x1=first, x2=last.
.globl ybank_sweep
ybank_sweep:
sub sp,sp,#96
stp x19,x30,[sp]
stp x0,x1,[sp,#16]
str x2,[sp,#32]
sub x8,x2,x1
add x8,x8,#1
mov x9,#0x20000000
cmp x9,x8
b.lo yfatal
str x8,[sp,#40]
str xzr,[sp,#48]
ybs_owner:
str xzr,[sp,#56]
ybs_bank:
mov x19,#0x7620000
movk x19,#0x2200,lsl #32
ldp x8,x9,[sp,#16]
cmp x8,#1
b.ne ybs_notbegin
ldr x10,[sp,#56]
lsl x10,x10,#19
ldr x11,[sp,#40]
sub x11,x11,x10
mov x12,#0x80000
cmp x11,x12
csel x11,x11,x12,lo
add x10,x9,x10
stp x8,x10,[x19]
stp x11,x9,[x19,#16]
b ybs_header
ybs_notbegin:
stp x8,x9,[x19]
stp xzr,xzr,[x19,#16]
ybs_header:
// Dedicated controller G-id 2^58: (low64,high20)=(2^58,0).
mov x8,#0x400000000000000
stp x8,xzr,[x19,#32]
stp xzr,xzr,[x19,#48]
ldr x8,[sp,#48]
ldr x9,[sp,#56]
lsl x1,x8,#10
orr x1,x1,x9
lsr x2,x8,#54
orr x2,x2,#0xc0000
mov x0,x19
bl ycmd_send
ldr x8,[sp,#56]
add x8,x8,#1
str x8,[sp,#56]
lsl x9,x8,#19
ldr x10,[sp,#40]
cmp x9,x10
b.lo ybs_bank
ldr x8,[sp,#48]
add x8,x8,#1
str x8,[sp,#48]
mov x9,#0x1000000000000
cmp x8,x9
b.lo ybs_owner
ldp x19,x30,[sp]
add sp,sp,#96
ret

// One entered group. Root invokes this only after the preceding whole-group
// stop decision continues. Return failure 1 exits the opened run with no output.
// Whole-group Hfinal/Phase4/stop is an unchanged root call AFTER successful RET;
// no future group's coins or archive are generated by this function.
.globl ygroup_run
ygroup_run:
sub sp,sp,#48
stp x19,x30,[sp]
stp x0,x1,[sp,#16]
mov x2,x1
mov x1,x0
mov x0,#1
bl ybank_sweep
mov x0,#1
ldp x1,x2,[sp,#16]
mov x3,#1
bl yworker_sweep
cbnz x0,ygr_return
mov x0,#2
ldp x1,x2,[sp,#16]
bl ybank_sweep
mov x0,#2
ldp x1,x2,[sp,#16]
mov x3,#1024
bl yworker_sweep
cbnz x0,ygr_return
mov x0,#3
ldp x1,x2,[sp,#16]
bl ybank_sweep
mov x0,#0
ygr_return:
ldp x19,x30,[sp]
add sp,sp,#48
ret

// One-time local scalar initialization; code/map delivery is separately
// charged in cold reserve. x0/1=G endpoint low/high (controller 2^58,0 allowed),
// x2=CHAIN_SET,x3=hooked yround,x4=fin1,x5=original persistent YSP.
.globl yboot_g
yboot_g:
cbnz x1,yfatal
mov x8,#0x400000000000000
cmp x8,x0
b.lo yfatal
ybg_valid:
lsr x8,x0,#10
orr x8,x8,x1,lsl #54
and x9,x0,#1023
mov x10,#0x7640000
movk x10,#0x2200,lsl #32
stp x8,x9,[x10]
stp x0,x1,[x10,#32]
stp x2,x3,[x10,#88]
stp x5,x4,[x10,#104]
mov x10,#0x7600000
movk x10,#0x2200,lsl #32
stp x0,x1,[x10,#64]
stp xzr,xzr,[x10,#16]
str xzr,[x10,#536]
str xzr,[x10,#592]
str x8,[x10,#616]
mov x10,#0x7300000
movk x10,#0x2200,lsl #32
and x8,x8,#0x7fff
stp x9,x8,[x10,#160]
stp xzr,xzr,[x10,#16]
str xzr,[x10,#48]
mov x10,#0x7500000
movk x10,#0x2200,lsl #32
str xzr,[x10,#24]
mov x8,#0x280000000000
str x8,[x10,#96]
str xzr,[x10,#104]
mov x10,#0x7520000
movk x10,#0x2200,lsl #32
stp xzr,xzr,[x10]
mov x10,#0x7200000
movk x10,#0x2200,lsl #32
stp xzr,xzr,[x10,#16]
mov x10,#0x7620000
movk x10,#0x2200,lsl #32
stp xzr,xzr,[x10]
stp xzr,xzr,[x10,#16]
stp xzr,xzr,[x10,#32]
stp xzr,xzr,[x10,#48]
stp xzr,xzr,[x10,#64]
stp xzr,xzr,[x10,#80]
stp xzr,xzr,[x10,#96]
stp xzr,xzr,[x10,#112]
dsb sy
ret
.globl yboot_bank
yboot_bank:
cbnz x1,yfatal
lsr x8,x0,#58
cbnz x8,yfatal
mov x16,#0x40000000
movk x16,#0x40,lsl #48
stp x0,x1,[x16,#224]
stp xzr,xzr,[x16,#16]
str xzr,[x16,#56]
str xzr,[x16,#96]
mov x8,#2
str x8,[x16,#104]
str xzr,[x16,#112]
str xzr,[x16,#536]
dsb sy
ret
.globl yboot_controller
yboot_controller:
mov x8,#0x7640000
movk x8,#0x2200,lsl #32
str x0,[x8,#168]
ret

// Cold Phase0 entry. yworker_once kind5 runs the fixed G mining primitive,
// then fin1+DSB before sending DONE; the next G is not started before DONE.
// All R/H partials remain live until the root's whole-Phase0 finalization.
.globl yboot_mining
yboot_mining:
mov x8,#0x7640000
movk x8,#0x2200,lsl #32
str x0,[x8,#176]
ret
.globl yphase0_mine
yphase0_mine:
mov x0,#5
mov x1,#0
mov x2,#0
mov x3,#1
b yworker_sweep

// Local H reset AFTER whole-Phase0 final/drain and original HHclear complete.
// Caller holds H receiver quiescent until this returns; no concurrent payload.
// x0=physical R2 start0; x1=physical sink/limit=32*CAP_H;
// x2=H scanner frame=LHR+64; x3=initial logical R2 cursor=-5460TiB.
// Existing frame[x2-8] names its E2 word at LHR+128.
// This reinitializes both coordinate systems; it never clears R/H mid-phase.
.globl yhome_reset
yhome_reset:
dsb sy
mov x5,#0x7ffffffc000000
str xzr,[x5]
stp x0,x1,[x5,#16]
str x3,[x2,#-24]
str xzr,[x2,#-16]
ldr x4,[x2,#-8]
str x3,[x4]
dsb sy
ret

// One-time R record-return selector reset; no EVENTI pulse. The fixed ingress
// source-PE shim writes this latch on RREADY commit; READY clear preserves it.
// Only matching completed G ACK clears valid during normal operation.
.globl rreturn_pe_boot
rreturn_pe_boot:
mov x0,#0x8020
str xzr,[x0]
dsb sy
ret

// Dedicated controller's fixed8TiB..16TiB aperture addressesRcoarse0x8018.
// ArithmeticG's same numeric VA class addressesRREADY instead. This is static
// endpoint-role wiring, not a runtime binding operation. Coarse writes do not
// pulse EVENTI. Root's8*N_REL+64/owner allowance includes this exact routine.
.globl ycontroller_coarse
ycontroller_coarse:
mov x1,#0x80000000000
mov x2,#0x10000000000
ycoarse_loop:
str x0,[x1]
dsb sy
ldr x3,[x1]
cmp x3,x0
b.ne ycoarse_error
add x1,x1,#8
subs x2,x2,#1
b.ne ycoarse_loop
ret
ycoarse_error:
udf #0

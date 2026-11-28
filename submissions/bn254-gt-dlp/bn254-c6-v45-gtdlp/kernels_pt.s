// kernels_pt.s -- variant B of the lattice band-2 tree: E3 hand-off, patch rows, macro-generated chains, builders.
// Geometry is parametrised by the .set constants below (one-line edits): CW cell width exponent of the level-1/0 centre fractions (24 -> 256 cells per axis,
// 23 -> 512), PB patch energy-bin shift, CW3/EB the E3 cell/bin shifts, NMAXPOS chord piece length, KKEY/OFFKEY key window.  Entry points: pt_handoff3, pt_entry_m2/m1/m0,
// pt_orbit_build, pt_frame_init, pt_global_init.  External symbol: pt_l3_next (= .Ltop_3 of the old tree; defined here only when assembled with -defsym PT_TEST=1).
.arch armv8.1-a
.text
.set CW,24
.set PB,28
.set CW3,22
.set EB,26
.set NMAXPOS,128
.set KKEY,96
.set OFFKEY,48
.set NF,32-CW
.set PBITS,38-PB
.set NF3,32-CW3
.set EBITS,38-EB
.set E3CELLSH,EBITS+1
.set DIRF1SH,NF+PBITS
.set SLOTN,NMAXPOS+1
.set NULLSLOT,NMAXPOS
.set LASTADD,SLOTN
.set SPANK,KKEY*2*SLOTN
.set REGBYTES,8*3*SPANK
.set NTWA,(NMAXPOS+2)/3+(NMAXPOS+1)/3
.set NTWS,(NMAXPOS+2)/3
.set NSA1,(NMAXPOS+1)/3+NMAXPOS/3
.set NSS1,(NMAXPOS+1)/3
.set NSS2,NMAXPOS/3
.if (KKEY-3) > (NMAXPOS-1)
.error "KKEY-3 must not exceed NMAXPOS-1"
.endif
.if (OFFKEY % 3) != 0
.error "OFFKEY must be a multiple of 3"
.endif
.if (KKEY % 6) != 0
.error "KKEY must be a multiple of 6"
.endif
.if ((KKEY*2*SLOTN) % 8) != 0
.error "KKEY*2*(NMAXPOS+1) must be a multiple of 8"
.endif
.if (NF+PBITS) > 30
.error "directory index field (NF+PBITS) must be <= 30 bits"
.endif
.if (1<<NF) > 1020
.error "2^NF must be <= 1020 (scratch bytemap)"
.endif
.if (KKEY*(NMAXPOS+SLOTN)) > 65535
.error "key constants must fit a 16-bit immediate"
.endif
.altmacro
.set _gid,0
.macro ROWDISP
ldp x2,x5,[x19],#16
ldr x6,[x19],#8
add x2,x2,x24
add x5,x5,x25
ldr x9,[x10,x6,lsl #3]
br x9
.endm
.macro EBODY
ldp x0,x6,[x2]
ldr x3,[x6,x5,lsl #3]
add x8,x17,x3,lsr #27
add x7,x0,x3,lsl #38
ldadd x1,x9,[x8]
str x7,[x9]
.endm
.macro ADV s
.if \s == 1
add x2,x2,x26
add x5,x5,x27
.elseif \s == 2
add x2,x2,x26,lsl #1
add x5,x5,x27,lsl #1
.elseif \s == 3
add x2,x2,x11
add x5,x5,x12
.endif
.endm
.macro E9CHK id
tbnz x9,#14,.Lpf_\id
.Lpc_\id:
.endm
.macro E9FLUSH id
.Lpf_\id:
sub x0,x9,#4,lsl #12
add x0,x0,#8
mov x3,#2046
ldr x7,[sp,#480]
stp x3,x7,[x0]
add x9,x0,#16
str x9,[x8]
b .Lpc_\id
.endm
.macro STEP f,s
EBODY
ADV \s
.if \f == 9
E9CHK %_gid
.set _gid,_gid + 1
.endif
.endm
.macro STUBS first,n
.set _c,\first
.rept \n
E9FLUSH %_c
.set _c,_c + 1
.endr
.endm
.macro TAIL ell
.if \ell == 0
ROWDISP
.else
b .Lp_top2
.endif
.endm
.macro LBLT f,ell,tw,k
.Lpt\f\ell\()_\tw\()_\k:
.endm
.macro TOWER f,ell,tw,n
.set _id0,_gid
.set _i,0
.rept \n
LBLT \f,\ell,\tw,%(\n - _i)
.if _i == \n - 1
STEP \f,0
.else
.if \tw == 0
.set _sd,1
.elseif \tw == 3
.set _sd,3
.elseif \tw == 1
.if ((\n - 2 - _i) & 1) == 0
.set _sd,1
.else
.set _sd,2
.endif
.else
.if ((\n - 2 - _i) & 1) == 0
.set _sd,2
.else
.set _sd,1
.endif
.endif
STEP \f,%_sd
.endif
.set _i,_i + 1
.endr
TAIL \ell
.if \f == 9
STUBS %_id0,\n
.endif
.endm
.macro LBLA f,ell,sk,fam,nk
.Lpa\f\ell\()_\sk\()_\fam\()_\nk:
.endm
.macro SA f,ell,sk,fam,nk
.set _id0,_gid
LBLA \f,\ell,\sk,\fam,\nk
ADV \sk
.set _i,0
.rept \nk
.if _i == \nk - 1
STEP \f,0
.else
.if \fam == 3
.set _sd,3
.elseif (_i & 1) == 0
.set _sd,1
.else
.set _sd,2
.endif
STEP \f,%_sd
.endif
.set _i,_i + 1
.endr
TAIL \ell
.if \f == 9
STUBS %_id0,\nk
.endif
.endm
.macro SAROW f,ell,sk,fam,nmax
.set _n,1
.rept \nmax
SA \f,\ell,\sk,\fam,%_n
.set _n,_n + 1
.endr
.endm
.macro LREFT f,ell,tw,k
.long .Lpt\f\ell\()_\tw\()_\k - ct_base
.endm
.macro LREFA f,ell,sk,fam,nk
.long .Lpa\f\ell\()_\sk\()_\fam\()_\nk - ct_base
.endm
.macro CHAFAM f,ell,fam
.set _n,1
.rept NMAXPOS
.if \fam == 0
LREFT \f,\ell,0,%_n
.elseif \fam == 3
.if _n <= NTWS
LREFT \f,\ell,3,%_n
.else
.long 0
.endif
.elseif \fam == 1
.if _n <= NTWA
.if (_n & 1) == 0
LREFT \f,\ell,1,%_n
.else
LREFT \f,\ell,2,%_n
.endif
.else
.long 0
.endif
.elseif \fam == 2
.if _n <= NTWA
.if (_n & 1) == 0
LREFT \f,\ell,2,%_n
.else
LREFT \f,\ell,1,%_n
.endif
.else
.long 0
.endif
.elseif \fam == 4
.if _n <= NSA1
LREFA \f,\ell,1,1,%_n
.else
.long 0
.endif
.elseif \fam == 5
.if _n <= NSS1
LREFA \f,\ell,1,3,%_n
.else
.long 0
.endif
.else
.if _n <= NSS2
LREFA \f,\ell,2,3,%_n
.else
.long 0
.endif
.endif
.set _n,_n + 1
.endr
.endm
.macro CHAFL f,ell
CHAFAM \f,\ell,0
CHAFAM \f,\ell,1
CHAFAM \f,\ell,2
CHAFAM \f,\ell,3
CHAFAM \f,\ell,4
CHAFAM \f,\ell,5
CHAFAM \f,\ell,6
.endm
.macro CHAINS f,ell
TOWER \f,\ell,0,NMAXPOS
TOWER \f,\ell,1,NTWA
TOWER \f,\ell,2,NTWA
TOWER \f,\ell,3,NTWS
SAROW \f,\ell,1,1,NSA1
SAROW \f,\ell,1,3,NSS1
SAROW \f,\ell,2,3,NSS2
.endm
.set O_K1,0
.set O_EK1,8
.set O_GW0,16
.set O_MU10,24
.set O_K0,32
.set O_DA1,40
.set O_DB1,48
.set O_LUT,56
.set O_DA0,64
.set O_DB0,72
.set O_MU21,80
.set O_MU20,88
.set O_DA2,96
.set O_DB2,104
.set O_K2,112
.set O_EK2,120
.set O_AE0,128
.set O_AE1,136
.set O_AE2,144
.set O_FLAV,152
.set O_A0,160
.set O_B0,168
.set O_PTRS,176
.set O_PTR0,184
.set O_DIR,192
.set O_E3,200
.set O_ARENA,208
.set O_ARENAEND,216
.set O_REC,224
.set O_P2,232
.set O_EMPTY,240
.set O_EMPTY3,248
.set O_SCR,256
.set O_GLB,264
.set O_STATUS,272
.set O_TNEG,280
.set O_TDIRP,288
.set O_TMP,296
.set O_SENT,304
.set O_SPSAVE,312
.set O_TPOS,320
.set O_TSLOW,328
.set O_FRAME,336
.set O_CELL0,344
.set O_CELL1,352
.set O_BIN0,360
.set O_BIN1,368
.set G_M3X,0
.set G_MB,8
.set G_TB8,16
.set G_TB9,24
.set G_LANE,32
.set G_SOFL,40
.set G_BAD,48
.set S_BM,0
.set S_BPL,2048
.set S_E3S,4096
.set S_X,262144
.set S_SIZE,524288
.set CELLW,1<<CW
.set NFV,1<<NF
.set PBMAX,(1<<PBITS)-1
.set E3BMAX,1<<EBITS
.set JMIN,-120
.set JMAX,120
.set MBT0,8192
.set MBS1,2048
.set M3XR,4096
.set M3XR0,(3-(M3XR%3))%3
.set MBR0,(3-(MBS1%3))%3
.set MBTM0,(3-(MBT0%3))%3
.set S_BMS,8192
.set G_SCR,56
.set CHK_K1,2*(MBS1-2)-120
.set CHK_KT,(4*(MBT0-3))/3-120
.if CHK_KT < CHK_K1
.set CHK_K1,CHK_KT
.endif
.set FR_3DA0,128
.set FR_MU,256
.set FR_DA1,272
.set FR_DA2,288
.set FR_E3,304
.set FR_M3X,312
.set FR_P2,320
.set FR_MB,328
.set FR_PTR0,416
.set FR_U3,3846
.set L3_AB,80
.set L3_D3,128
//@seg handoff3a
pt_handoff3:
ldp x8,x9,[x30,#L3_AB]
ldp x10,x13,[x30,#L3_D3]
madd x8,x0,x10,x8
madd x9,x0,x13,x9
lsr x10,x6,#EB
ubfx x13,x2,#CW3,#NF3
asr x3,x2,#32
bfi x10,x13,#E3CELLSH,#NF3
ldr x13,[sp,#FR_E3]
ldr x14,[x13,x10,lsl #3]
//@end
//@seg handoff3b
ldr x5,[x14],#8
add x5,x5,x3
//@end
//@seg handoff3c
pt_handoff3c:
ldp x10,x13,[sp,#FR_MU]
msub x21,x10,x5,x7
add x21,x21,x10
msub x22,x13,x5,x4
add x22,x22,x13
ldp x13,x15,[sp,#FR_DA2]
madd x20,x5,x13,x8
madd x23,x5,x15,x9
ldr x10,[sp,#FR_M3X]
ldrh w10,[x10,w5,sxtw #1]
ldrh w8,[sp,#FR_U3]
ldr x9,[sp,#FR_P2]
add x10,x10,x8
ldr x30,[x9,x10,lsl #3]
ldr x4,[sp,#FR_MB]
ldp x11,x12,[sp,#FR_3DA0]
mov x1,#8
//@end
//@seg node_L2
.Lp_top2:
ldr w6,[x14],#4
tbnz w6,#31,.Lp_exit
ldp x7,x8,[sp,#FR_MU]
sub x21,x21,x7
sub x22,x22,x8
asr x3,x21,#32
ubfx x7,x21,#CW,#NF
msub x8,x29,x3,x22
ubfx x9,x8,#CW,#NF
asr x5,x8,#32
bfi x6,x9,#0,#NF
bfi x6,x7,#DIRF1SH,#NF
ldr x19,[x28,x6,lsl #3]
ldp x8,x9,[sp,#FR_DA1]
madd x24,x3,x8,x20
madd x25,x3,x9,x23
madd x24,x5,x26,x24
madd x25,x5,x27,x25
add x20,x20,x13
add x23,x23,x15
add x9,x3,x5,lsl #12
ldrh w9,[x4,w9,sxtw #1]
ldr x10,[x30,x9]
ldr x30,[x30]
ldr x8,[sp,#FR_PTR0]
str x8,[x17]
//@end
//@seg row_disp
.Lp_row:
ldp x2,x5,[x19],#16
ldr x6,[x19],#8
add x2,x2,x24
add x5,x5,x25
ldr x9,[x10,x6,lsl #3]
br x9
//@end
//@seg l3_exit
.Lp_exit:
b pt_l3_next
//@end
//@seg null_nl
.Lp_null0:
ldp x2,x5,[x19],#16
ldr x6,[x19],#8
add x2,x2,x24
add x5,x5,x25
ldr x9,[x10,x6,lsl #3]
br x9
//@end
//@seg null_last
.Lp_null1:
b .Lp_top2
//@end
.ifdef PT_TEST
pt_l3_next:
nop
.endif
//@seg pair_init
pt_pair:
str x30,[x28,#O_SPSAVE]
str x1,[x28,#O_TMP]
ldp x20,x21,[x28,#O_K1]
ldp x22,x23,[x28,#O_GW0]
ldp x2,x25,[x28,#O_K0]
ldp x26,x24,[x28,#O_DB1]
ldr x15,[x28,#O_SCR]
add x3,x0,#1
lsl x3,x3,#PB
sub x3,x3,#1
mov x4,#(1<<38)
cmp x0,#PBMAX
csel x19,x4,x3,eq
lsl x3,x1,#DIRF1SH
orr x3,x3,x0,lsl #NF
ldr x4,[x28,#O_DIR]
add x4,x4,x3,lsl #3
str x4,[x28,#O_TDIRP]
mov x27,#1
lsl x27,x27,#32
mov x16,#NFV
mov x17,#1
mov x3,#S_X
add x9,x15,x3
lsl x0,x1,#CW
sub x5,x27,x2
sub x5,x5,#1
mov x3,#CELLW
add x6,x2,x3
sub x6,x6,#1
mov x7,#0
mov x8,#0
mov x12,#OFFKEY
mov x3,#SPANK
add x14,x12,x3
add x13,x14,x3
mov x10,#0x7fffffffffffffff
mov x11,#0x8000000000000000
//@end
//@seg gen_neg
.Lgn_top:
subs x2,x0,x20
csel x2,x2,xzr,hi
lsr x2,x2,#13
mul x2,x2,x2
umulh x2,x2,x21
subs x3,x19,x2
b.lo .Lgn_done
umulh x4,x3,x22
ldr x4,[x24,x4,lsl #3]
sub x1,x5,x4
add x2,x6,x4
cmp x1,x10
csel x10,x1,x10,lt
cmp x2,x11
csel x11,x2,x11,gt
ubfx x3,x1,#CW,#NF
sub x3,x16,x3
strb w17,[x15,x3]
ubfx x3,x2,#CW,#NF
sub x3,x16,x3
strb w17,[x15,x3]
stp x1,x2,[x9,#-40]!
stp x7,x8,[x9,#16]
str x12,[x9,#32]
add x0,x0,x27
add x5,x5,x23
add x6,x6,x23
sub x7,x7,x25
sub x8,x8,x26
mov x29,x12
mov x12,x13
mov x13,x14
mov x14,x29
b .Lgn_top
//@end
//@seg gn_to_gp
.Lgn_done:
str x9,[x28,#O_TNEG]
ldr x1,[x28,#O_TMP]
add x1,x1,#1
lsl x1,x1,#CW
sub x0,x27,x1
add x0,x0,#1
ldr x2,[x28,#O_K0]
sub x5,x27,x2
sub x5,x5,x23
sub x5,x5,#1
mov x3,#CELLW
add x6,x2,x3
sub x6,x6,x23
sub x6,x6,#1
mov x7,x25
mov x8,x26
mov x3,#S_X
add x9,x15,x3
mov x3,#SPANK
add x12,x3,#OFFKEY
add x13,x12,x3
mov x14,#OFFKEY
//@end
//@seg gen_pos
.Lgp_top:
subs x2,x0,x20
csel x2,x2,xzr,hi
lsr x2,x2,#13
mul x2,x2,x2
umulh x2,x2,x21
subs x3,x19,x2
b.lo .Lgp_done
umulh x4,x3,x22
ldr x4,[x24,x4,lsl #3]
sub x1,x5,x4
add x2,x6,x4
cmp x1,x10
csel x10,x1,x10,lt
cmp x2,x11
csel x11,x2,x11,gt
ubfx x3,x1,#CW,#NF
sub x3,x16,x3
strb w17,[x15,x3]
ubfx x3,x2,#CW,#NF
sub x3,x16,x3
strb w17,[x15,x3]
stp x1,x2,[x9],#40
stp x7,x8,[x9,#-24]
str x12,[x9,#-8]
add x0,x0,x27
sub x5,x5,x23
sub x6,x6,x23
add x7,x7,x25
add x8,x8,x26
mov x29,x12
mov x12,x13
mov x13,x14
mov x14,x29
b .Lgp_top
//@end
//@seg pair_post
.Lgp_done:
asr x2,x10,#32
asr x3,x11,#32
mov x5,#(KKEY-4-OFFKEY)
cmn x2,#OFFKEY
ccmp x3,x5,#0,ge
cset x25,gt
ldr x11,[x28,#O_TNEG]
ldr x20,[x28,#O_TDIRP]
ldr x24,[x28,#O_EMPTY]
mov x16,x9
str x16,[x28,#O_TPOS]
str x25,[x28,#O_TSLOW]
cmp x11,x16
b.eq .Lp_allempty
add x19,x15,#S_BPL
strh wzr,[x19],#2
mov x4,x15
mov x5,#NFV
add x6,x15,x5
//@end
//@seg bp_scan
.Lsc_w:
ldr x2,[x4]
cbz x2,.Lsc_n
str xzr,[x4]
sub x3,x4,x15
.Lsc_b:
rbit x7,x2
clz x7,x7
lsr x7,x7,#3
add x7,x7,x3
strh w7,[x19],#2
sub x8,x2,#1
ands x2,x2,x8
b.ne .Lsc_b
.Lsc_n:
add x4,x4,#8
cmp x4,x6
b.lo .Lsc_w
//@end
//@seg pair_arr_init
str xzr,[x4]
mov x5,#NFV
strh w5,[x19]
add x19,x15,#S_BPL
ldr x22,[x28,#O_ARENA]
ldp x26,x27,[x28,#O_DA0]
mov x12,#KKEY
mov x21,#(KKEY*LASTADD)
mov x13,#(NMAXPOS-1)
//@end
//@seg pair_arr
.Lp_arr:
ldrh w0,[x19],#2
ldrh w1,[x19]
lsl x10,x0,#CW
mov x17,x11
bl pt_rows
add x2,x20,x0,lsl #3
sub x3,x1,x0
tbz x3,#0,.Lf_even
str x7,[x2],#8
.Lf_even:
lsr x3,x3,#1
cbz x3,.Lf_done
.Lf_loop:
stp x7,x7,[x2],#16
subs x3,x3,#1
b.ne .Lf_loop
.Lf_done:
cmp x1,#NFV
b.ne .Lp_arr
str x22,[x28,#O_ARENA]
ldr x30,[x28,#O_SPSAVE]
ret
//@end
//@seg pair_empty
.Lp_allempty:
mov x3,#(NFV/2)
.Lp_ae:
stp x24,x24,[x20],#16
subs x3,x3,#1
b.ne .Lp_ae
ldr x30,[x28,#O_SPSAVE]
ret
//@end
//@seg rows_fast
pt_rows:
mov x23,x22
cbnz x25,.Lr_slow
.Lr_row:
ldp x2,x3,[x17],#16
ldp x4,x5,[x17],#16
ldr x6,[x17],#8
add x2,x2,x10
add x3,x3,x10
asr x2,x2,#32
asr x3,x3,#32
subs x7,x3,x2
b.lt .Lr_next
madd x4,x2,x26,x4
madd x5,x2,x27,x5
madd x6,x7,x12,x6
add x6,x6,x2
stp x4,x5,[x22],#16
str x6,[x22],#8
.Lr_next:
cmp x17,x16
b.lo .Lr_row
b .Lr_post
//@end
//@seg rows_slow
.Lr_slow:
.Ls_row:
ldp x2,x3,[x17],#16
ldp x4,x5,[x17],#16
ldr x6,[x17],#8
add x2,x2,x10
add x3,x3,x10
asr x2,x2,#32
asr x3,x3,#32
.Ls_piece:
subs x7,x3,x2
b.lt .Ls_next
cmp x7,x13
csel x7,x7,x13,ls
add x8,x2,#OFFKEY
cmp x8,#(KKEY-3)
csel x9,x2,xzr,ls
b.ls .Ls_key
mov x14,#3
sdiv x15,x2,x14
msub x9,x15,x14,x2
cmp x9,#0
add x15,x9,#3
csel x9,x15,x9,lt
.Ls_key:
madd x15,x7,x12,x6
add x15,x15,x9
madd x8,x2,x26,x4
madd x9,x2,x27,x5
stp x8,x9,[x22],#16
str x15,[x22],#8
add x2,x2,x7
add x2,x2,#1
b .Ls_piece
.Ls_next:
cmp x17,x16
b.lo .Ls_row
//@end
//@seg rows_post
.Lr_post:
cmp x22,x23
b.eq .Lr_none
ldr x6,[x22,#-8]
add x6,x6,x21
str x6,[x22,#-8]
mov x7,x23
ret
.Lr_none:
mov x7,x24
ret
//@end
//@seg e3_cell
pt_e3_cell:
ldp x2,x3,[x28,#O_K2]
ldr x15,[x28,#O_SCR]
add x15,x15,#(S_E3S>>12),lsl #12
mov x27,#1
lsl x27,x27,#32
lsl x4,x0,#CW3
add x5,x0,#1
lsl x5,x5,#CW3
sub x5,x5,#1
mov x6,#120
lsl x6,x6,#32
add x6,x4,x6
mvn x9,xzr
str x9,[x15,#56]
add x7,x15,#64
add x10,x15,#1032
//@end
//@seg e3_cell_neg
.Le3_neg:
subs x8,x6,x2
csel x8,x8,xzr,hi
lsr x8,x8,#13
mul x8,x8,x8
umulh x8,x8,x3
str x8,[x7],#8
sub x6,x6,x27
cmp x7,x10
b.lo .Le3_neg
//@end
//@seg e3_cell_mid
sub x6,x27,x5
add x10,x15,#1992
//@end
//@seg e3_cell_pos
.Le3_pos:
subs x8,x6,x2
csel x8,x8,xzr,hi
lsr x8,x8,#13
mul x8,x8,x8
umulh x8,x8,x3
str x8,[x7],#8
add x6,x6,x27
cmp x7,x10
b.lo .Le3_pos
//@end
//@seg e3_cell_end
str x9,[x7]
sub x8,x27,x5
cmp x4,x8
cset x0,hi
ret
//@end
//@seg e3_build_cell
pt_e3_build:
sub sp,sp,#48
str x30,[sp]
ldr x0,[x28,#O_CELL0]
str x0,[sp,#8]
.Le3_c:
bl pt_e3_cell
ldr x15,[x28,#O_SCR]
add x15,x15,#(S_E3S>>12),lsl #12
add x1,x15,#1024
add x1,x1,x0,lsl #3
sub x4,x1,#8
mov x5,x1
ldr x6,[sp,#8]
ldr x7,[x28,#O_E3]
lsl x8,x6,#E3CELLSH
add x7,x7,x8,lsl #3
ldr x22,[x28,#O_ARENA]
ldr x24,[x28,#O_EMPTY3]
mov x9,#0
mov x25,#PBMAX
mov x26,#1
lsl x26,x26,#38
//@end
//@seg e3_build_bin
.Le3_b:
add x10,x9,#1
lsl x10,x10,#EB
sub x10,x10,#1
cmp x9,#E3BMAX
csel x10,x26,x10,eq
.Le3_up:
ldr x11,[x5]
cmp x11,x10
b.hi .Le3_dn
add x5,x5,#8
b .Le3_up
.Le3_dn:
ldr x11,[x4]
cmp x11,x10
b.hi .Le3_em
sub x4,x4,#8
b .Le3_dn
.Le3_em:
add x12,x4,#8
sub x13,x5,#8
cmp x12,x13
b.hi .Le3_empty
mov x14,x22
sub x2,x12,x15
asr x2,x2,#3
sub x2,x2,#128
str x2,[x22],#8
//@end
//@seg e3_build_cand
.Le3_w:
ldr x2,[x12],#8
sub x2,x10,x2
lsr x2,x2,#PB
cmp x2,x25
csel x2,x2,x25,ls
lsl x2,x2,#NF
str w2,[x22],#4
cmp x12,x13
b.ls .Le3_w
//@end
//@seg e3_build_list_end
mov w2,#0x80000000
str w2,[x22],#4
add x22,x22,#7
and x22,x22,#-8
str x14,[x7],#8
b .Le3_nx
//@end
//@seg e3_build_empty
.Le3_empty:
str x24,[x7],#8
//@end
//@seg e3_build_next
.Le3_nx:
add x9,x9,#1
cmp x9,#E3BMAX
b.ls .Le3_b
str x22,[x28,#O_ARENA]
ldr x0,[sp,#8]
add x0,x0,#1
str x0,[sp,#8]
ldr x1,[x28,#O_CELL1]
cmp x0,x1
b.lo .Le3_c
ldr x30,[sp]
add sp,sp,#48
ret
//@end
//@seg ladd
pt_ladd:
mov x2,#0x9249
add x0,x0,x1
add x3,x0,x2
lsr x3,x3,#2
and x3,x3,x2
add x3,x3,x3,lsl #1
sub x0,x0,x3
ret
//@end
//@seg wheel_static
pt_wheel_static:
ldp x10,x11,[x0,#G_LANE]
ldr x12,[x0,#G_BAD]
mov x3,#3
mov x9,#0
.Lws_s:
mov x4,x9
mov x8,#0
mov x7,#0
.Lws_d:
udiv x5,x4,x3
msub x6,x5,x3,x4
lsl x6,x6,x7
orr x8,x8,x6
mov x4,x5
add x7,x7,#3
cmp x7,#18
b.lo .Lws_d
str w8,[x10,x9,lsl #2]
strh w9,[x11,x8,lsl #1]
lsr x5,x8,#9
orr x5,x8,x5
and x6,x5,#7
ubfx x7,x5,#3,#3
ubfx x4,x5,#6,#3
cmp x6,#0
ccmp x7,#0,#4,ne
ccmp x4,#0,#4,ne
cset x5,eq
strb w5,[x12,x9]
add x9,x9,#1
cmp x9,#729
b.lo .Lws_s
ret
//@end
//@seg bms_build
pt_bms_build:
str x30,[x28,#O_SPSAVE]
ldr x4,[x28,#O_GLB]
ldp x10,x11,[x4,#G_LANE]
ldr x12,[x4,#G_BAD]
ldr x5,[x28,#O_AE0]
ldr w6,[x10,x5,lsl #2]
mov x0,x6
mov x1,x6
bl pt_ladd
mov x7,x0
ldr x13,[x28,#O_SCR]
add x13,x13,#(S_BMS>>12),lsl #12
mov x9,#0
.Lb_s:
ldr w14,[x10,x9,lsl #2]
ldrb w15,[x12,x9]
mov x0,x14
mov x1,x6
bl pt_ladd
ldrh w0,[x11,x0,lsl #1]
ldrb w0,[x12,x0]
orr x15,x15,x0,lsl #1
mov x0,x14
mov x1,x7
bl pt_ladd
ldrh w0,[x11,x0,lsl #1]
ldrb w0,[x12,x0]
orr x15,x15,x0,lsl #2
strb w15,[x13,x9]
add x9,x9,#1
cmp x9,#729
b.lo .Lb_s
ldr x30,[x28,#O_SPSAVE]
ret
//@end
//@seg rec_build
pt_rec_build:
str x30,[x28,#O_SPSAVE]
ldr x4,[x28,#O_GLB]
ldp x10,x11,[x4,#G_LANE]
ldr x5,[x28,#O_AE1]
ldr w6,[x10,x5,lsl #2]
mov x0,x6
mov x1,x6
bl pt_ladd
mov x7,x0
ldr x5,[x28,#O_AE2]
ldr w8,[x10,x5,lsl #2]
ldr x13,[x28,#O_SCR]
add x13,x13,#(S_BMS>>12),lsl #12
ldr x14,[x28,#O_REC]
ldr x2,[x28,#O_FLAV]
add x3,x4,#G_TB8
ldr x15,[x3,x2,lsl #3]
movz x19,#(REGBYTES & 0xffff)
movk x19,#(REGBYTES >> 16),lsl #16
mov x9,#0
mov x20,x14
.Lr_u:
ldr w21,[x10,x9,lsl #2]
mov x0,x21
mov x1,x8
bl pt_ladd
ldrh w0,[x11,x0,lsl #1]
add x2,x0,x0,lsl #2
add x2,x14,x2,lsl #4
str x2,[x20]
ldrb w22,[x13,x9]
mov x0,x21
mov x1,x6
bl pt_ladd
ldrh w0,[x11,x0,lsl #1]
ldrb w23,[x13,x0]
mov x0,x21
mov x1,x7
bl pt_ladd
ldrh w0,[x11,x0,lsl #1]
ldrb w24,[x13,x0]
lsl x2,x22,#6
add x2,x2,x23,lsl #3
add x2,x2,x24
madd x3,x2,x19,x15
str x3,[x20,#8]
add x3,x3,#8
str x3,[x20,#16]
add x3,x3,#8
str x3,[x20,#24]
lsl x2,x23,#6
add x2,x2,x24,lsl #3
add x2,x2,x22
madd x3,x2,x19,x15
str x3,[x20,#32]
add x3,x3,#8
str x3,[x20,#40]
add x3,x3,#8
str x3,[x20,#48]
lsl x2,x24,#6
add x2,x2,x22,lsl #3
add x2,x2,x23
madd x3,x2,x19,x15
str x3,[x20,#56]
add x3,x3,#8
str x3,[x20,#64]
add x3,x3,#8
str x3,[x20,#72]
add x20,x20,#80
add x9,x9,#1
cmp x9,#729
b.lo .Lr_u
ldr x30,[x28,#O_SPSAVE]
ret
//@end
//@seg p2_build
pt_p2_build:
str x30,[x28,#O_SPSAVE]
ldr x4,[x28,#O_GLB]
ldp x10,x11,[x4,#G_LANE]
ldr x5,[x28,#O_AE2]
ldr w8,[x10,x5,lsl #2]
mov x0,x8
mov x1,x8
bl pt_ladd
mov x9,x0
ldr x14,[x28,#O_REC]
ldr x15,[x28,#O_P2]
mov x7,#0
mov x6,#0
.Lp2_c:
mov x5,#0
.Lp2_u:
ldr w21,[x10,x5,lsl #2]
mov x0,x21
mov x1,x7
bl pt_ladd
ldrh w0,[x11,x0,lsl #1]
add x2,x0,x0,lsl #2
add x2,x14,x2,lsl #4
str x2,[x15],#8
add x5,x5,#1
cmp x5,#729
b.lo .Lp2_u
add x6,x6,#1
cmp x6,#1
csel x7,x8,x9,eq
cmp x6,#3
b.lo .Lp2_c
ldr x30,[x28,#O_SPSAVE]
ret
//@end
//@seg m3x_build
pt_m3x_build:
ldr x0,[x0,#G_M3X]
sub x1,x0,#((2*M3XR)>>12),lsl #12
add x2,x0,#((2*M3XR)>>12),lsl #12
mov x5,#(M3XR0*729)
.Lm3_l:
strh w5,[x1],#2
add x5,x5,#729
cmp x5,#2187
csel x5,xzr,x5,eq
cmp x1,x2
b.lo .Lm3_l
ret
//@end
//@seg mb_blocks
pt_mb_build:
ldr x4,[x0,#G_MB]
ldr x1,[x0,#G_SCR]
movz x5,#((2*(4096*MBT0+2048)) & 0xffff)
movk x5,#((2*(4096*MBT0+2048)) >> 16),lsl #16
sub x0,x4,x5
mov x2,x1
mov x3,#0
.Lmb_blk:
mov x6,#(8+24*MBR0)
lsl x7,x3,#3
add x6,x6,x7
add x9,x7,#80
mov x5,#0
.Lmb_k:
strh w6,[x2],#2
add x6,x6,#24
sub x8,x6,#72
cmp x6,x9
csel x6,x8,x6,hs
add x5,x5,#1
cmp x5,#4096
b.lo .Lmb_k
add x3,x3,#1
cmp x3,#3
b.lo .Lmb_blk
//@end
//@seg mb_copy
mov x10,#(2*MBT0)
mov x3,#MBTM0
add x11,x1,#8192
add x12,x11,#8192
.Lmb_t:
cmp x3,#0
csel x13,x1,x11,eq
cmp x3,#2
csel x13,x12,x13,eq
add x14,x0,#8192
.Lmb_c:
ldp x15,x16,[x13],#16
stp x15,x16,[x0],#16
ldp x15,x16,[x13],#16
stp x15,x16,[x0],#16
ldp x15,x16,[x13],#16
stp x15,x16,[x0],#16
ldp x15,x16,[x13],#16
stp x15,x16,[x0],#16
ldp x15,x16,[x13],#16
stp x15,x16,[x0],#16
ldp x15,x16,[x13],#16
stp x15,x16,[x0],#16
ldp x15,x16,[x13],#16
stp x15,x16,[x0],#16
ldp x15,x16,[x13],#16
stp x15,x16,[x0],#16
cmp x0,x14
b.lo .Lmb_c
add x3,x3,#1
cmp x3,#3
csel x3,xzr,x3,eq
subs x10,x10,#1
b.ne .Lmb_t
ret
//@end
//@seg tb_chb3
pt_tb_chb3:
adr x12,famof_tab
mov x17,#255
mov x5,#0
.Lc3_m:
ldrb w7,[x12,x5]
eor x8,x5,#7
mov x6,#0
.Lc3_l:
cmp x6,#0
csel x9,x3,x4,eq
mov x10,#7
madd x10,x6,x10,x7
mov x11,#(NMAXPOS*4)
madd x10,x10,x11,x1
mov x13,#0
mov x14,#0
mov x15,#0
.Lc3_s:
lsr x16,x8,x14
and x16,x16,#1
add x13,x13,x16
add x14,x14,#1
cmp x14,#3
csel x14,xzr,x14,eq
cmp x13,#0
ccmp x7,x17,#4,ne
b.eq .Lc3_null
sub x16,x13,#1
ldrsw x16,[x10,x16,lsl #2]
add x16,x2,x16
str x16,[x0],#8
b .Lc3_n
.Lc3_null:
str x9,[x0],#8
.Lc3_n:
add x15,x15,#1
cmp x15,#NMAXPOS
b.lo .Lc3_s
str x9,[x0],#8
add x6,x6,#1
cmp x6,#2
b.lo .Lc3_l
add x5,x5,#1
cmp x5,#8
b.lo .Lc3_m
ret
//@end
//@seg tb_fb
pt_tb_fb:
mov x5,#0
mov x20,#(2*SLOTN*8)
.Lfb_b:
lsr x6,x5,#1
orr x6,x6,x5,lsl #2
and x6,x6,#7
lsr x7,x5,#2
orr x7,x7,x5,lsl #1
and x7,x7,#7
madd x8,x5,x20,x1
madd x9,x6,x20,x1
madd x10,x7,x20,x1
mov x11,#(2*SLOTN)
.Lfb_i:
ldr x12,[x8],#8
ldr x13,[x9],#8
ldr x14,[x10],#8
mov x15,#(KKEY/6)
.Lfb_k:
stp x12,x13,[x0],#16
stp x14,x12,[x0],#16
stp x13,x14,[x0],#16
subs x15,x15,#1
b.ne .Lfb_k
subs x11,x11,#1
b.ne .Lfb_i
add x5,x5,#1
cmp x5,#8
b.lo .Lfb_b
ret
//@end
//@seg tb_regions
pt_tb_regions:
movz x20,#((SPANK*8) & 0xffff)
movk x20,#((SPANK*8) >> 16),lsl #16
mov x5,#0
.Lrg_r:
lsr x6,x5,#6
ubfx x7,x5,#3,#3
and x8,x5,#7
mov x9,#0
.Lrg_b:
cmp x9,#0
csel x10,x6,x7,eq
cmp x9,#2
csel x10,x8,x10,eq
madd x11,x10,x20,x1
add x12,x0,x20
.Lrg_c:
ldp x13,x14,[x11],#16
stp x13,x14,[x0],#16
ldp x13,x14,[x11],#16
stp x13,x14,[x0],#16
ldp x13,x14,[x11],#16
stp x13,x14,[x0],#16
ldp x13,x14,[x11],#16
stp x13,x14,[x0],#16
cmp x0,x12
b.lo .Lrg_c
add x9,x9,#1
cmp x9,#3
b.lo .Lrg_b
add x5,x5,#1
cmp x5,#512
b.lo .Lrg_r
ret
//@end
//@seg tb_driver
pt_tb_build:
sub sp,sp,#48
str x30,[sp]
stp x0,x1,[sp,#8]
ldr x2,[x0,#G_SCR]
str x2,[sp,#24]
adr x3,ct_base
mov x4,#(14*NMAXPOS*4)
madd x1,x1,x4,x3
movz x4,#(28*NMAXPOS*4)
add x4,x3,x4
ldpsw x5,x6,[x4]
add x5,x3,x5
add x6,x3,x6
mov x0,x2
mov x2,x3
mov x3,x5
mov x4,x6
bl pt_tb_chb3
ldr x2,[sp,#24]
add x0,x2,#16,lsl #12
mov x1,x2
bl pt_tb_fb
ldp x3,x4,[sp,#8]
add x3,x3,#G_TB8
ldr x0,[x3,x4,lsl #3]
ldr x2,[sp,#24]
add x1,x2,#16,lsl #12
bl pt_tb_regions
ldr x30,[sp]
add sp,sp,#48
ret
//@end
//@seg global_init
pt_global_init:
sub sp,sp,#32
str x30,[sp]
str x0,[sp,#8]
bl pt_wheel_static
ldr x0,[sp,#8]
bl pt_m3x_build
ldr x0,[sp,#8]
bl pt_mb_build
ldr x0,[sp,#8]
mov x1,#0
bl pt_tb_build
ldr x0,[sp,#8]
mov x1,#1
bl pt_tb_build
ldr x30,[sp]
add sp,sp,#32
ret
//@end
//@seg orbit_check
pt_orbit_check:
ldr x1,[x28,#O_K1]
mov x2,#CHK_K1
cmp x1,x2
mov x0,#1
b.hi .Lck_ret
mov x2,#1
lsl x2,x2,#39
ldr x1,[x28,#O_EK1]
cmp x1,x2
mov x0,#3
b.lo .Lck_ret
ldr x1,[x28,#O_EK2]
cmp x1,x2
mov x0,#5
b.lo .Lck_ret
ldr x1,[x28,#O_GW0]
mov x0,#4
tbnz x1,#63,.Lck_ret
mov x0,#0
.Lck_ret:
ret
//@end
//@seg empty_init
pt_empty_init:
ldr x2,[x28,#O_EMPTY]
stp xzr,xzr,[x2]
mov x3,#(KKEY*(NULLSLOT+LASTADD))
str x3,[x2,#16]
ldr x2,[x28,#O_EMPTY3]
str xzr,[x2]
mov w3,#0x80000000
str w3,[x2,#8]
ldr x2,[x28,#O_SENT]
str w3,[x2]
ldr x2,[x28,#O_SCR]
mov x3,#(S_BPL/16)
.Lei_z:
stp xzr,xzr,[x2],#16
subs x3,x3,#1
b.ne .Lei_z
ret
//@end
//@seg build_dir
pt_build_dir:
sub sp,sp,#32
str x30,[sp]
ldr x2,[x28,#O_BIN0]
str x2,[sp,#8]
.Ld_bin:
mov x3,#0
str x3,[sp,#16]
.Ld_f1:
ldr x0,[sp,#8]
ldr x1,[sp,#16]
bl pt_pair
ldr x3,[sp,#16]
add x3,x3,#1
str x3,[sp,#16]
cmp x3,#NFV
b.lo .Ld_f1
ldr x2,[sp,#8]
add x2,x2,#1
str x2,[sp,#8]
ldr x4,[x28,#O_BIN1]
cmp x2,x4
b.lo .Ld_bin
ldr x30,[sp]
add sp,sp,#32
ret
//@end
//@seg orbit_build
pt_orbit_build:
sub sp,sp,#32
str x30,[sp]
bl pt_orbit_check
cbnz x0,.Lob_fail
bl pt_empty_init
bl pt_bms_build
bl pt_rec_build
bl pt_p2_build
bl pt_e3_build
bl pt_build_dir
mov x0,#0
.Lob_fail:
ldr x30,[sp]
add sp,sp,#32
ret
//@end
//@seg frame_init
pt_frame_init:
ldp x2,x3,[x28,#O_MU21]
stp x2,x3,[sp,#FR_MU]
ldp x2,x3,[x28,#O_DA1]
stp x2,x3,[sp,#FR_DA1]
ldp x2,x3,[x28,#O_DA2]
stp x2,x3,[sp,#FR_DA2]
ldr x2,[x28,#O_E3]
str x2,[sp,#FR_E3]
ldr x4,[x28,#O_GLB]
ldp x2,x3,[x4,#G_M3X]
str x2,[sp,#FR_M3X]
str x3,[sp,#FR_MB]
ldr x2,[x28,#O_P2]
str x2,[sp,#FR_P2]
ldr x2,[x28,#O_PTR0]
str x2,[sp,#FR_PTR0]
ldp x2,x3,[x28,#O_DA0]
add x2,x2,x2,lsl #1
add x3,x3,x3,lsl #1
stp x2,x3,[sp,#FR_3DA0]
strh wzr,[sp,#FR_U3]
//@end
//@seg pin_regs
pt_pin:
ldr x16,[x28,#O_LUT]
ldr x17,[x28,#O_PTRS]
ldp x26,x27,[x28,#O_DA0]
ldr x29,[x28,#O_MU10]
ldr x28,[x28,#O_DIR]
ret
//@end
//@seg entry_go
.Lent_go:
ldp x24,x25,[x28,#O_A0]
ldr x4,[x28,#O_REC]
ldr x10,[x4,#8]
ldr x14,[x28,#O_SENT]
mov x1,#8
ldp x2,x3,[x28,#O_DA0]
add x11,x2,x2,lsl #1
add x12,x3,x3,lsl #1
bl pt_pin
b .Lp_row
//@end
//@seg entry_m1
pt_entry_m1:
mov x28,x0
mov x0,#PBMAX
mov x1,#0
bl pt_pair
ldr x15,[x28,#O_SCR]
mov x3,#S_X
add x17,x15,x3
ldr x16,[x28,#O_TPOS]
ldr x25,[x28,#O_TSLOW]
mov x10,#0
ldr x22,[x28,#O_ARENA]
ldp x26,x27,[x28,#O_DA0]
mov x12,#KKEY
mov x21,#(KKEY*LASTADD)
mov x13,#(NMAXPOS-1)
ldr x24,[x28,#O_EMPTY]
cmp x17,x16
mov x7,x24
b.eq .Le1_skip
bl pt_rows
.Le1_skip:
str x22,[x28,#O_ARENA]
mov x19,x7
b .Lent_go
//@end
//@seg entry_m0
pt_entry_m0:
mov x28,x0
mov x0,#PBMAX
mov x1,#0
bl pt_pair
ldr x15,[x28,#O_SCR]
mov x3,#S_X
add x16,x15,x3
sub x17,x16,#40
ldr x2,[x28,#O_TNEG]
ldr x25,[x28,#O_TSLOW]
mov x10,#0
ldr x22,[x28,#O_ARENA]
ldp x26,x27,[x28,#O_DA0]
mov x12,#KKEY
mov x21,#(KKEY*LASTADD)
mov x13,#(NMAXPOS-1)
ldr x24,[x28,#O_EMPTY]
mov x4,#1
lsl x4,x4,#32
cmp x2,x16
mov x7,x24
b.eq .Le0_skip
str x4,[x17]
bl pt_rows
.Le0_skip:
str x22,[x28,#O_ARENA]
mov x19,x7
b .Lent_go
//@end
//@seg entry_m2
pt_entry_m2:
mov x28,x0
ldr x13,[x28,#O_E3]
mov x10,#E3BMAX
ldr x14,[x13,x10,lsl #3]
ldr x5,[x14],#8
mov x3,#1
subs x3,x3,x5
mov x6,#1
csel x5,x5,x6,le
b.le .Lm2_go
.Lm2_skip:
ldr w6,[x14]
tbnz w6,#31,.Lm2_go
add x14,x14,#4
subs x3,x3,#1
b.ne .Lm2_skip
.Lm2_go:
ldp x8,x9,[x28,#O_A0]
mov x7,#0
mov x4,#0
bl pt_pin
strh wzr,[sp,#FR_U3]
b pt_handoff3c
//@end
ct_base:
CHAFL 8,0
CHAFL 8,1
CHAFL 9,0
CHAFL 9,1
cha_null:
.long .Lp_null0 - ct_base
.long .Lp_null1 - ct_base
famof_tab:
.byte 0,4,2,6,1,5,3,255
.balign 4
CHAINS 8,0
CHAINS 8,1
CHAINS 9,0
CHAINS 9,1
.noaltmacro

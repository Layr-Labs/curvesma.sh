#define FE "adds x10,x10,x15\nadcs x11,x11,x16\nlsr x16,x9,#32\nadcs x12,x12,x17\n"
#define FS "subs x17,x9,x9,lsl #32\nsbc x1,x9,x16\nadds x9,x10,x9,lsl #32\n"
#define FQ FE "adcs x13,x13,x1\nadc x14,xzr,xzr\n" FS
#define FA(b) "mul x15,x3,x"#b"\nadcs x10,x11,x16\nmul x16,x4,x"#b"\nadcs x11,x12,x17\nmul x17,x5,x"#b"\nadcs x12,x13,x1\nmul x1,x6,x"#b"\n"
#define FB(b) "adds x9,x9,x15\numulh x15,x3,x"#b"\nadcs x10,x10,x16\numulh x16,x4,x"#b"\nadcs x11,x11,x17\numulh x17,x5,x"#b"\nadcs x12,x12,x1\n"
#define FR(b,c) FA(b) "adc x13,x"#c",xzr\n" FB(b) "umulh x1,x6,x"#b"\nadc x13,x13,xzr\n"
#define FC(b) FA(b) "umulh x14,x6,x"#b"\nadc x13,x14,xzr\n" FB(b) "adc x13,x13,xzr\n" FE "adc x13,x13,xzr\n" FS
#define FT(c) "adcs x10,x11,x16\nadcs x11,x12,x17\nadcs x12,x13,x1\nadc x13,x"#c",xzr\n"
#define F_M0(a) "ldp x7,x8," a "\nmul x9,x3,x7\numulh x15,x3,x7\nmul x10,x4,x7\numulh x16,x4,x7\nmul x11,x5,x7\numulh x17,x5,x7\nmul x12,x6,x7\n" \
"umulh x1,x6,x7\nadds x10,x10,x15\nadcs x11,x11,x16\nlsr x16,x9,#32\nadcs x12,x12,x17\nadc x13,xzr,x1\n" FS
#define F_MX(a,b) F_M0(a) FR(8,zr) "ldp x7,x8," b "\n" FQ FR(7,14) FQ FR(8,14) FQ FT(14)
#define F_MXC(a,b) F_M0(a) FC(8) "ldp x7,x8," b "\n" FC(7) FC(8) FT(zr)
#define FIN(a,b,c,d,e,k,l) "mov w"#k",#-1\nmov w"#l",#-2\n" FINB(a,b,c,d,e,k,l)
#define FINB(a,b,c,d,e,k,l) "adds x15,x"#a",#1\nsbcs x16,x"#b",x"#k"\nsbcs x17,x"#c",xzr\nadcs x1,x"#d",x"#l"\n" \
"sbcs xzr,x"#e",xzr\ncsel x"#a",x"#a",x15,lo\ncsel x"#b",x"#b",x16,lo\ncsel x"#c",x"#c",x17,lo\ncsel x"#d",x"#d",x1,lo\n"
#define F_MUL(m) F_MULX("[%["#m"]]","[%["#m"],#16]")
#define F_MULX(a,b) F_MX(a,b) FIN(9,10,11,12,13,7,8)
#define F_MAXC(a,b) F_MXC(a,b) F_MAXF
#define F_MAX(a,b) F_MX(a,b) F_MAXF
#define F_MAXF "neg x7,x13\nand x8,x7,#-0x100000000\nadds x9,x9,x13\nand x13,x7,#0xfffffffe\nadcs x10,x10,x8\nadcs x11,x11,x7\nadc x12,x12,x13\n"
#define SR "lsr x16,x7,#32\nadc x4,x1,xzr\nsubs x17,x7,x7,lsl #32\nsbc x1,x7,x16\nadds x7,x8,x7,lsl #32\nadcs x8,x3,x16\nadcs x3,x4,x17\n"
#define F_SQR "mul x8,x10,x9\numulh x16,x10,x9\nmul x3,x11,x9\numulh x17,x11,x9\nmul x4,x12,x9\numulh x5,x12,x9\nadds x3,x3,x16\nmul x15,x11,x10\n" \
"umulh x16,x11,x10\nadcs x4,x4,x17\nmul x17,x12,x10\numulh x1,x12,x10\nadc x5,x5,xzr\nmul x6,x12,x11\numulh x13,x12,x11\nadds x16,x16,x17\n" \
"mul x7,x9,x9\nadc x17,x1,xzr\nadds x4,x4,x15\numulh x9,x9,x9\nadcs x5,x5,x16\nmul x16,x10,x10\nadcs x6,x6,x17\numulh x10,x10,x10\n" \
"adc x13,x13,xzr\nadds x8,x8,x8\nmul x17,x11,x11\nadcs x3,x3,x3\numulh x11,x11,x11\nadcs x4,x4,x4\nmul x1,x12,x12\nadcs x5,x5,x5\n" \
"umulh x12,x12,x12\nadcs x6,x6,x6\nadcs x13,x13,x13\nadc x14,xzr,xzr\nadds x8,x8,x9\nadcs x3,x3,x16\nadcs x4,x4,x10\nadcs x5,x5,x17\n" \
"adcs x6,x6,x11\nadcs x13,x13,x1\nlsr x16,x7,#32\nadc x14,x14,x12\nsubs x17,x7,x7,lsl #32\nsbc x1,x7,x16\nadds x7,x8,x7,lsl #32\nadcs x8,x3,x16\n" \
"adcs x3,x4,x17\n" \
SR SR SR "adc x4,x1,xzr\nadds x7,x7,x5\nadcs x8,x8,x6\nadcs x3,x3,x13\nadcs x4,x4,x14\nadc x5,xzr,xzr\n" "ldp x9,x10,[%[s],#96]\n" FINB(7,8,3,4,5,9,10)
#define FK "ldp x7,x8,[%[s],#96]\n" FINB(9,10,11,12,13,7,8)
#define F_SUB(d0,d1,d2,d3,u0,u1,u2,u3,s0,s1,s2) "subs x"#d0",x"#d0",x"#u0"\nsbcs x"#d1",x"#d1",x"#u1"\n" \
"sbcs x"#d2",x"#d2",x"#u2"\nsbcs x"#d3",x"#d3",x"#u3"\nsbc x"#s0",xzr,xzr\nlsr x"#s1",x"#s0",#32\n" \
"and x"#s2",x"#s0",#-0xffffffff\nadds x"#d0",x"#d0",x"#s0"\nadcs x"#d1",x"#d1",x"#s1"\n" \
"adcs x"#d2",x"#d2",xzr\nadc x"#d3",x"#d3",x"#s2"\n"
#define F_LD(a,b,c,d,p) "ldp x"#a",x"#b",[%["#p"]]\nldp x"#c",x"#d",[%["#p"],#16]\n"
#define F_ST(a,b,c,d,p) "stp x"#a",x"#b",[%["#p"]]\nstp x"#c",x"#d",[%["#p"],#16]\n"
#define F_CL "x1","x3","x4","x5","x6","x7","x8","x9","x10","x11","x12","x13","x14","x15","x16","x17","cc","memory"
#define U(n) [n] "+r"(*n)
static inline void k_submul(uintptr_t *d, uintptr_t *r, const fe *b, uintptr_t *c) {
    __asm__ volatile(F_LD(3,4,5,6,b) "ldp x9,x10,[%[c],#16]\nldp x7,x8,[%[c]],#64\n" F_SUB(3,4,5,6,7,8,9,10,9,10,11)
        "stp x5,x6,[%[d],#16]\nstp x3,x4,[%[d]],#64\n" F_MAXC("[%[r],#-64]","[%[r],#-48]")
        "stp x11,x12,[%[r],#16]\nstp x9,x10,[%[r]],#64\n" : U(d), U(r), U(c) : [b] "r"(b->v) : F_CL);
}
static inline void k_subm(fe *r, const fe *b, const fe *c, const fe *m) {
    __asm__ volatile(F_LD(3,4,5,6,b) F_LD(7,8,9,10,c) F_SUB(3,4,5,6,7,8,9,10,9,10,11)
        F_MUL(m) F_ST(9,10,11,12,r) :: [r] "r"(r->v), [b] "r"(b->v), [c] "r"(c->v), [m] "r"(m->v) : F_CL);
}
static inline void k_mul2(fe *s, uintptr_t *m, uintptr_t *n) {
    __asm__ volatile(F_LD(3,4,5,6,s) F_MAX("[%[m],#-64]!","[%[m],#16]") "stp x9,x10,[%[s],#32]\nstp x11,x12,[%[s],#48]\n"
        "ldp x3,x4,[%[n],#-64]!\nldp x5,x6,[%[n],#16]\n" F_MAXC("[%[s]]","[%[s],#16]") F_ST(9,10,11,12,s) : U(m), U(n) : [s] "r"(s->v) : F_CL);
}
static inline void k_lxy(fe *s, const fe *R, const fe *c, fe *n) {
    __asm__ volatile("ldp x3,x4,[%[R],#32]\nldp x5,x6,[%[R],#48]\nldp x7,x8,[%[c],#32]\nldp x9,x10,[%[c],#48]\n"
        F_SUB(3,4,5,6,7,8,9,10,9,10,11) F_MXC("[%[s],#32]","[%[s],#48]") FK "stp x9,x10,[%[s],#64]\nstp x11,x12,[%[s],#80]\n"
        F_SQR F_LD(9,10,11,12,c)
        F_SUB(7,8,3,4,9,10,11,12,15,16,17) F_LD(5,6,13,14,R) F_SUB(7,8,3,4,5,6,13,14,15,16,17) F_ST(7,8,3,4,n)
        "subs x13,x9,x7\nsbcs x14,x10,x8\nsbcs x5,x11,x3\nsbcs x6,x12,x4\nsbc x15,xzr,xzr\nlsr x16,x15,#32\n"
        "and x17,x15,#-0xffffffff\nadds x3,x13,x15\nadcs x4,x14,x16\nadcs x5,x5,xzr\nadc x6,x6,x17\n"
        F_MXC("[%[s],#64]","[%[s],#80]") FK "ldp x13,x14,[%[c],#32]\nldp x15,x16,[%[c],#48]\n" F_SUB(9,10,11,12,13,14,15,16,3,4,5)
        "stp x9,x10,[%[n],#32]\nstp x11,x12,[%[n],#48]\n" :: [s] "r"(s->v), [R] "r"(R->v), [c] "r"(c->v), [n] "r"(n->v) : F_CL);
}
#undef U

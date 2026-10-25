/* Fused batch_step kernels: the field operations of the k_mul, k_sqr, k_subb_bl and k_sub2_bl calls
 * they replace, in order, with values in registers; no branch. F_MUL(m): x9..x12 = x3..x6*[m],
 * k_mul's 111 instructions, final constants in x7, x8 so x3..x6 survive; F_SQR: x7..x10 = x3..x6^2,
 * k_sqr's 87; F_SUB: d = d-u mod p, k_subb_bl's 12; FQ, FR, FIN, SR: repeated runs. */
#define FQ "adds x10,x10,x15\nadcs x11,x11,x16\nlsr x16,x9,#32\nadcs x12,x12,x17\nadcs x13,x13,x1\nadc x14,xzr,xzr\nsubs x17,x9,x9,lsl #32\nsbc x1,x9,x16\nadds x9,x10,x9,lsl #32\n"
#define FR(b,c) "mul x15,x3,x"#b"\nadcs x10,x11,x16\nmul x16,x4,x"#b"\nadcs x11,x12,x17\nmul x17,x5,x"#b"\nadcs x12,x13,x1\nmul x1,x6,x"#b"\nadc x13,x"#c",xzr\nadds x9,x9,x15\n" \
"umulh x15,x3,x"#b"\nadcs x10,x10,x16\numulh x16,x4,x"#b"\nadcs x11,x11,x17\numulh x17,x5,x"#b"\nadcs x12,x12,x1\numulh x1,x6,x"#b"\nadc x13,x13,xzr\n"
#define FIN(a,b,c,d,e,k,l) "mov w"#k",#-1\nmov w"#l",#-2\nadds x15,x"#a",#1\nsbcs x16,x"#b",x"#k"\nsbcs x17,x"#c",xzr\nadcs x1,x"#d",x"#l"\n" \
"sbcs xzr,x"#e",xzr\ncsel x"#a",x"#a",x15,lo\ncsel x"#b",x"#b",x16,lo\ncsel x"#c",x"#c",x17,lo\ncsel x"#d",x"#d",x1,lo\n"
#define F_MUL(m) "ldp x7,x8,[%["#m"]]\nmul x9,x3,x7\numulh x15,x3,x7\nmul x10,x4,x7\numulh x16,x4,x7\nmul x11,x5,x7\numulh x17,x5,x7\nmul x12,x6,x7\n" \
"umulh x1,x6,x7\nadds x10,x10,x15\nadcs x11,x11,x16\nlsr x16,x9,#32\nadcs x12,x12,x17\nadc x13,xzr,x1\nsubs x17,x9,x9,lsl #32\nsbc x1,x9,x16\n" \
"adds x9,x10,x9,lsl #32\n" \
FR(8,zr) "ldp x7,x8,[%["#m"],#16]\n" FQ FR(7,14) FQ FR(8,14) FQ "adcs x10,x11,x16\nadcs x11,x12,x17\nadcs x12,x13,x1\nadc x13,x14,xzr\n" \
FIN(9,10,11,12,13,7,8)
#define SR "lsr x16,x7,#32\nadc x10,x1,xzr\nsubs x17,x7,x7,lsl #32\nsbc x1,x7,x16\nadds x7,x8,x7,lsl #32\nadcs x8,x9,x16\nadcs x9,x10,x17\n"
#define F_SQR "mul x8,x4,x3\numulh x16,x4,x3\nmul x9,x5,x3\numulh x17,x5,x3\nmul x10,x6,x3\numulh x11,x6,x3\nadds x9,x9,x16\nmul x15,x5,x4\n" \
"umulh x16,x5,x4\nadcs x10,x10,x17\nmul x17,x6,x4\numulh x1,x6,x4\nadc x11,x11,xzr\nmul x12,x6,x5\numulh x13,x6,x5\nadds x16,x16,x17\n" \
"mul x7,x3,x3\nadc x17,x1,xzr\nadds x10,x10,x15\numulh x3,x3,x3\nadcs x11,x11,x16\nmul x16,x4,x4\nadcs x12,x12,x17\numulh x4,x4,x4\n" \
"adc x13,x13,xzr\nadds x8,x8,x8\nmul x17,x5,x5\nadcs x9,x9,x9\numulh x5,x5,x5\nadcs x10,x10,x10\nmul x1,x6,x6\nadcs x11,x11,x11\n" \
"umulh x6,x6,x6\nadcs x12,x12,x12\nadcs x13,x13,x13\nadc x14,xzr,xzr\nadds x8,x8,x3\nadcs x9,x9,x16\nadcs x10,x10,x4\nadcs x11,x11,x17\n" \
"adcs x12,x12,x5\nadcs x13,x13,x1\nlsr x16,x7,#32\nadc x14,x14,x6\nsubs x17,x7,x7,lsl #32\nsbc x1,x7,x16\nadds x7,x8,x7,lsl #32\nadcs x8,x9,x16\n" \
"adcs x9,x10,x17\n" \
SR SR SR "adc x10,x1,xzr\nadds x7,x7,x11\nadcs x8,x8,x12\nadcs x9,x9,x13\nadcs x10,x10,x14\nadc x11,xzr,xzr\n" FIN(7,8,9,10,11,3,4)
#define F_SUB(d0,d1,d2,d3,u0,u1,u2,u3,s0,s1,s2) "subs x"#d0",x"#d0",x"#u0"\nsbcs x"#d1",x"#d1",x"#u1"\n" \
"sbcs x"#d2",x"#d2",x"#u2"\nsbcs x"#d3",x"#d3",x"#u3"\nsbc x"#s0",xzr,xzr\nlsr x"#s1",x"#s0",#32\n" \
"lsl x"#s2",x"#s0",#32\nsub x"#s2",x"#s2",x"#s0"\nadds x"#d0",x"#d0",x"#s0"\nadcs x"#d1",x"#d1",x"#s1"\n" \
"adcs x"#d2",x"#d2",xzr\nadc x"#d3",x"#d3",x"#s2"\n"
#define F_LD(a,b,c,d,p) "ldp x"#a",x"#b",[%["#p"]]\nldp x"#c",x"#d",[%["#p"],#16]\n"
#define F_ST(a,b,c,d,p) "stp x"#a",x"#b",[%["#p"]]\nstp x"#c",x"#d",[%["#p"],#16]\n"
#define F_CL "x1","x3","x4","x5","x6","x7","x8","x9","x10","x11","x12","x13","x14","x15","x16","x17","cc","memory"
#define R(n) [n] "r"(n->v)
static inline void k_submul(fe *d, fe *r, const fe *b, const fe *c, const fe *m) {
    __asm__ volatile(F_LD(3,4,5,6,b) F_LD(7,8,9,10,c) F_SUB(3,4,5,6,7,8,9,10,9,10,11) F_ST(3,4,5,6,d)
        F_MUL(m) F_ST(9,10,11,12,r) :: R(d), R(r), R(b), R(c), R(m) : F_CL);
}
static inline void k_subm(fe *r, const fe *b, const fe *c, const fe *m) {
    __asm__ volatile(F_LD(3,4,5,6,b) F_LD(7,8,9,10,c) F_SUB(3,4,5,6,7,8,9,10,9,10,11)
        F_MUL(m) F_ST(9,10,11,12,r) :: R(r), R(b), R(c), R(m) : F_CL);
}
/* r = a*m; s = a*n (a is read before s is written) */
static inline void k_mul2(fe *r, fe *s, const fe *a, const fe *m, const fe *n) {
    __asm__ volatile(F_LD(3,4,5,6,a) F_MUL(m) F_ST(9,10,11,12,r) F_MUL(n) F_ST(9,10,11,12,s)
        :: R(r), R(s), R(a), R(m), R(n) : F_CL);
}
static inline void k_xy(fe *x, fe *y, const fe *m, const fe *b, const fe *c, const fe *e) {
    __asm__ volatile(F_LD(3,4,5,6,m) F_SQR F_LD(3,4,5,6,b) F_SUB(7,8,9,10,3,4,5,6,15,16,17)
        F_LD(11,12,13,14,c) F_SUB(7,8,9,10,11,12,13,14,15,16,17) F_ST(7,8,9,10,x)
        F_SUB(3,4,5,6,7,8,9,10,11,12,13) F_MUL(m) F_LD(13,14,15,16,e) F_SUB(9,10,11,12,13,14,15,16,3,4,5)
        F_ST(9,10,11,12,y) :: R(x), R(y), R(m), R(b), R(c), R(e) : F_CL);
}
#undef R

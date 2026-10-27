# Backward (BLOOP_FIELD_B, BLOOP_SLOW_B) and forward (FLOOP) blocks with branchy, metered corrections;
# inv canonical, PRE carry-only (the C code blinds each round's prefix chain with a fresh unit rho).
# One step of each block as a macro with a label suffix S; step.c composes eight copies per loop pass.
# Block constants: x28 = 0xffffffff (p1), x29 = 0xffffffff00000001 (p3); x27 = correction meter.
# usage: gen11.py <predecessor step.c> outdir [site]   (site: count only that site's events, weight 1)
# The predecessor (submission 778628fb, commit 1fcbf53a, measurement/step.c) defines BLOOP_FIELD.
import sys, re
_src = open(sys.argv[1]).read()
_m = re.search(r'#define BLOOP_FIELD \\\n(.*?)\n""\n', _src, re.S)
L = [l.strip() for l in "".join(re.findall(r'"((?:[^"\\]|\\.)*)"', _m.group(1))).split("\\n") if l.strip()]
assert len(L) == 574
outdir = sys.argv[2]
only = sys.argv[3] if len(sys.argv) > 3 else None
SITES = []
def wt(site, w):
    SITES.append((site, w))
    if only is None: return w
    return 1 if site == only else 0
lab = [0]
def nl():
    lab[0] += 1; return lab[0]
def sub_fix(r, site):  # after a borrowing 4-limb subtraction: r += p
    n = nl()
    return ["b.cs Lq%d_%%=" % n, "subs %s,%s,#1" % (r[0], r[0]), "adcs %s,%s,x28" % (r[1], r[1]),
            "adcs %s,%s,xzr" % (r[2], r[2]), "adc %s,%s,x29" % (r[3], r[3]), "add x27,x27,#%d" % wt(site, 5), "Lq%d_%%=:" % n]
def carry_fix(r, c, site):  # T = r + c*2^256, c in {0,1}, T < 2^256 + p: subtract c*p
    n = nl()
    return ["cbz %s,Lq%d_%%=" % (c, n), "adds %s,%s,#1" % (r[0], r[0]), "sbcs %s,%s,x28" % (r[1], r[1]),
            "sbcs %s,%s,xzr" % (r[2], r[2]), "sbc %s,%s,x29" % (r[3], r[3]), "add x27,x27,#%d" % wt(site, 5), "Lq%d_%%=:" % n]
slow = []
def slow_blk(n, r, s, site):
    slow.extend(["Ls%d_%%=:" % n, "adds %s,%s,#1" % (s[0], r[0]), "sbcs %s,%s,x28" % (s[1], r[1]),
                 "sbcs %s,%s,xzr" % (s[2], r[2]), "sbcs %s,%s,x29" % (s[3], r[3]),
                 "csel %s,%s,%s,hs" % (r[0], s[0], r[0]), "csel %s,%s,%s,hs" % (r[1], s[1], r[1]),
                 "csel %s,%s,%s,hs" % (r[2], s[2], r[2]), "csel %s,%s,%s,hs" % (r[3], s[3], r[3]),
                 "add x27,x27,#%d" % wt(site + "s", 10), "b Lq%d_%%=" % n])
def canon_fix(r, c, s, site):  # T = r + c*2^256 < 2p: canonical
    n = nl()
    code = ["cbz %s,Lp%d_%%=" % (c, n), "adds %s,%s,#1" % (r[0], r[0]), "sbcs %s,%s,x28" % (r[1], r[1]),
            "sbcs %s,%s,xzr" % (r[2], r[2]), "sbc %s,%s,x29" % (r[3], r[3]), "add x27,x27,#%d" % wt(site, 4),
            "b Lq%d_%%=" % n, "Lp%d_%%=:" % n, "cmn %s,x28" % r[3], "b.hs Ls%d_%%=" % n, "Lq%d_%%=:" % n]
    slow_blk(n, r, s, site)
    return code
def canon2_fix(r, c, s, site):  # T = r + c*2^256 < 2^256 + p: canonical (after a fix, re-test)
    n = nl()
    code = ["cbz %s,Lp%d_%%=" % (c, n), "adds %s,%s,#1" % (r[0], r[0]), "sbcs %s,%s,x28" % (r[1], r[1]),
            "sbcs %s,%s,xzr" % (r[2], r[2]), "sbc %s,%s,x29" % (r[3], r[3]), "add x27,x27,#%d" % wt(site, 5),
            "Lp%d_%%=:" % n, "cmn %s,x28" % r[3], "b.hs Ls%d_%%=" % n, "Lq%d_%%=:" % n]
    slow_blk(n, r, s, site)
    return code
def rename(lines, m):
    rx = re.compile(r'\bx(\d+)\b')
    return [rx.sub(lambda mo: m.get('x' + mo.group(1), 'x' + mo.group(1)), l) for l in lines]
# canonical-name FC product: a = x3..x6, window x9..x13, temps x14..x17, x1; b regs given
FS = ["subs x17,x9,x9,lsl #32", "sbc x1,x9,x16", "adds x9,x10,x9,lsl #32"]
def M0(b):
    return ["mul x9,x3,%s" % b, "umulh x15,x3,%s" % b, "mul x10,x4,%s" % b, "umulh x16,x4,%s" % b,
            "mul x11,x5,%s" % b, "umulh x17,x5,%s" % b, "mul x12,x6,%s" % b, "umulh x1,x6,%s" % b,
            "adds x10,x10,x15", "adcs x11,x11,x16", "lsr x16,x9,#32", "adcs x12,x12,x17", "adc x13,xzr,x1"] + FS
def FC(b): return (["mul x15,x3,%s" % b, "adcs x10,x11,x16", "mul x16,x4,%s" % b, "adcs x11,x12,x17",
                   "mul x17,x5,%s" % b, "adcs x12,x13,x1", "mul x1,x6,%s" % b, "umulh x14,x6,%s" % b, "adc x13,x14,xzr",
                   "adds x9,x9,x15", "umulh x15,x3,%s" % b, "adcs x10,x10,x16", "umulh x16,x4,%s" % b,
                   "adcs x11,x11,x17", "umulh x17,x5,%s" % b, "adcs x12,x12,x1", "adc x13,x13,xzr",
                   "adds x10,x10,x15", "adcs x11,x11,x16", "lsr x16,x9,#32", "adcs x12,x12,x17", "adc x13,x13,xzr"] + FS)
FT = ["adcs x10,x11,x16", "adcs x11,x12,x17", "adcs x12,x13,x1", "adc x13,xzr,xzr"]
def at(i, s): assert L[i-1] == s, (i, L[i-1], s)
R0 = ["x0", "x2", "x18", "x30"]; INV = ["x19", "x20", "x21", "x22"]; X36 = ["x3", "x4", "x5", "x6"]
TMP = ["x15", "x16", "x17", "x1"]
out = []
def seg(name, start):  # segment marks (used by the block tests)
    pass
_s = 0
# 1. iw = inv*PRE[w] with inv < p: five-word rows, window x0,x2,x18,x30,x13; PRE words loaded into x7,x8
iw = ["ldp x7,x8,[x23,#32]"] + M0("x7") + FC("x8") + ["ldp x7,x8,[x23,#48]"] + FC("x7") + FC("x8") + FT
out += rename(iw[:1], {}) + rename(iw[1:], {"x3": "x19", "x4": "x20", "x5": "x21", "x6": "x22",
                                           "x9": "x0", "x10": "x2", "x11": "x18", "x12": "x30"})
out += carry_fix(R0, "x13", "iw")
seg('iw', _s); _s = len(out)
# 2. inv = D[w]*inv (lines 108..199 of the measured text), canonical, written to x19..x22
at(108, "ldp x5,x6,[x23,#16]"); at(200, "adds x9,x10,x9,lsl #32")
out += L[107:199]
out += ["adds x19,x10,x9,lsl #32", "adcs x20,x11,x16", "adcs x21,x12,x17", "adcs x22,x13,x1", "adc x13,xzr,xzr"]
out += canon_fix(INV, "x13", TMP, "inv")
seg('inv', _s); _s = len(out)
# 3. dy
at(212, "ldr x7,[x26,#8]"); at(221, "sbc x13,xzr,xzr")
out += L[211:220] + sub_fix(X36, "dy")
seg('dy', _s); _s = len(out)
# 4. lambda = dy*iw, carry-only into x0,x2,x18,x30
at(228, "mul x9,x3,x0"); at(318, "adds x9,x10,x9,lsl #32")
out += L[227:317] + ["adds x0,x10,x9,lsl #32", "adcs x2,x11,x16", "adcs x18,x12,x17", "adcs x30,x13,x1", "adc x13,xzr,xzr"]
out += carry_fix(R0, "x13", "lam")
seg('lam', _s); _s = len(out)
# 5. lambda^2 (lambda < 2^256): canonical with re-test
at(332, "mul x4,x2,x0"); at(407, "adc x11,xzr,xzr")
out += L[331:407] + canon2_fix(X36, "x11", TMP, "sq")
seg('sq', _s); _s = len(out)
# 6. x3 = lambda^2 - x(Wc) - x(R)
at(417, "ldp x8,x9,[x24,#-32]"); at(423, "sbc x15,xzr,xzr")
out += L[416:422] + sub_fix(X36, "x3a")
at(430, "ldp x12,x13,[x7]"); at(436, "sbc x15,xzr,xzr")
out += L[429:435] + sub_fix(X36, "x3b")
seg('x3', _s); _s = len(out)
# 7. t = x(Wc) - x3
at(443, "stp x3,x4,[x25,#-48]"); at(449, "sbc x12,xzr,xzr")
out += L[442:448] + sub_fix(["x8", "x9", "x10", "x11"], "t")
seg('t', _s); _s = len(out)
# 8. t*lambda canonical into x0,x2,x18,x30
at(456, "mul x5,x8,x0"); at(546, "adds x5,x6,x5,lsl #32")
out += L[455:545] + ["adds x0,x6,x5,lsl #32", "adcs x2,x7,x16", "adcs x18,x12,x17", "adcs x30,x13,x1", "adc x13,xzr,xzr"]
out += canon_fix(R0, "x13", TMP, "y3m")
seg('y3m', _s); _s = len(out)
# 9. y3
at(560, "ldp x10,x11,[x24,#16]"); at(566, "sbc x12,xzr,xzr")
out += L[559:565] + sub_fix(R0, "y3")
at(573, "stp x0,x2,[x25,#-16]"); out += L[572:574]
seg('y3', _s)
def copy2(lines): return [re.sub(r'\bL([pqs])(\d+)_%=', r'L\1\2b_%=', l) for l in lines]
def cstr(lines): return '"' + "\\n".join(lines) + "\\n" + '"'
# forward block: x23 = &L[w].sel, x24 = &Wc[w], x25 = &DP[w], x26 = count, PRE in x19..x22
prod = M0("x19") + FC("x20") + FC("x21") + FC("x22")
prod = prod[:-1] + ["adds x19,x10,x9,lsl #32", "adcs x20,x11,x16", "adcs x21,x12,x17", "adcs x22,x13,x1", "adc x13,xzr,xzr"]
fw = (["1:", "ldr x7,[x23],#40", "ldp x3,x4,[x7]", "ldp x5,x6,[x7,#16]", "ldp x14,x15,[x24,#16]", "ldp x7,x8,[x24],#64",
       "subs x3,x3,x7", "sbcs x4,x4,x8", "sbcs x5,x5,x14", "sbcs x6,x6,x15"] + sub_fix(X36, "fd") +
      ["stp x5,x6,[x25,#16]", "stp x3,x4,[x25],#64"] + prod + carry_fix(INV, "x13", "fpre") +
      ["stp x19,x20,[x25,#32]", "stp x21,x22,[x25,#48]"])
FW1 = fw[1:]  # one step, after the "1:" label
with open(outdir + "/bloop_branchy.h", "w") as f:
    f.write("/* GENERATED by gen11.py: backward block with branchy corrections, meter in x27. */\n")
    def par(lines):  # one copy, labels suffixed by the macro argument S
        return cstr([re.sub(r'\bL([pqs])(\d+)_%=', r'L\1\2" S "_%=', l) for l in lines])
    f.write("#define BLOOP_FIELD_S(S) %s\n#define BLOOP_SLOW_S(S) %s\n" % (par(out), par(slow)))
with open(outdir + "/floop.h", "w") as f:
    fp = lambda lines: cstr([re.sub(r'\bL([pqs])(\d+)_%=', r'L\1\2" S "_%=', l) for l in lines])
    f.write("/* GENERATED by gen11.py: one forward step (prefix product), PRE in x19..x22. */\n#define FLOOP_S(S) %s\n" % fp(FW1))
if only is None:
    print("sites", SITES)

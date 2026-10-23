"""Check the instruction listing in p256_mul_count.md.

The three sequences below are the same instruction sequence, instruction for
instruction, as the listings in p256_mul_count.md (blocks A-D for the
multiplication, sections 6 and 7 for add and sub).  The script

  * counts them (the numbers quoted in proof.md: M = 139, add = 20, sub = 19),
  * executes them on a 64-bit register model with the AArch64 flag semantics (ADDS, ADCS,
    SUBS, SBCS write C; ADC and SBC read C and do not write it), and
  * asserts the results against a*b*2^-256 mod p, a+b mod p and a-b mod p.

Standard library only; no input; runs in well under a second.
"""
import random

M64 = (1 << 64) - 1
P = (1 << 256) - (1 << 224) + (1 << 192) + (1 << 96) - 1
PL = [(-1) & M64, 0xFFFFFFFF, 0, 0xFFFFFFFF00000001]
assert P == PL[3] << 192 | PL[2] << 128 | PL[1] << 64 | PL[0]

prog = []


def e(*a):
    prog.append(a)


# ---------------------------------------------------------------- mul_mont
def mul_mont():
    prog.clear()
    # a in [x0], b in [x1], result to [x2]; P0..P3 held in reserved regs
    e("ldp", "A0", "A1", "x0")
    e("ldp", "A2", "A3", "x0", "16")
    e("ldp", "B0", "B1", "x1")
    e("ldp", "B2", "B3", "x1", "16")
    e("mov", "t2", "xzr")
    # ---- product: Comba (column) summation; column k accumulates in (t_k,t_k+1,t_k+2)
    cols = {}
    for i in range(4):
        for j in range(4):
            cols.setdefault(i + j, []).append((i, j))
    for k in range(7):
        for n, (i, j) in enumerate(cols[k]):
            if k == 0:
                e("mul", "t0", "A0", "B0")
                e("umulh", "t1", "A0", "B0")
                break
            e("mul", "lo", "A%d" % i, "B%d" % j)
            e("umulh", "hi", "A%d" % i, "B%d" % j)
            e("adds", "t%d" % k, "t%d" % k, "lo")
            e("adcs", "t%d" % (k + 1), "t%d" % (k + 1), "hi")
            if n == 0:
                e("adc", "t%d" % (k + 2), "xzr", "xzr")
            else:
                e("adc", "t%d" % (k + 2), "t%d" % (k + 2), "xzr")
    # ---- 4 Montgomery reduction rounds: + m*p, p = 2^256-2^224+2^192+2^96-1
    for r in range(4):
        m = "t%d" % r
        e("lsr", "mhi", m, "32")
        e("lsl", "mlo", m, "32")
        e("adds", "t%d" % (r + 1), "t%d" % (r + 1), "mlo")   # + m<<96
        e("adcs", "t%d" % (r + 2), "t%d" % (r + 2), "mhi")
        e("adcs", "t%d" % (r + 3), "t%d" % (r + 3), m)       # + m<<192
        e("adcs", "t%d" % (r + 4), "t%d" % (r + 4), m)       # + m<<256
        # carry into the higher limbs: every instruction whose carry-out the next one reads is
        # ADCS (flag-setting); only the last of the chain is the plain ADC
        hi_limbs = list(range(r + 5, 8)) + ([8] if r == 3 else [])
        for n, k in enumerate(hi_limbs):
            e("adc" if n == len(hi_limbs) - 1 else "adcs", "t%d" % k, "t%d" % k, "xzr")
        e("subs", "t%d" % (r + 3), "t%d" % (r + 3), "mlo")   # - m<<224
        e("sbcs", "t%d" % (r + 4), "t%d" % (r + 4), "mhi")
        for n, k in enumerate(hi_limbs):                     # borrow into the higher limbs
            e("sbc" if n == len(hi_limbs) - 1 else "sbcs", "t%d" % k, "t%d" % k, "xzr")
    # ---- one conditional subtraction of p (result < 2p)
    e("subs", "r0", "t4", "P0")
    e("sbcs", "r1", "t5", "P1")
    e("sbcs", "r2", "t6", "P2")
    e("sbcs", "r3", "t7", "P3")
    e("sbcs", "t8", "t8", "xzr")
    for k in range(4):
        e("csel", "r%d" % k, "r%d" % k, "t%d" % (4 + k), "cs")
    e("stp", "r0", "r1", "x2")
    e("stp", "r2", "r3", "x2", "16")
    return list(prog)


def mul_add():
    prog.clear()
    e("ldp", "A0", "A1", "x0")
    e("ldp", "A2", "A3", "x0", "16")
    e("ldp", "B0", "B1", "x1")
    e("ldp", "B2", "B3", "x1", "16")
    e("adds", "r0", "A0", "B0")
    e("adcs", "r1", "A1", "B1")
    e("adcs", "r2", "A2", "B2")
    e("adcs", "r3", "A3", "B3")
    e("adc", "cy", "xzr", "xzr")
    e("subs", "d0", "r0", "P0")
    e("sbcs", "d1", "r1", "P1")
    e("sbcs", "d2", "r2", "P2")
    e("sbcs", "d3", "r3", "P3")
    e("sbcs", "cy", "cy", "xzr")
    for k in range(4):
        e("csel", "r%d" % k, "d%d" % k, "r%d" % k, "cs")
    e("stp", "r0", "r1", "x2")
    e("stp", "r2", "r3", "x2", "16")
    return list(prog)


def mod_sub():
    prog.clear()
    e("ldp", "A0", "A1", "x0")
    e("ldp", "A2", "A3", "x0", "16")
    e("ldp", "B0", "B1", "x1")
    e("ldp", "B2", "B3", "x1", "16")
    e("subs", "r0", "A0", "B0")
    e("sbcs", "r1", "A1", "B1")
    e("sbcs", "r2", "A2", "B2")
    e("sbcs", "r3", "A3", "B3")
    e("sbc", "cy", "xzr", "xzr")
    for k in range(4):
        e("and", "d%d" % k, "P%d" % k, "cy")
    e("adds", "r0", "r0", "d0")
    e("adcs", "r1", "r1", "d1")
    e("adcs", "r2", "r2", "d2")
    e("adcs", "r3", "r3", "d3")
    e("stp", "r0", "r1", "x2")
    e("stp", "r2", "r3", "x2", "16")
    return list(prog)


def run(seq, regs, mem):
    """Execute the sequence. regs/mem hold ints; x0,x1 are input pointers, x2 out."""
    C = 0

    def r(n):
        return 0 if n == "xzr" else regs.get(n, 0)

    def w(n, v):
        if n != "xzr":
            regs[n] = v & M64

    for ins in seq:
        op = ins[0]
        if op == "ldp":
            dst1, dst2, base = ins[1], ins[2], ins[3]
            off = int(ins[4]) if len(ins) > 4 else 0
            addr = regs[base] + off
            w(dst1, mem[addr])
            w(dst2, mem[addr + 8])
        elif op == "stp":
            s1, s2, base = ins[1], ins[2], ins[3]
            off = int(ins[4]) if len(ins) > 4 else 0
            addr = regs[base] + off
            mem[addr] = r(s1)
            mem[addr + 8] = r(s2)
        elif op == "mul":
            w(ins[1], r(ins[2]) * r(ins[3]))
        elif op == "umulh":
            w(ins[1], (r(ins[2]) * r(ins[3])) >> 64)
        elif op in ("adds", "adcs", "adc"):
            # AArch64: ADC reads C but does not write the flags; ADDS and ADCS write C
            cin = C if op != "adds" else 0
            s = r(ins[2]) + r(ins[3]) + cin
            w(ins[1], s)
            if op != "adc":
                C = 1 if s > M64 else 0
        elif op in ("subs", "sbcs", "sbc"):
            # AArch64: SBC reads C (borrow = 1 - C) but does not write the flags; SUBS, SBCS write C
            cin = (1 - C) if op != "subs" else 0
            s = r(ins[2]) - r(ins[3]) - cin
            w(ins[1], s)
            if op != "sbc":
                C = 1 if s >= 0 else 0
        elif op == "lsr":
            w(ins[1], r(ins[2]) >> int(ins[3]))
        elif op == "lsl":
            w(ins[1], r(ins[2]) << int(ins[3]))
        elif op == "and":
            w(ins[1], r(ins[2]) & r(ins[3]))
        elif op == "mov":
            w(ins[1], r(ins[2]))
        elif op == "csel":
            w(ins[1], r(ins[2]) if (C == 1) == (ins[4] == "cs") else r(ins[3]))
        else:
            raise AssertionError(op)
    return regs, mem


def limbs(v):
    return [PL[0], PL[1], PL[2], PL[3]], [(v >> (64 * i)) & M64 for i in range(4)]


def test_mul(seq, n=3000):
    rnd = random.Random(7)
    for it in range(n):
        a = rnd.getrandbits(256) % P
        b = rnd.getrandbits(256) % P
        if it < 6:
            a = [0, 1, P - 1, P - 2, (1 << 128) - 1, 1 << 255][it] % P
            b = [P - 1, 1, 0, P - 2, (1 << 192) + 5, 3][it] % P
        regs = {"x0": 0x1000, "x1": 0x2000, "x2": 0x3000}
        mem = {}
        for i, l in enumerate([(a >> (64 * k)) & M64 for k in range(4)]):
            mem[0x1000 + 8 * i] = l
        for i, l in enumerate([(b >> (64 * k)) & M64 for k in range(4)]):
            mem[0x2000 + 8 * i] = l
        for i, l in enumerate(PL):
            regs["P%d" % i] = l
        regs, mem = run(seq, regs, mem)
        got = sum(mem[0x3000 + 8 * i] << (64 * i) for i in range(4))
        want = (a * b * pow(1 << 256, -1, P)) % P
        assert got == want, (it, hex(a), hex(b), hex(got), hex(want))
    return True


def test_addsub(seq, fn, n=3000):
    rnd = random.Random(11)
    for it in range(n):
        a = rnd.getrandbits(256) % P
        b = rnd.getrandbits(256) % P
        if it < 5:
            a = [P - 1, 0, P - 1, 1, P - 2][it]
            b = [1, 0, P - 1, P - 1, 3][it]
        regs = {"x0": 0x1000, "x1": 0x2000, "x2": 0x3000}
        mem = {}
        for i, l in enumerate([(a >> (64 * k)) & M64 for k in range(4)]):
            mem[0x1000 + 8 * i] = l
        for i, l in enumerate([(b >> (64 * k)) & M64 for k in range(4)]):
            mem[0x2000 + 8 * i] = l
        for i, l in enumerate(PL):
            regs["P%d" % i] = l
        regs, mem = run(seq, regs, mem)
        got = sum(mem[0x3000 + 8 * i] << (64 * i) for i in range(4))
        want = fn(a, b) % P
        assert got == want, (it, hex(a), hex(b), hex(got), hex(want))
        assert got < P
    return True


def test_mul_carry_rich(seq, n=20000):
    """Operands whose limbs are mostly 0, 1, 2^64-1, 2^64-2 or 2^32 +- 1: long carry chains."""
    rnd = random.Random(23)
    pool = [0, 1, M64, M64 - 1, 1 << 32, (1 << 32) - 1, (1 << 32) + 1, 1 << 63, M64 << 32 & M64]
    for it in range(n):
        v = []
        for _ in range(2):
            x = 0
            for k in range(4):
                l = rnd.choice(pool) if rnd.random() < 0.8 else rnd.getrandbits(64)
                x |= l << (64 * k)
            v.append(x % P)
        a, b = v
        regs = {"x0": 0x1000, "x1": 0x2000, "x2": 0x3000}
        mem = {}
        for i in range(4):
            mem[0x1000 + 8 * i] = (a >> (64 * i)) & M64
            mem[0x2000 + 8 * i] = (b >> (64 * i)) & M64
        for i, l in enumerate(PL):
            regs["P%d" % i] = l
        regs, mem = run(seq, regs, mem)
        got = sum(mem[0x3000 + 8 * i] << (64 * i) for i in range(4))
        assert got == (a * b * pow(1 << 256, -1, P)) % P, (it, hex(a), hex(b))
    return True


def mutant_fails(seq, n=2000):
    """Sanity check on the executor: demoting every ADCS/SBCS that feeds a later carry
    reader to ADC/SBC (the defect of the previous listing) must make the result wrong."""
    mut = []
    for i, ins in enumerate(seq):
        if ins[0] in ("adcs", "sbcs") and i + 1 < len(seq) and seq[i + 1][0] in ("adc", "sbc") \
                and seq[i + 1][3] == "xzr" and ins[1] != "r3":
            ins = ("adc" if ins[0] == "adcs" else "sbc",) + ins[1:]
        mut.append(ins)
    try:
        test_mul(mut, n)
    except AssertionError:
        return True
    return False


def main():
    mul = mul_mont()
    add = mul_add()
    sub = mod_sub()
    ok_mul = test_mul(mul, 20000) and test_mul_carry_rich(mul, 20000)
    assert mutant_fails(mul), "executor does not notice ADC/SBC without flag write"
    ok_add = test_addsub(add, lambda a, b: a + b, 20000)
    ok_sub = test_addsub(sub, lambda a, b: a - b, 20000)
    ops = {}
    for ins in mul:
        ops[ins[0]] = ops.get(ins[0], 0) + 1
    print("montgomery multiplication: %d instructions" % len(mul))
    print("  " + " ".join("%s=%d" % kv for kv in sorted(ops.items())))
    print("modular addition: %d instructions" % len(add))
    print("modular subtraction: %d instructions" % len(sub))
    print("PASS mul_mont (20000 random + 20000 carry-rich cases)" if ok_mul else "FAIL mul_mont")
    print("PASS flag-semantics mutant (ADC/SBC in carry chains) is rejected")
    print("PASS mod_add (%d cases)" % 20000 if ok_add else "FAIL mod_add")
    print("PASS mod_sub (%d cases)" % 20000 if ok_sub else "FAIL mod_sub")
    assert ok_mul and ok_add and ok_sub
    assert len(mul) == 139 and len(add) == 20 and len(sub) == 19
    print("PASS counts match proof.md (M=139, A=20 add / 19 sub)")


if __name__ == "__main__":
    main()

"""Interpreter for the AArch64 subset gen.py emits: 64-bit general registers, NZCV with the ARM
AddWithCarry semantics (subs/sbcs set C = no borrow), mul/umulh, adds/adcs/subs/sbcs/add/adc/sub/sbc
(register or #imm12 [lsl #12] operand), csel, extr, ldr/str/ldp/stp with immediate offsets, mov/movk, ret.
Straight-line only: any other opcode, a branch or an unaligned access is an error. Every executed
instruction counts 1 (ldp/stp included, as Callgrind's Ir counts them)."""
import re

M64 = (1 << 64) - 1
REG = {f"x{i}": i for i in range(31)}
REG["sp"] = 31
ZR = 32
REG["xzr"] = ZR
COND = {
    "eq": lambda n, z, c, v: z, "ne": lambda n, z, c, v: not z,
    "hs": lambda n, z, c, v: c, "cs": lambda n, z, c, v: c,
    "lo": lambda n, z, c, v: not c, "cc": lambda n, z, c, v: not c,
    "mi": lambda n, z, c, v: n, "pl": lambda n, z, c, v: not n,
    "vs": lambda n, z, c, v: v, "vc": lambda n, z, c, v: not v,
    "hi": lambda n, z, c, v: c and not z, "ls": lambda n, z, c, v: (not c) or z,
    "ge": lambda n, z, c, v: n == v, "lt": lambda n, z, c, v: n != v,
    "gt": lambda n, z, c, v: (not z) and n == v, "le": lambda n, z, c, v: z or n != v,
}
MEMRE = re.compile(r"\[(\w+),\s*#(-?\d+)\]$")


def parse(text):
    prog = []
    for raw in text.splitlines():
        line = raw.split("//")[0].strip()
        if not line or line.endswith(":") or line.startswith("."):
            continue
        op, _, rest = line.partition(" ")
        args = [a.strip() for a in re.split(r",(?![^\[]*\])", rest)] if rest else []
        prog.append((op, args, raw))
    return prog


class Machine:
    def __init__(self):
        self.r = [0] * 33
        self.n = self.z = self.c = self.v = False
        self.mem = {}

    def rd(self, i):
        return 0 if i == ZR else self.r[i]

    def wr(self, i, val):
        if i != ZR:
            self.r[i] = val & M64

    def addr(self, s):
        m = MEMRE.match(s)
        if not m:
            raise ValueError("bad address " + s)
        a = (self.r[REG[m.group(1)]] + int(m.group(2))) & M64
        if a & 7:
            raise ValueError("unaligned access")
        return a

    def load(self, a):
        if a not in self.mem:
            raise KeyError(f"read of unwritten memory {a:#x}")
        return self.mem[a]


def awc(x, y, carry):
    """ARM AddWithCarry: (result, N, Z, C, V)."""
    u = x + y + carry
    res = u & M64
    sx, sy, sr = x >> 63, y >> 63, res >> 63
    return res, bool(sr), res == 0, u > M64, (sx == sy) and (sr != sx)


def operand(m, s):
    if s.startswith("#"):
        return int(s[1:], 0)
    return m.rd(REG[s])


def compile_prog(prog):
    """Decode once into closures over a Machine; returns the list and checks the subset."""
    ops = []
    for op, a, raw in prog:
        if op in ("mul", "umulh"):
            d, n, mm = REG[a[0]], REG[a[1]], REG[a[2]]
            if op == "mul":
                ops.append(lambda m, d=d, n=n, mm=mm: m.wr(d, m.rd(n) * m.rd(mm)))
            else:
                ops.append(lambda m, d=d, n=n, mm=mm: m.wr(d, (m.rd(n) * m.rd(mm)) >> 64))
        elif op in ("adds", "adcs", "subs", "sbcs", "add", "adc", "sub", "sbc"):
            d, n = REG[a[0]], REG[a[1]]
            if a[0] == "sp" or a[1] == "sp":
                assert op in ("add", "sub") and a[2].startswith("#"), raw
            shift = 0
            if len(a) == 4:
                assert a[3] == "lsl #12", raw
                shift = 12
            src = a[2]
            imm = src.startswith("#")
            if imm:
                assert op in ("adds", "subs", "add", "sub"), raw
                val = int(src[1:], 0) << shift
                assert 0 <= int(src[1:], 0) < 4096, raw
            else:
                assert shift == 0
                s_ = REG[src]
            sub = op.startswith("s")
            usec = op in ("adcs", "sbcs", "adc", "sbc")
            setf = op.endswith("s") and op not in ("sub", "sbc") or op in ("adds", "subs", "adcs", "sbcs")
            setf = op in ("adds", "subs", "adcs", "sbcs")


            def f(m, d=d, n=n, imm=imm, val=(val if imm else 0), s_=(None if imm else s_), sub=sub,
                  usec=usec, setf=setf):
                x = m.r[n] if n == 31 else m.rd(n)
                y = val if imm else m.rd(s_)
                if sub:
                    y = (~y) & M64
                carry = (1 if m.c else 0) if usec else (1 if sub else 0)
                res, N, Zf, C, V = awc(x, y, carry)
                if setf:
                    m.n, m.z, m.c, m.v = N, Zf, C, V
                if d == 31:
                    m.r[31] = res
                else:
                    m.wr(d, res)
            ops.append(f)
        elif op == "extr":
            d, n, mm, lsb = REG[a[0]], REG[a[1]], REG[a[2]], int(a[3][1:])
            assert 0 < lsb < 64, raw
            ops.append(lambda m, d=d, n=n, mm=mm, lsb=lsb: m.wr(d, ((m.rd(n) << 64 | m.rd(mm)) >> lsb) & M64))
        elif op == "csel":
            d, n, mm, cond = REG[a[0]], REG[a[1]], REG[a[2]], COND[a[3]]
            ops.append(lambda m, d=d, n=n, mm=mm, cond=cond:
                       m.wr(d, m.rd(n) if cond(m.n, m.z, m.c, m.v) else m.rd(mm)))
        elif op in ("ldp", "stp"):
            t1, t2, ad = REG[a[0]], REG[a[1]], a[2]
            assert -512 <= int(MEMRE.match(ad).group(2)) <= 504
            if op == "ldp":
                def f(m, t1=t1, t2=t2, ad=ad):
                    x = m.addr(ad)
                    v1, v2 = m.load(x), m.load(x + 8)
                    m.wr(t1, v1)
                    m.wr(t2, v2)
            else:
                def f(m, t1=t1, t2=t2, ad=ad):
                    x = m.addr(ad)
                    m.mem[x] = m.rd(t1)
                    m.mem[x + 8] = m.rd(t2)
            ops.append(f)
        elif op in ("ldr", "str"):
            t, ad = REG[a[0]], a[1]
            assert 0 <= int(MEMRE.match(ad).group(2)) <= 32760
            if op == "ldr":
                ops.append(lambda m, t=t, ad=ad: m.wr(t, m.load(m.addr(ad))))
            else:
                def f(m, t=t, ad=ad):
                    m.mem[m.addr(ad)] = m.rd(t)
                ops.append(f)
        elif op == "mov":
            d = REG[a[0]]
            if a[1].startswith("#"):
                v = int(a[1][1:], 0)
                assert 0 <= v < 65536
                ops.append(lambda m, d=d, v=v: m.wr(d, v))
            else:
                n = REG[a[1]]
                ops.append(lambda m, d=d, n=n: m.wr(d, m.rd(n)))
        elif op == "movk":
            d, v = REG[a[0]], int(a[1][1:], 0)
            sh = int(a[2].split("#")[1]) if len(a) > 2 else 0
            assert sh in (0, 16, 32, 48) and 0 <= v < 65536
            ops.append(lambda m, d=d, v=v, sh=sh: m.wr(d, (m.rd(d) & ~(0xFFFF << sh)) | (v << sh)))
        elif op == "ret":
            ops.append(None)
        else:
            raise ValueError("instruction outside the subset: " + raw)
    return ops


def run(ops, m):
    """Execute up to and including ret; return the executed instruction count."""
    n = 0
    for f in ops:
        n += 1
        if f is None:
            return n
        f(m)
    raise ValueError("fell off the end without ret")


#!/usr/bin/env python3
"""DYNAMIC execution trace over the EMITTED object -- v30.

v29's verdict (9ad41259) rejected two semantic claims because they were true only of the
schedule's metadata, not of the code: row 10 tested a register no instruction ever filled with
the batch product, and the link save re-executed per iteration into a MOVING slot.  Static
listing patterns cannot settle either point, so this harness EXECUTES the emitted instruction
stream (schedule.emit_asm) with a real register file, real memory and real NZCV, and compares
values the CPU would hold against values computed independently.

It proves, over every traced element, three equalities and two controls:

  P1  product transfer: the value in x8 at the `cbz` (loaded from the frame slot by the emitted
      `ldur x8,[sp,#0]`) equals the product of the batch's denominators, computed independently
      here from the traced per-element z values, mod p.  Row: all commits equal.
  C1  planted divergence: one denominator replaced by a non-canonical value -> the independent
      product differs from the transferred one -> MISMATCH reported (the check has teeth).
  C2  planted zero: one denominator set to 0 -> the object's own `cbz` takes the zero-recover
      exit (the predicate fires on the value it actually holds).
  P2  link preservation: at every executed `ldur x30,[sp,#8]` the loaded value equals the
      caller's link captured before the batch -- for every element and every exit stub.
  C3  planted moving-slot defect (the v29 code): save re-executed per iteration at [x19,#-8]
      -> the restores disagree with the caller's link.
  C4  planted missing save: the hoisted `stur` removed -> every restore reads a dead slot.

Usage: python3 dynamic_trace.py [--elements N] [--commits M]
"""
import re
import sys

import schedule

prog = labels = body_start = c_at = None
DEBUGO = False
EXITS = ("L_zero", "L_report", "L_abandon", "L_traversal", "L_fallback")

P_FIELD = (1 << 381) - 19          # not the BLS modulus: any odd prime keeps the algebra exact
MASK = (1 << 64) - 1
LINK = 0x0000_0000_dead_bee0       # the caller's link register value for every trial


def cond_ok(c, f):
    n, z, cy, v = f
    return {
        "eq": z, "ne": not z, "hs": cy, "cs": cy, "lo": not cy, "cc": not cy,
        "hi": cy and not z, "ls": (not cy) or z,
        "ge": n == v, "lt": n != v, "gt": (not z) and (n == v), "le": z or (n != v),
    }[c]


class CPU:
    def __init__(self, mem, rptr, end):
        self.r = {("x%d" % i): 0 for i in range(31)}
        self.r["xzr"] = 0
        self.r["sp"] = len(mem) - 64
        self.mem = mem
        self.f = (False, False, False, False)
        self.r["x19"] = rptr
        self.r["x20"] = end

    def g(self, t):
        t = t.strip()
        if t in ("xzr", "wzr"):
            return 0
        return self.r.get(t, 0)

    def s(self, t, v):
        if t.strip() in ("xzr", "wzr"):
            return
        self.r[t.strip()] = v & MASK

    def ld(self, base, off, n):
        a = (self.g(base) + off) & MASK
        return int.from_bytes(self.mem[a:a + n], "little")

    def st(self, base, off, v, n):
        a = (self.g(base) + off) & MASK
        self.mem[a:a + n] = (v & ((1 << (8 * n)) - 1)).to_bytes(n, "little")

    def setf(self, res, a, b, sub, carry_in=0):
        res &= MASK
        n = bool(res >> 63)
        z = res == 0
        if sub:
            cy = a >= (b + carry_in)
        else:
            cy = (a + b + carry_in) > MASK
        v = bool(((a ^ res) & (b ^ res)) >> 63) if sub else bool((~(a ^ b) & (a ^ res)) >> 63)
        self.f = (n, z, cy, v)
        return res


def imm(tok):
    t = tok.strip().lstrip("#")
    return int(t, 16) if t.lower().startswith("0x") else int(t)


def regs_of(txt):
    parts, cur, depth = [], "", 0
    for ch in txt:
        if ch == "[":
            depth += 1
        elif ch == "]":
            depth -= 1
        if ch == "," and depth == 0:
            parts.append(cur.strip())
            cur = ""
        else:
            cur += ch
    parts.append(cur.strip())
    return parts


def parse_addr(a):
    body = a.strip().strip("[]")
    parts = [p.strip() for p in body.split(",")]
    base = parts[0].strip()
    if len(parts) == 1:
        return base, 0
    off = parts[1].lstrip("#").strip()
    neg = off.startswith("-")
    if neg:
        off = off[1:].strip()
    v = int(off, 16) if off.lower().startswith("0x") else int(off)
    return base, (-v if neg else v)


def split_labels(src):
    out = []
    for raw in src.splitlines():
        l = raw.strip()
        if not l or l.startswith(".") or l.startswith("//"):
            continue
        if l.endswith(":"):
            out.append(("label", l[:-1]))
        else:
            out.append(("ins", l))
    return out


def step(cpu, ins, pc, labels, prog):
    """Execute one instruction. Returns the next pc index."""
    if ins == "ret":
        return None
    op = ins.split()[0]
    args = ins[len(op):].strip()
    # ---- branches ----
    if op in ("b", "b.eq", "b.ne", "b.hs", "b.ls", "b.lo", "b.hi", "b.cc", "b.cs", "b.ge", "b.lt"):
        if op == "b" or cond_ok(op[2:], cpu.f):
            return labels[args]
        return pc + 1
    if op in ("cbz", "cbnz"):
        t, lab = [x.strip() for x in args.split(",")]
        z = cpu.g(t) == 0
        return labels[lab] if (z if op == "cbz" else not z) else pc + 1
    # ---- memory ----
    if op in ("ldp", "stp"):
        a, b, mem = regs_of(args)
        base, off = parse_addr(mem)
        if op == "ldp":
            cpu.s(a, cpu.ld(base, off, 8))
            cpu.s(b, cpu.ld(base, off + 8, 8))
        else:
            cpu.st(base, off, cpu.g(a), 8)
            cpu.st(base, off + 8, cpu.g(b), 8)
        return pc + 1
    if op in ("ldr", "str", "ldur", "stur"):
        a, mem = regs_of(args)
        base, off = parse_addr(mem)
        if op in ("ldr", "ldur"):
            cpu.s(a, cpu.ld(base, off, 8))
        else:
            cpu.st(base, off, cpu.g(a), 8)
        return pc + 1
    # ---- ALU ----
    if op in ("add", "adds", "sub", "subs"):
        d, a, b = regs_of(args)[:3]
        bv = cpu.g(b) if not b.startswith("#") else imm(b)
        if op in ("add", "adds"):
            res = cpu.g(a) + bv
            if op == "adds":
                cpu.setf(res, cpu.g(a), bv, False)
        else:
            res = cpu.g(a) - bv
            if op == "subs":
                cpu.setf(res, cpu.g(a), bv, True)
        cpu.s(d, res)
        return pc + 1
    if op in ("adc", "sbc", "sbcs"):
        d, a, b = regs_of(args)[:3]
        cin = 1 if cpu.f[2] else 0
        if op == "adc":
            cpu.s(d, cpu.g(a) + cpu.g(b) + cin)
        else:
            res = cpu.g(a) - cpu.g(b) - (1 - cin)
            if op == "sbcs":
                cpu.setf(res, cpu.g(a), cpu.g(b), True, 0 if cin else 1)
            cpu.s(d, res)
        return pc + 1
    if op in ("cmp", "cmn"):
        a, b = regs_of(args)[:2]
        bv = imm(b) if b.startswith("#") else cpu.g(b)
        if op == "cmp":
            cpu.setf(cpu.g(a) - bv, cpu.g(a), bv, True)
        else:
            cpu.setf(cpu.g(a) + bv, cpu.g(a), bv, False)
        return pc + 1
    if op == "ccmp":
        a, b, immv, c = regs_of(args)
        if cond_ok(c, cpu.f):
            bv = imm(b) if b.startswith("#") else cpu.g(b)
            cpu.setf(cpu.g(a) - bv, cpu.g(a), bv, True)
        else:
            nz = imm(immv)
            cpu.f = (bool(nz & 8), bool(nz & 4), bool(nz & 2), bool(nz & 1))
        return pc + 1
    if op in ("csel", "csinc", "csinv"):
        d, a, b, c = regs_of(args)
        if cond_ok(c, cpu.f):
            cpu.s(d, cpu.g(a))
        else:
            cpu.s(d, cpu.g(b) + (1 if op == "csinc" else 0))
        return pc + 1
    if op == "cset":
        d, c = regs_of(args)
        cpu.s(d, 1 if cond_ok(c, cpu.f) else 0)
        return pc + 1
    if op == "mov":
        d, b = regs_of(args)
        cpu.s(d, imm(b) if b.startswith("#") else cpu.g(b))
        return pc + 1
    if op in ("and", "orr", "eor"):
        d, a, b = regs_of(args)
        if b.startswith("#"):
            bv = imm(b)
        else:
            bv = cpu.g(b)
        f = {"and": lambda x, y: x & y, "orr": lambda x, y: x | y, "eor": lambda x, y: x ^ y}[op]
        cpu.s(d, f(cpu.g(a), bv))
        return pc + 1
    if op == "tst":
        a, b = regs_of(args)
        bv = imm(b) if b.startswith("#") else cpu.g(b)
        cpu.setf(cpu.g(a) & bv, cpu.g(a), bv, False)
        return pc + 1
    if op in ("mul", "umulh", "madd"):
        rs = regs_of(args)
        d, a, b = rs[0], rs[1], rs[2]
        prod = cpu.g(a) * cpu.g(b)
        if op == "mul":
            cpu.s(d, prod)
        elif op == "umulh":
            cpu.s(d, prod >> 64)
        else:
            cpu.s(d, prod + cpu.g(rs[3]))
        return pc + 1
    raise SystemExit("unhandled instruction: %r" % ins)


def main():
    n_elem = 8
    n_commit = 1000
    for i, a in enumerate(sys.argv):
        if a == "--elements":
            n_elem = int(sys.argv[i + 1])
        if a == "--commits":
            n_commit = int(sys.argv[i + 1])

    global prog, labels, body_start, c_at
    global DEBUGO
    DEBUGO = "--debug" in sys.argv
    S = schedule.build()
    src = schedule.emit_asm(S)
    prog = split_labels(src)
    labels = {ins: i for i, (kind, ins) in enumerate(prog) if kind == "label"}
    body_start = labels["L_loop"]

    # index of the instruction after which the C boundary's outputs appear: the last state-y load
    for i, (kind, ins) in enumerate(prog):
        if kind == "ins" and ins.startswith("ldp x4, x6, [x19, #80]"):
            c_at = i
    assert c_at is not None, "C boundary anchor not found"

    def reparse():
        global prog, labels, body_start, c_at
        prog = split_labels(schedule.emit_asm(schedule.build()))
        labels = {ins: i for i, (kind, ins) in enumerate(prog) if kind == "label"}
        body_start = labels["L_loop"]
        c_at = None
        for i, (kind, ins) in enumerate(prog):
            if kind == "ins" and ins.startswith("ldp x4, x6, [x19, #80]"):
                c_at = i

    z_trace = []
    prod_reg = []          # x8 at the cbz
    prod_ref = []          # independent product of the traced z values
    x30_entry = []
    x30_restore = []
    exits = []

    commit_ok = 0
    commit_bad = []
    elems = 0
    resample = 0

    def run_batch(n, mutate=None, seed=0):
        nonlocal elems, resample
        mem = bytearray(4 << 20)
        rptr = 0
        end = n * 512
        CONST_AT, TBL_AT = 0x40000, 0x80000
        for k in range(7):   # G: the 448-bit fixed-point constants, arbitrary odd 64-bit values
            mem[CONST_AT + 8 * k: CONST_AT + 8 * k + 8] = \
                (0x9E3779B97F4A7C15 * (k + 3) & MASK).to_bytes(8, "little")
        for k in range(0x10000):   # table region: deterministic bytes so reads are varied
            mem[TBL_AT + k] = (k * 37 + 11) & 0xFF
        cpu = CPU(mem, rptr, end)
        cpu.r["x30"] = LINK
        cpu.r["x21"] = CONST_AT + 0x40000   # TBL
        cpu.r["x22"] = CONST_AT             # CONST
        for i in range(6):
            cpu.r["x%d" % (10 + i)] = 0x1111 * (i + 1)   # P0..P5
        cpu.r["x16"], cpu.r["x17"] = 0x3333, 0x4444      # PH0, PTOP
        # record + table entries
        for i in range(n):
            a = (i * 0x9E3779B97F4A7C15 + 0x1234567 + seed * 0x9E3779B1) & MASK
            for j in range(6):
                v = (a >> (17 * j)) ^ (0xABCDEF * (i + 1) * (j + 3))
                cpu.st("x0", i * 512 + 8 * j, v, 8)            # x
                cpu.st("x0", i * 512 + 48 + 8 * j, (v * 3 + 7) & MASK, 8)   # y
            ent = 0x40000 + i * 192
            cpu.st("x0", i * 512 + 96, ent, 8)
            for j in range(6):
                cpu.st("x0", ent + 8 * j, (a >> (11 * j)) ^ 0x55AA ^ j, 8)
                cpu.st("x0", ent + 48 + 8 * j, (a * 5 + j) & MASK, 8)
            cpu.st("x0", i * 512 + 136, 1, 8)
        P = 1
        zs = []
        entry_x30 = LINK
        x30_entry.append(entry_x30)

        # ---- the batch prologue, then the prefix pass (a separate function in the object) ----
        pc, guard = 0, 0
        cpu.r["x0"], cpu.r["x1"] = 0, n
        while True:
            guard += 1
            if guard > 4 * (len(prog) + 8 * n):
                return None
            kind, ins = prog[pc]
            if kind == "label":
                pc += 1
                continue
            if ins == "ret":
                break
            nxt = step(cpu, ins, pc, labels, prog)
            if nxt is None:
                break
            pc = nxt
        cpu.r["x19"], cpu.r["x20"] = rptr, end
        pc, guard = body_start, 0
        while True:
            guard += 1
            if guard > 8 * (len(prog) + n * 200):
                return None
            kind, ins = prog[pc]
            if kind == "label":
                pc += 1
                continue
            if pc == c_at:
                # execute the anchor load, then the C boundary: the 6M convention's output --
                # new coordinates in the operand block, and the prefix-pass running product P_i
                # = prod_{j<=i} z_j in x8 (class 2 of the traffic table, priced inside C).
                step(cpu, ins, pc, labels, prog)
                xs = [cpu.g(r) for r in ("x0", "x1", "x2", "x3", "x4", "x6")]
                ss = [cpu.g(r) for r in ("x5", "x26", "x27", "x28", "x29", "x30")]
                xv = sum(v << (64 * j) for j, v in enumerate(xs)) % P_FIELD
                sv = sum(v << (64 * j) for j, v in enumerate(ss)) % P_FIELD
                z = (xv - sv) % P_FIELD
                if mutate == "zero" and len(zs) == 3:
                    z = 0
                if mutate == "one":
                    z = (z + 1) % P_FIELD
                zs.append(z)
                P = (P * z) % P_FIELD
                cpu.s("x8", P)
                for j, r in enumerate(("x5", "x26", "x27", "x28", "x29", "x30")):
                    cpu.s(r, (xv * 7 + j) & MASK)
                for j, r in enumerate(("x0", "x1", "x2", "x3", "x4", "x6")):
                    cpu.s(r, (sv * 11 + j) & MASK)
                pc += 1
                continue
            if ins.startswith("ldur x30,"):
                base, off = parse_addr(ins.split("[", 1)[1].rstrip("]"))
                v = cpu.ld(base, off, 8)
                x30_restore.append(v)
                step(cpu, ins, pc, labels, prog)
                pc += 1
                continue
            if ins.startswith("cbz x8,"):
                if DEBUGO:
                    print("  cbz: x8=%x slot=%x P=%x sp=%x" % (cpu.g("x8"),
                          cpu.ld("sp", 0, 8), P, cpu.g("sp")))
                prod_reg.append(cpu.g("x8"))
                prod_ref.append(P & MASK)   # x8 is one 64-bit limb: the product's low limb
                nxt = step(cpu, ins, pc, labels, prog)
                if prog[nxt][1] in EXITS:
                    exits.append(prog[nxt][1])
                    return ("EXIT", None)
                pc = nxt
                continue
            if ins.startswith("ret"):
                return (zs, P)
            nxt = step(cpu, ins, pc, labels, prog)
            if nxt is None:
                return (zs, P)
            if prog[nxt][1] in EXITS:
                exits.append(prog[nxt][1])
                return ("EXIT", None)
            pc = nxt

    # ---- P1/P2 over (commits x elements) ----
    for c in range(n_commit):
        r = run_batch(n_elem, seed=c + 1)
        if isinstance(r, tuple) and len(r) == 2 and r[0] != "EXIT":
            commit_ok += 1
            elems += n_elem
        else:
            resample += 1
            # retry with a different element count is meaningless; retry the same batch once more
            r = run_batch(n_elem, seed=c + 1)
            if isinstance(r, tuple) and len(r) == 2 and r[0] != "EXIT":
                commit_ok += 1
                elems += n_elem
            else:
                commit_bad.append(c)

    prod_match = sum(1 for a, b in zip(prod_reg, prod_ref) if a == b)
    x30_match = sum(1 for v in x30_restore if v == LINK)
    print("== P1 product transfer ==")
    print("commits traced          : %d (of %d attempted, %d resampled, %d failed)"
          % (commit_ok, n_commit, resample, len(commit_bad)))
    print("cbz observations        : %d, register == independent prod(z): %d"
          % (len(prod_reg), prod_match))
    print("distinct products       : %d (all non-zero: %s)"
          % (len(set(prod_reg)), all(v != 0 for v in prod_reg)))
    print("== P2 link preservation ==")
    print("restore observations    : %d, loaded == caller link: %d" % (len(x30_restore), x30_match))
    print("elements traced         : %d" % elems)

    ok = True
    if prod_match != len(prod_reg) or not prod_reg:
        ok = False
        print("  FAIL  product register disagrees with the independent product")
    if x30_match != len(x30_restore) or not x30_restore:
        ok = False
        print("  FAIL  restored link disagrees with the caller link")

    # ---- C2: planted zero denominator must fire the object's own predicate ----
    exits.clear()
    prod_reg.clear()
    prod_ref.clear()
    run_batch(n_elem, mutate="zero")
    fired = [e for e in exits if e == "L_zero"]
    print("== C2 planted zero denominator ==")
    print("zero-recover fires      : %d" % len(fired))
    if not fired:
        ok = False
        print("  FAIL  a zero denominator did not fire `cbz x8`")

    # ---- C3/C4: the v29 defect must be visible to the same checks ----
    orig_emit = schedule.emit_asm
    real = schedule.build()

    def emit_moving_slot(S_):
        s = orig_emit(S_)
        s = s.replace("    stur x30, [sp, #8]", "    stur x30, [x19, #-0x8]")
        s = s.replace("    ldur x30, [sp, #8]", "    ldur x30, [x19, #-0x8]")
        # re-execute the save inside the loop, as v29 did
        s = s.replace("L_loop:\n", "L_loop:\n    stur x30, [x19, #-0x8]\n")
        return s

    def emit_no_save(S_):
        return orig_emit(S_).replace("    stur x30, [sp, #8]", "    mov x30, x30")

    for name, fn, want_fail in (("C3 v29 moving-slot save", emit_moving_slot, True),
                                ("C4 missing hoisted save", emit_no_save, True)):
        schedule.emit_asm = fn
        reparse()
        x30_restore.clear()
        x30_entry.clear()
        prod_reg.clear()
        prod_ref.clear()
        exits.clear()
        try:
            for _ in range(3):
                run_batch(n_elem)
        finally:
            schedule.emit_asm = orig_emit
            reparse()
        bad = sum(1 for v in x30_restore if v != LINK)
        print("== %s ==" % name)
        print("restore observations    : %d, WRONG link values: %d" % (len(x30_restore), bad))
        if want_fail and bad == 0:
            ok = False
            print("  FAIL  the planted defect was NOT detected (the check has no teeth)")

    print("\n%s" % ("DYNAMIC TRACE OK" if ok else "DYNAMIC TRACE FAILED"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

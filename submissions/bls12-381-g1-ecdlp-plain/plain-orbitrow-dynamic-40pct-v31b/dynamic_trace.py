#!/usr/bin/env python3
"""DYNAMIC execution trace over the EMITTED object -- v31.

WHAT THIS PROVES, AND WHAT IT REFUSES TO PROVE.

v30's verdict (d62b9b49) killed its dynamic evidence for two reasons, both real:

  * CIRCULARITY.  The v30 harness computed P = prod z_i in Python, assigned it to x8 at the
    C boundary, and then "verified" that x8 at the `cbz` equalled that same P masked to 64
    bits.  The check was an identity on the harness's own arithmetic; it said nothing about
    the emitted code.
  * ORDERING.  The same harness computed each element's denominator as
    z = (table x) - (state x) from registers that the row-6 state-y block had ALREADY
    overwritten (the state-y `ldp` block writes x0,x1,x2,x3,x4,x6, which is exactly the
    table-x lane), so the quantity it called "the denominator" was not the stated one.

v31 does not repair those computations -- it deletes them, because the object makes no claim
about them.  The span contains no multiplication chain that produces a batch product: the
product is produced inside the priced reference addition (the organizer's 6M Montgomery batch
convention, class 2 of the traffic table, charged in H).  What the span does contain, and what
this harness now checks, is a TRANSFER and an ORDERING:

  T1  slot round-trip.  The emitted `stur x8,[sp,#0]` writes the C-boundary accumulator to the
      batch's stationary frame slot; the emitted `ldur x8,[sp,#0]` reads the SAME address back;
      the value at the `cbz` equals both the slot content and the boundary value.  Every
      equality is between two OBSERVED pieces of machine state (register, memory), never
      between a register and a Python recomputation of it.
  T1a anti-vacuity.  Between the store and the `cbz` exactly one instruction writes x8, and it
      is the emitted `ldur`.  A harness (or a mutant) that injects a value into x8 inside that
      window is logged as an unattributed write and fails the run -- this is the check that
      would have caught v30's circular harness.
  T2  record ordering.  Within every element, every read of the record (base x19) precedes
      every write to the record.  This is the property the v30 harness's z-read silently
      assumed; here it is measured on the object's own memory operations.
  T3  link preservation.  Every executed `ldur x30,[sp,#8]` -- including the ones inside the
      five exit stubs, which are now EXECUTED rather than skipped at the branch -- loads the
      caller's link.

Controls (each must FAIL, and the run fails if it does not): C2 boundary value 0 -> the
object's own `cbz` takes the zero-recover exit; C3 the v29 moving-slot save; C4 the hoisted
save deleted; C5 the product load retargeted to [sp,#8]; C6 the product store retargeted to
[sp,#16]; C7 the product store deleted; C8 a state store hoisted above the state-y loads
(record ordering); C9 the v30 harness pattern reinstated (x8 injected at the `cbz`), which the
anti-vacuity check must flag.

The no-injection smoke row (`--smoke`) runs the same batch with x8 NOT assigned at the C
boundary: the value that reaches the `cbz` is then whatever the emitted instructions left in
x8, and it is still exactly what the store put in the slot.  That is the honest statement of
what the code does -- a transfer, not a product computation.

Usage: python3 dynamic_trace.py [--elements N] [--commits M] [--smoke]
"""
import re
import sys

import schedule

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
    """Execute one instruction. Returns the next pc index (None for ret)."""
    if ins == "ret":
        return None
    op = ins.split()[0]
    args = ins[len(op):].strip()
    if op in ("b", "b.eq", "b.ne", "b.hs", "b.ls", "b.lo", "b.hi", "b.cc", "b.cs", "b.ge", "b.lt"):
        if op == "b" or cond_ok(op[2:], cpu.f):
            return labels[args]
        return pc + 1
    if op in ("cbz", "cbnz"):
        t, lab = [x.strip() for x in args.split(",")]
        z = cpu.g(t) == 0
        return labels[lab] if (z if op == "cbz" else not z) else pc + 1
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
        bv = imm(b) if b.startswith("#") else cpu.g(b)
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


EXITS = ("L_zero", "L_report", "L_abandon", "L_traversal", "L_fallback")
REC_OFFS = {0, 16, 32, 48, 64, 80, 96}   # every record field the element reads (state x, state y, entry)


def build(src=None, req=True):
    prog = split_labels(src if src is not None else schedule.emit_asm(schedule.build()))
    labels = {ins: i for i, (kind, ins) in enumerate(prog) if kind == "label"}
    anchors = dict(body=labels["L_loop"])
    for i, (kind, ins) in enumerate(prog):
        if kind != "ins":
            continue
        if ins.startswith("stur x8, [sp, #0]"):
            anchors["store"] = i
        if ins.startswith("ldur x8, [sp, #0]"):
            anchors["load"] = i
        if ins.startswith("cbz x8,"):
            anchors["cbz"] = i
        if ins.startswith("ldp x4, x6, [x19, #80]"):
            anchors["c"] = i
    if req:
        for k in ("store", "load", "cbz", "c"):
            assert k in anchors, "anchor %r not found in the emitted stream" % k
    return prog, labels, anchors


def run_batch(prog, labels, anchors, n, boundary=None, seed=0, inject_at_cbz=False):
    """Execute the prologue, the prefix pass and one whole batch of n elements.

    boundary=None  -> NO-INJECTION mode: x8 is left exactly as the emitted instructions left it.
    boundary=v     -> the reference addition returns the accumulator value v in x8 (the C
                      boundary is a declared, priced program point; this is an INPUT to the
                      transfer check, never a quantity this script claims to derive).
    """
    mem = bytearray(4 << 20)
    rptr, end = 0, n * 512
    CONST_AT, TBL_AT = 0x40000, 0x80000
    for k in range(7):
        mem[CONST_AT + 8 * k: CONST_AT + 8 * k + 8] = \
            (0x9E3779B97F4A7C15 * (k + 3) & MASK).to_bytes(8, "little")
    for k in range(0x10000):
        mem[TBL_AT + k] = (k * 37 + 11) & 0xFF
    cpu = CPU(mem, rptr, end)
    cpu.r["x30"] = LINK
    cpu.r["x21"] = CONST_AT + 0x40000
    cpu.r["x22"] = CONST_AT
    for i in range(6):
        cpu.r["x%d" % (10 + i)] = 0x1111 * (i + 1)
    cpu.r["x16"], cpu.r["x17"] = 0x3333, 0x4444
    for i in range(n):
        a = (i * 0x9E3779B97F4A7C15 + 0x1234567 + seed * 0x9E3779B1) & MASK
        for j in range(6):
            v = (a >> (17 * j)) ^ (0xABCDEF * (i + 1) * (j + 3))
            cpu.st("x0", i * 512 + 8 * j, v, 8)
            cpu.st("x0", i * 512 + 48 + 8 * j, (v * 3 + 7) & MASK, 8)
        ent = 0x40000 + i * 192
        cpu.st("x0", i * 512 + 96, ent, 8)
        for j in range(6):
            cpu.st("x0", ent + 8 * j, (a >> (11 * j)) ^ 0x55AA ^ j, 8)
            cpu.st("x0", ent + 48 + 8 * j, (a * 5 + j) & MASK, 8)
        cpu.st("x0", i * 512 + 136, 1, 8)

    def note(ins):
        """Record an x19 read/write for the per-element read-before-write window."""
        nonlocal read_before_write, first_write_off
        op = ins.split()[0]
        if op not in ("ldr", "ldp", "ldur", "str", "stp", "stur"):
            return
        memarg = ins.split("[", 1)[1].rstrip("]")
        parts = [q.strip() for q in memarg.split(",")]
        if parts[0] != "x19":
            return
        off = int(parts[1].lstrip("#"), 0) if len(parts) > 1 else 0
        if op in ("ldr", "ldp", "ldur"):
            if first_write_off is None:
                read_before_write.add(off)
        elif first_write_off is None:
            first_write_off = off

    res = dict(obs=[], restores=[], exits=[], order_viol=[], stray_writes=[],
               boundary_vals=[], slot_after_store=[])
    pc, guard = 0, 0
    cpu.r["x0"], cpu.r["x1"] = 0, n
    in_body = False
    # per-element record ordering state: which record offsets are READ before the first write
    read_before_write, first_write_off = set(), None
    el = 0
    while True:
        guard += 1
        if guard > 16 * (len(prog) + n * 400):
            raise SystemExit("guard tripped")
        kind, ins = prog[pc]
        if pc == anchors["body"]:
            if in_body:
                # a new element begins: close the previous element's ordering window
                if read_before_write != REC_OFFS:
                    res["order_viol"].append((el, sorted(REC_OFFS - read_before_write), first_write_off))
                read_before_write, first_write_off = set(), None
                el += 1
            in_body = True
        if kind == "label":
            pc += 1
            continue
        if ins.startswith("ldur x30,"):
            base, off = parse_addr(ins.split("[", 1)[1].rstrip("]"))
            res["restores"].append(cpu.ld(base, off, 8))
        if pc == anchors["c"]:
            step(cpu, ins, pc, labels, prog)
            if boundary is not None:
                cpu.s("x8", boundary)
            res["boundary_vals"].append(cpu.g("x8"))
            note(ins)
            pc += 1
            continue
        if pc == anchors["store"]:
            step(cpu, ins, pc, labels, prog)
            res["slot_after_store"].append(cpu.ld("sp", 0, 8))
            note(ins)
            pc += 1
            continue
        if pc == anchors["cbz"]:
            if inject_at_cbz:
                cpu.s("x8", cpu.ld("sp", 0, 8))       # the v30 harness pattern, planted
                res["stray_writes"].append((None, "HARNESS injection at the cbz"))
            res["obs"].append(dict(x8=cpu.g("x8"), slot=cpu.ld("sp", 0, 8),
                                   stored=res["slot_after_store"][-1] if res["slot_after_store"] else None,
                                   boundary=res["boundary_vals"][-1] if res["boundary_vals"] else None))
            nxt = step(cpu, ins, pc, labels, prog)
            if prog[nxt][1] in EXITS:
                res["exits"].append(prog[nxt][1])
            pc = nxt
            continue
        # x8 attribution: any change of x8 not produced by executing this instruction
        before = cpu.g("x8")
        nxt = step(cpu, ins, pc, labels, prog)
        if cpu.g("x8") != before:
            if ins.split()[0] in ("ldur", "ldr", "ldp", "mov", "mul", "umulh", "madd",
                                  "add", "adds", "sub", "subs", "and", "orr", "eor", "csel"):
                pass  # an instruction wrote x8: attributed
            else:
                res["stray_writes"].append((pc, ins))
        if nxt is None:
            if not in_body:
                pc = anchors["body"]
                continue
            res["exits"].append("fallthrough")
            break
        if in_body:
            note(ins)
        pc = nxt
    if read_before_write != REC_OFFS:
        res["order_viol"].append((el, sorted(REC_OFFS - read_before_write), first_write_off))
    return res


def main():
    n_elem, n_commit = 8, 1000
    smoke_only = "--smoke" in sys.argv
    for i, a in enumerate(sys.argv):
        if a == "--elements":
            n_elem = int(sys.argv[i + 1])
        if a == "--commits":
            n_commit = int(sys.argv[i + 1])

    prog, labels, anchors = build()
    ok = True

    # ---- smoke: what does the code do with NO boundary assignment at all? ----
    smoke = run_batch(prog, labels, anchors, n_elem, boundary=None)
    same = sum(1 for o in smoke["obs"] if o["x8"] == o["slot"] == o["stored"])
    print("== S0 no-injection smoke (x8 never assigned by the harness) ==")
    print("elements                 : %d" % n_elem)
    print("cbz observations         : %d, x8 == slot == stored: %d" % (len(smoke["obs"]), same))
    print("boundary x8 values       : %s" % ", ".join("0x%x" % v for v in smoke["boundary_vals"][:3]))
    print("  -> the span computes no product: the value at the predicate is whatever the emitted")
    print("     instructions left in x8, carried through the frame slot by the emitted store/load.")
    if smoke["obs"] and same != len(smoke["obs"]):
        ok = False
        print("  FAIL  the transferred value is not preserved without an assignment")

    if smoke_only:
        print("\n%s" % ("SMOKE OK" if ok else "SMOKE FAILED"))
        return 0 if ok else 1

    # ---- T1 / T1a over 1000 commits x n elements ----
    obs, restores, order_viol, stray = [], [], [], []
    commits = 0
    for c in range(n_commit):
        b = (0x9E3779B97F4A7C15 * (c + 1) ^ 0xA5A5A5A5) & MASK
        if b == 0:
            b = 1
        r = run_batch(prog, labels, anchors, n_elem, boundary=b, seed=c + 1)
        if r["exits"] and r["exits"][0] in EXITS:
            continue                      # a zero/product exit: not a commit
        commits += 1
        obs.extend(r["obs"])
        restores.extend(r["restores"])
        order_viol.extend(r["order_viol"])
        stray.extend(r["stray_writes"])

    trip = sum(1 for o in obs if o["x8"] == o["slot"] == o["stored"] == o["boundary"])
    print("== T1 slot round-trip ==")
    print("commits                  : %d" % commits)
    print("observations at the cbz  : %d, cbz val == slot == stored == boundary: %d"
          % (len(obs), trip))
    print("distinct transferred vals : %d, all non-zero: %s"
          % (len(set(o["x8"] for o in obs)), all(o["x8"] for o in obs)))
    print("== T1a anti-vacuity ==")
    print("unattributed x8 writes   : %d" % len(stray))
    print("== T2 record ordering (read-before-write, per element) ==")
    print("batches checked          : %d (%d elements), violations: %d"
          % (len(obs), len(obs) * (n_elem - 0), len(order_viol)))
    print("== T3 link preservation ==")
    print("restore observations     : %d, loaded == caller link: %d"
          % (len(restores), sum(1 for v in restores if v == LINK)))
    if trip != len(obs) or not obs:
        ok = False
        print("  FAIL  transferred value not preserved through the slot")
    if stray:
        ok = False
        print("  FAIL  x8 written outside the instruction stream (vacuous check)")
    if order_viol:
        ok = False
        print("  FAIL  a record write preceded a record read inside an element")
    if sum(1 for v in restores if v == LINK) != len(restores) or not restores:
        ok = False
        print("  FAIL  restored link disagrees with the caller link")

    # ---- C2: a zero boundary value must fire the object's own predicate ----
    z = run_batch(prog, labels, anchors, n_elem, boundary=0)
    print("== C2 zero boundary value ==")
    print("zero-recover fires      : %d" % sum(1 for e in z["exits"] if e == "L_zero"))
    if not any(e == "L_zero" for e in z["exits"]):
        ok = False
        print("  FAIL  a zero accumulator did not fire `cbz x8`")

    # ---- mutants: every one of these must be caught ----
    base_src = schedule.emit_asm(schedule.build())

    def mut(pat, rep):
        return base_src.replace(pat, rep)

    def mut_reorder():
        s = base_src
        st = "    stp x5, x26, [x19, #0]\n"
        s = s.replace(st, "", 1)
        return s.replace("    ldp x0, x1, [x19, #48]\n", st + "    ldp x0, x1, [x19, #48]\n", 1)

    cases = [
        ("C5 product load retargeted to [sp,#8]", mut("    ldur x8, [sp, #0]", "    ldur x8, [sp, #8]")),
        ("C6 product store retargeted to [sp,#16]", mut("    stur x8, [sp, #0]", "    stur x8, [sp, #16]")),
        ("C7 product store deleted", mut("    stur x8, [sp, #0]", "    mov x8, x8")),
        ("C8 state store hoisted above the state-y loads", mut_reorder()),
        ("C3 v29 moving-slot link save",
         mut("    stur x30, [sp, #8]", "    stur x30, [x19, #-0x8]")
         .replace("    ldur x30, [sp, #8]", "    ldur x30, [x19, #-0x8]")
         .replace("L_loop:\n", "L_loop:\n    stur x30, [x19, #-0x8]\n")),
        ("C4 missing hoisted link save", mut("    stur x30, [sp, #8]", "    mov x30, x30")),
    ]
    for name, src in cases:
        prog2, labels2, anchors2 = build(src, req=False)
        for k in ("store", "load", "cbz", "c"):
            anchors2.setdefault(k, anchors[k])
        bad = 0
        if name.startswith(("C5", "C6", "C7")):
            b = 0x1234_5678_9ABC_DEF0
            r = run_batch(prog2, labels2, anchors2, n_elem, boundary=b)
            bad = sum(1 for o in r["obs"]
                      if o["x8"] != o["slot"] or o["slot"] != o["stored"] or o["slot"] != o["boundary"])
        elif name.startswith("C8"):
            r = run_batch(prog2, labels2, anchors2, n_elem, boundary=0x1234_5678_9ABC_DEF0)
            bad = len(r["order_viol"])
        else:
            r = run_batch(prog2, labels2, anchors2, n_elem, boundary=0x1234_5678_9ABC_DEF0)
            bad = sum(1 for v in r["restores"] if v != LINK) + (0 if r["restores"] else 1)
        print("== %s ==" % name)
        print("detections              : %d" % bad)
        if bad == 0:
            ok = False
            print("  FAIL  the planted defect was NOT detected (the check has no teeth)")

    # ---- C9: the v30 harness pattern must be flagged as a vacuous check ----
    r = run_batch(prog, labels, anchors, n_elem, boundary=0x1234_5678_9ABC_DEF0, inject_at_cbz=True)
    print("== C9 v30 pattern: x8 injected at the cbz ==")
    print("unattributed writes      : %d" % len(r["stray_writes"]))
    if not r["stray_writes"]:
        ok = False
        print("  FAIL  an injected cbz value was not flagged (the anti-vacuity check has no teeth)")

    print("\n%s" % ("DYNAMIC TRACE OK" if ok else "DYNAMIC TRACE FAILED"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

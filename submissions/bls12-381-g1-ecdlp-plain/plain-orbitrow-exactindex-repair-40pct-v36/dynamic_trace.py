#!/usr/bin/env python3
"""
DYNAMIC execution trace over the EMITTED listing -- v33.

(v33: see proof.md section 9 for the change log this file implements.)
"""
import sys

import schedule

MASK = (1 << 64) - 1
LINK = 0x0000_0000_dead_bee0       # the caller's link register value for every trial
FIELD = (0, 16, 32, 48, 64, 80, 96, 104, 120, 136)
FIELD_NAME = {0: "x limbs 0,1", 16: "x limbs 2,3", 32: "x limbs 4,5", 48: "y limbs 0,1",
              64: "y limbs 2,3", 80: "y limbs 4,5", 96: "entry pointer",
              104: "K_cp limbs 0,2", 120: "K_cp limb 3,valid", 136: "length (2 limbs)"}
SHIFTS = ("lsl", "lsr", "asr", "ror")
EXITS = ("L_zero", "L_report", "L_abandon", "L_traversal", "L_fallback")
SX_LANE = ("x5", "x26", "x27", "x28", "x29", "x30")
SY_LANE = ("x0", "x1", "x2", "x3", "x4", "x6")
TY_LANE = ("x7", "x8", "x9", "x23", "x24", "x25")


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

    def rd(self, a, n):
        a &= MASK
        return int.from_bytes(self.mem[a:a + n], "little")

    def ld(self, base, off, n):
        return self.rd((self.g(base) + off) & MASK, n)

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


def shifter(rs, ins):
    """The shifter on the last source operand, or None.  A form this parser does not model is a
    hard error: v32 silently ignored it, and that silence is what the verdict killed."""
    if len(rs) < 4:
        return None
    toks = rs[3].split()
    if toks and toks[0] in SHIFTS:
        if len(toks) != 2:
            raise SystemExit("unhandled shifter form %r in %r" % (rs[3], ins))
        return (toks[0], toks[1])
    return None


def apply_shift(v, sh):
    if sh is None:
        return v
    op, n = sh[0], imm(sh[1]) & 63
    if op == "lsl":
        return (v << n) & MASK
    if op == "lsr":
        return v >> n
    if op == "asr":
        return (((v ^ (1 << 64)) - (1 << 64)) >> n) & MASK
    return ((v >> n) | (v << ((64 - n) & 63))) & MASK


def split_labels(src):
    out = []
    for raw in src.splitlines():
        l = raw.strip()
        if not l or l.startswith(".") or l.startswith("//"):
            continue
        out.append(("label", l[:-1]) if l.endswith(":") else ("ins", l))
    return out


def step(cpu, ins, pc, labels, prog, shift_blind=False):
    """Execute one instruction.  Returns (next_pc, event) with event in (None, "ret")."""
    op = ins.split()[0]
    args = ins[len(op):].strip()
    if op == "ret":
        return None, "ret"
    if op in ("b", "b.eq", "b.ne", "b.hs", "b.ls", "b.lo", "b.hi", "b.cc", "b.cs", "b.ge", "b.lt"):
        if op == "b" or cond_ok(op[2:], cpu.f):
            return labels[args], None
        return pc + 1, None
    if op in ("cbz", "cbnz"):
        t, lab = [x.strip() for x in args.split(",")]
        z = cpu.g(t) == 0
        return (labels[lab] if (z if op == "cbz" else not z) else pc + 1), None
    if op in ("ldp", "stp"):
        a, b, mem = regs_of(args)
        base, off = parse_addr(mem)
        if op == "ldp":
            cpu.s(a, cpu.ld(base, off, 8))
            cpu.s(b, cpu.ld(base, off + 8, 8))
        else:
            cpu.st(base, off, cpu.g(a), 8)
            cpu.st(base, off + 8, cpu.g(b), 8)
        return pc + 1, None
    if op in ("ldr", "str", "ldur", "stur"):
        a, mem = regs_of(args)
        base, off = parse_addr(mem)
        if op in ("ldr", "ldur"):
            cpu.s(a, cpu.ld(base, off, 8))
        else:
            cpu.st(base, off, cpu.g(a), 8)
        return pc + 1, None
    if op in ("add", "adds", "sub", "subs"):
        rs = regs_of(args)
        d, a, b = rs[0], rs[1], rs[2]
        sh = None if shift_blind else shifter(rs, ins)
        bv = apply_shift(imm(b) if b.startswith("#") else cpu.g(b), sh)
        av = cpu.g(a)
        if op in ("add", "adds"):
            res = av + bv
            if op == "adds":
                cpu.setf(res, av, bv, False)
        else:
            res = av - bv
            if op == "subs":
                cpu.setf(res, av, bv, True)
        cpu.s(d, res)
        return pc + 1, None
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
        return pc + 1, None
    if op in ("cmp", "cmn"):
        a, b = regs_of(args)[:2]
        bv = imm(b) if b.startswith("#") else cpu.g(b)
        if op == "cmp":
            cpu.setf(cpu.g(a) - bv, cpu.g(a), bv, True)
        else:
            cpu.setf(cpu.g(a) + bv, cpu.g(a), bv, False)
        return pc + 1, None
    if op == "ccmp":
        a, b, immv, c = regs_of(args)
        if cond_ok(c, cpu.f):
            bv = imm(b) if b.startswith("#") else cpu.g(b)
            cpu.setf(cpu.g(a) - bv, cpu.g(a), bv, True)
        else:
            nz = imm(immv)
            cpu.f = (bool(nz & 8), bool(nz & 4), bool(nz & 2), bool(nz & 1))
        return pc + 1, None
    if op in ("csel", "csinc", "csinv"):
        d, a, b, c = regs_of(args)
        if cond_ok(c, cpu.f):
            cpu.s(d, cpu.g(a))
        else:
            cpu.s(d, cpu.g(b) + (1 if op == "csinc" else 0))
        return pc + 1, None
    if op == "cset":
        d, c = regs_of(args)
        cpu.s(d, 1 if cond_ok(c, cpu.f) else 0)
        return pc + 1, None
    if op == "mov":
        d, b = regs_of(args)
        cpu.s(d, imm(b) if b.startswith("#") else cpu.g(b))
        return pc + 1, None
    if op in ("and", "orr", "eor"):
        d, a, b = regs_of(args)
        bv = imm(b) if b.startswith("#") else cpu.g(b)
        f = {"and": lambda x, y: x & y, "orr": lambda x, y: x | y, "eor": lambda x, y: x ^ y}[op]
        cpu.s(d, f(cpu.g(a), bv))
        return pc + 1, None
    if op == "tst":
        a, b = regs_of(args)
        bv = imm(b) if b.startswith("#") else cpu.g(b)
        cpu.setf(cpu.g(a) & bv, cpu.g(a), bv, False)
        return pc + 1, None
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
        return pc + 1, None
    raise SystemExit("unhandled instruction: %r" % ins)


def build(src=None, req=True):
    prog = split_labels(src if src is not None else schedule.emit_asm(schedule.build()))
    labels = {ins: i for i, (kind, ins) in enumerate(prog) if kind == "label"}
    anchors = dict(body=labels["L_loop"], entry=labels["_batch_entry"])
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


def run_batch(prog, labels, anchors, n, boundary=None, seed=0, inject_at_cbz=False,
              shift_blind=False):
    """Execute the WHOLE function: `_batch_entry` -> prologue -> prefix pass -> leaf loop -> ret.

    boundary=None  -> NO-INJECTION mode: x8 is left exactly as the emitted instructions left it.
    boundary=v     -> the reference addition returns the accumulator value v in x8 (the C boundary
                      is a declared, priced program point; this is an INPUT to the transfer check,
                      never a quantity this script claims to derive).
    """
    mem = bytearray(4 << 20)
    rptr, end = 0, n * 512
    CONST_AT, ENT_AT = 0x40000, 0x100000
    for k in range(7):
        mem[CONST_AT + 8 * k: CONST_AT + 8 * k + 8] = \
            (0x9E3779B97F4A7C15 * (k + 3) & MASK).to_bytes(8, "little")
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
        ent = ENT_AT + i * 192
        cpu.st("x0", i * 512 + 96, ent, 8)
        for j in range(6):
            cpu.st("x0", ent + 8 * j, (a >> (11 * j)) ^ 0x55AA ^ j, 8)
            cpu.st("x0", ent + 48 + 8 * j, (a * 5 + j) & MASK, 8)
        cpu.st("x0", i * 512 + 104, (a ^ 0x104) & MASK, 8)
        cpu.st("x0", i * 512 + 120, (a ^ 0x120) & MASK, 8)
        # half the elements carry a length that is 0 mod 64, so the 1-in-64 refresh
        # path (the only writer of the checkpoint fields 104/120) is EXERCISED
        cpu.st("x0", i * 512 + 136, 64 if i % 2 == 0 else 1, 8)
    snapshot = {}
    for i in range(n):
        ent = ENT_AT + i * 192
        snapshot[ent] = bytes(mem[ent:ent + 192])

    res = dict(obs=[], restores=[], exits=[], order_viol=[], field_rd={}, field_wr={},
               stray_writes=[], boundary_vals=[], slot_after_store=[], bad_returns=[],
               returns=[], opmap=[], nonfield_writes=[], entw=[])

    rd_at, wr_at, el, counter = {}, {}, [0], [0]

    def note(ins):
        op = ins.split()[0]
        if op not in ("ldr", "ldp", "ldur", "str", "stp", "stur"):
            return
        parts = [q.strip() for q in ins.split("[", 1)[1].rstrip("]").split(",")]
        if parts[0] != "x19":
            return
        off = int(parts[1].lstrip("#"), 0) if len(parts) > 1 else 0
        if op in ("ldr", "ldp", "ldur"):
            rd_at.setdefault(off, counter[0])
        else:
            wr_at.setdefault(off, counter[0])

    def close_element():
        bad = []
        for off in FIELD:
            if off in wr_at and off not in rd_at:
                bad.append((off, "write with no read at all"))
            elif off in wr_at and rd_at[off] > wr_at[off]:
                bad.append((off, "read after write"))
        for off in sorted(wr_at):
            if off not in FIELD:
                res["nonfield_writes"].append((el[0], off))
        for off in FIELD:
            if off in rd_at:
                res["field_rd"][off] = res["field_rd"].get(off, 0) + 1
            if off in wr_at:
                res["field_wr"][off] = res["field_wr"].get(off, 0) + 1
        if bad:
            res["order_viol"].append((el[0], bad))
        rd_at.clear()
        wr_at.clear()

    pc, guard = anchors["entry"], 0
    cpu.r["x0"], cpu.r["x1"] = 0, n
    in_body = False
    while True:
        guard += 1
        if guard > 16 * (len(prog) + n * 400):
            raise SystemExit("guard tripped: the control flow did not terminate")
        kind, ins = prog[pc]
        if pc == anchors["body"]:
            if in_body:
                close_element()
                el[0] += 1
            in_body = True
        if kind == "label":
            pc += 1
            continue
        counter[0] += 1
        if ins.startswith("ldur x30,"):
            base, off = parse_addr(ins.split("[", 1)[1].rstrip("]"))
            res["restores"].append(cpu.ld(base, off, 8))
        if pc == anchors["c"]:
            step(cpu, ins, pc, labels, prog, shift_blind)
            ent = cpu.ld("x19", 96, 8)
            snap = snapshot.get(ent, b"\0" * 192)
            res["opmap"].append(dict(
                SX=sum(1 for k, r in enumerate(SX_LANE) if cpu.g(r) == cpu.ld("x19", 8 * k, 8)),
                SY=sum(1 for k, r in enumerate(SY_LANE) if cpu.g(r) == cpu.ld("x19", 48 + 8 * k, 8)),
                TY=sum(1 for k, r in enumerate(TY_LANE) if cpu.g(r) == cpu.rd(ent + 48 + 8 * k, 8)),
                TX=sum(1 for k in range(6) if cpu.rd(ent + 8 * k, 8) ==
                       int.from_bytes(snap[8 * k:8 * k + 8], "little"))))
            if boundary is not None:
                cpu.s("x8", boundary)
            res["boundary_vals"].append(cpu.g("x8"))
            note(ins)
            pc += 1
            continue
        if pc == anchors["store"]:
            step(cpu, ins, pc, labels, prog, shift_blind)
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
            nxt, _ = step(cpu, ins, pc, labels, prog, shift_blind)
            if nxt is not None and prog[nxt][1] in EXITS:
                res["exits"].append(prog[nxt][1])
            pc = nxt
            continue
        before = cpu.g("x8")
        nxt, ev = step(cpu, ins, pc, labels, prog, shift_blind)
        if ins.startswith("str x6, [x19, #96]"):
            res["entw"].append(cpu.g("x6"))           # row 5's shifted-add result, as stored
        if cpu.g("x8") != before and ins.split()[0] not in (
                "ldur", "ldr", "ldp", "mov", "mul", "umulh", "madd",
                "add", "adds", "sub", "subs", "and", "orr", "eor", "csel"):
            res["stray_writes"].append((pc, ins))
        if in_body:
            note(ins)
        if ev == "ret":
            res["returns"].append(cpu.g("x30"))
            if cpu.g("x30") != LINK:
                res["bad_returns"].append((pc, cpu.g("x30")))
            if in_body:
                close_element()
            break
        pc = nxt
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

    obs, restores, order_viol, stray, opmap, nonfield = [], [], [], [], [], []
    fld_rd, fld_wr = {}, {}
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
        opmap.extend(r["opmap"])
        nonfield.extend(r["nonfield_writes"])
        for k, v in r["field_rd"].items():
            fld_rd[k] = fld_rd.get(k, 0) + v
        for k, v in r["field_wr"].items():
            fld_wr[k] = fld_wr.get(k, 0) + v
        if r["bad_returns"] or len(r["returns"]) != 1:
            ok = False
            print("  FAIL  the batch did not return exactly once to the caller link")

    trip = sum(1 for o in obs if o["x8"] == o["slot"] == o["stored"] == o["boundary"])
    print("== T1 slot round-trip ==")
    print("commits                  : %d" % commits)
    print("observations at the cbz  : %d, cbz val == slot == stored == boundary: %d"
          % (len(obs), trip))
    print("distinct transferred vals : %d, all non-zero: %s"
          % (len(set(o["x8"] for o in obs)), all(o["x8"] for o in obs)))
    print("== T1a anti-vacuity ==")
    print("unattributed x8 writes   : %d" % len(stray))
    print("== T2 per-field read-before-write, ALL %d record fields ==" % len(FIELD))
    print("elements checked         : %d, violations: %d" % (commits * n_elem, len(order_viol)))
    print("  off  field                 reads    writes  read-before-own-write")
    for off in FIELD:
        r_, w_ = fld_rd.get(off, 0), fld_wr.get(off, 0)
        print("  %3d  %-20s %6d  %7d   %s"
              % (off, FIELD_NAME[off], r_, w_,
                 "yes" if r_ and w_ else ("read-only" if r_ else "no read/write")))
    print("  non-field x19 writes   : %d" % len(nonfield))
    for off in (104, 120):
        if not fld_wr.get(off):
            ok = False
            print("  FAIL  field %d is never written, so its read-before-write is untested" % off)
    print("== T3 link preservation ==")
    print("restore observations     : %d, loaded == caller link: %d"
          % (len(restores), sum(1 for v in restores if v == LINK)))
    print("== T4 operand map at the C boundary (24 limbs of the two points) ==")
    full = [m for m in opmap if m["SX"] == 6 and m["SY"] == 6 and m["TY"] == 6 and m["TX"] == 6]
    print("boundaries observed      : %d, all 24 limbs mapped: %d" % (len(opmap), len(full)))
    print("  state x (6) in registers : %d/%d" % (sum(m["SX"] == 6 for m in opmap), len(opmap)))
    print("  state y (6) in registers : %d/%d" % (sum(m["SY"] == 6 for m in opmap), len(opmap)))
    print("  table y (6) in registers : %d/%d" % (sum(m["TY"] == 6 for m in opmap), len(opmap)))
    print("  table x (6) at [x19,#96]->0..40, entry block unmodified: %d/%d"
          % (sum(m["TX"] == 6 for m in opmap), len(opmap)))
    if trip != len(obs) or not obs:
        ok = False
        print("  FAIL  transferred value not preserved through the slot")
    if stray:
        ok = False
        print("  FAIL  x8 written outside the instruction stream (vacuous check)")
    if order_viol:
        ok = False
        print("  FAIL  a record write preceded a read of the same field")
    if nonfield:
        ok = False
        print("  FAIL  a write at an offset that is not a declared record field")
    if sum(1 for v in restores if v == LINK) != len(restores) or not restores:
        ok = False
        print("  FAIL  restored link disagrees with the caller link")
    if len(full) != len(opmap) or not opmap:
        ok = False
        print("  FAIL  the operand map is not the declared one at some boundary")

    z = run_batch(prog, labels, anchors, n_elem, boundary=0)
    print("== C2 zero boundary value ==")
    print("zero-recover fires      : %d" % sum(1 for e in z["exits"] if e == "L_zero"))
    if not any(e == "L_zero" for e in z["exits"]):
        ok = False
        print("  FAIL  a zero accumulator did not fire `cbz x8`")

    base_src = schedule.emit_asm(schedule.build())

    def mut(pat, rep):
        s = base_src.replace(pat, rep)
        assert s != base_src, "mutation %r did not apply" % pat
        return s

    def mut_reorder():
        s = base_src
        st = "    stp x5, x26, [x19, #0]\n"
        s = s.replace(st, "", 1)
        # v28's defect: the result store of offset 0 ahead of the state-x read of offset 0
        return s.replace("    ldp x5, x26, [x19, #0]\n",
                         st + "    ldp x5, x26, [x19, #0]\n", 1)

    def mut_field(off_store, off_before):
        """a record store at off_store hoisted above the load at off_before"""
        s = base_src
        st = "    stp x0, x2, [x19, #%d]\n" % off_store
        s = s.replace(st, "", 1)
        return s.replace("    ldp x7, x8, [x19, #%d]\n" % off_before,
                         st + "    ldp x7, x8, [x19, #%d]\n" % off_before, 1)

    cases = [
        ("C5 product load retargeted to [sp,#8]", mut("    ldur x8, [sp, #0]", "    ldur x8, [sp, #8]")),
        ("C6 product store retargeted to [sp,#16]", mut("    stur x8, [sp, #0]", "    stur x8, [sp, #16]")),
        ("C7 product store deleted", mut("    stur x8, [sp, #0]", "    mov x8, x8")),
        ("C8 result store @0 hoisted above the state-x read @0", mut_reorder()),
        ("C3 v29 moving-slot link save",
         mut("    stur x30, [sp, #8]", "    stur x30, [x19, #-0x8]")
         .replace("    ldur x30, [sp, #8]", "    ldur x30, [x19, #-0x8]")
         .replace("L_loop:\n", "L_loop:\n    stur x30, [x19, #-0x8]\n")),
        ("C4 missing hoisted link save", mut("    stur x30, [sp, #8]", "    mov x30, x30")),
        ("C13 checkpoint store @104 hoisted above the @104 load", mut_field(104, 104)),
        ("C14 table-y block retargeted one field down",
         mut("    ldp x7, x8, [x5, #48]", "    ldp x7, x8, [x5, #64]")),
    ]
    for name, src in cases:
        prog2, labels2, anchors2 = build(src, req=False)
        for k in ("store", "load", "cbz", "c"):
            anchors2.setdefault(k, anchors[k])
        b = 0x1234_5678_9ABC_DEF0
        r = run_batch(prog2, labels2, anchors2, n_elem, boundary=b)
        if name.startswith(("C5", "C6", "C7")):
            bad = sum(1 for o in r["obs"]
                      if o["x8"] != o["slot"] or o["slot"] != o["stored"] or o["slot"] != o["boundary"])
        elif name.startswith(("C8", "C13")):
            bad = len(r["order_viol"])
        elif name.startswith("C14"):
            bad = sum(1 for m in r["opmap"] if m["TY"] != 6)
        else:
            bad = sum(1 for v in r["restores"] if v != LINK) + (0 if r["restores"] else 1)
        print("== %s ==" % name)
        print("detections              : %d" % bad)
        if bad == 0:
            ok = False
            print("  FAIL  the planted defect was NOT detected (the check has no teeth)")

    r = run_batch(prog, labels, anchors, n_elem, boundary=0x1234_5678_9ABC_DEF0, inject_at_cbz=True)
    print("== C9 v30 pattern: x8 injected at the cbz ==")
    print("unattributed writes      : %d" % len(r["stray_writes"]))
    if not r["stray_writes"]:
        ok = False
        print("  FAIL  an injected cbz value was not flagged (the anti-vacuity check has no teeth)")

    # ---- C10 shifted-add faithfulness: fixed sees the mutation, shift-blind does not ----
    b = 0x1234_5678_9ABC_DEF0
    mut_src = mut("    add x6, x21, x6, lsl #5", "    add x6, x21, x6, lsl #4")
    prog3, labels3, anchors3 = build(mut_src, req=False)
    for k in ("store", "load", "cbz", "c"):
        anchors3.setdefault(k, anchors[k])
    r0 = run_batch(prog, labels, anchors, n_elem, boundary=b)
    r1 = run_batch(prog3, labels3, anchors3, n_elem, boundary=b)
    r0b = run_batch(prog, labels, anchors, n_elem, boundary=b, shift_blind=True)
    r1b = run_batch(prog3, labels3, anchors3, n_elem, boundary=b, shift_blind=True)
    d_fixed = sum(1 for x, y in zip(r0["entw"], r1["entw"]) if x != y)
    d_blind = sum(1 for x, y in zip(r0b["entw"], r1b["entw"]) if x != y)
    print("== C10 shifted-add faithfulness (`lsl #5` -> `lsl #4`, table address) ==")
    print("stored table addresses that change, fixed interpreter : %d/%d"
          % (d_fixed, len(r0["entw"])))
    print("same comparison, shift-blind shadow (v32 parser)      : %d" % d_blind)
    if d_fixed == 0 or d_fixed != len(r0["entw"]) or d_blind != 0:
        ok = False
        print("  FAIL  the fix does not distinguish shifted from unshifted operands")

    # ---- C11 return address corrupted at the leaf epilogue ----
    head, tail = base_src.split("L_fallback:", 1)
    old = "    ldur x30, [sp, #8]"
    i = head.rfind(old)
    src11 = head[:i] + "    ldur x30, [sp, #0]" + head[i + len(old):] + "L_fallback:" + tail
    prog4, labels4, anchors4 = build(src11, req=False)
    for k in ("store", "load", "cbz", "c"):
        anchors4.setdefault(k, anchors[k])
    r = run_batch(prog4, labels4, anchors4, n_elem, boundary=b)
    print("== C11 return address corrupted at the leaf epilogue ==")
    print("bad returns             : %d" % len(r["bad_returns"]))
    if not r["bad_returns"]:
        ok = False
        print("  FAIL  a corrupted return address was not detected")

    # ---- C12 a write to far field 136 with no read of that field ----
    src12 = base_src.replace("L_loop:\n", "L_loop:\n    stp x27, x28, [x19, #136]\n", 1)
    prog5, labels5, anchors5 = build(src12, req=False)
    for k in ("store", "load", "cbz", "c"):
        anchors5.setdefault(k, anchors[k])
    r = run_batch(prog5, labels5, anchors5, n_elem, boundary=b)
    fire = [v for v in r["order_viol"] if any(o == 136 for o, _ in v[1])]
    print("== C12 write to record field 136 before any read of it ==")
    print("field-136 violations    : %d" % len(fire))
    if not fire:
        ok = False
        print("  FAIL  an unread field write at 136 was not detected")

    print("\n%s" % ("DYNAMIC TRACE OK" if ok else "DYNAMIC TRACE FAILED"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

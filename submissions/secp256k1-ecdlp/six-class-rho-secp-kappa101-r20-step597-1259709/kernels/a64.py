"""The AArch64 interpreter and instruction counter of experiment.py (proof.md section 2.6).

_asm_counts counts each function's instructions and the instructions per iteration of its single-block loops;
a64_load and a64_run execute the shipped kernels/secp_loops.s, kernels/secp_tests.s and kernels/arrive.s.
"""
M64 = (1 << 64) - 1


def _asm_counts(text):
    """Per function: instructions, and for single-block loops the instructions per iteration.
    Directive lines (.text, .globl, ...) are skipped; labels are `name:` (functions), `.Lname:` or `N:`
    (local); a branch to `Nb` is a backward branch to the nearest preceding label N."""
    funcs = {}
    cur = None
    for line in text.split("\n"):
        if not line.strip():
            continue
        if not line[0].isspace():
            name = line.split(":")[0]
            if name.startswith(".L") or name.isdigit():
                funcs[cur]["items"].append(("label", name))
            else:
                cur = name
                funcs[cur] = {"items": []}
            continue
        if line.strip().startswith("."):
            continue
        funcs[cur]["items"].append(("ins", line.strip()))
    res = {}
    for name, f in funcs.items():
        items = f["items"]
        pos = {}
        idx = 0
        for kind, val in items:
            if kind == "label":
                pos.setdefault(val, []).append(idx)
            else:
                idx += 1
        total = idx
        loops = []
        idx = 0
        for kind, val in items:
            if kind == "ins":
                mn = val.split()[0]
                tgt = val.split()[-1]
                key = tgt[:-1] if (tgt[:-1].isdigit() and tgt.endswith("b")) else tgt
                if (mn in ("b", "cbz", "cbnz", "tbz", "tbnz") or mn.startswith("b.")) and key in pos and pos[key][0] <= idx:
                    start = max(p for p in pos[key] if p <= idx)
                    body = [v for k, v in items if k == "ins"][start:idx + 1]
                    inner = [b for b in body[:-1] if b.split()[0] in ("b", "bl", "cbz", "cbnz", "tbz", "tbnz", "ret") or b.split()[0].startswith("b.")]
                    loops.append(dict(label=tgt, iteration=len(body), inner_branches=len(inner)))
                idx += 1
        res[name] = dict(instructions=total, loops=loops,
                         loop_iteration_max=max([l["iteration"] for l in loops]) if loops else 0,
                         calls=sum(1 for k, v in items if k == "ins" and v.split()[0] == "bl"))
    return res


# ---- a small AArch64 interpreter: executes the shipped kernels/secp_loops.s and kernels/secp_tests.s
# (proof.md section 2.6); flags C and Z only, which is all the kernels read
def _a64_ops(s):
    ops, cur, depth = [], "", 0
    for ch in s:
        if ch == "[":
            depth += 1
        if ch == "]":
            depth -= 1
        if ch == "," and depth == 0:
            ops.append(cur.strip())
            cur = ""
        else:
            cur += ch
    if cur.strip():
        ops.append(cur.strip())
    return ops


def a64_load(text):
    prog, labels = [], {}
    for raw in text.split("\n"):
        line = raw.split("//")[0].rstrip()
        if not line.strip():
            continue
        if not line[0].isspace():
            labels.setdefault(line.rstrip(":"), []).append(len(prog))
            continue
        ins = line.strip()
        if ins.startswith("."):
            continue
        mn, _, rest = ins.partition(" ")
        prog.append((mn, _a64_ops(rest)))
    return prog, labels


def a64_run(prog, labels, entry, regs, mem, limit=10**7):
    """Run from label `entry` until RET; returns the number of executed instructions. regs: dict x0..x30, sp."""
    M = M64
    R = {("x%d" % i): 0 for i in range(31)}
    R["sp"] = 0
    R.update(regs)
    flags = {"C": 0, "Z": 0}
    pc = labels[entry][0]
    steps = 0

    def rd(o):
        if o == "xzr":
            return 0
        return R[o]

    def wr(o, v):
        if o != "xzr":
            R[o] = v & M

    def imm(o):
        return int(o.lstrip("#"), 0)

    def ea(op, post):
        # returns (address, new_base_value or None, base_register)
        inner = op.strip()
        pre = inner.endswith("!")
        inner = inner.rstrip("!")[1:-1]
        parts = [p.strip() for p in inner.split(",")]
        base = parts[0]
        off = imm(parts[1]) if len(parts) > 1 else 0
        if post is not None:
            return rd(base), (rd(base) + imm(post)) & M, base
        if pre:
            return (rd(base) + off) & M, (rd(base) + off) & M, base
        return (rd(base) + off) & M, None, base

    def cond(c):
        return {"eq": flags["Z"] == 1, "ne": flags["Z"] == 0, "hs": flags["C"] == 1, "cs": flags["C"] == 1,
                "lo": flags["C"] == 0, "cc": flags["C"] == 0}[c]

    def tgt(t):           # a named label, or N b / N f: the nearest label N before / after this instruction
        if t[:-1].isdigit():
            c = labels[t[:-1]]
            return [p for p in c if p < pc][-1] if t[-1] == "b" else [p for p in c if p >= pc][0]
        return labels[t][0]

    while True:
        mn, o = prog[pc]
        steps += 1
        assert steps < limit
        pc += 1
        if mn in ("mul", "umulh", "smulh"):
            a, b = rd(o[1]), rd(o[2])
            if mn == "smulh":
                a, b = a - (a >> 63 << 64), b - (b >> 63 << 64)
            wr(o[0], (a * b) & M if mn == "mul" else (a * b) >> 64)
        elif mn == "asr":
            a = rd(o[1])
            wr(o[0], (a - (a >> 63 << 64)) >> imm(o[2]))
        elif mn in ("adds", "adc", "adcs", "add") and not (mn == "add" and (o[2].startswith("#") or len(o) == 4)):
            a, b = rd(o[1]), rd(o[2])
            if len(o) == 4:       # adds xd, xn, xm, lsl #k (shifted register; the bits shifted out are dropped)
                sh = o[3].split()
                assert mn == "adds" and sh[0] == "lsl"
                b = (b << imm(sh[1])) & M
            cin = flags["C"] if mn in ("adc", "adcs") else 0
            t = a + b + cin
            wr(o[0], t)
            if mn in ("adds", "adcs"):
                flags["C"] = t >> 64
                flags["Z"] = 1 if (t & M) == 0 else 0
        elif mn == "add":     # add xd, xn, #imm | add xd, xn, xm, lsl #k
            if o[2].startswith("#"):
                wr(o[0], rd(o[1]) + imm(o[2]))
            else:
                wr(o[0], rd(o[1]) + (rd(o[2]) << imm(o[3].split()[1])))
        elif mn in ("sub", "subs", "sbc", "sbcs"):
            a = rd(o[1])
            b = imm(o[2]) if o[2].startswith("#") else rd(o[2])
            cin = flags["C"] if mn in ("sbc", "sbcs") else 1
            t = a + ((~b) & M) + cin
            wr(o[0], t)
            if mn in ("subs", "sbcs"):
                flags["C"] = t >> 64
                flags["Z"] = 1 if (t & M) == 0 else 0
        elif mn == "csel":
            wr(o[0], rd(o[1]) if cond(o[3]) else rd(o[2]))
        elif mn == "madd":
            wr(o[0], rd(o[1]) * rd(o[2]) + rd(o[3]))
        elif mn in ("cmp", "cmn", "ccmp", "ccmn"):     # flags C and Z; ccmp/ccmn set #nzcv when cond fails
            if mn[1] == "c" and not cond(o[3]):
                flags["C"], flags["Z"] = imm(o[2]) >> 1 & 1, imm(o[2]) >> 2 & 1
            else:
                a, b = rd(o[0]), imm(o[1]) if o[1].startswith("#") else rd(o[1])
                t = a + b if mn[-1] == "n" else a + ((~b) & M) + 1
                flags["C"], flags["Z"] = t >> 64, 1 if (t & M) == 0 else 0
        elif mn == "orr":
            wr(o[0], rd(o[1]) | imm(o[2]))
        elif mn == "bfi":
            mask = ((1 << imm(o[3])) - 1) << imm(o[2])
            wr(o[0], (rd(o[0]) & ~mask) | ((rd(o[1]) << imm(o[2])) & mask))
        elif mn == "extr":
            wr(o[0], ((rd(o[1]) << 64 | rd(o[2])) >> imm(o[3])))
        elif mn == "lsr":
            wr(o[0], rd(o[1]) >> imm(o[2]))
        elif mn == "lsl":
            wr(o[0], rd(o[1]) << imm(o[2]))
        elif mn == "mov":
            wr(o[0], imm(o[1]) if o[1].startswith("#") else rd(o[1]))
        elif mn == "movk":
            sh = imm(o[2].split()[1]) if len(o) > 2 else 0
            wr(o[0], (rd(o[0]) & ~(0xFFFF << sh)) | ((imm(o[1]) & 0xFFFF) << sh))
        elif mn in ("ldp", "stp"):
            post = o[3] if len(o) > 3 else None
            addr, newb, base = ea(o[2], post)
            if mn == "ldp":
                wr(o[0], mem[addr])
                wr(o[1], mem[addr + 8])
            else:
                mem[addr] = rd(o[0])
                mem[addr + 8] = rd(o[1])
            if newb is not None:
                R[base] = newb
        elif mn in ("ldr", "str"):
            post = o[2] if len(o) > 2 else None
            addr, newb, base = ea(o[1], post)
            if mn == "ldr":
                wr(o[0], mem[addr])
            else:
                mem[addr] = rd(o[0])
            if newb is not None:
                R[base] = newb
        elif mn in ("b", "b.eq", "b.ne", "b.lo", "b.hs"):
            if mn == "b" or cond(mn[2:]):
                pc = tgt(o[0])
        elif mn in ("cbz", "tbz", "tbnz"):
            bit = rd(o[0]) if mn == "cbz" else rd(o[0]) >> imm(o[1]) & 1
            if (bit == 0) == (mn != "tbnz"):
                pc = tgt(o[-1])
        elif mn == "ret":
            return steps
        else:
            raise ValueError("unsupported instruction: " + mn)

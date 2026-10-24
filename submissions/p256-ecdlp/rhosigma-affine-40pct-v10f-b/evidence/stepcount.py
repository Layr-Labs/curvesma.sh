#!/usr/bin/env python3
"""stepcount.py — machine-generate the instruction rows of the counted walk step.

Reuses the v10 counter's method verbatim (clang CFG, backward edges dropped,
longest path from the loop head) and then classifies every instruction on the
counted path, so the ledger's rows ARE the counted block: they sum to the
reported longest_iteration by construction, not by hand.

usage: stepcount.py step_fix.s out.json
"""
import json
import sys

FUNC = __import__("re").compile(r"^([A-Za-z_][\w.$]*):")
LABEL = __import__("re").compile(r"^(\.[\w.$]+):")
BODY = __import__("re").compile(r"^\t([a-z][\w.]*)\s*(.*)$")
BR = {"b", "br"}


def parse(path):
    funcs, order = {}, []
    cur, labels, insns = None, {}, []
    for line in open(path):
        if not line.strip() or line.lstrip().startswith(("#", "//")):
            continue
        m = FUNC.match(line)
        if m and not line.startswith("\t"):
            if cur:
                funcs[cur] = (insns, labels)
            cur, labels, insns = m.group(1), {}, []
            order.append(cur)
            continue
        m = LABEL.match(line)
        if m and cur:
            labels[m.group(1)] = len(insns)
            continue
        m = BODY.match(line)
        if m and cur:
            insns.append((m.group(1), m.group(2).split("//")[0].strip()))
    if cur:
        funcs[cur] = (insns, labels)
    return funcs, order


def longest(insns, labels):
    n = len(insns)
    succ = []
    for i, (op, arg) in enumerate(insns):
        s = []
        tgt = arg.split(",")[-1].strip() if arg else ""
        is_branch = op in BR or op.startswith("b.") or op in ("cbz", "cbnz", "tbz", "tbnz")
        if is_branch and tgt in labels and labels[tgt] > i:
            s.append(labels[tgt])
        if op in BR:
            pass
        elif i + 1 < n:
            s.append(i + 1)
        succ.append(s)
    memo, path = {}, {}

    def best(i):
        if i in memo:
            return memo[i]
        memo[i], nxt = 1, None
        for j in succ[i]:
            if j >= n:
                continue
            v = 1 + best(j)
            if v > memo[i]:
                memo[i], nxt = v, j
        path[i] = nxt
        return memo[i]

    sys.setrecursionlimit(100000)
    heads = [labels[l] for l in labels if labels[l] < n]
    if not heads:
        return 0, []
    top = max(heads, key=lambda h: best(h))
    out, i = [], top
    while i is not None:
        out.append(insns[i])
        i = path.get(i)
    return len(out), out


MULS = {"mul", "umulh", "madd", "msub", "mneg", "smulh"}


def klass(op):
    if op in MULS:
        return "product"
    if op in ("add", "adc", "adcs", "sub", "sbc", "sbcs", "adds", "subs", "neg"):
        return "carry_chain"
    if op in ("csel", "csinc", "csinv", "cset", "csetm", "ccmp", "ccmn"):
        return "masked_select"
    if op in ("lsr", "lsl", "asr", "and", "orr", "eor", "bic", "ubfx", "sbfx", "extr", "rev", "rbit"):
        return "mix_key"
    if op in ("ldr", "ldur", "ldp", "ldrb", "ldrh", "ldrsb"):
        return "load"
    if op in ("str", "stur", "stp", "strb", "strh"):
        return "store"
    if op in ("cinc", "cinv", "cneg"):
        return "carry_propagate"
    if op in ("mrs", "msr"):
        return "guard_read"
    if op in ("cmp", "cmn", "tst"):
        return "flagtest"
    if op.startswith("mov") or op in ("sxtw", "uxtw", "sxtb", "uxtb", "sxth", "uxth", "nop"):
        return "move_extend"
    if op in BR or op.startswith("b.") or op in ("cbz", "cbnz", "tbz", "tbnz", "ret"):
        return "control"
    return "other"


def main():
    src, outp = sys.argv[1], sys.argv[2]
    funcs, order = parse(src)
    res = {}
    for fn in order:
        insns, labels = funcs[fn]
        if not insns:
            continue
        ln, blk = longest(insns, labels)
        if not ln:
            continue
        rows = {}
        for op, _ in blk:
            rows[klass(op)] = rows.get(klass(op), 0) + 1
        n_umulh = sum(1 for op, _ in blk if op == "umulh")
        n_mul = sum(1 for op, _ in blk if op == "mul")
        res[fn] = {"longest_iteration": ln, "rows": rows,
                   "umulh": n_umulh, "mul": n_mul,
                   "sum_of_rows": sum(rows.values())}
    json.dump(res, open(outp, "w"), indent=1, sort_keys=True)
    for fn, r in res.items():
        print(fn, r["longest_iteration"], "rows", r["sum_of_rows"], r["rows"])


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""count_iter.py — longest iteration of every looping function, from clang .s.

Method (the same in kind as the counts this track has been priced with):
  1. compile the shipped C for aarch64 with the reference flags (see run.sh);
  2. for each function, build the instruction CFG, drop every edge that jumps
     backwards (a loop back edge ends one iteration);
  3. report, per function, max over loop heads of the longest path from the
     head to a vertex with no outgoing in-iteration edge.

usage: count_iter.py p256_walk_v9f.s [--json out.json]
"""
import json
import re
import sys

FUNC = re.compile(r"^([A-Za-z_][\w.$]*):")
LABEL = re.compile(r"^(\.L[\w.$]+):")
BODY = re.compile(r"^\t([a-z][\w.]*)\s*(.*)$")
BR = {"b", "br"}


def parse(path):
    funcs, order = {}, []
    cur, labels, insns = None, {}, []
    for line in open(path):
        if not line.strip() or line.lstrip().startswith(("#", "//")):
            continue
        fm = FUNC.match(line)
        if fm and not line.startswith("\t"):
            if cur:
                funcs[cur] = (insns, labels)
            cur, labels, insns = fm.group(1), {}, []
            order.append(cur)
            continue
        lm = LABEL.match(line)
        if lm and cur:
            labels[lm.group(1)] = len(insns)
            continue
        bm = BODY.match(line)
        if bm and cur:
            insns.append((bm.group(1), bm.group(2).split("//")[0].strip()))
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
        if op in BR:                      # unconditional: no fallthrough
            pass
        elif i + 1 < n:
            s.append(i + 1)
        succ.append(s)
    # DAG (all backward edges dropped): longest path from each loop head
    memo = {}

    def best(i):
        if i in memo:
            return memo[i]
        memo[i] = 1
        for j in succ[i]:
            if j >= n:
                continue
            memo[i] = max(memo[i], 1 + best(j))
        return memo[i]

    sys.setrecursionlimit(100000)
    heads = [labels[l] for l in labels if labels[l] < n]
    if not heads:
        return 0, None
    top = max(heads, key=lambda h: best(h))
    return best(top), top


def main():
    funcs, order = parse(sys.argv[1])
    out = {}
    for fn in order:
        insns, labels = funcs[fn]
        if not insns:
            continue
        ln, head = longest(insns, labels)
        out[fn] = {"in_file": len(insns), "longest_iteration": ln}
        if ln:
            print("%-14s in_file %5d   longest_iteration %5d" % (fn, len(insns), ln))
    if "--json" in sys.argv:
        p = sys.argv[sys.argv.index("--json") + 1]
        json.dump(out, open(p, "w"), indent=1, sort_keys=True)
        print("wrote", p)


if __name__ == "__main__":
    main()

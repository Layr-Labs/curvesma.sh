#!/usr/bin/env python3
"""Static instruction count of one shipped walk step, from objdump -d.

The step and every function it calls compile to straight-line code: no
conditional branch of any kind (data dependent or not), no loop back-edge, no
early return.  The only control transfers in the step's call graph are the
fixed `bl` call sites.  Therefore the executed instruction count of one call is
exactly

    sum over functions f of  count(f) * calls(f)

where count(f) is the number of instructions in f's disassembly (prologue,
body and `ret`, all executed once per call) and calls(f) is the number of times
f is entered per walk step, which the call graph fixes.

This script:
  1. parses `objdump -d --no-show-raw-insn <binary>`;
  2. finds the step's call graph transitively from the root symbol;
  3. asserts every reached function contains no conditional branch and no
     backwards unconditional branch (i.e. no loop), and that every call site
     is a direct `bl`;
  4. reports count(f), calls(f) and the exact per-step total.

Usage: static_count.py <binary> [root-symbol]
Exit code 0 iff the straight-line assertion holds.
"""
import re
import subprocess
import sys

COND = re.compile(r'^(b\.[a-z]+|cbz|cbnz|tbz|tbnz)\b')
INSN = re.compile(r'^\s*([0-9a-f]+):\s+(\S+)\s*(.*)$')
SYM = re.compile(r'^([0-9a-f]+) <(.+)>:$')


def disasm(path):
    out = subprocess.run(['objdump', '-d', '--no-show-raw-insn', path],
                         check=True, capture_output=True, text=True).stdout
    funcs = {}
    cur = None
    for line in out.splitlines():
        m = SYM.match(line)
        if m:
            cur = m.group(2)
            funcs[cur] = []
            continue
        m = INSN.match(line)
        if m and cur is not None:
            funcs[cur].append((int(m.group(1), 16), m.group(2), m.group(3)))
    # Drop alignment padding after the last `ret`; it is never executed.
    for name, insns in funcs.items():
        last = max((i for i, (_, mn, _) in enumerate(insns) if mn == 'ret'), default=None)
        if last is not None and last + 1 < len(insns):
            funcs[name] = insns[:last + 1]
    return funcs


def analyse(funcs, root):
    calls = {}          # caller -> {callee: number of bl sites}
    extern = set()
    for name, insns in funcs.items():
        for addr, mn, ops in insns:
            if mn == 'bl':
                m = re.match(r'[0-9a-f]+ <([^>+]+)', ops)
                if not m:
                    continue
                callee = m.group(1)
                calls.setdefault(name, {}).setdefault(callee, 0)
                calls[name][callee] += 1
                if callee not in funcs:
                    extern.add((name, callee))
    mult = {root: 1}
    changed = True
    while changed:
        changed = False
        for f, n in list(mult.items()):
            for g, c in calls.get(f, {}).items():
                if g not in funcs:
                    continue
                if mult.get(g, 0) < n * c:
                    mult[g] = n * c
                    changed = True
    problems = ['%s: call to symbol outside the binary: %s' % (a, b) for a, b in sorted(extern)]
    for name in mult:
        for addr, mn, ops in funcs[name]:
            if COND.match(mn):
                problems.append('%s: conditional branch %s at %x' % (name, mn, addr))
            if mn == 'b':
                tgt = re.match(r'([0-9a-f]+)', ops)
                if tgt and int(tgt.group(1), 16) < addr:
                    problems.append('%s: backwards branch at %x' % (name, addr))
            if mn in ('br', 'blr'):
                problems.append('%s: indirect branch %s at %x' % (name, mn, addr))
    return mult, calls, problems


def main():
    path = sys.argv[1]
    root = sys.argv[2] if len(sys.argv) > 2 else 'rho_step.constprop.0'
    funcs = disasm(path)
    if root not in funcs:
        cands = [n for n in funcs if n.startswith('rho_step')]
        if not cands:
            print('root symbol not found', file=sys.stderr)
            return 2
        root = cands[0]
    mult, calls, problems = analyse(funcs, root)
    total = 0
    print('root: %s' % root)
    print('%-28s %8s %8s %10s' % ('function', 'insns', 'calls', 'insns/step'))
    for f in sorted(mult, key=lambda x: -mult[x]):
        n = len(funcs[f])
        total += n * mult[f]
        print('%-28s %8d %8d %10d' % (f, n, mult[f], n * mult[f]))
        for g, c in sorted(calls.get(f, {}).items()):
            print('    %-24s -> %-28s x%d' % (f, g, c))
    print('%-28s %8s %8s %10d' % ('TOTAL', '', '', total))
    if problems:
        print('\nSTRAIGHT-LINE ASSERTION FAILED:')
        for p in problems:
            print('  ' + p)
        return 1
    print('\nOK: every function in the call graph is straight-line '
          '(no conditional branch, no loop, no early return).')
    return 0


if __name__ == '__main__':
    sys.exit(main())

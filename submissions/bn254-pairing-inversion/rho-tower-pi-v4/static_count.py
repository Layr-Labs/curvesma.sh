#!/usr/bin/env python3
"""Static instruction count of one shipped walk step, from objdump -d.

Every function in the step's call graph must compile to straight-line code: no
conditional branch, no loop back-edge, no indirect branch, and exactly ONE
return, which must be the last instruction of the function (an early return
anywhere earlier is rejected explicitly).  The only control transfers allowed
are the fixed `bl` call sites.  Under those conditions the executed count of
one call is exactly

    sum over functions f of  count(f) * mult(f)

where count(f) is the number of instructions in f (prologue, body and the
single `ret`, all executed once per call) and mult(f) is the number of times
the call graph enters f per walk step.

mult(f) is the DAG path count: with the call graph taken as a DAG, a call site
of g in f contributes mult(f) entries of g, so mult(g) = sum over callers f of
mult(f) * (number of `bl g` sites in f).  Summing over converging paths (rather
than taking a per-edge maximum) is required; --selftest contains a converging
call graph where the two differ.

The call graph must be acyclic: a cycle would be a loop through calls and is
rejected.  Alignment padding after the single `ret` (nops only) is dropped;
a nop anywhere before the `ret` is rejected.

Usage: static_count.py <binary> [root-symbol]
       static_count.py --selftest
Exit code 0 iff every straight-line assertion holds.
"""
import re
import subprocess
import sys

COND = re.compile(r'^(b\.[a-z]+|cbz|cbnz|tbz|tbnz)$')
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
    return funcs


def analyse(funcs, root):
    calls = {}          # caller -> {callee: n_call_sites}
    external = []
    for name, insns in funcs.items():
        for addr, mn, ops in insns:
            if mn == 'bl':
                m = re.match(r'[0-9a-f]+ <([^>+]+)>$', ops.strip())
                if not m:
                    external.append('%s: call target is not a plain symbol (%s)' % (name, ops))
                    continue
                callee = m.group(1)
                calls.setdefault(name, {}).setdefault(callee, 0)
                calls[name][callee] += 1
                if callee not in funcs:
                    external.append('%s: call to symbol outside the binary: %s' % (name, callee))
    # Reachable set, cycle detection and topological order (Kahn from the root).
    reach = {root}
    stack = [root]
    while stack:
        f = stack.pop()
        for g in calls.get(f, {}):
            if g in funcs and g not in reach:
                reach.add(g)
                stack.append(g)
    indeg = {f: 0 for f in reach}
    for f in reach:
        for g in calls.get(f, {}):
            if g in reach:
                indeg[g] += 1
    order, ready = [], [f for f in reach if indeg[f] == 0]
    while ready:
        f = ready.pop()
        order.append(f)
        for g in calls.get(f, {}):
            if g in reach:
                indeg[g] -= 1
                if indeg[g] == 0:
                    ready.append(g)
    problems = list(dict.fromkeys(external))
    if len(order) != len(reach):
        problems.append('call graph has a cycle (loop through calls); straight-line path '
                        'counting is undefined')
    # DAG path multiplicity: sum over converging call paths.
    mult = {root: 1}
    for f in order:
        n = mult.get(f, 0)
        for g, c in sorted(calls.get(f, {}).items()):
            if g in funcs:
                mult[g] = mult.get(g, 0) + n * c
    for name in sorted(mult):
        insns = funcs[name]
        rets = [i for i, (_, mn, _) in enumerate(insns) if mn == 'ret']
        if len(rets) == 0:
            problems.append('%s: no return instruction' % name)
            continue
        if len(rets) > 1:
            problems.append('%s: %d returns; an early return is rejected (exactly one, last, '
                            'is required)' % (name, len(rets)))
        if rets[0] != len(insns) - 1:
            tail = insns[rets[-1] + 1:]
            if len(rets) == 1 and tail and all(mn == 'nop' for _, mn, _ in tail):
                # Alignment padding after the single `ret`: never executed, drop it.
                funcs[name] = insns[:rets[-1] + 1]
                insns = funcs[name]
            else:
                problems.append('%s: early return at instruction %d of %d'
                                % (name, rets[0], len(insns)))
        elif rets[-1] != len(insns) - 1:
            problems.append('%s: last instruction is not `ret`' % name)
        for i, (addr, mn, ops) in enumerate(insns):
            if COND.match(mn):
                problems.append('%s: conditional branch %s at %x' % (name, mn, addr))
            if mn in ('b', 'br', 'blr', 'braa', 'brab'):
                problems.append('%s: unconditional control transfer %s at %x' % (name, mn, addr))
            if mn == 'nop' and i < len(insns) - 1:
                problems.append('%s: nop at %x inside the counted range' % (name, addr))
    return mult, calls, problems, len(reach) != len(order)


def report(funcs, root, mult, calls, problems):
    total = 0
    print('root: %s' % root)
    print('%-24s %8s %8s %11s' % ('function', 'insns', 'calls', 'insns/step'))
    for f in sorted(mult, key=lambda x: -mult[x]):
        n = len(funcs[f])
        total += n * mult[f]
        print('%-24s %8d %8d %11d' % (f, n, mult[f], n * mult[f]))
        for g, c in sorted(calls.get(f, {}).items()):
            print('    %-20s -> %-24s x%d' % (f, g, c))
    print('%-24s %8s %8s %11d' % ('TOTAL', '', '', total))
    if problems:
        print('\nSTRAIGHT-LINE ASSERTION FAILED:')
        for p in problems:
            print('  ' + p)
        return total, 1
    print('\nOK: every reached function is straight-line (no conditional or unconditional '
          'branch, no call cycle, exactly one return and it is the last instruction, no nop '
          'in the counted range); multiplicities are DAG path sums.')
    return total, 0


def selftest():
    """Crafted counterexamples, one per audit rule.  Fails loudly if any are missed."""
    def ins(addr, mn, ops=''):
        return (addr, mn, ops)
    nop = ins(0, 'nop')
    ok = True

    # (a) converging call paths: root -> a, b(x2); a -> c; b -> c; c -> d.
    #     DAG path count of d is 1*1 + 2*1 = 3; a per-edge maximum would give 2.
    funcs = {
        'root': [ins(0, 'bl', '0 <a>'), ins(4, 'bl', '0 <b>'), ins(8, 'bl', '0 <b>'), ins(12, 'ret')],
        'a': [ins(0, 'bl', '0 <c>'), ins(4, 'ret')],
        'b': [ins(0, 'bl', '0 <c>'), ins(4, 'ret')],
        'c': [ins(0, 'bl', '0 <d>'), ins(4, 'ret')],
        'd': [ins(0, 'add', 'x0, x0, x0'), ins(4, 'ret')],
    }
    mult, _, probs, _ = analyse(funcs, 'root')
    if mult.get('d') != 3 or probs:
        print('SELFTEST FAIL (a): converging paths gave d=%r, problems=%r' % (mult.get('d'), probs))
        ok = False

    # (b) early return: a function with a `ret` in the middle and one at the end.
    funcs = {'root': [ins(0, 'bl', '0 <f>'), ins(4, 'ret')],
             'f': [ins(0, 'add', 'x0, x0, x0'), ins(4, 'ret'), ins(8, 'add', 'x0, x0, x0'), ins(12, 'ret')]}
    _, _, probs, _ = analyse(funcs, 'root')
    if not any('early return' in p or 'returns' in p for p in probs):
        print('SELFTEST FAIL (b): early return not rejected: %r' % probs)
        ok = False

    # (c) forward unconditional branch inside a counted function.
    funcs = {'root': [ins(0, 'b', '0 <root+8>'), ins(4, 'ret')]}
    _, _, probs, _ = analyse(funcs, 'root')
    if not any('unconditional control transfer' in p for p in probs):
        print('SELFTEST FAIL (c): forward branch not rejected: %r' % probs)
        ok = False

    # (d) conditional branch.
    funcs = {'root': [ins(0, 'cbz', 'x0, 0 <root+8>'), ins(4, 'ret')]}
    _, _, probs, _ = analyse(funcs, 'root')
    if not any('conditional branch' in p for p in probs):
        print('SELFTEST FAIL (d): conditional branch not rejected: %r' % probs)
        ok = False

    # (e) call cycle (a loop through calls).
    funcs = {'root': [ins(0, 'bl', '0 <f>'), ins(4, 'ret')],
             'f': [ins(0, 'bl', '0 <root>'), ins(4, 'ret')]}
    _, _, probs, _ = analyse(funcs, 'root')
    if not any('cycle' in p for p in probs):
        print('SELFTEST FAIL (e): call cycle not rejected: %r' % probs)
        ok = False

    # (g) alignment padding after the single ret is dropped, not counted.
    funcs = {'root': [ins(0, 'add', 'x0, x0, x0'), ins(4, 'ret'), ins(8, 'nop')]}
    mult, _, probs, _ = analyse(funcs, 'root')
    if probs or len(funcs['root']) != 2:
        print('SELFTEST FAIL (g): trailing padding not handled: %r %r' % (probs, funcs['root']))
        ok = False

    # (f) nop inside the counted range.
    funcs = {'root': [ins(0, 'nop'), ins(4, 'ret')]}
    _, _, probs, _ = analyse(funcs, 'root')
    if not any('nop' in p for p in probs):
        print('SELFTEST FAIL (f): nop not rejected: %r' % probs)
        ok = False

    print('selftest: %s (converging-path count, early return, forward branch, conditional '
          'branch, call cycle, nop, trailing padding)' % ('OK' if ok else 'FAILED'))
    return 0 if ok else 1


def main():
    if len(sys.argv) >= 2 and sys.argv[1] == '--selftest':
        return selftest()
    path = sys.argv[1]
    root = sys.argv[2] if len(sys.argv) > 2 else 'rho_step'
    funcs = disasm(path)
    if root not in funcs:
        cands = [n for n in funcs if n.startswith('rho_step')]
        if not cands:
            print('root symbol not found', file=sys.stderr)
            return 2
        root = cands[0]
    mult, calls, problems, _ = analyse(funcs, root)
    _, rc = report(funcs, root, mult, calls, problems)
    return rc


if __name__ == '__main__':
    sys.exit(main())

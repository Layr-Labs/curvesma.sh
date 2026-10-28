"""Paired reduction words. usage: gen48.py <predecessor measurement dir> <out dir> [--same]
Rewrites the six product and square sites of bloop_branchy.h and floop.h: m*p3 by mul/umulh against x29 = p3,
m_k>>32 and m_(k+1)<<32 joined by one extr, the last m>>32 a shifted adds; the square's first cross chain ends at limb
5. Each site (first instruction, length) occurs once; --same writes the sites back unchanged."""
import os, re, sys
Q, L, M = ['x0', 'x2', 'x18', 'x30'], ['x19', 'x20', 'x21', 'x22'], ['x3', 'x4', 'x5', 'x6']
SITES = {'bloop_branchy.h': [('ldp x7,x8,[x23,#32]', 96, L, ['x7', 'x8'] * 2, Q, {0: ['ldp x7,x8,[x23,#32]'], 2: ['ldp x7,x8,[x23,#48]']}),
                             ('ldp x5,x6,[x23,#16]', 96, M, L, L, {0: ['ldp x5,x6,[x23,#16]', 'ldp x3,x4,[x23],#-64']}),
                             ('mul x9,x3,x0', 94, M, Q, Q, {}), ('mul x4,x2,x0', 75, Q, None, M, {}),
                             ('mul x5,x8,x0', 94, ['x8', 'x9', 'x10', 'x11'], Q, Q, {})],
         'floop.h': [('mul x9,x3,x19', 94, M, L, L, {})]}
CUT = {'bloop_branchy.h': 21, 'floop.h': 4}

def emitter():
    ins = []
    def e(op, *s, sh=None, d=None):
        d = d if d is not None else '%%%d' % len(ins)
        ins.append((op, d, s, sh))
        return d
    return ins, e

def product(A, B, pre):
    ins, e = emitter()
    w = {}
    for j in range(4):
        for r in pre.get(j, []):
            e('raw', d=r)
        b = B[j]
        if j == 0:
            m = e('mul', A[0], b)
            for k in (1, 2, 3):
                w[k] = e('adcs' if k > 1 else 'adds', e('mul', A[k], b), e('umulh', A[k - 1], b))
            w[4] = e('adc', e('umulh', A[3], b), 'xzr')
        else:
            t = e('umulh', A[3], b)
            for k in range(4):
                w[j + k] = e('adcs' if k else 'adds', w[j + k], e('mul', A[k], b))
            t = e('adc', t, 'xzr')
            if j == 1:
                w[1] = e('adds', w[1], m, sh='lsl')
            for k, r in enumerate([e('extr', w[j], m), r17, r1]):
                w[j + 1 + k] = e('adcs' if k or j == 1 else 'adds', w[j + 1 + k], r)
            t = e('adc', t, 'xzr')
            for k in range(3):
                w[j + 1 + k] = e('adcs' if k else 'adds', w[j + 1 + k], e('umulh', A[k], b))
            w[j + 4] = e('adc', t, 'xzr')
            m = w[j]
        r17, r1 = e('mul', m, 'x29'), e('umulh', m, 'x29')
    w[5] = e('adds', w[5], m, sh='lsr')
    w[6], w[7] = e('adcs', w[6], r17), e('adcs', w[7], r1)
    return ins, [w[k] for k in (4, 5, 6, 7)]

def square(a):
    ins, e = emitter()
    s = {1: e('mul', a[0], a[1])}
    s[2] = e('adds', e('mul', a[0], a[2]), e('umulh', a[0], a[1]))
    s[3] = e('adcs', e('mul', a[0], a[3]), e('umulh', a[0], a[2]))
    s[4] = e('adcs', e('umulh', a[0], a[3]), e('mul', a[1], a[3]))
    s[5] = e('adc', e('umulh', a[1], a[3]), 'xzr')
    s[3] = e('adds', s[3], e('mul', a[1], a[2]))
    s[4] = e('adcs', s[4], e('umulh', a[1], a[2]))
    s[5] = e('adcs', s[5], e('mul', a[2], a[3]))
    s[6] = e('adc', e('umulh', a[2], a[3]), 'xzr')
    for k in range(1, 7):
        s[k] = e('adcs' if k > 1 else 'adds', s[k], s[k])
    s[7] = e('adc', 'xzr', 'xzr')
    T = {0: e('mul', a[0], a[0])}
    for k in range(1, 8):
        T[k] = e('adds' if k == 1 else 'adc' if k == 7 else 'adcs', s[k], e('umulh' if k % 2 else 'mul', a[k // 2], a[k // 2]))
    for k in range(4):
        m = T[k]
        r17, r1 = e('mul', m, 'x29'), e('umulh', m, 'x29')
        if k == 0:
            T[1] = e('adds', T[1], m, sh='lsl')
            T[2], T[3] = e('adcs', T[2], e('extr', T[1], m)), e('adcs', T[3], r17)
        elif k < 3:
            p, q = (3, 8) if k == 1 else (8, 9)
            T[p] = e('adds', T[p], e('extr', T[k + 1], m))
            T[q] = e('adcs', T[q], r17)
        else:
            T[9] = e('adds', T[9], m, sh='lsr')
            T[10] = e('adcs', T[10], r17)
        T[k + 8] = e('adc', r1, 'xzr')
    return ins, [e('adcs' if k else 'adds', T[8 + k], T[4 + k]) for k in range(4)]

def allocate(ins, outs, pins, pool):
    first, last = {}, {}
    for i, (op, d, s, sh) in enumerate(ins):
        for x in s:
            last[x] = i
        first[d] = i
    for o in outs:
        last[o] = len(ins)
    pin = dict(zip(outs, pins))
    held = {pin[o]: first[o] for o in outs}
    phys, busy = {}, {}
    for i, (op, d, s, sh) in enumerate(ins):
        for x in s:
            if last[x] == i and x in phys:
                busy.pop(phys[x], None)
        if op == 'raw':
            continue
        r = pin.get(d) or next(c for c in pool if c not in busy and last[d] <= held.get(c, 1 << 30))
        assert r not in busy
        phys[d] = r
        busy[r] = d
    R = lambda x: phys.get(x, x)
    return [d if op == 'raw' else '%s %s' % (op, ','.join([R(d)] + [R(x) for x in s]))
            + (',#32' if op == 'extr' else ',%s #32' % sh if sh else '') for op, d, s, sh in ins]

def written(code):
    w = set()
    for l in code:
        op, _, a = l.partition(' ')
        a = a.split(',')
        w |= set(a[:2]) if op == 'ldp' else {a[0]}
    return w

def rewrite(text, sites, same):
    ins = re.search(r'#define \w+_S\(S\) (.*)', text).group(1).replace('" S "', '').strip('"').split('\\n')
    at = 0
    for first, n, A, B, out, pre in sites:
        i = ins.index(first, at)
        old, at = ins[i:i + n], i + n
        assert not any(l.endswith(':') or l[0] == 'b' for l in old)
        pool = sorted(written(old) - set(A + (B or [])) - written(sum(pre.values(), [])) - {'x23'}, key=lambda r: int(r[1:]))
        new = old if same else allocate(*(product(A, B, pre) if B else square(A)), out, pool)
        o = '\\n'.join(old) + '\\n'
        assert text.count(o) == 1
        text = text.replace(o, '\\n'.join(new) + '\\n')
    return text

if __name__ == '__main__':
    src, dst = sys.argv[1:3]
    for f, sites in SITES.items():
        t = open(os.path.join(src, f)).read()
        u = rewrite(t, sites, '--same' in sys.argv)
        if '--same' not in sys.argv:
            u = u.replace('): ', '), then gen48.py (paired reduction words): ', 1)
            n = lambda x: re.search(r'#define \w+_S\(S\) (.*)', x).group(1).count('\\n')
            assert n(t) - n(u) == CUT[f]
        open(os.path.join(dst, f), 'w').write(u)
    print('ok')

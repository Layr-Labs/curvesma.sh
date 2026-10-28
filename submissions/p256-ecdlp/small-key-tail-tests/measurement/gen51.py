"""Deferred x3/y3' tests: derive bloop_branchy.h and step.c from the shared-row revision's (proof.md, "Deferred x3
and y3' tests"). usage: gen51.py <shared-row measurement dir> <output dir>; every edit must match once."""
import re, sys
L = lambda n: n + '" S "_%='
Y3, X3, YN = ['x0', 'x2', 'x18', 'x30'], ['x3', 'x4', 'x5', 'x6'], ['x13', 'x14', 'x15', 'x16']


def once(t, old, new):
    assert t.count(old) == 1, old[:60]
    return t.replace(old, new)


def seq(s, old, new):
    i = [k for k in range(len(s)) if s[k:k + len(old)] == old]
    assert len(i) == 1, old[0]
    return s[:i[0]] + new + s[i[0] + len(old):]


def fold(Y, n, ld):
    a, b, c, d = Y
    return ld + ['cmp x30,x11', 'b.hs ' + L('Lp' + n), 'subs %s,x8,x0' % a, 'sbcs %s,x9,x2' % b, 'sbcs %s,x10,x18' % c,
                 'sbc %s,x11,x30' % d, 'b ' + L('Lq' + n), L('Lp' + n) + ':', 'subs %s,x0,x8' % a, 'sbcs %s,x2,x9' % b,
                 'sbcs %s,x18,x10' % c, 'sbcs %s,x30,x11' % d, 'b.cs ' + L('Lq' + n)] + fix(Y, 'add', 5) + [L('Lq' + n) + ':']


def fix(r, op, w):
    s, a = ('subs', 'adc') if op == 'add' else ('adds', 'sbc')
    return ['%s %s,%s,#1' % (s, r[0], r[0]), '%ss %s,%s,x28' % (a, r[1], r[1]), '%ss %s,%s,xzr' % (a, r[2], r[2]),
            '%s %s,%s,x29' % (a, r[3], r[3]), 'add x27,x27,#%d' % w]


def canon(r, t):
    return (['adds %s,%s,#1' % (t[0], r[0]), 'sbcs %s,%s,x28' % (t[1], r[1]), 'sbcs %s,%s,xzr' % (t[2], r[2]),
             'sbcs %s,%s,x29' % (t[3], r[3])] + ['csel %s,%s,%s,hs' % (r[k], t[k], r[k]) for k in range(4)])


def macro(h, name, f):
    m = re.search(r'#define %s\(S\) "(.*)"\n' % name, h)
    return h[:m.start(1)] + '\\n'.join(f(m.group(1).split('\\n'))) + h[m.end(1):]


def body(b):
    b = seq(b, ['b.cc ' + L('Lp5')] + fix(X3, 'sub', 5) + [L('Lp5') + ':', 'cmn x6,x28', 'b.hs ' + L('Ls5'), L('Lq5') + ':'],
            ['b.cc ' + L('Lq5')] + fix(X3, 'sub', 5) + [L('Lq5') + ':'])
    old = ['sbcs x6,x6,x1', 'b.cs ' + L('Lq7'), 'subs x3,x3,#1', 'adcs x4,x4,x28', 'adcs x5,x5,xzr', 'adc x6,x6,x29',
           'add x27,x27,#5', L('Lq7') + ':', 'stp x5,x6,[x25,#-32]', 'subs x8,x3,x8', 'sbcs x9,x4,x9', 'sbcs x10,x5,x10',
           'sbcs x11,x6,x11']
    new = old[:]
    for k, v in ((0, 'sbcs x17,x6,x1'), (5, 'adc x17,x17,x29'), (8, 'stp x5,x17,[x25,#-32]'), (12, 'sbcs x11,x17,x11')):
        new[k] = v
    b = seq(b, old, new)
    b = seq(b, ['b.cc ' + L('Lp9')] + fix(Y3, 'sub', 4) + ['b ' + L('Lq9'), L('Lp9') + ':', 'cmn x30,x28',
                                                         'b.hs ' + L('Ls9'), L('Lq9') + ':'],
            ['b.cc ' + L('Lq9')] + fix(Y3, 'sub', 5) + [L('Lq9') + ':'])
    ld = ['ldp x10,x11,[x24,#16]', 'ldp x8,x9,[x24],#-64']
    b = seq(b, fold(Y3, '10', ld) + ['stp x0,x2,[x25,#-16]', 'stp x18,x30,[x25],#-64'],
            fold(YN, '10', ld) + ['stp x13,x14,[x25,#-16]', 'stp x15,x16,[x25],#-64'])
    assert sum('x17' in x for x in b[76:]) == 4
    return b


def slow(s):
    global n_fix
    T = ['x15', 'x16', 'x17', 'x1']
    z = s.index(L('Ls5') + ':')
    assert s[z:] == ([L('Ls5') + ':'] + canon(X3, T) + ['add x27,x27,#10', 'b ' + L('Lq5'), L('Ls9') + ':'] +
                     canon(Y3, T) + ['add x27,x27,#10', 'b ' + L('Lq9'), ''])
    assert s[:2] == [L('Lz0') + ':', 'b.mi 2f'] and s[z - 1] == 'b ' + L('Lq0')
    fx = (['cmn x17,x28', 'ccmn x30,x28,#2,lo', 'b.lo 2f', 'ldr x5,[x25,#32]'] +
          canon(['x3', 'x4', 'x5', 'x17'], ['x15', 'x16', 'x7', 'x1']) + ['stp x5,x17,[x25,#32]'] + canon(Y3, T))
    end = ['stp x13,x14,[x25,#48]', 'stp x15,x16,[x25,#64]', 'mov x17,#0', 'mov x30,#0', 'add x27,x27,#%c[wfix]',
           'b ' + L('Lt')]
    n_fix = 1 + len(fx) + 2 + 7 + len(end)
    return ([s[0], 'b.mi ' + L('Lx')] + s[2:z] + [L('Lx') + ':'] + fx +
            fold(YN, '11', ['ldp x10,x11,[x24,#80]', 'ldp x8,x9,[x24,#64]']) + end + [''])


src, dst = sys.argv[1], sys.argv[2]
h = macro(macro(open(src + '/bloop_branchy.h').read(), 'BLOOP_FIELD_S', body), 'BLOOP_SLOW_S', slow)
h = once(h, '#define BLOOP_SLOW_B', '#define BLOOP_FIX_N %d\n#define BLOOP_SLOW_B' % n_fix)
h = once(h, 'then gen49.py (shared rows)', 'then gen49.py (shared rows), then gen51.py (deferred x3/y3\' tests)')
open(dst + '/bloop_branchy.h', 'w').write(h)
c = open(src + '/step.c').read()
for old, new in (('ccmn x22,x28,#10,ne\\nb.hs Lz0', 'ccmn x17,x28,#10,ne\\nccmn x30,x28,#10,lo\\nccmn x22,x28,#10,lo\\nb.hs Lz0'),
                 ('BL_TAIL_P(S) "ldr x6', 'BL_TAIL_P(S) "Lt" S "_%=:\\nldr x6'),
                 ('BL_TAIL_A(S) "ldp x6', 'BL_TAIL_A(S) "Lt" S "_%=:\\nldp x6'),
                 ('BL_CAP BL_DP "ccmp " BL_F ",x6,#4,ne\\n" BL_INV(S) \\', 'BL_DP "ccmp " BL_F ",x6,#4,ne\\n" BL_INV(S) BL_CAP \\'),
                 ('tst x3,#%c[dpf]\\n"\n', 'tst x3,#%c[dpf]\\n"\n#define BL_DP_N 1\n'),
                 ('tst x8,#%c[dpf]\\n"\n', 'tst x8,#%c[dpf]\\n"\n#define BL_DP_N 2\n'),
                 ('#define BL_FM ""\n', '#define BL_FM ""\n#define BL_FM_N 0\n'),
                 ('#%c[fm]\\n"\n', '#%c[fm]\\n"\n#define BL_FM_N 1\n'),
                 ('[cap] "i"(CAP_PERIODS) \\', '[cap] "i"(CAP_PERIODS), \\\n                  [wfix] "i"(BLOOP_FIX_N + 7 + BL_FM_N + BL_DP_N) \\')):
    c = once(c, old, new)
assert c.count('bfi x9,x0,#6,#1') == 2
open(dst + '/step.c', 'w').write(c.replace('bfi x9,x0,#6,#1', 'bfi x9,x13,#6,#1'))
print('ok', n_fix)

"""Top-limb distinguished test: derive bloop_branchy.h and step.c from the deferred-tests revision's (proof.md,
"Distinguished test on the top limb"). usage: gen55.py <deferred-tests measurement dir> <output dir>; every edit must
match once."""
import re, sys


def once(t, old, new):
    assert t.count(old) == 1, old[:60]
    return t.replace(old, new)


src, dst = sys.argv[1], sys.argv[2]
h = open(src + '/bloop_branchy.h').read()
m = re.search(r'#define BLOOP_SLOW_S\(S\) "(.*)"\n', h)
s = m.group(1).split('\\n')
z = s.index('stp x5,x17,[x25,#32]')
y = ['adds x15,x0,#1', 'sbcs x16,x2,x28', 'sbcs x17,x18,xzr', 'sbcs x1,x30,x29', 'csel x0,x15,x0,hs', 'csel x2,x16,x2,hs',
     'csel x18,x17,x18,hs', 'csel x30,x1,x30,hs']
assert s[z + 1:z + 9] == y and s.count('mov x17,#0') == 1
s[z + 3], s[z + 7] = 'sbcs x5,x18,xzr', 'csel x18,x5,x18,hs'
s[s.index('mov x17,#0')] = 'and x17,x17,#%c[dpf]'
assert [l for l in s[z + 1:] if 'x17' in l] == ['and x17,x17,#%c[dpf]']
assert [l for l in s[z + 1:] if 'x5' in l] == ['sbcs x5,x18,xzr', 'csel x18,x5,x18,hs']
h = h[:m.start(1)] + '\\n'.join(s) + h[m.end(1):]
h = once(h, "then gen51.py (deferred x3/y3' tests)", "then gen51.py (deferred x3/y3' tests), then gen55.py (top-limb DP test)")
open(dst + '/bloop_branchy.h', 'w').write(h)
c = open(src + '/step.c').read()
i, j = c.index('#if DP_BITS + ADDER_BITS <= 64\n#define IS_DP'), c.index('typedef struct { fe x, y; } pt;')
c = c[:i] + ('_Static_assert(DP_BITS >= 1 && DP_BITS <= 63, "top-limb DP mask leaves bit 0 clear");\n'
             '#define IS_DP(kl, p) (((p)->x.v[3] >> (64 - DP_BITS)) == 0)\n#define IS_DP1 IS_DP\n\n') + c[j:]
i, j = c.index('#if DP_BITS + ADDER_BITS <= 64\n#define BL_DP'), c.index('#if FILTER_MASK == ~0ULL\n#define BL_F')
c = c[:i] + ('#define BL_DP "tst x17,#%c[dpf]\\n"\n#define BL_DP_N 1\n'
             '#define BL_DPF (((1ULL << DP_BITS) - 1) << (64 - DP_BITS))\n') + c[j:]
assert c.count('IS_DP') == 7 and c.count('extr') == 0
open(dst + '/step.c', 'w').write(c)
P3 = 0xffffffff00000001
for d in range(1, 64):
    f = ((1 << d) - 1) << (64 - d)
    assert (P3 & f) < P3 and f & 1 == 0
print('ok')

"""k5: lazy body -> side_k2.S, -26/lane (k2:11 k3:2 k4:12 T18:1; proof.md s1).
check() asserts byte equality; semantics: native_check.py. Stripped emit."""
import re
import sys
from pathlib import Path

root = Path(__file__).resolve().parent
if (root / 'files').is_dir():
    root = root / 'files'
sys.path.insert(0, str(root))
import lazy_witness

K2_SLOT = '[x29,#112]'


def _regions(text):
    m1 = text.index('// ---- lazy first Montgomery product')
    e1 = text.index('// ---- end __mul_mont_384', m1)
    while text[e1 - 1] in '\t ':
        e1 -= 1
    m2 = text.index('// ---- fused raw second Montgomery product')
    e2 = text.index('// ---- end __mul_mont_384', m2)
    while text[e2 - 1] in '\t ':
        e2 -= 1
    assert m1 < e1 < m2 < e2
    return (m1, e1), (m2, e2)


def _t1(text):
    out, n = [], 0
    lines = text.splitlines(keepends=True)
    i = 0
    while i < len(lines):
        if lines[i].split('//')[0].strip() == 'mov\tx17,xzr':
            j = i + 1
            while True:
                code = lines[j].split('//')[0].strip()
                assert not re.match(r'\S+\s+x17,', code), 'x17 clobbered'
                if code == 'adc\tx4,x17,xzr':
                    break
                j += 1
            lines[j] = lines[j].replace('adc\tx4,x17,xzr',
                                        'adc\tx4,xzr,xzr')
            n += 1
            i += 1
            continue
        out.append(lines[i])
        i += 1
    assert n == 2, n
    return ''.join(out)


def _t2(text):
    (m1, e1), (m2, e2) = _regions(text)
    parts = [text[:m1], text[m1:e1], text[e1:m2], text[m2:e2], text[e2:]]
    for k in (1, 3):
        body = parts[k]
        bare = '\n'.join(l.split('//')[0] for l in body.splitlines())
        assert not re.search(r'\bx18\b', bare), 'x18 used in product'
        assert body.count('\tmul\tx4,x4,') == 6
        body = body.replace('\tmul\tx4,x4,', '\tmul\tx4,x18,')
        assert body.count('\tldr\tx4,[x29,#12*__SIZEOF_POINTER__]') == 5
        body = body.replace('\tldr\tx4,[x29,#12*__SIZEOF_POINTER__]\n', '')
        old_ep = ('\tldp\tx4,x2,[x29,#12*__SIZEOF_POINTER__]'
                  '\t// pull r_ptr\n')
        assert body.count(old_ep) == 1
        body = body.replace(
            old_ep,
            '\tldr\tx18,[x29,#112]\t// k2: restore cursor (T4)\n')
        first = body.index('\tmul\t')
        spill = ('\tstr\tx18,[x29,#112]\t// k2: spill cursor\n'
                 '\tldr\tx18,[x29,#96]\t// k2: n0 pinned\n')
        body = body[:first] + spill + body[first:]
        body += ('\tldr\tx2,[x18,#-32]\t// k2: out via cursor (T4)\n'
                 '\tldr\tx4,[x29,#96]\t// k2: n0 live again\n')
        parts[k] = body
    return ''.join(parts)


def _t3(text):
    sbc = '\tsbc\tx0,xzr,xzr'
    assert text.count(sbc) == 1
    text, n = re.subn(r'\tsbc\tx0,xzr,xzr[ \t]*// all ones iff D < 0\n',
                      '', text)
    assert n == 1
    assert text.count('\tand\tx') == 6
    text = re.sub(r'\tand\tx\d+,x\d+,x0\n', '', text)
    anchor = '\tldp\tx17,x3,[x29,#160]\n'
    assert text.count(anchor) == 1
    text = text.replace(anchor, anchor + '\tb.hs\t.Lk2_x_nonneg\n')
    cond = '\tsubs\tx25,x11,x5\n'
    assert text.count(cond) == 1
    text = text.replace(cond, '.Lk2_x_nonneg:\n' + cond)
    return text


def _t4(text):
    text, n = re.subn(r'\tstr\tx0,\[x29,#104\][ \t]*'
                      r'// the inlined multiplication reloads \(n0, out\)\n',
                      '', text)
    assert n == 1
    bare = '\n'.join(l.split('//')[0] for l in text.splitlines())
    assert '[x29,#104]' not in bare, 'out slot still referenced'
    return text


TAIL_N0 = '\tldr\tx4,[x29,#96]\t// k2: n0 live again\n'


def _t5(text):
    """Drop second product tail n0 reload (x4 unread after; ldr sets no flags)."""
    assert text.count(TAIL_N0) == 2
    head, sep, tail = text.partition(TAIL_N0)
    assert sep and TAIL_N0 in tail
    tail = tail.replace(TAIL_N0, '', 1)
    assert TAIL_N0 not in tail
    return head + sep + tail


def _t6(text):
    """Fold mu-pointer add: 6 second-product mu loads go [sp,#96..136] (x2 dead)."""
    pre, marker, post = text.partition('// ---- fused raw second Montgomery product')
    assert marker
    unit = '\tadd\tx2,sp,#96\n\tldr\tx17,[x2]\n'
    assert pre.count(unit) == 1 and post.count('\tadd\tx2,sp,#96\n') == 0
    pre = pre.replace(unit, '\tldr\tx17,[sp,#96]\n', 1)
    body, endmark, rest = post.partition('// ---- end __mul_mont_384')
    assert endmark
    for n in range(1, 6):
        old = '\tldr\tx17,[x2,8*%d]\n' % n
        assert body.count(old) == 1, n
        body = body.replace(old, '\tldr\tx17,[sp,#%d]\n' % (96 + 8 * n), 1)
    assert '\tldr\tx17,[x2' not in body
    assert pre.count('\tldr\tx17,[x2') == 6  # first-product inv loads untouched
    return pre + marker + body + endmark + rest


TOPC = '\tadc\tx25,x17,xzr\n'
REST = '\tldr\tx18,[x29,#112]\t// k2: restore cursor (T4)\n'
TOP17 = '\tadc\tx17,x17,xzr\n'
ABT17 = '\tadc\tx17,xzr,xzr\n'


def _bl(s):
    return s.split('//')[0]
WR = {'adds', 'subs', 'cmp', 'tst'}
NEU = {'mov', 'add', 'sub', 'mul', 'umulh', 'madd', 'lsl', 'lsr', 'asr',
       'ubfx', 'bfi', 'and', 'orr', 'eor', 'ldr', 'str', 'ldp', 'stp',
       'adrp', 'adr', 'nop'}


def _flagdead(bare, i):
    while True:
        c = bare[i].strip()
        i += 1
        if not c or c.endswith(':'):
            continue
        op = c.split()[0]
        if op in WR:
            return True
        assert op in NEU, c


def _casc(text):
    (m1, e1), (m2, e2) = _regions(text)
    parts = [text[:m1], text[m1:e1], text[e1:m2], text[m2:e2], text[e2:]]
    for k in (1, 3):
        s = parts[k]
        assert s.count(TOPC) == 1  # T7/T8 top carry
        s = s.replace(TOPC, '')
        assert s.count(REST + TOP17) == 1  # T12/T13 above final adds
        s = s.replace(REST + TOP17, REST)
        assert s.count(TOP17 + '\tsubs') == 5  # T14/T16 round-6 q.p tail
        i = s.rfind(TOP17 + '\tsubs')
        s = s[:i] + s[i + len(TOP17):]
        assert s.count(ABT17) == 5  # T15/T17 round-6 a.b tail
        i = s.rfind(ABT17)
        s = s[:i] + s[i + len(ABT17):]
        parts[k] = s
    out = ''.join(parts)
    b = [l.split('//')[0] for l in out.splitlines(keepends=True)]
    r1 = b.index(_bl(REST))
    r2 = b.index(_bl(REST), r1 + 1)
    oc = _bl('\tldr\tx2,[x18,#-32]\t// k2: out via cursor (T4)\n')
    o1 = b.index(oc)
    o2 = b.index(oc, o1 + 1)
    w25 = b.index('\tadc\tx25,xzr,xzr\n')
    cmp_ = b.index('\tcmp\tx18,x30\n')
    f1 = max(i for i, l in enumerate(b) if l == '\tadcs\tx25,x25,xzr\n' and i < r1)
    f2 = max(i for i, l in enumerate(b) if l == '\tadcs\tx25,x25,xzr\n' and i < r2)
    for a, z in ((f1, r1), (r1, o1), (f2, r2), (r2, o2)):
        assert not any(re.search(r'\bx17\b', l) for l in b[a + 1:z]), (a, z)
    assert not any(re.search(r'\bx25\b', l) for l in b[o1 + 1:w25])
    assert not any(re.search(r'\bx2[57]\b', l) for l in b[o2 + 1:cmp_])
    for i in (f1 + 1, f2 + 1, r1 + 1, r2 + 1, o1, o2):
        _flagdead(b, i)
    s400 = max(i for i, l in enumerate(b) if l.startswith('\tsubs\txzr,') and i < r1)
    s1100 = max(i for i, l in enumerate(b) if l.startswith('\tsubs\txzr,') and i < r2)
    _flagdead(b, s400)
    _flagdead(b, s1100)
    return out


HEAD_LDP = '\tldp\tx2,x3,[x18,#16]\t\t// T, inv\n'
MOV_X2 = '\tmov\tx2,x3\n'


def _t9(text):
    assert text.count(HEAD_LDP) == 1 and text.count(MOV_X2) == 1
    y2 = ['\tldp\tx19,x20,[x2,#48]\n', '\tldp\tx21,x22,[x2,#64]\n',
          '\tldp\tx23,x24,[x2,#80]\n']
    for old in y2:
        assert text.count(old) == 1
    b = [l.split('//')[0] for l in text.splitlines(keepends=True)]
    h = b.index(_bl(HEAD_LDP))
    m = b.index(MOV_X2)
    got = [l.strip() for l in b[h + 1:m] if re.search(r'\bx[23]\b', l)]
    assert got == [l.strip() for l in y2], got
    for old in y2:
        text = text.replace(old, old.replace('[x2,', '[x3,'))
    text = text.replace(HEAD_LDP, HEAD_LDP.replace('x2,x3', 'x3,x2'))
    text = text.replace(MOV_X2, '')
    b = [l.split('//')[0] for l in text.splitlines(keepends=True)]
    h = b.index(_bl('\tldp\tx3,x2,[x18,#16]\t\t// T, inv\n'))
    nx3 = next(i for i in range(h + 1, len(b)) if re.search(r'\bx3\b', b[i]))
    assert b[nx3].strip() == 'ldp\tx19,x20,[x3,#48]', b[nx3]
    nx2 = next(i for i in range(h + 1, len(b)) if re.search(r'\bx2\b', b[i]))
    assert b[nx2].strip() == 'ldr\tx17,[x2]', b[nx2]
    y3 = b.index('\tldp\tx23,x24,[x3,#80]\n')
    nx3b = next(i for i in range(y3 + 1, len(b)) if re.search(r'\bx3\b', b[i]))
    assert b[nx3b].strip().startswith('umulh\tx3,'), b[nx3b]
    return text


def _t10(text):
    o1 = '\tldr\tx2,[x18,#-32]\t// k2: out via cursor (T4)\n'
    assert text.count(o1) == 2
    a424 = o1 + TAIL_N0
    assert text.count(a424) == 1
    e1131 = o1 + '\t// ---- end __mul_mont_384'
    assert text.count(e1131) == 1
    pt = '\tldp\tx0,x1,[x18,#-24]\t// P, T of this lane, reloaded from args\n'
    assert text.count(pt) == 1
    x1 = ['\tldp\tx19,x20,[x0]\n', '\tldp\tx21,x22,[x0,#16]\n',
          '\tldp\tx23,x24,[x0,#32]\t// x1 stays in x19..x24 for w\n']
    xs = ["\tstp\tx11,x12,[x2]\t\t// x2 = out: x'\n",
          '\tstp\tx13,x14,[x2,#16]\n', '\tstp\tx15,x16,[x2,#32]\n']
    y1 = ['\tldp x11,x12,[x0,#48]\n', '\tldp x13,x14,[x0,#64]\n',
          '\tldp x15,x16,[x0,#80]\n']
    ys = ["\tstp\tx11,x12,[x2,#48]\t// y'\n", '\tstp\tx13,x14,[x2,#64]\n',
          '\tstp\tx15,x16,[x2,#80]\n']
    ks = ['\tstp\tx12,x13,[x2,#96]\n', '\tstr\tx14,[x2,#112]\n']
    for old in x1 + xs + y1 + ys + ks:
        assert text.count(old) == 1, old
    assert text.count('\tldr x0,[x18,#-24]\n') == 1
    b = [l.split('//')[0] for l in text.splitlines(keepends=True)]
    i718 = b.index(_bl(pt))
    i762 = b.index('\tstp\tx15,x16,[x2,#32]\n')
    (m1, e1), (m2, e2) = _regions(text)
    ln = text[:m2].count('\n')
    lx = text[:e2].count('\n')
    oc = _bl(o1)
    seg = b[i718 + 1:]
    seg = seg[:next(k for k, l in enumerate(seg) if l == oc)]
    assert [l.strip() for l in seg if re.search(r'\bx2\b', l)] == \
        [_bl(l).strip() for l in xs], 'P-carry clobbered'
    n424 = next(k for k in range(i718) if b[k] == oc)
    assert not any(re.search(r'\bx0\b', l) for l in b[n424 + 1:i718])
    endo = next(kk for kk in range(i762 + 1, len(b)) if b[kk] == oc)
    x0s = [k for k in range(i762 + 1, endo) if re.search(r'\bx0\b', b[k])]
    assert all(ln <= k <= lx for k in x0s) and \
        b[min(x0s)].strip().startswith('umulh\tx0,'), x0s
    text = text.replace(a424, o1.replace('\tldr\tx2,', '\tldr\tx0,', 1) + TAIL_N0)
    text = text.replace(e1131, o1.replace('\tldr\tx2,', '\tldr\tx0,', 1) +
                        '\t// ---- end __mul_mont_384')
    text = text.replace(pt, pt.replace('\tldp\tx0,', '\tldp\tx2,', 1))
    for old in x1:
        text = text.replace(old, old.replace('[x0', '[x2'))
    for old in xs:
        text = text.replace(old, old.replace('[x2', '[x0').replace(
            '// x2 = out', '// x0 = out'))
    text = text.replace('\tldr x0,[x18,#-24]\n', '')
    for old in y1:
        text = text.replace(old, old.replace('[x0,', '[x2,'))
    for old in ys + ks:
        text = text.replace(old, old.replace('[x2,', '[x0,'))
    b = [l.split('//')[0] for l in text.splitlines(keepends=True)]
    j718 = next(k for k, l in enumerate(b) if l.startswith('\tldp\tx2,x1,[x18,#-24]'))
    j1136 = next(k for k, l in enumerate(b) if l == '\tldp x11,x12,[x2,#48]\n')
    want = {'ldp\tx2,x1,[x18,#-24]', 'ldp\tx19,x20,[x2]',
            'ldp\tx21,x22,[x2,#16]', 'ldp\tx23,x24,[x2,#32]',
            'ldp x11,x12,[x2,#48]', 'ldp x13,x14,[x2,#64]',
            'ldp x15,x16,[x2,#80]'}
    assert {b[k].strip() for k in range(j718, j1136 + 3)
            if re.search(r'\bx2\b', b[k])} == want
    j1131 = next(k for k, l in enumerate(b)
                 if l == _bl('\tldr\tx0,[x18,#-32]\t// k2: out via cursor (T4)\n')
                 and k > j718)
    j1167 = next(k for k, l in enumerate(b) if l.startswith("\tstp\tx11,x12,[x0,#48]"))
    assert not any(re.search(r'\bx0\b', l) for l in b[j1131 + 1:j1167])
    j1169 = next(k for k, l in enumerate(b) if l.startswith('\tstp\tx15,x16,[x0,#80]'))
    j1180 = next(k for k, l in enumerate(b) if l.startswith('\tstp\tx12,x13,[x0,#96]'))
    assert not any(re.search(r'\bx0\b', l) for l in b[j1169 + 1:j1180])
    return text


def _t11(text):
    a413 = REST + '\n\tadds\tx11,x12,x26\n'
    assert text.count(a413) == 1
    assert text.count(TAIL_N0) == 1
    pair = ('\tstr\tx18,[x29,#112]\t// k2: spill cursor\n'
            '\tldr\tx18,[x29,#96]\t// k2: n0 pinned\n')
    assert text.count(pair) == 2
    o424 = '\tldr\tx0,[x18,#-32]\t// k2: out via cursor (T4)\n'
    assert text.count(o424 + TAIL_N0) == 1
    pt = '\tldp\tx2,x1,[x18,#-24]\t// P, T of this lane, reloaded from args\n'
    assert text.count(pt) == 1
    assert text.count('\tmul\tx15,x4,x17\n') == 6
    assert text.count(K2_SLOT) == 4
    (m1, e1), (m2, e2) = _regions(text)
    assert text[m1:e1].count(pair) == 1 and text[m2:e2].count(pair) == 1
    b = [l.split('//')[0] for l in text.splitlines(keepends=True)]
    pin = _bl('\tldr\tx18,[x29,#96]\t// k2: n0 pinned\n')
    i77 = next(k for k, l in enumerate(b) if l == pin)
    i783 = next(k for k in range(i77 + 1, len(b)) if b[k].strip() == 'str\tx18,[x29,#112]')
    got = sorted(l.strip() for l in b[i77 + 1:i783] if re.search(r'\bx18\b', l))
    exp = sorted(['mul\tx4,x18,x11'] * 6 + ['ldr\tx18,[x29,#112]',
                                            'ldr\tx0,[x18,#-32]',
                                            'ldp\tx2,x1,[x18,#-24]'])
    assert got == exp, got
    i425 = b.index(_bl(TAIL_N0))
    i718 = b.index(_bl(pt))
    got = sorted(l.strip() for l in b[i425:i718] if re.search(r'\bx4\b', l))
    assert got == sorted([_bl(TAIL_N0).strip()] + ['mul\tx15,x4,x17'] * 6), got
    i421 = max(k for k, l in enumerate(b) if l == '\tadcs\tx16,x25,x3\n')
    nx3 = next(k for k in range(i421 + 1, len(b)) if re.search(r'\bx3\b', b[k]))
    assert b[nx3].strip() == 'mul\tx3,x7,x15', b[nx3]
    i714 = b.index('\tadc\tx16,x23,x16\n')
    i724 = b.index('\tldp\tx17,x3,[x1,#32]\n')
    got = [l.strip() for l in b[i714 + 1:i724 + 1] if re.search(r'\bx3\b', l)]
    assert got == ['ldp\tx17,x3,[x1,#32]'], got
    text = text.replace(a413, '\n\tadds\tx11,x12,x26\n')
    text = text.replace(o424 + TAIL_N0,
                        '\tldr\tx3,[x29,#112]\t// k4: cursor to tmp\n'
                        '\tldr\tx0,[x3,#-32]\t// k4: out via tmp\n')
    text = text.replace(pt, '\tldr\tx3,[x29,#112]\t// k4: cursor to tmp\n' + pt.replace(
        '[x18,#-24]', '[x3,#-24]'))
    parts = [text[:m1], text[m1:e1], text[e1:m2], text[m2:e2], text[e2:]]
    parts[3] = parts[3].replace(pair, '', 1)
    text = ''.join(parts)
    text = text.replace('\tmul\tx15,x4,x17\n', '\tmul\tx15,x18,x17\n')
    assert text.count(K2_SLOT) == 4
    b = [l.split('//')[0] for l in text.splitlines(keepends=True)]
    i77 = next(k for k, l in enumerate(b) if l == pin)
    i1120 = b.index(_bl(REST))
    for k in range(i77 + 1, i1120):
        if re.search(r'\bx18\b', b[k]):
            assert re.match(r'\tmul\tx\d+,x18,', b[k]), b[k]
    return text


def _t18(text):
    c = '\tldr\tx3,[x29,#112]\t// k4: cursor to tmp\n'
    o = '\tldr\tx0,[x3,#-32]\t// k4: out via tmp\n'
    p = '\tldp\tx2,x1,[x3,#-24]\t// P, T of this lane, reloaded from args\n'
    assert text.count(c) == 2 and text.count(c + o) == 1 and text.count(c + p) == 1
    b = [l.split('//')[0] for l in text.splitlines(keepends=True)]
    qo, pl = b.index(_bl(o)), b.index(_bl(p))
    q1 = max(k for k, l in enumerate(b[:qo]) if l == '\tumulh\tx3,x10,x4\n')
    q2 = next(k for k in range(pl + 1, len(b)) if re.search(r'\bx4\b', b[k]))
    assert b[q2].strip() == 'mul\tx4,x18,x19'
    assert not any(re.search(r'\bx4\b', l) for l in b[q1 + 1:pl] + b[pl + 1:q2])
    m3 = next(k for k in range(qo + 1, len(b)) if re.search(r'\bx3\b', b[k]))
    l3 = next(k for k in range(pl + 1, len(b)) if re.search(r'\bx3\b', b[k]))
    assert b[m3].strip() == 'mul\tx3,x7,x15'
    assert b[l3].strip() == 'ldp\tx17,x3,[x1,#32]'
    text = text.replace(c + o, '\tldr\tx4,[x29,#112]\n\tldr\tx0,[x4,#-32]\n')
    text = text.replace(c + p, '\tldp\tx2,x1,[x4,#-24]\n')
    assert c not in text and text.count(K2_SLOT) == 3
    return text


def _strip(text):
    out = []
    for l in text.splitlines(keepends=True):
        if '//' not in l or '// ----' in l:
            out.append(l)
        else:
            b = l.split('//')[0].rstrip()
            if b:
                out.append(b + '\n')
    return ''.join(out)


def generate(text):
    lazy, _ = lazy_witness.generate(text)
    assert K2_SLOT not in lazy
    out = _t18(_t11(_t10(_t9(_casc(_t6(_t5(_t4(_t3(_t2(_t1(lazy)))))))))))
    out = _strip(out)
    return out, {'t1_xzr_carry': 2, 't2_n0_pinned_net': 2,
                 't3_branchy_x_worst': 6, 't4_out_via_cursor': 1,
                 't5_dead_tail_n0': 1, 't6_mu_ptr_fold': 1,
                 't7_t8_dead_top_carry': 2, 't12_t13_dead_top17': 2,
                 't14_t16_dead_qptail': 2, 't15_t17_dead_abtail': 2,
                 't9_head_ptr_swap': 1, 't10_ptr_carry': 1,
                 't11_lane_pinned_n0': 2, 't18_cursor_via_x4': 1,
                 'worst_saved_per_lane': 26}


def check(root):
    text, meta = generate((root / 'measure/side_lanes.S').read_text())
    assert text == (root / 'side_k2.S').read_text()
    return meta

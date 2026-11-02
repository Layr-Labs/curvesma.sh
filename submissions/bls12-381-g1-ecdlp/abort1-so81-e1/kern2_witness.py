"""Second kernel transform: xzr carry, n0 pinned in x18, branchy x'.

Transforms the lazy side-output body (lazy_witness.generate output) into
side_k2.S, saving 10 instructions per lane on the worst path:
 T1 (-2): drop the two `mov x17,xzr` zeroings; the consuming
    `adc x4,x17,xzr` becomes `adc x4,xzr,xzr` (x17 is dead: reloaded next).
 T2 (-2): pin n0 in x18 across each Montgomery product (x18 is untouched
    inside both products), deleting 10 n0 reloads at the cost of 4 cursor
    spill/restore + 2 n0 loads + 2 x4 restores; epilogues load out only.
 T3 (-6): branchy x' normalization (b.hs over the +2p adds) replacing the
    mask sequence (sbc + 6 and); worst path (negative D) 28 vs 34.
 T4 (-1): out pointer via [x18,#-32] instead of the frame slot: the lane-head
    out spill is deleted and each epilogue restores the cursor first, then
    loads out cursor-relative (same count as the frame loads they replace).
Byte equality with side_k2.S is asserted by check(). Correctness of every
path is audited differentially by native_check.py (whole-integer reference).
"""
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


def generate(text):
    lazy, _ = lazy_witness.generate(text)
    assert K2_SLOT not in lazy
    out = _t4(_t3(_t2(_t1(lazy))))
    return out, {'t1_xzr_carry': 2, 't2_n0_pinned_net': 2,
                 't3_branchy_x_worst': 6, 't4_out_via_cursor': 1,
                 'worst_saved_per_lane': 11}


def check(root):
    text, meta = generate((root / 'measure/side_lanes.S').read_text())
    assert text == (root / 'side_k2.S').read_text()
    return meta

"""Derive measure/side_lanes.S (925 charge, 932 worst) from the validated k27r1 929 listing (SHA b4c7a7c2...).
Usage: python3 derive_k925.py <k27r1 side_lanes.S>. One exact text transform (proof.md s24), every edit asserted:
K3 (+13) inserts the frame-slot init, expensive-arm counter, node-local aggregate and counter BSS symbol."""
import hashlib, sys

INPUT_SHA = 'b4c7a7c2d17e6decfc8c0ee93f1def949602ad9d7e2ee0edff4a7ec3c8e4be72'
SHORT_HEADER = """// side_lanes.S -- counted-X kernel (925 mean charge, 932 worst path), derived by
// measure/derive_k925.py from the validated k27r1 929 listing (K3: X-expensive counting).
// Interface, arithmetic, K1/K2/K3: proof.md s1. args[0..n-1] = {out, P, T, inv},
// args[n][0] = NULL ends the lane loop; out[0..11] = P-T reduced, out[12..14] = key."""
PRO = ['\tstr\txzr,[sp,#240]\t// X-count init (s24)']
CNT = ['\tldr\tx0,[sp,#240]\t// s24: count X>=p lanes', '\tadc\tx0,x0,xzr', '\tstr\tx0,[sp,#240]']
EPI = ['\tadrp\tx0,x_count_words\t// s24: node-local 128-bit X counter', '\tadd\tx0,x0,:lo12:x_count_words',
 '\tldr\tx1,[sp,#240]', '\tldr\tx2,[x0]', '\tldr\tx3,[x0,#8]', '\tadds\tx2,x2,x1', '\tadc\tx3,x3,xzr',
 '\tstr\tx2,[x0]', '\tstr\tx3,[x0,#8]']


def k3(lines):
    bnd = [n for n, l in enumerate(lines) if l.startswith('#ifndef __SIZEOF_POINTER__')]
    assert len(bnd) == 1 and bnd[0] > 0
    old = '\n'.join(lines[:bnd[0]])
    assert 'K1:' in old and 'K2:' in old and 'derive_k929' in old   # replaced header is the 929 one
    out, st = SHORT_HEADER.split('\n'), dict(pro=0, cnt=0, epi=0, com=0)
    for raw in lines[bnd[0]:]:
        s = raw.split('//')[0].strip()
        out.append(raw)
        if s == 'stp\tx9,x10,[x29,#160]':
            out += PRO; st['pro'] += 1
        elif s == 'b.lo\t.Lx_done':
            out += CNT; st['cnt'] += 1
        elif raw.strip() == '.Lexit:':
            out += EPI; st['epi'] += 1
        elif raw.strip() == '.section .rodata':
            out.insert(len(out) - 1, '.comm\tx_count_words,16,8'); st['com'] += 1
    assert st == dict(pro=1, cnt=1, epi=1, com=1), st
    return out


def derive(text):
    assert hashlib.sha256(text.encode()).hexdigest() == INPUT_SHA, 'not the validated 929 listing'
    assert text.count('ldp\tx17,x29,') == 3 and text.count('ldp\tx17,x2,[sp,#') == 6
    assert text.count('x2 carried b words') == 2 and 'adc\tx17,xzr,xzr' not in text and 'adc\tx17,x17,xzr' not in text
    assert '[sp,#240]' not in text and 'x_count_words' not in text
    return '\n'.join(k3(text.split('\n')))


if __name__ == '__main__':
    sys.stdout.write(derive(open(sys.argv[1]).read()))

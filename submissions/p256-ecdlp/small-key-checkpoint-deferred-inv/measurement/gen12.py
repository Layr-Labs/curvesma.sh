"""Deferred-inv revision: derive this kernel's sources from the predecessor's.

usage: gen12.py <predecessor measurement dir> <output dir>

inv's canonical finish becomes carry-only (meter 4 -> 5); its r3 >= p3 test
moves into the step tail as `ccmn x22,x28,#10,ne; b.hs Lz` after the
rare-event compare chain, the old slow block becomes Lz (`b.mi 2f`, meter
10 -> 11, back to the join), and the rare exit `2:` canonicalises inv.
Each edit must match exactly the stated number of times.
"""
import os
import re
import sys

src, dst = sys.argv[1], sys.argv[2]
EX = r"2:\nadds x15,x19,#1\nsbcs x16,x20,x28\nsbcs x17,x21,xzr\nsbcs x1,x22,x29\ncsel x19,x15,x19,hs\n"
EDITS = {
    "bloop_branchy.h": [
        ("(reversed D and t)", "(reversed D and t), then gen12.py (deferred inv test)", 1),
        (r'b.cc Lp2" S "_%=', r'b.cc Lq2" S "_%=', 1),
        (r'add x27,x27,#4\nb Lq2" S "_%=\nLp2" S "_%=:\ncmn x22,x28\nb.hs Ls2" S "_%=\n', r"add x27,x27,#5\n", 1),
        (r'"Ls2" S "_%=:\n', r'"Lz0" S "_%=:\nb.mi 2f\n', 1),
        (r'add x27,x27,#10\nb Lq2" S', r'add x27,x27,#11\nb Lq0" S', 1),
    ],
    "step.c": [
        ('#define BL_TAIL_P "',
         '#define BL_INV(S) "ccmn x22,x28,#10,ne\\nb.hs Lz0" S "_%=\\nLq0" S "_%=:\\n"\n#define BL_EXIT "'
         + EX + '" \\\n    "csel x20,x16,x20,hs\\ncsel x21,x17,x21,hs\\ncsel x22,x1,x22,hs\\nmov x0,#1\\n"\n'
         '#define BL_TAIL_P(S) "', 1),
        ('#define BL_TAIL_A "', '#define BL_TAIL_A(S) "', 1),
        ('",x6,#4,ne\\nb.eq 2f\\n" \\', '",x6,#4,ne\\n" BL_INV(S) \\', 2),
        ('BLOOP_SLOW_B "2:\\nmov x0,#1\\n3:', 'BLOOP_SLOW_B BL_EXIT "3:', 2),
    ],
}
for name in ["p256_field.h", "p256_sub.h", "p256_fused.h", "bloop_branchy.h", "floop.h", "step.c"]:
    text = open(os.path.join(src, name)).read()
    for old, new, n in EDITS.get(name, []):
        assert text.count(old) == n, (name, old)
        text = text.replace(old, new)
    if name == "bloop_branchy.h":
        text, n = re.subn(r'(BLOOP_FIELD_S\("(\w?)"\) T)', r'\1("\2")', text)
        assert n == 8
    open(os.path.join(dst, name), "w").write(text)
print("ok")

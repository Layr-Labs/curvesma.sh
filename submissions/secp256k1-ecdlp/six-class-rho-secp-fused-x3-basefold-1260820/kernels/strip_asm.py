"""Split clang's assembly of secp_walk.c into the two shipped files.

usage: python3 strip_asm.py secp_walk.s OUTDIR
secp_fe.s : fe_mul_f, fe_sqr_f, fe_inv   (field kernels and the inversion bound)
Only labels and instructions are kept; directives and comments are dropped. The counted batch loops
pass1 / pass2 are not compiled: secp_loops.s is emitted by gen_loops.py. pass1_c, pass2_c and the
test drivers are not shipped.
"""
import os, re, sys

GROUPS = {"secp_fe.s": ("fe_mul_f", "fe_sqr_f", "fe_inv")}
src = open(sys.argv[1], encoding="utf-8").read().split("\n")
funcs, cur = {}, None
for line in src:
    m = re.match(r"^([A-Za-z_][\w]*):", line)
    if m:
        cur = m.group(1)
        funcs[cur] = []
        continue
    if cur is None:
        continue
    if line.startswith(".Lfunc_end"):
        cur = None
        continue
    s = line.split("//")[0].rstrip()
    if not s.strip():
        continue
    if re.match(r"^\s+\.", s):
        continue
    funcs[cur].append(s)
for fname, names in GROUPS.items():
    with open(os.path.join(sys.argv[2], fname), "w", encoding="utf-8") as fh:
        for n in names:
            fh.write(n + ":\n")
            fh.write("\n".join(funcs[n]) + "\n")

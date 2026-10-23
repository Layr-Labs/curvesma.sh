"""Split clang's AArch64 output for p256_walk.c into the two shipped files.

  clang -target aarch64-linux-gnu -O3 -fno-tree-vectorize -mgeneral-regs-only \
        -ffreestanding -fno-builtin -fPIC -S -o p256_walk.s p256_walk.c
  python3 strip_asm.py p256_walk.s .

Comments, data and assembler directives are dropped; only the listed functions
and their local labels are kept (the counts in proof.md are computed from them).
"""
import re, sys

GROUPS = {"p256_field.s": ["fe_mul", "fe_add", "fe_sub", "fe_inv"],
          "p256_loops.s": ["pass1", "pass2", "pass1_dbl", "pass2_dbl"],
          "p256_slow.s": ["move_to_dbl", "move_back", "trav_update"]}
src, outdir = sys.argv[1], sys.argv[2]
for fname, keep in GROUPS.items():
    out = ["// clang -target aarch64-linux-gnu -O3 -fno-tree-vectorize -mgeneral-regs-only -ffreestanding -fno-builtin -fPIC -S p256_walk.c",
           "// Apple clang 21.0.0 (LLVM). Comments and assembler directives stripped; functions: " + ", ".join(keep)]
    cur, emit = None, False
    for raw in open(src):
        t = raw.rstrip("\n").split("//")[0].rstrip()
        if not t.strip():
            continue
        m = re.match(r"^([A-Za-z_.$][\w.$]*):\s*$", t.strip())
        if m:
            name = m.group(1)
            if not name.startswith(".L") and not name.startswith("%bb"):
                cur, emit = name, name in keep
                if emit:
                    out.append(name + ":")
            elif emit and name.startswith(".L"):
                out.append(name + ":")
            continue
        if t.strip().startswith("."):
            continue
        if emit and cur:
            out.append("\t" + re.sub(r"\s+", " ", t.strip(), count=1))
    with open(outdir + "/" + fname, "w") as fh:
        fh.write("\n".join(out) + "\n")

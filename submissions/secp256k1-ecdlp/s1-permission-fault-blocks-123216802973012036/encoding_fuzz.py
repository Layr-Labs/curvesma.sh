"""Assemble each kernel as a64.py expands it (kernels/leaf.s, route.s, secp_fe.s) with llvm-mc and check that the object's
.text holds exactly 4 bytes for every instruction a64_load counted (proof.md 3.1). FEAT_CMPBR's CBHS/CBHI register forms are
assembled as written (+cmpbr). Needs llvm-mc and llvm-size on PATH; experiment.py does not run it."""
import os, re, subprocess, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "kernels"))
from a64 import expand, a64_load
ATTRS = "-mattr=+lse,+rand,+ls64,+rcpc-immo,+cmpbr"
ok = True
with tempfile.TemporaryDirectory() as tmp:
    for name in ("leaf.s", "route.s", "secp_fe.s"):
        text = open(os.path.join(HERE, "kernels", name)).read()
        insns = len(a64_load(text)[0])
        src = os.path.join(tmp, name)
        open(src, "w").write("\t.text\n" + expand(text) + "\n")
        asm = subprocess.run(["llvm-mc", "-triple=aarch64", ATTRS, "-filetype=obj", src, "-o", src + ".o"], capture_output=True, text=True)
        size = subprocess.run(["llvm-size", "-A", src + ".o"], capture_output=True, text=True).stdout
        text_bytes = sum(int(n) for n in re.findall(r"^\.text\S*\s+(\d+)", size, re.M))
        good = asm.returncode == 0 and text_bytes == 4 * insns
        print(name, "insns", insns, "llvm-mc rc", asm.returncode, ".text bytes", text_bytes, "ok" if good else "FAIL", asm.stderr[:300])
        ok &= good
print("ALL_OK" if ok else "FAILED")
sys.exit(0 if ok else 1)

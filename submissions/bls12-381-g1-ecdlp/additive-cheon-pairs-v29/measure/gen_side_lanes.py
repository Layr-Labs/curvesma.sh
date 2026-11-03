#!/usr/bin/env python3
"""Generate side_lanes.S: one call processes n lanes of the pair side output.

usage: python3 gen_side_lanes.py <blst checkout at f7a88067> > side_lanes.S

Bodies of blst's __mul_mont_384, __sqr_384, __mul_by_1_mont_384 (commit
f7a88067, build/elf/mul_mont_384-armv8.S, Apache-2.0, Supranational) are
inlined with these mechanical edits, each asserted below:
  all:   trailing ret removed.
  MULNR: (mu) final conditional subtraction removed; the final carry-propagating
         adds write x11..x16 directly; dead top carry chain of round 6 removed.
  MUL:   (z) dead top carry chain of round 6 removed and the final conditional
         subtraction removed (z unreduced < 1.25p in x19..x24; the key step reduces).
  SQR:   the six result stores dropped; temporaries x5..x10 renamed x17, x0, x2, x3, x4, x29
         (low word stays in x17, the modulus stays in x5..x10; x29 is restored after the loop).
  MUL:   b = mu read sp-relative (no pointer register).
  both:  round 1's mov x17,xzr dropped; its one consumer adc x4,x17,xzr becomes adc x4,xzr,xzr.
  BY1:   the three input loads dropped; registers renamed so that it reads the square's low
         half in place and leaves the high half (x24..x28, x16) intact; n0 is x18.
The dead-carry deletions are exact when the 7-limb result is < 2^384, which
holds because the results are < 1.21p and < 1.25p (lemma in proof.md).
"""
import re, sys

src = open(sys.argv[1].rstrip("/") + "/build/elf/mul_mont_384-armv8.S").read()

def body(name):
    e = re.escape(name)
    m = re.search(r"^%s:\n(.*?)^\.size\t%s,\.-%s\n" % (e, e, e), src, re.S | re.M)
    b = m.group(1)
    assert not re.search(r"^\.L|^[A-Za-z_]\w*:", b, re.M)
    assert b.count("\tret\n") == 1 and b.rstrip().endswith("ret")
    b = b[:b.rindex("\tret\n")]
    assert "x18" not in b and "x30" not in b
    return b

def sub1(s, old, new):
    assert s.count(old) == 1, old
    return s.replace(old, new)

def dead_round6(b):
    head, tail = b.split("\tldr\tx17,[x2,8*5]\n")
    n = len(re.findall(r"^\tadc\tx17,.*\n", tail, re.M))
    assert n == 3, n
    tail = re.sub(r"^\tadc\tx17,.*\n", "", tail, flags=re.M)
    tail = sub1(tail, "\tadc\tx25,x17,xzr\n", "")
    last = tail[tail.rindex("x17"):].split("\n")[0]; assert tail[:tail.rindex("x17")].rsplit("\n",1)[1].startswith("\tumulh\tx3,x16") and last == "x17", last
    return head + "\tldr\tx17,[x2,8*5]\n" + tail

def pin_n0(b, pull):
    # n0 lives in x18 for the whole call: per-round reloads become dead
    n = b.count("\tldr\tx4,[x29,#12*__SIZEOF_POINTER__]\n"); assert n == 5
    b = b.replace("\tldr\tx4,[x29,#12*__SIZEOF_POINTER__]\n", "")
    assert b.count("\tmul\tx4,x4,x19\n") == 6
    b = b.replace("\tmul\tx4,x4,x19\n", "\tmul\tx4,x18,x19\n")
    b = sub1(b, "\tldp\tx4,x2,[x29,#12*__SIZEOF_POINTER__]\t// pull r_ptr\n", pull)
    assert not re.search(r"\bx4\b(?!,)", b.replace("mul\tx4,x18", "").replace("x4,x17,xzr", "").replace(",x4,", ",").replace(",x4\n", "\n").replace("x25,x4,", "")), "x4 n0 use left"
    return b

def round1_zero(b):
    # round 1 materialises a zero in x17 only to add it as a carry source: use xzr instead
    b = sub1(b, "\tmov\tx17,xzr\n", "")
    i = b.index("\tadc\tx4,x17,xzr\n")
    assert "x17" not in b[b.index("\tsubs\txzr,x19,#1"):i].replace("x17,xzr\n", "")
    return b[:i] + "\tadc\tx4,xzr,xzr\n" + b[i + len("\tadc\tx4,x17,xzr\n"):]

mul = dead_round6(round1_zero(pin_n0(body("__mul_mont_384"), "\t// pull r_ptr (n0 pinned in x18; out pointer loaded by the driver)\n")))
MUL = mul[:mul.index("\n\tsubs\tx26,x19,x5")] + "\n"     # z unreduced (< 1.25p) in x19..x24
for k in range(1, 6):                                              # b = mu read sp-relative
    MUL = sub1(MUL, "\tldr\tx17,[x2,8*%d]\n" % k, "\tldr\tx17,[sp,#%d]\n" % (96 + 8 * k))
assert not re.search(r"\b(x2|x29)\b", MUL)   # out (x2) and T (x29) live across it
cut = mul.index("\n\tsubs\tx26,x19,x5")
MULNR = mul[:cut] + "\n"
fin = ["\tadds\tx19,x20,x26\n", "\tadcs\tx20,x21,x27\n", "\tadcs\tx21,x22,x28\n",
       "\tadcs\tx22,x23,x0\n", "\tadcs\tx23,x24,x1\n", "\tadcs\tx24,x25,x3\n"]
pre, post = MULNR.split("loaded by the driver)\n")
for i, l in enumerate(fin):
    post = sub1(post, l, l.replace("x%d," % (19 + i), "x%d," % (11 + i), 1))
MULNR = pre + "loaded by the driver)\n" + post
SQR = body("__sqr_384")
st = re.findall(r"^\tstp\t.*\[x0.*\n", SQR, re.M); assert len(st) == 6
for l in st: SQR = sub1(SQR, l, "")
assert not re.search(r"\bx0\b|\bx17\b", SQR)
assert not re.search(r"\b(x2|x3|x4|x29)\b", SQR)
sren = {"x5": "x17", "x6": "x0", "x7": "x2", "x8": "x3", "x9": "x4", "x10": "x29"}
SQR = re.sub(r"\bx\d+\b", lambda m: sren.get(m.group(0), m.group(0)), SQR)   # p stays in x5..x10
BY1 = body("__mul_by_1_mont_384")
ld = re.findall(r"^\tldp\tx1[1-6],x1[1-6],\[x1.*\n", BY1, re.M); assert len(ld) == 3
for l in ld: BY1 = sub1(BY1, l, "")
assert not re.search(r"\bx0\b|\bx1\b|\bx3\b|\bx17\b|\bx18\b", BY1)
ren = {"x11": "x17", "x12": "x19", "x13": "x20", "x14": "x21", "x15": "x22", "x16": "x23",
       "x19": "x0", "x20": "x1", "x21": "x3", "x22": "x11", "x23": "x12", "x24": "x13",
       "x25": "x14", "x26": "x15", "x4": "x18"}
BY1 = re.sub(r"\bx\d+\b", lambda m: ren.get(m.group(0), m.group(0)), BY1)
# square in registers: low (x17,x19..x23) feeds BY1 directly; high (x24..x28,x16) is untouched by BY1
assert not re.search(r"\bx(16|2[4-8])\b", BY1)

drv = open(sys.argv[0].replace("gen_side_lanes.py", "side_lanes.S.in")).read()
# The block after the x' reduction (second product, y', key, loop tail) is emitted twice:
# once after the X >= 2p block (labels suffixed 2), once after the other x' blocks, so no
# x' path needs a jump to a shared continuation.
lo, hi = drv.index(".Lw_done:\n"), drv.index(".Lexit:\n")
tail = drv[lo + len(".Lw_done:\n"):hi]
labs = set(re.findall(r"^(\.L\w+):", tail, re.M))
assert labs == {".Ly_nonneg", ".Lk_e", ".Ly_store", ".Lk_d"}, labs
copy = tail
for l in labs: copy = re.sub(re.escape(l) + r"\b", l + "2", copy)
assert copy.rstrip().endswith("cbnz\tx0,.Llane")
drv = sub1(drv, "\tb\t.Lw_done\n", "\t// continuation, copy for this path\n" + copy + "\tb\t.Lexit\n")
for k, v in (("MULNR", MULNR), ("MUL", MUL), ("SQR", SQR), ("BY1", BY1)):
    m = "\t//@%s\n" % k
    assert drv.count(m) == (2 if k == "MUL" else 1), k
    drv = drv.replace(m, "\t// ---- inlined %s\n%s\t// ---- end %s\n" % (k, v, k))
print(drv)

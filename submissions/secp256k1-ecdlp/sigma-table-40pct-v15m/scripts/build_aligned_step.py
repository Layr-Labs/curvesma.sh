#!/usr/bin/env python3
"""build_aligned_step.py — build the v15l ALIGNED step listing and COUNT it.

The aligned rule (A') is: class key = min-x over the three x-values of the
negation+GLV orbit {x, beta*x, p-x-beta*x}, table row = 3 alignment slots
tbl[j][u] = [u*e_j]P0, and the slot is carried in the state record.  Nothing
here is estimated: the beta-multiply + reduce is SLICED out of the shipped
v15g listing (so it is the same schoolbook+fold idiom the 748 was priced on),
the block is inserted at the end of the priced loop body, and the loop head's
two index instructions (`and x22,x0,#2047`, `add x22,x25,x22,lsl #6`) are
replaced by ONE reload of the address computed at the end of the previous
iteration (`ldr x22,[x28,#72]`, initialised in the prologue, outside the
priced iteration).  count_iter.py then counts the merged listing.

Writes evidence/step_secp_v15l_aligned.s and evidence/counts-step-v15l.json.
"""
import hashlib, json, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
EV = os.path.join(PKG, "evidence")
ASM = os.path.join(EV, "step_secp_v15g.s")

L = open(ASM).read().split("\n")

# ---- slice one complete register-to-register mul+reduce from the shipped asm
start = next(i for i, l in enumerate(L) if l.strip() == "mul\tx0, x8, x4")
end = next(i for i, l in enumerate(L) if i > start and "// y1 again" in l)
MACRO = [m for m in L[start:end] if m.strip() and not m.strip().startswith("//")]
assert all("\t" in m for m in MACRO), "macro lines must be instructions"

# ---- the alignment block (beta*x3, third orbit value, argmin, 3-slot address)
blk = []
blk.append("\tmov\tx4, x0                          // multiplicand = x3 (new point)")
blk.append("\tmov\tx5, x1")
blk.append("\tmov\tx6, x2")
blk.append("\tmov\tx7, x3")
blk.append("\tldrb\tw22, [x28, #64]                 // u_in = slot carried from previous step")
blk += MACRO                                      # beta*x3 mod p -> x0..x3
blk += ["\tsubs\tx12, x29, x0                    // p - beta*x3  (x29 = -C = p limb0)",
        "\tsbcs\tx13, x8, x1                     // x8 = ~0 (p limbs 1,3)",
        "\tsbcs\tx14, x8, x2",
        "\tsbcs\tx15, x8, x3",
        "\tcsel\tx16, x24, xzr, cc               // borrow -> subtract C",
        "\tsubs\tx12, x12, x16",
        "\tsbcs\tx13, x13, xzr",
        "\tsbcs\tx14, x14, xzr",
        "\tsbcs\tx15, x15, xzr",
        "\tsubs\tx12, x12, x4                    // third = p - beta*x3 - x3",
        "\tsbcs\tx13, x13, x5",
        "\tsbcs\tx14, x14, x6",
        "\tsbcs\tx15, x15, x7",
        "\tcsel\tx16, x24, xzr, cc",
        "\tsubs\tx12, x12, x16",
        "\tsbcs\tx13, x13, xzr",
        "\tsbcs\tx14, x14, xzr",
        "\tsbcs\tx15, x15, xzr"]
blk += ["\tsubs\tx16, x4, x0                     // argmin over (x3, beta*x3, third)",
        "\tsbcs\tx17, x5, x1",
        "\tsbcs\tx19, x6, x2",
        "\tsbcs\tx20, x7, x3",
        "\tcsel\tx4, x4, x0, cc                  // m = min(c0,c1)",
        "\tcsel\tx5, x5, x1, cc",
        "\tcsel\tx6, x6, x2, cc",
        "\tcsel\tx7, x7, x3, cc",
        "\tcsel\tx21, xzr, x22, cc                // idx = 0 if c0<=c1 else u_in",
        "\tsubs\tx16, x4, x12",
        "\tsbcs\tx17, x5, x13",
        "\tsbcs\tx19, x6, x14",
        "\tsbcs\tx20, x7, x15",
        "\tcsel\tx4, x4, x12, cc                  // min x -> m",
        "\tcsel\tx5, x5, x13, cc",
        "\tcsel\tx6, x6, x14, cc",
        "\tcsel\tx7, x7, x15, cc",
        "\tmov\tx22, #2",
        "\tcsel\tx22, x21, x22, cc                // slot = argmin index",
        "\tand\tx21, x4, #2047                  // j = (min x)[0] & 2047",
        "\tadd\tx21, x21, x21, lsl #1           // 3j",
        "\tadd\tx21, x21, x22                    // 3j + slot",
        "\tadd\tx22, x25, x21, lsl #5            // &TBL[j][slot]  (32 B per point)",
        "\tstrb\tw22, [x28, #64]                 // persist slot",
        "\tstr\tx22, [x28, #72]                 // persist &TBL[j][slot] for the next head"]
NB = len(blk)

# ---- merge: head's two index instructions -> one reload; block at loop end
ih = next(i for i, l in enumerate(L) if l.startswith("\tand\tx22, x0, #2047"))
assert "lsl #6" in L[ih + 1], L[ih + 1]
head = ["\tldr\tx22, [x28, #72]                 // &TBL[j][slot] from previous step",
        "\t// (superseded: and x22,x0,#2047 ; add x22,x25,x22,lsl #6)"]
out = L[:ih] + head + L[ih + 2:]
ie = next(i for i, l in enumerate(out) if l.startswith("\tb\t.Ltail"))
# prologue init of the two carried bytes (outside the priced iteration)
ip = next(i for i, l in enumerate(out) if l.startswith("\tb\t.Lhead"))
init = ["\tstrb\twzr, [x28, #64]                 // slot init (prologue, not priced)",
        "\tand\tx22, x0, #2047                   // first address (prologue, not priced)",
        "\tadd\tx22, x25, x22, lsl #5"]
out = out[:ie] + blk + out[ie:]
# the prologue init needs x0 = first point x; v15g has it live at .Lhead, so place
# the init immediately before the prologue's `b .Lhead` using the loaded x0.
ip = next(i for i, l in enumerate(out) if l.startswith("\tb\t.Lhead"))
out = out[:ip] + init + out[ip:]

HDR = ("// step_secp_v15l_aligned.s -- GENERATED by scripts/build_aligned_step.py from\n"
       "// evidence/step_secp_v15g.s (the 748-Ir shipped listing).  Aligned A' rule: min-x\n"
       "// class key over {x, beta*x, p-x-beta*x} + 3-slot orbit table tbl[j][u]=[u*e_j]P0.\n"
       "// The beta* multiply+reduce above is SLICED from the shipped listing (lines %d-%d),\n"
       "// so the multiplicand kernel priced here is the one the 748 was counted on.\n"
       "// Head: the two superseded index instructions are replaced by one reload of the\n"
       "// address computed at the end of the previous iteration (prologue init, unpriced).\n"
       % (start + 1, end))
open(os.path.join(EV, "step_secp_v15l_aligned.s"), "w").write(HDR + "\n".join(out) + "\n")

def count(path):
    r = subprocess.run([sys.executable, os.path.join(HERE, "count_iter.py"), path],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    m = re.search(r"longest_iteration\s+(\d+)", r.stdout)
    assert m, r.stdout
    return int(m.group(1)), r.stdout

base_ir, base_txt = count(ASM)
al_ir, al_txt = count(os.path.join(EV, "step_secp_v15l_aligned.s"))
doc = dict(kind="counts-step-v15l",
           base_listing="evidence/step_secp_v15g.s", base_ir=base_ir,
           aligned_listing="evidence/step_secp_v15l_aligned.s", aligned_ir=al_ir,
           alignment_block_ir=NB,
           superseded_insns=2, added_head_reload=1,
           identity=f"{base_ir} - 2 + 1 + {NB} = {base_ir - 1 + NB}",
           identity_holds=(base_ir - 2 + 1 + NB == al_ir),
           mul_macro_ir=len(MACRO),
           note="aligned_ir is the counted longest iteration of the merged listing; the "
                "record/carried-slot init sits in the prologue, outside the priced iteration",
           sha256={f: hashlib.sha256(open(os.path.join(PKG, f), "rb").read()).hexdigest()
                   for f in ("evidence/step_secp_v15g.s",
                             "evidence/step_secp_v15l_aligned.s")},
           base_stdout=base_txt.strip().splitlines()[-1],
           aligned_stdout=al_txt.strip().splitlines()[-1])
json.dump(doc, open(os.path.join(EV, "counts-step-v15l.json"), "w"), indent=1)
print(json.dumps({k: doc[k] for k in ("base_ir", "aligned_ir", "alignment_block_ir",
                                      "identity", "identity_holds")}, indent=1))

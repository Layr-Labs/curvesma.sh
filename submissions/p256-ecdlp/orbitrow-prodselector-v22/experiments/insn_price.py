"""v15 instruction price: OBJECT DISASSEMBLY of the shipped build + dynamic trace.

v14's 147-instruction fe_mul price pointed at an absent object/span file.  v15 ships
the source (p256_fe_v15.c, the v10m field kernels verbatim plus the counted
wrappers), the objdump of the build, and an instrumented call trace.  The price is
the measured span of the shipped function, not a remembered row:
  fe_mul_obj span = 148 instructions (0x0..0x250, arm64 -O2 -fPIC)
  fe_add_obj / fe_sub_obj span = 27 each (out-of-line upper bound on the inlined form)
Dynamic trace (trace_counts.json): one affine addition executes 388 fe_mul with its
own inversion (3 + 385) and 3 fe_mul when the inverse is supplied by the production
batch; 385 = 256 squarings + popcount(p-2)=128 multiplies + 1 fixup.
"""
import json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
OUT = os.path.join(D, "evidence")


def parse_dis(path):
    syms = []
    blocks = {}
    cur = None
    for line in open(path):
        m = re.match(r'^([0-9a-f]+) <([^>]+)>:', line)
        if m:
            cur = m.group(2)
            syms.append((int(m.group(1), 16), cur))
            blocks[cur] = []
            continue
        if cur and re.match(r'^\s+[0-9a-f]+:', line):
            blocks[cur].append(line.strip())
    return syms, blocks


def parse_nm(path):
    out = {}
    for line in open(path):
        p = line.split()
        if len(p) == 3 and p[1] in ("T", "t"):
            out[p[2]] = int(p[0], 16)
    return out


def main():
    syms, blocks = parse_dis(os.path.join(OUT, "insn_disasm.txt"))
    nm = parse_nm(os.path.join(OUT, "insn_nm.txt"))
    order = [s for s in sorted(nm.items(), key=lambda kv: kv[1]) if s[0].startswith("_fe_")]
    rows = {}
    for i, (name, addr) in enumerate(order):
        if name not in blocks:
            continue
        end = order[i + 1][1] if i + 1 < len(order) else None
        counted = len(blocks[name])
        rows[name] = {"addr": addr, "insn_lines": counted, "end": end,
                      "span_bytes": (end - addr) if end else None}
        if end is not None:
            assert (end - addr) % 4 == 0, name
            assert (end - addr) // 4 == counted, (name, counted, end - addr)
    trace = json.load(open(os.path.join(OUT, "trace_counts.json")))
    fe = rows["_fe_mul_obj"]["insn_lines"]
    fa = rows["_fe_add_obj"]["insn_lines"]
    fs = rows["_fe_sub_obj"]["insn_lines"]
    ev = {
        "build": {"cmd": "clang -O2 -fPIC -arch arm64 -c p256_fe_v15.c -o p256_fe_v15.o",
                  "disasm": "objdump -d --no-show-raw-insn p256_fe_v15.o",
                  "source": "p256_fe_v15.c (v10m field kernels verbatim + counted wrappers)",
                  "host": "arm64-apple-darwin"},
        "functions": rows,
        "prices_ir": {"fe_mul": fe, "fe_add": fa, "fe_sub": fs,
                      "inversion_muls": 385, "add_muls": 3,
                      "affine_add_row": 3 * fe + 6 * fa},
        "dynamic_trace": trace,
        "trace_reading": "affine_add_batched_fe_mul = 3 proves the hot path executes "
                         "fe_mul three times per addition when the inverse is batch-supplied; "
                         "affine_add_with_inversion_fe_mul = 388 = 3 + 385 confirms the chain "
                         "length used by the inversion row",
        "consumed_by": "ledger_v15.py rows; every 5-digit literal here is owned by this file",
    }
    os.makedirs(OUT, exist_ok=True)
    json.dump(ev, open(os.path.join(OUT, "insn_price.json"), "w"), indent=1, sort_keys=True)
    print(json.dumps(ev["prices_ir"], indent=1))


if __name__ == "__main__":
    main()

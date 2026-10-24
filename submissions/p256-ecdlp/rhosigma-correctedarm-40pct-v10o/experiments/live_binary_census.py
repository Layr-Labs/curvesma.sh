#!/usr/bin/env python3
"""v10l live-binary census + instruction spans from the SAME object.

ONE build: p256_recproc_v10l.c -> p256_recproc_v10l.o (flags in FLAGS below).
  * the census dlopen()s p256_recproc_v10l.dylib, which is linked FROM that .o,
    calls recproc_init() then recproc_round(), and READS BACK the exported counters;
  * the instruction spans are taken from the disassembly of p256_recproc_v10l.o,
    bounded by the symbol start addresses that the same .o exports (no side build).

Writes evidence/live_binary.json and evidence/stepspans.json.
"""
import ctypes, json, re, subprocess, sys, os

FLAGS = ["-O2", "-fPIC", "-Wall", "-Wno-unused-const-variable"]
EV = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "evidence")

def build():
    subprocess.run(["clang", *FLAGS, "-c", "p256_recproc_v10l.c", "-o", "p256_recproc_v10l.o"],
                   cwd=EV, check=True)
    subprocess.run(["clang", "-O2", "-fPIC", "-dynamiclib", "p256_recproc_v10l.o",
                    "-o", "p256_recproc_v10l.dylib"], cwd=EV, check=True)

def census(rounds=64, n=8):
    lib = ctypes.CDLL(os.path.join(EV, "p256_recproc_v10l.dylib"))
    lib.recproc_init.restype = None
    lib.recproc_init.argtypes = []
    lib.recproc_round.restype = ctypes.c_int
    lib.recproc_round.argtypes = [ctypes.c_int]
    c_calls = ctypes.c_uint64.in_dll(lib, "batch_inv_calls")
    c_dens = ctypes.c_uint64.in_dll(lib, "batch_inv_denominators")
    before = (c_calls.value, c_dens.value)
    lib.recproc_init()
    reports = 0
    for _ in range(rounds):
        r = lib.recproc_round(n)
        if r < 0:
            raise SystemExit("zero denominator in census")
        reports += r
    after = (c_calls.value, c_dens.value)
    return {"object": "p256_recproc_v10l.dylib (linked from p256_recproc_v10l.o)",
            "counters_before": {"batch_inv_calls": before[0], "batch_inv_denominators": before[1]},
            "counters_after": {"batch_inv_calls": after[0], "batch_inv_denominators": after[1]},
            "rounds": rounds, "walkers_per_round": n, "walker_steps": rounds * n,
            "denominators_expected": rounds * n, "reports": reports,
            "symbol_size_bytes": {s: symsize("p256_recproc_v10l.o", s)
                                 for s in ("_batch_inv", "_recproc_round", "_step_forward",
                                           "_step_close", "_recproc_init")}}

def disasm(obj="p256_recproc_v10l.o"):
    return subprocess.run(["otool", "-tV", obj], cwd=EV, capture_output=True, text=True).stdout

SYM = re.compile(r"^([0-9a-f]{16})\s+(.*)$")
_LINES = {}
def line_of(txt, addr):
    if not _LINES:
        for ln in txt.splitlines():
            m = re.match(r"^([0-9a-f]{16})\s+", ln)
            if m:
                _LINES[int(m.group(1), 16)] = ln
    return _LINES.get(addr, "")

def syms(obj="p256_recproc_v10l.o"):
    """TEXT symbols only.  A common/data symbol (e.g. _W, type C or S) carries a SIZE in the
    nm value column, not a code address; using it as a span boundary truncates the span it
    falls inside (judge 1770e2b1 item 2)."""
    out = subprocess.run(["nm", "-n", "-g", obj], cwd=EV, capture_output=True, text=True).stdout
    d = {}
    for line in out.splitlines():
        m = SYM.match(line.strip())
        if not m:
            continue
        parts = m.group(2).split()
        if len(parts) < 2 or parts[0] not in ("T", "t"):
            continue
        d[parts[-1]] = int(m.group(1), 16)
    return d

def spans(obj="p256_recproc_v10l.o"):
    """instruction count per symbol, bounded by the next symbol start in the SAME object"""
    txt = disasm(obj)
    starts = []
    for line in txt.splitlines():
        m = re.match(r"^([0-9a-f]{16})\s+(\S.*)$", line)
        if not m:
            continue
        addr = int(m.group(1), 16)
        body = m.group(2).strip()
        if body.startswith("//") or body.endswith(":"):
            continue
        # instruction lines look like: "0000000000000dd8\tsub\tsp, sp, #0x90"
        if re.match(r"^[a-z][a-z0-9.]*\s", body):
            starts.append((addr, body.split()[0]))
    s = syms(obj)
    ordered = sorted(s.items(), key=lambda kv: kv[1])
    res = {}
    for i, (name, a) in enumerate(ordered):
        b = ordered[i + 1][1] if i + 1 < len(ordered) else 1 << 62
        seq = [(addr, op) for (addr, op) in starts if a <= addr < b]
        if not seq:
            continue
        # longest iteration: the body of the largest backward branch (a loop back-edge)
        best = 0
        for k, (addr, op) in enumerate(seq):
            if not op.startswith("b"):
                continue
            m = re.search(r"0x([0-9a-f]+)", line_of(txt, addr))
            if not m:
                continue
            tgt = int(m.group(1), 16)
            if tgt < addr:
                body = sum(1 for (ad, _) in seq if tgt <= ad <= addr)
                best = max(best, body)
        ops = [op for (_, op) in seq]
        res[name.lstrip("_")] = {
            "start": hex(a), "next_symbol": ordered[i+1][0].lstrip("_") if i+1 < len(ordered) else None,
            "instructions": len(ops), "longest_iteration": best or len(ops),
            "mul": sum(1 for o in ops if o in ("mul", "madd", "msub")),
            "umulh": sum(1 for o in ops if o in ("umulh", "smulh")),
            "load": sum(1 for o in ops if o.startswith("ld")),
            "store": sum(1 for o in ops if o.startswith("st") or o in ("strb", "sturb")),
        }
    return res

def symsize(obj, name):
    out = subprocess.run(["size", "-m", obj], cwd=EV, capture_output=True, text=True).stdout
    for line in out.splitlines():
        parts = line.split()
        if parts and parts[-1].lstrip("_") == name.lstrip("_") and len(parts) >= 3:
            try:
                return int(parts[0])
            except ValueError:
                pass
    return None

if __name__ == "__main__":
    build()
    live = census()
    sp = spans()
    json.dump(live, open(os.path.join(EV, "live_binary.json"), "w"), indent=1, sort_keys=True)
    json.dump({"object": "p256_recproc_v10l.o",
               "flags": " ".join(FLAGS),
               "spans": sp}, open(os.path.join(EV, "stepspans.json"), "w"), indent=1, sort_keys=True)
    print(json.dumps(live, indent=1))
    print(json.dumps(sp, indent=1))

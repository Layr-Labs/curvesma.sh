"""Regenerate the four straight-line kernels (mul.s, sqr.s from gen.py; lazy.s, root.s from gen_lazy.py) in a
temporary directory, check that every immediate of every emitted instruction encodes in AArch64, and execute a
sample in emu.py against big-integer references (check.py, check_lazy.py).

Asserts each generator's emitted count equals the count the package charges (mul 12,735, sqr 6,786,
lazy 3,457, root 3,074), that every call executes exactly that many instructions, and that every output matches.
The .s listings are not shipped (above the package byte cap): the generators emit them byte-for-byte.
Usage: python3 -I check_kernels.py
"""
import json, os, re, shutil, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
WANT = {"mul": 12735, "sqr": 6786}
WANT_LAZY = {"lazy": 3457, "root": 3074}
MEM = re.compile(r"(ldp|stp|ldr|str)\s+(\w+),.*\[(\w+)(?:,\s*#(-?\d+))?\]$")
ARITH = re.compile(r"(add|adds|sub|subs|cmp|cmn)\s+.*#(-?\d+)(?:,\s*lsl\s*#(\d+))?$")
MOV = re.compile(r"mov\s+\w+,\s*#(-?\d+)$")
EXTR = re.compile(r"extr\s+x\w+,\s*x\w+,\s*x\w+,\s*#(\d+)$")


def encodable(ins):
  """None if the instruction has no immediate; else whether the immediate encodes (64-bit registers)."""
  m = MEM.match(ins)
  if m:
    off = int(m.group(4) or 0)
    if m.group(1) in ("ldp", "stp"):
      return -512 <= off <= 504 and off % 8 == 0
    return 0 <= off <= 32760 and off % 8 == 0
  m = ARITH.match(ins)
  if m:
    return 0 <= int(m.group(2)) <= 4095 and m.group(3) in (None, "0", "12")
  m = MOV.match(ins)
  if m:
    v = int(m.group(1)) % (1 << 64)
    nv = (~v) % (1 << 64)
    return any((x >> s) << s == x and x >> s < 1 << 16 for x in (v, nv) for s in (0, 16, 32, 48))
  m = EXTR.match(ins)
  if m:
    return 0 <= int(m.group(1)) <= 63
  assert "#" not in ins, ins
  return None


def immediates(d):
  """Every emitted instruction of the four listings: count of immediates checked per listing; fails on any."""
  out = {}
  for name in ("mul", "sqr", "lazy", "root"):
    n = 0
    for line in open(os.path.join(d, name + ".s")):
      ins = line.split("//")[0].strip()
      if not ins or ins.endswith(":") or ins.startswith("."):
        continue
      ok = encodable(ins)
      assert ok is not False, (name, ins)
      n += ok is True
    out[name] = n
  return out


def run(argv, cwd):
  p = subprocess.run([sys.executable, "-I"] + argv, cwd=cwd, capture_output=True, text=True)
  if p.returncode:
    sys.stdout.write(p.stdout)
    sys.stderr.write(p.stderr)
    raise SystemExit("%s failed with %d" % (argv[0], p.returncode))
  return p.stdout


def main():
  d = tempfile.mkdtemp(prefix="kernels-")
  try:
    for f in ("gen.py", "model.py", "emu.py", "check.py", "gen_lazy.py", "check_lazy.py"):
      shutil.copy(os.path.join(HERE, f), os.path.join(d, f))
    emitted = {}
    for line in (run(["gen.py"], d) + run(["gen_lazy.py"], d)).splitlines():
      parts = line.split()
      if len(parts) >= 2 and parts[0] in {**WANT, **WANT_LAZY}:
        emitted[parts[0]] = int(parts[1])
    assert emitted == {**WANT, **WANT_LAZY}, emitted
    imms = immediates(d)
    checks = {}
    for kind, argv in (("mul", ["check.py", "mul", "40", "6", "11"]),
             ("sqr", ["check.py", "sqr", "20", "0", "12"])):
      out = run(argv, d)
      row = next(json.loads(line) for line in out.splitlines() if line.startswith("{"))
      assert row["executed_counts"] == [WANT[kind]], row
      assert row["fails"] == 0, row
      assert out.strip().splitlines()[-1] == "PASS", out
      checks[kind] = {"cases": row["cases"], "fails": row["fails"],
              "executed_counts": row["executed_counts"]}
    out = run(["check_lazy.py"] + sys.argv[1:3] if len(sys.argv) > 2 else ["check_lazy.py", "8", "11"], d)
    row = next(json.loads(line) for line in out.splitlines() if line.startswith("{"))
    assert row["executed_counts"] == {k: [v] for k, v in WANT_LAZY.items()}, row
    assert not any(row["fails"].values()) and out.strip().splitlines()[-1] == "PASS", out
    checks["lazy+root"] = row
    print(json.dumps({"emitted_executed": emitted, "immediates_checked_all_encodable": imms, "checks": checks,
             "PASS": True}))
  finally:
    shutil.rmtree(d, ignore_errors=True)


if __name__ == "__main__":
  main()

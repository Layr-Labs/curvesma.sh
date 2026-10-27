"""Regenerate the four straight-line kernels (mul.s, sqr.s from gen.py; lazy.s, root.s from gen_lazy.py) in a
temporary directory and execute a reduced sample in emu.py against big-integer references (check.py,
check_lazy.py).

Asserts each generator's emitted count equals the count the package charges (mul 13,761, sqr 6,777,
lazy 4,964, root 1,561), that every call executes exactly that many instructions, and that every output matches.
The .s listings are not shipped (above the package byte cap): the generators emit them byte-for-byte.
Usage: python3 -I check_kernels.py
"""
import json, os, shutil, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
WANT = {"mul": 13761, "sqr": 6777}
WANT_LAZY = {"lazy": 4964, "root": 1561}


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
        gen = run(["gen.py"], d)
        emitted = {}
        for line in gen.splitlines():
            parts = line.split()
            if len(parts) >= 2 and parts[0] in WANT:
                emitted[parts[0]] = int(parts[1])
        assert emitted == WANT, (emitted, gen)
        for line in run(["gen_lazy.py"], d).splitlines():
            parts = line.split()
            if len(parts) >= 2 and parts[0] in WANT_LAZY:
                emitted[parts[0]] = int(parts[1])
        assert {k: emitted[k] for k in WANT_LAZY} == WANT_LAZY, emitted
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
        out = run(["check_lazy.py", "8", "11"], d)
        row = next(json.loads(line) for line in out.splitlines() if line.startswith("{"))
        assert row["executed_counts"] == {k: [v] for k, v in WANT_LAZY.items()}, row
        assert not any(row["fails"].values()) and out.strip().splitlines()[-1] == "PASS", out
        checks["lazy+root"] = row
        print(json.dumps({"emitted_executed": emitted, "checks": checks, "PASS": True}))
    finally:
        shutil.rmtree(d, ignore_errors=True)


if __name__ == "__main__":
    main()

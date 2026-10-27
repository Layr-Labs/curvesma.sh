"""Regenerate the two straight-line kernels from gen.py in a temporary directory and execute a reduced
sample in emu.py against check.py's big-integer reference.

Asserts the generator's emitted count equals the count the package charges (mul 14,410, sqr 6,773),
that every call executes exactly that many instructions, and that every output limb matches.
The full .s listings are not shipped (mul.s is 331 KB, sqr.s is 155 KB, above the package byte cap):
gen.py emits them byte-for-byte and check.py verifies the generated files.  Usage: python3 -I check_kernels.py
"""
import json, os, shutil, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
WANT = {"mul": 14410, "sqr": 6773}


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
        for f in ("gen.py", "model.py", "emu.py", "check.py"):
            shutil.copy(os.path.join(HERE, f), os.path.join(d, f))
        gen = run(["gen.py"], d)
        emitted = {}
        for line in gen.splitlines():
            parts = line.split()
            if len(parts) >= 2 and parts[0] in WANT:
                emitted[parts[0]] = int(parts[1])
        assert emitted == WANT, (emitted, gen)
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
        print(json.dumps({"emitted_executed": emitted, "checks": checks, "PASS": True}))
    finally:
        shutil.rmtree(d, ignore_errors=True)


if __name__ == "__main__":
    main()

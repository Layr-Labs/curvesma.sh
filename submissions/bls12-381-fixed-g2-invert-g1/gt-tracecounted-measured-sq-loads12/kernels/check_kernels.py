"""Regenerate both kernels from gen.py in a temp dir; execute a reduced sample in emu.py.
Asserts emitted == charged counts (mul 14,446, sqr 6,821), exact executed counts, 0 mismatches.
Listings not shipped (341 KB + 164 KB). Usage: python3 -I check_kernels.py
"""
import json, os, shutil, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
WANT = {"mul": 14446, "sqr": 6821}


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

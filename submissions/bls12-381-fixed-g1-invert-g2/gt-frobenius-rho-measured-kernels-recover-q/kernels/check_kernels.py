"""Regenerate both kernels from gen.py in a temporary directory and run a sample in emu.py against
check.py's big-integer reference: emitted count = charged count (mul 12,735, sqr 6,786) = executed count
on every call, every output limb equal. The .s listings (295 KB and 159 KB) are not shipped; gen.py emits
them. Usage: python3 -I check_kernels.py
"""
import contextlib, io, json, os, random, runpy, shutil, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
WANT = {"mul": 12735, "sqr": 6786}

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

def small():
	"""The smaller in-process check experiment.py runs: one cyclotomic element, 5 multiplications
    (one on random inputs) and 4 squarings, the same assertions."""
	d, path = tempfile.mkdtemp(), list(sys.path)
	try:
		for f in ("gen.py", "model.py", "emu.py", "check.py"):
			shutil.copy(os.path.join(HERE, f), d)
		with contextlib.redirect_stdout(io.StringIO()):
			got = {k: v[0] for k, v in runpy.run_path(d + "/gen.py")["main"]().items()}
		assert got == WANT, got
		c = runpy.run_path(d + "/check.py")
		P, mul, rnd, rng = c["P"], c["f12mul"], c["rand12"], random.Random(13)
		xs = [c["f12pow"](rnd(rng), (P ** 6 - 1) * (P ** 2 + 1))]
		for _ in range(4):
			xs.append(mul(xs[-1], xs[0]))
		n = {}
		for k, cases in (("mul", list(zip(xs, xs[1:])) + [(rnd(rng), rnd(rng))]), ("sqr", [(x, None) for x in xs[1:]])):
			lay, ops, pool = c["setup"](k)
			for x, m in cases:
				e, out = c["run_one"](ops, pool, lay, x, m, rng)
				assert e == WANT[k] and out == [v * c["RM"] % P for v in c["flat"](mul(x, m or x))], k
				n[k] = n.get(k, 0) + 1
		return {"emitted": got, "executed": WANT, "cases": n, "fails": 0}
	finally:
		sys.path[:] = path
		shutil.rmtree(d, ignore_errors=True)

if __name__ == "__main__":
	main()

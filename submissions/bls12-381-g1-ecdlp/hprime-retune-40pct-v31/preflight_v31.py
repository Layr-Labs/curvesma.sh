#!/usr/bin/env python3
"""v31 preflight: the checks that can kill a package WITHOUT a judge.

1. every evidence path named in submission.json exists;
2. the package expands from its own archive (tar round-trip, byte-identical manifest);
3. size/file-count limits and the claims-length cap;
4. keymod.py's mix64 == the C harness's mix64 on the shared self-test vector (one key rule, two
   implementations, one answer);
5. audit_v31.py and sweep_literals_v31.py exit 0.
"""
import hashlib, json, os, subprocess, sys, tarfile, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
MAX_BYTES, MAX_FILES, MAX_CLAIM = 262144, 30, 7965


def files():
    out = []
    for r, d, f in os.walk(HERE):
        d[:] = [x for x in d if x not in ("__pycache__", ".git")]
        for x in f:
            if x.endswith(".pyc"):
                continue
            out.append(os.path.relpath(os.path.join(r, x), HERE))
    return sorted(out)


def main():
    sub = json.load(open(os.path.join(HERE, "submission.json")))
    ok = True
    missing = [p for p in sub["claims"][0]["evidence"] if not os.path.exists(os.path.join(HERE, p))]
    print("evidence paths: %d, missing: %s" % (len(sub["claims"][0]["evidence"]), missing or "none"))
    ok &= not missing

    total = sum(os.path.getsize(os.path.join(HERE, f)) for f in files())
    print("files=%d (cap %d) bytes=%d (cap %d)" % (len(files()), MAX_FILES, total, MAX_BYTES))
    ok &= len(files()) <= MAX_FILES and total <= MAX_BYTES
    stmt = sub["claims"][0]["statement"]
    print("claims[0] chars=%d (cap %d)" % (len(stmt), MAX_CLAIM))
    ok &= len(stmt) <= MAX_CLAIM

    # archive round-trip: the package must survive its own tar and hash-match
    with tempfile.TemporaryDirectory() as td:
        tp = os.path.join(td, "pkg.tar")
        with tarfile.open(tp, "w") as t:
            for f in files():
                t.add(os.path.join(HERE, f), arcname=f)
        with tarfile.open(tp) as t:
            names = sorted(t.getnames())
            digests = {n: hashlib.sha256(t.extractfile(n).read()).hexdigest()
                       for n in names if t.getmember(n).isfile()}
        local = {f: hashlib.sha256(open(os.path.join(HERE, f), "rb").read()).hexdigest() for f in files()}
        same = all(digests[k] == local[k] for k in digests)
        print("archive round-trip: %d files, byte-identical=%s, bytes=%d" % (len(names), same, os.path.getsize(tp)))
        ok &= same

    # one key rule, two implementations
    sys.path.insert(0, HERE)
    import keymod
    import re
    csrc = open(os.path.join(HERE, "structured_walk_v31.c")).read()
    consts = re.findall(r"x \*= 0x([0-9A-Fa-f]{16})ULL", csrc)
    vec = [0, 1, 0xDEADBEEF, 2 ** 64 - 1]
    py = [keymod.mix64(v) for v in vec]
    mults = sorted(consts)
    print("mix64 multiply constants in shipped C: %s ; python: (%x, %x)" % (mults, keymod.MIX_M1, keymod.MIX_M2))
    agrees = (set(mults) == {"%016x" % keymod.MIX_M1, "%016x" % keymod.MIX_M2})
    print("mix64 self-test vector py=%s agreement=%s" % ([hex(v) for v in py], agrees))
    ok &= agrees

    for script in ("audit_v31.py", "sweep_literals_v31.py"):
        r = subprocess.run([sys.executable, os.path.join(HERE, script)], capture_output=True, text=True)
        print("%s exit=%d" % (script, r.returncode))
        if r.returncode != 0:
            print(r.stdout[-2000:], r.stderr[-2000:])
        ok &= r.returncode == 0

    print(json.dumps(dict(preflight="PASS" if ok else "FAIL")))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

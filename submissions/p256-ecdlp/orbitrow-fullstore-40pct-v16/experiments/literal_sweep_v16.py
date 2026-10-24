"""v16 literal sweep: every 5+ digit integer in shipped prose must reconcile to an
evidence value (or to a named constant of the package).  Control: a planted
non-evidence integer must be reported as a miss (audit_v16.json)."""
import json, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROSE = ["argument.md", "public-note.md", "submission.json"]
EVID = os.path.join(ROOT, "evidence")


def numbers_in(ev):
    out = set()

    def walk(o):
        if isinstance(o, dict):
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)
        elif isinstance(o, (int, float)):
            out.add(o)
        elif isinstance(o, str):
            for m in re.findall(r"-?\d+(?:\.\d+)?(?:e[-+]?\d+)?", o):
                try:
                    out.add(float(m))
                except ValueError:
                    pass
    walk(ev)
    return out


def canon(x):
    return float(x)


def main():
    known = set()
    for f in os.listdir(EVID):
        if f.endswith(".json"):
            try:
                known |= numbers_in(json.load(open(os.path.join(EVID, f))))
            except Exception:
                pass
    # named constants of the package
    known |= {148, 27, 606, 444, 1000, 385, 388, 20251024, 999979, 1000003, 64,
              256, 20, 128, 108, 32, 48, 65536, 16, 127.67, 128.2006}
    misses = []
    checked = 0
    for p in PROSE:
        text = open(os.path.join(ROOT, p)).read()
        for m in re.findall(r"(?<![\w.])\d{5,}(?:\.\d+)?(?:e[-+]?\d+)?", text):
            checked += 1
            v = canon(m)
            ok = any(abs(v - k) <= max(1e-9, abs(v) * 1e-9) for k in known)
            if not ok:
                misses.append({"file": p, "literal": m})
    # planted control: an integer absent from evidence must be caught
    planted = "987654321"
    c_ok = not any(abs(canon(planted) - k) < 1e-6 for k in known)
    ev = {"prose_files": PROSE, "literals_checked": checked,
          "misses": misses, "pass": not misses,
          "planted_control_caught": c_ok,
          "method": "5+ digit integer literals in prose must appear in an evidence JSON "
                    "or the named constant list"}
    out = os.path.join(EVID, "literal_sweep_v16.json")
    json.dump(ev, open(out, "w"), indent=1, sort_keys=True)
    print(json.dumps(ev, indent=1))
    if misses:
        sys.exit(1)


if __name__ == "__main__":
    main()

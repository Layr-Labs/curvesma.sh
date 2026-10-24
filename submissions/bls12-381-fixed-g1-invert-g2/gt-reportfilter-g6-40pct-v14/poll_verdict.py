#!/usr/bin/env python3
"""Poll one submission until it leaves review/validating, then write the verdict verbatim."""
import json, os, subprocess, sys, time

SID = sys.argv[1] if len(sys.argv) > 1 else "c92286c0-4586-4c95-a066-5e1639cebc5f"
OUT = sys.argv[2] if len(sys.argv) > 2 else "evidence/verdict-v13.md"
BENCH = "Holindauer/curvesmash/bls12-381-fixed-g1-invert-g2"
env = dict(os.environ); env.pop("YUKON_API_TOKEN", None)

deadline = time.time() + 60 * 14
last = None
while time.time() < deadline:
    r = subprocess.run(["yukon", "submissions", BENCH, "--all", "--json"],
                       capture_output=True, text=True, env=env, timeout=180)
    try:
        d = json.loads(r.stdout)
    except Exception:
        time.sleep(60); continue
    hit = [s for s in d["submissions"] if s["id"] == SID]
    if hit:
        s = hit[0]
        st = s.get("status"); sc = s.get("officialScore")
        if (st, sc) != last:
            print("%s status=%s score=%s" % (time.strftime("%H:%M:%SZ", time.gmtime()), st, sc), flush=True)
            last = (st, sc)
        if st not in ("review", "validating", "queued", "running"):
            reason = ""
            ev = s.get("evaluation") or {}
            reason = (ev.get("validity") or {}).get("reason") or s.get("rejectionReason") or ""
            body = ["# v13 verdict — submission %s" % SID, "",
                    "status: `%s`" % st, "officialScore: `%s`" % sc, "",
                    "The text below is the platform's `evaluation.validity.reason` field copied verbatim;",
                    "nothing is paraphrased and no number is retyped.", "", "---", "", reason, ""]
            os.makedirs(os.path.dirname(OUT), exist_ok=True)
            open(OUT, "w").write("\n".join(body))
            open(OUT.replace(".md", ".json"), "w").write(json.dumps(s, indent=1, sort_keys=True))
            print("VERDICT %s written to %s (%d chars)" % (st, OUT, len(reason)))
            sys.exit(0)
    time.sleep(60)
print("still in flight after the poll window")
sys.exit(3)

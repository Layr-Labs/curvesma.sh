#!/usr/bin/env python3
"""census_reduce_v15m.py — reduce the lane's raw v15m census to shipped evidence.

Input  : the lane census (census-v15m.json, 51,200 per-trial rows, ~13.5 MB) and
         the lane census script (census_v15m.py) for the planted-attribution
         selftest.  Neither fits the 262,144-byte archive cap, so the reduction is:
           * per-density EXACT INTEGER counts, re-summed here from the rows;
           * the dp1 rows themselves (the certified band) shipped zlib+base64,
             with the sha256 of the uncompressed JSON so the certificate can
             decompress, re-sum and byte-check them;
           * the sha256 + byte size of the raw census it came from.
Output : evidence/census-v15m-counts.json
         evidence/census-v15m-dp1-rows.json
         evidence/attribution-selftest-v15m.json

No number in the outputs is typed by hand: every count is a sum over the rows.
Usage: census_reduce_v15m.py <lane-dir> <package-dir>
"""
import base64, hashlib, json, os, subprocess, sys, zlib

lane = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else ".")
pkg = os.path.abspath(sys.argv[2] if len(sys.argv) > 2 else ".")
EV = os.path.join(pkg, "evidence")
os.makedirs(EV, exist_ok=True)

RAW = os.path.join(lane, "census-v15m.json")
raw_bytes = open(RAW, "rb").read()
cen = json.loads(raw_bytes)
rows = cen["rows"]
S = rows[0]["s_ops"]


def cp95(k, n):
    """one-sided Clopper-Pearson 95% lower bound, exact Decimal bisection."""
    from decimal import Decimal, getcontext
    getcontext().prec = 60
    k, n = int(k), int(n)
    if k == 0:
        return "0"
    lnC = (sum(Decimal(i).ln() for i in range(n - k + 1, n + 1))
           - sum(Decimal(i).ln() for i in range(1, k + 1)))

    def tail(p):
        p = Decimal(p)
        q = Decimal(1) - p
        t = (lnC + k * p.ln() + (n - k) * q.ln()).exp()
        terms, m = [t], t
        for j in range(k, n):
            t = t * (n - j) * p / ((j + 1) * q)
            terms.append(t)
            if t > m:
                m = t
        return m * sum(x / m for x in terms)

    lo, hi = Decimal("0.0"), Decimal("0.999999999")
    for _ in range(120):
        mid = (lo + hi) / 2
        if tail(mid) > Decimal("0.05"):
            hi = mid
        else:
            lo = mid
        if hi - lo < Decimal("1e-15"):
            break
    return str(lo)


dens = []
for dp in sorted({r["dp_bits"] for r in rows}):
    sub = [r for r in rows if r["dp_bits"] == dp]
    n = len(sub)
    ok = [r for r in sub if r["ok"]]
    strict = [r for r in ok if r["steps"] is not None and r["steps"] <= S]
    audit = [r for r in sub if r["audit_merge"]]
    drain = [r for r in ok if r["via"] == "drain"]
    drain_only = sum(1 for r in sub if r["drain_only"])
    dens.append(dict(
        dp_bits=dp, density="2^-%d" % dp, trials=n,
        verified_pass_bounded=len(ok), verified_strict=len(strict),
        audit_merge=len(audit), drain_success=len(drain),
        drain_only_pairs=drain_only,
        mean_records=round(sum(r["nreports"] for r in sub) / n, 3),
        p_uncond=round(len(ok) / n, 6),
        cp95_uncond=cp95(len(ok), n),
        cp95_strict=cp95(len(strict), n),
        cp95_audit=cp95(len(audit), n)))

doc = dict(
    kind="census-v15m-counts",
    statistic="UNCONDITIONAL verified recovery ([x]P0 == T) with no conditioning filter; "
              "the audit merge is published side by side, never used as a filter",
    budget_convention="a trial walks cap_mult * S_ops steps; the budget is tested once per "
                      "pass of ntrails walkers, so a recovery can land up to ntrails-1 steps "
                      "past S_ops (0.86%% of the toy budget). "
                      "'verified_strict' re-cuts the same rows at steps <= S_ops",
    channels=dict(budget="both records inside the search budget",
                  drain="one record delivered by the CHARGED drain flush",
                  drain_only="pair first meets inside the drain: counted, NEVER success"),
    audit="independent coverage audit: two distinct walkers in one class, in-budget",
    generator="census_reduce_v15m.py",
    lane=dict(census_script="census_v15m.py",
              raw_census_sha256=hashlib.sha256(raw_bytes).hexdigest(),
              raw_census_bytes=len(raw_bytes), raw_rows=len(rows),
              model_rate=cen["model_rate"], cap_mult=cen["cap_mult"],
              drain_spacings=cen["drain_spacings"], workers=cen["workers"],
              trials_per_density=cen["trials_per_density"],
              toy_curve=dict(name="q25", p_bits=25, R=2048, ntrails=16),
              s_ops_toy=S),
    densities=dens)

json.dump(doc, open(os.path.join(EV, "census-v15m-counts.json"), "w"),
          separators=(",", ":"))

# ---- dp1 rows, compact: [trial, ok, steps, audit_merge, drain_only, nreports]
dp1 = [r for r in rows if r["dp_bits"] == 1]
compact = [[r["trial"], r["ok"], r["steps"], r["audit_merge"],
            r["drain_only"], r["nreports"]] for r in dp1]
plain = json.dumps(compact, separators=(",", ":")).encode()
enc = base64.b64encode(zlib.compress(plain, 9)).decode()
json.dump(dict(
    kind="census-v15m-dp1-rows",
    note="exact per-trial rows of the CERTIFIED band (dp1). Fields: "
         "[trial, verified(0/1), steps_at_recovery_or_null, audit_merge(0/1), "
         "drain_only(0/1), reports_written]. Decompress with zlib+base64 and parse "
         "as JSON; the certificate re-sums every field and checks sha256(plain).",
    encoding="zlib+base64",
    fields=["trial", "verified", "steps", "audit_merge", "drain_only", "nreports"],
    trials=len(compact),
    sha256_plain=hashlib.sha256(plain).hexdigest(),
    plain_bytes=len(plain),
    s_ops_toy=S,
    generator="census_reduce_v15m.py",
    payload=enc),
    open(os.path.join(EV, "census-v15m-dp1-rows.json"), "w"), separators=(",", ":"))

# ---- planted attribution selftest, run now from the lane script
r = subprocess.run([sys.executable, os.path.join(lane, "census_v15m.py"), "--selftest"],
                   capture_output=True, text=True, cwd=lane,
                   env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"))
assert r.returncode == 0, r.stderr[-2000:]
st = json.loads(r.stdout[r.stdout.index("{"):])
json.dump(dict(kind="attribution-selftest-v15m",
               generator="census_v15m.py --selftest",
               checks="audit_key is the single derivation of a trial index; 40 planted "
                      "records are labelled through it and must come back with the "
                      "planted trial id and two distinct slots",
               planted=st["planted"], bad=st["bad"],
               end_to_end=dict(trials=len(st["trial_ids"]), trial_ids=st["trial_ids"],
                               verified=st["ok"], audit_merge=st["audit"],
                               note="rows come from a real 3-trial chunk at dp1; a verified "
                                    "recovery without an audit merge is legitimate (a "
                                    "same-walker self-collision across an epoch restart)"),
               verdict="PASS" if st["bad"] == 0 and st["planted"] == 40 else "FAIL"),
          open(os.path.join(EV, "attribution-selftest-v15m.json"), "w"), indent=1)

print(json.dumps({k: (v["verified_pass_bounded"], v["trials"], v["cp95_uncond"])
                  for k, v in ((d["density"], d) for d in dens)}, indent=1))
print("dp1 rows", len(compact), "plain", len(plain), "encoded", len(enc),
      "selftest bad", st["bad"])

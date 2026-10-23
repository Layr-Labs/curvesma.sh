"""Per-item instruction costs for the v3 kernel, from a measure.sh output directory.

  python3 -B summarize_v3.py <run-dir> > summary.json

Reference additions (ref_add, the walk's batched additions, and ref_proj_add, the table build)
are reference units and are removed from R. Everything else the kernel executes stays in R.
Per-call reference costs come from the union toggle (ref_add + ref_proj_add) with the call
counters printed by the kernel; the build-only ref_proj_add cost is the union minus the ref_add
part of the 512-step run.
"""
import json
import re
import sys

run = sys.argv[1]


def summary(name):
    for line in open("%s/%s" % (run, name)):
        if line.startswith("summary:"):
            return int(line.split()[1])
    raise SystemExit("no summary in " + name)


def counts(name):
    last = None
    for line in open("%s/%s" % (run, name)):
        if re.fullmatch(r"steps=\d+( \w+=\d+)*", line.strip()):
            last = line.strip()
    out = {}
    for tok in last.split():
        k, _, v = tok.partition("=")
        out[k] = int(v)
    return out


prod = {n: summary("step-%d.callgrind" % n) for n in (256, 512)}
ref = {n: summary("ref-%d.callgrind" % n) for n in (256, 512)}
cnt = {n: counts("step-%d.stderr" % n) for n in (256, 512)}
dp4 = {n: summary("mode-dp4-%d.callgrind" % n) for n in (256, 512)}
dp4c = counts("mode-dp4-512.stderr")
cap1 = {n: summary("mode-cap1-%d.callgrind" % n) for n in (256, 512)}
cap1c = counts("mode-cap1-512.stderr")
build = summary("mode-build-512.callgrind")

per_point_total = (prod[512] - prod[256]) / 256
adds_per_point = (cnt[512]["adds"] - cnt[256]["adds"]) / 256
ref_per_point = (ref[512] - ref[256]) / 256
ref_add_per_call = (ref[512] - ref[256]) / (cnt[512]["adds"] - cnt[256]["adds"])
R_point = per_point_total - ref_per_point
ref_proj_total = ref[512] - ref_add_per_call * cnt[512]["adds"]
R_build = build - ref_proj_total
dp_other = (dp4[512] - prod[512]) / dp4c["dps"]
restart_other = ((cap1[512] - prod[512]) - (cap1c["adds"] - cnt[512]["adds"]) * ref_add_per_call) / cap1c["restarts"]

out = {
    "run": run,
    "kernel": "cheon_step_v3.c (full walk step, 256 records, production tables)",
    "per_point_total_Ir": per_point_total,
    "adds_per_point_measured": adds_per_point,
    "ref_add_Ir_per_point": ref_per_point,
    "ref_add_Ir_per_call": ref_add_per_call,
    "R_point_instructions": R_point,
    "R_point_units": R_point / 2000,
    "dp_R_instructions_per_dp": dp_other,
    "restart_R_instructions_per_restart": restart_other,
    "restart_adds_per_restart": (cap1c["adds"] - cnt[512]["adds"]) / cap1c["restarts"],
    "build_total_Ir": build,
    "build_ref_proj_Ir": ref_proj_total,
    "build_proj_calls": cnt[512]["proj"],
    "R_build_instructions": R_build,
    "raw": {"prod": prod, "ref_union": ref, "dp4": dp4, "cap1": cap1, "build": build,
            "counts_prod": cnt, "counts_dp4_512": dp4c, "counts_cap1_512": cap1c},
}
print(json.dumps(out, indent=1))

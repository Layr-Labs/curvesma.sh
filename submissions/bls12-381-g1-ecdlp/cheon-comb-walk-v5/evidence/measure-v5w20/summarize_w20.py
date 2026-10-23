"""Worst-case per-item costs for the window-20 kernel.

  python3 -B summarize_w20.py <run1-dir> <run2-dir> > summary.json

run1: REF core (ref-4096/8192), natural production (step-4096/8192), forced and forcedcap1 (lite tables,
4096 and 8192 points), forcedfull (full tables, 4096 and 8192 points: cross-check of the lite tables).
run2: table build at 1024 and 2048 stages (build-only mode, lite tables).
"""
import json
import re
import sys

R1, R2 = sys.argv[1], sys.argv[2]


def summary(run, name):
    for line in open("%s/%s" % (run, name)):
        if line.startswith("summary:"):
            return int(line.split()[1])
    raise SystemExit("no summary in " + name)


def counts(run, name):
    last = None
    for line in open("%s/%s" % (run, name)):
        if re.match(r"^(steps|build done)", line.strip()):
            last = line.strip()
    return {k: int(v) for k, v in re.findall(r"(\w+)=(\d+)", last)}


ref = {n: summary(R1, "ref-%d.callgrind" % n) for n in (4096, 8192)}
nat = {n: counts(R1, "step-%d.stderr" % n) for n in (4096, 8192)}
core_adds = nat[8192]["adds"] - nat[4096]["adds"]
core_per_add = (ref[8192] - ref[4096]) / core_adds

f = {n: summary(R1, "mode-forced-%d.callgrind" % n) for n in (4096, 8192)}
fc = {n: summary(R1, "mode-forcedcap1-%d.callgrind" % n) for n in (4096, 8192)}
ff = {n: summary(R1, "mode-forcedfull-%d.callgrind" % n) for n in (4096, 8192)}
fcnt = {n: counts(R1, "mode-forced-%d.stderr" % n) for n in (4096, 8192)}
fccnt = {n: counts(R1, "mode-forcedcap1-%d.stderr" % n) for n in (4096, 8192)}
ffcnt = {n: counts(R1, "mode-forcedfull-%d.stderr" % n) for n in (4096, 8192)}
nat_round = summary(R1, "step-8192.callgrind") - summary(R1, "step-4096.callgrind")
f_round = f[8192] - f[4096]
ff_round = ff[8192] - ff[4096]
f_adds_round = fcnt[8192]["adds"] - fcnt[4096]["adds"]
ff_adds_round = ffcnt[8192]["adds"] - ffcnt[4096]["adds"]
forced_other_per_point = (f_round - f_adds_round * core_per_add) / 4096
natural_other_per_point = (nat_round - core_adds * core_per_add) / 4096
restart_other = ((fc[8192] - fc[4096]) - f_round) / (fccnt[8192]["restarts"] - fccnt[4096]["restarts"])
crosscheck_diff = ff_round - f_round

b1 = summary(R2, "mode-build-1024.callgrind")
b2 = summary(R2, "mode-build-2048.callgrind")
bc1 = counts(R2, "mode-build-1024.stderr")
bc2 = counts(R2, "mode-build-2048.stderr")
stage_adds = (bc2["adds"] - bc1["adds"]) / 1024
per_stage_total = (b2 - b1) / 1024
per_stage_other = per_stage_total - stage_adds * core_per_add
fixed_total = b1 - 1024 * per_stage_total
fixed_other = fixed_total - bc1["adds"] * core_per_add + 1024 * stage_adds * core_per_add - 1024 * per_stage_total + 1024 * per_stage_other + 0
# simpler: fixed other = build(1024 stages) other - 1024 * per-stage other
build1_other = b1 - bc1["adds"] * core_per_add
fixed_other = build1_other - 1024 * per_stage_other

out = {
    "window": 20,
    "windows": 13,
    "tables": "2 bases x 13 windows x 2^19 entries (full); lite runs use 4096 entries",
    "core_per_addition_Ir": core_per_add,
    "forced_adds_per_point": f_adds_round / 4096,
    "forced_other_per_point_Ir": forced_other_per_point,
    "natural_other_per_point_Ir": natural_other_per_point,
    "forced_round_total_Ir": f_round,
    "forcedfull_round_total_Ir": ff_round,
    "crosscheck_full_minus_lite_round_Ir": crosscheck_diff,
    "restart_other_per_restart_Ir": restart_other,
    "forced_lemma": fcnt[8192]["lemma"],
    "build_adds_per_stage": stage_adds,
    "build_other_per_stage_Ir": per_stage_other,
    "build_fixed_other_Ir": fixed_other,
    "forced_adds_per_round": f_adds_round,
}
print(json.dumps(out, indent=1))

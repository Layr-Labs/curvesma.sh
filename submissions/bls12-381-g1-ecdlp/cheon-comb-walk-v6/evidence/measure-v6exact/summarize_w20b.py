"""Exact per-item costs for the window-20 kernel (no averages of the reference core).

  python3 -B summarize_w20b.py <runA> <runB> <runC> > summary.json

runA: REF_FUNCS=core_mul (ref-4096/8192 totals and counts), natural production steps, forced, forcedcap1, forcedfull.
runB: REF_FUNCS=core_sqr (ref-4096/8192 totals and counts).
runC: build mode at 1024 and 2048 stages (lite tables).

Reference-core instructions are the measured wrapper totals of core_mul and core_sqr times their exact call counts (the
kernel counts every call). Other work is the remainder of each measured total. The per-call wrapper costs c_mul and c_sqr are
constant (constant-time field arithmetic), measured as the differenced total over the differenced call count of the same run.
"""
import json
import re
import sys
from fractions import Fraction as F

RA, RB, RC = sys.argv[1], sys.argv[2], sys.argv[3]


def summary(run, name):
    for line in open("%s/%s" % (run, name)):
        if line.startswith("summary:"):
            return int(line.split()[1])
    raise SystemExit("no summary in %s/%s" % (run, name))


def counts(run, name):
    last = None
    for line in open("%s/%s" % (run, name)):
        if re.match(r"^(steps|build done)", line.strip()):
            last = line.strip()
    return {k: int(v) for k, v in re.findall(r"(\w+)=(\d+)", last)}


refA = {n: summary(RA, "ref-%d.callgrind" % n) for n in (4096, 8192)}
cntA = {n: counts(RA, "ref-%d.stderr" % n) for n in (4096, 8192)}
refB = {n: summary(RB, "ref-%d.callgrind" % n) for n in (4096, 8192)}
cntB = {n: counts(RB, "ref-%d.stderr" % n) for n in (4096, 8192)}
c_mul = F(refA[8192] - refA[4096]) / (cntA[8192]["mul"] - cntA[4096]["mul"])
c_sqr = F(refB[8192] - refB[4096]) / (cntB[8192]["sqr"] - cntB[4096]["sqr"])

nat = {n: summary(RA, "step-%d.callgrind" % n) for n in (4096, 8192)}
natc = {n: counts(RA, "step-%d.stderr" % n) for n in (4096, 8192)}


def other_round(total_8192, total_4096, c8, c4):
    d_total = total_8192 - total_4096
    d_mul = c8["mul"] - c4["mul"]
    d_sqr = c8["sqr"] - c4["sqr"]
    core = d_mul * c_mul + d_sqr * c_sqr
    return d_total - core, d_mul, d_sqr


nat_other_round, nat_dm, nat_ds = other_round(nat[8192], nat[4096], natc[8192], natc[4096])
f = {n: summary(RA, "mode-forced-%d.callgrind" % n) for n in (4096, 8192)}
fc_ = {n: counts(RA, "mode-forced-%d.stderr" % n) for n in (4096, 8192)}
forced_other_round, f_dm, f_ds = other_round(f[8192], f[4096], fc_[8192], fc_[4096])
fc = {n: summary(RA, "mode-forcedcap1-%d.callgrind" % n) for n in (4096, 8192)}
fcc = {n: counts(RA, "mode-forcedcap1-%d.stderr" % n) for n in (4096, 8192)}
fcap_other_round, fcap_dm, fcap_ds = other_round(fc[8192], fc[4096], fcc[8192], fcc[4096])
restart_rounds = fcc[8192]["restarts"] - fcc[4096]["restarts"]
restart_other = F(fcap_other_round - forced_other_round, restart_rounds)
assert fcap_dm == f_dm and fcap_ds == f_ds, "restart path must not call the reference core"
ff = {n: summary(RA, "mode-forcedfull-%d.callgrind" % n) for n in (4096, 8192)}
ffc = {n: counts(RA, "mode-forcedfull-%d.stderr" % n) for n in (4096, 8192)}
ff_other_round, _, _ = other_round(ff[8192], ff[4096], ffc[8192], ffc[4096])

b = {n: summary(RC, "mode-build-%d.callgrind" % n) for n in (1024, 2048)}
bc = {n: counts(RC, "mode-build-%d.stderr" % n) for n in (1024, 2048)}
b1_other = F(b[1024]) - bc[1024]["mul"] * c_mul - bc[1024]["sqr"] * c_sqr
b2_other = F(b[2048]) - bc[2048]["mul"] * c_mul - bc[2048]["sqr"] * c_sqr
stage_mul = F(bc[2048]["mul"] - bc[1024]["mul"], 1024)
stage_sqr = F(bc[2048]["sqr"] - bc[1024]["sqr"], 1024)
stage_adds = F(bc[2048]["adds"] - bc[1024]["adds"], 1024)
per_stage_other = (b2_other - b1_other) / 1024
fixed_other = b1_other - 1024 * per_stage_other

out = {
    "window": 20,
    "c_mul_instructions_per_call": float(c_mul),
    "c_sqr_instructions_per_call": float(c_sqr),
    "c_mul_exact": str(c_mul),
    "c_sqr_exact": str(c_sqr),
    "walk_round_mul_calls": nat_dm, "walk_round_sqr_calls": nat_ds,
    "forced_round_mul_calls": f_dm, "forced_round_sqr_calls": f_ds,
    "forced_adds_per_point": F(f_ds, 4096).__float__(),
    "natural_other_per_point_Ir": float(F(nat_other_round, 4096)),
    "forced_other_per_point_Ir": float(F(forced_other_round, 4096)),
    "forced_other_per_point_exact": str(F(forced_other_round, 4096)),
    "forcedfull_minus_forced_round_other_Ir": float(ff_other_round - forced_other_round),
    "restart_other_per_restart_Ir": float(restart_other),
    "restart_other_per_restart_exact": str(restart_other),
    "build_stage_mul_calls": float(stage_mul),
    "build_stage_sqr_calls": float(stage_sqr),
    "build_adds_per_stage": float(stage_adds),
    "build_other_per_stage_Ir": float(per_stage_other),
    "build_other_per_stage_exact": str(per_stage_other),
    "build_fixed_other_Ir": float(fixed_other),
    "build_fixed_other_exact": str(fixed_other),
    "forced_lemma": fc_[8192]["lemma"],
}
print(json.dumps(out, indent=1))

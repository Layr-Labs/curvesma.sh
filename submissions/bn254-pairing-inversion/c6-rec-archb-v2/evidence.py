#!/usr/bin/env python3
"""
evidence.py — 100% self-contained pure-stdlib verification entry point (python -I -B).
Executes live:
  1. Big-int BN254 + C6-m184_w216-rec algebra & mod-r CRT factorization across all 6 roots;
  2. Exhaustive R^2=50 fundamental-domain bijection (107/107 orbits, 0 missed, 0 dups);
  3. Arch B-Hybrid 0-mismatch replay prototype + 64-bit scalar SWAR varying-T(a) identity;
  4. Self-contained Band-3 7D Fincke-Pohst coset-walk prototype (0 missed hits + nodes/leaves bounds);
  5. Exact 12D C6 special-q (2^65..2^125) & boot lattice det/membership + norm-envelope audit;
  6. Pure-Python hexagonal Z[omega] singleton + clique filtering simulator with exact kernel check;
  7. Live norm-histogram-weighted Euler-product Retention + Stage B + full de Bruijn yield;
  8. 4-parameter deterministic hard-cap LA + live 62-row LB=59 descent recurrence + exact rational gate K.
"""
import json, sys, time
from pathlib import Path

_BASE = Path(__file__).resolve().parent
HERE = (_BASE / "files") if (_BASE / "files" / "declared.json").is_file() else _BASE
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import model_bn254_c6 as M
import verify_bn254_c6 as V

def main():
    t0 = time.time()
    raw = sys.stdin.read().strip() if not sys.stdin.isatty() else ""
    req = json.loads(raw) if raw else {}
    part = req.get("part", "default")
    decl = json.loads((HERE / "declared.json").read_text())
    track = decl["track"]

    if part == "filtsim":
        lu = int(req.get("log2U", 13))
        sd = int(req.get("seed", 20261120))
        fs = V.verify_live_filtsim_z6(log2U=min(max(lu, 10), 15), seed=sd)
        print(json.dumps(dict(track=track, part=part, checks_ok=fs["ok"], live_filtsim_z6=fs, elapsed_s=round(time.time() - t0, 3))))
        return

    alg = M.verify_tower_algebra_and_crt()
    fd = M.verify_exhaustive_fd(R2=50.0)
    proto = M.verify_archB_and_scalar_swar()
    walk7 = V.verify_7d_walk_prototype()
    lat12 = V.verify_12d_descent_lattices()
    filt_live = V.verify_live_filtsim_z6(log2U=13, seed=20261120)
    led = M.compute_ledger(track=track, use_db=True, c_tail_z6=420.0)

    row_checks = {}
    all_rows_ok = True
    for k, dv in decl["rows"].items():
        cv = led["rows"][k]
        ok = (cv <= float(dv) + 1e-12)
        row_checks[k] = dict(computed=round(cv, 6), declared=dv, ok=ok)
        all_rows_ok &= ok

    checks = dict(
        algebra_and_mod_r_crt_ok=alg["ok"],
        exhaustive_fd_bijection_ok=fd["ok"],
        archB_and_scalar_swar_ok=proto["ok"],
        walk7_fp_prototype_ok=walk7["ok"],
        descent_12d_lattices_ok=lat12["ok"],
        live_filtsim_z6_ok=filt_live["ok"],
        filtsim_3scale_ladder_ok=(led["la"]["ladder_monotonic"] and led["la"]["all_runs_below_4caps"]),
        live_62row_descent_ok=(led["descent"]["ok"] and led["descent"]["nodes"] == 62),
        rows_ok=all_rows_ok,
        exp_run_ok=(led["exp_run"] <= float(decl["exp_run"])),
        p_cap_h3_ok=(led["la"]["p_cap_h3"] <= float(decl["p_cap_h3"])),
        success_ok=(led["succ_exact"] >= float(decl["success"])),
        gate_K_ok=(led["K_gate"] == decl["gate_K"] and led["gate_passes"] and led["gate_prev_fails"]),
        claim_ok=(led["claim_exact"] <= float(decl["claim"]) and int(led["score_numerator"]) <= int(decl["score_numerator"])),
        ram_ok=(led["geometry"]["ram_log2"] < 100.0 and led["geometry"]["score_arena_log2"] < 64.0 and led["la"]["fits_1own"]),
        geo_consistency_ok=(abs(led["pool"]["frac_reg_analytic"] - led["pool"]["frac_reg"]) < 0.005 and led["geometry"]["L_max"] < 64.0 * M.OP["w_sb"]),
        peel_ok=(not led["la"]["peel"]["stalls"] and led["la"]["peel"]["exc_crit"] < 0.80),
    )
    checks_ok = all(checks.values())

    out = dict(
        track=track,
        part=part,
        checks_ok=checks_ok,
        all_ok=checks_ok,
        checks=checks,
        algebra=alg,
        exhaustive_fd=fd,
        archB_prototype=proto,
        walk7_prototype=walk7,
        descent_12d_lattices=lat12,
        live_filtsim_z6=filt_live,
        ledger=dict(
            rows=row_checks,
            exp_run=round(led["exp_run"], 9),
            exp_run_decl=decl["exp_run"],
            success=round(led["succ_exact"], 8),
            success_decl=decl["success"],
            gate_K=led["K_gate"],
            claim=round(led["claim_exact"], 9),
            claim_decl=decl["claim"],
            score_numerator=decl["score_numerator"],
            fail_terms=led["fail_terms"],
        ),
        geometry=led["geometry"],
        per_cell=led["per_cell"],
        la=led["la"],
        descent=led["descent"],
    )

    if part in ("default", "sensitivity", "full"):
        led_plain = M.compute_ledger(track=track, use_db=False, c_tail_z6=420.0)
        led_tight = M.compute_ledger(track=track, use_db=True, c_tail_z6=420.0, cap_scale=0.80)
        out["sensitivity"] = dict(
            plain_dickman_4cap=dict(exp_run=round(led_plain["exp_run"], 6), claim=round(led_plain["claim_exact"], 6)),
            de_bruijn_tighter_cap=dict(exp_run=round(led_tight["exp_run"], 6), claim=round(led_tight["claim_exact"], 6)),
        )

    out["elapsed_s"] = round(time.time() - t0, 3)
    print(json.dumps(out))

if __name__ == "__main__":
    main()

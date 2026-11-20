#!/usr/bin/env python3
"""
evidence.py — Standalone pure-stdlib verification entry point (python -I -B).
Executes live big-int tower algebra, live mod-r CRT factorization of h/f_rec/g_rec,
live exhaustive R^2=50 fundamental domain bijection (107/107, 0 missed, 0 dups),
live Arch B-Hybrid + 64-bit scalar SWAR prototype, live Euler-product retention,
live full de Bruijn 200k-pool yield, live 8-seed Direct-2^21 LA & Cantelli bound,
and live 62-row LB=59 bottom-up descent recurrence!
"""
import json, sys, time
from pathlib import Path

_BASE = Path(__file__).resolve().parent
HERE = (_BASE / "files") if (_BASE / "files" / "declared.json").is_file() else _BASE
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import model_bn254_c6 as M

def main():
    t0 = time.time()
    raw = sys.stdin.read().strip() if not sys.stdin.isatty() else ""
    req = json.loads(raw) if raw else {}
    part = req.get("part", "default")
    decl = json.loads((HERE / "declared.json").read_text())
    track = decl["track"]

    alg = M.verify_tower_algebra_and_crt()
    fd = M.verify_exhaustive_fd(R2=50.0)
    proto = M.verify_archB_and_scalar_swar()
    led = M.compute_ledger(track=track, use_db=True, direct_21_mean=True, c_tail_z6=420.0)

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
        live_62row_descent_ok=(led["descent"]["ok"] and led["descent"]["nodes"] == 62),
        rows_ok=all_rows_ok,
        exp_run_ok=(led["exp_run"] <= float(decl["exp_run"])),
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
        led_plain = M.compute_ledger(track=track, use_db=False, direct_21_mean=True, c_tail_z6=420.0)
        led_ms = M.compute_ledger(track=track, use_db=True, direct_21_mean=False, c_tail_z6=420.0)
        led_plain_ms = M.compute_ledger(track=track, use_db=False, direct_21_mean=False, c_tail_z6=207.0)
        out["sensitivity"] = dict(
            plain_rho_direct21=dict(exp_run=round(led_plain["exp_run"], 6), claim=round(led_plain["claim_exact"], 6)),
            de_bruijn_matched_scale_la=dict(exp_run=round(led_ms["exp_run"], 6), claim=round(led_ms["claim_exact"], 6)),
            plain_rho_matched_scale_tail207=dict(exp_run=round(led_plain_ms["exp_run"], 6), claim=round(led_plain_ms["claim_exact"], 6)),
        )

    out["elapsed_s"] = round(time.time() - t0, 3)
    print(json.dumps(out))

if __name__ == "__main__":
    main()

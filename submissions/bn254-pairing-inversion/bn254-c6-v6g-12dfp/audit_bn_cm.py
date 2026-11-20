"""
audit_bn_cm.py — 100% self-contained pure-stdlib verification entry point (python3 -I -B) for BN254 v6g.
Executes live (<= 5s total):
1. Big-int BN254 + C6-m184_w216-rec algebra, exact rational sigma(w) == w, & mod-r CRT across all 6 roots;
2. Exhaustive R^2=50 fundamental-domain bijection (107/107 orbits, 0 missed, 0 dups);
3. 8-lane SWAR sb_u7 Band-2 prototype (1523/1523, 0 mismatches) + 10k uint16 SWAR scan+zero + 417-chunk P_14
 score-only Stage-B strip check (0 errors, ZERO Band-1/2/3 prime recording!);
4. Live sigma-orbit-folded 12D Optimal Bounding Ellipsoid E* Fincke-Pohst Band-3 prototype (4 real C6 orbits, exact active R, 408/408, 0 missed, 0 dups, chunked Q1.60 gate);
5. Live 12D C6 descent lattice audit (24 special-q 12f+12g 2^66..2^213 & 8 genuine G_T boot witnesses on exact FULL B_12(4), 6/6 sextiles + 5/5 quintiles);
6. Live pure-Python production hexagonal Z[omega] filtering simulator (U'=2^13 with kernel check & U'=2^14)
 + w_mid -> 0 monotone check + 5-scale ALL4 ladder verification;
7. Live UNPLANTED 2-sided 5-SM NFS Block-Wiedemann + F_1187^* DLP recovery + Z[omega] Schur check (hx_z6.py);
8. Full de Bruijn rho_dB 16-row ledger (89,151 all-active-cell cof_b, C_HIT_LO=13.25, C_FP12_CHUNK_INIT=120.0),
 fixed-instance H3 4-cap 1-owner LA, live 62-row LB=59 ZERO-RETRY DETERMINISTIC-DAG descent (k=1, boot_rels=1),
 strict FIX-U additive union bound, & exact rational run gate K.
"""
import json, sys, time
from pathlib import Path
_BASE = Path(__file__).resolve().parent
HERE = (_BASE / "files") if (_BASE / "files" / "declared.json").is_file() else _BASE
if str(HERE) not in sys.path:
 sys.path.insert(0, str(HERE))
import model_bn254_v6 as M
import verify_bn254_v6 as V
import hx_z6 as HX
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
  fs = V.verify_live_filtsim_z6(log2U=min(max(lu, 10), 15), seed=sd, W_target=160, check_kernel=(lu <= 13))
  print(json.dumps(dict(track=track, part=part, checks_ok=fs["ok"], live_filtsim_z6=fs, elapsed_s=round(time.time() - t0, 3))))
  return
 if part == "hx":
  lu = int(req.get("log2U", 12))
  sd = int(req.get("seed", 1))
  m = int(req.get("m", 10))
  n = int(req.get("n", 8))
  hx = HX.run_hx_1d_kernel_z6(log2U=min(max(lu, 11), 13), seed=sd, m=m, n=n)
  print(json.dumps(dict(track=track, part=part, checks_ok=hx["ok"], hx_z6=hx, elapsed_s=round(time.time() - t0, 3))))
  return
 alg = M.verify_tower_algebra_and_crt()
 fd = M.verify_exhaustive_fd(R2=50.0)
 proto = V.verify_swar_sb_u7_and_stageB()
 walk12 = V.verify_12d_sigma_folded_walk()
 lat12 = V.verify_12d_descent_lattices()
 filt_live = V.verify_live_filtsim_ladder(seed=20261120)
 hx = HX.run_hx_1d_kernel_z6(log2U=12, seed=1, m=10, n=8)
 led = M.compute_ledger(track=track, use_db=True, c_tail_z6=420.0)
 pq_bins = M.load_data()["pool_quad"]["bins"]
 pool_all_active_ok = (sum(x[2] for x in pq_bins) == 61463 and sum(x[3] for x in pq_bins) == 89151)
 row_checks = {}
 all_rows_ok = True
 for k, dv in decl["rows"].items():
  cv = led["rows"][k]
  ok = (cv <= float(dv) + 1e-12)
  row_checks[k] = dict(computed=round(cv, 6), declared=dv, ok=ok)
  all_rows_ok &= ok
 ladder5_ok = (
  filt_live["ok"] and led["la"]["ladder_all4_monotonic"] and led["la"]["all_runs_below_4caps"] and
  filt_live["u14"]["N_sw_U"] > led["la"]["filtsim_cur20"]["N_sw_U"] and
  filt_live["u14"]["N_l_U"] > led["la"]["filtsim_cur20"]["N_l_U"] and
  filt_live["u14"]["k5_exp"] > led["la"]["filtsim_cur20"]["k5_exp"] and
  filt_live["u14"]["tab_r"] > led["la"]["filtsim_cur20"]["tab_r"]
 )
 checks = dict(
  algebra_sigma_w_and_crt_ok=alg["ok"],
  exhaustive_fd_bijection_ok=fd["ok"],
  swar_sb_u7_and_stageB_ok=proto["ok"],
  walk12_sigma_folded_ok=walk12["ok"],
  descent_12d_lattices_24sq_8boot_ok=lat12["ok"],
  pool_all_active_89151_ok=pool_all_active_ok,
  live_filtsim_z6_and_wmid_ok=filt_live["ok"],
  filtsim_5scale_ladder_ok=ladder5_ok,
  unplanted_5sm_bw_dlp_and_schur_ok=hx["ok"],
  live_62row_zero_retry_dag_ok=(
   led["descent"]["ok"] and led["descent"]["nodes"] == 62
   and led["descent"]["mode"] == "zero_retry_deterministic_dag_k1_boot1"
   and led["descent"]["repeat_prob"] == 0.0
  ),
  rows_ok=all_rows_ok,
  exp_run_ok=(led["exp_run"] <= float(decl["exp_run"])),
  success_ok=(led["succ_exact"] >= float(decl["success"])),
  gate_K_ok=(led["K_gate"] == decl["gate_K"] and led["gate_passes"] and led["gate_prev_fails"]),
  claim_ok=(led["claim_exact"] <= float(decl["claim"]) and int(led["score_numerator"]) <= int(decl["score_numerator"])),
  ram_ok=(led["geometry"]["ram_log2"] < 100.0 and led["geometry"]["score_arena_log2"] < 64.0
    and led["geometry"]["b2_1p_ram_log2"] < 64.0 and led["la"]["fits_1own"]
    and led["descent"]["boot_root_bytes_log2"] < 64.0),
  geo_consistency_ok=(abs(led["pool"]["frac_reg_analytic"] - led["pool"]["frac_reg"]) < 0.005
       and led["geometry"]["max_sb"] <= 8),
  peel_ok=(not led["la"]["peel"]["stalls"] and led["la"]["peel"]["exc_crit"] < 0.75),
 )
 checks_ok = all(checks.values())
 out = dict(
  track=track, part=part, checks_ok=checks_ok, all_ok=checks_ok, checks=checks,
  algebra=alg, exhaustive_fd=fd, swar_sb_u7_stageB=proto,
  walk12_sigma_folded=walk12, descent_12d_lattices=lat12,
  live_filtsim_z6=filt_live, hx_unplanted_5sm_bw=hx,
  ledger=dict(
   rows=row_checks,
   exp_run=round(led["exp_run"], 9), exp_run_decl=decl["exp_run"],
   success=round(led["succ_exact"], 8), success_decl=decl["success"],
   gate_K=led["K_gate"],
   claim=round(led["claim_exact"], 9), claim_decl=decl["claim"],
   score_numerator=decl["score_numerator"], fail_terms=led["fail_terms"],
  ),
  geometry=led["geometry"], per_cell=led["per_cell"], la=led["la"], descent=led["descent"],
 )
 if part in ("sensitivity", "full"):
  led_plain = M.compute_ledger(track=track, use_db=False)
  led_bkt25 = M.compute_ledger(track=track, use_db=True, strict_per_owner=True, bkt_b3=25.0)
  led_flat = M.compute_ledger(track=track, use_db=True, strict_per_owner=False, bkt_b3=20.0)
  led_tight = M.compute_ledger(track=track, use_db=True, cap_scale=0.80)
  out["sensitivity"] = dict(
   plain_dickman_rho=dict(exp_run=round(led_plain["exp_run"], 6), claim=round(led_plain["claim_exact"], 6)),
   conservative_bkt25=dict(exp_run=round(led_bkt25["exp_run"], 6), claim=round(led_bkt25["claim_exact"], 6)),
   flat_shared_arena_ref=dict(exp_run=round(led_flat["exp_run"], 6), claim=round(led_flat["claim_exact"], 6)),
   de_bruijn_tighter_cap=dict(exp_run=round(led_tight["exp_run"], 6), claim=round(led_tight["claim_exact"], 6)),
  )
 out["elapsed_s"] = round(time.time() - t0, 3)
 print(json.dumps(out))
 if not checks_ok:
  raise SystemExit(1)
if __name__ == "__main__":
 main()

"""audit_bn_cm.py — Self-contained pure-stdlib live verifier (python3 -I -B, 18/18 checks) for BN254 v6p-noband3."""
import json,sys,time
from pathlib import Path
_BASE=Path(__file__).resolve().parent
HERE=(_BASE/"files")if(_BASE/"files"/"declared.json").is_file()else _BASE
if str(HERE)not in sys.path:
 sys.path.insert(0,str(HERE))
import model_bn254_v6 as M
import verify_bn254_v6 as V
import hx_z6 as HX
def main():
 t0=time.time()
 raw=sys.stdin.read().strip()if not sys.stdin.isatty()else""
 req=json.loads(raw)if raw else{}
 part=req.get("part","default")
 decl=json.loads((HERE/"declared.json").read_text())
 track=decl["track"]
 if part=="filtsim":
  lu=int(req.get("log2U",13))
  sd=int(req.get("seed",20261120))
  fs=V.verify_live_filtsim_z6(log2U=min(max(lu,10),15),seed=sd,W_target=160,check_kernel=(lu<=13))
  print(json.dumps(dict(track=track,part=part,checks_ok=fs["ok"],live_filtsim_z6=fs,elapsed_s=round(time.time()-t0,3))))
  return
 if part=="hx":
  lu=int(req.get("log2U",12))
  sd=int(req.get("seed",1))
  m=int(req.get("m",10))
  n=int(req.get("n",8))
  hx=HX.run_hx_1d_kernel_z6(log2U=min(max(lu,11),13),seed=sd,m=m,n=n)
  print(json.dumps(dict(track=track,part=part,checks_ok=hx["ok"],hx_z6=hx,elapsed_s=round(time.time()-t0,3))))
  return
 alg=M.verify_tower_algebra_and_crt()
 fd=M.verify_exhaustive_fd(R2=50.0)
 proto=V.verify_swar_sb_u7_and_stageB()
 nonprim_hnf=M.verify_nonprim_hnf_witnesses()
 lat12=V.verify_12d_descent_lattices()
 pgf_camp=V.verify_c6_pgf_and_hea_campaign()
 filt_live=V.verify_live_filtsim_ladder(seed=20261120)
 hx=HX.run_hx_1d_kernel_z6(log2U=12,seed=1,m=10,n=8)
 led=M.compute_ledger(track=track,use_db=True,c_tail_z6=420.0)
 pq=M.load_data()["pool_quad"]
 pq_bins=pq["bins"];np_bins=pq["nonprim_bins"]
 pool_all_active_ok=(sum(x[2]for x in pq_bins)==61051 and sum(x[1]for x in np_bins)==28100 and sum(x[3]for x in pq_bins)==89151 and nonprim_hnf["ok"])
 row_checks={}
 all_rows_ok=True
 for k,dv in decl["rows"].items():
  cv=led["rows"][k]
  ok=(cv<=float(dv)+1e-12)
  row_checks[k]=dict(computed=round(cv,6),declared=dv,ok=ok)
  all_rows_ok&=ok
 ladder5_ok=(
 filt_live["ok"]and led["la"]["ladder_all4_monotonic"]and led["la"]["all_runs_below_4caps"]and
 filt_live["u14"]["N_sw_U"]>led["la"]["filtsim_cur20"]["N_sw_U"]and
 filt_live["u14"]["N_l_U"]>led["la"]["filtsim_cur20"]["N_l_U"]and
 filt_live["u14"]["k5_exp"]>led["la"]["filtsim_cur20"]["k5_exp"]and
 filt_live["u14"]["tab_r"]>led["la"]["filtsim_cur20"]["tab_r"]
 )
 checks=dict(
 algebra_sigma_w_and_crt_ok=alg["ok"],
 exhaustive_fd_bijection_ok=fd["ok"],
 swar_sb_u7_and_stageB_ok=proto["ok"],
 descent_12d_lattices_24sq_8boot_ok=lat12["ok"],
 c6_pgf_usable_depth2_and_hea_campaign_ok=(pgf_camp["ok"]and pgf_camp["live_4_depth2_2hop_chains_exact"]),
 pool_all_active_89151_ok=pool_all_active_ok,
 live_filtsim_z6_and_wmid_ok=filt_live["ok"],
 filtsim_5scale_ladder_ok=ladder5_ok,
 unplanted_5sm_bw_dlp_and_schur_ok=hx["ok"],
 live_62row_k1_root_b10_metered_dag_ok=(
 led["descent"]["ok"]and led["descent"]["nodes"]==62
 and led["descent"]["mode"]=="zero_retry_dag_k1_root_iid_b10_metered"
 and led["descent"]["nb_r"]==1.0 and led["descent"]["nb_alpha"]==1.0
 and led["descent"]["boot_retries_max"]==10
 and led["descent"]["phi_run_one"]<=0.35
 and led["descent"]["phi_log_two"]<1e-4
 and led["descent"]["desc_cap_log2"]==112.0
 and led["descent"]["markov_cap_fail"]<1e-8
 and led["descent"]["repeat_prob"]==0.0
 ),
 rows_ok=all_rows_ok,
 exp_run_ok=(led["exp_run"]<=float(decl["exp_run"])),
 success_ok=(led["succ_exact"]>=float(decl["success"])),
 gate_K_ok=(led["K_gate"]==decl["gate_K"]and led["gate_passes"]and led["gate_prev_fails"]),
 claim_ok=(led["claim_exact"]<=float(decl["claim"])and int(led["score_numerator"])<=int(decl["score_numerator"])),
 ram_ok=(led["geometry"]["ram_log2"]<100.0 and led["geometry"]["score_arena_log2"]<64.0
 and led["geometry"]["b2_1p_ram_log2"]<64.0 and led["la"]["fits_1own"]
 and led["descent"]["boot_root_bytes_log2"]<64.0),
 geo_consistency_ok=(abs(led["pool"]["frac_reg_analytic"]-led["pool"]["frac_reg"])<0.005
 and led["geometry"]["max_sb"]<=8),
 peel_ok=(not led["la"]["peel"]["stalls"]and led["la"]["peel"]["exc_crit"]<0.75),
 )
 checks_ok=all(checks.values())
 out=dict(
 track=track,part=part,checks_ok=checks_ok,all_ok=checks_ok,checks=checks,
 algebra=alg,exhaustive_fd=fd,swar_sb_u7_stageB=proto,
 nonprim_hnf_witnesses=nonprim_hnf,descent_12d_lattices=lat12,
 c6_pgf_campaign=pgf_camp,
 live_filtsim_z6=filt_live,hx_unplanted_5sm_bw=hx,
 ledger=dict(
 rows=row_checks,
 exp_run=round(led["exp_run"],9),exp_run_decl=decl["exp_run"],
 success=round(led["succ_exact"],8),success_decl=decl["success"],
 gate_K=led["K_gate"],
 claim=round(led["claim_exact"],9),claim_decl=decl["claim"],
 score_numerator=decl["score_numerator"],fail_terms=led["fail_terms"],
 ),
 geometry=led["geometry"],per_cell=led["per_cell"],la=led["la"],descent=led["descent"],
 )
 if part in("sensitivity","full"):
  led_plain=M.compute_ledger(track=track,use_db=False)
  led_swar22=M.compute_ledger(track=track,use_db=True,strict_per_owner=True,swar_probe=22.0)
  led_flat=M.compute_ledger(track=track,use_db=True,strict_per_owner=False)
  led_tight=M.compute_ledger(track=track,use_db=True,cap_scale=0.80)
  out["sensitivity"]=dict(
  plain_dickman_rho=dict(exp_run=round(led_plain["exp_run"],6),claim=round(led_plain["claim_exact"],6)),
  conservative_swar22=dict(exp_run=round(led_swar22["exp_run"],6),claim=round(led_swar22["claim_exact"],6)),
  flat_shared_arena_ref=dict(exp_run=round(led_flat["exp_run"],6),claim=round(led_flat["claim_exact"],6)),
  de_bruijn_tighter_cap=dict(exp_run=round(led_tight["exp_run"],6),claim=round(led_tight["claim_exact"],6)),
  )
 out["elapsed_s"]=round(time.time()-t0,3)
 print(json.dumps(out))
 if not checks_ok:
  raise SystemExit(1)
if __name__=="__main__":
 main()

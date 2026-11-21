"""audit_bn_cm.py — BN254 C6 v19 live verifier (28/28 checks)."""
import json,sys,time
from pathlib import Path
_BASE=Path(__file__).resolve().parent
HERE=(_BASE/"files")if(_BASE/"files"/"declared.json").is_file()else _BASE
if str(HERE)not in sys.path:
 sys.path.insert(0,str(HERE))
import model_bn254_v6 as M
import verify_bn254_v6 as V
import descent_bn254_v6 as DESC
import hx_z6 as HX
import c6_adaptive as ADAP
import v19_lat_b2 as V19L
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
 pool_sub=V.verify_live_c6_pool_hnf_subsample(n_tot=1000,seed=20261121)
 gen_filt=V.verify_live_genuine_c6_relation_filter(seed=20261120)
 bareiss_lll=V.verify_live_bareiss_12d_lll(seed=20261120)
 mont5_ecm=DESC.verify_live_mont5_and_suyama12_ecm(seed=20261120)
 la_kern=DESC.verify_live_z6_la_kernels(seed=20261120)
 lat12=V.verify_12d_descent_lattices()
 pgf_camp=V.verify_c6_pgf_and_hea_campaign()
 deg8_camp=DESC.verify_deg8_dispersion_and_prod_yield(seed=20261120)
 adap_h5=ADAP.verify_c6_adaptive_h5_live()
 h3_scale=ADAP.verify_h3_scaling_and_nu_h()
 crt_mp=ADAP.verify_split_crt_norms_and_mp_rt(seed=20261120)
 filt_live=V.verify_live_filtsim_ladder(seed=20261120)
 hx=HX.run_hx_1d_kernel_z6(log2U=11,seed=1,m=10,n=8)
 lat_b2=V19L.verify_v19_lattice_b2_swar_and_supertile(seed=20261121)
 led=M.compute_ledger(track=track,use_db=True,c_tail_z6=740.0)
 pq=M.load_data()["pool_quad"]
 pq_bins=pq["bins"];np_bins=pq["nonprim_bins"]
 pool_all_active_ok=(sum(x[2]for x in pq_bins)==pq["n_prim_in"]and sum(x[1]for x in np_bins)==pq["n_nonprim_in"]and pq["n_in"]==200000 and pq["n_prim_in"]+pq["n_nonprim_in"]==200000 and nonprim_hnf["ok"])
 row_checks={}
 all_rows_ok=True
 for k,dv in decl["rows"].items():
  cv=led["rows"][k]
  ok=(cv<=float(dv)+1e-12)
  row_checks[k]=dict(computed=round(cv,6),declared=dv,ok=ok)
  all_rows_ok&=ok
 ladder7_ok=(filt_live["ok"]and led["la"]["ladder_all4_monotonic"]and led["la"]["all_runs_24_to_26_below_4caps"]and
 filt_live["u14"]["N_sw_U"]>led["la"]["filtsim_cur20"]["N_sw_U"]and
 filt_live["u14"]["N_l_U"]>led["la"]["filtsim_cur20"]["N_l_U"]and
 filt_live["u14"]["k5_exp"]>led["la"]["filtsim_cur20"]["k5_exp"]and
 filt_live["u14"]["tab_r"]<M.TAB_R_CAP
 )
 checks=dict(algebra_sigma_w_and_crt_ok=(alg["ok"]and alg["ntt3_shoup_butterfly_ok"]and alg["ntt3_garner_crt_ok"]and alg["ntt3_r4_r3_shoup_conv_ok"]),exhaustive_fd_pure_int_bijection_ok=(fd["ok"]and fd["zero_float_geom"]and fd["sigma_inv_exact"]and fd["max_order_dedekind_ok"]and fd["exact_frac_enum_ok"]and fd["certified_interval_ok"]and fd["global_hexagon_proof_ok"]and fd["unit_index1_proved"]),
 swar_sb_u7_pure_int_and_stageB_ok=(proto["ok"]and proto["det_occ_bound_ok"]and proto["p14_zero_gcd_gm_ok"]and proto["zero_float_geom"]and proto["flat_array_zero_dict_build"]and proto["flat_sbu7_exact_ok"]and proto["zero_sort_exact_ok"]and proto["exact_slot_cutoff_M_t"]and not proto["qa_stored_in_record"]),
 v19_lattice_b2_swar_and_supertile_ok=(lat_b2["ok"]and lat_b2["zero_float_fp_enum"]and lat_b2["strict_u64_i64_umulh64_enforced"]and lat_b2["fallback_128bit_triggered"]and lat_b2["fallback_128bit_mismatches"]==0 and lat_b2["isqrt64_5newton_errs"]==0 and lat_b2["lattice_b2_mismatches"]==0 and lat_b2["prod_51_to_61b_guards_passed"]==6 and lat_b2["prod_51_to_61b_actual_leaves_total"]>=500000 and lat_b2["prod_51_to_61b_actual_max_ir"]<=25.0 and lat_b2["prod_max_inv_b2_bits"]<=62 and lat_b2["prod_max_b2_q64_bits"]<=64 and lat_b2["prod_max_r2_q12_bits"]<64 and lat_b2["swar_step_errs"]==0 and lat_b2["swar_20op_canon_errs"]==0 and lat_b2["swar_tie_canon_errs"]==0 and lat_b2["supertile_replay_errs"]==0 and lat_b2["desc_ai_prefix_chain_errs"]==0 and lat_b2["swar_bounds_ok"]and lat_b2["supertile_structured_ht_avg_probes"]<=1.12),
 live_pool_hnf_subsample_ok=pool_sub["ok"],live_genuine_c6_rel_filter_ok=(gen_filt["ok"]and gen_filt["post_merge_cols_lost"]==0 and gen_filt["merge_backsub_errors"]==0),live_bareiss_12d_lll_limb_ok=(bareiss_lll["ok"]and bareiss_lll["inductive_sr_inv_ok"]and bareiss_lll["input_uniform_wc_proved"]and bareiss_lll["pair_setup_p59_proved"]and bareiss_lll["lat_b2_p61_proved"]),live_mont5_cios_and_suyama12_ecm_ok=(mont5_ecm["ok"]and mont5_ecm["ecm_ladder6_sep_ok"]and mont5_ecm["ecm_ladder6_min_sep_dag"]>=0.60 and mont5_ecm["ecm_ladder6_min_sep_boot"]>=0.60),
 live_z6_la_5limb_barrett_ok=(la_kern["ok"]and la_kern["unit4_v4_expand_ok"]),descent_12d_lattices_24sq_24boot_ok=(lat12["ok"]and lat12["prenorm_gate_uniform_ok"]and lat12["all_24_sq_reduced_live_by_bareiss_910"]and lat12["all_24_boot_bareiss_910_certified"]and lat12["boot_24_b_best_prim_ok"]),
 c6_pgf_usable_depth2_and_hea_campaign_ok=(pgf_camp["ok"]and pgf_camp["exact_ideal_primitivity_checked"]and pgf_camp["live_4_depth2_2hop_chains_exact"]and pgf_camp["unified_bareiss_910_used"]),
 deg8_dispersion_and_prod_yield_ok=(deg8_camp["ok"]and deg8_camp["tower_u5_p27631_r27481_deg8_frec_verified"]and deg8_camp["all_100_above_fb_1024"]and deg8_camp["live_6_deg8_ideals_and_2hop_d2_exact"]),
 c6_adaptive_h5_live_ok=(adap_h5["ok"]and adap_h5["strict_k1_zero_retry_protocol"]and adap_h5["live_k1_zero_retry_ok"]and adap_h5["live_witnesses_replayed"]==3 and adap_h5["after_fail_total"]==1856 and adap_h5["pooled_oe_fail"]<0.85),
 h3_scaling_and_nu_h_ok=(h3_scale["ok"]and h3_scale["fixed_instance_deterministic_h3"]and h3_scale["zero_loss_filtsim_ok"]and h3_scale["all4_monotone_20_26"]and h3_scale["hist_v4_exact_ok"]and h3_scale["fg_dens_ok"]and h3_scale["u_ram_all_lt_2p30"]),
 split_crt_norms_and_mp_rt_ok=(crt_mp["ok"]and crt_mp["chained_1pass_b2_errs"]==0),
 pool_all_active_200k_ok=pool_all_active_ok,live_filtsim_z6_and_wmid_ok=(filt_live["ok"]and filt_live["u12"]["post_merge_cols_lost"]==0 and filt_live["u12"]["merge_backsub_errors"]==0 and filt_live["u13"]["post_merge_cols_lost"]==0),filtsim_7scale_ladder_ok=ladder7_ok,unplanted_5sm_bw_dlp_and_schur_ok=hx["ok"],
 live_62row_pure_k1_boole_wald_dag_ok=(led["descent"]["ok"]and led["descent"]["nodes"]==62
 and led["descent"]["mode"]=="pure_k1_boole_wald_iid_boot"
 and led["descent"]["p0_floor"]==0.001
 and led["descent"]["phi_dag_one"]<=0.60
 and led["descent"]["e_runs_per_log"]<=2.0
 and led["descent"]["mu_min"]>=500.0
 and led["descent"]["p0_max"]<=0.005
 and led["descent"]["p_node_abort"]==0.0
 and led["descent"]["p_repeat_loss"]==0.0
 and led["descent"]["p_descent_fail"]==0.0
 and led["descent"]["desc_cap_log2"]==112.0
 ),rows_ok=all_rows_ok,exp_run_ok=(led["exp_run"]<=float(decl["exp_run"])),success_ok=(led["succ_exact"]>=float(decl["success"])),gate_K_ok=(led["K_gate"]==decl["gate_K"]and led["gate_passes"]and led["gate_prev_fails"]),claim_ok=(led["claim_exact"]<=float(decl["claim"])and int(led["score_numerator"])<=int(decl["score_numerator"])),
 ram_ok=(led["geometry"]["ram_log2"]<100.0 and led["geometry"]["score_arena_log2"]<64.0
 and led["geometry"]["b0_1own_ram_log2"]<64.0
 and led["geometry"]["b2_1p_ram_log2"]<64.0 and led["geometry"]["b2_build_1p_ram_log2"]<64.0
 and led["la"]["fits_1own"]and led["descent"]["boot_root_bytes_log2"]<64.0),geo_consistency_ok=(abs(led["pool"]["frac_reg_analytic"]-led["pool"]["frac_reg"])<0.005
 and led["geometry"]["max_sb"]<=8 and led["geometry"]["sup_nf_bits"]<448.0 and led["geometry"]["sup_ng_bits"]<=512.0 and led["geometry"]["log2_Na_shell"]<46.0),peel_ok=(not led["la"]["peel"]["stalls"]and led["la"]["peel"]["exc_crit"]<0.75),
 )
 checks_ok=all(checks.values())
 out=dict(track=track,part=part,checks_ok=checks_ok,all_ok=checks_ok,checks=checks,algebra=alg,exhaustive_fd=fd,swar_sb_u7_stageB=proto,v19_lattice_b2=lat_b2,live_pool_hnf_subsample=pool_sub,live_genuine_c6_rel_filter=gen_filt,
 live_bareiss_12d_lll=bareiss_lll,live_mont5_cios_and_suyama12_ecm=mont5_ecm,live_z6_la_kernels=la_kern,nonprim_hnf_witnesses=nonprim_hnf,descent_12d_lattices=lat12,
 c6_pgf_campaign=pgf_camp,deg8_dispersion_and_prod_yield=deg8_camp,c6_adaptive_h5_live=adap_h5,h3_scaling_and_nu_h=h3_scale,split_crt_norms_and_mp_rt=crt_mp,live_filtsim_z6=filt_live,hx_unplanted_5sm_bw=hx,ledger=dict(rows=row_checks,exp_run=round(led["exp_run"],9),exp_run_decl=decl["exp_run"],success=round(led["succ_exact"],8),success_decl=decl["success"],gate_K=led["K_gate"],claim=round(led["claim_exact"],9),claim_decl=decl["claim"],
 score_numerator=decl["score_numerator"],fail_terms=led["fail_terms"],),geometry=led["geometry"],per_cell=led["per_cell"],la=led["la"],descent=led["descent"],)
 if part in("sensitivity","full"):
  led_plain=M.compute_ledger(track=track,use_db=False)
  led_lat55=M.compute_ledger(track=track,use_db=True,strict_per_owner=True,swar_probe=10.0)
  led_flat=M.compute_ledger(track=track,use_db=True,strict_per_owner=False)
  led_tight=M.compute_ledger(track=track,use_db=True,cap_scale=0.80)
  out["sensitivity"]=dict(plain_dickman_rho=dict(exp_run=round(led_plain["exp_run"],6),claim=round(led_plain["claim_exact"],6)),conservative_lat_b2_plus10=dict(exp_run=round(led_lat55["exp_run"],6),claim=round(led_lat55["claim_exact"],6)),flat_shared_arena_ref=dict(exp_run=round(led_flat["exp_run"],6),claim=round(led_flat["claim_exact"],6)),
  de_bruijn_tighter_cap=dict(exp_run=round(led_tight["exp_run"],6),claim=round(led_tight["claim_exact"],6)),)
 out["elapsed_s"]=round(time.time()-t0,3)
 print(json.dumps(out))
 if not checks_ok:
  raise SystemExit(1)
if __name__=="__main__":
 main()

"""audit_bn_cm.py — BN254 C6 v29 live verifier (pure stdlib). Default run: the v26 checks (v20_lat_b2.py, v25_verify.py, unchanged v23b checks, 16-row ledger) plus the boot checks of ADAPTIVE.md section 1.2 (one randomly chosen stored production boot lattice replayed by boot_dp.py, the kept-population small-prime statistics and the miniature ratios recomputed from all stored records, and the boot ledger fields). Input {"part":"diagnostics"}: the 3 H5 scoped-diagnostic checks of ADAPTIVE.md section 2 (12D descent lattices incl. 24 stored boot lattices with all 144 candidates gated, adaptive h5, deg-8 dispersion), kept out of the default run by the board's per-run time limit. Inputs {"part":"boot_dp","j0":j,"n":<=8} and {"part":"boot_mini","u":-3|5,"j0":j,"n":<=40} replay any stored range; {"part":"dag_survivors"} prints the claim with twice the DAG survivor envelope (H5-W sensitivity). Other parts: sensitivity, full, filtsim, hx."""
import json,math,sys,time
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
import v20_lat_b2 as V20L
import v25_verify as V25V
import boot_dp as BD
import boot_mini as BM
import v32_verify as V32
import os
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
 if part=="dag_survivors":
  DESC.DAG_SURV_X=2.0;led2=M.compute_ledger(track=track,use_db=True,c_tail_z6=328.0)
  print(json.dumps(dict(track=track,part=part,checks_ok=True,dag_survivor_envelope_x2=dict(claim=round(led2["claim_exact"],6),exp_run=round(led2["exp_run"],6),two_dag=round(led2["descent"]["two_dag"],4),two_logs=round(led2["descent"]["two_log"],4)),elapsed_s=round(time.time()-t0,3))))
  return
 if part in("boot_dp","boot_mini"):
  u=int(req.get("u",-3));u=u if str(u)in BM.load()else -3
  N=256 if part=="boot_dp"else BM.load()[str(u)]["n"];j0=min(max(int(req.get("j0",0)),0),N-1)
  n=min(max(int(req.get("n",4)),1),8 if part=="boot_dp"else 40,N-j0)
  r=BD.replay(list(range(j0,j0+n)),full=True)if part=="boot_dp"else BM.replay(u,list(range(j0,j0+n)))
  print(json.dumps(dict(track=track,part=part,checks_ok=r["ok"],replay=r,elapsed_s=round(time.time()-t0,3))))
  return
 if part=="boot_evidence":
  rm=BD.rough_med_audit(range(4),1<<20);jm3=BD.joint_ea_mini_audit(-3,80)
  ok=bool(rm["ok"]and jm3["ok"])
  print(json.dumps(dict(track=track,part=part,checks_ok=ok,prod_rough_16_20=rm,mini_joint_ea_m3=jm3,elapsed_s=round(time.time()-t0,3))))
  if not ok:raise SystemExit(1)
  return
 if part=="hx":
  lu=int(req.get("log2U",12))
  sd=int(req.get("seed",1))
  m=int(req.get("m",10))
  n=int(req.get("n",8))
  hx=HX.run_hx_1d_kernel_z6(log2U=min(max(lu,11),13),seed=sd,m=m,n=n)
  print(json.dumps(dict(track=track,part=part,checks_ok=hx["ok"],hx_z6=hx,elapsed_s=round(time.time()-t0,3))))
  return
 diag=(part=="diagnostics")
 alg=M.verify_tower_algebra_and_crt() if not diag else None
 fd=M.verify_exhaustive_fd(R2=50.0) if not diag else None
 proto=V.verify_swar_sb_u7_and_stageB() if not diag else None
 nonprim_hnf=V25V.verify_nonprim_hnf_witnesses() if not diag else None
 pool_sub=V.verify_live_c6_pool_hnf_subsample(n_tot=1000,seed=20261121) if not diag else None
 gen_filt=V.verify_live_genuine_c6_relation_filter(seed=20261120) if not diag else None
 bareiss_lll=V.verify_live_bareiss_12d_lll(seed=20261120) if not diag else None
 mont5_ecm=DESC.verify_live_mont5_and_suyama12_ecm(seed=20261120) if not diag else None
 la_kern=DESC.verify_live_z6_la_kernels(seed=20261120) if not diag else None
 lat12=V.verify_12d_descent_lattices()if diag else None
 if diag:
  pgf_camp=V.verify_c6_pgf_and_hea_campaign()
  deg8_camp=DESC.verify_deg8_dispersion_and_prod_yield(seed=20261120)
  adap_h5=ADAP.verify_c6_adaptive_h5_live()
  checks=dict(descent_12d_lattices_24sq_24boot_ok=(lat12["ok"]and lat12["prenorm_gate_uniform_ok"]and lat12["all_24_sq_reduced_live_by_bareiss_910"]and lat12["all_24_boot_bareiss_910_certified"]and lat12["boot_candidates_gated"]==3456),deg8_dispersion_and_prod_yield_ok=(deg8_camp["ok"]and deg8_camp["tower_u5_p27631_r27481_deg8_frec_verified"]and deg8_camp["all_100_above_fb_1024"]and deg8_camp["live_6_deg8_ideals_and_2hop_d2_exact"]),c6_adaptive_h5_live_ok=(adap_h5["ok"]and adap_h5["strict_k1_zero_retry_protocol"]and adap_h5["live_k1_zero_retry_ok"]and adap_h5["live_witnesses_replayed"]==3 and adap_h5["after_fail_total"]==1856 and adap_h5["pooled_oe_fail"]<0.85))
  checks_ok=all(checks.values())
  print(json.dumps(dict(track=track,part=part,checks_ok=checks_ok,all_ok=checks_ok,checks=checks,descent_12d_lattices=lat12,deg8_dispersion_and_prod_yield=deg8_camp,c6_adaptive_h5_live=adap_h5,elapsed_s=round(time.time()-t0,3))))
  if not checks_ok:
   raise SystemExit(1)
  return
 pgf_camp=V.verify_c6_pgf_and_hea_campaign()
 h3_scale=ADAP.verify_h3_scaling_and_nu_h()
 crt_mp=ADAP.verify_split_crt_norms_and_mp_rt(seed=20261120)
 filt_live=V.verify_live_filtsim_ladder(seed=20261120)
 hx=HX.run_hx_1d_kernel_z6(log2U=11,seed=1,m=10,n=8)
 V20L.HT_SEED_BASE=int(req.get("ht_run",0))or 1+int.from_bytes(os.urandom(2),"big")%16
 lat_b2=V20L.verify_v20_lat_b2_and_srt(seed=20261121)
 ht_pool=V25V.verify_ht_pool(lat_b2)
 v25_upg=V25V.verify_v25_upgrades(seed=20261122)
 v32_upg=V32.verify_v32_all(M,str(HERE/"kernels_v32.s"),seed=20261123)
 jd=int.from_bytes(os.urandom(2),"big")%256;bd_rep=BD.replay([jd]);bd_a=BD.alpha16()
 bm_st=BM.stats();bm_rep=dict(status="not_run",how_to_run={"part":"boot_mini","u":5,"j0":0,"n":40},note="regenerates stored miniature lattices j0..j0+n-1 of tower u (-3 or 5) from scratch and compares them with the stored records")
 led=M.compute_ledger(track=track,use_db=True,c_tail_z6=328.0)
 pq=M.load_data()["pool_quad"]
 pq_bins=pq["bins"];np_bins=pq["nonprim_bins"]
 pool_all_active_ok=(sum(x[2]for x in pq_bins)==pq["n_prim_in"]and sum(x[1]for x in np_bins)==pq["n_nonprim_in"]and pq["n_in"]==2000000 and pq["n_prim_in"]+pq["n_nonprim_in"]==2000000 and nonprim_hnf["ok"])
 row_checks={}
 all_rows_ok=True
 for k,dv in decl["rows"].items():
  cv=led["rows"][k]
  ok=(cv<=float(dv)+1e-12)
  row_checks[k]=dict(computed=round(cv,6),declared=dv,ok=ok)
  all_rows_ok&=ok
 ladder8_ok=(filt_live["ok"]and led["la"]["ladder_all4_monotonic"]and led["la"]["all_runs_27_below_4caps"]and led["la"]["filtsim_cur27"]["n"]>=4 and
 filt_live["u14"]["N_sw_U"]>led["la"]["filtsim_cur20"]["N_sw_U"]and
 filt_live["u14"]["N_l_U"]>led["la"]["filtsim_cur20"]["N_l_U"]and
 filt_live["u14"]["k5_exp"]>led["la"]["filtsim_cur20"]["k5_exp"]and
 filt_live["u14"]["tab_r"]<M.TAB_R_CAP
 )
 fp_cnt,fp_segs=V20L.fp_prices_from_kernels(str(HERE/"kernels_v32.s"))
 fp_kernel_prices=dict(counted_k0_to_k10=fp_cnt,charged_FP_CL=V20L.FP_CL,root_k11_charged=V20L.FP_CL[11],match=(fp_cnt==V20L.FP_CL[:11]and V20L.FP_CL[11]>=max(fp_cnt[2:])))
 wheel_cnt=V20L.wheel_prices_from_kernels(str(HERE/"kernels_v32.s"))
 wheel_kernel_prices=dict(counted=wheel_cnt,charged=V20L.WHEEL_OPS,match=(wheel_cnt==V20L.WHEEL_OPS))
 checks=dict(fp_node_prices_counted_from_kernels_ok=fp_kernel_prices["match"],wheel_prices_counted_from_kernels_ok=wheel_kernel_prices["match"],algebra_sigma_w_and_crt_ok=(alg["ok"]and alg["ntt3_shoup_butterfly_ok"]and alg["ntt3_garner_crt_ok"]and alg["ntt3_r4_r3_shoup_conv_ok"]),exhaustive_fd_pure_int_bijection_ok=(fd["ok"]and fd["zero_float_geom"]and fd["sigma_inv_exact"]and fd["max_order_dedekind_ok"]and fd["exact_frac_enum_ok"]and fd["certified_interval_ok"]and fd["global_hexagon_proof_ok"]and fd["unit_index1_proved"]),
 swar_sb_u7_pure_int_and_stageB_ok=(proto["ok"]and proto["det_occ_bound_ok"]and proto["p14_zero_gcd_gm_ok"]and proto["zero_float_geom"]and proto["flat_array_zero_dict_build"]and proto["flat_sbu7_exact_ok"]and proto["zero_sort_exact_ok"]and proto["exact_slot_cutoff_M_t"]and not proto["qa_stored_in_record"]),
 lattice_b2_score_nodes_and_srt_ok=(lat_b2["ok"]and lat_b2["smat_cyclic_a345_ok"]and lat_b2["zero_float_fp_enum"]and lat_b2["strict_u64_i64_umulh_checked"]and lat_b2["zero_exception_64bit_q38_proved"]and lat_b2["cert_min_gsym_exact"]==753729 and lat_b2["isqrt_u40_lut7_errs"]==0 and lat_b2["isqrt_u40_lut7_max_over_pct"]<=0.12 and lat_b2["staddh_16b_semantic_errs"]==0 and lat_b2["b2_x4_prefix_staddh_errs"]==0 and lat_b2["lattice_b2_mismatches"]==0 and lat_b2["complete_root_to_leaf_leaves"]>=10000 and lat_b2["complete_root_to_leaf_mismatches"]==0 and lat_b2["complete_root_to_leaf_apply_errs"]==0 and lat_b2["complete_root_to_leaf_div_errs"]==0 and lat_b2["complete_root_to_leaf_max_sa_cones"]<=4 and lat_b2["wheel_tree_state_errs"]==0 and lat_b2["wheel_tree_leaves_checked"]>=10000 and lat_b2["swar21_cur_errs"]==0 and lat_b2["complete_scaled_max_rel_err_pct"]<=5.0 and lat_b2["complete_scaled_fp_ir_max_rel_err_pct"]<=5.0 and lat_b2["prod_exec_leaves"]>=2000 and lat_b2["prod_exec_mismatches"]==0 and lat_b2["prod_exec_apply_errs"]==0 and lat_b2["prod_exec_div_errs"]==0 and lat_b2["swar_tie_canon_errs"]==0 and lat_b2["barrett_kta_18op_errs"]==0 and lat_b2["barrett_kta_18op_tested"]>=3000 and lat_b2["desc_ai_prefix_chain_errs"]==0 and lat_b2["score_node_replay_errs"]==0 and lat_b2["fabric_1hop_addr_errs"]==0 and lat_b2["b_tab_packing_errs"]==0 and lat_b2["discrete_chord_bound_test_errs"]==0 and lat_b2["fabric_1hop_list_node_errs"]==0 and lat_b2["fabric_1hop_per_node_bytearray_zero_dict"]and lat_b2["fabric_1hop_zero_frag_exact_ok"]and lat_b2["greedy_node_bound_errs"]==0 and lat_b2["greedy_bound_nodes"]==16 and lat_b2["logp_bands_1_16b"]==224 and lat_b2["lat_8b_direct_ptr_page_cap"]==2046 and lat_b2["emit_full_pages_flushed"]>=20 and lat_b2["zero_hash_table_apply"]and lat_b2["inc_c0_exact_errs"]==0 and lat_b2["prod_256orb_exact_rejection_ht"]and lat_b2["prod_256orb_uniform_random_root_choice"] and lat_b2["prod_256orb_count"]==256 and lat_b2["prod_max_inv_b2_bits"]<64 and lat_b2["prod_max_r2_q12_bits"]<40 and lat_b2["prod_max_d2_bits"]<64 and lat_b2["prod_max_ck_bits"]<63 and lat_b2["swar_bounds_ok"]and lat_b2["srt_batch_mismatches"]==0 and lat_b2["e2e_factored_smooth_sqfree"]>=10 and lat_b2["e2e_factored_1lp"]>=100 and lat_b2["gate_0lp_fn"]==0 and lat_b2["gate_0lp_fp"]==0),
 v25_upgrades_ok=(v25_upg["ok"]and v25_upg["b1_e2e"]["ok"] and v25_upg["lat_apply_e2e"]["ok"] and v25_upg["a0_chord_kta_errs"]==0 and v25_upg["a0_chord_ai_order_errs"]==0 and v25_upg["cone3_s35_alg_ok"]and v25_upg["cone3_bijection_errs"]==0 and v25_upg["fd_rem0_u40_errs"]==0 and v25_upg["tbnz14_single_slot_errs"]==0 and v25_upg["sh17_x4_var_ks_errs"]==0 and v25_upg["b0_floor_div_nest_errs"]==0 and v25_upg["row_capacity_guard_errs"]==0 and v25_upg["dummy_sink_disjoint_ok"]and v25_upg["init_scan8_readonly_errs"]==0 and v25_upg["v28_wheel_ok"]),
 v32_upgrades_ok=(v32_upg["ok"]and v32_upg["adic3_exact"]["ok"]and v32_upg["ntt3_split_proto_errs"]==0 and v32_upg["ctail_lazy_errs"]==0 and v32_upg["bfi_b2_errs"]==0 and v32_upg["delta_b1_errs"]==0 and v32_upg["u16_range_ok"]and v32_upg["sf_uses_survivor_upper_ok"]and v32_upg["scan8_u16_errs"]==0 and v32_upg["b1_delta_ops_ok"]and v32_upg["boot_df_tree_ok"]and v32_upg["srt_g_leaf_ok"]and v32_upg["mul_kar_sb_ok"]and v32_upg["step1_tab_8192_inv_errs"]==0 and v32_upg["rho_tail_ok"]and v32_upg["ret2d_per_prime_errs"]==0 and v32_upg["ret2d_mc_ok"]),
 live_pool_hnf_subsample_ok=pool_sub["ok"],live_genuine_c6_rel_filter_ok=(gen_filt["ok"]and gen_filt["post_merge_cols_lost"]==0 and gen_filt["merge_backsub_errors"]==0),live_bareiss_12d_lll_limb_ok=(bareiss_lll["ok"]and bareiss_lll["inductive_sr_inv_ok"]and bareiss_lll["input_uniform_wc_proved"]and bareiss_lll["pair_setup_p59_proved"]and bareiss_lll["lat_b2_p61_proved"]),live_mont5_cios_and_suyama12_ecm_ok=(mont5_ecm["ok"]and mont5_ecm["ecm_ladder6_sep_ok"]and mont5_ecm["ecm_ladder6_min_sep_dag"]>=0.60 and mont5_ecm["ecm_ladder6_min_sep_boot"]>=0.60),
 live_z6_la_5limb_barrett_ok=(la_kern["ok"]and la_kern["unit4_v4_expand_ok"]and la_kern["unit4_8vec_ok"]and la_kern["unit4_8vec_u64_idx_43b_ok"]and la_kern["unit4_8vec_realistic_ir"]<=M.C_UNIT_5),
 c6_pgf_usable_depth2_and_hea_campaign_ok=(pgf_camp["ok"]and pgf_camp["exact_ideal_primitivity_checked"]and pgf_camp["live_4_depth2_2hop_chains_exact"]and pgf_camp["unified_bareiss_910_used"]),
 h3_scaling_and_nu_h_ok=(h3_scale["ok"]and h3_scale["fixed_instance_deterministic_h3"]and h3_scale["zero_loss_filtsim_ok"]and h3_scale["all4_monotone_20_27"]and h3_scale["hist_v4_exact_ok"]and h3_scale["fg_dens_ok"]and h3_scale["u_ram_all_lt_2p30"]),
 split_crt_norms_and_mp_rt_ok=(crt_mp["ok"]and crt_mp["chained_1pass_b2_errs"]==0),
 pool_all_active_200k_ok=pool_all_active_ok,live_filtsim_z6_and_wmid_ok=(filt_live["ok"]and filt_live["u12"]["post_merge_cols_lost"]==0 and filt_live["u12"]["merge_backsub_errors"]==0 and filt_live["u13"]["post_merge_cols_lost"]==0),filtsim_8scale_ladder_ok=ladder8_ok,unplanted_5sm_bw_dlp_and_schur_ok=hx["ok"],
 live_62row_pure_k1_boole_wald_dag_ok=(led["descent"]["ok"]and led["descent"]["nodes"]==62
 and led["descent"]["mode"]=="pure_k1_boole_wald_iid_dp_boot"
 and led["descent"]["p0_floor"]==0.001
 and led["descent"]["phi_dag_one"]<=0.50
 and led["descent"]["e_runs_per_log"]<=2.00
 and led["descent"]["mu_min"]>=500.0
 and led["descent"]["p0_max"]<=0.005
 and led["descent"]["p_node_abort"]==0.0
 and led["descent"]["p_repeat_loss"]==0.0
 and led["descent"]["p_descent_fail"]==0.0
 and led["descent"]["desc_cap_log2"]==112.0
 ),rows_ok=all_rows_ok,exp_run_ok=(led["exp_run"]<=float(decl["exp_run"])),success_ok=(led["succ_exact"]>=float(decl["success"])),gate_K_ok=(led["K_gate"]==decl["gate_K"]and led["gate_passes"]and led["gate_prev_fails"]),claim_ok=(led["claim_exact"]<=float(decl["claim"])and int(led["score_numerator"])<=int(decl["score_numerator"])),
 ram_ok=(led["geometry"]["ram_log2"]<=97.910 and led["geometry"]["ram_log2"]<100.0
 and led["geometry"]["peak_global_ram_cap_log2"]<=97.910
 and led["geometry"]["stor_bound"]["ok"]
 and led["geometry"]["n_sc_nodes_exp"]==34 and led["geometry"]["cs_blocks_max"]<(1<<26)
 and led["fail_terms"]["storage_cap_abort"]==0.0
 and led["geometry"]["list_1node_ram_log2"]<=60.176
 and led["geometry"]["emit_1own_ram_log2"]<=63.723
 and led["geometry"]["score_1node_ram_log2"]<=63.391
 and led["geometry"]["score_arena_log2"]<64.0
 and led["geometry"]["b2_1p_ram_log2"]<=57.58
 and led["geometry"]["b2_build_1p_ram_log2"]<=52.51
 and led["la"]["mem_1own_log2"]<=63.420 and led["la"]["fits_1own"]
 and led["descent"]["boot_root_retained_bytes_log2"]<=60.470 and abs(led["descent"]["boot_root_retained_bytes_log2"]-round(math.log2(M.pbits(0,led["descent"]["boot_x_ea"],1.0)/8.0),4))<1e-4
 and led["descent"]["boot_root_bytes_log2"]<=63.470 and abs(led["descent"]["boot_root_bytes_log2"]-round(math.log2(8.0*M.pbits(0,led["descent"]["boot_x_ea"],1.0)/8.0),4))<1e-4
 and led["geometry"]["stor_bound"]["checks"]["machine_lt_2_100"]
 and max(led["geometry"]["list_1node_ram_log2"],led["geometry"]["emit_1own_ram_log2"],led["geometry"]["score_1node_ram_log2"],led["geometry"]["b2_1p_ram_log2"],led["geometry"]["b2_build_1p_ram_log2"],led["la"]["mem_1own_log2"],led["descent"]["boot_root_bytes_log2"])<64.0),geo_consistency_ok=(abs(led["pool"]["frac_reg_analytic"]-led["pool"]["frac_reg"])<0.005
 and led["geometry"]["max_sb"]<=8 and led["geometry"]["sup_nf_bits"]<448.0 and led["geometry"]["sup_ng_bits"]<=512.0 and led["geometry"]["log2_Na_shell"]<47.0),manifest_subbounds_ok=(led["geometry"]["log2_V"]<=95.0920 and led["geometry"]["Ra"]<=575.55 and led["geometry"]["Rb_GA"]<=764.40
 and led["per_cell"]["rep_lo"]<=1.975 and led["per_cell"]["occ_sb"]<=0.179 and led["per_cell"]["hi_ap"]<=0.920
 and led["per_cell"]["lat_b2"]<=5.098 and led["per_cell"]["lat_lll"]<=0.00356 and led["per_cell"]["scan_zero"]<=0.426
 and led["per_cell"]["cof"]<=0.611 and led["per_cell"]["total"]<=9.209),peel_ok=(not led["la"]["peel"]["stalls"]and led["la"]["peel"]["exc_crit"]<0.75),
 )
 bdd=led["descent"];_r4=dict(led["rows"]);_r4["two_logs"]=bdd["boot_two_log_if_charge4"]
 bdd["claim_if_boot_charge_4"]=round(math.log2(1.0+led["K_gate"]/2.0**64*2.0**math.log2(sum(2.0**v for v in _r4.values()))),6)
 _r4["two_logs"]=bdd["boot_two_log_if_survivors_x4"]
 bdd["claim_if_boot_survivors_x4"]=round(math.log2(1.0+led["K_gate"]/2.0**64*2.0**math.log2(sum(2.0**v for v in _r4.values()))),6)
 checks.update(ht_pool_replay_and_pooled_bound_ok=ht_pool["ok"],boot_dp_random_lattice_replay_ok=bd_rep["ok"],
 boot_kept_rule_ledger_ok=(bdd["boot_dp_lattices"]==256 and bdd["boot_charge"]==2.0 and bdd["boot_dp_surv_monotone"]and bdd["boot_dp_pair_loss"]<1e-3
 and bdd["boot_dp_first_moment_factor"]-3.0*bdd["boot_dp_first_moment_factor_se"]>=1.0 and bdd["boot_dp_first_moment_exact_nl"]-3.0*bdd["boot_dp_first_moment_exact_nl_se"]>=1.0 and bdd["boot_dp_surv_exact_nl"]+3.0*bdd["boot_dp_surv_exact_nl_se"]<=0.5*bdd["boot_dp_surv_margin"] and bdd["boot_dp_surv_margin"]>=1.0),
 boot_mini_ratios_above_charge_ok=all(t["exhaustive"]and t["tower_ok"]and all(v["first_moment_ci95"][0]>=0.5 and v["lattice_ci95"][0]>=0.5 for v in t["by_L"].values())for t in bm_st.values()))
 checks_ok=all(checks.values())
 out=dict(ht_pool=ht_pool,boot_dp_replay=bd_rep,boot_dp_small_primes=bd_a,boot_mini_replay=bm_rep,boot_mini_stats=bm_st,track=track,part=part,checks_ok=checks_ok,all_ok=checks_ok,checks=checks,algebra=alg,exhaustive_fd=fd,swar_sb_u7_stageB=proto,v20_lat_b2_and_srt=lat_b2,fp_kernel_prices=fp_kernel_prices,wheel_kernel_prices=wheel_kernel_prices,v25_upgrades=v25_upg,v32_upgrades=v32_upg,live_pool_hnf_subsample=pool_sub,live_genuine_c6_rel_filter=gen_filt,
 live_bareiss_12d_lll=bareiss_lll,live_mont5_cios_and_suyama12_ecm=mont5_ecm,live_z6_la_kernels=la_kern,nonprim_hnf_witnesses=nonprim_hnf,
 c6_pgf_campaign=pgf_camp,diagnostics_part=dict(status="not_run",checks=["descent_12d_lattices_24sq_24boot_ok","deg8_dispersion_and_prod_yield_ok","c6_adaptive_h5_live_ok"],how_to_run={"part":"diagnostics"},note="H5 scoped diagnostics (ADAPTIVE.md section 2; v23b code, the boot-lattice part reduced to gating all 144 candidates) are not run in the default run because of the board's per-run time limit; the cost bound is conditional on declared heuristic H5; one run with part=diagnostics executes all three and reports their observations"),h3_scaling_and_nu_h=h3_scale,split_crt_norms_and_mp_rt=crt_mp,live_filtsim_z6=filt_live,hx_unplanted_5sm_bw=hx,ledger=dict(rows=row_checks,exp_run=round(led["exp_run"],9),exp_run_decl=decl["exp_run"],success=round(led["succ_exact"],8),success_decl=decl["success"],gate_K=led["K_gate"],claim=round(led["claim_exact"],9),claim_decl=decl["claim"],
 score_numerator=decl["score_numerator"],fail_terms=led["fail_terms"],),geometry=led["geometry"],per_cell=led["per_cell"],la=led["la"],descent=led["descent"],)
 if part in("sensitivity","full"):
  led_plain=M.compute_ledger(track=track,use_db=False)
  led_lat55=M.compute_ledger(track=track,use_db=True,strict_per_owner=True,swar_probe=10.0)
  led_flat=M.compute_ledger(track=track,use_db=True,strict_per_owner=False)
  led_tight=M.compute_ledger(track=track,use_db=True,cap_scale=0.80)
  led_srt30=M.compute_ledger(track=track,use_db=True,srt_g_cap=0.30)
  out["sensitivity"]=dict(plain_dickman_rho=dict(exp_run=round(led_plain["exp_run"],6),claim=round(led_plain["claim_exact"],6)),conservative_lat_b2_plus10=dict(exp_run=round(led_lat55["exp_run"],6),claim=round(led_lat55["claim_exact"],6)),conservative_srt_030_cap=dict(exp_run=round(led_srt30["exp_run"],6),claim=round(led_srt30["claim_exact"],6)),flat_shared_arena_ref=dict(exp_run=round(led_flat["exp_run"],6),claim=round(led_flat["claim_exact"],6)),
  de_bruijn_tighter_cap=dict(exp_run=round(led_tight["exp_run"],6),claim=round(led_tight["claim_exact"],6)),)
 out["elapsed_s"]=round(time.time()-t0,3)
 print(json.dumps(out))
 if not checks_ok:
  raise SystemExit(1)
if __name__=="__main__":
 main()

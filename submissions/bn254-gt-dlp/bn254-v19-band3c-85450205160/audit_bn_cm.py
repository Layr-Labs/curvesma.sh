"""audit_bn_cm.py — BN254 C6 v19-band3c live verifier (30/30 checks)."""
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
import band3_lbs as B3
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
 pool_sub=V.verify_live_c6_pool_hnf_subsample(n_tot=1000,seed=20261120)
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
 led=M.compute_ledger(track=track,use_db=True,c_tail_z6=740.0)
 b3=B3.verify_band3_tau_lbs(led=led,seed=20261121)
 pq=M.load_data()["pool_quad"]
 pq_bins=pq["bins"];np_bins=pq["nonprim_bins"]
 pool_all_active_ok=(sum(x[2]for x in pq_bins)==136563 and sum(x[1]for x in np_bins)==63437 and pq["n_in"]==200000 and nonprim_hnf["ok"])
 # --- Relation-conditioned Band-3 loss (organizer objection): independent live recomputation ---
 # loss = p_fb * ln(p_mid/p_list) * sup_{p in [2^p_list,2^p_mid), n in supp(mix)} p * P(p | usable g-norm)
 # with p*P = Num(n-p)/Num(n) = [D(n-p)/D(n)] * [at(n-p)/at(n)] * exp(tilt*p) from the same
 # Retention model that produces r_g. The mixture is the relation-weighted g-norm size distribution
 # (35 discrete norm sizes); the charged factor is the supremum over those sizes and the band, not
 # an average over the mixture.
 st3=0.25;RET=M.Retention(M.OP["lg"],M.OP["lx2"],prime_min=2**M.OP["p_min"],explicit_bits=min(M.OP["p_min"],14),step=st3,max_bits=530,tilt=0.05,odd=1)
 p_list_b3,pm_b3,sl_b3=M.OP["p_list"],M.OP["p_mid"],M.OP["slack"]
 d_b3=(led["geometry"]["log2_V"]-(pq["log2_V0"]+M.log2(pq["frac_reg"])))/12.0
 af_odd=M.ALPHA_F_RAW+2.0/3.0;ag_odd=M.ALPHA_G_RAW+1.0/3.0
 cl2=((7.0/9.0,0.0,0.0),(1.0/9.0,0.0,3.0),(1.0/9.0,6.0,0.0))
 loss_b3_led=led["geometry"]["b3_rel_loss"]
 def _rg_chk(loss_v):
  s0=s1=0.0
  for lnf0,lng0,cnt_prim,cnt_all in pq_bins:
   if not cnt_prim:continue
   nf_raw=lnf0+48.0*d_b3;ng_raw=lng0+12.0*d_b3;ng_q=round(ng_raw*2.0)/2.0
   a0=RET.at(ng_q,sl_b3);a3=RET.at(ng_q-3.0,sl_b3)
   y=rg=0.0
   for w2,v2f,v2g in cl2:
    nf_c=nf_raw-v2f+af_odd;ng_c=ng_raw-v2g+ag_odd;ac=(a3 if v2g else a0)*(1.0-loss_v)
    yc=w2*M.rho_db(nf_c/M.OP["lf"],nf_c)*M.rho_db(ng_c/M.OP["lg"],ng_c)*ac
    y+=yc;rg+=yc*ac
   s0+=cnt_prim*y;s1+=cnt_prim*y*(rg/y if y>0.0 else 0.0)
  return s0,s1
 s0n,s1n=_rg_chk(0.0);s0l,s1l=_rg_chk(loss_b3_led)
 r_g0=s1n/s0n;r_g_chk=s1l/s0l
 mix_g={};wtot=0.0
 for lnf0,lng0,cnt_prim,cnt_all in pq_bins:
  if not cnt_prim:continue
  nf_raw=lnf0+48.0*d_b3;ng_raw=lng0+12.0*d_b3;ng_q=round(ng_raw*2.0)/2.0
  a0=RET.at(ng_q,sl_b3);a3=RET.at(ng_q-3.0,sl_b3)
  for w2,v2f,v2g in cl2:
   nf_c=nf_raw-v2f+af_odd;ng_c=ng_raw-v2g+ag_odd
   wgt=cnt_prim*w2*M.rho_db(nf_c/M.OP["lf"],nf_c)*M.rho_db(ng_c/M.OP["lg"],ng_c)*(a3 if v2g else a0)
   ns=round(ng_c/st3);mix_g[ns]=mix_g.get(ns,0.0)+wgt;wtot+=wgt
 phi_b3_chk,b_b3_chk=0.0,0.0
 for i in range(13):
  b=p_list_b3+(pm_b3-p_list_b3)*i/12.0;kb=round(b/st3)
  for ns in mix_g:
   if ns-kb<1:continue
   phi=M.exp(RET.tilt*b)*(RET.denominator[ns-kb]/RET.denominator[ns])*(RET.at((ns-kb)*st3,sl_b3)/RET.at(ns*st3,sl_b3))
   if phi>phi_b3_chk:phi_b3_chk,b_b3_chk=phi,b
 loss_b3_chk=phi_b3_chk*M.P_SKIP*M.log(pm_b3/p_list_b3)
 ng_mean_b3=sum(ns*st3*w for ns,w in mix_g.items())/wtot
 mk1=M.band3_marked_retention(RET,ng_mean_b3,sl_b3,1.0,p_list_b3,14.0)
 mkp=M.band3_marked_retention(RET,ng_mean_b3,sl_b3,1.0-M.P_SKIP,p_list_b3,14.0)/mk1
 row_checks={}
 all_rows_ok=True
 for k,dv in decl["rows"].items():
  cv=led["rows"][k]
  ok=(cv<=float(dv)+1e-12)
  row_checks[k]=dict(computed=round(cv,6),declared=dv,ok=ok)
  all_rows_ok&=ok
 ladder5_ok=(filt_live["ok"]and led["la"]["ladder_all4_monotonic"]and led["la"]["all_runs_below_4caps"]and
 filt_live["u14"]["N_sw_U"]>led["la"]["filtsim_cur20"]["N_sw_U"]and
 filt_live["u14"]["N_l_U"]>led["la"]["filtsim_cur20"]["N_l_U"]and
 filt_live["u14"]["k5_exp"]>led["la"]["filtsim_cur20"]["k5_exp"]and
 filt_live["u14"]["tab_r"]>led["la"]["filtsim_cur20"]["tab_r"]
 )
 checks=dict(algebra_sigma_w_and_crt_ok=(alg["ok"]and alg["ntt3_shoup_butterfly_ok"]and alg["ntt3_garner_crt_ok"]and alg["ntt3_r4_r3_shoup_conv_ok"]),exhaustive_fd_pure_int_bijection_ok=(fd["ok"]and fd["zero_float_geom"]and fd["sigma_inv_exact"]and fd["max_order_dedekind_ok"]and fd["exact_frac_enum_ok"]and fd["certified_interval_ok"]and fd["global_hexagon_proof_ok"]and fd["unit_index1_proved"]),
 swar_sb_u7_pure_int_and_stageB_ok=(proto["ok"]and proto["det_occ_bound_ok"]and proto["p14_zero_gcd_gm_ok"]and proto["zero_float_geom"]and proto["flat_array_zero_dict_build"]and proto["flat_sbu7_exact_ok"]and proto["zero_sort_exact_ok"]and proto["exact_slot_cutoff_M_t"]and not proto["qa_stored_in_record"]),
 live_pool_hnf_subsample_ok=pool_sub["ok"],live_genuine_c6_rel_filter_ok=gen_filt["ok"],live_bareiss_12d_lll_limb_ok=(bareiss_lll["ok"]and bareiss_lll["inductive_sr_inv_ok"]and bareiss_lll["input_uniform_wc_proved"]and bareiss_lll["pair_setup_p59_proved"]),live_mont5_cios_and_suyama12_ecm_ok=(mont5_ecm["ok"]and mont5_ecm["ecm_ladder6_sep_ok"]and mont5_ecm["ecm_ladder6_min_sep_dag"]>=0.6174 and mont5_ecm["ecm_ladder6_min_sep_boot"]>=0.6236),
 live_z6_la_5limb_barrett_ok=(la_kern["ok"]and la_kern["unit4_v4_expand_ok"]),descent_12d_lattices_24sq_24boot_ok=(lat12["ok"]and lat12["prenorm_gate_uniform_ok"]and lat12["all_24_sq_reduced_live_by_bareiss_910"]and lat12["all_24_boot_bareiss_910_certified"]and lat12["boot_24_b_best_prim_ok"]),
 c6_pgf_usable_depth2_and_hea_campaign_ok=(pgf_camp["ok"]and pgf_camp["exact_ideal_primitivity_checked"]and pgf_camp["live_4_depth2_2hop_chains_exact"]and pgf_camp["unified_bareiss_910_used"]),
 deg8_dispersion_and_prod_yield_ok=(deg8_camp["ok"]and deg8_camp["tower_u5_p27631_r27481_deg8_frec_verified"]and deg8_camp["all_100_above_fb_1024"]and deg8_camp["live_6_deg8_ideals_and_2hop_d2_exact"]),
 c6_adaptive_h5_live_ok=(adap_h5["ok"]and adap_h5["strict_k1_zero_retry_protocol"]and adap_h5["live_k1_zero_retry_ok"]and adap_h5["live_witnesses_replayed"]==3 and adap_h5["after_fail_total"]==1856 and adap_h5["pooled_oe_fail"]<0.85),
 h3_scaling_and_nu_h_ok=(h3_scale["ok"]and h3_scale["fixed_instance_deterministic_h3"]and h3_scale["all4_monotone_20_24"]and h3_scale["hist_v4_exact_ok"]and h3_scale["u_ram_all_lt_2p30"]),
 split_crt_norms_and_mp_rt_ok=(crt_mp["ok"]and crt_mp["chained_1pass_b2_errs"]==0),
 pool_all_active_89200_ok=pool_all_active_ok,live_filtsim_z6_and_wmid_ok=filt_live["ok"],filtsim_5scale_ladder_ok=ladder5_ok,unplanted_5sm_bw_dlp_and_schur_ok=hx["ok"],
 live_62row_k1_boole_wald_dag_ok=(led["descent"]["ok"]and led["descent"]["nodes"]==62
 and led["descent"]["mode"]=="pure_k1_boole_wald_iid_boot"
 and led["descent"]["nb_r"]==1.0
 and led["descent"]["nb_alpha"]==0.5
 and led["descent"]["p0_floor"]==0.001
 and led["descent"]["p_boot_tree"]<=0.443500
 and led["descent"]["e_runs_per_log"]<=1.796900
 and led["descent"]["e_trials_per_log"]<=6.56e15
 and abs(led["descent"]["b3_skip"]-led["geometry"]["b3_rel_loss"])<1e-15
 and led["descent"]["mu_min"]>=725.46
 and led["descent"]["p0_max"]<=2.750e-3
 and led["descent"]["c_det_one"]==427.0
 and led["descent"]["p_node_abort"]==0.0
 and led["descent"]["p_repeat_loss"]==0.0
 and led["descent"]["p_descent_fail"]==0.0
 and led["descent"]["desc_cap_log2"]==112.0
 ),rows_ok=all_rows_ok,exp_run_ok=(led["exp_run"]<=float(decl["exp_run"])),success_ok=(led["succ_exact"]>=float(decl["success"])),gate_K_ok=(led["K_gate"]==decl["gate_K"]and led["gate_passes"]and led["gate_prev_fails"]),claim_ok=(led["claim_exact"]<=float(decl["claim"])and int(led["score_numerator"])<=int(decl["score_numerator"])),
 ram_ok=(led["geometry"]["ram_log2"]<100.0 and led["geometry"]["score_arena_log2"]<64.0
 and led["geometry"]["b0_1own_ram_log2"]<64.0
 and led["geometry"]["b2_1p_ram_log2"]<64.0 and led["geometry"]["b2_build_1p_ram_log2"]<64.0
 and led["geometry"]["b3_p1_1own_ram_log2"]<64.0 and led["geometry"]["b3_p2_1own_ram_log2"]<64.0
 and led["la"]["fits_1own"]and led["descent"]["boot_root_bytes_log2"]<64.0),geo_consistency_ok=(abs(led["pool"]["frac_reg_analytic"]-led["pool"]["frac_reg"])<0.005
 and led["geometry"]["max_sb"]<=8 and led["geometry"]["sup_nf_bits"]<448.0),peel_ok=(not led["la"]["peel"]["stalls"]and led["la"]["peel"]["exc_crit"]<0.75),band3_tau_lbs_ok=b3["ok"],
 band3_relation_loss_ok=(abs(loss_b3_chk-float(decl["b3_rel_loss"]))<1e-12 and abs(loss_b3_chk-loss_b3_led)<1e-12
 and abs(phi_b3_chk-float(decl["b3_phi_max"]))<1e-12 and abs(b_b3_chk-pm_b3)<1e-9 and phi_b3_chk>0.0
 and abs(mk1-RET.at(ng_mean_b3,sl_b3))<1e-12 and mkp>=1.0-loss_b3_chk and loss_b3_chk>M.P_SKIP),
 band3_retention_wiring_ok=(abs(r_g_chk-led["geometry"]["r_g"])<1e-12 and abs(r_g_chk-(1.0-loss_b3_led)*r_g0)<1e-12
 and r_g_chk<r_g0 and led["descent"]["b3_skip"]==loss_b3_led and led["geometry"]["b3_mix_sizes"]>=10),
 )
 checks_ok=all(checks.values())
 out=dict(track=track,part=part,checks_ok=checks_ok,all_ok=checks_ok,checks=checks,algebra=alg,exhaustive_fd=fd,swar_sb_u7_stageB=proto,live_pool_hnf_subsample=pool_sub,live_genuine_c6_rel_filter=gen_filt,
 live_bareiss_12d_lll=bareiss_lll,live_mont5_cios_and_suyama12_ecm=mont5_ecm,live_z6_la_kernels=la_kern,nonprim_hnf_witnesses=nonprim_hnf,descent_12d_lattices=lat12,
 c6_pgf_campaign=pgf_camp,deg8_dispersion_and_prod_yield=deg8_camp,c6_adaptive_h5_live=adap_h5,h3_scaling_and_nu_h=h3_scale,split_crt_norms_and_mp_rt=crt_mp,live_filtsim_z6=filt_live,hx_unplanted_5sm_bw=hx,band3_tau_lbs=b3,
 band3_relation_loss=dict(loss=round(loss_b3_chk,15),loss_src='phi_max*p_fb*ln(p_mid/p_list)',p_fb=M.P_SKIP,H=round(M.log(pm_b3/p_list_b3),15),
  phi_max=round(phi_b3_chk,12),phi_sup_at_b=round(b_b3_chk,4),ng_mixture_mean=round(ng_mean_b3,4),norm_sizes=len(mix_g),
  marked_at_z1=round(mk1,15),model_at=round(RET.at(ng_mean_b3,sl_b3),15),indep_thinning_retention=round(mkp,12),
  union_bound_loss=round(loss_b3_chk,15),indep_loss=round(1.0-mkp,12),
  r_g_no_loss=round(r_g0,12),r_g_charged=round(r_g_chk,12),r_g_ledger=round(led["geometry"]["r_g"],12),descent_b3_skip=led["descent"]["b3_skip"]),
 ledger=dict(rows=row_checks,exp_run=round(led["exp_run"],9),exp_run_decl=decl["exp_run"],success=round(led["succ_exact"],8),success_decl=decl["success"],gate_K=led["K_gate"],claim=round(led["claim_exact"],9),claim_decl=decl["claim"],
 score_numerator=decl["score_numerator"],fail_terms=led["fail_terms"],),geometry=led["geometry"],per_cell=led["per_cell"],la=led["la"],descent=led["descent"],)
 if part in("sensitivity","full"):
  led_plain=M.compute_ledger(track=track,use_db=False)
  led_swar22=M.compute_ledger(track=track,use_db=True,strict_per_owner=True,swar_probe=22.0)
  led_flat=M.compute_ledger(track=track,use_db=True,strict_per_owner=False)
  led_tight=M.compute_ledger(track=track,use_db=True,cap_scale=0.80)
  led_b3_85=M.compute_ledger(track=track,use_db=True,c_lbs=85.0)
  out["sensitivity"]=dict(plain_dickman_rho=dict(exp_run=round(led_plain["exp_run"],6),claim=round(led_plain["claim_exact"],6)),conservative_swar22=dict(exp_run=round(led_swar22["exp_run"],6),claim=round(led_swar22["claim_exact"],6)),flat_shared_arena_ref=dict(exp_run=round(led_flat["exp_run"],6),claim=round(led_flat["claim_exact"],6)),
  de_bruijn_tighter_cap=dict(exp_run=round(led_tight["exp_run"],6),claim=round(led_tight["claim_exact"],6)),band3_c_lbs_85=dict(exp_run=round(led_b3_85["exp_run"],6),claim=round(led_b3_85["claim_exact"],6)),)
 out["elapsed_s"]=round(time.time()-t0,3)
 print(json.dumps(out))
 if not checks_ok:
  raise SystemExit(1)
if __name__=="__main__":
 main()

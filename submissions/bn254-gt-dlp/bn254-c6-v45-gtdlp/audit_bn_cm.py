BH='sensitivity'
BG='boot_dp_surv_margin'
BF='f_node_log2'
BE='boot_x_ea'
BD='boot_root_retained_bytes_log2'
BC='mem_1own_log2'
BB='b2_build_1p_ram_log2'
BA='b2_1p_ram_log2'
B9='score_1node_ram_log2'
B8='emit_1own_ram_log2'
B7='list_1node_ram_log2'
B6='fail_terms'
B5='ram_log2'
B4='success'
B3='succ_exact'
B2='merge_backsub_errors'
B1='k5_exp'
B0='N_sw_U'
A_='n_nonprim_in'
Az='n_prim_in'
Ay='diagnostics'
Ax='boot_mini'
Aw='ht_run'
Av='zero_float_geom'
Au='all_off'
At='descbk'
As='lat_lll'
Ar='lat_b2'
Aq='phi_dag_one'
Ap='declared.json'
B='files'
Ah='checks'
Ag='boot_root_bytes_log2'
Af='score_numerator'
Ae='post_merge_cols_lost'
Ad='filtsim_cur20'
Ac='plan59'
Ab=range
AP="kernel-emulation and boot-replay checks; not in the default run because of the board's per-run time limit"
AO='K_gate'
AN='boot_dp'
AM='kernels_v32.s'
AL='two_logs'
AK=bool
A8='u13'
A7='kernels'
A6=SystemExit
A1='not_run'
A0='all_on'
z=float
y=max
s='part'
r=abs
o=1.
n='claim'
m='rows'
j=None
d='stor_bound'
c='claim_exact'
b=all
a=str
U='per_cell'
T=print
S='exp_run'
R=True
P='la'
N=int
G='descent'
F='geometry'
E=dict
D='ok'
C=round
import json,math as l,sys,time
from pathlib import Path
A=Path(__file__).resolve().parent
i=A/B if(A/B/Ap).is_file()else A
if a(i)not in sys.path:sys.path.insert(0,a(i))
import model_bn254_v6 as M,verify_bn254_v6 as V,descent_bn254_v6 as DESC,hx_z6 as An,c6_adaptive as AZ,v20_lat_b2 as Z,v25_verify as Aa,boot_dp as A5,boot_mini as AJ,v32_verify as Bn,pt_verify as Ao,os
def main():
	W=time.time();BI=sys.stdin.read().strip()if not sys.stdin.isatty()else'';O=json.loads(BI)if BI else{};H=O.get(s,'default');Q=json.loads((i/Ap).read_text());J=Q['track']
	if H=='filtsim':AQ=N(O.get('log2U',13));Ai=N(O.get('seed',20261120));BJ=V.verify_live_filtsim_z6(log2U=min(y(AQ,10),17),seed=Ai,W_target=100,check_kernel=AQ<=13,texc_abs=O.get('texc'));T(json.dumps(E(track=J,part=H,checks_ok=BJ[D],live_filtsim_z6=BJ,elapsed_s=C(time.time()-W,3))));return
	if H=='ablation':
		AR=M.load_data()[G]
		def Bo(lla,rp,pt,kn=1,pl=1,ew=1):
			B=M.w,E(M.OP),M.C_LAT_EMIT,M.C_LAT_LEAF,M.PT_FB_EXC,M.LLA_ON,M.EXACT_W2D,AR[Ac],M.PT_ON,DESC.KNAP_ON,DESC.descbk_swap(()if rp or pt else tuple(DESC.DESCBK_TYPED))
			if not ew:M.EXACT_W2D=False;M.OP['slack']=117.
			if not pl:M.w=((40.,42.5,8,'x4'),(42.5,44.5,16,'x4'),(44.5,46.,32,'x4'));M.OP['p_mid']=46.;M.OP['p_lat']=46.;M.C_LAT_EMIT=6.714;M.C_LAT_LEAF=6.714+M.C_LAT_APPLY;M.PT_FB_EXC=.03
			M.LLA_ON=AK(lla);M.PT_ON=AK(pt);DESC.KNAP_ON=AK(kn)
			if not rp:AR[Ac]=AR['plan59_pub']
			M._DESC_CACHE.clear();M.Am.clear()
			try:A=M.compute_ledger(track=J,use_db=R,c_tail_z6=364.)
			finally:M.w,op0,M.C_LAT_EMIT,M.C_LAT_LEAF,M.PT_FB_EXC,M.LLA_ON,M.EXACT_W2D,AR[Ac],M.PT_ON,DESC.KNAP_ON,D=B;M.OP.clear();M.OP.update(op0);DESC.descbk_swap(D);M._DESC_CACHE.clear();M.Am.clear()
			return E(claim=C(A[c],9),exp_run=C(A[S],9),la=C(A[m][P],6),gen=C(A[m]['gen'],6),collection=C(A[m]['collection'],6),lists=C(A[m]['lists'],6),two_logs=C(A[m][AL],6),phi_dag_one=C(A[G][Aq],6),lat_b2=C(A[U][Ar],6),lat_lll=C(A[U][As],6))
		A9=O.get('which','a');A2,BEw,BPl,BK,BL,BM,BN,A4p,A4a,AS,AA='shipped_all_on','EXACT_W2D_off','PLAT44_off','L_LA_off','REPLAN_off','PT_off','KNAP_off','anchor_09c6d291','anchor_4ae1ef46','anchor_lla','anchor_v43b';Bp={A2:(1,1,1,1,1,1),BEw:(1,1,1,1,1,0),BPl:(1,1,1,1,0,1),BK:(0,1,1,1,1,1),BL:(1,0,1,1,1,1),BM:(1,1,0,1,1,1),BN:(1,1,1,0,1,1),A4p:(1,1,1,1,1,0),A4a:(1,1,1,1,0,0),AS:(1,0,0,0,0,0),AA:(0,0,0,0,0,0)};AT=[A4p,A4a,AS,AA]if A9=='b'else[A2,BEw,BPl,BK,BL,BM,BN];X={aa:Bo(*Bp[aa])for aa in AT};t=lambda k:X[k][n];p={}
		if A2 in X:p.update(shipped_equals_declared=r(t(A2)-z(Q[n]))<=1.5e-09,each_switch_lowers_the_claim_from_all_on=b(t(A2)<t(A)for A in X if A!=A2),exact_w2d_off_equals_09c6d291=r(t(BEw)-82.302082328)<=1.5e-09)
		if AA in X:p.update(anchor_09c6d291_ok=r(t(A4p)-82.302082328)<=1.5e-09,anchor_4ae1ef46_ok=r(t(A4a)-82.327148880)<=1.5e-09,all_off_equals_v43b_82_671766821=r(t(AA)-82.671766821)<1e-08,replan_and_pt_off_equal_lla_package_82_527471810=r(t(AS)-82.52747181)<=1.5e-09,l_la_lowers_the_claim_with_replan_and_pt_off=t(AS)<t(AA))
		T(json.dumps(E(track=J,part=H,which=A9,checks_ok=b(p.values())and len(p)>0,ablation_checks=p,ablation=X,note="Pure 1-lever ExactScore2D ablation over 09c6d291: switches M.EXACT_W2D (integer-score Retention2D over exact W_int and bstep=0.5 R_bits grid, sl=112.6, [T-0.266,T+0.266] bin splitting), PLAT44, M.LLA_ON, stored plan59_pub, M.PT_ON, DESC.KNAP_ON. Variants (input which, default a): a = shipped_all_on, EXACT_W2D_off (=09c6d291), PLAT44_off, L_LA_off, REPLAN_off, PT_off, KNAP_off; b = anchor_09c6d291 (82.302082328), anchor_4ae1ef46 (82.327148880), anchor_lla (82.527471809), anchor_v43b (82.671766821).",elapsed_s=C(time.time()-W,3))));return
	if H==At:
		Bq=E(PAIR=('PAIR_SETUP',),FIXB=('FIX_LLL',),OMEGA=('OMEGA_BITS',),IDEAL=('IDEAL',),SCANZERO=('SCAN','ZERO'),HIT=('SIEVE_HIT',));Aj={ab:getattr(DESC,ab)for ab in DESC.DESCBK_TYPED}
		def Br(off):
			B=DESC.descbk_swap(tuple(off));M._DESC_CACHE.clear()
			try:A=M.compute_ledger(track=J,use_db=R,c_tail_z6=364.)
			finally:DESC.descbk_swap(B);M._DESC_CACHE.clear()
			return E(claim=C(A[c],9),exp_run=C(A[S],9),two_logs=C(A[m][AL],6),constants_off=list(off))
		A9=O.get('which','a');AT=[A0,Au,'PAIR_off','FIXB_off','OMEGA_off']if A9=='a'else[A0,'IDEAL_off','SCANZERO_off','HIT_off'];X={ac:Br(DESC.DESCBK_TYPED if ac==Au else()if ac==A0 else Bq[ac[:-4]])for ac in AT};p=E(shipped_equals_declared=r(X[A0][n]-z(Q[n]))<=1.5e-09,each_item_off_raises_the_claim=b(X[A0][n]<X[A][n]for A in AT if A!=A0 and any(DESC.DESCBK_NEW[A]!=DESC.DESCBK_TYPED[A]for A in X[A]['constants_off'])),constants_restored_after_the_states=b(getattr(DESC,A)==Aj[A]for A in Aj));T(json.dumps(E(track=J,part=H,which=A9,checks_ok=b(p.values()),checks=p,states=X,typed_values=DESC.DESCBK_TYPED,layer_values=DESC.DESCBK_NEW,constants_in_force=Aj,note='DESCBK layer (DESC.DESCBK_*): PAIR_SETUP, FIX_LLL, OMEGA_BITS (three of its seven constants) are the values the package proves, the others keep their typed values; descbk_constants_ok checks each against the count or worst case the package computes. States: all_on, all_off (all typed), X_off. Variants (input which): a = all_on, all_off, PAIR_off, FIXB_off, OMEGA_off; b = all_on, IDEAL_off, SCANZERO_off, HIT_off',elapsed_s=C(time.time()-W,3))));return
	if H==A7:
		e=V.verify_swar_sb_u7_and_stageB();I=Bn.verify_v32_all(M,a(i/AM),seed=20261123);K=Aa.verify_v25_upgrades(seed=20261122);Bs=N.from_bytes(os.urandom(2),'big')%256;BO=A5.replay([Bs]);f=E(swar_sb_u7_pure_int_and_stageB_ok=e[D]and e['det_occ_bound_ok']and e['p14_zero_gcd_gm_ok']and e[Av]and e['flat_array_zero_dict_build']and e['flat_sbu7_exact_ok']and e['zero_sort_exact_ok']and e['exact_slot_cutoff_M_t']and not e['qa_stored_in_record'],v32_upgrades_ok=I[D]and I['adic3_exact'][D]and I['ntt3_split_proto_errs']==0 and I['ctail_lazy_errs']==0 and I['bfi_b2_errs']==0 and I['delta_b1_errs']==0 and I['u16_range_ok']and I['sf_uses_survivor_upper_ok']and I['scan8_u16_errs']==0 and I['b1_delta_ops_ok']and I['boot_df_tree_ok']and I['srt_g_leaf_ok']and I['mul_kar_sb_ok']and I['v37_ntt_ok']and I['v39_emit8_ok']and I['v40_descent_ok']and I['step1_tab_8192_inv_errs']==0 and I['rho_tail_ok']and I['ret2d_per_prime_errs']==0 and I['ret2d_mc_ok']and I['exact_w2d_ok'],v25_upgrades_ok=K[D]and K['b1_e2e'][D]and K['lat_apply_e2e'][D]and K['a0_chord_kta_errs']==0 and K['a0_chord_ai_order_errs']==0 and K['cone3_s35_alg_ok']and K['cone3_bijection_errs']==0 and K['fd_rem0_u40_errs']==0 and K['tbnz14_single_slot_errs']==0 and K['sh17_x4_var_ks_errs']==0 and K['b0_floor_div_nest_errs']==0 and K['row_capacity_guard_errs']==0 and K['dummy_sink_disjoint_ok']and K['init_scan8_readonly_errs']==0 and K['v28_wheel_ok'],boot_dp_random_lattice_replay_ok=BO[D]);Y=b(f.values());T(json.dumps(E(track=J,part=H,checks_ok=Y,all_ok=Y,checks=f,swar_sb_u7_stageB=e,v32_upgrades=I,v25_upgrades=K,boot_dp_replay=BO,elapsed_s=C(time.time()-W,3))))
		if not Y:raise A6(1)
		return
	if H=='pt':
		u=Ao.pt_full(a(i),hb=N(O.get(Aw,0))or 5);T(json.dumps(E(track=J,part=H,checks_ok=u[D],all_ok=u[D],pt=u,elapsed_s=C(time.time()-W,3))))
		if not u[D]:raise A6(1)
		return
	if H=='dag_survivors':DESC.DAG_SURV_X=2.;AU=M.compute_ledger(track=J,use_db=R,c_tail_z6=364.);T(json.dumps(E(track=J,part=H,checks_ok=R,dag_survivor_envelope_x2=E(claim=C(AU[c],6),exp_run=C(AU[S],6),two_dag=C(AU[G]['two_dag'],4),two_logs=C(AU[G]['two_log'],4)),elapsed_s=C(time.time()-W,3))));return
	if H in(AN,Ax):AB=N(O.get('u',-3));AB=AB if a(AB)in AJ.load()else-3;BP=256 if H==AN else AJ.load()[a(AB)]['n'];AC=min(y(N(O.get('j0',0)),0),BP-1);AV=min(y(N(O.get('n',4)),1),8 if H==AN else 40,BP-AC);u=A5.replay(list(Ab(AC,AC+AV)),full=R)if H==AN else AJ.replay(AB,list(Ab(AC,AC+AV)));T(json.dumps(E(track=J,part=H,checks_ok=u[D],replay=u,elapsed_s=C(time.time()-W,3))));return
	if H=='boot_evidence':
		BQ=A5.rough_med_audit(Ab(4),1<<20);BR=A5.joint_ea_mini_audit(-3,80);AD=AK(BQ[D]and BR[D]);T(json.dumps(E(track=J,part=H,checks_ok=AD,prod_rough_16_20=BQ,mini_joint_ea_m3=BR,elapsed_s=C(time.time()-W,3))))
		if not AD:raise A6(1)
		return
	if H=='hx':AQ=N(O.get('log2U',12));Ai=N(O.get('seed',1));Bt=N(O.get('m',10));AV=N(O.get('n',8));AE=An.run_hx_1d_kernel_z6(log2U=min(y(AQ,11),13),seed=Ai,m=Bt,n=AV);T(json.dumps(E(track=J,part=H,checks_ok=AE[D],hx_z6=AE,elapsed_s=C(time.time()-W,3))));return
	g=H==Ay;AF=M.verify_tower_algebra_and_crt()if not g else j;k=M.verify_exhaustive_fd(R2=5e1)if not g else j;BS=Aa.verify_nonprim_hnf_witnesses()if not g else j;BT=V.verify_live_c6_pool_hnf_subsample(n_tot=1000,seed=20261121)if not g else j;AW=V.verify_live_genuine_c6_relation_filter(seed=20261120)if not g else j;v=V.verify_live_bareiss_12d_lll(seed=20261120)if not g else j;AG=DESC.verify_live_mont5_and_suyama12_ecm(seed=20261120)if not g else j;A3=DESC.verify_live_z6_la_kernels(seed=20261120)if not g else j;A4=V.verify_12d_descent_lattices()if g else j
	if g:
		AH=V.verify_c6_pgf_and_hea_campaign();AI=DESC.verify_deg8_dispersion_and_prod_yield(seed=20261120);w=AZ.verify_c6_adaptive_h5_live();BU=V.verify_live_trim_pair();f=E(descent_12d_lattices_24sq_24boot_ok=A4[D]and A4['prenorm_gate_uniform_ok']and A4['all_24_sq_reduced_live_by_bareiss_910']and A4['all_24_boot_bareiss_910_certified']and A4['boot_candidates_gated']==3456,deg8_dispersion_and_prod_yield_ok=AI[D]and AI['tower_u5_p27631_r27481_deg8_frec_verified']and AI['all_100_above_fb_1024']and AI['live_6_deg8_ideals_and_2hop_d2_exact'],c6_adaptive_h5_live_ok=w[D]and w['strict_k1_zero_retry_protocol']and w['live_k1_zero_retry_ok']and w['live_witnesses_replayed']==3 and w['after_fail_total']==1856 and w['pooled_oe_fail']<.85,c6_pgf_usable_depth2_and_hea_campaign_ok=AH[D]and AH['exact_ideal_primitivity_checked']and AH['live_4_depth2_2hop_chains_exact']and AH['unified_bareiss_910_used'],live_trim_pair_ok=BU[D]);Y=b(f.values());T(json.dumps(E(track=J,part=H,checks_ok=Y,all_ok=Y,checks=f,descent_12d_lattices=A4,deg8_dispersion_and_prod_yield=AI,c6_adaptive_h5_live=w,c6_pgf_campaign=AH,live_trim_pair=BU,elapsed_s=C(time.time()-W,3))))
		if not Y:raise A6(1)
		return
	Bu=V.verify_trim_pairs();q=AZ.verify_h3_scaling_and_nu_h();Ak=AZ.verify_split_crt_norms_and_mp_rt(seed=20261120);h=V.verify_live_filtsim_ladder(seed=20261120);AE=An.run_hx_1d_kernel_z6(log2U=11,seed=1,m=10,n=8);Z.HT_SEED_BASE=N(O.get(Aw,0))or 1+N.from_bytes(os.urandom(2),'big')%16;B=Z.verify_v20_lat_b2_and_srt(seed=20261121);Al=Ao.pt_light(a(i));BV=Aa.verify_ht_pool(B);Bv=A5.alpha16();BW=AJ.stats();Bw=E(status=A1,how_to_run={s:Ax,'u':5,'j0':0,'n':40},note='regenerates stored miniature lattices j0..j0+n-1 of tower u (-3 or 5) from scratch and compares them with the stored records');A=M.compute_ledger(track=J,use_db=R,c_tail_z6=364.);EwDef=Bn.verify_exact_w2d_default(M);BX=DESC.verify_cells_theta(M.load_data()[G]);BY=DESC.verify_boot_gate_knapsack(A[G],M.load_data()[G]);x=M.load_data()['pool_quad'];Bx=x['bins'];By=x['nonprim_bins'];Bz=sum(A[2]for A in Bx)==x[Az]and sum(A[1]for A in By)==x[A_]and x['n_in']==2000000 and x[Az]+x[A_]==2000000 and BS[D];BZ={};Ba=R
	for(Bb,Bc)in Q[m].items():Bd=A[m][Bb];AD=Bd<=z(Bc)+1e-12;BZ[Bb]=E(computed=C(Bd,6),declared=Bc,ok=AD);Ba&=AD
	B_=h[D]and A[P]['ladder_all4_monotonic']and A[P]['all_runs_28_below_4caps']and A[P]['filtsim_cur28']['n']>=4 and h[A8][B0]>A[P][Ad][B0]and h[A8]['N_l_U']>A[P][Ad]['N_l_U']and h[A8][B1]>A[P][Ad][B1]and h[A8]['tab_r']<M.TAB_R_CAP;Am,C0=Z.fp_prices_from_kernels(a(i/AM));Be=E(counted_k0_to_k10=Am,charged_FP_CL=Z.FP_CL,root_k11_charged=Z.FP_CL[11],match=Am==Z.FP_CL[:11]and Z.FP_CL[11]>=y(Am[2:]));Bf=Z.wheel_prices_from_kernels(a(i/AM));Bg=E(counted=Bf,charged=Z.WHEEL_OPS,match=Bf==Z.WHEEL_OPS);Bh=DESC.verify_descbk_constants(M,a(i/AM),v);f=E(fp_node_prices_counted_from_kernels_ok=Be['match']and Al['prices'][D],wheel_prices_counted_from_kernels_ok=Bg['match'],algebra_sigma_w_and_crt_ok=AF[D]and AF['ntt3_shoup_butterfly_ok']and AF['ntt3_garner_crt_ok']and AF['ntt3_r4_r3_shoup_conv_ok'],exhaustive_fd_pure_int_bijection_ok=k[D]and k[Av]and k['sigma_inv_exact']and k['max_order_dedekind_ok']and k['exact_frac_enum_ok']and k['certified_interval_ok']and k['global_hexagon_proof_ok']and k['unit_index1_proved'],lattice_b2_score_nodes_and_srt_ok=B[D]and B['smat_cyclic_a345_ok']and B['zero_float_fp_enum']and B['strict_u64_i64_umulh_checked']and B['zero_exception_64bit_q38_proved']and B['cert_min_gsym_exact']==753729 and B['isqrt_u40_lut7_errs']==0 and B['isqrt_u40_lut7_max_over_pct']<=.12 and B['staddh_16b_semantic_errs']==0 and B['b2_x4_prefix_staddh_errs']==0 and B['lattice_b2_mismatches']==0 and B['complete_root_to_leaf_leaves']>=10000 and B['complete_root_to_leaf_mismatches']==0 and B['complete_root_to_leaf_apply_errs']==0 and B['complete_root_to_leaf_div_errs']==0 and B['complete_root_to_leaf_max_sa_cones']<=4 and B['wheel_tree_state_errs']==0 and B['wheel_tree_leaves_checked']>=10000 and B['swar21_cur_errs']==0 and B['complete_scaled_max_rel_err_pct']<=5. and B['complete_scaled_fp_ir_max_rel_err_pct']<=5. and B['prod_exec_leaves']>=2000 and B['prod_exec_mismatches']==0 and B['prod_exec_apply_errs']==0 and B['prod_exec_div_errs']==0 and B['swar_tie_canon_errs']==0 and B['barrett_kta_18op_errs']==0 and B['barrett_kta_18op_tested']>=3000 and B['desc_ai_prefix_chain_errs']==0 and B['score_node_replay_errs']==0 and B['fabric_1hop_addr_errs']==0 and B['b_tab_packing_errs']==0 and B['discrete_chord_bound_test_errs']==0 and B['fabric_1hop_list_node_errs']==0 and B['fabric_1hop_per_node_bytearray_zero_dict']and B['fabric_1hop_zero_frag_exact_ok']and B['greedy_node_bound_errs']==0 and B['greedy_bound_nodes']==16 and B['logp_bands_1_16b']==224 and B['lat_8b_direct_ptr_page_cap']==2046 and B['emit_full_pages_flushed']>=20 and B['zero_hash_table_apply']and B['inc_c0_exact_errs']==0 and B['prod_256orb_exact_rejection_ht']and B['prod_256orb_uniform_random_root_choice']and B['prod_256orb_count']==256 and B['prod_max_inv_b2_bits']<64 and B['prod_max_r2_q12_bits']<40 and B['prod_max_d2_bits']<64 and B['prod_max_ck_bits']<63 and B['swar_bounds_ok']and B['srt_batch_mismatches']==0 and B['e2e_factored_smooth_sqfree']>=10 and B['e2e_factored_1lp']>=100 and B['gate_0lp_fn']==0 and B['gate_0lp_fp']==0,live_pool_hnf_subsample_ok=BT[D],live_genuine_c6_rel_filter_ok=AW[D]and AW[Ae]==0 and AW[B2]==0,live_bareiss_12d_lll_limb_ok=v[D]and v['inductive_sr_inv_ok']and v['input_uniform_wc_proved']and v['pair_setup_p59_proved']and v['lat_b2_p61_proved'],live_mont5_cios_and_suyama12_ecm_ok=AG[D]and AG['ecm_1rung_t2_ok']and AG['ecm_1rung_w_le_b1']and AG['ecm_1rung_min_iv_ratio']>=o,live_z6_la_5limb_barrett_ok=A3[D]and A3['unit4_v4_expand_ok']and A3['unit4_8vec_ok']and A3['unit4_8vec_u64_idx_43b_ok']and A3['unit4_8vec_realistic_ir']<=M.C_UNIT_5,h3_scaling_and_nu_h_ok=q[D]and q['fixed_instance_deterministic_h3']and q['zero_loss_filtsim_ok']and q['all4_monotone_20_27']and q['hist_v4_exact_ok']and q['fg_dens_ok']and q['u_ram_all_lt_2p30'],trim_pairs_ok=Bu[D],split_crt_norms_and_mp_rt_ok=Ak[D]and Ak['chained_1pass_b2_errs']==0,pool_all_active_200k_ok=Bz,live_filtsim_z6_and_wmid_ok=h[D]and h['u12'][Ae]==0 and h['u12'][B2]==0 and h[A8][Ae]==0,filtsim_8scale_ladder_ok=B_,unplanted_5sm_bw_dlp_and_schur_ok=AE[D],live_62row_pure_k1_boole_wald_dag_ok=A[G][D]and A[G]['nodes']==62 and A[G]['mode']=='pure_k1_boole_wald_iid_dp_boot'and A[G]['p0_floor']==.001 and A[G][Aq]<=.5 and A[G]['e_runs_per_log']<=2. and A[G]['mu_min']>=5e2 and A[G]['p0_max']<=.005 and A[G]['p_node_abort']==.0 and A[G]['p_repeat_loss']==.0 and A[G]['p_descent_fail']==.0 and A[G]['desc_cap_log2']==112.,rows_ok=Ba,exp_run_ok=A[S]<=z(Q[S]),success_ok=A[B3]>=z(Q[B4]),gate_K_ok=A[AO]==Q['gate_K']and A['gate_passes']and A['gate_prev_fails'],claim_ok=A[c]<=z(Q[n])and N(A[Af])<=N(Q[Af]),ram_ok=A[F][B5]<=97.5798 and A[F][B5]<1e2 and A[F]['peak_global_ram_cap_log2']<=97.5798 and A[F][d][D]and A[F]['n_sc_nodes_exp']==34 and A[F]['cs_blocks_max']<1<<26 and A[B6]['storage_cap_abort']==.0 and A[F][B7]<=60.176 and A[F][B8]<=63.973 and A[F][B9]<=63.393 and A[F]['score_arena_log2']<64. and A[F][BA]<=57.58 and A[F][BB]<=52.51 and A[P][BC]<=63.42 and A[P]['fits_1own']and A[G][BD]<=60.47 and r(A[G][BD]-C(l.log2(M.pbits(0,A[G][BE],o)/8.),4))<.0001 and A[G][Ag]<=63.47 and r(A[G][Ag]-C(l.log2(8.*M.pbits(0,A[G][BE],o)/8.),4))<.0001 and A[F][d][Ah]['stageb_g_fits_ok']and A[F][d][Ah]['f_node_lt_2_64']and A[F][d]['stageb_g_L']==13 and A[F][d]['f_node_L']==50 and __import__(d).STAGEB_G_SP>=l.ceil(M.pbits(0,M.OP['p_min'],M.KH_G_CAP)/64.)and __import__(d).F_NODE_SP>=l.ceil(M.pbits(0,M.OP['lf'],M.KH_F_CAP)/64.)and A[F][d][BF]<=62.8 and A[F][d][Ah]['machine_lt_2_100']and y(A[F][B7],A[F][B8],A[F][B9],A[F][d][BF],A[F][BA],A[F][BB],A[P][BC],A[G][Ag])<64.,geo_consistency_ok=r(A['pool']['frac_reg_analytic']-A['pool']['frac_reg'])<.005 and A[F]['max_sb']<=8 and A[F]['sup_nf_bits']<448. and A[F]['sup_ng_bits']<=512. and A[F]['log2_Na_shell']<47.,manifest_subbounds_ok=A[F]['log2_V']<=94.8905 and 568.8<=A[F]['Ra']<=568.9 and A[F]['Rb_GA']<=755.6 and A[U]['rep_lo']<=1.905 and A[U]['occ_sb']<=.142 and A[U]['hi_ap']<=.596 and A[U][Ar]<=4.192 and A[U][As]<=.0218 and A[U]['scan_zero']<=.427 and A[U]['cof']<=.490 and A[U]['total']<=7.770,peel_ok=not A[P]['peel']['stalls']and A[P]['peel']['exc_crit']<.72);L=A[G];AX=E(A[m]);AX[AL]=L['boot_two_log_if_charge4'];L['claim_if_boot_charge_4']=C(l.log2(o+A[AO]/2.**64*2.**l.log2(sum(2.**A for A in AX.values()))),6);AX[AL]=L['boot_two_log_if_survivors_x4'];L['claim_if_boot_survivors_x4']=C(l.log2(o+A[AO]/2.**64*2.**l.log2(sum(2.**A for A in AX.values()))),6);f.update(pt_light_ok=Al[D],ht_pool_replay_and_pooled_bound_ok=BV[D],cells_exact_ok=BX[D],boot_gate_knapsack_ok=BY[D],boot_kept_rule_ledger_ok=L['boot_dp_lattices']==256 and L['boot_charge']==2. and L['boot_dp_surv_monotone']and L['boot_dp_pair_loss']<.001 and L['boot_dp_first_moment_factor']-3.*L['boot_dp_first_moment_factor_se']>=o and L['boot_dp_first_moment_exact_nl']-3.*L['boot_dp_first_moment_exact_nl_se']>=o and L['boot_dp_surv_exact_nl']+3.*L['boot_dp_surv_exact_nl_se']<=.5*L[BG]and L[BG]>=o,boot_mini_ratios_above_charge_ok=b(A['exhaustive']and A['tower_ok']and b(A['first_moment_ci95'][0]>=.5 and A['lattice_ci95'][0]>=.5 for A in A['by_L'].values())for A in BW.values()));f.update(descbk_constants_ok=Bh[D],exact_w2d_default_ok=EwDef[D]);Y=b(f.values());AY=E(pt_light=Al,ht_pool=BV,exact_w2d_default=EwDef,descent_cells_theta=BX,boot_gate_knapsack=BY,boot_dp_replay=E(status=A1,how_to_run={s:A7},note=AP),boot_dp_small_primes=Bv,boot_mini_replay=Bw,boot_mini_stats=BW,track=J,part=H,checks_ok=Y,all_ok=Y,checks=f,algebra=AF,exhaustive_fd=k,swar_sb_u7_stageB=E(status=A1,how_to_run={s:A7},note=AP),v20_lat_b2_and_srt=B,fp_kernel_prices=Be,wheel_kernel_prices=Bg,v25_upgrades=E(status=A1,how_to_run={s:A7},note=AP),v32_upgrades=E(status=A1,how_to_run={s:A7},note=AP),live_pool_hnf_subsample=BT,live_genuine_c6_rel_filter=AW,live_bareiss_12d_lll=v,live_mont5_cios_and_suyama12_ecm=AG,live_z6_la_kernels=A3,nonprim_hnf_witnesses=BS,diagnostics_part=E(status=A1,checks=['descent_12d_lattices_24sq_24boot_ok','deg8_dispersion_and_prod_yield_ok','c6_adaptive_h5_live_ok','c6_pgf_usable_depth2_and_hea_campaign_ok','live_trim_pair_ok'],how_to_run={s:Ay},note='H5 scoped diagnostics (ADAPTIVE.md section 2) are not run in the default run (time limit); run them with part diagnostics'),h3_scaling_and_nu_h=q,split_crt_norms_and_mp_rt=Ak,live_filtsim_z6=h,hx_unplanted_5sm_bw=AE,ledger=E(rows=BZ,exp_run=C(A[S],9),exp_run_decl=Q[S],success=C(A[B3],8),success_decl=Q[B4],gate_K=A[AO],claim=C(A[c],9),claim_decl=Q[n],score_numerator=Q[Af],fail_terms=A[B6]),geometry=A[F],per_cell=A[U],la=A[P],descent=A[G]);AY[At]=Bh
	if H in(BH,'full'):Bi=M.compute_ledger(track=J,use_db=False);Bj=M.compute_ledger(track=J,use_db=R,strict_per_owner=R,swar_probe=1e1);Bk=M.compute_ledger(track=J,use_db=R,strict_per_owner=False);Bl=M.compute_ledger(track=J,use_db=R,cap_scale=.8);Bm=M.compute_ledger(track=J,use_db=R,srt_g_cap=.3);AY[BH]=E(plain_dickman_rho=E(exp_run=C(Bi[S],6),claim=C(Bi[c],6)),conservative_lat_b2_plus10=E(exp_run=C(Bj[S],6),claim=C(Bj[c],6)),conservative_srt_030_cap=E(exp_run=C(Bm[S],6),claim=C(Bm[c],6)),flat_shared_arena_ref=E(exp_run=C(Bk[S],6),claim=C(Bk[c],6)),de_bruijn_tighter_cap=E(exp_run=C(Bl[S],6),claim=C(Bl[c],6)))
	AY['elapsed_s']=C(time.time()-W,3);T(json.dumps(AY))
	if not Y:raise A6(1)
if __name__=='__main__':main()
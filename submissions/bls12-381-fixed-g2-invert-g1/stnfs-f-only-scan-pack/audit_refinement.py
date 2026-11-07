"""Stdlib ledger/kernel audit from recorded campaigns; optional descent replay."""
import base64,hashlib,json,lzma,math,re,sys
from fractions import Fraction as F
from decimal import Decimal,localcontext
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE/"files" if(HERE/"files").is_dir()else HERE
sys.path.insert(0,str(ROOT))
import descent16 as D,model16 as M,retain16 as R,owner10 as O
CLAIM="117.562576929701"
ORDER=Decimal("52435875175126190479447740508185965837690552500527637822603658699938581184513")
def derive(c,ck):
	# Recorded campaign totals are inputs; this recomputes their complete ledger.
	import lattice_count_bounds as B
	cfg=json.loads((ROOT/"config.json").read_text())
	pd=json.loads((ROOT/"pool_d.json").read_text(),parse_float=F)
	th=F(cfg["production_theta_region"]);tp=F(cfg["production_theta"])
	N=sum(pd["draws"][k]for k in pd["work"])
	q={k:sum(pd["sums"][s][k]for s in pd["work"])for k in M.SUM_KEYS}
	n=q["rows"];w=q["work"]
	se2=(n*q["work2"]-w*w)/((n-1)*w*w)+(1-F(n,N))/n
	allow=F(str(M.CANON_WORK_ALLOWANCE))
	ck["work_allowance_recomputed"]=(allow-1)**2>=F(str(M.WORK_SE_K))**2*se2 and se2<=F("0.004214")**2
	radius=math.isqrt(int(th*cfg["production_A_squared"]))
	corr=B.correction(cfg,radius);pi=B.pi_interval()[1]
	sd=B.sqrt_interval(B.determinant_rational(cfg["matrix_a"]))[0]
	va=pi**3/6*(th*cfg["production_A_squared"])**3/(2*sd)
	col=va/(N*(th/tp)**3)*w*allow/F(str(M.UNIT))
	ck["ledger_model_rows"]=[float(v)for v in c["ledger_upper_rows"]][3:]==[M.DESCENT_CAP,M.LOG2_CANON,M.LOG2_SM,M.LOG2_MISC,M.LOG2_FBC,M.LOG2_WIDE_LAT,M.LOG2_PREP7,50,M.LOG2_VERIFY,M.LOG2_GTREE,M.LOG2_GSPLIT,M.LOG2_SPLIT]and M.UNIT==cfg["reference_cost_instructions"]==20000
	if M.C_ROW==2710.0:O.install(M)
	with localcontext()as x:
		x.prec=80;d=lambda a:Decimal(a.numerator)/Decimal(a.denominator)if isinstance(a,F)else Decimal(str(a))
		ln2=Decimal(2).ln();lg=lambda a:a.ln()/ln2;p2=lambda a:(d(a)*ln2).exp()
		col_log=lg(d(col));z=Decimal(69)*ln2
		U=p2(69)/z*(1+1/z+2/z**2)+p2(31);nn=U*d(F(str(M.NU300)))
		m,nb=d(M.BW_M),d(M.BW_N)
		prod=d(F(str(M.W_ROW))*F(str(M.C_ENTRY_K))+F(str(M.C_ROW))+M.N_SM*F(str(M.C_FULL)))
		per=prod*(nb/m+1)+(1/m+1/nb)*(prod+nb*d(M.C_FULL)+d(M.C_RED9))+d(M.C_PROJ)*(nb/m+1)
		la=(nn**2*per+nn*(m+nb)*d(M.C_RED9))/d(M.UNIT)
		L=(nn+d(M.N_SM))*(1/m+1/nb)+2;k=m+nb
		gen=4*k**3*d(M.C_FULL+M.C_RED9)/d(M.UNIT)*(lg(L)*Decimal(3).ln()).exp()*lg(L)
		ck["collection_recomputed"]=col_log<=d(c["ledger_upper_rows"][0])
		ck["la_generator_recomputed"]=lg(la)<=d(c["ledger_upper_rows"][1])and lg(gen)<=d(c["ledger_upper_rows"][2])and nn<6*p2(56)
		ck["wrapper_recomputed"]=p2(80)/d(M.UNIT)<p2(c["wrapper_log2"])and c["wrapper_log2"]==66
		derived=dict(collection_log2_upper=str(col_log),la_log2=str(lg(la)),generator_log2=str(lg(gen)),work_se=str(d(se2).sqrt()))
	manifest=json.loads((ROOT/"submission.json").read_text())["claims"][0]["score"]
	ck["manifest_claim"]=F(int(manifest["numerator"]),int(manifest["denominator"]))==F(CLAIM)
	text=lzma.decompress(base64.b85decode("".join((ROOT/"kernels_asm.b85").read_text().split()))).decode()
	parts=re.split(r"(?m)^@@ ([^\n]+)\n",text);sections=dict(zip(parts[1::2],parts[2::2]));oh=sections["measurements-oh.txt"]
	scan=[(int(a),F(int(b)-int(c),int(e)))for a,b,c,e in re.findall(r"(?m)^OH scan (\d+) (\d+) (\d+) (\d+)$",oh)]
	cell,leaf=F(str(M.C_SCAN_CELL)),F(str(M.C_SCAN_LEAF));need=max(L*(v-cell)for L,v in scan)
	v=dict(scan);slope=(256*v[256]-64*v[64])/192
	ck["f_only_scan_exact"]=len(scan)==65 and {L for L,v in scan}==set(range(1,65))|{256} and slope<=cell and need<=leaf and (cell,leaf)==(F("1.376"),F("37.62"))and "19800 segments"in oh and "10119894 survivors, 0 mismatches"in oh
	cr={int(a):(int(b),int(c)-int(e)-int(f)-int(g))for a,b,c,e,f,g in re.findall(r"(?m)^CELLRES (\d+) (\d+) (\d+) (\d+) (\d+) (\d+)$",oh)}
	res=F(cr[26][1]-cr[24][1],cr[26][0]-cr[24][0]);ck["cell_residual_exact"]=res<=F(str(M.CELL_ROUTE))/4
	bm=sections["measurements12.txt"].split("== nbound31f (f bound only) measurements")[1]
	a,b,c0,e,n=map(int,re.search(r"bench 1/3 reps: (\d+) (\d+) ; benchsetup 1/3: (\d+) (\d+) ; evaluations/rep (\d+)",bm).groups())
	price=F((b-a)-(e-c0),2*n)
	ck["f_only_bound_exact"]=price<=F(str(M.C_NB_EVAL))==F("1154.6")and "f bound mismatches 0; f product bit mismatches 0"in sections["measurements12.txt"]
	derived.update(scan_cell_slope_exact=str(slope),scan_leaf_excess_exact=str(need),bound_price_exact=str(price),cell_residual_exact=str(res))
	return derived
def main():
	raw=sys.stdin.read();req=json.loads(raw)if raw.strip()else{}
	if req not in({},{"descent":True},{"descent":False}):
		raise ValueError('input: {} or {"descent": true}')
	c=json.loads((ROOT/"refinement.json").read_text());ck={}
	derived=derive(c,ck)
	r0=M.retentions()[0];sp=M.band_split(*c["split_inputs"],M.rmin(r0,M.OP['slack_f']),M.rmin(r0,M.OP['slack_g']))
	sb,sf=M.srt_toy(6,1);import base64,lzma,re;kt=re.split(r"(?m)^@@ fft.txt\n",lzma.decompress(base64.b85decode("".join((ROOT/"kernels_asm.b85").read_text().split()))).decode())[1];fp=[x.split(":")for l in kt.splitlines()if l.startswith("FFTPTS ")for x in l.split()[1:]];ck["srt"]=sb==0 and sf<=2 and len(fp)>=25 and all(float(c)<=M.fft_F(int(n))for n,c in fp)and M.gbatch_cost()<=M.C_GBATCH and M.gbatch_cost(pbits=1.000028*2**48/M.LN2,leaf=13,L=39,c_leaf=2**19)<=M.C_BATCH and (M.ecm_instr(69.0,576)*12*1.002+13*M.c_mr(576,M.MONT5))/M.UNIT*2**65.1456<2**M.LOG2_GSPLIT
	ck["band_split"]=sp["log2_units"]<=M.LOG2_SPLIT==float(c["ledger_upper_rows"][-1])
	ck["owner_schedule"]=O.ACTIVE_OWNERS==6 and O.row_tail()==1210 and O.C_ROW==1215.0 and O.LA_BOUND==float(c["ledger_upper_rows"][1])
	identity_bad=[n for n,h in c["code_sha256"].items()if hashlib.sha256((ROOT/n).read_bytes()).hexdigest()!=h]
	ck["code_identity"]=not identity_bad
	ck["parameters"]=D.P.MESH==384 and D.RSTEP==1/1024 and D.P.TRUNC==0 and M.DESCENT_CAP==c["descent_cap"]and D.P.LOG2_RELS==c["log2_rels_input"]<=c["log2_rels_2se"]and D.P.BOOT_RELS==18 and D.P.UDUP==0.65
	cf=json.loads((ROOT/"config.json").read_text());pa=json.loads((ROOT/"pool_d.json").read_text())["anchors"]
	ck["warm_start_anchors"]=all(list(M.pool_row_warm(cf,x[0],x[1],x[11])if M.RUN_T[0]<=x[12]-x[11]+1<=M.RUN_T[1]else(x[6],x[8]))==x[14:16]for x in pa)
	ck["windows"]=all(D.PLAN16[(b,s)]==(o[i][0],8*o[i][1],*o[i][2:])for b,o in D._PLAN.items()
	for i,s in enumerate("fg"))and D.BOOT16==(215,328,16,4)
	bad=D.P();bad.MESH=48
	try:
		D.evaluate(bad,D.PLAN16,D.BOOT16);ck["parameter_mismatch_rejected"]=False
	except ValueError as e:
		ck["parameter_mismatch_rejected"]="must agree" in str(e)
	if not req.get("descent"):
		ck["lazy_bound"]=all(D.lazy_mc(mu,r,psi,trials=2000)<=D.lazy_fraction(mu,r,psi)+0.01 for mu,r,psi in c["lazy_probes"])
	probes=[("kA0",(600,69.0,75.0,24,0,0,48,1.0,0.0)),("kA1",(600,69.0,75.0,42,1,0,48,1.0,0.0)),
	("kA4",(600,69.0,100.0,24,4,1,48,1.0,0.0)),("kA15",(1200,69.0,215.0,43,15,0,48,1.0,0.0))]
	ctl={x["name"]:x for x in c["probe_checks"]if x["rho_den"]==1024}
	for name,a in probes:
		lo=D.ss_lower2(*a);hi=R.ss_upper2(*a[:-1])
		ck[name]=0<lo<=hi and math.isclose(lo,ctl[name]["limited"],rel_tol=1e-12)and math.isclose(hi,ctl[name]["upper_limited"],rel_tol=1e-12)
	out={}
	if req.get("descent"):
		import counts16,base64,lzma
		R.load_table(json.loads(lzma.decompress(base64.b85decode((ROOT/"stats16.b85").read_text()))))
		res=D.evaluate(D.P(),D.PLAN16,D.BOOT16);EN,_=counts16.expected_nodes(D.P(),res);bt=res["boot"]
		rows=list(res["DET"].values())
		out=dict(two_log=res["two_log_units"],boot_phi=bt["phi"],boot_psi=bt["psi"],log2_EN=math.log2(EN),
		mu_min=min(d["mu"]for d in rows),phi_max=max(res["PHI"].values()),psi_max=max(max(d["psi"]for d in rows),bt["psi"]))
		ck["descent_recomputed"]=out["two_log"]<=c["two_log_upper"]and bt["phi"]<=c["boot_phi_upper"]and math.log2(EN)-M.NODE_CAP_LOG2<=c["counter_abort_log2_upper"]and out["mu_min"]>=D.P.MIN_RELS and out["phi_max"]<=0.05 and out["psi_max"]<=D.P.PSI_MAX and abs(out["two_log"]-c["two_log_measured"])<1e-9
	with localcontext()as x:
		x.prec=80;two=Decimal(2)
		p2=lambda v:(Decimal(str(v))*two.ln()).exp()
		core=sum(p2(v)for v in c["ledger_upper_rows"]);core_log2=core.ln()/two.ln()
		markov=p2(Decimal(str(c["two_log_upper"]))-Decimal(str(c["descent_cap"])))
		succ=(1-two**-66-Decimal(str(c["cap_abort"])))*(1-2*Decimal(str(c["boot_phi_upper"]))-markov-p2(c["counter_abort_log2_upper"])
		-p2(c["repeat_log2_upper"])-2/ORDER)
		tot=p2(c["ledger_upper_claim"])+p2(c["wrapper_log2"]);tot_log2=tot.ln()/two.ln()
		s_w=Decimal(str(c["success_lower_claim"]))-two**-48
		ck["mean_rounding"]=c["two_log_measured"]<c["two_log_upper"]and c["boot_phi_measured"]<=c["boot_phi_upper"]and c["node_log2_measured"]-M.NODE_CAP_LOG2<=c["counter_abort_log2_upper"]
		ck["ledger"]=core_log2<=Decimal(str(c["ledger_upper_claim"]))<122
		ck["ledger_rows_stated"]=core_log2<=Decimal(c["ledger_rows_sum_upper"])<=Decimal(str(c["ledger_upper_claim"]))
		ck["success"]=succ>=Decimal(str(c["success_lower_claim"]))>Decimal("0.39")
		ck["success_headroom"]=(1-two**-66-Decimal(str(c["cap_abort"])))*(1-2*Decimal(str(c["boot_phi_upper"]))-Decimal("1.1")*markov-p2(c["counter_abort_log2_upper"])-p2(c["repeat_log2_upper"])-2/ORDER)>=Decimal("0.39")+two**-48
		ck["claim_budget"]=tot_log2<=Decimal(CLAIM)and Decimal(CLAIM)==Decimal(c["claim"])
		ck["wrapped_success"]=s_w>Decimal("0.39")
	print(json.dumps({"checks_ok":all(ck.values()),"checks":ck,"identity_mismatches":identity_bad,"core_ledger_upper_log2":str(core_log2),
	"markov_upper":str(markov),"core_success_lower":str(succ),"log2_upper_with_wrapper":str(tot_log2),"declared":CLAIM,"success_lower_with_wrapper":str(s_w),"descent":out or "recorded (rerun with {\"descent\": true})","band_split":sp,"derived":derived,
	"scope":"Ledger and exact kernel accounting recomputed from inherited recorded pool sums and measurement logs; campaigns are not regenerated. evidence16.py runs the broader checks. Conditional on H1-H5, H-EA, H-X, H-FFT."}))
	if not all(ck.values()):
		raise SystemExit(1)
if __name__=="__main__":
	main()

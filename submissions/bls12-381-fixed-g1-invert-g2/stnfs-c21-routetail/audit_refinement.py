"""Hosted audit: {} recorded certificate; {"descent": "law_f"|"law_g"} recompute a stored omitted-mass law;
{"descent": "tree"|true} the descent on the stored laws; {"encountered": seed} [7e]."""
import hashlib,json,math,sys
from decimal import Decimal,localcontext
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE/"files" if(HERE/"files").is_dir()else HERE
sys.path.insert(0,str(ROOT))
import descent16 as D,model16 as M,retain16 as R,owner10 as O,common12 as C
CLAIM="115.658012484016"
ORDER=Decimal("52435875175126190479447740508185965837690552500527637822603658699938581184513")
def main():
	req=json.load(sys.stdin)
	if not(req=={}or list(req)==["descent"]and req["descent"]in(True,False,"tree","law_f","law_g")or list(req)==["encountered"]and type(req["encountered"])is int):raise ValueError("input: see docstring")
	c=json.loads((ROOT/"refinement.json").read_text());ck={}
	sp="with {}"
	if not req:
		r0=M.retentions()[0];sp=M.band_split(*c["split_inputs"],M.rmin(r0,M.OP['slack_f']),M.rmin(r0,M.OP['slack_g']))
		sb,sf=M.srt_toy(6,1);import base64,lzma,re;kt=re.split(r"(?m)^@@ fft.txt\n",lzma.decompress(base64.b85decode("".join((ROOT/"kernels_asm.b85").read_text().split()))).decode())[1];fp=[x.split(":")for l in kt.splitlines()if l.startswith("FFTPTS ")for x in l.split()[1:]];ck["srt"]=sb==0 and sf<=2 and len(fp)>=25 and all(float(c)<=M.fft_F(int(n))for n,c in fp)and M.gbatch_cost()<=M.C_GBATCH and M.gbatch_cost(pbits=1.000028*2**48/M.LN2,leaf=13,L=39,c_leaf=2**19)<=M.C_BATCH and (M.ecm_instr(69.0,576)*12*1.002+13*M.c_mr(576,M.MONT5)+M.ecm_instr(69.0,840)*17*1.002+17*2**26+18*M.c_mr(840,M.MONT5))/M.UNIT*2**65.1456<2**M.LOG2_GSPLIT and M.gbatch_cost(pbits=1.000028*2**64/M.LN2,leaf=13,L=54,c_leaf=2**19)<=M.C_BATCH2
		ck["band_split"]=sp["log2_units"]<=M.LOG2_SPLIT==float(c["ledger_upper_rows"][-1])
	ck["owner_schedule"]=O.ACTIVE_OWNERS==4 and O.row_tail()==857 and O.C_ROW==857.0 and O.LA_BOUND==float(c["ledger_upper_rows"][1])
	ck["code_identity"]=all(hashlib.sha256((ROOT/n).read_bytes()).hexdigest()==h for n,h in c["code_sha256"].items())
	ck["parameters"]=D.P.MESH==384 and D.RSTEP==1/1024 and D.P.TRUNC==0 and M.DESCENT_CAP==c["descent_cap"]and D.P.LOG2_RELS==c["log2_rels_input"]<=c["log2_rels_2se"]and D.P.BOOT_RELS==18 and D.P.UDUP==0.65
	st=C.stats();cf,pa=st["config"],st["pool_d"]["anchors"]
	if not req:ck["warm_start_anchors"]=all(list(M.pool_row_warm(cf,x[0],x[1],x[11])if M.RUN_T[0]<=x[12]-x[11]+1<=M.RUN_T[1]else(x[6],x[8]))==x[14:16]for x in pa)
	ck["windows"]=all(D.PLAN16[(b,s)]==(o[i][0],8*o[i][1],*o[i][2:])for b,o in D._PLAN.items()
	for i,s in enumerate("fg"))and D.BOOT16==(215,304,15,4)
	bad=D.P();bad.MESH=48
	try:
		D.evaluate(bad,D.PLAN16,D.BOOT16);ck["parameter_mismatch_rejected"]=False
	except ValueError as e:
		ck["parameter_mismatch_rejected"]="must agree" in str(e)
	probes=[("kA0",(600,69.0,75.0,24,0,0,48,1.0,0.0)),("kA1",(600,69.0,75.0,42,1,0,48,1.0,0.0)),
	("kA4",(600,69.0,100.0,24,4,1,48,1.0,0.0)),("kA15",(1200,69.0,215.0,43,15,0,48,1.0,0.0))]
	ctl={x["name"]:x for x in c["probe_checks"]if x["rho_den"]==1024}
	for name,a in probes:
		lo=D.ss_lower2(*a);hi=R.ss_upper2(*a[:-1])
		ck[name]=0<lo<=hi and math.isclose(lo,ctl[name]["limited"],rel_tol=1e-12)and math.isclose(hi,ctl[name]["upper_limited"],rel_tol=1e-12)
	out={}
	dv=req.get("descent");R.load_table(C.stats())
	if dv in("law_f","law_g"):
		sd=dv[-1];pdf,_=R.omitted_pdf(sd,True);mean=sum(x*v for x,v in pdf.items())/R.GRID;tail=R.cell_loss(sd,96)
		st=c["laws"][sd];out=dict(side=sd,mean_bits=mean,p_above_96=tail,stated=st)
		ck["law_"+sd]=R.law(sd)==R.LAWS[sd]and round(mean,2)==st[0]and tail<=st[1]
	if req.get("encountered")is not None:
		import descent_witness as W
		mm,ok=W.encountered_check(req["encountered"]);ck["encountered"]=ok
		out=dict(mean_cf_enc_ctl=mm[0][:2],mean_cg_enc_ctl=mm[1][:2],max_cf_cg_enc=[mm[0][2],mm[1][2]],rule="enc <= ctl + 10 (cf), + 4 (cg); max <= 97/14.5")
	if dv in(True,"tree"):
		import counts16
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
	print(json.dumps({"checks_ok":all(ck.values()),"checks":ck,"core_ledger_upper_log2":str(core_log2),
	"markov_upper":str(markov),"core_success_lower":str(succ),"log2_upper_with_wrapper":str(tot_log2),"declared":CLAIM,"success_lower_with_wrapper":str(s_w),"descent":out or "recorded (rerun with {\"descent\": true})","band_split":sp,
	"scope":"Recorded full-run certificate plus selected recomputation; evidence16.py reruns everything. Conditional on H1-H5, H-EA, H-X, H-FFT."}))
	if not all(ck.values()):
		raise SystemExit(1)
if __name__=="__main__":
	main()

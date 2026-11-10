import hashlib,json,math,sys
from decimal import Decimal,localcontext
from fractions import Fraction
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE/"files" if(HERE/"files").is_dir()else HERE
sys.path.insert(0,str(ROOT))
import descent16 as D,model16 as M,retain16 as R,owner10 as O,common12 as C
CLAIM="115.235972269343"
ORDER=Decimal("52435875175126190479447740508185965837690552500527637822603658699938581184513")
def main():
	req=json.load(sys.stdin)
	if not(set(req)<={"descent","ideal"}and req.get("descent","tree")in(True,"tree")and type(req.get("ideal",0))is int and len(req)<2):raise ValueError("one input key")
	c=json.loads((ROOT/"refinement.json").read_text());ck={}
	sp="with {}"
	if not req:
		r0=M.retentions()[0];sp=M.band_split(*c["split_inputs"],M.rmin(r0,M.OP['slack_f']),M.rmin(r0,M.OP['slack_g']))
		import base64,lzma,re;kt=re.split(r"(?m)^@@ fft.txt\n",lzma.decompress(base64.b85decode("".join((ROOT/"kernels_asm.b85").read_text().split()))).decode())[1];fp=[x.split(":")for l in kt.splitlines()if l.startswith("FFTPTS ")for x in l.split()[1:]];ck["fft"]=len(fp)>=25 and all(float(c)<=M.fft_F(int(n))for n,c in fp)and M.gbatch_cost()<=M.C_GBATCH_R and M.gbatch_cost(pbits=1.000028*2**48/M.LN2,leaf=13,L=39,c_leaf=2**19)<=M.C_BATCH_R and M.ssa_toy(256,5)==(True,True)and (M.ecm_instr(69.0,576)*12*3+13*M.c_mr(576,M.MONT5)+M.ecm_instr(69.0,840)*17*3+17*2**26+18*M.c_mr(840,M.MONT5))/M.UNIT*2**65.1456<2**M.LOG2_GSPLIT and M.gbatch_cost(pbits=1.000028*2**64/M.LN2,leaf=13,L=54,c_leaf=2**19)<=M.C_BATCH2_R
		ck["band_split"]=sp["log2_units"]<=M.LOG2_SPLIT==float(c["ledger_upper_rows"][-1])
	ck["owner_schedule"]=O.ACTIVE_OWNERS==3 and O.row_tail()==695 and O.C_ROW==700.0 and O.LA_BOUND==float(c["ledger_upper_rows"][1])
	ck["code_identity"]=hashlib.sha256(b"".join((ROOT/n).read_bytes()for n in c["code_files"])).hexdigest()==c["code_sha256"]
	ck["parameters"]=D.P.MESH==384 and D.RSTEP==1/1024 and D.P.TRUNC==0 and M.DESCENT_CAP==c["descent_cap"]and D.P.LOG2_RELS==c["log2_rels_input"]<=c["log2_rels_2se"]and D.P.BOOT_RELS==100 and D.P.UDUP==0.65 and D.P.CD==1.2 and D.P.YF==0.8
	st=C.stats();cf,pa=st["config"],st["pool_d"]["anchors"]
	if not req:ck["warm_start_anchors"]=all(list(M.pool_row_warm(cf,x[0],x[1],x[11])if M.RUN_T[0]<=x[12]-x[11]+1<=M.RUN_T[1]else(x[6],x[8]))==x[14:16]for x in pa)
	ck["windows"]=all(D.PLAN16[(b,s)]==(o[i][0],8*o[i][1],*o[i][2:])for b,o in D._PLAN.items()
	for i,s in enumerate("fg"))and D.BOOT16==(215,304,15,4)
	probes=[("kA0",(600,69.0,75.0,24,0,0,48,1.0,0.0)),("kA1",(600,69.0,75.0,42,1,0,48,1.0,0.0)),
	("kA4",(600,69.0,100.0,24,4,1,48,1.0,0.0)),("kA15",(1200,69.0,215.0,43,15,0,48,1.0,0.0))]
	for name,a in probes:
		lo=D.ss_lower2(*a);hi=R.ss_upper2(*a[:-1])
		ck[name]=0<lo<=hi
	out={}
	if"ideal"in req:
		import descent_witness as W;i=req["ideal"];r=D._D["disc"][i-1]if i else None
		x=W.camp_ideal(i,r[4]if i else 0);out.update(ideal=x,row=r)
		ck["ideal_replay"]=(x["yr_lo"]>=2*D.P.UDUP and x["raw_lo"]>=2)if i else x["yr"]>=D.P.UDUP
	dv=req.get("descent");R.load_table(C.stats())
	if dv in(True,"tree"):
		Pm=D.P();res=D.evaluate(Pm,D.PLAN16,D.BOOT16)
	if dv in(True,"tree"):
		import counts16
		EN,_=counts16.expected_nodes(D.P(),res);bt=res["boot"]
		rows=list(res["DET"].values())
		out.update(two_log=res["two_log_units"],boot_phi=bt["phi"],boot_psi=bt["psi"],log2_EN=math.log2(EN),
		mu_min=min(d["mu"]for d in rows),phi_max=max(res["PHI"].values()),psi_max=max(max(d["psi"]for d in rows),bt["psi"]))
		ck["descent_recomputed"]=out["two_log"]<=c["two_log_upper"]and bt["phi"]<=c["boot_phi_upper"]and math.log2(EN)-M.NODE_CAP_LOG2<=c["counter_abort_log2_upper"]and out["mu_min"]>=D.P.MIN_RELS and out["phi_max"]<=0.07 and out["psi_max"]<=D.P.PSI_MAX and abs(out["two_log"]-c["two_log_measured"])<1e-9
	with localcontext()as x:
		x.prec=80;two=Decimal(2)
		p2=lambda v:(Decimal(str(v))*two.ln()).exp()
		core=sum(p2(v)for v in c["ledger_upper_rows"]);core_log2=core.ln()/two.ln()
		markov=p2(Decimal(str(c["two_log_upper"]))-1-Decimal(str(M.ATT_CAP)));pm=Decimal(str(c["boot_phi_upper"]))+markov+p2(c["counter_abort_log2_upper"])+p2(c["repeat_log2_upper"])
		ck["attempt_model"]=pm*Decimal(str(M.P_TOL))<=Decimal(str(M.P_HI))and M.ATT_CAP>=c["two_log_upper"]-1+math.log2(M.CAP_MULT)and M.DESCENT_CAP==c["descent_cap"]
		succ=(1-two**-66)*(1-2*Decimal(str(M.P_HI))**M.ATT_K-2/ORDER)
		lo=succ-two**-48;K=M.gate_K(Fraction(c["lo_run"]))
		core_meter=p2(c["ledger_upper_claim"])+p2(c["wrapper_log2"]);ecore=core-p2(c["ledger_upper_rows"][3])+p2(M.ATT_CAP)*2/Decimal(str(1-M.P_HI))
		exp_w=Decimal(str(M.GATE_UNITS))+Decimal(K)/two**64*(p2(c["expected_ledger_claim"])+p2(c["wrapper_log2"]));exp_log2=exp_w.ln()/two.ln()
		ck["expected_ledger"]=ecore.ln()/two.ln()<=Decimal(c["expected_ledger_claim"])<=Decimal(c["ledger_upper_claim"])
		ck["mean_rounding"]=c["two_log_measured"]<c["two_log_upper"]and c["boot_phi_measured"]<=c["boot_phi_upper"]and c["node_log2_measured"]-M.NODE_CAP_LOG2<=c["counter_abort_log2_upper"]
		ck["ledger"]=core_log2<=Decimal(str(c["ledger_upper_claim"]))<122
		ck["opened_run_success"]=lo>=Decimal(c["lo_run"])and M.P_CAP==0
		ck["gate"]=K==c["gate_K"]and Fraction(c["lo_run"])*K>=Fraction(39,100)*2**64>Fraction(c["lo_run"])*(K-1)
		ck["claim_budget"]=exp_log2<=Decimal(CLAIM)and Decimal(CLAIM)==Decimal(c["claim"])
	print(json.dumps({"checks_ok":all(ck.values()),"checks":ck,"p_model":str(pm),"run_success":str(lo),"gate_K":K,"every_run_log2":str(core_meter.ln()/two.ln()),"expected_log2":str(exp_log2),"declared":CLAIM,"descent":out,"band_split":sp,"toy":None if req else __import__("toy_descent").run({})}))
	if not all(ck.values()):
		raise SystemExit(1)
if __name__=="__main__":
	main()

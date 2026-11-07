"""Hosted audit (stdlib, one JSON object): {} checks the recorded certificate; {"descent": true} also reruns the descent."""
import hashlib,json,math,sys
from decimal import Decimal,localcontext
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE/"files" if(HERE/"files").is_dir()else HERE
sys.path.insert(0,str(ROOT))
import descent16 as D,model16 as M,retain16 as R,owner10 as O
CLAIM="117.525153583201"          # core meter, dense-Fp12 units (C=20,000)
PLAIN_CLAIM="120.847081678088"    # plain bound, G1 units (C=2,000) = core + log2(10)
ORDER=Decimal("52435875175126190479447740508185965837690552500527637822603658699938581184513")
"""Paired record-work certificate. --replay regenerates every pool draw (expensive).
Default checks recorded sufficient statistics, bounds and eight regenerated anchors.
"""
import json,math,sys
from decimal import Decimal as DX,localcontext
from fractions import Fraction as F
from pathlib import Path
import model16 as M,common12 as C,lattice_count_bounds as L,runs16 as RE

def context():
 c=json.loads((ROOT/'config.json').read_text());t=F(c['production_theta_region'])
 b=L.correction(c,math.isqrt(int(t*c['production_A_squared'])))
 return c,t,M.Geometry(c,t,F(b['safe_yield_factor']),F(b['safe_work_population_factor']))

def delta(qa,T,lrs,g,dd):
 S=g.slice_cells(qa);lg=math.log2(g.vb*g.upper)+lrs;v=0.;used=0
 if M.RUN_T[0]<=T<=M.RUN_T[1]:
  for j in range(16):
   lo=52+j/4;hi=lo+.25;mid=(lo+hi)/2
   if lg-lo<=M.REC7_EXT_MAX and lg-hi>=M.REC7_TAB[-1][0]:
    dp=(M.NBAR/2)*(2**M.log2_li(hi)-2**M.log2_li(lo))
    v+=M.NBAR*math.log(hi/lo)*(M.step7(lg-hi)+M.EXTRA7-M.tab(M.REC_TAB,math.log2(S)-mid))+dp*M.C_SETUP6*(2**-lg-1/S)
    used+=1
 return (v+M.RUN_ALLOW*dd/S)*M.SIEVE_SIDES,v*M.SIEVE_SIDES,used

def extension_audit():
 c,t,g=context();e=json.loads((ROOT/'pool_d.json').read_text())['extension'];pd=json.loads((ROOT/'pool_d.json').read_text())
 d52=M.contam_sup(g);p=M.P7
 try:M.P7=56.;d56=M.contam_sup(g)
 finally:M.P7=p
 dd=d52-d56;B=e['B'];ck={}
 ck['schedule']=(M.P7,M.P7_FULL,M.X7,M.REC7_EXT_MAX)==(52.,56.,64.,21.99)
 ck['contamination']=B>=M.SIEVE_SIDES*M.RUN_ALLOW*dd*((1-float(t))/(g.lam-float(t)))**3 and math.isclose(dd,e['Ddiff'],rel_tol=1e-13)
 with localcontext()as x:
  x.prec=80
  for k in ('S','G'):
   old=e['sums'][k]['old'];new=pd['sums'][k];d=e['delta'][k]
   mean=DX(str(old['work']))+DX(str(d['work_delta']))
   second=DX(str(old['work2']))+2*DX(B)*DX(str(old['work']))+DX(old['rows'])*DX(B)**2
   ck[k+'_moments']=DX(str(new['work']))>=mean and DX(str(new['work2']))>=second and math.isclose(new['work'],float(mean),rel_tol=1e-11) and math.isclose(new['work2'],float(second),rel_tol=1e-11)
   ck[k+'_fixed']=all(new[n]==old[n]for n in old if n not in ('work','work2','ysum','y2','nfy','ngy')) and d['max_extension_delta_per_cell']==0
 ck['boundary']=sum(pd['sums'][k]['rows']-e['delta'][k]['rows']for k in ('S','G'))==4
 sub={k:{n:0. for n in ('ysum','y2','nfy','ngy')}for k in ('S','G')};ret=M.retentions()[0];outside=True
 for pool,seed,index in e['boundary']:
  a,b,aq,bq=C.regenerate_pair(c,seed,index);q,dq,_=C.integer_metric(c['matrix_a']);qa=round(F(aq,dq)/c['production_A_squared']*65536)/65536
  outside&=F(aq,dq)>t*c['production_A_squared']and qa<=float(t)
  nf,ng=[round(16*math.log2(z))/16 for z in C.norms(a,b)];uf=(nf+M.OP['alpha_f'])/M.OP['lb'];ug=(ng+M.OP['alpha_g'])/M.OP['lb']
  y=(1-qa)**3*C.rho(uf)*C.rho(ug)*ret.at(nf,M.OP['slack_f'])['retention']*M.ECM_FIND**M.lam_top(uf,M.OP['lx2']/M.OP['lb'])*(1-M.BAND_FMAX[0]*(1-M.BAND_Q)**M.BAND_ROUNDS)
  for n,v in dict(ysum=y,y2=y*y,nfy=nf*y,ngy=ng*y).items():sub[pool][n]+=v
 ck['boundary_yield']=outside and len(e['boundary'])==4 and all(math.isclose(pd['sums'][k][n],e['sums'][k]['old'][n]-v,rel_tol=1e-12)for k in sub for n,v in sub[k].items())
 qi,dq,_=C.integer_metric(c['matrix_a']);checks=[]
 for a,old in zip(pd['anchors'],e['anchors']):
  v=C.regenerate_pair(c,a[0],a[1])[0];lrs=RE.lrs(v,a[11],a[12],qi,dq,c['production_A_squared'],F(c['production_region_lam']),F(c['production_region_kap']))
  change,_,_=delta(a[2]/65536,a[12]-a[11]+1,lrs,g,dd)
  checks.append(abs(old['old']+change-a[13])<3e-6)
 ck['anchor_deltas']=all(checks)and len(checks)==8
 nd=sum(pd['draws'][k]for k in pd['work']);sc=float((t/F(c['production_theta']))**3)
 s=M.pool_totals(M.add_sums([pd['sums'][k]for k in pd['work']]),nd*sc,nd,g)
 col=s['log2_est']+math.log2(M.CANON_WORK_ALLOWANCE)
 ck['allowance']=M.CANON_WORK_ALLOWANCE>=1+2*s['work_se']
 ck['collection']=col<=116.865873450
 return dict(checks_ok=all(ck.values()),checks=ck,collection_log2=col,work_se_upper=s['work_se'],sampled_per_cell=s['per_cell'],D_sup_log2=math.log2(d52),B=B)

def replay():
 import ideal_canonical as I
 c,t,g=context();e=json.loads((ROOT/'pool_d.json').read_text())['extension'];pd=json.loads((ROOT/'pool_d.json').read_text())
 qi,dq,_=C.integer_metric(c['matrix_a']);qb,db,_=C.integer_metric(c['matrix_b']);A=c['production_A_squared'];units=RE.unit_list();lim=int(t*A*dq)
 out={}
 for pool,seeds in pd['pools'].items():
  terms=[];count=0;mx=-math.inf
  for seed in seeds:
   for i in range(pd['draws_per_seed']):
    a,b,aq,bq=C.regenerate_pair(c,seed,i);qa=round(F(aq,dq)/A*65536)/65536
    if aq>lim or qa>float(t)or F(bq)>(F(c['production_region_lam'])*A-F(aq,dq))*db:continue
    try:canonical=I.canonical_a(a,c['matrix_a'],100000)[0]
    except RuntimeError:canonical=I.canonical_a(a,c['matrix_a'],10**7)[0]
    if not canonical:continue
    lo,hi=RE.run(a,units,qi,lim);lrs=RE.lrs(a,lo,hi,qi,dq,A,F(c['production_region_lam']),F(c['production_region_kap']))
    d,ext,_=delta(qa,hi-lo+1,lrs,g,e['Ddiff']);terms.append(g.weight(qa)*d);count+=1;mx=max(mx,ext)
  v=math.fsum(terms);ref=e['delta'][pool]
  out[pool]=dict(rows=count,work_delta=v,extension_max=mx,ok=count==ref['rows']and mx<=0 and math.isclose(v,ref['work_delta'],rel_tol=1e-10))
 return out


def main():
	req=json.load(sys.stdin)
	if req not in({},{"descent":True},{"descent":False}):
		raise ValueError('input: {} or {"descent": true}')
	c=json.loads((ROOT/"refinement.json").read_text());ck={}
	ck["record_extension"]=extension_audit()["checks_ok"]
	r0=M.retentions()[0];sp=M.band_split(*c["split_inputs"],M.rmin(r0,M.OP['slack_f']),M.rmin(r0,M.OP['slack_g']))
	sb,sf=M.srt_toy(6,1);import base64,lzma,re;kt=re.split(r"(?m)^@@ fft.txt\n",lzma.decompress(base64.b85decode("".join((ROOT/"kernels_asm.b85").read_text().split()))).decode())[1];fp=[x.split(":")for l in kt.splitlines()if l.startswith("FFTPTS ")for x in l.split()[1:]];ck["srt"]=sb==0 and sf<=2 and len(fp)>=25 and all(float(c)<=M.fft_F(int(n))for n,c in fp)and M.gbatch_cost()<=M.C_GBATCH and M.gbatch_cost(pbits=1.000028*2**48/M.LN2,leaf=13,L=39,c_leaf=2**19)<=M.C_BATCH and (M.ecm_instr(69.0,576)*12*1.002+13*M.c_mr(576,M.MONT5))/M.UNIT*2**65.1456<2**M.LOG2_GSPLIT
	ck["band_split"]=sp["log2_units"]<=M.LOG2_SPLIT==float(c["ledger_upper_rows"][-1])
	ck["owner_schedule"]=O.ACTIVE_OWNERS==6 and O.row_tail()==1210 and O.C_ROW==1215.0 and O.LA_BOUND==float(c["ledger_upper_rows"][1])
	ck["code_identity"]=all(hashlib.sha256((ROOT/n).read_bytes()).hexdigest()==h for n,h in c["code_sha256"].items())
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
		ck["claim_budget"]=tot_log2<=Decimal(CLAIM)
		# plain bound = core meter + log2(10); declared value must be the least 12-dp value at or above it
		plain=tot_log2+Decimal(10).ln()/two.ln()
		# declared must be the least 12-dp value at or above the exact plain meter
		ck["plain_conversion"]=(plain<=Decimal(PLAIN_CLAIM)
			and Decimal(PLAIN_CLAIM)-Decimal("1e-12")<plain)
		ck["claim_declared"]=Decimal(PLAIN_CLAIM)==Decimal(c["claim"])
		ck["wrapped_success"]=s_w>Decimal("0.39")
	print(json.dumps({"checks_ok":all(ck.values()),"checks":ck,"core_ledger_upper_log2":str(core_log2),
	"markov_upper":str(markov),"core_success_lower":str(succ),"log2_upper_with_wrapper":str(tot_log2),"declared":CLAIM,"success_lower_with_wrapper":str(s_w),"descent":out or "recorded (rerun with {\"descent\": true})","band_split":sp,
	"scope":"Recorded full-run certificate plus selected recomputation; evidence16.py reruns the scientific checks; --replay-pools regenerates paired sample deltas. Conditional on H1-H5, H-EA, H-X, H-FFT."}))
	if not all(ck.values()):
		raise SystemExit(1)
if __name__=="__main__":
	if "--replay-pools" in sys.argv:
		r=replay();print(json.dumps(r))
		if not all(v["ok"]for v in r.values()):raise SystemExit(1)
	else:main()

"""Stdlib ledger/kernel audit from recorded campaigns; full descent replay."""
import base64,hashlib,json,lzma,math,re,sys
from fractions import Fraction as F
from decimal import Decimal,localcontext
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE/"files" if(HERE/"files").is_dir()else HERE
sys.path.insert(0,str(ROOT))
import descent16 as D,model16 as M,retain16 as R,owner10 as O,counts16 as T
CLAIM="117.617261684243"
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
	ck["ledger_model_rows"]=[float(v)for v in c["ledger_upper_rows"]][3:]==[M.DESCENT_CAP,M.LOG2_CANON,M.LOG2_SM,M.LOG2_MISC,M.LOG2_FBC,M.LOG2_WIDE_LAT,M.LOG2_PREP7,50,M.LOG2_VERIFY,M.LOG2_GTREE,M.LOG2_GSPLIT,M.LOG2_SPLIT,105]and M.UNIT==cfg["reference_cost_instructions"]==20000
	if M.C_ROW==2710.0:O.install(M)
	with localcontext()as x:
		x.prec=80;d=lambda a:Decimal(a.numerator)/Decimal(a.denominator)if isinstance(a,F)else Decimal(str(a))
		ln2=Decimal(2).ln();lg=lambda a:a.ln()/ln2;p2=lambda a:(d(a)*ln2).exp()
		col_log=lg(d(col));z=Decimal(69)*ln2
		U=d(T.factor_base_upper());nn=U*d(F(str(M.NU300)))
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
def _prefix_q_lower(D,L,U,d,alpha):
    base=D.P.LX;m=math.ceil(U)
    if d['n']<=base:return D.rho_lo((m+alpha)/base)*(1-2**-30)
    rows,h=D.window_rows(d['n'],base,d['iA'],D.P.MESH,d['kA'])
    tab=D.RHO_TAB;nt=len(tab);total=0.0
    for (ja,jb),row in rows.items():
        if ja>d['kA']or jb>d['kB']:continue
        j=ja+jb
        last=0 if j==0 else min(len(row)-1,math.floor((L-j*h-j*base)/h-2**-20))
        if last<0:continue
        k0=(m-j*base+alpha)/(base*D.RSTEP);hk=h/(base*D.RSTEP)
        for i in range(0,last+1,16):
            w=sum(row[i:min(i+16,last+1)])
            if w:
                k=math.ceil(k0-i*hk+2**-20)
                total+=w*(1.0 if k<=0 else tab[k]if k<nt else 0.0)
    return total*(1-2**-30)

def prefix_certificate(D,result):
    rows=[];used_u=0.0
    for (bin_,side),d in result['DET'].items():
        assert d['how']=='sieve'
        one=[];where=[]
        for name,alpha in [('nf',D.P.ALPHA_F),('ng',D.P.ALPHA_G)]:
            top=math.ceil(d[name]);best=(0.0,None)
            used_u=max(used_u,(top+alpha)/D.P.LX)
            for lo in range(0,top,64):
                hi=min(lo+64,top)
                q=_prefix_q_lower(D,lo,hi,d,alpha)
                sv=D.surv_upper(max(0,lo-1),D.P.LX,d['T'])
                value=sv/q if q else float('inf')
                if value>best[0]:best=(value,(lo,hi))
            one.append(best[0]);where.append(best[1])
        bound=(one[0]*one[1]+1)*(1+2**-30)
        charged=2**(d['surv_log2']-d['S']-d['pr_log2'])
        coverage=charged/bound
        work_upper=(1/(d['mu']*(1-d['psi']))+2**(-d['surv_log2']))*(1+2**-30)
        rows.append(dict(bin=bin_[0],side=side,coverage=coverage,
                         work_ratio=work_upper/d['lazy'],worst_intervals=where))
    end=min(len(D.RHO_TAB)-1,math.ceil(used_u/D.RSTEP)+2)
    checks=dict(interval_rate=all(r['coverage']>=42 for r in rows),
                prefix_work=all(r['work_ratio']<=1 for r in rows),
                prefix_truncation=math.log2(T.expected_nodes(D.P(),result)[0])-min(d['mu']*(1-d['psi'])for d in result['DET'].values())/math.log(2)<=-13,
                rho_monotone=all(a>=b for a,b in zip(D.RHO_TAB[:end],D.RHO_TAB[1:end+1])))
    return dict(checks=checks,cases=len(rows),raw_rate_threshold=42,distinct_fraction=1/64,good_mark_factor=1,
                minimum_coverage=min(r['coverage']for r in rows),
                maximum_work_ratio=max(r['work_ratio']for r in rows),rows=rows)

def tail(stats,ck):
	w=stats['norm_cdf']['minima'];o={}
	ck['cdf_witness_sizes']=len(w['node'])==2000 and len(w['boot'])==32
	for name,K,e in(('node',24,F(41,1000)),('boot',10,F(36,1000))):
		v=[F(c,2048)-F(19,20)*F(j,K)-F(1,16384)-e for j,c in w[name]]
		ck['cdf_'+name]=all(isinstance(j,int)and 1<=j<=K and isinstance(c,int)and 0<=c<=2048 for j,c in w[name])and min(v)>=0
		o[name]=dict(N=len(v),failures=sum(x<0 for x in v),min_margin=str(min(v)),index=v.index(min(v)))
	ck['cdf_test_levels']=29*(-Decimal(4096)*Decimal('.041')**2).exp()<Decimal('.05')and 10*(-Decimal(4096)*Decimal('.036')**2).exp()<Decimal('.05')
	ck['cdf_pooled_tail']=(1-F(19,4000))**2000<F(1,20)
	ck['cdf_boot_tail']=F(21,40)**32+F(21,40)**16<F(1,20)
	ck['uniform_ball']=all(2*D.CELLS[r]*(1<<16)>=(2*r+1)**12 for r in D.P.RADII)
	ck['audit_work']=F(2**28*(2**78+2**90+2**61)+2**47*2**69+2**110,20000)<2**105
	ck['audit_errors']=F(1,2**87)+F(1,2**88)+F(1,2**125)+F(1,2**160)+F(1,2**900)<F(1,2**80)
	br=stats['bound_calibration']['bounds62']
	ck['bound_calibration']=len(br)==62 and min(min(x[0],x[2])for x in br)>=0 and sum(sum(x[4:])for x in br)==0
	o['calibration']=dict(lattices=62,norms=1984,f_gap_max=max(x[1]for x in br),g_gap_max=max(x[3]for x in br),cdf_differences=sum(x[4]for x in br))
	ck['transport_parameters']=(F(str(D.P.NORM_MASS)),F(str(D.P.BAD_LATTICE)))==(F(19,20),F(1,200))
	return o
def replay(c):
	import counts16,time
	t=time.perf_counter()
	bad=[n for n,h in c['code_sha256'].items()if hashlib.sha256((ROOT/n).read_bytes()).hexdigest()!=h]
	ck={'code_identity':not bad,'parameters':D.P.MESH==384 and D.RSTEP==1/1024 and D.P.TRUNC==0 and D.P.LOG2_RELS==c['log2_rels_input']and D.P.BOOT_RELS==18}
	ck['windows']=all(D.PLAN16[(b,s)]==(o[i][0],8*o[i][1],*o[i][2:])for b,o in D._PLAN.items()for i,s in enumerate('fg'))and D.BOOT16==(215,328,16,4)
	manifest=json.loads((ROOT/'submission.json').read_text())['claims'][0]['score']
	ck['manifest_claim']=F(int(manifest['numerator']),int(manifest['denominator']))==F(CLAIM)==F(c['claim'])
	stats=json.loads(lzma.decompress(base64.b85decode((ROOT/'stats16.b85').read_text())));R.load_table(stats)
	derived=derive(c,ck);cdf=tail(stats,ck)
	fc=T.factor_base_certificate();rc=T.repeat_pool_certificate(tuple((int(a),int(b))for a,b in D.bins_of(D.P())))
	ck['prime_ideal_count']=all(fc['checks'].values())and fc['Ucap']==12632962157171371784
	ck['repeat_pools']=all(rc['checks'].values())and rc['pool_log2_lower']==68
	res=D.evaluate(D.P(),D.PLAN16,D.BOOT16);EN,_=counts16.expected_nodes(D.P(),res);bt=res['boot'];rows=list(res['DET'].values())
	pc=prefix_certificate(D,res);ck.update({'prefix_'+k:v for k,v in pc['checks'].items()});ck['prefix_memory']=max(d['surv_log2']for d in rows)<71
	out=dict(two_log=res['two_log_units'],boot_phi=bt['phi'],boot_psi=bt['psi'],log2_EN=math.log2(EN),mu_min=min(d['mu']for d in rows),phi_max=max(res['PHI'].values()),psi_max=max(max(d['psi']for d in rows),bt['psi']))
	ck['descent_recomputed']=out['two_log']<=c['two_log_upper']and bt['phi']<=c['boot_phi_upper']and math.log2(EN)-M.NODE_CAP_LOG2<=c['counter_abort_log2_upper']and out['mu_min']>=D.P.MIN_RELS and out['phi_max']<=0.05 and out['psi_max']<=D.P.PSI_MAX
	ck['measured_reproduced']=abs(out['two_log']-c['two_log_measured'])<1e-9 and abs(out['boot_phi']-c['boot_phi_measured'])<1e-15 and abs(out['log2_EN']-c['node_log2_measured'])<1e-9
	# Every denominator cutoff is below 512 bits; absorption does not remove
	# denominator mass. The numerator's hi(mass) is nondecreasing and <=1.
	ck['omission_absorption_safe']=all(0<D.P.SLACK-1-d['nf']/128<512 and 0<D.P.SLACK-1-d['ng']/128<512 for d in rows)
	with localcontext()as x:
		x.prec=80;two=Decimal(2);p2=lambda v:(Decimal(str(v))*two.ln()).exp()
		markov=p2(Decimal(str(c['two_log_upper']))-Decimal(str(c['descent_cap'])))
		main=1-two**-66-Decimal(str(c['cap_abort']))
		other=2*Decimal(str(c['boot_phi_upper']))+p2(c['counter_abort_log2_upper'])+p2(c['repeat_log2_upper'])+two**-80+p2(c['prefix_log2_upper'])+2/ORDER
		succ=main*(1-other-markov)-two**-48;head=main*(1-other-Decimal('1.1')*markov)-two**-48
		corelog=sum(p2(v)for v in c['ledger_upper_rows']).ln()/two.ln()
		log=(p2(c['ledger_upper_claim'])+p2(c['wrapper_log2'])).ln()/two.ln()
		ck['ledger']=corelog<=Decimal(c['ledger_rows_sum_upper'])<=Decimal(c['ledger_upper_claim'])and log<=Decimal(CLAIM)and c['descent_cap']==M.DESCENT_CAP
		ck['success']=succ>=Decimal(str(c['success_lower_claim']))-two**-48>Decimal('.39')
	print(json.dumps(dict(checks_ok=all(ck.values()),checks=ck,identity_mismatches=bad,derived=derived,cdf=cdf,prime_count=fc,repeat_pool=rc,descent=out,boot=bt,prefix=pc,nodes=[dict(bin=b[0],side=s,**d)for(b,s),d in res['DET'].items()],markov_upper=str(markov),success_lower_with_wrapper=str(succ),log2_upper_with_wrapper=str(log),elapsed_seconds=time.perf_counter()-t,scope='Full node/boot recurrence and exact numeric/count/kernel accounting recomputed. Recorded pool sums, prime structures and CDF minimum witnesses are empirical inputs; campaigns are not regenerated. Arrival transport remains H5n. No descent-output cache.')))
	if not all(ck.values()):raise SystemExit(1)

def main():
	raw=sys.stdin.read();req=json.loads(raw)if raw.strip()else{}
	if req not in({},{'descent':True},{'descent':False}):raise ValueError('input: {} or descent:boolean')
	return replay(json.loads((ROOT/'refinement.json').read_text()))
if __name__=='__main__':main()

"""rec7 run rule (proof.md section 3). {"mode":"dmax"}: d(qa,lg), slice walk minus run walk per cell, on the shipped grid;
{"mode":"replay","seed":s}: redraw a sample seed's first accepted row and test every a-value of its run."""
import json,math,os,sys
from fractions import Fraction as F
from random import Random
D=os.path.dirname(os.path.abspath(__file__));sys.path.insert(0,D)
import model16 as M
cfg=json.load(open(os.path.join(D,"config.json")));RR=json.load(open(os.path.join(D,"refinement.json")))["run_rule"]
def geo():
	import lattice_count_bounds as L
	th=F(cfg['production_theta_region']);c=L.correction(cfg,math.isqrt(int(th*cfg['production_A_squared'])))
	return M.Geometry(cfg,th,F(c['safe_yield_factor']),F(c['safe_work_population_factor']))
def rec(S,lg):
	lp=M.P_MID;t=0.0
	while lp<M.OP['lx2']:
		l2=min(lp+0.25,M.OP['lx2']);ds=math.log(l2/lp);h=math.log2(S)-0.5*(lp+l2);dp=(M.NBAR/2)*(2**M.log2_li(l2)-2**M.log2_li(lp))
		if lg is not None and M.P7<=lp<M.WIDE_FROM and h+lg>=M.REC7_TAB[-1][0]:
			t+=M.NBAR*ds*(M.step7(h+lg)+M.EXTRA7)/(1-M.RUN_ALLOW)+dp*M.C_SETUP6*(2**-lg-1)/S
		elif h>=M.REC_TAB[-1][0]:
			t+=M.NBAR*ds*(M.WIDE_REC*M.tab(M.REC_TAB_W,h)if lp>=M.WIDE_FROM else M.tab(M.REC_TAB,h))
		else:t+=dp*M.C_COSET_EMPTY/S
		lp=l2
	return t
def dmax(nq,nl,lgmax):
	g=geo();th=float(F(cfg['production_theta_region']));b=(-1,None)
	for i in range(nq+1):
		q=th*i/nq;S=g.slice_cells(q);s0=rec(S,None)
		for j in range(nl+1):
			lg=lgmax*j/nl;d=s0-rec(S,lg)
			if d>b[0]:b=(d,q,lg)
	return b
def canon(a):
	import ideal_canonical as IC
	try:return IC.canonical_a(a,cfg['matrix_a'],100000)[0]
	except RuntimeError:return IC.canonical_a(a,cfg['matrix_a'],10**7)[0]
def replay(seed,k):
	import common12 as C,runs16 as R
	pd=json.load(open(os.path.join(D,"pool_d.json")));u=R.unit_list();qi,dq,_=C.integer_metric(cfg['matrix_a']);A2=cfg['production_A_squared']
	lim=int(F(pd['theta_region'])*A2*dq);out=[];i=0
	while len(out)<k:
		a,_,aq,_=C.regenerate_pair(cfg,seed,i);x=i;i+=1
		if aq>lim or not canon(a)or Random(f"runtest-accept:{seed}:{x}").random()>=(1-aq/dq/A2)**3:continue
		lo,hi=R.run(a,u,qi,lim);bad=[t for t in range(lo,hi+1)if t and 32<=hi-lo+1<=1024 and not canon([a[0]+t]+a[1:])]
		out.append([seed,x,lo,hi,bad])
	return out
if __name__=="__main__":
	q=json.load(sys.stdin)
	if q.get("mode")=="replay":
		s,k=int(q["seed"]),int(q.get("rows",1));assert str(s)in RR["first"]and 1<=k<=2
		got=replay(s,k);ok=[g[:4]for g in got][:1]==[RR["first"][str(s)]]and all(not g[4]for g in got)
		print(json.dumps({"checks_ok":ok,"rows":got,"shipped_first":RR["first"][str(s)]}))
	else:
		d,qa,lg=dmax(*RR["grid"]);ok=abs(d-RR["d_max"])<=1e-9*RR["d_max"]and RR["d_max"]<=RR["d_charged"]
		print(json.dumps({"checks_ok":ok,"d_max":d,"at":[qa,lg],"recorded":RR["d_max"],"charged":RR["d_charged"]}))
	if not ok:raise SystemExit(1)

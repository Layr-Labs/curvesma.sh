"""stnfs GT-core experiment (recover-Q q12 = core of plain P10 9c904850; core of recover-Q 9332a74d). Stdlib only; stdout"""
import base64,hashlib,json,lzma,math,os,re,sys,time
from decimal import Decimal,localcontext
from fractions import Fraction as F
from random import Random
_HERE=os.path.dirname(os.path.abspath(__file__))
_DIR=os.path.join(_HERE,"files")if os.path.exists(os.path.join(_HERE,"files","model16.py"))else _HERE
sys.path.insert(0,_DIR)
import common12 as C
import model16 as M
import owner10 as O
O.install(M)
import descent16 as DS
import store16
from threshold_model import Retention
import ideal_canonical as IC
import lattice_count_bounds as LCB
import filtsim
P=4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787
U=C.U
H=C.H
LN2=math.log(2)
CLAIM=119.4568750535
import retain16
ALPHA_PMAX=2000
ALPHA_F,ALPHA_G=4.868096,0.843643
ALPHA_G_PMAX=50000
ALPHA_G_50K=0.731543
ALPHA_F_50K=5.085035
def alpha_contrib(p):
	vf=vg=0.0
	for d1,e,mf,mg in retain16.structure(p):
		q=float(p**d1)
		for side,mults in enumerate((mf,mg)):
			v=0.0
			for mu in mults:
				if mu==1:v+=q/ (q*q-1.0)
				elif p==2 and side==1 and mu==2:v+= (2.0*q/ (q*q-1.0)) * (1.0+q** -8)
				else:v+=1.0/ (q+1.0)
			if side==0:vf+=d1*v
			else:vg+=d1*v
	t=math.log(p) /math.log(2.0)
	return t* (1.0/ (p-1.0) -vf),t* (1.0/ (p-1.0) -vg)
def alpha_sweep(pmax):
	af=ag=0.0
	for p in retain16.primes(2,pmax):
		a,b=alpha_contrib(p);af+=a;ag+=b
	return af,ag
DECLARED= {'ledger': ('<=',119.4568750535),'relcol': ('<=',119.190851442),'la': ('<=',116.775787663),'gen': ('<=',107.389731954),'log2_rels': ('>=',63.6016),'rels_1se': ('>=',63.58379),'margin': ('>=',0.0798),'work_se': ('<=',0.002359),'two_log': ('<=',112.3879),'markov': ('<=',0.608726),'log2_EN': ('<=',15.6383),'over': ('<=', -11.66),'rep': ('<=', -10.7),'boot_phi': ('<=',9.363e-06),'success': ('>=',0.39034)}
_REPORT= []
T0=time.time()
def load_stats():
	import lzma,base64
	return json.loads(lzma.decompress(base64.b85decode(open(os.path.join(_DIR,"stats16.b85")).read())))
def say(*args):
	line=" ".join(str(a)for a in args)
	_REPORT.append(line)
	print(line,file=sys.stderr,flush=True)
def f6_mulp(A,B):
	return[c%P for c in C.mul([c%P for c in A], [c%P for c in B])]
def f6_pow(A,n):
	r,base= [1], [c%P for c in A]
	while n:
		if n&1:
			r=f6_mulp(r,base)
		base=f6_mulp(base,base)
		n>>=1
	return(r+ [0] *6)[:6]
def check_resultant():
	F0= [0] *6
	for d,c in C.TERMS:
		term= (c+ [0] *6)[:6]
		for _ in range(d//2):
			term=C.mul(term, [0, -U])
		F0= [x+y for x,y in zip(F0,term)]
	r=C.norm(F0)
	ok=r!=0 and r%P==0
	say(f"[0] Res_y(h, f(x^2 -> uy)): {r.bit_length()} bits, 0 mod p {'ok' if ok else 'FAIL'}")
	return ok
def check_tower():
	yp=f6_pow([0,1],P)
	cur,ok_irr= [0,1,0,0,0,0],True
	for k in range(1,7):
		acc= [0] *6
		for coef in reversed(cur):
			acc=f6_mulp(acc,yp)
			acc[0] = (acc[0] +coef) %P
		cur=acc
		is_y=cur== [0,1,0,0,0,0]
		if(k in(1,2,3)and is_y)or(k==6 and not is_y):
			ok_irr=False
	uy=f6_mulp([(-U) %P], [0,1])
	chi=f6_pow(uy, (P**6-1) //2)
	ok_sq=chi== [P-1,0,0,0,0,0]
	say(f"[0b] h irreducible mod p {'ok' if ok_irr else 'FAIL'}; (u y)^((p^6-1)/2) = -1 {'ok' if ok_sq else 'FAIL'}")
	cert=IC.field_certificate(P)
	ok_c=cert['squarefree']and cert['irreducible_mod_p']and cert['maximal_power_basis']and cert['discriminant'] ==58673
	say(f"[0c] Disc(h) = {cert['discriminant']} ({cert['discriminant_factors']}), O_h = Z[y] {'ok' if ok_c else 'FAIL'}")
	return ok_irr and ok_sq and ok_c
def check_rho_and_semismooth():
	ref= {2:0.30685281944,3:0.048608388,4:0.0049109256,6:1.9649696e-5,8:3.2320693e-8,10:2.7701718e-11}
	err=max(abs(DS.rho_series(u) -v) /v for u,v in ref.items())
	ok_br=all(DS.rho_lo(u) <=DS.rho_series(u) <=DS.rho_hi(u) +1e-300 for u in(1.5,2.3,7.7,15.2,33.3,41.0))
	rows,h=DS.window_rows(80.0,69.0,48,48)
	ok_w=abs(sum(rows[(1,0)]) -math.log(80/69)) <1e-12
	say(f"[1] rho: series vs reference {err:.1e} {'ok' if err < 1e-6 else 'FAIL'}; brackets {'ok' if ok_br else 'FAIL'}; weight sum = ln(n/base) "
	 f"{'ok' if ok_w else 'FAIL'}")
	return err<1e-6 and ok_br and ok_w
def canonical(cfg,a):
	try:
		return IC.canonical_a(a,cfg['matrix_a'],100000)[0]
	except RuntimeError:
		return IC.canonical_a(a,cfg['matrix_a'],10**7)[0]
def verify_pool(cfg,pd):
	t0=time.time();bad= [];rows= []
	import runs16 as R16
	units=R16.unit_list();qi,dq,_=C.integer_metric(cfg['matrix_a']);lim=int(F(pd['theta_region']) *cfg['production_A_squared'] *dq)
	for an in pd['anchors']:
		sd,i=an[0],an[1];r=M.pool_row(cfg,sd,i);a=C.regenerate_pair(cfg,sd,i)[0]
		f= [round(r['qa'] *65536),round(r['nf'] *16),round(r['ng'] *16),r['primitive'] +2*r['rejected'],r['evals'],r['bcells'],r['leaves'],r['blen'],r['rlen'] //16] +list(R16.run(a,units,qi,lim))
		if f!=an[2:]or not canonical(cfg,a):
			bad.append(('anchor',sd,i))
		rows.append(r)
	n_can=0
	for sd, (hi,idx)in sorted(pd['replay'].items()):
		got= [i for i in range(hi)if canonical(cfg,C.regenerate_pair(cfg,int(sd),i)[0])]
		n_can+=len(got);bad+= [('replay',sd)]if got!=idx else[]
	say(f"[2] {len(pd['anchors'])} judged-pool anchors (canonical, all fields, run ends) and draws 0..{hi - 1} of seeds {', '.join(sorted(pd['replay']))} "
	 f"({n_can} canonical) replay the shipped lists ({time.time()-t0:.1f}s): {'ok' if not bad else 'FAIL ' + str(bad[:3])}")
	return not bad,n_can,rows,len(units)
def box_check(cfg,pool,k=6,cells=8,seed=5):
	import embedding_bounds as EB
	rng=Random(f"box-check:{seed}")
	qa,da,_=C.integer_metric(cfg['matrix_a']);qb,db,_=C.integer_metric(cfg['matrix_b'])
	A2=cfg['production_A_squared'];disks=EB.root_disks();n=bad=nrow=0
	for r in rng.sample(pool,k):
		a,b,aq,bq=C.regenerate_pair(cfg,r['seed'],r['index'])
		est=EB.NormBounds(a,disks)
		box=M.box_rows(qb, (A2-F(aq,da)) *db,b, (b[0] //64) *64);nrow+=len(box)
		leaves=EB.box_leaves(est, [(est.row(list(s)),l,h)for s,l,h in box],min(x[1]for x in box),max(x[2]for x in box))
		for kk, (s,l,h)in enumerate(box):
			for b0 in rng.sample(range(l,h+1),min(cells,h-l+1)):
				lf= [x for x in leaves if x[0] <=b0<=x[1]and x[3]in(None,kk)][0]
				if isinstance(lf[2],dict):
					continue
				for(lo16,hi16),v in zip(lf[2],C.norms(a, [b0] +list(s))):
					n+=1;bad+=not(v!=0 and lo16-1e-6<=16*math.log2(abs(v)) <=hi16+1e-6)
	return bad==0,n,nrow
def poly_roots_count(poly,p):
	poly= [c%p for c in poly]
	while poly and poly[-1] ==0:
		poly.pop()
	if len(poly) <=1:
		return 0
	def pmod(a,m):
		a=a[:]
		inv=pow(m[-1], -1,p)
		while len(a) >=len(m):
			if a[-1]:
				z=a[-1] *inv%p;s=len(a) -len(m)
				for j,v in enumerate(m):
					a[s+j] = (a[s+j] -z*v) %p
			a.pop()
			while a and a[-1] ==0:
				a.pop()
		return a
	def pmulmod(a,b,m):
		r= [0] * (len(a) +len(b) -1)
		for i,x in enumerate(a):
			if x:
				for j,y in enumerate(b):
					r[i+j] = (r[i+j] +x*y) %p
		return pmod(r,m)
	r,base,e= [1], [0,1],p
	while e:
		if e&1:
			r=pmulmod(r,base,poly)
		base=pmulmod(base,base,poly);e>>=1
	r= (r+ [0,0])[:max(len(r),2)];r[1] = (r[1] -1) %p
	a,b=poly,pmod(r,poly)if len(r) >=len(poly)else[c%p for c in r]
	while b and any(b):
		a,b=b,pmod(a,b)
	return len(a) -1
def nbar(lo_bits,hi_bits):
	from threshold_model import primes_upto
	ps= [q for q in primes_upto(2**hi_bits)if q>=2**lo_bits]
	tot=0
	rng=Random(1)
	for p in ps:
		for s in sorted(set(roots_modp(H,p,rng))):
			tot+=poly_roots_count([1+s-s*s-s**3,0,s**5,0,0,0,2*s**3,0,0,0, -2*s,0,1],p) + \
    poly_roots_count([U*s,0,1],p)
	return tot/len(ps),len(ps)
def _asm(name):
	t=lzma.decompress(base64.b85decode("".join(open(os.path.join(_DIR,"kernels_asm.b85")).read().split()))).decode()
	return t.split("@@ "+name+"\n")[1].split("\n@@ ")[0].splitlines()
def kernel_rec7(ru,nunits):
	L=_asm("rec7.txt");T=M.REC7_TAB
	pts= [(math.log2(float(f[-1])), (int(f[9]) -int(f[8])) / (int(f[11]) -int(f[10])))for f in(ln.split()for ln in L if ln.startswith("CAMP "))]
	ok=len(pts) ==130 and all(a[0] >b[0]and a[1] <=b[1]for a,b in zip(T,T[1:]))and all(lh>=T[-1][0]and c<=M.step7(lh)for lh,c in pts)
	chk= [ln for ln in L if ln.startswith("CHECK ")];ok&=len(chk) ==1 and "ok (0 mismatches over 270 cosets" in chk[0]
	pm=max(int(ln.split()[4])for ln in L if ln.startswith("PREP "));ideals=M.NBAR* (2**M.log2_li(64) -2**M.log2_li(52))
	ok&=M.LOG2_PREP7>=math.log2(2*ideals*pm*2**0.05/M.UNIT)and M.P7>=52 and M.EXTRA7>=2
	Ts= [int(ln.split()[3])for ln in L if ln.startswith("CAMP ")];ok&=min(Ts) <=M.RUN_T[0]and M.RUN_T[1] <=max(Ts)
	ok&=nunits==ru['units']and ru['exact_bad'] ==0 and M.RUN_ALLOW>=1-0.05** (1/ru['exact_rows'])
	rn= {int(p):F(int(b) -int(a),int(d) -int(c))for p,a,b,c,d in re.findall(r"record p=(\d+): Ir (\d+) -> (\d+), hits (\d+) -> (\d+)","\n".join(_asm("measurements12.txt")))}
	ok&= [float(rn[p])for p in sorted(rn)] == [x[1]for x in M.REC_TAB]
	return ok, (f"runs of {M.RUN_T[0]}..{M.RUN_T[1]} a-values inside the campaign's {min(Ts)}..{max(Ts)}; {len(pts)} campaign points at or below the step table (worst ratio {max(c / M.step7(lh) for lh, c in pts):.4f}); {chk[0][6:]}; "
	 f"prepare7 <= {pm} instructions per ideal, 2 bases for 2^{math.log2(ideals):.2f} ideals: 2^{M.LOG2_PREP7}; runs: {nunits} units, exact sample "
	 f"{ru['exact_rows']} runs / {ru['exact_points']} a-values, {ru['exact_bad']} non-canonical, allowance {M.RUN_ALLOW} >= {1 - 0.05 ** (1 / ru['exact_rows']):.4f}; "
	 f"REC_TAB == measurements12.txt quotients {[round(float(rn[p]), 3) for p in sorted(rn)]}")
def contam_check(cfg,pd):
	import runs16 as R16
	c=pd['runs']['contam'];th=F(pd['theta_region']);A2=cfg['production_A_squared'];qi,dq,_=C.integer_metric(cfg['matrix_a'])
	U=R16.unit_list();lim=int(th*A2*dq);geo=M.Geometry(cfg,th);ok=True
	def band(r,S):
		lg=math.log2(r['gain'])if M.RUN_T[0] <=r['T'] <=M.RUN_T[1]else None;lp=M.P7;x=0.0
		while lp<M.WIDE_FROM:
			lp2=min(lp+.25,M.WIDE_FROM);ds=math.log(lp2/lp);lh=math.log2(S) -(lp+lp2) /2;dp=M.NBAR/2* (2**M.log2_li(lp2) -2**M.log2_li(lp))
			s_=M.NBAR*ds*M.tab(M.REC_TAB,lh)if lh>=M.REC_TAB[-1][0]else dp*M.C_COSET_EMPTY/S
			w_=M.NBAR*ds* (M.REC7_TAB[-1][1] +M.EXTRA7) / (1-M.RUN_ALLOW)
			m_=M.NBAR*ds* (M.step7(lh+lg) +M.EXTRA7) / (1-M.RUN_ALLOW) +dp*M.C_SETUP6* (2** -lg-1) /S if lg is not None and lh+lg>=M.REC7_TAB[-1][0]else s_
			x+=max(s_,w_) -m_;lp=lp2
		return x
	for sd,i,m0,c0 in c['replay']:
		r=M.pool_row(cfg,sd,i);a=C.regenerate_pair(cfg,sd,i)[0];lo,hi=R16.run(a,U,qi,lim);r['T'] =hi-lo+1;r['gain'] =R16.gain(a,lo,hi,qi,dq,A2)
		m=M.evaluate([r,dict(r,primitive=True,rejected=False,qa=max(r['qa'],.03))],1,geo,return_work=True)['work_contribs'][0];S=geo.slice_cells(r['qa'])
		ok&=math.isclose(m,m0,rel_tol=1e-9)and math.isclose(m+S*band(r,S),c0,rel_tol=1e-9)
	for sd,i,lo,hi in c['checked']:
		a=C.regenerate_pair(cfg,sd,i)[0];ok&=R16.run(a,U,qi,lim) == (lo,hi)and all(canonical(cfg, [a[0] +t] +a[1:])for t in range(lo,hi+1))
	ok&=c['bad'] ==0 and c['eps'] >=1-0.05** (1/c['draws'])and c['r_up'] >=c['R'] +2*c['se']and c['factor'] >=c['eps'] *c['r_up']
	return ok,c
def check_la_kernel(name="la_mac.s",start="LBB9_3:"):
	on,body=False, []
	for ln in _asm(name):
		s=ln.strip()
		if s.startswith(start):
			on=True;continue
		if on:
			if s.startswith(("LBB",".LBB"))or(s.startswith(";")and "bb." in s):
				break
			if s and not s.startswith(";")and not s.startswith("."):
				body.append(s.split()[0])
	return len(body), [m for m in body if m.startswith("b.")or m=="b"]
def check_la_func(fn):
	L= [x.strip()for x in _asm("la_mac.s")]
	i=L.index(fn+":") +1;n,seen,fwd,calls=0,set(),True, []
	while i<len(L)and not L[i].startswith(("_",";")):
		x=L[i];i+=1
		if x.startswith("LBB"):seen.add(x[:-1])
		elif x:
			n+=1;fwd=fwd and x.split()[-1]not in seen
			if x.split()[0] =="bl":calls.append(x.split()[-1])
	return n,fwd,calls
def la_row_tail():
	f= {k:check_la_func(k)for k in("_send_pair","_combine","_finish","_reduce9","_fold9","_reduce6")}
	ok=all(v[1]for v in f.values())and f["_finish"][2] == ["_fold9","_reduce6"]and f["_reduce9"][2] == ["_fold9","_reduce6"] \
  and not any(f[k][2]for k in("_send_pair","_combine","_fold9","_reduce6"))
	sub=f["_fold9"][0] +f["_reduce6"][0]
	row=O.row_tail(f["_send_pair"][0],f["_combine"][0],f["_finish"][0],f["_fold9"][0],f["_reduce6"][0])
	red9=f["_reduce9"][0] +sub+2
	return row,red9, {k:v[0]for k,v in f.items()},ok
def ball_counts_check():
	import functools
	@functools.lru_cache(None)
	def Nd(d,s):
		if s<0:
			return 0
		if d==0:
			return 1
		r=math.isqrt(s)
		return sum(Nd(d-1,s-x*x)for x in range(-r,r+1))
	ok=True
	for R in(4,6,8):
		n12=Nd(12,R*R)
		ok&= (n12-1) //2==DS.CELLS[R]
	V12=math.pi**6/720
	for R,half in DS.CELLS.items():
		n=2*half+1
		ok&=V12* (R-math.sqrt(3)) **12<=n<=V12* (R+math.sqrt(3)) **12
	return ok
def roots_modp(poly,p,rng):
	poly= [c%p for c in poly]
	while poly and poly[-1] ==0:
		poly.pop()
	if len(poly) <=1:
		return[]
	def pmod(a,m):
		a=a[:]
		inv=pow(m[-1], -1,p)
		while len(a) >=len(m):
		  if a[-1]:
		    z=a[-1] *inv%p;s=len(a) -len(m)
		    for j,v in enumerate(m):
		      a[s+j] = (a[s+j] -z*v) %p
		  a.pop()
		  while a and a[-1] ==0:
		    a.pop()
		return a
	def pmul(a,b,m):
		r= [0] * (len(a) +len(b) -1)
		for i,x in enumerate(a):
		  if x:
		    for j,y in enumerate(b):
		      r[i+j] = (r[i+j] +x*y) %p
		return pmod(r,m)
	def powmod(a,e,m):
		r= [1]
		while e:
		  if e&1:
		    r=pmul(r,a,m)
		  a=pmul(a,a,m);e>>=1
		return r
	xp=powmod([0,1],p,poly)
	xp=xp+ [0] *max(0,2-len(xp));xp[1] = (xp[1] -1) %p
	a,b=poly,pmod(xp,poly)
	while b and any(b):
		a,b=b,pmod(a,b)
	if not a or not any(a):
		return[]
	inv=pow(a[-1], -1,p);g= [x*inv%p for x in a]
	def split(g):
		if len(g) <=1:
		  return[]
		if len(g) ==2:
		  return[(-g[0] *pow(g[1], -1,p)) %p]
		while True:
		  c= [rng.randrange(p)for _ in range(len(g) -2)] + [rng.randrange(1,p),1]
		  t=powmod(c, (p-1) //2,g)
		  t=t+ [0] *max(0,1-len(t));t= [(t[0] -1) %p] +t[1:]
		  h=pmod(t,g)if len(t) >=len(g)else t
		  aa,bb=g,h
		  while bb and any(bb):
		    aa,bb=bb,pmod(aa,bb)
		  if aa and 1<len(aa) <len(g):
		    iv=pow(aa[-1], -1,p);hh= [x*iv%p for x in aa]
		    qq=g[:];v= [0] * (len(g) -len(hh) +1)
		    for k in range(len(v) -1, -1, -1):
		      v[k] =qq[k+len(hh) -1]
		      for j,x in enumerate(hh):
		        qq[k+j] = (qq[k+j] -v[k] *x) %p
		    assert not any(qq)
		    return split(hh) +split(v)
	return sorted(split(g))
def repeat_bound(dres):
	import counts16
	EN,_=counts16.expected_nodes(DS.P(),dres)
	pmin=min(math.log2(2* (2**DS.log2li(b[1]) -2**DS.log2li(b[0])))for b in dres['bins'])
	over=math.ceil((math.log2(EN) -M.NODE_CAP_LOG2) *100) /100
	coll=math.ceil((2*M.NODE_CAP_LOG2-1+6-pmin) *100) /100
	return math.log2(EN),over,coll,pmin
def audit(rng_seed=31001,npr=None,bpr=None):
	import descent_witness as W
	assert(W.P**6-1) %W.R!=0 and((W.P**12-1) //W.R) %W.R!=0
	rng=Random(rng_seed)
	while True:
		x= ([rng.randrange(W.P)for _ in range(6)], [rng.randrange(W.P)for _ in range(6)])
		if any(x[0] +x[1]):
			break
	w0=W.pow12(x, (W.P**12-1) //W.R)
	tt=rng.randrange(1,W.R)
	z=W.pow12(w0,tt)
	B,swaps=W.lll(W.boot_basis(z))
	bst=W.boot_ball_stats(z,B,rng_seed+1,points=6,radii=(96,))
	boot_bins_max=max(x[1]for x in DS.P.BOOT_BINS)
	sq= []
	Pm=DS.P()
	cfmax=max(b[1]for b in Pm.QBINS);cgmax=max(b[2]for b in Pm.QBINS)
	env_ok=True
	for bits,side,seed in((260,'f',13300),):
		r=W.specialq_ball_stats(bits,side,seed,points=6,radii=(96,))
		for st_ in r['stats']:
			env_ok&=st_['cf_max'] <=cfmax and st_['cg_max'] <=cgmax
		sq.append(dict(bits=bits,side=side,seed=seed,swaps=r['swaps'],stats=[{k: (round(v,2)if isinstance(v,float)else v)for k,v in s.items()}for s in r['stats']]))
	q_,_,_,_,nB,_=W.specialq_lattice(225,'g',Random(rng_seed+31))
	n_,iA_,kA_,kB_,R_,_=DS.PLAN16[(220,'g')]
	tr=W.orbit_triples(nB,rng_seed+32,R_,12)
	lq_=math.log2(q_);sn=0.0
	for lnf,lng,m,_h in tr:
		pf=DS.ss_lower2(lnf,Pm.LX,float(n_),iA_,kA_,kB_,Pm.MESH,Pm.ALPHA_F,Pm.TRUNC)
		pg=DS.ss_lower2(lng-lq_,Pm.LX,float(n_),iA_,kA_,kB_,Pm.MESH,Pm.ALPHA_G,Pm.TRUNC)
		sn+=pf*pg/m
	n1,iAb,kAb,_=DS.BOOT16
	tb=W.orbit_triples(B,rng_seed+33,4,24,boot=True)
	sb=sum(DS.ss_lower2(lng,Pm.LX,float(n1),iAb,kAb,1,Pm.MESH,Pm.ALPHA_G,Pm.TRUNC) /m for _,lng,m,_h in tb)
	orb=dict(node_m_max=max(t[2]for t in tr),node_height_max=max(t[3]for t in tr),
	     boot_m_max=max(t[2]for t in tb),boot_height_max=max(t[3]for t in tb))
	srng=Random(rng_seed+77)
	dp= [1] + [0] *16
	for _ in range(12):
		ndp= [0] *17
		for s in range(17):
			for v in range(-4,5):
				if s+v*v<=16:
					ndp[s+v*v] +=dp[s]
		dp=ndp
	tot=sum(dp[1:]);sh= [0] *17;okshape=True
	for _ in range(600):
		c=W.int_ball_point(srng,4)
		n2=sum(x*x for x in c)
		okshape&=0<n2<=16
		for x in reversed(c):
			if x:
				okshape&=x>0
				break
		sh[n2] +=1
	chi2=sum((sh[s] -600*dp[s] /tot) **2/ (600*dp[s] /tot)for s in range(1,17))
	return dict(boot=dict(swaps=swaps,stats=[{k: (round(v,2)if isinstance(v,float)else v)for k,v in s.items()}for s in bst],
	      within_model_bins=all(s['cg_max'] <=boot_bins_max for s in bst)),
	   specialq=sq,specialq_within_bins=env_ok,
	   relw=dict(node_yield_ratio=(sn/12) /2.0**npr,boot_yield_ratio=(sb/24) /2.0**bpr, **orb),
	   sampler=dict(chi2=chi2,ok=bool(okshape)and chi2<60),
	   storage=store16.check())
def check_fields():
	import retain16 as RT,descent_witness as W
	R=W.R;Hr= [x%R for x in H]
	def res(a,b):
		a=RT._trim([x%R for x in a]);b=RT._trim([x%R for x in b]);s=1
		while b and len(b) >1:
			da,db=len(a) -1,len(b) -1;rem=RT.pmod(a,b,R)
			if not rem:
				return 0
			s=s* (-1) ** (da*db) *pow(b[-1],da- (len(rem) -1),R) %R;a,b=b,rem
		return s*pow(b[0],len(a) -1,R) %R if b else 0
	def norm_poly(terms,deg):
		n=6*deg+1;xs=list(range(n))
		ys= [res(Hr, [sum(c[i] *pow(k,d,R)for d,c in terms if i<len(c)) %R for i in range(6)])for k in xs]
		for j in range(1,n):
			for i in range(n-1,j-1, -1):
				ys[i] = (ys[i] -ys[i-1]) *pow(xs[i] -xs[i-j], -1,R) %R
		poly= [ys[-1]]
		for i in range(n-2, -1, -1):
			poly= [(-xs[i] *poly[0] +ys[i]) %R] + [(poly[k-1] -xs[i] *poly[k]) %R for k in range(1,len(poly))] + [poly[-1]]
		return poly
	out= {}
	for name,terms,deg in(('f',C.TERMS,12), ('g',RT.G_TERMS,2)):
		poly=norm_poly(terms,deg)
		out['deg_'+name] =len(poly) -1
		out['sqf_'+name] =poly[-1] ==1 and len(poly) ==6*deg+1 and len(RT.pgcd(poly,RT.pderiv(poly,R),R)) ==1
	out['unramified'] =out['sqf_f']and out['sqf_g']
	out['disc_h'] =C.norm([-1,2,3,0,0,6])
	out['sig'] =signatures()
	return out
def signatures():
	def tr(p):
		while p and p[-1] ==0:p=p[:-1]
		return p
	def sturm(p):
		s= [tr([F(c)for c in p])];s.append([i*c for i,c in enumerate(s[0])][1:])
		while len(s[-1]) >1:
			a,b=s[-2][:],s[-1]
			while len(a) >=len(b):
				q=a[-1] /b[-1];d=len(a) -len(b);a=tr([x-q*b[i-d]if i>=d else x for i,x in enumerate(a)][:-1])
			if not a:break
			s.append([-x for x in a])
		return s
	def ch(s,x):
		v= [t for t in(sum(c*x**i for i,c in enumerate(q))for q in s)if t]
		return sum((a>0) != (b>0)for a,b in zip(v,v[1:]))
	nre=lambda p,B=F(10) **60: (lambda s:ch(s, -B) -ch(s,B))(sturm(p))
	sh=sturm(H);ys= []
	def iso(a,b):
		n=ch(sh,a) -ch(sh,b)
		if n==1 and b-a<F(1,10**30):ys.append((a+b) /2)
		elif n:iso(a, (a+b) /2);iso((a+b) /2,b)
	iso(F(-10),F(10))
	rf=rfp=rg=0
	for y in ys:
		fx= [0] *13
		for d,c in C.TERMS:fx[d] =sum(F(ci) *y**i for i,ci in enumerate(c))
		rf+=nre(fx);rfp+=nre(fx[0::2]);rg+=nre([U*y,0,1])
	rk=lambda n,r1:r1+ (n-r1) //2-1
	return dict(f=(rf,36-rf//2),fp=(rfp,18-rfp//2),g=(rg,6-rg//2),odd=rk(72,rf) -rk(36,rfp) +rk(12,rg) -rk(6,len(ys)))
def main():
	cfg_in=json.load(sys.stdin)if not sys.stdin.isatty()else{}
	seed=int(cfg_in.get("seed",1))
	say(f"evidence: stnfs GT core of bls12-381-fixed-g1-invert-g2 (q12); seed={seed}; unit {M.UNIT:.0f} instructions")
	ok0=check_resultant();ok0b=check_tower();ok1=check_rho_and_semismooth()
	cfg=json.load(open(os.path.join(_DIR,"config.json")))
	pd=json.load(open(os.path.join(_DIR,"pool_d.json")))
	ok2,n_can,anchors,nunits=verify_pool(cfg,pd)
	t0=time.time();ok2b,nb2,nr2=box_check(cfg,anchors,seed=seed)
	say(f"[2b] box bounds: 6 rebuilt boxes ({nr2} rows, hull of the row contexts, one bisection): {nb2} exact norms at sampled cells of every box row inside their leaf bounds {'ok' if ok2b else 'FAIL'} ({time.time()-t0:.1f}s)")
	theta_pool=F(cfg['production_theta']);theta_reg=F(cfg['production_theta_region'])
	assert theta_reg<=theta_pool==M.OP['theta']and F(pd['theta_region']) ==theta_reg
	ND=pd['draws_per_seed'];sc=float((theta_reg/theta_pool) **3)
	def region(th):
		rmin=math.isqrt(int(th*cfg['production_A_squared']));corr=LCB.correction(cfg,rmin)
		return rmin,corr,M.Geometry(cfg,th,F(corr['safe_yield_factor']),F(corr['safe_work_population_factor']))
	rmin,corr,geo=region(theta_reg)
	def tot(seeds):
		n=sum(pd['draws'][x]for x in seeds)
		return M.pool_totals(M.add_sums([pd['sums'][x]for x in seeds]),n*sc,n,geo)
	tD,tW,tB=tot(pd['judged']),tot(pd['work']),tot(pd['sizing'])
	th0=theta_reg-F(1,8000);nB=sum(pd['draws'][x]for x in pd['sizing'])
	tB0=M.pool_totals(pd['sizing_prev'],nB*float((th0/theta_pool) **3),nB,region(th0)[2])
	op=M.OP
	retain=Retention(int(op['lb']),int(op['lx2']),step=.25,max_bits=850,tilt=.05)
	cost=Retention(850,int(op['lx2']),step=.25,max_bits=851,tilt=.005)
	say(f"[3] collected region theta={cfg['production_theta_region']} (uniform sub-region of the pools sampled at theta={theta_pool}), a-radius {rmin}; "
	 f"eligible a 2^{math.log2(geo.va):.3f}, largest half-slice 2^{math.log2(geo.vb):.3f} cells; exact certified count products "
	 f"{float(corr['safe_yield_factor']):.9f} (yield), {float(corr['safe_work_population_factor']):.9f} (work) at every radius >= {rmin}; "
	 f"pools {pd['pools']} of {ND} draws per seed: sizing {pd['sizing']}, judged {pd['judged']}, work {pd['work']}; judged-pool region rows {sum(pd['sums'][x]['rows'] for x in pd['judged'])}, "
	 f"work region rows {sum(pd['sums'][x]['rows'] for x in pd['work'])}")
	log2_U=math.log2((M.NBAR/2) *2**M.log2_li(op['lb']) +2**31);need=log2_U+math.log2(op['excess'])
	e=dict(tD,log2_U=log2_U,log2_need=need,margin=tD['log2_rels'] -need,log2_relcol=tW['log2_est'] +math.log2(M.CANON_WORK_ALLOWANCE),
	   per_cell=tW['per_cell'] *M.CANON_WORK_ALLOWANCE,log2_population=tW['log2_population'])
	e['feasible'] =e['margin'] >=max(M.MARGIN,M.MARGIN_SE*tD['yield_se_log2'])
	rule=lambda t:M.MARGIN+math.sqrt(2) *t['yield_se_log2']
	ok_size=tB['log2_rels'] -need>=rule(tB)and tB0['log2_rels'] -need<rule(tB0)
	say(f"[4s] sizing on pool B: margin {tB['log2_rels'] - need:.6f} >= 0.06 + sqrt(2) SE_B = {rule(tB):.6f} at {cfg['production_theta_region']}; "
	 f"{tB0['log2_rels'] - need:.6f} < {rule(tB0):.6f} at {th0 * 8000}/8000: smallest grid point {'ok' if ok_size else 'FAIL'}")
	work_se=tW['work_se'];cf=pd['runs']['contam']['factor'];allow_need= (1+M.WORK_SE_K*work_se) * (1+cf)
	ok_allow=M.CANON_WORK_ALLOWANCE>=allow_need
	say(f"[4a] work estimator on pools {pd['work']}: conditional row SE {100*tW['cond_se']:.4f}%, total relative SE (incidence included) {100*work_se:.4f}%; "
	 f"charged allowance {M.CANON_WORK_ALLOWANCE:.4f} >= (1 + {M.WORK_SE_K:.2f} SE) x (1 + run contamination {cf}) = {allow_need:.6f} {'ok' if ok_allow else 'FAIL'}")
	ok_relmeans=e['nf_rel'] >=69*DS.P.MAIN_UF and e['ng_rel'] >=69*DS.P.MAIN_UG
	rels_lo=e['log2_rels'] +math.log2(1-e['yield_se_log2'] *math.log(2))
	ok_rels_input=rels_lo>=DS.P.LOG2_RELS
	rf=retain.at(e['nf_mean'],op['slack_f']);rg=retain.at(e['ng_mean'],op['slack_g'])
	say(f"[4] judged pool {pd['judged']}: B=2^{op['lb']:.0f}, X2=2^{op['lx2']:.0f}, slacks ({op['slack_f']}, {op['slack_g']}), excess {op['excess']}, "
	 f"width {op['width_bits']}, blocks of {op['L0']}, alpha ({op['alpha_f']}, {op['alpha_g']}); mean norms "
	 f"{e['nf_mean']:.1f}/{e['ng_mean']:.1f} (relation-weighted {e['nf_rel']:.2f}/{e['ng_rel']:.2f}); retention at the means {rf['retention']:.4f}/{rg['retention']:.4f}; "
	 f"rels 2^{e['log2_rels']:.5f} vs need 2^{e['log2_need']:.5f}, margin {math.floor(e['margin']*10000)/10000:.4f} >= max({M.MARGIN}, {M.MARGIN_SE:g} SE) "
	 f"({'feasible' if e['feasible'] else 'INFEASIBLE'}); pool SE {e['yield_se_log2']:.5f} bits; "
	 f"descent LOG2_RELS input {DS.P.LOG2_RELS:.3f} <= 1-SE lower count {rels_lo:.5f} {'ok' if ok_rels_input else 'FAIL'}")
	la=M.linear_algebra(e['log2_U'],op['excess'])
	say(f"[5] population 2^{e['log2_population']:.3f} cells at {e['per_cell']:.1f} instructions per cell (routing {M.REC_ROUTE:.0f}/record, "
	 f"{M.CELL_ROUTE:.0f}/cell; apply {M.C_APPLY:.2f}, recovery {M.C_RECOVER:.3f}): relation collection 2^{e['log2_relcol']:.4f} units; "
	 f"charged N(300)/U {M.NU300} -> N = 2^{la['log2_N']:.4f}; linear algebra 2^{la['log2_la']:.4f} units ({M.C_ENTRY:.0f} instr/entry)")
	t0=time.time()
	Us,res=filtsim.run(op['lb'],e['nf_mean'] /op['lb'],e['ng_mean'] /op['lb'],op['excess'],
	      seed,12,kmax=64,targets=(100,),verbose=False)
	r100=res[100]
	say(f"[6] filtering (scale 2^12, seed {seed}, {time.time()-t0:.1f}s): N/U = {1/r100[0]:.4f} at weight {r100[3]:.0f} (max |coefficient| {r100[2]}); "
	 f"shipped campaign runs per scale: "+", ".join(f"2^{k} {len(v)}" for k,v in sorted(M.FILT_RUNS.items())))
	t0=time.time()
	rec=filtsim.missing_recon(seed,excess=op['excess'])
	ok_rec=all(mf<=md for(lo,hi,n,mf,md,h3)in rec)
	say(f"[6b] missing logs after purge and back-substitution (scale 2^12, seed {seed}, {time.time()-t0:.1f}s) per lambda band "
	 "(ideals, measured, model P(Poisson <= 2), lost with >= 3 occurrences): "
	 +"; ".join(f"[{lo:g},{hi:g}): {n}, {mf:.4f}, {md:.4f}, {h3}" for(lo,hi,n,mf,md,h3)in rec) + f" {'ok' if ok_rec else 'FAIL'}")
	t0=time.time()
	ok_counts=ball_counts_check()
	import retain16 as RT0
	RT0.load_table(load_stats())
	Pm=DS.P()
	dres=DS.evaluate(Pm,DS.PLAN16,DS.BOOT16)
	two_log=dres['two_log_units']
	bt=dres['boot']
	lz=sorted((v for v in dres['DET'].values()if 'lazy' in v),key=lambda v: -v['survcost_log2'])[:3]
	lz_rows= []
	for i,v in enumerate(lz):
		mc=DS.lazy_mc(v['mu'],v['nb_r'],v['psi'],seed=seed+i)
		lz_rows.append((round(v['mu'],1),v['nb_r'],round(v['psi'],3),round(v['lazy'],4),round(mc,4),mc<=v['lazy'] +0.01))
	ok_lazy=all(x[-1]for x in lz_rows)and bt['own_log2'] <bt['sub_log2']
	say(f"[7l] lazy cofactorisation (mu, r, psi, bound, simulated, ok): {lz_rows}; lazy boot: {math.log2((1 - bt['psi'] ** Pm.BOOT_RELS) / (1 - bt['psi'])):.3f} bits of relations vs log2 {Pm.BOOT_RELS} {'ok' if ok_lazy else 'FAIL'}")
	cap=M.DESCENT_CAP
	markov=2** (two_log-cap)
	log2_EN,over,coll,pmin=repeat_bound(dres);OVER=2.0**over;REP=2.0**coll
	p_desc=1-2*bt['phi'] -markov-OVER-REP
	p_main=1.0-2.0** (la['log2_N'] +2-128)
	success=p_main*p_desc
	rows=DS.report(dres)
	sv= [dres['DET'][k]for k in dres['DET']if dres['DET'][k]['how'] =='sieve']
	hmax=max(DS.CELLS[d['R']]for d in sv)
	ok_width=hmax<1<<103 and(hmax-1) >>40<1<<63 and all(16*max(d['nf'],d['ng']) <1<<16 for d in sv)
	mu_min=min(r['mu']for r in rows);psi_max=max(max(r['psi']for r in rows),bt['psi'])
	say(f"[7] descent ({time.time()-t0:.1f}s): ball counts {'ok' if ok_counts else 'FAIL'}; {len(rows)} nodes, min mu {math.floor(mu_min*10)/10:.1f}, max psi {math.ceil(psi_max*1000)/1000:.3f}; "
	 f"boot n1={bt['n1']} nA={bt['nA']:.2f} kA={bt['kA']} R={bt['R']}: 2^{bt['trials_log2']:.2f} trials at 2^{bt['trial_log2']:.2f} (pr 2^{bt['pr_log2']:.1f}), "
	 f"own 2^{bt['own_log2']:.2f} + subtree 2^{bt['sub_log2']:.2f} x retry {bt['retry']:.3f}; one log 2^{bt['one_log_log2']:.3f}; two-log work 2^{two_log:.3f} units; "
	 f"cap 2^{cap:g}: Markov <= 2^{math.ceil((two_log-cap)*1000)/1000:.3f}; P(boot fails) < {math.ceil(bt['phi']*1e7)/1e7:.2e}; node counter: E[N] <= 2^{log2_EN:.2f} special-q nodes (counts16), N_MAX 2^{M.NODE_CAP_LOG2:g}, P(counter abort) <= 2^{over:.2f}; special-q repeats among <= N_MAX nodes in the smallest bin (2^{pmin:.2f} ideals): P <= 2^{coll:.2f}, charged; success = {math.floor(success*10000)/10000:.4f} (floor); record widths: half count <= "
	 f"2^{math.log2(hmax):.2f}, owner < 2^{math.log2((hmax - 1) >> 40):.2f}, max score {16 * max(max(d['nf'], d['ng']) for d in sv):.0f} {'ok' if ok_width else 'FAIL'}")
	nb,np_=nbar(8,12)
	say(f"[9] degree-1 ideals per prime, both sides, {np_} primes in [2^8, 2^12): {nb:.3f} (model NBAR = {M.NBAR})")
	ok9=nb<=M.NBAR*1.01
	la_n,la_br=check_la_kernel()
	tail,red9,tail_n,tail_ok=la_row_tail()
	say(f"[10] row tail (la_mac.s): {tail_n} instructions (loop-free, calls ok: {tail_ok}); {tail} per row <= charged C_ROW {M.C_ROW:.0f}; one 9-limb reduction {red9} <= charged C_RED9 {M.C_RED9:.0f}")
	def fst(k):
		x= [v[2] /2**k for v in M.FILT_RUNS[k]];n=len(x);m=sum(x) /n
		return m, (sum((y-m) **2 for y in x) / (n-1)) **.5/n**.5,n
	fs= {k:fst(k)for k in M.FILT_RUNS};K2=sorted(fs)[-2:];fneed=max(fs[k][0] +2*fs[k][1]for k in K2)
	nu_ok=K2==[18,19]and M.NU300>=fneed and all(fs[k][0] <=fs[15][0]for k in K2)and max(v[3]for vs in M.FILT_RUNS.values()for v in vs) <=255
	say(f"[10] N(300)/U by filtering scale (mean, SE, seeds): "+", ".join(f"2^{k} {m:.5f}/{s_:.5f}/{n}" for k, (m,s_,n)in sorted(fs.items()))
	 + f"; charged {M.NU300} >= max over 2^18, 2^19 (PREREG-G) of mean + 2 SE = {fneed:.6f}; means at 2^18 and 2^19 <= the 2^15 mean; "
	 f"largest merged coefficient bound {max(v[3] for vs in M.FILT_RUNS.values() for v in vs)} <= 255: {nu_ok}")
	say(f"[10] la_mac.s general group loop: {la_n} instructions, branches {la_br} {'ok' if la_n == 33 and la_br == ['b.ne'] else 'MISMATCH'}")
	bw_n,bw_br=check_la_kernel("bw_dot.s");bf_n,bf_br=check_la_kernel("bw_dot.s","LBB1_2:")
	say(f"[10] bw_dot.s bw_full (full-field multiply-accumulate) inner loop: {bf_n} instructions (charged {M.C_FULL:.0f}) {'ok' if bf_n == 97 and bf_n <= M.C_FULL and bf_br == ['b.ne'] else 'MISMATCH'}")
	wl="\n".join(_asm("wide12.txt"));ml="\n".join(_asm("measurements12x.txt"))
	fr=lambda k: {int(p):float(v)for p,v in re.findall(k+r" p=(\d+): Ir \d+ -> \d+, reps \d+ -> \d+, Ir per rep ([\d.]+)",wl)}
	ws,wr=fr("wsetup6w"),fr("wrec");nar=wr.get(2**63-25,0);big= [v for p,v in wr.items()if p>>63]
	pm=int(re.findall(r"Ir per pair min \d+ mean \d+ max (\d+)",wl)[0])
	recw= {int(p):F(int(b) -int(a),int(d) -int(c))for p,a,b,c,d in re.findall(r"TAB rec (\d+) (\d+) (\d+) (\d+) (\d+)",ml)}
	ok_wide= (len(ws) ==2 and max(ws.values()) <=406.87<=M.C_SETUP6*M.WIDE_SETUP and len(big) ==3 and 2.5*max(abs(v-nar)for v in big) <= (M.WIDE_REC-1) *nar
	 and wl.count("MATCH") ==6 and wl.count("certified 200") ==3 and "0 mismatches" in wl and "failures 0" in wl and "BYTE-IDENTICAL" in wl
	 and M.LOG2_WIDE_LAT>=math.log2(6*1.25506*2**64/ (64*math.log(2)) *2*pm/M.UNIT)
	 and[float(recw[p])for p in(257,65537,8388617,67108879,268435459)] == [x[1]for x in M.REC_TAB_W])
	say(f"[10w] primes in [2^63, 2^64) (wide12.c, wide12.txt, measurements12x.txt): setup6w {max(ws.values()):.2f} Ir <= 406.87 <= charged {M.C_SETUP6 * M.WIDE_SETUP:.2f}; "
	 f"record path per coset at three wide p within {max(abs(v - nar) for v in big):.2f} Ir of p = 2^63 - 25 ({nar:.2f}), 2.5x inside the {M.WIDE_REC - 1:.2f} factor; "
	 f"0bbca0b2 record path per hit {[round(float(recw[p]), 2) for p in sorted(recw)]} == REC_TAB_W; wide lattices <= 2^{math.log2(6 * 1.25506 * 2 ** 64 / (64 * math.log(2)) * 2 * pm / M.UNIT):.3f} <= 2^{M.LOG2_WIDE_LAT} units; "
	 f"checks and object identity in the logs {'ok' if ok_wide else 'FAIL'}")
	ol="\n".join(_asm("measurements-oh.txt"))
	br= [(int(e),int(c),int(sg),int(a.replace(',','')),float(pa),int(z.replace(',','')))for e,c,sg,a,pa,z in re.findall(r"entries=(\d+) rows=\d+ what=build cells=(\d+) K=\d+ L4=\d+ segments=(\d+).*\n *([\d,]+) \(([\d.]+)%\)  \?\?\?:build_layout'2.*\n *([\d,]+) \(",ol)]
	ok_build=len(br) ==4 and all(a/e<=M.C_BUILD_E and a/ (pa/100) -a-z<2**20 and(a/ (pa/100) -a-z) /sg<=M.C_BUILD_SEG for e,c,sg,a,pa,z in br)
	say(f"[10o] layout per entry "+", ".join(f"{a / e:.3f}" for e,c,sg,a,pa,z in br) + f" <= {M.C_BUILD_E}; run remainder < 2^20 "
	 f"and per segment <= 2^17: {'ok' if ok_build else 'FAIL'}")
	Lx=_asm("measurements12x.txt");q12= {};grp= {}
	for ln in Lx:
		f=ln.split()
		if ln.startswith("TAB line "):q12[int(f[2])] = (int(f[4]) -int(f[3])) / (int(f[6]) -int(f[5]))
		if ln.startswith("MP "):
			d=grp.setdefault((int(f[1]),int(f[2])), {}).setdefault(int(f[3]), [0,0,set(),0,0])
			d[0] +=int(f[6]) -int(f[5]);d[1] +=int(f[8]) -int(f[7]);a_,b_=map(int,f[4].split('-'));d[2] |=set(range(a_,b_+1));d[3] =int(f[10]);d[4] =int(f[11]) * (int(f[13]) -int(f[12]))
	ok_line=len(q12) ==4 and all(c>=q12[p]for(h,c),p in zip(M.LINE_TAB[1:], (257,4099,65537,1048583)));lmsg= []
	for(lp,pt),sds in sorted(grp.items()):
		hs= [math.log2(d[1] / (2*d[4] *d[3]))for d in sds.values()];cs= [d[0] /d[1]for d in sds.values()];n_=len(cs);mc=sum(cs) /n_
		ub=sum(d[0]for d in sds.values()) /sum(d[1]for d in sds.values()) +M.LINE_SE_K*math.sqrt(sum((c-mc) **2 for c in cs) / (n_-1) /n_)if n_>1 else float('inf')
		ok_line&=all(d[2] ==set(range(1,d[3] +1))for d in sds.values())and n_>=3 and ub<=M.tab(M.LINE_TAB,max(hs));lmsg.append(f"{lp}/{pt}: {ub:.2f} <= {M.tab(M.LINE_TAB, max(hs)):.2f}")
	sp= [ln.split()for ln in Lx if ln.startswith("SP ")]
	ok_line&=len(sp) >0 and all(M.tab(M.LINE_TAB,math.log2(int(f[8]) /int(f[9]) /int(f[2]))) >= (int(f[5]) -int(f[4])) / (int(f[7]) -int(f[6]))for f in sp)
	say(f"[10l] line path (measurements12x.txt, public slice12): base points {[round(q12[p], 3) for p in sorted(q12)]} <= LINE_TAB; random-pair groups (slice/p: pooled + 2 SE <= charged) "+", ".join(lmsg) + f"; {len(sp)} small-p runs below the table: {'ok' if ok_line else 'FAIL'}")
	okr,rmsg=kernel_rec7(pd['runs'],nunits)
	okc,cc=contam_check(cfg,pd);okr&=okc
	say(f"[10x] exact-canonical runs (PREREG-R): {cc['cand']} candidates, R {cc['R']:.4f} (SE {cc['se']:.4f}, charged {cc['r_up']}); {cc['draws']} Cc-weighted draws, {cc['points']} a-values, {cc['bad']} contaminated: eps {cc['eps']}, factor {cc['factor']}; {len(cc['replay'])} rows' m, Cc and {len(cc['checked'])} runs replayed {'ok' if okc else 'FAIL'}")
	say(f"[10r] rec7 run walk (rec7.c, rec7.txt, runs16.py): {rmsg} {'ok' if okr else 'FAIL'}")
	oh= [ln.split()for ln in ol.splitlines()if ln.startswith("OH ")]
	sc= [(int(f[2]), (int(f[3]) -int(f[4])) /int(f[5]))for f in oh if f[1] =="scan"]
	ap= [(int(f[2]) -int(f[3])) /int(f[4])for f in oh if f[1] =="apply"];rv= [int(f[2]) /int(f[3])for f in oh if f[1] =="recover"]
	need=max(v*b-M.C_SCAN_CELL*b for b,v in sc)
	ok_oh=set(range(1,65)) <= {b for b,v in sc}and need<=M.C_SCAN_LEAF and len(ap) ==len(rv) ==1 and ap[0] <=M.C_APPLY and rv[0] <=M.C_RECOVER and M.C_ZERO==0 and ol.count("MATCH | survivors") >=4
	say(f"[10s] scan2, apply, recover (slice12oh.c, measurements-oh.txt): {len(sc)} scan points incl. every leaf length 1..64: per leaf <= {M.C_SCAN_CELL} x length + {M.C_SCAN_LEAF} (largest need {need:.3f}); "
		f"apply {ap[0]:.3f} <= {M.C_APPLY}, recovery {rv[0]:.3f} <= {M.C_RECOVER} per record; zeroing fused; brute-force checks MATCH: {'ok' if ok_oh else 'FAIL'}")
	b2_n,b2_br=check_la_kernel("bw_dot2.s","LBB0_2:")
	k5_n,k5_br=check_la_kernel("krylov5.s",".LBB0_1:");g1_n,g1_br=check_la_kernel("krylov5.s",".LBB1_1:")
	k5_ok=k5_br== ['b.ne']and g1_br== ['b.ne']and M.C_ENTRY_K>=M.UNIT_FRAC_CHARGED* (4+3+k5_n) /5+ (1-M.UNIT_FRAC_CHARGED) * (4+3+g1_n)
	say(f"[10] krylov5.s (Krylov pass, five interleaved vectors): unit loop {k5_n} instructions per entry for five vectors, general loop {g1_n} per vector-entry; "
	 f"per vector-entry {M.UNIT_FRAC_CHARGED} x (4 + 3 + {k5_n})/5 + {1 - M.UNIT_FRAC_CHARGED:.2f} x (4 + 3 + {g1_n}) = {M.UNIT_FRAC_CHARGED * (4 + 3 + k5_n) / 5 + (1 - M.UNIT_FRAC_CHARGED) * (4 + 3 + g1_n):.3f} <= charged C_ENTRY_K {M.C_ENTRY_K} {'ok' if k5_ok else 'MISMATCH'}")
	say(f"[10] bw_dot2.s (128-bit X entry times field element, projections) inner loop: {b2_n} instructions (charged {M.C_PROJ:.0f}) {'ok' if b2_n <= M.C_PROJ and b2_br == ['b.ne'] else 'MISMATCH'}")
	t0=time.time()
	au=audit(npr=dres['DET'][((220,225),'g')]['pr_log2'],bpr=bt['pr_log2'])
	ok11=au['boot']['within_model_bins']and au['specialq_within_bins']
	say(f"[11] audit ({time.time()-t0:.1f}s): certified boot lattice ({au['boot']['swaps']} swaps), boot cg max {max(s['cg_max'] for s in au['boot']['stats'])} "
	 f"<= {max(x[1] for x in DS.P.BOOT_BINS):.1f} {'ok' if au['boot']['within_model_bins'] else 'FAIL'}; special-q lattices "
	 +"; ".join(f"{s['side']}@2^{s['bits']}: cf max {max(x['cf_max'] for x in s['stats'])}, cg max {max(x['cg_max'] for x in s['stats'])}" for s in au['specialq'])
	 + f" within bins (cf <= {max(b[1] for b in DS.P.QBINS):.1f}, cg <= {max(b[2] for b in DS.P.QBINS):.1f}) {'ok' if au['specialq_within_bins'] else 'FAIL'}; "
	 f"complete-orbit relw yield ratio E[s/m_true]/pr (live, uniform sampler; unit_orbit.py): node220g {au['relw']['node_yield_ratio']:.2f} (m_true max {au['relw']['node_m_max']}, unit height max {au['relw']['node_height_max']}), boot {au['relw']['boot_yield_ratio']:.3f} (m_true max {au['relw']['boot_m_max']}, height max {au['relw']['boot_height_max']}) vs charged UDUP {DS.P.UDUP}; "
	 f"sampler shell-chi2 {au['sampler']['chi2']:.1f} (< 60 {'ok' if au['sampler']['ok'] else 'FAIL'}); "
	 f"largest line <= 2^{au['storage']['line_cells_log2_max']} cells; chunk protocol {au['storage']['routing_stable']}")
	ok11r=au['relw']['node_yield_ratio'] >1.0 and au['relw']['boot_yield_ratio'] >0.70
	ok11s=au['sampler']['ok']
	import descent_witness as W
	import retain16 as RT
	import ecm16 as EC
	t0=time.time();st=load_stats()
	nb=st['nb_empty_window']
	rr= [(b,s_,R_) +W.retention(b,s_,R_,seed+b,40)for b,s_,R_ in[(b,s_,dres['DET'][(bb,s_)]['R'])for b,s_ in((69,'g'), (150,'g'))for bb in dres['bins']if bb[0] ==b]]
	nb_ok=all(z<=zp for b_,s_,R_,lb_,m_,r_,z,zp in nb['rows']if m_>=7.4)
	tab=st['prime_structure']['table']
	special= [p for p in RT.primes(2,512)if any(e!=1 or any(m!=1 for m in mf+mg)for d1,e,mf,mg in tab[str(p)])]
	spot=sorted(set([p for i,p in enumerate(RT.primes(2,512))if i%7==seed%7] +special))
	ok_tab=len(tab) ==len(RT.primes(2,512))and all(RT.structure(p) ==tab[str(p)]for p in spot)
	sample= {}
	for p in(2,23):
		rs=W.residue_valuations(p,400,seed+p)
		for side in 'fg':
			bound=1-RT.local_pdf(p,tab[str(p)],side,False,True).get(0,0.0)
			frac=1-rs[side].get(0,0) /rs['n']
			sample[f"{p}{side}"] = (round(frac,4),round(bound,4),frac<=bound+3*math.sqrt(bound* (1-bound) /rs['n']) +1e-12)
	ok_sample=all(v[2]for v in sample.values())
	pf,ordf=RT.omitted_pdf('f',True);pg,ordg=RT.omitted_pdf('g',True)
	ex=lambda pdf:sum(x*v for x,v in pdf.items()) /RT.GRID
	ret_min=min(v['ret']for v in dres['DET'].values()if v['how'] =='sieve');loss_max=max(max(v['loss_f'],v['loss_g'])for v in dres['DET'].values())
	ok11b=nb_ok and all(x[3]for x in rr)and ok_tab and ok_sample
	say(f"[11b] retention ({time.time()-t0:.1f}s): prime-ideal table (degree, exponent, root multiplicities of f, g) recomputed at {len(spot)} primes {'ok' if ok_tab else 'FAIL'}; "
	 "repeated roots or ramification at "+"; ".join(f"{p}: {tab[str(p)]}" for p in special)
	 +"; residue-class samples (400 base-primitive classes mod p^K, P(p | N) sample/bound): "+", ".join(f"{k} {v[0]}/{v[1]}" for k,v in sample.items()) + f" {'ok' if ok_sample else 'FAIL'}"
	 + f"; base-primitive {RT.base_primitive():.4f}, ordinary {ordf:.4f}/{ordg:.4f}; omitted mass means {ex(pf):.2f}/{ex(pg):.2f} bits, "
	 f"P(X > 96) {RT.cell_loss('f', 96):.2e}/{RT.cell_loss('g', 96):.2e}; loss <= {loss_max:.1e}, retention >= {math.floor(ret_min * 1e4) / 1e4:.4f}; bounds >= |N|: "+"; ".join(f"{b}{s_} R={R_}: {ok_}, max gap {g:.1e} bits" for b,s_,R_,ok_,g,k in rr)
	 +"; empty-window toys (bits, side, R, lb, mean, r, zeros, zeros at r=0.5): "+"; ".join(str(tuple(x))for x in nb['rows'])
	 + f" {'ok' if nb_ok else 'FAIL'}")
	t0=time.time()
	toyecm=EC.toy(60,seed,sizes=((36,1024,50), (40,2048,50)))
	sch= [EC.schedule(69.0,int(math.ceil(dres['DET'][((69,75),'f')]['nf'])) ),EC.schedule(bt['nA'],int(math.ceil(bt['m'])))]
	ok12=all(t['empirical'] >=t['model_corrected'] -2*math.sqrt(t['model_corrected'] * (1-t['model_corrected']) /t['trials'])for t in toyecm)
	say(f"[12] ECM schedules ({time.time()-t0:.1f}s): "+"; ".join(
	 f"{x['fbits']:.2f}-bit primes in {x['nbits']}-bit integers: B1 2^{math.log2(x['B1']):.0f}, B2 {x['B2'] // x['B1']} B1, wheel {x['wheel']}, "
	 f"{x['curves']} curves at p = {x['p_curve']:.2e}, {x['stage1_M']:.3g} + {x['stage2_M']:.3g} M per curve, 2^{x['log2_instructions']:.2f} instructions" for x in sch)
	 +"; toy detection (bits, trials, stage 1, stage 2, empirical, P1 without/with the power correction): "
	 +"; ".join(f"({t['bits']}, {t['trials']}, {t['stage1']}, {t['stage2']}, {t['empirical']:.3f}, {t['model']:.3f}, {t['model_corrected']:.3f})" for t in toyecm)
	 + f" {'ok' if ok12 else 'FAIL'}")
	t0=time.time();fl=check_fields()
	say(f"[13] number fields ({time.time()-t0:.1f}s): N(f(X)) of degree {fl['deg_f']} and N(g(X)) of degree {fl['deg_g']} squarefree mod r {'ok' if fl['unramified'] else 'FAIL'} (r unramified in K_f, K_g); "
	 f"disc h = {fl['disc_h']}; signatures K_f {fl['sig']['f']}, K_f+ {fl['sig']['fp']}, K_g {fl['sig']['g']}: odd unit rank {fl['sig']['odd']} = SM columns {M.N_SM} {'ok' if fl['sig']['odd'] == M.N_SM else 'FAIL'}")
	ok13=fl['unramified']and fl['sig']['odd'] ==M.N_SM
	import time as _t
	_t0=_t.time();al_s= {p:alpha_contrib(p)for p in(2,19,23)}
	ok_alpha=op['alpha_f'] >=ALPHA_F and op['alpha_f'] >=ALPHA_F_50K and op['alpha_g'] >=ALPHA_G and op['alpha_g'] >=ALPHA_G_50K and DS.P.ALPHA_G>=ALPHA_G and DS.P.ALPHA_G>=ALPHA_G_50K and DS.P.ALPHA_F>=ALPHA_F and DS.P.ALPHA_F>=ALPHA_F_50K and all(round(v[0],3) >=0 for v in al_s.values())
	say(f"[14] Murphy alpha (p <= {ALPHA_PMAX} both sides, p < 50000 both sides, exact Kummer-Dedekind roots): computed ({ALPHA_F:.6f}, {ALPHA_G:.6f}; 50k ({ALPHA_F_50K:.6f}, {ALPHA_G_50K:.6f})) <= charged ({op['alpha_f']}, {op['alpha_g']}) {'ok' if ok_alpha else 'FAIL'}; "
	 f"in-run exact checks at p=2/19/23 (f,g): "+", ".join(f"{p}: ({v[0]:.6f}, {v[1]:.6f})" for p,v in sorted(al_s.items())) + f" ({_t.time()-_t0:.1f}s)")
	tot=M.total_log2(e['log2_relcol'],la['log2_la'],cap,la['log2_gen'])
	say(f"ledger (log2 units): relation collection {e['log2_relcol']:.4f}, linear algebra {la['log2_la']:.4f}, generator {la['log2_gen']:.2f}, "
	 f"descent cap {cap:g} (expected {two_log:.3f}), canonical list {M.LOG2_CANON:.0f}, Schirokauer maps {M.LOG2_SM:.0f}, "
	 f"factor base {M.LOG2_FBC:.0f}, misc {M.LOG2_MISC:.0f}, verify {M.LOG2_VERIFY:.0f}")
	say(f"total_log2W_full_ledger = {tot:.12f} (claimed {CLAIM:.10f}: {'ok' if tot <= CLAIM else 'EXCEEDS CLAIM'}); "
	 f"reserve {CLAIM - tot:.3f} bits")
	vals=dict(ledger=tot,relcol=e['log2_relcol'],la=la['log2_la'],gen=la['log2_gen'],log2_rels=e['log2_rels'],rels_1se=rels_lo,margin=e['margin'],
	   work_se=work_se,two_log=two_log,markov=markov,log2_EN=log2_EN,over=math.log2(OVER),rep=math.log2(REP),
	   boot_phi=bt['phi'],success=success)
	dec_ok= {k: (vals[k] <=v if op=='<=' else vals[k] >=v)for k, (op,v)in DECLARED.items()}
	dec_ok['work_factor'] =782**12*10**9<=780**12*1031206887 and float(corr['safe_work_population_factor']) <=1.031206887
	ok_declared=all(dec_ok.values())
	say("declared bounds (unrounded value vs stated bound): "+"; ".join(f"{k} {vals[k]!r} {op} {v} {'ok' if dec_ok[k] else 'FAIL'}" for k, (op,v)in DECLARED.items()))
	checks=dict(owner_schedule=O.check(M,la),resultant=ok0,tower=ok0b,rho_semismooth=ok1,pool=ok2 and ok_size,box_bounds=ok2b,nbar=ok9,la_loop=k5_ok and la_n==33 and bf_n==97 and bf_br== ['b.ne']and b2_n<=M.C_PROJ and b2_br== ['b.ne']and check_la_kernel('la_mac.s','LBB9_1:') == (11, ['b.ne'])and tail_ok and tail<=M.C_ROW and red9<=M.C_RED9 and nu_ok and M.C_ENTRY>=4+3+M.UNIT_FRAC_CHARGED*11+ (1-M.UNIT_FRAC_CHARGED) *33 and ok_wide and ok_build and ok_line and ok_oh,rec7=okr,
	    ball_counts=ok_counts and ok_relmeans,lazy=ok_lazy,rels_input=ok_rels_input,work_allowance=ok_allow,declared_bounds=ok_declared,alpha_ok=ok_alpha,missing_logs=ok_rec,record_widths=ok_width,descent_plan_feasible=mu_min>=Pm.MIN_RELS and psi_max<=Pm.PSI_MAX and max(r['phi']for r in rows) <=0.05 and coll<= -4,
	    audit=ok11,retention=ok11b,ecm=ok12,fields=ok13,feasible=e['feasible'],within_claim=tot<=CLAIM,success_ok=success>=0.39,udup_relw=ok11r,sampler_uniform=ok11s)
	say(f"done in {time.time()-T0:.1f}s; checks: {checks}")
	print(json.dumps({"checks_ok":all(checks.values()),"checks":checks,
	     "total_log2W_full_ledger":math.ceil(tot*10**12) /10**12,"claimed":CLAIM,
	     "log2_relcol_ceil4":math.ceil(e['log2_relcol'] *10000) /10000,"log2_rels_floor4":math.floor(e['log2_rels'] *10000) /10000,
	     "margin_bits_floor4":math.floor(e['margin'] *10000) /10000,"work_se_ceil5":math.ceil(work_se*100000) /100000,"work_allowance":M.CANON_WORK_ALLOWANCE,"yield_se_log2":round(e['yield_se_log2'],4),"log2_N":round(la['log2_N'],4),
	     "log2_la_ceil4":math.ceil(la['log2_la'] *10000) /10000,"descent_two_logs_log2_units_ceil4":math.ceil(two_log*10000) /10000,
	     "descent_cap":cap,"markov_abort":markov,"log2_expected_nodes_ceil4":math.ceil(log2_EN*10000) /10000,"node_cap_log2":M.NODE_CAP_LOG2,"counter_abort":OVER,"repeat":REP,"success_exact":success,"modelled_success_floor4":math.floor(success*10000) /10000,
	     "descent_nodes":len(rows),"descent_boot_log2":round(bt['one_log_log2'],3)}),flush=True)
if __name__=="__main__":
	if "--alpha" in sys.argv:
		print(json.dumps(dict(zip(("alpha_f","alpha_g"),alpha_sweep(ALPHA_PMAX))),sort_keys=True));sys.exit(0)
	if "--alpha50k" in sys.argv:
		print(json.dumps(dict(zip(("alpha_f","alpha_g"),alpha_sweep(ALPHA_G_PMAX))),sort_keys=True));sys.exit(0)
	main()

""
import base64,hashlib,json,lzma,math,os,re,sys,time
from decimal import Decimal,localcontext
from fractions import Fraction as F
from math import log2,ceil,floor,sqrt
K=lambda b:'ok' if b else 'FAIL'
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
CLAIM=115.64681
import retain16
ALPHA_PMAX=2000
ALPHA_F,ALPHA_G=-0.047244,0.706684
ALPHA_G_PMAX=50000
ALPHA_G_50K=0.745775
ALPHA_F_50K=0.285135
def alpha_contrib(p):
	if p==2:
		import runs16 as local2
		return local2.alpha2()
	vf=vg=0.0
	for d1,e,mf,mg in retain16.structure(p):
		q=float(p**d1)
		for side,mults in enumerate((mf,mg)):
			v=0.0
			for mu in mults:
				v+=q/(q*q-1.0)if mu==1 else 1.0/(q+1.0)
			if side==0:vf+=d1*v
			else:vg+=d1*v
	t=math.log(p)/math.log(2.0)
	return t*(1.0/(p-1.0)-vf),t*(1.0/(p-1.0)-vg)
def toy_v2(T,k=7):
	import itertools as I
	V=[v for v in I.product((-1,0,1),repeat=6)if v>tuple(-x for x in v)];from threshold_model import primes_upto as pu;P=pu(1<<20);r=[0,0,0];ok=1;o=""
	while P[1:]:P=[math.prod(P[i:i+2])for i in range(0,len(P),2)]
	for j,(a,b)in enumerate(I.product(V,V)):
		if j%k==0 and M.coprime_pair(a,b):z=abs(C.norms(a,b)[0]);r[0]+=1;r[1]+=z%64==0;r[2]+=pow(P[0],256,z)==0
	for lb,x in sorted(T['lbs'].items()):
		c=x['cls'].values();t=[sum(d[j]for d in c)for j in(1,2,3)];q=lambda j:sum((d[1]-d[j])**2/d[j]for d in c);ok&=t[2]<=t[1]+t[0]**.5 and q(3)<q(2);o+=f" 2^{lb}: {[d[1:]for d in c]} chi2 {q(2):.0f}/{q(3):.0f};"
	say(f"[14t] toy census ({T['n']} pairs; v2 0/6: obs, avg, exact):{o} replay {r} {K(ok and r==T['replay'])}")
	return ok and r==T['replay']
def alpha_sweep(pmax):
	af=ag=0.0
	for p in retain16.primes(2,pmax):
		a,b=alpha_contrib(p);af+=a;ag+=b
	return af,ag
DECLARED={'ledger':('<=',115.64681),'relcol':('<=',115.300026350),'la':('<=',113.232777467),'gen':('<=',86.048),'log2_rels':('>=',62.5959),'rels_2se':('>=',62.5238),'margin':('>=',0.1335),'work_se':('<=',0.005049),'two_log':('<=',109.5279),'markov':('<=',0.553406),'log2_EN':('<=',15.2673),'over':('<=',-12.03),'rep':('<=',-10.7),'boot_phi':('<=',5.10e-06),'success':('>=',0.4456),'success_head':('>=',0.3903)}
_REPORT=[]
T0=time.time()
load_stats=C.stats
def say(*args):
	line=" ".join(str(a)for a in args)
	_REPORT.append(line)
	print(line,file=sys.stderr,flush=True)
def f6_mulp(A,B):
	return[c%P for c in C.pmul([c%P for c in A],[c%P for c in B])]
def f6_pow(A,n):
	r,base=[1],[c%P for c in A]
	while n:
		if n&1:
			r=f6_mulp(r,base)
		base=f6_mulp(base,base)
		n>>=1
	return(r+[0]*6)[:6]
def check_resultant():
	F0=[0]*6;uw=[-U*x for x in C.WP]
	for d,c in C.PTERMS:
		term=(c+[0]*6)[:6]
		for _ in range(d//2):
			term=C.pmul(term,uw)
		F0=[x+y for x,y in zip(F0,term)]
	r=C.norm(C.elt(F0))
	ok=r!=0 and r%P==0
	say(f"[0] Res: {r.bit_length()} bits, 0 mod p {K(ok)}")
	return ok
def dedekind(p):
	import retain16 as RT
	hp=[x%p for x in H];rep=[]
	for g in RT.factor_Fp(hp,p):
		if g not in rep:rep.append(g)
	prod_=[1];sq=[]
	for g in rep:
		e,rest=0,hp
		while True:
			q,rm=RT.pdiv(rest,g,p)
			if rm:break
			rest,e=q,e+1
		for _ in range(e):
			prod_=[sum(prod_[i]*g[k-i]for i in range(len(prod_))if 0<=k-i<len(g))for k in range(len(prod_)+len(g)-1)]
		if e>1:sq.append(g)
	assert len(prod_)==7 and all((a-b)%p==0 for a,b in zip(H,prod_))
	Fd=RT._trim([(a-b)//p%p for a,b in zip(H,prod_)])
	return bool(sq)and all(len(RT.pgcd(Fd,g,p))==1 for g in sq)
def check_tower():
	yp=f6_pow([0,1],P)
	cur,ok_irr=[0,1,0,0,0,0],True
	for k in range(1,7):
		acc=[0]*6
		for coef in reversed(cur):
			acc=f6_mulp(acc,yp)
			acc[0]=(acc[0]+coef)%P
		cur=acc
		is_y=cur==[0,1,0,0,0,0]
		if k==3:
			y3=cur
		if(k in(1,2,3)and is_y)or(k==6 and not is_y):
			ok_irr=False
	uw=[(-U*x)%P for x in C.WP]
	ok_sq=f6_pow(uw,(P**6-1)//2)==[P-1,0,0,0,0,0]
	ok_tau=y3==[x%P for x in C.SIGP]and f6_pow(uw,(P**3-1)//2)==[x*pow(2,-1,P)%P for x in C.C2P]
	import descent_witness as W
	chi=pow(P,3,W.R);ok_chi=(chi*chi+1)%W.R==0
	tf=lambda t:[C.sigma(c)for c in t]
	hs=[0]*6;pw=[1,0,0,0,0,0]
	for k in range(7):
		hs=[x+H[k]*y for x,y in zip(hs,pw)];pw=C.pmul(pw,C.SIGP)
	cc=C.CU;c2=C.mul(cc,cc);ok_aut=not any(hs)and C.mul(cc,C.sigma(cc))==[-1,0,0,0,0,0]and all(C.sigma(C.sigma(e))==e for e in([int(i==j)for j in range(6)]for i in range(6)))
	ok_aut&=all(C.mul(C.sigma(c),C._pw(cc,d))==C.mul(C._pw(c2,6),c)for d,c in C.TERMS)and C.mul(C.sigma(C.UW),[1,0,0,0,0,0])==C.mul(c2,C.UW)
	import runs16 as local2
	l2,_=local2.facts()
	d2=[[F(x,2)for x in r]for r in C.BAS2];det=F(IC.determinant([[2*x for x in r]for r in d2]),2**6)
	disc=C.norm(C.elt([0,-4,-6,8,0,6]))
	ok_c=disc==2957312==16*184832 and abs(det)==F(1,4)and l2 and dedekind(19)
	say(f"[0b] h irred {K(ok_irr)}; (u w)^((p^6-1)/2) = -1 {K(ok_sq)}; tau = Frob^3 {K(ok_tau)}; "
	f"chi^2 = -1 {K(ok_chi)}; sigma, c, tau(f), tau(g) {K(ok_aut)}")
	say(f"[0c] |Disc(h)| = |N(h'(y))| = {disc} = 4^2 x 184832; index 4; 2-maximal {K(l2)}, 19-maximal {K(dedekind(19))}: Disc(O_h) = 2^9 x 19^2 {K(ok_c)}")
	return ok_irr and ok_sq and ok_c and ok_tau and ok_chi and ok_aut
def check_rho_and_semismooth():
	ref={2:0.30685281944,3:0.048608388,4:0.0049109256,6:1.9649696e-5,8:3.2320693e-8,10:2.7701718e-11}
	err=max(abs(DS.rho_series(u)-v)/v for u,v in ref.items())
	ok_br=all(DS.rho_lo(u)<=DS.rho_series(u)<=DS.rho_hi(u)+1e-300 for u in(1.5,2.3,7.7,15.2,33.3,41.0))
	rows,h=DS.window_rows(80.0,69.0,48,48)
	ok_w=abs(sum(rows[(1,0)])-math.log(80/69))<1e-12
	say(f"[1] rho vs reference {err:.1e} {K(err<1e-6)}; brackets {K(ok_br)}; weights "
	f"{K(ok_w)}")
	return err<1e-6 and ok_br and ok_w
def canonical(cfg,a):
	try:
		return IC.canonical_a(a,cfg['matrix_a'],100000)[0]
	except RuntimeError:
		return IC.canonical_a(a,cfg['matrix_a'],10**7)[0]
def verify_pool(cfg,pd):
	bad=[];rows=[]
	import runs16 as R16
	units=R16.unit_list();qi,dq,_=C.integer_metric(cfg['matrix_a']);lim=int(F(pd['theta_region'])*cfg['production_A_squared']*dq)
	for an in pd['anchors']:
		sd,i=an[0],an[1];r=M.pool_row(cfg,sd,i);a=C.regenerate_pair(cfg,sd,i)[0]
		f=[round(r['qa']*65536),round(r['nf']*16),round(r['ng']*16),r['primitive']+2*r['rejected'],r['evals'],r['bcells'],r['leaves'],r['blen'],r['rlen']//16]+list(R16.run(a,units,qi,lim))
		lo,hi=f[-2:];r['lrs']=R16.lrs(a,lo,hi,qi,dq,cfg['production_A_squared'],F(cfg['production_region_lam']),F(cfg['production_region_kap']));r['T']=hi-lo+1;r['run']=(a,lo,hi)
		if f!=an[2:13]or not canonical(cfg,a):
			bad.append(('anchor',sd,i))
		rows.append(r)
	n_can=0
	for sd,(hi,idx)in sorted(pd['replay'].items()):
		got=[i for i in range(hi)if canonical(cfg,C.regenerate_pair(cfg,int(sd),i)[0])]
		n_can+=len(got);bad+=[('replay',sd)]if got!=idx else[]
	say(f"[2] {len(pd['anchors'])} anchors (all fields, run ends), draws 0..{hi-1} of seeds {', '.join(sorted(pd['replay']))} "
	f"({n_can} canonical) lists: {'ok' if not bad else 'FAIL '+str(bad[:3])}")
	return not bad,n_can,rows,len(units)
def box_check(cfg,pool,k=6,cells=8,seed=5):
	import embedding_bounds as EB
	rng=Random(f"box-check:{seed}")
	qa,da,_=C.integer_metric(cfg['matrix_a']);qb,db,_=C.integer_metric(cfg['matrix_b'])
	A2=cfg['production_A_squared'];disks=EB.root_disks();n=bad=nrow=0
	for r in rng.sample(pool,k):
		a,b,aq,bq=C.regenerate_pair(cfg,r['seed'],r['index'])
		est=EB.NormBounds(a,disks)
		box=M.box_rows(qb,(F(cfg['production_region_lam'])*A2-F(cfg['production_region_kap'])*F(aq,da))*db,b,(b[0]//64)*64);nrow+=len(box)
		leaves=EB.box_leaves(est,[(est.row(list(s)),l,h)for s,l,h in box],min(x[1]for x in box),max(x[2]for x in box))
		for kk,(s,l,h)in enumerate(box):
			for b0 in rng.sample(range(l,h+1),min(cells,h-l+1)):
				lf=[x for x in leaves if x[0]<=b0<=x[1]and x[3]in(None,kk)][0]
				if isinstance(lf[2],dict):
					continue
				for(lo16,hi16),v in zip(lf[2][:1],C.norms(a,[b0]+list(s))[:1]):
					n+=1;bad+=not(v!=0 and lo16-1e-6<=16*log2(abs(v))<=hi16+1e-6)
	return bad==0,n,nrow
def poly_roots_count(poly,p):
	import retain16 as RT
	g=RT._trim([c%p for c in poly])
	if len(g)<=1:return 0
	return len(RT.pgcd(g,RT.padd(RT.ppowmod([0,1],p,g,p),[0,p-1],p),p))-1
def nbar(lo_bits,hi_bits):
	from threshold_model import primes_upto
	ps=[q for q in primes_upto(2**hi_bits)if q>=2**lo_bits]
	tot=tf=0
	rng=Random(1)
	for p in ps:
		for s in sorted(set(roots_modp(H,p,rng))):
			w=s+1;c=poly_roots_count([w**6,0,w**5,0,0,0,2*w**3,0,0,0,-2*w,0,1],p);tf+=c;tot+=c+poly_roots_count([U*w,0,1],p)
	return tot/len(ps),len(ps),tf/len(ps)
def _asm(name):
	t=lzma.decompress(base64.b85decode("".join(open(os.path.join(_DIR,"kernels_asm.b85")).read().split()))).decode()
	parts=re.split(r"(?m)^@@ ([A-Za-z][A-Za-z0-9_.-]*)\n",t)
	return dict(zip(parts[1::2],parts[2::2]))[name].rstrip("\n").splitlines()
def kernel_rec7(ru,nunits):
	L=_asm("rec7.txt");T=M.REC7_TAB
	P_=lambda tg:[(log2(int(f[9])/int(f[10])),int(f[8])/int(f[9]),int(f[5]))for f in(ln.split()for ln in L if ln.startswith(tg))]
	pn,pw,ph=P_("CAMP "),P_("CAMPW "),P_("CAMPH ");pts=pn+ph+[x for x in pw if x[0]>=T[-1][0]]
	ok=len(pn)==130 and len(pw)==42 and len(ph)==96 and all(x[2]>>63 for x in pw)and all(a[0]>b[0]and a[1]<=b[1]for a,b in zip(T,T[1:]))and all(lh>=T[-1][0]and c<=M.step7(lh)for lh,c,_ in pts)and min(log2(x[2])for x in ph)<=M.P7
	chk=[ln for ln in L if ln.startswith(("CHECK ","WCHECK ","TCHECK "))];ok&=len(chk)==3 and "ok (0 mismatches over 270 cosets" in chk[0]and chk[1].endswith("mismatches 0")and chk[2].endswith("invalid 0: ok")
	pm=max(int(ln.split()[4])for ln in L if ln.startswith("PREP "));ideals=M.NBAR*(2**M.log2_li(64)-2**M.log2_li(M.P7))
	ok&=M.LOG2_PREP7>=log2(6*ideals*pm*2**0.05/M.UNIT)and M.P7>=48 and M.EXTRA7>=2
	cp={(f[1],f[5],f[7]):int(f[10])for f in(ln.split()for ln in L if ln.startswith("CAMP "))};cs=[ln.split()for ln in L if ln.startswith("COS ")]
	pc=[(log2(int(f[5+2*j])/cp[tuple(f[1:4])]),int(f[4+2*j])/int(f[5+2*j]))for f in cs for j in(0,1)];pc=[x for x in pc if x[0]>=T[-1][0]]
	cl=[(log2(int(f[4])/(2*int(f[7])*int(f[5])*int(f[6]))),int(f[3])/int(f[4]))for f in(ln.split()for ln in L if ln.startswith("COSL "))]
	E6=[[int(i==j)for j in range(6)]for i in range(6)]
	ok&=len(cs)==122 and all(c<=M.step7(lh)for lh,c in pc)and len(cl)==6 and all(lh>=M.LINE_TAB[-1][0]and c<=M.tab(M.LINE_TAB,lh)for lh,c in cl)and C.P2F[0]==1 and C.norm(C.W)==2 and all(sum(x*y for x,y in zip(C.mul(C.W,e),C.P2F))%2==0 for e in E6)
	Ts=[int(ln.split()[3])for ln in L if ln.startswith("CAMP ")];ok&=min(Ts)<=M.RUN_T[0]and M.RUN_T[1]<=max(Ts)
	ok&=nunits==ru['units']
	rn={int(p):F(int(b)-int(a),int(d)-int(c))for p,a,b,c,d in re.findall(r"record p=(\d+): Ir (\d+) -> (\d+), hits (\d+) -> (\d+)","\n".join(_asm("measurements12.txt")))}
	ok&=[float(rn[p])for p in sorted(rn)]==[x[1]for x in M.REC_TAB]
	return ok,(f"runs of {M.RUN_T[0]}..{M.RUN_T[1]} a-values in {min(Ts)}..{max(Ts)}; {len(pn)} + {len(ph)} (p >= 2^{min(log2(x[2])for x in ph):.2f}) + {len(pts)-len(pn)-len(ph)} wide points <= step table (worst {max(c/M.step7(lh)for lh,c,_ in pts):.4f}); {chk[0][6:]}; {chk[1][7:]}; {chk[2][7:]}; "
	f"prepare7 <= {pm} per ideal, 6 bases x 2^{log2(ideals):.2f} ideals: 2^{M.LOG2_PREP7}; coset: {len(cs)} COS ({len(pc)} sub-runs) <= step table, {len(cl)} COSL <= LINE_TAB, wO = ker P2F; runs: {nunits} units; REC_TAB == log {[round(float(rn[p]),3)for p in sorted(rn)]}")
def check_la_kernel(name,start):
	on,body=False,[]
	for ln in _asm(name):
		s=ln.strip()
		if s.startswith(start):
			on=True;continue
		if on:
			if s.startswith(("LBB",".LBB"))or(s.startswith(";")and "bb." in s):
				break
			if s and not s.startswith(";")and not s.startswith("."):
				body.append(s.split()[0])
	return len(body),[m for m in body if m.startswith("b.")or m=="b"]
def check_la_func(fn):
	L=[x.strip()for x in _asm("la_mac.s")]
	i=L.index(fn+":")+1;n,seen,fwd,calls=0,set(),True,[]
	while i<len(L)and not L[i].startswith(("_",";")):
		x=L[i];i+=1
		if x.startswith("LBB"):seen.add(x[:-1])
		elif x:
			n+=1;fwd=fwd and x.split()[-1]not in seen
			if x.split()[0]=="bl":calls.append(x.split()[-1])
	return n,fwd,calls
def la_part():
	c={};f="";ok=True
	for x in _asm("krylov5.s"):
		o=x.split()or[";"];g=f.split("/")[0]
		if x[-1:]==":":f=g+"/"*(x[-2]=="1")if x[0]=="."else x[:-1]
		elif o[0]!=";":
			c[f]=c.get(f,0)+1
			if f=="owner_part":ok&=o[0]in("stp","ldp","mov","add")or o[-1]in("unit_run5","gen_run1")
	return c["owner_part"]+2*c["unit_run5"]+10*c["gen_run1"],ok
def la_row_tail():
	f={k:check_la_func(k)for k in("_send_pair","_combine","_row_redc","_chi_mul")}
	d,dok=la_part();ok=all(v[1]and not v[2]for v in f.values())and dok
	f["owner_part"]=(d,)
	return O.row_tail(f["_send_pair"][0],f["_combine"][0],f["_row_redc"][0],d),f["_row_redc"][0]+6,f["_chi_mul"][0]+4,{k:v[0]for k,v in f.items()},ok
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
		ok&=(n12-1)//2==DS.CELLS[R]
	V12=math.pi**6/720
	for R,half in DS.CELLS.items():
		n=2*half+1
		ok&=V12*(R-sqrt(3))**12<=n<=V12*(R+sqrt(3))**12
	return ok
def roots_modp(poly,p,rng):
	import retain16 as RT
	return RT.roots_Fp(poly,p)if p<4096 else sorted((-f[0])%p for f in RT.factor_Fp(poly,p)if len(f)==2)
def repeat_bound(dres):
	import counts16
	EN,_=counts16.expected_nodes(DS.P(),dres)
	pmin=min(log2(2*(2**DS.log2li(b[1])-2**DS.log2li(b[0])))for b in dres['bins'])
	over=ceil((log2(EN)-M.NODE_CAP_LOG2)*100)/100
	coll=ceil((2*M.NODE_CAP_LOG2-1+6-pmin)*100)/100
	return log2(EN),over,coll,pmin
def audit(rng_seed=31001,npr=None,bpr=None):
	import descent_witness as W
	assert(W.P**6-1)%W.R!=0 and((W.P**12-1)//W.R)%W.R!=0
	rng=Random(rng_seed)
	while True:
		x=([rng.randrange(W.P)for _ in range(6)],[rng.randrange(W.P)for _ in range(6)])
		if any(x[0]+x[1]):
			break
	w0=W.pow12(x,(W.P**12-1)//W.R)
	tt=rng.randrange(1,W.R)
	z=W.pow12(w0,tt)
	B,swaps=W.reduce(W.boot_basis(z),W.QB)
	bst=W.boot_ball_stats(z,B,rng_seed+1,points=6,radii=(96,))
	boot_bins_max=max(x[1]for x in DS.P.BOOT_BINS)
	sq=[]
	Pm=DS.P()
	cfmax=max(b[1]for b in Pm.QBINS);cgmax=max(b[2]for b in Pm.QBINS)
	env_ok=True
	for bits,side,seed in((260,'f',13300),):
		r=W.specialq_ball_stats(bits,side,seed,points=6,radii=(96,))
		for st_ in r['stats']:
			env_ok&=st_['cf_max']<=cfmax and st_['cg_max']<=cgmax
		sq.append(dict(bits=bits,side=side,seed=seed,swaps=r['swaps'],stats=[{k:(round(v,2)if isinstance(v,float)else v)for k,v in s.items()}for s in r['stats']]))
	q_,_,_,_,nB,_=W.specialq_lattice(225,'g',Random(rng_seed+31))
	n_,iA_,kA_,kB_,R_,_=DS.PLAN16[(220,'g')]
	tr=W.orbit_triples(nB,rng_seed+32,R_,12)
	lq_=log2(q_);sn=0.0
	for lnf,lng,m,_h in tr:
		pf=DS.ss_lower2(lnf,Pm.LX,float(n_),iA_,kA_,kB_,Pm.MESH,Pm.ALPHA_F,Pm.TRUNC)
		pg=DS.ss_lower2(lng-lq_,Pm.LX,float(n_),iA_,kA_,kB_,Pm.MESH,Pm.ALPHA_G,Pm.TRUNC)
		sn+=pf*pg/m
	n1,iAb,kAb,_=DS.BOOT16
	tb=W.orbit_triples(B,rng_seed+33,4,24,boot=True)
	sb=sum(DS.ss_lower2(lng,Pm.LX,float(n1),iAb,kAb,1,Pm.MESH,Pm.ALPHA_G,Pm.TRUNC)/m for _,lng,m,_h in tb)
	orb=dict(node_m_max=max(t[2]for t in tr),node_height_max=max(t[3]for t in tr),
	boot_m_max=max(t[2]for t in tb),boot_height_max=max(t[3]for t in tb))
	srng=Random(rng_seed+77)
	dp=[1]+[0]*16
	for _ in range(12):
		ndp=[0]*17
		for s in range(17):
			for v in range(-4,5):
				if s+v*v<=16:
					ndp[s+v*v]+=dp[s]
		dp=ndp
	tot=sum(dp[1:]);sh=[0]*17;okshape=True
	for _ in range(600):
		c=W.int_ball_point(srng,4)
		n2=sum(x*x for x in c)
		okshape&=0<n2<=16
		for x in reversed(c):
			if x:
				okshape&=x>0
				break
		sh[n2]+=1
	chi2=sum((sh[s]-600*dp[s]/tot)**2/(600*dp[s]/tot)for s in range(1,17))
	return dict(boot=dict(swaps=swaps,stats=[{k:(round(v,2)if isinstance(v,float)else v)for k,v in s.items()}for s in bst],
	within_model_bins=all(s['cg_max']<=boot_bins_max for s in bst)),
	specialq=sq,specialq_within_bins=env_ok,
	relw=dict(node_yield_ratio=(sn/12)/2.0**npr,boot_yield_ratio=(sb/24)/2.0**bpr,**orb),
	sampler=dict(chi2=chi2,ok=bool(okshape)and chi2<60),
	storage=store16.check())
def check_fields():
	import retain16 as RT,descent_witness as W
	R=W.R;Hr=[x%R for x in H]
	def res(a,b):
		a=RT._trim([x%R for x in a]);b=RT._trim([x%R for x in b]);s=1
		while b and len(b)>1:
			da,db=len(a)-1,len(b)-1;rem=RT.pmod(a,b,R)
			if not rem:
				return 0
			s=s*(-1)**(da*db)*pow(b[-1],da-(len(rem)-1),R)%R;a,b=b,rem
		return s*pow(b[0],len(a)-1,R)%R if b else 0
	def norm_poly(terms,deg):
		n=6*deg+1;xs=list(range(n))
		ys=[res(Hr,[sum(c[i]*pow(k,d,R)for d,c in terms if i<len(c))%R for i in range(6)])for k in xs]
		for j in range(1,n):
			for i in range(n-1,j-1,-1):
				ys[i]=(ys[i]-ys[i-1])*pow(xs[i]-xs[i-j],-1,R)%R
		poly=[ys[-1]]
		for i in range(n-2,-1,-1):
			poly=[(-xs[i]*poly[0]+ys[i])%R]+[(poly[k-1]-xs[i]*poly[k])%R for k in range(1,len(poly))]+[poly[-1]]
		return poly
	out={}
	for name,terms,deg in(('f',C.PTERMS,12),('g',RT.G_TERMS,2)):
		poly=norm_poly(terms,deg)
		out['deg_'+name]=len(poly)-1
		out['sqf_'+name]=poly[-1]==1 and len(poly)==6*deg+1 and len(RT.pgcd(poly,RT.pderiv(poly,R),R))==1
	out['unramified']=out['sqf_f']and out['sqf_g']
	out['disc_h']=C.norm(C.elt([0,-4,-6,8,0,6]))
	out['sig']=signatures()
	return out
def signatures():
	def tr(p):
		while p and p[-1]==0:p=p[:-1]
		return p
	def sturm(p):
		s=[tr([F(c)for c in p])];s.append([i*c for i,c in enumerate(s[0])][1:])
		while len(s[-1])>1:
			a,b=s[-2][:],s[-1]
			while len(a)>=len(b):
				q=a[-1]/b[-1];d=len(a)-len(b);a=tr([x-q*b[i-d]if i>=d else x for i,x in enumerate(a)][:-1])
			if not a:break
			s.append([-x for x in a])
		return s
	def ch(s,x):
		v=[t for t in(sum(c*x**i for i,c in enumerate(q))for q in s)if t]
		return sum((a>0)!=(b>0)for a,b in zip(v,v[1:]))
	nre=lambda p,B=F(10)**60:(lambda s:ch(s,-B)-ch(s,B))(sturm(p))
	sh=sturm(H);ys=[]
	def iso(a,b):
		n=ch(sh,a)-ch(sh,b)
		if n==1 and b-a<F(1,10**30):ys.append((a+b)/2)
		elif n:iso(a,(a+b)/2);iso((a+b)/2,b)
	iso(F(-10),F(10))
	rf=rfp=rg=0
	for y in ys:
		fx=[0]*13
		for d,c in C.PTERMS:fx[d]=sum(F(ci)*y**i for i,ci in enumerate(c))
		rf+=nre(fx);rfp+=nre(fx[0::2]);rg+=nre([U*(y+1),0,1])
	rk=lambda n,r1:r1+(n-r1)//2-1
	return dict(f=(rf,36-rf//2),fp=(rfp,18-rfp//2),g=(rg,6-rg//2),odd=rk(72,rf)-rk(36,rfp)+rk(12,rg)-rk(6,len(ys)),sm=(rk(72,rf)-rk(36,rfp)+rk(12,rg)-rk(6,len(ys)))//2)
def hx_toy(seed,lu=12,P=2**61-31):
	import unit_entries as UE
	m,n,g=M.BW_M,M.BW_N,Random(seed);c,_,R,_=filtsim.build_system(69.0,(UE.MUF/69,UE.MUG/69),M.FILT_EXCESS,seed,lu);UE.MID.clear()
	filtsim.merge(filtsim.purge(R,len(c),M.N_SM,False)[0],len(c),(int(M.W_ROW)-1,),kmax=64,verbose=False,hook=UE)
	i4=next(x for x in(pow(a,(P-1)//4,P)for a in range(2,99))if x*x%P==P-1)
	R=[{k:(int(v.real)+int(v.imag)*i4)%P for k,v in r.items()if k>=0}for r in UE.LAST];ix={k:j for j,k in enumerate(sorted({k for r in R for k in r}))};D=len(R)
	A=[[(ix[k],v)for k,v in r.items()]+[(j,g.randrange(P))for j in range(len(ix),D-1)]for r in R]
	A=[r+[(D-1,-sum(x for _,x in r))]for r in A];T=[[]for _ in A]
	for i,r in enumerate(A):
		for j,x in r:T[j].append((i,x))
	def kr(B,V,d):
		o=[]
		for _ in range(d):o+=V;V=[[sum(x*v[j]for j,x in r)%P for r in B]for v in V]
		return o
	def rk(H,r=0):
		for c in range(len(H[0])):
			q=next((i for i in range(r,len(H))if H[i][c]),-1)
			if q<0:continue
			H[r],H[q]=H[q],H[r];iv=pow(H[r][c],-1,P);H[r]=[x*iv%P for x in H[r]]
			for i in range(r+1,len(H)):
				if H[i][c]:f=H[i][c];H[i]=[(a-f*b)%P for a,b in zip(H[i],H[r])]
			r+=1
		return r
	K=kr(A,[[g.randrange(P)for _ in A]for _ in range(n)],-(-D//n)+1);pm=list(range(D));g.shuffle(pm);Xb=[[0]*D for _ in range(m)]
	for t,j in enumerate(pm):Xb[t*m//D][j]=g.randrange(P)
	return D,rk([[sum(a*b for a,b in zip(u,v))%P for v in K]for u in kr(T,Xb,-(-D//m)+1)])
def main():
	cfg_in=json.load(sys.stdin)if not sys.stdin.isatty()else{}
	seed=int(cfg_in.get("seed",1))
	say(f"evidence: plain STNFS core; seed={seed}; unit {M.UNIT:.0f} instructions")
	ok0=check_resultant();ok0b=check_tower();ok1=check_rho_and_semismooth()
	cfg,pd=(lambda s:(s["config"],s["pool_d"]))(load_stats())
	ok2,n_can,anchors,nunits=verify_pool(cfg,pd)
	ok2b,nb2,nr2=box_check(cfg,anchors,seed=seed)
	say(f"[2b] 6 boxes ({nr2} rows): {nb2} exact f norms inside leaf bounds {K(ok2b)}")
	theta_pool=F(cfg['production_theta']);theta_reg=F(cfg['production_theta_region'])
	assert theta_reg<=theta_pool==M.OP['theta']and F(pd['theta_region'])==theta_reg
	ND=pd['draws_per_seed'];sc=float((theta_reg/theta_pool)**3)
	def region(th):
		rmin=math.isqrt(int(th*cfg['production_A_squared']));corr=LCB.correction(cfg,rmin)
		return rmin,corr,M.Geometry(cfg,th,F(corr['safe_yield_factor']),F(corr['safe_work_population_factor']))
	rmin,corr,geo=region(theta_reg)
	def tot(seeds,route=False):
		n=sum(pd['draws'][x]for x in seeds);parts=[]
		for x in seeds:
			sm=dict(pd['sums'][x])
			if route:sm['work']-=(0.5-M.CELL_ROUTE)*sm['cells']
			parts.append(sm)
		return M.pool_totals(M.add_sums(parts),n*sc,n,geo)
	tD,tW,tB=tot(pd['judged']),tot(pd['work'],True),tot(pd['sizing'])
	th0=theta_reg-F(1,8000);nB=sum(pd['draws'][x]for x in pd['sizing'])
	tB0=M.pool_totals(pd['sizing_prev'],nB*float((th0/theta_pool)**3),nB,region(th0)[2])
	op=M.OP
	retain,cost=M.retentions(op)
	for r,an in zip(anchors,pd['anchors']):r['evals_w'],r['leaves_w']=an[14],an[15]
	ws=[0.0 if M.skip(r)else M.evaluate([r],1.0,geo,op,retain,cost)['sums']['work']/geo.weight(r['qa'])for r in anchors]
	drop=(10.01-M.REC_ROUTE)*M.NBAR*math.log(M.OP['lx2']/M.P_MID)*M.SIEVE_SIDES
	ok_w=all(math.isclose(w,an[13]-(0 if M.skip(r) else 0.5-M.CELL_ROUTE+drop),rel_tol=1e-6)for w,an,r in zip(ws,pd['anchors'],anchors))
	ok_w=ok_w and all(list(M.pool_row_warm(cfg,an[0],an[1],an[11])if M.RUN_T[0]<=an[12]-an[11]+1<=M.RUN_T[1]else(an[6],an[8]))==an[14:16]for an in pd['anchors'])
	qq,dd,_=C.integer_metric(cfg['matrix_a']);lim=int(theta_reg*cfg['production_A_squared']*dd);ex=[]
	for r in sorted(anchors,key=lambda r:r['T'])[:3]:
		a,lo,hi=r['run'];can=[C.quadratic(qq,x)<=lim and x[5]>=1 and canonical(cfg,x)for x in([a[0]+t]+a[1:]for t in range(lo-1,hi+2))]
		ex.append((r['T'],all(can[1:-1]),can[0],can[-1]))
	ok_r=all(e[1]for e in ex)
	say(f"[2w] anchors' work {K(ok_w)}; [2r] 3 anchor runs "
	f"(T, canonical, ends) {ex} {K(ok_r)}")
	say(f"[3] theta={cfg['production_theta_region']} a-radius {rmin} (b-slices >= {math.isqrt(int((F(cfg['production_region_lam'])-F(cfg['production_region_kap'])*theta_reg)*cfg['production_A_squared']))}); "
	f"a 2^{log2(geo.va):.3f}, half-slice <= 2^{log2(geo.vb*geo.rad(0.0)**3):.3f}; certified "
	f"{float(corr['safe_yield_factor']):.9f} (yield), {float(corr['safe_work_population_factor']):.9f} (work) at R >= {rmin}; "
	f"rows {sum(pd['sums'][x]['rows']for x in pd['judged'])}/{sum(pd['sums'][x]['rows']for x in pd['work'])}")
	log2_U=M.log2_fb(op['lb']);need=log2_U+log2(op['excess'])
	e=dict(tD,log2_U=log2_U,log2_need=need,log2_rels_est=tD['log2_rels'],log2_rels=tD['log2_rels']-M.REL_BAND_LOSS,margin=tD['log2_rels']-M.REL_BAND_LOSS-need,log2_relcol=tW['log2_est']+log2(M.CANON_WORK_ALLOWANCE),
	per_cell=tW['per_cell']*M.CANON_WORK_ALLOWANCE,log2_population=tW['log2_population'])
	e['feasible']=e['margin']>=max(M.MARGIN,M.MARGIN_SE*tD['yield_se_log2'])
	rule=lambda t:M.MARGIN+sqrt(2)*t['yield_se_log2']
	ok_size=tB['log2_rels']-need>=rule(tB)and tB0['log2_rels']-need<rule(tB0)
	say(f"[4s] margin {tB['log2_rels']-need:.6f} >= 0.06 + sqrt(2) SE = {rule(tB):.6f} at {cfg['production_theta_region']}; "
	f"{tB0['log2_rels']-need:.6f} < {rule(tB0):.6f} at {th0*8000}/8000: smallest {K(ok_size)}")
	work_se=tW['work_se'];allow_need=1+M.WORK_SE_K*work_se
	ok_allow=M.CANON_WORK_ALLOWANCE>=allow_need
	say(f"[4a] work {pd['work']}: row SE {100*tW['cond_se']:.4f}%, total SE {100*work_se:.4f}%; "
	f"allowance {M.CANON_WORK_ALLOWANCE:.5f} >= 1 + {M.WORK_SE_K:.2f} SE = {allow_need:.6f} {K(ok_allow)}")
	ok_relmeans=e['nf_rel']>=69*DS.P.MAIN_UF and e['ng_rel']>=69*DS.P.MAIN_UG
	rels_lo=e['log2_rels']+log2(1-2*e['yield_se_log2']*math.log(2))
	m108=e['margin']-log2(M.FILT_EXCESS/op['excess']);ok_rels_input=rels_lo>=DS.P.LOG2_RELS and log2_U+log2(M.FILT_EXCESS)<=rels_lo and m108>=max(M.MARGIN,M.MARGIN_SE*e['yield_se_log2'])
	rf=retain.at(e['nf_mean'],op['slack_f']);rg=retain.at(e['ng_mean'],op['slack_g'])
	say(f"[4] judged {pd['judged']}: B=2^{op['lb']:.0f}, X2=2^{op['lx2']:.0f}, slacks ({op['slack_f']}, {op['slack_g']}), excess {op['excess']}, "
	f"alpha ({op['alpha_f']}, {op['alpha_g']}); mean norms "
	f"{e['nf_mean']:.1f}/{e['ng_mean']:.1f} (weighted {e['nf_rel']:.2f}/{e['ng_rel']:.2f}); retain f {rf['retention']:.4f} (g exact); "
	f"rels 2^{e['log2_rels_est']:.5f} - {M.REL_BAND_LOSS} = 2^{e['log2_rels']:.5f} vs need 2^{e['log2_need']:.5f}, margin {floor(e['margin']*10000)/10000:.4f} >= max({M.MARGIN}, {M.MARGIN_SE:g} SE) "
	f"({'feasible' if e['feasible']else 'INFEASIBLE'}); SE {e['yield_se_log2']:.5f} bits; "
	f"LOG2_RELS {DS.P.LOG2_RELS:.3f} and {M.FILT_EXCESS} U = 2^{log2_U+log2(M.FILT_EXCESS):.5f} <= 2-SE count {rels_lo:.5f} (margin {m108:.4f}) {K(ok_rels_input)}")
	la=M.linear_algebra(e['log2_U'],op['excess'])
	say(f"[5] 2^{e['log2_population']:.3f} cells at {e['per_cell']:.1f}/cell (routing {M.REC_ROUTE:.0f}/record, "
	f"{M.CELL_ROUTE:g}/cell; apply {M.C_APPLY:.2f}): collection 2^{e['log2_relcol']:.4f} units; "
	f"N(W)/U {M.NU} -> N = 2^{la['log2_N']:.4f}; LA 2^{la['log2_la']:.4f} units (k5 {M.C_ENTRY_K})")
	from threshold_model import primes_upto
	ntd=len(primes_upto(2**20));lf_,lg_=[-(-b//64)for b in M.BAND_NBITS];per_p=5*(lf_+lg_)+2*22+4;tdiv=ntd*per_p+12*(lf_*M.BAND_NBITS[0]+lg_*M.BAND_NBITS[1])
	ok_tdiv=ntd==82025 and M.PRIME_MIN<=2**20 and tdiv<=M.C_TDIV
	say(f"[5t] trial division: 2^{log2(tdiv):.3f} <= C_TDIV 2^{log2(M.C_TDIV):g} {K(ok_tdiv)}")
	Us,res=filtsim.run(op['lb'],e['nf_mean']/op['lb'],e['ng_mean']/op['lb'],M.FILT_EXCESS,
	seed,12,kmax=64,targets=(100,),verbose=False)
	r100=res[100]
	say(f"[6] (2^12, seed {seed}): N/U = {1/r100[0]:.4f} at weight {r100[3]:.0f} (max coef {r100[2]}); "
	f"FILT_RUNS seeds: "+", ".join(f"2^{k} {len(v)}"for k,v in sorted(M.FILT_RUNS.items())))
	rec=filtsim.missing_recon(seed,excess=M.FILT_EXCESS)
	ok_rec=all(mf<=md for(lo,hi,n,mf,md,h3)in rec)
	say(f"[6b] (2^12, seed {seed}) per lambda band "
	"(ideals, measured, model, lost >= 3): "
	+"; ".join(f"[{lo:g},{hi:g}): {n}, {mf:.4f}, {md:.4f}, {h3}"for(lo,hi,n,mf,md,h3)in rec)+f" {K(ok_rec)}")
	ok_counts=ball_counts_check()
	import retain16 as RT0
	RT0.load_table(load_stats())
	Pm=DS.P()
	dres=DS.evaluate(Pm,DS.PLAN16,DS.BOOT16)
	two_log=dres['two_log_units']
	bt=dres['boot']
	lz=sorted((v for v in dres['DET'].values()if 'lazy' in v),key=lambda v:-v['survcost_log2'])[:3]
	lz_rows=[(round(v['mu'],1),v['nb_r'],round(v['psi'],3),round(v['lazy'],4))for v in lz]
	ok_lazy=bt['own_log2']<bt['sub_log2']
	say(f"[7l] {lz_rows}; lazy boot: {log2((1-bt['psi']**Pm.BOOT_RELS)/(1-bt['psi'])):.3f} bits vs log2 {Pm.BOOT_RELS} {K(ok_lazy)}")
	cap=M.DESCENT_CAP
	markov=2**(two_log-cap)
	log2_EN,over,coll,pmin=repeat_bound(dres);OVER=2.0**over;REP=2.0**coll
	p_desc=1-2*bt['phi']-markov-OVER-REP
	p_main=1.0-2.0**(la['log2_N']+2-128)-M.P_CAP
	success=p_main*p_desc;success_head=p_main*(p_desc-0.1*markov)
	rows=DS.report(dres)
	sv=[dres['DET'][k]for k in dres['DET']if dres['DET'][k]['how']=='sieve']
	hmax=max(DS.CELLS[d['R']]for d in sv)
	ok_width=hmax<1<<103 and(hmax-1)>>40<1<<63 and all(16*max(d['nf'],d['ng'])<1<<16 for d in sv)
	mu_min=min(r['mu']for r in rows);psi_max=max(max(r['psi']for r in rows),bt['psi'])
	say(f"[7] descent: balls {K(ok_counts)}; {len(rows)} nodes, min mu {floor(mu_min*10)/10:.1f}, max psi {ceil(psi_max*1000)/1000:.3f}; "
	f"boot n1={bt['n1']} nA={bt['nA']:.2f} kA={bt['kA']} R={bt['R']}: "
	f"own 2^{bt['own_log2']:.2f} + subtree 2^{bt['sub_log2']:.2f} x retry {bt['retry']:.3f}; one log 2^{bt['one_log_log2']:.3f}; two logs 2^{two_log:.3f} units; "
	f"cap 2^{cap:g}: Markov <= 2^{ceil((two_log-cap)*1000)/1000:.3f}; boot fails < {ceil(bt['phi']*1e7)/1e7:.2e}; E[N] <= 2^{log2_EN:.2f} (N_MAX 2^{M.NODE_CAP_LOG2:g}), abort <= 2^{over:.2f}; repeats (smallest bin 2^{pmin:.2f}): 2^{coll:.2f}; success {floor(success*10000)/10000:.4f}, {floor(success_head*10000)/10000:.4f} at +10%; widths: half <= "
	f"2^{log2(hmax):.2f}, owner < 2^{log2((hmax-1)>>40):.2f}, score {16*max(max(d['nf'],d['ng'])for d in sv):.0f} {K(ok_width)}")
	nb,np_,nbf=nbar(14,15)
	say(f"[9] {np_} primes in [2^14, 2^15): {nb:.3f}, f only {nbf:.3f} (NBAR {M.NBAR}, sieved {M.NBAR*M.SIEVE_SIDES})")
	ok9=nb<=M.NBAR*1.01 and nbf<=M.NBAR*M.SIEVE_SIDES*1.01
	tail,red9,chi,tail_n,tail_ok=la_row_tail();ok_chi=chi<=M.C_CHI
	say(f"[10] tail: {tail_n} (straight-line {tail_ok}); {tail} <= C_ROW {M.C_ROW:.0f}; reduce {red9} <= C_RED9 {M.C_RED9:.0f}; chi x {chi} <= C_CHI {M.C_CHI:.0f} {K(ok_chi)}")
	def fst(k):
		x=[v[1]/2**k for v in M.FILT_RUNS[k]];n=len(x);m=sum(x)/n
		return m,(sum((y-m)**2 for y in x)/(n-1))**.5/n**.5,n
	fs={k:fst(k)for k in M.FILT_RUNS};fx={k:(lambda x:(sum(x)/len(x),(sum((y-sum(x)/len(x))**2 for y in x)/(len(x)-1)/len(x))**.5,len(x)))([v[1]/2**k for v in vs])for k,vs in M.FILT_X.items()};fxk=tuple(k for k in(19,20)if k in fx);fneed=max(fx[k][0]+2*fx[k][1]for k in fxk)
	sdu=max((sum((v[1]/2**k-sum(w[1]for w in vs)/len(vs)/2**k)**2 for v in vs)/(len(vs)-1))**.5*2**(k/2)for T in(M.FILT_RUNS,M.FILT_X)for k,vs in T.items()if len(vs)>=6)
	lu=la['log2_N']-log2(M.NU);chb=(M.SD_C*2**(-lu/2)/(M.NU-M.NU_H))**2;cap_ok=sdu<=M.SD_C and M.NU_H>=fneed and chb+1e-15<=M.P_CAP
	say(f"[10c] cap {M.NU} U: sd {sdu:.3f} <= {M.SD_C}; {M.NU_H} >= {fneed:.6f}; {chb:.2e} <= {M.P_CAP}: {cap_ok}")
	nu_ok=cap_ok and M.NU>=fneed and max([fs[k][0]for k in(19,20)if k in fs]+[fx[k][0]for k in fxk])<=fs[15][0]and max(v[2]for vs in(*M.FILT_RUNS.values(),*M.FILT_X.values())for v in vs)<=255
	say(f"[10] N(W)/U (mean/SE/seeds): "+", ".join(f"2^{k} {m:.5f}/{s_:.5f}/{n}"for k,(m,s_,n)in sorted(fs.items()))
	+f"; excess {M.FILT_EXCESS} (PREREG-J): "+", ".join(f"2^{k} {m:.5f}/{s_:.5f}/{n}"for k,(m,s_,n)in sorted(fx.items()))+f"; {M.NU} >= max 2^19, 2^20 mean + 2 SE {fneed:.6f}; 2^19, 2^20 means (both excesses) <= 2^15's; "
	f"max coef {max(v[2]for vs in(*M.FILT_RUNS.values(),*M.FILT_X.values())for v in vs)} <= 255: {nu_ok}")
	fy={k:[v[1]/2**k for v in vs]for k,vs in M.FILT_Y.items()};fyn={k:(sum(x)/len(x),(sum((y-sum(x)/len(x))**2 for y in x)/(len(x)-1)/len(x))**.5,len(x))for k,x in fy.items()}
	pr={k:[v[1]/v[2]for v in vs]for k,vs in M.FILT_Y.items()}
	prs={k:(sum(x)/len(x),(sum((y-sum(x)/len(x))**2 for y in x)/(len(x)-1)/len(x))**.5)for k,x in pr.items()};pw=sum(prs[k][1]**-2 for k in prs);pm=sum(prs[k][0]*prs[k][1]**-2 for k in prs)/pw;ps=pw**-.5
	ok_fy=all(len(pr[k])==len(M.FILT_Y[k])for k in pr)and pm+2*ps<=1 and max(v[3]for vs in M.FILT_Y.values()for v in vs)<=255
	say(f"[10y] excess {M.FILT_EXCESS}, means {M.FILT_Y_MEANS} vs campaign: "+", ".join(f"2^{k} {m:.5f}/{s_:.5f}/{n}, ratio {prs[k][0]:.4f} +- {prs[k][1]:.4f}"for k,(m,s_,n)in sorted(fyn.items()))+f"; pooled {pm:.4f} +- {ps:.4f}, + 2 SE {pm+2*ps:.4f} <= 1 {K(ok_fy)}")
	hxD,hxr=hx_toy(seed);ok_hx=hxr==hxD
	say(f"[10x] H-X rank {hxr} = {hxD} {K(ok_hx)}")
	bf_n,bf_br=check_la_kernel("bw_dot.s","LBB1_2:")
	say(f"[10] bw_full: {bf_n} (charged {M.C_FULL:.0f}) {K(bf_n==97 and bf_n<=M.C_FULL and bf_br==['b.ne'])}")
	wl="\n".join(_asm("wide12.txt"));ml="\n".join(_asm("measurements12x.txt"))
	fr=lambda k:{int(p):float(v)for p,v in re.findall(k+r" p=(\d+): Ir \d+ -> \d+, reps \d+ -> \d+, Ir per rep ([\d.]+)",wl)}
	ws,wr=fr("wsetup6w"),fr("wrec");nar=wr.get(2**63-25,0);big=[v for p,v in wr.items()if p>>63]
	pm=int(re.findall(r"Ir per pair min \d+ mean \d+ max (\d+)",wl)[0])
	recw={int(p):F(int(b)-int(a),int(d)-int(c))for p,a,b,c,d in re.findall(r"TAB rec (\d+) (\d+) (\d+) (\d+) (\d+)",ml)}
	ok_wide=(len(ws)==2 and max(ws.values())<=406.87<=M.C_SETUP6*M.WIDE_SETUP and len(big)==3 and 2.5*max(abs(v-nar)for v in big)<=(M.WIDE_REC-1)*nar
	and wl.count("MATCH")==6 and wl.count("certified 200")==3 and "0 mismatches" in wl and "failures 0" in wl and "BYTE-IDENTICAL" in wl
	and M.LOG2_WIDE_LAT>=log2(6*1.25506*2**64/(64*math.log(2))*4*pm/M.UNIT)
	and[float(recw[p])for p in(257,65537,8388617,67108879,268435459)]==[x[1]for x in M.REC_TAB_W])
	say(f"[10w] setup6w {max(ws.values()):.2f} <= 406.87 <= {M.C_SETUP6*M.WIDE_SETUP:.2f}; "
	f"wide p within {max(abs(v-nar)for v in big):.2f} of {nar:.2f}; "
	f"per hit {[round(float(recw[p]),2)for p in sorted(recw)]} == REC_TAB_W; lattices 2^{log2(6*1.25506*2**64/(64*math.log(2))*4*pm/M.UNIT):.3f} <= 2^{M.LOG2_WIDE_LAT} units; "
	f"log checks {K(ok_wide)}")
	ol="\n".join(_asm("measurements-oh.txt"))
	br=[(int(e),int(c),int(sg),int(a.replace(',','')),float(pa),int(z.replace(',','')))for e,c,sg,a,pa,z in re.findall(r"entries=(\d+) rows=\d+ what=build cells=(\d+) K=\d+ L4=\d+ segments=(\d+).*\n *([\d,]+) \(([\d.]+)%\)  \?\?\?:build_layout'2.*\n *([\d,]+) \(",ol)]
	ok_build=len(br)==4 and all(a/e<=M.C_BUILD_E and a/(pa/100)-a-z<2**20 and(a/(pa/100)-a-z)/sg<=M.C_BUILD_SEG for e,c,sg,a,pa,z in br)
	say(f"[10o] layout/entry "+", ".join(f"{a/e:.3f}"for e,c,sg,a,pa,z in br)+f" <= {M.C_BUILD_E}; run remainder < 2^20 "
	f"and per segment <= 2^17: {K(ok_build)}")
	Lx=_asm("measurements12x.txt");q12={};grp={}
	for ln in Lx:
		f=ln.split()
		if ln.startswith("TAB line "):q12[int(f[2])]=(int(f[4])-int(f[3]))/(int(f[6])-int(f[5]))
		if ln.startswith("MP "):
			d=grp.setdefault((int(f[1]),int(f[2])),{}).setdefault(int(f[3]),[0,0,set(),0,0])
			d[0]+=int(f[6])-int(f[5]);d[1]+=int(f[8])-int(f[7]);a_,b_=map(int,f[4].split('-'));d[2]|=set(range(a_,b_+1));d[3]=int(f[10]);d[4]=int(f[11])*(int(f[13])-int(f[12]))
	ok_line=len(q12)==4 and all(c>=q12[p]for(h,c),p in zip(M.LINE_TAB[1:],(257,4099,65537,1048583)));lmsg=[]
	for(lp,pt),sds in sorted(grp.items()):
		hs=[log2(d[1]/(2*d[4]*d[3]))for d in sds.values()];cs=[d[0]/d[1]for d in sds.values()];n_=len(cs);mc=sum(cs)/n_
		ub=sum(d[0]for d in sds.values())/sum(d[1]for d in sds.values())+M.LINE_SE_K*sqrt(sum((c-mc)**2 for c in cs)/(n_-1)/n_)if n_>1 else float('inf')
		ok_line&=all(d[2]==set(range(1,d[3]+1))for d in sds.values())and n_>=3 and ub<=M.tab(M.LINE_TAB,max(hs));lmsg.append(f"{lp}/{pt}: {ub:.2f} <= {M.tab(M.LINE_TAB,max(hs)):.2f}")
	sp=[ln.split()for ln in Lx if ln.startswith("SP ")]
	ok_line&=len(sp)>0 and all(M.tab(M.LINE_TAB,log2(int(f[8])/int(f[9])/int(f[2])))>=(int(f[5])-int(f[4]))/(int(f[7])-int(f[6]))for f in sp)
	say(f"[10l] line: base {[round(q12[p],3)for p in sorted(q12)]} <= LINE_TAB; groups "+", ".join(lmsg)+f"; {len(sp)} small-p runs: {K(ok_line)}")
	okr,rmsg=kernel_rec7(pd['runs'],nunits)
	import runs16 as R16r
	qq,dd,_=C.integer_metric(cfg['matrix_a']);lim=int(theta_reg*cfg['production_A_squared']*dd)
	say(f"[10r] rec7: {rmsg} {K(okr)}")
	kt="\n".join(_asm("fft.txt"));mm=[(int(a),float(c))for ln in kt.splitlines()if ln.startswith("FFTPTS ")for a,c in(x.split(":")for x in ln.split()[1:])];bc=M.gbatch_cost(pbits=1.000028*2**48/M.LN2,leaf=13,L=39,c_leaf=2**19);gc=M.gbatch_cost()
	sr=[x.split(":")for ln in kt.splitlines()if ln.startswith("SRTWAVES ")for x in ln.split()[1:]];okb=len(mm)>=25 and max(n for n,_ in mm)>=2**26 and all(c<=M.fft_F(n)for n,c in mm)and len(sr)>=10 and all(x[2]=="0"for x in sr)and bc<=M.C_BATCH and gc<=M.C_GBATCH and M.gbatch_cost(pbits=1.000028*2**64/M.LN2,leaf=13,L=54,c_leaf=2**19)<=M.C_BATCH2
	say(f"[10b] FFT: {len(mm)} products, charged/F <= {max(c/M.fft_F(n)for n,c in mm):.4f}; {len(sr)} SRT waves; band 2^{log2(bc):.4f} <= C_BATCH 2^{log2(M.C_BATCH):g}; [10g] g batch 2^{log2(gc):.4f} <= C_GBATCH 2^{log2(M.C_GBATCH):g} {K(okb)}")
	sb,sf=M.srt_toy(8,seed);ok_srt=sb==0 and sf<=2
	say(f"[10v] SRT: 2^8 leaves, Newton offsets -2/+2: {sb} mismatches vs %, {sf} corrections {K(ok_srt)}")
	rmf,rmg=M.rmin(retain,op['slack_f']),M.rmin(retain,op['slack_g']);sp=M.band_split(e['log2_rels_est'],e['yield_se_log2'],rmf,rmg);rj=json.load(open(os.path.join(_DIR,'refinement.json')));ok_split=sp['log2_units']<=M.LOG2_SPLIT==rj['ledger_upper_rows'][-1]and rj['split_inputs'][0]>=e['log2_rels_est']and rj['split_inputs'][1]>=e['yield_se_log2']
	say(f"[10t] band split: 2^{sp['c_rel_log2']:.4f} per relation; "
	f"relations <= 2^{sp['log2_count']:.4f}: 2^{sp['log2_units']:.4f} <= 2^{M.LOG2_SPLIT} {K(ok_split)}")
	D7=M.contam_sup(geo);g7=0.0
	for i in range(2001):
		S7=geo.slice_cells(float(theta_reg)*i/2000);c7=0.0;lp=M.P7
		while lp<M.X7:
			lp2=min(lp+0.25,M.X7);lh=log2(S7)-0.5*(lp+lp2);dp=(M.NBAR/2)*(2**M.log2_li(lp2)-2**M.log2_li(lp));wd=lp>=M.WIDE_FROM
			c7+=(S7*M.NBAR*math.log(lp2/lp)*(M.WIDE_REC*M.tab(M.REC_TAB_W,lh)if wd else M.tab(M.REC_TAB,lh))if lh>=M.REC_TAB[-1][0]else dp*M.C_COSET_EMPTY)+dp*M.C_SETUP6;lp=lp2
		g7=max(g7,c7)
	ct=pd['contam'];ex7=[]
	for sd,i,T_ in ct['replay']:
		a=C.regenerate_pair(cfg,sd,i)[0];lo,hi=R16r.run(a,R16r.unit_list(),qq,lim)
		ex7.append(hi-lo+1==T_ and all(C.quadratic(qq,x)<=lim and x[5]>=1 and canonical(cfg,x)for x in([a[0]+t]+a[1:]for t in range(lo,hi+1))))
	cr=float((F(ct['theta'])/theta_reg)**3);ok_c7=D7>=g7 and ct['bad']==0 and M.RUN_ALLOW>=(1-0.05**(1/ct['n']))*cr and all(ex7)
	say(f"[k7c] D_SUP 2^{log2(D7):.3f} >= 2^{log2(g7):.3f}; {ct['n']} runs ({ct['selection_sha']}): "
	f"{ct['points']} a, {ct['bad']} bad, RUN_ALLOW {M.RUN_ALLOW} >= {1-0.05**(1/ct['n']):.6f} x {cr:.4f}; replay {ex7} {K(ok_c7)}")
	oh=[ln.split()for ln in ol.splitlines()if ln.startswith("OH ")]
	sc=[(int(f[2]),(int(f[3])-int(f[4]))/int(f[5]))for f in oh if f[1]=="scan"]
	cr={int(f[1]):(int(f[2]),int(f[3])-int(f[4])-int(f[5])-int(f[6]))for f in(ln.split()for ln in ol.splitlines()if ln.startswith("CELLRES "))};cres=(cr[26][1]-cr[24][1])/(cr[26][0]-cr[24][0]);ok_oh_cr=cres<=M.CELL_ROUTE/4
	ap=[(int(f[2])-int(f[3]))/int(f[4])for f in oh if f[1]=="apply"]
	need=max(v*b-M.C_SCAN_CELL*b for b,v in sc)
	ok_oh=ok_oh_cr and set(range(1,65))<={b for b,v in sc}and need<=M.C_SCAN_LEAF and len(ap)==1 and ap[0]<=M.C_APPLY and M.C_ZERO==0 and ol.count("MATCH | survivors")>=4
	say(f"[10s] cell residual {cres:.3f}; scan2: {len(sc)} points: leaf <= {M.C_SCAN_CELL} x L + {M.C_SCAN_LEAF} (need {need:.3f}); "
	f"apply {ap[0]:.3f} <= {M.C_APPLY}; MATCH: {K(ok_oh)}")
	b2_n,b2_br=check_la_kernel("bw_dot2.s","LBB0_2:")
	k5_n,k5_br=check_la_kernel("krylov5.s",".LBB0_1:");g1_n,g1_br=check_la_kernel("krylov5.s",".LBB1_1:")
	u5=k5_n/5;UE=__import__('unit_entries');ku={k:(sum(v)/len(v),(sum((x-sum(v)/len(v))**2 for x in v)/(len(v)-1)/len(v))**.5)for k,v in M.K5_RUNS.items()};k5c=max(ku[k][0]+2*ku[k][1]for k in(16,17))
	k5_ok=k5_br==['b.ne']and g1_br==['b.ne']and u5==UE.U5==7.8 and g1_n==UE.G1 and M.C_ENTRY_K>=k5c and max(ku[16][0],ku[17][0])<=ku[15][0]and max(M.K5_RUNS[16]+M.K5_RUNS[17])<=12.54
	say(f"[10] krylov5: {k5_n}/entry, {u5:.1f}, {g1_n}; PREREG-U k5 "+", ".join(f"2^{k} {ku[k][0]:.4f}/{ku[k][1]:.4f}"for k in sorted(ku))+f"; max mean+2SE {k5c:.4f} <= {M.C_ENTRY_K} {K(k5_ok)}")
	say(f"[10] bw_dot2: {b2_n} (charged {M.C_PROJ:.0f}) {K(b2_n<=M.C_PROJ and b2_br==['b.ne'])}")
	au=audit(npr=dres['DET'][((220,225),'g')]['pr_log2'],bpr=bt['pr_log2'])
	ok11=au['boot']['within_model_bins']and au['specialq_within_bins']
	say(f"[11] audit: boot lattice ({au['boot']['swaps']} swaps), cg max {max(s['cg_max']for s in au['boot']['stats'])} "
	f"<= {max(x[1]for x in DS.P.BOOT_BINS):.1f} {K(au['boot']['within_model_bins'])}; special-q lattices "
	+"; ".join(f"{s['side']}@2^{s['bits']}: cf max {max(x['cf_max']for x in s['stats'])}, cg max {max(x['cg_max']for x in s['stats'])}"for s in au['specialq'])
	+f" in bins (cf <= {max(b[1]for b in DS.P.QBINS):.1f}, cg <= {max(b[2]for b in DS.P.QBINS):.1f}) {K(au['specialq_within_bins'])}; "
	f"UDUP yr: 220g {au['relw']['node_yield_ratio']:.2f}, boot {au['relw']['boot_yield_ratio']:.3f} vs {DS.P.UDUP}; "
	f"chi2 {au['sampler']['chi2']:.1f} (< 60 {K(au['sampler']['ok'])}); "
	f"line <= 2^{au['storage']['line_cells_log2_max']}; chunks {au['storage']['routing_stable']}")
	ok11r=au['relw']['node_yield_ratio']>1.0 and au['relw']['boot_yield_ratio']>0.70
	ok11s=au['sampler']['ok']
	import descent_witness as W
	import retain16 as RT
	import ecm16 as EC
	st=load_stats()
	nb=st['nb_empty_window']
	rr=[(b,s_,R_)+W.retention(b,s_,R_,seed+b,40)for b,s_,R_ in[(b,s_,dres['DET'][(bb,s_)]['R'])for b,s_ in((69,'g'),(150,'g'))for bb in dres['bins']if bb[0]==b]]
	nb_ok=all(z<=zp for b_,s_,R_,lb_,m_,r_,z,zp in nb['rows']if m_>=7.4)
	tab=st['prime_structure']['table']
	special=[p for p in RT.primes(2,512)if any(e!=1 or any(m!=1 for m in mf+mg)for d1,e,mf,mg in tab[str(p)])]
	spot=sorted(set([p for i,p in enumerate(RT.primes(2,512))if i%7==seed%7]+special))
	ok_tab=len(tab)==len(RT.primes(2,512))and all(RT.structure(p)==tab[str(p)]for p in spot)
	sample={}
	for p in special:
		rs=W.residue_valuations(p,400,seed+p)
		for side in 'fg':
			bound=1-RT.local_pdf(p,tab[str(p)],side,False,True).get(0,0.0)
			frac=1-rs[side].get(0,0)/rs['n']
			sample[f"{p}{side}"]=(round(frac,4),round(bound,4),frac<=bound+3*sqrt(bound*(1-bound)/rs['n'])+1e-12)
	ok_sample=all(v[2]for v in sample.values())
	pf,ordf=RT.omitted_pdf('f',True);pg,ordg=RT.omitted_pdf('g',True)
	ex=lambda pdf:sum(x*v for x,v in pdf.items())/RT.GRID
	ret_min=min(v['ret']for v in dres['DET'].values()if v['how']=='sieve');loss_max=max(max(v['loss_f'],v['loss_g'])for v in dres['DET'].values())
	ok11b=nb_ok and all(x[3]for x in rr)and ok_tab and ok_sample
	say(f"[11b]: table at {len(spot)} primes {K(ok_tab)}; "
	"special: "+"; ".join(f"{p}: {tab[str(p)]}"for p in special)
	+"; residue: "+", ".join(f"{k} {v[0]}/{v[1]}"for k,v in sample.items())+f" {K(ok_sample)}"
	+f"; primitive {RT.base_primitive():.4f}, ordinary {ordf:.4f}/{ordg:.4f}; omitted means {ex(pf):.2f}/{ex(pg):.2f} bits, "
	f"P(X > 96) {RT.cell_loss('f',96):.2e}/{RT.cell_loss('g',96):.2e}; loss <= {loss_max:.1e}, retention >= {floor(ret_min*1e4)/1e4:.4f}; bounds >= |N|: "+"; ".join(f"{b}{s_} R={R_}: {ok_}, gap {g:.1e} bits"for b,s_,R_,ok_,g,k in rr)
	+"; toys: "+"; ".join(str(tuple(x))for x in nb['rows'])
	+f" {K(nb_ok)}")
	mm,ok7e=W.encountered_check(seed)
	say(f"[7e] 8 encountered vs 8 sampled ideals: mean cf {mm[0][0]:.1f}/{mm[0][1]:.1f}, cg {mm[1][0]:.1f}/{mm[1][1]:.1f}; max {mm[0][2]:.1f}/{mm[1][2]:.1f} <= 97/14.5 {K(ok7e)}")
	toyecm=EC.toy(60,seed,sizes=((36,1024,50),(40,2048,50)))
	sch=[EC.schedule(69.0,int(ceil(dres['DET'][((69,75),'f')]['nf']))),EC.schedule(bt['nA'],int(ceil(bt['m'])))]
	ok12=all(t['empirical']>=t['model_corrected']-2*sqrt(t['model_corrected']*(1-t['model_corrected'])/t['trials'])for t in toyecm)
	say(f"[12] ECM: "+"; ".join(
	f"{x['fbits']:.2f}/{x['nbits']} bits: {x['curves']} curves, 2^{x['log2_instructions']:.2f}"for x in sch)
	+"; toys (bits, empirical, corrected): "
	+"; ".join(f"({t['bits']}, {t['empirical']:.3f}, {t['model_corrected']:.3f})"for t in toyecm)
	+f" {K(ok12)}")
	wt=M.toy(seed=seed);wn=M.window(cost,op['lb']);sh=[wn.ecm_share(r[k],op['slack_'+k[1]]+8)for k in('nf','ng')for r in anchors]
	se=lambda m:2*sqrt(m*(1-m)/wt['survivors']);ok12w=wt['measured']<=wt['model']+se(wt['model'])and wt['survivors']>=150
	say(f"[12w] anchors' ECM share f {min(sh[:8]):.4f}-{max(sh[:8]):.4f}; toy 2^{wt['n']}: {wt['survivors']} survivors, "
	f"ECM path {wt['measured']:.4f} <= {wt['model']:.4f} {K(ok12w)}")
	fl=check_fields()
	say(f"[13]: N(f(X)) deg {fl['deg_f']} and N(g(X)) deg {fl['deg_g']} squarefree {K(fl['unramified'])}; "
	f"disc h = {fl['disc_h']}; K_f {fl['sig']['f']}, K_f+ {fl['sig']['fp']}, K_g {fl['sig']['g']}: sigma-odd rank {fl['sig']['odd']}, chi-part {fl['sig']['sm']} = SM {M.N_SM} {K(fl['sig']['sm']==M.N_SM)}")
	ok13=fl['unramified']and fl['sig']['sm']==M.N_SM and fl['sig']['odd']==2*M.N_SM
	import runs16 as local2
	al_s={p:alpha_contrib(p)for p in(2,3,19)}
	ok_alpha=min(op['alpha_f'],DS.P.ALPHA_F)>=max(ALPHA_F,ALPHA_F_50K)and min(op['alpha_g'],DS.P.ALPHA_G)>=max(ALPHA_G,ALPHA_G_50K)and local2.facts()[0]
	say(f"[14] alpha: ({ALPHA_F:.6f}, {ALPHA_G:.6f}; 50k ({ALPHA_F_50K:.6f}, {ALPHA_G_50K:.6f})) <= ({op['alpha_f']}, {op['alpha_g']}) {K(ok_alpha)}; "
	f"in-run p=2 (exact, local2)/3/19: "+", ".join(f"{p}: ({v[0]:.6f}, {v[1]:.6f})"for p,v in sorted(al_s.items()))+f"")
	ev=local2.expectations();kp=(M.kappa2(840.29/69),M.kappa2(530.75/69));ok_v2=all(0<=M.V2E[i]-ev[i]<1e-6 and kp[i]<=M.KAP2[i]for i in(0,1));ok_toy=toy_v2(load_stats()['toy_v2'])
	say(f"[14v] v2: E {M.V2E} >= {[float(x)for x in ev]}; kappa2 (840, 530) {kp} <= {M.KAP2} {K(ok_v2)}")
	tot=M.total_log2(e['log2_relcol'],la['log2_la'],cap,la['log2_gen'])
	say(f"ledger: collection {e['log2_relcol']:.4f}, LA {la['log2_la']:.4f}, generator {la['log2_gen']:.2f}, "
	f"cap {cap:g} (expected {two_log:.3f}), canonical {M.LOG2_CANON:.0f}, SM {M.LOG2_SM:.0f}, "
	f"FB {M.LOG2_FBC:.0f}, misc {M.LOG2_MISC:.0f}, verify {M.LOG2_VERIFY:.0f}, band split {M.LOG2_SPLIT}")
	say(f"total_log2W_full_ledger = {tot:.12f} (claimed {CLAIM:.10f}: {'ok' if tot<=CLAIM else 'EXCEEDS CLAIM'}); "
	f"reserve {CLAIM-tot:.3f} bits")
	vals=dict(ledger=tot,relcol=e['log2_relcol'],la=la['log2_la'],gen=la['log2_gen'],log2_rels=e['log2_rels'],rels_2se=rels_lo,margin=e['margin'],
	work_se=work_se,two_log=two_log,markov=markov,log2_EN=log2_EN,over=log2(OVER),rep=log2(REP),
	boot_phi=bt['phi'],success=success,success_head=success_head)
	dec_ok={k:(vals[k]<=v if op=='<=' else vals[k]>=v)for k,(op,v)in DECLARED.items()}
	dec_ok['work_factor']=corr["safe_work_population_factor"]<=F(1032972622,10**9)and corr["safe_yield_factor"]>=F(966347073,10**9)
	ok_declared=all(dec_ok.values())
	say("declared: "+"; ".join(f"{k} {vals[k]!r} {op} {v} {K(dec_ok[k])}"for k,(op,v)in DECLARED.items()))
	checks=dict(owner_schedule=O.check(M,la),resultant=ok0,tower=ok0b,rho_semismooth=ok1,pool=ok2 and ok_size and ok_w and ok_r,box_bounds=ok2b,nbar=ok9,la_loop=k5_ok and bf_n==97 and bf_br==['b.ne']and b2_n<=M.C_PROJ and b2_br==['b.ne']and tail_ok and tail<=M.C_ROW and red9<=M.C_RED9 and ok_chi and nu_ok  and ok_wide and ok_build and ok_line and ok_oh,rec7=okr and ok_c7 and okb,band_split=ok_split,srt=ok_srt,tdiv=ok_tdiv,filt_y=ok_fy,hx=ok_hx,
	ball_counts=ok_counts and ok_relmeans,lazy=ok_lazy,rels_input=ok_rels_input,work_allowance=ok_allow,declared_bounds=ok_declared,alpha_ok=ok_alpha and ok_v2 and ok_toy,missing_logs=ok_rec,record_widths=ok_width,descent_plan_feasible=mu_min>=Pm.MIN_RELS and psi_max<=Pm.PSI_MAX and max(r['phi']for r in rows)<=0.05 and coll<=-4,
	audit=ok11,retention=ok11b,ecm=ok12 and ok12w,fields=ok13,feasible=e['feasible'],within_claim=tot<=CLAIM,success_ok=success>=0.39 and success_head>=0.39,udup_relw=ok11r,sampler_uniform=ok11s,encountered=ok7e)
	say(f"done in {time.time()-T0:.1f}s; checks: {checks}")
	print(json.dumps({"checks_ok":all(checks.values()),"checks":checks,"total_log2W_full_ledger":ceil(tot*10**12)/10**12,"claimed":CLAIM,"success_exact":success}),flush=True)
if __name__=="__main__":
	if "--alpha" in sys.argv:
		print(json.dumps(dict(zip(("alpha_f","alpha_g"),alpha_sweep(ALPHA_PMAX))),sort_keys=True));sys.exit(0)
	if "--alpha50k" in sys.argv:
		print(json.dumps(dict(zip(("alpha_f","alpha_g"),alpha_sweep(ALPHA_G_PMAX))),sort_keys=True));sys.exit(0)
	main()

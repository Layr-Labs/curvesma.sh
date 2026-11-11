"";import base64,hashlib,json,lzma,math,os,re,sys,time;from decimal import Decimal,localcontext;from fractions import Fraction as F;from math import log2,ceil,floor,sqrt;K=lambda b:'ok'if b else 'FAIL';from random import Random;_HERE=os.path.dirname(os.path.abspath(__file__));_DIR=os.path.join(_HERE,"files")if os.path.exists(os.path.join(_HERE,"files","model16.py"))else _HERE;sys.path.insert(0,_DIR);import common12 as C;import model16 as M;import owner10 as O;O.install(M);import descent16 as DS;import store16;import ideal_canonical as IC;import lattice_count_bounds as Xt;import filtsim;P=4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787;U=C.U;H=C.H;LN2=math.log(2);CLAIM=121.9088135851;import retain16;_P=2000;ALPHA_F,ALPHA_G=-.047244,.706684;_A=50000;_q=.745775;_r=.285135
def _m(p):
	if p==2:import runs16 as local2;return local2.alpha2()
	vf=vg=.0
	for d1,e,mf,mg in retain16.structure(p):
		q=float(p**d1)
		for side,mults in enumerate((mf,mg)):
			v=.0
			for mu in mults:v+=q/(q*q-1.)if mu==1 else 1./(q+1.)
			if side==0:vf+=d1*v
			else:vg+=d1*v
	t=math.log(p)/math.log(2.);return t*(1./(p-1.)-vf),t*(1./(p-1.)-vg)
def toy_v2(T,k=7):
	import itertools as I;V=[v for v in I.product((-1,0,1),repeat=6)if v>tuple(-x for x in v)];from threshold_model import primes_upto as pu;P=pu(1<<20);r=[0,0,0];ok=1;o=""
	while P[1:]:P=[math.prod(P[i:i+2])for i in range(0,len(P),2)]
	for j,(a,b)in enumerate(I.product(V,V)):
		if j%k==0 and M.coprime_pair(a,b):z=abs(C.norms(a,b)[0]);r[0]+=1;r[1]+=z%64==0;r[2]+=pow(P[0],256,z)==0
	for lb,x in sorted(T['lbs'].items()):c=x['cls'].values();t=[sum(d[j]for d in c)for j in(1,2,3)];q=lambda j:sum((d[1]-d[j])**2/d[j]for d in c);ok&=t[2]<=t[1]+t[0]**.5 and q(3)<q(2);o+=f" 2^{lb}: {[d[1:]for d in c]} chi2 {q(2):.0f}/{q(3):.0f};"
	Ya(f"[14t] toy census ({T['n']} pairs; v2 0/6):{o} replay {r} {K(ok and r==T['replay'])}");return ok and r==T['replay']
def _p(Zu):
	af=ag=.0
	for p in retain16.primes(2,Zu):a,b=_m(p);af+=a;ag+=b
	return af,ag
_W={'ledger':('<=',121.9088135851),'expected':('<=',113.229379668338),'relcol':('<=',114.157499000),'la':('<=',112.476927692),'gen':('<=',101.837),'log2_rels':('>=',63.1279),'rels_2se':('>=',63.0894),'margin':('>=',.1168),'work_se':('<=',.003321),'two_log':('<=',109.277),'p_cap':('<=',.0001586),'log2_EN':('<=',14.9427),'over':('<=',-12.35),'rep':('<=',-10.7),'boot_phi':('<=',7.2e-07),'success':('>=',.99904),'success_head':('>=',.999031)};Yz=[];T0=time.time();_k=C.stats
def Ya(*args):line=" ".join(str(a)for a in args);Yz.append(line);print(line,file=sys.stderr,flush=True)
def _U(A,B):return[c%P for c in C.pmul([c%P for c in A],[c%P for c in B])]
def _Y(A,n):
	r,base=[1],[c%P for c in A]
	while n:
		if n&1:r=_U(r,base)
		base=_U(base,base);n>>=1
	return(r+[0]*6)[:6]
def _o():
	F0=[0]*6;uw=[-U*x for x in C.WP]
	for d,c in C.PTERMS:
		term=(c+[0]*6)[:6]
		for _ in range(d//2):term=C.pmul(term,uw)
		F0=[x+y for x,y in zip(F0,term)]
	r=C.norm(C.elt(F0));ok=r!=0 and r%P==0;Ya(f"[0] Res: {r.bit_length()} bits, 0 mod p {K(ok)}");return ok
def _V(p):
	import retain16 as RT;hp=[x%p for x in H];rep=[]
	for g in RT.factor_Fp(hp,p):
		if g not in rep:rep.append(g)
	Ye=[1];sq=[]
	for g in rep:
		e,rest=0,hp
		while True:
			q,rm=RT.pdiv(rest,g,p)
			if rm:break
			rest,e=q,e+1
		for _ in range(e):Ye=[sum(Ye[i]*g[k-i]for i in range(len(Ye))if 0<=k-i<len(g))for k in range(len(Ye)+len(g)-1)]
		if e>1:sq.append(g)
	assert len(Ye)==7 and all((a-b)%p==0 for a,b in zip(H,Ye));Fd=RT._trim([(a-b)//p%p for a,b in zip(H,Ye)]);return bool(sq)and all(len(RT.pgcd(Fd,g,p))==1 for g in sq)
def _I():
	yp=_Y([0,1],P);cur,Yf=[0,1,0,0,0,0],True
	for k in range(1,7):
		acc=[0]*6
		for coef in reversed(cur):acc=_U(acc,yp);acc[0]=(acc[0]+coef)%P
		cur=acc;Zg=cur==[0,1,0,0,0,0]
		if k==3:y3=cur
		if(k in(1,2,3)and Zg)or k==6 and not Zg:Yf=False
	uw=[-U*x%P for x in C.WP];YC=_Y(uw,(P**6-1)//2)==[P-1,0,0,0,0,0];_8=y3==[x%P for x in C.SIGP]and _Y(uw,(P**3-1)//2)==[x*pow(2,-1,P)%P for x in C.C2P];import descent_witness as W;chi=pow(P,3,W.R);_T=(chi*chi+1)%W.R==0;tf=lambda t:[C.sigma(c)for c in t];hs=[0]*6;pw=[1,0,0,0,0,0]
	for k in range(7):hs=[x+H[k]*y for x,y in zip(hs,pw)];pw=C.pmul(pw,C.SIGP)
	cc=C.CU;c2=C.mul(cc,cc);_X=not any(hs)and C.mul(cc,C.sigma(cc))==[-1,0,0,0,0,0]and all(C.sigma(C.sigma(e))==e for e in([int(i==j)for j in range(6)]for i in range(6)));_X&=all(C.mul(C.sigma(c),C._pw(cc,d))==C.mul(C._pw(c2,6),c)for d,c in C.TERMS)and C.mul(C.sigma(C.UW),[1,0,0,0,0,0])==C.mul(c2,C.UW);import runs16 as local2;l2,_=local2.facts();d2=[[F(x,2)for x in r]for r in C.BAS2];det=F(IC.determinant([[2*x for x in r]for r in d2]),2**6);Zk=C.norm(C.elt([0,-4,-6,8,0,6]));Y4=Zk==2957312==16*184832 and abs(det)==F(1,4)and l2 and _V(19)
	Ya(f"[0b] h irred {K(Yf)}; (u w)^((p^6-1)/2) = -1 {K(YC)}; tau = Frob^3 {K(_8)}; "f"chi^2 = -1 {K(_T)}; sigma, c, tau(f), tau(g) {K(_X)}")
	Ya(f"[0c] |Disc(h)| = {Zk} = 4^2 x 184832; 2-maximal {K(l2)}, 19-maximal {K(_V(19))}: Disc(O_h) = 2^9 x 19^2 {K(Y4)}");return Yf and YC and Y4 and _8 and _T and _X
def _e():
	ref={2:.30685281944,3:.048608388,4:.0049109256,6:1.9649696e-5,8:3.2320693e-8,10:2.7701718e-11};ZG=max(abs(DS.rho_series(u)-v)/v for u,v in ref.items());YF=all(DS.rho_lo(u)<=DS.rho_series(u)<=DS.rho_hi(u)+1e-300 for u in(1.5,2.3,7.7,15.2,33.3,41.));rows,h=DS.window_rows(80.,69.,48,48);Yo=abs(sum(rows[1,0])-math.log(80/69))<1e-12
	Ya(f"[1] rho vs ref {ZG:.1e} {K(ZG<1e-6)}; brackets {K(YF)}; weights "f"{K(Yo)}")
	return ZG<1e-6 and YF and Yo
def canonical(cfg,a):
	try:return IC.canonical_a(a,cfg['matrix_a'],100000)[0]
	except RuntimeError:return IC.canonical_a(a,cfg['matrix_a'],10**7)[0]
def _B(cfg,pd):
	bad=[];rows=[];import runs16 as ZJ;units=ZJ.unit_list();qi,dq,_=C.integer_metric(cfg['matrix_a']);A2=cfg['production_A_squared']*dq
	for an in pd['anchors']:
		sd,i=an[0],an[1];r=M.pool_row(cfg,sd,i);a=C.regenerate_pair(cfg,sd,i)[0];la,th=M.reg(cfg,a);lim=int(th*A2);f=[round(r['qa']*65536),round(r['nf']*16),round(r['ng']*16),r['primitive']+2*r['rejected'],r['evals'],r['bcells'],r['leaves'],r['blen'],r['rlen']//16]+list(ZJ.run(a,units,qi,lim));lo,hi=f[-2:];r['lrs']=ZJ.lrs(a,lo,hi,qi,dq,cfg['production_A_squared'],la,F(cfg['production_region_kap']));r['T']=hi-lo+1;r['run']=a,lo,hi,lim
		if f!=an[2:13]or not canonical(cfg,a):bad.append(('anchor',sd,i))
		rows.append(r)
	_4=0
	for sd,(hi,idx)in sorted(pd['replay'].items()):got=[i for i in range(hi)if canonical(cfg,C.regenerate_pair(cfg,int(sd),i)[0])];_4+=len(got);bad+=[('replay',sd)]if got!=idx else[]
	Ya(f"[2] {len(pd['anchors'])} anchors, draws 0..{hi-1} of seeds {', '.join(sorted(pd['replay']))} "f"({_4} canonical) lists: {'ok' if not bad else 'FAIL '+str(bad[:3])}")
	return not bad,_4,rows,len(units)
def _Q(cfg,pool,l2,k=6,cells=8,seed=5):
	import embedding_bounds as EB;rng=Random(f"box-check:{seed}");qa,da,_=C.integer_metric(cfg['matrix_a']);qb,db,_=C.integer_metric(cfg['matrix_b']);A2=cfg['production_A_squared'];disks=EB.root_disks();n=bad=Zc=n2=0
	for r in rng.sample(pool,k)+[dict(seed=s,index=i)for s,i in l2]:
		a,b,aq,bq=C.regenerate_pair(cfg,r['seed'],r['index']);est=EB.NormBounds(a,disks);box=M.box_rows(qb,(M.reg(cfg,a)[0]*A2-F(cfg['production_region_kap'])*F(aq,da))*db,b);blo=b[0]//256*256;leaves=M.two_level(est,[(est.row(list(s)),l,h)for s,l,h,_ in box],[q for*_,q in box],blo);n2+=any(x[3]for x in leaves)
		for kk in rng.sample(range(len(box)),min(12,len(box))):
			s,l,h,_=box[kk];l,h=max(l,blo),min(h,blo+255);Zc+=l<=h
			for b0 in rng.sample(range(l,h+1),min(cells,max(0,h-l+1))):
				lf=[x for x in leaves if x[0]<=b0<=x[1]and(x[3]is None or kk in x[3])][0]
				if isinstance(lf[2],dict):continue
				for(ZA,ZC),v in zip(lf[2][:1],C.norms(a,[b0]+list(s))[:1]):n+=1;bad+=not(v!=0 and ZA-1e-6<=16*log2(abs(v))<=ZC+1e-6)
	return bad==0 and n2>=len(l2)>1,n,Zc
def _H(cfg,ct,Yh,th,m=1<<17):
	ok=ct['a'][0]<=1<=ct['a'][1]and M.COUNT[0]<=ct['a'][0]*ct['b'][0]and ct['a'][1]*ct['b'][1]<=M.COUNT[1]and Yh['safe_yield_factor']<=M.COUNT[0]and M.COUNT[1]<=Yh['safe_work_population_factor']and(lambda l:all(2*s*(l-th)>=1>=s*l for s in M.CLS))(F(cfg['production_region_lam']));o=[]
	for k,s in enumerate('ab'):
		q,d,_=C.integer_metric(cfg['matrix_'+s]);V,iq=C._frame(q);pts=ct['p'+s];cap=[math.isqrt(int(pts[-1][0]*iq[j][j]))for j in range(6)];g=Random(1900000000+300000000*k);h=[0]*len(pts);cr=Yh[('a_count_ratio','b_count_ratio_uniform')[k]];ok&=cr[0]<=ct[s][0]and ct[s][1]<=cr[1]
		for _ in range(m):
			y=[g.getrandbits(64)%(2*c+1)-c for c in cap];a=[sum(y[i]*V[i][j]for i in range(6))for j in range(6)]
			if a[5]>0:v=C.quadratic(q,a);h=[x+(v<=p[0])for x,p in zip(h,pts)]
		for(L,H,hp),hh in zip(pts,h):p=H/ct['n'];r=p*math.prod(2*c+1 for c in cap)/(math.pi**3/12*(L/d)**3)*sqrt(IC.determinant(q)/d**6);e=6*r/sqrt(H);ok&=hh==hp and abs(hh-m*p)<=4*sqrt(m*p)and ct[s][0]<=r-e and r+e<=ct[s][1];o.append(f"{r:.6f}+-{e/6:.6f}")
	return ok,", ".join(o)
def _g(Yb,p):
	import retain16 as RT;g=RT._trim([c%p for c in Yb])
	if len(g)<=1:return 0
	return len(RT.pgcd(g,RT.padd(RT.ppowmod([0,1],p,g,p),[0,p-1],p),p))-1
def nbar(Yu,Yw):
	from threshold_model import primes_upto;ps=[q for q in primes_upto(2**Yw)if q>=2**Yu];tot=tf=0;rng=Random(1)
	for p in ps:
		for s in sorted(set(_O(H,p,rng))):w=s+1;c=_g([w**6,0,w**5,0,0,0,2*w**3,0,0,0,-2*w,0,1],p);tf+=c;tot+=c+_g([U*w,0,1],p)
	return tot/len(ps),len(ps),tf/len(ps)
def Yc(name):t=lzma.decompress(base64.b85decode("".join(open(os.path.join(_DIR,"kernels_asm.b85")).read().split()))).decode();parts=re.split(r"(?m)^@@ ([A-Za-z][A-Za-z0-9_.-]*)\n",t);return dict(zip(parts[1::2],parts[2::2]))[name].rstrip("\n").splitlines()
def _G(ru,Yg):L=Yc("rec7.txt");T=M.REC7_TAB;KT=M.REC7_K;D={float(f[1]):[int(x)for x in f[2:]]for f in(ln.split()for ln in L if ln.startswith("K7 "))};nd=lambda v:(lambda m:m+2*(max(0,v[2]/1e6-v[0]*m*m)/(v[0]-1)/v[0])**.5)(v[1]/1e3/v[0]);ok=sorted(D)==sorted(e for e,_ in KT)and all(v[0]==512 and v[2]*512>=v[1]**2 for v in D.values())and all(a[0]>b[0]for a,b in zip(KT,KT[1:]))and KT[-1][0]==T[-1][0]and all(nd(D[e])<=k<=1 for e,k in KT)and all(a[0]>b[0]and a[1]<=b[1]for a,b in zip(T,T[1:]))and M.REC7_R2>=1;Zl=[ln for ln in L if ln.startswith(("CHECK ","WCHECK ","TCHECK "))];ok&=len(Zl)==3 and "ok (0 mismatches over 270 cosets"in Zl[0]and Zl[1].endswith("mismatches 0")and Zl[2].endswith("invalid 0: ok");pm=max(int(ln.split()[4])for ln in L if ln.startswith("PREP "));ideals=M.NBAR*(2**M.log2_li(64)-2**M.log2_li(M.P7));ok&=M.LOG2_PREP7>=log2(6*ideals*pm*2**.05/M.UNIT)and M.P7>=48 and(lambda L:"MATCH"in L[0]and max(int(f[6][4:])/int(f[3][8:])-(f[1]>"what=1")*M.C_APPLY7+(f[1]<"what=2")*(.0009-M.REC_ROUTE7)for f in L[1:5])<=0)([l.split()for l in Yc("lean8.txt")if l[0]!="#"])and"bit-20"in L[0];cl=[(log2(int(f[4])/(2*int(f[7])*int(f[5])*int(f[6]))),int(f[3])/int(f[4]))for f in(ln.split()for ln in L if ln.startswith("COSL "))];E6=[[int(i==j)for j in range(6)]for i in range(6)];ok&=len(cl)==6 and all(lh>=M.LINE_TAB[-1][0]and c<=M.tab(M.LINE_TAB,lh)for lh,c in cl)and C.P2F[0]==1 and C.norm(C.W)==2 and all(sum(x*y for x,y in zip(C.mul(C.W,e),C.P2F))%2==0 for e in E6);ok&=Yg==ru['units'];rn={int(p):F(int(b)-int(a),int(d)-int(c))for p,a,b,c,d in re.findall(r"record p=(\d+): Ir (\d+) -> (\d+), hits (\d+) -> (\d+)","\n".join(Yc("measurements12.txt")))};ok&=[float(rn[p])for p in sorted(rn)]==[x[1]for x in M.REC_TAB];return ok,f"{sum(v[0]for v in D.values())} production draws, mean + 2 SE of measured/reference <= K <= 1 in {len(D)} strata (largest need/K {max(nd(D[e])/k for e,k in KT if e in D):.4f}; largest draw {max(v[3]for v in D.values())/1e3:.3f}); {Zl[1][7:]}"
def _a(name,start):
	on,YV=False,[]
	for ln in Yc(name):
		s=ln.strip()
		if s.startswith(start):on=True;continue
		if on:
			if s.startswith(("LBB",".LBB"))or s.startswith(";")and "bb."in s:break
			if s and not s.startswith(";")and not s.startswith("."):YV.append(s.split()[0])
	return len(YV),[m for m in YV if m.startswith("b.")or m=="b"]
def _u(fn):
	L=[x.strip()for x in Yc("la_mac.s")];i=L.index(fn+":")+1;n,YY,ZE,calls=0,set(),True,[]
	while i<len(L)and not L[i].startswith(("_",";")):
		x=L[i];i+=1
		if x.startswith("LBB"):YY.add(x[:-1])
		elif x:
			n+=1;ZE=ZE and x.split()[-1]not in YY
			if x.split()[0]=="bl":calls.append(x.split()[-1])
	return n,ZE,calls
def Yv():
	c={};f="";ok=True
	for x in Yc("krylov5.s"):
		o=x.split()or[";"];g=f.split("/")[0]
		if x[-1:]==":":f=g+"/"*(x[-2]=="1")if x[0]=="."else x[:-1]
		elif o[0]!=";":
			c[f]=c.get(f,0)+1
			if f=="owner_part":ok&=o[0]in("stp","ldp","mov","add")or o[-1]in("unit_run5","gen_run1")
	return c["owner_part"]+2*c["unit_run5"]+10*c["gen_run1"],ok
def _F():f={k:_u(k)for k in("_send_pair","_combine","_row_red","_row_redc","_chi_mul")};d,Xr=Yv();ok=all(v[1]and not v[2]for v in f.values())and Xr;f["owner_part"]=d,;return O.row_tail(f["_send_pair"][0],f["_combine"][0],f["_row_red"][0],d),f["_row_redc"][0]+6,f["_chi_mul"][0]+4,{k:v[0]for k,v in f.items()},ok
def _n():
	import functools
	@functools.lru_cache(None)
	def Nd(d,s):
		if s<0:return 0
		if d==0:return 1
		r=math.isqrt(s);return sum(Nd(d-1,s-x*x)for x in range(-r,r+1))
	ok=True
	for R in 4,6,8:Xj=Nd(12,R*R);ok&=(Xj-1)//2==DS.CELLS[R]
	Z7=math.pi**6/720
	for R,half in DS.CELLS.items():n=2*half+1;ok&=Z7*(R-sqrt(3))**12<=n<=Z7*(R+sqrt(3))**12
	return ok
def _O(Yb,p,rng):import retain16 as RT;return RT.roots_Fp(Yb,p)if p<4096 else sorted(-f[0]%p for f in RT.factor_Fp(Yb,p)if len(f)==2)
def _w(dres):import counts16;EN,_=counts16.expected_nodes(DS.P(),dres);Xu=min(log2(2*(2**DS.log2li(b[1])-2**DS.log2li(b[0])))for b in dres['bins']);over=ceil((log2(EN)-M.NODE_CAP_LOG2)*100)/100;Yx=ceil((2*M.NODE_CAP_LOG2-1+6-Xu)*100)/100;return log2(EN),over,Yx,Xu
def audit(_f=31001,ZQ=None,Z6=None):
	import descent_witness as W;assert(W.P**6-1)%W.R!=0 and(W.P**12-1)//W.R%W.R!=0;rng=Random(_f)
	while True:
		x=[rng.randrange(W.P)for _ in range(6)],[rng.randrange(W.P)for _ in range(6)]
		if any(x[0]+x[1]):break
	w0=W.pow12(x,(W.P**12-1)//W.R);tt=rng.randrange(1,W.R);z=W.pow12(w0,tt);B,swaps=W.reduce(W.boot_basis(z),W.QB);Z5=W.boot_ball_stats(z,B,_f+1,points=6,radii=(96,));_v=max(x[1]for x in DS.P.BOOT_BINS);sq=[];Pm=DS.P();cfmax=max(b[1]for b in Pm.QBINS);cgmax=max(b[2]for b in Pm.QBINS);Qd=True
	for bits,side,seed in(260,'f',13300),:
		r=W.specialq_ball_stats(bits,side,seed,points=6,radii=(96,))
		for ZK in r['stats']:Qd&=ZK['cf_max']<=cfmax and ZK['cg_max']<=cgmax
		sq.append(dict(bits=bits,side=side,seed=seed,swaps=r['swaps'],stats=[{k:round(v,2)if isinstance(v,float)else v for k,v in s.items()}for s in r['stats']]))
	q_,_,_,_,nB,_=W.specialq_lattice(225,'g',Random(_f+31));n_,Z0,ZX,ZW,R_,_=DS.PLAN16[220,'g'];tr=W.orbit_triples(nB,_f+32,R_,12);Xk=log2(q_);sn=.0
	for Xm,Xl,m,_h in tr:pf=DS.ss_lower2(Xm,Pm.LX,float(n_),Z0,ZX,ZW,Pm.MESH,Pm.ALPHA_F,Pm.TRUNC);pg=DS.ss_lower2(Xl-Xk,Pm.LX,float(n_),Z0,ZX,ZW,Pm.MESH,Pm.ALPHA_G,Pm.TRUNC);sn+=pf*pg/m
	n1,Xo,Xn,_=DS.BOOT16;yb=lambda B,s:(lambda t:(sum(DS.ss_lower2(x[1],Pm.LX,float(n1),Xo,Xn,1,Pm.MESH,Pm.ALPHA_G,Pm.TRUNC)/x[2]for x in t),t))(W.orbit_triples(B,s,4,24,boot=True));sb,tb=yb(B,_f+33);ys=[sb]
	for i in range(7):z7=W.pow12(w0,rng.randrange(1,W.R));B7=W.reduce(W.boot_basis(z7),W.QB)[0];W.certify_boot(z7,B7);ys.append(yb(B7,_f+50+i)[0])
	Xd=dict(YR=max(t[2]for t in tr),Yl=max(t[3]for t in tr),YU=max(t[2]for t in tb),Ym=max(t[3]for t in tb))
	Zp=Random(_f+77);dp=[1]+[0]*16
	for _ in range(12):
		ZR=[0]*17
		for s in range(17):
			for v in range(-4,5):
				if s+v*v<=16:ZR[s+v*v]+=dp[s]
		dp=ZR
	tot=sum(dp[1:]);sh=[0]*17;_C=True
	for _ in range(600):
		c=W.int_ball_point(Zp,4);n2=sum(x*x for x in c);_C&=0<n2<=16
		for x in reversed(c):
			if x:_C&=x>0;break
		sh[n2]+=1
	chi2=sum((sh[s]-600*dp[s]/tot)**2/(600*dp[s]/tot)for s in range(1,17))
	return dict(boot=dict(swaps=swaps,stats=[{k:round(v,2)if isinstance(v,float)else v for k,v in s.items()}for s in Z5],_d=all(s['cg_max']<=_v for s in Z5)),YX=sq,_b=Qd,Yk=dict(_i=sn/12/2.**ZQ,_j=sb/24/2.**Z6,pair=8*sum(y*y for y in ys)/sum(ys)**2,**Xd),sampler=dict(chi2=chi2,ok=bool(_C)and chi2<60),storage=store16.check())
def _z():
	import retain16 as RT,descent_witness as W;R=W.R;Hr=[x%R for x in H]
	def res(a,b):
		a=RT._trim([x%R for x in a]);b=RT._trim([x%R for x in b]);s=1
		while b and len(b)>1:
			da,db=len(a)-1,len(b)-1;rem=RT.pmod(a,b,R)
			if not rem:return 0
			s=s*(-1)**(da*db)*pow(b[-1],da-(len(rem)-1),R)%R;a,b=b,rem
		return s*pow(b[0],len(a)-1,R)%R if b else 0
	def _R(terms,deg):
		n=6*deg+1;xs=list(range(n));ys=[res(Hr,[sum(c[i]*pow(k,d,R)for d,c in terms if i<len(c))%R for i in range(6)])for k in xs]
		for j in range(1,n):
			for i in range(n-1,j-1,-1):ys[i]=(ys[i]-ys[i-1])*pow(xs[i]-xs[i-j],-1,R)%R
		Yb=[ys[-1]]
		for i in range(n-2,-1,-1):Yb=[(-xs[i]*Yb[0]+ys[i])%R]+[(Yb[k-1]-xs[i]*Yb[k])%R for k in range(1,len(Yb))]+[Yb[-1]]
		return Yb
	out={}
	for name,terms,deg in('f',C.PTERMS,12),('g',RT.G_TERMS,2):Yb=_R(terms,deg);out['deg_'+name]=len(Yb)-1;out['sqf_'+name]=Yb[-1]==1 and len(Yb)==6*deg+1 and len(RT.pgcd(Yb,RT.pderiv(Yb,R),R))==1
	out['unramified']=out['sqf_f']and out['sqf_g'];out['sig']=signatures();return out
def signatures():
	def tr(p):
		while p and p[-1]==0:p=p[:-1]
		return p
	def YA(p):
		s=[tr([F(c)for c in p])];s.append([i*c for i,c in enumerate(s[0])][1:])
		while len(s[-1])>1:
			a,b=s[-2][:],s[-1]
			while len(a)>=len(b):q=a[-1]/b[-1];d=len(a)-len(b);a=tr([x-q*b[i-d]if i>=d else x for i,x in enumerate(a)][:-1])
			if not a:break
			s.append([-x for x in a])
		return s
	def ch(s,x):v=[t for t in(sum(c*x**i for i,c in enumerate(q))for q in s)if t];return sum((a>0)!=(b>0)for a,b in zip(v,v[1:]))
	Zz=lambda p,B=F(10)**60:(lambda s:ch(s,-B)-ch(s,B))(YA(p));sh=YA(H);ys=[]
	def ZB(a,b):
		n=ch(sh,a)-ch(sh,b)
		if n==1 and b-a<F(1,10**30):ys.append((a+b)/2)
		elif n:ZB(a,(a+b)/2);ZB((a+b)/2,b)
	ZB(F(-10),F(10));rf=Y0=rg=0
	for y in ys:
		fx=[0]*13
		for d,c in C.PTERMS:fx[d]=sum(F(ci)*y**i for i,ci in enumerate(c))
		rf+=Zz(fx);Y0+=Zz(fx[0::2]);rg+=Zz([U*(y+1),0,1])
	rk=lambda n,r1:r1+(n-r1)//2-1;return dict(f=(rf,36-rf//2),fp=(Y0,18-Y0//2),g=(rg,6-rg//2),odd=rk(72,rf)-rk(36,Y0)+rk(12,rg)-rk(6,len(ys)),sm=(rk(72,rf)-rk(36,Y0)+rk(12,rg)-rk(6,len(ys)))//2)
def YS(seed,lu=10,m=24,n=12):
	import unit_entries as UE,descent_witness as W;P=W.R;g=Random(seed);c,_,R,_=filtsim.build_system(69.,(UE.MUF/69,UE.MUG/69),M.FILT_EXCESS,seed,lu);UE.MID.clear();filtsim.merge(filtsim.purge(R,len(c),M.N_SM,False)[0],len(c),(int(M.W_ROW)-1,),kmax=64,verbose=False,hook=UE);i4=pow(W.P,3,P);R=[{k:(int(v.real)+int(v.imag)*i4)%P for k,v in r.items()if k>=0}for r in UE.LAST];ix={k:j for j,k in enumerate(sorted({k for r in R for k in r}))};D=len(R);N=len(ix);B=[[(ix[k],v)for k,v in r.items()]for r in R];v=[g.randrange(1,P)for _ in R];S=[[g.randrange(P)for _ in R]for _ in range(D-N-1)];S.append([-(sum(x*v[j]for j,x in r)+sum(s[i]*v[N+t]for t,s in enumerate(S)))*pow(v[-1],-1,P)%P for i,r in enumerate(B)]);mv=lambda x:[sum(a*x[j]for j,a in r)%P for r in B];cb=lambda V,f:[sum(a*b for a,b in zip(c,f))%P for c in zip(*V)];Z=[[g.randrange(P)for _ in R]for _ in range(n-D+N)];Y=S+[mv(z)for z in Z];co=list(range(N));g.shuffle(co);bl=[0]*N+list(range(D-N))
	for t,j in enumerate(co):bl[j]=t*m//N
	X=[g.randrange(1<<128)for _ in R];d=-(-D//n)+1;e=-(-D//m)+1;a=[];V=Y
	for i in range(d+e):
		q=[[0]*n for _ in range(m)]
		for j,y in enumerate(V):
			for b,x,z in zip(bl,X,y):q[b][j]+=x*z
		a.append(q);V=[mv(y)for y in V]
	H=[[a[i+k][l][j]%P for k in range(d+1)for j in range(n)]for i in range(e)for l in range(m)];H=[h[1:]+[-h[0]]for h in H];nc=n*d+n-1;pv=[]
	for c in range(nc):
		t=len(pv);q=next((i for i in range(t,len(H))if H[i][c]),-1)
		if q<0:continue
		H[t],H[q]=H[q],H[t];iv=pow(H[t][c],-1,P);p=H[t]=[x*iv%P for x in H[t]]
		for i in range(t+1,len(H)):
			if H[i][c]:f=H[i][c];H[i]=[(x-f*z)%P for x,z in zip(H[i],p)]
		pv.append(c)
	y=[0]*nc
	for t in range(len(pv)-1,-1,-1):y[pv[t]]=(H[t][nc]-sum(H[t][j]*y[j]for j in pv[t+1:]))%P
	f=[1]+y;w=cb(Y,f[-n:])
	for k in range(d-1,0,-1):w=[(p+q)%P for p,q in zip(mv(w),cb(Y,f[k*n:k*n+n]))]
	w=[(p+q)%P for p,q in zip(w,cb(Z,f[D-N:n]))][:N]+f[:D-N];l=w[0]*pow(v[0],-1,P);return D,len(pv),not any(h[nc]for h in H[len(pv):])and all((p-l*q)%P==0 for p,q in zip(w,v))and not any((sum(x*w[j]for j,x in r)+sum(s[i]*w[N+t]for t,s in enumerate(S)))%P for i,r in enumerate(B))
def main():
	YT=json.load(sys.stdin)if not sys.stdin.isatty()else{};seed=int(YT.get("seed",1));Ya(f"seed={seed}; unit {M.UNIT:.0f} instructions");Xh=_o();Zy=_I();Xg=_e();cfg,pd=(lambda s:(s["config"],s["pool_d"]))(_k());Xf,_4,anchors,Yg=_B(cfg,pd);Y6,nb2,Xi=_Q(cfg,anchors,pd['l2'],seed=seed);Ya(f"[2b] 6 + {len(pd['l2'])} level-2 blocks ({Xi} rows): {nb2} f norms in leaf bounds {K(Y6)}");_s=F(cfg['production_theta']);_c=F(cfg['production_theta_region']);assert _c<=_s==M.OP['theta']and F(pd['theta_region'])==_c;ND=pd['draws_per_seed'];rmin=math.isqrt(int(_c*M.CLS[0]*cfg['production_A_squared']));Yh=Xt.correction(cfg,rmin);geo=M.Geometry(cfg,_s)
	def tot(seeds):n=sum(pd['draws'][x]for x in seeds);return M.pool_totals(M.add_sums([pd['sums'][x]for x in seeds]),n,n,geo)
	tD,tW,tB=tot(pd['judged']),tot(pd['work']),tot(pd['sizing']);Z9=_c-F(1,8000);nB=sum(pd['draws'][x]for x in pd['sizing']);Zo=M.pool_totals(pd['sizing_prev'],nB,nB,geo);op=M.OP;retain,cost=M.retentions(op)
	ws=[.0 if M.skip(r)else M.evaluate([r],1.,geo,op,retain,cost)['sums']['work']/geo.weight(r['qa'])for r in anchors];Yo=all(math.isclose(w,an[13],rel_tol=1e-7)for w,an in zip(ws,pd['anchors']))and abs(retain.at(740,op['slack_f'],op['slack_b'])['retention']-0.854531)<1e-6;qq,dd,_=C.integer_metric(cfg['matrix_a']);ex=[]
	for r in sorted(anchors,key=lambda r:r['T'])[:3]:a,lo,hi,lim=r['run'];can=[C.quadratic(qq,x)<=lim and x[5]>=1 and canonical(cfg,x)for x in([a[0]+t]+a[1:]for t in range(lo-1,hi+2))];ex.append((r['T'],all(can[1:-1]),can[0],can[-1]))
	Y3=all(e[1]for e in ex)
	Ya(f"[2w] anchors' work {K(Yo)}; [2r] 3 anchor runs "f"(T, canonical, ends) {ex} {K(Y3)}")
	ZO,c3=_H(cfg,pd['count'],Yh,_c);Ya(f"[3] theta={cfg['production_theta_region']}; counts {c3}; {M.COUNT} in certified ({float(Yh['safe_yield_factor']):.6f}, {float(Yh['safe_work_population_factor']):.6f}), R >= {rmin} {K(ZO)}");log2_U=M.log2_fb(op['lb']);need=log2_U+log2(op['excess'])
	e=dict(tD,log2_U=log2_U,Xb=need,log2_rels_est=tD['log2_rels'],log2_rels=tD['log2_rels']-M.REL_BAND_LOSS,margin=tD['log2_rels']-M.REL_BAND_LOSS-need,log2_relcol=tW['log2_est']+log2(M.CANON_WORK_ALLOWANCE),per_cell=tW['per_cell']*M.CANON_WORK_ALLOWANCE,log2_population=tW['log2_population'])
	e['feasible']=e['margin']>=max(M.MARGIN,M.MARGIN_SE*tD['yield_se_log2']);rule=lambda t:M.MARGIN+sqrt(2)*t['yield_se_log2'];_2=tB['log2_rels']-need>=rule(tB)and Zo['log2_rels']-need<rule(Zo)
	Ya(f"[4s] margin {tB['log2_rels']-need:.6f} >= 0.06 + sqrt(2) SE = {rule(tB):.6f} at {cfg['production_theta_region']}; "f"{Zo['log2_rels']-need:.6f} < {rule(Zo):.6f} at {Z9*8000}/8000: smallest {K(_2)}")
	work_se=tW['work_se'];_t=1+M.WORK_SE_K*work_se;_N=M.CANON_WORK_ALLOWANCE>=_t
	Ya(f"[4a] work {pd['work']}: row SE {100*tW['cond_se']:.4f}%, total SE {100*work_se:.4f}%; "f"allowance {M.CANON_WORK_ALLOWANCE:.5f} >= 1 + {M.WORK_SE_K:.2f} SE = {_t:.6f} {K(_N)}")
	_D=e['nf_rel']>=69.5*DS.P.MAIN_UF and e['ng_rel']>=69.5*DS.P.MAIN_UG;Yd=e['log2_rels']+log2(1-2*e['yield_se_log2']*math.log(2));Zd=e['margin']-log2(M.FILT_EXCESS/op['excess']);_l=Yd>=DS.P.LOG2_RELS and log2_U+log2(M.FILT_EXCESS)<=Yd and Zd>=max(M.MARGIN,M.MARGIN_SE*e['yield_se_log2']);Ya(f"[4] judged {pd['judged']}: rels 2^{e['log2_rels']:.5f} vs need 2^{e['Xb']:.5f}, margin {floor(e['margin']*10000)/10000:.4f} ({K(e['feasible'])}); SE {e['yield_se_log2']:.5f}; 2-SE count {floor(Yd*1e5)/1e5} (margin {floor(Zd*1e4)/1e4}) {K(_l)}");la=M.linear_algebra(e['log2_U'],op['excess'])
	Ya(f"[5] 2^{e['log2_population']:.3f} cells at {e['per_cell']:.1f}/cell: collection 2^{e['log2_relcol']:.4f} units; "f"N(W)/U {M.NU} -> N = 2^{la['log2_N']:.4f}; LA 2^{la['log2_la']:.4f} units (k5 {M.C_ENTRY_K})")
	from threshold_model import primes_upto;ZP=len(primes_upto(2**20));ZV,ZU=[-(-b//64)for b in M.BAND_NBITS];Y2=5*(ZV+ZU)+2*22+4;tdiv=ZP*Y2+12*(ZV*M.BAND_NBITS[0]+ZU*M.BAND_NBITS[1]);nt="\n".join(Yc("norms.txt"));nk=max(map(int,re.findall(r"max (\d+)",nt)));_1=ZP==82025 and tdiv<=M.C_TDIV and nt.count("mismatches 0 L-violations 0")==2 and nk<=M.C_NORMS and max(geo.rb)<2**13;Ya(f"[5t] tdiv 2^{log2(tdiv):.3f} <= C_TDIV 2^{log2(M.C_TDIV):g}; norms {nk} {K(_1)}")
	Us,res=filtsim.run(op['lb'],e['nf_mean']/op['lb'],e['ng_mean']/op['lb'],M.FILT_EXCESS,seed,12,kmax=64,targets=(100,),verbose=False)
	YO=res[100];Ya(f"[6] (2^12, seed {seed}): N/U = {1/YO[0]:.4f} at weight {YO[3]:.0f} (max coef {YO[2]})");rec=filtsim.missing_recon(seed,excess=M.FILT_EXCESS);Qb=all(mf<=md for(lo,hi,n,mf,md,h3)in rec)
	Ya(f"[6b] bands (n, measured, model, lost >= 3): "+"; ".join(f"[{lo:g},{hi:g}): {n}, {mf:.4f}, {md:.4f}, {h3}"for(lo,hi,n,mf,md,h3)in rec)+f" {K(Qb)}")
	_y=_n();import retain16 as Xs;Xs.load_table(_k());Pm=DS.P();dres=DS.evaluate(Pm,DS.PLAN16,DS.BOOT16);two_log=dres['two_log_units'];bt=dres['boot'];lz=sorted((v for v in dres['DET'].values()if 'lazy'in v),key=lambda v:-v['survcost_log2'])[:3];Yt=[(round(v['mu'],1),v['nb_r'],round(v['psi'],3),round(v['lazy'],4))for v in lz];_3=bt['own_log2']<bt['sub_log2'];Ya(f"[7l] {Yt}; lazy boot: {log2((1-bt['psi']**Pm.BOOT_RELS)/(1-bt['psi'])):.3f} bits vs log2 {Pm.BOOT_RELS} {K(_3)}");cap=M.DESCENT_CAP;pcap=2**(two_log-cap);log2_EN,over,Yx,Xu=_w(dres);Zn=2.**over;Z8=2.**Yx;_7=1-2*bt['phi']-pcap-Zn-Z8;_6=1.-2.**(la['log2_N']+2-128)-M.P_CAP;success=_6*_7;success_head=_6*(_7+pcap-1.1*pcap);rows=DS.report(dres);sv=[dres['DET'][k]for k in dres['DET']if dres['DET'][k]['how']=='sieve'];Zh=max(DS.CELLS[d['R']]for d in sv);_J=Zh<1<<103 and Zh-1>>40<1<<63 and all(16*max(d['nf'],d['ng'])<1<<16 for d in sv);mu_min=min(r['mu']for r in rows);psi_max=max(max(r['psi']for r in rows),bt['psi'])
	Ya(f"[7] balls {K(_y)}; {len(rows)} nodes; one log 2^{bt['one_log_log2']:.3f}; two logs 2^{two_log:.3f} units; "f"cap 2^{cap:g}: P(cap) <= {ceil(pcap*1e4)/1e4:.4f}; success {floor(success*10000)/10000:.4f}, {floor(success_head*10000)/10000:.4f} at +10%; widths {K(_J)}")
	nb,np_,ZS=nbar(14,15);Ya(f"[9] {np_} primes in [2^14, 2^15): {nb:.3f}, f only {ZS:.3f} (NBAR {M.NBAR}, sieved {M.NBAR*M.SIEVE_SIDES})");Xe=nb<=M.NBAR*1.01 and ZS<=M.NBAR*M.SIEVE_SIDES*1.01;tail,Y1,chi,YN,Yi=_F();_T=chi<=M.C_CHI;Ya(f"[10] tail: {YN} ({Yi}); {tail} <= C_ROW {M.C_ROW:.0f}; reduce {Y1} <= C_RED9 {M.C_RED9:.0f}; chi x {chi} <= C_CHI {M.C_CHI:.0f} {K(_T)}")
	def Xq(T,k):x=[v[1]/2**k for v in T[k]];n=len(x);m=sum(x)/n;return m,(sum((y-m)**2 for y in x)/(n-1))**.5/n**.5,n
	fs={k:Xq(M.FILT_RUNS,k)for k in M.FILT_RUNS};fx={k:Xq(M.FILT_X,k)for k in M.FILT_X};Xp=sorted(fx);_5=fx[Xp[-1]][0]+2*fx[Xp[-1]][1];ZL=max((sum((v[1]/2**k-sum(w[1]for w in vs)/len(vs)/2**k)**2 for v in vs)/(len(vs)-1))**.5*2**(k/2)for T in(M.FILT_RUNS,M.FILT_X)for k,vs in T.items()if len(vs)>=6);lu=la['log2_N']-log2(M.NU);Z4=(M.SD_C*2**(-lu/2)/(M.NU-M.NU_H))**2;Qe=ZL<=M.SD_C and Z4+1e-15<=M.P_CAP;Ya(f"[10c] cap {M.NU} U: sd {ZL:.3f} <= {M.SD_C}; {Z4:.2e} <= {M.P_CAP}: {Qe}");mc=max(v[2]for vs in(*M.FILT_RUNS.values(),*M.FILT_X.values())for v in vs);YH=Qe and max(fs[19][0],fs[20][0])<=fs[15][0]and mc<=255
	Ya(f"[10] N(400)/U at excess 1.05, means 744/501 (mean/SE/seeds): "+", ".join(f"2^{k} {m:.5f}/{s_:.5f}/{n}"for k,(m,s_,n)in sorted(fs.items()))+f"; N({M.W_ROW-1:.0f})/U at excess {M.FILT_EXCESS}, means {M.FILT_MEANS}: "+", ".join(f"2^{k} {m:.6f}/{s_:.6f}/{n}"for k,(m,s_,n)in sorted(fx.items()))+f"; 2^21 mean + 2 SE {_5:.6f}; max coef {mc} <= 255: {YH}")
	Xw=M.FILT_W.get(M.W_ROW-1,[(0,0)])[0];YD=all(fx[a][0]>fx[b][0]for a,b in zip(Xp,Xp[1:]))and e['nf_rel']<=M.FILT_MEANS[0]and e['ng_rel']<=M.FILT_MEANS[1]and _5<=M.NU_H<M.NU and M.W_ROW-1<=Xw[0]<=Xw[1]<=M.W_ROW;Ya(f"[10y] falls at every step 2^{Xp[0]}..2^{Xp[-1]}; weights {Xw}; relation means {floor(e['nf_rel']*1e3)/1e3:.3f}/{floor(e['ng_rel']*1e3)/1e3:.3f} <= {M.FILT_MEANS}; {_5:.6f} <= {M.NU_H} {K(YD)}");Z2,Z1,Yq=YS(seed);Yq&=Z1==Z2-1;Ya(f"[10x] RHS BW toy: rank {Z1} = {Z2} - 1, w = planted logs x scalar {K(Yq)}");Qc=M.generator(la['log2_N'],1)<=2**la['log2_gen']*M.UNIT;Ya(f"[10m] generator sum <= model {K(Qc)}");Yy,YK=_a("bw_dot.s","LBB1_2:");Ya(f"[10] bw_full: {Yy} (charged {M.C_FULL:.0f}) {K(Yy==97 and Yy<=M.C_FULL and YK==['b.ne'])}");wl="\n".join(Yc("wide12.txt"));ml="\n".join(Yc("measurements12x.txt"));fr=lambda k:{int(p):float(v)for p,v in re.findall(k+r" p=(\d+): Ir \d+ -> \d+, reps \d+ -> \d+, Ir per rep ([\d.]+)",wl)};ws,wr=fr("wsetup6w"),fr("wrec");ZT=wr.get(2**63-25,0);big=[v for p,v in wr.items()if p>>63];pm=int(re.findall(r"Ir per pair min \d+ mean \d+ max (\d+)",wl)[0]);Zs={int(p):F(int(b)-int(a),int(d)-int(c))for p,a,b,c,d in re.findall(r"TAB rec (\d+) (\d+) (\d+) (\d+) (\d+)",ml)}
	_0=len(ws)==2 and max(ws.values())<=406.87<=M.C_SETUP6*M.WIDE_SETUP and len(big)==3 and 2.5*max(abs(v-ZT)for v in big)<=(M.WIDE_REC-1)*ZT and wl.count("MATCH")==6 and wl.count("certified 200")==3 and "0 mismatches"in wl and "failures 0"in wl and "BYTE-IDENTICAL"in wl and M.LOG2_WIDE_LAT>=log2(6*1.25506*2**64/(64*math.log(2))*4*pm/M.UNIT)and[float(Zs[p])for p in(257,65537,8388617,67108879,268435459)]==[x[1]for x in M.REC_TAB_W]
	Ya(f"[10w] wide path {K(_0)}");ol="\n".join(Yc("measurements-oh.txt"));br=[(int(e),int(c),int(sg),int(a.replace(',','')),float(pa),int(z.replace(',','')))for e,c,sg,a,pa,z in re.findall(r"entries=(\d+) rows=\d+ what=build cells=(\d+) K=\d+ L4=\d+ segments=(\d+).*\n *([\d,]+) \(([\d.]+)%\)  \?\?\?:build_layout'2.*\n *([\d,]+) \(",ol)];_L=len(br)==4 and all(a/e<=M.C_BUILD_E and a/(pa/100)-a-z<2**20 and(a/(pa/100)-a-z)/sg<=M.C_BUILD_SEG for e,c,sg,a,pa,z in br);Ya(f"[10o] layout {K(_L)}");Lx=Yc("measurements12x.txt");Zt={};Z3={}
	for ln in Lx:
		f=ln.split()
		if ln.startswith("TAB line "):Zt[int(f[2])]=(int(f[4])-int(f[3]))/(int(f[6])-int(f[5]))
		if ln.startswith("MP "):d=Z3.setdefault((int(f[1]),int(f[2])),{}).setdefault(int(f[3]),[0,0,set(),0,0]);d[0]+=int(f[6])-int(f[5]);d[1]+=int(f[8])-int(f[7]);a_,b_=map(int,f[4].split('-'));d[2]|=set(range(a_,b_+1));d[3]=int(f[10]);d[4]=int(f[11])*(int(f[13])-int(f[12]))
	_S=len(Zt)==4 and all(c>=Zt[p]for(h,c),p in zip(M.LINE_TAB[1:],(257,4099,65537,1048583)))
	for(lp,pt),YZ in sorted(Z3.items()):hs=[log2(d[1]/(2*d[4]*d[3]))for d in YZ.values()];cs=[d[0]/d[1]for d in YZ.values()];n_=len(cs);mc=sum(cs)/n_;ub=sum(d[0]for d in YZ.values())/sum(d[1]for d in YZ.values())+M.LINE_SE_K*sqrt(sum((c-mc)**2 for c in cs)/(n_-1)/n_)if n_>1 else float('inf');_S&=all(d[2]==set(range(1,d[3]+1))for d in YZ.values())and n_>=3 and ub<=M.tab(M.LINE_TAB,max(hs))
	sp=[ln.split()for ln in Lx if ln.startswith("SP ")];_S&=len(sp)>0 and all(M.tab(M.LINE_TAB,log2(int(f[8])/int(f[9])/int(f[2])))>=(int(f[5])-int(f[4]))/(int(f[7])-int(f[6]))for f in sp);Ya(f"[10l] line {K(_S)}");ZM,Zr=_G(pd['runs'],Yg);import runs16 as Zm;qq,dd,_=C.integer_metric(cfg['matrix_a']);Ya(f"[10r] rec7: {Zr} {K(ZM)}");kt="\n".join(Yc("fft.txt"));mm=[(int(a),float(c))for ln in kt.splitlines()if ln.startswith("FFTPTS ")for a,c in(x.split(":")for x in ln.split()[1:])];bc=M.gbatch_cost(pbits=1.000028*2**48/M.LN2,leaf=13,L=39,c_leaf=2**19);gc=M.gbatch_cost();sr=[x.split(":")for ln in kt.splitlines()if ln.startswith("SRTWAVES ")for x in ln.split()[1:]];ZN=len(mm)>=25 and max(n for n,_ in mm)>=2**26 and all(c<=M.fft_F(n)for n,c in mm)and len(sr)>=10 and all(x[2]=="0"for x in sr)and bc<=M.C_BATCH and gc<=M.C_GBATCH and M.gbatch_cost(pbits=1.000028*2**64/M.LN2,leaf=13,L=54,c_leaf=2**19)<=M.C_BATCH2;Ya(f"[10b] FFT: {len(mm)} products, charged/F <= {max(c/M.fft_F(n)for n,c in mm):.4f}; {len(sr)} SRT waves; band 2^{log2(bc):.4f} <= C_BATCH 2^{log2(M.C_BATCH):g}; [10g] g batch 2^{log2(gc):.4f} <= C_GBATCH 2^{log2(M.C_GBATCH):g} {K(ZN)}");sb,sf=M.srt_toy(8,seed);Qa=sb==0 and sf<=2;Ya(f"[10v] SRT: 2^8 leaves, {sb} mismatches, {sf} corrections {K(Qa)}");Xc=M.rmin(retain,op['slack_f'],op['slack_b']);sp=M.band_split(e['log2_rels_est'],e['yield_se_log2'],Xc);q=pd['sp'];n_=sum(q);ss,sg,ls,lg=[retain.co(x)for x in('sieve_sf','sieve_geo','large_sf','large_geo')];E=0;o4=""
	for nb in 700,740,780:
		n=4*nb;rb=(retain.cap(nb)+1)/16;a=floor((op['slack_f']-rb)*4);b=floor((op['slack_b']-rb-.125)*4);al=[];be=[]
		for j in range(len(q)):s=8*j+8;w=math.exp(retain.tilt*s/4);al.append(w*sum(ls[l]*ss[n-s-l]for l in range(min(b,a-s)+1)));be.append(w*sum(lg[l]*sg[n-s-l]for l in range(n-s+1)))
		B=sum(k*y for k,y in zip(q,be));r=sum(k*x for k,x in zip(q,al))/B;se=sqrt(sum(k*(x-r*y)**2 for k,x,y in zip(q,al,be)))/B;sl=retain.selection_normalizer;rl=retain.at(nb,op['slack_f'],op['slack_b'])['retention'];E=max(E,-log2(sl*(r-2*se)/rl));o4+=f" {nb}: {sl*r:.4f} +- {sl*se:.4f} vs {rl:.4f};"
	Xw=n_>=8000 and E<=M.REL_BAND_LOSS;Ya(f"[4b] r_f with the odd parts < 2^{M.LX0:g} of {n_} norms:{o4} loss at - 2 SE <= {ceil(E*1e4)/1e4} <= {M.REL_BAND_LOSS} {K(Xw)}");Xv=[[int(x)for x in f[2:]]for f in map(str.split,Yc("sb.txt"))if f[0]=="SB"];g2=[(lambda m:m+2*(max(0,v[2]-v[0]*m*m)/(v[0]-1)/v[0])**.5)(v[1]/v[0])for v in Xv];Xx=Yc("sb.txt")[0].count("MATCH")==3 and[v[0]for v in Xv]==[3,1024]and M.C_SB>=sum(g2)+200;Ya(f"[10k] stage B: K1 {ceil(g2[0])} + K2 {ceil(g2[1])} + 200 <= C_SB {M.C_SB:.0f} {K(Xx)}");rj=json.load(open(os.path.join(_DIR,'refinement.json')));_K=sp['log2_units']<=M.LOG2_SPLIT==rj['ledger_upper_rows'][-1]and rj['split_inputs'][0]>=e['log2_rels_est']and rj['split_inputs'][1]>=e['yield_se_log2']
	Ya(f"[10t] band split: 2^{sp['c_rel_log2']:.4f} per relation; "f"relations <= 2^{sp['log2_count']:.4f}: 2^{sp['log2_units']:.4f} <= 2^{M.LOG2_SPLIT} {K(_K)}")
	D7=M.contam_sup(geo);g7=.0
	for i in range(4002):
		S7=geo.slice_cells(float(_c*M.CLS[i&1])*(i>>1)/2000,i&1);c7=.0;lp=M.P7
		while lp<M.X7:lp2=min(lp+.25,M.X7);lh=log2(S7)-.5*(lp+lp2);dp=M.NBAR/2*(2**M.log2_li(lp2)-2**M.log2_li(lp));wd=lp>=M.WIDE_FROM;c7+=(S7*M.NBAR*math.log(lp2/lp)*(M.WIDE_REC*M.tab(M.REC_TAB_W,lh)if wd else M.tab(M.REC_TAB,lh))if lh>=M.REC_TAB[-1][0]else dp*M.C_COSET_EMPTY)+dp*M.C_SETUP6;lp=lp2
		g7=max(g7,c7)
	ct=pd['contam'];ZF=[]
	for sd,i,T_ in ct['replay']:a=C.regenerate_pair(cfg,sd,i)[0];lim=int(M.reg(cfg,a)[1]*cfg['production_A_squared']*dd);lo,hi=Zm.run(a,Zm.unit_list(),qq,lim);ZF.append(hi-lo+1==T_ and all(C.quadratic(qq,x)<=lim and x[5]>=1 and canonical(cfg,x)for x in([a[0]+t]+a[1:]for t in range(lo,hi+1))))
	cr=float((F(ct['theta'])/_c)**3);YE=D7>=g7 and ct['bad']==0 and M.RUN_ALLOW>=(1-.05**(1/ct['n']))*cr and all(ZF)
	Ya(f"[k7c] D_SUP 2^{log2(D7):.3f} >= 2^{log2(g7):.3f}; {ct['n']} runs ({ct['selection_sha']}): "f"{ct['points']} a, {ct['bad']} bad, RUN_ALLOW {M.RUN_ALLOW} >= {1-0.05**(1/ct['n']):.6f} x {cr:.4f}; replay {ZF} {K(YE)}")
	oh=[ln.split()for ln in ol.splitlines()if ln.startswith("OH ")];sc=[(int(f[2]),(int(f[3])-int(f[4]))/int(f[5]))for f in oh if f[1]=="scan"];cr={int(f[1]):(int(f[2]),int(f[3])-int(f[4])-int(f[5])-int(f[6]))for f in(ln.split()for ln in ol.splitlines()if ln.startswith("CELLRES "))};ZH=(cr[26][1]-cr[24][1])/(cr[26][0]-cr[24][0]);Yp=ZH<=M.CELL_ROUTE/4;ap=[(int(f[2])-int(f[3]))/int(f[4])for f in oh if f[1]=="apply"];need=max(v*b-M.C_SCAN_CELL*b for b,v in sc);Yj=Yp and set(range(1,65))<={b for b,v in sc}and need<=M.C_SCAN_LEAF and len(ap)==1 and ap[0]<=M.C_APPLY and M.C_ZERO==0 and ol.count("MATCH | survivors")>=4;nb=re.findall(r"([\d.]+)/(row|slice)","\n".join(Yc("nbrow.txt")));Yj&="mismatches=0"in Yc("nbrow.txt")[1]and float(nb[1][0])<=M.C_NB_SLICE;k2={f[1]:[int(x)for x in f[2:]]for f in map(str.split,Yc("nb2.txt"))if f[0]=="NB2"};n2=lambda v:(lambda m:m+2*(max(0,v[2]/100-v[0]*m*m)/(v[0]-1)/v[0])**.5)(v[1]/10/v[0]);Yj&="512 MATCH"in Yc("nb2.txt")[0]and all(v[0]==512 and v[2]*512>=v[1]**2 for v in k2.values())and M.C_NB_EVAL>=max(1154.6,n2(k2['kappa']))and M.C_NB_ROW>=float(nb[0][0])+n2(k2['eta'])and M.BOX+(M.S0,M.L1)==(8,8,8,64,16);Ya(f"[10s] scan, apply, bound rows {K(Yj)}");YW,YL=_a("bw_dot2.s","LBB0_2:");Ze,Zf=_a("krylov5.s",".LBB0_1:");Zi,Zj=_a("krylov5.s",".LBB1_1:");u5=Ze/5;UE=__import__('unit_entries');ku={k:(sum(v)/len(v),(sum((x-sum(v)/len(v))**2 for x in v)/(len(v)-1)/len(v))**.5)for k,v in M.K5_RUNS.items()};ZY=max(ku[k][0]+2*ku[k][1]for k in(16,17));YJ=Zf==['b.ne']and Zj==['b.ne']and u5==UE.U5==7.8 and Zi==UE.G1 and M.C_ENTRY_K>=ZY and max(ku[16][0],ku[17][0])<=ku[15][0]and max(M.K5_RUNS[16]+M.K5_RUNS[17])<=12.54;Ya(f"[10] krylov5: {Ze}/entry, {u5:.1f}, {Zi}; PREREG-U k5 "+", ".join(f"2^{k} {ku[k][0]:.4f}/{ku[k][1]:.4f}"for k in sorted(ku))+f"; max mean+2SE {ZY:.4f} <= {M.C_ENTRY_K} {K(YJ)}");Ya(f"[10] bw_dot2: {YW} (charged {M.C_PROJ:.0f}) {K(YW<=M.C_PROJ and YL==['b.ne'])}");au=audit(ZQ=dres['DET'][(220,225),'g']['pr_log2'],Z6=bt['pr_log2']);Zx=au['boot']['_d']and au['_b'];Ya(f"[11] boot {K(au['boot']['_d'])}; special-q {K(au['_b'])}; UDUP yr {au['Yk']['_i']:.2f}, {au['Yk']['_j']:.3f}; 8 lattices: pair ratio {au['Yk']['pair']:.4f} <= 2; chi2 {au['sampler']['chi2']:.1f} {K(au['sampler']['ok'])}");Y9=au['Yk']['_i']>1. and au['Yk']['_j']>.70 and au['Yk']['pair']*2**23<=Pm.BALL_PAIR;Y8=au['sampler']['ok'];import descent_witness as W;import retain16 as RT;import ecm16 as EC;st=_k();nb=st['nb_empty_window'];rr=[(b,s_,R_)+W.retention(b,s_,R_,seed+b,40)for b,s_,R_ in[(b,s_,dres['DET'][bb,s_]['R'])for b,s_ in((69,'g'),(150,'g'))for bb in dres['bins']if bb[0]==b]];YI=all(z<=zp for b_,s_,R_,lb_,m_,r_,z,zp in nb['rows']if m_>=7.4);tab=st['prime_structure']['table'];special=[p for p in RT.primes(2,512)if any(e!=1 or any(m!=1 for m in mf+mg)for d1,e,mf,mg in tab[str(p)])];Zq=sorted(set([p for i,p in enumerate(RT.primes(2,512))if i%7==seed%7]+special));_9=len(tab)==len(RT.primes(2,512))and all(RT.structure(p)==tab[str(p)]for p in Zq);sample={}
	for p in special:
		rs=W.residue_valuations(p,400,seed+p)
		for side in 'fg':bound=1-RT.local_pdf(p,tab[str(p)],side,False,True).get(0,.0);frac=1-rs[side].get(0,0)/rs['n'];sample[f"{p}{side}"]=round(frac,4),round(bound,4),frac<=bound+3*sqrt(bound*(1-bound)/rs['n'])+1e-12
	_x=all(v[2]for v in sample.values());pf,ordf=RT.omitted_pdf('f',True);pg,ordg=RT.omitted_pdf('g',True);ex=lambda pdf:sum(x*v for x,v in pdf.items())/RT.GRID;Ys=min(v['ret']for v in dres['DET'].values()if v['how']=='sieve');Yr=max(max(v['loss_f'],v['loss_g'])for v in dres['DET'].values());Zb=YI and all(x[3]for x in rr)and _9 and _x;Ya(f"[11b]: table {K(_9)}; residue {K(_x)}; loss <= {Yr:.1e}, retention >= {floor(Ys*1e4)/1e4:.4f}; toys {K(YI)}");mm,Y5=W.encountered_check(seed);Ya(f"[7e] enc/sampled: mean cf {mm[0][0]:.1f}/{mm[0][1]:.1f}, cg {mm[1][0]:.1f}/{mm[1][1]:.1f}; max {mm[0][2]:.1f}/{mm[1][2]:.1f} <= 97/14.5 {K(Y5)}");YM=EC.toy(60,seed,sizes=((36,1024,50),(40,2048,50)));Xa=[EC.schedule(69.,int(ceil(dres['DET'][(69,75),'f']['nf']))),EC.schedule(bt['nA'],int(ceil(bt['m'])))];Y7=all(t['empirical']>=t['model_corrected']-2*sqrt(t['model_corrected']*(1-t['model_corrected'])/t['trials'])for t in YM)
	Ya(f"[12] ECM: "+"; ".join(f"{x['fbits']:.2f}/{x['nbits']} bits: {x['curves']} curves, 2^{x['log2_instructions']:.2f}"for x in Xa)+f" {K(Y7)}")
	YG=True;wo=""
	for kw,pre,(n_,a_,b_,s_,k_)in pd['toy2']:wt=M.toy(**kw,trials=pre[0]);q=n_/pre[0];rt=wt['ret'];YG&=[wt[x]for x in('A','B','smooth','kept')]==pre[1:]and a_<=q*wt['sA']+2*(q*wt['sA'])**.5 and b_<=q*wt['sB']+2*(q*wt['sB'])**.5 and s_>=1000 and k_/s_>=rt-2*sqrt(rt*(1-rt)/s_);wo+=f" 2^{kw['n']}: A {a_} <= {q*wt['sA']:.0f}, B {b_} <= {q*wt['sB']:.0f}, kept {k_}/{s_} >= {rt:.4f};"
	Ya(f"[12w] toys ({n_} integers; replay {pre[0]}):{wo} {K(YG)}")
	fq=pd['fg'];fn=sum(map(sum,fq));fs=lambda n:cost.survivor_upper(n,op['slack_f'],.25,op['slack_b'],op['lb'],M.KMAX);ft=[math.log(fs(732)/fs(740))/8,math.log(M.rho((M.FILT_MEANS[1]-8)/op['lb'])/M.rho(M.FILT_MEANS[1]/op['lb']))/8]
	def fg(k):
		X,Y=[[math.exp(k*a*(8*i+4))for i in range(len(fq)+len(fq[0]))]for a in ft];A,B,G=[sum(v*w(i,j)for i,r in enumerate(fq)for j,v in enumerate(r))/fn for w in(lambda i,j:X[i]*Y[j],lambda i,j:X[i],lambda i,j:Y[j])];return A/B/G*(1+2*sqrt((sum(v*(X[i]*Y[j]/A-X[i]/B-Y[j]/G)**2 for i,r in enumerate(fq)for j,v in enumerate(r))/fn-1)/fn))
	fd=fg(1),fg(2);Yfg=fn>=8000 and min(ft)>0 and max(fd)<=M.Q_F_SAFETY;Ya(f"[12q] f/g small parts of {fn} pairs: joint/product at the events' tilts ({ft[0]:.4f}, {ft[1]:.4f} per bit) <= {fd[0]:.3f}, at twice <= {fd[1]:.3f} (+ 2 SE) <= {M.Q_F_SAFETY:g} {K(Yfg)}")
	fl=_z()
	Ya(f"[13]: N(f(X)) deg {fl['deg_f']} and N(g(X)) deg {fl['deg_g']} squarefree {K(fl['unramified'])}; "f"K_f {fl['sig']['f']}, K_f+ {fl['sig']['fp']}, K_g {fl['sig']['g']}: odd rank {fl['sig']['odd']}, chi-part {fl['sig']['sm']} = SM {M.N_SM} {K(fl['sig']['sm']==M.N_SM)}")
	Zw=fl['unramified']and fl['sig']['sm']==M.N_SM and fl['sig']['odd']==2*M.N_SM;import runs16 as local2;ZI={p:_m(p)for p in(2,3,19)};_M=min(op['alpha_f'],DS.P.ALPHA_F)>=max(ALPHA_F,_r)and min(op['alpha_g'],DS.P.ALPHA_G)>=max(ALPHA_G,_q)and local2.facts()[0]
	Ya(f"[14] alpha: ({ALPHA_F:.6f}, {ALPHA_G:.6f}; 50k ({_r:.6f}, {_q:.6f})) <= ({op['alpha_f']}, {op['alpha_g']}) {K(_M)}; "f"in-run p=2 (exact, local2)/3/19: "+", ".join(f"{p}: ({v[0]:.6f}, {v[1]:.6f})"for p,v in sorted(ZI.items()))+f"")
	ev=local2.expectations();kp=M.kappa2(840.29/69.5),M.kappa2(530.75/69.5);YB=all(0<=M.V2E[i]-ev[i]<1e-6 and kp[i]<=M.KAP2[i]for i in(0,1));YQ=toy_v2(_k()['toy_v2']);Ya(f"[14v] v2: E {M.V2E} >= {[float(x)for x in ev]}; kappa2 {kp} <= {M.KAP2} {K(YB)}");tot=M.total_log2(e['log2_relcol'],la['log2_la'],cap,la['log2_gen']);gK,glo=M.GATE;ZD=log2(1+gK/2**64*(2**M.total_log2(e['log2_relcol'],la['log2_la'],two_log,la['log2_gen'])+2**66))
	Ya(f"total_log2W_full_ledger = {tot:.12f} (every run <= {CLAIM:.10f} {K(tot<=CLAIM)})")
	vals=dict(ledger=tot,relcol=e['log2_relcol'],la=la['log2_la'],gen=la['log2_gen'],log2_rels=e['log2_rels'],rels_2se=Yd,margin=e['margin'],work_se=work_se,two_log=two_log,p_cap=pcap,log2_EN=log2_EN,over=log2(Zn),rep=log2(Z8),boot_phi=bt['phi'],success=success,success_head=success_head,expected=ZD)
	_Z={k:vals[k]<=v if op=='<='else vals[k]>=v for k,(op,v)in _W.items()};_Z['counts']=ZO;_E=all(_Z.values());Ya("declared: "+"; ".join(f"{k} {vals[k]!r} {op} {v} {K(_Z[k])}"for k,(op,v)in _W.items()))
	checks=dict(owner_schedule=O.check(M,la),resultant=Xh,tower=Zy,rho_semismooth=Xg,pool=Xf and _2 and Yo and Y3,box_bounds=Y6,nbar=Xe,la_loop=YJ and Yy==97 and YK==['b.ne']and YW<=M.C_PROJ and YL==['b.ne']and Yi and tail<=M.C_ROW and Y1<=M.C_RED9 and _T and YH and _0 and _L and _S and Yj,rec7=ZM and YE and ZN,band_split=_K,srt=Qa and Xx and Yfg,tdiv=_1,filt_y=YD,hx=Yq and Qc,ball_counts=_y and _D,lazy=_3,rels_input=_l,work_allowance=_N,declared_bounds=_E,alpha_ok=_M and YB and YQ,missing_logs=Qb,record_widths=_J,descent_plan_feasible=mu_min>=Pm.MIN_RELS and psi_max<=Pm.PSI_MAX and max(r['phi']for r in rows)<=.05 and Yx<=-4,audit=Zx,retention=Zb,ecm=Y7 and YG,fields=Zw,feasible=e['feasible']and Xw,within_claim=tot<=CLAIM,success_ok=success>=.39 and success_head>=.39,udup_relw=Y9,sampler_uniform=Y8,encountered=Y5)
	Ya(f"done in {time.time()-T0:.1f}s; checks: {checks}");print(json.dumps({"checks_ok":all(checks.values()),"checks":checks,"total_log2W_full_ledger":ceil(tot*10**12)/10**12,"claimed":CLAIM,"success_exact":success}),flush=True)
if __name__=="__main__":
	if "--alpha50k"in sys.argv:print(json.dumps(dict(zip(("alpha_f","alpha_g"),_p(_A))),sort_keys=True));sys.exit(0)
	main()

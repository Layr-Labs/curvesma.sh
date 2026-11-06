"""stnfs descent model and shipped plan (DESCENT.md)."""
import math,functools
from decimal import Decimal,localcontext
LN2=math.log(2);UNIT=20000.0
KMAX=64
RSTEP=1.0/1024
def _rho_coefficients(kmax=KMAX,nt=600):
	with localcontext()as ctx:
		ctx.prec=200
		zero,one=Decimal(0),Decimal(1)
		coefs=[None,[one]+[zero]*(nt-1)]
		for k in range(2,kmax+1):
			previous,c=coefs[k-1],[zero]*nt
			for m in range(nt-1):
				c[m+1]=(previous[m]+m*c[m])/(k*(m+1))
			c[0]=previous[0]-sum(c[1:])
			coefs.append(c)
		out=[None]
		for c in coefs[1:]:
			big=max(abs(x)for x in c)
			while len(c)>1 and abs(c[-1])<big*Decimal('1e-25'):
				c=c[:-1]
			out.append([float(x)for x in c])
		return out
RHO_COEFS=_rho_coefficients()
def rho_series(u):
	if u<=1:
		return 1.0
	k=math.ceil(u)
	if k>KMAX:
		return 0.0
	x,value=k-u,0.0
	for c in reversed(RHO_COEFS[k]):
		value=value*x+c
	return max(value,0.0)
RHO_TAB=[rho_series(i*RSTEP)for i in range(int(KMAX/RSTEP)+2)]
def rho_lo(u):
	if u<=1:
		return 1.0
	i=math.ceil(u/RSTEP)
	return RHO_TAB[i]if i<len(RHO_TAB)else 0.0
def rho_hi(u):
	if u<=1:
		return 1.0
	i=math.floor(u/RSTEP)
	if i+1>=len(RHO_TAB):
		return 1.0/math.factorial(math.floor(u))
	return RHO_TAB[i]
def kappa(u,t):
	return rho_series(u-t)/rho_series(u)if rho_series(u)>0 else 0.0
@functools.lru_cache(None)
def miss_side(u4,log2_rels,u_main):
	u=u4/4.0;tot=0.0
	for b10 in range(400,690):
		b=b10/10.0
		lam=2.0**(log2_rels-b)*2*kappa(u_main,b/69.0)
		tot+=kappa(u,b/69.0)*math.log((b+0.1)/b)*math.exp(-lam)*(1+lam+lam*lam/2)
	return tot
def missing(Pm,nf,ng):
	f=lambda m,um:miss_side(math.ceil(4*m/69.0),Pm.LOG2_RELS,um)
	return(f(nf,Pm.MAIN_UF)if nf else 0.0)+(f(ng,Pm.MAIN_UG)if ng else 0.0)
KA_MAX=16
def _conv(a,b):
	out=[0.0]*(len(a)+len(b)-1)
	for i,x in enumerate(a):
		if x:
			for j,y in enumerate(b):
				if y:
					out[i+j]+=x*y
	return out
@functools.lru_cache(None)
def window_rows(n,base,iA,mesh,kmax=KA_MAX):
	h=(n-base)/mesh
	wA=[math.log1p(h/(base+i*h))if i<iA else 0.0 for i in range(mesh)]
	wB=[math.log1p(h/(base+i*h))if i>=iA else 0.0 for i in range(mesh)]
	PA=[[1.0]]
	for j in range(1,kmax+1):
		PA.append([x/j for x in _conv(PA[-1],wA)])
	PB=[[1.0],wB]
	rows={}
	for jA in range(kmax+1):
		for jB in range(2):
			rows[(jA,jB)]=PA[jA]if jB==0 else _conv(PA[jA],PB[1])
	return rows,h
@functools.lru_cache(None)
def ss_cum2(m,base,n,iA,mesh,alpha,trunc,kmax=KA_MAX):
	rows,h=window_rows(n,base,iA,mesh,kmax)
	out=[[0.0,0.0]for _ in range(kmax+1)]
	for(jA,jB),row in rows.items():
		j=jA+jB
		imax=len(row)-1
		ifeas=min(imax,int(math.floor((m-j*h-j*base)/h)))if h>0 else imax
		if ifeas<0:
			continue
		imin=0 if not trunc else max(0,ifeas-int(trunc*base/h)-1)
		tab,ntab,ceil=RHO_TAB,len(RHO_TAB),math.ceil
		k0=(m-j*base+alpha)/(base*RSTEP);hk=h/(base*RSTEP)
		s=0.0
		for i in range(imin,ifeas+1):
			weight=row[i]
			if weight:
				k=ceil(k0-i*hk)
				s+=weight*(1.0 if k<=0 else(tab[k]if k<ntab else 0.0))
		out[jA][jB]=s
	cum=[[0.0,0.0]for _ in range(kmax+1)]
	for jB in range(2):
		acc=0.0
		for jA in range(kmax+1):
			acc+=out[jA][jB]
			cum[jA][jB]=acc
	return tuple(tuple(x)for x in cum)
def ss_lower2(m,base,n,iA,kA,kB,mesh,alpha,trunc=3.0):
	m=math.ceil(m)
	if n<=base:
		return rho_lo((m+alpha)/base)
	cum=ss_cum2(m,base,n,iA,mesh,alpha,trunc,kA)
	return cum[kA][0]+(cum[kA][1]if kB else 0.0)
def surv_upper(m,base,T):
	if T<=0:
		return 0.0
	if T>=m:
		return 1.0
	return min(1.0,max(1.0,T/base)*rho_hi((math.floor(m)-T)/base))
MONT5=367.76
def ecm_instr(fbits,nbits,success):
	import ecm16
	return ecm16.ecm_instr(fbits,nbits,success)
def log2li(l):
	x=l*LN2
	return l-math.log2(x)+math.log2(1+1/x+2/x/x)
CELLS={4:13315520,6:1573359464,8:47982017060,10:687193400292,12:6073556677752,14:38414493101628,
16:190058478539564,20:2754504316057176,24:24505531366261200,28:155610653675916932,
32:771909602591773000,40:11221331971811657968,48:99995402971660463340,
56:635616373364584622852,64:3155035383949694198996,80:45900097950907999286464,
96:409193506988359583551012,112:2601677688965738307448548,128:12916156661164168796080328,160:187942845067274337246511324,192:1675659512501707896389090868,
224:10654623198171176305341831664,256:52897551248440914796084637156,288:217401149762189752879411886484,
320:769748406668707279529221545272,384:6863095027560050329435265332068}
class P:
	LB=69.0;LX=69.0;PMIN=8.0;SLACK=96.0
	ALPHA_F=5.15;ALPHA_G=0.85;LOG2U6=6*63.714
	QBINS=[(1/6,-26.4+23,14.5+6),(1/6,-11.2+23,14.5+6),(1/6,-0.1+23,14.5+6),
	(1/6,10.5+23,14.5+6),(1/6,23.2+23,14.5+6),(1/6,86.6+23,14.5+6)]
	BOOT_BINS=[(0.5,1.4+4),(0.4,6.8+4),(0.1,12.7+4)]
	QBINS_LO=[(-140.05-23,-54.56-6),(-26.4-23,-54.56-6),(-11.2-23,-54.56-6),
	(-0.1-23,-54.56-6),(10.5-23,-54.56-6),(23.2-23,-54.56-6)]
	QBINS_SUB=[[-26.47,-0.76],[-15.7,1.36],[-8.64,2.41],[-3.29,4.1],[1.16,5.45],[5.07,4.44],[8.45,5.05],[11.67,5.09],[14.69,6.45],[17.52,5.98],[20.28,6.25],[22.96,6.93],[25.61,7.17],[28.26,8.98],[30.95,9.36],[33.69,9.07],[36.56,9.25],[39.6,10.43],[42.84,10.22],[46.5,11.73],[50.67,11.47],[55.75,12.05],[63.02,13.19],[109.6,20.5]]
	BOOT_SUB=[-13.24,-9.46,-7.11,-5.11,-3.15,-1.19,0.59,2.53,5.2,16.7]
	PRCAP_NODE=4.0;PRCAP_BOOT=1.5
	LOG2_RELS=63.575;MAIN_UF=739.1/69;MAIN_UG=498.2/69
	THIN_MAX=85
	QAUX=20
	UDUP=0.65
	MIN_RELS=18.0;PSI_MAX=0.8;PHI_FLOOR=0.02
	BOOT_RELS=18
	NB_R1=0.5;NB_R2=1.0
	BOUND_CELL=2400.0
	HITS_PER_CELL=2*math.log(69.0/8.0)
	SIEVE_HIT=1100.0;SCAN=10.0;ZERO=2.0
	PAIR_SETUP=2.0**25.1;LLL_NODE=2.0**36;NORMS=2.0**26;TDIV=2.0**24.5;IDEAL=2.0**34;PRIME=2.0**24
	ECM_SUCCESS=0.999
	RADII=(4,6,8,10,12,14,16,20,24,28,32,40,48,56,64,80,96,112,128,160,192,224,256,288,320,384)
	MESH=384;TRUNC=0.0
	LOG2P6=6*380.99700
	LMAX=225;BIN=5
def pairs_count(Pm):
	return 2.0*2**log2li(Pm.LX)
def node_cost(Pm,l_lo,l_hi,side,n,iA,kA,kB,R,thin=None):
	pairs=pairs_count(Pm)
	cfmax=max(b[1]for b in Pm.QBINS);cgmax=max(b[2]for b in Pm.QBINS)
	S=math.log2(CELLS[R]);c=(l_hi+S)/12
	h=(n-Pm.LX)/Pm.MESH if n>Pm.LX else 0.0
	nA=Pm.LX+iA*h if n>Pm.LX else Pm.LX
	T=kA*nA+kB*n+Pm.SLACK
	c_lo=(l_lo+S)/12
	pr=0.0;sv=0.0
	for(w,cf,cg),(cf_lo,cg_lo)in zip(Pm.QBINS,Pm.QBINS_LO):
		nf=72*c+cf-(l_lo if side=='f' else 0)
		ng=12*c+Pm.LOG2U6+cg-(l_lo if side=='g' else 0)
		pf=ss_lower2(nf,Pm.LX,n,iA,kA,kB,Pm.MESH,Pm.ALPHA_F,Pm.TRUNC)
		pg=ss_lower2(ng,Pm.LX,n,iA,kA,kB,Pm.MESH,Pm.ALPHA_G,Pm.TRUNC)
		pr+=w*pf*pg
		nf_lo=72*c_lo+cf_lo-(l_hi if side=='f' else 0)
		ng_lo=12*c_lo+Pm.LOG2U6+cg_lo-(l_hi if side=='g' else 0)
		sv+=w*surv_upper(nf_lo,Pm.LX,T)*surv_upper(ng_lo,Pm.LX,T)
	pri=0.0
	for cf,cg in Pm.QBINS_SUB:
		nf=72*c+cf-(l_lo if side=='f' else 0)
		ng=12*c+Pm.LOG2U6+cg-(l_lo if side=='g' else 0)
		pri+=ss_lower2(nf,Pm.LX,n,iA,kA,kB,Pm.MESH,Pm.ALPHA_F,Pm.TRUNC)*ss_lower2(ng,Pm.LX,n,iA,kA,kB,Pm.MESH,Pm.ALPHA_G,Pm.TRUNC)
	pr=min(pri/24,Pm.PRCAP_NODE*pr)
	if pr<=0:
		return None
	nf_max=72*c+cfmax-(l_lo if side=='f' else 0)
	ng_max=12*c+Pm.LOG2U6+cgmax-(l_lo if side=='g' else 0)
	det=1.0-(1-Pm.ECM_SUCCESS)*((nf_max+ng_max)/20.0+2*(kA+kB))
	miss=missing(Pm,nf_max,ng_max)
	if thin is None or l_lo<=Pm.THIN_MAX:
		thin=thin or l_lo<=Pm.THIN_MAX
	mu=CELLS[R]*pr*Pm.UDUP*det*(math.exp(-miss)if thin else 1.0)
	if mu<Pm.MIN_RELS:
		return None
	surv_cost=(Pm.NORMS+Pm.TDIV+ecm_instr(Pm.LX,nf_max,Pm.ECM_SUCCESS)+ecm_instr(Pm.LX,ng_max,Pm.ECM_SUCCESS)
	+2*(ecm_instr(nA,T,Pm.ECM_SUCCESS)+Pm.PRIME)+Pm.IDEAL)
	survivors=CELLS[R]*sv+CELLS[R]*pr
	sieve=CELLS[R]*(Pm.HITS_PER_CELL*Pm.SIEVE_HIT+Pm.SCAN+Pm.ZERO+Pm.BOUND_CELL)+pairs*Pm.PAIR_SETUP
	trials=CELLS[R]*surv_cost
	if trials<sieve+survivors*surv_cost:
		own=Pm.LLL_NODE+trials;how='ecm';survivors=CELLS[R]
		ret,loss_f,loss_g=1.0,0.0,0.0
	else:
		own=Pm.LLL_NODE+sieve+survivors*surv_cost;how='sieve'
		import retain16
		ret,loss_f,loss_g=retain16.retention(nf_max,ng_max,n,iA,kA,kB,Pm.SLACK)
		mu*=ret
		if mu<Pm.MIN_RELS:
			return None
	return dict(own=own,n=n,iA=iA,nA=nA,kA=kA,kB=kB,R=R,S=S,mu=mu,how=how,miss=miss,ret=ret,loss_f=loss_f,loss_g=loss_g,
	miss_children=0.0 if thin else miss,thin=bool(thin),pr_log2=math.log2(pr),
	surv_log2=math.log2(survivors),own_log2=math.log2(own),nf=nf_max,ng=ng_max,det=det,T=T,
	sieve_log2=math.log2(sieve),survcost_log2=math.log2(survivors*surv_cost))
def bins_of(Pm):
	bins=[(int(Pm.LB),int(Pm.LB)+Pm.BIN+1)]
	l=int(Pm.LB)+Pm.BIN+1
	while l<Pm.LMAX:
		bins.append((l,l+Pm.BIN));l+=Pm.BIN
	return bins
def phi_nb(mean,Pm,window):
	r=Pm.NB_R1 if window<=1 else Pm.NB_R2
	return max((1.0+mean/r)**(-r),Pm.PHI_FLOOR)
def lazy_fraction(mu,r,psi,tol=1e-13):
	retry=(1.0+1e-12)/(1.0-psi)
	lq=math.log(mu/(r+mu));lp0=r*math.log(r/(r+mu))
	tot=mass=0.0;m=0
	while m<=mu or 1.0-mass>tol:
		pm=math.exp(math.lgamma(m+r)-math.lgamma(r)-math.lgamma(m+1)+lp0+m*lq)
		tot+=pm*min(1.0,retry/(m+1)+psi**m);mass+=pm;m+=1
		if m>mu and pm<1e-20:break
	return min(1.0,tot+max(0.0,1.0-mass)+pm*(r+mu)/r)
def lazy_mc(mu,r,psi,s=10**6,trials=4000,seed=7):
	import random
	rng=random.Random(seed);tot=0.0
	for _ in range(trials):
		m=min(s,_poisson(rng,rng.gammavariate(r,mu/r)))
		pos=sorted(rng.sample(range(1,s+1),m));k=0
		while k<m and rng.random()<psi:k+=1
		tot+=(pos[k]if k<m else s)/s
	return tot/trials
def _poisson(rng,lam):
	if lam<50:
		L,k,p=math.exp(-lam),0,1.0
		while True:
			p*=rng.random()
			if p<=L:return k
			k+=1
	return max(0,int(round(rng.gauss(lam,math.sqrt(lam)))))
def evaluate(Pm,PLAN,BOOT):
	shared=('LX','MESH','TRUNC','ALPHA_F','ALPHA_G')
	if any(getattr(Pm,name)!=getattr(P,name)for name in shared):
		raise ValueError('descent and retention parameters must agree')
	import retain16
	retain16.loss.cache_clear()
	bins=bins_of(Pm)
	W={};PHI={};DET={}
	def wmax(side,n):
		return max([W[(b,side)]for b in bins if b[0]<n-1e-9]+[0.0])
	def phimax(side,n):
		return max([PHI[(b,side)]for b in bins if b[0]<n-1e-9]+[0.0])
	def miss_terms():
		low=[k for k in W if k[0][0]<=Pm.THIN_MAX]
		return(max([PHI[k]for k in low]+[0.0]),max([W[k]for k in low]+[0.0]))
	for b in bins:
		for side in('f','g'):
			n,iA,kA,kB,R,thin=PLAN[(b[0],side)]
			d=node_cost(Pm,b[0],b[1],side,float(n),iA,kA,kB,R,bool(thin))
			if d is None:
				raise ValueError('infeasible plan entry %s %s'%(b,side))
			nA=d['nA']
			pm,wm=miss_terms()
			psi=kA*(phimax('f',nA)+phimax('g',nA))+kB*(phimax('f',n)+phimax('g',n))+d['miss_children']*pm
			if psi>Pm.PSI_MAX:
				raise ValueError('psi too large at %s %s'%(b,side))
			retry=1.0/(1.0-psi)
			sub=kA*(wmax('f',nA)+wmax('g',nA))+kB*(wmax('f',n)+wmax('g',n))+d['miss_children']*wm
			window=(kA if nA>Pm.LX+1e-9 else 0)+(kB if n>nA+1e-9 else 0)
			if d['how']=='sieve':
				frac=lazy_fraction(d['mu'],Pm.NB_R1 if window<=1 else Pm.NB_R2,psi)
				d['own']=Pm.LLL_NODE+2**d['sieve_log2']+2**d['survcost_log2']*frac
				d['lazy']=frac;d['own_log2']=math.log2(d['own']);d['nb_r']=Pm.NB_R1 if window<=1 else Pm.NB_R2
			W[(b,side)]=d['own']+retry*sub
			PHI[(b,side)]=phi_nb(d['mu']*(1-psi),Pm,window)
			DET[(b,side)]=dict(d,psi=psi,retry=retry,sub_log2=math.log2(sub)if sub>0 else None,
			W_log2=math.log2(W[(b,side)]),phi=PHI[(b,side)])
	n1,iA,kA,Rb=BOOT
	if n1>Pm.LMAX:
		raise ValueError('boot window above the last bin')
	Sb=math.log2(CELLS[Rb]);cb=(Pm.LOG2P6+Sb)/12
	h=(n1-Pm.LX)/Pm.MESH;nA=Pm.LX+iA*h
	m_max=12*cb+Pm.LOG2U6+max(x[1]for x in Pm.BOOT_BINS)
	pr=sum(w*ss_lower2(12*cb+Pm.LOG2U6+cg,Pm.LX,float(n1),iA,kA,1,Pm.MESH,Pm.ALPHA_G,Pm.TRUNC)
	for w,cg in Pm.BOOT_BINS)
	pri=sum(ss_lower2(12*cb+Pm.LOG2U6+cg,Pm.LX,float(n1),iA,kA,1,Pm.MESH,Pm.ALPHA_G,Pm.TRUNC)
	for cg in Pm.BOOT_SUB)
	pr=min(pri/10,Pm.PRCAP_BOOT*pr)
	det=1.0-(1-Pm.ECM_SUCCESS)*(m_max/20.0+kA+1)
	import retain16
	per_rel=1.0/(pr*det*Pm.UDUP*retain16.counted_fraction())
	trial=Pm.NORMS+ecm_instr(nA,m_max,Pm.ECM_SUCCESS)+Pm.PRIME+Pm.IDEAL
	pm,wm=miss_terms();mb=missing(Pm,0,m_max)
	psi=kA*phimax('g',nA)+phimax('g',float(n1))+mb*pm
	if psi>Pm.PSI_MAX:
		raise ValueError('boot psi too large')
	retry=1/(1-psi)
	trials=per_rel*(1-psi**Pm.BOOT_RELS)/(1-psi)
	balls=trials # fresh exponent/lattice per cell trial
	sub=kA*wmax('g',nA)+wmax('g',float(n1))+mb*wm
	one_log=(2.0**36.1)*balls+trials*trial+retry*sub
	phi_boot=psi**Pm.BOOT_RELS
	boot=dict(n1=n1,nA=nA,kA=kA,R=Rb,S=Sb,m=m_max,pr_log2=math.log2(pr),trials_log2=math.log2(trials),
	trial_log2=math.log2(trial),balls=balls,miss=mb,own_log2=math.log2(trials*trial),sub_log2=math.log2(sub),
	psi=psi,retry=retry,phi=phi_boot,one_log_log2=math.log2(one_log))
	return dict(W=W,PHI=PHI,DET=DET,bins=bins,boot=boot,two_log_units=math.log2(2*one_log/UNIT),
	one_log_instr=one_log)
_PLAN={
69:((69,48,0,0,256,1),(69,48,0,0,288,1)),
75:((75,42,1,0,224,1),(75,42,1,0,256,1)),
80:((75,42,1,0,288,0),(75,47,2,0,192,0)),
85:((84,19,2,0,192,1),(75,44,2,0,192,1)),
90:((75,45,2,0,224,1),(78,32,2,0,256,1)),
95:((80,26,2,0,288,1),(80,26,2,1,224,1)),
100:((80,13,2,1,288,1),(80,21,2,1,384,1)),
105:((85,18,2,1,256,1),(85,18,2,1,384,1)),
110:((85,21,2,1,320,1),(90,22,2,1,320,1)),
115:((90,24,2,1,288,1),(95,24,2,1,320,1)),
120:((95,20,3,1,224,1),(96,45,3,0,320,1)),
125:((95,24,2,1,320,0),(99,20,3,1,320,0)),
130:((99,24,3,1,384,1),(100,24,4,1,384,0)),
135:((100,28,3,1,384,1),(110,28,3,1,320,1)),
140:((100,30,3,1,384,0),(110,30,3,1,320,0)),
145:((110,30,3,1,288,1),(109,36,4,1,288,1)),
150:((109,33,3,1,384,1),(105,41,4,1,384,0)),
155:((109,25,4,1,384,0),(125,26,4,1,384,1)),
160:((115,26,4,1,320,0),(124,22,5,1,384,0)),
165:((125,24,4,1,384,0),(128,36,4,1,384,0)),
170:((125,30,4,1,384,1),(135,36,4,1,384,0)),
175:((136,30,4,1,288,1),(145,29,4,1,384,0)),
180:((130,36,4,1,320,0),(124,36,6,1,320,0)),
185:((140,30,4,1,384,0),(133,30,6,1,384,0)),
190:((145,30,4,1,384,0),(150,24,6,1,384,0)),
195:((150,30,4,1,384,0),(159,37,4,1,384,0)),
200:((150,36,4,1,320,0),(144,36,6,1,288,0)),
205:((165,30,4,1,384,1),(155,30,6,1,384,0)),
210:((160,36,4,1,320,0),(165,28,6,1,384,0)),
215:((165,36,4,1,320,0),(170,36,6,1,224,0)),
220:((170,36,4,1,384,1),(170,36,6,1,288,0)),
}
PLAN16={(b,s):(v[i][0],8*v[i][1],*v[i][2:])
for b,v in _PLAN.items()for i,s in enumerate('fg')}
BOOT16=(215,328,16,4)
def report(res):
	rows=[]
	for b in res['bins']:
		for side in('f','g'):
			d=res['DET'][(b,side)]
			rows.append(dict(bin=b[0],side=side,n=int(d['n']),nA=round(d['nA'],1),kA=d['kA'],kB=d['kB'],R=d['R'],ret=round(d['ret'],4),
			loss_f=d['loss_f'],loss_g=d['loss_g'],how=d['how'],mu=round(d['mu'],2),own_log2=round(d['own_log2'],3),
			W_log2=round(d['W_log2'],3),psi=round(d['psi'],4),phi=round(d['phi'],6),nf=round(d['nf'],1),ng=round(d['ng'],1),T=round(d['T'],1)))
	return rows

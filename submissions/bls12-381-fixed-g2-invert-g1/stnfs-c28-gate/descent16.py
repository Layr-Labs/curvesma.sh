""
import math,functools
from math import exp,factorial,floor,lgamma,log,log1p,log2,sqrt
from decimal import Decimal,localcontext
LN2=log(2);UNIT=20000.0
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
	i=floor(u/RSTEP)
	if i+1>=len(RHO_TAB):
		return 1.0/factorial(floor(u))
	return RHO_TAB[i]
def kappa(u,t):
	return rho_series(u-t)/rho_series(u)if rho_series(u)>0 else 0.0
@functools.lru_cache(None)
def miss_side(u4,log2_rels,u_main,orbit=4):
	u=u4/4.0;tot=0.0
	for b10 in range(400,690):
		b=b10/10.0
		lam=2.0**(log2_rels-b)*orbit*kappa(u_main,b/69.0)
		tot+=kappa(u,b/69.0)*log((b+0.1)/b)*exp(-lam)*(1+lam+lam*lam/2)
	return tot
def missing(Pm,nf,ng):
	f=lambda m,um:miss_side(math.ceil(4*m/69.0),Pm.LOG2_RELS,um,Pm.ORBIT)
	return(f(nf,Pm.MAIN_UF)if nf else 0.0)+(f(ng,Pm.MAIN_UG)if ng else 0.0)
KA_MAX=16
def _conv(a,b):
	out=[0.0]*(len(a)+len(b)-1);nz=[j for j,y in enumerate(b)if y]
	if nz:
		j0,j1=nz[0],nz[-1]+1;bb=b[j0:j1]
		for i,x in enumerate(a):
			if x:out[i+j0:i+j1]=[o+x*y for o,y in zip(out[i+j0:i+j1],bb)]
	return out
@functools.lru_cache(None)
def window_rows(n,base,iA,mesh,kmax=KA_MAX):
	h=(n-base)/mesh
	wA=[log1p(h/(base+i*h))if i<iA else 0.0 for i in range(mesh)]
	wB=[log1p(h/(base+i*h))if i>=iA else 0.0 for i in range(mesh)]
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
		ifeas=min(imax,int(floor((m-j*h-j*base)/h)))if h>0 else imax
		if ifeas<0:
			continue
		imin=0 if not trunc else max(0,ifeas-int(trunc*base/h)-1)
		tab,ntab,ceil=RHO_TAB,len(RHO_TAB),math.ceil
		k0=(m-j*base+alpha)/(base*RSTEP);hk=h/(base*RSTEP)
		s=0.0
		for t in[row[i]*(1.0 if(k:=ceil(k0-i*hk))<=0 else(tab[k]if k<ntab else 0.0))for i in range(imin,ifeas+1)if row[i]]:s+=t
		out[jA][jB]=s
	cum=[[0.0,0.0]for _ in range(kmax+1)]
	for jB in range(2):
		acc=0.0
		for jA in range(kmax+1):
			acc+=out[jA][jB]
			cum[jA][jB]=acc
	return tuple(tuple(x)for x in cum)
def _rsum(row,m,base,h,j,alpha,up):
	return sum(w*(rho_hi((m-j*base-(i+j)*h+alpha)/base)if up else rho_lo((m-j*base-i*h+alpha)/base))
	for i,w in enumerate(row)if w and m-j*base-(i+(0 if up else j))*h>=0)
@functools.lru_cache(None)
def child_load(m,base,n,iA,kA,kB,mesh,alpha,g):
	"H1 mean of sum g(window primes) at norm m"
	rows,h=window_rows(n,base,iA,mesh,KA_MAX);num=den=0.0
	gA=[log1p(h/(base+i*h))*g[i]*(i<iA)for i in range(mesh)];gB=[log1p(h/(base+i*h))*g[i]*(i>=iA)for i in range(mesh)]
	for jA in range(kA+1):
		for jB in range(kB+1):
			j=jA+jB;den+=_rsum(rows[(jA,jB)],m+1,base,h,j,alpha,0)
			if jA:num+=_rsum(_conv(rows[(jA-1,jB)],gA),m,base,h,j,alpha,1)
			if jB:num+=_rsum(_conv(rows[(jA,0)],gB),m,base,h,j,alpha,1)
	return num/den
def cell_values(bins,V,side,base,n,mesh):
	h=(n-base)/mesh
	return tuple(max([V[(b,side)]for b in bins if b[0]<base+(i+1)*h and b[1]>base+i*h]+[0.0])for i in range(mesh))
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
	return min(1.0,max(1.0,T/base)*rho_hi((floor(m)-T)/base))
MONT5=367.76
def ecm_instr(fbits,nbits,success):
	import ecm16
	return ecm16.ecm_instr(fbits,nbits,success)
def log2li(l):
	x=l*LN2
	return l-log2(x)+log2(1+1/x+2/x/x)
import common12
_D=common12.stats('plan16.b85')
CELLS=dict(_D['cells'])
class P:
	LB=69.0;LX=69.0;PMIN=8.0;SLACK=96.0
	ALPHA_F=0.29;ALPHA_G=0.75;LOG2U6=6*63.714
	QBINS=[(1/6,-19.5+14,9.5+5),(1/6,-5.9+14,9.5+5),(1/6,4.3+14,9.5+5),(1/6,13.8+14,9.5+5),(1/6,25.1+14,9.5+5),(1/6,83.0+14,9.5+5)]
	BOOT_BINS=[(0.5,-185.0+2),(0.4,-182.1+2),(0.1,-179.6+2)]
	QBINS_LO=[(-148.3-14,-54.0-5),(-19.5-14,-54.0-5),(-5.9-14,-54.0-5),(4.3-14,-54.0-5),(13.8-14,-54.0-5),(25.1-14,-54.0-5)]
	QBINS_SUB=_D['qbins_sub']
	BOOT_SUB=_D['boot_sub']
	PRCAP_NODE=4.0;PRCAP_BOOT=1.5
	LOG2_RELS=62.583;MAIN_UF=740.0/69;MAIN_UG=498.4/69;ORBIT=4
	THIN_MAX=85
	UDUP=0.65
	MIN_RELS=18.0;LOAD=1.001;BOOT_DISC=0.9;CD=1.2;YF=0.43;OD=3.0
	def zeta(self,mu,key):return exp(-mu/self.OD)
	TAB=_D.get('loads');CHECK=()
	BOUND_CELL=2400.0
	HITS_PER_CELL=2*log(69.0/8.0)
	SIEVE_HIT=1100.0;SCAN=10.0;ZERO=2.0
	PAIR_SETUP=2.0**25.1;LLL_NODE=2.0**36;NORMS=2.0**26;TDIV=2.0**24.5;IDEAL=2.0**34;PRIME=2.0**24
	ECM_SUCCESS=0.999
	MESH=384;TRUNC=0.0
	LOG2P6=6*380.99700
	LMAX=225;BIN=5
def pairs_count(Pm):
	return 2.0*2**log2li(Pm.LX)
def node_cost(Pm,l_lo,l_hi,side,n,iA,kA,kB,R,thin=None):
	pairs=pairs_count(Pm)
	cfmax=max(b[1]for b in Pm.QBINS);cgmax=max(b[2]for b in Pm.QBINS)
	S=log2(CELLS[R]);c=(l_hi+S)/12
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
	mu=CELLS[R]*pr*Pm.UDUP*Pm.YF*det*(exp(-miss)if thin else 1.0)
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
	miss_children=0.0 if thin else miss,thin=bool(thin),pr_log2=log2(pr),
	surv_log2=log2(survivors),own_log2=log2(own),nf=nf_max,ng=ng_max,det=det,T=T,
	sieve_log2=log2(sieve),survcost_log2=log2(survivors*surv_cost))
def bins_of(Pm):
	bins=[(int(Pm.LB),int(Pm.LB)+Pm.BIN+1)]
	l=int(Pm.LB)+Pm.BIN+1
	while l<Pm.LMAX:
		bins.append((l,l+Pm.BIN));l+=Pm.BIN
	return bins
def p_short(Pm):
	"DESCENT.md (a); disc rows: b, side, n, shortfalls, worst j, w"
	return Pm.CD*log(20)*max(r[5]/r[2]for r in _D['disc'])
def tree_nodes():
	return 2*(1+sum(r[5]for r in _D['disc']))
def dom_c(u):
	"DESCENT.md, Encountered law"
	return max(max(w)/(sum(w)/len(w))for lo,hi in bins_of(P)for w in[[rho_series(u-b/P.LB)/b for b in[lo+(hi-lo)*i/64 for i in range(65)]]])
def disc_ok():
	return{(r[0],r[1])for r in _D['disc']}==set(PLAN16)and all(r[2]>=64 and r[3]==0 and r[5]>=0 for r in _D['disc'])
def _load(Pm,*a):
	"tabulated child load ([7t]); CHECK: indices recomputed"
	i=Pm._i;Pm._i+=1;t=Pm.TAB[i]if Pm.TAB else None
	if t is None or Pm.CHECK is True or i in Pm.CHECK:
		f=Pm.LOAD*child_load(*a);Pm._chk.append((i,f,t));return f if t is None else t
	return t
def boot_pr(Pm,BOOT):
	n1,iA,kA,Rb=BOOT;Sb=log2(CELLS[Rb]);cb=(Pm.LOG2P6+Sb)/12
	h=(n1-Pm.LX)/Pm.MESH;nA=Pm.LX+iA*h
	m_max=12*cb+Pm.LOG2U6+max(x[1]for x in Pm.BOOT_BINS)
	ss=lambda cg:ss_lower2(12*cb+Pm.LOG2U6+cg,Pm.LX,float(n1),iA,kA,1,Pm.MESH,Pm.ALPHA_G,Pm.TRUNC)
	pr=min(sum(ss(cg)for cg in Pm.BOOT_SUB)/10,Pm.PRCAP_BOOT*sum(w*ss(cg)for w,cg in Pm.BOOT_BINS))
	return Sb,cb,h,nA,m_max,pr
def evaluate(Pm,PLAN,BOOT):
	shared=('LX','MESH','TRUNC','ALPHA_F','ALPHA_G')
	if any(getattr(Pm,name)!=getattr(P,name)for name in shared):
		raise ValueError('descent and retention parameters must agree')
	import retain16
	retain16.loss.cache_clear();Pm._i=0;Pm._chk=[]
	bins=bins_of(Pm)
	W={};PHI={};DET={}
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
			if n>Pm.LX:
				L=[[_load(Pm,math.ceil(d['n'+sd]),Pm.LX,n,iA,kA,kB,Pm.MESH,al,cell_values(bins,V,sd,Pm.LX,n,Pm.MESH))
				for V in(PHI,W)]for sd,al in(('f',Pm.ALPHA_F),('g',Pm.ALPHA_G))]
			else:L=[[0.0,0.0],[0.0,0.0]]
			d['load']=L
			psi=L[0][0]+L[1][0]+d['miss_children']*pm
			sub=L[0][1]+L[1][1]+d['miss_children']*wm
			W[(b,side)]=d['own']+sub
			PHI[(b,side)]=Pm.zeta(d['mu'],(b[0],side))+psi
			DET[(b,side)]=dict(d,psi=psi,retry=1.0,sub_log2=log2(sub)if sub>0 else None,
			W_log2=log2(W[(b,side)]),phi=PHI[(b,side)])
	n1,iA,kA,Rb=BOOT
	if n1>Pm.LMAX:
		raise ValueError('boot window above the last bin')
	Sb,cb,h,nA,m_max,pr=boot_pr(Pm,BOOT)
	det=1.0-(1-Pm.ECM_SUCCESS)*(m_max/20.0+kA+1)
	import retain16
	per_rel=1.0/(pr*det*Pm.UDUP*Pm.YF*retain16.counted_fraction())
	trial=Pm.NORMS+ecm_instr(nA,m_max,Pm.ECM_SUCCESS)+Pm.PRIME+Pm.IDEAL
	pm,wm=miss_terms();mb=missing(Pm,0,m_max)
	Lb=[_load(Pm,math.ceil(m_max),Pm.LX,float(n1),iA,kA,1,Pm.MESH,Pm.ALPHA_G,cell_values(bins,V,'g',Pm.LX,float(n1),Pm.MESH))for V in(PHI,W)]
	psi=Lb[0]+mb*pm
	trials=per_rel
	balls=1.0+trials/CELLS[Rb]
	sub=Lb[1]+mb*wm
	one_log=(Pm.LLL_NODE+(Pm.LLL_NODE+CELLS[Rb]*Pm.NORMS)*Pm.BOOT_DISC/(1-Pm.BOOT_DISC))*balls+(trials+CELLS[Rb])*trial+sub
	boot=dict(n1=n1,nA=nA,kA=kA,R=Rb,S=Sb,m=m_max,pr_log2=log2(pr),trials_log2=log2(trials),
	trial_log2=log2(trial),balls=balls,miss=mb,own_log2=log2(trials*trial),sub_log2=log2(sub),
	psi=psi,retry=1.0,phi=psi,one_log_log2=log2(one_log))
	return dict(W=W,PHI=PHI,DET=DET,bins=bins,boot=boot,two_log_units=log2(2*one_log/UNIT),
	one_log_instr=one_log)
_PLAN={v[0]:tuple(map(tuple,v[1:]))for v in _D['plan']}
PLAN16={(b,s):(v[i][0],8*v[i][1],*v[i][2:])
for b,v in _PLAN.items()for i,s in enumerate('fg')}
BOOT16=(215,304,15,4)

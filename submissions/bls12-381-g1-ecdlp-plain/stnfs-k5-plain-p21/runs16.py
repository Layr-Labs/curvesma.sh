"""v10 canonical a-runs along e_0 (record path rec7, KERNELS.md). For a canonical a the run is the maximal integer interval
of t containing 0 with Q_a(eps (a + t e_0)) >= Q_a(a + t e_0) for every listed unit eps != +-1 (a quadratic in t per unit)
and Q_a(a + t e_0) <= theta A^2; the list is every u1^i u2^j u3^k (fundamental units of O_h) with log-radius <= 9.
Exact integer forms; float roots locate the ends, which are then verified exactly. Stdlib only."""
import cmath,math
import common12 as C
UF=[[0,-1,0,0,0,0],[1,0,-1,-1,1,-1],[-2,0,1,1,-1,1]]
K_EXP,LOG_RADIUS=34,9.0
def _roots():
	c=[complex(x)for x in C.H]
	z=[cmath.exp(2j*math.pi*k/6+0.4j)*1.1 for k in range(6)]
	ev=lambda x:sum(ck*x**k for k,ck in enumerate(c))
	for _ in range(500):
		z=[zi-ev(zi)/math.prod((zi-zj)for j,zj in enumerate(z)if j!=i)for i,zi in enumerate(z)]
	return z
def _inv(u):
	from fractions import Fraction as F
	rows=[C.mul([0]*j+[1],u)for j in range(6)]
	M=[[F(rows[j][i])for j in range(6)]+[F(int(i==0))]for i in range(6)]
	for k in range(6):
		p=next(i for i in range(k,6)if M[i][k]);M[k],M[p]=M[p],M[k]
		pv=M[k][k];M[k]=[x/pv for x in M[k]]
		for i in range(6):
			if i!=k and M[i][k]:
				f=M[i][k];M[i]=[x-f*y for x,y in zip(M[i],M[k])]
	inv=[int(M[i][6])for i in range(6)]
	assert C.mul(inv,u)==[1,0,0,0,0,0]
	return inv
def unit_list():
	r=_roots()
	lg=[[math.log(abs(sum(uk*x**k for k,uk in enumerate(u))))for x in r]for u in UF]
	pw=[]
	for u in UF:
		d={0:[1,0,0,0,0,0]};ui=_inv(u)
		for e in range(1,K_EXP+1):
			d[e]=C.mul(d[e-1],u);d[-e]=C.mul(d[-e+1],ui)
		pw.append(d)
	out=[]
	R=range(-K_EXP,K_EXP+1)
	for i in R:
		for j in R:
			for k in R:
				if max(abs(i*a+j*b+k*c)for a,b,c in zip(*lg))<=LOG_RADIUS:
					out.append(C.mul(C.mul(pw[0][i],pw[1][j]),pw[2][k]))
	return out
_QU={}
def run(a,units,q,lim):
	def form(x,y):return sum(x[i]*q[i][j]*y[j]for i in range(6)for j in range(6))
	E=[[int(i==j)for i in range(6)]for j in range(6)]
	cols=[C.mul(a,E[j])for j in range(6)]
	G=[[form(cols[i],cols[j])for j in range(6)]for i in range(6)]
	Bd=[[form(cols[i],E[j])for j in range(6)]for i in range(6)]
	qa0,bad0,qd0=form(a,a),form(a,E[0]),q[0][0]
	if qa0>lim:return 0,-1 # outside the region
	if id(units)not in _QU:_QU[id(units)]=[form(u,u)for u in units]
	coef=[]
	for u,quu in zip(units,_QU[id(units)]):
		if u in([1,0,0,0,0,0],[-1,0,0,0,0,0]):continue
		qu=sum(u[i]*G[i][j]*u[j]for i in range(6)for j in range(6))
		bu=sum(u[i]*Bd[i][j]*u[j]for i in range(6)for j in range(6))
		coef.append((qu-qa0,2*(bu-bad0),quu-qd0))
	if min(c[0]for c in coef)<0:return 0,0
	def ok(t):
		x=[a[0]+t]+a[1:]
		return form(x,x)<=lim and all(c0+c1*t+c2*t*t>=0 for c0,c1,c2 in coef)
	def first(sgn):
		best=math.inf
		for al,be,ga in coef:
			B_=be*sgn
			if al<0:return 0.0
			if ga==0:
				if B_<0:best=min(best,-al/B_)
				continue
			d=B_*B_-4*ga*al
			if d<0:continue
			s=math.sqrt(d)
			for rt in((-B_-s)/(2*ga),(-B_+s)/(2*ga)):
				if rt>0:best=min(best,rt)
		d=bad0*bad0-qd0*(qa0-lim);s=math.sqrt(max(d,0))
		return min(best,((-bad0+s)if sgn>0 else(bad0+s))/qd0)
	ends=[]
	for sgn in(1,-1):
		t=max(0,math.floor(first(sgn)+1e-9))
		while t>0 and not ok(sgn*t):t-=1
		while ok(sgn*(t+1)):t+=1
		ends.append(sgn*t)
	return ends[1],ends[0]
def lrs(a,t_lo,t_hi,q,den,A2,lam=1,kap=1):
	"""log2 sum over the run of (lam - kap qa)^3, exact"""
	from fractions import Fraction as F
	lam,kap=F(lam),F(kap);d=lam.denominator*kap.denominator;c=int(lam*d)*den*A2;k=int(kap*d)
	qa0=sum(a[i]*q[i][j]*a[j]for i in range(6)for j in range(6));b0=sum(a[i]*q[i][0]for i in range(6));q00=q[0][0]
	s=sum((c-k*(qa0+2*b0*t+q00*t*t))**3 for t in range(t_lo,t_hi+1))
	return math.log2(s)-3*math.log2(d*den*A2)

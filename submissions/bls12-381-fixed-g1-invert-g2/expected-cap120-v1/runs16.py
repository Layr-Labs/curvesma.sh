""
import cmath,math
import common12 as C
from fractions import Fraction as F
K_EXP,LOG_RADIUS=40,9.0
def _roots():
	c=[complex(x)for x in C.H];z=[cmath.exp(2j*math.pi*k/6+0.4j)*1.1 for k in range(6)];ev=lambda x:sum(ck*x**k for k,ck in enumerate(c))
	for _ in range(500):
		z=[zi-ev(zi)/math.prod((zi-zj)for j,zj in enumerate(z)if j!=i)for i,zi in enumerate(z)]
	return z
def _inv(u):
	rows=[C.mul([0]*j+[1],u)for j in range(6)];M=[[F(rows[j][i])for j in range(6)]+[F(int(i==0))]for i in range(6)]
	for k in range(6):
		p=next(i for i in range(k,6)if M[i][k]);M[k],M[p]=M[p],M[k];pv=M[k][k];M[k]=[x/pv for x in M[k]]
		for i in range(6):
			if i!=k and M[i][k]:
				f=M[i][k];M[i]=[x-f*y for x,y in zip(M[i],M[k])]
	inv=[int(M[i][6])for i in range(6)];assert C.mul(inv,u)==C.ONE;return inv
def _emb(u,z):
	p=C.power2(u);return sum(c*z**k for k,c in enumerate(p))/2
def unit_list():
	r=_roots();lg=[[math.log(abs(_emb(u,x)))for x in r]for u in C.UNITS];pw=[]
	for u in C.UNITS:
		d={0:C.ONE};ui=_inv(u)
		for e in range(1,K_EXP+1):
			d[e]=C.mul(d[e-1],u);d[-e]=C.mul(d[-e+1],ui)
		pw.append(d)
	out=[];R=range(-K_EXP,K_EXP+1)
	for i in R:
		for j in R:
			for k in R:
				if max(abs(i*a+j*b+k*c)for a,b,c in zip(*lg))<=LOG_RADIUS:
					assert max(abs(i),abs(j),abs(k))<K_EXP;out.append(C.mul(C.mul(pw[0][i],pw[1][j]),pw[2][k]))
	return out
_QU={}
def run(a,units,q,lim):
	def form(x,y):return sum(x[i]*q[i][j]*y[j]for i in range(6)for j in range(6))
	E=[[int(i==j)for i in range(6)]for j in range(6)];qa0,bad0,qd0=form(a,a),form(a,E[0]),q[0][0]
	if id(units)not in _QU:_QU[id(units)]=[form(u,u)for u in units]
	coef=[]
	for t in(a,C.sigma(a)):
		cols=[C.mul(t,E[j])for j in range(6)];G=[[form(cols[i],cols[j])for j in range(6)]for i in range(6)];Bd=[[form(cols[i],E[j])for j in range(6)]for i in range(6)]
		for u,quu in zip(units,_QU[id(units)]):
			if t is a and u in(C.ONE,[-x for x in C.ONE]):continue
			qu=sum(u[i]*G[i][j]*u[j]for i in range(6)for j in range(6));bu=sum(u[i]*Bd[i][j]*u[j]for i in range(6)for j in range(6));coef.append((qu-qa0,2*(bu-bad0),quu-qd0))
	if min(c[0]for c in coef)<0:return 0,0
	def ok(t):
		x=[a[0]+t]+a[1:];return form(x,x)<=lim and all(c0+c1*t+c2*t*t>=0 for c0,c1,c2 in coef)
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
		d=bad0*bad0-qd0*(qa0-lim);s=math.sqrt(max(d,0));return min(best,((-bad0+s)if sgn>0 else(bad0+s))/qd0)
	ends=[]
	for sgn in(1,-1):
		t=max(0,math.floor(first(sgn)+1e-9))
		while t>0 and not ok(sgn*t):t-=1
		while ok(sgn*(t+1)):t+=1
		ends.append(sgn*t)
	return ends[1],ends[0]
def lrs(a,t_lo,t_hi,q,den,A2,lam=1,kap=1):
	"";lam,kap=F(lam),F(kap);d=lam.denominator*kap.denominator;c=int(lam*d)*den*A2;k=int(kap*d);qa0=sum(a[i]*q[i][j]*a[j]for i in range(6)for j in range(6));b0=sum(a[i]*q[i][0]for i in range(6));q00=q[0][0];s=sum((c-k*(qa0+2*b0*t+q00*t*t))**3 for t in range(t_lo,t_hi+1));return math.log2(s)-3*math.log2(d*den*A2)
PI4=C.coords([-1,-5,2,-1,1,0])
def div(x,d):
	"";cols=[C.mul(d,[int(i==j)for j in range(6)])for i in range(6)];m=C.inverse([[cols[j][i]for j in range(6)]for i in range(6)]);q=[sum(m[i][j]*x[j]for j in range(6))for i in range(6)];return[int(v)for v in q]if all(F(v).denominator==1 for v in q)else None
def val(x,d,cap=99):
	k=0
	while k<cap and any(x):
		y=div(x,d)
		if y is None:
			break
		x,k=y,k+1
	return k
def facts():
	w,p4=C.W,PI4;one=C.ONE;two=[2,0,0,0,0,0];u2=div(two,C.mul(C.mul(w,w),C.mul(p4,p4)));ok=C.norm(w)==2 and C.norm(p4)==4 and div(p4,w)is None and u2 is not None and C.norm(u2)==1;reps=[[i,j,k,l,m,n]for i in range(2)for j in range(2)for k in range(2)for l in range(2)for m in range(2)for n in range(2)];al=[x for x in reps if div([a+b+c for a,b,c in zip(C.mul(x,x),x,one)],p4)is not None];ok&=len(al)>0;t4=[[0]*6,one,al[0],[a+b for a,b in zip(al[0],one)]];ff=lambda t:[sum(v)for v in zip(*[C.mul(c,C._pw(t,d))for d,c in C.TERMS])];ok&=all(div(ff(t),p4)is None for t in t4);vf=[val(c,w)for d,c in C.TERMS];vu2,vu4=val(C.UW,w),val(C.UW,p4);ok&=vf==[0,3,5,5,6]and vu2==33 and vu4==32;return ok,dict(vf=vf,vg=(vu2,vu4))
def expectations():
	"";ef=F(1,3)*6+2*0;eg2=F(2,3)*(sum(2*F(1,2**k)for k in range(1,17))+F(1,2**17));eg4=F(4,5)*sum(2*F(1,4**k)for k in range(1,17));return ef,eg2+2*eg4
def alpha2():
	ef,eg=expectations();return float(1-ef),float(1-eg)

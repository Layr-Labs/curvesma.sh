"";import sys as _s,os as _o
class _L:
	"""lib16.b85 modules (store16, lattice_count_bounds, unit_orbit, owner10; bytes)"""
	def find_spec(s,n,p=None,t=None):
		import importlib.util as u,lzma,base64,re
		if not hasattr(s,'m'):s.d=_o.path.dirname(_o.path.abspath(__file__));x=re.split(r"(?m)^@@ (\S+)\.py\n",lzma.decompress(base64.b85decode(open(_o.path.join(s.d,'lib16.b85')).read().strip())).decode());s.m=dict(zip(x[1::2],x[2::2]))
		if n in s.m:return u.spec_from_loader(n,s)
	def create_module(s,spec):return None
	def exec_module(s,m):m.__file__=_o.path.join(s.d,m.__name__+'.py');exec(compile(s.m[m.__name__],m.__file__,'exec'),m.__dict__)
_s.meta_path.append(_L())
from fractions import Fraction as F;from math import isqrt,ceil,lcm;H=[-1,0,-2,-2,2,0,1];U=15132376222941642752;WP=[1,1,0,0,0,0];P2F=[1,1,0,1,1,0];SIGP=[0,2,2,-2,0,-1];C2P=[-3,0,-3,3,0,1];BAS2=[[2,0,0,0,0,0],[-1,3,-2,1,-1,0],[1,-5,2,-1,1,0],[0,-4,-4,4,0,2],[3,0,3,-3,0,-1],[0,1,-5,2,-1,1]];UNITS2=[[0,-2,0,0,0,0],[1,4,3,-1,0,-1],[-2,-5,-1,2,1,1]]
def pmul(a,b):
	c=[0]*(len(a)+len(b)-1)
	for i,x in enumerate(a):
		for j,y in enumerate(b):c[i+j]+=x*y
	for i in range(len(c)-1,5,-1):
		for j in range(6):c[i-6+j]-=c[i]*H[j]
	return(c+[0]*6)[:6]
def inverse(matrix):
	n=len(matrix)
	a=[[F(x)for x in row]+[F(i==j)for j in range(n)]for i,row in enumerate(matrix)]
	for k in range(n):
		j=next(j for j in range(k,n)if a[j][k]);a[k],a[j]=a[j],a[k];p=a[k][k];a[k]=[x/p for x in a[k]]
		for j in range(n):
			if j!=k:p=a[j][k];a[j]=[x-p*y for x,y in zip(a[j],a[k])]
	return[row[n:]for row in a]
_BI=inverse(BAS2)
def coords(p2):"";v=[sum(F(p2[i])*_BI[i][j]for i in range(6))for j in range(6)];assert all(x.denominator==1 for x in v),'not in O_h';return[int(x)for x in v]
def power2(a):"";return[sum(a[i]*BAS2[i][j]for i in range(6))for j in range(6)]
T=[[coords([x//2 for x in pmul(BAS2[i],BAS2[j])])if all(x%2==0 for x in pmul(BAS2[i],BAS2[j]))else None for j in range(6)]for i in range(6)];assert all(t is not None for r in T for t in r);_TNZ=[[[(k,t)for k,t in enumerate(T[i][j])if t]for j in range(6)]for i in range(6)]
def mul(a,b):
	c=[0]*6
	for i,x in enumerate(a):
		if x:
			for j,y in enumerate(b):
				if y:
					for k,t in _TNZ[i][j]:c[k]+=x*y*t
	return c
ONE=[1,0,0,0,0,0]
def elt(p):"";return coords([2*x for x in(list(p)+[0]*6)[:6]])
W=elt(WP)
def _pw(x,n):
	r=ONE
	for _ in range(n):r=mul(r,x)
	return r
TERMS=[(12,ONE),(10,[-2*x for x in W]),(6,[2*x for x in _pw(W,3)]),(2,_pw(W,5)),(0,_pw(W,6))];PTERMS=[(12,[1]),(10,[-2*x for x in WP]),(6,[2*x for x in pmul(pmul(WP,WP),WP)]),(2,pmul(pmul(pmul(WP,WP),pmul(WP,WP)),WP)),(0,pmul(pmul(WP,WP),pmul(pmul(WP,WP),pmul(WP,WP))))];UW=[U*x for x in W];CU=coords(C2P);UNITS=[coords(u)for u in UNITS2]
def _sig_power2(p2):
	out=[0]*6;s=[1,0,0,0,0,0]
	for k in range(6):out=[x+p2[k]*y for x,y in zip(out,s)];s=pmul(s,SIGP)
	return out
S=[coords(_sig_power2(BAS2[i]))for i in range(6)]
def sigma(a):return[sum(a[i]*S[i][j]for i in range(6))for j in range(6)]
def norm(a):
	cols=[mul(a,[0]*i+[1]+[0]*(5-i))for i in range(6)];m=[[cols[j][i]for j in range(6)]for i in range(6)];prev,sign=1,1
	for k in range(5):
		if m[k][k]==0:
			row=next((j for j in range(k+1,6)if m[j][k]),None)
			if row is None:return 0
			m[k],m[row]=m[row],m[k];sign=-sign
		for i in range(k+1,6):
			for j in range(k+1,6):value=m[i][j]*m[k][k]-m[i][k]*m[k][j];assert value%prev==0;m[i][j]=value//prev
		prev=m[k][k]
	return abs(sign*m[-1][-1])
def norms(a,b):
	a2=mul(a,a);a4=mul(a2,a2);a6=mul(a4,a2);a10=mul(mul(a4,a4),a2);a12=mul(a6,a6);b2=mul(b,b);b4=mul(b2,b2);b6=mul(b4,b2);b10=mul(mul(b4,b4),b2);b12=mul(b6,b6);P={0:ONE,2:a2,6:a6,10:a10,12:a12};Q={0:ONE,2:b2,6:b6,10:b10,12:b12};f=[0]*6
	for d,c in TERMS:f=[x+y for x,y in zip(f,mul(c,mul(P[d],Q[12-d])))]
	g=[x+y for x,y in zip(a2,mul(UW,b2))];return norm(f),norm(g)
def integer_interval(a,t,c):
	a,t,c=F(a),F(t),F(c)
	if a<=0:raise ValueError('non-positive pivot')
	centre=-t/a;radius2=(t*t-a*c)/(a*a)
	if radius2<0:return range(0)
	n,d=centre.numerator,centre.denominator;z=isqrt(radius2.numerator*d*d//radius2.denominator);lo=-(-(n-z)//d);hi=(n+z)//d;return range(lo,hi+1)
def scaled_log(n,scale=16,upper=False):
	if n<1 or scale<1:raise ValueError('positive norm and scale required')
	x=n**scale;floor=x.bit_length()-1;return floor+int(upper and x!=1<<floor)
def integer_metric(matrix):den=lcm(*(F(x).denominator for row in matrix for x in row));q=[[int(F(x)*den)for x in row]for row in matrix];return q,den,inverse(q)
def quadratic(q,a):return sum(a[i]*q[i][j]*a[j]for i in range(6)for j in range(6))
def rho(u):import descent16;return 1. if u<=1 else descent16.rho_series(u)
_FR={}
def _frame(q):
	"";key=str(q)
	if key not in _FR:from collections import Counter;import ideal_canonical as IC;V=IC.lll([[int(i==j)for j in range(6)]for i in range(6)],[[F(x)for x in r]for r in q],Counter());assert abs(IC.determinant(V))==1;qq=[[sum(V[i][k]*q[k][l]*V[j][l]for k in range(6)for l in range(6))for j in range(6)]for i in range(6)];_FR[key]=V,inverse(qq)
	return _FR[key]
def uniform_half_ellipsoid(q,inverse_q,radius2,rng):
	"";V,iq=_frame(q);cap=[]
	for j in range(6):c=F(radius2)*iq[j][j];cap.append(isqrt(c.numerator//c.denominator))
	while True:
		y=[rng.randrange(-c,c+1)for c in cap];a=[sum(y[i]*V[i][j]for i in range(6))for j in range(6)]
		if a[5]>=1 and quadratic(q,a)<=radius2:return a,quadratic(q,a)
def regenerate_pair(cfg,seed,index):from random import Random;qa,da,ia=integer_metric(cfg['matrix_a']);qb,db,ib=integer_metric(cfg['matrix_b']);A2,theta=cfg['production_A_squared'],F(cfg['production_theta']);rng=Random(f"{cfg.get('sample_tag','v11-production-ideal-key')}:{seed}:{index}");a,aq=uniform_half_ellipsoid(qa,ia,int(A2*theta*da),rng);b,bq=uniform_half_ellipsoid(qb,ib,(A2-F(aq,da))*db,rng);return a,b,aq,bq
def stats():import json,lzma,base64,os;return json.loads(lzma.decompress(base64.b85decode(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'stats16.b85')).read())))

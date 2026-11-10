r"""Toy tower NFS vs the package's zero-retry first-choice attempt (tiny vs production). stdin {"seed",
"size","attempts"}+pb,lb,R,Rb,bw,top,KA. Base O=Z[s], s^2=2 (real quadratic, unit e=1+s, h=1). f=x^3+x+1,
m least >=2^(pb//3) with p=f(m) prime, p=3,5 mod 8; g=x-m. Relation (a0,a1,b0,b1), a=a0+a1s, b=b0+b1s:
Nf=|N(a^3+ab^2+b^3)|, Ng=|N(a-mb)|, N(x+ys)=x^2-2y^2. Node (q,u,t,side): q=+-1 mod 8, u^2=2, t a root of
the side polynomial mod q; rows (q,0,0,0),(-u,1,0,0),(t,0,1,0),(tu,0,0,1), W.lll; cells 0<|c|<=R (+-
once). Usable: gcd(N(a),N(b))=1; each norm (q removed) is 2^lb-smooth but for <=KA distinct split primes
in (2^lb,2^(lb+k bw)], k=(log2 q-lb)//bw: the children. Duplicates: (a,b), e^k(a,b) are one relation
(unit-balanced, sign-fixed representative). Attempt: z uniform in F_p*, boot rows (p,0),(z,1), cells
<=Rb, |a|,|b| semismooth over split primes in (2^lb,2^top] (children (l,u,m,g)), z redrawn until usable;
first usable relation everywhere, no retry; fails iff a node has none. Work: 1/cell, 8/usable cell.
Model (c28): S = sum over cells of prod ss_lower2(bits+alpha), alpha = sum_p (1/(p-1)-E_p)log2 p, E_p =
2n_p p/(p^2-1) split, 2n'_p p^2/(p^4-1) inert (n' roots in F_p^2), 2n_2/3 at 2; mu=YF S distinct/usable (D.P.YF);
z=exp(-mu/3); union of z over the tree."""
import sys,json,random
from math import gcd,log2,prod,exp,lgamma,log
from itertools import product
H=__import__('pathlib').Path(__file__).resolve().parent;sys.path.insert(0,str(H/'files'if(H/'files').is_dir()else H))
import descent_witness as W,descent16 as D
r=round;SIZES=[(48,11,4,5,4,24,6),(52,12,4,5,4,26,6),(56,13,4,5,4,28,6)]
def cp(k,n,lo=0,hi=1.):
	for _ in range(40*(k<n)):x=(lo+hi)/2;lo,hi=[(lo,x),(x,hi)][sum(exp(lgamma(n+1)-lgamma(i+1)-lgamma(n-i+1)+i*log(x)+(n-i)*log(1-x))for i in range(k+1))>.05]
	return r(hi,4)
def fac(n):
	if n<2:return[]
	if W.probable_prime(n):return[n]
	x=y=2;c=1
	while 1:
		x=(x*x+c)%n;y=(y*y+c)%n;y=(y*y+c)%n;d=gcd(x-y,n)
		if d==n:c+=1;x=y=2
		elif d>1:return fac(d)+fac(n//d)
ml=lambda a,b:(a[0]*b[0]+2*a[1]*b[1],a[0]*b[1]+a[1]*b[0])
nm=lambda a:abs(a[0]**2-2*a[1]**2)
sz=lambda x:x[0]**2+2*x[1]**2+x[2]**2+2*x[3]**2
def can(v):
	b=v
	for e in((1,1),(-1,1)):
		x=v
		while 1:
			y=ml(x[:2],e)+ml(x[2:],e)
			if sz(y)>=sz(x):break
			x=y
		if sz(x)<sz(b):b=x
	return b if next(c for c in b if c)>0 else tuple(-c for c in b)
def run(a):
	pb,lb,R,Rb,bw,top,KA=[a.get(k,v)for k,v in zip('pb lb R Rb bw top KA'.split(),SIZES[a.get('size',0)])]
	g=random.Random(a.get('seed',1));m=1<<pb//3
	while not(W.probable_prime(p:=m**3+m+1)and p%8 in(3,5)):m+=1
	B=1<<lb;pr=[x for x in range(2,B+1)if x<3 or W.probable_prime(x)];P=prod(pr);nb=(top-1-lb)//bw+1
	al=[sum((1/(x-1)-(n*2/3 if x<3 else 2*n*x/(x*x-1)if x%8 in(1,7)else 2*(n and 3 or len(c)<3)*x*x/(x**4-1)))*log2(x)for x in pr for n in[len(W.roots(c,x,g))])for c in([1,1,0,1],[-m,1])]
	bl=lambda R,d:[c for c in product(range(-R,R+1),repeat=d)if 0<sum(x*x for x in c)<=R*R and next(x for x in c if x)>0];C4,C2=bl(R,4),bl(Rb,2)
	M={};S={}
	def kids(N,w,sp):
		k=[];X=[]
		for x in N:
			d=gcd(x,P)
			while d>1:x//=d;d=gcd(x,d)
			if x>>w*KA:return
			X+=[x]
		for s,x in enumerate(X):
			f=fac(x)
			if len(f)>KA or len(set(f))<len(f)or any(not B<l<1<<w or sp and l%8 not in(1,7)for l in f):return
			k+=[(l,s)for l in f]
		return k
	def ch(l,s,a,b):
		for u in W.roots([-2,0,1],l,g):
			v=(b[0]+b[1]*u)%l
			if v and(((i:=(a[0]+a[1]*u)*pow(v,-1,l)%l)**3+i+1)*(s<1)+(i-m)*s)%l<1:return l,u,i,s
	def tree(F,w):
		o=F is not None
		for c in F or[]:
			x,y=nd(*c);o&=x;w+=y
		return o,w
	def nd(q,u,t,s):
		if q<=B:return 1,0
		I=q,u,t,s
		if I in M:return M[I]
		k=min((q.bit_length()-1-lb)//bw,nb-1);w=lb+k*bw;L,_=W.lll([[q,0,0,0],[-u%q,1,0,0],[t,0,1,0],[t*u%q,0,0,1]])
		T=n=0;F=None;Dd=set()
		for c in C4:
			v=tuple(sum(c[i]*L[i][j]for i in range(4))for j in range(4));A,Bb=v[:2],v[2:]
			Ef=ml(ml(A,A),A);Eb=ml(ml(A,Bb),Bb);Ec=ml(ml(Bb,Bb),Bb)
			N=[nm((Ef[0]+Eb[0]+Ec[0],Ef[1]+Eb[1]+Ec[1]))//q**(s<1),nm((A[0]-m*Bb[0],A[1]-m*Bb[1]))//q**s]
			if not all(N)or gcd(nm(A),nm(Bb))>1:continue
			T+=prod(D.ss_lower2(log2(x)+al[i],lb,w,64,KA,0,64,0,0)for i,x in enumerate(N));kk=kids(N,w,1)
			if kk is not None:
				n+=1;Dd.add(can(v))
				if F is None:F=[ch(l,x,A,Bb)for l,x in kk]
		S.setdefault((k,s),[]).append((n,len(Dd),T));M[I]=tree(F,len(C4)+8*n);return M[I]
	A_=[]
	for _ in range(a.get("attempts",40)):
		M.clear();F=None;w=0
		while F is None:
			(u,v),_=W.lll([[p,0],[g.randrange(1,p),1]])
			for i,j in C2:
				x=i*u[0]+j*v[0];y=i*u[1]+j*v[1];w+=1
				if gcd(x,y)<2 and x*y and(k:=kids([abs(x),abs(y)],top,1))is not None:w+=8;F=[(l,W.roots([-2,0,1],l,g)[0],m%l,1)for l,t in k];break
		A_+=[tree(F,w)]
	O={};Z=[0]*5;Y=[x for X in S.values()for x in X];ud=sum(x[1]for x in Y)/max(1,sum(x[0]for x in Y))
	for(k,s),X in sorted(S.items()):
		n=len(X);v=[sum(x[0]for x in X),sum(x[1]for x in X),sum(D.P.YF*ud*x[2]for x in X),sum(x[1]<1 for x in X),sum(exp(-D.P.YF*ud*x[2]/3)for x in X)];Z=[x+y for x,y in zip(Z,v)]
		O['%d%s'%(k,'fg'[s])]=dict(nodes=n,usable=v[0],distinct=v[1],mu=r(v[2],1),ratio=r(v[1]/v[2],3)if v[2]else None,zero=v[3],z_charged=r(v[4],3))
	n=len(A_);f=sum(1-o for o,_ in A_);wk=sorted(w for _,w in A_)
	return dict(p=p,udup=r(ud,4),attempts=n,fail=f,cp95=cp(f,n),union=r(Z[4]/n,4),nodes_per_tree=r(len(Y)/n,1),
	zero=Z[3],z_charged=r(Z[4],3),ratio=r(Z[1]/Z[2],3)if Z[2]else None,work=[r(wk[int(x*(n-1))])for x in(.5,.9,.99,1)],types=O)
if __name__=='__main__':print(json.dumps(run(json.loads(sys.stdin.read()or'{}')),separators=(',',':')))

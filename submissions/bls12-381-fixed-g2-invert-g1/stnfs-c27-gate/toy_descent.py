r"""Toy rational-NFS descent vs package laws; tiny vs production. stdin {"seed","size",
"attempts"}+pb,lb,R,Rb,bw,top. f=x^3+x+1, m least >=2^(pb//3) with p=f(m) prime, g=x-m. Attempt: z
uniform in F_p*; boot rows (p,0),(z,1), cells 0<|c|<=Rb (+- once, coprime); relation: |a|,|b|
2^lb-smooth but <=1 prime each in (2^lb,2^top], children (l,m mod l,g). Node (q,r,side): bin
k=(log2 q-lb)//bw, window (2^lb,2^(lb+k bw)], rows (q,0),(r,1), radius R; relation: |a^3+ab^2+b^3|,
|a-mb| (q removed) smooth but <=1 window prime each, children (l,a/b mod l,side). Cell order, <=32
tries; a relation succeeds iff its child subtrees do; outcome, work (cells) memoised. Model: S=sum
of prod ss_lower2(bits+alpha); psi,sub=child_load of phi,W at type max norms; phi=node mean
of max((1+2S(1-psi))^-.5+psi^32,.02); W=cells+sub/(1-psi). Attempt fails: dead end or work>c mean
model W (model phi+1/c); pair: consecutive fails, c=2.
Types: S=cert, nb/pois/phi = charged counts, retry [tries,fails] at relation j."""
import sys,json,random
from math import gcd,log2,prod,exp,lgamma,log,ceil
H=__import__('pathlib').Path(__file__).resolve().parent;sys.path.insert(0,str(H/'files'if(H/'files').is_dir()else H))
import descent_witness as W,descent16 as D
SIZES=[(36,10,4,3,2,19),(40,11,4,3,2,20),(44,11,4,3,2,21)]
r=round
def cp(k,n,lo=0,hi=1.):
	for _ in range(40*(k<n)):x=(lo+hi)/2;lo,hi=[(lo,x),(x,hi)][sum(exp(lgamma(n+1)-lgamma(i+1)-lgamma(n-i+1)+i*log(x)+(n-i)*log(1-x))for i in range(k+1))>.05]
	return r(hi,4)
def run(a):
	pb,lb,R,Rb,bw,top=[a.get(k,v)for k,v in zip('pb lb R Rb bw top'.split(),SIZES[a.get('size',0)])]
	g=random.Random(a.get('seed',1));m=1<<pb//3
	while not W.probable_prime(p:=m**3+m+1):m+=1
	B=1<<lb;pr=[x for x in range(2,B+1)if x<3 or W.probable_prime(x)];P=prod(pr);nb=(top-1-lb)//bw+1
	al=[sum((1-len(W.roots(c,x,g))*x/(x+1))*log2(x)/(x-1)for x in pr)for c in([1,1,0,1],[-m,1])]
	bl=lambda R:[(i,j)for i in range(-R,R+1)for j in range(-R,R+1)if 0<i*i+j*j<=R*R and(i>0 or i==0<j)]
	C,Cb=bl(R),bl(Rb);M={};S={};A=[];bins=[(x,x+bw)for x in range(lb,lb+nb*bw,bw)];Vp={};Vw={}
	def sc(L,Cs,w,Nf,sd):
		(u,v),_=W.lll(L);out=[];T=0;mx=[0,0]
		for i,j in Cs:
			a=i*u[0]+j*v[0];b=i*u[1]+j*v[1]
			if gcd(a,b)>1 or not all(N:=Nf(a,b)):continue
			y=[log2(x)+al[t]for x,t in zip(N,sd)];mx=list(map(max,mx,y));T+=prod(D.ss_lower2(z,lb,w,64,1,0,64,0,0)for z in y);c=[]
			for x,t in zip(N,sd):
				d=gcd(x,P)
				while d>1:x//=d;d=gcd(x,d)
				if x>1:
					if x>>w or not W.probable_prime(x):break
					c.append((x,t))
			else:out.append((a,b,c))
		return out,T,mx
	def tr(L,w):
		j=o=0
		for a,b,c in L[:32]:
			j+=1;o=1
			for x,t in c:
				y,v=d(x,m%x if t else a*pow(b,-1,x)%x,t);w+=v
				if not y:o=0;break
			if o:break
		return o,w,j
	def d(q,r,s):
		if q<=B:return 1,0
		I=q,r,s
		if I in M:return M[I]
		k=min((q.bit_length()-1-lb)//bw,nb-1);L,T,mx=sc([[q,0],[r,1]],C,lb+k*bw,lambda a,b:[abs(a**3+a*b*b+b**3)//q**(s<1),abs(a-m*b)//q**s],(0,1))
		o,w,j=tr(L,len(C));S.setdefault((k,s),[]).append((len(L),T,j-o,o,mx));M[I]=o,w;return o,w
	for _ in range(a.get('attempts',600)):
		L,T,mx=sc([[p,0],[g.randrange(1,p),1]],Cb,top,lambda a,b:[abs(a),abs(b)],(1,1));A.append((*tr(L,len(Cb))[:2],T,mx))
	def ld(mx,w,sd):
		r=[D.child_load(ceil(y),lb,w,64,1,0,64,0,D.cell_values(bins,V,t,lb,w,64))for V in(Vp,Vw)for y,t in zip(mx,sd)]
		return min(r[0]+r[1],1),r[2]+r[3]
	f=lambda T,ps:min(1,max((1+2*T*(1-ps))**-.5+ps**32,.02));O={}
	for k in range(nb):
		for s in(0,1):
			X=S.get((k,s));ps,sb=ld([*map(max,zip(*[x[4]for x in X]))],lb+k*bw,(0,1))if X else(1,0);K=bins[k],s
			Vp[K]=X and sum(f(x[1],ps)for x in X)/len(X)or 1.;Vw[K]=len(C)+sb/max(1-ps,1e-9)
			if X:
				n=len(X);c=list(zip(*X));v=[sum(c[0]),sum(c[1]),c[0].count(0),n-sum(c[3])]
				O['%d%s'%(k,'fg'[s])]=dict(nodes=n,rel=v[0],cert=r(v[1],2),ratio=r(v[0]/v[1],3),zeros=v[2],nb=r(sum((1+2*x[1])**-.5 for x in X),3),pois=r(sum(exp(-x[1])for x in X),3),
				fail=v[3],phi=r(Vp[K]*n,3),retry=[[sum(x[0]>=j and x[2]>=j-i for x in X)for i in(1,0)]for j in(1,2,3)])
	E=[(o,w,f(T,ps),len(Cb)+sb/max(1-ps,1e-9))for o,w,T,mx in A for ps,sb in[ld(mx,top,(1,1))]];n=len(E);Wm=sum(e[3]for e in E)/n;ph=sum(e[2]for e in E)/n
	rt=sorted(e[1]/Wm for e in E);F=[not o or w>2*Wm for o,w,*_ in E]
	att={c:dict(fail=k,cp95=cp(k,n),model=r(min(1,ph+1/c),4))for c in(2,4,8,16)for k in[sum(not o or w>c*Wm for o,w,*_ in E)]}
	return dict(p=p,attempts=n,phi=r(ph,4),W=r(Wm),cap=att,
	work=[r(rt[int(q*(n-1))],3)for q in(.5,.9,.99,1)]+[r(sum(rt)/n,3)],pair=[[sum((x,y)==(u,v)for x,y in zip(F,F[1:]))for v in(0,1)]for u in(0,1)],types=O)
if __name__=='__main__':print(json.dumps(run(json.loads(sys.stdin.read()or'{}')),separators=(',',':')))

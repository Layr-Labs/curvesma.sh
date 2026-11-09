"""H-X toy; construction and tested condition: proof.md section 5."""
from random import Random
import filtsim,model16 as M
def hx_toy(seed,lu=12,P=2**61-31,m=None,n=None):
	import unit_entries as UE
	m,n,g=m or M.BW_M,n or M.BW_N,Random(seed);c,_,R,_=filtsim.build_system(69.0,(UE.MUF/69,UE.MUG/69),M.FILT_EXCESS,seed,lu);UE.MID.clear()
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

"";import json,sys,os,random;_H=os.path.dirname(os.path.abspath(__file__));sys.path.insert(0,_H);import common12,filtsim as F
"""hxgen.py {"log2U":12..14,"seed":s,"n":11..24,"m":8..64}: H-X on this package's product-form live system (proof.md section 9): A = [M | S],
D = N + 10 live rows, S uniform, X as in production; at the charged d, e: rank [BY..B^dY] = N, rank O_e A = rank A, one-block pad control loses rank."""
Q=1000033;J=next(pow(g,(Q-1)//4,Q)for g in range(2,Q)if pow(g,(Q-1)//2,Q)==Q-1)
def cv(v):return(round(v.real)+round(v.imag)*J)%Q
def rank(rows):
	R=[r[:]for r in rows if any(r)];k=0;nc=len(R[0])if R else 0
	for c in range(nc):
		p=next((i for i in range(k,len(R))if R[i][c]),None)
		if p is None:continue
		R[k],R[p]=R[p],R[k];iv=pow(R[k][c],Q-2,Q);R[k]=[x*iv%Q for x in R[k]]
		for i in range(len(R)):
			if i!=k and R[i][c]:f=R[i][c];R[i]=[(a-f*b)%Q for a,b in zip(R[i],R[k])]
		k+=1
	return k
def main():
	q=json.loads(sys.stdin.read()or'{}');lu=int(q.get('log2U',13));sd=int(q.get('seed',1));n=int(q.get('n',13));m=int(q.get('m',32));assert 12<=lu<=14 and 11<=n<=24 and 8<=m<=64
	L=70.5;cols,w,rows,_=F.build_system(L,(740.3/L,498.7/L),1.09,sd,lu);Us=sum(1 for s,t in cols if t>=.8);rows2,_=F.purge(rows,len(cols),max(12,int(F.PEXC*Us)),False)
	H=type('H',(),dict(comb=lambda*a:0,stats=lambda*a:0,LAST=[]))();F.merge(rows2,len(cols),(200,),kmax=64,verbose=False,hook=H);R=[r for r in H.LAST if r]
	dc,lab=F.pivots(R,range(2**lu//2*2,len(cols)));live=[i for i in range(len(R))if i not in dc];lc=sorted({c for r in R for c in r if c>=0}-set(lab));N=len(lc);D=N+10
	order=sorted(dc,key=lambda p:lab[dc[p]]);ix={c:k for k,c in enumerate(lc)};lr=live[:D]
	if len(lr)<D:print(json.dumps(dict(log2U=lu,seed=sd,N=N,live_rows=len(live),ok=None,note='fewer live rows than N + 10')));return
	def op(x):
		v={c:x[ix[c]]for c in lc}
		for p in order:c=dc[p];v[c]=-cv(R[p][c].conjugate())*sum(cv(R[p][y])*v[y]for y in R[p]if y>=0 and y!=c)%Q
		return[sum(cv(R[i][y])*v[y]for y in R[i]if y>=0)%Q for i in lr]
	Mc=[op([int(j==k)for j in range(N)])for k in range(N)];g=random.Random(sd*7919+n*31+m)
	S=[[g.randrange(Q)for _ in range(D)]for _ in range(10)];A=[[Mc[k][i]for k in range(N)]+[S[t][i]for t in range(10)]for i in range(D)]
	Bv=lambda x:[sum(Mc[k][i]*x[k]for k in range(N))%Q for i in range(D)]
	d=-(-D//n)+1+-(-D//(128*n));e=-(-D//m)+1+-(-D//(128*m))
	Y=[S[t]for t in range(10)]+[Bv([g.randrange(Q)for _ in range(N)])for _ in range(n-10)];K=[];inc=[]
	cur=[Bv(y[:N])for y in Y]
	for s in range(d):K+=cur;inc.append(rank(K)if s<3 or s==d-1 else None);cur=[Bv(y[:N])for y in cur]
	rC=inc[-1]
	def Xof(onepad):
		idx=list(range(N));g.shuffle(idx);blk=[[]for _ in range(m)]
		for j,c in enumerate(idx):blk[j%m].append(c)
		for t in range(10):blk[0 if onepad else t%m].append(N+t)
		return[[(g.randrange(Q)if c in set(b)else 0)for c in range(D)]for b in blk]
	def obs(X):
		Ok=[];Xi=[r[:]for r in X]
		for _ in range(e):Ok+=Xi;Xi=[[sum(r[i]*A[i][k]for i in range(D))%Q if k<N else 0 for k in range(D)]for r in Xi]
		return rank([[sum(r[i]*A[i][k]for i in range(D))%Q for k in range(D)]for r in Ok])
	rO,rN=obs(Xof(False)),obs(Xof(True));rA=rank(A)
	print(json.dumps(dict(log2U=lu,seed=sd,N=N,D=D,n=n,m=m,d=d,e=e,rank_A=rA,rank_Cd=rC,first_ranks=inc[:3],rank_OeA=rO,control_one_block=rN,ok=rC==N and rO==rA and rN<rA)))
if __name__=="__main__":main()

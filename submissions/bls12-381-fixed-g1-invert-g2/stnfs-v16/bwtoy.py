import random
def mm(A,B,p):return[[sum(a*b for a,b in zip(r,c))%p for c in zip(*B)]for r in A]
def rk(M,p):
	M=[r[:]for r in M];R,C=len(M),len(M[0]);i=0
	for c in range(C):
		k=next((j for j in range(i,R)if M[j][c]),None)
		if k is None:continue
		M[i],M[k]=M[k],M[i];v=pow(M[i][c],p-2,p);M[i]=[x*v%p for x in M[i]]
		for j in range(R):
			if j!=i and M[j][c]:f=M[j][c];M[j]=[(a-f*b)%p for a,b in zip(M[j],M[i])]
		i+=1
		if i==R:break
	return i
def mbasis(H,dl,dr,m,n,p):
	K=n+m;P=[[([1]if r==c else[])for r in range(K)]for c in range(K)];sh=[0]*n+[1]*m
	for k in range(dl+dr+1):
		R=[[0]*K for _ in range(m)]
		for c in range(K):
			for u in range(m):
				s=sum(H[k-d][u][v]*cf for v in range(n)for d,cf in enumerate(P[c][v])if cf and 0<=k-d<len(H))
				pu=P[c][n+u];R[u][c]=(s-(pu[k]if k<len(pu)else 0))%p
		piv=[];used=set()
		for u in range(m):
			cand=[c for c in range(K)if c not in used and R[u][c]]
			if not cand:continue
			c0=min(cand,key=lambda c:(sh[c],c));used.add(c0);piv.append(c0);iv=pow(R[u][c0],p-2,p)
			for c in cand:
				if c==c0:continue
				f=R[u][c]*iv%p
				for r in range(m):R[r][c]=(R[r][c]-f*R[r][c0])%p
				for r in range(K):
					p1,p0=P[c][r],P[c0][r]
					if len(p1)<len(p0):p1.extend([0]*(len(p0)-len(p1)))
					for idx,val in enumerate(p0):p1[idx]=(p1[idx]-f*val)%p
		for c0 in piv:
			sh[c0]+=1
			for r in range(K):
				if P[c0][r]:P[c0][r]=[0]+P[c0][r]
	G=[]
	for c in sorted(range(K),key=lambda c:(sh[c],c))[:n]:
		dc=sh[c]
		if dc<=dr+1:
			vec=[(P[c][v][dc-j]if 0<=dc-j<len(P[c][v])else 0)%p for j in range(dr+2)for v in range(n)]
			if any(vec):G.append(vec)
	return G
def fam(k,N,p,g):
	if k=='sparse':
		while 1:
			A=[[g.randrange(p)if g.random()<.3 else 0 for _ in range(N)]for _ in range(N)];c=[g.randrange(p)for _ in range(N-1)]
			for r in A:r[-1]=sum(x*y for x,y in zip(c,r))%p
			if rk(A,p)==N-1:return A
	A=[[0]*N for _ in range(N)]
	for i in range(1,N):A[i][i]=(2 if i<=N//2 else 3)if k=='factors' else g.randrange(1,p)
	if k=='jordan':A[0][1],A[1][1]=1,0
	for _ in range(3*N):
		i,j=g.sample(range(N),2);c=g.randrange(p);A[i]=[(a+c*b)%p for a,b in zip(A[i],A[j])]
		for r in A:r[j]=(r[j]-c*r[i])%p
	return A
def trial(A,N,m,n,s,p,g):
	dl,dr=-(-N//m),-(-N//n);J=range(dr+2);X=[[g.randrange(s)for _ in range(N)]for _ in range(m)];Y=[[g.randrange(p)for _ in range(n)]for _ in range(N)];B=[mm(A,Y,p)]
	for _ in range(dl+dr):B.append(mm(A,B[-1],p))
	H=[mm(X,b,p)for b in B];R=[[B[j][i][v]for j in J for v in range(n)]for i in range(N)];T=[[H[i+j][u][v]for j in J for v in range(n)]for i in range(dl)for u in range(m)]
	rT=rk(T,p);rR=rk(R,p);K=mbasis(H,dl,dr,m,n,p)
	if any(sum(a*b for a,b in zip(r,col))%p for col in K for r in T):return rT,rR,0,1
	c=[g.randrange(p)for _ in K];x=[sum(a*b[t]for a,b in zip(c,K))%p for t in range(n*(dr+2))]
	if not any(x[:n]):return rT,rR,0,0
	w=[0]*N;Q=Y
	for k in J:w=[(w[i]+sum(Q[i][v]*x[k*n+v]for v in range(n)))%p for i in range(N)];Q=mm(A,Q,p)
	z=any(sum(a*b for a,b in zip(r,w))%p for r in A);return rT,rR,(not z)and any(w),z
def run(seed=1,t=40,N=16,m=4,n=3):
	g=random.Random(seed);rows=[];ok=1
	for p,s in((101,32),(7,7)):
		for k in('sparse','factors','jordan'):
			A=fam(k,N,p,g);f=bad=0
			for _ in range(t):rT,rR,good,z=trial(A,N,m,n,s,p,g);f+=not good;bad+=rT==rR and z
			b=(N-1)/s+(N+n+1)/p;ok&=bad==0 and(b>=1 or f/t<=b);rows.append((k,s,p,f,t,round(b,3),bad))
	return bool(ok),rows

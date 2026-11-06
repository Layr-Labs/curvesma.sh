import math,random
OPS=0
def mm(A,B,p):return[[sum(a*b for a,b in zip(r,c))%p for c in zip(*B)]for r in A]
def padd(a,b,p):return[((a[i]if i<len(a)else 0)+(b[i]if i<len(b)else 0))%p for i in range(max(len(a),len(b)))]
def kmul(a,b,p):
	global OPS
	if not a or not b:return[]
	n=max(len(a),len(b))
	if n<=2:
		o=[0]*(len(a)+len(b)-1)
		for i,x in enumerate(a):
			for j,y in enumerate(b):o[i+j]=(o[i+j]+x*y)%p;OPS+=1
		return o
	h=n//2;a0,a1=a[:h],a[h:];b0,b1=b[:h],b[h:]
	z0=kmul(a0,b0,p);z2=kmul(a1,b1,p);z1=kmul(padd(a0,a1,p),padd(b0,b1,p),p)
	mid=padd(z1,[(-x)%p for x in padd(z0,z2,p)],p)
	return padd(padd(z0,[0]*h+mid,p),[0]*(2*h)+z2,p)
def pmm(A,B,p):
	K=len(A);O=[[[]for _ in range(K)]for _ in range(K)]
	for c in range(K):
		for r in range(K):
			s=[]
			for k in range(K):s=padd(s,kmul(A[k][r],B[c][k],p),p)
			O[c][r]=s
	return O
def pmbasis(F,L,sh,m,p):
	global OPS;K=len(sh)
	if L==1:
		R=[r[:]for r in F[0]];P=[[([1]if r==c else[0])for r in range(K)]for c in range(K)];sh=sh[:];piv=[];used=set()
		for u in range(m):
			cand=[c for c in range(K)if c not in used and R[u][c]]
			if not cand:continue
			c0=min(cand,key=lambda c:(sh[c],c));used.add(c0);piv.append(c0);iv=pow(R[u][c0],p-2,p);OPS+=1
			for c in cand:
				if c==c0:continue
				f=R[u][c]*iv%p;OPS+=1
				for r in range(m):R[r][c]=(R[r][c]-f*R[r][c0])%p;OPS+=1
				for r in range(K):P[c][r][0]=(P[c][r][0]-f*P[c0][r][0])%p;OPS+=1
		for c0 in piv:
			sh[c0]+=1
			for r in range(K):P[c0][r]=[0]+P[c0][r]
		return P,sh
	h=L//2;P1,sh=pmbasis(F[:h],h,sh,m,p)
	Fm=[[[0]*K for _ in range(m)]for _ in range(L-h)]
	for c in range(K):
		for u in range(m):
			s=[]
			for k in range(K):s=padd(s,kmul([F[t][u][k]for t in range(len(F))],P1[c][k],p),p)
			for t in range(L-h):Fm[t][u][c]=s[h+t]if h+t<len(s)else 0
	P2,sh=pmbasis(Fm,L-h,sh,m,p)
	return pmm(P1,P2,p),sh
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
	global OPS;OPS=0;dl,dr=-(-N//m),-(-N//n);J=range(dr+2);Lg=dl+dr+1
	X=[[g.randrange(s)for _ in range(N)]for _ in range(m)];Y=[[g.randrange(p)for _ in range(n)]for _ in range(N)];B=[mm(A,Y,p)]
	for _ in range(dl+dr):B.append(mm(A,B[-1],p))
	H=[mm(X,b,p)for b in B];R=[[B[j][i][v]for j in J for v in range(n)]for i in range(N)];T=[[H[i+j][u][v]for j in J for v in range(n)]for i in range(dl)for u in range(m)]
	rT=rk(T,p);rR=rk(R,p);F=[[H[t][u]+[(-1 if t==0 and v==u else 0)%p for v in range(m)]for u in range(m)]for t in range(Lg)]
	P,sh=pmbasis(F,Lg,[0]*n+[1]*m,m,p);bnd=4*(m+n)**3*(Lg**math.log2(3))*math.log2(Lg)
	K=[(c,sh[c],[(P[c][v][sh[c]-j]if 0<=sh[c]-j<len(P[c][v])else 0)%p for j in J for v in range(n)])for c in sorted(range(n+m),key=lambda c:(sh[c],c))[:n]if sh[c]<=dr+1]
	K=[v for _,_,v in K if any(v)]
	if OPS>bnd or any(sum(a*b for a,b in zip(r,col))%p for col in K for r in T):return rT,rR,0,1
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

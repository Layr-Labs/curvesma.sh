import math,random,json,sys,time
from fractions import Fraction as F
Q32,Q38,M64=1<<32,1<<38,(1<<64)-1
def lut4(r2):
	B=r2-1
	if B==0:return math.isqrt(1<<52)+1
	A=64-B.bit_length();C=B<<A>>53;D=(C+1<<53-A)-1 if A<=53 else C>>A-53;return math.isqrt(D+1<<52)+1
def gs(G):
	C=len(G);D=[1]+[0]*C;E=[[0]*C for A in range(C)]
	for A in range(C):
		for B in range(A+1):
			F=G[A][B]
			for H in range(B):F=(D[H+1]*F-E[A][H]*E[B][H])//D[H]
			if B<A:E[A][B]=F
			else:D[A+1]=F
	return D,E
def tree(G,T,exact):
	B=12;A,P=gs(G);c=[[(P[C][B]<<32)//A[B+1]if C>B else 0 for B in range(B)]for C in range(B)];Q=[min(M64,(A[B+1]<<64)//(A[B]*T))for B in range(B)];d=[((A[B]*T<<38)+A[B+1]-1)//A[B+1]for B in range(B)];E=[0]*B
	for C in range(B-1,-1,-1):E[C]=math.isqrt((A[C]*T+A[C+1]-1)//A[C+1])+1+sum((E[A]+1)//2 for A in range(C+1,B))
	e=[sum(E[A]for A in range(A+1,B))+1 for A in range(B)];M=(A[2]<<38)//(A[1]*T);f=[exact and A[B+1]>=T*A[B]for B in range(B)];D=[0]*B;R=[];S=dict(exact_nodes=0)
	def L(k,rem,az):
		G=az;F=rem;H=-sum(c[A][k]*D[A]for A in range(k+1,B))
		if f[k]:
			S['exact_nodes']+=1;U=sum(P[A][k]*D[A]for A in range(k+1,B));V=F*T*A[k]*A[k+1];g=T*A[k]*A[k+1];N=-U//A[k+1]
			for C in range(max(N-1,0)if G else N-1,N+3):
				W=A[k+1]*C+U;X=Q38*W*W
				if X<=V:D[k]=C;L(k-1,-((X-V)//g),G and C==0)
			D[k]=0;return
		h=(F*d[k]>>64)+1;Y=lut4(h);I=e[k];E=(1 if k==0 else 0)if G else H-I-Y+Q32-1>>32;J=H+I+Y>>32
		if k==0:
			for C in range(E,J+1):D[0]=C;R.append(tuple(D))
		elif k==1:
			i=J-E+1;K=abs((E<<32)-H)>>13;O=F+(2*i*I*M>>32)+2+937*(J-E)+1-(K*K*Q[1]>>64);Z=((2*E+1)*Q32-2*H)*M>>32
			for C in range(E,J+1):
				if O>=0:D[1]=C;L(0,O,G and C==0)
				O-=Z;Z+=2*M
		else:
			for C in range(E,J+1):
				a=abs((C<<32)-H);K=(a-I if a>I else 0)>>13;b=K*K*Q[k]>>64
				if b<=F:D[k]=C;L(k-1,F-b,G and C==0)
		D[k]=0
	L(B-1,Q38,True);return R,S,A
def brute(G,T):
	A=12;C,D=gs(G);E=[F(C[A+1],C[A])for A in range(A)];K=[[F(D[B][A],C[A+1])for A in range(A)]for B in range(A)];H=[];B=[0]*A
	def I(k,rem):
		C=rem;D=-sum(K[A][k]*B[A]for A in range(k+1,A));F=math.isqrt(math.floor(C/E[k]))+1
		for G in range(math.floor(D)-F-1,math.ceil(D)+F+2):
			J=(G-D)**2*E[k]
			if J<=C:
				B[k]=G
				if k:I(k-1,C-J)
				else:H.append(tuple(B))
		B[k]=0
	I(A-1,F(T));return H
def canon(x):
	for A in x:
		if A:return x if A>0 else tuple(-A for A in x)
	return x
PROFS=[[1/30,1/25,1/28,1.3,1.,.8,2.5,1.9,1.5,1.15,.9,3.],[1/20,1/22,1/18,1.,.78,2.2,1.7,1.3,1.,4.,3.,2.3],[1/40,1/35,1/30,.9,1.6,1.25,.95,.75,1.8,1.4,1.05,2.6],[1/12,1/15,1/10,2.9,2.2,1.7,1.3,1.,.76,.6,.5,.45],[1/60,1/50,1/45,1.01,.76,.6,.5,1.2,.9,3.5,2.6,2.]]
def main(runs=16,seed=7):
	M=time.time();G=random.Random(seed);A=dict(runs=0,ref_points=0,exact_levels=0,exact_nodes=0,missing=0,extra=0,zero_or_dup=0,max_over=-1.,ctrl_runs_over_0024=0,ctrl_max_over=-1.)
	for N in range(runs):O=PROFS[N%len(PROFS)];B=G.randrange(1<<23,1<<25);E=[max(1,math.isqrt(int(A*B)))for A in O];H=[[E[B]if A==B else G.randint(-(E[A]//2),E[A]//2)if A<B else 0 for A in range(12)]for B in range(12)];D=[[sum(H[B][A]*H[C][A]for A in range(12))for C in range(12)]for B in range(12)];F=set(canon(A)for A in brute(D,B)if any(A));C,P,I=tree(D,B,True);J=set(C);K=lambda x:sum(D[A][B]*x[A]*x[B]for A in range(12)for B in range(12));A['runs']+=1;A['ref_points']+=len(F);A['exact_nodes']+=P['exact_nodes'];A['exact_levels']+=sum(I[A+1]>=B*I[A]for A in range(3,12));A['missing']+=len(F-set(canon(A)for A in J));A['extra']+=len(set(canon(A)for A in J)-F);A['zero_or_dup']+=sum(not any(A)for A in C)+len(C)-len(set(canon(A)for A in C));A['max_over']=max(A['max_over'],max((K(A)/B-1 for A in C),default=-1.));Q,R,R=tree(D,B,False);L=max((K(A)/B-1 for A in Q),default=-1.);A['ctrl_runs_over_0024']+=L>.0024;A['ctrl_max_over']=max(A['ctrl_max_over'],L)
	A['ok']=A['missing']==0 and A['zero_or_dup']==0 and A['max_over']<=.0024 and A['exact_levels']>0 and A['ctrl_runs_over_0024']>0;A['seconds']=round(time.time()-M,3);return A
if __name__=='__main__':o=main();print(json.dumps(o));sys.exit(0 if o['ok']else 1)
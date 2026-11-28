V='ctrl_max_over'
U='ctrl_runs_over_0024'
T='zero_or_dup'
S='missing'
R='exact_levels'
L='max_over'
Q='exact_nodes'
P=tuple
G=set
I=max
H=len
C=1.
E=sum
A=range
import math as D,random,json,sys,time
from fractions import Fraction as F
Q32,Q38,M=1<<32,1<<38,(1<<64)-1
def h(r2):
	B=r2-1
	if B==0:return D.isqrt(1<<52)+1
	A=64-B.bit_length();C=B<<A>>53;E=(C+1<<53-A)-1 if A<=53 else C>>A-53;return D.isqrt(E+1<<52)+1
def K(G):
	D=H(G);E=[1]+[0]*D;F=[[0]*D for a in A(D)]
	for B in A(D):
		for C in A(B+1):
			I=G[B][C]
			for J in A(C):I=(E[J+1]*I-F[B][J]*F[C][J])//E[J]
			if C<B:F[B][C]=I
			else:E[B+1]=I
	return E,F
def tree(G,T,exact):
	C=12;B,V=K(G);i=[[(V[q][v]<<32)//B[v+1]if q>v else 0 for v in A(C)]for q in A(C)];W=[min(M,(B[r+1]<<64)//(B[r]*T))for r in A(C)];j=[((B[s]*T<<38)+B[s+1]-1)//B[s+1]for s in A(C)];J=[0]*C
	for H in A(C-1,-1,-1):J[H]=D.isqrt((B[H]*T+B[H+1]-1)//B[H+1])+1+E((J[A]+1)//2 for A in A(H+1,C))
	l=[E(J[A]for A in A(t+1,C))+1 for t in A(C)];R=(B[2]<<38)//(B[1]*T);m=[exact and B[u+1]>=T*B[u]for u in A(C)];F=[0]*C;X=[];Y=dict(exact_nodes=0)
	def O(k,rem,az):
		J=az;H=rem;K=-E(i[A][k]*F[A]for A in A(k+1,C))
		if m[k]:
			Y[Q]+=1;Z=E(V[A][k]*F[A]for A in A(k+1,C));a=H*T*B[k]*B[k+1];n=T*B[k]*B[k+1];S=-Z//B[k+1]
			for D in A(I(S-1,0)if J else S-1,S+3):
				b=B[k+1]*D+Z;c=Q38*b*b
				if c<=a:F[k]=D;O(k-1,-((c-a)//n),J and D==0)
			F[k]=0;return
		o=(H*j[k]>>64)+1;d=h(o);L=l[k];G=(1 if k==0 else 0)if J else K-L-d+Q32-1>>32;M=K+L+d>>32
		if k==0:
			for D in A(G,M+1):F[0]=D;X.append(P(F))
		elif k==1:
			p=M-G+1;N=abs((G<<32)-K)>>13;U=H+(2*p*L*R>>32)+2+937*(M-G)+1-(N*N*W[1]>>64);e=((2*G+1)*Q32-2*K)*R>>32
			for D in A(G,M+1):
				if U>=0:F[1]=D;O(0,U,J and D==0)
				U-=e;e+=2*R
		else:
			for D in A(G,M+1):
				f=abs((D<<32)-K);N=(f-L if f>L else 0)>>13;g=N*N*W[k]>>64
				if g<=H:F[k]=D;O(k-1,H-g,J and D==0)
		F[k]=0
	O(C-1,Q38,True);return X,Y,B
def b(G,T):
	B=12;H,M=K(G);I=[F(H[a+1],H[a])for a in A(B)];N=[[F(M[c][d],H[d+1])for d in A(B)]for c in A(B)];J=[];C=[0]*B
	def L(k,rem):
		F=rem;G=-E(N[A][k]*C[A]for A in A(k+1,B));H=D.isqrt(D.floor(F/I[k]))+1
		for K in A(D.floor(G)-H-1,D.ceil(G)+H+2):
			M=(K-G)**2*I[k]
			if M<=F:
				C[k]=K
				if k:L(k-1,F-M)
				else:J.append(P(C))
		C[k]=0
	L(B-1,F(T));return J
def canon(x):
	for A in x:
		if A:return x if A>0 else P(-A for A in x)
	return x
O=[[1/30,1/25,1/28,1.3,C,.8,2.5,1.9,1.5,1.15,.9,3.],[1/20,1/22,1/18,C,.78,2.2,1.7,1.3,C,4.,3.,2.3],[1/40,1/35,1/30,.9,1.6,1.25,.95,.75,1.8,1.4,1.05,2.6],[1/12,1/15,1/10,2.9,2.2,1.7,1.3,C,.76,.6,.5,.45],[1/60,1/50,1/45,1.01,.76,.6,.5,1.2,.9,3.5,2.6,2.]]
def main(runs=16,seed=7):
	c=time.time();P=random.Random(seed);B=dict(runs=0,ref_points=0,exact_levels=0,exact_nodes=0,missing=0,extra=0,zero_or_dup=0,max_over=-C,ctrl_runs_over_0024=0,ctrl_max_over=-C)
	for d in A(runs):e=O[d%H(O)];F=P.randrange(1<<23,1<<25);M=[I(1,D.isqrt(int(i*F)))for i in e];W=[[M[j]if l==j else P.randint(-(M[l]//2),M[l]//2)if l<j else 0 for l in A(12)]for j in A(12)];K=[[E(W[k][A]*W[m][A]for A in A(12))for m in A(12)]for k in A(12)];N=G(canon(A)for A in b(K,F)if any(A));J,f,X=tree(K,F,True);Y=G(J);Z=lambda x:E(K[B][C]*x[B]*x[C]for B in A(12)for C in A(12));B['runs']+=1;B['ref_points']+=H(N);B[Q]+=f[Q];B[R]+=E(X[A+1]>=F*X[A]for A in A(3,12));B[S]+=H(N-G(canon(A)for A in Y));B['extra']+=H(G(canon(A)for A in Y)-N);B[T]+=E(not any(A)for A in J)+H(J)-H(G(canon(A)for A in J));B[L]=I(B[L],I((Z(A)/F-1 for A in J),default=-C));g,h,h=tree(K,F,False);a=I((Z(A)/F-1 for A in g),default=-C);B[U]+=a>.0024;B[V]=I(B[V],a)
	B['ok']=B[S]==0 and B[T]==0 and B[L]<=.0024 and B[R]>0 and B[U]>0;B['seconds']=round(time.time()-c,3);return B
if __name__=='__main__':B=main();print(json.dumps(B));sys.exit(0 if B['ok']else 1)
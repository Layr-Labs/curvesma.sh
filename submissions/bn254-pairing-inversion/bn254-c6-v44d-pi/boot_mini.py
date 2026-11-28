W='rec'
S=False
V=sorted
U=enumerate
R='r'
P=any
O=list
K=dict
J=.0
G=int
I=max
D=None
E=len
C=round
B=range
A=sum
import base64 as T,itertools as X,json,lzma,math as H
from math import log2
from pathlib import Path
import model_bn254_v6 as M,descent_bn254_v6 as DESC,boot_dp as L
LS=18,20,22
N=[D]
def d():
	if N[0]is D:
		C=1<<I(LS);F=bytearray(b'\x01')*(C+1);F[:2]=b'\x00\x00'
		for A in B(2,G(C**.5)+1):
			if F[A]:F[A*A::A]=b'\x00'*((C-A*A)//A+1)
		J=O(X.compress(B(C+1),F));K=[]
		for H in B(0,E(J),256):
			L=1
			for A in J[H:H+256]:L*=A
			K.append((L,J[H:H+256]))
		N[0]=K
	return N[0]
def tower_ok(u):
	C=36*u**4+36*u**3+24*u**2+6*u+1;H=M._pmul_h
	if C<5 or P(C%A==0 for A in B(2,G(C**.5)+1)):return S
	def I(a,e):
		A=[1,0,0,0,0,0]
		for B in bin(e)[2:]:
			A=H(A,A,C)
			if B=='1':A=H(A,a,C)
		return A
	F=[0,1,0,0,0,0]
	for O in B(3):
		F=I(F,C);K=[(F[b]-(1 if b==1 else 0))%C for b in B(6)];L=[c%C for c in M.H_COEFFS];A,D=L,K
		while P(D):
			while D and D[-1]==0:D=D[:-1]
			if not D:break
			while E(A)>=E(D)and P(A):
				N=A[-1]*pow(D[-1],-1,C)%C;J=E(A)-E(D);A=[(A[h]-N*D[h-J])%C if h>=J else A[h]for h in B(E(A))]
				while A and A[-1]==0:A=A[:-1]
			A,D=D,A
		if E(A)>1:return S
	return I([-((u+1)*g+(2 if f==0 else 0))%C for(f,g)in U(M.W_COEFFS)],(C**6-1)//2)==[C-1,0,0,0,0,0]
def record(j,u,sizes=D):
	Q=sizes;R=L.ctx(u);e=L.lattice_e(R,j+1);S=set();F=[]
	for K in L.cands(e):
		if I(abs(A)for A in K[:6]).bit_length()>=128 or I(abs(A)for A in K[6:]).bit_length()>=160:continue
		B=L.ng(R,K)
		if B==0:continue
		f=B not in S;S.add(B)
		if f and M.is_ideal_primitive_oh(K[:6],K[6:]):F.append(B)
	if Q is not D:Q+=[log2(i)for i in F]
	T=1
	for B in F:T*=B
	U=[]
	for(g,h)in d():
		V=H.gcd(T,g)
		if V>1:U+=[l for l in h if V%l==0]
	W=[];X=J
	for B in F:
		N=B;Y=1
		for G in U:
			if N%G==0:
				Y=G
				while N%G==0:
					N//=G
					if G<1<<16:X+=log2(G)
		W.append(Y if N==1 else D)
	Z=[]
	for a in LS:b=[M.rho((log2(k)+DESC.ALPHA_G)/a)for k in F];O=A(1 for A in W if A is not D and A<1<<a);P=A(b);c=(P*P-A(A*A for A in b))/2.;Z.append([O,O*(O-1)//2,C(100*log2(P)),C(100*log2(c))if c>0 else-10**9])
	return Z+[[E(F),C(10*X)]]
F=[D]
def Y(T):
	G=T['n'];F=E(LS);D=[[[0,0,0,0]for c in B(F)]+[[0,0]]for b in B(G)]
	for C in B(F):
		for(A,H)in T['M%d'%C]:D[A][C][0]=H
		for(A,H)in T['C%d'%C]:D[A][C][1]=H
	for A in B(G):
		I=T['e2'][A];D[A][F-1][2]=I
		for C in B(F-1):D[A][C][2]=I+T['d%d'%C][A]
		for C in B(F):D[A][C][3]=T['g%d'%C][A]+2*D[A][C][2]
		D[A][F]=[T['K'][A],T['s'][A]]
	return K(p=T['p'],r=T[R],n=G,rec=D)
def Z(raw):
	H,L,I=lzma.decompress(T.b85decode(raw)).partition(b'\n');E=json.loads(H);C=0
	for D in E.values():
		F=D.pop('_n')
		for(J,K,A)in D.pop('_c'):D[J]=[G.from_bytes(I[C+b:C+b+A],'little')+K for b in B(0,F*A,A)];C+=F*A
	return E
def load():
	if F[0]is D:
		A={}
		for G in V(Path(__file__).resolve().parent.glob('boot_mini*.b85')):
			for(B,C)in Z(G.read_bytes()).items():
				if B in A:
					for(H,E)in C.items():
						if isinstance(E,O):A[B][H]+=E
				else:A[B]=C
		F[0]={b:Y(c)for(b,c)in A.items()}
	return F[0]
def replay(u,js):F=load()[str(u)];B=[];D=[b for b in js if record(b,u,sizes=B)!=F[W][b]];return K(u=u,replayed=O(js),mismatches=D,kept_norm_bits_min_mean_max=[C(min(B),2),C(A(B)/E(B),2),C(I(B),2)],ok=not D)
def a(k):
	if k>300:return k-1.96*H.sqrt(k),k+1.96*H.sqrt(k)
	def A(n,m):
		A=C=H.exp(-m)
		for D in B(1,n+1):A*=m/D;C+=A
		return C
	def C(f,lo,hi):
		for C in B(80):
			A=(lo+hi)/2
			if f(A):hi=A
			else:lo=A
		return(lo+hi)/2
	D=J if k==0 else C(lambda m:1-A(k-1,m)>=.025,J,3.*k+50);return D,C(lambda m:A(k,m)<=.025,J,3.*k+50)
def Q(k,num,den,nl):B=A(den);D=A(num)/B;F,G=a(k);E=H.sqrt(A((A-D*B)**2 for(A,B)in zip(num,den))/(nl-1)/nl)/(B/nl);return[C(min(F/B,D-1.96*E),4),C(I(G/B,D+1.96*E),4)]
def stats():
	X={};f=A(log2(A)/(A-1)for A in M._primes_upto(1<<16))
	for(Y,I)in V(load().items(),key=lambda t:G(t[0])):
		B=I[W];H=E(B);Z={}
		for(F,g)in U(LS):L=[h[F][0]for h in B];N=[i[F][1]for i in B];O=[1 if j[F][0]>0 else 0 for j in B];S=[2.**(k[F][2]/1e2)for k in B];T=[2.**(l[F][3]/1e2)if l[F][3]>-10**8 else J for l in B];a=[m-n for(m,n)in zip(S,T)];b,P,c=A(S),A(T),A(a);Z[str(g)]=K(kept_smooth=A(L),model_sum=C(b,2),first_moment_ratio=C(A(L)/b,4),first_moment_ci95=Q(A(L),L,S,H),lattices_hit=A(O),bonferroni_indep=C(c,2),lattice_ratio=C(A(O)/c,4),lattice_ci95=Q(A(O),O,a,H),pairs=A(N),pairs_model=C(P,2),kappa_hat=C(A(N)/P,4)if P>0 else D,kappa_ci95=Q(A(N),N,T,H)if P>0 else D)
		d=A(A[3][0]for A in B);e=A(A[3][1]for A in B)/1e1/d;X[Y]=K(lattices=H,classes_of_e=(I[R]-1)//2,exhaustive=H==(I[R]-1)//2,tower_ok=tower_ok(G(Y)),p=I['p'],kept=d,mean_s16=C(e,4),alpha_emp=C(f-e,4),by_L=Z)
	return X
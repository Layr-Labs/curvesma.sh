import base64,itertools,json,lzma,math
from math import log2
from pathlib import Path
import model_bn254_v6 as M,descent_bn254_v6 as DESC,boot_dp as BD
TOWERS=-3,5
LS=18,20,22
_PR=[None]
def _chunks():
	if _PR[0]is None:
		B=1<<max(LS);C=bytearray(b'\x01')*(B+1);C[:2]=b'\x00\x00'
		for A in range(2,int(B**.5)+1):
			if C[A]:C[A*A::A]=b'\x00'*((B-A*A)//A+1)
		E=list(itertools.compress(range(B+1),C));F=[]
		for D in range(0,len(E),256):
			G=1
			for A in E[D:D+256]:G*=A
			F.append((G,E[D:D+256]))
		_PR[0]=F
	return _PR[0]
def tower_ok(u):
	B=36*u**4+36*u**3+24*u**2+6*u+1;E=M._pmul_h
	if B<5 or any(B%A==0 for A in range(2,int(B**.5)+1)):return False
	def F(a,e):
		A=[1,0,0,0,0,0]
		for C in bin(e)[2:]:
			A=E(A,A,B)
			if C=='1':A=E(A,a,B)
		return A
	D=[0,1,0,0,0,0]
	for K in range(3):
		D=F(D,B);H=[(D[A]-(1 if A==1 else 0))%B for A in range(6)];I=[A%B for A in M.H_COEFFS];A,C=I,H
		while any(C):
			while C and C[-1]==0:C=C[:-1]
			if not C:break
			while len(A)>=len(C)and any(A):
				J=A[-1]*pow(C[-1],-1,B)%B;G=len(A)-len(C);A=[(A[D]-J*C[D-G])%B if D>=G else A[D]for D in range(len(A))]
				while A and A[-1]==0:A=A[:-1]
			A,C=C,A
		if len(A)>1:return False
	return F([-((u+1)*C+(2 if A==0 else 0))%B for(A,C)in enumerate(M.W_COEFFS)],(B**6-1)//2)==[B-1,0,0,0,0,0]
def record(j,u,sizes=None):
	H=sizes;I=BD.ctx(u);V=BD.lattice_e(I,j+1);J=set();B=[]
	for D in BD.cands(V):
		if max(abs(A)for A in D[:6]).bit_length()>=128 or max(abs(A)for A in D[6:]).bit_length()>=160:continue
		A=BD.ng(I,D)
		if A==0:continue
		W=A not in J;J.add(A)
		if W and M.is_ideal_primitive_oh(D[:6],D[6:]):B.append(A)
	if H is not None:H+=[log2(A)for A in B]
	K=1
	for A in B:K*=A
	L=[]
	for(X,Y)in _chunks():
		N=math.gcd(K,X)
		if N>1:L+=[A for A in Y if N%A==0]
	O=[];P=.0
	for A in B:
		E=A;Q=1
		for C in L:
			if E%C==0:
				Q=C
				while E%C==0:
					E//=C
					if C<1<<16:P+=log2(C)
		O.append(Q if E==1 else None)
	R=[]
	for S in LS:T=[M.rho((log2(A)+DESC.ALPHA_G)/S)for A in B];F=sum(1 for A in O if A is not None and A<1<<S);G=sum(T);U=(G*G-sum(A*A for A in T))/2.;R.append([F,F*(F-1)//2,round(100*log2(G)),round(100*log2(U))if U>0 else-10**9])
	return R+[[len(B),round(10*P)]]
_S=[None]
def _dec(T):
	E=T['n'];D=len(LS);C=[[[0,0,0,0]for A in range(D)]+[[0,0]]for A in range(E)]
	for B in range(D):
		for(A,F)in T['M%d'%B]:C[A][B][0]=F
		for(A,F)in T['C%d'%B]:C[A][B][1]=F
	for A in range(E):
		G=T['e2'][A];C[A][D-1][2]=G
		for B in range(D-1):C[A][B][2]=G+T['d%d'%B][A]
		for B in range(D):C[A][B][3]=T['g%d'%B][A]+2*C[A][B][2]
		C[A][D]=[T['K'][A],T['s'][A]]
	return dict(p=T['p'],r=T['r'],n=E,rec=C)
def load():
	if _S[0]is None:
		A={}
		for E in sorted(Path(__file__).resolve().parent.glob('boot_mini*.b85')):
			for(B,C)in json.loads(lzma.decompress(base64.b85decode(E.read_bytes())).decode()).items():
				if B in A:
					for(F,D)in C.items():
						if isinstance(D,list):A[B][F]+=D
				else:A[B]=C
		_S[0]={A:_dec(B)for(A,B)in A.items()}
	return _S[0]
def replay(u,js):C=load()[str(u)];A=[];B=[B for B in js if record(B,u,sizes=A)!=C['rec'][B]];return dict(u=u,replayed=list(js),mismatches=B,kept_norm_bits_min_mean_max=[round(min(A),2),round(sum(A)/len(A),2),round(max(A),2)],ok=not B)
def _pois(k):
	if k>300:return k-1.96*math.sqrt(k),k+1.96*math.sqrt(k)
	def A(n,m):
		A=B=math.exp(-m)
		for C in range(1,n+1):A*=m/C;B+=A
		return B
	def B(f,lo,hi):
		for B in range(80):
			A=(lo+hi)/2
			if f(A):hi=A
			else:lo=A
		return(lo+hi)/2
	C=.0 if k==0 else B(lambda m:1-A(k-1,m)>=.025,.0,3.*k+50);return C,B(lambda m:A(k,m)<=.025,.0,3.*k+50)
def _wide(k,num,den,nl):A=sum(den);B=sum(num)/A;D,E=_pois(k);C=math.sqrt(sum((A-B*C)**2 for(A,C)in zip(num,den))/(nl-1)/nl)/(A/nl);return[round(min(D/A,B-1.96*C),4),round(max(E/A,B+1.96*C),4)]
def stats():
	K={};T=sum(log2(A)/(A-1)for A in M._primes_upto(1<<16))
	for(L,D)in sorted(load().items(),key=lambda t:int(t[0])):
		A=D['rec'];C=len(A);N={}
		for(B,U)in enumerate(LS):E=[A[B][0]for A in A];F=[A[B][1]for A in A];G=[1 if A[B][0]>0 else 0 for A in A];I=[2.**(A[B][2]/1e2)for A in A];J=[2.**(A[B][3]/1e2)if A[B][3]>-10**8 else .0 for A in A];O=[A-B for(A,B)in zip(I,J)];P,H,Q=sum(I),sum(J),sum(O);N[str(U)]=dict(kept_smooth=sum(E),model_sum=round(P,2),first_moment_ratio=round(sum(E)/P,4),first_moment_ci95=_wide(sum(E),E,I,C),lattices_hit=sum(G),bonferroni_indep=round(Q,2),lattice_ratio=round(sum(G)/Q,4),lattice_ci95=_wide(sum(G),G,O,C),pairs=sum(F),pairs_model=round(H,2),kappa_hat=round(sum(F)/H,4)if H>0 else None,kappa_ci95=_wide(sum(F),F,J,C)if H>0 else None)
		R=sum(A[3][0]for A in A);S=sum(A[3][1]for A in A)/1e1/R;K[L]=dict(lattices=C,classes_of_e=(D['r']-1)//2,exhaustive=C==(D['r']-1)//2,tower_ok=tower_ok(int(L)),p=D['p'],kept=R,mean_s16=round(S,4),alpha_emp=round(T-S,4),by_L=N)
	return K
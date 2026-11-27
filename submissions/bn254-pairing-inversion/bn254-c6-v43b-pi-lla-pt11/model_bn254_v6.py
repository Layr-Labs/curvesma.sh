Bs='mem_1own'
Br='log2_Ns_cap'
Bq='emitter_log2'
Bp='score_node_log2'
Bo='n_total'
Bn='nonprim_bins'
A8='bins_b62'
Bm='pool_quad'
Bl='descent'
A7=reversed
Ap='n_nonprim_in'
n='large'
m='sieve'
k='large_sf'
j='sieve_sf'
Ao='bins'
An=complex
AK='t4_geom'
AG=sorted
A1='n_prim_in'
A5='GA'
Y='large_geo'
X='sieve_geo'
z='EMB_ZK'
y=enumerate
S='small'
p=tuple
d='x4'
a=all
V=int
U=dict
Q=abs
O=True
N=round
M=float
K=False
J=pow
I=sum
H=min
G=None
E=max
D=.0
B=len
A=1.
C=range
import base64 as o,json,lzma as q,math as T,random,sys
from math import log,log2,exp,sqrt,ceil,floor,pi,log1p,gamma
from fractions import Fraction as F
from decimal import Decimal as e,localcontext as A9
from collections import Counter
from array import array
from pathlib import Path
import descent_bn254_v6 as DESC
L=T.log(2.)
D6=.5772156649015329
h=2e4
U_BN254=0x44e992b44a6909f1
P_BN254=36*U_BN254**4+36*U_BN254**3+24*U_BN254**2+6*U_BN254+1
R_BN254=36*U_BN254**4+36*U_BN254**3+18*U_BN254**2+6*U_BN254+1
H_COEFFS=[3,-4,0,2,1,-2,1]
W_COEFFS=[0,-2,-1,0,1,-1]
def AA(u=U_BN254):return[(u+1)*ab+(2 if aa==0 else 0)for(aa,ab)in y(W_COEFFS)]
G_LEAD_COEFFS=AA(U_BN254)
SIGMA_NUM=[9,-4,-5,1,3,-2]
ZK_NUM=[[3,0,0,0,0,0],[0,3,0,0,0,0],[9,-4,-5,1,3,-2],[0,2,1,-2,0,1],[6,1,-1,-1,3,-1],[-6,6,6,0,-3,3]]
F_REC_COEFFS={8:[1,0,0,0,0,0],6:[8,-10,-5,0,5,-5],4:[30,-12,-6,0,6,-6],2:[-4,14,7,0,-7,7],0:[-5,4,2,0,-2,2]}
CONTENT_BITS=D
D7=.683542
ALPHA_F_RAW=-.995
ALPHA_G_RAW=-1.026+.08
W=G
Z=G
_P17_G_CACHE=G
rand_below=DESC.rand_below
AB={ab:aa for(aa,ab)in y('0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ')}
def b62d(s):
	B=s.startswith('-');A=0
	for C in s[1:]if B else s:A=A*62+AB[C]
	return-A if B else A
def load_data():
	global W
	if W is G:
		J=Path(__file__).resolve().parent/'c6_data.b85';W=json.loads(q.decompress(o.b85decode(J.read_bytes())).decode());import boot_dp as D;B=D.load();B['off_hi']=E(B for A in B['rec']for B in A[3])/1e2+A;B['primes16']=list(_primes_upto(1<<16));B['ng_hi']=D.gate_norm_bits();W[Bl]['boot_dp']=B;C=W[Bm]
		if A8 in C and Ao not in C:
			F=0;H=[]
			for K in C[A8].split(';'):L,M,N,O=K.split();F+=b62d(L);I=b62d(N);H.append([F/1e2,(b62d(M)+42000)/1e2,I,I+b62d(O)])
			C[Ao]=H
	return W
def load_witnesses():
	global Z
	if Z is G:A=Path(__file__).resolve().parent/'c6_witnesses.b85';Z=json.loads(q.decompress(o.b85decode(A.read_bytes())).decode())
	return Z
C_HIT_LO=3.5704
C_B1_ROW=26.2
C_SB_ENTER_X4=42.6
D8=3.5
D9=7.05
C_B2_INT_X4=3.7
DA=6.55
C_B2_BND_HIT=7.75
C_B2_BND_MISS=5.75
C_B2_BND1_HIT=5.75
C_B2_BND1_MISS=3.75
C_SCAN_ZERO=.4
DB=8
C_SCAN_G8=272./1024.
C_LIST_B1=24.7
C_LIST_B2=40.75
Be=21.9
DC=4.75
DD=2.75
DE=.6
w=(4e1,42.5,8,d),(42.5,44.5,16,d),(44.5,46.,32,d)
B2_PREFIX_MULT=2.84
Al=D
C_LAT_EMIT=6.714
C_LAT_APPLY=1.53
C_LAT_LEAF=C_LAT_EMIT+C_LAT_APPLY
C_LAT_EMIT_V43B=9.3
PT_ON=O
PT_FB_EXC=.03
PT_GLOBAL_IR=2.**40
DF=D
AC=46
Bf=M(1<<AC+14)
Bg=34
l=M(1<<Bg)
Bh=38
DG=M(1<<Bh)
DH=64
SC_GUARD_DELTA=.005
DI=8
DJ=D
C_LAT_LLL=2.**34.5
C_LAT_EXC=2.**31.2
PT_CERT_IR=131072.
PT_BUILD_IR=2.**31.6
QA_MAX=.85
QB_MAX=.85
B_A_SHELLS=2
K_B1_PREF=256
B1_PREF_COV=1.0075
_WKF=512./729.
x=125./27.
Bi=8.
C_UNIT_5=6.71
C_GEN_1=42.
C_NORM_G=84e2
C_NORM_F=168e2
C_SRT_G_LEAF=46e3
C_SB=68e3
DK=2.**23.45
DL=2.**22.12
A0=1.001
DM=2.
AI,AL,AM,r,s,AD,AE=8.5,12e2,504.,67108864,2.,3.,32768.
NTT3=69.75,119.,9e1,151.5,371.,139.25,271.5
u=512.
AN,AO,AP,AQ,AR=70368744177664.,56,33024.,67108864.,576
DN,DO=112.,84.
W_ROW=100.5
NU_SW_CAP=.000207
NU_L_CAP=.000205
NC27=6.1e-05
DP=.999998
K5_LLA=1.25
TEXC_REL=.0001
TEXC_LA=6.
LLA_ON=O
K5_CAP=8.7
TAB_R_CAP=A
DQ=D
DR=112.
DS=28.
DT=97.
DU=4e1
Bj=27e1
Bk=.0078125
DV=2.**80
DW=2.**80
DX=2.**80
DY=2.**78
G2_OUTPUT_BOUND=69456896.
OP=U(lf=61.5,lg=61.5,exc=.72,W=100,p_min=24,p_cut=4e1,p_mid=46.,p_lat=46.,lx2=61.5,slack=117.,slack_b=G,w_sb=64,c_tail_z6=364.,margin=.005,work_allow=A)
R=G
def t(kmax=42,nt=280):
	with A9()as N:
		N.prec=220;H,O=e(0),e(1);F=[G,[O]+[H]*(nt-1)]
		for J in C(2,kmax+1):
			K,A=F[J-1],[H]*nt
			for D in C(nt-1):A[D+1]=(K[D]+D*A[D])/(J*(D+1))
			A[0]=K[0]-I(A[1:]);F.append(A)
		L=[G]
		for A in F[1:]:
			P=E(Q(A)for A in A)
			while B(A)>1 and Q(A[-1])<P*e('1e-18'):A=A[:-1]
			L.append(p(M(A)for A in A7(A)))
		return p(L)
def rho(u):
	global R
	if u<=A:return A
	if R is G:R=t()
	C=V(u)
	if u>C:C+=1
	if C>=B(R):return D
	F=C-u;E=D
	for H in R[C]:E=E*F+H
	return E if E>D else D
b=G
def rho_tab_1024():
	global R,b
	if b is not G:return b
	if R is G:R=t()
	F=R;B=[A]*1025
	for H in C(2,42):
		I=F[H]
		for J in C(1,1025):
			K=(1024-J)/1024.;E=D
			for L in I:E=E*K+L
			B.append(E if E>D else D)
	B.append(rho(41.+A/1024.));b=B;return B
def _rho_fast(u):
	if u<=A:return A
	D=rho_tab_1024();E=u*1024.;C=V(E)
	if C+1>=B(D):return rho(u)
	return D[C]+(E-C)*(D[C+1]-D[C])
def kappa_rho(u,t):A=rho(u);return rho(u-t)/A if A>D else D
def _primes_upto(n):
	B=bytearray(b'\x01')*(n+1);B[:2]=b'\x00\x00'
	for A in C(2,V(sqrt(n))+1):
		if B[A]:B[A*A::A]=b'\x00'*((n-A*A)//A+1)
	return[aa for aa in C(2,n+1)if B[aa]]
def v(weights,normalizer,length):
	H=weights;J=T.sumprod;A=length;K=[aa*H[aa]for aa in C(A,-1,-1)];G=[exp(-normalizer)]+[D]*A;B=1
	while B<=A and H[B]==D:B+=1
	F=A
	while F>=B and H[F]==D:F-=1
	if B>A:return G
	for E in C(B,A+1):I=0 if E<=F else E-F;G[E]=J(K[A-(E-I):A-B+1],G[I:E-B+1])/E
	return G
class Retention:
	def __init__(A,B_bits,X_bits,prime_min=2**14,step=.25,max_bits=600,tilt=.05,explicit_bits=14,power_terms=8,odd=1):
		W=tilt;U=prime_min;J=X_bits;G=B_bits;F=step;A.B,A.X,A.minimum,A.step=G,J,U,F;A.tilt,A.length=W,ceil(max_bits/F);A.names=S,X,Y,j,k;P={aa:[D]*(A.length+1)for aa in A.names};B=Counter();Q=H(M(explicit_bits),G)
		for c in _primes_upto(V(2**Q)-1)[odd:]:
			R=log2(c);d=N(R/F);T=exp(W*d*F-L*R);i=[(S,K)]if c<U else[((m if R<J else n)+'_geo',K),((m if R<J else n)+'_sf',O)]
			for(e,f)in i:
				B[e]+=log1p(T)if f else-log1p(-T);Z=T
				for a in C(1,A.length+1):
					g=a*d
					if g>A.length:break
					P[e][g]+=(-1 if f and a%2==0 else 1)*Z/a;Z*=T
					if Z<1e-25:break
		def b(left,right,pref):
			J=right;G=pref;E=left
			if E>=J:return
			H=F/4.;M=ceil((J-E)/H);H=(J-E)/M
			for U in C(M):
				K=E+(U+.5)*H;O=N(K/F)
				for D in C(1,power_terms+1):
					Q=D*O
					if Q>A.length:break
					I=exp(D*W*O*F-(D-1)*L*K)*H/(D*K)
					if I<1e-25:break
					V=-1 if D%2==0 else 1;X=((G,I),)if G==S else((G+'_geo',I),(G+'_sf',V*I))
					for(R,T)in X:B[R]+=T;P[R][Q]+=T
		h=log2(U);b(Q,H(h,G),S);b(E(Q,h),H(J,G),m);b(E(Q,J),G,n);A.weights,A.normalizers=P,{ac:B[ac]for ac in A.names};l=[I(P[B][ab]for B in(S,X,Y))for ab in C(A.length+1)];A.denominator=v(l,I(B[A]for A in(S,X,Y)),A.length);A.selection_normalizer=exp(B[j]+B[k]-B[X]-B[Y]);A._c={};A._ac={}
	def co(A,k):
		if k not in A._c:A._c[k]=v(A.weights[k],A.normalizers[k],A.length)
		return A._c[k]
	def two(A,n,a,b=9**9,win=G,sf=K):
		E=win;a=H(a,n)
		if a<0 or b<0:return D
		F=n,b,E,sf
		if F not in A._c:
			if F[1:]not in A._c:
				O=A.step;P=A.co(S);L=H(A.length,2048);I=[D]*(L+1);M=L
				while M>0 and P[M]<1e-22:M-=1
				for(N,Q)in y(A.co(k if sf else Y)[:H(b,L)+1]):
					if Q>1e-22 and(E is G or N*O<=E[2]or any(A.X*B-E[2]<=N*O<=E[1]*B+E[2]for B in C(1,E[0]+1))):
						for J in C(H(L-N,M)+1):I[N+J]+=Q*P[J]
				A._c[F[1:]]=I
			I=A._c[F[1:]];T=A.co(j if sf else X);K=array('d');J=D
			for R in C(H(n,B(I)-1)+1):J+=I[R]*T[n-R];K.append(J)
			A._c[F]=K
		K=A._c[F];return K[H(a,B(K)-1)]
	def cap(A,nb):return floor(ceil(nb)/log2(A.minimum))
	def at(B,nb,sl,sb=G):
		F=nb,sl,sb;C=B._ac.get(F)
		if C is not G:return C
		I=N(nb/B.step);J=(B.cap(nb)+1)/16.;K=B.selection_normalizer*B.two(I,floor((sl-J)/B.step),9**9 if sb is G else floor((sb-J-.125)/B.step),G,O)/B.denominator[I];C=E(D,H(A,K));B._ac[F]=C;return C
	def survivor_upper(A,nb,sl,ex=D,sb=G,lb=G,kmax=4):B=N(nb/A.step);C=floor((sl+ex)/A.step);return(A.two(B,C)if sb is G else A.two(B,C,floor((sb+.25)/A.step),(kmax,lb,(A.cap(nb)+2)/16.+.25)))/A.denominator[B]
def P(a):
	while B(a)>1 and a[-1]==0:a.pop()
	return a
def f(a,b,mod):
	D=B(a);E=B(b);A=[0]*(D+E-1)
	for F in C(D):
		G=a[F]
		if G:
			for H in C(E):
				I=b[H]
				if I:A[F+H]+=G*I
	A=[aa%mod for aa in A]
	while B(A)>1 and A[-1]==0:A.pop()
	return A
def _pmod(a,b,mod):
	D=mod;a=[aa%D for aa in a]
	while B(a)>1 and a[-1]==0:a.pop()
	A=B(b)
	while A>1 and b[A-1]%D==0:A-=1
	G=b[A-1]%D;H=1 if G==1 else J(G,-1,D);K=[b[ab]%D for ab in C(A-1)]
	while B(a)>=A and not(B(a)==1 and a[0]==0):
		I=B(a)-A;E=a.pop()
		if H!=1:E=E*H%D
		if E:
			for F in C(A-1):a[I+F]=(a[I+F]-E*K[F])%D
		while B(a)>1 and a[-1]==0:a.pop()
	return a if a else[0]
def g(a,b,mod):
	A=mod;a=P([aa%A for aa in a]);b=P([ab%A for ab in b]);G=[0]*E(1,B(a)-B(b)+1);I=J(b[-1],-1,A)
	while B(a)>=B(b)and not(B(a)==1 and a[0]==0):
		D=B(a)-B(b);H=a[-1]*I%A;G[D]=H
		for F in C(B(b)):a[D+F]=(a[D+F]-H*b[F])%A
		P(a)
	return P(G),a
def _pgcd(a,b,mod):
	A=mod;a=P([ab%A for ab in a]);b=P([ac%A for ac in b])
	while not(B(b)==1 and b[0]==0):a,b=b,_pmod(a,b,A)
	C=J(a[-1],-1,A);return[aa*C%A for aa in a]
def _ppowmod(base,e,m,mod):
	E=base;A=mod
	if B(m)==7 and B(E)<=6 and m[6]%A==1 and m[5]%A==A-2 and m[4]%A==1 and m[3]%A==2 and m[2]%A==0 and m[1]%A==A-4 and m[0]%A==3:return P(_ppow_h(E,e,A))
	if B(m)==3 and m[2]%A==1:
		I,J=-m[0]%A,-m[1]%A;F,G=1,0;C,D=E[0]%A if E else 0,E[1]%A if B(E)>1 else 0
		while e>0:
			if e&1:M=G*D;F,G=(F*C+M*I)%A,(F*D+G*C+M*J)%A
			M=D*D;C,D=(C*C+M*I)%A,(2*C*D+M*J)%A;e>>=1
		return[F,G]if G else[F]
	if B(m)==4 and m[3]%A==1:
		I,J,O=-m[0]%A,-m[1]%A,-m[2]%A;F,G,N=1,0,0;C,D,H=E[0]%A if E else 0,E[1]%A if B(E)>1 else 0,E[2]%A if B(E)>2 else 0
		while e>0:
			if e&1:K=N*H;L=G*H+N*D+K*O;F,G,N=(F*C+L*I)%A,(F*D+G*C+K*I+L*J)%A,(F*H+G*D+N*C+K*J+L*O)%A
			K=H*H;L=2*D*H+K*O;C,D,H=(C*C+L*I)%A,(2*C*D+K*I+L*J)%A,(2*C*H+D*D+K*J+L*O)%A;e>>=1
		return P([F,G,N])
	R=[1];Q=_pmod(E,m,A)
	while e>0:
		if e&1:R=_pmod(f(R,Q,A),m,A)
		Q=_pmod(f(Q,Q,A),m,A);e>>=1
	return R
def i(f,mod):
	A=mod;f=P([ac%A for ac in f]);M=J(f[-1],-1,A);f=[aa*M%A for aa in f];N=[ab*f[ab]%A for ab in C(1,B(f))]
	if B(_pgcd(f,N,A))!=1:return
	H=[];F=[0,1];D=f[:];G=1
	while B(D)-1>=2*G:
		F=_ppowmod(F,A,D,A);I=F[:]+[0]*E(0,2-B(F));I[1]=(I[1]-1)%A;K=_pgcd(I,D,A);L=B(K)-1
		if L>0:H.extend([G]*(L//G));D,O=g(D,K,A);F=_pmod(F,D,A)
		G+=1
	if B(D)-1>0:H.append(B(D)-1)
	return AG(H)
def _roots_mod(f,mod):
	A=mod;K=random.Random(12345);f=P([ab%A for ab in f]);L=J(f[-1],-1,A);f=[aa*L%A for aa in f];M=B(f)<=7 and _pmod(H_COEFFS,f,A)==[0];H,D=[],[f]
	while D:
		C=D.pop();E=B(C)-1
		if E==0:continue
		if E==1:H.append(-C[0]%A);continue
		while O:
			I=rand_below(K.getrandbits,A)
			if M and 4<=E<=5:F=_pmod(_ppow_h([I,1],(A-1)//2,A),C,A)
			else:F=_ppowmod([I,1],(A-1)//2,C,A)
			F[0]=(F[0]-1)%A;G=_pgcd(F,C,A)
			if 0<B(G)-1<E:N,Q=g(C,G,A);D.append(G);D.append(N);break
	return AG(H)
def _brent(n):
	if n%2==0:return 2
	for F in(1,2,3,5,7,11,13):
		A,J,B,D,G=2,64,1,1,1
		while B==1:
			K=A
			for L in C(D):A=(A*A+F)%n
			I=0
			while I<D and B==1:
				E=A
				for L in C(H(J,D-I)):A=(A*A+F)%n;G=G*Q(K-A)%n
				B=T.gcd(G,n);I+=J
			D<<=1
		if B==n:
			while O:
				E=(E*E+F)%n;B=T.gcd(Q(K-E),n)
				if B>1:break
		if 1<B<n:return B
	return n
def _isp(n):
	if n<2:return K
	for D in(2,3,5,7,11,13,17,19,23,29,31,37):
		if n==D:return O
		if n%D==0:return K
	B,E=n-1,0
	while B%2==0:B>>=1;E+=1
	for F in(2,7,61)if n<4759123141 else(2,325,9375,28178,450775,9780504,1795265022):
		if F%n==0:continue
		A=J(F,B,n)
		if A in(1,n-1):continue
		for G in C(E-1):
			A=A*A%n
			if A==n-1:break
		else:return K
	return O
AH={ac:(1<<ac,1<<ac-1,(1<<ac)-1,ac,2*ac,3*ac,4*ac,5*ac)for ac in(32,44,88)}
def _mh_k(a,b,S=44):
	B,F,G,C,Q,R,T,U=AH[S];A=(a[0]+(a[1]<<C)+(a[2]<<Q)+(a[3]<<R)+(a[4]<<T)+(a[5]<<U))*(b[0]+(b[1]<<C)+(b[2]<<Q)+(b[3]<<R)+(b[4]<<T)+(b[5]<<U));P=A&G
	if P>=F:P-=B;A+=B
	A>>=C;N=A&G
	if N>=F:N-=B;A+=B
	A>>=C;O=A&G
	if O>=F:O-=B;A+=B
	A>>=C;L=A&G
	if L>=F:L-=B;A+=B
	A>>=C;J=A&G
	if J>=F:J-=B;A+=B
	A>>=C;K=A&G
	if K>=F:K-=B;A+=B
	A>>=C;D=A&G
	if D>=F:D-=B;A+=B
	A>>=C;E=A&G
	if E>=F:E-=B;A+=B
	A>>=C;H=A&G
	if H>=F:H-=B;A+=B
	A>>=C;I=A&G
	if I>=F:I-=B;A+=B
	M=A>>C
	if M:I+=2*M;H-=M;E-=2*M;K+=4*M;J-=3*M
	if I:H+=2*I;E-=I;D-=2*I;J+=4*I;L-=3*I
	if H:E+=2*H;D-=H;K-=2*H;L+=4*H;O-=3*H
	if E:D+=2*E;K-=E;J-=2*E;O+=4*E;N-=3*E
	if D:K+=2*D;J-=D;L-=2*D;N+=4*D;P-=3*D
	return P,N,O,L,J,K
def _mh_int(a,b):
	E,F,G,H,I,J=a;K,L,M,N,O,P=b;W=E*K;U=E*L+F*K;V=E*M+F*L+G*K;T=E*N+F*M+G*L+H*K;R=E*O+F*N+G*M+H*L+I*K;S=E*P+F*O+G*N+H*M+I*L+J*K;A=F*P+G*O+H*N+I*M+J*L;B=G*P+H*O+I*N+J*M;C=H*P+I*O+J*N;D=I*P+J*O;Q=J*P
	if Q:D+=2*Q;C-=Q;B-=2*Q;S+=4*Q;R-=3*Q
	if D:C+=2*D;B-=D;A-=2*D;R+=4*D;T-=3*D
	if C:B+=2*C;A-=C;S-=2*C;T+=4*C;V-=3*C
	if B:A+=2*B;S-=B;R-=2*B;V+=4*B;U-=3*B
	if A:S+=2*A;R-=A;T-=2*A;U+=4*A;W-=3*A
	return W,U,V,T,R,S
def _pmul_h(a,b,mod):return[aa%mod for aa in _mh_int(a,b)]
def _ppow_h(a,e,mod):
	A=mod
	if e>0 and B(a)>=2 and a[0]==0 and a[1]==1 and not any(a[2:]):
		D,E,F,G,H,I=1,0,0,0,0,0
		for X in bin(e)[2:]:
			W=D*D;T=2*D*E;U=2*D*F+E*E;R=2*(D*G+E*F);O=2*(D*H+E*G)+F*F;P=2*(D*I+E*H+F*G);J=2*(E*I+F*H)+G*G;K=2*(F*I+G*H);L=2*G*I+H*H;M=2*H*I;N=I*I
			if N:M+=2*N;L-=N;K-=2*N;P+=4*N;O-=3*N
			if M:L+=2*M;K-=M;J-=2*M;O+=4*M;R-=3*M
			if L:K+=2*L;J-=L;P-=2*L;R+=4*L;U-=3*L
			if K:J+=2*K;P-=K;O-=2*K;U+=4*K;T-=3*K
			if J:P+=2*J;O-=J;R-=2*J;T+=4*J;W-=3*J
			D,E,F,G,H,I=W%A,T%A,U%A,R%A,O%A,P%A
			if X=='1':Q=I;D,E,F,G,H,I=-3*Q%A,(D+4*Q)%A,E,(F-2*Q)%A,(G-Q)%A,(H+2*Q)%A
		return[D,E,F,G,H,I]
	V=[1,0,0,0,0,0];S=[a[aa]%A if aa<B(a)else 0 for aa in C(6)]
	while e>0:
		if e&1:V=_pmul_h(V,S,A)
		S=_pmul_h(S,S,A);e>>=1
	return V
def _pinv_h(a,mod):
	A=mod;G,D=[aa%A for aa in H_COEFFS],P([ab%A for ab in a]);F,H=[0],[1]
	while not(B(D)==1 and D[0]==0):K,L=g(G,D,A);I=f(K,H,A);M=E(B(F),B(I));N=P([((F[ac]if ac<B(F)else 0)-(I[ac]if ac<B(I)else 0))%A for ac in C(M)]);G,D,F,H=D,L,H,N
	O=J(G[0],-1,A);return([ad*O%A for ad in F]+[0]*6)[:6]
def _norm_Kh_exact(coeffs_int):
	B=[[F(0)]*6 for aa in C(6)];G=[F(ab)for ab in coeffs_int]
	for A in C(6):
		for D in C(6):B[D][A]=G[D]
		L=G[5];J=[F(0)]+G[:5]
		for E in C(6):J[E]-=L*H_COEFFS[E]
		G=J
	H=F(1)
	for A in C(6):
		I=next(C for C in C(A,6)if B[C][A]!=0)
		if I!=A:B[A],B[I]=B[I],B[A];H=-H
		K=B[A][A];H*=K
		for D in C(A+1,6):
			if B[D][A]!=0:
				M=B[D][A]/K
				for E in C(A,6):B[D][E]-=M*B[A][E]
	return V(H)
def _pmul_h_q(a,b):
	A=[F(0)]*11
	for E in C(6):
		for B in C(6):A[E+B]+=a[E]*b[B]
	for D in C(10,5,-1):
		G=A[D]
		if G:
			A[D]=F(0)
			for B in C(6):A[D-6+B]-=G*H_COEFFS[B]
	return A[:6]
def verify_tower_algebra_and_crt():
	k=J(P_BN254,12,R_BN254)==1 and J(P_BN254,6,R_BN254)!=1 and J(P_BN254,4,R_BN254)!=1;P=[0,1,0,0,0,0];AU=J(3,-1,P_BN254);AV=[aa*AU%P_BN254 for aa in SIGMA_NUM];AW=_ppow_h(P,P_BN254,P_BN254);l=_ppow_h(AW,P_BN254,P_BN254);m=_ppow_h(l,P_BN254,P_BN254);n=_ppow_h(m,P_BN254,P_BN254);AX=_ppow_h(_ppow_h(n,P_BN254,P_BN254),P_BN254,P_BN254);o=AX==P and l!=P and m!=P;AY=n==AV;p=[F(ab,3)for ab in SIGMA_NUM];W=[[F(1)]+[F(0)]*5,p]
	for _ in C(4):W.append(_pmul_h_q(W[-1],p))
	AZ=[I(W_COEFFS[A]*W[A][ac]for A in C(6))for ac in C(6)];q=AZ==[F(ah)for ah in W_COEFFS];r=[-ad%P_BN254 for ad in G_LEAD_COEFFS];Aa=_ppow_h(r,(P_BN254-1)//2,P_BN254);X,R=Aa,[1,0,0,0,0,0];s=t=G
	for u in C(10):
		R=_pmul_h(R,X,P_BN254)
		if u==5:s=R[:]
		if u==9:t=R[:]
		X=_ppow_h(X,P_BN254,P_BN254)
	v=[P_BN254-1,0,0,0,0,0];w=s==v;x=AY and t==v;A=J(P_BN254,10,R_BN254);Ab=36*U_BN254**3+18*U_BN254**2+6*U_BN254+2;y=A==Ab and(A*A-A+1)%R_BN254==0 and(J(A,3,R_BN254)+1)%R_BN254==0;Ac=_pinv_h(r,P_BN254);Y=[0]*6;Z=[1,0,0,0,0,0]
	for b in(0,2,4,6,8):Ad=_pmul_h([ak%P_BN254 for ak in F_REC_COEFFS[b]],Z,P_BN254);Y=[(ai+aj)%P_BN254 for(ai,aj)in zip(Y,Ad)];Z=_pmul_h(Z,Ac,P_BN254)
	z=a(A==0 for A in Y);A0,A1,A2=Q(_norm_Kh_exact(W_COEFFS)),Q(_norm_Kh_exact(F_REC_COEFFS[8])),Q(_norm_Kh_exact(F_REC_COEFFS[0]));A3=i(H_COEFFS,R_BN254);A4=_roots_mod(H_COEFFS,R_BN254);c,d=[],[]
	for A6 in A4:
		Ae=I(W_COEFFS[A]*J(A6,A,R_BN254)for A in C(6))%R_BN254;c.append(i([1,0,((U_BN254+1)*Ae+2)%R_BN254],R_BN254));A8=[0]*9
		for(b,Af)in F_REC_COEFFS.items():A8[b]=I(Af[A]*J(A6,A,R_BN254)for A in C(6))%R_BN254
		d.append(i(A8,R_BN254))
	A9=A3==[1,1,1,1,1,1]and B(A4)==6 and T.gcd(Q(_norm_Kh_exact(G_LEAD_COEFFS)),R_BN254)==1 and a(A==[1,1]for A in c)and a(A==[2,2,2,2]for A in d);e=load_data()[AK];AA=e['SMAT']
	def AB(G):return E(Q(I(AA[D][A]*G[D][E]*AA[E][B]for D in C(6)for E in C(6))-G[A][B])for A in C(6)for B in C(6))
	AC,AD=AB(e[A5]),AB(e['GB']);AE=T.gcd(R_BN254,J(P_BN254,6,R_BN254)-1)==1;AF=[(17*ae+3)%P_BN254 for ae in C(6)];AG=[(29*af+11)%P_BN254 for af in C(6)];f=[(7*ag+5)%P_BN254 for ag in C(6)];AH=_pinv_h(AF,P_BN254);Ag=_pmul_h(AG,AH,P_BN254);Ah=_pmul_h(Ag,f,P_BN254);AI=_pmul_h(f,AH,P_BN254);AJ=_pmul_h(AI,AF,P_BN254)==f and _pmul_h(AI,AG,P_BN254)==Ah
	def Ai(n):
		B,D=n-1,0
		while B%2==0:B//=2;D+=1
		for E in(2,325,9375,28178,450775,9780504,1795265022):
			if E%n==0:continue
			A=J(E,B,n)
			if A in(1,n-1):continue
			for F in C(D-1):
				A=A*A%n
				if A==n-1:break
			else:return K
		return O
	S=1
	for Aj in NTT3_PRIMES:S*=Aj
	def Ak(p):
		B=(1<<64)-1;N=J(1<<64,2,p);I=-J(p,-1,1<<64)&B;D=random.Random(p&4294967295);L=[(p-1,p-1),(p-2,p-2),(2*p-1,2*p-1),(p-1,2*p-1),(0,p-1),(1,p-2)]+[(D.randrange(2*p),D.randrange(2*p))for al in C(32)]
		for(E,F)in L:
			G=E*F;M=(G&B)*I&B;H=G+M*p
			if H>=1<<128:return K
			A=H>>64
			if A>=2*p:return K
			if A>=p:A-=p
			if A!=E*F*J(1<<64,-1,p)%p:return K
		return O
	AL=a(Ak(A)for A in NTT3_PRIMES)
	def Al(p):
		I=(1<<64)-1;F=random.Random((p^21930)&4294967295)
		for M in C(64):
			B=F.randrange(p);D=F.randrange(p);E=F.randrange(p);J=(E<<64)//p;L=D*J>>64;A=D*E-L*p&I
			if A>=2*p:return K
			if A>=p:A-=p
			G=B+A
			if G>=p:G-=p
			H=B-A
			if H<0:H+=p
			if G!=(B+D*E)%p or H!=(B-D*E)%p:return K
		return O
	AM=a(Al(A)for A in NTT3_PRIMES);M,D,L=NTT3_PRIMES;AN=J(M,-1,D);AO=J(M*D%L,-1,L);g=random.Random(20261120);h=O
	for _ in C(32):
		V=g.randrange(S);AP,Am,An=V%M,V%D,V%L;Ao=(Am-AP)%D*AN%D;AQ=AP+Ao*M;Ap=(An-AQ%L)%L*AO%L;Aq=AQ+Ap*M*D
		if Aq!=V:h=K
	def j(arr,w_root,p,inv=K):
		r=w_root;M=B(arr);A=2*p;A3=(1<<64)-1;i=J(r,-1,p)if inv else r
		def N(v,w,wp):
			C=v*wp>>64;B=v*w-C*p&A3
			if B>=A:return-1
			return B
		T=M;a=[]
		if T%3==0:a.append(3);T//=3
		while T%4==0:a.append(4);T//=4
		if T==2:a.append(2);T//=2
		F=[0]*M
		for s in C(M):
			t=s;j=0
			for S in A7(a):j=j*S+t%S;t//=S
			F[j]=arr[s]%p
		u=J(i,M//4,p);A4=(u<<64)//p;b=J(i,M//3,p)if M%3==0 else 1;v=(b<<64)//p;k=b*b%p;w=(k<<64)//p;U=1
		for S in a:
			I=U;U*=S;l=J(i,M//U,p);D=1;V=[]
			if S==4:
				for x in C(I):P=D*D%p;h=P*D%p;V.append((D,(D<<64)//p,P,(P<<64)//p,h,(h<<64)//p));D=D*l%p
				for G in C(0,M,U):
					for E in C(I):
						D,c,P,m,h,A5=V[E];W=F[G+E];Q=N(F[G+E+I],D,c);R=N(F[G+E+2*I],P,m);n=N(F[G+E+3*I],h,A5)
						if Q<0 or R<0 or n<0:return
						X=W+R;X=X-A if X>=A else X;Y=W+A-R;Y=Y-A if Y>=A else Y;Z=Q+n;Z=Z-A if Z>=A else Z;o=N(Q+A-n,u,A4)
						if o<0:return
						K=X+Z;K=K-A if K>=A else K;O=X+A-Z;O=O-A if O>=A else O;L=Y+o;L=L-A if L>=A else L;d=Y+A-o;d=d-A if d>=A else d;F[G+E],F[G+E+I],F[G+E+2*I],F[G+E+3*I]=K,L,O,d
			elif S==2:
				for x in C(I):V.append((D,(D<<64)//p));D=D*l%p
				for G in C(0,M,U):
					for E in C(I):
						D,c=V[E];y=F[G+E];q=N(F[G+E+I],D,c)
						if q<0:return
						K=y+q;K=K-A if K>=A else K;L=y+A-q;L=L-A if L>=A else L;F[G+E],F[G+E+I]=K,L
			else:
				for x in C(I):P=D*D%p;V.append((D,(D<<64)//p,P,(P<<64)//p));D=D*l%p
				for G in C(0,M,U):
					for E in C(I):
						D,c,P,m=V[E];W=F[G+E];Q=N(F[G+E+I],D,c);R=N(F[G+E+2*I],P,m)
						if Q<0 or R<0:return
						e=Q+R;e=e-A if e>=A else e;K=W+e;K=K-A if K>=A else K;z=N(Q,b,v);A0=N(R,k,w);A1=N(Q,k,w);A2=N(R,b,v)
						if H(z,A0,A1,A2)<0:return
						f=z+A0;f=f-A if f>=A else f;L=W+f;L=L-A if L>=A else L;g=A1+A2;g=g-A if g>=A else g;O=W+g;O=O-A if O>=A else O;F[G+E],F[G+E+I],F[G+E+2*I]=K,L,O
		if inv:A6=J(M,-1,p);return[an*A6%p for an in F]
		return[am%p for am in F]
	def AR(n_len):
		A=n_len;N=[g.randrange(1<<64)for ao in C(A//2)]+[0]*(A-A//2);P=[g.randrange(1<<64)for ap in C(A//2)]+[0]*(A-A//2);Q=[0]*A
		for E in C(A//2):
			for R in C(A//2):Q[E+R]+=N[E]*P[R]
		F=[]
		for B in NTT3_PRIMES:
			H=2
			while J(H,(B-1)//2,B)==1 or J(H,(B-1)//3,B)==1:H+=1
			I=J(H,(B-1)//A,B);S=j(N,I,B,K);T=j(P,I,B,K)
			if S is G or T is G:return K
			X=[S[aq]*T[aq]%B for aq in C(A)];U=j(X,I,B,O)
			if U is G:return K
			F.append(U)
		for E in C(A):
			V,Y,Z=F[0][E],F[1][E],F[2][E];a=(Y-V)%D*AN%D;W=V+a*M;b=(Z-W%L)%L*AO%L
			if W+b*M*D!=Q[E]:return K
		return O
	AS=AR(16)and AR(24);AT=B(NTT3_PRIMES)==3 and a(A<1<<62 and(A-1)%(3*(1<<55))==0 and Ai(A)for A in NTT3_PRIMES)and S>(1<<54)*((1<<64)-1)**2 and AL and AM and h and AS;Ar=k and o and w and x and q and y and z and A9 and AE and AJ and AT and A0==216 and A1==1 and A2==729 and AC<1e-10 and AD<1e-10;return U(ok=Ar,ord_r_p_12=k,h_irred_mod_p=o,g_rec_irred_mod_p6=w,tau_is_frob10=x,sigma_fixes_w_exact=q,chi_192bit_ok=y,f_rec_vanishes_at_inv_uw=z,gcd_r_p6_minus_1_is_1=AE,boot_identity_ok=AJ,ntt3_primes_verified=AT,ntt3_montgomery_no_carry_ok=AL,ntt3_shoup_butterfly_ok=AM,ntt3_garner_crt_ok=h,ntt3_r4_r3_shoup_conv_ok=AS,ntt3_prod_bits=N(log2(S),3),norm_w=A0,norm_f_lead=A1,norm_f_const=A2,mod_r_h_degs=A3,mod_r_g_degs=c,mod_r_f_degs=d,mod_r_crt_all_deg2_ok=A9,sigma_iso_err_Ga=AC,sigma_iso_err_Gb=AD)
def build_exact_int_geom(T4=G):
	D=T4
	if D is G:D=load_data()[AK]
	P,L,B=D[A5],D['UMATS'],D['SMAT']
	def A(A,B):return[[I(A[ah][C]*B[C][ai]for C in C(6))for ai in C(6)]for ah in C(6)]
	def H(A):return[[A[ak][aj]for ak in C(6)]for aj in C(6)]
	M=1<<16;E=[[N(P[aa][ad]*M)for ad in C(6)]for aa in C(6)]
	for J in C(6):
		for O in C(J+1,6):E[O][J]=E[J][O]
	K=A(B,B);Q=A(H(B),A(E,B));R=A(H(K),A(E,K));F=[[E[ab][ae]+Q[ab][ae]+R[ab][ae]for ae in C(6)]for ab in C(6)];assert A(H(B),A(F,B))==F;S=[[[A(H(ac),A(F,ac))[af][ag]-F[af][ag]for ag in C(6)]for af in C(6)]for ac in L];return F,S,L,B,K,3*M
def verify_exhaustive_fd(R2=5e1):
	A6=load_data()[AK];A7,Be,Ad,Bf,Bg,A8=build_exact_int_geom(A6)
	def Ae(M):return p(p((A,M[B][A])for A in C(6)if M[B][A])for B in C(6))
	Af,B6=Ae(Bf),Ae(Bg);B7=p(Ae(A)for A in Ad)
	def A9(sp,v):return p(I(B*v[A]for(A,B)in A)for A in sp)
	def B8(Mat,v):
		D=0
		for A in C(6):
			B=v[A]
			if B:
				D+=Mat[A][A]*B*B
				for E in C(A+1,6):
					if v[E]:D+=2*Mat[A][E]*B*v[E]
		return D
	Ag=V(N(R2*A8));from fractions import Fraction as A0;from math import isqrt;AA=[A0(0)]*6;A1=[[A0(1 if aa==ai else 0)for ai in C(6)]for aa in C(6)]
	for P in C(6):
		g=A0(A7[P][P])-I(A1[P][A]*A1[P][A]*AA[A]for A in C(P));assert g>0;AA[P]=g
		for q in C(P+1,6):A1[q][P]=(A0(A7[q][P])-I(A1[q][A]*A1[P][A]*AA[A]for A in C(P)))/g
	AB,AL=[],[0]*6
	def B9(idx,rem):
		A=idx;B=I(A1[B][A]*AL[B]for B in C(A+1,6));D=rem/AA[A];E=A0(isqrt(D.numerator*D.denominator)+1,D.denominator);J=AA[A]
		for F in C((-B-E).__floor__(),(-B+E).__ceil__()+1):
			G=F+B;H=rem-J*G*G
			if H>=0:
				AL[A]=F
				if A==0:
					if any(AL):AB.append(p(AL))
				else:B9(A-1,H)
	B9(5,A0(Ag))
	def Bh(p):
		if(-p[0],-p[1],-p[2],-p[3],-p[4],-p[5])<p:return K
		for E in C(6):
			F=B8(Be[E],p)
			if F<0:return K
			if F==0:
				D=A9(B7[E],p)
				for A in(D,A9(Af,D),A9(B6,D)):
					G=-A[0],-A[1],-A[2],-A[3],-A[4],-A[5]
					if(G if G<A else A)<p:return K
		for I in(Af,B6):
			B=A9(I,p);H=-B[0],-B[1],-B[2],-B[3],-B[4],-B[5]
			if(H if H<B else B)<p:return K
		return O
	Bi=set(AB);Ah=set();AC=[]
	for AM in AB:
		if AM in Ah:continue
		Ai=[AM];Aj={AM};Ah.add(AM)
		while Ai:
			Bj=Ai.pop()
			for Bk in B7+(Af,):
				X=A9(Bk,Bj)
				for AD in(X,(-X[0],-X[1],-X[2],-X[3],-X[4],-X[5])):
					if AD in Bi and AD not in Aj:Aj.add(AD);Ah.add(AD);Ai.append(AD)
		AC.append(Aj)
	Ak=Al=Am=0
	for Bl in AC:
		P=I(1 for A in Bl if Bh(A))
		if P==1:Ak+=1
		elif P==0:Al+=1
		else:Am+=P-1
	BA=E(Q(A7[A][B]/M(A8)-A6[A5][A][B])for A in C(6)for B in C(6))
	def BB(a,b):
		B=[0]*11
		for D in C(6):
			for A in C(6):B[D+A]+=a[D]*b[A]
		for E in C(10,5,-1):
			F=B[E]
			for A in C(6):B[E-6+A]-=F*H_COEFFS[A]
		return B[:6]
	def Bm(p):
		B,A=[],list(p)
		for G in C(6):
			B.append(A)
			if G<5:
				D=[0]+A[:5];E=A[5]
				if E:
					for F in C(6):D[F]-=E*H_COEFFS[F]
				A=D
		return I(B[A][A]for A in C(6))
	if _OHM_CACHE is G:is_ideal_primitive_oh([1,0,0,0,0,0],[0,1,0,0,0,0])
	Bn=_OHM_CACHE;Bo=a(BB(ZK_NUM[A],ZK_NUM[B])==[3*I(Bn[A][B][C]*ZK_NUM[C][aD]for C in C(6))for aD in C(6)]for A in C(6)for B in C(6));Bp=[[1,0,1,-1,-3,-3],[0,1,-1,0,-1,-1],[0,0,-1,1,1,3],[0,0,-1,-1,0,2],[0,0,1,0,1,-1],[0,0,0,1,1,2]];Bq=[[I(Bp[A][aj]*ZK_NUM[A][an]for A in C(6))for an in C(6)]for aj in C(6)]==[[3 if ak==ao else 0 for ao in C(6)]for ak in C(6)];BC=AJ([[Bm(BB(ZK_NUM[al],ZK_NUM[ap]))//9 for ap in C(6)]for al in C(6)])
	def r(a,p):
		A=[ay%p for ay in a]
		while B(A)>1 and A[-1]==0:A.pop()
		return A
	def AN(a,b):
		A=[0]*(B(a)+B(b)-1)
		for(C,D)in y(a):
			for(E,F)in y(b):A[C+E]+=D*F
		return A
	def An(a,b,p):
		a,b=r(a,p),r(b,p)
		while b!=[0]:
			A=a[:];D=J(b[-1],-1,p)
			while B(A)>=B(b)and A!=[0]:
				G=A[-1]*D%p;F=B(A)-B(b)
				for E in C(B(b)):A[F+E]=(A[F+E]-G*b[E])%p
				A=r(A,p)
			a,b=b,A
		D=J(a[-1],-1,p);return[az*D%p for az in a]
	Ao=[1,0,1,1];BD=AN(Ao,Ao);Br=r(BD,2)==r(H_COEFFS,2)and An(Ao,[(H_COEFFS[ar]-BD[ar])//2 for ar in C(7)],2)==[1];AO=[14,6,18,1];AP=[6,1];BE=AN(AN(AP,AP),AN(AP,AO));Bs=r(BE,19)==r(H_COEFFS,19)and An(AO,[at*AO[at]for at in C(1,4)],19)==[1]and I(AO[A]*J(13,A,19)for A in C(4))%19!=0 and An(AP,[(H_COEFFS[au]-BE[au])//19 for au in C(7)],19)==[1];Ap=Bo and Bq and BC==-2**9*19**2 and Br and Bs;from decimal import Decimal as D,Context as BF,ROUND_FLOOR as Bt,ROUND_CEILING as BG;j,h=BF(prec=55,rounding=Bt),BF(prec=55,rounding=BG);k,c,Y,Bu,R,W,AQ,BH=j.add,h.add,j.subtract,h.subtract,j.multiply,h.multiply,j.divide,h.divide;A=lambda x:(D(x),D(x));S=lambda a,b:(k(a[0],b[0]),c(a[1],b[1]));Z=lambda a,b:(Y(a[0],b[1]),Bu(a[1],b[0]));L=lambda a,b:(H(R(A,B)for A in a for B in b),E(W(A,B)for A in a for B in b));d=lambda a:(D(0),E(W(a[0],a[0]),W(a[1],a[1])))if a[0]<=0<=a[1]else(R(a[0],a[0]),W(a[1],a[1]))if a[0]>0 else(R(a[1].copy_negate(),a[1].copy_negate()),W(a[0].copy_negate(),a[0].copy_negate()))
	def i(a,b):assert b[0]>0 or b[1]<0;return H(AQ(A,B)for A in a for B in b),E(BH(A,B)for A in a for B in b)
	def AR(a):
		A,B=j.sqrt(a[0]),h.sqrt(a[1])
		while W(A,A)>a[0]:A=j.next_minus(A)
		while R(B,B)<a[1]:B=h.next_plus(B)
		return A,B
	AE=D('1e-40');Bv,Bw=Y(D(1),AE),c(D(1),AE);BI=lambda a:(Y(j.ln(a[0]),AE),c(h.ln(a[1]),AE));BJ=lambda a:(R(j.exp(a[0]),Bv),W(h.exp(a[1]),Bw));BK=lambda a,b:(S(a[0],b[0]),S(a[1],b[1]));Aq=lambda a,b:(Z(L(a[0],b[0]),L(a[1],b[1])),S(L(a[0],b[1]),L(a[1],b[0])));BL=[D(ab)for ab in H_COEFFS];BM=[D(ac*H_COEFFS[ac])for ac in C(1,7)];l=D('1e-40');s=[];Ar=O
	for T in C(3):
		e,f=D(str(A6[z][T][1][0])),D(str(A6[z][T][1][1]))
		for _ in C(7):
			AF,m=BL[6],D(0);t,u=BM[5],D(0)
			for F in C(5,-1,-1):AF,m=k(Y(R(AF,e),R(m,f)),BL[F]),k(R(AF,f),R(m,e))
			for F in C(4,-1,-1):t,u=k(Y(R(t,e),R(u,f)),BM[F]),k(R(t,f),R(u,e))
			BN=k(R(t,t),R(u,u));e=Y(e,AQ(k(R(AF,t),R(m,u)),BN));f=Y(f,AQ(Y(R(m,t),R(AF,u)),BN))
		AS=(e,e),(f,f);b=[(A(0),A(0))]*7;b[0]=A(H_COEFFS[6]),A(0)
		for A2 in C(5,-1,-1):
			for F in C(6,0,-1):b[F]=BK(Aq(b[F],AS),b[F-1])
			b[0]=BK(Aq(b[0],AS),(A(H_COEFFS[A2]),A(0)))
		As=AR(S(d(b[0][0]),d(b[0][1])))[1];At=W(l,l)
		for F in C(2,7):As=c(As,W(AR(S(d(b[F][0]),d(b[F][1])))[1],At));At=W(At,l)
		if not As<R(AR(S(d(b[1][0]),d(b[1][1])))[0],l):Ar=K
		s.append(((Y(e,l),c(e,l)),(Y(f,l),c(f,l))))
	Bx=a(A[1][0]>D('0.05')for A in s)and a(E(Y(s[A][0][0],s[B][0][0]).copy_abs(),Y(s[A][1][0],s[B][1][0]).copy_abs())>D('0.05')for A in C(3)for B in C(A+1,3));v=[]
	for AS in s:
		AT=[(A(1),A(0))]
		for _ in C(5):AT.append(Aq(AT[-1],AS))
		n=[]
		for F in C(6):
			w,x=A(0),A(0)
			for A2 in C(6):
				if ZK_NUM[F][A2]:AG=A(ZK_NUM[F][A2]);w=S(w,L(AG,AT[A2][0]));x=S(x,L(AG,AT[A2][1]))
			n.append((i(w,A(3)),i(x,A(3))))
		v.append(n)
	BO=i(A(str(A6[A5][0][0])),A(3));AU=[]
	for o in C(6):
		n=[]
		for T in C(6):
			Au=A(0)
			for F in C(3):Au=S(Au,S(L(v[F][o][0],v[F][T][0]),L(v[F][o][1],v[F][T][1])))
			n.append(L(BO,Au))
		AU.append(n)
	BP=A('0.88');AV=[[A(1 if ad==am else 0)for am in C(6)]for ad in C(6)];Av=[A(0)]*6;BQ=O
	for P in C(6):
		g=Z(AU[P][P],BP)
		for F in C(P):g=Z(g,L(d(AV[P][F]),Av[F]))
		if g[0]<=0:BQ=K;break
		Av[P]=g
		for q in C(P+1,6):
			Aw=AU[q][P]
			for F in C(P):Aw=Z(Aw,L(L(AV[q][F],AV[P][F]),Av[F]))
			AV[q][P]=i(Aw,g)
	Ax=A(0)
	for o in C(6):
		for T in C(6):Ax=S(Ax,d(Z(i(A(A7[o][T]),A(A8)),AU[o][T])))
	AW=i(AR(Ax),BP)[1];AH=A('1e-4');A3=[]
	for m in C(6):
		n=[]
		for T in C(3):
			w,x=A(0),A(0)
			for F in C(6):
				if Ad[m][F][0]:AG=A(Ad[m][F][0]);w=S(w,L(AG,v[T][F][0]));x=S(x,L(AG,v[T][F][1]))
			n.append(S(d(w),d(x)))
		A3.append(n)
	def BR(m):return S(Z(L(m[0][0],Z(L(m[1][1],m[2][2]),L(m[1][2],m[2][1]))),L(m[0][1],Z(L(m[1][0],m[2][2]),L(m[1][2],m[2][0])))),L(m[0][2],Z(L(m[1][0],m[2][1]),L(m[1][1],m[2][0]))))
	AX=[]
	for o in C(6):
		for T in C(o+1,6):
			BS=[[Z(A3[o][av],A(1))for av in C(3)],[Z(A3[T][aw],A(1))for aw in C(3)],[A(1)]*3];Ay=BR(BS)
			if Ay[0]<=D('1e-20')and Ay[1]>=D('-1e-20'):continue
			X=[]
			for Az in C(3):AY=[ax[:]for ax in BS];AY[0][Az]=AY[1][Az]=AH[1].copy_negate(),AH[0].copy_negate();AY[2][Az]=A(1);X.append(i(BR(AY),Ay))
			if a(A[0]>0 for A in X):
				BT=O
				for By in C(6):
					A_=A(0)
					for F in C(3):A_=S(A_,L(X[F],Z(A3[By][F],A(1))))
					if A_[1]<Y(AH[0].copy_negate(),D('1e-20')):BT=K;break
				if BT:AX.append(X)
	B0=H(H(A[0]for A in A)for A in AX);B1=BH(c(D(1),AH[1]),B0);Bz=c(h.ln(B1),AE);AZ=[BI(ae)for ae in A3[0]];A4=[BI(af)for af in A3[2]];BU=Z(L(AZ[0],A4[1]),L(AZ[1],A4[0]));B_=E(A.copy_abs()for A in i(A4[1],BU));C0=E(A.copy_abs()for A in i(A4[0],BU));AI=V(ceil(M(W(W(c(B_,C0),D(2)),Bz))))+2;C1={(0,0),(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1)};B2=0;Aa=D('1e9')
	for BV in C(-AI,AI+1):
		for BW in C(-AI,AI+1):
			if(BV,BW)in C1:continue
			BX=[BJ(S(L(A(BV),AZ[aq]),L(A(BW),A4[aq])))for aq in C(3)]
			if E(A[0]for A in BX)>=B1:continue
			B2+=1
			for X in AX:
				Ab=A(0)
				for T in C(3):Ab=S(Ab,L(X[T],Z(BX[T],A(1))))
				if Ab[0]<Aa:Aa=Ab[0]
	B3=D(0)
	for BY in(D('-0.5'),D('0.5')):
		for BZ in(D('-0.5'),D('0.5')):
			Ac=A(0)
			for T in C(3):Ac=S(Ac,BJ(S(L((BY,BY),AZ[T]),L((BZ,BZ),A4[T]))))
			if Ac[1]>B3:B3=Ac[1]
	Ba=W(W(BO[1],B3),D('1.01'));B4=Ar and Bx and BQ and B(AX)==6 and AW<D('3.0e-5')and R(AQ(Y(D(1),AW),c(D(1),AW)),k(D(1),AH[0]))>D('1.00004')and B0>D('0.09')and B2==3 and Aa>D('0.40');Bb=V(W(Ba,D(A8)).to_integral_value(rounding=BG))
	def C2(vec):
		H=[I(vec[A]*ZK_NUM[A][aA]for A in C(6))for aA in C(6)];B,A=[],H[:]
		for J in C(6):
			B.append(A)
			if J<5:
				D=[0]+A[:5];E=A[5]
				if E:
					for F in C(6):D[F]-=E*H_COEFFS[F]
				A=D
		G=Q(AJ([[B[aC][aB]for aC in C(6)]for aB in C(6)]));assert G%729==0;return G//729
	Bc=[ag for ag in AB if B8(A7,ag)<=Bb];B5=[ah for ah in Bc if C2(ah)==1];Bd=Ap and Bb<Ag and set(B5)=={(1,0,0,0,0,0),(-1,0,0,0,0,0)};return U(ok=Ak==B(AC)and Al==0 and Am==0 and B(AC)==382 and BA<1e-05 and Ap and B4 and Bd,zero_float_geom=O,sigma_inv_exact=O,max_order_dedekind_ok=Ap,disc_Oh=BC,exact_frac_enum_ok=O,certified_interval_ok=B4,rouche_roots_ok=Ar,G_sym_scale=A8,G_sym_rel_err=M(f"{BA:.2e}"),unit_index1_proved=Bd,unit_box_R2_hi=N(M(Ba),4),unit_box_pts=B(Bc),unit_box_units=B(B5),roots_of_unity_exact=B(B5),global_hexagon_proof_ok=B4,eps_G_hi=N(M(AW),8),eps_inf=.0001,V_min_lo=N(M(B0),6),B_max_hi=N(M(B1),4),M_box_bound=AI,inner_nonfacet_units=B2,min_nonfacet_gap_lo=N(M(Aa),6),R2=R2,T_lim_int=Ag,nonzero_pts=B(AB),total_orbits=B(AC),exact1=Ak,missed_orbits=Al,duplicate_orbits=Am)
def log2_li(log2x):C=log2x;B=C*L;return C-log2(B)+log2(A+A/B+2./(B*B))
def log2_U_z6(lf,lg):A=2.**log2_li(lf)/6.;B=2.**log2_li(lg)/6.;C=(2.**(.5*lf)/(.5*lf*L)+2.**(.5*lg)/(.5*lg*L))/6.;return log2(A+B+C)
NTT3_PRIMES=69*(1<<55)+1,57*(1<<55)+1,54*(1<<55)+1
NTT2_52_PRIMES=4194240*(1<<40)+1,4194180*(1<<40)+1
def A2(m,t=A,d=A):
	if m<=1:return 22e1*d
	if m>1<<55:return d*ntt3_F(m)if d>1 else 3.*ntt3_F(ceil(m/2.))+2e1*m
	A=E(2,ceil(log2(m)));C=M(1<<A);F=t*46.5*C*A+d*(22e1*C+(u*C if C>r else D));B=3.*M(1<<A-2);return H(F,t*46.5*B*(A-2)+d*(364.*B+(u*B if B>r else D)))if B>=m and A>=3 else F
def ntt3_F(M_prod):return A2(ceil(E(M_prod,A)))
def fft_F(M_prod):return ntt3_F(M_prod)
def _kar_sb(n):return 6.*n*n+16.*n+2e1 if n<=16 else H(6.*n*n+16.*n+2e1,3.*_kar_sb((n+1)//2)+28.*n)
def A3(M_prod,two=K):
	F=ceil(E(M_prod,A));K=E(1,ceil(F/2.));G=5./3. if two else A;D=2. if two else A;C=A2(F,G,D)
	if K<=512:C=H(C,D*_kar_sb(K))
	I=ceil(F*64./52.)
	if I<=1<<18:
		B=E(2,ceil(log2(I)));J=16.*D*F;L=M(1<<B);C=H(C,2./3.*(G*46.5*L*B+D*22e1*L)+J)
		if B>=3 and 3*(1<<B-2)>=I:N=3.*M(1<<B-2);C=H(C,2./3.*(G*46.5*N*(B-2)+D*364.*N)+J)
		if B>=4 and 5*(1<<B-3)>=I:O=5.*M(1<<B-3);C=H(C,2./3.*(G*46.5*O*(B-3)+D*42e1*O)+J)
	return C
def mul_F(M_prod):return A3(M_prod,K)
def mul2_F(M_prod):return A3(M_prod,O)
KH_SPLIT_CAP=.3
KH_F_CAP=.08
KH_G_CAP=.17
KH_FG_CAP=.2
def A4(lo,hi):return 1.000028*(2.**hi-(2.**lo if lo else 0))/L if lo==0 else 1.000028*2.**hi/L
def pbits(lo,hi,cap=G):A=lo;C=KH_SPLIT_CAP if cap is G else cap;B=.5*hi;F=A4(A,B)if not A or A<B else D;return F+C*A4(E(A,B)if A else B,hi)
A6={}
def AF(pb,leaf=9,L=G,c_leaf=AE):
	F=c_leaf;D=leaf;J=ceil(pb/64.)
	if L is G:M=ceil(log2((J+1)/D));return H(AF(pb,D,A,F)for A in C(E(1,M-3),M+2))
	N=J,D,L,F;A=A6.get(N)
	if A is not G:return A
	K=L+4;B=lambda j:D*2**j;O=I(2**(L-A)*(mul_F(2*B(A-1))+s*2*B(A-1))for A in C(1,L+1));P=I(2**(L-A)*(mul2_F(B(A)+K)+2.*s*(B(A-1)+K))for A in C(1,L+1));A=(O+P+AD*mul_F(2*(E(J,B(L))+K)))/2**L+F;A6[N]=A;return A
def _det_mat(M):
	G=B(M);F=[aa[:]for aa in M];H=A
	for D in C(G):
		I=E(C(D,G),key=lambda r:Q(F[r][D]))
		if I!=D:F[D],F[I]=F[I],F[D];H=-H
		J=F[D][D];H*=J
		for K in C(D+1,G):
			N=F[K][D]/J
			for L in C(D,G):F[K][L]-=N*F[D][L]
	return H
def DZ(lf,lg,exc,nf=403.5,ng=468.2):
	O=log2_U_z6(lf,lg);H=[]
	for(P,B)in((nf,lf),(ng,lg)):
		R=P/B;S=V(N((B-3e1)*10))
		for T in C(S):I=3e1+T*.1;E=I+.1;J=kappa_rho(R,E/B);H.append((2.**(O-E)*6.*J,J*log(E/I)))
	def K(e_val):
		def F(h):
			B=D;E=A-h
			for(F,G)in H:
				C=F*e_val
				if C<1e2:B+=G*exp(-C*E)
			return A-exp(-B)
		B=A
		for G in C(100):
			E=F(B)
			if Q(E-B)<1e-10:return E
			B=E
		return B
	L=K(exc);M,F=.5,A
	for W in C(24):
		G=.5*(M+F)
		if K(G)>.75:M=G
		else:F=G
	return U(h=L,exc_crit=F,stalls=L>.75)
Am={}
_DESC_CACHE={}
def Da(W_k,Rb,g00):
	K=W_k*sqrt(g00)/Rb;M=N=D;L=120
	for S in C(L):
		B=(S+.5)/L;O=5.*B**4/L;F=sqrt(E(D,A-B*B));N+=O*(32./35.)*F**7;P=D
		for T in C(4):
			G=-F-T/4.*K
			while G<F:
				I=E(-F,G);J=H(F,G+K)
				if J>I:R=D if I<=D<=J else H(Q(I),Q(J));P+=(J-I)*E(D,A-B*B-R*R)**3
				G+=K
		M+=O*(P/4.)
	return M/N
def compute_ledger(track='bn254-gt-dlp',use_db=O,c_tail_z6=G,cap_scale=A,strict_per_owner=O,swar_probe=G,srt_g_cap=G):
	As=srt_g_cap;Ar=strict_per_owner;Aq=use_db;AH=c_tail_z6;A2=swar_probe;m=cap_scale;AI=load_data();J=AI[Bm];i=AI[AK];Db=AI['filtsim_raw'];j,e,At=OP['lf'],OP['lg'],OP['exc'];n,f,o,Au=OP['p_min'],OP['p_cut'],OP['p_lat'],OP['lx2'];assert o==OP['p_mid']and Au==e;Av,FK,Aw=OP['slack'],OP['slack_b'],OP['w_sb']
	if AH is G:AH=OP['c_tail_z6']
	if A2 is G:A2=D
	import v32_verify as AJ;Bt=AJ.get_ret2d_u8(sys.modules[__name__],pmin=n,B=e,odd=2);Dc,Bu=Bt['get_rg'],Bt['get_sf_s0'];p=AJ.get_exact_23_classes();FL=p['classes'];AL=ALPHA_F_RAW+p['ex_f_bits'];Dd=ALPHA_G_RAW+p['ex_g_bits'];De=H(A,D7/(J[A1]/M(J['n_in'])));A3=log2_U_z6(j,e);Ax=A3+log2(At)+OP['margin'];A4=J['frac_reg'];Df=J['frac_reg_trunc'];Dg=A/A4;A6=A/Df;Bv=J['log2_V0']+log2(A4);W=J['n_in'];Ay=J[Ao];assert I(A[2]for A in Ay)==J[A1]and I(A[1]for A in J[Bn])==J[Ap]and J[A1]+J[Ap]==W;k=rho_tab_1024();Bw=B(k)-1;Az=A/j;Bx=A/e;Dh=D6/L;A_=[(ab,ac,aa)for(ab,ac,aa,ad)in Ay if aa];By=AG({ai[0]for ai in A_});Bz=AG({aj[1]for aj in A_});B0,B_,C0=p['u_vf'],p['u_vg'],p['sparse_fg'];B1=B(B0);B2=T.sumprod
	def C1(u,ln):
		if u<=A:return A
		E=u*1024.;B=V(E);I=k[B]+(E-B)*(k[B+1]-k[B])if B<Bw else rho(u);D=u-A
		if D<=A:F=A
		else:G=D*1024.;C=V(G);F=k[C]+(G-C)*(k[C+1]-k[C])if C<Bw else rho(D)
		H=I-Dh*F/ln;return H if H>1e-300 else 1e-300
	def AM(d):
		C=N=O=P=Q=D;F=48.*d;I=12.*d;i=.9986*De;j={am:[C1((am+F-an+AL)*Az,am+F-an+AL)if Aq else rho((am+F-an+AL)*Az)for an in B0]for am in By};R={};S={}
		for A in Bz:
			T=A+I;J=[];K=[]
			for L in B_:
				U=Av-L
				if U<1e1:J.append(D);K.append(D);continue
				k=H(484,E(10,ceil(T-L)));V=Dc(k,U);M=T-L+Dd;X=(C1(M*Bx,M)if Aq else rho(M*Bx))*V;J.append(X);K.append(X*V)
			Y=[D]*B1;Z=[D]*B1
			for(a,b,c)in C0:Y[b]+=a*J[c];Z[b]+=a*K[c]
			R[A]=Y;S[A]=Z
		for(e,A,f)in A_:g=j[e];B=B2(g,R[A]);l=B2(g,S[A]);m=l/B if B>D else D;B*=i;G=f*B;C+=G;N+=f*(B*B);O+=G*(e+F);P+=G*(A+I);Q+=G*m
		h=C/W;return Bv+I+log2(E(h,1e-300)),h,N/W,O/C-CONTENT_BITS,P/C,Q/C
	AN,A7=-.1,.1;B3,AO=AM(AN)[0]-Ax,AM(A7)[0]-Ax
	for _ in C(9):
		AP=A7-AO*(A7-AN)/(AO-B3)if AO!=B3 else .5*(AN+A7);AQ=AM(AP)[0]-Ax
		if AQ<D:AN,B3=AP,AQ
		else:A7,AO=AP,AQ
		if Q(AQ)<1e-10:break
	q=AP;C2,B4,Di,B5,B6,Dj=AM(q);C3=sqrt(E(Di-B4*B4,D)/W)/E(B4,1e-300)/L;C4=sqrt((A-A4)/(J[Bo]*A4))/L;B7=sqrt(C3*C3+C4*C4);C5=Bv+12.*q;AR=2.**C5;B8=i[A5];r=i['GB'];B9=r[0][0];Dk=_det_mat(B8);Dl=_det_mat(r);Dm=r[0][0]/B8[0][0];C6=2.**q;s=J['ra0']*C6;t=J['rb0_ga']*C6;C7=t*sqrt(Dm);Dn=s*sqrt(QA_MAX);AS=C7*sqrt(QB_MAX);P=A/6.*(pi**3/6.)*Dn**6/sqrt(Dk);b=.5*(pi**3/6.)*AS**6/sqrt(Dl);C8=log2(P);Do=log2(b);C9=[[r[ae][ak]-r[ae][0]*r[0][ak]/B9 for ak in C(1,6)]for ae in C(1,6)];Dp=_det_mat(C9);BA=2.*AS/sqrt(B9);A8=4.*pi**2/15.*AS**5/sqrt(Dp);c=b/A8;R=B_A_SHELLS;BB=[QA_MAX*(af/R)**(A/3.)for af in C(R+1)];X=P/R;assert X+A<1<<47;Dq=1<<ceil(log2(X+A));CA=(BA+2.)*(2.*Dq)+X*4.;Dr=log2(CA);Ds=log2(b*2.);CB=E((2.**A/B for(D,A,B,C)in w if C=='ll'),default=D);CC=E((2.**A/B for(D,A,B,C)in w if C==d),default=D);Dt=log2(E(8.*(X+A)+32.*X+8.*CB,48.*X+8.*(CC+A))+8.25*P+CA);CD=log2(E(8.*(X+A)+32.*X+8.*CB,48.*X+8.*(CC+A))+8.25*P);AT=B8[0][0]/3.;CE=E(Q(An(A,B))for(A,B)in i['EMB_W']);BC=1.00003;Du=BC*s*s/AT;Dv=BC*t*t/AT;Dw=6.*log2(E(I(D*(B*Du)**C*((A-B)*Dv)**(4-C)for(C,D)in y((3.,14.072,19.9,7.349,A)))for B in(A/5e2 for A in C(501))))-24.*log2(3.);Dx=6.*log2(BC*E(((U_BN254+1)*CE+2.)*(s*s*QA_MAX/AT),t*t*QB_MAX/AT)/3.);BD=c/Aw+A;Dy=0;Dz=ceil(BD*Aw/c*1.005*1e3)/1e3;BE=U(cr2=BD*Aw/c,sb_row=BD,max_sb=Dy,frac_reg_analytic=.6586);CF=AF(pbits(0,j,KH_F_CAP),leaf=7,L=50);CG=AF(pbits(0,n,KH_G_CAP if As is G else As),leaf=8,L=13 if As is G else G,c_leaf=C_SRT_G_LEAF);CH=CF+64.,DK+DL;CI=J[Bn];BF=BG=CJ=D;AU=0;D_={ag:[H(A,DM*rho((ag+48.*q-al+AL)*Az))for al in B0]for ag in By};CK={}
	for AV in Bz:
		E0=AV+12.*q;CL=[];BH=[]
		for CM in B_:E1=H(484,E(10,floor(E0-CM)));E2,E3=Bu(E1,Av-CM);CL.append(E2);BH.append(E3)
		A9=AA=D;AW=[D]*B1
		for(BI,E4,BJ)in C0:A9+=BI*CL[BJ];AA+=BI*BH[BJ];AW[E4]+=BI*BH[BJ]
		CK[AV]=A9,AA,AW
	for(E5,AV,AB,FM)in Ay:
		AU+=AB
		if not AB:continue
		A9,AA,AW=CK[AV];E6=A9*(C_NORM_G+CG)+AA*(C_NORM_F+CH[0])+B2(D_[E5],AW)*CH[1];BF+=AB*E6;BG+=AB*A9;CJ+=AB*AA
	AX=0;BK=CN=D
	for(E7,BL)in CI:AX+=BL;E8=H(484,E(20,V(floor(E7+12.*q))));CO,_=Bu(E8,Av+A);BK+=BL*(CO*C_NORM_G);CN+=BL*CO
	assert AU==J[A1]and AX==J[Ap]and AU+AX==W;BF+=BK;BG+=CN;CP=BF/W;E9=BK/W;CQ=BG/W;EA=CJ/W;EB=C_HIT_LO if Ar else 2.65;BM=A0*(2.**log2_li(f)-2.**log2_li(n));FN=I(b*(H(QB_MAX,A-BB[B])/QB_MAX)**3 for B in C(R));CR=I(A8*(H(QB_MAX,A-BB[B])/QB_MAX)**2.5 for B in C(R));EC=AJ.b1_delta_save_per_cell(sys.modules[__name__],A6,n)if Ar else D;CS=_WKF*B1_PREF_COV*A6*log(f/n)*(EB*1.004+.06)+(BM*b*C_B1_ROW+BM*CR*64.)/AR-EC;ED=D;assert w[0][0]==f and w[-1][1]==o;BN=BO=CT=BP=BQ=CU=AY=D
	for(AZ,Aa,g,Ab)in w:
		AC=A0*(2.**log2_li(Aa)-2.**log2_li(AZ));EE=Da(g,AS,B9)*A6
		for EF in C(R):CV=sqrt(H(QB_MAX,A-BB[EF])/QB_MAX);EG=c*CV;EH=A8*CV**5;CW=EG/M(g)+A;BP+=AC*EH*(CW*C_SB_ENTER_X4 if Ab==d else CW*D8+D9)/AR
		EI=H(c,2.*g);EJ=E(D,c-2.*g);EK=1.005*g;EL=C_B2_INT_X4 if Ab==d else DA+Al;CX=H(A,2.*g/BA);Ac=A-(A-CX*CX)**2.5;EM=(A-Ac)*C_B2_BND1_HIT+Ac*C_B2_BND_HIT+Al;EN=(A-Ac)*C_B2_BND1_MISS+Ac*C_B2_BND_MISS+Al;EO=(EJ*EL+EI*EM+EK*EN)/c;BQ+=EE*log(Aa/AZ)*EO;BN+=AC;BO+=AC;CY=A0*.5*(2.**(2.*Aa)/(2.*Aa*L-A)-2.**(2.*AZ)/(2.*AZ*L-A));EP=C_LIST_B2 if Ab==d else Be;EQ=DD if Ab==d else DE;CU+=.5*AC*P*EP+R*EQ*(CY/g);AY+=B2_PREFIX_MULT*2.*(8.*(P+R)*(.5*AC)+R*8.*(CY/g))
	ER=BO;ES=log(Au/o);CZ=(C_LAT_LEAF+PT_FB_EXC if PT_ON else C_LAT_EMIT_V43B+C_LAT_APPLY)+(A2 if A2 else D);Ca=ES*CZ*Dg+DF;Cb=2.*KH_G_CAP*(2.**log2_li(Au)-2.**log2_li(o));Cc=Cb*(C_LAT_LLL+C_LAT_EXC+(PT_CERT_IR+PT_BUILD_IR if PT_ON else D))/AR;Cd=1.05*A6*(C_SCAN_G8+(DB-1)*(C_SCAN_ZERO-C_SCAN_G8)/c)+1024.*CQ*1.05*A6;Ce=CS+CT+BP+BQ+Ca+Cc+Cd+CP;ET=log2(AR*Ce*OP['work_allow']/h);Cf=A0*.5*(2.**log2_li(f)-2.**log2_li(n));Cg=A0*.5*2.**(2.*f)/(2.*f*L-A);EU=(Cf+A)*(C_LIST_B1*P*x)+DC*Cg*(R*K_B1_PREF)*Bi;Ad=8.*P*Cf*x+8.*Cg*(R*K_B1_PREF)*Bi;EV=Be;EW=M(79151674016786875);EX=M(0x182dba2118d3140c);Ae=P*(A+SC_GUARD_DELTA)**3;EY=Ae*(1e3+12.*log2(Ae))+4.*EW+2.*EX+16.*M(1<<37)+16.*M(1<<27)+5e1*(6.*Ae+8.*b)+2e1*P+6.4*M(1<<36)*P+128.*b+2**36*(4.4*b+8.4*CR);EZ=log2((EU+CU+EY)/h);import stor_bound as Ea;Y=Ea.main();assert Y['ok'];Eb=2.**Y['score_rows_total_log2'];Ec=2.**Y[Bp];Ch=V(ceil(Ec/DG));assert Ch<1<<26;Ed=Y[Bp];FO=2.**Y[Bq];Ee=Y[Bq];Ci=18.*P+256.*2.**f;Cj=2.**Y['SB1_log2'];assert 7.*Cj<=.2*b;Ck=log2(Bf+E(2.**CD,M(1<<56))+l*164e2+1064.*Cj);assert Ci<M(1<<56);Cl=64.*P+65536.*2.**f+8.*2.**o;assert(Ad*8/x+AY)/2**36+Cl<Bf;assert Ck<=Y['list_node_log2'];Ef=6;Cm=Y['machine_total_log2'];Eg=Cm;Eh=Y;Ei=I(51e1*E(Q(An(i[z][B][A][0],i[z][B][A][1]))for B in C(3))for A in C(6));Ej=I((1022. if A==1 else 51e1)*E(Q(An(i[z][B][A][0],i[z][B][A][1]))for B in C(3))for A in C(6));Ek=6.*log2(U_BN254*CE*Ei**2+Ej**2);El=V(floor(Ek/o));Cn=sqrt(E(I((.5 if A>>B&1 else-.5)*C9[B][D]*(.5 if A>>D&1 else-.5)for B in C(5)for D in C(5))for A in C(32)));Em=2.*A8*(1.0015+Cn/C7)**5
	def En(x):
		if B(x)>=12:C,D,E,F,G,H=x[6],x[7],x[8],x[9],x[10],x[11];assert Q(C+(D+E+F+G)/H-A)<.0005;return C*C_UNIT_5+(D*2.*C_UNIT_5+E*3.*C_UNIT_5+F*4.*C_UNIT_5+G*C_GEN_1)/H
		if B(x)==8:return x[7]
		return x[6]*C_UNIT_5+(A-x[6])*C_GEN_1
	def Eo(key):
		F=Db[key];A=B(F)
		def C(idx):B=[ap[idx]for ap in F];C=I(B)/A;G=I((A-C)**2 for A in B)/(A-1)if A>1 else D;return C,sqrt(G/A),E(B)
		K,L,M=C(1);N,O,P=C(2);Q,G,Y=C(3);R,G,S=C(4);T,G,G=C(5);V,G,G=C(6);J=[En(ao)for ao in F];W=I(J)/A;X=H(A[14]for A in F)if B(F[0])>=15 else 0;return U(n=A,N_sw_U=K,N_sw_se=L,N_sw_max=M,N_l_U=N,N_l_se=O,N_l_max=P,k5_ev=Q,k5_exp=W,k5_exp_max=E(J),tab_r=R,tab_r_max=S,unit_frac=T,exp_unit=V,min_ncols=X)
	Ep='m184_z6_61.50_61.50_e0.72_W100_';AD=[Eo(f"{Ep}{ah}")for ah in C(20,29)];assert a(A['n']>=4 for A in AD);Eq,Er,Es,Et,Eu,Ev,Ew,Ex,Af=AD;BR=NC27*m+TEXC_LA/2.**A3 if LLA_ON else NU_SW_CAP*m;Co=BR if LLA_ON else NU_L_CAP*m;BS=K5_CAP*(K5_LLA if LLA_ON else A);BT=TAB_R_CAP;Cp=NU_L_CAP/NU_SW_CAP;Cq=A3+log2(BR);Ag=A3+log2(Co);Ah=5;Ai=(W_ROW*BS*BT+AH)/Cp
	def Cr(lN,m,n):return(A+Bk)*(2.**lN+Ah)*(A/m+A/n)+6.
	def Ey(lN_cap,m,n):A=Cr(lN_cap,m,n);B=m+n;return log2(4.*B**3*ceil(log2(A))*(fft_F(18.*A)+64e1*A)/h)
	K=G
	for Z in(768,)if LLA_ON else(256,512,768,1024):
		for S in(32,)if LLA_ON else(24,32,40,48):
			AE=2.**Ag;Cs=Cr(Ag,Z,S);BU=AE*(A/Z+A/S);Ez=Cs/BU;E_=E(BU,(A+Bk)*(AE+Ah)/S+4.)/BU;F0=(Ai+DS)*(S/Z+A)*Ez+(A/Z+A/S)*(Ai+(Ah*DT+(S-Ah)*DU)+Bj)*E_;F1=(Z+S)*Bj;Ct=log2((AE*AE*F0+AE*F1+DV+DW+DX+DY+2.**20)/h);Cu=Ey(Ag,Z,S);Cv=log2(2.**Ct+2.**Cu)
			if K is G or Cv<K['both']:BV=2.**Cq;BW=ceil(S/8.);Cw=BV*(256.*BW+256.*BW+W_ROW*8.)+2.**55.04*BW+65.*2.**53;Cx=(3.*(Z+S)*Z+Z*S)*Cs*32.;F2=log2(Cw+Cx);K=U(la=Ct,gen=Cu,both=Cv,m=Z,n=S,mem_1own=F2,spmv_mem_log2=log2(Cw),pm_basis_mem_log2=log2(Cx),log2_N_cap=Ag,log2_Ns_cap=Cq,c_prod=Ai)
	BX=j,e,At
	if BX not in Am:Am[BX]=DZ(j,e,At,nf=B5,ng=B6)
	BY=Am[BX];Cy=C2-2.25*B7;BZ=N(Cy,3),N(BY['h'],4)
	if BZ not in _DESC_CACHE:_DESC_CACHE[BZ]=DESC.evaluate_descent(AI[Bl],rho,fft_F,AF,pbits,log2_rels=Cy,peel_h=BY['h'],nf_main=B5,ng_main=B6)
	Ba=_DESC_CACHE[BZ];u=C2+2.25*B7;v=log2(h);Aj=E(j,e);F3=C8+42.8-v;Cz=R_BN254.bit_length();F4=1695744+2**17+147456+184320+49152;F5=92160+2**15+36864+46080+12288;F6=u+log2(48*Cz*F4+12*Cz*F5)-v;F7=log2(2.**53+(PT_GLOBAL_IR if PT_ON else D)+3.*l*Ae+2.**(u+30.4)+2.**(u+14.)+(2.**(u+28.)if LLA_ON else D)+(Ef+1)*l*413e1+l*(257552.+8192.*CG)+(1<<50)*CF+M(DI)*M(1<<36)*l*413e1+BM*l*413e1)-v;F8=log2(1.25506*2.**Aj/(Aj*L));F9=F8+23.4-v;FA=AJ.split_ntt3_parts(sys.modules[__name__],2.**u,j,e);C_=log2(FA['total_ir'])-v+A;FB=ceil(Aj-6.+.53);FC=I(2.**(Aj-6.+.53-A-1)*(fft_F(2.**(A+1))+3.*2.**(A+1))for A in C(FB));FD=log2(2.*FC)-v;Bb=U(collection=ET,chi_setup=65.712288,la=K['la'],gen=K['gen'],lists=EZ,two_logs=Ba['two_log'],canon=F3,sm=F6,misc=F7,fb=F9,split=C_,gtree=FD,cert=50.000001,verify=11.000001,wrapper=53.712288)
	if track=='bn254-pairing-inversion':Bb['output']=log2(G2_OUTPUT_BOUND/h)
	Bc=log2(I(2.**A for A in Bb.values()));BV=ceil(2.**K[Br]);FE=2.**Bc/2.**DR;D0=U(main_matrix=4.*BV/2.**64,filt_cap_h3_given_H3=DQ,split_abort=2.**C_/(2.**DN-2.**DO/h),descent_wald_retry=D,descent_node_abort=D,descent_repeat_ideal=D,unified_run_meter=FE,storage_cap_abort=DJ,miller_rabin_64r=2.**69*2.**-128.,group_order_boot_and_sm=(8.+2.*Ba['e_trials_per_log'])/R_BN254,numerical_pure_int=D);D1=A-I(D0.values());D2=f"{H(floor(D1*1e6)/1e6,DP):.6f}";Bd=F(D2);Ak=V(-(-(F(39,100)*2**64)//Bd));FF=F(Ak,2**64)*Bd>=F(39,100);FG=F(Ak-1,2**64)*Bd<F(39,100);D3=ceil(Bc*1e9)/1e9;D4=log2(A+M(F(Ak,2**64))*2.**D3);D5=ceil(D4*1e9);FH=f"{D5/1e9:.9f}";FI=a(AD[B][A]>AD[B+1][A]for A in('N_sw_U','N_l_U','k5_exp')for B in C(8));FJ=Af['N_sw_max']<NU_SW_CAP*m and Af['N_l_max']<NU_L_CAP*m and Af['k5_exp_max']<BS and E(A['tab_r_max']for A in AD)<BT;return U(track=track,use_db=Aq,c_tail_z6=AH,cap_scale=m,strict_per_owner=Ar,swar_probe=A2,op=OP,pt=U(on=PT_ON,pt_fb_exc=PT_FB_EXC if PT_ON else D,pt_global_ir=PT_GLOBAL_IR if PT_ON else D,c_lat_use=CZ,c_lat_emit=C_LAT_EMIT if PT_ON else C_LAT_EMIT_V43B,c_lat_apply=C_LAT_APPLY,c_lat_leaf=C_LAT_LEAF if PT_ON else C_LAT_EMIT_V43B+C_LAT_APPLY,cert_ir=PT_CERT_IR if PT_ON else D,build_ir=PT_BUILD_IR if PT_ON else D,build_ir_log2=log2(PT_BUILD_IR),n_orb_lat=Cb,lat_lll_c_per_orbit=C_LAT_LLL+C_LAT_EXC+(PT_CERT_IR+PT_BUILD_IR if PT_ON else D)),pool=U(n_total=J[Bo],n_in=W,n_prim_in=J[A1],frac_reg=A4,frac_reg_analytic=BE['frac_reg_analytic']),geometry=U(log2_U=A3,log2_V=C5,se_bits=B7,log2_Vb=Do,log2_Na=C8,log2_Na_shell=log2(X),log2_N5=log2(A8),L_b_avg=c,L_max=BA,Ra=N(s,2),Rb_GA=N(t,2),sb_row=BE['sb_row'],max_sb=BE['max_sb'],cr2_tight=Dz,p_occ_b2=BO/E(BN,A),p_ent_b2=ER/E(BN,A),nf_rel=B5,ng_rel=B6,r_g=Dj,sf=CQ,s0=EA,score_arena_log2=Dr,vb_score_log2=Ds,b0_blocks=R,n_sc_nodes_exp=Bg,sc_block_exp=Bh,sc_w_subband=DH,sc_guard_delta=SC_GUARD_DELTA,stor_bound=Eh,live_score_total_log2=N(log2(Eb),4),cs_blocks_max=Ch,lat_8b_direct_ptr=O,delta5_blichfeldt=N(Cn,6),n5_disc_ub_log2=N(log2(Em),4),list_1node_ram_log2=Ck,b1_lists_log2=N(log2(Ad/2**36),4),lists_8x_log2=N(log2((Ad*8/x+AY)/2**36),4),lists_node_max_log2=N(log2((Ad*8/x+AY)/2**36+Cl),4),b1_build_scratch_log2=N(log2(Ci),4),emit_1own_ram_log2=Ee,score_1node_ram_log2=Ed,det_max_pf_lat=El,b2_1p_ram_log2=Dt,b2_build_1p_ram_log2=CD,avg_b2_list_ir=N(EV,3),ram_log2=Cm,peak_global_ram_cap_log2=Eg,sup_nf_bits=N(Dw,2),sup_ng_bits=N(Dx,2)),per_cell=U(total=Ce,rep_b0=ED,rep_lo=CS,probe=CT,occ_sb=BP,hi_ap=BQ,lat_b2=Ca,lat_lll=Cc,scan_zero=Cd,cof=CP,cof_nonprim=N(E9,6),n_prim_exact=AU,n_nonprim_exact=AX,nonprim_bins_n=B(CI)),la=U(bw=[K['m'],K['n']],log2_N_cap=K['log2_N_cap'],log2_Ns_cap=K[Br],W_row=W_ROW,nu_sw_cap=BR,nu_l_cap=Co,k5_cap=BS,tab_r_cap=BT,pf_rl_cap=Cp,c_prod=Ai,spmv_mem_log2=N(K['spmv_mem_log2'],4),pm_basis_mem_log2=N(K['pm_basis_mem_log2'],4),mem_1own_log2=K[Bs],fits_1own=K[Bs]<64.,h3_mode='fixed_instance_deterministic_4cap',ladder_all4_monotonic=FI,all_runs_28_below_4caps=FJ,filtsim_cur20=Eq,filtsim_cur21=Er,filtsim_cur22=Es,filtsim_cur23=Et,filtsim_cur24=Eu,filtsim_cur25=Ev,filtsim_cur26=Ew,filtsim_cur27=Ex,filtsim_cur28=Af,peel=BY),descent=Ba,rows=Bb,trim_ir_log2=u+28.,exp_run=Bc,exp_run_decl=f"{D3:.9f}",fail_terms=D0,succ_exact=D1,succ_decl=D2,K_gate=Ak,gate_passes=FF,gate_prev_fails=FG,claim_exact=D4,claim_decl=FH,score_numerator=str(D5))
_OH_U=(1,-1,-1,-1,1,0),(-1,1,1,1,0,0),(1,1,0,1,-1,0),(0,0,1,0,0,1),(-1,0,0,0,0,1),(1,0,0,-1,0,-1),(0,-1,0,-1,1,0),(1,0,0,0,1,-1),(0,1,1,1,0,1),(-3,0,2,1,0,1),(0,1,0,0,0,0),(1,0,-1,0,0,0),(-1,2,2,1,0,1),(-1,0,0,1,0,2),(-1,-1,-1,0,-1,1)
_OHM_CACHE=G
c=G
def AJ(m):
	B=[aa[:]for aa in m];H=1;D=1
	for A in C(6):
		if B[A][A]==0:
			for E in C(A+1,6):
				if B[E][A]:B[A],B[E]=B[E],B[A];D=-D;break
			else:return 0
		I=B[A][A];J=B[A]
		for K in C(A+1,6):
			F=B[K];L=F[A]
			for G in C(A+1,6):F[G]=(F[G]*I-L*J[G])//H
		H=I
	return D*B[5][5]
def is_ideal_primitive_oh(a,b):
	if T.gcd(*a,*b)>1:return K
	global _OHM_CACHE,c
	if c is G:
		if _OHM_CACHE is G:
			F=[[[1 if ap==ak else 0 for ap in C(6)]for ak in C(6)]for af in C(6)];R=0
			for E in C(1,6):
				F[E][0]=F[0][E]
				for S in C(E,6):F[E][S]=F[S][E]=_OH_U[R];R+=1
			_OHM_CACHE=F
		F=_OHM_CACHE;c=[[[F[al][ac][ag]for al in C(6)]for ag in C(6)]for ac in C(6)]
	U=c;V=T.sumprod;N=[[V(a,U[aa][ad])for ad in C(6)]for aa in C(6)];P=[[V(b,U[ab][ae])for ae in C(6)]for ab in C(6)];M=T.gcd(AJ(N),AJ(P))
	if M==1:return O
	A=[[am%M for am in ah]for ah in N+P]+[[M if an==ai else 0 for an in C(6)]for ai in C(6)]if M else N+P;Y=B(A);J=0
	for I in C(6):
		while O:
			L=[aj for aj in C(J,Y)if A[aj][I]]
			if not L:return K
			if B(L)==1:
				D=L[0]
				if D!=J:A[J],A[D]=A[D],A[J]
				if Q(A[J][I])!=1:return K
				J+=1;break
			D=H(L,key=lambda i:Q(A[i][I]))
			if A[D][I]<0:A[D]=[-ao for ao in A[D]]
			Z=A[D][I];d=A[D]
			for E in L:
				if E==D:continue
				W=A[E];e=W[I]//Z
				for X in C(I,6):W[X]-=e*d[X]
	return O
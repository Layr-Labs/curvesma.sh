A1=b'\x00\x00'
d='seed0'
c='s16'
b='little'
a='boot_s16.b85'
A0=bytearray
V='n'
U=False
n=set
X=dict
N=zip
T=enumerate
J='rec'
O=None
L=int
K=list
I=abs
G=.0
F=1.
E=max
D=round
C=len
B=sum
A=range
import base64 as e,json,lzma,math as H,random
from math import log2
from pathlib import Path
import model_bn254_v6 as M,descent_bn254_v6 as DESC
f=20261123
W=[[1,0,1,-1,-3,-3],[0,1,-1,0,-1,-1],[0,0,-1,1,1,3],[0,0,-1,-1,0,2],[0,0,1,0,1,-1],[0,0,0,1,1,2]]
P=[[0x13482b02e61085c8075c496b702953db8adb7d3406418007134e068338224da0,0x6f4f657a5d7d146a56a61a06011ee7443242e98d5dbde5a98c00952f77534c8,0x1a5a1b16367b977b2b01214f4018db8cedf9ebbf96c0bf4830bce5241e9af460,0x1f404dd1e5b65a1532934f7fcb09469589466a249aad1e9fc1e272dc09a5120e,0x2ef6ac4e331c03e9864bc9e7cebc0f061761dfef57875e074a8ad1c7ff506ac5,0x10ffb558754f7da7cd6107c7883f9ecf13ced8c611ff6d6e69d96d2ce4f8cd59],[0xa197ee65f48459da3dfdc12fae53f5965880ad4087e407ce118bac14f75fe38,0x16efb7210e98f80f88c5eae967b5491aa15931bbfee9e7cd8c6c03a57b51d218,0xb7ce411b9d0f2a9a9e636a7afe08d25625e807cc6de2ef147e9d1c8b0db2139,0x10761b9d6a5cfadd2f4e6f1db3afcaab25e9b2f5376ef8a44c81afd86b773bfe,0x17eb64652a8b4fd40a590d99908f4bcc97d050c934b447a1c9ce90f3730d46a0,0x18368a88b76285239e2274abe3f65a8410557f8de9928d5e846d1703c377fc69]]
R={}
def ctx(u=O,full=U):
	C=full;u=M.U_BN254 if u is O else u
	if u in R and not C:return R[u]
	A=36*u**4+36*u**3+24*u**2+6*u+1;E=36*u**4+36*u**3+18*u**2+6*u+1;B,H=M._pmul_h,M._pinv_h;U=M.W_COEFFS;J=H([-((u+1)*j+(2 if i==0 else 0))%A for(i,j)in T(U)],A)
	def F(x,y):return[(k+l)%A for(k,l)in N(B(x[0],y[0],A),B(B(x[1],y[1],A),J,A))],[(m+o)%A for(m,o)in N(B(x[0],y[1],A),B(x[1],y[0],A))]
	def V(x):C=H([(p-q)%A for(p,q)in N(B(x[0],x[0],A),B(B(x[1],x[1],A),J,A))],A);return B(x[0],C,A),B([-r%A for r in x[1]],C,A)
	def G(x,e):
		A=[1,0,0,0,0,0],[0]*6
		for B in bin(e)[2:]:
			A=F(A,A)
			if B=='1':A=F(A,x)
		return A
	L=[1,0,0,0,0,0],[0]*6;W=lambda:G(F(([3,1,0,0,0,0],[A-1,A-2,0,0,0,0]),V(([3,1,0,0,0,0],[1,2,0,0,0,0]))),(A**6+1)//E);D=(K(P[0]),K(P[1]))if u==M.U_BN254 and not C else W();assert D!=L and G(D,E)==L
	if C and u==M.U_BN254:assert D==(K(P[0]),K(P[1]))
	Y=M.build_exact_int_geom(M.load_data()['t4_geom'])[0];from fractions import Fraction as Z;Q=Z(0xa8ccdee3ccdbd36d*I(u),M.U_BN254);S=X(u=u,p=A,r=E,pw=G,gg=D,G6=Y,i3=pow(3,-1,A),ref=6*log2(A)+6*log2(I(u)),wa=Q.numerator,wb=Q.denominator)
	if not C:R[u]=S
	return S
def lattice(c,j):return lattice_e(c,random.Random(f+j).randrange(1,c['r']))
def lattice_e(c,e):
	C=c['p'];G,K=M._pmul_h,M._pinv_h;L=M.ZK_NUM;H=c['pw'](c['gg'],e);I=G(H[1],K(H[0],C),C);E=[[C if j==6+i else 0 for j in A(12)]for i in A(6)]
	for J in A(6):
		N=G(I,[l*c['i3']%C for l in L[J]],C);D=[1 if k==J else 0 for k in A(6)]
		for O in A(6):F=B(W[O][A]*N[A]for A in A(6))%C;D.append(F-C if F>C//2 else F)
		E.append(D)
	DESC.bareiss_lll_12d(E,c['G6'],c['wa'],c['wb'],U,11,40,(9,10))
	for D in E:assert M._pmul_h(I,Y(D[:6],C),C)==Y(D[6:],C)
	return E
def Y(v,p):return[B(v[A]*M.ZK_NUM[A][i]for A in A(6))*pow(3,-1,p)%p for i in A(6)]
def cands(B):return[B[i]for i in A(12)]+[[m+l*o for(m,o)in N(B[j],B[k])]for j in A(12)for k in A(j+1,12)for l in(1,-1)]
def ng(c,v):D,J,C=M.ZK_NUM,M.W_COEFFS,M._mh_int;E=[B(v[A]*D[A][i]for A in A(6))for i in A(6)];F=[B(v[6+A]*D[A][j]for A in A(6))for j in A(6)];G=C(E,E);K=C(J,G);L=C(F,F);H=I(M._norm_Kh_exact([(c['u']+1)*K[k]+2*G[k]+L[k]for k in A(6)]));assert H%9**6==0;return H//9**6
Q=[O]
def s16(N):
	if Q[0]is O:
		B=K(M._primes_upto(1<<16))
		while C(B)>1:B=[B[i]*B[i+1]if i+1<C(B)else B[i]for i in A(0,C(B),2)]
		Q[0]=B[0]
	E=N;D=H.gcd(Q[0],E);F=1
	while D>1:
		while E%D==0:E//=D;F*=D
		D=H.gcd(D,E)
	return log2(F)
def sample(j,u=O):
	K=ctx(u);L=n();B=[];O=0
	for(R,G)in T(cands(lattice(K,j))):
		if E(I(A)for A in G[:6]).bit_length()>=128 or E(I(A)for A in G[6:]).bit_length()>=160:continue
		F=ng(K,G)
		if F==0:continue
		S=F not in L;L.add(F)
		if not M.is_ideal_primitive_oh(G[:6],G[6:]):continue
		O+=1
		if S:B.append((R,D(100*(log2(F)-K['ref'])),D(100*s16(F)),F))
	P=0;N=0
	for Q in A(C(B)):
		for U in A(Q+1,C(B)):
			J=H.gcd(B[Q][3],B[U][3])
			if J>1:
				J=g(J)
				if J>1:P+=1;N=E(N,J.bit_length())
	return[C(L),O,[i[0]for i in B],[k[1]for k in B],[P,N],[l[2]for l in B]]
def g(g):
	s16(1);A=H.gcd(Q[0],g)
	while A>1:
		while g%A==0:g//=A
		A=H.gcd(A,g)
	return g
def gate_norm_bits():
	from fractions import Fraction as D;F=[[D(k,3)for k in i]for i in M.ZK_NUM]
	def J(a,b):
		B=[D(0)]*11
		for F in A(6):
			for C in A(6):B[F+C]+=a[F]*b[C]
		for E in A(10,5,-1):
			G=B[E]
			if G:
				for C in A(6):B[E-6+C]-=G*M.H_COEFFS[C]
				B[E]=D(0)
		return B[:6]
	G=lambda q:[L(B(W[o][A]*(3*q[A])for A in A(6)))//3 for o in A(6)];C=[[G(J(F[j],F[l]))for l in A(6)]for j in A(6)];K=E(B(I(C[B][E][D])for B in A(6)for E in A(6))for D in A(6));N=H.sqrt(E(B(B(I(C[A][D][E])for A in A(6))**2 for E in A(6))for D in A(6)));O=G([D(m)for m in M.W_COEFFS]);P=E(B(I(B(O[A]*C[A][E][D]for A in A(6)))for E in A(6))for D in A(6));return 6.*(log2(N)+log2(((M.U_BN254+1)*P+2.)*2.**256+2.**320)+log2(K))
def h(name,raw):
	R,X,F=lzma.decompress(e.b85decode(raw)).partition(b'\n');E=json.loads(R);C=E[V];K=F[:C];H=B(K)
	if name==a:
		M=F[C:C+2*H];D=0;N=[]
		for G in K:N.append([L.from_bytes(M[2*(D+i):2*(D+i)+2],b)+E['v0']for i in A(G)]);D+=G
		return{c:N}
	S=F[C:2*C];T=F[2*C:3*C];U=F[3*C:3*C+H];M=F[3*C+H:3*C+3*H];D=0;O=[]
	for I in A(C):
		G=K[I];P=0;Q=[]
		for W in U[D:D+G]:P+=W;Q.append(P)
		O.append([S[I]+E['a0'],T[I]+E['a1'],Q,[L.from_bytes(M[2*(D+j):2*(D+j)+2],b)+E['v0']for j in A(G)],E['xy'][I]]);D+=G
	return{d:E[d],V:C,J:O}
def Z(name):return h(name,(Path(__file__).resolve().parent/name).read_bytes())
S=[O]
def load():
	if S[0]is O:
		A=Z('boot_dp.b85');D=Z(a)[c];assert C(D)==C(A[J])==A[V]
		for(B,E)in N(A[J],D):assert C(E)==C(B[2])==C(B[3]);B.append(E)
		S[0]=A
	return S[0]
def replay(js,full=U):
	if full:ctx(full=True)
	B=load();A=[i for i in js if sample(i)!=B[J][i]];return X(replayed=K(js),mismatches=A,generator_recomputed=full,ok=not A)
def alpha16(S=O):S=S or load();L=[B(i[5])/1e2 for i in S[J]];F=[C(j[5])for j in S[J]];A=C(L);I=B(L)/B(F);R=B(F)/A;T=H.sqrt(B((A-I*B)**2 for(A,B)in N(L,F))/(A-1)/A)/R;O=B(log2(A)/(A-1)for A in M._primes_upto(1<<16));W=.08;G=[[2.**(W*o/1e2)for o in k[5]]for k in S[J]];K=B(F);P=B(map(B,G))/K;U=[(B(l)**2-B(A*A for A in l))/P**2 for l in G];Q=[C(m)*(C(m)-1)for m in G];V=B(U)/B(Q);Y=H.sqrt(B((A-V*B)**2 for(A,B)in N(U,Q))/(A-1)/A)/(B(Q)/A);Z=((K*P)**2-B(B(A)**2 for A in G))/P**2/(K*K-B(C(A)**2 for A in G));a=B(A[4][0]for A in S[J]);b=B(C(A[5])*(C(A[5])-1)//2 for A in S[J]);return X(kept_pairs=b,kept_pairs_sharing_a_prime_above_2_16=a,max_shared_rough_bits=E(A[4][1]for A in S[J]),kept=B(F),lattices=A,mean_kept=D(R,3),mean_s16=D(I,4),se_lattice=D(T,4),e_uniform=D(O,4),alpha_emp=D(O-I,4),alpha_emp_plus_3se=D(O-I+3*T,4),pair_moment_k008_within_lattice=D(V,4),pair_moment_k008_within_se=D(Y,4),pair_moment_k008_between_lattices=D(Z,4))
def rough_med_audit(js,max_p=1<<20):
	Q=max_p;T=A0(b'\x01')*(Q+1);T[:2]=A1
	for N in A(2,L(Q**.5)+1):
		if T[N]:T[N*N::N]=b'\x00'*((Q-N*N)//N+1)
	W=[o for o in A((1<<16)+1,Q+1)if T[o]];J=K(W)
	while C(J)>1:J=[J[p]*J[p+1]if p+1<C(J)else J[p]for p in A(0,C(J),2)]
	i=J[0];c=B(log2(A)/(A-F)for A in W);d=F
	for Y in W:d*=(F-F/Y)/(F-Y**.08/Y)
	e=ctx();R=[];Z=[];a=0;f=G
	for j in js:
		k=lattice(e,j);g=n();O=[]
		for S in cands(k):
			if E(I(A)for A in S[:6]).bit_length()>=128 or E(I(A)for A in S[6:]).bit_length()>=160:continue
			U=ng(e,S)
			if U==0:continue
			l=U not in g;g.add(U)
			if not l or not M.is_ideal_primitive_oh(S[:6],S[6:]):continue
			V=U;P=H.gcd(i,V);h=1
			while P>1:
				while V%P==0:V//=P;h*=P
				P=H.gcd(P,V)
			O.append(log2(h))
		a+=C(O);f+=B(O);R.append(B(O)/C(O));Z.append(B(2.**(.08*A)for A in O)/C(O))
	b=f/a;m=H.sqrt(B((A-b)**2 for A in R)/(C(R)-1)/C(R))if C(R)>1 else G;return X(lattices=K(js),kept_norms=a,band_bits=[16,L(D(log2(Q)))],mean=D(b,4),se=D(m,4),e_unif=D(c,4),alpha_med=D(c-b,4),moment_008_ratio=D(B(Z)/C(Z)/d,4),ok=True)
def joint_ea_mini_audit(u_val=-3,n_lat=60):
	c=u_val;import functools as A2;from math import log1p as A3,floor as o,ceil,sqrt;J,d,p,U,Y,q,W=(14.,2e1,72,4,96,16,64.)if c==-3 else(16.,21.,72,4,96,16,8e1)
	def Z(a,b):
		A=[G]*(C(a)+C(b)-1)
		for(D,B)in T(a):
			if B:
				for(E,F)in T(b):A[D+E]+=B*F
		return A
	P=(d-J)/Y;r=[A3(P/(J+aa*P))for aa in A(Y)];AK=[r[ab]if ab<q else G for ab in A(Y)];AL=[r[ac]if q<=ac<p else G for ac in A(Y)];AM=[r[ad]if ad>=p else G for ad in A(Y)];s=[[F]];e=[[F]]
	for V in A(1,U+1):s.append([ah/V for ah in Z(s[-1],AK)]);e.append([ai/V for ai in Z(e[-1],AL)])
	A4={}
	for N in A(U+2):
		t=o((W+1e-09-N*J)/P-N)
		if t<0:continue
		A5=e[N][:t+1]if N<=U else[];a=Z(e[N-1],AM)[:t+1]if N>=1 else[]
		for f in A(U-N+2):
			Q=K(A5)if A5 and f<=U-N else[]
			if a and f<=U-N+1:
				if not Q:Q=K(a)
				else:AN=E(C(Q),C(a));Q=[(Q[am]if am<C(Q)else G)+(a[am]if am<C(a)else G)for am in A(AN)]
			if not Q or not any(Q):continue
			u=Z(s[f],Q);g=A4.setdefault(f+N,[])
			if C(g)<C(u):g.extend([G]*(C(u)-C(g)))
			for(h,AO)in T(u):g[h]+=AO
	R=J+q*P;A6=J+p*P
	@A2.lru_cache(O)
	def AP(m):
		E=G
		for(D,H)in A4.items():
			I=min(C(H)-1,o((m-D*P-D*J)/P))+1
			if I<=0:continue
			E+=B(H[A]*(M.rho((m-D*J-A*P+DESC.ALPHA_G)/J)if m-D*J-A*P+DESC.ALPHA_G>0 else F)for A in A(I))
		return E
	A7=96;i=(W-R)/A7;AQ=[A3(i/(R+ae*i))for ae in A(A7)];v=[F];A8=[]
	for V in A(1,L(W//R)+1):v=[ag/V for ag in Z(v,AQ)];A8.append([(aj,ak)for(aj,ak)in T(v)if ak and V*R+aj*i<=W+1e-09])
	@A2.lru_cache(O)
	def AR(m):
		A=M.rho(m/R)if m>0 else F
		for(B,H)in T(A8,1):
			C=G
			for(I,J)in H:D=E(m-B*R-(I+B)*i,G)/R;C+=J*(F if D<=F else M.rho(D))
			A+=C
		return min(F,A)
	j=1<<L(d);k=A0(b'\x01')*(j+1);k[:2]=A1
	for S in A(2,L(j**.5)+1):
		if k[S]:k[S*S::S]=b'\x00'*((j-S*S)//S+1)
	A9=[af for af in A(2,j+1)if k[af]];AA=[]
	for h in A(0,C(A9),256):
		AB=A9[h:h+256];AC=1
		for S in AB:AC*=S
		AA.append((AC,AB))
	AD=ctx(c);AS=random.Random(20261123+c);AT,AU,AV,AW,AX=2.**R,2.**J,2.**A6,2.**d,2.**W
	def AE(nums):
		D=1
		for B in nums:D*=B
		E=[]
		for(Q,R)in AA:
			F=H.gcd(D,Q)
			if F>1:E+=[an for an in R if F%an==0]
		I=J=0;K=L=G
		for B in nums:
			M=log2(B);K+=AR(o(M));L+=AP(ceil(M));N=C=B;O=P=0
			for A in E:
				if C%A==0:
					while C%A==0:
						C//=A
						if A<=AT:N//=A
						if A>AU:
							if A<=AV:O+=1
							elif A<=AW:P+=1
			if N<=AX:
				I+=1
				if C==1 and O<=U and P<=1:J+=1
		return I,K,J,L
	w=[];x=[];y=[];z=[];AF=AG=AH=0
	for V in A(n_lat):
		AY=lattice_e(AD,V+1);AI=n();l=[]
		for b in cands(AY):
			if E(I(A)for A in b[:6]).bit_length()>=128 or E(I(A)for A in b[6:]).bit_length()>=160:continue
			m=ng(AD,b)
			if m==0:continue
			AZ=m not in AI;AI.add(m)
			if AZ and M.is_ideal_primitive_oh(b[:6],b[6:]):l.append(m)
		AH+=C(l);Aa,Ab,Ac,Ad=AE(l);w.append(Aa);x.append(Ab);y.append(Ac);z.append(Ad);Ae,_,Af,_=AE([AS.randrange(E(2,L(al*.99)),L(al*1.01)+2)for al in l]);AF+=Ae;AG+=Af
	AJ=B(w)/B(x);Ag=B(y)/B(z);return X(u=c,lattices=n_lat,kept=AH,x_ea=D(R,2),tau=W,LB=J,L_mid=D(A6,2),L_top=d,surv_obs=B(w),surv_ctrl=AF,surv_mod=D(B(x),1),surv_vs_mod=D(AJ,4),use_obs=B(y),use_ctrl=AG,use_mod=D(B(z),1),use_vs_mod=D(Ag,4),ok=AJ<=F)
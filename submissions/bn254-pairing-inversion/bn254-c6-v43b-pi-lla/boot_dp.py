import base64,json,lzma,math,random
from math import log2
from pathlib import Path
import model_bn254_v6 as M,descent_bn254_v6 as DESC
SEED0=20261123
ZK3I=[[1,0,1,-1,-3,-3],[0,1,-1,0,-1,-1],[0,0,-1,1,1,3],[0,0,-1,-1,0,2],[0,0,1,0,1,-1],[0,0,0,1,1,2]]
GG_PROD=[[0x13482b02e61085c8075c496b702953db8adb7d3406418007134e068338224da0,0x6f4f657a5d7d146a56a61a06011ee7443242e98d5dbde5a98c00952f77534c8,0x1a5a1b16367b977b2b01214f4018db8cedf9ebbf96c0bf4830bce5241e9af460,0x1f404dd1e5b65a1532934f7fcb09469589466a249aad1e9fc1e272dc09a5120e,0x2ef6ac4e331c03e9864bc9e7cebc0f061761dfef57875e074a8ad1c7ff506ac5,0x10ffb558754f7da7cd6107c7883f9ecf13ced8c611ff6d6e69d96d2ce4f8cd59],[0xa197ee65f48459da3dfdc12fae53f5965880ad4087e407ce118bac14f75fe38,0x16efb7210e98f80f88c5eae967b5491aa15931bbfee9e7cd8c6c03a57b51d218,0xb7ce411b9d0f2a9a9e636a7afe08d25625e807cc6de2ef147e9d1c8b0db2139,0x10761b9d6a5cfadd2f4e6f1db3afcaab25e9b2f5376ef8a44c81afd86b773bfe,0x17eb64652a8b4fd40a590d99908f4bcc97d050c934b447a1c9ce90f3730d46a0,0x18368a88b76285239e2274abe3f65a8410557f8de9928d5e846d1703c377fc69]]
_CTX={}
def ctx(u=None,full=False):
	C=full;u=M.U_BN254 if u is None else u
	if u in _CTX and not C:return _CTX[u]
	A=36*u**4+36*u**3+24*u**2+6*u+1;E=36*u**4+36*u**3+18*u**2+6*u+1;B,H=M._pmul_h,M._pinv_h;N=M.W_COEFFS;I=H([-((u+1)*C+(2 if B==0 else 0))%A for(B,C)in enumerate(N)],A)
	def F(x,y):return[(B+C)%A for(B,C)in zip(B(x[0],y[0],A),B(B(x[1],y[1],A),I,A))],[(B+C)%A for(B,C)in zip(B(x[0],y[1],A),B(x[1],y[0],A))]
	def O(x):C=H([(B-C)%A for(B,C)in zip(B(x[0],x[0],A),B(B(x[1],x[1],A),I,A))],A);return B(x[0],C,A),B([-B%A for B in x[1]],C,A)
	def G(x,e):
		A=[1,0,0,0,0,0],[0]*6
		for B in bin(e)[2:]:
			A=F(A,A)
			if B=='1':A=F(A,x)
		return A
	J=[1,0,0,0,0,0],[0]*6;P=lambda:G(F(([3,1,0,0,0,0],[A-1,A-2,0,0,0,0]),O(([3,1,0,0,0,0],[1,2,0,0,0,0]))),(A**6+1)//E);D=(list(GG_PROD[0]),list(GG_PROD[1]))if u==M.U_BN254 and not C else P();assert D!=J and G(D,E)==J
	if C and u==M.U_BN254:assert D==(list(GG_PROD[0]),list(GG_PROD[1]))
	Q=M.build_exact_int_geom(M.load_data()['t4_geom'])[0];from fractions import Fraction as R;K=R(0xa8ccdee3ccdbd36d*abs(u),M.U_BN254);L=dict(u=u,p=A,r=E,pw=G,gg=D,G6=Q,i3=pow(3,-1,A),ref=6*log2(A)+6*log2(abs(u)),wa=K.numerator,wb=K.denominator)
	if not C:_CTX[u]=L
	return L
def lattice(c,j):return lattice_e(c,random.Random(SEED0+j).randrange(1,c['r']))
def lattice_e(c,e):
	A=c['p'];E,I=M._pmul_h,M._pinv_h;J=M.ZK_NUM;F=c['pw'](c['gg'],e);G=E(F[1],I(F[0],A),A);C=[[A if C==6+B else 0 for C in range(12)]for B in range(6)]
	for H in range(6):
		K=E(G,[B*c['i3']%A for B in J[H]],A);B=[1 if A==H else 0 for A in range(6)]
		for L in range(6):D=sum(ZK3I[L][A]*K[A]for A in range(6))%A;B.append(D-A if D>A//2 else D)
		C.append(B)
	DESC.bareiss_lll_12d(C,c['G6'],c['wa'],c['wb'],False,11,40,(9,10))
	for B in C:assert M._pmul_h(G,_zkp(B[:6],A),A)==_zkp(B[6:],A)
	return C
def _zkp(v,p):return[sum(v[A]*M.ZK_NUM[A][B]for A in range(6))*pow(3,-1,p)%p for B in range(6)]
def cands(B):return[B[A]for A in range(12)]+[[A+D*B for(A,B)in zip(B[A],B[C])]for A in range(12)for C in range(A+1,12)for D in(1,-1)]
def ng(c,v):B,G,A=M.ZK_NUM,M.W_COEFFS,M._mh_int;C=[sum(v[A]*B[A][C]for A in range(6))for C in range(6)];D=[sum(v[6+A]*B[A][C]for A in range(6))for C in range(6)];E=A(C,C);H=A(G,E);I=A(D,D);F=abs(M._norm_Kh_exact([(c['u']+1)*H[A]+2*E[A]+I[A]for A in range(6)]));assert F%9**6==0;return F//9**6
_PP=[None]
def s16(N):
	if _PP[0]is None:
		A=list(M._primes_upto(1<<16))
		while len(A)>1:A=[A[B]*A[B+1]if B+1<len(A)else A[B]for B in range(0,len(A),2)]
		_PP[0]=A[0]
	C=N;B=math.gcd(_PP[0],C);D=1
	while B>1:
		while C%B==0:C//=B;D*=B
		B=math.gcd(B,C)
	return log2(D)
def sample(j,u=None):
	E=ctx(u);F=set();A=[];H=0
	for(K,C)in enumerate(cands(lattice(E,j))):
		if max(abs(A)for A in C[:6]).bit_length()>=128 or max(abs(A)for A in C[6:]).bit_length()>=160:continue
		B=ng(E,C)
		if B==0:continue
		L=B not in F;F.add(B)
		if not M.is_ideal_primitive_oh(C[:6],C[6:]):continue
		H+=1
		if L:A.append((K,round(100*(log2(B)-E['ref'])),round(100*s16(B)),B))
	I=0;G=0
	for J in range(len(A)):
		for N in range(J+1,len(A)):
			D=math.gcd(A[J][3],A[N][3])
			if D>1:
				D=_rough16(D)
				if D>1:I+=1;G=max(G,D.bit_length())
	return[len(F),H,[A[0]for A in A],[A[1]for A in A],[I,G],[A[2]for A in A]]
def _rough16(g):
	s16(1);A=math.gcd(_PP[0],g)
	while A>1:
		while g%A==0:g//=A
		A=math.gcd(A,g)
	return g
def gate_norm_bits():
	from fractions import Fraction as C;B=[[C(A,3)for A in A]for A in M.ZK_NUM]
	def E(a,b):
		A=[C(0)]*11
		for E in range(6):
			for B in range(6):A[E+B]+=a[E]*b[B]
		for D in range(10,5,-1):
			F=A[D]
			if F:
				for B in range(6):A[D-6+B]-=F*M.H_COEFFS[B]
				A[D]=C(0)
		return A[:6]
	D=lambda q:[int(sum(ZK3I[B][A]*(3*q[A])for A in range(6)))//3 for B in range(6)];A=[[D(E(B[A],B[C]))for C in range(6)]for A in range(6)];F=max(sum(abs(A[C][D][B])for C in range(6)for D in range(6))for B in range(6));G=math.sqrt(max(sum(sum(abs(A[D][B][C])for D in range(6))**2 for C in range(6))for B in range(6)));H=D([C(A)for A in M.W_COEFFS]);I=max(sum(abs(sum(H[B]*A[B][D][C]for B in range(6)))for D in range(6))for C in range(6));return 6.*(log2(G)+log2(((M.U_BN254+1)*I+2.)*2.**256+2.**320)+log2(F))
def _ld(name):return json.loads(lzma.decompress(base64.b85decode((Path(__file__).resolve().parent/name).read_bytes())).decode())
_S=[None]
def load():
	if _S[0]is None:
		A=_ld('boot_dp.b85');C=_ld('boot_s16.b85')['s16'];assert len(C)==len(A['rec'])==A['n']
		for(B,D)in zip(A['rec'],C):assert len(D)==len(B[2])==len(B[3]);B.append(D)
		_S[0]=A
	return _S[0]
def replay(js,full=False):
	if full:ctx(full=True)
	B=load();A=[A for A in js if sample(A)!=B['rec'][A]];return dict(replayed=list(js),mismatches=A,generator_recomputed=full,ok=not A)
def alpha16(S=None):S=S or load();F=[sum(A[5])/1e2 for A in S['rec']];B=[len(A[5])for A in S['rec']];A=len(F);D=sum(F)/sum(B);J=sum(B)/A;K=math.sqrt(sum((A-D*B)**2 for(A,B)in zip(F,B))/(A-1)/A)/J;G=sum(log2(A)/(A-1)for A in M._primes_upto(1<<16));O=.08;C=[[2.**(O*A/1e2)for A in A[5]]for A in S['rec']];E=sum(B);H=sum(map(sum,C))/E;L=[(sum(A)**2-sum(A*A for A in A))/H**2 for A in C];I=[len(A)*(len(A)-1)for A in C];N=sum(L)/sum(I);P=math.sqrt(sum((A-N*B)**2 for(A,B)in zip(L,I))/(A-1)/A)/(sum(I)/A);Q=((E*H)**2-sum(sum(A)**2 for A in C))/H**2/(E*E-sum(len(A)**2 for A in C));R=sum(A[4][0]for A in S['rec']);T=sum(len(A[5])*(len(A[5])-1)//2 for A in S['rec']);return dict(kept_pairs=T,kept_pairs_sharing_a_prime_above_2_16=R,max_shared_rough_bits=max(A[4][1]for A in S['rec']),kept=sum(B),lattices=A,mean_kept=round(J,3),mean_s16=round(D,4),se_lattice=round(K,4),e_uniform=round(G,4),alpha_emp=round(G-D,4),alpha_emp_plus_3se=round(G-D+3*K,4),pair_moment_k008_within_lattice=round(N,4),pair_moment_k008_within_se=round(P,4),pair_moment_k008_between_lattices=round(Q,4))
def rough_med_audit(js,max_p=1<<20):
	E=max_p;H=bytearray(b'\x01')*(E+1);H[:2]=b'\x00\x00'
	for B in range(2,int(E**.5)+1):
		if H[B]:H[B*B::B]=b'\x00'*((E-B*B)//B+1)
	K=[A for A in range((1<<16)+1,E+1)if H[A]];A=list(K)
	while len(A)>1:A=[A[B]*A[B+1]if B+1<len(A)else A[B]for B in range(0,len(A),2)]
	W=A[0];Q=sum(log2(A)/(A-1.)for A in K);R=1.
	for L in K:R*=(1.-1./L)/(1.-L**.08/L)
	S=ctx();F=[];N=[];O=0;T=.0
	for X in js:
		Y=lattice(S,X);U=set();C=[]
		for G in cands(Y):
			if max(abs(A)for A in G[:6]).bit_length()>=128 or max(abs(A)for A in G[6:]).bit_length()>=160:continue
			I=ng(S,G)
			if I==0:continue
			Z=I not in U;U.add(I)
			if not Z or not M.is_ideal_primitive_oh(G[:6],G[6:]):continue
			J=I;D=math.gcd(W,J);V=1
			while D>1:
				while J%D==0:J//=D;V*=D
				D=math.gcd(D,J)
			C.append(log2(V))
		O+=len(C);T+=sum(C);F.append(sum(C)/len(C));N.append(sum(2.**(.08*A)for A in C)/len(C))
	P=T/O;a=math.sqrt(sum((A-P)**2 for A in F)/(len(F)-1)/len(F))if len(F)>1 else .0;return dict(lattices=list(js),kept_norms=O,band_bits=[16,int(round(log2(E)))],mean=round(P,4),se=round(a,4),e_unif=round(Q,4),alpha_med=round(Q-P,4),moment_008_ratio=round(sum(N)/len(N)/R,4),ok=True)
def joint_ea_mini_audit(u_val=-3,n_lat=60):
	l=n_lat;O=u_val;import functools as m;from math import log1p as n,floor as Z,ceil,sqrt;A,P,a,G,J,b,I=(14.,2e1,72,4,96,16,64.)if O==-3 else(16.,21.,72,4,96,16,8e1)
	def K(a,b):
		A=[.0]*(len(a)+len(b)-1)
		for(C,B)in enumerate(a):
			if B:
				for(D,E)in enumerate(b):A[C+D]+=B*E
		return A
	C=(P-A)/J;c=[n(C/(A+B*C))for B in range(J)];A4=[c[A]if A<b else .0 for A in range(J)];A5=[c[A]if b<=A<a else .0 for A in range(J)];A6=[c[A]if A>=a else .0 for A in range(J)];d=[[1.]];Q=[[1.]]
	for H in range(1,G+1):d.append([A/H for A in K(d[-1],A4)]);Q.append([A/H for A in K(Q[-1],A5)])
	o={}
	for B in range(G+2):
		e=Z((I+1e-09-B*A)/C-B)
		if e<0:continue
		p=Q[B][:e+1]if B<=G else[];L=K(Q[B-1],A6)[:e+1]if B>=1 else[]
		for R in range(G-B+2):
			D=list(p)if p and R<=G-B else[]
			if L and R<=G-B+1:
				if not D:D=list(L)
				else:A7=max(len(D),len(L));D=[(D[A]if A<len(D)else .0)+(L[A]if A<len(L)else .0)for A in range(A7)]
			if not D or not any(D):continue
			f=K(d[R],D);S=o.setdefault(R+B,[])
			if len(S)<len(f):S.extend([.0]*(len(f)-len(S)))
			for(T,A8)in enumerate(f):S[T]+=A8
	E=A+b*C;q=A+a*C
	@m.lru_cache(None)
	def A9(m):
		D=.0
		for(B,E)in o.items():
			F=min(len(E)-1,Z((m-B*C-B*A)/C))+1
			if F<=0:continue
			D+=sum(E[D]*(M.rho((m-B*A-D*C+DESC.ALPHA_G)/A)if m-B*A-D*C+DESC.ALPHA_G>0 else 1.)for D in range(F))
		return D
	r=96;U=(I-E)/r;AA=[n(U/(E+A*U))for A in range(r)];g=[1.];s=[]
	for H in range(1,int(I//E)+1):g=[A/H for A in K(g,AA)];s.append([(A,B)for(A,B)in enumerate(g)if B and H*E+A*U<=I+1e-09])
	@m.lru_cache(None)
	def AB(m):
		A=M.rho(m/E)if m>0 else 1.
		for(B,F)in enumerate(s,1):
			C=.0
			for(G,H)in F:D=max(m-B*E-(G+B)*U,.0)/E;C+=H*(1. if D<=1. else M.rho(D))
			A+=C
		return min(1.,A)
	V=1<<int(P);W=bytearray(b'\x01')*(V+1);W[:2]=b'\x00\x00'
	for F in range(2,int(V**.5)+1):
		if W[F]:W[F*F::F]=b'\x00'*((V-F*F)//F+1)
	t=[A for A in range(2,V+1)if W[A]];u=[]
	for T in range(0,len(t),256):
		v=t[T:T+256];w=1
		for F in v:w*=F
		u.append((w,v))
	x=ctx(O);AC=random.Random(20261123+O);AD,AE,AF,AG,AH=2.**E,2.**A,2.**q,2.**P,2.**I
	def y(nums):
		D=1
		for B in nums:D*=B
		E=[]
		for(P,Q)in u:
			F=math.gcd(D,P)
			if F>1:E+=[A for A in Q if F%A==0]
		H=I=0;J=K=.0
		for B in nums:
			L=log2(B);J+=AB(Z(L));K+=A9(ceil(L));M=C=B;N=O=0
			for A in E:
				if C%A==0:
					while C%A==0:
						C//=A
						if A<=AD:M//=A
						if A>AE:
							if A<=AF:N+=1
							elif A<=AG:O+=1
			if M<=AH:
				H+=1
				if C==1 and N<=G and O<=1:I+=1
		return H,J,I,K
	h=[];i=[];j=[];k=[];z=A0=A1=0
	for H in range(l):
		AI=lattice_e(x,H+1);A2=set();X=[]
		for N in cands(AI):
			if max(abs(A)for A in N[:6]).bit_length()>=128 or max(abs(A)for A in N[6:]).bit_length()>=160:continue
			Y=ng(x,N)
			if Y==0:continue
			AJ=Y not in A2;A2.add(Y)
			if AJ and M.is_ideal_primitive_oh(N[:6],N[6:]):X.append(Y)
		A1+=len(X);AK,AL,AM,AN=y(X);h.append(AK);i.append(AL);j.append(AM);k.append(AN);AO,_,AP,_=y([AC.randrange(max(2,int(A*.99)),int(A*1.01)+2)for A in X]);z+=AO;A0+=AP
	A3=sum(h)/sum(i);AQ=sum(j)/sum(k);return dict(u=O,lattices=l,kept=A1,x_ea=round(E,2),tau=I,LB=A,L_mid=round(q,2),L_top=P,surv_obs=sum(h),surv_ctrl=z,surv_mod=round(sum(i),1),surv_vs_mod=round(A3,4),use_obs=sum(j),use_ctrl=A0,use_mod=round(sum(k),1),use_vs_mod=round(AQ,4),ok=A3<=1.)
AW='phi_dag_one'
AS='w_dag_one'
AR='c_det_one_wc'
AP='univ_boot_wc_ir_int'
AO='univ_pair_setup_p59_ir_int'
AV='n_zero'
AN='two_stage_measured_ir_log2'
AU='special_q'
AM='measured_limb_ir_log2'
AL='s_max_det'
AT='t4_geom'
Bz='c_vals'
By='phi_vals'
Bx='w_vals'
Bw='boot_dp'
AK='OMEGA_BITS'
AJ='FIX_LLL'
AI='PAIR_SETUP'
AH='SIEVE_HIT'
AD='phi_dag_one_wc'
AC='w_dag_one_wc'
AB='c_tail_z6'
AA='obs_within_uniform'
A9='exact_div_errors'
A5='inductive_sr_inv_ok'
A4='monotone_d_ok'
A3='two_stage_ok'
BC='primes16'
A2=divmod
n='little'
u='fano_v_over_m'
l='swaps'
Ai='boot'
Ah='inf'
j=complex
o='descent'
AQ='cells'
m='B'
h='mean'
A8='rec'
y='plan59'
q=list
k=bool
d=enumerate
a=sorted
Z=abs
Y='f'
X=zip
W=pow
U=all
R='g'
Q=True
P=int
O=dict
M=tuple
N=float
J=False
K=min
I=None
H=round
G=sum
E=max
D=len
B=.0
A=1.
C=range
import functools as V,math as S,stor_bound as AX
from math import exp,factorial,floor,ceil,log,log1p,log2,sqrt
from fractions import Fraction as F
import c6_adaptive as L
psurv_iv_tab=L.psurv_iv_tab
psurv_iv=L.psurv_iv
psurv512=L.psurv512
def rand_below(getrandbits,n):
	B=getrandbits;C=n.bit_length();A=B(C)
	while A>=n:A=B(C)
	return A
AY=log(2.)
A6=2e4
T=A/1024.
Bn=16
MONT5=286.
TORSION=12
WHEELS=30,210,2310,30030,510510
ECM_T2_PAD=.05
SEPARATION=.6
MR_ROUNDS=110
NEWTON_PP_IR=2.**29
B_MIN_ECM=32
b=2
D1=.9
D2=2.
SURV_DM=1e1
Af=A
B9=.5
Bo=.001
Bp=112.
Bq=44
Br=6.*log2(0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47)
LOG2U6=6.*log2(0x44e992b44a6909f1)
QBINS=[(1/6,-.1,-2.1),(1/6,9.5,1.2),(1/6,15.85,3.75),(1/6,22.2,6.35),(1/6,29.,9.),(1/6,50.5,18.5)]
QBINS_LO=[(-95.,-55.),(-1e1,-55.),(B,-55.),(8.,-55.),(14.,-55.),(2e1,-55.)]
BA=2.
DAG_SURV_X=A
Bs=3e1
D3=2.**27
D4=2.**26
D5=2.**70
BOOT_KEEP_IR=8.
D6=64.
DF=48
ALPHA_F=-.995
ALPHA_G=-.946
BB=96.
Bt=.6869
c=96
A7=.999
BOUND_CELL=24.
D7=2.*log(59./8.)
SIEVE_HIT=2e1
SCAN=2.
ZERO=.25
PAIR_SETUP=2.**26.2
Bu=2.**33
NORMS=2.**20
TDIV=2.**20
IDEAL=2.**26
PRIME=2.**23
p=r=s=J
t=v=w=Q
TDIV20_PRIMES=82025.
FIX_LLL=2.**40
OMEGA_BITS=E(1,B_MIN_ECM//2)
_ECM_OMEGA={}
e,f,g,z,i,A0,A1=SIEVE_HIT,SCAN,ZERO,PAIR_SETUP,IDEAL,FIX_LLL,OMEGA_BITS
DESCBK_TYPED=O(SIEVE_HIT=e,SCAN=f,ZERO=g,PAIR_SETUP=z,IDEAL=i,FIX_LLL=A0,OMEGA_BITS=A1)
DESCBK_NEW=O(SIEVE_HIT=e,SCAN=f,ZERO=g,PAIR_SETUP=65888196.,IDEAL=i,FIX_LLL=960960331872.,OMEGA_BITS=20)
if p:SIEVE_HIT=DESCBK_NEW[AH]
if r:SCAN=DESCBK_NEW['SCAN'];ZERO=DESCBK_NEW['ZERO']
if t:PAIR_SETUP=DESCBK_NEW[AI]
if s:IDEAL=DESCBK_NEW['IDEAL']
if v:FIX_LLL=DESCBK_NEW[AJ]
if w:OMEGA_BITS=DESCBK_NEW[AK]
D8=4
D9=328,380
KNAP_ON=Q
def DA(l):B=l*AY;return l-log2(B)+log2(A+A/B+2./(B*B))
def mul_L(L):return 48. if L<=2 else MONT5 if L<=5 else 22.*L*L+32.*L+21.
def Bv(nbits):A=E(b,ceil(nbits/64.));return mul_L(A),A
def lehmer_gcd_ir(nbits):A=nbits;B=E(b,ceil(A/64.));return 2.*A*(23.*B+23.)+4.*B+1e1
DB=lehmer_gcd_ir
def AZ(x):return 1.26*x/log(x)if x>=17 else 7.
def Aa(x):return AZ(x)+2.*sqrt(x)
@V.lru_cache(I)
def _stage_ops(b1,b2,wheel):
	A=wheel;E=ceil(1.5*b1)*11;B=A
	for C in(2,3,5,7,11,13,17):
		if A%C==0:B=B//C*(C-1)
	F=B//2;D=ceil(b2/A)+1;G=6*(A//2+D)+3*F*D if b2>b1 else 0;return 1.1*(E+G),E,G,F,D
def _conv(a,b):
	I,J=D(a),D(b);K=[B]*(I+J-1);E=0
	while E<I and a[E]==B:E+=1
	if E==I:return K
	L=I-1
	while L>E and a[L]==B:L-=1
	F=0
	while F<J and b[F]==B:F+=1
	if F==J:return K
	M=J-1
	while M>F and b[M]==B:M-=1
	N=a[E:L+1];O=b[F:M+1];A,H=D(N),D(O)
	if A>H:N,O=O,N;A,H=H,A
	R=N[::-1];T=S.sumprod;U=E+F
	for G in C(A+H-1):P=0 if G<H else G-H+1;Q=G if G<A-1 else A-1;K[U+G]=T(R[A-1-Q:A-P],O[G-Q:G-P+1])
	return K
Ab=M(a({P(H(2.**(aa/2.)))for aa in C(12,74)}))
AE={}
_ECM_BEST_C={}
_SQ_NORMS_CACHE={}
_SQ_BASIS_CACHE={}
@V.lru_cache(I)
def Ac(b1,b2,steps=256):
	D=steps;A,E=log2(b1),log2(b2)
	if E<=A:return A,()
	B=(E-A)/D;return A,M((A+B*C,A+B*(C+1),log((A+B*(C+1))/(A+B*C)))for C in C(D))
def x(fbits,b1,b2,rt,hi=J,steps=256):
	N=fbits;C=hi;O=H(N,6),b1,b2,C;E=AE.get(O)
	if E is not I:return E
	P=N-log2(TORSION);F,Q=Ac(b1,b2,steps);J=D(rt);S=(lambda u:A if u<=A else rt[floor(u/T)]if floor(u/T)+1<J else A/factorial(floor(u)))if C else lambda u:A if u<=A else rt[ceil(u/T)]if ceil(u/T)<J else B;G=S(P/F)
	if Q:
		U=A/(F*T);V=floor if C else ceil
		for(W,X,L)in Q:
			M=P-(X if C else W)
			if M<=F:G+=L
			else:
				R=V(M*U)
				if R<J:G+=rt[R]*L
				elif C:G+=A/factorial(floor(M/F))*L
	E=K(G,A-1e-15);AE[O]=E;return E
def Ad(fbits,b1,rt):I=fbits-log2(TORSION);C=I/log2(b1);G=D(rt);F=E(C-2.,A);J=A if F<=A else rt[floor(F/T)]if floor(F/T)+1<G else A/factorial(floor(F));K=A if C<=A else rt[ceil(C/T)]if ceil(C/T)<G else B;H=sqrt(b1);return E(B,A-J/E(K,1e-300)*2./(H*log(H)))
def modinv_ir(nbits):A=nbits;B=E(1,ceil(A/64.));return 2.*A*(48.*B+41.)+1e1*B+1e1
AF={}
@V.lru_cache(I)
def branch_caps_all(nbits,b_min=B_MIN_ECM,l_min=b):
	a=l_min;T=P(ceil(nbits))-1;H=E(1,P(b_min)-1);b=H,a;J=AF.get(b)
	if J is I or D(J)<=T:
		G=E(T,2048)
		def f(b):A=E(a,b+63>>6);return mul_L(A),lehmer_gcd_ir(b),modinv_ir(b),2*b+20+63>>6,A
		U=[f(aa)for aa in C(G+2)];M=[B]*(G+1);N=[B]*(G+1);O=[B]*(G+1);Q=[0]*(G+1);R=[0]*(G+1)
		for A in C(G+1):
			if A<H:M[A],N[A],O[A],Q[A],R[A]=U[E(1,A+1)];continue
			c=A+1;V,W,X,Y,Z=U[c];d=A-H
			if d>=H:
				e=[H]
				for g in C(1,(c>>5)+2):
					S=g<<5
					for K in(S-11,S-10,S-1,S):
						if H<K<=d:e.append(K)
				for K in e:
					h,i,j,k,l=U[K+1];L=A-K;F=h+M[L]
					if F>V:V=F
					F=i+N[L]
					if F>W:W=F
					F=j+O[L]
					if F>X:X=F
					F=k+Q[L]
					if F>Y:Y=F
					F=l+R[L]
					if F>Z:Z=F
			M[A],N[A],O[A],Q[A],R[A]=V,W,X,Y,Z
		J=[(M[ab],N[ab],O[ab],Q[ab],R[ab])for ab in C(G+1)];AF[b]=J
	return J[E(0,T)]
@V.lru_cache(I)
def AG(d,nb,mul_tree,s_sum_max,l_sum_max):
	D=l_sum_max;B=mul_tree;import model_bn254_v6 as J;L=N((d+1)*(d+1))*(B+8.*D+12.);A=2*(d+1)*ceil((2*nb+20)/64.);F=2*(d+1)*ceil((2*B_MIN_ECM+20)/64.);assert 1<<E(2,ceil(log2(A)))<1<<54 and d+2<=1<<20;G=[F,A]
	for H in C(2,54):
		for I in(1<<H,3<<H-1):
			if F<=I+1<=A:G.append(I+1)
	M=E(J.ntt3_F(A)/N(A)for A in G);O=2*(d+1)*s_sum_max*(M+2e1)+(2*d+1)*(2.*B+2e1*D);return K(L,O)
@V.lru_cache(I)
def _ptree_tree_ir(d,nb,mul_tree,s_sum_max,l_sum_max):
	D=l_sum_max;C=s_sum_max;A=mul_tree
	if d<=1:return B
	E=d//2;return _ptree_tree_ir(E,nb,A,C,D)+_ptree_tree_ir(d-E,nb,A,C,D)+AG(d-E,nb,A,C,D)
def _ecm_best(fbits,nbits,success,rt):
	M=H(N(fbits),6),P(ceil(nbits)),N(success),OMEGA_BITS;g=_ECM_BEST_C.get(M)
	if g is not I:return g
	Z=H(M[0]+ECM_T2_PAD,4);G=M[1];D,V,n,h,W=branch_caps_all(G,B_MIN_ECM,b);o=E(b,ceil(G/64.));p=mul_L(o);F=G//OMEGA_BITS;_ECM_OMEGA[M]=F;q=(2.*MR_ROUNDS*(G+F)+8./3.*F*G)*p;O=2*F*NEWTON_PP_IR+q;X=12.*W+12.*F;i=log1p(-M[2]);a=I;T=N(Ah)
	for A in Ab:
		r=Aa(A);Q=1.1*ceil(1.5*A)*(11.*D+1e1*X)+r*V
		if Q+O>=T:continue
		c=Ad(Z,A,rt)
		if c<=B:continue
		j=x(Z,A,A*3000,rt,J)*c*SEPARATION
		if j<=1e-12 or ceil(i/log1p(-j))*Q+O>=T:continue
		for k in(3000,1000,500,300,200,100,50,30,10,1):
			d=A*k;e=x(Z,A,d,rt,J)*c*SEPARATION
			if e<=1e-12:break
			Y=ceil(i/log1p(-e))
			if Y*Q+O>=T:break
			for L in WHEELS if k>1 else(0,):
				if L and L>A:break
				if L==0:R=Y*Q+O
				else:
					A0,s,s,S,C=_stage_ops(A,d,L);f=ceil(log2(E(2,S)));R=Y*(Q+1.15*((L//2+C)*(6.*D+6.*X)+S*C*(3.*D+X))+(C+2*f*F)*V)+O;U=K(S,C);t=ceil(E(S,C)/N(U));u=1+ceil(C/N(U));v=1.15*((L//2+C)*(6.*D+6.*X)+4.*(S+C)*D)+u*n+2.*f*F*V;l=_ptree_tree_ir(U,G,D,h,W);w=AG(U,G,D,h,W);y=1.15*(l+t*(3.*l+4.*w+U*(D+4.*W)));z=ceil(log2(E(2,C)));m=Y*(Q+v+y+(C+3.*F*S)*D+(1+(2*z+1+2*f)*F)*V)+O
					if m<R:R=m
				if R<T:T=R;a=R,A,d,L,e
	_ECM_BEST_C[M]=a;return a
def boot_gate_knapsack(rows,vals,LB,x_ea,nA,nB,kA,tau):
	V=vals;Q=kA;P=nA;B=rows;W=D(B);assert W and B[0][0]==LB and LB<=x_ea<=P<nB and U(B[A][1]==B[A+1][0]for A in C(W-1))and B[-1][1]>=nB and D(V)==W and U(A>=0 for A in V);c=[N(aa).as_integer_ratio()for aa in V];e=E(A.bit_length()-1 for(B,A)in c);R=[ab<<e-ac.bit_length()+1 for(ab,ac)in c];Y=[];f=-1
	for(A,H,J)in a((0 if A<x_ea else floor(A),-R[B],B)for(B,(A,C))in d(B)if A<P):
		if-H>f:Y.append((A,-H,J));f=-H
	L=floor(tau);K=[[0]*(L+1)]
	for k in C(Q):
		g=K[-1];Z=g[:]
		for(A,H,J)in Y:Z[A:]=[ad if ad>ae+H else ae+H for(ad,ae)in X(Z[A:],g)]
		K.append(Z)
	S=K[Q][L];M=I
	for(J,(h,j))in d(B):
		if j>P and h<nB:
			A=floor(E(h,P))
			if A<=L and K[Q][L-A]+R[J]>S:S=K[Q][L-A]+R[J];M=J
	T=L-(floor(E(B[M][0],P))if M is not I else 0);b=[]
	for i in C(Q,0,-1):
		for(A,H,J)in Y:
			if A<=T and K[i-1][T-A]+H==K[i][T]:b.append(J);T-=A;break
	assert G(R[A]for A in b)+(R[M]if M is not I else 0)==S;return F(S,1<<e),O(A=b,B=M)
def Ag(x):A=N(x);return S.nextafter(A,S.inf)if F(A)<x else A
def evaluate_descent(desc_data,rho_fn,fft_F_fn,gbatch_cost_fn,pbits_fn,log2_rels,peel_h,nf_main=405.,ng_main=469.,b3_skip=B):
	BD=pbits_fn;Aj=rho_fn;z=desc_data;J=59.;DC=6;h={P(aa):ab for(aa,ab)in z[AQ].items()};BE=z['omitted_law'];B_=z['base_prim_c6'];DD=z[y];C0=M(z[Ai]);DE=M(z['ea']);l=getattr(Aj,'__globals__',{}).get('rho_tab_1024',lambda:[Aj(aV*T)for aV in C(32770)])()
	def DG(u):
		if u<=A:return A
		C=ceil(u/T);return l[C]if C<D(l)else B
	def C1(u):
		if u<=A:return A
		B=floor(u/T)
		if B+1>=D(l):return A/factorial(floor(u))
		return l[B]
	def DH(rem,base):C=ceil(rem/(base*T));return A if C<=0 else l[C]if C<D(l)else B
	def C2(u,t):A=Aj(u);return Aj(u-t)/A if A>B else B
	def Ak(fbits,nbits,success=.999):return _ecm_best(fbits,nbits,success,l)[0]
	def BF(lo_bits,hi_bits,nbits,success=.999,step=6.):A=hi_bits;return Ak(A,nbits,success)if A>lo_bits+1e-09 else B
	@V.lru_cache(I)
	def DI(fbits,nbits,logS):
		D=nbits;F=E(1,ceil(D/64.));G=BD(0,fbits,A);M=E(0,K(ceil(log2((ceil(G/64.)+1)/F)),P(logS)));B,H=I,0
		for J in C(0,M+1):
			L=gbatch_cost_fn(G,leaf=F,L=J)
			if B is I or L<B:B,H=L,J
		return B+DB(D)+12.*Bv(D)[0],H
	def BG(fbits,nbits,Ns):B=E(Ns,A);C,D=DI(fbits,nbits,log2(B));return C*(A+2.**D/E(Ns,1e-300))
	def DJ(x):return G(2.**(x-6.+.53-A-1)*(fft_F_fn(2.**(A+1))+3.*2.**(A+1))for A in C(ceil(x-6.+.53)))
	@V.lru_cache(I)
	def C3(n,base,iA,mesh,kmax=Bn):
		G=mesh;F=base;D=(n-F)/G;K=[log1p(D/(F+ax*D))if ax<iA else B for ax in C(G)];L=[log1p(D/(F+ay*D))if ay>=iA else B for ay in C(G)];E=[[A]]
		for M in C(1,kmax+1):E.append([az/M for az in _conv(E[-1],K)])
		N=[[A],L];I={}
		for H in C(kmax+1):
			for J in C(2):I[H,J]=E[H]if J==0 else _conv(E[H],N[1])
		return I,D
	@V.lru_cache(I)
	def DK(m,base,n,iA,mesh,alpha,kmax=Bn):
		H=kmax;E=base;W,I=C3(n,E,iA,mesh,H);Q=[[B,B]for aA in C(H+1)];X=D(l)
		for((F,G),J)in W.items():
			L=F+G;R=K(D(J)-1,P(floor((m-L*I-L*E)/I)))if I>0 else D(J)-1
			if R<0:continue
			Y=(m-L*E+alpha)/(E*T);Z=I/(E*T);S=B
			for N in C(R+1):
				if J[N]:O=ceil(Y-N*Z);S+=J[N]*(A if O<=0 else l[O]if O<X else B)
			Q[F][G]=S
		U=[[B,B]for aB in C(H+1)]
		for G in C(2):
			V=B
			for F in C(H+1):V+=Q[F][G];U[F][G]=V
		return M(M(A)for A in U)
	def BH(m,base,n,iA,kA,kB,mesh,alpha):
		C=alpha;A=base;m=ceil(m)
		if n<=A:return DG((m+C)/A)
		D=DK(m,A,n,iA,mesh,C,kA);return D[kA][0]+(D[kA][1]if kB else B)
	def DL(m,base,n,iA,kA,kB,mesh,alpha):
		F=alpha;C=base;m=floor(m)
		if n<=C:return C1((m+F)/C)
		N,D=C3(n,C,iA,mesh,kA);G=B
		for((H,I),O)in N.items():
			if H>kA or I and not kB:continue
			E=H+I
			for(J,L)in d(O):
				if L and E*C+J*D<=m:M=(m-(E*C+J*D)-E*D+F)/C;G+=L*(A if M<=A else C1(M))
		return K(G,A)
	def Al(m,base,T):
		if T<=0:return B
		if T>=m:return A
		return K(A,2.*psurv512(floor(m),base,T,512))
	DM=E(B,BE[Y][3]+BE[R][3]-A);Am=DM*B_
	def C4(side,m,n,iA,kA,kB,alpha,slack):B=alpha;D=P(floor((slack-A-m/128.)*D8));F,H,I,L=BE[side];C=lambda xg:DL(E(m-xg,A),J,n,iA,kA,kB,c,B)if m-xg>A else A;K=G(H[E(0,D-D9[0]):])*C(128.)+G(B*C(64.*A)for(A,B)in d(I,2));return K/(BH(m,J,n,iA,kA,kB,c,B)*F)
	def DN(nf,ng,n,iA,kA,kB,slack):C=slack;D=C4(Y,ceil(nf),n,iA,kA,kB,ALPHA_F,C);F=C4(R,ceil(ng),n,iA,kA,kB,ALPHA_G,C);return E(B,A-D-F)*Am
	@V.lru_cache(I)
	def C5(u4,u_main,bskip=B):
		G=bskip;L=u4/4.;H=B
		for E in C(300,P(10*J)):
			I=E/1e1;D=I+.1;M=2.**(log2_rels-D)*DC*C2(u_main,D/J);F=exp(-M*(A-(peel_h if E>=544 else .4)))
			if G>B and E>=550:F=K(A,F+G)
			H+=C2(L,D/J)*log(D/I)*F
		return H
	DO,DP=nf_main/J,ng_main/J
	def BI(nf,ng):return(C5(ceil(4.*nf/J),DO,B)if nf else B)+(C5(ceil(4.*ng/J),DP,b3_skip)if ng else B)
	DQ=2.*2.**DA(J);DR=E(A[1]for A in QBINS);DS=E(A[2]for A in QBINS);An={};BJ={};BK={};Ao={}
	def A9(s,lim):return E([aE for((aC,aD),aE)in An.items()if aD==s and aC[0]<lim-1e-09]+[B])
	def AA(s,lim):return E([aH for((aF,aG),aH)in BJ.items()if aG==s and aF[0]<lim-1e-09]+[B])
	def AB(s,lim):return E([aK for((aI,aJ),aK)in BK.items()if aJ==s and aI[0]<lim-1e-09]+[B])
	for(m,Ap,e,AR,AS,W,Z,f,_)in DD:
		assert AR<=m and(AR<m or W==0 and Z==0),f"non-decreasing DAG edge at ({m},{e})";Q=N(AR);C6=log2(h[f]);Aq=(Ap+C6)/12.;C7=(m+C6)/12.;DT=(Q-J)/c if Q>J else B;n=J+AS*DT if Q>J else J;BL=W*n+Z*Q+BB;BM=ceil(W*(F(P(J))+(F(AS*(AR-P(J)),c)if Q>J else 0))+Z*AR);r=AC=C8=B
		for((BN,DU,DV),(DW,DX))in X(QBINS,QBINS_LO):DY=DF*Aq+DU-(m if e==Y else B);DZ=12.*Aq+LOG2U6+DV-(m if e==R else B);r+=BN*BH(DY,J,Q,AS,W,Z,c,ALPHA_F)*BH(DZ,J,Q,AS,W,Z,c,ALPHA_G);C9=DF*C7+DW-(Ap if e==Y else B);CA=12.*C7+LOG2U6+DX-(Ap if e==R else B);AC+=DAG_SURV_X*BN*Al(C9,J,BL)*Al(CA,J,BL);C8+=DAG_SURV_X*BN*Al(C9,J,BM)*Al(CA,J,BM)
		AD=DF*Aq+DR-(m if e==Y else B);AE=12.*Aq+LOG2U6+DS-(m if e==R else B);Da=A-(A-A7)*((AD+AE)/2e1+2.*(W+Z));Db=BI(AD,AE);AT=h[f]*r*Bt*Da*E(B,A-Db);CB=Ak(J,AD,A7);CC=Ak(J,AE,A7);CD=NORMS+TDIV+CB+CC+2.*(BF(J,n,BL,A7)+(W+Z)*PRIME)+IDEAL;CE=h[f]*(D7*SIEVE_HIT+SCAN+ZERO+BOUND_CELL)+DQ*PAIR_SETUP
		if h[f]*CD<CE+h[f]*(AC+r)*CD:BO='ecm';AT*=Am;Ar=h[f]
		else:BO='sieve';AT*=DN(AD,AE,Q,AS,W,Z,BB);Ar=h[f]*(AC+r)
		CF=NORMS+TDIV+BG(J,AD,Ar)+BG(J,AE,Ar)+IDEAL;Dc=2.*TDIV20_PRIMES*(ceil(E(AD,AE)/64.)*5.+1e1);CG=Dc+2.*(BF(J,n,BM,A7)+(W+Z)*PRIME)+CB+CC;CH=K(h[f],h[f]*(C8+r));Dd=(Bu+h[f]*CF+CH*CG if BO=='ecm'else Bu+CE+Ar*CF+CH*CG)+2**40;CI=E(Bo,(A+B9*AT/Af)**-Af)if AT>B else A;De=W*(AA(Y,n)+AA(R,n))+Z*(AA(Y,Q)+AA(R,Q));CJ=CI+De;CK=Dd+W*(A9(Y,n)+A9(R,n))+Z*(A9(Y,Q)+A9(R,Q));CL=A+W*(AB(Y,n)+AB(R,n))+Z*(AB(Y,Q)+AB(R,Q));As=m,Ap;An[As,e]=CK;BJ[As,e]=CJ;BK[As,e]=CL;Ao[As,e]=O(how=BO,R=f,mu=AT,p0=CI,phi_dag=CJ,c_det=CL,W_log2=log2(CK/A6))
	A0,CM,s,EC=C0;CN,AF=DE;CO=(A0-J)/c;AG=J+CM*CO;g=J+CN*CO
	@V.lru_cache(I)
	def Df(base,n,iA,kA,mesh,kx_i,tau_i):
		O=base;I=mesh;H=kA;K=(n-O)/I;P=[log1p(K/(O+aL*K))for aL in C(I)];W=[P[aM]if aM<kx_i else B for aM in C(I)];X=[P[aN]if kx_i<=aN<iA else B for aN in C(I)];Y=[P[aO]if aO>=iA else B for aO in C(I)];Q=[[A]];L=[[A]]
		for T in C(1,H+1):Q.append([aP/T for aP in _conv(Q[-1],W)]);L.append([aQ/T for aQ in _conv(L[-1],X)])
		U={}
		for F in C(H+2):
			R=floor((tau_i+1e-09-F*O)/K-F)
			if R<0:continue
			V=L[F][:R+1]if F<=H else[];J=_conv(L[F-1],Y)[:R+1]if F>=1 else[]
			for M in C(H-F+2):
				G=[]
				if V and M<=H-F:G=q(V)
				if J and M<=H-F+1:
					if not G:G=q(J)
					else:Z=E(D(G),D(J));G=[(G[aR]if aR<D(G)else B)+(J[aR]if aR<D(J)else B)for aR in C(Z)]
				if not G or not any(G):continue
				S=_conv(Q[M],G);N=U.setdefault(M+F,[])
				if D(N)<D(S):N.extend([B]*(D(S)-D(N)))
				for(a,b)in d(S):N[a]+=b
		return U,K
	@V.lru_cache(I)
	def A1(m):
		L,A=Df(J,N(A0),CM,s,c,CN,AF);F=B
		for(E,H)in L.items():
			I=K(D(H)-1,floor((m-E*A-E*J)/A))+1
			if I<=0:continue
			F+=G(H[B]*DH(m-E*J-B*A+ALPHA_G,J)for B in C(I))
		return F
	i=z[Bw];t=Br+LOG2U6;AU=ceil(i['ng_hi'])
	@V.lru_cache(I)
	def BP(m):return D1*(A-.001*(m/2e1+2.*(s+1)))*Bt*Am*.5*E(B,A-BI(B,m))
	@V.lru_cache(I)
	def At(m):return K(A,D2*psurv512(m-SURV_DM,g,AF,512))
	import model_bn254_v6 as Dg;Au=Bv(1975)[0];Dh=(1<<31)*Au;CP=lambda da,db:Dg.ntt3_F(63*(da+db+1))+(da+db+1)*4.*Au;Di=G((1<<31-A)*6.*CP(1<<A,1<<A)for A in C(32));Dj=G((1<<32-A)*G(2.*CP(2<<A,2<<A)for A in C(A))for A in C(1,32));Av=Dh+Di+Dj+(1<<32)*(1224500+2.*Au)+61.*(1<<32)*Au+1975.*IDEAL;Aw=N(1<<59);assert Av<Aw and log2(Av)<=58.324;AH=64;CQ=10611996;A2=CQ;CR=0
	for CS in C(2,AH+1):
		while A2%CS==0:A2//=CS
	for j in C(1,AH):
		BQ=1
		for BR in C(1,AH+1):BQ=BQ*(j*AH+BR)%CQ
		if S.gcd(A2,BQ%A2)>1:
			CR+=1
			for CT in C(j*AH+1,(j+1)*AH+1):
				while A2%CT==0:A2//=CT
	assert A2==1 and CR==2;o=lambda a,b,n:[G(a[B]*b[aU-B]for B in C(E(0,aU-D(b)+1),K(D(a),aU+1)))%n for aU in C(D(a)+D(b)-1)]
	def Dk(m,d,n):
		B=[1];A=1;H=m[::-1]
		while A<d:A<<=1;E=H[:A];assert D(E)+D(B)-1<=2*A+1;F=o(E,B,n)[:A];G=[((2 if aS==0 else 0)-F[aS])%n for aS in C(D(F))];assert D(B)+D(G)-1<=2*A+1;B=o(B,G,n)[:A]
		return B[:d]
	def Ax(a,b,n):
		E,A=D(a)-1,D(b)-1
		if E<A:return a[:]
		B=E-A+1;F=Dk(b,B,n);G=a[::-1][:B];assert D(G)+D(F)-1<=2*A+1;H=o(G,F,n)[:B][::-1];assert D(H)+D(b)-1<=2*A+1;I=o(H,b,n);return[(a[aT]-(I[aT]if aT<D(I)else 0))%n for aT in C(A)]
	AI=13199;BS=[1]
	for BR in C(1,9):BS=o(BS,[BR,1],AI)
	A3=[[[-ak*8%AI,1]for ak in C(8)]]
	while D(A3[-1])>1:A3.append([o(A3[-1][2*at],A3[-1][2*at+1],AI)for at in C(0,D(A3[-1])//2)])
	BT=[Ax(BS,A3[-1][0],AI)]
	for CU in reversed(A3[:-1]):BT=[Ax(BT[al//2],CU[al],AI)for al in C(D(CU))]
	assert[am[0]for am in BT]==[S.prod(an*8+B for B in C(1,9))%AI for an in C(8)];BU=BG(g,AU,2.**Bq);Ay=Ak(g,AU,.999)+BF(g,AG,AF,.999)+(s+1)*PRIME+Aw;BV=FIX_LLL+D3+144.*(NORMS+D4)+D6;b=[];AJ=[];A4=[];BW=[];u=[];BX=K(BP(A)for A in C(1,AU+1))
	for AV in i[A8]:Az=[A1(ceil(t+ao/1e2+.005))for ao in AV[3]];BY=G(Az);AC=G(At(floor(t+A/1e2-.005))for A in AV[3]);b.append(BY);BW.append((BY*BY-G(A*A for A in Az))/2.);u.append(AC);AJ.append(BV+D(Az)*(BU+IDEAL+BOOT_KEEP_IR)+AC*Ay);A4.append(D(Az))
	L=D(b);v=G(b)/L;p=G(AJ)/L;A_=G(BW)/L;CV=G(A4)/L;B0=sqrt(G((A-v)**2 for A in b)/(L-1));B1=sqrt(G((A-p)**2 for A in AJ)/(L-1));CW=B0/sqrt(L);CX=B1/sqrt(L);Dl=G(((A-v)/B0)**3 for A in b)/L;Dm=G(((A-v)/B0)**4 for A in b)/L;Dn=G(((A-p)/B1)**3 for A in AJ)/L;Do=G(((A-p)/B1)**4 for A in AJ)/L;B2=v-3.*CW;Dp=p+3.*CX;CY=CV-3.*sqrt(G((A-CV)**2 for A in A4)/(L-1)/L);CZ=sqrt(G((A-A_)**2 for A in BW)/(L-1)/L);Ca=2.**Bs*(A_+3.*CZ)/(B2/BA);w=[ceil(t+ad/1e2+.005)for ac in i[A8]for ad in ac[3]];Cb=[(log2(A1(ae))-log2(A1(ae+8)))/8. for ae in C(K(w)-80,E(w)+1,8)];Cc=K(Cb);Cd=E(Cb);B3=BZ=I
	for AW in(Cc,Cd):
		AX=A
		for AY in i[BC]:AX*=(A-A/AY)/(A-AY**(AW-A))
		Ce=[G(A1(ceil(t+A/1e2+.005))*2.**(AW*(B/1e2+ALPHA_G))for(A,B)in X(ap[3],ap[5]))/AX for ap in i[A8]];Ba=G(Ce)/G(b);Cf=sqrt(G((A-Ba*B)**2 for(A,B)in X(Ce,b))/(L-1)/L)/v
		if B3 is I or Ba-3.*Cf<B3-3.*BZ:B3,BZ=Ba,Cf
	Dq=U(At(A)>=At(A+1)for A in C(K(w)-2,E(w)+1))and U(A1(A)>=A1(A+1)for A in C(K(w),E(w)+1));AK=V.lru_cache(I)(lambda m:psurv512(m,g,AF,512));Cg=[floor(t+ag/1e2-.005)for af in i[A8]for ag in af[3]];Bb=[(log2(AK(ah))-log2(AK(ah+8)))/8. for ah in C(K(Cg)-90,E(Cg)+1,8)];Ch=K(Bb);B4=B
	for AW in(K(Bb),E(Bb)):
		AX=A
		for AY in i[BC]:AX*=(A-A/AY)/(A-AY**(AW-A))
		Ci=[G(At(floor(t+A/1e2-.005))*2.**(AW*B/1e2)for(A,B)in X(aq[3],aq[5]))/AX for aq in i[A8]];Cj=G(Ci)/G(u);Dr=sqrt(G((A-Cj*B)**2 for(A,B)in X(Ci,u))/(L-1)/L)/(G(u)/L);B4=E(B4,Cj+3.*Dr)
	Ds=2.**(A+1e1*Ch)/B4;AL=[B]*161;AL[0]=A;Bc={}
	for AZ in i[BC]:
		if AZ<=256:
			Aa=[B]*161;Bd=A;Ck=A/AZ;Dt=A-Ck;j=0;Du=log2(AZ)
			while Bd>1e-13:
				A5=P(H(j*Du))
				if A5>160:break
				Ax=Dt*Bd
				for Ab in C(161-A5):Aa[Ab+A5]+=AL[Ab]*Ax
				Bd*=Ck;j+=1
			AL=Aa
		else:B5=P(H(log2(AZ)));Bc[B5]=Bc.get(B5,B)+A/AZ
	for(B5,Cl)in a(Bc.items()):
		Aa=[B]*161;Be=exp(-Cl);j=0
		while Be>1e-14:
			A5=j*B5
			if A5>160:break
			for Ab in C(161-A5):Aa[Ab+A5]+=AL[Ab]*Be
			j+=1;Be*=Cl/j
		AL=Aa
	Bf={}
	for(Dv,r)in d(AL):
		if r>1e-08:Cm=Dv+1&~1;Bf[Cm]=Bf.get(Cm,B)+r
	Cn=q(Bf.items());Bg=V.lru_cache(I)(lambda m:A1(m-ALPHA_G));Dw=V.lru_cache(I)(lambda m:G(B*Bg(m-A)for(A,B)in Cn));Dx=V.lru_cache(I)(lambda m:G(B*AK(m-A)for(A,B)in Cn));Bh=[];Bi=[];B6=[]
	for AV in i[A8]:
		Co=Cp=o=B
		for(Cq,Dy)in X(AV[3],AV[5]):Bj=ceil(t+Cq/1e2+.005);B7=floor(t+Cq/1e2-.005);Cr=P(H(Dy/1e2))+1&~1;Co+=Bg(Bj-Cr)*(Bg(Bj)/Dw(Bj));Cp+=AK(B7-Cr)*(AK(B7)/Dx(B7));o+=AK(B7)
		Bh.append(Co);Bi.append(Cp);B6.append(o)
	Cs=G(Bh)/G(b);Dz=sqrt(G((A-Cs*B)**2 for(A,B)in X(Bh,b))/(L-1)/L)/(G(b)/L);Ct=G(Bi)/G(B6);D_=sqrt(G((A-Ct*B)**2 for(A,B)in X(Bi,B6))/(L-1)/L)/(G(B6)/L);Cu=[BV+ai*(BU+IDEAL+BOOT_KEEP_IR)+4.*aj*Ay for(ai,aj)in X(A4,u)];Cv=G(Cu)/L;Cw=Cv+3.*sqrt(G((A-Cv)**2 for A in Cu)/(L-1)/L);AM=BX*(B2/BA)*(A-Ca);ED=B2;E0=BI(B,AU);EE=G(u)/G(A4);AN=Dp;Ac=(2.**Bq+144.)/CY*AN;B8=2.*DJ(g)+D5;Bk=s*A9(R,AG)+A9(R,N(A0));Ad=s*AA(R,AG)+AA(R,N(A0));Bl=s*AB(R,AG)+AB(R,N(A0));Ae=a(A[0]for A in An if A[1]==R);Bm=O(on=k(KNAP_ON),rows=Ae,w_vals=[An[au,R]for au in Ae],phi_vals=[BJ[av,R]for av in Ae],c_vals=[BK[aw,R]for aw in Ae],params=O(LB=J,x_ea=g,nA=AG,nB=N(A0),kA=s,tau=AF));AO,AP,Cx=Bk,Ad,Bl
	if KNAP_ON:Cy,Cz,C_=[boot_gate_knapsack(Ae,Bm[ar],J,g,AG,N(A0),s,AF)for ar in(Bx,By,Bz)];AO=K(Ag(Cy[0]),Bk);AP=K(Ag(Cz[0]),Ad);Cx=K(Ag(C_[0]),Bl);Bm.update(w_cfg=Cy[1],phi_cfg=Cz[1],c_cfg=C_[1])
	x=A/(A-AP)if AP<A else N(Ah);E1=A/(A-Ad)if Ad<A else N(Ah);E2=x/AM;E3=AN/AM+Ac+AO;E4=B8+x*E3;D0=2.*E4/A6;E5=BX*(B2/4.-2.**Bs*(A_+3.*CZ));E6=log2(2.*(B8+x*(AN/E5+Ac+AO))/A6);E7=log2(2.*(B8+x*(Cw/AM+Ac*Cw/AN+AO))/A6);E8=log2(D0);E9=D0/2.**Bp;EA=log2(2.*x*AO/A6);EB=log2(2.*(B8+x*(AN/AM+Ac))/A6);return O(ok=AP<=.6 and x<5. and Av<Aw,mode='pure_k1_boole_wald_iid_dp_boot',nodes=D(Ao),lb=J,nb_r=Af,nb_alpha=B9,p0_floor=Bo,base_prim=B_,counted_frac=Am,two_log=E8,desc_cap_log2=Bp,markov_cap_fail=E9,two_dag=EA,two_boot_exp=EB,boot_split_ps_ir_log2=H(log2(Av),3),boot_split_det_charge_log2=H(log2(Aw),3),boot_cfg=q(C0),boot_dp_lattices=L,boot_dp_mean_k=G(A4)/L,boot_dp_E_mean=v,boot_dp_E_se=CW,boot_dp_E_cv=H(B0/v,4),boot_dp_E_min_ratio=H(K(b)/v,4),boot_dp_E_skew=H(Dl,3),boot_dp_E_kurt=H(Dm,2),boot_dp_W_mean_log2=log2(p),boot_dp_W_se_rel=CX/p,boot_dp_W_cv=H(B1/p,4),boot_dp_W_max_ratio=H(E(AJ)/p,4),boot_dp_W_skew=H(Dn,3),boot_dp_W_kurt=H(Do,2),boot_dp_K_min=K(A4),boot_dp_K_lo=CY,boot_dp_Q_mean=A_,boot_dp_pair_loss=Ca,boot_dp_strip_ir=BU,boot_dp_cst_min=BX,boot_dp_cst_used_min_max=[K(BP(A)for A in set(w)),E(BP(A)for A in set(w))],boot_dp_slope_lo=Cc,boot_dp_slope_hi=Cd,boot_dp_first_moment_factor=B3,boot_dp_first_moment_factor_se=BZ,boot_dp_surv_monotone=Dq,boot_dp_ecm_ir=Ay,boot_dp_fix_ir=BV,boot_dp_surv_per_lattice=G(u)/L,boot_charge=BA,boot_two_log_if_charge4=E6,boot_two_log_if_survivors_x4=E7,boot_dp_overshoot_share=Ac/(AN/AM),boot_dp_ecm_share_of_W=G(A*Ay for A in u)/L/p,boot_dp_surv_slope_lo=Ch,boot_dp_surv_smallprime_factor=B4,boot_dp_surv_margin=Ds,boot_dp_first_moment_exact_nl=H(Cs,4),boot_dp_first_moment_exact_nl_se=H(Dz,4),boot_dp_surv_exact_nl=H(Ct,4),boot_dp_surv_exact_nl_se=H(D_,4),boot_dp_ng_hi=AU,boot_x_ea=H(g,4),boot_root_retained_bytes_log2=H(log2(BD(0,g,A)/8.),4),boot_root_bytes_log2=H(log2(8.*BD(0,g,A)/8.),4),boot_prob=AM,boot_miss=E0,phi_dag_one=AP,p_boot_tree=AP,e_runs_per_log=x,e_trials_per_log=E2,c_det_one=Cx,w_dag_one=AO,w_dag_one_wc=Bk,phi_dag_one_wc=Ad,e_runs_per_log_wc=E1,c_det_one_wc=Bl,boot_knap=Bm,p_node_abort=B,p_repeat_loss=B,p_descent_fail=B,repeat_prob=B,mu_min=K(A['mu']for A in Ao.values()),p0_max=E(A['p0']for A in Ao.values()))
def bareiss_lll_12d(B,G6,sa_num=1,sa_den=1,track_stats=J,dnum=9,dden=10,stage2=I,ret_dl=J):
	p=S.sumprod;d=12;M=[1]+[0]*d;V=[[0]*d for aa in C(d)];AG=[[p(ad,ab[:6])for ad in G6]for ab in B];AH=[[p(ae,ac[6:])for ae in G6]for ac in B]
	for e in C(d):
		AI=AG[e];AJ=AH[e];X=V[e]
		for L in C(e+1):
			A8=B[L];AK=V[L];r=p(AI,A8[:6])*sa_num+p(AJ,A8[6:])*sa_den
			for F in C(L):r=(M[F+1]*r-X[F]*AK[F])//M[F]
			if L<e:X[L]=r
			else:M[e+1]=r
	A9=((dnum,dden),stage2)if stage2 is not I else((dnum,dden),)
	if not track_stats:
		for(z,A0)in A9:
			F=1
			while F<d:
				N=V[F];A=B[F]
				for P in C(F-1,-1,-1):
					T=N[P]
					if not T:continue
					R=M[P+1]
					if T>0:
						if T<<1>R:K=((T<<1)+R)//(R<<1)
						else:continue
					elif-T<<1>R:K=-(((-T<<1)+R)//(R<<1))
					else:continue
					D=B[P];i=V[P]
					if K==1:
						A[0]-=D[0];A[1]-=D[1];A[2]-=D[2];A[3]-=D[3];A[4]-=D[4];A[5]-=D[5];A[6]-=D[6];A[7]-=D[7];A[8]-=D[8];A[9]-=D[9];A[10]-=D[10];A[11]-=D[11]
						for L in C(P):N[L]-=i[L]
					elif K==-1:
						A[0]+=D[0];A[1]+=D[1];A[2]+=D[2];A[3]+=D[3];A[4]+=D[4];A[5]+=D[5];A[6]+=D[6];A[7]+=D[7];A[8]+=D[8];A[9]+=D[9];A[10]+=D[10];A[11]+=D[11]
						for L in C(P):N[L]+=i[L]
					else:
						A[0]-=K*D[0];A[1]-=K*D[1];A[2]-=K*D[2];A[3]-=K*D[3];A[4]-=K*D[4];A[5]-=K*D[5];A[6]-=K*D[6];A[7]-=K*D[7];A[8]-=K*D[8];A[9]-=K*D[9];A[10]-=K*D[10];A[11]-=K*D[11]
						for L in C(P):N[L]-=K*i[L]
					N[P]=T-K*R
				Y=M[F];f=M[F-1];g=M[F+1];c=N[F-1];s=f*g+c*c
				if A0*s<z*Y*Y:
					B[F-1],B[F]=A,B[F-1];j=V[F-1]
					if F>1:j[:F-1],N[:F-1]=N[:F-1],j[:F-1]
					A1=F-1
					for e in C(F+1,d):X=V[e];k=X[A1];l=X[F];X[A1]=(k*c+l*f)//Y;X[F]=(k*g-l*c)//Y
					M[F]=s//Y;F=A1 if F>1 else 1
				else:F+=1
		return(B,M,V)if ret_dl else B
	AA=q(M);A3=E(A.bit_length()for A in M[1:]);a=E(Z(V[A][B]).bit_length()for A in C(d)for B in C(A));b=1<<a;AL=G(A.bit_length()for A in M[1:12]);m=AB=A4=0;AC=A5=Q
	for(z,A0)in A9:
		F=1
		while F<d:
			N=V[F];A=B[F]
			for P in C(F-1,-1,-1):
				T=N[P]
				if not T:continue
				R=M[P+1]
				if T>0:
					if T<<1>R:K=((T<<1)+R)//(R<<1)
					else:continue
				elif-T<<1>R:K=-(((-T<<1)+R)//(R<<1))
				else:continue
				D=B[P];i=V[P]
				if K==1:A[0]-=D[0];A[1]-=D[1];A[2]-=D[2];A[3]-=D[3];A[4]-=D[4];A[5]-=D[5];A[6]-=D[6];A[7]-=D[7];A[8]-=D[8];A[9]-=D[9];A[10]-=D[10];A[11]-=D[11]
				elif K==-1:A[0]+=D[0];A[1]+=D[1];A[2]+=D[2];A[3]+=D[3];A[4]+=D[4];A[5]+=D[5];A[6]+=D[6];A[7]+=D[7];A[8]+=D[8];A[9]+=D[9];A[10]+=D[10];A[11]+=D[11]
				else:A[0]-=K*D[0];A[1]-=K*D[1];A[2]-=K*D[2];A[3]-=K*D[3];A[4]-=K*D[4];A[5]-=K*D[5];A[6]-=K*D[6];A[7]-=K*D[7];A[8]-=K*D[8];A[9]-=K*D[9];A[10]-=K*D[10];A[11]-=K*D[11]
				t=T-K*R;N[P]=t
				if t>=b or t<=-b:a=Z(t).bit_length();b=1<<a
				for L in C(P):
					u=N[L]-K*i[L];N[L]=u
					if u>=b or u<=-b:a=Z(u).bit_length();b=1<<a
				AD=(Z(K).bit_length()+63)//64;A6=(R.bit_length()+63)//64;m+=12*A6*A6+(P+1)*5*AD*A6+240*AD
			for L in C(F):
				if 2*Z(N[L])>M[L+1]:A5=J
			Y,f,g,c=M[F],M[F-1],M[F+1],N[F-1];h=(Y.bit_length()+63)//64;m+=15*h*h
			if A0*(f*g+c*c)<z*Y*Y:
				B[F-1],B[F]=A,B[F-1];j=V[F-1]
				if F>1:j[:F-1],N[:F-1]=N[:F-1],j[:F-1]
				for e in C(F+1,d):
					X=V[e];k,l=X[F-1],X[F];AM=k*c+l*f;AN=k*g-l*c;v,AO=A2(AM,Y);w,AP=A2(AN,Y)
					if AO or AP:A4+=1
					X[F-1]=v;X[F]=w;m+=44*h*h
					if v>=b or v<=-b:a=Z(v).bit_length();b=1<<a
					if w>=b or w<=-b:a=Z(w).bit_length();b=1<<a
				s=f*g+c*c;AE,AQ=A2(s,Y)
				if AQ:A4+=1
				M[F]=AE;m+=12*h*h
				if AE>AA[F]:AC=J
				AB+=1;F=F-1 if F>1 else 1
			else:F+=1
	if not U(2*Z(V[A][B])<=M[B+1]for A in C(1,d)for B in C(A)):A5=J
	x=E(A.bit_length()for A in AA[1:]);W=(x+63)//64;n=(ceil(1.5*x+9)+63)//64;y=(x+9+63)//64;A7=(ceil(.5*x+9)+63)//64;o=ceil(AL/log2(1e1/9.));AR=11*(12*n*W+60*y*A7)+275*y*W;AS=12*n*W+60*y*A7+22*W;AT=15*W*W;AU=12*W*W+440*n*W;AF=(o+11)*AR+o*AS+(2*o+11)*AT+o*AU+6188*W*W;return B,O(swaps=AB,s_max_det=o,monotone_d_ok=AC,inductive_sr_inv_ok=A5,exact_div_errors=A4,max_d_bits=A3,max_d_limbs=(A3+63)//64,max_lam_bits=a,max_lam_limbs=(a+63)//64,uniform_limbs=[W,n,y,A7],obs_within_uniform=(a+63)//64<=n and(A3+63)//64<=W,measured_limb_ir_log2=H(log2(E(1,m)),2),wc_det_ir_int=AF,wc_det_ir_log2=ceil(log2(AF)*1e2)/1e2)
def verify_bareiss_910_cert(B,G6,sa_num=1,sa_den=1):
	I=S.sumprod;D=12;A=[1]+[0]*D;F=[[0]*D for aa in C(D)];O=[[I(ad,ab[:6])for ad in G6]for ab in B];P=[[I(ae,ac[6:])for ae in G6]for ac in B]
	for E in C(D):
		Q=O[E];R=P[E];L=F[E]
		for G in C(E+1):
			M=B[G];T=F[G];K=I(Q,M[:6])*sa_num+I(R,M[6:])*sa_den
			for H in C(G):
				N=A[H+1]*K-L[H]*T[H]
				if N%A[H]!=0:return J
				K=N//A[H]
			if G<E:L[G]=K
			else:A[E+1]=K
	V=U(2*Z(F[B][D])<=A[D+1]for B in C(1,D)for D in C(B));W=U(10*(A[B-1]*A[B+1]+F[B][B-1]*F[B][B-1])>=9*A[B]*A[B]for B in C(1,D));return V and W
def verify_live_bareiss_12d_lll(seed=20261120):
	import bisect,json,random;from pathlib import Path;import model_bn254_v6 as A;AB=A.load_data();s=AB[AT];Y=A.build_exact_int_geom(s)[0];D=A.P_BN254;f=W(3,-1,D);AC=A.load_witnesses();Ak=[[1,0,1,-1,-3,-3],[0,1,-1,0,-1,-1],[0,0,-1,1,1,3],[0,0,-1,-1,0,2],[0,0,1,0,1,-1],[0,0,0,1,1,2]];t=0xa8ccdee3ccdbd36d;AD=[A.b62d(aa)for aa in AC[Ai][0][0].split()];AE=[[D if al==6+ab else 0 for al in C(12)]for ab in C(6)]
	for L in C(6):AF=[1 if am==L else 0 for am in C(6)];u=[G(AF[B]*A.ZK_NUM[B][an]for B in C(6))*f%D for an in C(6)];Al=[3*ao%D for ao in A._pmul_h(AD,u,D)];Am=[G(Ak[ap][A]*Al[A]for A in C(6))*f%D for ap in C(6)];AE.append(AF+[aJ-D if aJ>D//2 else aJ for aJ in Am])
	def v(B_in,sn,sd,lim):G,A=bareiss_lll_12d(B_in,Y,sn,sd,Q,11,40);I,B=bareiss_lll_12d(G,Y,sn,sd,Q,9,10);D=A[l]+B[l];C,J,E,K=A['uniform_limbs'];F=A['wc_det_ir_int']+11*(12*J*C+60*E*K+275*E*C)+6188*C*C;B.update(stage1_swaps=A[l],stage2_swaps=B[l],two_stage_total_swaps=D,two_stage_s_max_d0=A[AL],two_stage_joint_wc_ir_int=F,two_stage_joint_wc_ir_log2=ceil(log2(F)*1e2)/1e2,two_stage_measured_ir_log2=H(log2(2.**A[AM]+2.**B[AM]),2));B[A3]=A[A4]and A[A5]and A[A9]==0 and A[AA]and D<=A[AL]and B['two_stage_joint_wc_ir_log2']<lim;return I,B
	V,I=v(AE,t,1,4e1);w=verify_bareiss_910_cert(V,Y,t,1)
	for AG in V:
		u=[G(AG[B]*A.ZK_NUM[B][aq]for B in C(6))*f%D for aq in C(6)];An=[G(AG[6+B]*A.ZK_NUM[B][ar]for B in C(6))*f%D for ar in C(6)]
		if A._pmul_h(AD,u,D)!=An:w=J
	AH=[[j(at,au)for(at,au)in ac]for ac in s['EMB_ZK']];Ao=[j(ad,ae)for(ad,ae)in s['EMB_W']];AI=M(M(A.real for A in A)for A in AH);AJ=M(M(A.imag for A in A)for A in AH);BI=log2(AB[o][AQ]['4']);Ap=random.Random(seed);X=[factorial(af)for af in C(13)];AK=[]
	for a in C(2):
		for b in C(2):
			for c in C(5):
				if(Aq:=16-16*a-9*b-4*c)<0:continue
				for x in C(K(Aq,12-a-b-c)+1):
					if(AO:=x+c+b+a)==0:continue
					AK.append((X[12]//(X[12-AO]*X[x]*X[c]*X[b]*X[a])*(1<<AO),(x,c,b,a)))
	Ar=[];AP=0
	for(As,g)in AK:AP+=As;Ar.append(AP)
	d=[];y=q(V)
	for z in C(12):
		for L in C(z+1,12):y.append([V[z][aK]+V[L][aK]for aK in C(12)]);y.append([V[z][aL]-V[L][aL]for aL in C(12)])
	for AR in y:
		AS,AV=AR[:6],AR[6:];AW=E(0,E(E(Z(A)for A in AS),E(Z(A)for A in AV),1).bit_length()-40);AX=N(1<<AW);AY=[av/AX for av in AS];AZ=[aw/AX for aw in AV];Aa=B
		for L in C(3):Ab=j(S.sumprod(AI[L],AY),S.sumprod(AJ[L],AY));Ac=j(S.sumprod(AI[L],AZ),S.sumprod(AJ[L],AZ));h=A.U_BN254*Ao[L]*Ab*Ab-Ac*Ac;Aa+=log2(h.real*h.real+h.imag*h.imag)+4.*AW
		d.append(Aa-Br-LOG2U6)
	d.sort();I['live_boot0_144_offsets_min_med_max']=[H(d[0],2),H(d[72],2),H(d[-1],2)];I['boot0_exact_match']=w;i={}
	for(At,Au)in((18,'sq162_g'),(23,'sq213_f')):
		g,Av,Aw,Ax=AC[AU][At];R,Ay,Az=A.b62d(Av),A.b62d(Aw),A.b62d(Ax);A_=W(3,-1,R);A0=[1]
		for g in C(5):A0.append(A0[-1]*Ay%R)
		Ad=[G(A.ZK_NUM[ax][B]*A0[B]for B in C(6))*A_%R for ax in C(6)];A1=Ad+[-Az*aF%R for aF in Ad];Ae=W(A1[0],-1,R);B0=[[R if aM==0 else 0 for aM in C(12)]]+[[-A1[aG]*Ae%R-(R if-A1[aG]*Ae%R>R//2 else 0)if aN==0 else 1 if aN==aG else 0 for aN in C(12)]for aG in C(1,12)];g,i[Au]=v(B0,1,1,34.)
	e=(1<<220)-1;B1=[Ap.randrange(-e//2,e//2)for ag in C(11)];B2=[[e if aH==0 else 0 for aH in C(12)]]+[[B1[ay]if aI==0 else 1 if aI==ay+1 else 0 for aI in C(12)]for ay in C(11)];g,i['sq220_wc']=v(B2,1,1,34.)
	def Af(mat,k):
		B=[[F(mat[aO][aP])for aP in C(k)]for aO in C(k)];D=F(1)
		for A in C(k):
			E=next(C for C in C(A,k)if B[C][A]!=0)
			if E!=A:B[A],B[E]=B[E],B[A];D=-D
			G=B[A][A];D*=G
			for H in C(A+1,k):
				J=B[H][A]/G
				for I in C(A,k):B[H][I]-=J*B[A][I]
		return D
	B3=[[Y[ah][az]-(50000 if ah==az else 0)for az in C(6)]for ah in C(6)];Ag=U(Af(B3,A)>0 for A in C(1,7));T=[P(Af(Y,ai))for ai in C(1,7)];B4=T+[T[5]*T[aA-7]for aA in C(7,13)]
	def k(d_arr,extra=0,den=10,num=9):I=d_arr;B=E(A.bit_length()for A in I);L=G(A.bit_length()for A in I[:11]);A,D,F,H=(B+63)//64,(ceil(1.5*B+9)+63)//64,(B+9+63)//64,(ceil(.5*B+9)+63)//64;C=ceil(L/log2(N(den)/N(num)));J=11*(12*D*A+60*F*H)+275*F*A;M=12*D*A+60*F*H+22*A;O=(C+11)*J+C*M+(2*C+11)*(15*A*A)+C*(12*A*A+440*D*A)+6188*A*A;K=O+11*J+11*(15*A*A)+6188*A*A+extra;return K,ceil(log2(K)*1e2)/1e2,B,C,[A,D,F,H]
	B5=[D**(2*aB)*T[aB-1]for aB in C(1,7)]+[D**12*T[5]*t**(aC-6)*T[aC-7]for aC in C(7,13)];B6=[e*e*B4[aj-1]for aj in C(1,13)];B7=[((1<<59)-1)**2]*12;B8=[1764**aD*T[aD-1]for aD in C(1,7)]+[1764**6*T[5]*1000**(aE-6)*T[aE-7]for aE in C(7,13)];Ah=(1<<62)-1;B9=[Ah*Ah*B8[ak-1]for ak in C(1,13)];A2,A6,BA,BB,BC=k(B5,0);A7,A8,BJ,BD,BE=k(B6,0);m,n,BK,BL,BM=k(B7,13828);p,r,BN,BF,BG=k(B9,0,100,99);Aj=Ag and BA==I['max_d_bits']and A2<=2.**A6 and A6<=39.81 and A2<=FIX_LLL and A7<=2.**A8 and A8<=32.89 and A7<2**36 and m<=2.**n and n<=25.98 and m<=PAIR_SETUP and p<=2.**r and r<=34.45 and p<2**34.5;BH=Aj and I[A4]and I[A5]and I[A9]==0 and I[AA]and I[A3]and I['wc_det_ir_log2']<4e1 and I[AN]<32. and w and U(A[A4]and A[A5]and A[A9]==0 and A[AA]and A[A3]and A[AN]<28. for A in i.values());return O(ok=BH,zero_float_lll=Q,unified_delta='9/10',input_uniform_wc=Q,input_uniform_wc_proved=Aj,g6_lambda_min_gt_50k=Ag,univ_boot_wc_ir_int=A2,univ_boot_wc_ir_log2=A6,univ_boot_s_max=BB,univ_boot_limbs=BC,univ_sq220_wc_ir_int=A7,univ_sq220_wc_ir_log2=A8,univ_sq220_s_max=BD,univ_sq220_limbs=BE,univ_pair_setup_p59_ir_int=m,univ_pair_setup_p59_ir_log2=n,pair_setup_budget_log2=H(log2(PAIR_SETUP),3),pair_setup_p59_proved=m<=PAIR_SETUP and n<=25.98,univ_lat_b2_p61_ir_int=p,univ_lat_b2_p61_ir_log2=r,univ_lat_b2_p61_s_max=BF,univ_lat_b2_p61_limbs=BG,lat_b2_p61_proved=p<2**34.5 and r<=34.45,inductive_sr_inv_ok=Q,all_step_lam_tracked=Q,boot_p6=I,special_q_live=i)
def verify_live_mont5_and_suyama12_ecm(seed=20261120):
	import random;V=random.Random(seed);I=(1<<64)-1;b=S.gcd
	def B7(a_int,b_int,N_int,n0p):
		F=N_int;N=[a_int>>64*af&I for af in C(5)];O=[b_int>>64*ag&I for ag in C(5)];L=[F>>64*ah&I for ah in C(5)];A=[0]*7
		for P in C(5):
			Q=N[P];H=0
			for B in C(5):D=Q*O[B]+A[B]+H;A[B]=D&I;H=D>>64
			E=A[5]+H;A[5]=E&I;A[6]=E>>64;M=A[0]*n0p&I;J=A[0]+M*L[0]>>64
			for B in C(1,5):D=M*L[B]+A[B]+J;A[B-1]=D&I;J=D>>64
			E=A[5]+J;A[4]=E&I;A[5]=A[6]+(E>>64);A[6]=0
		K=G(A[B]<<64*B for B in C(6))
		if K>=F:K-=F
		return K
	def B8(a_int,b_int,N_int,n0p):
		F=b_int;E=a_int;B=N_int;G=[E&I,E>>64&I];H=[F&I,F>>64&I];J=[B&I,B>>64&I];A=[0,0,0]
		for K in C(2):L=G[K]*H[0]+A[0];M=L&I;R=L>>64;N=G[K]*H[1]+A[1]+R;S=N&I;T=A[2]+(N>>64);O=M*n0p&I;U=M+O*J[0]>>64;P=O*J[1]+S+U;A[0]=P&I;Q=T+(P>>64);A[1]=Q&I;A[2]=Q>>64
		D=A[0]|A[1]<<64|A[2]<<128;return D-B if D>=B else D
	AI=AJ=AK=Q
	for e in C(32):
		F=V.getrandbits(318)|1<<318|1;B9=-W(F,-1,1<<64)&I;z=V.randrange(F);A0=V.randrange(F)
		if B7(z,A0,F,B9)!=z*A0*W(1<<320,-1,F)%F:AI=J
		h=(1<<128)-1 if e==0 else(1<<127)+1 if e==1 else V.getrandbits(126)|1<<126|1;BA=-W(h,-1,1<<64)&I;Aa=V.randrange(h);n=V.randrange(h)
		if B8(Aa,n,h,BA)!=Aa*n*W(1<<128,-1,h)%h:AJ=J
		f,M=z|1,A0;AL=0
		while M>0:
			AL+=1
			if M&1==0:M>>=1
			else:
				if f<M:f,M=M,f
				if f==M:break
				f,M=M,f-M>>1
		if f!=b(z|1,A0)or AL>640 or AL*138.+3e1>lehmer_gcd_ir(320):AK=J
	Ab=279;Ac=47<=mul_L(2)and Ab<=mul_L(5)and U(22*A*A+32*A+21<=mul_L(A)for A in C(6,101))
	def A1(n):
		if n<2 or n%2==0:return n==2
		B,D=n-1,0
		while B%2==0:B>>=1;D+=1
		for E in(2,3,5,7,11,13,17):
			if E>=n:continue
			A=W(E,B,n)
			if A in(1,n-1):continue
			for F in C(D-1):
				A=A*A%n
				if A==n-1:break
			else:return J
		return Q
	def A2(X,Z,a24,N):B=(X+Z)%N;C=(X-Z)%N;D=B*B%N;A=C*C%N;E=(D-A)%N;return D*A%N,E*(A+a24*E)%N
	def A3(X1,Z1,X2,Z2,Xd,Zd,N):A=(X1-Z1)*(X2+Z2)%N;B=(X1+Z1)*(X2-Z2)%N;C=(A+B)%N;D=(A-B)%N;return Zd*C*C%N,Xd*D*D%N
	def A4(k,X,Z,a24,N):
		E=a24;A,B=X,Z;C,D=A2(X,Z,E,N)
		for F in bin(k)[3:]:
			if F=='1':A,B=A3(A,B,C,D,X,Z,N);C,D=A2(C,D,E,N)
			else:C,D=A3(A,B,C,D,X,Z,N);A,B=A2(A,B,E,N)
		return A,B
	BC=[aa for aa in C(2,180)if A1(aa)];BD=[ab**P(log(180)/log(ab))for ab in BC];A5=30;BE=[ac for ac in C(1,A5//2)if b(ac,A5)==1]
	def Ad(arr,modN):
		C=arr;B=modN;A=b(C[-1],B)
		if A==1:return 1
		I=A if 1<A<B else B;E,F=0,D(C)-1
		while E+1<F:
			G=(E+F)//2;H=b(C[G],I)
			if 1<H<I:return H
			if H==1:E=G
			else:F=G
		return A if 1<A<B else 1
	Ae=p=q=AM=AN=AO=Af=r=s=0
	for e in C(8):
		L=V.randrange(1<<16,(1<<17)-1)|1
		while not A1(L):L+=2
		X=V.randrange(1<<16,(1<<17)-1)|1
		while not A1(X)or X==L:X+=2
		AP=V.randrange(1<<221,(1<<222)-1)|1
		while not A1(AP):AP+=2
		F=L*X*AP
		for AR in C(6,14):
			Ae+=1;i=(AR*AR-5)%F;A6=4*AR%F;BF=(A6-i)%F;Ag=16*i*i*i*A6%F;AS=b(Ag,F)
			if AS>1:
				if AS in(L,X):
					p+=1;q+=1;AO+=1
					if AS==L:r+=1;s+=1
				continue
			t=W(BF,3,F)*(3*i+A6)*W(Ag,-1,F)%F;j,d=W(i,3,F),W(A6,3,F);k=[d]
			for BG in BD:j,d=A4(BG,j,d,t,F);k.append(k[-1]*d%F)
			BH=k[-1]%L==0;BI=k[-1]%X==0;BJ,BK=A2(j,d,t,F);A8={1:(j,d),2:(BJ,BK)}
			for u in C(3,A5//2):BL,BM=A8[u-1];BN,BO=A8[u-2];A8[u]=A3(BL,BM,j,d,BN,BO,F)
			AT,AU=A4(A5,j,d,t,F);A9,AA=A4(6,AT,AU,t,F);Aj,Ak=A4(5,AT,AU,t,F);l=1;AV=[1]
			for e in C(6,55):
				for u in BE:BP,BQ=A8[u];l=l*(A9*BQ-AA*BP)%F
				AV.append(l);BR,BS=A3(A9,AA,AT,AU,Aj,Ak,F);Aj,Ak,A9,AA=A9,AA,BR,BS
			g=BH or l%L==0;AB=BI or l%X==0
			if not(g or AB):continue
			p+=1
			if g:r+=1
			if g and AB:AM+=1
			AC=Ad(k,F)
			if AC in(L,X):
				q+=1;AO+=1
				if g and(AC==L or b(k[-1],F//AC)==L or b(l,F//AC)==L):s+=1
				if g and AB:AN+=1
			else:
				AW=Ad(AV,F)
				if AW in(L,X):
					q+=1;Af+=1
					if g and(AW==L or b(AV[-1],F//AW)==L):s+=1
					if g and AB:AN+=1
	Al=q/N(E(1,p));Am=s/N(E(1,r));import model_bn254_v6 as An;m=An.rho_tab_1024();v=An.load_data()[o];Z=set();BT=E(A[1]for A in QBINS);BU=E(A[2]for A in QBINS);BV={P(ad):ae for(ad,ae)in v[AQ].items()}
	for(Ao,BW,Ap,BX,BY,BZ,Ba,Bb,e)in v[y]:
		AX=N(BX);AD=59.+BY*(AX-59.)/c if AX>59. else 59.;Bc=log2(BV[Bb]);Aq=(BW+Bc)/12.;Bd=ceil(DF*Aq+BT-(Ao if Ap==Y else B));Be=ceil(12.*Aq+LOG2U6+BU-(Ao if Ap==R else B));Ar=ceil(BZ*AD+Ba*AX);Bf=ceil(Ar+BB);Z.add((59.,Bd));Z.add((59.,Be))
		if AD>59.+1e-09:Z.add((H(AD,4),Ar));Z.add((H(AD,4),Bf))
	Bg,Bh,Bu,e=v[Ai];Bi,Bj=v['ea'];As=(Bg-59.)/c;Bk=H(59.+Bi*As,4);Bl=H(59.+Bh*As,4);Z.add((Bk,ceil(v[Bw]['ng_hi'])));Z.add((Bl,Bj));AE=N(Ah);AY=0;AZ=Q;At=log2(TORSION);AF=.05;Au=D(m)
	def Bm(u):return A if u<=A else m[floor(u/T)]if floor(u/T)+1<Au else A/factorial(floor(u))
	def Bn(u):return A if u<=A else m[ceil(u/T)]if ceil(u/T)<Au else B
	for(AG,Bo)in a(Z):
		e,w,n,Av,Aw=_ecm_best(AG,Bo,A7,m)
		if Av and Av>w:AZ=J
		Ax=sqrt(w);Bp=2./(Ax*log(Ax));Ay=A/log2(w);Bv=A/Aw;Az=E(1,P(ceil((AG-16.)/AF)));A_=SEPARATION/Aw
		for B0 in C(0,Az,20):
			B1=K(Az,B0+20);Bq=x(K(AG,16.+B1*AF),w,n,m,J)
			for Br in C(B0,B1):
				B2=16.+Br*AF;B3=K(AG,B2+AF);B4=(B2-At)*Ay;Bs=(B3-At)*Ay;B5=E(B,A-Bm(B4-2. if B4>3. else A)/E(Bn(Bs),1e-300)*Bp);AH=Bq*B5*A_
				if AH<1.05:AH=x(B3,w,n,m,J)*B5*A_
				if AH<AE:AE=AH
				AY+=1
	B6=AZ and AE>=A and D(Z)>=40 and AY>=50000;Bt=AI and AJ and AK and Ac and p>=45 and AM>=15 and r>=40 and Al>=.85 and Am>=.85 and B6;return O(ok=Bt,mont5_32_exact_ok=AI,mont2_32_exact_ok=AJ,lehmer32_exact_ok=AK,cios5_exact_ir=Ab,mont5_budget_ir=MONT5,all_L_5_100_below_budget=Ac,suyama12_curves=Ae,suyama12_smooth_hits=p,suyama12_isolated_single_prime=q,suyama12_multi_factor_hits=AM,suyama12_multi_factor_split=AN,suyama12_p1_hits=r,suyama12_p1_isolated=s,suyama12_p1_cond_iso_rate=H(Am,4),suyama12_s1_iso=AO,suyama12_s2_iso=Af,suyama12_per_curve_cond_iso_rate=H(Al,4),ecm_1rung_t2_ok=B6,ecm_1rung_targets=D(Z),ecm_1rung_intervals_checked=AY,ecm_1rung_min_iv_ratio=H(AE,6),ecm_1rung_w_le_b1=AZ)
def verify_live_z6_la_kernels(seed=20261120):
	import random,model_bn254_v6 as L;B=L.R_BN254;A=W(L.P_BN254,10,B);e=(1<<317)//B;M=(1<<64)-1;F=random.Random(seed);assert e<1<<64 and(A*A-A+1)%B==0
	def f(A):
		F=A<0;E=-A if F else A;G=[E>>64*an&M for an in C(5)];I=(G[4]<<3|G[3]>>61)&M;H=I*e>>64;D=E-H*B
		if D>=B:D-=B
		if D>=B:D-=B
		if D>=B or D<0:return-1,-1
		if F and D:D=B-D
		return D,E//B-H
	g=Q;h=0
	for v in[0,1,-1,B-1,-(B-1),B,-B,(1<<59)*(B-1),-((1<<59)*(B-1)),(1<<315)-1,-((1<<315)-1)]+[F.randrange(-(1<<315),1<<315)for ac in C(50)]:
		A7,a=f(v)
		if A7!=v%B or a<0 or a>2:g=J
		if a>h:h=a
	b=Q
	for T in C(8):
		R=[0,0,0];N=0
		for V in(1,-1):
			for S in(0,1,2):I=F.randrange(B);R[S]+=V*I;i=1 if S==0 else A if S==1 else(A-1)%B;N=(N+V*i*I)%B
		for(j,k)in[(1,1),(2,-3),(-2,5),(-255,255)]+[(F.randrange(-255,256),F.randrange(-255,256))for aj in C(8)]:I=F.randrange(B);A8,A9=1 if j>=0 else-1,1 if k>=0 else-1;R[0]+=A8*(Z(j)*I);R[1]+=A9*(Z(k)*I);N=(N+(j+k*A)*I)%B
		P=[f(ad)[0]for ad in R];l=(P[0]-P[2]+A*(P[1]+P[2]))%B
		if l!=N:b=J
	AA=5*(B-1)<1<<256;Y=128;m=[[F.randrange(B)for ae in C(16)]for aa in C(Y)];m[0]=[B-1]*16;c=[[[m[ab][af]>>64*ak&M for ak in C(4)]for af in C(16)]for ab in C(Y)];w=(1<<43)-128;x=U(w+A<<9<=(1<<64)-1 and w+A>>32>0 for A in C(Y))
	def X(cols):
		Q=cols;R=[[0,0,0,0,0]for ao in C(16)];S=0;E=0;V=D(Q);L=0
		while L+15<=V:
			N=Q[L:L+15];E+=23
			for I in C(16):
				H=R[I];E+=3
				for W in(0,5,10):
					J=q(c[N[W]][I]);E+=2
					for T in C(1,5):
						U=c[N[W+T]][I];E+=6;A=0
						for B in C(4):F=J[B]+U[B]+A;J[B]=F&M;A=F>>64
						if A:S+=1
					A=0
					for B in C(4):F=H[B]+J[B]+A;H[B]=F&M;A=F>>64
					H[4]=H[4]+A&M;E+=5
				E+=3
			E+=2;L+=15
		O=V-L
		if O>0:
			N=Q[L:];E+=(O+1)//2+O
			for I in C(16):
				H=R[I];E+=3;P=0
				while P<O:
					X=K(5,O-P);J=q(c[N[P]][I]);E+=2
					for T in C(1,X):
						U=c[N[P+T]][I];E+=6;A=0
						for B in C(4):F=J[B]+U[B]+A;J[B]=F&M;A=F>>64
						if A:S+=1
					A=0
					for B in C(4):F=H[B]+J[B]+A;H[B]=F&M;A=F>>64
					H[4]=H[4]+A&M;E+=5;P+=X
				E+=3
		return[G(R[ap][A]<<64*A for A in C(5))for ap in C(16)],S,E
	y=0;n=0
	for z in[[0]*200]+[[F.randrange(Y)for al in C(ag)]for ag in C(0,35)]:
		AC,o,T=X(z);n+=o
		if AC!=[G(m[B][am]for B in z)for am in C(16)]:y+=1
	A0=X([0]*15)[2];p=A0/24e1;AD=U(X([0]*(15+A))[2]-X([0]*A)[2]==A0 for A in C(15));AE=E(X([0]*A)[2]-L.C_UNIT_5*16.*A for A in C(15));r=6.*AE/16.;s=323.+r;A1=0;AF=0;A2=F.getrandbits
	for T in C(64):
		A3=[[]for ah in C(6)]
		for T in C(100):A3[rand_below(A2,6)].append(rand_below(A2,Y))
		for A4 in A3:T,o,AG=X(A4);n+=o;A1+=AG;AF+=D(A4)*16
	A5=A1/1024.;A6=H(p,4);t=AA and x and y==0 and n==0 and AD and p<=L.C_UNIT_5 and s<=L.OP[AB]and A5<=1e2*L.C_UNIT_5+r;AH=t;u=Q
	for T in C(12):
		R=[0,0,0];N=0
		for d in(1,2,3,4,5,17):
			for S in(0,1,2):
				V=F.choice((1,-1));I=F.randrange(B);i=1 if S==0 else A if S==1 else(A-1)%B
				if d<=4:
					for T in C(d):R[S]+=V*I
				else:R[S]+=V*(d*I)
				N=(N+V*d*i*I)%B
		P=[f(ai)[0]for ai in R];l=(P[0]-P[2]+A*(P[1]+P[2]))%B
		if l!=N:u=J
	return O(ok=g and b and t and u and s<=L.OP[AB],mu_r_hex=hex(e),barrett_5limb_exact_ok=g,barrett_max_q_err=h,z6_3bucket_2coord_row_ok=b,unit4_5entry_subblock_max=5,unit4_5limb_carry_ok=AH,unit4_8vec_ok=t,unit4_8vec_u64_idx_43b_ok=x,unit4_8vec_10blk_ir=H(p,4),unit4_8vec_realistic_ir=A6,unit4_8vec_w100_row_ir=H(A5,2),unit4_8vec_6bkt_tail_wc_ir=H(r,2),row_tail_total_wc_ir=H(s,2),unit4_v4_expand_ok=u,unit4_reg_exact_ir=A6,unit4_reg_ir=L.C_UNIT_5,gen1_reg_2coord_ir=42.,gen1_sign_split_42op_ok=b,row_tail_exact_ir=323,c_tail_budget_ir=L.OP[AB],zero_tabm_plain_sparse=Q)
def verify_deg8_dispersion_and_prod_yield(seed=20261120):
	import json;from pathlib import Path;import model_bn254_v6 as L,verify_bn254_v6 as AW;A7=L.load_data();BA=A7[AT];BB=L.build_exact_int_geom(BA)[0];AX=L.ZK_NUM;v=L.W_COEFFS;A8,BC=AW._eval_zk_mod,AW._eval_poly_mod;AY=L.load_witnesses();p=AY['deg8_campaign'];c=p['tower'];V,N,AZ,BD,Aa,BE,BF,BG=c['u_t'],c['p_t'],c['r_t'],c['lb'],c['nA'],c['nB'],c['kA'],c['kB'];i,Ab,w,Ac=1<<BD,1<<Aa,1<<BE,1<<2*Aa;q=L._isp;BH=L._ppow_h([0,1,0,0,0,0],N,N);Ad=L._ppow_h(BH,N,N);Ae=L._ppow_h(Ad,N,N);BI=L._ppow_h(Ae,N**3,N);r=[-((V+1)*ab+(2 if aa==0 else 0))%N for(aa,ab)in d(v)];A9=r;AA=[1,0,0,0,0,0]
	for x in C(6):AA=L._pmul_h(AA,A9,N);A9=L._ppow_h(A9,N,N)
	z=L._pmul_h(r,r,N);BJ=L._pmul_h(z,r,N);BK=L._pmul_h(z,z,N);BL,BM,BN,BO,BP=[L.F_REC_COEFFS[ac]for ac in(8,6,4,2,0)];AB=[ad%N for ad in BL]
	for(BQ,BR)in((BM,r),(BN,z),(BO,BJ),(BP,BK)):BS=L._pmul_h([ay%N for ay in BQ],BR,N);AB=[(ai+aj)%N for(ai,aj)in X(AB,BS)]
	Ag=q(N)and q(AZ)and N==36*V**4+36*V**3+24*V**2+6*V+1 and AZ==36*V**4+36*V**3+18*V**2+6*V+1 and BI==[0,1,0,0,0,0]and Ad!=[0,1,0,0,0,0]and Ae!=[0,1,0,0,0,0]and W(AA[0],(N-1)//2,N)==N-1 and AB==[0]*6;Ah=[F(ae,3)for ae in L.SIGMA_NUM];AC=[[F(1)]+[F(0)]*5,Ah]
	for x in C(4):AC.append(L._pmul_h_q(AC[-1],Ah))
	Ai=M(M(P(3*AC[B][A])for B in C(6))for A in C(6));T=L._mh_int;Aj=S.sumprod
	def Ak(p):B=M(Aj(A,p)for A in Ai);E=M(Aj(A,B)for A in Ai);C=T(T(p,B),E);A=-C[2]//27;D=C[0]//27-2*A;return D*D+2*A*A
	def BT(pa,pb):A=M(A//3 for A in T(pa,pa));D=M(A//3 for A in T(pb,pb));E=T(v,A);return A,D,Ak(M((V+1)*E[B]+2*A[B]+D[B]for B in C(6)))//729
	def BU(pa2,pb2):E=pb2;D=pa2;A=M(A//3 for A in T(D,D));B=M(A//3 for A in T(E,E));F=M(A//3 for A in T(D,E));K=M(A//3 for A in T(A,A));G=M(A//3 for A in T(B,B));H=M(A//3 for A in T(A,F));I=M(A//3 for A in T(A,B));J=M(A//3 for A in T(F,B));L=M(K[A]+8*H[A]+30*I[A]-4*J[A]-5*G[A]for A in C(6));N=M(5*H[A]+6*I[A]-7*J[A]-2*G[A]for A in C(6));O=T(v,N);return Ak(M(L[A]+O[A]for A in C(6)))//729
	def AD(xs):
		A=xs
		while D(A)>1:A=[A[aL]*A[aL+1]if aL+1<D(A)else A[aL]for aL in C(0,D(A),2)]
		return A[0]if A else 1
	AE=L._primes_upto(Ab);BV=AD([ak for ak in AE if ak<=i]);BW,BX=AD([az for az in AE if i<az<=8192]),AD([aA for aA in AE if 8192<aA<=Ab])
	def Al(n):return L._brent(n)
	def Am(n):
		A=S.gcd(n,BV)
		while A>1:n//=A;A=S.gcd(n,A)
		if n==1:return Q
		F=0
		for H in(BW,BX):
			A=H
			while n>1:
				A=S.gcd(n,A)
				if A==1:break
				n//=A;B=[A]
				while B:
					E=B.pop()
					if E<1<<20 or q(E):
						F+=1
						if F>BF:return J
					else:D=Al(E);B.append(D);B.append(E//D)
		if n==1:return Q
		if n<Ac:return n<=w
		if n>w*w or q(n):return J
		B=[n];G=0
		while B:
			C=B.pop()
			if C==1:continue
			if C<Ac or q(C):
				if C>w:return J
				G+=1
				if G>BG:return J
				continue
			D=Al(C)
			if D==C:return J
			B.append(D);B.append(C//D)
		return Q
	BY=[((aB,1),)for aB in C(12)]+[((aC,1),(aD,aE))for aC in C(12)for aD in C(aC+1,12)for aE in(1,-1)]+[((am,1),(an,ap),(ao,aq))for am in C(12)for an in C(am+1,12)for ao in C(an+1,12)for ap in(1,-1)for aq in(1,-1)]
	def A0(side,q,s,rho,first=J):
		b=first;c=[A8([1 if aT==aM else 0 for aT in C(6)],s,q)for aM in C(6)];g=W(c[0],-1,q);N=[[0]*12 for aN in C(12)];N[0][0]=q
		for K in C(1,6):N[K][K]=1;N[K][0]=-c[K]*g%q
		for K in C(6):N[6+K][6+K]=1;N[6+K][0]=rho*c[K]*g%q
		bareiss_lll_12d(N,BB,1,1,J);Y=[([G(aO[B]*AX[B][aU]for B in C(6))for aU in C(6)],[G(aO[6+B]*AX[B][aV]for B in C(6))for aV in C(6)])for aO in N];S=[aP[0]for aP in Y]+[[-aW for aW in aQ[0]]for aQ in Y];T=[aR[1]for aR in Y]+[[-aX for aX in aS[1]]for aS in Y];d=0;h=set();e=I
		for H in BY:
			if D(H)==1:O=H[0][0];Z=M(S[O]);a=M(T[O])
			elif D(H)==2:O=H[0][0];X=H[1][0]+(0 if H[1][1]==1 else 12);A,B,E,F=S[O],S[X],T[O],T[X];Z=A[0]+B[0],A[1]+B[1],A[2]+B[2],A[3]+B[3],A[4]+B[4],A[5]+B[5];a=E[0]+F[0],E[1]+F[1],E[2]+F[2],E[3]+F[3],E[4]+F[4],E[5]+F[5]
			else:O=H[0][0];X=H[1][0]+(0 if H[1][1]==1 else 12);i=H[2][0]+(0 if H[2][1]==1 else 12);A,B,U,E,F,V=S[O],S[X],S[i],T[O],T[X],T[i];Z=A[0]+B[0]+U[0],A[1]+B[1]+U[1],A[2]+B[2]+U[2],A[3]+B[3]+U[3],A[4]+B[4]+U[4],A[5]+B[5]+U[5];a=E[0]+F[0]+V[0],E[1]+F[1]+V[1],E[2]+F[2]+V[2],E[3]+F[3]+V[3],E[4]+F[4]+V[4],E[5]+F[5]+V[5]
			if not any(Z)or not any(a):continue
			k,l,P=BT(Z,a)
			if side==0:
				if P%q!=0:continue
				P//=q
			if P<=0 or not Am(P):continue
			Q=BU(k,l)
			if side==1:
				if Q%q!=0:continue
				Q//=q
			if Q<=0 or(Q,P)in h or not Am(Q):continue
			R=[0]*12
			for(m,n)in H:
				j=N[m]
				for f in C(12):R[f]+=j[f]if n==1 else-j[f]
			if not L.is_ideal_primitive_oh(R[:6],R[6:]):continue
			if b==1:return R[:6],R[6:],Q,P
			if e is I:e=R[:6],R[6:],Q,P
			h.add((Q,P));d+=1
		return I if b==1 else(d,e)if b==2 else d
	def An(cside,cq,av,bv):
		I=cside;A=cq;C=L._ppowmod([0,1],A,L.H_COEFFS,A)
		if D(C)<2:C+=[0]*(2-D(C))
		C[1]=(C[1]-1)%A
		for E in L._roots_mod(L._pgcd(C,L.H_COEFFS,A),A):
			K,J=A8(av,E,A),A8(bv,E,A)
			if J==0:continue
			B=K*W(J,-1,A)%A
			if B==0:continue
			F=BC(v,E,A)
			if I==0 and(((V+1)*F+2)*B*B+1)%A==0:return 0,A,E,B
			if I==1:
				G=B*B%A;H=G*G%A;M=H*G%A;N=H*H%A
				if(N+(8+5*F)*M+(30+6*F)*H+(-4-7*F)*G+(-5-2*F))%A==0:return 1,A,E,B
	e=Q
	for(BZ,Ao,Ba,Bb,Bc)in p['live_ctrl']:
		if Ao<=i or A0(BZ,Ao,Ba,Bb,J)!=Bc:e=J
	for(Ap,Aq,AF,Bd,Ar,AG,Be)in p['live_d2']:
		s=A0(*Ap,first=1)
		if s is I or not i<AG<AF<Ap[1]or(s[2]if Aq==1 else s[3])%AF!=0:e=J;continue
		As=An(Aq,AF,s[0],s[1])
		if As is I:e=J;continue
		Bf,t=A0(*As,first=2)
		if Bf!=Bd or t is I or(t[2]if Ar==1 else t[3])%AG!=0:e=J;continue
		At=An(Ar,AG,t[0],t[1])
		if At is I or A0(*At,first=0)!=Be:e=J
	Z=[[P(ar)for ar in af.split()]for af in p['compact_100'].split(';')];Au=U(A[1]>i for A in Z)
	def j(rs):A=[aY[3]for aY in rs];F=D(A);C=G(A)/F;I=G((A-C)**2 for A in A)/(F-1);J=E(B,(I-C)/E(C*C,1e-12));return O(n=F,min=K(A),max=E(A),mean=H(C,2),fano_v_over_m=H(I/C,3),k_hat=H(J,4),n_zero=G(1 for A in A if A==0))
	AH,AI,AJ=j([aF for aF in Z if aF[0]==1 and aF[2]==0]),j([aG for aG in Z if aG[0]==1 and aG[2]==1]),j([aH for aH in Z if aH[0]==1]);AK,AL,AM=j([aI for aI in Z if aI[0]==0 and aI[2]==0]),j([aJ for aJ in Z if aJ[0]==0 and aJ[2]==1]),j([aK for aK in Z if aK[0]==0]);AN=L.rho_tab_1024();Av=D(AN);Bg=A7[o][y];Aw={P(ag):ah for(ag,ah)in A7[o][AQ].items()};Bh=log2(Aw[4]);AO={}
	def k(m,nA,nB,kA,kB,al):
		if nB<=59.:
			J=(m+al)/59.
			if J<=A:return A
			L=ceil(J*1024.);return AN[L]if L<Av else B
		G=nA,nB,kA,kB
		if G not in AO:
			D=(nB-59.)/24.;M=P(H(24.*(nA-59.)/(nB-59.)));R=[log1p(D/(59.+aZ*D))if aZ<M else B for aZ in C(24)];S=[log1p(D/(59.+a0*D))if a0>=M else B for a0 in C(24)];F=[[A]]
			for E in C(1,kA+1):F.append([a3/E for a3 in _conv(F[-1],R)])
			AO[G]=D,[(a1+a2,[(a5,a4)for(a5,a4)in d(F[a1]if a2==0 else _conv(F[a1],S))if a4])for a1 in C(kA+1)for a2 in C(kB+1)]
		D,T=AO[G];N=B;O=1024./59.;U=D*O
		for(E,V)in T:
			W=(m-E*59.-E*D)/D if D>0 else 1e9;X=(m-E*59.+al)*O
			for(Q,Y)in V:
				if Q<=W:I=ceil(X-Q*U);N+=Y*(A if I<=0 else AN[I]if I<Av else B)
		return K(A,N)
	b=[];f=[]
	for(Ax,Bi,x,x)in AY[AU]:
		Ay=L.b62d(Bi);AP=log2(Ay);A1=_SQ_NORMS_CACHE.get((Ax,Ay));Az=next((A for A in Bg if A[2]==Ax and A[0]<=AP<=A[1]),I)
		if Az is I or not A1:continue
		A_,Bj,A2,g,Bk,l,m,Bl,x=Az;n=59.+Bk*(g-59.)/96.;B0=log2(Aw[Bl]);B1=(B0-Bh)/12.;B2=(Bj+B0)/12.;B3=G(A*k(48*B2+B-(A_ if A2==Y else 0),n,g,l,m,ALPHA_F)*k(12*B2+LOG2U6+C-(A_ if A2==R else 0),n,g,l,m,ALPHA_G)for(A,B,C)in QBINS);B4=B5=B
		for(Bm,Bn)in A1:B6=Bm+48.*B1-(AP if A2==Y else B);B7=Bn+12.*B1-(AP if A2==R else B);B4+=k(B6,n,g,l,m,ALPHA_F)*k(B7,n,g,l,m,ALPHA_G);B5+=k(B6,n,g,l,m,L.ALPHA_F_RAW)*k(B7,n,g,l,m,L.ALPHA_G_RAW)
		b.append(B4/(D(A1)*B3));f.append(B5/(D(A1)*B3))
	AR=G(b)/D(b);Bo=G((A-AR)**2 for A in b)/E(1,D(b)-1);Bs=Bo/(AR*AR);AS=G(f)/D(f);Bp=G((A-AS)**2 for A in f)/E(1,D(f)-1);Bt=Bp/(AS*AS);A3=a(b);B8=a(f);Bq=G(A/(A+.5*(AJ[h]if B[0]==1 else AM[h]))for B in Z);Br=Ag and e and Au and D(Z)==100 and AJ[AV]==0 and AM[AV]==0 and E(AH[u],AI[u],AK[u],AL[u])<=1.5 and D(b)==24 and K(b)>=1.5 and K(f)>=2.5;A4=p.get('surv_100',[]);A5=[A4[at]for(at,au)in d(Z)if au[2]==0]if A4 else[];A6=[A4[aw]for(aw,ax)in d(Z)if ax[2]==1]if A4 else[];return O(ok=Br,tower_u5_p27631_r27481_deg8_frec_verified=Ag,monic_m1_deg8_h5w_surv_ctrl_mean=H(G(A5)/E(1,D(A5)),2),monic_m1_deg8_h5w_surv_d1_mean=H(G(A6)/E(1,D(A6)),2),monic_m1_deg8_h5w_surv_d1_over_ctrl=H(G(A6)/E(1,D(A6))/E(1e-09,G(A5)/E(1,D(A5))),4),all_100_above_fb_1024=Au,live_6_deg8_ideals_and_2hop_d2_exact=e,nb_r=Af,nb_alpha=B9,b_fb=i,deg8_f_ctrl=AH,deg8_f_d1_k1=AI,deg8_f_all=AJ,deg2_g_ctrl=AK,deg2_g_d1_k1=AL,deg2_g_all=AM,deg8_exp_p0_sum=H(Bq,2),deg8_direct_oe_p0=B,f_d1_over_ctrl=H(AI[h]/AH[h],4),g_d1_over_ctrl=H(AL[h]/AK[h],4),prod_24_lattices_n=D(b),prod_mu_cons_over_mod_min=H(A3[0],2),prod_mu_cons_over_mod_p10=H(A3[2],2),prod_mu_cons_over_mod_p25=H(A3[6],2),prod_mu_cons_over_mod_med=H(A3[12],2),prod_mu_true_over_mod_min=H(B8[0],1),prod_mu_true_over_mod_med=H(B8[12],1))
def verify_v40_descent(M,segs):return L.verify_v43_descent(__import__(__name__),M,segs)
def descbk_swap(s):global SIEVE_HIT,SCAN,ZERO,PAIR_SETUP,IDEAL,FIX_LLL,OMEGA_BITS;B=O(SIEVE_HIT=SIEVE_HIT,SCAN=SCAN,ZERO=ZERO,PAIR_SETUP=PAIR_SETUP,IDEAL=IDEAL,FIX_LLL=FIX_LLL,OMEGA_BITS=OMEGA_BITS);A=s if isinstance(s,O)else{aa:(DESCBK_TYPED if aa in s else DESCBK_NEW)[aa]for aa in DESCBK_TYPED};SIEVE_HIT,SCAN,ZERO,PAIR_SETUP,IDEAL,FIX_LLL,OMEGA_BITS=A[AH],A['SCAN'],A['ZERO'],A[AI],A['IDEAL'],A[AJ],A[AK];return B
def verify_descbk_constants(M,path,bare):
	R=bare;import re;x=re.compile('^\\s+[a-z][a-z0-9.]*\\b');J={};L=I;Q=1
	for a in AX.ko(path):
		V=a.strip()
		if V.startswith('//@seg'):L=V.split()[1];J.setdefault(L,0);Q=1;continue
		if V.startswith('//@end'):L=I;Q=1;continue
		A=a.split('//',1)[0].rstrip()
		if L is I or not A.strip()or A.endswith(':'):continue
		if A.strip().startswith('.rept'):Q=P(A.strip().split()[1]);continue
		if A.strip().startswith('.endr'):Q=1;continue
		if A.strip().startswith('.'):continue
		if x.match(A):J[L]+=Q
	W,X,Y=J.get('desc_zero16',-1),J.get('desc_scan8_u8',-1),J.get('desc_sieve_hit',-1);b=M.C_LAT_EMIT_V43B+M.C_LAT_APPLY;c=W>0 and X>0 and W/256.<=ZERO and X/64.<=SCAN;d=Y>0 and Y+M.C_LAT_APPLY<=SIEVE_HIT and b<=SIEVE_HIT;h=E(N(A[1])for A in M.load_data()[o][y]);j=6.*h*45.*mul_L(E(4,ceil(h/64.)));l=j<IDEAL;m=R[AO]<=PAIR_SETUP;n=R[AP]<=FIX_LLL;T=1<<K(OMEGA_BITS,24);F=bytearray([1])*(T+1);F[0]=F[1]=0
	for B in C(2,S.isqrt(T)+1):
		if F[B]:F[B*B::B]=bytes(D(C(B*B,T+1,B)))
	H=T+1
	while any(H%A==0 for A in C(2,S.isqrt(H)+1)):H+=1
	def A2(nb):
		A=0;B=H
		while B<1<<nb:B*=H;A+=1
		return A
	Z=O(_ECM_OMEGA);q=OMEGA_BITS<=24 and G(F)<=TDIV20_PRIMES and D(Z)>=40 and U(B>=A2(A[1])for(A,B)in Z.items());u=p==(SIEVE_HIT!=e)and r==(SCAN!=f and ZERO!=g)and t==(PAIR_SETUP!=z)and s==(IDEAL!=i)and v==(FIX_LLL!=A0)and w==(OMEGA_BITS!=A1);return O(ok=k(c and d and l and m and n and q and u),switches_on=O(HIT=p,SCANZERO=r,PAIR=t,IDEAL=s,FIXB=v,OMEGA=w),switches_match_constants=u,scan_zero=O(ok=k(c),scan=SCAN,scan_typed=f,counted_scan_ops_per_64_cells=X,zero=ZERO,zero_typed=g,counted_zero_ops_per_256_cells=W),sieve_hit=O(ok=k(d),used=SIEVE_HIT,typed=e,counted_emit_ops=Y,c_lat_apply=M.C_LAT_APPLY,unaccelerated_leaf_price=b),ideal=O(ok=k(l),used=IDEAL,typed=i,worst_case=j),pair_setup=O(ok=k(m),used=PAIR_SETUP,typed=z,proved_ir=R[AO]),boot_lll=O(ok=k(n),used=FIX_LLL,typed=A0,proved_ir=R[AP]),omega=O(ok=k(q),bits=OMEGA_BITS,typed_bits=A1,primes_up_to_2_bits=G(F),tdiv_primes_charged=TDIV20_PRIMES,smallest_prime_entering_ecm=H,ecm_targets_checked=D(Z)))
def verify_cells_theta(desc_data,rmax=117):
	L=desc_data;B=rmax;from math import sumprod as b;A=B*B;J=[0]*(A+1);F=0
	while F*F<=A:J[F*F]+=1 if F==0 else 2;F+=1
	Q=[aa for(aa,ab)in d(J)if ab];K=[0]*(A+1)
	for G in Q:
		for M in Q:
			if G+M>A:break
			K[G+M]+=J[G]*J[M]
	def R(a,b):B=lambda a:P.from_bytes(b''.join(A.to_bytes(5,n)for A in a),n);D=(B(a)*B(b)&(1<<40*(A+1))-1).to_bytes(5*(A+1),n);return[P.from_bytes(D[5*ah:5*ah+5],n)for ah in C(A+1)]
	S=R(R(K,K),K);T=[0]*(A+1);V=0
	for G in C(A+1):V+=S[G];T[G]=V
	N={ac:(b(S[:ac*ac+1],T[ac*ac::-1])-1)//2 for ac in C(2,B+1)};H={P(ad):ae for(ad,ae)in L[AQ].items()};E=a(A for A in H if 2<=A<=B);W=U(H[A]==N[A]for A in E);I=a({af[7]for af in L[y]});X=U(type(A)is P and A in H and 2<=A<=B for A in I);Y=O(H);Y[E[D(E)//2]]+=1;Z=not U(Y[A]==N[A]for A in E);return O(ok=W and X and Z and D(E)>=D(I),keys_checked=D(E),rmax=B,largest_key_checked=E[-1],all_keys_exact=W,plan_rows=D(L[y]),plan_radii=D(I),plan_radius_min=I[0],plan_radius_max=I[-1],plan_radii_all_checked_keys=X,control_off_by_one_detected=Z,unchecked_keys_above_rmax=[ag for ag in a(H)if ag>B],cells4=H[4],cells4_theta=N[4])
def verify_boot_gate_knapsack(desc,desc_data):
	x=desc_data;H=desc;from itertools import combinations_with_replacement as AZ;j=H['boot_knap'];Aa=j['params'];V=[M(aa)for aa in j['rows']];AJ,Ab,n,_=x[Ai];Ac,r=x['ea'];AK=x[y];l=N(K(A[0]for A in AK));AL=(AJ-l)/c;s=l+Ac*AL;i=l+Ab*AL;e=N(AJ);AM=Aa==O(LB=l,x_ea=s,nA=i,nB=e,kA=n,tau=r)and H['lb']==l;AN=V==a((A[0],A[1])for A in AK if A[2]==R)and V[0][0]==l and U(V[A][1]==V[A+1][0]for A in C(D(V)-1))and V[-1][1]>=e;t=lambda vs,lim:E([ax for((aw,ay),ax)in X(V,vs)if aw<lim-1e-09]+[B]);z=lambda vs:n*t(vs,i)+t(vs,e);A0,A1,AO=j[Bx],j[By],j[Bz];AP=z(A0)==H[AC]and z(A1)==H[AD]and z(AO)==H[AR]
	def A2(rs,cfg,vs,val,xe,na,nb_,k,ta,fl):B=cfg;C=floor if fl else N;H=[rs[al]for al in B['A']];A=rs[B[m]]if B[m]is not I else I;return D(H)<=k and U(A<na for(A,B)in H)and(A is I or A[1]>na and A[0]<nb_)and G(0 if A<xe else C(A)for(A,B)in H)+(C(E(A[0],na))if A else 0)<=C(ta)and G(F(vs[A])for A in B['A'])+(F(vs[B[m]])if A else 0)==val
	def A3(rs,vs,xe,na,nb_,k,ta,fl):
		H=[F(am)for am in vs];A=1
		for I in H:A=A*I.denominator//S.gcd(A,I.denominator)
		J=[P(an*A)for an in H];B=floor if fl else N;Q=[(0 if ao<xe else B(ao),ap)for((ao,aq),ap)in X(rs,J)if ao<na];R=[(B(E(ar,na)),au)for((ar,at),au)in X(rs,J)if at>na and ar<nb_]+[(0,0)];ta=B(ta);D=0
		for T in C(k+1):
			for K in AZ(Q,T):
				L=G(A for(A,B)in K)
				if L<=ta:
					M=G(A for(B,A)in K)
					for(U,O)in R:
						if L+U<=ta and M+O>D:D=M+O
		return F(D,A)
	A4=A5=A6=Q;AQ={}
	for(Ad,u,Ae,v)in(('w',A0,H[AS],H[AC]),('phi',A1,H[AW],H[AD]),('c_det',AO,H['c_det_one'],H[AR])):W,f=boot_gate_knapsack(V,u,l,s,i,e,n,r);A4&=Ae==(K(Ag(W),v)if j['on']else v);A5&=A2(V,f,u,W,s,i,e,n,r,J);A6&=W<=F(n)*F(t(u,i))+F(t(u,e));AQ[Ad]=O(bound=N(W),worst_row=v,ratio=N(W)/v,A_rows=[V[ak]for ak in a(f['A'])],B_row=V[f[m]]if f[m]is not I else I)
	A7=[ab for ab in C(D(V))if V[ab][0]>=120];AT=[V[ac]for ac in A7];p=[(5*ad,5*ad+5)for ad in C(7)];AU=[4.,2e1,7.,3e1,1e2,3.,4e1];AV=[(AT,[A0[ai]for ai in A7],122.,i,e,4,500),(AT,[A1[aj]for aj in A7],122.,i,e,4,430),(p,[4.,A,7.,2.,9.,3.,6.],3.,18.,3e1,3,41),(p,[4.,A,7.,5e1,9.,3.,6.],3.,18.,3e1,3,41),(p,AU,5.,2e1,3e1,3,41),(p,AU,5.,2e1,3e1,3,33),(p,[.25]*7,3.,18.,3e1,3,41)];A8=A9=AA=Q
	for(Y,b,d,L,T,Z,g)in AV:W,f=boot_gate_knapsack(Y,b,Y[0][0],d,L,T,Z,g);A8&=W==A3(Y,b,d,L,T,Z,g,J);A9&=A2(Y,f,b,W,d,L,T,Z,g,J);AA&=W<boot_gate_knapsack(Y,b,Y[0][0],d,L,T,Z,Z*L+T)[0]or b[0]==.25
	q=48;AB=[20261126]
	def h(n):AB[0]=AB[0]*0x5851f42d4c957f2d+0x14057b7ef767814f&(1<<64)-1;return(AB[0]>>33)%n
	def AE(bd,a,z,q):
		A=[av for av in bd if a<=av<=z]
		if A and h(2):return A[h(D(A))]
		B=ceil(a*q);return(B+h(floor(z*q)-B+1))/q
	w=Q;AF=AG=AH=0
	for _ in C(q):
		AX=[1+h(4)for ae in C(5+h(4))];AI=[0]
		for Af in AX:AI.append(AI[-1]+Af)
		o=[af/2. for af in AI];Y=[(o[ag],o[ag+1])for ag in C(D(AX))];AY=o[-1];T=AE(o,AY/2.+.5,AY,2);L=AE(o,.5,T-.5,2);d=AE(o,B,L,4);Z=1+h(3);g=(1+h(P(Z*L+T)))/2.;b=[(h(9)if h(5)else 0)/4. for ah in Y];W,f=boot_gate_knapsack(Y,b,B,d,L,T,Z,g);w&=W==A3(Y,b,d,L,T,Z,g,Q)and W>=A3(Y,b,d,L,T,Z,g,J)and A2(Y,f,b,W,d,L,T,Z,g,Q);AF+=W<boot_gate_knapsack(Y,b,B,d,L,T,Z,Z*L+T)[0];AG+=f[m]is not I;AH+=any(A%1!=0 for A in(d,L,T,g))
	w&=AF>=q//3 and AG>=q//4 and AH>=q//4;return O(ok=k(AM and AN and AP and A4 and A5 and A6 and A8 and A9 and AA and w),params_from_stored_data_ok=AM,rows_tile_ok=AN,worst_row_rederived_ok=AP,dp_reproduces_exported_ok=A4,argmax_certificate_ok=A5,bound_le_worst_row_exact_ok=A6,equals_exhaustive_enumeration_ok=A8,reduced_argmax_certificates_ok=A9,reduced_budgets_bind_ok=AA,random_instances_ok=w,random_instances=q,random_budget_binding=AF,random_with_B_child=AG,random_fractional_edges=AH,reduced_instances=D(AV),knap_on=j['on'],gate_binds=H[AS]<H[AC]and H[AW]<H[AD],bounds=AQ,tau=r,kA=n,x_ea=s,nA=i,nB=e)
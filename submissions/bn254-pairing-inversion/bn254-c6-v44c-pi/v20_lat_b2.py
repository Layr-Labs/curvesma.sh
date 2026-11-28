EP='emit9_step1'
EO='chord_enter'
EN='max_sa'
EM='div_err'
EL='t4_geom'
T='top_%d'
CU='maxb'
CT='maxb2'
CS='lvs'
CR='wheel_disp4'
CQ=bytearray
BY=list
BX=enumerate
AW=OverflowError
AI=dict
A1=False
A0=zip
s=.0
p=abs
n=int
o=min
h=any
g=tuple
a=float
Q=True
O=1.
N=max
H=pow
F=round
E=None
D=len
C=sum
A=range
import math as P,random,struct as BS,stor_bound as U
FP_CL=[20,72,97,79,77,75,73,71,69,67,65,97]
PHI_F3=(1,0,2,2,0,1),(1,2,0,2,1,2),(1,2,2,1,2,2)
WHEEL_OPS=AI(chord_enter=0,emit9_step1=9,emit9_step2=9,emit9_step3=9,emit9_last=7,emit8_step1=8,emit8_step2=8,emit8_step3=8,emit8_last=6,emit8_epoch_flush=10,wheel_disp4=4,wheel_skip1_stub=2,wheel_skip2_stub=2)
EMIT8_P_MIN=6369051672525773
EMIT8_FLUSH_AMORT=.225
BT=405,377,352,269,269,269
BU=537,500,468,357,357,357
R=g(2*A+1 for A in BT)
S=g(2*A+1 for A in BU)
I=[1]
for L in R[:-1]:I.append(I[-1]*L)
WB_MR=[1]
for L in S[:-1]:WB_MR.append(WB_MR[-1]*L)
PROD_MA_MR=I[-1]*R[-1]
PROD_MB_MR=WB_MR[-1]*S[-1]
A_OFF_MR=C(BT[A]*I[A]for A in A(6))
B_OFF_MR=C(BU[A]*WB_MR[A]for A in A(6))
def _mr_da(a):return C(a[A]*I[A]for A in A(6))&(1<<64)-1
def _mr_db5(b):return C(b[A]*WB_MR[A]for A in A(6))<<5&(1<<64)-1
def _mr_pack_a(a):return A_OFF_MR+C(a[A]*I[A]for A in A(6))&(1<<64)-1
def _mr_pack_b(b):return B_OFF_MR+C(b[A]*WB_MR[A]for A in A(6))
HT_SEED_BASE=1
PT_ORB=[]
PT_OLD=[]
PT_OLD_HAND=30
HT_NEW=AI(PRICE_ROW=6,PRICE_CAND=27,PRICE_HAND=33,PT_EXTRA=.0054,PT_ROWS=1.0065,PT_CANDS=1.015)
def Gd(v):A=HT_NEW;return C(A*B for(A,B)in A0([A['PRICE_ROW']-WHEEL_OPS[CR],A['PRICE_CAND'],FP_CL[2]-PT_OLD_HAND+A['PRICE_HAND']]+FP_CL[3:],v))
def _st_f3(v):A,G,H,B,I,C,D,J,K,E,L,F=v;M=(A+2*H+2*B+C)%3;N=(A+2*G+2*B+I+2*C)%3;O=(A+2*G+2*H+B+2*I+2*C)%3;P=(D+2*K+2*E+F)%3;Q=(D+2*J+2*E+L+2*F)%3;R=(D+2*J+2*K+E+2*L+2*F)%3;return M+3*(N+3*(O+3*(P+3*(Q+3*R))))
def _add_st_f3(s1,s2):
	B=0;C=1
	for D in A(6):B+=(s1%3+s2%3)%3*C;s1//=3;s2//=3;C*=3
	return B
def _bad_st_f3(st):A=st;B=A%3;A//=3;C=A%3;A//=3;D=A%3;A//=3;E=A%3;A//=3;F=A%3;G=A//3;return B==0 and E==0 or C==0 and F==0 or D==0 and G==0
def _bad_f3(v):A,G,H,B,I,C,D,J,K,E,L,F=v;return(A+2*H+2*B+C)%3==0 and(D+2*K+2*E+F)%3==0 or(A+2*G+2*B+I+2*C)%3==0 and(D+2*J+2*E+L+2*F)%3==0 or(A+2*G+2*H+B+2*I+2*C)%3==0 and(D+2*J+2*K+E+2*L+2*F)%3==0
def wheel_prices_from_kernels(path):B=count_kernel_segments(path);return{j:B.get(j,-1)for j in WHEEL_OPS}
def BV(v,o=0):return C(FP_CL[A+o]*B for(A,B)in BX(v))
from fractions import Fraction as z
from math import sqrt,ceil,floor,log2
import model_bn254_v6 as M,descent_bn254_v6 as DESC
K=(1<<64)-1
EF=-(1<<63)
EG=(1<<63)-1
Q32=1<<32
Q38=1<<38
def G(x):
	if x<0 or x>K:raise AW(f"u64 overflow: {x}")
	return x
def J(x):
	if x<EF or x>EG:raise AW(f"i64 overflow: {x}")
	return x
def count_kernel_segments(path):
	import re;G=re.compile('^\\s+[a-z][a-z0-9.]*\\b');C={};A=E
	for F in U.ko(path):
		D=F.strip()
		if D.startswith('//@seg'):A=D.split()[1];C.setdefault(A,0);continue
		if D.startswith('//@end'):A=E;continue
		B=F.split('//',1)[0].rstrip()
		if A is E or not B.strip()or B.endswith(':')or B.strip().startswith('.'):continue
		if G.match(B):C[A]+=1
	return C
def fp_prices_from_kernels(path):
	B=count_kernel_segments(path);D=[B['k0_body']+B['k0_back'],B['top_2']+B['body_2']+B['k1_exit']]
	for C in A(2,11):D.append(B[T%(C+1)]+B['body_%d'%(C+1)]+B[T%C])
	return D,B
def Ar(a,b):G(a);G(b);return a*b>>64
V=[0]*512+[ceil(sqrt(aa+1)*2.**52.5)+(1<<31)for aa in A(512,2048)]
def Ge(r2):
	A=r2
	if A<1 or A>=1<<40:raise AW(f"u40 range error: {A}")
	C=64-A.bit_length();B=C&~1;D=A<<B>>53;return G(V[D]>>(B>>1))
def EH():
	D=[]
	for B in A(24,65):
		if B==64:D.append((64<<11,0,0));continue
		for C in A(1024,2048):
			if B<=53:E,F=C<<53-B,(C+1<<53-B)-1
			else:
				if C&(1<<B-53)-1:continue
				E=F=C>>B-53
			D.append((B<<11|C,E,F))
	return D
f=AI(n=0,fail=0,maxb=-999.,n2=0,fail2=0,maxb2=-999.)
AE=[0]*((64<<11)+1)
for(W,d,X)in EH():AE[W]=P.isqrt(X+1<<52)+1
def As(r2):
	if r2<1 or r2>=(1<<40)+1:raise AW(f"u40 range error: {r2}")
	A=r2-1;B=64-A.bit_length();C=A<<(B&63)&K;return AE[B<<11|C>>53]
def At(mem,off,w_lp):mem[off]=mem[off]+(w_lp&65535)&65535
B=P.sumprod
def _mv6(Mat,v):A=Mat;return B(A[0],v),B(A[1],v),B(A[2],v),B(A[3],v),B(A[4],v),B(A[5],v)
def AV(G,v):return B(v,(B(G[0],v),B(G[1],v),B(G[2],v),B(G[3],v),B(G[4],v),B(G[5],v)))
def _eval_poly(cs,s,p):
	A=0
	for B in reversed(cs):A=(A*s+B)%p
	return A
def EI(v,s,p):
	J=H(3,-1,p);D=s*s%p;E=D*s%p;F=E*s%p;K=F*s%p;C=0;L=M.ZK_NUM
	for G in A(6):
		I=v[G]
		if I:B=L[G];N=(B[0]+B[1]*s+B[2]*D+B[3]*E+B[4]*F+B[5]*K)%p;C=(C+I*N)%p
	return C*J%p
def _sqrt_mod(n,p):
	n%=p
	if n==0:return 0
	if H(n,(p-1)//2,p)!=1:return-1
	if p%4==3:A=H(n,(p+1)//4,p);return o(A,p-A)
	B,J=p-1,0
	while B&1==0:B>>=1;J+=1
	E=2
	while H(E,(p-1)//2,p)!=p-1:E+=1
	K,F,C,A=J,H(E,B,p),H(n,B,p),H(n,(B+1)//2,p)
	while C!=1:
		D,G=C,0
		while D!=1:D=D*D%p;G+=1
		I=H(F,1<<K-G-1,p);K=G;F=I*I%p;C=C*F%p;A=A*I%p
	return o(A,p-A)
def Y(p):
	def L(a,b):
		B,C=D(a),D(b);E=[0]*(B+C-1)
		for F in A(B):
			G=a[F]
			if G:
				for H in A(C):E[F+H]+=G*b[H]
		return[k%p for k in E]
	def I(a,m):
		a=a[:];B=D(m);E=1 if m[-1]==1 else H(m[-1],-1,p)
		while D(a)>=B:
			G=a[-1]if E==1 else a[-1]*E%p;F=D(a)-B
			for C in A(B):a[F+C]=(a[F+C]-G*m[C])%p
			while a and a[-1]==0:a.pop()
		return a
	def T(a,b):
		a=a[:];C=D(b);I=H(b[-1],-1,p);B=[0]*(D(a)-C+1)
		while D(a)>=C:
			G=a[-1]*I%p;E=D(a)-C;B[E]=G
			for F in A(C):a[E+F]=(a[E+F]-G*b[F])%p
			while a and a[-1]==0:a.pop()
		while B and B[-1]==0:B.pop()
		return B
	def N(b,e,m):
		if D(m)==7 and m==J:
			A=M._ppow_h(b,e,p)
			while A and A[-1]==0:A.pop()
			return A
		B=[1]
		while e:
			if e&1:B=I(L(B,b),m)
			b=I(L(b,b),m);e>>=1
		return B
	def O(a,b):
		while b:a,b=b,I(a,b)
		if a:A=H(a[-1],-1,p);a=[l*A%p for l in a]
		return a
	J=[j%p for j in M.H_COEFFS];C=N([0,1],p,J)
	if D(C)<2:C+=[0]*(2-D(C))
	C[1]=(C[1]-1)%p
	while C and C[-1]==0:C.pop()
	P=O(J,C)
	if D(P)<4:return[]
	U=random.Random(p);G=[P];F=[]
	while G:
		B=G.pop()
		if D(B)==2:F.append(-B[0]*H(B[1],-1,p)%p)
		elif D(B)==3:R=H(2*B[2],-1,p);S=_sqrt_mod((B[1]*B[1]-4*B[2]*B[0])%p,p);F.append((-B[1]+S)*R%p);F.append((-B[1]-S)*R%p)
		else:
			while Q:
				V=U.randrange(p);E=N([V,1],(p-1)//2,B)
				if E:E[0]=(E[0]-1)%p
				while E and E[-1]==0:E.pop()
				K=O(B,E)
				if 1<D(K)<D(B):G.append(K);G.append(T(B,K));break
	F.sort();return F
def _find_root_h(p):A=Y(p);return A[0]if A else E
EJ=20,10,0,50,40,30
Z=512,512,512,512,512,512
Gf=512,512,512,0,512,512
b=51,30,20,10,0,41
c=1024,1024,512,512,512,512
EK=0,5,1,2,3,4
e=C(2*Z[A]<<EJ[A]for A in A(6))
i=C(2*c[A]<<b[A]for A in A(6))
def Gg(v,sh,bi):
	C=0
	for B in A(6):C|=v[B]+bi[B]<<sh[B]
	return C
def _cone_of(v):
	A,B,C=v[3],-v[4],v[5];F=A+B+C
	if F==0:return-1
	if F>0:G=1
	else:G=0;A,B,C=-A,-B,-C
	if A>B and A>C:E,D=0,v[6]
	elif B>A and B>C:E,D=1,v[6]+v[8]+v[10]
	elif C>A and C>B:E,D=2,v[6]+v[7]+v[10]
	else:return-1
	if D==0:return-1
	return E<<2|G<<1|n(D>0)
def _canon_cone_ref(a,b,S,S2):
	D=a[3]-a[4]+a[5];M=(1,)if D>0 else(-1,)if D<0 else(1,-1);F=E;H=E
	for(B,G)in((a,b),(_mv6(S,a),_mv6(S,b)),(_mv6(S2,a),_mv6(S2,b))):
		I=1
		for J in EK:
			if G[J]!=0:I=1 if G[J]>0 else-1;break
		C=g(I*A for A in G);N=C[0],C[1],C[2],C[3],C[4],C[5]
		for K in M:
			A=B if K==1 else(-B[0],-B[1],-B[2],-B[3],-B[4],-B[5]);L=(K*D,A[3],A[3],A[4],A[5],A[0],A[1],A[2])+N
			if F is E or L>F:F=L;H=A,C
	return H
def Gh(G_sym,G_SC):
	G=G_sym;E=[z(0)]*6;D=[[z(1 if j==k else 0)for k in A(6)]for j in A(6)]
	for B in A(6):
		H=z(G[B][B])-C(D[B][A]*D[B][A]*E[A]for A in A(B));assert H>0;E[B]=H
		for I in A(B+1,6):D[I][B]=(z(G[I][B])-C(D[I][A]*D[B][A]*E[A]for A in A(B)))/H
	M=5*G_SC;F=[0]*6;J=[10**18];K=[0]
	def L(idx,rem):
		B=idx;H=C(D[A][B]*F[A]for A in A(B+1,6));I=rem/E[B];M=z(P.isqrt(I.numerator*I.denominator)+1,I.denominator);S=E[B]
		for N in A((-H-M).__floor__(),(-H+M).__ceil__()+1):
			O=N+H;Q=rem-S*O*O
			if Q>=0:
				F[B]=N
				if B==0:
					if h(F):
						R=AV(G,F);K[0]+=1
						if R<J[0]:J[0]=R
				else:L(B-1,Q)
	L(5,z(M));return J[0],K[0]
BW=32
HT_K=8
def ht_stratum(s_idx):
	C=s_idx;A=M.OP['p_lat'];D=M.OP['lx2'];F=A+(D-A)*C/BW;G=A+(D-A)*(C+1)/BW;E=n(2.**F)|1;B=n(2.**G)
	if B&1==0:B-=1
	return E,(B-E)//2+1
def ht_stream(hb):
	K,B,B,B,B,L=M.build_exact_int_geom(M.load_data()[EL]);N=n(F(1764*566.26**2*L));C=[];D=[]
	for E in A(BW):
		O,G=ht_stratum(E);H=0;I=[]
		for P in A(HT_K):Q,J,R,S,T,U,B=ht_orbit(hb,E,P,O,G,K);H+=Q;I.append((J,R));D.append((J,[j[:]for j in S],T,U,N))
		C.append((G,H,I))
	return C,D
def ht_orbit(hb,s_idx,kd,L_i,M_s,G_sym):
	F=random.Random(hb<<20|s_idx<<8|kd);N=0
	while Q:
		N+=1;B=L_i+2*M.rand_below(F.getrandbits,M_s)
		if B&7 not in(1,3)or not M._isp(B):continue
		G=M._ppow_h([0,1,0,0,0,0],B,B);G[1]=(G[1]-1)%B;O=M.H_COEFFS if not h(G)else M._pgcd(G,M.H_COEFFS,B)
		if D(O)<=1:continue
		P=M._roots_mod(O,B)
		if not P:continue
		I=[];f=M.U_BN254%B
		for R in P:
			W=_eval_poly(M.W_COEFFS,R,B);S=-((M.U_BN254+1)%B*W+2)%B
			if S!=0:
				E=_sqrt_mod(H(S,-1,B),B)
				if E>0:I.append((R,E))
		if I:X=D(I)//3;Y,E=F.choice(I);Z=E if F.randrange(2)==0 else B-E;break
	a=H(3,-1,B);J=[1]
	for b in A(5):J.append(J[-1]*Y%B)
	T=[C(M.ZK_NUM[j][A]*J[A]for A in A(6))*a%B for j in A(6)];U=T+[-Z*k%B for k in T];c=H(U[0],-1,B);K=[[B if l==0 else 0 for l in A(12)]]
	for V in A(1,12):
		L=-U[V]*c%B
		if L>B//2:L-=B
		K.append([L if m==0 else 1 if m==V else 0 for m in A(12)])
	b,d,e=DESC.bareiss_lll_12d(K,G_sym,1764,1000,A1,9,10,(99,100),Q);return N,B,X,K,d,e,F
def verify_v20_lat_b2_and_srt(seed=20261121):
	BZ=M.load_data()[EL];q,EQ,ER,R,i,AX=M.build_exact_int_geom(BZ);ES=R[3]==[0,0,0,0,-1,0]and R[4]==[0,0,0,0,0,-1]and R[5]==[0,0,0,1,0,0];CV,Gi=Gh(q,AX);Gj=1000*CV;Au=z(Gj)*z(13,20)**11;ET=n(F(1764*595.**2*AX));CW=n((z(Q38*ET)/Au).__ceil__());CX=(Q38*CW>>64)+1;Gk=P.isqrt(n((z(ET)/Au).__ceil__()))+1;EU=n((z(3,2)**11*Gk).__ceil__());CY=11*((1<<31)+1)*EU;EV=As(CX);Gl=CY+11*EU+1+EV+Q32>>32;CZ=EV+Q32>>13;EW=CZ*CZ;t=sorted(EH(),key=lambda e:e[1]);Gm=all(AE[t[A][0]]<=AE[t[A+1][0]]for A in A(D(t)-1))and all(t[A][2]+1==t[A+1][1]for A in A(D(t)-1))and t[0][1]==0 and t[-1][2]==(1<<40)-1;EX=C(1 for(A,C,B)in t if AE[A]*AE[A]<=B+1<<52);EY=C(1 for(A,B,C)in t if AE[A]>Ge(B+1));Gn=N(AE[A]/(sqrt(B+O)*67108864.)-O for(A,B,C)in t);EZ=CV==q[0][0]==753729 and Gm and Au>32*AX and CW<=K and CX<1<<40 and CY<1<<63 and Gl<1<<22 and CZ<1<<32 and EW<=K;Ca=random.Random(seed^20903);Cb=0;Ba=s
	for Ea in BY(A(1,2049))+[(1<<bH)-1 for bH in A(8,41)]+[1<<bh for bh in A(8,40)]+[1+M.rand_below(Ca.getrandbits,(1<<40)-1)for au in A(19000)]:
		Eb=As(Ea);Cc=sqrt(Ea)*(1<<26)
		if Eb<=Cc:Cb+=1
		Ec=(Eb-Cc)/Cc
		if Ec>Ba:Ba=Ec
	Cd=0;Ed=[0]*256;Ce=[0]*256
	for Ee in A(4096):
		Av=Ca.randrange(256);Ef=65535 if Ee<16 else 32768 if Ee<32 else Ca.randrange(65536);At(Ed,Av,Ef);Cf=Ce[Av];Cf=Cf+(Ef&65535)&65535;Ce[Av]=Cf
		if Ed[Av]!=Ce[Av]:Cd+=1
	Aw=[];X=101
	while D(Aw)<6:
		if M._isp(X)and X not in(2,3,19):
			Eg=[bi for bi in A(X)if _eval_poly(M.H_COEFFS,bi,X)==0]
			if D(Eg)>=3:
				Eh=set();Bb=H(3,-1,X)
				for AB in Eg:
					if AB in Eh:continue
					AJ=_eval_poly(M.SIGMA_NUM,AB,X)*Bb%X;AY=_eval_poly(M.SIGMA_NUM,AJ,X)*Bb%X
					if D({AB,AJ,AY})!=3:continue
					Eh.update((AB,AJ,AY));Cg=_eval_poly(M.W_COEFFS,AB,X);assert _eval_poly(M.W_COEFFS,AJ,X)==Cg and _eval_poly(M.W_COEFFS,AY,X)==Cg;Ei=-((M.U_BN254+1)%X*Cg+2)%X
					if Ei==0:continue
					Ax=_sqrt_mod(H(Ei,-1,X),X)
					if Ax>0:
						Aw.append((X,(AB,AJ,AY),Ax))
						if D(Aw)==6:break
		X+=2
	def Ej(a,b,Wa,Wb):
		for B in A(6):
			C=Wa*AV(EQ[B],a)+Wb*AV(EQ[B],b)
			if C<0:return A1
			if C==0:
				D,E=_mv6(ER[B],a),_mv6(ER[B],b)
				if _canon_cone_ref(D,E,R,i)<_canon_cone_ref(a,b,R,i):return A1
		return Q
	def Go(a,b):
		D=_canon_cone_ref(a,b,R,i);A=0
		for(B,C)in((a,b),(_mv6(R,a),_mv6(R,b)),(_mv6(i,a),_mv6(i,b))):
			if _canon_cone_ref(B,C,R,i)==D:
				if B in(a,g(-A for A in a))and C in(b,g(-A for A in b)):A+=1
		return A
	Ch,Ci=1,1;AZ=18*AX;Cj=[]
	for Gp in A(-2,3):
		for Gq in A(-2,3):
			for Gr in A(-2,3):
				for Gs in A(-2,3):
					for Gt in A(-2,3):
						for Gu in A(-3,4):
							Ck=Gu,Gt,Gs,Gr,Gq,Gp
							if not h(Ck):continue
							Ek=AV(q,Ck)
							if Ek<=AZ:Cj.append((Ck,Ek))
	Cl={};El={}
	for(Gv,Gw)in Cj:
		for(Gx,Gy)in Cj:
			if Gw+Gy>AZ:continue
			A2,A3=_canon_cone_ref(Gv,Gx,R,i)
			if(A2,A3)not in Cl:
				Bc=Go(A2,A3);Cl[A2,A3]=Bc
				if Ej(A2,A3,Ch,Ci):El[A2,A3]=Bc
	Cm=0;Em=0;En=0
	for(j,(AB,AJ,AY),Ax)in Aw:
		Eo={};Ep={}
		for((A2,A3),Bc)in Cl.items():
			Aa=0
			for Eq in(AB,AJ,AY):
				u=EI(A2,Eq,j);Ab=Ax*EI(A3,Eq,j)%j
				if u==Ab:Aa+=1
				if u==(j-Ab)%j:Aa+=1
			if Bc==3:Aa//=3
			if Aa:
				Eo[A2,A3]=Aa
				if(A2,A3)in El:Ep[A2,A3]=Aa
		Bb=H(3,-1,j);Cn=[1]
		for _ in A(5):Cn.append(Cn[-1]*AB%j)
		Er=[C(M.ZK_NUM[aw][A]*Cn[A]for A in A(6))*Bb%j for aw in A(6)];Es=Er+[-Ax*bj%j for bj in Er];Gz=H(Es[0],-1,j);AK=[[j if bk==0 else 0 for bk in A(12)]]
		for U in A(1,12):
			AL=-Es[U]*Gz%j
			if AL>j//2:AL-=j
			AK.append([AL if bI==0 else 1 if bI==U else 0 for bI in A(12)])
		DESC.bareiss_lll_12d(AK,q,Ch,Ci,A1,9,10);G_=[[C(AK[ax][C]*q[C][E]*AK[bl][E]for C in A(6)for E in A(6))+C(AK[ax][6+C]*q[C][E]*AK[bl][6+E]for C in A(6)for E in A(6))for bl in A(12)]for ax in A(12)];v=[1]+[0]*12;Bd=[[0]*12 for ay in A(12)]
		for U in A(12):
			for b in A(U+1):
				AC=G_[U][b]
				for B in A(b):AC=(v[B+1]*AC-Bd[U][B]*Bd[b][B])//v[B]
				if b<U:Bd[U][b]=AC
				else:v[U+1]=AC
		H0=[[J((Bd[aA][bm]<<32)//v[bm+1])if aA>bm else 0 for bm in A(12)]for aA in A(12)];H1=[G(o(K,(v[aB+1]<<64)//(v[aB]*AZ)))for aB in A(12)];H2=[G(((v[aC]*AZ<<38)+v[aC+1]-1)//v[aC+1])for aC in A(12)];Co=[0]*12
		for B in A(11,-1,-1):Ay=P.isqrt((v[B]*AZ+v[B+1]-1)//v[B+1])+1;Co[B]=Ay+C((Co[A]+1)//2 for A in A(B+1,12))
		H3=[C(Co[A]for A in A(aD+1,12))+1 for aD in A(12)];Be={};Bf={};Et=[0]*12
		def Eu(k,rem,az,cur):
			H=rem;B=0
			for L in A(k+1,12):B=J(B-J(H0[L][k]*Et[L]))
			T=G(Ar(H,H2[k])+1);M=As(T);F=H3[k];U=(1 if k==0 else 0)if az else J(B-F-M+Q32-1)>>32;V=J(B+F+M)>>32;N=AK[k];W=H1[k]
			for C in A(U,V+1):
				if k==0:
					O=[cur[cc]+C*N[cc]for cc in A(12)];I,K=g(O[:6]),g(O[6:])
					if h(I)and h(K)and AV(q,I)+AV(q,K)<=AZ:
						D,E=_canon_cone_ref(I,K,R,i);Be[D,E]=Be.get((D,E),0)+1
						if Ej(D,E,Ch,Ci):Bf[D,E]=Bf.get((D,E),0)+1
				else:
					P=G(p(J((C<<32)-B)));X=P-F if P>F else 0;Q=X>>13;S=Ar(G(Q*Q),W)
					if S<=H:Et[k]=C;Eu(k-1,H-S,az and C==0,[cur[cd]+C*N[cd]for cd in A(12)])
		Eu(11,Q38,Q,[0]*12)
		if Be!=Eo or Bf!=Ep:Cm+=1
		Em+=C(Be.values());En+=C(Bf.values())
	def Ev(tb):
		A=1
		for B in EK:
			if tb[B]!=0:A=1 if tb[B]>0 else-1;break
		return g(A*B for B in tb)
	def H4(a,b):
		G,L=_mv6(R,a),_mv6(i,a);C=a[3]-a[4]+a[5];H=(1,)if C>0 else(-1,)if C<0 else(1,-1)
		if G==a:M=N(((A*C,A*a[3],A*a[3],A*a[4],A*a[5],A*a[0],A*a[1],A*a[2]),g(A*B for B in a))for A in H)[1];O=N(Ev(A)for A in(b,_mv6(R,b),_mv6(i,b)));return 3,M,O
		F=E;D=0;I=E
		for(P,B)in((0,a),(1,G),(2,L)):
			for J in H:
				A=B if J==1 else(-B[0],-B[1],-B[2],-B[3],-B[4],-B[5]);K=J*C,A[3],A[3],A[4],A[5],A[0],A[1],A[2]
				if F is E or K>F:F=K;D=P;I=A
		Q=b if D==0 else _mv6(R,b)if D==1 else _mv6(i,b);return D,I,Ev(Q)
	def Cp(base,B0,l0,h0,T12_chk=E,verify_ref=Q,p_div=E,s0_div=E,rho_div=E):
		U=verify_ref;B=p_div
		if l0>h0:return AI(lvs=0,fast=0,bnd=0,segs=0,bspl=0,fspl=0,cross=0,clips=0,mm=0,ae=0,div_err=0,max_sa=0)
		O=h0-l0+1;Q=[base[bY]+l0*B0[bY]for bY in A(12)];F=_mr_pack_a(Q[:6])&K;D=_mr_pack_b(Q[6:])<<5&K;Y=_mr_da(B0[:6]);Z=_mr_db5(B0[6:]);V=G=W=0;I=E
		if B is not E:
			a=H(3,-1,B);S=[1]
			for o in A(5):S.append(S[-1]*s0_div%B)
			I=[C(M.ZK_NUM[bZ][A]*S[A]for A in A(6))*a%B for bZ in A(6)];b=[-rho_div*b0%B for b0 in I];c=I+b;d=P.sumprod
		for e in A(O):
			if U or B is not E:
				T=[Q[b1]+e*B0[b1]for b1 in A(12)]
				if I is not E and d(T,c)%B!=0:W+=1
				J,L=g(T[:6]),g(T[6:])
				if h(p(J[A])>BT[A]for A in A(6))or h(p(L[A])>BU[A]for A in A(6)):G+=1
				else:
					f=_mr_pack_a(J);N=_mr_pack_b(L)
					if F!=f or D&31!=0 or D>>5!=N or F>=PROD_MA_MR or N>=PROD_MB_MR or N>=1<<59:G+=1
					if U:
						X,j,k=H4(J,L);l=(X<<64|D)>>5&K
						if l!=X<<59|N:G+=1
						m,n=_canon_cone_ref(J,L,R,i)
						if j!=m or k!=n:V+=1
			F=F+Y&K;D=D+Z&K
		return AI(lvs=O,fast=O,bnd=0,segs=1,bspl=0,fspl=0,cross=0,clips=0,mm=V,ae=G,div_err=W,max_sa=1)
	Cq=[((1,0,0,0,0,0),(0,1,0,0,0,0)),((5,0,0,0,0,0),(1,2,-3,4,-1,2)),((-7,0,0,0,0,0),(0,0,1,-2,3,-4)),((1,1,1,0,0,0),(0,1,-1,2,0,-1))]
	for(Cr,H5)in[((0,1,0,0,0,0),(0,1,0,0,0,0)),((0,0,0,1,0,0),(1,-2,3,0,1,-1)),((1,0,1,0,1,0),(0,0,1,0,0,-2))]:H6=_mv6(R,Cr);H7=_mv6(i,Cr);Cs=g(Cr[A]+H6[A]+H7[A]for A in A(6));assert _mv6(R,Cs)==Cs;Cq.append((Cs,H5))
	Bg=0
	for(Az,Ew)in Cq:
		Ex=set();Ct,Cu=Az,Ew
		for _ in A(3):
			for H8 in(1,-1):
				for H9 in(1,-1):Ex.add(_canon_cone_ref(g(H8*A for A in Ct),g(H9*A for A in Cu),R,i))
			Ct,Cu=_mv6(R,Ct),_mv6(R,Cu)
		if D(Ex)!=1:Bg+=1
		if Az[3]==0 and Az[4]==0 and Az[5]==0:
			Ey=Cp(BY(Az)+BY(Ew),[1,0,0,0,0,0,0,1,0,0,0,0],0,4,verify_ref=Q)
			if Ey['mm']or Ey['ae']:Bg+=1
	A4=random.Random(seed);c=n(F(1764*566.26**2*AX));AM=[P.pi**(aa/2.)/P.gamma(aa/2.+O)for aa in A(13)];HA=sqrt(1764.**6*1e3**6*M._det_mat(BZ['GA'])**2)*AX**6;Bh=49,51,53,55,57,59,61;Cv=0;Ez=0;Bi=0;Bj=0;Bk=0;A_=0;Bl=0;Bm=0;Bn=0;Cw=[];Cx=[];Cy=[];Cz=[];C_=0;D0=0;D1=0;D2=0;E_=0;Bo=[];D3=999.;Bp=0;Bq=0;Br=0;Bs=0;Bt=0
	for HB in Bh:
		S=(1<<HB)+1
		while Q:
			if M._isp(S):
				B0=_find_root_h(S)
				if B0 is not E:
					HC=_eval_poly(M.W_COEFFS,B0,S);F0=-((M.U_BN254+1)%S*HC+2)%S
					if F0!=0:
						Bu=_sqrt_mod(H(F0,-1,S),S)
						if Bu>0:break
			S+=2
		HD=H(3,-1,S);D4=[1]
		for _ in A(5):D4.append(D4[-1]*B0%S)
		F1=[C(M.ZK_NUM[aE][A]*D4[A]for A in A(6))*HD%S for aE in A(6)];F2=F1+[-Bu*bn%S for bn in F1];HE=H(F2[0],-1,S);A5=[[S if bo==0 else 0 for bo in A(12)]]
		for U in A(1,12):
			AL=-F2[U]*HE%S
			if AL>S//2:AL-=S
			A5.append([AL if bJ==0 else 1 if bJ==U else 0 for bJ in A(12)])
		DESC.bareiss_lll_12d(A5,q,1764,1000,A1,9,10,(99,100));F3=[(_mv6(q,aF[:6]),_mv6(q,aF[6:]))for aF in A5];HF=[[1764*P.sumprod(F3[aG][0],A5[bp][:6])+1000*P.sumprod(F3[aG][1],A5[bp][6:])for bp in A(12)]for aG in A(12)];I=[1]+[0]*12;B1=[[0]*12 for aH in A(12)]
		for U in A(12):
			for b in A(U+1):
				AC=HF[U][b]
				for B in A(b):AC=(I[B+1]*AC-B1[U][B]*B1[b][B])//I[B]
				if b<U:B1[U][b]=AC
				else:I[U+1]=AC
		F4=o(I[A+1]/I[A]for A in A(12));D3=o(D3,log2(F4));assert F4>=a(Au)*.999999;Ac=n((42e2*S*HA/AM[12])**(O/6.));D5=[[J((B1[aI][bq]<<32)//I[bq+1])if aI>bq else 0 for bq in A(12)]for aI in A(12)];F5=[G(o(K,(I[aJ+1]<<64)//(I[aJ]*Ac)))for aJ in A(12)];HG=[G(((I[aK]*Ac<<38)+I[aK+1]-1)//I[aK+1])for aK in A(12)];D6=[0]*12
		for B in A(11,-1,-1):Ay=P.isqrt((I[B]*Ac+I[B+1]-1)//I[B+1])+1;D6[B]=Ay+C((D6[A]+1)//2 for A in A(B+1,12))
		HH=[C(D6[A]for A in A(aL+1,12))+1 for aL in A(12)];Bv=[0]*12;D7=[0]*12;Bw=[0]*12;B2=[0];F6=[[D5[br][aM]for br in A(aM+1,12)]for aM in A(12)];D8=(I[2]<<38)//(I[1]*Ac);HI=2*D8;F7=P.sumprod;HJ=[(_mr_da(aN[:6]),_mr_db5(aN[6:]))for aN in A5];HK=A_OFF_MR&K,B_OFF_MR<<5&K;Ad=[_st_f3(A5[aO])for aO in A(12)];Bx=[E]+[[[_add_st_f3(bS,_add_st_f3(Ad[bs],Ad[bs])if bK==2 else Ad[bs]if bK==1 else 0)<<(3 if bs==1 else 0)for bS in A(729)]for bK in A(3)]for bs in A(1,12)];HL=[Bx[1][1][aP>>3]|((1 if _bad_st_f3(aP>>3)else 0)|(2 if _bad_st_f3(_add_st_f3(aP>>3,Ad[0]))else 0)|(4 if _bad_st_f3(_add_st_f3(aP>>3,_add_st_f3(Ad[0],Ad[0])))else 0))if aP>>3<729 else 0 for aP in A(8192)]
		def D9(k,rem,az,cur,c0_passed=E,st_cur=0,bmask_passed=E,c21=E):
			Y=bmask_passed;X=st_cur;W=c0_passed;L=cur;J=rem;F=c21;nonlocal Bi,Bj,Bk,A_,Bt,Bl,Bm,Bn;Bv[k]+=1;D=-F7(F6[k],Bw[k+1:])if k<11 else 0
			if D<EF or D>EG:raise AW(f"i64 ck overflow: {D}")
			if k==0 and W is not E and W!=D:Bt+=1
			O=(J*HG[k]>>64)+1
			if O<1 or O>=1<<40:raise AW(f"u40 r2 overflow: {O}")
			Z=As(O);M=HH[k];B=(1 if k==0 else 0)if az else D-M-Z+Q32-1>>32;G=D+M+Z>>32
			if B>G:return
			N=A5[k];H=HJ[k]
			if k==0:
				if F is not E:
					i=F[0]+B*H[0]&K;a=F[1]+B*H[1]&K;b=[L[ce]+B*N[ce]for ce in A(12)]
					if i!=_mr_pack_a(b[:6])or a&31!=0 or a>>5!=_mr_pack_b(b[6:]):Bn+=1
				if Y is not E:
					for c in A(B,G+1):
						Bm+=1
						if bool(Y>>c%3&1)!=_bad_f3([L[cf]+c*N[cf]for cf in A(12)]):Bl+=1
				I=Cp(L,N,B,G,T12_chk=Ac,verify_ref=B2[0]<400,p_div=S,s0_div=B0,rho_div=Bu);B2[0]+=I[CS];Bi+=I['mm'];Bj+=I['ae'];Bk+=I[EM]
				if I[EN]>A_:A_=I[EN]
			elif k==1:
				Q=G-B+1;j=(2*Q*M*D8>>32)+2+Q*Q;l=(B<<32)-D;d=p(l)>>13;m=d*d*F5[1]>>64;P=J+j-m;e=((2*B+1)*Q32-2*D)*D8>>32;f=-F7(F6[0][1:],Bw[2:])-D5[1][0]*B;R=Bx[1][B%3][X]
				for C in A(B,G+1):
					R=HL[R];D7[0]+=P<0
					if P>=0:Bw[1]=C;D9(0,P,az and C==0,[L[cg]+C*N[cg]for cg in A(12)],f,0,R&7,(F[0]+C*H[0]&K,F[1]+C*H[1]&K))
					P-=e;e+=HI;f-=D5[1][0]
			else:
				n=F5[k];T=Bx[k][B%3][X]
				for C in A(B,G+1):
					U=(C<<32)-D;g=-U if U<0 else U;o=g-M if g>M else 0;h=o>>13;V=h*h*n>>64;D7[k-1]+=V>J
					if V<=J:Bw[k]=C;D9(k-1,J-V,az and C==0,[L[ch]+C*N[ch]for ch in A(12)],E,T,E,(F[0]+C*H[0]&K,F[1]+C*H[1]&K))
					T=Bx[k][1][T]
		D9(11,Q38,Q,[0]*12,E,0,E,HK);Cv+=B2[0];Ez+=C(Bv);Cw.append(F(C(Bv)/N(1,B2[0]),3));F8=[2.**(log2(AM[12-aQ]/AM[12])+.5*(log2(I[aQ])-aQ*log2(Ac)))for aQ in A(1,12)];Cx.append(F(C(F8),3));Cy.append(F(BV([bU+bV for(bU,bV)in A0(Bv,D7)])/N(1,B2[0]),2));Cz.append(F(BV(F8),2));F9=[[J((B1[aR][bt]<<32)//I[bt+1])if aR>bt else 0 for bt in A(12)]for aR in A(12)];HM=[G(o(K,(I[aS+1]<<64)//(I[aS]*c)))for aS in A(12)];B3=N(log2(I[A+1])-log2(I[A])-log2(c)for A in A(12));f['n2']+=1;f['fail2']+=B3>=0;f[CT]=N(f[CT],B3);FA=[G(((I[aT]*c<<38)+I[aT+1]-1)//I[aT+1])for aT in A(12)]
		for HN in FA:Bp=N(Bp,HN.bit_length())
		DA=[0]*12
		for B in A(11,-1,-1):Ay=P.isqrt((I[B]*c+I[B+1]-1)//I[B+1])+1;DA[B]=Ay+C((DA[A]+1)//2 for A in A(B+1,12))
		By=[C(DA[A]for A in A(aU+1,12))+1 for aU in A(12)];HO=[I[aV+1]/I[aV]for aV in A(12)];HP=sqrt(a(c));DB=0;HQ=0;HR=0;HS=0;HT=0;HU=0;HV=0;HW=0;FB=0
		while DB<300:
			FC=[A4.gauss(s,O)for bu in A(11)];HX=C(A*A for A in FC);HY=HP*A4.random()**(O/11.)/sqrt(HX);B4=[0]*12;DC=Q38;FD=Q
			for B in A(11,0,-1):
				AF=0
				for b in A(B+1,12):AF=J(AF-J(F9[b][B]*B4[b]))
				Bs=N(Bs,p(AF).bit_length());B4[B]=n(F(AF/a(Q32)+FC[B-1]*HY/sqrt(HO[B])));B5=G(p(J((B4[B]<<32)-AF)));HZ=B5-By[B]if B5>By[B]else 0;B6=HZ>>13;FE=G(B6*B6);Br=N(Br,FE.bit_length());B7=Ar(FE,HM[B])
				if B7>DC:FD=A1;break
				DC-=B7
			if not FD:continue
			AN=0
			for b in A(1,12):AN=J(AN-J(F9[b][0]*B4[b]))
			Bz=G(Ar(DC,FA[0])+1);Bq=N(Bq,Bz.bit_length());B8=As(Bz);Ae=J(AN-By[0]-B8+Q32-1)>>32;B9=J(AN+By[0]+B8)>>32
			if Ae>B9:continue
			Ha=[C(B4[A]*A5[A][bv]for A in A(1,12))for bv in A(12)];w=Cp(Ha,A5[0],Ae,B9,T12_chk=c,verify_ref=Q,p_div=S,s0_div=B0,rho_div=Bu)
			if w[CS]==0:continue
			HQ+=1;DB+=w[CS];HR+=w['fast'];HS+=w['bnd'];HT+=w['segs'];HU+=w['bspl'];HV+=w['fspl'];HW+=w['cross'];FB+=w['clips'];D0+=w['mm'];D1+=w['ae'];D2+=w[EM]
		C_+=DB;E_+=FB;DD=[2.**(log2(AM[12-aW]/AM[12])+.5*(log2(I[aW])-aW*log2(c)))for aW in A(1,12)];Hb=DD[0];JI=C(DD[1:]);Hc=BV(DD);Hd=a(WHEEL_OPS[EO]+WHEEL_OPS[CR]);He=Hb*Hd;Hf=a(WHEEL_OPS[EP])+.015;Hg=Hc+He+Hf+M.C_LAT_APPLY;Bo.append(F(Hg,2))
	Hh=C(O/A for A in Bh);Hi=F(C(O/A*B for(A,B)in A0(Bh,Bo))/Hh,2);Hj=N(Bo);FF=F(1e2*N(p(A-B)/o(A,B)for(A,B)in A0(Cw,Cx)),2);DE=F(1e2*N(p(A-B)/o(A,B)for(A,B)in A0(Cy,Cz)),2);DF=BW;d=HT_K;FG=DF*d;DG=[];FH=0;B_=[];DH=[];FI=0;del PT_ORB[:];del PT_OLD[:];import pt_verify as Hk;FJ=[0];FK=HT_SEED_BASE;C0=P.sumprod
	for FL in A(DF):
		Hl,AO=ht_stratum(FL);AP=0;AQ=[];FM=[]
		for Hm in A(d):
			FN,C1,FO,DI,V,DJ,FP=ht_orbit(FK,FL,Hm,Hl,AO,q);AP+=FN;FI+=FN;PT_ORB.append((C1,[bT[:]for bT in DI],V,DJ,c));Af=[2.**(log2(AM[12-bw]/AM[12])+.5*(log2(V[bw])-bw*log2(c)))for bw in A(1,12)];Hn,JJ=Af[0],C(Af[1:]);Ho=[[(DJ[bx][bL]<<32)//V[bL+1]if bx>bL else 0 for bL in A(12)]for bx in A(12)];Hp=[o(K,(V[by+1]<<64)//(V[by]*c))for by in A(12)];B3=N(log2(V[A+1])-log2(V[A])-log2(c)for A in A(12));f['n']+=1;f['fail']+=B3>=0;f[CU]=N(f[CU],B3);Hq=[((V[bz]*c<<38)+V[bz+1]-1)//V[bz+1]for bz in A(12)];DK=[0]*12
			for B in A(11,-1,-1):Hr=P.isqrt((V[B]*c+V[B+1]-1)//V[B+1])+1;DK[B]=Hr+C((DK[A]+1)//2 for A in A(B+1,12))
			FQ=[C(DK[A]for A in A(bA+1,12))+1 for bA in A(12)];Hs=[V[bB+1]/V[bB]for bB in A(12)];Ht=sqrt(a(c));Hu=[O/sqrt(Hs[bC])for bC in A(12)];FR=[[Ho[bM][bD]for bM in A(bD+1,12)]for bD in A(12)];Hv=[[DI[bN][bE]for bN in A(1,12)]for bE in A(12)];DL=BA=0;BB=BC=C2=s;Hw=FP.gauss;Hx=FP.random;Hy=O/a(Q32);FS=DI[0];Hz=Hq[0];FT=FQ[0];H_=FR[0];Ag=WHEEL_OPS;DM=C1>=EMIT8_P_MIN;BD=EMIT8_FLUSH_AMORT if DM else .015;DN=a(Ag[EO]);DO=a(Ag['emit8_step1'if DM else EP]);FU=a(Ag['emit8_last'if DM else'emit9_last']);BE=a(Ag[CR]);I0=a(Ag['wheel_skip1_stub']);I1=a(Ag['wheel_skip2_stub']);I2=P.isqrt((V[0]*c+V[1]-1)//V[1])+1;BF=s if 2*I2+3<=128 else 2.
			while BA<160:
				DP=[Hw(s,O)for bO in A(11)];I3=Ht*Hx()**(O/11.)/sqrt(C0(DP,DP));C3=[0]*12;DQ=Q38;FV=Q
				for B in A(11,0,-1):
					AF=-C0(FR[B],C3[B+1:])if B<11 else 0;FW=n(F(AF*Hy+DP[B-1]*I3*Hu[B]));C3[B]=FW;B5=p((FW<<32)-AF);FX=FQ[B]
					if B5>FX:
						B6=B5-FX>>13;B7=B6*B6*Hp[B]>>64
						if B7>DQ:FV=A1;break
						DQ-=B7
				if not FV:continue
				AN=-C0(H_,C3[1:]);Bz=(DQ*Hz>>64)+1;B8=P.isqrt(Bz<<52)+1;Ae=AN-FT-B8+Q32-1>>32;B9=AN+FT+B8>>32
				if Ae>B9:continue
				I4=C3[1:];DR=[C0(I4,Hv[bP])for bP in A(12)];C4=(1 if _bad_f3(DR)else 0)|(2 if _bad_f3([DR[bX]+FS[bX]for bX in A(12)])else 0)|(4 if _bad_f3([DR[bW]+2*FS[bW]for bW in A(12)])else 0);Ah=B9-Ae+1;DL+=1;BA+=Ah
				if C4==7:BB+=DN+BE+BF
				elif C4==0:DS=(Ah+127)//128;BB+=DN+BE+BF+(DS-1)*(BE+BF+2e1);BC+=(Ah-1)*(DO+BD)+(FU+BD);C2+=Ah*M.C_LAT_APPLY
				else:
					BB+=DN+BE+BF;C5=(C4|C4<<3)>>Ae%3&7;I5=C5 in(3,5,6);I6,I7=divmod(Ah,3);C6=I6*(1 if I5 else 2)+C(1 for A in A(I7)if not C5>>A&1);DS=(Ah+127)//128;BG=s
					if C6>=1:
						if C5&1==0:BG=s
						elif C5&2==0:BG=I0
						else:BG=I1
					BC+=BG+(DS-1)*(BE+BF+BG+2e1)
					if C6>=1:BC+=(C6-1)*(DO+BD)+(FU+BD);C2+=C6*M.C_LAT_APPLY
			FH+=BA;C7=BV(Af)+Hn*(BB/DL);PT_OLD.append(C7+(BC+C2)/BA)
			if Hk.pt_cert(V,DJ,c)[0]:FY=[Af[0]*HT_NEW['PT_ROWS'],Af[1]*HT_NEW['PT_CANDS']]+Af[2:];C7=Gd(FY)+FY[0]*(BB/DL)+HT_NEW['PT_EXTRA']*(DO+BD+M.C_LAT_APPLY)
			else:FJ[0]+=1
			I8=(BC+C2)/BA;FZ=C7+I8;DG.append(F(FZ,2));AQ.append((C1,FO,FZ));FM.append((C1,FO,C7))
		B_.append((AO,AP,AQ));DH.append((AO,AP,FM))
	Fa=BH=Fb=s
	for((AO,AP,AQ),(_,_,DT))in A0(B_,DH):A6=AO*(d-O)/(AP-O);Fa+=A6*C(B*C/A for(A,B,C)in AQ)/d;Fb+=A6*C(B*C/A for(A,B,C)in DT)/d;BH+=A6*C(B/A for(A,B,C)in AQ)/d
	BI=Fa/BH;C8=Fb/BH;Fc=Fd=s
	for((AO,AP,AQ),(_,_,DT))in A0(B_,DH):A6=AO*(d-O)/(AP-O);Fe=[aY/aX*(aZ-BI)for(aX,aY,aZ)in AQ];DU=C(Fe)/d;I9=C((A-DU)**2 for A in Fe)/(d-1);C9=A6*A6/a(d-2);Fc+=(A6*A6+C9)*(I9/d)+C9*(DU*DU);Ff=[a1/a0*(a2-C8)for(a0,a1,a2)in DT];DV=C(Ff)/d;IA=C((A-DV)**2 for A in Ff)/(d-1);Fd+=(A6*A6+C9)*(IA/d)+C9*(DV*DV)
	DW=sqrt(Fc)/BH;Fg=sqrt(Fd)/BH;IB=F(BI+3.*DW,3);DX=O+N(.03,DE/1e2);IC=F(BI+(DX-O)*(C8+3.*Fg)+DX*3.*DW,3);ID=BZ['GA'];k=[ID[ab][:]+[O if ab==bF else s for bF in A(6)]for ab in A(6)]
	for x in A(6):
		Fh=N(A(x,6),key=lambda r:p(k[r][x]));k[x],k[Fh]=k[Fh],k[x];IE=k[x][x];k[x]=[a3/IE for a3 in k[x]]
		for BJ in A(6):
			if BJ!=x and k[BJ][x]:IF=k[BJ][x];k[BJ]=[k[BJ][bQ]-IF*k[x][bQ]for bQ in A(12)]
	Fi=[sqrt(k[ac][6+ac])for ac in A(6)];Fj=[F(572.*sqrt(O+M.SC_GUARD_DELTA)*ad,1)for ad in Fi];Fk=[F(572.*sqrt(1.764*(O+M.SC_GUARD_DELTA))*ae,1)for ae in Fi];Fl=all(Fj[A]<BT[A]+O for A in A(6))and all(Fk[A]<BU[A]+O for A in A(6))and PROD_MB_MR<1<<59;Y=1024;IG=16;IH=16;JK=6;AR=[]
	for II in A(0,4):
		for IJ in A(-4,5):
			for IK in A(-4,5):
				for IL in A(-2,3):
					CA=IL,0,0,II,IJ,IK
					if h(CA):IM=C(CA[B]*q[B][C]*CA[C]for B in A(6)for C in A(6));AR.append((IM,CA))
					if D(AR)>=Y:break
				if D(AR)>=Y:break
			if D(AR)>=Y:break
		if D(AR)>=Y:break
	AR.sort(key=lambda x:x[0]);DY=[af for(ag,af)in AR];Fm=[Gg(ah,EJ,Gf)for ah in DY]
	def Fn(pa_v):
		A=pa_v;B=A>>50&511;C=(A>>40&1023)-508;D=(A>>30&1023)-508;E=(A>>20&1023)-510
		if A>>10&1023!=512 or A&1023!=512 or not(0<=B<4 and 0<=C<9 and 0<=D<9 and 0<=E<5):return-1
		return B<<11|C<<7|D<<3|E
	CB=[-1]*8192
	for(IN,IO)in BX(Fm,1):DZ=Fn(IO);assert DZ>=0 and CB[DZ]==-1;CB[DZ]=IN
	def IP(av):A=av;B=A[0]&65535|(A[1]&65535)<<16|(A[2]&65535)<<32|(A[3]&65535)<<48;C=A[4]&65535|(A[5]&65535)<<16;return G(B),G(C)
	def Ai(v):v&=65535;return v-65536 if v>=32768 else v
	Fo=[E]+[IP(a4)for a4 in DY]
	def Fp(w0,w1,c_k,p_val,bias_val,inv_p):C=p_val;B=c_k;F=Ai(w0);H=Ai(w0>>16);I=Ai(w0>>32);K=Ai(w0>>48);L=Ai(w1);M=Ai(w1>>16);A=J(F*B[0]);A=J(A+J(H*B[1]));A=J(A+J(I*B[2]));A=J(A+J(K*B[3]));A=J(A+J(L*B[4]));A=J(A+J(M*B[5]));E=G(A+bias_val);N=Ar(E,inv_p);D=G(E-G(N*C));assert D<2*C;return D-C if D>=C else D
	Da=0;Db=0;CC=0
	for IQ in(42,46,48):
		Z=(1<<IQ)+1
		while Q:
			if M._isp(Z):
				Dc=_find_root_h(Z)
				if Dc is not E:
					IR=_eval_poly(M.W_COEFFS,Dc,Z);Fq=-((M.U_BN254+1)%Z*IR+2)%Z
					if Fq!=0:
						Fr=_sqrt_mod(H(Fq,-1,Z),Z)
						if Fr>0:break
			Z+=2
		IS=H(3,-1,Z);Dd=[1]
		for _ in A(5):Dd.append(Dd[-1]*Dc%Z)
		Fs=[C(M.ZK_NUM[a5][A]*Dd[A]for A in A(6))*IS%Z for a5 in A(6)];Ft=-H(Fr,-1,Z)%Z;IT=[a6*Ft%Z for a6 in Fs];Fu=G(3072*Z);IU=G((1<<64)//Z);CC=N(CC,(2*Fu).bit_length())
		for L in A(Y,0,-1):
			De,Df=Fo[L];IV=Fp(De,Df,IT,Z,Fu,IU);u=DY[L-1];IW=C(u[A]*Fs[A]for A in A(6))*Ft%Z
			if IV!=IW:Da+=1
			Db+=1
	W=127
	while Q:
		if M._isp(W):
			Dg=_find_root_h(W)
			if Dg is not E:
				IX=_eval_poly(M.W_COEFFS,Dg,W);Fv=-((M.U_BN254+1)%W*IX+2)%W
				if Fv!=0:
					Fw=_sqrt_mod(H(Fv,-1,W),W)
					if Fw>0:break
		W+=2
	IY=H(3,-1,W);Dh=[1]
	for _ in A(5):Dh.append(Dh[-1]*Dg%W)
	IZ=[C(M.ZK_NUM[ai][A]*Dh[A]for A in A(6))*IY%W for ai in A(6)];Ia=-H(Fw,-1,W)%W;Ib=[aj*Ia%W for aj in IZ];Ic=G(3072*W);Id=G((1<<64)//W);T=16;BK=(W+T-1)//T;A7=[0]*(Y+1);Di=[0]*BK;Fx=[0]*(Y+1);Ie=(1<<47)-1
	for L in A(Y,0,-1):De,Df=Fo[L];B=Fp(De,Df,Ib,W,Ic,Id);A7[L]=B;e=B//T;Fx[L]=Di[e];Di[e]=(B&T-1)<<47|L
	Dj=0
	for e in A(BK):
		BL=A4.randint(40,Y);Fy=set();CD=Di[e]
		while CD:
			L=CD&Ie
			if L>BL:break
			Fy.add((CD>>47,L));CD=Fx[L]
		If={(A7[a7]&T-1,a7)for a7 in A(1,BL+1)if A7[a7]//T==e}
		if Fy!=If:Dj+=1
	A8=[0]*(BK+1)
	for L in A(1,Y+1):A8[A7[L]//T+1]+=1
	for e in A(BK):A8[e+1]+=A8[e]
	Aj=[0]*Y;Ig=A8[-1]
	for L in A(Y,0,-1):B=A7[L];e=B//T;A8[e+1]-=1;Aj[A8[e+1]]=L<<17|B&T-1
	A8=A8[1:]+[Ig];Dk=0;AG=2048
	for e in A(BK):
		BL=A4.randint(40,Y);Fz=BL<<17|15
		for CE in(1,-1):
			Dl=[0]*(T*AG);Dm=[0]*(T*AG);Ih,F_=A8[e],A8[e+1];AD=Ih;CF=0 if CE>0 else(T-1)*AG
			while AD+4<=F_:
				Ii,Ij,Ik,G0=Aj[AD],Aj[AD+1],Aj[AD+2],Aj[AD+3]
				if G0>Fz:break
				for Ak in(Ii,Ij,Ik,G0):A9=Ak&15;Dn=CF+(A9<<11)if CE>0 else CF-(A9<<11);At(Dl,Dn+(Ak>>16>>1),720)
				AD+=4
			while AD<F_:
				Ak=Aj[AD];AD+=1
				if Ak>Fz:break
				A9=Ak&15;Dn=CF+(A9<<11)if CE>0 else CF-(A9<<11);At(Dl,Dn+(Ak>>16>>1),720)
			for L in A(1,BL+1):
				if A7[L]//T==e:
					A9=A7[L]&T-1
					if CE<0:A9=T-1-A9
					Dm[A9*AG+L]=Dm[A9*AG+L]+720&65535
			if Dl!=Dm:Dk+=1
	l=IG;Al=4;AA=0
	for G1 in(0,1,(1<<34)+7,(1<<35)-1):
		if G(G1<<15)>=1<<50:AA+=1
		if G(G1<<14)>=1<<49:AA+=1
	for CG in(0,1,32768,(1<<28)-1):
		Do=G((1<<35)-1<<29|CG)
		if Do>>29!=(1<<35)-1 or Do<<35&K!=G(CG<<35):AA+=1
		if G((Do<<35&K)+((1<<49)-2)&K)!=(CG<<35)+((1<<49)-2):AA+=1
		if G((CG<<35)+((1<<49)-2))>=1<<64:AA+=1
	Dp=0;AS=l*IH*8;BM=[0]*(AS+1);Am=[1]+[0]*(l-1)
	for BN in A(AS):
		if BN&7==0:Il=o(A(l),key=lambda s:Am[s]);Dq=Il
		G2=Am[Dq];Am[Dq]=G2+1;BM[BN]=G(Dq<<29|G2)
	BM[AS]=0
	if h(BM[A]==0 for A in A(AS)):Dp+=1
	G3=[[[0]*(Y+1)for a8 in A(Am[ak])]for ak in A(l)];Dr=[[[0]*(Y+1)for a9 in A(Am[al])]for al in A(l)];Ds=2046;Im=(1<<35)-1;AT=[16]*l;Dt=[CQ(16384)for am in A(l)];BO=0;CH=0;BP=0;G4=BS.Struct('<Q').pack_into;G5=BS.Struct('<Q').unpack_from;In=BS.Struct('<QQ').pack_into;Io=BS.Struct('<QQ').unpack_from
	def G6(sid_v,w_lp_v,cnt_v):
		C=sid_v;nonlocal BP;B=Dt[C];In(B,0,cnt_v|w_lp_v<<16,0);D,H=Io(B,0);E=D&65535;I=D>>16&65535;assert 1<=E<=Ds;BP+=1;F=16;J=G3[C]
		for H in A(E):G=G5(B,F)[0];F+=8;At(J[G>>35],(G&Im)>>1,I)
	G7=M.rand_below;G8=A4.getrandbits
	for G9 in A(224):
		Du=752+G9;Ip=37000 if G9 in(0,223)else 18
		for Iq in A(Ip):
			L=1+G7(G8,Y);BN=AS if Iq%500==0 else G7(G8,AS);Ir=Fm[L-1];Dv=BM[BN];r=Dv>>29;GA=Dv&(1<<29)-1
			if BN<AS:Dr[r][GA][L]=Dr[r][GA][L]+Du&65535
			Is=CB[Fn(Ir)];GB=G((Dv<<35&K)+(Is<<1))
			if GB>=1<<64:AA+=1
			G4(Dt[r],AT[r],GB);AT[r]+=8
			if AT[r]&16384:
				if AT[r]!=16384:AA+=1
				G6(r,Du,Ds);AT[r]=16;BO+=1
		for r in A(l):
			GC=AT[r]-16>>3
			if GC:G6(r,Du,GC);AT[r]=16;CH+=1
	Dt=E;CB=E;BM=E;CI=0
	for CJ in A(l):
		It=1 if CJ==0 else 0
		for GD in A(It,Am[CJ]):
			if G3[CJ][GD]!=Dr[CJ][GD]:CI+=1
	GE=CI==0 and AA==0 and BP==BO+CH;Iu=l;GF=[CQ(16384)for an in A(Al)];GG=[CQ(16384)for ao in A(Al)];Dw=[0]*Al
	for An in A(1,Y+1):AH=A7[An]//T%Al;GH=Dw[AH];G4(GF[AH],GH*8,An<<17|A7[An]&T-1);Dw[AH]=GH+1
	Dx=0;GI=[0]*(T*AG);Dy=[0]*(T*AG)
	for An in A(1,Y+1):GJ=(A7[An]&T-1)*AG+An;Dy[GJ]=Dy[GJ]+720&65535
	for AH in A(Al):
		GG[AH][:]=GF[AH]
		for Iv in A(Dw[AH]):GK=G5(GG[AH],Iv*8)[0];At(GI,((GK&15)<<11)+(GK>>16>>1),720)
	if GI!=Dy:Dx+=1
	GL=256;Iw=l;Dz=0;D_=[0]*l;E0=[A4.randint(1,64)for ap in A(GL)]
	for Ix in E0:Iy=o(A(l),key=lambda s:D_[s]);D_[Iy]+=Ix
	if N(D_)>C(E0)/a(l)+N(E0)+1e-09:Dz+=1
	CK=BZ['GB'];Iz=CK[0][0];CL=[[CK[aq][a_]-CK[aq][0]*CK[0][a_]/Iz for a_ in A(1,6)]for aq in A(1,6)];I_=sqrt(N(C((.5 if B>>C&1 else-.5)*CL[C][D]*(.5 if B>>D&1 else-.5)for C in A(5)for D in A(5))for B in A(32)));J0=8.*P.pi*P.pi/15./sqrt(M._det_mat(CL));m=[CL[ar][:]+[O if ar==bG else s for bG in A(5)]for ar in A(5)]
	for y in A(5):
		GM=N(A(y,5),key=lambda r:p(m[r][y]));m[y],m[GM]=m[GM],m[y];J1=m[y][y];m[y]=[ba/J1 for ba in m[y]]
		for BQ in A(5):
			if BQ!=y and m[BQ][y]:J2=m[BQ][y];m[BQ]=[m[BQ][bR]-J2*m[y][bR]for bR in A(10)]
	J3=[sqrt(m[at][5+at])for at in A(5)];E1=0
	for CM in(4.,8.,12.):
		GN=[n(floor(CM*J3[bb]))for bb in A(5)];J4=CM*CM;GO=[0]
		def GP(ix,v5):
			B=ix
			if B==5:
				if C(v5[B]*CL[B][C]*v5[C]for B in A(5)for C in A(5))<=J4+1e-12:GO[0]+=1
				return
			for D in A(-GN[B],GN[B]+1):v5[B]=D;GP(B+1,v5)
		GP(0,[0]*5)
		if GO[0]>J0*(CM+I_)**5:E1+=1
	def E2(a,b):
		C=[0]*11
		for D in A(6):
			for B in A(6):C[D+B]+=a[D]*b[B]
		for E in A(10,5,-1):
			F=C[E]
			if F:
				for B in A(6):C[E-6+B]-=F*M.H_COEFFS[B]
		return C[:6]
	def J5(M6):
		B=[b2[:]for b2 in M6];D=1;H=1
		for C in A(5):
			if B[C][C]==0:
				E=next((A for A in A(C+1,6)if B[A][C]!=0),-1)
				if E<0:return 0
				B[C],B[E]=B[E],B[C];D=-D
			I=B[C][C]
			for F in A(C+1,6):
				for G in A(C+1,6):B[F][G]=(B[F][G]*I-B[F][C]*B[C][G])//H
			H=I
		return D*B[5][5]
	def GQ(a,b,u_val=M.U_BN254):
		D=[C(a[A]*M.ZK_NUM[A][b3]for A in A(6))for b3 in A(6)];E=[C(b[A]*M.ZK_NUM[A][b4]for A in A(6))for b4 in A(6)];F=E2(D,D);L=[(u_val+1)*b5+2*b6+b7 for(b5,b6,b7)in A0(E2(M.W_COEFFS,F),F,E2(E,E))];G=[];B=L[:]
		for N in A(6):
			G.append(B)
			if N<5:
				H=[0]+B[:5];I=B[5]
				if I:
					for J in A(6):H[J]-=I*M.H_COEFFS[J]
				B=H
		K=p(J5([[G[b9][b8]for b9 in A(6)]for b8 in A(6)]));assert K%9**6==0;return K//9**6
	def GR(P0,xs):
		B=[BY(xs)]
		while D(B[-1])>1:C=B[-1];B.append([C[cb]*C[cb+1]if cb+1<D(C)else C[cb]for cb in A(0,D(C),2)])
		G=[P0%B[-1][0]]
		for J in A(D(B)-2,-1,-1):G=[G[b_>>1]%ca for(b_,ca)in BX(B[J])]
		H=[]
		for(I,K)in A0(xs,G):
			E=I;F=P.gcd(E,K)
			while F>1:E//=F;F=P.gcd(E,F)
			H.append((E,I//E))
		return H
	Ao=[]
	while D(Ao)<64:
		u=[A4.randint(-220,220)for bc in A(6)];AU=[A4.randint(-300,300)for bd in A(6)]
		if h(u)and h(AU)and M.is_ideal_primitive_oh(u,AU):Ao.append(GQ(u,AU))
	GS=M._primes_upto(1<<14);GT=1
	for Ap in GS:GT*=Ap
	GU=[]
	for BR in A(0,64,16):GU.extend(GR(GT,Ao[BR:BR+16]))
	E3=0
	for(GV,(J6,J7))in BX(GU):
		CN=Ao[GV]
		for Ap in GS:
			while CN%Ap==0:CN//=Ap
		if CN!=J6 or Ao[GV]//CN!=J7:E3+=1
	def J8(n):
		if n<2:return A1
		for D in(2,3,5,7,11,13,17,19,23,29,31,37):
			if n==D:return Q
			if n%D==0:return A1
		C=n-1;E=0
		while C&1==0:C>>=1;E+=1
		for F in(2,3,5,7,11,13,17,19,23,29,31,37):
			if F>=n:continue
			B=H(F,C,n)
			if B in(1,n-1):continue
			for G in A(E-1):
				B=B*B%n
				if B==n-1:break
			else:return A1
		return Q
	def J9(n):
		if n&1==0:return 2
		for E in A(1,50):
			C=2;B=2;D=1
			while D==1:C=(C*C+E)%n;B=(B*B+E)%n;B=(B*B+E)%n;D=P.gcd(p(C-B),n)
			if D!=n:return D
		return n
	def JA(n,small_ps):
		C={}
		for A in small_ps:
			if A*A>n:break
			while n%A==0:C[A]=C.get(A,0)+1;n//=A
		D=[n]if n>1 else[]
		while D:
			B=D.pop()
			if B==1:continue
			if J8(B):C[B]=C.get(B,0)+1
			else:E=J9(B);D.append(E);D.append(B//E)
		return C
	E4=256;E5=1<<20;JB=M._primes_upto(1024);JC=M._primes_upto(E4);GW=1
	for Ap in JC:GW*=Ap
	Aq=[]
	while D(Aq)<320:
		u=[A4.randint(-4,4)for be in A(6)];AU=[A4.randint(-6,6)for bf in A(6)]
		if h(u)and h(AU)and M.is_ideal_primitive_oh(u,AU):
			CO=GQ(u,AU,u_val=5)
			if CO>=1<<25:Aq.append(CO)
	GX=[]
	for BR in A(0,D(Aq),32):GX.extend(GR(GW,Aq[BR:BR+32]))
	E6=0;GY=0;E7=0;GZ=0;E8=0;E9=0;EA=1e9;EB=1e9
	for(CO,(EC,_))in A0(Aq,GX):
		JD=log2(CO);CP=JA(EC,JB);JE=C(floor(16.*log2(A))for A in CP if E4<=A<=E5);JF=h(B>=2 for(A,B)in CP.items()if E4<=A<=E5);JG=floor(ceil(JD)/8.);Ab=(JG+1)/16.;ED=(log2(EC)if EC>1 else s)-JE/16.;EE=[bg for bg in CP if bg>E5];Ga=ED<=Ab+.25
		if not EE:
			if not JF:
				E6+=1;Gb=Ab+.25-ED
				if Gb<EA:EA=Gb
				if not Ga:E8+=1
			else:GY+=1
		else:
			if D(EE)==1 and CP[EE[0]]==1:E7+=1
			else:GZ+=1
			Gc=ED-(Ab+.25)
			if Gc<EB:EB=Gc
			if Ga:E9+=1
	JH=ES and EZ and Cb==0 and Ba<=.0012 and EX==0 and EY==0 and Cd==0 and Dk==0 and Cm==0 and Bg==0 and Da==0 and Db>=3000 and CC<64 and Dj==0 and CI==0 and AA==0 and Dp==0 and E1==0 and Dx==0 and GE and Dz==0 and BO>=20 and BP==BO+CH and Bt==0 and Fl and Cv>=10000 and Bi==0 and Bj==0 and Bk==0 and A_<=4 and Bl==0 and Bm>=10000 and Bn==0 and FF<=5. and DE<=5. and C_>=2000 and D0==0 and D1==0 and D2==0 and Bp<64 and Bq<40 and Br<64 and Bs<63 and E3==0 and E6>=10 and E7>=100 and E8==0 and E9==0;return AI(ok=JH,smat_cyclic_a345_ok=ES,zero_float_fp_enum=Q,strict_u64_i64_umulh_checked=Q,zero_exception_64bit_q38_proved=EZ,cert_min_gsym_exact=CV,cert_min_gsym_pts_le_5gsc=Gi,uncond_min_bk_star_bits=F(log2(a(Au)),3),uncond_inv_max_bits=CW.bit_length(),uncond_r2_q12_max_bits=CX.bit_length(),uncond_ck_max_bits=CY.bit_length(),uncond_d2_max_bits=EW.bit_length(),isqrt_u40_lut7_errs=Cb,lut4_buckets=D(t),lut4_bucket_violations=EX,lut4_buckets_above_lut7=EY,lut4_max_overshoot_pct=F(Gn*1e2,4),isqrt_u40_lut7_max_over_pct=F(Ba*1e2,4),staddh_16b_semantic_errs=Cd,b2_x4_prefix_staddh_errs=Dk,orbits_6ideal_checked=D(Aw),lattice_b2_mismatches=Cm,lattice_b2_raw_hits=Em,lattice_b2_jvor_hits=En,complete_root_to_leaf_orbits=D(Bh),complete_root_to_leaf_leaves=Cv,complete_root_to_leaf_nodes=Ez,complete_root_to_leaf_mismatches=Bi,complete_root_to_leaf_apply_errs=Bj,complete_root_to_leaf_div_errs=Bk,complete_root_to_leaf_max_sa_cones=A_,wheel_tree_leaves_checked=Bm,wheel_tree_state_errs=Bl,swar21_cur_errs=Bn,complete_scaled_nodes_per_leaf_obs=Cw,complete_scaled_nodes_per_leaf_pred=Cx,complete_scaled_max_rel_err_pct=FF,complete_scaled_fp_ir_obs=Cy,complete_scaled_fp_ir_pred=Cz,complete_scaled_fp_ir_max_rel_err_pct=DE,prod_exec_leaves=C_,prod_exec_mismatches=D0,prod_exec_apply_errs=D1,prod_exec_div_errs=D2,prod_exec_endpoint_clips=E_,prod_min_bk_star_bits=F(D3,2),prod_max_inv_b2_bits=Bp,prod_max_r2_q12_bits=Bq,prod_max_d2_bits=Br,prod_max_ck_bits=Bs,prod_12lvl_analytic_ir_per_leaf=Bo,prod_inv_bit_weighted_mean_ir=Hi,prod_max_ir=Hj,prod_256orb_strata=DF,prod_256orb_draws_per_stratum=d,prod_256orb_count=FG,prod_256orb_leaves_sampled=FH,prod_256orb_min_ir=o(DG),prod_256orb_max_ir=N(DG),prod_256orb_wmean_ht_ir=F(BI,3),prod_256orb_se_strat_ir=F(DW,3),prod_256orb_gauss_wmean_ht_ir=F(C8,4),prod_256orb_gauss_se_strat_ir=F(Fg,4),prod_256orb_exact_wmean_ht_ir=F(BI-C8,4),prod_256orb_wmean_plus_3se_ir=IB,prod_256orb_extrap_factor=F(DX,4),prod_256orb_extrap_upper_bound_ir=IC,prod_256orb_exact_rejection_ht=Q,prod_256orb_uniform_random_root_choice=Q,prod_256orb_pt_fallback_orbits=FJ[0],prod_256orb_run=FK,prod_256orb_strata_ht=B_,beta_ge1_census=AI(orbits_256=f['n'],ge1_256=f['fail'],max_log2_beta_256=F(f[CU],3),orbits_prod=f['n2'],ge1_prod=f['fail2'],max_log2_beta_prod=F(f[CT],3)),prod_256orb_avg_rej_tries=F(FI/a(FG),1),prod_budget_c_lat_leaf=M.C_LAT_LEAF,swar_tie_cases_checked=D(Cq),swar_tie_canon_errs=Bg,barrett_kta_18op_tested=Db,barrett_kta_18op_errs=Da,barrett_max_sp_bits=CC,desc_ai_prefix_chain_prime=W,desc_ai_prefix_chain_errs=Dj,score_node_replay_errs=CI,fabric_1hop_addr_errs=AA,b_tab_packing_errs=Dp,discrete_chord_bound_test_errs=E1,fabric_1hop_list_node_errs=Dx,fabric_1hop_per_node_bytearray_zero_dict=Q,fabric_1hop_zero_frag_exact_ok=GE,score_nodes_sim=Iu,greedy_bound_subbands=GL,greedy_bound_nodes=Iw,greedy_node_bound_errs=Dz,fabric_1hop_list_nodes_sim=Al,logp_bands_1_16b=224,lat_8b_direct_ptr_page_cap=Ds,emit_full_pages_flushed=BO,emit_tail_pages_flushed=CH,emit_pages_streamed=BP,zero_hash_table_apply=Q,inc_c0_exact_errs=Bt,swar_bounds_ok=Fl,max_a_prod_coords=Fj,max_b_prod_coords=Fk,srt_p24_genuine_norms_checked=D(Ao),srt_batch_mismatches=E3,e2e_factored_c6_norms_checked=D(Aq),e2e_factored_smooth_sqfree=E6,e2e_factored_smooth_p2=GY,e2e_factored_1lp=E7,e2e_factored_multi_lp=GZ,gate_0lp_fn=E8,gate_0lp_fp=E9,gate_0lp_min_smooth_margin_b=F(EA,4),gate_0lp_min_1lp_gap_b=F(EB,2))
import math,random,struct
FP_C0=21.
FP_C1=82.
FP_C2=99.
FP_CL=[20,72,97,79,77,75,73,71,69,67,65,97]
PHI_F3=(1,0,2,2,0,1),(1,2,0,2,1,2),(1,2,2,1,2,2)
WHEEL_OPS=dict(chord_enter=0,emit9_step1=9,emit9_step2=9,emit9_step3=9,emit9_last=7,emit8_step1=8,emit8_step2=8,emit8_step3=8,emit8_last=6,emit8_epoch_flush=10,wheel_disp4=4,wheel_skip1_stub=2,wheel_skip2_stub=2)
EMIT8_P_MIN=6369051672525773
EMIT8_FLUSH_AMORT=.165
AMAX_MR=405,377,352,269,269,269
BMAX_MR=537,500,468,357,357,357
MA_MR=tuple(2*A+1 for A in AMAX_MR)
MB_MR=tuple(2*A+1 for A in BMAX_MR)
WA_MR=[1]
for _m in MA_MR[:-1]:WA_MR.append(WA_MR[-1]*_m)
WB_MR=[1]
for _m in MB_MR[:-1]:WB_MR.append(WB_MR[-1]*_m)
PROD_MA_MR=WA_MR[-1]*MA_MR[-1]
PROD_MB_MR=WB_MR[-1]*MB_MR[-1]
A_OFF_MR=sum(AMAX_MR[A]*WA_MR[A]for A in range(6))
B_OFF_MR=sum(BMAX_MR[A]*WB_MR[A]for A in range(6))
def _mr_da(a):return sum(a[A]*WA_MR[A]for A in range(6))&(1<<64)-1
def _mr_db5(b):return sum(b[A]*WB_MR[A]for A in range(6))<<5&(1<<64)-1
def _mr_pack_a(a):return A_OFF_MR+sum(a[A]*WA_MR[A]for A in range(6))&(1<<64)-1
def _mr_pack_b(b):return B_OFF_MR+sum(b[A]*WB_MR[A]for A in range(6))
HT_SEED_BASE=1
def _st_f3(v):A,G,H,B,I,C,D,J,K,E,L,F=v;M=(A+2*H+2*B+C)%3;N=(A+2*G+2*B+I+2*C)%3;O=(A+2*G+2*H+B+2*I+2*C)%3;P=(D+2*K+2*E+F)%3;Q=(D+2*J+2*E+L+2*F)%3;R=(D+2*J+2*K+E+2*L+2*F)%3;return M+3*(N+3*(O+3*(P+3*(Q+3*R))))
def _add_st_f3(s1,s2):
	A=0;B=1
	for C in range(6):A+=(s1%3+s2%3)%3*B;s1//=3;s2//=3;B*=3
	return A
def _bad_st_f3(st):A=st;B=A%3;A//=3;C=A%3;A//=3;D=A%3;A//=3;E=A%3;A//=3;F=A%3;G=A//3;return B==0 and E==0 or C==0 and F==0 or D==0 and G==0
def _bad_f3(v):A,G,H,B,I,C,D,J,K,E,L,F=v;return(A+2*H+2*B+C)%3==0 and(D+2*K+2*E+F)%3==0 or(A+2*G+2*B+I+2*C)%3==0 and(D+2*J+2*E+L+2*F)%3==0 or(A+2*G+2*H+B+2*I+2*C)%3==0 and(D+2*J+2*K+E+2*L+2*F)%3==0
def wheel_prices_from_kernels(path):B=count_kernel_segments(path);return{A:B.get(A,-1)for A in WHEEL_OPS}
def _fpw(v,o=0):return sum(FP_CL[A+o]*B for(A,B)in enumerate(v))
from fractions import Fraction as Frac
from math import sqrt,ceil,floor,log2
import model_bn254_v6 as M
import descent_bn254_v6 as DESC
M64=(1<<64)-1
I64_MIN=-(1<<63)
I64_MAX=(1<<63)-1
Q32=1<<32
Q38=1<<38
def _chk_u64(x):
	if x<0 or x>M64:raise OverflowError(f"u64 overflow: {x}")
	return x
def _chk_i64(x):
	if x<I64_MIN or x>I64_MAX:raise OverflowError(f"i64 overflow: {x}")
	return x
def count_kernel_segments(path):
	import re;F=re.compile('^\\s+[a-z][a-z0-9.]*\\b');C={};A=None
	for E in open(path,encoding='ascii'):
		D=E.strip()
		if D.startswith('//@seg'):A=D.split()[1];C.setdefault(A,0);continue
		if D.startswith('//@end'):A=None;continue
		B=E.split('//',1)[0].rstrip()
		if A is None or not B.strip()or B.endswith(':')or B.strip().startswith('.'):continue
		if F.match(B):C[A]+=1
	return C
def fp_prices_from_kernels(path):
	A=count_kernel_segments(path);C=[A['k0_body']+A['k0_back'],A['top_2']+A['body_2']+A['k1_exit']]
	for B in range(2,11):C.append(A['top_%d'%(B+1)]+A['body_%d'%(B+1)]+A['top_%d'%B])
	return C,A
def _umulh64(a,b):_chk_u64(a);_chk_u64(b);return a*b>>64
_SQRT_LUT=[0]*512+[ceil(sqrt(A+1)*2.**52.5)+(1<<31)for A in range(512,2048)]
def _isqrt_u40_lut7_q32(r2):
	A=r2
	if A<1 or A>=1<<40:raise OverflowError(f"u40 range error: {A}")
	C=64-A.bit_length();B=C&~1;D=A<<B>>53;return _chk_u64(_SQRT_LUT[D]>>(B>>1))
def _lut4_buckets():
	C=[]
	for A in range(24,65):
		if A==64:C.append((64<<11,0,0));continue
		for B in range(1024,2048):
			if A<=53:D,E=B<<53-A,(B+1<<53-A)-1
			else:
				if B&(1<<A-53)-1:continue
				D=E=B>>A-53
			C.append((A<<11|B,D,E))
	return C
_PB=dict(n=0,fail=0,maxb=-999.,n2=0,fail2=0,maxb2=-999.)
_LUT4=[0]*((64<<11)+1)
for(_i,_lo,_hi)in _lut4_buckets():_LUT4[_i]=math.isqrt(_hi+1<<52)+1
def _isqrt_lut4(r2):
	if r2<1 or r2>=(1<<40)+1:raise OverflowError(f"u40 range error: {r2}")
	A=r2-1;B=64-A.bit_length();C=A<<(B&63)&M64;return _LUT4[B<<11|C>>53]
def _staddh_u16(mem,off,w_lp):mem[off]=mem[off]+(w_lp&65535)&65535
_SP=math.sumprod
def _mv6(Mat,v):A=Mat;return _SP(A[0],v),_SP(A[1],v),_SP(A[2],v),_SP(A[3],v),_SP(A[4],v),_SP(A[5],v)
def _qf6(G,v):return _SP(v,(_SP(G[0],v),_SP(G[1],v),_SP(G[2],v),_SP(G[3],v),_SP(G[4],v),_SP(G[5],v)))
def _eval_poly(cs,s,p):
	A=0
	for B in reversed(cs):A=(A*s+B)%p
	return A
def _eval_zk(v,s,p):
	H=pow(3,-1,p);C=s*s%p;D=C*s%p;E=D*s%p;I=E*s%p;B=0;J=M.ZK_NUM
	for F in range(6):
		G=v[F]
		if G:A=J[F];K=(A[0]+A[1]*s+A[2]*C+A[3]*D+A[4]*E+A[5]*I)%p;B=(B+G*K)%p
	return B*H%p
def _sqrt_mod(n,p):
	n%=p
	if n==0:return 0
	if pow(n,(p-1)//2,p)!=1:return-1
	if p%4==3:A=pow(n,(p+1)//4,p);return min(A,p-A)
	B,I=p-1,0
	while B&1==0:B>>=1;I+=1
	E=2
	while pow(E,(p-1)//2,p)!=p-1:E+=1
	J,F,C,A=I,pow(E,B,p),pow(n,B,p),pow(n,(B+1)//2,p)
	while C!=1:
		D,G=C,0
		while D!=1:D=D*D%p;G+=1
		H=pow(F,1<<J-G-1,p);J=G;F=H*H%p;C=C*F%p;A=A*H%p
	return min(A,p-A)
def _find_all_roots_h(p):
	def I(a,b):
		A,B=len(a),len(b);C=[0]*(A+B-1)
		for D in range(A):
			E=a[D]
			if E:
				for F in range(B):C[D+F]+=E*b[F]
		return[A%p for A in C]
	def F(a,m):
		a=a[:];A=len(m);C=1 if m[-1]==1 else pow(m[-1],-1,p)
		while len(a)>=A:
			E=a[-1]if C==1 else a[-1]*C%p;D=len(a)-A
			for B in range(A):a[D+B]=(a[D+B]-E*m[B])%p
			while a and a[-1]==0:a.pop()
		return a
	def P(a,b):
		a=a[:];B=len(b);F=pow(b[-1],-1,p);A=[0]*(len(a)-B+1)
		while len(a)>=B:
			E=a[-1]*F%p;C=len(a)-B;A[C]=E
			for D in range(B):a[C+D]=(a[C+D]-E*b[D])%p
			while a and a[-1]==0:a.pop()
		while A and A[-1]==0:A.pop()
		return A
	def J(b,e,m):
		if len(m)==7 and m==G:
			A=M._ppow_h(b,e,p)
			while A and A[-1]==0:A.pop()
			return A
		B=[1]
		while e:
			if e&1:B=F(I(B,b),m)
			b=F(I(b,b),m);e>>=1
		return B
	def K(a,b):
		while b:a,b=b,F(a,b)
		if a:A=pow(a[-1],-1,p);a=[B*A%p for B in a]
		return a
	G=[A%p for A in M.H_COEFFS];B=J([0,1],p,G)
	if len(B)<2:B+=[0]*(2-len(B))
	B[1]=(B[1]-1)%p
	while B and B[-1]==0:B.pop()
	L=K(G,B)
	if len(L)<4:return[]
	Q=random.Random(p);E=[L];D=[]
	while E:
		A=E.pop()
		if len(A)==2:D.append(-A[0]*pow(A[1],-1,p)%p)
		elif len(A)==3:N=pow(2*A[2],-1,p);O=_sqrt_mod((A[1]*A[1]-4*A[2]*A[0])%p,p);D.append((-A[1]+O)*N%p);D.append((-A[1]-O)*N%p)
		else:
			while True:
				R=Q.randrange(p);C=J([R,1],(p-1)//2,A)
				if C:C[0]=(C[0]-1)%p
				while C and C[-1]==0:C.pop()
				H=K(A,C)
				if 1<len(H)<len(A):E.append(H);E.append(P(A,H));break
	D.sort();return D
def _find_root_h(p):A=_find_all_roots_h(p);return A[0]if A else None
SH_A=20,10,0,50,40,30
BI_A_STEP=512,512,512,512,512,512
BI_A_CAN=512,512,512,0,512,512
MK_A=(1023,)*6
SH_B=51,30,20,10,0,41
BI_B_STEP=1024,1024,512,512,512,512
BI_B_CAN=0,1024,512,512,512,512
MK_B=2047,2047,1023,1023,1023,1023
A_ORD=3,4,5,0,1,2
B_ORD=0,5,1,2,3,4
BIAS_A_STEP_X2=sum(2*BI_A_STEP[A]<<SH_A[A]for A in range(6))
A3_BIAS_SHIFT=512<<50
BIAS_B_STEP_X2=sum(2*BI_B_STEP[A]<<SH_B[A]for A in range(6))
B0_BIAS_SHIFT=1024<<51
def _pack_v(v,sh,bi):
	B=0
	for A in range(6):B|=v[A]+bi[A]<<sh[A]
	return B
def _pack_d(d,sh):
	A=0
	for B in range(6):A=A+(d[B]<<sh[B])&M64
	return A
def _unpack_v(p,sh,bi,mk):return tuple((p>>sh[A]&mk[A])-bi[A]for A in range(6))
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
	return E<<2|G<<1|int(D>0)
def _canon_cone_ref(a,b,S,S2):
	D=a[3]-a[4]+a[5];L=(1,)if D>0 else(-1,)if D<0 else(1,-1);E=None;G=None
	for(B,F)in((a,b),(_mv6(S,a),_mv6(S,b)),(_mv6(S2,a),_mv6(S2,b))):
		H=1
		for I in B_ORD:
			if F[I]!=0:H=1 if F[I]>0 else-1;break
		C=tuple(H*A for A in F);M=C[0],C[1],C[2],C[3],C[4],C[5]
		for J in L:
			A=B if J==1 else(-B[0],-B[1],-B[2],-B[3],-B[4],-B[5]);K=(J*D,A[3],A[3],A[4],A[5],A[0],A[1],A[2])+M
			if E is None or K>E:E=K;G=A,C
	return G
def _floor_div(num,den):
	B=num;A=den
	if A<0:B,A=-B,-A
	return B//A
def _ceil_div(num,den):
	B=num;A=den
	if A<0:B,A=-B,-A
	return-(-B//A)
def _certify_min_gsym(G_sym,G_SC):
	E=G_sym;C=[Frac(0)]*6;B=[[Frac(1 if A==B else 0)for B in range(6)]for A in range(6)]
	for A in range(6):
		F=Frac(E[A][A])-sum(B[A][D]*B[A][D]*C[D]for D in range(A));assert F>0;C[A]=F
		for G in range(A+1,6):B[G][A]=(Frac(E[G][A])-sum(B[G][D]*B[A][D]*C[D]for D in range(A)))/F
	K=5*G_SC;D=[0]*6;H=[10**18];I=[0]
	def J(idx,rem):
		A=idx;F=sum(B[C][A]*D[C]for C in range(A+1,6));G=rem/C[A];K=Frac(math.isqrt(G.numerator*G.denominator)+1,G.denominator);P=C[A]
		for L in range((-F-K).__floor__(),(-F+K).__ceil__()+1):
			M=L+F;N=rem-P*M*M
			if N>=0:
				D[A]=L
				if A==0:
					if any(D):
						O=_qf6(E,D);I[0]+=1
						if O<H[0]:H[0]=O
				else:J(A-1,N)
	J(5,Frac(K));return H[0],I[0]
def verify_v20_lat_b2_and_srt(seed=20261121):
	Az=M.load_data()['t4_geom'];X,Df,Dg,E,R,A8=M.build_exact_int_geom(Az);Dh=E[3]==[0,0,0,0,-1,0]and E[4]==[0,0,0,0,0,-1]and E[5]==[0,0,0,1,0,0];Bu,Fy=_certify_min_gsym(X,A8);Fz=1000*Bu;AR=Frac(Fz)*Frac(13,20)**11;Di=int(round(1764*595.**2*A8));Bv=int((Frac(Q38*Di)/AR).__ceil__());Bw=(Q38*Bv>>64)+1;F_=math.isqrt(int((Frac(Di)/AR).__ceil__()))+1;Dj=int((Frac(3,2)**11*F_).__ceil__());Bx=11*((1<<31)+1)*Dj;Dk=_isqrt_lut4(Bw);G0=Bx+11*Dj+1+Dk+Q32>>32;By=Dk+Q32>>13;Dl=By*By;Z=sorted(_lut4_buckets(),key=lambda e:e[1]);G1=all(_LUT4[Z[A][0]]<=_LUT4[Z[A+1][0]]for A in range(len(Z)-1))and all(Z[A][2]+1==Z[A+1][1]for A in range(len(Z)-1))and Z[0][1]==0 and Z[-1][2]==(1<<40)-1;Dm=sum(1 for(A,C,B)in Z if _LUT4[A]*_LUT4[A]<=B+1<<52);Dn=sum(1 for(A,B,C)in Z if _LUT4[A]>_isqrt_u40_lut7_q32(B+1));G2=max(_LUT4[A]/(sqrt(B+1.)*67108864.)-1. for(A,B,C)in Z);Do=Bu==X[0][0]==753729 and G1 and AR>32*A8 and Bv<=M64 and Bw<1<<40 and Bx<1<<63 and G0<1<<22 and By<1<<32 and Dl<=M64;Bz=random.Random(seed^20903);B_=0;A_=.0
	for Dp in list(range(1,2049))+[(1<<A)-1 for A in range(8,41)]+[1<<A for A in range(8,40)]+[1+M.rand_below(Bz.getrandbits,(1<<40)-1)for A in range(19000)]:
		Dq=_isqrt_lut4(Dp);C0=sqrt(Dp)*(1<<26)
		if Dq<=C0:B_+=1
		Dr=(Dq-C0)/C0
		if Dr>A_:A_=Dr
	C1=0;Ds=[0]*256;C2=[0]*256
	for Dt in range(4096):
		AS=Bz.randrange(256);Du=65535 if Dt<16 else 32768 if Dt<32 else Bz.randrange(65536);_staddh_u16(Ds,AS,Du);C3=C2[AS];C3=C3+(Du&65535)&65535;C2[AS]=C3
		if Ds[AS]!=C2[AS]:C1+=1
	AT=[];J=101
	while len(AT)<6:
		if M._isp(J)and J not in(2,3,19):
			Dv=[A for A in range(J)if _eval_poly(M.H_COEFFS,A,J)==0]
			if len(Dv)>=3:
				Dw=set();B0=pow(3,-1,J)
				for q in Dv:
					if q in Dw:continue
					x=_eval_poly(M.SIGMA_NUM,q,J)*B0%J;A9=_eval_poly(M.SIGMA_NUM,x,J)*B0%J
					if len({q,x,A9})!=3:continue
					Dw.update((q,x,A9));C4=_eval_poly(M.W_COEFFS,q,J);assert _eval_poly(M.W_COEFFS,x,J)==C4 and _eval_poly(M.W_COEFFS,A9,J)==C4;Dx=-((M.U_BN254+1)%J*C4+2)%J
					if Dx==0:continue
					AU=_sqrt_mod(pow(Dx,-1,J),J)
					if AU>0:
						AT.append((J,(q,x,A9),AU))
						if len(AT)==6:break
		J+=2
	def Dy(a,b,Wa,Wb):
		for A in range(6):
			B=Wa*_qf6(Df[A],a)+Wb*_qf6(Df[A],b)
			if B<0:return False
			if B==0:
				C,D=_mv6(Dg[A],a),_mv6(Dg[A],b)
				if _canon_cone_ref(C,D,E,R)<_canon_cone_ref(a,b,E,R):return False
		return True
	def G3(a,b):
		D=_canon_cone_ref(a,b,E,R);A=0
		for(B,C)in((a,b),(_mv6(E,a),_mv6(E,b)),(_mv6(R,a),_mv6(R,b))):
			if _canon_cone_ref(B,C,E,R)==D:
				if B in(a,tuple(-A for A in a))and C in(b,tuple(-A for A in b)):A+=1
		return A
	C5,C6=1,1;AA=18*A8;C7=[]
	for G4 in range(-2,3):
		for G5 in range(-2,3):
			for G6 in range(-2,3):
				for G7 in range(-2,3):
					for G8 in range(-2,3):
						for G9 in range(-3,4):
							C8=G9,G8,G7,G6,G5,G4
							if not any(C8):continue
							Dz=_qf6(X,C8)
							if Dz<=AA:C7.append((C8,Dz))
	C9={};D_={}
	for(GA,GB)in C7:
		for(GC,GD)in C7:
			if GB+GD>AA:continue
			h,i=_canon_cone_ref(GA,GC,E,R)
			if(h,i)not in C9:
				B1=G3(h,i);C9[h,i]=B1
				if Dy(h,i,C5,C6):D_[h,i]=B1
	CA=0;E0=0;E1=0
	for(S,(q,x,A9),AU)in AT:
		E2={};E3={}
		for((h,i),B1)in C9.items():
			AB=0
			for E4 in(q,x,A9):
				a=_eval_zk(h,E4,S);AC=AU*_eval_zk(i,E4,S)%S
				if a==AC:AB+=1
				if a==(S-AC)%S:AB+=1
			if B1==3:AB//=3
			if AB:
				E2[h,i]=AB
				if(h,i)in D_:E3[h,i]=AB
		B0=pow(3,-1,S);CB=[1]
		for f in range(5):CB.append(CB[-1]*q%S)
		E5=[sum(M.ZK_NUM[B][A]*CB[A]for A in range(6))*B0%S for B in range(6)];E6=E5+[-AU*A%S for A in E5];GE=pow(E6[0],-1,S);y=[[S if A==0 else 0 for A in range(12)]]
		for H in range(1,12):
			g=-E6[H]*GE%S
			if g>S//2:g-=S
			y.append([g if A==0 else 1 if A==H else 0 for A in range(12)])
		DESC.bareiss_lll_12d(y,X,C5,C6,False,9,10);GF=[[sum(y[A][C]*X[C][D]*y[B][D]for C in range(6)for D in range(6))+sum(y[A][6+C]*X[C][D]*y[B][6+D]for C in range(6)for D in range(6))for B in range(12)]for A in range(12)];b=[1]+[0]*12;B2=[[0]*12 for A in range(12)]
		for H in range(12):
			for N in range(H+1):
				r=GF[H][N]
				for A in range(N):r=(b[A+1]*r-B2[H][A]*B2[N][A])//b[A]
				if N<H:B2[H][N]=r
				else:b[H+1]=r
		GG=[[_chk_i64((B2[B][A]<<32)//b[A+1])if B>A else 0 for A in range(12)]for B in range(12)];GH=[_chk_u64(min(M64,(b[A+1]<<64)//(b[A]*AA)))for A in range(12)];GI=[_chk_u64(((b[A]*AA<<38)+b[A+1]-1)//b[A+1])for A in range(12)];CC=[0]*12
		for A in range(11,-1,-1):AV=math.isqrt((b[A]*AA+b[A+1]-1)//b[A+1])+1;CC[A]=AV+sum((CC[A]+1)//2 for A in range(A+1,12))
		GJ=[sum(CC[A]for A in range(A+1,12))+1 for A in range(12)];B3={};B4={};E7=[0]*12
		def E8(k,rem,az,cur):
			G=rem;A=0
			for J in range(k+1,12):A=_chk_i64(A-_chk_i64(GG[J][k]*E7[J]))
			Q=_chk_u64(_umulh64(G,GI[k])+1);K=_isqrt_lut4(Q);F=GJ[k];S=(1 if k==0 else 0)if az else _chk_i64(A-F-K+Q32-1)>>32;T=_chk_i64(A+F+K)>>32;L=y[k];U=GH[k]
			for B in range(S,T+1):
				if k==0:
					M=[cur[A]+B*L[A]for A in range(12)];H,I=tuple(M[:6]),tuple(M[6:])
					if any(H)and any(I)and _qf6(X,H)+_qf6(X,I)<=AA:
						C,D=_canon_cone_ref(H,I,E,R);B3[C,D]=B3.get((C,D),0)+1
						if Dy(C,D,C5,C6):B4[C,D]=B4.get((C,D),0)+1
				else:
					N=_chk_u64(abs(_chk_i64((B<<32)-A)));V=N-F if N>F else 0;O=V>>13;P=_umulh64(_chk_u64(O*O),U)
					if P<=G:E7[k]=B;E8(k-1,G-P,az and B==0,[cur[A]+B*L[A]for A in range(12)])
		E8(11,Q38,True,[0]*12)
		if B3!=E2 or B4!=E3:CA+=1
		E0+=sum(B3.values());E1+=sum(B4.values())
	def E9(tb):
		A=1
		for B in B_ORD:
			if tb[B]!=0:A=1 if tb[B]>0 else-1;break
		return tuple(A*B for B in tb)
	def GK(a,b):
		G,L=_mv6(E,a),_mv6(R,a);C=a[3]-a[4]+a[5];H=(1,)if C>0 else(-1,)if C<0 else(1,-1)
		if G==a:M=max(((A*C,A*a[3],A*a[3],A*a[4],A*a[5],A*a[0],A*a[1],A*a[2]),tuple(A*B for B in a))for A in H)[1];N=max(E9(A)for A in(b,_mv6(E,b),_mv6(R,b)));return 3,M,N
		F=None;D=0;I=None
		for(O,B)in((0,a),(1,G),(2,L)):
			for J in H:
				A=B if J==1 else(-B[0],-B[1],-B[2],-B[3],-B[4],-B[5]);K=J*C,A[3],A[3],A[4],A[5],A[0],A[1],A[2]
				if F is None or K>F:F=K;D=O;I=A
		P=b if D==0 else _mv6(E,b)if D==1 else _mv6(R,b);return D,I,E9(P)
	def CD(base,B0,l0,h0,T12_chk=None,verify_ref=True,p_div=None,s0_div=None,rho_div=None):
		O=verify_ref;A=p_div
		if l0>h0:return dict(lvs=0,fast=0,bnd=0,segs=0,bspl=0,fspl=0,cross=0,clips=0,mm=0,ae=0,div_err=0,max_sa=0)
		J=h0-l0+1;K=[base[A]+l0*B0[A]for A in range(12)];C=_mr_pack_a(K[:6])&M64;B=_mr_pack_b(K[6:])<<5&M64;T=_mr_da(B0[:6]);U=_mr_db5(B0[6:]);P=D=Q=0;F=None
		if A is not None:
			V=pow(3,-1,A);L=[1]
			for g in range(5):L.append(L[-1]*s0_div%A)
			F=[sum(M.ZK_NUM[B][A]*L[A]for A in range(6))*V%A for B in range(6)];W=[-rho_div*B%A for B in F];X=F+W;Y=math.sumprod
		for Z in range(J):
			if O or A is not None:
				N=[K[A]+Z*B0[A]for A in range(12)]
				if F is not None and Y(N,X)%A!=0:Q+=1
				G,H=tuple(N[:6]),tuple(N[6:])
				if any(abs(G[A])>AMAX_MR[A]for A in range(6))or any(abs(H[A])>BMAX_MR[A]for A in range(6)):D+=1
				else:
					a=_mr_pack_a(G);I=_mr_pack_b(H)
					if C!=a or B&31!=0 or B>>5!=I or C>=PROD_MA_MR or I>=PROD_MB_MR or I>=1<<59:D+=1
					if O:
						S,b,c=GK(G,H);d=(S<<64|B)>>5&M64
						if d!=S<<59|I:D+=1
						e,f=_canon_cone_ref(G,H,E,R)
						if b!=e or c!=f:P+=1
			C=C+T&M64;B=B+U&M64
		return dict(lvs=J,fast=J,bnd=0,segs=1,bspl=0,fspl=0,cross=0,clips=0,mm=P,ae=D,div_err=Q,max_sa=1)
	CE=[((1,0,0,0,0,0),(0,1,0,0,0,0)),((5,0,0,0,0,0),(1,2,-3,4,-1,2)),((-7,0,0,0,0,0),(0,0,1,-2,3,-4)),((1,1,1,0,0,0),(0,1,-1,2,0,-1))]
	for(CF,GL)in[((0,1,0,0,0,0),(0,1,0,0,0,0)),((0,0,0,1,0,0),(1,-2,3,0,1,-1)),((1,0,1,0,1,0),(0,0,1,0,0,-2))]:GM=_mv6(E,CF);GN=_mv6(R,CF);CG=tuple(CF[A]+GM[A]+GN[A]for A in range(6));assert _mv6(E,CG)==CG;CE.append((CG,GL))
	B5=0
	for(AW,EA)in CE:
		EB=set();CH,CI=AW,EA
		for f in range(3):
			for GO in(1,-1):
				for GP in(1,-1):EB.add(_canon_cone_ref(tuple(GO*A for A in CH),tuple(GP*A for A in CI),E,R))
			CH,CI=_mv6(E,CH),_mv6(E,CI)
		if len(EB)!=1:B5+=1
		if AW[3]==0 and AW[4]==0 and AW[5]==0:
			EC=CD(list(AW)+list(EA),[1,0,0,0,0,0,0,1,0,0,0,0],0,4,verify_ref=True)
			if EC['mm']or EC['ae']:B5+=1
	j=random.Random(seed);T=int(round(1764*568.83**2*A8));z=[math.pi**(A/2.)/math.gamma(A/2.+1.)for A in range(13)];GQ=sqrt(1764.**6*1e3**6*M._det_mat(Az['GA'])**2)*A8**6;B6=49,51,53,55,57,59,61;CJ=0;ED=0;B7=0;B8=0;B9=0;AX=0;BA=0;BB=0;BC=0;CK=[];CL=[];CM=[];CN=[];CO=0;CP=0;CQ=0;CR=0;EE=0;BD=[];CS=999.;BE=0;BF=0;BG=0;BH=0;BI=0
	for GR in B6:
		F=(1<<GR)+1
		while True:
			if M._isp(F):
				AY=_find_root_h(F)
				if AY is not None:
					GS=_eval_poly(M.W_COEFFS,AY,F);EF=-((M.U_BN254+1)%F*GS+2)%F
					if EF!=0:
						BJ=_sqrt_mod(pow(EF,-1,F),F)
						if BJ>0:break
			F+=2
		GT=pow(3,-1,F);CT=[1]
		for f in range(5):CT.append(CT[-1]*AY%F)
		EG=[sum(M.ZK_NUM[B][A]*CT[A]for A in range(6))*GT%F for B in range(6)];EH=EG+[-BJ*A%F for A in EG];GU=pow(EH[0],-1,F);l=[[F if A==0 else 0 for A in range(12)]]
		for H in range(1,12):
			g=-EH[H]*GU%F
			if g>F//2:g-=F
			l.append([g if A==0 else 1 if A==H else 0 for A in range(12)])
		DESC.bareiss_lll_12d(l,X,1764,1000,False,9,10,(99,100));EI=[(_mv6(X,A[:6]),_mv6(X,A[6:]))for A in l];GV=[[1764*math.sumprod(EI[A][0],l[B][:6])+1000*math.sumprod(EI[A][1],l[B][6:])for B in range(12)]for A in range(12)];B=[1]+[0]*12;AZ=[[0]*12 for A in range(12)]
		for H in range(12):
			for N in range(H+1):
				r=GV[H][N]
				for A in range(N):r=(B[A+1]*r-AZ[H][A]*AZ[N][A])//B[A]
				if N<H:AZ[H][N]=r
				else:B[H+1]=r
		EJ=min(B[A+1]/B[A]for A in range(12));CS=min(CS,log2(EJ));assert EJ>=float(AR)*.999999;AD=int((42e2*F*GQ/z[12])**(1./6.));CU=[[_chk_i64((AZ[C][A]<<32)//B[A+1])if C>A else 0 for A in range(12)]for C in range(12)];EK=[_chk_u64(min(M64,(B[A+1]<<64)//(B[A]*AD)))for A in range(12)];GW=[_chk_u64(((B[A]*AD<<38)+B[A+1]-1)//B[A+1])for A in range(12)];CV=[0]*12
		for A in range(11,-1,-1):AV=math.isqrt((B[A]*AD+B[A+1]-1)//B[A+1])+1;CV[A]=AV+sum((CV[A]+1)//2 for A in range(A+1,12))
		GX=[sum(CV[A]for A in range(A+1,12))+1 for A in range(12)];BK=[0]*12;CW=[0]*12;BL=[0]*12;Aa=[0];EL=[[CU[B][A]for B in range(A+1,12)]for A in range(12)];CX=(B[2]<<38)//(B[1]*AD);GY=2*CX;EM=math.sumprod;GZ=[(_mr_da(A[:6]),_mr_db5(A[6:]))for A in l];Ga=A_OFF_MR&M64,B_OFF_MR<<5&M64;AE=[_st_f3(l[A])for A in range(12)];BM=[None]+[[[_add_st_f3(C,_add_st_f3(AE[A],AE[A])if B==2 else AE[A]if B==1 else 0)<<(3 if A==1 else 0)for C in range(729)]for B in range(3)]for A in range(1,12)];Gb=[BM[1][1][A>>3]|((1 if _bad_st_f3(A>>3)else 0)|(2 if _bad_st_f3(_add_st_f3(A>>3,AE[0]))else 0)|(4 if _bad_st_f3(_add_st_f3(A>>3,_add_st_f3(AE[0],AE[0])))else 0))if A>>3<729 else 0 for A in range(8192)]
		def CY(k,rem,az,cur,c0_passed=None,st_cur=0,bmask_passed=None,c21=None):
			V=bmask_passed;U=st_cur;T=c0_passed;J=cur;I=rem;D=c21;nonlocal B7,B8,B9,AX,BI,BA,BB,BC;BK[k]+=1;C=-EM(EL[k],BL[k+1:])if k<11 else 0
			if C<I64_MIN or C>I64_MAX:raise OverflowError(f"i64 ck overflow: {C}")
			if k==0 and T is not None and T!=C:BI+=1
			M=(I*GW[k]>>64)+1
			if M<1 or M>=1<<40:raise OverflowError(f"u40 r2 overflow: {M}")
			W=_isqrt_lut4(M);K=GX[k];A=(1 if k==0 else 0)if az else C-K-W+Q32-1>>32;E=C+K+W>>32
			if A>E:return
			L=l[k];G=GZ[k]
			if k==0:
				if D is not None:
					f=D[0]+A*G[0]&M64;X=D[1]+A*G[1]&M64;Y=[J[B]+A*L[B]for B in range(12)]
					if f!=_mr_pack_a(Y[:6])or X&31!=0 or X>>5!=_mr_pack_b(Y[6:]):BC+=1
				if V is not None:
					for Z in range(A,E+1):
						BB+=1
						if bool(V>>Z%3&1)!=_bad_f3([J[A]+Z*L[A]for A in range(12)]):BA+=1
				H=CD(J,L,A,E,T12_chk=AD,verify_ref=Aa[0]<400,p_div=F,s0_div=AY,rho_div=BJ);Aa[0]+=H['lvs'];B7+=H['mm'];B8+=H['ae'];B9+=H['div_err']
				if H['max_sa']>AX:AX=H['max_sa']
			elif k==1:
				O=E-A+1;g=(2*O*K*CX>>32)+2+O*O;h=(A<<32)-C;a=abs(h)>>13;i=a*a*EK[1]>>64;N=I+g-i;b=((2*A+1)*Q32-2*C)*CX>>32;c=-EM(EL[0][1:],BL[2:])-CU[1][0]*A;P=BM[1][A%3][U]
				for B in range(A,E+1):
					P=Gb[P];CW[0]+=N<0
					if N>=0:BL[1]=B;CY(0,N,az and B==0,[J[A]+B*L[A]for A in range(12)],c,0,P&7,(D[0]+B*G[0]&M64,D[1]+B*G[1]&M64))
					N-=b;b+=GY;c-=CU[1][0]
			else:
				j=EK[k];Q=BM[k][A%3][U]
				for B in range(A,E+1):
					R=(B<<32)-C;d=-R if R<0 else R;m=d-K if d>K else 0;e=m>>13;S=e*e*j>>64;CW[k-1]+=S>I
					if S<=I:BL[k]=B;CY(k-1,I-S,az and B==0,[J[A]+B*L[A]for A in range(12)],None,Q,None,(D[0]+B*G[0]&M64,D[1]+B*G[1]&M64))
					Q=BM[k][1][Q]
		CY(11,Q38,True,[0]*12,None,0,None,Ga);CJ+=Aa[0];ED+=sum(BK);CK.append(round(sum(BK)/max(1,Aa[0]),3));EN=[2.**(log2(z[12-A]/z[12])+.5*(log2(B[A])-A*log2(AD)))for A in range(1,12)];CL.append(round(sum(EN),3));CM.append(round(_fpw([A+B for(A,B)in zip(BK,CW)])/max(1,Aa[0]),2));CN.append(round(_fpw(EN),2));EO=[[_chk_i64((AZ[C][A]<<32)//B[A+1])if C>A else 0 for A in range(12)]for C in range(12)];Gc=[_chk_u64(min(M64,(B[A+1]<<64)//(B[A]*T)))for A in range(12)];Ab=max(log2(B[A+1])-log2(B[A])-log2(T)for A in range(12));_PB['n2']+=1;_PB['fail2']+=Ab>=0;_PB['maxb2']=max(_PB['maxb2'],Ab);EP=[_chk_u64(((B[A]*T<<38)+B[A+1]-1)//B[A+1])for A in range(12)]
		for Gd in EP:BE=max(BE,Gd.bit_length())
		CZ=[0]*12
		for A in range(11,-1,-1):AV=math.isqrt((B[A]*T+B[A+1]-1)//B[A+1])+1;CZ[A]=AV+sum((CZ[A]+1)//2 for A in range(A+1,12))
		BN=[sum(CZ[A]for A in range(A+1,12))+1 for A in range(12)];Ge=[B[A+1]/B[A]for A in range(12)];Gf=sqrt(float(T));Ca=0;Gg=0;Gh=0;Gi=0;Gj=0;Gk=0;Gl=0;Gm=0;EQ=0
		while Ca<300:
			ER=[j.gauss(.0,1.)for A in range(11)];Gn=sum(A*A for A in ER);Go=Gf*j.random()**(1./11.)/sqrt(Gn);Ac=[0]*12;Cb=Q38;ES=True
			for A in range(11,0,-1):
				u=0
				for N in range(A+1,12):u=_chk_i64(u-_chk_i64(EO[N][A]*Ac[N]))
				BH=max(BH,abs(u).bit_length());Ac[A]=int(round(u/float(Q32)+ER[A-1]*Go/sqrt(Ge[A])));Ad=_chk_u64(abs(_chk_i64((Ac[A]<<32)-u)));Gp=Ad-BN[A]if Ad>BN[A]else 0;Ae=Gp>>13;ET=_chk_u64(Ae*Ae);BG=max(BG,ET.bit_length());Af=_umulh64(ET,Gc[A])
				if Af>Cb:ES=False;break
				Cb-=Af
			if not ES:continue
			A0=0
			for N in range(1,12):A0=_chk_i64(A0-_chk_i64(EO[N][0]*Ac[N]))
			BO=_chk_u64(_umulh64(Cb,EP[0])+1);BF=max(BF,BO.bit_length());Ag=_isqrt_lut4(BO);AF=_chk_i64(A0-BN[0]-Ag+Q32-1)>>32;Ah=_chk_i64(A0+BN[0]+Ag)>>32
			if AF>Ah:continue
			Gq=[sum(Ac[A]*l[A][B]for A in range(1,12))for B in range(12)];c=CD(Gq,l[0],AF,Ah,T12_chk=T,verify_ref=True,p_div=F,s0_div=AY,rho_div=BJ)
			if c['lvs']==0:continue
			Gg+=1;Ca+=c['lvs'];Gh+=c['fast'];Gi+=c['bnd'];Gj+=c['segs'];Gk+=c['bspl'];Gl+=c['fspl'];Gm+=c['cross'];EQ+=c['clips'];CP+=c['mm'];CQ+=c['ae'];CR+=c['div_err']
		CO+=Ca;EE+=EQ;Cc=[2.**(log2(z[12-A]/z[12])+.5*(log2(B[A])-A*log2(T)))for A in range(1,12)];Gr=Cc[0];Id=sum(Cc[1:]);Gs=_fpw(Cc);Gt=float(WHEEL_OPS['chord_enter']+WHEEL_OPS['wheel_disp4']);Gu=Gr*Gt;Gv=float(WHEEL_OPS['emit9_step1'])+.015;Gw=Gs+Gu+Gv+M.C_LAT_APPLY;BD.append(round(Gw,2))
	Gx=sum(1./A for A in B6);Gy=round(sum(1./A*B for(A,B)in zip(B6,BD))/Gx,2);Gz=max(BD);EU=round(1e2*max(abs(A-B)/min(A,B)for(A,B)in zip(CK,CL)),2);Cd=round(1e2*max(abs(A-B)/min(A,B)for(A,B)in zip(CM,CN)),2);Ai=32;P=8;EV=Ai*P;Ce=[];EW=0;BP=[];Cf=[];EX=0;EY=HT_SEED_BASE;BQ=M.OP['p_lat'];EZ=M.OP['lx2'];BR=math.sumprod
	for Cg in range(Ai):
		G_=BQ+(EZ-BQ)*Cg/Ai;H0=BQ+(EZ-BQ)*(Cg+1)/Ai;Ea=int(2.**G_)|1;Ch=int(2.**H0)
		if Ch&1==0:Ch-=1
		A1=(Ch-Ea)//2+1;A2=0;A3=[];Eb=[]
		for H1 in range(P):
			Aj=random.Random(EY<<20|Cg<<8|H1);Ci=0
			while True:
				Ci+=1;C=Ea+2*M.rand_below(Aj.getrandbits,A1)
				if C&7 not in(1,3)or not M._isp(C):continue
				s=M._ppow_h([0,1,0,0,0,0],C,C);s[1]=(s[1]-1)%C;Ec=M.H_COEFFS if not any(s)else M._pgcd(s,M.H_COEFFS,C)
				if len(Ec)<=1:continue
				Ed=M._roots_mod(Ec,C)
				if not Ed:continue
				BS=[];Ie=M.U_BN254%C
				for Ee in Ed:
					H2=_eval_poly(M.W_COEFFS,Ee,C);Ef=-((M.U_BN254+1)%C*H2+2)%C
					if Ef!=0:
						Ak=_sqrt_mod(pow(Ef,-1,C),C)
						if Ak>0:BS.append((Ee,Ak))
				if BS:Eg=len(BS)//3;H3,Ak=Aj.choice(BS);H4=Ak if Aj.randrange(2)==0 else C-Ak;break
			A2+=Ci;EX+=Ci;H5=pow(3,-1,C);Cj=[1]
			for f in range(5):Cj.append(Cj[-1]*H3%C)
			Eh=[sum(M.ZK_NUM[B][A]*Cj[A]for A in range(6))*H5%C for B in range(6)];Ei=Eh+[-H4*A%C for A in Eh];H6=pow(Ei[0],-1,C);BT=[[C if A==0 else 0 for A in range(12)]]
			for BU in range(1,12):
				g=-Ei[BU]*H6%C
				if g>C//2:g-=C
				BT.append([g if A==0 else 1 if A==BU else 0 for A in range(12)])
			f,O,H7=DESC.bareiss_lll_12d(BT,X,1764,1000,False,9,10,(99,100),True);Ck=[2.**(log2(z[12-A]/z[12])+.5*(log2(O[A])-A*log2(T)))for A in range(1,12)];H8,If=Ck[0],sum(Ck[1:]);H9=[[(H7[B][A]<<32)//O[A+1]if B>A else 0 for A in range(12)]for B in range(12)];HA=[min(M64,(O[A+1]<<64)//(O[A]*T))for A in range(12)];Ab=max(log2(O[A+1])-log2(O[A])-log2(T)for A in range(12));_PB['n']+=1;_PB['fail']+=Ab>=0;_PB['maxb']=max(_PB['maxb'],Ab);HB=[((O[A]*T<<38)+O[A+1]-1)//O[A+1]for A in range(12)];Cl=[0]*12
			for A in range(11,-1,-1):HC=math.isqrt((O[A]*T+O[A+1]-1)//O[A+1])+1;Cl[A]=HC+sum((Cl[A]+1)//2 for A in range(A+1,12))
			Ej=[sum(Cl[A]for A in range(A+1,12))+1 for A in range(12)];HD=[O[A+1]/O[A]for A in range(12)];HE=sqrt(float(T));HF=[1./sqrt(HD[A])for A in range(12)];Ek=[[H9[B][A]for B in range(A+1,12)]for A in range(12)];HG=[[BT[B][A]for B in range(1,12)]for A in range(12)];El=BV=0;BW=BX=Cm=.0;HH=Aj.gauss;HI=Aj.random;HJ=1./float(Q32);Em=BT[0];HK=HB[0];En=Ej[0];HL=Ek[0];AG=WHEEL_OPS;Cn=C>=EMIT8_P_MIN;BY=EMIT8_FLUSH_AMORT if Cn else .015;Co=float(AG['chord_enter']);Eo=float(AG['emit8_step1'if Cn else'emit9_step1']);Ep=float(AG['emit8_last'if Cn else'emit9_last']);Al=float(AG['wheel_disp4']);HM=float(AG['wheel_skip1_stub']);HN=float(AG['wheel_skip2_stub']);HO=math.isqrt((O[0]*T+O[1]-1)//O[1])+1;Am=.0 if 2*HO+3<=128 else 2.
			while BV<160:
				Cp=[HH(.0,1.)for A in range(11)];HP=HE*HI()**(1./11.)/sqrt(BR(Cp,Cp));s=[0]*12;Cq=Q38;Eq=True
				for A in range(11,0,-1):
					u=-BR(Ek[A],s[A+1:])if A<11 else 0;Er=int(round(u*HJ+Cp[A-1]*HP*HF[A]));s[A]=Er;Ad=abs((Er<<32)-u);Es=Ej[A]
					if Ad>Es:
						Ae=Ad-Es>>13;Af=Ae*Ae*HA[A]>>64
						if Af>Cq:Eq=False;break
						Cq-=Af
				if not Eq:continue
				A0=-BR(HL,s[1:]);BO=(Cq*HK>>64)+1;Ag=math.isqrt(BO<<52)+1;AF=A0-En-Ag+Q32-1>>32;Ah=A0+En+Ag>>32
				if AF>Ah:continue
				HQ=s[1:];Cr=[BR(HQ,HG[A])for A in range(12)];BZ=(1 if _bad_f3(Cr)else 0)|(2 if _bad_f3([Cr[A]+Em[A]for A in range(12)])else 0)|(4 if _bad_f3([Cr[A]+2*Em[A]for A in range(12)])else 0);AH=Ah-AF+1;El+=1;BV+=AH
				if BZ==7:BW+=Co+Al+Am
				elif BZ==0:Cs=(AH+127)//128;BW+=Co+Al+Am+(Cs-1)*(Al+Am+2e1);BX+=(AH-1)*(Eo+BY)+(Ep+BY);Cm+=AH*M.C_LAT_APPLY
				else:
					BW+=Co+Al+Am;Ba=(BZ|BZ<<3)>>AF%3&7;HR=Ba in(3,5,6);HS,HT=divmod(AH,3);Bb=HS*(1 if HR else 2)+sum(1 for A in range(HT)if not Ba>>A&1);Cs=(AH+127)//128;An=.0
					if Bb>=1:
						if Ba&1==0:An=.0
						elif Ba&2==0:An=HM
						else:An=HN
					BX+=An+(Cs-1)*(Al+Am+An+2e1)
					if Bb>=1:BX+=(Bb-1)*(Eo+BY)+(Ep+BY);Cm+=Bb*M.C_LAT_APPLY
			EW+=BV;Et=_fpw(Ck)+H8*(BW/El);HU=(BX+Cm)/BV;Eu=Et+HU;Ce.append(round(Eu,2));A3.append((C,Eg,Eu));Eb.append((C,Eg,Et))
		BP.append((A1,A2,A3));Cf.append((A1,A2,Eb))
	Ev=Ao=Ew=.0
	for((A1,A2,A3),(f,f,Ct))in zip(BP,Cf):k=A1*(P-1.)/(A2-1.);Ev+=k*sum(B*C/A for(A,B,C)in A3)/P;Ew+=k*sum(B*C/A for(A,B,C)in Ct)/P;Ao+=k*sum(B/A for(A,B,C)in A3)/P
	Ap=Ev/Ao;Bc=Ew/Ao;Ex=Ey=.0
	for((A1,A2,A3),(f,f,Ct))in zip(BP,Cf):k=A1*(P-1.)/(A2-1.);Ez=[B/A*(C-Ap)for(A,B,C)in A3];Cu=sum(Ez)/P;HV=sum((A-Cu)**2 for A in Ez)/(P-1);Bd=k*k/float(P-2);Ex+=(k*k+Bd)*(HV/P)+Bd*(Cu*Cu);E_=[B/A*(C-Bc)for(A,B,C)in Ct];Cv=sum(E_)/P;HW=sum((A-Cv)**2 for A in E_)/(P-1);Ey+=(k*k+Bd)*(HW/P)+Bd*(Cv*Cv)
	Cw=sqrt(Ex)/Ao;F0=sqrt(Ey)/Ao;HX=round(Ap+3.*Cw,3);Cx=1.+max(.03,Cd/1e2);HY=round(Ap+(Cx-1.)*(Bc+3.*F0)+Cx*3.*Cw,3);HZ=Az['GA'];U=[HZ[A][:]+[1. if A==B else .0 for B in range(6)]for A in range(6)]
	for d in range(6):
		F1=max(range(d,6),key=lambda r:abs(U[r][d]));U[d],U[F1]=U[F1],U[d];Ha=U[d][d];U[d]=[A/Ha for A in U[d]]
		for Aq in range(6):
			if Aq!=d and U[Aq][d]:Hb=U[Aq][d];U[Aq]=[U[Aq][A]-Hb*U[d][A]for A in range(12)]
	F2=[sqrt(U[A][6+A])for A in range(6)];F3=[round(572.*sqrt(1.+M.SC_GUARD_DELTA)*A,1)for A in F2];F4=[round(572.*sqrt(1.764*(1.+M.SC_GUARD_DELTA))*A,1)for A in F2];F5=all(F3[A]<AMAX_MR[A]+1. for A in range(6))and all(F4[A]<BMAX_MR[A]+1. for A in range(6))and PROD_MB_MR<1<<59;K=1024;Hc=16;Hd=16;Ig=6;A4=[]
	for He in range(0,4):
		for Hf in range(-4,5):
			for Hg in range(-4,5):
				for Hh in range(-2,3):
					Be=Hh,0,0,He,Hf,Hg
					if any(Be):Hi=sum(Be[A]*X[A][B]*Be[B]for A in range(6)for B in range(6));A4.append((Hi,Be))
					if len(A4)>=K:break
				if len(A4)>=K:break
			if len(A4)>=K:break
		if len(A4)>=K:break
	A4.sort(key=lambda x:x[0]);Cy=[A for(B,A)in A4];F6=[_pack_v(A,SH_A,BI_A_CAN)for A in Cy]
	def F7(pa_v):
		A=pa_v;B=A>>50&511;C=(A>>40&1023)-508;D=(A>>30&1023)-508;E=(A>>20&1023)-510
		if A>>10&1023!=512 or A&1023!=512 or not(0<=B<4 and 0<=C<9 and 0<=D<9 and 0<=E<5):return-1
		return B<<11|C<<7|D<<3|E
	Bf=[-1]*8192
	for(Hj,Hk)in enumerate(F6,1):Cz=F7(Hk);assert Cz>=0 and Bf[Cz]==-1;Bf[Cz]=Hj
	def Hl(av):A=av;B=A[0]&65535|(A[1]&65535)<<16|(A[2]&65535)<<32|(A[3]&65535)<<48;C=A[4]&65535|(A[5]&65535)<<16;return _chk_u64(B),_chk_u64(C)
	def AI(v):v&=65535;return v-65536 if v>=32768 else v
	F8=[None]+[Hl(A)for A in Cy]
	def F9(w0,w1,c_k,p_val,bias_val,inv_p):C=p_val;B=c_k;F=AI(w0);G=AI(w0>>16);H=AI(w0>>32);I=AI(w0>>48);J=AI(w1);K=AI(w1>>16);A=_chk_i64(F*B[0]);A=_chk_i64(A+_chk_i64(G*B[1]));A=_chk_i64(A+_chk_i64(H*B[2]));A=_chk_i64(A+_chk_i64(I*B[3]));A=_chk_i64(A+_chk_i64(J*B[4]));A=_chk_i64(A+_chk_i64(K*B[5]));E=_chk_u64(A+bias_val);L=_umulh64(E,inv_p);D=_chk_u64(E-_chk_u64(L*C));assert D<2*C;return D-C if D>=C else D
	C_=0;D0=0;Bg=0
	for Hm in(42,46,48):
		L=(1<<Hm)+1
		while True:
			if M._isp(L):
				D1=_find_root_h(L)
				if D1 is not None:
					Hn=_eval_poly(M.W_COEFFS,D1,L);FA=-((M.U_BN254+1)%L*Hn+2)%L
					if FA!=0:
						FB=_sqrt_mod(pow(FA,-1,L),L)
						if FB>0:break
			L+=2
		Ho=pow(3,-1,L);D2=[1]
		for f in range(5):D2.append(D2[-1]*D1%L)
		FC=[sum(M.ZK_NUM[B][A]*D2[A]for A in range(6))*Ho%L for B in range(6)];FD=-pow(FB,-1,L)%L;Hp=[A*FD%L for A in FC];FE=_chk_u64(3072*L);Hq=_chk_u64((1<<64)//L);Bg=max(Bg,(2*FE).bit_length())
		for D in range(K,0,-1):
			D3,D4=F8[D];Hr=F9(D3,D4,Hp,L,FE,Hq);a=Cy[D-1];Hs=sum(a[A]*FC[A]for A in range(6))*FD%L
			if Hr!=Hs:C_+=1
			D0+=1
	I=127
	while True:
		if M._isp(I):
			D5=_find_root_h(I)
			if D5 is not None:
				Ht=_eval_poly(M.W_COEFFS,D5,I);FF=-((M.U_BN254+1)%I*Ht+2)%I
				if FF!=0:
					FG=_sqrt_mod(pow(FF,-1,I),I)
					if FG>0:break
		I+=2
	Hu=pow(3,-1,I);D6=[1]
	for f in range(5):D6.append(D6[-1]*D5%I)
	Hv=[sum(M.ZK_NUM[B][A]*D6[A]for A in range(6))*Hu%I for B in range(6)];Hw=-pow(FG,-1,I)%I;Hx=[A*Hw%I for A in Hv];Hy=_chk_u64(3072*I);Hz=_chk_u64((1<<64)//I);G=16;Ar=(I+G-1)//G;m=[0]*(K+1);D7=[0]*Ar;FH=[0]*(K+1);H_=(1<<47)-1
	for D in range(K,0,-1):D3,D4=F8[D];A=F9(D3,D4,Hx,I,Hy,Hz);m[D]=A;Q=A//G;FH[D]=D7[Q];D7[Q]=(A&G-1)<<47|D
	D8=0
	for Q in range(Ar):
		As=j.randint(40,K);FI=set();Bh=D7[Q]
		while Bh:
			D=Bh&H_
			if D>As:break
			FI.add((Bh>>47,D));Bh=FH[D]
		I0={(m[A]&G-1,A)for A in range(1,As+1)if m[A]//G==Q}
		if FI!=I0:D8+=1
	n=[0]*(Ar+1)
	for D in range(1,K+1):n[m[D]//G+1]+=1
	for Q in range(Ar):n[Q+1]+=n[Q]
	AJ=[0]*K;I1=n[-1]
	for D in range(K,0,-1):A=m[D];Q=A//G;n[Q+1]-=1;AJ[n[Q+1]]=D<<17|A&G-1
	n=n[1:]+[I1];D9=0;v=1<<11
	for Q in range(Ar):
		As=j.randint(40,K);FJ=As<<17|15
		for Bi in(1,-1):
			DA=[0]*(G*v);DB=[0]*(G*v);I2,FK=n[Q],n[Q+1];t=I2;Bj=0 if Bi>0 else(G-1)*v
			while t+4<=FK:
				I3,BU,I4,FL=AJ[t],AJ[t+1],AJ[t+2],AJ[t+3]
				if FL>FJ:break
				for AK in(I3,BU,I4,FL):o=AK&15;DC=Bj+(o<<11)if Bi>0 else Bj-(o<<11);_staddh_u16(DA,DC+(AK>>16>>1),720)
				t+=4
			while t<FK:
				AK=AJ[t];t+=1
				if AK>FJ:break
				o=AK&15;DC=Bj+(o<<11)if Bi>0 else Bj-(o<<11);_staddh_u16(DA,DC+(AK>>16>>1),720)
			for D in range(1,As+1):
				if m[D]//G==Q:
					o=m[D]&G-1
					if Bi<0:o=G-1-o
					DB[o*v+D]=DB[o*v+D]+720&65535
			if DA!=DB:D9+=1
	V=Hc;AL=4;p=0
	for FM in(0,1,(1<<34)+7,(1<<35)-1):
		if _chk_u64(FM<<15)>=1<<50:p+=1
		if _chk_u64(FM<<14)>=1<<49:p+=1
	for Bk in(0,1,32768,(1<<28)-1):
		DD=_chk_u64((1<<35)-1<<29|Bk)
		if DD>>29!=(1<<35)-1 or DD<<35&M64!=_chk_u64(Bk<<35):p+=1
		if _chk_u64((DD<<35&M64)+((1<<49)-2)&M64)!=(Bk<<35)+((1<<49)-2):p+=1
		if _chk_u64((Bk<<35)+((1<<49)-2))>=1<<64:p+=1
	DE=0;A5=V*Hd*8;At=[0]*(A5+1);AM=[1]+[0]*(V-1)
	for Au in range(A5):
		if Au&7==0:I5=min(range(V),key=lambda s:AM[s]);DF=I5
		FN=AM[DF];AM[DF]=FN+1;At[Au]=_chk_u64(DF<<29|FN)
	At[A5]=0
	if any(At[A]==0 for A in range(A5)):DE+=1
	FO=[[[0]*(K+1)for A in range(AM[A])]for A in range(V)];DG=[[[0]*(K+1)for A in range(AM[A])]for A in range(V)];DH=2046;I6=(1<<35)-1;A6=[16]*V;DI=[bytearray(16384)for A in range(V)];Av=0;Bl=0;Aw=0;FP=struct.Struct('<Q').pack_into;FQ=struct.Struct('<Q').unpack_from;I7=struct.Struct('<QQ').pack_into;I8=struct.Struct('<QQ').unpack_from
	def FR(sid_v,w_lp_v,cnt_v):
		B=sid_v;nonlocal Aw;A=DI[B];I7(A,0,cnt_v|w_lp_v<<16,0);C,G=I8(A,0);D=C&65535;H=C>>16&65535;assert 1<=D<=DH;Aw+=1;E=16;I=FO[B]
		for G in range(D):F=FQ(A,E)[0];E+=8;_staddh_u16(I[F>>35],(F&I6)>>1,H)
	FS=M.rand_below;FT=j.getrandbits
	for FU in range(224):
		DJ=752+FU;I9=37000 if FU in(0,223)else 18
		for IA in range(I9):
			D=1+FS(FT,K);Au=A5 if IA%500==0 else FS(FT,A5);IB=F6[D-1];DK=At[Au];Y=DK>>29;FV=DK&(1<<29)-1
			if Au<A5:DG[Y][FV][D]=DG[Y][FV][D]+DJ&65535
			IC=Bf[F7(IB)];FW=_chk_u64((DK<<35&M64)+(IC<<1))
			if FW>=1<<64:p+=1
			FP(DI[Y],A6[Y],FW);A6[Y]+=8
			if A6[Y]&16384:
				if A6[Y]!=16384:p+=1
				FR(Y,DJ,DH);A6[Y]=16;Av+=1
		for Y in range(V):
			FX=A6[Y]-16>>3
			if FX:FR(Y,DJ,FX);A6[Y]=16;Bl+=1
	DI=None;Bf=None;At=None;Bm=0
	for Bn in range(V):
		ID=1 if Bn==0 else 0
		for FY in range(ID,AM[Bn]):
			if FO[Bn][FY]!=DG[Bn][FY]:Bm+=1
	FZ=Bm==0 and p==0 and Aw==Av+Bl;IE=V;Fa=[bytearray(16384)for A in range(AL)];Fb=[bytearray(16384)for A in range(AL)];DL=[0]*AL
	for AN in range(1,K+1):w=m[AN]//G%AL;Fc=DL[w];FP(Fa[w],Fc*8,AN<<17|m[AN]&G-1);DL[w]=Fc+1
	DM=0;Fd=[0]*(G*v);DN=[0]*(G*v)
	for AN in range(1,K+1):Fe=(m[AN]&G-1)*v+AN;DN[Fe]=DN[Fe]+720&65535
	for w in range(AL):
		Fb[w][:]=Fa[w]
		for IF in range(DL[w]):Ff=FQ(Fb[w],IF*8)[0];_staddh_u16(Fd,((Ff&15)<<11)+(Ff>>16>>1),720)
	if Fd!=DN:DM+=1
	Fg=256;IG=V;DO=0;DP=[0]*V;DQ=[j.randint(1,64)for A in range(Fg)]
	for IH in DQ:II=min(range(V),key=lambda s:DP[s]);DP[II]+=IH
	if max(DP)>sum(DQ)/float(V)+max(DQ)+1e-09:DO+=1
	Bo=Az['GB'];IJ=Bo[0][0];Bp=[[Bo[A][B]-Bo[A][0]*Bo[0][B]/IJ for B in range(1,6)]for A in range(1,6)];IK=sqrt(max(sum((.5 if A>>B&1 else-.5)*Bp[B][C]*(.5 if A>>C&1 else-.5)for B in range(5)for C in range(5))for A in range(32)));IL=8.*math.pi*math.pi/15./sqrt(M._det_mat(Bp));W=[Bp[A][:]+[1. if A==B else .0 for B in range(5)]for A in range(5)]
	for e in range(5):
		Fh=max(range(e,5),key=lambda r:abs(W[r][e]));W[e],W[Fh]=W[Fh],W[e];IM=W[e][e];W[e]=[A/IM for A in W[e]]
		for Ax in range(5):
			if Ax!=e and W[Ax][e]:IN=W[Ax][e];W[Ax]=[W[Ax][A]-IN*W[e][A]for A in range(10)]
	IO=[sqrt(W[A][5+A])for A in range(5)];DR=0
	for Bq in(4.,8.,12.):
		Fi=[int(floor(Bq*IO[A]))for A in range(5)];IP=Bq*Bq;Fj=[0]
		def Fk(ix,v5):
			A=ix
			if A==5:
				if sum(v5[A]*Bp[A][B]*v5[B]for A in range(5)for B in range(5))<=IP+1e-12:Fj[0]+=1
				return
			for B in range(-Fi[A],Fi[A]+1):v5[A]=B;Fk(A+1,v5)
		Fk(0,[0]*5)
		if Fj[0]>IL*(Bq+IK)**5:DR+=1
	def DS(a,b):
		B=[0]*11
		for C in range(6):
			for A in range(6):B[C+A]+=a[C]*b[A]
		for D in range(10,5,-1):
			E=B[D]
			if E:
				for A in range(6):B[D-6+A]-=E*M.H_COEFFS[A]
		return B[:6]
	def IQ(M6):
		A=[A[:]for A in M6];C=1;G=1
		for B in range(5):
			if A[B][B]==0:
				D=next((C for C in range(B+1,6)if A[C][B]!=0),-1)
				if D<0:return 0
				A[B],A[D]=A[D],A[B];C=-C
			H=A[B][B]
			for E in range(B+1,6):
				for F in range(B+1,6):A[E][F]=(A[E][F]*H-A[E][B]*A[B][F])//G
			G=H
		return C*A[5][5]
	def Fl(a,b,u_val=M.U_BN254):
		B=[sum(a[A]*M.ZK_NUM[A][B]for A in range(6))for B in range(6)];C=[sum(b[A]*M.ZK_NUM[A][B]for A in range(6))for B in range(6)];D=DS(B,B);J=[(u_val+1)*A+2*B+C for(A,B,C)in zip(DS(M.W_COEFFS,D),D,DS(C,C))];E=[];A=J[:]
		for K in range(6):
			E.append(A)
			if K<5:
				F=[0]+A[:5];G=A[5]
				if G:
					for H in range(6):F[H]-=G*M.H_COEFFS[H]
				A=F
		I=abs(IQ([[E[B][A]for B in range(6)]for A in range(6)]));assert I%9**6==0;return I//9**6
	def Fm(P0,xs):
		A=[list(xs)]
		while len(A[-1])>1:B=A[-1];A.append([B[A]*B[A+1]if A+1<len(B)else B[A]for A in range(0,len(B),2)])
		E=[P0%A[-1][0]]
		for H in range(len(A)-2,-1,-1):E=[E[A>>1]%B for(A,B)in enumerate(A[H])]
		F=[]
		for(G,I)in zip(xs,E):
			C=G;D=math.gcd(C,I)
			while D>1:C//=D;D=math.gcd(C,D)
			F.append((C,G//C))
		return F
	AO=[]
	while len(AO)<64:
		a=[j.randint(-220,220)for A in range(6)];A7=[j.randint(-300,300)for A in range(6)]
		if any(a)and any(A7)and M.is_ideal_primitive_oh(a,A7):AO.append(Fl(a,A7))
	Fn=M._primes_upto(1<<14);Fo=1
	for AP in Fn:Fo*=AP
	Fp=[]
	for Ay in range(0,64,16):Fp.extend(Fm(Fo,AO[Ay:Ay+16]))
	DT=0
	for(Fq,(IR,IS))in enumerate(Fp):
		Br=AO[Fq]
		for AP in Fn:
			while Br%AP==0:Br//=AP
		if Br!=IR or AO[Fq]//Br!=IS:DT+=1
	def IT(n):
		if n<2:return False
		for C in(2,3,5,7,11,13,17,19,23,29,31,37):
			if n==C:return True
			if n%C==0:return False
		B=n-1;D=0
		while B&1==0:B>>=1;D+=1
		for E in(2,3,5,7,11,13,17,19,23,29,31,37):
			if E>=n:continue
			A=pow(E,B,n)
			if A in(1,n-1):continue
			for F in range(D-1):
				A=A*A%n
				if A==n-1:break
			else:return False
		return True
	def IU(n):
		if n&1==0:return 2
		for D in range(1,50):
			B=2;A=2;C=1
			while C==1:B=(B*B+D)%n;A=(A*A+D)%n;A=(A*A+D)%n;C=math.gcd(abs(B-A),n)
			if C!=n:return C
		return n
	def IV(n,small_ps):
		C={}
		for A in small_ps:
			if A*A>n:break
			while n%A==0:C[A]=C.get(A,0)+1;n//=A
		D=[n]if n>1 else[]
		while D:
			B=D.pop()
			if B==1:continue
			if IT(B):C[B]=C.get(B,0)+1
			else:E=IU(B);D.append(E);D.append(B//E)
		return C
	DU=1<<8;DV=1<<20;IW=M._primes_upto(1<<10);IX=M._primes_upto(DU);Fr=1
	for AP in IX:Fr*=AP
	AQ=[]
	while len(AQ)<320:
		a=[j.randint(-4,4)for A in range(6)];A7=[j.randint(-6,6)for A in range(6)]
		if any(a)and any(A7)and M.is_ideal_primitive_oh(a,A7):
			Bs=Fl(a,A7,u_val=5)
			if Bs>=1<<25:AQ.append(Bs)
	Fs=[]
	for Ay in range(0,len(AQ),32):Fs.extend(Fm(Fr,AQ[Ay:Ay+32]))
	DW=0;Ft=0;DX=0;Fu=0;DY=0;DZ=0;Da=1e9;Db=1e9
	for(Bs,(Dc,f))in zip(AQ,Fs):
		IY=log2(Bs);Bt=IV(Dc,IW);IZ=sum(floor(16.*log2(A))for A in Bt if DU<=A<=DV);Ia=any(B>=2 for(A,B)in Bt.items()if DU<=A<=DV);Ib=floor(ceil(IY)/8.);AC=(Ib+1)/16.;Dd=(log2(Dc)if Dc>1 else .0)-IZ/16.;De=[A for A in Bt if A>DV];Fv=Dd<=AC+.25
		if not De:
			if not Ia:
				DW+=1;Fw=AC+.25-Dd
				if Fw<Da:Da=Fw
				if not Fv:DY+=1
			else:Ft+=1
		else:
			if len(De)==1 and Bt[De[0]]==1:DX+=1
			else:Fu+=1
			Fx=Dd-(AC+.25)
			if Fx<Db:Db=Fx
			if Fv:DZ+=1
	Ic=Dh and Do and B_==0 and A_<=.0012 and Dm==0 and Dn==0 and C1==0 and D9==0 and CA==0 and B5==0 and C_==0 and D0>=3000 and Bg<64 and D8==0 and Bm==0 and p==0 and DE==0 and DR==0 and DM==0 and FZ and DO==0 and Av>=20 and Aw==Av+Bl and BI==0 and F5 and CJ>=10000 and B7==0 and B8==0 and B9==0 and AX<=4 and BA==0 and BB>=10000 and BC==0 and EU<=5. and Cd<=5. and CO>=2000 and CP==0 and CQ==0 and CR==0 and BE<64 and BF<40 and BG<64 and BH<63 and DT==0 and DW>=10 and DX>=100 and DY==0 and DZ==0;return dict(ok=Ic,smat_cyclic_a345_ok=Dh,zero_float_fp_enum=True,strict_u64_i64_umulh_checked=True,zero_exception_64bit_q38_proved=Do,cert_min_gsym_exact=Bu,cert_min_gsym_pts_le_5gsc=Fy,uncond_min_bk_star_bits=round(log2(float(AR)),3),uncond_inv_max_bits=Bv.bit_length(),uncond_r2_q12_max_bits=Bw.bit_length(),uncond_ck_max_bits=Bx.bit_length(),uncond_d2_max_bits=Dl.bit_length(),isqrt_u40_lut7_errs=B_,lut4_buckets=len(Z),lut4_bucket_violations=Dm,lut4_buckets_above_lut7=Dn,lut4_max_overshoot_pct=round(G2*1e2,4),isqrt_u40_lut7_max_over_pct=round(A_*1e2,4),staddh_16b_semantic_errs=C1,b2_x4_prefix_staddh_errs=D9,orbits_6ideal_checked=len(AT),lattice_b2_mismatches=CA,lattice_b2_raw_hits=E0,lattice_b2_jvor_hits=E1,complete_root_to_leaf_orbits=len(B6),complete_root_to_leaf_leaves=CJ,complete_root_to_leaf_nodes=ED,complete_root_to_leaf_mismatches=B7,complete_root_to_leaf_apply_errs=B8,complete_root_to_leaf_div_errs=B9,complete_root_to_leaf_max_sa_cones=AX,wheel_tree_leaves_checked=BB,wheel_tree_state_errs=BA,swar21_cur_errs=BC,complete_scaled_nodes_per_leaf_obs=CK,complete_scaled_nodes_per_leaf_pred=CL,complete_scaled_max_rel_err_pct=EU,complete_scaled_fp_ir_obs=CM,complete_scaled_fp_ir_pred=CN,complete_scaled_fp_ir_max_rel_err_pct=Cd,prod_exec_leaves=CO,prod_exec_mismatches=CP,prod_exec_apply_errs=CQ,prod_exec_div_errs=CR,prod_exec_endpoint_clips=EE,prod_min_bk_star_bits=round(CS,2),prod_max_inv_b2_bits=BE,prod_max_r2_q12_bits=BF,prod_max_d2_bits=BG,prod_max_ck_bits=BH,prod_12lvl_analytic_ir_per_leaf=BD,prod_inv_bit_weighted_mean_ir=Gy,prod_max_ir=Gz,prod_256orb_strata=Ai,prod_256orb_draws_per_stratum=P,prod_256orb_count=EV,prod_256orb_leaves_sampled=EW,prod_256orb_min_ir=min(Ce),prod_256orb_max_ir=max(Ce),prod_256orb_wmean_ht_ir=round(Ap,3),prod_256orb_se_strat_ir=round(Cw,3),prod_256orb_gauss_wmean_ht_ir=round(Bc,4),prod_256orb_gauss_se_strat_ir=round(F0,4),prod_256orb_exact_wmean_ht_ir=round(Ap-Bc,4),prod_256orb_wmean_plus_3se_ir=HX,prod_256orb_extrap_factor=round(Cx,4),prod_256orb_extrap_upper_bound_ir=HY,prod_256orb_exact_rejection_ht=True,prod_256orb_uniform_random_root_choice=True,prod_256orb_run=EY,prod_256orb_strata_ht=BP,beta_ge1_census=dict(orbits_256=_PB['n'],ge1_256=_PB['fail'],max_log2_beta_256=round(_PB['maxb'],3),orbits_prod=_PB['n2'],ge1_prod=_PB['fail2'],max_log2_beta_prod=round(_PB['maxb2'],3)),prod_256orb_avg_rej_tries=round(EX/float(EV),1),prod_budget_c_lat_leaf=M.C_LAT_LEAF,swar_tie_cases_checked=len(CE),swar_tie_canon_errs=B5,barrett_kta_18op_tested=D0,barrett_kta_18op_errs=C_,barrett_max_sp_bits=Bg,desc_ai_prefix_chain_prime=I,desc_ai_prefix_chain_errs=D8,score_node_replay_errs=Bm,fabric_1hop_addr_errs=p,b_tab_packing_errs=DE,discrete_chord_bound_test_errs=DR,fabric_1hop_list_node_errs=DM,fabric_1hop_per_node_bytearray_zero_dict=True,fabric_1hop_zero_frag_exact_ok=FZ,score_nodes_sim=IE,greedy_bound_subbands=Fg,greedy_bound_nodes=IG,greedy_node_bound_errs=DO,fabric_1hop_list_nodes_sim=AL,logp_bands_1_16b=224,lat_8b_direct_ptr_page_cap=DH,emit_full_pages_flushed=Av,emit_tail_pages_flushed=Bl,emit_pages_streamed=Aw,zero_hash_table_apply=True,inc_c0_exact_errs=BI,swar_bounds_ok=F5,max_a_prod_coords=F3,max_b_prod_coords=F4,srt_p24_genuine_norms_checked=len(AO),srt_batch_mismatches=DT,e2e_factored_c6_norms_checked=len(AQ),e2e_factored_smooth_sqfree=DW,e2e_factored_smooth_p2=Ft,e2e_factored_1lp=DX,e2e_factored_multi_lp=Fu,gate_0lp_fn=DY,gate_0lp_fp=DZ,gate_0lp_min_smooth_margin_b=round(Da,4),gate_0lp_min_1lp_gap_b=round(Db,2))
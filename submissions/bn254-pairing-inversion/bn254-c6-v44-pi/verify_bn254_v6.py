Ad='two_sided_surv_max_over_up'
b='N_l_U'
Z='N_sw_U'
V='N_sw'
U='ncols'
Aj='descent'
Be='pool_quad'
Bd=bytearray
AF='boot'
AE='special_q'
Ap=reversed
a='frac_reg'
A7='res_fano_v_over_m'
T='texc'
N='surplus'
AQ='EMB_ZK'
AP='t4_geom'
AY=next
t='mean_M_over_mu'
s='pgf_ok'
r='ok'
A6=list
o=complex
l=sorted
h=float
j=all
g=enumerate
c=zip
Y=dict
X=pow
W=any
S=None
O=tuple
R=min
P=int
L=True
K=.0
J=abs
I=False
H=round
G=1.
B=max
D=sum
A=len
C=range
import bisect as AX,functools as CK,heapq as z,json,math as Q,random
from fractions import Fraction as F
from math import ceil,exp,floor,log2,pi,sqrt
from pathlib import Path
import model_bn254_v6 as M,descent_bn254_v6 as DESC
E=Q.sumprod
def _mv6(Mat,v):A=Mat;return[E(A[0],v),E(A[1],v),E(A[2],v),E(A[3],v),E(A[4],v),E(A[5],v)]
def Cx(A,B):return[[D(A[e][C]*B[C][f]for C in C(6))for f in C(6)]for e in C(6)]
def Bc(G,v):return E(v,(E(G[0],v),E(G[1],v),E(G[2],v),E(G[3],v),E(G[4],v),E(G[5],v)))
def _det_int(Mat):
	E=A(Mat);B=[e[:]for e in Mat];G=J=1
	for D in C(E-1):
		if B[D][D]==0:
			F=AY((A for A in C(D+1,E)if B[A][D]),-1)
			if F<0:return 0
			B[D],B[F]=B[F],B[D];G=-G
		F=B[D][D]
		for H in C(D+1,E):
			for I in C(D+1,E):B[H][I]=(B[H][I]*F-B[H][D]*B[D][I])//J
		J=F
	return G*B[E-1][E-1]
def _eval_poly_mod(cs,s,p):
	A=0
	for B in Ap(cs):A=(A*s+B)%p
	return A
def _eval_zk_mod(vec,s,p):
	B=X(3,-1,p);A=[1]
	for E in C(5):A.append(A[-1]*s%p)
	return D(vec[B]*(D(M.ZK_NUM[B][C]*A[C]for C in C(6))%p)for B in C(6)if vec[B])*B%p
def _is_prime(n):return n>=2 and j(n%A for A in C(2,P(n**.5)+1))
def verify_swar_sb_u7_and_stageB():
	Bf=M.load_data();s=Bf[AP];Du=Bf[Be]['kap'];Cy,Cz,AR,Dv,Bg=s['GA'],s['GB'],s['SMAT'],s['FAC6'],s['UMATS'];Bh=[[o(aq,ar)for(aq,ar)in aa]for aa in s[AQ]];C_=[o(ab,ac)for(ab,ac)in s['EMB_W']];Bi=Cx(AR,AR);b=random.Random(20261120);AS=0;D0=0x8000800080008000;Bj=M.rand_below;Bk=b.getrandbits
	for k in C(2000):
		D1=b.random()<G/9.;Bl=48 if D1 else 0;Bm=0;Bn=I
		for Dw in C(16):
			Bo=[3520+Bj(Bk,2401)for aD in C(4)];Bp=[Bj(Bk,8481)for aE in C(4)];D2=[32768-aF+Bl for aF in Bo];A4=[D2[aG]+Bp[aG]for aG in C(4)]
			if W(A>=65536 for A in A4):AS+=1
			D3=A4[0]|A4[1]<<16|A4[2]<<32|A4[3]<<48;Bm|=D3
			if W(Bp[A]+Bl>=Bo[A]for A in C(4)):Bn=L
		if(Bm&D0!=0)!=Bn:AS+=1
	AT=K
	for k in C(500):
		Bq=[b.randint(-450,450)for at in C(6)];Br=[b.randint(-750,750)for au in C(6)]
		if not W(Bq)or not W(Br):continue
		Bs=Bt=K
		for z in C(3):Bu=D(Bh[z][A]*Bq[A]for A in C(6));Bv=D(Bh[z][A]*Br[A]for A in C(6));An=M.U_BN254*C_[z]*(Bu*Bu);Bw=An-Bv*Bv;Bs+=log2(Bw.real**2+Bw.imag**2);Bt+=log2(An.real**2+An.imag**2)
		AT=B(AT,J(Bs-Bt))
	Bx=M._primes_upto(16383);Ao=Bx[1:];V=(1<<64)-1
	if M._P17_G_CACHE is S:
		D4=M._primes_upto((1<<17)-1);m=[]
		for E in D4:
			if E>=16384:
				if M.U_BN254%E==0:m.append(E);continue
				if E&7 not in(1,3):continue
				Ap=M._ppow_h([0,1,0,0,0,0],E,E);Ap[1]=(Ap[1]-1)%E;A5=M._pgcd(Ap,M.H_COEFFS,E)
				if A(A5)>1:
					By=M._pmod([-aK%E for aK in M.G_LEAD_COEFFS],A5,E)
					if A(M._pgcd(By,A5,E))>1:m.append(E);continue
					p=M._ppowmod(By,(E-1)//2,A5,E);p[0]=(p[0]-1)%E
					if A(M._pgcd(p,A5,E))>1:m.append(E)
		M._P17_G_CACHE=m
	m=M._P17_G_CACHE;A7=l(Ao+m,reverse=L);AU,Aq=[],[];D5=(1<<64)-1
	for E in A7:
		Bz=I
		for(B_,C0)in g(AU):
			if C0*E<=D5:AU[B_]=C0*E;Aq[B_].append(E);Bz=L;break
		if not Bz:AU.append(E);Aq.append([E])
	AV=A6(c(AU,Aq));C1=[]
	for(C2,D6)in AV:n=C2<<64-C2.bit_length();Ar=((1<<128)-1)//n-(1<<64);C1.append((n,Ar,[X(1<<64,aI,n)for aI in C(1,8)],[(aJ,X(aJ,-1,1<<64),V//aJ)for aJ in D6]))
	As=[2]+A7;AW=0;AX=0;D7=Q.sumprod
	for k in C(20):
		t=b.getrandbits(390)|1<<389
		for k in C(3):t*=As[b.randrange(A(As))]
		if m:t*=m[b.randrange(A(m))]
		C3=Ao[b.randrange(A(Ao))]
		if t.bit_length()+3*C3.bit_length()<=511:t*=C3**3
		At=1;C4=(t&-t).bit_length()-1;D8=t>>C4;At<<=C4;A8=[D8>>64*av&V for av in C(8)];Au=30
		for(n,Ar,AY,C5)in C1:
			Au+=66+4*A(C5);Av=A8[0]+D7(A8[1:],AY);A9=Av&V;D9=Av>>64&V;DA=Av>>128;Aw=A9+D9*AY[0]+DA*AY[1];AA=Aw&V;AB=Aw>>64&V;DB=Aw>>128
			if DB:AA+=AY[1];AB+=AA>>64;AA&=V
			if AB>=n:AB-=n
			C6=Ar*AB;Ax=(C6&V)+AA;DC=(C6>>64)+AB+(Ax>>64)+1&V;Ax&=V;x=AA-DC*n&V
			if x>Ax:x=x+n&V
			if x>=n:x=x-n&V
			for(E,C7,DD)in C5:
				if x*C7&V<=DD:
					while L:
						Au+=50;C8=[0]*8;Ay=0
						for AC in C(8):C9=A8[AC]-Ay;Az=1 if C9<0 else 0;CA=(C9&V)*C7&V;C8[AC]=CA;Ay=(CA*E>>64)+Az
						if Ay==0:A8=C8;At*=E
						else:break
		AX=B(AX,Au);DE=D(A8[A]<<64*A for A in C(8));A_=t;CB=1
		for E in As:
			while A_%E==0:A_//=E;CB*=E
		if At!=CB or DE!=A_:AW+=1
	CC,DF,Bg,AR,Bi,CD=M.build_exact_int_geom(s);B0=Bc
	def DG(a):
		if O(-A for A in a)<O(a):return I
		for D in C(6):
			E=B0(DF[D],a)
			if E<0:return I
			if E==0:
				A=_mv6(Bg[D],a);F=[-aP for aP in A]
				if O(F if O(F)<O(A)else A)<O(a):return I
		for H in(AR,Bi):
			B=_mv6(H,a);G=[-aO for aO in B]
			if O(G if O(G)<O(B)else B)<O(a):return I
		return L
	A0=[]
	for E in C(101,1800,2):
		if not _is_prime(E):continue
		for h in C(E):
			if _eval_poly_mod(M.H_COEFFS,h,E)==0:
				CE=-((M.U_BN254+1)%E*_eval_poly_mod(M.W_COEFFS,h,E)+2)%E
				if CE==0:continue
				CF=X(CE,-1,E)
				if X(CF,(E-1)//2,E)==1:
					for B1 in C(1,(E+1)//2):
						if B1*B1%E==CF:A0.append((E,h,B1));break
		if A(A0)>=24:break
	DH=36*CD;DI=P(H(11e1/(Cz[0][0]/Cy[0][0])*CD));Z=[]
	for DJ in C(-2,3):
		for DK in C(-2,3):
			for DL in C(-2,3):
				for DM in C(-2,3):
					for DN in C(-2,3):
						for AZ in C(-2,3):
							T=[AZ,DN,DM,DL,DK,DJ]
							if W(T)and(B2:=B0(CC,T))<=DH and DG(T):Z.append((B2,T))
	Z.sort(key=lambda x:x[0]);Z=Z[:80];u,AD=[],[]
	for DO in C(1,3):
		for DP in C(-2,3):
			for DQ in C(-2,3):
				for DR in C(-2,3):
					for DS in C(-2,3):
						AE=[(aL,B3,B2)for aL in C(-6,7)if(B2:=B0(CC,(B3:=[aL,DS,DR,DQ,DP,DO])))<=DI]
						if AE:
							u.append((A(AE),[aM[1]for aM in AE],[aN[2]for aN in AE]))
							for(k,B3,DT)in AE:AD.append((B3,DT,A(AD)))
						if A(u)>=50:break
					if A(u)>=50:break
				if A(u)>=50:break
			if A(u)>=50:break
		if A(u)>=50:break
	v={(ad,ae):0 for ad in C(A(Z))for ae in C(A(AD))};Aa=0
	for(E,h,AF)in A0:
		DU=[_eval_zk_mod(aw,h,E)for(ax,aw)in Z];DV=[_eval_zk_mod(ay,h,E)for(ay,az,az)in AD];B4={}
		for(N,DW)in g(DU):B4.setdefault(DW,[]).append(N)
		for(Ab,DX)in g(DV):
			AG=AF*DX%E;DY=(E-AG)%E
			for N in B4.get(AG,()):v[N,Ab]+=1;Aa+=1
			for N in B4.get(DY,()):v[N,Ab]+=1;Aa+=1
	a=4;CG={af:0 for af in v};CH={O(ag):ah for(ag,ai,ah)in AD};CI=[aj for(aj,ak)in Z];CJ=A(set(CI))<A(CI);CK=B5=B6=0;B7=S;B8=L
	for(E,h,AF)in A0:
		B9=X(AF,-1,E);w=(E+a-1)//a;d=A(Z);Ac=[0]*w;CL=[0]*d;CM=[0]*d
		for(N,(k,T))in g(Z):e=_eval_zk_mod(T,h,E)*B9%E;F=e//a;Ac[F]+=1;CL[N]=F;CM[N]=N<<8|e&a-1
		Ad=[0]*(w+1);BA=[0]*w;Ae=0;Af=Bd((w+15)//8+1);CN=0
		for F in C(w):
			Ad[F]=Ae;BA[F]=Ae
			if Ac[F]:Af[F>>3]|=1<<(F&7);CN+=1
			Ae+=Ac[F];Ac[F]=0
		if CN>R(w,d):B8=I
		Ad[w]=Ae;CO=[0]*d
		if B7 is S:B7=bytes(Af)
		for N in C(d):F=CL[N];CP=BA[F];CO[CP]=CM[N];BA[F]=CP+1
		for(AH,AI,k)in u:
			A1=_eval_zk_mod(AI[0],h,E)
			for BB in(1,-1):
				BC=A1 if BB==1 else-(A1+AH-1)%E;BD=(A1+AH-1)%E if BB==1 else-A1%E;DZ=[(BC//a,BD//a)]if BC<=BD else[(BC//a,w-1),(0,BD//a)]
				for(Ag,Da)in DZ:
					AJ=Da-Ag+1;CK+=AJ;assert AJ<=8,f"1-byte occ_bm overflow: n_span={AJ} > 8";CQ=Ag>>3;Db=Ag&7;Dc=Af[CQ]|Af[CQ+1]<<8;CR=Dc>>Db&(1<<AJ)-1
					if CR==0:continue
					for z in C(AJ):
						if CR>>z&1:
							B5+=1;B6+=1;F=Ag+z
							for Dd in C(Ad[F],Ad[F+1]):
								CS=CO[Dd];e=F*a|CS&a-1;N=CS>>8;p=(e-A1)%E if BB==1 else(-e-A1)%E
								if p<AH:CG[N,CH[O(AI[p])]]+=1
	CT=D(1 for A in v if CG[A]!=v[A]);BE=BF=BG=0;AK=(1<<64)-1;y=835
	for U in((1<<44)+7,(1<<55)-55):
		for q in(0,1,255,256,U-1,U-y,U-y+5,U-100):
			for AL in(1,-1):
				BH=q if AL==1 else-(q+y-1)%U;A2=(q+y-1)%U if AL==1 else-q%U;De=A6(C(BH>>8,(A2>>8)+1))if BH<=A2 else A6(C(BH>>8,U+255>>8))+A6(C(0,(A2>>8)+1))
				for BI in De:
					AM=BI<<8;CU=((AM-q if BI>=q>>8 else U-q+AM)if AL==1 else A2-AM if BI<=A2>>8 else A2+U-AM)&AK
					for BJ in(0,1,7,127,128,254,255):
						BK=AM+BJ
						if BK>=U:continue
						BL=(BK-q)%U if AL==1 else(-BK-q)%U;CV=BJ+CU&AK if AL==1 else CU-BJ&AK
						if(CV<y)!=(BL<y)or BL<y and CV!=BL:BE+=1
		BM=[b.randrange(U)for aA in C(6)];Df=[aB&268435455 for aB in BM];Dg=[aC>>28 for aC in BM];Dh=(1<<64)%U
		for k in C(32):
			T=[b.randrange(0,8192)for aH in C(6)];CW=D(T[A]*Df[A]for A in C(6));BN=D(T[A]*Dg[A]for A in C(6));CX=BN>>36;Di=(BN<<28)+CW&AK;Dj=(BN<<28&AK)+CW>>64;CX+=Dj;Az=(Di+CX*Dh)%U
			if Az!=D(T[A]*BM[A]for A in C(6))%U:BG+=1
	CY=[];CZ={}
	for(N,(Dk,T))in g(Z):CZ.setdefault(O(T[1:]),[]).append((T[0],N,T))
	for Ah in CZ.values():Ah.sort();CY.append(Ah)
	r,A9,Dl=A0[0];BO=X(Dl,-1,r);Dm=_eval_zk_mod([1,0,0,0,0,0],A9,r)*BO%r;i=(r+a-1)//a;d=A(Z);Ca=[0]*i;Cb=[0]*d;Cc=[0]*d
	for(N,(Dk,T))in g(Z):e=_eval_zk_mod(T,A9,r)*BO%r;F=e//a;Ca[F]+=1;Cb[N]=F;Cc[N]=N<<8|e&a-1
	AN=[0]*(i+1);BP=[0]*i;f=0
	for F in C(i):AN[F]=f;BP[F]=f;f+=Ca[F]
	AN[i]=f;Cd=[0]*d
	for N in C(d):F=Cb[N];Cd[BP[F]]=Cc[N];BP[F]+=1
	Ai=[0]*i;Ce=[0]*d;Cf=[0]*d;BQ=0
	for Ah in CY:
		BR=S;Aj=0
		for(AZ,N,T)in Ah:Aj=(Aj+Dm)%r if BR is not S and AZ==BR+1 else _eval_zk_mod(T,A9,r)*BO%r;BR=AZ;F=Aj//a;Ai[F]+=1;Ce[BQ]=F;Cf[BQ]=N<<8|Aj&a-1;BQ+=1
	AO=[0]*(i+1);BS=[0]*i;f=0;Cg=Bd((i+15)//8+1)
	for F in C(i):
		AO[F]=f;BS[F]=f
		if Ai[F]:Cg[F>>3]|=1<<(F&7)
		f+=Ai[F];Ai[F]=0
	AO[i]=f;Ch=[0]*d
	for AC in C(d):F=Ce[AC];Ch[BS[F]]=Cf[AC];BS[F]+=1
	BT=bytes(Cg)==B7;BU=AO==AN and j(l(Ch[AO[A]:AO[A+1]])==l(Cd[AN[A]:AN[A+1]])for A in C(i))and BT
	if not BU:BF=1
	Ci={al:0 for al in v}
	for(E,h,AF)in A0:
		B9=X(AF,-1,E);Cj=[0]*E;Ck=[0]*A(Z)
		for(N,(k,T))in g(Z):e=_eval_zk_mod(T,h,E)*B9%E;Cj[e]+=1;Ck[N]=e
		Ak=[0]*(E+1);BV=[0]*E;f=0
		for BW in C(E):Ak[BW]=f;BV[BW]=f;f+=Cj[BW]
		Ak[E]=f;Cl=[0]*A(Z)
		for N in C(A(Z)):e=Ck[N];Cl[BV[e]]=N;BV[e]+=1
		for(AH,AI,k)in u:
			for p in C(AH):
				AG=_eval_zk_mod(AI[p],h,E);Ab=CH[O(AI[p])]
				for Cm in(AG,-AG%E):
					for Dn in C(Ak[Cm],Ak[Cm+1]):Ci[Cl[Dn],Ab]+=1
	Cn=D(1 for A in v if Ci[A]!=v[A]);BX=0;Al=floor(468/14);Am=[(0,Al+1),(883.,960+Al+1),(1767.,1920+Al+1)];Co=j(Am[A][1]<Am[A+1][0]for A in C(A(Am)-1))
	for Do in(0,1,2):
		for k in C(50):
			BY=D(b.uniform(55.25,6e1)for A in C(Do))
			if BY>14e1:continue
			Cp=b.uniform(K,B(.5,R(11e1,168.-BY)));Cq=468.-BY-Cp;Cr=b.randint(8,R(Al,P(Cq/14.)));Dp=[Cq/Cr]*Cr;Dq=D(floor(16.*A)for A in Dp);Cs=floor(7488.)-Dq-floor(16.*Cp)
			if Cs<=2240:
				if not W(A<=Cs<=B for(A,B)in Am):BX+=1
	BZ=0
	for k in C(1000):
		Ct=b.randint(6400,8000);Cu=Ct-2720;Dr=32768-Cu;Cv=b.randint(0,Cu+200);Ds=Cv+48+Dr&32768!=0;Dt=Ct-48-Cv<=2720
		if Ds!=Dt:BZ+=1
	A3=P(2.**45.352);Ba=[(am,an)for am in(0,400,795)for an in(0,1,(1<<32)-1,1<<32,(1<<32)+1,(1<<44)+12345,A3-1)];Bb=[ao*A3+ap for(ao,ap)in Ba];Cw=A(set(Bb))==A(Ba)and j(0<=A<1<<56 for A in Bb)and j((A//A3,A%A3)==(B,C)for(A,(B,C))in c(Bb,Ba))and 0*A3+0!=0*A3+(1<<32);return Y(ok=CT==0 and B8 and BE==0 and BF==0 and BT and BG==0 and BU and Cn==0 and CJ and AS==0 and AW==0 and Co and BX==0 and BZ==0 and Cw and AT<=1e-12 and Aa>=1400,zero_float_geom=L,det_occ_bound_ok=B8,flat_array_zero_dict_build=L,flat_sbu7_exact_ok=BT,zero_sort_exact_ok=BU,exact_slot_cutoff_M_t=L,packed_8b_u64_wrap_errors=BE,band2_a0_line_build_errors=BF,split28_6coord_res_errors=BG,band0_block_direct_apply_mismatches=Cn,qa_stored_in_record=I,qa_ties_tested=CJ,scalar_swar_10k_errors=AS,p2_shift_threshold_errors=BZ,a_slot_64bit_linear_ok=Cw,max_ng_invariance_bits=AT,p14_primes=A(Bx),p17_g_split_primes=A(m),p17_strip_odd_primes=A(A7),p17_word_chunks=A(AV),p14_zero_gcd_gm_ok=AW==0 and AX<=M.C_SB and A(AV)*66+A(A7)*4+30+1600<=M.C_SB,p14_strip_ops_max=AX,p14_strip_ops_e16=A(AV)*66+A(A7)*4+30+800,p14_strip_errors=AW,stageB_window_disjoint=Co,stageB_window_false_neg=BX,sb_u7_hits=Aa,sb_u7_mismatches=CT,sb_total=CK,sb_occ=B5,sb_ent=B6,sb_occ_skip_pct=H(1e2*(G-B6/B(1,B5)),2))
def verify_12d_descent_lattices():
	AG=M.load_data();AH=AG[AP];f=M.load_witnesses();AI=[[o(am,an)for(am,an)in aa]for aa in AH[AQ]];BQ=[o(ab,ac)for(ab,ac)in AH['EMB_W']];q={ad:[D(ae[A]*AI[ao][1]**A for A in C(6))for ao in C(3)]for(ad,ae)in M.F_REC_COEFFS.items()};AJ,Br=6.*log2(M.U_BN254),6.*log2(M.P_BN254);Ak=[[1,0,1,-1,-3,-3],[0,1,-1,0,-1,-1],[0,0,-1,1,1,3],[0,0,-1,-1,0,2],[0,0,1,0,1,-1],[0,0,0,1,1,2]];AK=_det_int
	def d(a,b):
		A=[0]*11
		for D in C(6):
			if a[D]:
				for B in C(6):A[D+B]+=a[D]*b[B]
		for E in C(10,5,-1):
			F=A[E]
			if F:
				A[E]=0
				for B in C(6):A[E-6+B]-=F*M.H_COEFFS[B]
		return A[:6]
	def Al(poly_num,denom):
		B=denom;D,A=[],poly_num[:]
		for I in C(6):
			D.append(A)
			if I<5:
				E=[0]+A[:5];F=A[5]
				if F:
					for G in C(6):E[G]-=F*M.H_COEFFS[G]
				A=E
		H=J(AK([[D[a6][a5]for a6 in C(6)]for a5 in C(6)]));assert H%B**6==0;return H//B**6
	def BR(a,b):
		G=[D(a[A]*M.ZK_NUM[A][a7]for A in C(6))for a7 in C(6)];H=[D(b[A]*M.ZK_NUM[A][a8]for A in C(6))for a8 in C(6)];A,B=d(G,G),d(H,H);L=d(M.W_COEFFS,A);N=Al([(M.U_BN254+1)*L[a9]+2*A[a9]+B[a9]for a9 in C(6)],9);E=[[1,0,0,0,0,0],A];F=[[1,0,0,0,0,0],B]
		for Q in C(3):E.append(d(E[-1],A));F.append(d(F[-1],B))
		I=[0]*6
		for(J,O)in M.F_REC_COEFFS.items():
			P=d(d(O,E[J//2]),F[4-J//2])
			for K in C(6):I[K]+=P[K]
		return Al(I,3**8),N
	def Am(vec,mod):A=X(3,-1,mod);return[D(vec[A]*M.ZK_NUM[A][a_]for A in C(6))*A%mod for a_ in C(6)]
	An=O(O(A.real for A in A)for A in AI);Ao=O(O(A.imag for A in A)for A in AI);g=Q.sumprod
	def AL(a,b,need_f=L):
		N=O=K;U=B(B(J(A)for A in a),B(J(A)for A in b),1);L=B(0,U.bit_length()-40);P=h(1<<L);Q=[ba/P for ba in a];R=[bb/P for bb in b]
		for A in C(3):
			S=o(g(An[A],Q),g(Ao[A],Q));T=o(g(An[A],R),g(Ao[A],R));D,E=S*S,T*T;F=((M.U_BN254+1)*BQ[A]+2.)*D+E;O+=log2(F.real*F.real+F.imag*F.imag)+4.*L
			if need_f:G,H=D*D,E*E;I=q[8][A]*G*G+q[6][A]*G*D*E+q[4][A]*G*H+q[2][A]*D*H*E+q[0][A]*H*H;N+=log2(I.real*I.real+I.imag*I.imag)+16.*L
		return N,O
	BS=log2(AG[Aj]['cells']['4']);i=[Q.factorial(af)for af in C(13)];AM=[]
	for r in C(2):
		for s in C(2):
			for t in C(5):
				Ap=16-16*r-9*s-4*t
				if Ap<0:continue
				for AN in C(R(Ap,12-r-s-t)+1):
					AO=AN+t+s+r
					if AO==0:continue
					AR=i[12]//(i[12-AO]*i[AN]*i[t]*i[s]*i[r])*(1<<AO);AM.append((AR,(AN,t,s,r)))
	Aq,A1=[],0
	for(AR,_)in AM:A1+=AR;Aq.append(A1)
	Ar=A1//2==AG[Aj]['cells']['4']
	def BT(rng_loc):A=rng_loc;G=AX.bisect_left(Aq,A.randrange(A1)+1);B,C,D,E=AM[G][1];F=[1]*B+[2]*C+[3]*D+[4]*E+[0]*(12-B-C-D-E);A.shuffle(F);return[bc*(1 if A.getrandbits(1)else-1)if bc else 0 for bc in F]
	As=M.build_exact_int_geom(AH)[0];BU=random.Random(42);At=random.Random(20261121);A2=AS=L;A3=K;A4=0;k=0;Au,Av,AT,AU,a,b=[],[],[],[],[],[]
	for e in f[AE]:
		l,F,AV,Aw=e[0],M.b62d(e[1]),M.b62d(e[2]),M.b62d(e[3]);AW=F.bit_length();Au.append(AW);Av.append(l);BV=X(3,-1,F);AY=[1]
		for _ in C(5):AY.append(AY[-1]*AV%F)
		Ax=[D(M.ZK_NUM[ap][A]*AY[A]for A in C(6))*BV%F for ap in C(6)];V=Ax+[-Aw*aB%F for aB in Ax];Ay=X(V[0],-1,F);AZ=[[F if aV==0 else 0 for aV in C(12)]]+[[-V[aC]*Ay%F-(F if-V[aC]*Ay%F>F//2 else 0)if aW==0 else 1 if aW==aC else 0 for aW in C(12)]for aC in C(1,12)];DESC.bareiss_lll_12d(AZ,As,1,1,I,11,40,(9,10));N=[[AZ[aD][aq]for aD in C(12)]for aq in C(12)]
		if J(AK(N))!=F:A2=I
		for Z in C(12):
			m=[N[aE][Z]for aE in C(6)];A5=[N[6+aF][Z]for aF in C(6)]
			if(_eval_zk_mod(m,AV,F)-Aw*_eval_zk_mod(A5,AV,F))%F!=0:A2=I
		A6=[N[ar][0]+N[ar][1]for ar in C(12)];Az,A_=BR(A6[:6],A6[6:])
		if Az%F!=0 if l=='f'else A_%F!=0:AS=I
		BW,BX=AL(A6[:6],A6[6:],L);A3=B(A3,J(log2(Az)-BW),J(log2(A_)-BX));A7=(log2(F)+BS)/12.;A8=AW if l=='f'else K;A9=AW if l=='g'else K;BY=D(A*M.rho((48.*A7+B-A8)/59.)*M.rho((12.*A7+AJ+C-A9)/59.)for(A,B,C)in DESC.QBINS);Aa=[];B0=[]
		for _ in C(60):
			u=BT(BU);U=[g(N[aG],u)for aG in C(12)];S=B(J(A).bit_length()for A in U);k=B(k,S)
			if S>=48:A4+=1
			V,v=AL(U[:6],U[6:],L);B0.append((V,v));a.append(V-48.*A7);b.append(v-12.*A7-AJ);Aa.append(M._rho_fast((V-A8)/59.)*M._rho_fast((v-A9)/59.))
		DESC._SQ_NORMS_CACHE[l,F]=B0;DESC._SQ_BASIS_CACHE[l,F]=AZ;AT.append(H(D(Aa)/A(Aa)/B(BY,1e-300),2))
		for B1 in(16,64,128):
			BZ=log2(.5*(pi**6/72e1)*B1**12);B2=(log2(F)+BZ)/12.;Ba=D(A*M._rho_fast((48.*B2+B-A8)/59.)*M._rho_fast((12.*B2+AJ+C-A9)/59.)for(A,B,C)in DESC.QBINS);Ab=[]
			for _ in C(32):
				B3=[At.gauss(K,G)for aX in C(12)];Bb=sqrt(D(A*A for A in B3));Bc=B1*At.random()**(G/12.);u=[P(H(aY*Bc/Bb))for aY in B3]
				if not W(u):u[0]=1
				U=[g(N[aZ],u)for aZ in C(12)];S=B(J(A).bit_length()for A in U);k=B(k,S)
				if S>=48:A4+=1
				V,v=AL(U[:6],U[6:],L);Ab.append(M._rho_fast((V-A8)/59.)*M._rho_fast((v-A9)/59.))
			AU.append(D(Ab)/A(Ab)/B(Ba,1e-300))
	a.sort();b.sort()
	def AA(arr,qv):C=arr;return H(C[R(A(C)-1,B(0,ceil(qv*A(C))-1))],2)
	B4=j(a[R(A(a)-1,floor(B/6.*A(a)))]>=DESC.QBINS_LO[B][0]and b[R(A(b)-1,floor(B/6.*A(b)))]>=DESC.QBINS_LO[B][1]for B in C(6));B5=[AA(a,(ag+1)/6.)for ag in C(6)];B6=[AA(b,(ah+1)/6.)for ah in C(6)];B7=j(B5[A]<=DESC.QBINS[A][1]and B6[A]<=DESC.QBINS[A][2]for A in C(6));B8=.1,.25,.5,.75,.9,G;n=Ac=L;Ad=0;AB=0;T=[0,0];E=M.P_BN254;Bd=X(3,-1,E);Be=[-ai%E for ai in M.G_LEAD_COEFFS];B9=M._pinv_h(Be,E)
	def Bf(a12,b12):B=b12;A=a12;C=[(bd+be)%E for(bd,be)in c(M._pmul_h(A[0],B[0],E),M._pmul_h(M._pmul_h(A[1],B[1],E),B9,E))];D=[(bf+bg)%E for(bf,bg)in c(M._pmul_h(A[0],B[1],E),M._pmul_h(A[1],B[0],E))];return C,D
	def BA(a):
		D,F,G,H,I,J=a;R=D*D;P=2*D*F;Q=2*D*G+F*F;O=2*(D*H+F*G);M=2*(D*I+F*H)+G*G;N=2*(D*J+F*I+G*H);A=2*(F*J+G*I)+H*H;B=2*(G*J+H*I);C=2*H*J+I*I;K=2*I*J;L=J*J
		if L:K+=2*L;C-=L;B-=2*L;N+=4*L;M-=3*L
		if K:C+=2*K;B-=K;A-=2*K;M+=4*K;O-=3*K
		if C:B+=2*C;A-=C;N-=2*C;O+=4*C;Q-=3*C
		if B:A+=2*B;N-=B;M-=2*B;Q+=4*B;P-=3*B
		if A:N+=2*A;M-=A;O-=2*A;P+=4*A;R-=3*A
		return[R%E,P%E,Q%E,O%E,M%E,N%E]
	def Bg(a12):A=a12;B=[(bh+bi)%E for(bh,bi)in c(BA(A[0]),M._pmul_h(BA(A[1]),B9,E))];C=[2*bj%E for bj in M._pmul_h(A[0],A[1],E)];return B,C
	def Bh(base12,exp_v):
		A=[1,0,0,0,0,0],[0]*6
		for B in bin(exp_v)[2:]:
			A=Bg(A)
			if B=='1':A=Bf(A,base12)
		return A
	for e in f[AF]:
		AC=[M.b62d(at)for at in e[0].split()];Bi=[M.b62d(au)for au in e[1].split()];BB=[Bi[av*12:(av+1)*12]for av in C(6)];Bj=Bh(([1,0,0,0,0,0],[a0%E for a0 in AC]),M.R_BN254)
		if Bj[1]!=[0]*6 or not W(AC):Ac=I
		BC=[[0]*12 for aw in C(6)]
		for Z in C(12):
			m=[BB[aH][Z]for aH in C(6)];Bk=[D(m[A]*M.ZK_NUM[A][aI]for A in C(6))*Bd%E for aI in C(6)];Bl=M._pmul_h(AC,Bk,E);Ae=[3*aJ%E for aJ in Bl];Ae=[aK-E if aK>E//2 else aK for aK in Ae]
			for BD in C(6):
				BE=D(Ak[BD][A]*Ae[A]for A in C(6))
				if BE%3!=0:n=I
				BC[BD][Z]=BE//3
		N=BB+BC
		if J(AK(N))!=E**6:n=I
		if not DESC.verify_bareiss_910_cert([[N[a4][a1]for a4 in C(12)]for a1 in C(12)],As,0xa8ccdee3ccdbd36d,1):n=I
		for Z in C(12):
			m=[N[aL][Z]for aL in C(6)];A5=[N[6+aM][Z]for aM in C(6)]
			if M._pmul_h(AC,Am(m,E),E)!=Am(A5,E):n=I
			w=B(J(A).bit_length()for A in m);S=B(J(A).bit_length()for A in A5);T[0]=B(T[0],w);T[1]=B(T[1],S)
			if w>=128 or S>=160:AB+=1
		Af=[[N[aN][ax]for aN in C(12)]for ax in C(12)]
		for U in[Af[aO]+[]for aO in C(12)]+[[a2+aR*a3 for(a2,a3)in c(Af[aP],Af[aQ])]for aP in C(12)for aQ in C(aP+1,12)for aR in(1,-1)]:
			w=B(J(A).bit_length()for A in U[:6]);S=B(J(A).bit_length()for A in U[6:]);T[0]=B(T[0],w);T[1]=B(T[1],S)
			if w>=128 or S>=160:AB+=1
		Ad+=1
	from fractions import Fraction as x;BF=[[x(ay,3)for ay in aj]for aj in M.ZK_NUM]
	def Bm(a,b):
		A=[x(0)]*11
		for E in C(6):
			for B in C(6):A[E+B]+=a[E]*b[B]
		for D in C(10,5,-1):
			F=A[D]
			if F:
				for B in C(6):A[D-6+B]-=F*M.H_COEFFS[B]
				A[D]=x(0)
		return A[:6]
	def Ag(p):return[P(D(Ak[bk][A]*(3*p[A])for A in C(6)))//3 for bk in C(6)]
	AD=[[Ag(Bm(BF[ak],BF[az]))for az in C(6)]for ak in C(6)];Ah=B(D(J(AD[B][D][A])for B in C(6)for D in C(6))for A in C(6));p=sqrt(B(D(D(J(AD[C][A][B])for C in C(6))**2 for B in C(6))for A in C(6)));Bn=Ag([x(aA)for aA in M.W_COEFFS]);Bo={al:Ag([x(aS)for aS in M.F_REC_COEFFS[al]])for al in(0,2,4,6,8)};Bp=D(B(D(J(D(Bo[E][A]*AD[A][F][B]for A in C(6)))for F in C(6))for B in C(6))for E in(0,2,4,6,8));BG=B(D(J(D(Bn[A]*AD[A][E][B]for A in C(6)))for E in C(6))for B in C(6));BH=log2(Bp)+7.*log2(Ah)+384.;BI=log2((M.U_BN254+1)*BG+3.)+log2(Ah)+96.;BJ=log2(((M.U_BN254+1)*BG+2.)*2.**256+2.**320)+log2(Ah);BK=H(6.*(log2(p)+BH),2);BL=H(6.*(log2(p)+BI),2);BM=H(6.*(log2(p)+BJ),2)
	def Ai(rb):return D((5-A)**2*(24*ceil((A+1)*rb/64.)**2+16*ceil((A+1)*rb/64.)*(ceil(A*rb/64.)if A else 1)+60*ceil((A+1)*rb/64.))for A in C(5))
	y=ceil(BK/64.);z=ceil(BL/64.);A0=ceil(BM/64.);BN=B(Ai(log2(p)+BH)+Ai(log2(p)+BI)+450000,Ai(log2(p)+BJ)+150000);BO=B(4*(y+z)+4170*(y+z)+18990+41*y*(6*y+10)+41*z*(6*z+10),4*A0+4170*A0+9495+41*A0*(6*A0+10));BP=A4==0 and AB==0 and k<48 and T[0]<128 and T[1]<160 and y<=40 and z<=16 and A0<=32 and BN<DESC.NORMS and BO<DESC.TDIV;Bq=Ar and AS and A3<1e-08 and j(A>=G for A in AT)and B4 and B7 and R(AU)>=G and Ad==24 and A(f[AE])==24 and A(f[AF])==24;return Y(ok=A2 and n and Ac and Bq and BP,prenorm_gate_uniform_ok=BP,sq_gate_max_cb=k,bt_gate_max_ab=T,sq_gate_fails=A4,bt_gate_fails=AB,sq_nf_hi_bits=BK,sq_ng_hi_bits=BL,bt_ng_hi_bits=BM,gate_norms_wc=BN,gate_tdiv_wc=BO,all_24_sq_reduced_live_by_bareiss_910=A2,all_24_boot_bareiss_910_certified=n,boot_gt_subgroup_ok=Ac,cells4_exact_cardinality_ok=Ar,special_q_count=A(f[AE]),special_q_sides=Av,special_q_bits=Au,exact_side_divisibility_ok=AS,boot_count=A(f[AF]),boot_candidates_gated=144*Ad,exact_int_vs_fp_max_err_bits=H(A3,12),paired_joint_smooth_ratios_r4=AT,hi_radius_R16_128_min_paired_ratio=H(R(AU),2),cf_sextiles_obs=B5,cg_sextiles_obs=B6,cf_quantiles_p10_25_50_75_90_max=[AA(a,aT)for aT in B8],cg_quantiles_p10_25_50_75_90_max=[AA(b,aU)for aU in B8],qbins_hi_6_of_6_ok=B7,qbins_lo_6_of_6_ok=B4)
def d(log2U,Lf=62.5,Lg=62.,Nf=390.64,Ng=466.6,excess=.687,t_dense=.8):
	A=t_dense;J=excess*2.**log2U;E=R(A*R(Lf,Lg),B(8.,log2(J/25.)+log2(3.)));F=K
	for(L,N)in((Nf/Lf,Lf),(Ng/Lg,Lg)):
		D=E/N
		if D>=A:continue
		G=(A-D)/2e2
		for O in C(200):I=D+G*(O+.5);F+=M.kappa_rho(L,I)*G/I
	return H(E,4),H(F,6)
def verify_live_filtsim_z6(log2U=13,seed=20261120,W_target=60,check_kernel=I,export_hx=I,texc_abs=S):
	BO=texc_abs;BN=export_hx;BM=W_target;A7=check_kernel;o=random.Random(seed);BP=Q.log(2.);p,q,CL,CM,Aq=62.5,62.,390.64,466.6,.687;Ar=(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1);BQ={(1,0):(1,0),(0,1):(1,-1),(-1,1):(0,-1),(-1,0):(-1,0),(0,-1):(-1,1),(1,-1):(0,1)}
	def AZ(x,y):return x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1]
	def Aa(a,b):
		if a>=0 and b>=0 or a<=0 and b<=0:return a,b,0
		if J(a)<=J(b):return 0,a+b,-a
		return a+b,0,b
	def r(v):A,B,C=Aa(v[0],v[1]);return(A!=0)+(B!=0)+(C!=0)
	def BR(v):A,B,C=Aa(v[0],v[1]);return J(A)+J(B)+J(C)
	def CN(v):A,C,D=Aa(v[0],v[1]);return B(J(A),J(C),J(D))
	def Ab(v):return v in BQ
	As,At=CL/p,CM/q;BS=2.**p/(p*BP);CO=2.**q/(q*BP);Au=BS+CO;CP=BS/Au;k=1<<log2U;BT=P(H(k*CP));CQ=k-BT;d=P(H(Aq*k));AI=.8;Av=R(AI*R(p,q),B(8.,log2(d/25.)+log2(3.)));Aw=M._primes_upto(P(floor(2.**Av)));CR=(A(Aw)+5)//6;A8=k+2*CR;AJ=400;AK=[AI+(G-AI)*(aa+.5)/AJ for aa in C(AJ)];Ac,Ad=[],[];Ax=Ay=K
	for b in AK:Ax+=2.**(p*(b-G))/b;Ay+=2.**(q*(b-G))/b;Ac.append(Ax);Ad.append(Ay)
	Ac=[ab/Ax for ab in Ac];Ad=[ac/Ay for ac in Ad];CS=[M.kappa_rho(As,ad)for ad in AK];CT=[M.kappa_rho(At,ae)for ae in AK];A0=[K]*A8;Ae=[K]*A8;m=0;BU=(G-AI)/AJ
	for _ in C(BT):AL=R(AJ-1,AX.bisect_left(Ac,o.random()));b=AK[AL]+(o.random()-.5)*BU;A0[m]=R(G,Aq*Au*CS[AL]/2.**(p*b)/d);Ae[m]=p*b;m+=1
	for _ in C(CQ):AL=R(AJ-1,AX.bisect_left(Ad,o.random()));b=AK[AL]+(o.random()-.5)*BU;A0[m]=R(G,Aq*Au*CT[AL]/2.**(q*b)/d);Ae[m]=q*b;m+=1
	BV=[]
	for(Az,Af)in((As,p),(At,q)):
		for T in C(0,A(Aw),6):A9=Aw[T];U=R(G,6.*M.kappa_rho(Az,log2(A9)/Af)/A9);A0[m]=U;Ae[m]=log2(A9);BV.append((m,U,A9));m+=1
	AA=BW=K
	for(Az,Af)in((As,p),(At,q)):
		BX=Av/Af;BY=(AI-BX)/2e2
		for T in C(200):
			b=BX+BY*(T+.5);BZ=M.kappa_rho(Az,b)*BY/b;AA+=BZ
			if Af*b<=28.:BW+=BZ
	CU=[af for af in C(k)if A0[af]>=.2];A_=[ag for ag in C(k)if A0[ag]<.2];Ba=[];Ag=K
	for T in A_:Ag+=A0[T];Ba.append(Ag)
	def Bb(lam):
		A=lam
		if A<3e1:
			E=exp(-A);C=G;D=0
			while L:
				C*=o.random()
				if C<=E:return D
				D+=1
		return B(0,P(H(A+sqrt(A)*o.gauss(K,G))))
	e=1000000009;CV=AY(A for B in C(2,50)for A in[X(B,(e-1)//6,e)]if(A*A-A+1)%e==0)
	def u(v):return(v[0]+v[1]*CV)%e
	AB=[o.randrange(1,e)for ap in C(A8)]if A7 else S;E=[];N=[set()for ah in C(A8)];Ah=o.random;B0=CK.partial(M.rand_below,o.getrandbits);CW=AX.bisect_left;Bc=A(A_)-1;CX=Q.log;CY=[(ai,A0[ai])for ai in CU];Bd=[{}for aj in C(d)]
	def B1(d,c,mv):
		if c in d:
			A=d[c];B=A[0]+mv[0],A[1]+mv[1]
			if B==(0,0):del d[c]
			else:d[c]=B
		else:d[c]=mv
	for(F,U,A9)in BV:
		if U<=K:continue
		CZ=G/A9;Be=G/Q.log1p(-U)if U<.999999 else K;O=-1
		while L:
			O=O+1 if Be==K else O+1+P(CX(Ah()or 1e-300)*Be)
			if O>=d:break
			Z=1
			while Ah()<CZ:Z+=1
			Bf=Ar[B0(6)];B1(Bd[O],F,(Z*Bf[0],Z*Bf[1]))
	for O in C(d):
		AC=Bd[O]
		for(F,Ca)in CY:
			if Ah()<Ca:B1(AC,F,Ar[B0(6)])
		for _ in C(Bb(Ag)):Bg=CW(Ba,Ah()*Ag);B1(AC,A_[Bc if Bg>Bc else Bg],Ar[B0(6)])
		Cb=AA+D(r(A)for A in AC.values());Cc=-D(u(B)*AB[A]for(A,B)in AC.items())%e if A7 else 0;E.append([AC,1,AA,AA,Cb,L,Cc])
		for F in AC:N[F].add(O)
	f={ak for ak in C(A8)if N[ak]}
	def B2(nalive,stack=S,seeds=S):
		K=seeds;J=stack;H=nalive;F=[aS for aS in f if A(N[aS])==1]if K is S else A6(K)
		while F:
			C=F.pop();L=N[C]
			if A(L)!=1:continue
			M=AY(iter(L));D=E[M]
			if not D[5]:continue
			B=D[0]
			if J is not S:
				if W(A(N[B])<2 for B in B if B!=C):continue
				if u(B[C])==0:continue
				J.append((C,Y(B),D[6]))
			D[5]=I;H-=1
			for G in B:
				O=N[G];O.discard(M);P=A(O)
				if P==1:F.append(G)
				elif P==0:f.discard(G)
			B.clear()
		return H
	g=B2(d);Ai=A(f);Cd=set(f);AM=[];Bh=0;Aj=B(6 if BN else 2,P(.0001*k))if BO is S else B(2,P(BO));Ce=g-Ai
	while g-Ai>Aj:
		Bh+=1;AD=A6(C(d))
		def Ak(x):
			while AD[x]!=x:AD[x]=AD[AD[x]];x=AD[x]
			return x
		for F in f:
			Bi=N[F]
			if A(Bi)==2:
				Cf,Cg=Bi;Bj,Bk=Ak(Cf),Ak(Cg)
				if Bj!=Bk:AD[Bj]=Bk
		B3={}
		for T in C(d):
			O=E[T]
			if O[5]:Bl=Ak(T);B3[Bl]=B3.get(Bl,K)+O[4]
		Bm=l(B3.items(),key=lambda kv:-kv[1]);Bn=g-Ai-Aj;Ch=B(1,R(Bn//2+1,B(1,A(Bm)//100)));Ci={aq for(aq,ar)in Bm[:Ch]};Bo=I
		for T in C(d):
			O=E[T]
			if O[5]and Ak(T)in Ci:
				A1=O[0]
				if W(A(N[B])<2 for B in A1):continue
				O[5]=I;g-=1;Bo=L;A2=[]
				for a in A1:
					AN=N[a];AN.discard(T);AO=A(AN)
					if AO==1:A2.append(a)
					elif AO==0:f.discard(a)
				A1.clear()
				if A2:g=B2(g,stack=AM,seeds=A2)
		if not Bo:
			Bp=l([(E[aL][4],aL)for aL in C(d)if E[aL][5]and j(A(N[B])>=2 for B in E[aL][0])],reverse=L)
			if not Bp:break
			for(_,T)in Bp[:B(1,Bn//2)]:
				O=E[T]
				if not O[5]:continue
				A1=O[0]
				if W(A(N[B])<2 for B in A1):continue
				O[5]=I;g-=1;A2=[]
				for a in A1:
					AN=N[a];AN.discard(T);AO=A(AN)
					if AO==1:A2.append(a)
					elif AO==0:f.discard(a)
				A1.clear()
				if A2:g=B2(g,stack=AM,seeds=A2)
				if g-A(f)<=Aj:break
		Ai=A(f)
	Bq=A(Cd-(f|{aM for(aM,aN,aN)in AM}));B4=0
	if A7:
		Br={at:AB[at]for at in f}
		for(Al,Bs,Cj)in Ap(AM):
			B5=u(Bs[Al]);B6=D(u(B)*Br[A]for(A,B)in Bs.items()if A!=Al)%e;Bt=-(Cj+B6)*X(B5,-1,e)%e;Br[Al]=Bt
			if Bt!=AB[Al]:B4+=1
	AE,Bu=64,255;v=D(E[A][4]for A in C(d)if E[A][5])
	def AF(c):
		B=A(N[c])
		if B<=2:return h(B-3)
		C=1e18
		for D in N[c]:
			if Ab(E[D][0][c])and E[D][4]<C:C=E[D][4]
		return h(B-2)*C
	w=[(AF(al),al)for al in C(A8)if 1<=A(N[al])<=AE];z.heapify(w);B7=0;AP=[];Ck=set(f)
	while w and v/B(1,g)<BM:
		Cl,F=z.heappop(w);B8=A(N[F])
		if B8==0 or B8>AE:continue
		Bv=AF(F)
		if Bv>Cl+1e-09:z.heappush(w,(Bv,F));continue
		if B8==1:
			AQ=AY(iter(N[F]));AR=E[AQ][0]
			if not Ab(AR[F])or W(A(N[B])<2 for B in AR if B!=F):continue
			AP.append((F,Y(AR),E[AQ][6]));E[AQ][5]=I;g-=1;v-=E[AQ][4]
			for a in AR:
				N[a].discard(AQ)
				if 1<=A(N[a])<=AE:z.heappush(w,(AF(a),a))
			AR.clear();continue
		AS=-1;Bw=1e18
		for c in N[F]:
			if Ab(E[c][0][F])and E[c][4]<Bw:Bw=E[c][4];AS=c
		if AS<0:continue
		U=E[AS];A3=U[0];Bx=BQ[A3[F]];By=[au for au in N[F]if au!=AS];AG=L;B9=-U[4];AT={}
		for T in By:
			O=E[T];x=AZ(O[0][F],Bx);AU=r(x)
			if O[1]+BR(x)*U[1]>Bu:AG=I;break
			B9+=U[2]if AU==1 else AU*U[3];s=O[0]
			for(V,BA)in A3.items():
				i=s.get(V);t=AZ(x,BA);n=(-t[0],-t[1])if i is S else(i[0]-t[0],i[1]-t[1])
				if CN(n)>Bu:AG=I;break
				B9+=r(n)-(0 if i is S else r(i))
				if V!=F:
					if i is S:AT[V]=AT.get(V,0)+1
					elif n==(0,0):AT[V]=AT.get(V,0)-1
			if not AG:break
		if AG and W(A(N[B])-1+AT.get(B,0)<=0 for B in A3 if B!=F):AG=I
		if not AG:continue
		if v+B9>M.W_ROW*B(1,g-1):break
		AP.append((F,Y(A3),U[6]))
		for T in By:
			O=E[T];s=O[0];x=AZ(s[F],Bx);AU=r(x);O[1]+=BR(x)*U[1];Bz=U[2]if AU==1 else AU*U[3];O[2]+=Bz;O[3]+=U[3];v+=Bz
			if A7:O[6]=(O[6]-u(x)*U[6])%e
			for(V,BA)in A3.items():
				i=s.get(V);t=AZ(x,BA)
				if i is S:n=-t[0],-t[1];s[V]=n;N[V].add(T);v+=r(n)
				else:
					n=i[0]-t[0],i[1]-t[1];v+=r(n)-r(i)
					if n==(0,0):
						del s[V];N[V].discard(T)
						if 1<=A(N[V])<=AE:z.heappush(w,(AF(V),V))
					else:s[V]=n
			O[4]=O[2]+D(r(A)for A in s.values())
			for V in s:
				if A(N[V])<=AE:z.heappush(w,(AF(V),V))
			B7+=1
		U[5]=I;g-=1;v-=U[4]
		for a in A3:
			N[a].discard(AS)
			if 1<=A(N[a])<=AE:z.heappush(w,(AF(a),a))
		A3.clear()
	A4=[am for am in C(d)if E[am][5]];Am={ao for an in A4 for ao in E[an][0]};B_=A(Ck-(Am|{aO for(aO,aP,aP)in AP}));BB=0
	if A7:
		C0={av:AB[av]for av in Am}
		for(An,C1,Cm)in Ap(AP):
			B5=u(C1[An]);B6=D(u(B)*C0[A]for(A,B)in C1.items()if A!=An)%e;C2=-(Cm+B6)*X(B5,-1,e)%e;C0[An]=C2
			if C2!=AB[An]:BB+=1
	AH=A(A4);BC=0
	if A7:
		for c in A4:
			if(D(u(B)*AB[A]for(A,B)in E[c][0].items())+E[c][6])%e!=0:BC+=1
	Cn=l([(A(N[aw]),aw)for aw in C(k)if N[aw]]);A5={};AV={};BD={};BE=0
	for(_,F)in Cn:
		C3=BD.get(F,G);Co=l([(E[aF][4],aF)for aF in N[F]if aF not in A5 and Ab(E[aF][0][F])])
		for(_,c)in Co:
			AW=K
			for y in E[c][0]:
				if y!=F and y in AV and AV[y]>AW:AW=AV[y]
			Z=AW+.5*(C3-AW)
			if AW<Z<C3:
				A5[c]=F;AV[F]=Z;BE+=1
				for y in E[c][0]:
					if y!=F and BD.get(y,G)>Z:BD[y]=Z
				break
	BF=AH-BE;BG,Cp=M.C_UNIT_5,M.C_GEN_1;BH=C4=C5=C6=C7=K;C8=BW/B(AA,1e-12)
	for c in A4:
		Cq=A5.get(c,-1);BI=C9=BJ=0;Ao=K
		for(F,CA)in E[c][0].items():
			BK=Ae[F];Cr=F==Cq
			for b in Aa(CA[0],CA[1]):
				Z=J(b)
				if Z==0:continue
				BH+=G
				if Z==1:C4+=G
				CB=BG*Z if Z<=4 else Cp;C5+=CB
				if not Cr:
					Ao+=CB
					if Z<=4 and BK<=28.:
						if BK<=14.:BI+=Z
						elif BK<=18.25:C9+=Z
						else:BJ+=Z
		BL=E[c][3];Ao+=BL*BG
		if C8>K and BL>K:BJ+=Bb(BL*C8)
		CC=BI//4;Cs=BI-4*CC;CD=C9+Cs;CE=CD//3;Ct=CD-3*CE;Cu=(BJ+Ct)//2;Cv=(3*CC+2*CE+Cu)*BG;C6+=Ao;C7+=Ao-Cv
	Cw=C4/B(BH,G);Cx=C5/B(BH,G);CF=C7/B(C6,1e-12);CG=v/B(1,AH);CH=Y(ok=BC==0 and Bq==0 and B_==0 and B4==0 and BB==0 and AH>0 and BF>0 and B7>0 and CG<=M.W_ROW and CF<M.TAB_R_CAP,passes=Bh,s0=Ce,log2U=log2U,seed=seed,W_target=BM,initial_rows=d,texc=Aj,ncols=A(Am),surplus=AH-A(Am),nrec=BE,N_sw=AH,N_l=BF,N_sw_U=H(AH/k,6),N_l_U=H(BF/k,6),avg_w=H(CG,2),merges=B7,w_mid=H(AA,4),lsmall=H(Av,3),exp_unit=H(Cw,5),k5_exp=H(Cx,4),tab_r=H(CF,5),kernel_errors=BC,post_peel_cols_lost=Bq,post_merge_cols_lost=B_,repurge_stack_len=A(AM),merge_stack_len=A(AP),repurge_backsub_errors=B4,merge_backsub_errors=BB)
	if BN:CI=l([(AV[A5[aG]],A5[aG],aG)for aG in A4 if aG in A5]);CJ=[ay for(ax,ay,ax)in CI];Cy=set(CJ);Cz=[az for az in A4 if az not in A5];f=l({aI for aH in A4 for aI in E[aH][0]});C_=[aA for aA in f if aA not in Cy];D0=[(aB,E[aC][0][aB],[(aQ,aR)for(aQ,aR)in E[aC][0].items()if aQ!=aB])for(aD,aB,aC)in CI];D1=[[(aJ,aK)for(aJ,aK)in E[aE][0].items()]for aE in Cz];CH['hx_export']=Y(rec_order=CJ,recorded_pivots=D0,live_col_set=C_,live_rows_full=D1)
	return CH
def verify_trim_pairs():
	L=M.load_witnesses()['trim_pairs'];j=L['texc_new'];N=L['all_ok'];C=L['rows'];O=[];P=[];U=[];V=[];W=[];F=[]
	for(q,r,k,l,s,X,m,Z,G,n,a,b,c,o)in C:
		d=G>=100;N=N and G<=X and n<=m-(k-j)and(not d or b<M.TAB_R_CAP and c<=M.W_ROW and a<=M.K5_LLA*Z)
		if d:O.append(a/Z);V.append(b);W.append(c)
		else:F.append(G)
		P.append(G/X);U.append(o/Q.log(l))
	p=B(A[3]/2.**A[0]for A in C);e=Q.ceil(B(U)*100)/100;f=Q.ceil(p*1e4)/1e4;S=Q.ceil(e*(Q.log(f)+M.log2_U_z6(M.OP['lf'],M.OP['lg'])*Q.log(2.)));I=Q.ceil(M.log2_U_z6(M.OP['lf'],M.OP['lg']));T=8*I;J=512;g=Q.ceil(926./Q.log2(5.));assert g<=J;h=Y(parent_size_reset=8,col_unions=(3+2*T+12)/M.OP['exc'],accumulate=2+T+10,component_sort=64*I,kill_sweep=2+T+10,fallback_sweep=2+8*J+64*I);E=D(h.values());K=16*J;i=lambda x:Q.ceil(x*100)/100;return Y(ok=N and S<=4096 and E<=2.**16 and 4096*E+K<=2.**28,pass_ir_terms={u:H(v,1)for(u,v)in h.items()},pass_ir_per_relation=H(E,1),drop_once_ir=K,find_steps_max=I,row_entries_max=J,row_entries_bound=g,budget_log2_at_pass_cap=i(Q.log2(4096*E+K)),budget_log2_at_predicted_passes=i(Q.log2(S*E+K)),tab_r_max=Q.ceil(B(V)*1e3)/1e3,avg_w_max=Q.ceil(B(W)*100)/100,ncols_lt_100=[R(F),B(F)],pairs_ncols_lt_100=A(F),pairs=A(C),ncols_ratio_range=[Q.floor(R(P)*1000)/1000,Q.ceil(B(P)*1000)/1000],pairs_ncols_ge_100=A(O),k5_ratio_max=Q.ceil(B(O)*1e4)/1e4,passes_max=B(A[13]for A in C),extra_passes_max=B(A[13]-A[4]for A in C),passes_per_ln_s0_max=e,s0_over_U_max=f,passes_at_prod=S)
def verify_live_trim_pair(seed=20261120):A=verify_live_filtsim_z6(log2U=11,seed=seed,W_target=60,check_kernel=L,texc_abs=32);B=verify_live_filtsim_z6(log2U=11,seed=seed,W_target=60,check_kernel=L,texc_abs=M.TEXC_LA);return Y(ok=A[r]and B[r]and A[N]==A[T]==32 and B[N]==B[T]==P(M.TEXC_LA)and B[U]<=A[U]and B[V]<=A[V],old_texc32=A,new_texc6=B)
def verify_live_filtsim_ladder(seed=20261120):B=verify_live_filtsim_z6(log2U=12,seed=seed,W_target=60,check_kernel=L);D=verify_live_filtsim_z6(log2U=13,seed=seed,W_target=60,check_kernel=I);F=D[Z]<B[Z]and D[b]<B[b];G=[12,13,20,21,22,23,55.29];E=[d(e)for e in G];H=j(E[A][1]>E[A+1][1]for A in C(A(E)-1))and E[-1][1]<.25;return Y(ok=B[r]and D[r]and F and H and B[N]==B[T]and D[N]==D[T],dim_monotonic_12_to_13=F,w_mid_monotone_to_zero_ok=H,w_mid_prod_per_side_zero=L,w_mid_ladder={str(f):Y(lsmall=i,w_mid=k)for(f,(i,k))in c(G,E)},u12=B,u13=D)
def verify_c6_pgf_and_hea_campaign():
	from fractions import Fraction as b;BH=M.load_data();Ae=BH[AP];Z=M.load_witnesses()['pgf_campaign'];F=Z['tower'];U,V,Bg=P(F['u_t']),P(F['p_t']),h(F['wa_lll']);E,AM,A6,u,v=P(F['lb']),P(F['nA']),P(F['nB']),P(F['b_d2']),P(F['T_surv']);AN,BI=h(F['p_ecm']),P(F['miss_fb_step']);A8,d,e,Af=1<<E,1<<AM,1<<A6,1<<v;Ag,Ah,BJ,BK=F['qbins_hi_f'],F['qbins_hi_g'],F['qbins_lo_f'],F['qbins_lo_g'];BL,BM,Ai,A9=h(F['alpha_f']),h(F['alpha_g']),M.ZK_NUM,M.W_COEFFS;Aj=[-4-aa if ab==0 else-aa for(ab,aa)in g(A9)];Bh=[[o(ak,al)for(ak,al)in ac]for ac in Ae[AQ]];AO=sqrt(2.);BN=M._pmul_h_q([b(am)for am in A9],[b(an)for an in Aj]);BO=M._ppow_h([0,1,0,0,0,0],V,V);Ak=M._ppow_h(BO,V,V);Al=M._ppow_h(Ak,V,V);BP=M._ppow_h(Al,V**3,V);Am=BN==[b(6)]+[b(0)]*5 and 1+U*6==V and BP==[0,1,0,0,0,0]and Ak!=[0,1,0,0,0,0]and Al!=[0,1,0,0,0,0];An=[b(ad,3)for ad in M.SIGMA_NUM];AR=[[b(1)]+[b(0)]*5,An]
	for _ in C(4):AR.append(M._pmul_h_q(AR[-1],An))
	Ao=O(O(P(3*AR[B][A])for B in C(6))for A in C(6));m=M._mh_k;Ap=Q.sumprod
	def AA(p):B=O(Ap(A,p)for A in Ao);E=O(Ap(A,B)for A in Ao);C=m(m(p,B,44),E,88);A=-C[2]//27;D=C[0]//27-2*A;return D*D+2*A*A
	def Aq(pa2,pb2):A=pb2;B=m(A9,pa2,44);return AA((U*B[0]-A[0],U*B[1]-A[1],U*B[2]-A[2],U*B[3]-A[3],U*B[4]-A[4],U*B[5]-A[5]))//729
	def Ar(pa2,pb2):A=pa2;B=m(Aj,pb2,44);return AA((A[0]+B[0],A[1]+B[1],A[2]+B[2],A[3]+B[3],A[4]+B[4],A[5]+B[5]))//729
	def AB(xs):
		B=xs
		while A(B)>1:B=[B[aF]*B[aF+1]if aF+1<A(B)else B[aF]for aF in C(0,A(B),2)]
		return B[0]if B else 1
	AC=M._primes_upto(A8);BQ=[ae for ae in AC if ae>2048];AS=[ag for(af,ag)in g(BQ)if af%BI==0];As=AB(AS);BR,At,Au=AB([au for au in AC if au<=256]),AB([aw for aw in AC if 256<aw<=2048]),AB([ax for ax in AC if 2048<ax<=A8]);BS=At*Au;AD=M._isp;Av=M._brent;w=Q.gcd;Aw=1<<2*E;Bi=e*e;BT=H(AN*1000)
	def Ax(n):
		for B in(BR,At,Au):
			while n>1:
				A=w(n,B)
				if A==1:break
				n//=A
			if n==1:return 1
		return n
	def Ay(n):
		D=0
		while L:
			A=w(n,BS)
			if A==1:break
			n//=A
			if A<65536 or AD(A):D+=1
			else:
				B=[A]
				while B:
					C=B.pop()
					if C<65536 or AD(C):D+=1
					else:E=Av(C);B.append(E);B.append(C//E)
		return D
	def BU(side,q,sp_idx,n_fac):
		B=side<<56^q<<24^sp_idx<<8^1371938854
		for D in C(n_fac):
			A=(B^D)+0x9e3779b97f4a7c15&0xffffffffffffffff;A=(A^A>>30)*0xbf58476d1ce4e5b9&0xffffffffffffffff;A=(A^A>>27)*0x94d049bb133111eb&0xffffffffffffffff
			if(A^A>>31)%1000>=BT:return I
		return L
	def Az(cof,b_a=d,b_b=e):
		E=b_b;B=cof
		if B==1:return L,[]
		if B<Aw:return B<=E,[B]if B<=E else[]
		if B>E*E or AD(B):return I,[]
		F,G=[],[B]
		while G:
			C=G.pop()
			if C==1:continue
			if C<Aw or AD(C):
				if C>E:return I,[]
				F.append(C)
				if A(F)>2:return I,[]
				continue
			H=Av(C)
			if H==C:return I,[]
			G.append(H);G.append(C//H)
		J=D(1 for A in F if A>b_a);return A(F)-J<=1 and J<=1,F
	BV=M.build_exact_int_geom(Ae)[0]
	def BW(q,s,rho):
		E=[_eval_zk_mod([1 if aJ==aG else 0 for aJ in C(6)],s,q)for aG in C(6)];F=X(E[0],-1,q);B=[[0]*12 for aH in C(12)];B[0][0]=q
		for A in C(1,6):B[A][A]=1;B[A][0]=-E[A]*F%q
		for A in C(6):B[6+A][6+A]=1;B[6+A][0]=rho*E[A]*F%q
		DESC.bareiss_lll_12d(B,BV,121,25,I);return B,[([D(aI[B]*Ai[B][aK]for B in C(6))for aK in C(6)],[D(aI[6+B]*Ai[B][aL]for B in C(6))for aL in C(6)])for aI in B]
	BX=[((ay,1),)for ay in C(12)]+[((az,1),(aA,aB))for az in C(12)for aA in C(az+1,12)for aB in(1,-1)]+[((ao,1),(ap,ar),(aq,at))for ao in C(12)for ap in C(ao+1,12)for aq in C(ap+1,12)for ar in(1,-1)for at in(1,-1)]
	def AT(cside,cq,av,bv):
		F=cside;B=cq;C=M._ppowmod([0,1],B,M.H_COEFFS,B)
		if A(C)<2:C+=[0]*(2-A(C))
		C[1]=(C[1]-1)%B
		for E in M._roots_mod(M._pgcd(C,M.H_COEFFS,B),B):
			I,G=_eval_zk_mod(av,E,B),_eval_zk_mod(bv,E,B)
			if G==0:continue
			D=I*X(G,-1,B)%B
			if D==0:continue
			H=_eval_poly_mod(A9,E,B)
			if F==0 and(U*H*D*D-1)%B==0 or F==1 and(D*D-4-H)%B==0:return F,B,E,D
	def A_(side,q,s,rho,b_a=d,b_b=e,first=I):
		Y=first;X=side;o,U=BW(q,s,rho);c=d=0;I=[aM[0]for aM in U]+[[-aQ for aQ in aN[0]]for aN in U];J=[aO[1]for aO in U]+[[-aR for aR in aP[1]]for aP in U];p,r=(Ar,Aq)if X==0 else(Aq,Ar)
		for(t,B)in g(BX):
			if A(B)==3:R,T,e=B[0][0],B[1][0]+(0 if B[1][1]==1 else 12),B[2][0]+(0 if B[2][1]==1 else 12);D,E,K=I[R],I[T],I[e];F,G,L=J[R],J[T],J[e];N=D[0]+E[0]+K[0],D[1]+E[1]+K[1],D[2]+E[2]+K[2],D[3]+E[3]+K[3],D[4]+E[4]+K[4],D[5]+E[5]+K[5];P=F[0]+G[0]+L[0],F[1]+G[1]+L[1],F[2]+G[2]+L[2],F[3]+G[3]+L[3],F[4]+G[4]+L[4],F[5]+G[5]+L[5]
			elif A(B)==2:R,T=B[0][0],B[1][0]+(0 if B[1][1]==1 else 12);D,E=I[R],I[T];F,G=J[R],J[T];N=D[0]+E[0],D[1]+E[1],D[2]+E[2],D[3]+E[3],D[4]+E[4],D[5]+E[5];P=F[0]+G[0],F[1]+G[1],F[2]+G[2],F[3]+G[3],F[4]+G[4],F[5]+G[5]
			else:N,P=O(I[B[0][0]]),O(J[B[0][0]])
			if not W(N)or not W(P):continue
			f=O(A//3 for A in m(N,N,32));h=O(A//3 for A in m(P,P,32));V=p(f,h)
			if V<=0:continue
			i=Ax(V);j=not Y and i<Af;k,Z=Az(i,b_a,b_b)
			if not j and not k:continue
			Q=r(f,h)
			if Q%q!=0:continue
			Q//=q
			if Q<=0:continue
			l=Ax(Q)
			if j and l<Af:d+=1
			if k:
				u,a=Az(l,b_a,b_b)
				if u and w(V,As)==1 and w(Q,As)==1 and BU(X,q,t,A(Z)+A(a)+Ay(V)+Ay(Q)):
					H=[0]*12
					for(v,x)in B:
						n=o[v]
						for b in C(12):H[b]+=n[b]if x==1 else-n[b]
					if w(AA(N)//729,AA(P)//729)==1 or M.is_ideal_primitive_oh(H[:6],H[6:]):
						if Y:return(H[:6],H[6:],a,Z)if X==0 else(H[:6],H[6:],Z,a)
						c+=1
		return S if Y else(c,d)
	AU=lambda ps,qp,spy,rhp,b_a=d,b_b=e:A_(ps,qp,spy,rhp,b_a,b_b,L);AE=lambda side,q,s,rho:A_(side,q,s,rho);B0={}
	def B1(u):
		A=B0.get(u)
		if A is S:A=M.rho(u);B0[u]=A
		return A
	AV=lambda u:G if u<=G else B1(ceil(u*1024.)/1024.);AW=lambda u:G if u<=G else B1(floor(u*1024.)/1024.);B2=P(H(48*(AM-E)/(A6-E)));f=(A6-E)/48.;B3=[Q.log1p(f/(E+ah*f))if ah<B2 else K for ah in C(48)];B4=[Q.log1p(f/(E+ai*f))if ai>=B2 else K for ai in C(48)];BY=[(0,[G]),(1,B3),(1,B4),(2,DESC._conv(B3,B4))]
	def B5(m,alp):m=ceil(m);return R(G,D(C*AV((m-A*E-B*f+alp)/E)for(A,D)in BY for(B,C)in g(D)if C and A*E+B*f+A*f<=m))
	B6=lambda m:G if v>=m else R(G,B(G,v/E)*AW((floor(m)-v)/E));AX=6.*log2(U)
	def AF(side,q):F=log2(q);A=(F+1e1)/12.;I=F if side==1 else K;J=F if side==0 else K;L=D(G/6.*B5(12.*A+Ag[B]-I,BL)*B5(12.*A+AX+Ah[B]-J,BM)for B in C(6));Q=D(G/6.*B6(12.*A+BJ[B]-I)*B6(12.*A+AX+BK[B]-J)for B in C(6));M=12.*A+Ag[-1]-I;N=12.*A+AX+Ah[-1]-J;S=B(K,G-(G-AN)*((M+N)/8.+4.));O,P=M/E,N/E;T,U=B(AV(O),1e-300),B(AV(P),1e-300);V=D((AW(B(K,O-log2(A)/E))/T+AW(B(K,P-log2(A)/E))/U)/A for A in AS);return H(699.947008*L*S*B(K,G-V),2),H(1024.*R(G,Q+L),2)
	N=L;n=[];a=L
	for((x,B7,AG,B8,AY),(AH,AI))in c(Z['live_d2_chains'][:2],Z['live_d2_expected'][:2]):
		T=AU(*x,d,e)
		if T is S or AG not in(T[2]if B7==0 else T[3])or not u<AG<=d<x[1]:N=a=I;continue
		i=AT(B7,AG,T[0],T[1])
		if i is S:N=a=I;continue
		B9,BA=AE(*i);BB,BC=AF(i[0],i[1])
		if B9!=AH[3]or BA!=AH[5]or J(BB-AH[2])>.02 or J(BC-AH[4])>.02:N=a=I
		n.append([i[0],i[1],1,BB,B9,BC,BA]);y=AU(*i,u,u)
		if y is S or AY not in(y[2]if B8==0 else y[3])or not A8<AY<=u<AG:N=a=I;continue
		p=AT(B8,AY,y[0],y[1])
		if p is S:N=a=I;continue
		BD,AO=AE(*p);BE,BF=AF(p[0],p[1])
		if BD!=AI[3]or AO!=AI[5]or J(BE-AI[2])>.02 or J(BF-AI[4])>.02:N=a=I
		n.append([p[0],p[1],2,BE,BD,BF,AO])
	for((x,BG,AZ),AJ)in c(Z['live_d1_parents'][:2],Z['live_d1_expected'][:2]):
		T=AU(*x,d,e)
		if T is S or AZ not in(T[2]if BG==0 else T[3])or not A8<AZ<=d<x[1]:N=I;continue
		q=AT(BG,AZ,T[0],T[1])
		if q is S:N=I;continue
		z,A0=AE(*q);A1,A2=AF(q[0],q[1])
		if z!=AJ[3]or A0!=AJ[5]or J(A1-AJ[2])>.02 or J(A2-AJ[4])>.02:N=I
		n.append([q[0],q[1],1,A1,z,A2,A0])
	for(Aa,Ab,BZ,Ba,Bb,Bc,Bd,Be)in Z['live_ctrl_ideals'][:2]:
		z,A0=AE(Aa,Ab,BZ,Ba);A1,A2=AF(Aa,Ab)
		if z!=Bc or A0!=Be or J(A1-Bb)>.02 or J(A2-Bd)>.02:N=I
		n.append([Aa,Ab,0,A1,z,A2,A0])
	A3=[[P(aj[0]),M.b62d(aj[1]),P(aj[2]),M.b62d(aj[3])/1e2,M.b62d(aj[4]),M.b62d(aj[5])/1e2,M.b62d(aj[6])]for aj in(A.split()for A in Z['compact_224_b62'].split(';'))]
	def A4(rs):I=rs;C=A(I);E=[aS[3]for aS in I];F=[aT[4]for aT in I];L=[aU[5]for aU in I];M=[aV[6]for aV in I];Q=D(F)/B(D(E),1e-12);S=D(F)/C;T=H(D((A-Q*B)**2 for(A,B)in c(F,E))/B(1,C-1)/B(S,1e-12),3);U=H(D((G+.5*A)**-G for A in E)/C,6);N=H(D((G+.5*A)**-G for A in E),3);J=l(A/B for(A,B)in c(F,E));O={str(aW):Y(emp=H(D(aW**B for B in F)/C,8),h_pgf_1_r1=H(D((G+B*(G-aW))**-G for B in E)/C,8),h_pgf_05_r1=H(D((G+.5*B*(G-aW))**-G for B in E)/C,8),ok=D(aW**B for B in F)/C<=D((G+B*(G-aW))**-G for B in E)/C and D(aW**B for B in F)/C<=D((G+.5*B*(G-aW))**-G for B in E)/C)for aW in(K,.2,.5,.8)};return Y(n=C,mu_min=H(R(E),2),mu_max=H(B(E),2),M_min=R(F),M_max=B(F),n_zero=D(1 for A in F if A==0),nb_r1_p0=U,exp_p0_sum=N,direct_oe_p0=H(D(1 for A in F if A==0)/B(N,1e-12),4),min_M_over_mu=H(J[0],4),p05_M_over_mu=H(J[B(0,P(.05*C))],4),below_half_mu=D(1 for A in J if A<=.5),below_08_mu=D(1 for A in J if A<=.8),res_fano_v_over_m=T,mean_M_over_mu=H(D(J)/C,4),two_sided_surv_obs_over_up=H(D(M)/D(L),4),two_sided_surv_max_over_up=H(B(A/B for(A,B)in c(M,L)),4),pgf_ok=j(A[r]for A in O.values()),pgf=O)
	k,AK,AL,A5,Ac=A4(A3),A4([aC for aC in A3 if aC[2]==1]),A4([aD for aD in A3 if aD[2]==2]),A4([aE for aE in A3 if aE[2]==0]),A4(n);Bf=Am and N and a and A(A3)==224 and A(n)==8 and k['n_zero']==0 and k['mu_min']>=25. and B(k[A7],AK[A7],AL[A7],A5[A7])<3. and k[s]and AK[s]and AL[s]and A5[s]and Ac[s]and k[t]>=G and k[Ad]<G and Ac[Ad]<G;return Y(ok=Bf,unified_bareiss_910_used=L,exact_ideal_primitivity_checked=L,tower_verified=Am,live_16_usable_and_depth2_exact=N,live_4_depth2_2hop_chains_exact=a,tower_params=Y(u_t=U,p_t=V,lb=E,nA=AM,nB=A6,b_d2=u,T_surv=v,ball_pts=1024,p_ecm=AN,miss_fb_primes=A(AS)),usable_thinning_totals=Z['usable_thinning_totals'],slot_d1_vs_ctrl_M_over_mu=H(AK[t]/A5[t],4),slot_d2_vs_ctrl_M_over_mu=H(AL[t]/A5[t],4),all_224=k,slot_d1_100=AK,slot_d2_24=AL,ctrl_100=A5,live_16=Ac)
def verify_live_c6_pool_hnf_subsample(n_tot=1000,seed=20261121):
	w=M.load_data();j=w[AP];E=w[Be];AO=j['GA'];k=[[o(af,ag)for(af,ag)in aa]for aa in j[AQ]];x,y,An,Ao,Ap,AR=M.build_exact_int_geom(j);d=Bc;l=E['wa'];z=E['wb'];m=E['ra0'];AS=E['rb0_ga'];n=l*P(H(m*m))*AR
	def AT(a,b):
		for A in C(6):
			if l*d(y[A],a)+z*d(y[A],b)<0:return I
		return L
	AU=M.is_ideal_primitive_oh
	def AV(z):return((((z-2)*z+1)*z+2)*z*z-4)*z+3
	V=[o(.4+0.9j)**ab for ab in C(6)]
	for _ in C(150):
		for b in C(6):
			A0=G+0j
			for O in C(6):
				if b!=O:A0*=V[b]-V[O]
			V[b]-=AV(V[b])/A0
	A1=[R(V,key=lambda r:J(r-k[ac][1]))for ac in C(3)];A2=lambda cs,z:D(B*z**A for(A,B)in g(cs));AW=[A2([0,-2,-1,0,1,-1],ad)for ad in A1];AX=[[A2(M.F_REC_COEFFS[ah],ae)for ah in(0,2,4,6,8)]for ae in A1];AY=h(M.U_BN254)
	def AZ(G):
		B=[[K]*6 for al in C(6)]
		for E in C(6):
			for A in C(E+1):F=G[E][A]-D(B[E][C]*B[A][C]for C in C(A));B[E][A]=sqrt(F)if E==A else F/B[A][A]
		return B
	def Aa(L):
		B=[[K]*6 for am in C(6)]
		for A in C(6):
			B[A][A]=G/L[A][A]
			for E in C(A):B[A][E]=-D(L[A][C]*B[C][E]for C in C(E,A))/L[A][A]
		return B
	A3=Aa(AZ(AO));A4=random.Random(seed);X=N=e=p=0;Q=[];S=[];Ab=(1,0,2,2,0,1),(1,2,0,2,1,2),(1,2,2,1,2,2);q=r=0
	def Ac(a,b):return W(D(A*B for(A,B)in c(A,a))%3==0 and D(A*B for(A,B)in c(A,b))%3==0 for A in Ab)
	while X<n_tot:
		s=[A4.gauss(0,1)for ai in C(12)];A5=A4.random()**(G/12.)/sqrt(D(A*A for A in s));T=[P(H(D(A3[B][aj]*s[B]*A5*m for B in C(aj,6))))for aj in C(6)];U=[P(H(D(A3[B][ak]*s[6+B]*A5*AS for B in C(ak,6))))for ak in C(6)]
		if not W(T)or not W(U):continue
		A6=l*d(x,T);A7=z*d(x,U)
		if A6+A7>n:continue
		X+=1
		if 20*A6>17*n or 20*A7>17*n:continue
		if not AT(T,U):continue
		N+=1;A8=Ac(T,U);q+=A8
		if AU(T,U):
			r+=A8;e+=1;A9=AA=K
			for O in C(3):AB=D(k[O][A]*T[A]for A in C(6));AC=D(k[O][A]*U[A]for A in C(6));Z,F=AB*AB,AC*AC;Ad,Ae,Af,Ag,Ah=AX[O];Ai=(((Ah*Z+Ag*F)*Z+Af*F*F)*Z+Ae*F*F*F)*Z+Ad*F*F*F*F;Aj=F+((AY+1)*AW[O]+2.)*Z;A9+=2.*log2(J(Ai));AA+=2.*log2(J(Aj))
			Q.append(A9+M.CONTENT_BITS);S.append(AA)
		else:p+=1
	t=E['bins'];AD=D(A[2]for A in t);AE=D(A[0]*A[2]for A in t)/AD;AF=D(A[1]*A[2]for A in t)/AD;AG=N/h(X);AH=e/B(1,N);f=E['n_prim_in']/h(E['n_in']);u=D(Q)/A(Q);v=D(S)/A(S);AI=(AG-E[a])/sqrt(B(1e-12,E[a]*(G-E[a])/X));AJ=(AH-f)/sqrt(f*(G-f)/B(1,N));Ak=sqrt(D((A-u)**2 for A in Q)/B(1,A(Q)-1)/A(Q));Al=sqrt(D((A-v)**2 for A in S)/B(1,A(S)-1)/A(S));AK=(u-AE)/Ak;AL=(v-AF)/Al;i=217./729.;AM=q/h(B(1,N));AN=(AM-i)/sqrt(i*(G-i)/B(1,N));Am=J(AI)<=2.5 and J(AJ)<=2.5 and J(AK)<=2.5 and J(AL)<=2.5 and e+p==N and r==0 and J(AN)<=2.5;return Y(ok=Am,mod3_bad_cells=q,mod3_bad_frac=H(AM,4),mod3_bad_frac_exact=H(i,4),z_mod3_bad=H(AN,2),mod3_bad_but_primitive=r,zero_float_fd_and_hnf=L,region='E12_JointVor',n_samples=X,n_active=N,n_prim=e,n_nonprim=p,frac_reg=H(AG,4),frac_reg_200k=E[a],z_frac_reg=H(AI,2),prim_rate=H(AH,4),prim_rate_200k=H(f,4),z_prim_rate=H(AJ,2),mean_lnf=H(u,2),mean_lnf_200k=H(AE,2),z_lnf=H(AK,2),mean_lng=H(v,2),mean_lng_200k=H(AF,2),z_lng=H(AL,2))
def verify_live_genuine_c6_relation_filter(seed=20261120):import c6_adaptive as A;return A.verify_live_genuine_c6_relation_filter(seed=seed)
verify_live_bareiss_12d_lll=DESC.verify_live_bareiss_12d_lll
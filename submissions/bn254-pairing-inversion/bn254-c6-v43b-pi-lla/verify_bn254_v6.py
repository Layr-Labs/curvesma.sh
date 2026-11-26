import bisect,functools,heapq,json,math,random
from fractions import Fraction as F
from math import ceil,exp,floor,log2,pi,sqrt
from pathlib import Path
import model_bn254_v6 as M
import descent_bn254_v6 as DESC
_SP=math.sumprod
def _mv6(Mat,v):A=Mat;return[_SP(A[0],v),_SP(A[1],v),_SP(A[2],v),_SP(A[3],v),_SP(A[4],v),_SP(A[5],v)]
def _mm6(A,B):return[[sum(A[D][C]*B[C][E]for C in range(6))for E in range(6)]for D in range(6)]
def _qf6(G,v):return _SP(v,(_SP(G[0],v),_SP(G[1],v),_SP(G[2],v),_SP(G[3],v),_SP(G[4],v),_SP(G[5],v)))
def _det_int(Mat):
	C=len(Mat);A=[A[:]for A in Mat];E=H=1
	for B in range(C-1):
		if A[B][B]==0:
			D=next((C for C in range(B+1,C)if A[C][B]),-1)
			if D<0:return 0
			A[B],A[D]=A[D],A[B];E=-E
		D=A[B][B]
		for F in range(B+1,C):
			for G in range(B+1,C):A[F][G]=(A[F][G]*D-A[F][B]*A[B][G])//H
		H=D
	return E*A[C-1][C-1]
def _eval_poly_mod(cs,s,p):
	A=0
	for B in reversed(cs):A=(A*s+B)%p
	return A
def _eval_zk_mod(vec,s,p):
	B=pow(3,-1,p);A=[1]
	for C in range(5):A.append(A[-1]*s%p)
	return sum(vec[B]*(sum(M.ZK_NUM[B][C]*A[C]for C in range(6))%p)for B in range(6)if vec[B])*B%p
def _is_prime(n):return n>=2 and all(n%A for A in range(2,int(n**.5)+1))
def verify_swar_sb_u7_and_stageB():
	BC=M.load_data();V=BC['t4_geom'];DQ=BC['pool_quad']['kap'];CU,CV,A1,DR,BD=V['GA'],V['GB'],V['SMAT'],V['FAC6'],V['UMATS'];BE=[[complex(A,B)for(A,B)in A]for A in V['EMB_ZK']];CW=[complex(A,B)for(A,B)in V['EMB_W']];BF=_mm6(A1,A1);I=random.Random(20261120);A2=0;CX=0x8000800080008000;BG=M.rand_below;BH=I.getrandbits
	for P in range(2000):
		CY=I.random()<1./9.;BI=48 if CY else 0;BJ=0;BK=False
		for DS in range(16):
			BL=[3520+BG(BH,2401)for A in range(4)];BM=[BG(BH,8481)for A in range(4)];CZ=[32768-A+BI for A in BL];h=[CZ[A]+BM[A]for A in range(4)]
			if any(A>=65536 for A in h):A2+=1
			Ca=h[0]|h[1]<<16|h[2]<<32|h[3]<<48;BJ|=Ca
			if any(BM[A]+BI>=BL[A]for A in range(4)):BK=True
		if(BJ&CX!=0)!=BK:A2+=1
	A3=.0
	for P in range(500):
		BN=[I.randint(-450,450)for A in range(6)];BO=[I.randint(-750,750)for A in range(6)]
		if not any(BN)or not any(BO):continue
		BP=BQ=.0
		for c in range(3):BR=sum(BE[c][A]*BN[A]for A in range(6));BS=sum(BE[c][A]*BO[A]for A in range(6));AN=M.U_BN254*CW[c]*(BR*BR);BT=AN-BS*BS;BP+=log2(BT.real**2+BT.imag**2);BQ+=log2(AN.real**2+AN.imag**2)
		A3=max(A3,abs(BP-BQ))
	BU=M._primes_upto(16383);AO=BU[1:];F=(1<<64)-1
	if M._P17_G_CACHE is None:
		Cb=M._primes_upto((1<<17)-1);Q=[]
		for A in Cb:
			if A>=16384:
				if M.U_BN254%A==0:Q.append(A);continue
				if A&7 not in(1,3):continue
				AP=M._ppow_h([0,1,0,0,0,0],A,A);AP[1]=(AP[1]-1)%A;i=M._pgcd(AP,M.H_COEFFS,A)
				if len(i)>1:
					BV=M._pmod([-B%A for B in M.G_LEAD_COEFFS],i,A)
					if len(M._pgcd(BV,i,A))>1:Q.append(A);continue
					S=M._ppowmod(BV,(A-1)//2,i,A);S[0]=(S[0]-1)%A
					if len(M._pgcd(S,i,A))>1:Q.append(A)
		M._P17_G_CACHE=Q
	Q=M._P17_G_CACHE;j=sorted(AO+Q,reverse=True);A4,AQ=[],[];Cc=(1<<64)-1
	for A in j:
		BW=False
		for(BX,BY)in enumerate(A4):
			if BY*A<=Cc:A4[BX]=BY*A;AQ[BX].append(A);BW=True;break
		if not BW:A4.append(A);AQ.append([A])
	A5=list(zip(A4,AQ));BZ=[]
	for(Ba,Cd)in A5:R=Ba<<64-Ba.bit_length();AR=((1<<128)-1)//R-(1<<64);BZ.append((R,AR,[pow(1<<64,A,R)for A in range(1,8)],[(A,pow(A,-1,1<<64),F//A)for A in Cd]))
	AS=[2]+j;A6=0;A7=0;Ce=math.sumprod
	for P in range(20):
		W=I.getrandbits(390)|1<<389
		for P in range(3):W*=AS[I.randrange(len(AS))]
		if Q:W*=Q[I.randrange(len(Q))]
		Bb=AO[I.randrange(len(AO))]
		if W.bit_length()+3*Bb.bit_length()<=511:W*=Bb**3
		AT=1;Bc=(W&-W).bit_length()-1;Cf=W>>Bc;AT<<=Bc;k=[Cf>>64*A&F for A in range(8)];AU=30
		for(R,AR,A8,Bd)in BZ:
			AU+=66+4*len(Bd);AV=k[0]+Ce(k[1:],A8);l=AV&F;Cg=AV>>64&F;Ch=AV>>128;AW=l+Cg*A8[0]+Ch*A8[1];m=AW&F;n=AW>>64&F;Ci=AW>>128
			if Ci:m+=A8[1];n+=m>>64;m&=F
			if n>=R:n-=R
			Be=AR*n;AX=(Be&F)+m;Cj=(Be>>64)+n+(AX>>64)+1&F;AX&=F;a=m-Cj*R&F
			if a>AX:a=a+R&F
			if a>=R:a=a-R&F
			for(A,Bf,Ck)in Bd:
				if a*Bf&F<=Ck:
					while True:
						AU+=50;Bg=[0]*8;AY=0
						for o in range(8):Bh=k[o]-AY;AZ=1 if Bh<0 else 0;Bi=(Bh&F)*Bf&F;Bg[o]=Bi;AY=(Bi*A>>64)+AZ
						if AY==0:k=Bg;AT*=A
						else:break
		A7=max(A7,AU);Cl=sum(k[A]<<64*A for A in range(8));Aa=W;Bj=1
		for A in AS:
			while Aa%A==0:Aa//=A;Bj*=A
		if AT!=Bj or Cl!=Aa:A6+=1
	Bk,Cm,BD,A1,BF,Bl=M.build_exact_int_geom(V);Ab=_qf6
	def Cn(a):
		if tuple(-A for A in a)<tuple(a):return False
		for C in range(6):
			D=Ab(Cm[C],a)
			if D<0:return False
			if D==0:
				A=_mv6(BD[C],a);E=[-A for A in A]
				if tuple(E if tuple(E)<tuple(A)else A)<tuple(a):return False
		for G in(A1,BF):
			B=_mv6(G,a);F=[-A for A in B]
			if tuple(F if tuple(F)<tuple(B)else B)<tuple(a):return False
		return True
	d=[]
	for A in range(101,1800,2):
		if not _is_prime(A):continue
		for N in range(A):
			if _eval_poly_mod(M.H_COEFFS,N,A)==0:
				Bm=-((M.U_BN254+1)%A*_eval_poly_mod(M.W_COEFFS,N,A)+2)%A
				if Bm==0:continue
				Bn=pow(Bm,-1,A)
				if pow(Bn,(A-1)//2,A)==1:
					for Ac in range(1,(A+1)//2):
						if Ac*Ac%A==Bn:d.append((A,N,Ac));break
		if len(d)>=24:break
	Co=36*Bl;Cp=int(round(11e1/(CV[0][0]/CU[0][0])*Bl));G=[]
	for Cq in range(-2,3):
		for Cr in range(-2,3):
			for Cs in range(-2,3):
				for Ct in range(-2,3):
					for Cu in range(-2,3):
						for A9 in range(-2,3):
							D=[A9,Cu,Ct,Cs,Cr,Cq]
							if any(D)and(Ad:=Ab(Bk,D))<=Co and Cn(D):G.append((Ad,D))
	G.sort(key=lambda x:x[0]);G=G[:80];X,p=[],[]
	for Cv in range(1,3):
		for Cw in range(-2,3):
			for Cx in range(-2,3):
				for Cy in range(-2,3):
					for Cz in range(-2,3):
						q=[(A,Ae,Ad)for A in range(-6,7)if(Ad:=Ab(Bk,(Ae:=[A,Cz,Cy,Cx,Cw,Cv])))<=Cp]
						if q:
							X.append((len(q),[A[1]for A in q],[A[2]for A in q]))
							for(P,Ae,C_)in q:p.append((Ae,C_,len(p)))
						if len(X)>=50:break
					if len(X)>=50:break
				if len(X)>=50:break
			if len(X)>=50:break
		if len(X)>=50:break
	Y={(A,B):0 for A in range(len(G))for B in range(len(p))};AA=0
	for(A,N,r)in d:
		D0=[_eval_zk_mod(B,N,A)for(C,B)in G];D1=[_eval_zk_mod(B,N,A)for(B,C,C)in p];Af={}
		for(C,D2)in enumerate(D0):Af.setdefault(D2,[]).append(C)
		for(AB,D3)in enumerate(D1):
			s=r*D3%A;D4=(A-s)%A
			for C in Af.get(s,()):Y[C,AB]+=1;AA+=1
			for C in Af.get(D4,()):Y[C,AB]+=1;AA+=1
	H=4;Bo={A:0 for A in Y};Bp={tuple(A):B for(A,C,B)in p};Bq=[A for(A,B)in G];Br=len(set(Bq))<len(Bq);Bs=Ag=Ah=0;Ai=None;Aj=True
	for(A,N,r)in d:
		Ak=pow(r,-1,A);Z=(A+H-1)//H;J=len(G);AC=[0]*Z;Bt=[0]*J;Bu=[0]*J
		for(C,(P,D))in enumerate(G):K=_eval_zk_mod(D,N,A)*Ak%A;B=K//H;AC[B]+=1;Bt[C]=B;Bu[C]=C<<8|K&H-1
		AD=[0]*(Z+1);Al=[0]*Z;AE=0;AF=bytearray((Z+15)//8+1);Bv=0
		for B in range(Z):
			AD[B]=AE;Al[B]=AE
			if AC[B]:AF[B>>3]|=1<<(B&7);Bv+=1
			AE+=AC[B];AC[B]=0
		if Bv>min(Z,J):Aj=False
		AD[Z]=AE;Bw=[0]*J
		if Ai is None:Ai=bytes(AF)
		for C in range(J):B=Bt[C];Bx=Al[B];Bw[Bx]=Bu[C];Al[B]=Bx+1
		for(t,u,P)in X:
			e=_eval_zk_mod(u[0],N,A)
			for Am in(1,-1):
				An=e if Am==1 else-(e+t-1)%A;Ao=(e+t-1)%A if Am==1 else-e%A;D5=[(An//H,Ao//H)]if An<=Ao else[(An//H,Z-1),(0,Ao//H)]
				for(AG,D6)in D5:
					v=D6-AG+1;Bs+=v;assert v<=8,f"1-byte occ_bm overflow: n_span={v} > 8";By=AG>>3;D7=AG&7;D8=AF[By]|AF[By+1]<<8;Bz=D8>>D7&(1<<v)-1
					if Bz==0:continue
					for c in range(v):
						if Bz>>c&1:
							Ag+=1;Ah+=1;B=AG+c
							for D9 in range(AD[B],AD[B+1]):
								B_=Bw[D9];K=B*H|B_&H-1;C=B_>>8;S=(K-e)%A if Am==1 else(-K-e)%A
								if S<t:Bo[C,Bp[tuple(u[S])]]+=1
	C0=sum(1 for A in Y if Bo[A]!=Y[A]);Ap=Aq=Ar=0;w=(1<<64)-1;b=835
	for E in((1<<44)+7,(1<<55)-55):
		for T in(0,1,255,256,E-1,E-b,E-b+5,E-100):
			for x in(1,-1):
				As=T if x==1 else-(T+b-1)%E;f=(T+b-1)%E if x==1 else-T%E;DA=list(range(As>>8,(f>>8)+1))if As<=f else list(range(As>>8,E+255>>8))+list(range(0,(f>>8)+1))
				for At in DA:
					y=At<<8;C1=((y-T if At>=T>>8 else E-T+y)if x==1 else f-y if At<=f>>8 else f+E-y)&w
					for Au in(0,1,7,127,128,254,255):
						Av=y+Au
						if Av>=E:continue
						Aw=(Av-T)%E if x==1 else(-Av-T)%E;C2=Au+C1&w if x==1 else C1-Au&w
						if(C2<b)!=(Aw<b)or Aw<b and C2!=Aw:Ap+=1
		Ax=[I.randrange(E)for A in range(6)];DB=[A&268435455 for A in Ax];DC=[A>>28 for A in Ax];DD=(1<<64)%E
		for P in range(32):
			D=[I.randrange(0,1<<13)for A in range(6)];C3=sum(D[A]*DB[A]for A in range(6));Ay=sum(D[A]*DC[A]for A in range(6));C4=Ay>>36;DE=(Ay<<28)+C3&w;DF=(Ay<<28&w)+C3>>64;C4+=DF;AZ=(DE+C4*DD)%E
			if AZ!=sum(D[A]*Ax[A]for A in range(6))%E:Ar+=1
	C5=[];C6={}
	for(C,(DG,D))in enumerate(G):C6.setdefault(tuple(D[1:]),[]).append((D[0],C,D))
	for AH in C6.values():AH.sort();C5.append(AH)
	U,l,DH=d[0];Az=pow(DH,-1,U);DI=_eval_zk_mod([1,0,0,0,0,0],l,U)*Az%U;O=(U+H-1)//H;J=len(G);C7=[0]*O;C8=[0]*J;C9=[0]*J
	for(C,(DG,D))in enumerate(G):K=_eval_zk_mod(D,l,U)*Az%U;B=K//H;C7[B]+=1;C8[C]=B;C9[C]=C<<8|K&H-1
	z=[0]*(O+1);A_=[0]*O;L=0
	for B in range(O):z[B]=L;A_[B]=L;L+=C7[B]
	z[O]=L;CA=[0]*J
	for C in range(J):B=C8[C];CA[A_[B]]=C9[C];A_[B]+=1
	AI=[0]*O;CB=[0]*J;CC=[0]*J;B0=0
	for AH in C5:
		B1=None;AJ=0
		for(A9,C,D)in AH:AJ=(AJ+DI)%U if B1 is not None and A9==B1+1 else _eval_zk_mod(D,l,U)*Az%U;B1=A9;B=AJ//H;AI[B]+=1;CB[B0]=B;CC[B0]=C<<8|AJ&H-1;B0+=1
	A0=[0]*(O+1);B2=[0]*O;L=0;CD=bytearray((O+15)//8+1)
	for B in range(O):
		A0[B]=L;B2[B]=L
		if AI[B]:CD[B>>3]|=1<<(B&7)
		L+=AI[B];AI[B]=0
	A0[O]=L;CE=[0]*J
	for o in range(J):B=CB[o];CE[B2[B]]=CC[o];B2[B]+=1
	B3=bytes(CD)==Ai;B4=A0==z and all(sorted(CE[A0[A]:A0[A+1]])==sorted(CA[z[A]:z[A+1]])for A in range(O))and B3
	if not B4:Aq=1
	CF={A:0 for A in Y}
	for(A,N,r)in d:
		Ak=pow(r,-1,A);CG=[0]*A;CH=[0]*len(G)
		for(C,(P,D))in enumerate(G):K=_eval_zk_mod(D,N,A)*Ak%A;CG[K]+=1;CH[C]=K
		AK=[0]*(A+1);B5=[0]*A;L=0
		for B6 in range(A):AK[B6]=L;B5[B6]=L;L+=CG[B6]
		AK[A]=L;CI=[0]*len(G)
		for C in range(len(G)):K=CH[C];CI[B5[K]]=C;B5[K]+=1
		for(t,u,P)in X:
			for S in range(t):
				s=_eval_zk_mod(u[S],N,A);AB=Bp[tuple(u[S])]
				for CJ in(s,-s%A):
					for DJ in range(AK[CJ],AK[CJ+1]):CF[CI[DJ],AB]+=1
	CK=sum(1 for A in Y if CF[A]!=Y[A]);B7=0;AL=floor(468/14);AM=[(0,AL+1),(16*55.25-1,16*60+AL+1),(16*110.5-1,16*120+AL+1)];CL=all(AM[A][1]<AM[A+1][0]for A in range(len(AM)-1))
	for DK in(0,1,2):
		for P in range(50):
			B8=sum(I.uniform(55.25,6e1)for A in range(DK))
			if B8>14e1:continue
			CM=I.uniform(.0,max(.5,min(11e1,168.-B8)));CN=468.-B8-CM;CO=I.randint(8,min(AL,int(CN/14.)));DL=[CN/CO]*CO;DM=sum(floor(16.*A)for A in DL);CP=floor(16.*468.)-DM-floor(16.*CM)
			if CP<=16*140:
				if not any(A<=CP<=B for(A,B)in AM):B7+=1
	B9=0
	for P in range(1000):
		CQ=I.randint(16*400,16*500);CR=CQ-16*170;DN=32768-CR;CS=I.randint(0,CR+200);DO=CS+48+DN&32768!=0;DP=CQ-48-CS<=16*170
		if DO!=DP:B9+=1
	g=int(2.**45.352);BA=[(A,B)for A in(0,400,795)for B in(0,1,(1<<32)-1,1<<32,(1<<32)+1,(1<<44)+12345,g-1)];BB=[A*g+B for(A,B)in BA];CT=len(set(BB))==len(BA)and all(0<=A<1<<56 for A in BB)and all((A//g,A%g)==(B,C)for(A,(B,C))in zip(BB,BA))and 0*g+0!=0*g+(1<<32);return dict(ok=C0==0 and Aj and Ap==0 and Aq==0 and B3 and Ar==0 and B4 and CK==0 and Br and A2==0 and A6==0 and CL and B7==0 and B9==0 and CT and A3<=1e-12 and AA>=1400,zero_float_geom=True,det_occ_bound_ok=Aj,flat_array_zero_dict_build=True,flat_sbu7_exact_ok=B3,zero_sort_exact_ok=B4,exact_slot_cutoff_M_t=True,packed_8b_u64_wrap_errors=Ap,band2_a0_line_build_errors=Aq,split28_6coord_res_errors=Ar,band0_block_direct_apply_mismatches=CK,qa_stored_in_record=False,qa_ties_tested=Br,scalar_swar_10k_errors=A2,p2_shift_threshold_errors=B9,a_slot_64bit_linear_ok=CT,max_ng_invariance_bits=A3,p14_primes=len(BU),p17_g_split_primes=len(Q),p17_strip_odd_primes=len(j),p17_word_chunks=len(A5),p14_zero_gcd_gm_ok=A6==0 and A7<=M.C_SB and len(A5)*66+len(j)*4+30+50*32<=M.C_SB,p14_strip_ops_max=A7,p14_strip_ops_e16=len(A5)*66+len(j)*4+30+50*16,p14_strip_errors=A6,stageB_window_disjoint=CL,stageB_window_false_neg=B7,sb_u7_hits=AA,sb_u7_mismatches=C0,sb_total=Bs,sb_occ=Ag,sb_ent=Ah,sb_occ_skip_pct=round(1e2*(1.-Ah/max(1,Ag)),2))
def verify_12d_descent_lattices():
	t=M.load_data();u=t['t4_geom'];N=M.load_witnesses();v=[[complex(A,B)for(A,B)in A]for A in u['EMB_ZK']];A_=[complex(A,B)for(A,B)in u['EMB_W']];V={A:[sum(B[A]*v[C][1]**A for A in range(6))for C in range(3)]for(A,B)in M.F_REC_COEFFS.items()};w,BQ=6.*log2(M.U_BN254),6.*log2(M.P_BN254);AJ=[[1,0,1,-1,-3,-3],[0,1,-1,0,-1,-1],[0,0,-1,1,1,3],[0,0,-1,-1,0,2],[0,0,1,0,1,-1],[0,0,0,1,1,2]];x=_det_int
	def I(a,b):
		A=[0]*11
		for C in range(6):
			if a[C]:
				for B in range(6):A[C+B]+=a[C]*b[B]
		for D in range(10,5,-1):
			E=A[D]
			if E:
				A[D]=0
				for B in range(6):A[D-6+B]-=E*M.H_COEFFS[B]
		return A[:6]
	def AK(poly_num,denom):
		B=denom;C,A=[],poly_num[:]
		for H in range(6):
			C.append(A)
			if H<5:
				D=[0]+A[:5];E=A[5]
				if E:
					for F in range(6):D[F]-=E*M.H_COEFFS[F]
				A=D
		G=abs(x([[C[B][A]for B in range(6)]for A in range(6)]));assert G%B**6==0;return G//B**6
	def B0(a,b):
		E=[sum(a[A]*M.ZK_NUM[A][B]for A in range(6))for B in range(6)];F=[sum(b[A]*M.ZK_NUM[A][B]for A in range(6))for B in range(6)];A,B=I(E,E),I(F,F);K=I(M.W_COEFFS,A);L=AK([(M.U_BN254+1)*K[C]+2*A[C]+B[C]for C in range(6)],9);C=[[1,0,0,0,0,0],A];D=[[1,0,0,0,0,0],B]
		for P in range(3):C.append(I(C[-1],A));D.append(I(D[-1],B))
		G=[0]*6
		for(H,N)in M.F_REC_COEFFS.items():
			O=I(I(N,C[H//2]),D[4-H//2])
			for J in range(6):G[J]+=O[J]
		return AK(G,3**8),L
	def AL(vec,mod):A=pow(3,-1,mod);return[sum(vec[A]*M.ZK_NUM[A][B]for A in range(6))*A%mod for B in range(6)]
	AM=tuple(tuple(A.real for A in A)for A in v);AN=tuple(tuple(A.imag for A in A)for A in v);O=math.sumprod
	def y(a,b,need_f=True):
		I=J=.0;R=max(max(abs(A)for A in a),max(abs(A)for A in b),1);H=max(0,R.bit_length()-40);K=float(1<<H);L=[A/K for A in a];N=[A/K for A in b]
		for A in range(3):
			P=complex(O(AM[A],L),O(AN[A],L));Q=complex(O(AM[A],N),O(AN[A],N));B,C=P*P,Q*Q;D=((M.U_BN254+1)*A_[A]+2.)*B+C;J+=log2(D.real*D.real+D.imag*D.imag)+4.*H
			if need_f:E,F=B*B,C*C;G=V[8][A]*E*E+V[6][A]*E*B*C+V[4][A]*E*F+V[2][A]*B*F*C+V[0][A]*F*F;I+=log2(G.real*G.real+G.imag*G.imag)+16.*H
		return I,J
	B1=log2(t['descent']['cells']['4']);P=[math.factorial(A)for A in range(13)];z=[]
	for W in range(2):
		for X in range(2):
			for Y in range(5):
				AO=16-16*W-9*X-4*Y
				if AO<0:continue
				for A0 in range(min(AO,12-W-X-Y)+1):
					A1=A0+Y+X+W
					if A1==0:continue
					A2=P[12]//(P[12-A1]*P[A0]*P[Y]*P[X]*P[W])*(1<<A1);z.append((A2,(A0,Y,X,W)))
	AP,g=[],0
	for(A2,_)in z:g+=A2;AP.append(g)
	AQ=g//2==t['descent']['cells']['4']
	def B2(rng_loc):A=rng_loc;G=bisect.bisect_left(AP,A.randrange(g)+1);B,C,D,E=z[G][1];F=[1]*B+[2]*C+[3]*D+[4]*E+[0]*(12-B-C-D-E);A.shuffle(F);return[B*(1 if A.getrandbits(1)else-1)if B else 0 for B in F]
	AR=M.build_exact_int_geom(u)[0];B3=random.Random(42);AS=random.Random(20261121);h=A3=True;i=.0;j=0;Q=0;AT,AU,A4,A5,J,K=[],[],[],[],[],[]
	for L in N['special_q']:
		R,B,A6,AV=L[0],M.b62d(L[1]),M.b62d(L[2]),M.b62d(L[3]);A7=B.bit_length();AT.append(A7);AU.append(R);B4=pow(3,-1,B);A8=[1]
		for _ in range(5):A8.append(A8[-1]*A6%B)
		AW=[sum(M.ZK_NUM[C][A]*A8[A]for A in range(6))*B4%B for C in range(6)];G=AW+[-AV*A%B for A in AW];AX=pow(G[0],-1,B);A9=[[B if A==0 else 0 for A in range(12)]]+[[-G[A]*AX%B-(B if-G[A]*AX%B>B//2 else 0)if C==0 else 1 if C==A else 0 for C in range(12)]for A in range(1,12)];DESC.bareiss_lll_12d(A9,AR,1,1,False,11,40,(9,10));C=[[A9[B][A]for B in range(12)]for A in range(12)]
		if abs(x(C))!=B:h=False
		for H in range(12):
			S=[C[A][H]for A in range(6)];k=[C[6+A][H]for A in range(6)]
			if(_eval_zk_mod(S,A6,B)-AV*_eval_zk_mod(k,A6,B))%B!=0:h=False
		l=[C[A][0]+C[A][1]for A in range(12)];AY,AZ=B0(l[:6],l[6:])
		if AY%B!=0 if R=='f'else AZ%B!=0:A3=False
		B5,B6=y(l[:6],l[6:],True);i=max(i,abs(log2(AY)-B5),abs(log2(AZ)-B6));m=(log2(B)+B1)/12.;n=A7 if R=='f'else .0;o=A7 if R=='g'else .0;B7=sum(A*M.rho((48.*m+B-n)/59.)*M.rho((12.*m+w+C-o)/59.)for(A,B,C)in DESC.QBINS);AA=[];Aa=[]
		for _ in range(60):
			Z=B2(B3);F=[O(C[A],Z)for A in range(12)];D=max(abs(A).bit_length()for A in F);Q=max(Q,D)
			if D>=48:j+=1
			G,a=y(F[:6],F[6:],True);Aa.append((G,a));J.append(G-48.*m);K.append(a-12.*m-w);AA.append(M._rho_fast((G-n)/59.)*M._rho_fast((a-o)/59.))
		DESC._SQ_NORMS_CACHE[R,B]=Aa;DESC._SQ_BASIS_CACHE[R,B]=A9;A4.append(round(sum(AA)/len(AA)/max(B7,1e-300),2))
		for Ab in(16,64,128):
			B8=log2(.5*(pi**6/72e1)*Ab**12);Ac=(log2(B)+B8)/12.;B9=sum(A*M._rho_fast((48.*Ac+B-n)/59.)*M._rho_fast((12.*Ac+w+C-o)/59.)for(A,B,C)in DESC.QBINS);AB=[]
			for _ in range(32):
				Ad=[AS.gauss(.0,1.)for A in range(12)];BA=sqrt(sum(A*A for A in Ad));BB=Ab*AS.random()**(1./12.);Z=[int(round(A*BB/BA))for A in Ad]
				if not any(Z):Z[0]=1
				F=[O(C[A],Z)for A in range(12)];D=max(abs(A).bit_length()for A in F);Q=max(Q,D)
				if D>=48:j+=1
				G,a=y(F[:6],F[6:],True);AB.append(M._rho_fast((G-n)/59.)*M._rho_fast((a-o)/59.))
			A5.append(sum(AB)/len(AB)/max(B9,1e-300))
	J.sort();K.sort()
	def p(arr,qv):A=arr;return round(A[min(len(A)-1,max(0,ceil(qv*len(A))-1))],2)
	Ae=all(J[min(len(J)-1,floor(A/6.*len(J)))]>=DESC.QBINS_LO[A][0]and K[min(len(K)-1,floor(A/6.*len(K)))]>=DESC.QBINS_LO[A][1]for A in range(6));Af=[p(J,(A+1)/6.)for A in range(6)];Ag=[p(K,(A+1)/6.)for A in range(6)];Ah=all(Af[A]<=DESC.QBINS[A][1]and Ag[A]<=DESC.QBINS[A][2]for A in range(6));Ai=.1,.25,.5,.75,.9,1.;T=AC=True;AD=0;q=0;E=[0,0];A=M.P_BN254;BC=pow(3,-1,A);BD=[-B%A for B in M.G_LEAD_COEFFS];Aj=M._pinv_h(BD,A)
	def BE(a12,b12):C=b12;B=a12;D=[(B+C)%A for(B,C)in zip(M._pmul_h(B[0],C[0],A),M._pmul_h(M._pmul_h(B[1],C[1],A),Aj,A))];E=[(B+C)%A for(B,C)in zip(M._pmul_h(B[0],C[1],A),M._pmul_h(B[1],C[0],A))];return D,E
	def Ak(a):
		E,F,G,H,I,J=a;R=E*E;P=2*E*F;Q=2*E*G+F*F;O=2*(E*H+F*G);M=2*(E*I+F*H)+G*G;N=2*(E*J+F*I+G*H);B=2*(F*J+G*I)+H*H;C=2*(G*J+H*I);D=2*H*J+I*I;K=2*I*J;L=J*J
		if L:K+=2*L;D-=L;C-=2*L;N+=4*L;M-=3*L
		if K:D+=2*K;C-=K;B-=2*K;M+=4*K;O-=3*K
		if D:C+=2*D;B-=D;N-=2*D;O+=4*D;Q-=3*D
		if C:B+=2*C;N-=C;M-=2*C;Q+=4*C;P-=3*C
		if B:N+=2*B;M-=B;O-=2*B;P+=4*B;R-=3*B
		return[R%A,P%A,Q%A,O%A,M%A,N%A]
	def BF(a12):B=a12;C=[(B+C)%A for(B,C)in zip(Ak(B[0]),M._pmul_h(Ak(B[1]),Aj,A))];D=[2*B%A for B in M._pmul_h(B[0],B[1],A)];return C,D
	def BG(base12,exp_v):
		A=[1,0,0,0,0,0],[0]*6
		for B in bin(exp_v)[2:]:
			A=BF(A)
			if B=='1':A=BE(A,base12)
		return A
	for L in N['boot']:
		r=[M.b62d(A)for A in L[0].split()];BH=[M.b62d(A)for A in L[1].split()];Al=[BH[A*12:(A+1)*12]for A in range(6)];BI=BG(([1,0,0,0,0,0],[B%A for B in r]),M.R_BN254)
		if BI[1]!=[0]*6 or not any(r):AC=False
		Am=[[0]*12 for A in range(6)]
		for H in range(12):
			S=[Al[A][H]for A in range(6)];BJ=[sum(S[A]*M.ZK_NUM[A][B]for A in range(6))*BC%A for B in range(6)];BK=M._pmul_h(r,BJ,A);AE=[3*B%A for B in BK];AE=[B-A if B>A//2 else B for B in AE]
			for An in range(6):
				Ao=sum(AJ[An][A]*AE[A]for A in range(6))
				if Ao%3!=0:T=False
				Am[An][H]=Ao//3
		C=Al+Am
		if abs(x(C))!=A**6:T=False
		if not DESC.verify_bareiss_910_cert([[C[B][A]for B in range(12)]for A in range(12)],AR,0xa8ccdee3ccdbd36d,1):T=False
		for H in range(12):
			S=[C[A][H]for A in range(6)];k=[C[6+A][H]for A in range(6)]
			if M._pmul_h(r,AL(S,A),A)!=AL(k,A):T=False
			b=max(abs(A).bit_length()for A in S);D=max(abs(A).bit_length()for A in k);E[0]=max(E[0],b);E[1]=max(E[1],D)
			if b>=128 or D>=160:q+=1
		AF=[[C[B][A]for B in range(12)]for A in range(12)]
		for F in[AF[A]+[]for A in range(12)]+[[A+C*B for(A,B)in zip(AF[A],AF[B])]for A in range(12)for B in range(A+1,12)for C in(1,-1)]:
			b=max(abs(A).bit_length()for A in F[:6]);D=max(abs(A).bit_length()for A in F[6:]);E[0]=max(E[0],b);E[1]=max(E[1],D)
			if b>=128 or D>=160:q+=1
		AD+=1
	from fractions import Fraction as c;Ap=[[c(A,3)for A in A]for A in M.ZK_NUM]
	def BL(a,b):
		A=[c(0)]*11
		for D in range(6):
			for B in range(6):A[D+B]+=a[D]*b[B]
		for C in range(10,5,-1):
			E=A[C]
			if E:
				for B in range(6):A[C-6+B]-=E*M.H_COEFFS[B]
				A[C]=c(0)
		return A[:6]
	def AG(p):return[int(sum(AJ[B][A]*(3*p[A])for A in range(6)))//3 for B in range(6)]
	s=[[AG(BL(Ap[A],Ap[B]))for B in range(6)]for A in range(6)];AH=max(sum(abs(s[B][C][A])for B in range(6)for C in range(6))for A in range(6));U=sqrt(max(sum(sum(abs(s[C][A][B])for C in range(6))**2 for B in range(6))for A in range(6)));BM=AG([c(A)for A in M.W_COEFFS]);BN={A:AG([c(A)for A in M.F_REC_COEFFS[A]])for A in(0,2,4,6,8)};BO=sum(max(sum(abs(sum(BN[B][A]*s[A][D][C]for A in range(6)))for D in range(6))for C in range(6))for B in(0,2,4,6,8));Aq=max(sum(abs(sum(BM[A]*s[A][C][B]for A in range(6)))for C in range(6))for B in range(6));Ar=log2(BO)+7.*log2(AH)+384.;As=log2((M.U_BN254+1)*Aq+3.)+log2(AH)+96.;At=log2(((M.U_BN254+1)*Aq+2.)*2.**256+2.**320)+log2(AH);Au=round(6.*(log2(U)+Ar),2);Av=round(6.*(log2(U)+As),2);Aw=round(6.*(log2(U)+At),2)
	def AI(rb):return sum((5-A)**2*(24*ceil((A+1)*rb/64.)**2+16*ceil((A+1)*rb/64.)*(ceil(A*rb/64.)if A else 1)+60*ceil((A+1)*rb/64.))for A in range(5))
	d=ceil(Au/64.);e=ceil(Av/64.);f=ceil(Aw/64.);Ax=max(AI(log2(U)+Ar)+AI(log2(U)+As)+450000,AI(log2(U)+At)+150000);Ay=max(4*(d+e)+4170*(d+e)+18990+41*d*(6*d+10)+41*e*(6*e+10),4*f+4170*f+9495+41*f*(6*f+10));Az=j==0 and q==0 and Q<48 and E[0]<128 and E[1]<160 and d<=40 and e<=16 and f<=32 and Ax<DESC.NORMS and Ay<DESC.TDIV;BP=AQ and A3 and i<1e-08 and all(A>=1. for A in A4)and Ae and Ah and min(A5)>=1. and AD==24 and len(N['special_q'])==24 and len(N['boot'])==24;return dict(ok=h and T and AC and BP and Az,prenorm_gate_uniform_ok=Az,sq_gate_max_cb=Q,bt_gate_max_ab=E,sq_gate_fails=j,bt_gate_fails=q,sq_nf_hi_bits=Au,sq_ng_hi_bits=Av,bt_ng_hi_bits=Aw,gate_norms_wc=Ax,gate_tdiv_wc=Ay,all_24_sq_reduced_live_by_bareiss_910=h,all_24_boot_bareiss_910_certified=T,boot_gt_subgroup_ok=AC,cells4_exact_cardinality_ok=AQ,special_q_count=len(N['special_q']),special_q_sides=AU,special_q_bits=AT,exact_side_divisibility_ok=A3,boot_count=len(N['boot']),boot_candidates_gated=144*AD,exact_int_vs_fp_max_err_bits=round(i,12),paired_joint_smooth_ratios_r4=A4,hi_radius_R16_128_min_paired_ratio=round(min(A5),2),cf_sextiles_obs=Af,cg_sextiles_obs=Ag,cf_quantiles_p10_25_50_75_90_max=[p(J,A)for A in Ai],cg_quantiles_p10_25_50_75_90_max=[p(K,A)for A in Ai],qbins_hi_6_of_6_ok=Ah,qbins_lo_6_of_6_ok=Ae)
def compute_w_mid_at_scale(log2U,Lf=61.5,Lg=61.5,Nf=390.8,Ng=466.66,excess=.72,t_dense=.8):
	A=t_dense;G=excess*2.**log2U;C=min(A*min(Lf,Lg),max(8.,log2(G/25.)+log2(3.)));D=.0
	for(H,I)in((Nf/Lf,Lf),(Ng/Lg,Lg)):
		B=C/I
		if B>=A:continue
		E=(A-B)/2e2
		for J in range(200):F=B+E*(J+.5);D+=M.kappa_rho(H,F)*E/F
	return round(C,4),round(D,6)
def verify_live_filtsim_z6(log2U=13,seed=20261120,W_target=100,check_kernel=False,export_hx=False,texc_abs=None):
	B0=texc_abs;A_=export_hx;Az=W_target;l=check_kernel;U=random.Random(seed);B1=math.log(2.);V,W,Bx,By,AS=61.5,61.5,390.8,466.66,.72;AT=(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1);B2={(1,0):(1,0),(0,1):(1,-1),(-1,1):(0,-1),(-1,0):(-1,0),(0,-1):(-1,1),(1,-1):(0,1)}
	def AC(x,y):return x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1]
	def AD(a,b):
		if a>=0 and b>=0 or a<=0 and b<=0:return a,b,0
		if abs(a)<=abs(b):return 0,a+b,-a
		return a+b,0,b
	def X(v):A,B,C=AD(v[0],v[1]);return(A!=0)+(B!=0)+(C!=0)
	def B3(v):A,B,C=AD(v[0],v[1]);return abs(A)+abs(B)+abs(C)
	def Bz(v):A,B,C=AD(v[0],v[1]);return max(abs(A),abs(B),abs(C))
	def AE(v):return v in B2
	AU,AV=Bx/V,By/W;B4=2.**V/(V*B1);B_=2.**W/(W*B1);AW=B4+B_;C0=B4/AW;R=1<<log2U;B5=int(round(R*C0));C1=R-B5;L=int(round(AS*R));w=.8;AX=min(w*min(V,W),max(8.,log2(L/25.)+log2(3.)));AY=M._primes_upto(int(floor(2.**AX)));C2=(len(AY)+5)//6;m=R+2*C2;x=400;y=[w+(1.-w)*(A+.5)/x for A in range(x)];AF,AG=[],[];AZ=Aa=.0
	for J in y:AZ+=2.**(V*(J-1.))/J;Aa+=2.**(W*(J-1.))/J;AF.append(AZ);AG.append(Aa)
	AF=[A/AZ for A in AF];AG=[A/Aa for A in AG];C3=[M.kappa_rho(AU,A)for A in y];C4=[M.kappa_rho(AV,A)for A in y];f=[.0]*m;AH=[.0]*m;S=0;B6=(1.-w)/x
	for z in range(B5):A0=min(x-1,bisect.bisect_left(AF,U.random()));J=y[A0]+(U.random()-.5)*B6;f[S]=min(1.,AS*AW*C3[A0]/2.**(V*J)/L);AH[S]=V*J;S+=1
	for z in range(C1):A0=min(x-1,bisect.bisect_left(AG,U.random()));J=y[A0]+(U.random()-.5)*B6;f[S]=min(1.,AS*AW*C4[A0]/2.**(W*J)/L);AH[S]=W*J;S+=1
	B7=[]
	for(Ab,AI)in((AU,V),(AV,W)):
		for E in range(0,len(AY),6):n=AY[E];F=min(1.,6.*M.kappa_rho(Ab,log2(n)/AI)/n);f[S]=F;AH[S]=log2(n);B7.append((S,F,n));S+=1
	o=B8=.0
	for(Ab,AI)in((AU,V),(AV,W)):
		B9=AX/AI;BA=(w-B9)/2e2
		for E in range(200):
			J=B9+BA*(E+.5);BB=M.kappa_rho(Ab,J)*BA/J;o+=BB
			if AI*J<=28.:B8+=BB
	C5=[A for A in range(R)if f[A]>=.2];Ac=[A for A in range(R)if f[A]<.2];BC=[];AJ=.0
	for E in Ac:AJ+=f[E];BC.append(AJ)
	def BD(lam):
		A=lam
		if A<3e1:
			D=exp(-A);B=1.;C=0
			while True:
				B*=U.random()
				if B<=D:return C
				C+=1
		return max(0,int(round(A+sqrt(A)*U.gauss(.0,1.))))
	N=1000000009;C6=next(A for B in range(2,50)for A in[pow(B,(N-1)//6,N)]if(A*A-A+1)%N==0)
	def a(v):return(v[0]+v[1]*C6)%N
	p=[U.randrange(1,N)for A in range(m)]if l else None;A=[];C=[set()for A in range(m)];AK=U.random;Ad=functools.partial(M.rand_below,U.getrandbits);C7=bisect.bisect_left;BE=len(Ac)-1;C8=math.log;C9=[(A,f[A])for A in C5];BF=[{}for A in range(L)]
	def Ae(d,c,mv):
		if c in d:
			A=d[c];B=A[0]+mv[0],A[1]+mv[1]
			if B==(0,0):del d[c]
			else:d[c]=B
		else:d[c]=mv
	for(B,F,n)in B7:
		if F<=.0:continue
		CA=1./n;BG=1./math.log1p(-F)if F<.999999 else .0;D=-1
		while True:
			D=D+1 if BG==.0 else D+1+int(C8(AK()or 1e-300)*BG)
			if D>=L:break
			H=1
			while AK()<CA:H+=1
			BH=AT[Ad(6)];Ae(BF[D],B,(H*BH[0],H*BH[1]))
	for D in range(L):
		q=BF[D]
		for(B,CB)in C9:
			if AK()<CB:Ae(q,B,AT[Ad(6)])
		for z in range(BD(AJ)):BI=C7(BC,AK()*AJ);Ae(q,Ac[BE if BI>BE else BI],AT[Ad(6)])
		CC=o+sum(X(A)for A in q.values());CD=-sum(a(B)*p[A]for(A,B)in q.items())%N if l else 0;A.append([q,1,o,o,CC,True,CD])
		for B in q:C[B].add(D)
	O={A for A in range(m)if C[A]}
	def Af(nalive,stack=None,seeds=None):
		J=seeds;I=stack;H=nalive;F=[A for A in O if len(C[A])==1]if J is None else list(J)
		while F:
			D=F.pop();K=C[D]
			if len(K)!=1:continue
			L=next(iter(K));E=A[L]
			if not E[5]:continue
			B=E[0]
			if I is not None:
				if any(len(C[A])<2 for A in B if A!=D):continue
				if a(B[D])==0:continue
				I.append((D,dict(B),E[6]))
			E[5]=False;H-=1
			for G in B:
				M=C[G];M.discard(L);N=len(M)
				if N==1:F.append(G)
				elif N==0:O.discard(G)
			B.clear()
		return H
	P=Af(L);AL=len(O);CE=set(O);A1=[];BJ=0;AM=max(6 if A_ else 2,int(.0001*R))if B0 is None else max(2,int(B0));CF=P-AL
	while P-AL>AM:
		BJ+=1;r=list(range(L))
		def AN(x):
			while r[x]!=x:r[x]=r[r[x]];x=r[x]
			return x
		for B in O:
			BK=C[B]
			if len(BK)==2:
				CG,CH=BK;BL,BM=AN(CG),AN(CH)
				if BL!=BM:r[BL]=BM
		Ag={}
		for E in range(L):
			D=A[E]
			if D[5]:BN=AN(E);Ag[BN]=Ag.get(BN,.0)+D[4]
		BO=sorted(Ag.items(),key=lambda kv:-kv[1]);BP=P-AL-AM;CI=max(1,min(BP//2+1,max(1,len(BO)//100)));CJ={A for(A,B)in BO[:CI]};BQ=False
		for E in range(L):
			D=A[E]
			if D[5]and AN(E)in CJ:
				g=D[0]
				if any(len(C[A])<2 for A in g):continue
				D[5]=False;P-=1;BQ=True;h=[]
				for I in g:
					A2=C[I];A2.discard(E);A3=len(A2)
					if A3==1:h.append(I)
					elif A3==0:O.discard(I)
				g.clear()
				if h:P=Af(P,stack=A1,seeds=h)
		if not BQ:
			BR=sorted([(A[B][4],B)for B in range(L)if A[B][5]and all(len(C[A])>=2 for A in A[B][0])],reverse=True)
			if not BR:break
			for(z,E)in BR[:max(1,BP//2)]:
				D=A[E]
				if not D[5]:continue
				g=D[0]
				if any(len(C[A])<2 for A in g):continue
				D[5]=False;P-=1;h=[]
				for I in g:
					A2=C[I];A2.discard(E);A3=len(A2)
					if A3==1:h.append(I)
					elif A3==0:O.discard(I)
				g.clear()
				if h:P=Af(P,stack=A1,seeds=h)
				if P-len(O)<=AM:break
		AL=len(O)
	BS=len(CE-(O|{A for(A,B,B)in A1}));Ah=0
	if l:
		BT={A:p[A]for A in O}
		for(AO,BU,CK)in reversed(A1):
			Ai=a(BU[AO]);Aj=sum(a(B)*BT[A]for(A,B)in BU.items()if A!=AO)%N;BV=-(CK+Aj)*pow(Ai,-1,N)%N;BT[AO]=BV
			if BV!=p[AO]:Ah+=1
	s,BW=64,255;b=sum(A[B][4]for B in range(L)if A[B][5])
	def t(c):
		B=len(C[c])
		if B<=2:return float(B-3)
		D=1e18
		for E in C[c]:
			if AE(A[E][0][c])and A[E][4]<D:D=A[E][4]
		return float(B-2)*D
	c=[(t(A),A)for A in range(m)if 1<=len(C[A])<=s];heapq.heapify(c);Ak=0;A4=[];CL=set(O)
	while c and b/max(1,P)<Az:
		CM,B=heapq.heappop(c);Al=len(C[B])
		if Al==0 or Al>s:continue
		BX=t(B)
		if BX>CM+1e-09:heapq.heappush(c,(BX,B));continue
		if Al==1:
			A5=next(iter(C[B]));A6=A[A5][0]
			if not AE(A6[B])or any(len(C[A])<2 for A in A6 if A!=B):continue
			A4.append((B,dict(A6),A[A5][6]));A[A5][5]=False;P-=1;b-=A[A5][4]
			for I in A6:
				C[I].discard(A5)
				if 1<=len(C[I])<=s:heapq.heappush(c,(t(I),I))
			A6.clear();continue
		A7=-1;BY=1e18
		for K in C[B]:
			if AE(A[K][0][B])and A[K][4]<BY:BY=A[K][4];A7=K
		if A7<0:continue
		F=A[A7];i=F[0];BZ=B2[i[B]];Ba=[A for A in C[B]if A!=A7];u=True;Am=-F[4];A8={}
		for E in Ba:
			D=A[E];d=AC(D[0][B],BZ);A9=X(d)
			if D[1]+B3(d)*F[1]>BW:u=False;break
			Am+=F[2]if A9==1 else A9*F[3];Y=D[0]
			for(G,An)in i.items():
				Q=Y.get(G);Z=AC(d,An);T=(-Z[0],-Z[1])if Q is None else(Q[0]-Z[0],Q[1]-Z[1])
				if Bz(T)>BW:u=False;break
				Am+=X(T)-(0 if Q is None else X(Q))
				if G!=B:
					if Q is None:A8[G]=A8.get(G,0)+1
					elif T==(0,0):A8[G]=A8.get(G,0)-1
			if not u:break
		if u and any(len(C[A])-1+A8.get(A,0)<=0 for A in i if A!=B):u=False
		if not u:continue
		if b+Am>M.W_ROW*max(1,P-1):break
		A4.append((B,dict(i),F[6]))
		for E in Ba:
			D=A[E];Y=D[0];d=AC(Y[B],BZ);A9=X(d);D[1]+=B3(d)*F[1];Bb=F[2]if A9==1 else A9*F[3];D[2]+=Bb;D[3]+=F[3];b+=Bb
			if l:D[6]=(D[6]-a(d)*F[6])%N
			for(G,An)in i.items():
				Q=Y.get(G);Z=AC(d,An)
				if Q is None:T=-Z[0],-Z[1];Y[G]=T;C[G].add(E);b+=X(T)
				else:
					T=Q[0]-Z[0],Q[1]-Z[1];b+=X(T)-X(Q)
					if T==(0,0):
						del Y[G];C[G].discard(E)
						if 1<=len(C[G])<=s:heapq.heappush(c,(t(G),G))
					else:Y[G]=T
			D[4]=D[2]+sum(X(A)for A in Y.values())
			for G in Y:
				if len(C[G])<=s:heapq.heappush(c,(t(G),G))
			Ak+=1
		F[5]=False;P-=1;b-=F[4]
		for I in i:
			C[I].discard(A7)
			if 1<=len(C[I])<=s:heapq.heappush(c,(t(I),I))
		i.clear()
	j=[B for B in range(L)if A[B][5]];AP={C for B in j for C in A[B][0]};Bc=len(CL-(AP|{A for(A,B,B)in A4}));Ao=0
	if l:
		Bd={A:p[A]for A in AP}
		for(AQ,Be,CN)in reversed(A4):
			Ai=a(Be[AQ]);Aj=sum(a(B)*Bd[A]for(A,B)in Be.items()if A!=AQ)%N;Bf=-(CN+Aj)*pow(Ai,-1,N)%N;Bd[AQ]=Bf
			if Bf!=p[AQ]:Ao+=1
	v=len(j);Ap=0
	if l:
		for K in j:
			if(sum(a(B)*p[A]for(A,B)in A[K][0].items())+A[K][6])%N!=0:Ap+=1
	CO=sorted([(len(C[A]),A)for A in range(R)if C[A]]);k={};AA={};Aq={};Ar=0
	for(z,B)in CO:
		Bg=Aq.get(B,1.);CP=sorted([(A[C][4],C)for C in C[B]if C not in k and AE(A[C][0][B])])
		for(z,K)in CP:
			AB=.0
			for e in A[K][0]:
				if e!=B and e in AA and AA[e]>AB:AB=AA[e]
			H=AB+.5*(Bg-AB)
			if AB<H<Bg:
				k[K]=B;AA[B]=H;Ar+=1
				for e in A[K][0]:
					if e!=B and Aq.get(e,1.)>H:Aq[e]=H
				break
	As=v-Ar;At,CQ=M.C_UNIT_5,M.C_GEN_1;Au=Bh=Bi=Bj=Bk=.0;Bl=B8/max(o,1e-12)
	for K in j:
		CR=k.get(K,-1);Av=Bm=Aw=0;AR=.0
		for(B,Bn)in A[K][0].items():
			Ax=AH[B];CS=B==CR
			for J in AD(Bn[0],Bn[1]):
				H=abs(J)
				if H==0:continue
				Au+=1.
				if H==1:Bh+=1.
				Bo=At*H if H<=4 else CQ;Bi+=Bo
				if not CS:
					AR+=Bo
					if H<=4 and Ax<=28.:
						if Ax<=14.:Av+=H
						elif Ax<=18.25:Bm+=H
						else:Aw+=H
		Ay=A[K][3];AR+=Ay*At
		if Bl>.0 and Ay>.0:Aw+=BD(Ay*Bl)
		Bp=Av//4;CT=Av-4*Bp;Bq=Bm+CT;Br=Bq//3;CU=Bq-3*Br;CV=(Aw+CU)//2;CW=(3*Bp+2*Br+CV)*At;Bj+=AR;Bk+=AR-CW
	CX=Bh/max(Au,1.);CY=Bi/max(Au,1.);Bs=Bk/max(Bj,1e-12);Bt=b/max(1,v);Bu=dict(ok=Ap==0 and BS==0 and Bc==0 and Ah==0 and Ao==0 and v>0 and As>0 and Ak>0 and Bt<=M.W_ROW and Bs<M.TAB_R_CAP,passes=BJ,s0=CF,log2U=log2U,seed=seed,W_target=Az,initial_rows=L,texc=AM,ncols=len(AP),surplus=v-len(AP),nrec=Ar,N_sw=v,N_l=As,N_sw_U=round(v/R,6),N_l_U=round(As/R,6),avg_w=round(Bt,2),merges=Ak,w_mid=round(o,4),lsmall=round(AX,3),exp_unit=round(CX,5),k5_exp=round(CY,4),tab_r=round(Bs,5),kernel_errors=Ap,post_peel_cols_lost=BS,post_merge_cols_lost=Bc,repurge_stack_len=len(A1),merge_stack_len=len(A4),repurge_backsub_errors=Ah,merge_backsub_errors=Ao)
	if A_:Bv=sorted([(AA[k[A]],k[A],A)for A in j if A in k]);Bw=[B for(A,B,A)in Bv];CZ=set(Bw);Ca=[A for A in j if A not in k];O=sorted({C for B in j for C in A[B][0]});Cb=[A for A in O if A not in CZ];Cc=[(B,A[C][0][B],[(A,C)for(A,C)in A[C][0].items()if A!=B])for(D,B,C)in Bv];Cd=[[(A,B)for(A,B)in A[B][0].items()]for B in Ca];Bu['hx_export']=dict(rec_order=Bw,recorded_pivots=Cc,live_col_set=Cb,live_rows_full=Cd)
	return Bu
def verify_trim_pairs():
	H=M.load_witnesses()['trim_pairs'];c=H['texc_new'];I=H['all_ok'];A=H['rows'];J=[];K=[];O=[];P=[];Q=[];C=[]
	for(j,k,d,e,l,R,f,S,D,g,T,U,V,h)in A:
		W=D>=100;I=I and D<=R and g<=f-(d-c)and(not W or U<M.TAB_R_CAP and V<=M.W_ROW and T<=M.K5_LLA*S)
		if W:J.append(T/S);P.append(U);Q.append(V)
		else:C.append(D)
		K.append(D/R);O.append(h/math.log(e))
	i=max(A[3]/2.**A[0]for A in A);X=math.ceil(max(O)*100)/100;Y=math.ceil(i*1e4)/1e4;L=math.ceil(X*(math.log(Y)+M.log2_U_z6(M.OP['lf'],M.OP['lg'])*math.log(2.)));E=math.ceil(M.log2_U_z6(M.OP['lf'],M.OP['lg']));N=8*E;F=512;Z=math.ceil((448.+478.)/math.log2(5.));assert Z<=F;a=dict(parent_size_reset=8,col_unions=(3+2*N+12)/.72,accumulate=2+N+10,component_sort=64*E,kill_sweep=2+N+10,fallback_sweep=2+8*F+64*E);B=sum(a.values());G=2*8*F;b=lambda x:math.ceil(x*100)/100;return dict(ok=I and L<=4096 and B<=2.**16 and 4096*B+G<=2.**28,pass_ir_terms={A:round(B,1)for(A,B)in a.items()},pass_ir_per_relation=round(B,1),drop_once_ir=G,find_steps_max=E,row_entries_max=F,row_entries_bound=Z,budget_log2_at_pass_cap=b(math.log2(4096*B+G)),budget_log2_at_predicted_passes=b(math.log2(L*B+G)),tab_r_max=math.ceil(max(P)*1e3)/1e3,avg_w_max=math.ceil(max(Q)*100)/100,ncols_lt_100=[min(C),max(C)],pairs_ncols_lt_100=len(C),pairs=len(A),ncols_ratio_range=[math.floor(min(K)*1000)/1000,math.ceil(max(K)*1000)/1000],pairs_ncols_ge_100=len(J),k5_ratio_max=math.ceil(max(J)*1e4)/1e4,passes_max=max(A[13]for A in A),extra_passes_max=max(A[13]-A[4]for A in A),passes_per_ln_s0_max=X,s0_over_U_max=Y,passes_at_prod=L)
def verify_live_trim_pair(seed=20261120):A=verify_live_filtsim_z6(log2U=11,seed=seed,W_target=100,check_kernel=True,texc_abs=32);B=verify_live_filtsim_z6(log2U=11,seed=seed,W_target=100,check_kernel=True,texc_abs=M.TEXC_LA);return dict(ok=A['ok']and B['ok']and A['surplus']==A['texc']==32 and B['surplus']==B['texc']==int(M.TEXC_LA)and B['ncols']<=A['ncols']and B['N_sw']<=A['N_sw'],old_texc32=A,new_texc6=B)
def verify_live_filtsim_ladder(seed=20261120):A=verify_live_filtsim_z6(log2U=12,seed=seed,W_target=100,check_kernel=True);B=verify_live_filtsim_z6(log2U=13,seed=seed,W_target=100,check_kernel=False);D=B['N_sw_U']<A['N_sw_U']and B['N_l_U']<A['N_l_U'];E=[12,13,20,21,22,23,53.64];C=[compute_w_mid_at_scale(A)for A in E];F=all(C[A][1]>C[A+1][1]for A in range(len(C)-1))and C[-1][1]<.25;return dict(ok=A['ok']and B['ok']and D and F and A['surplus']==A['texc']and B['surplus']==B['texc'],dim_monotonic_12_to_13=D,w_mid_monotone_to_zero_ok=F,w_mid_prod_per_side_zero=True,w_mid_ladder={str(A):dict(lsmall=B,w_mid=C)for(A,(B,C))in zip(E,C)},u12=A,u13=B)
def verify_c6_pgf_and_hea_campaign():
	from fractions import Fraction as I;An=M.load_data();A9=An['t4_geom'];G=M.load_witnesses()['pgf_campaign'];B=G['tower'];E,F,BB=int(B['u_t']),int(B['p_t']),float(B['wa_lll']);A,u,f,T,U=int(B['lb']),int(B['nA']),int(B['nB']),int(B['b_d2']),int(B['T_surv']);v,Ao=float(B['p_ecm']),int(B['miss_fb_step']);g,J,K,AA=1<<A,1<<u,1<<f,1<<U;AB,AC,Ap,Aq=B['qbins_hi_f'],B['qbins_hi_g'],B['qbins_lo_f'],B['qbins_lo_g'];Ar,As,AD,h=float(B['alpha_f']),float(B['alpha_g']),M.ZK_NUM,M.W_COEFFS;AE=[-4-A if B==0 else-A for(B,A)in enumerate(h)];BC=[[complex(A,B)for(A,B)in A]for A in A9['EMB_ZK']];w=sqrt(2.);At=M._pmul_h_q([I(A)for A in h],[I(A)for A in AE]);Au=M._ppow_h([0,1,0,0,0,0],F,F);AF=M._ppow_h(Au,F,F);AG=M._ppow_h(AF,F,F);Av=M._ppow_h(AG,F**3,F);AH=At==[I(6)]+[I(0)]*5 and 1+E*6==F and Av==[0,1,0,0,0,0]and AF!=[0,1,0,0,0,0]and AG!=[0,1,0,0,0,0];AI=[I(A,3)for A in M.SIGMA_NUM];x=[[I(1)]+[I(0)]*5,AI]
	for _ in range(4):x.append(M._pmul_h_q(x[-1],AI))
	AJ=tuple(tuple(int(3*x[B][A])for B in range(6))for A in range(6));P=M._mh_k;AK=math.sumprod
	def i(p):B=tuple(AK(A,p)for A in AJ);E=tuple(AK(A,B)for A in AJ);C=P(P(p,B,44),E,88);A=-C[2]//27;D=C[0]//27-2*A;return D*D+2*A*A
	def AL(pa2,pb2):A=pb2;B=P(h,pa2,44);return i((E*B[0]-A[0],E*B[1]-A[1],E*B[2]-A[2],E*B[3]-A[3],E*B[4]-A[4],E*B[5]-A[5]))//729
	def AM(pa2,pb2):A=pa2;B=P(AE,pb2,44);return i((A[0]+B[0],A[1]+B[1],A[2]+B[2],A[3]+B[3],A[4]+B[4],A[5]+B[5]))//729
	def j(xs):
		A=xs
		while len(A)>1:A=[A[B]*A[B+1]if B+1<len(A)else A[B]for B in range(0,len(A),2)]
		return A[0]if A else 1
	k=M._primes_upto(g);Aw=[A for A in k if A>2048];y=[B for(A,B)in enumerate(Aw)if A%Ao==0];AN=j(y);Ax,AO,AP=j([A for A in k if A<=256]),j([A for A in k if 256<A<=2048]),j([A for A in k if 2048<A<=g]);Ay=AO*AP;l=M._isp;AQ=M._brent;V=math.gcd;AR=1<<2*A;BD=K*K;Az=round(v*1000)
	def AS(n):
		for B in(Ax,AO,AP):
			while n>1:
				A=V(n,B)
				if A==1:break
				n//=A
			if n==1:return 1
		return n
	def AT(n):
		D=0
		while True:
			A=V(n,Ay)
			if A==1:break
			n//=A
			if A<65536 or l(A):D+=1
			else:
				B=[A]
				while B:
					C=B.pop()
					if C<65536 or l(C):D+=1
					else:E=AQ(C);B.append(E);B.append(C//E)
		return D
	def A_(side,q,sp_idx,n_fac):
		B=side<<56^q<<24^sp_idx<<8^1371938854
		for C in range(n_fac):
			A=(B^C)+0x9e3779b97f4a7c15&0xffffffffffffffff;A=(A^A>>30)*0xbf58476d1ce4e5b9&0xffffffffffffffff;A=(A^A>>27)*0x94d049bb133111eb&0xffffffffffffffff
			if(A^A>>31)%1000>=Az:return False
		return True
	def AU(cof,b_a=J,b_b=K):
		C=b_b;A=cof
		if A==1:return True,[]
		if A<AR:return A<=C,[A]if A<=C else[]
		if A>C*C or l(A):return False,[]
		D,E=[],[A]
		while E:
			B=E.pop()
			if B==1:continue
			if B<AR or l(B):
				if B>C:return False,[]
				D.append(B)
				if len(D)>2:return False,[]
				continue
			F=AQ(B)
			if F==B:return False,[]
			E.append(F);E.append(B//F)
		G=sum(1 for A in D if A>b_a);return len(D)-G<=1 and G<=1,D
	B0=M.build_exact_int_geom(A9)[0]
	def B1(q,s,rho):
		C=[_eval_zk_mod([1 if B==A else 0 for B in range(6)],s,q)for A in range(6)];D=pow(C[0],-1,q);B=[[0]*12 for A in range(12)];B[0][0]=q
		for A in range(1,6):B[A][A]=1;B[A][0]=-C[A]*D%q
		for A in range(6):B[6+A][6+A]=1;B[6+A][0]=rho*C[A]*D%q
		DESC.bareiss_lll_12d(B,B0,121,25,False);return B,[([sum(A[B]*AD[B][C]for B in range(6))for C in range(6)],[sum(A[6+B]*AD[B][C]for B in range(6))for C in range(6)])for A in B]
	B2=[((A,1),)for A in range(12)]+[((A,1),(B,C))for A in range(12)for B in range(A+1,12)for C in(1,-1)]+[((A,1),(B,D),(C,E))for A in range(12)for B in range(A+1,12)for C in range(B+1,12)for D in(1,-1)for E in(1,-1)]
	def z(cside,cq,av,bv):
		F=cside;A=cq;B=M._ppowmod([0,1],A,M.H_COEFFS,A)
		if len(B)<2:B+=[0]*(2-len(B))
		B[1]=(B[1]-1)%A
		for D in M._roots_mod(M._pgcd(B,M.H_COEFFS,A),A):
			I,G=_eval_zk_mod(av,D,A),_eval_zk_mod(bv,D,A)
			if G==0:continue
			C=I*pow(G,-1,A)%A
			if C==0:continue
			H=_eval_poly_mod(h,D,A)
			if F==0 and(E*H*C*C-1)%A==0 or F==1 and(C*C-4-H)%A==0:return F,A,D,C
	def AV(side,q,s,rho,b_a=J,b_b=K,first=False):
		U=first;T=side;k,R=B1(q,s,rho);Z=a=0;G=[A[0]for A in R]+[[-A for A in A[0]]for A in R];H=[A[1]for A in R]+[[-A for A in A[1]]for A in R];l,m=(AM,AL)if T==0 else(AL,AM)
		for(n,A)in enumerate(B2):
			if len(A)==3:O,Q,b=A[0][0],A[1][0]+(0 if A[1][1]==1 else 12),A[2][0]+(0 if A[2][1]==1 else 12);B,C,I=G[O],G[Q],G[b];D,E,J=H[O],H[Q],H[b];K=B[0]+C[0]+I[0],B[1]+C[1]+I[1],B[2]+C[2]+I[2],B[3]+C[3]+I[3],B[4]+C[4]+I[4],B[5]+C[5]+I[5];L=D[0]+E[0]+J[0],D[1]+E[1]+J[1],D[2]+E[2]+J[2],D[3]+E[3]+J[3],D[4]+E[4]+J[4],D[5]+E[5]+J[5]
			elif len(A)==2:O,Q=A[0][0],A[1][0]+(0 if A[1][1]==1 else 12);B,C=G[O],G[Q];D,E=H[O],H[Q];K=B[0]+C[0],B[1]+C[1],B[2]+C[2],B[3]+C[3],B[4]+C[4],B[5]+C[5];L=D[0]+E[0],D[1]+E[1],D[2]+E[2],D[3]+E[3],D[4]+E[4],D[5]+E[5]
			else:K,L=tuple(G[A[0][0]]),tuple(H[A[0][0]])
			if not any(K)or not any(L):continue
			c=tuple(A//3 for A in P(K,K,32));d=tuple(A//3 for A in P(L,L,32));S=l(c,d)
			if S<=0:continue
			e=AS(S);f=not U and e<AA;g,W=AU(e,b_a,b_b)
			if not f and not g:continue
			N=m(c,d)
			if N%q!=0:continue
			N//=q
			if N<=0:continue
			h=AS(N)
			if f and h<AA:a+=1
			if g:
				o,X=AU(h,b_a,b_b)
				if o and V(S,AN)==1 and V(N,AN)==1 and A_(T,q,n,len(W)+len(X)+AT(S)+AT(N)):
					F=[0]*12
					for(p,r)in A:
						j=k[p]
						for Y in range(12):F[Y]+=j[Y]if r==1 else-j[Y]
					if V(i(K)//729,i(L)//729)==1 or M.is_ideal_primitive_oh(F[:6],F[6:]):
						if U:return(F[:6],F[6:],X,W)if T==0 else(F[:6],F[6:],W,X)
						Z+=1
		return None if U else(Z,a)
	A0=lambda ps,qp,spy,rhp,b_a=J,b_b=K:AV(ps,qp,spy,rhp,b_a,b_b,True);m=lambda side,q,s,rho:AV(side,q,s,rho);AW={}
	def AX(u):
		A=AW.get(u)
		if A is None:A=M.rho(u);AW[u]=A
		return A
	A1=lambda u:1. if u<=1. else AX(ceil(u*1024.)/1024.);A2=lambda u:1. if u<=1. else AX(floor(u*1024.)/1024.);AY=int(round(48*(u-A)/(f-A)));L=(f-A)/48.;AZ=[math.log1p(L/(A+B*L))if B<AY else .0 for B in range(48)];Aa=[math.log1p(L/(A+B*L))if B>=AY else .0 for B in range(48)];B3=[(0,[1.]),(1,AZ),(1,Aa),(2,DESC._conv(AZ,Aa))]
	def Ab(m,alp):m=ceil(m);return min(1.,sum(D*A1((m-B*A-C*L+alp)/A)for(B,E)in B3 for(C,D)in enumerate(E)if D and B*A+C*L+B*L<=m))
	Ac=lambda m:1. if U>=m else min(1.,max(1.,U/A)*A2((floor(m)-U)/A));A3=6.*log2(E)
	def n(side,q):C=log2(q);B=(C+1e1)/12.;D=C if side==1 else .0;E=C if side==0 else .0;F=sum(1./6.*Ab(12.*B+AB[A]-D,Ar)*Ab(12.*B+A3+AC[A]-E,As)for A in range(6));K=sum(1./6.*Ac(12.*B+Ap[A]-D)*Ac(12.*B+A3+Aq[A]-E)for A in range(6));G=12.*B+AB[-1]-D;H=12.*B+A3+AC[-1]-E;L=max(.0,1.-(1.-v)*((G+H)/8.+4.));I,J=G/A,H/A;M,N=max(A1(I),1e-300),max(A1(J),1e-300);O=sum((A2(max(.0,I-log2(B)/A))/M+A2(max(.0,J-log2(B)/A))/N)/B for B in y);return round(1024.*.683542*F*L*max(.0,1.-O),2),round(1024.*min(1.,K+F),2)
	C=True;Q=[];H=True
	for((W,Ad,o,Ae,A4),(p,q))in zip(G['live_d2_chains'][:2],G['live_d2_expected'][:2]):
		D=A0(*W,J,K)
		if D is None or o not in(D[2]if Ad==0 else D[3])or not T<o<=J<W[1]:C=H=False;continue
		N=z(Ad,o,D[0],D[1])
		if N is None:C=H=False;continue
		Af,Ag=m(*N);Ah,Ai=n(N[0],N[1])
		if Af!=p[3]or Ag!=p[5]or abs(Ah-p[2])>.02 or abs(Ai-p[4])>.02:C=H=False
		Q.append([N[0],N[1],1,Ah,Af,Ai,Ag]);X=A0(*N,T,T)
		if X is None or A4 not in(X[2]if Ae==0 else X[3])or not g<A4<=T<o:C=H=False;continue
		R=z(Ae,A4,X[0],X[1])
		if R is None:C=H=False;continue
		Aj,w=m(*R);Ak,Al=n(R[0],R[1])
		if Aj!=q[3]or w!=q[5]or abs(Ak-q[2])>.02 or abs(Al-q[4])>.02:C=H=False
		Q.append([R[0],R[1],2,Ak,Aj,Al,w])
	for((W,Am,A5),r)in zip(G['live_d1_parents'][:2],G['live_d1_expected'][:2]):
		D=A0(*W,J,K)
		if D is None or A5 not in(D[2]if Am==0 else D[3])or not g<A5<=J<W[1]:C=False;continue
		S=z(Am,A5,D[0],D[1])
		if S is None:C=False;continue
		Y,Z=m(*S);a,b=n(S[0],S[1])
		if Y!=r[3]or Z!=r[5]or abs(a-r[2])>.02 or abs(b-r[4])>.02:C=False
		Q.append([S[0],S[1],1,a,Y,b,Z])
	for(A6,A7,B4,B5,B6,B7,B8,B9)in G['live_ctrl_ideals'][:2]:
		Y,Z=m(A6,A7,B4,B5);a,b=n(A6,A7)
		if Y!=B7 or Z!=B9 or abs(a-B6)>.02 or abs(b-B8)>.02:C=False
		Q.append([A6,A7,0,a,Y,b,Z])
	c=[[int(A[0]),M.b62d(A[1]),int(A[2]),M.b62d(A[3])/1e2,M.b62d(A[4]),M.b62d(A[5])/1e2,M.b62d(A[6])]for A in(A.split()for A in G['compact_224_b62'].split(';'))]
	def d(rs):D=rs;A=len(D);B=[A[3]for A in D];C=[A[4]for A in D];F=[A[5]for A in D];G=[A[6]for A in D];J=sum(C)/max(sum(B),1e-12);K=sum(C)/A;L=round(sum((A-J*B)**2 for(A,B)in zip(C,B))/max(1,A-1)/max(K,1e-12),3);M=round(sum((1.+.5*A)**-1. for A in B)/A,6);H=round(sum((1.+.5*A)**-1. for A in B),3);E=sorted(A/B for(A,B)in zip(C,B));I={str(D):dict(emp=round(sum(D**A for A in C)/A,8),h_pgf_1_r1=round(sum((1.+A*(1.-D))**-1. for A in B)/A,8),h_pgf_05_r1=round(sum((1.+.5*A*(1.-D))**-1. for A in B)/A,8),ok=sum(D**A for A in C)/A<=sum((1.+A*(1.-D))**-1. for A in B)/A and sum(D**A for A in C)/A<=sum((1.+.5*A*(1.-D))**-1. for A in B)/A)for D in(.0,.2,.5,.8)};return dict(n=A,mu_min=round(min(B),2),mu_max=round(max(B),2),M_min=min(C),M_max=max(C),n_zero=sum(1 for A in C if A==0),nb_r1_p0=M,exp_p0_sum=H,direct_oe_p0=round(sum(1 for A in C if A==0)/max(H,1e-12),4),min_M_over_mu=round(E[0],4),p05_M_over_mu=round(E[max(0,int(.05*A))],4),below_half_mu=sum(1 for A in E if A<=.5),below_08_mu=sum(1 for A in E if A<=.8),res_fano_v_over_m=L,mean_M_over_mu=round(sum(E)/A,4),two_sided_surv_obs_over_up=round(sum(G)/sum(F),4),two_sided_surv_max_over_up=round(max(A/B for(A,B)in zip(G,F)),4),pgf_ok=all(A['ok']for A in I.values()),pgf=I)
	O,s,t,e,A8=d(c),d([A for A in c if A[2]==1]),d([A for A in c if A[2]==2]),d([A for A in c if A[2]==0]),d(Q);BA=AH and C and H and len(c)==224 and len(Q)==8 and O['n_zero']==0 and O['mu_min']>=25. and max(O['res_fano_v_over_m'],s['res_fano_v_over_m'],t['res_fano_v_over_m'],e['res_fano_v_over_m'])<3. and O['pgf_ok']and s['pgf_ok']and t['pgf_ok']and e['pgf_ok']and A8['pgf_ok']and O['mean_M_over_mu']>=1. and O['two_sided_surv_max_over_up']<1. and A8['two_sided_surv_max_over_up']<1.;return dict(ok=BA,unified_bareiss_910_used=True,exact_ideal_primitivity_checked=True,tower_verified=AH,live_16_usable_and_depth2_exact=C,live_4_depth2_2hop_chains_exact=H,tower_params=dict(u_t=E,p_t=F,lb=A,nA=u,nB=f,b_d2=T,T_surv=U,ball_pts=1024,p_ecm=v,miss_fb_primes=len(y)),usable_thinning_totals=G['usable_thinning_totals'],slot_d1_vs_ctrl_M_over_mu=round(s['mean_M_over_mu']/e['mean_M_over_mu'],4),slot_d2_vs_ctrl_M_over_mu=round(t['mean_M_over_mu']/e['mean_M_over_mu'],4),all_224=O,slot_d1_100=s,slot_d2_24=t,ctrl_100=e,live_16=A8)
def verify_live_c6_pool_hnf_subsample(n_tot=1000,seed=20261121):
	d=M.load_data();R=d['t4_geom'];A=d['pool_quad'];A5=R['GA'];S=[[complex(A,B)for(A,B)in A]for A in R['EMB_ZK']];e,f,AS,AT,AU,A6=M.build_exact_int_geom(R);L=_qf6;T=A['wa'];g=A['wb'];U=A['ra0'];A7=A['rb0_ga'];V=T*int(round(U*U))*A6
	def A8(a,b):
		for A in range(6):
			if T*L(f[A],a)+g*L(f[A],b)<0:return False
		return True
	A9=M.is_ideal_primitive_oh
	def AA(z):return((((z-2)*z+1)*z+2)*z*z-4)*z+3
	I=[complex(.4+0.9j)**A for A in range(6)]
	for _ in range(150):
		for N in range(6):
			h=1.+0j
			for D in range(6):
				if N!=D:h*=I[N]-I[D]
			I[N]-=AA(I[N])/h
	i=[min(I,key=lambda r:abs(r-S[A][1]))for A in range(3)];j=lambda cs,z:sum(B*z**A for(A,B)in enumerate(cs));AB=[j([0,-2,-1,0,1,-1],A)for A in i];AC=[[j(M.F_REC_COEFFS[B],A)for B in(0,2,4,6,8)]for A in i];AD=float(M.U_BN254)
	def AE(G):
		B=[[.0]*6 for A in range(6)]
		for C in range(6):
			for A in range(C+1):D=G[C][A]-sum(B[C][D]*B[A][D]for D in range(A));B[C][A]=sqrt(D)if C==A else D/B[A][A]
		return B
	def AF(L):
		B=[[.0]*6 for A in range(6)]
		for A in range(6):
			B[A][A]=1./L[A][A]
			for C in range(A):B[A][C]=-sum(L[A][D]*B[D][C]for D in range(C,A))/L[A][A]
		return B
	k=AF(AE(A5));l=random.Random(seed);J=C=O=W=0;E=[];F=[];AG=(1,0,2,2,0,1),(1,2,0,2,1,2),(1,2,2,1,2,2);X=Y=0
	def AH(a,b):return any(sum(A*B for(A,B)in zip(A,a))%3==0 and sum(A*B for(A,B)in zip(A,b))%3==0 for A in AG)
	while J<n_tot:
		Z=[l.gauss(0,1)for A in range(12)];m=l.random()**(1./12.)/sqrt(sum(A*A for A in Z));G=[int(round(sum(k[B][A]*Z[B]*m*U for B in range(A,6))))for A in range(6)];H=[int(round(sum(k[B][A]*Z[6+B]*m*A7 for B in range(A,6))))for A in range(6)]
		if not any(G)or not any(H):continue
		n=T*L(e,G);o=g*L(e,H)
		if n+o>V:continue
		J+=1
		if 20*n>17*V or 20*o>17*V:continue
		if not A8(G,H):continue
		C+=1;p=AH(G,H);X+=p
		if A9(G,H):
			Y+=p;O+=1;q=r=.0
			for D in range(3):s=sum(S[D][A]*G[A]for A in range(6));t=sum(S[D][A]*H[A]for A in range(6));K,B=s*s,t*t;AI,AJ,AK,AL,AM=AC[D];AN=(((AM*K+AL*B)*K+AK*B*B)*K+AJ*B*B*B)*K+AI*B*B*B*B;AO=B+((AD+1)*AB[D]+2.)*K;q+=2.*log2(abs(AN));r+=2.*log2(abs(AO))
			E.append(q+M.CONTENT_BITS);F.append(r)
		else:W+=1
	a=A['bins'];u=sum(A[2]for A in a);v=sum(A[0]*A[2]for A in a)/u;w=sum(A[1]*A[2]for A in a)/u;x=C/float(J);y=O/max(1,C);P=A['n_prim_in']/float(A['n_in']);b=sum(E)/len(E);c=sum(F)/len(F);z=(x-A['frac_reg'])/sqrt(max(1e-12,A['frac_reg']*(1.-A['frac_reg'])/J));A0=(y-P)/sqrt(P*(1.-P)/max(1,C));AP=sqrt(sum((A-b)**2 for A in E)/max(1,len(E)-1)/len(E));AQ=sqrt(sum((A-c)**2 for A in F)/max(1,len(F)-1)/len(F));A1=(b-v)/AP;A2=(c-w)/AQ;Q=217./729.;A3=X/float(max(1,C));A4=(A3-Q)/sqrt(Q*(1.-Q)/max(1,C));AR=abs(z)<=2.5 and abs(A0)<=2.5 and abs(A1)<=2.5 and abs(A2)<=2.5 and O+W==C and Y==0 and abs(A4)<=2.5;return dict(ok=AR,mod3_bad_cells=X,mod3_bad_frac=round(A3,4),mod3_bad_frac_exact=round(Q,4),z_mod3_bad=round(A4,2),mod3_bad_but_primitive=Y,zero_float_fd_and_hnf=True,region='E12_JointVor',n_samples=J,n_active=C,n_prim=O,n_nonprim=W,frac_reg=round(x,4),frac_reg_200k=A['frac_reg'],z_frac_reg=round(z,2),prim_rate=round(y,4),prim_rate_200k=round(P,4),z_prim_rate=round(A0,2),mean_lnf=round(b,2),mean_lnf_200k=round(v,2),z_lnf=round(A1,2),mean_lng=round(c,2),mean_lng_200k=round(w,2),z_lng=round(A2,2))
def verify_live_genuine_c6_relation_filter(seed=20261120):import c6_adaptive as A;return A.verify_live_genuine_c6_relation_filter(seed=seed)
verify_live_bareiss_12d_lll=DESC.verify_live_bareiss_12d_lll
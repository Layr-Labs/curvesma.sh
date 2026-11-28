F='files'
AA=dict
A9=list
A8=sorted
w=int
l=True
k=False
e=all
W=None
U=enumerate
H=pow
G=sum
C=len
A=range
import json,math as AM,random,sys
from pathlib import Path
B=Path(__file__).resolve().parent
J=B/F if(B/F).is_dir()else B
for D in(str(J),str(B)):
	if D not in sys.path:sys.path.insert(0,D)
import verify_bn254_v6 as V
I=1000033
S=next(A for B in A(2,50)for A in[H(B,(I-1)//6,I)]if(A*A-A+1)%I==0)
def O(v,chi=S,q=I):return(v[0]+v[1]*chi)%q
def v(rows,ncols,q):
	F=ncols;B=[a[:]for a in rows];J=C(B);K=[];D=0
	for E in A(F):
		L=next((A for A in A(D,J)if B[A][E]),W)
		if L is W:continue
		B[D],B[L]=B[L],B[D];S=H(B[D][E],q-2,q);M=[c*S%q for c in B[D]];B[D]=M
		for G in A(J):
			if G!=D and B[G][E]:
				T=B[G][E];P=B[G]
				for I in A(E,F):
					if M[I]:P[I]=(P[I]-T*M[I])%q
		K.append(E);D+=1
		if D==J:break
	V=set(K);X=[b for b in A(F)if b not in V];Q=[]
	for N in X:
		O=[0]*F;O[N]=1
		for(R,Y)in U(K):
			if B[R][N]:O[Y]=-B[R][N]%q
		Q.append(O)
	return D,Q
def Ac(rows,q):
	A=rows
	if not A:return 0
	return v(A,C(A[0]),q)[0]
def unplanted_two_sided_5sm_bw_check(seed=1,m_blk=10,n_blk=8):
	L=m_blk;J=n_blk;O,B,h=1187,593,534;AB=B*B;m=[3,1,-7,0,1];x=[-1,-6,-2,1];R=5
	def n(c,x,mod):
		A=0
		for B in reversed(c):A=(A*x+B)%mod
		return A
	assert n(m,h,O)==0 and n(x,h,O)==0
	def Ad(a,b,mod_p,q):
		G=mod_p;E=C(G)-1;B=[0]*(C(a)+C(b)-1)
		for(H,I)in U(a):
			if I:
				for(D,K)in U(b):B[H+D]=(B[H+D]+I*K)%q
		for F in A(C(B)-1,E-1,-1):
			J=B[F]
			if J:
				for D in A(E):B[F-E+D]=(B[F-E+D]-J*G[D])%q
		return B[:E]
	def BK(base,exp,mod_p,q):
		D=mod_p;F=C(D)-1;E=[1]+[0]*(F-1);A=[aX%q for aX in base]+[0]*(F-C(base));B=exp
		while B:
			if B&1:E=Ad(E,A,D,q)
			B>>=1
			if B:A=Ad(A,A,D,q)
		return E
	def Ae(c,q):
		for B in A(q):
			if n(c,B,q)==0:
				if G(A*c[A]*H(B,A-1,q)for A in A(1,C(c)))%q==0:return k
		return l
	AN,BL=B*B-1,B-1
	def AC(a,b,poly_lo,eps):A=BK([a%AB,-b%AB],eps,poly_lo,AB);A[0]=(A[0]-1)%AB;assert e(A%B==0 for A in A);return[aY//B%B for aY in A]
	Af=229;BM=[aa for aa in A(2,Af+1)if e(aa%A for A in A(2,w(aa**.5)+1))];y,z={},{}
	for V in BM:
		if Ae(m,V):
			for i in A(V):
				if n(m,i,V)==0:y[V,i]=C(y)
		if Ae(x,V):
			for i in A(V):
				if n(x,i,V)==0:z[V,i]=C(z)
	BN=A8({ap for(ap,aq)in y});BO=A8({ar for(ar,at)in z})
	def Ag(val,plist):
		A=abs(val)
		if A==0:return
		C={}
		for B in plist:
			if A%B==0:
				D=0
				while A%B==0:A//=B;D+=1
				C[B]=D
				if A==1:return C
		return C if A==1 else W
	def BP(fb_d,poly_lo):
		F=poly_lo;I=[]
		for(B,J)in fb_d:
			K=w(round(AM.log2(B)*4));I.append((B,J,K));D,E=B*B,J
			while D<=128:L=n(F,E,D);M=G(A*F[A]*H(E,A-1,B)for A in A(1,C(F)))%B;E=(E-L*H(M,-1,B))%D;I.append((D,E,K));D*=B
		return I
	BQ=BP(z,x);K=[]
	for D in A(1,500):
		AO,BR=-w(3.2*D)-10,w(4.2*D)+10;AP=BR-AO+1;A0=bytearray(AP)
		for(V,i,BS)in BQ:
			BT=(D*i-AO)%V
			for S in A(BT,AP,V):Ah=A0[S]+BS;A0[S]=Ah if Ah<255 else 255
		BU=0 if D<=20 else max(0,w(12.*AM.log2(D))-20)
		for S in A(AP):
			if A0[S]<BU:continue
			F=AO+S
			if F==0 or AM.gcd(abs(F),D)!=1 or(F-D*h)%O==0:continue
			Ai=Ag(F*F*F-2*F*F*D-6*F*D*D-D*D*D,BO)
			if Ai is W:continue
			Aj=Ag(F**4-7*F*F*D*D+F*D**3+3*D**4,BN)
			if Aj is W:continue
			A1={y[aI,F*H(D,-1,aI)%aI]:aJ for(aI,aJ)in Aj.items()};AD={z[aK,F*H(D,-1,aK)%aK]:aL for(aK,aL)in Ai.items()};K.append((F,D,A1,AD))
	o=C(y);Q=A9(A(C(K)))
	while l:
		p={}
		for P in Q:
			for M in K[P][2]:p[M]=p.get(M,0)+1
			for M in K[P][3]:p[o+M]=p.get(o+M,0)+1
		Ak=[au for au in Q if e(p[A]>=2 for A in K[au][2])and e(p[o+A]>=2 for A in K[au][3])]
		if C(Ak)==C(Q):break
		Q=Ak
	Al=A8({aN for aM in Q for aN in K[aM][2]}|{o+aP for aO in Q for aP in K[aO][3]});q={ac:ab for(ab,ac)in U(Al)};I=C(Al);E=I+R;AQ=[0]*I;AE=[]
	for P in Q:
		Am=[q[aQ]for aQ in K[P][2]]+[q[o+aR]for aR in K[P][3]];AE.append(Am)
		for T in Am:AQ[T]+=1
	AF=A9(A(C(Q)))
	while C(AF)>E:
		An=[av for av in AF if e(AQ[A]>=2 for A in AE[av])]
		if not An:break
		Ao=max(An,key=lambda i:(C(AE[i]),i));AF.remove(Ao)
		for T in AE[Ao]:AQ[T]-=1
	Ap=[Q[ad]for ad in AF];assert C(Ap)==E;A2=random.Random(2026+seed);AR=[[A2.randrange(B)for aw in A(4)]for ae in A(3)];BV=[[A2.randrange(B)for ax in A(3)]for af in A(2)];r=[[0]*I for ag in A(E)];A3=[[0]*E for ah in A(R)]
	for(N,P)in U(Ap):
		F,D,A1,AD=K[P]
		for(M,a)in A1.items():r[N][q[M]]=a%B
		for(M,a)in AD.items():r[N][q[o+M]]=-a%B
		AS=AC(F,D,m,AN);A0=AC(F,D,x,BL)
		for X in A(3):A3[X][N]=G(AR[X][A]*AS[A]for A in A(4))%B
		for X in A(2):A3[3+X][N]=G(BV[X][A]*A0[A]for A in A(3))%B
	Aq,_=v(r,I,B);AG=[r[ai]+[A3[aS][ai]for aS in A(R)]for ai in A(E)];AH,Ar=v(AG,E,B);BW=[r[aj]+[A3[aT][aj]for aT in A(R-1)]for aj in A(E)];_,As=v(BW,E-1,B)
	def At(x_live):return[G(r[aZ][A]*x_live[A]for A in A(I))%B for aZ in A(E)]
	BX=lambda x:At(x[:I]);Y=-(-E//J)+1+-(-E//(128*J));a=-(-E//L)+1+-(-E//(128*L));Au=[[A2.randrange(B)for ay in A(I)]for ak in A(J-R)];BY=[A3[az]for az in A(R)]+[At(aA)for aA in Au];s=[BY]
	for _ in A(Y+a-1):s.append([BX(aU)for aU in s[-1]])
	BZ=A9(s[0][R:])+[aC for aB in A(1,Y+1)for aC in s[aB]];Av=Ac(BZ,B)
	def Aw(onepad):
		D=A9(A(I));A2.shuffle(D);C=[[]for a0 in A(L)]
		for(G,H)in U(D):C[G%L].append(H)
		for F in A(R):C[0 if onepad else F%L].append(I+F)
		return[[A2.randrange(1,B)if a2 in set(a1)else 0 for a2 in A(E)]for a1 in C]
	Ax=Aw(k);Ba=Aw(l)
	def Ay(X):
		D=[];C=[a3[:]for a3 in X]
		for F in A(a):D.extend(C);C=[[G(a4[A]*AG[A][a6]for A in A(E))%B if a6<I else 0 for a6 in A(E)]for a4 in C]
		return Ac([[G(a5[A]*AG[A][a7]for A in A(E))%B for a7 in A(E)]for a5 in D],B)
	Az=Ay(Ax);A_=Ay(Ba);B0=[[[G(Ax[aD][A]*s[al][aV][A]for A in A(E))%B for aV in A(J)]for aD in A(L)]for al in A(Y+a)];AI=[]
	for N in A(a):
		for Bb in A(L):
			f=[]
			for Z in A(Y+1):f.extend(B0[N+Z][Bb])
			AI.append(f)
	Bc,_=v(AI,(Y+1)*J,B);t=Y+a;c=J+L;j=[[[1 if aE==am else 0]for aE in A(c)]for am in A(c)];b=[[[0]*t for aF in A(L)]for an in A(c)]
	for M in A(t):
		for N in A(L):
			for Z in A(J):b[Z][N][M]=B0[M][N][Z]%B
	for N in A(L):b[J+N][N][0]=-1%B
	AJ=[0]*J+[1]*L
	for M in A(t):
		for N in A(L):
			AT=[aW for aW in A(c)if b[aW][N][M]]
			if not AT:continue
			d=min(AT,key=lambda c:(AJ[c],c));Bd=H(b[d][N][M],-1,B)
			for T in AT:
				if T==d:continue
				B1=b[T][N][M]*Bd%B
				for P in A(L):
					B2,B3=b[T][P],b[d][P]
					for AK in A(M,t):
						if B3[AK]:B2[AK]=(B2[AK]-B1*B3[AK])%B
				for f in A(c):
					A4,AU=j[T][f],j[d][f]
					if C(A4)<C(AU):A4.extend([0]*(C(AU)-C(A4)))
					for(B4,B5)in U(AU):
						if B5:A4[B4]=(A4[B4]-B1*B5)%B
			for P in A(L):b[d][P]=[0]+b[d][P][:-1]
			for f in A(c):j[d][f]=[0]+j[d][f]
			AJ[d]+=1
	B6=e(b[B][C][D]==0 for B in A(c)for C in A(L)for D in A(t));B7=[ao for ao in A(c)if AJ[ao]<=Y];AV=0;AW=l;g=W;B8=W
	for T in B7:
		A5=AJ[T];Be=[j[T][aG][A5]if A5<C(j[T][aG])else 0 for aG in A(R)]
		if not any(Be):continue
		AV+=1;u=[0]*((Y+1)*J)
		for Z in A(A5+1):
			B9=A5-Z
			for BA in A(J):BB=j[T][BA];u[Z*J+BA]=BB[B9]if B9<C(BB)else 0
		if not e(G(D[A]*u[A]for A in A(C(u)))%B==0 for D in AI):AW=k
		if g is W:
			Bf=u[:R];A6=[0]*I
			for X in A(J-R):
				A7=u[R+X]
				if A7:
					Bg=Au[X]
					for S in A(I):A6[S]=(A6[S]+A7*Bg[S])%B
			for Z in A(1,Y+1):
				Bh=u[Z*J:(Z+1)*J];Bi=s[Z-1]
				for X in A(J):
					A7=Bh[X]
					if A7:
						Bj=Bi[X]
						for S in A(I):A6[S]=(A6[S]+A7*Bj[S])%B
			BC=A6+Bf
			if e(G(AG[C][A]*BC[A]for A in A(E))%B==0 for C in A(E)):g=BC;B8=A5
	BD=k;BE=[];BF={};BG=0
	if g is not W:
		BH=g[:I];Bk=next(A for A in A(2,O)if H(A,B,O)==1);Bl={H(Bk,aH,O):aH for aH in A(B)};Bm=H(2,-1,B);AL=set()
		for P in Q:
			F,D,A1,AD=K[P];Bn=H((F-D*h)%O,2,O);BI=Bl[Bn]*Bm%B
			if BI==0:continue
			BG+=1;AS=AC(F,D,m,AN);Bo=(G(B*BH[q[A]]for(A,B)in A1.items())+G(g[I+B]*G(AR[B][A]*AS[A]for A in A(4))for B in A(3)))%B;AL.add(Bo*H(BI,-1,B)%B)
		BE=A8(AL)
		def BJ(ri):C,D,E,H=K[ri];F=AC(C,D,m,AN);return(G(B*BH[q[A]]for(A,B)in E.items())+G(g[I+B]*G(AR[B][A]*F[A]for A in A(4))for B in A(3)))%B
		AX,AY=Q[0],Q[1];AZ=H((K[AX][0]-K[AX][1]*h)%O,2,O);Aa=H((K[AY][0]-K[AY][1]*h)%O,2,O);Ab=BJ(AY)*H(BJ(AX),-1,B)%B;BD=C(AL)==1 and 0 not in AL and H(AZ,Ab,O)==Aa;BF=AA(g=AZ,z=Aa,recovered_x=Ab,g_pow_x_eq_z=H(AZ,Ab,O)==Aa)
	Bp=Aq==I and AH==E-1 and C(Ar)==1 and C(As)==0 and Av==I and Az==AH and A_<AH and B6 and AW and AV>=1 and g is not W and BD;return AA(ok=Bp,planted=k,p=O,r=B,m_root=h,s_sm=R,B=Af,raw_rels=C(K),purged_rels=C(Q),held_out_rels=C(Q)-E,N=I,D=E,m_blk=L,n_blk=J,d=Y,e=a,rank_M=Aq,rank_A=AH,ker_dim_5sm=C(Ar),neg_ctrl_4sm_ker_dim=C(As),rank_Cd=Av,rank_OeA=Az,control_one_block=A_,hankel_shape=[C(AI),(Y+1)*J],hankel_rank=Bc,sigma_basis_order=t,sigma_basis_dim=c,sigma_basis_residual_zero=B6,sigma_basis_retained_cols=C(B7),sigma_basis_sm_nonzero_cols=AV,sigma_basis_chosen_deg=B8,sigma_basis_all_hf_zero=AW,hankel_recovered=g is not W,dlp_checked_rels=BG,dlp_all_purged_ratios=BE,dlp_live_check=BF)
def z6_schur_operator_check(log2U=12,seed=1):
	B,J=I,S;assert(J*J-J+1)%B==0 and(H(J,3,B)+1)%B==0;P=V.verify_live_filtsim_z6(log2U=log2U,seed=seed,W_target=160,check_kernel=k,export_hx=l);N=P['hx_export'];T=N['rec_order'];F={p:o for(o,p)in U(T)};W=N['live_col_set'];X={r:q for(q,r)in U(W)};Y,K=C(T),C(W);Z=N['live_rows_full'][:K+5];g=C(Z);a=[]
	for(b,(h,i,Q))in U(N['recorded_pivots']):
		assert F[h]==b;j=-H(O(i),-1,B)%B;L,M=[],[]
		for(D,R)in Q:
			E=O(R)
			if not E:continue
			if D in F:assert F[D]<b;L.append((F[D],E))
			else:M.append((X[D],E))
		a.append((j,L,M))
	c=[]
	for Q in Z:
		L,M=[],[]
		for(D,R)in Q:
			E=O(R)
			if not E:continue
			if D in F:L.append((F[D],E))
			else:M.append((X[D],E))
		c.append((L,M))
	def d(x_live):
		C=x_live;A=[0]*Y
		for(D,(E,F,H))in U(a):I=G(C*A[B]for(B,C)in F)+G(B*C[A]for(A,B)in H);A[D]=E*(I%B)%B
		return[(G(C*A[B]for(B,C)in y)+G(B*C[A]for(A,B)in z))%B for(y,z)in c]
	m=[d([1 if x==s else 0 for x in A(K)])for s in A(K)];n=random.Random(seed*7919+31);e=[n.randrange(B)for t in A(K)];f=d(e)==[G(m[A][u]*e[A]for A in A(K))%B for u in A(g)];return AA(ok=f,Q=B,chi=J,filt_N_sw=P['N_sw'],filt_N_l=P['N_l'],filt_N_rec=Y,l_rec_strict_lower_tri=l,prod_eq_schur=f)
def run_hx_1d_kernel_z6(log2U=12,seed=1,m=10,n=8,s_sm=5):A=unplanted_two_sided_5sm_bw_check(seed=seed,m_blk=m,n_blk=n);B=z6_schur_operator_check(log2U=log2U,seed=seed);C=bool(A['ok']and B['ok']);return AA(checks_ok=C,ok=C,source='unplanted_2sided_5sm_bw_F1187_dlp_plus_z6_schur',unplanted_5sm_bw=A,z6_schur_op=B)
if __name__=='__main__':
	E=run_hx_1d_kernel_z6();print(json.dumps(E,separators=(',',':')))
	if not E['checks_ok']:raise SystemExit(1)
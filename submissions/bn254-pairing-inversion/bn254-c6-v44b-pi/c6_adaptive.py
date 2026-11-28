BE='MONT5_ROUND'
BM=reversed
w='ret'
v=':'
r=','
AE=all
AD=next
A5=any
y=abs
u=sorted
k=float
n=list
l=zip
h=dict
X=enumerate
W=set
U=min
T=True
R=False
O=tuple
P=1.
N=None
L=max
J=round
I=.0
H=int
E=pow
B=sum
C=len
A=range
import base64,json,lzma,math as G,stor_bound as CJ
from pathlib import Path
from fractions import Fraction as F
def verify_threshold_counterexample():
	L,F,C=6e1,P,.1;H=.5*L/F;B=lambda t:(P+H*(P-t))**-F
	for M in A(81):G=M/1e2;assert .1*G+.9*G**100<=B(G)+1e-12
	assert .1>B(C)+.05;D=B(I)+C
	for E in A(2,65):K=U(.8,(E-1)*(P+P/H)/(E-1+F));D=U(D,U(P,B(K)/K**(E-1))+C**E)
	assert D>=.1-1e-12;return h(pgf_psi_invalid=J(B(C),5),true_fail=.1,tr_bound=J(D,5),tr_valid=T)
def verify_c6_adaptive_h5_live():
	import model_bn254_v6 as D,verify_bn254_v6 as w,descent_bn254_v6 as B2;AO=verify_threshold_counterexample();K=json.loads(lzma.decompress(base64.b85decode((Path(__file__).resolve().parent/'c6_data.b85').read_bytes())).decode())['adaptive_h5'];B3,x,B4,g=K['names'],K['cal_tbl'],K['after_by_set'],K['first_stats']
	def z(si,cs,cb):
		D=f"{si}:{cs}:{cb}"
		if D in x:A,B,C=x[D];return B/A,C/A
		A=B=C=0
		for(E,(F,G,I))in x.items():
			L,J,K=map(H,E.split(v))
			if J==cs and K==cb:A+=F;B+=G;C+=I
		return B/A,C/A
	AP={};AQ=AR=AS=AT=A0=AU=AV=AW=0
	for(V,B5)in X(B3):
		i=B4[V];AX=C(i);AY=B(A[1]for A in i);AZ=B(A[3]for A in i);Aa=B(A[4]for A in i);A1=A2=A3=A4=I
		for(Y,B6,B7,Y,Y,B8,Ab)in i:Ac=[B8]+([Ab]if Ab else[]);A1+=B6;A2+=B(z(V,A>>5,A&31)[0]for A in Ac);A3+=B7;A4+=B(z(V,A>>5,A&31)[1]for A in Ac)
		Ad,Ae=A1/A2,A3/A4;assert Ad<.9 and .9<Ae<1.1;AP[B5]=h(roots=g[V][2],first_n=g[V][0],first_fail=g[V][1],after_n=AX,after_fail=AY,oe_fail=J(Ad,3),oe_work=J(Ae,3),sh_prev=AZ,sh_run=Aa);A0+=AX;AU+=AY;AQ+=A1;AR+=A2;AS+=A3;AT+=A4;AV+=AZ;AW+=Aa
	A6,Af=AQ/AR,AS/AT;assert A0==1856 and A6<.85 and .95<Af<1.05;b=K['conc_rows'];assert C(b)==10 and L(A[5]for A in b)<15. and L(A[6]for A in b)<3.;c,M,j,Ag=H(K['u_t']),H(K['p_t']),H(K['lx_b']),H(K['ns_max']);Ah=1<<j;d,e=D.ZK_NUM,D.W_COEFFS;Ai=[-4-aa if ab==0 else-aa for(ab,aa)in X(e)];B9=D._ppow_h([0,1,0,0,0,0],M,M);Aj=D._ppow_h(B9,M,M);Ak=D._ppow_h(Aj,M,M);BA=D._ppow_h(Ak,M**3,M);BB=[c*ac%M for ac in e];Z=BB;A7=[1,0,0,0,0,0]
	for Y in A(6):A7=D._pmul_h(A7,Z,M);Z=D._ppow_h(Z,M,M)
	A8=M==1+6*c and BA==[0,1,0,0,0,0]and Aj!=[0,1,0,0,0,0]and Ak!=[0,1,0,0,0,0]and E(A7[0],(M-1)//2,M)==M-1;assert A8;Al=[F(ad,3)for ad in D.SIGMA_NUM];A9=[[F(1)]+[F(0)]*5,Al]
	for Y in A(4):A9.append(D._pmul_h_q(A9[-1],Al))
	Am=O(O(H(3*A9[A][B])for A in A(6))for B in A(6));BC=[[1 if ae==ao else 0 for ao in A(6)]for ae in A(6)];An=D._isp;Ao=G.sumprod
	def f(n,p):
		n%=p
		if n==0:return 0
		if E(n,(p-1)//2,p)!=1:return-1
		if p%4==3:A=E(n,(p+1)//4,p);return U(A,p-A)
		B,I=p-1,0
		while B%2==0:B//=2;I+=1
		G=2
		while E(G,(p-1)//2,p)!=p-1:G+=1
		J,K,C,A=I,E(G,B,p),E(n,B,p),E(n,(B+1)//2,p)
		while C!=1:
			H,F=1,C*C%p
			while F!=1:F=F*F%p;H+=1
			D=E(K,1<<J-H-1,p);J,K,C,A=H,D*D%p,C*D*D%p,A*D%p
		return U(A,p-A)
	P=D._mh_int
	def m(p):C=O(Ao(A,p)for A in Am);E=O(Ao(A,C)for A in Am);D=P(P(p,C),E);A=-(D[2]//729);B=D[0]//729-2*A;return B,A,B*B+2*A*A
	AA={}
	def S(p):
		A=AA.get(p)
		if A is not N:return A
		if p<=3 or p&7 not in(1,3)or D._ppow_h([0,1,0,0,0,0],p,p)!=[0,1,0,0,0,0]:AA[p]=N;return
		B=f(-2,p);C=p-B;F=D._roots_mod(D.H_COEFFS,p);G=AD(A for A in F if(2+w._eval_poly_mod(e,A,p))%p==B);H=AD(A for A in F if(2+w._eval_poly_mod(e,A,p))%p==C);A=B,C,(G,H),(f(2+B,p),f(2+C,p)),(f(E(c*(B-2)%p,-1,p),p),f(E(c*(C-2)%p,-1,p),p));AA[p]=A;return A
	Ap={(H(af),H(ag))for(af,ag)in K['ret_f_split']};Aq={(H(ah),H(ai))for(ah,ai)in K['ret_g_split']};BD={(H(aj),H(ak))for(aj,ak)in K['ret_f_ns']};BE={(H(al),H(am))for(al,am)in K['ret_g_ns']};k=u({ap for(ap,aq)in Ap}|{ar for(ar,at)in Aq});Ar=D._primes_upto(Ah);AB=W(k);AC=[];As={}
	for Q in Ar:
		if Q<=Ag and Q not in AB:
			if S(Q)is N:AC.append(Q);As[Q]=f(-2,Q)if Q>3 and Q&7 in(1,3)else-1
			else:k.append(Q);AB.add(Q)
	for Q in k:S(Q)
	def AE(xs):
		B=xs;B=n(B)
		while C(B)>1:B=[B[av]*B[av+1]if av+1<C(B)else B[av]for av in A(0,C(B),2)]
		return B[0]
	BF,BG,BH=AE(AC),AE(k),AE([au for au in Ar if au>Ag and au not in AB]);AF=[];Z=[0]*6
	def At(idx,rem):
		B=idx
		if B==6:
			if A5(Z):AF.append(O(Z))
			return
		D=G.isqrt(rem)
		for C in A(-D,D+1):Z[B]=C;At(B+1,rem-C*C)
	At(0,H(K['ball_r2']));assert C(AF)==H(K['ball_pts'])==484;BI=[O((B,A)for(B,A)in X(an)if A)for an in AF]
	def Au(n,c0,c1,sd):
		B=n;F=BD if sd==0 else BE;H=Ap if sd==0 else Aq;I=G.gcd(B,BF)
		if I>1:
			for A in AC:
				if I%A==0:
					while B%A==0:B//=A
					C=As[A]
					if C>=0:
						D=(c0+c1*C)%A==0;E=(c0+c1*(A-C))%A==0
						if not D and not E or D and(A,0)not in F or E and(A,1)not in F:return
					elif(A,0)not in F:return
		J=G.gcd(B,BG)
		if J>1:
			for A in k:
				if J%A==0:
					while B%A==0:B//=A
					C,K=S(A)[0],S(A)[1];D=(c0+c1*C)%A==0;E=(c0+c1*K)%A==0
					if not D and not E or D and(A,0)not in H or E and(A,1)not in H:return
		if G.gcd(n,BH)>1:return
		return B
	def Av(cof,nw,km):
		K=cof
		if K==1:return[]
		if km<=0 or K.bit_length()>nw*km:return
		O=1<<nw;I=[];L=[K]
		while L:
			B=L.pop()
			if B==1:continue
			if An(B):
				if B>O or B<=Ah or S(B)is N:return
				if B not in I:
					if C(I)>=km:return
					I.append(B)
			else:
				J=N
				for H in A(1,32):
					E=D=2;F=1
					while F==1:
						M=1;P,Q=E,D
						for R in A(16):E=(E*E+H)%B;D=(D*D+H)%B;D=(D*D+H)%B;M=M*(E-D)%B
						F=G.gcd(M,B)
						if F==B:
							E,D,F=P,Q,1
							while F==1:E=(E*E+H)%B;D=(D*D+H)%B;D=(D*D+H)%B;F=G.gcd(y(E-D),B)
					if 1<F<B:J=F;break
				if J is N:return
				L.extend((J,B//J))
		return I
	o={}
	def AG(sd,q,oi,nw,km):
		J=sd;L=J,q,oi,nw,km
		if L in o:return o[L]
		Ad,Ae,AH,AI,AJ=S(q);AK=AH[oi];AL=AI[oi]if J==0 else AJ[oi];s=[w._eval_zk_mod([1 if aA==aw else 0 for aA in A(6)],AK,q)for aw in A(6)];A6=E(s[0],-1,q);I=[[0]*12 for ax in A(12)];I[0][0]=q
		for H in A(1,6):I[H][H]=1;I[H][0]=-s[H]*A6%q
		for H in A(6):I[6+H][6+H]=1;I[6+H][0]=AL*s[H]*A6%q
		B2.bareiss_lll_12d(I,BC,1,1,R);A7=O(B(I[0][A]*d[A][C]for A in A(6))for C in A(6));A8=O(B(I[0][6+A]*d[A][C]for A in A(6))for C in A(6));t=A9=N
		for AA in A(1,12):
			AB=O(B(I[AA][A]*d[A][C]for A in A(6))for C in A(6));AC=O(B(I[AA][6+A]*d[A][C]for A in A(6))for C in A(6))
			if P(A7,AC)!=P(A8,AB):t,A9=AB,AC;break
		if t is N:o[L]=[];return[]
		AM=[O(A//3 for A in P(d[ay],t))for ay in A(6)];AN=[O(A//3 for A in P(d[az],A9))for az in A(6)];u=[]
		for AO in BI:
			M,Q,U,V,W,X=A7;Y,Z,a,b,f,g=A8
			for(H,G)in AO:
				D,F=AM[H],AN[H]
				if G==1:M+=D[0];Q+=D[1];U+=D[2];V+=D[3];W+=D[4];X+=D[5];Y+=F[0];Z+=F[1];a+=F[2];b+=F[3];f+=F[4];g+=F[5]
				elif G==-1:M-=D[0];Q-=D[1];U-=D[2];V-=D[3];W-=D[4];X-=D[5];Y-=F[0];Z-=F[1];a-=F[2];b-=F[3];f-=F[4];g-=F[5]
				else:M+=G*D[0];Q+=G*D[1];U+=G*D[2];V+=G*D[3];W+=G*D[4];X+=G*D[5];Y+=G*F[0];Z+=G*F[1];a+=G*F[2];b+=G*F[3];f+=G*F[4];g+=G*F[5]
			v=M,Q,U,V,W,X;x=Y,Z,a,b,f,g
			if not A5(v)or not A5(x):continue
			h=O(A//3 for A in P(v,v));i=O(A//3 for A in P(x,x))
			if J==0:y=P(Ai,i);j,k,l=m(O(h[A]+y[A]for A in A(6)))
			else:z=P(e,h);j,k,l=m(O(c*z[A]-i[A]for A in A(6)))
			if l%q!=0 or l==0:continue
			AD=Au(l//q,j,k,J);n=Av(AD,nw,km)if AD is not N else N
			if n is N:continue
			if J==0:z=P(e,h);p,r,A0=m(O(c*z[A]-i[A]for A in A(6)))
			else:y=P(Ai,i);p,r,A0=m(O(h[A]+y[A]for A in A(6)))
			AE=Au(A0,p,r,1-J)if A0 else N;A1=Av(AE,nw,km-C(n))if AE is not N else N
			if A1 is N:continue
			(AP,AQ,AR),(AS,AT,AU)=((n,j,k),(A1,p,r))if J==0 else((A1,p,r),(n,j,k));A2=[];K=T
			for(AF,AV,AW,AX)in((0,AP,AQ,AR),(1,AS,AT,AU)):
				for A3 in AV:
					AY,AZ,_,Aa,Ab=S(A3);AG=0
					for(A4,Ac)in((0,AY),(1,AZ)):
						if(AW+AX*Ac)%A3==0:
							if(Aa[A4]if AF==0 else Ab[A4])<=0 or C(A2)>=km:K=R;break
							A2.append((AF,A3,A4));AG+=1
					if not K or AG==0:K=R;break
				if not K:break
			if K:u.append(A2);break
		o[L]=u;return u
	def Aw(sd,q,oi,nw,km,dsh,rec=R):
		G=AG(sd,q,oi,nw,km);A=1
		if not G:return R,A,([],0,0)
		B=G[0]
		if not B:return T,A,([[0,0,0,0,0]]if rec else[],0,0)
		D=T;E=0;F=0;I=0
		for(N,J,O)in B:
			F+=1;P=J.bit_length();K=L(j,P-dsh);Q=0 if K==j else 1;M,S,(X,U,V)=Aw(N,J,O,K,Q,dsh,R);E+=S;F+=U;I+=(0 if M else 1)+V
			if not M:D=R;break
		A+=E;W=[[0,H(not D),E,C(B),0]]if rec else[];return D,A,(W,F,I)
	AH=0
	for(V,AI,AJ,AK,BJ,BK,BL,Ax)in K['live_witnesses']:Ay,Y,(p,Y,Y)=Aw(AI,AJ,AK,BJ,BK,BL,T);assert p and p[0][0]==0 and p[0][1]==Ax[0][1]and p[0][3]==Ax[0][3];AH+=1
	a=[];l=(1<<16)+1
	while C(a)<64:
		if An(l):
			AL=S(l)
			if AL is not N:
				for q in(0,1):
					if AL[3][q]>0:a.append((0,l,q))
					if AL[4][q]>0:a.append((1,l,q))
		l+=2
	a=a[:64];r=s=AM=Az=t=0;AN=I
	def A_(sd,q,oi,nw,km,dsh):
		M=AG(sd,q,oi,nw,km);A=1
		if not M:return R,R,A,0,0,0,I
		N=M[0]
		if not N:return T,T,A,0,0,0,I
		B=T;C=0;D=0;E=0;F=0;G=I
		for(H,J,O)in N:
			D+=1;P=J.bit_length();K=L(j,P-dsh);Q=0 if K==j else 1;G+=z(0,H,P)[0];U=AG(H,J,O,K,Q)
			if not U:C+=1;F+=1;E+=1;B=R;break
			S,a,V,W,X,Y,Z=A_(H,J,O,K,Q,dsh);C+=V;D+=W;E+=(0 if S else 1)+X;F+=Y;G+=Z
			if not S:B=R;break
		A+=C;return B,T,A,D,E,F,G
	for(AI,AJ,AK)in a:
		Ay,BM,BR,BN,BO,BP,BQ=A_(AI,AJ,AK,16,2,1)
		if Ay:s+=1
		if BM:r+=1
		AM+=BN;Az+=BO;t+=BP;AN+=BQ
	B0=t/L(AN,1e-12);B1=r==s+t and r>=25 and s>=18 and AM>=15 and B0<=.85;return h(ok=AH==3 and B1 and A8 and AO['tr_valid'],strict_k1_zero_retry_protocol=T,live_k1_zero_retry_ok=B1,k1_17b_roots=C(a),k1_17b_nonempty=r,k1_17b_completed=s,k1_17b_reached_children=AM,k1_17b_slot_zeros=t,k1_17b_exp_p0_sum=J(AN,3),k1_17b_direct_slot_oe_p0=J(B0,4),k1_17b_child_fails=Az,threshold_counterexample=AO,tower_u2_p13_verified=A8,roots_total=B(A[2]for A in g),first_attempts_total=B(A[0]for A in g),after_fail_total=A0,after_fail_failures=AU,pooled_oe_fail=J(A6,4),boole_union_attempt_oe=J(A6,4),pooled_oe_work=J(Af,4),shared_prev_total=AV,shared_run_total=AW,conc_bins_checked=C(b),max_C_emp=L(A[5]for A in b),max_C_bday=L(A[6]for A in b),live_witnesses_replayed=AH,per_setting=AP)
def verify_h3_scaling_and_nu_h():
	import model_bn254_v6 as D;O=D.load_data()['filtsim_raw'];Q='m184_z6_62.50_62.00_e0.687_W60_';K=[aa for aa in(20,21,22,23,24,25,26)if f"{Q}{aa}"in O];S=[];g=[];i=[];A9=[];u=T;v=T;j=10**9
	for Aq in K:
		Z=O[f"{Q}{Aq}"];m=C(Z);S.append(B(A[1]for A in Z)/m);g.append(B(A[2]for A in Z)/m);n=[]
		for H in Z:
			if C(H)>=12:
				AA,AB,AC,AD,AF,AG=H[6],H[7],H[8],H[9],H[10],H[11]
				if y(AA+(AB+AC+AD+AF)/AG-P)>=.0005:u=R
				n.append(AA*D.C_UNIT_5+(AB*2.*D.C_UNIT_5+AC*3.*D.C_UNIT_5+AD*4.*D.C_UNIT_5+AF*D.C_GEN_1)/AG)
				if C(H)>=14 and(H[12]!=0 or H[13]!=0):v=R
				if C(H)>=15:j=U(j,H[14])
			elif C(H)==8:n.append(H[7])
			else:n.append(H[6]*D.C_UNIT_5+(P-H[6])*D.C_GEN_1)
		i.append(B(n)/m);A9.append(B(A[4]for A in Z)/m)
	w=AE(B[D]>B[D+1]for B in(S,g)for D in A(C(K)-1));o=[G.log(ab)for ab in K];x=[G.log(ac)for ac in S];z,AH=B(o)/C(o),B(x)/C(x);a=-B((A-z)*(B-AH)for(A,B)in l(o,x))/B((A-z)**2 for A in o);A0=G.exp(AH+a*z);AI=L(y(A0/B**a-A)/A for(B,A)in l(K,S));AJ=A0/D.log2_U_z6(D.OP['lf'],D.OP['lg'])**a;A1=D.log2_U_z6(D.OP['lf'],D.OP['lg']);V=C(K);M=[B(A[14]for A in O[f"{Q}{ad}"])/C(O[f"{Q}{ad}"])/2.**ad for ad in K];Ar=[B(B[1]*2.**ae-B[14]for B in O[f"{Q}{ae}"])/C(O[f"{Q}{ae}"])/G.floor(D.TEXC_REL*2.**ae)for ae in K];A2=D.NC27;A3=[G.log(ai)for ai in M]if U(M)>0 else[I]*V;W=B(K)/V;p=B(A3)/V;AK=G.exp(p+B((A-W)*(B-p)for(A,B)in l(K,A3))/B((A-W)**2 for A in K)*(A1-W));A4=[G.log(af)for af in K];A6=B(A4)/V;AL=-B((A-A6)*(B-p)for(A,B)in l(A4,A3))/B((A-A6)**2 for A in A4);AM=G.exp(p+AL*A6)/A1**AL;AN=[P/ag for ag in M];AO=B(AN)/V;As=B((A-W)*(B-AO)for(A,B)in l(K,AN))/B((A-W)**2 for A in K);AP=AO+As*(A1-W);AQ=P/AP if AP>0 else P;A7=h(ncols_over_U=[J(aj,8)for aj in M],surplus_over_texc=[J(ak,4)for ak in Ar],max_ratio=J(L(M[A+1]/M[A]for A in A(V-1)),4),pred_ypow=AM,pred_exp=AK,pred_invlin=AQ,cap=A2,cap_over_max_pred=J(A2/L(AM,AK,AQ),3),cap_over_c26=J(A2/M[-1],3));AR=U(M)>0 and A7['max_ratio']<P and M[-2]<=D.NC27 and A7['cap_over_max_pred']>=P;AS=D._primes_upto(16384);q=b=c=0;AT=AU=AV=0;d=r=AW=AX=AY=I;At=[D.F_REC_COEFFS[ah]for ah in(0,2,4,6,8)]
	for F in AS:
		if F in(2,3,19):continue
		c+=1;N=G.log2(F);d+=N;s=D._ppow_h([0,1,0,0,0,0],F,F);s[1]=(s[1]-1)%F
		if not A5(s):b+=1;r+=N;A8=D.H_COEFFS
		else:
			A8=D._pgcd(s,D.H_COEFFS,F);t=C(A8)-1
			if t==3:q+=1;r+=N
			elif t==6:b+=1;r+=N
			elif t!=0:raise AssertionError(f"d={t} p={F}")
			else:continue
		e=f=R
		for H in D._roots_mod(A8,F):
			Au=B(D.W_COEFFS[A]*E(H,A,F)for A in A(6))%F;AZ=-((D.U_BN254+1)%F*Au+2)%F
			if not f and(AZ==0 or E(AZ,(F-1)//2,F)==1):f=T
			if not e:
				X=[B(al[A]*E(H,A,F)for A in A(6))%F for al in At];Aa=[X[0],0,X[1],0,X[2],0,X[3],0,X[4]];Y=D._ppowmod([0,1],F,Aa,F)
				if C(Y)<2:Y+=[0]*(2-C(Y))
				Y[1]=(Y[1]-1)%F
				if C(D._pgcd(Y,Aa,F))>1 or X[0]==0:e=T
			if e and f:break
		if e:AT+=1;AW+=N
		if f:AU+=1;AX+=N
		if e or f:AV+=1;AY+=N
	Ab=3,13,29,983,11003,405928799;Ac=1
	for Av in Ab:Ac*=Av
	Ad=Ac==D.U_BN254 and AE(D._isp(A)and A<1<<30 for A in Ab);Ae=q+2*b;Af=Ae/(C(AS)/3.);Ag=(q+b)/k(c);Ah=r/d;Ai=AT/k(c);Aj=AW/d;Ak=AU/k(c);Al=AX/d;Am=AV/k(c);An=AY/d;Ao=Ai<=D.KH_F_CAP and Aj<=D.KH_F_CAP and Ak<=D.KH_G_CAP and Al<=D.KH_G_CAP and Am<=D.KH_FG_CAP and An<=D.KH_FG_CAP;Ap=1.15;Aw=K==[20,21,22,23,24,25,26]and j>=150 and w and u and v and a>2.5 and AJ<D.NU_SW_CAP/1e1 and D.NU_SW_CAP/S[-1]>=Ap and D.NU_L_CAP/g[-1]>=Ap and L(i)<D.K5_CAP and L(A9)<D.TAB_R_CAP and .97<=Af<=1.03 and Ag<=D.KH_SPLIT_CAP and Ah<=D.KH_SPLIT_CAP and Ao and Ad and AR;return h(ok=Aw,fixed_instance_deterministic_h3=T,zero_loss_filtsim_ok=v,min_ncols_all28=j,scales=K,dim_monotone_20_26=w,dim_and_ncols_monotone_20_26=w and A7['max_ratio']<P,hist_v4_exact_ok=u,inv_log_A=J(A0,3),inv_log_beta=J(a,3),max_rel_fit_err=J(AI,4),pred_nsw_prod=J(AJ,6),nu_sw_cap=D.NU_SW_CAP,cap_over_top_scale=J(D.NU_SW_CAP/S[-1],3),nu_l_cap_over_top_scale=J(D.NU_L_CAP/g[-1],3),k5_v8_top_scale=J(i[-1],4),k5_v8_max_scale=J(L(i),4),nu_h_upto_2p14=J(Af,4),kh_split_frac_2p14=J(Ag,4),kh_log_wt_frac_2p14=J(Ah,4),f_rec_split_2p14=J(Ai,4),f_rec_log_2p14=J(Aj,4),g_rec_split_2p14=J(Ak,4),g_rec_log_2p14=J(Al,4),fg_split_2p14=J(Am,4),fg_log_2p14=J(An,4),fg_dens_ok=Ao,u_ram_all_lt_2p30=Ad,c6_d3_primes_2p14=q,c6_d6_primes_2p14=b,c6_deg1_orbits_2p14=Ae,h3_scale=A7,h3_scale_ok=AR)
def verify_split_crt_norms_and_mp_rt(seed=20261120):
	import random;from fractions import Fraction as d;import model_bn254_v6 as F;CA=143,495,663,1701,3023,3551,4485,4557,5055;Aw=[(1<<62)-aa for aa in CA];Ac=[]
	for D in Aw:
		Ad=F._roots_mod(F.H_COEFFS,D)
		if C(Ad)!=6:return h(ok=R,err='not_split')
		CB=E(3,-1,D);A1=[[B(F.ZK_NUM[aE][A]*E(ah,A,D)for A in A(6))*CB%D for aE in A(6)]for ah in Ad];AK=[((F.U_BN254+1)*B(F.W_COEFFS[A]*E(ai,A,D)for A in A(6))+2)%D for ai in Ad];Ae=[[B(F.F_REC_COEFFS[aF][A]*E(aj,A,D)for A in A(6))%D for aF in(0,2,4,6,8)]for aj in Ad];Ac.append((D,A1,AK,Ae))
	Ba=[[d(ak,3)for ak in ab]for ab in F.ZK_NUM]
	def e(a,b):
		B=[d(0)]*11
		for E in A(6):
			for C in A(6):B[E+C]+=a[E]*b[C]
		for D in A(10,5,-1):
			G=B[D]
			if G:
				for C in A(6):B[D-6+C]-=G*F.H_COEFFS[C]
				B[D]=d(0)
		return B[:6]
	def Bb(p_vec):
		F=[p_vec]
		for K in A(5):F.append(e(F[-1],[d(0),d(1),d(0),d(0),d(0),d(0)]))
		C=[[F[aM][aL]for aM in A(6)]for aL in A(6)];D=d(1)
		for B in A(6):
			E=AD((A for A in A(B,6)if C[A][B]!=0),N)
			if E is N:return 0
			if E!=B:C[B],C[E]=C[E],C[B];D=-D
			D*=C[B][B]
			for G in A(B+1,6):
				J=C[G][B]/C[B][B]
				for I in A(B,6):C[G][I]-=J*C[B][I]
		assert D.denominator==1;return y(H(D.numerator))
	CC=[d(ac)for ac in F.W_COEFFS];CD={ad:[d(al)for al in F.F_REC_COEFFS[ad]]for ad in(0,2,4,6,8)};Z=random.Random(seed);Ax=0;Ay=Az=0
	for _ in A(20):
		Af=[Z.randint(-400,400)for am in A(6)];Ag=[Z.randint(-500,500)for an in A(6)];AL=[];Bc=[]
		for(D,A1,AK,Ae)in Ac:
			A2=A_=1
			for V in A(6):AM=B(Af[A]*A1[V][A]for A in A(6))%D;q=B(Ag[A]*A1[V][A]for A in A(6))%D;o=AM*AM%D;r=q*q%D;A2=A2*(AK[V]*o+r)%D;s=o*o%D;t=r*r%D;AN=o*r%D;AO=Ae[V];A_=A_*(AO[0]*t*t+AO[1]*AN*t+AO[2]*s*t+AO[3]*s*AN+AO[4]*s*s)%D
			AL.append(A2);Bc.append(A_)
		def B0(rs,cnt):
			A=0;B=1
			for((C,*G),D)in l(Ac[:cnt],rs[:cnt]):F=(D-A)%C*E(B%C,-1,C)%C;A+=F*B;B*=C
			return y(A-B if A>B//2 else A)
		CE=B0(AL,9);CF=B0(Bc,9);Bd=[B(Af[A]*Ba[A][ao]for A in A(6))for ao in A(6)];Be=[B(Ag[A]*Ba[A][ap]for A in A(6))for ap in A(6)];o=e(Bd,Bd);r=e(Be,Be);CG=[(F.U_BN254+1)*aq+2*ar+at for(aq,ar,at)in l(e(CC,o),o,r)];s=e(o,o);t=e(r,r);AN=e(o,r);CH=[e(t,t),e(AN,t),e(s,t),e(s,AN),e(s,s)];B1=[d(0)]*6
		for(Ah,CI)in X((0,2,4,6,8)):CJ=e(CD[CI],CH[Ah]);B1=[B1[aG]+CJ[aG]for aG in A(6)]
		Bf=Bb(CG);Bg=Bb(B1)
		if CE!=Bf or CF!=Bg:Ax+=1
		Ay=L(Ay,Bf.bit_length());Az=L(Az,Bg.bit_length())
	a=1<<64;Ai=0;B2=0
	def B3(R_int,P_int,m,N):
		G=[R_int>>64*aN&a-1 for aN in A(N)];H=[P_int>>64*aO&a-1 for aO in A(m)];C=[0]*N
		for D in A(N):
			E=G[D]
			if E:
				for F in A(m):C[(D+F)%N]+=E*H[F]
		B=0
		for I in A(N-1,m-1,-1):B=(B<<64)+C[I]
		return B&(1<<64*(N-m))-1
	for A3 in(2,4,8):
		CK=8;Aj=2*A3+CK;Bh=a**(Aj-1)+(a**A3-1)//(a+1);Bi=a+1;CL=Bh*Bi//a**A3%a**(Aj-A3);CM=B3(Bh,Bi,A3,Aj)
		if not 0<=(CL-CM)%a**(Aj-A3)<A3*a:Ai+=1
	AP=4;AE=2;Ak=AP+4;CN=1<<64*(AE+Ak)-1
	for Bj in(0,1):
		AF=[a**AE-1-Z.randrange(1<<20)if Bj else Z.randrange(a,a**AE)for au in A(1<<AP)];Bk=a**40-1-Z.randrange(1<<30)if Bj else Z.randrange(1,a**40);AQ=[AF[:]]
		for _ in A(AP):B4=AQ[-1];AQ.append([B4[2*aJ]*B4[2*aJ+1]for aJ in A(C(B4)//2)])
		Bl=AQ[-1][0];CO=AE*(1<<AP);B5=[Bk%Bl*a**(CO+Ak)//Bl]
		for V in A(AP,0,-1):
			B6=AE*(1<<V-1);Bm=2*B6+Ak;A4=[]
			for(Ah,Al)in X(B5):A4.append(B3(Al,AQ[V-1][2*Ah+1],B6,Bm));A4.append(B3(Al,AQ[V-1][2*Ah],B6,Bm))
			B5=A4
		for(Al,B7)in l(B5,AF):
			if(Al*B7+CN>>64*(AE+Ak))%B7==Bk%B7:B2+=1
			else:Ai+=1
	Bn=J(B(G.log2(A)for A in Aw),3);A6=F.H_COEFFS;AG=F.W_COEFFS;CP=F.U_BN254;CQ=[F.F_REC_COEFFS[ae]for ae in(0,2,4,6,8)];B8=B9=I
	for D in F._primes_upto(500):
		if D<5 or D==19:continue
		CR=[av for av in A(D)if((((((A6[6]*av+A6[5])*av+A6[4])*av+A6[3])*av+A6[2])*av+A6[1])*av+A6[0])%D==0];Am=BA=0
		for p in CR:
			AH=[(((((aH[5]*p+aH[4])*p+aH[3])*p+aH[2])*p+aH[1])*p+aH[0])%D for aH in CQ]
			for x in A(D):
				if((((AH[4]*x+AH[3])*x+AH[2])*x+AH[1])*x+AH[0])%D==0:
					if x==0:Am+=1
					elif E(x,(D-1)//2,D)==1:Am+=2
			if AH[4]==0:Am+=1
			Bo=-((CP+1)*(((((AG[5]*p+AG[4])*p+AG[3])*p+AG[2])*p+AG[1])*p+AG[0])+2)%D
			if Bo==0:BA+=1
			elif E(Bo,(D-1)//2,D)==1:BA+=2
		B8+=(P/(D-1)-Am*D/k(D*D-1))*G.log2(D);B9+=(P/(D-1)-BA*D/k(D*D-1))*G.log2(D)
	CS=[0 if af%3 else af+1 for af in A(64)];An=0;BB=0
	for(z,Bp)in X(CS):
		CT=1 if Bp>=1 else 0;An=An<<1|CT
		if Bp>0:BB|=1<<63-z
	BC=0;Ao=0;c=47;BD=(1<<c)-1;BE=[(Z.randrange(1,500),Z.randint(12,25))for ag in A(8)];BF=B(A for(B,A)in BE)
	for g in(8,16,32,64,128,256,512,1024):
		K=2003;Bq=Z.randrange(1,K);Ap=(K+g-1)//g;CU=g-1<<c;BG=[0]*BF;Br=[0]*Ap;u=0
		for(BH,BI)in BE:
			f=BH%K
			for _ in A(BI):
				BG[u]=f;Br[f//g]+=1;u+=1;f+=Bq
				if f>=K:f-=K
		AR=[0]*(Ap+1)
		for AI in A(Ap):AR[AI+1]=AR[AI]+Br[AI]
		Bs=AR[:Ap];Bt=[0]*BF;AS=0;CV=g-1<<c
		for z in A(BF):
			f=BG[z];AI=f//g;AS=AS&~CV|(f&g-1)<<c;Bt[Bs[AI]]=AS;Bs[AI]+=1;AS+=1
			if AS&BD!=z+1:Ao+=1
		BJ=[]
		for(BH,BI)in BE:
			for x in A(BI):BJ.append((BH+x*Bq)%K)
		if BG!=BJ:Ao+=1
		for A7 in(20,390):
			for M in(0,g-3,K-15,Z.randrange(K)):
				for BK in(1,-1):
					Bu=W()
					for(z,f)in X(BJ):
						V=(f-M)*BK%K
						if V<A7:Bu.add((V,z))
					v=W();A8=[]
					if BK==1:
						if M+A7<=K:A8.append((0,M,M+A7))
						else:A8.append((0,M,K));A8.append((K-M,0,M+A7-K))
					elif M+1>=A7:A8.append((0,M-A7+1,M+1))
					else:A8.append((0,0,M+1));A8.append((M+1,K-(A7-M-1),K))
					for(AT,S,T)in A8:
						A9=S//g;AA=(T-1)//g
						for Q in A(A9,AA+1):
							i=Q*g;A0=i+g;Bv=i>=S and A0<=T
							for Y in Bt[AR[Q]:AR[Q+1]]:
								if BK==1:
									j=Y;Aq=j>>c;w=j&BD
									if Bv:v.add((AT+(i-S)+Aq,w))
									elif A9!=AA and Q==A9:
										if j>=S-i<<c:v.add((AT+(i-S)+Aq,w))
									elif A9!=AA and Q==AA:
										if j<T-i<<c:v.add((AT+(i-S)+Aq,w))
									else:
										BL=L(S,i)-i<<c;BM=U(T,A0)-L(S,i)<<c
										if 0<=j-BL<BM:v.add((AT+(i-S)+Aq,w))
								else:
									AU=Y^CU;CW=AU>>c;w=AU&BD;Ar=AT+(T-A0)+CW
									if Bv:v.add((Ar,w))
									elif A9!=AA and Q==AA:
										if AU>=A0-T<<c:v.add((Ar,w))
									elif A9!=AA and Q==A9:
										if AU<A0-S<<c:v.add((Ar,w))
									else:
										BL=A0-U(T,A0)<<c;BM=U(T,A0)-L(S,i)<<c
										if 0<=AU-BL<BM:v.add((Ar,w))
					if v!=Bu:BC+=1
	BN=0;Bw=0;As=(1<<47)-1;Bx=(1<<64)-1
	for At in(6,7,8,9):
		O=1<<At;BO=O-1;CX=BO<<47
		for K in(O*3+17,O*5+3):
			CY=(K+O-1)//O;Au=180;By=[0]+[Z.randrange(K)for aK in A(Au)];BP=[0]*CY;A4=[0]*(Au+1)
			for u in A(1,Au+1):f=By[u];Q=f>>At;Y=(f&BO)<<47|u;A4[u]=BP[Q];BP[Q]=Y
			for m in(19,O//2+5,O+11,U(K-1,2*O+37)):
				if m>=K:continue
				for BQ in(1,-1):
					for M in(0,O-3,K-1,K-m//2):
						Bw+=1;Bz=W()
						for u in A(1,Au+1):
							V=(By[u]-M)*BQ%K
							if V<m:Bz.add((V,u))
						BR=W()
						if BQ==1:B_=[(M,M+m-1,0)]if M+m<=K else[(M,K-1,0),(0,M+m-1-K,K-M)]
						else:B_=[(M-m+1,M,m-1)]if M>=m-1 else[(0,M,M),(K-(m-1-M),K-1,m-1)]
						for(AB,AV,C0)in B_:
							AW=AB>>At;AX=AV>>At
							for Q in A(AW,AX+1):
								Y=BP[Q]
								if Y==0:continue
								if BQ==1:
									q=C0+(Q*O-AB)
									if AW==AX:S=AB-Q*O;T=AV-Q*O+1;b=3
									elif Q==AW:S=AB-Q*O;T=O;b=1
									elif Q==AX:S=0;T=AV-Q*O+1;b=2
									else:S=0;T=O;b=0
									AY=S<<47;BS=T<<47;AF=T-S<<47
									while Y:
										if b==0 or b==1 and Y>=AY or b==2 and Y<BS or b==3 and Y-AY&Bx<AF:BR.add((q+(Y>>47),Y&As))
										Y=A4[Y&As]
								else:
									q=C0+AB-(Q*O+BO)
									if AW==AX:S=O-(AV-Q*O+1);T=O-(AB-Q*O);b=3
									elif Q==AW:S=0;T=O-(AB-Q*O);b=2
									elif Q==AX:S=O-(AV-Q*O+1);T=O;b=1
									else:S=0;T=O;b=0
									AY=S<<47;BS=T<<47;AF=T-S<<47
									while Y:
										j=Y^CX
										if b==0 or b==1 and j>=AY or b==2 and j<BS or b==3 and j-AY&Bx<AF:BR.add((q+(j>>47),j&As))
										Y=A4[j&As]
						if BR!=Bz:BN+=1
	AZ=F._primes_upto(16383)[1:]
	if F._P17_G_CACHE is not N:AZ=AZ+F._P17_G_CACHE
	else:
		for D in F._primes_upto((1<<17)-1):
			if D>=16384:
				if F.U_BN254%D==0:AZ.append(D);continue
				BT=F._ppow_h([0,1,0,0,0,0],D,D);BT[1]=(BT[1]-1)%D;C1=F._pgcd(BT,F.H_COEFFS,D)
				if C(C1)>1:
					for Aa in F._roots_mod(C1,D):
						C2=F.U_BN254*B(F.W_COEFFS[A]*E(Aa,A,D)for A in A(6))%D
						if C2==0 or E(C2,(D-1)//2,D)==1:AZ.append(D);break
	Ab=[]
	for _ in A(200):
		Af=[Z.randint(-2**15,2**15)for aw in A(6)];Ag=[Z.randint(-2**17,2**17)for ax in A(6)];AL=[]
		for(D,A1,AK,Ae)in Ac:
			A2=1
			for V in A(6):AM=B(Af[A]*A1[V][A]for A in A(6))%D;q=B(Ag[A]*A1[V][A]for A in A(6))%D;A2=A2*(AK[V]*AM*AM-q*q)%D
			AL.append(A2)
		AC=B0(AL,9)
		if AC==0:continue
		AC>>=(AC&-AC).bit_length()-1;BU=0
		for D in AZ:
			if AC%D==0:
				while AC%D==0:AC//=D;BU+=1
				BU+=1
		Ab.append(BU)
	C3=B(Ab)/C(Ab);CZ=L(Ab);BV=0;Ca=0x8000800080008000;Ci=65535;BW=F.rand_below;BX=Z.getrandbits
	for Cb in A(48):
		AJ=[[BW(BX,32768)for aI in A(128)]for ay in A(8)]
		if Cb%3==0:AJ[Z.randrange(8)][Z.randrange(128)]|=32768
		BY=[32768-(3000+BW(BX,600))for az in A(128)];C4=[aA+48 for aA in BY];C5=[BW(BX,9)==0 for aB in A(8)];Cc=lambda r,i:r[4*i]|r[4*i+1]<<16|r[4*i+2]<<32|r[4*i+3]<<48;C6=0;BZ=0
		for Aa in AJ:
			for z in A(32):C6|=Cc(Aa,z)
		BZ+=385;C7=W()
		if C6&Ca:
			for(Av,Aa)in X(AJ):
				for(Cd,Ce)in X(Aa):
					if Ce&32768:C7.add((Av,Cd))
		Cf={(aC,aD)for aC in A(8)for aD in A(128)if AJ[aC][aD]>=32768}
		for Av in A(8):AJ[Av]=n(C4 if C5[Av]else BY)
		BZ+=164
		if C7!=Cf or BZ!=549 or A5(AJ[A]!=(C4 if C5[A]else BY)for A in A(8)):BV+=1
	Cg=3.<B8<5. and 1.5<B9<3. and An==BB and BC==0 and Ao==0 and BN==0 and BV==0 and C3<=16.;C8=8122;C9=16654;Ch=Ax==0 and Ai==0 and B2==32 and Cg and Bn>557. and C8<=F.C_NORM_G and C9<=F.C_NORM_F;return h(ok=Ch,split_primes_9=C(Aw),prod_bits_9=Bn,crt_norm_samples=20,crt_norm_errs=Ax,max_nf_bits=Az,max_ng_bits=Ay,op_ng_wc=C8,op_nf_wc=C9,cyclic_mp_rt_errs=Ai,srt_5trans_leaves_exact=B2,murphy_af_unram_500=J(B8,4),murphy_ag_unram_500=J(B9,4),subs_adc_occbm_ok=An==BB,flat_b2_subband_errs=BC,lex_tmp_kta_errs=Ao,chained_1pass_b2_tests=Bw,chained_1pass_b2_errs=BN,scan8_blocks=48,scan8_errs=BV,scan8_ops_per_1024=549,strip_att_n=C(Ab),strip_att_mean=J(C3,3),strip_att_max=CZ)
def verify_live_genuine_c6_relation_filter(seed=20261120):
	import random as BN,model_bn254_v6 as g,verify_bn254_v6 as Ao;Az,F=13,28393;AZ=E(Az,10,F);i=F*F;assert g._isp(F)and F==Az**4-Az**2+1 and(AZ*AZ-AZ+1)%F==0;Aa=g.H_COEFFS;BO=[aa*Aa[aa]for aa in A(1,7)];o,AO=Ao._eval_poly_mod,Ao._eval_zk_mod;AP=[[[0]*6 for aA in A(6)]for ab in A(6)]
	for AF in A(6):AP[0][AF][AF]=1
	BP=0
	for a in A(1,6):
		AP[a][0]=AP[0][a]
		for AF in A(a,6):AP[a][AF]=AP[AF][a]=g._OH_U[BP];BP+=1
	A_=G.sumprod;BQ=[[[AP[aM][ac][aB]for aM in A(6)]for aB in A(6)]for ac in A(6)]
	def Ap(u,v):B=[[A_(u,BQ[aY][a0])for a0 in A(6)]for aY in A(6)];return[A_(v,[B[a1][aZ]for a1 in A(6)])for aZ in A(6)]
	def BR(v):return y(Ao._det_int([[A_(v,BQ[a2][a3])for a3 in A(6)]for a2 in A(6)]))
	CD=[-1,-1,-1,-1,1,-1];CE=[2*ad for ad in CD];CF=[3,-1,-1,-1,1,-1];BS=(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1);Aq={};Ar={};B0=W();Ab={}
	def BT(k):
		if k not in Ab:Ab[k]=C(Ab)
		return Ab[k]
	for D in A(5,601):
		if not Ao._is_prime(D)or D in(13,19):continue
		e=[aC for aC in A(D)if o(Aa,aC,D)==0]
		if C(e)in(3,6):
			BU=E(3,-1,D);BV=W();Aq[D]=u(e)
			for f in u(e):
				if f in BV:continue
				p=o(g.SIGMA_NUM,f,D)*BU%D;z=o(g.SIGMA_NUM,p,D)*BU%D
				if C({f,p,z})!=3:continue
				BV.update((f,p,z));BW=o(g.W_COEFFS,f,D);B1=2*BW%D
				if B1 and E(B1,(D-1)//2,D)==1:
					CG=E(B1,-1,D);Z=AD(A for A in A(1,D)if A*A%D==CG);AQ=BT(('tau6',0,D,f,Z))
					for(Ac,(AR,B2))in X(((f,Z),(z,D-Z),(p,Z),(f,D-Z),(z,Z),(p,D-Z))):Ar[0,D,AR,B2]=AQ,Ac
				B3=(BW+4)%D
				if B3 and E(B3,(D-1)//2,D)==1:
					Z=AD(A for A in A(1,D)if A*A%D==B3);AQ=BT(('tau6',1,D,f,Z))
					for(Ac,(AR,B2))in X(((f,Z),(z,D-Z),(p,Z),(f,D-Z),(z,Z),(p,D-Z))):Ar[1,D,AR,B2]=AQ,Ac
		elif C(e)==0:
			if g._ppow_h([0,1,0,0,0,0],D**3,D)==[0,1,0,0,0,0]:B0.add(D)
	BX=u(n(Aq.keys())+n(B0));CH=W(BX)
	def BY(n):
		while n%2==0:n//=2
		while n%3==0:n//=3
		B={}
		for A in BX:
			if A*A>n:break
			if n%A==0:
				C=0
				while n%A==0:n//=A;C+=1
				B[A]=C
		if n>1:
			if n in CH:B[n]=B.get(n,0)+1
			else:return
		return B
	CI=g._roots_mod(Aa,F);BZ=E(3,-1,F)
	def CJ(s):return(s-o(Aa,s,i)*E(o(BO,s,i),-1,i))%i
	Ba=[];Bb=W()
	for f in CI:
		if f in Bb:continue
		p=o(g.SIGMA_NUM,f,F)*BZ%F;z=o(g.SIGMA_NUM,p,F)*BZ%F;Bb.update((f,p,z));Ba.append([CJ(aN)for aN in(f,p,z)])
	def Bc(u,v,c,md):return(u[0]*v[0]+c*u[1]*v[1])%md,(u[0]*v[1]+u[1]*v[0])%md
	def CK(u,e,c,md):
		B=1,0;A=u
		while e>0:
			if e&1:B=Bc(B,A,c,md)
			A=Bc(A,A,c,md);e>>=1
		return B
	As=[]
	for Ad in Ba:
		Bd=o(g.W_COEFFS,Ad[0],i);Be=E(2*Bd%i,-1,i);Bf=(Bd+4)%i;CL=F-1 if E(Be%F,(F-1)//2,F)==1 else i-1;CM=F-1 if E(Bf%F,(F-1)//2,F)==1 else i-1
		if not As:As.append((Ad,Be,CL))
		As.append((Ad,Bf,CM))
	AG=[]
	for a in A(6):
		b=[0]*6;b[a]=1;AG.append(O(b))
		for AF in A(a+1,6):
			for CN in(1,-1):b=[0]*6;b[a]=1;b[AF]=CN;AG.append(O(b))
	CO=BN.Random(seed)
	for _ in A(300):
		b=O(CO.choice((-1,0,0,0,1))for A in A(6))
		if A5(b):AG.append(b)
	AG=n(h.fromkeys(AG));Y=[];B4=[];Bg=W();Bh=0;CP=E(AZ,2,F)
	for AS in AG:
		A6=Ap(AS,AS);CQ=Ap(CE,A6)
		for AT in AG:
			if AS==AT:continue
			v=Ap(AT,AT);Bi=BR([CQ[aS]-v[aS]for aS in A(6)])
			if Bi<=1:continue
			Bj=BY(Bi)
			if Bj is N:continue
			CR=Ap(CF,v);Bk=BR([A6[aT]-CR[aT]for aT in A(6)])
			if Bk<=1:continue
			Bl=BY(Bk)
			if Bl is N:continue
			AH={};Ae=T;Bm=0
			for(At,CS)in((0,Bj),(1,Bl)):
				Bn=1 if At==0 else-1
				for(D,AU)in CS.items():
					if D in B0:continue
					Bo=0
					for Af in Aq[D]:
						Bp=AO(AS,Af,D);Bq=AO(AT,Af,D)
						if Bq==0:
							if Bp==0:Ae=R;break
							continue
						Z=Bp*E(Bq,-1,D)%D
						if(At,D,Af,Z)in Ar:
							b=1
							if AU>1:
								m=D**AU;AI=Af;AV=Z
								for _ in A(AU):AI=(AI-o(Aa,AI,m)*E(o(BO,AI,m),-1,m))%m
								Br=o(g.W_COEFFS,AI,m);CT=E(2*Br%m,-1,m)if At==0 else(Br+4)%m
								for _ in A(AU):AV=(AV-(AV*AV-CT)*E(2*AV,-1,m))%m
								Bs=(AO(AS,AI,m)-AO(AT,AI,m)*AV)%m;b=0
								while b<AU and Bs%D==0:b+=1;Bs//=D
							AQ,Ac=Ar[At,D,Af,Z];CU,CV=BS[Ac];CW,CX=AH.get(AQ,(0,0));AH[AQ]=CW+Bn*b*CU,CX+Bn*b*CV;Bo+=b;Bm+=b
					if not Ae:break
					if C(Aq[D])==6 and Bo!=AU:Ae=R;break
				if not Ae:break
			if not Ae:continue
			AH={aP:aO for(aP,aO)in AH.items()if aO!=(0,0)}
			if C(AH)>=2:
				Bt=O(u(AH.items()))
				if Bt not in Bg:
					Bg.add(Bt);Y.append(AH);Bh+=Bm;Bu=[]
					for(Ad,CY,CZ)in As:
						AW=0;Ag=1
						for AR in Ad:Ca=CK((AO(AS,AR,i),-AO(AT,AR,i)%i),CZ,CY,i);AW=(AW+Ag*(Ca[1]//F))%F;Ag=Ag*CP%F
						Bu.append(AW)
					B4.append(Bu)
			if C(Y)>=160:break
		if C(Y)>=160:break
	A7=lambda a,b:(a+b*AZ)%F;Ah=u({aE for aD in Y for aE in aD});Cb={af:ae for(ae,af)in X(Ah)};AJ=C(Ah);Au=AJ+3;c=[]
	for(q,Cc)in l(Y,B4):
		B5=[0]*Au
		for(K,Cd)in q.items():B5[Cb[K]]=A7(*Cd)
		for(Ce,Cf)in X(Cc):B5[AJ+Ce]=Cf
		c.append(B5)
	A0=0;B6=[]
	for Ai in A(Au):
		B7=AD((A for A in A(A0,C(c))if c[A][Ai]!=0),N)
		if B7 is N:continue
		c[A0],c[B7]=c[B7],c[A0];Cg=E(c[A0][Ai],-1,F);c[A0]=[aF*Cg%F for aF in c[A0]]
		for a in A(C(c)):
			if a!=A0 and c[a][Ai]:Ch=c[a][Ai];c[a]=[(aU-Ch*aV)%F for(aU,aV)in l(c[a],c[A0])]
		B6.append(Ai);A0+=1
	Bv=[ag for ag in A(Au)if ag not in B6];B8=C(Bv);Av=[0]*Au
	if B8==1:
		Bw=Bv[0];Av[Bw]=1
		for(a,Ci)in X(B6):Av[Ci]=-c[a][Bw]%F
	AX={Ah[ah]:Av[ah]for ah in A(AJ)};Cj=Av[AJ:];AY={ai:-B(B4[ai][A]*Cj[A]for A in A(3))%F for ai in A(C(Y))};Bx=B(1 for(A,C)in X(Y)if B(A7(*B)*AX[A]for(A,B)in C.items())%F!=AY[A]);B9=C(Y);j={}
	for(Q,q)in X(Y):
		for K in q:j.setdefault(K,W()).add(Q)
	A1=W(A(B9));AK=[aj for(aj,ak)in j.items()if C(ak)==1]
	while AK:
		K=AK.pop();e=[aG for aG in j.get(K,())if aG in A1]
		if C(e)==1:
			Q=e[0];A1.remove(Q)
			for M in Y[Q]:
				j[M].discard(Q)
				if C(j[M])==1:AK.append(M)
	Aw={al for(al,am)in j.items()if C(am)>=2};BA=[];Ck=U(C(A1),H(1.2*C(Aw))+4)
	for Q in u(A1,key=lambda i:C(Y[i]),reverse=T):
		if C(A1)<=Ck:break
		if AE(C(j[A])>=2 for A in Y[Q]):
			A1.remove(Q)
			for K in Y[Q]:j[K].discard(Q)
			BB=[aQ for aQ in Y[Q]if C(j[aQ])==1]
			while BB:
				K=BB.pop()
				if C(j[K])==1:
					S=AD(iter(j[K]))
					if AE(C(j[A])>=2 for A in Y[S]if A!=K):
						A1.remove(S);BA.append((K,S,h(Y[S])))
						for M in Y[S]:
							j[M].discard(S)
							if C(j[M])==1:BB.append(M)
	BC={ao for an in A1 for ao in Y[an]};By={ap for(ap,aq,aq)in BA};Bz=C(Aw-(BC|By));BD={ar:AX[ar]for ar in BC};Ax=0
	for(K,Q,B_)in BM(BA):
		Ay=A7(*B_[K])
		if Ay==0:Ax+=1;continue
		BE=B(A7(*B)*BD[A]for(A,B)in B_.items()if A!=K)%F;BD[K]=(AY[Q]-BE)*E(Ay,-1,F)%F
		if BD[K]!=AX[K]:Ax+=1
	V={at:h(Y[at])for at in A1};A2={}
	for(Q,q)in V.items():
		for K in q:A2.setdefault(K,W()).add(Q)
	Cl=W(BC);BF=[];C0=0;BG=0
	for K in u(n(A2.keys()),key=lambda cc:C(A2[cc])):
		e=n(A2.get(K,()))
		if C(e)!=2:continue
		d,S=e
		if d not in V or S not in V:continue
		A3,r=V[d][K]
		if A3*A3+A3*r+r*r!=1:continue
		Aj,AL=A3+r,-r;A6,v=V[S][K];A8,s=Aj*A6-AL*v,Aj*v+AL*A6+AL*v;BH={}
		for(M,(A9,t))in V[d].items():
			if M==K:continue
			Ak,Al=A8*A9-s*t,A8*t+s*A9+s*t;AA,AB=V[S].get(M,(0,0));w,x=AA-Ak,AB-Al
			if(AA,AB)==(0,0):BH[M]=1
			elif(w,x)==(0,0):BH[M]=-1
		if A5(C(A2.get(A,()))-1+BH.get(A,0)<=0 for A in V[d]if A!=K):continue
		BF.append((K,h(V[d]),AY[d]))
		for(M,(A9,t))in V[d].items():
			Ak,Al=A8*A9-s*t,A8*t+s*A9+s*t;AA,AB=V[S].get(M,(0,0));w,x=AA-Ak,AB-Al
			if w==0 and x==0:V[S].pop(M,N);A2[M].discard(S)
			else:V[S][M]=w,x;A2.setdefault(M,W()).add(S)
		Cm=B(A7(B,C)*AX[A]for(A,(B,C))in V[S].items())%F;C1=(AY[S]-A7(A8,s)*AY[d])%F
		if Cm!=C1:BG+=1
		AY[S]=C1;V.pop(d);C0+=1
		for M in n(A2):A2[M].discard(d)
	BI={av for au in V.values()for av in au};C2=C(Cl-(BI|{aW for(aW,aX,aX)in BF}));BJ={aw:AX[aw]for aw in BI};BK=0
	for(Am,C3,Cn)in BM(BF):
		Ay=A7(*C3[Am]);BE=B(A7(*B)*BJ[A]for(A,B)in C3.items()if A!=Am)%F;BJ[Am]=(Cn-BE)*E(Ay,-1,F)%F
		if BJ[Am]!=AX[Am]:BK+=1
	C4=B(C(A)for A in V.values());Co=B(1 for C in V.values()for(A,B)in C.values()if A*A+A*B+B*B==1)
	def Cp(a,b):return U(y(a+A)+y(b-A)+y(A)for A in(0,-a,b))
	def C5(rmap):
		A=0;B=I
		for D in rmap.values():
			for(E,F)in D.values():A+=1;C=Cp(E,F);B+=C*g.C_UNIT_5 if C<=4 else g.C_GEN_1
		return B/L(1,A)
	C6=C5(V);Cq={ay:ax for(ax,ay)in Ab.items()};C7=[P/Cq[az][2]for az in Ah];Cr=B(C7);C8=[];AW=I
	for Ag in C7:AW+=Ag/Cr;C8.append(AW)
	import bisect;C9=[];CA=[];CB=[]
	for Cs in A(4):
		CC=BN.Random(seed+100+Cs);A4=[]
		for q in Y:
			An={}
			for _ in A(C(q)):
				BL=Ah[U(AJ-1,bisect.bisect_left(C8,CC.random()))];Ct,Cu=BS[CC.randrange(6)];AA,AB=An.get(BL,(0,0));w,x=AA+Ct,AB+Cu
				if w==0 and x==0:An.pop(BL,N)
				else:An[BL]=w,x
			if C(An)>=2:A4.append(An)
		AC={};AM=W(A(C(A4)))
		for(Q,q)in X(A4):
			for K in q:AC.setdefault(K,W()).add(Q)
		AK=[aH for(aH,aI)in AC.items()if C(aI)==1]
		while AK:
			K=AK.pop();e=[aR for aR in AC.get(K,())if aR in AM]
			if C(e)==1:
				Q=e[0];AM.remove(Q)
				for M in A4[Q]:
					AC[M].discard(Q)
					if C(AC[M])==1:AK.append(M)
		Cv={aJ for(aJ,aK)in AC.items()if C(aK)>=2};Cw=U(C(AM),H(1.2*C(Cv))+4)
		for Q in u(AM,key=lambda i:C(A4[i]),reverse=T):
			if C(AM)<=Cw:break
			if AE(C(AC[A])>=2 for A in A4[Q]):
				AM.remove(Q)
				for K in A4[Q]:AC[K].discard(Q)
		k={aL:h(A4[aL])for aL in AM};AN={}
		for(Q,q)in k.items():
			for K in q:AN.setdefault(K,W()).add(Q)
		for K in u(n(AN.keys()),key=lambda cc:C(AN[cc])):
			e=n(AN.get(K,()))
			if C(e)!=2:continue
			d,S=e
			if d not in k or S not in k:continue
			A3,r=k[d][K]
			if A3*A3+A3*r+r*r!=1:continue
			Aj,AL=A3+r,-r;A6,v=k[S][K];A8,s=Aj*A6-AL*v,Aj*v+AL*A6+AL*v
			for(M,(A9,t))in k[d].items():
				Ak,Al=A8*A9-s*t,A8*t+s*A9+s*t;AA,AB=k[S].get(M,(0,0));w,x=AA-Ak,AB-Al
				if w==0 and x==0:k[S].pop(M,N);AN[M].discard(S)
				else:k[S][M]=w,x;AN.setdefault(M,W()).add(S)
			k.pop(d)
			for M in A4[d]:AN.get(M,W()).discard(d)
		Cx=B(C(A)for A in k.values());Cy=B(1 for C in k.values()for(A,B)in C.values()if A*A+A*B+B*B==1);C9.append(C(k));CA.append(Cy/L(1,Cx));CB.append(C5(k))
	Cz=B(C9)/4.;C_=B(CA)/4.;D0=B(CB)/4.;D1=B9>=150 and AJ>=25 and B8==1 and Bx==0 and Bz==0 and C2==0 and Ax==0 and BK==0 and BG==0 and C(V)>=C(BI)and C6<=g.K5_CAP;return h(ok=D1,synthetic_incidences=R,zero_nonsplit_inert_aggregation=T,unplanted_f13p12_r28393_chi_kernel=T,unplanted_nullity=B8,unplanted_rel_errs=Bx,genuine_c6_rels=B9,pure_tau6_orbit_cols=AJ,sm_chi_cols=3,post_peel_cols=C(Aw),active_tau6_orbits=C(Aw),tau6_unit_hits=Bh,repurged_cols=C(By),post_peel_cols_lost=Bz,post_merge_cols_lost=C2,repurge_backsub_errors=Ax,merge_backsub_errors=BK,z6_unit_merges=C0,final_rows=C(V),final_avg_w=J(C4/L(1,C(V)),2),unit_entry_frac=J(Co/L(1,C4),4),gen_k5_v4=J(C6,4),syn_matched_final_rows=J(Cz,2),syn_matched_unit_frac=J(C_,4),syn_matched_k5_v4=J(D0,4),kernel_errors=BG)
import functools as D,math as G,random
from math import ceil,floor,log,log2
@D.lru_cache(N)
def psurv_iv_tab(u_max=45.,m=512):
	B=G.nextafter;K=H(ceil(u_max*m))+2;b=B(P/m,G.inf);c=B(P/m,I);Q=[P]*(m+1);R=[P]*(m+1);J=32;V=[];T=[];d=[I]*(K+1);W=[I]*(K+1)
	def e(bi):
		F=bi*J;C=D=I
		for E in A(F+J-1,F-1,-1):C=B(C+Q[E],G.inf);D=B(D+R[E],I);d[E]=C;W[E]=D
		V.append(C);T.append(D)
	for o in A(m//J):e(o)
	U=Q[m];S=R[m]
	for D in A(m+1,K+1):
		if D<=2*m:f=P-log(D/m);F=B(B(f,G.inf),G.inf);E=L(I,B(B(f,I),I))
		else:
			M=D//J;N=(D-m)//J
			if D==M*J:g=G.fsum(V[N:M]);h=G.fsum(T[N+1:M])+W[D-m+1]
			else:g=U+d[D-m]+G.fsum(V[N+1:M]);i=D-m+1;h=S+G.fsum(T[N+1:M])if i==(N+1)*J else S+W[i]+G.fsum(T[N+1:M])
			F=g/D;E=h/(D-1)
			for j in A(16):F=B(F,G.inf);E=B(E,I)
			E=L(I,E)
		Q.append(F);R.append(E)
		if(D+1)%J==0:e((D+1)//J-1);U=S=I
		else:U=B(U+F,G.inf);S=B(S+E,I)
	l=[I]*K;n=[I]*K;X=[P]*(K+1);Y=[P]*(K+1)
	for C in A(m,K):
		F=m/k(C)if C<2*m else X[C+1-m]*m/k(C);E=m/k(C+1)if C<2*m else Y[C-m]*m/k(C+1)
		for j in A(4):F=B(F,G.inf);E=B(E,I)
		l[C]=F;n[C]=L(I,E);Z=X[C]+b*F;a=Y[C]+c*E
		for j in A(4):Z=B(Z,G.inf);a=B(a,I)
		X[C+1]=Z;Y[C+1]=a
	return m,b,c,O(Q),O(R),O(l),O(n),O(Q[::-1]),O(R[::-1])
@D.lru_cache(N)
def psurv_iv(m_bits,B,T,mesh=512):
	Y=m_bits
	if T<I:return I,I
	if Y<=T:return P,P
	S=Y/B;f=T/B;g=L(45.,ceil(S)+2.);D,h,i,Z,j,a,b,k,l=psurv_iv_tab(g,mesh);O=H(floor(S*D));E=H(ceil(S*D));X=U(S,f);F=C(Z);c=Z[O];d=j[E]if E<F else I
	if X<=P:return U(P,d),U(P,c)
	V=H(ceil(X*D));W=H(floor(X*D));e=G.sumprod;J=O-1 if O-1<V else V;Q=I
	if J>D:Q+=e(k[F-O+D:F-O+J],a[D:J])
	if V>(J if J>D else D):Q+=G.fsum(a[L(D,J):V])
	K=D if E-D<F else E-F+1;M=E if E<W else W;R=I
	if M>K:R+=e(l[F-1-E+K:F-1-E+M],b[K:M])
	if W>(M if M>K else K):R+=G.fsum(b[L(K,M):W])
	N=G.nextafter
	for m in A(8):Q=N(Q,G.inf);R=N(R,I)
	return U(P,L(I,N(d+N(i*L(I,R),I),I))),U(P,N(c+N(h*Q,G.inf),G.inf))
def psurv512(m_bits,B,T,mesh=512):return psurv_iv(m_bits,B,T,mesh)[1]
def verify_v43_descent(D,M,segs):
	b=segs;import re;BF=b.get('desc_zero16',-1);BG=b.get('desc_scan8_u8',-1);CK=b.get('desc_sieve_hit',-1);BH=b.get('mont5_cios_prologue',-1);BI=b.get('mont5_cios_round',-1);BJ=b.get('mont5_cios_epilogue',-1);BK=b.get('mont2_cios',-1);BL=b.get('bingcd_gen_loop',-1);BM=b.get('modinv_gen_loop',-1);Af=BH+5*BI+BJ;CL=b.get('cios_gen_mul_inner',-1);CM=b.get('cios_gen_red_inner',-1);Ag=D.MONT5;Ah=25.;d=3*Ag+5*Ah;x=2*Ag+3*Ah;g=5e1;AJ=M.P_BN254-2;BN=(4+(AJ.bit_length()-1)+(AJ.bit_count()-1))*Ag+2*Ah;BO=254*(3*d+4*x+10*g)+254*(12*d+4*x+10*g)+BN+3*d+x;BP=8*d+15*g;Ai=3*BP+d+30*g;CN=12*d+3*x+15*g+BN;CO=4*BP+d+CN+10*g;BQ=M.P_BN254**4-M.P_BN254**2+1;AK=BQ//M.R_BN254;CP=AK.bit_length()-1+(AK.bit_count()-1);BR=65*(2*Ai+6*d+6*x+20*g)+67*(Ai+12*d+4*x+20*g)+((CP+4)*Ai+CO+2e5);CQ=AJ.bit_length()==254 and AJ.bit_count()==110 and BQ%M.R_BN254==0 and AK.bit_length()==761 and AK.bit_count()==389 and BO<=5422256. and BO<=M.G2_OUTPUT_BOUND and BR<=40297005. and BR<=2e4*2.**11;Aj=BF==18 and BG==22 and CK==9 and BH==18 and BI==48 and BJ==21 and Af==279 and Af+7<=D.MONT5 and BK==47 and BL==49 and BM==91 and CL==11 and CM==11 and BF/256.<=D.ZERO and BG/64.<=D.SCAN and CQ;CR=Path(__file__).resolve().parent/'kernels_v32.s';CS=CJ.ko(CR).read().splitlines();Ak=[];Al=R;AL=[];Am=[];AM=[];An=[];AN=[];BS={};BT={};BU={};a=N;AO=0;Ao=[]
	for CT in CS:
		O=CT.split('//',1)[0].strip()
		if not O:continue
		if O.startswith('.macro MONT5_ROUND'):Al=T;continue
		if Al:
			if O.startswith('.endm'):Al=R
			else:Ak.append(O)
			continue
		if O=='seg_mont5_cios:':a='m5';continue
		if O=='seg_mont2_cios:':a='m2';continue
		if O=='seg_bingcd_gen:':a='bg';continue
		if O=='seg_cios_gen:':a='cg';continue
		if O=='seg_modinv_gen:':a='mi';continue
		if a=='m5':
			if O.startswith('.rept '):AO=H(O.split()[1]);Ao=[];continue
			if O=='.endr':
				for y in A(AO):
					for BV in Ao:AL.extend(Ak if BV==BE else[BV])
				AO=0;continue
			if AO:Ao.append(O);continue
			if O.endswith(v)or O.startswith('.'):continue
			AL.extend(Ak if O==BE else[O])
			if O==w:a=N
		elif a=='m2':
			if O.endswith(v)or O.startswith('.'):continue
			Am.append(O)
			if O==w:a=N
		elif a=='bg':
			if O.endswith(v):BS[O[:-1]]=C(AM);continue
			if O.startswith('.'):continue
			AM.append(O)
			if O==w:a=N
		elif a=='cg':
			if O.endswith(v):BT[O[:-1]]=C(An);continue
			if O.startswith('.'):continue
			An.append(O)
			if O==w:a=N
		elif a=='mi':
			if O.endswith(v):BU[O[:-1]]=C(AN);continue
			if O.startswith('.'):continue
			AN.append(O)
			if O==w:a=N
	F=(1<<64)-1
	def S(s):s=s.strip().rstrip(r);return 32 if s=='xzr'else 31 if s=='sp'else H(s[1:])
	Ap={};CU={'mul':26,'umulh':27,'orr':28,'add':29,'sub':30,'adds':31,'adcs':32,'adc':33,'subs':34,'sbcs':35,'csel':36};BW={'b':1,'b.ne':2,'b.lo':3,'b.hs':4,'b.pl':5};BX={'add':13,'sub':14,'lsr':15,'lsl':16,'subs':17}
	def CV(prog,labels):
		G=labels;I=id(prog)
		if I in Ap:return Ap[I]
		C=[]
		for O in prog:
			B,U,D=O.partition(' ');D=D.strip()
			if B==w:C.append((0,0,0,0,0))
			elif B in BW:C.append((BW[B],G[D],0,0,0))
			elif B in('cbz','cbnz'):J,K=D.split(r);C.append((6 if B=='cbz'else 7,S(J),G[K.strip()],0,0))
			elif B in('tbz','tbnz'):J,Q,K=D.split(r);C.append((8 if B=='tbz'else 9,S(J),H(Q.strip().lstrip('#')),G[K.strip()],0))
			elif B=='cmp':E,F=D.split(r);C.append((10,S(E),S(F),0,0))
			elif B=='mov':E,F=D.split(r);L=F.strip();C.append((11,S(E),H(L.lstrip('#')),0,0)if L.startswith('#')else(12,S(E),S(L),0,0))
			elif B in BX and'#'in D:E,F,R=D.split(r);C.append((BX[B],S(E),S(F),H(R.strip().lstrip('#')),0))
			elif B in('stp','ldp'):A=re.match('(\\w+),\\s*(\\w+),\\s*\\[(\\w+),\\s*#(\\d+)\\]$',D);C.append((18 if B=='stp'else 19,S(A.group(1)),S(A.group(2)),S(A.group(3)),H(A.group(4))))
			elif B in('str','ldr'):
				M=20 if B=='str'else 23;A=re.match('(\\w+),\\s*\\[(\\w+),\\s*(\\w+),\\s*lsl\\s*#(\\d+)\\]$',D)
				if A:C.append((M,S(A.group(1)),S(A.group(2)),S(A.group(3)),H(A.group(4))))
				else:
					A=re.match('(\\w+),\\s*\\[(\\w+)\\],\\s*#(\\d+)$',D)
					if A:C.append((M+1,S(A.group(1)),S(A.group(2)),H(A.group(3)),0))
					else:A=re.match('(\\w+),\\s*\\[(\\w+)(?:,\\s*#(\\d+))?\\]$',D);C.append((M+2,S(A.group(1)),S(A.group(2)),H(A.group(3)or 0),0))
			else:
				E,F,T=[ap.strip()for ap in D.split(r)[:3]];P=CU.get(B)
				if P is N:raise ValueError(f"unsupported op: {O}")
				C.append((P,S(E),S(F),S(T),0))
		Ap[I]=C;return C
	def A1(prog,x_init,mem,labels=N,max_ops=2000000):
		L=mem;Q=CV(prog,labels or{});A=n(x_init)+[0];I=M=N=J=P=0
		while P<max_ops:
			B,C,D,G,O=Q[J];J+=1;P+=1
			if B==0:return A[:32],P
			elif B==1:J=C
			elif B==2:
				if not M:J=C
			elif B==3:
				if not I:J=C
			elif B==4:
				if I:J=C
			elif B==5:
				if not N:J=C
			elif B==6:
				if A[C]==0:J=D
			elif B==7:
				if A[C]!=0:J=D
			elif B==8:
				if not A[C]>>D&1:J=G
			elif B==9:
				if A[C]>>D&1:J=G
			elif B==10:E=A[C]+(A[D]^F)+1;I=E>>64;M=H(E&F==0);N=E>>63&1
			elif B==11:A[C]=D&F;A[32]=0
			elif B==12:A[C]=A[D];A[32]=0
			elif B==13:A[C]=A[D]+G&F;A[32]=0
			elif B==14:A[C]=A[D]-G&F;A[32]=0
			elif B==15:A[C]=A[D]>>G&F;A[32]=0
			elif B==16:A[C]=A[D]<<G&F;A[32]=0
			elif B==17:E=A[D]+(G^F)+1;I=E>>64;M=H(E&F==0);N=E>>63&1;A[C]=E&F;A[32]=0
			elif B==18:K=A[G]+O;L[K]=A[C];L[K+8]=A[D]
			elif B==19:K=A[G]+O;A[C]=L[K];A[D]=L[K+8];A[32]=0
			elif B==20:L[A[D]+(A[G]<<O)]=A[C]
			elif B==21:K=A[D];L[K]=A[C];A[D]=K+G&F;A[32]=0
			elif B==22:L[A[D]+G]=A[C]
			elif B==23:A[C]=L[A[D]+(A[G]<<O)];A[32]=0
			elif B==24:K=A[D];A[C]=L[K];A[D]=K+G&F;A[32]=0
			elif B==25:A[C]=L[A[D]+G];A[32]=0
			elif B==26:A[C]=A[D]*A[G]&F;A[32]=0
			elif B==27:A[C]=A[D]*A[G]>>64&F;A[32]=0
			elif B==28:A[C]=A[D]|A[G];A[32]=0
			elif B==29:A[C]=A[D]+A[G]&F;A[32]=0
			elif B==30:A[C]=A[D]-A[G]&F;A[32]=0
			elif B==31:E=A[D]+A[G];I=E>>64;M=H(E&F==0);N=E>>63&1;A[C]=E&F;A[32]=0
			elif B==32:E=A[D]+A[G]+I;I=E>>64;M=H(E&F==0);N=E>>63&1;A[C]=E&F;A[32]=0
			elif B==33:A[C]=A[D]+A[G]+I&F;A[32]=0
			elif B==34:E=A[D]+(A[G]^F)+1;I=E>>64;M=H(E&F==0);N=E>>63&1;A[C]=E&F;A[32]=0
			elif B==35:E=A[D]+(A[G]^F)+I;I=E>>64;M=H(E&F==0);N=E>>63&1;A[C]=E&F;A[32]=0
			elif B==36:A[C]=A[D]if I==1 else A[G];A[32]=0
		raise RuntimeError('max_ops exceeded')
	V=random.Random(20261143);BY=C(AL)==279;BZ=C(Am)==47;Ba=C(AM)==49;Bb=C(AN)==BM;Bc=T
	for Y in A(32):
		K=V.getrandbits(318)|1<<318|1;CW=-E(K,-1,1<<64)&F;Bd=K-1 if Y==0 else V.randrange(K);AP=K-1 if Y==0 else V.randrange(K);A2={};A3,A4,A6,Be,Bf=4096,8192,12288,16384,32768
		for Q in A(5):A2[A3+8*Q]=Bd>>64*Q&F;A2[A4+8*Q]=AP>>64*Q&F;A2[Be+8*Q]=K>>64*Q&F
		e=[0]*32;e[5],e[6],e[7],e[8],e[9],e[31]=A3,A4,CW,A6,Be,Bf
		for Bg in A(19,28):e[Bg]=0xdeadbeef0000+Bg
		Bh,CX=A1(AL,e,A2);CY=B(A2[A6+8*A]<<64*A for A in A(5))
		if CX!=279 or CY!=Bd*AP*E(1<<320,-1,K)%K or Bh[31]!=Bf or A5(Bh[A]!=e[A]for A in A(19,28)):BY=R
		f=(1<<128)-1 if Y==0 else(1<<127)+1 if Y==1 else V.getrandbits(126)|1<<126|1;CZ=-E(f,-1,1<<64)&F;Aq=f-1 if Y==0 else V.randrange(f);Ar=f-1 if Y==0 else V.randrange(f);i=[0]*32;i[0],i[1],i[2],i[3],i[4],i[5],i[6]=Aq&F,Aq>>64&F,Ar&F,Ar>>64&F,f&F,f>>64&F,CZ;Bi,Ca=A1(Am,i,{})
		if Ca!=47 or Bi[0]|Bi[1]<<64!=Aq*Ar*E(1<<128,-1,f)%f:BZ=R
	for W in(2,5,6,8,12,20):
		j=64*W;Cb=H(D.lehmer_gcd_ir(j))
		for Y in A(3):
			s=V.getrandbits(j)|1;AQ=0 if Y==0 else s if Y==1 else V.getrandbits(j)
			if Y==2:
				AR=V.getrandbits(U(60,j//2))|1;s=s//AR*AR;AQ=AQ//AR*AR
				if s&1==0:s|=1
			z={};A7,A8=4096,12288
			for Q in A(W):z[A7+8*Q]=s>>64*Q&F;z[A8+8*Q]=AQ>>64*Q&F
			AS=[0]*32;AS[0],AS[1],AS[2]=A7,A8,W;Cc,Cd=A1(AM,AS,z,BS);Ce=B(z[Cc[0]+8*A]<<64*A for A in A(W))
			if Ce!=G.gcd(s,AQ)or Cd>Cb:Ba=R
	for W in(6,8,12,20):
		As=64*W;Cf=H(D.mul_L(W))
		for Y in A(4):
			m=V.getrandbits(As-2)|1<<As-2|1;Cg=-E(m,-1,1<<64)&F;Bj=m-1 if Y==0 else V.randrange(m);Bk=m-1 if Y==0 else V.randrange(m);A9={};A3,A4,AA,Ch,A6=4096,12288,20480,28672,36864
			for Q in A(W):A9[A3+8*Q]=Bj>>64*Q&F;A9[A4+8*Q]=Bk>>64*Q&F;A9[AA+8*Q]=m>>64*Q&F
			o=[0]*32;o[0],o[1],o[2],o[3],o[4],o[5],o[6]=A3,A4,AA,Ch,W,Cg,A6;y,Ci=A1(An,o,A9,BT);At=B(A9[A6+8*A]<<64*A for A in A(W))
			if At!=Bj*Bk*E(1<<As,-1,m)%m or Ci>Cf:Bc=R
	for W in(2,5,6,8,12,20):
		j=64*W;Cj=H(D.modinv_ir(j))
		for Y in(0,2)if W==20 else(0,1,2):
			Z=V.getrandbits(j-1)|1<<j-1|1;t=V.randrange(1,Z)if Y<2 else V.randrange(1,Z//3)*3%Z or 3
			if Y==2 and Z%3!=0:Z=Z//3*3;Z=Z if Z&1 else Z+3
			if t>=Z:t%=Z
			if t==0:t=1
			p={};A7,A8,Bl,Bm,AA=4096,12288,20480,28672,36864
			for Q in A(W):p[A7+8*Q]=Z>>64*Q&F;p[A8+8*Q]=t>>64*Q&F;p[Bl+8*Q]=0;p[Bm+8*Q]=1 if Q==0 else 0;p[AA+8*Q]=Z>>64*Q&F
			u=[0]*32;u[0],u[1],u[2],u[3],u[4],u[5]=A7,A8,W,Bl,Bm,AA;Bn,Ck=A1(AN,u,p,BU);At=B(p[Bn[0]+8*A]<<64*A for A in A(W));Cl=B(p[Bn[3]+8*A]<<64*A for A in A(W));Bo=G.gcd(t,Z)
			if At!=Bo or Ck>Cj or Bo==1 and Cl*t%Z!=1:Bb=R
	Aj=Aj and BY and BZ and Ba and Bb and Bc;Cm=BL;Au=M.load_data()['descent'];Cn=Au['plan59'];Co={H(aa):ab for(aa,ab)in Au['cells'].items()};Cp=L(A[1]for A in D.QBINS);Cq=L(A[2]for A in D.QBINS);Av=Aw=AB=Ax=I
	for(Bp,Bq,Br,y,y,y,y,Cr,y)in Cn:Bs=log2(Co[Cr]);Bt=(Bq+Bs)/12.;Av=L(Av,D.DF*Bt+Cp-(Bp if Br=='f'else I));Aw=L(Aw,12.*Bt+D.LOG2U6+Cq-(Bp if Br=='g'else I));AB=L(AB,k(Bq));Ax=L(Ax,Bs)
	Ay=L(ceil(Av),ceil(Aw),ceil(Au['boot_dp']['ng_hi']));Bu=L(Ay,2560);AT=ceil(Bu/6e1)+1;Cs=AT/8.*(M.C_NORM_G+M.C_NORM_F)+6e1*AT*AT;Ct=309*(ceil(Bu/64.)*5+10);Bv=D.mul_L(L(4,ceil(AB/64.)));Cu=AT>=44 and Cs<D.NORMS and Ct<D.TDIV and 64.*AB*Bv<D.PRIME and 6.*AB*45.*Bv<D.IDEAL and Ax<=82.;Az=T;Bw=P
	for z in A(400,1900,100):
		Cv,A_=psurv_iv(z,59.,5e2,512)
		if not(I<=Cv<=A_<=P and A_<=Bw+1e-15):Az=R
		Bw=A_
	Bx=V.randrange(1<<592,1<<593)|1;By=V.randrange(1<<591,1<<592)|1;K=Bx*By
	def q(A,B):
		D=[0]*(C(A)+C(B)-1)
		for(F,E)in X(A):
			if E:
				for(G,H)in X(B):D[F+G]+=E*H
		return[aq%K for aq in D]
	def B0(arr,w,p,inv=R):
		B=C(arr);Q=E(w,-1,p)if inv else w;R=B.bit_length()-1;D=[arr[H(f"{ar:0{R}b}"[::-1],2)]%p for ar in A(B)];I=1
		while I<B:
			K=I;I<<=1;S=E(Q,B//I,p);G=1;M=[]
			for U in A(K):M.append((G,(G<<64)//p));G=G*S%p
			for L in A(0,B,I):
				for J in A(K):G,T=M[J];N=D[L+J];O=D[L+J+K];P=(O*G-(O*T>>64)*p&F)%p;D[L+J]=(N+P)%p;D[L+J+K]=(N-P)%p
		return[at*E(B,-1,p)%p for at in D]if inv else D
	Bz=[V.randrange(K)for ac in A(16)];B_=[V.randrange(K)for ad in A(16)];AC=ceil((2*K.bit_length()+5)/64.);B1=32*AC;B2=1<<ceil(log2(B1));AU,AD,AF=M.NTT3_PRIMES
	def C0(P):
		B=[0]*B2
		for(D,E)in X(P):
			C=E%K;G=D*AC
			for H in A(AC):B[G+H]=C&F;C>>=64
		return B
	Cw,Cx=C0(Bz),C0(B_);Cy=E(AU,-1,AD);Cz=E(AU*AD%AF,-1,AF);AV=[]
	for c in M.NTT3_PRIMES:
		AW=2
		while E(AW,(c-1)//2,c)==1 or E(AW,(c-1)//3,c)==1:AW+=1
		B3=E(AW,(c-1)//B2,c);C_,D0=B0(Cw,B3,c),B0(Cx,B3,c);AV.append(B0([C_[ao]*D0[ao]%c for ao in A(B2)],B3,c,T))
	C1=[0]*B1;B4=0
	for AX in A(B1):C2,D1,D2=AV[0][AX],AV[1][AX],AV[2][AX];C3=C2+(D1-C2)%AD*Cy%AD*AU;C4=C3+(D2-C3%AF)%AF*Cz%AF*AU*AD+B4;C1[AX]=C4&F;B4=C4>>64
	C5=[]
	for D3 in A(31):
		B5=0
		for D4 in A(AC-1,-1,-1):B5=B5<<64|C1[D3*AC+D4]
		C5.append(B5%K)
	C6=B4==0 and C5==q(Bz,B_)
	def C7():
		while T:
			A=V.randrange(1,K)
			if G.gcd(A,K)==1:return A
	D5=[V.randrange(1,K)for ae in A(32)];D6=[C7()for af in A(32)];D7=[V.randrange(1,K)for ag in A(32)];D8=[C7()for ah in A(32)];AY=D6+D8;AZ=[AY[0]];AG=0
	for B6 in AY[1:]:AZ.append(AZ[-1]*B6%K);AG+=1
	Aa=E(AZ[-1],-1,K);Ab=[0]*64
	for B7 in A(63,0,-1):Ab[B7]=Aa*AZ[B7-1]%K;Aa=Aa*AY[B7]%K;AG+=2
	Ab[0]=Aa;C8=[ai*aj%K for(ai,aj)in l(D5,Ab[:32])];C9=[ak*al%K for(ak,al)in l(D7,Ab[32:])];AG+=64;B8=[1];A0=n(AY);A0[3]=A0[3]*Bx%K;A0[47]=A0[47]*By%K;AH=[A0[0]]
	for B6 in A0[1:]:AH.append(AH[-1]*B6%K)
	Ac,B9,CA=0,63,0
	while Ac<B9:
		BA=(Ac+B9)//2;CA+=1
		if G.gcd(AH[BA],K)>1:B9=BA
		else:Ac=BA+1
	D9=G.gcd(AH[-1],K)==K and CA<=6 and 1<G.gcd(AH[Ac],K)<K
	for AP in C8:B8=q(B8,[-AP%K,1])
	def BB(pts):
		A=pts
		if C(A)==1:return[-A[0]%K,1],N,N,1
		B=C(A)//2;D=BB(A[:B]);E=BB(A[B:]);return q(D[0],E[0]),D,E,C(A)
	CB=BB(C9);CC=CB[0];DA=[(B8[am]-CC[am])%K for am in A(32)];DB=CC[::-1];Ad=[1];BC=1
	while BC<32:Ae=U(2*BC,32);CD=q(DB[:Ae],Ad)[:Ae];CE=[-an%K for an in CD];CE[0]=(2-CD[0])%K;Ad=q(Ad,CE)[:Ae];BC=Ae
	DC=q(DA[::-1],Ad)[:32]
	def BD(Mc,nd):
		D,A,B,C=nd
		if A is N:return[Mc[0]%K]
		return BD(q(Mc,B[0][::-1])[B[3]:C],A)+BD(q(Mc,A[0][::-1])[A[3]:C],B)
	CF=D9 and AG==253 and AG<=256 and AE(B==G.prod((A-B)%K for B in C8)%K for(A,B)in l(C9,BD(DC,CB)));AI=L(D._stage_ops(A,A,A)[3]for A in D.WHEELS);CG=ceil(L(Ay,1225)/64.);DD=ceil((2*L(Ay,1225)+20)/64.);CH=1<<ceil(log2(2*(AI+1)*DD));CI=6*CH*8+4*AI*CG*8*4+4*AI*CG*8;DE=AI==46080 and AE(D._stage_ops(C,E,A)[3]<=AI for(B,C,E,A,B)in D._ECM_BEST_C.values()if A>0);DF=Aj and Cu and Az and C6 and CF and DE and CH<1<<54 and CI<1<<30;return h(ok=DF,cios5=Af,mont2=BK,bingcd_gen=Cm,psurv_iv_ok=Az,kron_ntt3_ok=C6,bls_mpe_ok=CF,ecm_ram_mib=J(CI/(1<<20),2))
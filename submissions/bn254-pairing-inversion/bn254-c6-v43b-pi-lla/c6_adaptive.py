import base64,json,lzma,math
from pathlib import Path
from fractions import Fraction as F
def verify_threshold_counterexample():
	I,E,B=6e1,1.,.1;G=.5*I/E;A=lambda t:(1.+G*(1.-t))**-E
	for J in range(81):F=J/1e2;assert .1*F+.9*F**100<=A(F)+1e-12
	assert .1>A(B)+.05;C=A(.0)+B
	for D in range(2,65):H=min(.8,(D-1)*(1.+1./G)/(D-1+E));C=min(C,min(1.,A(H)/H**(D-1))+B**D)
	assert C>=.1-1e-12;return dict(pgf_psi_invalid=round(A(B),5),true_fail=.1,tr_bound=round(C,5),tr_valid=True)
def verify_c6_adaptive_h5_live():
	import model_bn254_v6 as A,verify_bn254_v6 as j,descent_bn254_v6 as Ag;A1=verify_threshold_counterexample();B=json.loads(lzma.decompress(base64.b85decode((Path(__file__).resolve().parent/'c6_adaptive_h5.b85').read_bytes())).decode());Ah,c,Ai,Q=B['names'],B['cal_tbl'],B['after_by_set'],B['first_stats']
	def d(si,cs,cb):
		D=f"{si}:{cs}:{cb}"
		if D in c:A,B,C=c[D];return B/A,C/A
		A=B=C=0
		for(E,(F,G,H))in c.items():
			K,I,J=map(int,E.split(':'))
			if I==cs and J==cb:A+=F;B+=G;C+=H
		return B/A,C/A
	A2={};A3=A4=A5=A6=e=A7=A8=A9=0
	for(H,Aj)in enumerate(Ah):
		R=Ai[H];AA=len(R);AB=sum(A[1]for A in R);AC=sum(A[3]for A in R);AD=sum(A[4]for A in R);f=g=h=i=.0
		for(I,Ak,Al,I,I,Am,AE)in R:AF=[Am]+([AE]if AE else[]);f+=Ak;g+=sum(d(H,A>>5,A&31)[0]for A in AF);h+=Al;i+=sum(d(H,A>>5,A&31)[1]for A in AF)
		AG,AH=f/g,h/i;assert AG<.9 and .9<AH<1.1;A2[Aj]=dict(roots=Q[H][2],first_n=Q[H][0],first_fail=Q[H][1],after_n=AA,after_fail=AB,oe_fail=round(AG,3),oe_work=round(AH,3),sh_prev=AC,sh_run=AD);e+=AA;A7+=AB;A3+=f;A4+=g;A5+=h;A6+=i;A8+=AC;A9+=AD
	k,AI=A3/A4,A5/A6;assert e==1856 and k<.85 and .95<AI<1.05;L=B['conc_rows'];assert len(L)==10 and max(A[5]for A in L)<15. and max(A[6]for A in L)<3.;M,C,S,AJ=int(B['u_t']),int(B['p_t']),int(B['lx_b']),int(B['ns_max']);AK=1<<S;N,O=A.ZK_NUM,A.W_COEFFS;AL=[-4-A if B==0 else-A for(B,A)in enumerate(O)];An=A._ppow_h([0,1,0,0,0,0],C,C);AM=A._ppow_h(An,C,C);AN=A._ppow_h(AM,C,C);Ao=A._ppow_h(AN,C**3,C);Ap=[M*A%C for A in O];J=Ap;l=[1,0,0,0,0,0]
	for I in range(6):l=A._pmul_h(l,J,C);J=A._ppow_h(J,C,C)
	m=C==1+6*M and Ao==[0,1,0,0,0,0]and AM!=[0,1,0,0,0,0]and AN!=[0,1,0,0,0,0]and pow(l[0],(C-1)//2,C)==C-1;assert m;AO=[F(A,3)for A in A.SIGMA_NUM];n=[[F(1)]+[F(0)]*5,AO]
	for I in range(4):n.append(A._pmul_h_q(n[-1],AO))
	AP=tuple(tuple(int(3*n[B][A])for B in range(6))for A in range(6));Aq=[[1 if A==B else 0 for B in range(6)]for A in range(6)];AQ=A._isp;AR=math.sumprod
	def P(n,p):
		n%=p
		if n==0:return 0
		if pow(n,(p-1)//2,p)!=1:return-1
		if p%4==3:A=pow(n,(p+1)//4,p);return min(A,p-A)
		B,H=p-1,0
		while B%2==0:B//=2;H+=1
		F=2
		while pow(F,(p-1)//2,p)!=p-1:F+=1
		I,J,C,A=H,pow(F,B,p),pow(n,B,p),pow(n,(B+1)//2,p)
		while C!=1:
			G,E=1,C*C%p
			while E!=1:E=E*E%p;G+=1
			D=pow(J,1<<I-G-1,p);I,J,C,A=G,D*D%p,C*D*D%p,A*D%p
		return min(A,p-A)
	D=A._mh_int
	def V(p):C=tuple(AR(A,p)for A in AP);F=tuple(AR(A,C)for A in AP);E=D(D(p,C),F);A=-(E[2]//729);B=E[0]//729-2*A;return B,A,B*B+2*A*A
	o={}
	def G(p):
		B=o.get(p)
		if B is not None:return B
		if p<=3 or p&7 not in(1,3)or A._ppow_h([0,1,0,0,0,0],p,p)!=[0,1,0,0,0,0]:o[p]=None;return
		C=P(-2,p);D=p-C;E=A._roots_mod(A.H_COEFFS,p);F=next(A for A in E if(2+j._eval_poly_mod(O,A,p))%p==C);G=next(A for A in E if(2+j._eval_poly_mod(O,A,p))%p==D);B=C,D,(F,G),(P(2+C,p),P(2+D,p)),(P(pow(M*(C-2)%p,-1,p),p),P(pow(M*(D-2)%p,-1,p),p));o[p]=B;return B
	AS={(int(A),int(B))for(A,B)in B['ret_f_split']};AT={(int(A),int(B))for(A,B)in B['ret_g_split']};Ar={(int(A),int(B))for(A,B)in B['ret_f_ns']};As={(int(A),int(B))for(A,B)in B['ret_g_ns']};T=sorted({A for(A,B)in AS}|{A for(A,B)in AT});AU=A._primes_upto(AK);p=set(T);q=[];AV={}
	for E in AU:
		if E<=AJ and E not in p:
			if G(E)is None:q.append(E);AV[E]=P(-2,E)if E>3 and E&7 in(1,3)else-1
			else:T.append(E);p.add(E)
	for E in T:G(E)
	def r(xs):
		A=xs;A=list(A)
		while len(A)>1:A=[A[B]*A[B+1]if B+1<len(A)else A[B]for B in range(0,len(A),2)]
		return A[0]
	At,Au,Av=r(q),r(T),r([A for A in AU if A>AJ and A not in p]);s=[];J=[0]*6
	def AW(idx,rem):
		A=idx
		if A==6:
			if any(J):s.append(tuple(J))
			return
		C=math.isqrt(rem)
		for B in range(-C,C+1):J[A]=B;AW(A+1,rem-B*B)
	AW(0,int(B['ball_r2']));assert len(s)==int(B['ball_pts'])==484;Aw=[tuple((B,A)for(B,A)in enumerate(A)if A)for A in s]
	def AX(n,c0,c1,sd):
		B=n;F=Ar if sd==0 else As;H=AS if sd==0 else AT;I=math.gcd(B,At)
		if I>1:
			for A in q:
				if I%A==0:
					while B%A==0:B//=A
					C=AV[A]
					if C>=0:
						D=(c0+c1*C)%A==0;E=(c0+c1*(A-C))%A==0
						if not D and not E or D and(A,0)not in F or E and(A,1)not in F:return
					elif(A,0)not in F:return
		J=math.gcd(B,Au)
		if J>1:
			for A in T:
				if J%A==0:
					while B%A==0:B//=A
					C,K=G(A)[0],G(A)[1];D=(c0+c1*C)%A==0;E=(c0+c1*K)%A==0
					if not D and not E or D and(A,0)not in H or E and(A,1)not in H:return
		if math.gcd(n,Av)>1:return
		return B
	def AY(cof,nw,km):
		I=cof
		if I==1:return[]
		if km<=0 or I.bit_length()>nw*km:return
		L=1<<nw;F=[];J=[I]
		while J:
			A=J.pop()
			if A==1:continue
			if AQ(A):
				if A>L or A<=AK or G(A)is None:return
				if A not in F:
					if len(F)>=km:return
					F.append(A)
			else:
				H=None
				for E in range(1,32):
					C=B=2;D=1
					while D==1:
						K=1;M,N=C,B
						for O in range(16):C=(C*C+E)%A;B=(B*B+E)%A;B=(B*B+E)%A;K=K*(C-B)%A
						D=math.gcd(K,A)
						if D==A:
							C,B,D=M,N,1
							while D==1:C=(C*C+E)%A;B=(B*B+E)%A;B=(B*B+E)%A;D=math.gcd(abs(C-B),A)
					if 1<D<A:H=D;break
				if H is None:return
				J.extend((H,A//H))
		return F
	W={}
	def t(sd,q,oi,nw,km):
		H=sd;J=H,q,oi,nw,km
		if J in W:return W[J]
		AV,AW,A8,A9,AA=G(q);AB=A8[oi];AC=A9[oi]if H==0 else AA[oi];k=[j._eval_zk_mod([1 if B==A else 0 for B in range(6)],AB,q)for A in range(6)];x=pow(k[0],-1,q);F=[[0]*12 for A in range(12)];F[0][0]=q
		for E in range(1,6):F[E][E]=1;F[E][0]=-k[E]*x%q
		for E in range(6):F[6+E][6+E]=1;F[6+E][0]=AC*k[E]*x%q
		Ag.bareiss_lll_12d(F,Aq,1,1,False);y=tuple(sum(F[0][A]*N[A][B]for A in range(6))for B in range(6));z=tuple(sum(F[0][6+A]*N[A][B]for A in range(6))for B in range(6));l=A0=None
		for A1 in range(1,12):
			A2=tuple(sum(F[A1][A]*N[A][B]for A in range(6))for B in range(6));A3=tuple(sum(F[A1][6+A]*N[A][B]for A in range(6))for B in range(6))
			if D(y,A3)!=D(z,A2):l,A0=A2,A3;break
		if l is None:W[J]=[];return[]
		AD=[tuple(A//3 for A in D(N[A],l))for A in range(6)];AE=[tuple(A//3 for A in D(N[A],A0))for A in range(6)];m=[]
		for AF in Aw:
			K,L,P,Q,R,S=y;T,U,X,Y,Z,a=z
			for(E,C)in AF:
				A,B=AD[E],AE[E]
				if C==1:K+=A[0];L+=A[1];P+=A[2];Q+=A[3];R+=A[4];S+=A[5];T+=B[0];U+=B[1];X+=B[2];Y+=B[3];Z+=B[4];a+=B[5]
				elif C==-1:K-=A[0];L-=A[1];P-=A[2];Q-=A[3];R-=A[4];S-=A[5];T-=B[0];U-=B[1];X-=B[2];Y-=B[3];Z-=B[4];a-=B[5]
				else:K+=C*A[0];L+=C*A[1];P+=C*A[2];Q+=C*A[3];R+=C*A[4];S+=C*A[5];T+=C*B[0];U+=C*B[1];X+=C*B[2];Y+=C*B[3];Z+=C*B[4];a+=C*B[5]
			n=K,L,P,Q,R,S;o=T,U,X,Y,Z,a
			if not any(n)or not any(o):continue
			b=tuple(A//3 for A in D(n,n));c=tuple(A//3 for A in D(o,o))
			if H==0:p=D(AL,c);d,e,f=V(tuple(b[A]+p[A]for A in range(6)))
			else:r=D(O,b);d,e,f=V(tuple(M*r[A]-c[A]for A in range(6)))
			if f%q!=0 or f==0:continue
			A4=AX(f//q,d,e,H);g=AY(A4,nw,km)if A4 is not None else None
			if g is None:continue
			if H==0:r=D(O,b);h,i,s=V(tuple(M*r[A]-c[A]for A in range(6)))
			else:p=D(AL,c);h,i,s=V(tuple(b[A]+p[A]for A in range(6)))
			A5=AX(s,h,i,1-H)if s else None;t=AY(A5,nw,km-len(g))if A5 is not None else None
			if t is None:continue
			(AG,AH,AI),(AJ,AK,AM)=((g,d,e),(t,h,i))if H==0 else((t,h,i),(g,d,e));u=[];I=True
			for(A6,AN,AO,AP)in((0,AG,AH,AI),(1,AJ,AK,AM)):
				for v in AN:
					AQ,AR,_,AS,AT=G(v);A7=0
					for(w,AU)in((0,AQ),(1,AR)):
						if(AO+AP*AU)%v==0:
							if(AS[w]if A6==0 else AT[w])<=0 or len(u)>=km:I=False;break
							u.append((A6,v,w));A7+=1
					if not I or A7==0:I=False;break
				if not I:break
			if I:m.append(u);break
		W[J]=m;return m
	def AZ(sd,q,oi,nw,km,dsh,rec=False):
		F=t(sd,q,oi,nw,km);A=1
		if not F:return False,A,([],0,0)
		B=F[0]
		if not B:return True,A,([[0,0,0,0,0]]if rec else[],0,0)
		C=True;D=0;E=0;G=0
		for(K,H,L)in B:
			E+=1;M=H.bit_length();I=max(S,M-dsh);N=0 if I==S else 1;J,O,(T,P,Q)=AZ(K,H,L,I,N,dsh,False);D+=O;E+=P;G+=(0 if J else 1)+Q
			if not J:C=False;break
		A+=D;R=[[0,int(not C),D,len(B),0]]if rec else[];return C,A,(R,E,G)
	u=0
	for(H,v,w,x,Ax,Ay,Az,Aa)in B['live_witnesses']:Ab,I,(X,I,I)=AZ(v,w,x,Ax,Ay,Az,True);assert X and X[0][0]==0 and X[0][1]==Aa[0][1]and X[0][3]==Aa[0][3];u+=1
	K=[];U=(1<<16)+1
	while len(K)<64:
		if AQ(U):
			y=G(U)
			if y is not None:
				for Y in(0,1):
					if y[3][Y]>0:K.append((0,U,Y))
					if y[4][Y]>0:K.append((1,U,Y))
		U+=2
	K=K[:64];Z=a=z=Ac=b=0;A0=.0
	def Ad(sd,q,oi,nw,km,dsh):
		K=t(sd,q,oi,nw,km);A=1
		if not K:return False,False,A,0,0,0,.0
		L=K[0]
		if not L:return True,True,A,0,0,0,.0
		B=True;C=0;D=0;E=0;F=0;G=.0
		for(H,I,M)in L:
			D+=1;N=I.bit_length();J=max(S,N-dsh);O=0 if J==S else 1;G+=d(0,H,N)[0];Q=t(H,I,M,J,O)
			if not Q:C+=1;F+=1;E+=1;B=False;break
			P,X,R,T,U,V,W=Ad(H,I,M,J,O,dsh);C+=R;D+=T;E+=(0 if P else 1)+U;F+=V;G+=W
			if not P:B=False;break
		A+=C;return B,True,A,D,E,F,G
	for(v,w,x)in K:
		Ab,A_,B4,B0,B1,B2,B3=Ad(v,w,x,16,2,1)
		if Ab:a+=1
		if A_:Z+=1
		z+=B0;Ac+=B1;b+=B2;A0+=B3
	Ae=b/max(A0,1e-12);Af=Z==a+b and Z>=25 and a>=18 and z>=15 and Ae<=.85;return dict(ok=u==3 and Af and m and A1['tr_valid'],strict_k1_zero_retry_protocol=True,live_k1_zero_retry_ok=Af,k1_17b_roots=len(K),k1_17b_nonempty=Z,k1_17b_completed=a,k1_17b_reached_children=z,k1_17b_slot_zeros=b,k1_17b_exp_p0_sum=round(A0,3),k1_17b_direct_slot_oe_p0=round(Ae,4),k1_17b_child_fails=Ac,threshold_counterexample=A1,tower_u2_p13_verified=m,roots_total=sum(A[2]for A in Q),first_attempts_total=sum(A[0]for A in Q),after_fail_total=e,after_fail_failures=A7,pooled_oe_fail=round(k,4),boole_union_attempt_oe=round(k,4),pooled_oe_work=round(AI,4),shared_prev_total=A8,shared_run_total=A9,conc_bins_checked=len(L),max_C_emp=max(A[5]for A in L),max_C_bday=max(A[6]for A in L),live_witnesses_replayed=u,per_setting=A2)
def verify_h3_scaling_and_nu_h():
	import model_bn254_v6 as A;G=A.load_data()['filtsim_raw'];H='m184_z6_61.50_61.50_e0.72_W100_';D=[A for A in(20,21,22,23,24,25,26,27,28)if f"{H}{A}"in G];I=[];U=[];V=[];s=[];f=True;g=True;W=10**9
	for AY in D:
		N=G[f"{H}{AY}"];X=len(N);I.append(sum(A[1]for A in N)/X);U.append(sum(A[2]for A in N)/X);Y=[]
		for C in N:
			if len(C)>=12:
				t,u,v,w,x,y=C[6],C[7],C[8],C[9],C[10],C[11]
				if abs(t+(u+v+w+x)/y-1.)>=.0005:f=False
				Y.append(t*A.C_UNIT_5+(u*2.*A.C_UNIT_5+v*3.*A.C_UNIT_5+w*4.*A.C_UNIT_5+x*A.C_GEN_1)/y)
				if len(C)>=14 and(C[12]!=0 or C[13]!=0):g=False
				if len(C)>=15:W=min(W,C[14])
			elif len(C)==8:Y.append(C[7])
			else:Y.append(C[6]*A.C_UNIT_5+(1.-C[6])*A.C_GEN_1)
		V.append(sum(Y)/X);s.append(sum(A[4]for A in N)/X)
	h=all(A[B]>A[B+1]for A in(I,U,V)for B in range(len(D)-1));Z=[math.log(A)for A in D];i=[math.log(A)for A in I];j,z=sum(Z)/len(Z),sum(i)/len(i);O=-sum((A-j)*(B-z)for(A,B)in zip(Z,i))/sum((A-j)**2 for A in Z);k=math.exp(z+O*j);A0=max(abs(k/B**O-A)/A for(B,A)in zip(D,I));A1=k/A.log2_U_z6(A.OP['lf'],A.OP['lg'])**O;l=A.log2_U_z6(A.OP['lf'],A.OP['lg']);J=len(D);E=[sum(A[14]for A in G[f"{H}{A}"])/len(G[f"{H}{A}"])/2.**A for A in D];AZ=[sum(A[1]*2.**B-A[14]for A in G[f"{H}{B}"])/len(G[f"{H}{B}"])/math.floor(A.TEXC_REL*2.**B)for B in D];m=A.NC27;n=[math.log(A)for A in E]if min(E)>0 else[.0]*J;K=sum(D)/J;a=sum(n)/J;A2=math.exp(a+sum((A-K)*(B-a)for(A,B)in zip(D,n))/sum((A-K)**2 for A in D)*(l-K));o=[math.log(A)for A in D];p=sum(o)/J;A3=-sum((A-p)*(B-a)for(A,B)in zip(o,n))/sum((A-p)**2 for A in o);A4=math.exp(a+A3*p)/l**A3;A5=[1./A for A in E];A6=sum(A5)/J;Aa=sum((A-K)*(B-A6)for(A,B)in zip(D,A5))/sum((A-K)**2 for A in D);A7=A6+Aa*(l-K);A8=1./A7 if A7>0 else 1.;q=dict(ncols_over_U=[round(A,8)for A in E],surplus_over_texc=[round(A,4)for A in AZ],max_ratio=round(max(E[A+1]/E[A]for A in range(J-1)),4),pred_ypow=A4,pred_exp=A2,pred_invlin=A8,cap=m,cap_over_max_pred=round(m/max(A4,A2,A8),3),cap_over_c28=round(m/E[-1],3));A9=min(E)>0 and q['max_ratio']<1. and E[-2]<=A.NC27 and q['cap_over_max_pred']>=1.;AA=A._primes_upto(16384);b=P=Q=0;AB=AC=AD=0;R=c=AE=AF=AG=.0;Ab=[A.F_REC_COEFFS[B]for B in(0,2,4,6,8)]
	for B in AA:
		if B in(2,3,19):continue
		Q+=1;F=math.log2(B);R+=F;d=A._ppow_h([0,1,0,0,0,0],B,B);d[1]=(d[1]-1)%B
		if not any(d):P+=1;c+=F;r=A.H_COEFFS
		else:
			r=A._pgcd(d,A.H_COEFFS,B);e=len(r)-1
			if e==3:b+=1;c+=F
			elif e==6:P+=1;c+=F
			elif e!=0:raise AssertionError(f"d={e} p={B}")
			else:continue
		S=T=False
		for C in A._roots_mod(r,B):
			Ac=sum(A.W_COEFFS[D]*pow(C,D,B)for D in range(6))%B;AH=-((A.U_BN254+1)%B*Ac+2)%B
			if not T and(AH==0 or pow(AH,(B-1)//2,B)==1):T=True
			if not S:
				L=[sum(D[A]*pow(C,A,B)for A in range(6))%B for D in Ab];AI=[L[0],0,L[1],0,L[2],0,L[3],0,L[4]];M=A._ppowmod([0,1],B,AI,B)
				if len(M)<2:M+=[0]*(2-len(M))
				M[1]=(M[1]-1)%B
				if len(A._pgcd(M,AI,B))>1 or L[0]==0:S=True
			if S and T:break
		if S:AB+=1;AE+=F
		if T:AC+=1;AF+=F
		if S or T:AD+=1;AG+=F
	AJ=3,13,29,983,11003,405928799;AK=1
	for Ad in AJ:AK*=Ad
	AL=AK==A.U_BN254 and all(A._isp(B)and B<1<<30 for B in AJ);AM=b+2*P;AN=AM/(len(AA)/3.);AO=(b+P)/float(Q);AP=c/R;AQ=AB/float(Q);AR=AE/R;AS=AC/float(Q);AT=AF/R;AU=AD/float(Q);AV=AG/R;AW=AQ<=A.KH_F_CAP and AR<=A.KH_F_CAP and AS<=A.KH_G_CAP and AT<=A.KH_G_CAP and AU<=A.KH_FG_CAP and AV<=A.KH_FG_CAP;AX=1.15;Ae=D==[20,21,22,23,24,25,26,27,28]and W>=1000 and h and f and g and O>2.5 and A0<.15 and A1<A.NU_SW_CAP/1e1 and A.NU_SW_CAP/I[-1]>=AX and A.NU_L_CAP/U[-1]>=AX and V[-1]<A.K5_CAP and s[-1]<A.TAB_R_CAP and .97<=AN<=1.03 and AO<=A.KH_SPLIT_CAP and AP<=A.KH_SPLIT_CAP and AW and AL and A9;return dict(ok=Ae,fixed_instance_deterministic_h3=True,zero_loss_filtsim_ok=g,min_ncols_all44=W,scales=D,all4_monotone_20_28=h,all4_monotone_20_27=h,hist_v4_exact_ok=f,inv_log_A=round(k,3),inv_log_beta=round(O,3),max_rel_fit_err=round(A0,4),pred_nsw_prod=round(A1,6),nu_sw_cap=A.NU_SW_CAP,cap_over_top_scale=round(A.NU_SW_CAP/I[-1],3),nu_l_cap_over_top_scale=round(A.NU_L_CAP/U[-1],3),k5_v8_top_scale=round(V[-1],4),nu_h_upto_2p14=round(AN,4),kh_split_frac_2p14=round(AO,4),kh_log_wt_frac_2p14=round(AP,4),f_rec_split_2p14=round(AQ,4),f_rec_log_2p14=round(AR,4),g_rec_split_2p14=round(AS,4),g_rec_log_2p14=round(AT,4),fg_split_2p14=round(AU,4),fg_log_2p14=round(AV,4),fg_dens_ok=AW,u_ram_all_lt_2p30=AL,c6_d3_primes_2p14=b,c6_d6_primes_2p14=P,c6_deg1_orbits_2p14=AM,h3_scale=q,h3_scale_ok=A9)
def verify_split_crt_norms_and_mp_rt(seed=20261120):
	import random;from fractions import Fraction as O;import model_bn254_v6 as B;Bp=143,495,663,1701,3023,3551,4485,4557,5055;Aa=[(1<<62)-A for A in Bp];AG=[]
	for A in Aa:
		AH=B._roots_mod(B.H_COEFFS,A)
		if len(AH)!=6:return dict(ok=False,err='not_split')
		Bq=pow(3,-1,A);h=[[sum(B.ZK_NUM[E][C]*pow(D,C,A)for C in range(6))*Bq%A for E in range(6)]for D in AH];y=[((B.U_BN254+1)*sum(B.W_COEFFS[C]*pow(D,C,A)for C in range(6))+2)%A for D in AH];AI=[[sum(B.F_REC_COEFFS[E][C]*pow(D,C,A)for C in range(6))%A for E in(0,2,4,6,8)]for D in AH];AG.append((A,h,y,AI))
	BE=[[O(A,3)for A in A]for A in B.ZK_NUM]
	def P(a,b):
		A=[O(0)]*11
		for E in range(6):
			for C in range(6):A[E+C]+=a[E]*b[C]
		for D in range(10,5,-1):
			F=A[D]
			if F:
				for C in range(6):A[D-6+C]-=F*B.H_COEFFS[C]
				A[D]=O(0)
		return A[:6]
	def BF(p_vec):
		E=[p_vec]
		for I in range(5):E.append(P(E[-1],[O(0),O(1),O(0),O(0),O(0),O(0)]))
		B=[[E[B][A]for B in range(6)]for A in range(6)];C=O(1)
		for A in range(6):
			D=next((C for C in range(A,6)if B[C][A]!=0),None)
			if D is None:return 0
			if D!=A:B[A],B[D]=B[D],B[A];C=-C
			C*=B[A][A]
			for F in range(A+1,6):
				H=B[F][A]/B[A][A]
				for G in range(A,6):B[F][G]-=H*B[A][G]
		assert C.denominator==1;return abs(int(C.numerator))
	Br=[O(A)for A in B.W_COEFFS];Bs={A:[O(A)for A in B.F_REC_COEFFS[A]]for A in(0,2,4,6,8)};K=random.Random(seed);Ab=0;Ac=Ad=0
	for _ in range(20):
		AJ=[K.randint(-400,400)for A in range(6)];AK=[K.randint(-500,500)for A in range(6)];z=[];BG=[]
		for(A,h,y,AI)in AG:
			i=Ae=1
			for I in range(6):A0=sum(AJ[A]*h[I][A]for A in range(6))%A;X=sum(AK[A]*h[I][A]for A in range(6))%A;V=A0*A0%A;Y=X*X%A;i=i*(y[I]*V+Y)%A;Z=V*V%A;a=Y*Y%A;A1=V*Y%A;A2=AI[I];Ae=Ae*(A2[0]*a*a+A2[1]*A1*a+A2[2]*Z*a+A2[3]*Z*A1+A2[4]*Z*Z)%A
			z.append(i);BG.append(Ae)
		def Af(rs,cnt):
			A=0;B=1
			for((C,*F),D)in zip(AG[:cnt],rs[:cnt]):E=(D-A)%C*pow(B%C,-1,C)%C;A+=E*B;B*=C
			return abs(A-B if A>B//2 else A)
		Bt=Af(z,9);Bu=Af(BG,9);BH=[sum(AJ[A]*BE[A][B]for A in range(6))for B in range(6)];BI=[sum(AK[A]*BE[A][B]for A in range(6))for B in range(6)];V=P(BH,BH);Y=P(BI,BI);Bv=[(B.U_BN254+1)*A+2*C+D for(A,C,D)in zip(P(Br,V),V,Y)];Z=P(V,V);a=P(Y,Y);A1=P(V,Y);Bw=[P(a,a),P(A1,a),P(Z,a),P(Z,A1),P(Z,Z)];Ag=[O(0)]*6
		for(AL,Bx)in enumerate((0,2,4,6,8)):By=P(Bs[Bx],Bw[AL]);Ag=[Ag[A]+By[A]for A in range(6)]
		BJ=BF(Bv);BK=BF(Ag)
		if Bt!=BJ or Bu!=BK:Ab+=1
		Ac=max(Ac,BJ.bit_length());Ad=max(Ad,BK.bit_length())
	L=1<<64;AM=0;Ah=0
	def Ai(R_int,P_int,m,N):
		F=[R_int>>64*A&L-1 for A in range(N)];G=[P_int>>64*A&L-1 for A in range(m)];B=[0]*N
		for C in range(N):
			D=F[C]
			if D:
				for E in range(m):B[(C+E)%N]+=D*G[E]
		A=0
		for H in range(N-1,m-1,-1):A=(A<<64)+B[H]
		return A&(1<<64*(N-m))-1
	for j in(2,4,8):
		Bz=8;AN=2*j+Bz;BL=L**(AN-1)+(L**j-1)//(L+1);BM=L+1;B_=BL*BM//L**j%L**(AN-j);C0=Ai(BL,BM,j,AN)
		if not 0<=(B_-C0)%L**(AN-j)<j*L:AM+=1
	A3=4;s=2;AO=A3+4;C1=1<<64*(s+AO)-1
	for BN in(0,1):
		t=[L**s-1-K.randrange(1<<20)if BN else K.randrange(L,L**s)for A in range(1<<A3)];BO=L**40-1-K.randrange(1<<30)if BN else K.randrange(1,L**40);A4=[t[:]]
		for _ in range(A3):Aj=A4[-1];A4.append([Aj[2*A]*Aj[2*A+1]for A in range(len(Aj)//2)])
		BP=A4[-1][0];C2=s*(1<<A3);Ak=[BO%BP*L**(C2+AO)//BP]
		for I in range(A3,0,-1):
			Al=s*(1<<I-1);BQ=2*Al+AO;k=[]
			for(AL,AP)in enumerate(Ak):k.append(Ai(AP,A4[I-1][2*AL+1],Al,BQ));k.append(Ai(AP,A4[I-1][2*AL],Al,BQ))
			Ak=k
		for(AP,Am)in zip(Ak,t):
			if(AP*Am+C1>>64*(s+AO))%Am==BO%Am:Ah+=1
			else:AM+=1
	BR=round(sum(math.log2(A)for A in Aa),3);l=B.H_COEFFS;u=B.W_COEFFS;C3=B.U_BN254;C4=[B.F_REC_COEFFS[A]for A in(0,2,4,6,8)];An=Ao=.0
	for A in B._primes_upto(500):
		if A<5 or A==19:continue
		C5=[B for B in range(A)if((((((l[6]*B+l[5])*B+l[4])*B+l[3])*B+l[2])*B+l[1])*B+l[0])%A==0];AQ=Ap=0
		for W in C5:
			v=[(((((B[5]*W+B[4])*W+B[3])*W+B[2])*W+B[1])*W+B[0])%A for B in C4]
			for e in range(A):
				if((((v[4]*e+v[3])*e+v[2])*e+v[1])*e+v[0])%A==0:
					if e==0:AQ+=1
					elif pow(e,(A-1)//2,A)==1:AQ+=2
			if v[4]==0:AQ+=1
			BS=-((C3+1)*(((((u[5]*W+u[4])*W+u[3])*W+u[2])*W+u[1])*W+u[0])+2)%A
			if BS==0:Ap+=1
			elif pow(BS,(A-1)//2,A)==1:Ap+=2
		An+=(1./(A-1)-AQ*A/float(A*A-1))*math.log2(A);Ao+=(1./(A-1)-Ap*A/float(A*A-1))*math.log2(A)
	C6=[0 if A%3 else A+1 for A in range(64)];AR=0;Aq=0
	for(f,BT)in enumerate(C6):
		C7=1 if BT>=1 else 0;AR=AR<<1|C7
		if BT>0:Aq|=1<<63-f
	Ar=0;AS=0;N=47;As=(1<<N)-1;At=[(K.randrange(1,500),K.randint(12,25))for A in range(8)];Au=sum(A for(B,A)in At)
	for R in(8,16,32,64,128,256,512,1024):
		C=2003;BU=K.randrange(1,C);AT=(C+R-1)//R;C8=R-1<<N;Av=[0]*Au;BV=[0]*AT;b=0
		for(Aw,Ax)in At:
			Q=Aw%C
			for _ in range(Ax):
				Av[b]=Q;BV[Q//R]+=1;b+=1;Q+=BU
				if Q>=C:Q-=C
		A5=[0]*(AT+1)
		for w in range(AT):A5[w+1]=A5[w]+BV[w]
		BW=A5[:AT];BX=[0]*Au;A6=0;C9=R-1<<N
		for f in range(Au):
			Q=Av[f];w=Q//R;A6=A6&~C9|(Q&R-1)<<N;BX[BW[w]]=A6;BW[w]+=1;A6+=1
			if A6&As!=f+1:AS+=1
		Ay=[]
		for(Aw,Ax)in At:
			for e in range(Ax):Ay.append((Aw+e*BU)%C)
		if Av!=Ay:AS+=1
		for m in(20,390):
			for D in(0,R-3,C-15,K.randrange(C)):
				for Az in(1,-1):
					BY=set()
					for(f,Q)in enumerate(Ay):
						I=(Q-D)*Az%C
						if I<m:BY.add((I,f))
					c=set();n=[]
					if Az==1:
						if D+m<=C:n.append((0,D,D+m))
						else:n.append((0,D,C));n.append((C-D,0,D+m-C))
					elif D+1>=m:n.append((0,D-m+1,D+1))
					else:n.append((0,0,D+1));n.append((D+1,C-(m-D-1),C))
					for(A7,G,H)in n:
						o=G//R;p=(H-1)//R
						for F in range(o,p+1):
							S=F*R;g=S+R;BZ=S>=G and g<=H
							for J in BX[A5[F]:A5[F+1]]:
								if Az==1:
									T=J;AU=T>>N;d=T&As
									if BZ:c.add((A7+(S-G)+AU,d))
									elif o!=p and F==o:
										if T>=G-S<<N:c.add((A7+(S-G)+AU,d))
									elif o!=p and F==p:
										if T<H-S<<N:c.add((A7+(S-G)+AU,d))
									else:
										A_=max(G,S)-S<<N;B0=min(H,g)-max(G,S)<<N
										if 0<=T-A_<B0:c.add((A7+(S-G)+AU,d))
								else:
									A8=J^C8;CA=A8>>N;d=A8&As;AV=A7+(H-g)+CA
									if BZ:c.add((AV,d))
									elif o!=p and F==p:
										if A8>=g-H<<N:c.add((AV,d))
									elif o!=p and F==o:
										if A8<g-G<<N:c.add((AV,d))
									else:
										A_=g-min(H,g)<<N;B0=min(H,g)-max(G,S)<<N
										if 0<=A8-A_<B0:c.add((AV,d))
					if c!=BY:Ar+=1
	B1=0;Ba=0;AW=(1<<47)-1;Bb=(1<<64)-1
	for AX in(6,7,8,9):
		E=1<<AX;B2=E-1;CB=B2<<47
		for C in(E*3+17,E*5+3):
			CC=(C+E-1)//E;AY=180;Bc=[0]+[K.randrange(C)for A in range(AY)];B3=[0]*CC;k=[0]*(AY+1)
			for b in range(1,AY+1):Q=Bc[b];F=Q>>AX;J=(Q&B2)<<47|b;k[b]=B3[F];B3[F]=J
			for U in(19,E//2+5,E+11,min(C-1,2*E+37)):
				if U>=C:continue
				for B4 in(1,-1):
					for D in(0,E-3,C-1,C-U//2):
						Ba+=1;Bd=set()
						for b in range(1,AY+1):
							I=(Bc[b]-D)*B4%C
							if I<U:Bd.add((I,b))
						B5=set()
						if B4==1:Be=[(D,D+U-1,0)]if D+U<=C else[(D,C-1,0),(0,D+U-1-C,C-D)]
						else:Be=[(D-U+1,D,U-1)]if D>=U-1 else[(0,D,D),(C-(U-1-D),C-1,U-1)]
						for(q,A9,Bf)in Be:
							AA=q>>AX;AB=A9>>AX
							for F in range(AA,AB+1):
								J=B3[F]
								if J==0:continue
								if B4==1:
									X=Bf+(F*E-q)
									if AA==AB:G=q-F*E;H=A9-F*E+1;M=3
									elif F==AA:G=q-F*E;H=E;M=1
									elif F==AB:G=0;H=A9-F*E+1;M=2
									else:G=0;H=E;M=0
									AC=G<<47;B6=H<<47;t=H-G<<47
									while J:
										if M==0 or M==1 and J>=AC or M==2 and J<B6 or M==3 and J-AC&Bb<t:B5.add((X+(J>>47),J&AW))
										J=k[J&AW]
								else:
									X=Bf+q-(F*E+B2)
									if AA==AB:G=E-(A9-F*E+1);H=E-(q-F*E);M=3
									elif F==AA:G=0;H=E-(q-F*E);M=2
									elif F==AB:G=E-(A9-F*E+1);H=E;M=1
									else:G=0;H=E;M=0
									AC=G<<47;B6=H<<47;t=H-G<<47
									while J:
										T=J^CB
										if M==0 or M==1 and T>=AC or M==2 and T<B6 or M==3 and T-AC&Bb<t:B5.add((X+(T>>47),T&AW))
										J=k[T&AW]
						if B5!=Bd:B1+=1
	AD=B._primes_upto(16383)[1:]
	if B._P17_G_CACHE is not None:AD=AD+B._P17_G_CACHE
	else:
		for A in B._primes_upto((1<<17)-1):
			if A>=16384:
				if B.U_BN254%A==0:AD.append(A);continue
				B7=B._ppow_h([0,1,0,0,0,0],A,A);B7[1]=(B7[1]-1)%A;Bg=B._pgcd(B7,B.H_COEFFS,A)
				if len(Bg)>1:
					for AE in B._roots_mod(Bg,A):
						Bh=B.U_BN254*sum(B.W_COEFFS[C]*pow(AE,C,A)for C in range(6))%A
						if Bh==0 or pow(Bh,(A-1)//2,A)==1:AD.append(A);break
	AF=[]
	for _ in range(200):
		AJ=[K.randint(-2**15,2**15)for A in range(6)];AK=[K.randint(-2**17,2**17)for A in range(6)];z=[]
		for(A,h,y,AI)in AG:
			i=1
			for I in range(6):A0=sum(AJ[A]*h[I][A]for A in range(6))%A;X=sum(AK[A]*h[I][A]for A in range(6))%A;i=i*(y[I]*A0*A0-X*X)%A
			z.append(i)
		r=Af(z,9)
		if r==0:continue
		r>>=(r&-r).bit_length()-1;B8=0
		for A in AD:
			if r%A==0:
				while r%A==0:r//=A;B8+=1
				B8+=1
		AF.append(B8)
	Bi=sum(AF)/len(AF);CD=max(AF);B9=0;CE=0x8000800080008000;CM=65535;BA=B.rand_below;BB=K.getrandbits
	for CF in range(48):
		x=[[BA(BB,32768)for A in range(128)]for A in range(8)]
		if CF%3==0:x[K.randrange(8)][K.randrange(128)]|=32768
		BC=[32768-(3000+BA(BB,600))for A in range(128)];Bj=[A+48 for A in BC];Bk=[BA(BB,9)==0 for A in range(8)];CG=lambda r,i:r[4*i]|r[4*i+1]<<16|r[4*i+2]<<32|r[4*i+3]<<48;Bl=0;BD=0
		for AE in x:
			for f in range(32):Bl|=CG(AE,f)
		BD+=128+255+2;Bm=set()
		if Bl&CE:
			for(AZ,AE)in enumerate(x):
				for(CH,CI)in enumerate(AE):
					if CI&32768:Bm.add((AZ,CH))
		CJ={(A,B)for A in range(8)for B in range(128)if x[A][B]>=32768}
		for AZ in range(8):x[AZ]=list(Bj if Bk[AZ]else BC)
		BD+=16*(2+8)+4
		if Bm!=CJ or BD!=549 or any(x[A]!=(Bj if Bk[A]else BC)for A in range(8)):B9+=1
	CK=3.<An<5. and 1.5<Ao<3. and AR==Aq and Ar==0 and AS==0 and B1==0 and B9==0 and Bi<=16.;Bn=9*(6*122+70)+904;Bo=9*(6*280+70)+904;CL=Ab==0 and AM==0 and Ah==32 and CK and BR>557. and Bn<=B.C_NORM_G and Bo<=B.C_NORM_F;return dict(ok=CL,split_primes_9=len(Aa),prod_bits_9=BR,crt_norm_samples=20,crt_norm_errs=Ab,max_nf_bits=Ad,max_ng_bits=Ac,op_ng_wc=Bn,op_nf_wc=Bo,cyclic_mp_rt_errs=AM,srt_5trans_leaves_exact=Ah,murphy_af_unram_500=round(An,4),murphy_ag_unram_500=round(Ao,4),subs_adc_occbm_ok=AR==Aq,flat_b2_subband_errs=Ar,lex_tmp_kta_errs=AS,chained_1pass_b2_tests=Ba,chained_1pass_b2_errs=B1,scan8_blocks=48,scan8_errs=B9,scan8_ops_per_1024=549,strip_att_n=len(AF),strip_att_mean=round(Bi,3),strip_att_max=CD)
def verify_live_genuine_c6_relation_filter(seed=20261120):
	import random as Ay;import model_bn254_v6 as P,verify_bn254_v6 as AP;Aa,B=13,28393;AA=pow(Aa,10,B);Q=B*B;assert P._isp(B)and B==Aa**4-Aa**2+1 and(AA*AA-AA+1)%B==0;AB=P.H_COEFFS;Az=[A*AB[A]for A in range(1,7)];U,z=AP._eval_poly_mod,AP._eval_zk_mod;A0=[[[0]*6 for A in range(6)]for A in range(6)]
	for q in range(6):A0[0][q][q]=1
	A_=0
	for J in range(1,6):
		A0[J][0]=A0[0][J]
		for q in range(J,6):A0[J][q]=A0[q][J]=P._OH_U[A_];A_+=1
	Ab=math.sumprod;B0=[[[A0[C][A][B]for C in range(6)]for B in range(6)]for A in range(6)]
	def AQ(u,v):A=[[Ab(u,B0[A][B])for B in range(6)]for A in range(6)];return[Ab(v,[A[C][B]for C in range(6)])for B in range(6)]
	def B1(v):return abs(AP._det_int([[Ab(v,B0[A][B])for B in range(6)]for A in range(6)]))
	Bo=[-1,-1,-1,-1,1,-1];Bp=[2*A for A in Bo];Bq=[3,-1,-1,-1,1,-1];B2=(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1);AR={};AS={};Ac=set();AC={}
	def B3(k):
		if k not in AC:AC[k]=len(AC)
		return AC[k]
	for A in range(5,601):
		if not AP._is_prime(A)or A in(13,19):continue
		N=[B for B in range(A)if U(AB,B,A)==0]
		if len(N)in(3,6):
			B4=pow(3,-1,A);B5=set();AR[A]=sorted(N)
			for O in sorted(N):
				if O in B5:continue
				V=U(P.SIGMA_NUM,O,A)*B4%A;d=U(P.SIGMA_NUM,V,A)*B4%A
				if len({O,V,d})!=3:continue
				B5.update((O,V,d));B6=U(P.W_COEFFS,O,A);Ad=2*B6%A
				if Ad and pow(Ad,(A-1)//2,A)==1:
					Br=pow(Ad,-1,A);I=next(B for B in range(1,A)if B*B%A==Br);A1=B3(('tau6',0,A,O,I))
					for(AD,(A2,Ae))in enumerate(((O,I),(d,A-I),(V,I),(O,A-I),(d,I),(V,A-I))):AS[0,A,A2,Ae]=A1,AD
				Af=(B6+4)%A
				if Af and pow(Af,(A-1)//2,A)==1:
					I=next(B for B in range(1,A)if B*B%A==Af);A1=B3(('tau6',1,A,O,I))
					for(AD,(A2,Ae))in enumerate(((O,I),(d,A-I),(V,I),(O,A-I),(d,I),(V,A-I))):AS[1,A,A2,Ae]=A1,AD
		elif len(N)==0:
			if P._ppow_h([0,1,0,0,0,0],A**3,A)==[0,1,0,0,0,0]:Ac.add(A)
	B7=sorted(list(AR.keys())+list(Ac));Bs=set(B7)
	def B8(n):
		while n%2==0:n//=2
		while n%3==0:n//=3
		B={}
		for A in B7:
			if A*A>n:break
			if n%A==0:
				C=0
				while n%A==0:n//=A;C+=1
				B[A]=C
		if n>1:
			if n in Bs:B[n]=B.get(n,0)+1
			else:return
		return B
	Bt=P._roots_mod(AB,B);B9=pow(3,-1,B)
	def Bu(s):return(s-U(AB,s,Q)*pow(U(Az,s,Q),-1,Q))%Q
	BA=[];BB=set()
	for O in Bt:
		if O in BB:continue
		V=U(P.SIGMA_NUM,O,B)*B9%B;d=U(P.SIGMA_NUM,V,B)*B9%B;BB.update((O,V,d));BA.append([Bu(A)for A in(O,V,d)])
	def BC(u,v,c,md):return(u[0]*v[0]+c*u[1]*v[1])%md,(u[0]*v[1]+u[1]*v[0])%md
	def Bv(u,e,c,md):
		B=1,0;A=u
		while e>0:
			if e&1:B=BC(B,A,c,md)
			A=BC(A,A,c,md);e>>=1
		return B
	AT=[]
	for AE in BA:
		BD=U(P.W_COEFFS,AE[0],Q);BE=pow(2*BD%Q,-1,Q);BF=(BD+4)%Q;Bw=B-1 if pow(BE%B,(B-1)//2,B)==1 else Q-1;Bx=B-1 if pow(BF%B,(B-1)//2,B)==1 else Q-1
		if not AT:AT.append((AE,BE,Bw))
		AT.append((AE,BF,Bx))
	r=[]
	for J in range(6):
		K=[0]*6;K[J]=1;r.append(tuple(K))
		for q in range(J+1,6):
			for By in(1,-1):K=[0]*6;K[J]=1;K[q]=By;r.append(tuple(K))
	Bz=Ay.Random(seed)
	for _ in range(300):
		K=tuple(Bz.choice((-1,0,0,0,1))for A in range(6))
		if any(K):r.append(K)
	r=list(dict.fromkeys(r));H=[];Ag=[];BG=set();BH=0;B_=pow(AA,2,B)
	for A3 in r:
		i=AQ(A3,A3);C0=AQ(Bp,i)
		for A4 in r:
			if A3==A4:continue
			a=AQ(A4,A4);BI=B1([C0[A]-a[A]for A in range(6)])
			if BI<=1:continue
			BJ=B8(BI)
			if BJ is None:continue
			C1=AQ(Bq,a);BK=B1([i[A]-C1[A]for A in range(6)])
			if BK<=1:continue
			BL=B8(BK)
			if BL is None:continue
			s={};AF=True;BM=0
			for(AU,C2)in((0,BJ),(1,BL)):
				BN=1 if AU==0 else-1
				for(A,A5)in C2.items():
					if A in Ac:continue
					BO=0
					for AG in AR[A]:
						BP=z(A3,AG,A);BQ=z(A4,AG,A)
						if BQ==0:
							if BP==0:AF=False;break
							continue
						I=BP*pow(BQ,-1,A)%A
						if(AU,A,AG,I)in AS:
							K=1
							if A5>1:
								T=A**A5;t=AG;A6=I
								for _ in range(A5):t=(t-U(AB,t,T)*pow(U(Az,t,T),-1,T))%T
								BR=U(P.W_COEFFS,t,T);C3=pow(2*BR%T,-1,T)if AU==0 else(BR+4)%T
								for _ in range(A5):A6=(A6-(A6*A6-C3)*pow(2*A6,-1,T))%T
								BS=(z(A3,t,T)-z(A4,t,T)*A6)%T;K=0
								while K<A5 and BS%A==0:K+=1;BS//=A
							A1,AD=AS[AU,A,AG,I];C4,C5=B2[AD];C6,C7=s.get(A1,(0,0));s[A1]=C6+BN*K*C4,C7+BN*K*C5;BO+=K;BM+=K
					if not AF:break
					if len(AR[A])==6 and BO!=A5:AF=False;break
				if not AF:break
			if not AF:continue
			s={B:A for(B,A)in s.items()if A!=(0,0)}
			if len(s)>=2:
				BT=tuple(sorted(s.items()))
				if BT not in BG:
					BG.add(BT);H.append(s);BH+=BM;BU=[]
					for(AE,C8,C9)in AT:
						A7=0;AH=1
						for A2 in AE:CA=Bv((z(A3,A2,Q),-z(A4,A2,Q)%Q),C9,C8,Q);A7=(A7+AH*(CA[1]//B))%B;AH=AH*B_%B
						BU.append(A7)
					Ag.append(BU)
			if len(H)>=160:break
		if len(H)>=160:break
	k=lambda a,b:(a+b*AA)%B;AI=sorted({B for A in H for B in A});CB={B:A for(A,B)in enumerate(AI)};u=len(AI);AV=u+3;L=[]
	for(W,CC)in zip(H,Ag):
		Ah=[0]*AV
		for(C,CD)in W.items():Ah[CB[C]]=k(*CD)
		for(CE,CF)in enumerate(CC):Ah[u+CE]=CF
		L.append(Ah)
	e=0;Ai=[]
	for AJ in range(AV):
		Aj=next((A for A in range(e,len(L))if L[A][AJ]!=0),None)
		if Aj is None:continue
		L[e],L[Aj]=L[Aj],L[e];CG=pow(L[e][AJ],-1,B);L[e]=[A*CG%B for A in L[e]]
		for J in range(len(L)):
			if J!=e and L[J][AJ]:CH=L[J][AJ];L[J]=[(A-CH*C)%B for(A,C)in zip(L[J],L[e])]
		Ai.append(AJ);e+=1
	BV=[A for A in range(AV)if A not in Ai];Ak=len(BV);AW=[0]*AV
	if Ak==1:
		BW=BV[0];AW[BW]=1
		for(J,CI)in enumerate(Ai):AW[CI]=-L[J][BW]%B
	A8={AI[A]:AW[A]for A in range(u)};CJ=AW[u:];A9={A:-sum(Ag[A][B]*CJ[B]for B in range(3))%B for A in range(len(H))};BX=sum(1 for(A,C)in enumerate(H)if sum(k(*B)*A8[A]for(A,B)in C.items())%B!=A9[A]);Al=len(H);R={}
	for(E,W)in enumerate(H):
		for C in W:R.setdefault(C,set()).add(E)
	f=set(range(Al));v=[A for(A,B)in R.items()if len(B)==1]
	while v:
		C=v.pop();N=[A for A in R.get(C,())if A in f]
		if len(N)==1:
			E=N[0];f.remove(E)
			for D in H[E]:
				R[D].discard(E)
				if len(R[D])==1:v.append(D)
	AX={A for(A,B)in R.items()if len(B)>=2};Am=[];CK=min(len(f),int(1.2*len(AX))+4)
	for E in sorted(f,key=lambda i:len(H[i]),reverse=True):
		if len(f)<=CK:break
		if all(len(R[A])>=2 for A in H[E]):
			f.remove(E)
			for C in H[E]:R[C].discard(E)
			An=[A for A in H[E]if len(R[A])==1]
			while An:
				C=An.pop()
				if len(R[C])==1:
					F=next(iter(R[C]))
					if all(len(R[A])>=2 for A in H[F]if A!=C):
						f.remove(F);Am.append((C,F,dict(H[F])))
						for D in H[F]:
							R[D].discard(F)
							if len(R[D])==1:An.append(D)
	Ao={B for A in f for B in H[A]};BY={A for(A,B,B)in Am};BZ=len(AX-(Ao|BY));Ap={A:A8[A]for A in Ao};AY=0
	for(C,E,Ba)in reversed(Am):
		AZ=k(*Ba[C])
		if AZ==0:AY+=1;continue
		Aq=sum(k(*B)*Ap[A]for(A,B)in Ba.items()if A!=C)%B;Ap[C]=(A9[E]-Aq)*pow(AZ,-1,B)%B
		if Ap[C]!=A8[C]:AY+=1
	G={A:dict(H[A])for A in f};g={}
	for(E,W)in G.items():
		for C in W:g.setdefault(C,set()).add(E)
	CL=set(Ao);Ar=[];Bb=0;As=0
	for C in sorted(list(g.keys()),key=lambda cc:len(g[cc])):
		N=list(g.get(C,()))
		if len(N)!=2:continue
		M,F=N
		if M not in G or F not in G:continue
		h,X=G[M][C]
		if h*h+h*X+X*X!=1:continue
		AK,w=h+X,-X;i,a=G[F][C];l,Y=AK*i-w*a,AK*a+w*i+w*a;At={}
		for(D,(m,Z))in G[M].items():
			if D==C:continue
			AL,AM=l*m-Y*Z,l*Z+Y*m+Y*Z;n,o=G[F].get(D,(0,0));b,c=n-AL,o-AM
			if(n,o)==(0,0):At[D]=1
			elif(b,c)==(0,0):At[D]=-1
		if any(len(g.get(A,()))-1+At.get(A,0)<=0 for A in G[M]if A!=C):continue
		Ar.append((C,dict(G[M]),A9[M]))
		for(D,(m,Z))in G[M].items():
			AL,AM=l*m-Y*Z,l*Z+Y*m+Y*Z;n,o=G[F].get(D,(0,0));b,c=n-AL,o-AM
			if b==0 and c==0:G[F].pop(D,None);g[D].discard(F)
			else:G[F][D]=b,c;g.setdefault(D,set()).add(F)
		CM=sum(k(B,C)*A8[A]for(A,(B,C))in G[F].items())%B;Bc=(A9[F]-k(l,Y)*A9[M])%B
		if CM!=Bc:As+=1
		A9[F]=Bc;G.pop(M);Bb+=1
		for D in list(g):g[D].discard(M)
	Au={B for A in G.values()for B in A};Bd=len(CL-(Au|{A for(A,B,B)in Ar}));Av={A:A8[A]for A in Au};Aw=0
	for(AN,Be,CN)in reversed(Ar):
		AZ=k(*Be[AN]);Aq=sum(k(*B)*Av[A]for(A,B)in Be.items()if A!=AN)%B;Av[AN]=(CN-Aq)*pow(AZ,-1,B)%B
		if Av[AN]!=A8[AN]:Aw+=1
	Bf=sum(len(A)for A in G.values());CO=sum(1 for C in G.values()for(A,B)in C.values()if A*A+A*B+B*B==1)
	def CP(a,b):return min(abs(a+A)+abs(b-A)+abs(A)for A in(0,-a,b))
	def Bg(rmap):
		A=0;B=.0
		for D in rmap.values():
			for(E,F)in D.values():A+=1;C=CP(E,F);B+=C*P.C_UNIT_5 if C<=4 else P.C_GEN_1
		return B/max(1,A)
	Bh=Bg(G);CQ={B:A for(A,B)in AC.items()};Bi=[1./CQ[A][2]for A in AI];CR=sum(Bi);Bj=[];A7=.0
	for AH in Bi:A7+=AH/CR;Bj.append(A7)
	import bisect;Bk=[];Bl=[];Bm=[]
	for CS in range(4):
		Bn=Ay.Random(seed+100+CS);j=[]
		for W in H:
			AO={}
			for _ in range(len(W)):
				Ax=AI[min(u-1,bisect.bisect_left(Bj,Bn.random()))];CT,CU=B2[Bn.randrange(6)];n,o=AO.get(Ax,(0,0));b,c=n+CT,o+CU
				if b==0 and c==0:AO.pop(Ax,None)
				else:AO[Ax]=b,c
			if len(AO)>=2:j.append(AO)
		p={};x=set(range(len(j)))
		for(E,W)in enumerate(j):
			for C in W:p.setdefault(C,set()).add(E)
		v=[A for(A,B)in p.items()if len(B)==1]
		while v:
			C=v.pop();N=[A for A in p.get(C,())if A in x]
			if len(N)==1:
				E=N[0];x.remove(E)
				for D in j[E]:
					p[D].discard(E)
					if len(p[D])==1:v.append(D)
		CV={A for(A,B)in p.items()if len(B)>=2};CW=min(len(x),int(1.2*len(CV))+4)
		for E in sorted(x,key=lambda i:len(j[i]),reverse=True):
			if len(x)<=CW:break
			if all(len(p[A])>=2 for A in j[E]):
				x.remove(E)
				for C in j[E]:p[C].discard(E)
		S={A:dict(j[A])for A in x};y={}
		for(E,W)in S.items():
			for C in W:y.setdefault(C,set()).add(E)
		for C in sorted(list(y.keys()),key=lambda cc:len(y[cc])):
			N=list(y.get(C,()))
			if len(N)!=2:continue
			M,F=N
			if M not in S or F not in S:continue
			h,X=S[M][C]
			if h*h+h*X+X*X!=1:continue
			AK,w=h+X,-X;i,a=S[F][C];l,Y=AK*i-w*a,AK*a+w*i+w*a
			for(D,(m,Z))in S[M].items():
				AL,AM=l*m-Y*Z,l*Z+Y*m+Y*Z;n,o=S[F].get(D,(0,0));b,c=n-AL,o-AM
				if b==0 and c==0:S[F].pop(D,None);y[D].discard(F)
				else:S[F][D]=b,c;y.setdefault(D,set()).add(F)
			S.pop(M)
			for D in j[M]:y.get(D,set()).discard(M)
		CX=sum(len(A)for A in S.values());CY=sum(1 for C in S.values()for(A,B)in C.values()if A*A+A*B+B*B==1);Bk.append(len(S));Bl.append(CY/max(1,CX));Bm.append(Bg(S))
	CZ=sum(Bk)/4.;Ca=sum(Bl)/4.;Cb=sum(Bm)/4.;Cc=Al>=150 and u>=25 and Ak==1 and BX==0 and BZ==0 and Bd==0 and AY==0 and Aw==0 and As==0 and len(G)>=len(Au)and Bh<=P.K5_CAP;return dict(ok=Cc,synthetic_incidences=False,zero_nonsplit_inert_aggregation=True,unplanted_f13p12_r28393_chi_kernel=True,unplanted_nullity=Ak,unplanted_rel_errs=BX,genuine_c6_rels=Al,pure_tau6_orbit_cols=u,sm_chi_cols=3,post_peel_cols=len(AX),active_tau6_orbits=len(AX),tau6_unit_hits=BH,repurged_cols=len(BY),post_peel_cols_lost=BZ,post_merge_cols_lost=Bd,repurge_backsub_errors=AY,merge_backsub_errors=Aw,z6_unit_merges=Bb,final_rows=len(G),final_avg_w=round(Bf/max(1,len(G)),2),unit_entry_frac=round(CO/max(1,Bf),4),gen_k5_v4=round(Bh,4),syn_matched_final_rows=round(CZ,2),syn_matched_unit_frac=round(Ca,4),syn_matched_k5_v4=round(Cb,4),kernel_errors=As)
import functools,math,random
from math import ceil,floor,log,log2
@functools.lru_cache(None)
def psurv_iv_tab(u_max=45.,m=512):
	A=math.nextafter;G=int(ceil(u_max*m))+2;U=A(1./m,math.inf);V=A(1./m,.0);J=[1.]*(m+1);K=[1.]*(m+1);F=32;O=[];M=[];W=[.0]*(G+1);P=[.0]*(G+1)
	def X(bi):
		E=bi*F;B=C=.0
		for D in range(E+F-1,E-1,-1):B=A(B+J[D],math.inf);C=A(C+K[D],.0);W[D]=B;P[D]=C
		O.append(B);M.append(C)
	for f in range(m//F):X(f)
	N=J[m];L=K[m]
	for C in range(m+1,G+1):
		if C<=2*m:Y=1.-log(C/m);E=A(A(Y,math.inf),math.inf);D=max(.0,A(A(Y,.0),.0))
		else:
			H=C//F;I=(C-m)//F
			if C==H*F:Z=math.fsum(O[I:H]);a=math.fsum(M[I+1:H])+P[C-m+1]
			else:Z=N+W[C-m]+math.fsum(O[I+1:H]);b=C-m+1;a=L+math.fsum(M[I+1:H])if b==(I+1)*F else L+P[b]+math.fsum(M[I+1:H])
			E=Z/C;D=a/(C-1)
			for c in range(16):E=A(E,math.inf);D=A(D,.0)
			D=max(.0,D)
		J.append(E);K.append(D)
		if(C+1)%F==0:X((C+1)//F-1);N=L=.0
		else:N=A(N+E,math.inf);L=A(L+D,.0)
	d=[.0]*G;e=[.0]*G;Q=[1.]*(G+1);R=[1.]*(G+1)
	for B in range(m,G):
		E=m/float(B)if B<2*m else Q[B+1-m]*m/float(B);D=m/float(B+1)if B<2*m else R[B-m]*m/float(B+1)
		for c in range(4):E=A(E,math.inf);D=A(D,.0)
		d[B]=E;e[B]=max(.0,D);S=Q[B]+U*E;T=R[B]+V*D
		for c in range(4):S=A(S,math.inf);T=A(T,.0)
		Q[B+1]=S;R[B+1]=T
	return m,U,V,tuple(J),tuple(K),tuple(d),tuple(e),tuple(J[::-1]),tuple(K[::-1])
@functools.lru_cache(None)
def psurv_iv(m_bits,B,T,mesh=512):
	P=m_bits
	if T<.0:return .0,.0
	if P<=T:return 1.,1.
	L=P/B;X=T/B;Y=max(45.,ceil(L)+2.);A,Z,a,Q,b,R,S,c,d=psurv_iv_tab(Y,mesh);I=int(floor(L*A));C=int(ceil(L*A));O=min(L,X);D=len(Q);U=Q[I];V=b[C]if C<D else .0
	if O<=1.:return min(1.,V),min(1.,U)
	M=int(ceil(O*A));N=int(floor(O*A));W=math.sumprod;E=I-1 if I-1<M else M;J=.0
	if E>A:J+=W(c[D-I+A:D-I+E],R[A:E])
	if M>(E if E>A else A):J+=math.fsum(R[max(A,E):M])
	F=A if C-A<D else C-D+1;G=C if C<N else N;K=.0
	if G>F:K+=W(d[D-1-C+F:D-1-C+G],S[F:G])
	if N>(G if G>F else F):K+=math.fsum(S[max(F,G):N])
	H=math.nextafter
	for e in range(8):J=H(J,math.inf);K=H(K,.0)
	return min(1.,max(.0,H(V+H(a*max(.0,K),.0),.0))),min(1.,H(U+H(Z*J,math.inf),math.inf))
def psurv512(m_bits,B,T,mesh=512):return psurv_iv(m_bits,B,T,mesh)[1]
def verify_v43_descent(D,M,segs):
	L=segs;import re;Ar=L.get('desc_zero16',-1);As=L.get('desc_scan8_u8',-1);Bv=L.get('desc_sieve_hit',-1);At=L.get('mont5_cios_prologue',-1);Au=L.get('mont5_cios_round',-1);Av=L.get('mont5_cios_epilogue',-1);Aw=L.get('mont2_cios',-1);Ax=L.get('bingcd_gen_loop',-1);Ay=L.get('modinv_gen_loop',-1);AH=At+5*Au+Av;Bw=L.get('cios_gen_mul_inner',-1);Bx=L.get('cios_gen_red_inner',-1);AI=D.MONT5;AJ=25.;O=3*AI+5*AJ;b=2*AI+3*AJ;R=5e1;v=M.P_BN254-2;Az=(4+(v.bit_length()-1)+(v.bit_count()-1))*AI+2*AJ;A_=254*(3*O+4*b+10*R)+254*(12*O+4*b+10*R)+Az+3*O+b;B0=8*O+15*R;AK=3*B0+O+30*R;By=12*O+3*b+15*R+Az;Bz=4*B0+O+By+10*R;B1=M.P_BN254**4-M.P_BN254**2+1;w=B1//M.R_BN254;B_=w.bit_length()-1+(w.bit_count()-1);B2=65*(2*AK+6*O+6*b+20*R)+67*(AK+12*O+4*b+20*R)+((B_+4)*AK+Bz+2e5);C0=v.bit_length()==254 and v.bit_count()==110 and B1%M.R_BN254==0 and w.bit_length()==761 and w.bit_count()==389 and A_<=5422256. and A_<=M.G2_OUTPUT_BOUND and B2<=40297005. and B2<=2e4*2.**11;AL=Ar==18 and As==22 and Bv==9 and At==18 and Au==48 and Av==21 and AH==279 and AH+7<=D.MONT5 and Aw==47 and Ax==49 and Ay==91 and Bw==11 and Bx==11 and Ar/256.<=D.ZERO and As/64.<=D.SCAN and C0;C1=Path(__file__).resolve().parent/'kernels_v32.s';C2=C1.read_text(encoding='ascii').splitlines();AM=[];AN=False;x=[];AO=[];y=[];AP=[];z=[];B3={};B4={};B5={};K=None;A0=0;AQ=[]
	for C3 in C2:
		A=C3.split('//',1)[0].strip()
		if not A:continue
		if A.startswith('.macro MONT5_ROUND'):AN=True;continue
		if AN:
			if A.startswith('.endm'):AN=False
			else:AM.append(A)
			continue
		if A=='seg_mont5_cios:':K='m5';continue
		if A=='seg_mont2_cios:':K='m2';continue
		if A=='seg_bingcd_gen:':K='bg';continue
		if A=='seg_cios_gen:':K='cg';continue
		if A=='seg_modinv_gen:':K='mi';continue
		if K=='m5':
			if A.startswith('.rept '):A0=int(A.split()[1]);AQ=[];continue
			if A=='.endr':
				for c in range(A0):
					for B6 in AQ:x.extend(AM if B6=='MONT5_ROUND'else[B6])
				A0=0;continue
			if A0:AQ.append(A);continue
			if A.endswith(':')or A.startswith('.'):continue
			x.extend(AM if A=='MONT5_ROUND'else[A])
			if A=='ret':K=None
		elif K=='m2':
			if A.endswith(':')or A.startswith('.'):continue
			AO.append(A)
			if A=='ret':K=None
		elif K=='bg':
			if A.endswith(':'):B3[A[:-1]]=len(y);continue
			if A.startswith('.'):continue
			y.append(A)
			if A=='ret':K=None
		elif K=='cg':
			if A.endswith(':'):B4[A[:-1]]=len(AP);continue
			if A.startswith('.'):continue
			AP.append(A)
			if A=='ret':K=None
		elif K=='mi':
			if A.endswith(':'):B5[A[:-1]]=len(z);continue
			if A.startswith('.'):continue
			z.append(A)
			if A=='ret':K=None
	B=(1<<64)-1
	def F(s):s=s.strip().rstrip(',');return 32 if s=='xzr'else 31 if s=='sp'else int(s[1:])
	AR={};C4={'mul':26,'umulh':27,'orr':28,'add':29,'sub':30,'adds':31,'adcs':32,'adc':33,'subs':34,'sbcs':35,'csel':36};B7={'b':1,'b.ne':2,'b.lo':3,'b.hs':4,'b.pl':5};B8={'add':13,'sub':14,'lsr':15,'lsl':16,'subs':17}
	def C5(prog,labels):
		H=labels;I=id(prog)
		if I in AR:return AR[I]
		C=[]
		for N in prog:
			B,S,D=N.partition(' ');D=D.strip()
			if B=='ret':C.append((0,0,0,0,0))
			elif B in B7:C.append((B7[B],H[D],0,0,0))
			elif B in('cbz','cbnz'):J,K=D.split(',');C.append((6 if B=='cbz'else 7,F(J),H[K.strip()],0,0))
			elif B in('tbz','tbnz'):J,P,K=D.split(',');C.append((8 if B=='tbz'else 9,F(J),int(P.strip().lstrip('#')),H[K.strip()],0))
			elif B=='cmp':E,G=D.split(',');C.append((10,F(E),F(G),0,0))
			elif B=='mov':E,G=D.split(',');L=G.strip();C.append((11,F(E),int(L.lstrip('#')),0,0)if L.startswith('#')else(12,F(E),F(L),0,0))
			elif B in B8 and'#'in D:E,G,Q=D.split(',');C.append((B8[B],F(E),F(G),int(Q.strip().lstrip('#')),0))
			elif B in('stp','ldp'):A=re.match('(\\w+),\\s*(\\w+),\\s*\\[(\\w+),\\s*#(\\d+)\\]$',D);C.append((18 if B=='stp'else 19,F(A.group(1)),F(A.group(2)),F(A.group(3)),int(A.group(4))))
			elif B in('str','ldr'):
				M=20 if B=='str'else 23;A=re.match('(\\w+),\\s*\\[(\\w+),\\s*(\\w+),\\s*lsl\\s*#(\\d+)\\]$',D)
				if A:C.append((M,F(A.group(1)),F(A.group(2)),F(A.group(3)),int(A.group(4))))
				else:
					A=re.match('(\\w+),\\s*\\[(\\w+)\\],\\s*#(\\d+)$',D)
					if A:C.append((M+1,F(A.group(1)),F(A.group(2)),int(A.group(3)),0))
					else:A=re.match('(\\w+),\\s*\\[(\\w+)(?:,\\s*#(\\d+))?\\]$',D);C.append((M+2,F(A.group(1)),F(A.group(2)),int(A.group(3)or 0),0))
			else:
				E,G,R=[A.strip()for A in D.split(',')[:3]];O=C4.get(B)
				if O is None:raise ValueError(f"unsupported op: {N}")
				C.append((O,F(E),F(G),F(R),0))
		AR[I]=C;return C
	def f(prog,x_init,mem,labels=None,max_ops=2000000):
		K=mem;P=C5(prog,labels or{});A=list(x_init)+[0];H=L=M=I=O=0
		while O<max_ops:
			C,D,E,G,N=P[I];I+=1;O+=1
			if C==0:return A[:32],O
			elif C==1:I=D
			elif C==2:
				if not L:I=D
			elif C==3:
				if not H:I=D
			elif C==4:
				if H:I=D
			elif C==5:
				if not M:I=D
			elif C==6:
				if A[D]==0:I=E
			elif C==7:
				if A[D]!=0:I=E
			elif C==8:
				if not A[D]>>E&1:I=G
			elif C==9:
				if A[D]>>E&1:I=G
			elif C==10:F=A[D]+(A[E]^B)+1;H=F>>64;L=int(F&B==0);M=F>>63&1
			elif C==11:A[D]=E&B;A[32]=0
			elif C==12:A[D]=A[E];A[32]=0
			elif C==13:A[D]=A[E]+G&B;A[32]=0
			elif C==14:A[D]=A[E]-G&B;A[32]=0
			elif C==15:A[D]=A[E]>>G&B;A[32]=0
			elif C==16:A[D]=A[E]<<G&B;A[32]=0
			elif C==17:F=A[E]+(G^B)+1;H=F>>64;L=int(F&B==0);M=F>>63&1;A[D]=F&B;A[32]=0
			elif C==18:J=A[G]+N;K[J]=A[D];K[J+8]=A[E]
			elif C==19:J=A[G]+N;A[D]=K[J];A[E]=K[J+8];A[32]=0
			elif C==20:K[A[E]+(A[G]<<N)]=A[D]
			elif C==21:J=A[E];K[J]=A[D];A[E]=J+G&B;A[32]=0
			elif C==22:K[A[E]+G]=A[D]
			elif C==23:A[D]=K[A[E]+(A[G]<<N)];A[32]=0
			elif C==24:J=A[E];A[D]=K[J];A[E]=J+G&B;A[32]=0
			elif C==25:A[D]=K[A[E]+G];A[32]=0
			elif C==26:A[D]=A[E]*A[G]&B;A[32]=0
			elif C==27:A[D]=A[E]*A[G]>>64&B;A[32]=0
			elif C==28:A[D]=A[E]|A[G];A[32]=0
			elif C==29:A[D]=A[E]+A[G]&B;A[32]=0
			elif C==30:A[D]=A[E]-A[G]&B;A[32]=0
			elif C==31:F=A[E]+A[G];H=F>>64;L=int(F&B==0);M=F>>63&1;A[D]=F&B;A[32]=0
			elif C==32:F=A[E]+A[G]+H;H=F>>64;L=int(F&B==0);M=F>>63&1;A[D]=F&B;A[32]=0
			elif C==33:A[D]=A[E]+A[G]+H&B;A[32]=0
			elif C==34:F=A[E]+(A[G]^B)+1;H=F>>64;L=int(F&B==0);M=F>>63&1;A[D]=F&B;A[32]=0
			elif C==35:F=A[E]+(A[G]^B)+H;H=F>>64;L=int(F&B==0);M=F>>63&1;A[D]=F&B;A[32]=0
			elif C==36:A[D]=A[E]if H==1 else A[G];A[32]=0
		raise RuntimeError('max_ops exceeded')
	G=random.Random(20261143);B9=len(x)==279;BA=len(AO)==47;BB=len(y)==49;BC=len(z)==Ay;BD=True
	for I in range(32):
		C=G.getrandbits(318)|1<<318|1;C6=-pow(C,-1,1<<64)&B;BE=C-1 if I==0 else G.randrange(C);A1=C-1 if I==0 else G.randrange(C);g={};h,i,j,BF,BG=4096,8192,12288,16384,32768
		for E in range(5):g[h+8*E]=BE>>64*E&B;g[i+8*E]=A1>>64*E&B;g[BF+8*E]=C>>64*E&B
		P=[0]*32;P[5],P[6],P[7],P[8],P[9],P[31]=h,i,C6,j,BF,BG
		for BH in range(19,28):P[BH]=0xdeadbeef0000+BH
		BI,C7=f(x,P,g);C8=sum(g[j+8*A]<<64*A for A in range(5))
		if C7!=279 or C8!=BE*A1*pow(1<<320,-1,C)%C or BI[31]!=BG or any(BI[A]!=P[A]for A in range(19,28)):B9=False
		Q=(1<<128)-1 if I==0 else(1<<127)+1 if I==1 else G.getrandbits(126)|1<<126|1;C9=-pow(Q,-1,1<<64)&B;AS=Q-1 if I==0 else G.randrange(Q);AT=Q-1 if I==0 else G.randrange(Q);S=[0]*32;S[0],S[1],S[2],S[3],S[4],S[5],S[6]=AS&B,AS>>64&B,AT&B,AT>>64&B,Q&B,Q>>64&B,C9;BJ,CA=f(AO,S,{})
		if CA!=47 or BJ[0]|BJ[1]<<64!=AS*AT*pow(1<<128,-1,Q)%Q:BA=False
	for H in(2,5,6,8,12,20):
		T=64*H;CB=int(D.lehmer_gcd_ir(T))
		for I in range(3):
			Y=G.getrandbits(T)|1;A2=0 if I==0 else Y if I==1 else G.getrandbits(T)
			if I==2:
				A3=G.getrandbits(min(60,T//2))|1;Y=Y//A3*A3;A2=A2//A3*A3
				if Y&1==0:Y|=1
			d={};k,l=4096,12288
			for E in range(H):d[k+8*E]=Y>>64*E&B;d[l+8*E]=A2>>64*E&B
			A4=[0]*32;A4[0],A4[1],A4[2]=k,l,H;CC,CD=f(y,A4,d,B3);CE=sum(d[CC[0]+8*A]<<64*A for A in range(H))
			if CE!=math.gcd(Y,A2)or CD>CB:BB=False
	for H in(6,8,12,20):
		AU=64*H;CF=int(D.mul_L(H))
		for I in range(4):
			U=G.getrandbits(AU-2)|1<<AU-2|1;CG=-pow(U,-1,1<<64)&B;BK=U-1 if I==0 else G.randrange(U);BL=U-1 if I==0 else G.randrange(U);m={};h,i,n,CH,j=4096,12288,20480,28672,36864
			for E in range(H):m[h+8*E]=BK>>64*E&B;m[i+8*E]=BL>>64*E&B;m[n+8*E]=U>>64*E&B
			V=[0]*32;V[0],V[1],V[2],V[3],V[4],V[5],V[6]=h,i,n,CH,H,CG,j;c,CI=f(AP,V,m,B4);AV=sum(m[j+8*A]<<64*A for A in range(H))
			if AV!=BK*BL*pow(1<<AU,-1,U)%U or CI>CF:BD=False
	for H in(2,5,6,8,12,20):
		T=64*H;CJ=int(D.modinv_ir(T))
		for I in(0,2)if H==20 else(0,1,2):
			J=G.getrandbits(T-1)|1<<T-1|1;Z=G.randrange(1,J)if I<2 else G.randrange(1,J//3)*3%J or 3
			if I==2 and J%3!=0:J=J//3*3;J=J if J&1 else J+3
			if Z>=J:Z%=J
			if Z==0:Z=1
			W={};k,l,BM,BN,n=4096,12288,20480,28672,36864
			for E in range(H):W[k+8*E]=J>>64*E&B;W[l+8*E]=Z>>64*E&B;W[BM+8*E]=0;W[BN+8*E]=1 if E==0 else 0;W[n+8*E]=J>>64*E&B
			a=[0]*32;a[0],a[1],a[2],a[3],a[4],a[5]=k,l,H,BM,BN,n;BO,CK=f(z,a,W,B5);AV=sum(W[BO[0]+8*A]<<64*A for A in range(H));CL=sum(W[BO[3]+8*A]<<64*A for A in range(H));BP=math.gcd(Z,J)
			if AV!=BP or CK>CJ or BP==1 and CL*Z%J!=1:BC=False
	AL=AL and B9 and BA and BB and BC and BD;CM=Ax;AW=M.load_data()['descent'];CN=AW['plan59'];CO={int(A):B for(A,B)in AW['cells'].items()};CP=max(A[1]for A in D.QBINS);CQ=max(A[2]for A in D.QBINS);AX=AY=o=AZ=.0
	for(BQ,BR,BS,c,c,c,c,CR,c)in CN:BT=log2(CO[CR]);BU=(BR+BT)/12.;AX=max(AX,D.DF*BU+CP-(BQ if BS=='f'else .0));AY=max(AY,12.*BU+D.LOG2U6+CQ-(BQ if BS=='g'else .0));o=max(o,float(BR));AZ=max(AZ,BT)
	Aa=max(ceil(AX),ceil(AY),ceil(AW['boot_dp']['ng_hi']));BV=max(Aa,2560);A5=ceil(BV/6e1)+1;CS=A5/8.*(M.C_NORM_G+M.C_NORM_F)+6e1*A5*A5;CT=309*(ceil(BV/64.)*5+10);BW=D.mul_L(max(4,ceil(o/64.)));CU=A5>=44 and CS<D.NORMS and CT<D.TDIV and 64.*o*BW<D.PRIME and 6.*o*45.*BW<D.IDEAL and AZ<=82.;Ab=True;BX=1.
	for d in range(400,1900,100):
		CV,Ac=psurv_iv(d,59.,5e2,512)
		if not(.0<=CV<=Ac<=1. and Ac<=BX+1e-15):Ab=False
		BX=Ac
	BY=G.randrange(1<<592,1<<593)|1;BZ=G.randrange(1<<591,1<<592)|1;C=BY*BZ
	def X(A,B):
		D=[0]*(len(A)+len(B)-1)
		for(F,E)in enumerate(A):
			if E:
				for(G,H)in enumerate(B):D[F+G]+=E*H
		return[A%C for A in D]
	def Ad(arr,w,p,inv=False):
		A=len(arr);M=pow(w,-1,p)if inv else w;N=A.bit_length()-1;C=[arr[int(f"{A:0{N}b}"[::-1],2)]%p for A in range(A)];E=1
		while E<A:
			G=E;E<<=1;O=pow(M,A//E,p);D=1;I=[]
			for Q in range(G):I.append((D,(D<<64)//p));D=D*O%p
			for H in range(0,A,E):
				for F in range(G):D,P=I[F];J=C[H+F];K=C[H+F+G];L=(K*D-(K*P>>64)*p&B)%p;C[H+F]=(J+L)%p;C[H+F+G]=(J-L)%p
		return[B*pow(A,-1,p)%p for B in C]if inv else C
	Ba=[G.randrange(C)for A in range(16)];Bb=[G.randrange(C)for A in range(16)];p=ceil((2*C.bit_length()+5)/64.);Ae=32*p;Af=1<<ceil(log2(Ae));A6,q,r=M.NTT3_PRIMES
	def Bc(P):
		A=[0]*Af
		for(E,F)in enumerate(P):
			D=F%C;G=E*p
			for H in range(p):A[G+H]=D&B;D>>=64
		return A
	CW,CX=Bc(Ba),Bc(Bb);CY=pow(A6,-1,q);CZ=pow(A6*q%r,-1,r);A7=[]
	for N in M.NTT3_PRIMES:
		A8=2
		while pow(A8,(N-1)//2,N)==1 or pow(A8,(N-1)//3,N)==1:A8+=1
		Ag=pow(A8,(N-1)//Af,N);Ca,Cb=Ad(CW,Ag,N),Ad(CX,Ag,N);A7.append(Ad([Ca[A]*Cb[A]%N for A in range(Af)],Ag,N,True))
	Bd=[0]*Ae;Ah=0
	for A9 in range(Ae):Be,Cc,Cd=A7[0][A9],A7[1][A9],A7[2][A9];Bf=Be+(Cc-Be)%q*CY%q*A6;Bg=Bf+(Cd-Bf%r)%r*CZ%r*A6*q+Ah;Bd[A9]=Bg&B;Ah=Bg>>64
	Bh=[]
	for Ce in range(31):
		Ai=0
		for Cf in range(p-1,-1,-1):Ai=Ai<<64|Bd[Ce*p+Cf]
		Bh.append(Ai%C)
	Bi=Ah==0 and Bh==X(Ba,Bb)
	def Bj():
		while True:
			A=G.randrange(1,C)
			if math.gcd(A,C)==1:return A
	Cg=[G.randrange(1,C)for A in range(32)];Ch=[Bj()for A in range(32)];Ci=[G.randrange(1,C)for A in range(32)];Cj=[Bj()for A in range(32)];AA=Ch+Cj;AB=[AA[0]];s=0
	for Aj in AA[1:]:AB.append(AB[-1]*Aj%C);s+=1
	AC=pow(AB[-1],-1,C);AD=[0]*64
	for Ak in range(63,0,-1):AD[Ak]=AC*AB[Ak-1]%C;AC=AC*AA[Ak]%C;s+=2
	AD[0]=AC;Bk=[A*B%C for(A,B)in zip(Cg,AD[:32])];Bl=[A*B%C for(A,B)in zip(Ci,AD[32:])];s+=64;Al=[1];e=list(AA);e[3]=e[3]*BY%C;e[47]=e[47]*BZ%C;t=[e[0]]
	for Aj in e[1:]:t.append(t[-1]*Aj%C)
	AE,Am,Bm=0,63,0
	while AE<Am:
		An=(AE+Am)//2;Bm+=1
		if math.gcd(t[An],C)>1:Am=An
		else:AE=An+1
	Ck=math.gcd(t[-1],C)==C and Bm<=6 and 1<math.gcd(t[AE],C)<C
	for A1 in Bk:Al=X(Al,[-A1%C,1])
	def Ao(pts):
		A=pts
		if len(A)==1:return[-A[0]%C,1],None,None,1
		B=len(A)//2;D=Ao(A[:B]);E=Ao(A[B:]);return X(D[0],E[0]),D,E,len(A)
	Bn=Ao(Bl);Bo=Bn[0];Cl=[(Al[A]-Bo[A])%C for A in range(32)];Cm=Bo[::-1];AF=[1];Ap=1
	while Ap<32:AG=min(2*Ap,32);Bp=X(Cm[:AG],AF)[:AG];Bq=[-A%C for A in Bp];Bq[0]=(2-Bp[0])%C;AF=X(AF,Bq)[:AG];Ap=AG
	Cn=X(Cl[::-1],AF)[:32]
	def Aq(Mc,nd):
		E,A,B,D=nd
		if A is None:return[Mc[0]%C]
		return Aq(X(Mc,B[0][::-1])[B[3]:D],A)+Aq(X(Mc,A[0][::-1])[A[3]:D],B)
	Br=Ck and s==4*64-3 and s<=4*64 and all(B==math.prod((A-B)%C for B in Bk)%C for(A,B)in zip(Bl,Aq(Cn,Bn)));u=max(D._stage_ops(A,A,A)[3]for A in D.WHEELS);Bs=ceil(max(Aa,1225)/64.);Co=ceil((2*max(Aa,1225)+20)/64.);Bt=1<<ceil(log2(2*(u+1)*Co));Bu=6*Bt*8+4*u*Bs*8*4+4*u*Bs*8;Cp=u==46080 and all(D._stage_ops(C,E,A)[3]<=u for(B,C,E,A,B)in D._ECM_BEST_C.values()if A>0);Cq=AL and CU and Ab and Bi and Br and Cp and Bt<1<<54 and Bu<1<<30;return dict(ok=Cq,cios5=AH,mont2=Aw,bingcd_gen=CM,psurv_iv_ok=Ab,kron_ntt3_ok=Bi,bls_mpe_ok=Br,ecm_ram_mib=round(Bu/(1<<20),2))
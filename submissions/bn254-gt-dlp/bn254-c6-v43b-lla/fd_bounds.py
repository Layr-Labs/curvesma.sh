import math,random,sys,time
from fractions import Fraction as Fr
Q32=1<<32
Q38=1<<38
M64=(1<<64)-1
BMIN=Fr(753729000)*Fr(13,20)**11
TMAX=int(2.**46.803)+2
KAPPA=Fr(99,100)-Fr(1,4)
BETA_LO=BMIN/TMAX
LUT=[0]*512+[math.ceil(math.sqrt(A+1)*2.**52.5)+(1<<31)for A in range(512,2048)]
MD0=2*0xe8d4a51027
def sqrt_up(q,p=60):q=Fr(q);A=q.numerator<<2*p;B=-(-A//q.denominator);return Fr(math.isqrt(B)+1,1<<p)
def sqrt_dn(q,p=60):q=Fr(q);return Fr(math.isqrt((q.numerator<<2*p)//q.denominator),1<<p)
def isqrt_lut(r2):
	if not 1<=r2<1<<40:raise OverflowError('r2 outside [1, 2^40): %d'%r2)
	A=64-r2.bit_length()&~1;return LUT[r2<<A>>53]>>(A>>1)
def lut_certificate():
	D=0;E=Fr(0);I=0
	for A in range(0,63,2):
		F,J=max(1,1<<max(0,62-A)),min(1<<40,1<<64-A)
		if A==62:F=1
		if F>=J:continue
		for B in range(512,2048):
			N=-(-(B<<53)>>A)if A<=53 else-(-(B<<53)//(1<<A));O=-(-(B+1<<53)//(1<<A))-1;G,H=max(N,F),min(O,J-1)
			if G>H:continue
			for K in(G,H):
				L=64-K.bit_length()&~1
				if L!=A or K<<L>>53!=B:D+=1
			I+=1;C=LUT[B]>>(A>>1)
			if C*C<=H<<52:D+=1
			M=Fr(C*C,G<<52)
			if M>E:E=M
	return D,E,I
def uniform_bounds(LAM,kappa=KAPPA,beta_lo=BETA_LO,grid=400):
	R=kappa;N=grid;F=beta_lo;H=sqrt_up(LAM);B=[None]*12
	for D in range(0,11):
		S={A:sqrt_up(Fr(1)/R**(A-D))for A in range(D+1,12)};I={}
		for E in range(11,D,-1):I[E]=S[E]+Fr(3,2)+sum((I[A]+1)/2 for A in range(E+1,12))
		O={}
		for E in range(11,D,-1):W=sum(I[A]for A in range(E+1,12))+1;O[E]=sum(O[A]for A in range(E+1,12))/2+W/Q32+H*(S[E]+Fr(1,64)+Fr(1,1<<20))
		X=sum(I.values())+1;Y=sum(O.values());B[D]=(X+Y)/Q32
	B[11]=Fr(1,Q32);P={}
	for D in range(2,12):Q=1+(1<<19)*B[D];P[D]=((1<<20)*Q*(1+Fr(1,1<<37))+2*Q*Q+2)/Q38
	T=sum(P.values());A=H*sqrt_up(1+Fr(1,1<<12)+Fr(1,1<<38))+B[1]
	def J(b):C=sqrt_up(Fr(1)/b);D=Fr(1)/b;E=(2*A*C+1)**2+937*(2*A*C)+1+(1+Fr(1,1<<26))*A*A*D+2*A*C+4+(1<<20)*A;return E/Q38+B[1]*(6*A+2+B[1])
	def U(Fv,b0):A=(1+Fv)*(1+Fr(1,Q38))+b0/4096;return T+(H*sqrt_up(A)+B[0])**2-1
	C=[F*(Fr(1)/F)**Fr(A,N)if False else None for A in range(0)];C=[];a=sqrt_dn(sqrt_dn(sqrt_dn(Fr(1)/F)));Z=Fr(1)+Fr(1,N//8 if N>=16 else 2);K=F
	while K<1:C.append(K);K=K*Z
	C.append(Fr(1));G=None
	for L in range(len(C)-1):
		V=U(J(C[L]),min(Fr(1),C[L+1]/R))
		if G is None or V>G[0]:G=V,C[L],C[L+1]
	M={};M['R1K1_max_log2']=math.log2((1<<39)*(2*A+1)*Q32*B[1]+(1<<33));M['r2_0_max_log2']=math.log2((1+J(F))*(Fr(4096)/F+Fr(1,1<<26))+1);M['dr1_max_log2']=math.log2(Q38*(2*A+2)+Q38);return dict(sL=H,z=B,e=P,D1=T,W=A,Fb=J,eps_at=U,eps=G[0],eps_where=G[1:],F_lo=J(F),npts=len(C),WR=M)
G6=[[753729,251243,251243,251244,251242,251244],[251243,1163069,-455913,-455913,517661,517662],[251243,-455913,1163069,189495,189494,-455913],[251244,-455913,189495,1542060,127747,-394165],[251242,517661,189494,127747,1542058,455914],[251244,517662,-455913,-394165,455914,1542060]]
KAP_L=Fr(74,100)
P_E_CHARGED=3<<12
def _det(M):
	M=[[Fr(A)for A in A]for A in M];C=Fr(1)
	for A in range(len(M)):
		C*=M[A][A]
		for B in range(A+1,len(M)):D=M[B][A]/M[A][A];M[B]=[A-D*B for(A,B)in zip(M[B],M[A])]
	return C
def lemma_e():
	D=Fr(1<<124)*Fr(1764)**6*Fr(1000)**6*_det(G6)**2;E=Fr(1764)*Fr(56880,100)**2*196608-1;B=[]
	for A in range(12):
		C=BMIN**A
		for F in range(A,12):C*=max(BMIN,KAP_L**(F-A)*E)
		B.append(math.log2(C/D))
	return B
def exact_nodes_bound(UB,eta_max,lem):
	B=1+UB['D1']+Fr(12,1<<38)
	def E(b):A=2*sqrt_up(UB['sL']**2*(B/b+Fr(1,4096)+B/(1<<38)))+2*eta_max+Fr(1,1<<31);return A.numerator//A.denominator+1
	C=0
	for A in range(12):
		if lem[A]>0:continue
		D=1
		for F in range(A+1,12):D*=E(KAP_L**(F-A))
		C+=D
	return C
def fp_nodes_above_bound(UB,eta_max,lem):
	B=1+UB['D1']+Fr(12,1<<38)
	def E(b):A=2*sqrt_up(UB['sL']**2*(B/b+Fr(1,4096)+B/(1<<38)))+2*eta_max+Fr(1,1<<31);return A.numerator//A.denominator+1
	C=0
	for A in range(12):
		if lem[A]>0:continue
		for F in range(A+1,12):
			D=1
			for G in range(F+1,12):D*=E(KAP_L**(G-A))
			C+=D
	return C
FP_NODE_MAX_IR=99
C_LAT_EXC_IR=int(2.**31.2)
def direct_q12_lut_certificate():
	C=(Fr(1)+Fr(1,1<<26))**2;D=0
	for E in(0,1,2,3,7,15,63,1024,(1<<20)-1,(1<<36)-1,(1<<37)-1):
		B=E+1;A=math.isqrt(B<<52)+1
		if A*A<=B<<52 or Fr(A*A,B<<52)>C:D+=1
	return D,C
def uncond_xmax_bound(eps_or):
	C=sqrt_up(1+eps_or);D=[sqrt_up(Fr(1)/KAPPA**A)for A in range(12)];A=[Fr(0)]*12
	for B in range(11,-1,-1):A[B]=C+Fr(1,2)+sum(Fr(1,2)*A[C]*D[C-B]for C in range(B+1,12))
	E=Fr(1)/sqrt_dn(BETA_LO);return[A[B]*E for B in range(12)]
def swar21_max_coord_bound(eps_or):
	E=eps_or;C=6;A=[[Fr(G6[A][B])for B in range(C)]+[Fr(int(A==B))for B in range(C)]for A in range(C)]
	for B in range(C):
		I=next(C for C in range(B,C)if A[C][B]!=0);A[B],A[I]=A[I],A[B];P=A[B][B];A[B]=[A/P for A in A[B]]
		for D in range(C):
			if D!=B and A[D][B]:Q=A[D][B];A[D]=[A[D][C]-Q*A[B][C]for C in range(2*C)]
	J=max(A[B][C+B]for B in range(C));R=sqrt_up((1+E)*Fr(TMAX,1000)*J);S=sqrt_up(1+E);K=[sqrt_up(Fr(1)/KAPPA**A)for A in range(12)];F=[Fr(0)]*12
	for G in range(11,-1,-1):F[G]=S+Fr(1,2)+sum(Fr(1,2)*F[A]*K[A-G]for A in range(G+1,12))
	H=Fr(1)+E
	for L in range(1,12):
		M=Fr(0)
		for T in range(L):N=sum(Fr(1,2)*F[A]*K[A-T]for A in range(L,12));M+=N*N
		O=1+E+M
		if O>H:H=O
	U=(sqrt_up(H)+sqrt_up(Fr(5,4)))*sqrt_up(Fr(TMAX,1000)*J);return R,U
def exact_node_ops():A=lambda a,b:4*a*b+2*(a+b);B=lambda n:A(1,n)+9*n+10;C=8*(A(8,1)+3*9)+B(9);D=A(1,8)+3*9+A(9,9)+3*19+3*19+10;E=3*19+B(16)+100;return C+A(1,16)+4*D+3*E+20
def completeness_margin():
	B=None;I=sqrt_up(Fr(1)/BETA_LO)+Fr(3,2);C={}
	for E in range(11,1,-1):C[E]=I+sum((C[A]+1)/2 for A in range(E+1,12))
	F=(sum(C.values())+1)/Q32
	for A in range(1,65):
		for G in(Fr(4),Fr(1<<38)):
			for D in(Fr(0),F):
				H=G*D*(A-1)-2*A*D+1+Fr((3*A+1)*(A-1),4)
				if B is None or H<B[0]:B=H,G,D,A
	return B+(F,)
def main():N=time.time();J,O,P=lut_certificate();K,L=direct_q12_lut_certificate();A=uniform_bounds(L);D=completeness_margin();F=math.log2((1+A['D1'])*(Fr(4096)/BETA_LO+Fr(1,1<<26))+1);A['WR']['r2_k1_max_log2']=F;G=uncond_xmax_bound(A['eps']);Q=max(math.log2(A)for A in G);E=dict(lut7_buckets=P,lut7_outward_errs=J,lut7_lambda=float(O),lut_direct_q12_errs=K,lut_direct_q12_lambda=float(L),lut_direct_q12_in_bounds_ok=A['WR']['r2_0_max_log2']<37. and F<37.,uncond_xmax_max_log2=round(Q,3),uncond_xmax_lt_2p22=max(G)<1<<22,beta_lo_log2=math.log2(BETA_LO),z_max_log2=math.log2(max(A['z'][:11])),delta1_levels_11_to_2=float(A['D1']),F_beta_lo=float(A['F_lo']),eps_or_bound=float(A['eps']),word_widths=A['WR'],pad1_min_slack=float(D[0]));B=lemma_e();C=exact_nodes_bound(A,D[4],B);H=fp_nodes_above_bound(A,D[4],B);M=exact_node_ops();R,I=swar21_max_coord_bound(A['eps']);E.update(lemma_e_log2_margin_k0_k3=[round(A,3)for A in B[:4]],exact_levels_possible=[A for A in range(12)if B[A]<=0],exact_nodes_per_orbit_max=C,exact_nodes_log2=round(math.log2(C),3),exact_node_ops=M,exact_node_price_charged=P_E_CHARGED,exact_per_orbit_ir_log2=round(math.log2(C*P_E_CHARGED),3),fp_nodes_above_exact_max=H,exact_level_work_ir=C*P_E_CHARGED+H*FP_NODE_MAX_IR,exact_level_allowance_ir=C_LAT_EXC_IR,swar21_leaf_coord_ub=round(float(R),2),swar21_max_coord_ub=round(float(I),2),swar21_zero_overflow_ok=I<1<<20);E['ok']=J==0 and K==0 and A['WR']['r2_0_max_log2']<37. and F<37. and max(G)<1<<22 and A['eps']<=Fr(23389,10000000)and D[0]>0 and min(B[:3])>0 and 2*M<=P_E_CHARGED and C*P_E_CHARGED+H*FP_NODE_MAX_IR<=C_LAT_EXC_IR and I<1<<20;E['seconds']=round(time.time()-N,3);return E
if __name__=='__main__':import json;o=main();print(json.dumps(o,default=str));sys.exit(0 if o['ok']else 1)
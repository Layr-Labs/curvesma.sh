T='r2_0_max_log2'
R=int
O='eps'
M='WR'
K='D1'
N=len
L=min
J=None
G=round
F=float
E=sum
D=max
A=range
import math as B,random,sys,time
from fractions import Fraction as Fr
import stor_bound as Q
Q32=1<<32
Q38=1<<38
V=(1<<64)-1
C=Fr(753729000)*Fr(13,20)**11
TMAX=R(2.**46.803)+2
KAPPA=Fr(99,100)-Fr(1,4)
H=C/TMAX
LUT=[0]*512+[B.ceil(B.sqrt(o+1)*2.**52.5)+(1<<31)for o in A(512,2048)]
def sqrt_up(q,p=60):q=Fr(q);A=q.numerator<<2*p;C=-(-A//q.denominator);return Fr(B.isqrt(C)+1,1<<p)
def sqrt_dn(q,p=60):q=Fr(q);return Fr(B.isqrt((q.numerator<<2*p)//q.denominator),1<<p)
def c():
	F=0;G=Fr(0);K=0
	for B in A(0,63,2):
		H,M=D(1,1<<D(0,62-B)),L(1<<40,1<<64-B)
		if B==62:H=1
		if H>=M:continue
		for C in A(512,2048):
			Q=-(-(C<<53)>>B)if B<=53 else-(-(C<<53)//(1<<B));R=-(-(C+1<<53)//(1<<B))-1;I,J=D(Q,H),L(R,M-1)
			if I>J:continue
			for N in(I,J):
				O=64-N.bit_length()&~1
				if O!=B or N<<O>>53!=C:F+=1
			K+=1;E=LUT[C]>>(B>>1)
			if E*E<=J<<52:F+=1
			P=Fr(E*E,I<<52)
			if P>G:G=P
	return F,G,K
def d(LAM,kappa=KAPPA,beta_lo=H,grid=400):
	Y=kappa;U=grid;I=beta_lo;M=sqrt_up(LAM);D=[J]*12
	for G in A(0,11):
		Z={k:sqrt_up(Fr(1)/Y**(k-G))for k in A(G+1,12)};O={}
		for H in A(11,G,-1):O[H]=Z[H]+Fr(3,2)+E((O[A]+1)/2 for A in A(H+1,12))
		V={}
		for H in A(11,G,-1):d=E(O[A]for A in A(H+1,12))+1;V[H]=E(V[A]for A in A(H+1,12))/2+d/Q32+M*(Z[H]+Fr(1,64)+Fr(1,1<<20))
		e=E(O.values())+1;f=E(V.values());D[G]=(e+f)/Q32
	D[11]=Fr(1,Q32);W={}
	for G in A(2,12):X=1+(1<<19)*D[G];W[G]=((1<<20)*X*(1+Fr(1,1<<37))+2*X*X+2)/Q38
	a=E(W.values());C=M*sqrt_up(1+Fr(1,4096)+Fr(1,1<<38))+D[1]
	def P(b):A=sqrt_up(Fr(1)/b);B=Fr(1)/b;E=(2*C*A+1)**2+937*(2*C*A)+1+(1+Fr(1,1<<26))*C*C*B+2*C*A+4+(1<<20)*C;return E/Q38+D[1]*(6*C+2+D[1])
	def b(Fv,b0):A=(1+Fv)*(1+Fr(1,Q38))+b0/4096;return a+(M*sqrt_up(A)+D[0])**2-1
	F=[I*(Fr(1)/I)**Fr(j,U)if False else J for j in A(0)];F=[];h=sqrt_dn(sqrt_dn(sqrt_dn(Fr(1)/I)));g=Fr(1)+Fr(1,U//8 if U>=16 else 2);Q=I
	while Q<1:F.append(Q);Q=Q*g
	F.append(Fr(1));K=J
	for R in A(N(F)-1):
		c=b(P(F[R]),L(Fr(1),F[R+1]/Y))
		if K is J or c>K[0]:K=c,F[R],F[R+1]
	S={};S['R1K1_max_log2']=B.log2((1<<39)*(2*C+1)*Q32*D[1]+(1<<33));S[T]=B.log2((1+P(I))*(Fr(4096)/I+Fr(1,1<<26))+1);S['dr1_max_log2']=B.log2(Q38*(2*C+2)+Q38);return dict(sL=M,z=D,e=W,D1=a,W=C,Fb=P,eps_at=b,eps=K[0],eps_where=K[1:],F_lo=P(I),npts=N(F),WR=S)
G6=[[753729,251243,251243,251244,251242,251244],[251243,1163069,-455913,-455913,517661,517662],[251243,-455913,1163069,189495,189494,-455913],[251244,-455913,189495,1542060,127747,-394165],[251242,517661,189494,127747,1542058,455914],[251244,517662,-455913,-394165,455914,1542060]]
P=Fr(74,100)
I=3<<12
def U(M):
	M=[[Fr(b)for b in a]for a in M];D=Fr(1)
	for B in A(N(M)):
		D*=M[B][B]
		for C in A(B+1,N(M)):E=M[C][B]/M[B][B];M[C]=[j-E*k for(j,k)in zip(M[C],M[B])]
	return D
def e():
	H=Fr(1<<124)*Fr(1764)**6*Fr(1000)**6*U(G6)**2;I=Fr(1764)*Fr(56880,100)**2*196608-1;F=[]
	for E in A(12):
		G=C**E
		for J in A(E,12):G*=D(C,P**(J-E)*I)
		F.append(B.log2(G/H))
	return F
def exact_nodes_bound(UB,eta_max,lem):
	C=1+UB[K]+Fr(12,1<<38)
	def F(b):A=2*sqrt_up(UB['sL']**2*(C/b+Fr(1,4096)+C/(1<<38)))+2*eta_max+Fr(1,1<<31);return A.numerator//A.denominator+1
	D=0
	for B in A(12):
		if lem[B]>0:continue
		E=1
		for G in A(B+1,12):E*=F(P**(G-B))
		D+=E
	return D
def f(UB,eta_max,lem):
	C=1+UB[K]+Fr(12,1<<38)
	def F(b):A=2*sqrt_up(UB['sL']**2*(C/b+Fr(1,4096)+C/(1<<38)))+2*eta_max+Fr(1,1<<31);return A.numerator//A.denominator+1
	D=0
	for B in A(12):
		if lem[B]>0:continue
		for G in A(B+1,12):
			E=1
			for H in A(G+1,12):E*=F(P**(H-B))
			D+=E
	return D
W=99
X=R(2.**31.2)
def direct_q12_lut_certificate():
	D=(Fr(1)+Fr(1,1<<26))**2;E=0
	for F in(0,1,2,3,7,15,63,1024,(1<<20)-1,(1<<36)-1,(1<<37)-1):
		C=F+1;A=B.isqrt(C<<52)+1
		if A*A<=C<<52 or Fr(A*A,C<<52)>D:E+=1
	return E,D
def g(eps_or):
	D=sqrt_up(1+eps_or);F=[sqrt_up(Fr(1)/KAPPA**a)for a in A(12)];B=[Fr(0)]*12
	for C in A(11,-1,-1):B[C]=D+Fr(1,2)+E(Fr(1,2)*B[A]*F[A-C]for A in A(C+1,12))
	G=Fr(1)/sqrt_dn(H);return[B[b]*G for b in A(12)]
def h(eps_or):
	H=eps_or;F=6;B=[[Fr(G6[a][k])for k in A(F)]+[Fr(R(a==l))for l in A(F)]for a in A(F)]
	for C in A(F):
		L=next(A for A in A(C,F)if B[A][C]!=0);B[C],B[L]=B[L],B[C];T=B[C][C];B[C]=[j/T for j in B[C]]
		for G in A(F):
			if G!=C and B[G][C]:U=B[G][C];B[G]=[B[G][m]-U*B[C][m]for m in A(2*F)]
	M=D(B[A][F+A]for A in A(F));V=sqrt_up((1+H)*Fr(TMAX,1000)*M);W=sqrt_up(1+H);N=[sqrt_up(Fr(1)/KAPPA**b)for b in A(12)];I=[Fr(0)]*12
	for J in A(11,-1,-1):I[J]=W+Fr(1,2)+E(Fr(1,2)*I[A]*N[A-J]for A in A(J+1,12))
	K=Fr(1)+H
	for O in A(1,12):
		P=Fr(0)
		for X in A(O):Q=E(Fr(1,2)*I[A]*N[A-X]for A in A(O,12));P+=Q*Q
		S=1+H+P
		if S>K:K=S
	Y=(sqrt_up(K)+sqrt_up(Fr(5,4)))*sqrt_up(Fr(TMAX,1000)*M);return V,Y
def exact_node_ops():A=lambda a,b:4*a*b+2*(a+b);B=lambda n:A(1,n)+9*n+10;C=8*(A(8,1)+27)+B(9);D=A(1,8)+27+A(9,9)+57+57+10;E=57+B(16)+100;return C+A(1,16)+4*D+3*E+20
def i():
	C=J;M=sqrt_up(Fr(1)/H)+Fr(3,2);D={}
	for G in A(11,1,-1):D[G]=M+E((D[A]+1)/2 for A in A(G+1,12))
	I=(E(D.values())+1)/Q32
	for B in A(1,65):
		for K in(Fr(4),Fr(1<<38)):
			for F in(Fr(0),I):
				L=K*F*(B-1)-2*B*F+1+Fr((3*B+1)*(B-1),4)
				if C is J or L<C[0]:C=L,K,F,B
	return C+(I,)
def main():j=time.time();Y,k,l=c();Z,a=direct_q12_lut_certificate();C=d(a);N=i();R=B.log2((1+C[K])*(Fr(4096)/H+Fr(1,1<<26))+1);C[M]['r2_k1_max_log2']=R;S=g(D(C[O],Q.EPS_OR_MAX));m=D(B.log2(A)for A in S);P=dict(lut7_buckets=l,lut7_outward_errs=Y,lut7_lambda=F(k),lut_direct_q12_errs=Z,lut_direct_q12_lambda=F(a),lut_direct_q12_in_bounds_ok=C[M][T]<37. and R<37.,uncond_xmax_max_log2=G(m,3),uncond_xmax_lt_2p22=D(S)<1<<22,beta_lo_log2=B.log2(H),z_max_log2=B.log2(D(C['z'][:11])),delta1_levels_11_to_2=F(C[K]),F_beta_lo=F(C['F_lo']),eps_or_bound=F(C[O]),word_widths=C[M],pad1_min_slack=F(N[0]));E=e();J=exact_nodes_bound(C,N[4],E);U=f(C,N[4],E);b=exact_node_ops();n,V=h(D(C[O],Q.EPS_OR_MAX));P.update(lemma_e_log2_margin_k0_k3=[G(o,3)for o in E[:4]],exact_levels_possible=[p for p in A(12)if E[p]<=0],exact_nodes_per_orbit_max=J,exact_nodes_log2=G(B.log2(J),3),exact_node_ops=b,exact_node_price_charged=I,exact_per_orbit_ir_log2=G(B.log2(J*I),3),fp_nodes_above_exact_max=U,exact_level_work_ir=J*I+U*W,exact_level_allowance_ir=X,swar21_leaf_coord_ub=G(F(n),2),swar21_max_coord_ub=G(F(V),2),swar21_zero_overflow_ok=V<1<<20);P['ok']=Y==0 and Z==0 and C[M][T]<37. and R<37. and D(S)<1<<22 and C[O]<=Q.EPS_OR_MAX and N[0]>0 and L(E[:3])>0 and 2*b<=I and J*I+U*W<=X and V<1<<20;P['seconds']=G(time.time()-j,3);return P
if __name__=='__main__':import json;S=main();print(json.dumps(S,default=str));sys.exit(0 if S['ok']else 1)
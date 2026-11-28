r=tuple
j=next
V=None
U=all
T=zip
M=round
N=max
K=float
I=int
H=len
D=sum
C=range
import sys
from fractions import Fraction as F
from itertools import product
from math import isqrt,log2
G=[[753729,251243,251243,251244,251242,251244],[251243,1163069,-455913,-455913,517661,517662],[251243,-455913,1163069,189495,189494,-455913],[251244,-455913,189495,1542060,127747,-394165],[251242,517661,189494,127747,1542058,455914],[251244,517662,-455913,-394165,455914,1542060]]
S=[[1,0,1,0,1,0],[0,0,-1,0,-1,1],[0,1,-1,0,0,1],[0,0,0,0,-1,0],[0,0,0,0,0,-1],[0,0,0,1,0,0]]
G_SC=196608
R_A_MAX=F(572,1)
W_A,s=1764,1000
DELTA=F(1,200)
EPS_OR_MAX=F(35,10000)
t=16
u,K_B1=2,256
QA_MAX,QB_MAX=F(17,20),F(17,20)
BLK=38
K_MIN_STRIDE=41
J=64
N_SC=1<<34
Y=1<<36
P=1<<14
k=405,377,352,269,269,269
l=537,500,468,357,357,357
v=r(2*A+1 for A in k)
w=r(2*A+1 for A in l)
Q=1
for W in v:Q*=W
Z=1
for W in w:Z*=W
m=4*Z
a=1<<49
y=(1<<37)*8
WHEEL_TABS=1<<31
STAGEB=1<<25
B=13
STAGEB_G_SP=64388
E=50
F_NODE_SP=5880868228802566
z=1<<60
A0=52.6
A9,AA,AB=63.42,62.5,59.5
A=96
X=F(355,113)
def up(q):q=F(q);return F(-((-q.numerator<<A)//q.denominator),1<<A)
def dn(q):q=F(q);return F((q.numerator<<A)//q.denominator,1<<A)
def sqrt_up(q):q=F(q);assert q>=0;C=(q.numerator<<2*A)//q.denominator+1;B=isqrt(C);B+=B*B<C;return F(B,1<<A)
def sqrt_dn(q):q=F(q);return F(isqrt((q.numerator<<2*A)//q.denominator),1<<A)
def det(M):
	M=[[F(ab)for ab in aa]for aa in M];B=H(M);D=F(1)
	for A in C(B):
		E=j(B for B in C(A,B)if M[B][A]!=0)
		if E!=A:M[A],M[E]=M[E],M[A];D=-D
		D*=M[A][A]
		for G in C(A+1,B):I=M[G][A]/M[A][A];M[G]=[M[G][ac]-I*M[A][ac]for ac in C(B)]
	return D
def inv(M):
	D=H(M);A=[[F(M[aa][ad])for ad in C(D)]+[F(I(aa==ae))for ae in C(D)]for aa in C(D)]
	for B in C(D):
		G=j(C for C in C(B,D)if A[C][B]!=0);A[B],A[G]=A[G],A[B];J=A[B][B];A[B]=[ac/J for ac in A[B]]
		for E in C(D):
			if E!=B and A[E][B]:K=A[E][B];A[E]=[A[E][af]-K*A[B][af]for af in C(2*D)]
	return[ab[D:]for ab in A]
mm=lambda A,B:[[D(A[aa][C]*B[C][ab]for C in C(H(B)))for ab in C(H(B[0]))]for aa in C(H(A))]
L=lambda A:[list(aa)for aa in T(*A)]
def x(M):A=H(M);return N(F(D(B[D]*M[D][E]*B[E]for D in C(A)for E in C(A)),4)for B in product((1,-1),repeat=A))
b={6:up(X**3/6),5:up(8*X**2/15),3:up(4*X/3)}
class c:
	def __init__(A,M):A.n=H(M);A.h=sqrt_up(x(M));A.isd=up(1/sqrt_dn(det(M)))
	def vol(A,R2,pad=True):B=sqrt_up(R2)+(A.h if pad else 0);return up(b[A.n]*B**A.n*A.isd)
assert mm(L(S),mm(G,S))==G,'G not S-invariant'
A1=mm(S,S)
assert mm(S,A1)==[[I(ab==af)for af in C(6)]for ab in C(6)]
for O in([1,2,3,4,5,6],[0,0,0,7,-2,5]):d=O[3],-O[4],O[5];f=[D(S[ac][A]*O[A]for A in C(6))for ac in C(6)];assert(f[3],-f[4],f[5])==(d[1],d[2],d[0])
assert DELTA>=EPS_OR_MAX
T12_MAX=W_A*G_SC*R_A_MAX**2
assert(1+EPS_OR_MAX)*W_A*G_SC*F(56883,100)**2<=(1+DELTA)*T12_MAX
R=T12_MAX/W_A
B2=T12_MAX/s
E6=c(G)
def e(i):return[I(i==aa)for aa in C(6)]
def add(*A):return[D(aa)for aa in T(*A)]
def neg(v):return[-aa for aa in v]
A2={'s_x=0':([0,0,0,1,-1,1],[e(0),e(1),e(2),add(e(3),e(4)),add(e(4),e(5))]),'x0=x1':([0,0,0,1,1,0],[e(0),e(1),e(2),add(e(3),neg(e(4))),e(5)]),'x1=x2':([0,0,0,0,1,1],[e(0),e(1),e(2),e(3),add(e(4),neg(e(5)))]),'x2=x0':([0,0,0,1,0,-1],[e(0),e(1),e(2),add(e(3),e(5)),e(4)])}
n=[]
for(AC,(o,g))in A2.items():assert U(D(A*B for(A,B)in T(o,A))==0 for A in g);A3=j(e(A)for A in C(6)if abs(D(A*B for(A,B)in T(o,e(A))))==1);assert abs(det(L(g+[A3])))==1;p=L(g);n.append(c(mm(L(p),mm(G,p))))
A4=c([ad[:3]for ad in G[:3]])
def cnt_a_up(s):A=s*R;return up((E6.vol(A)+D(B.vol(A)for B in n)+A4.vol(A))/6)
A5=I(cnt_a_up(QA_MAX))+1
D_B1=-(--(-A5//u)//K_B1)
def CA(u):return I(cnt_a_up(min(u,QA_MAX)+DELTA))+1+D_B1+t
S5=[[F(G[aa][ae])-F(G[aa][0]*G[0][ae],G[0][0])for ae in C(1,6)]for aa in C(1,6)]
A6=up(1/sqrt_dn(det(S5)))
def SB(s):A=sqrt_up(s*B2)+E6.h;B=up(b[6]*A**6*E6.isd);C=up(b[5]*A**5*A6);return up(B/(2*J)+C)
h=1<<24
def A7(k):
	B=1<<k-1
	if CA(F(0))>B:return
	A,C=0,h
	if CA(QA_MAX)<=B:return QA_MAX
	while C-A>1:
		D=(A+C)//2
		if CA(F(D,h)*QA_MAX)<=B:A=D
		else:C=D
	return F(A,h)*QA_MAX
def main():
	i=CA(QA_MAX);b=N(K_MIN_STRIDE,(i-1).bit_length());O=SB(QB_MAX);L=J*(F(1<<K_MIN_STRIDE)*O);j=[]
	for c in C(K_MIN_STRIDE+1,b+1):d=A7(c);n=O if d is V else SB(min(QB_MAX,F(1)-d));j.append((c,d,n));L+=J*F(1<<c-1)*n
	L=up(L);o=J<<b;H=I(cnt_a_up(QA_MAX+DELTA))+1;S=inv(G);A1=[sqrt_up((QA_MAX+DELTA)*R*S[aa][aa])for aa in C(6)];A2=[sqrt_up(QB_MAX*B2*S[ab][ab])for ab in C(6)];A3=[sqrt_up((F(1)+DELTA)*R*S[ac][ac])for ac in C(6)];A4=[sqrt_up((F(1)+DELTA)*B2*S[ad][ad])for ad in C(6)];p=U(A3[A]<F(k[A]+1)for A in C(6))and U(A4[A]<F(l[A]+1)for A in C(6))and 16*Q<1<<64 and Z<1<<59;q=H<a and H<1<<49;T=256*(1<<20);r=N_SC-1;A5=I(J*O//r)+1+8*J;s=H+16*H;t=8*(1<<B);u=STAGEB_G_SP*8+(B+1)*t*8+D((1<<B-A)*(8*(1<<A)+B+4)for A in C(B+1))*8+8*N(STAGEB_G_SP,t)*8+24*(1<<B);v=7*(1<<E);e=F_NODE_SP*8+(E+1)*v*8+D((1<<E-A)*(7*(1<<A)+E+4)for A in C(E+1))*8+8*N(F_NODE_SP,v)*8+24*(1<<E)+N_SC*P;W=up(L/r+2*o+a+(Y+1)*P+T+STAGEB+s);X=-(-W.numerator//W.denominator>>BLK)+1;f=16*Q+8*m+N_SC*(T+8)+32*H+(1<<14)+y+WHEEL_TABS;A6=F(2)**56 if A0<=56 else V;g=z+A6+N_SC*(P+16)+1072*O;w=up(L+N_SC*(2*o+(Y+1)*P+T+STAGEB+a+s)+Y*g+f+e+F(2)**64+F(2)**64+F(2)**60);A=lambda x:K(log2(F(x).numerator)-log2(F(x).denominator))+5e-07;h=dict(D_B1_log2=A(N(D_B1,1)),T12_MAX_log2=A(T12_MAX),corner_h_a=K(E6.h/sqrt_dn(R)),corner_h_b=K(E6.h/sqrt_dn(B2)),CA1_log2=A(i),kappa_max=b,Ea_plus_log2=A(H),SB1_log2=A(O),max_rows_per_Cs=A5,epoch_buf_bytes=T,stageb_g_L=B,stageb_g_bytes=u,stageb_cap_bytes=STAGEB,f_node_L=E,f_node_log2=A(e),levels=[(af,V if ae is V else M(K(ae),7),M(A(ag),4))for(af,ae,ag)in j],score_rows_total_log2=A(L),score_node_log2=A(W),score_node_blocks_log2=A(X),emitter_log2=A(f),list_node_log2=A(g),machine_total_log2=A(w),box_ok=p,sink_ok=q,amax=[M(K(ah),2)for ah in A1],bmax=[M(K(ai),2)for ai in A2]);x=dict(delta_ge_eps_or=DELTA>=EPS_OR_MAX,box_ok=p,sink_ok=q,stageb_g_fits_ok=u<=STAGEB,f_node_lt_2_64=e<1<<64,score_node_lt_2_64=W<F(2)**64,row_block_fits_26_bits=X<1<<26,row_block_fits_29_bits=X<1<<29,max_rec_addr_lt_2_64=(X-1<<BLK)+H+1<1<<64,a_tab_16b_ptr_ok=16*Q+8*m<1<<64,emitter_lt_2_64=f<2**64,list_node_lt_2_64=g<2**64,machine_lt_2_100=w<F(2)**100);h['checks']=x;h['ok']=U(x.values());return h
if __name__=='__main__':import json,time as q;A8=q.time();i=main();i['seconds']=M(q.time()-A8,3);print(json.dumps(i,indent=1,default=str));sys.exit(0 if i['ok']else 1)
def ko(p):import io;return io.StringIO('\n'.join(A if A==''or A[:2]=='//'or A[-1]==':'or A[0]==' 'else' '+A.replace(',',', ')for A in open(p,encoding='ascii').read().split('\n')))
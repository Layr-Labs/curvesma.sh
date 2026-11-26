import sys
from fractions import Fraction as F
from itertools import product
from math import isqrt,log2
G=[[753729,251243,251243,251244,251242,251244],[251243,1163069,-455913,-455913,517661,517662],[251243,-455913,1163069,189495,189494,-455913],[251244,-455913,189495,1542060,127747,-394165],[251242,517661,189494,127747,1542058,455914],[251244,517662,-455913,-394165,455914,1542060]]
S=[[1,0,1,0,1,0],[0,0,-1,0,-1,1],[0,1,-1,0,0,1],[0,0,0,0,-1,0],[0,0,0,0,0,-1],[0,0,0,1,0,0]]
G_SC=196608
R_A_MAX=F(572,1)
W_A,W_B=1764,1000
DELTA=F(1,200)
EPS_OR_MAX=F(23389,10000000)
GUARD=16
B_A,K_B1=2,256
QA_MAX,QB_MAX=F(17,20),F(17,20)
BLK=38
K_MIN_STRIDE=41
W_SB=64
N_SC=1<<34
N_LIST=1<<36
PAGE=1<<14
A_TAB_AMAX_FULL=405,377,352,269,269,269
B_TAB_BMAX_FULL=537,500,468,357,357,357
A_TAB_MA=tuple(2*A+1 for A in A_TAB_AMAX_FULL)
B_TAB_MB=tuple(2*A+1 for A in B_TAB_BMAX_FULL)
A_TAB_ENTRIES=1
for _m in A_TAB_MA:A_TAB_ENTRIES*=_m
B_TAB_PLANE_ENTRIES=1
for _m in B_TAB_MB:B_TAB_PLANE_ENTRIES*=_m
B_TAB_ENTRIES=4*B_TAB_PLANE_ENTRIES
SINK=1<<49
LUT_DIRECT_Q12=(1<<37)*8
WHEEL_TABS=1<<31
STAGEB=1<<25
STAGEB_G_L=13
STAGEB_G_SP=64388
F_NODE_L=50
F_NODE_SP=5880868228802566
LIST_POOL=1<<60
LIST_SCRATCH_LOG2=52.6
LA_LOG2,BOOT_LOG2,SPLIT_LOG2=63.42,62.5,59.5
PREC=96
PI_UP=F(355,113)
def up(q):q=F(q);return F(-((-q.numerator<<PREC)//q.denominator),1<<PREC)
def dn(q):q=F(q);return F((q.numerator<<PREC)//q.denominator,1<<PREC)
def sqrt_up(q):q=F(q);assert q>=0;B=(q.numerator<<2*PREC)//q.denominator+1;A=isqrt(B);A+=A*A<B;return F(A,1<<PREC)
def sqrt_dn(q):q=F(q);return F(isqrt((q.numerator<<2*PREC)//q.denominator),1<<PREC)
def det(M):
	M=[[F(A)for A in A]for A in M];B=len(M);C=F(1)
	for A in range(B):
		D=next(B for B in range(A,B)if M[B][A]!=0)
		if D!=A:M[A],M[D]=M[D],M[A];C=-C
		C*=M[A][A]
		for E in range(A+1,B):G=M[E][A]/M[A][A];M[E]=[M[E][B]-G*M[A][B]for B in range(B)]
	return C
def inv(M):
	C=len(M);A=[[F(M[A][B])for B in range(C)]+[F(int(A==B))for B in range(C)]for A in range(C)]
	for B in range(C):
		E=next(C for C in range(B,C)if A[C][B]!=0);A[B],A[E]=A[E],A[B];G=A[B][B];A[B]=[A/G for A in A[B]]
		for D in range(C):
			if D!=B and A[D][B]:H=A[D][B];A[D]=[A[D][C]-H*A[B][C]for C in range(2*C)]
	return[A[C:]for A in A]
mm=lambda A,B:[[sum(A[D][C]*B[C][E]for C in range(len(B)))for E in range(len(B[0]))]for D in range(len(A))]
tr=lambda A:[list(A)for A in zip(*A)]
def corner_h2(M):A=len(M);return max(F(sum(B[C]*M[C][D]*B[D]for C in range(A)for D in range(A)),4)for B in product((1,-1),repeat=A))
OMEGA_UP={6:up(PI_UP**3/6),5:up(8*PI_UP**2/15),3:up(4*PI_UP/3)}
class Ell:
	def __init__(A,M):A.n=len(M);A.h=sqrt_up(corner_h2(M));A.isd=up(1/sqrt_dn(det(M)))
	def vol(A,R2,pad=True):B=sqrt_up(R2)+(A.h if pad else 0);return up(OMEGA_UP[A.n]*B**A.n*A.isd)
assert mm(tr(S),mm(G,S))==G,'G not S-invariant'
S2=mm(S,S)
assert mm(S,S2)==[[int(A==B)for B in range(6)]for A in range(6)]
for v in([1,2,3,4,5,6],[0,0,0,7,-2,5]):x=v[3],-v[4],v[5];w=[sum(S[B][A]*v[A]for A in range(6))for B in range(6)];assert(w[3],-w[4],w[5])==(x[1],x[2],x[0])
assert DELTA>=EPS_OR_MAX
T12_MAX=W_A*G_SC*R_A_MAX**2
A2=T12_MAX/W_A
B2=T12_MAX/W_B
E6=Ell(G)
def e(i):return[int(i==A)for A in range(6)]
def add(*A):return[sum(A)for A in zip(*A)]
def neg(v):return[-A for A in v]
WALLS={'s_x=0':([0,0,0,1,-1,1],[e(0),e(1),e(2),add(e(3),e(4)),add(e(4),e(5))]),'x0=x1':([0,0,0,1,1,0],[e(0),e(1),e(2),add(e(3),neg(e(4))),e(5)]),'x1=x2':([0,0,0,0,1,1],[e(0),e(1),e(2),e(3),add(e(4),neg(e(5)))]),'x2=x0':([0,0,0,1,0,-1],[e(0),e(1),e(2),add(e(3),e(5)),e(4)])}
WALL_ELL=[]
for(name,(nv,Bs))in WALLS.items():assert all(sum(A*B for(A,B)in zip(nv,A))==0 for A in Bs);comp=next(e(A)for A in range(6)if abs(sum(A*B for(A,B)in zip(nv,e(A))))==1);assert abs(det(tr(Bs+[comp])))==1;Bm=tr(Bs);WALL_ELL.append(Ell(mm(tr(Bm),mm(G,Bm))))
L0_ELL=Ell([A[:3]for A in G[:3]])
def cnt_a_up(s):A=s*A2;return up((E6.vol(A)+sum(B.vol(A)for B in WALL_ELL)+L0_ELL.vol(A))/6)
NA_REAL_UP=int(cnt_a_up(QA_MAX))+1
D_B1=-(--(-NA_REAL_UP//B_A)//K_B1)
def CA(u):return int(cnt_a_up(min(u,QA_MAX)+DELTA))+1+D_B1+GUARD
S5=[[F(G[A][B])-F(G[A][0]*G[0][B],G[0][0])for B in range(1,6)]for A in range(1,6)]
isd5=up(1/sqrt_dn(det(S5)))
def SB(s):A=sqrt_up(s*B2)+E6.h;B=up(OMEGA_UP[6]*A**6*E6.isd);C=up(OMEGA_UP[5]*A**5*isd5);return up(B/(2*W_SB)+C)
GRID=1<<24
def t_level(k):
	B=1<<k-1
	if CA(F(0))>B:return
	A,C=0,GRID
	if CA(QA_MAX)<=B:return QA_MAX
	while C-A>1:
		D=(A+C)//2
		if CA(F(D,GRID)*QA_MAX)<=B:A=D
		else:C=D
	return F(A,GRID)*QA_MAX
def main():
	R=CA(QA_MAX);K=max(K_MIN_STRIDE,(R-1).bit_length());D=SB(QB_MAX);C=W_SB*(F(1<<K_MIN_STRIDE)*D);S=[]
	for L in range(K_MIN_STRIDE+1,K+1):M=t_level(L);T=D if M is None else SB(min(QB_MAX,F(1)-M));S.append((L,M,T));C+=W_SB*F(1<<L-1)*T
	C=up(C);U=W_SB<<K;B=int(cnt_a_up(QA_MAX+DELTA))+1;E=inv(G);e=[sqrt_up((QA_MAX+DELTA)*A2*E[A][A])for A in range(6)];f=[sqrt_up(QB_MAX*B2*E[A][A])for A in range(6)];g=[sqrt_up((F(1)+DELTA)*A2*E[A][A])for A in range(6)];h=[sqrt_up((F(1)+DELTA)*B2*E[A][A])for A in range(6)];V=all(g[A]<F(A_TAB_AMAX_FULL[A]+1)for A in range(6))and all(h[A]<F(B_TAB_BMAX_FULL[A]+1)for A in range(6))and 16*A_TAB_ENTRIES<1<<64 and B_TAB_PLANE_ENTRIES<1<<59;W=B<SINK and B<1<<49;H=256*(1<<20);X=N_SC-1;i=int(W_SB*D//X)+1+8*W_SB;Y=B+16*B;Z=8*(1<<STAGEB_G_L);a=STAGEB_G_SP*8+(STAGEB_G_L+1)*Z*8+sum((1<<STAGEB_G_L-A)*(8*(1<<A)+STAGEB_G_L+4)for A in range(STAGEB_G_L+1))*8+8*max(STAGEB_G_SP,Z)*8+24*(1<<STAGEB_G_L);b=7*(1<<F_NODE_L);N=F_NODE_SP*8+(F_NODE_L+1)*b*8+sum((1<<F_NODE_L-A)*(7*(1<<A)+F_NODE_L+4)for A in range(F_NODE_L+1))*8+8*max(F_NODE_SP,b)*8+24*(1<<F_NODE_L)+N_SC*PAGE;I=up(C/X+2*U+SINK+(N_LIST+1)*PAGE+H+STAGEB+Y);J=-(-I.numerator//I.denominator>>BLK)+1;O=16*A_TAB_ENTRIES+8*B_TAB_ENTRIES+N_SC*(H+8)+32*B+(1<<14)+LUT_DIRECT_Q12+WHEEL_TABS;j=F(2)**56 if LIST_SCRATCH_LOG2<=56 else None;P=LIST_POOL+j+N_SC*(PAGE+16)+1072*D;c=up(C+N_SC*(2*U+(N_LIST+1)*PAGE+H+STAGEB+SINK+Y)+N_LIST*P+O+N+F(2)**64+F(2)**64+F(2)**60);A=lambda x:float(log2(F(x).numerator)-log2(F(x).denominator))+5e-07;Q=dict(D_B1_log2=A(max(D_B1,1)),T12_MAX_log2=A(T12_MAX),corner_h_a=float(E6.h/sqrt_dn(A2)),corner_h_b=float(E6.h/sqrt_dn(B2)),CA1_log2=A(R),kappa_max=K,Ea_plus_log2=A(B),SB1_log2=A(D),max_rows_per_Cs=i,epoch_buf_bytes=H,stageb_g_L=STAGEB_G_L,stageb_g_bytes=a,stageb_cap_bytes=STAGEB,f_node_L=F_NODE_L,f_node_log2=A(N),levels=[(C,None if B is None else round(float(B),7),round(A(D),4))for(C,B,D)in S],score_rows_total_log2=A(C),score_node_log2=A(I),score_node_blocks_log2=A(J),emitter_log2=A(O),list_node_log2=A(P),machine_total_log2=A(c),box_ok=V,sink_ok=W,amax=[round(float(A),2)for A in e],bmax=[round(float(A),2)for A in f]);d=dict(delta_ge_eps_or=DELTA>=EPS_OR_MAX,box_ok=V,sink_ok=W,stageb_g_fits_ok=a<=STAGEB,f_node_lt_2_64=N<1<<64,score_node_lt_2_64=I<F(2)**64,row_block_fits_26_bits=J<1<<26,row_block_fits_29_bits=J<1<<29,max_rec_addr_lt_2_64=(J-1<<BLK)+B+1<1<<64,a_tab_16b_ptr_ok=16*A_TAB_ENTRIES+8*B_TAB_ENTRIES<1<<64,emitter_lt_2_64=O<2**64,list_node_lt_2_64=P<2**64,machine_lt_2_100=c<F(2)**100);Q['checks']=d;Q['ok']=all(d.values());return Q
if __name__=='__main__':import json,time;t0=time.time();o=main();o['seconds']=round(time.time()-t0,3);print(json.dumps(o,indent=1,default=str));sys.exit(0 if o['ok']else 1)
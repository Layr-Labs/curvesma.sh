"""stor_bound.py -- outward, input-independent storage bound for the v26 direct live Score Nodes (stdlib only, < 1 s).

Every quantity below is an exact rational or an outward-rounded dyadic rational; no float enters a bound.
Floats are used only to print log2 values (the printed log2 is rounded up in the 6th decimal).

Geometry (exact):  E12 = {(a,b) in Z^12 : 1764 a^T G a + 1000 b^T G b <= T12},  G = G_sym (integer, S-invariant).
  Q(a) = 1764 a^T G a / T12,  Q(b) = 1000 b^T G b / T12   (E12 boundary = 1).
  T12 <= T12_MAX = 1764 * G_SC * R_A_MAX^2 with R_A_MAX = 572.0 (single-ball T=1 production R_a <= 566.0).  Every count used
  below is nondecreasing in T12, so the bound holds for every T12 <= T12_MAX.

Index sets:
  K_a   = {a : s_x >= 0, x_0 >= x_1, x_0 >= x_2},  x = (a_3, -a_4, a_5), s_x = x_0 + x_1 + x_2  (closed cone)
  E_a+  = K_a cap Z^6 cap {0 < Q(a) <= 1 + DELTA}, ranked 1..|E_a+| by (Q(a), fixed tie order); ranks <= |E_a| are
          the real a-set {Q <= 1}.  a_tab[x] = (ai, &B_TAB[k][0]) (16 B) for x = pack(a), a in E_a+, and 0 otherwise.
  E_b   = canonical b (first nonzero of (b_0, b_5, b_1, b_2, b_3, b_4) positive) with 0 < Q(b) <= 1; one row each.
  Row capacity (exact, uint8_t 1 B/cell): cap(sb) = #{a in E_a+ : Q(a) <= u_sb + DELTA} + D_B1 + 16, u_sb = max_{b in sb cap E_b} (1 - Q(b)),
          D_B1 = ceil(ceil(|E_a| / B_a) / K_B1) = Band-1 prefix-checkpoint spacing (Band 1 applies every list entry
          up to the first checkpoint >= N_a(b), i.e. up to D_B1 cells past the row's logical end).
  Stride 2^k_s, k_s = max(K_MIN_STRIDE, ceil(log2(cap))) with explicit minimum stride K_MIN_STRIDE = 41 (8 contiguous
  2^BLK = 2^38 B blocks, uint8_t: 1 B/cell); every sub-band therefore belongs to one of the 8 stride classes
  k_s in {41..48}, and a sub-band (64 b_0 values b_0 >> 6 = sb on one (b_1..b_5) chord) occupies 64 * 2^k_s bytes.

Lemmas used (all elementary):
  L1 (cube packing).  For convex K in R^n: #(K cap Z^n) <= Vol(K + [-1/2,1/2]^n).  For K = {x^T M x <= R^2},
      K + cube is inside {x^T M x <= (R + h)^2}, h = max over the 2^n cube corners c of sqrt(c^T M c) (triangle
      inequality for the norm sqrt(x^T M x)); volume = omega_n (R + h)^n / sqrt(det M).
  L2 (cone).  g in {+-I, +-S, +-S^2} are Q-preserving unimodular maps; S cycles (x_0, x_1, x_2).  The six cones gK_a
      overlap only on the 4 walls H: s_x = 0, x_i = x_j.  With m(x) = #{g : x in gK_a}: m(x) = #argmax x (s_x > 0),
      6 #(K_a cap E) = sum_x m(x) <= #E + sum_H #(H cap E) + #(L_0 cap E), L_0 = {a_3 = a_4 = a_5 = 0}.  Each
      H cap Z^6 is a rank-5 lattice with an explicit unimodular basis, L_0 cap Z^6 = Z^3 (L1 in 5 and 3 dimensions).
  L3 (sub-bands).  A sub-band active at level s (contains some b in E_b with Q(b) <= s) is a point
      y = (64 sb + 31.5, l) of 64Z x Z^5 (+ shift) in (E_b(s) + [-31.5, 31.5] e_0) cap {y_0 >= 31.5}.  Translates
      y + [-32,32] x [-1/2,1/2]^5 are disjoint and lie in W = (E(r) + [-63.5, 63.5] e_0) cap {y_0 >= -1/2} with
      r = sqrt(s) + h_b.  Pairing l with -l (chord centre c(-l) = -c(l)) gives Vol(W) <= Vol6(E(r))/2 + 64 Vol5(P(r)),
      P = projection of E(r) along e_0 (Schur complement S5, det S5 = det G / G_00).  So
      SB(s) <= Vol6(E(r))/128 + Vol5(P(r)).  This includes the worst-case 64-alignment waste at both chord ends.
  L4 (levels).  kappa(u) = max(K_MIN_STRIDE, ceil(log2(CA(u)))) with K_MIN_STRIDE = 41 (1 B/cell, BLK = 38 addressing),
      CA(u) >= cap, CA nondecreasing.  For k >= K_MIN_STRIDE + 1 pick t_k with CA(t_k) <= 2^(k-1) (largest dyadic grid point);
      then 1[k_s >= k] <= 1[u_sb >= t_k] and 2^k_s = 2^K_MIN_STRIDE + sum_{k=K_MIN_STRIDE+1}^{K} 2^(k-1) 1[k_s >= k].
      Summing over sub-bands: TOTAL_ROWS <= 64 (2^K_MIN_STRIDE SB(QB_MAX) + sum_{k=42}^{K} 2^(k-1) SB(min(QB_MAX, 1 - t_k))).
  L5 (greedy).  Least-loaded assignment of sub-bands to N_SC nodes: each node <= TOTAL/N_SC + largest sub-band.
"""
import sys
from fractions import Fraction as F
from itertools import product
from math import isqrt,log2
G= [[753729,251243,251243,251244,251242,251244],
[251243,1163069,-455913,-455913,517661,517662],
[251243,-455913,1163069,189495,189494,-455913],
[251244,-455913,189495,1542060,127747,-394165],
[251242,517661,189494,127747,1542058,455914],
[251244,517662,-455913,-394165,455914,1542060]]
S= [[1,0,1,0,1,0],[0,0,-1,0,-1,1],[0,1,-1,0,0,1],
[0,0,0,0,-1,0],[0,0,0,0,0,-1],[0,0,0,1,0,0]]
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
A_TAB_AMAX_FULL= (405,377,352,269,269,269)
B_TAB_BMAX_FULL= (537,500,468,357,357,357)
A_TAB_MA=tuple(2*m+1 for m in A_TAB_AMAX_FULL)
B_TAB_MB=tuple(2*m+1 for m in B_TAB_BMAX_FULL)
A_TAB_ENTRIES=1
for _m in A_TAB_MA:A_TAB_ENTRIES*=_m
B_TAB_PLANE_ENTRIES=1
for _m in B_TAB_MB:B_TAB_PLANE_ENTRIES*=_m
B_TAB_ENTRIES=4*B_TAB_PLANE_ENTRIES
SINK=1<<49
LUT_DIRECT_Q12= (1<<37) *8
WHEEL_TABS=1<<31
STAGEB=1<<25
STAGEB_G_L=13
STAGEB_G_SP=64388
F_NODE_L=50
F_NODE_SP=5880868228802566
LIST_POOL=1<<60
LIST_SCRATCH_LOG2=52.6
LA_LOG2,BOOT_LOG2,SPLIT_LOG2=63.42,62.50,59.50
PREC=96
PI_UP=F(355,113)
def up(q):
 q=F(q); return F(-((-q.numerator<<PREC) //q.denominator),1<<PREC)
def dn(q):
 q=F(q); return F((q.numerator<<PREC) //q.denominator,1<<PREC)
def sqrt_up(q):
 q=F(q); assert q>=0
 n= (q.numerator<< (2*PREC)) //q.denominator+1
 r=isqrt(n);r+= (r*r<n)
 return F(r,1<<PREC)
def sqrt_dn(q):
 q=F(q); return F(isqrt((q.numerator<< (2*PREC)) //q.denominator),1<<PREC)
def det(M):
 M= [[F(x) for x in r] for r in M];n=len(M);d=F(1)
 for c in range(n):
  p=next(r for r in range(c,n) if M[r][c] !=0)
  if p!=c:M[c],M[p] =M[p],M[c];d= -d
  d*=M[c][c]
  for r in range(c+1,n):
   f=M[r][c] /M[c][c]
   M[r] = [M[r][j] -f*M[c][j] for j in range(n)]
 return d
def inv(M):
 n=len(M);A= [[F(M[i][j]) for j in range(n)] + [F(int(i==j)) for j in range(n)] for i in range(n)]
 for c in range(n):
  p=next(r for r in range(c,n) if A[r][c] !=0);A[c],A[p] =A[p],A[c]
  pv=A[c][c];A[c] = [x/pv for x in A[c]]
  for r in range(n):
   if r!=c and A[r][c]:
    f=A[r][c];A[r] = [A[r][j] -f*A[c][j] for j in range(2*n)]
 return [r[n:] for r in A]
mm= lambda A,B: [[sum(A[i][k] *B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
tr= lambda A: [list(r) for r in zip(*A)]
def corner_h2(M):
 n=len(M)
 return max(F(sum(s[i] *M[i][j] *s[j] for i in range(n) for j in range(n)),4) for s in product((1,-1),repeat=n))
OMEGA_UP= {6:up(PI_UP**3/6),5:up(8*PI_UP**2/15),3:up(4*PI_UP/3)}
class Ell:
 """Lattice-count / volume upper bounds for {x in Z^n : x^T M x <= R^2} (M integer, in G units)."""
 def __init__(self,M):
  self.n=len(M);self.h=sqrt_up(corner_h2(M));self.isd=up(1/sqrt_dn(det(M)))
 def vol(self,R2,pad=True):
  r=sqrt_up(R2) + (self.h if pad else 0)
  return up(OMEGA_UP[self.n] *r**self.n*self.isd)
assert mm(tr(S),mm(G,S)) ==G,"G not S-invariant"
S2=mm(S,S); assert mm(S,S2) == [[int(i==j) for j in range(6)] for i in range(6)]
for v in ([1,2,3,4,5,6],[0,0,0,7,-2,5]):
 x= (v[3],-v[4],v[5]);w= [sum(S[i][j] *v[j] for j in range(6)) for i in range(6)]
 assert (w[3],-w[4],w[5]) == (x[1],x[2],x[0])
assert DELTA>=EPS_OR_MAX
T12_MAX=W_A*G_SC*R_A_MAX**2
A2=T12_MAX/W_A
B2=T12_MAX/W_B
E6=Ell(G)
def e(i): return [int(i==j) for j in range(6)]
def add(*vs): return [sum(c) for c in zip(*vs)]
def neg(v): return [-x for x in v]
WALLS= {
 "s_x=0":  ([0,0,0,1,-1,1],[e(0),e(1),e(2),add(e(3),e(4)),add(e(4),e(5))]),
 "x0=x1":  ([0,0,0,1,1,0],[e(0),e(1),e(2),add(e(3),neg(e(4))),e(5)]),
 "x1=x2":  ([0,0,0,0,1,1],[e(0),e(1),e(2),e(3),add(e(4),neg(e(5)))]),
 "x2=x0":  ([0,0,0,1,0,-1],[e(0),e(1),e(2),add(e(3),e(5)),e(4)]),
}
WALL_ELL= []
for name,(nv,Bs) in WALLS.items():
 assert all(sum(a*b for a,b in zip(nv,v)) ==0 for v in Bs)
 comp=next(e(i) for i in range(6) if abs(sum(a*b for a,b in zip(nv,e(i)))) ==1)
 assert abs(det(tr(Bs+ [comp]))) ==1
 Bm=tr(Bs);WALL_ELL.append(Ell(mm(tr(Bm),mm(G,Bm))))
L0_ELL=Ell([row[:3] for row in G[:3]])
def cnt_a_up(s):
 R2=s*A2
 return up((E6.vol(R2) +sum(we.vol(R2) for we in WALL_ELL) +L0_ELL.vol(R2)) /6)
NA_REAL_UP=int(cnt_a_up(QA_MAX)) +1
D_B1= -(-(-(-NA_REAL_UP//B_A)) //K_B1)
def CA(u):
 return int(cnt_a_up(min(u,QA_MAX) +DELTA)) +1+D_B1+GUARD
S5= [[F(G[i][j]) -F(G[i][0] *G[0][j],G[0][0]) for j in range(1,6)] for i in range(1,6)]
isd5=up(1/sqrt_dn(det(S5)))
def SB(s):
 r=sqrt_up(s*B2) +E6.h
 v6=up(OMEGA_UP[6] *r**6*E6.isd)
 v5=up(OMEGA_UP[5] *r**5*isd5)
 return up(v6/ (2*W_SB) +v5)
GRID=1<<24
def t_level(k):
 lim=1<< (k-1)
 if CA(F(0)) >lim: return None
 lo,hi=0,GRID
 if CA(QA_MAX) <=lim: return QA_MAX
 while hi-lo>1:
  mid= (lo+hi) //2
  if CA(F(mid,GRID) *QA_MAX) <=lim:lo=mid
  else:hi=mid
 return F(lo,GRID) *QA_MAX
def main():
 ca1=CA(QA_MAX);K=max(K_MIN_STRIDE,(ca1-1).bit_length())
 SB1=SB(QB_MAX)
 total=W_SB* (F(1<<K_MIN_STRIDE) *SB1)
 levels= []
 for k in range(K_MIN_STRIDE+1,K+1):
  t=t_level(k)
  sbk=SB1 if t is None else SB(min(QB_MAX,F(1) -t))
  levels.append((k,t,sbk));total+=W_SB*F(1<< (k-1)) *sbk
 total=up(total)
 sb_max=W_SB<<K
 na_plus=int(cnt_a_up(QA_MAX+DELTA)) +1
 Gi=inv(G)
 amax= [sqrt_up((QA_MAX+DELTA) *A2*Gi[i][i]) for i in range(6)]
 bmax= [sqrt_up(QB_MAX*B2*Gi[i][i]) for i in range(6)]
 amax_full= [sqrt_up((F(1) +DELTA) *A2*Gi[i][i]) for i in range(6)]
 bmax_full= [sqrt_up((F(1) +DELTA) *B2*Gi[i][i]) for i in range(6)]
 box_ok= (all(amax_full[i] <F(A_TAB_AMAX_FULL[i] +1) for i in range(6))
and all(bmax_full[i] <F(B_TAB_BMAX_FULL[i] +1) for i in range(6))
and 16*A_TAB_ENTRIES< (1<<64) and B_TAB_PLANE_ENTRIES< (1<<59))
 sink_ok=na_plus<SINK and na_plus< (1<<49)
 EPOCH_BUF=256* (1<<20)
 N_SC_REAL=N_SC-1
 max_rows_per_Cs=int((W_SB*SB1) //N_SC_REAL) +1+8*W_SB
 tabs_c=na_plus+16*na_plus
 g_nroot=8* (1<<STAGEB_G_L)
 g_stageb_bytes= (STAGEB_G_SP*8+ (STAGEB_G_L+1) *g_nroot*8
+sum((1<< (STAGEB_G_L-j)) * (8* (1<<j) +STAGEB_G_L+4) for j in range(STAGEB_G_L+1)) *8
+8*max(STAGEB_G_SP,g_nroot) *8+24* (1<<STAGEB_G_L))
 f_nroot=7* (1<<F_NODE_L)
 f0_node= (F_NODE_SP*8+ (F_NODE_L+1) *f_nroot*8
+sum((1<< (F_NODE_L-j)) * (7* (1<<j) +F_NODE_L+4) for j in range(F_NODE_L+1)) *8
+8*max(F_NODE_SP,f_nroot) *8+24* (1<<F_NODE_L) +N_SC*PAGE)
 c_node=up(total/N_SC_REAL+2*sb_max+SINK+ (N_LIST+1) *PAGE+EPOCH_BUF+STAGEB+tabs_c)
 blocks= -((-c_node.numerator//c_node.denominator) >>BLK) +1
 e0=16*A_TAB_ENTRIES+8*B_TAB_ENTRIES+N_SC* (EPOCH_BUF+8) +32*na_plus+ (1<<14) +LUT_DIRECT_Q12+WHEEL_TABS
 lscr=F(2) **56 if LIST_SCRATCH_LOG2<=56 else None
 l_node=LIST_POOL+lscr+N_SC* (PAGE+16) +1072*SB1
 machine=up(total+N_SC* (2*sb_max+ (N_LIST+1) *PAGE+EPOCH_BUF+STAGEB+SINK+tabs_c)
+N_LIST*l_node+e0+f0_node+F(2) **64+F(2) **64+F(2) **60)
 L= lambda x:float(log2(F(x).numerator) -log2(F(x).denominator)) +5e-7
 out=dict(
D_B1_log2=L(max(D_B1,1)),T12_MAX_log2=L(T12_MAX),corner_h_a=float(E6.h/sqrt_dn(A2)),corner_h_b=float(E6.h/sqrt_dn(B2)),
CA1_log2=L(ca1),kappa_max=K,Ea_plus_log2=L(na_plus),SB1_log2=L(SB1),
max_rows_per_Cs=max_rows_per_Cs,epoch_buf_bytes=EPOCH_BUF,
stageb_g_L=STAGEB_G_L,stageb_g_bytes=g_stageb_bytes,stageb_cap_bytes=STAGEB,
f_node_L=F_NODE_L,f_node_log2=L(f0_node),
levels=[(k,None if t is None else round(float(t),7),round(L(sb),4)) for k,t,sb in levels],
score_rows_total_log2=L(total),score_node_log2=L(c_node),score_node_blocks_log2=L(blocks),
emitter_log2=L(e0),list_node_log2=L(l_node),machine_total_log2=L(machine),
box_ok=box_ok,sink_ok=sink_ok,amax=[round(float(x),2) for x in amax],bmax=[round(float(x),2) for x in bmax])
 checks=dict(
delta_ge_eps_or=DELTA>=EPS_OR_MAX,box_ok=box_ok,sink_ok=sink_ok,
stageb_g_fits_ok=(g_stageb_bytes<=STAGEB),f_node_lt_2_64=(f0_node< (1<<64)),
score_node_lt_2_64=c_node<F(2) **64,row_block_fits_26_bits=blocks< (1<<26),row_block_fits_29_bits=blocks< (1<<29),
max_rec_addr_lt_2_64=(((blocks-1) <<BLK) +na_plus+1) < (1<<64),
a_tab_16b_ptr_ok=(16*A_TAB_ENTRIES+8*B_TAB_ENTRIES< (1<<64)),
emitter_lt_2_64=e0<2**64,list_node_lt_2_64=l_node<2**64,machine_lt_2_100=machine<F(2) **100)
 out["checks"] =checks;out["ok"] =all(checks.values())
 return out
if __name__=="__main__":
 import json,time
 t0=time.time();o=main();o["seconds"] =round(time.time() -t0,3)
 print(json.dumps(o,indent=1,default=str))
 sys.exit(0 if o["ok"] else 1)

import functools,itertools,json,math,os,re,sys
from decimal import Decimal,ROUND_CEILING,getcontext
from fractions import Fraction as Fr
FAILED=[]
def check(name,ok,detail=""):
 print("%s %s%s"%(name,"ok"if ok else"FAIL",("  "+detail)if detail else""))
 if not ok:FAILED.append(name)
def lg(x):
 f=Fr(x)
 return math.log2(f.numerator)-math.log2(f.denominator)
def cdiv(a,b):return-((-a)//b)
def exp_interval(x,terms=60):
 assert abs(x)<=1
 n,d,S=x.numerator,x.denominator,1<<400
 lo=hi=0;tl=th=S
 for k in range(terms+1):
  if k:a,b=sorted((tl*n,th*n));tl,th=a//(d*k),cdiv(b,d*k)
  lo,hi=lo+tl,hi+th
 rem=cdiv(max(-tl,th)*abs(n)*3,d*(terms+1))
 return Fr(lo-rem,S),Fr(hi+rem,S)
def ln2_interval(bits=300):
 z,z2=Fr(1,3),Fr(1,9)
 tot,pw,k=Fr(0),z,0
 while True:
  tot+=pw/(2*k+1)
  rest=pw*z2/((2*k+3)*(1-z2))
  if rest<Fr(1,2**bits):return 2*tot,2*(tot+rest)
  pw*=z2;k+=1
LN2=ln2_interval()
@functools.lru_cache(maxsize=None)
def P2(c):
 n=c.numerator//c.denominator
 f=c-n
 lo,hi=exp_interval(f*LN2[0])[0],exp_interval(f*LN2[1])[1]
 sc=Fr(2)**n
 return lo*sc,hi*sc
def exp_neg_hi(m_val):
 if not m_val:return Fr(1)
 k_h,mv=0,m_val
 while mv>1:mv/=2;k_h+=1
 mlo=Fr((mv.numerator<<200)//mv.denominator,2**200)
 eh=exp_interval(-mlo,45)[1]
 for _ in range(k_h):
  eh=eh*eh
  eh=Fr(cdiv(eh.numerator<<300,eh.denominator),1<<300)
 return eh
p=2**256-2**32-977
r=int("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141",16)
C0=2**32+977
BETA_P=int("7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee",16)
LAM_R=int("5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72",16)
C=1000
NODES=2**66
M_G=1024
N_CLUST=NODES//M_G
D1,D2,D3,D4=41,41,64,23
F1,F2=1<<D1,1<<D2
N_REL,N_HOME=F1,F1*F2
ALLN=NODES+N_REL+N_HOME
MX=8
YPC=33
YROUNDS=63520
L_CHAIN=YROUNDS*YPC
Q2=MX*L_CHAIN
KEYBITS=64+D4
PART=64+D1
KBF=6310
K1,K2B=64,128
K=K1
HB=50
THRESH=Fr(39,100)
LEAVES=6
FWD_INC,FWD_FIX=Fr(2042,32),7
def fwd_count(q2):
 assert q2%32==0
 return FWD_FIX+(2042*q2)//32
BWD_INC,BWD_FIX=Fr(14422,16),7
def bwd_count(q2):
 assert q2%16==0
 return BWD_FIX+(14422*q2)//16
FIX_SITES_QUAD=dict(SHX=54,FDX=12,F5=8,FS=6,FADDW=8,BLK=28)
FIX_EXTRA=dict(SHX=3,FDX=3,F5=3,FS=4,FADDW=4,BLK=6)
EMIT=7
EMIT_INC=LEAVES*EMIT
KEY_W3=Fr(1,1)
KEY_INC,KEY_FIX=Fr(533,4),5
FLUSH=10
FLUSH1_X,FIN1_X=8,14
FLUSH2_X,FLUSH2_D=10,8
XEMIT_HOT=14
JPATH,SJOIN=3,1
INV=22549
INIT3_BKT,TINIT_ENTRY=4,3
FIN_DEST,FIN_SENT,FIN_HELD,FIN_FIX=15,4,2,2
FINH_DEST,FINH_SENT=16,5
H2S_ARR,H2S_ONE=6,1
H2S_FULL=1
H2S_EMPTY=6
H2S_ARR_FIX=16+H2S_ONE
def h2s_count(n):
 if n==K:return H2S_FULL+H2S_ARR*K
 if n==0:return H2S_EMPTY
 return 16+H2S_ONE*(n%2)+H2S_ARR*n
check("c1",h2s_count(K)==385 and max(h2s_count(n)-H2S_ARR*n for n in range(K))==H2S_ARR_FIX and H2S_EMPTY<=H2S_ARR_FIX)
check("c2",KEY_INC==8+52+KEY_W3+6*(5+EMIT)+Fr(1,4)and KEY_FIX==5)
check("c3",BWD_INC==Fr(14422,16)and Q2%32==0 and BWD_FIX==7 and bwd_count(Q2)==BWD_FIX+(14422*Q2)//16)
check("c4",FWD_INC==Fr(2042,32)and Q2%32==0 and FWD_FIX==7 and fwd_count(Q2)==FWD_FIX+(2042*Q2)//32)
check("c5",(3+3+2+EMIT)-(5+EMIT)==JPATH and SJOIN==1)
check("c6",EMIT==1+1+2+2+1 and XEMIT_HOT==2+4+4+4)
DONE=2
RCV,WAKE,WY,EXPECT,SWAP,CTX=1,5,20,39,4+4*15,0
DONE_C=DONE+RCV+WAKE
RCV1,FIN_CNT,H2S_DONE=4,6,4
EXPECTH=78
RDONE=1+H2S_DONE
DONE_C1=DONE+RCV1+WAKE
DLV=6
MSG_LINK=DONE_C+DLV+CTX
FINAL=FIN_DEST+FIN_SENT+DONE_C1+H2S_ARR_FIX+FIN_CNT
check("c7",FLUSH1_X<=FLUSH+DONE and FIN1_X<=FIN_DEST+FIN_SENT+DONE and FLUSH2_X<=FLUSH and FLUSH2_D<=FLUSH and FIN_HELD<=FINH_SENT and DONE==2 and(DONE_C,DONE_C1,MSG_LINK,FINAL,RDONE)==(8,11,14,53,5)and H2S_ARR_FIX==17)
check("c8",5<=WAKE and 9+2*5<=WY and(SWAP,CTX)==(64,0)and 3+5+5+26<=EXPECT and EXPECTH==2*EXPECT and FIN_FIX==2)
GEN_MUL=645410
GEN_MULX=GEN_MUL+4
YR_FIX,YR_BLK,YB=22556,42763,256
CHAIN_ADD=2**15
CHAIN_SET=2*GEN_MUL+17*CHAIN_ADD
YSND,YDLV,YDLVE=24,32,23
YRING_TOT=YSND+3+(M_G-2)*(YDLV+WY)+YDLVE+WY
def y_group(q):
 q_gen=cdiv(q,M_G)
 gen_y=YROUNDS*(cdiv(q_gen,YB)*YR_FIX+q_gen*YR_BLK)+q_gen*CHAIN_SET
 ring_y=cdiv(q*L_CHAIN*YRING_TOT,M_G)
 return gen_y+ring_y
XT_SLOTS=6
XT_W,XT_WHI=43,42
XTBLD_PAIRED=6*(4+2*11)+20*2*(4*2**XT_W+2*2**XT_WHI)
XTBLD_ALL=20*XTBLD_PAIRED
XSETUP=2**17
ENTRY_X=18
COIN_SCALAR=64
BATCH_SETUP=1352
INV_CALL=16
NODE_INIT=2**30
HEADER=80
BT_BASE,BT_W=0xfff8<<48,43
TBLD=4+7*11+20*(6*2**BT_W+2**(256-5*BT_W))
MSG=128
CAND_CAP=2**29
MULX,SQRX,ADDX=105+4,80+3,9+2*9
DBL,MADD_=2*MULX+5*SQRX+10*ADDX,8*MULX+3*SQRX+10*ADDX
CAND_LEAVES=6*CHAIN_ADD+12*MULX
CAND_ITEM=4*8*2*MSG+2*CAND_LEAVES+36*600+2*(384*200+256*(DBL+MADD_)+INV)
RESOLVE=2*CAND_ITEM+1
REQ_CAP=6*2**41
REQ_SERVE=2**10+GEN_MUL
RELAY_CAP=2**49
BARRIER=2**12*ALLN
OUTPUT=2**20
DELTA_F=2**60
QP,QQ=Q2//2,Q2//4
FIX_BATCH=sum(FIX_SITES_QUAD.values())*QQ
P1=Fr(6*(r-1)+1,r*r)
TWO_BOUND=Fr(2,r-1)
SHX_BAD=2**64*2**33
FOLD_BAD=2**66*C0
SQ_BAD=C0*C0+2*C0
FS_P1_BAD=C0
BDAY63_BAD=2**63*C0
BDAY64_BAD=2**64*C0
SITES=[(3*QP,SHX_BAD,1),(QP,FOLD_BAD,1),(3*QP,SHX_BAD,1),(QP,FOLD_BAD,12),(QP,FOLD_BAD,1),(3*QP,SHX_BAD,6),(QP,FOLD_BAD,6),(Q2,SQ_BAD,12),(Q2,FOLD_BAD,12),(2*Q2,FS_P1_BAD,6),(QP,BDAY63_BAD,12),(QP,FS_P1_BAD,6),(QP,BDAY64_BAD,2),(QP,FS_P1_BAD,8),(9*Q2,SHX_BAD,8)]
SITES_QUAD=dict(SHX=54,FDX=12,F5=4,FS_FA=12,BDAY63=2,BDAY64=2)
D_SITES=sum(Fr(n*bad*mult,r-1)for n,bad,mult in SITES)
D_YGEN=Fr(YROUNDS*(YPC+2*16+2),r-1)
D_BATCH=D_SITES+Q2*Fr(6,r-1)+D_YGEN
FIX_ROW=C0*2**192
FIX_FOLD=(2**66+1)*2**128*C0
FIX_SQ=2**193
FIX_LK6=15*C0*2**192
FIX_LK7=21*C0*2**192
FIX_BLK=30*C0*2**192
FIX_SUB=C0*2**192
FIX_BDAY63=((1<<128)+1)*(1<<63)*C0
FIX_BDAY64=((1<<128)+1)*(1<<64)*C0
FIXUPS=[("BLK",QP,FIX_LK7,12),("SHX",3*QP,FIX_ROW,1),("FDX",QP,FIX_FOLD,1),("SHX",3*QP,FIX_ROW,1),("FDX",QP,FIX_FOLD,12),("FDX",QP,FIX_FOLD,1),("FS",QP,FIX_SUB,6),("SHX",3*QP,FIX_ROW,6),("FDX",QP,FIX_FOLD,6),("F5",QP,FIX_BDAY63,12),("FADDW",QP,FIX_SUB,6),("F5",Q2,FIX_SQ,12),("BLK",Q2,FIX_LK7,12),("BLK",Q2,FIX_LK6,12),("FDX",Q2,FIX_FOLD,12),("FS",Q2,FIX_SUB,6),("FADDW",Q2,FIX_SUB,6),("BLK",2*Q2,FIX_BLK,6),("BLK",2*QQ,FIX_LK7,8),("F5",2*QQ,FIX_BDAY64,2),("FADDW",2*QQ,FIX_SUB,8),("BLK",2*Q2,FIX_LK6,8),("SHX",9*Q2,FIX_ROW,8)]
E_FIX=sum(Fr(n*bad*mult,r-1)for _,n,bad,mult in FIXUPS)+D_BATCH*FIX_BATCH
KB=functools.lru_cache(None)(lambda fixed:Fr(6*2**(256-fixed)*(KBF+1),(r-1)*KBF))
RHO=Fr(r,r-1)**2
CLUST_BATCH=M_G*LEAVES*Q2
FIXED_LEAVES=2*CLUST_BATCH
PARTNER_PAIRS=CLUST_BATCH*(CLUST_BATCH-1)
FA,FB=int("114ca50f7a8e2f3f657c1108d9d44cfd9",16),int("e4437ed6010e88286f547fa90abfe4c3",16)
MH_A,MH_B,ML=20,21,23
def isum(*rs):return sum(x[0]for x in rs),sum(x[1]for x in rs)
def _table(raw):
 cor=isum(raw,(-1,0))
 u3,uc3=(raw,raw,(-2-2*raw[1],-2*raw[0])),(cor,cor,(-2-2*cor[1],-2*cor[0]))
 wr=[[isum(uc3[m],uc3[(m+t)%3],u3[(m+2*t)%3],(0,2))for t in range(3)]for m in(0,1)]
 nr=[(-wr[0][t][1]-wr[1][t][1]-1,-wr[0][t][0]-wr[1][t][0]+1)for t in range(3)]
 return[(min(x[0]for x in rs),max(x[1]for x in rs))for rs in(wr[0],wr[1],nr)]
TFP_H,TFP_L=((0,1),(0,7)),((0,1),(0,7),(0,1))
ERRH,ERRL=_table(isum(*TFP_H)),_table(isum(*TFP_L))
WDH,WDL=sum(hi-lo for lo,hi in ERRH),sum(hi-lo for lo,hi in ERRL)
WH,WL=max(hi-lo for lo,hi in ERRH)+1,max(hi-lo for lo,hi in ERRL)+1
(ELO,EHI)=ERRH[0]
NWIN=2*(EHI-ELO)+(1-2*ELO)+(3+2*EHI)+2*(EHI+1-ELO)+2*(3*EHI+3-3*ELO)
D_WJ=TWO_BOUND*NWIN*2**192
D_LAB=TWO_BOUND*(WDH*2**(256-MH_A)+WDL*2**(256-ML))
D_LAB_SYM=TWO_BOUND*3*((WH-1)*2**(256-MH_A)+(WL-1)*2**(256-ML))
B_HA,B_HB,B_L=2**MH_A+WH-1,2**MH_B+WH-1,2**ML+WL-1
check("c9",FA*FA-FA*FB+FB*FB==p and(FB*BETA_P-(FA-FB))%p==0 and((FA-FB)*BETA_P+FA)%p==0 and B_HA*B_HB*(B_L*2**128+2*FA)*(B_L*2**128+4*FA-2*FB)<p*2**(MH_A+MH_B+2*ML)*(1+Fr(1,KBF))and(2*B_HA+1)*(2*B_HB+1)*(2*B_L*2**128+FA)*(2*B_L*2**128+2*FA-FB)<16*p*2**(MH_A+MH_B+2*ML)*(1+Fr(1,KBF))and(B_L*2**128+2*FA)*(B_L*2**128+4*FA-2*FB)<p*2**(2*ML)*(1+Fr(1,KBF))and B_L*2**233<p*(1+Fr(1,KBF))and(B_HA*2**192+2*FA)*(B_HB*2**192+4*FA-2*FB)<p*2**(128+MH_A+MH_B)*(1+Fr(1,KBF)))
check("c10",isum(*TFP_H)==(0,8)and isum(*TFP_L)==(0,9)and NWIN==564 and(1<<322)*C0<(1<<99)*p)
check("c11",ERRH==[(-20,26),(-20,26),(-53,40)]and ERRL==[(-22,29),(-22,29),(-59,44)]and(WDH,WDL,WH,WL)==(185,205,94,104))
check("c12",D_WJ<P2(Fr(-5386,100))[0],"2^%.4f"%lg(D_WJ))
check("c13",sum(n for n,b,_ in SITES if b==SHX_BAD)==SITES_QUAD["SHX"]*QQ and sum(n for n,b,_ in SITES if b==FOLD_BAD)==SITES_QUAD["FDX"]*QQ and sum(n for n,b,_ in SITES if b==SQ_BAD)==SITES_QUAD["F5"]*QQ and sum(n for n,b,_ in SITES if b==FS_P1_BAD)==SITES_QUAD["FS_FA"]*QQ and sum(n for n,b,_ in SITES if b==BDAY63_BAD)==SITES_QUAD["BDAY63"]*QQ and sum(n for n,b,_ in SITES if b==BDAY64_BAD)==SITES_QUAD["BDAY64"]*QQ and len(SITES)==15)
FIX_RUNS={m:sum(n for ms,n,_,_ in FIXUPS if ms==m)for m in FIX_SITES_QUAD}
check("c14",FIX_RUNS=={m:FIX_SITES_QUAD[m]*QQ for m in FIX_SITES_QUAD}and sum(FIX_SITES_QUAD.values())==116)
check("c15",FIX_BATCH*C0<2**63 and FIX_BATCH*C0>>32==FIX_BATCH+(FIX_BATCH*977>>32),"2^%.3f"%lg(FIX_BATCH*C0))
check("c16",D_LAB<P2(Fr(-1088,100))[0]and 3*D_WJ+D_LAB_SYM<P2(Fr(-1068,100))[0],"2^%.3f 2^%.4f"%(lg(D_LAB),lg(3*D_WJ+D_LAB_SYM)))
check("c17",D_BATCH<P2(Fr(-1268,10))[0],"2^%.3f"%lg(D_BATCH))
G_GROUPS=2048
G_EXP=25
def group_sizes(q1):
 ms=[0]
 g_pow=G_GROUPS**G_EXP
 q32=q1**32
 for g in range(1,G_GROUPS):
  rhs=q32*(g**G_EXP)
  m=int(q1*((g/G_GROUPS)**(G_EXP/32.0)))+2
  while(m**32)*g_pow>rhs:m-=1
  while((m+1)**32)*g_pow<=rhs:m+=1
  ms.append(m)
 ms.append(q1)
 qs=[ms[_g+1]-ms[_g]for _g in range(G_GROUPS)]
 return qs,ms
HCHG_NEW,HCHG_VISIT,HCHG_END,HCHG_VER=4,4,0,5
HCHG_LS0,HCHG_LSA=2,1
HCHG_LS1=5+HCHG_LSA
check("c18",HCHG_NEW==4 and HCHG_VISIT==4 and HCHG_END==0 and HCHG_VER==5 and HCHG_LS1==6)
DLVHB=256
def dlvh_seg(n,m):
 if m==0:return 1
 head=n%2;rest=m-head
 return 5+head*5+(rest//2)*6+(rest%2)*4-(rest==0)
def dlvh_copy(m1,m2,K_=K2B):
 if m1==0 and m2==K_:return DLVHB
 return 23-(m1>0)+dlvh_seg(0,m1)+dlvh_seg(m1,m2)+dlvh_seg(m1+m2,K_-m1-m2)
DLVHS=max(dlvh_copy(a,b)for a in range(K2B+1)for b in range(K2B+1-a)if(a,b)!=(0,K2B))
BLOCKH=FLUSH+DONE_C+DLVHB
FINALH=FINH_DEST+FINH_SENT+DONE_C+DLVHS+MSG_LINK+(DLVHS-DLVHB)
check("c19",DLVHB==256 and DLVHS==427 and BLOCKH==FLUSH+DONE_C+DLVHB and FINALH==641,"DLVHS %d FINALH %d"%(DLVHS,FINALH))
ZW=16384
ZERO_PASS=2*(ZW//256)+2
DLVH0=15
STOP_CHECK=5
NBAR_COUNT=6+19*G_GROUPS+1
NBAR_ALLOW=40000
PHASE_CTRL_G=2**20
KB_BUCKET=KB(D1+D2+HB)
KB64=KB(64)
R2_UP_TIB=682
R2_TIB=38220
CAP_H=((R2_TIB*(2**40))//(32*K2B))*K2B
P1_GOOD=Fr(6*(r-1),r*r)*(1-3*D_WJ-D_LAB_SYM)
K_U_CLUST=LEAVES*M_G*MX
SAMEX_FAC=Fr(CLUST_BATCH-K_U_CLUST,CLUST_BATCH-1)*Fr(r,r-2)+Fr(96*(K_U_CLUST-1),CLUST_BATCH-1)
EF_Q=E_FIX*C0/(2**32*(DELTA_F-FIX_BATCH))
KB_M=KB(KEYBITS+D1+D2)
KB_D1=KB(D1)
SIP=240*RHO*KB_M
LOG2_POOLS=15
N_POOLS=1<<LOG2_POOLS
R_POOL=N_POOLS//M_G
K_GEN=NODES//N_POOLS
def gf1024_mul(a,b):
 res=0
 for _i in range(10):
  if b>>_i&1:res^=a<<_i
 for _i in range(18,9,-1):
  if res>>_i&1:res^=0x409<<(_i-10)
 return res
def pool_of_gen(g):
 c,m=divmod(g,M_G)
 return R_POOL*m+(((c&31)^gf1024_mul((c>>5)&0x3FF,m))&31)
def _vand_gf1024_ok():
 b=2
 for _ in range(10):b=gf1024_mul(b,b)
 if b!=2:return False
 for _d in(2,5):
  bd=2
  for _ in range(_d):bd=gf1024_mul(bd,bd)
  if bd==2 or any(gf1024_mul(bd^2,z)==0 for z in range(1,1024)):return False
 ns_tot=0
 for d in range(1024):
  ws=[gf1024_mul(1<<j,d)&31 for j in range(10)]
  vals=[0]
  for w in ws:vals+=[x^w for x in vals]
  if d and len(set(vals))!=32:return False
  cnt=[0]*32
  for v in vals:cnt[v]+=1
  for du in range(32):
   ks=cnt[du]
   ns=1024*1024+2*ks*max(0,ks-1)
   ns_tot+=ns*((1<<41)-(d==0 and du==0))
 return ns_tot*2304<9216*1024*1024*((1<<56)-1)
def _bx_tuples_ok():
 if(LAM_R*LAM_R+LAM_R+1)%r!=0:return False
 mu6=((1,0),(-1,0),(0,1),(0,-1),(-1,-1),(1,1))
 norms,n2,counts=set(),set(),{}
 for u0 in mu6:
  n2.add(1)
  for u1 in mu6:
   a,b=u0[0]+u1[0],u0[1]+u1[1]
   if a or b:n2.add(a*a-a*b+b*b)
 _inv={}
 for k in(1,2,3,4):
  seen=set()
  for asg in itertools.product(range(k),repeat=4):
   if len(set(asg))!=k:continue
   for u0 in mu6:
    for u1 in mu6:
     for u2 in mu6:
      sums=[[0,0]for _ in range(k)]
      for pos,u in enumerate((u0,u1,u2,(1,0))):sums[asg[pos]][0]+=u[0];sums[asg[pos]][1]+=u[1]
      if any(a==0 and b==0 for a,b in sums):continue
      for a,b in sums:
       norms.add(a*a-a*b+b*b)
       if not(a+b*LAM_R)%r:return False
      den=(sums[-1][0]+sums[-1][1]*LAM_R)%r
      inv_l=_inv.get(den)or _inv.setdefault(den,pow(den,-1,r))
      seen.add(tuple(((a+b*LAM_R)*inv_l)%r for a,b in sums[:-1]))
  counts[k]=len(seen)
 return norms=={1,3,4,7,9,12,13,16}and n2=={1,3,4}and counts=={1:1,2:78,3:252,4:216}and 96*(78+252+216)<65536
SAME_M_NODE_PAIRS=M_G*(N_CLUST*N_CLUST*(R_POOL-1)//R_POOL)//2
SHARED_BOTH_CLUST_PAIRS=(N_CLUST*N_CLUST//(R_POOL*R_POOL)-N_CLUST)//2
ALL_CLUST_PAIRS=N_CLUST*(N_CLUST-1)//2
GOOD_NODE_PAIRS=SAME_M_NODE_PAIRS+SHARED_BOTH_CLUST_PAIRS*(M_G*(M_G-1)//2)+(ALL_CLUST_PAIRS-SHARED_BOTH_CLUST_PAIRS)*(M_G*(M_G-1))
YR_MINE=992
PTS_MINE_CHAIN=YR_MINE*YPC
NBATCH_MINE=62835969*YB
PTS_PER_GEN_MINE=NBATCH_MINE*PTS_MINE_CHAIN
M_POOL=K_GEN*PTS_PER_GEN_MINE
S_MINE_TOT=NODES*PTS_PER_GEN_MINE
N_CHAINS_POOL=K_GEN*NBATCH_MINE
PAIRS_PER_MINE_CHAIN=PTS_MINE_CHAIN*(PTS_MINE_CHAIN-1)//2
P_BUCKET_MAX=Fr(6*(1<<63),r-1)
P_SAME_CHAIN_POOL=N_CHAINS_POOL*PAIRS_PER_MINE_CHAIN*P_BUCKET_MAX
P_CHAIN_COUPLE_POOL=N_CHAINS_POOL*PAIRS_PER_MINE_CHAIN*(M_POOL*M_POOL+(N_CHAINS_POOL-1)*PAIRS_PER_MINE_CHAIN)*(P_BUCKET_MAX**2)
GAMMA_MINE_CROSS=M_POOL*(M_POOL-1)//2-N_CHAINS_POOL*PAIRS_PER_MINE_CHAIN
P_BDAY_HIT_MIN=Fr((r-1)*(r-1),(1<<193)*r*r)-Fr(6*(r-1),r*r)
MU_POOL_LO=GAMMA_MINE_CROSS*P_BDAY_HIT_MIN
CS_POOL=GAMMA_MINE_CROSS*(2*PTS_MINE_CHAIN*M_POOL)*(2*P_BUCKET_MAX*P_BUCKET_MAX)
_P_A102=Fr(6*(1<<102),r-1)
_P_A174=Fr(6*(1<<(256-D1-D2)),r-1)
_P_A215=Fr(6*(1<<(256-D1)),r-1)
P_CELL154_MASK_POOL=(Fr(M_POOL*(M_POOL-1),2)*M_POOL+2*N_CHAINS_POOL*PAIRS_PER_MINE_CHAIN*M_POOL)*P_BUCKET_MAX*_P_A102
POIS_LE3=exp_neg_hi(MU_POOL_LO)*sum((MU_POOL_LO**_k)/math.factorial(_k)for _k in range(4))
P_FAIL_POOL=POIS_LE3+CS_POOL+P_CELL154_MASK_POOL+P_SAME_CHAIN_POOL+P_CHAIN_COUPLE_POOL
P_FAIL_ALL_POOLS=N_POOLS*P_FAIL_POOL
LAM_ORD=N_POOLS*(M_POOL*M_POOL)*P_BUCKET_MAX
POINTWISE_X_MAX=M_POOL*M_POOL*P_BUCKET_MAX
SAME_OFF_VAND=(2*LEAVES*MX)**2+(M_G//R_POOL)*(LEAVES*MX)**2+Fr(K_U_CLUST*K_U_CLUST,N_CLUST)
P_BAD_X=65536*sum((LAM_ORD**(_d-1))*(LAM_ORD/r)for _d in(2,3,4))+65536*sum((LAM_ORD**(_d-1))*N_POOLS*(M_POOL*M_POOL)*Fr(30*(1<<64),r*r)for _d in(2,3))
P_X_ORBIT_ALL_POOLS=LAM_ORD*(6*LAM_ORD/r)
P_MINE_ZERO_ALL=Fr(S_MINE_TOT,r)
P_SAME_A_ALL=N_POOLS*Fr(M_POOL*(M_POOL-1)//2*6,r-1)
D_XMINE_ALL=NODES*NBATCH_MINE*Fr(YR_MINE*(YPC+34),r-1)
_MINE_Y_MEAN=Fr(S_MINE_TOT,PTS_MINE_CHAIN)*Fr(6*(1<<(256-D1-D2)),r-1)
_MINE_K_CH=CAP_H//PTS_MINE_CHAIN
P_MINE_HOME_CAP=Fr(1,1<<200)
E_CAND_MINE_EX=N_POOLS*Fr(M_POOL*(M_POOL-1),2)*_P_A102
E_CAND_MINE=cdiv(E_CAND_MINE_EX.numerator,E_CAND_MINE_EX.denominator)
E_VIS_MINE=cdiv((S_MINE_TOT*(S_MINE_TOT-1)//2)*6*(1<<(256-131)),r-1)
E_K2_MINE=cdiv(N_POOLS*(M_POOL*(M_POOL-1)//2)*6*(1<<(256-131)),r-1)
_OTH_H_OUT=N_POOLS*Fr(M_POOL*(M_POOL-1),2)*_P_A174*_P_A102+2*PTS_MINE_CHAIN*M_POOL*_P_A102
_OTH_G_OUT=2*PTS_PER_GEN_MINE*M_POOL*_P_A102
_M_R_OUT=18*N_POOLS*Fr(M_POOL*(M_POOL-1),2)*_P_A215*_P_A102
_RL_CR_CAP_MINE=1<<20
_RL_CR_MINE=36*PTS_MINE_CHAIN*M_POOL*_P_A102
_VAR_R_OUT=324*N_POOLS*M_POOL*_P_A215*(1+RHO*RHO*PTS_MINE_CHAIN*_P_A215)*(M_POOL*_P_A102)**2+18*PTS_MINE_CHAIN*_M_R_OUT
P_MINE_CAPS=P_MINE_HOME_CAP+LAM_ORD*(_OTH_H_OUT/Fr(CAND_CAP-2*PTS_MINE_CHAIN)+2*_OTH_G_OUT/Fr(REQ_CAP-4*PTS_MINE_CHAIN)+_RL_CR_MINE/_RL_CR_CAP_MINE+_VAR_R_OUT/Fr(RELAY_CAP-36*PTS_MINE_CHAIN-_RL_CR_CAP_MINE-_M_R_OUT)**2)
check("c19b",_vand_gf1024_ok()and _bx_tuples_ok()and all(pool_of_gen(1024*(_p&31)+(_p>>5))==_p for _p in range(N_POOLS))and NBATCH_MINE%YB==0 and _MINE_K_CH-2*_MINE_Y_MEAN-82>200 and MU_POOL_LO>56 and P_FAIL_ALL_POOLS<P2(Fr(-28,1))[0]and P_BAD_X<P2(Fr(-140,1))[0]and P_MINE_CAPS<P2(Fr(-19,1))[0],"mu_p %.2f p_fail 2^%.2f p_caps 2^%.2f"%(float(MU_POOL_LO),lg(P_FAIL_ALL_POOLS),lg(P_MINE_CAPS)))
EXCH_ROW=ALLN*4*2**10
HR_INIT=N_HOME*(INIT3_BKT+1+DLVH0)
VISIT=HCHG_VISIT+Fr(HCHG_END,2)
GEN_MINE_Y=YR_MINE*((NBATCH_MINE//YB)*YR_FIX+NBATCH_MINE*YR_BLK)+NBATCH_MINE*(CHAIN_SET+2*COIN_SCALAR+64)
GEN_MINE_NODE=64+GEN_MINE_Y+PTS_PER_GEN_MINE*XEMIT_HOT
_DEF_MINE_EX=D_XMINE_ALL*(PTS_MINE_CHAIN*CAP_H*VISIT+PTS_MINE_CHAIN*S_MINE_TOT*(VISIT+HCHG_VER+HCHG_LS1+RESOLVE+2*REQ_SERVE+8*MSG))
DEF_MINE=cdiv(_DEF_MINE_EX.numerator,_DEF_MINE_EX.denominator)
PSEL_POOL=8+1024*256
ROUTE_MINE=(S_MINE_TOT*(FLUSH+DONE_C1)//K1+S_MINE_TOT*H2S_ARR+cdiv(S_MINE_TOT,K1)*H2S_FULL+NODES*(FINAL*F1+FIN_FIX)+cdiv(S_MINE_TOT*BLOCKH,K2B)+N_REL*(FINALH*F2+FIN_FIX)+N_HOME*EXPECTH+N_REL*(EXPECT+RDONE)+S_MINE_TOT*(HCHG_NEW+RCV)+VISIT*E_VIS_MINE+HCHG_VER*E_K2_MINE+HCHG_LS1*E_CAND_MINE+E_CAND_MINE*(RESOLVE+2*REQ_SERVE+8*MSG)+DEF_MINE+EXCH_ROW+N_POOLS*(PSEL_POOL+MSG+4)+ALLN*(MSG+STOP_CHECK)+NODES*MX*(4*MSG)+3*N_HOME+6*S_MINE_TOT+HR_INIT+NODES*(INIT3_BKT*F1+2)+N_REL*((INIT3_BKT+TINIT_ENTRY)*F2+64)+5*BARRIER)
E_MINE=NODES*GEN_MINE_NODE+ROUTE_MINE
_GLOB_BASE=P_FAIL_ALL_POOLS+P_SAME_A_ALL+P_BAD_X+D_XMINE_ALL+P_MINE_ZERO_ALL+P_X_ORBIT_ALL_POOLS+P_MINE_CAPS
_CS_K=P1*P1+Fr(36,(r-1)**2)
_KBD12=KB(D1+D2)
_COV_K=Fr(K_U_CLUST*L_CHAIN,r-1)*(_KBD12**2)
_ACROSS0=SAMEX_FAC*PARTNER_PAIRS+CLUST_BATCH*CLUST_BATCH+96*K_U_CLUST*K_U_CLUST
_FDIFF0=SAME_OFF_VAND+160*K_U_CLUST*K_U_CLUST*((POINTWISE_X_MAX/MX)**4)*KB_M+2*CLUST_BATCH*CLUST_BATCH*RHO*KB_M
_ROW_K1=1+Fr(1,r)
def success_prefix(q1p):
 n0p=LEAVES*q1p*Q2
 Sp=NODES*n0p
 gamma=GOOD_NODE_PAIRS*(n0p*n0p)
 mu_raw=gamma*P1
 mu=gamma*P1_GOOD
 cs=gamma*(2*CLUST_BATCH*Sp)*_CS_K
 e_hi=exp_neg_hi(mu)
 d_f=q1p*(EF_Q+D_BATCH)
 m_p2=Sp*_KBD12+1
 v_p2=Sp*(K_U_CLUST*_KBD12+_COV_K)+1
 d_p2=v_p2/Fr(CAP_H-FIXED_LEAVES-m_p2)**2
 lam_int=RHO*(Sp+FIXED_LEAVES)*KB_M
 alpha_cross=(2*K_U_CLUST+SAMEX_FAC*(CLUST_BATCH-1))*Sp+_ACROSS0
 other_h=RHO*(Fr(Sp)*Sp*_KBD12*KB_M+alpha_cross*KB_M+2*L_CHAIN*_KBD12)+SIP*Sp
 d_cl=other_h/(CAND_CAP-76)
 other_g=RHO*KB_M*(Fr(n0p)*Sp*2+2*LEAVES*MX*Sp+96*K_U_CLUST*n0p+96*LEAVES*MX*K_U_CLUST)
 d_rq=other_g/(REQ_CAP-2*L_CHAIN-600)
 rl_cap_cr=RELAY_CAP//4
 rl_cr=18*FIXED_LEAVES*Fr(Sp)*KB_M
 rl_ediff=RHO*Fr(Sp)*Sp*KB_D1*KB_M
 rl_fdiff=_FDIFF0+2*CLUST_BATCH*Fr(Sp)*KB_M
 v_rl=648*(rl_ediff*rl_fdiff+CLUST_BATCH*Fr(Sp)*RL_SAME_PL)
 d_rl=rl_cr/rl_cap_cr+v_rl/_RL_DEN2
 per_leaf=D_BATCH+d_f+d_p2+d_cl+d_rq+d_rl
 rows_=mu_raw*_ROW_K1*(Fr(1,r)+lam_int+2*per_leaf)
 g_def=NODES*q1p*D_BATCH
 glob=_GLOB_BASE+NODES*(3*q1p+2*NBATCH_MINE+1)*Fr(1,2**254)+g_def
 lo=1-e_hi-cs-rows_-glob
 return dict(mu=mu,lo=lo,cs=cs,rows=rows_,glob=glob,Sp=Sp,n0p=n0p,d_p2=d_p2,m_p2=m_p2,d_f=d_f,lam_int=lam_int,d_cl=d_cl,d_rq=d_rq,d_rl=d_rl,g_def=g_def)
def grouped_full_lo(q):return success_prefix(q)["lo"]
Q1GD=70000000000
Q1G,MQ1=group_sizes(Q1GD)
RL_SAME_PL=RHO*Fr(CLUST_BATCH-1,2)*KB_D1*KB_M+96*K_U_CLUST*KB_M+SIP
_SG_G=[NODES*LEAVES*q*Q2 for q in Q1G]
SG_MAX=max(_SG_G)
SPROD_MAX=max(_SG_G[g]*(NODES*LEAVES*MQ1[g+1]*Q2)for g in range(G_GROUPS))
M_VG_MAX=18*(RHO*Fr(SPROD_MAX)*KB_D1*KB_M+Fr(SG_MAX)*RL_SAME_PL)
_RL_DEN2=(Fr(RELAY_CAP-18*FIXED_LEAVES-RELAY_CAP//4)-M_VG_MAX)**2
N0G=LEAVES*Q1GD*Q2
SG_FULL=NODES*N0G
LO_RUN=grouped_full_lo(Q1GD)
GATE_BITS,GATE_INSNS=128,256
THRESH_GATE=Fr(3905,10**4)
K_GATE=cdiv(THRESH_GATE.numerator*LO_RUN.denominator*(1<<GATE_BITS),THRESH_GATE.denominator*LO_RUN.numerator)
P_GATE=Fr(K_GATE,1<<GATE_BITS)
GATED_LO=P_GATE*LO_RUN
check("c20 %d"%Q1GD,LO_RUN>=Fr(9985,10000),"lo_run %.12f"%(math.floor(LO_RUN*10**12)/10**12))
check("c21",0<K_GATE<(1<<GATE_BITS)and GATED_LO>=THRESH_GATE>THRESH and Fr(K_GATE-1,1<<GATE_BITS)*LO_RUN<THRESH_GATE and 2*COIN_SCALAR+8+3+4<=GATE_INSNS,"p_gate %.9f gated_lo %.12f"%(float(P_GATE),math.floor(GATED_LO*10**12)/10**12))
check("c22",all((MQ1[g]**32)*(G_GROUPS**G_EXP)<=(Q1GD**32)*(g**G_EXP)<((MQ1[g]+1)**32)*(G_GROUPS**G_EXP)for g in range(1,G_GROUPS)))
TPRE=[success_prefix(m)for m in MQ1[:-1]]
def ceil_256(x):
 x=min(Fr(1),max(Fr(0),x))
 return Fr(cdiv(x.numerator<<256,x.denominator),1<<256)
PHI=[Fr(1)]+[ceil_256(1-t["lo"])for t in TPRE[1:]]
TFULL=success_prefix(Q1GD)
check("c23 %.6f"%lg(SG_FULL),SG_FULL==NODES*6*Q1GD*Q2==TFULL["Sp"])
check("c24 %.9f"%float(TFULL["mu"]),TFULL["mu"]==GOOD_NODE_PAIRS*(N0G*N0G)*P1_GOOD and GATED_LO>=THRESH)
check("c25",TFULL["cs"]<Fr(1,2**82),"2^%.4f"%lg(TFULL["cs"]))
check("c26",TFULL["lam_int"]<P2(Fr(-377,10))[0],"2^%.3f"%lg(TFULL["lam_int"]))
check("c27 %.3f"%lg(Q1GD*E_FIX),Q1GD*E_FIX<P2(Fr(39,1))[0]and DELTA_F>2**20*Q1GD*E_FIX+FIX_BATCH)
check("c28",TFULL["d_f"]<P2(Fr(-214,10))[0],"2^%.3f"%lg(TFULL["d_f"]))
check("c29",TFULL["d_p2"]<P2(Fr(-336,10))[0],"2^%.3f"%lg(TFULL["d_p2"]))
check("c30",TFULL["d_cl"]<P2(Fr(-175,10))[0]and TFULL["d_rq"]<P2(Fr(-175,10))[0]and TFULL["d_rl"]<P2(Fr(-185,10))[0],"2^%.3f 2^%.3f 2^%.3f"%(lg(TFULL["d_cl"]),lg(TFULL["d_rq"]),lg(TFULL["d_rl"])))
check("c31",TFULL["g_def"]<P2(Fr(-255,10))[0]and TFULL["glob"]<P2(Fr(-19,1))[0],"2^%.3f"%lg(TFULL["glob"]))
check("c32",PHI[0]==1 and all(0<=x<=1 for x in PHI)and all(PHI[g]>=PHI[g+1]for g in range(G_GROUPS-1))and all(PHI[g]>=min(Fr(1),max(Fr(0),1-TPRE[g]["lo"]))for g in range(1,G_GROUPS)),"PHI[1] %.4f PHI[16] %.4f PHI[31] %.4f"%(float(PHI[1]),float(PHI[16]),float(PHI[31])))
HOME_MAX=CAP_H
check("c33",HOME_MAX==CAP_H<2**59,"2^%.4f"%lg(HOME_MAX))
PAIRM=KB_BUCKET+Fr(1,r)
PAIRK2=KB(D1+D2+64)+Fr(1,r)
PJ=Fr(2,3)+2*D_WJ+Fr(1,r-1)
check("c34",VISIT==4)
SG=[NODES*LEAVES*q*Q2 for q in Q1G]
WG,SAMEX,WGJ,WGK2,WGL,WGR,WGG=[],[],[],[],[],[],[]
for _g,_q in enumerate(Q1G):
 _s=SG[_g]
 _leaf=NODES*_q*(fwd_count(Q2)+INV+INV_CALL+BATCH_SETUP)+NODES*y_group(_q)+int(NODES*_q*bwd_count(Q2))+int(NODES*_q*((KEY_INC-EMIT_INC)*Q2+KEY_FIX))
 _route=_s*(EMIT*K+FLUSH+DONE_C1)//K+_s*H2S_ARR+(_s//K)*H2S_FULL+NODES*(FINAL*F1+FIN_FIX)+_s*BLOCKH//K2B+N_REL*(FINALH*F2+FIN_FIX)+N_HOME*EXPECTH+N_REL*(EXPECT+RDONE)+NODES*(INIT3_BKT*F1+1)
 _gp=NODES*LEAVES*MQ1[_g]*Q2
 _samex=N_CLUST*_q*CLUST_BATCH*(CLUST_BATCH-1)//2
 _main=_s*_gp+_s*(_s-1)//2-_samex
 _evis=_main*PAIRM
 _gfix=_s*HCHG_NEW+ALLN*STOP_CHECK
 WGK2.append(HCHG_VER*_main*PAIRK2)
 _group=_gfix+VISIT*_evis+WGK2[-1]
 _jkey=JPATH*_s*PJ+_s//3*SJOIN
 WGJ.append(_jkey);WGL.append(_leaf);WGR.append(_route);WGG.append(_gfix+VISIT*_evis)
 WG.append(_leaf+_route+_group+_jkey+EXCH_ROW)
 SAMEX.append(_samex)
check("c35",sum(Q1G)==Q1GD and min(Q1G)>0 and all(s%K==0 for s in SG))
check("c36",NODES*2**50<=2**116 and NODES<2**67 and max(Q1G)*LEAVES*Q2//K+1<2**50)
ZERO_INIT=N_HOME*(2**HB//ZW*ZERO_PASS+1)
check("c37",ZERO_PASS==130 and 2**HB%ZW==0 and(2**64-2**56)%2048==0)
E_PRE=ALLN*NODE_INIT+E_MINE+NODES*(3*cdiv(Q1GD,M_G)+1)*COIN_SCALAR+NODES*MX*XSETUP+NODES*Q2*ENTRY_X+NODES*HEADER+N_REL*(INIT3_BKT*F2+1+TINIT_ENTRY*F2+1)+NODES*(INIT3_BKT*F1+2)+ZERO_INIT+HR_INIT+NODES*XTBLD_ALL+NODES*TBLD
_SIPT=240*RHO
P_COIN_FAIL_ALL=NODES*(3*Q1GD+2*NBATCH_MINE+1)*Fr(1,2**254)
P_EXC=P_X_ORBIT_ALL_POOLS+P_BAD_X+D_XMINE_ALL+P_MINE_ZERO_ALL+P_COIN_FAIL_ALL
W_DET_MAX=(NODES*Q1GD*FIX_BATCH*max(FIX_EXTRA.values())+SG_FULL*JPATH+SG_FULL*CAP_H*(VISIT+HCHG_VER+HCHG_LS1)+N_HOME*CAND_CAP*(RESOLVE+2*REQ_SERVE))
_EXC_W_EX=P_EXC*W_DET_MAX
EXC_WORK=cdiv(_EXC_W_EX.numerator,_EXC_W_EX.denominator)
_ACROSS_FULL=2*K_U_CLUST*SG_FULL+SAMEX_FAC*SG_FULL*(CLUST_BATCH-1)+SAMEX_FAC*PARTNER_PAIRS+CLUST_BATCH*CLUST_BATCH+96*K_U_CLUST*K_U_CLUST
ELIST=RHO*(Fr(SG_FULL)*SG_FULL*KB_M+_ACROSS_FULL*KB_M)+_SIPT*KB_M*SG_FULL+72
EK2=RHO*(_ACROSS_FULL*KB64)+_SIPT*KB64*SG_FULL+72
SXK2=SAMEX_FAC*PAIRM*sum(SAMEX)
SXV=SAMEX_FAC*VISIT*PAIRM*sum(SAMEX)+EXC_WORK
BRACKET=RHO*(Fr(SG_FULL)*SG_FULL*KB(D1+D2)+_ACROSS_FULL)+_SIPT*SG_FULL+72
DEFV=VISIT*D_BATCH*(HOME_MAX*SG_FULL+sum(SAMEX))
DEFL_PAIRS=2*D_BATCH*BRACKET+D_BATCH*HOME_MAX*SG_FULL
DEFK2=HCHG_VER*DEFL_PAIRS
DEFL=HCHG_LS1*DEFL_PAIRS
DEFJ=JPATH*D_BATCH*SG_FULL
E_FIX_WORK=NODES*Q1GD*E_FIX*max(FIX_EXTRA.values())
E_CAND_WORK=(ELIST+DEFL_PAIRS)*RESOLVE
E_REQ_WORK=(ELIST+DEFL_PAIRS)*(2*REQ_SERVE)
E_SHARED_REST=E_PRE+NODES*GEN_MUL+(2+2*G_GROUPS)*ALLN*MSG+NBAR_ALLOW*BARRIER+ALLN*PHASE_CTRL_G+OUTPUT
E_RUN=E_SHARED_REST+E_FIX_WORK+E_CAND_WORK+E_REQ_WORK+sum(PHI[_g]*WG[_g]for _g in range(G_GROUPS))+SXV+EK2*HCHG_VER+ELIST*HCHG_LS1+DEFV+DEFK2+DEFL+DEFJ+SXK2*HCHG_VER
E_TOTAL=GATE_INSNS+P_GATE*E_RUN
def claim6_of(x):
 c=Fr(math.ceil(lg(x)*10**6),10**6)
 for _ in range(50):
  if P2(c)[0]>=x and P2(c-Fr(1,10**6))[1]<x:return c
  c+=Fr(1,10**6)if P2(c)[0]<x else-Fr(1,10**6)
 raise AssertionError("claim6 not certified")
E6=claim6_of(E_TOTAL)
W_E=E_TOTAL/C
def claim15_of(w):
 getcontext().prec=80
 start=(Decimal(w.numerator).ln()-Decimal(w.denominator).ln())/Decimal(2).ln()*10**15
 c=Fr(int(start.to_integral_value(rounding=ROUND_CEILING)),10**15)
 for _ in range(400):
  if P2(c)[0]>=w and P2(c-Fr(1,10**15))[1]<w:return c
  c+=Fr(1,10**15)if P2(c)[0]<w else-Fr(1,10**15)
 raise AssertionError("E-claim not certified")
CLAIM_E=claim15_of(W_E)
JROW_E=sum(PHI[_g]*WGJ[_g]for _g in range(G_GROUPS))+DEFJ
K2ROW_E=sum(PHI[_g]*WGK2[_g]for _g in range(G_GROUPS))+SXK2*HCHG_VER+EK2*HCHG_VER
PER_LEAF_OPEN=E_RUN/SG_FULL
PER_LEAF_E=PER_LEAF_OPEN
LOW,UP_BASE=2**56,2**64-2**56
def valid_va(a):return 0<=a<2**64 and(a>>56)in(0,255)
TIB=2**40
XT_PAIR6_BYTES=12*2560*TIB
AXYX_BASE=4096*TIB
AXYX_HI=20480*TIB
XT_BASE=24576*TIB
XT_END=XT_BASE+XT_PAIR6_BYTES
CN_BASE=84*TIB
CNM_BASE=90*TIB
H1L_BASE=57344*TIB
H1U_BASE=2**64-6144*TIB
H2U_BASE=2**64-16384*TIB
check("c38",BT_BASE==2**64-2**51 and AXYX_BASE%(2**49)==0 and AXYX_HI%(2**48)==0 and XT_BASE%(2**48)==0 and valid_va(XT_END-1)and all(valid_va(BT_BASE|((w<<8)%2**64))for w in(0,1,2**43-1,2**64-1)))
YR_SLOT=2**31
YR_BASE=32*TIB+YR_SLOT
EB_END=33*TIB+2**34
EB_START=EB_END-Q2*512
T34=34*TIB
LINK_GEN,LINK_GEN_MSG,PORTY_GEN=T34,T34+0x1000,T34+0xB000
G34_END=38*TIB
REQ_LO=96*TIB-64
REQ_END=REQ_LO+REQ_CAP*64+64
LINK_HR=2**56-2**40
LINK_HR_MSG=LINK_HR+0x8000
NEXT_LO=2**64-57344*TIB
NEXT_BASE=NEXT_LO+341*TIB
R2_UP_BYTES=R2_UP_TIB*TIB
R2_LOW_BYTES=CAP_H*32-R2_UP_BYTES
HCC_BYTES=2**20+8*2048+4096+CAND_CAP*528
G34_CONTENT=2**20+8*2048+MX*6*XT_SLOTS*64+2048+MX*64+(2*MX+2*cdiv(Q1GD,M_G)+Q1GD+1)*32+2**30+YB*2304+1024*128+512
REGISTRY=[("G","LC1",0,F1*8,"init3"),("G","YR",YR_BASE-64,YR_BASE-64+(YB+1)*YR_SLOT,"yround"),("G","EB",EB_START,EB_END,"setup:ENTRY_X fwd"),("G","LG",LINK_GEN,LINK_GEN+4096,"MMIO"),("G","LGM",LINK_GEN_MSG,LINK_GEN_MSG+4096,"MMIO"),("G","PY",PORTY_GEN,PORTY_GEN+512,"MMIO"),("G","G34",T34+0x10000,G34_END,"setup:NODE_INIT setup:XSETUP setup:COIN_SCALAR ydlv0 ydsp wcnt dlvq psel"),("G","CN",CN_BASE,CN_BASE+64*Q1GD,"setup:BATCH_SETUP"),("G","CNM",CNM_BASE,CNM_BASE+64*NBATCH_MINE,"setup:COIN_SCALAR"),("G","REQ",REQ_LO,REQ_END,"dlvq"),("G","AXY",AXYX_BASE,XT_BASE,"tbld"),("G","XT6",XT_BASE,XT_END,"tbld"),("G","H1L",H1L_BASE,H1L_BASE+4096*TIB,"init3"),("G","H1U",H1U_BASE,BT_BASE,"init3"),("G","BT",BT_BASE,2**64,"tbld"),("R","RCS",0x10000,0x10000+2**20+8*2048,"setup:NODE_INIT h2s0"),("R","LC2",2**45,2**45+F2*8,"init3"),("R","RC2",3*2**45,3*2**45+F2*8,"tinit"),("R","E2P",2**47,2**47+2*F2*8,"analytic"),("R","DPW",2**48-4096,2**48,"h2s0"),("R","RCB",2**52,2**52+0x10000+RELAY_CAP*64,"setup:NODE_INIT setup:EXCH_ROW dlvq"),("R","LR",LINK_HR,LINK_HR+0x8000,"MMIO"),("R","LRM",LINK_HR_MSG,LINK_HR_MSG+0x8000,"MMIO"),("R","RCT",2**56-8192,2**56-4096,"dlvq"),("R","H2U",H2U_BASE,2**64,"init3"),("H","R2L",0,R2_LOW_BYTES,"dlvh"),("H","LH",LINK_HR,LINK_HR+0x8000,"MMIO"),("H","LHM",LINK_HR_MSG,LINK_HR_MSG+0x8000,"MMIO"),("H","HCC",LINK_HR+0x10000,LINK_HR+0x10000+HCC_BYTES,"setup:NODE_INIT dlvh"),("H","HCS",2**56-16384,2**56,"dlvh0"),("H","HH",UP_BASE,UP_BASE+2**HB*8,"zero16384 hhclr"),("H","NX",NEXT_LO,NEXT_LO+CAP_H*16,"dlvh"),("H","R2U",2**64-R2_UP_BYTES,2**64,"dlvh")]
_HERE=os.path.dirname(os.path.abspath(__file__))
_SYMS=set(re.findall(r"\.globl\s+(\w+)","".join(open(os.path.join(_HERE,"kernels",f)).read()for f in("leaf.s","route.s","secp_fe.s"))))
def _bld_ok(b):
 t=b.split()
 return bool(t)and all(x in _SYMS or x in("MMIO","analytic")or(x[:6]=="setup:"and isinstance(globals().get(x[6:]),int))for x in t)and(len(t)==1 or not{"MMIO","analytic"}&set(t))
def disjoint_va(iv):
 s_iv=sorted(iv)
 return(all(0<=a<b<=2**64 and valid_va(a)and valid_va(b-1)and((a<LOW)==(b-1<LOW))for a,b in s_iv)and all(s_iv[i][1]<=s_iv[i+1][0]for i in range(len(s_iv)-1)))
ROLE_N=dict(G=NODES,R=N_REL,H=N_HOME)
ROLE_B={k:sum(hi-lo for r_,_,lo,hi,_ in REGISTRY if r_==k)for k in ROLE_N}
PG,BLK=2**16,2**28
def _pt(k):
 iv=sorted((lo//PG*PG,cdiv(hi,PG)*PG,b=="MMIO",hi-lo)for r_,_,lo,hi,b in REGISTRY if r_==k);pg={};ok=True;l1=set();l2=set();l3=set()
 for lo,hi,m,_ in iv:
  for a in(lo,hi-PG):ok&=pg.setdefault(a,m)==m
  l1.update(range(lo>>52,((hi-1)>>52)+1));l2.update(range(lo>>40,((hi-1)>>40)+1))
  if lo%BLK:l3.add(lo>>28)
  if hi%BLK:l3.add((hi-1)>>28)
 T=2+len(l1)+len(l2)+len(l3)
 return ok,T,PG*T+sum(hi-lo-n for lo,hi,_,n in iv)
PT={k:_pt(k)for k in ROLE_N}
PROG_DELIV,NODE_STATE,MMU_ON,PT_DESC=(2**18//32)*2*MSG,2**12,64,3
NI_ITEM={k:PROG_DELIV+NODE_STATE+MMU_ON+8*2**10+PT_DESC*4096*PT[k][1]for k in ROLE_N}
check("c39",all(PT[k][0]for k in ROLE_N)and max(NI_ITEM.values())<=NODE_INIT,"tables "+"/".join(str(PT[k][1])for k in"GRH"))
QUEUE_MINE=cdiv(S_MINE_TOT,K)*(32*K+32)+cdiv(S_MINE_TOT,K2B)*(32*K2B+32)+NODES*F1*(32*K+32)+N_REL*F2*(32*K2B+32+64)+(36*N_HOME*CAND_CAP+4*NODES*MX+5*2*(N_HOME+NODES+2*(N_REL-1)))*64
QUEUE=max(NODES*q*(L_CHAIN*512+64)+(s//K)*(32*K+32)+(s//K2B)*(32*K2B+32)+NODES*F1*(32*K+32)+N_REL*F2*(32*K2B+32+64)for q,s in zip(Q1G,SG))+(36*N_HOME*CAND_CAP+4*ALLN+19*2*(N_HOME+NODES+2*(N_REL-1)))*64
PEAK=sum(ROLE_N[k]*(ROLE_B[k]+PT[k][2])for k in ROLE_N)+QUEUE
check("c40",QUEUE_MINE<=QUEUE and len(REGISTRY)==33 and[sum(r_==k for r_,_,_,_,_ in REGISTRY)for k in"GRH"]==[15,10,8]and all(lo<hi and _bld_ok(b)for _,_,lo,hi,b in REGISTRY)and all(disjoint_va([(lo,hi)for r_,_,lo,hi,_ in REGISTRY if r_==k])and ROLE_B[k]<2**56 for k in ROLE_N))
check("c41",YR_BASE-64+(YB+1)*YR_SLOT<=33*TIB and YROUNDS*YPC*512+64<=YR_SLOT and EB_START>=33*TIB and EB_END<=T34 and G34_END-T34-0x10000>=G34_CONTENT and CN_BASE+64*Q1GD<=CNM_BASE and CNM_BASE+64*NBATCH_MINE<=REQ_LO and REQ_END<=AXYX_BASE and XT_END<=H1L_BASE and H1L_BASE+4096*TIB<2**56 and R2_LOW_BYTES<LINK_HR and UP_BASE+2**HB*8==NEXT_LO and NEXT_LO+CAP_H*16<2**64-R2_UP_BYTES and HCC_BYTES<2**40-0x10000-16384)
check("c42",Q1GD<=2**38 and Q2<2**24 and YROUNDS*YPC==L_CHAIN and PTS_MINE_CHAIN<2**18 and NBATCH_MINE<2**44)
def least_at_or_above(d,x,scale):
 return P2(d)[0]>=x and P2(d-Fr(1,scale))[1]<x
def least_grid_up(x,scale):
 k=math.ceil(lg(x)*scale)
 for _ in range(20):
  d=Fr(k,scale)
  if least_at_or_above(d,x,scale):return d
  k+=1 if P2(d)[0]<x else-1
 raise AssertionError("grid bound not certified")
DECLARED=dict(time_log2=CLAIM_E,memory_log2_bytes=least_grid_up(PEAK,10**4),preprocessing_log2=least_grid_up(Fr(E_PRE,C),10**4),success=THRESH,per_leaf=Fr(cdiv(PER_LEAF_E.numerator*10**4,PER_LEAF_E.denominator),10**4))
GEN_ITEM=256*(DBL+MADD_)+2*(INV+3*MULX+2*SQRX+10*ADDX)+64
GEN_TEST,GEN_CTRL,NEGT,RED=256*2*13+64,256*6+64,30,9
X_ENTRY=MADD_+3*SQRX+11*MULX+2*ADDX+2*NEGT+17*RED+32
check("c43",GEN_ITEM+GEN_TEST+GEN_CTRL+X_ENTRY==GEN_MUL and GEN_MUL+4==GEN_MULX)
ENTRY_W=[14+(1 if j&1 else 0)+(3 if j%4==0 else 0)+(2 if(j&1)==0 else 0)+(2 if j==0 else 0)for j in range(MX)]
ENTRY_ITEM=sum(2*cdiv(w,2)for w in ENTRY_W)+3
TEMPLATE_ITEM=3*MULX+SQRX+24+7+20
check("c44",INV+3*MULX+2*SQRX+10*ADDX+4<CHAIN_ADD and 56*(2*MULX+ADDX+4)+TEMPLATE_ITEM+X_ENTRY<XSETUP and ENTRY_ITEM<=MX*ENTRY_X,"entries %d/8"%ENTRY_ITEM)
BARRIER_MSGS=2*(N_HOME+NODES+2*(N_REL-1))
BATCH_REST=11+3*9+2+3+6+3+MSG+4+SWAP+20+8+7+6
PH_STEP=8+8+2+15
BARRIER_MSG=7+DONE+DLV+RCV+WAKE+CTX
BARRIER_ITEM=ALLN*(3*13+1)+BARRIER_MSGS*BARRIER_MSG
check("c45",COIN_SCALAR+cdiv(4*YROUNDS,YB)+BATCH_REST<=BATCH_SETUP and G_GROUPS*(12*PH_STEP+16+4)<PHASE_CTRL_G and BARRIER>=7*BARRIER_ITEM and NBAR_COUNT<=NBAR_ALLOW and 9+MSG_LINK<=MSG)
DLVQ,DLVQ_DROP,REQH_C,FWDQ,SENDR,FWDR=10,7,17,11,10,15
PH4_HOPS={k:v+DONE+RCV+WAKE+CTX for k,v in dict(req_home_relay=cdiv(REQH_C,2)+DLVQ,req_relay_gen=FWDQ+DLVQ,rep_gen_relay=SENDR+DLVQ,rep_relay_home=FWDR+DLV).items()}
_h=list(PH4_HOPS.values())
check("c46",max(_h)<=MSG and RESOLVE==2*CAND_ITEM+1 and PROG_DELIV<NODE_INIT)
COIN_ITEM=2*(4*2+5+4)+16+2+1+2
P_COIN_FAIL=Fr(2**256-r+1,2**256)**2
INVC_ITEM,STOP_ITEM,HEADER_ITEM,REQ_ITEM,PH4_FIX,OUT_ITEM=2+1+6+6,5,17*4+9,16+12+600+16,EXPECT+DONE_C+CTX+32+24,16
check("c47",COIN_ITEM<=COIN_SCALAR and INVC_ITEM<=INV_CALL and STOP_ITEM<=STOP_CHECK and HEADER_ITEM<=HEADER and REQ_ITEM+GEN_MUL<=REQ_SERVE and 4*PH4_FIX*ALLN<=EXCH_ROW and OUT_ITEM<=OUTPUT and P_COIN_FAIL<Fr(1,2**254))
print("q1 %d  q2 %d  S 2^%.6f  mu %.9f  lo_run %.12f  p_gate %.9f  success >= %.12f"%(Q1GD,Q2,lg(SG_FULL),float(TFULL["mu"]),math.floor(LO_RUN*10**12)/10**12,float(P_GATE),math.floor(GATED_LO*10**12)/10**12))
print("E_mine 2^%.6f  E_run 2^%.6f insns (2^%.6f units)  E_total 2^%.6f  claim %d/10^15  per leaf %.4f  mem 2^%.7f  pre 2^%.7f  grids %s"%(lg(E_MINE),lg(E_RUN),lg(Fr(E_RUN,C)),float(E6),CLAIM_E*10**15,float(PER_LEAF_E),lg(PEAK),lg(Fr(E_PRE,C)),"/".join(format(float(least_grid_up(ROLE_B[k],10**4)),".4f")for k in"GRH")))
def _ps(a):return sum(PHI[g]*a[g]for g in range(G_GROUPS))
LEDGER=[ALLN*NODE_INIT,NODES*GEN_MUL,(2+2*G_GROUPS)*ALLN*MSG,E_MINE,NODES*(3*cdiv(Q1GD,M_G)+1)*COIN_SCALAR,NODES*MX*XSETUP,NODES*Q2*ENTRY_X+NODES*HEADER,NODES*XTBLD_ALL,NODES*TBLD,N_REL*((INIT3_BKT+TINIT_ENTRY)*F2+2)+NODES*(INIT3_BKT*F1+2),ZERO_INIT,HR_INIT,_ps(WGL),E_FIX_WORK,JROW_E-DEFJ,_ps(WGR),_ps(WGG),SXV,K2ROW_E,ELIST*HCHG_LS1,DEFV+DEFK2+DEFL+DEFJ,E_CAND_WORK,E_REQ_WORK,sum(PHI)*EXCH_ROW,NBAR_ALLOW*BARRIER,ALLN*PHASE_CTRL_G,OUTPUT]
LEDGER_CELLS=[format(float(least_grid_up(v,10**4)),".4f")for v in LEDGER]
E_RUN_6=least_grid_up(E_RUN,10**6)
E_RUN_U6=least_grid_up(Fr(E_RUN,C),10**6)
print("E_RUN_6 %.6f  E_RUN_U6 %.6f  E6 %.6f  QUEUE 2^%.4f  tables %s"%(float(E_RUN_6),float(E_RUN_U6),float(E6),lg(QUEUE),"/".join(str(PT[k][1])for k in"GRH")))
print("LEDGER_CELLS:",LEDGER_CELLS)
check("c48",len(LEDGER)==27 and sum(LEDGER)==E_RUN and P2(CLAIM_E)[0]>=Fr(E_TOTAL,C)and P2(E_RUN_6)[0]>=E_RUN and P2(E_RUN_U6)[0]>=Fr(E_RUN,C)and P2(E6)[0]>=E_TOTAL,"cells "+" ".join(LEDGER_CELLS))
STATED=dict(time_log2=Fr(123358839330617711,10**15),memory_log2_bytes=Fr(1380002,10**4),preprocessing_log2=Fr(1162652,10**4),success=Fr(39,100),per_leaf=Fr(671232,10**4))
_PF=open(os.path.join(_HERE,"proof.md")).read()
_SF=open(os.path.join(_HERE,"success.md")).read()
_SJT=open(os.path.join(_HERE,"submission.json")).read()
_SC=json.loads(_SJT)["claims"][0]["score"]
def cnt(text,tok):
 return len(re.findall(r"(?<![0-9])"+re.escape(tok)+r"(?![0-9])",text))
EC=[("123.358839330617711",3,1,2),("123358839330617711",1,0,1),("138.0002",2,0,1),("116.2652",2,0,1),("67.1232",2,1,1),("0.39",2,3,2),("0.3905",1,1,1),("0.390926575",1,1,1),("0.998908811125",1,1,1),("6.989415493",0,1,0),("134.679655",2,1,0),("124.713870",2,1,1),("133.324624",2,1,0),("128.610916",1,0,1),("1,313,229,200,424,960",1,0,1),("70,000,000,000",1,1,1),("16,769,280",1,0,0),("55.9255",1,0,0),("55.5869",1,0,0),("55.9997",1,0,0),("2048",10,1,2),("133025420141064816783379887723529591501",1,1,1),("1,185,778,776,658,177,683,654,854,049,792",1,1,0),("16,086,008,064",1,0,0),("112.0001",1,0,0),("85.2999",1,0,0),("101.0008",1,0,0),("125.3432",3,0,1),("99.6116",1,0,0),("86.0000",1,0,0),("94.1693",1,0,0),("120.9658",1,0,0),("115.9658",1,0,0),("109.0001",1,0,0),("125.0224",1,0,0),("86.3220",1,0,0),("134.5043",1,0,0),("107.1028",1,0,0),("128.2576",1,0,0),("130.9906",1,0,0),("129.2016",1,0,0),("43.4356",1,0,0),("112.3440",1,0,0),("93.3920",1,0,0),("55.5443",1,0,0),("112.5857",1,0,0),("111.1092",1,0,0),("103.0517",1,0,0),("109.2878",1,0,0),("102.0001",1,0,0),("20.0000",1,0,0)]
_cm=[k for(k,p,s,j)in EC if(cnt(_PF,k),cnt(_SF,k),cnt(_SJT,k))!=(p,s,j)]
check("c49",not _cm,"count mismatch: "+",".join(_cm))
check("c50",set(LEDGER_CELLS)<={k for(k,p,s,j)in EC},"unenforced: "+",".join(sorted(set(LEDGER_CELLS)-{k for(k,p,s,j)in EC})))
check("c51",Fr(int(_SC["numerator"]),int(_SC["denominator"]))==STATED["time_log2"]==DECLARED["time_log2"]and math.gcd(int(_SC["numerator"]),int(_SC["denominator"]))==1)
RATES={"2^-53.86":(D_WJ,1,1,0),"2^-10.88":(D_LAB,1,1,0),"2^-10.68":(3*D_WJ+D_LAB_SYM,1,0,0),"2^-93":(Fr((1<<322)*C0,(1<<192)*p),1,0,0),"2^-127.9":(D_BATCH,1,1,0),"2^-82":(TFULL["cs"],0,2,0),"2^-21.4":(TFULL["d_f"],0,2,0),"2^-37.7":(TFULL["lam_int"],0,2,0),"2^-33.6":(TFULL["d_p2"],0,2,0),"2^-17.5":(max(TFULL["d_cl"],TFULL["d_rq"]),0,4,0),"2^-18.5":(TFULL["d_rl"],0,2,0),"2^-25.5":(TFULL["g_def"],0,1,0),"2^-19":(max(TFULL["glob"],P_MINE_CAPS),0,3,0),"2^-28":(P_FAIL_ALL_POOLS,0,2,0),"2^-39":(P_SAME_A_ALL,0,1,0),"2^-140":(max(P_BAD_X,D_XMINE_ALL,P_MINE_ZERO_ALL),1,4,0),"2^-139":(P_EXC,1,1,0),"2^-204":(P_X_ORBIT_ALL_POOLS,0,1,0),"2^-254":(P_COIN_FAIL,1,1,0),"2^-102":(160*K_U_CLUST*K_U_CLUST*((POINTWISE_X_MAX/MX)**4)*KB_M,0,1,0)}
ANALYTIC={"2^-200"}
def rcnt(text,tok):
 return len(re.findall(re.escape(tok)+r"(?!\.?[0-9])",text))
_RT=[set(re.findall(r"2\^-[0-9.]*[0-9]",t))for t in(_PF,_SF,_SJT)]
_rbad=[k for k,(v,p_,s_,j_)in RATES.items()if not(v<P2(-Fr(k[3:]))[0]and(rcnt(_PF,k),rcnt(_SF,k),rcnt(_SJT,k))==(p_,s_,j_))]
check("c52",not _rbad and set().union(*_RT)<=set(RATES)|ANALYTIC,"bad: "+",".join(_rbad+sorted(set().union(*_RT)-set(RATES)-ANALYTIC)))
check("c53",K_GATE==133025420141064816783379887723529591501 and CAP_H==1313229200424960 and PSEL_POOL==262152 and EXC_WORK==11778206816208 and E_RUN_6==Fr(134679655,10**6)and E_RUN_U6==Fr(124713870,10**6)and E6==Fr(133324624,10**6)and STATED==DECLARED)
print("FAIL: "+", ".join(FAILED)if FAILED else"PASS")
if __name__=="__main__":
 sys.exit(1 if FAILED else 0)

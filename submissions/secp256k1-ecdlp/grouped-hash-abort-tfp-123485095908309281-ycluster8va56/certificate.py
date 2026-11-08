import functools,json,math,os,re,sys
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
p=2**256-2**32-977
r=int("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141",16)
C0=2**32+977
BETA_P=int("7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee",16)
C=1000
NODES=2**68
M_G=1024
N_CLUST=NODES//M_G
F1,F2=2**42,2**42
N_REL,N_HOME=F1,F1*F2
ALLN=NODES+N_REL+N_HOME
SIDE=F1
MX=8
L_CHAIN=2**21
Q2=MX*L_CHAIN
D1,D2,D3,D4=42,42,64,22
KEYBITS=64+D4
PART=64+D1
KBF=3195
K1,K2B=16,16
THRESH=Fr(39,100)
LEAVES=6
FWD_INC,FWD_FIX=Fr(259,4),7
def fwd_count(q2):
 assert q2%MX==0
 return FWD_FIX+(259*q2)//4
BWD_INC,BWD_FIX,BWD_YSW=Fr(2013,2),7,0
def bwd_count(q2):
 assert q2%MX==0
 return BWD_FIX+(2013*q2)//2
FIXSP=dict(SHX=27,FDX=6,F5=2,FS=3,FADDW=2,BLK=16)
FIX_EXTRA=dict(SHX=3,FDX=3,F5=3,FS=4,FADDW=4,BLK=5)
EMIT=5
EMIT_INC=LEAVES*EMIT
KEY_INC,KEY_FIX=Fr(257,2),5
KEY_W3=Fr(2,1)
K=K1
FLUSH=10
FLUSH2_X,FLUSH2_D=10,8
JPATH,SJOIN=1,1
INV=22549
TINIT_ENTRY=3
FIN_DEST,FIN_SENT,FIN_HELD,FIN_FIX=15,4,2,2
FINH_DEST,FINH_SENT=16,5
H2S_ARR,H2S_ONE=6,1
H2S_FULL= -1
FLUSH1_X,FIN1_X=8,14
H2S_EMPTY=6
H2S_ARR_FIX=14+H2S_ONE
def h2s_count(n):
 if n==K:return H2S_FULL+H2S_ARR*K
 if n==0:return H2S_EMPTY
 return 14+H2S_ONE*(n%2)+H2S_ARR*n
check("c1",h2s_count(K)==95 and max(h2s_count(n)-H2S_ARR*n for n in range(K))==H2S_ARR_FIX and H2S_EMPTY<=H2S_ARR_FIX)
INIT3_BKT=4
check("c2",KEY_INC==8+52+KEY_W3+6*(6+EMIT)+Fr(1,2)and KEY_FIX==5)
check("c3",BWD_INC==Fr(2013,2)and Q2%MX==0 and BWD_FIX==7 and bwd_count(Q2)==BWD_FIX+(2013*Q2)//2)
check("c4",FWD_INC==Fr(259,4)and Q2%MX==0 and FWD_FIX==7 and fwd_count(Q2)==FWD_FIX+(259*Q2)//4)
check("c5",(4+2+1+EMIT)-(6+EMIT)==JPATH and SJOIN==1)
check("c6",EMIT==1+1+2+1)
DONE=2
RCV,WAKE,WY,EXPECT,SWAP,CTX=1,5,21,34,4+4*15,0
DONE_C=DONE+RCV+WAKE
RCV1,FIN_CNT,H2S_DONE=4,6,4
EXPECTH=68
RDONE=1+H2S_DONE
DONE_C1=DONE+RCV1+WAKE
DLV=6
MSG_LINK=DONE_C+DLV+CTX
FINAL=FIN_DEST+FIN_SENT+DONE_C1+H2S_ARR_FIX+FIN_CNT
check("c7",FLUSH1_X<=FLUSH+DONE and FIN1_X-H2S_FULL<=FIN_DEST+FIN_SENT+DONE and FLUSH2_X<=FLUSH and FLUSH2_D<=FLUSH and FIN_HELD<=FINH_SENT and DONE==2 and(DONE_C,DONE_C1,MSG_LINK,FINAL,RDONE)==(8,11,14,51,5)and H2S_ARR_FIX==15)
check("c8",5<=WAKE and 11+2*5<=WY and(SWAP,CTX)==(64,0)and 3+5+26<=EXPECT and EXPECTH==2*EXPECT and FIN_FIX==2)
GEN_MUL=645410
GEN_MULX=GEN_MUL+4
YR_FIX,YR_BLK,YB=22556,42697,256
YPC=33
YROUNDS=cdiv(L_CHAIN,YPC)
CHAIN_ADD=2**15
CHAIN_SET=2*GEN_MUL+17*CHAIN_ADD
YSND,YDLV,YDLVE=24,32,23
YRING_TOT=YSND+3+(M_G-2)*(YDLV+WY)+YDLVE+WY
def y_group(q):
 q_gen=cdiv(q,M_G)
 gen_y=YROUNDS*(cdiv(q_gen,YB)*YR_FIX+q_gen*YR_BLK)+q_gen*CHAIN_SET
 ring_y=cdiv(q*L_CHAIN*YRING_TOT,M_G)
 return gen_y+ring_y
XT_SLOTS=7
XT_W,XT_WTOP=37,34
XT_ENTRIES_SINGLE=6*2**XT_W+2**XT_WTOP
XTBLD_SINGLE=4+XT_SLOTS*11+20*XT_ENTRIES_SINGLE
XTBLD_PAIRED=6*(4+2*11)+20*2*(4*2**43+2*2**42)
XTBLD_ALL=4*XTBLD_SINGLE+22*XTBLD_PAIRED
XSETUP=2**17
ENTRY_X=18
COIN_SCALAR=64
BATCH_SETUP=1344
INV_CALL=16
NODE_INIT=2**30
HEADER=80
BT_BASE,BT_W=0xfff8<<48,43
TBLD=4+7*11+20*(6*2**BT_W+2**(256-5*BT_W))
MSG=128
CAND_CAP=2**27
RESOLVE=1502275
REQ_CAP=96*2**33
REQ_SERVE=2**10+GEN_MUL
RELAY_CAP=2**49
BARRIER=2**12*ALLN
OUTPUT=2**20
DELTA_F=2**52
FIX_BATCH=sum(FIXSP.values())*(Q2//2)
P1=Fr(6*(r-1)+1,r*r)
TWO_BOUND=Fr(2,r-1)
SHX_BAD=2**64*2**33
FOLD_BAD=2**66*C0
SQ_BAD=C0*C0+2*C0
FS_P1_BAD=C0
QP=Q2//2
SITES=[(3*QP,SHX_BAD,1),(QP,FOLD_BAD,1),(3*QP,SHX_BAD,1),(QP,FOLD_BAD,12),(QP,FOLD_BAD,1),(3*QP,SHX_BAD,6),(QP,FOLD_BAD,6),(Q2,SQ_BAD,12),(Q2,FOLD_BAD,12),(2*Q2,FS_P1_BAD,6),(9*Q2,SHX_BAD,8)]
SITES_PAIR=dict(SHX=27,FDX=6,F5=2,FS_FA=4)
D_SITES=sum(Fr(n*bad*mult,r-1)for n,bad,mult in SITES)
FIX_ROW=C0*2**192
FIX_FOLD=(2**66+1)*2**128*C0
FIX_SQ=2**193
FIX_LK6=15*C0*2**128
FIX_LK7=21*C0*2**128
FIX_LK8=28*C0*2**128
FIX_BLK=30*C0*2**128
FIX_SUB=C0*2**192
FIXUPS=[("BLK",QP,FIX_LK7,12),("SHX",3*QP,FIX_ROW,1),("FDX",QP,FIX_FOLD,1),("SHX",3*QP,FIX_ROW,1),("FDX",QP,FIX_FOLD,12),("FDX",QP,FIX_FOLD,1),("FS",QP,FIX_SUB,6),("SHX",3*QP,FIX_ROW,6),("FDX",QP,FIX_FOLD,6),("BLK",QP,FIX_LK8,6),("F5",Q2,FIX_SQ,12),("BLK",Q2,FIX_LK7,12),("BLK",Q2,FIX_LK6,12),("FDX",Q2,FIX_FOLD,12),("FS",Q2,FIX_SUB,6),("FADDW",Q2,FIX_SUB,6),("BLK",2*Q2,FIX_BLK,6),("BLK",Q2,FIX_LK7,8),("BLK",2*Q2,FIX_LK6,8),("SHX",9*Q2,FIX_ROW,8)]
E_FIX=sum(Fr(n*bad*mult,r-1)for _,n,bad,mult in FIXUPS)
KB=functools.lru_cache(None)(lambda fixed:Fr(6*2**(256-fixed)*(KBF+1),(r-1)*KBF))
RHO=Fr(r,r-1)**2
CLUST_BATCH=M_G*LEAVES*Q2
FIXED_LEAVES=2*CLUST_BATCH
PARTNER_PAIRS=CLUST_BATCH*(CLUST_BATCH-1)
D_YGEN=Fr(YROUNDS*(YPC+2*16+2),r-1)
D_BATCH=D_SITES+Q2*Fr(6,r-1)+D_YGEN
D_WJ=TWO_BOUND*23*2**192
D_LAB=TWO_BOUND*1305*2**(128+PART)
FA,FB=int("114ca50f7a8e2f3f657c1108d9d44cfd9",16),int("e4437ed6010e88286f547fa90abfe4c3",16)
BOX=(p>>PART)+1
WIDE=654*(p>>128)+4
check("c9",FA*FA-FA*FB+FB*FB==p and(FB*BETA_P-(FA-FB))%p==0 and((FA-FB)*BETA_P+FA)%p==0 and(BOX+WIDE+FA)*(BOX+WIDE+2*FA-FB)*2**(2*PART)<p*p*(1+Fr(1,KBF))and(2*(BOX+WIDE)+FA)*(2*(BOX+WIDE)+2*FA-FB)*2**(2*PART)<4*p*p*(1+Fr(1,KBF)))
def isum(*rs):return sum(x[0]for x in rs),sum(x[1]for x in rs)
TFP_TERMS=((0,1),(0,7),(0,56))
RAW=isum(*TFP_TERMS);COR=isum(RAW,(-1,0))
U3,UC3=(RAW,RAW,(-2-2*RAW[1],-2*RAW[0])),(COR,COR,(-2-2*COR[1],-2*COR[0]))
check("c10",RAW==(0,64)and COR==(-1,64)and U3[2]==(-130,0)and UC3[2]==(-130,2))
WR=[[isum(UC3[m],UC3[(m+t)%3],U3[(m+2*t)%3],(0,2))for t in range(3)]for m in(0,1)]
NR=[(-WR[0][t][1]-WR[1][t][1]-1,-WR[0][t][0]-WR[1][t][0]+1)for t in range(3)]
LABEL_ERR=[(min(x[0]for x in rs),max(x[1]for x in rs))for rs in(WR[0],WR[1],NR)]
check("c11",LABEL_ERR==[(-132,194),(-132,194),(-389,264)]and sum(hi-lo for lo,hi in LABEL_ERR)==1305 and max(hi-lo for lo,hi in LABEL_ERR)+1==654 and WIDE==654*(p>>128)+4)
check("c12",D_WJ<P2(Fr(-5847,100))[0],"2^%.4f"%lg(D_WJ))
check("c13",sum(n for n,b,_ in SITES if b==SHX_BAD)==SITES_PAIR["SHX"]*QP and sum(n for n,b,_ in SITES if b==FOLD_BAD)==SITES_PAIR["FDX"]*QP and sum(n for n,b,_ in SITES if b==SQ_BAD)==SITES_PAIR["F5"]*QP and sum(n for n,b,_ in SITES if b==FS_P1_BAD)==SITES_PAIR["FS_FA"]*QP and len(SITES)==11)
FIX_RUNS={m:sum(n for ms,n,_,_ in FIXUPS if ms==m)for m in FIXSP}
check("c14",FIX_RUNS=={m:FIXSP[m]*QP for m in FIXSP}and sum(FIXSP.values())==56)
check("c15",FIX_BATCH*C0<2**63 and FIX_BATCH*C0>>32==FIX_BATCH+106,"2^%.3f"%lg(FIX_BATCH*C0))
check("c16",D_LAB<P2(Fr(-1065,100))[0]and 3*D_WJ+TWO_BOUND*1959*2**(128+PART)<P2(Fr(-1006,100))[0],"2^%.3f 2^%.4f"%(lg(D_LAB),lg(3*D_WJ+TWO_BOUND*1959*2**(128+PART))))
check("c17",D_BATCH<P2(Fr(-1274,10))[0],"2^%.3f"%lg(D_BATCH))
G_GROUPS=2048
G_EXP=24
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
HB=47
HCHG_NEW,HCHG_VISIT,HCHG_END,HCHG_VER=4,4,0,5
HCHG_LS0,HCHG_LSA=2,1
HCHG_LS1=5+HCHG_LSA
check("c18",HCHG_NEW==4 and HCHG_VISIT==4 and HCHG_END==0 and HCHG_VER==5 and HCHG_LS1==6)
DLVC,DLVS1=2,2
DLVHB=31
DLV_RESTORE=5
DLVH_SEG,DLVH_SLOW_FIX=7,17+DLV_RESTORE
def dlvh_seg(n,m):
 if m==0:return 1
 head=n%2;rest=m-head
 return DLVH_SEG+head+DLVC*m+DLVS1*(head+rest%2)
def dlvh_copy(m1,m2,K_=K2B):
 if m1==0 and m2==K_:return DLVHB
 return DLVH_SLOW_FIX-(m1>0)+dlvh_seg(0,m1)+dlvh_seg(m1,m2)+dlvh_seg(m1+m2,K_-m1-m2)
DLVHS=max(dlvh_copy(a,b)for a in range(K2B+1)for b in range(K2B+1-a)if(a,b)!=(0,K2B))
BLOCKH=FLUSH+DONE_C+DLVHB
FINALH=FINH_DEST+FINH_SENT+DONE_C+DLVHS+MSG_LINK+(DLVHS-DLVHB)
check("c19",DLVHB==31 and DLVHS==84 and BLOCKH==FLUSH+DONE_C+DLVHB and FINALH==180)
ZW=16384
ZERO_PASS=2*(ZW//32)+2
STOP_CHECK=5
NBAR_COUNT=4+19*G_GROUPS+1
NBAR_ALLOW=40000
PHASE_CTRL_G=2**20
KB_BUCKET=KB(D1+D2+HB)
KB64=KB(64)
R2_UP_TIB=682
CAP_H=(((R2_UP_TIB+4095)*(2**40))//512)*16
D_LAB_SYM=TWO_BOUND*(3*653)*2**(128+PART)
P1_GOOD=Fr(6*(r-1),r*r)*(1-3*D_WJ-D_LAB_SYM)
K_U_CLUST=LEAVES*M_G*MX
SAMEX_FAC=Fr(CLUST_BATCH-K_U_CLUST,CLUST_BATCH-1)*Fr(r,r-1)+Fr(24*(K_U_CLUST-1),CLUST_BATCH-1)
EF_Q=E_FIX*C0/(2**32*(DELTA_F-FIX_BATCH))
KB_M=KB(KEYBITS+D1+D2)
SIP=RHO*Fr(5,2)*(1+11*(D_WJ+D_LAB))*KB_M
def success_prefix(q1p,qs=None):
 n0p=LEAVES*q1p*Q2
 Sp=NODES*n0p
 clust_n0=M_G*n0p
 gamma=Fr(Sp*Sp-N_CLUST*clust_n0*clust_n0,2)
 mu_raw=gamma*P1
 mu=gamma*P1_GOOD
 n_y=q1p*L_CHAIN
 deg=2*LEAVES*(n_y+M_G*Q2)*Sp
 cs=gamma*deg*(P1*P1+Fr(36,(r-1)**2))
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
 e_hi=exp_neg_hi(mu)
 d_batch,d_wj,d_lab=D_BATCH,D_WJ,D_LAB
 d_f=q1p*EF_Q+q1p*d_batch
 m_p2=Sp*KB(D1+D2)+1
 cov_p2=Fr(Sp*K_U_CLUST*L_CHAIN,r-1)*(KB(D1+D2)**2)
 v_p2=K_U_CLUST*Sp*KB(D1+D2)+cov_p2+1
 d_p2=v_p2/Fr(CAP_H-FIXED_LEAVES-m_p2)**2
 kb_m=KB_M
 lam_int=RHO*(Sp+FIXED_LEAVES)*KB(D1+D2+D3+D4)
 sip=SIP
 alpha_cross=2*K_U_CLUST*Sp+SAMEX_FAC*Sp*(CLUST_BATCH-1)+PARTNER_PAIRS+CLUST_BATCH*CLUST_BATCH
 other_h=RHO*(Fr(Sp)*Sp*KB(D1+D2)*kb_m+alpha_cross*kb_m+2*L_CHAIN*KB(D1+D2))+sip*Sp
 d_cl=other_h/(CAND_CAP-4)
 other_g=RHO*(Fr(n0p)*Sp*kb_m*2+2*LEAVES*MX*Sp*kb_m)+24*LEAVES*MX*n0p*kb_m
 d_rq=other_g/(REQ_CAP-2*L_CHAIN-8)
 rl_cap_cr=RELAY_CAP//4
 rl_cr=18*FIXED_LEAVES*Fr(Sp)*kb_m
 rl_ediff=RHO*Fr(Sp)*Sp*KB_D1*kb_m
 rl_fdiff=LEAVES*LEAVES+Fr(2*K_U_CLUST*K_U_CLUST,N_CLUST)+2*K_U_CLUST*K_U_CLUST*kb_m+2*CLUST_BATCH*Fr(Sp)*kb_m+2*CLUST_BATCH*CLUST_BATCH*RHO*kb_m
 v_rl=648*(rl_ediff*rl_fdiff+CLUST_BATCH*Fr(Sp)*RL_SAME_PL)
 d_rl=rl_cr/rl_cap_cr+v_rl/(Fr(RELAY_CAP-18*FIXED_LEAVES-rl_cap_cr)-M_VG_MAX)**2
 per_leaf=d_batch+d_f+d_p2+d_cl+d_rq+d_rl
 rows_=mu_raw*(1+Fr(1,r))*(Fr(1,r)+lam_int+2*per_leaf)
 g_def=NODES*q1p*d_batch
 glob=(NODES*(3*q1p+MX+1)*Fr(1,2**254)+NODES*MX*Fr(1,r)+N_CLUST*((M_G*MX)*(M_G*MX-1)//2)*Fr(6,r-1)+g_def)
 lo=1-e_hi-cs-rows_-glob
 return dict(mu=mu,lo=lo,cs=cs,rows=rows_,glob=glob,Sp=Sp,n0p=n0p,d_p2=d_p2,m_p2=m_p2,d_f=d_f,lam_int=lam_int,d_cl=d_cl,d_rq=d_rq,d_rl=d_rl,g_def=g_def)
def grouped_full_lo(q):
 return success_prefix(q)["lo"]
Q1GD=17170000000
Q1G,MQ1=group_sizes(Q1GD)
KB_D1=KB(D1)
RL_SAME_PL=RHO*Fr(CLUST_BATCH-1,2)*KB_D1*KB_M+12*K_U_CLUST*KB_M+SIP
_SG_G=[NODES*LEAVES*q*Q2 for q in Q1G]
SG_MAX=max(_SG_G)
SPROD_MAX=max(_SG_G[g]*(NODES*LEAVES*MQ1[g+1]*Q2)for g in range(G_GROUPS))
M_VG_MAX=18*(RHO*Fr(SPROD_MAX)*KB_D1*KB_M+Fr(SG_MAX)*RL_SAME_PL)
N0G=LEAVES*Q1GD*Q2
SG_FULL=NODES*N0G
LO_RUN=grouped_full_lo(Q1GD)
GATE_BITS=128
GATE_INSNS=256
THRESH_GATE=Fr(39000001,10**8)
K_GATE=cdiv(THRESH_GATE.numerator*LO_RUN.denominator*(1<<GATE_BITS),THRESH_GATE.denominator*LO_RUN.numerator)
P_GATE=Fr(K_GATE,1<<GATE_BITS)
GATED_LO=P_GATE*LO_RUN
check("c20 %d"%Q1GD,LO_RUN>=Fr(9983,10000),"lo_run %.12f"%(math.floor(LO_RUN*10**12)/10**12))
check("c21",0<K_GATE<(1<<GATE_BITS)and GATED_LO>=THRESH_GATE>THRESH and Fr(K_GATE-1,1<<GATE_BITS)*LO_RUN<THRESH_GATE and 2*COIN_SCALAR+8+3+4<=GATE_INSNS,"p_gate %.9f gated_lo %.12f"%(float(P_GATE),math.floor(GATED_LO*10**12)/10**12))
check("c22",all((MQ1[g]**32)*(G_GROUPS**G_EXP)<=(Q1GD**32)*(g**G_EXP)<((MQ1[g]+1)**32)*(G_GROUPS**G_EXP)for g in range(1,G_GROUPS)))
TPRE=[success_prefix(m)for m in MQ1[: -1]]
def ceil_256(x):
 x=min(Fr(1),max(Fr(0),x))
 return Fr(cdiv(x.numerator<<256,x.denominator),1<<256)
PHI=[Fr(1)]+[ceil_256(1-t["lo"])for t in TPRE[1:]]
TFULL=success_prefix(Q1GD)
check("c23 %.6f"%lg(SG_FULL),SG_FULL==NODES*6*Q1GD*Q2==TFULL["Sp"])
check("c24 %.9f"%float(TFULL["mu"]),TFULL["mu"]==Fr(SG_FULL*SG_FULL-N_CLUST*(M_G*N0G)**2,2)*Fr(6*(r-1),r*r)*(1-3*D_WJ-TWO_BOUND*1959*2**(128+PART))and GATED_LO>=THRESH)
check("c25",TFULL["cs"]<Fr(1,2**62),"2^%.4f"%lg(TFULL["cs"]))
check("c26",TFULL["lam_int"]<P2(Fr(-388,10))[0],"2^%.3f"%lg(TFULL["lam_int"]))
check("c27 %.3f"%lg(Q1GD*E_FIX),Q1GD*E_FIX<P2(Fr(35,1))[0]and DELTA_F>2**17*Q1GD*E_FIX+FIX_BATCH)
check("c28",TFULL["d_f"]<P2(Fr(-184,10))[0],"2^%.3f"%lg(TFULL["d_f"]))
check("c29",TFULL["d_p2"]<P2(Fr(-199,10))[0],"2^%.3f"%lg(TFULL["d_p2"]))
check("c30",TFULL["d_cl"]<P2(Fr(-186,10))[0]and TFULL["d_rq"]<P2(Fr(-158,10))[0]and TFULL["d_rl"]<P2(Fr(-274,10))[0],"2^%.3f 2^%.3f 2^%.3f"%(lg(TFULL["d_cl"]),lg(TFULL["d_rq"]),lg(TFULL["d_rl"])))
check("c31",TFULL["g_def"]<P2(Fr(-255,10))[0]and TFULL["glob"]<P2(Fr(-255,10))[0],"2^%.3f"%lg(TFULL["glob"]))
check("c32",PHI[0]==1 and all(0<=x<=1 for x in PHI)and all(PHI[g]>=PHI[g+1]for g in range(G_GROUPS-1))and all(PHI[g]>=min(Fr(1),max(Fr(0),1-TPRE[g]["lo"]))for g in range(1,G_GROUPS)),"PHI[1] %.4f PHI[16] %.4f PHI[31] %.4f"%(float(PHI[1]),float(PHI[16]),float(PHI[31])))
HOME_MAX=CAP_H
check("c33",HOME_MAX==CAP_H<2**59,"2^%.4f"%lg(HOME_MAX))
PAIRM=KB_BUCKET+Fr(1,r)
PAIRK2=KB(D1+D2+64)+Fr(1,r)
VISIT=HCHG_VISIT+Fr(HCHG_END,2)
PJ=Fr(2,3)+2*D_WJ+Fr(1,r-1)
check("c34",VISIT==4)
SG=[NODES*LEAVES*q*Q2 for q in Q1G]
WG,SAMEX,WGJ,WGK2,WGL,WGR,WGG=[],[],[],[],[],[],[]
EXCH_ROW=ALLN*4*2**10
for _g,_q in enumerate(Q1G):
 _s=SG[_g]
 _leaf=NODES*_q*(fwd_count(Q2)+INV+INV_CALL+BATCH_SETUP)+NODES*y_group(_q)+int(NODES*_q*bwd_count(Q2))+int(NODES*_q*((KEY_INC-EMIT_INC)*Q2+KEY_FIX))
 _route=_s*(EMIT*K+FLUSH+DONE_C1)//K+_s*H2S_ARR+(_s//K)*H2S_FULL+NODES*(FINAL*F1+FIN_FIX)+_s*BLOCKH//K2B+N_REL*(FINALH*F2+FIN_FIX)+N_HOME*EXPECTH+N_REL*(EXPECT+RDONE)+NODES*(INIT3_BKT*F1+1)
 _gp=NODES*LEAVES*MQ1[_g]*Q2
 _samex=N_CLUST*_q*CLUST_BATCH*(CLUST_BATCH-1)//2
 _main=_s*_gp+_s*(_s-1)//2-_samex
 _evis=_main*PAIRM
 _gfix=_s*HCHG_NEW+ALLN*STOP_CHECK
 WGK2 .append(HCHG_VER*_main*PAIRK2)
 _group=_gfix+VISIT*_evis+WGK2[-1]
 _jkey=JPATH*_s*PJ+_s//3*SJOIN
 WGJ.append(_jkey);WGL.append(_leaf);WGR.append(_route);WGG.append(_gfix+VISIT*_evis)
 WG.append(_leaf+_route+_group+_jkey+EXCH_ROW)
 SAMEX.append(_samex)
check("c35",sum(Q1G)==Q1GD and min(Q1G)>0 and all(s%K==0 for s in SG))
check("c36",NODES*2**50<=2**118 and 2**126+2**118<2**127 and NODES<2**69 and max(Q1G)*LEAVES*Q2//K+1<2**50)
ZERO_INIT=N_HOME*(2**HB//ZW*ZERO_PASS+1)
check("c37",ZERO_PASS==1026 and 2**HB%ZW==0 and(2**64-2**52)%256==0)
DLVH0=15
HR_INIT=N_HOME*(INIT3_BKT+1+DLVH0)
E_PRE=ALLN*NODE_INIT+NODES*(3*cdiv(Q1GD,M_G)+MX+1)*COIN_SCALAR+NODES*MX*(GEN_MULX+XSETUP)+NODES*XTBLD_ALL+NODES*Q2*ENTRY_X+NODES*HEADER+N_REL*(INIT3_BKT*F2+1+TINIT_ENTRY*F2+1)+NODES*(INIT3_BKT*F1+2)+ZERO_INIT+HR_INIT+NODES*TBLD
_KB170=KB(KEYBITS+D1+D2)
_SIPT=RHO*Fr(5,2)*(1+11*(D_WJ+D_LAB))
_ACROSS_FULL=2*K_U_CLUST*SG_FULL+SAMEX_FAC*SG_FULL*(CLUST_BATCH-1)+PARTNER_PAIRS+CLUST_BATCH*CLUST_BATCH
ELIST=RHO*(Fr(SG_FULL)*SG_FULL*_KB170+_ACROSS_FULL*_KB170)+_SIPT*_KB170*SG_FULL
EK2=RHO*(_ACROSS_FULL*KB64)+_SIPT*KB64*SG_FULL
SXK2=SAMEX_FAC*PAIRM*sum(SAMEX)
SXV=SAMEX_FAC*VISIT*PAIRM*sum(SAMEX)
BRACKET=RHO*(Fr(SG_FULL)*SG_FULL*KB(D1+D2)+_ACROSS_FULL)+_SIPT*SG_FULL
DEFV=VISIT*D_BATCH*(HOME_MAX*SG_FULL+sum(SAMEX))
DEFK2=HCHG_VER*2*D_BATCH*BRACKET+HCHG_VER*D_BATCH*HOME_MAX*SG_FULL
DEFL_PAIRS=2*D_BATCH*BRACKET+D_BATCH*HOME_MAX*SG_FULL
DEFL=HCHG_LS1*DEFL_PAIRS
DEFJ=JPATH*D_BATCH*SG_FULL
E_FIX_WORK=NODES*Q1GD*E_FIX*max(FIX_EXTRA.values())
E_CAND_WORK=(ELIST+DEFL_PAIRS)*RESOLVE
E_REQ_WORK=(ELIST+DEFL_PAIRS)*(2*REQ_SERVE)
E_SHARED_REST=E_PRE+NODES*GEN_MUL+2*ALLN*MSG+NBAR_ALLOW*BARRIER+ALLN*PHASE_CTRL_G+OUTPUT+G_GROUPS*2*ALLN*MSG
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
LOW,UP_BASE=2**56,2**64-2**52
def valid_va(a):return 0<=a<2**64 and(a>>56)in(0,255)
TIB=2**40
XT_PAIR6_BYTES=14*2560*TIB
XT_DEL32_BYTES=4*(2**XT_W<<8)
AXYX_BASE=4096*TIB
AXYX_HI=20480*TIB
XT_BASE=24576*TIB
XT_END=XT_BASE+XT_PAIR6_BYTES
DEL_BASE=96*TIB
CN_BASE=84*TIB
check("c38",BT_BASE==2**64-2**51 and DEL_BASE%(2**45)==0 and AXYX_BASE%(2**49)==0 and AXYX_HI%(2**48)==0 and XT_BASE%(2**48)==0 and valid_va(XT_END-1)and valid_va(DEL_BASE)and all(valid_va(BT_BASE|((w<<8)%2**64))for w in(0,1,2**43-1,2**64-1)))
YR_SLOT=2**31
YR_BASE=32*TIB+YR_SLOT
EB_END=33*TIB+2**34
EB_START=EB_END-Q2*512
T34=34*TIB
LINK_GEN,LINK_GEN_MSG,PORTY_GEN=T34,T34+0x1000,T34+0xB000
REQ_LO=T34+2**40-64
REQ_END=REQ_LO+REQ_CAP*64+64
LINK_HR=2**52-2**40
LINK_HR_MSG=LINK_HR+0x1000
NEXT_BASE=2**64-3*2**50
R2_UP_BYTES=R2_UP_TIB*TIB
R2_LOW_BYTES=CAP_H*32-R2_UP_BYTES
HCC_BYTES=2**20+8*2048+4096+CAND_CAP*528
G34_CONTENT=2**20+8*2048+MX*6*XT_SLOTS*64+2048+MX*64+(MX+2*cdiv(Q1GD,M_G)+Q1GD+1)*32+2**30+YB*2304
REGISTRY=[("G","LC1",0,F1*8,"init3"),("G","YR",YR_BASE-64,YR_BASE-64+(YB+1)*YR_SLOT,"yround"),("G","EB",EB_START,EB_END,"setup:ENTRY_X fwd"),("G","LG",LINK_GEN,LINK_GEN+4096,"MMIO"),("G","LGM",LINK_GEN_MSG,LINK_GEN_MSG+4096,"MMIO"),("G","PY",PORTY_GEN,PORTY_GEN+512,"MMIO"),("G","G34",T34+0x10000,REQ_LO,"setup:NODE_INIT setup:XSETUP setup:COIN_SCALAR ydlv0"),("G","REQ",REQ_LO,REQ_END,"dlvq"),("G","CN",CN_BASE,CN_BASE+64*Q1GD,"setup:BATCH_SETUP"),("G","DEL",DEL_BASE,DEL_BASE+XT_DEL32_BYTES,"tbld"),("G","H1L",2**51,2**52,"init3"),("G","AXY",AXYX_BASE,XT_BASE,"tbld"),("G","XT6",XT_BASE,XT_END,"tbld"),("G","H1U",UP_BASE,UP_BASE+2**51,"init3"),("G","BT",BT_BASE,2**64,"tbld"),("R","RCS",0x10000,0x10000+2**20+8*2048,"setup:NODE_INIT h2s0"),("R","LC2",2**45,2**45+F2*8,"init3"),("R","RC2",3*2**45,3*2**45+F2*8,"tinit"),("R","E2P",2**47,2**47+2*F2*8,"analytic"),("R","DPW",2**48-4096,2**48,"h2s0"),("R","RCB",2**52,2**52+0x10000+RELAY_CAP*64,"setup:NODE_INIT setup:EXCH_ROW dlvq"),("R","LR",LINK_HR,LINK_HR+4096,"MMIO"),("R","LRM",LINK_HR_MSG,LINK_HR_MSG+4096,"MMIO"),("R","RCT",2**52-8192,2**52-4096,"dlvq"),("R","H2U",UP_BASE,2**64,"init3"),("H","R2L",0,R2_LOW_BYTES,"dlvh"),("H","LH",LINK_HR,LINK_HR+4096,"MMIO"),("H","LHM",LINK_HR_MSG,LINK_HR_MSG+4096,"MMIO"),("H","HCC",LINK_HR+0x10000,LINK_HR+0x10000+HCC_BYTES,"setup:NODE_INIT dlvh"),("H","HCS",2**52-16384,2**52,"dlvh0"),("H","HH",UP_BASE,UP_BASE+2**HB*8,"zero16384"),("H","NX",NEXT_BASE,NEXT_BASE+CAP_H*16,"dlvh"),("H","R2U",2**64-R2_UP_BYTES,2**64,"dlvh")]
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
  l1 .update(range(lo>>52,((hi-1)>>52)+1));l2 .update(range(lo>>40,((hi-1)>>40)+1))
  if lo%BLK:l3 .add(lo>>28)
  if hi%BLK:l3 .add((hi-1)>>28)
 T=2+len(l1)+len(l2)+len(l3)
 return ok,T,PG*T+sum(hi-lo-n for lo,hi,_,n in iv)
PT={k:_pt(k)for k in ROLE_N}
PROG_DELIV,NODE_STATE,MMU_ON,PT_DESC=(2**17//32)*2*MSG,2**12,64,3
NI_ITEM={k:PROG_DELIV+NODE_STATE+MMU_ON+8*2**10+PT_DESC*4096*PT[k][1]for k in ROLE_N}
check("c39",all(PT[k][0]for k in ROLE_N)and max(NI_ITEM.values())<=NODE_INIT,"tables "+"/".join(str(PT[k][1])for k in"GRH"))
QUEUE=max(NODES*q*(L_CHAIN*512+64)+(2*(s//K)+NODES*F1+2*N_REL*F2)*544 for q,s in zip(Q1G,SG))+(36*N_HOME*CAND_CAP+4*ALLN+19*2*(N_HOME+NODES+2*(N_REL-1)))*64
PEAK=sum(ROLE_N[k]*(ROLE_B[k]+PT[k][2])for k in ROLE_N)+QUEUE
check("c40",len(REGISTRY)==33 and[sum(r_==k for r_,_,_,_,_ in REGISTRY)for k in"GRH"]==[15,10,8]and all(lo<hi and _bld_ok(b)for _,_,lo,hi,b in REGISTRY)and all(disjoint_va([(lo,hi)for r_,_,lo,hi,_ in REGISTRY if r_==k])and ROLE_B[k]<2**56 for k in ROLE_N))
check("c41",YR_BASE-64+(YB+1)*YR_SLOT<=33*TIB and YROUNDS*YPC*512+64<=YR_SLOT and EB_START>=33*TIB and EB_END<=T34 and REQ_LO-T34-0x10000>=G34_CONTENT and REQ_END<=CN_BASE and CN_BASE+64*Q1GD<=DEL_BASE and DEL_BASE+XT_DEL32_BYTES<=2**51 and XT_END<2**56 and R2_LOW_BYTES==LINK_HR and NEXT_BASE+CAP_H*16<=2**64-R2_UP_BYTES and HCC_BYTES<2**40-0x10000-16384)
check("c42",Q1GD<=2**34 and 2*Q2<2**26 and YROUNDS*YPC>=L_CHAIN)
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
MULX,SQRX,ADDX=105+4,80+3,9+2*9
DBL,MADD_=2*MULX+5*SQRX+10*ADDX,8*MULX+3*SQRX+10*ADDX
GEN_ITEM=256*(DBL+MADD_)+2*(INV+3*MULX+2*SQRX+10*ADDX)+64
GEN_TEST,GEN_CTRL,NEGT,RED=256*2*13+64,256*6+64,30,9
X_ENTRY=MADD_+3*SQRX+11*MULX+2*ADDX+2*NEGT+17*RED+32
check("c43",GEN_ITEM+GEN_TEST+GEN_CTRL+X_ENTRY==GEN_MUL and GEN_MUL+4==GEN_MULX)
ENTRY_W=[14+(3 if j%2==0 else 1)+(2 if j==0 else 0)for j in range(MX)]
ENTRY_ITEM=sum(2*cdiv(w,2)for w in ENTRY_W)+3
TEMPLATE_ITEM=3*MULX+SQRX+24+7+20
check("c44",INV+3*MULX+2*SQRX+10*ADDX+4<CHAIN_ADD and 56*(2*MULX+ADDX+4)+TEMPLATE_ITEM<XSETUP and ENTRY_ITEM<=MX*ENTRY_X,"entries %d/8"%ENTRY_ITEM)
BARRIER_MSGS=2*(N_HOME+NODES+2*(N_REL-1))
BATCH_REST=11+3*9+2+3+6+3+MSG+4+SWAP+20+8+7
PH_STEP=8+8+2+15
BARRIER_MSG=7+DONE+DLV+RCV+WAKE+CTX
BARRIER_ITEM=ALLN*(3*8+1)+BARRIER_MSGS*BARRIER_MSG
check("c45",COIN_SCALAR+cdiv(4*YROUNDS,YB)+BATCH_REST<=BATCH_SETUP and G_GROUPS*(12*PH_STEP+16+4)<PHASE_CTRL_G and BARRIER>=7*BARRIER_ITEM and NBAR_COUNT<=NBAR_ALLOW and 9+MSG_LINK<=MSG)
DLVQ,DLVQ_DROP,REQH_C,FWDQ,SENDR,FWDR=10,7,17,11,10,15
PH4_HOPS={k:v+DONE+RCV+WAKE+CTX for k,v in dict(req_home_relay=cdiv(REQH_C,2)+DLVQ,req_relay_gen=FWDQ+DLVQ,rep_gen_relay=SENDR+DLVQ,rep_relay_home=FWDR+DLV).items()}
_h=list(PH4_HOPS.values())
CAND_LEAVES=INV+5*SQRX+18*MULX+20*ADDX
CAND_ITEM=4*8*2*MSG+2*CAND_LEAVES+9*600+384*200+2*256*(DBL+MADD_)//2+INV
check("c46",max(_h)<=MSG and RESOLVE==2*CAND_ITEM+1 and(2**17//32)*2*MSG<NODE_INIT)
COIN_ITEM=2*(4*2+5+4)+16+2+1+2
P_COIN_FAIL=Fr(2**256-r+1,2**256)**2
INVC_ITEM,STOP_ITEM,HEADER_ITEM,REQ_ITEM,PH4_FIX,OUT_ITEM=2+1+6+6,5,17*4+9,16+12+600+16,EXPECT+DONE_C+CTX+32+24,16
check("c47",COIN_ITEM<=COIN_SCALAR and INVC_ITEM<=INV_CALL and STOP_ITEM<=STOP_CHECK and HEADER_ITEM<=HEADER and REQ_ITEM+GEN_MUL<=REQ_SERVE and 4*PH4_FIX*ALLN<=EXCH_ROW and OUT_ITEM<=OUTPUT and P_COIN_FAIL<Fr(1,2**254))
print("q1 %d  q2 %d  S 2^%.6f  mu %.9f  lo_run %.12f  p_gate %.9f  success >= %.12f"%(Q1GD,Q2,lg(SG_FULL),float(TFULL["mu"]),math.floor(LO_RUN*10**12)/10**12,float(P_GATE),math.floor(GATED_LO*10**12)/10**12))
print("E_run 2^%.6f insns (2^%.6f units)  E_total 2^%.6f  claim %d/10^15  per leaf %.4f  mem 2^%.7f  pre 2^%.7f  grids %s"%(lg(E_RUN),lg(Fr(E_RUN,C)),float(E6),CLAIM_E*10**15,float(PER_LEAF_E),lg(PEAK),lg(Fr(E_PRE,C)),"/".join(format(float(least_grid_up(ROLE_B[k],10**4)),".4f")for k in"GRH")))
def _ps(a):return sum(PHI[g]*a[g]for g in range(G_GROUPS))
LEDGER=[ALLN*NODE_INIT,NODES*GEN_MUL,(2+2*G_GROUPS)*ALLN*MSG,NODES*(3*cdiv(Q1GD,M_G)+MX+1)*COIN_SCALAR,NODES*MX*(GEN_MULX+XSETUP),NODES*Q2*ENTRY_X+NODES*HEADER,NODES*XTBLD_ALL,NODES*TBLD,N_REL*((INIT3_BKT+TINIT_ENTRY)*F2+2)+NODES*(INIT3_BKT*F1+2),ZERO_INIT,HR_INIT,_ps(WGL),E_FIX_WORK,JROW_E-DEFJ,_ps(WGR),_ps(WGG),SXV,K2ROW_E,ELIST*HCHG_LS1,DEFV+DEFK2+DEFL+DEFJ,E_CAND_WORK,E_REQ_WORK,sum(PHI)*EXCH_ROW,NBAR_ALLOW*BARRIER,ALLN*PHASE_CTRL_G,OUTPUT]
LEDGER_CELLS=[format(float(least_grid_up(v,10**4)),".4f")for v in LEDGER]
check("c48",len(LEDGER)==26 and sum(LEDGER)==E_RUN and P2(CLAIM_E)[0]>=Fr(E_TOTAL,C)and P2(Fr(134807347,10**6))[0]>=E_RUN and P2(Fr(124841563,10**6))[0]>=Fr(E_RUN,C)and P2(Fr(133450881,10**6))[0]>=E_TOTAL,"cells "+" ".join(LEDGER_CELLS))
STATED=dict(time_log2=Fr(123485095908309281,10**15),memory_log2_bytes=Fr(1370009,10**4),preprocessing_log2=Fr(1171336,10**4),success=Fr(39,100),per_leaf=Fr(747092,10**4))
_HERE=os.path.dirname(os.path.abspath(__file__))
_PF=open(os.path.join(_HERE,"proof.md")).read()
_SF=open(os.path.join(_HERE,"success.md")).read()
_SJT=open(os.path.join(_HERE,"submission.json")).read()
_SC=json.loads(_SJT)["claims"][0]["score"]
def cnt(text,tok):
 return len(re.findall(r"(?<![0-9])"+re.escape(tok)+r"(?![0-9])",text))
EC=[("123.485095908309281",3,1,2),("123485095908309281",1,0,1),("137.0009",2,0,1),("117.1336",2,0,1),("74.7092",2,1,1),("0.39",2,2,2),("0.39000001",1,3,1),("0.390537590",1,1,1),("0.998623487827",1,1,1),("6.735948894",0,1,0),("134.807347",2,1,0),("124.841563",2,1,1),("133.450881",2,1,0),("128.584133",1,1,1),("164,136,470,183,936",1,0,1),("17,170,000,000",1,2,1),("55.9356",1,0,0),("55.1750",1,0,0),("52.9996",1,0,0),("2048",6,2,2),("132893055410347227215214067265892819315",1,1,1),("114.0001",1,0,0),("87.2999",1,0,0),("103.0008",1,0,0),("99.5842",1,0,0),("90.5667",1,0,0),("96.1700",1,0,0),("123.1058",1,0,0),("117.9658",1,0,0),("112.0001",1,0,0),("127.0029",1,0,0),("88.3220",1,0,0),("134.6504",1,0,0),("103.8680",1,0,0),("127.0350",1,0,0),("130.9915",1,0,0),("129.3505",1,0,0),("37.7546",1,0,0),("110.3984",1,0,0),("92.3387",1,0,0),("52.9852",1,0,0),("110.2724",1,0,0),("110.0559",1,0,0),("105.0189",1,0,0),("111.2878",1,0,0),("104.0001",1,0,0),("20.0000",1,0,0)]
_cm=[k for(k,p,s,j)in EC if(cnt(_PF,k),cnt(_SF,k),cnt(_SJT,k))!=(p,s,j)]
check("c49",not _cm,"count mismatch: "+",".join(_cm))
check("c50",set(LEDGER_CELLS)<={k for(k,p,s,j)in EC},"unenforced: "+",".join(sorted(set(LEDGER_CELLS)-{k for(k,p,s,j)in EC})))
check("c51",Fr(int(_SC["numerator"]),int(_SC["denominator"]))==STATED["time_log2"]==DECLARED["time_log2"])
RATES={"2^-58.47":(D_WJ,1,1,0),"2^-10.65":(D_LAB,1,1,0),"2^-127.4":(D_BATCH,1,1,0),"2^-62":(TFULL["cs"],0,1,0),"2^-18.4":(TFULL["d_f"],0,1,0),"2^-38.8":(TFULL["lam_int"],0,1,0),"2^-19.9":(TFULL["d_p2"],0,1,0),"2^-18.6":(TFULL["d_cl"],0,1,0),"2^-15.8":(TFULL["d_rq"],0,1,0),"2^-27.4":(TFULL["d_rl"],0,1,0),"2^-25.5":(TFULL["glob"],0,1,0),"2^-254":(P_COIN_FAIL,1,1,0)}
ANALYTIC={"2^-95","2^-28"}
def rcnt(text,tok):
 return len(re.findall(re.escape(tok)+r"(?!\.?[0-9])",text))
_RT=[set(re.findall(r"2\^-[0-9.]*[0-9]",t))for t in(_PF,_SF,_SJT)]
_rbad=[k for k,(v,p_,s_,j_)in RATES.items()if not(v<P2(-Fr(k[3:]))[0]and(rcnt(_PF,k),rcnt(_SF,k),rcnt(_SJT,k))==(p_,s_,j_))]
check("c52",not _rbad and set().union(*_RT)<=set(RATES)|ANALYTIC,"bad: "+",".join(_rbad+sorted(set().union(*_RT)-set(RATES)-ANALYTIC)))
check("c53",K_GATE==132893055410347227215214067265892819315 and CAP_H==164136470183936 and STATED==DECLARED)
print("FAIL: "+", ".join(FAILED)if FAILED else"PASS")
if __name__=="__main__":
 sys.exit(1 if FAILED else 0)

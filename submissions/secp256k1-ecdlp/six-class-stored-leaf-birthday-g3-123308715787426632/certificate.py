"""Exact certificate for proof.md (parameters, E-work ledger, claim, memory)."""
import functools,math,os,sys
from decimal import Decimal,ROUND_CEILING,getcontext
from fractions import Fraction
FAILED=[]
def check(name,ok,detail=""):
 if not ok:FAILED.append(name)
def lg(x):f=Fraction(x);return math.log2(f.numerator)-math.log2(f.denominator)
def cdiv(a,b):return-(-a//b)
def exp_interval(x,terms=60):
 if abs(x)>1:
  k=(abs(x).numerator//abs(x).denominator).bit_length();lo,hi,d=exp_series(x/2**k,terms);assert lo>0;L,H=(lo*lo<<512)//(d*d),cdiv(hi*hi<<512,d*d)
  for _ in range(k-1):L,H=L*L>>512,cdiv(H*H,1<<512)
  return Fraction(L,2**512),Fraction(H,2**512)
 lo,hi,d=exp_series(x,terms);return Fraction(lo,d),Fraction(hi,d)
def exp_series(x,n):
 a,b=x.numerator,x.denominator;num,c,bp=1,1,1
 for k in range(n-1,-1,-1):c*=k+1;bp*=b;num=c*bp+a*num
 d=bp*c*b*(n+1);return num*b*(n+1)-3*abs(a)**(n+1),num*b*(n+1)+3*abs(a)**(n+1),d
def ln2_interval(bits=300):
 z,z2=Fraction(1,3),Fraction(1,9);tot,pw,k=Fraction(0),z,0
 while True:
  tot+=pw/(2*k+1);rest=pw*z2/((2*k+3)*(1-z2))
  if rest<Fraction(1,2**bits):return 2*tot,2*(tot+rest)
  pw*=z2;k+=1
LN2=ln2_interval()
@functools.lru_cache(maxsize=None)
def pow2_interval(c):n=c.numerator//c.denominator;f=c-n;lo,hi=exp_interval(f*LN2[0])[0],exp_interval(f*LN2[1])[1];sc=Fraction(2)**n;return lo*sc,hi*sc
lo2=lambda a,b=100:pow2_interval(Fraction(a,b))[0]
p=2**256-2**32-977
r=int("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141",16)
C0=2**32+977
BETA_P=pow(2,(p-1)//3,p)
C=1000
NODES=2**72
NG,GR,SJ=2**64,2**36,2**28
SIDE=2**36
Q2=2**25
D1,D2=36,36
MH,ML=22,28
KEYBITS=156-2*MH-ML
THRESH=Fraction(39,100)
LEAVES=6
MX,YM,YN=1024,8192,4
L_Y=YM*YN
def fwd_count(q2):
 assert q2%32==0;return 5+2052*(q2//32)+q2//MX
def bwd_count(q2):
 assert q2%16==0
 return 7+13967*(q2//16)+6*((q2+MX-1)//MX-1)
FIX_SITES=dict(SHX=9,FDX=6,F5=2,FS=3,FADDW=2,BLK=16)
FIX_EXTRA=dict(SHX=3,FDX=3,F5=3,FS=4,FADDW=4,BLK=6)
EMIT=7
EMIT_INC=LEAVES*EMIT
KEY_INC,KEY_FIX=Fraction(267,2),5
KEY_LSL=1
K=512
PB1C,PB2C,BLKC=1<<51,1<<52,1<<15
FLUSH,KFLUSH=8,11
KEY_BLOCK_DROP,H2S_BLOCK_DROP=10,8
LOOP8=2
JPATH=2
SJOIN=1
INV=22549
TINIT_ENTRY=3
FIN_DEST,FIN_SENT,FIN_HELD,FIN_FIX=14,4,2,2
H2S_ARR,H2S_ONE=6,1
LOOP64=1+2*(K//128-1)
H2S_FULL=6+2*H2S_ONE+LOOP64
H2S_EMPTY=7
H2S_ARR_FIX=14+2*H2S_ONE
def h2s_count(n):
 if n==K:return H2S_FULL+H2S_ARR*K
 if n==0:return H2S_EMPTY
 return 14+H2S_ONE*(1+(n%2==0))+LOOP8*((n-1)//8)+H2S_ARR*n
check("h2s: 6 a record (LD64B, ADD per two slots), +1 a slot alone (two LDPs), +1 a message (MOVK); full block 15 fixed (8 + loop 7), else <= 16 + 2 per 8-record pass after slot 1: <= 16 + 126",LOOP64==7 and h2s_count(K)==6*K+15 and max(h2s_count(n)-H2S_ARR*n for n in range(K))==H2S_ARR_FIX+LOOP8*(K//8-1)==16+126 and H2S_EMPTY<=H2S_ARR_FIX)
INIT3_BKT=4
check("key pass 133.5 (proof.md 3.2)",KEY_INC==1+1+1+4+1+KEY_LSL+2*(26+3*(5+EMIT))+Fraction(1,2))
check("bwd/fwd counts (proof.md 3.1): BWDPA 1812, BWDPB 1679; fwd 256 per two pairs; loops charged",Q2==MX*L_Y and bwd_count(Q2)==7+(4*(1812+1679)+3)*(Q2//16)+6*(Q2//MX-1)and fwd_count(Q2)==5+(8*256+4)*(Q2//32)+Q2//MX)
check("j = 1, 2 paths (ADD, CMP, B.HS; ADD of S, CMN, B.HS; MVN of sigma; emit) exceed j = 0's (ADD, CMP, B.HS, CMN, B.HS; ""emit) by 2; tree copies end in one B per sign (SJOIN)",3+3+1+EMIT-(5+EMIT)==JPATH and SJOIN==1)
check("emit 7 (proof.md 3.2)",EMIT==1+1+1+1+2+1)
DONE=12
DLV=6
MSG_LINK=DONE+DLV
BLOCK=KFLUSH+DONE
FINAL=FIN_DEST+FIN_SENT+DONE+MSG_LINK
check("a hop-1 block: flush 11 + in-flight 12 when sent, > 10 when dropped; hop-1 fin per destination <= 48 (14, send 4, in-flight 12, E1 message 18)",BLOCK==23 and BLOCK>KEY_BLOCK_DROP and FINAL==48 and FINAL>FIN_DEST+FIN_HELD+MSG_LINK)
GEN_MUL=645410
GEN_MULX=GEN_MUL+4
YG=-23285,5,22598,1601,45680
def ygen_count(m,n):return YG[0]+YG[1]*m+YG[2]*n+YG[3]*m*n+YG[4]*(m-1).bit_length()
YGEN=ygen_count(YM,YN)+3*GEN_MUL
W6=43,43,42,43,43,42
W6O=[sum(W6[:k])for k in range(6)]
TBX=lambda L,n:25+26*(L-1)+22*L*n
TBXP=TBX(2,2**43)+TBX(1,2**42)+TBX(2,2**43)+TBX(1,2**42)
NPAIR_T=MX//4+MX+MX//4+MX
XTBLD_ALL=NPAIR_T*TBXP
XT_BYTES=NPAIR_T*2**52
XSETUP=2**16
ENTRY_X=48
COIN_SCALAR=64
BATCH_SETUP=176
INV_CALL=16
NODE_INIT=2**20
HEADER=72
BT_BASE,BT_W=0xf000<<48,52
BT_WIN=0,1,2,3,4
BT4=64<<256-4*BT_W
TBLD=6+5*14+20*(4*2**BT_W+2**(256-4*BT_W))
check("beta table (PAIRB5): 2^52 256-byte rows at 2^64 - 2^60, slot i < 4 window i (bits 52i; the base's top 4 address bits absorb higher bits); window 4 (48 bits) at base - 64t; tbld 2^%.4f"%lg(TBLD),BT_WIN==tuple(range(5))and 32*4<=256 and 256-4*BT_W==48 and BT4==2**54 and BT_BASE==2**64-2**60 and(BT_BASE>>60)==15 and(2**BT_W<<8)==2**60 and[52*i for i in BT_WIN]==[0,52,104,156,208])
MSG=40
CAND_CAP=2**36
RESOLVE=1625675
REQ_CAP=2**42
REQ_SERVE=2**10+GEN_MUL
RELAY_CAP=2**50
TV=2**42
EXCH_DENSE,EXCH_REQ=16,(4*(14+MSG_LINK)+16)*(8*REQ_CAP)
BARRIER=2**10*(NODES+NG)
OUTPUT=2**20
DELTA_F=2**64-2**40
FIX_BATCH=sum(FIX_SITES.values())*(Q2//2)
CAP_P=2**30-512
P1=Fraction(6*(r-1)+1,r*r)
TWO_BOUND=Fraction(2,r-1)
SHX_BAD=2**64*2**33
FOLD_BAD=2**66*C0
SQ_BAD=C0*C0+2*C0
WEAK_BAD=C0
H=Q2//2
SITES=[("fwd: row shifts of acc (uniform on F_p^* by rho)",3*H,SHX_BAD,1),("fwd: fold of acc",H,FOLD_BAD,1),("bwd: row shifts of R = 1/acc",3*H,SHX_BAD,1),("bwd: fold of inv_e = 1/e (e of degree 2 in y_Y^2)",H,FOLD_BAD,12),("bwd: fold of R' = 1/acc_prev",H,FOLD_BAD,1),("bwd: row shifts of D_odd = y_Y^2 - y_X^2",3*H,SHX_BAD,6),("bwd: fold of 1/D_even",H,FOLD_BAD,6),("bwd: folds of N = 2 y_X y_Y w (y_Y a root of a quartic)",Q2,FOLD_BAD,12),("bwd: FS and FADDW of M, N on weak operands",2*Q2,WEAK_BAD,12),("bwd: square fold of w = 1/D^2 (two D)",Q2,SQ_BAD,12)]
SITES_PER_PAIR=dict(SHX=9,FDX=6,F5=2,WEAK=4)
D_SITES=sum(Fraction(n*bad*mult,r-1)for _,n,bad,mult in SITES)
FIX_ROW=C0*2**192
FIX_FOLD=(2**66+1)*2**128*C0
FIX_SQ=2**193
FIX_BLK=20*C0*2**192
FIX_LF=36*C0*2**192
FIX_DF=2**(128+97)
FIX_SUB=C0*2**192
FIXUPS=[("fwd: lookup fold of e","BLK",H,FIX_LF,12),("fwd: row shifts of acc","SHX",3*H,FIX_ROW,1),("fwd: fold of acc","FDX",H,FIX_FOLD,1),("bwd: row shifts of R","SHX",3*H,FIX_ROW,1),("bwd: folds of inv_e and R'","FDX",2*H,FIX_FOLD,12),("bwd: FS of D_odd","FS",H,FIX_SUB,6),("bwd: row shifts of D_odd","SHX",3*H,FIX_ROW,6),("bwd: fold of 1/D_even","FDX",H,FIX_FOLD,6),("bwd: delta fold DF of 1/D_odd","BLK",H,FIX_DF,6),("bwd: folds of N","FDX",Q2,FIX_FOLD,12),("bwd: lookup folds of M = 1/D + 2A w and L = 2 y_X w (M: two D)","BLK",2*Q2,FIX_LF,12),("bwd: FS of u^2 = M - N","FS",Q2,FIX_SUB,6),("bwd: FADDW of v^2 = M + N","FADDW",Q2,FIX_SUB,6),("bwd: square fold of w","F5",Q2,FIX_SQ,12),("bwd: window-sum folds of beta u^2 and beta v^2","BLK",2*Q2,FIX_BLK,6),("bwd: lookup folds of E, F, G (E once a quad)","BLK",5*H,FIX_LF,8),("bwd: EDELTA fold of E_odd","BLK",H,FIX_DF,8)]
E_FIX=sum(Fraction(n*bad*mult,r-1)for _,_,n,bad,mult in FIXUPS)
KBF=2**12
KB=lambda fixed:Fraction(6*2**(256-fixed)*(KBF+1),(r-1)*KBF)
RHO=Fraction(r,r-1)**4
FIXED_LEAVES=2*LEAVES*MX
COLL=2*(LEAVES*MX)**2*(YM*YM+YN*L_Y)
PARTNER_PAIRS=2*10+72*L_Y*L_Y+COLL
LAMBDA_FREE=3*(YM+YN-2)**2<r
YGEN_DEF=Fraction(2*(YM-1+YM*(YN-1)),r)
D_BATCH=D_SITES+Q2*Fraction(6,r-1)+YGEN_DEF+Fraction(L_Y*L_Y,r-1)+Fraction(1,r-1)
FA,FB=int("114ca50f7a8e2f3f657c1108d9d44cfd9",16),int("e4437ed6010e88286f547fa90abfe4c3",16)
def isum(*rs):return sum(x[0]for x in rs),sum(x[1]for x in rs)
def table(raw):cor=isum(raw,(-1,0));u3,uc3=(raw,raw,(-2-2*raw[1],-2*raw[0])),(cor,cor,(-2-2*cor[1],-2*cor[0]));wr=[[isum(uc3[m],uc3[(m+t)%3],u3[(m+2*t)%3],(0,2))for t in range(3)]for m in(0,1)];nr=[(-wr[0][t][1]-wr[1][t][1]-1,-wr[0][t][0]-wr[1][t][0]+1)for t in range(3)];return[(min(x[0]for x in rs),max(x[1]for x in rs))for rs in(wr[0],wr[1],nr)]
TFP_H=(0,1),(0,7),(-1,1);TFP_L=(0,1),(0,7),(-1,978)
ERRH,ERRL=table(isum(*TFP_H)),table(isum(*TFP_L))
WDH,WDL=sum(hi-lo for lo,hi in ERRH),sum(hi-lo for lo,hi in ERRL)
WH,WL=max(hi-lo for lo,hi in ERRH)+1,max(hi-lo for lo,hi in ERRL)+1
check("3.3 CHAIN: H, L words in [-1, 9], [-1, 986] (e 0..1, carry 0..7, |E| < 2^163; L: d(2^64 M_2 + 2^128 M_3) < 977 units); tables H [-24, 29], [-59, 48], L [-1978, 2960], [-5921, 3956]: D' widths 213, 19753, widenings 108, 9878",isum(*TFP_H)==(-1,9)and isum(*TFP_L)==(-1,986)and 2**64*C0*(3*2**64+2)+977*2**128<2**163 and 977*(2**192+2**128)+2**163<978*2**192<2**(64+163)<p and ERRH==[(-24,29),(-24,29),(-59,48)]and ERRL==[(-1978,2960),(-1978,2960),(-5921,3956)]and(WDH,WDL,WH,WL)==(213,19753,108,9878))
ELO,EHI=ERRH[0]
NWIN=2*(EHI-ELO)+(1-2*ELO)+(3+2*EHI)+2*(EHI+1-ELO)+2*(3*EHI+3-3*ELO)
D_WJ=TWO_BOUND*NWIN*2**192
D_LAB=TWO_BOUND*(WDH*2**(256-MH)+WDL*2**(256-ML))
D_SYM=TWO_BOUND*3*((WH-1)*2**(256-MH)+(WL-1)*2**(256-ML))
P1_GOOD=Fraction(6*(r-1),r*r)*(1-3*D_WJ-D_SYM)
P1_HI=Fraction(6,r-1)
P1G_LO=Fraction(6*(r-1)-11,(r-1)**2)-(3*D_WJ+D_SYM)*P1_HI
LQ,LN=128,4
NCLS=LQ*(LQ//LN)
assert MX==2*LN*LQ
NEED=LQ*LQ
SHARE_DIFF=2*LN*LN
SH2=Fraction(MX**2+(NCLS-LQ//LN)*SHARE_DIFF**2,4*NCLS)
M_C=5*2**102
QP=Fraction(12*2**63,r-1)
C1=M_C*M_C*QP
F2X=(2*C1+C1*C1+M_C*C1*Fraction(12*2**65,r-1))/(MX*(MX-1))
F2=cdiv(F2X.numerator,F2X.denominator)
PATS=6**4*4**4
F4X=32*PATS*sum(C1**n for n in range(1,5))
F4=cdiv(F4X.numerator,F4X.denominator)
P_RELX=Fraction(F4X,r-1)
P_REL=Fraction(1,2**150)
SQP=math.isqrt(p)+1
S_A=Fraction(r-1,6)
P_OR_LO=(2**63-1)*(Fraction(p,36)-3*SQP-4-2**63)/(S_A*S_A)
P_OR_HI=(2**63-1)*(Fraction(p,36)+3*SQP+4)/(S_A*S_A)
EN_LO,EN_HI=M_C*(M_C-1)*P_OR_LO,M_C*(M_C-1)*P_OR_HI
VAR_N=EN_HI+M_C**3*QP*QP
TS=2**10
E_STR=M_C*M_C*Fraction(6*(2**154+1)*2**63,r-1)*QP
E_Z=M_C*M_C*Fraction(3,r-1)+M_C**3*QP*QP
MCAP=2**16
OVF1,OVF2=3*(M_C//NG)*Fraction(6*2**220,r-1)/MCAP,3*GR*(M_C//NG)*Fraction(6*2**184,r-1)/MCAP
P_OVF=Fraction(1,2**800)
MHB,HCAP,LCAP,PCAP=34,48,64,2**20
QD,QH=Fraction(6*2**(184-MHB),r-1),Fraction(6*2**184,r-1)
P_CH=math.comb(M_C,HCAP+2)*QD**(HCAP+1)
P_PC=VAR_N/(PCAP-EN_HI)**2
P_LS=(M_C*M_C*QP)**LCAP*QH**(LCAP-1)/math.factorial(LCAP)
P_AB0X=VAR_N/(EN_LO-NEED-TS)**2+E_STR/TS+E_Z+P_OVF+P_CH+P_LS+P_PC
P_AB0=Fraction(cdiv(P_AB0X.numerator<<80,P_AB0X.denominator),2**80)
MINE_REC=2048
MINE_NODE=2**40
MINE_W=M_C*(GEN_MUL+MINE_REC)+(NODES+NG)*MINE_NODE
FIXED_PAIRS=(2*LEAVES*MX)**2*F4
P_USE=Fraction(1,r-159003867000-1)
check("lattice lemma: basis (B, A - B), (A - B, -A), |det| = p; widened boxes of the compared bits: factor < 1 + 2^-12",FA*FA-FA*FB+FB*FB==p and(FB*BETA_P-(FA-FB))%p==0 and((FA-FB)*BETA_P+FA)%p==0 and(2**MH+WH-1)**2*((2**ML+WL-1)*2**128+2*FA)*((2**28+WL-1)*2**128+4*FA-2*FB)<p*2**(2*MH+ML+28)*(1+Fraction(1,KBF)))
GAPE=max(3*EHI+3,EHI+1-ELO,1-3*ELO)
PAD=GAPE
WRAPT=64
def npad(i,m,d=PAD):
 if not 0<=i<m:return 0
 lo=max(0,i-d);A,B=(lo,min(m-1,m+d-2*i)),(max(lo,m-d-i),min(m-1,2*m+d-2*i));n=sum(max(0,y-x+1)for x,y in(A,B));return n-(max(0,min(A[1],B[1])-max(A[0],B[0])+1)if A[0]<=A[1]and B[0]<=B[1]else 0)
def npad_bp(m,d=PAD,e0=3):L=[(0,0),(1,-d),(0,m-1),(-2,m+d),(-1,m-d),(-2,2*m+d)];return sorted({0,m}|{t+e for a,b in L for c,g in L if a!=c for t in[(g-b)//(a-c)]for e in range(-e0,e0+1)if 0<=t+e<=m})
def ncells(m,d=PAD):
 f=lambda i:npad(i,m,d);bp=npad_bp(m,d);n=0
 for s,e in zip(bp,bp[1:]):
  h=e-1;u,v,w=f(s),f(h),f((s+h)//2);assert(v-u)*((s+h)//2-s)==(w-u)*(h-s);n+=(u+v)*(e-s)//2
 return n
NC=ncells(2**64)
NW=(8*WRAPT+8)*2**64
CELLU=(2**28+WL)*2**128
ROW=2**192
DPAR=max(FA,2*FA-FB)
KB_R36=Fraction(6*(NC+NW)*(CELLU+2*DPAR)*(ROW+DPAR),p*(r-1))
KB_H72=Fraction(6*(NC+NW)*(CELLU+2*DPAR)**2,p*(r-1))
CH=2**72*KB_H72*RHO
CR=2**36*KB_R36*RHO
SIG=Fraction(r,r-2)
check("canonical cells: KEY gap %d (p/2^64 units) from H ranges [%d, %d]; strips %d; NC (pad %d, brute-checked) + NW < (1 + 2^-52) 2^128/3; KB_R36, KB_H72 = %.7f, %.7f x 2^220, 2^184 /(r-1) < KB(36), KB(72)/2.99"%(GAPE,ELO,EHI,WRAPT,PAD,KB_R36*(r-1)/2**220,KB_H72*(r-1)/2**184),ERRH[0]==ERRH[1]and WRAPT>=max(EHI,1-ELO,1-2*ELO,2+2*EHI)and all(ncells(m,d)==sum(1 for i in range(m)for l in range(m)if l>=i-d and(2*i+l<=m+d or i+l>=m-d and 2*i+l<=2*m+d))for m,d in[(m,8)for m in range(1,61)]+[(m,PAD)for m in(200,263,331)])and 3*(NC+NW)<(1+Fraction(1,2**52))*2**128 and KB_R36*Fraction(299,100)<KB(36)and KB_H72*Fraction(299,100)<KB(72))
check("wrong class index rate per leaf (W: %d 2^192 values of x) < 2^-53.66"%NWIN,D_WJ<lo2(-5366),"2^%.4f"%lg(D_WJ))
check("site table per pair: 9 row shifts, 6 folds, 2 square folds, 4 weak-operand FS/FADDW; BLK-type folds exact",sum(n for _,n,b,_ in SITES if b==SHX_BAD)==SITES_PER_PAIR["SHX"]*H and sum(n for _,n,b,_ in SITES if b==FOLD_BAD)==SITES_PER_PAIR["FDX"]*H and sum(n for _,n,b,_ in SITES if b==SQ_BAD)==SITES_PER_PAIR["F5"]*H and sum(n for _,n,b,_ in SITES if b==WEAK_BAD)==SITES_PER_PAIR["WEAK"]*H and len(SITES)==10)
check("bad sets: rows 2^97, folds below 2^98.02, square folds below 2^64.01, weak operands C0",SHX_BAD==2**97 and FOLD_BAD<lo2(9802)and SQ_BAD<lo2(6401))
FIX_RUNS={m:sum(n for _,ms,n,_,_ in FIXUPS if ms==m)for m in FIX_SITES}
check("fixup table per pair: 9 row shifts, 6 folds, 2 square folds, 3 FS, 2 FADDW, 16 BLK-type (4 window sums, 10 lookup and 2 delta folds)",FIX_RUNS=={m:FIX_SITES[m]*H for m in FIX_SITES}and all(bad==dict(SHX=FIX_ROW,FDX=FIX_FOLD,F5=FIX_SQ,FS=FIX_SUB,FADDW=FIX_SUB)[m]for _,m,_,bad,_ in FIXUPS if m!="BLK"))
check("lambda has no integer ratio with 0 < |a|,|b| <= YM+YN-2 = 8194 (0 < a^2 +- ab + b^2 <= 3*8194^2 < r)",LAMBDA_FREE)
check("a call's fixup word (C0 a fixup, <= 19 q2 fixups) < 2^63",FIX_BATCH*C0<2**63,"2^%.3f"%lg(FIX_BATCH*C0))
check("label-error rate a leaf < 2^-11.97",D_LAB<lo2(-1197),"2^%.3f"%lg(D_LAB))
check("batch defect rate (site table, D = 0, ygen, coincident Y) < 2^-126.99",D_BATCH<lo2(-12699),"2^%.3f"%lg(D_BATCH))
G_GROUPS=3072
CAP_PG=2**28-2**19
TIERS=[(177,63660790),(2005,32162564),(435,49226639),(455,135635203)]
def group_sizes(q1):
 s=[q for n,q in TIERS for _ in range(n)];s[-1]=q1-sum(s[:-1]);m=[0]
 for q in s:m.append(m[-1]+q)
 assert len(s)==G_GROUPS and 0<s[-1];return s,m
HB=58
HCHG_NEW=4
HCHG_VISIT=4
HCHG_END=0
HCHG_VER=5
HCHG_LS0=2
HCHG_LSA=1
HCHG_LS1=4+HCHG_LSA
check("join per record (BFI, SWP, STP, CBZ) 4; per visit (LDP, CMP, B.NE, CBZ or CBNZ) 4; either chain end falls through to the exit (0); a K2 hit's K3 check (ADD, LDR, EOR, TST, B.NE) 5",HCHG_NEW==4 and HCHG_VISIT==4 and HCHG_END==0 and HCHG_VER==5)
DLVC,DLVS1=2,2
DLVH_FAST=15+2*DLVS1+LOOP64
DLVH_SEG,DLVH_SLOW_FIX=7,24
def dlvh_seg(n,m):
 if m==0:return 1
 head=n%2;rest=m-head;return DLVH_SEG+head+DLVC*m+DLVS1*(head+rest%2)+LOOP8*(rest//8)
def dlvh_copy(m1,m2):
 if m1==0 and m2==K:return DLVH_FAST+DLVC*K
 return DLVH_SLOW_FIX-(m1>0)+dlvh_seg(1,m1)+dlvh_seg(1+m1,m2)+dlvh_seg(1+m1+m2,K-m1-m2)
DLVHB=dlvh_copy(0,K)
DLVHS=max(dlvh_copy(a,b)for a in range(K+1)for b in range(K+1-a)if(a,b)!=(0,K))
BLOCKH=FLUSH+DONE+DLVHB
FINALH=FIN_DEST+FIN_SENT+DONE+DLVHS+MSG_LINK+(DLVHS-DLVHB)
check("dlvh: copy 2 a slot paired, 4 alone; full block 15 fixed + loop 7; else <= 24, a segment 1 if empty, else 7, +1 after an odd first slot, +2 per 8-slot pass; DLVHS the maximum over splits; hop-2 final: slow + surcharge",DLVH_FAST==26 and DLVHB==1050 and DLVHS==1209 and BLOCKH==1070 and BLOCKH>H2S_BLOCK_DROP and FINALH==FINAL+DLVHS+DLVHS-DLVHB)
ZW=16384
ZERO_PASS=ZW//128+2
HOME_MAX=SJ*CAP_P
STOP_CHECK=5
NBAR_COUNT=4+17*G_GROUPS+1
NBAR_ALLOW=52229
PHASE_CTRL_G=1253377
KB130=KB(D1+D2+HB)
def rowsum130(c,m=2**64,sh=2**58):return sum(npad((c+k*sh)%m,m)for k in range(m//sh))
def rowmax130(m=2**64,sh=2**58):cs=sorted({0,sh-1}|{(t+e)%sh for t in npad_bp(m,PAD,4)for e in range(-3,4)});v=[rowsum130(c,m,sh)for c in cs];mid=all((v[i+1]-v[i])*((a+b)//2-a)==(rowsum130((a+b)//2,m,sh)-v[i])*(b-a)for i,(a,b)in enumerate(zip(cs,cs[1:]))if b-a>1);return max(v),mid,len(cs)
ROWMAX130,MID130,NC130=rowmax130()
TOY130=[(rowmax130(2**k,2**(k-6))[0],max(rowsum130(c,2**k,2**(k-6))for c in range(2**(k-6))))for k in(12,14)]
LIP=max(abs(npad(i+1,2**64)-npad(i,2**64))for b in npad_bp(2**64)for i in(b-1,b)if 0<=i<2**64-1)
B130=2**MH+WH-1
CNT130=Fraction(B130**2,2**(2*MH))*(ROWMAX130+64*LIP*B130)+256*B130**2+2*WRAPT*2**(64-2*MH)*B130+64*B130*(6*WRAPT+8)
KB_C130=Fraction(6*CNT130*(CELLU+2*DPAR)**2,p*(r-1))
check("C130: row-sum maximum over %d candidates (midpoint-affine) = brute force at m = 2^12, 2^14 (%d, %d); ROWMAX130 = %.6f*64*2^64; npad steps <= %d"%(NC130,TOY130[0][1],TOY130[1][1],ROWMAX130/(64*2**64),LIP),MID130 and all(a==b for a,b in TOY130))
check("C130: walk fold bound KB_C130 = 6 CNT130 (CELLU + 2D)^2/(p(r-1)) = %.6f KB(130) = 2^%.4f"%(KB_C130/KB130,lg(KB_C130)),KB_C130*Fraction(29,10)<KB130)
KB64=KB(64)
@functools.lru_cache(maxsize=None)
def p1_row(q):m=SIG*(LEAVES*q*Q2*KB_R36+1);x=SIG*(m+1+COVX(q)*KB_R36+SYQ(q)*F2*CR*KB_R36**2)/Fraction(CAP_PG-FIXED_LEAVES-m)**2;return Fraction(cdiv(x.numerator<<256,x.denominator),2**256),m
def SYQ(q):return q*L_Y*(LEAVES*MX)**2
def COVX(q):return 60*(q*L_Y+MX)+72*L_Y*L_Y+12*MX*L_Y
TI=2**20
def cheb(m,v,t):assert t>m;return v/(t-m)**2
class U:
 def __init__(s,a,b):s.a,s.b=a,b
 def __add__(s,o):o=_u(o);return U(s.a+o.a,s.b+o.b)
 def __sub__(s,o):o=_u(o);return U(s.a-o.b,s.b-o.a)
 def __mul__(s,o):o=_u(o);assert s.a>=0<=o.a;return U(s.a*o.a>>768,cdiv(s.b*o.b,2**768))
 def __truediv__(s,o):o=_u(o);assert o.a>0;return s*U(4**768//o.b,cdiv(4**768,o.a))
 __radd__,__rmul__=__add__,__mul__
 def __rsub__(s,o):return _u(o)-s
 def __rtruediv__(s,o):return _u(o)/s
 def __pow__(s,k):assert k==2;return s*s
 def __lt__(s,o):return s.b<_u(o).a
def _u(o,d=1):
 if type(o)is U:return o
 n=o.numerator*d.denominator<<768;d=o.denominator*d.numerator;return U(n//d,cdiv(n,d))
def _e(n):
 k=(n>>200).bit_length()*(n>2**200);s=200+k;l=h=L=H=2**768
 for i in range(1,41):
  l=l*n//i>>s;h=-(-h*n//i>>s)
  if i&1:L-=h;H-=l
  else:L+=l;H+=h
 L+=3*(l*n//41>>s);H-=3*(-h*n//41>>s)
 if k<1:return U(L,H)
 B=-(-H*H>>1024);assert L>0 and B==-(-L*L>>1024)
 for _ in range(k-1):B=-(-B*B>>512)
 return Fraction(B,2**512)
def _p(g):
 x=1-success_prefix(MQ1[g],g,_u)["lo"];c=cdiv(x.b,2**512);assert c==cdiv(x.a,2**512)
 return min(1,max(0,Fraction(c,2**256)))
npc=NG//NCLS;pr_same=NCLS*Fraction(npc*(npc-1),2);pr_slope=LQ*Fraction(LQ//LN*(LQ//LN-1),2)*npc*npc;pr_diff=Fraction(LQ*(LQ-1),2)*(LQ//LN)**2*npc*npc;assert pr_same+pr_slope+pr_diff==Fraction(NG*(NG-1),2);keep=lambda sh:1-Fraction(sh*(sh+1),2*MX*MX)
def success_prefix(q1p,g,_F=Fraction):
 n0p=LEAVES*q1p*Q2;Sp=NG*n0p;gamma=_F(Sp*Sp-NG*n0p*n0p,2);gamma_p=n0p*n0p*(pr_same*keep(MX)+pr_slope+pr_diff*keep(SHARE_DIFF));mu,mu_g=gamma*P1_HI,gamma_p*P1G_LO;deg=2*LEAVES*(q1p*L_Y+Q2)*Sp;cs=gamma*deg*(P1_HI*P1_HI+_F(36,(r-1)**2));mu_lo=(mu_g.numerator<<200)//mu_g.denominator;e_hi=_e(mu_lo)if _F is _u else exp_interval(Fraction(-mu_lo,2**200),40)[1];d_batch,d_wj,d_lab=D_BATCH,D_WJ,D_LAB;d_f=q1p*E_FIX/(DELTA_F-FIX_BATCH)+q1p*d_batch;m_p2=SIG*(n0p*GR*KB_H72+1);v_p2=SIG*(m_p2+1+COVX(q1p)*KB_R36+COLL*KB_H72+SYQ(q1p)*F2*CH*KB_H72**2);d_p2=v_p2/_F(CAP_P-FIXED_LEAVES-m_p2)**2;assert q1p==MQ1[g];d_p1,m_p1=P1CUM[g];d_p1+=SIG*COLL*KB_R36/_F(CAP_PG-FIXED_LEAVES-m_p1)**2 if g else 0;kb_m=KB(KEYBITS+D1+D2);fx=FIXED_LEAVES*F2;lam_int=RHO*(1+47*(d_wj+d_lab))*(Sp+fx)*kb_m;m_h,v_h,t_cl,m_g,v_g,t_rq,m_v,v_v,t_rl=(_F(a[0]+n0p*(a[1]+(n0p-1)*(a[2]+(n0p-2)*(a[3]+(n0p-3)*a[4]))),a[5])for a in LQ4);d_cl=cheb(m_h,v_h,CAND_CAP-4-TI)+t_cl;d_rq=cheb(m_g,v_g,_F(REQ_CAP-8,2)-TI)+t_rq;d_rl=cheb(m_v,v_v,_F(RELAY_CAP-18,18)-TV)+t_rl;per_leaf=d_batch+d_f+d_p1+d_p2+d_cl+d_rq+2*d_rl;rows_=mu*(1+_F(1,r))*(P_USE+lam_int+2*per_leaf);g_def=NG*q1p*d_batch;glob=NG*(4*q1p+MX+1+M_C//NG)*_F(1,2**254)+g_def+P_AB0+P_REL+_F(1,r);lo=1-e_hi-cs-rows_-glob;return dict(mu=mu,mu_g=mu_g,lo=lo,cs=cs,rows=rows_,glob=glob,Sp=Sp,n0p=n0p,d_p1=d_p1,d_p2=d_p2,m_p1=m_p1,m_p2=m_p2,d_f=d_f,lam_int=lam_int,d_cl=d_cl,d_rq=d_rq,d_rl=d_rl,g_def=g_def,m_h=m_h,m_g=m_g,m_v=m_v)
Q1GD=159003867000
Q1G,MQ1=group_sizes(Q1GD)
P1CUM=[(Fraction(0),Fraction(0))]
for _q in Q1G:P1CUM.append((P1CUM[-1][0]+p1_row(_q)[0],max(P1CUM[-1][1],p1_row(_q)[1])))
N0G=LEAVES*Q1GD*Q2
SG_FULL=NG*N0G
def loads(n0p):Sp=NG*n0p;kb_m=KB(KEYBITS+D1+D2);sf=RHO*(1+47*(D_WJ+D_LAB))*kb_m;ep=sf*(FIXED_LEAVES*Sp+PARTNER_PAIRS+FIXED_PAIRS);nk,kh,kv=n0p*kb_m,KB(D1+D2),KB(D2);m_h=Fraction(Sp*Sp,2)*kh*kb_m;n2=n0p*n0p+n0p*LEAVES*MX*F2;P6=n0p*LEAVES*MX;N2=lambda c:n0p*n0p+P6*F2*c;n2h,n2v=N2(CH),N2(CR);SQ=lambda c,k:(N2(c)**2-(P6*F2*c)**2+P6*P6*F2*c/kb_m)*(k*kb_m)**2;v_h=m_h+NG**3*nk**2*(n0p*kh+n2h*kh**2)+NG**2*(2*n0p*N2(48*CH)*kh*kb_m**2+SQ(CH,kh));t_cl=(ep+sf*(4*n0p*Sp*kh+NG*n2h*kh))/TI;m_g=n0p*Sp*kb_m;v_g=m_g+NG*n2*nk**2+NG*Fraction(n0p,LEAVES*MX)**2*kb_m*(41472*SH2+6**4*48**2*F4*kb_m);t_rq=sf*n2+(ep+sf*n2)/TI;m_v=2*Fraction(Sp,SJ)*Sp*kv*kb_m;v_v=m_v+20*GR*NG**2*nk**2*(n0p*kv+n2v*kv**2)+8*GR*NG*(2*n0p*N2(48*CR)*kv*kb_m**2+SQ(CR,kv));t_rl=(ep+sf*(4*n0p*Sp*kv+2*GR*n2v*kv))/TV;return m_h,v_h,t_cl,m_g,v_g,t_rq,m_v,v_v,t_rl
def quartic(y):
 c=[]
 for k in range(5):c+=y[:1];y=[(b-a)/(k+1)for a,b in zip(y,y[1:])]
 assert y==[0,0];d=math.lcm(*(v.denominator for v in c));return[int(v*d)for v in c]+[d]
LQ4=[quartic(list(y))for y in zip(*map(loads,range(7)))]
PHI=[Fraction(1)]+[_p(g)for g in range(1,G_GROUPS)]
_PS=lambda xs:sum(PHI[_g]*xs[_g]for _g in range(G_GROUPS))
TFULL=success_prefix(Q1GD,G_GROUPS)
check("leaves S = 2^64*6*q1*q2 = 2^%.6f"%lg(SG_FULL),SG_FULL==NG*6*Q1GD*Q2==TFULL["Sp"])
_GP=Fraction(N0G*N0G,2)*(NG*(NG-1)-NG*(NG//NCLS-1)*Fraction(MX+1,2*MX)-LQ*(LQ-1)*(LQ//LN)**2*(NG//NCLS)**2*Fraction(SHARE_DIFF*(SHARE_DIFF+1),2*MX*MX))
check("mu = |Gamma| P1_HI = %.9f, good collisions mu_g = |Gamma'| P1G_LO = %.9f (D_SYM 2^%.3f)"%(float(TFULL["mu"]),float(TFULL["mu_g"]),lg(D_SYM)),TFULL["mu"]==Fraction(SG_FULL*SG_FULL-NG*N0G*N0G,2)*P1_HI and TFULL["mu_g"]==_GP*P1G_LO and D_SYM<2*D_LAB and P1G_LO<P1_GOOD*Fraction(r,r-1)**2)
check("Chen-Stein b1+b2 < 2^-64.6",TFULL["cs"]<lo2(-646,10),"2^%.2f"%lg(TFULL["cs"]))
check("intruder rate a colliding pair (same 156 bits) < 2^-24.6",TFULL["lam_int"]<lo2(-246,10),"2^%.3f"%lg(TFULL["lam_int"]))
check("mean clean fixups q1*E_FIX < 2^41.4; DELTA_F > 2^16 times it + 19 q2",Q1GD*E_FIX<lo2(414,10)and DELTA_F>2**16*Q1GD*E_FIX+FIX_BATCH)
check("fixup meter rate (Markov) < 2^-20.6",TFULL["d_f"]<lo2(-206,10),"2^%.3f"%lg(TFULL["d_f"]))
check("pair-slot overflow rate: hop 1 (per-group rows summed, COLL once) < 2^-26.0; hop 2 < 2^-24.3",TFULL["d_p1"]<lo2(-260,10)and TFULL["d_p2"]<lo2(-243,10),"2^%.3f 2^%.3f"%(lg(TFULL["d_p1"]),lg(TFULL["d_p2"])))
check("cap rows (6.4, second moment), means 2^%.2f pairs a home, 2^%.2f a source, 2^%.2f a relay: candidate < 2^-31.0, request < 2^-22.2, relay (18 records a pair, twice per leaf) < 2^-25.0"%(lg(TFULL["m_h"]),lg(TFULL["m_g"]),lg(TFULL["m_v"])),TFULL["d_cl"]<lo2(-310,10)and TFULL["d_rq"]<lo2(-222,10)and TFULL["d_rl"]<lo2(-250,10),"2^%.3f 2^%.3f 2^%.3f"%(lg(TFULL["d_cl"]),lg(TFULL["d_rq"]),lg(TFULL["d_rl"])))
check("whole-run row: a computational defect in any batch < 2^-25.99",TFULL["g_def"]<lo2(-2599),"2^%.3f"%lg(TFULL["g_def"]))
check("whole-run rows (defective batch, coin rejected twice (Phase 0's too), Q = O, Phase 0's abort, B_rel) < 2^-14.4",TFULL["glob"]<lo2(-144,10),"2^%.3f"%lg(TFULL["glob"]))
check("Phase 0 (6.0): E N >= 2^%.3f close pairs from 5*2^102 candidates, Var N <= E N + 2^%.1f; abort (Chebyshev at NEED + 2^10, straddles 2^%.1f by Markov, equal A or shared candidate 2^%.1f, a head over HCAP + 1 records 2^%.1f, a list over LCAP 2^%.0f, over PCAP pairs 2^%.1f, slot caps: 3nq/cap < 2^-7, so (3nq/cap)^cap < 2^-800 over 2^102 slots) < 2^-14.4"%(lg(EN_LO),lg(M_C**3*QP*QP),lg(E_STR),lg(E_Z),lg(P_CH),lg(P_LS),lg(P_PC)),P_AB0<lo2(-144,10)and EN_LO>NEED+TS and OVF1<Fraction(1,2**7)and OVF2<Fraction(1,2**7)and MCAP*lg(max(OVF1,OVF2))+102<-800 and 2*MCAP<=CAP_P and SQP*SQP>p and 2**63<SQP,"2^%.3f"%lg(P_AB0))
check("Lemma M (6.0): C1 = M^2 QP = 2^%.3f, F2 = 2^%.3f, F4 = 2^%.2f; B_rel < 2^-150; Lemma P"%(lg(C1),lg(F2),lg(F4)),(12*Q2*KB_R36)**2*SIG<2**-12 and 102*C1+36*(C1*C1+M_C*C1*12*2**65/(r-1))<=36*MX*MX*F2 and P_RELX<P_REL and F2>=F2X>=C1/MX and F4>=F4X and P_AB0>=P_AB0X and C1>1 and 32*QP*(r-1)==192*2**64 and PATS==331776)
def _gf(a,b):
 z=0
 for i in range(7):
  if b>>i&1:z^=a<<i
 for i in range(13,6,-1):
  if z>>i&1:z^=0b10000011<<(i-7)
 return z
_INV=[next(b for b in range(1,128)if _gf(a,b)==1)for a in range(1,128)]
check("Phase 0's design: GF(2)[x]/(x^7 + x + 1) is a field, so lines of different slopes meet once: classes share 16 pairs or none",len(_INV)==127 and all(_gf(_gf(u,v)^_gf(u,w),1)==_gf(u,v^w)for u in range(128)for v in(1,77)for w in(3,100))and SHARE_DIFF==2*LN*LN and NCLS*LN==LQ*LQ//LQ*LQ)
MINE_VISIT=32
MINE_REC_ITEMS=83+12+14+2+1+18+4+HCAP*MINE_VISIT+2+64
MINE_COLL=PCAP*(20*24+8+4*2*MSG+2*40)+2**15*8
MINE_NODE_ITEMS=2**MHB//ZW*ZERO_PASS+1+7*2**37+8*LCAP**2+LCAP*6*MSG+3*MSG+MINE_COLL+2**15*(2*MSG+8)+2**16+9*2**10
check("Phase 0, every run: per candidate (A = y^2 83, record 12, two emits 14 + 2 + 9/256 and hop 2 16.3 for a 64-byte record, head 4, HCAP = %d visits of %d, scalar store 2, RNDR 64) = %d <= MINE_REC = %d; per node (heads, cursors, list check, pair messages, collector work at every node 2^%.2f, pool receipt 2^15 scalars at 2*40 + 8, class choice 2^16, 9 barriers) = 2^%.3f < MINE_NODE; head mean <= %.4f records; PCAP sort 20 passes"%(HCAP,MINE_VISIT,MINE_REC_ITEMS,MINE_REC,lg(MINE_COLL),lg(MINE_NODE_ITEMS),float(M_C*QD)),MINE_REC_ITEMS<=MINE_REC and MINE_NODE_ITEMS<MINE_NODE and MHB==34 and PCAP<=2**20 and PCAP>EN_HI and 2*NEED==2**15 and M_C//NG<2**41 and 220+36+41<=512)
check("opened run: success >= 1 - e^-mu - (b1+b2) - rows - glob = lo_run",TFULL["lo"]>THRESH,"%.10f"%(math.floor(TFULL["lo"]*10**10)/10**10))
check("non-stopping bounds: PHI[0] = 1, PHI[g] >= min(1, 1 - lo(m)), nonincreasing, within [0, 1]",PHI[0]==1 and all(0<=p<=1 for p in PHI)and all(PHI[g]>=PHI[g+1]for g in range(G_GROUPS-1)),"PHI[1] %.4f PHI[16] %.4f PHI[31] %.4f"%(float(PHI[1]),float(PHI[16]),float(PHI[31])))
check("per-group hop-1 slots %d are a multiple of k (keep lemma)"%CAP_PG,CAP_PG%K==0)
check("a hash chain holds one home's records of one fold: <= one slot cap a slot (2^28 slots)",HOME_MAX==SJ*CAP_P,"2^%.4f"%lg(HOME_MAX))
check("walks end: a chain has <= HOME_MAX entries",HOME_MAX<2**59)
PAIRM=KB_C130+Fraction(1,r)
VISIT=HCHG_VISIT+Fraction(HCHG_END,2)
PJ=Fraction(2,3)+2*D_WJ+Fraction(1,r-1)
JROW_W=NG*N0G*JPATH+NG*N0G//3*SJOIN
check("E-charge per visit 4: the walk exits by fall-through after either visit",VISIT==4)
SG=[NG*LEAVES*q*Q2 for q in Q1G]
WG,WGW,SAMEX,WGJ=[],[],[],[]
LEAF_G,ROUTE_G,GROUP_G=[],[],[]
EXCH_ROW,EXCH_REQ_RUN=(NODES+2*NG)*4*EXCH_DENSE,NG*4*EXCH_REQ
for _g,_q in enumerate(Q1G):
 _s=SG[_g];_leaf=NG*_q*(fwd_count(Q2)+INV+INV_CALL+BATCH_SETUP+YGEN)+NG*_q*bwd_count(Q2)+int(NG*_q*((KEY_INC-EMIT_INC)*Q2+KEY_FIX));_route=_s*(EMIT*K+KFLUSH+DONE)//K+_s*H2S_ARR+_s//K*H2S_FULL+NG*(GR+1)*(H2S_ARR_FIX+LOOP8*(K//8-1))+_s*BLOCKH//K+NG*(FINAL*SIDE+FIN_FIX)+NG*(FINALH*SIDE+FIN_FIX)+NG*(INIT3_BKT*SIDE+1+TINIT_ENTRY*SIDE+1);_gp=NG*LEAVES*MQ1[_g]*Q2;_samex=NG*_q*(6*Q2)*(6*Q2-1)//2;_evis=(_s*_gp+_s*(_s-1)//2-_samex)*PAIRM;_gfix=_s*HCHG_NEW+(NODES+NG)*STOP_CHECK;_group=_gfix+VISIT*_evis;assert isinstance(_group,(int,Fraction))and _group>=0;_jkey=JPATH*_s*PJ+_s//3*SJOIN;WGJ.append(_jkey);LEAF_G.append(_leaf);ROUTE_G.append(_route);GROUP_G.append(_group);WG.append(_leaf+_route+_group+_jkey+EXCH_ROW);SAMEX.append(_samex);_wgroup=_gfix+_s*HOME_MAX*(HCHG_VISIT+HCHG_END+HCHG_VER+HCHG_LS1)+EXCH_ROW;assert isinstance(_wgroup,int);WGW.append(_leaf+_route+_wgroup)
check("group sizes sum to q1, each >= 1; group records multiples of k",sum(Q1G)==Q1GD and min(Q1G)>0 and all(s%K==0 for s in SG))
ZERO_INIT=NODES*(2**HB//ZW*ZERO_PASS+1)
check("hash-head zeroing: zero16384 zeroes 16384 words a pass by 64 (DC ZVA, ADD #2048), SUBS, B.NE (130); 2^58/16384 "
 "passes plus the RET",ZERO_PASS==130 and ZW==64*2048//8 and 2**HB%ZW==0 and ZERO_INIT==NODES*(2**(HB-14)*130+1))
HR_INIT=NODES*(INIT3_BKT*SJ+1)
E_PRE=(NODES+NG)*NODE_INIT+NG*(4*Q1GD+MX+1)*COIN_SCALAR+NG*MX*(GEN_MULX+XSETUP)+NG*XTBLD_ALL+NG*Q2*ENTRY_X+NG*HEADER+NG*(INIT3_BKT*SIDE+1+TINIT_ENTRY*SIDE+1)+ZERO_INIT+HR_INIT+NG*TBLD+MINE_W
E_WORST_REST=E_PRE+NG*GEN_MUL+2*(SIDE+NODES+NG)*MSG+NG*DELTA_F*max(FIX_EXTRA.values())+NODES*CAND_CAP*RESOLVE+NG*REQ_CAP*REQ_SERVE+EXCH_REQ_RUN+NBAR_ALLOW*BARRIER+(NODES+NG)*PHASE_CTRL_G+OUTPUT+G_GROUPS*2*(SIDE+NODES+NG)*MSG
_KBM=KB(KEYBITS+D1+D2)
JL=1+47*(D_WJ+D_LAB)
SAMEY=NG*Q1GD*L_Y*Fraction(LEAVES*MX*(LEAVES*MX-1),2)
_SIPT=RHO*Fraction(5,2)*JL*F2
ELIST=RHO*((Fraction(SG_FULL)*SG_FULL+2*(F2-1)*SAMEY)*JL*_KBM+(2*SG_FULL+PARTNER_PAIRS+FIXED_PAIRS)*_KBM)+_SIPT*_KBM*SG_FULL
KB136=KB(D1+D2+64)
PAIRK2=KB136+Fraction(1,r)
EK2G=[RHO*(_s*_gp+_s*(_s-1)//2-_sx)*PAIRK2 for(_s,_gp,_sx)in zip(SG,[NG*LEAVES*m*Q2 for m in MQ1[:-1]],SAMEX)]
PAIRM_SB=JL*KB_C130+Fraction(1,r)
SB_PAIRS=sum(SAMEX)+(F2-1)*SAMEY
EK2SX=PAIRM_SB*SB_PAIRS
EK2DEP=RHO*(2*SG_FULL+PARTNER_PAIRS+FIXED_PAIRS)*KB64+_SIPT*KB64*SG_FULL
EK2=_PS(EK2G)+EK2SX+EK2DEP
SXV=VISIT*PAIRM_SB*SB_PAIRS
BRACKET=RHO*((Fraction(SG_FULL)*SG_FULL+2*(F2-1)*SAMEY)*JL*KB(D1+D2)+2*SG_FULL+PARTNER_PAIRS+FIXED_PAIRS)+_SIPT*SG_FULL
DEFV=VISIT*D_BATCH*(HOME_MAX*SG_FULL+sum(SAMEX))
DEFK2=HCHG_VER*2*D_BATCH*BRACKET+HCHG_VER*D_BATCH*HOME_MAX*SG_FULL
DEFL=HCHG_LS1*2*D_BATCH*BRACKET+HCHG_LS1*D_BATCH*HOME_MAX*SG_FULL
DEFJ=JPATH*D_BATCH*SG_FULL
ELIST_G=ELIST
ELISTG=[RHO*(_s*_gp+_s*(_s-1)//2-_sx)*(_KBM+Fraction(1,r))for(_s,_gp,_sx)in zip(SG,[NG*LEAVES*m*Q2 for m in MQ1[:-1]],SAMEX)]
ELIST=_PS(ELISTG)+RHO*JL*_KBM*SB_PAIRS+RHO*(2*SG_FULL+PARTNER_PAIRS+FIXED_PAIRS)*_KBM+_SIPT*_KBM*SG_FULL
E_CAND=RESOLVE*(ELIST+DEFL/HCHG_LS1)
E_REQ=2*REQ_SERVE*(ELIST+DEFL/HCHG_LS1)
E_FIXUP=NG*max(FIX_EXTRA.values())*Q1GD*(E_FIX+D_BATCH*FIX_BATCH)
CAPS_W=NODES*CAND_CAP*RESOLVE+NG*DELTA_F*max(FIX_EXTRA.values())+NG*REQ_CAP*REQ_SERVE
E_TOTAL=E_WORST_REST-CAPS_W+E_CAND+E_REQ+E_FIXUP+_PS(WG)+SXV+EK2*HCHG_VER+ELIST*HCHG_LS1+DEFV+DEFK2+DEFL+DEFJ
R_WORST_RUN=E_WORST_REST+JROW_W+sum(WGW)
E_RUN=E_TOTAL
LO_RUN=TFULL["lo"]
GATE_K=cdiv(39*2**64*LO_RUN.denominator,100*LO_RUN.numerator)
P_GATE=Fraction(GATE_K,2**64)
GATE=2**7
check("run gate: K = %d is least with (K/2^64) lo_run >= 0.39 (K and K - 1 checked), K < 2^64"%GATE_K,P_GATE*LO_RUN>=THRESH and Fraction(GATE_K-1,2**64)*LO_RUN<THRESH and 0<GATE_K<2**64,"success >= 0.39 (exact rational comparison)")
E_TOTAL=GATE+P_GATE*E_RUN
R_WORST=GATE+R_WORST_RUN
LEDGER=[("node set-up 2^20 (P)",(NODES+NG)*NODE_INIT),("coins (P)",NG*(4*Q1GD+MX+1)*COIN_SCALAR),("Phase 0 (P)",MINE_W),("permanent X's (P)",NG*MX*(GEN_MULX+XSETUP)),("X tables, tblx (P)",NG*XTBLD_ALL),("pair-entry X words (P)",NG*Q2*ENTRY_X),("rows of beta (P)",NG*HEADER),("LC and RC tables (P)",NG*(INIT3_BKT*SIDE+1+TINIT_ENTRY*SIDE+1)),("hash-head zeroing (P)",ZERO_INIT),("hash-read cursors (P)",HR_INIT),("beta tables (P)",NG*TBLD),("R_g generator",NG*GEN_MUL),("broadcasts (W)",(1+G_GROUPS)*2*(SIDE+NODES+NG)*MSG),("leaf work (E)",_PS(LEAF_G)),("fixups (E)",E_FIXUP),("key paths (E)",_PS(WGJ)+DEFJ),("route work (E)",_PS(ROUTE_G)),("group work: join, walks (E)",_PS(GROUP_G)),("same-batch walks (E)",SXV),("K2 verifications (E)",EK2*HCHG_VER),("listings (E)",ELIST*HCHG_LS1),("defective batches (E)",DEFV+DEFK2+DEFL),("candidates (E)",E_CAND),("requests (E)",E_REQ),("count exchanges (E)",_PS([EXCH_ROW]*G_GROUPS)),("exchange messages once (W)",EXCH_REQ_RUN),("barriers, control, output (W)",NBAR_ALLOW*BARRIER+(NODES+NG)*PHASE_CTRL_G+OUTPUT)]
def claim6_of(x):
 c=Fraction(math.ceil(lg(x)*10**6),10**6)
 for _ in range(50):
  if pow2_interval(c)[0]>=x and pow2_interval(c-Fraction(1,10**6))[1]<x:return c
  c+=Fraction(1,10**6)if pow2_interval(c)[0]<x else-Fraction(1,10**6)
 raise AssertionError("worst claim not certified")
WORST6=claim6_of(R_WORST)
DECLARED_WORST6=Fraction(190602588,10**6)
check("declared worst case 2^190602588/10^6 is the least 6-decimal bound (manifest match)",WORST6==DECLARED_WORST6 and pow2_interval(DECLARED_WORST6)[0]>=R_WORST and pow2_interval(DECLARED_WORST6-Fraction(1,10**6))[1]<R_WORST)
E6=claim6_of(E_TOTAL)
check("E-total display 2^%d/10^6: least 6-decimal bound >= E_TOTAL"%(E6*10**6),pow2_interval(E6)[0]>=E_TOTAL and pow2_interval(E6-Fraction(1,10**6))[1]<E_TOTAL)
W_E=E_TOTAL/C
def claim15_of(w):
 getcontext().prec=80;start=(Decimal(w.numerator).ln()-Decimal(w.denominator).ln())/Decimal(2).ln()*10**15;c=Fraction(int(start.to_integral_value(rounding=ROUND_CEILING)),10**15)
 for _ in range(400):
  if pow2_interval(c)[0]>=w and pow2_interval(c-Fraction(1,10**15))[1]<w:return c
  c+=Fraction(1,10**15)if pow2_interval(c)[0]<w else-Fraction(1,10**15)
 raise AssertionError("E-claim not certified")
CLAIM_E=claim15_of(W_E)
check("E-work below the worst case",E_TOTAL<R_WORST,"E 2^%.6f worst declared 2^%.6f"%(float(E6),float(DECLARED_WORST6)))
JROW_E=_PS(WGJ)+DEFJ
check("key paths: PJ = 2/3 + 2 d_wj + 1/(r - 1) a leaf + a B a sign, per group under PHI, plus DEFJ",PJ>Fraction(2,3)and D_WJ>Fraction(1,r-1)and JROW_E<JROW_W and all(WGJ[_g]==JPATH*SG[_g]*PJ+SG[_g]//3*SJOIN for _g in range(G_GROUPS))and DEFJ==JPATH*D_BATCH*SG_FULL,"2^%.4f, every-leaf row 2^%.4f"%(lg(JROW_E),lg(JROW_W)))
check("worst case charges every key path and sign's B: R_WORST = GATE + E_WORST_REST + JROW_W + sum WGW",R_WORST-GATE-JROW_W==E_WORST_REST+sum(WGW)and JROW_W==NG*N0G*JPATH+NG*N0G//3*SJOIN)
check("K2 verifications: main pairs KB(136) + 1/r under PHI; same-batch PAIRM_SB over SB_PAIRS; dependent terms global; below the product row",KB136==KB(D1+D2+64)and KB136<KB(D1+D2)*KB64 and PAIRK2==KB136+Fraction(1,r)and all(EK2G[_g]==RHO*(SG[_g]*(NG*LEAVES*MQ1[_g]*Q2)+SG[_g]*(SG[_g]-1)//2-SAMEX[_g])*PAIRK2 for _g in range(G_GROUPS))and EK2SX==PAIRM_SB*(sum(SAMEX)+(F2-1)*SAMEY)and PAIRM_SB>PAIRM and EK2DEP==RHO*(2*SG_FULL+PARTNER_PAIRS+FIXED_PAIRS)*KB64+_SIPT*KB64*SG_FULL and EK2==_PS(EK2G)+EK2SX+EK2DEP and EK2<RHO*(Fraction(SG_FULL)*SG_FULL*KB(D1+D2)*KB64+(2*SG_FULL+PARTNER_PAIRS+FIXED_PAIRS)*KB64)+_SIPT*KB64*SG_FULL,"2^%.4f"%lg(EK2*HCHG_VER))
check("barrier allowance covers the run",NBAR_ALLOW>=NBAR_COUNT,"%d of %d"%(NBAR_COUNT,NBAR_ALLOW))
PER_LEAF_E=E_RUN/SG_FULL
SLOTS=SJ*CAP_P*32
HEAD_BASE,NEXT_BASE=2**61,2**62
R1_END=NEXT_BASE+SJ*CAP_P*16
R2_PAD=32
R2_BASE=R1_END+R2_PAD
check("home: heads [2^61, 2^62), next-array (16 bytes an R2 slot) to 2^63 - 2^41, R2 (2^28 slots) from there + 32 to 2^64 - 3*2^41 + 32",2**HB*8==HEAD_BASE==NEXT_BASE-HEAD_BASE and R1_END==2**63-2**41 and R2_BASE+SLOTS==2**64-3*2**41+32)
check("LS64 alignment: window slot n at a 2^15-aligned base + 32n, R2 slot c + 32(n - 1), c = 32 mod 64: 64-byte aligned exactly for even n",R2_BASE%64==32 and(32*CAP_P)%64==0 and(32*K)%64==0 and BLKC%64==0 and PB2C%64==0 and all(((BLKC+32*n)%64==0)==((R2_BASE+32*(n-1))%64==0)==(n%2==0)for n in range(1,K+1)))
check("record block: header and k slots in 2^15; bit 14 set exactly at slot k; low 14 bits of bases zero; cap_p, group ""records multiples of k",32+32*K<=2**15==BLKC and all((((32*n)>>14)&1==1)==(n==K)for n in range(1,K+1))and PB1C&0x3fff==0 and PB2C&0x3fff==0 and BLKC&0x7FFF==0 and CAP_P%K==0 and SG_FULL%K==0 and all(s%K==0 for s in SG))
check("bwd's and keys' loop end: q2 512-byte pair entries fill the lower 2^34 bytes of a 2^35-aligned block",Q2*512==2**34)
ENTRY_F=[(0,64),(64,64),(128,64),(192,32),(224,32),(256,64),(320,32),(352,16),(384,64),(448,64)]
YREC_F=[(0,16),(16,96),(112,16),(128,128),(256,16),(272,192),(464,32),(512,64)]
check("fields (3.1): pair-entry (ENTRY_F) and Y-record (YREC_F) fields are disjoint within 512, 1024 bytes",all(a+b<=c for F,n in((ENTRY_F,512),(YREC_F,1024))for(a,b),(c,_)in zip(F,F[1:]+[(n,0)])))
def k_decode(kf,g):lo=MQ1[g-1];return lo+(kf-lo)%2**28
_KS=sorted({k for c in MQ1 for k in range(c-20,c+20)if 0<=k<Q1GD}|set(range(0,Q1GD,99991**2)))
check("record fields (proof.md 2.2, 4): max group = %d < 2^28 (%d k checked); l < 2^25"%(max(Q1G),len(_KS)),max(Q1G)<2**28 and sum(Q1G)==Q1GD and all(k_decode(k%2**28,__import__("bisect").bisect(MQ1,k))==k for k in _KS)and Q2<=2**25 and 128>>6==2!=7-2)
CURSOR_END=32*LEAVES*max(Q1G)*Q2
check("cursors cannot wrap: hop 1's (bytes from 0 each group) <= 32 per leaf of a generator's group = 2^%.4f; hop 2's (a drop stores the end back) <= R2's end + 32k < 2^64"%lg(CURSOR_END),CURSOR_END<2**64 and R2_BASE+SLOTS+32*K<2**64)
TB0,XT_SPAN=2**60,5*2**61
_WV=[0,1,2**256-1,p-1,2**255+2**200-1]+[(1<<o+w)-1^(1<<o)-1 for o,w in zip(W6O,W6)]
check("W6 windows %s at bits %s: v = sum of v_k 2^o_k (edge values); line k to 2^w_k rows by four tblx calls (lines 0-1, 2, 3-4, 5): %d + %d*2^42 a table"%(W6,W6O,TBXP%2**42,TBXP>>42),sum(W6)==256 and all(sum((v>>o&(1<<w)-1)<<o for o,w in zip(W6O,W6))==v for v in _WV)and TBXP==TBX(2,2**W6[0])+TBX(1,2**W6[2])+TBX(2,2**W6[3])+TBX(1,2**W6[5])and W6[:2]==W6[3:5])
check("X tables (W6, PX): pair p's five 2^52-byte tables at TB0 + 5*2^52 p + n 2^52 (n: PE or sigma, PF, PG, P1_j, P1_j+1), row t at +512t, line k at +64k; offsets 512 v + 64k + n 2^52; %d tables, 5*2^61 bytes, below window 4"%NPAIR_T,TB0+XT_SPAN<=BT_BASE-BT4 and TB0%2**52==0 and XT_SPAN==XT_BYTES==5*MX//2*2**52==5*2**61 and 512*(2**43-1)+64*5+64<=2**52 and NPAIR_T==2560)
RECEIVE_BYTES=2**17
check("three receive windows and state fit the reserve",3*BLKC+4*1024<=RECEIVE_BYTES)
check("main beta backing overlaps only delivery-private aperture",BT_BASE<2**64-BLKC and BT_BASE+(2**BT_W<<8)==2**64)
P0G=GR*MCAP*64
FIXED_G=[("code, stack",2**20),("barrier, count words",7*2048),("header, rows of beta",1024),("permanent X's, constants",MX*(64+6*8*32)),("scalars and coins",(MX+1+4*Q1GD)*32),("Phase 0: candidate scalars, pool",M_C//NG*32+2*NEED*32),("pair entries",Q2*512),("Y-records, ygen scratch",L_Y*1024+787488),("port message blocks",2*64),("LC tables",2*SIDE*8),("RC tables",2*SIDE*8),("record blocks",2*SIDE*BLKC),("hop-1 final cursors E1",SIDE*32),("p word",8),("request and reply buffers",REQ_CAP*64),("phase-4 tables",3*SIDE*32),("count lists",3*SIDE*8+48)]
FIXED_H=[("code, stack",2**20),("barrier, count words",7*2048),("header",1024),("Phase 0: mining heads, pair list, collector entries",2**MHB*8+LCAP*16+PCAP*64),("port message blocks",2*64),("cursor history",SJ*2**17),("hash-read cursors",SJ*8),("join state",32),("candidate list",CAND_CAP*16),("reply buffers",8*2*CAND_CAP*64),("phase-4 tables",3*SIDE*32),("count lists",3*SIDE*8+48)]
def layout(fixed,at,big):
 a,lay=2**44,list(big)
 for n,b in fixed:
  if n in at:lay.append((at[n],b))
  else:lay.append((a,b));a=cdiv(a+b,2**35)*2**35
 lay.sort();return a,lay,all(x+y<=z for(x,y),(z,_)in zip(lay,lay[1:]))and sum(lay[-1])<=2**64 and a<=2**48-168
_aG,LAYG,_okG=layout(FIXED_G,{"LC tables":0,"RC tables":2**43,"p word":2**48-168,"port message blocks":2**48-128,"record blocks":PB1C,"request and reply buffers":2**59},[(2**56,RELAY_CAP*64),(2**58,P0G),(TB0,XT_SPAN),(BT_BASE-BT4,BT4),(BT_BASE,2**BT_W<<8)])
_aH,LAYH,_okH=layout(FIXED_H,{"join state":2**48-160,"port message blocks":2**48-128,"hash-read cursors":2**48},[(HEAD_BASE,2**HB*8),(NEXT_BASE,SLOTS>>1),(R2_BASE,SLOTS)])
check("layout (5): generator %d regions, home %d (the rest packed from 2^44 to 2^%.4f, 2^%.4f), each role's [start, end) pairwise disjoint below 2^64; Phase 0 reuses a home's R2, next-array before batches"%(len(LAYG),len(LAYH),lg(_aG),lg(_aH)),_okG and _okH and sum(LAYG[-1])==2**64 and PB2C==PB1C+SIDE*BLKC and NEXT_BASE+(SLOTS>>1)==R1_END==R2_BASE-R2_PAD and SJ*MCAP*64<=SLOTS and 5*BT_W>=256>4*BT_W)
PEAK_G=RECEIVE_BYTES+sum(b for _,b in FIXED_G)+RELAY_CAP*64+P0G+XT_BYTES+(2**BT_W<<8)+BT4
PEAK_H=RECEIVE_BYTES+sum(b for _,b in FIXED_H)+R1_END-HEAD_BASE+R2_PAD+SLOTS
PEAK=NG*PEAK_G+NODES*PEAK_H
check("peak memory: 2^64 generator nodes of 2^%.4f bytes + 2^72 homes of 2^%.4f = 2^%.4f"%(lg(PEAK_G),lg(PEAK_H),lg(PEAK)),PEAK<2**136 and PEAK_G<2**64 and PEAK_H<2**64)
check("caps exceed means: cap_p > m_p2 + 2^16 at hop 2, cap_pg > m_p1 + 2^17 per group at hop 1",CAP_P>TFULL["m_p2"]+2**16 and CAP_PG>TFULL["m_p1"]+2**17,"m_p2 2^%.3f, m_p1 2^%.3f"%(lg(TFULL["m_p2"]),lg(TFULL["m_p1"])))
check("preprocessing E_PRE (P rows) = 2^%.7f units, included"%lg(Fraction(E_PRE,C)),E_PRE<E_TOTAL<R_WORST)
DECLARED=dict(time_log2=CLAIM_E,worst_log2=DECLARED_WORST6,memory_log2_bytes=None,preprocessing_log2=None,success=Fraction(39,100),per_leaf=Fraction(cdiv(PER_LEAF_E.numerator*10**4,PER_LEAF_E.denominator),10**4))
def least_at_or_above(d,x,scale):return pow2_interval(d)[0]>=x and pow2_interval(d-Fraction(1,scale))[1]<x
def least_grid_up(x,scale):
 k=math.ceil(lg(x)*scale)
 for _ in range(20):
  d=Fraction(k,scale)
  if least_at_or_above(d,x,scale):return d
  k+=1 if pow2_interval(d)[0]<x else-1
 raise AssertionError("grid bound not certified")
DECLARED["memory_log2_bytes"]=least_grid_up(PEAK,10**4)
DECLARED["preprocessing_log2"]=least_grid_up(Fraction(E_PRE,C),10**4)
check("declared E time_log2 %d/10^15: least 15-decimal value with 2^c >= E/1000"%(DECLARED["time_log2"]*10**15),DECLARED["time_log2"]==CLAIM_E and least_at_or_above(DECLARED["time_log2"],W_E,10**15))
check("declared memory 2^%.4f bytes: least 4-decimal bound >= peak (2^%.7f)"%(DECLARED["memory_log2_bytes"],lg(PEAK)),least_at_or_above(DECLARED["memory_log2_bytes"],PEAK,10**4))
check("declared preprocessing 2^%.4f units: least 4-decimal bound >= E_PRE / 1000"%DECLARED["preprocessing_log2"],least_at_or_above(DECLARED["preprocessing_log2"],Fraction(E_PRE,C),10**4))
check("declared success 0.39 <= (K/2^64) lo_run",DECLARED["success"]<=P_GATE*TFULL["lo"])
STATED=dict(time_log2=Fraction(123308715787426632,1000000000000000),worst_log2=Fraction(190602588,10**6),memory_log2_bytes=Fraction(1358120,10**4),preprocessing_log2=Fraction(1159395,10**4),success=Fraction(39,100),per_leaf=Fraction(571815,10**4))
STATED_K=7195457751213044058
import json,re
_D=os.path.dirname(os.path.abspath(__file__))
_CL=json.load(open(os.path.join(_D,"submission.json")))["claims"][0]
_ST,_SC=_CL["statement"],_CL["score"]
_PM=open(os.path.join(_D,"proof.md"),encoding="utf-8").read()
_T7=[l for l in _PM[_PM.index("## 7."):].split("\n")if l.startswith("| ")and"|"in l[2:]][2:]
ROWS7=["| %s | %.4f |"%(n,lg(v))for n,v in LEDGER]
check("proof.md 7: %d log2 (4 decimals); sum E_RUN; total E6"%len(LEDGER),_T7[:len(LEDGER)]==ROWS7 and sum(v for _,v in LEDGER)==E_RUN and _T7[len(LEDGER)]=="| **opened run E_run** | **%.6f** |"%claim6_of(E_RUN)and _T7[len(LEDGER)+1]=="| **expected total E = 2^7 + (K/2^64) E_run** | **%.6f** |"%E6)
E_NOABORT=E_RUN+sum((1-PHI[_g])*(WG[_g]+HCHG_VER*EK2G[_g]+(HCHG_LS1+RESOLVE+2*REQ_SERVE)*ELISTG[_g])for _g in range(G_GROUPS))
check("proof.md 7: a run that never stops early (all PHI = 1) costs 2^%.4f units, as printed"%lg(Fraction(E_NOABORT,C)),"costs about 2^%.4f units (the certificate's E with all PHI[g] = 1)"%lg(Fraction(E_NOABORT,C))in _PM and"an opened run costs about 2^%.4f units"%lg(Fraction(E_RUN,C))in _PM)
check("stated figures (literals, submission.json, score) equal the certified ones",STATED==DECLARED and Fraction(int(_SC["numerator"]),int(_SC["denominator"]))==STATED["time_log2"]and all(re.search(re.escape(t)+r"(?![0-9])",_ST)for t in("2^123.308715787426632","2^190.602588","2^135.8120","2^115.9395","57.1815","0.39","K = %d"%STATED_K))and STATED_K==GATE_K and"K = %d"%STATED_K in _PM and"K = %d"%STATED_K in open(os.path.join(_D,"success.md"),encoding="utf-8").read())
check("declared E work per leaf %.4f: least 4-decimal bound >= E/S (%.7f)"%(DECLARED["per_leaf"],PER_LEAF_E),DECLARED["per_leaf"]>=PER_LEAF_E>DECLARED["per_leaf"]-Fraction(1,10**4))
MULX,SQRX,ADDX=105+4,80+3,9+2*9
DBL,MADD_=2*MULX+5*SQRX+10*ADDX,8*MULX+3*SQRX+10*ADDX
GEN_ITEM=256*(DBL+MADD_)+2*(INV+3*MULX+2*SQRX+10*ADDX)+64
GEN_TEST=256*2*(9+4)+64
GEN_CTRL=1600
NEGT=30
RED=9
Y_ENTRY=3*SQRX+2*MULX+2*ADDX+2*NEGT+7*RED+14
X_ENTRY=MADD_+3*SQRX+11*MULX+2*ADDX+2*NEGT+17*RED+32
check("generator: doubling 903, mixed addition 1,391; with tests, control, max(Y_ENTRY, X_ENTRY) <= %d, the allowance"%(GEN_ITEM+GEN_TEST+GEN_CTRL+max(Y_ENTRY,X_ENTRY)),DBL==903 and MADD_==1391 and GEN_ITEM+GEN_TEST+GEN_CTRL+max(Y_ENTRY,X_ENTRY)<=GEN_MUL,"%.4fx"%(math.floor(10**4*GEN_MUL/(GEN_ITEM+GEN_TEST+GEN_CTRL+max(Y_ENTRY,X_ENTRY)))/10**4))
check("permanent X generator: [a]P + R_g, tests, control, four affine-store STPs fit GEN_MULX",GEN_ITEM+GEN_TEST+GEN_CTRL+4<=GEN_MULX,"%d <= %d"%(GEN_ITEM+GEN_TEST+GEN_CTRL+4,GEN_MULX))
XSET_ITEM=4*MULX+3*SQRX+10*ADDX+2*NEGT+6*6*MULX+64
check("X setup: 2a^3, 2a, a^2, a^4, y^2, 2y^2, 2y, sigma, pi - 7 sigma, delta, 2 delta, -T(a), -T(beta a), 6 steps a table fit XSETUP",XSET_ITEM<=XSETUP,"%d <= %d"%(XSET_ITEM,XSETUP))
ENTRY_ITEM=4+4+2+2+11+7+5+3+9
check("pair entry X words (corrections, y^2: 2 LDP + 2 STP each; table, sigma bases; P1's block of 6 line bases, 2 YLS pointers 7 + 4 + 4 (LSR, ADD shifted, ADD 128, ADD 64), delta's STR, 2 delta 4 + 3; loop 3; address words at most 9: odd TBZ, SUB 448, TBZ, SUB 640, ADD 384, STP +232, SUB 64, STR +40, B; even TBZ, ADD 1088, STR +360) fit ENTRY_X",ENTRY_ITEM<=ENTRY_X)
BARRIER_ITEM=7*(8+5)+5*(6+MSG_LINK)+1+20
BATCH_ITEM=2+4+60+(1+4+4+4+2+1)+(2*4+1)+(4+1+2+1+1+1)+2*(4+1)+2*(2+1+3)+2*(1+4+1+1)+(2+2+1+1)+2
check("batch set-up: base %d + 16 group-controller instructions fit 176"%BATCH_ITEM,BATCH_ITEM+16<=BATCH_SETUP)
check("phase control: 12 kernel calls a group <= 33 each, stop test and loop 8, h2s's p word 4, over %d groups < %d"%(G_GROUPS,PHASE_CTRL_G),G_GROUPS*(12*33+8+4)<PHASE_CTRL_G)
check("barrier: at most 232 per node: 7 waits (5 barrier words, 2 count words) of LDAXR, CMP, B.EQ, WFE, B, re-check, one stale-event pass of 5; 5 messages of 24, DSB, 20 set-up (2^10 is > 4x)",BARRIER_ITEM==232 and BARRIER>4*(NODES+NG)*BARRIER_ITEM)
EXCH_ITEM=16+(4*(14+MSG_LINK)+16)*(8*REQ_CAP)
check("phase-4 count exchange 16, 144 per message (proof.md 7)",EXCH_ITEM==EXCH_DENSE+EXCH_REQ==16+144*2**45 and NODES*2*CAND_CAP==NG*8*REQ_CAP==2**109 and EXCH_ROW==(NODES+2*NG)*4*EXCH_DENSE and EXCH_REQ_RUN==NG*4*144*8*REQ_CAP==NODES*4*144*2*CAND_CAP,"2^%.4f"%lg(EXCH_ITEM))
check("listings, candidates and requests (2 each, REQ_SERVE): main pairs per group under PHI at rho (KB(156) + 1/r); same-batch JL KB(156) over SB_PAIRS; dependent terms global; R_worst charges each source's cap once",sum(_s*_gp+_s*(_s-1)//2 for _s,_gp in zip(SG,[NG*LEAVES*m*Q2 for m in MQ1[:-1]]))==SG_FULL*(SG_FULL-1)//2 and all(ELISTG[_g]==RHO*(SG[_g]*NG*LEAVES*MQ1[_g]*Q2+SG[_g]*(SG[_g]-1)//2-SAMEX[_g])*(KB(156)+Fraction(1,r))for _g in range(G_GROUPS))and ELIST==_PS(ELISTG)+RHO*(JL*SB_PAIRS+2*SG_FULL+PARTNER_PAIRS+FIXED_PAIRS)*KB(156)+_SIPT*KB(156)*SG_FULL<ELIST_G and E_REQ==2*REQ_SERVE*(ELIST+DEFL/HCHG_LS1)<NG*REQ_CAP*REQ_SERVE and CAPS_W-NG*REQ_CAP*REQ_SERVE==NODES*CAND_CAP*RESOLVE+NG*DELTA_F*6,"2^%.4f of 2^%.4f"%(lg(ELIST),lg(ELIST_G)))
check("a broadcast message: <= 26 with in-flight count and delivery, within its price 40",2+3+3+MSG_LINK<=MSG)
CAND_LEAVES=INV+5*SQRX+18*MULX+20*ADDX
EPOCH_DEC=100
assert EPOCH_DEC==8+12*6+8+12 and 2**12>G_GROUPS
CAND_ITEM=4*8*2*MSG+4*CAND_LEAVES+36*600+384*200+2*256*(DBL+MADD_)//2+INV+2*EPOCH_DEC
check("candidate: messages, both records' decode, six leaves each, 36 unit tests, solve mod r, [x]P = Q = %d; charged above twice that"%CAND_ITEM,RESOLVE==2*CAND_ITEM+1,"%.3fx"%(math.floor(10**3*RESOLVE/CAND_ITEM)/10**3))
PROGRAM_BYTES=2**18
KERNEL_INSTR=21368+13802+1519
check("program delivery: 2^18 bytes (%d kernel instructions, %d bytes, + <= 2^16 + 64 generator, resolver, control) = 2^13 messages, 3 hops at 40 < 2^20"%(KERNEL_INSTR,4*KERNEL_INSTR),4*KERNEL_INSTR+2**16+64<=PROGRAM_BYTES and PROGRAM_BYTES//32*3*MSG+128<NODE_INIT and 256*(MULX+SQRX)+2**17*48+2*GEN_MUL+2**16<2**23<NODES)
ANALYTIC_E=NG*(4*Q1GD+MX+1)*COIN_SCALAR+NG*MX*(GEN_MULX+XSETUP)+NG*Q2*ENTRY_X+NG*HEADER+(NODES+NG)*NODE_INIT+NG*Q1GD*3*GEN_MUL+NG*GEN_MUL+NG*Q1GD*(BATCH_SETUP+INV_CALL)+(1+G_GROUPS)*2*(SIDE+NODES+NG)*MSG+E_CAND+E_REQ+EXCH_REQ_RUN+sum(PHI[_g]*EXCH_ROW for _g in range(G_GROUPS))+NBAR_ALLOW*BARRIER+(NODES+NG)*PHASE_CTRL_G+OUTPUT+MINE_W
check("analytic allowances (C, T, S, requests, exchange messages at lifetime caps, Phase 0) are below 2^-9.88 of E_run",ANALYTIC_E<lo2(-988)*E_RUN,"2^%.3f"%lg(Fraction(ANALYTIC_E,E_RUN)))
print("RESULT q1=%d G=%d q=%d mu=%.6f lo=%.10f score=%.15f E_run_log2=%.6f d_f_log2=%.3f d_p1_log2=%.3f d_p2_log2=%.3f d_rq_log2=%.3f d_rl_log2=%.3f d_cl_log2=%.3f lam_log2=%.3f q1EFIX_log2=%.4f K=%d worst6=%.6f preproc_log2=%.7f perleaf=%.7f FAILED=%s"%(Q1GD,G_GROUPS,Q1G[0],float(TFULL["mu"]),int(TFULL["lo"]*10**10)/1e10,float(CLAIM_E),lg(E_RUN),lg(TFULL["d_f"]),lg(TFULL["d_p1"]),lg(TFULL["d_p2"]),lg(TFULL["d_rq"]),lg(TFULL["d_rl"]),lg(TFULL["d_cl"]),lg(TFULL["lam_int"]),lg(Q1GD*E_FIX),GATE_K,float(DECLARED_WORST6),lg(Fraction(E_PRE,C)),float(PER_LEAF_E),";".join(FAILED)if FAILED else"none"))
if not FAILED:print("PASS")
if __name__=="__main__":sys.exit(1 if FAILED else 0)

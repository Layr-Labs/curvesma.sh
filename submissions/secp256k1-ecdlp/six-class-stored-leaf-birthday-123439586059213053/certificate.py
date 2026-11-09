"""Exact certificate for proof.md (parameters, E-work ledger, claim, memory)."""
import functools
import math
import os
import sys
from decimal import Decimal,ROUND_CEILING,getcontext
from fractions import Fraction
FAILED=[]
def check(name,ok,detail=""):
 print("%-66s %s%s"%(name,"ok" if ok else "FAIL",("  "+detail)if detail else ""))
 if not ok:
  FAILED.append(name)
def lg(x):
 f=Fraction(x)
 return math.log2(f.numerator)-math.log2(f.denominator)
def cdiv(a,b):
 return-((-a)//b)
def exp_interval(x,terms=60):
 if abs(x)>1:  # e^x = (e^(x/2^k))^(2^k), squared and rounded outward to 2^-512
  k=(abs(x).numerator//abs(x).denominator).bit_length()
  lo,hi=exp_interval(x/2**k,terms)
  assert lo>0
  for _ in range(k):
   lo,hi=lo*lo,hi*hi
   lo,hi=Fraction((lo.numerator<<512)//lo.denominator,2**512),Fraction(cdiv(hi.numerator<<512,hi.denominator),2**512)
  return lo,hi
 a,b,n=x.numerator,x.denominator,terms  # s = sum_(k<=n) x^k/k! over one denominator b^n n!, t = x^n/n!
 num,ak,bk,fk=0,1,b**n,math.factorial(n)
 for k in range(n+1):
  num+=ak*bk*fk
  ak,bk,fk=ak*a,bk//b,fk//(k+1)
 den=b**n*math.factorial(n)
 s,t=Fraction(num,den),Fraction(a**n,den)
 rem=abs(t*x)/(terms+1)*3
 return s-rem,s+rem
def ln2_interval(bits=300):
 z,z2=Fraction(1,3),Fraction(1,9)
 tot,pw,k=Fraction(0),z,0
 while True:
  tot+=pw/(2*k+1)
  rest=pw*z2/((2*k+3)*(1-z2))
  if rest<Fraction(1,2**bits):
   return 2*tot,2*(tot+rest)
  pw*=z2
  k+=1
LN2=ln2_interval()
@functools.lru_cache(maxsize=None)  # exact interval for 2^c
def pow2_interval(c):
 n=c.numerator//c.denominator
 f=c-n
 lo,hi=exp_interval(f*LN2[0])[0],exp_interval(f*LN2[1])[1]
 sc=Fraction(2)**n
 return lo*sc,hi*sc
# 1. constants
p=2**256-2**32-977
r=int("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141",16)
C0=2**32+977
BETA_P=int("7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee",16)
C=1000
NODES=2**72
SIDE=2**36
Q2=2**25
D1,D2=36,36
MH,ML=22,28  # H words compared from bit 22, L words (d1, d2) from bit 28 (proof.md 3.3)
KEYBITS=156-2*MH-ML  # compared bits 228 - 2MH - ML = 156 (emit 7), less d1, d2
THRESH=Fraction(39,100)  # the rules' minimum
LEAVES=6
MX,YM,YN=1024,8192,4  # X's per node; a batch's Y = C + [m t + n s]P, m < YM, n < YN
L_Y=YM*YN
def fwd_count(q2):  # proof.md 3.1 (the check below)
 assert q2%32==0
 return 5+2180*(q2//32)+q2//MX
def bwd_count(q2):
 assert q2%16==0
 return 6+15395*(q2//16)+5*((q2+MX-1)//MX-1)
FIX_SITES=dict(SHX=9,FDX=6,F5=2,FS=3,FADDW=2,BLK=16)  # per pair
FIX_EXTRA=dict(SHX=3,FDX=3,F5=3,FS=4,FADDW=4,BLK=6)
EMIT=7
EMIT_INC=LEAVES*EMIT
KEY_INC,KEY_FIX=Fraction(267,2),4
KEY_LSL=1  # per incidence: LSL x18 = x3<<30 (signless w3)
K=256
PB1C,PB2C,BLKC=1<<50,1<<51,1<<14
FLUSH=9
KEY_BLOCK_DROP,H2S_BLOCK_DROP=7,9
LOOP8=2
JPATH=2
SJOIN=1  # a sign with a j = 1 or 2 leaf ends in one B, charged on every sign
INV=22549
TINIT_ENTRY=3
ECOPY_ENTRY,KCALL=4,33
FIN_DEST,FIN_SENT,FIN_HELD,FIN_FIX=14,4,2,2
H2S_ARR,H2S_ONE=6,1
H2S_FULL=6+2*H2S_ONE
H2S_EMPTY=7  # count 0: MOVN, LDR, LDR, MOV, TBNZ, CBZ, RET
H2S_ARR_FIX=14+2*H2S_ONE  # + CBZ, SUB, ADD, 3 TBZ, CMP, B.HS; slot 1, a last slot alone
def h2s_count(n):  # h2s on an n-record block
 if n==K:
  return H2S_FULL+H2S_ARR*K
 if n==0:
  return H2S_EMPTY
 return 14+H2S_ONE*(1+(n%2==0))+LOOP8*((n-1)//8)+H2S_ARR*n
check("h2s: 6 per record (LD64B, ADD per two slots), +1 per slot alone (two LDPs); full block 8 fixed, else <= 16 + 2 per "
   "8-record pass after slot 1: <= 16 + 62",
   h2s_count(K)==6*K+8 and max(h2s_count(n)-H2S_ARR*n for n in range(K))==H2S_ARR_FIX+LOOP8*(K//8-1)
   and H2S_EMPTY<=H2S_ARR_FIX)
INIT3_BKT=4
check("key pass: 133.5 = 2 signs * (LD64B of the gamma, beta'' block, 2 LDP of alpha', SUB of x3, 26 word sums, 3 leaves * "
   "(index ADD, CMP, B.HS, CMN, B.HS + emit 7)) + LSL (signless w3) + TBZ per two incidences",
   KEY_INC==2*(1+2+1+26+3*(5+EMIT))+KEY_LSL+Fraction(1,2))
check("bwd: 1924 per pair (BWDP: LK2 E, F, G of both incidences 98 + 86 + 86, LD64B of e, acc, W2S, C1BS 157, two INC with P1 248, "
   "3 TFP4S 97), 3 per 16 incidences, 5 per inner Y switch, 6 fixed; fwd: 272 per two pairs (PFWD2: EFWD2 97, P1X, P1XB), 4 per 32 "
   "incidences, 1 per Y switch, 5 fixed",
   Q2==MX*L_Y and bwd_count(Q2)==6+(8*1924+3)*(Q2//16)+5*(Q2//MX-1)
   and fwd_count(Q2)==5+(8*272+4)*(Q2//32)+Q2//MX)
check("j = 1, 2 paths (ADD, CMP, B.HS; ADD of S, CMN, B.HS; MVN of sigma; emit) exceed j = 0's (ADD, CMP, B.HS, CMN, B.HS; "
   "emit) by 2; tree copies end in one B per sign (SJOIN)",
   (3+3+1+EMIT)-(5+EMIT)==JPATH and SJOIN==1)
check("emit: AND or BIC (8*d1), LDADD of its LC entry (the slot), w3 in one ADD x18,W,LSR #28 (j = 0, 2) or EON "
   "x18,S,LSR #28 (j = 1) on x18 = x3<<30, K2 by BFXIL, "
   "K3 the second part's H word, two STPs, TBZ of bit 13 for slot k",EMIT==1+1+1+1+2+1)
DONE=2
DLV=6
MSG_LINK=DONE+DLV
BLOCK=FLUSH+DONE
FINAL=FIN_DEST+FIN_SENT+DONE+MSG_LINK
check("a hop-1 block: flush 9 + in-flight 2 when sent, > 7 when dropped; hop-1 fin per destination <= 28 (14, send 4, "
   "in-flight 2, E1 message 8)",
   BLOCK==11 and BLOCK>KEY_BLOCK_DROP and FINAL==28 and FINAL>FIN_DEST+FIN_HELD+MSG_LINK)
GEN_MUL=645410
GEN_MULX=GEN_MUL+4  # + 4 STPs of X's affine store
YG=(-23285,5,22598,1606,45680)  # ygen: A + B m + C n + D m n + E ceil(log2 m)
def ygen_count(m,n):
 return YG[0]+YG[1]*m+YG[2]*n+YG[3]*m*n+YG[4]*(m-1).bit_length()
YGEN=ygen_count(YM,YN)+3*GEN_MUL  # with C, T, S
W7=(37,)*4+(36,)*3  # W7: seven windows of v, window k bits W7O[k].. (proof.md 3.1)
W7O=[sum(W7[:k])for k in range(7)]
TBX=lambda L,n:25+26*(L-1)+22*L*n  # one tblx call: L lines of n rows
TBXP,TBXD=TBX(4,2**37)+TBX(3,2**36),TBX(2,2**37)+TBX(2,2**36)  # pair: lines 0-3, 4-6; delta: lines 0-1 (slots 0-3), 2-3 (4-7)
NPAIR_T,NSINGLE_T=3*MX//2+MX//4+MX,MX//2
XTBLD_ALL=NPAIR_T*TBXP+NSINGLE_T*TBXD
XT_BYTES=NPAIR_T*2**46+NSINGLE_T*2**45
XSETUP=2**16
ENTRY_X=40
COIN_SCALAR=64
BATCH_SETUP=176
INV_CALL=16
NODE_INIT=2**20
HEADER=72
BT_BASE,BT_W=0xfff8<<48,43
BT_WIN=(0,1,3,2,4,5,3)  # PAIRB6: window of each beta slot
TBLD=4+7*12+20*(6*2**BT_W+2**(256-5*BT_W))
check("beta table (PAIRB6): a 256-byte row's 7 slots hold windows %s (window 3 in slots 2, 6: bit 128 at address bit 7; "
   "slot 5 41 bits); tbld 2^%.4f"
   %(BT_WIN,lg(TBLD)),sorted(set(BT_WIN))==list(range(6))and BT_WIN.count(3)==2 and BT_WIN.index(3)==2 and BT_WIN[6]==3
   and 32*len(BT_WIN)<=256 and (6-2)*32==1<<7 and 256-5*BT_W==41 and BT_WIN[5]==5 and 43*3==129)
MSG=40
CAND_CAP=2**35
RESOLVE=1625371
REQ_CAP=2**33
REQ_SERVE=2**10+GEN_MUL  # regenerates Y
RELAY_CAP=2**40
EXCH_G=12*SIDE+4*(14+MSG_LINK)*(8*REQ_CAP)
BARRIER=2**82
OUTPUT=2**20
DELTA_F=2**62
FIX_BATCH=sum(FIX_SITES.values())*(Q2//2)
CAP_P=2**22-2**17
P1=Fraction(6*(r-1)+1,r*r)
TWO_BOUND=Fraction(2,r-1)
SHX_BAD=2**64*2**33
FOLD_BAD=2**66*C0
SQ_BAD=C0*C0+2*C0
WEAK_BAD=C0  # P1: N's true value below C0
H=Q2//2
SITES=[
 ("fwd: row shifts of acc (uniform on F_p^* by rho)",3*H,SHX_BAD,1),
 ("fwd: fold of acc",H,FOLD_BAD,1),
 ("bwd: row shifts of R = 1/acc",3*H,SHX_BAD,1),
 ("bwd: fold of inv_e = 1/e (e of degree 2 in y_Y^2)",H,FOLD_BAD,12),
 ("bwd: fold of R' = 1/acc_prev",H,FOLD_BAD,1),
 ("bwd: row shifts of D_odd = y_Y^2 - y_X^2",3*H,SHX_BAD,6),
 ("bwd: fold of 1/D_even",H,FOLD_BAD,6),
 ("bwd: folds of N = 2 y_X y_Y w (y_Y a root of a quartic)",Q2,FOLD_BAD,12),
 ("bwd: FS and FADDW of M, N on weak operands",2*Q2,WEAK_BAD,12),
 ("bwd: square fold of w = 1/D^2 (two D)",Q2,SQ_BAD,12),
]
SITES_PER_PAIR=dict(SHX=9,FDX=6,F5=2,WEAK=4)
D_SITES=sum(Fraction(n*bad*mult,r-1)for _,n,bad,mult in SITES)
FIX_ROW=C0*2**192
FIX_FOLD=(2**66+1)*2**128*C0
FIX_SQ=2**193
FIX_BLK=30*C0*2**192  # BLK6 (LF3): t <= 5, low 64 bits >= 2^64 - t C0, 2 words w
FIX_LF=36*C0*2**192  # LF (LF3): c <= 8, low 64 bits >= 2^64 - c C0
FIX_SUB=C0*2**192
FIXUPS=[
 ("fwd: lookup fold of e","BLK",H,FIX_LF,12),
 ("fwd: row shifts of acc","SHX",3*H,FIX_ROW,1),
 ("fwd: fold of acc","FDX",H,FIX_FOLD,1),
 ("bwd: row shifts of R","SHX",3*H,FIX_ROW,1),
 ("bwd: folds of inv_e and R'","FDX",2*H,FIX_FOLD,12),
 ("bwd: FS of D_odd","FS",H,FIX_SUB,6),
 ("bwd: row shifts of D_odd","SHX",3*H,FIX_ROW,6),
 ("bwd: fold of 1/D_even","FDX",H,FIX_FOLD,6),
 ("bwd: lookup fold of 1/D_odd","BLK",H,FIX_LF,6),
 ("bwd: folds of N","FDX",Q2,FIX_FOLD,12),
 ("bwd: lookup folds of M = 1/D + 2A w and L = 2 y_X w (M: two D)","BLK",2*Q2,FIX_LF,12),
 ("bwd: FS of u^2 = M - N","FS",Q2,FIX_SUB,6),
 ("bwd: FADDW of v^2 = M + N","FADDW",Q2,FIX_SUB,6),
 ("bwd: square fold of w","F5",Q2,FIX_SQ,12),
 ("bwd: window-sum folds of beta u^2 and beta v^2","BLK",2*Q2,FIX_BLK,6),
 ("bwd: lookup folds of E, F, G","BLK",3*Q2,FIX_LF,8),
]
E_FIX=sum(Fraction(n*bad*mult,r-1)for _,_,n,bad,mult in FIXUPS)
KBF=2**13  # KB's widening factor 1 + 1/KBF (3.3)
KB=lambda fixed:Fraction(6*2**(256-fixed)*(KBF+1),(r-1)*KBF)
RHO=Fraction(r,r-1)**4
FIXED_LEAVES=2*LEAVES
PARTNER_PAIRS=2*10+72*L_Y*L_Y  # 6.4: pairs without a home factor
LAMBDA_FREE=3*(YM+YN-2)**2<r
YGEN_DEF=Fraction(2*(YM-1+YM*(YN-1)),r)
D_BATCH=D_SITES+Q2*Fraction(6,r-1)+YGEN_DEF+Fraction(L_Y*L_Y,r-1)  # sites, D = 0, ygen, Y_v = +-Y_w
FA,FB=int("114ca50f7a8e2f3f657c1108d9d44cfd9",16),int("e4437ed6010e88286f547fa90abfe4c3",16)
def isum(*rs):
 return sum(x[0]for x in rs),sum(x[1]for x in rs)
def table(raw):  # 3.3's error table (W0, W1, N' over t) from a word's range
 cor=isum(raw,(-1,0))
 u3,uc3=(raw,raw,(-2-2*raw[1],-2*raw[0])),(cor,cor,(-2-2*cor[1],-2*cor[0]))
 wr=[[isum(uc3[m],uc3[(m+t)%3],u3[(m+2*t)%3],(0,2))for t in range(3)]for m in (0,1)]
 nr=[(-wr[0][t][1]-wr[1][t][1]-1,-wr[0][t][0]-wr[1][t][0]+1)for t in range(3)]
 return [(min(x[0]for x in rs),max(x[1]for x in rs))for rs in (wr[0],wr[1],nr)]
TFP_H=TFP_L=((0,1),(0,7),(-1,1))  # e, carry, the unit the rows' errors E (|E| < 2^163) can move H by
ERRH,ERRL=table(isum(*TFP_H)),table(isum(*TFP_L))
WDH,WDL=sum(hi-lo for lo,hi in ERRH),sum(hi-lo for lo,hi in ERRL)
WH,WL=max(hi-lo for lo,hi in ERRH)+1,max(hi-lo for lo,hi in ERRL)+1
check("3.3's word ranges (CHAIN): H, L in [-1, 9] (e 0..1, carry 0..7, rows' error |E| < 2^64 C0 (3*2^64 + 2) < 2^163: one "
   "unit); error tables W0, W1 [-24, 29], N' [-59, 48]: D' widths 213, widenings 108",
   isum(*TFP_H)==isum(*TFP_L)==(-1,9)and 2**64*C0*(3*2**64+2)<2**163 and 2**(64+163)<p
   and ERRH==ERRL==[(-24,29),(-24,29),(-59,48)]and (WDH,WDL,WH,WL)==(213,213,108,108))
(ELO,EHI)=ERRH[0]
NWIN=2*(EHI-ELO)+(1-2*ELO)+(3+2*EHI)+2*(EHI+1-ELO)+2*(3*EHI+3-3*ELO)  # 3.3's W1-W4, units of 2^192 values of x
D_WJ=TWO_BOUND*NWIN*2**192
D_LAB=TWO_BOUND*(WDH*2**(256-MH)+WDL*2**(256-ML))
D_SYM=TWO_BOUND*3*((WH-1)*2**(256-MH)+(WL-1)*2**(256-ML))  # D' closed under the class (6.1)
P1_GOOD=Fraction(6*(r-1),r*r)*(1-3*D_WJ-D_SYM)
check("lattice lemma: basis (B, A - B), (A - B, -A) of {(v, u): u = beta v mod p}, |det| = p, parallelogram 2A - B wide "
   "in u and A in v; boxes of the compared bits widened by 3.3's ranges: factor below 1 + 2^-13",
   FA*FA-FA*FB+FB*FB==p and (FB*BETA_P-(FA-FB))%p==0 and ((FA-FB)*BETA_P+FA)%p==0
   and (2**MH+WH-1)**2*((2**ML+WL-1)*2**128+2*FA)*((2**28+WL-1)*2**128+4*FA-2*FB)<p*2**(2*MH+ML+28)*(1+Fraction(1,KBF)))
# 6.4 canonical routing cells
GAPE=max(3*EHI+3,EHI+1-ELO,1-3*ELO)  # KEY's gap, units of p/2^64
PAD=GAPE
WRAPT=64
def npad(i,m,d=PAD):  # canonical cells in macrorow i
 if not 0<=i<m:
  return 0
 lo=max(0,i-d)
 A,B=(lo,min(m-1,m+d-2*i)),(max(lo,m-d-i),min(m-1,2*m+d-2*i))
 n=sum(max(0,y-x+1)for x,y in (A,B))
 return n-(max(0,min(A[1],B[1])-max(A[0],B[0])+1)if A[0]<=A[1]and B[0]<=B[1]else 0)
def npad_bp(m,d=PAD,e0=3):  # npad is affine between these (its lines' crossings +-e0)
 L=[(0,0),(1,-d),(0,m-1),(-2,m+d),(-1,m-d),(-2,2*m+d)]
 return sorted({0,m}|{t+e for a,b in L for c,g in L if a!=c for t in [(g-b)//(a-c)]for e in range(-e0,e0+1)if 0<=t+e<=m})
def ncells(m,d=PAD):  # #{(i,l): l >= i-d, 2i+l <= m+d or (i+l >= m-d, 2i+l <= 2m+d)}
 f=lambda i:npad(i,m,d)
 bp=npad_bp(m,d)
 n=0
 for s,e in zip(bp,bp[1:]):
  h=e-1
  u,v,w=f(s),f(h),f((s+h)//2)
  assert (v-u)*((s+h)//2-s)==(w-u)*(h-s)
  n+=(u+v)*(e-s)//2
 return n
NC=ncells(2**64)
NW=(8*WRAPT+8)*2**64  # cells with a strip row, column or z-diagonal (6.4)
CELLU=(2**28+WL)*2**128  # v's of a widened digit cell in a macrocolumn
ROW=2**192
DPAR=max(FA,2*FA-FB)
KB_R36=Fraction(6*(NC+NW)*(CELLU+2*DPAR)*(ROW+DPAR),p*(r-1))
KB_H72=Fraction(6*(NC+NW)*(CELLU+2*DPAR)**2,p*(r-1))
check("canonical cells: KEY gap %d (p/2^64 units) from H ranges [%d, %d]; strips %d; NC (pad %d, brute-checked) + NW < "
   "(1 + 2^-52) 2^128/3; KB_R36, KB_H72 = %.7f, %.7f x 2^220, 2^184 /(r-1) < KB(36), KB(72)/2.99"%(GAPE,ELO,EHI,WRAPT,PAD,KB_R36*(r-1)/2**220,KB_H72*(r-1)/2**184),
   ERRH[0]==ERRH[1]and WRAPT>=max(EHI,1-ELO,1-2*ELO,2+2*EHI)
   and all(ncells(m,d)==sum(1 for i in range(m)for l in range(m)if l>=i-d and (2*i+l<=m+d or (i+l>=m-d and 2*i+l<=2*m+d)))
   for m,d in [(m,8)for m in range(1,61)]+[(m,PAD)for m in (200,263,331)])and 3*(NC+NW)<(1+Fraction(1,2**52))*2**128
   and KB_R36*Fraction(299,100)<KB(36)and KB_H72*Fraction(299,100)<KB(72))
check("wrong class index rate per leaf (W: %d 2^192 values of x) < 2^-53.66"%NWIN,D_WJ<pow2_interval(Fraction(-5366,100))[0],"2^%.4f"%lg(D_WJ))
check("site table per pair: 9 row shifts, 6 folds, 2 square folds, 4 weak-operand FS/FADDW; lookup and window-sum folds exact (BLK-type)",
   sum(n for _,n,b,_ in SITES if b==SHX_BAD)==SITES_PER_PAIR["SHX"]*H
   and sum(n for _,n,b,_ in SITES if b==FOLD_BAD)==SITES_PER_PAIR["FDX"]*H
   and sum(n for _,n,b,_ in SITES if b==SQ_BAD)==SITES_PER_PAIR["F5"]*H
   and sum(n for _,n,b,_ in SITES if b==WEAK_BAD)==SITES_PER_PAIR["WEAK"]*H and len(SITES)==10)
check("bad sets: rows 2^97, folds below 2^98.02, square folds below 2^64.01, weak operands C0",
   SHX_BAD==2**97 and FOLD_BAD<pow2_interval(Fraction(9802,100))[0]and SQ_BAD<pow2_interval(Fraction(6401,100))[0])
FIX_RUNS={m:sum(n for _,ms,n,_,_ in FIXUPS if ms==m)for m in FIX_SITES}
check("fixup table per pair: 9 row shifts, 6 folds, 2 square folds, 3 FS, 2 FADDW, 16 BLK-type (4 window sums, 12 lookup folds)",
   FIX_RUNS=={m:FIX_SITES[m]*H for m in FIX_SITES}
   and all(bad==dict(SHX=FIX_ROW,FDX=FIX_FOLD,F5=FIX_SQ,FS=FIX_SUB,FADDW=FIX_SUB)[m]for _,m,_,bad,_ in FIXUPS if m!="BLK"))
check("lambda has no integer ratio with 0 < |a|,|b| <= YM+YN-2 = 8194 (0 < a^2 +- ab + b^2 <= 3*8194^2 < r)",LAMBDA_FREE)
check("a call's fixup word (C0 per fixup, at most 19 q2 fixups) stays below 2^63",FIX_BATCH*C0<2**63,
   "2^%.3f"%lg(FIX_BATCH*C0))
check("label-error rate per leaf < 2^-13.24",D_LAB<pow2_interval(Fraction(-1324,100))[0],"2^%.3f"%lg(D_LAB))
check("batch defect rate (site table, D = 0, ygen, coincident Y) < 2^-126.99",D_BATCH<pow2_interval(Fraction(-12699,100))[0],
   "2^%.3f"%lg(D_BATCH))
G_GROUPS=780  # equal groups
CAP_PG=2**20-2**11
def group_sizes(q1):
 assert q1%G_GROUPS==0
 q=q1//G_GROUPS
 return [q]*G_GROUPS,[g*q for g in range(G_GROUPS+1)]
HB=58
HCHG_NEW=4
HCHG_VISIT=4
HCHG_END=1
HCHG_VER=5
HCHG_LS0=2
HCHG_LSA=1
HCHG_LS1=4+HCHG_LSA
check("join per record (BFI, SWP, STP, CBZ) 4; per visit (LDP, CMP, B.EQ, CBZ or CBNZ) 4; B after an even visit count 1; "
   "a K2 hit's K3 check (ADD, LDR, EOR, TST, B.NE) 5",
   HCHG_NEW==4 and HCHG_VISIT==4 and HCHG_END==1 and HCHG_VER==5)
DLVC,DLVS1=2,2
DLVH_FAST=15+2*DLVS1
DLVH_SEG,DLVH_SLOW_FIX=7,24
def dlvh_seg(n,m):  # m slots from slot n
 if m==0:
  return 1
 head=n%2
 rest=m-head
 return DLVH_SEG+head+DLVC*m+DLVS1*(head+rest%2)+LOOP8*(rest//8)
def dlvh_copy(m1,m2):
 if m1==0 and m2==K:
  return DLVH_FAST+DLVC*K
 return DLVH_SLOW_FIX-(m1>0)+dlvh_seg(1,m1)+dlvh_seg(1+m1,m2)+dlvh_seg(1+m1+m2,K-m1-m2)
DLVHB=dlvh_copy(0,K)
DLVHS=max(dlvh_copy(a,b)for a in range(K+1)for b in range(K+1-a)if (a,b)!=(0,K))
BLOCKH=FLUSH+DONE+DLVHB
FINALH=FIN_DEST+FIN_SENT+DONE+DLVHS+MSG_LINK+(DLVHS-DLVHB)
check("dlvh: copy 2 per slot paired, 4 alone; full block 15 fixed; else <= 24, per segment 1 if empty, else 7, +1 after an "
   "odd first slot, +2 per 8-slot pass; DLVHS the maximum over splits; hop-2 final: slow + surcharge",
   DLVHB==531 and DLVHS==633 and BLOCKH==542 and BLOCKH>H2S_BLOCK_DROP and FINALH==FINAL+DLVHS+DLVHS-DLVHB)
ZW=16384
ZERO_PASS=ZW//16+2
HOME_MAX=SIDE*CAP_P
STOP_CHECK=5
NBAR_COUNT=4+17*G_GROUPS+1
NBAR_ALLOW=13265
PHASE_CTRL_G=318241
KB130=KB(D1+D2+HB)
# C130 (success.md 6.6): ROWMAX130 = max_c sum_k npad(c + k 2^58), affine between candidates, so peaking at one
def rowsum130(c,m=2**64,sh=2**58):
 return sum(npad((c+k*sh)%m,m)for k in range(m//sh))
def rowmax130(m=2**64,sh=2**58):
 cs=sorted({0,sh-1}|{(t+e)%sh for t in npad_bp(m,PAD,4)for e in range(-3,4)})
 v=[rowsum130(c,m,sh)for c in cs]
 mid=all((v[i+1]-v[i])*((a+b)//2-a)==(rowsum130((a+b)//2,m,sh)-v[i])*(b-a)for i,(a,b)in enumerate(zip(cs,cs[1:]))if b-a>1)
 return max(v),mid,len(cs)
ROWMAX130,MID130,NC130=rowmax130()
TOY130=[(rowmax130(2**k,2**(k-6))[0],max(rowsum130(c,2**k,2**(k-6))for c in range(2**(k-6))))for k in (12,14)]
LIP=max(abs(npad(i+1,2**64)-npad(i,2**64))for b in npad_bp(2**64)for i in (b-1,b)if 0<=i<2**64-1)  # npad's largest step
B130=2**MH+WH-1  # exact rows of one fixed H block
CNT130=Fraction(B130**2,2**(2*MH))*(ROWMAX130+64*LIP*B130)+256*B130**2+2*WRAPT*2**(64-2*MH)*B130+64*B130*(6*WRAPT+8)
KB_C130=Fraction(6*CNT130*(CELLU+2*DPAR)**2,p*(r-1))
check("C130: the row-sum maximum over %d candidates (midpoint-affine) = brute force at m = 2^12, 2^14 (%d, %d); "
   "ROWMAX130 = %.6f*64*2^64; npad steps <= %d"%(NC130,TOY130[0][1],TOY130[1][1],ROWMAX130/(64*2**64),LIP),
   MID130 and all(a==b for a,b in TOY130))
check("C130: walk fold bound KB_C130 = 6 CNT130 (CELLU + 2D)^2/(p(r-1)) = %.6f KB(130) = 2^%.4f"
   %(KB_C130/KB130,lg(KB_C130)),KB_C130*Fraction(29,10)<KB130)
KB64=KB(64)
@functools.lru_cache(maxsize=None)  # a group's hop-1 row (6.4, rounded up to 2^-256) and mean
def p1_row(q):
 m=LEAVES*q*Q2*KB_R36+1
 x=(m+1+COVX(q)*KB_R36)/Fraction(CAP_PG-FIXED_LEAVES-m)**2
 return Fraction(cdiv(x.numerator<<256,x.denominator),2**256),m
def COVX(q):  # 6.4: dependent ordered pairs
 return 60*(q*L_Y+MX)+72*L_Y*L_Y+12*MX*L_Y
TI=2**20
def cheb(m,v,t):  # P(Z >= t) <= v/(t - m)^2, E Z <= m < t, Var Z <= v
 assert t>m
 return v/(t-m)**2
def success_prefix(q1p,qs):
 n0p=LEAVES*q1p*Q2
 Sp=NODES*n0p
 gamma=Fraction(Sp*Sp-NODES*n0p*n0p,2)
 mu,mu_g=gamma*P1,gamma*P1_GOOD
 deg=2*LEAVES*(q1p*L_Y+Q2)*Sp
 cs=gamma*deg*(P1*P1+Fraction(36,(r-1)**2))
 mu_lo=Fraction((mu_g.numerator<<200)//mu_g.denominator,2**200)if mu_g else Fraction(0)
 e_hi=exp_interval(-mu_lo,40)[1]
 d_batch,d_wj,d_lab=D_BATCH,D_WJ,D_LAB
 d_f=q1p*E_FIX/(DELTA_F-FIX_BATCH)+q1p*d_batch
 m_p2=n0p*SIDE*KB_H72+1
 v_p2=m_p2+1+COVX(q1p)*KB_R36
 d_p2=v_p2/Fraction(CAP_P-FIXED_LEAVES-m_p2)**2
 assert list(qs)==Q1G[:len(qs)]
 d_p1,m_p1=P1CUM[len(qs)]  # the prefix's per-group rows, summed
 kb_m=KB(KEYBITS+D1+D2)
 lam_int=RHO*(Sp+FIXED_LEAVES)*kb_m
 sf=RHO*(1+47*(d_wj+d_lab))*kb_m  # 6.4's loads: the other pairs by Markov at TI
 ep=sf*(FIXED_LEAVES*Sp+PARTNER_PAIRS)
 nk,kh,kv=n0p*kb_m,KB(D1+D2),KB(D2)
 m_h=Fraction(Sp*Sp,2)*kh*kb_m  # the home
 v_h=m_h+NODES**3*nk**2*(n0p*kh+n0p**2*kh**2)+NODES**2*(2*n0p**3*kh+n0p**4*kh**2)*kb_m**2
 d_cl=cheb(m_h,v_h,CAND_CAP-4-TI)+(ep+sf*(4*n0p*Sp*kh+NODES*n0p**2*kh))/TI
 m_g=n0p*Sp*kb_m  # the source (requests <= 2 per record)
 v_g=m_g+NODES*n0p**2*nk**2
 d_rq=cheb(m_g,v_g,Fraction(REQ_CAP-8,2)-TI)+sf*n0p**2+(ep+sf*n0p**2)/TI
 m_v=2*Fraction(Sp,SIDE)*Sp*kv*kb_m  # the relay (<= 18 records per pair)
 v_v=m_v+20*SIDE*NODES**2*nk**2*(n0p*kv+n0p**2*kv**2)+8*SIDE*NODES*(2*n0p**3*kv+n0p**4*kv**2)*kb_m**2
 d_rl=cheb(m_v,v_v,Fraction(RELAY_CAP-18,18)-2**25)+(ep+sf*(4*n0p*Sp*kv+2*SIDE*n0p**2*kv))/2**25
 per_leaf=d_batch+d_f+d_p1+d_p2+d_cl+d_rq+2*d_rl  # d_wj, d_lab: in mu_g
 rows_=mu*(1+Fraction(1,r))*(Fraction(1,r)+lam_int+2*per_leaf)
 g_def=NODES*q1p*d_batch
 glob=NODES*(4*q1p+MX+1)*Fraction(1,2**254)+NODES*Fraction(MX,r)+g_def
 lo=1-e_hi-cs-rows_-glob
 return dict(mu=mu,mu_g=mu_g,lo=lo,cs=cs,rows=rows_,glob=glob,Sp=Sp,n0p=n0p,d_p1=d_p1,d_p2=d_p2,m_p1=m_p1,m_p2=m_p2,
    d_f=d_f,lam_int=lam_int,d_cl=d_cl,d_rq=d_rq,d_rl=d_rl,g_def=g_def,m_h=m_h,m_g=m_g,m_v=m_v)
Q1GD=536870880  # 780*688296
Q1G,MQ1=group_sizes(Q1GD)
P1CUM=[(Fraction(0),Fraction(0))]  # (sum, max) of p1_row over g groups
for _q in Q1G:
 P1CUM.append((P1CUM[-1][0]+p1_row(_q)[0],max(P1CUM[-1][1],p1_row(_q)[1])))
N0G=LEAVES*Q1GD*Q2
SG_FULL=NODES*N0G
TPRE=[success_prefix(m,Q1G[:g])for g,m in enumerate(MQ1[:-1])]
PHI=[Fraction(1)]+[min(Fraction(1),max(Fraction(0),Fraction(cdiv((1-t["lo"]).numerator<<256,(1-t["lo"]).denominator),2**256)))
   for t in TPRE[1:]]  # 6.6: 1 - lo(m) rounded up to 2^-256
TFULL=success_prefix(Q1GD,Q1G)
check("leaves S = 2^72*6*q1*q2 = 2^%.6f"%lg(SG_FULL),SG_FULL==NODES*6*Q1GD*Q2==TFULL["Sp"])
check("mu = |Gamma|*P1 = %.9f (cross-node pairs), good collisions mu_g = |Gamma|*P1_GOOD = %.9f (D_SYM 2^%.3f)"%(float(TFULL["mu"]),float(TFULL["mu_g"]),lg(D_SYM)),TFULL["mu"]==Fraction(SG_FULL*SG_FULL-NODES*N0G*N0G,2)*P1
   and TFULL["mu_g"]==TFULL["mu"]/P1*P1_GOOD and D_SYM<2*D_LAB)
check("Chen-Stein b1+b2 < 2^-73",TFULL["cs"]<Fraction(1,2**73),"2^%.2f"%lg(TFULL["cs"]))
check("intruder rate per colliding pair (same 156 bits) < 2^-24.8",TFULL["lam_int"]<pow2_interval(Fraction(-248,10))[0],
   "2^%.3f"%lg(TFULL["lam_int"]))
check("mean fixups per node (clean batches) q1*E_FIX = 2^%.3f < 2^33.36; cap 2^62 > 2^16 times it plus 19 q2"
   %lg(Q1GD*E_FIX),Q1GD*E_FIX<pow2_interval(Fraction(3336,100))[0]
   and DELTA_F>2**16*Q1GD*E_FIX+FIX_BATCH)
check("fixup meter rate (Markov) < 2^-28.64",TFULL["d_f"]<pow2_interval(Fraction(-2864,100))[0],
   "2^%.3f"%lg(TFULL["d_f"]))
check("pair-slot overflow rate: hop 1, the per-group rows summed, < 2^-18.3; hop 2 < 2^-18.0",
   TFULL["d_p1"]<pow2_interval(Fraction(-183,10))[0]and TFULL["d_p2"]<pow2_interval(Fraction(-180,10))[0],
   "2^%.3f 2^%.3f"%(lg(TFULL["d_p1"]),lg(TFULL["d_p2"])))
check("cap rows by the second moment (6.4), means 2^%.2f pairs at a home, 2^%.2f at a source, 2^%.2f at a relay: "
   "candidate < 2^-35.5, request < 2^-26.8, relay (18 records a pair, twice per leaf) < 2^-24.6"%(lg(TFULL["m_h"]),lg(TFULL["m_g"]),lg(TFULL["m_v"])),
   TFULL["d_cl"]<pow2_interval(Fraction(-355,10))[0]and TFULL["d_rq"]<pow2_interval(Fraction(-268,10))[0]and TFULL["d_rl"]<pow2_interval(Fraction(-246,10))[0],
   "2^%.3f 2^%.3f 2^%.3f"%(lg(TFULL["d_cl"]),lg(TFULL["d_rq"]),lg(TFULL["d_rl"])))
check("whole-run row: a computational defect in any batch < 2^-25.99",TFULL["g_def"]<pow2_interval(Fraction(-2599,100))[0],
   "2^%.3f"%lg(TFULL["g_def"]))
check("whole-run rows (defective batch anywhere, coin rejection twice) < 2^-25.99",
   TFULL["glob"]<pow2_interval(Fraction(-2599,100))[0],"2^%.3f"%lg(TFULL["glob"]))
check("opened run: success >= 1 - e^-mu - (b1+b2) - rows - global = lo_run",TFULL["lo"]>THRESH,
   "%.10f"%(math.floor(TFULL["lo"]*10**10)/10**10))
check("non-stopping bounds: PHI[0] = 1, nonincreasing, within [0, 1]",
   PHI[0]==1 and all(0<=p<=1 for p in PHI)and all(PHI[g]>=PHI[g+1]for g in range(G_GROUPS-1)),
   "PHI[1] %.4f PHI[16] %.4f PHI[31] %.4f"%(float(PHI[1]),float(PHI[16]),float(PHI[31])))
check("per-group hop-1 slots %d are a multiple of k (keep lemma)"%CAP_PG,CAP_PG%K==0)
check("a hash chain holds one home's records of one fold: at most one slot cap per slot (2^36 slots)",
   HOME_MAX==SIDE*CAP_P,"2^%.4f"%lg(HOME_MAX))
check("uncapped walks always terminate: every chain is finite (at most HOME_MAX entries)",
   HOME_MAX<2**59)
PAIRM=KB_C130+Fraction(1,r)  # C130 (success.md 6.6)
VISIT=HCHG_VISIT+Fraction(HCHG_END,2)
PJ=Fraction(2,3)+2*D_WJ+Fraction(1,r-1)
JROW_W=NODES*N0G*JPATH+NODES*N0G//3*SJOIN
check("E-charge per visit 4.5: 4 per visit plus the B after an even visit count, at most one per two visits",
   VISIT==Fraction(9,2))
SG=[NODES*LEAVES*q*Q2 for q in Q1G]
WG,WGW,SAMEX,WGJ=[],[],[],[]
LEAF_G,ROUTE_G,GROUP_G=[],[],[]
EXCH_ROW=NODES*4*EXCH_G
for _g,_q in enumerate(Q1G):
 _s=SG[_g]
 _leaf=NODES*_q*(fwd_count(Q2)+INV+INV_CALL+BATCH_SETUP+YGEN) \
  +NODES*_q*bwd_count(Q2)+int(NODES*_q*((KEY_INC-EMIT_INC)*Q2+KEY_FIX))
 _route=_s*(EMIT*K+FLUSH+DONE)//K+_s*H2S_ARR+(_s//K)*H2S_FULL \
  +NODES*(SIDE+1)*(H2S_ARR_FIX+LOOP8*(K//8-1)) \
  +_s*BLOCKH//K+NODES*(FINAL*SIDE+FIN_FIX)+NODES*(FINALH*SIDE+FIN_FIX) \
  +NODES*(INIT3_BKT*SIDE+1+TINIT_ENTRY*SIDE+1)
 _gp=NODES*LEAVES*MQ1[_g]*Q2
 _samex=NODES*_q*(6*Q2)*(6*Q2-1)//2
 _evis=(_s*_gp+_s*(_s-1)//2-_samex)*PAIRM
 _gfix=_s*HCHG_NEW+NODES*STOP_CHECK
 _group=_gfix+VISIT*_evis
 assert isinstance(_group,(int,Fraction))and _group>=0
 _jkey=JPATH*_s*PJ+_s//3*SJOIN
 WGJ.append(_jkey)
 LEAF_G.append(_leaf);ROUTE_G.append(_route);GROUP_G.append(_group)
 WG.append(_leaf+_route+_group+_jkey+EXCH_ROW)
 SAMEX.append(_samex)
 _wgroup=_gfix+_s*HOME_MAX*(HCHG_VISIT+HCHG_END+HCHG_VER+HCHG_LS1)+EXCH_ROW
 assert isinstance(_wgroup,int)
 WGW.append(_leaf+_route+_wgroup)
check("group sizes sum to the grouped q1, each >= 1; group records multiples of k",
   sum(Q1G)==Q1GD and min(Q1G)>0 and all(s%K==0 for s in SG))
ZERO_INIT=NODES*(2**HB//ZW*ZERO_PASS+1)
check("hash-head zeroing: zero16384 zeroes 16384 words a pass by 512 (DC ZVA, ADD #256), SUBS, B.NE (1026); 2^58/16384 "
   "passes plus the RET",ZERO_PASS==1026
   and ZW==512*256//8 and 2**HB%ZW==0 and ZERO_INIT==NODES*(2**(HB-14)*1026+1))
HR_INIT=NODES*(INIT3_BKT*SIDE+1)
E3_INIT=NODES*(TINIT_ENTRY*SIDE+1+KCALL)
E3_COPY=NODES*(ECOPY_ENTRY*SIDE+1+KCALL)  # charged on every run
E_PRE=NODES*NODE_INIT+NODES*(4*Q1GD+MX+1)*COIN_SCALAR+NODES*MX*(GEN_MULX+XSETUP) \
 +NODES*XTBLD_ALL+NODES*Q2*ENTRY_X+NODES*HEADER \
 +NODES*(INIT3_BKT*SIDE+1+TINIT_ENTRY*SIDE+1)+ZERO_INIT+HR_INIT+NODES*TBLD+E3_INIT
E_WORST_REST=E_PRE+E3_COPY+NODES*GEN_MUL+2*(SIDE+NODES)*MSG \
 +NODES*DELTA_F*max(FIX_EXTRA.values()) \
 +NODES*CAND_CAP*RESOLVE+NODES*REQ_CAP*REQ_SERVE \
 +NBAR_ALLOW*BARRIER+NODES*PHASE_CTRL_G+OUTPUT+G_GROUPS*2*(SIDE+NODES)*MSG
_KBM=KB(KEYBITS+D1+D2)
_SIPT=RHO*Fraction(5,2)*(1+47*(D_WJ+D_LAB))
ELIST=RHO*(Fraction(SG_FULL)*SG_FULL*_KBM
      +(2*SG_FULL+PARTNER_PAIRS)*_KBM)+_SIPT*_KBM*SG_FULL
KB136=KB(D1+D2+64)
PAIRK2=KB136+Fraction(1,r)
EK2G=[RHO*(_s*_gp+_s*(_s-1)//2-_sx)*PAIRK2 for (_s,_gp,_sx)in
    zip(SG,[NODES*LEAVES*m*Q2 for m in MQ1[:-1]],SAMEX)]
EK2SX=PAIRM*sum(SAMEX)
EK2DEP=RHO*(2*SG_FULL+PARTNER_PAIRS)*KB64+_SIPT*KB64*SG_FULL
EK2=sum(PHI[_g]*EK2G[_g]for _g in range(G_GROUPS))+EK2SX+EK2DEP
SXV=VISIT*PAIRM*sum(SAMEX)
BRACKET=RHO*(Fraction(SG_FULL)*SG_FULL*KB(D1+D2)+2*SG_FULL+PARTNER_PAIRS) \
 +_SIPT*SG_FULL
DEFV=VISIT*D_BATCH*(HOME_MAX*SG_FULL+sum(SAMEX))
DEFK2=HCHG_VER*2*D_BATCH*BRACKET+HCHG_VER*D_BATCH*HOME_MAX*SG_FULL
DEFL=HCHG_LS1*2*D_BATCH*BRACKET+HCHG_LS1*D_BATCH*HOME_MAX*SG_FULL
DEFJ=JPATH*D_BATCH*SG_FULL
E_CAND=RESOLVE*(ELIST+DEFL/HCHG_LS1)
E_FIXUP=NODES*max(FIX_EXTRA.values())*Q1GD*(E_FIX+D_BATCH*FIX_BATCH)
CAPS_W=NODES*CAND_CAP*RESOLVE+NODES*DELTA_F*max(FIX_EXTRA.values())
E_TOTAL=E_WORST_REST-CAPS_W+E_CAND+E_FIXUP+sum(PHI[_g]*WG[_g]for _g in range(G_GROUPS))+SXV \
 +EK2*HCHG_VER+ELIST*HCHG_LS1+DEFV+DEFK2+DEFL+DEFJ
R_WORST_RUN=E_WORST_REST+JROW_W+sum(WGW)
E_RUN=E_TOTAL
LO_RUN=TFULL["lo"]
GATE_K=cdiv(39*2**64*LO_RUN.denominator,100*LO_RUN.numerator)  # least K with (K/2^64) LO_RUN >= 0.39
P_GATE=Fraction(GATE_K,2**64)
GATE=2**7  # proof.md 2.4
check("run gate: K = %d is least with (K/2^64) lo_run >= 0.39 (K and K - 1 checked), K < 2^64"%GATE_K,
   P_GATE*LO_RUN>=THRESH and Fraction(GATE_K-1,2**64)*LO_RUN<THRESH and 0<GATE_K<2**64,
   "success >= 0.39 (exact rational comparison)")
E_TOTAL=GATE+P_GATE*E_RUN
R_WORST=GATE+R_WORST_RUN
# proof.md 7's rows, summing exactly to E_RUN
_PS=lambda xs:sum(PHI[_g]*xs[_g]for _g in range(G_GROUPS))
LEDGER=[
 ("node set-up: 2^20 per node (P)",NODES*NODE_INIT),
 ("coins: (4 q1 + MX + 1)*64 (P)",NODES*(4*Q1GD+MX+1)*COIN_SCALAR),
 ("permanent X's: MX*(GEN_MULX + XSETUP) (P)",NODES*MX*(GEN_MULX+XSETUP)),
 ("X tables: tblx, %d pair and %d delta tables (P)"%(NPAIR_T,NSINGLE_T),NODES*XTBLD_ALL),
 ("pair-entry X words: q2*ENTRY_X (P)",NODES*Q2*ENTRY_X),
 ("rows of beta: 72 (P)",NODES*HEADER),
 ("row LC and RC tables (P)",NODES*(INIT3_BKT*SIDE+1+TINIT_ENTRY*SIDE+1)),
 ("hash-head zeroing: zero16384 (P)",ZERO_INIT),
 ("hash-read-cursor init: init3 (P)",HR_INIT),
 ("beta tables: tbld (P)",NODES*TBLD),
 ("E3: tinit to 2^64 - 1 (P)",E3_INIT),
 ("E3 copy at group 391: ecopy (W)",E3_COPY),
 ("R_g generator",NODES*GEN_MUL),
 ("Q and stop broadcasts at 40 per message (W)",(1+G_GROUPS)*2*(SIDE+NODES)*MSG),
 ("leaf work: ygen, C, T, S, fwd, inv, bwd, keys w/o emits, set-up (E)",_PS(LEAF_G)),
 ("fixups: 6 per expected fixup (E)",E_FIXUP),
 ("key paths: s_g (PJ + 1/3) + DEFJ (E)",_PS(WGJ)+DEFJ),
 ("route work: emits, arrival, row blocks, fin, tables (E)",_PS(ROUTE_G)),
 ("group work: join, main-pair walks, stop test (E)",_PS(GROUP_G)),
 ("same-batch walks (E)",SXV),
 ("K2 verifications (E)",EK2*HCHG_VER),
 ("listings (E)",ELIST*HCHG_LS1),
 ("defective runs: DEFV + DEFK2 + DEFL (E)",DEFV+DEFK2+DEFL),
 ("candidates: 1625371 per listing (E)",E_CAND),
 ("requests, each regenerating Y: 2^33 (2^10 + GEN_MUL) (W)",NODES*REQ_CAP*REQ_SERVE),
 ("phase-4 count exchanges (E)",_PS([EXCH_ROW]*G_GROUPS)),
 ("barriers, phase control, output (W)",NBAR_ALLOW*BARRIER+NODES*PHASE_CTRL_G+OUTPUT),
]
def claim6_of(x):
 c=Fraction(math.ceil(lg(x)*10**6),10**6)
 for _ in range(50):
  if pow2_interval(c)[0]>=x and pow2_interval(c-Fraction(1,10**6))[1]<x:
   return c
  c+=Fraction(1,10**6)if pow2_interval(c)[0]<x else-Fraction(1,10**6)
 raise AssertionError("worst claim not certified")
WORST6=claim6_of(R_WORST)
DECLARED_WORST6=Fraction(190446050,10**6)
check("disclosed worst case 2^%d/10^6: least 6-decimal bound >= R_WORST"%(WORST6*10**6),
   pow2_interval(WORST6)[0]>=R_WORST and pow2_interval(WORST6-Fraction(1,10**6))[1]<R_WORST)
check("declared worst case 2^190446050/10^6 equals the least 6-decimal bound (manifest match)",
   WORST6==DECLARED_WORST6 and pow2_interval(DECLARED_WORST6)[0]>=R_WORST
   and pow2_interval(DECLARED_WORST6-Fraction(1,10**6))[1]<R_WORST)
E6=claim6_of(E_TOTAL)
check("E-total display 2^%d/10^6: least 6-decimal bound >= E_TOTAL"%(E6*10**6),
   pow2_interval(E6)[0]>=E_TOTAL and pow2_interval(E6-Fraction(1,10**6))[1]<E_TOTAL)
W_E=E_TOTAL/C
def claim15_of(w):
 getcontext().prec=80
 start=(Decimal(w.numerator).ln()-Decimal(w.denominator).ln())/Decimal(2).ln()*10**15
 c=Fraction(int(start.to_integral_value(rounding=ROUND_CEILING)),10**15)
 for _ in range(400):
  if pow2_interval(c)[0]>=w and pow2_interval(c-Fraction(1,10**15))[1]<w:
   return c
  c+=Fraction(1,10**15)if pow2_interval(c)[0]<w else-Fraction(1,10**15)
 raise AssertionError("E-claim not certified")
CLAIM_E=claim15_of(W_E)
check("E-work below the worst case",E_TOTAL<R_WORST,
   "E 2^%.6f worst declared 2^%.6f"%(float(E6),float(DECLARED_WORST6)))
check("E-claim = least 15-decimal c with 2^c >= E/1000",pow2_interval(CLAIM_E)[0]>=W_E
   and pow2_interval(CLAIM_E-Fraction(1,10**15))[1]<W_E,"%d/10^15"%(CLAIM_E*10**15))
JROW_E=sum(PHI[_g]*WGJ[_g]for _g in range(G_GROUPS))+DEFJ
check("key paths E-priced: PJ = 2/3 + 2 d_wj + 1/(r - 1) per leaf + a B per sign, per group under PHI, plus DEFJ",
   PJ>Fraction(2,3)and D_WJ>Fraction(1,r-1)and JROW_E<JROW_W
   and all(WGJ[_g]==JPATH*SG[_g]*PJ+SG[_g]//3*SJOIN for _g in range(G_GROUPS))and DEFJ==JPATH*D_BATCH*SG_FULL,
   "2^%.4f, every-leaf row 2^%.4f"%(lg(JROW_E),lg(JROW_W)))
check("worst case charges every leaf's key path and sign's B: R_WORST = GATE + E_WORST_REST + JROW_W + sum WGW",
   R_WORST-GATE-JROW_W==E_WORST_REST+sum(WGW)and JROW_W==NODES*N0G*JPATH+NODES*N0G//3*SJOIN)
check("K2 verifications: main pairs at KB(136) + 1/r per group under PHI; same-X pairs at their visits (PAIRM); "
   "dependent terms global; below the product row",
   KB136==KB(D1+D2+64)and KB136<KB(D1+D2)*KB64 and PAIRK2==KB136+Fraction(1,r)
   and all(EK2G[_g]==RHO*(SG[_g]*(NODES*LEAVES*MQ1[_g]*Q2)+SG[_g]*(SG[_g]-1)//2-SAMEX[_g])
      *PAIRK2 for _g in range(G_GROUPS))and EK2SX==PAIRM*sum(SAMEX)
   and EK2DEP==RHO*(2*SG_FULL+PARTNER_PAIRS)*KB64+_SIPT*KB64*SG_FULL
   and EK2==sum(PHI[_g]*EK2G[_g]for _g in range(G_GROUPS))+EK2SX+EK2DEP
   and EK2<RHO*(Fraction(SG_FULL)*SG_FULL*KB(D1+D2)*KB64+(2*SG_FULL+PARTNER_PAIRS)*KB64)
      +_SIPT*KB64*SG_FULL,
   "2^%.4f"%lg(EK2*HCHG_VER))
check("barrier allowance covers the grouped run",NBAR_ALLOW>=NBAR_COUNT,"%d of %d"%(NBAR_COUNT,NBAR_ALLOW))
PER_LEAF_E=E_RUN/SG_FULL
SLOTS=SIDE*CAP_P*32
FIXED=[("code, stack",2**20),("barrier words",5*2048),("header, rows of beta",1024),("permanent X's, constants",MX*(64+6*8*32)),
   ("scalars and coins",(MX+1+4*Q1GD)*32),("pair entries",Q2*512),("Y-records, ygen scratch",L_Y*512+787488),
   ("port message blocks",2*64),("LC tables",2*SIDE*8),("RC tables",2*SIDE*8),("record blocks",2*SIDE*2**14),
   ("final-cursor tables E1, E2",3*SIDE*32),("epoch cursors E3",SIDE*8),("hash-read cursors",SIDE*8),("join state, p word",40),
   ("candidate list",CAND_CAP*16),("request and reply buffers",REQ_CAP*64+8*2*CAND_CAP*64),("relay buffer",RELAY_CAP*64),
   ("phase-4 tables",3*SIDE*32)]
FIXED_BYTES=sum(b for _,b in FIXED)
SLOTS_G=2**D1*CAP_PG*32
PHASE3G=[("per-group hop-1 slots",SLOTS_G),("hash heads",2**HB*8),("next-array (next address, K2 per R2 slot)",SIDE*CAP_P*16)]
PHASE3G_BYTES=sum(b for _,b in PHASE3G)
R1_BASE=2**52
HEAD_BASE,NEXT_BASE=2**61,2**62
R1_END=NEXT_BASE+SIDE*CAP_P*16  # R1 sized to its contents
check("R1 = [2^52, 2^63 - 2^57) holds exactly the per-group slots, hash heads and next-array",
   PHASE3G_BYTES==R1_END-R1_BASE==2**63-2**57-2**52,"2^%.4f"%lg(PHASE3G_BYTES))
R2_PAD=32
R2_BASE=R1_END+R2_PAD
check("LS64 alignment: window slot n at a 2^14-aligned base + 32n, R2 slot c + 32(n - 1), c = 32 mod 64: 64-byte aligned "
   "exactly for even n",
   R2_BASE%64==32 and (32*CAP_P)%64==0 and (32*K)%64==0 and BLKC%64==0 and PB2C%64==0
   and all(((BLKC+32*n)%64==0)==((R2_BASE+32*(n-1))%64==0)==(n%2==0)for n in range(1,K+1)))
check("R1 layout: hop-1 slots at 2^52, heads at 2^61 (2^61 bytes, bucket in bits 3..60: one BFI), next-array at 2^62 "
   "to R1's end; disjoint, in R1",R2_BASE%2==0 and R1_BASE+SLOTS_G<=HEAD_BASE and HEAD_BASE+2**HB*8<=NEXT_BASE
   and NEXT_BASE+(SLOTS>>1)==R1_END==R2_BASE-R2_PAD and HEAD_BASE%2**61==0 and 2**HB*8==2**61
   and (SLOTS>>1)==PHASE3G[2][1],"heads end 2^%.4f, next-array end 2^%.4f, R1 end 2^%.4f"
   %(lg(HEAD_BASE+2**HB*8),lg(NEXT_BASE+(SLOTS>>1)),lg(R2_BASE-R2_PAD)))
for name,b in FIXED+[("R1 at 2^52",R1_END-R1_BASE)]+PHASE3G+[("R2 at R1's end + 32",SLOTS)]:
 print("  %-42s 2^%.4f bytes"%(name,lg(b)))
check("fixed regions lie below 2^52",FIXED_BYTES<=2**52,
   "2^%.4f"%lg(FIXED_BYTES))
ALIGNED=4*SIDE*8+2*SIDE*2**14+SIDE*8
check("aligned tables: LC at 0 and 2^39, RC at 2^43, hash-read cursors at 2^48, column blocks at 2^50, row blocks at "
   "2^51; the rest fits between 2^44 and 2^48",
   SIDE*8==2**39 and SIDE*2**14==2**50 and PB1C==1<<50 and PB2C==1<<51
   and PB2C+(1<<50)<=1<<52 and FIXED_BYTES-ALIGNED<=(1<<48)-(1<<44),
   "2^%.4f"%lg(FIXED_BYTES-ALIGNED))
check("record block: header and k slots in 2^14; bit 13 set exactly at slot k; low 13 bits of bases zero; cap_p, group "
   "records multiples of k",
   32+32*K<=2**14 and all((((32*n)>>13)&1==1)==(n==K)for n in range(1,K+1))
   and PB1C&0x1FFF==0 and PB2C&0x1FFF==0 and BLKC&0x1FFF==0 and CAP_P%K==0
   and SG_FULL%K==0 and all(s%K==0 for s in SG))
check("bwd's and keys' loop end: q2 512-byte pair entries fill the lower 2^34 bytes of a 2^35-aligned block",
   Q2*512==2**34)
EPOCH_G=390
EPOCH_B=EPOCH_G*Q1G[0]
def k_decode(kf,late):  # kf = k mod 2^28; late: R2 address >= E3[j]
 return kf+2**28*(late and kf<EPOCH_B)
_KS=sorted({k for c in (0,EPOCH_B,2**28,Q1GD)for k in range(c-300,c+300)if 0<=k<Q1GD}|set(range(0,Q1GD,99991)))
check("record fields: w2 bits 36-63 hold kf = k mod 2^28; k = kf + 2^28 [late and kf < B] exact, B = %d*%d = %d <= 2^28, "
   "q1 - 2^28 < B (%d k checked); w3's code 2 (5 on j = 1) and l < 2^25 in bits 39-63"%(EPOCH_G,Q1G[0],EPOCH_B,len(_KS)),Q1GD<2**29 and EPOCH_B<=2**28<Q1GD<2**28+EPOCH_B
   and all(Q1G[0]*G_GROUPS==Q1GD and k_decode(k%2**28,k>=EPOCH_B)==k for k in _KS)
   and Q2<=2**25 and 128>>6==2!=7-2)
check("R1 and R2 end below the beta tables at 2^64 - 2^51 (2^43 entries of 256 bytes)",
   R2_BASE+SLOTS<=BT_BASE and BT_BASE+(2**BT_W<<8)==2**64 and 6*BT_W>=256>5*BT_W,
   "2^%.6f"%lg(R2_BASE+SLOTS))
CURSOR_END=2**52+32*CAP_P*(SIDE-1)+32*N0G
check("cursors cannot wrap: the last slot's base + 32 per leaf of a node (all groups) < 2^64",CAP_PG<=CAP_P and CURSOR_END<2**64,"2^%.4f"%lg(CURSOR_END))
check("hash-read cursors occupy [2^48, 2^48 + 2^39), below the column blocks at 2^50",(1<<48)+SIDE*8<=PB1C)
TB0,XT_SPAN=2**64-3*2**57+2**49,2**58
_WV=[0,1,2**256-1,p-1,2**255+2**200-1]+[(1<<o+w)-1^(1<<o)-1 for o,w in zip(W7O,W7)]
check("W7 windows %s at bits %s: v = sum of v_k 2^o_k (edge values); lines 0-3 built to 2^37 rows, 4-6 to 2^36: "
   "tblx %d + %d*2^36 a pair table, %d + %d*2^36 a delta table"%(W7,W7O,TBXP%2**36,TBXP>>36,TBXD%2**36,TBXD>>36),
   sum(W7)==256 and all(sum((v>>o&(1<<w)-1)<<o for o,w in zip(W7O,W7))==v for v in _WV)
   and TBXP==TBX(4,2**W7[0])+TBX(3,2**W7[4])and TBXD==TBX(2,2**W7[0])+TBX(2,2**W7[4]))
check("X tables (W7, PX): pair p's 2^49 bytes at TB0 + 2^49 p, 2^46-byte slots PE, PF, PG, sigma, P1_j, P1_j+1 (+512t + 64k), delta "
   "(2^45: +256t + 32k); offsets 512 v + 64k + n 2^46; R2 < TB0, TB0 + 2^58 <= beta; %d + %d tables, 3*2^56 bytes"%(NPAIR_T,NSINGLE_T),
   R2_BASE+SLOTS<=TB0 and TB0+XT_SPAN<=BT_BASE and TB0%2**49==0 and MX//2*2**49==XT_SPAN and 512*(2**37-1)+64*6+64<=2**46
   and 256*(2**37-1)+64*3+64<=2**45 and 6*2**46+2**45<=2**49 and XT_BYTES==3*2**56 and NPAIR_T==3*MX//2+MX//4+MX)
RECEIVE_BYTES=2**16
check("three finite receive latches and processor state fit reserve",3*2**14+4*1024<=RECEIVE_BYTES)
check("main beta backing overlaps only delivery-private aperture",BT_BASE<2**64-2**14 and BT_BASE+(2**BT_W<<8)==2**64)
PEAK=NODES*(RECEIVE_BYTES+FIXED_BYTES+R1_END-R1_BASE+R2_PAD+SLOTS+(2**BT_W<<8)+XT_BYTES)
check("peak memory 2^72 * (fixed + R1 + 32 + R2 + beta tables + 3*2^56 bytes of X tables + receive reserve) = 2^%.4f bytes"%lg(PEAK),PEAK<2**136)
check("caps exceed crude means: cap_p > m_p2 + 2^16 at hop 2, cap_pg > m_p1 + 2^17 per group at hop 1",
   CAP_P>TFULL["m_p2"]+2**16 and CAP_PG>TFULL["m_p1"]+2**17,
   "m_p2 2^%.3f, m_p1 2^%.3f"%(lg(TFULL["m_p2"]),lg(TFULL["m_p1"])))
check("preprocessing E_PRE (the ledger's P rows) = 2^%.7f units, included"
   %lg(Fraction(E_PRE,C)),E_PRE<E_TOTAL<R_WORST)
DECLARED=dict(time_log2=CLAIM_E,worst_log2=DECLARED_WORST6,memory_log2_bytes=None,preprocessing_log2=None,
    success=Fraction(39,100),
    per_leaf=Fraction(cdiv(PER_LEAF_E.numerator*10**4,PER_LEAF_E.denominator),10**4))
def least_at_or_above(d,x,scale):
 return pow2_interval(d)[0]>=x and pow2_interval(d-Fraction(1,scale))[1]<x
def least_grid_up(x,scale):
 k=math.ceil(lg(x)*scale)
 for _ in range(20):
  d=Fraction(k,scale)
  if least_at_or_above(d,x,scale):
   return d
  k+=1 if pow2_interval(d)[0]<x else-1
 raise AssertionError("grid bound not certified")
DECLARED["memory_log2_bytes"]=least_grid_up(PEAK,10**4)
DECLARED["preprocessing_log2"]=least_grid_up(Fraction(E_PRE,C),10**4)
check("declared E time_log2 %d/10^15: least 15-decimal value with 2^c >= E/1000"%(DECLARED["time_log2"]*10**15),
   DECLARED["time_log2"]==CLAIM_E and least_at_or_above(DECLARED["time_log2"],W_E,10**15))
check("declared memory 2^%s bytes: least 4-decimal bound >= peak (2^%.7f)"
   %(format(float(DECLARED["memory_log2_bytes"]),".4f"),lg(PEAK)),
   least_at_or_above(DECLARED["memory_log2_bytes"],PEAK,10**4))
check("declared preprocessing 2^%s units: least 4-decimal bound >= E_PRE / 1000"
   %format(float(DECLARED["preprocessing_log2"]),".4f"),
   least_at_or_above(DECLARED["preprocessing_log2"],Fraction(E_PRE,C),10**4))
check("declared success 0.39 <= (K/2^64) lo_run",DECLARED["success"]<=P_GATE*TFULL["lo"])
STATED=dict(time_log2=Fraction(123439586059213053,10**15),worst_log2=Fraction(190446050,10**6),memory_log2_bytes=Fraction(1359831,10**4),preprocessing_log2=Fraction(1180053,10**4),success=Fraction(39,100),per_leaf=Fraction(723559,10**4))
STATED_K=7203334947470835147
import json,os,re
_CL=json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"submission.json")))["claims"][0]
_ST,_SC=_CL["statement"],_CL["score"]
_PM=open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"proof.md"),encoding="utf-8").read()
_T7=[l for l in _PM[_PM.index("## 7."):].split("\n")if l.startswith("| ")and "|" in l[2:]][2:]
ROWS7=["| %s | %.4f |"%(n,lg(v))for n,v in LEDGER]
check("proof.md 7: each of the %d printed rows equals its log2 here (4 decimals), the rows sum to E_TOTAL, and the "
   "printed total is E6"%len(LEDGER),_T7[:len(LEDGER)]==ROWS7 and sum(v for _,v in LEDGER)==E_RUN
   and _T7[len(LEDGER)]=="| **opened run E_run** | **%.6f** |"%lg(E_RUN)
   and _T7[len(LEDGER)+1]=="| **expected total E = 2^7 + (K/2^64) E_run** | **%s** |"%format(float(E6),".6f"))
E_NOABORT=E_RUN+sum((1-PHI[_g])*(WG[_g]+HCHG_VER*EK2G[_g])for _g in range(G_GROUPS))
check("proof.md 7: a run that never stops early (all PHI = 1) costs 2^%.4f units, as printed"%lg(Fraction(E_NOABORT,C)),
   "costs about 2^%.4f units (the certificate's E with all PHI[g] = 1)"%lg(Fraction(E_NOABORT,C))in _PM
   and "an opened run costs about 2^%.4f units"%lg(Fraction(E_RUN,C))in _PM)
check("stated figures (literals here, submission.json's statement and score) equal the certified ones",
   STATED==DECLARED and Fraction(int(_SC["numerator"]),int(_SC["denominator"]))==STATED["time_log2"]and all(re.search(re.escape(t)+r"(?![0-9])",_ST)for t in ("2^123.439586059213053","2^190.446050","2^135.9831","2^118.0053",
              "72.3559","0.39","K = %d"%STATED_K))and STATED_K==GATE_K
   and "K = %d"%STATED_K in _PM and "K = %d"%STATED_K in open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"success.md"),encoding="utf-8").read())
check("declared E work per leaf %s: least 4-decimal bound >= E/S (%.7f)"
   %(format(float(DECLARED["per_leaf"]),".4f"),float(PER_LEAF_E)),
   DECLARED["per_leaf"]>=PER_LEAF_E>DECLARED["per_leaf"]-Fraction(1,10**4))
MULX,SQRX,ADDX=105+4,80+3,9+2*9
DBL,MADD_=2*MULX+5*SQRX+10*ADDX,8*MULX+3*SQRX+10*ADDX
GEN_ITEM=256*(DBL+MADD_)+2*(INV+3*MULX+2*SQRX+10*ADDX)+64
GEN_TEST=256*2*(9+4)+64
GEN_CTRL=1600
NEGT=30
RED=9
Y_ENTRY=3*SQRX+2*MULX+2*ADDX+2*NEGT+7*RED+14
X_ENTRY=MADD_+3*SQRX+11*MULX+2*ADDX+2*NEGT+17*RED+32
check("generator: doubling 903, mixed addition 1,391; with tests, control and max(Y_ENTRY, X_ENTRY) at most %d, the allowance"%(GEN_ITEM+GEN_TEST+GEN_CTRL+max(Y_ENTRY,
   X_ENTRY)),DBL==903 and MADD_==1391 and GEN_ITEM+GEN_TEST+GEN_CTRL+max(Y_ENTRY,X_ENTRY)<=GEN_MUL,
   "%.4fx"%(math.floor(10**4*GEN_MUL/(GEN_ITEM+GEN_TEST+GEN_CTRL+max(Y_ENTRY,X_ENTRY)))/10**4))
check("permanent X generator: [a]P + R_g with tests, control and four affine-store STPs fits GEN_MULX",
   GEN_ITEM+GEN_TEST+GEN_CTRL+4<=GEN_MULX,"%d <= %d"%(GEN_ITEM+GEN_TEST+GEN_CTRL+4,GEN_MULX))
XSET_ITEM=4*MULX+3*SQRX+10*ADDX+2*NEGT+6*7*MULX+64
check("X setup: 2a^3, 2a, a^2, a^4, y^2, 2y^2, 2y, sigma, pi - 7 sigma, delta, -T(a), -T(beta a), 7 steps a table fit XSETUP",
   XSET_ITEM<=XSETUP,"%d <= %d"%(XSET_ITEM,XSETUP))
ENTRY_ITEM=4+4+2+2+11+6+3
check("pair entry X words (corrections, y^2: 2 LDP + 2 STP each; PE, sigma bases; P1, delta line bases 7 + 4, 4 + 2; loop 3) "
   "fit ENTRY_X",ENTRY_ITEM<=ENTRY_X)
BARRIER_ITEM=5*8+5*(6+MSG_LINK)+1+20
BATCH_ITEM=(2+4)+60+(1+4+4+4+2+1)+(2*4+1)+(4+1+2+1+1+1)+2*(4+1) \
 +2*(2+1+3)+2*(1+4+1+1)+(2+2+1+1)+2
check("batch set-up: base %d plus16 equal-group controller instructions fit176"%BATCH_ITEM,BATCH_ITEM+16<=BATCH_SETUP)
check("phase control: 12 kernel calls a group <= 33 each, stop test and loop 8, h2s's p word 4, over %d groups < %d"%(G_GROUPS,PHASE_CTRL_G),G_GROUPS*(12*33+8+4)<PHASE_CTRL_G)
check("barrier: at most 131 instructions per node with in-flight counts and deliveries (2^10 is > 7x)",
   BARRIER_ITEM==131 and BARRIER>7*NODES*BARRIER_ITEM)
EXCH_ITEM=12*SIDE+4*(14+MSG_LINK)*(8*REQ_CAP)
check("phase-4 count exchange per node and sub-phase: 12 per entry of the 2^36-entry table, 88 per message for 8*2^33",EXCH_ITEM==EXCH_G==12*2**36+88*2**36
   and EXCH_ROW==NODES*4*EXCH_G,"2^%.4f"%lg(EXCH_ITEM))
check("a broadcast message costs <= 16 with in-flight count and delivery, within its price of 40",
   2+3+3+MSG_LINK<=MSG)
CAND_LEAVES=INV+5*SQRX+18*MULX+20*ADDX
EPOCH_DEC=24  # j (SUB, UDIV, 6 MOV), E3[j] (2 MOV, LDR), CMP, CSET, LSR, B (2 MOV), CMP, CSET, AND, ADD, STR
CAND_ITEM=4*8*2*MSG+4*CAND_LEAVES+36*600+384*200+2*256*(DBL+MADD_)//2+INV+2*EPOCH_DEC
check("candidate: messages, both records' decode, six leaves each, 36 unit tests, solve mod r, [x]P = Q = %d; charged "
   "above twice that"%CAND_ITEM,RESOLVE==2*CAND_ITEM+1,
   "%.3fx"%(math.floor(10**3*RESOLVE/CAND_ITEM)/10**3))
PROGRAM_BYTES=2**18
KERNEL_INSTR=22610+14963+1524  # leaf.s, route.s, secp_fe.s (experiment.py checks)
check("program delivery: 2^18 bytes (%d kernel instructions, %d bytes, + <= 2^16 + 64 for generator, resolver, control) "
   "= 2^13 messages over 2 hops at 40 < 2^20"%(KERNEL_INSTR,4*KERNEL_INSTR),
   4*KERNEL_INSTR+2**16+64<=PROGRAM_BYTES and PROGRAM_BYTES//32*2*MSG+128<NODE_INIT)
ANALYTIC_E=NODES*(4*Q1GD+MX+1)*COIN_SCALAR+NODES*MX*(GEN_MULX+XSETUP)+NODES*Q2*ENTRY_X \
 +NODES*HEADER+NODES*NODE_INIT+NODES*Q1GD*3*GEN_MUL+NODES*GEN_MUL \
 +NODES*Q1GD*(BATCH_SETUP+INV_CALL)+(1+G_GROUPS)*2*(SIDE+NODES)*MSG \
 +E_CAND+NODES*REQ_CAP*REQ_SERVE+sum(PHI[_g]*EXCH_ROW for _g in range(G_GROUPS)) \
 +NBAR_ALLOW*BARRIER+NODES*PHASE_CTRL_G+OUTPUT
check("analytic allowances (with C, T, S, request regeneration, the exchange rows at the caps) are below 2^-8.6 of E_run",
   ANALYTIC_E<pow2_interval(Fraction(-86,10))[0]*E_RUN,"2^%.3f"%lg(Fraction(ANALYTIC_E,E_RUN)))
print()
print("q1 %d  q2 2^25  n0 %d  S 2^%.6f  mu %.9f  run success >= %.10f"%(Q1GD,N0G,lg(SG_FULL),float(TFULL["mu"]),
                 math.floor(TFULL["lo"]*10**10)/10**10))
print("E_total 2^%.6f instructions  W_E 2^%.10f  E-claim %d/10^15  per leaf %.7f declared %s  R_worst 2^%.6f"
   %(float(E6),lg(W_E),CLAIM_E*10**15,float(PER_LEAF_E),format(float(DECLARED["per_leaf"]),".4f"),
   float(DECLARED_WORST6)))
print("key paths: E-priced 2^%.6f (every-leaf row 2^%.6f in the worst case, JPATH %d); Y generation %d per batch (GEN_MUL %d); join "
   "%d fixed, %s per visit"%(lg(JROW_E),lg(JROW_W),JPATH,YGEN,GEN_MUL,HCHG_NEW,VISIT))
print("declared: memory_log2_bytes %s (exact %.7f)  preprocessing_log2 %s (exact %.7f)"%(
 format(float(DECLARED["memory_log2_bytes"]),".4f"),lg(PEAK),
 format(float(DECLARED["preprocessing_log2"]),".4f"),lg(Fraction(E_PRE,C))))
print("FAIL: "+", ".join(FAILED)if FAILED else "PASS")
if __name__=="__main__":
 sys.exit(1 if FAILED else 0)

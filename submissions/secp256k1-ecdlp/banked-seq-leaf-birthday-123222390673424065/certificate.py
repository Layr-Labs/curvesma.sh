import functools,itertools,json,math,os,re,sys;from decimal import Decimal,ROUND_CEILING,getcontext;from fractions import Fraction as Fr;FAILED=[];rn,en,dc=range,enumerate,dict
def check(I,ok,o=""):
 print("%s %s%s"%(I,"ok"if ok else"FAIL",("  "+o)if o else""))
 if not ok:FAILED.append(I)
def lg(x):f=Fr(x);return math.log2(f.numerator)-math.log2(f.denominator)
def cd(a,b):return-((-a)//b)
def expi(x,terms=60):
 assert abs(x)<=1;n,d,S=x.numerator,x.denominator,1<<400;lo=hi=0;I=J=S
 for k in rn(terms+1):
  if k:a,b=sorted((I*n,J*n));I,J=a//(d*k),cd(b,d*k)
  lo,hi=lo+I,hi+J
 o=cd(max(-I,J)*abs(n)*3,d*(terms+1));return Fr(lo-o,S),Fr(hi+o,S)
def ln2_interval(bits=300):
 z,J=Fr(1,3),Fr(1,9);I,pw,k=Fr(0),z,0
 while True:
  I+=pw/(2*k+1);o=pw*J/((2*k+3)*(1-J))
  if o<Fr(1,2**bits):return 2*I,2*(I+o)
  pw*=J;k+=1
LN2=ln2_interval()
@functools.lru_cache(maxsize=None)
def P2(c):n=c.numerator//c.denominator;f=c-n;lo,hi=expi(f*LN2[0])[0],expi(f*LN2[1])[1];o=Fr(2)**n;return lo*o,hi*o
def Z2(n,d):return P2(Fr(n,d))[0]
def enh(o):
 if not o:return Fr(1)
 J,mv=0,o
 while mv>1:mv/=2;J+=1
 mlo=Fr((mv.numerator<<200)//mv.denominator,2**200);I=expi(-mlo,45)[1]
 for _ in rn(J):I=I*I;I=Fr(cd(I.numerator<<300,I.denominator),1<<300)
 return I
p=2**256-2**32-977;r=int("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141",16);C0=2**32+977;BETA_P=int("7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee",16);LAM_R=int("5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72",16);C=1000;NGEN=2**66;M_G=1024;N_CLUST=NGEN//M_G;D1,D2,D3,D4=41,41,64,23;BANK_K2,SLOT_BITS,SLOT_M=44,47,34;N_HBANK=1<<(D1+BANK_K2);HA_SP,HB_SP=44,47;MX=6;YPC=33;YROUNDS=63520;L_CHAIN=YROUNDS*YPC;Q2=MX*L_CHAIN;KEYBITS=64+D4;PART=64+D1;KBF=700;HB=50;THRESH=Fr(39,100);LEAVES=6;FWD_INC,FWD_FIX,FWD_REM=Fr(1494,24),7,(8,381,753,1128)
def fwd_count(q2):assert q2%6==0;m,c=divmod(q2//6,4);return FWD_FIX+FWD_REM[c]+1494*m
BWD_INC,BWD_FIX,BWD_ODD=Fr(10226,12),9,5115
def bwd_count(q2):assert q2%6==0;n=q2//6;return BWD_FIX+10226*(n//2)+BWD_ODD*(n%2)
CH,B_CAP,CHK=2**22,2**26-1,2;FB_FIX=max(fwd_count(6*c)-6*c*FWD_INC for c in rn(4))+BWD_FIX+max(0,BWD_ODD-6*BWD_INC)
FIX_SITES_QUAD=dc(SHX=18,FDX=12,F5=4,FS=6,FADDW=4,BLK=10,DF=2,RW=8);FIX_EXTRA=dc(SHX=5,FDX=5,F5=5,FS=6,FADDW=6,BLK=8,DF=8,RW=5);INS=11;INS_INC=LEAVES*INS;KEY_W3=Fr(1,1);KFC=10;KEY_INC,KEY_FIX=Fr(945,6),5;XEMIT_HOT,XEMIT_FIX=17,3;SELECT=3;OCC_FIX,OCC_WALK,OCC_VISIT,OCC_VER,OCC_LS1,OCC_LS0,OCC_DROP=23,4,8,3,10,7,8;OCC_PAIR=23;OCC_M=21;YSTORE,YSTORE_FIX=41,1;BSEL=SELECT+1;HZ_BANK=SELECT+6;JPATH,JPATH2,SJOIN=3,2,0;JX=((0,JPATH-1,JPATH2-1),)*2+((1,JPATH,JPATH2-1),);INV=22549;TINIT_ENTRY=3
check("c1",OCC_FIX==(9+4+21)-INS and OCC_M==21 and OCC_PAIR>=max(OCC_FIX,OCC_M)and OCC_WALK+OCC_VISIT<=OCC_PAIR and OCC_VISIT==6+2 and OCC_VER==3 and OCC_LS1==max(4+6+2-2,1+6+2-2-2)and OCC_LS0==OCC_LS1-3 and OCC_DROP==12-4);check("c2",KEY_INC==8+48+KEY_W3+6*(4+INS)+KFC+Fr(1,2)and KEY_FIX==5);check("c3",all(fwd_count(6*n)+bwd_count(6*n)<=FB_FIX+6*n*(FWD_INC+BWD_INC)for n in rn(1,9))and CH*512==2**31 and(B_CAP+1)*512==2**35);check("c4",CHK>=1+Fr(L_CHAIN,CH)and 6*CH*512<=2**37);check("c5",(3+3+INS+1)-(4+INS)==JPATH and(4+1+INS+1)-(4+INS)==JPATH2 and SJOIN==0 and JX==((0,2,1),(0,2,1),(1,3,1)));check("c6",SELECT==3 and INS==1+2+SELECT+1+2+2 and XEMIT_HOT==2+3+SELECT+1+2+2+4 and XEMIT_FIX==3);check("c7",YSTORE==2+SELECT+1+7*4+1+3+1+2 and YSTORE_FIX==1 and BSEL==SELECT+1 and HZ_BANK==SELECT+2+4);GEN_MUL=645410;GEN_MULX=GEN_MUL+4;YR_FIX,YR_BLK,YR_BLK_MINE,YB=22556,62398,40552,256;CHAIN_ADD=2**15;CHAIN_SET=2*GEN_MUL+17*CHAIN_ADD;Y_PROD_FIX=64
def y_group(q):o=cd(q,M_G);J=YROUNDS*(cd(o,YB)*YR_FIX+o*YR_BLK)+o*CHAIN_SET;return J+o*(L_CHAIN*YSTORE+YSTORE_FIX+Y_PROD_FIX)
NCAND=dc(A=8,B=6,EB=7);EISR={k:357+31*v for k,v in NCAND.items()};KCL=dc(E=("A",)*6,F=("A",)*3+("B",)*3,G=("A",)*3+("B",)*3)
def tbldk_count(n,cl):return 5+n*(197+6*EISR[cl])
def tbldf_count(n,v):return 4+n*(6+46*v)
XT_SLOTS=6;XT_W,XT_WHI=43,42;XTBLD_PAIRED=6*(4+2*11)+20*2*(4*2**XT_W+2*2**XT_WHI);XTBLD_ALL=6*XTBLD_PAIRED+2*(5*26+200*2**XT_W)+sum(tbldk_count(2**(XT_WHI if k in(2,5)and c!="F"else XT_W),KCL[c][k])for c in"EFG"for k in rn(5 if c=="F"else 6));XSETUP=2**17;ENTRY_X=16;COIN_SCALAR=64;BATCH_SETUP=2048;INV_CALL=16;MACH_INIT=2**31;GEN_INIT=2**12;YQ_INIT=64;TB_BANK=64;HZ_D1=16;GROUP_CTRL=2**10;CHUNK_SAVE=16;XDIST=64;REPORT=2**10;HEADER=80;BT_BASE,BT_W=0xfff8<<48,43;TBLD=tbldf_count(2**BT_W,8);CAND_CTL=2**13;LIST_CAP=36*CH;MULX,SQRX,ADDX=105+4,80+3,9+2*9;DBL,MADD_=2*MULX+5*SQRX+10*ADDX,8*MULX+3*SQRX+10*ADDX;CAND_LEAVES=6*CHAIN_ADD+12*MULX;CAND_ITEM=CAND_CTL+2*CAND_LEAVES+36*600+2*(384*200+256*(DBL+MADD_)+INV);RESOLVE=2*CAND_ITEM+1;REQ_SERVE=2**10+GEN_MUL;OUTPUT=2**20;DELTA_F=2**63-2**34;QP,QQ=Q2//2,Q2//4;FIX_BATCH=sum(FIX_SITES_QUAD.values())*CH*6//4;P1=Fr(6*(r-1)+1,r*r);TWB=Fr(2,r-1);SHX_BAD=2**64*2**33;FOLD_BAD=2**66*C0;SQ_BAD=C0*C0+2*C0;FS_P1_BAD=C0;SITES=[(3*QP,SHX_BAD,36),(QP,FOLD_BAD,36),(3*QP,SHX_BAD,36),(QP,FOLD_BAD,36),(QP,FOLD_BAD,36),(3*QP,SHX_BAD,36),(QP,FOLD_BAD,36),(Q2,SQ_BAD,12),(Q2,FOLD_BAD,12),(2*Q2,FS_P1_BAD,6)];SITES_QUAD=dc(SHX=18,FDX=12,F5=4,FS_FA=8);D_SITES=sum(Fr(n*bad*mult,r-1)for n,bad,mult in SITES);D_YGEN=Fr(YROUNDS*(YPC+2*16+2),r-1);D_BATCH=Fr(CH,L_CHAIN)*(D_SITES+Q2*Fr(6,r-1))+CH*D_YGEN;FIX_ROW=C0*2**192;FIX_FOLD=(2**66+1)*2**128*C0;FIX_SQ=2**193;FIX_LK6=15*C0*2**192;FIX_LK7=21*C0*2**192;FIX_SUB=C0*2**192;FIX_DF63=2**128*(2**63+2)*C0;FIX_RW=C0*2**128;FIXUPS=[("BLK",QP,FIX_LK7,36),("SHX",3*QP,FIX_ROW,36),("FDX",QP,FIX_FOLD,36),("SHX",3*QP,FIX_ROW,36),("FDX",QP,FIX_FOLD,36),("FDX",QP,FIX_FOLD,36),("FS",QP,FIX_SUB,36),("SHX",3*QP,FIX_ROW,36),("FDX",QP,FIX_FOLD,36),("DF",QP,FIX_DF63,6),("F5",Q2,FIX_SQ,12),("BLK",Q2,FIX_LK7,12),("BLK",Q2,FIX_LK6,12),("FDX",Q2,FIX_FOLD,12),("FS",Q2,FIX_SUB,6),("FADDW",Q2,FIX_SUB,6),("RW",2*Q2,FIX_RW,6)];E_FIX=sum(Fr(n*bad*mult,r-1)for _,n,bad,mult in FIXUPS)+CHK*D_BATCH*FIX_BATCH;KB=functools.lru_cache(None)(lambda fixed:Fr(6*2**(256-fixed)*(KBF+1),(r-1)*KBF));RHO=Fr(r,r-1)**2;SLF=6*RHO*Fr(2**256,r);SL6=6*L_CHAIN*Fr(2**256,r);D_DROP=2*SLF*L_CHAIN/(B_CAP-2*L_CHAIN-SL6)**2;CLUST_BATCH=M_G*LEAVES*Q2;FIXED_LEAVES=2*CLUST_BATCH;PARTNER_PAIRS=CLUST_BATCH*(CLUST_BATCH-1);FA,FB=int("114ca50f7a8e2f3f657c1108d9d44cfd9",16),int("e4437ed6010e88286f547fa90abfe4c3",16);MH_A,MH_B,ML=20,21,23
def isum(*rs):return sum(x[0]for x in rs),sum(x[1]for x in rs)
def _cvg():
 ds,u,v=[],BETA_P,p
 while u:ds+=[v//u];u,v=v%u,u
 o,qm,qk,pm,pk=[],0,1,1,0
 for d in ds:qn,pn=d*qk+qm,d*pk+pm;o+=[(qn,qn*BETA_P-pn*p,d)];qm,qk,pm,pk=qk,qn,pk,pn
 return o
CVG=_cvg()
def eis_cls(X2):
 for k in rn(1,len(CVG)-1):
  (qk,dk,_),(qm,dm,_),(qn,_,a)=CVG[k],CVG[k-1],CVG[k+1]
  if qk+qm<=X2-1<qn+qk:g=max(abs(dk),abs(dm)-(min((X2-1-qm)//qk,a)-1)*abs(dk))
 bb=None
 for k in rn(1,len(CVG)):
  (qk,dk,_),(qm,dm,_)=CVG[k],CVG[k-1]
  if qk+qm>X2:break
  if(dk>0)==(dm>0):continue
  for j in rn((X2-qk-qm)//qk+1):
   dw,qw=dm+j*dk,qm+j*qk
   if(dw>0)==(dk>0)and dw:break
   if bb is None or abs(dk)+abs(dw)<bb[0]:bb=(abs(dk)+abs(dw),(-dk,qk),(-dw,qw))
 _,U,W=bb
 if U[0]*W[1]<U[1]*W[0]:U,W=W,U
 lo,hi=min(0,U[0])+min(0,W[0]),max(0,U[0])+max(0,W[0]);o=(-(lo+hi)//2,(X2-U[1]-W[1])//2);C=(g+1)//2
 cs=sorted(((i*U[0]+j*W[0],i*U[1]+j*W[1])for i in rn(-4,5)for j in rn(-4,5)if o[0]+lo+i*U[0]+j*W[0]<=C and o[0]+hi+i*U[0]+j*W[0]>=-C and o[1]+i*U[1]+j*W[1]<X2 and o[1]+U[1]+W[1]+i*U[1]+j*W[1]>0),key=lambda t:(abs(t[0])+64*abs(t[1]),t))
 return dc(X2=X2,C=C,U=U,W=W,o=o,c=cs,kS=(W[0]*o[1]-W[1]*o[0])%p,kT=(U[1]*o[0]-U[0]*o[1])%p)
def eis_red(c,v):
 U,W,o=c["U"],c["W"],c["o"];S=(W[1]*v+c["kS"])%p;T=(c["kT"]-U[1]*v)%p;f1,f2=o[0]+(S*U[0]+T*W[0])//p,o[1]+(S*U[1]+T*W[1])//p
 for l1,l2 in c["c"]:
  if 0<=f2+l2<c["X2"]and abs(f1+l1)<=c["C"]:return f1+l1,f2+l2
 return f1,f2
def eis_cover(c):
 U,W,(o1,o2),X2,C=c["U"],c["W"],c["o"],c["X2"],c["C"];V=[(o1,o2),(o1+U[0],o2+U[1]),(o1+U[0]+W[0],o2+U[1]+W[1]),(o1+W[0],o2+W[1])];E_=list(zip(V,V[1:]+V[:1]))
 ys=sorted({y for y in[v[1]for v in V]+[e for _,l in c["c"]for e in(-l,X2-l)]if o2<=y<=o2+U[1]+W[1]})
 def xs(y):z=[x0+Fr(y-y0)*(x1-x0)/(y1-y0)for(x0,y0),(x1,y1)in E_ if y0!=y1 and min(y0,y1)<=y<=max(y0,y1)];return min(z),max(z)
 for ya,yb in zip(ys,ys[1:]):
  cm=[]
  for a,b in sorted((-C-l1,C-l1)for l1,l2 in c["c"]if-l2<=ya and yb<=X2-l2):
   if cm and a<=cm[-1][1]+1:cm[-1]=(cm[-1][0],max(cm[-1][1],b))
   else:cm+=[(a,b)]
  (La,Ra),(Lb,Rb)=xs(ya),xs(yb)
  if not any(a<=min(La,Lb)and max(Ra,Rb)<=b for a,b in cm):return False
 return True
X2A,X2B=45402545270420317080910594221470922860,64502973549206556628585045361533709079;EIS=dc(A=eis_cls(X2A),B=eis_cls(X2B),EB=eis_cls(2**128-6*(X2A-1)));C_E,C_FG=sum(EIS[c]["C"]for c in KCL["E"])+EIS["EB"]["C"],sum(EIS[c]["C"]for c in KCL["F"]);E0R,E1R=(0,1+Fr(1,2**30)),(Fr(-3),Fr(3)+Fr(1,2**180))
def _half(lo,hi,ti,e):t0,t1=(lo>>128,hi>>128)if ti else(0,0);s=[m*x/2**128 for m in(lo,hi)for x in e];return(min(t0,0)-1,max(t1-1,0)+5+ti),(math.floor(min(s)),2+math.ceil(max(s)))
def _kt(C):es=(E0R[0]+E1R[0],E0R[1]+E1R[1]);p1,p2,p3=_half(-C,C,1,E0R),_half(0,2**128-1,0,E1R),_half(-C-2**128,C,1,es);return(isum(p1[0],p2[0]),isum(p1[1],p2[1])),((p3[0][0]-p1[0][1],p3[0][1]-p1[0][0]),(p3[1][0]-p1[1][1],p3[1][1]-p1[1][0]))
KT=dc(E=_kt(C_E),F=_kt(C_FG),G=_kt(C_FG));KB2=[[tuple(isum(KT[c][r][w],(-1,0))if c=="E"else KT[c][r][w]for w in(0,1))for r in(0,1)]for c in"EFG"]
def _table3(rs):
 rot=lambda a,b:(a,b,(-2-a[1]-b[1],-a[0]-b[0]));(e0,e1),(f0,f1),(g0,g1)=rs;al=rot(isum(e0,(-1,0)),isum(e1,(-1,0)));bp=rot(*[(-f[1]-2,-f[0]-1)for f in(f0,f1)]);ga=rot(g0,g1)
 W=lambda m,t:isum(al[m%3],bp[(m+t)%3],ga[(m+2*t+1)%3],(0,2));wr=[[W(o+i,t)for t,o in en((2,2,0))]for i in(0,1)];nr=[(-wr[0][t][1]-wr[1][t][1]-1,-wr[0][t][0]-wr[1][t][0]+1)for t in rn(3)];return[(min(x[0]for x in r),max(x[1]for x in r))for r in(wr[0],wr[1],nr)]
ERRH,ERRL=_table3([(KT[c][0][0],KT[c][1][0])for c in"EFG"]),_table3([(KT[c][0][1],KT[c][1][1])for c in"EFG"]);NWIN_K,ERR_K,WW_K,DSYM_K=5748,([(-241,237),(-136,133),(-361,360)],[(-399,401),(-234,228),(-592,588)]),(1468,2442,722,1181),-765;WDH,WDL=sum(hi-lo for lo,hi in ERRH),sum(hi-lo for lo,hi in ERRL);WH,WL=max(hi-lo for lo,hi in ERRH)+1,max(hi-lo for lo,hi in ERRL)+1;(ELO,EHI)=(min(ERRH[0][0],ERRH[1][0]),max(ERRH[0][1],ERRH[1][1]));NWIN=2*(EHI-ELO)+(1-2*ELO)+(3+2*EHI)+2*(EHI+1-ELO)+2*(3*EHI+3-3*ELO);D_WJ=TWB*NWIN*2**192;D_LAB=TWB*(WDH*2**(256-MH_A)+WDL*2**(256-ML));D_LAB_SYM=TWB*3*((WH-1)*2**(256-MH_A)+(WL-1)*2**(256-ML));B_HA,B_HB,B_L=2**MH_A+WH-1,2**MH_B+WH-1,2**ML+WL-1;check("c9",FA*FA-FA*FB+FB*FB==p and(FB*BETA_P-(FA-FB))%p==0 and((FA-FB)*BETA_P+FA)%p==0 and B_HA*B_HB*(B_L*2**128+2*FA)*(B_L*2**128+4*FA-2*FB)<p*2**(MH_A+MH_B+2*ML)*(1+Fr(1,KBF))and(2*B_HA+1)*(2*B_HB+1)*(2*B_L*2**128+FA)*(2*B_L*2**128+2*FA-FB)<16*p*2**(MH_A+MH_B+2*ML)*(1+Fr(1,KBF))and(B_L*2**128+2*FA)*(B_L*2**128+4*FA-2*FB)<p*2**(2*ML)*(1+Fr(1,KBF))and B_L*2**233<p*(1+Fr(1,KBF))and(B_HA*2**192+2*FA)*(B_HB*2**192+4*FA-2*FB)<p*2**(128+MH_A+MH_B)*(1+Fr(1,KBF))and 2**(64-HA_SP)*(2**(HA_SP-BANK_K2+D2)+WH-1)*(B_L*2**128+2*FA)*(B_L*2**128+4*FA-2*FB)<p*2**(64-BANK_K2+D2+2*ML)*(1+Fr(1,KBF))and D1+D2+(HA_SP-MH_A)+(HB_SP-MH_B)==D1+BANK_K2+SLOT_BITS==D1+D2+HB and HA_SP==BANK_K2 and HB_SP==SLOT_BITS);check("c10",NCAND=={k:len(c["c"])for k,c in EIS.items()}and all(eis_cover(c)and all((u[0]+u[1]*BETA_P)%p==0 for u in(c["U"],c["W"]))and c["U"][0]*c["W"][1]-c["U"][1]*c["W"][0]==p for c in EIS.values())and all(sum(EIS[c]["X2"]-1 for c in KCL[k]+(("EB",)if k=="E"else()))<2**128 for k in KCL)and NWIN==NWIN_K,"C_E 2^%.4f C_FG 2^%.4f"%(lg(C_E),lg(C_FG)));check("c11",(ERRH,ERRL)==ERR_K and(WDH,WDL,WH,WL)==WW_K);check("c12",D_WJ<Z2(-5051,100),"2^%.4f"%lg(D_WJ));check("c13",all(sum(n for n,b,_ in SITES if b==B)==SITES_QUAD[k]*QQ for k,B in(("SHX",SHX_BAD),("FDX",FOLD_BAD),("F5",SQ_BAD),("FS_FA",FS_P1_BAD)))and len(SITES)==10);FIX_RUNS={m:sum(n for ms,n,_,_ in FIXUPS if ms==m)for m in FIX_SITES_QUAD};check("c14",FIX_RUNS=={m:FIX_SITES_QUAD[m]*QQ for m in FIX_SITES_QUAD}and sum(FIX_SITES_QUAD.values())==64);check("c15",FIX_BATCH*C0<2**63 and FIX_BATCH*C0>>32==FIX_BATCH+(FIX_BATCH*977>>32),"2^%.3f"%lg(FIX_BATCH*C0));check("c16",D_LAB<Z2(-820,100)and 3*D_WJ+D_LAB_SYM<P2(Fr(DSYM_K,100))[0],"2^%.3f 2^%.4f"%(lg(D_LAB),lg(3*D_WJ+D_LAB_SYM)));check("c17",D_BATCH<Z2(-1258,10)and D_DROP<Z2(-26,1),"2^%.3f"%lg(D_BATCH));G_GROUPS=6656;G_EXP=21
def group_sizes(q1):
 ms=[0];o=G_GROUPS**G_EXP;J=q1**32
 for g in rn(1,G_GROUPS):
  I=J*(g**G_EXP);m=int(q1*((g/G_GROUPS)**(G_EXP/32.0)))+2
  while(m**32)*o>I:m-=1
  while((m+1)**32)*o<=I:m+=1
  ms+=[m]
 ms+=[q1];qs=[ms[_g+1]-ms[_g]for _g in rn(G_GROUPS)];return qs,ms
ZW=16384;ZERO_PASS=2*(ZW//256)+2;STOP_CHECK=5;KB_BUCKET=KB(D1+D2+HB);PAD_C=max(3*EHI+3,EHI+1-ELO,1-3*ELO);WRAPT_C=483
def _npad(i,m,d=PAD_C):
 if not 0<=i<m:return 0
 lo=max(0,i-d);A,B=(lo,min(m-1,m+d-2*i)),(max(lo,m-d-i),min(m-1,2*m+d-2*i));return sum(max(0,y-x+1)for x,y in(A,B))-(max(0,min(A[1],B[1])-max(A[0],B[0])+1)if A[0]<=A[1]and B[0]<=B[1]else 0)
_LC=[(0,0),(1,-PAD_C),(0,2**64-1),(-2,2**64+PAD_C),(-1,2**64-PAD_C),(-2,2*2**64+PAD_C)];co=sorted({0,2**64}|{t+e for a,b in _LC for cg,g in _LC if a!=cg for t in[(g-b)//(a-cg)]for e in rn(-3,4)if 0<=t+e<=2**64});NC_C=sum((_npad(s,2**64)+_npad(e-1,2**64))*(e-s)//2 for s,e in zip(co,co[1:]));NRB,NCB=1<<(64-HA_SP),1<<(64-HB_SP);ROWMAX132=cd(NC_C,1<<HA_SP)+(len(co)-1)*2**64;CNT132=cd(B_HA*B_HB*(ROWMAX132+NRB*3*B_HA),1<<HB_SP)+NRB*4*B_HA*B_HB+2*WRAPT_C*NCB*B_HB+NRB*B_HA*(6*WRAPT_C+8);DPAR=max(FA,2*FA-FB);KB_C132=Fr(6*CNT132*((B_L*(1<<128)+2*DPAR)**2),p*(r-1));check("c19c",WRAPT_C>=max(EHI,1-ELO,1-2*ELO,2+2*EHI)and all((_npad(e-1,2**64)-_npad(s,2**64))*((s+e-1)//2-s)==(_npad((s+e-1)//2,2**64)-_npad(s,2**64))*(e-1-s)for s,e in zip(co,co[1:]))and max(sum(_npad(h*16+r_,256,12)for h in rn(16))for r_ in rn(16))<=sum(_npad(i,256,12)for i in rn(256))//16+6*256 and KB_C132<Fr(34,100)*KB_BUCKET);KB64=KB(64);OVF_CAP,OVF_CAP_M=2**SLOT_BITS-1,2**SLOT_M-1;P1_GOOD=Fr(6*(r-1),r*r)*(1-3*D_WJ-D_LAB_SYM);K_U_CLUST=LEAVES*M_G*MX;W_ELL=1+47*(3*D_WJ+D_LAB_SYM);SAMEX_FAC=Fr(CLUST_BATCH-K_U_CLUST,CLUST_BATCH-1)*Fr(r,r-2)*W_ELL+Fr(96*(K_U_CLUST-1),CLUST_BATCH-1);EF_Q=E_FIX*C0/(2**32*(DELTA_F-33*FIX_BATCH));KB_M=KB(KEYBITS+D1+D2);SIP=240*RHO*KB_M;LOG2_POOLS=14;N_POOLS=1<<LOG2_POOLS;R_POOL=N_POOLS//M_G;K_GEN=NGEN//N_POOLS
def gf1024_mul(a,b):
 o=0
 for I in rn(10):
  if b>>I&1:o^=a<<I
 for I in rn(18,9,-1):
  if o>>I&1:o^=0x409<<(I-10)
 return o
def pool_of_gen(g):c,m=divmod(g,M_G);return R_POOL*m+(((c&15)^gf1024_mul((c>>4)&0x3FF,m))&15)
def _vand_gf1024_ok():
 b=2
 for _ in rn(10):b=gf1024_mul(b,b)
 if b!=2:return False
 for _d in(2,5):
  J=2
  for _ in rn(_d):J=gf1024_mul(J,J)
  if J==2 or any(gf1024_mul(J^2,z)==0 for z in rn(1,1024)):return False
 I=0
 for d in rn(1024):
  ws=[gf1024_mul(1<<j,d)&15 for j in rn(10)];o=[0]
  for w in ws:o+=[x^w for x in o]
  if d and len(set(o))!=16:return False
  cnt=[0]*16
  for v in o:cnt[v]+=1
  for du in rn(16):ks=cnt[du];ns=1024*1024+2*ks*max(0,ks-1);I+=ns*((1<<42)-(d==0 and du==0))
 return I*2304<9216*1024*1024*((1<<56)-1)
def _bx_tuples_ok():
 if(LAM_R*LAM_R+LAM_R+1)%r!=0:return False
 o=((1,0),(-1,0),(0,1),(0,-1),(-1,-1),(1,1));norms,n2,counts=set(),set(),{}
 for u0 in o:
  n2.add(1)
  for u1 in o:
   a,b=u0[0]+u1[0],u0[1]+u1[1]
   if a or b:n2.add(a*a-a*b+b*b)
 I={}
 for k in(1,2,3,4):
  seen=set()
  for asg in itertools.product(rn(k),repeat=4):
   if len(set(asg))!=k:continue
   for u0 in o:
    for u1 in o:
     for u2 in o:
      sums=[[0,0]for _ in rn(k)]
      for pos,u in en((u0,u1,u2,(1,0))):sums[asg[pos]][0]+=u[0];sums[asg[pos]][1]+=u[1]
      if any(a==0 and b==0 for a,b in sums):continue
      for a,b in sums:
       norms.add(a*a-a*b+b*b)
       if not(a+b*LAM_R)%r:return False
      den=(sums[-1][0]+sums[-1][1]*LAM_R)%r;J=I.get(den)or I.setdefault(den,pow(den,-1,r));seen.add(tuple(((a+b*LAM_R)*J)%r for a,b in sums[:-1]))
  counts[k]=len(seen)
 pm={k:sum(1 for p_ in itertools.product(rn(k),repeat=k)for mb in itertools.product((0,1),repeat=k)if len(set(p_))>=2 and list(dc.fromkeys(p_))==list(rn(len(set(p_))))and len(set(zip(p_,mb)))==k)for k in(2,3,4)};return norms=={1,3,4,7,9,12,13,16}and n2=={1,3,4}and counts=={1:1,2:78,3:252,4:216}and pm=={2:4,3:20,4:76}and max(pm.values())<=96 and 96*(78+252+216)<65536
SAME_M_GEN_PAIRS=M_G*(N_CLUST*N_CLUST*(R_POOL-1)//R_POOL)//2;SHARED_BOTH_CLUST_PAIRS=(N_CLUST*N_CLUST//(R_POOL*R_POOL)-N_CLUST)//2;ACP=N_CLUST*(N_CLUST-1)//2;GOOD_GEN_PAIRS=SAME_M_GEN_PAIRS+SHARED_BOTH_CLUST_PAIRS*(M_G*(M_G-1)//2)+(ACP-SHARED_BOTH_CLUST_PAIRS)*(M_G*(M_G-1));YR_MINE=992;PTS_MINE_CHAIN=YR_MINE*YPC;NBATCH_MINE=23750000*YB;PTS_PER_GEN_MINE=NBATCH_MINE*PTS_MINE_CHAIN;M_POOL=K_GEN*PTS_PER_GEN_MINE;S_MINE_TOT=NGEN*PTS_PER_GEN_MINE;N_CHAINS_POOL=K_GEN*NBATCH_MINE;PPMC=PTS_MINE_CHAIN*(PTS_MINE_CHAIN-1)//2;P_BUCKET_MAX=Fr(6*(1<<63),r-1);P_SAME_CHAIN_POOL=N_CHAINS_POOL*PPMC*P_BUCKET_MAX;P_CHAIN_COUPLE_POOL=N_CHAINS_POOL*PPMC*(M_POOL*M_POOL+(N_CHAINS_POOL-1)*PPMC)*(P_BUCKET_MAX**2);GAMMA_MINE_CROSS=M_POOL*(M_POOL-1)//2-N_CHAINS_POOL*PPMC;P_BDAY_HIT_MIN=Fr((r-1)*(r-1),(1<<193)*r*r)-Fr(6*(r-1),r*r);MU_POOL_LO=GAMMA_MINE_CROSS*P_BDAY_HIT_MIN;CS_POOL=GAMMA_MINE_CROSS*(2*PTS_MINE_CHAIN*M_POOL)*(2*P_BUCKET_MAX*P_BUCKET_MAX);_P_A102=Fr(6*(1<<102),r-1);_P_A171=Fr(6*(1<<(256-D1-BANK_K2)),r-1);P_CELL154_MASK_POOL=(Fr(M_POOL*(M_POOL-1),2)*M_POOL+2*N_CHAINS_POOL*PPMC*M_POOL)*P_BUCKET_MAX*_P_A102;POIS_LE3=enh(MU_POOL_LO)*sum((MU_POOL_LO**_k)/math.factorial(_k)for _k in rn(4));P_FAIL_POOL=POIS_LE3+CS_POOL+P_CELL154_MASK_POOL+P_SAME_CHAIN_POOL+P_CHAIN_COUPLE_POOL;P_FAIL_ALL_POOLS=N_POOLS*P_FAIL_POOL;LAM_ORD=N_POOLS*(M_POOL*M_POOL)*P_BUCKET_MAX;POINTWISE_X_MAX=M_POOL*M_POOL*P_BUCKET_MAX;P_BAD_X=65536*sum((LAM_ORD**(_d-1))*(2*LAM_ORD/r)for _d in(2,3,4))+65536*sum((LAM_ORD**(_d-1))*N_POOLS*(M_POOL*M_POOL)*Fr(30*(1<<63),r*r)for _d in(2,3));P_X_ORBIT_ALL_POOLS=4*LAM_ORD*(6*LAM_ORD/r);P_MINE_ZERO_ALL=Fr(S_MINE_TOT,r);P_SAME_A_ALL=N_POOLS*Fr(M_POOL*(M_POOL-1)//2*6,r-1);D_XMINE_ALL=NGEN*NBATCH_MINE*Fr(YR_MINE*(YPC+34),r-1);_MINE_Y_MEAN=S_MINE_TOT*_P_A171;_MINE_EXP=D1+BANK_K2+Fr(3,2)*_MINE_Y_MEAN/PTS_MINE_CHAIN-Fr(OVF_CAP_M,PTS_MINE_CHAIN);P_MINE_BANK_CAP=Fr(1,1<<200);ECMX=N_POOLS*Fr(M_POOL*(M_POOL-1),2)*_P_A102;E_CAND_MINE=cd(ECMX.numerator,ECMX.denominator);MINE_SLOT_BITS=D1+BANK_K2+SLOT_M;MINE_K3_BITS=D1+BANK_K2+49;E_VIS_MINE=cd((S_MINE_TOT*(S_MINE_TOT-1)//2)*6*(1<<(256-MINE_SLOT_BITS)),r-1);E_K2_MINE=cd(N_POOLS*(M_POOL*(M_POOL-1)//2)*6*(1<<(256-MINE_K3_BITS)),r-1);P_MINE_CAPS=P_MINE_BANK_CAP;check("c19b",_vand_gf1024_ok()and _bx_tuples_ok()and all(pool_of_gen(1024*(_p&15)+(_p>>4))==_p for _p in rn(N_POOLS))and NBATCH_MINE%YB==0 and LN2[0]>Fr(2,3)and _MINE_EXP<-200 and SLOT_M+5+1<=SLOT_BITS+5 and MU_POOL_LO>31 and P_FAIL_ALL_POOLS<Z2(-19,1)and P_BAD_X<Z2(-140,1)and P_MINE_CAPS<Z2(-19,1),"mu_p %.2f p_fail 2^%.2f p_caps 2^%.2f"%(float(MU_POOL_LO),lg(P_FAIL_ALL_POOLS),lg(P_MINE_CAPS)));GEN_MINE_Y=YR_MINE*((NBATCH_MINE//YB)*YR_FIX+NBATCH_MINE*YR_BLK_MINE)+NBATCH_MINE*(CHAIN_SET+2*COIN_SCALAR+64);GEN_MINE=64+GEN_MINE_Y+PTS_PER_GEN_MINE*XEMIT_HOT;_P_MINE_DEF=D_XMINE_ALL+P_MINE_ZERO_ALL+NGEN*(2*NBATCH_MINE+1)*Fr(1,2**254);W_MINE_DET_MAX=S_MINE_TOT*(OVF_CAP_M+1)*(OCC_PAIR+OCC_VER)+S_MINE_TOT*(OCC_LS1+RESOLVE+2*REQ_SERVE+REPORT);ce=_P_MINE_DEF*W_MINE_DET_MAX;DEF_MINE=cd(ce.numerator,ce.denominator);PSEL_POOL=16+1024*256;ZERO_MINE=N_HBANK*((2**(SLOT_M+5)//(8*ZW))*ZERO_PASS+HZ_BANK)+2**D1*HZ_D1;JOIN_MINE=(ZERO_MINE+OCC_PAIR*E_VIS_MINE+OCC_VER*E_K2_MINE+OCC_LS1*E_CAND_MINE+E_CAND_MINE*(RESOLVE+2*REQ_SERVE+REPORT)+DEF_MINE+N_POOLS*(PSEL_POOL+1024*32+STOP_CHECK+8*REQ_SERVE)+NGEN*MX*XDIST);E_MINE=NGEN*GEN_MINE+JOIN_MINE;ck=P_FAIL_ALL_POOLS+P_SAME_A_ALL+P_BAD_X+D_XMINE_ALL+P_MINE_ZERO_ALL+P_X_ORBIT_ALL_POOLS+P_MINE_CAPS;d5=P1*P1+Fr(36,(r-1)**2);cf=KB(D1+BANK_K2);cz=Fr(K_U_CLUST*L_CHAIN,r-1)*(cf**2);cr=1+Fr(1,r)
def spx(q1p):n0p=LEAVES*q1p*Q2;Sp=NGEN*n0p;gamma=GOOD_GEN_PAIRS*(n0p*n0p);I=gamma*P1;mu=gamma*P1_GOOD;cs=gamma*(2*CLUST_BATCH*Sp)*d5;e_hi=enh(mu);d_f=SLF*q1p*(EF_Q+CHK*D_BATCH);m_p2=Sp*cf+1;v_p2=Sp*(K_U_CLUST*cf+cz)+1;d_p2=v_p2/Fr(OVF_CAP-FIXED_LEAVES-m_p2)**2;lam_int=RHO*(Sp+FIXED_LEAVES)*KB_M;per_leaf=D_BATCH+d_f+d_p2+D_DROP;rows_=I*cr*(Fr(1,r)+lam_int+2*per_leaf);g_def=NGEN*q1p*CHK*D_BATCH;glob=ck+NGEN*(3*q1p+2*NBATCH_MINE+1)*Fr(1,2**254)+g_def;lo=1-e_hi-cs-rows_-glob;return dc(mu=mu,lo=lo,cs=cs,rows=rows_,glob=glob,Sp=Sp,n0p=n0p,d_p2=d_p2,m_p2=m_p2,d_f=d_f,lam_int=lam_int,g_def=g_def)
def gfl(q):return spx(q)["lo"]
Q1GD=103079062500;Q1G,MQ1=group_sizes(Q1GD);N0G=LEAVES*Q1GD*Q2;SG_FULL=NGEN*N0G;LO_RUN=gfl(Q1GD);GB,GATE_INSNS=128,256;THRESH_GATE=THRESH;K_GATE=cd(THRESH_GATE.numerator*LO_RUN.denominator*(1<<GB),THRESH_GATE.denominator*LO_RUN.numerator);P_GATE=Fr(K_GATE,1<<GB);GATED_LO=P_GATE*LO_RUN;check("c20 %d"%Q1GD,LO_RUN>=Fr(9985,10000),"lo_run %.12f"%(math.floor(LO_RUN*10**12)/10**12));check("c21",0<K_GATE<(1<<GB)and GATED_LO>=THRESH_GATE and THRESH_GATE>=THRESH and Fr(K_GATE-1,1<<GB)*LO_RUN<THRESH_GATE and 2*COIN_SCALAR+8+3+4+16<=GATE_INSNS,"p_gate %.9f gated_lo %.12f"%(float(P_GATE),math.floor(GATED_LO*10**12)/10**12));check("c22",all((MQ1[g]**32)*(G_GROUPS**G_EXP)<=(Q1GD**32)*(g**G_EXP)<((MQ1[g]+1)**32)*(G_GROUPS**G_EXP)for g in rn(1,G_GROUPS)));TPRE=[spx(m)for m in MQ1[:-1]]
def ceil_256(x):x=min(Fr(1),max(Fr(0),x));return Fr(cd(x.numerator<<256,x.denominator),1<<256)
PHI=[Fr(1)]+[ceil_256(1-t["lo"])for t in TPRE[1:]];TFULL=spx(Q1GD);check("c23 %.6f"%lg(SG_FULL),SG_FULL==NGEN*6*Q1GD*Q2==TFULL["Sp"]);check("c24 %.9f"%float(TFULL["mu"]),TFULL["mu"]==GOOD_GEN_PAIRS*(N0G*N0G)*P1_GOOD and GATED_LO>=THRESH);check("c25",TFULL["cs"]<Fr(1,2**82),"2^%.4f"%lg(TFULL["cs"]));check("c26",TFULL["lam_int"]<Z2(-376,10),"2^%.3f"%lg(TFULL["lam_int"]));check("c27 %.3f"%lg(Q1GD*E_FIX),Q1GD*E_FIX<Z2(39,1)and DELTA_F>2**20*Q1GD*E_FIX+FIX_BATCH);check("c28",TFULL["d_f"]<Z2(-188,10),"2^%.3f"%lg(TFULL["d_f"]));check("c29",TFULL["d_p2"]<Z2(-295,10),"2^%.3f"%lg(TFULL["d_p2"]));check("c30",OVF_CAP-FIXED_LEAVES-TFULL["m_p2"]>0 and LIST_CAP>=LEAVES*MX*CH and LIST_CAP>=PTS_MINE_CHAIN,"m_p2 2^%.3f"%lg(TFULL["m_p2"]));check("c31",TFULL["g_def"]<Z2(-217,10)and TFULL["glob"]<Z2(-19,1),"2^%.3f"%lg(TFULL["glob"]));check("c32",PHI[0]==1 and all(0<=x<=1 for x in PHI)and all(PHI[g]>=PHI[g+1]for g in rn(G_GROUPS-1))and all(PHI[g]>=min(Fr(1),max(Fr(0),1-TPRE[g]["lo"]))for g in rn(1,G_GROUPS)),"PHI[1] %.4f PHI[16] %.4f PHI[31] %.4f"%(float(PHI[1]),float(PHI[16]),float(PHI[31])));HMX=OVF_CAP+1;check("c33",HMX==2**SLOT_BITS and 32*HMX==2**52,"2^%.4f"%lg(HMX));PAIRM=KB_C132+Fr(1,r);PAIRK2=KB(D1+D2+64)+Fr(1,r);PJ_MEAN=Fr(sum(map(sum,JX)),9)+JPATH*(2*D_WJ+Fr(1,r-1));check("c34",OCC_PAIR==23);SG=[NGEN*LEAVES*q*Q2 for q in Q1G];WG,SAMEX,WGJ,WGK2,WGL,WGI,WGG=[],[],[],[],[],[],[]
for _g,_q in en(Q1G):_s=SG[_g];_leaf=NGEN*_q*(CHK*(FB_FIX+COIN_SCALAR+41+CHUNK_SAVE+INV+INV_CALL+KEY_FIX)+BATCH_SETUP)+NGEN*y_group(_q)+cd(NGEN*_q*(FWD_INC+BWD_INC+KEY_INC-INS_INC)*Q2,1);_ins=_s*INS;_gp=NGEN*LEAVES*MQ1[_g]*Q2;_samex=N_CLUST*_q*CLUST_BATCH*(CLUST_BATCH-1)//2;_main=_s*_gp+_s*(_s-1)//2-_samex;_evis=_main*PAIRM;_gfix=GROUP_CTRL;WGK2+=[OCC_VER*_main*PAIRK2];_group=_gfix+OCC_PAIR*_evis+WGK2[-1];_jkey=_s*PJ_MEAN+_s//3*SJOIN;WGJ+=[_jkey];WGL+=[_leaf];WGI+=[_ins];WGG+=[_gfix+OCC_PAIR*_evis];WG+=[_leaf+_ins+_group+_jkey];SAMEX+=[_samex]
check("c35",sum(Q1G)==Q1GD and min(Q1G)>0);N_TBANK,N_YBANK=N_POOLS<<41,N_CLUST;check("c36",NGEN==2**66 and N_YBANK==2**56 and N_TBANK==2**55 and N_HBANK==2**85 and L_CHAIN<2**21 and N_HBANK<<SLOT_BITS==2**(D1+D2+HB));ZERO_INIT=N_HBANK*((2**(SLOT_BITS+5)//(8*ZW))*ZERO_PASS+HZ_BANK)+2**D1*HZ_D1;check("c37",ZERO_PASS==130 and 2**(SLOT_BITS+5)%(8*ZW)==0 and 2**(SLOT_M+5)%(8*ZW)==0 and(2**64-2**55)%2048==0);E_PRE=MACH_INIT+NGEN*GEN_INIT+N_YBANK*YQ_INIT+E_MINE+NGEN*(3*cd(Q1GD,M_G)+1)*COIN_SCALAR+NGEN*MX*XSETUP+NGEN*6*CH*ENTRY_X+NGEN*HEADER+ZERO_INIT+N_TBANK*(XTBLD_ALL+TB_BANK)+TBLD;_SIPT=240*RHO;PCFA=NGEN*(3*Q1GD+2*NBATCH_MINE+1)*Fr(1,2**254);P_EXC=P_X_ORBIT_ALL_POOLS+P_BAD_X+D_XMINE_ALL+P_MINE_ZERO_ALL+PCFA;W_CORE_DET=sum(WGL[_g]+WGI[_g]+GROUP_CTRL for _g in rn(G_GROUPS))+(SG_FULL//3)*SJOIN;W_DET_MAX=(W_CORE_DET+NGEN*Q1GD*(B_CAP//CH+1)*FIX_BATCH*max(FIX_EXTRA.values())+SG_FULL*JPATH+SG_FULL*HMX*(OCC_PAIR+OCC_VER)+SG_FULL*(OCC_LS1+RESOLVE+2*REQ_SERVE));cj=P_EXC*W_DET_MAX;EXC_WORK=cd(cj.numerator,cj.denominator);c4=2*K_U_CLUST*SG_FULL+SAMEX_FAC*SG_FULL*(CLUST_BATCH-1)+SAMEX_FAC*PARTNER_PAIRS+CLUST_BATCH*CLUST_BATCH+96*K_U_CLUST*K_U_CLUST;ELIST=RHO*(Fr(SG_FULL)*SG_FULL*KB_M+c4*KB_M)+_SIPT*KB_M*SG_FULL+72;EK2=RHO*(c4*KB64)+_SIPT*KB64*SG_FULL+72;PAIR_SAMEX=Fr(CLUST_BATCH-K_U_CLUST,CLUST_BATCH-1)*Fr(r,r-2)*W_ELL*PAIRM+Fr(96*(K_U_CLUST-1),CLUST_BATCH-1)*(KB_BUCKET+Fr(1,r));SXK2=PAIR_SAMEX*sum(SAMEX);SXV=PAIR_SAMEX*OCC_PAIR*sum(SAMEX)+EXC_WORK;BRACKET=RHO*(Fr(SG_FULL)*SG_FULL*KB(D1+D2)+c4)+_SIPT*SG_FULL+72;DEFV=OCC_PAIR*D_BATCH*(HMX*SG_FULL+sum(SAMEX));DEFL_PAIRS=2*D_BATCH*BRACKET+D_BATCH*HMX*SG_FULL;DEFK2=OCC_VER*DEFL_PAIRS;DEFL=OCC_LS1*DEFL_PAIRS;DEFJ=JPATH*D_BATCH*SG_FULL;E_FIX_WORK=NGEN*Q1GD*E_FIX*max(FIX_EXTRA.values());E_CAND_WORK=(ELIST+DEFL_PAIRS)*RESOLVE;E_REQ_WORK=(ELIST+DEFL_PAIRS)*(2*REQ_SERVE);E_SHARED_REST=E_PRE+NGEN*GEN_MUL+OUTPUT;E_RUN=E_SHARED_REST+E_FIX_WORK+E_CAND_WORK+E_REQ_WORK+sum(PHI[_g]*WG[_g]for _g in rn(G_GROUPS))+SXV+EK2*OCC_VER+ELIST*OCC_LS1+DEFV+DEFK2+DEFL+DEFJ+SXK2*OCC_VER;E_TOTAL=GATE_INSNS+P_GATE*E_RUN;W_E=E_TOTAL/C
def claim15_of(w):
 getcontext().prec=80;start=(Decimal(w.numerator).ln()-Decimal(w.denominator).ln())/Decimal(2).ln()*10**15;c=Fr(int(start.to_integral_value(rounding=ROUND_CEILING)),10**15)
 for _ in rn(400):
  if P2(c)[0]>=w and P2(c-Fr(1,10**15))[1]<w:return c
  c+=Fr(1,10**15)if P2(c)[0]<w else-Fr(1,10**15)
 raise AssertionError("E-claim not certified")
CLAIM_E=claim15_of(W_E);JROW_E=sum(PHI[_g]*WGJ[_g]for _g in rn(G_GROUPS))+DEFJ;K2ROW_E=sum(PHI[_g]*WGK2[_g]for _g in rn(G_GROUPS))+SXK2*OCC_VER+EK2*OCC_VER;PER_LEAF_OPEN=E_RUN/SG_FULL;PLE=PER_LEAF_OPEN;LOW,UP_BASE=2**55,2**64-2**55
def valid_va(a):return 0<=a<2**64 and(a<2**55 or a>=2**64-2**55)
TB=2**40;XTK_BYTES=3*10240*TB;AXB=4096*TB;AXYX_HI=16384*TB;FWB=19456*TB;XTB=2**64-32768*TB;AXY_END=24576*TB;XT_END=XTB+XTK_BYTES;CN_BASE=84*TB;CMB=90*TB;H_BASE=24576*TB;OV_BASE=H_BASE|2**52;T34=34*TB;SEL=dc(T=T34+2**32,Y=T34+2*2**32,S=T34+3*2**32,H=T34+4*2**32);_KS=[open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"kernels",f)).read()for f in("leaf.s","secp_fe.s")];check("c38",_KS[0].count("mov x1,#%s\nmov x4,#%s\n"%(hex(SEL["H"]),hex(H_BASE)))==1 and _KS[0].count("orr x6,x4,#%s\n"%hex(OV_BASE-H_BASE-32))==2 and H_BASE%2**53==0 and _KS[0].count("HOCC hocc,%d,%s\n"%(SLOT_BITS,hex(2**64-2**SLOT_BITS)))==1 and _KS[0].count("HOCC hoccm,%d,%s\n"%(SLOT_M,hex(2**64-2**SLOT_M)))==1 and len(re.findall(r"bfi x4,[\\\w]+,#5,#%d\n"%SLOT_BITS,_KS[0]))==3 and _KS[0].count("bfi x4,x25,#5,#%d\n"%SLOT_M)==1 and _KS[0].count("tst x0,#%s\n"%hex(2**64-2**D1))==2 and _KS[0].count("ldr x0,[x4,#24]\ncbnz x0,")==4 and _KS[0].count("dsb sy\n")==2 and "stadd"not in _KS[0]and set(re.findall(r"(?:LK1|LKP6_PAIR) (0x[0-9a-f]+)",_KS[0]))=={hex(FWB),hex(FWB+2560*TB)}and[tuple(int(x,16)for x in m.split(","))for m in re.findall(r"OFFS7 (?:\w+,){10}(0x[0-9a-f,x]+)",_KS[1])]==[tuple(XTB+k*10240*TB+d*2**50 for d in(0,2,8,4,6,9))for k in(0,2)]and[tuple(int(x,0)for x in m.split(","))for m in re.findall(r"\nOFF5 (?:\w+,){10}([\w,]+)",_KS[1])]==[tuple(d<<48 for d in(0,2,4,6,8)),tuple(XTB+10240*TB+d*2**50 for d in(0,2,4,6,8))]and _KS[1].count("stp x15,xzr,[x8,#352]\n")==1 and XTB%2**50==0 and BT_BASE==2**64-2**51 and AXB%(2**49)==0 and AXYX_HI%(2**48)==0 and XTB%(2**48)==0 and valid_va(XT_END-1)and all(valid_va(BT_BASE|((w<<8)%2**64))for w in(0,1,2**43-1,2**64-1)));YRSL=2**31;YR_BASE=32*TB+YRSL;EB_END=42*TB+2**37;EB_START=EB_END-6*CH*512;YQ_BASE=40*TB;YQ_LIM=YQ_BASE+64+512*B_CAP;GS_BASE=36*TB;CL_BASE=35*TB;RP_BASE=37*TB;GS_CONTENT=MX*6*XT_SLOTS*64+MX*64+4096+YB*2304+2**20;NEAR1=(T34+0x10000,T34+2**32);REGISTRY=[("N","N1",NEAR1[0],NEAR1[1],"setup:MACH_INIT psel"),("N","CL",CL_BASE,CL_BASE+32*LIST_CAP+64,"hocc hoccm"),("N","RP",RP_BASE,RP_BASE+N_POOLS*(1024*128+256),"setup:REPORT psel"),("N","BT",BT_BASE,2**64,"tbldf")]+[("D","SEL"+w,SEL[w],SEL[w]+2**16,"MMIO")for w in"TYSH"]+[("T","AXY",AXB,AXY_END,"tbld"),("T","XT6",XTB,XT_END,"tbldk"),("Y","YQ",YQ_BASE-64,YQ_BASE+512*(B_CAP+1),"ystore"),("S","YR",YR_BASE-64,YR_BASE-64+(YB+1)*YRSL,"yround"),("S","GS",GS_BASE,GS_BASE+2**30,"setup:GEN_INIT setup:XSETUP yround"),("S","EB",EB_START,EB_END,"setup:ENTRY_X fwd"),("S","CN",CN_BASE,CN_BASE+64*Q1GD,"setup:COIN_SCALAR"),("S","CNM",CMB,CMB+64*NBATCH_MINE,"setup:COIN_SCALAR"),("H","HP",H_BASE,OV_BASE,"hzero keys xemit"),("H","HO",OV_BASE,H_BASE+2**53,"hzero hocc hoccm")];WBITS=dc(T=(LOG2_POOLS,41),Y=(56,0),S=(64,2),H=(D1,BANK_K2));WINS=[(w,SEL[w],WBITS[w],[(lo,hi)for k_,_,lo,hi,_ in REGISTRY if k_==w])for w in"TYSH"];N_BANK=dc(T=N_TBANK,Y=N_YBANK,S=NGEN,H=N_HBANK);cl=os.path.dirname(os.path.abspath(__file__));d3=set(re.findall(r"\.globl\s+(\w+)","".join(open(os.path.join(cl,"kernels",f)).read()for f in("leaf.s","bank.s","secp_fe.s"))))|{"hocc","hoccm"}
def _mvz(a):return any(a&~(0xFFFF<<q)==0 for q in(0,16,32,48))
def _bld_ok(b):t=b.split();return bool(t)and all(x in d3 or x in("MMIO","analytic")or(x[:6]=="setup:"and isinstance(globals().get(x[6:]),int))for x in t)and(len(t)==1 or not{"MMIO","analytic"}&set(t))
def disjoint_va(I):o=sorted(I);return(all(0<=a<b<=2**64 and valid_va(a)and valid_va(b-1)and((a<LOW)==(b-1<LOW))for a,b in o)and all(o[i][1]<=o[i+1][0]for i in rn(len(o)-1)))
KINDS="NDTYSH";ROLE_B={k:sum(hi-lo for r_,_,lo,hi,_ in REGISTRY if r_==k)for k in KINDS};PG,BLK=2**16,2**28
def _pt():
 iv=sorted((lo//PG*PG,cd(hi,PG)*PG,b=="MMIO",hi-lo)for r_,_,lo,hi,b in REGISTRY);pg={};ok=True;l1=set();l2=set();l3=set();I=[];J=[]
 for lo,hi,m,_ in iv:
  for a in(lo,hi-PG):ok&=pg.setdefault(a,m)==m
  t=J if m else I
  if t and lo<=t[-1][1]:t[-1]=(t[-1][0],max(t[-1][1],hi))
  else:t+=[(lo,hi)]
  l1.update(rn(lo>>52,((hi-1)>>52)+1));l2.update(rn(lo>>40,((hi-1)>>40)+1))
  if m:l3.update(rn(lo>>28,((hi-1)>>28)+1))
  else:
   if lo%BLK:l3.add(lo>>28)
   if hi%BLK:l3.add((hi-1)>>28)
 ok&=all(b1<=a2 or b2<=a1 for a1,b1 in I for a2,b2 in J);T=2+len(l1)+len(l2)+len(l3);o=PG*T+sum(b-a for a,b in I)+sum(b-a for a,b in J);return ok,T,o-sum(ROLE_B.values()),o
PT=_pt();MACH_STATE,MMU_ON,PT_DESC=2**12,64,4;NI_ITEM=MACH_STATE+MMU_ON+8*2**10+PT_DESC*4096*PT[1]+2**21;check("c39",PT[0]and NI_ITEM<=MACH_INIT and BSEL+5+3*1024<=GEN_INIT and(BSEL+5)+(6*TINIT_ENTRY+1+4)+8<=YQ_INIT and BSEL+5+8<=TB_BANK,"tables %d"%PT[1]);NEAR_B=ROLE_B["N"]+ROLE_B["D"]+PT[3]-sum(ROLE_B.values());PEAK=NEAR_B+sum(N_BANK[w]*ROLE_B[w]for w in"TYSH");check("c40",len(REGISTRY)==18 and[sum(r_==k for r_,_,_,_,_ in REGISTRY)for k in KINDS]==[4,4,2,1,5,2]and all(lo<hi and _bld_ok(b)for _,_,lo,hi,b in REGISTRY)and disjoint_va([(lo,hi)for _,_,lo,hi,_ in REGISTRY])and sum(ROLE_B.values())<PT[3]<2**56 and ROLE_B["H"]==2**53 and all(2**b0*2**b1==N_BANK[w]for w,_,(b0,b1),_ in WINS)and all(a%PG==0 and _mvz(a)for a in list(SEL.values())+[H_BASE]),"near 2^%.4f PA 2^%.4f"%(lg(NEAR_B),lg(PT[3])));check("c41",YR_BASE-64+(YB+1)*YRSL<=33*TB and YROUNDS*YPC*512+64<=YRSL and NEAR1[1]<=SEL["T"]and SEL["H"]+PG<=CL_BASE and CL_BASE+32*LIST_CAP+64<=GS_BASE and GS_CONTENT<=2**30 and GS_BASE+2**30<=RP_BASE and RP_BASE+N_POOLS*(1024*128+256)<=YQ_BASE-64 and YQ_BASE+512*(B_CAP+1)<=EB_START and EB_END<=CN_BASE and CN_BASE+64*Q1GD<=CMB and CMB+64*NBATCH_MINE<=AXB and AXYX_HI+12*2**48==FWB and FWB+2*2560*TB==AXY_END<=H_BASE and H_BASE+2**53==LOW and OV_BASE+32*(OVF_CAP+1)==H_BASE+2**53 and YQ_LIM+448<=YQ_BASE+512*(B_CAP+1));check("c42",Q1GD<=2**37 and Q2<2**24 and YROUNDS*YPC==L_CHAIN and PTS_MINE_CHAIN<2**18 and NBATCH_MINE<2**44)
def laa(d,x,o):return P2(d)[0]>=x and P2(d-Fr(1,o))[1]<x
def least_grid_up(x,o):
 k=math.ceil(lg(x)*o)
 for _ in rn(20):
  d=Fr(k,o)
  if laa(d,x,o):return d
  k+=1 if P2(d)[0]<x else-1
 raise AssertionError("grid bound not certified")
DECLARED=dc(time_log2=CLAIM_E,memory_log2_bytes=least_grid_up(PEAK,10**4),preprocessing_log2=least_grid_up(Fr(E_PRE,C),10**4),success=THRESH,per_leaf=Fr(cd(PLE.numerator*10**4,PLE.denominator),10**4));GEN_ITEM=256*(DBL+MADD_)+2*(INV+3*MULX+2*SQRX+10*ADDX)+64;GEN_TEST,GEN_CTRL,NEGT,RED=256*2*13+64,256*6+64,30,9;X_ENTRY=MADD_+3*SQRX+11*MULX+2*ADDX+2*NEGT+17*RED+32;check("c43",GEN_ITEM+GEN_TEST+GEN_CTRL+X_ENTRY==GEN_MUL and GEN_MUL+4==GEN_MULX);EIT=MX*(2*4+3+1)+3+MX//2;TEMPLATE_ITEM=3*MULX+SQRX+24+7+20;check("c44",INV+3*MULX+2*SQRX+10*ADDX+4<CHAIN_ADD and 56*(2*MULX+ADDX+4)+TEMPLATE_ITEM+X_ENTRY<XSETUP and EIT<=MX*ENTRY_X,"entries %d/%d"%(EIT,MX*ENTRY_X));BS_ITEM=dc(rho=COIN_SCALAR,ycalls=cd(4*YROUNDS,YB),meter=11,calls=41,selects=2*(BSEL+5),pool=96,ycursor=BSEL+5+16,cand=2*(4+3),loop=32);check("c45",sum(BS_ITEM.values())+(BS_ITEM["meter"]+BS_ITEM["cand"])*(B_CAP//CH)<=BATCH_SETUP and BSEL+5+24<=Y_PROD_FIX and 2*STOP_CHECK+64<=GROUP_CTRL and 16<=HZ_D1 and 8+BSEL+5+16<=XDIST and XEMIT_FIX+5+4+3<=64 and BSEL+5<=64 and 2*(3+5)<=CHUNK_SAVE,"batch items %d"%sum(BS_ITEM.values()));CAND_SEL=4*(BSEL+5)+2*(BSEL+5)+2*64+4*64;check("c46",CAND_SEL<=CAND_CTL and RESOLVE==2*CAND_ITEM+1 and 128+64<=REPORT);COIN_ITEM=2*(4*2+5+4)+16+2+1+2;P_COIN_FAIL=Fr(2**256-r+1,2**256)**2;INVC_ITEM,STOP_ITEM,HEADER_ITEM,REQ_ITEM,OUT_ITEM=2+1+6+6,5,17*4+9,16+12+600+16,16;check("c47",COIN_ITEM<=COIN_SCALAR and INVC_ITEM<=INV_CALL and STOP_ITEM<=STOP_CHECK and HEADER_ITEM<=HEADER and REQ_ITEM+GEN_MUL<=REQ_SERVE and OUT_ITEM<=OUTPUT and P_COIN_FAIL<Fr(1,2**254));print("q1 %d  q2 %d  S 2^%.6f  mu %.9f  lo_run %.12f  p_gate %.9f  success >= %.12f"%(Q1GD,Q2,lg(SG_FULL),float(TFULL["mu"]),math.floor(LO_RUN*10**12)/10**12,float(P_GATE),math.floor(GATED_LO*10**12)/10**12));print("E_mine 2^%.6f  E_run 2^%.6f insns (2^%.6f units)  claim %d/10^15  per leaf %.4f  mem 2^%.7f  pre 2^%.7f  grids %s"%(lg(E_MINE),lg(E_RUN),lg(Fr(E_RUN,C)),CLAIM_E*10**15,float(PLE),lg(PEAK),lg(Fr(E_PRE,C)),"/".join(format(float(least_grid_up(ROLE_B[k],10**4)),".4f")for k in"NTYSH")))
def _ps(a):return sum(PHI[g]*a[g]for g in rn(G_GROUPS))
LEDGER=[MACH_INIT+NGEN*GEN_INIT+N_YBANK*YQ_INIT,NGEN*GEN_MUL,E_MINE,NGEN*(3*cd(Q1GD,M_G)+1)*COIN_SCALAR,NGEN*MX*XSETUP,NGEN*6*CH*ENTRY_X+NGEN*HEADER,N_TBANK*(XTBLD_ALL+TB_BANK),TBLD,ZERO_INIT,_ps(WGL),_ps(WGI),E_FIX_WORK,JROW_E-DEFJ,_ps(WGG),SXV,K2ROW_E,ELIST*OCC_LS1,DEFV+DEFK2+DEFL+DEFJ,E_CAND_WORK,E_REQ_WORK,OUTPUT];LEDGER_CELLS=[format(float(least_grid_up(v,10**4)),".4f")for v in LEDGER];E_RUN_6=least_grid_up(E_RUN,10**6);E_RUN_U6=least_grid_up(Fr(E_RUN,C),10**6);E6=least_grid_up(E_TOTAL,10**6);print("E_total %.6f  E_RUN_6 %.6f  E_RUN_U6 %.6f  tables %d"%(float(E6),float(E_RUN_6),float(E_RUN_U6),PT[1]));check("c48",len(LEDGER)==21 and sum(LEDGER)==E_RUN and P2(CLAIM_E)[0]>=Fr(E_TOTAL,C)and P2(E_RUN_6)[0]>=E_RUN and P2(E_RUN_U6)[0]>=Fr(E_RUN,C)and P2(E6)[0]>=E_TOTAL,"cells "+" ".join(LEDGER_CELLS));d8=open(os.path.join(cl,"proof.md")).read();d7=open(os.path.join(cl,"success.md")).read();cq=open(os.path.join(cl,"submission.json")).read();d6=json.loads(cq)["claims"][0]["score"]
INS_N=sum(PHI[_g]*SG[_g]for _g in rn(G_GROUPS))+S_MINE_TOT;SENS_INS=lg(Fr(E_TOTAL+P_GATE*INS_N,E_TOTAL));print("BM sensitivity +%.6f bits per instruction-equivalent added to every hash select (2^%.4f expected selects)"%(SENS_INS,lg(INS_N)))
# DCZID_EL0.BS=9 (2,048-byte DC ZVA blocks) is assumed by hzero; with 64-byte blocks the same loop needs ZP64
# instructions for 131,072 B in place of ZERO_PASS, in both zeroing passes
ZP64=2*(8*ZW//64)+2;DZ64=(ZP64-ZERO_PASS)*N_HBANK*((2**(SLOT_BITS+5)+2**(SLOT_M+5))//(8*ZW));CLAIM_Z64=least_grid_up(Fr(E_TOTAL+P_GATE*DZ64,C),10**4);print("claim with 64-byte DC ZVA blocks %.4f (%d instructions a pass)"%(float(CLAIM_Z64),ZP64));check("c56",ZP64==4098 and ZERO_PASS==2*(8*ZW//2048)+2 and ZERO_INIT+ZERO_MINE+DZ64==N_HBANK*(((2**(SLOT_BITS+5)+2**(SLOT_M+5))//(8*ZW))*ZP64+2*HZ_BANK)+2*2**D1*HZ_D1 and CLAIM_Z64>CLAIM_E)
W_EVERY=GATE_INSNS+E_SHARED_REST+W_MINE_DET_MAX+W_DET_MAX;EVERY_U=least_grid_up(Fr(W_EVERY,C),10**4);print("every-run ceiling 2^%.4f units (W_EVERY 2^%.4f insns)"%(float(EVERY_U),lg(W_EVERY)));check("c55",W_EVERY>E_RUN and P2(EVERY_U)[0]>=Fr(W_EVERY,C))
# selects executed, per window: a round of Phases 1-3 (hash: one a leaf; Y: one a stored record and one a member
# visit; state: one a member visit and one a producing cluster; tables: one a member visit) and set-up (hash: one a
# mined point and one a bank in each of the two zeroing passes; Y: one a sub-ring; state: three visits a
# generator; tables: one a bank)
SELN=dc(H=(NGEN*LEAVES*Q2,S_MINE_TOT+2*N_HBANK),Y=(N_CLUST*L_CHAIN+NGEN,N_YBANK),S=(NGEN+N_CLUST,3*NGEN),T=(NGEN,N_TBANK));check("c54",SELECT==3 and BSEL==SELECT+1 and all(a>0 and b>0 for a,b in SELN.values())and SELN["H"][0]*Q1GD==SG_FULL)
def g4(x):return format(float(least_grid_up(x,10**4)),".4f")
def grp(n):return format(n,",")
def dfl(x,n):return"%d.%0*d"%((int(x*10**n)//10**n,n,int(x*10**n)%10**n))
FIGV=dc(claim=dfl(CLAIM_E,15),num=str(CLAIM_E.numerator),den=str(CLAIM_E.denominator),kgate=str(K_GATE),pgate="%.9f"%float(P_GATE),lorun=dfl(LO_RUN,12),mu="%.9f"%float(TFULL["mu"]),thresh="0.39",erun=dfl(E_RUN_6,6),erunu=dfl(E_RUN_U6,6),etot=dfl(E6,6),slog=dfl(least_grid_up(SG_FULL,10**6),6),perleaf=dfl(DECLARED["per_leaf"],4),mem=dfl(DECLARED["memory_log2_bytes"],4),pre=dfl(DECLARED["preprocessing_log2"],4),every=dfl(EVERY_U,4),sens="%.6f"%SENS_INS,z64=dfl(CLAIM_Z64,4),zp64=grp(ZP64),xentry=grp(X_ENTRY),q1=grp(Q1GD),q2=grp(Q2),groups=str(G_GROUPS),ovf=grp(OVF_CAP),lcap=grp(LIST_CAP),mpool=grp(M_POOL),nbm=grp(NBATCH_MINE),tables=grp(PT[1]),bsum=grp(sum(BS_ITEM.values())),insn="%.4f"%lg(INS_N),candsel=str(CAND_SEL),gap="%.1f"%(float(EVERY_U)-float(CLAIM_E)),gnear=g4(NEAR_B),gT=g4(ROLE_B["T"]),gY=g4(ROLE_B["Y"]),gS=g4(ROLE_B["S"]),gH=g4(ROLE_B["H"]),pa=g4(PT[3]),**{"sel"+w+"01"[i]:g4(SELN[w][i])for w in"TYSH"for i in(0,1)},**{"row%d"%(i+1):c for i,c in en(LEDGER_CELLS)})
FIGN=dc(**{'claim':(2,0,2),'num':(1,0,1),'den':(1,0,1),'kgate':(1,0,1),'pgate':(1,0,1),'lorun':(1,1,1),'mu':(0,1,0),'thresh':(2,2,2),'erun':(1,0,0),'erunu':(1,0,1),'etot':(1,0,0),'slog':(1,0,1),'perleaf':(1,0,1),'mem':(1,0,1),'pre':(1,0,1),'every':(1,0,1),'sens':(1,0,0),'z64':(1,0,0),'zp64':(1,0,0),'xentry':(1,0,0),'q1':(1,1,1),'q2':(1,0,0),'groups':(3,0,1),'ovf':(0,1,0),'lcap':(1,0,0),'mpool':(1,0,0),'nbm':(1,0,0),'tables':(2,0,0),'bsum':(1,0,0),'insn':(1,0,0),'candsel':(1,0,0),'gap':(1,0,0),'gnear':(1,0,0),'gT':(1,0,0),'gY':(1,0,0),'gS':(1,0,0),'gH':(1,0,0),'pa':(1,0,0),'selT0':(1,0,0),'selT1':(1,0,0),'selY0':(1,0,0),'selY1':(1,0,0),'selS0':(1,0,0),'selS1':(1,0,0),'selH0':(1,0,0),'selH1':(1,0,0),'row1':(1,0,0),'row2':(1,0,0),'row3':(2,0,1),'row4':(1,0,0),'row5':(1,0,0),'row6':(1,0,0),'row7':(1,0,0),'row8':(1,0,0),'row9':(1,0,0),'row10':(1,0,0),'row11':(1,0,0),'row12':(1,0,0),'row13':(1,0,0),'row14':(1,0,0),'row15':(1,0,0),'row16':(1,0,0),'row17':(1,0,0),'row18':(1,0,0),'row19':(1,0,0),'row20':(1,0,0),'row21':(1,0,0)})
def cnt(text,o):return len(re.findall(r"(?<![0-9.])"+re.escape(o)+r"(?![0-9])",text))
FIGC={k:(cnt(d8,v),cnt(d7,v),cnt(cq,v))for k,v in FIGV.items()};d9=[k for k in FIGV if FIGC[k]!=FIGN.get(k)];check("c49",not d9,"count mismatch: "+",".join(d9));check("c50",set(FIGV)==set(FIGN)and all(sum(v)>0 for v in FIGN.values())and all(FIGN["row%d"%(i+1)][0]>=1 for i in rn(len(LEDGER))),"unstated figure");check("c51",Fr(int(d6["numerator"]),int(d6["denominator"]))==DECLARED["time_log2"]and math.gcd(int(d6["numerator"]),int(d6["denominator"]))==1 and DECLARED["success"]==THRESH)
RATES={"2^-50.51":D_WJ,"2^-8.20":D_LAB,"2^-7.65":3*D_WJ+D_LAB_SYM,"2^-126":D_BATCH,"2^-26":D_DROP,"2^-82":TFULL["cs"],"2^-18.8":TFULL["d_f"],"2^-37.6":TFULL["lam_int"],"2^-29.5":TFULL["d_p2"],"2^-21.7":TFULL["g_def"],"2^-19":max(TFULL["glob"],P_FAIL_ALL_POOLS),"2^-39":P_SAME_A_ALL,"2^-140":max(P_BAD_X,D_XMINE_ALL,P_MINE_ZERO_ALL),"2^-139":P_EXC,"2^-202":P_X_ORBIT_ALL_POOLS,"2^-254":P_COIN_FAIL,"2^-200":P_MINE_CAPS};ANALYTIC={"2^-30","2^-187","2^-190"}
RATEN=dc(**{'2^-50.51':(1,1,0),'2^-8.20':(1,1,0),'2^-7.65':(1,0,0),'2^-126':(1,1,0),'2^-26':(0,2,0),'2^-82':(0,2,0),'2^-18.8':(0,2,0),'2^-37.6':(0,2,0),'2^-29.5':(0,2,0),'2^-21.7':(0,1,0),'2^-19':(0,3,0),'2^-39':(0,1,0),'2^-140':(1,4,0),'2^-139':(0,1,0),'2^-202':(0,1,0),'2^-254':(1,1,0),'2^-200':(0,2,0)})
def rcnt(text,o):return len(re.findall(re.escape(o)+r"(?!\.?[0-9])",text))
RATEC={k:(rcnt(d8,k),rcnt(d7,k),rcnt(cq,k))for k in RATES};da=[set(re.findall(r"2\^-[0-9.]*[0-9]",t))for t in(d8,d7,cq)];cx=[k for k,v in RATES.items()if not(v<=P2(-Fr(k[3:]))[0]and RATEC[k]==RATEN.get(k))];check("c52",not cx and set().union(*da)<=set(RATES)|ANALYTIC,"bad: "+",".join(cx+sorted(set().union(*da)-set(RATES)-ANALYTIC)));check("c53",PSEL_POOL==262160 and OVF_CAP==2**47-1 and OVF_CAP_M==2**34-1 and LIST_CAP==150994944 and DECLARED["time_log2"]==CLAIM_E);print("FAIL: "+", ".join(FAILED)if FAILED else"PASS")
if __name__=="__main__":sys.exit(1 if FAILED else 0)

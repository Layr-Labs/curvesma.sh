import functools,itertools,json,math,os,re,sys;from decimal import Decimal,ROUND_CEILING,getcontext;from fractions import Fraction as Fr;FAILED=[];rn,en,dc=range,enumerate,dict
def check(I,ok,o=""):
 print("%s %s%s"%(I,"ok"if ok else"FAIL","  "+o if o else""))
 if not ok:FAILED.append(I)
def lg(x):f=Fr(x);return math.log2(f.numerator)-math.log2(f.denominator)
def cd(a,b):return-(-a//b)
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
p=2**256-2**32-977;r=int("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141",16);C0=2**32+977;BETA_P=int("7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee",16);LAM_R=int("5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72",16);C=1000;NODES=2**66;M_G=1024;N_CLUST=NODES//M_G;D1,D2,D3,D4=40,42,64,22;F1,F2=1<<D1,1<<D2;N_REL,N_HOME=F1,F1*F2;ALLN=NODES+N_REL+N_HOME;MX=6;YPC=33;YROUNDS=63520;L_CHAIN=YROUNDS*YPC;Q2=MX*L_CHAIN;KEYBITS=64+D4;KBF=700;K1,K2B=128,128;K=K1;HB=50;THRESH=Fr(39,100);LEAVES=6;FWD_INC,FWD_FIX,FWD_REM=Fr(1494,24),7,(8,381,753,1128)
def fwd_count(q2):assert q2%6==0;m,c=divmod(q2//6,4);return FWD_FIX+FWD_REM[c]+1494*m
BWD_INC,BWD_FIX,BWD_ODD=Fr(10226,12),9,5115
def bwd_count(q2):assert q2%6==0;n=q2//6;return BWD_FIX+10226*(n//2)+BWD_ODD*(n%2)
CH,B_CAP,CHK=2**22,2**26-1,2;FB_FIX=max(fwd_count(6*c)-6*c*FWD_INC for c in rn(4))+BWD_FIX+max(0,BWD_ODD-6*BWD_INC);FIX_SITES_QUAD=dc(SHX=18,FDX=12,F5=4,FS=6,FADDW=4,BLK=10,DF=2,RW=8);FIX_EXTRA=dc(SHX=3,FDX=3,F5=3,FS=4,FADDW=4,BLK=6,DF=6,RW=3);EMIT=7;EMIT_INC=LEAVES*EMIT;KEY_W3=Fr(1,1);KFC=10;KEY_INC,KEY_FIX=Fr(801,6),5;FLUSH=10;FLUSH1_X,FIN1_X=8,14;FLUSH2_X,FLUSH2_D=11,9;FLUSH2=FLUSH+1;XEMIT_HOT=14;JPATH,JPATH2,SJOIN=3,2,0;JX=((0,JPATH,JPATH2),)*2+((1,JPATH,JPATH2-1),);INV=22549;INIT3_BKT,TINIT_ENTRY=4,3;FIN_DEST,FIN_SENT,FIN_HELD,FIN_FIX=15,4,2,2;FINH_DEST,FINH_SENT=16,5;H2S_ARR,H2S_ONE=6,1;H2S_FULL=1;H2S_EMPTY=6;H2S_ARR_FIX=17+H2S_ONE
def h2s_count(n):
 if n==K:return H2S_FULL+H2S_ARR*K
 if n==0:return H2S_EMPTY
 return 17+H2S_ONE*(n%2)+H2S_ARR*n
check("c1",h2s_count(K)==769 and max(h2s_count(n)-H2S_ARR*n for n in rn(K))==H2S_ARR_FIX and H2S_EMPTY<=H2S_ARR_FIX);check("c2",KEY_INC==8+48+KEY_W3+6*(4+EMIT)+KFC+Fr(1,2)and KEY_FIX==5);check("c3",all(fwd_count(6*n)+bwd_count(6*n)<=FB_FIX+6*n*(FWD_INC+BWD_INC)for n in rn(1,9))and CH*512==2**31 and(B_CAP+1)*512==2**35);check("c4",CHK>=1+Fr(L_CHAIN,CH)and 6*CH*512<=2**37);check("c5",(3+3+EMIT+1)-(4+EMIT)==JPATH and(4+1+EMIT+1)-(4+EMIT)==JPATH2 and SJOIN==0 and all(sum(x)==JPATH+JPATH2 and max(x)==JPATH for x in JX));check("c6",EMIT==1+1+2+2+1 and XEMIT_HOT==2+4+4+4);DONE=2;RCV,WAKE,WY,EXPECT,SWAP,CTX=1,5,23,39,4+4*15,0;DONE_C=DONE+RCV+WAKE;RCV1,FIN_CNT,H2S_DONE=4,6,4;EXPECTH=78;RDONE=1+H2S_DONE;DONE_C1=DONE+RCV1+WAKE;DLV=6;MSG_LINK=DONE_C+DLV+CTX;FINAL=FIN_DEST+FIN_SENT+DONE_C1+H2S_ARR_FIX+FIN_CNT;check("c7",FLUSH1_X<=FLUSH+DONE and FIN1_X<=FIN_DEST+FIN_SENT+DONE and FLUSH2_X<=FLUSH2 and FLUSH2_D<=FLUSH2 and FIN_HELD<=FINH_SENT and DONE==2 and(DONE_C,DONE_C1,MSG_LINK,FINAL,RDONE)==(8,11,14,54,5)and H2S_ARR_FIX==18);check("c8",5<=WAKE and 13+2*5<=WY and(SWAP,CTX)==(64,0)and 3+5+5+26<=EXPECT and EXPECTH==2*EXPECT and FIN_FIX==2);GEN_MUL=645410;GEN_MULX=GEN_MUL+4;YR_FIX,YR_BLK,YR_BLK_MINE,YB=22556,62398,40552,256;CHAIN_ADD=2**15;CHAIN_SET=2*GEN_MUL+17*CHAIN_ADD;YSND,YDLV,YDLVE=30,32,23;YRING_TOT=YSND+3+(M_G-1)*(YDLV+WY)+YDLVE+WY
def y_group(q):o=cd(q,M_G);J=YROUNDS*(cd(o,YB)*YR_FIX+o*YR_BLK)+o*CHAIN_SET;I=cd(q*L_CHAIN*YRING_TOT,M_G);return J+I
NCAND=dc(A=8,B=6,EB=7);EISR={k:357+31*v for k,v in NCAND.items()};KCL=dc(E=("A",)*6,F=("A",)*3+("B",)*3,G=("A",)*3+("B",)*3)
def tbldk_count(n,cl):return 5+n*(197+6*EISR[cl])
def tbldf_count(n,v):return 4+n*(6+46*v)
XT_SLOTS=6;XT_W,XT_WHI=43,42;XTBLD_PAIRED=6*(4+2*11)+20*2*(4*2**XT_W+2*2**XT_WHI);XTBLD_ALL=6*XTBLD_PAIRED+2*(5*26+200*2**XT_W)+sum(tbldk_count(2**(XT_WHI if k in(2,5)and c!="F"else XT_W),KCL[c][k])for c in"EFG"for k in rn(5 if c=="F"else 6));XSETUP=2**17;ENTRY_X=16;COIN_SCALAR=64;BATCH_SETUP=9544;INV_CALL=16;NODE_INIT=2**31;HEADER=80;BT_BASE,BT_W=0xfff8<<48,43;TBLD=tbldf_count(2**BT_W,7);MSG=128;CAND_CAP=31*2**26;MULX,SQRX,ADDX=105+4,80+3,9+2*9;DBL,MADD_=2*MULX+5*SQRX+10*ADDX,8*MULX+3*SQRX+10*ADDX;CAND_LEAVES=6*CHAIN_ADD+12*MULX;CAND_ITEM=4*8*2*MSG+2*CAND_LEAVES+36*600+2*(384*200+256*(DBL+MADD_)+INV);RESOLVE=2*CAND_ITEM+1;REQ_CAP=125*2**39;REQ_SERVE=2**10+GEN_MUL;RELAY_CAP=2**49-2**43+2**41;BARRIER=2**12*ALLN;OUTPUT=2**20;DELTA_F=2**62;QP,QQ=Q2//2,Q2//4;FIX_BATCH=sum(FIX_SITES_QUAD.values())*CH*6//4;P1=Fr(6*(r-1)+1,r*r);TWB=Fr(2,r-1);SHX_BAD=2**64*2**33;FOLD_BAD=2**66*C0;SQ_BAD=C0*C0+2*C0;FS_P1_BAD=C0;SITES=[(3*QP,SHX_BAD,36),(QP,FOLD_BAD,36),(3*QP,SHX_BAD,36),(QP,FOLD_BAD,36),(QP,FOLD_BAD,36),(3*QP,SHX_BAD,36),(QP,FOLD_BAD,36),(Q2,SQ_BAD,12),(Q2,FOLD_BAD,12),(2*Q2,FS_P1_BAD,6)];SITES_QUAD=dc(SHX=18,FDX=12,F5=4,FS_FA=8);D_SITES=sum(Fr(n*bad*mult,r-1)for n,bad,mult in SITES);D_YGEN=Fr(YROUNDS*(YPC+2*16+2),r-1);D_BATCH=Fr(CH,L_CHAIN)*(D_SITES+Q2*Fr(6,r-1))+CH*D_YGEN;FIX_ROW=C0*2**192;FIX_FOLD=(2**66+1)*2**128*C0;FIX_SQ=2**193;FIX_LK6=15*C0*2**192;FIX_LK7=21*C0*2**192;FIX_SUB=C0*2**192;FIX_DF63=2**128*(2**63+2)*C0;FIX_RW=C0*2**128;FIXUPS=[("BLK",QP,FIX_LK7,36),("SHX",3*QP,FIX_ROW,36),("FDX",QP,FIX_FOLD,36),("SHX",3*QP,FIX_ROW,36),("FDX",QP,FIX_FOLD,36),("FDX",QP,FIX_FOLD,36),("FS",QP,FIX_SUB,36),("SHX",3*QP,FIX_ROW,36),("FDX",QP,FIX_FOLD,36),("DF",QP,FIX_DF63,6),("F5",Q2,FIX_SQ,12),("BLK",Q2,FIX_LK7,12),("BLK",Q2,FIX_LK6,12),("FDX",Q2,FIX_FOLD,12),("FS",Q2,FIX_SUB,6),("FADDW",Q2,FIX_SUB,6),("RW",2*Q2,FIX_RW,6)];E_FIX=sum(Fr(n*bad*mult,r-1)for _,n,bad,mult in FIXUPS)+CHK*D_BATCH*FIX_BATCH;KB=functools.lru_cache(None)(lambda fixed:Fr(6*2**(256-fixed)*(KBF+1),(r-1)*KBF));RHO=Fr(r,r-1)**2;SLF=6*RHO*Fr(2**256,r);SL6=6*L_CHAIN*Fr(2**256,r);D_DROP=2*SLF*L_CHAIN/(B_CAP-2*L_CHAIN-SL6)**2;CLUST_BATCH=M_G*LEAVES*Q2;FIXED_LEAVES=2*CLUST_BATCH;PARTNER_PAIRS=CLUST_BATCH*(CLUST_BATCH-1);FA,FB=int("114ca50f7a8e2f3f657c1108d9d44cfd9",16),int("e4437ed6010e88286f547fa90abfe4c3",16);MH_A,MH_B,ML=20,21,23
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
 lo,hi=min(0,U[0])+min(0,W[0]),max(0,U[0])+max(0,W[0]);o=(-(lo+hi)//2,(X2-U[1]-W[1])//2);C=(g+1)//2;cs=sorted(((i*U[0]+j*W[0],i*U[1]+j*W[1])for i in rn(-4,5)for j in rn(-4,5)if o[0]+lo+i*U[0]+j*W[0]<=C and o[0]+hi+i*U[0]+j*W[0]>=-C and o[1]+i*U[1]+j*W[1]<X2 and o[1]+U[1]+W[1]+i*U[1]+j*W[1]>0),key=lambda t:(abs(t[0])+64*abs(t[1]),t));return dc(X2=X2,C=C,U=U,W=W,o=o,c=cs,kS=(W[0]*o[1]-W[1]*o[0])%p,kT=(U[1]*o[0]-U[0]*o[1])%p)
def eis_red(c,v):
 U,W,o=c["U"],c["W"],c["o"];S=(W[1]*v+c["kS"])%p;T=(c["kT"]-U[1]*v)%p;f1,f2=o[0]+(S*U[0]+T*W[0])//p,o[1]+(S*U[1]+T*W[1])//p
 for l1,l2 in c["c"]:
  if 0<=f2+l2<c["X2"]and abs(f1+l1)<=c["C"]:return f1+l1,f2+l2
 return f1,f2
def eis_cover(c):
 U,W,(o1,o2),X2,C=c["U"],c["W"],c["o"],c["X2"],c["C"];V=[(o1,o2),(o1+U[0],o2+U[1]),(o1+U[0]+W[0],o2+U[1]+W[1]),(o1+W[0],o2+W[1])];E_=list(zip(V,V[1:]+V[:1]));ys=sorted({y for y in[v[1]for v in V]+[e for _,l in c["c"]for e in(-l,X2-l)]if o2<=y<=o2+U[1]+W[1]})
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
def _table3(rs):rot=lambda a,b:(a,b,(-2-a[1]-b[1],-a[0]-b[0]));(e0,e1),(f0,f1),(g0,g1)=rs;al=rot(isum(e0,(-1,0)),isum(e1,(-1,0)));bp=rot(*[(-f[1]-2,-f[0]-1)for f in(f0,f1)]);ga=rot(g0,g1);W=lambda m,t:isum(al[m%3],bp[(m+t)%3],ga[(m+2*t+1)%3],(0,2));wr=[[W(o+i,t)for t,o in en((2,2,0))]for i in(0,1)];nr=[(-wr[0][t][1]-wr[1][t][1]-1,-wr[0][t][0]-wr[1][t][0]+1)for t in rn(3)];return[(min(x[0]for x in r),max(x[1]for x in r))for r in(wr[0],wr[1],nr)]
ERRH,ERRL=_table3([(KT[c][0][0],KT[c][1][0])for c in"EFG"]),_table3([(KT[c][0][1],KT[c][1][1])for c in"EFG"]);NWIN_K,ERR_K,WW_K,DSYM_K=5748,([(-241,237),(-136,133),(-361,360)],[(-399,401),(-234,228),(-592,588)]),(1468,2442,722,1181),-765;WDH,WDL=sum(hi-lo for lo,hi in ERRH),sum(hi-lo for lo,hi in ERRL);WH,WL=max(hi-lo for lo,hi in ERRH)+1,max(hi-lo for lo,hi in ERRL)+1;(ELO,EHI)=(min(ERRH[0][0],ERRH[1][0]),max(ERRH[0][1],ERRH[1][1]));NWIN=2*(EHI-ELO)+(1-2*ELO)+(3+2*EHI)+2*(EHI+1-ELO)+2*(3*EHI+3-3*ELO);D_WJ=TWB*NWIN*2**192;D_LAB=TWB*(WDH*2**(256-MH_A)+WDL*2**(256-ML));D_LAB_SYM=TWB*3*((WH-1)*2**(256-MH_A)+(WL-1)*2**(256-ML));B_HA,B_HB,B_L=2**MH_A+WH-1,2**MH_B+WH-1,2**ML+WL-1;B_LA=2**(ML+1)+WL-1;check("c9",FA*FA-FA*FB+FB*FB==p and(FB*BETA_P-(FA-FB))%p==0 and((FA-FB)*BETA_P+FA)%p==0 and B_HA*B_HB*(B_LA*2**128+2*FA)*(B_L*2**128+4*FA-2*FB)<p*2**(MH_A+MH_B+2*ML+1)*(1+Fr(1,KBF))and(2*B_HA+1)*(2*B_HB+1)*(2*B_LA*2**128+FA)*(2*B_L*2**128+2*FA-FB)<16*p*2**(MH_A+MH_B+2*ML+1)*(1+Fr(1,KBF))and(2**41+WH-1)*(B_LA*2**128+2*FA)*(B_L*2**128+4*FA-2*FB)<p*2**(2*ML+42)*(1+Fr(1,KBF))and B_LA*2**232<p*(1+Fr(1,KBF))and(B_HA*2**192+2*FA)*(B_HB*2**192+4*FA-2*FB)<p*2**(128+MH_A+MH_B)*(1+Fr(1,KBF)));check("c10",NCAND=={k:len(c["c"])for k,c in EIS.items()}and all(eis_cover(c)and all((u[0]+u[1]*BETA_P)%p==0 for u in(c["U"],c["W"]))and c["U"][0]*c["W"][1]-c["U"][1]*c["W"][0]==p for c in EIS.values())and all(sum(EIS[c]["X2"]-1 for c in KCL[k]+(("EB",)if k=="E"else()))<2**128 for k in KCL)and NWIN==NWIN_K,"C_E 2^%.4f C_FG 2^%.4f"%(lg(C_E),lg(C_FG)));check("c11",(ERRH,ERRL)==ERR_K and(WDH,WDL,WH,WL)==WW_K);check("c12",D_WJ<Z2(-5051,100),"2^%.4f"%lg(D_WJ));check("c13",all(sum(n for n,b,_ in SITES if b==B)==SITES_QUAD[k]*QQ for k,B in(("SHX",SHX_BAD),("FDX",FOLD_BAD),("F5",SQ_BAD),("FS_FA",FS_P1_BAD)))and len(SITES)==10);FIX_RUNS={m:sum(n for ms,n,_,_ in FIXUPS if ms==m)for m in FIX_SITES_QUAD};check("c14",FIX_RUNS=={m:FIX_SITES_QUAD[m]*QQ for m in FIX_SITES_QUAD}and sum(FIX_SITES_QUAD.values())==64);check("c15",FIX_BATCH*C0<2**63 and FIX_BATCH*C0>>32==FIX_BATCH+(FIX_BATCH*977>>32),"2^%.3f"%lg(FIX_BATCH*C0));check("c16",D_LAB<Z2(-820,100)and 3*D_WJ+D_LAB_SYM<P2(Fr(DSYM_K,100))[0],"2^%.3f 2^%.4f"%(lg(D_LAB),lg(3*D_WJ+D_LAB_SYM)));check("c17",D_BATCH<Z2(-1258,10)and D_DROP<Z2(-26,1),"2^%.3f"%lg(D_BATCH));G_GROUPS=10752;G_EXP=21
def group_sizes(q1):
 ms=[0];o=G_GROUPS**G_EXP;J=q1**32
 for g in rn(1,G_GROUPS):
  I=J*(g**G_EXP);m=int(q1*((g/G_GROUPS)**(G_EXP/32.0)))+2
  while m**32*o>I:m-=1
  while((m+1)**32)*o<=I:m+=1
  ms+=[m]
 ms+=[q1];qs=[ms[_g+1]-ms[_g]for _g in rn(G_GROUPS)];return qs,ms
HCHG_NEW,HCHG_VISIT,HCHG_END,HCHG_VER=4,4,0,5;HCHG_LS0,HCHG_LSA=2,1;HCHG_LS1=5+HCHG_LSA;check("c18",HCHG_NEW==4 and HCHG_VISIT==4 and HCHG_END==0 and HCHG_VER==5 and HCHG_LS1==6);DLVHB=256
def dsg(n,m):
 if m==0:return 1
 I=n%2;o=m-I;return 5+I*5+o//2*6+o%2*4-(o==0)
def dlvh_copy(o,I,J=K2B):
 if o==0 and I==J:return DLVHB
 return 23-(o>0)+dsg(0,o)+dsg(o,I)+dsg(o+I,J-o-I)
DLVHS=max(dlvh_copy(a,b)for a in rn(K2B+1)for b in rn(K2B+1-a)if(a,b)!=(0,K2B));BLOCKH=FLUSH2+DONE_C+DLVHB;FINALH=FINH_DEST+FINH_SENT+DONE_C+DLVHS+MSG_LINK+(DLVHS-DLVHB);check("c19",DLVHB==256 and DLVHS==427 and BLOCKH==FLUSH2+DONE_C+DLVHB and FINALH==641,"DLVHS %d FINALH %d"%(DLVHS,FINALH));ZW=16384;ZERO_PASS=2*(ZW//256)+2;DLVH0=15;STOP_CHECK=5;NBAR_COUNT=8+19*G_GROUPS+1;NBAR_ALLOW=204297;PHASE_CTRL_G=2**23;KB_BUCKET=KB(D1+D2+HB);PAD_C=max(3*EHI+3,EHI+1-ELO,1-3*ELO);WRAPT_C=483
def _npad(i,m,d=PAD_C):
 if not 0<=i<m:return 0
 lo=max(0,i-d);A,B=(lo,min(m-1,m+d-2*i)),(max(lo,m-d-i),min(m-1,2*m+d-2*i));return sum(max(0,y-x+1)for x,y in(A,B))-(max(0,min(A[1],B[1])-max(A[0],B[0])+1)if A[0]<=A[1]and B[0]<=B[1]else 0)
_LC=[(0,0),(1,-PAD_C),(0,2**64-1),(-2,2**64+PAD_C),(-1,2**64-PAD_C),(-2,2*2**64+PAD_C)];co=sorted({0,2**64}|{t+e for a,b in _LC for cg,g in _LC if a!=cg for t in[(g-b)//(a-cg)]for e in rn(-3,4)if 0<=t+e<=2**64});NC_C=sum((_npad(s,2**64)+_npad(e-1,2**64))*(e-s)//2 for s,e in zip(co,co[1:]));ROWMAX132=cd(NC_C,1<<42)+(len(co)-1)*2**64;CNT132=cd(B_HA*B_HB*(ROWMAX132+(1<<22)*3*B_HA),1<<50)+(1<<22)*4*B_HA*B_HB+2*WRAPT_C*(1<<14)*B_HB+(1<<22)*B_HA*(6*WRAPT_C+8);DPAR=max(FA,2*FA-FB);KB_C132=Fr(6*CNT132*((B_LA*(1<<128)+2*DPAR)*(B_L*(1<<128)+2*DPAR)),p*(r-1));check("c19c",WRAPT_C>=max(EHI,1-ELO,1-2*ELO,2+2*EHI)and all((_npad(e-1,2**64)-_npad(s,2**64))*((s+e-1)//2-s)==(_npad((s+e-1)//2,2**64)-_npad(s,2**64))*(e-1-s)for s,e in zip(co,co[1:]))and max(sum(_npad(h*16+r_,256,12)for h in rn(16))for r_ in rn(16))<=sum(_npad(i,256,12)for i in rn(256))//16+6*256 and KB_C132<Fr(34,100)*KB_BUCKET);KB64=KB(64);R2_UP_TIB=5460;R2_TIB=38220;CAP_H=((R2_TIB*(2**40))//(32*K2B))*K2B;P1_GOOD=Fr(6*(r-1),r*r)*(1-3*D_WJ-D_LAB_SYM);K_U_CLUST=LEAVES*M_G*MX;W_ELL=1+47*(3*D_WJ+D_LAB_SYM);SAMEX_FAC=Fr(CLUST_BATCH-K_U_CLUST,CLUST_BATCH-1)*Fr(r,r-2)*W_ELL+Fr(96*(K_U_CLUST-1),CLUST_BATCH-1);EF_Q=E_FIX*C0/(2**32*(DELTA_F-33*FIX_BATCH));KB_M=KB(KEYBITS+D1+D2);KB_D1=KB(D1);SIP=240*RHO*KB_M;LOG2_POOLS=14;N_POOLS=1<<LOG2_POOLS;R_POOL=N_POOLS//M_G;K_GEN=NODES//N_POOLS
def gf1024_mul(a,b):
 o=0
 for I in rn(10):
  if b>>I&1:o^=a<<I
 for I in rn(18,9,-1):
  if o>>I&1:o^=0x409<<I-10
 return o
def pool_of_gen(g):c,m=divmod(g,M_G);return R_POOL*m+((c&15^gf1024_mul(c>>4&0x3FF,m))&15)
def _vand_gf1024_ok():
 b=2
 for _ in rn(10):b=gf1024_mul(b,b)
 if b!=2:return False
 for _d in 2,5:
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
 for k in 1,2,3,4:
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
SAME_M_NODE_PAIRS=M_G*(N_CLUST*N_CLUST*(R_POOL-1)//R_POOL)//2;SHARED_BOTH_CLUST_PAIRS=(N_CLUST*N_CLUST//(R_POOL*R_POOL)-N_CLUST)//2;ACP=N_CLUST*(N_CLUST-1)//2;GOOD_NODE_PAIRS=SAME_M_NODE_PAIRS+SHARED_BOTH_CLUST_PAIRS*(M_G*(M_G-1)//2)+(ACP-SHARED_BOTH_CLUST_PAIRS)*(M_G*(M_G-1));YR_MINE=992;PTS_MINE_CHAIN=YR_MINE*YPC;NBATCH_MINE=31417985*YB;PTS_PER_GEN_MINE=NBATCH_MINE*PTS_MINE_CHAIN;M_POOL=K_GEN*PTS_PER_GEN_MINE;S_MINE_TOT=NODES*PTS_PER_GEN_MINE;N_CHAINS_POOL=K_GEN*NBATCH_MINE;PPMC=PTS_MINE_CHAIN*(PTS_MINE_CHAIN-1)//2;P_BUCKET_MAX=Fr(6*(1<<63),r-1);P_SAME_CHAIN_POOL=N_CHAINS_POOL*PPMC*P_BUCKET_MAX;P_CHAIN_COUPLE_POOL=N_CHAINS_POOL*PPMC*(M_POOL*M_POOL+(N_CHAINS_POOL-1)*PPMC)*(P_BUCKET_MAX**2);GAMMA_MINE_CROSS=M_POOL*(M_POOL-1)//2-N_CHAINS_POOL*PPMC;P_BDAY_HIT_MIN=Fr((r-1)*(r-1),(1<<193)*r*r)-Fr(6*(r-1),r*r);MU_POOL_LO=GAMMA_MINE_CROSS*P_BDAY_HIT_MIN;CS_POOL=GAMMA_MINE_CROSS*(2*PTS_MINE_CHAIN*M_POOL)*(2*P_BUCKET_MAX*P_BUCKET_MAX);_P_A102=Fr(6*(1<<103),r-1);_P_A174=Fr(6*(1<<(256-D1-D2)),r-1);_P_A215=Fr(6*(1<<256-D1),r-1);P_CELL154_MASK_POOL=(Fr(M_POOL*(M_POOL-1),2)*M_POOL+2*N_CHAINS_POOL*PPMC*M_POOL)*P_BUCKET_MAX*_P_A102;POIS_LE3=enh(MU_POOL_LO)*sum((MU_POOL_LO**_k)/math.factorial(_k)for _k in rn(4));P_FAIL_POOL=POIS_LE3+CS_POOL+P_CELL154_MASK_POOL+P_SAME_CHAIN_POOL+P_CHAIN_COUPLE_POOL;P_FAIL_ALL_POOLS=N_POOLS*P_FAIL_POOL;LAM_ORD=N_POOLS*(M_POOL*M_POOL)*P_BUCKET_MAX;POINTWISE_X_MAX=M_POOL*M_POOL*P_BUCKET_MAX;SAME_OFF_VAND=(2*LEAVES*MX)**2+(M_G//R_POOL)*(LEAVES*MX)**2+Fr(K_U_CLUST*K_U_CLUST,N_CLUST);P_BAD_X=65536*sum((LAM_ORD**(_d-1))*(2*LAM_ORD/r)for _d in(2,3,4))+65536*sum((LAM_ORD**(_d-1))*N_POOLS*(M_POOL*M_POOL)*Fr(30*(1<<63),r*r)for _d in(2,3));P_X_ORBIT_ALL_POOLS=4*LAM_ORD*(6*LAM_ORD/r);P_MINE_ZERO_ALL=Fr(S_MINE_TOT,r);P_SAME_A_ALL=N_POOLS*Fr(M_POOL*(M_POOL-1)//2*6,r-1);D_XMINE_ALL=NODES*NBATCH_MINE*Fr(YR_MINE*(YPC+34),r-1);_MINE_Y_MEAN=S_MINE_TOT*_P_A174;_MINE_K_CH=CAP_H//PTS_MINE_CHAIN;P_MINE_HOME_CAP=Fr(1,1<<200);ECMX=N_POOLS*Fr(M_POOL*(M_POOL-1),2)*_P_A102;E_CAND_MINE=cd(ECMX.numerator,ECMX.denominator);E_VIS_MINE=cd((S_MINE_TOT*(S_MINE_TOT-1)//2)*6*(1<<(256-131)),r-1);E_K2_MINE=cd(N_POOLS*(M_POOL*(M_POOL-1)//2)*6*(1<<(256-131)),r-1);_OTH_H_OUT=N_POOLS*Fr(M_POOL*(M_POOL-1),2)*_P_A174*_P_A102+2*PTS_MINE_CHAIN*M_POOL*_P_A102;_OTH_G_OUT=2*PTS_PER_GEN_MINE*M_POOL*_P_A102;_M_R_OUT=18*N_POOLS*Fr(M_POOL*(M_POOL-1),2)*_P_A215*_P_A102;cc=1<<20;_RL_CR_MINE=36*PTS_MINE_CHAIN*M_POOL*_P_A102;_VAR_R_OUT=324*N_POOLS*M_POOL*_P_A215*(1+RHO*RHO*PTS_MINE_CHAIN*_P_A215)*(M_POOL*_P_A102)**2+18*PTS_MINE_CHAIN*_M_R_OUT;P_MINE_CAPS=P_MINE_HOME_CAP+LAM_ORD*(_OTH_H_OUT/Fr(CAND_CAP-2*PTS_MINE_CHAIN)+2*_OTH_G_OUT/Fr(REQ_CAP-4*PTS_MINE_CHAIN)+_RL_CR_MINE/cc+_VAR_R_OUT/Fr(RELAY_CAP-36*PTS_MINE_CHAIN-cc-_M_R_OUT)**2);check("c19b",_vand_gf1024_ok()and _bx_tuples_ok()and all(pool_of_gen(1024*(_p&15)+(_p>>4))==_p for _p in rn(N_POOLS))and NBATCH_MINE%YB==0 and LN2[0]>Fr(2,3)and _MINE_K_CH-Fr(3,2)*_MINE_Y_MEAN-82>200 and MU_POOL_LO>56 and P_FAIL_ALL_POOLS<Z2(-28,1)and P_BAD_X<Z2(-140,1)and P_MINE_CAPS<Z2(-19,1),"mu_p %.2f p_fail 2^%.2f p_caps 2^%.2f"%(float(MU_POOL_LO),lg(P_FAIL_ALL_POOLS),lg(P_MINE_CAPS)));EXCH_ROW=ALLN*4*2**10;HR_INIT=N_HOME*(INIT3_BKT+1+DLVH0);VISIT=HCHG_VISIT+Fr(HCHG_END,2);GEN_MINE_Y=YR_MINE*((NBATCH_MINE//YB)*YR_FIX+NBATCH_MINE*YR_BLK_MINE)+NBATCH_MINE*(CHAIN_SET+2*COIN_SCALAR+64);GEN_MINE_NODE=64+GEN_MINE_Y+PTS_PER_GEN_MINE*XEMIT_HOT;_P_MINE_DEF=D_XMINE_ALL+P_MINE_ZERO_ALL+NODES*(2*NBATCH_MINE+1)*Fr(1,2**254);W_MINE_DET_MAX=S_MINE_TOT*CAP_H*(VISIT+HCHG_VER)+N_HOME*CAND_CAP*(HCHG_LS1+RESOLVE+2*REQ_SERVE+8*MSG);ce=_P_MINE_DEF*W_MINE_DET_MAX;DEF_MINE=cd(ce.numerator,ce.denominator);PSEL_POOL=16+1024*256;ROUTE_MINE=(S_MINE_TOT*(FLUSH+DONE_C1)//K1+S_MINE_TOT*H2S_ARR+cd(S_MINE_TOT,K1)*H2S_FULL+NODES*(FINAL*F1+FIN_FIX)+cd(S_MINE_TOT*BLOCKH,K2B)+N_REL*(FINALH*F2+FIN_FIX)+N_HOME*EXPECTH+N_REL*(EXPECT+RDONE)+S_MINE_TOT*(HCHG_NEW+RCV)+VISIT*E_VIS_MINE+HCHG_VER*E_K2_MINE+HCHG_LS1*E_CAND_MINE+E_CAND_MINE*(RESOLVE+2*REQ_SERVE+8*MSG)+DEF_MINE+EXCH_ROW+N_POOLS*(PSEL_POOL+1024*(5*MSG+32)+MSG+4+8*(REQ_SERVE+4*MSG))+ALLN*(MSG+STOP_CHECK)+NODES*MX*(4*MSG)+3*N_HOME+6*S_MINE_TOT+HR_INIT+NODES*(INIT3_BKT*F1+2)+N_REL*((INIT3_BKT+TINIT_ENTRY)*F2+64)+5*BARRIER);E_MINE=NODES*GEN_MINE_NODE+ROUTE_MINE;ck=P_FAIL_ALL_POOLS+P_SAME_A_ALL+P_BAD_X+D_XMINE_ALL+P_MINE_ZERO_ALL+P_X_ORBIT_ALL_POOLS+P_MINE_CAPS;d5=P1*P1+Fr(36,(r-1)**2);cf=KB(D1+D2);cz=Fr(K_U_CLUST*L_CHAIN,r-1)*cf**2;cp=SAMEX_FAC*PARTNER_PAIRS+CLUST_BATCH*CLUST_BATCH+96*K_U_CLUST*K_U_CLUST;_FDIFF0=SAME_OFF_VAND+320*K_U_CLUST*K_U_CLUST*((POINTWISE_X_MAX/MX)**4)*KB_M+2*CLUST_BATCH*CLUST_BATCH*RHO*KB_M;cr=1+Fr(1,r)
def spx(q1p):n0p=LEAVES*q1p*Q2;Sp=NODES*n0p;gamma=GOOD_NODE_PAIRS*(n0p*n0p);I=gamma*P1;mu=gamma*P1_GOOD;cs=gamma*(2*CLUST_BATCH*Sp)*d5;e_hi=enh(mu);d_f=SLF*q1p*(EF_Q+CHK*D_BATCH);m_p2=Sp*cf+1;v_p2=Sp*(K_U_CLUST*cf+cz)+1;d_p2=v_p2/Fr(CAP_H-FIXED_LEAVES-m_p2)**2;lam_int=RHO*(Sp+FIXED_LEAVES)*KB_M;alpha_cross=(2*K_U_CLUST+SAMEX_FAC*(CLUST_BATCH-1))*Sp+cp;other_h=RHO*(Fr(Sp)*Sp*cf*KB_M+alpha_cross*KB_M+2*L_CHAIN*cf)+SIP*Sp;d_cl=other_h/(CAND_CAP-76);other_g=RHO*KB_M*(Fr(n0p)*Sp*2+2*LEAVES*MX*Sp+96*K_U_CLUST*n0p+96*LEAVES*MX*K_U_CLUST);d_rq=other_g/(REQ_CAP-2*L_CHAIN-600);rl_cap_cr=RELAY_CAP//4;rl_cr=18*FIXED_LEAVES*Fr(Sp)*KB_M;o=RHO*Fr(Sp)*Sp*KB_D1*KB_M;rl_fdiff=_FDIFF0+2*CLUST_BATCH*Fr(Sp)*KB_M;v_rl=648*(o*rl_fdiff+CLUST_BATCH*Fr(Sp)*RL_SAME_PL);d_rl=rl_cr/rl_cap_cr+v_rl/cn;per_leaf=D_BATCH+d_f+d_p2+d_cl+d_rq+d_rl+D_DROP;rows_=I*cr*(Fr(1,r)+lam_int+2*per_leaf);g_def=NODES*q1p*CHK*D_BATCH;glob=ck+NODES*(3*q1p+2*NBATCH_MINE+1)*Fr(1,2**254)+g_def;lo=1-e_hi-cs-rows_-glob;return dc(mu=mu,lo=lo,cs=cs,rows=rows_,glob=glob,Sp=Sp,n0p=n0p,d_p2=d_p2,m_p2=m_p2,d_f=d_f,lam_int=lam_int,d_cl=d_cl,d_rq=d_rq,d_rl=d_rl,g_def=g_def)
def gfl(q):return spx(q)["lo"]
Q1GD=103000000000;Q1G,MQ1=group_sizes(Q1GD);RL_SAME_PL=W_ELL*Fr(r,r-2)*(RHO*Fr(CLUST_BATCH-1,2)*KB_D1*KB_M+96*K_U_CLUST*KB_M+SIP);cy=[NODES*LEAVES*q*Q2 for q in Q1G];SG_MAX=max(cy);SPROD_MAX=max(cy[g]*(NODES*LEAVES*MQ1[g+1]*Q2)for g in rn(G_GROUPS));M_VG_MAX=18*(RHO*Fr(SPROD_MAX)*KB_D1*KB_M+Fr(SG_MAX)*RL_SAME_PL);cn=(Fr(RELAY_CAP-18*FIXED_LEAVES-RELAY_CAP//4)-M_VG_MAX)**2;N0G=LEAVES*Q1GD*Q2;SG_FULL=NODES*N0G;LO_RUN=gfl(Q1GD);GB,GATE_INSNS=128,256;THRESH_GATE=THRESH;K_GATE=cd(THRESH_GATE.numerator*LO_RUN.denominator*(1<<GB),THRESH_GATE.denominator*LO_RUN.numerator);P_GATE=Fr(K_GATE,1<<GB);GATED_LO=P_GATE*LO_RUN;check("c20 %d"%Q1GD,LO_RUN>=Fr(9985,10000),"lo_run %.12f"%(math.floor(LO_RUN*10**12)/10**12));check("c21",0<K_GATE<(1<<GB)and GATED_LO>=THRESH_GATE and THRESH_GATE>=THRESH and Fr(K_GATE-1,1<<GB)*LO_RUN<THRESH_GATE and 2*COIN_SCALAR+8+3+4+16<=GATE_INSNS,"p_gate %.9f gated_lo %.12f"%(float(P_GATE),math.floor(GATED_LO*10**12)/10**12));check("c22",all((MQ1[g]**32)*(G_GROUPS**G_EXP)<=(Q1GD**32)*(g**G_EXP)<((MQ1[g]+1)**32)*(G_GROUPS**G_EXP)for g in rn(1,G_GROUPS)))
def pq(q):n=LEAVES*q*Q2;g=GOOD_NODE_PAIRS*n*n;return g*P1_GOOD,2*cr*g*P1*(NODES*n*(K_U_CLUST*cf+cz)+1)/Fr(CAP_H-FIXED_LEAVES-NODES*n*cf-1)**2
V=[1-spx(k)["lo"]-enh(pq(k)[0])-pq(k)[1]for k in rn(6)];DF=[]
while V:DF+=[V[0]];V=[b-a for a,b in zip(V,V[1:])]
PDR=math.lcm(*[x.denominator for x in DF]);PRN=[int(x*PDR)for x in DF]
def phig(q):m,t=pq(q);e=enh(m)+t;return e.numerator*PDR+sum(c*math.comb(q,k)for k,c in en(PRN))*e.denominator,e.denominator*PDR
TPN=[phig(m)for m in MQ1[:-1]]
PHI=[Fr(1)]+[Fr(1)if N>=D else Fr(cd(N<<256,D),1<<256)for N,D in TPN[1:]];TFULL=spx(Q1GD);check("c23 %.6f"%lg(SG_FULL),SG_FULL==NODES*6*Q1GD*Q2==TFULL["Sp"]);check("c24 %.9f"%float(TFULL["mu"]),TFULL["mu"]==GOOD_NODE_PAIRS*(N0G*N0G)*P1_GOOD and GATED_LO>=THRESH);check("c25",TFULL["cs"]<Fr(1,2**82),"2^%.4f"%lg(TFULL["cs"]));check("c26",TFULL["lam_int"]<Z2(-366,10),"2^%.3f"%lg(TFULL["lam_int"]));check("c27 %.3f"%lg(Q1GD*E_FIX),Q1GD*E_FIX<Z2(39,1)and DELTA_F>2**20*Q1GD*E_FIX+FIX_BATCH);check("c28",TFULL["d_f"]<Z2(-188,10),"2^%.3f"%lg(TFULL["d_f"]));check("c29",TFULL["d_p2"]<Z2(-336,10),"2^%.3f"%lg(TFULL["d_p2"]));check("c30",TFULL["d_cl"]<Z2(-158,10)and TFULL["d_rq"]<Z2(-158,10)and TFULL["d_rl"]<Z2(-158,10),"2^%.3f 2^%.3f 2^%.3f"%(lg(TFULL["d_cl"]),lg(TFULL["d_rq"]),lg(TFULL["d_rl"])));check("c31",TFULL["g_def"]<Z2(-217,10)and TFULL["glob"]<Z2(-19,1),"2^%.3f"%lg(TFULL["glob"]));check("c32",PHI[0]==1 and all(0<=x<=1 for x in PHI)and all(PHI[g]>=PHI[g+1]for g in rn(G_GROUPS-1))and all(PHI[g].numerator*TPN[g][1]>=min(TPN[g])*PHI[g].denominator for g in rn(1,G_GROUPS))and all(Fr(*phig(q))==1-spx(q)["lo"]for q in(MQ1[1],MQ1[G_GROUPS//2],MQ1[-2],Q1GD)),"PHI[1] %.4f PHI[16] %.4f PHI[31] %.4f"%(float(PHI[1]),float(PHI[16]),float(PHI[31])));HMX=CAP_H;check("c33",HMX==CAP_H<2**59,"2^%.4f"%lg(HMX));PAIRM=KB_C132+Fr(1,r);PAIRK2=KB(D1+D2+64)+Fr(1,r);PJ_MEAN=Fr(JPATH+JPATH2,3)+JPATH*(2*D_WJ+Fr(1,r-1));check("c34",VISIT==4);SG=[NODES*LEAVES*q*Q2 for q in Q1G];WG,SAMEX,WGJ,WGK2,WGL,WGR,WGG=[],[],[],[],[],[],[]
for _g,_q in en(Q1G):_s=SG[_g];_leaf=NODES*_q*(CHK*(FB_FIX+COIN_SCALAR+41+INV+INV_CALL+KEY_FIX)+BATCH_SETUP)+NODES*y_group(_q)+cd(NODES*_q*(FWD_INC+BWD_INC+KEY_INC-EMIT_INC)*Q2,1);_route=_s*(EMIT*K+FLUSH+DONE_C1)//K+_s*H2S_ARR+_s//K*H2S_FULL+NODES*(FINAL*F1+FIN_FIX)+_s*BLOCKH//K2B+N_REL*(FINALH*F2+FIN_FIX)+N_HOME*EXPECTH+N_REL*(EXPECT+RDONE)+NODES*(INIT3_BKT*F1+1);_gp=NODES*LEAVES*MQ1[_g]*Q2;_samex=N_CLUST*_q*CLUST_BATCH*(CLUST_BATCH-1)//2;_main=_s*_gp+_s*(_s-1)//2-_samex;_evis=_main*PAIRM;_gfix=_s*HCHG_NEW+ALLN*STOP_CHECK;WGK2+=[HCHG_VER*_main*PAIRK2];_group=_gfix+VISIT*_evis+WGK2[-1];_jkey=_s*PJ_MEAN+_s//3*SJOIN;WGJ+=[_jkey];WGL+=[_leaf];WGR+=[_route];WGG+=[_gfix+VISIT*_evis];WG+=[_leaf+_route+_group+_jkey+EXCH_ROW];SAMEX+=[_samex]
check("c35",sum(Q1G)==Q1GD and min(Q1G)>0 and all(s%K==0 for s in SG));check("c36",NODES*2**50<=2**116 and NODES<2**67 and max(Q1G)*LEAVES*Q2//K+1<2**50);ZERO_INIT=N_HOME*(2**HB//ZW*ZERO_PASS+1);check("c37",ZERO_PASS==130 and 2**HB%ZW==0 and(2**64-2**55)%2048==0);E_PRE=ALLN*NODE_INIT+E_MINE+NODES*(3*cd(Q1GD,M_G)+1)*COIN_SCALAR+NODES*MX*XSETUP+NODES*6*CH*ENTRY_X+NODES*HEADER+N_REL*(INIT3_BKT*F2+1+TINIT_ENTRY*F2+1)+NODES*(INIT3_BKT*F1+2)+ZERO_INIT+HR_INIT+NODES*XTBLD_ALL+NODES*TBLD;_SIPT=240*RHO;PCFA=NODES*(3*Q1GD+2*NBATCH_MINE+1)*Fr(1,2**254);P_EXC=P_X_ORBIT_ALL_POOLS+P_BAD_X+D_XMINE_ALL+P_MINE_ZERO_ALL+PCFA;W_CORE_DET=sum(WGL[_g]+WGR[_g]+SG[_g]*HCHG_NEW+ALLN*STOP_CHECK+EXCH_ROW for _g in rn(G_GROUPS))+(SG_FULL//3)*SJOIN;W_DET_MAX=(W_CORE_DET+NODES*Q1GD*(B_CAP//CH+1)*FIX_BATCH*max(FIX_EXTRA.values())+SG_FULL*JPATH+SG_FULL*CAP_H*(VISIT+HCHG_VER+HCHG_LS1)+N_HOME*CAND_CAP*(RESOLVE+2*REQ_SERVE));cj=P_EXC*W_DET_MAX;EXC_WORK=cd(cj.numerator,cj.denominator);c4=2*K_U_CLUST*SG_FULL+SAMEX_FAC*SG_FULL*(CLUST_BATCH-1)+SAMEX_FAC*PARTNER_PAIRS+CLUST_BATCH*CLUST_BATCH+96*K_U_CLUST*K_U_CLUST;ELIST=RHO*(Fr(SG_FULL)*SG_FULL*KB_M+c4*KB_M)+_SIPT*KB_M*SG_FULL+72;EK2=RHO*(c4*KB64)+_SIPT*KB64*SG_FULL+72;PAIR_SAMEX=Fr(CLUST_BATCH-K_U_CLUST,CLUST_BATCH-1)*Fr(r,r-2)*W_ELL*PAIRM+Fr(96*(K_U_CLUST-1),CLUST_BATCH-1)*(KB_BUCKET+Fr(1,r));SXK2=PAIR_SAMEX*sum(SAMEX);SXV=PAIR_SAMEX*VISIT*sum(SAMEX)+EXC_WORK;BRACKET=RHO*(Fr(SG_FULL)*SG_FULL*KB(D1+D2)+c4)+_SIPT*SG_FULL+72;DEFV=VISIT*D_BATCH*(HMX*SG_FULL+sum(SAMEX));DEFL_PAIRS=2*D_BATCH*BRACKET+D_BATCH*HMX*SG_FULL;DEFK2=HCHG_VER*DEFL_PAIRS;DEFL=HCHG_LS1*DEFL_PAIRS;DEFJ=JPATH*D_BATCH*SG_FULL;E_FIX_WORK=NODES*Q1GD*E_FIX*max(FIX_EXTRA.values());E_CAND_WORK=(ELIST+DEFL_PAIRS)*RESOLVE;E_REQ_WORK=(ELIST+DEFL_PAIRS)*(2*REQ_SERVE);E_SHARED_REST=E_PRE+NODES*GEN_MUL+(2+2*G_GROUPS)*ALLN*MSG+NBAR_ALLOW*BARRIER+ALLN*PHASE_CTRL_G+OUTPUT;E_RUN=E_SHARED_REST+E_FIX_WORK+E_CAND_WORK+E_REQ_WORK+sum(PHI[_g]*WG[_g]for _g in rn(G_GROUPS))+SXV+EK2*HCHG_VER+ELIST*HCHG_LS1+DEFV+DEFK2+DEFL+DEFJ+SXK2*HCHG_VER;E_TOTAL=GATE_INSNS+P_GATE*E_RUN;W_E=E_TOTAL/C
def claim15_of(w):
 getcontext().prec=80;start=(Decimal(w.numerator).ln()-Decimal(w.denominator).ln())/Decimal(2).ln()*10**15;c=Fr(int(start.to_integral_value(rounding=ROUND_CEILING)),10**15)
 for _ in rn(400):
  if P2(c)[0]>=w and P2(c-Fr(1,10**15))[1]<w:return c
  c+=Fr(1,10**15)if P2(c)[0]<w else-Fr(1,10**15)
 raise AssertionError("E-claim not certified")
CLAIM_E=claim15_of(W_E);JROW_E=sum(PHI[_g]*WGJ[_g]for _g in rn(G_GROUPS))+DEFJ;K2ROW_E=sum(PHI[_g]*WGK2[_g]for _g in rn(G_GROUPS))+SXK2*HCHG_VER+EK2*HCHG_VER;PER_LEAF_OPEN=E_RUN/SG_FULL;PLE=PER_LEAF_OPEN;LOW,UP_BASE=2**55,2**64-2**55
def valid_va(a):return 0<=a<2**64 and(a<2**55 or a>=2**64-2**55)
TB=2**40;XTK_BYTES=3*10240*TB;AXB=4096*TB;AXYX_HI=16384*TB;FWB=19456*TB;XTB=2**64-32768*TB;AXY_END=24576*TB;XT_END=XTB+XTK_BYTES;CN_BASE=84*TB;CMB=90*TB;H1B=24576*TB;H1U_BASE=28672*TB;H2U_BASE=2**64-32768*TB;_KS=[open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"kernels",f)).read()for f in("leaf.s","secp_fe.s")];check("c38",set(re.findall(r"(?:LK1|LKP6_PAIR) (0x[0-9a-f]+)",_KS[0]))=={hex(FWB),hex(FWB+2560*TB)}and[tuple(int(x,16)for x in m.split(","))for m in re.findall(r"OFFS7 (?:\w+,){10}(0x[0-9a-f,x]+)",_KS[1])]==[tuple(XTB+k*10240*TB+d*2**50 for d in(0,2,8,4,6,9))for k in(0,2)]and[tuple(int(x,0)for x in m.split(","))for m in re.findall(r"\nOFF5 (?:\w+,){11}([\w,]+)",_KS[1])]==[tuple(d<<48 for d in(0,2,4,6,8)),tuple(XTB+10240*TB+d*2**50 for d in(0,2,4,6,8))]and XTB%2**50==0 and BT_BASE==2**64-2**51 and AXB%2**49==0 and AXYX_HI%2**48==0 and XTB%2**48==0 and valid_va(XT_END-1)and all(valid_va(BT_BASE|((w<<8)%2**64))for w in(0,1,2**43-1,2**64-1)));YRSL=2**31;YR_BASE=32*TB+YRSL;EB_END=42*TB+2**37;EB_START=EB_END-6*CH*512;YQ_BASE=40*TB;T34=34*TB;LINK_GEN,LINK_GEN_MSG,PORTY_GEN=T34,T34+0x1000,T34+0xB000;G34_END=38*TB;REQ_LO=96*TB-64;REQ_END=REQ_LO+REQ_CAP*64+64;LHR=2**55-2**40;LHRM=LHR+0x8000;NEXT_LO=2**64-24576*TB;NEXT_BASE=NEXT_LO+(R2_UP_TIB//2)*TB;R2UB=R2_UP_TIB*TB;R2LB=CAP_H*32-R2UB;HCC_BYTES=2**20+8*2048+4096+CAND_CAP*528;G34_CONTENT=2**20+8*2048+MX*6*XT_SLOTS*64+2048+MX*64+(2*MX+2*cd(Q1GD,M_G)+Q1GD+1)*32+2**30+YB*2304+1024*128+512+2048;REGISTRY=[("G","LC1",0,F1*8,"init1 init3"),("G","YR",YR_BASE-64,YR_BASE-64+(YB+1)*YRSL,"yround"),("G","EB",EB_START,EB_END,"setup:ENTRY_X fwd"),("G","YQ",YQ_BASE-64,YQ_BASE+2**36+2**31,"ydlv"),("G","LG",LINK_GEN,LINK_GEN+4096,"MMIO"),("G","LGM",LINK_GEN_MSG,LINK_GEN_MSG+4096,"MMIO"),("G","PY",PORTY_GEN,PORTY_GEN+1024,"MMIO"),("G","G34",T34+0x10000,G34_END,"setup:NODE_INIT setup:XSETUP setup:COIN_SCALAR ydlv0 ydsp wcnt dlvq psel setup:CHAIN_SET yround"),("G","CN",CN_BASE,CN_BASE+64*Q1GD,"setup:BATCH_SETUP"),("G","CNM",CMB,CMB+64*NBATCH_MINE,"setup:COIN_SCALAR"),("G","REQ",REQ_LO,REQ_END,"setup:NODE_INIT dlvq"),("G","AXY",AXB,AXY_END,"tbld"),("G","XT6",XTB,XT_END,"tbldk"),("G","H1L",H1B,H1B+4096*TB,"init1"),("G","H1U",H1U_BASE,H1U_BASE+4096*TB,"init1"),("G","BT",BT_BASE,2**64,"tbldf"),("R","RCS",0x10000,0x10000+2**20+8*2048,"setup:NODE_INIT h2s0"),("R","LC2",2**45,2**45+F2*8,"init3"),("R","RC2",3*2**45,3*2**45+F2*8,"tinit"),("R","E2P",2**47,2**47+2*F2*8,"analytic"),("R","DPW",2**48-4096,2**48,"setup:NODE_INIT h2s0 h2s"),("R","RCB",2**48,2**48+0x10000+RELAY_CAP*64,"setup:NODE_INIT setup:EXCH_ROW dlvq"),("R","LR",LHR,LHR+0x8000,"MMIO"),("R","LRM",LHRM,LHRM+0x8000,"MMIO"),("R","RCT",2**55-8192,2**55-4096,"setup:NODE_INIT dlvq"),("R","H2U",H2U_BASE,2**64,"init3"),("H","R2L",0,R2LB,"dlvh"),("H","LH",LHR,LHR+0x8000,"MMIO"),("H","LHM",LHRM,LHRM+0x8000,"MMIO"),("H","HCC",LHR+0x10000,LHR+0x10000+HCC_BYTES,"setup:NODE_INIT init3 dlvh"),("H","HCS",2**55-16384,2**55,"setup:NODE_INIT dlvh0"),("H","HH",UP_BASE,UP_BASE+2**HB*8,"zero16384 hhclr dlvh"),("H","NX",NEXT_LO,NEXT_LO+CAP_H*16,"dlvh"),("H","R2U",2**64-R2UB,2**64,"dlvh")];cl=os.path.dirname(os.path.abspath(__file__));d3=set(re.findall(r"\.globl\s+(\w+)","".join(open(os.path.join(cl,"kernels",f)).read()for f in("leaf.s","route.s","secp_fe.s"))))
def _bld_ok(b):t=b.split();return bool(t)and all(x in d3 or x in("MMIO","analytic")or(x[:6]=="setup:"and isinstance(globals().get(x[6:]),int))for x in t)and(len(t)==1 or not{"MMIO","analytic"}&set(t))
def disjoint_va(I):o=sorted(I);return(all(0<=a<b<=2**64 and valid_va(a)and valid_va(b-1)and((a<LOW)==(b-1<LOW))for a,b in o)and all(o[i][1]<=o[i+1][0]for i in rn(len(o)-1)))
RN=dc(G=NODES,R=N_REL,H=N_HOME);ROLE_B={k:sum(hi-lo for r_,_,lo,hi,_ in REGISTRY if r_==k)for k in RN};PG,BLK=2**16,2**28
def _pt(k):
 iv=sorted((lo//PG*PG,cd(hi,PG)*PG,b=="MMIO",hi-lo)for r_,_,lo,hi,b in REGISTRY if r_==k);pg={};ok=True;l1=set();l2=set();l3=set();I=[];J=[]
 for lo,hi,m,_ in iv:
  for a in(lo,hi-PG):ok&=pg.setdefault(a,m)==m
  t=J if m else I
  if t and lo<=t[-1][1]:t[-1]=(t[-1][0],max(t[-1][1],hi))
  else:t+=[(lo,hi)]
  l1.update(rn(lo>>52,(hi-1>>52)+1));l2.update(rn(lo>>40,(hi-1>>40)+1))
  if m:l3.update(rn(lo>>28,(hi-1>>28)+1))
  else:
   if lo%BLK:l3.add(lo>>28)
   if hi%BLK:l3.add(hi-1>>28)
 ok&=all(b1<=a2 or b2<=a1 for a1,b1 in I for a2,b2 in J);T=2+len(l1)+len(l2)+len(l3);o=PG*T+sum(b-a for a,b in I)+sum(b-a for a,b in J);return ok,T,o-ROLE_B[k],o
PT={k:_pt(k)for k in RN};PROG_DELIV,NODE_STATE,MMU_ON,PT_DESC=(2**18//32)*2*MSG,2**12,64,4;NI_ITEM={k:PROG_DELIV+NODE_STATE+MMU_ON+8*2**10+PT_DESC*4096*PT[k][1]+2**21 for k in RN};check("c39",all(PT[k][0]for k in RN)and PT["G"][1]==65488 and PT["R"][1]==65307 and PT["H"][1]==65544 and max(NI_ITEM.values())<=NODE_INIT and sum(RN[k]*NI_ITEM[k]for k in RN)+NODES*F1<=ALLN*NODE_INIT,"tables "+"/".join(str(PT[k][1])for k in"GRH"));PEAK=sum(RN[k]*PT[k][3]for k in RN);check("c40",len(REGISTRY)==34 and[sum(r_==k for r_,_,_,_,_ in REGISTRY)for k in"GRH"]==[16,10,8]and all(lo<hi and _bld_ok(b)for _,_,lo,hi,b in REGISTRY)and all(disjoint_va([(lo,hi)for r_,_,lo,hi,_ in REGISTRY if r_==k])and ROLE_B[k]<PT[k][3]<2**56 for k in RN));check("c41",YR_BASE-64+(YB+1)*YRSL<=33*TB and YROUNDS*YPC*512+64<=YRSL and T34<YQ_BASE-2**37 and YQ_BASE+2**36+2**31<=EB_START and EB_END<=CN_BASE and G34_END-T34-0x10000>=G34_CONTENT and CN_BASE+64*Q1GD<=CMB and CMB+64*NBATCH_MINE<=REQ_LO and REQ_END<=AXB and AXYX_HI+12*2**48==FWB and FWB+2*2560*TB==AXY_END<=H1B and R2LB<LHR and UP_BASE+2**HB*8==NEXT_LO and NEXT_LO+CAP_H*16<2**64-R2UB and HCC_BYTES<2**40-0x10000-16384);check("c42",Q1GD<=2**37 and Q2<2**24 and YROUNDS*YPC==L_CHAIN and PTS_MINE_CHAIN<2**18 and NBATCH_MINE<2**44)
def laa(d,x,o):return P2(d)[0]>=x and P2(d-Fr(1,o))[1]<x
def least_grid_up(x,o):
 k=math.ceil(lg(x)*o)
 for _ in rn(20):
  d=Fr(k,o)
  if laa(d,x,o):return d
  k+=1 if P2(d)[0]<x else-1
 raise AssertionError("grid bound not certified")
DECLARED=dc(time_log2=CLAIM_E,memory_log2_bytes=least_grid_up(PEAK,10**4),preprocessing_log2=least_grid_up(Fr(E_PRE,C),10**4),success=THRESH,per_leaf=Fr(cd(PLE.numerator*10**4,PLE.denominator),10**4));GEN_ITEM=256*(DBL+MADD_)+2*(INV+3*MULX+2*SQRX+10*ADDX)+64;GEN_TEST,GEN_CTRL,NEGT,RED=256*2*13+64,256*6+64,30,9;X_ENTRY=MADD_+3*SQRX+11*MULX+2*ADDX+2*NEGT+17*RED+32;check("c43",GEN_ITEM+GEN_TEST+GEN_CTRL+X_ENTRY==GEN_MUL and GEN_MUL+4==GEN_MULX);ENTRY_W=[14+2*(j&1)for j in rn(MX)];EIT=sum(2*cd(w,2)for w in ENTRY_W)+3+MX//2;TEMPLATE_ITEM=3*MULX+SQRX+24+7+20;check("c44",INV+3*MULX+2*SQRX+10*ADDX+4<CHAIN_ADD and 56*(2*MULX+ADDX+4)+TEMPLATE_ITEM+X_ENTRY<XSETUP and EIT<=MX*ENTRY_X,"entries %d/8"%EIT);BARRIER_MSGS=2*(N_HOME+NODES+2*(N_REL-1));BATCH_REST=11+3*9+2+3+6+3+MSG+4+SWAP+20+8+7+6+2**13;PH_STEP=8+8+2+15;BARRIER_MSG=7+DONE+DLV+RCV+WAKE+CTX;BARRIER_ITEM=ALLN*(3*13+1)+BARRIER_MSGS*BARRIER_MSG;check("c45",COIN_SCALAR+cd(4*YROUNDS,YB)+BATCH_REST<=BATCH_SETUP and G_GROUPS*(12*PH_STEP+16+4)<PHASE_CTRL_G and BARRIER>=7*BARRIER_ITEM and NBAR_COUNT<=NBAR_ALLOW and 9+MSG_LINK<=MSG);DLVQ,DLVQ_DROP,DLVQ_CARRY,REQH_C,FWDQ,SENDR,FWDR=14,11,4,17,11,10,15;PH4_HOPS={k:v+DONE+RCV+2*WAKE+CTX for k,v in dc(req_home_relay=cd(REQH_C,2)+DLVQ+DLVQ_CARRY,req_relay_gen=FWDQ+3+DLVQ+DLVQ_CARRY,rep_gen_relay=SENDR+DLVQ+DLVQ_CARRY,rep_relay_home=FWDR+DLV).items()};_h=list(PH4_HOPS.values());check("c46",max(_h)<=MSG and RESOLVE==2*CAND_ITEM+1 and PROG_DELIV<NODE_INIT);COIN_ITEM=2*(4*2+5+4)+16+2+1+2;P_COIN_FAIL=Fr(2**256-r+1,2**256)**2;INVC_ITEM,STOP_ITEM,HEADER_ITEM,REQ_ITEM,PH4_FIX,OUT_ITEM=2+1+6+6,5,17*4+9,16+12+600+16,2*EXPECT+DONE_C+CTX+32+24,16;check("c47",COIN_ITEM<=COIN_SCALAR and INVC_ITEM<=INV_CALL and STOP_ITEM<=STOP_CHECK and HEADER_ITEM<=HEADER and REQ_ITEM+GEN_MUL<=REQ_SERVE and 4*PH4_FIX*ALLN<=EXCH_ROW and OUT_ITEM<=OUTPUT and P_COIN_FAIL<Fr(1,2**254));print("q1 %d  q2 %d  S 2^%.6f  mu %.9f  lo_run %.12f  p_gate %.9f  success >= %.12f"%(Q1GD,Q2,lg(SG_FULL),float(TFULL["mu"]),math.floor(LO_RUN*10**12)/10**12,float(P_GATE),math.floor(GATED_LO*10**12)/10**12));print("E_mine 2^%.6f  E_run 2^%.6f insns (2^%.6f units)  claim %d/10^15  per leaf %.4f  mem 2^%.7f  pre 2^%.7f  grids %s"%(lg(E_MINE),lg(E_RUN),lg(Fr(E_RUN,C)),CLAIM_E*10**15,float(PLE),lg(PEAK),lg(Fr(E_PRE,C)),"/".join(format(float(least_grid_up(ROLE_B[k],10**4)),".4f")for k in"GRH")))
def _ps(a):return sum(PHI[g]*a[g]for g in rn(G_GROUPS))
LEDGER=[ALLN*NODE_INIT,NODES*GEN_MUL,(2+2*G_GROUPS)*ALLN*MSG,E_MINE,NODES*(3*cd(Q1GD,M_G)+1)*COIN_SCALAR,NODES*MX*XSETUP,NODES*6*CH*ENTRY_X+NODES*HEADER,NODES*XTBLD_ALL,NODES*TBLD,N_REL*((INIT3_BKT+TINIT_ENTRY)*F2+2)+NODES*(INIT3_BKT*F1+2),ZERO_INIT,HR_INIT,_ps(WGL),E_FIX_WORK,JROW_E-DEFJ,_ps(WGR),_ps(WGG),SXV,K2ROW_E,ELIST*HCHG_LS1,DEFV+DEFK2+DEFL+DEFJ,E_CAND_WORK,E_REQ_WORK,sum(PHI)*EXCH_ROW,NBAR_ALLOW*BARRIER,ALLN*PHASE_CTRL_G,OUTPUT];LEDGER_CELLS=[format(float(least_grid_up(v,10**4)),".4f")for v in LEDGER];E_RUN_6=least_grid_up(E_RUN,10**6);E_RUN_U6=least_grid_up(Fr(E_RUN,C),10**6);E6=least_grid_up(E_TOTAL,10**6);print("E_total %.6f  E_RUN_6 %.6f  E_RUN_U6 %.6f  tables %s"%(float(E6),float(E_RUN_6),float(E_RUN_U6),"/".join(str(PT[k][1])for k in"GRH")));check("c48",len(LEDGER)==27 and sum(LEDGER)==E_RUN and P2(CLAIM_E)[0]>=Fr(E_TOTAL,C)and P2(E_RUN_6)[0]>=E_RUN and P2(E_RUN_U6)[0]>=Fr(E_RUN,C)and P2(E6)[0]>=E_TOTAL,"cells "+" ".join(LEDGER_CELLS));STATED=dc(time_log2=Fr(123293102848764316,1000000000000000),memory_log2_bytes=Fr(1379998,10**4),preprocessing_log2=Fr(1163558,10**4),success=Fr(39,100),per_leaf=Fr(582232,10**4));d8=open(os.path.join(cl,"proof.md")).read();d7=open(os.path.join(cl,"success.md")).read();cq=open(os.path.join(cl,"submission.json")).read();d6=json.loads(cq)["claims"][0]["score"]
def cnt(text,o):return len(re.findall(r"(?<![0-9])"+re.escape(o)+r"(?![0-9])",text))
EC=[("123.293102848764316",2,0,2),("123293102848764316",1,0,0),("1000000000000000",1,0,0),("137.9998",1,0,1),("116.3558",1,0,1),("58.2232",1,0,1),("0.39",2,2,2),("0.390195572",1,0,1),("0.999498785848",1,1,1),("8.462363218",0,1,0),("134.616618",1,0,0),("124.650834",1,0,1),("133.258888",1,0,0),("128.753096",1,0,1),("1,313,229,200,424,960",1,0,1),("103,000,000,000",1,1,1),("12,576,960",1,0,0),("55.9984",1,0,0),("55.9944",1,0,0),("55.9998",1,0,0),("2048",4,0,0),("10752",3,0,1),("132776672646451629524037084662993969194",1,0,1),("1,185,778,795,529,196,871,059,725,352,960",1,0,0),("8,043,004,160",1,0,0)]+[(c,2 if c=="124.2717"else 1,0,1 if c=="124.2717"else 0)for c in LEDGER_CELLS];d9=[k for(k,p,s,j)in EC if(cnt(d8,k),cnt(d7,k),cnt(cq,k))!=(p,s,j)];check("c49",not d9,"count mismatch: "+",".join(d9));check("c50",set(LEDGER_CELLS)<={k for(k,p,s,j)in EC},"unenforced");check("c51",Fr(int(d6["numerator"]),int(d6["denominator"]))==STATED["time_log2"]==DECLARED["time_log2"]and math.gcd(int(d6["numerator"]),int(d6["denominator"]))==1);RATES={"2^-50.51":(D_WJ,1,1,0),"2^-8.20":(D_LAB,1,1,0),"2^-7.65":(3*D_WJ+D_LAB_SYM,1,0,0),"2^-126":(D_BATCH,1,1,0),"2^-26":(D_DROP,0,2,0),"2^-82":(TFULL["cs"],0,2,0),"2^-18.8":(TFULL["d_f"],0,2,0),"2^-36.6":(TFULL["lam_int"],0,2,0),"2^-33.6":(TFULL["d_p2"],0,2,0),"2^-15.8":(max(TFULL["d_cl"],TFULL["d_rq"],TFULL["d_rl"]),0,6,0),"2^-21.7":(TFULL["g_def"],0,1,0),"2^-19":(max(TFULL["glob"],P_MINE_CAPS),0,3,0),"2^-28":(P_FAIL_ALL_POOLS,0,2,0),"2^-39":(P_SAME_A_ALL,0,1,0),"2^-140":(max(P_BAD_X,D_XMINE_ALL,P_MINE_ZERO_ALL),1,4,0),"2^-139":(P_EXC,0,1,0),"2^-202":(P_X_ORBIT_ALL_POOLS,0,1,0),"2^-254":(P_COIN_FAIL,1,1,0),"2^-99":(320*K_U_CLUST*K_U_CLUST*((POINTWISE_X_MAX/MX)**4)*KB_M,0,1,0)};ANALYTIC={"2^-200","2^-30","2^-187","2^-190"}
def rcnt(text,o):return len(re.findall(re.escape(o)+r"(?!\.?[0-9])",text))
da=[set(re.findall(r"2\^-[0-9.]*[0-9]",t))for t in(d8,d7,cq)];cx=[k for k,(v,p_,s_,j_)in RATES.items()if not(v<P2(-Fr(k[3:]))[0]and(rcnt(d8,k),rcnt(d7,k),rcnt(cq,k))==(p_,s_,j_))];check("c52",not cx and set().union(*da)<=set(RATES)|ANALYTIC,"bad: "+",".join(cx+sorted(set().union(*da)-set(RATES)-ANALYTIC)));check("c53",K_GATE==132776672646451629524037084662993969194 and CAP_H==1313229200424960 and PSEL_POOL==262160 and DEF_MINE==116754189 and EXC_WORK==5941171950738 and E_RUN_6==Fr(134616618,10**6)and E_RUN_U6==Fr(124650834,10**6)and E6==Fr(133258888,10**6)and STATED==DECLARED);print("FAIL: "+", ".join(FAILED)if FAILED else"PASS")
if __name__=="__main__":sys.exit(1 if FAILED else 0)

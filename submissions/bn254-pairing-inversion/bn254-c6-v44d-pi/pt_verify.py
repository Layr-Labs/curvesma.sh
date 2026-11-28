D3='chain_table_E9'
D2='chain_table_E8'
D1='global'
D0='orbits_with_runs'
C_='flush_E9'
Cz='ok_runs'
Cy='exact_ball_selfcheck_bad'
Cx='orbits'
Cw='seconds'
Cv='mnemonics_executed'
Cu='nmaxpos'
Ct='geometry_matches_file'
Cs='pkg_counter_agrees'
Cr='l3_counted'
Cq='/kernels_v32.s'
Cp=' b .Lp_top2\n .endif\n .endm\n .macro LBLT'
Co='tails_are_the_segments'
Cn='^\\s+[a-z][a-z0-9.]*\\b'
Cm='region'
Cl='pt_tb_regions'
Ck='pt_orbit_check'
Cj='CHK_K1'
Ci='e3_build_bin'
Ch='slow_pairs'
Cg='dev_ex'
Cf='pt_empty_init'
Ce='first_bad'
Cd='handoff3c'
Cc='handoff3b'
Cb='handoff3a'
Ca='null_nl'
CZ='.Lp_top2'
CY='FR_3DA0'
CX='FR_DA2'
CW='FR_DA1'
CV='famof_tab'
CU='OFFKEY'
CT='CHAFAM'
CS='.ifndef'
CR='staddb'
CQ=ValueError
CB='pop'
CA='t2_chunked'
C9='prices'
C8='cert'
C7='caught'
C6='t3w'
C5='t2f'
C4='removed_old_handoff_lines'
C3='.Lmb_t'
C2='cases'
C1='G_SCR'
C0='O_K1'
B_='rows_skip'
Bz='rows_emitting'
By='dyn'
Bx='hand'
Bw='l3_exit'
Bv='node'
Bu='last'
Bt='FR_U3'
Bs='cha_null'
Br='NULLSLOT'
Bq='KKEY'
Bp='.ifdef'
Bo='.rept'
Bn='.macro'
Bm='neg'
Bl='csneg'
Bk='ands'
Bj='cmp'
Bi='blr'
Bh='tbz'
Bg='cbz'
Bf='wzr'
Be='xzr'
Bd='little'
Bc='empty_nodes'
Bb='unweighted_size_biased_sample'
Ba='CW3'
BM='fallback'
BL='/kernels_pt.s'
BK='t6'
BJ='flushes'
BI='leaves'
BH='fault'
BG='PBITS'
BF='OFF'
BE='REGBYTES'
BD='//@end'
BC='//@seg'
BB='\n'
BA='ldp'
B9='csinc'
B8='tbnz'
B7='uxtw'
B6='cands_P'
B5='rows_P'
B4='rows_exec'
B3='BMS'
B2='c4'
B1='GW'
B0='EK'
A_='C2'
Az=next
A=Exception
Ay=list
An='t3p'
Am='count_dev'
Al='row_disp'
Ak='node_L2'
Aj='.if'
Ai='//'
Ah='sxtw'
Ag='asr'
Af='lsr'
Ae='lam'
Ad='CW'
AW='PRICE_ROW'
AV='PRICE_CAND'
AU='row'
AT='null_last'
AS='NF3'
AR='SPANK'
AQ='LASTADD'
AP='lsl'
AO=tuple
AH='t3'
AG='.Ltop_3'
AF='PRICE_HAND'
AE='EBITS'
AD='over'
AC='miss'
AB='c3'
AA='GX'
A2='t1'
A1='pt_l3_next'
A0='NF'
z='NMAXPOS'
y=False
x='rows'
w=abs
v='nodes'
t='maxover_ppm'
u='x'
s='cands'
o='T'
r='AE'
q='DB'
p='DA'
n='gate'
i='per_orbit'
h=enumerate
f='E8'
e='bad'
d='ret'
c='MU'
b=divmod
a=sorted
Y='#'
X=set
U='E'
T=str
S=min
Q=bool
P=max
N=all
L='t2'
I=round
H=sum
K=int
G=None
E=dict
D='ok'
C=len
B=range
import random,re,time,stor_bound as AI
from fractions import Fraction as Fr
from math import isqrt,log2
import v20_lat_b2 as W,stor_bound as SB
M=(1<<64)-1
j=(1<<32)-1
Q32=1<<32
V=E(CW=24,PB=28,CW3=22,EB=26)
PT_D1=39*(1<<64)//10**6+1
PT_GATES=25,35,50
PT_G=1
def pt_cert(d,lam,T,gate=G,geo=G,isq=isqrt):
	Q=isq;L=gate;M=geo or V;Z,o,R,p=M[Ad],M['PB'],M[Ba],M['EB'];L=L if L is not G else PT_GATES[PT_G]*(1<<64)//10000;a=lambda n:Q(n-1)+1;I=lambda a,b:-(-a//b);b=(1<<64)-1;U=1<<32;C=12;W=[S(b,(d[aa+1]<<64)//(d[aa]*T))for aa in B(C)];q=[((d[ab]*T<<38)+d[ab+1]-1)//d[ab+1]for ab in B(C)];X=[0]*C
	for A in B(C-1,-1,-1):X[A]=Q((d[A]*T+d[A+1]-1)//d[A+1])+1+H((A+1)//2 for A in X[A+1:])
	J=[H(X[ac+1:])+1 for ac in B(C)];c=lambda k,D:(P(D-J[k],0)>>13)**2*W[k]>>64;Y=E(C1=N(2*w(lam[A][C])<=d[C+1]for A in B(C)for C in B(A)),C2=N(1<<39<=W[A]<b and 0<d[A+1]<d[A]*T for A in B(3)),C3=c(1,120*U+1)>1<<38 and c(2,120*U+1)>1<<38);K=[0]*C;e=0
	for A in B(C-1,0,-1):K[A]=I(J[A]+Q((q[A]>>26)+1<<52)+1+(1<<25),U)+e;e+=K[A]+1>>1
	f=H(K[3:]);g=f+K[2];r=g+K[1];h=(1<<Z)-1+g;i=(1<<Z)-1+J[0]+r;D=[ad+1 for ad in W[:3]];F=[a(ae)for ae in D];j=I(D[1]*h*h+D[0]*i*i,1<<64);k=a(j);O=8192+J[1];l=(1<<27)+2*O*F[1]+(2*O*F[1]>>30)+I(D[1]*O*O,1<<64)+4+(D[0]>>12)+(D[0]>>37)+2*F[0]+(2*F[0]>>12)+5;s=l+j+(k<<33)+I(k*l,1<<32)+1;t=((1<<p)+(1<<o)<<26)+(F[2]<<R+1)+(F[2]>>38-R)+(D[2]>>64-2*R)+(1<<26)+3;u=I(2*(8192+J[2]+f)*F[2]*1003,1000)+(1<<27)+2;m=PT_D1+u+t+s;Y[n]=m<=L;return N(Y.values()),(m,1<<64),Y
CERT=[]
def CC(orbs,gates=PT_GATES):
	A=gates;M=W.M.OP['p_lat'];B=[];N={};D=E(C1=0,C2=0,C3=0);G={aa:0 for aa in A};del CERT[:]
	for(U,d,V,X,Y)in orbs:
		Z,(J,e),F=pt_cert(V,X,Y,gate=1<<64);B.append(J);CERT.append((Z,J))
		for O in D:D[O]+=not F[O]
		b=K((log2(U)-M)//1.5);Q=N.setdefault(b,[0]*(1+C(A)));Q[0]+=1
		for(c,R)in h(A):T=not(F['C1']and F[A_]and F['C3'])or J>R*(1<<64)//10000;G[R]+=T;Q[1+c]+=T
	L=lambda e:I(e*100/(1<<64),5);return E(n=C(B),eps_min_pct=L(S(B)),eps_mean_pct=L(H(B)//C(B)),eps_max_pct=L(P(B)),fallbacks_at_gate_pct={ab/100:G[ab]for ab in A},gate_pct=A[PT_G]/100,failed_checks=D,bands_1p5bit=[(M+1.5*ad,ac[0],ac[1:])for(ad,ac)in a(N.items())],ok=D['C1']+D[A_]+D['C3']==0 and G[A[PT_G]]==0)
def pt_orb(o,full=1):
	O,G,C,N,I=o;D=12;A=E(p=O,B=G,d=C,lam=N,T=I);A[c]=[[(N[aa][ai]<<32)//C[ai+1]if aa>ai else 0 for ai in B(D)]for aa in B(D)];A[B0]=[S(M,(C[ab+1]<<64)//(C[ab]*I))for ab in B(D)];A[B1]=[((C[ac]*I<<38)+C[ac+1]-1)//C[ac+1]for ac in B(D)];J=[0]*D
	for F in B(D-1,-1,-1):J[F]=isqrt((C[F]*I+C[F+1]-1)//C[F+1])+1+H((A+1)//2 for A in J[F+1:])
	A[AA]=[H(J[ad+1:])+1 for ad in B(D)];A[AB]={};A[B2]={};A['c5']={}
	if not full:return A
	A[p]=[16*W._mr_da(G[ae][:6])&M for ae in B(D)];A[q]=[H(G[af][6+A]*W.WB_MR[A]for A in B(6))&M for af in B(D)];A['A0']=BT+16*W.A_OFF_MR&M;A['B0']=W.B_OFF_MR&M;A[r]=[W._st_f3(G[ag])for ag in B(D)];K=W._add_st_f3;L=A[r][0];P=K(L,L);A[B3]=[(1 if W._bad_st_f3(ah)else 0)|(2 if W._bad_st_f3(K(ah,L))else 0)|(4 if W._bad_st_f3(K(ah,P))else 0)for ah in B(729)];return A
def BN(O,k,D):return(P(D-O[AA][k],0)>>13)**2*O[B0][k]>>64
def BO(M):return isqrt(M<<52)+1
def pt_geo(O,f1,bp,mut=0):
	C=mut;F='g',f1,bp,C
	if F in O[AB]:return O[AB][F]
	G,H=V[Ad],V['PB'];R=1<<38-H;T=((bp+1<<H)-1 if C!=2 else bp<<H)if bp<R-1 else 1<<38;I=f1<<G;L=I if C==1 else(f1+1<<G)-1;M=O[AA][0];U=O[B1][0];N=O[c][1][0];A=[]
	for D in(-1,1):
		J=[]
		for E in B(0 if D<0 else 1,-121 if D<0 else 121,D):
			K=E<<32;W=0 if I<=K<=L else S(w(I-K),w(L-K));P=T-BN(O,1,W)
			if P<0:break
			Q=BO((P*U>>64)+1);J.append((E,-N*E-M-Q+Q32-1,-N*E+M+Q-1+(0 if C==1 else 1<<G)))
		A=J[::-1]+A if D<0 else A+J
	O[AB][F]=A;return A
def pt_rows(O,f1,f0,bp,mut=0):
	B=f1,f0,bp,mut
	if B in O[AB]:return O[AB][B]
	C=V[Ad];A=[(aa,ab+(f0<<C)>>32,ac+(f0<<C)>>32)for(aa,ab,ac)in pt_geo(O,f1,bp,mut)];A=[ad for ad in A if ad[1]<=ad[2]];O[AB][B]=A;return A
def pt_e3(O,cell,b3,mut=0):
	E=cell;F=E,b3,mut
	if F in O[B2]:return O[B2][F]
	J,K,L=V[Ba],V['EB'],V['PB'];X=1<<38-L;M=E<<J;Q=(E+1<<J)-1;R=(b3+1<<K)-1 if b3<1<<38-K else 1<<38;A=[]
	for D in(-1,1):
		H=[]
		for T in B(0 if D<0 else 1,-121 if D<0 else 121,D):
			I=T<<32;Y=0 if M<=I<=Q else S(w(M-I),w(Q-I));U=BN(O,2,Y)
			if U>R:break
			H.append((T,P(S(R-U>>L,X-1)-(mut==3),0)))
		A=H[::-1]+A if D<0 else A+H
	assert N(A[B][0]+1==A[B+1][0]for B in B(C(A)-1));W=(A[0][0],[aa for(ab,aa)in A])if A else G;O[B2][F]=W;return W
def pt_node_data(O,x,rem3):C=O[c];A=lambda k:-H(C[A][k]*x[A]for A in B(k+1,12));return E(x=Ay(x),C2=A(2),c1b=A(1),c0b=A(0),rem3=rem3)
def pt_cands(O,nd,mut=0):
	B=mut;A,H,Q=V[Ad],V[Ba],V['EB'];C=O[c];I=nd[A_];R=I>>32;D=pt_e3(O,I>>H&(1<<32-H)-1,nd['rem3']>>Q,B==3 and 3);E=[];S=B if B<3 else 0
	if D is G:return E
	for(T,J)in h(D[1]):F=R+D[0]+T;K=nd['c1b']-C[2][1]*F;U=nd['c0b']-C[2][0]*F;L=K>>32;M=U-C[1][0]*L;N=K>>A&(1<<32-A)-1;P=M>>A&(1<<32-A)-1;E.append((F,(N,P,J),M>>32,L,pt_rows(O,N,P,J,S)))
	return E
def D4(O,nd,mut=0):
	A=X()
	for(C,K,D,E,F)in pt_cands(O,nd,mut):
		for(G,H,I)in F:
			for J in B(H,I+1):A.add((C,E+G,D+J))
	return A
def pt_ball(O,x,e=0):
	D,A,V=O['d'],O[Ae],O[o];E,F,C=D[1],D[2],D[3];W=H(Fr((D[C+1]*x[C]+H(A[B][C]*x[B]for B in B(C+1,12)))**2,D[C]*D[C+1])for C in B(3,12));Y=E*F*C;G=(Fr(V)*(1+Fr(e,1<<64))-W)*Y//1;J=X()
	if G<0:return J
	G=K(G);L=H(A[B][2]*x[B]for B in B(3,12));Q=isqrt(G//E)
	for I in B(-((L+Q)//C),(Q-L)//C+1):
		M=G-(C*I+L)**2*E
		if M<0:continue
		N=A[2][1]*I+H(A[B][1]*x[B]for B in B(3,12));R=isqrt(M//C)
		for P in B(-((N+R)//F),(R-N)//F+1):
			S=M-(F*P+N)**2*C
			if S<0:continue
			T=A[1][0]*P+A[2][0]*I+H(A[B][0]*x[B]for B in B(3,12));U=isqrt(S//(F*C))
			for Z in B(-((T+U)//E),(U-T)//E+1):J.add((I,P,Z))
	return J
def CD(O,x,t):
	C,F=O['d'],O[Ae];D=[t[2],t[1],t[0]]+Ay(x[3:]);E=0
	for A in B(12):E+=Fr((C[A+1]*D[A]+H(F[B][A]*D[B]for B in B(A+1,12)))**2,C[A]*C[A+1])
	return E/O[o]
def pt_ch(O,x,k,rem):
	D=rem;A,R,G=O['d'],O[Ae],O[o];S,T,I=O[c],O[B1],O[AA];E=[]
	if A[k+1]>=G*A[k]:
		J=H(R[A][k]*x[A]for A in B(k+1,12));K=-J//A[k+1];L=G*A[k]*A[k+1];M=D*L
		for C in B(K-1,K+3):
			N=(1<<38)*(A[k+1]*C+J)**2
			if N<=M:E.append((C,-((N-M)//L)))
	else:
		F=-H(S[A][k]*x[A]for A in B(k+1,12));P=BO((D*T[k]>>64)+1)
		for C in B(F-I[k]-P+Q32-1>>32,(F+I[k]+P>>32)+1):
			Q=BN(O,k,w((C<<32)-F))
			if Q<=D:E.append((C,D-Q))
	return E
def CE(O,n,seed,tries=300):
	J=random.Random(seed);D=[];G=X()
	for L in B(tries):
		A=[0]*12;E=1<<38;H=1
		for I in B(11,2,-1):
			F=pt_ch(O,A,I,E)
			if not F:H=0;break
			K,E=J.choices(F,[aa+1 for(ab,aa)in F])[0];A[I]=K
		if H and AO(A)not in G:
			G.add(AO(A));D.append(pt_node_data(O,A,E))
			if C(D)>=n:break
	return D
def BP(r):A=r[B4]/P(r[B5],1);B=r[s]/P(r[B6],1);C=(r['N']-r[U])/P(r[U],1);r[Bb]=E(rows=I(A,5),cands=I(B,5),extra=I(C,5),note='NOT production means (nodes drawn with child weight budget+1, not by their share of the level-3 nodes): the production-weighted factors are derived by pt_pop in part pt');return r
def Ao(O,nodes,e,mut=0,kv=G):
	J=nodes;A=E(nodes=C(J),E=0,N=0,miss=0,over=0,rows=0,cands=0,maxover_ppm=0,empty_nodes=0,cands_P=0,rows_P=0,rows_exec=0)
	for B in J:
		F=pt_ball(O,B[u]);G=D4(O,B,mut);L=pt_ball(O,B[u],e);I=pt_cands(O,B,mut);A[s]+=C(I);A[x]+=H(C(A[4])for A in I);A[B4]+=H(P(1,C(pt_pieces(A[4],kv)))for A in I);A[B6]+=C({aa[0]for aa in F});A[B5]+=C({ab[:2]for ab in F});A[U]+=C(F);A['N']+=C(G);A[AC]+=C(F-G);A[AD]+=C(G-L);A[Bc]+=not G
		for M in G-F:A[t]=P(A[t],K((CD(O,B[u],M)-1)*10**6))
	A[D]=A[AC]==0 and A[AD]==0 and A[U]>0;return BP(A)
def D5(O,nd):
	A=pt_ball(O,nd[u]);J=nd[u];C=0
	if not A:return 0
	D=[aa[0]for aa in A];E=[ab[1]for ab in A];F=[ac[2]for ac in A]
	for G in B(S(D)-1,P(D)+2):
		for H in B(S(E)-1,P(E)+2):
			for I in B(S(F)-1,P(F)+2):C+=((G,H,I)in A)!=(CD(O,J,(G,H,I))<=1)
	return C
R=[0]*40
F=[0,0,0,0]
CODEB=0x600000000000
STOPA=0x5ffffffffff0
A3=31232<<32
class Fault(A):0
class AX(A):0
def Ap(v):v&=M;return v-(1<<64)if v>>63 else v
def l(v):v&=j;return v-(1<<32)if v>>31 else v
class D6:
	def __init__(A):A.reg=[];A.pg={}
	def add(A,lo,hi,fn=G):A.reg.append((lo,hi,fn))
	def chk(B,a,n,w):
		for(C,D,A)in B.reg:
			if C<=a and a+n<=D:
				if w and A:raise Fault('write to table %x'%a)
				return A
		raise Fault(('write %x'if w else'read %x')%a)
	def rd(B,a,n):
		A=a&4095;C=B.pg.get(a>>12)
		if C is not G and A+n<=4096:return K.from_bytes(C[A:A+n],Bd)
		D=B.chk(a,n,0)
		if D:return D(a)
		if C is G and A+n<=4096:return 0
		return B.rd(a,4096-A)|B.rd(a+4096-A,n-4096+A)<<8*(4096-A)
	def wr(D,a,n,v):
		A=a&4095;C=D.pg.get(a>>12);v&=(1<<8*n)-1
		if C is not G and A+n<=4096:C[A:A+n]=v.to_bytes(n,Bd);return
		D.chk(a,n,1)
		while n:A=a&4095;B=S(n,4096-A);C=D.pg.setdefault(a>>12,bytearray(4096));C[A:A+B]=(v&(1<<8*B)-1).to_bytes(B,Bd);v>>=8*B;a+=B;n-=B
MEM=D6()
EM=[]
DR={}
D7={'||':lambda a,b:K(Q(a or b)),'&&':lambda a,b:K(Q(a and b)),'|':lambda a,b:a|b,'^':lambda a,b:a^b,'&':lambda a,b:a&b,'==':lambda a,b:K(a==b),'!=':lambda a,b:K(a!=b),'<':lambda a,b:K(a<b),'>':lambda a,b:K(a>b),'<=':lambda a,b:K(a<=b),'>=':lambda a,b:K(a>=b),'<<':lambda a,b:a<<b,'>>':lambda a,b:a>>b,'+':lambda a,b:a+b,'-':lambda a,b:a-b,'*':lambda a,b:a*b,'/':lambda a,b:a//b,'%':lambda a,b:a%b}
BQ={'||':1,'&&':1,'==':2,'!=':2,'<':2,'>':2,'<=':2,'>=':2,'+':3,'-':3,'|':4,'&':4,'^':4,'*':5,'/':5,'%':5,'<<':5,'>>':5}
CF={}
def O(e,sym):
	B=CF.get(e)
	if B is G:
		D=re.findall('0x[0-9a-fA-F]+|\\d+|[A-Za-z_.$][\\w.$]*|<<|>>|<=|>=|==|!=|&&|\\|\\||[-+*/%&|^~!<>()]',e);A=[0]
		def F():
			B=D[A[0]];A[0]+=1
			if B=='(':C=E(0);A[0]+=1;return C
			if B in'-~!':C=F();return(lambda S:-C(S))if B=='-'else(lambda S:~C(S))if B=='~'else lambda S:K(not C(S))
			if B[0].isdigit():G=K(B,0);return lambda S:G
			return lambda S:S[B]
		def E(m):
			B=F()
			while A[0]<C(D)and D[A[0]]in BQ and BQ[D[A[0]]]>m:G=D[A[0]];A[0]+=1;H=E(BQ[G]);B=(lambda a,b,g:lambda S:g(a(S),b(S)))(B,H,D7[G])
			return B
		B=CF[e]=E(0)
	return B(sym)
def split(r):
	if'['not in r:return[aa.strip()for aa in r.split(',')]if r.strip()else[]
	C=[];A='';D=0
	for B in r:
		D+=(B=='[')-(B==']')
		if B==','and D==0:C.append(A.strip());A=''
		else:A+=B
	if A.strip():C.append(A.strip())
	return C
def J(t,d=0):
	if t in(Be,Bf):return 33 if d else 31
	return 32 if t in('sp','wsp')else K(t[1:])
AJ={'eq':lambda:F[1],'ne':lambda:not F[1],'cs':lambda:F[2],'hs':lambda:F[2],'cc':lambda:not F[2],'lo':lambda:not F[2],'mi':lambda:F[0],'pl':lambda:not F[0],'vs':lambda:F[3],'vc':lambda:not F[3],'hi':lambda:F[2]and not F[1],'ls':lambda:not F[2]or F[1],'ge':lambda:F[0]==F[3],'lt':lambda:F[0]!=F[3],'gt':lambda:not F[1]and F[0]==F[3],'le':lambda:F[1]or F[0]!=F[3],'al':lambda:1}
BR={'eq':'ne','ne':'eq','cs':'cc','cc':'cs','hs':'lo','lo':'hs','mi':'pl','pl':'mi','hi':'ls','ls':'hi','ge':'lt','lt':'ge','gt':'le','le':'gt'}
CG='ldp x0, x6, [x2]','ldr x3, [x6, x5, lsl #3]'
def Aq(o,i,m):
	E=o[i]
	if E[0]==Y:
		F=K(E[1:],0)
		if i+1<C(o):F<<=K(o[i+1].split(Y)[1])
		F&=m;return lambda:F
	B=J(E)
	if i+1<C(o):
		D,A=(o[i+1].split()+['#0'])[:2];A=K(A[1:])
		if D==AP:return lambda:R[B]<<A&m
		if D==Af:return lambda:(R[B]&m)>>A
		if D==Ag:G=l if m==j else Ap;return lambda:G(R[B])>>A&m
		if D==Ah:return lambda:l(R[B])<<A&m
		if D==B7:return lambda:(R[B]&j)<<A&m
	return lambda:R[B]&m
def AY(o,k):
	D=o[k];A=J(D[1:].split(',')[0].split(']')[0].strip());B=[aa.strip()for aa in D[1:D.index(']')].split(',')]
	if C(o)>k+1:
		E=K(o[k+1][1:],0)
		def F():B=R[A];R[A]=B+E&M;return B
		return F
	if C(B)==1:return lambda:R[A]
	if B[1][0]==Y:
		E=K(B[1][1:],0)
		if D.endswith('!'):
			def F():B=R[A]+E&M;R[A]=B;return B
			return F
		return lambda:R[A]+E&M
	G=J(B[1]);H=K(B[2].split(Y)[1])if C(B)>2 and Y in B[2]else 0;I=B[2].split()[0]if C(B)>2 else AP
	if I==Ah:return lambda:R[A]+(l(R[G])<<H)&M
	if I==B7:return lambda:R[A]+((R[G]&j)<<H)&M
	return lambda:R[A]+(R[G]<<H)&M
def D8(t,pc,A):
	V=pc;D,x,p=t.partition(' ');B=split(p)
	if t in CG and A.sink:
		if t==CG[0]:
			def L():
				A=R[2]-BT
				if A<0 or A&15 or A>>4>=W.PROD_MA_MR or R[5]>=W.PROD_MB_MR:raise Fault('emit address')
				EM.append((R[2],R[5]));R[0]=0;R[6]=0
		else:
			def L():R[3]=0
		return L
	if D in('b','bl'):
		def L():
			if D=='bl':R[30]=CODEB+4*(V+1)
			return A.pcof(B[0],V)
		return L
	if D[:2]=='b.':
		T=AJ[D[2:]]
		def L():
			if T():return A.pcof(B[0],V)
		return L
	if D in(Bg,'cbnz'):
		H=J(B[0]);E=j if B[0][0]=='w'else M;f=D==Bg
		def L():
			if(R[H]&E==0)==f:return A.pcof(B[1],V)
		return L
	if D in(Bh,B8):
		H=J(B[0]);Q=K(B[1][1:],0);f=D==Bh
		def L():
			if(not R[H]>>Q&1)==f:return A.pcof(B[2],V)
		return L
	if D in('br',Bi,d):
		H=J(B[0])if B else 30
		def L():
			B=R[H]
			if D==Bi:R[30]=CODEB+4*(V+1)
			return A.ja(B)
		return L
	if D=='adr':
		I=J(B[0],1)
		def L():C=A.dadr(B[1]);R[I]=CODEB+4*A.pcof(B[1],V)if C is G else C
		return L
	if D=='nop':return lambda:G
	X=B[0][0]=='w'if B else 0;E=j if X else M;O=l if X else Ap;Q=32 if X else 64
	if D in('add','sub','adds','subs',Bj,'cmn'):
		if D in(Bj,'cmn'):B=[Bf if X else Be]+B;D='subs'if D==Bj else'adds'
		I=J(B[0],1);H=J(B[1]);P=Aq(B,2,E);k=D[:3]=='sub'
		if D[-1]!='s':
			if k:
				def L():R[I]=R[H]-P()&E
			else:
				def L():R[I]=R[H]+P()&E
			return L
		def L():
			A=R[H]&E;B=P()
			if k:C=A-B&E;F[2]=K(A>=B);F[3]=K(O(A)-O(B)!=O(C))
			else:C=A+B&E;F[2]=K(A+B>E);F[3]=K(O(A)+O(B)!=O(C))
			F[0]=C>>Q-1;F[1]=K(C==0);R[I]=C
		return L
	if D in('and','orr','eor',Bk,'tst','bic'):
		if D=='tst':B=[Bf if X else Be]+B;D=Bk
		I=J(B[0],1);H=J(B[1]);P=Aq(B,2,E);q=D==Bk
		def L():
			A=R[H]&E;B=P();C=(A&B if D[:3]=='and'else A|B if D=='orr'else A^B if D=='eor'else A&~B)&E
			if q:F[0]=C>>Q-1;F[1]=K(C==0);F[2]=F[3]=0
			R[I]=C
		return L
	if D=='mov':
		I=J(B[0],1)
		if B[1][0]==Y:b=K(B[1][1:],0)&E;return lambda:R.__setitem__(I,b)
		H=J(B[1]);return lambda:R.__setitem__(I,R[H]&E)
	if D in(AP,Af,Ag,'ror')and B[2][0]==Y:
		I=J(B[0],1);H=J(B[1]);U=K(B[2][1:],0)
		def L():A=R[H]&E;R[I]=A<<U&E if D==AP else A>>U if D==Af else O(A)>>U&E if D==Ag else(A>>U|A<<Q-U)&E
		return L
	if D in('ubfx','sbfx','bfi','bfxil','ubfiz',Ah,B7):
		I=J(B[0],1);H=J(B[1])
		if D in(Ah,B7):return lambda:R.__setitem__(I,(l(R[H])if D==Ah else R[H]&j)&M)
		a=K(B[2][1:]);g=K(B[3][1:]);Z=(1<<g)-1
		def L():
			A=R[H]&E
			if D=='ubfx':R[I]=A>>a&Z
			elif D=='sbfx':B=A>>a&Z;R[I]=(B-(1<<g)if B>>g-1 else B)&E
			elif D=='bfi':R[I]=(R[I]&E&~(Z<<a)|(A&Z)<<a)&E
			elif D=='bfxil':R[I]=(R[I]&E&~Z|A>>a&Z)&E
			else:R[I]=(A&Z)<<a&E
		return L
	if D in('mul','madd','msub','mneg','umulh','smulh','umull','smull'):
		I=J(B[0],1);H=J(B[1]);N=J(B[2]);U=J(B[3])if C(B)>3 else 31
		def L():
			if D=='umulh':R[I]=(R[H]&M)*(R[N]&M)>>64
			elif D=='smulh':R[I]=Ap(R[H])*Ap(R[N])>>64&M
			elif D=='umull':R[I]=(R[H]&j)*(R[N]&j)&M
			elif D=='smull':R[I]=l(R[H])*l(R[N])&M
			elif D=='madd':R[I]=R[U]+(R[H]&E)*(R[N]&E)&E
			elif D=='msub':R[I]=R[U]-(R[H]&E)*(R[N]&E)&E
			elif D=='mul':R[I]=(R[H]&E)*(R[N]&E)&E
			else:R[I]=-(R[H]&E)*(R[N]&E)&E
		return L
	if D in('csel',B9,Bl,'cneg','cset','cinc'):
		I=J(B[0],1)
		if D=='cset':H=N=31;T=AJ[BR[B[1]]];D=B9
		elif D=='cinc':H=N=J(B[1]);T=AJ[BR[B[2]]];D=B9
		elif D=='cneg':H=N=J(B[1]);T=AJ[BR[B[2]]];D=Bl
		else:H=J(B[1]);N=J(B[2]);T=AJ[B[3]]
		def L():
			if T():R[I]=R[H]&E
			else:R[I]=(R[N]+1 if D==B9 else-R[N]if D==Bl else R[N])&E
		return L
	if D in('sdiv','udiv'):
		I=J(B[0],1);H=J(B[1]);N=J(B[2])
		def L():
			C=R[H]&E;A=R[N]&E
			if D=='udiv':R[I]=C//A if A else 0
			else:F=O(C);B=O(A);R[I]=w(F)//w(B)*(1 if(F<0)==(B<0)else-1)&E if B else 0
		return L
	if D in(Bm,'mvn','orn','eon','negs'):
		I=J(B[0],1)
		if D in(Bm,'mvn'):P=Aq(B,1,E);return(lambda:R.__setitem__(I,-P()&E))if D==Bm else lambda:R.__setitem__(I,~P()&E)
		H=J(B[1]);P=Aq(B,2,E);return lambda:R.__setitem__(I,(R[H]&E|~P()&E if D=='orn'else R[H]&E^~P()&E)&E)
	if D in(AP,Af,Ag,'ror'):
		I=J(B[0],1);H=J(B[1]);N=J(B[2])
		def L():A=R[H]&E;B=R[N]&Q-1;R[I]=A<<B&E if D==AP else A>>B if D==Af else O(A)>>B&E if D==Ag else(A>>B|A<<Q-B)&E
		return L
	if D in('movz','movk'):
		I=J(B[0],1);b=K(B[1][1:],0);h=K(B[2].split(Y)[1])if C(B)>2 else 0
		if D=='movz':return lambda:R.__setitem__(I,b<<h&E)
		return lambda:R.__setitem__(I,R[I]&E&~(65535<<h)|b<<h)
	if D in('ccmp','ccmn'):
		H=J(B[0]);T=AJ[B[3]];c=K(B[2][1:],0);r=J(B[1])if B[1][0]!=Y else 0;s=K(B[1][1:],0)&E if B[1][0]==Y else 0;v=B[1][0]==Y
		def L():
			if T():
				A=R[H]&E;B=s if v else R[r]&E
				if D=='ccmp':C=A-B&E;F[2]=K(A>=B);F[3]=K(O(A)-O(B)!=O(C))
				else:C=A+B&E;F[2]=K(A+B>E);F[3]=K(O(A)+O(B)!=O(C))
				F[0]=C>>Q-1;F[1]=K(C==0)
			else:F[0]=c>>3&1;F[1]=c>>2&1;F[2]=c>>1&1;F[3]=c&1
		return L
	if D in('rbit','clz'):
		I=J(B[0],1);H=J(B[1])
		if D=='clz':return lambda:R.__setitem__(I,Q-(R[H]&E).bit_length())
		return lambda:R.__setitem__(I,K(format(R[H]&E,'0%db'%Q)[::-1],2))
	if D in('ldr','ldrh','ldrb','ldrsw','ldrsh','ldrsb','str','strh','strb'):
		H=4 if D=='ldrsw'else 2 if D[3:]in('h','sh')else 1 if D[3:]in('b','sb')else 4 if X else 8;I=J(B[0],D[0]=='l');S=AY(B,1)
		if D[0]=='l':
			m=D[3:5]=='sw'or D[3:5]=='sh'or D[3:5]=='sb'
			def L():
				A=MEM.rd(S(),H)
				if m and A>>8*H-1:A-=1<<8*H
				R[I]=A&E if not m else A&M
			return L
		def L():MEM.wr(S(),H,R[I])
		return L
	if D=='ldpsw':
		I=J(B[0],1);e=J(B[1],1);S=AY(B,2)
		def L():A=S();B=l(MEM.rd(A,4))&M;C=l(MEM.rd(A+4,4))&M;R[I]=B;R[e]=C
		return L
	if D in(BA,'stp'):
		H=4 if X else 8;I=J(B[0],D==BA);e=J(B[1],D==BA);S=AY(B,2)
		if D==BA:
			def L():A=S();B=MEM.rd(A,H);C=MEM.rd(A+H,H);R[I]=B;R[e]=C
			return L
		def L():A=S();MEM.wr(A,H,R[I]);MEM.wr(A+H,H,R[e])
		return L
	if D in('ldadd',CR,'stadd'):
		if D=='ldadd':n=J(B[0]);o=J(B[1],1);S=AY(B,2);i=8 if B[0][0]==u else 4
		else:n=J(B[0]);o=33;S=AY(B,1);i=1 if D==CR else 8 if B[0][0]==u else 4
		def L():A=S();B=MEM.rd(A,i);MEM.wr(A,i,B+R[n]);R[o]=B
		return L
	raise Fault('unsupported instruction: '+t)
D9='b','bl','br',Bi,d,'adr',Bg,'cbnz',Bh,B8
DA='add sub adds subs cmp cmn ccmp ccmn rbit clz sdiv udiv neg mvn orn eon movz movk and orr eor ands tst bic mov lsl lsr asr ubfx sbfx bfi bfxil ubfiz sxtw uxtw mul madd msub mneg umulh smulh umull smull csel csinc csneg cneg cset cinc ldr ldrh ldrb ldrsw ldrsh ldrsb str strh strb ldp ldpsw stp ldadd stadd staddb b bl br blr ret b.cond cbz cbnz tbz tbnz adr nop'.split()
CH={}
class A4:
	def __init__(A,text,ext=(),sink=0,over=G,lazy=('TOWER','SA'),dl=('CHAFL',)):
		A.mac={};A.sym={};A.code=[];A.txt=[];A.sg=[];A.hits=[];A.lab={};A.alt=0;A.cs=G;A.ext=X(ext);A.nl={};A.dc={};A.cls={};A.sink=sink;A.over=over or{};A.lazy=lazy;A.dl=dl;A.lz={};A.lw=[];A.by=[];A.lp={};A.bp={};A.sc=0;A.nm={};A.nn={};A.labs=[];A.ct={};A.dq=[];A.used=X();A.mm={};A.fx=0;A.da={};A.pd=G;B=[]
		for D in text.split(BB):
			C=D.strip()
			if C[:6]==BC:B.append('@seg '+C.split()[1]);continue
			if C[:6]==BD:B.append('@end');continue
			E=D.split(Ai,1)[0].strip()
			if E:B.append(E)
		F=A.parse(B);A.sig=hash((repr(a(A.mac.items())),AO(A for A in B if A[0]=='.')));A.ex(F,{})
	def parse(K,L):
		F=[];B=0
		while B<C(L):
			M=L[B];B+=1;E=M.split(G,1)[0];H=M[C(E):].strip()
			if E in(Bn,Bo,'.irp'):
				I=B;A=0
				while 1:
					D=L[I].split(G,1)[0]
					if D in(Bn,Bo,'.irp','.irpc'):A+=1
					elif D in('.endm','.endr'):
						if A==0:break
						A-=1
					I+=1
				N=K.parse(L[B:I]);B=I+1
				if E==Bn:P=split(H.replace(' ',',',1));K.mac[P[0]]=[(aa.split('=')[0].split(':')[0],aa.split('=')[1]if'='in aa else G,aa.endswith(':vararg'))for aa in P[1:]],N
				elif E==Bo:F.append(('r',H,N))
				else:Q=split(H);F.append(('i',Q[0],Q[1:],N))
			elif E in(Aj,Bp,CS):
				J=[(E,H,[])];A=0
				while 1:
					O=L[B];B+=1;D=O.split(G,1)[0]
					if D in(Aj,Bp,CS):A+=1
					elif D=='.endif':
						if A==0:break
						A-=1
					elif A==0 and D=='.elseif':J.append((Aj,O[7:].strip(),[]));continue
					elif A==0 and D=='.else':J.append((Aj,'1',[]));continue
					J[-1][2].append(O)
				F.append(('c',[(ab,ac,K.parse(ad))for(ab,ac,ad)in J]))
			else:F.append(('l',M))
		return F
	def sub(A,l,env):
		if'\\'not in l:return l
		return re.sub('\\\\(\\w+)(\\\\\\(\\))?',lambda m:T(env.get(m.group(1),m.group(0))),l)
	def pure(A,d):
		D=A.nm.get(d)
		if D is G:A.nm[d]=0;B=A.mac[d][1];D=A.nm[d]=2 if C(B)==1 and B[0][0]=='l'and re.fullmatch('[\\w.$\\\\()]+:',B[0][1])else K(any(A.need(B)for B in B))
		return D
	def need(A,n):
		B=A.nn.get(id(n))
		if B is G:
			C=n[0]
			if C=='l':D=n[1].split(G,1)[0];B=D in A.mac and A.pure(D)>0
			elif C=='c':B=any(A.need(D)for(B,B,C)in n[1]for D in C)
			else:B=any(A.need(B)for B in n[2 if C=='r'else 3])
			A.nn[id(n)]=B
		return B
	def ex(A,N,e):
		for C in N:
			F=C[0]
			if F=='l':A.line(C[1],e)
			elif A.sc and not A.need(C):continue
			elif F=='c':
				for(G,D,I)in C[1]:
					D=A.sub(D,e)
					if O(D,A.sym)if G==Aj else(D in A.sym or D in A.lab)==(G==Bp):A.ex(I,e);break
			elif F=='r':
				for K in B(O(A.sub(C[1],e),A.sym)):A.ex(C[2],e)
			else:
				for J in C[2]:H=E(e);H[C[1]]=J;A.ex(C[3],H)
	def line(A,l,e):
		if'\\'in l:l=A.sub(l,e)
		H=A.cls.get(l)
		if H:
			if not A.sc:A.addn(H)
			return
		if l[0]=='@':
			if not A.sc:A.cs=l.split()[1]if l[1]=='s'else G
			return
		B=l.split(G,1)[0];D=l[C(B):].strip()
		if B=='.error':raise CQ('assembler .error: '+D)
		if B in('.set','.equ'):E,J=split(D);A.sym[E]=A.over[E]if E in A.over else O(J,A.sym);A.fx+=1;return
		if B=='.altmacro':A.alt=1;return
		if B=='.noaltmacro':A.alt=0;return
		if B in('.long','.byte')and A.pd and not A.sc:A.da[A.pd]=C(A.lw),C(A.by),C(A.dq);A.pd=G
		if B=='.long':
			if not A.sc:I=re.fullmatch('(\\S+)\\s*-\\s*(\\S+)',D);A.lw.append(I.groups()if I else O(D,A.sym));A.fx+=1
			return
		if B=='.byte':
			if not A.sc:A.by+=[O(ae,A.sym)for ae in split(D)];A.fx+=1
			return
		if B[0]=='.'and B[-1]!=':'and not re.match('\\.[A-Za-z_]\\w*:',l):return
		F=re.match('([A-Za-z_.$][\\w.$]*|\\d+):\\s*(.*)$',l)
		if F:
			if not A.sc:A.defl(F.group(1))
			l=F.group(2)
			if not l:return
			B=l.split(G,1)[0];D=l[C(B):].strip()
		if B in A.mac:A.call(B,D)
		elif not A.sc:A.add(l)
	def call(A,d,a):
		Q,H=A.mac[d];J=split(a);B={}
		for(D,(K,R,S))in h(Q):
			if S:B[K]=','.join(J[D:]);break
			F=J[D]if D<C(J)else R
			if F and F[0]=='%'and A.alt:F=T(O(F[1:],A.sym))
			B[K]=F
		if A.sc:
			L=A.pure(d)
			if L==2:A.labs.append(A.sub(H[0][1],B)[:-1])
			elif L:A.ex(H,B)
			return
		if d in A.dl:
			if A.pd:A.da[A.pd]=C(A.lw),C(A.by),C(A.dq);A.pd=G
			A.dq.append((d,AO(B.values()),A.alt));return
		if d in A.lazy:
			M=E(A.sym);D=A.sig,d,AO(B.items());I=CH.get(D)
			if I is G:
				A.sc=1;A.labs=[]
				try:A.ex(H,B)
				finally:A.sc=0;A.sym=E(M)
				I=CH[D]=A.labs
			U=[d,B,M,I,A.alt]
			for F in I:A.lz[F]=U
			return
		D=d,AO(B.values())
		if A.sc:A.ex(H,B);return
		N=A.mm.get(D)
		if N is not G:
			for P in N:A.addn((P,P.split(G,1)[0],y))
			return
		V=A.fx;W=C(A.code);A.ex(H,B)
		if A.fx==V:A.mm[D]=A.txt[W:]
	def defl(A,n):
		A.fx+=1;A.pd=n
		if n.isdigit():A.nl.setdefault(K(n),[]).append(C(A.code))
		else:A.lab[n]=C(A.code);A.lp[n]=C(A.lw);A.bp[n]=C(A.by)
	def add(B,l):A=re.sub('\\s*,\\s*',', ',re.sub('\\s+',' ',l));C=B.cls[l]=A,A.split()[0],any(not re.fullmatch('-?(0x[0-9a-fA-F]+|\\d+)',A.strip())for A in re.findall('#([^,\\]!]*)',A));B.addn(C)
	def addn(A,c):
		B,D,E=c
		if E:B=re.sub('#([^,\\]!\\s][^,\\]!]*)',lambda m:Y+T(O(m.group(1).strip(),A.sym)),B)
		if B not in A.dc or D in D9 or D[:2]=='b.':A.dc[B]=D8(B,C(A.code),A)
		A.code.append(A.dc[B]);A.txt.append(B);A.sg.append(A.cs);A.hits.append(0);A.pd=G
	def unit(A,n):
		B,D,F,H,I=A.lz[n]
		for J in H:A.lz.pop(J,G)
		K,L,M=A.cs,A.sym,A.alt;A.cs='chain';A.alt=I;A.sym=E(F);A.sym['_gid']=C(A.code)
		try:A.ex(A.mac[B][1],D)
		finally:A.cs,A.sym,A.alt=K,L,M
	def pcof(A,n,at=G):
		if n in A.lab:return A.lab[n]
		B=re.fullmatch('(\\d+)([bf])',n)
		if B:D=A.nl[K(B.group(1))];return P(A for A in D if A<=at)if B.group(2)=='b'else S(A for A in D if A>at)
		if n in A.lz:A.unit(n);return A.lab[n]
		if n in A.ext:A.lab[n]=C(A.code);A.code.append(lambda:(A for A in()).throw(AX(n)));A.txt.append('@halt '+n);A.sg.append(G);A.hits.append(0);return A.lab[n]
		raise Fault('undefined label '+n)
	def ja(B,a):
		if a==STOPA:raise AX(d)
		A=a-CODEB
		if A<0 or A&3 or A>>2>=C(B.code):raise Fault('branch to non-code %x'%a)
		return A>>2
	def dadr(B,l):A=B.da.get(l);return G if A is G else A3+4*(A[0]+7*B.sym[z]*A[2])+A[1]
	def cte(A,f,ell,fam,n):
		B=f,ell,fam,n
		if B not in A.ct:
			F=Az(A for A in A.mac[CT][1]if A[0]=='r');H=A.lw;A.lw=[];I=E(A.sym);J=A.alt;A.alt=A.dq[0][2];A.sym['_n']=n
			try:A.ex(F[2],E(f=T(f),ell=T(ell),fam=T(fam)))
			finally:D=A.lw;A.lw=H;A.sym=I;A.alt=J
			assert C(D)==1;A.ct[B]=D[0][0]if D[0]else G
		return A.ct[B]
	def run(D,start,stops=(),mx=4*10**6,hooks=G):
		I=stops;E=D.pcof(start);C=D.code;J=D.hits;H=[]
		for F in I:A=D.pcof(F);H.append((A,C[A]));C[A]=lambda n=F:(A for A in()).throw(AX(n))
		for(F,N)in(hooks or{}).items():A=D.pcof(F);K=C[A];H.append((A,K));C[A]=lambda fn=N,o=K:(fn(),o())[1]
		try:
			try:
				for P in B(mx):J[E]+=1;L=C[E]();E=E+1 if L is G else L
			except AX as M:
				if T(M)in I:J[E]-=1
				return T(M)
			raise Fault('step limit')
		finally:
			for(A,O)in H:C[A]=O
	def counts(C):
		A={}
		for(B,D)in zip(C.sg,C.hits):
			if B:A[B]=A.get(B,0)+D
		return A
	def clear(A):
		A.use()
		for D in B(C(A.hits)):A.hits[D]=0
	def use(A):A.used|={af.split()[0]for(af,ag)in zip(A.txt,A.hits)if ag and af[0]!='@'}
SPB=1048576
AZ=SPB+512
PTB=10485760
EBB=11534336
Aa=28672<<32
ORBA=65536
GCA=69632
BS=536870912
DS=131072
CI=73728
Ar=73984
AK=74240
DIRB,A5,AL,A6,g,A7,k,As=1<<44,2<<44,3<<44,4<<44,5<<44,30720<<32,31744<<32,32256<<32
BT=1<<48
A8=26624<<32,27136<<32
At=1
Ab={}
def pt_kv(A):B=A.sym;C=E(K=B[Bq],OFF=B[CU],NMAXPOS=B[z],NULLSLOT=B[Br],LASTADD=B[AQ],SLOTN=B['SLOTN'],SPANK=B[AR],REGBYTES=B[BE]);C.update(Ab);return C
def pt_pieces(rows,kv):
	B=[]
	for(E,A,C)in rows:
		while A<=C:D=S(C,A+kv[z]-1);B.append((E,A,D));A=D+1
	return B
def CJ(O,rows,kv):
	A=kv;F,G=A['K'],A[BF];H=pt_pieces(rows,A);I=[]
	for(J,(D,B,K))in h(H):
		E=B+G
		if not 0<=E<=F-3:E=G+B%3
		I+=[D*O[p][1]+B*O[p][0]&M,D*O[q][1]+B*O[q][0]&M,E+F*(K-B+(A[AQ]if J==C(H)-1 else 0))+A[AR]*(D%3)]
	return I
def Au(A,fl,b,ph,sl,la,mut=0):
	C=fl,b,ph,sl,la,mut;D=A.__dict__.setdefault('pc',{})
	if C in D:return D[C]
	F=A.sym[z];E=(b>>ph|b<<3-ph)&7;J=A.lw[A.lp[Bs]+la][0]
	if mut==1:E=0
	G=A.by[A.bp[CV]+E];I=H((~E&7)>>A%3&1 for A in B(sl+1))if sl<F else 0;D[C]=J if sl>=F or G==255 or not I else A.cte(fl,la,G,I);return D[C]
def pt_mem(O,A,fl=8,tmut=0):
	F=fl;D=pt_kv(A);J=A.sym;K=W._add_st_f3;N=O[r];H=[AL+8];E=O['c5'];E.clear();I,P,R,L=J[A0],J[BG],J[AS],J[AE];Q,U,S,V=D[BE],D[AR],D['K'],D[AQ];MEM.reg=[];MEM.pg={};MEM.add(SPB-4096,SPB+8192);MEM.add(AL,AL+(1<<36));MEM.add(PTB,PTB+4096);MEM.add(EBB,EBB+(1<<24));MEM.add(As,As+64);MEM.wr(As+16,8,S*(D[Br]+D[AQ]))
	def X(ws):
		A=H[0]
		for(B,D)in h(ws):MEM.wr(A+8*B,8,D)
		H[0]+=8*C(ws);return A
	def T(s,n):
		A=0
		for C in B(n):A=K(A,s)
		return A
	def Y(u,c):A=[O[B3][K(u,T(N[1],(c+aa)%3))]for aa in B(3)];return A[0]*64+A[1]*8+A[2]
	def Z(a):
		A=a-DIRB>>3;B='d',A
		if B not in E:F=A&(1<<I)-1;G=A>>I&(1<<P)-1;H=A>>I+P;C=pt_rows(O,H,F,G);E[B]=X(CJ(O,C,D))if C else As
		return E[B]
	def a(a):
		D=a-A5>>3;F='e',D
		if F not in E:
			J,K=D>>L+1,D&(1<<L+1)-1;B=pt_e3(O,J,K)if J<1<<R and K<=1<<L else G;A=H[0]
			if B is G:MEM.wr(A,8,0);MEM.wr(A+8,4,1<<31);H[0]+=16
			else:
				MEM.wr(A,8,B[0]&M)
				for(N,P)in h(B[1]):MEM.wr(A+8+4*N,4,P<<I)
				MEM.wr(A+8+4*C(B[1]),4,1<<31);H[0]+=8+4*C(B[1])+4+7&~7
			E[F]=A
		return E[F]
	def c(a):
		B,A=b(a-g,80);A>>=3
		if A==0:return g+80*K(B,N[2])
		C,D=b(A-1,3);return A8[F-8]+Y(B,C)*Q+8*D
	def d(a):
		C,D=b(a-A8[F-8],Q);B,E=b(D>>3,U)
		if B>2:raise Fault('tb index')
		G,H=b(E,S);I,J=b(G,V);return CODEB+4*A.pcof(Au(A,F,C>>3*(2-B)&7,H%3,J,I,tmut))
	def e(a):A,B=b(a-A7>>3,729);return g+80*K(B,T(N[2],A))
	def f(a):A=a-A6>>1;B=(A+2048)//4096;C=A-4096*B;return 8+8*(3*(C%3)+B%3)
	MEM.add(DIRB,DIRB+(8<<2*I+P),Z);MEM.add(A5,A5+(8<<R+L+1),a);MEM.add(A8[F-8],A8[F-8]+512*Q,d);MEM.add(g,g+58320,c);MEM.add(A7,A7+17496,e);MEM.add(k-(1<<20),k+(1<<20),lambda a:(a-k>>1)%3*729);MEM.add(A6-(1<<28),A6+(1<<28),f)
def BU(O,S,p0=EBB):C,A,B=O[c],O[p],O[q];return[(S['FR_MU'],C[2][1]),(S['FR_MU']+8,C[2][0]),(S[CW],A[1]),(S[CW]+8,B[1]),(S[CX],A[2]),(S[CX]+8,B[2]),(S['FR_E3'],A5),(S['FR_M3X'],k),(S['FR_P2'],A7),(S['FR_MB'],A6),(S['FR_PTR0'],p0),(S[CY],3*A[0]),(S[CY]+8,3*B[0])]
def DB(O,x):
	A=0
	for C in B(3,12):
		for D in B(x[C]%3):A=W._add_st_f3(A,O[r][C])
	return A
def pt_pk(O,base,x2,x1,x0):C=O['B'];A=[base[aa]+x2*C[2][aa]+x1*C[1][aa]+x0*C[0][aa]for aa in B(12)];return G if W._bad_f3(A)else(BT+16*W._mr_pack_a(A[:6])&M,W._mr_pack_b(A[6:]))
def DC(O,nd,kv,fl,p0=EBB):
	V=nd[u];W=[H(V[A]*O['B'][A][aa]for A in B(3,12))for aa in B(12)];F={};K=L=G=M=N=0;Q=pt_cands(O,nd);A=p0;R=[]
	for(X,d,Y,Z,a)in Q:
		A=p0;S=pt_pieces(a,kv);K+=P(1,C(S))
		for(b,T,c)in S:
			D=0;E=-1
			for U in B(T,c+1):
				I=pt_pk(O,W,X,Z+b,Y+U)
				if I:
					F[I]=F.get(I,0)+1;D+=1;J=A;A+=8
					if E<0:E=U-T
					if fl==9 and J>>14&1:G+=1;R.append(J-16384+8);A=J-16384+24
			if D:L+=(2 if E>0 else 0)+6*D+2*(D-1)+(D if fl==9 else 0);M+=1;N+=E>0
	return F,C(Q),K,L+8*G,G,A,R,M,N
def DD(O,A,nd,fl=8,p0=EBB):
	I=nd;A0=pt_kv(A);S=I[u];A2,a,b=O[c],O[p],O[q];P=A.sym
	for A3 in[aa for aa in MEM.pg if aa<<12<AL]:del MEM.pg[A3]
	for A4 in B(C(R)):R[A4]=0
	EM.clear();A.clear();J=lambda a,v:MEM.wr(a,8,v&M);J(AZ+P['L3_AB'],O['A0']+H(S[A]*a[A]for A in B(4,12)));J(AZ+P['L3_AB']+8,O['B0']+H(S[A]*b[A]for A in B(4,12)));J(AZ+P['L3_D3'],a[3]);J(AZ+P['L3_D3']+8,b[3])
	for(A5,A6)in BU(O,P,p0):J(SPB+A5,A6)
	J(SPB+480,19088743);MEM.wr(SPB+P[Bt],2,DB(O,S));J(PTB,p0);R[17]=PTB;R[26]=a[0];R[27]=b[0];R[28]=DIRB;R[29]=A2[1][0]&M;R[0]=S[3]&M;R[2]=I[A_]&M;R[6]=I['rem3'];R[7]=I['c1b']&M;R[4]=I['c0b']&M;R[30]=AZ;R[32]=SPB;V,D,K,U,j,A7,A8,k,l=DC(O,I,A0,fl,p0)
	try:Y=A.run('pt_handoff3',mx=3*(34+27*D+6*K+U)+2000);m=''
	except Fault as A9:Y=BH;m=T(A9)
	L={}
	for n in EM:L[n]=L.get(n,0)+1
	Z=A.hits;d=A.txt;e=H(B for(A,B)in zip(d,Z)if A[0]!='@');o=Z[A.lab[CZ]]-1;r=H(B for(A,B)in zip(d,Z)if A=='br x9');F=DE(A);f=A.counts();s=f.get('chain',0)+f.get(Ca,0)+f.get(AT,0);g=W.HT_NEW;h=t=0
	if Y!=BH and I.get('ball')is not G:AA=[H(S[A]*O['B'][A][ab]for A in B(3,12))for ab in B(12)];v={AB for ac in I['ball']if(AB:=pt_pk(O,AA,ac[0],ac[1],ac[2]))};h=C(v-X(L));t=C(X(L)-v)
	w=C(EM);AC=s-F[AU]*(K-D)-F[Bu]*D;x=H(B for(A,B)in zip(d,Z)if A=='sub x0, x9, #4, lsl #12');z=AC==8*w-2*k+2*l+(w if fl==9 else 0)+8*x
	try:i=MEM.rd(PTB,8)==A7 and N(MEM.rd(A,8)==2046 and MEM.rd(A+8,8)==19088743 for A in A8)
	except Fault:i=y
	AD=Y==A1 and L==V and h==0 and i and z and x==j and o==D and r==K and e==F['pre']+(F[Bv]+F[AU])*D+F[AU]*(K-D)+F[Bu]*D+U and s==U+F[AU]*(K-D)+F[Bu]*D and e==g[AF]-At+g[AV]*D+g[AW]*K+U;return E(ok=Q(AD),end=Y,fault=m,leaves=C(EM),want=H(V.values()),cands=D,rows=K,E=U,flushes=j,rows_emitting=k,rows_skip=l,ident=Q(z),dyn=e,c_run=o,R_run=r,ball_missing=h,state_ok=Q(i),extras=t,mismatch=0 if L==V else C(X(L.items())^X(V.items())))
def DE(A):
	if hasattr(A,'sh'):return A.sh
	B=lambda g:H(1 for A in A.sg if A==g);D=A.lab[CZ];C=D
	while not A.txt[C].startswith(B8):C+=1
	A.sh=E(hand=B(Cb)+B(Cc)+B(Cd),sent=C-D+1,exit=B(Bw),node=B(Ak),row=B(Al),last=B(AT));A.sh['pre']=A.sh[Bx]+A.sh['sent']+A.sh['exit'];return A.sh
def Z(O,A,nodes,fl=8,tmut=0,p0=EBB,keep=0,stop=0):
	pt_mem(O,A,fl,tmut);B=E(nodes=0,ok=0,leaves=0,cands=0,rows=0,E=0,dyn=0,faults=0,bad=0,flushes=0,ball_missing=0,rows_emitting=0,rows_skip=0,flavour=fl)
	for G in nodes:
		C=DD(O,A,G,fl,p0);B[v]+=1;B[D]+=C[D];B[e]+=not C[D];B['faults']+=C['end']==BH
		for F in(BI,s,x,U,By,BJ,'ball_missing',Bz,B_):B[F]+=C[F]
		if not C[D]and Ce not in B:B[Ce]=C
		if keep:B.setdefault('runs',[]).append(C)
		if stop and not C[D]:break
	A.use();B[D]=B[e]==0 and B[v]>0;return B
CK,CL,CM,m,STKB=1879048192,1895825408,1912602624,1929379840,1946157056
def DF(O,u,c):
	D=W._add_st_f3;E=O[r];A=[]
	for F in B(3):
		C=0
		for G in B((c+F)%3):C=D(C,E[1])
		A.append(O[B3][D(u,C)])
	return A[0]*64+A[1]*8+A[2]
def pt_bmem(O,A,fl=8,cell=(0,1),bins=(0,1),add=0):
	B=A.sym;F,G,H,I=B[A0],B[BG],B[AS],B[AE]
	if not add:MEM.reg=[];MEM.pg={}
	for(E,C)in((ORBA,16384),(STKB,8192),(BS,B['S_SIZE']+8192),(AL,1<<36),(DIRB,8<<2*F+G),(A5,8<<H+I+1),(g,58320),(A7,17496),(k-(1<<20),2<<20),(CK,65536),(CL,1<<19),(CM,4096),(m,4194304)):MEM.add(E,E+C)
	MEM.add(Aa,Aa+(8<<40),lambda a:BO((a-Aa>>3)+1));J=lambda n,v:MEM.wr(ORBA+B[n],8,v&M)
	for(C,D)in((C0,O[AA][1]),('O_EK1',O[B0][1]),('O_GW0',O[B1][0]),('O_MU10',O[c][1][0]),('O_K0',O[AA][0]),('O_DA1',O[p][1]),('O_DB1',O[q][1]),('O_LUT',Aa),('O_DA0',O[p][0]),('O_DB0',O[q][0]),('O_MU21',O[c][2][1]),('O_MU20',O[c][2][0]),('O_DA2',O[p][2]),('O_DB2',O[q][2]),('O_K2',O[AA][2]),('O_EK2',O[B0][2]),('O_AE0',O[r][0]),('O_AE1',O[r][1]),('O_AE2',O[r][2]),('O_FLAV',fl-8),('O_A0',O['A0']),('O_B0',O['B0']),('O_PTRS',PTB),('O_PTR0',EBB),('O_DIR',DIRB),('O_E3',A5),('O_ARENA',AL),('O_REC',g),('O_P2',A7),('O_EMPTY',AK),('O_EMPTY3',Ar),('O_SCR',BS),('O_GLB',GCA),('O_SENT',CI),('O_CELL0',cell[0]),('O_CELL1',cell[1]),('O_BIN0',bins[0]),('O_BIN1',bins[1])):J(C,D)
	for(C,D)in(('G_M3X',k),('G_MB',A6),('G_TB8',A8[0]),('G_TB9',A8[1]),('G_LANE',CK),('G_SOFL',CL),('G_BAD',CM),(C1,m)):MEM.wr(GCA+B[C],8,D)
def A9(A,lab,x0=0,x1=0,stops=(),hooks=G):
	for D in B(C(R)):R[D]=0
	R[28]=ORBA;R[30]=STOPA;R[32]=STKB+4096;R[0]=x0;R[1]=x1;return A.run(lab,stops=stops,hooks=hooks)
def DG(kv,p):
	A=[]
	for C in B(1200):
		A+=[MEM.rd(p,8),MEM.rd(p+8,8),MEM.rd(p+16,8)];p+=24
		if A[-1]%kv[AR]//kv['K']>=kv[AQ]:return A
def pt_pair_cost(O,A,f1,bp):
	F=pt_kv(A);Q=A.sym;G=1<<Q[A0];e=Q[Ad];B=pt_geo(O,f1,bp);R=H(1 for(A,B,B)in B if A<=0);f=C(B)-R;T=C(B);I=E(pair_init=38,gen_neg=34*R+7,gn_to_gp=22,gen_pos=34*f+7)
	if not B:I.update(pair_post=14,pair_empty=1+3*(G//2)+2);return I
	g=not(S(B for(A,B,A)in B)>>32>=-F[BF]and P(B for(A,A,B)in B)>>32<=F['K']-4-F[BF]);U={0}
	for(l,D,i)in B:
		for j in(D,i):
			V=G-(j>>e&G-1)
			if V<=G-1:U.add(V)
	D=a(U);k=C({aa>>3 for aa in D if aa});W=J=X=Y=0
	for(Z,b)in h(D):
		K=(D[Z+1]if Z+1<C(D)else G)-b;L=pt_rows(O,f1,b,bp);c=pt_pieces(L,F);M=C(L);W+=8+(K&1)+2+(3*(K>>1)if K>>1 else 0)+2;J+=2
		if not g:J+=1+17*M+11*(T-M)
		else:
			N={}
			for(d,m,n)in c:N[d]=N.get(d,0)+1
			X+=H(11+17*A for A in N.values())+11*(T-M)+6*H(1 for(B,A,C)in c if not 0<=A+F[BF]<=F['K']-3)
		Y+=7 if L else 4
	I.update(pair_post=19,bp_scan=5*(G//8)+2*k+8*(C(D)-1),pair_arr_init=9,pair_arr=W+3,rows_fast=J,rows_slow=X,rows_post=Y);return I
def Av(O,A,keys,fl=8,stop=0):
	P=A.sym;I=pt_kv(A);L=P[A0];U=P[BG];pt_bmem(O,A,fl);C=E(cases=0,bad=0,arrays=0,instr=0,count_dev=0,slow_pairs=0,empty_pairs=0)
	if A9(A,Cf)!=d:C[e]+=1
	if[MEM.rd(AK,8),MEM.rd(AK+8,8),MEM.rd(AK+16,8)]!=[0,0,I['K']*(I[Br]+I[AQ])]or MEM.rd(Ar,8)!=0 or MEM.rd(Ar+8,4)!=1<<31 or MEM.rd(CI,4)!=1<<31:C[e]+=1;C['empty_init_bad']=1
	for(F,G)in keys:
		A.clear();V=A9(A,'pt_pair',F,G);J=A.counts();C[C2]+=1;K=V==d;C['instr']+=H(J.values())
		try:
			for Q in B(1<<L):
				R=pt_rows(O,G,Q,F);M=MEM.rd(DIRB+8*(G<<L+U|F<<L|Q),8);W=CJ(O,R,I);C['arrays']+=1;K&=M==AK if not R else M!=AK and DG(I,M)==W
				if stop and not K:break
			N=pt_pair_cost(O,A,G,F);S={aa:(J.get(aa,0),N.get(aa,0))for aa in X(J)|X(N)if J.get(aa,0)!=N.get(aa,0)}
			if S:C[Am]+=1;C.setdefault(Cg,(F,G,S))
			C['empty_pairs']+=not pt_geo(O,G,F);C[Ch]+=J.get('rows_slow',0)>0
		except Fault as Y:K=y;C[BH]=T(Y)
		C[e]+=not K
		if stop and(not K or C[Am]):break
	A.use();C[D]=C[e]==0 and C[Am]==0 and C[C2]>0;return C
def BV(O,A,cells,bins,fl=8):
	J=bins;V=A.sym;b=V[A0];W=V[AE];F=E(cells=0,bad=0,lists=0,empty=0,instr=0,count_dev=0);c=(1<<W)+1;J=a(X(J)|{0,c-1});i=pt_kv(A)
	for L in cells:
		pt_bmem(O,A,fl,cell=(L,L+1));A9(A,Cf);f=iter(J[1:]);j=[0]
		def g():
			A=Az(f,G)
			if A is not G:R[7]+=8*(A-1-R[9]);R[9]=A-1
		A.clear();h=A9(A,'pt_e3_build',hooks={'.Le3_nx':g});P=A.counts();F['cells']+=1;S=h==d;F['instr']+=H(P.values());Y=0;I=E(e3_build_cell=21,e3_cell=16,e3_cell_neg=1089,e3_cell_mid=2,e3_cell_pos=1080,e3_cell_end=5,e3_build_bin=0,e3_build_cand=0,e3_build_list_end=0,e3_build_empty=0,e3_build_next=3*C(J)+7+3)
		for Z in J:
			K=pt_e3(O,L,Z);N=C(K[1])if K else 0;T=MEM.rd(A5+8*(L<<W+1|Z),8);F['lists']+=1
			if not K:S&=T==Ar;I[Ci]+=15;I['e3_build_empty']+=1;F['empty']+=1
			else:S&=MEM.rd(T,8)==K[0]&M and[MEM.rd(T+8+4*ac,4)for ac in B(N+1)]==[ad<<b for ad in K[1]]+[1<<31];I[Ci]+=20+5*(N-Y);I['e3_build_cand']+=9*N;I['e3_build_list_end']+=6
			Y=N
		U={aa:(P.get(aa,0),ab)for(aa,ab)in I.items()if P.get(aa,0)!=ab};F[Am]+=Q(U)
		if U:F.setdefault(Cg,U)
		F[e]+=not S
	A.use();F[D]=F[e]==0 and F[Am]==0 and F['cells']>0;return F
def BW(O,A,fl=8,only=3):
	L=A.sym;G=W._add_st_f3;M=O[r];pt_bmem(O,A,fl);A.clear();P=E(bad=0)
	for X in B(C(R)):R[X]=0
	R[30]=STOPA;R[32]=STKB+4096;R[0]=GCA;D=A.run('pt_wheel_static')==d
	for Y in('pt_bms_build','pt_rec_build','pt_p2_build'):D&=A9(A,Y)==d
	D&=N(MEM.rd(BS+L['S_BMS']+A,1)==O[B3][A]for A in B(729));I=[]
	for F in B(729 if only&1 else 0):
		S=[MEM.rd(g+80*F+8*aa,8)for aa in B(10)];I.append(S[0]==g+80*G(F,M[2]))
		for T in B(3):I+=[S[1+3*T+ab]==A8[fl-8]+DF(O,F,T)*L[BE]+8*ab for ab in B(3)]
	D&=N(I)
	for U in B(3 if only&2 else 0):
		for F in B(729):
			J=0
			for Z in B(U):J=G(J,M[2])
			D&=MEM.rd(A7+8*(U*729+F),8)==g+80*G(F,J)
	V=A.counts();A.use();P.update(ok=Q(D),instr=H(V.values()),counts=E((A,V.get(A,0))for A in('wheel_static','ladd','bms_build','rec_build','p2_build')),bad=K(not D));return P
def BX(O,A):
	F=A.sym;pt_bmem(O,A);D=F[Cj];B=[]
	for(G,H,I)in((C0,D,0),(C0,D+1,1),('O_EK1',(1<<39)-1,3),('O_EK2',(1<<39)-1,5),('O_GW0',1<<63,4)):pt_bmem(O,A);MEM.wr(ORBA+F[G],8,H);A9(A,Ck);B.append(R[0]==I)
	pt_bmem(O,A);A9(A,Ck);B.append(R[0]==0);return E(ok=N(B),cases=C(B),CHK_K1=D)
def BY(A,rows=6,only=3):
	I=rows;F=only;J=A.sym;G=J['M3XR'];M=J['MBT0'];MEM.reg=[];MEM.pg={};O=A6-2*(4096*M+2048);D=1;A.clear()
	for(P,K)in((GCA-4096,16384),(STKB,8192),(k-2*G-256,4*G+512),(m,4194304),(O,I*8192+4096)):MEM.add(P,P+K)
	for(K,S)in(('G_M3X',k),('G_MB',A6),(C1,m)):MEM.wr(GCA+J[K],8,S)
	for L in B(C(R)):R[L]=0
	R[30]=STOPA;R[32]=STKB+4096;R[0]=GCA
	if F&1:D&=A.run('pt_m3x_build')==d and N(MEM.rd(k+2*A,2)==A%3*729 for A in B(-G,G))
	if F&2:
		for L in B(C(R)):R[L]=0
		R[30]=STOPA;R[32]=STKB+4096;R[0]=GCA;D&=A.run('pt_mb_build',stops=(C3,))==C3;D&=N(MEM.rd(m+2*(4096*A+C),2)==8+8*(3*((C-2048)%3)+A)for A in B(3)for C in B(4096));R[10]=I;D&=A.run(C3)==d;D&=N(MEM.rd(O+2*(4096*A+C),2)==8+8*(3*((C-2048)%3)+(A-M)%3)for A in B(I)for C in B(4096))
	T=A.counts();A.use();return E(ok=Q(D),instr=H(T.values()),m3x_entries=2*G if F&1 else 0,mb_blocks_entries=12288 if F&2 else 0,mb_copy_rows=I if F&2 else 0)
def DH(O,A,npairs=10,ncells=2,nbins=24,seed=1):
	U=npairs;J=A.sym;V,P,W,Q=J[A0],J[BG],J[AS],J[AE];K=random.Random(seed);F=[]
	for a in B(U):b=K.randrange(1<<V);c=K.randrange(1<<P);F.append(H(pt_pair_cost(O,A,b,c).values())+8)
	L=H(F)/C(F);d=H((A-L)**2 for A in F)/(C(F)-1)/C(F);D=[]
	for a in B(ncells):
		X=K.randrange(1<<W);e=[K.randrange((1<<Q)+1)for aa in B(nbins)];R=[]
		for f in e:G=pt_e3(O,X,f);Y=C(G[1])if G else 0;R.append(19 if not Y else 29+9*Y)
		G=pt_e3(O,X,1<<Q);D.append(2216+(1<<Q)*H(R)/C(R)+5*(C(G[1])if G else 0)+0)
	S=H(D)/C(D);g=H((A-S)**2 for A in D)/(C(D)-1)/C(D)if C(D)>1 else .0;M=1<<V+P;N=1<<W;T=120811+8*(1<<P)+7+M*L+N*S;Z=(M*M*d+N*N*g)**.5;return E(pairs_sampled=U,mean_pair=I(L),est=I(T),se=I(Z),pair_part=I(M*L),e3_part=I(N*S),log2_est=I(log2(T),3),log2_bound=I(log2(T+3*Z),3))
def DI(O,A):
	D=A.sym;pt_bmem(O,A);MEM.add(SPB-4096,SPB+8192)
	for G in B(C(R)):R[G]=0
	R[28]=ORBA;R[30]=STOPA;R[32]=SPB;MEM.wr(SPB+D[Bt],2,4660);F=A.run('pt_frame_init')==d;F&=N(MEM.rd(SPB+A,8)==B&M for(A,B)in BU(O,D))and MEM.rd(SPB+D[Bt],2)==0;F&=[R[16],R[17],R[26],R[27],R[29],R[28]]==[Aa,PTB,O[p][0],O[q][0],O[c][1][0]&M,DIRB];A.use();return E(ok=Q(F),slots=C(BU(O,D))+1)
Ac=29184<<32,29440<<32
def Aw(A,fl=8,regs=(5,300),full=1):
	V=regs;N=full;G=fl;F=A.sym;O=F[z];J=F[Bq];D=F['SLOTN'];P=F[BE];Y={};MEM.reg=[];MEM.pg={};l=28*O+2;S=A3+4*l;I=1;K=lambda l:Y.setdefault(l,C(Y)+1)
	def n(a):
		B=a-A3>>2
		if B>=28*O:return 4*K(A.lw[A.lp[Bs]+B-28*O][0])
		C,E=b(B,O);C,F=b(C,7);G,H=b(C,2);D=A.cte(8+G,H,F,E+1);return 4*K(D)if D else 0
	MEM.add(A3,S,n);MEM.add(S,S+16,lambda a:A.by[a-S])
	for(Z,W)in((GCA-4096,16384),(STKB,8192),(m,4194304),(Ac[0],512*P),(Ac[1],512*P)):MEM.add(Z,Z+W)
	for(W,o)in((C1,m),('G_TB8',Ac[0]),('G_TB9',Ac[1])):MEM.wr(GCA+F[W],8,o)
	for L in B(C(R)):R[L]=0
	R[30]=STOPA;R[32]=STKB+4096;R[0]=GCA;R[1]=G-8;A.clear();a=Cl if N else'pt_tb_fb';I&=A.run('pt_tb_build',stops=(a,))==a;X=random.Random(G)
	for M in B(8):
		for c in B(2):
			for d in B(D):I&=MEM.rd(m+8*((M*2+c)*D+d),8)==A3+4*K(Au(A,G,M,0,d,c))
	e=0
	if N:
		for M in B(8):
			for L in B(2*D):
				for T in(0,X.randrange(J),X.randrange(J),J-1):I&=MEM.rd(m+65536+8*((M*2*D+L)*J+T),8)==A3+4*K(Au(A,G,M,T%3,L%D,L//D));e+=1
		f,p=R[0],R[1]
		for U in V:
			g=[0]
			def q():
				g[0]+=1
				if g[0]==1:R[5]=U;R[0]=f+U*P
				else:raise AX(Cm)
			R[0]=f;R[1]=p;I&=A.run(Cl,hooks={'.Lrg_r':q})==Cm
			for h in B(3):
				r=U>>3*(2-h)&7
				for s in B(800):i=X.randrange(F[AR]);j,T=b(i,J);I&=MEM.rd(Ac[0]+U*P+8*(h*F[AR]+i),8)==A3+4*K(Au(A,G,r,T%3,j%D,j//D))
	k=A.counts();A.use();return E(ok=Q(I),instr=H(k.values()),counts=E((A,k.get(A,0))for A in('tb_driver','tb_chb3','tb_fb','tb_regions')),chb3_entries=16*D,fb_entries_checked=e,regions=Ay(V)if N else[],region_entries_checked=2400*C(V)if N else 0,flavour=G)
BZ=K(I(1764*569.23**2*SB.G_SC))
def AM(txt,seg=G,mac=G):
	H=re.compile(Cn);E=[];C=G;D=0
	for F in txt.split(BB):
		A=F.strip()
		if A[:6]==BC:C=A.split()[1];continue
		if A[:6]==BD:C=G;continue
		if A.startswith('.macro '):D=A.split()[1]
		elif A=='.endm':D=0
		B=F.split(Ai,1)[0].rstrip()
		if(seg and C==seg or mac and D==mac)and B.strip()and not B.endswith(':')and not B.strip().startswith('.')and H.match(B):E.append(B.strip())
	return E
def CN(txt):
	F=re.compile(Cn);C={};A=G
	for E in txt.split(BB):
		D=E.strip()
		if D[:6]==BC:A=D.split()[1];C.setdefault(A,0);continue
		if D[:6]==BD:A=G;continue
		B=E.split(Ai,1)[0].rstrip()
		if A is G or not B.strip()or B.endswith(':')or B.strip().startswith('.'):continue
		if F.match(B):C[A]+=1
	return C
def Ax(d,fn='kernels_pt.s',txt=G,H=G):
	K=txt;H=H or W.HT_NEW;C=K if K is not G else AI.ko(d+'/'+fn).read();B=CN(C);L=W.count_kernel_segments(d+'/'+fn)if K is G else G;R=AM(C,Ak);M=1+Az(A for(A,B)in h(R)if B.startswith(B8));O=B[Cb]+B[Cc]+B[Cd];A=E(row=B[Al],node_L2=B[Ak],last_tail=B[AT],node=B[Ak]+B[AT],hand_segments=O,sentinel=M,exit=B[Bw],jump_in=At,hand=At+O+M+B[Bw],empty_patch=B[Ak]+B[Al]+B[AT]);A[Co]=AM(C,mac='ROWDISP')==AM(C,Al)and AM(C,AT)==['b .Lp_top2']and AM(C,Ca)==AM(C,Al)and Cp in C;F=I=0
	for P in AI.ko(d+Cq):
		J=P.strip()
		if J.startswith(BC):I=J.split()[1]=='body_3';continue
		if J.startswith(BD):I=0
		if I and J.startswith('stp x6, x7, [sp, #336]'):F=1
		if I and F and re.match('\\s+[a-z][a-z0-9.]*\\b',P.split(Ai,1)[0].rstrip()):F+=1
	A[C4]=F-1 if F else 0;A[Cr]=W.fp_prices_from_kernels(d+Cq)[0][2]-A[C4]+A[Bx];A[Cs]=L is G or N(L.get(A,0)==B for(A,B)in B.items());A['charged']=E(row=H[AW],cand=H[AV],hand=H[AF],source='v20_lat_b2.HT_NEW');A[D]=Q(A[AU]==H[AW]and A[Bv]==H[AV]and A[Bx]==H[AF]and A['empty_patch']==A[Bv]+A[AU]and A[C4]==W.PT_OLD_HAND and A[Cr]==W.FP_CL[2]-W.PT_OLD_HAND+H[AF]and A[Cs]and A[Co]);A['segment_ops']={aa:ab for(aa,ab)in a(B.items())if aa[:2]in('ha','no','ro','l3','nu')};return A
def DJ(t,old,new,n=1):A=old;assert t.count(A)==n,A;return t.replace(A,new)
DK=('node_L2: wrong ubfx field of f1 (directory key)',L,' ubfx x7, x21, #CW, #NF\n',' ubfx x7, x21, #(CW+1), #NF\n'),('hand-off: msub -> madd (sign of the level-1 centre)',L,' msub x21, x10, x5, x7\n',' madd x21, x10, x5, x7\n'),('node_L2: C1 recurrence with the wrong sign',L,' sub x21, x21, x7\n sub x22, x22, x8\n',' add x21, x21, x7\n sub x22, x22, x8\n'),('node_L2: c0 recurrence uses MU21 instead of MU20 (wrong recurrence register)',L,' sub x22, x22, x8\n',' sub x22, x22, x7\n'),('node_L2: node base built from the wrong register',L,' madd x24, x3, x8, x20\n',' madd x24, x3, x8, x23\n'),('node_L2: directory index f1 field one bit too high (wrong key bit)',L,' bfi x6, x7, #DIRF1SH, #NF\n',' bfi x6, x7, #(DIRF1SH+1), #NF\n'),('node_L2: t0 = d0 >> 31 (wrong asr)',L,' asr x5, x8, #32\n',' asr x5, x8, #31\n'),('hand-off: wheel state of the prefix (U3) ignored',L,' ldrh w8, [sp, #FR_U3]\n',' mov x8, xzr\n'),('chain tails swapped (last <-> non-last)',L,' .if \\ell == 0\n ROWDISP\n',' .if \\ell == 1\n ROWDISP\n'),('wrong chain tail: the last chain branches to the row dispatch',L,Cp,' b .Lp_row\n .endif\n .endm\n .macro LBLT'),('chain advance ADV 2: shift by 2 replaced by 4 (wrong stride)',L,' add x2, x2, x26, lsl #1\n add x5, x5, x27, lsl #1\n .elseif \\s == 3',' add x2, x2, x26, lsl #2\n add x5, x5, x27, lsl #1\n .elseif \\s == 3'),('chain family map (famof_tab): masks 4 and 5 swapped (wheel pattern)',L,' .byte 0, 4, 2, 6, 1, 5, 3, 255\n',' .byte 0, 4, 2, 6, 5, 1, 3, 255\n'),('chain-start table CHAFAM: family-1 tower choice flipped',L,' .elseif \\fam == 1\n .if _n <= NTWA\n .if (_n & 1) == 0\n',' .elseif \\fam == 1\n .if _n <= NTWA\n .if (_n & 1) == 1\n'),('E9 flush stub: pointer reset one word short (E9, forced flush)',C5,' add x9, x0, #16\n str x9, [x8]\n',' add x9, x0, #8\n str x9, [x8]\n'),('builder: LUT index +1 (T[i] convention off by one)',An,' umulh x4, x3, x22\n ldr x4, [x24, x4, lsl #3]\n',' umulh x4, x3, x22\n add x4, x4, #1\n ldr x4, [x24, x4, lsl #3]\n',2),('builder: BL initial value loses its -1 (outward rounding)',An,' sub x5, x27, x2\n sub x5, x5, #1\n mov x3, #CELLW',' sub x5, x27, x2\n mov x3, #CELLW'),('builder: row key adds n0 instead of lo',An,' madd x6, x7, x12, x6\n add x6, x6, x2\n stp x4, x5, [x22], #16\n',' madd x6, x7, x12, x6\n add x6, x6, x7\n stp x4, x5, [x22], #16\n'),('builder: dropped null row (shared empty entry loses the last flag)',An,' mov x3, #(KKEY*(NULLSLOT+LASTADD))\n str x3, [x2, #16]\n',' mov x3, #(KKEY*NULLSLOT)\n str x3, [x2, #16]\n'),('builder: E3 list bin index shift PB+1','t3e',' lsr x2, x2, #PB\n cmp x2, x25\n',' lsr x2, x2, #(PB+1)\n cmp x2, x25\n'),('builder: wheel record block c=1 mixes the masks',C6,' lsl x2, x23, #6\n add x2, x2, x24, lsl #3\n add x2, x2, x22\n',' lsl x2, x23, #6\n add x2, x2, x22, lsl #3\n add x2, x2, x24\n'),('builder: P2 multiple of AE2 for c2=1,2 swapped',C6,' csel x7, x8, x9, eq\n',' csel x7, x9, x8, eq\n'),('builder: rec_build region base always from the emit8 table (flavour ignored)','t3w9',' ldr x15, [x3, x2, lsl #3]\n',' ldr x15, [x3]\n'),('gate: EK2 comparison inverted',n,' mov x0, #5\n b.lo .Lck_ret\n',' mov x0, #5\n b.hs .Lck_ret\n'),('row dispatch: key load advances 16 bytes',L,'.Lp_row:\n ldp x2, x5, [x19], #16\n ldr x6, [x19], #8\n','.Lp_row:\n ldp x2, x5, [x19], #16\n ldr x6, [x19], #16\n'),('global: M3X table value step 728 instead of 729','t3g1',' add x5, x5, #729\n cmp x5, #2187\n',' add x5, x5, #728\n cmp x5, #2187\n'),('global: MB table increment 8 instead of 24','t3g2',' add x6, x6, #24\n sub x8, x6, #72\n',' add x6, x6, #8\n sub x8, x6, #72\n'),('global: chain table region block r=2 copies the b1 block','t3t2',' csel x10, x8, x10, eq\n',' csel x10, x7, x10, eq\n'),('gate: K1 limit raised by 64',n,' mov x2, #CHK_K1\n',' mov x2, #(CHK_K1+64)\n'),('hand-off count off by one (an extra instruction in handoff3a)',BK,' ldr x13, [sp, #FR_E3]\n',' nop\n ldr x13, [sp, #FR_E3]\n'),('candidate price 27 -> 28 (an extra instruction in node_L2)',BK,' ldr w6, [x14], #4\n tbnz w6, #31, .Lp_exit',' nop\n ldr w6, [x14], #4\n tbnz w6, #31, .Lp_exit')
def DL(d,O,nodes,keys,cells,bins,e,n1,fl_nodes=G):
	P=nodes;V=AI.ko(d+BL).read();K=[];k=G
	for(R,B,a,b,*c)in DK:
		f=time.time();M=DJ(V,a,b,*c);J=E(name=R,test=B)
		try:
			if B==BK:A=Ax(d,txt=M)[D]
			else:
				F=A4(M,ext=(AG,A1),sink=B in(L,C5))
				if B==L:A=Z(O,F,P,8,stop=1)[D]
				elif B==C5:A=Z(O,F,P[:2],9,p0=EBB+16384-24,stop=1)[D]
				elif B==An:A=Av(O,F,keys,stop=1)[D]
				elif B=='t3e':A=BV(O,F,cells,bins)[D]
				elif B==C6:A=BW(O,F,only=1 if'block c=1'in R else 2)[D]
				elif B=='t3t':A=Aw(F,8,())[D]
				elif B=='t3t2':A=Aw(F,8,(300,))[D]
				elif B=='t3w9':A=BW(O,F,9,only=1)[D]
				elif B=='t3g1':A=BY(F,2,1)[D]
				elif B=='t3g2':A=BY(F,2,2)[D]
				else:A=BX(O,F)[D]
		except(Fault,KeyError,CQ,IndexError,TypeError,StopIteration,AssertionError)as g:A=y;J['exception']=type(g).__name__
		J[C7]=not A;J['sec']=I(time.time()-f,3);K.append(J)
	S=A4(V,ext=(AG,A1),sink=1);A=Z(O,S,P,8,tmut=1,stop=1)[D];K.append(E(name='chain table ignores the wheel rule (bad_f3): test-sensitivity control on the table spec',test=L,caught=not A));X=E(W.HT_NEW)
	for(T,M)in((AW,7),(AV,28),(AF,32)):Y=E(X);Y[T]=M;K.append(E(name='charged %s %d -> %d'%(T,X[T],M),test=BK,caught=not Ax(d,H=Y)[D]))
	for(R,h,i)in(('spec: rho-bin top of the patch rows replaced by the bin bottom (outward rounding dropped)',2,0),('spec: E3 patch-bin index one bin too low (one-bin margin dropped)',3,0),('spec: cell extents of the patch rows dropped (cells as points): INFORMATIONAL: not resolved by sampled nodes, the exhaustive key-space scans are the evidence for the cell extents',1,1)):J=Ao(O,n1,e,h,kv=pt_kv(S));K.append(E(name=R,test=A2,caught=not J[D],missed=J[AC],informational=Q(i)))
	J=Ao(O,n1,e//4,kv=pt_kv(S));K.append(E(name='certificate weakened to eps/4',test=A2,caught=not J[D],over_eps=J[AD]));j=E(O,DA=[W._mr_da(ab[:6])for ab in O['B']],DB=[W._mr_db5(ac[6:])for ac in O['B']],A0=W.A_OFF_MR,B0=W.B_OFF_MR<<5);K.append(E(name='descriptor in the units of v20_lat_b2.py',test=L,caught=not Z(j,S,P,8,stop=1)[D]));U=[aa for aa in K if not aa.get('informational')];return E(mutants=K,n=C(U),n_caught=H(A[C7]for A in U),all_caught=N(A[C7]for A in U))
def pt_gen(pb):
	import fd_bounds as O;C=W.M;A=(1<<pb)+1
	while 1:
		if C._isp(A):
			D=W._find_root_h(A)
			if D is not G:
				P=W._eval_poly(C.W_COEFFS,D,A);J=-((C.U_BN254+1)%A*P+2)%A
				if J!=0:
					K=W._sqrt_mod(pow(J,-1,A),A)
					if K>0:break
		A+=2
	Q=pow(3,-1,A);E=[1]
	for R in B(5):E.append(E[-1]*D%A)
	L=[H(C.ZK_NUM[aa][A]*E[A]for A in B(6))*Q%A for aa in B(6)];M=L+[-K*ab%A for ab in L];S=pow(M[0],-1,A);F=[[A]+[0]*11]
	for N in B(1,12):
		I=-M[N]*S%A
		if I>A//2:I-=A
		F.append([I if ac==0 else 1 if ac==N else 0 for ac in B(12)])
	R,T,U=W.DESC.bareiss_lll_12d(F,O.G6,1764,1000,y,9,10,(99,100),True);return A,F,T,U,BZ
def DM(orbs,k):A=orbs;D=a(B(C(A)),key=lambda i:A[i][0]);E=C(D);return[A[D[S(E-1,(2*aa+1)*E//(2*k))]]for aa in B(k)]if k>1 else[A[D[E//2]]]
def DN(A):D=A.sym[z];E=A.mac['CHAFL'][1];F=Az(A for A in A.mac[CT][1]if A[0]=='r');return Q([aa[1]for aa in A.dq]==[('8','0'),('8','1'),('9','0'),('9','1')]and[ab[1]for ab in E]==['CHAFAM \\f, \\ell, %d'%ac for ac in B(7)]and O(F[1],A.sym)==D and A.da.get('ct_base')==(0,0,0)and A.da.get(Bs)==(0,0,4)and A.da.get(CV)==(2,0,4)and C(A.by)==8 and A.by[7]==255)
def CO(A):return N(A.sym.get(B)==C for(B,C)in V.items())and DN(A)
def AN(t,q,ks):
	for A in ks:t[A]=t.get(A,0)+q[A]
def CP(orbs):A=E(pt_gate_pct=PT_GATES[PT_G]/100,stor_bound_eps_or_max_pct=float(SB.EPS_OR_MAX)*100,ok=Fr(PT_GATES[PT_G],10000)==SB.EPS_OR_MAX);return A,E(orbit_T_is_production=N(A[4]==BZ for A in orbs),ok=N(A[4]==BZ for A in orbs))
def pt_light(d,orbs=G):
	try:return DO(d,orbs)
	except A as B:return E(ok=y,error=type(B).__name__+': '+T(B)[:300],prices=E(ok=y))
def DO(d,orbs=G,n1=15,n2=4):
	M=orbs;import fd_bounds as l;e=l.main()['delta1_levels_11_to_2'];assert Fr(PT_D1,1<<64)>=Fr(e)*(1+Fr(1,1<<50)),'PT_D1 below fd_bounds delta1';F=time.time();d=T(d);m=M is G and C(W.PT_ORB)==256;M=M or Ay(W.PT_ORB)or[pt_gen(ab)for ab in(47,54,60)];A=E(geometry=E(V),orbits=C(M));J={};A['fd_delta1']=e;A[n],A[o]=CP(M);A[C8]=p=CC(M);J[C8]=time.time()-F;F=time.time();A[C9]=q=Ax(d);J[C9]=time.time()-F;F=time.time();K=A4(AI.ko(d+BL).read(),ext=(AG,A1),sink=1);r=pt_kv(K);J['asm']=time.time()-F;F=time.time();H=E(nodes=0,E=0,N=0,miss=0,over=0,rows=0,cands=0,maxover_ppm=0,empty_nodes=0,cands_P=0,rows_P=0,rows_exec=0,per_orbit=[])
	for(f,u)in h(DM(M,3)):
		O=pt_orb(u);A9,(b,g),w=pt_cert(O['d'],O[Ae],O[o],gate=1<<64);j=CE(O,n1,1+f);X=Ao(O,j,b,kv=r);AN(H,X,(v,U,'N',AC,AD,x,s,Bc,B6,B5,B4));H[t]=P(H[t],X[t]);H[i].append(E(log2p=I(log2(O['p']),2),eps_cert_pct=I(b*100/g,4),E=X[U],size_biased=[X[Bb][ad]for ad in(x,s,'extra')],ok=X[D]))
		if f==1:R,y,c,A0=O,b,j,w
	H[D]=H[AC]==0 and H[AD]==0 and H[U]>0 and N(A[D]for A in H[i]);A[A2]=BP(H);A['orbit']=E(log2p=I(log2(R['p']),3),eps_cert_pct=I(y*100/g,5),checks=A0);J[A2]=time.time()-F;F=time.time();A[L]=A3=Z(R,K,c[:n2],8);A[L][Ct]=A5=CO(K);A['t2_e9']=A6=Z(R,K,c[:2],9,p0=EBB+16384-24);J[L]=time.time()-F;F=time.time();Ab[z]=7
	try:A[CA]=k=Z(R,K,c[:2],8);k[Cu]=7
	finally:Ab.clear()
	J[CA]=time.time()-F;F=time.time();Y=random.Random(5);A7=[(Y.choice([0,1,Y.randrange(1024),1023]),Y.randrange(1<<DP(K)))for aa in B(2)];S=A4(AI.ko(d+BL).read(),ext=(AG,A1));A[AH]=E(gate=BX(R,S),pairs=Av(R,S,A7),e3=BV(R,S,[Y.randrange(1<<S.sym[AS])],[Y.randrange((1<<S.sym[AE])+1)for ac in B(8)]));A[AH][D]=N(A[D]for A in A[AH].values());J[AH]=time.time()-F
	if m:import pt_pop as A8;F=time.time();A[BM]=A8.pp_light(d,K.sym[Cj]);J[BM]=time.time()-F
	A[Cv]=a(K.used|S.used);A[Cw]=E((A,I(B,3))for(A,B)in J.items());A[D]=Q(p[D]and q[D]and H[D]and A3[D]and A6[D]and k[D]and A5 and A[n][D]and A[o][D]and A[AH][D]and(A[BM][D]if A.get(BM)else 1));return A
def DP(A):return A.sym[A0]
def pt_full(d,orbs=G,hb=5):
	try:return DQ(d,orbs,hb=hb)
	except A as B:import traceback as C;return E(ok=y,error=type(B).__name__+': '+T(B)[:300],where=C.format_exc()[-400:])
def DQ(d,orbs=G,seed=7,hb=5):
	r=hb;q=orbs;O=time.time();d=T(d);Y={};AO=q is G;A5=G;u=q;r=S(P(r,1),16)
	if AO:q=[pt_gen(ac)for ac in(46,47,49,51,53,55,57,59,61)];Y['gen_orbits']=time.time()-O;O=time.time();A5,u=W.ht_stream(r);Y['stream']=time.time()-O;O=time.time()
	F=E(geometry=E(V),orbits=C(u),orbits_generated=AO,ht_run=r,orbits_in_T1_T3_T5=C(q),orbits_T1_T3_T5='generated at 2^46..2^61 (pt_gen), as in the first version of this part; the replayed orbits of the run enter T0 and pt_pop');F[n],F[o]=CP(u);F[C8]=Af=CC(u);F[C9]=Ag=Ax(d);Y['cert_prices']=time.time()-O;O=time.time();w=AI.ko(d+BL).read();c=A4(w,ext=(AG,A1),sink=1);A9=pt_kv(c);F[Ct]=Ah=CO(c);F['key_layout']=E(A9);AP=v,U,'N',AC,AD,x,s,Bc,B6,B5,B4;R=E(maxover_ppm=0,per_orbit=[],**{af:0 for af in AP});J=E(E8=E(nodes=0,ok=0,leaves=0,cands=0,rows=0,E=0,dyn=0,bad=0,faults=0,ball_missing=0,flushes=0,rows_emitting=0,rows_skip=0),E9=G,orbits=0);g=E(J[f]);j=E(J[f]);X=[];A6=E(t1=0,t2=0)
	for(AA,Aj)in h(q):
		AB=time.time();b=pt_orb(Aj);Au,(AJ,Ak),Ay=pt_cert(b['d'],b[Ae],b[o],gate=1<<64);y=CE(b,14,seed+AA);A=Ao(b,y,AJ,kv=A9);X.append((b,y,AJ));AN(R,A,AP);R[t]=P(R[t],A[t]);R[i].append(E(log2p=I(log2(b['p']),2),eps_cert_pct=I(AJ*100/Ak,4),E=A[U],size_biased=[A[Bb][al]for al in(x,s,'extra')],maxover_ppm=A[t],ok=A[D]));A6[A2]+=time.time()-AB;AB=time.time();A=Z(b,c,y[:4],8);AN(J[f],A,[ag for ag in J[f]if ag!=D]);J[f][D]+=A[D]and 1 or 0;J[Cx]+=A[v]>0
		if AA%2==0:A=Z(b,c,y[:2],9);AN(g,A,[aj for aj in g if aj!=D]);g[D]+=A[D]and 1 or 0
		if AA in(2,6):A=Z(b,c,y[:2],9,p0=EBB+16384-24);AN(j,A,[ak for ak in j if ak!=D]);j[D]+=A[D]and 1 or 0
		A6[L]+=time.time()-AB
	BP(R);R[Cy]=H(D5(X[0][0],A)for A in X[0][1][:2]);R[D]=R[AC]==0 and R[AD]==0 and R[U]>0 and R[Cy]==0 and N(A[D]for A in R[i]);F[A2]=R;Y[A2]=A6[A2];Y['t2_runs']=A6[L];O=time.time();Ab[z]=7
	try:k=E(J[f]);k=E((A,0)for A in k);A=Z(X[1][0],c,X[1][1][:4],8);AN(k,A,[ah for ah in k if ah!=D]);k[D]=A[D];k[Cu]=7
	finally:Ab.clear()
	J[f]['ok_orbits']=J[f].pop(D);g[Cz]=g.pop(D);j[Cz]=j.pop(D);J['E9']=g;J[C_]=j;J['chunked']=k;J[D0]=J.pop(Cx);AK=J[f];J[D]=AK[e]==0 and AK[v]>=30 and J[D0]>=8 and g[e]==0 and g[v]>0 and j[e]==0 and j[BJ]>0 and k[D]and Ah;p=W.HT_NEW;A3=E(identity='ops = (PRICE_HAND-1) + PRICE_CAND c + PRICE_ROW R + E ;  E = 8 L - 2 R_g + 2 R_s + (L if E9) + 8 F',charged=E(PRICE_ROW=p[AW],PRICE_CAND=p[AV],PRICE_HAND=p[AF]),sets={})
	for(AQ,A)in((f,AK),('E9',g),(C_,j)):Al=AQ!=f;A3['sets'][AQ]=E(hand_offs=A[v],cands_c=A[s],rows_R=A[x],records_L=A[BI],rows_emitting_Rg=A[Bz],rows_skip_Rs=A[B_],flushes_F=A[BJ],E=A[U],ops=A[By],ok=Q(A[By]==A[v]*(p[AF]-At)+p[AV]*A[s]+p[AW]*A[x]+A[U]and A[U]==8*A[BI]-2*A[Bz]+2*A[B_]+(A[BI]if Al else 0)+8*A[BJ]))
	A3[D]=N(A[D]for A in A3['sets'].values());F['accounting']=A3;J[D]=J[D]and A3[D];F[L]=J;Y[CA]=time.time()-O;O=time.time();A7,AR=X[0][0],X[-1][0];Am=X[C(X)//2][0];l=random.Random(seed);M=A4(w,ext=(AG,A1));K=E(per_orbit=[],gate=G);AT=M.sym[A0]
	for m in(A7,AR):AU=[(l.choice([0,1,2,l.randrange(1024),l.randrange(1024),1023]),l.randrange(1<<AT))for ad in B(14)];An=a({0,(1<<M.sym[AS])-1,l.randrange(1<<M.sym[AS])});AX=[l.randrange((1<<M.sym[AE])+1)for ai in B(24)]+[1,2,(1<<M.sym[AE])-1];K[i].append(E(log2p=I(log2(m['p']),2),pairs=Av(m,M,AU),e3=BV(m,M,An,AX),wheel=BW(m,M,8+C(K[i]))))
	K[n]=BX(A7,M);K[D1]=BY(M);K[D2]=Aw(M,8,(300,));K[D3]=Aw(M,9,(),0);K['frame']=DI(A7,M);AL=A4(w,ext=(AG,A1),over=E(KKEY=12,OFFKEY=6));K['slow_path_variant']=A8=Av(Am,AL,[(l.randrange(1024),l.randrange(1<<AT))for ae in B(3)]);A8[Bq]=12;A8[CU]=6;K[D]=K[n][D]and K[D1][D]and K[D2][D]and K[D3][D]and K['frame'][D]and A8[D]and A8[Ch]>0 and N(A['pairs'][D]and A['e3'][D]and A['wheel'][D]for A in K[i])and N(A['pairs'][C2]>=14 for A in K[i]);F[AH]=K;Y[AH]=time.time()-O;O=time.time();F['build_estimate']=AY=E(declared_PT_BUILD_IR=W.M.PT_BUILD_IR,declared_log2=I(log2(W.M.PT_BUILD_IR),3),per_orbit=[]);Ap=X[C(X)//2][0]
	for m in(A7,Ap,AR):AZ=DH(m,M);AZ['log2p']=I(log2(m['p']),2);AY[i].append(AZ)
	AY['note']='INFORMATIONAL: build estimates of single orbits (2^46, 2^53, 2^61); PT_BUILD_IR bounds the production-weighted MEAN per orbit (pt_pop.pp_acct); full-orbit builds are not executed in this package';Y['build_est']=time.time()-O;O=time.time();AM,Aa,Aq=X[C(X)//2];Ar=a(Aa,key=lambda n:C(pt_cands(AM,n))if C(pt_cands(AM,n))>=3 else 99)[:2];F['t5']=DL(d,AM,Ar,AU[:1],[0],AX[:5],Aq,Aa[:8]);Y['t5']=time.time()-O;O=time.time()
	if A5:import pt_pop as As;F[CB]=As.pp_part(d,A5,u,r,M,A9,K[i][0]['wheel']['counts'],CN(w),c);Y[CB]=time.time()-O
	c.use();M.use();AL.use();Ac=c.used|M.used|AL.used;Ad={ab.group(1)for aa in w.split(BB)for ab in[re.match('^\\s+([a-z][a-z0-9.]*)\\b',aa.split(Ai,1)[0].rstrip())]if ab and not aa.strip().startswith('.')};F['mnemonics_supported']=DA;F[Cv]=a(Ac);F['mnemonics_in_file']=C(Ad);F['mnemonics_not_executed']=a(Ad-Ac);F[Cw]=E((A,I(B,3))for(A,B)in Y.items());F[D]=Q(Af[D]and Ag[D]and R[D]and J[D]and K[D]and F['t5']['all_caught']and F[n][D]and F[o][D]and(F[CB][D]if A5 else 1));return F
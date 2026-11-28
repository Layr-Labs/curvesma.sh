Ap='switches'
Ao='max_page'
An='b1_page'
Am='b1_term'
AZ='b1_switch'
AY='b1_row'
BT='.Lchunk_loop'
BS='.Lfl_g_'
BR='.rept '
BQ='apply_page'
BP='wheel_disp4'
BO='wheel_skip2_stub'
BN='emit8_step1:'
BM='emit9_step1:'
BL='.Lchain_entry'
BK='.Lchord_split'
BJ='//@end'
BI='//@seg'
AX='ldr_uxtw3'
AW='ldr_sxtw3'
AV='lsl_reg'
AU='mov_imm'
AT='sub_imm12'
AS='add_lsr'
AR='str_sp_off'
AQ='ldr_sp_off'
AP='ldr_lsl3'
AO='ldr_off'
AN='ldp_off'
AM='stp_sp_off'
AL='ldp_sp_off'
AB='staddb'
AA='staddh'
A9='csel_hs'
A8='sub_reg'
A7='sub_imm'
A6='sub_lsl'
A5='add_reg'
A3='add_lsl'
u='str_post'
s='ldr_reg'
r='ldr_post'
q='stp'
m='stp_post'
f='ldp'
e='ldp_post'
BE='kernels_v32.s'
d='L'
A2='term_pages'
A1='b1_line'
A0='entries'
z='pages'
y='min_page'
Al='.Ltop_2'
Ak='.Lk0_top'
Aj='b.pl .Lk0_top'
Ai='wheel_skip1_stub'
Ah='.Lcont_s3'
Ag='.Lcont_s2'
Af='.Lcont_s1'
Ae='.Lfl_last'
b='bcc'
a='tbnz'
Z='tbz'
W='b'
U='adr'
Q='and_imm'
P='cmp_imm'
x='add_imm'
T=abs
Ad=IndexError
Ac=KeyError
Aa=sorted
AK='other'
AJ='b .Ltop_2'
AI='.Lcont_last'
AH='.Lfl_s3'
AG='.Lfl_s2'
AE='.Lfl_s1'
AD='ret'
O='subs_reg'
p='add_imm12'
AF=float
A4='.Lmask0'
t='ok'
w=':'
v=zip
o=tuple
n=pow
l=dict
k=all
j=set
i=bool
h=list
g=str
Y=.0
X=max
c=enumerate
V=True
S=ValueError
R=int
N=False
J=1.
I=round
G=None
E=sum
C=len
A=range
from builtins import eval as C9
import stor_bound as Ab,math,random,re
from bisect import bisect_left
from pathlib import Path
from math import log,log2,sqrt,ceil,exp
import model_bn254_v6 as M,v20_lat_b2 as L,descent_bn254_v6 as DESC
F=(1<<64)-1
def verify_v25_upgrades(seed=20261122):
	O=random.Random(seed);CR=M.load_data()['t4_geom'];BF,S,S,U,BG,CS=M.build_exact_int_geom(CR)
	def BH(K):
		B,E=M.QA_MAX,M.QB_MAX;F=lambda t:20*t**3-45*t**4+36*t**5-10*t**6;N=.05*(F(B)+F(E)-J);O=X(Y,J-B);G=O**3*B**3;P=J-E;H=M.B_A_SHELLS*K;C=J/H
		for D in A(1,H+1):
			I=X(P,B*((D-1)*C)**(J/3.));L=min(B,B*(D*C)**(J/3.))
			if L>I:G+=D*C*B**3*((J-I)**3-(J-L)**3)
		return G/N
	BI=BH(64);BJ=BH(M.K_B1_PREF);CT=BI<=1.031 and BJ<=M.B1_PREF_COV
	for S in A(5):O.sample(A(1,1024),320);[O.randrange(0,32768)for ah in A(1024)]
	u=CQ(g(Path(__file__).resolve().parent/BE),seed=seed+2);BK=lat_apply_e2e(g(Path(__file__).resolve().parent/BE),seed=seed+3);CU=0 if u[t]and BK[t]else 1;BL=u['c_hit_lo_exact'];CV=BL<=M.C_HIT_LO and u['row_static']+u['switch_allow']<=M.C_B1_ROW;AG=0;Ab=0;Ac=0
	for Ad in(42,46,48):
		D=(1<<Ad)+1
		while V:
			if D&7 in(1,3)and M._isp(D):
				Ae=L._find_root_h(D)
				if Ae is not G:
					CW=L._eval_poly(M.W_COEFFS,Ae,D);BM=-((M.U_BN254+1)%D*CW+2)%D
					if BM!=0:
						BN=L._sqrt_mod(n(BM,-1,D),D)
						if BN>0:break
			D+=2
		CX=n(3,-1,D);Af=[1]
		for S in A(5):Af.append(Af[-1]*Ae%D)
		CY=[E(M.ZK_NUM[ai][A]*Af[A]for A in A(6))*CX%D for ai in A(6)];CZ=-n(BN,-1,D)%D;AH=[aj*CZ%D for aj in CY];T=AH[0];m=[];BO=[]
		for DL in A(32):Ca=[O.randint(-300,300)for ar in A(5)];Cb=O.randint(-200,150);v=[[Cb+at]+Ca for at in A(32)];BO.append(v);m.extend(v)
		Cc=lambda v:math.sumprod(L._mv6(BF,v),v);BP=Aa(A(C(m)),key=lambda idx:(Cc(m[idx]),m[idx]));BQ=[0]*C(m)
		for(Z,Ag)in c(BP):BQ[Ag]=Z
		BR=[-1]*C(m);BS=0
		for v in BO:
			Cd=v[0];AI=E(Cd[A]*AH[A]for A in A(6))%D
			for W in v:
				Z=BQ[BS];BS+=1;BR[Z]=AI;w=E(W[A]*AH[A]for A in A(6))%D
				if AI!=w:AG+=1
				Ab+=1;p=AI+T;AI=p-D if p>=D else p
		for(Z,Ag)in c(BP):
			W=m[Ag];w=E(W[A]*AH[A]for A in A(6))%D
			if BR[Z]!=w:Ac+=1
	Ce=569.23*sqrt(M.QA_MAX)*sqrt(196608./753729.);AJ=20000;Cf=E(5.*((A+.5)/AJ)**4*sqrt(J-((A+.5)/AJ)**2)for A in A(AJ))/AJ;BT=2.*Ce*Cf;BU=J+2.*E(((A/AF(M.B_A_SHELLS))**(J/3.))**2.5 for A in A(1,M.B_A_SHELLS));BV=BT/BU
	if not BV>=34.:AG+=1
	BW=U[3]==[0,0,0,0,-1,0]and U[4]==[0,0,0,0,0,-1]and U[5]==[0,0,0,1,0,0];x=0;Ah=0
	for S in A(450):
		q=[O.randint(-300,300)for ak in A(6)];y=[O.randint(-400,400)for al in A(6)];W=q+y;Ai=[];Aj,Ak=q[:],y[:]
		for S in A(3):
			for Cg in(1,-1):
				for Ch in(1,-1):Ai.append([Cg*aI for aI in Aj]+[Ch*aJ for aJ in Ak])
			Aj=L._mv6(U,Aj);Ak=L._mv6(U,Ak)
		Ci=q[3]-q[4]+q[5];Al=L._mv6(U,q)
		if Al[3]-Al[4]+Al[5]!=Ci:x+=1
		T=L._cone_of(W)
		if T>=0:
			BX,BY=L._canon_cone_ref(o(W[:6]),o(W[6:]),U,BG);Cj=h(BX)+h(BY);BZ=[au for au in Ai if au[3]-au[4]+au[5]>0 and au[3]>-au[4]and au[3]>au[5]and au[6]>0]
			if C(BZ)!=1 or BZ[0]!=Cj:x+=1
			for z in Ai:
				Ck,Cl=L._canon_cone_ref(o(z[:6]),o(z[6:]),U,BG)
				if(Ck,Cl)!=(BX,BY):x+=1
			Ah+=1
	Am=R(I(1764*572.**2*CS));Cm=1<<38;AK=Y;An=0;Ao=1<<64
	for Ad in(47,52,56,61):
		H=(1<<Ad)+1
		while V:
			if H&7 in(1,3)and M._isp(H):
				Ap=L._find_root_h(H)
				if Ap is not G:
					Cn=L._eval_poly(M.W_COEFFS,Ap,H);Ba=-((M.U_BN254+1)%H*Cn+2)%H
					if Ba!=0:
						Bb=L._sqrt_mod(n(Ba,-1,H),H)
						if Bb>0:break
			H+=2
		Co=n(3,-1,H);Aq=[1]
		for S in A(5):Aq.append(Aq[-1]*Ap%H)
		Bc=[E(M.ZK_NUM[am][A]*Aq[A]for A in A(6))*Co%H for am in A(6)];Bd=Bc+[-Bb*av%H for av in Bc];Cp=n(Bd[0],-1,H);Be=[[H if aw==0 else 0 for aw in A(12)]]
		for Bf in A(1,12):
			Ar=-Bd[Bf]*Cp%H
			if Ar>H//2:Ar-=H
			Be.append([Ar if aE==0 else 1 if aE==Bf else 0 for aE in A(12)])
		S,a,S=DESC.bareiss_lll_12d(Be,BF,1764,1000,N,9,10,(99,100),V);Cq=((a[0]*Am<<38)+a[1]-1)//a[1];Bg=(a[2]<<38)//(a[1]*Am);Ao=min(Ao,Bg);AL=[0]*12
		for A0 in A(11,-1,-1):Cr=math.isqrt((a[A0]*Am+a[A0+1]-1)//a[A0+1])+1;AL[A0]=Cr+E((AL[A]+1)//2 for A in A(A0+1,12))
		Cs=E(AL[A]for A in A(2,12))+1;As=2*AL[1]+2;Ct=(2*As*Cs*Bg>>32)+2;Cu=Cm+Ct+As*As;Bh=(Cu*Cq>>64)+1
		if Bh>=1<<40:An+=1
		AK=X(AK,log2(Bh))
	r=0;Cv=1<<40;Bi=0
	for Cw in(0,1,17,(1<<34)-1):
		AM=Cv+(Cw<<15)+8;At=AM+16
		for Cx in A(1,4100):
			AN=At;At+=8
			if not(AM+16<=AN and AN+8<=AM+16384):r+=1
			Bj=AN>>14&1;Cy=Cx%2046==0
			if i(Bj)!=Cy:r+=1
			if Bj:
				Bk=AN-16384+8
				if Bk!=AM:r+=1
				At=Bk+16;Bi+=1
	if Bi!=8:r+=1
	AO=0
	for Au in(8,16):
		for AP in(35,42,48,50):
			Cz=1<<AP;AQ=4000;Av=AQ<<17|15;Bl=Aa((O.randint(1,5000),O.randrange(Au))for A in A(200));C_=[ax<<17|ay for(ax,ay)in Bl]
			for Aw in(1,-1):
				Ax={};Ay={};Bm=1<<55 if Aw>0 else(1<<55)+(Au-1)*Cz
				for Az in C_:
					if Az>Av:break
					A1=Az&15;D0=Bm+(A1<<AP)if Aw>0 else Bm-(A1<<AP);A2=D0+(Az>>16);Ax[A2]=Ax.get(A2,0)+720&65535
				for(Z,A1)in Bl:
					if Z>AQ:break
					D1=A1 if Aw>0 else Au-1-A1;A2=(1<<55)+(D1<<AP)+(Z<<1);Ay[A2]=Ay.get(A2,0)+720&65535
				if Ax!=Ay:AO+=1
			if AQ<<17|15>Av or AQ+1<<17|0<=Av:AO+=1
	A_=0
	for Bn in(-64,-63,-1,0,1,63,64,559):
		Bo=Bn>>6
		for Bp in(8,16,32,64):
			Bq=R(I(log2(Bp)));Br=Bn>>Bq<<Bq;D2=Br+Bp-1
			if Br>>6!=Bo or D2>>6!=Bo:A_+=1
	import fd_bounds as D3;AR=D3.main();import exact_levels as D4;Bs=D4.main();B0=AR['eps_or_bound'];A3=0
	if not B0<=M.SC_GUARD_DELTA:A3+=1
	import stor_bound as AS;A4=AS.main()
	if not(A4[t]and AR[t]and B0<=AF(AS.EPS_OR_MAX)and AF(AS.EPS_OR_MAX)<=AF(AS.DELTA)and AR['exact_per_orbit_ir_log2']<=log2(M.C_LAT_EXC)and Bs[t]):A3+=1
	Bt=i(A4['sink_ok']);Bu=Bv=Y;Bw=551.94
	for D5 in A(120):
		Bx=((D5+.5)/12e1)**.4;By=sqrt(J-Bx)*Bw
		for D6 in A(8):
			B1=8.*(D6+.5)/8.;Bz={}
			for y in A(ceil(B1-By),R(B1+By)+1):W=X(Y,J-Bx-((y-B1)/Bw)**2)**3;Bz.setdefault(y>>3,[]).append(W)
			for B_ in Bz.values():Bu+=8.*X(B_);Bv+=E(B_)
	C0=Bu/Bv
	if not C0<=1.05:A3+=1
	B2=0;B3=256;D7=[32768-O.randint(2000,8000)for aa in A(B3)];b=[D7[:]for ab in A(8)];C1=j()
	for P in A(8):
		for AT in A(B3):
			if O.random()<.03:b[P][AT]=b[P][AT]+9000&65535
			if b[P][AT]&32768:C1.add((P,AT))
	C2=j()
	for T in A(0,B3,4):
		C3=0
		for P in A(8):D8=b[P][T]|b[P][T+1]<<16|b[P][T+2]<<32|b[P][T+3]<<48;C3|=D8
		if C3&0x8000800080008000:
			for P in A(8):
				for C4 in A(4):
					if b[P][T+C4]&32768:C2.add((P,T+C4))
	if C2!=C1:B2+=1
	from fractions import Fraction as A5;D9=[A5(ac)for ac in M.H_COEFFS];B4=[[A5(an,3)for an in ad]for ad in M.ZK_NUM];B=L.PHI_F3
	def DA(a,b):
		B=[A5(0)]*(C(a)+C(b)-1)
		for(F,G)in c(a):
			for(H,I)in c(b):B[F+H]+=G*I
		while C(B)>6:
			D=B.pop()
			if D:
				J=C(B)-6
				for E in A(6):B[J+E]-=D*D9[E]
		return B+[A5(0)]*(6-C(B))
	K=[B4[ae][:]+[A5(1 if ae==az else 0)for az in A(6)]for ae in A(6)]
	for Q in A(6):
		C5=next(A for A in A(Q,6)if K[A][Q]!=0);K[Q],K[C5]=K[C5],K[Q];DB=K[Q][Q];K[Q]=[ao/DB for ao in K[Q]]
		for A6 in A(6):
			if A6!=Q and K[A6][Q]:DC=K[A6][Q];K[A6]=[K[A6][aF]-DC*K[Q][aF]for aF in A(12)]
	DD=[af[6:]for af in K];AU=0
	for d in A(6):
		for e in A(6):
			DE=DA(B4[d],B4[e]);C6=[E(DE[A]*DD[A][aA]for A in A(6))for aA in A(6)]
			if any(A.denominator!=1 for A in C6):AU+=1
			else:
				DF=[aG.numerator for aG in C6]
				for B5 in A(3):
					if E(B[B5][A]*DF[A]for A in A(6))%3!=B[B5][d]*B[B5][e]%3:AU+=1
	C7=(B[0][0]*(B[1][1]*B[2][2]-B[1][2]*B[2][1])-B[0][1]*(B[1][0]*B[2][2]-B[1][2]*B[2][0])+B[0][2]*(B[1][0]*B[2][1]-B[1][1]*B[2][0]))%3;DG=C7!=0 and k(B[A][0]==1 for A in A(3));B6=[B.index(o(E(B[ag][A]*U[A][D]for A in A(6))%3 for D in A(6)))for ag in A(3)];DH=Aa(B6)==[0,1,2]and B6!=[0,1,2];C8=E(1 for B in A(729)if any(B//3**A%3==0 and B//3**(3+A)%3==0 for A in A(3)));C9=E(1<<3-E(1 for A in A(3)if B//3**A%3==0)for B in A(27));B7=0;CA=0
	for S in A(250):
		B8=[O.randint(-40,40)for ap in A(6)];B9=[O.randint(-40,40)for aq in A(6)]
		if not any(B8)or not any(B9):continue
		if L._bad_f3(B8+B9):
			CA+=1
			if M.is_ideal_primitive_oh(B8,B9):B7+=1
	BA=0;AV=[G]*448
	for A7 in A(1,7):
		for A8 in A(1,61):
			BB=A7<<6|A8;f=[aB for aB in A(A8)if not A7>>aB%3&1]
			if not f:AV[BB]='NONE',0,0
			elif C(f)==1:AV[BB]='TAIL1',f[0],1
			else:CB=[f[aH+1]-f[aH]for aH in A(C(f)-1)];A9='C3'if k(A==3 for A in CB)else'C12'if CB[-1]==1 else'C21';AV[BB]=A9,f[0],C(f)
	AW=[(0,2880),(3840,3864)]+[(8192*aC,8192*aC+5554)for aC in A(1,12)]+[(131072,139264)];AW.sort();CC=E(1 for A in A(C(AW)-1)if AW[A][1]>AW[A+1][0])
	for A7 in A(1,7):
		for CD in h(A(0,64))+[120,121,125]:
			AX=[];CE,CF,AA=1000,2000,CD
			while AA>0:
				A8=60 if AA>60 else AA;A9,CG,e=AV[A7<<6|A8]
				if e>0:
					AB,AC=CE+CG*7,CF+CG*13
					if A9=='TAIL1':AX.append((AB,AC))
					elif A9=='C3':
						for d in A(e):
							AX.append((AB,AC))
							if d+1<e:AB+=21;AC+=39
					else:
						BC=1 if A9=='C12'else 2;AD=[]
						for S in A(e-1):AD.append(BC);BC=3-BC
						AD.reverse()
						for d in A(e):
							AX.append((AB,AC))
							if d<C(AD):AB+=AD[d]*7;AC+=AD[d]*13
				AA-=A8
				if AA>0:CE+=420;CF+=780
			w=[(1000+aD*7,2000+aD*13)for aD in A(CD)if not A7>>aD%3&1]
			if AX!=w:BA+=1
	AE=0
	for CH in(-3145728,-100000,-17,-1,0,1,2,3,17,100000,3145728):
		z=CH+(768<<12)&4294967295;DI=z*1431655766>>32;DJ=z-3*DI
		if DJ!=CH%3:AE+=1
	for CI in(0,1,2,50,400):
		AY=2369;s=CI<<12|AY;CJ=0
		while V:
			if s&4095!=AY:AE+=1
			p=2048+((AY-2048)*5+7)%729;s=s&~4095|p;s=s-4096&F;AY=p;CJ+=1
			if s>>63&1:break
		if CJ!=CI+1:AE+=1
	AZ=0;BD=ceil(2.**A4['score_node_log2']/2.**38);CK=R(2.**A4['Ea_plus_log2'])+1
	if not(BD<1<<26 and(BD-1<<38)+CK+1<1<<64):AZ+=1
	for CL in(0,1,(1<<34)-1,(1<<35)-1):
		for CM in(0,1,BD-1):
			for CN in(1,CK):
				CO=CL<<29|CM
				if CO>>26!=8*CL or CN+(CO<<38&F)&F!=(CM<<38)+CN:AZ+=1
	CP=AU==0 and DG and DH and C8==217 and C9==125 and B7==0 and CA>=40 and BA==0 and CC==0 and AE==0 and AZ==0;DK=CP and CT and CV and CU==0 and AG==0 and Ac==0 and Ab==3072 and BW and x==0 and Ah>=300 and An==0 and AK<25. and r==0 and AO==0 and A_==0 and A3==0 and Bt and B2==0;return l(ok=DK,b1_cov_k64_exact=I(BI,5),b1_cov_k256_exact=I(BJ,5),b1_pref_cov_budget=M.B1_PREF_COV,b1_c_hit_lo_exact=I(BL,6),b1_c_hit_lo_budget=M.C_HIT_LO,b1_c_b1_row_budget=M.C_B1_ROW,b1_e2e=u,lat_apply_e2e=BK,a0_chord_kta_tested=Ab,a0_chord_kta_errs=AG,a0_chord_ai_order_errs=Ac,a0_chord_mean_points=I(BT,2),a0_chord_mean_pieces=I(BU,4),a0_chord_piece_mean_points=I(BV,2),cone3_s35_alg_ok=BW,cone3_bijection_tested=Ah,cone3_bijection_errs=x,fd_rem0_max_r2_bits=I(AK,2),fd_rem0_u40_errs=An,tbnz14_single_slot_errs=r,sh17_x4_var_ks_errs=AO,b0_floor_div_nest_errs=A_,over_radius_eps_bound=I(B0,7),fd_bounds=AR,exact_levels=Bs,scan_group_overhead=I(C0,4),stor_bound=A4,fd_A1_min_log2=I(log2(Ao),2),row_capacity_guard_errs=A3,dummy_sink_blocks_c0=1<<15,dummy_sink_disjoint_ok=Bt,init_scan8_readonly_errs=B2,v28_wheel_ok=CP,v28_wheel_hom_errs=AU,v28_wheel_det3_f3=C7,v28_wheel_sigma_perm=B6,v28_wheel_bad_729=C8,v28_wheel_rep_27=C9,v28_wheel_false_discards=B7,v28_wheel_wtab_errs=BA,v28_sp_alias_errs=CC,v28_wheel_x23_barrett_errs=AE,v28_emit10_errs=AZ)
D=re.compile
B='x(\\d+)'
H='#(\\d+)'
Aq=[(e,D(rf"ldp {B}, {B}, \[{B}\], {H}$")),(f,D(rf"ldp {B}, {B}, \[{B}\]$")),(m,D(rf"stp {B}, {B}, \[{B}\], {H}$")),(q,D(rf"stp {B}, {B}, \[{B}\]$")),(r,D(rf"ldr {B}, \[{B}\], {H}$")),(s,D(rf"ldr {B}, \[{B}, {B}\]$")),(u,D(rf"str {B}, \[{B}\], {H}$")),(A3,D(f"add {B}, {B}, {B}, lsl {H}$")),(A5,D(f"add {B}, {B}, {B}$")),(p,D(f"add {B}, {B}, {H}, lsl #12$")),(x,D(f"add {B}, {B}, {H}$")),(A6,D(f"sub {B}, {B}, {B}, lsl {H}$")),(A7,D(f"sub {B}, {B}, {H}$")),(A8,D(f"sub {B}, {B}, {B}$")),(O,D(f"subs {B}, {B}, {B}$")),(P,D(f"cmp {B}, {H}$")),('cmp',D(f"cmp {B}, {B}$")),(A9,D(f"csel {B}, {B}, {B}, hs$")),(Q,D(f"and {B}, {B}, #0x([0-9a-f]+)$")),('lsr',D(f"lsr {B}, {B}, {H}$")),('ubfx',D(f"ubfx {B}, {B}, {H}, {H}$")),(U,D(f"adr {B}, ([.\\w]+)$")),('blr',D(f"blr {B}$")),('br',D(f"br {B}$")),(W,D('b ([.\\w]+)$')),(AD,D('ret$')),(Z,D(f"tbz {B}, {H}, ([.\\w]+)$")),(a,D(f"tbnz {B}, {H}, ([.\\w]+)$")),(b,D('b\\.(eq|hi|lo|ls|pl) ([.\\w]+)$')),(AA,D(rf"staddh w(\d+), \[{B}\]$")),(AB,D(rf"staddb w(\d+), \[{B}\]$")),(AL,D(rf"ldp {B}, {B}, \[sp, {H}\]$")),(AM,D(rf"stp {B}, {B}, \[sp, {H}\]$")),(AN,D(rf"ldp {B}, {B}, \[{B}, {H}\]$")),(AO,D(rf"ldr {B}, \[{B}, {H}\]$")),('ldr',D(rf"ldr {B}, \[{B}\]$")),(AP,D(rf"ldr {B}, \[{B}, {B}, lsl #3\]$")),(AQ,D(rf"ldr {B}, \[sp, {H}\]$")),(AR,D(rf"str {B}, \[sp, {H}\]$")),('str',D(rf"str {B}, \[{B}\]$")),('ldadd',D(rf"ldadd {B}, {B}, \[{B}\]$")),(AS,D(f"add {B}, {B}, {B}, lsr {H}$")),(AT,D(f"sub {B}, {B}, {H}, lsl #12$")),(AU,D(f"mov {B}, {H}$")),('extr',D(f"extr {B}, {B}, {B}, {H}$")),('sbfx',D(f"sbfx {B}, {B}, {H}, {H}$")),('mul',D(f"mul {B}, {B}, {B}$")),('madd',D(f"madd {B}, {B}, {B}, {B}$")),('msub',D(f"msub {B}, {B}, {B}, {B}$")),('udiv',D(f"udiv {B}, {B}, {B}$")),(AV,D(f"lsl {B}, {B}, {B}$")),('bfi_w',D('bfi w(\\d+), w(\\d+), #(\\d+), #(\\d+)$')),(AW,D(rf"ldr {B}, \[{B}, w(\d+), sxtw #3\]$")),(AX,D(rf"ldr {B}, \[{B}, w(\d+), uxtw #3\]$"))]
def AC(path,prefix):
	D,J,F,B,E,H=[],{},[],G,G,[]
	for K in Ab.ko(path):
		I=K.strip()
		if I.startswith(BI):B=I.split()[1];continue
		if I.startswith(BJ):B=G;continue
		A=K.split('//',1)[0].strip()
		if B is G or not B.startswith(prefix)or not A:continue
		if A.startswith('.rept'):E,H=R(A.split()[1]),[];continue
		if A.startswith('.endr'):
			for L in H*E:D.append(L);F.append(B)
			E=G;continue
		if A.endswith(w):J[A[:-1]]=C(D);continue
		if A.startswith('.'):continue
		if E is not G:H.append(A)
		else:D.append(A);F.append(B)
	return D,J,F
def As(path):return AC(path,'b1_')
def lat_apply_e2e(path,seed=20261125):
	x=path;CA,A5,CB=AC(x,'lat_apply');Am=_b1_decode(CA,A5);CC=L.count_kernel_segments(x);An=CC.get('emit_direct_flush',-1);AF=E(1 for A in CB if A=='lat_apply_page');y=A5['apply_end']-A5['apply'];BU=(An+2+20)/2046.;BV=(y+2)/(y//3*2)+(AF+3+20)/2046.;CD=random.Random(seed);CE,Cm,_=AC(x,'emit');Cn=[aa for aa in CE];f,T,A6,Y,Ao={},{},G,G,G;a={AE:[],AG:[],AH:[],Ae:[]};AL=G
	for BW in Ab.ko(x):
		Ap=BW.strip()
		if Ap.startswith(BI):A6=Ap.split()[1];continue
		if Ap.startswith(BJ):A6=G;continue
		D=BW.split('//',1)[0].strip()
		if not D or D.startswith('.')and not D.endswith(w):continue
		if D.endswith(w):
			if Y is not G and Y!=A4:T.setdefault(Y,[]).append(D)
			if Y!=A4:Y=D[:-1]if D[:-1]in(BK,Af,Ag,Ah,AI,BL,A4)else G
			AL=D[:-1]if D[:-1]in a else G
			if D in(BM,BN)and Ao in(Ai,BO):T[Ao]=D
			continue
		if Y is not G:
			T.setdefault(Y,[]).append(D)
			if Y not in(AI,BL,A4)or C(T[Y])>=2:Y=G
		if A6:f.setdefault(A6,[]).append(D);Ao=A6
		if AL:
			a[AL].append(D)
			if D.startswith('b '):AL=G
	AM=T.get(A4,[]);CF=T.get(BK)==['wheel_disp4:']and T.get(Ai)==BM and T.get(BO)==BN and T.get(Af)==['emit9_step2:']and T.get(Ag)==['emit9_step3:']and T.get(Ah)==['emit9_last:']and T.get(AI)==[Aj,AJ]and AM==[Aj,AJ]and AD not in AM;z=f.get(BP,[]);Aq=f.get(Ai,[]);AN=f.get('emit9_step1',[]);Ar=f.get('emit9_step2',[]);As=f.get('emit9_step3',[]);At=f.get('emit9_last',[]);CG=[Aj,AJ];BX=z+Aq+AN+Ar+As+At+CG;BY=C(z);BZ=BY+C(Aq);Ba=BZ+C(AN);Bb=Ba+C(Ar);Bc=Bb+C(As);CH=Bc+C(At);Bd=C(BX);Au=Bd+C(a[AE]);Av=Au+C(a[AG]);Be=Av+C(a[AH]);Aw=BX+a[AE]+a[AG]+a[AH]+a[Ae];Ax=C(Aw);Bf=Ax+1;CI={Af:Ba,Ag:Bb,Ah:Bc,AI:CH,AE:Bd,AG:Au,AH:Av,Ae:Be,Ak:Ax,Al:Bf};Bg=_b1_decode(Aw,CI);J,Ay,Az={},{},{};m=1<<48;A_=1<<50;U=1<<52;B0=1<<46;n=1<<45;A0=n+(1<<20);CJ=1<<40;B1=19;J[B0+480]=B1;o=(1<<34)-1;Bh=0,1,o;p={ab:(1<<54)+32768*ab+8 for ab in Bh}
	for Bi in Bh:J[U+8*Bi]=p[Bi]+16
	B2=4294967291&~(7<<23)|1<<23;CK=B2-(1<<32)if B2>>31 else B2;J[n+(CK<<3)&F]=A0;J[A0+48]=4*BZ;J[A0+56]=4*BY;AO=0;O=V
	try:
		B3=_B1CPU(Bg,J,Ay);H=B3.x;H[1]=8;H[10]=n;H[11]=48;H[12]=3;H[26]=16;H[27]=1;H[17]=U;H[31]=B0;AP=1,(1<<37)-1,1<<38,(1<<37)+1,(1<<47)+12345,(1<<48)-1;Bj=43963475;CL=0,1,3,6;AQ=525
		for d in A(AQ):
			for(CM,CN)in c(CL):
				q=d*4+CM;B4=d*7+CN;A7=1 if q<2050 else 0 if q<2080 else o;Bk=0 if A7==0 else Bj if A7==o else 1+q%3;Bl=AP[(q-2080)%C(AP)]if A7==o else 1+q%256;Bm=A_+(q&3)*CJ;J[m+16*B4]=Bl;J[m+16*B4+8]=Bm;J[Bm+8*B4]=A7<<30|Bk
				if A7==1 and q<2046:A1=(Bk<<38)+Bl;Az[A1]=Az.get(A1,0)+B1&255
		def AR(c):
			nonlocal AO;AO+=1;A=c.x[0]
			if A!=p[1]or J[A]!=2046 or J[A+8]!=B1:raise S('bad flush header')
			if J[U+8]!=p[1]+16 or c.x[1]!=8:raise S('bad ptr reset')
			B=_B1CPU(Am,J,Ay);B.x[28]=A;B.run(A5[BQ],j(),{},{})
		CO={Au-1:AR,Av-1:AR,Be-1:AR,C(Aw)-1:AR};Bn={};H[23]=AQ-1
		for d in A(AQ):H[2]=m+16*(d*7-(1 if d==0 else 0));H[5]=d*7-(1 if d==0 else 0);H[6]=4294967291;H[14]=1;H[0]=7 if d==0 else 6;H[23]=H[23]-1&F;B3.N=i(H[23]>>63);B3.run(0,{Ax,Bf},CO,Bn)
		CP=k(J[p[o]+16+8*A]==(Bj<<38)+AP[A%C(AP)]for A in A(20));O=CF and C(z)==4 and C(Aq)==2 and C(AN)==9 and C(Ar)==9 and C(As)==9 and C(At)==7 and AO==1 and E(Bn.values())==AQ*39+2+1+8 and Ay==Az and CP and H[1]==8 and J[U+8]==p[1]+16+32 and J[U+0]==p[0]+16+240 and J[U+8*o]==p[o]+16+160 and k(A[0]!=AK for A in Bg)
		if O:
			CQ=[ac.split('//',1)[0].strip()for ac in Ab.ko(x)];AS,Bo,B5,AT={},[],N,G
			for D in CQ:
				if not D:continue
				if D=='.altmacro':B5=V;continue
				if D=='.noaltmacro':B5=N;continue
				if not B5:continue
				if D.startswith('.macro '):B6=D[7:].strip().split(G,1);Bp=B6[0];CR=[af.strip()for af in B6[1].split(',')]if C(B6)>1 else[];AT=Bp;AS[Bp]=CR,[];continue
				if D=='.endm':AT=G;continue
				if AT is not G:AS[AT][1].append(D)
				else:Bo.append(D)
			B7={};b=[];r={}
			def B8(expr):
				A=expr
				if A.isdigit():return R(A)
				return R(C9(A,{'__builtins__':{}},B7))
			def Bq(ln,env):
				A=ln
				for(B,C)in env.items():A=A.replace(rf"\{B}\()",g(C)).replace(f"\\{B}",g(C))
				return A
			def AU(lines,env):
				K=env;J=lines;E,M=0,C(J)
				while E<M:
					B=Bq(J[E],K);E+=1
					if not B:continue
					if B.startswith('.set '):O,P=B[5:].split(',',1);B7[O.strip()]=B8(P.strip());continue
					if B.startswith('.if '):
						H=[(B[4:].strip(),[])];F=0
						while E<M:
							D=J[E];E+=1
							if D.startswith('.if '):F+=1;H[-1][1].append(D)
							elif D=='.endif':
								if F==0:break
								F-=1;H[-1][1].append(D)
							elif F==0 and D.startswith('.elseif '):H.append((D[8:].strip(),[]))
							elif F==0 and D=='.else':H.append(('1',[]))
							else:H[-1][1].append(D)
						for(Q,I)in H:
							if B8(Bq(Q,K)):AU(I,K);break
						continue
					if B.startswith(BR):
						R=B8(B[6:].strip());I=[];F=0
						while E<M:
							D=J[E];E+=1
							if D.startswith(BR):F+=1;I.append(D)
							elif D=='.endr':
								if F==0:break
								F-=1;I.append(D)
							else:I.append(D)
						for X in A(R):AU(I,K)
						continue
					if w in B and not B.startswith('.'):
						N,S=B.split(w,1)
						if N.strip()and' 'not in N.strip():
							r[N.strip()]=C(b);B=S.strip()
							if not B:continue
					elif B.endswith(w):r[B[:-1].strip()]=C(b);continue
					L=B.split(G,1)
					if L[0]in AS:T,U=AS[L[0]];V=[ak.strip()for ak in(L[1].split(',')if C(L)>1 else[])];W={am:g(B7[al[1:]])if al.startswith('%')else al for(am,al)in v(T,V)};AU(U,W);continue
					b.append(B)
			AU(Bo,{})
			def AV(start_idx):
				B=[]
				for D in A(start_idx,C(b)):
					B.append(b[D])
					if b[D]==AJ:break
				return B
			K={};A2=(1<<54)+8
			for AW in A(790):K[m+16*AW]=AW+1;K[m+16*AW+8]=A_;K[A_+8*AW]=1<<30|1
			def Br(mode,s_op,bm,r0,Lc):
				E=mode;B=[an for an in A(Lc)if not bm>>(r0+an)%3&1]
				if not B:return B,h(AM)
				if C(B)==1:
					D=B[0]
					if D==0:return B,AV(r[f"e{E}_s0_p0"]+127*s_op)
					return B,AV(r[f"wheel_skip{D}_emit{E}_last"])
				D=B[0];H=[ap-ao for(ao,ap)in v(B,B[1:])];G=bin(bm).count('1');F=0 if G==0 else 3 if G==2 else 1 if H[-1]==1 else 2
				if D==0:I=128 if F==0 else 43 if F==3 else 86;return B,AV(r[f"e{E}_s0_p{F}"]+(I-C(B))*s_op)
				return B,AV(r[f"e{E}_s{D}_p{F}_k{C(B)}"])
			for(AX,A8,B9)in((8,8,6),(9,9,7)):
				for P in A(8):
					for e in A(3):
						for Z in(1,2,3,4,5,6,7,8,9,10,11,12,128):
							W,Bs=Br(AX,A8,P,e,Z);Bt=h(z)+Bs;CS=C(z);Q=C(Bt);A9={Ak:Q,Al:Q+1}
							for s in Bs:
								if s.startswith('tbnz ')and BS in s:A9[s.rsplit(',',1)[1].strip()]=Q+2
							CT=_b1_decode(Bt,A9);t=e&4294967295&~(7<<23)|P<<23;BA=t-(1<<32)if t>>31 else t;K[n+(BA<<3)&F]=A0;K[A0+(Z-1)*8]=4*CS;K[U+8]=A2+16;A3=_B1CPU(CT,K,{});B=A3.x;B[1],B[2],B[5],B[6],B[0],B[14]=8,m,0,e,Z-1,P;B[10],B[11],B[12],B[17],B[23],B[26],B[27]=n,48,3,U,1,16,1;AA={};A3.run(0,{Q,Q+1},{},AA);BB=(K[U+8]-(A2+16))//8;BC=[K[A2+16+8*ag]-(1<<38)-1 for ag in A(BB)];CU=4+(2 if W and W[0]>0 else 0)+(A8*(C(W)-1)+B9 if W else 0)+1
							if BC!=W or E(AA.values())!=CU:O=N;break
						if not O:break
					if not O:break
				if not O:break
			if O:
				CV,CW,_=AC(x,'chord_chunk')
				def CX(start_idx):
					B=[]
					for D in A(start_idx,C(b)):
						B.append(b[D])
						if b[D]==AD:break
					return B
				Bu=(0,0,128),(0,1,86),(0,1,85),(1,2,85),(0,3,43),(1,3,43),(2,3,42);Bv=[A0+8192*(ad+1)for ad in A(3)]
				for(AX,A8,B9)in((8,8,6),(9,9,7)):
					CY={ae:CX(r[f"c{AX}_s{ae[0]}_p{ae[1]}_k{ae[2]}"])for ae in Bu}
					for P in A(8):
						for e in(-5,-4,-3):
							for Z in(129,256,385,781):
								W=[ah for ah in A(Z)if not P>>(e+ah)%3&1];AY=(Z+127)//128;Bw=Z-128*(AY-1);CZ=(e+128*(AY-1))%3;_,Ca=Br(AX,A8,P,CZ,Bw);X=h(z);Bx=C(X);X.extend(CV);Cb=C(X);X.extend(AM);Cc=C(X);X.extend(Ca);By={}
								for Bz in Bu:By[Bz]=C(X);X.extend(CY[Bz])
								Q=C(X);A9={BP:0,BT:Bx+CW[BT],Ak:Q,Al:Q+1}
								for s in X:
									if s.startswith('tbnz ')and BS in s:A9[s.rsplit(',',1)[1].strip()]=Q+2
								Cd=_b1_decode(X,A9)
								for B_ in A(3):
									BD=Bv[B_];K[BD+(Z-1)*8]=4*(Cb if P==7 else Bx);K[BD+(Bw-1)*8]=4*Cc
									if P<7:AZ=[aj for aj in A(128)if not P>>(B_+aj)%3&1];Ce=AZ[0];Cf=AZ[-1]-AZ[-2];C0=bin(P).count('1');Cg=0 if C0==0 else 3 if C0==2 else 1 if Cf==1 else 2;K[BD+6248]=4*By[Ce,Cg,C(AZ)]
								for Ch in A(AY):C1=e+128*Ch;t=C1&4294967295&~(7<<23)|P<<23;BA=t-(1<<32)if t>>31 else t;K[n+(BA<<3)&F]=Bv[C1%3]
								K[U+8]=A2+16;A3=_B1CPU(Cd,K,{});B=A3.x;BE=0 if Z==781 else F;A3.N=i(BE>>63);B[1],B[2],B[5],B[6],B[0],B[14]=8,m,0,e&4294967295,Z-1,P;B[10],B[11],B[12],B[17],B[23],B[26],B[27]=n,48,3,U,BE,16,1;B[30],B[31]=0xcafebabe12345678,B0;AA={};Ci=A3.run(0,{Q,Q+1},{},AA);BB=(K[U+8]-(A2+16))//8;BC=[K[A2+16+8*ai]-(1<<38)-1 for ai in A(BB)];C2=2 if W and W[0]>0 else 0;Cj=8 if P==7 else 6+C2+(AY-1)*(6+C2+20)+(A8*(C(W)-1)+B9 if W else 0)+2;Ck=Q if BE==0 else Q+1
								if BC!=W or B[30]!=0xcafebabe12345678 or Ci!=Ck or E(AA.values())>Cj:O=N;break
							if not O:break
						if not O:break
					if not O:break
	except(Ac,S,Ad):O=N
	Aa,C3,BF={},{},{};AB=1<<44;C4={};u=O
	try:
		for(BG,C5)in((2046,19),(257,17),(1,15)):
			Aa[AB]=BG;Aa[AB+8]=C5
			for Cl in A(BG):A1=(1<<40)+CD.randrange(1,512);Aa[AB+16+8*Cl]=A1;BF[A1]=BF.get(A1,0)+C5&255
			BH=_B1CPU(Am,Aa,C3);BH.x[28]=AB;BH.run(A5[BQ],j(),{},C4)
			if BH.x[4]!=AB+16+8*BG:u=N
		u=u and C3==BF
	except(Ac,S,Ad):u=N
	C6=E(C4.values());C7=2304;C8=C7*(y+2)/(y//3*2)+3*(AF+3);u=i(u and An==8 and AF==11 and y==384 and BU<=.015 and BV<=M.C_LAT_APPLY and C6<=C8 and k(A[0]!=AK for A in Am));return l(ok=u,emit9_e2e_ok=O,emit9_flushes=AO,emit9_ops=C(AN),emit_direct_flush_ops=An,flush_ir_per_rec=I(BU,6),apply_page_ops=AF,apply_body_ops=y,apply_ir_per_rec=I(BV,6),c_lat_apply_budget=M.C_LAT_APPLY,pages_tested=3,recs_tested=C7,executed_ir=C6,charged_ir=I(C8,2))
K={}
def _b1_decode(prog,labels):
	F=labels;D=[]
	for B in prog:
		if B in K and not F:D.append(K[B]);continue
		A=G
		if B in K:A=K[B]
		else:
			for(E,I)in Aq:
				H=I.match(B)
				if H:
					J=h(H.groups());A=(E,)+o(R(A)if isinstance(A,g)and A.isdigit()else A for A in J)
					if E==Q:A=(E,A[1],A[2],R(A[3],16))if C(A)>3 else(E,A[1],R(A[2],16))
					K[B]=A;break
			else:A=AK,B;K[B]=A
		if A[0]in(U,W,Z,a,b):D.append(A[:-1]+(F[A[-1]],))
		else:D.append(A)
	return D
class _B1CPU:
	def __init__(A,P,m8,m16):A.P,A.x,A.m8,A.m16,A.Z,A.C,A.N=P,[0]*32,m8,m16,N,N,N
	def run(E,pc,stop,hook,tally):
		M=tally;G=pc;c,B,D,K=E.P,E.x,E.m8,E.m16
		while V:
			if G in stop:return G
			if G in hook:hook[G](E)
			A=c[G];C=A[0];H=G+1;M[G]=M.get(G,0)+1
			if C==e:B[A[1]]=D[B[A[3]]];B[A[2]]=D[B[A[3]]+8];B[A[3]]=B[A[3]]+A[4]&F
			elif C==f:B[A[1]]=D[B[A[3]]];B[A[2]]=D[B[A[3]]+8]
			elif C==m:D[B[A[3]]]=B[A[1]];D[B[A[3]]+8]=B[A[2]];B[A[3]]=B[A[3]]+A[4]&F
			elif C==q:D[B[A[3]]]=B[A[1]];D[B[A[3]]+8]=B[A[2]]
			elif C==r:B[A[1]]=D[B[A[2]]];B[A[2]]=B[A[2]]+A[3]&F
			elif C==s:B[A[1]]=D[B[A[2]]+B[A[3]]&F]
			elif C==u:D[B[A[2]]]=B[A[1]];B[A[2]]=B[A[2]]+A[3]&F
			elif C==A3:B[A[1]]=B[A[2]]+(B[A[3]]<<A[4])&F
			elif C==A5:B[A[1]]=B[A[2]]+B[A[3]]&F
			elif C==p:B[A[1]]=B[A[2]]+(A[3]<<12)&F
			elif C==x:B[A[1]]=B[A[2]]+A[3]&F
			elif C==A6:B[A[1]]=B[A[2]]-(B[A[3]]<<A[4])&F
			elif C==A7:B[A[1]]=B[A[2]]-A[3]&F
			elif C==A8:B[A[1]]=B[A[2]]-B[A[3]]&F
			elif C in(O,'cmp',P):
				N,R=(B[A[2]],B[A[3]])if C==O else(B[A[1]],A[2])if C==P else(B[A[1]],B[A[2]]);L=N-R&F;E.Z,E.C,E.N=L==0,N>=R,i(L>>63)
				if C==O:B[A[1]]=L
			elif C==A9:B[A[1]]=B[A[2]]if E.C else B[A[3]]
			elif C==Q:B[A[1]]=B[A[2]]&A[3]
			elif C=='lsr':B[A[1]]=B[A[2]]>>A[3]
			elif C=='ubfx':B[A[1]]=B[A[2]]>>A[3]&(1<<A[4])-1
			elif C==U:B[A[1]]=4*A[2]
			elif C=='blr':
				if B[A[1]]&3:raise S('misaligned blr')
				B[30]=4*H;H=B[A[1]]>>2
			elif C=='br':
				if B[A[1]]&3:raise S('misaligned br')
				H=B[A[1]]>>2
			elif C==W:H=A[1]
			elif C==AD:
				if B[30]==0:return
				if B[30]&3:raise S('misaligned ret')
				H=B[30]>>2
			elif C==Z:
				if not B[A[1]]>>A[2]&1:H=A[3]
			elif C==a:
				if B[A[1]]>>A[2]&1:H=A[3]
			elif C==b:
				d={'eq':E.Z,'hi':E.C and not E.Z,'lo':not E.C,'ls':not E.C or E.Z,'pl':not E.N}[A[1]]
				if d:H=A[2]
			elif C==AA:
				J=B[A[2]]
				if J&1:raise S('unaligned staddh')
				K[J]=K.get(J,0)+(B[A[1]]&65535)&65535
			elif C==AB:J=B[A[2]];K[J]=K.get(J,0)+(B[A[1]]&255)&255
			elif C==AL:B[A[1]]=D[B[31]+A[3]&F];B[A[2]]=D[B[31]+A[3]+8&F]
			elif C==AM:D[B[31]+A[3]&F]=B[A[1]];D[B[31]+A[3]+8&F]=B[A[2]]
			elif C==AN:B[A[1]]=D[B[A[3]]+A[4]];B[A[2]]=D[B[A[3]]+A[4]+8]
			elif C==AO:B[A[1]]=D[B[A[2]]+A[3]&F]
			elif C=='ldr':B[A[1]]=D[B[A[2]]]
			elif C==AP:B[A[1]]=D[B[A[2]]+(B[A[3]]<<3)&F]
			elif C==AQ:B[A[1]]=D[B[31]+A[2]&F]
			elif C==AR:D[B[31]+A[2]&F]=B[A[1]]
			elif C=='str':D[B[A[2]]]=B[A[1]]
			elif C=='ldadd':T=B[A[3]];X=D[T];D[T]=X+B[A[1]]&F;B[A[2]]=X
			elif C==AS:B[A[1]]=B[A[2]]+(B[A[3]]>>A[4])&F
			elif C==AT:B[A[1]]=B[A[2]]-(A[3]<<12)&F
			elif C==AU:B[A[1]]=A[2]&F
			elif C=='extr':B[A[1]]=(B[A[2]]<<64|B[A[3]])>>A[4]&F
			elif C=='sbfx':I=B[A[2]]>>A[3]&(1<<A[4])-1;B[A[1]]=I-(1<<A[4])&F if I>>A[4]-1 else I
			elif C=='mul':B[A[1]]=B[A[2]]*B[A[3]]&F
			elif C=='madd':B[A[1]]=B[A[2]]*B[A[3]]+B[A[4]]&F
			elif C=='msub':B[A[1]]=B[A[4]]-B[A[2]]*B[A[3]]&F
			elif C=='udiv':B[A[1]]=B[A[2]]//B[A[3]]if B[A[3]]else 0
			elif C==AV:B[A[1]]=B[A[2]]<<(B[A[3]]&63)&F
			elif C=='bfi_w':Y=(1<<A[4])-1<<A[3];B[A[1]]=(B[A[1]]&~Y|B[A[2]]<<A[3]&Y)&4294967295
			elif C==AW:I=B[A[3]]&4294967295;B[A[1]]=D[B[A[2]]+((I-(1<<32)if I>>31 else I)<<3)&F]
			elif C==AX:B[A[1]]=D[B[A[2]]+((B[A[3]]&4294967295)<<3)&F]
			else:raise S('instruction outside the executed set: '+g(A[1]))
			G=H
def CQ(path,seed=20261124,NA=2000,primes=((7,45,8,256),(61,94,24,8),(89,103,24,8)),nstreams=3):
	w=nstreams;U=NA;At,K,AG=As(path);o=_b1_decode(At,K);G={}
	for AH in AG:G[AH]=G.get(AH,0)+1
	W=K['b1_copy_end']-K['b1_copy'];Y=K['b1_apply_end']-K['b1_apply'];AI=G[AY]+4+2;A3=W/(W//2)+4./(W//2);AJ=[aa for(ab,aa)in v(AG,o)if ab==AZ and aa[0]in(p,x)and aa[1]==26];Z=E(A[3]<<12 if A[0]==p else A[3]for A in AJ);AL=16*(W//4-1);Au=C(AJ)==2 and 16<=Z and Z+X(8+AL,AL+16)<=16384;q=(Z-16)//8;AM,Av,Aw=G['b1_flush'],G[Am],G[AZ];AN=G[An]+3;Ax=A3+(AM+20)/q+(Y+2)/(Y//3*2)+(AN+20)/q;H=random.Random(seed);B,A4,r={},{},{};L=[1+ad*U//4 for ad in A(4)]+[U+1];AO={ac:E(1<<A for A in A(3)if H.random()<1/3)for ac in A(1,U+1)};D=l(rows=0,empty_rows=0,entries=0,lines=0,pages=0,term_pages=0,min_page=1<<30,max_page=0,switches=0);Ay,AP,a,b,A5=1<<40,1<<36,1<<44,1<<32,1<<34;Az=j();A6,AQ={},{};A7=[N]
	def AR(cpu):
		C=cpu.x[28];A=B[C]
		if 16+8*A>16384:raise S('page overflow')
		if not A7[0]:D[y]=min(D[y],A)
		D[Ao]=X(D[Ao],A);D[z]+=1;E=_B1CPU(o,B,A4);E.x[28]=C;E.run(K[An],Az,{},AQ)
		if E.x[4]!=C+16+8*A:raise S('receiver did not consume exactly its page')
	AS={K['b1_handoff']:AR,K['b1_term_handoff']:AR};d=V;s=0
	try:
		for(J,AT,A_,P)in primes:
			e=P*8;A8=32*e;B0={ae:H.randrange(J)for ae in A(1,U+1)};AU=[-(-(L[af+1]-L[af])//P)for af in A(4)];A9={}
			for M in A(1,U+1):A9.setdefault(B0[M],[]).append(M)
			t=Ay
			for f in A(J+1):
				for AV in A(8):
					for g in A(4):
						AW=AP+f*A8+(4*AV+g)*e
						if f==J:B[AW]=t;continue
						AX=[ai for ai in A9.get(f,())if L[g]<=ai<L[g+1]and not AO[ai]&AV]
						for Aa in A(P):B[AW+8*Aa]=t+8*bisect_left(AX,L[g]+Aa*AU[g])
						for M in AX:B[t]=M;t+=8
			AA=[H.randrange(J)for ag in A(6)];Ab=1<<30
			for(B1,B2)in c(AA):B[Ab+8*B1]=B2
			AB=_B1CPU(o,B,A4);F=AB.x;F[23]=J*A8;F[24]=AP;F[3]=b;F[29]=AT
			for Q in A(w):B[b+16*Q]=a+(Q<<14);B[b+16*Q+8]=a+(Q<<14)+16
			F[19]=b;F[28]=a;F[25]=a+16;F[26]=a+Z;Ae=1+(J<<8)
			for AC in A(A_):
				Af=H.randrange(1,40);AD=[H.randrange(-559,560-Af)]+[H.randrange(-510,511)for ah in A(5)];B3=E(A*B for(A,B)in v(AD,AA))%J;B4=AA[0];T=A5+AC*4096;h=A5-(1<<24)+32*AC
				for AE in A(Af):
					if AE==0 and J<16:m,O,u=0,3,P
					else:m,O,u=E(1<<A for A in A(3)if H.random()<1/3),H.randrange(4),H.randrange(1,P+1)
					Ag=AE==0 or H.random()<.08;Ah=Ae<<38;Ae+=1;B[T]=Ah;B[T+8]=(4*m+O)*e+8*u<<32|R(Ag)<<31|4*m*e;T+=16;s=X(s,4*m*e)
					if Ag:B[T]=H.randrange(w);T+=8;D[Ap]+=1
					f=(B3+AE*B4)%J;B5=L[O]+(u*AU[O]if u<P else L[O+1]-L[O]);AF=0
					for M in A9.get(f,()):
						if not AO[M]&m and M<B5 and M<L[O+1]:Ai=Ah+M;r[Ai]=r.get(Ai,0)+AT&255;AF+=1
					D['rows']+=1;D[A0]+=AF;D['empty_rows']+=AF==0
				Aj=lambda vs:E((B&65535)<<16*A for(A,B)in c(vs));B[h]=A5+AC*4096;B[h+8]=T;B[h+16]=Aj(AD[:4]);B[h+24]=Aj(AD[4:]);F[1]=h;F[2]=Ab;F[5]=J;F[6]=A8.bit_length()-1;F[7]=J*(6<<15);D['lines']+=1;AB.run(K[A1],j(),AS,A6)
			B[F[19]]=F[28];B[F[19]+8]=F[25];A7[0]=V
			for Q in A(w):
				n=b+16*Q
				if B[n+8]!=B[n]+16:F[28]=B[n];F[25]=B[n+8];AB.run(K[Am],j(),AS,A6);D[A2]+=1;B[n+8]=F[25]
			A7[0]=N
		d=A4==r
	except(Ac,S,Ad):d=N
	Ak=E(A6.values());Al=E(AQ.values());B6=D[z]-D[A2];Aq=D['rows']*AI+D[A0]*A3+B6*AM+D[A2]*Av+D[Ap]*Aw+D['lines']*(G[A1]+1);Ar=D[A0]*(Y+2)/(Y//3*2)+D[z]*AN;d=i(d and D[y]>=q and Ak<=Aq and Al<=Ar and G[AY]==20 and G[A1]+1<=64 and W==256 and Y==384 and Au and s>=1<<15 and k(A[0]!=AK for A in o));return l(ok=d,cells=C(r),stats=D,seg_counts=G,row_static=AI,ent_static=A3,flush_threshold=Z,max_start_offset=s,page_min_records=q,c_hit_lo_exact=Ax,switch_allow=.2,sender_ir=Ak,sender_charged=I(Aq,2),receiver_ir=Al,receiver_charged=I(Ar,2))
def verify_nonprim_hnf_witnesses():
	X=M.load_data()['pool_quad']
	if M._OHM_CACHE is G:M.is_ideal_primitive_oh([1,0,0,0,0,0],[0,1,0,0,0,0])
	O=M._OHM_CACHE
	def D(u,v):return[E(u[B]*v[C]*O[B][C][ac]for B in A(6)for C in A(6))for ac in A(6)]
	def F(v):return[[E(v[A]*O[A][ad][ae]for A in A(6))for ae in A(6)]for ad in A(6)]
	def P(rows):
		B=[af[:]for af in rows];I=1;E=0
		for F in A(6):
			while V:
				G=[ag for ag in A(E,C(B))if B[ag][F]!=0]
				if C(G)==1:
					D=G[0]
					if D!=E:B[E],B[D]=B[D],B[E]
					if B[E][F]<0:B[E]=[-ai for ai in B[E]]
					I*=B[E][F];E+=1;break
				D=min(G,key=lambda i:T(B[i][F]))
				if B[D][F]<0:B[D]=[-ah for ah in B[D]]
				K=B[D][F];L=B[D]
				for H in G:
					if H==D:continue
					M=B[H][F]//K;N=B[H]
					for J in A(F,6):N[J]-=M*L[J]
		return I
	Q=[-1,-1,-1,-1,1,-1];assert[E(Q[A]*M.ZK_NUM[A][aa]for A in A(6))//3 for aa in A(6)]==M.W_COEFFS;R=V;H=[]
	for(Y,I,J)in X['nonprim_witnesses']:
		B=P(F(I)+F(J));S=D(I,I);Z=D(J,J);a=D(Q,S);U=T(P(F([Z[ab]+(M.U_BN254+1)*a[ab]+2*S[ab]for ab in A(6)])));K=U%(B*B)==0;b=U//(B*B)if K else 0;L=0;N=b
		while N>0 and N&1==0:L+=1;N>>=1
		W=B==Y and K and L in(0,3,6);R&=W;H.append(l(N_I=B,divides_Ng_sq=K,v2_NJ=L,ok=W))
	return l(ok=R,n_checked=C(H),witnesses=H)
def _ht_load(raw):
	import json,lzma,base64 as Q;S,Y,K=lzma.decompress(Q.b85decode(raw)).partition(b'\n');F=json.loads(S);H,G,I,L,O=F['R'],F['G'],F['kw'],F['stw'],F['iuw'];M=C(I);B=0;J=H*M;D=H*G;T=[R.from_bytes(K[B+aa:B+aa+L],'little')for aa in A(0,J*L,L)];B+=J*L;P=[]
	for N in A(M):U=[K[B+ad*D:B+ad*D+D]for ad in A(I[N])];B+=I[N]*D;P.append([E(U[A][af]<<8*(I[N]-1-A)for A in A(I[N]))for af in A(D)])
	V=K[B:B+J*G];B+=J*G;D=J*G;W=[K[B+ab*D:B+ab*D+D]for ab in A(O)];X=[E(W[A][ac]<<8*(O-1-A)for A in A(O))+F['iu0']for ac in A(D)];return{d:F[d],'runs':[[[F['M'][ag],T[ae*M+ag],[[P[ag][ae*G+ah],V[(ag*H+ae)*G+ah],X[(ag*H+ae)*G+ah]]for ah in A(G)]]for ag in A(M)]for ae in A(H)]}
def verify_ht_pool(lat,extrap_factor=G):
	W=extrap_factor;H=lat;import json,lzma,base64;X=_ht_load((Path(__file__).resolve().parent/'ht_pool.b85').read_bytes());B,K=X['runs'],X[d];Z=H['prod_256orb_run'];s=H['prod_256orb_strata_ht'];t=B[Z-1];P=0
	for(u,((L,N,a),(w,x,b)))in c(v(t,s)):
		if L!=w or N!=x or C(a)!=C(b):P+=1;continue
		for((y,z,A0),(A1,A2,A3))in v(a,b):
			if K[u]+2*y!=A1 or z!=A2 or T(A0-A3*1e6)>J:P+=1
	D=8*C(B);e=[(B[0][aa][0],E(B[aa][1]for B in B),[(K[aa]+2*af,ag,ah/1e6)for ae in B for(af,ag,ah)in ae[aa][2]])for aa in A(C(K))];f=S=Y
	for(L,N,Q)in e:F=L*(D-J)/(N-J);f+=F*E(B*C/A for(A,B,C)in Q)/D;S+=F*E(B/A for(A,B,C)in Q)/D
	R=f/S;g=Y
	for(L,N,Q)in e:F=L*(D-J)/(N-J);h=[ac/ab*(ad-R)for(ab,ac,ad)in Q];U=E(h)/D;A4=E((A-U)**2 for A in h)/(D-1);i=F*F/(D-2.);g+=(F*F+i)*(A4/D)+i*U*U
	V=sqrt(g)/S;O=W if W is not G else H['prod_256orb_extrap_factor'];j=H['prod_256orb_gauss_wmean_ht_ir'];m=H['prod_256orb_gauss_se_strat_ir'];n=j+3.*m;A5=O*(R+3.*V);o=R+(O-J)*n+O*3.*V
	def A6(strata,Kd):
		I=strata;A=Kd;K=D=Y
		for(F,G,C)in I:B=F*(A-J)/(G-J);K+=B*E(B*C/A for(A,B,C)in C)/A;D+=B*E(B/A for(A,B,C)in C)/A
		L=K/D;M=Y
		for(F,G,C)in I:B=F*(A-J)/(G-J);N=[aq/ap*(ar-L)for(ap,aq,ar)in C];H=E(N)/A;M+=(B*B+B*B/(A-2.))*(E((A-H)**2 for A in N)/(A-1)/A)+B*B/(A-2.)*H*H
		return L,sqrt(M)/D
	p=[];q=[]
	for A7 in B:r,A8=A6([(aj,ak,[(K[ai]+2*am,an,ao/1e6)for(am,an,ao)in al])for(ai,(aj,ak,al))in c(A7)],8);p.append(I(r,3));q.append(I(O*(r+3.*A8),3))
	A9=P==0 and o<=M.C_LAT_LEAF and C(B)==16 and k(C(A)==32 for A in B);return l(ok=A9,runs=C(B),orbits=E(C(B[D][2])for B in B for D in A(C(K))),replay_run=Z,replay_errs=P,pooled_R=I(R,4),pooled_SE=I(V,4),extrap_factor=O,live_gauss_wmean=I(j,4),live_gauss_se=I(m,4),live_gauss_ub=I(n,4),pooled_bound_uniform_103=I(A5,4),pooled_bound=I(o,4),c_lat_leaf=M.C_LAT_LEAF,per_run_R=p,per_run_bound_uniform_103=q)
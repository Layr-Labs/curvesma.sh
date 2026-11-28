AC='emit8_step3'
AB='emit8_step2'
AA='emit8_step1'
A7='wheel_skip2_stub'
A6='wheel_disp4'
Cw='c_tail_z6'
Cv='hensel_div8_inner'
Cu='hensel_div8_outer'
Ct='reg_bingcd8_u_even'
Cs='reg_bingcd8_v_odd_lt'
Cr='reg_bingcd8_v_odd_ge'
Cq='reg_bingcd8_v_odd_chk'
Cp='reg_bingcd8_io'
Co='srt_odd_norm8_shift'
Cn='srt_odd_norm8_mid'
Cm='srt_odd_norm8_find'
Cl='cios8_csel_init'
Ck='cios8_csub_loop'
Cj='cios8_csub_init'
Ci='cios8_outer_post'
Ch='cios8_inner_red'
Cg='cios8_outer_mid'
Cf='cios8_inner_mul'
Ce='cios8_outer_pre'
Cd='cios8_setup'
Cc='//@seg'
Cb='get_rg'
Ca='get_sf_s0'
CZ='sv_sf_ref_350'
CY='sv_rg_ref_350'
CX='//@seg '
l='ex_g_bits'
e='trunc_g_bits'
d='ex_f_bits'
Y='trunc_f_bits'
s='small'
r='sieve_geo'
k='emit8_last'
Ba='reg_bingcd8_done_v'
BZ='reg_bingcd8_done_u'
BY='reg_bingcd8_v_even'
BX='reg_bingcd8_u_odd_lt'
BW='reg_bingcd8_u_odd_ge'
BV='inv64_odd'
BU='cios8_csel_loop'
BT='b.lo'
B0='x0'
A_='.endr'
Az='x16_nobfi_bnd2'
Ay='x16_nobfi_bnd1_lo'
Ax='x16_nobfi_bnd1_hi'
Aw=ValueError
h=zip
c=sorted
AY='reg_bingcd8_u_odd_chk'
AX='//'
AL='ok'
AK='x22'
AJ='//@end'
AI='x16_nobfi_pos'
A9=list
A8=all
X=tuple
q=any
n=min
g=dict
f=enumerate
b=True
a=1.
Z=float
U=','
T=False
P=.0
O=sum
K=pow
I=None
H=max
G=round
F=int
D=len
C=range
import math as AH,random,re,stor_bound as AW
from fractions import Fraction as Fr
from math import log,log2,exp,sqrt,ceil,floor
from pathlib import Path
o=AH.log(2.)
L3=log2(3.)
E=(1<<64)-1
j=.35
A=29.5
L=2.82
i={}
def get_ret2d_u8(M,pmin=24,B=61.5,odd=2,scale=j,bstep=a,dstep=.5,sl_step=.25):
	j=odd;J=sl_step;S=pmin;N=dstep;L=scale;A=bstep;k=S,B,j,L,A,N,J
	if k in i:return i[k]
	W=488;l,A8=-4e1,4e1;E=ceil(W/A);Q=F(G(l/N));p=F(G(A8/N));D=p-Q+1;t=.04;U=M.Retention(B,B,prime_min=1<<S,explicit_bits=14,step=A,max_bits=W,tilt=t,odd=j);K=M.Retention(Z(W),B,prime_min=1<<S,explicit_bits=14,step=A,max_bits=W+1,tilt=.005,odd=j)
	def u(tilt_v,is_sf):
		y=tilt_v;q=is_sf;d=[g()for ac in C(E+1)];Y=P;AD=F(floor(L*S+.5));AE=F(ceil(L*B+.5))+1;z=c(set([Z(S),Z(B)]+[(af-.5)/L for af in C(AD,AE+1)if S<(af-.5)/L<B]))
		for(e,r)in h(z[:-1],z[1:]):
			A0=F(floor(L*.5*(e+r)+.5));A1=H(4,F(ceil((r-e)/.05)));s=(r-e)/A1
			for AF in C(A1):
				f=e+(AF+.5)*s;R=G(f/A)
				if R>E:continue
				A2=R*A-A0/L;t=ceil(A2/N)if q else floor(A2/N);A3=exp(y*R*A)*s/f;Y+=A3;A4=d[R];A4[t]=A4.get(t,P)+A3
				for M in C(2,5):
					A5=M*R
					if A5>E:break
					b=exp(M*y*R*A-(M-1)*o*f)*s/(M*f)
					if b<=1e-25:break
					i=d[A5]
					if q:A6=-a if M%2==0 else a;Y+=A6*b;A7=H(Q,n(p,M*t));i[A7]=i.get(A7,P)+A6*b
					else:
						Y+=b
						for u in C(1,M+1):AG=AH.comb(M,u)*(a if u%2==1 else-a);A8=H(Q,n(p,floor((M*R*A-u*A0/L)/N)));i[A8]=i.get(A8,P)+AG*b
		A9=[(ad,X((B,ad*C)for(B,C)in d[ad].items()))for ad in C(1,E+1)if d[ad]];j=[[P]*D for ae in C(E+1)];j[0][-Q]=exp(-Y);AA=[(-Q,-Q)]+[(D,-1)]*E
		for T in C(A9[0][0],E+1):
			U=j[T];AI=a/T;J=D;K=-1
			for(v,AJ)in A9:
				if v>T:break
				w=j[T-v];V,W=AA[T-v]
				if V>W:continue
				for(I,x)in AJ:
					k=V if V+I>=0 else-I;l=W if W+I<D else D-1-I
					if k<=l:
						for m in C(k,l+1):U[m+I]+=x*w[m]
						if k+I<J:J=k+I
						if l+I>K:K=l+I
					if not q and V+I<0:
						AB=O(w[A]for A in C(V,n(W+1,-I)))
						if AB:
							U[0]+=x*AB;J=0
							if 0>K:K=0
					if W+I>=D:
						AC=O(w[A]for A in C(H(V,D-I),W+1))
						if AC:
							U[D-1]+=x*AC
							if D-1<J:J=D-1
							K=D-1
			if J<=K:
				for m in C(J,K+1):U[m]*=AI
				while J<=K and abs(U[J])<1e-26:J+=1
				while K>=J and abs(U[K])<1e-26:K-=1
				AA[T]=J,K
		return j,Y
	v,A9=u(t,b);w,AG=u(.005,T);AA=exp(A9-U.normalizers[r]);AB=U.co(s);x=K.co(s);AC=K.co('large_geo');y=[P]*(E+1);Y=E
	while Y>0 and x[Y]<1e-22:Y-=1
	for(z,R)in f(AC):
		if R>1e-22:
			for A0 in C(n(E-z,Y)+1):y[z+A0]+=R*x[A0]
	A1=[[P]*D for aa in C(E+1)];A2=[[P]*D for ab in C(E+1)]
	for d in C(E+1):
		V=P
		for(m,R)in f(v[d]):V+=R;A1[d][m]=V
		V=P
		for(m,R)in f(w[d]):V+=R;A2[d][m]=V
	q=.5/L+.016+.5*A;AD=U.denominator;A3=K.denominator;A4={};A5={};A6={};e={}
	def AE(n_int,sl_c):
		O=n_int;E=F(floor(sl_c/J))
		if E<0:return P
		R=O,E;B=A4.get(R)
		if B is not I:return B
		K=G(O/A);S=E*J-q;V=n(K,ceil((S-l)/A)+1);T=P
		for L in C(H(0,V+1)):
			U=AB[L]
			if U>1e-22:
				M=floor((S-L*A)/N)-Q
				if M>=0:T+=U*A1[K-L][M if M<D else D-1]
		B=H(P,n(a,T*AA/AD[K]));A4[R]=B;return B
	def AF(n_int,sl_c):
		R=n_int;S=F(ceil(sl_c/J));E=R,S;U=A5.get(E)
		if U is not I:return U,A6[E]
		B=G(R/A)
		if B not in e:e[B]=K.two(B,B,0,I,T)/A3[B]
		V=S*J+q;Z=n(B,ceil((V-l)/A)+1);W=P
		for M in C(H(0,Z+1)):
			X=y[M]
			if X>1e-22:
				O=floor((V-M*A)/N)-Q
				if O>=0:W+=X*A2[B-M][O if O<D else D-1]
		L=W/A3[B];Y=e[B]if e[B]<L else L;A5[E]=L;A6[E]=Y;return L,Y
	A7=g(get_rg=AE,get_sf_s0=AF,ret_1d_rg=U,ret_1d_sf=K,scale=L,bstep=A,dstep=N,sl_step=J,m_margin=q,sv_rg_sum_350=O(v[G(350/A)]),sv_rg_ref_350=U.co('sieve_sf')[G(350/A)],sv_sf_sum_350=O(w[G(350/A)]),sv_sf_ref_350=K.co(r)[G(350/A)]);i[k]=A7;return A7
_GSRT2_CACHE={}
def get_gsrt_2stage(M,B1=15,T1=86.0,B=61.5,odd=2,step=0.5):
	k=(B1,T1,round(B,3),odd,step)
	if k in _GSRT2_CACHE:return _GSRT2_CACHE[k]
	m_mar=0.5/j+0.016+0.5;max_b=496;nmax=int(round(max_b/step))
	ret=M.Retention(B,24.0,prime_min=(1<<B1),explicit_bits=min(B1,14),step=step,max_bits=max_b,tilt=0.04,odd=odd)
	c1,c2,c3=ret.co('small'),ret.co('sieve_geo'),ret.co('large_geo')
	conv12=[0.0]*(nmax+1)
	for ia in range(min(len(c1),nmax+1)):
		v1=c1[ia]
		if v1<1e-30:continue
		for ib in range(min(len(c2),nmax+1-ia)):conv12[ia+ib]+=v1*c2[ib]
	ret_fp=M.Retention(float(max_b-4),B,prime_min=(1<<24),explicit_bits=14,step=step,max_bits=max_b,tilt=0.005,odd=odd)
	cs_fp,c_sm_fp,c_lg_fp=ret_fp.co('sieve_geo'),ret_fp.co('small'),ret_fp.co('large_geo')
	ret_s1=M.Retention(float(max_b-4),24.0,prime_min=(1<<B1),explicit_bits=min(B1,14),step=step,max_bits=260,tilt=0.005,odd=odd)
	c1_fp,c_m2=ret_s1.co('small'),ret_s1.co('sieve_geo')
	r_lim=int(ceil((130.0+m_mar)/step))+2
	c_rem=[0.0]*(r_lim+1);c_comp=[0.0]*(r_lim+1)
	for ia in range(min(len(c_sm_fp),r_lim+1)):
		v1=c_sm_fp[ia]
		if v1>=1e-30:
			for ib in range(min(len(c_lg_fp),r_lim+1-ia)):c_rem[ia+ib]+=v1*c_lg_fp[ib]
	for ia in range(min(len(c_m2),r_lim+1)):
		v1=c_m2[ia]
		if v1>=1e-30:
			for ib in range(min(len(c_lg_fp),r_lim+1-ia)):c_comp[ia+ib]+=v1*c_lg_fp[ib]
	t1_b=int(floor((T1-m_mar)/step));t1_hi=int(ceil((T1+m_mar)/step))
	conv12_t1=[0.0]*(nmax+1)
	for jb in range(min(t1_b+1,len(c2))):
		v2=c2[jb]
		if v2>=1e-30:
			for ia in range(min(len(c1),nmax+1-jb)):conv12_t1[ia+jb]+=c1[ia]*v2
	c_rem_st2=[sum(c1_fp[ia]*c_comp[rv-ia]for ia in range(max(0,rv-t1_hi),min(rv+1,len(c1_fp))))for rv in range(r_lim+1)]
	c3_rev=c3[::-1];cs_rev=cs_fp[::-1];nc3=len(c3)-1;ncs=len(cs_fp)-1;sp=AH.sumprod
	rg_c={};fp_c={}
	def get_rg_rat(n_int,sl_c):
		sb=int(floor((sl_c-m_mar)/step))
		if sb<=0:return 1.0
		key=(n_int,sb)
		if key in rg_c:return rg_c[key]
		nb=min(nmax,int(round(n_int/step)));lim=min(sb,nb)+1;o=nc3-nb
		w3=c3_rev[o:o+lim]
		p1=sp(conv12[:lim],w3);p2=sp(conv12_t1[:lim],w3)
		res=min(1.0,max(0.0,p2/p1))if p1>1e-300 else 1.0
		rg_c[key]=res;return res
	def get_fp_st2(n_int,sl_c):
		sh=min(r_lim,int(ceil((sl_c+m_mar)/step)))
		if sh<=0:return 1.0
		key=(n_int,sh)
		if key in fp_c:return fp_c[key]
		nb=min(ncs,int(round(n_int/step)));lim=min(sh,nb)+1;o=ncs-nb
		w3=cs_rev[o:o+lim]
		sf_tot=sp(c_rem[:lim],w3);sf_st2=sp(c_rem_st2[:lim],w3)
		res=min(1.0,max(0.0,sf_st2/sf_tot))if sf_tot>1e-300 else 1.0
		fp_c[key]=res;return res
	out=dict(get_rg_rat=get_rg_rat,get_fp_st2=get_fp_st2,B1=B1,T1=T1,B=B,m_margin=m_mar)
	_GSRT2_CACHE[k]=out;return out
def verify_gsrt_2stage(M,seed=20261123):
	rng=random.Random(seed);B1=getattr(M,'GSRT2_B1',15);T1=getattr(M,'GSRT2_T1',86.0);B=M.OP['lg']
	ps1=M._primes_upto(1<<B1)[2:];ps2=[p for p in M._primes_upto(1<<17)[2:]if p>(1<<B1)]
	# 1. Exact algebraic partition check: stripping [0,B1] then (B1,24] removes identical smooth part as [0,24]
	def prod_tree(ps):
		a=list(ps)
		while len(a)>1:a=[a[i]*a[i+1]if i+1<len(a)else a[i]for i in range(0,len(a),2)]
		return a[0]
	P1=prod_tree(ps1);P2=prod_tree(ps2);P12=P1*P2;alg_errs=0
	def strip_all(n,P):
		while True:
			g=AH.gcd(n,P)
			if g==1:break
			n//=g
		return n
	for _ in range(100):
		n=rng.getrandbits(460)|1
		for _ in range(4):n*=rng.choice(ps1)
		for _ in range(2):n*=rng.choice(ps2)
		if strip_all(strip_all(n,P1),P2)!=strip_all(n,P12):alg_errs+=1
	# 2. Live Monte Carlo check of L1 tail (for fp_st2) against exact Retention model
	g2=get_gsrt_2stage(M,B1=B1,T1=T1,B=B,odd=2)
	nt=4000;hits30=0
	for _ in range(nt):
		s=0.0
		for p in ps1:
			if rng.random()<1.0/p:
				lp=log2(p);s+=lp
				while rng.random()<1.0/p:s+=lp
		if s>=30.0:hits30+=1
	ret_s1=M.Retention(492.0,24.0,prime_min=(1<<B1),explicit_bits=min(B1,14),step=0.5,max_bits=160,tilt=0.005,odd=2)
	c1=ret_s1.co('small');pc=[c1[i]*exp(-0.005*0.5*i)for i in range(len(c1))];zc=sum(pc);p_mod30=sum(pc[60:])/zc
	p_emp30=hits30/nt;se30=sqrt(max(1e-12,p_mod30*(1.0-p_mod30)/nt))
	mc_ok=p_emp30<=p_mod30+3.5*se30 and abs(p_emp30-p_mod30)<=4.0*se30
	rg467=g2['get_rg_rat'](467,M.OP['slack']);fp467=g2['get_fp_st2'](467,M.OP['slack'])
	return dict(ok=(alg_errs==0 and mc_ok and 0.98<rg467<1.0 and 0.2<fp467<0.6),alg_errs=alg_errs,mc_ok=mc_ok,p_emp30=round(p_emp30,5),p_mod30=round(p_mod30,5),rg_rat_467=round(rg467,5),fp_st2_467=round(fp467,5))
J=I
def get_exact_23_classes():
	global J
	if J is not I:return J
	R=6;B={(0,0):Fr(1,2)};S=Fr(1,2)
	for E in C(1,R):T=Fr(1,3**E);B[2,E]=T;S-=T
	B[2,R]=S;assert O(B.values())==1
	def L(d1,d2):
		A={}
		for((C,D),E)in d1.items():
			for((F,G),H)in d2.items():B=C+F,D+G;A[B]=A.get(B,Fr(0))+E*H
		return A
	k=L(L(B,B),B);l={(0,0):Fr(27,28),(6,0):Fr(1,28)};U=L(k,l);assert O(U.values())==1;M={}
	for((G,H),K)in U.items():E=n(G,12),n(H,12);M[E]=M.get(E,Fr(0))+K
	m=(Fr(7,9),0,0),(Fr(1,9),6,0),(Fr(7,72),0,6),(Fr(1,72),0,9);o,p=Fr(2,3),Fr(17,24);q,r=Fr(45,14),Fr(9,4);V=Z(o)+Z(q)*L3;W=Z(p)+Z(r)*L3;A=[];Y=P;a=Fr(0);N=P;Q=P
	for(s,b,d)in m:
		for((G,H),t)in c(M.items()):
			K=s*t;a+=K;F=Z(K)
			if b==0 and d==0 and G==0 and H==0 or F<1e-05:Y+=F
			else:e=b+G*L3;h=d+H*L3;N+=F*e;Q+=F*h;A.append((F,e,h))
	A.insert(0,(Y,P,P));i=c({ab[1]for ab in A});j=c({ac[2]for ac in A});u={y:x for(x,y)in f(i)};v={aa:z for(z,aa)in f(j)};w=X((A,u[B],v[C])for(A,B,C)in A);assert a==1;assert N<=V+1e-12 and Q<=W+1e-12;J=g(classes=X(A),u_vf=X(i),u_vg=X(j),sparse_fg=w,ex_f_bits=V,ex_g_bits=W,trunc_f_bits=N,trunc_g_bits=Q,n_classes=D(A));return J
def ET(M,seed=20261123):
	E=get_exact_23_classes();L=(1,473,344,58,71,434),(1,642,473,17,681,58),(1,344,642,434,199,17)
	if M._OHM_CACHE is I:M.is_ideal_primitive_oh([1,0,0,0,0,0],[0,1,0,0,0,0])
	n=M._OHM_CACHE;R=0
	for F in L:
		for N in C(6):
			for Z in C(6):
				o=F[N]*F[Z]%729;p=O(n[N][Z][A]*F[A]for A in C(6))%729
				if o!=p:R+=1
	a=D({X(A%3 for A in al)for al in L})==3;from fractions import Fraction as S;r=[[S(ae,3)for ae in aa]for aa in M.ZK_NUM];A=[r[ab][:]+[S(1 if ab==am else 0)for am in C(6)]for ab in C(6)]
	for B in C(6):
		b=next(C for C in C(B,6)if A[C][B]!=0);A[B],A[b]=A[b],A[B];s=A[B][B];A[B]=[af/s for af in A[B]]
		for P in C(6):
			if P!=B and A[P][B]:t=A[P][B];A[P]=[A[P][an]-t*A[B][an]for an in C(12)]
	u=[ac[6:]for ac in A]
	def c(poly):A=[O(S(poly[A])*u[A][ao]for A in C(6))for ao in C(6)];assert A8(A.denominator==1 for A in A);return[ap.numerator for ap in A]
	v={ad:c(M.F_REC_COEFFS[ad])for ad in(8,6,4,2,0)};w=c(M.G_LEAD_COEFFS);x=[1,0,0,0,0,0];T=0
	for F in L:
		y=[O(v[ag][A]*F[A]for A in C(6))%729 for ag in(8,6,4,2,0)];z=[O(w[A]*F[A]for A in C(6))%729,O(x[A]*F[A]for A in C(6))%729]
		if y!=[1,374,615,650,286]or z!=[68,1]:T+=1
	f=random.Random(seed);h=lambda x:6 if x%729==0 else next(A for A in C(6)if x//3**A%3!=0);Q=0;U=0
	while U<1500:
		A0=[f.randrange(729)for ah in C(6)];A1=[f.randrange(729)for ai in C(6)];i=[O(A0[A]*aj[A]for A in C(6))%729 for aj in L];j=[O(A1[A]*ak[A]for A in C(6))%729 for ak in L]
		if q(i[A]%3==0 and j[A]%3==0 for A in C(3)):continue
		U+=1
		for N in C(3):
			V,W=i[N],j[N];H,J=V*V%729,W*W%729;A2=(K(H,4,729)+374*K(H,3,729)*J+615*H*H*J*J+650*H*K(J,3,729)+286*K(J,4,729))%729;A3=(68*H+J)%729;k,m=h(A2),h(A3)
			if V%3==0 or W%3==0:
				if k!=0 or m!=0:Q+=1
			elif k!=2 or m<1:Q+=1
	A4=R==0 and a and T==0 and Q==0 and E[Y]<=E[d]+1e-12 and E[e]<=E[l]+1e-12;return g(ok=A4,hom_errs=R,distinct_3_primes_ok=a,coeff_mod729_errs=T,law_errs=Q,tested=U,n_classes=E['n_classes'],ex_f_bits=G(E[d],5),ex_g_bits=G(E[l],5),trunc_f_bits=G(E[Y],5),trunc_g_bits=G(E[e],5))
def split_ntt3_parts(M,R,lf,lg):
	C=Z(1<<46);T=ceil(R/C);E={};J=P;F=P
	for(U,K,D,L,N)in(('f',M.KH_F_CAP,lf,448,7),('g',M.KH_G_CAP,lg,512,8)):V=M.pbits(0,D,K);A=ceil(V/64.);W=1.25506*2.**(.5*D)/(.5*D*o)+K*1.25506*2.**D/(D*o);O=ceil(log2(W));G=L//24;B=(N+G)*C;X=ceil(log2(B));Y=3.*(M.ntt3_F(A)+M.ntt3_F(B));a=X*(M.ntt3_F(B)+M.mul2_F(B)+16.*B);b=G*C*(2.*(2.*L*(12.*N+16.))+2e4);c=G*C*1e2;d=O*(M.ntt3_F(A)+16.*A);e=R*6e5;f=R*35.*2.**19;Q=T*O*(Y+a+b+c)+d+e+f;E[U]=Q;J+=Q;I=H(2,ceil(log2(A)));g=Z(1<<I);S=3.*Z(1<<I-2);h=S if I>=3 and S>=A else g;F=H(F,28.*h+24.*A+32.*B+8.*C)
	E['total_ir']=J;E['ram_log2']=log2(F);return E
def b1_delta_save_per_cell(M,raw_mult,p_min=24.,p_hi=A):
	B=p_min;J=[AH.comb(3,m)*(2./9.)**m*(2./3.)**(3-m)/(512./729.)for m in C(4)];E=P;H=F(G((p_hi-B)*4));A=(p_hi-B)/H
	for K in C(H):
		D=B+(K+.5)*A;N=2.**(D+.5*A);Q=O(J[A]*exp(-2.**32*(2./3.)**A/N)for A in C(4));I=L+24.*Q
		if I<M.C_HIT_LO:E+=(M.C_HIT_LO-I)*log((D+.5*A)/(D-.5*A))
	return M._WKF*M.B1_PREF_COV*raw_mult*E
def verify_v32_all(M,kernels_path,seed=20261123):
	B1=kernels_path;N=random.Random(seed);Cx=ET(M,seed=seed);Cy=M._primes_upto(200)
	def Bb(arr):
		A=arr
		if D(A)==1:return A[0],I,I
		B=D(A)//2;C,E=Bb(A[:B]),Bb(A[B:]);return C[0]*E[0],C,E
	EW=Bb(Cy);Bc=[]
	for _ in C(64):
		B2=1
		for _ in C(N.randint(3,9)):B2*=N.choice(Cy)
		Bc.append(B2)
	def Bd(node,frags):
		E=frags
		if not E:return[]
		F,G,P=node
		if G is I:
			N=[]
			for(B,A)in E:
				H=0;J=A
				while J%F==0:J//=F;H+=1
				assert J==1 and H>=1;N.append((B,F,H))
			return N
		Q=G[0];K,L=[],[]
		for(B,A)in E:
			C=AH.gcd(A,Q%A)
			if C==1:L.append((B,A))
			elif C==A:K.append((B,A))
			else:
				O=C;D=A//C
				while b:
					M=AH.gcd(D,C)
					if M==1:break
					O*=M;D//=M
				K.append((B,O))
				if D>1:L.append((B,D))
		return Bd(G,K)+Bd(P,L)
	EX=Bd(EW,A9(f(Bc)));Cz={}
	for(y,z,Be)in EX:Cz.setdefault(y,[]).append((z,Be))
	Bf=0
	for(y,B2)in f(Bc):
		C_=1
		for(z,Be)in Cz[y]:C_*=z**Be
		if C_!=B2:Bf+=1
	c=M.R_BN254;D0=K(M.P_BN254,10,c);EY=(1<<317)//c;EZ=D0*(1<<256)%c;Ea=-K(c,-1,1<<64)&E
	def Eb(a,b):
		A=0
		for B in C(4):D=a>>64*B&E;A+=D*b;F=A*Ea&E;A=A+F*c>>64
		if A>=c:A-=c
		return A
	def D1(A):
		D=A<0;C=-A if D else A;F=(C>>256<<3|C>>192>>61)&E;G=F*EY>>64;B=C-G*c
		if B>=c:B-=c
		if B>=c:B-=c
		assert 0<=B<c;return c-B if D and B else B
	Bg=0
	for _ in C(200):
		D2=N.randrange(-(1<<313),1<<313);D3=N.randrange(-(1<<313),1<<313);B3=N.randrange(-(1<<313),1<<313);Ec=(D1(D2-B3)+Eb(D1(D3+B3),EZ))%c;Ed=(D2-B3+D0*(D3+B3))%c
		if Ec!=Ed:Bg+=1
	D4=323;import stor_bound as X;Ee=X.main();Ft=X.SB(X.QB_MAX);Fu=F(X.cnt_a_up(X.QA_MAX+X.DELTA))+1;B4=A9(C(X.K_MIN_STRIDE,Ee['kappa_max']+1));D5=X.F;Ef=X.dn(D5(314159265,100000000)**3/6);Eg=X.dn(1/X.sqrt_up(X.det(X.G)));Eh=X.sqrt_dn(X.QA_MAX*X.G_SC*D5(56880,100)**2)-X.E6.h;Ei=F(Ef*Eh**6*Eg/6)-1;Bh=O(n(a,(1<<A)/Z(Ei))for A in B4);AM={};AZ=I;Ej=AI,Ax,Ay,Az;D6={aa+ab for aa in Ej for ab in('_pre','','_tail')}
	for D7 in AW.ko(B1):
		Bi=D7.strip()
		if Bi.startswith(CX):AZ=Bi.split()[1];AM[AZ]=[];continue
		if Bi.startswith(AJ):AZ=I;continue
		if AZ in D6:
			AN=D7.split(AX,1)[0].strip()
			if AN:AM[AZ].append(AN)
	def AA(sn,ents,bp,wa,r24=0,r25=0,r26=0):
		M=ents;J=AM.get(sn+'_pre',[])+AM.get(sn,[])+AM.get(sn+'_tail',[]);K,H=[],0
		while H<D(J):
			if J[H].startswith('.rept '):
				b=F(J[H].split()[1]);H+=1;R=[]
				while H<D(J)and J[H]!=A_:R.append(J[H]);H+=1
				K.extend(R*b)
			elif J[H]!=A_:K.append(J[H])
			H+=1
		I={f"x{aL}":0 for aL in C(32)};I.update(x10=bp,x19=wa&255,x22=0,x23=D(M)*8,x24=r24,x25=r25,x26=r26);N=S=T;Q={};L=V=0;B=lambda s:I['x'+s[1:]]&(4294967295 if s[0]=='w'else E)
		def O(s,v):I['x'+s[1:]]=v&(4294967295 if s[0]=='w'else E)
		def W(tg):return next(A for A in(C(L+1,D(K))if tg.endswith('f')else C(D(K)))if K[A]==(tg[:-1]if tg.endswith('f')else tg)+':')
		while L<D(K):
			P=K[L]
			if P.endswith(':'):L+=1;continue
			V+=1;A=[aM for aM in P.replace('[',' ').replace(']',' ').replace(U,' ').split()if aM];G=A[0]
			if G=='ldp'and A[1:5]==[B0,'x1',AK,'#16']:X=I[AK]//8;I[B0],I['x1']=M[X],M[X+1];I[AK]+=16
			elif G=='ldr'and A[1:4]==[B0,AK,'#8']:I[B0]=M[I[AK]//8];I[AK]+=8
			elif G=='cmp'and D(A)==3:Y,Z=B(A[1]),B(A[2]);N=Y>=Z;S=Y==Z
			elif G=='lsr'and D(A)==4:
				c=F(A[3][1:])if A[3].startswith('#')else B(A[3])&63
				if A[3]=='#48':raise Aw(P)
				O(A[1],B(A[2])>>c)
			elif G=='bic'and D(A)==4:O(A[1],B(A[2])&~F(A[3][1:]))
			elif G=='sub':O(A[1],B(A[2])-B(A[3]))
			elif G=='add':d=B(A[3])<<F(A[5][1:])if D(A)==6 and A[4]=='lsl'else B(A[3]);O(A[1],B(A[2])+d)
			elif G=='staddb':a=B(A[2]);Q[a]=Q.get(a,0)+(B(A[1])&255)&255
			elif G=='tbz':
				if B(A[1])>>F(A[2][1:])&1==0:L=W(A[3]);continue
			elif G in(BT,'b.hi','b.hs'):
				if not N if G==BT else N if G=='b.hs'else N and not S:L=W(A[1]);continue
			elif G=='ret':break
			else:raise Aw(P)
			L+=1
		return Q,V
	Ek=A8('#48'not in'\n'.join(AM.get(A,[]))for A in D6);A0=0 if Ek and B4==A9(C(41,49))and Bh<=M.B2_PREFIX_MULT and M.C_LIST_B2>=29.7+2.*(Bh-a)*3. else 1;Bj=0
	for k in B4:
		AO=1<<k
		for A1 in(8,16,32):
			for D8 in(1,-1):
				p=N.randrange(1<<40,(1<<63)-A1*AO-16);AP=[]
				for _ in C(16):D9=N.randrange(A1);d=N.randrange(1,AO);l=D9 if D8>0 else A1-1-D9;AP.append((l,d,l<<k|d))
				DA=k==41 and A1==8 and D8==1
				if DA:AP[0],AP[1],AP[2]=(1,1,1<<41|1),(0,(1<<41)-1,(1<<41)-1),(1,0,1<<41)
				Aa,Ab=(0,1)if DA else(N.randrange(A1),N.randrange(A1));A2=N.randrange(A1);Ac=N.randrange(A2,A1);B5=[ax[2]for ax in AP];(h,Fv),(El,Fw),(Em,Fx),(En,Fy)=AA(AI,B5,p,13),AA(Ax,B5,p,13,r25=(Aa+1<<k)-1),AA(Ay,B5,p,13,r25=Ab<<k),AA(Az,B5,p,13,r24=k,r25=A2,r26=Ac-A2);Bk,Bl,Bm,Bn={},{},{},{}
				for(l,d,DB)in AP:
					Bj+=1;i=p+l*AO+d
					if DB>=1<<53 or p+DB!=i or i>=1<<64:A0+=1
					Bk[i]=Bk.get(i,0)+13&255
					if l<=Aa:Bl[i]=Bl.get(i,0)+13&255
					if l>=Ab:Bm[i]=Bm.get(i,0)+13&255
					if A2<=l<=Ac:Bn[i]=Bn.get(i,0)+13&255
				if h!=Bk or El!=Bl or Em!=Bm or En!=Bn:A0+=1
	for k in C(41,49):
		AO=1<<k;p=1<<50;Aa,Ab,A2,Ac=3,2,1,5
		for DC in(0,1,7,15,16,17,19,31,33):
			DD=[(N.randrange(8),N.randrange(1,AO))for aq in C(DC)];Eo=[ar<<k|at for(ar,at)in DD]
			for(Ep,Eq,Er,Es)in((AI,{},lambda j:b,42./16.),(Ax,g(r25=(Aa+1<<k)-1),lambda j:j<=Aa,74./16.),(Ay,g(r25=Ab<<k),lambda j:j>=Ab,74./16.),(Az,g(r24=k,r25=A2,r26=Ac-A2),lambda j:A2<=j<=Ac,106./16.)):
				Et,Eu=AA(Ep,Eo,p,13,**Eq);Bo={}
				for(l,d)in DD:
					if Er(l):i=p+l*AO+d;Bo[i]=Bo.get(i,0)+13&255
				if Et!=Bo or Eu>Es*DC+11.55:A0+=1
	Ad=[];Ae={};Bp=T;Af=I
	for AN in AW.ko(B1):
		S=AN.strip()
		if S in('//@seg x16_bfi_copy_page','//@seg x16_bfi_copy'):Bp=b;continue
		if S==AJ:Bp=T;continue
		if not Bp or not S or S.startswith(AX):continue
		if S.startswith('.rept'):Ev=F(S.split()[1]);Af=[];continue
		if S.startswith(A_):
			for _ in C(Ev):Ad.extend(Af)
			Af=I;continue
		if S.endswith(':'):Ae[S[:-1]]=D(Ad);continue
		if Af is not I:Af.append(S)
		else:Ad.append(S)
	for r in(0,1,2,3,15,16,127,128,129,255,256,2045,2046):
		Bq=1<<40;DE=16;B6={Bq+8*al:3405643776+al for al in C(r)};AQ={};J={11:r,15:Bq,25:DE};AB=0;Br=0;DF=T
		while AB<D(Ad):
			R=Ad[AB];Br+=1;AB+=1
			if R.startswith('cbz x11,'):
				if J[11]==0:AB=Ae[R.split(U)[1].strip()]
			elif R=='add x16, x15, x11, lsl #3':J[16]=J[15]+(J[11]<<3)
			elif R.startswith('tbz x11, #0,'):
				if J[11]&1==0:AB=Ae[R.split(U)[2].strip()]
			elif R=='ldr x0, [x15], #8':J[0]=B6[J[15]];J[15]+=8
			elif R=='str x0, [x25], #8':AQ[J[25]]=J[0];J[25]+=8
			elif R=='ubfx x12, x11, #1, #6':J[12]=J[11]>>1&63
			elif R=='adr x13, copy_end':J[13]=Ae['copy_end']*4
			elif R=='sub x13, x13, x12, lsl #3':J[13]-=J[12]<<3
			elif R=='br x13':AB=J[13]//4
			elif R=='ldp x0, x1, [x15], #16':J[0]=B6[J[15]];J[1]=B6[J[15]+8];J[15]+=16
			elif R=='stp x0, x1, [x25], #16':AQ[J[25]]=J[0];AQ[J[25]+8]=J[1];J[25]+=16
			elif R=='cmp x15, x16':DF=J[15]<J[16]
			elif R.startswith('b.lo '):
				if DF:AB=Ae[R.split()[1]]
			elif R=='ret':break
			else:raise Aw(R)
		Ew={DE+8*am:B6[Bq+8*am]for am in C(r)}
		if AQ!=Ew or AQ and H(AQ)+8>16384 or Br>r*(13e1/128.)+11.+1e-09:A0+=1
		if r==2046 and(Br+59.)/2046.>1.05:A0+=1
	DG=[(N.randrange(8),N.randrange(1,1<<41))for ac in C(33)];DH=[ad<<41|ae for(ad,ae)in DG];Ex,Ey=AA(AI,DH[:19],1<<50,13);Ez,E_=AA(AI,DH[19:],1<<50,13);Bs=g(Ex)
	for(Ag,F0)in Ez.items():Bs[Ag]=Bs.get(Ag,0)+F0&255
	Bt={}
	for(l,d)in DG:Ag=(1<<50)+(l<<41)+d;Bt[Ag]=Bt.get(Ag,0)+13&255
	if Bs!=Bt or Ey+E_>33*(42./16.)+23.1:A0+=1
	Bu=0;B7=[10,250,(1<<32)+500,(1<<32)+900,(1<<33)+1200,(1<<33)+1205,(1<<33)+1299,(1<<33)+1400,(1<<33)+1555];B8=[];Bv=B7[0];B9=[0];Bw=B7[0]
	for d in B7[1:]:
		L=d-Bw
		if L>=1<<32:B8.append((Bv,B9));Bv=d;B9=[0];Bw=d
		else:B9.append(L);Bw=d
	B8.append((Bv,B9));Bx=1<<50;By=17;Ah={};Bz={}
	for d in B7:Bz[Bx+d]=Bz.get(Bx+d,0)+By&255
	for(F1,Ai)in B8:
		A3=Bx+F1
		for BA in C(0,D(Ai),2):
			DI=Ai[BA]|(Ai[BA+1]<<32 if BA+1<D(Ai)else 0);A3=A3+(DI&4294967295)&E;Ah[A3]=Ah.get(A3,0)+By&255
			if BA+1<D(Ai):A3=A3+(DI>>32)&E;Ah[A3]=Ah.get(A3,0)+By&255
	if Ah!=Bz or D(B8)!=3:Bu+=1
	DJ=M.load_data()['t4_geom'];DK=DJ['GA'][0][0]/3.;F2=n(abs(Z(M.U_BN254+1)*complex(A,B)+2.)for(A,B)in DJ['EMB_W']);DL=572.**2*M.QA_MAX/DK;F3=782.76**2*M.QB_MAX/DK;F4=a/(DL*DL);DM=F3/(F2*F4);DN=H(6.*log2(a+DM),-6.*log2(a-DM));Aj=M.OP['slack'];DO=478.;B_=F(floor(j*(DO-Aj)+.5));DP=128-B_;F5=H(F(floor(j*(24.+A*.01)+.5))/(24.+A*.01)for A in C(3701));DQ=F(ceil(128.+j*Aj+.5+(F5-j)*DO+.02));m=get_ret2d_u8(M,pmin=M.OP['p_min'],B=M.OP['lg'],odd=2);DR=abs(m['sv_rg_sum_350']-m[CY])/H(m[CY],1e-30);DS=abs(m['sv_sf_sum_350']-m[CZ])/H(m[CZ],1e-30);DT,F6=m[Ca](468,Aj);F7=m[Cb](468,Aj);F8=m['ret_1d_sf'].survivor_upper(468.,Aj,ex=P);BB,BC=m['bstep'],m['dstep'];BD=0
	for DU in[24.51]+[24.+an*.01 for an in C(3701)]:
		C0=F(floor(j*DU+.5));C1=G(DU/BB);DV=C1*BB-C0/j
		if C1*BB-ceil(DV/BC)*BC>C0/j+1e-12:BD+=1
		if C1*BB-floor(DV/BC)*BC<C0/j-1e-12:BD+=1
	DW=M.rho(33.5);DX=M.rho(4e1);DY=1e-60<DW<1e-55 and 1e-75<DX<1e-70 and D(M.rho_tab_1024())>=41984;DZ=get_ret2d_u8(M,pmin=14,B=3e1,odd=2);F9=M._primes_upto(4096)
	def FA(n_v):
		A=n_v;D={}
		for B in F9:
			if B*B>A:break
			if A%B==0:
				F=0
				while A%B==0:A//=B;F+=1
				D[B]=F
		E=[A]if A>1 else[]
		while E:
			C=E.pop()
			if C==1:continue
			if M._isp(C):D[C]=D.get(C,0)+1
			else:G=M._brent(C);E.append(G);E.append(C//G)
		return D
	C2=4000;C3=22.;BE=Da=Db=0;FB=F(2.**49.5);FC=1<<50;FD=random.Random(seed^23100)
	for _ in C(C2):
		while b:
			C4=FD.randrange(FB,FC)
			if C4%6 in(1,5):break
		Dc=FA(C4);FE=O(F(floor(j*log2(A)+.5))for A in Dc if 1<<14<=A<1<<30);FF=F(floor(j*(log2(C4)-C3)+.5));Dd=FE>=FF
		if Dd:Db+=1
		if A8(A<1<<30 for A in Dc):
			BE+=1
			if Dd:Da+=1
	BF=Da/H(1,BE);Fz=sqrt(BF*(a-BF)/H(1,BE));BG=Db/Z(C2);F_=sqrt(BG*(a-BG)/Z(C2));De=DZ[Cb](50,C3);Df,_=DZ[Ca](49,C3);Dg=BE>=1400 and BF>=De and BG<=Df;BH=0;C5=0;Dh=[(af*7+13)%729<<3|af*3+5&7 for af in C(729)for ag in C(8)]
	for BI in C(729):
		FG=Dh[BI<<3]
		for FH in C(8):
			if Dh[BI<<3|FH]!=FG:C5+=1
	Di={}
	for Ak in C(8):
		for BJ in(-4194304,-4194303,-2100000,-1000001,-3,-2,-1,0,1,2,3,1234567,2100000,4194303):
			FI=BJ&4294967295;BI=(Ak*89+abs(BJ))%729;FJ=BI<<3|Ak;C6=FI&~(7<<23)|(FJ&7)<<23;FK=C6-(1<<32)if C6&1<<31 else C6;y=(1<<26)+FK;FL=(Ak|Ak<<3)>>BJ%3&7
			if not 0<=y<1<<27:BH+=1;continue
			C7=y-(1<<26);C8=C7>>23&7;Dj=C7&8388607;Dk=Dj-(1<<23)if C7<0 else Dj;Dl=(C8|C8<<3)>>Dk%3&7
			if C8!=Ak or Dk!=BJ or Dl!=FL or y in Di:BH+=1
			Di[y]=Dl
	Dm=DR<.001 and DS<.001 and DT>=F8 and .9<F7<a and F6<=DT and BD==0 and Dg and DY;Dn=DN<=.02 and B_<=128 and DP>=0 and DQ<=255 and Dm and BH==0 and C5==0;C9=0;r=128;FM=[128-N.randint(90,127)for ah in C(r)];Al=[FM[:]for ai in C(8)];Do=set()
	for s in C(8):
		for AR in C(r):
			if N.random()<.03:Al[s][AR]=Al[s][AR]+N.randint(100,130)&255
			if Al[s][AR]&128:Do.add((s,AR))
	Dp=set()
	for CA in C(0,r,8):
		Dq=0
		for s in C(8):
			Dr=0
			for Am in C(8):Dr|=Al[s][CA+Am]<<8*Am
			Dq|=Dr
		if Dq&0x8080808080808080:
			for s in C(8):
				for Am in C(8):
					if Al[s][CA+Am]&128:Dp.add((s,CA+Am))
	if Dp!=Do:C9+=1
	FN=re.compile('^\\s+[a-z][a-z0-9.]*\\b');A={};t=I;An=1
	for Ds in AW.ko(B1):
		R=Ds.strip()
		if R.startswith(Cc):t=R.split()[1];A.setdefault(t,0);An=1;continue
		if R.startswith(AJ):t=I;An=1;continue
		AC=Ds.split(AX,1)[0].rstrip()
		if t is I or not AC.strip()or AC.endswith(':'):continue
		if AC.strip().startswith('.rept'):An=F(AC.strip().split()[1]);continue
		if AC.strip().startswith(A_):An=1;continue
		if AC.strip().startswith('.'):continue
		if FN.match(AC):A[t]+=An
	Ao=13e1/128.+92./2046.;Dt=A.get('delta_b1_copy')==35 and A.get('delta_b1_apply')==291 and A.get('delta_b1_apply_tail')==8 and A.get('x16_bfi_copy_page')==9 and A.get('x16_bfi_copy')==131 and 13e1/128.+7e1/2046.<=1.05 and M.C_SB_ENTER_X4>=42.6 and A.get(AI)==42 and A.get(Ax)==74 and A.get(Ay)==74 and A.get(Az)==106 and 42./16.+Ao<=M.C_B2_INT_X4 and 74./16.+Ao<=M.C_B2_BND1_HIT and 42./16.+Ao<=M.C_B2_BND1_MISS and 106./16.+Ao<=M.C_B2_BND_HIT and 74./16.+Ao<=M.C_B2_BND_MISS and A.get('init8_u8')==75 and A.get('scan8_u8')==199 and 271./1024.<=M.C_SCAN_G8 and .375<=M.C_SCAN_ZERO;FO=M.pbits(0,62.9375,a);CB=ceil(FO/64.);CC=8*CB
	def FP(m_v):A=H(2,ceil(log2(m_v)));B=1<<A;C=3*(1<<A-2)if A>=3 else B;return C if A>=3 and C>=m_v else B
	def CD(m_v,m_base=1<<55):
		B=m_base;A=m_v
		if A<=B:return 2*(8*FP(A))+2*(8*A)
		C=(A+1)//2;return H(CD(C,B),16*(C+2)+CD(C+2,B))
	Ap=0
	for Du in C(61):
		FQ=O(ceil(CB/(1<<A))for A in C(1,Du+1));Dv=H(1,ceil(CB/(1<<Du)));Dw=8*FQ+8*(2*Dv)+CD(Dv)
		if Dw>Ap:Ap=Dw
	u,v,A4=M.NTT3_PRIMES;Dx=(1<<55)//2*(E*E)<u*v*A4;FR=K(u,-1,v);FS=K(u*v%A4,-1,A4);Dy={}
	for CE in M.NTT3_PRIMES:
		CF=3
		while K(CF,(CE-1)//2,CE)==1:CF+=1
		Dy[CE]=CF
	def w(A_v,N_v,p_v,inv_f=T):
		O=inv_f;E=N_v;A=p_v;L=K(Dy[A],(A-1)//E,A)
		if O:L=K(L,-1,A)
		B=[aN%A for aN in A_v]+[0]*(E-D(A_v));G=0
		for F in C(1,E):
			I=E>>1
			while G&I:G^=I;I>>=1
			G^=I
			if F<G:B[F],B[G]=B[G],B[F]
		H=2
		while H<=E:
			R=K(L,E//H,A);M=H>>1
			for F in C(0,E,H):
				N=1
				for J in C(M):P=B[F+J];Q=B[F+J+M]*N%A;B[F+J]=(P+Q)%A;B[F+J+M]=(P-Q)%A;N=N*R%A
			H<<=1
		if O:S=K(E,-1,A);B=[aO*S%A for aO in B]
		return B
	Aq=[0];BK=[0]
	def Ar(nb):
		Aq[0]+=nb
		if Aq[0]>BK[0]:BK[0]=Aq[0]
	def CG(nb):Aq[0]-=nb
	def FT(A_v,B_v,m_o):
		J=B_v;I=A_v;B=m_o;A=1<<H(2,ceil(log2(H(1,B))));Ar(16*A);D=w(I,A,u);F=w(J,A,u);P=w([D[aP]*F[aP]%u for aP in C(A)],A,u,b)[:B];Ar(8*B);D=w(I,A,v);F=w(J,A,v);Q=w([D[aQ]*F[aQ]%v for aQ in C(A)],A,v,b)[:B];Ar(8*B);D=w(I,A,A4);F=w(J,A,A4);R=w([D[aR]*F[aR]%A4 for aR in C(A)],A,A4,b)[:B];K=[0]*B;L=0
		for G in C(B):M=P[G];S=(Q[G]-M)*FR%v;N=M+u*S;T=(R[G]-N)*FS%A4;O=N+u*v*T+L;K[G]=O&E;L=O>>64
		CG(16*B+16*A);return K
	def BL(A_v,B_v,m_o,m_base=8):
		W=B_v;V=A_v;O=m_base;I=m_o
		if I<=O:return FT(V,W,I)
		e=(I+1)//2;B=(e+1)//2;P,Q=V[:B],V[B:];R,S=W[:B],W[B:];Z=BL(P,R,2*B,O);a=BL(Q,S,H(1,I-2*B),O);K=H(D(P),D(Q))+1;L=H(D(R),D(S))+1;M=K+L;Ar(8*(K+L+M));b=[0]*K;F=0
		for A in C(K):G=(P[A]if A<D(P)else 0)+(Q[A]if A<D(Q)else 0)+F;b[A]=G&E;F=G>>64
		c=[0]*L;F=0
		for A in C(L):G=(R[A]if A<D(R)else 0)+(S[A]if A<D(S)else 0)+F;c[A]=G&E;F=G>>64
		T=BL(b,c,M,O)
		for d in(Z,a):
			X=0
			for A in C(M):
				Y=T[A]-(d[A]if A<D(d)else 0)-X
				if Y<0:T[A]=Y+(1<<64);X=1
				else:T[A]=Y;X=0
		J=[0]*I
		for(A,N)in f(Z):J[A]=N
		for(A,N)in f(a):
			if 2*B+A<I:J[2*B+A]=N
		F=0
		for(A,N)in f(T):
			if B+A<I:G=J[B+A]+N+F;J[B+A]=G&E;F=G>>64
		U=B+M
		while F and U<I:G=J[U]+F;J[U]=G&E;F=G>>64;U+=1
		CG(8*(K+L+M));return J
	Dz=[N.randrange(1<<64)for aj in C(32)];D_=[N.randrange(1<<64)for ak in C(32)];Ar(1024);FU=BL(Dz,D_,64,m_base=8);CG(1024)
	def CH(L_v):
		A=0
		for B in C(D(L_v)-1,-1,-1):A=A<<64|L_v[B]
		return A
	FV=CH(FU)==CH(Dz)*CH(D_)and Aq[0]==0 and BK[0]<4096;CI=M._primes_upto(8200)[:1024];assert D(CI)==1024;A5=A9(CI);CJ=O(A.bit_length()for A in A5)
	while D(A5)>1:A5=[A5[ao]*A5[ao+1]for ao in C(0,D(A5),2)];CJ+=O(A.bit_length()for A in A5)
	E0=A5[0];BM=E0.bit_length();AS=0;As=0
	def CK(arr):
		A=arr;nonlocal AS,As
		if D(A)==1:return A[0]
		B=D(A)//2;C=CK(A[:B]);E=C.bit_length();AS+=E
		if AS>As:As=AS
		F=CK(A[B:]);AS-=E;return C*F
	FW=CK(CI);FX=Ap/CC;E1=FW==E0 and AS==0 and As<1.02*BM and Dx and FV and Ap<5.*CC and 8*CC<1<<64 and CJ>9.*BM
	def AT(x):return x&E
	def FY(a,b):A=(a&E)*(b&E);return A&E,A>>64&E
	def AD(x,n):x&=(1<<64*n)-1;return[x>>64*aS&E for aS in C(n)]
	def BN(L):
		A=0
		for B in C(D(L)-1,-1,-1):A=A<<64|L[B]&E
		return A
	def FZ(A):
		B=0
		while B<8 and A[B]==0:B+=1
		G=A[B];F=(G&-G).bit_length()-1;H=[0]*8
		for D in C(8):
			if D+B<8:I=A[D+B]>>F;J=A[D+B+1]<<64-F&E if F>0 and D+B+1<8 else 0;H[D]=I|J
		return H
	def CL(a):
		a&=E;A=1
		for B in C(6):A=AT(A*AT(2-AT(a*A)))
		return A
	E2=A.get(Cd,0)+8*(A.get(Ce,0)+8*A.get(Cf,0)+A.get(Cg,0)+7*A.get(Ch,0)+A.get(Ci,0))+A.get(Cj,0)+8*A.get(Ck,0)+A.get(Cl,0)+8*A.get(BU,0)+20;CM=[];BO={};E3=T
	for AN in AW.ko(Path(__file__).with_name('kernels_v32.s')).read().splitlines():
		S=AN.strip()
		if S=='//@seg cios8_setup':E3=b;continue
		if E3:
			if S.startswith(CX)or S.startswith(AJ)or not S or S.startswith(AX):
				if S==AJ and BO.get(BU)is not I:break
				continue
			if S.endswith(':'):BO[S[:-1]]=D(CM)
			else:V,Y=S.split(I,1);CM.append((V,Y))
	def B(k):return 32 if k=='xzr'else F(k[1:])
	W=[]
	for(V,Y)in CM:
		if V=='stp':Q,e,CN=Y.split(U,2);h=[au.strip()for au in CN.strip()[1:-1].split(U)];W.append((0,B(Q.strip()),B(e.strip()),B(h[0]),F(h[1][1:])//8 if D(h)>1 else 0))
		elif V in('ldr','str'):
			L,CN=Y.split(U,1);h=[ay.strip()for ay in CN.strip()[1:-1].split(U)];CO=V=='str'
			if D(h)==1:W.append((1 if CO else 2,B(L.strip()),B(h[0]),0,0))
			elif h[1].startswith('#'):W.append((1 if CO else 2,B(L.strip()),B(h[0]),F(h[1][1:])//8,0))
			else:W.append((3 if CO else 4,B(L.strip()),B(h[0]),B(h[1]),0))
		elif V=='mov':L,AE=[az.strip()for az in Y.split(U)];W.append((5,B(L),F(AE[1:]),0,0))
		elif V=='mul':L,Q,e=[aA.strip()for aA in Y.split(U)];W.append((6,B(L),B(Q),B(e),0))
		elif V=='umulh':L,Q,e=[aB.strip()for aB in Y.split(U)];W.append((7,B(L),B(Q),B(e),0))
		elif V=='adds':L,Q,e=[aC.strip()for aC in Y.split(U)];W.append((8,B(L),B(Q),B(e),0))
		elif V=='adc':L,Q,e=[aD.strip()for aD in Y.split(U)];W.append((9,B(L),B(Q),B(e),0))
		elif V=='add':L,Q,AE=[aE.strip()for aE in Y.split(U)];W.append((10,B(L),B(Q),F(AE[1:]),0))
		elif V=='sub':L,Q,AE=[aF.strip()for aF in Y.split(U)];W.append((11,B(L),B(Q),F(AE[1:]),0))
		elif V=='subs':L,Q,e=[aG.strip()for aG in Y.split(U)];W.append((12,B(L),B(Q),B(e),0))
		elif V=='sbcs':L,Q,e=[aH.strip()for aH in Y.split(U)];W.append((13,B(L),B(Q),B(e),0))
		elif V=='cmp':Q,AE=[aI.strip()for aI in Y.split(U)];W.append((14,B(Q),F(AE[1:]),0,0))
		elif V==BT:W.append((15,BO[Y.strip()],0,0,0))
		elif V=='cbnz':Q,Fa=[aJ.strip()for aJ in Y.split(U)];W.append((16,B(Q),BO[Fa],0,0))
		elif V=='csel':L,Q,e,AR=[aK.strip()for aK in Y.split(U)];assert AR=='hs';W.append((17,B(L),B(Q),B(e),0))
		else:raise Aw(V)
	def E4(A,B,N,n0inv):
		L={0:A9(A),1:A9(B),2:A9(N),4:[0]*10,17:[0]*8};C=[0]*33;C[0]=0;C[1]=1;C[2]=2;C[3]=n0inv&E;C[4]=4;C[7]=3735928559;C[8]=99;C[17]=17;K=0;M=0;Q=0;U=D(W)
		while M<U:
			G,F,H,I,R=W[M];M+=1;Q+=1
			if G==0:J=L[I];J[R]=C[F];J[R+1]=C[H]
			elif G==1:J=L[H];assert 0<=I<D(J);J[I]=C[F]
			elif G==2:J=L[H];assert 0<=I<D(J);C[F]=J[I]if F<32 else 0
			elif G==3:J=L[H];O=C[I];assert 0<=O<D(J);J[O]=C[F]
			elif G==4:J=L[H];O=C[I];assert 0<=O<D(J);C[F]=J[O]if F<32 else 0
			elif G==5:C[F]=H
			elif G==6:C[F]=C[H]*C[I]&E
			elif G==7:C[F]=C[H]*C[I]>>64&E
			elif G==8:P=C[H]+C[I];K=1 if P>E else 0;C[F]=P&E if F<32 else 0
			elif G==9:C[F]=C[H]+C[I]+K&E if F<32 else 0
			elif G==10:C[F]=C[H]+I&E
			elif G==11:C[F]=C[H]-I&E
			elif G==12:S,T=C[H],C[I];K=1 if S>=T else 0;C[F]=S-T&E if F<32 else 0
			elif G==13:P=C[H]+(C[I]^E)+K;K=1 if P>E else 0;C[F]=P&E if F<32 else 0
			elif G==14:K=1 if C[F]>=H else 0
			elif G==15:
				if K==0:M=F
			elif G==16:
				if C[F]!=0:M=H
			elif G==17:C[F]=C[H]if K else C[I]
		assert Q<=E2;return L[17]
	def BP(A,B):
		F=[0]*8;D=0
		for E in C(8):
			G=A[E]-B[E]-D
			if G<0:F[E]=G+(1<<64);D=1
			else:F[E]=G;D=0
		return F,D
	def AU(X):return[X[aT]>>1|(X[aT+1]&1)<<63 for aT in C(7)]+[X[7]>>1]
	def E5(U_in,V_in):
		A=U_in[:];B=V_in[:]
		if not q(A!=0 for A in B):return A,0
		C=0;E=0
		while b:
			C+=1
			if C>1024:return I,C
			if E==0:
				if B[0]&1==0:
					B=AU(B)
					if not q(A!=0 for A in B):return I,C
				else:
					D,F=BP(B,A)
					if not F:
						B=AU(D)
						if not q(A!=0 for A in B):return A,C
					else:
						D,G=BP(A,B)
						if G:return I,C
						A=AU(D)
						if not q(A!=0 for A in A):return I,C
						E=1
			elif A[0]&1==0:
				A=AU(A)
				if not q(A!=0 for A in A):return I,C
			else:
				D,F=BP(A,B)
				if not F:
					A=AU(D)
					if not q(A!=0 for A in A):return B,C
				else:
					D,G=BP(B,A)
					if G:return I,C
					B=AU(D)
					if not q(A!=0 for A in B):return I,C
					E=0
	def Fb(N_odd_L,g_L):
		L=CL(g_L[0]);B=N_odd_L[:];I=[0]*8
		for A in C(8):
			J=AT(B[A]*L);I[A]=J;K=0;F=0
			for D in C(8-A):
				M,N=FY(J,g_L[D]);G=M+K;K=N+(G>>64);G&=E;H=B[A+D]-G-F
				if H<0:B[A+D]=H+(1<<64);F=1
				else:B[A+D]=H;F=0
			if B[A]!=0:return
		return I
	AF=0;AV=0;E6=M._primes_upto(64);E7=1
	for z in E6:E7*=z
	At,AG=1,1
	while At+AG<1<<512:At,AG=AG,At+AG
	A6,Fc=(AG,At)if AG&1 else(At,AG);Fd=[1<<480,1<<511,3**302,2**199*3**120*5**40*1000000007,(1<<479)-1,(1<<512)-2,A6,A6*2 if A6*2<1<<512 else A6-1]+[N.randrange(1<<470,1<<512)for ap in C(56)]
	for E8 in Fd:
		Au=FZ(AD(E8,8));CP=BN(Au)
		if CP&1!=1:AF+=1;continue
		Fe=E7%CP
		if CP==1:E9=1
		else:
			Ff=AT(-CL(Au[0]));BQ=AD(Fe,8)
			for _ in C(9):BQ=E4(BQ,BQ,Au,Ff)
			EA,CQ=E5(Au,BQ)
			if CQ>AV:AV=CQ
			if EA is I or CQ>1024:AF+=1;continue
			EB=Fb(Au,EA)
			if EB is I:AF+=1;continue
			E9=BN(EB)
		t=E8
		for z in E6:
			while t%z==0:t//=z
		if E9!=t:AF+=1
	for(EC,ED)in((A6,0),(A6,A6),(A6,Fc),((1<<512)-1,(1<<512)-3),((1<<511)+1,1),(((1<<255)-19)*((1<<255)-31),((1<<255)-19)*((1<<255)-101))):
		EE,CR=E5(AD(EC,8),AD(ED,8))
		if CR>AV:AV=CR
		if EE is I or CR>1024 or BN(EE)!=AH.gcd(EC,ED):AF+=1
	Fg=Cm,Cn,Co,BV,Cd,Ce,Cf,Cg,Ch,Ci,Cj,Ck,Cl,BU,Cp,AY,BW,BX,BY,Cq,Cr,Cs,Ct,BZ,Ba,Cu,Cv;Fh=A8(A.get(B,0)>0 for B in Fg)and(A.get(AY)==A.get(Cq)==1 and A.get(BW)==A.get(Cr)==26 and A.get(BX)==A.get(Cs)==17 and A.get(BY)==A.get(Ct)==9 and A.get(BZ)==A.get(Ba)==4);Fi=8*A.get(Cm,0)+A.get(Cn,0)+8*A.get(Co,0)+10;Fj=A.get(BV,0)+2+9*E2
	for(CS,EF,EG)in(((1<<512)-1,(1<<512)-2,(1<<512)-2),((1<<512)-3,(1<<512)-4,(1<<512)-4),((1<<511)+1,1<<511,1<<511)):
		EH=AD(CS,8);Fk=AT(-CL(EH[0]));Fl=BN(E4(AD(EF,8),AD(EG,8),EH,Fk))
		if Fl!=EF*EG*K(1<<512,-1,CS)%CS:AF+=1
	EI=H(A.get(AY,0)+A.get(BW,0),A.get(AY,0)+9+A.get(BX,0),A.get(AY,0)+A.get(BY,0));Fm=A.get(Cp,0)+1024*EI+H(A.get(BZ,0),A.get(Ba,0))+1;Fn=A.get(BV,0)+8*A.get(Cu,0)+36*A.get(Cv,0)+20;EJ=Fi+Fj+Fm+Fn;EK=Fh and AF==0 and AV<=1024 and EJ<=M.C_SRT_G_LEAF
	def CT(L):
		A=0
		for B in C(D(L)-1,-1,-1):A=A<<64|L[B]
		return A
	def Fo(A,B):
		K,H=D(A),D(B);F=[0]*(K+H)
		for G in C(K):
			M=A[G];I=0
			for J in C(H):L=M*B[J]+F[G+J]+I;F[G+J]=L&E;I=L>>64
			F[G+H]=I
		return F
	def EL(A,B):
		I=H(D(A),D(B));F=[0]*(I+1);J=0
		for G in C(I):K=(A[G]if G<D(A)else 0)+(B[G]if G<D(B)else 0)+J;F[G]=K&E;J=K>>64
		F[I]=J
		while D(F)>1 and F[-1]==0:F.pop()
		return F
	def Av(A,B):
		E=[0]*D(A);G=0
		for F in C(D(A)):
			H=A[F]-(B[F]if F<D(B)else 0)-G
			if H<0:E[F]=H+(1<<64);G=1
			else:E[F]=H;G=0
		assert G==0
		while D(E)>1 and E[-1]==0:E.pop()
		return E
	def EM(A,B):
		E,F=D(A),D(B)
		while E>0 and A[E-1]==0:E-=1
		while F>0 and B[F-1]==0:F-=1
		if E!=F:return 1 if E>F else-1
		for G in C(E-1,-1,-1):
			if A[G]!=B[G]:return 1 if A[G]>B[G]else-1
		return 0
	def BR(A,B):
		Q=H(D(A),D(B))
		if Q<=16:return Fo(A,B)
		F=Q//2;G,I=A[:F],A[F:];J,K=B[:F],B[F:];R=BR(G,J);S=BR(I,K);T=EM(I,G);U=EM(K,J);V=Av(I,G)if T>=0 else Av(G,I);W=Av(K,J)if U>=0 else Av(J,K);assert D(V)<=H(D(G),D(I))and D(W)<=H(D(J),D(K));X=BR(V,W);Y=EL(R,S);a=Av(Y,X)if T*U>=0 else EL(Y,X);C=[0]*(D(A)+D(B))
		for(N,Z)in((0,R),(F,a),(2*F,S)):
			L=0
			for(P,b)in f(Z):
				if N+P<D(C):M=C[N+P]+b+L;C[N+P]=M&E;L=M>>64
			O=N+D(Z)
			while L and O<D(C):M=C[O]+L;C[O]=M&E;L=M>>64;O+=1
		return C
	CU=0
	for CV in C(200):
		BS=(4,8,12,16,20,24,32,48,64)[CV%9];EN=[E]*BS if CV<9 else[N.randrange(1<<64)for av in C(BS)];EO=[E]*BS if CV<9 else[N.randrange(1<<64)for aw in C(BS)];Fp=CT(BR(EN,EO));Fq=CT(EN)*CT(EO)
		if Fp!=Fq:CU+=1
	CW=A.get('mul4_sb',-1);EP=CU==0 and CW==105 and CW<=F(M._kar_sb(4));o=EU(M,A,N);EQ=o[AL];x=EV(M,B1);ER=x[AL];import descent_bn254_v6 as Fr;A7=Fr.verify_v40_descent(M,A);ES=A7[AL];G2S=verify_gsrt_2stage(M,seed=seed);Fs=Cx[AL]and Bf==0 and Bg==0 and D4<=M.OP[Cw]and A0==0 and Bj==768 and Bu==0 and Dn and C9==0 and Dt and E1 and EK and EP and EQ and ER and ES and G2S[AL];return g(ok=Fs,gsrt_2st=G2S,gsrt_2st_ok=G2S[AL],adic3_exact=Cx,ntt3_split_proto_errs=Bf,ctail_lazy_errs=Bg,ctail_lazy_ops=D4,ctail_budget=M.OP[Cw],bfi_b2_errs=A0,bfi_b2_cases=Bj,b2_prefix_mult=G(Bh,4),b2_active_ks=B4,delta_b1_errs=Bu,ng_diff_ub_bits=G(DN,6),u8_scale=j,u8_thr_max=B_,u8_init_min=DP,u8_max_val=DQ,ret2d_m_margin_bits=G(m['m_margin'],6),ret2d_rg_marg_err=G(DR,8),ret2d_sf_marg_err=G(DS,8),ret2d_per_prime_errs=BD,ret2d_mc_ok=Dg,ret2d_mc_emp_rg=G(BF,4),ret2d_mc_mod_rg=G(De,4),ret2d_mc_emp_sf=G(BG,4),ret2d_mc_mod_sf=G(Df,4),rho_tail_ok=DY,rho_335=f"{DW:.6e}",rho_400=f"{DX:.6e}",rot_tab_2s_errs=BH,step1_tab_8192_inv_errs=C5,sf_uses_survivor_upper_ok=Dm,u16_range_ok=Dn,scan8_u16_errs=C9,b1_delta_ops_ok=Dt,boot_df_tree_ok=E1,boot_df_stack_ratio=G(As/BM,4),boot_df_peak_ratio=G(FX,4),boot_df_prod_peak_log2=G(log2(Ap),4),boot_kar_ntt3_live_peak_ratio=G(BK[0]/512,4),ntt3_crt_no_ovf=Dx,boot_full_tree_ratio=G(CJ/BM,2),srt_g_leaf_ok=EK,srt_g_leaf_9sq_wc_ops=F(EJ),srt_g_leaf_budget=M.C_SRT_G_LEAF,bingcd8_max_iters=AV,bingcd8_iter_wc_ops=EI,mul_kar_sb_ok=EP,mul_kar_sb_errs=CU,mul4_sb_ops=CW,v37_ntt_ok=EQ,v37_r4_inner=o['r4_in'],v37_r4_hoist_x2=o['r4_ho'],v37_r3=o['r3'],v37_r5=o['r5'],v37_r4_blend=o['c_bl'],v37_r3_lin=o['c_r3'],v37_r5_lin=o['c_r5'],v37_bf_errs=o['bf_errs'],v37_conv_errs=o['conv_errs'],v39_emit8_ok=ER,v39_emit8_flush_amort=x['flush_amort'],v39_emit8_b1cpu_errs=x['cpu_errs'],v39_max_rows_per_Cs=x['max_rows'],v39_kmax_525=x['k525'],v39_kmax_53=x['k53'],v39_kmax_555=x['k555'],v39_sink_ub=x['sink_ub'],v39_sink_reset_ok=x['sink_reset_ok'],v40_descent_ok=ES,v43_cios5_ir=A7['cios5'],v43_mont2_ir=A7['mont2'],v43_bingcd_gen_ir=A7['bingcd_gen'],v43_psurv_iv_ok=A7['psurv_iv_ok'],v43_kron_ntt3_ok=A7['kron_ntt3_ok'],v43_bls_mpe_ok=A7['bls_mpe_ok'],v43_ecm_ram_mib=A7['ecm_ram_mib'])
def EU(M,segs,rng):
	f=segs;N=rng;A7=(1<<64)-1;l,m=f.get('r4_4clamp_inner_tw',-1),f.get('r4_4clamp_hoist_tw_x2',-1);n,o=f.get('r3_shoup_bfly',-1),f.get('r5_shoup_bfly',-1);t=.5*(9.*l/8.+9.*m/16.);u,v=9.*n/3.+22e1,9.*o/5.+22e1
	def w(n):
		if n<2:return T
		for D in(2,3,5,7,11,13,17,19,23,29,31,37):
			if n==D:return b
			if n%D==0:return T
		B,E=n-1,0
		while B&1==0:B>>=1;E+=1
		for F in(2,3,5,7,11,13,17,19,23,29,31,37):
			A=K(F,B,n)
			if A in(1,n-1):continue
			for G in C(E-1):
				A=A*A%n
				if A==n-1:break
			else:return T
		return b
	def A9(p):
		F,C,A=p-1,set(),p-1
		for D in(2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61):
			while A%D==0:C.add(D);A//=D
		B=67
		while B*B<=A:
			while A%B==0:C.add(B);A//=B
			B+=2
		if A>1:C.add(A)
		E=2
		while q(K(E,F//A,p)==1 for A in C):E+=1
		return E
	def F(v,w,wp,p):return v*w-(v*wp>>64)*p&A7
	def E(x,p2):return x-p2 if x>=p2 else x
	def x(a,ws,J,Jp,p):B=ws;A=2*p;C=E(a[0],A);D=F(a[1],B[0],(B[0]<<64)//p,p);G=F(a[2],B[1],(B[1]<<64)//p,p);H=F(a[3],B[2],(B[2]<<64)//p,p);I,K=E(C+G,A),E(C+A-G,A);L,M=E(D+H,A),F(D+A-H,J,Jp,p);return I+L,K+M,I+A-L,K+A-M
	def y(a,ws,c1,c1p,c2v,c2p,p):A=2*p;B=E(a[0],A);C=F(a[1],ws[0],(ws[0]<<64)//p,p);D=F(a[2],ws[1],(ws[1]<<64)//p,p);G=E(C+D,A);J=F(G,c1,c1p,p);H=F(C+A-D,c2v,c2p,p);I=E(B+J,A);return B+G,I+H,I+A-H
	def z(w3,p):A=(p+1)//2;B=(p-A)%p;C=(w3-w3*w3%p)*A%p;return B,(B<<64)//p,C,(C<<64)//p
	def A0(w5,p):A=w5;C=A*A%p;D=C*A%p;F=D*A%p;B=(p+1)//2;G=B*B%p;I=-G%p;J=(A+F-C-D)*G%p;H=(A-F)*B%p;E=(C-D)*B%p;K,L,M=E,(H-E)%p,-(H+E)%p;return X((A,(A<<64)//p)for A in(I,J,K,L,M))
	def A1(a,ws,cs,p):A=2*p;(S,T),(U,V),(W,X),(Y,Z),(b,c)=cs;D=E(a[0],A);B=[F(a[aD+1],ws[aD],(ws[aD]<<64)//p,p)for aD in C(4)];G,H=E(B[0]+B[3],A),E(B[1]+B[2],A);I,J=E(B[0]+A-B[3],A),E(B[1]+A-B[2],A);K=E(G+H,A);L=E(D+F(K,S,T,p),A);M=F(G+A-H,U,V,p);N,O=E(L+M,A),E(L+A-M,A);P=F(I+J,W,X,p);Q,R=E(P+F(I,Y,Z,p),A),E(P+F(J,b,c,p),A);return D+K,N+Q,O+R,O+A-R,N+A-Q
	Y,R=M.NTT2_52_PRIMES;AA=A8(w(A)and(A-1)%(3*(1<<55))==0 and 4*A<1<<64 for A in M.NTT3_PRIMES)and A8(w(A)and(A-1)%(15*(1<<40))==0 and 4*A<1<<64 for A in M.NTT2_52_PRIMES)and Y*R>(1<<18)*((1<<52)-1)**2*3;A2={aa:A9(aa)for aa in M.NTT3_PRIMES+M.NTT2_52_PRIMES};d=0
	for A in M.NTT3_PRIMES+M.NTT2_52_PRIMES:
		P,p=4*A,A2[A];i=K(p,(A-1)//4,A);AB=(i<<64)//A
		for S in C(516):
			J=[P-1 if S>>ad&1 else 0 for ad in C(4)]if S<16 else[N.randrange(P)for ae in C(4)];U=[N.randrange(A)for ab in C(3)];V=x(J,U,i,AB,A);B=[J[0]%A]+[J[af+1]*U[af]%A for af in C(3)]
			if H(V)>=P or[au%A for au in V]!=[(B[0]+B[2]+B[1]+B[3])%A,(B[0]-B[2]+i*(B[1]-B[3]))%A,(B[0]+B[2]-B[1]-B[3])%A,(B[0]-B[2]-i*(B[1]-B[3]))%A]:d+=1
		Z=K(p,(A-1)//3,A);AC=z(Z,A)
		for S in C(258):
			J=[P-1 if S>>ag&1 else 0 for ag in C(3)]if S<8 else[N.randrange(P)for ah in C(3)];U=[N.randrange(A)for ac in C(2)];V=y(J,U,*AC,A);B=[J[0]%A]+[J[ai+1]*U[ai]%A for ai in C(2)]
			if H(V)>=P or[av%A for av in V]!=[(B[0]+B[1]+B[2])%A,(B[0]+Z*B[1]+Z*Z*B[2])%A,(B[0]+Z*Z*B[1]+Z*B[2])%A]:d+=1
		if(A-1)%5==0:
			A3=K(p,(A-1)//5,A);AD=A0(A3,A)
			for S in C(282):
				J=[P-1 if S>>aw&1 else 0 for aw in C(5)]if S<32 else[N.randrange(P)for ax in C(5)];U=[N.randrange(A)for aj in C(4)];V=A1(J,U,AD,A);B=[J[0]%A]+[J[ay+1]*U[ay]%A for ay in C(4)]
				if H(V)>=P or[aB%A for aB in V]!=[O(B[C]*K(A3,aC*C,A)for C in C(5))%A for aC in C(5)]:d+=1
	def Q(a,p,inv=T):
		B=D(a);A=K(A2[p],(p-1)//B,p)
		if inv:A=K(A,p-2,p)
		def O(v,ws):
			B=D(v)
			if B==1:return[v[0]]
			A=5 if B%5==0 else 3 if B%3==0 else 4 if B%4==0 else 2;N=B//A;U=[O(v[aH::A],K(ws,A,p))for aH in C(A)];P=[0]*B;G=1;L=K(ws,N,p);V=z(L,p)if A==3 else I;W=A0(L,p)if A==5 else I;X=(L<<64)//p if A==4 else 0
			for Q in C(N):
				H=[G]
				for Y in C(A-2):H.append(H[-1]*G%p)
				J=[U[aI][Q]for aI in C(A)]
				if A==5:M=A1(J,H,W,p)
				elif A==3:M=y(J,H,*V,p)
				elif A==4:M=x(J,H,L,X,p)
				else:R,S=E(J[0],2*p),F(J[1],G,(G<<64)//p,p);M=R+S,R+2*p-S
				for T in C(A):P[Q+T*N]=M[T]
				G=G*ws%p
			return P
		G=O([aF%p for aF in a],A)
		if inv:H=K(B,p-2,p);return[aG*H%p for aG in G]
		return[aE%p for aE in G]
	j=0;AE=K(Y,R-2,R);e=1<<52
	for L in(12,20,24,30,48,60):
		W=L//2
		for A4 in(0,1):
			a=([e-1]*W if A4==0 else[N.randrange(e)for az in C(W)])+[0]*(L-W);c=([e-1]*W if A4==0 else[N.randrange(e)for aA in C(W)])+[0]*(L-W);A5=Q([ak*al%Y for(ak,al)in h(Q(a,Y),Q(c,Y))],Y,b);AF=Q([am*an%R for(am,an)in h(Q(a,R),Q(c,R))],R,b);r=s=0
			for k in C(L):A6=A5[k]+Y*((AF[k]-A5[k])*AE%R)+s;r|=(A6&e-1)<<52*k;s=A6>>52
			r|=s<<52*L
			if r!=O(a[A]<<52*A for A in C(W))*O(c[A]<<52*A for A in C(W)):j+=1
	for A in M.NTT3_PRIMES:
		for L in(12,24):
			a,c=[N.randrange(A)for ao in C(L)],[N.randrange(A)for ap in C(L)];AG=Q([aq*ar%A for(aq,ar)in h(Q(a,A),Q(c,A))],A,b)
			if AG!=[O(a[A]*c[(at-A)%L]for A in C(L))%A for at in C(L)]:j+=1
	AH=l==45 and m==74 and n==36 and o==93 and t<=46.5 and u<=364. and v<=42e1 and AA and d==0 and j==0;return g(ok=AH,r4_in=l,r4_ho=m,r3=n,r5=o,c_bl=G(t,4),c_r3=G(u,2),c_r5=G(v,2),bf_errs=d,conv_errs=j)
def EV(M,kernels_path):
	import stor_bound as B,v20_lat_b2 as Y,v25_verify as A9;from fractions import Fraction as E;from math import isqrt as AD;AE=M.verify_exhaustive_fd(R2=5e1);j=E(3,100000);Ab=E(B.G[0][0],B.G_SC)-j*E(88,100);Ac=(E(1)-j)*Ab*B.G_SC;l=(E(1)+B.DELTA)*B.G_SC*B.R_A_MAX**2;V={}
	for(AF,AG)in((105,2),(53,1),(107,2),(54,1),(109,2),(111,2),(123,2)):
		m=3*AG;AH=AF+m*64;a,n=0,1<<AH//m+2
		while n-a>1:
			o=(a+n)//2
			if o**m<=1<<AH:a=o
			else:n=o
		c=Ac*E(a,1<<64)
		if c<=2*l:return g(ok=T,flush_amort=99.,cpu_errs=1,max_rows=0,k525=0,k53=0,k555=0,sink_ub=0,sink_reset_ok=T)
		N=F(c//(c-2*l))
		if c>4*l and N!=1:return g(ok=T,flush_amort=99.,cpu_errs=1,max_rows=0,k525=0,k53=0,k555=0,sink_ub=0,sink_reset_ok=T)
		V[G(AF/AG,4)]=N
	p=B.main();q=p['max_rows_per_Cs'];r=p['epoch_buf_bytes']//8;Ad=(52.5,53.,13),(53.,53.5,5),(53.5,54.,3),(54.,55.5,2),(55.5,61.5,1);AI=P;AK=Y.EMIT8_P_MIN==6369051672525773 and Y.EMIT8_P_MIN*Y.EMIT8_P_MIN>1<<105
	for(AM,Ae,N)in Ad:
		s=r//(6*N*q)
		if V[AM]>N or s<1 or s*6*N*q>r:AK=T
		AI+=ceil(2.*M.KH_G_CAP*(2.**M.log2_li(Ae)-2.**M.log2_li(AM))/s)
	Af=2.**94.889;AN=AI*Z(B.N_SC)*5e1/(Af*1.07*log(61.5/52.5));L={};d=[];Q=I;e=[];t=T
	for AO in AW.ko(kernels_path):
		u=AO.strip()
		if u.startswith(Cc):Q=u.split()[1];t=T;continue
		if u.startswith(AJ):
			if Q==k:t=b
			Q=I;continue
		R=AO.split(AX,1)[0].strip()
		if not R or R.endswith(':')or R.startswith('.'):continue
		if t and D(e)<2:e.append(R)
		if Q in(A6,A7,AA,AB,AC,k):L.setdefault(Q,[]).append(R)
		if Q=='body_2':d.append(R)
	AP=(E(1)+B.DELTA)*B.W_A*B.G_SC*B.R_A_MAX**2;AQ=E(753729000);Ag=E(74,100)*AQ;AR=2*(AD(F(AP//AQ))+1)+3;S=2*(AD(F(AP//Ag))+1)+3;AS=D(d)>=2 and d[-2]=='ldr x8, [sp, #416]'and d[-1]=='str x8, [x17]'and AR*S<=r;Ah=e==['b.pl .Lk0_top','b .Ltop_2'];v=L.get(A6,[]);w=L.get(A7,[]);AT=L.get(AA,[]);W=v+w+AT+L.get(AB,[])+L.get(AC,[])+L.get(k,[])+e;Ai=A9._b1_decode(W,{'.Lk0_top':D(W),'.Ltop_2':D(W)+1});H={};x=1<<48;Aj=1<<50;y=1<<52;z=1<<45;A0=z+(1<<20);AU=(1<<34)-1;X=0,1,AU;A1={ab:(1<<54)+(1<<27)*(aa+1)for(aa,ab)in f(X)}
	for J in X:H[y+8*J]=A1[J]
	A2=4294967291&~(7<<23)|1<<23;Ak=A2-(1<<32)if A2>>31 else A2;H[z+(Ak<<3)&(1<<64)-1]=A0;H[A0+48]=4*(D(v)+D(w));H[A0+64]=4*D(v);A3=A9._B1CPU(Ai,H,{});A=A3.x;A[1],A[10],A[11],A[12],A[26],A[27],A[17]=8,z,48,3,16,1,y;A[4],A[13]=0xdeadbeefcafe1234,0x123456789abcdef0;Al=1,(1<<37)-1,1<<38,(1<<47)+12345,(1<<48)-1;A4={ac:[]for ac in X};Am=0,1,3,6
	for K in C(128):
		for(An,Ao)in f(Am):h=K*4+An;A5=K*7+Ao;i=X[h%3];AV=0 if i==0 else 43963475 if i==AU else 1+h%7;AY=Al[h%5];AZ=Aj+(h&3)*(1<<40);H[x+16*A5]=AY;H[x+16*A5+8]=AZ;H[AZ+8*A5]=i<<30|AV;A4[i].append((AV<<38)+AY)
	Aa={};U=0;A[23]=127
	for K in C(128):A[2],A[5],A[6],A[14],A[0]=x+16*(K*7-(2 if K==0 else 0)),K*7-(2 if K==0 else 0),4294967291,1,8 if K==0 else 6;A[23]=A[23]-1&(1<<64)-1;A3.N=bool(A[23]>>63);A3.run(0,{D(W),D(W)+1},{},Aa)
	if O(Aa.values())!=4483 or not Ah or D(w)!=2:U+=1
	if A[4]!=0xdeadbeefcafe1234 or A[13]!=0x123456789abcdef0 or A[1]!=8:U+=1
	for J in X:
		if(H[y+8*J]-A1[J])//8!=D(A4[J]):U+=1
		for(Ap,Aq)in f(A4[J]):
			if H.get(A1[J]+8*Ap)!=Aq:U+=1
	Ar=A8((A+1)**2<=937*A+1 for A in C(S+1))and(S+2)*S+1<<32<1<<64 and S==909;As=AE[AL]and E(str(AE['eps_G_hi']))<j and p[AL]and AK and AN<=Y.EMIT8_FLUSH_AMORT and D(AT)==8 and U==0 and Ar and AS;return g(ok=As,flush_amort=G(AN,5),cpu_errs=U,max_rows=q,k525=V[52.5],k53=V[53.],k555=V[55.5],sink_ub=AR*S,sink_reset_ok=AS)

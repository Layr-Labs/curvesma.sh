"";import math;from math import ceil,exp,floor,gcd,log,log2,prod,sqrt;from fractions import Fraction as F;from common12 import P2F,rho,integer_metric,regenerate_pair,norms,mul,norm,integer_interval;import ideal_canonical as IC;import embedding_bounds as EB;import filtsim;from threshold_model import Retention;LN2=log(2);UNIT=20000.
OP=dict(lb=69.5,lx2=64.,excess=1.05,slack_f=346.5,slack_b=276,width_bits=32,theta=F(3,8),alpha_f=.29,alpha_g=.75)
MARGIN=.06;CLS=F(47,50),F(103,100);V2E,KAP2=(2.,2.399985),(.06,.05)
def reg(cfg,a):s=CLS[sum(x*y for x,y in zip(a,P2F))%2==0];return F(cfg['production_region_lam'])*s,F(cfg['production_theta_region'])*s
def a2(r,s,op=OP):return op['alpha_'+'fg'[s]]+V2E[s]+KAP2[s]-r['v2'+'fg'[s]]
def kappa2(u,lb=69.5):
	lo,hi,r0=.0,1.,rho(u)
	for _ in range(40):
		m=(lo+hi)/2
		if sum(2.**-(k+1)*rho(u-(k-1-m)/lb)for k in range(60))>r0:lo=m
		else:hi=m
	return hi
MARGIN_SE=2.;Q_CORE=.02;DESCENT_CAP=121.9;GATE=7201208159453235316,'0.999031';NODE_CAP_LOG2=27.3;CANON_WORK_ALLOWANCE=1.00666;WORK_SE_K=2.;Q_F_SAFETY=10.;ECM_FIND=.999;REC_ROUTE,REC_ROUTE7=10.01,1.07;CELL_ROUTE=.02;LINE_TAB=[(18.78,4.91),(13.82,7.78),(9.82,12.83),(5.82,34.04),(1.83,167.01)];LINE_SE_K=2.;REC_TAB=[(18.91,23331674/522256),(10.91,3183663/41121),(3.91,1638258/6539),(.91,3216517/3967),(-1.09,4964184/2901)];C_COSET_EMPTY=1600.;C_COS=10.;C_APPLY,C_APPLY7=5.76,3.13;WIDE_FROM,WIDE_REC,WIDE_SETUP=63.,1.01,1.;REC_TAB_W=[(18.91,23449665/522256),(10.91,3216127/41121),(3.91,1745730/6539),(.91,3721383/3967),(-1.09,6474056/2901)];LOG2_WIDE_LAT=71.6;REC7_TAB=[(25.,33.17),(24.01,37.65),(19.,37.65),(18.01,38.33),(16.02,38.33),(16.,43.76),(14.49,49.66),(14.01,54.79),(13.01,54.95),(12.98,60.33),(12.03,60.33),(12.,62.88),(11.02,62.88),(11.01,73.74),(10.55,73.74),(10.31,82.85),(10.,86.39),(9.7,87.26),(9.55,87.26),(8.65,127.16),(7.47,127.16),(7.,167.75),(5.5,243.83),(4.99,243.83)];REC7_R2=1.29;REC7_K=[(24,.6102),(22,.6144),(20,.6744),(18,.7537),(16,.853),(14,.7948),(12,.8072),(10,.8279),(8,.7379),(4.99,.7398)];P7,X7,RUN_ALLOW,RUN_T=48.,64.,.00091,(32,2048);REL_BAND_LOSS=.009;COUNT=.998097,.999595
def hc_sup(t,xm):
	if xm<t[-1][0]:return.0
	cand=[xm]+[h for h,c in t if h<=xm]
	for(h1,c1),(h2,c2)in zip(t,t[1:]):
		s=(c1-c2)/(h1-h2);x=h2+(-s/log(2)-c2)/s
		if h2<=x<=min(h1,xm):cand.append(x)
	return max(2**x*tab(t,x)for x in cand)
def contam_sup(geo):
	S=geo.slice_cells(.0,1);lp=P7;d=.0
	while lp<X7:lp2=min(lp+.25,X7);m=.5*(lp+lp2);dp=NBAR/2*(2**log2_li(lp2)-2**log2_li(lp));wd=lp>=WIDE_FROM;d+=NBAR*log(lp2/lp)*2**m*(WIDE_REC*hc_sup(REC_TAB_W,log2(S)-m)if wd else hc_sup(REC_TAB,log2(S)-m))+dp*(C_COSET_EMPTY+C_SETUP6);lp=lp2
	return d
LOG2_PREP7=75.09
def step7(l,T=0,m=1):return tab(REC7_TAB,l)*(REC7_R2 if T>1024 else 1)*(next(k for e,k in REC7_K if l>=e)if m else 1)
C_SCAN_CELL,C_SCAN_LEAF=1.125,42.88;C_BUILD_E,C_BUILD_SLICE,C_BUILD_SEG=28.57,2.**21,2.**17;C_ZERO=.0;C_SETUP6=406.87+72.+21.13;C_NB_EVAL=1178.;C_NB_ROW=378.;C_NB_SLICE=18500.;MONT5=367.76;C_NORMS=48300.;C_TDIV=2.**24.5;NBAR=2.;LX0,P_MID=24.,48.;PRIME_MIN,EXPLICIT_BITS=2**24,14;C_SB=820000.;BAND_NBITS=832,576;BAND_FMAX=41,28;BAND_Q=.999;BAND_ROUNDS=2
def ecm_instr(fbits,nbits,success=ECM_FIND):import ecm16;return ecm16.ecm_instr(fbits,nbits,success)
def lam_top(u,t0,m=40):
	r=filtsim.rho(u)
	if r<=0 or t0>=1:return.0
	return sum(filtsim.rho(u-t)/r/t*((1-t0)/m)for t in(t0+(1-t0)*(i+.5)/m for i in range(m)))
def log2_fb(lb):return log2(NBAR/2*2**log2_li(lb)+2**31)-1
def log2_li(log2x):x=log2x*LN2;return log2x-log2(x)+log2(1+1/x+2/x/x)
def tab(t,log2h):
	if log2h>=t[0][0]:return t[0][1]
	for(h1,c1),(h2,c2)in zip(t,t[1:]):
		if log2h>=h2:return c2+(c1-c2)*(log2h-h2)/(h1-h2)
	return t[-1][1]
def skip(r):return r.get('p2a')and r.get('p2b')
def coprime_pair(a,b):
	d=gcd(norm(a),norm(b))
	if d==1:return True
	rows=[mul(v,[0]*k+[1])for v in(a,b)for k in range(6)];key=IC.row_hnf(rows,d);return all(key[k][k]==1 for k in range(6))
_CTX={}
def pool_row(cfg,sd,i):
	if not _CTX:_CTX.update(qa=integer_metric(cfg['matrix_a']),qb=integer_metric(cfg['matrix_b']),disks=EB.root_disks())
	(qa,da,ia),(qb,db,ib),disks=_CTX['qa'],_CTX['qb'],_CTX['disks'];A2=cfg['production_A_squared'];a,b,aq,bq=regenerate_pair(cfg,sd,i);p2={k:sum(x*y for x,y in zip(v,P2F))%2==0 for k,v in(('p2a',a),('p2b',b))};qa16=round(F(aq,da)/A2*65536);nf,ng=norms(a,b);nf16=round(16*log2(nf));ng16=round(16*log2(ng));p2.update(v2f=(nf&-nf).bit_length()-1,v2g=(ng&-ng).bit_length()-1);prim=coprime_pair(a,b);lam,kap=reg(cfg,a)[0],F(cfg['production_region_kap']);radius=(lam*A2-kap*F(aq,da))*db
	if F(bq)>radius:return dict(**p2,seed=sd,index=i,qa=qa16/65536.,nf=nf16/16.,ng=ng16/16.,primitive=prim,rejected=False,nonfinite=False,inreg=False)
	t=sum(qb[0][j]*b[j]for j in range(1,6));c=sum(qb[j][k]*b[j]*b[k]for j in range(1,6)for k in range(1,6))-radius;iv=integer_interval(qb[0][0],t,c);blo=b[0]//256*256;bhi=blo+255;lo,hi=max(iv.start,blo),min(iv.stop-1,bhi);rlen=16*max(1,(iv.stop-iv.start)//16);blen=hi-lo+1;est=EB.NormBounds(a,disks);box=box_rows(qb,radius,b);leaves=[l for l in two_level(est,[(est.row(list(s)),l,h)for s,l,h,_ in box],[q for*_,q in box],blo)if l[0]<=hi and l[1]>=lo and(l[3]is None or 0 in l[3])];at=[l for l in leaves if l[0]<=b[0]<=l[1]][0]
	return dict(**p2,seed=sd,index=i,qa=qa16/65536.,nf=nf16/16.,ng=ng16/16.,primitive=prim,rejected=isinstance(at[2],dict),nonfinite=False,evals=est.counts['block_evaluations'],bcells=sum(min(h,bhi)-max(l,blo)+1 for _,l,h,_ in box if l<=bhi and h>=blo),leaves=sum(1 if y-x>254 else y//64-x//64+1 for x,y in((max(l[0],lo),min(l[1],hi))for l in leaves)),blen=blen,rlen=rlen,inreg=True)
S0,L1=64,16
def two_level(est,rows,keys,blo):
	W=16*OP['width_bits'];ok=lambda bd:bd[0]is not None and bd[0][1]-bd[0][0]<=W;A=min(a for _,a,_ in rows);B=max(b for*_,b in rows);sub={};P=[(x,x+S0-1)for x in range(blo+256-S0,blo-1,-S0)];L=[];l2=0
	for k,q in enumerate(keys):sub.setdefault(q,[]).append(k)
	H=EB.hull([EB.hull([rows[k][0]for k in ks])for ks in sub.values()])
	while P:
		lo,hi=P.pop()
		if hi<A or lo>B:continue
		bd=est.block(H,lo,hi)
		if ok(bd):L.append((lo,hi,bd,None))
		elif hi-lo>=L1:m=(lo+hi)//2;P+=[(m+1,hi),(lo,m)]
		else:
			l2=1
			for ks in sub.values():L+=[(x,y,b2,tuple(ks)if k2 is None else(ks[k2],))for x,y,b2,k2 in EB.box_leaves(est,[rows[k]for k in ks],lo,hi,W)]
	if not l2:
		lo,hi=min(l[2][0][0]for l in L),max(l[2][0][1]for l in L)
		if hi-lo<=W:L=[(blo,blo+255,[[lo,hi],None],None)]
	return L
def pool_totals(s,n_eff,n_draws,geo):
	n,ys,y2,w,w2=s['rows'],s['ysum'],s['y2'],s['work'],s['work2'];mw=w/n;cond=sqrt(max(.0,(w2-n*mw*mw)/(n-1))/n)/mw
	return dict(log2_rels=log2(geo.va/n_eff*geo.vb*geo.lower*ys),yield_se_log2=sqrt(max(.0,y2/n_eff-(ys/n_eff)**2)/n_eff)/(ys/n_eff)/LN2,nf_rel=s['nfy']/ys,ng_rel=s['ngy']/ys,nf_mean=s['nfw']/s['w'],ng_mean=s['ngw']/s['w'],log2_population=log2(geo.va/n_eff*s['cells']),per_cell=w/s['cells'],log2_est=log2(geo.va/n_eff*w/UNIT),cond_se=cond,work_se=sqrt(cond**2+(1-n/n_draws)/n))
SUM_KEYS='rows','ysum','y2','nfy','ngy','work','work2','w','nfw','ngw','cells'
def add_sums(parts):return{k:sum(p[k]for p in parts)for k in SUM_KEYS}
BOX=8,8,8
def box_rows(qb,radius,b):
	import itertools;base=[b[j]//s*s for j,s in zip((2,3,4),BOX)];out=[]
	for d in itertools.product(*(range(s)for s in BOX)):
		x=b[1],base[0]+d[0],base[1]+d[1],base[2]+d[2],b[5];t=sum(qb[0][j]*x[j-1]for j in range(1,6));c=sum(qb[j][k]*x[j-1]*x[k-1]for j in range(1,6)for k in range(1,6))-radius;iv=integer_interval(qb[0][0],t,c)
		if iv:out.append((x,iv.start,iv.stop-1,(d[0]//2,d[1]//2,d[2]//4)))
	out.sort(key=lambda r:r[0]!=tuple(b[1:]));return out
class Geometry:
	def __init__(self,cfg,theta,lower=COUNT[0],upper=COUNT[1]):self.A2=cfg['production_A_squared'];self.lam,self.kap=float(F(cfg['production_region_lam'])),float(F(cfg['production_region_kap']));self.qa,self.da,self.ia=integer_metric(cfg['matrix_a']);self.qb,self.db,self.ib=integer_metric(cfg['matrix_b']);from ideal_canonical import determinant;det_a=F(determinant(self.qa),self.da**6);det_b=F(determinant(self.qb),self.db**6);v6=math.pi**3/6;self.va=v6*(float(theta)*self.A2)**3/(2*sqrt(float(det_a)));self.vb=v6*self.A2**3/(2*sqrt(float(det_b)));self.lower,self.upper=float(lower),float(upper);self.rb=[sqrt(float(self.ib[j][j])*self.A2*self.db)for j in range(6)];self.inv55=float(self.ib[5][5])*self.db
	def rad(self,q,y=0):return self.lam*float(CLS[y])-self.kap*q
	def slice_cells(self,q,y=0):return self.vb*self.rad(q,y)**3*self.upper
	def weight(self,qa_frac):return self.vb*(1-qa_frac)**3*self.upper
	def lines(self,q,y=0):return max(1.,floor(sqrt(self.rad(q,y)*self.A2*self.inv55)))
	def entries_upper(self,q,y=0):
		s=sqrt(self.rad(q,y));L4=1.
		for j in range(1,5):L4*=2*floor(self.rb[j]*s)+3
		return(floor(sqrt(self.rad(q,y)*self.A2*self.inv55))+2)*L4
def retentions(op=OP):
	return Retention(op['lb'],int(op['lx2']),prime_min=PRIME_MIN,explicit_bits=EXPLICIT_BITS,step=.25,max_bits=850,tilt=.05,odd=1),Retention(850,int(op['lx2']),prime_min=PRIME_MIN,explicit_bits=EXPLICIT_BITS,step=.25,max_bits=851,tilt=.005,odd=1)
def evaluate(rows,n_total,geo,op=OP,retain=None,cost=None):
	lb,lx2=op['lb'],op['lx2'];Sf,Sb,W=op['slack_f'],op['slack_b'],op['width_bits'];retain=retain or retentions(op)[0];cost=cost or retentions(op)[1];wsum=.0;ysum=.0;y2=.0;nf_m=ng_m=.0;nf_r=ng_r=.0;work_sum=work2=.0;c_ecm_f=ecm_instr(lb,Sb+8);c_mr_f=c_mr(Sb+8,MONT5);band_margin=1-BAND_FMAX[0]*(1-BAND_Q)**BAND_ROUNDS;pairs_all=NBAR/2*2**log2_li(lx2)+((WIDE_SETUP-1)*(NBAR/2)*(2**log2_li(lx2)-2**log2_li(WIDE_FROM))if lx2>WIDE_FROM else.0);rows=[r for r in rows if r.get('inreg',True)];D7=contam_sup(geo)
	for r in rows:
		w=(1-r['qa'])**3;wsum+=w
		if r['primitive']and not r['rejected']and r['qa']>=Q_CORE:rf=retain.at(r['nf']-r['v2f'],Sf,Sb)['retention'];rg=1.;uf_r,ug_r=(r['nf']+a2(r,0,op))/lb,(r['ng']+a2(r,1,op))/lb;miss=ECM_FIND**lam_top(uf_r,lx2/lb);y=rho(uf_r)*rho(ug_r)*rf*rg*miss*band_margin;ysum+=w*y;y2+=w*w*y*y;nf_r+=w*y*r['nf'];ng_r+=w*y*r['ng']
		nf_m+=w*r['nf'];ng_m+=w*r['ng']
		if skip(r):continue
		Y=bool(r.get('p2a'));fk=.5 if Y else 1.;S=geo.slice_cells(r['qa'],Y)*fk*geo.lower/geo.upper;K=geo.lines(r['qa'],Y)
		p_all=dict(build=(C_BUILD_E*geo.entries_upper(r['qa'],Y)+C_BUILD_SLICE)/S+C_BUILD_SEG*4/2**32,zero=.0,setup6=pairs_all*(C_SETUP6+C_COS*(fk<1))/S,thr=(C_NB_EVAL*r['evals']/r['bcells']+C_NB_ROW/r['rlen'])/fk+C_NB_SLICE/S)
		lp=LX0;line=.0
		while lp<P_MID:lp2=min(lp+.25,P_MID);dsum=log(lp2/lp);line+=NBAR*dsum*tab(LINE_TAB,log2(S/K)-.5*(lp+lp2));lp=lp2
		lp=P_MID;rec=.0;empty=.0;u7=h7=.0;lg=log2(geo.vb*geo.upper*fk/2)+r['lrs']if RUN_T[0]<=r.get('T',0)<=RUN_T[1]else None;l7=lg and lg+1+log2(geo.lower/geo.upper)
		while lp<lx2:
			lp2=min(lp+.25,lx2);dsum=log(lp2/lp);log2h=log2(S)-.5*(lp+lp2)
			if lg is not None and P7<=lp<X7 and l7-.5*(lp+lp2)>=REC7_TAB[-1][0]:dp=NBAR/2*(2**log2_li(lp2)-2**log2_li(lp));rec+=NBAR*dsum*step7(l7-.5*(lp+lp2),r['T']);h7+=NBAR*dsum;u7+=dp*((C_SETUP6+C_COS)*2**-l7-(C_SETUP6+C_COS*(fk<1))/S)
			elif log2h>=REC_TAB[-1][0]:rec+=NBAR*dsum*(WIDE_REC*tab(REC_TAB_W,log2h)if lp>=WIDE_FROM else tab(REC_TAB,log2h))
			else:dpairs=NBAR/2*(2**log2_li(lp2)-2**log2_li(lp));empty+=dpairs*C_COSET_EMPTY/S
			lp=lp2
		rec+=RUN_ALLOW*D7/S;p_all['setup6']=max(.0,p_all['setup6']+u7);hits=NBAR*log(lx2/P_MID);no=r['nf']-r['v2f'];sf=cost.survivor_upper(no,Sf,W);sa=cost.survivor_upper(no,Sf,.25);s0=cost.survivor_upper(no,Sf,.25,Sb,lb,KMAX);q_g=min(1.,Q_F_SAFETY*rho((r['ng']+a2(r,1,op))/lb));c_surv=C_GBATCH+q_g*(C_TDIV+C_BATCH+C_BATCH2+c_mr_f+c_ecm_f);line*=SIEVE_SIDES;rec*=SIEVE_SIDES;empty*=SIEVE_SIDES;hits=hits*SIEVE_SIDES
		p_keep=dict(line=line,rec=rec,empty=empty,apply=hits*C_APPLY-(C_APPLY-C_APPLY7)*h7*SIEVE_SIDES,route=CELL_ROUTE+REC_ROUTE*hits-(REC_ROUTE-REC_ROUTE7)*h7*SIEVE_SIDES,scan=C_SCAN_CELL+C_SCAN_LEAF*r['leaves']/r['blen']/fk,zero=C_ZERO,cof=sf*C_NORMS+sa*C_SB+s0*c_surv)
		per_cell=sum(p_all.values())+sum(p_keep.values());contrib=geo.weight(r['qa'])*per_cell;work_sum+=contrib;work2+=contrib*contrib
	nf_m/=wsum;ng_m/=wsum;return dict(sums=dict(rows=len(rows),ysum=ysum,y2=y2,nfy=nf_r,ngy=ng_r,work=work_sum,work2=work2,w=wsum,nfw=nf_m*wsum,ngw=ng_m*wsum,cells=sum(geo.weight(r['qa'])for r in rows if not skip(r))))
def _filt():import base64,lzma,json,os;d=json.loads(lzma.decompress(base64.b85decode(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'filt16.b85')).read().strip())));return tuple({int(k):tuple(map(tuple,v))for k,v in d[n].items()}for n in('FILT_RUNS','FILT_X','FILT_W'))
FILT_RUNS,FILT_X,FILT_W=_filt();NU_H=.02282;NU=.02284;SD_C,P_CAP=1.,1e-09;FILT_EXCESS=1.09;FILT_MEANS=740.3,498.7;W_ROW=251.;C_ENTRY_K=8.4;K5_RUNS={15:(8.426,8.44,8.401,8.419),16:(8.374,8.405,8.373,8.375,8.368,8.393,8.386,8.376),17:(8.341,8.345,8.341,8.349,8.355,8.363,8.348,8.34)};C_ROW=2710.;C_RED9=270.;BW_M=8192;BW_N=130;C_PROJ=60.;N_SM=10;C_FULL=97.;C_CHI=200.;LOG2_SM=88.;LOG2_FBC=73.1;LOG2_CANON=92.;LOG2_MISC=80.;LOG2_VERIFY=11.;LOG2_GTREE,LOG2_GSPLIT=74.,96.5
def linear_algebra(log2_U,excess=None):ratio=NU;log2_N=log2_U+log2(NU);c_prod=W_ROW*C_ENTRY_K+C_ROW+C_CHI;c_kry=c_prod;per_n2=c_kry*(BW_N/BW_M+1)+(1/BW_M+1/BW_N)*(c_prod+BW_N*C_FULL+C_RED9)+C_PROJ*(BW_N/BW_M+1);per_n1=(BW_M+BW_N)*C_RED9;log2_la=log2((2**(2*log2_N)*per_n2+2**log2_N*per_n1)/UNIT);return dict(ratio=ratio,log2_N=log2_N,log2_la=log2_la,log2_gen=generator(log2_N))
def generator(log2_N,x=0):
	m,k=BW_M,BW_M+BW_N;L=(2**log2_N+N_SM)*(1/m+1/BW_N)+2;J=ceil(log2(L))
	if x:return sum(2**j*(m+k)*k*k*(fft_F(18*l+18)+640*l)for j in range(J)for l in[-(-L//2**j)])+2**J*m*k*k*(C_FULL+C_RED9)
	return log2(4*k**3*log2(L)*(fft_F(18*L)+640*L)/UNIT)
def total_log2(log2_relcol,log2_la,log2_desc=DESCENT_CAP,log2_gen=.0):
	return log2(2**log2_relcol+2**log2_la+2**log2_desc+2**log2_gen+2**LOG2_SM+2**LOG2_FBC+2**LOG2_CANON+2**LOG2_MISC+2**LOG2_VERIFY+2**LOG2_WIDE_LAT+2**LOG2_PREP7+2**LOG2_SPLIT+2**LOG2_GTREE+2**LOG2_GSPLIT+2**50)
KMAX,MR_ROUNDS=4,40;SIEVE_SIDES=.5;C_GBATCH=2.**24.98;C_BATCH=2.**24;C_BATCH2=2.**25.25
def srt_toy(k,seed,dxs=(-2,2),G=None):
	import random;from threshold_model import primes_upto;ps=[p for p in primes_upto(2**21)if p>2**20];P=ptree(ps)[-1][0];g=random.Random(seed);zs=[]
	for _ in range(2**k):
		z=g.getrandbits(832)|1
		for _ in range(g.randrange(3)):z*=g.choice(ps)**g.randrange(1,3)
		zs.append(z)
	R=[srt(P,zs,G,dx)for dx in dxs];z=3**500*7*(2**31-1);return sum(a!=P%z for r in R for a,z in zip(r[0],zs))+sum(k2(z,a)!=part(z,P)for a,z in zip(R[0][0],zs))+(k2(z,105%z)!=3**500*7),max(r[1]for r in R)
def k2(z,r):v=(z&-z).bit_length()-1;z>>=v;return gcd(z,pow(r,2**10,z))<<v
def part(z,P):
	m=1;g=gcd(z,P)
	while g>1:z//=g;m*=g;g=gcd(z,g)
	return m
def ptree(zs):
	T=[list(zs)]
	while len(T[-1])>1:a=T[-1];T.append([a[i]*a[i+1]if i+1<len(a)else a[i]for i in range(0,len(a),2)])
	return T
def srt(P,zs,G=None,dx=0):
	T=ptree(zs);G=G if G is not None else len(T)+1;w=lambda Z:64*(-(-Z.bit_length()//64)+G);Z=T[-1][0];b=w(Z);s=b+P.bit_length()+8;x=(1<<s)//Z+dx;q=P*x>>s-b;r=(P<<b)-q*Z;fix=0
	while r<0:q-=1;r+=Z;fix+=1
	while r>=Z:q+=1;r-=Z;fix+=1
	Y=[q%(1<<b)]
	for l in range(len(T)-2,-1,-1):
		C=T[l];Yn=[]
		for i,z in enumerate(C):v=i//2;Zv=T[l+1][v];sib=C[i^1]if i^1<len(C)else 1;bv,bc=w(Zv),w(z);Yn.append(Y[v]*sib%(1<<bv)>>bv-bc)
		Y=Yn
	out=[]
	for z,y in zip(zs,Y):b=w(z);r=z*y+(1<<b)-1>>b;out.append(0 if r==z else r)
	return out,fix

FFT_A,FFT_B,FFT_ROOT,FFT_COPY,FFT_LEAF,FFT_XFER,FFT_LOCAL=8.5,1200.,6.,3.,2.**15,504.,2**26
def fft_F(N):return FFT_A*N*log2(max(N,2))**2+FFT_B+(FFT_XFER*N if N>FFT_LOCAL else.0)
def gbatch_cost(pbits=1.000028*2**69.5/LN2,leaf=9,L=None,c_leaf=FFT_LEAF):SP=ceil(pbits/64);L=L if L is not None else ceil(log2((SP+1)/leaf));G=L+2;n=lambda j:leaf*2**j;PT=sum(2**(L-j)*(fft_F(2*n(j-1))+FFT_COPY*2*n(j-1))for j in range(1,L+1));RT=sum(2**(L-j+1)*(fft_F(n(j)+G+n(j-1))+FFT_COPY*(n(j-1)+G))for j in range(1,L+1));return(PT+RT+FFT_ROOT*fft_F(2*(max(SP,n(L))+G)))/2**L+c_leaf
LOG2_SPLIT=101.27;SPLIT_LIMB,SPLIT_SIDE,SPLIT_LOOP=6,22,4;SPLIT_DIV,SPLIT_IDEAL=2.**10,2.**26
def band_split(log2_rels_est,yield_se,rmin_f):lf,lg=[-(-b//64)for b in BAND_NBITS];np_=1.25506*2**48/(48*LN2);c_rel=np_*(SPLIT_LIMB*(lf+lg)+2*SPLIT_SIDE+SPLIT_LOOP)+sum(BAND_FMAX)*SPLIT_DIV+SPLIT_IDEAL;bm=(1-BAND_FMAX[0]*(1-BAND_Q)**BAND_ROUNDS)*(1-BAND_FMAX[1]*(1-BAND_Q)**BAND_ROUNDS);lc=log2_rels_est+log2(COUNT[1]/COUNT[0])-(lf+lg)*log2(ECM_FIND)-log2(bm*rmin_f)+2*yield_se;return dict(log2_primes=log2(np_),c_rel_log2=log2(c_rel),log2_count=lc,log2_units=lc+log2(c_rel/UNIT))
def rmin(ret,slack,sb=None,top=850):
	v=[]
	for k in range(4,4*top+1):
		try:v.append(ret.at(k/4,slack,sb)['retention'])
		except ArithmeticError:pass
	return min(v)
def c_mr(nbits,mul):return MR_ROUNDS*2*nbits*mul
def toy(n=48,X=10.5,lb=11.5,l0=5,SA=30.,SB=22.5,kmax=2,trials=150000,seed=1):
	"";import random;from threshold_model import Retention,primes_upto;ps=primes_upto(int(2**lb));q0=2**l0;sv=[p for p in ps if q0<=p<2**X];Ps=prod(p for p in ps if p<q0);Pv=prod(sv);Pl=prod(p for p in ps if p>=2**X);kw=dict(prime_min=q0,explicit_bits=ceil(lb)+1,step=1/16);R=Retention(lb,X,max_bits=n+3,tilt=.05,**kw);C=Retention(n+2,X,max_bits=n+3,tilt=.005,**kw);c=C.cap(n);rng=random.Random(seed);f16=lambda v:floor(16*log2(v));A=Bw=sm=kp=0
	for _ in range(trials):
		x=rng.randrange(int(2**(n-.125)),int(2**(n+.125)));s=x;m=1;g=gcd(s,Ps)
		while g>1:s//=g;m*=g;g=gcd(s,g)
		v=gcd(s,Pv);r=s;g=gcd(r,Pv*Pl)
		while g>1:r//=g;g=gcd(r,g)
		sm+=r==1
		if log2(x/v)<=SA+1:
			d=f16(x)-sum(f16(p)for p in sv if v%p==0)
			if d<=16*SA:
				A+=1;O=d-f16(m)
				if O<=16*SB and any(16*X*k-1<=O<=16*lb*k+c+1 for k in range(kmax+1)):Bw+=1;kp+=r==1
	return dict(n=n,trials=trials,A=A,B=Bw,smooth=sm,kept=kp,sA=C.survivor_upper(n,SA,.25)*trials,sB=C.survivor_upper(n,SA,.25,SB,lb,kmax)*trials,ret=R.at(n,SA,SB)['retention'])

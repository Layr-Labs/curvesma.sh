"""stnfs ledger model: relation collection (chunk-routed records, every block"""
import math
from fractions import Fraction as F
from common12 import rho,integer_metric,quadratic,regenerate_pair,norms,mul,norm,integer_interval
import ideal_canonical as IC
import embedding_bounds as EB
import filtsim
from threshold_model import Retention
LN2=math.log(2)
UNIT=20000.0
OP=dict(lb=69.0,lx2=64.0,excess=1.05,slack_f=302,slack_g=233,width_bits=8,L0=64,
theta=F(3,8),alpha_f=5.15,alpha_g=0.85)
MARGIN=0.06
MARGIN_SE=2.0
Q_CORE=0.02
DESCENT_CAP=113.279
NODE_CAP_LOG2=27.3
CANON_WORK_ALLOWANCE=1.00806
WORK_SE_K=2.0
Q_F_SAFETY=100.0
ECM_FIND=0.999
REC_ROUTE=10.01
CELL_ROUTE=8.0
LINE_TAB=[(18.81,17.2),(13.82,24.1),(9.82,38.99),(5.82,88.32),(1.83,352.1)]
LINE_SE_K=2.0
REC_TAB=[(18.91,23331674/522256),(10.91,3183663/41121),(3.91,1638258/6539),(0.91,3216517/3967),(-1.09,4964184/2901)]
C_COSET_EMPTY=1600.0
C_APPLY=5.76
C_RECOVER=5.76
WIDE_FROM,WIDE_REC,WIDE_SETUP=63.0,1.01,1.0
REC_TAB_W=[(18.91,23449665/522256),(10.91,3216127/41121),(3.91,1745730/6539),(0.91,3721383/3967),(-1.09,6474056/2901)]
LOG2_WIDE_LAT=70.6
REC7_TAB=[(21.99,66.02),(19.99,72.59),(17.99,82.58),(15.99,98.16),(14.99,107.63),(13.99,121.21),(12.99,121.21),(11.99,136.87),(10.99,164.07),(9.99,198.04),(8.99,236.60),(7.99,278.37),(6.99,380.56),(5.99,417.68),(4.99,541.42)]
P7,X7,EXTRA7,RUN_ALLOW,RUN_T=56.0,64.0,2.0,0.00075,(32,1024)
REL_BAND_LOSS=0.009 # [11c]
def hc_sup(t,xm):
	if xm<t[-1][0]:
		return 0.0
	cand=[xm]+[h for h,c in t if h<=xm]
	for(h1,c1),(h2,c2)in zip(t,t[1:]):
		s=(c1-c2)/(h1-h2);x=h2+(-s/math.log(2)-c2)/s
		if h2<=x<=min(h1,xm):cand.append(x)
	return max(2**x*tab(t,x)for x in cand)
def contam_sup(geo):
	S=geo.slice_cells(0.0);lp=P7;d=0.0
	while lp<X7:
		lp2=min(lp+0.25,X7);m=0.5*(lp+lp2);dp=(NBAR/2)*(2**log2_li(lp2)-2**log2_li(lp));wd=lp>=WIDE_FROM
		d+=NBAR*math.log(lp2/lp)*2**m*(WIDE_REC*hc_sup(REC_TAB_W,math.log2(S)-m)if wd else hc_sup(REC_TAB,math.log2(S)-m))+dp*(C_COSET_EMPTY+C_SETUP6);lp=lp2
	return d
LOG2_PREP7=73.5
def step7(l):
	return next(c for h,c in REC7_TAB if l>=h)
C_SCAN_CELL,C_SCAN_LEAF=3.002,90.98
C_BUILD_E,C_BUILD_SLICE,C_BUILD_SEG=28.17,2.0**21,2.0**17
C_ZERO=0.0
C_SETUP6=406.87+72.0+21.13
C_NB_EVAL=1508.2
C_NB_ROW=4600.0
C_NB_SLICE=17500.0
MONT5=367.76
C_NORMS=2.0**22
C_TDIV=2.0**24.5
NBAR=2.0
LX0,P_MID=14.0,48.0
PRIME_MIN,EXPLICIT_BITS=16384,14 # unsieved below PRIME_MIN: omitted mass
BAND_FBITS=48.0
BAND_NBITS=(832,576)
BAND_FMAX=(41,28)
BAND_Q=0.999
BAND_ROUNDS=2
_WIN={}
def window(cost,lb):
	if id(cost)not in _WIN:_WIN[id(cost)]=Window(cost,lb)
	return _WIN[id(cost)]
def ecm_instr(fbits,nbits,success=ECM_FIND):
	import ecm16
	return ecm16.ecm_instr(fbits,nbits,success)
def lam_top(u,t0,m=40):
	r=filtsim.rho(u)
	if r<=0 or t0>=1:
		return 0.0
	return sum(filtsim.rho(u-t)/r/t*((1-t0)/m)for t in(t0+(1-t0)*(i+0.5)/m for i in range(m)))
def log2_li(log2x):
	x=log2x*LN2
	return log2x-math.log2(x)+math.log2(1+1/x+2/x/x)
def tab(t,log2h):
	if log2h>=t[0][0]:
		return t[0][1]
	for(h1,c1),(h2,c2)in zip(t,t[1:]):
		if log2h>=h2:
			return c2+(c1-c2)*(log2h-h2)/(h1-h2)
	return t[-1][1]
def coprime_pair(a,b):
	d=math.gcd(norm(a),norm(b))
	if d==1:
		return True
	rows=[mul(v,[0]*k+[1])for v in(a,b)for k in range(6)]
	key=IC.row_hnf(rows,d)
	return all(key[k][k]==1 for k in range(6))
_CTX={}
def pool_row(cfg,sd,i):
	if not _CTX:
		_CTX.update(qa=integer_metric(cfg['matrix_a']),qb=integer_metric(cfg['matrix_b']),disks=EB.root_disks())
	(qa,da,ia),(qb,db,ib),disks=_CTX['qa'],_CTX['qb'],_CTX['disks']
	A2=cfg['production_A_squared']
	a,b,aq,bq=regenerate_pair(cfg,sd,i)
	qa16=round(F(aq,da)/A2*65536)
	nf,ng=norms(a,b)
	nf16=round(16*math.log2(nf));ng16=round(16*math.log2(ng))
	prim=coprime_pair(a,b)
	lam,kap=F(cfg.get('production_region_lam','1')),F(cfg.get('production_region_kap','1'))
	radius=(lam*A2-kap*F(aq,da))*db
	if F(bq)>radius: # not a region row
		return dict(seed=sd,index=i,qa=qa16/65536.0,nf=nf16/16.0,ng=ng16/16.0,primitive=prim,rejected=False,nonfinite=False,inreg=False)
	t=sum(qb[0][j]*b[j]for j in range(1,6))
	c=sum(qb[j][k]*b[j]*b[k]for j in range(1,6)for k in range(1,6))-radius
	iv=integer_interval(qb[0][0],t,c)
	blo=(b[0]//64)*64;lo,hi=max(iv.start,blo),min(iv.stop-1,blo+63)
	rlen=16*max(1,(iv.stop-iv.start)//16);blen=hi-lo+1
	est=EB.NormBounds(a,disks)
	box=box_rows(qb,radius,b,blo)
	leaves=EB.box_leaves(est,[(est.row(list(s)),l,h)for s,l,h in box],min(x[1]for x in box),max(x[2]for x in box))
	leaves=[l for l in leaves if l[0]<=hi and l[1]>=lo and l[3]in(None,0)]
	at=[l for l in leaves if l[0]<=b[0]<=l[1]][0]
	return dict(seed=sd,index=i,qa=qa16/65536.0,nf=nf16/16.0,ng=ng16/16.0,primitive=prim,rejected=isinstance(at[2],dict),
	nonfinite=False,evals=est.counts['block_evaluations'],bcells=sum(h-l+1 for _,l,h in box),leaves=len(leaves),blen=blen,rlen=rlen,inreg=True)
def pool_totals(s,n_eff,n_draws,geo):
	n,ys,y2,w,w2=s['rows'],s['ysum'],s['y2'],s['work'],s['work2']
	mw=w/n
	cond=math.sqrt(max(0.0,(w2-n*mw*mw)/(n-1))/n)/mw
	return dict(log2_rels=math.log2(geo.va/n_eff*geo.vb*geo.lower*ys),
	yield_se_log2=math.sqrt(max(0.0,y2/n_eff-(ys/n_eff)**2)/n_eff)/(ys/n_eff)/LN2,
	nf_rel=s['nfy']/ys,ng_rel=s['ngy']/ys,nf_mean=s['nfw']/s['w'],ng_mean=s['ngw']/s['w'],
	log2_population=math.log2(geo.va/n_eff*s['cells']),per_cell=w/s['cells'],log2_est=math.log2(geo.va/n_eff*w/UNIT),
	cond_se=cond,work_se=math.sqrt(cond**2+(1-n/n_draws)/n))
SUM_KEYS=('rows','ysum','y2','nfy','ngy','work','work2','w','nfw','ngw','cells')
def add_sums(parts):
	return{k:sum(p[k]for p in parts)for k in SUM_KEYS}
BOX=(2,2,4)
def box_rows(qb,radius,b,blo):
	import itertools
	base=[(b[j]//s)*s for j,s in zip((2,3,4),BOX)];out=[]
	for d in itertools.product(*(range(s)for s in BOX)):
		x=(b[1],base[0]+d[0],base[1]+d[1],base[2]+d[2],b[5])
		t=sum(qb[0][j]*x[j-1]for j in range(1,6))
		c=sum(qb[j][k]*x[j-1]*x[k-1]for j in range(1,6)for k in range(1,6))-radius
		iv=integer_interval(qb[0][0],t,c);lo,hi=max(iv.start,blo),min(iv.stop-1,blo+63)
		if lo<=hi:
			out.append((x,lo,hi))
	out.sort(key=lambda r:r[0]!=tuple(b[1:]))
	return out
class Geometry:
	def __init__(self,cfg,theta,lower=F(2401,2500),upper=F(2601,2500)):
		self.A2=cfg['production_A_squared']
		self.lam,self.kap=float(F(cfg.get('production_region_lam','1'))),float(F(cfg.get('production_region_kap','1')))
		self.qa,self.da,self.ia=integer_metric(cfg['matrix_a'])
		self.qb,self.db,self.ib=integer_metric(cfg['matrix_b'])
		from ideal_canonical import determinant
		det_a=F(determinant(self.qa),self.da**6);det_b=F(determinant(self.qb),self.db**6)
		v6=math.pi**3/6
		self.va=v6*(float(theta)*self.A2)**3/(2*math.sqrt(float(det_a)))
		self.vb=v6*self.A2**3/(2*math.sqrt(float(det_b)))
		self.lower,self.upper=float(lower),float(upper)
		self.rb=[math.sqrt(float(self.ib[j][j])*self.A2*self.db)for j in range(6)]
		self.inv55=float(self.ib[5][5])*self.db
	def rad(self,qa_frac):
		return self.lam-self.kap*qa_frac
	def slice_cells(self,qa_frac):
		return self.vb*self.rad(qa_frac)**3*self.upper
	def weight(self,qa_frac):
		return self.vb*(1-qa_frac)**3*self.upper # sampling density
	def lines(self,qa_frac):
		return max(1.0,math.floor(math.sqrt(self.rad(qa_frac)*self.A2*self.inv55)))
	def entries_upper(self,qa_frac):
		s=math.sqrt(self.rad(qa_frac));L4=1.0
		for j in range(1,5):
			L4*=2*math.floor(self.rb[j]*s)+3
		return(math.floor(math.sqrt(self.rad(qa_frac)*self.A2*self.inv55))+2)*L4
def retentions(op=OP):
	return(Retention(int(op['lb']),int(op['lx2']),prime_min=PRIME_MIN,explicit_bits=EXPLICIT_BITS,step=.25,max_bits=850,tilt=.05),
	Retention(850,int(op['lx2']),prime_min=PRIME_MIN,explicit_bits=EXPLICIT_BITS,step=.25,max_bits=851,tilt=.005))
def evaluate(rows,n_total,geo,op=OP,retain=None,cost=None,scale_consts=1.0,verbose=False,return_work=False):
	lb,lx2,excess=op['lb'],op['lx2'],op['excess']
	Sf,Sg,W=op['slack_f'],op['slack_g'],op['width_bits']
	retain=retain or retentions(op)[0]
	cost=cost or retentions(op)[1]
	wsum=0.0;ysum=0.0;y2=0.0;rej_w=0.0;nf_m=ng_m=0.0;nf_r=ng_r=0.0
	work_sum=work2=0.0
	work_list=[]if return_work else None
	c_ecm_f,c_ecm_g=ecm_instr(lb,Sf+8),ecm_instr(lb,Sg+8)
	win=window(cost,lb);c_mr_f,c_mr_g=c_mr(Sf+8,MONT5),c_mr(Sg+8,MONT5)
	c_band_f=BAND_ROUNDS*ecm_instr(BAND_FBITS,BAND_NBITS[0],BAND_Q)
	c_band_g=BAND_ROUNDS*ecm_instr(BAND_FBITS,BAND_NBITS[1],BAND_Q)
	band_keep=(1-BAND_FMAX[0]*(1-BAND_Q)**BAND_ROUNDS)*(1-BAND_FMAX[1]*(1-BAND_Q)**BAND_ROUNDS)
	pairs_all=(NBAR/2)*2**log2_li(lx2)+((WIDE_SETUP-1)*(NBAR/2)*(2**log2_li(lx2)-2**log2_li(WIDE_FROM))if lx2>WIDE_FROM else 0.0)
	rows=[r for r in rows if r.get('inreg',True)]
	D7=contam_sup(geo)
	for r in rows:
		w=(1-r['qa'])**3
		wsum+=w
		if r['primitive']and not r['rejected']and r['qa']>=Q_CORE:
			rf=retain.at(r['nf'],Sf)['retention'];rg=retain.at(r['ng'],Sg)['retention']
			uf_r,ug_r=(r['nf']+op['alpha_f'])/lb,(r['ng']+op['alpha_g'])/lb
			miss=ECM_FIND**(lam_top(uf_r,lx2/lb)+lam_top(ug_r,lx2/lb))
			y=rho(uf_r)*rho(ug_r)*rf*rg*miss*band_keep
			ysum+=w*y;y2+=w*w*y*y
			nf_r+=w*y*r['nf'];ng_r+=w*y*r['ng']
		if r['rejected']:
			rej_w+=w
		nf_m+=w*r['nf'];ng_m+=w*r['ng']
		S=geo.slice_cells(r['qa']);K=geo.lines(r['qa'])
		p_all=dict(build=(C_BUILD_E*geo.entries_upper(r['qa'])+C_BUILD_SLICE)/S+C_BUILD_SEG*4/2**32,zero=0.0,setup6=pairs_all*C_SETUP6/S,
		thr=C_NB_EVAL*r['evals']/r['bcells']+C_NB_ROW/r['rlen']+C_NB_SLICE/S)
		lp=LX0;line=0.0
		while lp<P_MID:
			lp2=min(lp+0.25,P_MID);dsum=math.log(lp2/lp)
			line+=NBAR*dsum*tab(LINE_TAB,math.log2(S/K)-0.5*(lp+lp2));lp=lp2
		lp=P_MID;rec=0.0;empty=0.0;lg=math.log2(geo.vb*geo.upper)+r['lrs']if RUN_T[0]<=r.get('T',0)<=RUN_T[1]else None
		while lp<lx2:
			lp2=min(lp+0.25,lx2);dsum=math.log(lp2/lp)
			log2h=math.log2(S)-0.5*(lp+lp2)
			if lg is not None and P7<=lp<X7 and lg-0.5*(lp+lp2)>=REC7_TAB[-1][0]:
				dp=(NBAR/2)*(2**log2_li(lp2)-2**log2_li(lp))
				rec+=NBAR*dsum*(step7(lg-0.5*(lp+lp2))+EXTRA7)+dp*C_SETUP6*(2**-lg-1/S)
			elif log2h>=REC_TAB[-1][0]:
				rec+=NBAR*dsum*(WIDE_REC*tab(REC_TAB_W,log2h)if lp>=WIDE_FROM else tab(REC_TAB,log2h))
			else:
				dpairs=(NBAR/2)*(2**log2_li(lp2)-2**log2_li(lp))
				empty+=dpairs*C_COSET_EMPTY/S
			lp=lp2
		rec+=RUN_ALLOW*D7/S # clean-run rule
		hits=NBAR*math.log(lx2/P_MID)
		sf=cost.survivor_upper(r['nf'],Sf,W);sg=cost.survivor_upper(r['ng'],Sg,W)
		q_f=min(1.0,Q_F_SAFETY*rho((r['nf']+op['alpha_f'])/lb)/cost.survivor_upper(r['nf'],Sf,0))
		cbf=BAND_ROUNDS*ecm_instr(BAND_FBITS,min(BAND_NBITS[0],math.ceil(r['nf'])+1),BAND_Q)
		cbg=BAND_ROUNDS*ecm_instr(BAND_FBITS,min(BAND_NBITS[1],math.ceil(r['ng'])+1),BAND_Q)
		ef=win.ecm_share(r['nf'],Sf+W);eg=win.ecm_share(r['ng'],Sg+W);hf=min(1,win.ecm_share(r['nf'],Sf+W,2)+R2_MARGIN);hg=min(1,win.ecm_share(r['ng'],Sg+W,2)+R2_MARGIN)
		c_surv=C_NORMS+C_TDIV+cbf*(1+hf)/2+2*c_mr_f+ef*c_ecm_f+q_f*(cbg*(1+hg)/2+2*c_mr_g+eg*c_ecm_g)
		p_keep=dict(line=line,rec=rec,empty=empty,
		apply=hits*C_APPLY,recover=hits*C_RECOVER,
		route=CELL_ROUTE+REC_ROUTE*hits,
		scan=C_SCAN_CELL+C_SCAN_LEAF*r['leaves']/r['blen'],
		zero=C_ZERO,cof=sf*sg*c_surv)
		per_cell=sum(p_all.values())+sum(p_keep.values())
		per_cell=per_cell*scale_consts if scale_consts!=1.0 else per_cell
		contrib=geo.weight(r['qa'])*per_cell
		work_sum+=contrib;work2+=contrib*contrib
		if return_work:work_list.append(contrib)
	n_a_upper=geo.va/n_total
	population=n_a_upper*sum(geo.weight(r['qa'])for r in rows)
	relcol_instr=n_a_upper*work_sum*CANON_WORK_ALLOWANCE
	per_cell_mean=relcol_instr/population
	rels=n_a_upper*geo.vb*geo.lower*ysum
	yield_se=math.sqrt(max(0.0,y2/n_total-(ysum/n_total)**2)/n_total)/(ysum/n_total)/LN2 if ysum>0 else float('inf')
	nf_m/=wsum;ng_m/=wsum
	log2_U=math.log2((NBAR/2)*2**log2_li(lb)+2**31)
	log2_need=log2_U+math.log2(excess)
	log2_rels=math.log2(rels)if rels>0 else float('-inf')
	feasible=log2_rels>=log2_need+max(MARGIN,MARGIN_SE*yield_se)
	log2_relcol=math.log2(relcol_instr/UNIT)
	out=dict(op={k:(str(v)if isinstance(v,F)else v)for k,v in op.items()},
	log2_population=math.log2(population),per_cell=per_cell_mean,
	log2_relcol=log2_relcol,
	log2_rels=log2_rels,log2_U=log2_U,log2_need=log2_need,
	margin=log2_rels-log2_need,yield_se_log2=yield_se,feasible=feasible,
	rejected_weight=rej_w/wsum,nf_mean=nf_m,ng_mean=ng_m,nf_rel=nf_r/ysum if ysum else 0.0,ng_rel=ng_r/ysum if ysum else 0.0,
	records_per_cell=NBAR*math.log(lx2/P_MID),
	c_surv_log2=math.log2(C_NORMS+C_TDIV+c_band_f+c_band_g+c_ecm_f+c_ecm_g),mean_slice_log2=math.log2(population/(n_a_upper*len(rows))))
	out['sums']=dict(rows=len(rows),ysum=ysum,y2=y2,nfy=nf_r,ngy=ng_r,work=work_sum,work2=work2,w=wsum,nfw=nf_m*wsum,ngw=ng_m*wsum,
	cells=sum(geo.weight(r['qa'])for r in rows))
	if return_work:out['work_contribs']=work_list
	return out
def _filt(): # every filtering run (filt16.b85): FILT_RUNS/FILT_X (seed, N300, max coef) at excess 1.05/1.08; FILT_Y [10y]
	import base64,lzma,json,os
	d=json.loads(lzma.decompress(base64.b85decode(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'filt16.b85')).read().strip())))
	return tuple({int(k):tuple(map(tuple,v))for k,v in d[n].items()}for n in('FILT_RUNS','FILT_X','FILT_Y'))
FILT_RUNS,FILT_X,FILT_Y=_filt()
NU300_H=0.03329 # H3: E[N(300)/U] at excess >= 1.08 (all retained relations, [4])
NU300=0.0333 # cap: abort if N > NU300 U
SD_C,P_CAP=1.0,1e-09 # sd(N(300)/U) <= SD_C/sqrt(U) ([10c]); Chebyshev abort bound, charged
FILT_EXCESS=1.08
FILT_Y_MEANS=(739.1,498.2) # [10y]
W_ROW=301.0
C_ENTRY=21.1
C_ENTRY_K=9.97
K5_RUNS={15:(10.083,10.049,10.197,10.057),16:(9.925,9.938,9.919,9.916,9.958,9.998,9.906,9.978),17:(9.857,9.882,9.855,9.844,9.859,9.881,9.897,9.927)}
C_ROW=2710.0
C_RED9=420.0
UNIT_FRAC_CHARGED=0.86
BW_M=72
BW_N=5
C_PROJ=60.0
N_SM=21
C_FULL=97.0
LOG2_SM=88.0
LOG2_FBC=72.0
LOG2_CANON=91.0
LOG2_MISC=80.0
LOG2_VERIFY=11.0
def linear_algebra(log2_U,excess=None):
	ratio=NU300
	log2_N=log2_U+math.log2(NU300)
	c_prod=W_ROW*C_ENTRY_K+C_ROW+N_SM*C_FULL
	c_kry=W_ROW*C_ENTRY_K+C_ROW+N_SM*C_FULL
	per_n2=c_kry*(BW_N/BW_M+1)+(1/BW_M+1/BW_N)*(c_prod+BW_N*C_FULL+C_RED9)+C_PROJ*(BW_N/BW_M+1)
	per_n1=(BW_M+BW_N)*C_RED9
	log2_la=math.log2((2**(2*log2_N)*per_n2+2**log2_N*per_n1)/UNIT)
	return dict(ratio=ratio,log2_N=log2_N,log2_la=log2_la,log2_gen=generator(log2_N))
def generator(log2_N):
	k=BW_M+BW_N
	log2_L=math.log2((2**log2_N+N_SM)*(1/BW_M+1/BW_N)+2)
	return math.log2(4*k**3*(C_FULL+C_RED9)/UNIT)+math.log2(3)*log2_L+math.log2(log2_L)
def total_log2(log2_relcol,log2_la,log2_desc=DESCENT_CAP,log2_gen=0.0):
	return math.log2(2**log2_relcol+2**log2_la+2**log2_desc+2**log2_gen+2**LOG2_SM+2**LOG2_FBC
	+2**LOG2_CANON+2**LOG2_MISC+2**LOG2_VERIFY+2**LOG2_WIDE_LAT+2**LOG2_PREP7+2**50)
# Cofactor windows (H-EA, proof.md section 3): ECM-path share of the survivor law; toy check ([12w]).
from functools import lru_cache
from threshold_model import coefficients
WIDEN,KMAX,MR_ROUNDS,R2_MARGIN=1.0,4,40,0.08
def c_mr(nbits,mul):
	return MR_ROUNDS*2*nbits*mul
class Window:
	def __init__(self,ret,lb,widen=WIDEN,primes=False):
		r=self.r=ret;st=r.step;n=r.length+1
		self.small=coefficients(r.weights['small'],r.normalizers['small'],r.length)
		self.large=coefficients(r.weights['large_geo'],r.normalizers['large_geo'],r.length)
		w1=[0.0]*n;left=max(r.explicit_cutoff_bits,r.X);cells=math.ceil((r.B-left)/(st/4));width=(r.B-left)/cells
		for i in range(cells):
			b=left+(i+.5)*width;j=round(b/st)
			if j<n:w1[j]+=math.exp(r.tilt*j*st)*width/b
		e=math.exp(-r.normalizers['large_geo'])
		self.comp=self.large if primes else[max(0.0,x-e*y)for x,y in zip(self.large,w1)]
		self.win=[any(r.X*k-widen<=j*st<=lb*k+widen for k in range(2,KMAX+1))for j in range(n)]
	@lru_cache(None)
	def fraction(self,n,stop):
		sv=self.r.pdf['sieve_geo'];sm=self.small;v=[sv[n-t]if t<=n else 0.0 for t in range(stop+1)];tot=ecm=cp=0.0
		for L in range(min(stop,n)+1):
			a=sum(sm[s]*v[L+s]for s in range(stop-L+1));tot+=self.large[L]*a
			if self.win[L]:ecm+=self.comp[L]*a
			cp+=self.comp[L]*a
		return tot,ecm,cp
	def ecm_share(self,norm_bits,slack_bits,i=1):
		st=self.r.step;n=round(norm_bits/st);f=self.fraction(n,min(n,math.floor(slack_bits/st)))
		return min(1.0,f[i]/f[0])if f[0]>0 else 1.0
def composite(n,rng,rounds=MR_ROUNDS):
	d=n-1;s=0
	while d%2==0:d//=2;s+=1
	for _ in range(rounds):
		x=pow(rng.randrange(2,n-1),d,n)
		if x in(1,n-1):continue
		for _ in range(s-1):
			x=x*x%n
			if x==n-1:break
		else:return True
	return False
def toy(n=120,slack=77.5,trials=20000,seed=1,X=16,pmin=16):
	"""The rule on random integers near 2^n (sieve [pmin, 2^X), LB = 69X/64) vs this model there."""
	import random
	from threshold_model import Retention,primes_upto
	lb=X*69/64;R=Retention(n+2,X,prime_min=pmin,explicit_bits=10,step=.25,max_bits=n+3,tilt=.005)
	wm=Window(R,lb);model=wm.ecm_share(n,slack);cm=wm.ecm_share(n,slack,2);law=Window(R,lb,0.0,True).ecm_share(n,slack)
	rng=random.Random(seed);P=math.prod(primes_upto(2**X-1));Ps=math.prod(primes_upto(pmin-1));sv=ec=wn=cc=0
	for _ in range(trials):
		x=rng.randrange(int(2**(n-.125)),int(2**(n+.125)));r=x;g=math.gcd(r,P)
		while g>1:r//=g;g=math.gcd(r,g)
		s=x;g=math.gcd(s,Ps);m=1
		while g>1:s//=g;m*=g;g=math.gcd(s,g)
		if math.log2(m)+math.log2(r)<=slack:
			c=r>3 and composite(r,rng);iw=any(2**(X*k)<r<=2**(lb*k)for k in range(2,KMAX+1));sv+=1;wn+=iw;ec+=iw and c;cc+=c
	return dict(n=n,survivors=sv,ecm=ec,measured=ec/max(sv,1),model=model,in_window=wn/max(sv,1),law=law,comp=cc/max(sv,1),comp_model=cm)

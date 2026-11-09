""
import math
from math import ceil,exp,floor,gcd,log,log2,prod,sqrt
from fractions import Fraction as F
from common12 import P2F,rho,integer_metric,quadratic,regenerate_pair,norms,mul,norm,integer_interval
import ideal_canonical as IC
import embedding_bounds as EB
import filtsim
from threshold_model import Retention
LN2=log(2)
UNIT=20000.0
OP=dict(lb=69.0,lx2=64.0,excess=1.05,slack_f=302,slack_g=233,width_bits=8,L0=64,
theta=F(3,8),alpha_f=0.29,alpha_g=0.75)
MARGIN=0.06
V2E,KAP2=(2.0,2.399985),(0.06,0.05)
def a2(r,s,op=OP):
	return op['alpha_'+'fg'[s]]+V2E[s]+KAP2[s]-r['v2'+'fg'[s]]
def kappa2(u,lb=69.0):
	lo,hi,r0=0.0,1.0,rho(u)
	for _ in range(40):
		m=(lo+hi)/2
		if sum(2.0**-(k+1)*rho(u-(k-1-m)/lb)for k in range(60))>r0:lo=m
		else:hi=m
	return hi
MARGIN_SE=2.0
Q_CORE=0.02
DESCENT_CAP=110.0
NODE_CAP_LOG2=27.3
CANON_WORK_ALLOWANCE=1.01033
WORK_SE_K=2.0
Q_F_SAFETY=100.0
ECM_FIND=0.999
REC_ROUTE=10.01
REC_ROUTE_STORED=10.01
CELL_ROUTE=0.02
ROUTE_CUT=0.48
LINE_TAB=[(18.78,6.06),(13.82,9.04),(9.82,14.26),(5.82,36.0),(1.83,174.93)]
LINE_SE_K=2.0
REC_TAB=[(18.91,23331674/522256),(10.91,3183663/41121),(3.91,1638258/6539),(0.91,3216517/3967),(-1.09,4964184/2901)]
C_COSET_EMPTY=1600.0
C_COS=10.0
C_APPLY=5.76
WIDE_FROM,WIDE_REC,WIDE_SETUP=63.0,1.01,1.0
REC_TAB_W=[(18.91,23449665/522256),(10.91,3216127/41121),(3.91,1745730/6539),(0.91,3721383/3967),(-1.09,6474056/2901)]
LOG2_WIDE_LAT=71.6
REC7_TAB=[(24.99,34.19),(17.99,39.33),(15.99,44.76),(14.99,48.2),(12.99,57.85),(11.99,63.88),(10.99,74.74),(9.99,87.39),(8.99,104.4),(7.99,128.16),(6.99,168.75),(5.99,178.61),(4.99,244.82)]
P7,X7,EXTRA7,RUN_ALLOW,RUN_T=48.0,64.0,2.0,0.00091,(32,1024)
REL_BAND_LOSS=0.009
def hc_sup(t,xm):
	if xm<t[-1][0]:
		return 0.0
	cand=[xm]+[h for h,c in t if h<=xm]
	for(h1,c1),(h2,c2)in zip(t,t[1:]):
		s=(c1-c2)/(h1-h2);x=h2+(-s/log(2)-c2)/s
		if h2<=x<=min(h1,xm):cand.append(x)
	return max(2**x*tab(t,x)for x in cand)
def contam_sup(geo):
	S=geo.slice_cells(0.0);lp=P7;d=0.0
	while lp<X7:
		lp2=min(lp+0.25,X7);m=0.5*(lp+lp2);dp=(NBAR/2)*(2**log2_li(lp2)-2**log2_li(lp));wd=lp>=WIDE_FROM
		d+=NBAR*log(lp2/lp)*2**m*(WIDE_REC*hc_sup(REC_TAB_W,log2(S)-m)if wd else hc_sup(REC_TAB,log2(S)-m))+dp*(C_COSET_EMPTY+C_SETUP6);lp=lp2
	return d
LOG2_PREP7=75.09
def step7(l):
	return next(c for h,c in REC7_TAB if l>=h)
C_SCAN_CELL,C_SCAN_LEAF=1.125,42.88
C_BUILD_E,C_BUILD_SLICE,C_BUILD_SEG=28.57,2.0**21,2.0**17
C_ZERO=0.0
C_SETUP6=406.87+72.0+21.13
C_NB_EVAL=1154.6
C_NB_ROW=535.0
C_NB_SLICE=18500.0
MONT5=367.76
C_NORMS=2.0**22
C_TDIV=2.0**24.5
NBAR=2.0
LX0,P_MID=14.0,48.0
PRIME_MIN,EXPLICIT_BITS=16384,14
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
def log2_fb(lb):
	return log2((NBAR/2)*2**log2_li(lb)+2**31)-1
def log2_li(log2x):
	x=log2x*LN2
	return log2x-log2(x)+log2(1+1/x+2/x/x)
def tab(t,log2h):
	if log2h>=t[0][0]:
		return t[0][1]
	for(h1,c1),(h2,c2)in zip(t,t[1:]):
		if log2h>=h2:
			return c2+(c1-c2)*(log2h-h2)/(h1-h2)
	return t[-1][1]
def skip(r):
	return r.get('p2a')and r.get('p2b')
def coprime_pair(a,b):
	d=gcd(norm(a),norm(b))
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
	p2={k:sum(x*y for x,y in zip(v,P2F))%2==0 for k,v in(('p2a',a),('p2b',b))}
	qa16=round(F(aq,da)/A2*65536)
	nf,ng=norms(a,b)
	nf16=round(16*log2(nf));ng16=round(16*log2(ng));p2.update(v2f=(nf&-nf).bit_length()-1,v2g=(ng&-ng).bit_length()-1)
	prim=coprime_pair(a,b)
	lam,kap=F(cfg['production_region_lam']),F(cfg['production_region_kap'])
	radius=(lam*A2-kap*F(aq,da))*db
	if F(bq)>radius:
		return dict(**p2,seed=sd,index=i,qa=qa16/65536.0,nf=nf16/16.0,ng=ng16/16.0,primitive=prim,rejected=False,nonfinite=False,inreg=False)
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
	return dict(**p2,seed=sd,index=i,qa=qa16/65536.0,nf=nf16/16.0,ng=ng16/16.0,primitive=prim,rejected=isinstance(at[2],dict),
	nonfinite=False,evals=est.counts['block_evaluations'],bcells=sum(h-l+1 for _,l,h in box),leaves=len(leaves),blen=blen,rlen=rlen,inreg=True)
def box_leaves_from(est,rows,start,width=128):
	H,P,L=EB.hull([c for c,_,_ in rows]),list(start),[];ok=lambda bd:bd[0]is not None and bd[0][1]-bd[0][0]<=width
	while P:
		lo,hi=P.pop();act=[k for k,(_,a,b)in enumerate(rows)if a<=hi and b>=lo]
		if not act:continue
		bd=est.block(H,lo,hi)
		if ok(bd):L.append((lo,hi,bd,None))
		elif lo==hi:
			for k in act:bk=est.block(rows[k][0],lo,hi);L.append((lo,hi,bk if ok(bk)else{'rejected':True},k))
		else:m=(lo+hi)//2;P.extend(((m+1,hi),(lo,m)))
	return L
def pool_row_warm(cfg,sd,i,t_lo):
	pool_row(cfg,sd,i);(qa,da,_),(qb,db,_),dk=_CTX['qa'],_CTX['qb'],_CTX['disks'];A2=cfg['production_A_squared']
	lam,kap=F(cfg['production_region_lam']),F(cfg['production_region_kap']);a,b,aq,bq=regenerate_pair(cfg,sd,i);blo=b[0]//64*64
	rad=lambda x:(lam*A2-kap*F(quadratic(qa,x),da))*db;box_of=lambda x:box_rows(qb,rad(x),b,blo)
	def cert(x,bx,st):
		e=EB.NormBounds(x,dk);rw=[(e.row(list(s)),l,h)for s,l,h in bx];return e,(EB.box_leaves(e,rw,blo,blo+63)if st is None else box_leaves_from(e,rw,st))
	st=None
	for t in range(t_lo,0):
		x=[a[0]+t]+a[1:];bx=box_of(x)
		if bx:st=[(l[0],l[1])for l in cert(x,bx,None)[1]];break
	box=box_of(a);lo0,hi0=min(r[1]for r in box),max(r[2]for r in box)
	if st is not None:
		st=sorted({(max(l,lo0),min(h,hi0))for l,h in st if max(l,lo0)<=min(h,hi0)});g,c=[],lo0
		for l,h in st:
			if l>c:g.append((c,l-1))
			c=max(c,h+1)
		st+=g+([(c,hi0)]if c<=hi0 else[])
	e,lv=cert(a,box,st);iv=integer_interval(qb[0][0],sum(qb[0][j]*b[j]for j in range(1,6)),sum(qb[j][k]*b[j]*b[k]for j in range(1,6)for k in range(1,6))-rad(a))
	lo,hi=max(iv.start,blo),min(iv.stop-1,blo+63)
	return e.counts['block_evaluations'],len([l for l in lv if l[0]<=hi and l[1]>=lo and l[3]in(None,0)])
def pool_totals(s,n_eff,n_draws,geo):
	n,ys,y2,w,w2=s['rows'],s['ysum'],s['y2'],s['work']-(ROUTE_CUT+ROUTE_REC_CUT)*s['cells'],s['work2']
	mw=w/n
	cond=sqrt(max(0.0,(w2-n*mw*mw)/(n-1))/n)/mw
	return dict(log2_rels=log2(geo.va/n_eff*geo.vb*geo.lower*ys),
	yield_se_log2=sqrt(max(0.0,y2/n_eff-(ys/n_eff)**2)/n_eff)/(ys/n_eff)/LN2,
	nf_rel=s['nfy']/ys,ng_rel=s['ngy']/ys,nf_mean=s['nfw']/s['w'],ng_mean=s['ngw']/s['w'],
	log2_population=log2(geo.va/n_eff*s['cells']),per_cell=w/s['cells'],log2_est=log2(geo.va/n_eff*w/UNIT),
	cond_se=cond,work_se=sqrt(cond**2+(1-n/n_draws)/n))
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
		self.lam,self.kap=float(F(cfg['production_region_lam'])),float(F(cfg['production_region_kap']))
		self.qa,self.da,self.ia=integer_metric(cfg['matrix_a'])
		self.qb,self.db,self.ib=integer_metric(cfg['matrix_b'])
		from ideal_canonical import determinant
		det_a=F(determinant(self.qa),self.da**6);det_b=F(determinant(self.qb),self.db**6)
		v6=math.pi**3/6
		self.va=v6*(float(theta)*self.A2)**3/(2*sqrt(float(det_a)))
		self.vb=v6*self.A2**3/(2*sqrt(float(det_b)))
		self.lower,self.upper=float(lower),float(upper)
		self.rb=[sqrt(float(self.ib[j][j])*self.A2*self.db)for j in range(6)]
		self.inv55=float(self.ib[5][5])*self.db
	def rad(self,qa_frac):
		return self.lam-self.kap*qa_frac
	def slice_cells(self,qa_frac):
		return self.vb*self.rad(qa_frac)**3*self.upper
	def weight(self,qa_frac):
		return self.vb*(1-qa_frac)**3*self.upper
	def lines(self,qa_frac):
		return max(1.0,floor(sqrt(self.rad(qa_frac)*self.A2*self.inv55)))
	def entries_upper(self,qa_frac):
		s=sqrt(self.rad(qa_frac));L4=1.0
		for j in range(1,5):
			L4*=2*floor(self.rb[j]*s)+3
		return(floor(sqrt(self.rad(qa_frac)*self.A2*self.inv55))+2)*L4
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
	band_margin=(1-BAND_FMAX[0]*(1-BAND_Q)**BAND_ROUNDS)
	pairs_all=(NBAR/2)*2**log2_li(lx2)+((WIDE_SETUP-1)*(NBAR/2)*(2**log2_li(lx2)-2**log2_li(WIDE_FROM))if lx2>WIDE_FROM else 0.0)
	rows=[r for r in rows if r.get('inreg',True)]
	D7=contam_sup(geo)
	for r in rows:
		w=(1-r['qa'])**3
		wsum+=w
		if r['primitive']and not r['rejected']and r['qa']>=Q_CORE:
			rf=retain.at(r['nf'],Sf)['retention'];rg=1.0
			uf_r,ug_r=(r['nf']+a2(r,0,op))/lb,(r['ng']+a2(r,1,op))/lb
			miss=ECM_FIND**lam_top(uf_r,lx2/lb)
			y=rho(uf_r)*rho(ug_r)*rf*rg*miss*band_margin
			ysum+=w*y;y2+=w*w*y*y
			nf_r+=w*y*r['nf'];ng_r+=w*y*r['ng']
		if r['rejected']:
			rej_w+=w
		nf_m+=w*r['nf'];ng_m+=w*r['ng']
		if skip(r):continue
		fk=.5 if r.get('p2a')else 1.;S=geo.slice_cells(r['qa'])*fk;K=geo.lines(r['qa'])
		p_all=dict(build=(C_BUILD_E*geo.entries_upper(r['qa'])+C_BUILD_SLICE)/S+C_BUILD_SEG*4/2**32,zero=0.0,setup6=pairs_all*(C_SETUP6+C_COS*(fk<1))/S,
		thr=(C_NB_EVAL*(lambda w:w+RUN_ALLOW*max(0,r['evals']-w))(r.get('evals_w',r['evals']))/r['bcells']+C_NB_ROW/r['rlen'])/fk+C_NB_SLICE/S)
		lp=LX0;line=0.0
		while lp<P_MID:
			lp2=min(lp+0.25,P_MID);dsum=log(lp2/lp)
			line+=NBAR*dsum*tab(LINE_TAB,log2(S/K)-0.5*(lp+lp2));lp=lp2
		lp=P_MID;rec=0.0;empty=0.0;u7=0.0;lg=log2(geo.vb*geo.upper*fk/2)+r['lrs']if RUN_T[0]<=r.get('T',0)<=RUN_T[1]else None
		while lp<lx2:
			lp2=min(lp+0.25,lx2);dsum=log(lp2/lp)
			log2h=log2(S)-0.5*(lp+lp2)
			if lg is not None and P7<=lp<X7 and lg-0.5*(lp+lp2)>=REC7_TAB[-1][0]:
				dp=(NBAR/2)*(2**log2_li(lp2)-2**log2_li(lp))
				rec+=NBAR*dsum*(step7(lg-0.5*(lp+lp2))+EXTRA7);u7+=dp*((C_SETUP6+C_COS)*2**-lg-(C_SETUP6+C_COS*(fk<1))/S)
			elif log2h>=REC_TAB[-1][0]:
				rec+=NBAR*dsum*(WIDE_REC*tab(REC_TAB_W,log2h)if lp>=WIDE_FROM else tab(REC_TAB,log2h))
			else:
				dpairs=(NBAR/2)*(2**log2_li(lp2)-2**log2_li(lp))
				empty+=dpairs*C_COSET_EMPTY/S
			lp=lp2
		rec+=RUN_ALLOW*D7/S;p_all['setup6']=max(0.0,p_all['setup6']+u7)
		hits=NBAR*log(lx2/P_MID)
		sf=cost.survivor_upper(r['nf'],Sf,W)
		q_g=min(1.0,Q_F_SAFETY*rho((r['ng']+a2(r,1,op))/lb))
		ef=win.ecm_share(r['nf'],Sf+W)
		c_surv=C_NORMS+C_GBATCH+q_g*(C_TDIV+C_BATCH+C_BATCH2+c_mr_f+ef*c_ecm_f)
		line*=SIEVE_SIDES;rec*=SIEVE_SIDES;empty*=SIEVE_SIDES;hits=hits*SIEVE_SIDES
		p_keep=dict(line=line,rec=rec,empty=empty,
		apply=hits*C_APPLY,
		route=CELL_ROUTE+REC_ROUTE*hits,
		scan=C_SCAN_CELL+C_SCAN_LEAF*max(r['leaves'],r.get('leaves_w',0))/r['blen']/fk,
		zero=C_ZERO,cof=sf*c_surv)
		per_cell=sum(p_all.values())+sum(p_keep.values())
		per_cell=per_cell*scale_consts if scale_consts!=1.0 else per_cell
		contrib=geo.weight(r['qa'])*per_cell
		work_sum+=contrib;work2+=contrib*contrib
		if return_work:work_list.append(contrib)
	n_a_upper=geo.va/n_total
	population=n_a_upper*sum(geo.weight(r['qa'])for r in rows if not skip(r))
	relcol_instr=n_a_upper*work_sum*CANON_WORK_ALLOWANCE
	per_cell_mean=relcol_instr/population
	rels=n_a_upper*geo.vb*geo.lower*ysum
	yield_se=sqrt(max(0.0,y2/n_total-(ysum/n_total)**2)/n_total)/(ysum/n_total)/LN2 if ysum>0 else float('inf')
	nf_m/=wsum;ng_m/=wsum
	log2_U=log2_fb(lb)
	log2_need=log2_U+log2(excess)
	log2_rels=log2(rels)if rels>0 else float('-inf')
	feasible=log2_rels>=log2_need+max(MARGIN,MARGIN_SE*yield_se)
	log2_relcol=log2(relcol_instr/UNIT)
	out=dict(op={k:(str(v)if isinstance(v,F)else v)for k,v in op.items()},
	log2_population=log2(population),per_cell=per_cell_mean,
	log2_relcol=log2_relcol,
	log2_rels=log2_rels,log2_U=log2_U,log2_need=log2_need,
	margin=log2_rels-log2_need,yield_se_log2=yield_se,feasible=feasible,
	rejected_weight=rej_w/wsum,nf_mean=nf_m,ng_mean=ng_m,nf_rel=nf_r/ysum if ysum else 0.0,ng_rel=ng_r/ysum if ysum else 0.0,
	records_per_cell=NBAR*log(lx2/P_MID),
	c_surv_log2=log2(C_NORMS+C_GBATCH+C_TDIV+C_BATCH+c_ecm_f),mean_slice_log2=log2(population/(n_a_upper*len(rows))))
	out['sums']=dict(rows=len(rows),ysum=ysum,y2=y2,nfy=nf_r,ngy=ng_r,work=work_sum,work2=work2,w=wsum,nfw=nf_m*wsum,ngw=ng_m*wsum,
	cells=sum(geo.weight(r['qa'])for r in rows if not skip(r)))
	if return_work:out['work_contribs']=work_list
	return out
def _filt():
	import base64,lzma,json,os
	d=json.loads(lzma.decompress(base64.b85decode(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'filt16.b85')).read().strip())))
	return tuple({int(k):tuple(map(tuple,v))for k,v in d[n].items()}for n in('FILT_RUNS','FILT_X','FILT_Y'))
FILT_RUNS,FILT_X,FILT_Y=_filt()
NU_H=0.02687
NU=0.0274
SD_C,P_CAP=1.0,0.0
FILT_EXCESS=1.095
FILT_Y_MEANS=(740.0,498.4)
W_ROW=401.0
C_ENTRY_K=8.38
K5_RUNS={15:(8.425,8.443,8.385,8.382),16:(8.365,8.373,8.348,8.34,8.364,8.379,8.365,8.366),17:(8.316,8.327,8.317,8.32,8.323,8.334,8.323,8.314)}
C_ROW=2710.0
C_RED9=266.0
BW_M=512
BW_N=25
C_PROJ=60.0
N_SM=10
C_FULL=97.0
C_CHI=198.0
LOG2_SM=88.0
LOG2_FBC=73.0
LOG2_CANON=92.0
LOG2_MISC=80.0
LOG2_VERIFY=11.0
LOG2_GTREE,LOG2_GSPLIT=74.0,98.0
def linear_algebra(log2_U,excess=None):
	ratio=NU
	log2_N=log2_U+log2(NU)
	c_prod=W_ROW*C_ENTRY_K+C_ROW+N_SM*C_FULL+C_CHI
	c_kry=c_prod
	per_n2=c_kry*(BW_N/BW_M+1)+(1/BW_M+1/BW_N)*(c_prod+BW_N*C_FULL+C_RED9)+C_PROJ*(BW_N/BW_M+1)
	per_n1=(BW_M+BW_N)*C_RED9
	log2_la=log2((2**(2*log2_N)*per_n2+2**log2_N*per_n1)/UNIT)
	return dict(ratio=ratio,log2_N=log2_N,log2_la=log2_la,log2_gen=generator(log2_N))
def generator(log2_N):
	k=BW_M+BW_N
	log2_L=log2((2**log2_N+N_SM)*(1/BW_M+1/BW_N)+2)
	L=2**log2_L;return log2(4*k**3*(fft_F(18*L)+L*(18*FFT_COPY+4*(C_RED9+C_FULL)))/UNIT)+log2(log2_L)
def total_log2(log2_relcol,log2_la,log2_desc=DESCENT_CAP,log2_gen=0.0):
	return log2(2**log2_relcol+2**log2_la+2**log2_desc+2**log2_gen+2**LOG2_SM+2**LOG2_FBC
	+2**LOG2_CANON+2**LOG2_MISC+2**LOG2_VERIFY+2**LOG2_WIDE_LAT+2**LOG2_PREP7+2**LOG2_SPLIT+2**LOG2_GTREE+2**LOG2_GSPLIT+2**50)
from functools import lru_cache
from threshold_model import coefficients
WIDEN,KMAX,MR_ROUNDS=1.0,4,40
SIEVE_SIDES=0.5
ROUTE_REC_CUT=(REC_ROUTE_STORED-REC_ROUTE)*NBAR*log(OP['lx2']/P_MID)*SIEVE_SIDES
C_GBATCH=2.0**25
C_BATCH=2.0**24
C_BATCH2=2.0**25.25
def srt_toy(k,seed,dxs=(-2,2),G=None):
	import random
	from threshold_model import primes_upto
	ps=[p for p in primes_upto(2**21)if p>2**20];P=ptree(ps)[-1][0];g=random.Random(seed);zs=[]
	for _ in range(2**k):
		z=g.getrandbits(832)|1
		for _ in range(g.randrange(3)):z*=g.choice(ps)**g.randrange(1,3)
		zs.append(z)
	R=[srt(P,zs,G,dx)for dx in dxs]
	return sum(a!=P%z for r in R for a,z in zip(r[0],zs)),max(r[1]for r in R)
def ptree(zs):
	T=[list(zs)]
	while len(T[-1])>1:
		a=T[-1];T.append([a[i]*a[i+1]if i+1<len(a)else a[i]for i in range(0,len(a),2)])
	return T
def srt(P,zs,G=None,dx=0):
	T=ptree(zs)
	G=G if G is not None else len(T)+1;w=lambda Z:64*(-(-Z.bit_length()//64)+G);Z=T[-1][0];b=w(Z)
	s=b+P.bit_length()+8;x=(1<<s)//Z+dx
	q=(P*x)>>(s-b);r=(P<<b)-q*Z;fix=0
	while r<0:q-=1;r+=Z;fix+=1
	while r>=Z:q+=1;r-=Z;fix+=1
	Y=[q%(1<<b)]
	for l in range(len(T)-2,-1,-1):
		C=T[l];Yn=[]
		for i,z in enumerate(C):
			v=i//2;Zv=T[l+1][v];sib=C[i^1]if(i^1)<len(C)else 1;bv,bc=w(Zv),w(z)
			Yn.append(((Y[v]*sib)%(1<<bv))>>(bv-bc))
		Y=Yn
	out=[]
	for z,y in zip(zs,Y):
		b=w(z);r=(z*y+(1<<b)-1)>>b;out.append(0 if r==z else r)
	return out,fix

FFT_A,FFT_B,FFT_ROOT,FFT_COPY,FFT_LEAF,FFT_XFER,FFT_LOCAL=8.5,1200.0,6.0,3.0,2.0**15,504.0,2**26
def fft_F(N):return FFT_A*N*log2(max(N,2))**2+FFT_B+(FFT_XFER*N if N>FFT_LOCAL else 0.0)
def gbatch_cost(pbits=1.000028*2**69/LN2,leaf=9,L=None,c_leaf=FFT_LEAF):
	SP=ceil(pbits/64);L=L if L is not None else ceil(log2((SP+1)/leaf));G=L+2;n=lambda j:leaf*2**j
	PT=sum(2**(L-j)*(fft_F(2*n(j-1))+FFT_COPY*2*n(j-1))for j in range(1,L+1))
	RT=sum(2**(L-j+1)*(fft_F(n(j)+G+n(j-1))+FFT_COPY*(n(j-1)+G))for j in range(1,L+1))
	return(PT+RT+FFT_ROOT*fft_F(2*(max(SP,n(L))+G)))/2**L+c_leaf
LOG2_SPLIT=100.62
SPLIT_LIMB,SPLIT_SIDE,SPLIT_LOOP=6,22,4
SPLIT_DIV,SPLIT_IDEAL=2.0**10,2.0**26
def band_split(log2_rels_est,yield_se,rmin_f,rmin_g):
	lf,lg=[-(-b//64)for b in BAND_NBITS];np_=1.25506*2**48/(48*LN2)
	c_rel=np_*(SPLIT_LIMB*(lf+lg)+2*SPLIT_SIDE+SPLIT_LOOP)+sum(BAND_FMAX)*SPLIT_DIV+SPLIT_IDEAL
	bm=(1-BAND_FMAX[0]*(1-BAND_Q)**BAND_ROUNDS)*(1-BAND_FMAX[1]*(1-BAND_Q)**BAND_ROUNDS)
	lc=log2_rels_est+log2(1.032972622/0.966347073)-(lf+lg)*log2(ECM_FIND)-log2(bm*rmin_f*rmin_g)+2*yield_se
	return dict(log2_primes=log2(np_),c_rel_log2=log2(c_rel),log2_count=lc,log2_units=lc+log2(c_rel/UNIT))
def rmin(ret,slack,top=850):
	v=[]
	for k in range(4,4*top+1):
		try:v.append(ret.at(k/4,slack)['retention'])
		except ArithmeticError:pass
	return min(v)
def c_mr(nbits,mul):
	return MR_ROUNDS*2*nbits*mul
class Window:
	def __init__(self,ret,lb,widen=WIDEN,primes=False):
		r=self.r=ret;st=r.step;n=r.length+1
		self.small=coefficients(r.weights['small'],r.normalizers['small'],r.length)
		self.large=coefficients(r.weights['large_geo'],r.normalizers['large_geo'],r.length)
		w1=[0.0]*n;left=max(r.explicit_cutoff_bits,r.X);cells=ceil((r.B-left)/(st/4));width=(r.B-left)/cells
		for i in range(cells):
			b=left+(i+.5)*width;j=round(b/st)
			if j<n:w1[j]+=exp(r.tilt*j*st)*width/b
		e=exp(-r.normalizers['large_geo'])
		self.comp=self.large if primes else[max(0.0,x-e*y)for x,y in zip(self.large,w1)]
		self.win=[any(r.X*k-widen<=j*st<=lb*k+widen for k in range(2,KMAX+1))for j in range(n)]
	@lru_cache(None)
	def fraction(self,n,stop):
		sv=self.r.pdf['sieve_geo'];sm=self.small;v=[sv[n-t]if t<=n else 0.0 for t in range(stop+1)];tot=ecm=0.0
		for L in range(min(stop,n)+1):
			a=sum(sm[s]*v[L+s]for s in range(stop-L+1));tot+=self.large[L]*a
			if self.win[L]:ecm+=self.comp[L]*a
		return tot,ecm
	def ecm_share(self,norm_bits,slack_bits,i=1):
		st=self.r.step;n=round(norm_bits/st);f=self.fraction(n,min(n,floor(slack_bits/st)))
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
	""
	import random
	from threshold_model import Retention,primes_upto
	lb=X*69/64;R=Retention(n+2,X,prime_min=pmin,explicit_bits=10,step=.25,max_bits=n+3,tilt=.005)
	wm=Window(R,lb);model=wm.ecm_share(n,slack);law=Window(R,lb,0.0,True).ecm_share(n,slack)
	rng=random.Random(seed);P=prod(primes_upto(2**X-1));Ps=prod(primes_upto(pmin-1));sv=ec=wn=cc=0
	for _ in range(trials):
		x=rng.randrange(int(2**(n-.125)),int(2**(n+.125)));r=x;g=gcd(r,P)
		while g>1:r//=g;g=gcd(r,g)
		s=x;g=gcd(s,Ps);m=1
		while g>1:s//=g;m*=g;g=gcd(s,g)
		if log2(m)+log2(r)<=slack:
			c=r>3 and composite(r,rng);iw=any(2**(X*k)<r<=2**(lb*k)for k in range(2,KMAX+1));sv+=1;wn+=iw;ec+=iw and c;cc+=c
	return dict(n=n,survivors=sv,ecm=ec,measured=ec/max(sv,1),model=model,in_window=wn/max(sv,1),law=law)

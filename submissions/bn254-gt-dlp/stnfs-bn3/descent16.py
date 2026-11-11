import math,functools;from math import exp,factorial,floor,lgamma,log,log1p,log2,sqrt;from decimal import Decimal,localcontext;LN2=log(2);UNIT=20000.;KMAX=64;RSTEP=1./1024
def _rho_coefficients(kmax=KMAX,nt=600):
	with localcontext()as ctx:
		ctx.prec=200;zero,one=Decimal(0),Decimal(1);coefs=[None,[one]+[zero]*(nt-1)]
		for k in range(2,kmax+1):
			previous,c=coefs[k-1],[zero]*nt
			for m in range(nt-1):c[m+1]=(previous[m]+m*c[m])/(k*(m+1))
			c[0]=previous[0]-sum(c[1:]);coefs.append(c)
		out=[None]
		for c in coefs[1:]:
			big=max(abs(x)for x in c)
			while len(c)>1 and abs(c[-1])<big*Decimal('1e-25'):c=c[:-1]
			out.append([float(x)for x in c])
		return out
RHO_COEFS=_rho_coefficients()
def rho_series(u):
	if u<=1:return 1.
	k=math.ceil(u)
	if k>KMAX:return.0
	x,value=k-u,.0
	for c in reversed(RHO_COEFS[k]):value=value*x+c
	return max(value,.0)
RHO_TAB=[rho_series(i*RSTEP)for i in range(int(KMAX/RSTEP)+2)]
def rho_lo(u):
	if u<=1:return 1.
	i=math.ceil(u/RSTEP);return RHO_TAB[i]if i<len(RHO_TAB)else.0
def rho_hi(u):
	if u<=1:return 1.
	i=floor(u/RSTEP)
	if i+1>=len(RHO_TAB):return 1./factorial(floor(u))
	return RHO_TAB[i]
def kappa(u,t):return rho_series(u-t)/rho_series(u)if rho_series(u)>0 else.0
@functools.lru_cache(None)
def miss_side(u4,log2_rels,u_main,orbit=4,lb=56.):
	u=u4/4.;tot=.0
	for b10 in range(300,int(10*lb)):b=b10/10.;lam=2.**(log2_rels-b)*orbit*kappa(u_main,b/lb);tot+=kappa(u,b/lb)*log((b+.1)/b)*exp(-lam)*(1+lam+lam*lam/2)
	return tot
def missing(Pm,nf,ng):f=lambda m,um:miss_side(math.ceil(4*m/Pm.LX),Pm.LOG2_RELS,um,Pm.ORBIT,Pm.LX);return(f(nf,Pm.MAIN_UF)if nf else.0)+(f(ng,Pm.MAIN_UG)if ng else.0)
KA_MAX=16
def _conv(a,b):
	out=[.0]*(len(a)+len(b)-1);nz=[j for j,y in enumerate(b)if y]
	if nz:
		j0,j1=nz[0],nz[-1]+1;bb=b[j0:j1]
		for i,x in enumerate(a):
			if x:out[i+j0:i+j1]=[o+x*y for o,y in zip(out[i+j0:i+j1],bb)]
	return out
@functools.lru_cache(None)
def window_rows(n,base,iA,mesh,kmax=KA_MAX):
	h=(n-base)/mesh;wA=[log1p(h/(base+i*h))if i<iA else.0 for i in range(mesh)];wB=[log1p(h/(base+i*h))if i>=iA else.0 for i in range(mesh)];PA=[[1.]]
	for j in range(1,kmax+1):PA.append([x/j for x in _conv(PA[-1],wA)])
	PB=[[1.],wB];rows={}
	for jA in range(kmax+1):
		for jB in range(2):rows[jA,jB]=PA[jA]if jB==0 else _conv(PA[jA],PB[1])
	return rows,h
@functools.lru_cache(None)
def ss_cum2(m,base,n,iA,mesh,alpha,trunc,kmax=KA_MAX):
	rows,h=window_rows(n,base,iA,mesh,kmax);out=[[.0,.0]for _ in range(kmax+1)]
	for(jA,jB),row in rows.items():
		j=jA+jB;imax=len(row)-1;ifeas=min(imax,int(floor((m-j*h-j*base)/h)))if h>0 else imax
		if ifeas<0:continue
		imin=0 if not trunc else max(0,ifeas-int(trunc*base/h)-1);tab,ntab,ceil=RHO_TAB,len(RHO_TAB),math.ceil;k0=(m-j*base+alpha)/(base*RSTEP);hk=h/(base*RSTEP);s=.0
		for t in[row[i]*(1. if(k:=ceil(k0-i*hk))<=0 else tab[k]if k<ntab else.0)for i in range(imin,ifeas+1)if row[i]]:s+=t
		out[jA][jB]=s
	cum=[[.0,.0]for _ in range(kmax+1)]
	for jB in range(2):
		acc=.0
		for jA in range(kmax+1):acc+=out[jA][jB];cum[jA][jB]=acc
	return tuple(tuple(x)for x in cum)
def ss_lower2(m,base,n,iA,kA,kB,mesh,alpha,trunc=3.):
	m=math.ceil(m)
	if n<=base:return rho_lo((m+alpha)/base)
	cum=ss_cum2(m,base,n,iA,mesh,alpha,trunc,kA);return cum[kA][0]+(cum[kA][1]if kB else.0)
def surv_upper(m,base,T):
	if T<=0:return.0
	if T>=m:return 1.
	return min(1.,max(1.,T/base)*rho_hi((floor(m)-T)/base))
MONT5=367.76
def ecm_instr(fbits,nbits,success):import ecm16;return ecm16.ecm_instr(fbits,nbits,success)
def log2li(l):x=l*LN2;return l-log2(x)+log2(1+1/x+2/x/x)
CELLS={4:13315520,6:1573359464,8:47982017060,10:687193400292,12:6073556677752,14:38414493101628,16:190058478539564,20:2754504316057176,24:24505531366261200,28:155610653675916932,32:771909602591773000,40:11221331971811657968,48:99995402971660463340,56:635616373364584622852,64:3155035383949694198996,80:45900097950907999286464,96:409193506988359583551012,112:2601677688965738307448548,128:12916156661164168796080328,160:187942845067274337246511324,192:1675659512501707896389090868,224:10654623198171176305341831664,256:52897551248440914796084637156,288:217401149762189752879411886484,320:769748406668707279529221545272,384:6863095027560050329435265332068}
class P:LB=56.;LX=56.;PMIN=8.;SLACK=96.;DF=48;ALPHA_F=.0;ALPHA_G=3.4;LOG2U6=6*log2(4965661367192848881);QBINS=[(1/6,13.6,14.6),(1/6,23.7,14.6),(1/6,31.5,14.6),(1/6,38.7,14.6),(1/6,47.6,14.6),(1/6,98.3,14.6)];BOOT_BINS=[(0.5,-178.78),(0.4,-175.78),(0.1,-171.93)];QBINS_LO=[(-92.7,-53.4),(-6.5,-53.4),(3.6,-53.4),(11.4,-53.4),(18.6,-53.4),(27.5,-53.4)];QBINS_SUB=[[-1.9,5.4],[5.4,4.7],[10.0,6.0],[13.6,6.4],[16.6,6.7],[19.2,7.6],[21.6,8.5],[23.7,8.3],[25.8,9.4],[27.8,10.0],[29.6,9.4],[31.5,10.2],[33.2,10.4],[35.0,10.2],[36.8,11.3],[38.7,11.3],[40.7,11.5],[42.8,12.6],[45.0,11.2],[47.6,12.2],[50.4,13.6],[53.9,12.7],[59.2,11.6],[98.3,14.6]];BOOT_SUB=[-189.12,-186.97,-185.57,-184.44,-183.37,-182.43,-181.41,-180.27,-178.78,-171.93];PRCAP_NODE=4.;PRCAP_BOOT=1.5;LOG2_RELS=49.85;MAIN_UF=417./56.;MAIN_UG=471./56.;ORBIT=4;THIN_MAX=72;UDUP=.65;MIN_RELS=18.;PSI_MAX=.8;PHI_FLOOR=.02;DISC=.01;BOOT_RELS=18;BALL_PAIR=2.**24;NB_R1=.5;NB_R2=1.;BOUND_CELL=2400.;HITS_PER_CELL=2*log(56./8.);SIEVE_HIT=1100.;SCAN=10.;ZERO=2.;PAIR_SETUP=2.**25.1;LLL_NODE=2.**36;NORMS=2.**26;TDIV=2.**24.5;IDEAL=2.**34;PRIME=2.**24;ECM_SUCCESS=.999;MESH=96;TRUNC=.0;LOG2P6=6*log2(21888242871839275222246405745257275088696311157297823662689037894645226208583);LMAX=262;BIN=5
def pairs_count(Pm):return 2.*2**log2li(Pm.LX)
def node_cost(Pm,l_lo,l_hi,side,n,iA,kA,kB,R,thin=None):
	pairs=pairs_count(Pm);cfmax=max(b[1]for b in Pm.QBINS);cgmax=max(b[2]for b in Pm.QBINS);S=log2(CELLS[R]);c=(l_hi+S)/12;h=(n-Pm.LX)/Pm.MESH if n>Pm.LX else.0;nA=Pm.LX+iA*h if n>Pm.LX else Pm.LX;T=kA*nA+kB*n+Pm.SLACK;c_lo=(l_lo+S)/12;pr=.0;sv=.0
	for(w,cf,cg),(cf_lo,cg_lo)in zip(Pm.QBINS,Pm.QBINS_LO):nf=Pm.DF*c+cf-(l_lo if side=='f'else 0);ng=12*c+Pm.LOG2U6+cg-(l_lo if side=='g'else 0);pf=ss_lower2(nf,Pm.LX,n,iA,kA,kB,Pm.MESH,Pm.ALPHA_F,Pm.TRUNC);pg=ss_lower2(ng,Pm.LX,n,iA,kA,kB,Pm.MESH,Pm.ALPHA_G,Pm.TRUNC);pr+=w*pf*pg;nf_lo=Pm.DF*c_lo+cf_lo-(l_hi if side=='f'else 0);ng_lo=12*c_lo+Pm.LOG2U6+cg_lo-(l_hi if side=='g'else 0);sv+=w*surv_upper(nf_lo,Pm.LX,T)*surv_upper(ng_lo,Pm.LX,T)
	pri=.0
	for cf,cg in Pm.QBINS_SUB:nf=Pm.DF*c+cf-(l_lo if side=='f'else 0);ng=12*c+Pm.LOG2U6+cg-(l_lo if side=='g'else 0);pri+=ss_lower2(nf,Pm.LX,n,iA,kA,kB,Pm.MESH,Pm.ALPHA_F,Pm.TRUNC)*ss_lower2(ng,Pm.LX,n,iA,kA,kB,Pm.MESH,Pm.ALPHA_G,Pm.TRUNC)
	pr=min(pri/24,Pm.PRCAP_NODE*pr)
	if pr<=0:return None
	nf_max=Pm.DF*c+cfmax-(l_lo if side=='f'else 0);ng_max=12*c+Pm.LOG2U6+cgmax-(l_lo if side=='g'else 0);det=1.-(1-Pm.ECM_SUCCESS)*((nf_max+ng_max)/20.+2*(kA+kB));miss=missing(Pm,nf_max,ng_max)
	if thin is None or l_lo<=Pm.THIN_MAX:thin=thin or l_lo<=Pm.THIN_MAX
	mu=CELLS[R]*pr*Pm.UDUP*det*(exp(-miss)if thin else 1.)
	if mu<Pm.MIN_RELS:return None
	surv_cost=Pm.NORMS+Pm.TDIV+ecm_instr(Pm.LX,nf_max,Pm.ECM_SUCCESS)+ecm_instr(Pm.LX,ng_max,Pm.ECM_SUCCESS)+2*(ecm_instr(nA,T,Pm.ECM_SUCCESS)+Pm.PRIME)+Pm.IDEAL
	survivors=CELLS[R]*sv+CELLS[R]*pr;sieve=CELLS[R]*(Pm.HITS_PER_CELL*Pm.SIEVE_HIT+Pm.SCAN+Pm.ZERO+Pm.BOUND_CELL)+pairs*Pm.PAIR_SETUP;trials=CELLS[R]*surv_cost
	if trials<sieve+survivors*surv_cost:own=Pm.LLL_NODE+trials;how='ecm';survivors=CELLS[R];ret,loss_f,loss_g=1.,.0,.0
	else:
		own=Pm.LLL_NODE+sieve+survivors*surv_cost;how='sieve';import retain16;ret,loss_f,loss_g=retain16.retention(nf_max,ng_max,n,iA,kA,kB,Pm.SLACK);mu*=ret
		if mu<Pm.MIN_RELS:return None
	return dict(own=own,n=n,iA=iA,nA=nA,kA=kA,kB=kB,R=R,S=S,mu=mu,how=how,miss=miss,ret=ret,loss_f=loss_f,loss_g=loss_g,miss_children=.0 if thin else miss,thin=bool(thin),pr_log2=log2(pr),surv_log2=log2(survivors),own_log2=log2(own),nf=nf_max,ng=ng_max,det=det,T=T,sieve_log2=log2(sieve),survcost_log2=log2(survivors*surv_cost))
def bins_of(Pm):
	bins=[(int(Pm.LB),int(Pm.LB)+Pm.BIN+1)];l=int(Pm.LB)+Pm.BIN+1
	while l<Pm.LMAX:bins.append((l,l+Pm.BIN));l+=Pm.BIN
	return bins
_PLAN={56:((56,0,0,0,80,1),(56,0,0,0,80,1)),62:((56,0,0,0,80,1),(62,72,1,0,56,1)),67:((67,48,1,0,64,1),(67,48,1,0,64,1)),72:((64,72,1,0,80,1),(64,72,1,0,80,1)),77:((69,72,1,0,80,1),(77,48,1,0,80,1)),82:((66,48,1,1,48,1),(82,72,1,0,80,1)),87:((71,24,1,1,56,1),(71,24,1,1,56,1)),92:((76,24,1,1,56,1),(76,24,1,1,56,1)),97:((73,48,1,1,64,1),(73,24,1,1,80,1)),102:((70,24,1,1,80,0),(78,24,1,1,80,1)),107:((75,24,1,1,80,0),(75,48,1,1,80,0)),112:((80,24,1,1,80,0),(80,48,1,1,80,0)),117:((85,24,1,1,80,0),(85,24,1,1,96,0)),122:((90,48,1,1,80,0),(106,48,2,0,96,0)),127:((95,48,1,1,80,0),(95,48,1,1,96,0)),132:((100,24,1,1,96,0),(88,48,2,1,80,0)),137:((105,48,1,1,80,0),(93,48,2,1,80,0)),142:((110,48,1,1,96,0),(98,48,2,1,96,0)),147:((115,48,1,1,96,0),(103,48,2,1,96,0)),152:((120,48,1,1,96,0),(108,48,2,1,96,0)),157:((141,24,1,1,96,0),(113,48,2,1,96,0)),162:((130,48,3,1,20,0),(118,48,2,1,96,0)),167:((135,48,3,1,20,0),(123,48,2,1,112,0)),172:((128,48,4,1,20,0),(128,48,2,1,112,0)),177:((121,72,4,1,14,0),(133,48,2,1,112,0)),182:((126,72,4,1,14,0),(126,48,3,1,96,0)),187:((131,72,4,1,12,0),(131,48,3,1,96,0)),192:((136,72,4,1,12,0),(136,48,3,1,96,0)),197:((141,72,4,1,12,0),(141,48,3,1,96,0)),202:((146,72,4,1,12,1),(146,48,3,1,112,0)),207:((151,72,4,1,12,1),(151,48,3,1,112,0)),212:((156,48,6,1,14,1),(156,48,3,1,112,0)),217:((161,48,6,1,14,1),(161,48,3,1,128,0)),222:((166,72,4,1,10,1),(166,48,3,1,128,0)),227:((219,48,6,0,10,0),(183,48,3,1,96,0)),232:((224,48,6,0,10,0),(188,48,3,1,112,0)),237:((221,48,6,0,12,0),(193,48,3,1,112,0)),242:((234,48,6,0,10,0),(186,72,6,1,10,1)),247:((239,48,6,0,10,0),(191,72,6,1,10,1)),252:((244,48,6,0,10,0),(196,72,6,1,10,1)),257:((249,48,6,0,10,0),(201,72,6,1,10,1))}
PLAN16={(b,s):v[i]for b,v in _PLAN.items()for i,s in enumerate('fg')}
BOOT16=(260,36,12,4)

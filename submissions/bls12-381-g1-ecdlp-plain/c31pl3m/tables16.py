"""Dense-band subset-sum tables (proof.md section 4, [10d]): ratio() = row entry cost with tables / without, from filtsim row
compositions (COMP in filt16.b85; band terms above the explicit primes split by filtsim's kappa(u, t)/t law; Poisson band counts, exact
greedy groups, densest band first; builds per row at C_PAIR (pairs) and C_RED (g >= 3), 4^g/2 patterns per multiset; nbytes(): one slot)."""
import math,json,lzma,base64,os,filtsim
L=69.5;US=(740.3/L,498.7/L);EDGES={2:31.5,3:21.,4:16.};U5,G1,C_PAIR,C_RED=7.8,33.,8.4,18.4;LN2=math.log(2)
def li(x):y=x*LN2;return 2**x/y*(1+1/y+2/y**2)
def ncols(x):return li(x)/2
def kap_t(t):return sum(filtsim.kappa(u,t)for u in US)/t
def comp():return json.loads(lzma.decompress(base64.b85decode(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'filt16.b85')).read().strip())))['COMP']
def bins(a,step=.25):
	b=[(k,min(k+1,a['lsmall']),x)for k,x in sorted((int(k),x)for k,x in a['uts'].items())];t0,t1,m=a['lsmall']/L,.8,400;tot=sum(kap_t(t0+(t1-t0)*(i+.5)/m)for i in range(m));x=a['lsmall']
	while x<.8*L-1e-9:y=min(x+step,.8*L);b.append((x,y,a['mid']*sum(kap_t((x+(y-x)*(i+.5)/40)/L)for i in range(40))*((y-x)/L/40)/((t1-t0)/m)/tot));x=y
	return b
def pmf(l):n=int(l+12*l**.5+40);p=[math.exp(-l)];[p.append(p[-1]*l/k)for k in range(1,n+1)];return p
def adds(lam,G):
	c={0:1.};s=.0
	for g in range(G,1,-1):
		nc={}
		for r,pr in c.items():
			for k,pk in enumerate(pmf(lam.get(g,.0))):
				w=pr*pk
				if w<1e-18:continue
				q=r+k;s+=w*(q//g);nc[q%g]=nc.get(q%g,.0)+w
		c=nc
	return s+sum(r*p for r,p in c.items())+lam.get(1,.0)
def nbytes(e=EDGES):return 160*sum(4**g/2*math.comb(int(ncols(x))+g-1,g)for g,x in e.items())
def tab_cost(e):return 8*math.comb(int(ncols(e[2]))+1,2)*C_PAIR+4*ncols(e[2])*C_RED+8*math.comb(int(ncols(e[3]))+1,2)*C_RED+sum(4**g/2*math.comb(int(ncols(x))+g-1,g)*C_RED for g,x in e.items()if g>=3)
def ratio(a,log2N,e=EDGES):
	lam={}
	for lo,hi,x in bins(a):
		cs=sorted({lo,hi,*[v for v in e.values()if lo<v<hi]})
		for p,q in zip(cs,cs[1:]):g=max([h for h,v in e.items()if(p+q)/2<=v],default=1);lam[g]=lam.get(g,0)+x*(q-p)/(hi-lo)
	units=sum(x for _,_,x in bins(a))+a['large'];old=U5*units+G1*a['gen'];return(U5*(adds(lam,max(e))+a['large'])+G1*a['gen']+tab_cost(e)/2**log2N)/old
def scale(rows):a=dict(lsmall=sum(r['lsmall']for r in rows)/len(rows),uts={},gen=0,mid=0,large=0);n=len(rows);[a['uts'].__setitem__(k,a['uts'].get(k,0)+x/n)for r in rows for k,x in r['uts'].items()];[a.__setitem__(k,a[k]+r[k]/n)for r in rows for k in('gen','mid','large')];return a

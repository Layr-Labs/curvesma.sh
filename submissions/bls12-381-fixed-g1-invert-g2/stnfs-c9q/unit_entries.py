""
import json,sys,os
_HERE=os.path.dirname(os.path.abspath(__file__))
_DIR=os.path.join(_HERE,"files")if os.path.exists(os.path.join(_HERE,"files","filtsim.py"))else _HERE
sys.path.insert(0,_DIR)
import filtsim
from filtsim import npt,U4
MID={};LAST=[]
MUF,MUG,EXCESS=744.0,501.0,1.095
CHARGED,CK=0.86,9.89
def cn(v):
	return next(w for w in(v*u for u in U4)if w.real>0 and w.imag>=0)
def cost(a):
	return 9.2*a if a<=4 else 40.0
def comb(i,p,m,pr):
	q=MID.get(p,{});c=cn(m)
	if c==1 and not q:return
	d=MID.setdefault(i,{})
	for k,x in([(1,pr[-2]-sum(q.values()))]if c!=1 else[])+list(q.items()):
		k=cn(c*k);d[k]=d.get(k,0.0)+x
def stats(rows,alive):
	h={};mu=mp=mc=0.0
	for i,r in enumerate(rows):
		if not alive[i]:continue
		for c,v in r.items():
			if c>=0:
				for a in(abs(v.real),abs(v.imag)):
					if a:h[a]=h.get(a,0)+1
		q=MID.get(i,{});u=r[-2]-sum(q.values());mu+=u;mp+=u;mc+=9.2*u
		for k,x in q.items():mp+=x*npt(k);mc+=x*sum(cost(abs(z))for z in(k.real,k.imag)if z)
	return dict(explicit=sum(h.values()),ones=h.get(1,0),xcost=sum(cost(a)*n for a,n in h.items()),mid=mp,mid_unit=mu,mid_cost=mc,hist={str(int(k)):h[k]for k in sorted(h)[:8]})
def main():
	req=json.load(sys.stdin)
	log2U=int(req.get("log2U",12));seed=int(req.get("seed",1));W=int(req.get("weight",400))
	assert 10<=log2U<=17
	MID.clear()
	L,ex=69.0,float(req.get("excess",EXCESS))
	cols,w,rows,w_mid=filtsim.build_system(L,(MUF/L,MUG/L),ex,seed,log2U)
	Usamp=sum(1 for(s_,t)in cols if t>=0.8)
	rows2,nc=filtsim.purge(rows,len(cols),max(2,int(filtsim.PEXC*Usamp)),False)
	n,ncl,avg,maxc,cw,st=filtsim.merge(rows2,len(cols),(W,),kmax=64,verbose=False,hook=sys.modules[__name__])[W]
	tot=st['explicit']+st['mid']
	unit=(st['ones']+st['mid_unit'])/tot
	k5=(st['xcost']+st['mid_cost'])/tot
	print(json.dumps(dict(log2U=log2U,seed=seed,U_over_N=round(Usamp/n,3),mean_row_weight=round(avg,2),max_coef=int(maxc),
	explicit_fraction=round(st['explicit']/tot,4),unit_of_explicit=round(st['ones']/st['explicit'],4),
	implicit_nonunit_of_implicit=round(1-st['mid_unit']/st['mid'],4),unit_fraction=round(unit,4),
	charged_unit_fraction=CHARGED,small_shares=[round(st['hist'].get(str(c),0)/tot,5)for c in(2,3,4)],k5=round(k5,3),charged_k5=CK,ok=unit>=CHARGED,explicit_hist=st['hist'])))
if __name__=="__main__":
	sys.modules.setdefault("unit_entries",sys.modules[__name__])
	main()

import json,sys,os,heapq,math
_HERE=os.path.dirname(os.path.abspath(__file__))
_DIR=os.path.join(_HERE,"files")if os.path.exists(os.path.join(_HERE,"files","filtsim.py"))else _HERE
sys.path.insert(0,_DIR)
import filtsim
from filtsim import rowweight
MIDNU={}
MUF,MUG,EXCESS=739.1,498.2,1.05
def coef_stats(rows,alive):
		ex=one=0;mid=mid_nu=0.0;hist={}
		for i,r in enumerate(rows):
				if not alive[i]:continue
				for c,v in r.items():
						if c>=0:
								ex+=1;a=abs(v);hist[a]=hist.get(a,0)+1
								if a==1:one+=1
				mid+=r.get(-1,0.0);mid_nu+=MIDNU.get(i,0.0)
		return dict(explicit=ex,ones=one,mid=mid,mid_nonunit=mid_nu,hist={str(k):hist[k]for k in sorted(hist)[:8]})
def merge_stats(rows,U,targets,kmax=32,cmax=255,verbose=True):
	rows=[dict(r)for r in rows];alive=[True]*len(rows);imul=[1]*len(rows);colrows={}
	for i,r in enumerate(rows):
		for c in r:
			if c>=0:colrows.setdefault(c,set()).add(i)
	nrows=len(rows);ncols=len(colrows);totw=sum(rowweight(r)for r in rows);heap=[(len(s),c)for c,s in colrows.items()];heapq.heapify(heap);results={};tlist=sorted(targets);ti=0;maxc=1;eliminated_cols=0
	while heap and ti<len(tlist):
		wgt,c=heapq.heappop(heap)
		if c not in colrows:continue
		s=colrows[c]
		if len(s)!=wgt:heapq.heappush(heap,(len(s),c));continue
		if wgt==0:del colrows[c];ncols-=1;continue
		if wgt>kmax:break
		if wgt==1:
			i=next(iter(s));alive[i]=False;nrows-=1;totw-=rowweight(rows[i])
			for cc in rows[i]:
				if cc<0:continue
				colrows[cc].discard(i)
				if not colrows[cc]:del colrows[cc];ncols-=1
				else:heapq.heappush(heap,(len(colrows[cc]),cc))
			rows[i]={};continue
		piv=[i for i in s if abs(rows[i][c])==1]
		if not piv:continue
		p=min(piv,key=lambda i:rowweight(rows[i]));pr=rows[p];pc=pr[c];ok=True
		for i in s:
			if i==p:continue
			r=rows[i];m=r[c]*pc
			if imul[i]+abs(m)*imul[p]>cmax:ok=False;break
			for cc,v in pr.items():
				if cc<0:continue
				if abs(r.get(cc,0)-m*v)>cmax:ok=False;break
			if not ok:break
		if not ok:continue
		for i in list(s):
			if i==p:continue
			r=rows[i];m=r[c]*pc;imul[i]+=abs(m)*imul[p]
			if imul[i]>maxc:maxc=imul[i]
			r[-1]=r.get(-1,0.0)+pr.get(-1,0.0);totw+=pr.get(-1,0.0);MIDNU[i]=MIDNU.get(i,0.0)+MIDNU.get(p,0.0)+((pr.get(-1,0.0)-MIDNU.get(p,0.0))if abs(m)>1 else 0.0)
			for cc,v in pr.items():
				if cc<0:continue
				nv=r.get(cc,0)-m*v
				if nv==0:
					if cc in r:
						del r[cc];colrows[cc].discard(i);totw-=1
						if not colrows[cc]:del colrows[cc];ncols-=1
						else:heapq.heappush(heap,(len(colrows[cc]),cc))
				else:
					if cc not in r:colrows[cc].add(i);totw+=1;heapq.heappush(heap,(len(colrows[cc]),cc))
					r[cc]=nv
					if abs(nv)>maxc:maxc=abs(nv)
		alive[p]=False;nrows-=1;totw-=rowweight(pr)
		for cc in pr:
			if cc==c or cc<0:continue
			colrows[cc].discard(p)
			if not colrows[cc]:del colrows[cc];ncols-=1
			else:heapq.heappush(heap,(len(colrows[cc]),cc))
		rows[p]={};del colrows[c];ncols-=1;eliminated_cols+=1;avg=totw/max(nrows,1)
		while ti<len(tlist)and avg>=tlist[ti]:results[tlist[ti]]=(nrows,ncols,avg,maxc,wgt,coef_stats(rows,alive));ti+=1
	if ti<len(tlist):
		avg=totw/max(nrows,1)
		for t in tlist[ti:]:results[t]=(nrows,ncols,avg,maxc,-1)
	return results
def main():
		req=json.load(sys.stdin);log2U=int(req.get("log2U",12));seed=int(req.get("seed",1));W=int(req.get("weight",300));assert 10<=log2U<=15;MIDNU.clear();L,uf,ug,ex=69.0,MUF/69,MUG/69,float(req.get("excess",EXCESS));cols,w,rows,w_mid=filtsim.build_system(L,(uf,ug),ex,seed,log2U);Usamp=sum(1 for(s_,t)in cols if t>=0.8);rows2,nc=filtsim.purge(rows,len(cols),max(2,int(filtsim.PEXC*Usamp)),False);res=merge_stats(rows2,len(cols),(W,),kmax=64,verbose=False);n,ncl,avg,maxc,cw,st=res[W];tot=st['explicit']+st['mid'];unit=(st['ones']+st['mid']-st['mid_nonunit'])/tot;charged=0.87;CK=12.47;sh=[st['hist'].get(str(c),0)/tot for c in(2,3,4)];k5=9.2*unit+sum(9.2*c*x for c,x in zip((2,3,4),sh))+40*(1-unit-sum(sh));print(json.dumps(dict(log2U=log2U,seed=seed,U_over_N=round(Usamp/n,3),mean_row_weight=round(avg,2),max_coef=maxc,
		           explicit_fraction=round(st['explicit']/tot,4),unit_of_explicit=round(st['ones']/st['explicit'],4),
		           implicit_nonunit_of_implicit=round(st['mid_nonunit']/st['mid'],4),unit_fraction=round(unit,4),
		           charged_unit_fraction=charged,small_shares=[round(x,5)for x in sh],k5=round(k5,3),charged_k5=CK,ok=unit>=charged and k5<=CK,explicit_hist=st['hist'])))
if __name__=="__main__":
		main()

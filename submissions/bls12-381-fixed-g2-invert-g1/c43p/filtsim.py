"";import math,random,bisect,heapq;LN2=math.log(2);from decimal import Decimal,getcontext;Vi=40;Vs=600;Vf=.001
def Vh():
	getcontext().prec=150;zero,one=Decimal(0),Decimal(1);coefs=[None,[one]+[zero]*(Vs-1)]
	for k in range(2,Vi+1):
		prev=coefs[k-1];c=[zero]*Vs
		for m in range(Vs-1):c[m+1]=(prev[m]+m*c[m])/(k*(m+1))
		c[0]=prev[0]-sum(c[1:]);coefs.append(c)
	out=[None]
	for c in coefs[1:]:
		big=max(abs(x)for x in c)
		while len(c)>1 and abs(c[-1])<big*Decimal('1e-25'):c=c[:-1]
		out.append([float(x)for x in c])
	return out
Vy=Vh()
def rho_series(u):
	if u<=1.:return 1.
	k=math.ceil(u)
	if k>Vi or k<1:return.0
	xi=k-u;s=.0
	for ci in reversed(Vy[k]):s=s*xi+ci
	return s if s>.0 else.0
Vq=int((Vi-.001)/Vf);Vw=[1. if i*Vf<=1. else rho_series(i*Vf)for i in range(Vq+2)]
def rho(u):
	if u<=1.:return 1.
	x=u/Vf
	if x>=Vq:return.0
	i=int(x);f=x-i;return Vw[i]*(1.-f)+Vw[i+1]*f
def primes_upto(n):
	s=bytearray([1])*(n+1);s[0]=s[1]=0
	for i in range(2,int(n**.5)+1):
		if s[i]:s[i*i::i]=bytearray(len(s[i*i::i]))
	return[i for i in range(n+1)if s[i]]
U4=1,1j,-1,-1j
def n1(v):return abs(v.real)+abs(v.imag)
def mx(v):return max(abs(v.real),abs(v.imag))
def npt(v):return(v.real!=0)+(v.imag!=0)
def kappa(u,t):return rho(u-t)/rho(u)
def lsm(L,Rp,c4,t=.8):return min(max(8.,math.log2(Rp/25.))+c4,t*L)
def build_system(L,us,excess,seed,log2U,t_dense=.8,nt=4000,lsmall=None,c4=True):
	rng=random.Random(seed);VT=2.**L/(L*LN2);VM=excess*VT;Up=2**log2U;Rp=int(round(excess*Up))
	if lsmall is None:lsmall=lsm(L,Rp,c4,t_dense)
	ts=[t_dense+(1-t_dense)*(i+.5)/nt for i in range(nt)];VR=[2.**(L*(t-1))/t for t in ts];cum=[];tot=.0
	for d in VR:tot+=d;cum.append(tot)
	cum=[c/tot for c in cum];cols=[];w=[]
	for s_,u in enumerate(us):
		for _ in range(Up//2):i=bisect.bisect_left(cum,rng.random());t=ts[min(i,nt-1)]+(rng.random()-.5)*(1-t_dense)/nt;lam=VM*2*kappa(u,t)/2.**(L*t);cols.append((s_,t));w.append(min(1.,lam/Rp))
	VF=[]
	for s_,u in enumerate(us):
		for q in primes_upto(int(2**lsmall))[::1+c4]:t=math.log(q)/(L*LN2);pr=min(1.,(2+2*c4)*kappa(u,t)/q);cols.append((s_,t));w.append(pr);VF.append((len(cols)-1,pr,q))
	VI=lsmall/L;w_mid=.0;m=2000
	for s_,u in enumerate(us):
		for i in range(m):t=VI+(t_dense-VI)*(i+.5)/m;w_mid+=kappa(u,t)*((t_dense-VI)/m)/t
	U=len(cols);VG=Up//2*len(us);VL=[i for i in range(VG)if w[i]>=.2];VC=[i for i in range(VG)if w[i]<.2];VK=[];VD=.0
	for i in VC:VD+=w[i];VK.append(VD)
	sg=(lambda:U4[int(4*rng.random())])if c4 else lambda:1 if rng.random()<.5 else-1;rows=[]
	for _ in range(Rp):
		row={}
		for i in VL:
			if rng.random()<w[i]:row[i]=sg()
		k=Vt(rng,VD)
		for _ in range(k):
			x=rng.random()*VD;j=VC[bisect.bisect_left(VK,x)];row[j]=row.get(j,0)+sg()
			if row[j]==0:del row[j]
		for(j,pr,q)in VF:
			if rng.random()<pr:
				v=1
				while rng.random()<1./q:v+=1
				row[j]=v*sg()
		row[-1]=row[-2]=w_mid;rows.append(row)
	return cols,w,rows,w_mid
def Vt(rng,lam):
	if lam<30:
		L=math.exp(-lam);k=0;p=1.
		while True:
			p*=rng.random()
			if p<=L:return k
			k+=1
	return max(0,int(round(rng.gauss(lam,math.sqrt(lam)))))
def Vj(rows,alive,Vm,Ve):
	Vd=True
	while Vd:
		Vd=False
		for i,r in enumerate(rows):
			if not alive[i]:continue
			for c in r:
				if c>=0 and Vm[c]==1:
					alive[i]=False;Ve-=1;Vd=True
					for cc in r:
						if cc>=0:Vm[cc]-=1
					break
	return Ve
def Vc(r):return sum(npt(v)for c,v in r.items()if c>=0)+r.get(-1,.0)
def purge(rows,U,Vg,verbose=True):
	Vm=[0]*U
	for r in rows:
		for c in r:
			if c>=0:Vm[c]+=1
	alive=[True]*len(rows);Ve=Vj(rows,alive,Vm,len(rows));ncols=sum(1 for x in Vm if x>0)
	while Ve-ncols>Vg:
		parent=list(range(len(rows)))
		def VJ(x):
			while parent[x]!=x:parent[x]=parent[parent[x]];x=parent[x]
			return x
		Vv={}
		for i,r in enumerate(rows):
			if not alive[i]:continue
			for c in r:
				if c>=0 and Vm[c]==2:Vv.setdefault(c,[]).append(i)
		for c,rs in Vv.items():
			if len(rs)==2:
				a,b=VJ(rs[0]),VJ(rs[1])
				if a!=b:parent[a]=b
		comp={}
		for i in range(len(rows)):
			if alive[i]:comp.setdefault(VJ(i),[]).append(i)
		VH=sorted(comp.values(),key=lambda rs:-sum(Vc(rows[i])for i in rs));VP=Ve-ncols-Vg;nb=max(1,min(VP//2+1,max(1,len(VH)//100)));Vz=0
		for rs in VH:
			if Vz>=nb:break
			for i in rs:
				if alive[i]:
					alive[i]=False;Ve-=1
					for c in rows[i]:
						if c>=0:Vm[c]-=1
			Vz+=1
		Ve=Vj(rows,alive,Vm,Ve);ncols=sum(1 for x in Vm if x>0)
	return[rows[i]for i in range(len(rows))if alive[i]],ncols
def merge(rows,U,targets,kmax=32,cmax=255,verbose=True,hook=None):
	rows=[dict(r)for r in rows];alive=[True]*len(rows);VB=[1]*len(rows);Va={}
	for i,r in enumerate(rows):
		for c in r:
			if c>=0:Va.setdefault(c,set()).add(i)
	Vo=len(rows);ncols=len(Va);VA=sum(Vc(r)for r in rows);rw=[Vc(r)for r in rows];K=lambda c:(lambda s:(len(s)-2)*min([rw[i]for i in s if rows[i][c]in U4]or[1e9])if len(s)>2 else len(s)-3)(Va[c]);P=lambda c:len(Va[c])>kmax or heapq.heappush(VE,(K(c),c));VE=[(K(c),c)for c in Va if len(Va[c])<=kmax];heapq.heapify(VE);Vr={};Vn=sorted(targets);ti=0;maxc=1
	while VE and ti<len(Vn):
		key,c=heapq.heappop(VE)
		if c not in Va:continue
		s=Va[c];k=K(c);VN=len(s)
		if k>key+1e-9:P(c);continue
		if VN==0:del Va[c];ncols-=1;continue
		if VN>kmax:continue
		if VN==1:
			i=next(iter(s));alive[i]=False;Vo-=1;VA-=Vc(rows[i])
			for cc in rows[i]:
				if cc<0:continue
				Va[cc].discard(i)
				if not Va[cc]:del Va[cc];ncols-=1
				else:P(cc)
			rows[i]={};continue
		piv=[i for i in s if abs(rows[i][c])==1]
		if not piv:continue
		p=min(piv,key=lambda i:Vc(rows[i]));pr=rows[p];pc=pr[c];ok=True
		for i in s:
			if i==p:continue
			r=rows[i];m=r[c]*pc.conjugate()
			if VB[i]+n1(m)*VB[p]>cmax:ok=False;break
			for cc,v in pr.items():
				if cc<0:continue
				if mx(r.get(cc,0)-m*v)>cmax:ok=False;break
			if not ok:break
		if not ok:continue
		for i in list(s):
			if i==p:continue
			r=rows[i];m=r[c]*pc.conjugate();VB[i]+=n1(m)*VB[p]
			if hook:hook.comb(i,p,m,pr)
			if VB[i]>maxc:maxc=VB[i]
			dm=pr[-1]if npt(m)==1 else 2*pr[-2];r[-1]+=dm;r[-2]+=pr[-2];VA+=dm
			for cc,v in pr.items():
				if cc<0:continue
				ov=r.get(cc,0);nv=ov-m*v;VA+=npt(nv)-npt(ov)
				if nv==0:
					if cc in r:
						del r[cc];Va[cc].discard(i)
						if not Va[cc]:del Va[cc];ncols-=1
						else:P(cc)
				else:
					if cc not in r:Va[cc].add(i)
					r[cc]=nv
					if mx(nv)>maxc:maxc=int(mx(nv))
			rw[i]=Vc(r)
			for cc in r:
				if cc>=0:P(cc)
		alive[p]=False;Vo-=1;VA-=Vc(pr)
		for cc in pr:
			if cc==c or cc<0:continue
			Va[cc].discard(p)
			if not Va[cc]:del Va[cc];ncols-=1
			else:P(cc)
		rows[p]={};del Va[c];ncols-=1;avg=VA/max(Vo,1)
		while ti<len(Vn)and avg>=Vn[ti]:Vr[Vn[ti]]=(Vo,ncols,avg,maxc,VN)+((hook.stats(rows,alive),)if hook else());ti+=1
	if ti<len(Vn):
		avg=VA/max(Vo,1)
		for t in Vn[ti:]:Vr[t]=Vo,ncols,avg,maxc,-1
	if hook:hook.LAST[:]=[r for i,r in enumerate(rows)if alive[i]]
	return Vr
PEXC=.0001
def run(L,uf,ug,excess,seed,log2U,pexc=PEXC,kmax=32,targets=(100,150,200,300,400,600),verbose=True,c4=True):
	cols,w,rows,w_mid=build_system(L,(uf,ug),excess,seed,log2U,c4=c4);Usamp=sum(1 for(s_,t)in cols if t>=.8);U=len(cols);rows2,ncols=purge(rows,U,max(2,int(pexc*Usamp)),verbose);res=merge(rows2,U,targets,kmax=kmax,verbose=verbose);out={}
	for W,(n,nc,avg,maxc,cw)in sorted(res.items()):out[W]=Usamp/n,n,maxc,avg
	return Usamp,out
def missing_recon(seed,log2U=12,L=69.,uf=750./69,ug=498.5/69,excess=1.3,pexc=PEXC,c4=True):
	cols,w,rows,_=build_system(L,(uf,ug),excess,seed,log2U,c4=c4);U=len(cols);Rp=len(rows);VO=[0]*U
	for r in rows:
		for c in r:
			if c>=0:VO[c]+=1
	kept=purge([dict(r)for r in rows],U,max(2,int(pexc*U)),verbose=False)[0];Vu=[False]*U
	for r in kept:
		for c in r:
			if c>=0:Vu[c]=True
	Vd=True
	while Vd:
		Vd=False
		for r in rows:
			VS=[c for c in r if c>=0 and not Vu[c]]
			if len(VS)==1:Vu[VS[0]]=True;Vd=True
	out=[]
	for lo,hi in(1,2),(2,3),(3,5),(5,8),(8,1e9):
		idx=[i for i in range(2*(2**log2U//2))if lo<=w[i]*Rp<hi and VO[i]>0]
		if not idx:continue
		VQ=sum(1 for i in idx if not Vu[i])/len(idx);model=sum(min(1.,math.exp(-w[i]*Rp)*(1+w[i]*Rp+(w[i]*Rp)**2/2))for i in idx)/len(idx);out.append((lo,hi,len(idx),VQ,model,sum(1 for i in idx if not Vu[i]and VO[i]>=3)))
	return out

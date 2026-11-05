"""H3x toy (proof.md 5): left Krylov fill of unit-vector X over F_r. Stdin {"dim","m","trials","seed","w"}. A: production-
shaped (Zipf columns, 87% +-1, column weight >= 2, 21 dense F_r columns, rank dim-1). Fill: rank [X..(A^T)^(k-1) X] = dim."""
import json,random,sys
R=52435875175126190479447740508185965837690552500527637822603658699938581184513
def mat(D,w,g):
	n=D-21;cs=list(range(n));cw=[1/(i+8)for i in cs];rows=[]
	for _ in range(D):
		r={}
		for c in g.choices(cs,cw,k=w):r[c]=g.choice((1,-1))*(1 if g.random()<0.87 else g.randint(2,8))
		rows.append(r)
	for c in cs:
		while sum(c in r for r in rows)<2:rows[g.randrange(D)][c]=g.choice((1,-1))
	for r in rows:
		for c in range(n,D):r[c]=g.randrange(R)
	v=[g.randrange(1,R)for _ in range(D)]
	for r in rows:r[n]=-sum(a*v[c]for c,a in r.items()if c!=n)*pow(v[n],-1,R)%R
	return rows
def add(B,v,D):
	p=next((i for i,a in enumerate(v)if a),None)
	while p is not None and p in B:
		f=v[p];v=[(a-f*c)%R for a,c in zip(v,B[p])];p=next((i for i in range(p,D)if v[i]),None)
	if p is not None:B[p]=[a*pow(v[p],-1,R)%R for a in v]
def fill(rows,D,J,k):
	B={};X=[[int(i==j)for i in range(D)]for j in J]
	for s in range(k):
		for x in X:add(B,x,D)
		if len(B)==D:return s+1
		Y=[[0]*D for _ in X]
		for x,y in zip(X,Y):
			for i,xi in enumerate(x):
				if xi:
					for c,a in rows[i].items():y[c]+=xi*a
		X=[[t%R for t in y]for y in Y]
def run(D,m,T,s,w):
	g=random.Random(s);k=-(-D//m);res=[];bad=0
	while len(res)<T:
		rows=mat(D,w,g);B={}
		for r in rows:add(B,[r.get(c,0)%R for c in range(D)],D)
		if len(B)!=D-1:bad+=1;continue
		res.append(fill(rows,D,g.sample(range(D),m),k+1))
	ok=res.count(k)
	return dict(dim=D,m=m,w=w,k=k,steps={str(x):res.count(x)for x in set(res)},filled_at_k=ok,trials=T,h3_redrawn=bad,ok=ok==T)
if __name__=="__main__":
	q=json.load(sys.stdin);D=int(q.get("dim",96));m=int(q.get("m",8));T=int(q.get("trials",10));w=int(q.get("w",24))
	assert 48<=D<=160 and 2<=m<=32 and 1<=T<=40 and 4<=w<=40
	print(json.dumps(run(D,m,T,int(q.get("seed",1)),w)))

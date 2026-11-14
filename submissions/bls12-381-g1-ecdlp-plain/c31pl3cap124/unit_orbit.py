"";import math,itertools,cmath;from fractions import Fraction as Fr;import common12 as C;H=C.H
def rows_of(v):return[C.mul([0]*j+[1],v[:6])+C.mul([0]*j+[1],v[6:])for j in range(6)]
def lattice_adj(B):
	n=len(B);a=[[Fr(x)for x in row]+[Fr(int(i==j))for j in range(n)]for i,row in enumerate(B)]
	for k in range(n):
		piv=next(i for i in range(k,n)if a[i][k]);a[k],a[piv]=a[piv],a[k];p=a[k][k];a[k]=[x/p for x in a[k]]
		for i in range(n):
			if i!=k and a[i][k]:f=a[i][k];a[i]=[x-f*y for x,y in zip(a[i],a[k])]
	den=1
	for row in a:
		for x in row[n:]:den=den*x.denominator//math.gcd(den,x.denominator)
	return[[int(x*den)for x in row[n:]]for row in a],den
def ball_coords(w,adj,den):cm=[sum(w[i]*adj[i][j]for i in range(len(adj)))for j in range(len(adj[0]))];return None if any(x%den for x in cm)else[x//den for x in cm]
def enum_form(G,bound):
	n=len(G);Gf=[[Fr(x)for x in row]for row in G];d=[Fr(0)]*n;L=[[Fr(0)]*n for _ in range(n)]
	for i in range(n):
		s=Gf[i][i]-sum(L[i][k]**2*d[k]for k in range(i))
		if s<=0:raise ValueError('form not positive definite')
		d[i]=s
		for j in range(i+1,n):L[j][i]=(Gf[j][i]-sum(L[j][k]*L[i][k]*d[k]for k in range(i)))/d[i]
	out=[];e=[0]*n
	def rec(i,rem):
		c=sum(L[j][i]*e[j]for j in range(i+1,n));r2=rem/d[i];s=math.isqrt(r2.numerator//r2.denominator)if r2>0 else 0
		for x in range(math.floor(-c)-s-1,math.ceil(-c)+s+2):
			t=d[i]*(x+c)**2
			if t<=rem:
				e[i]=x
				if i==0:out.append(e[:])
				else:rec(i-1,rem-t)
		e[i]=0
	rec(n-1,Fr(bound));return out
def cell_form(v,adj,den):rows=rows_of(v);M=[[sum(r[k]*adj[k][l]for k in range(12))for l in range(12)]for r in rows];return[[sum(M[i][l]*M[j][l]for l in range(12))for j in range(6)]for i in range(6)]
def orbit(v,adj,den,R):
	pts=enum_form(cell_form(v,adj,den),R*R*den*den);units=[e for e in pts if abs(C.norm(e))==1]
	for e in units:c=ball_coords(C.mul(e,v[:6])+C.mul(e,v[6:]),adj,den);assert c is not None and sum(x*x for x in c)<=R*R
	assert[1,0,0,0,0,0]in units and[-1,0,0,0,0,0]in units;return units,len(pts)
def height(e):return max(abs(x)for x in e)
def m_true(units):return len(units)//2

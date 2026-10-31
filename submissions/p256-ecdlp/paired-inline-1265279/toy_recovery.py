"""Actual small a=-3 curves, only evidence for the recovery procedure."""
import random
def prime(n):return n>=2 and all(n%d for d in range(2,int(n**.5)+1))
def check():
    p=1009;a=p-3
    for b in range(1,p):
        if (4*a**3+27*b*b)%p==0:continue
        points=[]
        for x in range(p):
            rhs=(x*x*x+a*x+b)%p
            if rhs==0:points.append((x,0))
            elif pow(rhs,(p-1)//2,p)==1:
                y=next(y for y in range(1,p) if y*y%p==rhs)
                points.extend([(x,y),(x,-y%p)])
        r=len(points)+1
        if prime(r):break
    else:raise AssertionError('no test curve')
    G=points[0]
    def add(A,B):
        if A is None:return B
        if B is None:return A
        x,y=A;u,v=B
        if x==u:
            if (y+v)%p==0:return None
            t=(3*x*x+a)*pow(2*y,-1,p)%p
        else:t=(v-y)*pow(u-x,-1,p)%p
        xx=(t*t-x-u)%p
        return xx,(t*(x-xx)-y)%p
    def mul(k,A):
        out=None
        while k:
            if k&1:out=add(out,A)
            A=add(A,A);k>>=1
        return out
    assert mul(r,G) is None
    rng=random.Random(271828);wins=aborts=miss=0;adjacency_checks=0
    for _ in range(120):
        secret=rng.randrange(r);Q=mul(secret,G)
        if Q is None:wins+=1;continue
        groups={};bad=False
        for node in range(4):
            bg=rng.randrange(r);cs=[rng.randrange(1,r) for _ in range(2)]
            for k in range(2):
                ag=rng.randrange(r);X=add(mul(ag,G),mul(bg,Q))
                for c in cs:
                    Y=mul(c,G)
                    if X is None or X[0]==Y[0]:bad=True;continue
                    for sign in [1,-1]:
                        L=add(X,(Y[0],sign*Y[1]%p))
                        assert L is not None
                        groups.setdefault(L[0],[]).append((L,(ag+sign*c)%r,bg))
        if bad:aborts+=1;continue
        answer=None
        for group in groups.values():
            useful_any=False;useful_adj=False
            for i in range(len(group)):
                for j in range(i+1,len(group)):
                    A,m,n=group[i];B,mm,nn=group[j]
                    u=1 if A==B else -1
                    assert B==(A[0],u*A[1]%p)
                    if (nn-u*n)%r:
                        useful_any=True
                        if j==i+1:useful_adj=True
            assert not useful_any or useful_adj
            adjacency_checks+=1
            for left,right in zip(group,group[1:]):
                A,m,n=left;B,mm,nn=right;u=1 if A==B else -1
                den=(nn-u*n)%r
                if den:
                    x=(u*m-mm)*pow(den,-1,r)%r
                    assert mul(x,G)==Q and x==secret
                    answer=x
        if answer is None:miss+=1
        else:wins+=1
    return {'curve_p':p,'a':a,'b':b,'prime_order':r,'trials':120,
            'recoveries':wins,'denominator_or_identity_aborts':aborts,
            'no_useful_candidate':miss,'adjacency_checks':adjacency_checks,
            'all_recovered_scalars_correct':True,
            'scope':'Tiny actual curve; tests recovery and adjacency, not production success.'}

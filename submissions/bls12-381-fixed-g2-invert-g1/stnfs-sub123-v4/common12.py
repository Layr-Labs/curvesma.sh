from fractions import Fraction as F
from math import isqrt, ceil, lcm
from decimal import Decimal, localcontext
H = [-1,-1,1,1,0,0,1]
U = 15132376222941642752
TERMS = [(12,[1]),(10,[0,-2]),(6,[0,0,0,2]),(2,[0,0,0,0,0,1]),(0,[1,1,-1,-1])]
def mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            c[i+j]+=x*y
    for i in range(len(c)-1,5,-1):
        for j in range(6):
            c[i-6+j]-=c[i]*H[j]
    return (c+[0]*6)[:6]
def norm(a):
    cols=[mul(a,[0]*i+[1]) for i in range(6)]
    m=[[cols[j][i] for j in range(6)] for i in range(6)]
    prev, sign=1,1
    for k in range(5):
        if m[k][k]==0:
            row=next((j for j in range(k+1,6) if m[j][k]),None)
            if row is None:
                return 0
            m[k],m[row]=m[row],m[k]; sign=-sign
        for i in range(k+1,6):
            for j in range(k+1,6):
                value=m[i][j]*m[k][k]-m[i][k]*m[k][j]
                assert value%prev==0
                m[i][j]=value//prev
        prev=m[k][k]
    return abs(sign*m[-1][-1])
def norms(a,b):
    ap,bp=[[1]],[[1]]
    for _ in range(12):
        ap.append(mul(ap[-1],a)); bp.append(mul(bp[-1],b))
    f=[0]*6
    for d,c in TERMS:
        term=mul(c,mul(ap[d],bp[12-d]))
        f=[x+y for x,y in zip(f,term)]
    g=[x+y for x,y in zip(ap[2],mul([0,U],bp[2]))]
    return norm(f),norm(g)
def integer_interval(a, t, c):
    a, t, c = F(a), F(t), F(c)
    if a <= 0:
        raise ValueError('non-positive pivot')
    centre = -t/a
    radius2 = (t*t-a*c)/(a*a)
    if radius2 < 0:
        return range(0)
    n, d = centre.numerator, centre.denominator
    z = isqrt((radius2.numerator*d*d)//radius2.denominator)
    lo = -((-(n-z))//d)
    hi = (n+z)//d
    return range(lo, hi+1)
def scaled_log(n, scale=16, upper=False):
    if n < 1 or scale < 1:
        raise ValueError('positive norm and scale required')
    x = n**scale
    floor = x.bit_length()-1
    return floor + int(upper and x != 1 << floor)
def inverse(matrix):
    n = len(matrix)
    a = [[F(x) for x in row]+[F(i == j) for j in range(n)]
         for i, row in enumerate(matrix)]
    for k in range(n):
        j = next(j for j in range(k, n) if a[j][k])
        a[k], a[j] = a[j], a[k]
        p = a[k][k]
        a[k] = [x/p for x in a[k]]
        for j in range(n):
            if j != k:
                p = a[j][k]
                a[j] = [x-p*y for x, y in zip(a[j], a[k])]
    return [row[n:] for row in a]
def integer_metric(matrix):
    den = lcm(*(F(x).denominator for row in matrix for x in row))
    q = [[int(F(x)*den) for x in row] for row in matrix]
    return q, den, inverse(q)
def quadratic(q, a):
    return sum(a[i]*q[i][j]*a[j] for i in range(6) for j in range(6))
def rho_coefficients(kmax=40, nt=600):
    with localcontext() as ctx:
        ctx.prec = 150
        zero, one = Decimal(0), Decimal(1)
        coefs = [None, [one]+[zero]*(nt-1)]
        for k in range(2, kmax+1):
            previous, c = coefs[k-1], [zero]*nt
            for m in range(nt-1):
                c[m+1] = (previous[m]+m*c[m])/(k*(m+1))
            c[0] = previous[0]-sum(c[1:])
            coefs.append(c)
        out = [None]
        for c in coefs[1:]:
            big = max(abs(x) for x in c)
            while len(c) > 1 and abs(c[-1]) < big*Decimal('1e-25'):
                c = c[:-1]
            out.append([float(x) for x in c])
        return out
COEFS = rho_coefficients()
def rho(u):
    if u <= 1:
        return 1.0
    k = ceil(u)
    if k >= len(COEFS):
        return 0.0
    x, value = k-u, 0.0
    for c in reversed(COEFS[k]):
        value = value*x+c
    return max(value, 0.0)

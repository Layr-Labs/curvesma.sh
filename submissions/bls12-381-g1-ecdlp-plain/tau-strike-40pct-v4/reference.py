"""Same orbit6-v3 map; Jacobian scalar action with one affine normalization.

Python big integers are reference arithmetic, not unit-cost native operations.
The independent affine routines below are also the toy correctness oracle.
"""

def add(P, Q, p):
    if P is None:
        return Q
    if Q is None:
        return P
    x, y = P
    X, Y = Q
    if x == X:
        if (y+Y) % p == 0:
            return None
        num, den = 3*x*x % p, 2*y % p
    else:
        num, den = (Y-y) % p, (X-x) % p
    slope = num*pow(den,p-2,p) % p
    xx = (slope*slope-x-X) % p
    return xx, (slope*(x-xx)-y) % p


def mul(k, P, p):
    Q = None
    for i in range(k.bit_length()-1,-1,-1):
        Q = add(Q,Q,p)
        if (k >> i) & 1:
            Q = add(Q,P,p)
    return Q


def projective_action(Q, k, G, p, census=False):
    """Return Q+[k]G. All mixed-add exceptions precede any division.

    Doubling: 7 products. Mixed addition: 11, including fallback 4+7.
    Normalization: <=2*bit_length(p)+4 products, including Fermat inverse.
    """
    count = {'doublings':0,'mixed_calls':0,'field_products':0,
             'normalizations':0,'mixed_double':0,'mixed_opposite':0}
    def f(a,b):
        count['field_products'] += 1
        return a*b % p
    def double(J):
        count['doublings'] += 1
        if J is None or J[1] == 0:
            return None
        X,Y,Z = J
        A = f(X,X)
        B = f(Y,Y)
        C = f(B,B)
        S = 4*f(X,B) % p
        M = 3*A % p
        XX = (f(M,M)-2*S) % p
        YY = (f(M,(S-XX)%p)-8*C) % p
        ZZ = 2*f(Y,Z) % p
        return XX,YY,ZZ
    def mixed(J,P):
        count['mixed_calls'] += 1
        if P is None:
            return J
        if J is None:
            return P[0],P[1],1
        X,Y,Z = J
        x,y = P
        Z2 = f(Z,Z)
        U = f(x,Z2)
        V = f(f(y,Z),Z2)
        H = (U-X) % p
        D = (V-Y) % p
        if H == 0:
            if D == 0:
                count['mixed_double'] += 1
                return double(J)
            count['mixed_opposite'] += 1
            return None
        HH = f(H,H)
        HHH = f(H,HH)
        XHH = f(X,HH)
        XX = (f(D,D)-HHH-2*XHH) % p
        YY = (f(D,(XHH-XX)%p)-f(Y,HHH)) % p
        ZZ = f(Z,H)
        return XX,YY,ZZ
    J = None
    for i in range(k.bit_length()-1,-1,-1):
        J = double(J)
        if (k >> i) & 1:
            J = mixed(J,G)
    J = mixed(J,Q)
    if J is None:
        out = None
    else:
        X,Y,Z = J
        count['normalizations'] += 1
        # Literal binary Fermat path, not an unpriced arbitrary inverse.
        inv = 1
        for i in range((p-2).bit_length()-1,-1,-1):
            inv = f(inv,inv)
            if ((p-2) >> i) & 1:
                inv = f(inv,Z)
        inv2 = f(inv,inv)
        out = f(X,inv2), f(Y,f(inv2,inv))
    return (out,count) if census else out


def canonical(P,c):
    if P is None:
        return bytes(49),None,1
    p,r = c['p'],c['r']
    x,y = P
    bx = c['beta']*x % p
    xs = [x,bx,(-x-bx)%p]
    j = xs.index(min(xs))
    sign = 1 if 2*y < p else -1
    yy = sign*y % p
    u = sign*c['units'][j] % r
    return b'\x01'+yy.to_bytes(48,'big'),(xs[j],yy),u


def toy_context(p,r,cofactor):
    beta = next(t for g in range(2,p) if (t:=pow(g,(p-1)//3,p)) != 1)
    # Only toy setup: exact residue-to-root enumeration, never a production table.
    roots = {}
    for y in range(1,p):
        roots.setdefault(y*y%p,y)
    G = None
    for x in range(p):
        y = roots.get((x*x*x+4)%p)
        if y is not None:
            G = mul(cofactor,(x,y),p)
            if G is not None:
                break
    assert G is not None and mul(r,G,p) is None
    lam = next(k for k in range(2,r) if (k*k+k+1)%r == 0
               and mul(k,G,p) == (beta*G[0]%p,G[1]))
    return {'p':p,'r':r,'beta':beta,'lam':lam,'G':G,
            'units':(1,lam,lam*lam%r)}


def budget(r):
    from math import isqrt
    N = (r-1)//6
    B = isqrt(169*N//100)
    return B+(100*B*B < 169*N)


def solve(T,c,seed,check=False,affine_steps=False,observe=None):
    """One attempt, same domains and state map as accepted orbit6 v3.

    affine_steps is an independent experimental comparison, not the scored path.
    observe only records toy histograms; neither diagnostic alters decisions.
    """
    from hashlib import shake_256
    p,r,G = c['p'],c['r'],c['G']
    B = budget(r)
    count = 0
    def done(x,reason,mu=None,period=None):
        if x is not None:
            assert 0 <= x < r and mul(x,G,p) == T
        return {'x':x,'reason':reason,'evaluations':count,
                'mu':mu,'period':period,'budget':B}
    if T is None:
        return done(0,'identity target')
    assert len(seed) == 32
    raw = shake_256(b'initial-v3'+seed).digest(128)
    a,b = int.from_bytes(raw[:64],'big')%r,int.from_bytes(raw[64:],'big')%r
    if b == 0:
        return done(None,'zero initial target coefficient')
    P = add(mul(a,G,p),mul(b,T,p),p)
    _,P,u = canonical(P,c)
    initial = (P,u*a%r,u*b%r)
    if P is None:
        return done(-a*pow(b,r-2,r)%r,'initial identity relation')
    def key(state):
        Q = state[0]
        return bytes(49) if Q is None else b'\x01'+Q[1].to_bytes(48,'big')
    def step(state):
        nonlocal count
        assert count < 6*B
        count += 1
        Q,aa,bb = state
        h = int.from_bytes(shake_256(b'orbit6-v3'+seed+key(state)).digest(64),'big')%r
        Q = add(Q,mul(h,G,p),p) if affine_steps else projective_action(Q,h,G,p)
        _,Q,u = canonical(Q,c)
        if observe is not None:
            observe(h,u)
        new = (Q,u*(aa+h)%r,u*bb%r)
        if check:
            assert add(mul(new[1],G,p),mul(new[2],T,p),p) == Q
        return new
    slow = fast = initial
    met = False
    for _ in range(B):
        slow = step(slow)
        fast = step(step(fast))
        if key(slow) == key(fast):
            met = True
            break
    if not met:
        return done(None,'Floyd iteration cap')
    left,right,mu = initial,slow,0
    while key(left) != key(right):
        if mu >= B:
            return done(None,'entry cap')
        left,right = step(left),step(right)
        mu += 1
    entry = left
    lap,period = step(entry),1
    while key(lap) != key(entry):
        if period >= B:
            return done(None,'period cap')
        lap = step(lap)
        period += 1
    den = (lap[2]-entry[2])%r
    if den == 0:
        return done(None,'one-lap orbit product is identity',mu,period)
    return done((entry[1]-lap[1])*pow(den,r-2,r)%r,
                'verified one-lap relation',mu,period)


def wire_output(x,r):
    if not isinstance(x,int) or not 0 <= x < r:
        raise ValueError('noncanonical scalar')
    return {'x':format(x,'064x')}

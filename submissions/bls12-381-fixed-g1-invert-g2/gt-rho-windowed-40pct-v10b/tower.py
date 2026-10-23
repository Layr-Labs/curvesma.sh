#!/usr/bin/env python3
"""Tower arithmetic Fp2 -> Fp6 -> Fp12 with BLS12-style xi=1+u; used for Part A and the toy."""
ZERO2=(0,0); ONE2=(1,0)

def is_prime(n):
    if n<2: return False
    for sp in (2,3,5,7,11,13,17,19,23,29,31,37):
        if n%sp==0: return n==sp
    d=n-1; s=0
    while d%2==0: d//=2; s+=1
    for a in (2,3,5,7,11,13,17,19,23,29,31,37):
        x=pow(a,d,n)
        if x in (1,n-1): continue
        for _ in range(s-1):
            x=x*x%n
            if x==n-1: break
        else: return False
    return True

class Tower:
    def __init__(self,p):
        assert p%12==7, p%12
        self.p=p
        self.xi=(1,1)
        self.xi3=self.f2pow(self.xi,(p-1)//3)          # v^p = xi3*v
        self.gamma=self.f2pow(self.xi,(p-1)//6)        # w^p = gamma*w
        self.zeta=(self.gamma[0]*self.gamma[0]+self.gamma[1]*self.gamma[1])%p  # |gamma|^2, 6th root of 1
        self.zero6=(ZERO2,ZERO2,ZERO2)
        self.one6=(ONE2,ZERO2,ZERO2)
    # ---- Fp2
    def f2mul(self,a,b):
        p=self.p; return ((a[0]*b[0]-a[1]*b[1])%p,(a[0]*b[1]+a[1]*b[0])%p)
    def f2sqr(self,a):
        p=self.p; return (((a[0]+a[1])*(a[0]-a[1]))%p,(2*a[0]*a[1])%p)
    def f2add(self,a,b): return ((a[0]+b[0])%self.p,(a[1]+b[1])%self.p)
    def f2sub(self,a,b): return ((a[0]-b[0])%self.p,(a[1]-b[1])%self.p)
    def f2neg(self,a): return ((-a[0])%self.p,(-a[1])%self.p)
    def f2conj(self,a): return (a[0],(-a[1])%self.p)
    def f2mulr(self,a,s): return (a[0]*s%self.p,a[1]*s%self.p)
    def f2pow(self,a,n):
        o=ONE2
        while n:
            if n&1: o=self.f2mul(o,a)
            a=self.f2sqr(a); n>>=1
        return o
    # ---- Fp6
    def f6mul(self,a,b):
        xi=self.xi
        t0=self.f2mul(a[0],b[0]); t1=self.f2mul(a[1],b[1]); t2=self.f2mul(a[2],b[2])
        c0=self.f2add(t0,self.f2mul(xi,self.f2sub(self.f2mul(self.f2add(a[1],a[2]),self.f2add(b[1],b[2])),self.f2add(t1,t2))))
        c1=self.f2add(self.f2sub(self.f2sub(self.f2mul(self.f2add(a[0],a[1]),self.f2add(b[0],b[1])),t0),t1),self.f2mul(xi,t2))
        c2=self.f2add(self.f2sub(self.f2sub(self.f2mul(self.f2add(a[0],a[2]),self.f2add(b[0],b[2])),t0),t2),t1)
        return (c0,c1,c2)
    def f6sqr(self,a):
        xi=self.xi
        s0=self.f2sqr(a[0]); s1=self.f2sqr(a[1]); s2=self.f2sqr(a[2])
        p01=self.f2mul(a[0],a[1]); p02=self.f2mul(a[0],a[2]); p12=self.f2mul(a[1],a[2])
        return (self.f2add(s0,self.f2mul(xi,self.f2add(p12,p12))),
                self.f2add(self.f2add(p01,p01),self.f2mul(xi,s2)),
                self.f2add(s1,self.f2add(p02,p02)))
    def f6sqr_fast(self,a):
        """Fp6 squaring with six Fp2 squarings (no Fp2 multiplication)."""
        xi=self.xi
        t0=self.f2sqr(a[0]); t1=self.f2sqr(a[1]); t2=self.f2sqr(a[2])
        p12=self.f2sub(self.f2sub(self.f2sqr(self.f2add(a[1],a[2])),t1),t2)
        p02=self.f2sub(self.f2sub(self.f2sqr(self.f2add(a[0],a[2])),t0),t2)
        p01=self.f2sub(self.f2sub(self.f2sqr(self.f2add(a[0],a[1])),t0),t1)
        # p12 = 2 a1 a2, p02 = 2 a0 a2, p01 = 2 a0 a1 already
        return (self.f2add(t0,self.f2mul(xi,p12)),
                self.f2add(p01,self.f2mul(xi,t2)),
                self.f2add(t1,p02))
    def f6add(self,a,b): return (self.f2add(a[0],b[0]),self.f2add(a[1],b[1]),self.f2add(a[2],b[2]))
    def f6sub(self,a,b): return (self.f2sub(a[0],b[0]),self.f2sub(a[1],b[1]),self.f2sub(a[2],b[2]))
    def f6neg(self,a): return (self.f2neg(a[0]),self.f2neg(a[1]),self.f2neg(a[2]))
    def f6scale(self,a,s): return (self.f2mul(a[0],s),self.f2mul(a[1],s),self.f2mul(a[2],s))
    def f6muxi(self,a): return (self.f2mul(self.xi,a[2]),a[0],a[1])      # multiply by v
    def f6pow(self,a,n):
        o=self.one6
        while n:
            if n&1: o=self.f6mul(o,a)
            a=self.f6sqr(a); n>>=1
        return o
    # ---- Fp12 = Fp6[w]/(w^2-v)
    def mul(self,a,b):
        t0=self.f6mul(a[0],b[0]); t1=self.f6mul(a[1],b[1])
        t2=self.f6mul(self.f6add(a[0],a[1]),self.f6add(b[0],b[1]))
        return (self.f6add(t0,self.f6muxi(t1)), self.f6sub(self.f6sub(t2,t0),t1))
    def sqr(self,a): return self.mul(a,a)
    def one12(self): return (self.one6,self.zero6)
    def pow(self,a,n):
        o=self.one12()
        while n:
            if n&1: o=self.mul(o,a)
            a=self.mul(a,a); n>>=1
        return o
    def conj12(self,a): return (a[0],self.f6neg(a[1]))
    def cyclo_sqr(self,a):
        """g^2 = (1+2 v c1^2) + ((c0+c1)^2 - (1+v c1^2) - c1^2) w, valid for g^{p^4-p^2+1}=1."""
        c0,c1=a
        t=self.f6sqr(c1)                     # c1^2
        vt=self.f6muxi(t)                    # v c1^2
        c0sq=self.f6add(vt,self.one6)        # 1 + v c1^2 = c0^2
        s=self.f6sqr(self.f6add(c0,c1))
        out0=self.f6add(self.f6add(vt,vt),self.one6)   # 1 + 2 v c1^2
        out1=self.f6sub(self.f6sub(s,c0sq),t)          # 2 c0 c1
        return (out0,out1)
    def cyclo_sqr_fast(self,a):
        """Granger-Scott cyclotomic squaring with the optimal Fp6 squaring."""
        c0,c1=a
        t=self.f6sqr_fast(c1)
        vt=self.f6muxi(t)
        c0sq=self.f6add(vt,self.one6)
        s=self.f6sqr_fast(self.f6add(c0,c1))
        out0=self.f6add(self.f6add(vt,vt),self.one6)
        out1=self.f6sub(self.f6sub(s,c0sq),t)
        return (out0,out1)
    def frob6(self,c):
        return (self.f2conj(c[0]), self.f2mul(self.f2conj(c[1]),self.xi3),
                self.f2mul(self.f2conj(c[2]),self.f2mul(self.xi3,self.xi3)))
    def frob(self,a):
        return (self.frob6(a[0]), self.f6scale(self.frob6(a[1]),self.gamma))
    # ---- class keys
    def keys12(self,a):
        """12 Frobenius keys: {a*zeta^j} U {B*zeta^j}; returns list aligned k=0..11."""
        d0=a[1][0]; g=self.gamma; p=self.p; z=self.zeta
        A=d0[0]; B=(d0[0]*g[0]+d0[1]*g[1])%p
        out=[0]*12
        zj=1
        for j in range(6):
            out[2*j]=(A*zj)%p
            out[2*j+1]=(B*zj)%p
            zj=zj*z%p
        return out
    def kappa(self,a):
        ks=self.keys12(a)
        m=min(ks)
        return m, ks.index(m)

def random_gt(T,r,rng):
    """uniform nonidentity element of the order-r subgroup of GT (r prime, r | Phi12(p))."""
    p=T.p; e=(p**12-1)//r
    while True:
        x=tuple(tuple((rng.randrange(p),rng.randrange(p)) for _ in range(3)) for _ in range(2))
        g=T.pow(x,e)
        if g!=T.one12() and T.pow(g,r)==T.one12(): return g

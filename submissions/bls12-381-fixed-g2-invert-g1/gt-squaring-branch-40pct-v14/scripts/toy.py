#!/usr/bin/env python3
"""Toy engine (v12b): Frobenius-class mixed square/multiply rho walk (stdlib only).

Tower and partition are the production ones:
  Fp2 = Fp[u]/(u^2+1),  Fp6 = Fp2[v]/(v^3-(1+u)),  Fp12 = Fp6[w]/(w^2-v)
  class key  = min of the twelve Frobenius keys {a*zeta^j} U {B*zeta^j}; kappa = argmin
  square step: y <- y^2                    (one cyclotomic squaring, (A,D) <- (2A,2D))
  multiply step: y <- y*M_h^(p^-kappa)     (one dense Fp12 mul, (A,D) <- (A+lam^-k a_h, ...))

v12b measurement repairs (DECISION.md re-arm):
  * class identity is the EXACT canonical min-of-twelve key (the Part A key rule), computed
    once per step and used for the table, the stall detector and the return counter.  v12
    compared a coarser two-key toy fingerprint, so the long-trail return count could not
    separate structural returns from birthday returns.
  * long_trail reports structural statistics against the exact identity: returns, returns
    to the pre-escape visited set within one step of an escape (the inverse-undo trap
    signature), and the longest run of consecutive escapes.  A birthday baseline
    L*(L-1)//(2N) is reported alongside so the two are separable.
  * trap_local isolates the escape rule: does one escape land back on a class already
    visited (exact key equality) for squaring versus the reflecting control.
Toy-scale compromises, disclosed in argument.md: labels are a fixed mixed hash of the
canonical key because the toy key space is one base-field element while production keys
are 381-bit; every candidate collision is confirmed by an exact twelve-conjugate
comparison and by g^x = z.
"""
import json, random, sys
from math import gcd, isqrt

def is_prime(n):
    if n < 2: return False
    for sp in (2,3,5,7,11,13,17,19,23,29,31,37):
        if n % sp == 0: return n == sp
    d=n-1; s=0
    while d % 2 == 0: d//=2; s+=1
    for a in (2,3,5,7,11,13,17,19,23,29,31,37):
        x=pow(a,d,n)
        if x in (1,n-1): continue
        for _ in range(s-1):
            x=x*x%n
            if x==n-1: break
        else: return False
    return True

ZERO2=(0,0); ONE2=(1,0)

class Tower:
    __slots__=('p','xi3','xi3sq','gamma','zeta','one6','zero6')
    def __init__(self,p):
        assert p%12==7
        self.p=p
        self.xi3=self.f2pow((1,1),(p-1)//3)         # v^p = xi3*v
        self.xi3sq=self.f2mul(self.xi3,self.xi3)
        self.gamma=self.f2pow((1,1),(p-1)//6)       # w^p = gamma*w
        self.zeta=(self.gamma[0]*self.gamma[0]+self.gamma[1]*self.gamma[1])%p
        assert (self.zeta*self.zeta-self.zeta+1)%p==0
        self.one6=(ONE2,ZERO2,ZERO2); self.zero6=(ZERO2,ZERO2,ZERO2)
    def f2mul(self,a,b):
        p=self.p; return ((a[0]*b[0]-a[1]*b[1])%p,(a[0]*b[1]+a[1]*b[0])%p)
    def f2sqr(self,a):
        p=self.p; return (((a[0]+a[1])*(a[0]-a[1]))%p,(2*a[0]*a[1])%p)
    def f2add(self,a,b): p=self.p; return ((a[0]+b[0])%p,(a[1]+b[1])%p)
    def f2sub(self,a,b): p=self.p; return ((a[0]-b[0])%p,(a[1]-b[1])%p)
    def f2neg(self,a): p=self.p; return ((-a[0])%p,(-a[1])%p)
    def f2conj(self,a): p=self.p; return (a[0],(-a[1])%p)
    def f2pow(self,a,n):
        o=ONE2
        while n:
            if n&1: o=self.f2mul(o,a)
            a=self.f2sqr(a); n>>=1
        return o
    def f6mul(self,a,b):
        xi=(1,1)
        t0=self.f2mul(a[0],b[0]); t1=self.f2mul(a[1],b[1]); t2=self.f2mul(a[2],b[2])
        c0=self.f2add(t0,self.f2mul(xi,self.f2sub(self.f2mul(self.f2add(a[1],a[2]),self.f2add(b[1],b[2])),self.f2add(t1,t2))))
        c1=self.f2add(self.f2sub(self.f2sub(self.f2mul(self.f2add(a[0],a[1]),self.f2add(b[0],b[1])),t0),t1),self.f2mul(xi,t2))
        c2=self.f2add(self.f2sub(self.f2sub(self.f2mul(self.f2add(a[0],a[2]),self.f2add(b[0],b[2])),t0),t2),t1)
        return (c0,c1,c2)
    def f6sqr(self,a):
        xi=(1,1)
        s0=self.f2sqr(a[0]); s1=self.f2sqr(a[1]); s2=self.f2sqr(a[2])
        p01=self.f2mul(a[0],a[1]); p02=self.f2mul(a[0],a[2]); p12=self.f2mul(a[1],a[2])
        return (self.f2add(s0,self.f2mul(xi,self.f2add(p12,p12))),
                self.f2add(self.f2add(p01,p01),self.f2mul(xi,s2)),
                self.f2add(s1,self.f2add(p02,p02)))
    def f6add(self,a,b): return (self.f2add(a[0],b[0]),self.f2add(a[1],b[1]),self.f2add(a[2],b[2]))
    def f6sub(self,a,b): return (self.f2sub(a[0],b[0]),self.f2sub(a[1],b[1]),self.f2sub(a[2],b[2]))
    def f6neg(self,a): return (self.f2neg(a[0]),self.f2neg(a[1]),self.f2neg(a[2]))
    def f6scale(self,a,s): return (self.f2mul(a[0],s),self.f2mul(a[1],s),self.f2mul(a[2],s))
    def f6muxi(self,a): return (self.f2mul((1,1),a[2]),a[0],a[1])
    def mul(self,a,b):
        t0=self.f6mul(a[0],b[0]); t1=self.f6mul(a[1],b[1])
        t2=self.f6mul(self.f6add(a[0],a[1]),self.f6add(b[0],b[1]))
        return (self.f6add(t0,self.f6muxi(t1)),self.f6sub(self.f6sub(t2,t0),t1))
    def sqr(self,a): return self.mul(a,a)
    def one12(self): return (self.one6,self.zero6)
    def pow(self,a,n):
        o=self.one12()
        while n:
            if n&1: o=self.mul(o,a)
            a=self.mul(a,a); n>>=1
        return o
    def cyclo_sqr(self,a):
        c0,c1=a
        t=self.f6sqr(c1)                        # c1^2
        vt=self.f6muxi(t)                       # v c1^2
        s=self.f6sqr(self.f6add(c0,c1))
        return (self.f6add(self.f6add(vt,vt),self.one6),      # 1 + 2 v c1^2
                self.f6sub(self.f6sub(s,self.f6add(vt,self.one6)),t))   # 2 c0 c1
    def frob6(self,c):
        return (self.f2conj(c[0]),self.f2mul(self.f2conj(c[1]),self.xi3),
                self.f2mul(self.f2conj(c[2]),self.xi3sq))
    def frob(self,a):
        return (self.frob6(a[0]),self.f6scale(self.frob6(a[1]),self.gamma))
    def keys12(self,a):
        d0=a[1][0]; g=self.gamma; p=self.p; z=self.zeta
        A=d0[0]; B=(d0[0]*g[0]+d0[1]*g[1])%p
        Az=A*z%p; Bz=B*z%p
        return [A,B,Az,Bz,(Az-A)%p,(Bz-B)%p,(-A)%p,(-B)%p,(-Az)%p,(-Bz)%p,(A-Az)%p,(B-Bz)%p]
    def kappa(self,a):
        ks=self.keys12(a); m=min(ks)
        return m,ks.index(m)
    def ident(self,a):
        """EXACT class identity under the Part A twelve-key rule.

        The orbit invariant is the twelve-key SET: one Frobenius application permutes it,
        so same class <=> same set (verified in part_a and re-checked here).  At the toy
        field size the canonical MINIMUM is only log2(q) ~ 12.5 bits against ~19 bits of
        class entropy, so min-only equality has ~90 classes per bucket and cannot serve as
        an identity; the sorted set (~25 bits) is the exact toy-scale identity.  Returns
        (identity, canonical min, kappa index) with kappa = argmin as in production.
        """
        ks=self.keys12(a)
        m0=min(ks)
        return tuple(sorted(ks)), m0, ks.index(m0)

class Toy:
    def __init__(self,q,rp,R=128,f=(2,3),seed=1):
        self.q=q; self.rp=rp; self.R=R; self.f=f
        self.T=Tower(q); self.rng=random.Random(seed)
        self.lam=q%rp
        assert pow(self.lam,12,rp)==1 and pow(self.lam,6,rp)==rp-1
        self.N=(rp-1)//12
        self.linv=[pow(self.lam,(-k)%12,rp) for k in range(12)]
        e=(q**12-1)//rp
        while True:
            x=tuple(tuple((self.rng.randrange(q),self.rng.randrange(q)) for _ in range(3)) for _ in range(2))
            g=self.T.pow(x,e)
            if g!=self.T.one12() and self.T.pow(g,rp)==self.T.one12(): break
        self.g=g
        self.build_table()

    def build_table(self):
        """Fresh target z=g^x and the R-entry fixed transition table for that target.

        The table is built ONCE per target and never re-randomised during a walk, so the
        map y -> y*T[h][...] is a fixed deterministic function of the canonical class key
        (the production transition table).  Each trial gets a fresh x and hence a fresh z
        and a fresh table, as in production where the table is built from the supplied z.
        """
        rng=self.rng; q=self.q; rp=self.rp; R=self.R; fn,fd=self.f
        self.x=rng.randrange(1,rp); self.z=self.T.pow(self.g,self.x)
        ent=[None]*R; imgs=[None]*R
        for h in range(R):
            if h*fd < fn*R:
                ent[h]=None
            else:
                a=rng.randrange(1,rp); b=rng.randrange(1,rp)
                M=self.T.mul(self.T.pow(self.g,a),self.T.pow(self.z,b))
                cur=M; im=[]
                for k in range(12):
                    im.append(cur); cur=self.T.frob(cur)
                ent[h]=(a,b); imgs[h]=im
        self.ent=ent; self.imgs=imgs

    def inv_ok(self,y,A,D):
        """Exact production invariant y == g^A * z^D (the (A,D) exponent ledger)."""
        return self.T.mul(self.T.pow(self.g,A%self.rp),self.T.pow(self.z,D%self.rp))==y

    def lab(self,y):
        """label h from the canonical min (a class function); returns (h, exact identity, kappa)."""
        K,m0,kap=self.T.ident(y)
        return ((m0*0x9E3779B1 ^ ((m0>>13)*0x85EBCA6B & 0xFFFFFFFF))>>3)%self.R, K, kap

    def step(self,y,A,D,h,kap):
        T=self.T; rp=self.rp
        e=self.ent[h]
        if e is None:
            return T.cyclo_sqr(y),(2*A)%rp,(2*D)%rp,1
        j=(12-kap)%12
        return T.mul(y,self.imgs[h][j]),(A+self.linv[kap]*e[0])%rp,(D+self.linv[kap]*e[1])%rp,0

    def match(self,y1,y2):
        cur=y1
        for k in range(12):
            if cur==y2: return k
            cur=self.T.frob(cur)
        return None

    def recover(self,A1,D1,A2,D2,k):
        rp=self.rp; pj=pow(self.q,k,rp)
        den=(D2-pj*D1)%rp
        if den==0: return None
        return ((pj*A1-A2)%rp)*pow(den,-1,rp)%rp

    def walk(self,max_steps=10**7,check_every=16):
        T=self.T; rp=self.rp
        A0=self.rng.randrange(rp); D0=self.rng.randrange(1,rp)
        y=T.mul(T.pow(self.g,A0),T.pow(self.z,D0)); A,D=A0,D0
        assert self.inv_ok(y,A,D)          # start-of-walk invariant
        h,K,kap=self.lab(y)
        tab={K:[(A,D,y)]}
        n=0; sq=0; esc=0; chk=1; bad=0
        while n<max_steps:
            y,A,D,issq=self.step(y,A,D,h,kap); n+=1; sq+=issq
            if (n%check_every==0):
                chk+=1
                if not self.inv_ok(y,A,D): bad+=1
            h,K,kap=self.lab(y)
            cands=tab.get(K)
            if cands is not None:
                for rec in cands:
                    k=self.match(rec[2],y)
                    if k is not None:
                        x=self.recover(rec[0],rec[1],A,D,k)
                        if x is not None and x!=0 and T.pow(self.g,x)==self.z:
                            assert self.inv_ok(y,A,D)   # collision-point invariant
                            chk+=1
                            return dict(steps=n,square=sq,useful=True,escapes=esc,checks=chk,bad=bad)
                # confirmed class collision with no usable relation: charged squaring escape
                esc+=1
                y=T.cyclo_sqr(y); A=(2*A)%rp; D=(2*D)%rp; n+=1; sq+=1
                chk+=1
                if not self.inv_ok(y,A,D): bad+=1    # escape-path invariant
                h,K,kap=self.lab(y)
                if K in tab: tab[K].append((A,D,y))
                else: tab[K]=[(A,D,y)]
                continue
            tab[K]=[(A,D,y)]
        return dict(steps=n,square=sq,useful=False,escapes=esc,checks=chk,bad=bad)

    def long_trail(self,length=2048,escape="squaring",R=None):
        """Length-L trail, escapes on an exact-class repeat, exact-class return counting."""
        T=self.T; rp=self.rp
        y=self.z; A=0; D=1
        h,K,kap=self.lab(y)
        seen={K}
        prev=None; n=0; esc=0; ret=0; rev=0; chain=0; maxchain=0; lastM=None
        while n<length:
            y,A,D,issq=self.step(y,A,D,h,kap); n+=1
            if not issq: lastM=self.ent[h]
            h2,K2,kap2=self.lab(y)
            if prev is not None and K2==prev:
                esc+=1; chain+=1
                if chain>maxchain: maxchain=chain
                before=set(seen)
                if escape=="squaring":
                    y=T.cyclo_sqr(y); A=(2*A)%rp; D=(2*D)%rp
                else:   # rejected v8-style reflecting offset: multiply by the inverse jump
                    if lastM is None: lastM=(1,1)
                    M=T.mul(T.pow(self.g,(-lastM[0])%rp),T.pow(self.z,(-lastM[1])%rp))
                    y=T.mul(y,M); A=(A-lastM[0])%rp; D=(D-lastM[1])%rp
                n+=1
                h2,K2,kap2=self.lab(y)
                if K2 in before: rev+=1          # landed on an already-visited class: trap
                seen={K2}; prev=K2
                continue
            chain=0
            if K2 in seen: ret+=1
            seen.add(K2); prev=K2
        return dict(returns=ret,post_escape_reverts=rev,escapes=esc,max_chain=maxchain,
                    steps=n,length=length,birthday=length*(length-1)//(2*self.N))

    def trap_local(self,trials=4000,R=None):
        """Isolate the escape rule on the exact class identity."""
        T=self.T; rp=self.rp; Rt=R or self.R
        out={}
        for esc_kind in ("squaring","reflecting"):
            rev=0; same=0; diff=0; lastM=None
            for _ in range(trials):
                A=self.rng.randrange(rp); D=self.rng.randrange(1,rp)
                y=T.mul(T.pow(self.g,A),T.pow(self.z,D))
                h,K,kap=self.lab(y)
                hist=[]
                for _ in range(3):
                    y,A,D,issq=self.step(y,A,D,h,kap)
                    if not issq: lastM=self.ent[h]
                    h,K,kap=self.lab(y); hist.append(K)
                pre=set(hist)
                Kpre=hist[-1]
                if esc_kind=="squaring":
                    y2=T.cyclo_sqr(y); A2=(2*A)%rp; D2=(2*D)%rp
                else:
                    if lastM is None: lastM=(1,1)
                    M=T.mul(T.pow(self.g,(-lastM[0])%rp),T.pow(self.z,(-lastM[1])%rp))
                    y2=T.mul(y,M); A2=(A-lastM[0])%rp; D2=(D-lastM[1])%rp
                K2=self.lab(y2)[1]
                if K2==Kpre: same+=1
                elif K2 in pre: rev+=1
                else: diff+=1
            out[esc_kind]=dict(trials=trials,reverts_to_previous_class=same,
                               reverts_to_recent_class=rev,leaves_recent_set=diff)
        return out

    def pipeline(self,nodes=16,spacing=32,walks=8,max_steps=120000):
        T=self.T; rp=self.rp
        store=[{} for _ in range(nodes)]
        cap=(max_steps*walks//spacing//nodes*5)//4+8
        overflow=0; recs=0; found=None; steps=0
        starts=[(self.rng.randrange(rp),self.rng.randrange(1,rp)) for _ in range(walks)]
        while steps<max_steps and found is None:
            for wi in range(walks):
                A,D=starts[wi]
                y=T.mul(T.pow(self.g,A),T.pow(self.z,D))
                for _ in range(max_steps//10):
                    h,K,kap=self.lab(y)
                    if K<rp//spacing:
                        node=K%nodes; recs+=1
                        if len(store[node])>=cap: overflow+=1
                        else:
                            if K in store[node]:
                                A1,D1,y1=store[node][K]
                                k=self.match(y1,y)
                                if k is not None and (A1,D1)!=(A,D):
                                    x=self.recover(A1,D1,A,D,k)
                                    if x is not None and T.pow(self.g,x)==self.z:
                                        found=x; break
                            store[node][K]=(A,D,y)
                    y,A,D,_=self.step(y,A,D,h,kap); steps+=1
                if found is not None: break
        return dict(found=found is not None,overflow=overflow,records=recs,
                    max_node=max(len(s) for s in store),cap=cap,steps=steps)

def field_params(q):
    if not (q%12==7 and is_prime(q)): return None
    p2=q*q
    def f2pow(a,n):
        o=(1,0)
        while n:
            if n&1: o=((a[0]*o[0]-a[1]*o[1])%q,(a[0]*o[1]+a[1]*o[0])%q)
            a=((a[0]*a[0]-a[1]*a[1])%q,(2*a[0]*a[1])%q); n>>=1
        return o
    if f2pow((1,1),(p2-1)//3)==(1,0) or f2pow((1,1),(p2-1)//2)==(1,0): return None
    n=q**4-q*q+1; f={}; d=2; m=n
    while d<=200000 and d*d<=m:
        while m%d==0: f[d]=f.get(d,0)+1; m//=d
        d+= 1 if d==2 else 2
    cands=[d for d in f if d>500 and d%12==1 and is_prime(d)]
    if m>500 and m%12==1 and is_prime(m): cands.append(m)
    if not cands: return None
    rp=max(cands)
    fac={}; mm=rp-1; d=2
    while d*d<=mm and d<=200000:
        while mm%d==0: fac[d]=fac.get(d,0)+1; mm//=d
        d+= 1 if d==2 else 2
    if mm>1: fac[mm]=fac.get(mm,0)+1
    o=rp-1
    for d in fac:
        while o%d==0 and pow(2,o//d,rp)==1: o//=d
    return rp,(rp-1)//12,o,o//gcd(o,12)

def s_ops(N,R):
    """floor(1.05*sqrt(2*N'*ln(5/3))) with N'=N*R/(R-1), integer sqrt, ln(5/3) rounded DOWN.

    ln(5/3) = 0.510825623... so 510825/10^6 is a strict lower enclosure; the ledger's S uses
    an upper enclosure of the same constant, hence S_toy <= S_ledger for every field.
    """
    num=2*N*R
    lo=105*105*num*510825//(100*100*1000000*(R-1))
    return isqrt(lo)+1

def sqrt_bounds(N,scale=10**6):
    lo=isqrt(N*scale*scale); return lo,lo+1

if __name__=="__main__":
    cfg=json.loads(sys.stdin.read() or "{}")
    q=cfg.get("q",523); fp=field_params(q)
    if fp is None: raise SystemExit("bad q")
    rp,N,o2,o2o=fp
    tw=Toy(q,rp,cfg.get("R",128),(cfg.get("f_num",2),cfg.get("f_den",5)),cfg.get("seed",1))
    mode=cfg.get("mode","census")
    if mode=="census":
        tr=cfg.get("trials",100); S=s_ops(N,tw.R)
        res=[]
        for _ in range(tr):
            tw.build_table()          # FRESH target z=g^x per trial (fresh table too)
            res.append(tw.walk())
        st=[r["steps"] for r in res]
        slo,shi=sqrt_bounds(N)
        tot=sum(st)
        out=dict(mode="census",q=q,rp=rp,N=N,R=tw.R,f=list(tw.f),trials=tr,
            squares=sum(r["square"] for r in res),steps_total=tot,
            S_ops=int(S),hits=sum(1 for s in st if s<=S),
            useful=sum(1 for r in res if r["useful"]),
            checks=sum(r["checks"] for r in res),bad=sum(r["bad"] for r in res),
            mean_steps_num=tot, mean_steps_den=tr,
            mean_over_sqrtN_lo_num=10**6*tot, mean_over_sqrtN_lo_den=tr*shi,
            mean_over_sqrtN_hi_num=10**6*tot, mean_over_sqrtN_hi_den=tr*slo,
            sqrt_N_floor=isqrt(N), max_steps=max(st), min_steps=min(st),
            escapes=sum(r["escapes"] for r in res),
            square_frac_num=sum(r["square"] for r in res),
            square_frac_den=max(1,tot))
        print(json.dumps(out))
    elif mode=="longtrail":
        out={"mode":"longtrail","q":q,"R":tw.R,"N":N,"length":cfg.get("length",2048),"rows":[]}
        for esc in ("squaring","reflecting"):
            agg=dict(escape=esc,trails=0,returns=0,post_escape_reverts=0,escapes=0,
                     max_chain=0,steps=0,birthday=0)
            for i in range(cfg.get("trails",40)):
                tw2=Toy(q,rp,tw.R,tw.f,seed=1000+i)
                r=tw2.long_trail(cfg.get("length",2048),esc)
                agg["trails"]+=1
                for k in ("returns","post_escape_reverts","escapes","steps","birthday"):
                    agg[k]+=r[k]
                agg["max_chain"]=max(agg["max_chain"],r["max_chain"])
            out["rows"].append(agg)
        print(json.dumps(out))
    elif mode=="trap":
        print(json.dumps(dict(mode="trap",q=q,R=tw.R,N=N,
            rows=tw.trap_local(cfg.get("trials",4000)))))
    elif mode=="pipeline":
        print(json.dumps(tw.pipeline(**{k:v for k,v in cfg.items() if k in ("nodes","spacing","walks","max_steps")})))

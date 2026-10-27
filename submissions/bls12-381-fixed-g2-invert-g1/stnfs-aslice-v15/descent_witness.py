"""Production-field boot and special-q lattice witnesses; stdlib only.
All arithmetic/membership/norm checks are exact. Decimal LLL proposes a basis;
independent exact Gram-Schmidt verifies its LLL conditions and determinant.
"""
import os,sys,json,math,time,random
from fractions import Fraction as F
from decimal import Decimal as D, localcontext
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,os.path.join(HERE,'files') if os.path.isdir(os.path.join(HERE,'files')) else HERE)
import common12 as C
P=4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787
R=0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
U=C.U
Z=[0]*6; ONE=([1]+[0]*5,Z[:])
def add(a,b): return [(x+y)%P for x,y in zip(a,b)]
def neg(a): return [(-x)%P for x in a]
def mul(a,b): return [x%P for x in C.mul(a,b)]
def pow6(a,n):
    v=ONE[0]
    while n:
        if n&1:v=mul(v,a)
        a=mul(a,a);n>>=1
    return v
def mul12(a,b):
    ac=mul(a[0],b[0]);bd=mul(a[1],b[1]);cross=mul(add(a[0],a[1]),add(b[0],b[1]))
    return add(ac,neg(mul([0,U,0,0,0,0],bd))),add(cross,neg(add(ac,bd)))
def pow12(a,n):
    v=ONE
    while n:
        if n&1:v=mul12(v,a)
        a=mul12(a,a);n>>=1
    return v
def uniform_gt(seed):
    rng=random.Random(seed)
    while True:
        x=([rng.randrange(P) for _ in range(6)],[rng.randrange(P) for _ in range(6)])
        if any(x[0]+x[1]):break
    # Surjective homomorphism Fp12* -> mu_r; every fibre has equal cardinality.
    return pow12(x,(P**12-1)//R)
def boot_basis(z):
    if not any(z[1]):raise ValueError('subfield/identity: solve separately')
    ratio=mul(neg(z[0]),pow6(z[1],P**6-2))
    return [[P if i==j else 0 for j in range(12)] for i in range(6)]+[mul(ratio,[0]*j+[1])+[int(i==j) for i in range(6)] for j in range(6)]
def gs_exact(B):
    stars=[];mu=[];norm=[]
    for v in B:
        row=[]; s=list(map(F,v))
        for j,b in enumerate(stars):
            m=sum(x*y for x,y in zip(v,b))/norm[j];row.append(m)
            s=[x-m*y for x,y in zip(s,b)]
        stars.append(s);mu.append(row);norm.append(sum(x*x for x in s))
    return mu,norm
def lll(B):
    # Incremental GSO updates; precision high enough for 381-bit q-ary input.
    B=[v[:] for v in B];n=len(B);swaps=0
    with localcontext() as ctx:
        ctx.prec=300
        mu=[];norm=[];star=[]
        for v in B:
            s=list(map(D,v)); row=[]
            for j,b in enumerate(star):
                x=sum(D(a)*c for a,c in zip(v,b))/norm[j];row.append(x)
                s=[a-x*c for a,c in zip(s,b)]
            star.append(s);mu.append(row);norm.append(sum(x*x for x in s))
        k=1
        while k<n:
            for j in range(k-1,-1,-1):
                q=int(mu[k][j].to_integral_value(rounding='ROUND_HALF_EVEN'))
                if q:
                    B[k]=[x-q*y for x,y in zip(B[k],B[j])]
                    for i in range(j):mu[k][i]-=q*mu[j][i]
                    mu[k][j]-=q
            m=mu[k][k-1]
            if norm[k]>=(D('0.99')-m*m)*norm[k-1]:k+=1;continue
            A=norm[k]+m*m*norm[k-1];newm=m*norm[k-1]/A
            norm[k]=norm[k]*norm[k-1]/A;norm[k-1]=A;mu[k][k-1]=newm
            B[k],B[k-1]=B[k-1],B[k]
            for j in range(k-1):mu[k][j],mu[k-1][j]=mu[k-1][j],mu[k][j]
            for i in range(k+1,n):
                t=mu[i][k];mu[i][k]=mu[i][k-1]-m*t;mu[i][k-1]=t+newm*mu[i][k]
            k=max(1,k-1);swaps+=1
    return B,swaps

def certify_boot(z,B):
    for v in B:
        assert mul(v[:6],z[1])==neg(mul(v[6:],z[0]))
    mu,norm=gs_exact(B)
    for i in range(12):
        assert all(abs(x)<=F(1,2) for x in mu[i])
        if i:assert norm[i]>=(F(98,100)-mu[i][-1]**2)*norm[i-1]
    det2=math.prod(norm);assert det2==P**12
    return float(sum(norm)),mu,norm

def boot_stats(z,B,seed,points=64,scales=(2,4,8,16)):
    cert,mu,norm=certify_boot(z,B);rng=random.Random(seed)
    out=[]
    for scale in scales:
        vals=[];radii=[]
        for _ in range(points):
            c=[rng.randint(-scale,scale) for _ in B]
            if not any(c):c[0]=1
            v=[sum(c[i]*B[i][j] for i in range(12)) for j in range(12)]
            assert mul(v[:6],z[1])==neg(mul(v[6:],z[0]))
            ng=C.norm([a+b for a,b in zip(C.mul(v[:6],v[:6]),C.mul([0,U],C.mul(v[6:],v[6:])) )])
            assert ng>0
            vals.append(math.log2(ng));radii.append(math.log2(sum(x*x for x in v))/2)
        out.append({'scale':scale,'points':points,'mean_bits':sum(vals)/points,'max_bits':max(vals),'min_bits':min(vals),'radius_log2_max':max(radii)})
    return out

def main():
    inp=json.load(sys.stdin);seed=inp.get('seed',13001);t=time.time()
    z=uniform_gt(seed); print('GT generated',time.time()-t,file=sys.stderr,flush=True)
    B,swaps=lll(boot_basis(z));print('LLL done',swaps,time.time()-t,file=sys.stderr,flush=True)
    stats=boot_stats(z,B,seed+1,points=inp.get('points',16))
    print(json.dumps({'seed':seed,'z':z,'basis':B,'swaps':swaps,'stats':stats,'elapsed':time.time()-t}))
if __name__=='__main__':main()

# Finite-field polynomial tools for actual relative-degree-one special-q ideals.
def trim(a):
    while a and a[-1]==0:a.pop()
    return a
def psub(a,b,q):return trim([((a[i] if i<len(a) else 0)-(b[i] if i<len(b) else 0))%q for i in range(max(len(a),len(b)))])
def pmod(a,b,q):
    a=trim([x%q for x in a]);iv=pow(b[-1],-1,q)
    while len(a)>=len(b):
        v=a[-1]*iv%q;k=len(a)-len(b)
        for i,x in enumerate(b):a[k+i]=(a[k+i]-v*x)%q
        trim(a)
    return a
def pmul(a,b,f,q):
    v=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):v[i+j]=(v[i+j]+x*y)%q
    return pmod(v,f,q)
def ppow(a,n,f,q):
    v=[1]
    while n:
        if n&1:v=pmul(v,a,f,q)
        a=pmul(a,a,f,q);n>>=1
    return v
def pgcd(a,b,q):
    while b:a,b=b,pmod(a,b,q)
    if not a:return []
    iv=pow(a[-1],-1,q);return [x*iv%q for x in a]
def roots(f,q,rng):
    f=trim([x%q for x in f]);f=pgcd(f,psub(ppow([0,1],q,f,q),[0,1],q),q)
    def split(g):
        if len(g)<=1:return []
        if len(g)==2:return [-g[0]*pow(g[1],-1,q)%q]
        while True:
            a=[rng.randrange(q) for _ in range(len(g)-1)]
            h=pgcd(g,psub(ppow(a,(q-1)//2,g,q),[1],q),q)
            if 1<len(h)<len(g):
                # monic exact quotient
                a=g[:];v=[0]*(len(g)-len(h)+1)
                for k in range(len(v)-1,-1,-1):
                    v[k]=a[k+len(h)-1]
                    for j,x in enumerate(h):a[k+j]=(a[k+j]-v[k]*x)%q
                assert not any(a)
                return split(h)+split(trim(v))
    return sorted(split(f))
def probable_prime(n):
    if n%2==0:return False
    d=n-1;s=0
    while d%2==0:d//=2;s+=1
    for a in (2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53):
        if n==a:return True
        x=pow(a,d,n)
        if x in (1,n-1):continue
        for _ in range(s-1):
            x=x*x%n
            if x==n-1:break
        else:return False
    return True

def specialq_witness(bits,side,seed,points=32,translated=False,scales=(2,3,4,5,6,8,10,12,16,24,32,48,64,96,128)):
    rng=random.Random(seed)
    while True:
        q=rng.getrandbits(bits)|(1<<(bits-1))|1
        if not probable_prime(q):continue
        ss=roots(C.H,q,rng)
        opts=[]
        for s in ss:
            if side=='g':ff=[U*s,0,1]
            else:ff=[1+s-s*s-s**3,0,s**5,0,0,0,2*s**3,0,0,0,-2*s,0,1]
            opts.extend((s,t) for t in roots(ff,q,rng))
        if opts:break
    s,t=rng.choice(opts)
    w=[pow(s,j,q) for j in range(6)]+[(-t*pow(s,j,q))%q for j in range(6)]
    B=[[q]+[0]*11]+[[(-w[j])%q]+[int(k==j) for k in range(1,12)] for j in range(1,12)]
    B,swaps=lll(B)
    assert all(sum(a*b for a,b in zip(v,w))%q==0 for v in B)
    mu,nn=gs_exact(B);assert math.prod(nn)==q*q
    rows=[]
    for scale in scales:
        nf=[];ng=[];cf=[];cg=[];cb=[]
        # Exactly (2*scale+1)^12 coefficient tuples: S = log2(cardinality).
        S=12*math.log2(2*scale+1);c=(math.log2(q)+S)/12
        for _ in range(points):
            k=[rng.randint(-scale,scale)+(rng.randint(-scale,scale) if translated else 0) for _ in range(12)]
            if not any(k):k[0]=1
            v=[sum(k[i]*B[i][j] for i in range(12)) for j in range(12)]
            a,b=C.norms(v[:6],v[6:]);assert (a if side=='f' else b)%q==0
            af=math.log2(a);ag=math.log2(b)
            nf.append(af);ng.append(ag);cf.append(af-72*c);cg.append(ag-12*c-6*math.log2(U));cb.append(math.log2(max(abs(x) for x in v))-c)
        
        import descent14 as DP
        ratios=[]
        for n in (69,79,99,129,159,199,229):
            if n>=math.log2(q):continue
            df=math.log2(q) if side=='f' else 0; dg=math.log2(q) if side=='g' else 0
            expected=DP.proba(72*c+64-df,5,n)*DP.proba(12*c+6*math.log2(U)+12-dg,1,n)
            actual=[DP.proba(a-df,5,n)*DP.proba(b-dg,1,n) for a,b in zip(nf,ng)]
            ratios.append(sum(actual)/points/expected)
        rows.append({'joint_model_ratio_min':min(ratios),'scale':scale,'S':S,'c':c,'points':points,'cf_max':max(cf),'cf_mean':sum(cf)/points,'cg_max':max(cg),'cg_mean':sum(cg)/points,'coord_excess_max':max(cb)})
    return {'bits':bits,'side':side,'seed':seed,'q':q,'s':s,'t':t,'basis':B,'stats':rows,'swaps':swaps}

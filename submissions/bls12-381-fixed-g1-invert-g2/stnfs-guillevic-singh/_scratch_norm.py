import cmath, math, random
U_PARAM = -0xD201000000010000

def roots6():
    d_rad = math.sqrt(23.0/108.0)
    u_c = (0.5+d_rad)**(1/3.); v_c=(0.5-d_rad)**(1/3.)
    z1=-(u_c+v_c)
    z2=complex(0.5*(u_c+v_c), 0.5*math.sqrt(3.0)*(u_c-v_c))
    z3=complex(0.5*(u_c+v_c), -0.5*math.sqrt(3.0)*(u_c-v_c))
    out=[]
    for z in (z1,z2,z3):
        r=cmath.sqrt(z)
        out += [r,-r]
    return out

def norms(a,b):
    lf=0.0; lg=0.0
    for th in roots6():
        aw=sum(c*th**i for i,c in enumerate(a))
        bw=sum(c*th**i for i,c in enumerate(b))
        F=aw**12-2*th*aw**10*bw**2+2*th**3*aw**6*bw**6+th**5*aw**2*bw**10+(th**2-1)*bw**12
        G=aw**2-U_PARAM*th*bw**2
        lf+=math.log(abs(F)); lg+=math.log(abs(G))
    return lf,lg

def buggy(a,b):
    d_rad = math.sqrt(23.0/108.0)
    u_c = (0.5+d_rad)**(1/3.); v_c=(0.5-d_rad)**(1/3.)
    z1=-(u_c+v_c)
    z2=complex(0.5*(u_c+v_c), 0.5*math.sqrt(3.0)*(u_c-v_c))
    z3=complex(0.5*(u_c+v_c), -0.5*math.sqrt(3.0)*(u_c-v_c))
    lf=0.;lg=0.
    for z in (z1,z2,z3):
        th=cmath.sqrt(z)
        aw=sum(c*th**i for i,c in enumerate(a))
        bw=sum(c*th**i for i,c in enumerate(b))
        F=aw**12-2*th*aw**10*bw**2+2*th**3*aw**6*bw**6+th**5*aw**2*bw**10+(th**2-1)*bw**12
        G=aw**2-U_PARAM*th*bw**2
        lf+=2*math.log(abs(F)); lg+=2*math.log(abs(G))
    return lf,lg

# reviewer example
print("reviewer example buggy", buggy([1,1,0,0,0,0],[0]*6), "expect 0,0 (N_g=1, N_f=h(-1)^12? )")
print("reviewer example correct", norms([1,1,0,0,0,0],[0]*6))

rng=random.Random(7); A=686
sf=sg=bf=bg=0.0; n=0
df=[];dg=[]
while n<2000:
    a=[rng.randint(-A,A) for _ in range(6)]
    b=[rng.randint(-A,A) for _ in range(5)]+[rng.randint(0,A)]
    g=0
    for c in a+b:
        g=math.gcd(g,abs(c))
    if g!=1: continue
    lf,lg=norms(a,b)
    bf_,bg_=buggy(a,b)
    sf+=lf; sg+=lg; bf+=bf_; bg+=bg_
    df.append((lf-bf_)/math.log(2)); dg.append((lg-bg_)/math.log(2))
    n+=1
print("avg log2 Nf correct %.3f buggy %.3f"%(sf/n/math.log(2), bf/n/math.log(2)))
print("avg log2 Ng correct %.3f buggy %.3f"%(sg/n/math.log(2), bg/n/math.log(2)))
print("mean diff f %.3f bits, g %.3f bits"%(sum(df)/n, sum(dg)/n))
print("ratio Nf correct/buggy in bits range", min(df), max(df))

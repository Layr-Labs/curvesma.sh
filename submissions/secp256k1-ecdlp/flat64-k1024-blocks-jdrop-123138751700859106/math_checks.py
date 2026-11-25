"""Finite field and selected-X coefficient facts for the supplied proofs.
The raw Eisenstein inventory follows the denominator-clearing check accompanying
submission283cdacf. These checks establish finite algebra, not the production
pole-degree, collision, or performance arguments.
"""
def check(K=None,full=False):
    import itertools,math,time
    from fractions import Fraction as F
    from types import SimpleNamespace
    started=time.monotonic()
    if K is None:
        K=SimpleNamespace(p=int('fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f',16),r=int('fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141',16),LAM_R=int('5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72',16),NODES=1<<58,N_CLUST=1<<48,M_G=1024,LEAVES=6,MX=6)
    def rem(a,b):
        while a and a.bit_length()>=b.bit_length():a^=b<<(a.bit_length()-b.bit_length())
        return a
    def gcd(a,b):
        while b:a,b=b,rem(a,b)
        return a
    def polymul(a,b):
        z=0
        while b:
            if b&1:z^=a
            a<<=1;b>>=1
        return z
    def irreducible(f):
        powers=[2]
        for _ in range(10):powers.append(rem(polymul(powers[-1],powers[-1]),f))
        return powers[10]==2 and gcd(f,powers[5]^2)==1 and gcd(f,powers[2]^2)==1
    assert irreducible(0x409) and not irreducible(0x400)
    def gf(a,b):return rem(polymul(a,b),0x409)
    def rank(vs):
        piv={}
        for a in vs:
            while a:
                i=a.bit_length()-1
                if i in piv:a^=piv[i]
                else:piv[i]=a;break
        return len(piv)
    assert all(rank([gf(a,1<<i)&15 for i in range(10)])==4 for a in range(1,1024))
    assert rank([gf(1,1<<i)&7 for i in range(10)])==3
    # Low4(gf(v,m)) has64preimages whenv!=0; v=0 is constant.
    # Count every(delta_u,delta_v) residue with the exact linear-map rank.
    counts=[]
    enumerated=0
    for v in range(1024):
        hist=[64]*16 if v else[1024]+[0]*15
        # Rank4 was checked for every nonzero v, hence each fiber has2^6
        # members. Enumerate four nonzero maps and zero as a separate check;
        # optional full mode enumerates all multiplication-table entries.
        if full or v in(0,1,2,512,1023):
            observed=[0]*16
            for m in range(1024):observed[gf(v,m)&15]+=1
            assert observed==hist;enumerated+=1024
        counts.extend(hist)
    mean=F(sum(a*(a-1) for a in counts),16384)
    avg=(K.LEAVES*K.MX)**2*(1+2*mean/K.M_G**2)
    # Removing identical fullcluster pairs costsatmostN/(N-1).
    avg_distinct=avg*F(K.N_CLUST,K.N_CLUST-1)
    assert avg_distinct<9216 and K.NODES//16384==1<<44
    units=((1,0),(-1,0),(0,1),(0,-1),(-1,-1),(1,1))
    def add(a,b):return(a[0]+b[0],a[1]+b[1])
    def mul(a,b):return(a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]-a[1]*b[1])
    def conj(a):return(a[0]-a[1],-a[1])
    def norm(a):return a[0]*a[0]-a[0]*a[1]+a[1]*a[1]
    def val(a):return(a[0]+a[1]*K.LAM_R)%K.r
    def ratio(a,b):
        c=mul(a,conj(b));d=norm(b);g=math.gcd(math.gcd(c[0],c[1]),d)
        return(c[0]//g,c[1]//g,d//g)
    assert(K.LAM_R*K.LAM_R+K.LAM_R+1)%K.r==0
    assert all(norm(u)==1 and mul(u,conj(u))==(1,0)for u in units)
    assert all(mul(u,v)in units for u,v in itertools.product(units,repeat=2))
    def partitions(t=(0,)):
        if len(t)==4:yield t;return
        for a in range(max(t)+2):yield from partitions(t+(a,))
    parts=list(partitions());assert len(parts)==15
    norms=set();eta=set();inventory={};invs={};zero_cases=0;eligible=0
    for k in range(1,5):
        seenQ=set();seenF=set()
        for part in parts:
            if max(part)+1!=k:continue
            for perm in itertools.permutations(range(k)):
                asg=tuple(perm[z]for z in part);sizes=[asg.count(i)for i in range(k)]
                for uu in itertools.product(units,repeat=3):
                    co=[(0,0)]*k
                    for slot,u in zip(asg,uu+((1,0),)):co[slot]=add(co[slot],u)
                    if(0,0)in co:zero_cases+=1;continue
                    for c in co:
                        n=norm(c);norms.add(n)
                        assert 1<=n<=16 and math.gcd(n,K.r)==math.gcd(n,K.p)==1 and val(c)!=0
                        assert mul(c,conj(c))==(n,0)
                    den=co[-1]
                    if den not in invs:invs[den]=pow(val(den),-1,K.r)
                    qs=tuple(ratio(c,den)for c in co[:-1]);fs=tuple(val(c)*invs[den]%K.r for c in co[:-1])
                    assert tuple((a+b*K.LAM_R)*pow(d,-1,K.r)%K.r for a,b,d in qs)==fs
                    seenQ.add(qs);seenF.add(fs)
                    if k>=3:
                        for i,j in itertools.combinations(range(k),2):
                            assert min(sizes[i],sizes[j])==1
                            a,b=co[i],co[j];unit,other=(a,b)if norm(a)==1 else(b,a)
                            assert norm(unit)==1 and norm(other)in(1,3,4)
                            eta.add(norm(mul(conj(unit),other)));eligible+=1
        assert len(seenQ)==len(seenF)=={1:1,2:78,3:252,4:216}[k]
        inventory[k]=len(seenQ)
    assert norms=={1,3,4,7,9,12,13,16} and eta=={1,3,4} and zero_cases>0
    # Restrictedgrowth partitions for2..4slots into>=2selectedpair blocks.
    assignments={}
    for k in (2,3,4):
        ps=set(tuple(p[:k])for p in parts)
        assignments[k]=sum(2**(max(t)+1)for t in ps if max(t)>=1)
    assert assignments=={2:4,3:20,4:92}
    assert max(assignments.values())<=96 and 96*sum(inventory[k]for k in(2,3,4))<65536
    fiber12=sum(1 for x in range(61) for y in range(61) if (y*y-x*x*x-7)%61==0 and 22*(y*y-15)*(y*y-34)%61==1)
    assert fiber12==12 and fiber12>6
    return dict(all_checks_pass=True,pair_inverse_fiber_counterexample=dict(p=61,A0=15,A1=34,inverse=22,fiber=fiber12,wrong_six_bound_rejected=True),gf_polynomial=0x409,gf_nonzero_multiplier_rank_checks=1023,gf_residue_pairs=len(counts),gf_projection_enumerated_products=enumerated,GF_pool_size=K.NODES//16384,nonsingleton_bound_numerator=str(avg_distinct.numerator),nonsingleton_bound_denominator=str(avg_distinct.denominator),nonsingleton_bound_lt_9216=True,normalized_inventory=inventory,raw_norms=sorted(norms),eta_norms=sorted(eta),eligible_two_slot_choices=eligible,zero_coefficient_cases_excluded=zero_cases,partition_member_counts=assignments,union_assignment_allowance=65536,degree_bound=6*(1+max(eta)),controls=dict(reducible_polynomial_detected=True,rank3_projection_detected=True),scope='Finitefield/ring/inventory facts only; production divisor and probability arguments remain analytic.',seconds=time.monotonic()-started)

if __name__=='__main__':
    import json,sys
    cfg=json.loads(sys.stdin.read()or'{}')
    print(json.dumps(check(full=bool(cfg.get('full',False))),sort_keys=True))

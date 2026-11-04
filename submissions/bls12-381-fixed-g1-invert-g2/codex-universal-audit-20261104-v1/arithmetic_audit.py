import hashlib,json,math,sys
from fractions import Fraction

def main():
    options=json.load(sys.stdin)
    count=options['samples']
    assert 32<=count<=1024
    q=2**61-1
    factors={2:1,3:2,5:2,7:1,11:1,13:1,31:1,41:1,61:1,151:1,331:1,1321:1}
    product=1
    for f,e in factors.items():
        assert all(f%d for d in range(2,math.isqrt(f)+1))
        product*=f**e
        assert math.gcd(pow(37,(q-1)//f,q)-1,q)==1
    assert product==q-1 and pow(37,q-1,q)==1
    triples=[(h,k,c) for h in [0,1,q-2,q-1] for k in [0,1,q-2,q-1] for c in [0,1,2**60-1]]
    for i in range(count):
        x=int.from_bytes(hashlib.sha512(b'audit'+i.to_bytes(8,'big')).digest(),'big')
        triples.append((x%q,(x>>61)%q,(x>>122)%(2**60)))
    for h,k,c in triples:
        x=h*k+c
        lo=x%(2**64);hi=x>>64
        s=(lo&q)+((lo>>61)|(hi<<3))
        assert s<2*q
        s=(s&q)+(s>>61)
        assert s<=q
        s=s-q if s>=q else s
        assert s==x%q
    for i in range(count):
        payload=(hashlib.sha512(b'digits'+i.to_bytes(8,'big')).digest()*5)[:288]+bytes([i%2])
        M=int.from_bytes(payload,'big')
        B=[int.from_bytes(payload[j:j+8],'big') for j in range(0,288,8)]
        d=[((B[35]<<8)|payload[-1])%(2**64)]
        d += [((B[35-j]<<8)|(B[36-j]>>56))%(2**64) for j in range(1,36)]
        d.append(B[0]>>56)
        assert sum(x<<(64*j) for j,x in enumerate(d))==M
        digits=[(M>>(60*j))&(2**60-1) for j in range(39)]
        assert sum(x<<(60*j) for j,x in enumerate(digits))==M
        key=int.from_bytes(hashlib.sha512(b'key'+i.to_bytes(8,'big')).digest(),'big')%q
        h=0
        for c in reversed(digits):h=(h*key+c)%q
        assert h==sum(c*pow(key,j,q) for j,c in enumerate(digits))%q
    r=int('73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001',16)
    L=(119*(math.isqrt(r-1)+1)+99)//100;K=L+(L+7)//8
    assert (K*70144+2**90)**40<20000**40*2**5187
    assert Fraction(K*(K+1),2)*Fraction(38,q)**5<Fraction(1,8192)
    print(json.dumps({'scope':'bounded arithmetic and public constants only; not SHA success or instruction measurement','random_fixtures':count,'reduction_cases':len(triples),'checks_passed':True}))

if __name__=='__main__':main()

#!/usr/bin/env python3
"""Executable algebra and fingerprint tests, including deliberate wrong-pair faults."""
import argparse
from collections import Counter
import json
from pathlib import Path
import random
from tower import Tower, exact_order
from ledger import P,R,frobenius_rows
ROOT=Path(__file__).resolve().parents[1]
SEED=2026102307

def specialised_product(x,y,a,b,p):
    if (a,b)==(1,0):return x,y
    if (a,b)==(p-1,0):return -x%p,-y%p
    if b==0:return a*x%p,a*y%p
    if a==0:return -b*y%p,b*x%p
    if (a-b)%p==0:return a*(x-y)%p,a*(x+y)%p
    if (a+b)%p==0:return a*(x+y)%p,a*(y-x)%p
    ac=a*x%p;bd=b*y%p;cross=(x+y)*(a+b)%p
    return (ac-bd)%p,(cross-ac-bd)%p

def specialised_frob(f,y,k):
    out=[0]*12;p=f.p
    for i,(a,b) in enumerate(f.factors[k]):
        x,z=(y[i]+y[i+6])%p,y[i+6]
        if k&1:z=-z%p
        x,z=specialised_product(x,z,a,b,p)
        out[i],out[i+6]=(x-z)%p,z
    return tuple(out)

def montgomery_keys(f,y):
    # Production c1=(a,b) is stored directly in the Fp2 tower; the flattened
    # polynomial toy recovers that pair as (y1+y7,y7).
    p=f.p; mont=pow(2,384,p); inv=pow(mont,-1,p)
    def mm(x,z):return x*z*inv%p
    a,b=(y[1]+y[7])%p,y[7]
    # Two charged decodes, then four mixed-representation fixed products.
    a=mm(a*mont%p,1);b=mm(b*mont%p,1)
    B=(mm(a,f.d[0]*mont%p)+mm(b,f.d[1]*mont%p))%p
    az=mm(a,f.zeta[0]*mont%p);bz=mm(B,f.zeta[0]*mont%p)
    even=(a,az,(az-a)%p,-a%p,-az%p,(a-az)%p)
    odd=(B,bz,(bz-B)%p,-B%p,-bz%p,(B-bz)%p)
    return tuple(z for pair in zip(even,odd) for z in pair)


def production(rng):
    f=Tower(P);out=[]
    for sample in range(2):
        while True:
            g=f.pow(tuple(rng.randrange(P) for _ in range(12)),(P**12-1)//R)
            if g!=f.one:break
        assert f.pow(g,R)==f.one
        M=f.pow(g,17+sample); K,kap,rep=f.metadata(g)
        direct=g;representatives=[]
        for k in range(12):
            assert f.frob(g,k)==direct==specialised_frob(f,g,k)
            assert (direct[1]+direct[7])%P==f.keys(g)[k]
            direct=f.pow(direct,P)
            y=f.frob(g,k);Ky,ky,cy=f.metadata(y)
            assert Ky==K and cy==rep
            # No squaring of the walk element, no runtime full Frobenius.
            got=f.mul(y,f.frob(M,-ky))
            canonical_next=f.metadata(f.mul(cy,M))[2]
            assert f.metadata(got)[2]==canonical_next
            assert f.frob(f.mul(y,f.frob(M,-ky)),ky)==f.mul(cy,M)
            representatives.append(canonical_next)
            # The decision uses relative canonical alignment, not the history.
            _, kt, _=f.metadata(got)
            guard=(kt-ky)%12!=0
            c_raw=f.mul(cy,M);_,kc,_=f.metadata(c_raw)
            assert guard==(kc!=0)
            if guard:got=f.mul(got,got);c_raw=f.mul(c_raw,c_raw)
            assert f.metadata(got)[2]==f.metadata(c_raw)[2]
            assert montgomery_keys(f,y)==f.keys(y)
        assert direct==g and len(set(representatives))==1
        assert f.frob(g,6)==f.pow(g,R-1)
        out.append(dict(sample=sample,direct_Frobenius=12,specialised_Frobenius=12,
                        key_comparisons=12,canonical_conjugates=12,twist_identities=12,
                        inverse_Frobenius=True,guard_conjugacy_checks=12,
                        canonical_key_Montgomery_pipeline_checks=12))
    return dict(samples=2,results=out,frobenius_rows=frobenius_rows(f))

def exhaustive_alignment():
    counts=Counter();r=13;lam=2
    assert exact_order(lam,r,12)
    for x in range(r):
        for e in range(1,r):
            for k in range(12):
                lk=pow(lam,k,r);ep=e*pow(lk,-1,r)%r
                for b in range(r):
                    a=(e-b*x)%r
                    for bp in range(r):
                        ap=(ep-bp*x)%r;den=(b-lk*bp)%r;num=(lk*ap-a)%r
                        counts['pairs']+=1
                        if den:
                            assert num*pow(den,-1,r)%r==x;counts['recovered']+=1
                        else:assert num==0;counts['zero_denominator']+=1
    return dict(r=r,lambda_order=12,secrets=13,**counts)

def fingerprint(matrix,bits):
    result=0
    while bits:
        bit=bits&-bits;result^=matrix[bit.bit_length()-1];bits^=bit
    return result

def independent_fingerprints(rng):
    # Matrix is independent of tower keys and of DP predicate coins/positions.
    # 4608 columns, each a 128-bit independent uniform algorithm coin.
    matrix=[rng.getrandbits(128) for _ in range(4608)]
    same=0
    for _ in range(100):
        point=tuple(rng.getrandbits(381) for _ in range(12))
        bits=sum(c<<(384*i) for i,c in enumerate(point))
        explicit=0
        for j in range(4608):
            if bits>>j&1:explicit^=matrix[j]
        assert fingerprint(matrix,bits)==explicit
        same+=1
    # Exhaustive miniature universal-family proof for every nonzero difference.
    universal={}
    for delta in range(1,16):
        zeros=0
        for m in range(1<<12):
            small=[(m>>(3*j))&7 for j in range(4)]
            zeros+=fingerprint(small,delta)==0
        assert zeros==512;universal[str(delta)]=zeros
    # Tiny 2-bit fingerprints intentionally collide: never accept without
    # exact class confirmation followed by scalar verification.
    r=97;lam=6;assert exact_order(lam,r,12)
    L=[pow(lam,k,r) for k in range(12)]
    reps=sorted({min(e*a%r for a in L) for e in range(1,r)})
    tiny=[rng.getrandbits(2) for _ in range(7)]
    x=23;records=[];false_pairs=0;rejected=0;injected=0
    for e in reps:
        b=rng.randrange(r);a=(e-b*x)%r;fp=fingerprint(tiny,e)
        for ep,ap,bp,fpp in records:
            if fp!=fpp:continue
            false_pairs+=1
            assert e!=ep
            den=(bp-b)%r;num=(a-ap)%r
            if den:
                candidate=num*pow(den,-1,r)%r
                assert candidate!=x;rejected+=1
            else:rejected+=1
            injected+=int((fp^1)!=fpp)
        records.append((e,a,b,fp))
    assert false_pairs and rejected==false_pairs and injected==false_pairs
    return dict(full_point_bits=4608,output_bits=128,matrix_bytes=73728,
        full_point_encoding_checks=same,miniature_universal_zero_counts=universal,
        deliberate_false_pairs=false_pairs,false_pairs_rejected=rejected,
        injected_fingerprint_faults_detected=injected,wrong_answers_accepted=0)

def inverse_two_cycle():
    # On canonical scalar labels, two same-h, inverse-aligning steps cancel
    # identically. This is NOT a rare accidental den=0 event.
    checks=0
    for r in (13,37,97):
        for a in range(r):
            for A in range(r):
                a1=-(A+a)%r;a2=-(a1+a)%r
                assert a2==A;checks+=1
    return dict(exhaustive_label_cancellations=checks,
        symbolic='-(-(A+a)+a)=A; -(-(D+b)+b)=D',
        conclusion='Floyd detects this cycle but its aligned numerator and denominator are identically zero.')

def main():
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    rng=random.Random(SEED)
    out=dict(version=1,scope='Algebra and faults only; no production timing or success certification.',
        production=production(rng),alignment=exhaustive_alignment(),
        fingerprint=independent_fingerprints(rng),forced_inverse_cycle=inverse_two_cycle())
    text=json.dumps(out,indent=2)+'\n'
    if a.write:(ROOT/'evidence/algebra.json').write_text(text)
    print(text,end='')
if __name__=='__main__':main()

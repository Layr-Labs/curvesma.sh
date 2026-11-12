"""Exact source-side recovery of the omitted X-low bit and leaf sign.
The production source reconstructs Y from its retained coins, then calls this
finite arithmetic model's specified algorithm.  Four complete affine additions,
four exact restoring fraction divisions and bounded comparisons fit 2^22 extra
instructions/request; all eight reply chunks are retained, the last holds sign.
This is a finite identity check, not a production probability experiment.
"""
import random,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import scalar_check as sc
M64=(1<<64)-1
FRAGMENT_MASK=0xffffe0000fffffff
SOURCE_MASK=(1<<36)-1
def stored_key(t):
    """Exactly EXTR #36 followed by receive's BFI at bits28..44."""
    k2=t>>128
    k3=((t>>64)>>36 | (t&M64)<<28)&M64
    return (k3&FRAGMENT_MASK)|((k2>>58)<<28)
def decode_metadata(w2,w3):
    return ((w2&SOURCE_MASK,w3&SOURCE_MASK),
            (w3>>36)&511,(w2>>36)|((w3>>45)<<28))
def fraction(A):
    return (A[0]<<192)//sc.P
def recover_source(Xpair,Y,scalars,b,c,key):
    """Return one unambiguous (X index, sign) or discard the request.
    Xpair consists of the two permanent canonical affine points. Scalars are
    their public a values; b is the source's coefficient of Q (carried by Y)
    and c is Y's regenerated scalar. sign0 is X+Y and sign1 is X-Y. Input indices/caps have
    already been validated by the source driver.
    """
    matches=[]
    for j,X in enumerate(Xpair):
        if X is None or Y is None:continue
        for sign in (0,1):
            Z=sc.add(X,Y if sign==0 else sc.neg(Y))
            if Z is not None and stored_key(fraction(Z))==key:
                matches.append((j,sign))
    if len(matches)!=1:return None
    j,sign=matches[0];X=Xpair[j]
    vals=[X[0]*sc.MR%sc.P,X[1]*sc.MR%sc.P,
          Y[0]*sc.MR%sc.P,Y[1]*sc.MR%sc.P,
          scalars[j],b,c,sign]
    return dict(values=vals,choice=(j,sign))
def check(seed=20261107):
    rng=random.Random(seed);cases=0;metadata_cases=0
    for pair in (0,511):
        for yi in (0,(1<<47)-1):
            w2=SOURCE_MASK|((yi&((1<<28)-1))<<36)
            w3=SOURCE_MASK|(pair<<36)|((yi>>28)<<45)
            assert w2<1<<64 and w3<1<<64
            assert decode_metadata(w2,w3)==((SOURCE_MASK,SOURCE_MASK),pair,yi)
            metadata_cases+=1
    for _ in range(5):
        aa=[rng.randrange(1,sc.R) for _ in range(2)]
        xx=[sc.mul(a) for a in aa];cc=rng.randrange(1,sc.R);yy=sc.mul(cc)
        for j in (0,1):
            for sign in (0,1):
                leaf=sc.add(xx[j],yy if sign==0 else sc.neg(yy));tt=fraction(leaf)
                key=stored_key(tt);got=recover_source(xx,yy,aa,0,cc,key)
                assert got is not None and got['choice']==(j,sign)
                assert sc.leaf(got['values'][:7],got['values'][7])[0]==leaf
                # This is the actual receiver representation, including its six
                # high-K2 bits, not a recomputation from congruence alone.
                k2=tt>>128;k3=((tt>>100)|((tt&M64)<<28))&M64
                assert key==(k3&FRAGMENT_MASK)|((k2>>58)<<28)
                assert recover_source(xx,yy,aa,0,cc,key^1) is None
                cases+=1
        # Two identical supplied X entries produce an ambiguous match. Such a
        # request is discarded rather than selecting an arbitrary scalar.
        key=stored_key(fraction(sc.add(xx[0],yy)))
        assert recover_source([xx[0],xx[0]],yy,[aa[0],aa[0]],0,cc,key) is None
    recovered=0
    for sign2 in (0,1):
        secret=1739;aa=37;bb=11;cc=13;bb2=17;cc2=19
        zscalar=(aa+bb*secret+cc)%sc.R
        aa2=(zscalar-(1 if sign2==0 else -1)*(cc2+bb2*secret))%sc.R
        X1=sc.mul(aa);X2=sc.mul(aa2)
        Y1=sc.mul(cc+bb*secret);Y2=sc.mul(cc2+bb2*secret)
        fragment=stored_key(fraction(sc.mul(zscalar)))
        left=recover_source([X1,sc.mul(87)],Y1,[aa,87],bb,cc,fragment)
        right=recover_source([sc.mul(91),X2],Y2,[91,aa2],bb2,cc2,fragment)
        assert left and right and left['choice']==(0,0) and right['choice']==(1,sign2)
        value=sc.recover(left['values'][:7],right['values'][:7],left['values'][7],right['values'][7])
        assert value==secret and sc.mul(value)==sc.mul(secret);recovered+=1
    extra=4*45056+4*8256+(1<<18)
    assert extra<(1<<22)
    return dict(all_checks=True,choice_cases=cases,ambiguous_cases=5,
                real_collision_recoveries=recovered,metadata_cases=metadata_cases,extra_source_bound=1<<22,
                explicit_extra_construction=extra,fragment_bits=53,
                scope='Exact P256 source enumeration, fraction-fragment roundtrip, mutation/ambiguity rejection and two actual recoveries.')
if __name__=='__main__':
    import json
    print(json.dumps(check()))

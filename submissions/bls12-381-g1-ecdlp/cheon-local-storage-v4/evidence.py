"""Execute only using the pinned llm-verifier credential-free Docker runner."""
import json
import random
import runpy
import sys
from pathlib import Path

ROOT = Path('/snapshot/files')
core = runpy.run_path(str(ROOT/'herd_core.py'),run_name='core')
word = runpy.run_path(str(ROOT/'word_kernels.py'),run_name='word')
ledger = runpy.run_path(str(ROOT/'ledger.py'),run_name='ledger')


def addresses(seed):
    rng = random.Random(seed)
    maximum = 0
    seen = set()
    for kind in range(2):
        for row in range(8):
            digits = [0,1,(1 << 32)-1]+[rng.getrandbits(32) for _ in range(100)]
            for digit in digits:
                address = 4096+(((kind*8+row)*(1 << 32))+digit)*96
                assert 4096 <= address and address+96 <= (1 << 43)
                assert address < 1 << 64
                identity = (kind,row,digit)
                if identity not in seen:
                    seen.add(identity)
                maximum = max(maximum,address+95)
    for index in [0,1,(1 << 36)-1]+[rng.getrandbits(36) for _ in range(1000)]:
        address = (1 << 43)+index*128
        assert (1 << 43) <= address and address+128 <= (1 << 44)
        maximum = max(maximum,address+127)
    return {'local_address_bits_used':maximum.bit_length(),'maximum_checked_address':str(maximum),
            'point_entry_vectors':16*103,'record_entry_vectors':1003,
            'scope':'Production layout arithmetic, not a production allocation or hardware timing.'}


def production_algebra(seed,count=16):
    profile = json.loads((ROOT/'target-profile.json').read_text())
    p,r = int(profile['field_modulus'],16),int(profile['subgroup_order'],16)
    E = core['helpers']['Curve'](p,4)
    P = tuple(int(profile['generator_P0'][k],16) for k in ('x','y'))
    x = -int(profile['parameter_x'][1:],16)
    lam = (x*x-1) % r
    image = E.mul(lam,P)
    beta = image[0]*pow(P[0],-1,p) % p
    assert image[1]==P[1] and pow(beta,3,p)==1 and beta!=1
    units=[1,r-1,lam,r-lam,lam*lam % r,(-lam*lam) % r]
    assert len(set(units))==6 and all(pow(u,(r-1)//(1 << 20),r)==1 for u in units)
    rng=random.Random(seed)
    key=lambda A:min(A[1]*pow(2,384,p) % p,(-A[1]*pow(2,384,p)) % p)
    for _ in range(count):
        coefficient=rng.randrange(1,r)
        A=E.mul(coefficient,P)
        orbit=[E.mul(u,A) for u in units]
        assert len(set(orbit))==6 and len({key(S) for S in orbit})==1
        partial,total=None,0
        # Production 32-bit decomposition; compute individual table entries on
        # demand for checking, never claim this is production table construction.
        for j in range(8):
            digit=(coefficient >> (32*j)) & ((1 << 32)-1)
            if digit:
                chunk=digit << (32*j)
                Q=E.mul(chunk,P)
                if partial is not None:
                    assert total<chunk and total+chunk<=coefficient<r
                    assert partial[0]!=Q[0]
                partial=E.add(partial,Q)
                total+=chunk
        assert partial==A and total==coefficient
        h=pow(rng.randrange(1,r),1 << 20,r)
        c_mont=coefficient*(1 << 256) % r
        h_mont=h*(1 << 256) % r
        assert word['montgomery4'](c_mont,1,r)[0]==coefficient
        next_mont=word['montgomery4'](c_mont,h_mont,r)[0]
        assert word['montgomery4'](next_mont,1,r)[0]==coefficient*h % r
        assert key(E.mul(h,A))==key(E.mul(h,orbit[1]))
    return {'production_points':count,'Montgomery_key_orbit_checks':6*count,
            'scalar_Montgomery_decode_and_transition_checks':count,
            'fixed_base_32bit_decompositions':count,'nondegeneracy_checks_passed':True,
            'scope':'Production curve identities and sampled decompositions; not a production search or probability proof.'}


def cap_checks(seed):
    rng=random.Random(seed)
    r=core['helpers']['toy_prime'](22,16,rng)
    lam=next(pow(a,(r-1)//3,r) for a in range(2,r) if pow(a,(r-1)//3,r)!=1)
    units=[1,r-1,lam,r-lam,lam*lam % r,(-lam*lam) % r]
    canonical=lambda c:min(u*c % r for u in units)
    label=lambda kind,c:(canonical(c),0)
    recover=lambda a,b:None
    cases=[('table_cap',{'table_cap':0,'slots':8,'dp_bits':0,'endpoint_cap':100}),
           ('endpoint_cap',{'table_cap':100,'slots':512,'dp_bits':0,'endpoint_cap':0}),
           ('probe_cap',{'table_cap':100,'slots':512,'dp_bits':0,'probe_cap':0}),
           ('chain_cap',{'table_cap':100,'slots':512,'dp_bits':88,'chain_cap':0})]
    rows=[]
    for reason,options in cases:
        row=core['one_trial'](r,16,units,random.Random(seed),label,recover,**options)
        assert not row['success'] and row['reason']==reason
        rows.append({'expected_failure':reason,'observed_failure':row['reason']})
    return rows


def main():
    cfg=json.load(sys.stdin)
    seed=cfg.get('seed',1); phase=cfg.get('phase','full')
    out={'phase':phase,'seed':seed,'ledger':ledger['compute']()}
    profile=json.loads((ROOT/'target-profile.json').read_text())
    if phase in ('full','algebra'):
        out['addresses']=addresses(seed)
        out['production_algebra']=production_algebra(seed,cfg.get('production_points',12))
        out['scalar_word_kernel']=word['check'](int(profile['subgroup_order'],16),seed,256,4)
        out['base_field_word_kernel']=word['check'](int(profile['field_modulus'],16),seed+1,256,6)
        out['forced_cap_checks']=cap_checks(seed)
    if phase in ('full','curve'):
        out['curve_trials']=core['curve_trials'](seed,cfg.get('trials',100),cfg.get('bits',24),
                                                  cfg.get('logd',4),cfg.get('dp_bits',2),True)
    if phase in ('full','scalar'):
        out['scalar_trials']=core['scalar_trials'](seed+100,cfg.get('scalar_trials',500),
                                                  cfg.get('bits',28),cfg.get('logd',6),cfg.get('dp_bits',2))
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()

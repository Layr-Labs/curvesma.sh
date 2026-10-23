"""Captured real-curve toy pipeline plus independent arithmetic controls.
This file does not simulate production cost or infer production success rates.
"""
import hashlib
import itertools
import json
import math
import random
from mechanism import (Field, affine_add, affine_scalar, jac, key, make_tree,
                       mixed, normalize, progression, recover, twice)

class ToyField:
    """Small prime arithmetic; SAME Jacobian formulas, no scalar/exponent oracle."""
    def __init__(self, p):
        self.p = p
    def enc(self, a):
        return a % self.p
    def dec(self, a):
        return a % self.p
    def mul(self, a, b):
        return a * b % self.p
    def add(self, a, b):
        return (a + b) % self.p
    def sub(self, a, b):
        return (a - b) % self.p
    def inv(self, a):
        return pow(a, -1, self.p)

def captured(check, profile, m, length, h, r, p):
    gen = tuple(int(profile['generator_P'][c], 16) for c in ('x', 'y'))
    bcurve = int(profile['coefficient_b'], 16)
    check('production_generator_curve', (gen[1]**2-gen[0]**3+3*gen[0]-bcurve) % p == 0)
    check('production_generator_order', affine_scalar(r, gen, p) is None)
    lattice_counts = []
    for digits in range(1, 6):
        vals = [sum(v * 4**i for i,v in enumerate(ds))
                for ds in itertools.product((-3,-1,1,3), repeat=digits)]
        want = list(range(1-4**digits, 4**digits, 2))
        check(f'radix4_bijection_{digits}', sorted(vals) == want)
        positive = sorted(v for v in vals if v > 0)
        check(f'positive_half_exact_{digits}', positive == list(range(1,4**digits,2)))
        lattice_counts.append({'digits':digits,'strings':len(vals),'distinct':len(set(vals)),
                               'positive':len(positive)})
    f = Field(p)
    operands = [0, 1, 2, 3, p-1, p-2, (1<<64)-1, 1<<64,
                (1<<64)+1, (1<<128)-1, 1<<128,
                (1<<192)-1, 1<<192, (1<<255)-1, 1<<255, *gen, bcurve]
    rng = random.Random(20261023)
    operands.extend(rng.randrange(p) for _ in range(8))
    pairs, good = 0, True
    for a in operands:
        for b in operands:
            before = f.pairs
            good &= f.cios(a,b) == a*b*f.radix_inv % p and f.pairs-before == 32
            pairs += 1
    check('CIOS_fullwidth_pairs', good and pairs == 676, pairs)
    branch_cases = 0
    for scalar in [0,1,2,3,17]:
        point = affine_scalar(scalar, gen, p)
        for scale in [1,2,p-1]:
            jp = jac(None,f) if point is None else (
                f.enc(point[0]*scale**2 % p),f.enc(point[1]*scale**3 % p),f.enc(scale))
            check(f'P256_double_{scalar}_{scale}', normalize([twice(jp,f)],f)[0] == affine_add(point,point,p))
            for label,addend in [('O',None),('P',gen),('negP',(gen[0],-gen[1]%p)),('self',point)]:
                before = f.calls
                got = mixed(jp,addend,f)
                calls = f.calls-before
                check(f'P256_mixed_{scalar}_{scale}_{label}',
                      calls <= 16 and normalize([got],f)[0] == affine_add(point,addend,p))
                branch_cases += 1
    check('Y_zero_branch', twice((f.enc(1),0,f.enc(1)),f)[2] == 0)
    # Production-field chain starts at P, not O, and stores [2j+1]P.
    baby = progression(gen, affine_scalar(2,gen,p),17,f,batch=4)
    check('P256_odd_chain_all_labels', baby == [affine_scalar(2*j+1,gen,p) for j in range(17)])
    check('negative_control_j_label_is_wrong', baby[0] != affine_scalar(0,gen,p))
    start = affine_scalar(77,gen,p)
    stride = affine_scalar(64,gen,p)
    assert stride is not None
    giant = progression(start,(stride[0],-stride[1]%p),17,f,batch=4)
    check('P256_giant_chain_all_labels', giant == [affine_scalar(77-64*i,gen,p) for i in range(17)])
    records = [(gen,7),((gen[0],-gen[1]%p),8),(None,0),(affine_scalar(2,gen,p),9),(gen,10)]
    tree, drained = make_tree(records,8), []
    while True:
        item = tree.pull()
        if item is None:
            break
        drained.append(item)
    check('sorted_tree_duplicates_and_identity', drained == sorted(records,key=key))
    check('END_sticky', tree.pull() is None and tree.pull() is None)
    tp, tb, tr, tg, tm = 97, 6, 107, (0,43), 8
    points = [(x,y) for x in range(tp) for y in range(tp)
              if (y*y-x**3+3*x-tb) % tp == 0]
    check('toy_counted_prime_order', len(points)+1 == tr and all(tr%d for d in range(2,math.isqrt(tr)+1)))
    check('toy_generator_on_curve_order', tg in points and affine_scalar(tr,tg,tp) is None)
    th = (tr+2*tm-1)//(2*tm)
    tf = ToyField(tp)
    # Cross-check optimized tiny-field arithmetic against inherited full CIOS.
    got_fast, nv_fast = recover(affine_scalar(80,tg,tp),tg,tp,tr,tm,th,tf)
    got_cios, nv_cios = recover(affine_scalar(80,tg,tp),tg,tp,tr,tm,th,Field(tp))
    check('toy_CIOS_and_smallfield_same_full_pipeline', (got_fast,nv_fast) == (got_cios,nv_cios) == (80,1))
    observed = []
    max_verify = 0
    for secret in range(tr):
        got, nv = recover(affine_scalar(secret,tg,tp),tg,tp,tr,tm,th,tf)
        observed.append(got)
        max_verify = max(max_verify,nv)
    check('toy_exhaustive_recovered', observed == list(range(tr)), f'{sum(a==b for a,b in enumerate(observed))}/{tr}')
    check('toy_zero_recovered_without_shortcut', observed[0] == 0)
    seed, trials = 2026102308, 1000
    rng = random.Random(seed)
    transcript, successes = [], 0
    for _ in range(trials):
        secret = rng.randrange(tr)
        got, nv = recover(affine_scalar(secret,tg,tp),tg,tp,tr,tm,th,tf)
        successes += got == secret
        max_verify = max(max_verify,nv)
        transcript.append([secret,got,nv])
    check('toy_random_1000_zero_misses', successes == trials, f'{successes}/{trials}')
    check('toy_max_two_verifications', max_verify <= 2, max_verify)
    # Finite capture is split into its own file by the certificate generator.
    edge_secrets = [0,1,2,3,7,8,15,16,17,31,32,33,105,106]
    edges = [{'x':x,'recovered':observed[x]} for x in edge_secrets]
    # End cell can overhang 2r: 208+15=223=9+2*107. Reduction must precede filtering.
    raw, tc = 223, 9
    check('overhang_match_is_real', affine_scalar(raw,tg,tp) == affine_scalar(tc,tg,tp))
    check('overhang_two_subtractions_required', raw >= 2*tr and raw-2*tr == tc)
    check('raw_filter_would_discard_correct_match', not (1 <= raw < 2*tr) and tc & 1)
    fixtures = []
    secrets = [0,1,2,3,length-2,length-1,length,length+1,2*length-1,r//2,r-2,r-1]
    for x in secrets:
        z = x if x & 1 else x+r
        i = (z-1)//(2*length)
        center = (2*i+1)*length
        delta = z-center
        odd = abs(delta)
        j = (odd-1)//2
        check('production_witness_range_'+str(x), 0 <= i < h and 0 <= j < m and delta & 1 and 1 <= odd < length)
        q = affine_scalar(x,gen,p)
        bp = affine_scalar(odd,gen,p)
        cp = affine_scalar(center,gen,p)
        assert cp is not None and bp is not None
        gp = affine_add(q,(cp[0],-cp[1]%p),p)
        want = bp if delta > 0 else (bp[0],-bp[1]%p)
        check('production_witness_point_'+str(x), gp == want)
        signed = center+odd if delta > 0 else center-odd
        c = signed-r if signed >= r else signed
        c = c-r if c >= r else c
        check('production_witness_recovery_'+str(x), c == x and affine_scalar(c,gen,p) == q)
        fixtures.append({'x':str(x),'odd_lift':str(z),'i':str(i),'j':str(j),'sign':1 if delta>0 else -1,'recovered':str(c)})
    summary = {'curve':{'p':tp,'a':tp-3,'b':tb,'r':tr,'P':list(tg)},'m':tm,'L':2*tm,'h':th,
               'exhaustive_trials':tr,'exhaustive_successes':sum(x==v for x,v in enumerate(observed)),
               'random_seed':seed,'random_trials':trials,'random_successes':successes,'max_verifications':max_verify,
               'edge_cases':edges,'random_transcript_sha256':hashlib.sha256(json.dumps(transcript,separators=(',',':')).encode()).hexdigest(),
               'CIOS_operand_pairs':pairs,'mixed_cases':branch_cases,'production_witness_fixtures':len(fixtures),
               'scope':'Toy executes full ECC pipeline; P256 fixtures check witness relations only, not a production search.'}
    print('CAPTURED_TOY '+json.dumps(summary,sort_keys=True))
    return {'summary':summary,'radix4_enumerations':lattice_counts,'P256_witnesses':fixtures}, transcript

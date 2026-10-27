"""Two labelled herds on the implicit multiplicative Cheon subgroup.

The production specification uses addressable fixed-base tables and batched
affine additions. These scaled tests execute the same scalar transitions,
full-key collision test and recovery; they do not allocate production tables.
"""
import hashlib
import math
import random
import runpy
from pathlib import Path

helpers = runpy.run_path('/snapshot/files/curve_support.py', run_name='helpers')


def units_for_curve(E, r, P):
    beta = next(pow(a,(E.p-1)//3,E.p) for a in range(2,100)
                if pow(a,(E.p-1)//3,E.p) != 1)
    roots = {pow(a,(r-1)//3,r) for a in range(2,100)}-{1}
    lam = next(x for x in roots if E.mul(x,P)==(beta*P[0] % E.p,P[1]))
    return [1,r-1,lam,r-lam,lam*lam % r,(-lam*lam) % r]


def sample_nonzero(rng, r, cap=128):
    for _ in range(cap):
        a = rng.getrandbits(r.bit_length())
        if 0 < a < r:
            return a
    raise RuntimeError('random_draw_cap')


def jumps(rng, r, d, units, bits):
    values, used = [], set()
    canonical = lambda x: min(u*x % r for u in units)
    for _ in range(1 << bits):
        for attempt in range(128):
            h = pow(sample_nonzero(rng,r),d,r)
            a,b = canonical(h),canonical(pow(h,-1,r))
            if a != 1 and a != b and a not in used and b not in used:
                values.append(h)
                used.add(a)
                break
        else:
            raise RuntimeError('jump_draw_cap')
    return values


def point_digest(key, salt):
    raw = key.to_bytes(48,'big')
    salted = bytes(a^b for a,b in zip(raw[:32],salt))+raw[32:]
    return hashlib.sha256(salted).digest()


def table_build(E, base, r, window):
    tables = []
    shift = base
    for _ in range((r.bit_length()+window-1)//window):
        row = [None, shift]
        for k in range(2,1 << window):
            row.append(E.add(row[-1],shift))
        tables.append(row)
        shift = E.mul(1 << window,shift)
    return tables


def fixed_multiply(E, coefficient, tables, r, window):
    assert 0 < coefficient < r
    point, partial, additions = None, 0, 0
    for j,row in enumerate(tables):
        digit = (coefficient >> (window*j)) & ((1 << window)-1)
        if not digit:
            continue
        chunk = digit << (window*j)
        if point is not None:
            assert partial < chunk and partial+chunk <= coefficient < r
            assert point[0] != row[digit][0]
            additions += 1
        point = E.add(point,row[digit])
        partial += chunk
    assert partial == coefficient and point is not None
    return point, additions


def one_trial(r, d, units, rng, label, recover, lanes=8, dp_bits=2, jump_bits=6,
              table_cap=None, slots=None, probe_cap=128, chain_cap=None, endpoint_cap=None):
    assert lanes % 2 == 0 and len(units) == 6
    m, order = (r-1)//d, (r-1)//(6*d)
    assert all(pow(u,m,r)==1 for u in units)
    factor = jumps(rng,r,d,units,jump_bits)
    salt = rng.getrandbits(256).to_bytes(32,'big')
    maximum = lanes*((4*(math.isqrt(order)+1)+lanes-1)//lanes)
    if table_cap is None:
        table_cap = 8*((maximum >> dp_bits)+1)
    if slots is None:
        slots = 1 << (4*table_cap).bit_length()
    if endpoint_cap is None:
        endpoint_cap = table_cap
    assert slots & (slots-1) == 0
    if chain_cap is None:
        chain_cap = 256*(1 << dp_bits)
    records = [None]*slots
    coefficients = [pow(sample_nonzero(rng,r),d,r) for _ in range(lanes)]
    age = [0]*lanes
    occupied = states = endpoints = restarts = duplicate_endpoints = 0
    additions = 0
    reason = 'budget'
    for _ in range(maximum//lanes):
        labels = []
        for i,c in enumerate(coefficients):
            key,count = label(i & 1,c)
            labels.append((key,point_digest(key,salt)))
            additions += count
        states += lanes
        for i,(key,digest) in enumerate(labels):
            is_dp = int.from_bytes(digest[:12],'big') >> (96-dp_bits) == 0
            if is_dp:
                endpoints += 1
                if endpoints > endpoint_cap:
                    return {'success':False,'reason':'endpoint_cap','states':states,'endpoints':endpoints}
                index = int.from_bytes(digest[-8:],'little') & (slots-1)
                for probe in range(probe_cap):
                    old = records[index]
                    if old is None:
                        if occupied >= table_cap:
                            reason = 'table_cap'
                            return {'success':False,'reason':reason,'states':states,'endpoints':endpoints}
                        records[index] = (key,coefficients[i],i & 1)
                        occupied += 1
                        break
                    if old[0] == key:
                        duplicate_endpoints += 1
                        if old[2] != (i & 1):
                            tame,wild = ((old[1],coefficients[i]) if old[2]==0
                                         else (coefficients[i],old[1]))
                            solution = recover(tame,wild)
                            return {'success':True,'reason':'cross_herd','solution':solution,
                                    'states':states,'endpoints':endpoints,'records':occupied,
                                    'duplicate_endpoints':duplicate_endpoints,
                                    'fixed_base_additions':additions,'maximum_states':maximum}
                        break
                    index = (index+1) & (slots-1)
                else:
                    return {'success':False,'reason':'probe_cap','states':states,'endpoints':endpoints}
                coefficients[i] = pow(sample_nonzero(rng,r),d,r)
                age[i] = 0
                restarts += 1
            else:
                coefficients[i] = coefficients[i]*factor[int.from_bytes(digest[12:14],'little')
                                                         & ((1 << jump_bits)-1)] % r
                age[i] += 1
                if age[i] > chain_cap:
                    return {'success':False,'reason':'chain_cap','states':states,'endpoints':endpoints}
    return {'success':False,'reason':reason,'states':states,'endpoints':endpoints,
            'duplicate_endpoints':duplicate_endpoints,
            'records':occupied,'fixed_base_additions':additions,'maximum_states':maximum}


def summarize(rows):
    successes = sum(row['success'] for row in rows)
    failures = {}
    for row in rows:
        if not row['success']:
            failures[row['reason']] = failures.get(row['reason'],0)+1
    return {'trials':len(rows),'successes':successes,'success_rate':successes/len(rows),
            'failure_reasons':failures,'average_states':sum(x['states'] for x in rows)/len(rows),
            'maximum_duplicate_endpoints':max(x.get('duplicate_endpoints',0) for x in rows),
            'mean_duplicate_endpoints':sum(x.get('duplicate_endpoints',0) for x in rows)/len(rows),
            'all_returned_scalars_correct':True}


def curve_trials(seed, trials=60, bits=24, logd=4, dp_bits=2, montgomery_keys=True):
    rng = random.Random(seed)
    d = 1 << logd
    E,r,P = helpers['find_curve'](bits,d,rng)
    units = units_for_curve(E,r,P)
    window = 6
    tame_tables = table_build(E,P,r,window)
    rows = []
    for _ in range(trials):
        secret = rng.randrange(1,r)
        target,t = E.mul(secret,P),pow(secret,d,r)
        Ad = E.mul(t,P)
        wild_tables = table_build(E,Ad,r,window)

        def label(kind,c):
            out,count = fixed_multiply(E,c,wild_tables if kind else tame_tables,r,window)
            stored_y = out[1]*pow(2,384,E.p) % E.p if montgomery_keys else out[1]
            return min(stored_y,E.p-stored_y),count

        def recover(tame,wild):
            ratio = tame*pow(wild,-1,r) % r
            found = next(u*ratio % r for u in units if E.mul(u*ratio % r,P)==Ad)
            root = helpers['dth_root'](found,d,r)
            # The same d roots as production, exhaustively enumerated at toy size.
            generator = next(a for a in range(2,r) if pow(pow(a,(r-1)//d,r),d//2,r)!=1)
            omega = pow(generator,(r-1)//d,r)
            candidate = root
            for j in range(d):
                if E.mul(candidate,P)==target:
                    assert candidate == secret
                    return candidate
                candidate = candidate*omega % r
            raise AssertionError('root recovery failed')

        row = one_trial(r,d,units,rng,label,recover,dp_bits=dp_bits)
        if row['success']:
            assert row['solution'] == secret
        rows.append(row)
    result = summarize(rows)
    result.update({'p':E.p,'r':r,'d':d,'quotient_order':(r-1)//(6*d),
                   'dp_bits':dp_bits,'lanes':8,'jump_bits':6,'fixed_base_window':window,
                   'key_representation':'Montgomery y*2^384 mod p' if montgomery_keys else 'ordinary y',
                   'scope':'Scaled end-to-end curve attack; includes failed capped trials and exact scalar verification, not a production probability proof'})
    return result


def scalar_trials(seed,trials=200,bits=28,logd=6,dp_bits=2):
    rng = random.Random(seed)
    d = 1 << logd
    r = helpers['toy_prime'](bits,d,rng)
    lam = next(pow(a,(r-1)//3,r) for a in range(2,r) if pow(a,(r-1)//3,r)!=1)
    units = [1,r-1,lam,r-lam,lam*lam % r,(-lam*lam) % r]
    rows=[]
    for _ in range(trials):
        secret=rng.randrange(1,r)
        t=pow(secret,d,r)
        def label(kind,c):
            value = c*t % r if kind else c
            return min(u*value % r for u in units),0
        def recover(tame,wild):
            ratio=tame*pow(wild,-1,r) % r
            assert t in {u*ratio % r for u in units}
            return secret
        rows.append(one_trial(r,d,units,rng,label,recover,dp_bits=dp_bits))
    result=summarize(rows)
    result.update({'r':r,'d':d,'quotient_order':(r-1)//(6*d),'dp_bits':dp_bits,
                   'scope':'Finite scalar quotient with SHA-labelled orbits; success-model evidence only, no curve-operation cost evidence'})
    return result

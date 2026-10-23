"""Own, bounded toy-only GLV DP pipeline. No production parameters accepted.
Run: python3 -B toy_pipeline.py > toy_pipeline_results.json
The legacy mean-repeat-time rho is NOT the exponent rho measured here.
"""
import json
import math
import random
from collections import defaultdict

P, B, Q = 524347, 13, 523969
N, DP_BITS, TABLE = (Q - 1) // 6, 4, 256
TRIALS, TABLES, KEY_BITS = 2000, 8, 12

def add(a, b):
    if a is None: return b
    if b is None: return a
    x, y = a; X, Y = b
    if x == X and (y + Y) % P == 0: return None
    s = ((3*x*x)*pow(2*y, -1, P) if x == X else (Y-y)*pow((X-x) % P, -1, P)) % P
    z = (s*s-x-X) % P
    return z, (s*(x-z)-y) % P

def cube_root(mod):
    g = 2
    while pow(g, (mod-1)//3, mod) == 1: g += 1
    return pow(g, (mod-1)//3, mod)

assert all(n % d for n in (P, Q) for d in range(2, math.isqrt(n)+1))
beta, lam = cube_root(P), cube_root(Q)
beta2 = beta*beta % P
ys = {}
for x in range(P):
    v = (x*x % P*x+B) % P
    if pow(v, (P-1)//2, P) == 1:
        y = pow(v, (P+1)//4, P)
        ys[x] = y
assert 2*len(ys)+1 == Q and 0 not in ys
G = next(iter(ys.items()))
point_log, scalar_point = {}, [None]*Q
cur = None
for k in range(1, Q):
    cur = add(cur, G); point_log[cur] = k; scalar_point[k] = cur
assert add(cur, G) is None
if point_log[(beta*G[0] % P, G[1])] != lam:
    lam = lam*lam % Q
assert point_log[(beta*G[0] % P, G[1])] == lam
units = (1, lam, lam*lam % Q)

def canon(a):
    if a is None: return None, 1
    x, y = a
    xs = (x, beta*x % P, beta2*x % P)
    j = min(range(3), key=xs.__getitem__)
    u = units[j]
    if y & 1: y, u = P-y, -u % Q
    return (xs[j], y), u

reps = {}
for x, y in ys.items():
    a, _ = canon((x, y)); reps[a[0]] = a
xs = sorted(reps)
assert len(xs) == N
T = math.isqrt((98859*108*108*N)//(100000*100*100))
while T*T*100000*100*100 < 98859*108*108*N: T += 1
horizon = 40*(1 << DP_BITS)
rows = []
for table_id in range(TABLES):
    rr = random.Random(8007+table_id)
    tlogs = [rr.randrange(1, Q) for _ in range(TABLE)]
    points = [scalar_point[k] for k in tlogs]
    trans, doubles = {}, {}
    for x in xs:
        j = (x >> DP_BITS) & (TABLE-1)
        a, u = canon(add(reps[x], points[j]))
        trans[x] = (None if a is None else a[0], u, j)
        a, u = canon(add(reps[x], reps[x]))
        doubles[x] = (a[0], u)
    successes = core_collisions = recoveries = spurious = escapes = aborts = work = starts = 0
    for trial in range(TRIALS):
        rng = random.Random(1000+table_id*100003+trial)
        secret = rng.randrange(1, Q)
        bs = [rng.randrange(Q) for _ in range(TABLE)]
        aa = [(k-b*secret) % Q for k, b in zip(tlogs, bs)]
        table = defaultdict(list)
        first = {}
        used = 0; trail = 0; any_core = False; good = False; recovered = False
        while used < T:
            trail += 1; starts += 1
            x = xs[rng.randrange(N)]
            b = rng.randrange(Q); a = (point_log[reps[x]]-b*secret) % Q
            ring = []; checkpoint = (x, ()); power = 1; age = 0; length = 0
            core = set()
            while length < horizon:
                if used < T:
                    core.add(x)
                    if x in first and first[x] != trail: any_core = True
                    first.setdefault(x, trail)
                xn, u, j = trans[x]
                a, b = u*(a+aa[j]) % Q, u*(b+bs[j]) % Q
                used += 1; length += 1
                if xn is None:
                    # An actual identity yields the logarithm directly if b != 0.
                    if b:
                        assert (-a*pow(b, -1, Q)) % Q == secret
                        good = recovered = True
                    break
                x = xn
                if (x & ((1 << DP_BITS)-1)) == 0:
                    key = x & ((1 << KEY_BITS)-1)
                    for oldx, olda, oldb, oldcore in table[key]:
                        if oldx != x:
                            spurious += 1; continue
                        if oldb != b:
                            candidate = (olda-a)*pow((b-oldb) % Q, -1, Q) % Q
                            assert candidate == secret and scalar_point[candidate] == scalar_point[secret]
                            recovered = True
                            if core.intersection(oldcore): good = True
                    table[key].append((x, a, b, core))
                    break
                if x in ring:
                    cycle = ring[ring.index(x):]
                    minimum = min(cycle)
                    while x != minimum and length < horizon:
                        xn, u, j = trans[x]
                        a, b = u*(a+aa[j]) % Q, u*(b+bs[j]) % Q
                        x = xn; used += 1; length += 1
                    if length >= horizon: break
                    if length + 2 > horizon: break
                    x, u = doubles[x]
                    a, b = 2*u*a % Q, 2*u*b % Q
                    used += 2; length += 2; escapes += 1; ring = []
                    if (x & ((1 << DP_BITS)-1)) == 0:
                        key = x & ((1 << KEY_BITS)-1)
                        for oldx, olda, oldb, oldcore in table[key]:
                            if oldx != x: spurious += 1; continue
                            if oldb != b:
                                candidate = (olda-a)*pow((b-oldb) % Q, -1, Q) % Q
                                assert candidate == secret
                                recovered = True
                                if core.intersection(oldcore): good = True
                        table[key].append((x, a, b, core)); break
                ring.append(x)
                if len(ring) > 8: ring.pop(0)
                # Brent detects repeated expanded state; it does not alter a transition.
                state = (x, tuple(ring)); age += 1
                if state == checkpoint:
                    aborts += 1; break
                if age == power:
                    checkpoint = state; power *= 2; age = 0
            # Finish only this in-flight trail after T; no new starts after T.
        successes += bool(good)
        core_collisions += bool(any_core)
        recoveries += bool(recovered)
        work += used
    prob = successes/TRIALS
    exponent = T*(T-1)/(2*N)
    rows.append(dict(table_id=table_id, seed=8007+table_id, trials=TRIALS,
        success_with_core_overlap=successes, core_collision_trials=core_collisions,
        all_dp_recovery_trials=recoveries, probability=prob,
        effective_exponent_rho=-math.log1p(-prob)/exponent,
        spurious_key_matches=spurious, short_cycle_escapes=escapes,
        brent_cycle_aborts=aborts, total_main_evaluations=work, total_starts=starts))
s = sum(r['success_with_core_overlap'] for r in rows)
p = s/(TABLES*TRIALS)
# Wilson interval: descriptive pooled Bernoulli interval, NOT independent-table uncertainty.
z = 1.96; n = TABLES*TRIALS
center = (p+z*z/(2*n))/(1+z*z/n)
radius = z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)
print(json.dumps(dict(version=1, curve=dict(p=P,b=B,order=Q,classes=N,beta=beta,lambda_scalar=lam),
    parameters=dict(table_size=TABLE,dp_bits=DP_BITS,truncated_key_bits=KEY_BITS,
      tables=TABLES,trials_per_table=TRIALS,core_budget=T,trail_cap=horizon,
      cycle_ring=8,detector='Brent on expanded (point,ring) state; no detector escapes'),
    table_results=rows,total_trials=n,successes=s,pooled_probability=p,
    pooled_effective_exponent_rho=-math.log1p(-p)/(T*(T-1)/(2*N)),
    table_rho_min=min(r['effective_exponent_rho'] for r in rows),
    table_rho_max=max(r['effective_exponent_rho'] for r in rows),
    pooled_wilson95_descriptive=[center-radius,center+radius],
    scope='Toy-only complete DP pipeline, truncated-key checks and verified recovery. '
      'Toy partition/table/density differ from production. No production rate, confidence, '
      'cycle-tail, Poisson model, or instruction-cost certification.'),indent=2))

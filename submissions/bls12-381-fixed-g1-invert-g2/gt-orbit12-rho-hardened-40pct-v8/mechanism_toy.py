#!/usr/bin/env python3
"""Actual Fp12 tower diagnostic; no hash replaces either coordinate partition or DP.
Own standard-library code. Visible exponents are equality/recovery instrumentation
only, never partition inputs. Field p and subgroup r differ from production.
"""
import argparse
import json
from math import isqrt
import random
import time

PARAMETERS = [
    {'p': 70368891220219, 'r': 1549741},
    {'p': 70377839719483, 'r': 54298861},
    {'p': 70381155765259, 'r': 923346397},
    {'p': 73986475718707, 'r': 106302407641},
    {'p': 1044763899747811, 'r': 4836869906233},
]


def is_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
        if a % n == 0:
            continue
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


class Tower:
    """Fp2[w]/(w^6-(1+u)), u^2=-1, exactly the production tower form."""
    def __init__(self, parameter):
        self.p, self.r = parameter['p'], parameter['r']
        p, r = self.p, self.r
        assert p < 2**64 and r < 2**64 and is_prime(p) and is_prime(r)
        assert p % 24 == 19 and (p**4-p*p+1) % r == 0
        self.lp = [pow(p, k, r) for k in range(12)]
        assert len(set(self.lp)) == 12 and self.lp[6] == r-1
        assert pow(2, 12, r) != 1
        # Binomial irreducibility: xi is neither a square nor a cube in Fp2.
        assert self.cpow((1, 1), (p*p-1)//2) != (1, 0)
        assert self.cpow((1, 1), (p*p-1)//3) != (1, 0)
        self.one = (1,) + (0,)*11
        self.mont = (1 << 64) % p
        self.d1 = self.cpow((1, 1), (p-1)//6)
        self.zeta = self.cpow((1, 1), (p*p-1)//6)[0]
        assert self.cpow((1, 1), (p*p-1)//6)[1] == 0
        assert self.zeta*self.zeta % p == (self.zeta-1) % p
        self.gamma = []
        for k in range(12):
            d = self.cpow((1, 1), (p**k-1)//6)
            self.gamma.append([self.cpow(d, i) for i in range(6)])
        self.g = self.power((1, 0, 1, 0)+(0,)*8, (p**12-1)//r)
        assert self.g != self.one and self.power(self.g, r) == self.one
        # Fixed-base byte windows accelerate setup only. No walk lookup oracle.
        self.windows = []
        base = self.g
        for _ in range((r.bit_length()+7)//8):
            row = [self.one]
            for j in range(1, 256):
                row.append(self.mul(row[-1], base))
            self.windows.append(row)
            for _ in range(8):
                base = self.square(base)

    def cmul(self, a, b):
        p = self.p
        return ((a[0]*b[0]-a[1]*b[1]) % p,
                (a[0]*b[1]+a[1]*b[0]) % p)

    def cpow(self, a, n):
        z = (1, 0)
        while n:
            if n & 1:
                z = self.cmul(z, a)
            a = self.cmul(a, a)
            n >>= 1
        return z

    def mul(self, a, b):
        re, im = [0]*11, [0]*11
        for i in range(6):
            ar, ai = a[2*i], a[2*i+1]
            for j in range(6):
                br, bi = b[2*j], b[2*j+1]
                re[i+j] += ar*br-ai*bi
                im[i+j] += ar*bi+ai*br
        for k in range(6, 11):
            re[k-6] += re[k]-im[k]
            im[k-6] += re[k]+im[k]
        return tuple(v % self.p for k in range(6) for v in (re[k], im[k]))

    def square(self, a):
        re, im = [0]*11, [0]*11
        for i in range(6):
            ar, ai = a[2*i], a[2*i+1]
            re[2*i] += ar*ar-ai*ai
            im[2*i] += 2*ar*ai
            for j in range(i+1, 6):
                br, bi = a[2*j], a[2*j+1]
                re[i+j] += 2*(ar*br-ai*bi)
                im[i+j] += 2*(ar*bi+ai*br)
        for k in range(6, 11):
            re[k-6] += re[k]-im[k]
            im[k-6] += re[k]+im[k]
        return tuple(v % self.p for k in range(6) for v in (re[k], im[k]))

    def power(self, a, n):
        z = self.one
        while n:
            if n & 1:
                z = self.mul(z, a)
            a = self.square(a)
            n >>= 1
        return z

    def element(self, e):
        z, i = self.one, 0
        while e:
            z = self.mul(z, self.windows[i][e & 255])
            e >>= 8
            i += 1
        return z

    def frobenius(self, a, k):
        out = []
        for i, gamma in enumerate(self.gamma[k % 12]):
            c = (a[2*i], -a[2*i+1] if k & 1 else a[2*i+1])
            out.extend(self.cmul(c, gamma))
        return tuple(out)

    def canonical(self, a):
        if a == self.one:
            return 0, 0
        p, zeta = self.p, self.zeta
        x, y = a[2]*self.mont % p, a[3]*self.mont % p
        B = (x*self.d1[0]+y*self.d1[1]) % p
        xz, bz = x*zeta % p, B*zeta % p
        even = (x, xz, (xz-x) % p, -x % p, -xz % p, (x-xz) % p)
        odd = (B, bz, (bz-B) % p, -B % p, -bz % p, (B-bz) % p)
        keys = tuple(v for pair in zip(even, odd) for v in pair)
        K = min(keys)
        ties = [k for k, v in enumerate(keys) if v == K]
        if len(ties) == 1:
            return K, ties[0]
        k = min(ties, key=lambda j: tuple(v*self.mont % p for v in self.frobenius(a, j)))
        return K, k

    @staticmethod
    def mark(K, bits):
        return ((K >> 10) & ((1 << bits)-1)) == 0


class SharedMap:
    """One fixed R=1024 table shared by every walker and density comparison."""
    def __init__(self, field, rng):
        self.f, self.rng = field, rng
        self.x = rng.randrange(1, field.r)
        self.table = []
        r = field.r
        for _ in range(1024):
            a, b = rng.randrange(r), rng.randrange(r)
            y = field.element((a+self.x*b) % r)
            self.table.append([(field.frobenius(y, -k % 12),
                                a*field.lp[-k % 12] % r,
                                b*field.lp[-k % 12] % r) for k in range(12)])

    def fresh(self):
        f, r = self.f, self.f.r
        a, b = self.rng.randrange(r), self.rng.randrange(r)
        y = f.element((a+self.x*b) % r)
        K, k = f.canonical(y)
        return y, a, b, K, k

    def step(self, s):
        y, a, b, K, k = s
        f = self.f
        jump, aa, bb = self.table[K & 1023][k]
        y = f.mul(f.square(y), jump)
        K, k = f.canonical(y)
        return y, (2*a+aa) % f.r, (2*b+bb) % f.r, K, k

    def normalized(self, s):
        _, a, b, _, k = s
        f = self.f
        A, B = f.lp[k]*a % f.r, f.lp[k]*b % f.r
        # Visible exponent is instrumentation only, not an input to step/mark.
        return (A+self.x*B) % f.r, A, B

    def recovered(self, one, two):
        if one[0] != two[0]:
            return False
        den = (one[2]-two[2]) % self.f.r
        if not den:
            return False
        guess = (two[1]-one[1])*pow(den, -1, self.f.r) % self.f.r
        assert guess == self.x, 'incorrect nondegenerate aligned recovery'
        return True


def census(mapping, deadline):
    s, seen, keys = mapping.fresh(), {}, []
    while True:
        c = mapping.normalized(s)[0]
        if c in seen:
            mu, length = seen[c], len(keys)-seen[c]
            return {'mu': mu, 'lambda': length, 'rho': len(keys),
                    'marks': [sum(mapping.f.mark(K, d) for K in keys) for d in (8, 10, 12)],
                    'cycle_marks': [sum(mapping.f.mark(K, d) for K in keys[mu:]) for d in (8, 10, 12)]}
        seen[c] = len(keys)
        keys.append(s[3])
        s = mapping.step(s)
        if len(keys) % 4096 == 0 and time.monotonic() > deadline:
            raise TimeoutError('toy wall cap reached; no partial result may be called complete')


def reference_cycle(N, rng):
    """Exact lazy random-function path law, not an asymptotic rho sampler."""
    seen, y = {}, rng.randrange(N)
    while y not in seen:
        seen[y] = len(seen)
        y = rng.randrange(N)
    mu = seen[y]
    return {'mu': mu, 'lambda': len(seen)-mu, 'rho': len(seen)}


def pipeline(mapping, bits, windows, deadline):
    f, N = mapping.f, (mapping.f.r-1)//12+1
    S = isqrt(6*N//5)
    if 5*S*S < 6*N:
        S += 1
    out = dict(density_bits=bits, windows=windows, budget_per_window=S,
               tortoise_iterations=0, evaluations=0, dense_multiplication_analogues=0,
               starts=0, reports=0, floyd_detections=0, timeout_restarts=0,
               cycle_basin_iterations=0, verified_floyd_relations=0,
               verified_report_relations=0, degenerate_floyd_relations=0,
               successful_windows=0, dp_segment_lengths=[], cycle_segment_lengths=[])
    for _ in range(windows):
        reports, walkers = {}, []
        success = False
        for _ in range(2):
            s = mapping.fresh()
            walkers.append((s, s, 0))
            out['starts'] += 1
        for t in range(S):
            wi = t % 2
            tort, hare, age = walkers[wi]
            tort, hare = mapping.step(tort), mapping.step(mapping.step(hare))
            age += 1
            out['tortoise_iterations'] += 1
            out['evaluations'] += 3
            out['dense_multiplication_analogues'] += 6
            if f.mark(tort[3], bits):
                norm = mapping.normalized(tort)
                for previous in reports.get(norm[0], []):
                    if mapping.recovered(norm, previous):
                        out['verified_report_relations'] += 1
                        success = True
                reports.setdefault(norm[0], []).append(norm)
                out['reports'] += 1
                out['dp_segment_lengths'].append(age)
            elif tort[3] == hare[3] and mapping.normalized(tort)[0] == mapping.normalized(hare)[0]:
                out['floyd_detections'] += 1
                out['cycle_basin_iterations'] += age
                out['cycle_segment_lengths'].append(age)
                if mapping.recovered(mapping.normalized(tort), mapping.normalized(hare)):
                    out['verified_floyd_relations'] += 1
                    success = True
                else:
                    out['degenerate_floyd_relations'] += 1
            elif age >= 128*(1 << bits):
                out['timeout_restarts'] += 1
            else:
                walkers[wi] = tort, hare, age
                continue
            s = mapping.fresh()
            walkers[wi] = s, s, 0
            out['starts'] += 1
        out['successful_windows'] += success
        if time.monotonic() > deadline:
            raise TimeoutError('toy pipeline wall cap reached')
    return out


def experiment(seed=83, profile='full', wall_seconds=780, progress=None):
    start = time.monotonic()
    deadline = start+wall_seconds
    rng, reference_rng = random.Random(seed), random.Random(seed ^ 0x5A17)
    counts = [128, 64, 32, 12, 4] if profile == 'full' else [1]
    result = {'schema_version':1, 'seed':seed, 'profile':profile,
              'partition':'minimum real w coefficient over all 12 actual Frobenius conjugates; low 10 bits',
              'DP_predicate':'((K >> 10) & ((1 << density_bits)-1)) == 0',
              'densities_bits':[8,10,12], 'R':1024, 'fields':[], 'pipeline':[],
              'equivariance_failures':0, 'incorrect_recoveries':0,
              'cycles_observed':0, 'completed':False,
              'scope':'Actual order-12 Frobenius on the production-form Fp12 tower with small prime-order subgroups. Visible exponents certify equality and aligned recovery only. No production-parameter execution or production probability theorem.'}
    for fi, trials in enumerate(counts):
        f = Tower(PARAMETERS[fi])
        N = (f.r-1)//12+1
        row = {'p':f.p, 'r':f.r, 'N':N, 'lambda':f.lp[1],
               'independent_tables':trials, 'first_cycles':[], 'reference_cycles':[],
               'basin_trials':[], 'uniform_coordinate_samples':0,
               'uniform_mark_counts':[0,0,0], 'partition_counts':[0]*1024}
        nsamples = 32768 if profile == 'full' else 1024
        for _ in range(nsamples):
            K, k = f.canonical(f.element(rng.randrange(f.r)))
            row['uniform_coordinate_samples'] += 1
            row['partition_counts'][K & 1023] += 1
            for i, d in enumerate((8, 10, 12)):
                row['uniform_mark_counts'][i] += f.mark(K, d)
        for ti in range(trials):
            mapping = SharedMap(f, rng)
            probe = mapping.fresh()
            expected = mapping.normalized(mapping.step(probe))[0]
            for j in range(12):
                y, a, b, _, _ = probe
                y = f.frobenius(y, j)
                K, k = f.canonical(y)
                shifted = y, a*f.lp[j] % f.r, b*f.lp[j] % f.r, K, k
                result['equivariance_failures'] += mapping.normalized(mapping.step(shifted))[0] != expected
            first = census(mapping, deadline)
            result['cycles_observed'] += 1
            row['first_cycles'].append(first)
            row['reference_cycles'].append(reference_cycle(N, reference_rng))
            # Long restart/basin probes: >=32 actual cycles for each such trial.
            long_trial = fi <= 1 and ti < (4 if profile == 'full' else 1)
            if long_trial:
                cycles = [first]
                for _ in range(31):
                    cycles.append(census(mapping, deadline))
                    result['cycles_observed'] += 1
                row['basin_trials'].append({'table':ti, 'cycles':cycles})
            if ti < (4 if fi == 0 else 2 if fi == 1 else 0):
                for d in (8,10,12):
                    p = pipeline(mapping, d, 128 if profile == 'full' else 8, deadline)
                    p.update(field_index=fi, table=ti, N=N)
                    result['pipeline'].append(p)
            if progress is not None:
                progress({'field':fi,'table':ti,'tables':trials,'elapsed_seconds':time.monotonic()-start})
        result['fields'].append(row)
    result['completed'] = True
    result['elapsed_seconds'] = time.monotonic()-start
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--profile', choices=['smoke','full'], default='full')
    parser.add_argument('--seed', type=int, default=83)
    parser.add_argument('--wall-seconds', type=int, default=780)
    parser.add_argument('--output')
    args = parser.parse_args()
    import sys
    out = experiment(args.seed, args.profile, args.wall_seconds,
                     lambda row: print(json.dumps(row), file=sys.stderr, flush=True))
    text = json.dumps(out, separators=(',', ':'))+'\n'
    if args.output:
        from pathlib import Path
        target = Path(args.output)
        # Two JSON files preserve the same capture while respecting 64KiB/file.
        companion = target.with_name('mechanism-pipeline.json')
        companion.write_text(json.dumps({'pipeline':out['pipeline']}, separators=(',', ':'))+'\n')
        main = {k:v for k,v in out.items() if k != 'pipeline'}
        main['pipeline_file'] = companion.name
        target.write_text(json.dumps(main, separators=(',', ':'))+'\n')
        print(json.dumps({'completed':out['completed'], 'cycles':out['cycles_observed'],
                          'elapsed_seconds':out['elapsed_seconds'], 'output':args.output}))
    else:
        print(text)

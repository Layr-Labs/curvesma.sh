""
from math import exp, log, log2, log1p, ceil, floor, sqrt
from random import Random
from collections import Counter
from array import array
from common12 import scaled_log

LN2 = log(2)

def primes_upto(n):
    sieve = bytearray(b'\1')*(n+1)
    sieve[:2] = b'\0\0'
    for p in range(2, int(sqrt(n))+1):
        if sieve[p]:
            sieve[p*p::p] = b'\0'*((n-p*p)//p+1)
    return [p for p in range(2, n+1) if sieve[p]]

def coefficients(weights, normalizer, length):
    ""
    terms = [(j, j*v) for j, v in enumerate(weights) if j and abs(v)>=1e-24]
    out = [exp(-normalizer)]+[0.0]*length
    for n in range(1, length+1):
        value = 0.0
        for j, w in terms:
            if j > n:
                break
            value += w*out[n-j]
        out[n] = value/n
        if value < -1e-13:
            raise ArithmeticError('negative coefficient beyond roundoff tolerance')
    return out

def convolution(a, b, length):
    out = [0.0]*(length+1)
    for i, x in enumerate(a[:length+1]):
        if x:
            for j in range(min(len(b)-1, length-i)+1):
                out[i+j] += x*b[j]
    return out

class Retention:
    def __init__(self, B_bits, X_bits, *, prime_min=256, step=.25,
                 max_bits=1024, tilt=.05, explicit_bits=12, power_terms=8):
        if not 1 < prime_min <= 2**X_bits <= 2**B_bits:
            raise ValueError('invalid sieve/factor-base range')
        if not 0<step<=1 or max_bits <= B_bits or tilt >= LN2:
            raise ValueError('invalid grid/tilt')
        if prime_min > 2**min(float(explicit_bits), B_bits):

            raise ValueError('explicit prime range must cover prime_min')
        self.B, self.X, self.minimum, self.step = B_bits, X_bits, prime_min, step
        self.tilt, self.length = tilt, ceil(max_bits/step)

        self.names = ('small','sieve_geo','large_geo','sieve_sf','large_sf')
        weights = {k:[0.0]*(self.length+1) for k in self.names}
        totals = Counter()
        cutoff = min(float(explicit_bits), B_bits)
        self.explicit_cutoff_bits = cutoff

        def add_prime(p):
            b = log2(p)
            base_index = round(b/step)
            q = exp(tilt*base_index*step-LN2*b)
            if not 0<q<1:
                raise ValueError('tilt makes prime series divergent')
            if p < prime_min:
                groups = [('small',False)]
            else:
                prefix = 'sieve' if b < X_bits else 'large'
                groups = [(prefix+'_geo',False),(prefix+'_sf',True)]
            for name, squarefree in groups:
                totals[name] += log1p(q) if squarefree else -log1p(-q)
                qk = q
                for k in range(1, self.length+1):
                    index = k*base_index
                    if index > self.length:
                        break
                    sign = -1 if squarefree and k%2 == 0 else 1
                    weights[name][index] += sign*qk/k
                    qk *= q
                    if qk < 1e-25:
                        break
        for p in primes_upto(int(2**cutoff)-1):
            add_prime(p)

        def continuum(left, right, prefix):
            if left >= right:
                return

            width = step/4
            cells = ceil((right-left)/width)
            width = (right-left)/cells
            for i in range(cells):
                b = left+(i+.5)*width
                base_index = round(b/step)
                for k in range(1,power_terms+1):
                    amount = exp(k*tilt*base_index*step-(k-1)*LN2*b)*width/(k*b)
                    index = k*base_index
                    geo, sf = prefix+'_geo', prefix+'_sf'
                    sign = -1 if k%2 == 0 else 1
                    totals[geo] += amount; totals[sf] += sign*amount
                    if index <= self.length:
                        weights[geo][index] += amount
                        weights[sf][index] += sign*amount
        continuum(cutoff,min(X_bits,B_bits),'sieve')
        continuum(max(cutoff,X_bits),B_bits,'large')
        self.weights, self.normalizers = weights, dict(totals)
        self.pdf = {name:coefficients(weights[name],totals[name],self.length)
                    for name in ('sieve_sf','sieve_geo')}

        missing_weights = [x+y for x,y in zip(weights['small'],weights['large_sf'])]
        self.missing_sf = coefficients(missing_weights,totals['small']+totals['large_sf'],self.length)
        missing_geo_weights = [x+y for x,y in zip(weights['small'],weights['large_geo'])]
        self.missing_geo = coefficients(missing_geo_weights,totals['small']+totals['large_geo'],self.length)
        all_weights = [sum(weights[k][j] for k in ('small','sieve_geo','large_geo'))
                       for j in range(self.length+1)]
        self.denominator = coefficients(all_weights,
                                        sum(totals[k] for k in ('small','sieve_geo','large_geo')),
                                        self.length)
        self.selection_normalizer = exp(totals['sieve_sf']+totals['large_sf']
                                        -totals['sieve_geo']-totals['large_geo'])
        self.log_partition = totals['small']+totals['sieve_geo']+totals['large_geo']
        self._conditional_cdfs = {}

    def _mass(self,n,stop,squarefree):
        if stop<0:
            return 0.0
        key=(n,squarefree)
        if key not in self._conditional_cdfs:
            missing=self.missing_sf if squarefree else self.missing_geo
            sieve=self.pdf['sieve_sf' if squarefree else 'sieve_geo']
            values,total=array('d'),0.0
            for d in range(n+1):
                total+=missing[d]*sieve[n-d]
                values.append(total)
            self._conditional_cdfs[key]=values
        return self._conditional_cdfs[key][min(n,stop)]

    def at(self, norm_bits, slack_bits, *, scale=16, guard_bits=0.0):
        if norm_bits <= 0 or slack_bits < 0 or norm_bits/self.step > self.length:
            raise ValueError('norm or slack outside model grid')
        n = round(norm_bits/self.step)

        count_cap = floor(ceil(norm_bits)/log2(self.minimum))
        rounding_budget = (count_cap+1)/scale
        allowed = slack_bits-rounding_budget-guard_bits
        stop = min(n,floor(allowed/self.step))
        numerator = self._mass(n,stop,True)
        denominator = self.denominator[n]
        if denominator <= 0:
            raise ArithmeticError('zero model denominator')
        retained = self.selection_normalizer*numerator/denominator
        if not -1e-10 <= retained <= 1+1e-9:
            raise ArithmeticError('retention outside probability range')
        smooth_probability = exp(self.log_partition-self.tilt*n*self.step)*denominator/self.step/LN2
        return dict(retention=max(0.0,min(1.0,retained)),
                    rounding_budget_bits=rounding_budget,
                    omitted_mass_limit_bits=allowed,
                    smooth_probability_model=smooth_probability)

    def survivor_upper(self,norm_bits,slack_bits,preliminary_extra_bits=0):
        ""
        n=round(norm_bits/self.step)
        if not 0<n<=self.length or slack_bits<0 or preliminary_extra_bits<0:
            raise ValueError('invalid survivor range')
        stop=min(n,floor((slack_bits+preliminary_extra_bits)/self.step))
        numerator=self._mass(n,stop,False)
        return numerator/self.denominator[n]

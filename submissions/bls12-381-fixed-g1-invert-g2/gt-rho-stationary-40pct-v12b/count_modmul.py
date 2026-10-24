#!/usr/bin/env python3
"""Differential test + instruction count for the counted renormalisation routine.

1. builds a small harness that calls mulmod_r on random (a, cs) pairs and prints
   them, so the C result can be compared against bigint arithmetic in Python;
2. compiles mulmod_r.c for arm64 with clang -O2 and counts the instructions of
   the emitted body with otool -tvV, writing both numbers to modmul-count.json.
"""
import json
import os
import random
import subprocess
import sys

D = os.path.dirname(os.path.abspath(__file__))
R = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
P = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab

HARNESS = r'''
#include <stdint.h>
#include <stdio.h>
void mulmod_r(uint64_t out[4], const uint64_t a[4], const uint64_t cs[4]);
int main(void){
    uint64_t a[4], c[4], o[4];
    while (scanf("%llx %llx %llx %llx %llx %llx %llx %llx",
                 (unsigned long long*)&a[0],(unsigned long long*)&a[1],
                 (unsigned long long*)&a[2],(unsigned long long*)&a[3],
                 (unsigned long long*)&c[0],(unsigned long long*)&c[1],
                 (unsigned long long*)&c[2],(unsigned long long*)&c[3]) == 8) {
        mulmod_r(o, a, c);
        printf("%016llx %016llx %016llx %016llx\n",(unsigned long long)o[0],(unsigned long long)o[1],
               (unsigned long long)o[2],(unsigned long long)o[3]);
    }
    return 0;
}
'''


def limbs(x, n=4):
    return [(x >> (64 * i)) & ((1 << 64) - 1) for i in range(n)]


def main():
    open(os.path.join(D, 'harness.c'), 'w').write(HARNESS)
    subprocess.run(['clang', '-O2', '-arch', 'arm64', '-c',
                    os.path.join(D, 'mulmod_r.c'), '-o', os.path.join(D, 'mulmod_r.o')],
                   check=True)
    subprocess.run(['clang', '-O2', '-arch', 'arm64',
                    os.path.join(D, 'harness.c'), os.path.join(D, 'mulmod_r.o'),
                    '-o', os.path.join(D, 'harness')], check=True)

    rng = random.Random(20261024)
    cases = []
    # p^k * 2^256 mod r, the twelve pre-scaled constants
    psc = [pow(P, k, R) * pow(2, 256, R) % R for k in range(12)]
    for k in range(12):
        for _ in range(20):
            a = rng.randrange(0, R) if _ else (R - 1 if k % 3 == 0 else 0)
            cases.append((a, psc[k]))
    # edge cases
    for a in (0, 1, R - 1, R - 2, 2):
        for k in (0, 6, 11):
            cases.append((a, psc[k]))
    inp = '\n'.join(' '.join('%016x' % v for v in limbs(a) + limbs(c)) for a, c in cases)
    out = subprocess.run([os.path.join(D, 'harness')], input=inp,
                         capture_output=True, text=True, check=True).stdout.strip().split('\n')
    bad = 0
    for (a, c), line in zip(cases, out):
        got = int(''.join(reversed(line.split())), 16)
        want = a * c * pow(2, -256, R) % R
        if got != want:
            bad += 1
    n_cases = len(cases)

    txt = subprocess.run(['otool', '-tvV', os.path.join(D, 'mulmod_r.o')],
                         capture_output=True, text=True, check=True).stdout
    body = []
    inside = False
    for ln in txt.split('\n'):
        if '_mulmod_r:' in ln:
            inside = True
            continue
        if inside:
            if ln.strip() == '':
                break
            body.append(ln)
    n_instr = len(body)

    # also compile at -O0 to show the priced row is the -O2 scheduled body
    subprocess.run(['clang', '-O0', '-arch', 'arm64', '-c',
                    os.path.join(D, 'mulmod_r.c'), '-o', os.path.join(D, 'mulmod_r_O0.o')],
                   check=True)
    txt0 = subprocess.run(['otool', '-tvV', os.path.join(D, 'mulmod_r_O0.o')],
                          capture_output=True, text=True, check=True).stdout
    body0, inside = [], False
    for ln in txt0.split('\n'):
        if '_mulmod_r:' in ln:
            inside = True
            continue
        if inside:
            if ln.strip() == '':
                break
            body0.append(ln)

    res = dict(
        routine='mulmod_r(a, cs) = REDC(a * (p^k 2^256 mod r)) = a * p^k mod r',
        source='mulmod_r.c', compiler='clang -O2 -arch arm64',
        instructions_O2=n_instr, instructions_O0=len(body0),
        differential_cases=n_cases, differential_mismatches=bad,
        covers='twelve k values 0..11 x 20 random labels, plus a in {0,1,2,r-1,r-2}',
        renorm_row_per_step=2 * n_instr,
        note='labelled row = 2 x this count, one multiply per label (a and b)',
    )
    json.dump(res, open(os.path.join(D, 'modmul-count.json'), 'w'), indent=1)
    print(json.dumps(res, indent=1))
    return 0 if bad == 0 else 1


if __name__ == '__main__':
    sys.exit(main())

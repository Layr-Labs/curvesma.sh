"""Python replica of cheon_step.c: runs the same non-reference step on standard integers.

  python3 -I replica.py <target-profile.json> <inputs.txt> <n> [kernel-dump.txt]

Prints the final state in the kernel's dump format. With a kernel dump, it also checks
equality and exits non-zero on a mismatch. The Montgomery limbs (R = 2^384) are reproduced
where the kernel's choices depend on them: the multiplier index, the least-of-three
selection and the distinguished-point test.
"""
import json
import sys

prof = json.load(open(sys.argv[1]))
p = int(prof["field_modulus"], 16)
r = int(prof["subgroup_order"], 16)
lines = [ln.strip() for ln in open(sys.argv[2]) if ln.strip()]
n = int(sys.argv[3])
WINDOWS, TBL, NMULT, NMU = 16, 64, 64, 6
MASK64 = (1 << 64) - 1
pos = 0


def take(k):
    global pos
    v = [int(t, 16) for t in lines[pos:pos + k]]
    pos += k
    return v


beta = take(1)[0]
x, y = take(2)
c = take(1)[0]
mult = take(NMULT)
mu = take(NMU)
comb = [[None] * TBL for _ in range(WINDOWS)]
for k in range(WINDOWS):
    for i in range(TBL):
        comb[k][i] = take(2)


def mont_limbs(v):
    m = v * (1 << 384) % p
    return [(m >> (64 * i)) & MASK64 for i in range(6)]


def less(a, b):
    for i in range(5, -1, -1):
        if a[i] != b[i]:
            return a[i] < b[i]
    return False


dps = 0
for _ in range(n):
    j = mont_limbs(x)[1] & (NMULT - 1)
    c2 = c * mult[j] % r
    carry = 0
    acc = None
    for k in range(WINDOWS):
        v = (c2 >> (16 * k)) & 0xFFFF
        v += carry
        carry = 1 if v >= 0x8000 else 0
        neg = carry
        if neg:
            v = 0x10000 - v
        ex, ey = comb[k][v & (TBL - 1)]
        if neg:
            ey = (p - ey) % p
        acc = (ex, ey)
    ax, ay = acc
    bx = beta * ax % p
    nx = (-(ax + bx)) % p
    best, sel = ax, 0
    if less(mont_limbs(bx), mont_limbs(best)):
        best, sel = bx, 1
    if less(mont_limbs(nx), mont_limbs(best)):
        best, sel = nx, 2
    sign = ay & 1
    if sign:
        ay = (p - ay) % p
    c = c2 * mu[sel * 2 + sign] % r
    x, y = best, ay
    bl = mont_limbs(x)
    if (bl[0] | (bl[1] & ((1 << 26) - 1))) == 0:
        dps += 1

dump = "%096x\n%096x\n%064x\n" % (x, y, c)
sys.stderr.write("steps=%d dps=%d\n" % (n, dps))
sys.stdout.write(dump)
if len(sys.argv) > 4:
    kern = open(sys.argv[4]).read()
    if kern != dump:
        sys.stderr.write("MISMATCH with kernel dump\n")
        sys.exit(1)
    sys.stderr.write("kernel dump matches replica\n")

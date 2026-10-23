"""Python replica of cheon_step_v3.c: the same walk on standard integers and the BLS12-381 group law.

  python3 -B replica_v3.py <target-profile.json> <inputs.txt> <n> <prod|dump|dp4|cap1> [kernel-dump]

Runs n walk points round-robin over 256 records (even records tame, odd wild), with the same
digit, table, canonical-rotation, restart and distinguished-point rules as the kernel. Prints
the final records in the kernel's dump format. With a kernel dump, it checks equality.
"""
import sys

import ec

prof_path, inputs_path, n, cfg = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
kern_path = sys.argv[5] if len(sys.argv) > 5 else None
p, r, _P0 = ec.load_profile(prof_path)
lines = [ln.strip() for ln in open(inputs_path) if ln.strip()]
pos = 0


def take(k):
    global pos
    v = [int(t, 16) for t in lines[pos:pos + k]]
    pos += k
    return v


WINDOWS, TBLN, NMULT, NMU, NSTART, WALKS = 16, 1 << 15, 64, 6, 16, 256
beta = take(1)[0]
P0 = tuple(take(2))
Q = tuple(take(2))
mult = take(NMULT)
mu = take(NMU)
starts = take(NSTART)
dp_bits = 4 if cfg == "dp4" else 90
cap = 1 if cfg == "cap1" else 1 << 100
_base_cache = {}


def base_k(fam, k):
    """2^(16k) * B for the family's base B (P0 for tame, Q for wild)."""
    key = (fam, k)
    if key not in _base_cache:
        B = Q if fam else P0
        _base_cache[key] = ec.mul(1 << (16 * k), B, p)
    return _base_cache[key]


def entry(fam, k, m):
    """Table entry for digit magnitude m >= 1: m * 2^(16k) * B (the kernel's entry m-1)."""
    return ec.mul(m, base_k(fam, k), p)


def comb_eval(c, fam):
    carry = 0
    acc = None
    for k in range(WINDOWS):
        v = (c >> (16 * k)) & 0xFFFF
        v += carry
        carry = 1 if v >= 0x8000 else 0
        neg = carry
        if neg:
            v = 0x10000 - v
        if v == 0:
            continue
        ex, ey = entry(fam, k, v)
        if neg:
            ey = (p - ey) % p
        acc = (ex, ey) if acc is None else ec.add(acc, (ex, ey), p)
    return acc


def canon(pt):
    x, y = pt
    bx = beta * x % p
    nx = (-(x + bx)) % p
    best, sel = x, 0
    if ec.less_mont(bx, best, p):
        best, sel = bx, 1
    if ec.less_mont(nx, best, p):
        best, sel = nx, 2
    sign = y & 1
    if sign:
        y = (p - y) % p
    return (best, y), sel * 2 + sign


recs = []
for w in range(WALKS):
    fam = w & 1
    c = starts[w & (NSTART - 1)] * mult[(w >> 4) & (NMULT - 1)] % r
    newpt, idx = canon(comb_eval(c, fam))
    recs.append([newpt, c * mu[idx] % r, 0, fam])
dps = 0
restarts = 0
for i in range(n):
    rec = recs[i % WALKS]
    pt, c, since, fam = rec
    j = (ec.mont_limbs(pt[0], p)[1] >> 26) & (NMULT - 1)
    c2 = c * mult[j] % r
    newpt, idx = canon(comb_eval(c2, fam))
    c3 = c2 * mu[idx] % r
    since += 1
    if since >= cap:
        cs = starts[restarts % NSTART]
        restarts += 1
        newpt, idx = canon(comb_eval(cs, fam))
        c3 = cs * mu[idx] % r
        since = 0
    rec[0], rec[1], rec[2] = newpt, c3, since
    m = ec.mont_limbs(newpt[0], p)
    if ((m[0] | (m[1] << 64)) & ((1 << dp_bits) - 1)) == 0:
        dps += 1

sys.stderr.write("steps=%d dps=%d restarts=%d\n" % (n, dps, restarts))
dump = "".join("%096x\n%096x\n%064x\n%032x\n" % (rec[0][0], rec[0][1], rec[1], rec[2]) for rec in recs)
sys.stdout.write(dump)
if kern_path:
    if open(kern_path).read() != dump:
        sys.stderr.write("MISMATCH with kernel dump\n")
        sys.exit(1)
    sys.stderr.write("kernel dump matches replica\n")

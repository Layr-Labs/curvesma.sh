"""Input file for cheon_step_full.c (hex, one value per line; order matches the C reader).

  python3 -I gen_inputs_full.py <target-profile.json> <seed> > inputs.txt

Order: beta, x0, y0, c0, 64 multipliers, 6 mu6 elements, comb table for P0 (16 x 64 points),
comb table for A_d (16 x 64 points). Values are arithmetic test data: the kernel never checks
curve membership, and the instruction count does not depend on the values. Multipliers and mu6
come from the group orders so that they lie in the right ranges. The output is deterministic
in the seed.
"""
import json
import random
import sys

prof = json.load(open(sys.argv[1]))
p = int(prof["field_modulus"], 16)
r = int(prof["subgroup_order"], 16)
D = 1 << 20
rng = random.Random(int(sys.argv[2]))
WINDOWS, TBL, NMULT, NMU = 16, 64, 64, 6


def cube_root_of_unity(m, rng_):
    while True:
        g = rng_.randrange(2, m - 1)
        c = pow(g, (m - 1) // 3, m)
        if c != 1:
            return c


beta = cube_root_of_unity(p, rng)              # Fp cube root; the endomorphism is (x,y)->(beta x, y)
lam = cube_root_of_unity(r, rng)               # matching eigenvalue in Fr
mu = [1, r - 1, lam, r - lam, lam * lam % r, r - lam * lam % r]   # index sel*2+sign -> lam^sel (-1)^sign
x0 = rng.randrange(1, p)
y0 = rng.randrange(1, p)
c0 = pow(rng.randrange(2, r - 1), D, r)        # element of H (d-th power)
mult = [pow(rng.randrange(2, r - 1), D, r) for _ in range(NMULT)]   # walk multipliers in H

out = ["%096x" % beta, "%096x" % x0, "%096x" % y0, "%064x" % c0]
out += ["%064x" % g for g in mult]
out += ["%064x" % m for m in mu]
for _ in range(2 * WINDOWS * TBL):             # P0 table, then A_d table
    out.append("%096x" % rng.randrange(1, p))
    out.append("%096x" % rng.randrange(1, p))
print("\n".join(out))

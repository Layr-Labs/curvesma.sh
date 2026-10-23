"""Input file for cheon_step_v3.c (hex, one value per line; order matches the C reader).

  python3 -B gen_inputs_v3.py <target-profile.json> <seed> > inputs.txt

Order: beta; P0 (x, y); Q = [e]P0 (x, y) for a random e; 64 walk multipliers (elements of H);
6 mu6 elements; 16 restart exponents (elements of H). The kernel builds the comb tables itself.
beta and lambda are matched so that (beta x, y) = [lambda](x, y) on G1.
"""
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ec  # noqa: E402

p, r, P0 = ec.load_profile(sys.argv[1])
rng = random.Random(int(sys.argv[2]))
D = 1 << 20
NMULT, NSTART = 64, 16

beta, lam = ec.matched_beta_lambda(P0, p, r)
mu = [1, r - 1, lam, r - lam, lam * lam % r, r - lam * lam % r]    # index sel*2+sign
e = rng.randrange(2, r - 1)
Q = ec.mul(e, P0, p)
mult = [pow(rng.randrange(2, r - 1), D, r) for _ in range(NMULT)]
starts = [pow(rng.randrange(2, r - 1), D, r) for _ in range(NSTART)]

out = ["%096x" % beta, "%096x" % P0[0], "%096x" % P0[1], "%096x" % Q[0], "%096x" % Q[1]]
out += ["%064x" % m for m in mult]
out += ["%064x" % m for m in mu]
out += ["%064x" % s for s in starts]
print("\n".join(out))

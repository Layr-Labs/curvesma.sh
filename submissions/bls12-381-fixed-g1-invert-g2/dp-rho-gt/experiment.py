"""Toy parallel distinguished-point Pollard rho for discrete log (no automorphism).

Standard van Oorschot-Wiener DP collision search: many trails walk the
deterministic map F(X)=X*M_{h(X)}; each trail runs to a distinguished point (DP),
reports (element, A, B), and restarts; a collision between two trails' DPs solves
the log. ONE group multiplication per step (the reference unit) -- no tortoise/hare
3x factor. In production the DP reports are routed to many storage nodes (each a
separate <2^64-byte machine); the toy stores them in one dict and only measures
the walk's collision behaviour within a hard step budget.

Group: prime-order-l subgroup of F_q* (cyclic, standing in for GT). Success =
usable collision found within the hard budget T=ceil(budget_factor*k0) TOTAL
multiplications, k0=sqrt(2*l*ln(1/0.61)). Toy scale only.
"""
import json, math, random, sys


def is_prime(n):
    if n < 2: return False
    for sp in (2,3,5,7,11,13,17,19,23,29,31,37):
        if n % sp == 0: return n == sp
    d,s=n-1,0
    while d%2==0: d//=2; s+=1
    for a in (2,3,5,7,11,13,17,19,23,29,31,37):
        x=pow(a,d,n)
        if x in (1,n-1): continue
        for _ in range(s-1):
            x=x*x%n
            if x==n-1: break
        else: return False
    return True


def setup(bits, rng):
    while True:
        l = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if not is_prime(l): continue
        for c in range(2, 4000, 2):
            q = c*l + 1
            if is_prime(q):
                g = 1
                while g == 1: g = pow(rng.randrange(2, q-1), (q-1)//l, q)
                return l, q, g


def trial(l, q, g, rng, R, theta_bits, walks, seg_cap, max_steps, counters):
    x = rng.randrange(1, l); h = pow(g, x, q)
    table = [(pow(g, (a:=rng.randrange(l)), q) * pow(h, (b:=rng.randrange(l)), q) % q, a, b) for _ in range(R)]
    mask = (1 << theta_bits) - 1
    seg_limit = seg_cap << theta_bits
    seen = {}

    def fresh(prev):
        if prev is None:
            a, b = rng.randrange(l), rng.randrange(l); u = pow(g, a, q)*pow(h, b, q) % q
        else:
            pu, pa, pb = prev; jj = rng.randrange(R); m, ma, mb = table[jj]
            u = pu*m % q; a, b = (pa+ma)%l, (pb+mb)%l; counters['mults'] += 1
        return [u, a, b, 0]

    state = [fresh(None) for _ in range(walks)]
    while counters['mults'] < max_steps:
        for w in state:
            if counters['mults'] >= max_steps: return False
            u, a, b, n = w
            j = (u * 0x9E3779B1 >> 7) % R
            m, ma, mb = table[j]
            nu, na, nb = u*m % q, (a+ma)%l, (b+mb)%l
            counters['mults'] += 1; n += 1
            if (nu >> 3) & mask == 0:                 # distinguished point
                if nu in seen:
                    pa, pb = seen[nu]
                    if (pb-nb) % l:
                        xr = (na-pa)*pow(pb-nb, l-2, l) % l
                        if pow(g, xr, q) == h: return True
                else: seen[nu] = (na, nb)
                w[:] = fresh((nu, na, nb)); continue
            if n >= seg_limit:                        # abandon stuck trail
                counters['abandon'] += 1; w[:] = fresh((nu, na, nb)); continue
            w[:] = [nu, na, nb, n]
    return False


def main():
    args = json.load(sys.stdin)
    bits = int(args.get("bits", 24)); trials = int(args.get("trials", 150))
    seed = int(args.get("seed", 1)); R = int(args.get("R", 2048))
    theta_bits = int(args.get("theta_bits", 4)); walks = int(args.get("walks", 8))
    seg_cap = int(args.get("seg_cap", 256)); bf = float(args.get("budget_factor", 1.3))
    rng = random.Random(seed); l, q, g = setup(bits, rng)
    k0 = math.sqrt(2*l*math.log(1/0.61)); max_steps = int(math.ceil(bf*k0))
    ok = abn = 0
    for _ in range(trials):
        c = {'mults': 0, 'abandon': 0}
        if trial(l, q, g, rng, R, theta_bits, walks, seg_cap, max_steps, c): ok += 1
        abn += c['abandon']
    print(json.dumps({
        "l": l, "trials": trials, "R": R, "theta_bits": theta_bits, "walks": walks,
        "budget_factor": bf, "max_steps": max_steps,
        "max_steps_over_sqrt_l": max_steps/math.sqrt(l),
        "success_within_budget": ok/trials,
        "mean_abandons_per_trial": abn/trials,
        "random_mapping_lb": 1 - 0.61**(bf*bf),
    }, indent=1))


if __name__ == "__main__":
    main()

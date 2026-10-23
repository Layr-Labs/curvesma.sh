"""Toy SINGLE-TRAIL distinguished-point rho (self-collision) for discrete log.

One continuous trajectory X_{i+1}=F(X_i)=X_i*M_{h(X_i)} walks; at each
distinguished point it stores (element -> (A,B)); when a DP element REPEATS (the
trajectory has re-entered its own cycle), the two exponent pairs solve the log.
ONE group multiplication per step -- no tortoise/hare factor and no cross-trail
merge analysis. In production the DP store is distributed over many <2^64-byte
nodes (routing charged); the toy uses one dict and measures the self-collision
behaviour within a hard step budget.

Group: prime-order-l subgroup of F_q* (cyclic, standing in for GT). Success =
usable self-collision detected within T=ceil(budget_factor*k0)+L_cap total
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


def trial(l, q, g, rng, R, theta_bits, max_steps, counters):
    """One single-trail DP rho; return True iff a usable self-collision is found."""
    x = rng.randrange(1, l); h = pow(g, x, q)
    table = [(pow(g, (a:=rng.randrange(l)), q) * pow(h, (b:=rng.randrange(l)), q) % q, a, b) for _ in range(R)]
    mask = (1 << theta_bits) - 1
    seen = {}
    A, B = rng.randrange(l), rng.randrange(l)
    u = pow(g, A, q) * pow(h, B, q) % q
    steps = 0
    while steps < max_steps:
        j = (u * 0x9E3779B1 >> 7) % R
        m, ma, mb = table[j]
        u = u*m % q; A = (A+ma)%l; B = (B+mb)%l
        steps += 1; counters['mults'] += 1
        if (u >> 3) & mask == 0:                    # distinguished point
            if u in seen:                           # self-collision: trajectory cycled
                pa, pb = seen[u]
                if (pb-B) % l:
                    xr = (A-pa)*pow(pb-B, l-2, l) % l
                    if pow(g, xr, q) == h: return True
                return False                        # unusable (B==pb); rare
            seen[u] = (A, B)
    return False


def main():
    args = json.load(sys.stdin)
    bits = int(args.get("bits", 24)); trials = int(args.get("trials", 150))
    seed = int(args.get("seed", 1)); R = int(args.get("R", 2048))
    theta_bits = int(args.get("theta_bits", 4)); bf = float(args.get("budget_factor", 1.3))
    rng = random.Random(seed); l, q, g = setup(bits, rng)
    k0 = math.sqrt(2*l*math.log(1/0.61))
    L_cap = 20 << theta_bits
    max_steps = int(math.ceil(bf*k0)) + L_cap        # hard budget T + L_cap
    ok = 0
    for _ in range(trials):
        c = {'mults': 0}
        if trial(l, q, g, rng, R, theta_bits, max_steps, c): ok += 1
    print(json.dumps({
        "l": l, "trials": trials, "R": R, "theta_bits": theta_bits,
        "budget_factor": bf, "max_steps": max_steps,
        "success_within_budget": ok/trials,
        "random_mapping_lb_P_rho_le_T": 1 - 0.61**(bf*bf),
    }, indent=1))


if __name__ == "__main__":
    main()

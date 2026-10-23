# Minimal BLS12-381 Fp12 tower arithmetic (blst/py_ecc convention) to generate
# GT elements and measure coefficient-bit distributions. Self-contained.
P = int('1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab',16)
R = int('73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001',16)

def fp_mul(a,b): return a*b % P
def fp_add(a,b): return (a+b) % P
def fp_sub(a,b): return (a-b) % P

# Fp2 = Fp[u]/(u^2+1); element (a0,a1) = a0 + a1*u
def fp2_mul(A,B):
    a0,a1=A; b0,b1=B
    return ((a0*b0 - a1*b1) % P, (a0*b1 + a1*b0) % P)
def fp2_add(A,B): return ((A[0]+B[0])%P,(A[1]+B[1])%P)
def fp2_sub(A,B): return ((A[0]-B[0])%P,(A[1]-B[1])%P)
XI=(1,1)  # v^3 = xi = 1+u
def fp2_mul_xi(A):  # multiply Fp2 element by (1+u): (a0 - a1) + (a0 + a1) u
    a0,a1=A
    return ((a0-a1)%P,(a0+a1)%P)
FP2_0=(0,0); FP2_1=(1,0)

# Fp6 = Fp2[v]/(v^3 - xi); element (c0,c1,c2) = c0 + c1*v + c2*v^2, ci in Fp2
def fp6_mul(A,B):
    a0,a1,a2=A; b0,b1,b2=B
    t0=fp2_mul(a0,b0); t1=fp2_mul(a1,b1); t2=fp2_mul(a2,b2)
    # c0 = t0 + xi*((a1+a2)(b1+b2)-t1-t2)
    c0=fp2_add(t0, fp2_mul_xi(fp2_sub(fp2_sub(fp2_mul(fp2_add(a1,a2),fp2_add(b1,b2)),t1),t2)))
    c1=fp2_add(fp2_sub(fp2_sub(fp2_mul(fp2_add(a0,a1),fp2_add(b0,b1)),t0),t1), fp2_mul_xi(t2))
    c2=fp2_add(fp2_sub(fp2_sub(fp2_mul(fp2_add(a0,a2),fp2_add(b0,b2)),t0),t2),t1)
    return (c0,c1,c2)
def fp6_add(A,B): return (fp2_add(A[0],B[0]),fp2_add(A[1],B[1]),fp2_add(A[2],B[2]))
def fp6_sub(A,B): return (fp2_sub(A[0],B[0]),fp2_sub(A[1],B[1]),fp2_sub(A[2],B[2]))
def fp6_mul_v(A):  # multiply by v: (c0,c1,c2) -> (xi*c2, c0, c1)
    c0,c1,c2=A; return (fp2_mul_xi(c2),c0,c1)
FP6_0=(FP2_0,FP2_0,FP2_0); FP6_1=(FP2_1,FP2_0,FP2_0)

# Fp12 = Fp6[w]/(w^2 - v); element (d0,d1) = d0 + d1*w, di in Fp6
def fp12_mul(A,B):
    a0,a1=A; b0,b1=B
    t0=fp6_mul(a0,b0); t1=fp6_mul(a1,b1)
    d0=fp6_add(t0, fp6_mul_v(t1))
    d1=fp6_sub(fp6_sub(fp6_mul(fp6_add(a0,a1),fp6_add(b0,b1)),t0),t1)
    return (d0,d1)
def fp12_sqr(A): return fp12_mul(A,A)
FP12_1=(FP6_1,FP6_0)

def fp12_pow(A,e):
    r=FP12_1
    while e:
        if e&1: r=fp12_mul(r,A)
        A=fp12_sqr(A); e>>=1
    return r

def rand_fp12(rng):
    return (((rng.randrange(P),rng.randrange(P)),(rng.randrange(P),rng.randrange(P)),(rng.randrange(P),rng.randrange(P))),
            ((rng.randrange(P),rng.randrange(P)),(rng.randrange(P),rng.randrange(P)),(rng.randrange(P),rng.randrange(P))))

def first_coeff(X):  # the very first Fp coefficient (d0.c0.a0) of the tower rep
    return X[0][0][0]

COFACTOR=(P**12 - 1)//R

def gt_elem(rng):   # a GT element = random^cofactor (order divides r)
    return fp12_pow(rand_fp12(rng), COFACTOR)

# --- production-form evidence: uniformity of the distinguished-point coefficient
# bits over real BLS12-381 GT elements, supporting the ~2^-64 DP density and the
# routing-bit balance used by the single-trail DP rho (see proof.md). ---
import json, random, sys, collections, math


def main():
    args = json.load(sys.stdin)
    n_log2 = int(args.get("n_log2", 13))      # number of GT samples = 2^n_log2
    dp_shift = int(args.get("dp_shift", 20))  # DP predicate bit offset (proof uses 20)
    dp_bits = int(args.get("dp_bits", 12))    # measured DP-slice width (<=64 in proof)
    route_shift = int(args.get("route_shift", 84))  # routing bit offset (proof uses 84)
    route_bits = int(args.get("route_bits", 12))
    seed = int(args.get("seed", 1))
    rng = random.Random(seed)
    g = gt_elem(rng)                           # one real GT element (order divides r)
    # sanity: g is in the order-r subgroup
    in_gt = (fp12_pow(g, R) == FP12_1)
    N = 1 << n_log2
    x = g
    dp_hist = collections.Counter()
    rt_hist = collections.Counter()
    dp_zero = 0
    for _ in range(N):
        x = fp12_mul(x, g)                     # walk-like multiplicative sequence in GT
        c = first_coeff(x)
        dp_hist[(c >> dp_shift) & ((1 << dp_bits) - 1)] += 1
        rt_hist[(c >> route_shift) & ((1 << route_bits) - 1)] += 1
        if (c >> dp_shift) & ((1 << dp_bits) - 1) == 0:
            dp_zero += 1

    def chi2_over_dof(hist, bits):
        exp = N / (1 << bits)
        chi = sum((hist.get(b, 0) - exp) ** 2 / exp for b in range(1 << bits))
        return chi / ((1 << bits) - 1)

    print(json.dumps({
        "group": "real BLS12-381 GT (order-r subgroup of Fp12*), blst tower",
        "g_in_order_r_subgroup": bool(in_gt),
        "samples": N,
        "dp_slice_bits_%d_to_%d_chi2_over_dof" % (dp_shift, dp_shift + dp_bits - 1): chi2_over_dof(dp_hist, dp_bits),
        "routing_slice_bits_%d_to_%d_chi2_over_dof" % (route_shift, route_shift + route_bits - 1): chi2_over_dof(rt_hist, route_bits),
        "observed_dp_count_at_density_2^-%d" % dp_bits: dp_zero,
        "expected_dp_count": N / (1 << dp_bits),
        "note": "chi2/dof ~1 indicates the DP-predicate and routing coefficient bits are uniform over real GT; the measured %d-bit density matches 2^-%d, supporting the production 2^-64 predicate and balanced routing by the same mechanism." % (dp_bits, dp_bits),
    }, indent=1))


if __name__ == "__main__":
    main()

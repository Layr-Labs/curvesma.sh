"""Matrix generator (proof.md section 4): order basis by the recursive PM-basis (divide and conquer on the order)
and by the iterative M-basis; live toy checks that both return the same basis, that F P = 0 mod z^d, and that
the recursion's polynomial work stays within the charged mass 4 k^3 d log2(d) (k = m + n): the sum, over all
scalar polynomial products of the recursion, of the product length (one FFT product of length l costs fft_F(18 l)).

Polynomials: lists of k x k (or s x k) coefficient matrices over F_p, index = power of z.
"""
import math, random


def _mm(A, B, p, cnt):
    n, m, l = len(A), len(B), len(B[0])
    cnt[0] += n * m * l
    return [[sum(A[i][t] * B[t][j] for t in range(m)) % p for j in range(l)] for i in range(n)]


def _pmul(F, G, d, p, cnt):
    """F * G mod z^d (schoolbook in z; the production uses FFT products of the same entries)."""
    rows, cols = len(F[0]), len(G[0][0])
    cnt[1] += rows * len(G[0]) * cols * min(d, len(F) + len(G) - 1)
    out = [[[0] * cols for _ in range(rows)] for _ in range(min(d, len(F) + len(G) - 1))]
    for i, Fi in enumerate(F):
        for j, Gj in enumerate(G):
            if i + j < len(out):
                P = _mm(Fi, Gj, p, cnt)
                for r in range(rows):
                    for c in range(cols):
                        out[i + j][r][c] = (out[i + j][r][c] + P[r][c]) % p
    return out


def _step(D, shift, p):
    """One order step on the constant residual D (s x k): returns k x k basis of degree <= 1 and new shift."""
    s, k = len(D), len(D[0])
    order = sorted(range(k), key=lambda j: (shift[j], j))
    C = [[D[r][c] for c in range(k)] for r in range(s)]
    T0 = [[int(i == j) for j in range(k)] for i in range(k)]
    piv = []
    for r in range(s):
        jp = next((j for j in order if j not in piv and C[r][j] % p), None)
        if jp is None:
            continue
        inv = pow(C[r][jp], p - 2, p)
        for j in order:
            if j != jp and j not in piv and C[r][j] % p:
                f = C[r][j] * inv % p
                for rr in range(s):
                    C[rr][j] = (C[rr][j] - f * C[rr][jp]) % p
                for i in range(k):
                    T0[i][j] = (T0[i][j] - f * T0[i][jp]) % p
        piv.append(jp)
    T1 = [[0] * k for _ in range(k)]
    for jp in piv:
        for i in range(k):
            T1[i][jp] = T0[i][jp]; T0[i][jp] = 0
    sh = [shift[j] + (j in piv) for j in range(k)]
    return [T0, T1], sh


def mbasis(F, d, shift, p, cnt):
    k = len(F[0][0]); P = [[[int(i == j) for j in range(k)] for i in range(k)]]
    for t in range(d):
        R = _pmul(F, P, t + 1, p, cnt)
        D = R[t] if t < len(R) else [[0] * k for _ in range(len(F[0]))]
        S, shift = _step(D, shift, p)
        P = _pmul(P, S, len(P) + 1, p, cnt)
    return P, shift


def pmbasis(F, d, shift, p, cnt):
    if d == 1:
        return mbasis(F, 1, shift, p, cnt)
    h = d // 2
    P1, s1 = pmbasis(F, h, shift, p, cnt)
    R = _pmul(F, P1, d, p, cnt)
    F2 = R[h:d] + [[[0] * len(R[0][0]) for _ in R[0]] for _ in range(d - h - len(R[h:d]))]
    P2, s2 = pmbasis(F2, d - h, s1, p, cnt)
    return _pmul(P1, P2, len(P1) + len(P2) - 1, p, cnt), s2


def toy(seed=3, s=2, k=4, d=16, p=10007):
    rng = random.Random(seed)
    F = [[[rng.randrange(p) for _ in range(k)] for _ in range(s)] for _ in range(d)]
    cnt_i = [0, 0]; cnt_r = [0, 0]
    Pi, si = mbasis(F, d, [0] * k, p, cnt_i)
    Pr, sr = pmbasis(F, d, [0] * k, p, cnt_r)
    strip = lambda P: P[:max(i for i, c in enumerate(P) if any(any(x for x in row) for row in c)) + 1]
    same = strip(Pi) == strip(Pr) and si == sr
    R = _pmul(F, Pr, d, p, [0, 0])
    zero = all(not any(any(x for x in row) for row in R[t]) for t in range(min(d, len(R))))
    sum_deg = sum(sr)
    charged = 4 * k ** 3 * d * math.log2(d)
    return dict(s=s, k=k, d=d, same_as_iterative=same, order_ok=zero, sum_shift=sum_deg, sum_shift_expected=s * d,
                poly_mass_recursive=cnt_r[1], charged_mass=round(charged), within=cnt_r[1] <= charged)

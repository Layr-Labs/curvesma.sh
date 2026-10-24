"""Finite checks for the indexed-baby Cheon argument. Not a production attack."""
import json
import sys


R = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001
P = 0x1A0111EA397FE69A4B1BA7B6434BACD764774B84F38512BF6730D2A0F6B0F6241EABFFFEB153FFFFB9FEFFFFFFFFAAAB
FACTORS = (
    (2, 32),
    (3, 1),
    (11, 1),
    (19, 1),
    (10177, 1),
    (125527, 1),
    (859267, 1),
    (906349, 2),
    (2508409, 1),
    (2529403, 1),
    (52437899, 1),
    (254760293, 2),
)


def is_prime(n):
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return n == d
        d += 1 if d == 2 else 2
    return True


def check_arithmetic():
    n = R - 1
    prod = 1
    primes = []
    for q, e in FACTORS:
        assert is_prime(q)
        prod *= q ** e
        primes.append(q)
    assert prod == n
    for q in primes:
        assert pow(7, n // q, R) != 1

    d = 1 << 20
    assert n % d == 0
    m = n // d
    # Integer square root by Newton; the proof's N is this ceil.
    s = 1 << ((m.bit_length() + 1) // 2)
    while True:
        nxt = (s + m // s) // 2
        if nxt >= s:
            break
        s = nxt
    N = s if s * s >= m else s + 1
    while (N - 1) * (N - 1) >= m:
        N -= 1
    while N * N < m:
        N += 1
    assert (N - 1) * (N - 1) < m <= N * N
    assert N == 223621884914696853578958571250384896
    assert m == 50006747412801924209068050869165387952509453297164571592906626415194112
    delta = N * N - m
    assert delta == 218380746987008646073201729736704
    assert delta < N
    assert (1 << 117) < N < (1 << 118)
    assert N < 11 * (1 << 114)
    G = (m + N - 1) // N
    assert G == N
    # Every canonical lambda = j + N*l with 0 <= j < N lands in range, and the
    # giant exponents do not wrap.
    assert N * (N - 1) < m

    sched = 16 * 512 + 24 * 128 + 2048
    dbl = 8 * 512 + 20 * 128 + 1024
    norm = 12 * 512 + 1024
    assert sched == 13312 and dbl == 7680 and norm == 7168
    coeffs = (
        2 * 3 * sched,  # query
        2 * norm,  # normalisation of babies and giants
        2 * 2048,  # Fr geometric updates
        2 * 3 * 512,  # digit extraction
        256,  # record write
        118 * 256,  # sort comparisons
        118 * 256,  # sort moves
        120 * 256,  # binary search, G <= N
    )
    assert sum(coeffs) == 192768
    assert 192768 * 11 == 2120448
    assert 2120448 * (1 << 114) == 2070 * (1 << 124) + 768 * (1 << 114)
    assert 768 * (1 << 114) == 3 * (1 << 122)
    others = (
        6 * sched * (1 << 85)
        + 512 * dbl
        + 6 * norm * (1 << 85)
        + 8 * (611 * 512 + 4096)
        + 4 * (611 * 512)
        + 4096 * 256 * (dbl + sched)
        + 64 * 512 * 2048
        + (1 << 20) * 4096
        + (1 << 48)
    )
    assert others < (1 << 104)
    # Saturated linear term plus the other-term cap.
    assert 192768 * 11 * (1 << 114) + (1 << 104) == 2070 * (1 << 124) + 3 * (1 << 122) + (1 << 104)
    assert 3 * (1 << 122) + (1 << 104) < (1 << 124)
    assert 2071 * 100 < 107 * 2000
    assert 107 ** 10 < 2 * (100 ** 10)
    # Memory: 342 * 2^118 < 384 * 2^118 = (3/2) * 2^126, and (3/2)^5 = 243/32 < 8.
    assert 342 < 384
    assert 243 < 256
    assert P.bit_length() == 381

    # Straight-line envelopes dominate the derived programs in proof.md.
    derived_sched = 15 * 383 + 14 * 128 + 128 + 512
    derived_dbl = 7 * 383 + 14 * 128 + 256
    assert derived_sched == 8177 and sched > derived_sched
    assert derived_dbl == 4729 and dbl > derived_dbl


def toy_indexed_cheon():
    order = 233
    d = 8
    assert (order - 1) % d == 0
    m = (order - 1) // d
    primes = []
    x = order - 1
    pdiv = 2
    while pdiv * pdiv <= x:
        if x % pdiv == 0:
            primes.append(pdiv)
            while x % pdiv == 0:
                x //= pdiv
        pdiv += 1 if pdiv == 2 else 2
    if x > 1:
        primes.append(x)
    g = None
    for cand in range(2, order):
        if all(pow(cand, (order - 1) // q, order) != 1 for q in primes):
            g = cand
            break
    assert g is not None
    omega = pow(g, d, order)
    assert pow(omega, m, order) == 1
    kappa = pow(g, m, order)
    assert pow(kappa, d, order) == 1
    N = 1
    while N * N < m:
        N += 1
    P0 = 5
    rank_differs = False
    for secret in range(1, order):
        Ad = (pow(secret, d, order) * P0) % order
        omega_inv = pow(omega, m - 1, order)
        records = []
        scal = 1
        for j in range(N):
            point = (scal * Ad) % order
            records.append((point, j))
            scal = (scal * omega_inv) % order
        records.sort()
        for rank, (_point, j) in enumerate(records):
            if rank != j:
                rank_differs = True
        gamma = pow(omega, N, order)
        found = None
        gscal = 1
        for lam1 in range(N):
            giant = (gscal * P0) % order
            # Binary search on the point key; the returned value is the stored index.
            lo, hi = 0, len(records)
            while lo < hi:
                mid = (lo + hi) // 2
                if records[mid][0] < giant:
                    lo = mid + 1
                else:
                    hi = mid
            if lo < len(records) and records[lo][0] == giant:
                found = (records[lo][1], lam1)
                break
            gscal = (gscal * gamma) % order
        assert found is not None
        lam0, lam1 = found
        lam = (lam0 + N * lam1) % m
        assert pow(omega, lam, order) == pow(secret, d, order)
        sigma = pow(g, lam, order)
        U = (sigma * P0) % order
        # Order-d coset search, again storing j0 in the record.
        B = 1
        while B * B < d:
            B += 1
        kappa_inv = pow(kappa, d - 1, order)
        babies = []
        ks = 1
        for j0 in range(B):
            babies.append(((ks * U) % order, j0))
            ks = (ks * kappa) % order
        babies.sort()
        step = pow(kappa_inv, B, order)
        hs = 1
        j = None
        for j1 in range(B):
            giant = (hs * (secret * P0)) % order
            lo, hi = 0, len(babies)
            while lo < hi:
                mid = (lo + hi) // 2
                if babies[mid][0] < giant:
                    lo = mid + 1
                else:
                    hi = mid
            if lo < len(babies) and babies[lo][0] == giant:
                j = babies[lo][1] + B * j1
                break
            hs = (hs * step) % order
        assert j is not None and 0 <= j < d
        recovered = (sigma * pow(kappa, j, order)) % order
        assert recovered == secret
        assert (recovered * P0) % order == (secret * P0) % order
    assert rank_differs


def main():
    payload = json.load(sys.stdin)
    assert payload["mode"] == "all"
    check_arithmetic()
    toy_indexed_cheon()
    print(json.dumps({
        "arithmetic": True,
        "toy_order": 233,
        "toy_recovered_all": True,
        "stored_index_used": True,
        "scope": "identities and a 233-element toy group; not a BLS12-381 execution",
    }))


if __name__ == "__main__":
    main()

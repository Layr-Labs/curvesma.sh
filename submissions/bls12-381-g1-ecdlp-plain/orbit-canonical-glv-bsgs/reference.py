"""Orbit-canonical GLV baby-step/giant-step on BLS12-381 G1.

Baby steps are a*P0 for a in [0, m). Giant steps subtract b*phi(P0), and
phi(P0) = [lambda]P0 with lambda^2 + lambda + 1 = r as integers. Each stored
key is the canonical compressed encoding of an order-6 orbit, so each chain
point contributes one record. The scored claim is the ledger in proof.md.
"""

from __future__ import annotations

BLS12_381_R = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001
# lambda^2 + lambda + 1 = r, and [lambda] is the eigenvalue of phi(x, y) = (beta*x, y).
LAMBDA = 0xAC45A4010001A40200000000FFFFFFFF
BETA = 0x1A0111EA397FE699EC02408663D4DE85AA0D857D89759AD4897D29650FB85F9B409427EB4F49FFFD8BFD00000000AAAC
BETA2 = 0x5F19672FDF76CE51BA69C6076A0F77EADDB3A93BE6F89688DE17D813620A00022E01FFFFFFFEFFFE

M_PROD = 82563089651286698697488428308387906491
H_PROD = 82563089651286698697488428308387906490

KEY_BYTES = 48
PER_POINT = 1750  # two field products, orbit-min selection, store, 12-pass sort, merge
PER_IDENTITY = 920
VERIFY_OPS = 2048
BUCKET_INSTR = 12 * (1 << 32) * 6
FIXED_INSTR = 8192


def _ledger_numerator(m: int, h: int) -> int:
    """Numerator of W = H + R/2000. Denominator is 2000."""
    group_ops = (m - 1) + (h - 1) + VERIFY_OPS
    instructions = (
        PER_IDENTITY
        + PER_POINT * (m - 1)
        + PER_POINT * h
        + 2 * BUCKET_INSTR
        + FIXED_INSTR
    )
    return 2000 * group_ops + instructions


def coverage(m: int, h: int) -> int:
    """Cardinality of the explicit four-multiplier set at these parameters."""
    return h * (2 * m - 1) + m * (m - 1)


def select_parameters(r: int = BLS12_381_R):
    """Production (m, h) with h = m-1 and |S| = (m-1)*(3*m-1) >= ceil(39r/100)."""
    if r != BLS12_381_R:
        raise ValueError("toy orders are constructed in experiment.py")
    m, h = M_PROD, H_PROD
    lam = LAMBDA
    assert lam * lam + lam + 1 == r
    assert h == m - 1
    assert 2 * m - 1 <= lam
    assert (m - 1) + (h - 1) <= lam - 1
    target = (39 * r + 99) // 100
    got = coverage(m, h)
    assert got == (m - 1) * (3 * m - 1)
    assert got >= target
    assert (m - 2) * (3 * (m - 1) - 1) < target
    assert m < (1 << 126) and h < (1 << 126)
    assert m + h < (1 << 128)
    return m, h


def orbit_min_scalar(a: int, mod: int, lam: int) -> int:
    """Smallest residue in the order-6 orbit of a. Toy stand-in for the canonical key."""
    if a % mod == 0:
        return 0
    sigmas = (1, -1, lam, -lam, 1 + lam, -(1 + lam))
    return min((s * a) % mod for s in sigmas)


def explicit_set(m: int, h: int, mod: int, lam: int) -> set[int]:
    """Residues b*lam + sigma*a for sigma in {±1, ±(1+lam)}, a=0 only with +1."""
    out = set()
    sigmas = (1, -1, 1 + lam, -(1 + lam))
    for b in range(h):
        for a in range(m):
            for sigma in sigmas:
                if a == 0 and sigma != 1:
                    continue
                out.add((b * lam + sigma * a) % mod)
    return out


def solve_orbit(target: int, m: int, h: int, mod: int, lam: int):
    """Return target if its orbit meets a baby orbit under stride lam, else None."""
    babies = {}
    for a in range(m):
        key = orbit_min_scalar(a, mod, lam)
        babies.setdefault(key, []).append(a)
    sigmas = (1, -1, lam, -lam, 1 + lam, -(1 + lam))
    found = None
    for b in range(h):
        key = orbit_min_scalar((target - b * lam) % mod, mod, lam)
        for a in babies.get(key, ()):
            for sigma in sigmas:
                cand = (b * lam + sigma * a) % mod
                if cand == target:
                    found = cand
    return found


if __name__ == "__main__":
    m, h = select_parameters()
    num = _ledger_numerator(m, h)
    assert num == 619223172384650240231163212313527778068736
    print("m", m)
    print("h", h)
    print("coverage", coverage(m, h))
    print("ledger_numerator", num)

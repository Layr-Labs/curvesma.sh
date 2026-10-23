"""Toy checks for the Frobenius-scaled collision equation.

The BLS12-381 constants are the demontgomery values of blst coeffs2[] at
commit f7a88067ed98385ebc6cf41ccecf21ead706a5b7. The group check is a small
cyclic group, not GT.
"""

import json
import sys


P = 0x1A0111EA397FE69A4B1BA7B6434BACD764774B84F38512BF6730D2A0F6B0F6241EABFFFEB153FFFFB9FEFFFFFFFFAAAB

# demontgomery(coeffs2[n-1]) for n = 1, 2, 3
C1 = 0x1A0111EA397FE699EC02408663D4DE85AA0D857D89759AD4897D29650FB85F9B409427EB4F49FFFD8BFD00000000AAAD
C2 = 0x1A0111EA397FE699EC02408663D4DE85AA0D857D89759AD4897D29650FB85F9B409427EB4F49FFFD8BFD00000000AAAC
C3 = 0x1A0111EA397FE69A4B1BA7B6434BACD764774B84F38512BF6730D2A0F6B0F6241EABFFFEB153FFFFB9FEFFFFFFFFAAAA


def check_constants():
    if pow(C1, 6, P) != 1:
        raise SystemExit("c1^6 != 1")
    if pow(C1, 3, P) == 1:
        raise SystemExit("c1 has order dividing 3")
    if pow(C1, 2, P) != C2:
        raise SystemExit("c2 != c1^2")
    if pow(C1, 3, P) != C3:
        raise SystemExit("c3 != c1^3")


def check_collision(trials):
    # Prime field of scalars. 2 is a generator mod 13 and has order 12.
    q = 13
    zeta = 2
    if pow(zeta, 12, q) != 1 or pow(zeta, 6, q) == 1:
        raise SystemExit("toy automorphism is not order 12")
    for i in range(trials):
        y = (2 + i) % q
        alpha1 = (5 + 3 * i) % q
        beta1 = (7 + i) % q
        delta = i % 12
        scale = pow(zeta, delta, q)
        # Force a nonzero denominator, then the matching alpha2.
        beta2 = (beta1 * scale + 1 + (i % 5)) % q
        alpha2 = (alpha1 * scale + y * (beta1 * scale - beta2)) % q
        numer = (alpha1 * scale - alpha2) % q
        denom = (beta2 - beta1 * scale) % q
        if denom == 0:
            raise SystemExit("degenerate toy denominator")
        y_hat = (numer * pow(denom, -1, q)) % q
        if y_hat != y:
            raise SystemExit("collision equation failed")


def main():
    raw = sys.stdin.read() if sys.stdin.isatty() is False else ""
    trials = 20
    if raw.strip():
        data = json.loads(raw)
        trials = int(data["trials"])
    check_constants()
    check_collision(trials)
    json.dump({"ok": True, "trials": trials}, sys.stdout)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()

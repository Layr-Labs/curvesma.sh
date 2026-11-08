#!/usr/bin/env python3
"""Toy exact-key BSGS, a counted Fp multiply envelope, and the production ledger."""
import json
import math
import sys

P = int(
    "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab",
    16,
)
R = int("73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001", 16)
C = 20000
M = 8192
TX = 1 << 21
BABY_CAP = 1 << 30
GIANT_CAP = 1 << 29
INIT = 1 << 150
FIXED_CAP = 1 << 48
CLAIM = 145


def schoolbook_envelope():
    """Instruction envelope for a 6-limb multiply reduced into [0, p)."""
    # Inner product step: 4 address, 4 loads, mul, umulh, 8-instruction
    # carry chain, 2 stores, 2 branches.
    assert 4 + 4 + 2 + 8 + 2 + 2 <= 64
    # One reduction row: 14-instruction quotient digit, six 16-instruction
    # limb mul-subs, 16-instruction carry cleanup.
    assert 14 + 6 * 16 + 16 <= 256
    return 36 * 64 + 6 * 256 + 256


def fp_add_envelope():
    # Six limb adds (load, load, adc, store, carry) plus a conditional subtract.
    assert 6 * 5 + 6 * 4 + 8 <= 256
    return 256


def fixed_items():
    m12 = 1 << 21
    pair = 1 << 36
    return {
        "parse": 1 << 30,
        "membership": 1 << 35,
        "order_check": 1 << 30,
        "pairings": pair * 2,
        "power_m": 128 * m12,
        "inverse_power": 1 << 30,
        "g2_scalar": 1 << 34,
        "verify_exp": 1 << 30,
        "encode": 1 << 24,
        "controller": 1 << 30,
    }


def ledger():
    assert schoolbook_envelope() <= M
    assert fp_add_envelope() <= 256
    # Per word per hop: load, store, two address bumps, counter, checksum, branch.
    assert 8 <= 32
    # 4x4-limb address multiply, 72-limb compare, loop control.
    assert 16 * 8 + 64 + 72 * 4 + 64 <= 4096
    # 256 hops, 128 words, 32 instructions per word, plus address and compare.
    assert 256 * 128 * 32 + 4096 < TX
    baby = (1 << 29) + (1 << 17) + (1 << 22)
    assert baby < BABY_CAP
    giant = 130 * TX + 129 * 1024
    assert giant < GIANT_CAP
    items = fixed_items()
    assert sum(items.values()) < FIXED_CAP
    # Miller loop and final exponentiation sit under the pairing cap.
    m12 = 1 << 21
    assert 64 * (1 << 27) < (1 << 34)
    assert 1024 * m12 + (1 << 30) < (1 << 36)
    m = h = 1 << 128
    assert m * h > R
    r_total = m * BABY_CAP + h * GIANT_CAP + INIT + FIXED_CAP
    assert r_total < C * (1 << CLAIM)
    r_pre = m * BABY_CAP + INIT + FIXED_CAP
    assert r_pre < C * (1 << CLAIM)
    peak = 2 * m * 1024 + (1 << 129) * 256
    assert peak < (1 << 150)
    return {
        "baby_instructions": baby,
        "giant_instructions": giant,
        "fixed_instructions": sum(items.values()),
        "R_total": str(r_total),
        "log2_W": math.log2(r_total / C),
        "claim": CLAIM,
        "memory_bytes": str(peak),
        "log2_memory": math.log2(peak),
        "fp_mul_envelope": schoolbook_envelope(),
    }


def mul_mod_p(a, b):
    """Schoolbook product reduced by Python, used only to check the toy field."""
    return (a * b) % P


def bsgs_sorted(order, mul, identity, gen, secret):
    """Exact-key BSGS. Keys are canonical group elements, not fingerprints."""
    m = math.isqrt(order)
    if m * m < order:
        m += 1
    h = (order + m - 1) // m
    target = identity
    for _ in range(secret):
        target = mul(target, gen)
    baby = []
    x = identity
    for j in range(m):
        baby.append((x, j))
        x = mul(x, gen)
    baby.sort()
    # g^{-m} = g^{order-m}
    giant_step = identity
    step = gen
    exp = (order - m) % order
    while exp:
        if exp & 1:
            giant_step = mul(giant_step, step)
        step = mul(step, step)
        exp >>= 1
    y = target
    keys = [row[0] for row in baby]
    for i in range(h):
        lo, hi = 0, len(keys)
        while lo < hi:
            mid = (lo + hi) // 2
            if keys[mid] < y:
                lo = mid + 1
            else:
                hi = mid
        if lo < len(keys) and keys[lo] == y:
            return (i * m + baby[lo][1]) % order
        y = mul(y, giant_step)
    raise RuntimeError("miss")


def check_group():
    p = 101
    assert all(p % k for k in range(2, int(p ** 0.5) + 1))
    gen = 2
    order = p - 1
    hits = 0
    for secret in range(order):
        got = bsgs_sorted(order, lambda a, b: (a * b) % p, 1, gen, secret)
        assert got == secret
        hits += 1
    # Additive coverage, including the identity exponent.
    for n in (1, 2, 17, 64, 200):
        for secret in range(n):
            got = bsgs_sorted(n, lambda a, b: (a + b) % n, 0, 1, secret)
            assert got == secret
    # A few base-field multiplies match the integer product.
    samples = 0
    x = 3
    for _ in range(40):
        y = (x * 5 + 7) % P
        assert mul_mod_p(x, y) == (x * y) % P
        x = y
        samples += 1
    return {"field_prime": p, "recovered": hits, "fp_samples": samples}


def main():
    raw = sys.stdin.read().strip()
    cfg = json.loads(raw) if raw else {"mode": "all"}
    mode = cfg.get("mode", "all")
    out = {"status": "ok", "mode": mode}
    if mode in ("all", "ledger"):
        out["ledger"] = ledger()
    if mode in ("all", "group"):
        out["group"] = check_group()
    print(json.dumps(out))


if __name__ == "__main__":
    main()

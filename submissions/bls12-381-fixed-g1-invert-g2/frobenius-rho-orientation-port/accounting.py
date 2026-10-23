#!/usr/bin/env python3
"""Conservative fixed-G1 orientation ledger."""

from math import log2


# The public fixed-G2 analysis gives this conservative S0 exponent for the
# same order-12 GT quotient.  R=2048 can only lower it, so we retain it.
S0 = 2.0**125.6522
ARRIVALS = 1.05 * S0
STEP = 1.16
CYCLE_AND_ROUTING_FACTOR = 1.001
FIXED_AND_FINISHING = 2.0**100

TOTAL = ARRIVALS * STEP * CYCLE_AND_ROUTING_FACTOR + FIXED_AND_FINISHING
TOTAL_LOG2 = log2(TOTAL)
CLAIM = 126.0

assert TOTAL_LOG2 < CLAIM
print(f"arrival_log2={log2(ARRIVALS):.9f}")
print(f"total_log2={TOTAL_LOG2:.9f}")
print(f"claim={CLAIM:.2f}")

#!/usr/bin/env python3
"""Conservative ledger transformation for the 22-bit comb candidate."""

from math import log2


# Published upper components of Randomblock1's accepted 20-bit package.  One
# ulp at six decimals is added before exponentiation so this script never
# treats the printed decimals as exact lower values.
OLD_H = 2.0**120.460293
OLD_OTHER = 2.0**118.987424

# The old hot path has twelve additions per evaluated point.  The 22-bit
# recoding has twelve digits and therefore eleven additions.  OLD_H may also
# contain fixed work, so add the whole new fixed-work reserve rather than
# assuming every old reference operation scales.
NEW_H = OLD_H * (11.0 / 12.0) + 2.0**28

# Keep every published non-reference charge unchanged.  The extra reserve is
# deliberately far larger than the enlarged comb construction and addressing
# work, so no reduction in per-point ancillary work is claimed.
NEW_OTHER = OLD_OTHER + 2.0**112
NEW_TOTAL = NEW_H + NEW_OTHER
NEW_LOG2 = log2(NEW_TOTAL)

CLAIM = 120.83
assert NEW_LOG2 < CLAIM

print(f"new_H_log2={log2(NEW_H):.9f}")
print(f"new_other_log2={log2(NEW_OTHER):.9f}")
print(f"new_total_log2={NEW_LOG2:.9f}")
print(f"claim={CLAIM:.2f}")

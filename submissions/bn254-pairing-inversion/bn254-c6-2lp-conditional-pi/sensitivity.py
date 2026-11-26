"""Arithmetic-only sensitivity calculation; does not establish H-LP."""

import json
import math


ROWS = {
    "collection": 83.773593,
    "chi_setup": 65.712288,
    "la": 77.500018,
    "gen": 73.170471,
    "lists": 78.536663,
    "two_logs": 78.941351,
    "canon": 75.695321,
    "sm": 74.465934,
    "misc": 71.799087,
    "fb": 65.526296,
    "split": 71.808841,
    "gtree": 59.508896,
    "cert": 50.000001,
    "verify": 11.000001,
    "wrapper": 53.712288,
    "output": 11.761903,
}
GATE_K = 7194244577235879602
FACTOR = 47821 / 16874


def logsum2(values):
    top = max(values)
    return top + math.log2(sum(2 ** (v - top) for v in values))


def main():
    changed = dict(ROWS)
    changed["collection"] -= math.log2(FACTOR)
    opened = logsum2(changed.values())
    score = math.log2(1 + (GATE_K / 2**64) * 2**opened)
    assert score <= 81.210350
    print(json.dumps({
        "scope": "arithmetic sensitivity only; H-LP untested",
        "c4_row_ratio": FACTOR,
        "hypothetical_collection_log2": changed["collection"],
        "hypothetical_opened_run_log2": opened,
        "hypothetical_expected_work_log2": score,
        "declared_upper_bound_log2": 81.210350,
    }, sort_keys=True))


if __name__ == "__main__":
    main()

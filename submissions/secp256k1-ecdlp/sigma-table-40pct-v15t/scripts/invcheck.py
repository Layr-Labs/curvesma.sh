"""invcheck.py -- the endpoint comparison, DERIVED from invariant-row.json.

No caller writes the formula: the residual is evaluated from the row's term list,
so a change of convention is a one-line change of the JSON and every consumer
(harness, census, certificate check) moves with it.
"""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
ROW_PATH = os.path.join(HERE, "invariant-row.json")


def load(path=ROW_PATH):
    return json.load(open(path))


def coef_value(c, a, cval, r):
    if c == "1":
        return 1 % r
    if c == "-1":
        return (-1) % r
    if c == "a":
        return a % r
    if c == "-a":
        return (-a) % r
    if c == "c":
        return cval % r
    if c == "-c":
        return (-cval) % r
    raise ValueError("unknown coefficient %r" % c)


def residual(cur, pts, a, cval, row=None):
    """pts maps the row's point names to curve points; returns the residual."""
    row = row or load()
    acc = None
    for t in row["terms"]:
        k = coef_value(t["coef"], a, cval, cur.r)
        if k == 0:
            continue
        acc = cur.add(acc, cur.mul(k, pts[t["point"]]))
    return acc


def holds(cur, pts, a, cval, row=None):
    return residual(cur, pts, a, cval, row) is None

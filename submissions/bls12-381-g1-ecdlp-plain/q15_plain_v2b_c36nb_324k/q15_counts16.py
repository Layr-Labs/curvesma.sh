"""q15_counts16.py — Exact expected special-q node counter for multi-retry NB descent (matching c36p counts16.py)."""


def expected_nodes(Pm, dres):
    bins = dres["bins"]
    DET = dres["DET"]
    C = {}
    def cmax(side, lim):
        return max([C[b, side] for b in bins if (b, side) in C and b[0] < lim - 1e-9] + [0.0])
    for b in bins:
        for side in ("f", "g"):
            d = DET[b, side]
            sub = (d["kA"] * (cmax("f", d["nA"]) + cmax("g", d["nA"])) if d["kA"] else 0.0) + (
                d["kB"] * (cmax("f", float(d["n"])) + cmax("g", float(d["n"]))) if d["kB"] else 0.0
            )
            C[b, side] = 1.0 + d["retry"] * sub
    bt = dres["boot"]
    sub_b = bt["kA"] * cmax("g", bt["nA"]) + cmax("g", float(bt["n1"]))
    return 2.0 * bt["retry"] * sub_b, C

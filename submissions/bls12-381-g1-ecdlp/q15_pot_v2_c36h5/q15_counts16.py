"""q15_counts16.py — Expected descent node count across both individual logs."""

def expected_nodes(Pm, res):
    return float(res["boot"]["e_nodes_2log"]), dict(res["EN_NODE"])

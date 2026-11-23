"""q15_counts16.py — Exact tree-slot counter for Pure k=1 Zero-Retry Boole + Wald i.i.d. Boot."""


def expected_nodes(Pm, dres):
    bt = dres["boot"]
    return 2.0 * bt["e_runs"] * bt["n_slots_tree"], dres["NSL"]

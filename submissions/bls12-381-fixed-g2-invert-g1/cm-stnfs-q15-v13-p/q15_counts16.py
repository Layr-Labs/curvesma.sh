"""q15_counts16.py — Expected node-trial visit upper bound for 2-State Wald-Boole + Wald i.i.d. Boot."""


def expected_nodes(Pm, dres):
    bt = dres["boot"]
    return 2.0 * bt["e_runs"] * bt["e_node_trials_1run"], dres["NSL"]

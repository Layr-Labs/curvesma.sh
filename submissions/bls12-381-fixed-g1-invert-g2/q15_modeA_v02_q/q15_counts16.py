"""q15_counts16.py — Exact deterministic slot counter for Mode-A pure k=1 Boole + Wald i.i.d. fresh-e boot."""


def expected_nodes(Pm, dres):
    bt = dres["boot"]
    return 2.0 * bt["e_boot_runs"] * bt["c_det_one"], dres["C_DET"]

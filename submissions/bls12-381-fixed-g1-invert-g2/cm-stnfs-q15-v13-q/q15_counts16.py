"""q15_counts16.py — Exact 2-state expected slot-trial counter for v13 Wald-Boole + Fresh-e Boot."""


def expected_nodes(Pm, dres):
    bt = dres["boot"]
    return 2.0 * bt["e_boot_runs"] * bt["en_ch_tree"], dres["EN_CH"]

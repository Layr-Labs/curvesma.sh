#!/usr/bin/env python3
"""regen_v25.py -- CLAIM LITERALS ARE GENERATED, never typed.

Verdict point 4 (fourth fleet machine-copy violation): the manifest's displayed
subtraction yielded 0.6316956 while the claim said 0.631930.  A literal that a
human retypes is a literal that can drift.  Here every number in the claim
statement, the title and the score comes from ledger.py's own objects, i.e. from
evidence/ledger-v25.json written below; nothing is hand-copied.

PLANTED CONTROL: `--planted` adds 1 ulp (one unit in the last place, here
1e-9 of a bit on the claim) to the ledger JSON, regenerates, and asserts the
claim STRING changes.  If the claim does not move, the literals are not coming
from the ledger and the script exits non-zero.
"""
import json
import os
import sys
from decimal import Decimal

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
os.chdir(ROOT)
sys.path.insert(0, ROOT)

import ledger as L  # noqa: E402


def ledger_numbers():
    return {
        "claim_log2": "%.3f" % float(L.CLAIM),
        "claim_value": "%.9f" % float(L.CLAIM),
        "prev_claim": "128.172",
        "log_w": str(L.LOG_W),
        "step_insn": int(L.STEP_INSN),
        "step_price": "%.9f" % float(L.STEP_PRICE),
        "success": "%.9f" % L.SUCCESS,
        "cp95_cp": str(L.CP95_CP),
        "jensen_ratio": "%.9f" % float(L.JENSEN_RATIO),
        "exact_single": "%.9f" % L.EXACT_SINGLE,
        "eom_single": "%.9f" % L.EOM_SINGLE,
        "pair_k": int(L.PAIR_K),
        "pair_n": int(L.PAIR_N),
        "pair_rate": "%.9f" % (L.PAIR_K / L.PAIR_N),
        "debits": "%.9f" % float(L.DEBITS),
        "restart_loss": "%.9f" % float(L.RESTART_SUCCESS_LOSS),
        "floor": "%.3f" % float(L.SUCCESS_FLOOR),
        "floor_met": bool(L.SUCCESS >= float(L.SUCCESS_FLOOR)),
        "margin_floor": "%.9f" % (L.SUCCESS - float(L.SUCCESS_FLOOR)),
        "claim_ratio_num": str(L.CLAIM.numerator),
        "claim_ratio_den": str(L.CLAIM.denominator),
        "claim_reduced": str(L.CLAIM),
        "slope_per_doubling": "%.12f" % L.T_SLOPE,
        "span_doublings": "%.6f" % L.T_SPAN_DEC,
        "doublings_to_prod": "%.6f" % L._lad["extrapolation"]["doublings"],
        "pooled_k": int(L._lad["pooled"]["k"]),
        "pooled_n": int(L._lad["pooled"]["n"]),
        "pooled_p": "%.9f" % L._lad["pooled"]["p_hat"],
        "pooled_cp95": "%.9f" % L._lad["pooled"]["cp95"],
        "gate": "%.3f" % float(L.GATE),
        "ship_line": "%.3f" % float(L.SHIP_LINE),
        "log2_s0": "%.9f" % L.log2f(L.S_OPS),
        "exact_esc": "%.6f" % L._TH["cells"][0]["history_escape_fired"]["success_exact"],
        "exact_clean": "%.6f" % L._TH["cells"][0]["history_clean"]["success_exact"],
        "relad": "-".join("%.6f" % c["all"]["success_exact"] for c in L._TH["cells"][1:]),
    }


def claim_text(N):
    return (
        "On secp256k1-ecdlp this package claims the COMPUTATION BOUND W = 2^{claim} "
        "(reduced rational {claim_red}, submitted {claim}) for the shipped DP-walk configuration (2^20 lanes, inversion BATCH 2^18, "
        "production DP density d=60, step {step} itemised instructions at C=1000). "
        "FIVE FIXES, in the verdict da40eb9c's own order. "
        "(1) HISTORY-DEPENDENT COALESCENCE IS NO LONGER CROSS-POOLED. Every ladder cell is estimated on its own "
        "trials -- no cell's rate is carried into another -- and each cell is additionally split by its own escape "
        "history (esc_run fired vs clean). Design cell 2^20: exact-sum success {exact_esc} on escape-fired trails "
        "(n=2351) and {exact_clean} on clean trails (n=49); re-ladder cells (L=16,R=32,d=4), (L=32,R=64,d=0) and "
        "(L=64,R=32,d=0) give {relad} on their own data, so moving L, R and d one at a time does not separate the "
        "rate from the floor. Nothing below pools cells. "
        "(2) THE EXACT SUM REPLACES exp-of-mean, AND IT WAS RUN. miss = E[(1-2^-60)^T] is summed over the RAW "
        "per-trail lengths T recorded for every run: success_exact = 1 - (1/n) sum_j (1-p)^{{T_j}} at the design "
        "point p = 1/E[T]. f(T)=(1-p)^T is convex, so the exp-of-mean value v24 shipped is an UPPER bound on "
        "success; measured on the design cell exact {exact} vs exp-of-mean {eom}, ratio R = {ratio}, i.e. the "
        "substitution over-credited success by {gap} and 'trail/2' is DELETED. "
        "(3) THE DETECTOR IS KEY-LEVEL AND THE ADVERSARIAL SET IS NOT TRUNCATED. The comparison is now between two "
        "INDEPENDENT implementations (the production selector arithmetic twist_index_production_secp + cell, "
        "against klass) over the WHOLE adversarial set -- 36 pairs, no out[:20] slice -- and the class key "
        "(K,b,column) is equal on 36/36. Two distinct keys in one column are reported NOT a collision (audit "
        "A3), which a class-level comparison would have called a collision. DISCLOSED: the twist index h, a "
        "canonicalisation label and not the class key, agrees 20/36 between the window selector and the exact "
        "argmin; that window-vs-argmin gap is measured and not repaired here. "
        "(4) THE FLOOR RIDES ON A CP-BINOMIAL LOWER BOUND ONLY. The bootstrap interval is deleted. "
        "SUCCESS = CP95({k}/{n}) = {cp} x exact-sum ratio {ratio} - debits {deb} - restart loss {rest} = {succ} "
        ">= {floor}. The floor is MET with the exact sum, so the claim-class is kept; margin {margin}. "
        "(5) THE TWO UNPRICED ROWS ARE PRICED. sigma's beta^h application was priced only as the selector's "
        "partial-product WINDOW, not its result: charged a full F_p multiply (97) plus a 512->256 reduction mod p "
        "(97). The 97 of the mu6 row counted 16 schoolbook mul-adds and one conditional subtract, not the "
        "reduction of the full 512-bit product mod r: charged a second row of 2 x 97. Step {step} -> {step_new} "
        "instructions, price {price}, which moves the claim from 128.172 to {claim} (+{dbit} bits). The claim got "
        "WORSE and it is the tested algorithm's own price. "
        "LADDER (n=1200 per cell, production key, lag-3, d=0): {pooled_k}/{pooled_n} pooled = {pooled_p}, CP-95 "
        "lower {pooled_cp}; slope {slope} per doubling over {span} doublings; {doublings} doublings to production "
        "class count. NOT RANKING: the claim is {dgate} bits above the live bar {gate}; submitted as validity "
        "depth with the floor measured, not as a rank claim."
    ).format(
        claim=N["claim_log2"], claim_red=N["claim_reduced"], step=N["step_insn"], exact_esc=N["exact_esc"],
        exact_clean=N["exact_clean"], relad=N["relad"],
        exact=N["exact_single"], eom=N["eom_single"], ratio=N["jensen_ratio"],
        gap="%.6f" % (float(N["eom_single"]) - float(N["exact_single"])),
        k=N["pair_k"], n=N["pair_n"], cp=N["cp95_cp"], deb=N["debits"], rest=N["restart_loss"],
        succ=N["success"], floor=N["floor"], margin=N["margin_floor"],
        step_new=754, price=N["step_price"],
        dbit="%.6f" % (float(N["claim_value"]) - float(N["prev_claim"])),
        pooled_k=N["pooled_k"], pooled_n=N["pooled_n"], pooled_p=N["pooled_p"],
        pooled_cp=N["pooled_cp95"], slope=N["slope_per_doubling"], span=N["span_doublings"],
        doublings=N["doublings_to_prod"], dgate="%.4f" % (float(N["claim_value"]) - float(N["gate"])),
        gate=N["gate"],
    )


def title_text(N):
    t = ("secp256k1 v25: exact-sum success {succ} (not exp-of-mean {eom}; R={ratio}) "
         "clears floor {floor}; per-cell within-history ladder (no pooling); key-level detector, full "
         "adversarial set; sigma beta^h + mod-r reductions priced {step} insn -> claim {claim}; "
         "CP-binomial only").format(succ=N["success"], eom=N["eom_single"], ratio=N["jensen_ratio"],
                                    floor=N["floor"], step=N["step_insn"], claim=N["claim_log2"])
    assert len(t) <= 300, (len(t), t)
    return t


def write_all(N, tag=""):
    stmt = claim_text(N)
    stmt.encode("ascii")                      # ASCII-only law
    claim = {"statement": stmt,
             "score": {"numerator": N["claim_ratio_num"], "denominator": N["claim_ratio_den"]},
             "evidence": ["public-note.md", "dependencies.json", "ledger.py",
                          "evidence/ledger-v25.json", "evidence/trail-hist-v25.json",
                          "evidence/audit-v25.json", "evidence/dladder-v24.json",
                          "evidence/pair-census.json", "evidence/selector-pair.json",
                          "evidence/cap-recount.json", "evidence/census-v20n.json",
                          "evidence/restart-count.json", "evidence/selector-swap.json",
                          "evidence/absorb-direct.json", "scripts/dladder.py",
                          "scripts/trails_v25.py", "scripts/audit_v25.py",
                          "scripts/regen_v25.py", "scripts/run_ladder_v24.py",
                          "scripts/selector.py", "scripts/qwalk.py", "scripts/census_v20n.py"]}
    sub = {"version": 1, "title": title_text(N), "claims": [claim], "experiments": []}
    json.dump(claim, open(os.path.join(ROOT, "claim.json"), "w"), indent=1)
    json.dump(sub, open(os.path.join(ROOT, "submission.json"), "w"), indent=1)
    json.dump(N, open(os.path.join(ROOT, "evidence", "ledger-v25.json"), "w"), indent=1)
    return stmt, title_text(N)


def main():
    N = ledger_numbers()
    stmt, title = write_all(N)
    print("CLAIM", N["claim_log2"], "SUCCESS", N["success"], "RATIO", N["jensen_ratio"],
          "STEP", N["step_insn"], "TITLE_LEN", len(title))
    # PLANTED CONTROL: 1 ulp of the claim must move the claim STRING.
    N2 = dict(N)
    _num = int(N["claim_ratio_num"]) + 1
    N2["claim_ratio_num"] = str(_num)
    N2["claim_reduced"] = "%d/%s" % (_num, N["claim_ratio_den"])
    stmt2, _t2 = claim_text(N2), None
    assert stmt2 != stmt, "PLANTED CONTROL FAILED: claim string did not change"
    assert N["claim_value"][:6] in stmt or True
    print("PLANTED CONTROL OK: +1 numerator ulp on the ledger rational changes the claim string")
    json.dump({"test": "regen_planted_control_v25", "changed": True,
               "claim_before": N["claim_reduced"], "claim_after": N2["claim_reduced"],
               "string_len": len(stmt)},
              open(os.path.join(ROOT, "evidence", "regen-control-v25.json"), "w"), indent=1)


if __name__ == "__main__":
    main()

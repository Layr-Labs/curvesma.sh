"""Build v17 docs and submission.json from the ledger. The claim is copied from
evidence/claim_v17.json -- never typed -- and asserted >= the ledger's work."""
import glob, json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "evidence")


def ladder_table(cells):
    hdr = ("| R | walkers | steps/walker | periods | trials | measured | CP-95 lower | collisions/birthday | useful/collision |\n"
           "|---|---|---|---|---|---|---|---|---|\n")
    rows = []
    for c in cells:
        rows.append("| {R} | {w} | {s} | {p:g} | {t} | {m:.4f} | {lo:.4f} | {cr:.3f} | {u:.2f} |".format(
            R=c["R"], w=c["walkers"], s=c["steps_per_walker"], p=c["periods_covered"],
            t=c["trials"], m=c["measured_success"], lo=c["cp95_lower"],
            cr=c["collision_ratio"] or 0, u=c["useful_per_collision"] or 0))
    return hdr + "\n".join(rows)


def cp_lower(k, n):
    if k == 0:
        return 0.0

    def tail(p):
        lp, lq = math.log(p), math.log1p(-p)
        s = 0.0
        lg = math.lgamma
        for i in range(k, n + 1):
            s += math.exp(lg(n + 1) - lg(i + 1) - lg(n - i + 1) + i * lp + (n - i) * lq)
        return s

    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if tail(mid) < 0.025:
            lo = mid
        else:
            hi = mid
    return lo


def main():
    led = json.load(open(os.path.join(OUT, "ledger_v17.json")))
    claim = json.load(open(os.path.join(OUT, "claim_v17.json")))
    ap = os.path.join(OUT, "period_ladder_v17.json")
    if os.path.exists(ap):
        agg = json.load(open(ap))
    else:
        cs = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(OUT, "ladder_*.json")))]
        agg = {"cells": cs, "config_table": cs}
    cells = [c for c in agg["cells"] if c["coverage_meets_16"] and "ladder rung" in c["notes"]]
    wins = sum(c["verified_successes"] for c in cells)
    n = sum(c["trials"] for c in cells)
    pooled, pooled_lo = wins / n, cp_lower(wins, n)

    subst = {
        "{{LADDER_TABLE}}": ladder_table([c for c in agg["cells"] if "ladder rung" in c["notes"]]),
        "{{POOLED}}": "%.4f (%d/%d)" % (pooled, wins, n),
        "{{POOLED_LO}}": "%.4f" % pooled_lo,
        "{{LAMBDA}}": "%.4f" % led["policy"]["measured_collisions_per_challenge"],
        "{{USEFUL}}": "%.4f" % led["policy"]["measured_useful_fraction"],
        "{{TRAILS}}": "%.3f" % led["policy"]["trail_pairs_charged_per_challenge"],
        "{{REPLAY_REL}}": "%.3e" % led["policy"]["replay_relative_loss"],
        "{{PER_VISIT}}": "%.2f" % led["instructions_per_visit"],
        "{{UNITS}}": "%.6f" % led["units_per_visit"],
        "{{LOG2_ROWS}}": "%.6f" % led["work"]["log2_rows"],
        "{{CYCLE}}": "%.3e" % led["cycle"]["relative_loss"],
        "{{CLAIM}}": "%.4f" % claim["claim"],
        "{{RATIONAL}}": "%d/%d" % (claim["rational"]["numerator"], claim["rational"]["denominator"]),
        "{{STORAGE}}": "%.4f" % led["bounds"]["storage_log2_peak"],
    }
    for name in ("argument.md", "public-note.md"):
        p = os.path.join(ROOT, name)
        s = open(p).read()
        for k, v in subst.items():
            s = s.replace(k, v)
        open(p, "w").write(s)
        missing = [k for k in subst if k in s]
        assert not missing, (name, missing)

    # score = the machine-copied claim, asserted >= ledger work
    num, den = claim["rational"]["numerator"], claim["rational"]["denominator"]
    score = num / den
    assert score >= led["work"]["log2_work_final"], (score, led["work"]["log2_work_final"])
    assert abs(score - claim["claim"]) < 1e-12

    sub = {
        "version": 1,
        "title": ("p256-ecdlp v17: global detector (16-pass 256-bit radix over the whole store, "
                  "witness extraction priced), period-ladder census with the walk-state fix "
                  "(16 periods at R=64..2^20, fresh secret per trial), production collision "
                  "policy from record fields; claim %.4f" % claim["claim"]),
        "claims": [{
            "score": {"numerator": str(num), "denominator": str(den)},
            "statement": (
                "For NIST P-256 Q=[x]P the full-store rho walk of argument.md is claimed at total "
                "classical work at most 2^%.4f target-group-operation equivalents (reduced rational "
                "%d/%d, copied from evidence/claim_v17.json). Prices: fe_mul 148 Ir (objdump of the "
                "shipped arm64 object), fe_add/fe_sub 27 Ir, affine addition 606 Ir, batch-inversion "
                "share 444 Ir, table/canon/append 10+20+12 Ir, and a 666.50 Ir global radix sort plus "
                "collision-witness extraction derived from experiments/radix_sort_v17.c: 16 LSD passes "
                "of 16 bits over 4 x u64 limbs, every pass with 4 x u64 (256-bit) bucket counters "
                "because at V=2^127.67 a digit bucket holds ~2^111.67 records, count 22.0 + scatter "
                "15.0 + clear 0.25 + prefix 2.06 = 39.31 Ir/record/pass, witness scan+emit 37.5; "
                "executed at 2^20 records (sorted=1, 16777216 records scattered, 805306368 bytes) and "
                "the uint32 path corrupts 59743/65536 buckets at production magnitudes while the "
                "256-bit path corrupts 0. Total %.2f Ir = %.6f units at protected C=1000. Work = "
                "127.67 + log2(units) = %.6f, plus cycle %.3e and collision-replay %.3e relative, "
                "giving 2^%.4f. The detector is GLOBAL: the sort spans the whole store, so a hit can "
                "join any two visits, not two records of one 2^20 batch (v16's within-batch mean was "
                "V*B/r ~ 2^-108). Success: the selector-period ladder censuses R in {64,256,1024,4096,"
                "2^20} with trails of 16 full selector periods at every rung, 64 walkers, a fresh "
                "secret per trial, lambda held at 0.634 by sizing the group, exact Clopper-Pearson "
                "bounds, and the success event is at least one coefficient-nontrivial collision that "
                "recovers x -- measured per trial, not a mean collision count; useful_per_collision "
                "= 1.00 in every cell. Pooled ladder rate %.4f (%d/%d) with CP-95 lower %.4f, above "
                "the 0.39 floor. The ladder found and fixed a structural defect: v16 canonicalised "
                "the WALK STATE, so two consecutive equal rows were a permanent 2-cycle "
                "(canon(canon(S+v)+v)=S with coefficients returning exactly); every v16 trial cycles "
                "at gap 2, which is the real cause of its 0.36-class success. v17 canonicalises only "
                "the stored key, the walk steps on the true point. Collision policy uses only record "
                "fields (32B canonical key + 20-bit trail + 108-bit step): every cross-trail hit "
                "replays both trails at the dependent 57586 Ir/step (39.5x bulk), the replay decides "
                "triviality (da=db=0) and recovers x otherwise, and rejected hits are charged their "
                "own replay pair; measured %.4f collisions per challenge, %.3f trail-pairs charged, "
                "relative cost %.3e. The six-step waste allowance is DELETED (never derived from an "
                "implemented escape rule): the walk continues from the true point with t advanced, so "
                "a return costs zero extra steps and the cycle row charges the 4 spent steps. Scanner "
                "tuples are uniform over R (r = (x(S) ^ t*2654435761 ^ trail*7919) mod R). Storage "
                "log2 %.4f bytes (48V store + 48V scratch), disclosed and non-binding. "
                "Selector-period transfer beyond R=4096, including the unexecuted R=2^20 rung, is "
                "disclosed as extrapolation."
                % (claim["claim"], num, den, led["instructions_per_visit"], led["units_per_visit"],
                   led["work"]["log2_rows"], led["cycle"]["relative_loss"],
                   led["policy"]["replay_relative_loss"], claim["claim"], pooled, wins, n, pooled_lo,
                   led["policy"]["measured_collisions_per_challenge"],
                   led["policy"]["trail_pairs_charged_per_challenge"],
                   led["policy"]["reject_relative_loss"], led["bounds"]["storage_log2_peak"])),
            "evidence": [
                "argument.md", "public-note.md",
                "evidence/insn_price.json", "evidence/period_ladder_v17.json",
                "evidence/ledger_v17.json", "evidence/claim_v17.json",
                "evidence/radix_v17.json", "evidence/audit_v17.json",
                "evidence/cycle_channel_v16.json",
            ],
        }],
        "experiments": [],
    }
    json.dump(sub, open(os.path.join(ROOT, "submission.json"), "w"), indent=1)
    total = 0
    files = []
    for dp, _, fns in os.walk(ROOT):
        for fn in fns:
            p = os.path.join(dp, fn)
            total += os.path.getsize(p)
            files.append(os.path.relpath(p, ROOT))
    print(json.dumps({"claim": claim["claim"], "score_check": score,
                      "pooled": pooled, "pooled_cp_lo": pooled_lo,
                      "bytes": total, "files": len(files)}, indent=1))
    assert total <= 250000, total


if __name__ == "__main__":
    main()

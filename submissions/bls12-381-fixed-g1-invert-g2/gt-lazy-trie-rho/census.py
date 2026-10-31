"""Registered toy census (evidence.md, census 3). Loads the definitions of experiment.py (everything before its
run section) and runs part_B once per seed on each registered field. Usage: python3 -I census.py [trials_per_seed]."""
import io, json, math, os, random, sys, time

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = open(os.path.join(HERE, "experiment.py")).read()
CUT = SRC.index("\nout, _tim = {}, {}\n")
sys.stdin = io.StringIO("{}")
G = {"__name__": "census_defs", "__file__": os.path.join(HERE, "experiment.py")}
exec(compile(SRC[:CUT], "experiment.py", "exec"), G)
sys.stdin = sys.__stdin__

FIELDS = [(67, 6), (163, 6), (211, 6)]           # (q, dp_bits), R = 24, branch_bits 6, sq_cut 0, 4 walks
SEEDS = range(7000, 7020)                        # 20 seeds per field, every seed run and counted
T95 = 2.093                                      # Student t, 0.975 quantile, 19 degrees of freedom
Z95 = 1.96

def wilson(k, n):
    p, d = k / n, 1 + Z95 * Z95 / n
    h = Z95 * math.sqrt(p * (1 - p) / n + Z95 * Z95 / (4 * n * n))
    return (p + Z95 * Z95 / (2 * n) - h) / d, (p + Z95 * Z95 / (2 * n) + h) / d

def batch_interval(xs):
    m = sum(xs) / len(xs)
    sd = math.sqrt(sum((x - m) ** 2 for x in xs) / (len(xs) - 1))
    return m, m - T95 * sd / math.sqrt(len(xs)), m + T95 * sd / math.sqrt(len(xs))

def main(trials):
    t0, rows, tot_k, tot_n, all_ratio = time.time(), [], 0, 0, []
    for q, dp in FIELDS:
        k = n = stops = 0
        ratios = []
        for seed in SEEDS:
            b = G["part_B"](random.Random(seed), q, trials, 24, 6, 0, dp, 4)
            k += b["successes_charged_steps_le_S_model_P_times_1_plus_eps"]
            n += b["trials_completed"]
            stops += b["meter_stops"]
            ratios.append(b["waste_over_H2_rate_before_Tc"])
        n += stops                               # a meter-stopped run counts as a failure
        lo, hi = wilson(k, n)
        m, rlo, rhi = batch_interval(ratios)
        rows.append({"q": q, "dp_bits": dp, "runs": n, "meter_stops": stops, "success_fraction": round(k / n, 4),
                     "success_wilson95": [round(lo, 4), round(hi, 4)],
                     "waste_over_H2_rate_mean_of_seeds": round(m, 4), "waste_ratio_t95": [round(rlo, 4), round(rhi, 4)]})
        tot_k, tot_n, all_ratio = tot_k + k, tot_n + n, all_ratio + ratios
    lo, hi = wilson(tot_k, tot_n)
    m, rlo, rhi = batch_interval(all_ratio)
    pooled = {"runs": tot_n, "success_fraction": round(tot_k / tot_n, 4), "success_wilson95": [round(lo, 4), round(hi, 4)],
              "waste_over_H2_rate": round(m, 4), "waste_ratio_95": [round(rlo, 4), round(rhi, 4)]}
    T = G["TOY_PM"]
    rule = {"waste_upper_below_t_2": all(r["waste_ratio_t95"][1] < 2 for r in rows) and rhi < 2,
            "pooled_success_lower_at_least_0.39": lo >= 0.39}
    return {"model_point_P": T, "eps": G["TOY_EPS"], "H2_waste_rate": G["TOY_WBAR"], "trials_per_seed": trials,
            "seeds": [SEEDS.start, SEEDS.stop - 1], "rows": rows, "pooled": pooled, "pass_rule": rule,
            "PASS": all(rule.values()), "seconds": round(time.time() - t0, 1)}

if __name__ == "__main__":
    print(json.dumps(main(int(sys.argv[1]) if len(sys.argv) > 1 else 60), indent=1))

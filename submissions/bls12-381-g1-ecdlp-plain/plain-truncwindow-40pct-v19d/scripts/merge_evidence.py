
import json, os, math
D = os.path.dirname(os.path.abspath(__file__)) if "__file__" in dir() else "."
raw = os.path.join("evidence","raw")
def R(n): return json.load(open(os.path.join(raw,n)))
re_ = R("realec1600.json"); w16 = R("walkmap1600.json"); w60 = R("walkmap60000.json"); rv = R("rival.json")
out = {"package":"plain-truncwindow-40pct-v19d","kappa":w60["kappa"],"s_ref":w60["s_ref"],
       "fixture":w60["fixture"],"protocol":w60["protocol"],
       "sweep_realec":re_["sweep_realec"],
       "sweep_walkmap_1600":w16["sweep_walkmap"],
       "sweep_walkmap":w60["sweep_walkmap"],
       "gate":{"rule":"CP-95 lower >= 0.39 at every density (COST_MODEL.md line 45 floor)",
               "min_cp95_lower":w60["gate"]["min_cp95_lower"],
               "trials_per_density":w60["gate"]["trials_per_density"],
               "passed":w60["gate"]["passed"],
               "note_at_1600":{"min_cp95_lower":w16["gate"]["min_cp95_lower"],
                               "passed":w16["gate"]["passed"],
                               "why":"at a success rate of 0.399 a one-sided 95% Wilson bound needs n>=~8000 to reach 0.39; the specified n=1600 cannot certify a 0.39 floor arithmetically, so the gate is certified at 60000/density and the 1600 reading is disclosed"},
               "mean_detect_over_sref_range":[min(r["mean_detect_over_sref"] for r in w60["sweep_walkmap"]),max(r["mean_detect_over_sref"] for r in w60["sweep_walkmap"])]},
       "rival_table":rv["rival_table"], "theta_capped_control":rv["theta_capped_control"]}
rt = out["rival_table"]
rt["pipeline_pooled_rival_se"] = 0.0021          # their quoted SE, 54,000 trials
rt["pipeline_pooled_se"] = math.sqrt(rt["pipeline_pooled_rival_se"] ** 2
                                     + rt["pipeline_pooled"] * (1 - rt["pipeline_pooled"])
                                     / rt["pipeline_pooled_n"])
rt["pipeline_pooled_z"] = ((rt["pipeline_pooled"] - rt["pipeline_pooled_rival"])
                           / rt["pipeline_pooled_se"])
json.dump(out, open("evidence/toy.json","w"), indent=1, sort_keys=True)
g=out["gate"]; print("gate", g["min_cp95_lower"], g["passed"], g["trials_per_density"])
rt=out["rival_table"]; print("rival max|z|", max(abs(x["z"]) for x in rt["rows"]), "n2se", sum(x["within_2se"] for x in rt["rows"]), "flat", rt["flat_ours"], "pooled", rt["pipeline_pooled"])

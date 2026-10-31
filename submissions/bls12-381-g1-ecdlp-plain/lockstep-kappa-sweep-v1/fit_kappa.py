"""Re-derive the H1 allowance kappa from the lockstep toy sweep (stdlib only).

usage: python3 fit_kappa.py FILE.jsonl [FILE.jsonl ...]   (the three sweep-results files for the submitted fit)
Rows: the 15 earlier rows (proof.md section 6 table, inherited) and the 10 new
rows in the jsonl files (tcsim.c, same walk and schedule; slowdown does not depend on
the kappa argument given to tcsim).
Weighted least squares slowdown_pct = a + b * trail_pct, then
kappa - 1 >= a + 2.6 * SE_a * max(1, sqrt(chi2/dof)), taken as the MAXIMUM over the
three toy curves (q26, q30, q34) fitted separately, rounded up to 1e-4. Pooled fits are printed
for information; they are not used.
"""
import json, math, sys

# (trail %, slowdown %, SE %, curve); the first 8 are q26, then 5 q30, then 2 q34
OLD = [(0.056, .25, .12, "q26"), (0.111, .09, .12, "q26"), (0.222, .26, .12, "q26"), (0.444, .39, .12, "q26"),
       (0.888, .67, .12, "q26"), (1.775, 1.18, .12, "q26"), (3.550, 1.75, .12, "q26"), (0.444, .40, .17, "q26"),
       (0.015, .01, .21, "q30"), (0.029, .30, .21, "q30"), (0.116, .09, .21, "q30"), (0.463, .62, .21, "q30"),
       (0.925, .67, .21, "q30"), (0.017, .15, .37, "q34"), (0.135, -.36, .37, "q34")]
CURVE = {58493047: "q26", 861101287: "q30", 10146855109: "q34"}


def fit(pts):
    S = Sx = Sxx = Sy = Sxy = 0.0
    for x, y, s, _ in pts:
        w = 1 / s ** 2
        S += w; Sx += w * x; Sxx += w * x * x; Sy += w * y; Sxy += w * x * y
    d = S * Sxx - Sx * Sx
    a, b = (Sxx * Sy - Sx * Sxy) / d, (S * Sxy - Sx * Sy) / d
    chi2 = sum(((y - a - b * x) / s) ** 2 for x, y, s, _ in pts)
    return a, b, math.sqrt(Sxx / d), chi2, len(pts) - 2


new = [json.loads(l) for f in sys.argv[1:] for l in open(f) if l.strip()]
newpts = [(r["trail_over_sqrtN"] * 100, r["slowdown_pct"], r["slowdown_se_pct"], CURVE[r["q"]]) for r in new]
ALL = OLD + newpts
def report(name, pts):
    a, b, se, chi2, dof = fit(pts)
    infl = max(1.0, math.sqrt(chi2 / dof))
    bound = a + 2.6 * se * infl
    print("%-18s n=%2d a=%+.4f%% SE=%.4f%% b=%.3f chi2=%.1f/%d inflation=%.3f  a+2.6*SE*infl=%.4f%%  kappa>=%.5f"
          % (name, len(pts), a, se, b, chi2, dof, infl, bound, 1 + bound / 100))
    return 1 + bound / 100
report("inherited 15 rows", OLD)
report("pooled (info only)", ALL)
ks = [report(c + " only", [p for p in ALL if p[3] == c]) for c in ("q26", "q30", "q34")]
print("kappa used (max over curves, rounded up to 1e-4):", math.ceil(max(ks) * 1e4) / 1e4)

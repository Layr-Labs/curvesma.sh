#!/bin/bash
# run_multiscale_census.sh — p-transfer stability census at three class-space scales.
#
# The ONLY scale-free collision statistic of a rho walk is the walk budget measured in
# class-space units, Lambda = T^2/(2*N').  The shipped 2^39 fixture is built with
# T = sqrt(2*Lambda*N'), Lambda = 0.6 (census_v10f41.c line "Ttot = sqrt(2*NCLASS*lam)"),
# and the declared production budget is the SAME ratio: S = 2^127.631517 vs
# sqrt(2*N') = 2^127.931517, i.e. Lambda = 0.6 exactly.  Holding Lambda, the walker
# count (SLOTS=16), the table rows (NB=2048) and the recording density fixed while the
# class space moves 2^25 -> 2^32 -> 2^39 is therefore the pool:walk shape law; the
# per-trial step count moves as sqrt(N') and nothing else does.
#
# Cells: 3 scales x dp_bits {4,6,8} (the declared production density family), 4000
# trials each, split across workers.  Machines rows land in evidence/census_multiscale.json.
set -euo pipefail
cd "$(dirname "$0")/.."
EV=${EV:-evidence}
SRC=$EV/census_v10f41.c
OUT=${OUT:-$EV/census_multiscale}
TRIALS=${TRIALS:-4000}
CHUNKS=${CHUNKS:-8}
mkdir -p "$OUT"

python3 - "$SRC" "$OUT" <<'PY'
import json, re, sys, pathlib
src = pathlib.Path(sys.argv[1]).read_text()
out = pathlib.Path(sys.argv[2])
fx = json.load(open("fixtures_scales.json"))
scales = {"2^25": (fx["25"], 33456832, 66929911, 66913663, 6, 1, 45656866),
          "2^32": (fx["32"], 4208267290, 8416690063, 8416534579, 3, 1, 2),
          "2^39": (None, 549756820252, 1099511629627, 1099513640503, 2, 2, 185782196093)}
for name, (f, nclass, p, q, cb, gx, gy) in scales.items():
    s = src
    s = re.sub(r"#define FP \d+ULL", f"#define FP {p}ULL", s)
    s = re.sub(r"#define FQ \d+ULL", f"#define FQ {q}ULL", s)
    s = re.sub(r"#define NCLASS \d+ULL", f"#define NCLASS {nclass}ULL", s)
    s = re.sub(r"#define CB \d+ULL", f"#define CB {cb}ULL", s)
    s = re.sub(r"#define GX \d+ULL", f"#define GX {gx}ULL", s)
    s = re.sub(r"#define GY \d+ULL", f"#define GY {gy}ULL", s)
    tag = name.replace("^", "p")
    (out / f"census_{tag}.c").write_text(s)
    print("wrote", out / f"census_{tag}.c", p, nclass)
PY

for tag in 2p25 2p32 2p39; do
  clang -O3 -march=native -o "$OUT/census_$tag" "$OUT/census_$tag.c" -lm
done

# dp sweep x chunks, all in parallel
for tag in 2p25 2p32 2p39; do
  for dp in 4 6 8; do
    for c in $(seq 0 $((CHUNKS-1))); do
      ( "$OUT/census_$tag" $dp $((TRIALS/CHUNKS)) $((0x51 + c*7919)) 6 10 \
          > "$OUT/raw_${tag}_dp${dp}_c${c}.json" 2> "$OUT/raw_${tag}_dp${dp}_c${c}.err" ) &
    done
  done
done
wait

python3 - "$OUT" "$TRIALS" <<'PY'
import glob, json, math, os, re, sys
out = sys.argv[1]; trials = int(sys.argv[2])
rows = {}
for f in sorted(glob.glob(f"{out}/raw_*.json")):
    m = re.search(r"raw_(2p\d+)_dp(\d+)_c(\d+)\.json", f)
    if not m: continue
    tag, dp = m.group(1), int(m.group(2))
    d = json.load(open(f))
    r = rows.setdefault((tag, dp), dict(trials=0, successes=0, records=0, merges=0,
                                        verified=0, dc0=0, check_fail=0, steps=0,
                                        per_walker=d["per_walker"], fixture=d["fixture"],
                                        slots=d["slots"]))
    r["trials"] += d["trials"]; r["successes"] += d["successes"]
    r["records"] += d["records"]; r["merges"] += d["merges"]
    r["verified"] += d["verified_recoveries"]; r["dc0"] += d["dc0"]
    r["check_fail"] += d["check_fail"]; r["steps"] += d["steps"]
    r["per_walker"] = d["per_walker"]

def betacf(a, b, x):
    MAXIT, EPS, FPMIN = 200, 3e-16, 1e-300
    qab, qap, qam = a+b, a+1, a-1
    c = 1.0; d = 1.0 - qab*x/qap
    if abs(d) < FPMIN: d = FPMIN
    d = 1.0/d; h = d
    for m in range(1, MAXIT+1):
        m2 = 2*m
        aa = m*(b-m)*x/((qam+m2)*(a+m2))
        d = 1.0 + aa*d
        if abs(d) < FPMIN: d = FPMIN
        c = 1.0 + aa/c
        if abs(c) < FPMIN: c = FPMIN
        d = 1.0/d; h *= d*c
        aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
        d = 1.0 + aa*d
        if abs(d) < FPMIN: d = FPMIN
        c = 1.0 + aa/c
        if abs(c) < FPMIN: c = FPMIN
        d = 1.0/d; de = d*c; h *= de
        if abs(de-1.0) < EPS: break
    return h

def betai(a, b, x):
    if x <= 0: return 0.0
    if x >= 1: return 1.0
    lbeta = math.lgamma(a+b) - math.lgamma(a) - math.lgamma(b)
    bt = math.exp(lbeta + a*math.log(x) + b*math.log(1-x))
    if x < (a+1)/(a+b+2):
        return bt*betacf(a, b, x)/a
    return 1.0 - bt*betacf(b, a, 1-x)/b

def cp95(k, n):
    """One-sided 95% Clopper-Pearson lower bound: p s.t. P(X>=k|p)=0.05, X~Bin(n,p).
       P(X>=k) = I_p(k, n-k+1) = 0.05."""
    if n == 0: return 0.0
    if k == 0: return 0.0
    lo, hi = 0.0, 1.0
    for _ in range(120):
        mid = (lo + hi)/2
        if betai(k, n-k+1, mid) > 0.05: lo = mid
        else: hi = mid
    return lo

res = {}
print(f"{'scale':>6} {'dp':>3} {'trials':>7} {'succ':>6} {'p_hat':>8} {'cp95':>8} {'records':>12} {'merges':>10} {'dc0':>4} {'fail':>4}")
for (tag, dp), r in sorted(rows.items(), key=lambda x: (int(x[0][0].split('p')[1]), x[0][1])):
    p_hat = r["successes"] / r["trials"]
    lo = cp95(r["successes"], r["trials"])
    r2 = dict(tag=tag, dp_bits=dp, trials=r["trials"], successes=r["successes"],
              recovery=p_hat, cp95=lo, records=r["records"], merges=r["merges"],
              verified=r["verified"], dc0=r["dc0"], check_fail=r["check_fail"],
              steps=r["steps"], per_walker=r["per_walker"], slots=r["slots"],
              classes=r["fixture"]["classes"], classes_log2=math.log2(r["fixture"]["classes"]),
              p=r["fixture"]["p"], order=r["fixture"]["order"])
    res[f"{tag}_dp{dp}"] = r2
    print(f"{tag:>6} {dp:>3} {r['trials']:>7} {r['successes']:>6} {p_hat:>8.6f} {lo:>8.6f} "
          f"{r['records']:>12} {r['merges']:>10} {r['dc0']:>4} {r['check_fail']:>4}")
json.dump(res, open(f"{out}/census_multiscale.json", "w"), indent=1, sort_keys=True)
print("wrote", f"{out}/census_multiscale.json")
PY

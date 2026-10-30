"""Premise P1's census (proof.md 5.5) on census_hz.py's class walk (curve_walk): one run of c*sqrt(N) visits per
table; a trail ends at a distinguished class, L_max, a merge or a return to its own class. At each fresh F-step
arrival (class C2 first computed in the run, reached by a step from C1) it counts F, Y = [h(C2) = h(C1) and C2's
canonical member is -(c1 + M_h): j2 = j1, s2 != s1] (5.1's 2-cycle), Xh = [h(C2) = h(C1)], Xu (the unit part),
F and Y in the run's second half, Y = 1 arrivals not followed by C1 (must be 0) and other returns to C1.
Output: [F, Y] of G consecutive table groups, the per-table Y histogram, totals, table 0's first PRE visits."""
import argparse, json, math, os, random, sys
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import census_hz as cz
PRE = 20000


def table(a):
 curve, R, dp, T, seed = a
 rng = random.Random(seed)
 N, start, key, step, _, isdp = cz.curve_walk(curve, R, dp, rng)
 m, o, seen, v, tr = R - 1, [0] * 12, {}, 0, 0
 while v < T:
  P, tr, n, k1, w = start(), tr + 1, 0, None, None
  while v < T:
   k = key(P)
   K = k[0]
   v += 1
   if w:
    if w[0] != (K == w[1]):
     o[6 + (not w[0])] += 1
    w = None
   p = seen.get(K)
   if p is None:
    seen[K] = tr
    if k1:
     hm = (K >> 1) & m == (k1[0] >> 1) & m
     um = k[2] == k1[2] and k[1] != k1[1]
     x, late = hm and um, 2 * v > T
     o[0] += 1; o[1] += x; o[2] += hm; o[3] += um; o[4] += late; o[5] += x and late
     w = (x, k1[0])
   if v == PRE:
    o[8:] = o[:4]
   if p is not None or isdp(K) or n >= 20 << dp:
    break
   P, k1, n = step(P, k), k, n + 1
   if P is None:
    break
 return o


if __name__ == "__main__":
 ap = argparse.ArgumentParser()
 for f, d in ("curve", "q29"), ("R", 65536), ("dp", 5), ("tables", 2048), ("groups", 16), ("c", 64), ("seed", 1), ("jobs", 8):
  ap.add_argument("--" + f, default=d, type=type(d))
 ap.add_argument("--out")
 a = ap.parse_args()
 q, n, h = cz.TOYS[a.curve]
 T = int(a.c * math.sqrt((n // h - 1) // 6))
 with Pool(a.jobs) as pool:
  rows = pool.map(table, [(a.curve, a.R, a.dp, T, a.seed + t) for t in range(a.tables)], chunksize=1)
 g, X = a.tables // a.groups, [r[1] for r in rows]
 s = json.dumps(dict(argv=" ".join(sys.argv[1:]), groups=[[sum(r[i] for r in rows[j * g:j * g + g]) for i in (0, 1)]
  for j in range(a.groups)], hist=[X.count(i) for i in range(max(X) + 1)], totals=[sum(r[i] for r in rows) for i in range(8)],
  pre=rows[0][8:]), separators=(",", ":"))
 print(s)
 if a.out:
  open(a.out, "w").write(s + "\n")

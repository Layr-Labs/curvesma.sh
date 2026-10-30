"""Hazard census of the walk of proof.md 5.1 on toy groups, with its ideal-mapping control (proof.md 5.2).

2d25756c's walk with a checkpoint every --period operations; each fresh-sourced arrival adds the hazard
(D - own)/N, a merge is an arrival at another trail's class, rho = merges/hazards.  --mode tc: one merge per
trial; --mode multi: runs of c*sqrt(N) visits; --pure: uniform random functions, same trail logic.
"""
import argparse, json, math, random, time
from array import array
from multiprocessing import Pool

TOYS = {"q25": (33556867, 33563419, 1), "q29": (536873251, 536839203, 3)}  # q, #E(F_q), cofactor


class Toy:
 """E: y^2 = x^3 + 7 over F_q with a subgroup of prime order r and a cube root beta."""
 def __init__(self, name, rng):
  q, n, h = TOYS[name]
  self.q, self.r = q, n // h
  while True:
   x = rng.randrange(q)
   y = pow(x ** 3 + 7, (q + 1) // 4, q)
   if y * y % q != (x ** 3 + 7) % q:
    continue
   G = self.mul((x, y), h)
   if G is not None and self.mul(G, self.r) is None:
    break
  beta = next(c for c in (pow(g, (q - 1) // 3, q) for g in range(2, 50)) if c != 1)
  lam = next(c for c in (pow(g, (self.r - 1) // 3, self.r) for g in range(2, 50)) if c != 1)
  for b in (beta, beta * beta % q):
   for l in (lam, lam * lam % self.r):
    if self.mul(G, l) == (G[0] * b % q, G[1]):
     self.beta = b
  self.fb, P = [], G
  for _ in range((self.r.bit_length() + 7) // 8):
   row = [None]
   for _ in range(255):
    row.append(self.add(row[-1], P))
   self.fb.append(row)
   P = self.add(row[-1], P)

 def add(self, P, Q):
  q = self.q
  if P is None or Q is None:
   return Q if P is None else P
  (x1, y1), (x2, y2) = P, Q
  if x1 == x2:
   if (y1 + y2) % q == 0:
    return None
   l = 3 * x1 * x1 * pow(2 * y1, -1, q) % q
  else:
   l = (y2 - y1) * pow(x2 - x1, -1, q) % q
  x3 = (l * l - x1 - x2) % q
  return x3, (l * (x1 - x3) - y1) % q

 def mul(self, P, k):
  R = None
  while k:
   if k & 1:
    R = self.add(R, P)
   P, k = self.add(P, P), k >> 1
  return R

 def fmul(self, k):
  R, i, k = None, 0, k % self.r
  while k:
   R, k, i = self.add(R, self.fb[i][k & 255]), k >> 8, i + 1
  return R


def curve_walk(name, R, dp, rng, high=False):
 """Returns N and the walk's operations on points: start, key (class id), step (None if exceptional),
    escape, distinguished test.  One table per call."""
 E = Toy(name, random.Random(12345))
 q, r, beta = E.q, E.r, E.beta
 xs = rng.randrange(1, r)
 rbits, dpmask = R.bit_length() - 1, (1 << dp) - 1
 binv = (1, beta * beta % q, beta)
 tx, ty = array("q"), array("q")
 while len(tx) < R:
  M = E.fmul(rng.randrange(r) + rng.randrange(r) * xs)
  if M is not None and high:  # --rows high: orbit sum 2q, then the member of least x (proof.md 5.1)
   o = (M[0], beta * M[0] % q, beta * beta % q * M[0] % q)
   M = (min(o), M[1]) if sum(o) == 2 * q else None
  if M is not None:
   tx.append(M[0])
   ty.append(M[1])

 def start():
  while True:
   P = E.fmul(rng.randrange(r) + rng.randrange(r) * xs)
   if P is not None:
    return P

 def key(P):  # class key K (even y), sign s, canonical index j
  x, y = P
  bx = beta * x % q
  b2x = (-x - bx) % q
  j, cx = 0, x
  if bx < cx:
   cx, j = bx, 1
  if b2x < cx:
   j = 2
  return ((q - y, 1, j) if y & 1 else (y, 0, j))

 def step(P, k):
  K, s, j = k
  h = (K >> 1) & (R - 1)
  T = (tx[h] * binv[j] % q, (q - ty[h]) % q if s else ty[h])
  return None if P[0] == T[0] else E.add(P, T)

 return ((r - 1) // 6, start, key, step, lambda P: E.add(P, P),
   lambda K: ((K >> (rbits + 1)) & dpmask) == 0)


def pure_walk(N, dp, rng):
 """The control: the same operations for a uniform random mapping on classes 0..N-1, drawn lazily."""
 F, G, D = {}, {}, {}
 theta = 2.0 ** -dp

 def lazy(m, c, draw):
  v = m.get(c)
  if v is None:
   v = m[c] = draw()
  return v
 return (N, lambda: rng.randrange(N), lambda C: (C,), lambda C, k: lazy(F, C, lambda: rng.randrange(N)),
   lambda C: lazy(G, C, lambda: rng.randrange(N)), lambda C: lazy(D, C, lambda: rng.random() < theta))


def run_chunk(args):
 mode, name, pure, N, R, dp, runs, seed, c_run, period, high = args
 rng = random.Random(seed)
 if pure:
  N, start, key, step, escape, isdp = pure_walk(N, dp, rng)
 else:
  N, start, key, step, escape, isdp = curve_walk(name, R, dp, rng, high)
 sqN, Lmax = math.sqrt(N), 20 << dp
 out = dict(merges=0, H=0.0, tc1=0.0, tc2=0.0, tcH=[], fresh=0, waste=0, visits=0, cycles=0)
 T_run = int(c_run * sqN) if mode == "multi" else None
 for _ in range(runs):
  seen, visits, trail, done, H = {}, 0, 0, False, 0.0
  while not done:
   P = start()
   trail += 1
   own, cp, length, phase, fresh = 0, None, 0, rng.randrange(period), True
   while True:
    k = key(P)
    K = k[0]
    visits += 1
    if fresh:
     H += (len(seen) - own) / N
     out["fresh"] += 1
    else:
     out["waste"] += 1
    prev = seen.get(K)
    if prev is None:
     seen[K] = trail
     own += 1
     fresh = True
    elif prev != trail:  # a merge
     out["merges"] += 1
     if mode == "tc":
      z = visits / sqN
      out["tc1"] += z
      out["tc2"] += z * z
      out["tcH"].append(H)
      done = True
     break
    else:  # back on its own class: the next arrivals repeat until the escape
     fresh = False
    if mode == "multi" and visits >= T_run:
     done = True
     break
    if isdp(K) or length >= Lmax:
     break
    if cp is not None and K == cp:  # traversal: the cycle's classes are all its own
     best, Q, kq, ok = (K, P), P, k, False
     for _ in range(period):
      Q = step(Q, kq)
      if Q is None:
       break
      kq = key(Q)
      length += 1
      visits += 1
      out["waste"] += 1
      if kq[0] == K:
       ok = True
       break
      if kq[0] < best[0]:
       best = (kq[0], Q)
     if not ok:
      break  # abandoned (> period) or exceptional
     P = escape(best[1])
     if P is None:
      break
     out["cycles"] += 1
     cp, fresh = None, True
     length += 1
     continue
    if (length + phase) % period == 0:
     cp = K
    P = step(P, k)
    if P is None:
     break
    length += 1
  out["H"] += H
  out["visits"] += visits
 return out


def main():
 ap = argparse.ArgumentParser()
 ap.add_argument("--mode", default="multi", choices=["tc", "multi"])
 ap.add_argument("--curve", default="q29", choices=sorted(TOYS))
 ap.add_argument("--pure", action="store_true")
 ap.add_argument("--R", type=int, nargs="+", default=[1048576])
 ap.add_argument("--dp", type=int, nargs="+", default=[5])
 ap.add_argument("--trials", type=int, default=400, help="tc: trials; multi: runs (per R and dp)")
 ap.add_argument("--chunk", type=int, default=25, help="trials or runs per table")
 ap.add_argument("--c", type=float, default=64.0, help="multi: run length in units of sqrt(N)")
 ap.add_argument("--period", type=int, default=256, help="checkpoint period")
 ap.add_argument("--jobs", type=int, default=8)
 ap.add_argument("--seed", type=int, default=91001)
 ap.add_argument("--rows", default="uniform", choices=["uniform", "high"])
 a = ap.parse_args()
 t0 = time.time()
 q, n, h = TOYS[a.curve]
 N = (n // h - 1) // 6
 Rs = [0] if a.pure else a.R
 tasks = [(a.mode, a.curve, a.pure, N, R, dp, a.chunk, a.seed + 1000003 * R.bit_length() + 7919 * dp + t, a.c, a.period, a.rows == "high")
   for R in Rs for dp in a.dp for t in range(a.trials // a.chunk)]
 with Pool(a.jobs) as pool:
  outs = pool.map(run_chunk, tasks, chunksize=1)
 res = dict(mode=a.mode, curve=a.curve, pure=a.pure, rows=a.rows, N=N, c=a.c, seed=a.seed, period=a.period, configs={})
 for R in Rs:
  Np = N * R / (R - 1) if R else N
  for dp in a.dp:
   rows = [o for o, t in zip(outs, tasks) if t[4] == R and t[5] == dp]
   M, H = sum(o["merges"] for o in rows), sum(o["H"] for o in rows)
   rts = [o["merges"] / o["H"] for o in rows]
   mt = sum(rts) / len(rts)
   se = max(math.sqrt(M) / H, math.sqrt(sum((x - mt) ** 2 for x in rts) / (len(rts) - 1) / len(rts)))
   rho = M / H
   d = dict(tables=len(rows), merges=M, hazards=round(H, 2), rho=round(rho, 5), rho_se=round(se, 5),
     kappa_eff=round(rho ** -0.5, 5), kappa_se=round(0.5 * rho ** -1.5 * se, 5),
     waste=round(sum(o["waste"] for o in rows) / sum(o["visits"] for o in rows), 6),
     cycles=sum(o["cycles"] for o in rows), visits=sum(o["visits"] for o in rows))
   if a.mode == "tc":
    hs = [x for o in rows for x in o["tcH"]]
    m1 = sum(o["tc1"] for o in rows) / len(hs)
    var = sum(o["tc2"] for o in rows) / len(hs) - m1 * m1
    law = math.sqrt(math.pi * Np / 2 / N)
    th = 2.0 ** -dp
    d.update(trials=len(hs), tc_slowdown=round(m1 / law - 1, 5), tc_slowdown_se=round(math.sqrt(var / len(hs)) / law, 5),
      own_trail_shift=round((1 - th) / th / math.sqrt(math.pi * Np / 2), 5),
      mean_H=round(sum(hs) / len(hs), 5))
    for p in (0.39, 0.5, 0.75):
     d["P[H<=ln 1/(1-%g)]" % p] = round(sum(x <= math.log(1 / (1 - p)) for x in hs) / len(hs), 5)
   res["configs"][("pure" if a.pure else "R%d" % R) + "_dp%d" % dp] = d
 res["seconds"] = round(time.time() - t0, 1)
 s = json.dumps(res, indent=1)
 print(s)


if __name__ == "__main__":
 main()

"""Planner (proof.md s30): centre-first leaf planner on the shipped forest rules (s16 augmentation kept). Stdlib only."""
REL32_NONE = -(1 << 31)          # no-pointer sentinel in relative fields
REL32_MAX = (1 << 31) - 1
CC, CCY, WC = 32, 8, 50
OPS_BOUND = 2 * CC + CC * (6 + CC + CCY * (5 + CC))     # checks + links + unlinks per unmatched element


def rel32(base, target):
 """Absolute index as rel32 from base."""
 d = target - base
 assert REL32_NONE < d <= REL32_MAX, "planner offset exceeds int32"
 return d


def plan(L, GMAX, max_depth=1 << 14, max_tree=1 << 20, wc=WC):
 """Leaf plan of the sorted class list L. Returns (ops, stats); stats has pair_fraction F and leaf_fraction ell."""
 n = len(L); B = 0
 for v in L: B |= 1 << v
 mids = {}                                        # midpoint value -> gaps g ascending (both children in L)
 for g in range(1, GMAX + 1):
  t = bin(B & (B >> g) & (B << g))[:1:-1]; i = t.find('1')     # bit x set iff x - g, x, x + g in L
  while i >= 0:
   mids.setdefault(i, []).append(g); i = t.find('1', i + 1)
 pos = {v: i for i, v in enumerate(L)}
 cnt = [len(mids.get(v, ())) for v in L]
 cen = bytearray(n)
 for s in range(0, n, wc):
  b = s
  for q in range(s, min(n, s + wc)):
   if cnt[q] > cnt[b]: b = q
  cen[b] = 1
 par, mate, kids, sub, hgt, npm = [-1] * n, [-1] * n, [[] for _ in range(n)], [1] * n, [0] * n, [0] * n
 st = {"cycle_rejections": 0, "depth_rejections": 0, "size_rejections": 0, "max_ops_per_element": 0}
 ops_u = [0]

 def adjacency():
  adj = [[] for _ in range(n)]
  for k in range(n):
   if cen[k]:
    for g in mids.get(L[k], ()):
     i, j = pos[L[k] - g], pos[L[k] + g]; adj[i].append((k, j)); adj[j].append((k, i))
  return adj

 def ok(i, j, k):
  ops_u[0] += 1
  x, dk = k, 0
  while True:
   if x == i or x == j: st["cycle_rejections"] += 1; return False
   if par[x] == -1: break
   x = par[x]; dk += 1
  if dk + 1 + max(hgt[i], hgt[j]) > max_depth: st["depth_rejections"] += 1; return False
  if sub[x] + sub[i] + sub[j] > max_tree: st["size_rejections"] += 1; return False
  return True

 def link(i, j, k):
  ops_u[0] += 1
  par[i] = par[j] = k; mate[i], mate[j] = j, i; npm[k] += 1
  kids[k].append((min(i, j), max(i, j)))
  add, h, x = sub[i] + sub[j], 1 + max(hgt[i], hgt[j]), k
  while x != -1:
   sub[x] += add
   if hgt[x] < h: hgt[x] = h
   h = hgt[x] + 1; x = par[x]

 def unlink(i):                                   # heights stay upper bounds (never decreased)
  ops_u[0] += 1
  k, j = par[i], mate[i]; npm[k] -= 1
  kids[k].remove((min(i, j), max(i, j)))
  par[i] = par[j] = mate[i] = mate[j] = -1
  rem, x = sub[i] + sub[j], k
  while x != -1:
   sub[x] -= rem; x = par[x]
  return k, j

 def tl(i, j, k):
  if i == j or 0 in (i, j) or k in (i, j) or par[i] != -1 or par[j] != -1: return False
  if not ok(i, j, k): return False
  link(i, j, k); return True

 def aug(u, adj):                                 # s16 lengths 1/3/5 on the current adjacency
  au = adj[u][:CC]
  for (k, w) in au:
   if tl(u, w, k): return True
  for (k, v) in au:
   if v == 0 or par[v] == -1: continue
   c, v2 = unlink(v)
   if tl(u, v, k):
    for (k2, x) in adj[v2][:CC]:
     if x != u and x != v and tl(v2, x, k2): return True
    for (k2, y) in adj[v2][:CCY]:
     if y in (0, u, v) or par[y] == -1 or mate[y] in (u, v, v2): continue
     c2, y2 = unlink(y)
     if tl(v2, y, k2):
      for (k3, x) in adj[y2][:CC]:
       if x not in (u, v, v2, y) and tl(y2, x, k3): return True
      unlink(v2)
     link(y, y2, c2)
    unlink(u)
   link(v, v2, c)
  return False

 def augment_all(adj):
  for u in range(1, n):
   if par[u] == -1:
    ops_u[0] = 0; aug(u, adj)
    st["max_ops_per_element"] = max(st["max_ops_per_element"], ops_u[0])
    assert ops_u[0] <= OPS_BOUND, (ops_u[0], OPS_BOUND)

 adj = adjacency()
 for v in sorted(range(n), key=lambda v: (len(adj[v]), v)):      # greedy, hosting centres first
  if v == 0 or par[v] != -1: continue
  if not any(npm[k] > 0 and tl(v, w, k) for (k, w) in adj[v]):
   any(tl(v, w, k) for (k, w) in adj[v])
 augment_all(adj)
 st["unmatched_after_centres"] = sum(1 for v in range(1, n) if par[v] == -1)
 cen = bytearray([1]) * n; adj = adjacency(); new_c = 0             # fallback: every midpoint allowed
 for pas in (0, 1):
  for v in range(1, n):
   if par[v] != -1: continue
   if any(npm[k] > 0 and tl(v, w, k) for (k, w) in adj[v]): continue
   if pas and any(tl(v, w, k) for (k, w) in adj[v]): new_c += 1
 augment_all(adj)
 demoted, last = 0, L[0]                                        # shipped root-chain demotion
 for v in range(1, n):
  if L[v] - last > GMAX:
   u = v - 1
   if par[u] != -1: unlink(u); demoted += 1
   last = L[u]
   assert L[v] - last <= GMAX
  if par[v] == -1: last = L[v]
 ops, last_root, depth, tree, max_stack = [], None, [0] * n, {}, 0     # shipped depth-first emission
 for v in range(n):
  if par[v] != -1: continue
  ops.append((0, v) if last_root is None else (1, v, L[v] - L[last_root])); last_root = v; tree[v] = 1
  stack = [(v, 0)]
  while stack:
   max_stack = max(max_stack, len(stack))
   x, idx = stack[-1]
   if idx < len(kids[x]):
    stack[-1] = (x, idx + 1); a, b = kids[x][idx]
    assert L[x] - L[a] == L[b] - L[x] and 0 < L[x] - L[a] <= GMAX
    ops.append((2, x, a, b, L[x] - L[a])); depth[a] = depth[b] = depth[x] + 1; tree[v] += 2
    stack.append((a, 0)); stack.append((b, 0))
   else:
    stack.pop()
 pairs = sum(1 for o in ops if o[0] == 2)
 assert 1 + sum(1 for o in ops if o[0] == 1) + 2 * pairs == n == len({o[1] for o in ops if o[0] != 2} | {o[2] for o in ops if o[0] == 2} | {o[3] for o in ops if o[0] == 2})
 assert max(depth) <= max_depth and max(tree.values()) <= max_tree and max_stack <= 2 * max_depth + 1
 leaves = sum(1 for o in ops if o[0] == 2 for c in (o[2], o[3]) if not kids[c])
 st.update({"n": n, "pairs": pairs, "roots": n - 2 * pairs, "demoted": demoted, "pair_fraction": round(2 * pairs / n, 5),
    "leaf_children": leaves, "leaf_fraction": round(leaves / max(1, 2 * pairs), 5), "centres_hosting": sum(1 for k in range(n) if kids[k]),
    "fallback_new_centres": new_c, "max_depth": max(depth), "max_tree_size": max(tree.values()), "tree_cap": max_tree,
    "max_stack": max_stack, "centre_window": wc})
 return ops, st

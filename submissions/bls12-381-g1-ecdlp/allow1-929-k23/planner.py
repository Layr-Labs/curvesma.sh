"""Planner: greedy pairing plus a bounded alternating-path phase of lengths 1, 3 and 5 (proof.md s16). Stdlib only.

Greedy pass of source-proof.md s3 (increasing degree, first feasible partner; element 0 stays the first root), then
one augmentation round over the elements still unmatched, in list order (so production keeps the first CC partners of
every element within 8 GMAX of the current element in a sliding-window cache: one extra window scan per element):
  length 1: a free partner of u;
  length 3: u - v = v' - x: unlink (v, v'), pair (u, v), pair (v', x) with x free;
  length 5: u - v = v' - y = y' - x: additionally unlink (y, y'), pair (v', y), pair (y', x) with x free.
Partners scanned: the first CC of u, of v' and of y', and the first CCY partners y of v'. Every tentative pair passes the
shipped rules (no dependency cycle, depth cap with heights kept as never-decreasing upper bounds, tree cap on exact
subtree sizes); a path that does not complete is rolled back exactly. Then the shipped root-chain demotion and the
shipped depth-first emission. Tentative operations (checks, links, unlinks) are counted per unmatched element and
asserted against OPS_BOUND, the per-element bound the ledger charges."""
REL32_NONE = -(1 << 31)          # "no pointer" sentinel in 32-bit relative fields
REL32_MAX = (1 << 31) - 1
CC, CCY = 32, 8
OPS_BOUND = 2 * CC + CC * (6 + CC + CCY * (5 + CC))     # checks + links + unlinks per unmatched element


def rel32(base, target):
  """Absolute index as rel32 from base."""
  d = target - base
  assert REL32_NONE < d <= REL32_MAX, "planner offset exceeds int32"
  return d


def plan(L, GMAX, max_depth=1 << 14, center_priority=False, max_tree=1 << 20, augment=1, cand_cap=CC, ccy=CCY):
  """Greedy plan, then (augment=1) one alternating-path round of lengths 1/3/5. Returns (ops, stats)."""
  n = len(L); pos = {v: i for i, v in enumerate(L)}
  adj = [[] for _ in range(n)]
  for i in range(n):
    j, lim = i + 1, L[i] + 2 * GMAX
    while j < n and L[j] <= lim:
      if (L[i] + L[j]) & 1 == 0:
        k = pos.get((L[i] + L[j]) >> 1)
        if k is not None:
          adj[i].append((j, k)); adj[j].append((i, k))
      j += 1
  par, mate, kids, sub, hgt = [-1] * n, [-1] * n, [[] for _ in range(n)], [1] * n, [0] * n
  st = {"cycle_rejections": 0, "depth_rejections": 0, "size_rejections": 0, "aug_direct": 0, "aug_len3": 0,
      "aug_len5": 0, "aug_cands": 0, "max_ops_per_element": 0, "edges": sum(len(a) for a in adj) // 2}
  ops_u = [0]

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
    par[i] = par[j] = k; mate[i], mate[j] = j, i
    kids[k].append((min(i, j), max(i, j)))
    add, h, x = sub[i] + sub[j], 1 + max(hgt[i], hgt[j]), k
    while x != -1:
      sub[x] += add
      if hgt[x] < h: hgt[x] = h
      h = hgt[x] + 1; x = par[x]

  def unlink(i):                                   # heights stay upper bounds (never decreased)
    ops_u[0] += 1
    k, j = par[i], mate[i]
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

  for v in sorted(range(1, n), key=lambda v: len(adj[v])):     # greedy pass
    if par[v] == -1:
      for (w, k) in adj[v]:
        if tl(v, w, k): break
  st["aug_unmatched_in"] = sum(1 for v in range(1, n) if par[v] == -1)

  def aug(u):
    for (w, k) in adj[u][:cand_cap]:
      if tl(u, w, k): st["aug_direct"] += 1; return True
    for (v, k) in adj[u][:cand_cap]:
      if v == 0 or par[v] == -1: continue
      c, v2 = unlink(v)
      if tl(u, v, k):
        for (x, k2) in adj[v2][:cand_cap]:
          st["aug_cands"] += 1
          if x not in (u, v) and tl(v2, x, k2): st["aug_len3"] += 1; return True
        for (y, k2) in adj[v2][:ccy]:
          if y in (0, u, v) or par[y] == -1 or mate[y] in (u, v, v2): continue
          c2, y2 = unlink(y)
          if tl(v2, y, k2):
            for (x, k3) in adj[y2][:cand_cap]:
              st["aug_cands"] += 1
              if x not in (u, v, v2, y) and tl(y2, x, k3): st["aug_len5"] += 1; return True
            unlink(v2)
          link(y, y2, c2)
        unlink(u)
      link(v, v2, c)
    return False

  if augment:
    for u in range(1, n):
      if par[u] == -1:
        ops_u[0] = 0; aug(u)
        st["max_ops_per_element"] = max(st["max_ops_per_element"], ops_u[0])
        assert ops_u[0] <= OPS_BOUND, (ops_u[0], OPS_BOUND)
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
  st.update({"n": n, "pairs": pairs, "roots": n - 2 * pairs, "demoted": demoted, "pair_fraction": round(2 * pairs / n, 5),
         "max_depth": max(depth), "max_tree_size": max(tree.values()), "tree_cap": max_tree, "max_stack": max_stack,
         "augmenting_paths": bool(augment), "aug_cand_cap": cand_cap, "aug_len5_cap": ccy})
  return ops, st

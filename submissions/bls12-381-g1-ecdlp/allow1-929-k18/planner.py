"""Planner: greedy pairing plus a bounded augmenting-path phase (proof.md sections 3 and 16). Stdlib only."""
REL32_NONE = -(1 << 31)          # "no pointer" sentinel in 32-bit relative fields
REL32_MAX = (1 << 31) - 1
def rel32(base, target):
    """Absolute index as rel32 from base."""
    d = target - base
    assert REL32_NONE < d <= REL32_MAX, "planner offset exceeds int32"
    return d
def plan(L, GMAX, max_depth=1 << 14, center_priority=False, max_tree=1 << 20, augment=1, cand_cap=32):
    """Greedy pairing plan with rel32 pointers, then one bounded alternating-path augmentation pass (proof.md s16)."""
    n = len(L)
    pos = {v: i for i, v in enumerate(L)}
    adj = [[] for _ in range(n)]
    n_edges = 0
    for i in range(n):
        vi = L[i]
        lim = vi + 2 * GMAX
        j = i + 1
        while j < n and L[j] <= lim:
            vj = L[j]
            if (vi + vj) & 1 == 0:
                k = pos.get((vi + vj) >> 1)
                if k is not None:
                    adj[i].append((j, k))
                    adj[j].append((i, k))
                    n_edges += 1
            j += 1
    parent = [REL32_NONE] * n
    partner = [REL32_NONE] * n
    head = [REL32_NONE] * n
    nextsib = [REL32_NONE] * n
    depth = [0] * n
    pair_slot = [-1] * n
    height = [0] * n
    stats = {"cycle_rejections": 0, "depth_rejections": 0, "demoted": 0, "size_rejections": 0,
             "depth_increments": 0, "depth_decrements": 0, "height_increments": 0}
    size = [1] * n          # subtree size at every node (the root entry is the tree size); kept by ancestor walks
    def add_size(x, delta):
        while x != -1:
            size[x] += delta
            o = parent[x]
            x = -1 if o == REL32_NONE else x + o
    def in_subtree(x, roots):
        while x != -1:
            if x in roots: return True
            o = parent[x]
            x = -1 if o == REL32_NONE else x + o
        return False
    def parent_of(x):
        o = parent[x]
        return -1 if o == REL32_NONE else x + o
    def partner_of(x):
        o = partner[x]
        return -1 if o == REL32_NONE else x + o
    def first_child(x):
        o = head[x]
        return -1 if o == REL32_NONE else x + o
    def next_pair(x):
        o = nextsib[x]
        return -1 if o == REL32_NONE else x + o
    def set_next(x, target):
        nextsib[x] = REL32_NONE if target < 0 else rel32(x, target)
    def root_of(v):
        while parent[v] != REL32_NONE:
            v += parent[v]
        return v
    def set_depth(v):
        stack = [v]
        while stack:
            x = stack.pop()
            px = parent[x]
            p = -1 if px == REL32_NONE else x + px
            next_depth = depth[p] + 1 if p != -1 else 0
            delta = next_depth - depth[x]
            stats["depth_increments"] += max(0, delta)
            stats["depth_decrements"] += max(0, -delta)
            depth[x] = next_depth
            a = first_child(x)
            while a != -1:
                b = a + partner[a]
                stack.append(a)
                stack.append(b)
                a = next_pair(a)
    def creates_cycle(i, j, k):
        x = k
        while x != -1:
            if x == i or x == j:
                return True
            o = parent[x]
            x = -1 if o == REL32_NONE else x + o
        return False
    def try_match(i, j, k, force=False):
        if parent[i] != REL32_NONE or parent[j] != REL32_NONE:
            return False
        if not force:
            if depth[k] + 1 + max(height[i], height[j]) > max_depth:
                stats["depth_rejections"] += 1
                return False
            if creates_cycle(i, j, k):
                stats["cycle_rejections"] += 1
                return False
            rk = root_of(k)
            if size[rk] + size[i] + size[j] > max_tree:
                stats["size_rejections"] += 1
                return False
        add_size(k, size[i] + size[j])
        # append pair to k's child list, preserving arrival order
        deg, last = 0, -1
        a = first_child(k)
        while a != -1:
            deg += 1
            last = a
            a = next_pair(a)
        lo, hi = (i, j) if i < j else (j, i)
        parent[lo] = rel32(lo, k)
        parent[hi] = rel32(hi, k)
        partner[lo] = rel32(lo, hi)
        partner[hi] = rel32(hi, lo)
        pair_slot[lo] = pair_slot[hi] = deg
        nextsib[lo] = REL32_NONE
        if last == -1:
            head[k] = rel32(k, lo)
        else:
            nextsib[last] = rel32(last, lo)
        set_depth(i)
        set_depth(j)
        h, x = 1 + max(height[i], height[j]), k
        while x != -1 and height[x] < h:
            stats["height_increments"] += h - height[x]
            height[x] = h
            h += 1
            o = parent[x]
            x = -1 if o == REL32_NONE else x + o
        return True
    def unmatch(u):
        k, b = u + parent[u], u + partner[u]
        add_size(k, -(size[u] + size[b]))
        lo, hi = (u, b) if u < b else (b, u)
        position = pair_slot[u]
        prev, cur, idx = -1, first_child(k), 0
        while cur != -1 and idx < position:
            prev, cur, idx = cur, next_pair(cur), idx + 1
        assert idx == position and cur == lo and partner_of(cur) == hi
        last, lastprev = cur, prev
        nxt = next_pair(cur)
        while nxt != -1:
            lastprev, last, nxt = last, nxt, next_pair(nxt)
        if last == cur:
            if prev == -1:
                head[k] = REL32_NONE
            else:
                set_next(prev, -1)
        else:
            after = next_pair(cur)
            set_next(lastprev, -1)
            if after == last:
                set_next(last, -1)
            else:
                set_next(last, after)
            if prev == -1:
                head[k] = rel32(k, last)
            else:
                set_next(prev, last)
            pair_slot[last] = pair_slot[partner_of(last)] = position
        parent[lo] = parent[hi] = REL32_NONE
        partner[lo] = partner[hi] = REL32_NONE
        pair_slot[lo] = pair_slot[hi] = -1
        set_depth(lo)
        set_depth(hi)
    if center_priority:
        by_center = [[] for _ in range(n)]
        for i, edges in enumerate(adj):
            for j, k in edges:
                if i < j:
                    by_center[k].append((i, j))
        for k in sorted(range(n), key=lambda k: -len(by_center[k])):
            for i, j in by_center[k]:
                if parent[i] == REL32_NONE and parent[j] == REL32_NONE:
                    try_match(i, j, k)
    for v in sorted(range(n), key=lambda v: len(adj[v])):
        if parent[v] != REL32_NONE:
            continue
        for (w, k) in adj[v]:
            if parent[w] == REL32_NONE and try_match(v, w, k):
                break
    # s16 augmentation (one pass): for each element u still unmatched after the greedy pass, in increasing degree,
    # take the first cand_cap partners of u in list order and try a direct pair with a free partner; failing that,
    # for each matched partner w with pair (w, w') under kw, take the first cand_cap partners v of w' and try the
    # length-3 alternating path: dissolve (w, w'), form (u, w) and (w', v). A candidate is checked conservatively on
    # the current forest (feasible) and committed through the same unmatch/try_match operations; a rejected candidate
    # mutates nothing. The ledger charges the pass for every element of the class (no count assumption).
    stats.update({"aug_len3": 0, "aug_direct": 0, "aug_cands": 0, "aug_unmatched_in": 0})
    def feasible(new_pairs, detached):
        moved = set()
        for (i, j, k) in new_pairs: moved.add(i); moved.add(j)
        for w in detached: moved.add(w); moved.add(w + partner[w])
        gain = {}
        for (i, j, k) in new_pairs:
            if in_subtree(k, moved): return False                                   # acyclicity and depth independence
            if depth[k] + 1 + max(height[i], height[j]) > max_depth: return False  # depth cap (heights never lowered)
            rk = root_of(k); gain[rk] = gain.get(rk, 0) + size[i] + size[j]
        return all(size[rk] + g <= max_tree for rk, g in gain.items())            # tree cap, detachment ignored
    def commit(new_pairs, detached):
        for w in detached: unmatch(w)
        for (i, j, k) in new_pairs:
            assert try_match(i, j, k), "committed pair failed the live checks"
    for _round in range(augment):
        free_in = [v for v in range(n) if parent[v] == REL32_NONE]
        stats["aug_unmatched_in"] = max(stats["aug_unmatched_in"], len(free_in))
        for u in sorted(free_in, key=lambda v: len(adj[v])):
            if parent[u] != REL32_NONE: continue
            done = False
            for (w, k) in adj[u][:cand_cap]:
                if parent[w] == REL32_NONE and feasible([(u, w, k)], []):
                    commit([(u, w, k)], []); stats["aug_direct"] += 1; done = True; break
            if not done:
                for (w, k) in adj[u][:cand_cap]:
                    if parent[w] == REL32_NONE: continue
                    wp = w + partner[w]
                    for (v, k2) in adj[wp][:cand_cap]:
                        if v == u or v == w or parent[v] != REL32_NONE: continue
                        stats["aug_cands"] += 1
                        if feasible([(u, w, k), (wp, v, k2)], [w]):
                            commit([(u, w, k), (wp, v, k2)], [w]); stats["aug_len3"] += 1; done = True; break
                    if done: break
    if parent[0] != REL32_NONE:
        unmatch(0)
        stats["demoted"] += 1
    last = L[0]
    for v in range(1, n):
        if L[v] - last > GMAX:
            u = v - 1
            if parent[u] != REL32_NONE:
                unmatch(u)
                stats["demoted"] += 1
            last = L[u]
            assert L[v] - last <= GMAX
        if parent[v] == REL32_NONE:
            last = L[v]
    emitted = [False] * n
    ops = []
    max_stack = n_pairs = n_single = 0
    last_root = None
    for v in range(n):
        if emitted[v] or parent[v] != REL32_NONE:
            continue
        if last_root is None:
            ops.append((0, v))
        else:
            ops.append((1, v, L[v] - L[last_root]))
            n_single += 1
        emitted[v] = True
        last_root = v
        stack = [(v, 0)]
        while stack:
            max_stack = max(max_stack, len(stack))
            x, idx = stack[-1]
            a, seen = first_child(x), 0
            while a != -1 and seen < idx:
                a = next_pair(a)
                seen += 1
            if a != -1:
                stack[-1] = (x, idx + 1)
                b = a + partner[a]
                g = L[x] - L[a]
                assert L[b] - L[x] == g and 0 < g <= GMAX
                ops.append((2, x, a, b, g))
                n_pairs += 1
                emitted[a] = emitted[b] = True
                stack.append((a, 0))
                stack.append((b, 0))
            else:
                stack.pop()
    assert all(emitted), "dependency cycle remains"
    if not augment:
        assert stats["depth_increments"] <= n * max_depth
        assert stats["depth_decrements"] <= stats["depth_increments"]
        assert stats["height_increments"] <= n * max_depth
    for v in range(n):   # subtree sizes consistent at every root
        if parent[v] == REL32_NONE:
            st, tot = [v], 0
            while st:
                y = st.pop(); tot += 1; a = first_child(y)
                while a != -1:
                    st.append(a); st.append(a + partner[a]); a = next_pair(a)
            assert size[v] == tot
    assert max(depth) <= max_depth and max_stack <= 2 * max_depth + 1
    assert n_single + 1 + 2 * n_pairs == n
    tree_sizes = {}
    for v in range(n):
        rv = root_of(v)
        tree_sizes[rv] = tree_sizes.get(rv, 0) + 1
    assert max(tree_sizes.values()) <= max_tree
    degrees = []
    for x in range(n):
        d, a = 0, first_child(x)
        while a != -1:
            d += 1
            a = next_pair(a)
        degrees.append(d)
    stats.update({"n": n, "edges": n_edges, "pairs": n_pairs, "roots": n_single + 1,
                  "max_tree_size": max(tree_sizes.values()), "tree_cap": max_tree,
                  "pair_fraction": round(2 * n_pairs / n, 5), "max_depth": max(depth),
                  "max_stack": max_stack, "center_priority": center_priority,
                  "augmenting_paths": bool(augment), "aug_cand_cap": cand_cap,
                  "used_centers": sum(1 for d in degrees if d),
                  "center_fraction": sum(1 for d in degrees if d) / n,
                  "max_pairs_per_center": max(degrees),
                  "depth_update_bound_per_vertex": 2 * max_depth})
    return ops, stats

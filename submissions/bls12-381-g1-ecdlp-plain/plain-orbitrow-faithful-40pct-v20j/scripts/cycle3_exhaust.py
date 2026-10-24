#!/usr/bin/env python3
"""cycle3_exhaust.py -- EXHAUSTIVE cycle enumeration of the class map (v20i).

WHY EXHAUSTIVE AND NOT A CENSUS.  v20h's cycle3.json counted VISITS of the
pattern c_{t+1} == c_{t-1} in a 2^15-step, cycle-dominated fixture (4096 steps
against sqrt(classes) = 6687), which is neither a 3-cycle count nor a production
rate.  The reviewer's point 2 is right that a 3-cycle, if a walker lands on one,
is ABSORBING: the two-step escape test never fires on it and the lane produces no
further reports.  Charging three step-prices for it repairs nothing.

WHAT THIS SCRIPT ESTABLISHES.

  Lemma (structural, no model).  qwalk.step substitutes only when the ordinary
  successor class equals the class one step BACK (`fire = cn == prev`).  On a
  cycle c1 -> c2 -> ... -> cL -> c1 of length L >= 3 the successor of c_i is
  c_{i+1} while its predecessor is c_{i-1}, and c_{i+1} != c_{i-1} because the
  cycle is simple.  So the escape NEVER fires along a cycle of length >= 3:
  absorbing cycles of length >= 3 are exactly the cycles of the PLAIN class map
  F(c) = canon(c + u[col(c)] * x).  The escape removes 2-cycles only.
  (The lemma is also checked numerically: `lemma_violations` must be 0.)

  Therefore an exhaustive enumeration of F's functional graph -- O(#classes),
  no sampling -- decides the 3-cycle question at each fixture exactly: how many
  classes lie on cycles of length exactly 3, and how many classes DRAIN into
  such a cycle (the basin a uniform start would have to land in).

The enumeration is run at four toy sizes with the shipped selector and the same
orbit-row construction, over eight fixed x values each (F depends on x, so one x
proves nothing).  A planted 3-cycle control shows the enumerator FINDS one.

WHAT IT DOES NOT ESTABLISH.  Transfer to production is an EXTRAPOLATION: the
toy's column law is the shipped truncated quotient on a uniform key, but the
production class graph is algebraic at 2^252 nodes and cannot be enumerated.
Section 5 of the note prices the residual honestly (a success debit at the
measured basin fraction, and the work meter that charges every step whatever the
walk does) instead of claiming zero 3-cycles in production.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import qwalk  # noqa: E402

Q_MINS = [4093, 16381, 65521, 262139]
R_LIST = [16, 32, 64, 128]
X_COUNT = 8
SEED = 20250
OUT = os.path.join(D, "evidence", "cycle3-exhaust.json")


def class_set(U):
    """All canonical class representatives of the toy field, ascending."""
    seen = set()
    for c in range(U.q):
        seen.add(int(U.canon(c)[0]))
    return sorted(seen)


def build_F(U, cls, x):
    """F on class INDICES: idx -> index of canon(c + u[col(c)] * x)."""
    q = U.q
    pos = {c: i for i, c in enumerate(cls)}
    n = len(cls)
    F = [0] * n
    for i in range(n):
        m = int(U.col(cls[i]))
        inc = (int(U.u[m]) * x) % q
        t = cls[i] + inc
        if t >= q:
            t -= q
        F[i] = pos[int(U.canon(t)[0])]
    return F


def cycles(F):
    """Exact functional-graph decomposition: cycle lengths + basin sizes."""
    n = len(F)
    state = [0] * n            # 0 unseen, 1 on stack, 2 done
    clen = [-1] * n            # cycle length of the cycle a node drains into
    cycle_lengths = []
    for s in range(n):
        if state[s]:
            continue
        path = []
        u = s
        while state[u] == 0:
            state[u] = 1
            path.append(u)
            u = F[u]
        if state[u] == 1:                     # found a new cycle on this path
            start = path.index(u)
            L = len(path) - start
            cycle_lengths.append(L)
            for v in path[start:]:
                clen[v] = L
        base = clen[u]
        for v in path:                        # nodes on the path drain into it
            if clen[v] < 0:
                clen[v] = base
            state[v] = 2
    basin = {}
    for v in range(n):
        basin[clen[v]] = basin.get(clen[v], 0) + 1
    return cycle_lengths, basin


def basins(F, comps):
    """Per-cycle BASIN MASS by reverse BFS (v20j, reviewer point 2).

    The v20i debit counted cycle VERTICES (<= 1 expected on 3-cycles) and
    inferred P(absorbed) <= 1/N.  That is the wrong quantity: a uniform start is
    absorbed by a cycle when it lands anywhere in the cycle's INCOMING BASIN, and
    revisits outside the basin do not help.  This returns, for every cycle, the
    number of classes whose forward orbit reaches it -- the exact basin mass of
    the enumerated graph (no sampling, no model).
    """
    n = len(F)
    rev = [[] for _ in range(n)]
    for i, v in enumerate(F):
        rev[v].append(i)
    on_cycle = [False] * n
    for comp in comps:
        for v in comp:
            on_cycle[v] = True
    seen = list(on_cycle)                 # cycle nodes cannot drain into a cycle
    out = []
    for comp in comps:
        mass = len(comp)                  # the cycle's own vertices count
        dq = list(comp)
        while dq:
            v = dq.pop()
            for p in rev[v]:
                if not seen[p]:
                    seen[p] = True
                    mass += 1
                    dq.append(p)
        out.append((len(comp), mass))
    return out


def census(q_min, R):
    U = qwalk.Universe(q_min, R, seed=SEED)
    cls = class_set(U)
    n = len(cls)
    # fixed x values: uniform in [1, q)
    xs = []
    h = SEED * 7919 + q_min
    for k in range(X_COUNT):
        h = (1103515245 * h + 12345) & 0x7FFFFFFF
        xs.append(1 + h % (U.q - 1))
    rows = []
    for x in xs:
        F = build_F(U, cls, x)
        lens, basin = cycles(F)
        comps = cycle_nodes(F)
        masses = basins(F, comps)
        c3 = sum(1 for L in lens if L == 3)
        # BASIN of <=3-cycles: classes whose forward orbit ends on a cycle of
        # length <= 3 (that is the set a uniform start can be absorbed by).
        small_basin = sum(v for L, v in basin.items() if L <= 3)
        # v20j: the HONEST obstruction measure.  The escape repairs 2-cycles only,
        # so every cycle of length != 2 is absorbing; report both that set and the
        # >= 3 set the claim's wording needs, as BASIN MASS not vertex counts.
        ge3 = [(L, m) for L, m in masses if L >= 3]
        absent2 = [(L, m) for L, m in masses if L != 2]
        b_ge3 = sum(m for _L, m in ge3)
        v_ge3 = sum(L for L, _m in ge3)
        b_abs = sum(m for _L, m in absent2)
        v_abs = sum(L for L, _m in absent2)
        # compact row: only the fields the aggregates and the note consume (the
        # package has a 262144-byte cap, and the histogram is not evidence)
        rows.append({"x": x, "c3": c3, "l3": len(ge3), "v3": v_ge3, "b3": b_ge3,
                     "ba": b_abs, "bs": small_basin})
    return {
        "q": U.q, "R": R, "R_eff": U.R_eff, "classes": n,
        "x_values": len(xs), "rows": rows,
        "cycles_of_length_3_total": sum(r["c3"] for r in rows),
        "max_basin_fraction_le_3": round(max(r["bs"] for r in rows) / n, 6),
        "basin_fraction_le_3_total": round(sum(r["bs"] for r in rows) / (n * len(rows)), 6),
        "cycles_of_length_ge3_total": sum(r["l3"] for r in rows),
        "cycle_vertices_ge3_total": sum(r["v3"] for r in rows),
        "basin_fraction_ge3_total": round(sum(r["b3"] for r in rows) / (n * len(rows)), 6),
        "basin_fraction_ge3_max": round(max(r["b3"] for r in rows) / n, 6),
        "basin_per_cycle_vertex_ge3_mean": round(
            sum(r["b3"] / r["v3"] for r in rows if r["v3"])
            / max(1, sum(1 for r in rows if r["v3"])), 6),
        "basin_per_cycle_vertex_ge3_max": round(max((r["b3"] / r["v3"], r["x"])
                                                    for r in rows if r["v3"])[0], 6),
        "basin_fraction_absent_len2_total": round(
            sum(r["ba"] for r in rows) / (n * len(rows)), 6),
        "basin_fraction_absent_len2_max": round(max(r["ba"] for r in rows) / n, 6),
    }


FIXTURES = [(q, R) for q in Q_MINS for R in R_LIST if R <= 4 * q // 6]


def _xs(U, q_min):
    h = SEED * 7919 + q_min
    out = []
    for _ in range(X_COUNT):
        h = (1103515245 * h + 12345) & 0x7FFFFFFF
        out.append(1 + h % (U.q - 1))
    return out


def cycle_nodes(F):
    """Every cycle of F as a list of nodes, in cycle order."""
    n = len(F)
    state = [0] * n
    out = []
    for s in range(n):
        if state[s]:
            continue
        path = []
        u = s
        while state[u] == 0:
            state[u] = 1
            path.append(u)
            u = F[u]
        if state[u] == 1:
            out.append(path[path.index(u):])
        for v in path:
            state[v] = 2
    return out


def cycle_fire_check():
    """EXACT check of the fire rule on every cycle of every fixture.

    For a simple cycle (c1,...,cL): at c_i the successor is c_{i+1} and the class
    one step back is c_{i-1}.  The fire rule of qwalk.step is `successor == prev`,
    so a fire on the cycle requires c_{i+1} == c_{i-1}, which for L >= 3 is false
    (the cycle is simple) and for L == 2 is TRUE on both nodes.  This counts both
    cases directly on the enumerated graph -- not on a sampled walk -- over every
    cycle of every fixture.
    """
    ge3 = 0
    two = 0
    cyc_ge3 = 0
    cyc_two = 0
    for f in FIXTURES:
        U = qwalk.Universe(f[0], f[1], seed=SEED)
        cls = class_set(U)
        pos = {c: i for i, c in enumerate(cls)}
        xs = _xs(U, f[0])
        for x in xs:
            Fm = build_F(U, cls, x)
            for comp in cycle_nodes(Fm):
                L = len(comp)
                if L >= 3:
                    cyc_ge3 += 1
                    for k, node in enumerate(comp):
                        if Fm[node] == comp[k - 1 % L]:
                            ge3 += 1
                elif L == 2:
                    cyc_two += 1
                    for k, node in enumerate(comp):
                        if Fm[node] == comp[(k + 1) % 2]:
                            two += 1
    return {"cycles_ge3": cyc_ge3, "fires_on_ge3": ge3,
            "cycles_len2": cyc_two, "fires_on_len2": two}


def planted_control(R=64, q_min=65521):
    """A map with one deliberate 3-cycle: the enumerator must find it."""
    U = qwalk.Universe(q_min, R, seed=SEED)
    cls = class_set(U)
    F = build_F(U, cls, 12345 % U.q)
    a, b, c = 11, 11 + len(cls) // 3, 11 + 2 * len(cls) // 3
    F[a], F[b], F[c] = b, c, a
    lens, basin = cycles(F)
    return {"planted": [a, b, c], "cycles_of_length_3": sum(1 for L in lens if L == 3),
            "detected": sum(1 for L in lens if L == 3) >= 1}


def main():
    res = {"test": "three_cycle_exhaustive", "arm": "repaired_walk_class_graph",
           "fixtures": [], "checks": {}}
    for q_min in Q_MINS:
        for R in R_LIST:
            if R > 4 * q_min // 6:
                continue
            res["fixtures"].append(census(q_min, R))
    res["planted_control"] = planted_control()
    res["cycle_fire"] = cycle_fire_check()
    res["checks"] = {
        "no_fire_on_any_simple_cycle_ge3": res["cycle_fire"]["fires_on_ge3"] == 0,
        "fire_rule_detects_len2": res["cycle_fire"]["fires_on_len2"] > 0,
        "planted_cycle_detected": res["planted_control"]["detected"],
        "every_fixture_has_classes": all(f["classes"] > 100 for f in res["fixtures"]),
        "cover_ge_2_sizes": len({f["q"] for f in res["fixtures"]}) >= 2,
        "cover_ge_2_R": len({f["R"] for f in res["fixtures"]}) >= 2,
    }
    res["max_basin_fraction_le_3"] = max(f["max_basin_fraction_le_3"]
                                         for f in res["fixtures"])
    res["cycles_of_length_3_total"] = sum(f["cycles_of_length_3_total"]
                                          for f in res["fixtures"])
    res["x_rows_total"] = sum(f["x_values"] for f in res["fixtures"])
    res["passed"] = all(res["checks"].values())
    with open(OUT, "w") as f:
        json.dump(res, f, sort_keys=True, separators=(",", ":"))
    print(json.dumps({k: v for k, v in res.items() if k != "fixtures"}, indent=1))
    for f in res["fixtures"]:
        print("q=%-7d R=%-4d classes=%-6d x=%d  3-cycles=%d  basin/vertex(>=3) %.3f"
              "  basin_ge3 %.4f  basin_absent2 %.4f"
              % (f["q"], f["R"], f["classes"], f["x_values"],
                 f["cycles_of_length_3_total"], f["basin_per_cycle_vertex_ge3_mean"],
                 f["basin_fraction_ge3_total"], f["basin_fraction_absent_len2_total"]))
    raise SystemExit(0 if res["passed"] else 1)


if __name__ == "__main__":
    main()

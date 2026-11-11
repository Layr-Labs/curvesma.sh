import importlib.util
import json
import math
import random
import re
import sys
from decimal import Decimal, ROUND_CEILING, getcontext
from pathlib import Path

base = Path("/snapshot/files")
if not base.is_dir():
    base = Path(__file__).parent
spec = importlib.util.spec_from_file_location("tt_ref", base / "reference.py")
ref = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ref)
request = json.load(sys.stdin)
seed, trials = request["seed"], request["trials"]
coins = random.Random(seed)
profile = json.loads((base / "target-profile.json").read_text())
p = int(profile["field_modulus"], 16)
r = int(profile["subgroup_order"], 16)
g = tuple(int(v, 16) for v in profile["generator_g"]["coefficients"])
G = ref.BN254GT(p, r, g)
M64 = ref.M64

# ---- 1. Group and tower checks at full width --------------------------------
assert G.power(g, r) == G.identity and g != G.identity

def f2mul(a, b):
    v0, v1 = a[0] * b[0] % p, a[1] * b[1] % p
    return ((v0 - v1) % p, ((a[0] + a[1]) * (b[0] + b[1]) - v0 - v1) % p)
def f2add(a, b): return ((a[0] + b[0]) % p, (a[1] + b[1]) % p)
def f2sub(a, b): return ((a[0] - b[0]) % p, (a[1] - b[1]) % p)
def f2xi(a): return ((9 * a[0] - a[1]) % p, (9 * a[1] + a[0]) % p)
def f6mul(a, b):
    v0, v1, v2 = f2mul(a[0], b[0]), f2mul(a[1], b[1]), f2mul(a[2], b[2])
    c0 = f2add(v0, f2xi(f2sub(f2sub(f2mul(f2add(a[1], a[2]), f2add(b[1], b[2])), v1), v2)))
    c1 = f2add(f2sub(f2sub(f2mul(f2add(a[0], a[1]), f2add(b[0], b[1])), v0), v1), f2xi(v2))
    c2 = f2add(f2sub(f2sub(f2mul(f2add(a[0], a[2]), f2add(b[0], b[2])), v0), v2), v1)
    return (c0, c1, c2)
def f6add(a, b): return tuple(f2add(x, y) for x, y in zip(a, b))
def f6sub(a, b): return tuple(f2sub(x, y) for x, y in zip(a, b))
def f12mul(a, b):                       # mirror of tower.c fp12_mul
    v0, v1 = f6mul(a[0], b[0]), f6mul(a[1], b[1])
    c1 = f6sub(f6sub(f6mul(f6add(a[0], a[1]), f6add(b[0], b[1])), v0), v1)
    c0 = (f2add(v0[0], f2xi(v1[2])), f2add(v0[1], v1[0]), f2add(v0[2], v1[1]))
    return (c0, c1)
def to_tower(Y):
    c = G.coords(Y)
    return tuple(tuple((c[2 * (3 * h + m)], c[2 * (3 * h + m) + 1]) for m in range(3)) for h in range(2))
def flat(t): return [x for h in t for m in h for x in m]

def rand_gt():
    return G.power(g, coins.randrange(1, r))
def rand_fp12():
    return tuple(coins.randrange(p) for _ in range(12))

tower_checks = row_checks = 0
for _ in range(6):
    Y, T = rand_fp12(), rand_fp12()
    assert flat(f12mul(to_tower(Y), to_tower(T))) == G.coords(G.mul(Y, T))
    tower_checks += 1
for _ in range(40):
    Y, T = (rand_gt(), rand_gt()) if _ < 3 else (rand_fp12(), rand_fp12())
    assert G.coord0(G.mul(Y, T)) == sum(a * b for a, b in zip(G.coords(Y), G.row(T))) % p
    row_checks += 1

# ---- 2. The tag certificate: 2^64*t/p lies in [A, A+96) mod 2^64 -------------
def F_of(T):
    rinv = pow(1 << 256, -1, p)
    out = []
    for vk in G.row(T):
        for i in range(4):
            f = ((vk << (64 * i)) * rinv % p << 128) // p
            out.append((f & M64, f >> 64))
    return out
bound_checks, worst = 0, 0
edge = [G.identity, tuple([p - 1] * 12), (0,) * 11 + (p - 1,), (p - 1,) + (0,) * 11]
for i in range(260):
    if i < len(edge) ** 2:
        Y, T = edge[i // len(edge)], edge[i % len(edge)]
    elif i < 22:
        Y, T = rand_gt(), rand_gt()
    else:
        Y, T = rand_fp12(), rand_fp12()
    A = ref.tag_sum(ref.mont_limbs(G, Y), F_of(T))
    t = G.coord0(G.mul(Y, T))
    D = ((t << 64) - A * p) % (p << 64)          # (2^64 t/p - A) mod 2^64, times p
    assert D < 96 * p
    worst = max(worst, D // p)
    bound_checks += 1

# ---- 3. Full-width walk through the real solve loop --------------------------
def traced_run(cbits, steps_wanted, run_seed):
    rng = random.Random(run_seed)
    x = rng.randrange(1, r)
    z = G.power(g, x)
    log = {"t": [], "s": [], "prev": None, "label_checks": 0, "n": 0}
    def trace(table, Y0, e, s, dp, a, b):
        Ycur = G.mul(Y0, e.T)
        t = G.coord0(Ycur)
        assert s == ref.exact_selector(G, t) and dp == ((t << 80) < p)
        if log["prev"] is not None:               # one multiplier per step
            Yp, sp = log["prev"]
            assert Ycur == G.mul(Yp, table.mults[sp][0])
        if log["n"] % 40 == 0:                    # Y = g^a z^b with the entry's increments
            aa, bb = (a + e.alpha) % r, (b + e.beta) % r
            assert Ycur == G.mul(G.power(g, aa), G.power(z, bb))
            log["label_checks"] += 1
        log["prev"] = (Ycur, s)
        log["t"].append(t); log["s"].append(s); log["n"] += 1
    got, steps, outcome, stats = ref.solve(G, z, 5, 80, cbits, steps_wanted, 1, 1 << 30, rng, trace)
    assert outcome == "cap" and steps == steps_wanted
    return log, stats
log3, stats3 = traced_run(3, 120, seed + 1)        # production certificate width
log54, stats54 = traced_run(54, 120, seed + 1)     # forces many exact paths
assert log3["t"] == log54["t"] and log3["s"] == log54["s"]
assert stats3["exact_paths"] == 0 and stats54["exact_paths"] > 0 and stats3["rebases"] >= 15

# ---- 3b. Scalar output at full width (gt-dlp track: return x itself) -------
x_id, _, outcome_id, _ = ref.solve(G, G.identity, 5, 80, 3, 1, 1, 1, None)
assert outcome_id == "identity" and x_id == 0 and G.power(g, x_id) == G.identity
def encode_x(x):
    assert 0 <= x < r
    return format(x, "064x")
enc_checks = 0
for xv in (0, 1, r - 1, random.Random(seed + 2).randrange(r)):
    h = encode_x(xv)
    assert len(h) == 64 and h == h.lower() and int(h, 16) == xv
    enc_checks += 1

# ---- 4. Complete algorithm on toy prime-field subgroups ----------------------
def isprime(n):
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for q in small:
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in small:                               # deterministic below 3.3e24
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True
def toy_group(nbits):
    while True:
        n = coins.getrandbits(nbits) | (1 << (nbits - 1)) | 1
        if isprime(n):
            break
    c = (1 << 60) // n
    while not isprime(c * n + 1):
        c += 1
    q = c * n + 1
    while True:
        gg = pow(coins.randrange(2, q), c, q)
        if gg != 1:
            return ref.ToyField(q, n, gg)
toy = []
for nbits, dbits, share in ((16, 3, 1.0), (20, 4, 1.0), (24, 5, 0.5)):
    n_trials = max(16, int(trials * share))
    counts = {k: 0 for k in ("success", "identity", "degenerate", "cap", "walk_cap")}
    for _ in range(n_trials):
        T = toy_group(nbits)                      # fresh group, target and coins each trial
        secret = coins.randrange(T.r)
        got, n, outcome, st = ref.solve(T, T.power(T.g, secret), 6, dbits, dbits,
                                        ref.step_cap(T.r), 1 << 20, 64 << dbits, coins)
        counts[outcome] += 1
        if got is not None:
            assert got == secret
    toy.append({"group_order_bits": nbits, "field_bits": 61, "distinguished_bits": dbits,
                "trials": n_trials, "outcomes_at_cap": counts,
                "success_fraction_at_cap": (counts["success"] + counts["identity"]) / n_trials})

# ---- 5. Instruction recount from the supplied listings -----------------------
def parse(name):
    fn, cur = {}, None
    for line in (base / name).read_text().splitlines():
        head = re.match(r"^([0-9a-f]+) <_?(\w+)>:", line)
        if head:
            cur = head.group(2); fn[cur] = []
            continue
        ins = re.match(r"^\s*([0-9a-f]+):\s+(\S+)\s*(.*)$", line)
        if ins and cur:
            fn[cur].append((int(ins.group(1), 16), ins.group(2), ins.group(3)))
    return fn
SIMD = re.compile(r"\b[vqdsz]\d+\b")
tower = parse("tower.disasm.txt")
memo = {}
def executed(name):
    if name not in memo:
        total = 0
        for _, op, arg in tower[name]:
            total += 1
            callee = re.search(r"<_?(\w+)>", arg)
            if op == "bl" or (op == "b" and callee and callee.group(1) != name):
                total += executed(callee.group(1))
            else:
                assert not (op == "b" or op.startswith("b.") or op in ("cbz", "cbnz", "tbz", "tbnz", "br", "blr"))
                assert not SIMD.search(arg)
        memo[name] = total
    return memo[name]
tower_counts = {n: executed(n) for n in tower if n != "main"}
run = parse("ttstep.disasm.txt")["tt_run"]
addr = [a for a, _, _ in run]
back = [(a, int(re.match(r"0x([0-9a-f]+)", arg).group(1), 16)) for a, op, arg in run
        if op.startswith("b.") and int(re.match(r"0x([0-9a-f]+)", arg).group(1), 16) <= a]
assert len(back) == 1                              # exactly one backward branch: the loop
tail, head = back[0]
body = [(a, op, arg) for a, op, arg in run if head <= a <= tail]
for a, op, arg in body:
    assert not SIMD.search(arg) and op not in ("bl", "blr", "br", "ret", "b")
    if (op.startswith("b.") or op in ("cbz", "cbnz", "tbz", "tbnz")) and a != tail:
        target = int(re.search(r"0x([0-9a-f]+)", arg).group(1), 16)
        assert target > tail                       # the only other branches leave the loop
loop_instructions = len(body)
call_overhead = len(run) - loop_instructions
ops = [op for _, op, _ in body]
kernel = {"loop_instructions": loop_instructions, "call_overhead": call_overhead,
          "umulh": ops.count("umulh"), "madd_or_mul": ops.count("madd") + ops.count("mul"),
          "loads": sum(2 if op == "ldp" else 1 for op in ops if op in ("ldr", "ldp")),
          "stores": sum(1 for op in ops if op in ("str", "stp"))}
assert kernel["umulh"] == 48 and kernel["madd_or_mul"] == 48 and kernel["stores"] == 0

# ---- 6. Ledger -----------------------------------------------------------------
STEP, WINDOW, ELL, FIXED = 273, 24576, 256, 1 << 110
assert loop_instructions <= STEP
assert tower_counts["rebase"] + call_overhead + 128 <= WINDOW
S = ref.step_cap(r)
windows = S // ELL + 1
R_total = S * STEP + windows * WINDOW + FIXED
getcontext().prec = 90
score = (Decimal(R_total) / 20000).ln() / Decimal(2).ln()
score_up = str(score.quantize(Decimal("0.000001"), rounding=ROUND_CEILING))
entries = math.comb(ELL + 8, 8)
table_bytes = entries * 1280
store_bytes = (1 << 52) * 96
print(json.dumps({
    "seed": seed,
    "tower_vs_profile_basis_cases": tower_checks, "row_identity_cases": row_checks,
    "tag_certificate_cases": bound_checks, "tag_certificate_worst_offset_of_96": worst,
    "full_width_walk": {"steps": log3["n"], "selectors": log3["s"][:24],
                        "label_checks": log3["label_checks"], "window_rebases": stats3["rebases"],
                        "exact_paths_cbits3": stats3["exact_paths"],
                        "exact_paths_cbits54_same_walk": stats54["exact_paths"]},
    "scalar_output": {"identity_x": encode_x(x_id), "encoding_checks": enc_checks},
    "toy_runs": toy, "model_success_at_cap": 1 - math.exp(-1409 * 7 / 16000),
    "tower_instruction_counts": tower_counts, "tag_kernel": kernel,
    "step_cap_S": str(S), "windows": str(windows),
    "instructions_per_step": STEP, "instructions_per_window": WINDOW,
    "fixed_instruction_allowance": "2^110", "total_instruction_bound": str(R_total),
    "score_rounded_up": score_up,
    "table_entries": entries, "table_bytes_log2": round(math.log2(table_bytes), 3),
    "peak_bytes_log2_bound": round(math.log2(table_bytes + store_bytes + (1 << 40)), 3),
    "scope": "Full-width BN254 checks of the tower formulas, the tag row identity, the tag "
             "certificate, a 120-step walk through the real solve loop, the identity case and the "
             "scalar encoding; complete-algorithm "
             "runs on toy prime-field subgroups only; static instruction recount of the supplied "
             "listings. No production DLP, profiler measurement or acceptance is claimed.",
}))

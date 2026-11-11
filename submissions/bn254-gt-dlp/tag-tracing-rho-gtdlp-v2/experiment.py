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
kspec = importlib.util.spec_from_file_location("tt_kernels", base / "kernels.py")
kern = importlib.util.module_from_spec(kspec)
kspec.loader.exec_module(kern)
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
def f12mul(a, b):      # formulas of kernels.listing_rebase (b-side sums are the stored ones)
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

# ---- 5. Instruction counts ------------------------------------------------------
# 5a. The step loop: static recount of the clang listing of ttstep.c.
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
run = parse("ttstep.disasm.txt")["tt_run"]
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
assert all(not SIMD.search(arg) for _, _, arg in run)
loop_instructions = len(body)
call_overhead = len(run) - loop_instructions
ops = [op for _, op, _ in body]
kernel = {"loop_instructions": loop_instructions, "call_overhead": call_overhead,
          "umulh": ops.count("umulh"), "madd_or_mul": ops.count("madd") + ops.count("mul"),
          "loads": sum(2 if op == "ldp" else 1 for op in ops if op in ("ldr", "ldp")),
          "stores": sum(1 for op in ops if op in ("str", "stp"))}
assert kernel["umulh"] == 48 and kernel["madd_or_mul"] == 48 and kernel["stores"] == 0
# tt_run reads only F (entry offsets 0..767) and next[s] (768 + 8s) of an entry (pointer in x12),
# so the fields appended after them in v2 (kernels.py layout) do not change the listing
ent_ends = [int(m or "0", 16) + (16 if op == "ldp" else 8) for _, op, arg in body
            for m in re.findall(r"\[x12(?:, #0x([0-9a-f]+))?\]", arg)]
assert max(ent_ends) == 768 and sum("#0x300]" in arg for _, _, arg in body) == 1
assert kern.ENTRY_T == 96 * 8 + 8 * 8

# 5b. The window re-base: generated AArch64 listings executed on the simulator.
prog = kern.program()
labels = {"K": 0x4000, "KR": 0x4040}
names = set(prog) | set(labels)
for L in prog.values():
    for line in L:
        kern.check_encodable(line, names)
listing_text = kern.render(prog)
shipped = (base / "listings.txt").read_text()
assert shipped.rstrip("\n") == listing_text.rstrip("\n"), "listings.txt differs from the generator"
M = kern.Machine(prog, labels)
krng = random.Random(seed ^ 0x6B65726E)            # separate stream: sections 1-4 are unchanged
def limbs4(v): return [(v >> (64 * i)) & M64 for i in range(4)]
def put(mem, a, ws):
    for i, w in enumerate(ws):
        mem[a + 8 * i] = w
def fresh_mem():
    mem = {}
    put(mem, labels["K"], limbs4(p) + [(-pow(p, -1, 1 << 64)) % (1 << 64)])
    put(mem, labels["KR"], limbs4(r) + [(-pow(r, -1, 1 << 64)) % (1 << 64)])
    return mem
SP0 = 0x90000
def call(fn, x0, x1, x2, mem):
    regs = {"x%d" % i: krng.getrandbits(64) for i in range(29)}
    regs.update({"x0": x0, "x1": x1, "x2": x2, "sp": SP0})
    saved = {k: regs[k] for k in kern.CALLEE + ["sp"]}
    before = dict(mem)
    n = M.run(fn, regs, mem)
    assert all(regs[k] == saved[k] for k in saved), "callee-saved register or sp not restored"
    return n, before
def entry_words(T, al, be):
    """T_e, alpha, beta and the stored sums, as the entry holds them (Montgomery limbs)."""
    tl = ref.mont_limbs(G, T)
    v = [sum(tl[4 * j + i] << (64 * i) for i in range(4)) for j in range(12)]
    f2 = [(v[2 * k], v[2 * k + 1]) for k in range(6)]          # memory order a0,a2,a4,a1,a3,a5
    add2 = lambda x, y: ((x[0] + y[0]) % p, (x[1] + y[1]) % p)
    s3 = lambda B: [add2(B[1], B[2]), add2(B[0], B[1]), add2(B[0], B[2])]
    Y0, Y1 = f2[0:3], f2[3:6]
    Ys = [add2(x, y) for x, y in zip(Y0, Y1)]
    pre = [w for e in Ys + s3(Y0) + s3(Y1) + s3(Ys) for c in e for w in limbs4(c)]
    return tl, limbs4(al) + limbs4(be), pre

# leaf tests (edge, random, all-limbs-large; separate and aliased output)
leaf_cases = leaf_fail = 0
Rinv = pow(1 << 256, -1, p)
small = [0, 1, p - 1, (1 << 192) - 1, (1 << 64) - 1, p >> 1]
quads = [(a0, a1, b0, b1) for a0 in small for a1 in small for b0 in (0, p - 1) for b1 in (0, p - 1, p >> 1)]
quads += [tuple(krng.randrange(p) for _ in range(4)) for _ in range(40)]
quads += [tuple(p - 1 - krng.getrandbits(32) for _ in range(4)) for _ in range(20)]
for a0, a1, b0, b1 in quads:
    want = {"f2mul": [(a0 * b0 - a1 * b1) * Rinv % p, (a0 * b1 + a1 * b0) * Rinv % p],
            "f2mulxi": [(9 * a0 - a1) % p, (9 * a1 + a0) % p],
            "fp_add": [(a0 + b0) % p], "fp_sub": [(a0 - b0) % p]}
    for fn, w in want.items():
        for alias in (False, True):
            mem = fresh_mem()
            A, B = 0x10000, 0x20000
            O = A if alias else 0x30000
            put(mem, A, limbs4(a0) + limbs4(a1)); put(mem, B, limbs4(b0) + limbs4(b1))
            n, before = call(fn, O, A, B, mem)
            got = [sum(mem[O + 32 * k + 8 * i] << (64 * i) for i in range(4)) for k in range(len(w))]
            changed = {ad for ad in mem if mem[ad] != before.get(ad)}
            assert all(O <= ad < O + 32 * len(w) or SP0 - 256 <= ad < SP0 for ad in changed)
            leaf_cases += 1; leaf_fail += got != w
for x, y in [(0, 0), (r - 1, r - 1), (r - 1, 1)] + [(krng.randrange(r), krng.randrange(r)) for _ in range(20)]:
    mem = fresh_mem(); put(mem, 0x10000, limbs4(x)); put(mem, 0x20000, limbs4(y))
    call("fr_add", 0x30000, 0x10000, 0x20000, mem)
    leaf_cases += 1; leaf_fail += sum(mem[0x30000 + 8 * i] << (64 * i) for i in range(4)) != (x + y) % r
assert leaf_fail == 0
# exact path: one f2mul by the integer pair (1, 0) takes a_0 out of Montgomery form
Yx = G.power(g, krng.randrange(1, r))
yl = ref.mont_limbs(G, Yx)
mem = fresh_mem(); put(mem, 0x10000, yl[:8]); put(mem, 0x20000, limbs4(1) + limbs4(0))
call("f2mul", 0x30000, 0x10000, 0x20000, mem)
out_mont = [sum(mem[0x30000 + 32 * k + 8 * i] << (64 * i) for i in range(4)) for k in range(2)]
assert out_mont == G.coords(Yx)[:2] and (out_mont[0] - 9 * out_mont[1]) % p == G.coord0(Yx)

# the full re-base on real GT elements, arbitrary Fp12 elements and edge values
STATE, ENTRY = 0x100000, 0x200000
edge12 = [G.identity, tuple([p - 1] * 12), (0,) * 12, (p - 1,) + (0,) * 11]
rebase_cases, rebase_counts, cert_checks = 0, set(), 0
pairs = [(Y, T) for Y in edge12 for T in edge12]
pairs += [(G.power(g, krng.randrange(1, r)), G.power(g, krng.randrange(1, r))) for _ in range(8)]
pairs += [(tuple(krng.randrange(p) for _ in range(12)), tuple(krng.randrange(p) for _ in range(12))) for _ in range(8)]
F_root = F_of(G.identity)
for idx, (Y, T) in enumerate(pairs):
    a, b, al, be = (krng.randrange(r) for _ in range(4))
    if idx < 4:
        a, b, al, be = r - 1, r - 1, r - 1, r - 1
    mem = fresh_mem()
    tl, ab, pre = entry_words(T, al, be)
    put(mem, STATE, ref.mont_limbs(G, Y) + limbs4(a) + limbs4(b))
    put(mem, ENTRY + kern.ENTRY_T, tl + ab); put(mem, ENTRY + kern.ENTRY_PRE, pre)
    OUT = STATE if idx % 2 == 0 else 0x300000                     # in place, or a second buffer
    n, before = call("rebase", OUT, STATE, ENTRY, mem)
    rebase_counts.add(n)
    changed = {ad for ad in mem if mem[ad] != before.get(ad)}
    assert all(OUT <= ad < OUT + 448 or SP0 - kern.FRAME - 256 <= ad < SP0 for ad in changed)
    got = [mem[OUT + 8 * i] for i in range(56)]
    YT = G.mul(Y, T)
    assert got == ref.mont_limbs(G, YT) + limbs4((a + al) % r) + limbs4((b + be) % r)
    # the executed output is exactly what tt_run reads: its tag certifies coordinate 0 of Y*T
    A = ref.tag_sum(got[:48], F_root)
    D = ((G.coord0(YT) << 64) - A * p) % (p << 64)
    assert D < 96 * p
    cert_checks += 1
    rebase_cases += 1
assert len(rebase_counts) == 1                       # straight line: same count on every input
REBASE = rebase_counts.pop()
calls = kern.count_by_function(prog)
assert REBASE == len(prog["rebase"]) + sum(n * len(prog[f]) for f, n in calls.items())
window_end = {"rebase_executed": REBASE, "driver_listing": len(prog["rebase"]),
              "leaf_calls": calls, "leaf_lengths": {f: len(prog[f]) for f in calls},
              "rebase_cases": rebase_cases, "tag_certificate_on_executed_output": cert_checks,
              "leaf_cases": leaf_cases, "leaf_failures": leaf_fail,
              "listing_lines_checked_encodable": sum(len(L) for L in prog.values())}

# ---- 6. Ledger -----------------------------------------------------------------
DRIVER = 128
STEP, WINDOW, ELL, FIXED = 273, REBASE + call_overhead + DRIVER, 256, 1 << 110
assert loop_instructions <= STEP and WINDOW == 11923
S = ref.step_cap(r)
windows = S // ELL + 1
R_total = S * STEP + windows * WINDOW + FIXED
getcontext().prec = 90
score = (Decimal(R_total) / 20000).ln() / Decimal(2).ln()
score_up = str(score.quantize(Decimal("0.000001"), rounding=ROUND_CEILING))
entries = math.comb(ELL + 8, 8)
table_bytes = entries * kern.ENTRY_BYTES
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
    "tag_kernel": kernel, "window_end": window_end,
    "step_cap_S": str(S), "windows": str(windows),
    "instructions_per_step": STEP, "instructions_per_window": WINDOW,
    "fixed_instruction_allowance": "2^110", "total_instruction_bound": str(R_total),
    "score_rounded_up": score_up,
    "table_entries": entries, "table_bytes_log2": round(math.log2(table_bytes), 3),
    "peak_bytes_log2_bound": round(math.log2(table_bytes + store_bytes + (1 << 40)), 3),
    "scope": "Full-width BN254 checks of the tower formulas, the tag row identity, the tag "
             "certificate, a 120-step walk through the real solve loop, the identity case and the "
             "scalar encoding; the window re-base executed instruction by instruction on a "
             "simulator of the AArch64 subset used, on real GT elements, arbitrary Fp12 elements "
             "and edge values; complete-algorithm runs on toy prime-field subgroups only; static "
             "recount of the step-loop listing. No production DLP, profiler measurement or "
             "acceptance is claimed.",
}))

"""Audit of the measured side-output kernel measure/side_lanes.S (proof.md section 1).

(1) Static path enumeration of one lane (.Llane .. cbnz): every path and its instruction count;
    the worst path is the charged per-lane price.
(2) Word-level interpretation with a64.py (inherited interpreter) of whole calls over groups of
    lanes against exact integer field/curve arithmetic: x', y' fully reduced, key = limbs 1..3 of
    the even element of {y', p - y'}; edge values, randoms, genuine curve points, and inputs biased
    to the rare branches (x1 = x2 = 0 for X >= 2p; y2 near p for D >= p). Every instruction of the
    lane body is visit-counted: both successors of every data-dependent conditional branch are asserted reached,
    each call's lane-body count is asserted <= WORST * n + 1 (the final lane's b .Lexit), and the remainder (prologue, epilogue,
    terminator) is asserted identical for every call and <= 512 (the ledger's fixed-call row).
    Input memory and callee-saved registers are asserted unchanged.
No new Callgrind run happens here; measure/measurement.log is the Callgrind record."""
import sys, json, re, random, time
from pathlib import Path
root = Path(__file__).resolve().parent
if (root / 'files').is_dir(): root = root / 'files'
sys.path.insert(0, str(root))
import a64
params = json.loads((root / 'params.json').read_text()); p = int(params['p'], 16)
R = pow(2, 384, p); RI = pow(R, -1, p); MASK = (1 << 64) - 1; BASE = 0x80000
WORST = 971

def lane_paths(text):
    lines = text.split("\n"); st = lines.index(".Llane:"); code, labels = [], {}
    for l in lines[st:]:
        l = l.split("//")[0].rstrip()
        if l.startswith(".Lexit:"): break
        if re.match(r"^\.L\w+:", l): labels[l[:-1]] = len(code); continue
        if l.startswith("\t") and l.strip(): code.append(l.strip())
    paths = []
    def walk(i, n, tr):
        while True:
            ins = code[i]; n += 1; op = ins.split()[0]
            if op == "cbnz": paths.append((n, tr)); return
            m = re.match(r"(b\.\w+|tbz|tbnz|cbz|b)\b.*?(\.L\w+)$", ins)
            if m and op == "b": assert labels[m.group(2)] > i; i = labels[m.group(2)]; continue
            if m: assert labels[m.group(2)] > i; walk(labels[m.group(2)], n, tr + [op + "+"]); tr = tr + [op + "-"]
            i += 1
    walk(0, 0, [])
    return paths

def prepare(text):
    code, data = text.split('.section .rodata', 1)
    vals = []
    for line in data.splitlines():
        if '.quad' in line: vals += [int(x.strip(), 0) for x in line.split('.quad', 1)[1].split(',')]
    assert sum(vals[i] << (64 * i) for i in range(6)) == p and vals[6] == (-pow(p, -1, 1 << 64)) & MASK
    code = code.replace('__SIZEOF_POINTER__', '8')
    code = re.sub(r'(?<![A-Za-z_])(-?\d+)\*(\d+)', lambda m: str(int(m[1]) * int(m[2])), code)
    code, n = re.subn(r'\badrp\s+x0,\.Lsk_consts', 'mov x0,#' + str(BASE), code); assert n == 1
    code, n = re.subn(r'\badd\s+x0,x0,:lo12:\.Lsk_consts', 'add x0,x0,#0', code); assert n == 1
    code = '\n'.join(l for l in code.splitlines() if not l.lstrip().startswith('#'))
    pr, labels = a64.a64_load(code)
    return pr, labels, {BASE + 8 * i: v for i, v in enumerate(vals)}

def put(m, a, v):
    for i in range(6): m[a + 8 * i] = (v >> (64 * i)) & MASK
def get(m, a): return sum(m.get(a + 8 * i, 0) << (64 * i) for i in range(6))

def add(A, B):
    x, y = A; u, v = B
    t = (3 * x * x * pow(2 * y, -1, p) if x == u else (v - y) * pow(u - x, -1, p)) % p
    xx = (t * t - x - u) % p
    return xx, (t * (x - xx) - y) % p

def expect(x1, y1, x2, y2, iv):   # Montgomery residues in, residues out (formula of the kernel)
    mu = (y1 + y2) * iv * RI % p; xx = (mu * mu * RI - x1 - x2) % p
    yy = (mu * ((xx - x2) % p) * RI + y2) % p
    e = yy if yy % 2 == 0 else p - yy
    return xx, yy, [(e >> (64 * j)) & MASK for j in (1, 2, 3)]

def run(cfg):
    t0 = time.time(); rng = random.Random(cfg.get('seed', 1))
    text = (root / 'measure/side_lanes.S').read_text()
    paths = lane_paths(text); worst = max(n for n, _ in paths)
    assert worst == WORST, worst
    prog, labels, cm = prepare(text)
    vals = [0, 1, p - 1, p - 2, R]; cases = []
    for i, a in enumerate(vals):
        for j, b in enumerate(vals):
            for iv in vals: cases.append((vals[(i + j) % 5], a, vals[(2 * i + j) % 5], b or 1, iv))
    cases += [tuple(rng.randrange(1 if k == 3 else 0, p) for k in range(5)) for _ in range(cfg.get('samples', 200))]
    cases += [(0, rng.randrange(p), 0, rng.randrange(1, p), rng.randrange(p)) for _ in range(60)]
    cases += [(rng.randrange(p), rng.randrange(p), rng.randrange(p), p - rng.randrange(1, 1 << 300), rng.randrange(p)) for _ in range(60)]
    G = (int(params['gx'], 16), int(params['gy'], 16)); A, B = G, add(G, G); curve = []
    for _ in range(24):
        A = add(A, G); B = add(B, add(G, G))
        if A[0] == B[0]: continue
        D = add(A, (B[0], -B[1] % p)); curve.append(len(cases))
        cases.append((A[0] * R % p, A[1] * R % p, B[0] * R % p, B[1] * R % p, pow((B[0] - A[0]) % p, -1, p) * R % p))
        assert expect(*cases[-1])[:2] == (D[0] * R % p, D[1] * R % p)
    fixed, lanes, per_call = set(), 0, []
    lo = labels['.Llane'][0]; hi = labels['.Lexit'][0]
    body_visits = {i: 0 for i in range(lo, hi)}
    for begin in range(0, len(cases), 8):
        grp = cases[begin:begin + 8]; n = len(grp); m = dict(cm); want = []
        for i, (x1, y1, x2, y2, iv) in enumerate(grp):
            out, pp, tt, ii = 0x20000 + 256 * i, 0x30000 + 256 * i, 0x40000 + 256 * i, 0x50000 + 256 * i
            for j, v in enumerate((out, pp, tt, ii)): m[0x10000 + 32 * i + 8 * j] = v
            put(m, pp, x1); put(m, pp + 48, y1); put(m, tt, x2); put(m, tt + 48, y2)
            put(m, tt + 96, p - x2); put(m, tt + 144, p - y2); put(m, ii, iv)
            want.append(expect(x1, y1, x2, y2, iv))
        m[0x10000 + 32 * n] = m[0x10000 + 32 * n + 8] = 0         # terminator args[n] = {NULL, NULL}
        regs = {f'x{i}': rng.getrandbits(64) for i in range(31)}; regs.update(x0=0x10000, x1=n, sp=0x7000000)
        mem = dict(m); ro = {}; vis = {i: 0 for i in range(lo, hi)}
        c = a64.a64_run(prog, labels, 'side_lanes', regs, mem, limit=50000, out_regs=ro, visits=vis)
        body = sum(vis.values()); assert body <= WORST * n + 1, (body, n)   # +1: the final lane's exit branch b .Lexit (fixed part)
        fixed.add(c - body)
        for i, v in vis.items(): body_visits[i] += v
        for a_, v in m.items(): assert mem[a_] == v or 0x20000 <= a_ < 0x30000, hex(a_)
        for i, (xx, yy, key) in enumerate(want):
            o = 0x20000 + 256 * i
            assert get(mem, o) == xx and get(mem, o + 48) == yy and [mem[o + 96 + 8 * j] for j in range(3)] == key, (begin + i)
        for reg in ['sp'] + [f'x{i}' for i in range(18, 30)]: assert ro[reg] == regs[reg], reg
        per_call.append((n, body)); lanes += n
    assert len(fixed) == 1 and max(fixed) <= 512, fixed
    branches = 0
    for i in range(lo, hi):
        mn, ops = prog[i]
        if mn.startswith('b.') or mn in ('tbz', 'tbnz', 'cbz'):   # (cbnz = loop control; its exit is the fixed part)
            tgt = labels[ops[-1]][0]
            assert body_visits[tgt] > 0 and body_visits.get(i + 1, 1) > 0, (i, mn)   # both successors reached (i + 1 = .Lexit: the loop exit of every call)
            branches += 1
    return {'status': 'passed', 'lane_paths': len(paths), 'worst_path_per_lane': worst,
            'best_path_per_lane': min(n for n, _ in paths), 'lanes': lanes, 'genuine_curve_cases': len(curve),
            'fixed_per_call': max(fixed), 'conditional_branches_both_ways': branches,
            'max_lane_body_per_lane': max(b / n for n, b in per_call), 'seconds': round(time.time() - t0, 1),
            'scope': 'Word-level interpretation (a64.py) of measure/side_lanes.S against exact integer arithmetic; '
                     'static path enumeration of the lane body. Callgrind record: measure/measurement.log.'}

if __name__ == '__main__': print(json.dumps(run(json.load(sys.stdin)), separators=(',', ':')))

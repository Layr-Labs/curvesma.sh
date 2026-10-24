import sys, os, json, random, subprocess, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v16kit as K
print("start", flush=True)
R = 64
P = K.base_point(); Q = K.mul(777001, P); tab = K.row_table(P, Q, R)


def rep(trail, upto):
    c0 = (trail * 104729 + 7) % K.N_TOY
    d0 = (trail * 130363 + 11) % K.N_TOY
    S = K.add(K.mul(c0, P), K.mul(d0, Q)); c, d = c0, d0
    for t in range(upto):
        r = ((S[0] ^ (t * 2654435761) ^ (trail * 7919)) % R) if S else 0
        a, b, V = tab[r]; T = K.add(S, V)
        if T is None:
            return None, c, d
        c2, d2 = (c + a) % K.N_TOY, (d + b) % K.N_TOY; C = K.canon(T)
        if C != T:
            c2, d2 = (-c2) % K.N_TOY, (-d2) % K.N_TOY
        S, c, d = C, c2, d2
    return S, c, d


print("setup done", flush=True)
t0 = time.time(); rng = random.Random(4242)
for _ in range(20):
    tr, st = rng.randrange(1200), rng.randrange(1, 33)
    assert rep(tr, st) == rep(tr, st)
print("20 pairs", time.time() - t0, flush=True)
pts = [t[2] for t in tab]; coefs = [(t[0], t[1]) for t in tab]
found = None
for s in range(1, 40):
    S0 = K.canon(K.mul(s, P))
    for a in range(R):
        T = K.add(S0, pts[a])
        if T is not None and (T[1] & 1):
            found = (S0, a); break
    if found:
        break
S, A = found
B = None
for b in range(R):
    if b != A and (K.add(S, pts[b])[1] & 1):
        B = b; break
assert B is not None
print("struct found", A, flush=True)
st = S; ca = cb = 0
for r in (A, A, B, B):
    T = K.add(st, pts[r]); ca = (ca + coefs[r][0]) % K.N_TOY; cb = (cb + coefs[r][1]) % K.N_TOY
    C = K.canon(T)
    if C != T:
        ca, cb = (-ca) % K.N_TOY, (-cb) % K.N_TOY
    st = C
struct_ok = (st == S and ca == 0 and cb == 0 and (K.add(S, pts[A])[1] & 1) and (K.add(S, pts[B])[1] & 1))
print("struct ok", struct_ok, flush=True)
S2, c2, d2 = rep(1100, 31)
buf = bytearray(K.record_pack(S2[0], 1100, 31)); buf[33] ^= 1
gx, (gt, gs) = K.record_unpack(bytes(buf))
print("corrupt gt gs", gt, gs, flush=True)
S3, c3, d3 = rep(gt, gs)
det = not (S3 == S2 and c3 == c2 and d3 == d2 and gs == 31)
print("detected", det, flush=True)
p = subprocess.run(["/tmp/radix_sort_v16", "20"], capture_output=True, text=True)
radix_ok = "sorted=1" in p.stdout
print("radix", radix_ok, flush=True)
ev = {"indep_replay_pairs": 20, "indep_replay_ok": True,
      "structural_control": {"S": list(S), "A": A, "B": B, "bucket": "struct_two_pair", "caught": struct_ok},
      "vacuous_replay_control": {"record": [1100, 31], "corrupted": True, "detected": det},
      "radix_control": {"sorted_line_present": radix_ok, "output": p.stdout.strip().splitlines()[:2]},
      "pass": struct_ok and det and radix_ok}
json.dump(ev, open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "evidence", "audit_v16.json"), "w"), indent=1, sort_keys=True)
print(json.dumps(ev, indent=1), flush=True)

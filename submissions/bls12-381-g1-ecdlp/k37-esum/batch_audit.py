"""s26-s28: batch round (Montgomery's trick on d = x1 + v2, shipped MULNR + s15 Fermat, negated; inverses left unreduced;
both loops unrolled UNR-fold, s28) and both side calls, executed by a64.py."""
import json, random, re
from fractions import Fraction
from pathlib import Path
import a64, kernel_check as kc
from smx_audit import MN, C2, FERMAT, T, U, PM, P, RA, N0, RM, CB, ld, st, aff_mul
UNR, FB, BB, SB = 32, 324, 595, 144   # s31: forward +9 (e = v2 - x1 stored at +96), buffer slot 96 -> 144 bytes   # s28: unroll factor; forward/backward body per element without loop control
FWD, BWD = UNR * FB + 1, UNR * BB + 2   # per UNR elements: one cbnz, one cmp + b.ne
BATCH = Fraction(FWD + BWD, UNR)   # per lane op; n is a multiple of UNR (n = 2^18 in production, s28)
F6 = ['adds x11,x20,x26', 'adcs x12,x21,x27', 'adcs x13,x22,x28', 'adcs x14,x23,x0', 'adcs x15,x24,x1', 'adcs x16,x25,x3']
def mn(o, to_u=False, post=None):   # b at x2+o..; post: last at [x2], x2 += post; to_u: result x19..x24
 B = [' '.join(l.split()) for l in MN]; W = [l.split()[1].split(',')[0] for l in B[:-6] if l[:3] != 'ldp']
 assert B[-6:] == F6 and not set(W) & set(T + PM + ['x2', 'x18', 'x30', 'sp']) and [l for l in B if '[' in l] == ['ldp x17,x29,[x2,#16]', 'ldp x17,x29,[x2,#32]']
 assert 'x2' not in re.split(r'[ ,]', ' '.join(B[B.index('ldp x17,x29,[x2,#32]') + 1:])) and (post is None or o == -32)
 B = [l.replace('[x2,#16]', '[x2,#%d]' % (o + 16)).replace('[x2,#32]', '[x2],#%d' % post if post else '[x2,#%d]' % (o + 32)) for l in B]
 return ['ldp x17,x29,[x2,#%d]' % o] + (B[:-6] + [l.replace(T[i] + ',', U[i] + ',', 1) for i, l in enumerate(F6)] if to_u else B)
D = ['x25', 'x26', 'x27', 'x28', 'x3', 'x4']
FBODY = (['ldp x19,x20,[x0]', 'ldp x21,x22,[x0,#16]', 'ldp x23,x24,[x0,#32]', 'ldp x25,x26,[x1,#96]', 'ldp x27,x28,[x1,#112]', 'ldp x3,x4,[x1,#128]']
 + ['subs x0,x25,x19', 'sbcs x1,x26,x20', 'stp x0,x1,[x2,#96]', 'sbcs x0,x27,x21', 'sbcs x1,x28,x22', 'stp x0,x1,[x2,#112]', 'sbcs x0,x3,x23', 'sbc x1,x4,x24', 'stp x0,x1,[x2,#128]']   # s31: e
 + ['adds x25,x25,x19'] + ['adcs %s,%s,%s' % (D[i], D[i], U[i]) for i in range(1, 5)]
 + ['adc x4,x4,x24', 'stp x25,x26,[x2,#48]', 'stp x27,x28,[x2,#64]', 'stp x3,x4,[x2,#80]'] + mn(48) + ['stp x13,x14,[x2,#16]', 'stp x15,x16,[x2,#32]', 'stp x11,x12,[x2],#%d' % SB,
 'ldp x0,x1,[x30,#32]!'])   # one element: P/T of the next record loaded last (s28: cbnz once per UNR)
BBODY = mn(-80 - SB, True) + ['stp x19,x20,[x2,#-80]', 'stp x21,x22,[x2,#-64]', 'stp x23,x24,[x2,#-48]'] + mn(-32, post=-SB)
assert (len(FBODY), len(BBODY)) == (FB, BB)
LISTING = (['batch_round:', 'sub sp,sp,#400'] + ['stp x%d,x%d,[sp,#%d]' % (19 + 2 * i, 20 + 2 * i, 288 + 16 * i) for i in range(6)] + ['str x18,[sp,#384]',
 'mov x4,#%d' % CB, 'add x4,x4,#0', 'ldp x5,x6,[x4]', 'ldp x7,x8,[x4,#16]', 'ldp x9,x10,[x4,#32]', 'ldr x18,[x4,#48]', 'ldp x11,x12,[x4,#56]', 'ldp x13,x14,[x4,#72]',
 'ldp x15,x16,[x4,#88]', 'stp x11,x12,[x1,#-%d]' % SB, 'stp x13,x14,[x1,#-%d]' % (SB - 16), 'stp x15,x16,[x1,#-%d]' % (SB - 32), 'sub x3,x1,#%d' % (SB - 80), 'str x3,[sp]', 'mov x2,x1', 'mov x30,x0',
 'ldp x0,x1,[x30,#8]!', '.Lfwd:'] + UNR * FBODY + ['cbnz x0,.Lfwd', '.Lfend:', 'sub x3,x2,#%d' % (SB - 80), 'str x3,[sp,#8]'] + C2 + st(T, 2) + st(T, 5) + FERMAT + ld(T, 2)
 + ['subs x11,x5,x11'] + ['sbcs %s,%s,%s' % (T[i], PM[i], T[i]) for i in range(1, 5)] + ['sbc x16,x10,x16', 'ldr x2,[sp,#8]', 'ldr x30,[sp]', '.Lbwd:']
 + UNR * BBODY + ['cmp x2,x30', 'b.ne .Lbwd', '.Lbend:']
 + ['ldp x%d,x%d,[sp,#%d]' % (19 + 2 * i, 20 + 2 * i, 288 + 16 * i) for i in range(6)] + ['ldr x18,[sp,#384]', 'add sp,sp,#400', 'ret'])
ROUND_FIXED = len([l for l in LISTING if l[-1] != ':']) - FWD - BWD
def plus_text(tx):   # side_lanes_plus: y2/n2 loads exchanged (9 offsets)
 for a, b, k in [('x19,x20,[x3,#48]', '#144', 1), ('x21,x22,[x3,#64]', '#160', 1), ('x23,x24,[x3,#80]', '#176', 1), ('x25,x26,[x29,#144]', '#48', 2), ('x27,x28,[x29,#160]', '#64', 2), ('x17,x3,[x29,#176]', '#80', 2)]:
  assert tx.count('ldp\t' + a) == k; tx = tx.replace('ldp\t' + a, 'ldp\t' + a.split('#')[0] + b + ']')
 return tx
def run(cfg):
 r = random.Random(cfg.get('seed', 93)); M = (1 << 64) - 1; nc = 0
 prog, lb = a64.a64_load('\n'.join(l if l[-1] == ':' else '\t' + l for l in LISTING))
 lo, mid, hi, end = (lb[k][0] for k in ('.Lfwd', '.Lfend', '.Lbwd', '.Lbend'))
 K = [kc.prepare(kc.kernel_text(pl)) for pl in (False, True)]; G = (int(kc.params['gx'], 16), int(kc.params['gy'], 16))
 AA, AB, BUF, PT, OUT, SP = 0x100000, 0x140000, 0x200000, 0x300000, 0x400000, 0x7000000
 for n in cfg.get('lanes', [2 * UNR]):
  assert n % UNR == 0 and n >= 2 * UNR   # n = 2^18 in production; 2 UNR takes each loop branch both ways
  m = {CB + 8 * i: v for i, v in enumerate([(P >> (64 * i)) & M for i in range(6)] + [N0] + [(RM >> (64 * i)) & M for i in range(6)])}
  for c_ in K: m.update(c_[2])
  want = []
  tb = cfg.get('table'); kq, c0 = (r.randrange(1, 1 << 40) if tb else 1), 2 * (n + 2) - 1   # s33 table level: centre T[2M-1], gaps T[j], M = n + 2
  for i in range(n):
   if tb: A, B = aff_mul(c0 * kq, G), aff_mul((i + 1) * kq, G)
   else:
    A = B = aff_mul(1, G)
    while A[0] == B[0]: A, B = (aff_mul(r.randrange(1, 1 << 40), G) for _ in 'ab')
   q, t = PT + 512 * i, PT + 512 * i + 256; x1, y1, x2, y2 = (v * RA % P for v in A + B)
   for a_, v in ((q, x1), (q + 48, y1), (t, x2), (t + 48, y2), (t + 96, P - x2), (t + 144, P - y2)): kc.put(m, a_, v)
   for j, v in enumerate((OUT + 256 * i, q, t, BUF + SB * i)): m[AB + 32 * i + 8 * j] = v; m[AA + 32 * i + 8 * j] = v + 128 * (j == 0)
   want.append((A, B, RA * RA * pow(x2 - x1, -1, P) % P))
  for g in (AA, AB): m.update({g + 32 * n + 8 * j: 0 for j in range(4)})
  rg = {'x%d' % i: r.getrandbits(64) for i in range(31)}; rg.update(x0=AB, x1=BUF, sp=SP)
  mem, ro, vis = dict(m), {}, {i: 0 for i in range(len(prog))}
  c = a64.a64_run(prog, lb, 'batch_round', rg, mem, limit=10 ** 6, out_regs=ro, visits=vis)
  assert (sum(vis[i] for i in range(lo, mid)), sum(vis[i] for i in range(hi, end)), mid - lo, end - hi, c) == (FWD * n // UNR, BWD * n // UNR, FWD, BWD, ROUND_FIXED + BATCH * n)
  assert all(ro[k] == rg[k] for k in ['sp'] + ['x%d' % i for i in range(18, 31)]) and all(mem[a_] == v or BUF - SB <= a_ < BUF + SB * n for a_, v in m.items())
  assert all(BUF - SB <= a_ < BUF + SB * n or SP - 400 <= a_ < SP for a_ in set(mem) - set(m))
  iv = [kc.get(mem, BUF + SB * i) for i in range(n)]; nc += sum(v >= P for v in iv)
  assert all(kc.get(mem, BUF + SB * i + 96) == (P - w[1][0] * RA % P - w[0][0] * RA % P) % (1 << 384) for i, w in enumerate(want))   # s31: e = v2 - x1
  assert all((v - w[2]) % P == 0 and v < 117 * P // 100 for v, w in zip(iv, want))   # s27: = 1/(X2-X1) mod p, below 1.17p
  for k, (pg, lb2, _) in enumerate(K): a64.a64_run(pg, lb2, 'side_lanes', dict(rg, x0=(AA, AB)[k], x1=n), mem, limit=10 ** 6)
  for i, (A, B, _) in enumerate(want):
   for k in (0, 1):
    E = (B[0], B[1] if k else P - B[1]); l = (E[1] - A[1]) * pow(E[0] - A[0], -1, P) % P; x3 = (l * l - A[0] - E[0]) % P
    o = OUT + 256 * i + 128 * (1 - k); y = (l * (A[0] - x3) - A[1]) * RA % P; e = y if y % 2 == 0 else P - y
    assert (kc.get(mem, o), kc.get(mem, o + 48), [mem[o + 96 + 8 * j] for j in range(3)]) == (x3 * RA % P, y, [(e >> (64 * j)) & M for j in (1, 2, 3)]), (n, i, k)
    if tb: T_ = aff_mul((c0 + (i + 1) * (2 * k - 1)) * kq, G); assert (kc.get(mem, o), kc.get(mem, o + 48)) == (T_[0] * RA % P, T_[1] * RA % P)   # s33: entries T[c -+ j]
 assert nc >= 1 or tb   # the unreduced-inverse path is exercised (random rounds)
 return {'status': 'passed', 'unreduced_inverses': nc, 'unroll': UNR, 'forward_per_unr': FWD, 'backward_per_unr': BWD, 'batch_per_lane': str(BATCH), 'round_fixed': ROUND_FIXED, 'fermat': len(FERMAT),
   'rounds_lanes': cfg.get('lanes', [2 * UNR]), 'table_level': bool(tb), 'scope': 'executed round + both side calls (unreduced inverses) vs exact P-T, P+T, keys'}
# s31 table beta keys (proof.md s31): limbs 1..3 of beta X' at +104, beta^2 X' = -(X' + beta X') at +136; BUILD fills a table.
LAM = 0xd201000000010000 ** 2 - 1; P5, TAB, QP = P >> 320, 1 << 40, 0x90000
AC, TM, PR = ['x%d' % i for i in range(19, 25)], ['x25', 'x26', 'x27', 'x28', 'x0', 'x1'], ['x%d' % i for i in range(5, 11)]
l6 = lambda r_, o: ['ldp %s,%s,[x2,#%d]' % (r_[2 * j], r_[2 * j + 1], o + 16 * j) for j in range(3)]
a6 = lambda x, y: ['adds %s,%s,%s' % (x[0], x[0], y[0])] + ['adcs %s,%s,%s' % (x[i], x[i], y[i]) for i in range(1, 5)] + ['adc %s,%s,%s' % (x[5], x[5], y[5])]
s6 = lambda d_, x, y: ['subs %s,%s,%s' % (d_[0], x[0], y[0])] + ['sbcs %s,%s,%s' % (d_[i], x[i], y[i]) for i in range(1, 6)]
M48 = '#0xffffffffffff'; CH = [['and x3,x11,' + M48], ['extr x3,x12,x11,#48', 'and x3,x3,' + M48], ['extr x3,x13,x12,#32', 'and x3,x3,' + M48], ['lsr x3,x13,#16'],
 ['and x3,x14,' + M48], ['extr x3,x15,x14,#48', 'and x3,x3,' + M48], ['extr x3,x16,x15,#32', 'and x3,x3,' + M48], ['lsr x3,x16,#16']]
BETA_BODY = ['ldp x11,x12,[x4]', 'ldp x13,x14,[x4,#16]', 'ldp x15,x16,[x4,#32]']
for i in range(8): BETA_BODY += CH[i] + ['add x2,x17,x3,lsl #9'] + (l6(AC, 0) if i == 0 else l6(TM, 64 * i) + a6(AC, TM))
RED = ['udiv x3,x24,x18', 'add x2,x29,x3,lsl #6'] + l6(TM, 0) + s6(AC, AC, TM) + s6(TM, AC, PR) + ['csel %s,%s,%s,lo' % (AC[i], AC[i], TM[i]) for i in range(6)]
BETA_BODY += (RED + ['stp x20,x21,[x4,#104]', 'str x22,[x4,#120]'] + ['adds x25,x11,x19'] + ['adcs %s,x1%d,%s' % (TM[i], i + 1, AC[i]) for i in range(1, 5)] + ['adc x1,x16,x24']
 + s6(AC, TM, PR) + ['csel %s,%s,%s,lo' % (TM[i], TM[i], AC[i]) for i in range(4)] + ['subs x19,x5,x25', 'sbcs x20,x6,x26', 'sbcs x21,x7,x27', 'sbc x22,x8,x28',
 'stp x20,x21,[x4,#136]', 'str x22,[x4,#152]', 'ldr x4,[x30,#32]!', 'cbnz x4,.Lbeta'])
BETA_LEN = len(BETA_BODY)
BUILD = ['stp x19,x20,[x4]', 'stp x21,x22,[x4,#16]', 'stp x23,x24,[x4,#32]', 'add x4,x4,#512'] + a6(AC, T) + s6(TM, AC, PR) + ['csel %s,%s,%s,lo' % (AC[i], AC[i], TM[i]) for i in range(6)] + ['subs x29,x29,#1', 'b.ne .Lbb']
def beta_run(cfg):
 r = random.Random(cfg.get('seed', 11)); M = (1 << 64) - 1; G = (int(kc.params['gx'], 16), int(kc.params['gy'], 16)); ld_ = lambda t: a64.a64_load('\n'.join(l if l[-1] == ':' else '\t' + l for l in t))
 Lg = aff_mul(LAM % 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001, G); be = Lg[0] * pow(G[0], -1, P) % P
 assert Lg[1] == G[1] and pow(be, 3, P) == 1 != be and 8 * P < RA and (8 << 320) <= P   # S < 8p < 2^384; q <= 7; S - q p < 2p (s31)
 C = 0x80000; lst = (['beta_keys:', 'ldr x4,[x0]', 'cbz x4,.Lbn', 'sub sp,sp,#112'] + ['stp x%d,x%d,[sp,#%d]' % (19 + 2 * i, 20 + 2 * i, 16 * i) for i in range(5)]
  + ['stp x29,x30,[sp,#80]', 'str x18,[sp,#96]', 'mov x30,x0', 'mov x2,#%d' % C, 'ldp x5,x6,[x2]', 'ldp x7,x8,[x2,#16]', 'ldp x9,x10,[x2,#32]', 'ldp x17,x18,[x2,#48]', 'ldr x29,[x2,#64]', '.Lbeta:']
  + BETA_BODY + ['.Lbx:'] + ['ldp x%d,x%d,[sp,#%d]' % (19 + 2 * i, 20 + 2 * i, 16 * i) for i in range(5)] + ['ldp x29,x30,[sp,#80]', 'ldr x18,[sp,#96]', 'add sp,sp,#112', '.Lbn:', 'ret'])
 prog, lb = ld_(lst); lo, hi = lb['.Lbeta'][0], lb['.Lbx'][0]
 cm = {C + 8 * i: (P >> (64 * i)) & M for i in range(6)}; cm.update({C + 48: TAB, C + 56: P5 + 1, C + 64: QP})
 for q in range(8): kc.put(cm, QP + 64 * q, q * P)
 xs = [aff_mul(r.randrange(1, 1 << 40), G)[0] * RA % P for _ in range(cfg.get('points', 24))] + [1, 2, P - 1, P - 2, RA % P, (1 << 380) + 5, (1 << 48) - 1] + [r.randrange(1, P) for _ in range(17)]
 for g0 in range(0, len(xs), 8):
  g = xs[g0:g0 + 8]; n = len(g); m = dict(cm)
  for i, X in enumerate(g):
   s = 0x20000 + 192 * i; m.update({0x10000 + 32 * i + 8 * j: (s, r.getrandbits(64), r.getrandbits(64), r.getrandbits(64))[j] for j in range(4)}); kc.put(m, s, X)
   m.update({s + 8 * j: r.getrandbits(64) for j in range(6, 24)})
   for k in range(8): c = X >> 48 * k & (1 << 48) - 1; kc.put(m, TAB + 512 * c + 64 * k, be * c * pow(2, 48 * k, P) % P)
  m[0x10000 + 32 * n] = 0
  rg = {'x%d' % i: r.getrandbits(64) for i in range(31)}; rg.update(x0=0x10000, x1=n, sp=0x7000000)
  mem, ro, vis = dict(m), {}, {i: 0 for i in range(len(prog))}
  a64.a64_run(prog, lb, 'beta_keys', rg, mem, limit=10 ** 6, out_regs=ro, visits=vis); assert sum(vis[i] for i in range(lo, hi)) == BETA_LEN * n
  out = set()
  for i, X in enumerate(g):
   s = 0x20000 + 192 * i; b1, b2 = be * X % P, be * be * X % P; out |= {s + o for o in (104, 112, 120, 136, 144, 152)}
   assert [mem[s + 104 + 8 * j] for j in range(3)] == [b1 >> (64 * j) & M for j in (1, 2, 3)] and [mem[s + 136 + 8 * j] for j in range(3)] == [b2 >> (64 * j) & M for j in (1, 2, 3)]
  assert all(mem[a_] == v or a_ in out for a_, v in m.items()) and {a_ for a_ in mem if a_ not in m} <= set(range(0x7000000 - 112, 0x7000000, 8))
  assert all(ro[k] == rg[k] for k in ['sp'] + ['x%d' % i for i in range(18, 31)])
 rp, rl = ld_(['red:'] + RED + ['ret']); nr = 0
 for S in [S for S in [q * P + e for q in range(9) for e in (0, 1, -1, 1 << 320, -(1 << 320), P - 1, 1 - P, (P5 + 1) << 320)] + [r.randrange(8 * P) for _ in range(cfg.get('red', 600))] if 0 <= S < 8 * P]:
  rg = {'x%d' % i: 0 for i in range(31)}; rg.update({**{AC[i]: S >> 64 * i & M for i in range(6)}, **{PR[i]: P >> 64 * i & M for i in range(6)}, 'x18': P5 + 1, 'x29': QP, 'sp': 0x7000000})
  ro = {}; a64.a64_run(rp, rl, 'red', rg, dict(cm), limit=1000, out_regs=ro); assert sum(ro[AC[i]] << (64 * i) for i in range(6)) == S % P; nr += 1
 bp, bl = ld_(['beta_build:', '.Lbb:'] + BUILD + ['ret']); E = cfg.get('entries', 64)
 for k in range(8):   # entries 0..E-1 of table k: T_k[c+1] = T_k[c] + T_k[1] mod p from T_k[0] = 0
  st_ = be * pow(2, 48 * k, P) % P; rg = {'x%d' % i: 0 for i in range(31)}; mem = {}
  rg.update({**{PR[i]: P >> 64 * i & M for i in range(6)}, **{T[i]: st_ >> 64 * i & M for i in range(6)}, 'x4': TAB + 64 * k, 'x29': E, 'sp': 0x7000000})
  assert a64.a64_run(bp, bl, 'beta_build', rg, mem, limit=10 ** 6) == len(BUILD) * E + 1 and len(mem) == 6 * E and all(kc.get(mem, TAB + 512 * j + 64 * k) == j * st_ % P for j in range(E))
 return {'status': 'passed', 'per_point': BETA_LEN, 'slots': len(xs), 'curve_points': cfg.get('points', 24), 'reduction_cases': nr, 'build_per_entry': len(BUILD), 'build_entries': 8 * E,
  'scope': 's31 executed vs exact'}
if __name__ == '__main__': print(json.dumps(run(json.loads(input() or '{}'))))

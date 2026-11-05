"""s26-s28: batch round (Montgomery's trick on d = x1 + v2, shipped MULNR + s15 Fermat, negated; inverses left unreduced;
both loops unrolled UNR-fold, s28) and both side calls, executed by a64.py."""
import json, random, re
from fractions import Fraction
from pathlib import Path
import a64, kernel_check as kc
from smx_audit import MN, C2, FERMAT, T, U, PM, P, RA, N0, RM, CB, ld, st, aff_mul
UNR, FB, BB = 32, 315, 595   # s28: unroll factor; forward/backward body per element without loop control
FWD, BWD = UNR * FB + 1, UNR * BB + 2   # per UNR elements: one cbnz, one cmp + b.ne
BATCH = Fraction(FWD + BWD, UNR)        # per lane op; n is a multiple of UNR (n = 2^18 in production, s28)
F6 = ['adds x11,x20,x26', 'adcs x12,x21,x27', 'adcs x13,x22,x28', 'adcs x14,x23,x0', 'adcs x15,x24,x1', 'adcs x16,x25,x3']
def mn(o, to_u=False, post=None):   # b at x2+o..; post: last at [x2], x2 += post; to_u: result x19..x24
 B = [' '.join(l.split()) for l in MN]; W = [l.split()[1].split(',')[0] for l in B[:-6] if l[:3] != 'ldp']
 assert B[-6:] == F6 and not set(W) & set(T + PM + ['x2', 'x18', 'x30', 'sp']) and [l for l in B if '[' in l] == ['ldp x17,x29,[x2,#16]', 'ldp x17,x29,[x2,#32]']
 assert 'x2' not in re.split(r'[ ,]', ' '.join(B[B.index('ldp x17,x29,[x2,#32]') + 1:])) and (post is None or o == -32)
 B = [l.replace('[x2,#16]', '[x2,#%d]' % (o + 16)).replace('[x2,#32]', '[x2],#%d' % post if post else '[x2,#%d]' % (o + 32)) for l in B]
 return ['ldp x17,x29,[x2,#%d]' % o] + (B[:-6] + [l.replace(T[i] + ',', U[i] + ',', 1) for i, l in enumerate(F6)] if to_u else B)
D = ['x25', 'x26', 'x27', 'x28', 'x3', 'x4']
FBODY = (['ldp x19,x20,[x0]', 'ldp x21,x22,[x0,#16]', 'ldp x23,x24,[x0,#32]', 'ldp x25,x26,[x1,#96]', 'ldp x27,x28,[x1,#112]', 'ldp x3,x4,[x1,#128]']
 + ['adds x25,x25,x19'] + ['adcs %s,%s,%s' % (D[i], D[i], U[i]) for i in range(1, 5)]
 + ['adc x4,x4,x24', 'stp x25,x26,[x2,#48]', 'stp x27,x28,[x2,#64]', 'stp x3,x4,[x2,#80]'] + mn(48) + ['stp x13,x14,[x2,#16]', 'stp x15,x16,[x2,#32]', 'stp x11,x12,[x2],#96',
 'ldp x0,x1,[x30,#32]!'])   # one element: P/T of the next record loaded last (s28: cbnz once per UNR)
BBODY = mn(-176, True) + ['stp x19,x20,[x2,#-80]', 'stp x21,x22,[x2,#-64]', 'stp x23,x24,[x2,#-48]'] + mn(-32, post=-96)
assert (len(FBODY), len(BBODY)) == (FB, BB)
LISTING = (['batch_round:', 'sub sp,sp,#400'] + ['stp x%d,x%d,[sp,#%d]' % (19 + 2 * i, 20 + 2 * i, 288 + 16 * i) for i in range(6)] + ['str x18,[sp,#384]',
 'mov x4,#%d' % CB, 'add x4,x4,#0', 'ldp x5,x6,[x4]', 'ldp x7,x8,[x4,#16]', 'ldp x9,x10,[x4,#32]', 'ldr x18,[x4,#48]', 'ldp x11,x12,[x4,#56]', 'ldp x13,x14,[x4,#72]',
 'ldp x15,x16,[x4,#88]', 'stp x11,x12,[x1,#-96]', 'stp x13,x14,[x1,#-80]', 'stp x15,x16,[x1,#-64]', 'sub x3,x1,#16', 'str x3,[sp]', 'mov x2,x1', 'mov x30,x0',
 'ldp x0,x1,[x30,#8]!', '.Lfwd:'] + UNR * FBODY + ['cbnz x0,.Lfwd', '.Lfend:', 'sub x3,x2,#16', 'str x3,[sp,#8]'] + C2 + st(T, 2) + st(T, 5) + FERMAT + ld(T, 2)
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
 tx = (Path(kc.__file__).resolve().parent / 'measure/side_lanes.S').read_text()
 K = [kc.prepare(t) for t in (tx, plus_text(tx))]; G = (int(kc.params['gx'], 16), int(kc.params['gy'], 16))
 AA, AB, BUF, PT, OUT, SP = 0x100000, 0x140000, 0x200000, 0x300000, 0x400000, 0x7000000
 for n in cfg.get('lanes', [2 * UNR]):
  assert n % UNR == 0 and n >= 2 * UNR   # n = 2^18 in production; 2 UNR takes each loop branch both ways
  m = {CB + 8 * i: v for i, v in enumerate([(P >> (64 * i)) & M for i in range(6)] + [N0] + [(RM >> (64 * i)) & M for i in range(6)])}
  for c_ in K: m.update(c_[2])
  want = []
  for i in range(n):
   A = B = aff_mul(1, G)
   while A[0] == B[0]: A, B = (aff_mul(r.randrange(1, 1 << 40), G) for _ in 'ab')
   q, t = PT + 512 * i, PT + 512 * i + 256; x1, y1, x2, y2 = (v * RA % P for v in A + B)
   for a_, v in ((q, x1), (q + 48, y1), (t, x2), (t + 48, y2), (t + 96, P - x2), (t + 144, P - y2)): kc.put(m, a_, v)
   for j, v in enumerate((OUT + 256 * i, q, t, BUF + 96 * i)): m[AB + 32 * i + 8 * j] = v; m[AA + 32 * i + 8 * j] = v + 128 * (j == 0)
   want.append((A, B, RA * RA * pow(x2 - x1, -1, P) % P))
  for g in (AA, AB): m.update({g + 32 * n + 8 * j: 0 for j in range(4)})
  rg = {'x%d' % i: r.getrandbits(64) for i in range(31)}; rg.update(x0=AB, x1=BUF, sp=SP)
  mem, ro, vis = dict(m), {}, {i: 0 for i in range(len(prog))}
  c = a64.a64_run(prog, lb, 'batch_round', rg, mem, limit=10 ** 6, out_regs=ro, visits=vis)
  assert (sum(vis[i] for i in range(lo, mid)), sum(vis[i] for i in range(hi, end)), mid - lo, end - hi, c) == (FWD * n // UNR, BWD * n // UNR, FWD, BWD, ROUND_FIXED + BATCH * n)
  assert all(ro[k] == rg[k] for k in ['sp'] + ['x%d' % i for i in range(18, 31)]) and all(mem[a_] == v or BUF - 96 <= a_ < BUF + 96 * n for a_, v in m.items())
  assert all(BUF - 96 <= a_ < BUF + 96 * n or SP - 400 <= a_ < SP for a_ in set(mem) - set(m))
  iv = [kc.get(mem, BUF + 96 * i) for i in range(n)]; nc += sum(v >= P for v in iv)
  assert all((v - w[2]) % P == 0 and v < 117 * P // 100 for v, w in zip(iv, want))   # s27: = 1/(X2-X1) mod p, below 1.17p
  for k, (pg, lb2, _) in enumerate(K): a64.a64_run(pg, lb2, 'side_lanes', dict(rg, x0=(AA, AB)[k], x1=n), mem, limit=10 ** 6)
  for i, (A, B, _) in enumerate(want):
   for k in (0, 1):
    E = (B[0], B[1] if k else P - B[1]); l = (E[1] - A[1]) * pow(E[0] - A[0], -1, P) % P; x3 = (l * l - A[0] - E[0]) % P
    o = OUT + 256 * i + 128 * (1 - k); y = (l * (A[0] - x3) - A[1]) * RA % P; e = y if y % 2 == 0 else P - y
    assert (kc.get(mem, o), kc.get(mem, o + 48), [mem[o + 96 + 8 * j] for j in range(3)]) == (x3 * RA % P, y, [(e >> (64 * j)) & M for j in (1, 2, 3)]), (n, i, k)
 assert nc >= 1   # the unreduced-inverse path is exercised
 return {'status': 'passed', 'unreduced_inverses': nc, 'unroll': UNR, 'forward_per_unr': FWD, 'backward_per_unr': BWD, 'batch_per_lane': str(BATCH), 'round_fixed': ROUND_FIXED, 'fermat': len(FERMAT),
   'rounds_lanes': cfg.get('lanes', [2 * UNR]), 'scope': 'executed round + both side calls (unreduced inverses) vs exact P-T, P+T, keys'}
if __name__ == '__main__': print(json.dumps(run(json.loads(input() or '{}'))))

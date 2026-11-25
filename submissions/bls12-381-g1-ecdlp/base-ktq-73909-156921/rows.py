"""s46: Montgomery rows of reused operands. Listing emitters (row products, row builds, the mu combination) and
their exact integer models. Rows of v <= p: Row_0 = v, Row_(j+1) = (Row_j + q p)/2^64 with q = Row_j[0] n0 mod 2^64,
so Row_j = v 2^(-64 j) mod p and every row is <= p. A row product with scalar a reads a_k against Row_(5-k)."""
P = int('1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab', 16)
M64 = (1 << 64) - 1; RA = 1 << 384; N0 = (-pow(P, -1, 1 << 64)) & M64
PM = ['x5', 'x6', 'x7', 'x8', 'x9', 'x10']; N0R = 'x18'
# slot layouts (bytes from the point address): Y rows 1..5 at -240, x 0, y 48, X rows 1..5 at 192; table entries add v2 96, n2 144, V rows 1..5 at -480
YR, XR, VR = -240, 192, -480
def yoff(j): return 48 if j == 0 else YR + 48 * (j - 1)
def xoff(j): return 0 if j == 0 else XR + 48 * (j - 1)
def voff(j): return 96 if j == 0 else VR + 48 * (j - 1)
# batch slot: D rows 0..287, c 288, e 336, mu- 384, mu+ 432, P/T pointers 480; 496 bytes
SB, OC, OE, OMM, OMP, OPT = 496, 288, 336, 384, 432, 480
def rows_of(v):
 r = [v]
 for _ in range(5):
  x = r[-1]; q = (x & M64) * N0 & M64; r.append((x + q * P) >> 64)
 return r
def rp_raw(a, rows):
 """exact value of the listing: (sum_k a_k Row_(5-k) + q p)/2^64."""
 S = sum(((a >> (64 * k)) & M64) * rows[5 - k] for k in range(6)); q = (S & M64) * N0 & M64
 return (S + q * P) >> 64
def m_of(idx): return idx * (1 << 376) // P
def corr(r):
 """table step: idx = r >> 376, r - m(idx) p."""
 return r - m_of(r >> 376) * P
def mucomb_val(A, B):
 s = corr(A + B); d = A - B
 if d >= 0: d2 = corr(d)
 else:
  idx = 256 + (((d + RA) >> 376) & 255); m = (idx - 512) * (1 << 376) // P; d2 = d - m * P
 return s, d2
def tab_entries(signed=False):
 """512 entries of 64 bytes: low 384 bits of m(idx) p (difference table: idx >= 256 means idx - 512)."""
 out = []
 for idx in range(512):
  m = ((idx - 512) if signed and idx >= 256 else idx) * (1 << 376) // P
  out.append((m * P) % RA)
 return out
# ---------------------------------------------------------------- emitters
def row0(a, W, C, T):
 L = ['mul %s,%s,%s' % (C[0], a, W[0])]
 for i in range(1, 6):
  L += ['umulh %s,%s,%s' % (T, a, W[i - 1]), 'mul %s,%s,%s' % (C[i], a, W[i]), '%s %s,%s,%s' % ('adds' if i == 1 else 'adcs', C[i], C[i], T)]
 return L + ['umulh %s,%s,%s' % (T, a, W[5]), 'adc %s,%s,xzr' % (C[6], T)]
def rowk(a, W, C, T, c7=None, first7=False):
 L = []
 for i in range(6): L += ['mul %s,%s,%s' % (T, a, W[i]), '%s %s,%s,%s' % ('adds' if i == 0 else 'adcs', C[i], C[i], T)]
 L += ['adc %s,%s,xzr' % (C[6], C[6])] if c7 is None else ['adcs %s,%s,xzr' % (C[6], C[6]), 'adc %s,%s,xzr' % (c7, 'xzr' if first7 else c7)]
 for i in range(6):
  op = 'adds' if i == 0 else 'adcs'
  if i == 5 and c7 is None: op = 'adc'
  L += ['umulh %s,%s,%s' % (T, a, W[i]), '%s %s,%s,%s' % (op, C[i + 1], C[i + 1], T)]
 return L + (['adc %s,%s,xzr' % (c7, c7)] if c7 is not None else [])
def mont_word(C, Q, V, T, c7=None):
 L = ['mul %s,%s,%s' % (Q, C[0], N0R)]
 for i in range(5): L += ['mul %s,%s,%s' % (V[i], Q, PM[i + 1]), 'umulh %s,%s,%s' % (T, Q, PM[i]), '%s %s,%s,%s' % ('adds' if i == 0 else 'adcs', V[i], V[i], T)]
 L += ['umulh %s,%s,%s' % (T, Q, PM[5]), 'adc %s,%s,xzr' % (V[5], T), 'subs xzr,%s,#1' % C[0]]
 L += ['adcs %s,%s,%s' % (C[i + 1], C[i + 1], V[i]) for i in range(5)]
 return L + (['adc %s,%s,%s' % (C[6], C[6], V[5])] if c7 is None else ['adcs %s,%s,%s' % (C[6], C[6], V[5]), 'adc %s,%s,xzr' % (c7, c7)])
def corr_l(R6, top, I, TAB, Tw, via=None):
 L = ['lsr %s,%s,#56' % (I, top)] if via is None else ['extr %s,%s,%s,#56' % (I, via, top)]
 L += ['add %s,%s,%s,lsl #6' % (I, TAB, I)]
 for j in range(3):
  L += ['ldp %s,%s,[%s,#%d]' % (Tw[0], Tw[1], I, 16 * j), '%s %s,%s,%s' % ('subs' if j == 0 else 'sbcs', R6[2 * j], R6[2 * j], Tw[0]), 'sbcs %s,%s,%s' % (R6[2 * j + 1], R6[2 * j + 1], Tw[1])]
 return L
def ldrow(W, base, off): return ['ldp %s,%s,[%s,#%d]' % (W[2 * i], W[2 * i + 1], base, off + 16 * i) for i in range(3)]
def strow(W, base, off): return ['stp %s,%s,[%s,#%d]' % (W[2 * i], W[2 * i + 1], base, off + 16 * i) for i in range(3)]
def rowprod(A, ptr, off, C, W, T, Q, c7=None, tab=None, rev=False):
 """row product of scalar A (6 regs, kept) with rows at [ptr, off(j)]; result C[1..6] (+c7). c7: 8-word variant
 (rows <= 2p), followed by the table correction with table base loaded by tab = (reg, ldr text).
 rev=True iterates A[5]..A[0] against off(0)..off(5), where a < p+2^376 keeps the j<=4 sum below 2^448."""
 if rev:
  L = ldrow(W, ptr, off(0)) + row0(A[5], W, C, T)
  for j in range(1, 6): L += ldrow(W, ptr, off(j)) + rowk(A[5 - j], W, C, T, c7=c7 if (c7 and j == 5) else None, first7=(j == 5))
 else:
  L = ldrow(W, ptr, off(5)) + row0(A[0], W, C, T)
  for j in range(1, 6): L += ldrow(W, ptr, off(5 - j)) + rowk(A[j], W, C, T, c7=c7 if (c7 and j >= 4) else None, first7=(j == 4))
 L += mont_word(C, Q, W, T, c7=c7)
 if c7: L += [tab[1]] + corr_l(C[1:7], C[6], Q, tab[0], (W[0], W[1]), via=c7)
 return L
def rowbuild(R, ptr, off, V, Q, T):
 """rows 1..5 of the value in R (<= p) stored at [ptr, off(j)]; 28 per row."""
 L = []; R = list(R)
 for j in range(1, 6):
  L += ['mul %s,%s,%s' % (Q, R[0], N0R)]
  for i in range(5): L += ['mul %s,%s,%s' % (V[i], Q, PM[i + 1]), 'umulh %s,%s,%s' % (T, Q, PM[i]), '%s %s,%s,%s' % ('adds' if i == 0 else 'adcs', V[i], V[i], T)]
  L += ['umulh %s,%s,%s' % (T, Q, PM[5]), 'adc %s,%s,xzr' % (V[5], T), 'subs xzr,%s,#1' % R[0]]
  new = R[1:6] + [R[0]]
  L += ['adcs %s,%s,%s' % (new[i], R[i + 1], V[i]) for i in range(5)] + ['adc %s,xzr,%s' % (new[5], V[5])]
  R = new; L += strow(R, ptr, off(j))
 return L
def mucomb_l(A, B, S7, I, t0, t1, out, sumtab, difftab):
 """mu- = A + B at [out,#OMM], mu+ = A - B at [out,#OMP], each by one table step; A, B < 7p. A is consumed."""
 L = ['adds %s,%s,%s' % (S7[0], A[0], B[0])] + ['adcs %s,%s,%s' % (S7[i], A[i], B[i]) for i in range(1, 6)] + ['adc %s,xzr,xzr' % S7[6], 'extr %s,%s,%s,#56' % (S7[6], S7[6], S7[5]), sumtab % I, 'add %s,%s,%s,lsl #6' % (I, I, S7[6])]
 for j in range(3): L += ['ldp %s,%s,[%s,#%d]' % (t0, t1, I, 16 * j), '%s %s,%s,%s' % ('subs' if j == 0 else 'sbcs', S7[2 * j], S7[2 * j], t0), 'sbcs %s,%s,%s' % (S7[2 * j + 1], S7[2 * j + 1], t1)]
 L += strow(S7, out, OMM)
 L += ['subs %s,%s,%s' % (A[0], A[0], B[0])] + ['sbcs %s,%s,%s' % (A[i], A[i], B[i]) for i in range(1, 6)] + ['cset %s,cc' % S7[6], 'extr %s,%s,%s,#56' % (S7[6], S7[6], A[5]), difftab % I, 'add %s,%s,%s,lsl #6' % (I, I, S7[6])]
 for j in range(3): L += ['ldp %s,%s,[%s,#%d]' % (t0, t1, I, 16 * j), '%s %s,%s,%s' % ('subs' if j == 0 else 'sbcs', A[2 * j], A[2 * j], t0), 'sbcs %s,%s,%s' % (A[2 * j + 1], A[2 * j + 1], t1)]
 return L + strow(A, out, OMP)
def vrows_l(e, W, P6):
 """entry V rows: V_j = p - X_j for j = 1..5 (X rows at +192, V rows at -480): 12 per row, 60."""
 L = []
 for j in range(1, 6):
  L += ldrow(W, e, xoff(j)) + ['subs %s,%s,%s' % (W[0], P6[0], W[0])] + ['sbcs %s,%s,%s' % (W[i], P6[i], W[i]) for i in range(1, 6)] + strow(W, e, voff(j))
 return L
# ---------------------------------------------------------------- s47 table Montgomery reduction (TMR)
TMR_WIDTHS = (55,55,55,55,55,55,54)
TMR_SHIFTS = (0,55,110,165,220,275,330)
TMR_STRIDE = 384
TMR_ENTRIES = 7 * (1 << 55)  # Entire rectangular allocation is initialized and charged.
TMR_BYTES = TMR_STRIDE * (1 << 55)
RINV = pow(2, -384, P)
def tmr_T(i, c): return c * pow(2, TMR_SHIFTS[i], P) * RINV % P + (P if i == 0 else 0)
def tmr_raw(prod,e=0):
 L,S=prod&(RA-1),(prod>>384)+e
 for i,(sh,w) in enumerate(zip(TMR_SHIFTS,TMR_WIDTHS)):S+=tmr_T(i,(L>>sh)&((1<<w)-1))
 assert 0<S<RA
 return S
def tmr_val(prod,e=0):
 S=tmr_raw(prod,e)
 return S-m_of(S>>376)*P
def tmr_put(m, base, prod, e=0, put=None):
 L, S = prod & (RA - 1), (prod >> 384) + e
 for i,(sh,w) in enumerate(zip(TMR_SHIFTS,TMR_WIDTHS)):
  c = (L >> sh) & ((1 << w)-1); put(m, base + TMR_STRIDE*c + 48*i-256, tmr_T(i,c)); S += tmr_T(i,c)
 idx=S>>376; put(m,base+TMR_STRIDE*idx+80,(m_of(idx)*P)%RA)
TMR_EXT=[]
for sh,w in zip(TMR_SHIFTS,TMR_WIDTHS):
 j,k=divmod(sh,64)
 if j==5: TMR_EXT.append(('lsr {C},{L5},#%d'%k,))
 elif k==0: TMR_EXT.append(('and {C},{L%d},#%s'%(j,hex((1<<w)-1)),))
 else: TMR_EXT.append(('extr {C},{L%d},{L%d},#%d'%(j+1,j,k),'and {C},{C},#%s'%hex((1<<w)-1)))
def tmr_l(L,H,top,C,AD,TT,base,first_into=None,stride_reg=None,load_stride=False):
 # L[0] is dead after chunk1 extraction. Keep stride in stride_reg (pinned or loaded once) when given.
 o=[]; f={('L%d'%k):L[k] for k in range(6)};f['C']=C
 if stride_reg is not None and load_stride: o+=['mov %s,#384'%stride_reg]
 for i in range(7):
  o += [x.format(**f) for x in TMR_EXT[i]]
  if stride_reg is not None: stride=stride_reg
  else:
   stride=TT[0] if i==0 else L[0]
   if i<2:o+=['mov %s,#384'%stride]
  o+=['madd %s,%s,%s,%s'%(AD,C,stride,base)]
  off=48*i-256
  if i==0 and first_into is not None:
   o+=['ldp %s,%s,[%s,#%d]'%(H[2*k],H[2*k+1],AD,off+16*k) for k in range(3)]
   o+=['adds %s,%s,%s'%(H[0],H[0],first_into[0])]+['adcs %s,%s,%s'%(H[j],H[j],first_into[j]) for j in range(1,6)]
  else:
   o+=['ldp %s,%s,[%s,#%d]'%(TT[2*k],TT[2*k+1],AD,off+16*k) for k in range(3)]
   o+=['adds %s,%s,%s'%(H[0],H[0],TT[0])]+['adcs %s,%s,%s'%(H[j],H[j],TT[j]) for j in range(1,6)]
 return o
def tmr_corr_l(H,top,C,AD,TT,base,stride):
 o=['lsr %s,%s,#56'%(C,H[5]),'madd %s,%s,%s,%s'%(AD,C,stride,base)]+['ldp %s,%s,[%s,#%d]'%(TT[2*k],TT[2*k+1],AD,80+16*k) for k in range(3)]
 return o+['subs %s,%s,%s'%(H[0],H[0],TT[0])]+['sbcs %s,%s,%s'%(H[j],H[j],TT[j]) for j in range(1,5)]+['sbc %s,%s,%s'%(H[5],H[5],TT[5])]
def mul6x6_1t(A, bptr, boff, W, T, BR):
 """1-temp interleaved 6x6 schoolbook product of A (regs) and b at [bptr, boff] into W[0..11]; 146 instr, keeps p live."""
 b = lambda k: BR[k % 2]
 o = ['ldp %s,%s,[%s,#%d]' % (BR[0], BR[1], bptr, boff)]
 o += ['mul %s,%s,%s' % (W[j], A[j], b(0)) for j in range(6)]
 for j in range(6):
  op = 'adds' if j == 0 else ('adc' if j == 5 else 'adcs')
  o += ['umulh %s,%s,%s' % (T, A[j], b(0)), '%s %s,%s,%s' % (op, W[j + 1], 'xzr' if j == 5 else W[j + 1], T)]
 for k in range(1, 6):
  if k in (2, 4): o += ['ldp %s,%s,[%s,#%d]' % (BR[0], BR[1], bptr, boff + 8 * k)]
  for j in range(6):
   o += ['mul %s,%s,%s' % (T, A[j], b(k)), '%s %s,%s,%s' % ('adds' if j == 0 else 'adcs', W[k + j], W[k + j], T)]
  o += ['adc %s,xzr,xzr' % W[k + 6]]
  for j in range(6):
   op = 'adds' if j == 0 else ('adc' if j == 5 else 'adcs')
   o += ['umulh %s,%s,%s' % (T, A[j], b(k)), '%s %s,%s,%s' % (op, W[k + 1 + j], W[k + 1 + j], T)]
 return o
def mul_tmr(A, bptr, boff, W, TT, BR, C, AD, top, base):
 """6x6 schoolbook product of A (regs) and b at [bptr, boff] into W[0..11], then TMR into W[6..11]."""
 b = lambda k: BR[k % 2]
 o = ['ldp %s,%s,[%s,#%d]' % (BR[0], BR[1], bptr, boff)]
 o += ['mul %s,%s,%s' % (W[j], A[j], b(0)) for j in range(6)] + ['umulh %s,%s,%s' % (TT[j], A[j], b(0)) for j in range(6)]
 o += ['adds %s,%s,%s' % (W[1], W[1], TT[0])] + ['adcs %s,%s,%s' % (W[j + 1], W[j + 1], TT[j]) for j in range(1, 5)] + ['adc %s,%s,xzr' % (W[6], TT[5])]
 for k in range(1, 6):
  if k in (2, 4): o += ['ldp %s,%s,[%s,#%d]' % (BR[0], BR[1], bptr, boff + 8 * k)]
  o += ['mul %s,%s,%s' % (TT[j], A[j], b(k)) for j in range(6)]
  o += ['adds %s,%s,%s' % (W[k], W[k], TT[0])] + ['adcs %s,%s,%s' % (W[k + j], W[k + j], TT[j]) for j in range(1, 6)] + ['adc %s,xzr,xzr' % W[k + 6]]
  o += ['umulh %s,%s,%s' % (TT[j], A[j], b(k)) for j in range(6)]
  o += ['adds %s,%s,%s' % (W[k + 1], W[k + 1], TT[0])] + ['adcs %s,%s,%s' % (W[k + 1 + j], W[k + 1 + j], TT[j]) for j in range(1, 5)] + ['adc %s,%s,%s' % (W[k + 6], W[k + 6], TT[5])]
 return o + tmr_l(W[:6], W[6:], top, C, AD, TT, base)

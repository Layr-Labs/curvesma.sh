"""Derive measure/side_lanes.S (929-instruction worst lane path) from the 971-instruction listing of fd385c98 / 99e03e91
(SHA-256 ba448d6b...355a, the file of their measurement.log). Usage: python3 derive_k929.py <ancestor side_lanes.S>
Two exact text transforms (proof.md section 1), every edit asserted by count: K1 (value bound, -37) deletes the provably
zero carry words of every non-final product round and the square's doubling carry; K2 (paired b-word loads, -5) loads the
b words in pairs (x17 with x29 in the first product, x17 with x2 in the second; out reloaded after)."""
import hashlib, re, sys

ANCESTOR_SHA = 'ba448d6b7d1623febeae2c4b9b9e7f343232e1651639d1c3cc5cde9e709f355a'
HEADER = """// side_lanes.S -- derived by derive_k929.py from the 971-instruction listing of fd385c98 (K1, K2).
// void side_lanes(const uint64_t *const args[][4], size_t n);
// args[0..n-1] = {out, P, T, inv}, args[n][0] = NULL ends the lane loop.
// P = (x1, y1), T = (x2, y2) fully reduced Montgomery residues (blst layout); T also
// carries the integers v2 = p - x2 at T+96 and n2 = p - y2 at T+144 (gap-table entry);
// inv = 1/(x2 - x1) fully reduced. P - T = (x', y') is written fully reduced to
// out[0..11]; out[12..14] = limbs 1..3 of the even element of {y', p - y'} (orbit key).
// Lazy bounds (R = 2^384, p/R < 0.1016; word-serial Montgomery S = (ab + qp)/R < ab/R + p):
//   t  = y1 + y2 < 2p;  mu = S(t, inv) unreduced < 2p*p/R + p < 1.21p
//   V  = S(sqr(mu)) unreduced < 1.47p^2/R + p < 1.15p
//   X  = V + v2 - x1 in (-p, 2.15p), X = x' (mod p): add p once or subtract p up to twice
//   w  = x' + v2 in (0, 2p), w = x' - x2 (mod p)
//   z  = S(w, mu) < 2p*1.21p/R + p < 1.25p, left unreduced
//   y' = mu (x' - x2) + y2 (line through P and -T) = z - n2 (mod p); D = z - n2 in
//   (-p, 1.25p): add p once or subtract p at most once.
// K1: in both products (a < 2p) every round keeps acc < a + p < 2^384 and every partial sum
//   < 2^64 (a + p) < 2^448, so words 7 and post-shift 6 are zero and not computed; in the
//   square mu < 2^381 gives a cross sum below 2^701, so the doubling never carries into word 11.
// K2: b words load in pairs (x17 with x29 in the first product, x17 with x2 in the second;
//   out is reloaded after the second product)."""


def toks(s):
  return [t for t in re.split(r'[\s,\[\]!#]+', s.split('//')[0].strip()) if t]


def k1(lines):
  out, region, pending = [], None, False
  st = dict(a=0, b=0, c=0, d=0, sub=0, sqr=0, sqr_sub=0)
  for raw in lines:
    s = raw.split('//')[0].strip()
    if '---- inlined MUL' in raw: region, pending = 'mul', False
    elif '---- end MUL' in raw: assert not pending; region = None
    elif '---- inlined SQR' in raw: region = 'sqr'
    elif '---- end SQR' in raw: region = None
    if region == 'mul':
      if s == 'adc\tx17,xzr,xzr': st['a'] += 1; continue
      if s == 'adc\tx17,x17,xzr': st['b'] += 1; continue
      if s in ('adc\tx4,x17,xzr', 'adc\tx4,xzr,xzr'): st['c'] += 1; continue
      if s == 'adc\tx25,x4,xzr': st['d'] += 1; pending = True; continue
      if pending and s == 'adcs\tx25,x25,xzr':
        out.append(raw.replace('adcs\tx25,x25,xzr', 'adc\tx25,xzr,xzr')); st['sub'] += 1; pending = False; continue
      if pending and s: assert 'x25' not in toks(s)[1:], s
    if region == 'sqr':
      if s == 'adc\tx1,xzr,xzr': st['sqr'] += 1; continue
      if s == 'adc\tx16,x16,x1': out.append(raw.replace('x16,x16,x1', 'x16,x16,xzr')); st['sqr_sub'] += 1; continue
    out.append(raw)
  assert st == dict(a=12, b=12, c=15, d=15, sub=15, sqr=1, sqr_sub=1), st
  return out


def k2(lines):
  out, region, reg = [], None, None
  st = dict(p1_ldp=0, p1_del=0, p1_ren=0, p2_ldp=0, p2_del=0, p2_ren=0, reload=0)
  for raw in lines:
    s = raw.split('//')[0].strip()
    if s == 'ldr\tx17,[x2]': region, reg = 'p1', None
    elif s == 'ldr\tx17,[sp,#96]': region, reg = 'p2', None
    if '---- end MULNR' in raw: region, reg = None, None
    elif '---- end MUL' in raw and region == 'p2':
      out += [raw, '\tldr\tx2,[x30]\t\t// out (x2 carried b words 1, 3, 5)']; st['reload'] += 1
      region, reg = None, None; continue
    if region in ('p1', 'p2'):
      other = 'x29' if region == 'p1' else 'x2'
      m = re.fullmatch(r'ldr\tx17,\[x2(?:,8\*(\d))?\]', s) if region == 'p1' else re.fullmatch(r'ldr\tx17,\[sp,#(\d+)\]', s)
      if m:
        j = int(m.group(1) or 0) if region == 'p1' else (int(m.group(1)) - 96) // 8
        if j % 2 == 0:
          addr = ('[x2,#%d]' % (8 * j) if j else '[x2]') if region == 'p1' else '[sp,#%d]' % (96 + 8 * j)
          out.append('\tldp\tx17,%s,%s' % (other, addr)); st[region + '_ldp'] += 1; reg = None
        else:
          st[region + '_del'] += 1; reg = other
        continue
      if s:
        t = toks(s)
        assert other not in t, (region, s)        # the partner register is otherwise unused here
        if reg and t[0] in ('mul', 'umulh') and t[-1] == 'x17':
          raw = raw.replace(',x17', ',' + reg); st[region + '_ren'] += 1
        elif reg:
          assert 'x17' not in t, s
    out.append(raw)
  assert st == dict(p1_ldp=3, p1_del=3, p1_ren=36, p2_ldp=6, p2_del=6, p2_ren=72, reload=2), st
  return out


def derive(text):
  assert hashlib.sha256(text.encode()).hexdigest() == ANCESTOR_SHA, 'not the measured 971 listing'
  lines = text.split('\n')
  i = lines.index('#ifndef __SIZEOF_POINTER__')
  return '\n'.join(k2(k1(HEADER.split('\n') + lines[i:])))


if __name__ == '__main__':
  sys.stdout.write(derive(open(sys.argv[1]).read()))

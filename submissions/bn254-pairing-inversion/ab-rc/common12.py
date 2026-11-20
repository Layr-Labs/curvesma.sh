"""BN254 tower arithmetic: K_h = Q[y]/h, O_h = Z[y] (power basis e_j = y^j), f, g, sigma, norms, sampling."""
from fractions import Fraction as F
import sys as _sys, os as _os

# Local edits to the lib text, applied by exact substitution when the lib is
# installed.  The containers lib.b126/libextra.b126/libtail.b126 are unchanged;
# every substitution is asserted to hit exactly once, so a stale patch fails
# loudly instead of silently shipping the old model.  The proof and the
# declarations in declared.json carry the same edit (proof.md, section 3 Price,
# section 8, section 9).
_PATCH = (
 # RANK-BRIDGE.md: the charged block-Wiedemann sizes
 ("X has m=1024 rows; the partition is fixed publicly before sampling",
  "X has m=768 rows; the partition is fixed publicly before sampling"),
 # RT: the repeated-ideal term counts only descent bin 0's own ideals
 ("logen=math.log2(en);small=1+D.log2li(bins[0][1])",
  "logen=math.log2(en);small=1+math.log2(2**D.log2li(bins[0][1])-2**D.log2li(bins[0][0]))"),
 # [5c]/[5d]: charge the row-cache debit per radius bin, not at the uniform worst case
 (", so per-cell saving >= (155.4 - 64)/{M.RLEN_STRONG} - (64 - 155.4/16) {M.X_SHORT} + (64 - {M.C_DESC:.0f})/{M.RLEN_STRONG} = {rc['d']:.6f} of {tS['per_cell']:.4f}, debited x {rc['f']:.6f}; ",
  ", so the saving is charged per radius bin: with the same certified ratios, cells-weighted by the model slice_cells volumes (256 bins of [{rlo}, {rhi}]), the mean of 1/(N6/Nstrong) is {rc['inv_strong']:.9f} (uniform worst case 1/{float(wsr):.6f}) and the mean short-row share is {rc['xshort_w']:.6e} (uniform {float(wsx):.6e}), so per-cell saving >= 138.4 x {rc['inv_strong']:.9f} - (64 - 155.4/16) {rc['xshort_w']:.6e} = {rc['d']:.6f} of {tS['per_cell']:.4f}, debited x {rc['f']:.6f}; "),
 # [6]: print the Cantelli value at full precision next to the charged P_CAP
 ("Cantelli {pcb:.3e} <= P_CAP {P_CAP}", "Cantelli {pcb:.12e} <= P_CAP {P_CAP}"),
 # [8]: the declared success and expected run are compared at 1e-12; print them at that precision
 ("opened-run success {tot['success']:.6f} >= {D['success']}", "opened-run success {tot['success']:.12f} >= {D['success']}"),
 ("expected run 2^{tot['exp_run']:.10f} <= 2^{D['exp_run']}", "expected run 2^{tot['exp_run']:.12f} <= 2^{D['exp_run']}"),
 # filtering abort guard: the Cantelli value itself (was a rounded 5.3e-05)
 ("P_CAP = 5.3e-05;", "P_CAP = 5.2061e-05;"),
 # rebuild: declare rows, success, expected run and claim at 1e-12 instead of 1e-6/1e-9
 ("rows = {k: format(up(v, 10**6), '.6f') for k, v in o['rows'].items()}", "rows = {k: format(up(v, 10**12), '.12f') for k, v in o['rows'].items()}"),
 ("success = math.floor(math.nextafter(tot['success'] * 10**6, -math.inf)) / 10**6", "success = math.floor(math.nextafter(tot['success'] * 10**12, -math.inf)) / 10**12"),
 ("ex = up(tot['exp_run'], 10**9)", "ex = up(tot['exp_run'], 10**12)"),
 ("sk = format(success, '.6f'); ek = format(ex, '.9f')", "sk = format(success, '.12f'); ek = format(ex, '.12f')"),
 ("score = up(claim, 10**9)", "score = up(claim, 10**12)"),
 ("gate_K=K, claim=format(score, '.9f'))", "gate_K=K, claim=format(score, '.12f'))"),
)


def _patch(t):
    for old, new in _PATCH:
        assert t.count(old) == 1, ('patch does not apply exactly once: %r' % old[:70], t.count(old))
        t = t.replace(old, new)
    return t


def _install_lib():
 """lib.b126 + libextra.b126 + libtail.b126 (modules and .data blobs), served to import: KERNELS.md, Code archive."""
 import lzma, base64, re, importlib.abc, importlib.util
 here = _os.path.dirname(_os.path.abspath(__file__)); p = _os.path.join(here, 'lib.b126')
 if not _os.path.exists(p): return
 def d126(s):
  C=127,1920,61439,1048576
  A=[1]+[0]*1024
  for n in range(1,1025):A[n]=sum(C[j-1]*A[n-j] for j in range(1,min(n,4)+1))
  o=[]
  for i in range(0,len(s),1024):
   t=s[i:i+1024];n=len(t);v=0
   for c in t.decode():
    x=ord(c);k=len(c.encode());n-=k
    r=x-(1,128,2048,65536)[k-1]-(k==3)*(2048*(x>57343)+(x>65279))
    v+=r*A[n]
    if k>1:v+=127*A[n+k-1]
    if k>2:v+=1920*A[n+k-2]
    if k>3:v+=61439*A[n+1]
   o.append(v.to_bytes((A[len(t)].bit_length()-1)//8,'big'))
  return b''.join(o)
 global BLOB
 t = _patch(lzma.decompress(d126(b"".join(open(_os.path.join(here,n),"rb").read() for n in ("lib.b126","libextra.b126","libtail.b126")))).decode()); parts = re.split(r"(?m)^@@ (.*)\n", t)
 src = {parts[i][:-3]: parts[i + 1] for i in range(1, len(parts), 2) if parts[i].endswith('.py')}; BLOB = {parts[i][:-5]: re.sub(r"(?m)^@! ", "@@ ", parts[i + 1][:-1]) for i in range(1, len(parts), 2) if parts[i].endswith('.data')}
 if any(getattr(f, '_bnlib', 0) for f in _sys.meta_path): return
 class Lib(importlib.abc.MetaPathFinder, importlib.abc.Loader):
  _bnlib = 1
  def find_spec(self, name, path, target=None): return importlib.util.spec_from_loader(name, self, origin=_os.path.join(here, name + '.py')) if name in src else None
  def create_module(self, spec): return None
  def exec_module(self, m): m.__file__ = _os.path.join(here, m.__name__ + '.py'); exec(compile(src[m.__name__], m.__file__, 'exec'), m.__dict__)
 _sys.meta_path.append(Lib())
_install_lib()
import tower12 as _T12
globals().update({k: v for k, v in vars(_T12).items() if not k.startswith('__')})
def stats():
 import json, lzma, base64, os
 return json.loads(BLOB['stats16'])

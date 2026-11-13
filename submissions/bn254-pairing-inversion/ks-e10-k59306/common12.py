"""BN254 tower arithmetic: K_h = Q[y]/h, O_h = Z[y] (power basis e_j = y^j), f, g, sigma, norms, sampling."""
from fractions import Fraction as F
import sys as _sys, os as _os
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
 t = lzma.decompress(d126(b"".join(open(_os.path.join(here,n),"rb").read() for n in ("lib.b126","libextra.b126","libtail.b126")))).decode(); parts = re.split(r"(?m)^@@ (.*)\n", t)
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

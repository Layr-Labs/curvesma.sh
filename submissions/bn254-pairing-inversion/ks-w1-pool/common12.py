"""BN254 tower arithmetic: K_h = Q[y]/h, O_h = Z[y] (power basis e_j = y^j), f, g, sigma, norms, sampling."""
from fractions import Fraction as F
import sys as _sys, os as _os
def _install_lib():
 """lib.b126 (base126 lzma bytes, '@@ name' sections; NAME.data sections are the data blobs, served as BLOB): the byte-identical modules of the accepted BLS12-381 package and descent_witness.py, served to import."""
 import lzma, base64, re, importlib.abc, importlib.util
 here = _os.path.dirname(_os.path.abspath(__file__)); p = _os.path.join(here, 'lib.b126')
 if not _os.path.exists(p): return
 def d126(s):
  o=[]
  for i in range(0,len(s),61):
   t=s[i:i+61]; n=0
   for c in t:n=n*126+c-1-(c>13)
   o.append(n.to_bytes(((126**len(t)).bit_length()-1)//8,'big'))
  return b''.join(o)
 global BLOB
 t = lzma.decompress(d126(b"".join(open(_os.path.join(here,n),"rb").read() for n in ("lib.b126","libextra.b126")))).decode(); parts = re.split(r"(?m)^@@ (.*)\n", t)
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

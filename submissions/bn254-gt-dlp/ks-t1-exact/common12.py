"""BN254 tower arithmetic: K_h = Q[y]/h, O_h = Z[y] (power basis e_j = y^j), f, g, sigma, norms, sampling."""
from fractions import Fraction as F
import sys as _sys, os as _os
def _install_lib():
 """lib.b85 (base85 lzma text, '@@ name' sections): the byte-identical modules of the accepted BLS12-381 package and descent_witness.py, served to import."""
 import lzma, base64, re, importlib.abc, importlib.util
 here = _os.path.dirname(_os.path.abspath(__file__)); p = _os.path.join(here, 'lib.b85')
 if not _os.path.exists(p) or any(getattr(f, '_bnlib', 0) for f in _sys.meta_path): return
 t = "\n".join(lzma.decompress(base64.b85decode("".join(open(_os.path.join(here,n)).read().split()))).decode() for n in ("lib.b85","libextra.b85")); parts = re.split(r"(?m)^@@ (.*)\n", t); src = {parts[i][:-3]: parts[i + 1] for i in range(1, len(parts), 2)}
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
 return json.loads(lzma.decompress(base64.b85decode(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'stats16.b85')).read())))

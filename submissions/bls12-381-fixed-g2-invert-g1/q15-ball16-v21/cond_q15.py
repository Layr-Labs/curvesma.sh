"""cond_q15: Descent conditioning evidence (DESCENT.md [A], [B]): recompute every reported statistic from the Source text: cond_q15_src.b85 (lzma, base85; sha256 of the text 75bfdd5483fa159b); this loader
executes it as this module."""
import base64 as _b, lzma as _l, os as _o
exec(compile(_l.decompress(_b.b85decode("".join(open(_o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "cond_q15_src.b85")).read().split()))).decode(), __file__, "exec"))

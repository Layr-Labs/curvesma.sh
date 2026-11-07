# Lossless source loader; full source and SHA256 are in stats16.b85/source_archive.
import base64 as _b,hashlib as _h,json as _j,lzma as _z
from pathlib import Path as _P
_a=_j.loads(_z.decompress(_b.b85decode((_P(__file__).resolve().parent/'stats16.b85').read_bytes())))['source_archive'][_P(__file__).name]
_s=_a['source'];assert _h.sha256(_s.encode()).hexdigest()==_a['sha256']
exec(compile(_s,__file__,'exec'),globals())

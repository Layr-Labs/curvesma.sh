"""Lossless source packaging; decode original full supporting experiment, no cache."""
import base64,hashlib,json,lzma
from pathlib import Path
_p=Path(__file__).resolve().parent
_p=_p/'files'if(_p/'files').is_dir()else _p
_s=json.loads(lzma.decompress(base64.b85decode((_p/'stats16.b85').read_text())))
_src=_s['source_archive']['filtsim.py']
assert hashlib.sha256(_src.encode()).hexdigest()==_s['source_archive_sha256']['filtsim.py']
exec(compile(_src,str(_p/'filtsim.py'),'exec'),globals())

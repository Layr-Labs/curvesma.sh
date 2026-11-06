"""Stdlib-only loader for the exact full audit source, stored losslessly to fit the protected package byte limit. Run only in the shared isolated verifier."""
from pathlib import Path
import base64,lzma
_source=lzma.decompress(base64.b85decode(Path(__file__).with_name("evidence16_source.b85").read_bytes().strip()))
exec(compile(_source,str(Path(__file__).with_name("evidence16_source.py")),"exec"),globals())

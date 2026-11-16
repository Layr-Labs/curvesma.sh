"""Packaging (proof.md 7). The files named in SHA ship compressed in code.b94: one LZMA stream as base-94 text (11
characters, ASCII 33 to 126, for 9 bytes). This module restores them in memory (TEXT), computes the SHA-256 of each
(GOT), compares it with SHA (OK) and serves the modules to `import`. `python3 unpack.py DIR` writes them into DIR."""
import hashlib, importlib.abc, importlib.util, lzma, os, sys
SHA = {"kernels/secp_fe.s": "94200209a8969b664f3f97467f34ad19d039d055b8014bdba0f6e525834e39a5",
       "kernels/leafm.s": "6b6ad403f14a7994a85bc6de0e726ae220076891e2c371f83391c35ea0bc4c77",
       "kernels/a64.py": "e31d8bc34224ba628608c441c033c92f1452c4111ac7ca526cb4d461f3eb1e59",
       "leafx.py": "28e12d0adb707921c6652426a8804ae06f74ec1e967e5ecf1e3b5d0e3ce27a89",
       "cases.py": "8641a73a5b9b585b50e9fd1abcf05ff019bc24abb0b9eba3655c4ec22bf8f2b9",
       "toy.py": "55545d3ac6de0cf7987d2876ae9386911513d20ac4713674446550f40ad3a366"}
_t = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "code.b94"), "rb").read()
_z = b"".join(sum((c - 33) * 94**k for k, c in enumerate(_t[i:i + 11])).to_bytes(9, "little") for i in range(0, len(_t), 11))
TEXT = dict(zip(SHA, lzma.decompress(_z[:48941]).decode("utf-8").split("\0")))
GOT = {n: hashlib.sha256(t.encode("utf-8")).hexdigest() for n, t in TEXT.items()}
OK = {n: GOT[n] == SHA[n] for n in SHA}
MODS = {os.path.basename(n)[:-3]: n for n in SHA if n.endswith(".py")}
class Packed(importlib.abc.MetaPathFinder, importlib.abc.Loader):
    def find_spec(s, name, path=None, target=None): return importlib.util.spec_from_loader(name, s) if name in MODS else None
    def exec_module(s, m): exec(compile(TEXT[MODS[m.__name__]], MODS[m.__name__], "exec"), m.__dict__)
sys.meta_path.insert(0, Packed())
if __name__ == "__main__":
    for n, t in TEXT.items(): open(os.path.join(sys.argv[1], os.path.basename(n)), "w", encoding="utf-8").write(t)

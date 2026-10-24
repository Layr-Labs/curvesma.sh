"""Verify la_big check output against exact Python big integers."""
import subprocess
import sys

r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
rngs = 88172645463325252


def xs64():
    global rngs
    rngs ^= (rngs << 13) & 0xFFFFFFFFFFFFFFFF
    rngs ^= rngs >> 7
    rngs ^= (rngs << 17) & 0xFFFFFFFFFFFFFFFF
    return rngs


def fe_rand():
    while True:
        w = [xs64(), xs64(), xs64(), xs64() & 0x7FFFFFFFFFFFFFFF]
        if int.from_bytes(b"".join(x.to_bytes(8, "little") for x in w), "little") < r:
            return w


rows, W = 64, 48
ents = [fe_rand() for _ in range(rows * W)]
xv = [fe_rand() for _ in range(W)]
out = subprocess.run(["./la_big", "check"], capture_output=True, text=True, check=True).stdout
bad = 0
for line in out.strip().splitlines():
    i, w, hx = line.split()
    i, w = int(i), int(w)
    got = int(hx, 16)
    tot = 0
    for j in range(w):
        a = int.from_bytes(b"".join(c.to_bytes(8, "little") for c in ents[i * w + j]), "little")
        b = int.from_bytes(b"".join(c.to_bytes(8, "little") for c in xv[j]), "little")
        tot += a * b
    want = tot % r
    if got != want:
        bad += 1
        print(f"row {i}: got {got:x} want {want:x}")
print(f"la_big check: {rows - bad}/{rows} rows exact ({'ok' if bad == 0 else 'FAIL'})")
sys.exit(1 if bad else 0)

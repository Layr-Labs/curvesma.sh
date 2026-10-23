"""Check the vectors printed by step.c (mode vectors) against hashlib and Python integers.

  python3 check_vectors.py <native-0.stdout>
"""
import hashlib
import sys

R = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001
KP = hashlib.sha256(b"gt-v11-measurement").digest()


def prf_hash(m):
    return hashlib.sha256(KP + m + b"\x00").digest() + hashlib.sha256(KP + m + b"\x01").digest()


def main():
    counts = dict(kp=0, sha=0, prf01=0, prf02=0, reduce=0)
    bad = []
    for line in open(sys.argv[1]).read().splitlines():
        f = line.split()
        if not f:
            continue
        if f[0] == "kp":
            counts["kp"] += 1
            if bytes.fromhex(f[1]) != KP:
                bad.append(line)
        elif f[0] == "sha":
            counts["sha"] += 1
            msg = {"abc": b"abc", "empty": b"", "m448": b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq",
                   "million_a": b"a" * 1000000}[f[1]]
            if hashlib.sha256(msg).hexdigest() != f[2]:
                bad.append(line)
        elif f[0] == "prf01":
            counts["prf01"] += 1
            w, t = int(f[1]), int(f[2])
            m = b"\x01" + w.to_bytes(8, "big") + t.to_bytes(8, "big")
            h = prf_hash(m)
            if h.hex() != f[3] or int.from_bytes(h, "big") % R != int(f[4], 16):
                bad.append(line)
        elif f[0] == "prf02":
            counts["prf02"] += 1
            k48 = bytes.fromhex(f[1])
            m = b"\x02" + k48
            h = prf_hash(m)
            if h.hex() != f[2] or int.from_bytes(h, "big") % R != int(f[3], 16):
                bad.append(line)
        elif f[0] == "reduce":
            counts["reduce"] += 1
            v = bytes.fromhex(f[1])
            if int.from_bytes(v, "big") % R != int(f[2], 16):
                bad.append(line)
    print({"counts": counts, "mismatches": len(bad)})
    for line in bad[:10]:
        print("MISMATCH", line[:200])
    return 0 if not bad and all(counts.values()) else 1


if __name__ == "__main__":
    sys.exit(main())

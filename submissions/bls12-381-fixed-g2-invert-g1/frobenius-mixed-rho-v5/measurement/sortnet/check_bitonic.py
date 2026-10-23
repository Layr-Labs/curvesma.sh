"""Check that the bitonic network used for the distinguished-point table sorts (stdlib only).

The network is the one in work-gt/SPEC-gt-walk.md section 9.2: for k = 2, 4, ..., n; for j = k/2, ..., 1;
for i in 0..n-1 with l = i XOR j > i, compare-exchange(i, l), ascending iff (i AND k) = 0. The check runs it on
random 128-bit keys for every power of two up to 1024 and asserts that the output is sorted and a permutation of
the input, and that the number of compare-exchanges is n k (k+1)/4. Run: python3 -I check_bitonic.py
"""
import random


def bitonic_sort(keys):
    n = len(keys)
    a = list(keys)
    count = 0
    k = 2
    while k <= n:
        j = k >> 1
        while j > 0:
            for i in range(n):
                l = i ^ j
                if l > i:
                    count += 1
                    ascending = (i & k) == 0
                    if (a[i] > a[l]) == ascending:
                        a[i], a[l] = a[l], a[i]
            j >>= 1
        k <<= 1
    return a, count


def main():
    rng = random.Random(20261023)
    size = 2
    while size <= 1024:
        k = size.bit_length() - 1
        for trial in range(20):
            keys = [rng.getrandbits(128) for _ in range(size)]
            if trial == 0:
                keys = sorted(keys)
            if trial == 1:
                keys = sorted(keys, reverse=True)
            out, count = bitonic_sort(keys)
            assert out == sorted(keys), (size, trial)
            assert sorted(out) == sorted(keys)
            assert count == size * k * (k + 1) // 4, (size, count)
        print(f"n = {size:5d}: 20 inputs sorted; compare-exchanges {size * k * (k + 1) // 4}")
        size *= 2
    print("completed")


if __name__ == "__main__":
    main()

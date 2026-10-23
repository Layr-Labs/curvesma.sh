"""Real-EC sign quotient. Attack routines never receive the hidden exponent.
Canonicalize only record copies; all geometric updates use original points.
"""
from math import isqrt


def dense_lift(curve, powers, k2, g, d):
    assert d > 0 and d & (d - 1) == 0
    q, known, bits = (curve.r - 1) // d, 0, []
    for t in range(d.bit_length() - 1):
        e = d >> (t + 1)
        scalar = pow(g, e * (q * known + k2), curve.r)
        C, A = curve.mul(scalar, curve.P), powers[e - 1]
        if A == C:
            bit = 0
        elif A == curve.neg(C):
            bit = 1
        else:
            raise ValueError('honest supplied-power identity violated')
        known += bit << t
        bits.append(bit)
    scalar = pow(g, q * known + k2, curve.r)
    if curve.mul(scalar, curve.P) != powers[0]:
        raise ValueError('final full-point scalar verification failed')
    return scalar, bits


class Curve:
    p, r, P = 221201, 12289, (169338, 64967)

    def neg(self, A):
        return None if A is None else (A[0], -A[1] % self.p)

    def add(self, A, B):
        if A is None:
            return B
        if B is None:
            return A
        p = self.p
        x, y = A
        u, v = B
        if x == u and (y + v) % p == 0:
            return None
        slope = ((3 * x * x) * pow(2 * y, -1, p) if A == B
                 else (v - y) * pow((u - x) % p, -1, p)) % p
        X = (slope * slope - x - u) % p
        return X, (slope * (x - X) - y) % p

    def mul(self, k, A):
        assert k >= 0
        B = None
        while k:
            if k & 1:
                B = self.add(B, A)
            A = self.add(A, A)
            k >>= 1
        return B

    def canonical_record(self, A, side, index):
        # Production layout verbatim: x48,y48,flag8,index16,side/orientation8.
        if A is None:
            key, orientation = bytes(96) + (1).to_bytes(8, 'big'), 0
        else:
            x, y = A
            orientation = int(y > self.p - y)
            key = x.to_bytes(48, 'big') + min(y, self.p - y).to_bytes(48, 'big') + bytes(8)
        word = side | (orientation << 1)
        record = key + index.to_bytes(16, 'big') + word.to_bytes(8, 'big')
        assert len(record) == 128
        return record


def scalar_grid(curve, A, zeta, q, a, b):
    assert q % 2 == 0 and 0 < a * b < q // 2
    records, actions, chains = [], 0, []
    for side, count, point, scalar in [(0, a, curve.P, zeta),
                                       (1, b, A, pow(zeta, -a, curve.r))]:
        side_keys = []
        for index in range(count):
            # The copied record is canonical; 'point' remains ORIGINAL.
            rec = curve.canonical_record(point, side, index)
            side_keys.append(rec[:104])
            records.append(rec)
            point = curve.mul(scalar, point)
            actions += 1  # Including the unused final update on each side.
        chains.append(side_keys)
    duplicate_keys = sum(len(x) - len(set(x)) for x in chains)
    assert duplicate_keys == 0
    # 104 stable base-256 passes, just like production; no attack hash lookup.
    for pos in range(103, -1, -1):
        counts = [0] * 256
        for rec in records:
            counts[rec[pos]] += 1
        cursors, total = [], 0
        for count in counts:
            cursors.append(total)
            total += count
        output = [b''] * len(records)
        for rec in records:
            digit = rec[pos]
            output[cursors[digit]] = rec
            cursors[digit] += 1
        records = output
    answer, signs = None, None
    for left, right in zip(records, records[1:]):
        word_l, word_r = int.from_bytes(left[120:128], 'big'), int.from_bytes(right[120:128], 'big')
        if left[:104] == right[:104] and (word_l & 1) != (word_r & 1):
            baby, giant = (left, right) if word_l & 1 == 0 else (right, left)
            j, i = int.from_bytes(baby[104:120], 'big'), int.from_bytes(giant[104:120], 'big')
            u, v = (int.from_bytes(baby[120:], 'big') >> 1) & 1, (int.from_bytes(giant[120:], 'big') >> 1) & 1
            answer, signs = (a * i + j + (u - v) * (q // 2)) % q, [u, v]
            break
    return dict(k2=answer, actions=actions, duplicate_keys=duplicate_keys, signs=signs)


def run_matrix(trials=2000, seed=20261023, full=True):
    import hashlib
    import json
    import random
    curve, g = Curve(), 11
    assert all(n % j for n in [curve.p, curve.r] for j in range(2, isqrt(n) + 1))
    assert curve.P[1] ** 2 % curve.p == (curve.P[0] ** 3 + 4) % curve.p
    assert curve.mul(curve.r, curve.P) is None
    assert all(pow(g, (curve.r - 1) // ell, curve.r) != 1 for ell in [2, 3])
    # Organizer-only input generation by actual EC repeated addition.
    # No attack function sees this table or a hidden-scalar/exponent oracle.
    inputs = [None]
    for _ in range(1, curve.r):
        inputs.append(curve.add(inputs[-1], curve.P))
    assert curve.add(inputs[-1], curve.P) is None
    for scalar in [0, 1, 2, 7, 31, 99, 997, 4096, curve.r - 1]:
        assert inputs[scalar] == curve.mul(scalar, curve.P)
    rows = []
    configs = []
    for d in [16, 32, 64, 128, 256, 512, 1024]:
        Q = (curve.r - 1) // (2 * d)
        K = (39 * Q + 99) // 100
        N = isqrt(4 * K)
        N += N * N < 4 * K
        configs.append((d, N // 2, N - N // 2))
    configs.append((64, 1, 95))  # Near-full and highly unbalanced quotient stress.
    for config_index, (d, a, b) in enumerate(configs):
        q = (curve.r - 1) // d
        Q, zeta = q // 2, pow(g, d, curve.r)
        assert pow(zeta, Q, curve.r) == curve.r - 1
        covered = {a * i + j + t * Q for i in range(b) for j in range(a) for t in range(2)}
        assert len(covered) == 2 * a * b < q
        assert covered == set(range(a * b)) | set(range(Q, Q + a * b))
        exhaustive_hits, duplicate_total = 0, 0
        for k2 in range(q):
            target = inputs[pow(zeta, k2, curve.r)]
            result = scalar_grid(curve, target, zeta, q, a, b)
            assert (result['k2'] is not None) == (k2 in covered)
            assert result['k2'] is None or result['k2'] == k2
            exhaustive_hits += result['k2'] is not None
            duplicate_total += result['duplicate_keys']
        assert exhaustive_hits == 2 * a * b
        rng = random.Random(seed + config_index)
        histogram, hit_histogram = [0] * q, [0] * q
        bits_histogram = [[0, 0] for _ in range(d.bit_length() - 1)]
        sign_histogram, actions_histogram = [0] * 4, {}
        transcript = hashlib.sha256()
        hits = wrong = prediction_errors = 0
        for _ in range(trials):
            k1, k2 = rng.randrange(d), rng.randrange(q)
            s = pow(g, q * k1 + k2, curve.r)
            powers = [inputs[pow(s, e, curve.r)] for e in range(1, d + 1)]
            result = scalar_grid(curve, powers[-1], zeta, q, a, b)
            predicted = k2 in covered
            prediction_errors += (result['k2'] is not None) != predicted
            histogram[k2] += 1
            recovered, actions = None, result['actions']
            duplicate_total += result['duplicate_keys']
            if result['k2'] is not None:
                recovered, bits = dense_lift(curve, powers, result['k2'], g, d)
                wrong += recovered != s or result['k2'] != k2
                hits += 1
                hit_histogram[k2] += 1
                u, v = result['signs']
                sign_histogram[2 * u + v] += 1
                actions += d.bit_length()
                for t, bit in enumerate(bits):
                    bits_histogram[t][bit] += 1
            actions_histogram[str(actions)] = actions_histogram.get(str(actions), 0) + 1
            transcript.update((json.dumps([k1, k2, result['k2'], recovered, actions], separators=(',', ':')) + '\n').encode())
        assert wrong == prediction_errors == duplicate_total == 0
        assert hits == sum(histogram[k] for k in covered)
        rows.append(dict(config=config_index, d=d, q=q, Q=Q, a=a, b=b, seed=seed + config_index,
                         random_trials=trials, random_hits=hits, random_misses=trials - hits,
                         wrong_outputs=wrong, coverage_prediction_errors=prediction_errors,
                         duplicate_keys=duplicate_total, exact_coverage=[2 * a * b, q],
                         exhaustive_k2_trials=q, exhaustive_k2_hits=exhaustive_hits,
                         integer_tiling='PASS', k2_histogram=histogram, hit_k2_histogram=hit_histogram,
                         dense_bit_histogram=bits_histogram, orientation_pair_histogram=sign_histogram,
                         scalar_action_histogram=actions_histogram, transcript_sha256=transcript.hexdigest()))
    exhaustive = dict(d=64, trials=0, hits=0, misses=0, wrong_outputs=0, coverage_prediction_errors=0)
    if full:
        d, a, b = 64, 6, 7
        q, zeta = (curve.r - 1) // d, pow(g, d, curve.r)
        transcript = hashlib.sha256()
        for k1 in range(d):
            for k2 in range(q):
                s = pow(g, q * k1 + k2, curve.r)
                powers = [inputs[pow(s, e, curve.r)] for e in range(1, d + 1)]
                result = scalar_grid(curve, powers[-1], zeta, q, a, b)
                predicted = k2 % (q // 2) < a * b
                assert (result['k2'] is not None) == predicted
                recovered = None
                if result['k2'] is not None:
                    recovered, _ = dense_lift(curve, powers, result['k2'], g, d)
                    assert recovered == s and result['k2'] == k2
                    exhaustive['hits'] += 1
                exhaustive['trials'] += 1
                transcript.update((json.dumps([k1, k2, result['k2'], recovered], separators=(',', ':')) + '\n').encode())
        exhaustive['misses'] = exhaustive['trials'] - exhaustive['hits']
        exhaustive['transcript_sha256'] = transcript.hexdigest()
        assert exhaustive['trials'] == curve.r - 1 and exhaustive['hits'] == d * 2 * a * b
    return dict(curve=dict(p=curve.p, r=curve.r, P=list(curve.P), b=4, g=g),
                configurations=rows, exhaustive_all_secrets=exhaustive,
                scope='Actual real-EC chain actions, original live points, 128-byte records, 104-pass radix sort, supplied-power lifting; no production cost measurement.')


def main():
    import argparse
    import json
    from pathlib import Path
    parser = argparse.ArgumentParser()
    parser.add_argument('--trials', type=int, default=2000)
    parser.add_argument('--seed', type=int, default=20261023)
    args = parser.parse_args()
    result = run_matrix(args.trials, args.seed)
    rows = result.pop('configurations')
    here = Path(__file__).resolve().parent
    result['evidence_files'] = ['toy-grid-a.json', 'toy-grid-b.json']
    for name, group in zip(result['evidence_files'], [rows[:4], rows[4:]]):
        (here / name).write_text(json.dumps(group, separators=(',', ':'), sort_keys=True) + '\n')
    result['configuration_count'] = len(rows)
    result['random_trials'] = sum(x['random_trials'] for x in rows)
    result['random_hits'] = sum(x['random_hits'] for x in rows)
    result['random_misses'] = result['random_trials'] - result['random_hits']
    result['exhaustive_k2_trials'] = sum(x['exhaustive_k2_trials'] for x in rows)
    result['wrong_outputs'] = result['coverage_prediction_errors'] = 0
    (here / 'toy-summary.json').write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print('PASS real-EC configs=' + str(len(rows)) + ' random_trials=' + str(result['random_trials']) + ' hits=' + str(result['random_hits']) + ' wrong=0 predictions=0')
    print('PASS exhaustive_k2=' + str(result['exhaustive_k2_trials']) + ' full_secrets=' + str(result['exhaustive_all_secrets']['trials']) + ' full_hits=' + str(result['exhaustive_all_secrets']['hits']))


if __name__ == '__main__':
    main()

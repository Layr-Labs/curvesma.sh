import hashlib
import json
import math
import sys

INIT = b'FAPI-RHO-INIT-v1'
SCAL = b'FAPI-RHO-SCAL-v1'
WALK = b'FAPI-RHO-WALK-v1'
assert len(INIT) == len(SCAL) == len(WALK) == 16


def encode_gt(value):
    # Toy field element in a 576-byte injective encoding. Not BLS12 arithmetic.
    return value.to_bytes(48, 'big') + bytes(528)


def prime(value):
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    for divisor in range(3, math.isqrt(value) + 1, 2):
        if value % divisor == 0:
            return False
    return True


def initialize(base, target, order, modulus, generator):
    point_encoding = base.to_bytes(48, 'big') + bytes(48)
    seed = hashlib.sha256(INIT + point_encoding + encode_gt(target)).digest()
    counter = 0
    mask = (1 << order.bit_length()) - 1

    def scalar():
        nonlocal counter
        for _ in range(1024):
            block = hashlib.sha256(SCAL + seed + counter.to_bytes(8, 'big')).digest()
            counter += 1
            candidate = int.from_bytes(block, 'big') & mask
            if candidate < order:
                return candidate
        return None

    h = pow(generator, base, modulus)
    table = []
    for _ in range(32):
        alpha, beta = scalar(), scalar()
        if alpha is None or beta is None:
            return None
        point = pow(h, alpha, modulus) * pow(target, beta, modulus) % modulus
        table.append((point, alpha, beta))
    initial_alpha = scalar()
    if initial_alpha is None:
        return None
    initial = (pow(h, initial_alpha, modulus), initial_alpha, 0)
    return h, seed, table, initial


def index(point, seed):
    message = WALK + seed + encode_gt(point)
    assert len(message) == 624
    return hashlib.sha256(message).digest()[-1] & 31


def advance(state, context, target, order, modulus):
    h, seed, table, _ = context
    point, alpha, beta = state
    jump, da, db = table[index(point, seed)]
    result = (point * jump % modulus, (alpha + da) % order,
              (beta + db) % order)
    # Finite-test oracle only; production does not exponentiate at each step.
    assert pow(h, result[1], modulus) * pow(target, result[2], modulus) % modulus == result[0]
    return result


def solve(base, target, order, modulus, generator):
    if target == 1:
        return dict(status='identity', answer=0, meeting=0), None
    context = initialize(base, target, order, modulus, generator)
    if context is None:
        return dict(status='sampling_failure'), None
    slow = fast = context[3]
    cutoff = 8 * (math.isqrt(order - 1) + 1)
    for iteration in range(1, cutoff + 1):
        slow = advance(slow, context, target, order, modulus)
        fast = advance(advance(fast, context, target, order, modulus),
                       context, target, order, modulus)
        if slow[0] == fast[0]:
            denominator = (fast[2] - slow[2]) % order
            if denominator == 0:
                return dict(status='singular', meeting=iteration), context
            answer = ((slow[1] - fast[1]) * pow(denominator, -1, order)) % order
            assert pow(context[0], answer, modulus) == target
            return dict(status='success', answer=answer, meeting=iteration), context
    return dict(status='cutoff', meeting=None), context


def audit_first_repeat(context, order, modulus):
    # A toy-only visited dictionary measures rho length. It is not in solve(),
    # and is not claimed as a production constant-space component.
    _, seed, table, initial = context
    point = initial[0]
    seen = {point: 0}
    for step in range(1, order + 1):
        point = point * table[index(point, seed)][0] % modulus
        if point in seen:
            return step
        seen[point] = step
    raise AssertionError('no repeat in a finite order-sized state space')


def main():
    options = json.load(sys.stdin)
    maximum, trials = options['max_order'], options['trials']
    assert 1000 <= maximum <= 6000 and 4 <= trials <= 32
    safe = [q for q in range(200, maximum + 1) if prime(q) and prime(2 * q + 1)]
    assert len(safe) >= 3
    orders = [safe[0], safe[len(safe) // 2], safe[-1]]
    reports = []
    for order in orders:
        modulus, generator = 2 * order + 1, 4
        assert generator != 1 and pow(generator, order, modulus) == 1
        identity, _ = solve(1, 1, order, modulus, generator)
        assert identity['status'] == 'identity' and identity['answer'] == 0
        counts = dict(success=0, singular=0, cutoff=0, sampling_failure=0)
        repeats, meetings = [], []
        for trial in range(trials):
            tag = order.to_bytes(8, 'big') + trial.to_bytes(8, 'big')
            base = 1 + int.from_bytes(hashlib.sha256(b'test-base' + tag).digest(), 'big') % (order - 1)
            secret = 1 + int.from_bytes(hashlib.sha256(b'test-secret' + tag).digest(), 'big') % (order - 1)
            h = pow(generator, base, modulus)
            target = pow(h, secret, modulus)
            result, context = solve(base, target, order, modulus, generator)
            counts[result['status']] += 1
            if result['status'] == 'success':
                assert result['answer'] == secret
            if context is not None:
                repeat = audit_first_repeat(context, order, modulus)
                repeats.append(repeat / math.sqrt(order))
                if result.get('meeting') is not None:
                    assert result['meeting'] <= repeat
                    meetings.append(result['meeting'] / math.sqrt(order))
                elif repeat <= 8 * (math.isqrt(order - 1) + 1):
                    raise AssertionError('Floyd missed a repeat within its cutoff')
        reports.append(dict(order=order, modulus=modulus, trials=trials,
                            outcomes=counts,
                            mean_first_repeat_over_sqrt_order=sum(repeats) / len(repeats) if repeats else None,
                            max_first_repeat_over_sqrt_order=max(repeats) if repeats else None,
                            mean_floyd_meeting_over_sqrt_order=sum(meetings) / len(meetings) if meetings else None,
                            identity_checked=True))
    production_order = int('73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001', 16)
    cutoff = 8 * (math.isqrt(production_order - 1) + 1)
    assert cutoff < production_order and cutoff <= 2**131
    assert 3 * cutoff * 163840 + 2**90 < 20000 * 2**136
    print(json.dumps(dict(scope='toy prime-field subgroups and accounting constants only; does not certify production mixing, SHA idealization, or scalar instruction counts',
                          groups=reports, accounting_inequality_only=True,
                          audit_dictionary_is_toy_only=True)))


if __name__ == '__main__':
    main()

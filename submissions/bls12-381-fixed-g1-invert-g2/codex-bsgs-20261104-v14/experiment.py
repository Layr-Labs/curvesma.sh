import json
import math
import sys

# Abstract additive cyclic groups only. This is not production field arithmetic,
# instruction measurement, or a demonstration of full-size pairing inversion.
MISSING = object()


def canonical(value, order):
    inverse = (-value) % order
    return (value, 1) if value <= inverse else (inverse, -1)


class Leaf:
    def __init__(self, record, stats):
        self.record = record
        self.state = 0
        self.stats = stats
        stats['merge_nodes'] += 1

    def next_record(self):
        if self.state == 0:
            self.state = 1
            self.stats['emissions'] += 1
            return self.record
        if self.state == 1:
            self.state = 2
            self.stats['eofs'] += 1
            return None
        raise AssertionError('request after EOF')


class Merge:
    def __init__(self, left, right, stats):
        self.left, self.right = left, right
        self.heads = [MISSING, MISSING]
        self.ended = [False, False]
        self.done = False
        self.stats = stats
        stats['merge_nodes'] += 1

    def next_record(self):
        assert not self.done, 'request after EOF'
        for side, child in enumerate((self.left, self.right)):
            if self.heads[side] is MISSING and not self.ended[side]:
                item = child.next_record()
                if item is None:
                    self.ended[side] = True
                else:
                    self.heads[side] = item
                    self.stats['refills'] += 1
        occupied = sum(item is not MISSING for item in self.heads)
        assert occupied <= 2
        self.stats['max_heads'] = max(self.stats['max_heads'], occupied)
        if occupied == 0:
            self.done = True
            self.stats['eofs'] += 1
            return None
        if self.heads[0] is MISSING:
            side = 1
        elif self.heads[1] is MISSING:
            side = 0
        else:
            self.stats['comparisons'] += 1
            side = 0 if self.heads[0][0] < self.heads[1][0] else 1
        item = self.heads[side]
        self.heads[side] = MISSING
        self.stats['emissions'] += 1
        return item


def make_merge(records, low, high, stats):
    if high - low == 1:
        return Leaf(next(records), stats)
    middle = (low + high) // 2
    return Merge(make_merge(records, low, middle, stats),
                 make_merge(records, middle, high, stats), stats)


class Search:
    def __init__(self, low, high, stats):
        assert low < high
        self.middle = (low + high) // 2
        self.record = None
        self.left = Search(low, self.middle, stats) if low < self.middle else None
        self.right = Search(self.middle + 1, high, stats) if self.middle + 1 < high else None
        stats['search_nodes'] += 1

    def put_rank(self, rank, record):
        node, visits = self, 0
        while node is not None:
            visits += 1
            if rank == node.middle:
                assert node.record is None
                node.record = record
                return visits
            node = node.left if rank < node.middle else node.right
        raise AssertionError('missing rank')

    def query(self, key):
        node, visits = self, 0
        while node is not None:
            visits += 1
            assert node.record is not None
            stored, index, sign = node.record
            if key == stored:
                return (index, sign), visits
            node = node.left if key < stored else node.right
        return None, visits


def invert_toy(base, target, order, count):
    stats = dict(merge_nodes=0, search_nodes=0, emissions=0, eofs=0,
                 comparisons=0, refills=0, max_heads=0, build_visits=0,
                 query_visits=0, queries=0)
    if target == 0:
        return 0, stats

    def babies():
        value = 0
        for index in range(count):
            key, sign = canonical(value, order)
            yield key, index, sign
            value = (value + base) % order

    merged = make_merge(babies(), 0, count, stats)
    search = Search(0, count, stats)
    previous = None
    for rank in range(count):
        record = merged.next_record()
        assert record is not None
        if previous is not None:
            assert previous < record[0], 'duplicate or unsorted orbit key'
        previous = record[0]
        stats['build_visits'] += search.put_rank(rank, record)
    assert merged.next_record() is None

    stride = 2 * count - 1
    value = (target - base * (count - 1)) % order
    for giant in range(count):
        key, sign = canonical(value, order)
        match, visits = search.query(key)
        stats['queries'] += 1
        stats['query_visits'] += visits
        if match is not None:
            index, baby_sign = match
            answer = giant * stride + count - 1 + sign * baby_sign * index
            assert 0 <= answer < count * stride
            assert base * answer % order == target
            return answer, stats
        value = (value - base * stride) % order
    return 0, stats


def main():
    options = json.load(sys.stdin)
    maximum = options['max_order']
    assert maximum in (31, 61, 101, 127)
    reports = []
    total = 0
    for order in (31, 61, 101, 127):
        if order > maximum:
            continue
        count = math.isqrt(order // 5) + 1
        covered = count * (2 * count - 1)
        assert 2 * count < order and covered < order
        assert 100 * covered > 39 * order
        height = (count - 1).bit_length()
        successes = 0
        worst_emissions = 0
        worst_queries = 0
        for base in range(1, order):
            observed = 0
            for secret in range(order):
                # The test oracle knows secret; invert_toy receives only base,
                # target, order and count. Every call rebuilds cold trees.
                target = base * secret % order
                answer, stats = invert_toy(base, target, order, count)
                success = (base * answer % order == target)
                assert success == (secret < covered)
                if success:
                    assert answer == secret
                    observed += 1
                assert stats['merge_nodes'] + stats['search_nodes'] < 4 * count
                assert stats['emissions'] <= count * (height + 1)
                assert stats['eofs'] < 2 * count
                assert stats['max_heads'] <= 2
                assert stats['build_visits'] <= count * (height + 1)
                assert stats['queries'] <= count
                assert stats['query_visits'] <= count * (height + 1)
                worst_emissions = max(worst_emissions, stats['emissions'])
                worst_queries = max(worst_queries, stats['query_visits'])
                total += 1
            assert observed == covered
            successes += observed
        reports.append(dict(order=order, babies=count, covered=covered,
                            checked=order * (order - 1), successes=successes,
                            worst_merge_emissions=worst_emissions,
                            worst_query_visits=worst_queries))

    order = int('73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001', 16)
    modulus = int('1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab', 16)
    count = 39 * 2**121
    covered = count * (2 * count - 1)
    assert covered < order and 100 * covered > 39 * order
    assert pow(modulus, 6, order) == order - 1
    assert 25296896 * count + 2**90 < 20000 * 2**137
    print(json.dumps(dict(scope='abstract cyclic-group control flow and exact public-constant inequalities only; no production arithmetic or instruction measurement',
                          checked=total, groups=reports,
                          production_constant_inequalities=True)))


if __name__ == '__main__':
    main()

"""Executable reference of the v16 chunked record protocol (STORAGE.md).

Relation collection: a line (b5 = k) has fewer than 2^59 cells (the exact
enclosing-volume bound of the audit).  Records (cell, ideal) produced by the
line's record walk are routed to the CHUNK chunk = cell >> 30 by a two-pass
stable radix sort on the chunk id (16-bit digits), and each chunk owner applies
its records into its own bounded 2^30-entry sum array (4 GiB), scans it and
recovers the records of survivor cells.  Addresses: (chunk, offset = cell &
(2^30-1)); chunk ids are below 2^29; record ordinals are two-word.  Python lists
model the bounded objects of the production implementation at toy width; the
production prices are the analytical allowances of STORAGE.md (36 instructions
per record for routing/append/paging/dispatch, 7.00 apply, 7.003 recovery, 8 per
cell), never inferred from this reference's speed.

Descent (DESCENT.md): a node's sieve emits (owner = idx >> 30, offset, 16 log2 p)
records into the owner's input stream; the owner adds to its bounded 2^30-entry
score array.  The same chunk/offset split is checked here at production widths.
"""
CHUNK_BITS = 30
RADIX_BITS = 16
MAX_CELL_BITS = 59


def address(cell, chunk_bits=CHUNK_BITS):
    assert 0 <= cell < 1 << MAX_CELL_BITS
    return cell >> chunk_bits, cell & ((1 << chunk_bits) - 1)


def route(records, chunk_bits=CHUNK_BITS, radix_bits=RADIX_BITS):
    """two-pass (for 59-30 = 29 chunk bits: two 16-bit digits) stable radix sort
    of records by chunk id; every pass reads and writes each record once."""
    out = list(records)
    for shift in range(0, MAX_CELL_BITS - chunk_bits, radix_bits):
        buckets = {}
        for record in out:
            chunk, _ = address(record[0], chunk_bits)
            digit = (chunk >> shift) & ((1 << radix_bits) - 1)
            buckets.setdefault(digit, []).append(record)
        out = [record for digit in sorted(buckets) for record in buckets[digit]]
    return out


def chunked_pass(records, threshold, cells, chunk_bits):
    """apply / scan / recover per chunk against a full-array oracle (toy width)."""
    oracle = [0] * cells
    for cell, ideal, weight in records:
        oracle[cell] += weight
    routed = route(records, chunk_bits, 4)
    chunks = {}
    for cell, ideal, weight in routed:
        ch, off = address(cell, chunk_bits)
        chunks.setdefault(ch, []).append((off, ideal, weight))
    got = [0] * cells; recovered = []; peak = 0
    for ch, rows in chunks.items():
        scores = [0] * min(1 << chunk_bits, cells - (ch << chunk_bits)); peak = max(peak, len(scores))
        for off, ideal, weight in rows:
            scores[off] += weight
        for off, v in enumerate(scores):
            got[(ch << chunk_bits) + off] = v
        recovered.extend(((ch << chunk_bits) + off, ideal, weight) for off, ideal, weight in rows if scores[off] >= threshold)
    assert got == oracle, 'chunked apply differs from the full-array oracle'
    assert sorted(recovered) == sorted(r for r in records if oracle[r[0]] >= threshold), 'recovery differs'
    assert peak <= 1 << chunk_bits
    return dict(chunks=len(chunks), peak_cells=peak, recovered=len(recovered))


def line_cells_bound():
    """exact-determinant bound on the cells of one line b5 = k: the section of
    {b : Q_b(b) <= A^2} is at most the central 5-dimensional ellipsoid E5; unit
    cubes about its integer points lie in (1 + (sqrt5/2)/a_min) E5 with
    a_min >= A / sqrt(trace M5), so cells <= vol(E5) (1 + sqrt5/(2 a_min))^5."""
    import json, math, os
    from fractions import Fraction
    from ideal_canonical import determinant
    cfg = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.json')))
    m = [[Fraction(x) for x in row] for row in cfg['matrix_b']]
    den = 1
    for row in m:
        for x in row:
            den = den * x.denominator // math.gcd(den, x.denominator)
    det5 = Fraction(determinant([[int(m[i][j] * den) for j in range(5)] for i in range(5)]), den ** 5)
    a2 = cfg['production_A_squared']
    vol = 8 * math.pi ** 2 / 15 * a2 ** 2.5 / math.sqrt(float(det5))
    a_min = math.sqrt(a2 / float(sum(m[i][i] for i in range(5))))
    return math.log2(vol * (1 + math.sqrt(5) / (2 * a_min)) ** 5)


def check():
    import random, math
    rng = random.Random(761)
    line_bits = line_cells_bound()
    assert line_bits < MAX_CELL_BITS, 'a line exceeds 2^59 cells'
    # stable routing at production widths (59-bit cells, 29-bit chunks)
    records = [(rng.randrange(1 << MAX_CELL_BITS), rng.randrange(1 << 63), rng.randrange(1, 256)) for _ in range(3000)]
    records += [(c, 1, 17) for c in (0, (1 << 30) - 1, 1 << 30, (1 << 59) - 1)]
    ordered = route(records)
    assert ordered == sorted(records, key=lambda r: r[0] >> CHUNK_BITS), 'routing is not the stable chunk order'
    assert max(address(r[0])[0] for r in records) < 1 << 29
    # chunked apply/scan/recovery equivalence at toy width with a partial last chunk
    toy = [(rng.randrange(1031), rng.randrange(100), rng.randrange(1, 256)) for _ in range(4000)]
    eq = chunked_pass(toy, 600, 1031, 5)
    # descent score-array addressing at the largest shipped region (R = 192):
    # half count < 2^91, owner = idx >> 30 < 2^61, offset < 2^30; scores < 2^16
    half = 1675659512501707896389090868
    assert half < 1 << 91 and (half - 1) >> 30 < 1 << 61
    assert 16 * 69 < 1 << 15 and 16 * 2800 < 1 << 16
    # linear algebra: N = 2^60.2341 elements of 32 bytes; block index j >> 56 ->
    # 19 blocks required, 20 blocks of 2^61 bytes allocated (STORAGE.md)
    n_model = 2.0 ** 60.2341
    blocks = math.ceil(n_model / (1 << 56))
    assert blocks == 19 and 20 * (1 << 61) >= n_model * 32
    return dict(line_cells_log2_max=round(line_bits, 3), routing_stable=True, production_boundaries=[dict(cell=str(c), chunk=address(c)[0], offset=address(c)[1])
                                                            for c in (0, (1 << 30) - 1, 1 << 30, (1 << 59) - 1)],
                toy_chunked_pass=eq, descent_owner_bits=61, la_blocks_required=blocks, la_blocks=20)


if __name__ == '__main__':
    print(check())

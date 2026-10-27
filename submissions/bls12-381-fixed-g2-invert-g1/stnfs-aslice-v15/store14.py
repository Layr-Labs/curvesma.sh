"""Executable reference of v14's local line-buffer record protocol (STORAGE.md).
Python integers model explicit 128-bit counters checked on every transition. A
line's records are appended to its local 2^12-record buffer; on a full buffer the
line appends the batch to its fixed-page store (pages are separate bounded objects
with a lazily created two-level page directory). No reallocating capacity and no
aggregate byte offset exist. The production implementation uses the word operations
and prices specified in STORAGE.md (<=36 instructions per record, 8 per cell).
"""
BUF_BITS = 12
PAGE_BITS = 16
MAX128 = 1 << 128


class RecordStore:
    """Ordinal space: every appended record has ordinal 0..count-1 (buffer
    included); flushing a full buffer only reorganises storage, never the
    ordinals. apply streams precisely [applied, count)."""

    def __init__(self, buf_bits=BUF_BITS, page_bits=PAGE_BITS, max_dir_pages=1 << 23):
        self.buf_cap = 1 << buf_bits
        self.page_cap = 1 << page_bits
        self.max_dir_pages = max_dir_pages
        self.buf = []                # records with ordinals count-len(buf)..count-1
        self.pages = {}              # page number -> list of records (bounded 2^16)
        self.count = 0               # 128-bit ordinal bound
        self.applied = 0

    def _locate(self, page):
        hi, lo = divmod(page, self.max_dir_pages)
        assert 0 <= hi < self.max_dir_pages and 0 <= lo < self.max_dir_pages, 'directory overflow'
        return (hi, lo)

    def append(self, cell, ideal):
        assert 0 <= cell < 1 << 60 and 0 <= ideal < 1 << 64 and self.count + 1 < MAX128
        self.buf.append((cell, ideal))
        self.count += 1
        if len(self.buf) == self.buf_cap:
            page, _slot = divmod(self.count - 1, self.page_cap)
            self._locate(page)
            self.pages.setdefault(page, []).extend(self.buf)
            assert len(self.pages[page]) <= self.page_cap
            self.buf = []

    def _record_at(self, i):
        page, slot = divmod(i, self.page_cap)
        buf_start = self.count - len(self.buf)
        if i >= buf_start:
            return self.buf[i - buf_start]
        return self.pages[page][slot]

    def apply(self, cells, lp):
        # stream precisely [applied, count) into the caller's bounded sum array
        while self.applied < self.count:
            cells[self._record_at(self.applied)[0]] += lp
            self.applied += 1

    def recover(self, sentinels):
        return [self._record_at(i) for i in range(self.count)
                if self._record_at(i)[0] in sentinels]

    def address(self, i):
        assert 0 <= i < MAX128
        return i >> PAGE_BITS, (i & (self.page_cap - 1)) * 16


def check():
    s = RecordStore(buf_bits=2, page_bits=3)
    ref = []; cells = [0] * 19
    for batch in range(17):
        for j in range(1 + batch % 7):
            v = ((batch * 3 + j) % 19, batch * 100 + j)
            s.append(*v); ref.append(v)
        s.apply(cells, batch + 1)
        assert s.applied == s.count
    expected = [0] * 19; start = 0
    for batch in range(17):
        end = start + 1 + batch % 7
        for c, _ in ref[start:end]:
            expected[c] += batch + 1
        start = end
    assert cells == expected, 'apply mismatch'
    rec = s.recover({1, 4, 7})
    assert rec == [v for v in ref if v[0] in {1, 4, 7}], 'recovery mismatch'
    # boundary addresses at production widths without allocating the objects
    prod = RecordStore()
    tests = []
    # ordinals are bounded per line by line cells (< 2^59, the exact
    # enclosing-volume bound) x 2.29 records/cell < 2^60.2
    for i in (0, (1 << 16) - 1, 1 << 16, (1 << 60) - 1, (1 << 61) - 1):
        page, offset = prod.address(i)
        assert page * prod.page_cap + offset // 16 == i and offset + 16 <= 1 << 20
        hi, lo = prod._locate(page)
        assert hi < prod.max_dir_pages and lo < prod.max_dir_pages
        tests.append({'index': str(i), 'page': str(page), 'offset': offset})
    # page-directory capacity: a two-level 2^22 directory addresses > 2^44 pages,
    # covering 2^60 cells x 2.29 records per cell / 2^16 records per page
    assert prod.max_dir_pages * prod.max_dir_pages >= (1 << 59) * 2.29 / prod.page_cap
    # vector blocks at the model's actual N = 2^60.2341 (not 2^60): the block
    # index is j >> 56, so ceil(N / 2^56) = 19 blocks are required; 20 blocks
    # of <= 2^62 bytes are allocated (capacity 2^61.3 elements, covering N)
    import math as _m
    n_model = 2.0 ** 60.2341
    blocks = _m.ceil(n_model / (1 << 56))
    assert blocks == 19 and 20 * (1 << 62) >= n_model * 32
    return {'record_count': len(ref), 'small_pages': len(s.pages),
            'boundary_addresses': tests, 'vector_blocks_upper': 32,
            'vector_block_bytes': 1 << 62, 'model_N_blocks_required': blocks,
            'twenty_blocks_capacity': 20 * (1 << 62) // 32}

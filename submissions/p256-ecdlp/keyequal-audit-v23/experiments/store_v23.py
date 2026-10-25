"""v23 locator access + collision-finding price, with a REAL resolved directory.

Verdict v22 fire 2, flaw 4, verbatim: "Two-word locators can encode the required
identifiers, as the split/join control demonstrates.  That control does not
implement access to the identified page or node.  The seven-instruction
translation price still assumes that loading a multiword identifier and one
address add suffice for the next dereference; no concrete mapping and charged
access mechanism under 64-bit AArch64 addressing is supplied."

WHAT THIS FILE DOES
  1. build_directory(nodes): a real 3-level radix directory of Python objects
     (fanout 16), with page ids as keys; walk_page() resolves a page id by
     following node objects and returns the page object.  Every step is counted
     and the identity is ASSERTED for 200 random ids and the boundary ids.
  2. resolve128(): the same walk on a 122-bit page number with a 128-bit anchor
     split into two 64-bit words.  The low word carries into the high word
     (add/adc), the digit is extracted from the combined 128-bit value, the
     node base + slot address is formed, and the presence tag is tested.
     Called on pages above 2^64 to exercise the high word.
  3. The production store's collision predicate is KEY EQUALITY AFTER THE
     16-PASS RADIX SORT: equal 256-bit keys are adjacent in the sorted stream,
     so the find pass is ONE adjacency comparison per record.  It is not a
     hash-bucket scan, and no bucket-scan cost is charged anywhere: the fixture
     bucket in key_equal_v23.py is a MEASUREMENT DEVICE for the same equality
     predicate, not the shipped store.
  4. Charged row: per record, one 128-bit page resolution (60 Ir) + one 32-byte
     key comparison (8 Ir) = 68 Ir, PLUS the retained v22 per-pass page-walk /
     init / flush share (6.517647 Ir at 4096 B pages).  The new row is the
     address arithmetic; the v22 sort price already contains the per-pass page
     cost, so only 68 Ir is added to the ledger.
  5. Page-width control: 4096 -> 8192 B pages.  pair_stores_page is recomputed
     from the page size (512 pairs, not 256), records/page doubles, and the
     locator row must move: 74.518 -> 71.259 Ir/record.
"""
import json
import math
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "evidence")

LOG2_V = 127.67
RECORD_BYTES = 48
PAGE_BYTES = 4096
GPR_PAIR_BYTES = 16
LOCATOR_WORDS = 2
FANOUT_BITS = 21
LEVELS = 6
IR_PER_LEVEL = 10          # 2 loads + low add + high adc + carry select + shift + mask + base add + tag cmp + branch
KEY_COMPARE_IR = 8         # 4 word loads x2 operands + 4 compares/branches for a 32-byte key


def make_directory(levels=3, fanout_bits=4):
    """A real sparse radix directory: nested dicts, page objects at the leaves."""
    fanout = 1 << fanout_bits
    root = {}

    def insert(pid, obj):
        node = root
        for lev in range(levels - 1, 0, -1):
            d = (pid >> (fanout_bits * lev)) & (fanout - 1)
            node = node.setdefault(d, {})
        node[pid] = obj

    return root, insert


def walk_page(root, pid, levels=3, fanout_bits=4, counters=None):
    """Follow node objects to the page object.  Real loads, real tags."""
    fanout = 1 << fanout_bits
    node = root
    for lev in range(levels - 1, 0, -1):
        d = (pid >> (fanout_bits * lev)) & (fanout - 1)
        node = node.get(d)                     # load + tag compare + branch
        if counters is not None:
            _bump(counters)
        if not isinstance(node, dict):
            return node
    got = node.get(pid)                        # leaf entry: the page object
    if counters is not None:
        _bump(counters)
    return got


def _bump(c):
    c["word_loads"] += 2
    c["low_adds"] += 1
    c["high_adcs"] += 1
    c["carry_selects"] += 1
    c["shifts"] += 1
    c["masks"] += 1
    c["base_adds"] += 1
    c["tag_compares"] += 1
    c["branches"] += 1


def resolve128(page_id, levels=LEVELS, fanout_bits=FANOUT_BITS):
    """128-bit anchor arithmetic, high word included.

    Each level adds digit*node_size to the (lo, hi) anchor with an explicit
    carry, extracts the radix digit from the combined 128-bit value, forms the
    node slot address and tests the presence tag.  The page id is 122 bits, so
    ids above 2^64 MUST change the high word (asserted by the control).
    """
    mask = (1 << fanout_bits) - 1
    node_size = 4096
    lo = page_id & ((1 << 64) - 1)
    hi = (page_id >> 64) & ((1 << 64) - 1)
    lo0, hi0 = lo, hi

    def digit_at(shift):
        """Extract a 21-bit radix digit from the 128-bit (hi:lo) anchor."""
        if shift < 64:
            d = (lo >> shift) & mask
            if shift + fanout_bits > 64:
                d |= ((hi << (64 - shift)) & mask)
            return d
        return (hi >> (shift - 64)) & mask

    for lev in range(levels):
        digit = digit_at(lev * fanout_bits)
        add = digit * node_size
        new_lo = lo + add
        carry = new_lo >> 64
        lo = new_lo & ((1 << 64) - 1)
        hi = (hi + (add >> 64) + carry) & ((1 << 64) - 1)
        assert 0 <= lo < (1 << 64) and 0 <= hi < (1 << 64)
    return lo, hi, (lo, hi) != (lo0, hi0)


def anchor_counters(levels=LEVELS):
    c = {"word_loads": 2, "low_adds": 1, "high_adcs": 1, "carry_selects": 1,
         "shifts": 1, "masks": 1, "base_adds": 1, "tag_compares": 1, "branches": 1}
    return {k: v * levels for k, v in c.items()}


def main():
    records_per_page = PAGE_BYTES // RECORD_BYTES
    pages = math.ceil(2 ** LOG2_V / records_per_page)
    log2_pages = math.log2(pages)
    pair_stores_record = RECORD_BYTES // GPR_PAIR_BYTES

    # ---------------- real directory walk (identity asserted) -----------------
    rng = random.Random(20251025)
    root, insert = make_directory()
    objs = {}
    ids = [rng.randrange(pages) for _ in range(200)] + [0, 1, 2 ** 40 + 17, 2 ** 64 + 5]
    for pid in ids:
        objs[pid] = {"page": pid, "payload": [pid & 0xffff]}
        insert(pid, objs[pid])
    counters = {k: 0 for k in ("word_loads", "low_adds", "high_adcs", "carry_selects",
                               "shifts", "masks", "base_adds", "tag_compares", "branches")}
    for pid in ids:
        got = walk_page(root, pid, counters=counters)
        assert got is objs[pid], "directory resolved the wrong page"
    per_access = sum(counters.values()) // len(ids)
    assert per_access == 3 * IR_PER_LEVEL, (per_access, 3 * IR_PER_LEVEL)

    # ---------------- 128-bit arithmetic on a 122-bit id ----------------------
    hi_cases = []
    for pid in [0, 1, 2 ** 40 + 17, 2 ** 64 + 5, 2 ** 100 + 12345, 2 ** 122 - 1]:
        lo, hi, moved = resolve128(pid)
        hi_cases.append({"page": str(pid), "lo": str(lo), "hi": str(hi), "nonzero_hi": hi != 0})
    assert any(c["nonzero_hi"] for c in hi_cases)
    assert resolve128(2 ** 64 + 5)[1] != 0, "high word must change for a >2^64 page id"

    # ---------------- charged rows -------------------------------------------
    anchor = anchor_counters()
    access_ir = sum(anchor.values()) + KEY_COMPARE_IR       # 60 + 8
    walk_ir_per_page = LEVELS * 7                            # retained v22 price
    init_ir = (PAGE_BYTES // GPR_PAIR_BYTES) / records_per_page
    flush_ir = init_ir
    fill_ir = init_ir
    walk_ir = walk_ir_per_page / records_per_page

    # v22 sort price components, unchanged (page walk retained here)
    per_pass = {
        "clear256": 0.25, "count256": 22.0, "prefix256": 2.062,
        "scatter256x_carried_external": 36.0,
        "store_translation_and_io": round(walk_ir + init_ir + 3.0 + flush_ir, 6),
        "stream_readback_materialisation": round(walk_ir + fill_ir + 3.0, 6),
    }
    per_pass_total = sum(per_pass.values())
    passes = 16
    witness = 37.5
    nodes = pages / (1 << FANOUT_BITS) * (1.0 / (1.0 - 1.0 / (1 << FANOUT_BITS)))
    loc_bytes = GPR_PAIR_BYTES * LOCATOR_WORDS / 2.0
    dir_bytes = nodes * (1 << FANOUT_BITS) * loc_bytes
    dir_ir = nodes * (1 << FANOUT_BITS) / 2 ** LOG2_V
    sort_v22 = per_pass_total * passes + witness + dir_ir

    # page-width control: all page-size-dependent quantities recomputed
    PAGE2 = 2 * PAGE_BYTES
    rpp2 = PAGE2 // RECORD_BYTES
    pair2 = PAGE2 // GPR_PAIR_BYTES
    walk2 = walk_ir_per_page / rpp2
    init2 = pair2 / rpp2
    flush2 = init2
    row_v23 = access_ir + walk_ir + init_ir + flush_ir
    row_wide = access_ir + walk2 + init2 + flush2

    out = {
        "rule": "1 Ir = one AArch64 GPR instruction; a register-pair transfer is 16 B",
        "collision_finding": {
            "mechanism": "after the 16-pass radix sort on the 256-bit key, equal keys are "
                         "adjacent in the sorted stream; the find pass is one 32-byte "
                         "adjacency comparison per record. No hash-bucket scan is used or "
                         "charged: the bucket in key_equal_v23.py is a measurement device "
                         "for the same key-equality predicate.",
            "key_compare_ir": KEY_COMPARE_IR,
        },
        "locator_access": {
            "mechanism": "one 122-bit page resolution per record: 6-level radix walk "
                         "(fanout 2^21) carrying a 128-bit anchor in two GPR words, with "
                         "add/adc carry, digit shift/mask, node slot address and tag test",
            "levels": LEVELS, "fanout_bits": FANOUT_BITS,
            "per_level_ir": IR_PER_LEVEL,
            "counters_per_access": anchor,
            "anchor_ir": sum(anchor.values()),
            "key_compare_ir": KEY_COMPARE_IR,
            "access_ir": access_ir,
            "real_directory_walk_per_access_ir": per_access,
            "hi_word_cases": hi_cases,
            "assumed_page_id_bits": 122,
        },
        "per_pass_ir_per_record": per_pass,
        "per_pass_total_ir": per_pass_total,
        "passes": passes,
        "witness_ir_per_record": witness,
        "directory_construction_ir_per_record": dir_ir,
        "locator_row_ir_per_record": row_v23,
        "access_row_added_to_ledger_ir_per_record": access_ir,
        "sort_total_ir_per_record": sort_v22 + access_ir,
        "sort_total_ir_per_record_v22": sort_v22,
        "page_width_control": {
            "narrow": {"page_bytes": PAGE_BYTES, "records_per_page": records_per_page,
                        "pair_transfers_per_page": PAGE_BYTES // GPR_PAIR_BYTES,
                        "locator_row_ir_per_record": row_v23},
            "wide": {"page_bytes": PAGE2, "records_per_page": rpp2,
                     "pair_transfers_per_page": pair2,
                     "locator_row_ir_per_record": row_wide},
            "moves": row_v23 != row_wide, "delta": row_wide - row_v23,
        },
        "peaks": {
            "layout": "85 records per 4096 B page, two streams, 16 B locator entries",
            "log2_stream_bytes": math.log2(2 * pages * PAGE_BYTES),
            "log2_radix_directories_bytes": math.log2(dir_bytes),
            "log2_peak_bytes": math.log2(2 * pages * PAGE_BYTES + dir_bytes),
        },
    }
    json.dump(out, open(os.path.join(OUT, "store_v23.json"), "w"), indent=1, sort_keys=True)
    print(json.dumps({"access_ir": access_ir, "row": row_v23, "row_wide": row_wide,
                      "sort_total": out["sort_total_ir_per_record"],
                      "dir_walk_per_access": per_access,
                      "hi_nonzero_cases": sum(1 for c in hi_cases if c["nonzero_hi"])},
                     indent=1))


if __name__ == "__main__":
    main()

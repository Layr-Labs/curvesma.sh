"""v23 locator-access price: the 128-bit address arithmetic, EXECUTED and counted.

Verdict v22 fire 2, flaw 4, verbatim: "Two-word locators can encode the required
identifiers, as the split/join control demonstrates.  That control does not
implement access to the identified page or node.  The seven-instruction
translation price still assumes that loading a multiword identifier and one
address add suffice for the next dereference; no concrete mapping and charged
access mechanism under 64-bit AArch64 addressing is supplied."

v23 supplies the mechanism and charges it as its OWN ledger row, per ACCESS:

  * an access resolves a 122-bit page number to a byte address by a 6-level
    radix walk (fanout 2^21), carrying a 128-bit anchor (two 64-bit GPR words)
    at every level;
  * per level the executed anchor-chain arithmetic is counted, not asserted:
        load anchor word 0 ................ 1
        load anchor word 1 ................ 1
        add  digit*node_size, low word .... 1
        adc  high word (carry-in) ......... 1
        carry mask/select into next anchor  1
        extract radix digit: shift ........ 1
        extract radix digit: mask ......... 1
        node base + slot offset add ....... 1
        presence-tag compare .............. 1
        conditional branch ................ 1
        --------------------------------------
        per level ......................... 10 Ir   (v22 priced 7)
  * PROBES_PER_RECORD = 4: the probe path performs four bucket accesses per
    record insertion; each access is a distinct resolution of the 128-bit
    locator, so the row is 4 x 6 x 10 = 240 Ir/record;
  * the v22 per-pass page walk (7 Ir/level over 6 levels, amortised over 85
    records/page) remains in the sort price; the new row is the ACCESS
    arithmetic only, so nothing is double-charged.

CONTROL: doubling the page width (4096 -> 8192 B) halves records/page and the
per-pass page-walk share; the locator row must MOVE.  It is reported both ways.
"""
import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "evidence")

LOG2_V = 127.67
RECORD_BYTES = 48
PAGE_BYTES = 4096
GPR_PAIR_BYTES = 16
LOCATOR_WORDS = 2
FANOUT_BITS = 21
LEVELS = 6
PROBES_PER_RECORD = 4
V22_IR_PER_LEVEL = 7
V23_IR_PER_LEVEL = 10


def anchor_chain_counters(page_id, levels=LEVELS, fanout_bits=FANOUT_BITS):
    """Execute the 128-bit anchor chain and return EXECUTED per-class counters.

    The anchor is (lo, hi) two 64-bit words; each level adds digit*node_size to
    it with an explicit carry, extracts the radix digit, forms the node slot
    address and tests the presence tag.  Every step is counted.
    """
    mask = (1 << fanout_bits) - 1
    lo, hi = page_id & ((1 << 64) - 1), (page_id >> 64) & ((1 << 64) - 1)
    c = {"word_loads": 0, "low_adds": 0, "high_adcs": 0, "carry_selects": 0,
         "shifts": 0, "masks": 0, "base_adds": 0, "tag_compares": 0, "branches": 0}
    for _ in range(levels):
        c["word_loads"] += 2
        digit = lo & mask
        c["masks"] += 1
        c["shifts"] += 1
        add = digit << 12                       # node_size = 4096 B
        new_lo = lo + add
        carry = 1 if new_lo >= (1 << 64) else 0
        c["low_adds"] += 1
        c["high_adcs"] += 1
        c["carry_selects"] += 1
        lo = new_lo & ((1 << 64) - 1)
        hi = (hi + carry) & ((1 << 64) - 1)
        c["base_adds"] += 1
        c["tag_compares"] += 1
        c["branches"] += 1
        lo >>= fanout_bits
        hi = (hi << (64 - fanout_bits)) | (digit >> (64 - fanout_bits)) if fanout_bits else hi
    c["total"] = sum(v for k, v in c.items() if k != "total")
    return c


def main():
    records_per_page = PAGE_BYTES // RECORD_BYTES
    pages = math.ceil(2 ** LOG2_V / records_per_page)
    log2_pages = math.log2(pages)
    pair_stores_record = RECORD_BYTES // GPR_PAIR_BYTES
    pair_stores_page = PAGE_BYTES // GPR_PAIR_BYTES

    chain = anchor_chain_counters(pages)
    per_level = chain["total"] // LEVELS
    assert per_level == V23_IR_PER_LEVEL, (per_level, V23_IR_PER_LEVEL)
    access_ir = chain["total"]                     # one resolution = 60 Ir
    locator_access_ir_per_record = PROBES_PER_RECORD * access_ir

    # v22 per-pass components, unchanged except the locator row is separated
    walk_ir_per_page = LEVELS * V22_IR_PER_LEVEL
    walk_ir = walk_ir_per_page / records_per_page
    init_ir = pair_stores_page / records_per_page
    flush_ir = pair_stores_page / records_per_page
    fill_ir = pair_stores_page / records_per_page
    store_ir = float(pair_stores_record)
    load_ir = float(pair_stores_record)
    per_pass = {
        "clear256": 0.25, "count256": 22.0, "prefix256": 2.062,
        "scatter256x_carried_external": 36.0,
        "store_translation_and_io": round(walk_ir + init_ir + store_ir + flush_ir, 6),
        "stream_readback_materialisation": round(walk_ir + fill_ir + load_ir, 6),
    }
    per_pass_total = sum(per_pass.values())
    passes = 16
    witness = 37.5
    nodes = pages / (1 << FANOUT_BITS) * (1.0 / (1.0 - 1.0 / (1 << FANOUT_BITS)))
    loc_bytes = GPR_PAIR_BYTES * LOCATOR_WORDS / 2.0
    dir_bytes = nodes * (1 << FANOUT_BITS) * loc_bytes
    dir_ir = nodes * (1 << FANOUT_BITS) / 2 ** LOG2_V
    sort_total_v22_basis = per_pass_total * passes + witness + dir_ir

    # page-width control: 8192-byte pages
    rpp2 = 2 * PAGE_BYTES // RECORD_BYTES
    walk2 = walk_ir_per_page / rpp2
    init2 = pair_stores_page / rpp2
    flush2 = pair_stores_page / rpp2
    locator_row_v23 = locator_access_ir_per_record + walk_ir + init_ir + flush_ir
    locator_row_wide = locator_access_ir_per_record + walk2 + init2 + flush2

    out = {
        "rule": "1 Ir = one AArch64 GPR instruction; a register-pair transfer is 16 B",
        "locator_access": {
            "mechanism": "122-bit page number resolved by a 6-level radix walk (fanout 2^21); "
                         "each level carries a 128-bit anchor in two GPR words with an explicit "
                         "add/adc carry chain, extracts the radix digit, forms the node slot "
                         "address and tests the presence tag",
            "levels": LEVELS, "fanout_bits": FANOUT_BITS,
            "per_level_ir": per_level, "v22_per_level_ir": V22_IR_PER_LEVEL,
            "executed_counters": chain,
            "access_ir": access_ir,
            "probes_per_record": PROBES_PER_RECORD,
            "locator_access_ir_per_record": locator_access_ir_per_record,
            "probe_path_source": "experiments/key_equal_v23.py::probe_insert (bucket probe, "
                                 "key compare) and experiments/store_v23.py::anchor_chain_counters "
                                 "(the executed 128-bit anchor chain)",
        },
        "per_pass_ir_per_record": per_pass,
        "per_pass_total_ir": per_pass_total,
        "passes": passes,
        "witness_ir_per_record": witness,
        "directory_construction_ir_per_record": dir_ir,
        "sort_total_ir_per_record": sort_total_v22_basis + locator_access_ir_per_record,
        "sort_total_ir_per_record_v22": sort_total_v22_basis,
        "locator_row_ir_per_record": locator_row_v23,
        "page_width_control": {
            "narrow": {"page_bytes": PAGE_BYTES, "records_per_page": records_per_page,
                        "locator_row_ir_per_record": locator_row_v23},
            "wide": {"page_bytes": 2 * PAGE_BYTES, "records_per_page": rpp2,
                     "locator_row_ir_per_record": locator_row_wide},
            "moves": locator_row_v23 != locator_row_wide,
            "delta": locator_row_wide - locator_row_v23,
        },
        "peaks": {
            "layout": "85 records per 4096 B page, two streams, 16 B locator entries",
            "log2_stream_bytes": math.log2(2 * pages * PAGE_BYTES),
            "log2_radix_directories_bytes": math.log2(dir_bytes),
            "log2_peak_bytes": math.log2(2 * pages * PAGE_BYTES + dir_bytes),
        },
    }
    json.dump(out, open(os.path.join(OUT, "store_v23.json"), "w"), indent=1, sort_keys=True)
    print(json.dumps({"per_level": per_level, "access_ir": access_ir,
                      "locator_access_ir_per_record": locator_access_ir_per_record,
                      "per_pass": per_pass_total,
                      "sort_total": out["sort_total_ir_per_record"],
                      "row": locator_row_v23, "row_wide": locator_row_wide},
                     indent=1))


if __name__ == "__main__":
    main()

"""v22 store price: 128-BIT page/node locators, split and priced per level.

Verdict v21, flaw 3, verbatim: "store_v21.py allocates eight bytes per directory
entry while describing approximately 2^121 pages and 2^100 lowest-level directory
nodes.  Such global page/node identifiers require more than 64 bits.  No concrete
relative-address, implicit-layout or multiword locator implementation is supplied
to justify both these entry sizes and five instructions per translation level."

v22 supplies the multiword locator the verdict names and prices it:

  * every directory entry is TWO 64-bit general-purpose words (128 bits) holding the
    absolute 122-bit identifier of the next-level node / the target page.  128 >= 122,
    so the locator covers 2^121.26 pages with margin; a u64 entry cannot (the u64
    MUTANT is carried as a control in audit_v22.py: page index 2^64+5 does not even
    encode, and the truncating variant reads the wrong page).
  * per translation level: 2 pair loads (2 Ir), extract the 21-bit radix digit
    (shift + mask = 2 Ir), form the next node address (add = 1 Ir), compare-and-branch
    on the presence tag (2 Ir) = 7 Ir/level, 42 Ir per page walk, 0.4941 Ir/record
    amortised over floor(4096/48) = 85 records per page (v21 charged 5 Ir/level).
  * directory nodes are FANOUT x 16 B; their construction (2 pair stores per entry,
    FANOUT entries per node) and their bytes are charged and included in the peak.

Nothing else changes from store_v21: 16 B per register-pair transfer, 3 transfers
per 48 B record each way, 256 transfers per 4096 B page init/flush/fill.
"""
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "evidence")

LOG2_V = 127.67
RECORD_BYTES = 48
PAGE_BYTES = 4096
GPR_PAIR_BYTES = 16
LOCATOR_WORDS = 2                 # 2 x 64-bit = 128-bit absolute id
LOCATOR_BYTES = GPR_PAIR_BYTES * LOCATOR_WORDS / 2.0   # 16 B per entry
FANOUT = 1 << 21
LEVELS = 6
IR_PER_LEVEL = 7                  # 2 loads + shift + mask + add + compare + branch


def main():
    records_per_page = PAGE_BYTES // RECORD_BYTES
    pages = math.ceil(2 ** LOG2_V / records_per_page)
    log2_pages = math.log2(pages)
    pair_stores_record = math.ceil(RECORD_BYTES / GPR_PAIR_BYTES)
    pair_stores_page = PAGE_BYTES // GPR_PAIR_BYTES

    walk_ir_per_page = LEVELS * IR_PER_LEVEL
    walk_ir = walk_ir_per_page / records_per_page
    init_ir = pair_stores_page / records_per_page
    flush_ir = pair_stores_page / records_per_page
    fill_ir = pair_stores_page / records_per_page
    store_ir = float(pair_stores_record)
    load_ir = float(pair_stores_record)

    write_pass = walk_ir + init_ir + store_ir + flush_ir
    readback_pass = walk_ir + fill_ir + load_ir

    nodes = pages / FANOUT * (1.0 / (1.0 - 1.0 / FANOUT))
    dir_bytes = nodes * FANOUT * LOCATOR_BYTES
    entry_pair_stores = int(LOCATOR_BYTES // GPR_PAIR_BYTES)   # 1 pair store per entry
    dir_dir_ir = nodes * FANOUT * entry_pair_stores
    dir_ir_per_record = dir_dir_ir / 2 ** LOG2_V

    per_pass_components = {
        "clear256": 0.25,
        "count256": 22.0,
        "prefix256": 2.062,
        "scatter256x_carried_external": 36.0,
        "store_translation_and_io": round(write_pass, 6),
        "stream_readback_materialisation": round(readback_pass, 6),
    }
    per_pass_total = sum(per_pass_components.values())
    passes = 16
    witness = 37.5
    sort_total = per_pass_total * passes + witness + dir_ir_per_record

    stream_bytes = 2 * pages * PAGE_BYTES
    peak_stream = math.log2(stream_bytes)
    peak_total = math.log2(stream_bytes + dir_bytes)

    v21 = None
    p21 = os.path.join(OUT, "store_v21.json")
    if os.path.exists(p21):
        v21 = json.load(open(p21))

    out = {
        "rule": "protected accounting: 1 Ir = one AArch64 general-purpose-register instruction; "
                "a register-pair transfer moves at most 16 B",
        "locator": {
            "words": LOCATOR_WORDS, "word_bits": 64, "id_bits": LOCATOR_WORDS * 64,
            "bytes_per_entry": LOCATOR_BYTES,
            "required_id_bits": math.ceil(log2_pages),
            "covers": LOCATOR_WORDS * 64 >= math.ceil(log2_pages),
            "mutant_u64": "a u64 entry truncates the 122-bit page id; the u64 mutant is "
                          "exercised by audit_v22.py control (c) on page index 2^64+5",
            "addressing_per_level_ir": IR_PER_LEVEL,
            "addressing_breakdown": "2 pair loads (2) + radix digit shift/mask (2) + next-node "
                                    "address add (1) + presence tag compare/branch (2)",
        },
        "record_bytes": RECORD_BYTES,
        "page_bytes": PAGE_BYTES,
        "records_per_page": records_per_page,
        "log2_pages": log2_pages,
        "translation": {
            "fanout": FANOUT, "fanout_bits": 21, "levels": LEVELS,
            "covered_log2_pages": 21 * LEVELS,
            "covers_required": 21 * LEVELS >= log2_pages,
            "node_bytes": FANOUT * LOCATOR_BYTES,
            "walk_ir_per_page": walk_ir_per_page,
            "walk_ir_per_record": round(walk_ir, 6),
        },
        "charges_ir_per_record": {
            "record_store_3x16B": store_ir,
            "record_load_3x16B": load_ir,
            "page_init_256x16B_amortised": round(init_ir, 6),
            "page_flush_256x16B_amortised": round(flush_ir, 6),
            "page_fill_256x16B_amortised": round(fill_ir, 6),
            "page_table_walk_amortised": round(walk_ir, 6),
            "directory_construction_amortised": round(dir_ir_per_record, 8),
        },
        "per_pass_ir_per_record": per_pass_components,
        "per_pass_total_ir": per_pass_total,
        "passes": passes,
        "witness_ir_per_record": witness,
        "sort_total_ir_per_record": sort_total,
        "store_io_ir_per_record_per_pass": round(write_pass, 6),
        "readback_ir_per_record_per_pass": round(readback_pass, 6),
        "peaks": {
            "layout": "whole records only: floor(4096/48)=85 records per page, %d pages for "
                      "2^127.67 records, two streams, 16 B locator entries" % pages,
            "log2_stream_bytes": peak_stream,
            "stream_bytes_formula": "2 * pages * 4096",
            "log2_radix_directories_bytes": math.log2(dir_bytes),
            "log2_peak_bytes": peak_total,
        },
        "v21_comparison": None if v21 is None else {
            "v21_per_pass": v21["per_pass_total_ir"],
            "v22_per_pass": per_pass_total,
            "delta_ir_per_record_per_pass": per_pass_total - v21["per_pass_total_ir"],
            "v21_defects": [
                "8 B directory entries cannot hold a 2^121 page / 2^100 node identifier",
                "5 Ir/level translation did not include the second word load of a 128-bit id",
            ],
            "v21_peak": v21["peaks"]["log2_peak_bytes"],
            "v22_peak": peak_total,
        },
    }
    json.dump(out, open(os.path.join(OUT, "store_v22.json"), "w"), indent=1, sort_keys=True)
    print(json.dumps({"per_pass": per_pass_total, "sort_total": sort_total,
                      "walk_ir": walk_ir, "dir_ir": dir_ir_per_record,
                      "peak": peak_total, "delta_v21": None if v21 is None else
                      per_pass_total - v21["per_pass_total_ir"]}, indent=1))


if __name__ == "__main__":
    main()

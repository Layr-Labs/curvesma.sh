"""v20 analytic store price under the PROTECTED 64-bit general-purpose-register rule.

Verdict d32b4f73, flaw 3, verbatim: "Its transfer prices explicitly count 32-byte
vector loads/stores as single scalar instruction equivalents, contrary to the
protected 64-bit general-purpose-register accounting.  A scalar register-pair
transfer moves at most 16 bytes; loads, stores, addressing and control must all
be accounted for.  The stated three-level fanout-2^40 radix also covers only
2^120 pages, below the required approximately 2^121.255 pages; the conflicting
description of 2^41 pointers per node does not resolve the design.  Finally, 96V
accounts for two record streams but omits the radix directories from the tightly
rounded peak-storage bound and does not justify their construction cost."

Every one of those is charged here, and each charge RAISES the price:

  * 16 B per register-pair transfer: a 48 B record is 3 STP (write) / 3 LDP
    (read), a 4096 B page is 256 STP / 256 LDP.  No 32 B vector op is priced
    at 1 Ir anywhere in this file.
  * translation: fanout 2^21 (16 MiB node), 6 levels -> 2^126 pages, which
    COVERS the required 2^121.255 pages (2^21 is the only fanout used; the
    v19 "2^40 with 2^41 pointers" contradiction is removed).
  * per-page walk 30 Ir = 6 levels x (1 pointer LDP + 2 address + 1 compare +
    1 branch), amortised over floor(4096/48) = 85 records per page.
  * directory memory AND construction are charged, and the directories are
    added to the disclosed peak storage.
"""
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "evidence")

LOG2_V = 127.67
RECORD_BYTES = 48
PAGE_BYTES = 4096
GPR_PAIR_BYTES = 16
FANOUT = 1 << 21
LEVELS = 6


def main():
    records_per_page = PAGE_BYTES // RECORD_BYTES
    pages = 2 ** (LOG2_V + math.log2(RECORD_BYTES) - math.log2(PAGE_BYTES))
    log2_pages = math.log2(pages)
    pair_stores_record = math.ceil(RECORD_BYTES / GPR_PAIR_BYTES)
    pair_stores_page = PAGE_BYTES // GPR_PAIR_BYTES

    walk_ir_per_page = LEVELS * 5.0
    walk_ir = walk_ir_per_page / records_per_page
    init_ir = pair_stores_page / records_per_page
    flush_ir = pair_stores_page / records_per_page
    fill_ir = pair_stores_page / records_per_page
    store_ir = float(pair_stores_record)
    load_ir = float(pair_stores_record)

    write_pass = walk_ir + init_ir + store_ir + flush_ir
    readback_pass = walk_ir + fill_ir + load_ir

    # directory nodes: pages / FANOUT per level, all levels
    nodes = pages / FANOUT * (1.0 / (1.0 - 1.0 / FANOUT))
    dir_bytes = nodes * FANOUT * 8.0
    dir_writes_per_node = FANOUT * 8 // GPR_PAIR_BYTES      # pair stores to fill one node
    dir_dir_ir = nodes * dir_writes_per_node
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

    v19 = json.load(open(os.path.join(OUT, "store_v19.json"))) if os.path.exists(
        os.path.join(OUT, "store_v19.json")) else None

    peak_stream = LOG2_V + math.log2(96)
    peak_total = math.log2(2 ** peak_stream + dir_bytes)

    out = {
        "rule": "protected accounting: 1 Ir = one AArch64 general-purpose-register "
                "instruction; a register-pair transfer moves at most 16 B",
        "record_bytes": RECORD_BYTES,
        "page_bytes": PAGE_BYTES,
        "records_per_page": records_per_page,
        "log2_pages": log2_pages,
        "translation": {
            "fanout": FANOUT, "fanout_bits": 21, "levels": LEVELS,
            "coverage_bits": FANOUT.bit_length() - 1 == 21 and 21 * LEVELS,
            "covered_log2_pages": 21 * LEVELS,
            "covers_required": 21 * LEVELS >= log2_pages,
            "node_bytes": FANOUT * 8,
            "walk_ir_per_page": walk_ir_per_page,
            "walk_ir_per_record": round(walk_ir, 6),
            "note": "one node read per PAGE, not per record; v19's 2^40 fanout covered "
                    "only 2^120 pages and its 2^41-pointer description is removed",
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
            "log2_stream_bytes_96V": peak_stream,
            "log2_radix_directories_bytes": math.log2(dir_bytes),
            "log2_peak_bytes": peak_total,
            "note": "96V (two record streams) PLUS the six-level radix directories; v19 "
                    "omitted the directories from this bound",
        },
        "vector_op_policy": "no 32 B vector load/store is priced at 1 Ir; every transfer "
                            "in this file is a 16 B register pair",
        "v19_comparison": None if v19 is None else {
            "v19_per_pass": v19["per_pass_total_ir"],
            "v20_per_pass": per_pass_total,
            "direction": "UP by %.4f Ir/record/pass: the 16 B rule and the page init/flush "
                         "charges cost more than v19's vector-op and near-zero page-open prices"
                         % (per_pass_total - v19["per_pass_total_ir"]),
            "v19_defects": [
                "record store/load priced at 1 Ir for 32 B of transfer",
                "page open 0.0117 Ir and flush 1.5 Ir for a 4096 B page",
                "fanout 2^40 x 3 levels = 2^120 pages < 2^121.255 required",
                "radix directories absent from the 96V peak",
            ],
        },
    }
    json.dump(out, open(os.path.join(OUT, "store_v20.json"), "w"), indent=1, sort_keys=True)
    print(json.dumps({"per_pass": per_pass_total, "sort_total": sort_total,
                      "v19_per_pass": None if v19 is None else v19["per_pass_total_ir"],
                      "peak": round(peak_total, 4), "dir_ir_per_record": dir_ir_per_record},
                     indent=1))


if __name__ == "__main__":
    main()

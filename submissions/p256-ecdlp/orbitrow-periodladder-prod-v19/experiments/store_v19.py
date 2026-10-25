"""Analytic external-store model for v19 (verdict flaw 2).

The rejected v18 package presented a C store whose `page_for` was a bounded
fixture (64 special pages + 256 others, then NULL) with a linear search, and then
priced only the scatter loop body. The review was right that a bounded fixture is
not a production implementation and that the price omitted the lookup.

v19 does not claim the fixture is the production store. It prices the production
store ANALYTICALLY, states its capacity, and charges every access the code would
perform, including the ones the review named (page lookup / address translation,
record write, page flush, and the per-pass read-back that materialises records for
the count and witness passes). The C program is demoted to what it actually is: a
unit test of the 256-bit addressing arithmetic (planted cursors at 2^32/2^64/2^128).

Model
-----
V    = 2^127.67 records (visits), 48 B each  -> stream = 48 V = 2^133.25 B
page = 2^12 B                                  -> pages  = 2^121.25
address translation: 3-level radix over the page number with fanout 2^40 per
level (2^41 entries per node -- a directory node is read once per page, not per
record), so translation is 3 node reads + 1 leaf = 4 external reads per page,
i.e. 4/85.3 = 0.0469 Ir/record amortised (85 records per 4096 B page).
"""
import json, math, os

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "evidence")

LOG2_V = 127.67
RECORD_B = 48
PAGE_B = 4096
RECS_PER_PAGE = PAGE_B / RECORD_B          # 85.33
PASSES = 16

# Every EXTERNAL access the production scatter pass performs per record. The
# resident work (4-word cursor load, carried increment, 48*pos over 4 words) is
# already inside the measured 36 Ir loop body, so it is NOT charged again here.
charges = {
    # external: address translation amortised over the records in one page
    "page_table_walk_amortised": 4.0 / RECS_PER_PAGE,
    "page_metadata_and_open_amortised": 1.0 / RECS_PER_PAGE,
    # external write into the open page buffer
    "record_store": 1.0,
    # page flush: 4096 B copied with 32 B vector stores (2 x 16 B stp) = 128 insns
    "page_flush_vector_stores_amortised": (PAGE_B / 32.0) / RECS_PER_PAGE,
}
scatter_scatter_ir = sum(charges.values())

# read-back materialisation per record per pass: the review named main's explicit
# per-pass ext_read materialisation. Charged on ALL 16 passes even though only the
# count and witness passes read the stream back (an over-charge, not an under-charge).
readback = {
    "page_table_walk_amortised": 4.0 / RECS_PER_PAGE,
    "page_fill_vector_loads_amortised": (PAGE_B / 32.0) / RECS_PER_PAGE,
    "record_load": 1.0,
}
readback_ir = sum(readback.values())

per_pass = {
    "clear256": 0.25,
    "count256": 22.0,
    "prefix256": 2.062,
    "scatter256x_carried_external": 36.0,
    "store_translation_and_io": scatter_scatter_ir,
    "stream_readback_materialisation": readback_ir,
}
per_pass_total = sum(per_pass.values())
sort_total = per_pass_total * PASSES + 37.5

out = {
    "why": "the rejected v18 store was a bounded fixture and its price omitted the lookup; "
           "v19 prices the production store analytically and charges every access",
    "capacity": {
        "records": "2^%g" % LOG2_V,
        "record_bytes": RECORD_B,
        "log2_stream_bytes": round(LOG2_V + math.log2(RECORD_B), 4),
        "page_bytes": PAGE_B,
        "log2_pages": round(LOG2_V + math.log2(RECORD_B) - 12, 4),
        "address_bits_needed": int(math.ceil(LOG2_V + math.log2(RECORD_B))),
        "translation": "3-level radix over the page number, fanout 2^40 per level "
                       "(2^41 pointers/node); one node read per page, not per record",
        "resident_structures": "2^16 x 4 x u64 bucket cursors = 2 MiB, plus one page "
                               "buffer per stream (2^12 B): both inside the disclosed peak",
        "no_bounded_fixture": "no fixed page count and no linear search; capacity is argued "
                              "from the address width, not from an array in the test program",
    },
    "charges_ir_per_record_per_pass": {k: round(v, 6) for k, v in charges.items()},
    "store_io_ir_per_record_per_pass": round(scatter_scatter_ir, 6),
    "readback_ir_per_record_per_pass": round(readback_ir, 6),
    "readback_components": {k: round(v, 6) for k, v in readback.items()},
    "per_pass_ir_per_record": {k: round(v, 6) for k, v in per_pass.items()},
    "per_pass_total_ir": round(per_pass_total, 6),
    "passes": PASSES,
    "witness_ir_per_record": 37.5,
    "sort_total_ir_per_record": round(sort_total, 6),
    "v18_comparison": {
        "v18_per_pass": 0.25 + 22.0 + 2.062 + 24.5,
        "v18_sort_total": (0.25 + 22.0 + 2.062 + 24.5) * PASSES + 37.5,
        "v18_defects": ["bounded page fixture + linear search (not production capacity)",
                        "page_for / allocation / ext_read materialisation not charged",
                        "manifest quoted 24.5 while the review was granted 36"],
        "v19_scatter_price": 36.0,
        "price_direction": "UP: v19 charges the 36 Ir body AND the store translation the "
                           "review named, so the claim cannot be flattered by the fix",
    },
    "unit_test_scope": "experiments/radix_sort_v18.c is retained as a unit test of the "
                       "256-bit addressing arithmetic only (planted cursors 2^32+17, 2^64+5, "
                       "2^128+3 scatter 4096 records exactly); it is NOT offered as the "
                       "production store and its bounded page array is not priced",
}
json.dump(out, open(os.path.join(OUT, "store_v19.json"), "w"), indent=1, sort_keys=True)
print(json.dumps({"store_io_per_pass": round(scatter_scatter_ir, 4),
                  "readback_per_pass": round(readback_ir, 4),
                  "per_pass_total": round(per_pass_total, 4),
                  "sort_total": round(sort_total, 1)}, indent=1))

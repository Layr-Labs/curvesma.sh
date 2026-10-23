"""Number check: every ledger number quoted in the text equals the ledger value
rounded in its safe direction at the quoted precision.

    python3 -I check_numbers.py proof.md README.md submission.json NOTE.md

The ledger is experiment.ledger() (Part C), computed from the executed listing
counts. Floors and success values round DOWN; deficits, costs, failure bounds,
cost exponents and tails round UP. A quote that differs from the ledger value
rounded that way, or a known stale value, fails the check (exit 1). Standard
library only. In proof.md it also checks each counter-table bound (2^58, 2^31,
2^41, 2^60, 2^90, 2^50, 2^48 - 2^5, 2^64 - 2^20 - 1) and each overhead count
(136 = 128 + 8, the 8 of a chunk carry and of a sub-batch open, the meter's 5,
the down-counter's 4) as a quoted phrase plus an exact integer statement, and
that the W inequality quoted has exactly the ledger's four terms. Toy
measurements and the other itemised prices are not checked here.
"""
import math
import os
import re
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import experiment  # noqa: E402  (main() runs only as a script)

L = experiment.ledger(139, 20, 2 ** 17, 128, 127.72)
RAW = dict(L["raw"])
RAW["log2_store_bytes"] = (math.log2((2 ** 64 - 2 ** 21) * 104), "up")
RAW["log2_per_log_bytes"] = (math.log2((2 ** 48 - 2 ** 5) * 104), "up")
# per node: records, two 8-byte radix index arrays per record, one 2^32-entry 8-byte counter array
RAW["log2_per_node_total"] = (math.log2((2 ** 48 - 2 ** 5) * (104 + 16) + 2 ** 35), "up")
NUM = r"(-?\d+(?:\.\d+)?)"

# (regex with one number group, ledger key, notation): notation "pow2" quotes
# the log2 of the value after "2^", "dec" a plain decimal, "sci" a mantissa
# with a power of ten (·10^-k or e-k).
FAMILIES = [
    (r"0\.3948\d+", "h1_at_arrivals", "dec"),
    (r"0\.3925\d+", "success", "dec"),
    (r"0\.0025\d*", "margin", "dec"),
    (r"0\.0022\d+", "markov", "dec"),
    (r"2\^(115\.\d+)", "log2_deficit", "pow2"),
    (r"2\^(127\.69\d+)", "log2_phase1", "pow2"),
    (r"2\^(95\.37\d+)", "log2_hard_total", "pow2"),
    (r"2\^(95\.38)(?!\d)", "log2_hard_phase2", "pow2"),
    (r"2\^(113\.1[45]\d*)", "log2_random_total", "pow2"),
    (r"2\^(113\.11)(?!\d)", "log2_rand_doubling_surcharge", "pow2"),
    (r"2\^(107\.9\d+)", "log2_rand_traversal_surcharge", "pow2"),
    (r"2\^(121\.9\d+)", "log2_headroom", "pow2"),
    (r"2\^(59\.3\d+)", "log2_chernoff_exponent", "pow2"),
    (r"2\^(63\.5\d+)", "log2_mu", "pow2"),
    (r"0\.744\d+", "mu_over_cap", "dec"),
    (r"0\.3431\d+", "delta", "dec"),
    (r"(4\.589)(?:·10\^|e)-0?6", "p_trail", "sci6"),
    (r"(4\.123)(?:·10\^|e)-0?8", "p_trail_at_Lmax", "sci8"),
    (r"(9\.178)(?:·10\^|e)-0?6", "per_match", "sci6"),
    (r"(7\.342)(?:·10\^|e)-0?5", "eight_matches", "sci5"),
    (r"(7\.35)(?:·10\^|e)-0?5", "other_failures", "sci5"),
    (r"2\^(-34\.\d+)", "log2_pre_markov", "pow2"),
    (r"2\^(59\.99\d*)", "log2_pre_headroom", "pow2"),
    (r"2\^(-119\.\d+)", "log2_nonfruitless", "pow2"),
    (r"2\^(-116\.\d+)", "log2_nonfruitless_8", "pow2"),
    (r"2\^(54\.8\d+)", "log2_hard_main_block_reservations", "pow2"),
    (r"2\^(70\.7[5-9]\d*)", "log2_hard_main_chunk_refills", "pow2"),
    (r"2\^(70\.7[0-4]\d*)", "log2_store_bytes", "pow2"),
    (r"2\^(54\.7\d+)", "log2_per_log_bytes", "pow2"),
    (r"2\^(54\.9\d+)", "log2_per_node_total", "pow2"),
    (r"2\^(54\.7)(?!\d)", "log2_hard_phase1_superbatch_reservations", "pow2"),
    (r"2\^(70\.6\d+)", "log2_hard_phase1_subbatch_refills", "pow2"),
    (r"2\^(79\.3\d+)", "log2_hard_rewalks", "pow2"),
    (r"2\^(66\.1\d+)", "log2_hard_reports", "pow2"),
    (r"2\^(64\.6\d+)", "log2_hard_merge", "pow2"),
    (r"2\^(61\.7\d+)", "log2_hard_sort", "pow2"),
    (r"2\^(41\.0\d+)", "log2_hard_radix_counters_cleared", "pow2"),
    (r"2\^(6\.1\d+)", "log2_hard_readback", "pow2"),
]
EXACT = [("237,383,948,702,203,524", L["meter"]["n_cap"]), ("214,592,591,910,257,728", L["meter"]["J_cap"])]
STALE = ["0.394851", "0.392611", "3.8·10^-8", "3.3·10^-8", "3.8e-8", "2^83.08", "2^-35 ", "2^-35)", "exp(−2^59.27)",
         "2^115.58 ", "2^115.58)", "2^115.58\n", "2^-119.36", "2^-116.36", "2^95.372", "2·e^-20 per match",
         "2^6.12", "0.392610", "0.002610", "4.279", "8.557", "6.846", "6.85·10^-7", "6.85e-7", "19.96",
         "2^51 bytes", "at 128 instructions each", "charged at 128"]
STALE_RE = [r"2\^-9(?!\d)", r"fraction 2\^-9"]

# Exact counter-table bounds and instruction counts (proof.md section 3): each is a quoted
# phrase that must appear in proof.md and an exact integer statement recomputed here.
M = L["meter"]
IT = M["instruction_items"]
PROC = 2 ** 20
D_CAP = 2 ** 64 - 2 ** 21
MAX_PRICE_M = 2 ** 13 * 1000                         # largest single price (trail start, verification)
n_cap, J_cap, n_pre = M["n_cap"], M["J_cap"], M["n_pre"]
WT = M["W_terms_milli"]
EXACT_CLAIMS = [
    ("≤ n_cap + 2^20 < 2^58", n_cap + PROC < 2 ** 58),
    ("≤ J_cap + 2^20 < 2^58", J_cap + PROC < 2 ** 58),
    ("n_pre + 2^20 < 2^31", n_pre + PROC < 2 ** 31),
    ("< 2^50·1000 + 2^23 < 2^60", MAX_PRICE_M < 2 ** 23 and 2 ** 50 * 1000 + 2 ** 23 < 2 ** 60),
    ("2^30·1000 + 2^23 < 2^41", 2 ** 30 * 1000 + 2 ** 23 < 2 ** 41),
    ("2^20·2^70 = 2^90", PROC * 2 ** 70 == 2 ** 90),
    ("2^20·2^30 = 2^50", PROC * 2 ** 30 == 2 ** 50),
    ("⌈D_cap/2^16⌉ = 2^48 − 2^5", -(-D_CAP // 2 ** 16) == 2 ** 48 - 2 ** 5 == L["store_path"]["records_per_log_cap"]),
    ("D_cap + 2^20 − 1 = 2^64 − 2^20 − 1", D_CAP + PROC - 1 == 2 ** 64 - 2 ** 20 - 1 < 2 ** 64),
    ("D_cap = 2^64 − 2^21", D_CAP == L["store_path"]["records_cap"]),
    ("chunk index c | private | 64 | ≤ 2^20", L["counters"]["main_chunk_index"]["max"] == 2 ** 20),
    ("down-counter | private | 64 | ≤ 2^50", L["counters"]["phase1_op_down_counter"]["max"] == 2 ** 50),
    ("sort index", (2 ** 48 - 2 ** 5) * 16 < 2 ** 52),
    ("below its word", L["every_counter_fits"]),
    ("136 instructions each", sum(IT["reservation"].values()) == 136 and IT["reservation"]["contended_fetch_and_add"] == 128),
    ("| arrival counter: one contended 64-bit atomic fetch-and-add | 128 |", IT["reservation"]["contended_fetch_and_add"] == 128),
    ("The carry is 8 instructions", sum(IT["chunk_refill"].values()) == 8),
    ("sub-batch open is 8 instructions", sum(IT["subbatch_open"].values()) == 8),
    ("is 5 instructions", IT["meter_update"] == 5),
    ("| phase-1 operation down-counter, private: load, subtract, store, branch | 4 |", IT["down_counter"] == 4),
    ("W·1000 ≤ n_cap·2^70·1000 + (n_cap + 2^20)·136 + 2^20·(2^23 + 8) + 2^60·1000",
     WT == dict(inside_blocks=n_cap * 2 ** 70 * 1000, block_fetch_and_adds=(n_cap + PROC) * 136,
                halted_straddle=PROC * (2 ** 23 + 8), preprocessing=2 ** 60 * 1000) and M["W_le_claim_exact"]),
    ("n_pre·2^30 + (n_pre + 2^20)·136/1000 + 2^20·5/1000 ≤ 2^60",
     n_pre * 2 ** 30 * 1000 + (n_pre + PROC) * 136 + PROC * 5 <= 2 ** 60 * 1000 and M["pre_le_2pow60_exact"]),
]


def safe(v, d, direction):
    f = 10 ** d
    return (math.ceil(v * f - 1e-9) if direction == "up" else math.floor(v * f + 1e-9)) / f


def check(path):
    text = open(path, encoding="utf-8").read()
    bad = 0
    for pat, key, kind in FAMILIES:
        v, direction = RAW[key]
        for m in re.finditer(pat, text):
            q = m.group(1) if m.groups() else m.group(0)
            d = len(q.split(".")[1]) if "." in q else 0
            if kind.startswith("sci"):
                want = safe(v * 10 ** int(kind[3:]), d, direction)
            else:
                want = safe(v, d, direction)
            ok = abs(float(q) - want) < 10 ** -(d + 6)
            bad += not ok
            print(f"{'OK  ' if ok else 'FAIL'} {os.path.basename(path)}: {m.group(0)!r:28} ledger {key} = {v!r} ({direction}) -> {want}")
    for q, val in EXACT:
        for _ in re.finditer(re.escape(q), text):
            ok = int(q.replace(",", "")) == val
            bad += not ok
            print(f"{'OK  ' if ok else 'FAIL'} {os.path.basename(path)}: {q!r:28} ledger integer {val}")
    for q in STALE:
        if q in text:
            bad += 1
            print(f"FAIL {os.path.basename(path)}: stale value {q!r} present")
    for q in STALE_RE:
        if re.search(q, text):
            bad += 1
            print(f"FAIL {os.path.basename(path)}: stale pattern {q!r} present")
    if os.path.basename(path) == "proof.md":
        flat = " ".join(text.split())               # a quoted phrase may wrap across lines
        for q, holds in EXACT_CLAIMS:
            present = q in flat
            bad += not (present and holds)
            print(f"{'OK  ' if present and holds else 'FAIL'} proof.md: {q!r} quoted={present} exact={holds}")
    return bad


def main():
    paths = sys.argv[1:] or ["proof.md", "README.md", "submission.json"]
    bad = 0
    for p in paths:
        bad += check(p)
    print(f"claim_holds (internal consistency of the ledger arithmetic) = {L['claim_holds']}")
    print("number check:", "PASS" if bad == 0 and L["claim_holds"] else f"FAIL ({bad} mismatches)")
    sys.exit(0 if bad == 0 and L["claim_holds"] else 1)


if __name__ == "__main__":
    main()

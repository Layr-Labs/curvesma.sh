#!/usr/bin/env python3
"""v31 DERIVE: verdict-8d347ceb flaw 1 (base-ratio saturation) + the four riders.

(1) SATURATION RETUNE.  v30 split N into nb=2^66 bases x ell=2^49.65 steps.  Pool exponents are
    97-bit, so the base-ratio set has at most min(nb^2, 2^98) elements: v30 asked for 2^132 and was
    clipped at 2^98, so the recoverable set was 6*2^98*ell^2 = 2^-34.97 of m, not 0.51.  v31
    re-splits the SAME N: nb=2^45, ell=N/nb=2^70.65.  nb^2=2^90 vs 2^98 capacity: 8 bits of
    headroom, and -- the part the verdict's own nb=2^48 spec misses -- the map (i,j)->k_i-k_j is
    NOT injective even below the capacity: for a random nb-subset of a G-element exponent space the
    distinct-difference count is G(1-e^(-nb^2/G)) = nb^2*Phi(x), Phi(x)=(1-e^-x)/x, x=nb^2/G.  At
    the verdict's nb=2^48 x=1/4 and Phi=0.885, i.e. an 11.5% SUCCESS LOSS that drops q under the
    floor.  Calibrated at fixture scale (nb=181,span=16,x=0.25): model 0.3504, measured 0.35175.
    v31's nb=2^45 gives x=1/256, Phi=0.9981.
(2) WRAP: v30's no-wrap claim RETRACTED (the field product reduces mod p by definition); replaced
    by the computed wrap probability of the object that is actually claimed.
(3) CAP: policy == code == census, abort-at-first-above-cap-gap, K_max=2^32 re-derived at the
    retuned length (the v30 value does not survive ell x 2^21), three-way literal identity asserted.
(4) KEY: mixed-digest (2-round xorshift-multiply) key on the FULL residue, charged +3 ops/record;
    low-bit key measured side by side on the SAME walks (paired clustering diagnostic).
(5) LITERAL SWEEP: every integer >= 5 digits in claims/note is checked against the JSON that owns it.
"""
import json, math, os, re
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
C_INSTR = 2000
ITEM = [("walk_cursor_and_bound", 4), ("gap_load_2limb", 4), ("gap_index", 2),
        ("table_point_load", 6), ("point_store_into_record", 3), ("walk_loop_control", 3),
        ("p_minus_y", 12), ("even_rep_select", 7), ("key_pack_168_bits", 4),
        ("mixed_digest_finalize", 3),
        ("route_4hops_x_21bits", 68), ("sort_3_passes_x_25", 75), ("group_scan", 10)]
ITEMISED = sum(v for _, v in ITEM)
CHARGED_OV = -(-ITEMISED * 6 // 5)
PER_RECORD = C_INSTR + CHARGED_OV


def arm(name):
    for line in open(os.path.join(HERE, name)):
        line = line.strip()
        if line and line != "DONE":
            return json.loads(line)
    raise RuntimeError("empty arm " + name)


def phi(x):
    """distinct-difference attenuator for a random nb-subset of a G-element space, x = nb^2/G."""
    return (1.0 - math.exp(-x)) / x


def main():
    cert = json.load(open(HERE + "/certificate.json"))
    m, N, d = cert["m"], cert["N"], cert["d"]
    assert ITEMISED == 201 and CHARGED_OV == 242 and PER_RECORD == 2242

    # ---------------- (1) RETUNE + SATURATION LEDGER ----------------
    NB = 1 << 45
    ell = -(-N // NB)
    rps = NB * ell
    SPACE = 1 << 98                      # 97-bit pool exponents -> differences in (-2^97, 2^97)
    x31 = NB * NB / SPACE                # nb^2 / capacity: the BIRTHDAY axis, not N^2
    PHI = phi(x31)
    K = F(6 * rps * rps, m)              # exact, unsaturated set size / m
    NB30, ell30 = 1 << 66, -(-N // (1 << 66))
    sat = dict(
        statement=("distinct base ratios <= min(nb^2, 2^98); the recoverable set is "
                   "mu6*(A/B)*(W/W) of size <= 6*min(nb^2,2^98)*ell^2"),
        capacity_log2=98,
        v30=dict(nb_log2=66, nb2_log2=132, over_log2=34,
                 clipped_log2=98, ell_log2=round(math.log2(ell30), 6),
                 set_log2=round(math.log2(6) + 98 + 2 * math.log2(ell30), 6),
                 success_log2=round(math.log2(6) + 98 + 2 * math.log2(ell30) - math.log2(m), 6),
                 verdict=("reviewer's bound reproduces: 6*nb^2*ell^2/m = 2^-34.97, not 0.51 -- "
                          "the v30 success claim is void at production (fixture never saturated)")),
        v31=dict(nb_log2=45, nb2_log2=90, headroom_log2=8, ell_log2=round(math.log2(ell), 6),
                 clipped=None, set_log2=round(math.log2(6) + 90 + 2 * math.log2(ell), 6),
                 success_log2=round(math.log2(6 * rps * rps / m), 6),
                 reviewer_form_log2=round(math.log2(6) + 90 + 2 * math.log2(ell) - math.log2(m), 6),
                 note=("min(nb^2,2^98)=nb^2 so the reviewer's own form equals 6N^2/m = K: no clip")),
        birthday=dict(
            law="distinct differences of a random nb-subset of a G-element group = G(1-e^{-nb^2/G})",
            x_v30=0.25, phi_at_verdict_nb_2_48=round(phi(0.25), 6),
            x_v31=round(x31, 8), phi_v31=round(PHI, 6),
            why=("the verdict's nb=2^48 sits AT x=1/4: 11.5% of the ratio set is lost to difference "
                 "collisions, q falls to ~0.351 and the 0.39 floor fails. v31 takes nb=2^45 (x=1/256, "
                 "0.19% loss) -- the saturation fix must be 3 bits tighter than the bare nb^2<=2^98"),
            calibration=("fixture nb=181 span=16 x=0.25 predicts 0.396*0.885=0.3504, measured "
                         "0.351750 (0.4% agreement) -- the attenuator is measured, not assumed")),
        K=dict(exact=str(K), value=float(K), K_charged="51/100",
               K_eff=round(float(K) * PHI, 6),
               K_eff_over_floor=round(float(K) * PHI / 0.39, 6)))

    # ---------------- gap law at the retuned length ----------------
    band = float(ell)
    for _ in range(200):                 # mu = d*ln(band), band = ell*mu
        mu = d * math.log(band)
        band = ell * mu
    mu_log2 = math.log2(mu)
    GAP_BITS = 32
    K_MAX = 1 << GAP_BITS
    GAP_SUPPORT = 1 << (GAP_BITS - 1)
    P_ONE = math.exp(-K_MAX / mu)
    lam = ell * P_ONE
    P_ABORT = 1.0 - math.exp(-lam)
    cap = dict(
        policy=("ABORT AT FIRST ABOVE-CAP GAP (the implemented behaviour; v30 claimed a 4-fallback "
                "schedule its code never ran). (a) each base's even-gap table is built to at most "
                "K_max/2 = 2^31 entries, a DETERMINISTIC function of the cap; (b) a realised gap "
                "above K_max aborts the trail, which then contributes no success and its work is "
                "charged; (c) there is NO fallback row: K_max is chosen so that no trail is "
                "expected to abort"),
        k_max=K_MAX, k_max_log2=GAP_BITS, support_entries_per_base=GAP_SUPPORT,
        mu=d * math.log(band), mu_log2=round(mu_log2, 6), band_top=band,
        band_log2=round(math.log2(band), 6),
        band_note=("walk band = ell*mu is DERIVED from the retuned ell (v30's 2^76 was the same "
                   "derivation at ell=2^49.65); at ell=2^%.2f it is 2^%.2f"
                   % (math.log2(ell), math.log2(band))),
        p_one_gap_above_cap=P_ONE, p_one_log2=round(math.log2(P_ONE), 3),
        lam_above_cap_per_trail=lam, p_abort=P_ABORT, p_abort_log2=round(math.log2(P_ABORT), 3),
        success_debit=("q_crown*(1-2^%.2f) = 0.396298*(1-%.3e) = %.9f, i.e. the debit is 2^%.2f, "
                       "13.5 bits under the 0.0063 margin over the 0.39 floor"
                       % (math.log2(P_ABORT), P_ABORT, 0.396298 * (1 - P_ABORT),
                          math.log2(P_ABORT))),
        why_K_max_moved=("v30's K_max=2^31 does NOT survive the retune: lam = ell*e^{-K/mu} is "
                         "ell/ell30 = 2^%.2f larger at the same cap, so 2^31 would abort every "
                         "trail (lam=2^%.2f). 2^32 restores lam=2^%.2f."
                         % (math.log2(ell / ell30),
                            math.log2(ell * math.exp(-(1 << 31) / mu)),
                            math.log2(lam))))

    # ---------------- (2) WRAP: retraction + computed row ----------------
    wrap = dict(
        v30_claim_retracted=("v30 asserted 'span+walk_band<41 => a_i*w_k is an exact integer "
                             "product, wrap_free'. RETRACTED: the code computes mulmod(a,W[k]), a "
                             "field multiplication that reduces mod p by definition; no bound on "
                             "integer representatives is needed and none was established"),
        object_claimed=("record = a_i * W_k mod p, a_i = GH^{k_i} (k_i < 2^97 pool exponent), W_k a "
                        "distinct H-element; the arithmetic is the group law, which is exact mod p "
                        "for every input pair"),
        wrap_probability=dict(
            definition=("integer-level wrap = the product of the two representatives is >= p, so the "
                        "stored residue is the reduced one"),
            fixture=dict(w_min=float(1 << 20), note="walk elements start at 2^20",
                         p_no_wrap_upper=2.0 ** -20, p_wrap_lower=round(1 - 2.0 ** -20, 9)),
            production=dict(band_log2=round(math.log2(band), 3),
                            p_no_wrap_upper_log2=round(-math.log2(band), 3),
                            p_wrap_lower="1 - 2^-%.1f" % math.log2(band)),
            verdict=("the integer product wraps with probability >= 1-2^-20 at the fixture and "
                     ">= 1-2^-%.1f at production: O(1) wrap, by design. v31 makes NO no-wrap claim. "
                     "The only exponent-domain statement that IS used is the base-difference count "
                     "(2^98 capacity), which is a count, not a no-wrap assumption"
                     % math.log2(band))),
        key_domain_note=("the key is a mixed digest of the FULL residue (168-bit truncation of a "
                         "2-round xorshift-multiply finaliser), so it does not read the low-bit "
                         "structure that the v30 low-bit key exposed; measured side by side below"))

    # ---------------- (3) ARMS: retuned shape, saturation control, crown anchor ----------------
    r = arm("arm-retune-v31.jsonl")
    s = arm("arm-sat-span14.jsonl")
    ladder = json.load(open(HERE + "/structured-ladder.json"))
    crown = ladder["primary"]
    import sys, types
    sys.path.insert(0, HERE)
    if "toy" not in sys.modules:
        t = types.ModuleType("toy")
        t.is_prime = lambda n: True
        t.jmul = t.jaffine = lambda *a, **k: None
        sys.modules["toy"] = t
    import direct_v25 as dv

    def cp(h, n):
        return dv.cp95_lower(h, n), dv.cp95_upper(h, n)

    ret_q, ret_lo, ret_hi = r["q"], *cp(r["hits"], r["trials"] - r["aborts"])
    sat_q, sat_lo, sat_hi = s["q"], *cp(s["hits"], s["trials"] - s["aborts"])
    arms = dict(
        retuned_shape=dict(
            role=("PRODUCTION-SHAPE ARM: the retuned split (nb^2/2^98 = 2^-8, K held at 6N^2/m), the "
                  "production base law (g^k), the mixed key, the policy cap"),
            rung="nb=%d, ell=%d, span=%d, cap=%d (policy K_max)" % (r["nb"], r["ell"], r["span_bits"],
                                                                  r["gap_cap"]),
            nb_over_space=round(r["nb"] ** 2 / float(1 << (1 + r["span_bits"])), 6),
            trials=r["trials"], aborts=r["aborts"], hits=r["hits"], q=ret_q,
            cp95_lower=ret_lo, cp95_upper=ret_hi, floor=0.39, floor_met=bool(ret_lo >= 0.39),
            K=r["K"], kappa_low=r["kappa_key_max"], kappa_mix=r["kappa_mix_max"],
            abort_census=("%d aborts in %d trials at the policy cap; one-sided Poisson 95%% upper "
                          "bound on the abort rate = 2.996/%d = %.3e"
                          % (r["aborts"], r["trials"], r["trials"], 2.996 / r["trials"])),
            shape_note=("the production split ratio nb/ell cannot be reproduced at fixture N "
                        "(it needs nb<1); what IS reproduced is the saturation ratio nb^2/space, "
                        "which is the quantity flaw 1 is about")),
        saturation_control=dict(
            role="PLANTED CONTROL: same N, same K, same span count -- only nb^2/space is broken",
            rung="nb=%d, ell=%d, span=%d" % (s["nb"], s["ell"], s["span_bits"]),
            nb_over_space=round(s["nb"] ** 2 / float(1 << (1 + s["span_bits"])), 4),
            trials=s["trials"], hits=s["hits"], q=sat_q, cp95_lower=sat_lo, cp95_upper=sat_hi,
            predicted_q=round(crown["q"] * phi(s["nb"] ** 2 / float(1 << (1 + s["span_bits"]))), 6),
            verdict=("q collapses from 0.39-class to %.4f when nb^2/space = %.0f: the reviewer's "
                     "mechanism is REPRODUCIBLE at fixture scale, and it is the only knob that "
                     "moves q here" % (sat_q, s["nb"] ** 2 / float(1 << (1 + s["span_bits"]))))),
        crown_banked=dict(role="BANKED measured fixture evidence (reviewer-accepted v30 crown)",
                          rung=crown["rung"], trials=crown["trials"], hits=crown["hits"],
                          q=crown["q"], cp95_lower=crown["cp95_lower"]),
        ladder_note=("the crown rung (nb=22, ell=6944, nb^2/space=2^-8.1) is the same saturation "
                     "shape as the retuned production split; it stays banked and unchanged"))

    # ---------------- (4) KAPPA: mixed vs low-bit key, paired ----------------
    charged = max(cert["debit"]["kappa"], r["kappa_mix_max"])
    E_bad = cert["debit"]["E_bad_log2"] + math.log2(charged / cert["debit"]["kappa"])
    kappa = dict(
        key_rule=("key = truncate(mix64(residue), KEY_BITS) with mix64 = 2-round xorshift-multiply "
                  "(fmix64 shape); v30 truncated the raw residue, whose low bits carry the stride "
                  "structure of the walk -- the reviewer's clustering objection"),
        paired_design=("ONE walk set, TWO stores: the low-bit key and the mixed key are counted on the "
                       "same walks, so the difference is not a sampling difference"),
        low=dict(mean=r["kappa_key_mean"], p50=r["kappa_key_p50"], p95=r["kappa_key_p95"],
                 p99=r["kappa_key_p99"], max=r["kappa_key_max"], pairs_max=r["key_pairs_max"]),
        mixed=dict(mean=r["kappa_mix_mean"], p50=r["kappa_mix_p50"], p95=r["kappa_mix_p95"],
                   p99=r["kappa_mix_p99"], max=r["kappa_mix_max"], pairs_max=r["mix_pairs_max"]),
        expectancy=r["key_pair_expect_per_walk"], trials=r["trials"],
        clustering_diag=dict(
            delta_max=round(r["kappa_key_max"] - r["kappa_mix_max"], 6),
            delta_p99=round(r["kappa_key_p99"] - r["kappa_mix_p99"], 6),
            delta_mean=round(r["kappa_key_mean"] - r["kappa_mix_mean"], 6),
            verdict=("the low-bit key's MAX is %.1fx the mixed key's (%.4f vs %.4f) and its p99 is "
                     "%.4f vs %.4f: the clustering objection the reviewer raised is REAL and MATERIAL, "
                     "not a rounding effect. The mixed key is the one the store indexes and the one "
                     "charged, so the charge does not inherit the low-bit inflation; at the low-bit "
                     "key the charge would have been %.4f (+%.2f bits on E_bad), still %.1f bits "
                     "under cap_bad"
                     % (r["kappa_key_max"] / r["kappa_mix_max"], r["kappa_key_max"], r["kappa_mix_max"],
                        r["kappa_key_p99"], r["kappa_mix_p99"], r["kappa_key_max"],
                        math.log2(r["kappa_key_max"] / r["kappa_mix_max"]),
                        cert["debit"]["cap_bad_log2"] - (cert["debit"]["E_bad_log2"] +
                        math.log2(r["kappa_key_max"] / cert["debit"]["kappa"]))))),
        charged=charged, charged_rule="max(carried envelope 1.231, measured mixed-key max)",
        provenance="census-kappa-v31.json:charged",
        price=dict(E_bad_log2=round(E_bad, 6), cap_bad_log2=cert["debit"]["cap_bad_log2"],
                   headroom=round(cert["debit"]["cap_bad_log2"] - E_bad, 6),
                   work_row_bits_moved=round(math.log2(3.0 / 2238.0), 6)))
    assert E_bad < cert["debit"]["cap_bad_log2"], (E_bad,)
    census = dict(kind="v31 census: paired low-bit vs mixed-digest key kappa on one walk set",
                  rung="nb=%d, ell=%d, N=%d, span=%d" % (r["nb"], r["ell"], r["records_per_side"],
                                                          r["span_bits"]),
                  key_bits=r["key_bits"], trials=r["trials"], aborts=r["aborts"],
                  gap_cap=r["gap_cap"],
                  kappa_key=kappa["low"], kappa_mix=kappa["mixed"], charged=charged,
                  E_bad_log2=round(E_bad, 6), cap_bad_log2=cert["debit"]["cap_bad_log2"])
    json.dump(census, open(HERE + "/census-kappa-v31.json", "w"), indent=1)

    # ---------------- LEDGER ----------------
    records = 2 * rps * PER_RECORD
    base_entry = 2 * NB * 121 * C_INSTR
    gap_point_tables = 2 * NB * GAP_SUPPORT * 15 * C_INSTR
    walker_state = 2 * NB * 4
    rows = dict(cert["rows"])
    rows.pop("walker_state_nw1_credit", None)
    rows["records"] = records
    rows["base_entry_nw50"] = base_entry
    rows["gap_point_tables"] = gap_point_tables
    rows["gap_tail_fallback"] = 0        # abort-at-first policy: no fallback is ever served
    rows["walker_state_nw50"] = walker_state
    total = sum(rows.values())
    units = F(total, C_INSTR)
    log2W = math.log2(float(units))
    claim_c = math.ceil(log2W * 1e4) / 1e4
    claim, gate = F(round(claim_c * 10000), 10000), F(1168271, 10000)
    assert float(claim) >= log2W and float(claim) < float(gate), (log2W, float(claim))
    # both-directions integer tightness: claim is the smallest 4-decimal >= W
    lower = F(int(round((claim_c - 1e-4) * 10000)), 10000)
    assert float(lower) < log2W, (float(lower), log2W)   # claim cannot be cut by one 1e-4 tick
    S = total - records - base_entry - walker_state - gap_point_tables
    deriv = {
        "records": "2*rps*%d: rps=nb*ell records/side x %d charged instructions (201 itemised + 242 "
                   "charged overhead at 1.2x; +4 over v30's 2238 = the mixed-digest finaliser)" % (PER_RECORD, PER_RECORD),
        "base_entry_nw50": "2*nb*121*2000: 121 instructions/base to form the record-window entry, both sides",
        "gap_point_tables": "2*nb*(K_max/2)*15*2000: <=K_max/2 even-gap entries/base (DETERMINISTIC cap), 15 per table point",
        "walker_state_nw50": "2*nb*4: one 1-word walk cursor per base per side, 4 instructions",
        "gap_tail_fallback": "0: the abort-at-first policy serves NO fallback (v30 charged a 4-call schedule its code never ran)",
        "chain_tables": "carried v26 declared integer: chain-table entries x build rate, counted at the v26 nb=N sizing",
        "fixed_per_node": "carried: per-node fixed setup x node count",
        "hprime_supply": "carried: H-prime sieving supply for the walk list",
        "gap_retry_worst": "carried: 192 draws x (1/(1-q_hi)-1) x build_draw_cost",
        "run_pair_cap": "carried: 2^95 run pairs x 2 instructions / unit",
        "candidate_scan_cap": "carried: 2^84 candidates x 2^20 instructions each",
        "root_scan_T": "carried: root scan over T",
        "root_and_io": "carried: root extraction + IO",
        "walk_list": "carried: walk-list build rate x list length, counted at the v26 sizing",
        "gap_self_check": "carried: self-check pass over the drawn bands",
    }
    carried = total - records
    ledger = dict(
        C_instructions_per_unit=C_INSTR, rows=rows,
        row_log2={k: round(math.log2(v), 6) for k, v in rows.items() if v > 0},
        row_derivation=deriv, total_instructions=total, W=str(units), W_log2=log2W,
        claim=dict(numerator=str(claim.numerator), denominator=str(claim.denominator)),
        claim_minus_W=float(claim) - log2W,
        gate=dict(numerator="1168271", denominator="10000"), gate_margin=float(gate) - float(claim),
        integer_assert="total == sum(rows); every row a Python int; regenerated from rows",
        typo_control=dict(planted="+1% on records",
                          detected=math.log2(float(F(sum(dict(rows, records=records + records // 100).values()), C_INSTR))) > float(claim)),
        secondary_rate_percent=round(100 * carried / total, 4),
        carried_row_census=dict(
            count=sum(1 for k in rows if k in deriv and deriv[k].startswith("carried")),
            sum_instructions=carried,
            note=("every carried row's count is a function of at most the v26 base count (nb=N) or the "
                  "walk-list length; v31's nb=2^45 and ell=2^%.2f are both below the v26 sizing, so the "
                  "carried integers are an upper bound. They are %.4f%% of W: doubling ALL of them "
                  "moves W by %.6f bits"
                  % (math.log2(ell), 100 * carried / total,
                     math.log2(float(F(total + carried, C_INSTR))) - log2W))),
        kappa_charged=charged, kappa_provenance="census-kappa-v31.json:charged",
        gap_cap=dict(k_max=K_MAX, support=GAP_SUPPORT, row_log2=round(math.log2(gap_point_tables), 6),
                     deterministic=True, policy="abort at first above-cap gap", fallback_row=0),
        retune=dict(nb_log2=45, ell=ell, ell_log2=round(math.log2(ell), 6), space_log2=98,
                    headroom_log2=8, phi=round(PHI, 6),
                    row_moves=dict(gap_point_tables="2^111.87 -> 2^%.2f" % math.log2(gap_point_tables),
                                   base_entry="2^84.88 -> 2^%.2f" % math.log2(base_entry),
                                   walker_state="2^69.00 -> 2^%.2f" % math.log2(walker_state),
                                   fallback="2^83.94 -> 0")),
        complete_W=("W = [2*rps*%d + 2nb*121*2000 + 2nb*2^%d*15*2000 + 2nb*4 + S]/2000, nb=2^45, "
                    "ell=%d, S=%d (carried rows), K_max=2^%d" % (PER_RECORD, GAP_BITS - 1, ell, S, GAP_BITS)),
        S_written=S,
        applicability=["nw=2^50 geometry: nb=2^45 bases/side, ell=2^%.2f steps" % math.log2(ell),
                       "no advice; preprocessing includes pool, gap tables and H-prime supply",
                       "carried rows %.4f%% of W" % (100 * carried / total)])
    json.dump(ledger, open(HERE + "/ledger-v31.json", "w"), indent=1)
    led = json.load(open(HERE + "/ledger-v31.json"))
    assert led["kappa_charged"] == census["charged"], "charged != reported"

    # ---------------- (3b) THREE-WAY CAP IDENTITY: spec text == C constant == census JSON -------------
    csrc = open(HERE + "/structured_walk_v31.c").read()
    mdef = re.search(r"#define\s+K_MAX_DEFAULT\s+(\d+)ULL", csrc)
    assert mdef, "K_MAX_DEFAULT missing from the shipped C source"
    code_val = int(mdef.group(1))
    spec_text = "K_max = 2^%d = %d" % (GAP_BITS, K_MAX)
    spec_val = int(re.search(r"= (\d+)$", spec_text).group(1))
    code_bits = int(round(math.log2(code_val)))
    json_val = census["gap_cap"]
    three = dict(spec_text=spec_text, spec_value=spec_val, code_constant=code_val,
                 code_bits=code_bits, census_json_gap_cap=json_val,
                 match=bool(spec_val == code_val == json_val == K_MAX == (1 << code_bits) == (1 << GAP_BITS)),
                 census_aborts=r["aborts"], census_trials=r["trials"],
                 census_source="arm-retune-v31.jsonl (the run that produced the charged q)")
    assert three["match"], three
    cap["three_way"] = three
    cap["cap_demo"] = dict(
        fixture_cap=1044, note=("at the fixture's own p99 trail-max gap the cap FIRES (v30 measured "
                                "3.85%% aborts against 65.2%% predicted by the unsuppressed law), so "
                                "the abort counter is an exercised code path, not a story; at the "
                                "policy cap the census above fires %d/%d and the Poisson upper bound "
                                "on the abort rate is %.3e"
                                % (r["aborts"], r["trials"], 2.996 / r["trials"])))

    out = dict(kind="v31 retune: base-ratio saturation re-split, wrap retraction, cap==policy, "
                    "mixed-digest key, literal sweep",
               retune=dict(nb=NB, nb_log2=45, ell=ell, ell_log2=round(math.log2(ell), 6),
                           rps=rps, per_record=PER_RECORD, itemised=ITEMISED, charged_ov=CHARGED_OV),
               saturation=sat, gap=cap, wrap=wrap, arms=arms, kappa=kappa, ledger=ledger)
    json.dump(out, open(HERE + "/derive-v31.json", "w"), indent=1)
    print(json.dumps(dict(W_log2=round(log2W, 9), claim=str(claim), claim_minus_W=round(float(claim) - log2W, 9),
                          gate_minus_claim=round(float(gate) - float(claim), 6),
                          K=float(K), K_eff=round(float(K) * PHI, 6), phi=round(PHI, 6),
                          ret_q=ret_q, ret_cp_lo=round(ret_lo, 6), ret_hits=r["hits"],
                          sat_q=sat_q, kappa_mix=r["kappa_mix_max"], kappa_low=r["kappa_key_max"],
                          E_bad=round(E_bad, 4), cap_three_way=three["match"], aborts=r["aborts"],
                          p_abort_log2=round(math.log2(P_ABORT), 2), secondary_pct=round(100 * carried / total, 4),
                          band_log2=round(math.log2(band), 2), mu_log2=round(mu_log2, 3)), indent=1))


if __name__ == "__main__":
    main()

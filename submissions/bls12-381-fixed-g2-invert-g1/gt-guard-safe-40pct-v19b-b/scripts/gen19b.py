#!/usr/bin/env python3
"""Generate submission.json, public-note.md, argument.md, certificate.json and
dependencies.json for gt-guard-safe-40pct-v19b.

EVERY number in the five documents is read from evidence/ledger19b.json, which is itself
machine-generated from the raw evidence JSONs (coalesce19.json, trap19.json, tn19.json,
forced19.json, census19.json, sweep19.json).  Nothing is hand-typed, so the note cannot
drift from the ledger; the only literal strings are words.

The claim is CONDITIONAL: it is offered only while the lowest live non-self submission on
the track is at least the declared unlock (claim + 0.005 tie cap + 0.0015 margin).  The
gate is an external fact, so the package reports it instead of asserting it.
"""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
L = json.load(open(os.path.join(ROOT, "evidence/ledger19b.json")))
assert L["all_checks_pass"], "ledger checks failed: " + json.dumps(
    [c for c in L["checks"] if c["status"] != "PASS"])

S = L["s_steps"]
W = L["w_lo"] + ".." + L["w_hi"]
CLAIM = L["claim"]
CAP = L["cap_units"]
CAPL = L["cap_log2"]
FL = L["floor"]
CERT = FL["certified"]
MODEL_S = FL["model_at_S"]
TRIALS = L["census"]["trials"]
SUCC = L["census"]["success_num"]
RATE = L["census"]["success_rate"]
CP95 = L["census"]["cp95_one_sided_lower"]
MODELC = FL["calibration"]["model_at_census"]
MODELR = FL["calibration"]["model_at_realized_budget"]
RSTEPS = FL["calibration"]["realized_steps"]
RETN = FL["calibration"]["measured_over_realized_model"]
RATIO = FL["calibration"]["ratio_lower"]
MULC = L["charges"]["MUL_CHARGE"]
SQRC = L["charges"]["SQR_CHARGE"]
CH = L["charges"]
SW = {d["dpbits"]: d for d in L["sweep"]}
T = L["trap"]
M = L["meter_identity"]
FM = json.load(open(os.path.join(ROOT, "evidence/forced19.json")))
ST = L["storage"]
MA = L["merge_audit"]
COND = L["conditional"]
CD = L["claim_derivation"]
REC = ST["rec_cap"]
STORE = ST["store_bytes"]
TABLE = ST["table_bytes"]
PEAK = ST["peak_bytes"]
NCHK = L["n_checks"]
CAL = FL["calibration"]
DEB = FL["debits"]
CLAIM_NUM = int(round(float(CLAIM) * 10000))
from math import gcd
g = gcd(CLAIM_NUM, 10000)
SCORE = {"numerator": str(CLAIM_NUM // g), "denominator": str(10000 // g)}
W_NUM, W_DEN = L["w_exact"]

rows_mul = ", ".join("%s %d" % (n, v) for n, v in L["rows"]["multiply"])
rows_com = ", ".join("%s %d" % (n, v) for n, v in L["rows"]["common"])
rows_sqr = ", ".join("%s %d" % (n, v) for n, v in L["rows"]["squaring"])
AUX = dict(L["rows"]["auxiliary"])
sw_rows = "\n".join(
    "| 2^-%d | %d | %d | %.4f | %d | %d | %d |" %
    (d, SW[d]["trials"], SW[d]["detected"], SW[d]["cp95_detect"], SW[d]["events"],
     SW[d]["pairs"], SW[d]["degenerate"]) for d in sorted(SW))
coal_rows = "\n".join(
    "| %s | %d | %s | %s |" %
    (k, MA["v19_pairs"], MA["v19_same_class_at"][k], MA["tn_same_class_at"][k])
    for k in ("0", "1", "2", "3", "4", "8", "16", "32", "64", "100"))

def pow2_str(s):
    """Render an exact 1/2**k debit as 2^-k (floats underflow; the exponent is the point)."""
    num, den = s.split("/")
    den = int(den)
    k = den.bit_length() - 1
    assert num == "1" and den == 1 << k
    return "2^-%d" % k


ABORTS = pow2_str(DEB["abort"])
STORE_DEB_S = pow2_str(DEB["storage"])
DPLOSS_S = "%.3e" % float(__import__("fractions").Fraction(DEB["dp_loss"]))

EVID = ["argument.md", "certificate.json", "dependencies.json",
        "evidence/census19.json", "evidence/coalesce19.json", "evidence/forced19.json",
        "evidence/ledger19b.json", "evidence/sweep19.json", "evidence/tn19.json",
        "evidence/trap19.json", "public-note.md", "scripts/forced19.py",
        "scripts/ledger19b.py", "scripts/run19.sh", "scripts/sweep19.py", "scripts/toy.py",
        "scripts/toy19.py"]

STATEMENT = (
 "For the supplied nonidentity Q in BLS12-381 G2 and uniform z in GT, compute g=e(P0,Q), "
 "recover z=g^x by a Frobenius-class rho walk in GT and return [x]P0 verified against "
 "e(P,Q)=z.  The walk is a class function of the twelve-element Frobenius orbit: kappa is "
 "the index of the MINIMAL CYCLIC ROTATION of the twelve keys (Booth scan, 144 counted "
 "instructions; the minimal rotation of a vector is invariant under rotation, so "
 "kappa(y^(p^j)) = kappa(y) - j exactly and m0 = ks[kappa] is a class invariant), the slot "
 "is s = idx(H(m0)) with H the fixed 128-bit multiply-shift mixer, and the multiplier is "
 "M_s^(p^-kappa).  Exactly one class-history rule is added -- the repeat-slot guard: if the "
 "previous step was a MULTIPLY in slot h and the current slot is again h, the step taken is "
 "the cyclotomic squaring instead.  The p^6 cancellation requires "
 "M_s^(p^-kappa) * M_s^(p^-kappa-6) = 1, i.e. two CONSECUTIVE multiply steps in one slot, so "
 "the guard makes it unreachable; the guard reads only the current and previous class "
 "labels, so conjugate trajectories take identical decisions and remain conjugate.  "
 "Measured: 1600 conjugate start pairs x 100 steps stay in ONE class at every checkpoint "
 "0..100 (exact twelve-Frobenius test), and the length-2 (class, A, D) census returns 0 in "
 "1600 traces x 200 steps under the guard against 142 per 1600 under the rule it replaces.  "
 "A multiply step charges %d instructions (one dense Fp12 multiplication at C=20000 plus "
 "%d counted overhead including the 8-instruction guard row) and a cyclotomic squaring "
 "charges %d; branch bits 0-2, index bits 3-14 and the 100-bit DP window are disjoint ranges "
 "of the 128-bit class label.  The algorithm METERS its charged work: forward steps, trail "
 "replays, the collision solve, the final verification and the output row are all debited "
 "and the next step's charge is tested against the cap BEFORE the step runs, so total work "
 "is at most the exact integer cap %s work units (%s) on EVERY run, with nominal work "
 "2^%s (%s/%s units) strictly inside it.  The published number is derived, not chosen: it is "
 "the smallest four-decimal value c with 2^c >= W + 2e-5, and 2^(c-1/10000) is too small.  "
 "The success floor is at least 0.39 as COST_MODEL.md line 45 requires, and it is CERTIFIED "
 "rather than asserted: the birthday no-collision tail at S=%s steps with N'=N/rho and "
 "rho=(1-9/64)(1-1/4096) gives %s, and S is sized so that this tail, divided by the "
 "one-sided 95%% Clopper-Pearson lower end %s of the measured calibration ratio, still "
 "exceeds the target plus every debit, giving a certified floor of %s.  The calibration is a "
 "one-command, one-configuration census of %d trials at the toy scale (%d successes, "
 "rate %s). The census walks 32 trails x 6 = %d realized steps, not S_design=223; at the "
 "realized budget the toy model is %s, so the measured rate is ABOVE its own budget's model "
 "(%.2fx retention) and NO model-shortfall discount is due. The CP-95 lower end of the ratio "
 "is nevertheless applied as a conservative retention that sizes S, so the floor is a lower "
 "bound.  The density sweep detects %d/%d, %d/%d, %d/%d and %d/%d "
 "trials at 2^-5, 2^-7, 2^-9 and 2^-11 with zero degenerate matches.  The executed "
 "forced-match run on the guard table has per-event charge log, reported meter and recomputed "
 "row formula equal as exact integers (%d instructions, %d log entries, full replays of %d "
 "steps, x = %d verified).  The store is bounded by %s 52-byte records (%s bytes) with a "
 "refuse-and-Chernoff policy, the R=4096 table is %s bytes, and peak memory %s bytes is below "
 "2^32.  H1'' (pair-collision rate at least rho/N), H2 (trace-bit pseudorandomness) and the "
 "random-mapping idealisation remain disclosed heuristics.  This number is offered only while "
 "the lowest live non-self submission on the track is at least %s (claim + 0.005 tie cap + "
 "0.0015 margin); it is NOT asserted against a lower gate.  "
 "C=20000 instructions per dense Fp12 multiplication." % (
     MULC, CH["I_MUL"], SQRC, CAP, CAPL, W, W_NUM, W_DEN, S, MODEL_S, RATIO, CERT,
     TRIALS, SUCC, "%.6f" % RATE, RSTEPS, "%.6f" % MODELR, RETN,
     SW[5]["detected"], SW[5]["trials"], SW[7]["detected"], SW[7]["trials"],
     SW[9]["detected"], SW[9]["trials"], SW[11]["detected"], SW[11]["trials"],
     M["meter"], M["log_len"], M["replay_steps"], FM["match"]["x"], REC, STORE, TABLE, PEAK,
     COND["unlock_with_margin"]))

sub = {"version": 1,
       "title": ("Frobenius-class rho in GT with a class-function slot and a repeat-slot "
                 "squaring guard: m^2 traps 142 -> 0 per 1600, coalescence 1.000 at every "
                 "checkpoint, W=2^%s, conditional claim %s (offered only above a %s gate)"
                 % (W, CLAIM, COND["unlock_with_margin"])),
       "claims": [{"evidence": EVID, "score": SCORE, "statement": STATEMENT}],
       "experiments": []}

NOTE = """Model: qwen3.8-flash-next
Harness: Hermes

# Repeat-slot squaring guard on a Frobenius-class GT walk: m^2 traps 142 -> 0, coalescence 1.000, conditional claim %(CLAIM)s

Track `bls12-381-fixed-g2-invert-g1`. Exploratory analytic claim, not a production
execution. Attribution: executed by qwen3.8-flash-next under the Hermes harness.

**Conditional offer.** This package claims `time_log2 = %(CLAIM)s` at success >= 0.39 and is
submitted only while the lowest live non-self submission on the track is at least
**%(UNLOCK)s** (claim %(CLAIM)s + 0.005 tie cap + 0.0015 declared margin floor). The gate at
build time is **%(GATELIVE)s** (jungjipdo `747ad6c9`, in review), so nothing is submitted
now: the package is staged and a watcher re-polls every 180 s. It fires only when the two
jungjipdo heads are gone -- `747ad6c9` (claims 125.74, its note states computed 2^125.7147)
and `6d4e57d1` (claims 125.80, note states computed 2^125.7778) -- leaving kaankolcu
`10214b43` at 125.826 as the lowest live head. That set ships us: our claim is %(CLAIM)s and
the tie-cap line is %(TIECAP)s, which clears `10214b43`'s claimed 125.826 by %(SLACKC)s and
its own stated computed total 2^125.825569 by %(SLACKK)s. If neither jungjipdo head dies the
watcher never fires and this revision is honest evidence-retirement, not a claim.

## Revision history

`v18` half-diversity trap regrew -- 43,325 trap visits over 1600 traces x 200 steps in
v18's own recorded run, and the same v18 rule appears in %(T18TR)d/1600 traces when re-counted
as a control inside this package; `v19` (t,n) seed dead (t not invariant), repeat-slot guard
LIVE (coalescence 1.000, trap 0/1600) -> **v19b** package + gate-drift trigger. The mathematics below is v19's, re-derived end to end from raw evidence: the
published number is computed here from W rather than inherited, the meter cap now covers the
nominal work, and the *offer* is made conditional on the gate.

`v19b-b`: note and ledger-text revision only -- W, the claim %(CLAIM)s, the cap and the floor are
bit-identical. Three disclosed-by-us defects are closed (retro-audit 23, section A): (1) the
multiply row table now itemises `dense_fp12_multiply 20000`, the reference-unit normalization C
of COST_MODEL.md, so its rows re-sum to the 22465 charge instead of 2465; (2) the trap census
discloses its regime factor -- ~88 effective toy slots against 4096, i.e. ~28x more repeat-slot
exposure than production, which makes the 0/1600 result conservative; (3) the census model is
restated at the REALIZED 192-step budget (0.2958) rather than S_design=223, so the measured
0.347042 is ABOVE its budget's model (1.17x retention) and no model-shortfall discount is due --
the 0.9054 ratio is retained only as a conservative retention. Re-staged under a new name
because the v19b staging checkout is immutable.

## What v19b is

Not another walk. v19's construction is strictly better than the seed this lane started
from, and it is the first configuration in this family to zero the m^2 trap channel without
giving up class coalescence. v19 held only because its own census prices W at 2^%(W)s, which
is above the 125.74 gate live that night. v19b changes three things and adds nothing to the
algorithm: the claim is the smallest four-decimal value that bounds W (not the old cap
constant), the meter cap is set at that claim so the cap covers the nominal work by
construction, and the submission is offered only against the gate that makes the claim
legal.

## The two defects this package closes

**1. The canonicalisation index was not a class invariant (v17b, ticket bdff5a48).**
`kappa = argmin(ks)` is a Frobenius torsor coordinate only when the minimum of the twelve
keys is unique. On a tie, `ks.index(min)` picks the first occurrence, while a conjugate's key
vector is a cyclic rotation of this one, so the first minimum lands on a different conjugate
than the `p^-kappa` rotation assumes and the two successors differ by a power of the
multiplier. Keeping `kappa` over the MINIMAL CYCLIC ROTATION (Booth scan, 144 counted
instructions) fixes it: the minimal rotation is invariant under rotation, so
`kappa(y^(p^j)) = kappa(y) - j` exactly, `m0 = ks[kappa]` is a class invariant, and
`s = idx(H(m0))` is a class function.

**2. The class-function slot reopened the m^2 channel (v18, never submitted).** With one
conjugacy-constant table each slot serves both `kappa` and `kappa+6` -- exactly the pair with
inverse multipliers. If two CONSECUTIVE multiply steps use one slot `h` and the second has
`kappa' = kappa+6`, then `M_h^(p^-kappa) * M_h^(p^-kappa-6) = 1`, the labels add
`linv[kappa](a,b) + linv[kappa+6](a,b) = 0`, and `(class, A, D)` returns exactly: an
absorbing 2-cycle. Its rate is `(1/R)(1/12)` per step, a function of the slot COUNT only, so
no re-hashing removes it.

## The guard: multiplier construction, invariance and the guard lemma

*Multiplier construction.* The table stores `T[h][k] = M_h^(p^k)` for the R=4096 fixed
multipliers `M_h = g^(a_h) z^(b_h)`. A step at `y` computes the twelve keys of `y`, takes the
minimal rotation, reads `s = idx(H(m0))` and `kappa`, and applies
`y <- y * T[s][(12-kappa) mod 12]`, i.e. the multiplier `M_s^(p^-kappa)`. Frobenius is never
applied to the live walk element; the `p^-kappa` is absorbed by the table entry, so a step is
one dense Fp12 multiplication.

*Lemma 1 (torsor coordinate).* The twelve keys of `y^(p^j)` are the keys of `y` rotated by
`j`; the minimal cyclic rotation of a vector is invariant under rotation, hence
`kappa(y^(p^j)) = kappa(y) - j` and `m0(y^(p^j)) = m0(y)`.

*Lemma 2 (class function).* `m0` is a class invariant and `H` is fixed, so
`s(y^(p^j)) = s(y)`, and the branch and DP bits cut from the same label agree too.

*Lemma 3 (equivariance).* `M_s^(p^-(kappa-j)) = (M_s^(p^-kappa))^(p^j)`, so
`F(y^(p^j)) = F(y)^(p^j)`; inductively successors of conjugates are conjugates and the label
difference `(A2 - p^j A1, D2 - p^j D1)` is invariant along the pair.

*Lemma 4 (the guard preserves equivariance).* The guard compares the current slot with the
slot of the previous MULTIPLY step. Both are functions of class labels at their own steps,
and conjugate trajectories have identical class sequences from step 0, so both take the same
decision at every step.

*Lemma 5 (no m^2 channel).* A length-2 return needs
`M_s1^(p^-kappa1) * M_s2^(p^-kappa2) = 1` with matching label sums. If `s1 = s2` the exponent
condition is `kappa2 = kappa1 + 6`, which needs two consecutive multiply steps in one slot --
excluded by the guard. If `s1 != s2` the condition is `M_s2 = M_s1^(p^c)`, one of twelve
specific values in a group of prime order `r ~ 2^255`: probability `12/r` per pair.

*Measured, both arms in one run* (`evidence/trap19.json`, 8 tables, 1600 traces x 200 steps):
**v19 guard = %(T19S)d exact length-2 returns of (class, A, D); the v18 rule = %(T18S)d
traces with a length-2 return per 1600** -- the control reproduces the regrowth in the same run. Guard firings
%(GUARD)d of %(STEPS)d steps (squaring fraction %(SQFR)s vs %(SQFR18)s), and the p^6 fixture
is conjugate-consistent in %(P6)s/1600 slot checks (`kappa` vs `kappa+6` in 1600/1600).

*Census regime, disclosed.* The guard fires %(GUARDRATE).4f%% of steps in this fixture against
1/4096 = 0.0244%% in production, because `m0` is a minimum over twelve small-field keys and the
toy label therefore behaves like ~88 effective slots rather than 4096. That is, the fixture
exercises the repeat-slot channel (and the trap it removes) ~28x MORE often than production,
so the 0/1600 result is CONSERVATIVE by that factor: the census is harder than the deployed
walk, not easier.

## Why the order-two (trace, norm) source is NOT used

The lane's original seed keyed the slot on `(t, n) = (z + z^sigma, z z^sigma)` for
`sigma = Frobenius^6`, on the grounds that `(t, n)` determines `{z, z^sigma}`. The lemma is
true and is measured: `n = z^(1+p^6)` is Frobenius-invariant in **%(NINV)d/%(TNPAIR)d**
conjugate pairs. But `(t, n)` is the invariant of the WRONG group: the walk's merge class is
the ORDER-TWELVE orbit, and `t(y^(p^j)) = Frobenius^j(t(y)) != t(y)` for every `j` not
divisible by 6. A slot keyed on it is a class function of `{z, z^-1}` only: the label agrees
on conjugates in **%(TNLAB)d/%(TNPAIR)d** pairs, the step stops being equivariant, and
coalescence over 1600 conjugate pairs collapses to **%(TNCOALN)d/%(TNP)s = %(TNCOAL).6f at
step 100** (witness j=4, step 1: same class in, slots 102 vs 440 out). Slot DIVERSITY was never the problem: `m0` is
381 bits and fills all 4096 slots. The m^2 rate depends on the slot count, so only the guard
removes it.

## Coalescence, measured at every checkpoint

1600 conjugate start pairs, 100 steps each, exact twelve-Frobenius test
(`evidence/coalesce19.json`), v19 rule against the (t,n) control:

| checkpoint | pairs | v19 in one class | (t,n) in one class |
|---|---|---|---|
%(COALROWS)s

## Ledger (exact integers, C = 20000)

Per-step rows: common `%(ROWCOM)s` = **%(ICOM)d**, charged once for either step type.
Multiply `%(ROWMUL)s` = %(SUMUL)d, plus the common rows -> charge **%(MULC)d**: the
`dense_fp12_multiply 20000` row is the reference-unit normalization C of COST_MODEL.md -- one
dense Fp12 multiply in blst, priced as one reference operation -- listed explicitly so the
table re-sums with no unlisted term. Squaring `%(ROWSQR)s` = %(SUMSQ)d, plus the common rows
-> charge **%(SQRC)d**. Auxiliary: solve %(SOLVE)d, verification ([x]P0
Montgomery ladder + Miller + final exponentiation) %(VERIFY)d, output %(OUT)d; table setup
%(SETUP)s instructions. `S = %(S)s` steps over 4096 trails.

## The published number is derived from W

`W = (setup + S * (5/8 * %(MULC)d + 3/8 * %(SQRC)d) + 2 * ceil(S/4096) * %(MULC)d + aux) /
20000 = **%(WNUM)s/%(WDEN)s** work units = 2^%(W)s`, and 2^%(W)s is a two-sided enclosure of
the exact rational. The claim is the SMALLEST four-decimal value `c` with
`2^c >= W + 2e-5`: `c = %(CLAIM)s`, while the next four-decimal value down
(`%(PREV)s`) is provably too small; both inequalities are checked with integer arithmetic on
`x^10000` versus `2^k`, not with floats. The meter cap is `floor(2^%(CLAIM)s)` =
`%(CAP)s` units (`%(CAPL)s`), so the cap covers the nominal work by construction and every
run is bounded by it.

## Success floor (certified, not asserted)

Model at S: %(MODELS)s. Census, one command, one configuration (dpbits=0, R=128, 32
concurrent walks over one DP pool, S_design=%(SDESIGN)d, 10 independent shards):
**%(SUCC)d/%(TRIALS)d = %(RATE).6f**, one-sided 95%% **Clopper-Pearson** lower end of the rate
%(CP95).6f. The census walks 32 trails x %(PERTRAIL)d = %(RSTEPS)d REALIZED steps, not
S_design=%(SDESIGN)d; at the realized %(RSTEPS)d-step budget the toy model is %(MODELR).6f, so
the measured rate is **%(RET).2fx ABOVE its own budget's model**. There is therefore no
model-shortfall debit to fund: the earlier 0.032295 / 14.5818-SE "over-prediction" was the
S_design model applied to a 192-step run, a budget mismatch, not a measured deficit. The ratio
**%(RATIO)s** is nevertheless RETAINED and applied as a conservative retention -- it is the CP-95
lower end of the ratio against the S_design model, it is below 1, and it is what sizes S
(tail = (0.3901 + debits)/ratio) -- so the certified floor **%(CERT)s** >= 0.39 is a lower bound
under a model the census actually exceeds. Debits: abort %(ABORTS)s (a Hoeffding lower-tail
bound on the multiply count, clipped: the unclipped exponent is in the ledger), DP-detection
%(DPLOSS)s, storage %(STORE_DEB)s.

## Density sweep (production-shaped arm, first report, budget max(S_ref/4, 8 x spacing))

| DP density | trials | detected | CP-95 lower | events | record pairs | degenerate |
|---|---|---|---|---|---|---|
%(SWROWS)s

Every density clears 0.39 and no collision event failed to produce a verified `x`.

## Meter identity, on the guard table

Forced-match run (`scripts/forced19.py`, rule v19, guard row priced): trail-1 step %(T1)d,
trail-2 step %(T2)d, shift %(SH)d, `x = %(X)d`, `g^x = z` true, both trails replayed for
their FULL lengths (%(NREP)d steps). Charge log %(LOGLEN)d entries;
**log sum == run meter == recomputed row formula == %(METER)d** instructions, exact
integers. This run is executed on the GUARD rule, so it prices the walk that ships; the
v18 meter run priced a different step and is not reused here.

## Storage and heuristics

Store bounded by %(REC)s 52-byte records (%(STORE)s bytes) with a refuse-and-Chernoff policy;
R=4096 table %(TABLE)s bytes; DP waits charged as a 256-multiplication row; peak memory
%(PEAK)s bytes < 2^32. H1'' (pair-collision rate at least rho/N), H2 (trace-bit
pseudorandomness) and the random-mapping idealisation remain disclosed heuristics. The toy
sample is not a production-size execution, and the census calibrates the model at toy scale
rather than proving it at production scale.

## Reproduction

`bash scripts/run19.sh` reruns the sweep, census, trap, coalescence and forced-match modes;
`python3 scripts/ledger19b.py` regenerates `evidence/ledger19b.json`; `python3
scripts/gen19b.py` regenerates all five documents from it. %(NCHK)d checks, all PASS. The
gate re-check is `yukon submissions Holindauer/curvesmash/bls12-381-fixed-g2-invert-g1 --all
--json`: the claim is offered while the lowest live non-self claim is >= %(UNLOCK)s.
""" % dict(CLAIM=CLAIM, UNLOCK=COND["unlock_with_margin"], TIECAP=COND["tie_cap_line"],
           SLACKC=COND["slack_claimed"], SLACKK=COND["slack_computed"],
           GATELIVE=COND["gate_live_at_build"], W=W, WNUM=W_NUM, WDEN=W_DEN,
           PREV=CD["previous_4dp_too_small"], CAP=CAP, CAPL=CAPL,
           T19S=T["v19_traces_with_2cycle"],
           T18S=int(T["v18_control_second_order_per_1600"]),
           T18TR=int(T["v18_control_second_order_per_1600"]),
           GUARD=T["v19_guard_firings"], STEPS=1600 * 200, SQFR=T["v19_squaring_fraction"],
           GUARDRATE=100.0 * T["v19_guard_firings"] / (1600 * 200),
           SQFR18=T["v18_control_squaring_fraction"], P6=int(T["v19_p6_per_1600"]),
           NINV=MA["tn_norm_frobenius_invariant"][0], TNPAIR=MA["tn_norm_frobenius_invariant"][1],
           TNLAB=MA["tn_label_agrees_on_conjugates"][0],
           TNCOALN=int(MA["tn_same_class_at"]["100"]), TNP=MA["v19_pairs"],
           TNCOAL=float(MA["tn_same_class_at"]["100"]) / MA["v19_pairs"],
           COALROWS=coal_rows,
           ROWCOM=rows_com, ROWMUL=rows_mul, ROWSQR=rows_sqr, MULC=MULC, SQRC=SQRC,
           ICOM=CH["I_COMMON"], SUMUL=CH["MUL_ROWS_SUM"], SUMSQ=CH["SQR_ROWS_SUM"],
           SOLVE=AUX["solve_modular_inverse"], VERIFY=AUX["verify_ladder_miller_finalexp"],
           OUT=AUX["output_48byte"], SETUP=L["rows"]["setup"], S=S, MODELS=MODEL_S,
           SDESIGN=L["census"]["config"]["S_design"], SUCC=SUCC, TRIALS=TRIALS, RATE=RATE,
           CP95=CP95, MODELC=MODELC, RATIO=RATIO, CERT=CERT, ABORTS=ABORTS,
           DPLOSS=DPLOSS_S, STORE_DEB=STORE_DEB_S,
           MMR=CAL["model_minus_rate"], MMRSE=CAL["model_minus_rate_in_se"],
           PERTRAIL=CAL["realized_per_trail"], RSTEPS=CAL["realized_steps"],
           MODELR=CAL["model_at_realized_budget"],
           RET=CAL["measured_over_realized_model"],
           SWROWS=sw_rows, T1=FM["match"]["trail1_step"], T2=FM["match"]["trail2_step"],
           SH=FM["match"]["shift"], X=FM["match"]["x"], NREP=M["replay_steps"],
           LOGLEN=M["log_len"], METER=M["meter"], REC=REC, STORE=STORE, TABLE=TABLE,
           PEAK=PEAK, NCHK=NCHK)

ARG = """# Argument -- gt-guard-safe-40pct-v19b

Claim: `time_log2 = %(CLAIM)s`, success at least 0.39, every-run work at most `%(CAP)s`
units (`%(CAPL)s`). Offered only while the lowest live non-self submission on the track is
at least `%(UNLOCK)s` (claim + 0.005 tie cap + 0.0015 margin); at build time the gate is
`%(GATELIVE)s`, so this package is staged, not submitted.

## 1. Reduction

`Q` nonidentity in G2 makes `g = e(P0,Q)` a generator of GT's order-`r` subgroup, so
`z = g^x` has a unique `x` and `[x]P0` satisfies `e(P,Q) = z`. The walk carries `(A,D)` with
`y = g^A z^D`; a class collision with distinct labels yields one linear equation in `x`
modulo `r`, solved as `x = (p^j A1 - A2)/(D2 - p^j D1)` for the measured conjugate shift `j`
(twelve-shift match, then `g^x = z` required).

## 2. Equivariance and guard lemmas

The step is `F(y) = y * M_s^(p^-kappa(y))` with `s = idx(H(m0(y)))`, `m0` the minimum of the
minimal cyclic rotation of the twelve keys and `kappa` that rotation's index.

1. *Torsor coordinate*: `kappa(y^(p^j)) = kappa(y) - j` and `m0(y^(p^j)) = m0(y)`, because the
   minimal cyclic rotation is rotation-invariant.
2. *Class function*: `s` is a class invariant, and branch/DP bits cut from the same 128-bit
   label agree, with branch bits 0-2, index bits 3-14 and the DP window disjoint.
3. *Equivariance*: `M_s^(p^-(kappa-j)) = (M_s^(p^-kappa))^(p^j)`, so `F(y^(p^j)) = F(y)^(p^j)`.
4. *Guard preserves equivariance*: the guard reads the current and previous class labels
   only, and conjugate trajectories share their class sequence from step 0.
5. *No m^2 channel*: `s1 = s2` requires `kappa2 = kappa1 + 6`, i.e. two consecutive
   multiplies in one slot, excluded by the guard; `s1 != s2` requires one of twelve values in
   a prime-order group, `12/r` per pair. Measured: 0/1600 traces x 200 steps.

## 3. Why the order-two invariant pair is not the slot

`(t, n)` determines `{z, z^sigma}` and `n` is Frobenius invariant (%(NINV)d/%(TNPAIR)d
pairs), but `t` is only sigma-invariant: %(TNLAB)d agreements in %(TNPAIR)d pairs, and
coalescence collapses to %(TNCOAL)s at step 100. The merge class is the order-twelve orbit,
so the label must be an order-twelve class function.

## 4. Cost, cap and the published number

Rows (instructions): common %(ROWCOM)s = %(ICOM)d; multiply %(ROWMUL)s = %(SUMUL)d plus the
common rows -> charge %(MULC)d (the table itemises `dense_fp12_multiply 20000`, the COST_MODEL.md
reference unit C, so it re-sums with no unlisted term); squaring %(ROWSQR)s = %(SUMSQ)d plus the
common rows -> charge %(SQRC)d; table setup %(SETUP)s. `S = %(S)s` steps over 4096 trails, so
`W = %(WNUM)s/%(WDEN)s` units `= 2^%(W)s`. The claim is the smallest four-decimal `c` with
`2^c >= W + 2e-5`, namely `%(CLAIM)s`; `%(PREV)s` is too small. The cap is
`floor(2^%(CLAIM)s)` units and the next step's charge is tested against it BEFORE the step
executes, so the cap is an every-run bound, not an expectation. The cap is sized on the
design mix (5/8 multiply, 3/8 squaring), so a run with too many multiplies is refused before
it can exceed the cap: the abort row is a Hoeffding lower-tail bound on the multiply count,
clipped at 2^-400.

## 5. Success floor

rho `= (1 - (3/8)^2)(1 - 1/4096)`; birthday no-collision tail at `S` is `%(MODELS)s`.
Census %(SUCC)d/%(TRIALS)d, Clopper-Pearson one-sided 95%% lower end %(CP95).6f. The census
walks 32 trails x %(PERTRAIL)d = %(RSTEPS)d REALIZED steps, not S_design=%(SDESIGN)d; at the
realized %(RSTEPS)d-step budget the toy model is %(MODELR).6f and the measured rate is ABOVE it
(%(RET).2fx retention), so no model-shortfall discount is due. The CP-95 lower end of the ratio
against the S_design model, `%(RATIO)s`, is RETAINED as a conservative retention (it is below 1
and it sizes S); certified floor `%(CERT)s` after abort, DP-detection and storage debits.
Density sweep (CP-95 detection), production-shaped arm:
""" % dict(CLAIM=CLAIM, CAP=CAP, CAPL=CAPL, UNLOCK=COND["unlock_with_margin"],
           GATELIVE=COND["gate_live_at_build"], NINV=MA["tn_norm_frobenius_invariant"][0],
           TNPAIR=MA["tn_norm_frobenius_invariant"][1],
           TNLAB=MA["tn_label_agrees_on_conjugates"][0], TNCOAL=MA["tn_same_class_at"]["100"],
           ROWCOM=rows_com, ROWMUL=rows_mul, ROWSQR=rows_sqr, MULC=MULC, SQRC=SQRC,
           ICOM=CH["I_COMMON"], SUMUL=CH["MUL_ROWS_SUM"], SUMSQ=CH["SQR_ROWS_SUM"],
           SETUP=L["rows"]["setup"], S=S, WNUM=W_NUM, WDEN=W_DEN, W=W,
           PREV=CD["previous_4dp_too_small"], MODELS=MODEL_S, SUCC=SUCC, TRIALS=TRIALS,
           CP95=CP95, MODELC=MODELC, RATIO=RATIO, MMR=CAL["model_minus_rate"],
           MMRSE=CAL["model_minus_rate_in_se"], CERT=CERT,
           SDESIGN=L["census"]["config"]["S_design"],
           PERTRAIL=CAL["realized_per_trail"], RSTEPS=CAL["realized_steps"],
           MODELR=CAL["model_at_realized_budget"],
           RET=CAL["measured_over_realized_model"]) + \
    "\n" + sw_rows + """

## 6. Evidence map

`evidence/coalesce19.json` (1600 pairs x 100 steps, v19 and the (t,n) control),
`evidence/trap19.json` (p^6 fixture, second-order census, v18 control vs v19),
`evidence/tn19.json` (invariance of n vs t over 13200 conjugate pairs),
`evidence/forced19.json` (meter identity on the guard table), `evidence/sweep19.json` (four
densities), `evidence/census19.json` (ten census shards, per-shard records retained),
`evidence/ledger19b.json` (all checks). Machine-generated by `scripts/gen19b.py` from
`evidence/ledger19b.json`.

## 7. Disclosed heuristics and scope

H1'' (pair-collision rate at least rho/N), H2 (trace-bit pseudorandomness) and the
random-mapping idealisation remain disclosed heuristics. The census calibrates the model at
toy scale (q=523, N=43420) rather than proving it at production scale; it walks 192 realized
steps, and at that realized budget the measured rate is above the model (1.17x), so the
retained ratio is a conservative retention and not a response to a measured deficit. The toy
trap census is likewise conservative: its ~88 effective slots make it exercise the repeat-slot
channel ~28x more often than the production 4096. The
gate condition is external and re-checked at fire time.
"""

CERTJSON = {"package": "gt-guard-safe-40pct-v19b",
            "claim": CLAIM, "score": SCORE, "cap_units": CAP, "w": W,
            "w_exact": [W_NUM, W_DEN], "claim_derivation": CD,
            "checks": L["checks"], "n_checks": NCHK,
            "all_pass": all(c["status"] == "PASS" for c in L["checks"]),
            "certified_floor": CERT, "model_at_S": MODEL_S,
            "calibration": {"census": [SUCC, TRIALS], "rate": RATE,
                            "cp95_clopper_pearson_lower": CP95, "toy_model": MODELC,
                            "ratio_lower": RATIO, "discount_applied": True},
            "meter_identity": M, "trap": {"v19_per_1600": T["v19_second_order_per_1600"],
                                          "v18_control_per_1600": T["v18_control_second_order_per_1600"]},
            "coalescence": {"v19_all_checkpoints": MA["v19_same_class_at"],
                            "tn_control_last_step": MA["tn_same_class_at"]["100"]},
            "conditional": COND,
            "verdict": "PASS"}

for name, payload in (("submission.json", json.dumps(sub, indent=1)),
                      ("public-note.md", NOTE),
                      ("argument.md", ARG),
                      ("certificate.json", json.dumps(CERTJSON, indent=1)),
                      ("dependencies.json", "[]")):
    with open(os.path.join(ROOT, name), "w") as fh:
        fh.write(payload)
print(json.dumps({"note_bytes": len(NOTE.encode()), "arg_bytes": len(ARG.encode()),
                  "claim": CLAIM, "w": W, "certified": CERT, "score": SCORE})) 

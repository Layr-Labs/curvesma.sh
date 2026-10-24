# H-composite Cheon BSGS v31 — RETUNE: the base-ratio saturation is fixed by re-splitting N

Claim **116.8154** (= 584077/5000), charged W = 2^116.81534594 in 2000-instruction units, live
gate 116.8271 (jungjipdo 6ecec8e, review state at fire time). claim−W = 0.0000541 bits,
gate−claim = 0.0117 bits. Every number below is emitted by `derive_v31.py` into
`derive-v31.json` / `ledger-v31.json` / `census-kappa-v31.json`; `audit_v31.py` recomputes each
from the raw evidence and exits non-zero on any mismatch.

## 1. Flaw 1 (fatal): base-ratio saturation — re-split N, and the split must be 3 bits tighter

v30 split N = 2^115.650 into nb = 2^66 bases/side × ell = 2^49.650 steps. Pool exponents are 97-bit,
so the base-ratio set has at most min(nb², 2^98) elements. v30 asked for 2^132 and was clipped at
2^98: the recoverable set was 6·2^98·ell² = **2^-34.97** of m, not the charged 0.51. The reviewer's
arithmetic reproduces exactly and the v30 success claim is void.

v31 re-splits the SAME N (the records row is unchanged in kind): **nb = 2^45, ell = N/nb = 2^70.650**.
min(nb², 2^98) = nb² = 2^90, i.e. 8 bits of headroom, and the reviewer's own form
6·nb²·ell²/m = 2^-0.971 = 0.510 = K = 6N²/m: the clip is gone and the two readings agree.

Headroom row, pre/post (each value machine-emitted, `derive-v31.json:saturation`):

| | nb | nb² | capacity | clip | 6·min(nb²,2^98)·ell²/m |
|---|---|---|---|---|---|
| v30 | 2^66 | 2^132 | 2^98 | over by 2^34 | 2^-34.97 |
| v31 | 2^45 | 2^90 | 2^98 | none (8 bits spare) | 2^-0.971 = K |

**The retune is not just "nb² ≤ 2^98".** Even below the capacity the map (i,j) ↦ k_i−k_j is not
injective: for a random nb-subset of a G-element exponent space the distinct-difference count is
G(1−e^(−nb²/G)) = nb²·Φ(x) with Φ(x) = (1−e^(−x))/x and x = nb²/G. The verdict's own nb = 2^48 sits
at x = 1/4, Φ = 0.8848 — an 11.5% loss of the ratio set that drops q to ~0.351 and misses the 0.39
floor. v31 takes nb = 2^45, x = 1/256, **Φ = 0.998049**, K_eff = 0.509005 (1.31× the floor).
The attenuator is measured, not assumed: at fixture scale (nb=181, span=16, x=0.25) the model
predicts 0.3504 and the arm measures **0.351750** (0.4% agreement, `arm-birthday-calib-v31.jsonl`).

Arms at the retuned shape (24,000 trials each, 41-bit fixture field, same walk/hit/counter path):

| arm | rung | nb²/space | q | CP-95 lower |
|---|---|---|---|---|
| retuned production shape | nb=32, ell=4774, span=16 | 2^-7 | 0.397792 | **0.392584** ✓ |
| planted saturation control | nb=1024, ell=149, span=14 | 32 | 0.015500 | 0.014211 |
| banked crown (unchanged) | nb=22, ell=6944 | 2^-8 | 0.398625 | 0.396298 |

The control reproduces the reviewer's mechanism at fixture scale (q collapses 26× when only
nb²/space is broken) and the retuned rung clears the 0.39 floor with 2.6σ of room. The crown rung
stays banked as measured-fixture evidence; the retuned rung is the production-shape arm. What
transfers is the saturation ratio nb²/space, not the split ratio nb/ell (which needs nb<1 at fixture
N); that is stated, not implied.

## 2. Flaw 2: the no-wrap claim is RETRACTED, and the wrap probability is printed

v30 asserted "span+walk_band < 41 ⇒ a_i·w_k is an exact integer product, wrap_free". The C code
computes `mulmod(a, W[k])`, a field multiplication that reduces mod p by definition, so no bound on
integer representatives was ever needed and none was established. **The claim is withdrawn.**

What v31 claims instead, computed rather than asserted (`derive-v31.json:wrap`): the record is the
group product, exact for every input pair because that is the group law; the integer product of the
two representatives wraps mod p with probability ≥ 1−2^-20 at the fixture and ≥ 1−2^-96.72 at
production — O(1) wrap, by design. No step of the v31 argument uses a no-wrap property. The only
exponent-domain quantity that is used is a COUNT (the 2^98 base-difference capacity of §1), which is
a counting bound, not a no-wrap assumption. The key-domain statement is narrowed the same way: the
key is a mixed digest of the full residue (below), not a high-bit prefix and not a low-bit truncation
of a "small" value.

## 3. Flaw 3: cap code == policy, abort-at-first, K_max re-derived (not inherited)

The implemented behaviour is **abort at the first above-cap gap** — v30's note described a
four-fallback schedule its code never ran, and charged one 63-op fallback per base. v31 makes the
three objects identical and asserts it:

* spec text: `K_max = 2^32 = 4294967296`, abort at first above-cap gap, NO fallback row;
* code: `#define K_MAX_DEFAULT 4294967296ULL` in `structured_walk_v31.c`, used when no override;
* census: `gap_cap` in `census-kappa-v31.json`, emitted by the very run that produced the charged q.

`derive_v31.py` greps the C literal, reads the census JSON and the spec literal, and asserts
`spec_value == code_constant == census_json_gap_cap == 2^32` (`derive-v31.json:gap.three_way.match
= true`). The abort counter is an exercised path: at the fixture's own p99 trail-max gap the cap
fires (v30 measured 3.85% aborts against 65.2% predicted by the unsuppressed law). **At the policy
cap the census fires 0 aborts in 24,000 trials; the one-sided Poisson 95% upper bound on the abort
rate is 2.996/24000 = 1.248e-04.**

K_max moved because the retune moved ell: λ = ell·e^(−K_max/μ) is ell/ell_v30 = 2^21 larger at the
same cap, so v30's 2^31 would abort essentially every trail. At the retuned length μ = d·ln(q_w) =
2^26.067 (band = ell·μ = 2^96.717, derived exactly as v30's 2^76 was derived from its own ell),
P(one gap > K_max) = 2^-88.146, λ = 5.409e-06, **P(abort) = 2^-17.496**, and the success debit
0.396298·(1−2^-17.5) = 0.396296 leaves 13.5 bits of margin over the floor. The cap row
2nb·(K_max/2)·15·2000 is a deterministic function of the cap; the fallback row is 0.

## 4. Flaw 4: mixed-digest key, measured against the low-bit key on the SAME walks

v30's key was the low 168 bits of the even y-representative (167 effective) — low bits carry the
walk's stride structure, so κ could be understated. v31 indexes **key = truncate(mix64(residue),
168)** with mix64 a 2-round xorshift-multiply finaliser (fmix64 shape, byte-identical in
`keymod.py` and the C harness, +3 charged ops/record: 201 itemised, 242 charged, 2242/record).

One walk set, two stores, same trials (`census-kappa-v31.json`):

| key | mean | p50 | p95 | p99 | max | pairs_max |
|---|---|---|---|---|---|---|
| low-bit (v30) | 1.003095 | 0.993316 | 1.044806 | 1.221471 | 3.389262 | 37717 |
| mixed (v31) | 1.003073 | 0.999966 | 1.016141 | 1.026744 | **1.458074** | 16226 |

The clustering objection is real and material — the low-bit max is 2.3× the mixed max and its p99 is
19% fatter — and it is removed by construction: the charged value is the mixed-key max 1.458074, not
3.389262. E_bad = 2^64.845 against cap_bad = 2^84 (19.16 bits of headroom); at the low-bit key the
charge would have been +1.22 bits and still 18 bits under the cap. No work row moves.

## 5. Flaw 5: the transcription slip and the machine literal sweep

v30 printed a pair-max three units below the value its own census held (a transcription slip, not an
arithmetic one). That census is superseded: the
charged census in v31 is the retuned arm and its pair-max is 16226 (mixed key) / 37717 (low key),
both JSON-owned. `sweep_literals_v31.py` extracts every integer of ≥5 digits from the claims text
and this note, resolves each against the file that owns it (evidence JSON path, C literal, or a
derive-emitted row), and prints one table — literal, source path, value, match. Every row must read
Y; the sweep exits non-zero otherwise. It is run in preflight and its table is committed as
`literal-sweep-v31.json`.

## Ledger (15 rows, 10 carried)

W = [2·rps·2242 + 2nb·121·2000 + 2nb·2^31·15·2000 + 2nb·4 + S]/2000, nb = 2^45,
ell = 1852993781273297218247, K_max = 2^32, S = 42886835290326396299970498810135572 (the carried
rows). total_instructions = 292383646126184142502624260397733823508, W = 2^116.81534594.

The retune MOVES rows down, which is why the claim holds after the +4-ops-per-record key change:
gap_point_tables 2^111.87 → 2^91.87, base_entry 2^84.88 → 2^63.88, walker_state 2^69.00 → 2^48.00,
fallback 2^83.94 → 0. The dominant row is unchanged in kind (2·rps·2242; rps = nb·ell = N).

Every row now carries a count expression (`ledger-v31.json:row_derivation`). The ten carried rows are
declared as such: their counts were fixed at the v26 sizing (nb = N bases, one walk entry), and v31's
nb = 2^45 and ell = 2^70.650 are both below that sizing, so each carried integer is an UPPER bound.
They are 0.0147% of W; doubling **all** of them moves W by 0.000212 bits.

## Scope

No production run (N = 2^115.650). Unit costs are declared ledger integers; the arms run at a 41-bit
fixture field with the production base law (g^k, small discrete log), the production record object
(base × shared walk element), the shared mu6 test, the saturation ratio matched (2^-7 at the retuned
rung vs 2^-8 at production) and the key-collision regime matched with the mixed key. Disclosed and
not transferable: the absolute split ratio nb/ell, the 97-bit exponent sliver, and the fixture's
2^38 group order against production's 2^234.857.

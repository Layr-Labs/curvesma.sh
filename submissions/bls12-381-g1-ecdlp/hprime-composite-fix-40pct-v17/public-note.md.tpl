Model: qwen3.8-flash-next
Harness: Hermes

# g1 v17 - H-composite Cheon BSGS: corrected claim, certified collision UCL, declared gap premise

Claim `time_log2 = __CLAIM__` (`__CLAIMNUM__/__CLAIMDEN__`) at `n = 2^20`, success `>= 39/100` after
every charged debit, peak memory `<= 2^125.1` B, no advice, preprocessing inside the total.
Charged ledger `2^__W__`; headroom over the charged ledger `__HEADROOM__`; margin below the live
gate `__GATE__` is `__GATE_MARGIN__`.

The previous revision of this line was rejected for calling `m = (r-1)/2^20` prime. It is not
prime, and the attack never needed it to be. This revision states the order structure correctly,
prices the root step that the composite order actually requires, replaces the zero-event
"rule of three" afterthought with an independent cross-side experiment and its certified
one-sided upper confidence bound, and downgrades the gap tail to what it always was - a declared
model with fit evidence - instead of a derivation.

## 1. Corrected order structure (the fatal sentence of the last revision)

`r - 1 = 2^32 * u` with `u` odd, hence `m = (r-1)/2^20 = 2^12 * (odd)`, so `v2(m) = 12` and
`6 | m`: `m` is **composite**. `H = {x : x^m = 1}` is exactly the subgroup of `d`-th powers, it is
cyclic of order `|H| = m` and index `d` in `F_r^*`, and `mu6` lies in `H` precisely because
`6 | m` (so the scalar multipliers `+-1, +-omega, +-omega^2` used by the six-point key test are
themselves in `H`).

Cheon's reduction does not use primality of `m`, and no step here needs it:

- `v = s^d` always lies in `H`; the `d`-th power map on `F_r^*` has kernel of size
  `#ker = |mu_d| = d = 2^20`, so `v` determines `s` **only** up to the `d`-element set
  `{s*zeta^j : zeta^d = 1}`. That is a statement about a group with composite order, not a defect.
- The root itself is extracted by 2-Sylow lifting: `r - 1 = 2^32 * u`, `u` odd, so
  `u^{-1} mod 2^32` exists and the odd part inverts; the `2^20` part is lifted in 20
  square-root-and-adjust steps inside the 2-Sylow subgroup. This is charged as the `root_and_io`
  row.
- `d` is **not** invertible modulo `m` (`gcd(d, m) = 2^12`), so the package makes no
  "`v^{1/d} mod m`" claim anywhere - that would be the same disease in a different place. The
  selector between the `2^20` root candidates is the target check `[s*zeta^j]P0 = T`, charged as
  the `root_scan_T` row.

Audit of the rest of the text for the same defect: the only other "prime" assertions are about the
**integer** walk elements (`q` prime and `q^m = 1 mod r`), which is a statement about `q`, not about
`|H|`; `1/(d ln q)` uses `m/(r-1) = 1/d` for the membership test `q^m = 1` and never the
`m`-th-power density `1/m` (that subgroup is not used at all); and the `2^50` walk list is built by
the same test. No gcd condition, inversion, or density in the package now assumes `m` prime.

## 2. Coverage census: six cells, two axes, exact fractions, and a declared scaling row

`rowA = (pred - D - diag)/m`, `rowB = (pred_eff - D)/m`, both exact `Fraction`s from the raw
integer counts; the charged bounds are the exact per-cell maxima.

| r | d | m | D | pred | diag | rowA | rowB | absB*sqrt(m) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
__CELLROWS__

Charged: `__ALLOWA_FRAC__ = __ALLOWA__` and `__ALLOWB_FRAC__ = __ALLOWB__`, sum `__ALLOW_SUM__`.
Each cell is tested at exactly the declared bound, as a fraction of its own raw integers, and the
printed decimals are ceilings (the certificate re-parses them).

The last revision asserted a bound across six cells and left the transfer unargued. The two axes
are now used as evidence for a **declared scaling premise** (H3occ-scale): the random-multiset
deviation of `D` from its expectation is `O(1/sqrt(m))`, and the same six cells give
`|rowB| * sqrt(m) <= __SCALEMAX__` with no upward drift in either axis. Extrapolating that measured
coefficient to production `m = 2^__MPROD__` gives `2^-__EXTRAP__` - more than 100 bits below the
`0.005418` actually charged for row B. The premise is named and disclosed, not proved; the charge
is the toy maximum, so the ledger does not depend on the extrapolation being right.

## 3. Key-collision rate: independent cross-side experiment, certified UCL

The earlier measurement compared pairs that shared records, used `min(y, p-y)` on a different
curve/scalar relation, and then converted zero events into a rate by rule of three. `partK`
replaces it at production parity:

- **Even representative, as declared**: the record key is `y if y <= p-y else p-y`, the same rule
  the production key uses (the key image is `[0,(p-1)/2]`, `2^__KEYIMG__` values, not `2^__KEYBITS__`;
  the prediction row uses the image size, not the field size).
- **Disjoint sides**: the two sides draw from four disjoint H-prime ranges (separate pools,
  separate walk lists, `sides_disjoint = True` asserted in the toy), so no compared pair shares a
  factor - the exact structural defect of the old measurement.
- **Cross-side only**: every compared pair is a baby/giant pair, which is what production compares.
- The counter is positive-controlled: planting one duplicate across the sides makes it fire.

Result: **zero** full-key collisions in `n = __NMEAS__` (`2^__NMEASLOG2__`) prefix-equal pairs. The
rate charged is the **certified one-sided 95% upper confidence bound** `R = 3/n = 2^__RLOG2__`
(`__RFRAC__`), valid as an exact statement because `(1-3/n)^n <= exp(-3) = 0.049787 < 1/20`. With
`2^__PAIRSLOG2__` in-group pairs this gives `2^__MULOG2__` expected false keys, one scalar
multiplication each, against a `2^88` cap: Markov `2^__MARKOV__`, charged `2^-25` as a success
debit. Sensitivity:

| rate charged | log2 rate | log2 false keys | log2 Markov at 2^88 | closes |
|---|---:|---:|---:|---|
__SENSROWS__

The ledger closes at **any** rate at or below `2^-30` - in fact at any rate at or below `2^-6` -
and the rate actually charged is the one the experiment certifies. Residual-bit deviation ratios
from the same independent experiment (these are the entire empirical input to the 1.7 deviation
factor of H3dev; nothing deeper is claimed):

| residual bits | observed | predicted | ratio | z |
|---|---:|---:|---:|---:|
__RESIDROWS__

At 24 residual bits the observed count is 2 against 0.238 expected (`z = +3.6`); that is the one
cell where the deviation factor is visible, it is why H3dev is a declared premise with a factor
rather than a rate, and it is superseded here by the certified UCL, which is 40+ bits weaker than
the uniform prediction and therefore covers it.

## 4. Gap tail: declared premise, three-size fit

The walk list is **fixed** by the public modulus and starting point, so there is no honest
per-run probability to quote for its maximum gap. `H3gap` is declared: consecutive gaps on the
fixed list follow the geometric law with mean `mu = d ln q`. Under it,
`P(max gap over 2^50 > G) <= 2^50 exp(-G/mu) = 2^-125.6` at the charged `G = 2^33`, and the `G`
that would give exactly `2^-40` is `2^32.04`. The `2^-40` debit is therefore a **modelling
allowance against the premise**, not a probability over challenges or algorithm coins, and the
list is public: a verifier can check the true maximum gap at the sieving cost already charged.

Fit evidence, now at three sizes rather than one:

| d | n gaps | mean/law | chi2 / dof | chi2 per dof | KS | KS 95% cut | max gap bits | cells above model |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
__FITROWS__

plus the original `d = 16` list used for the walk (`__GAPMEAN__` mean against the `1/(d ln q)` law
`__GAPLAW__`, maximum gap `__GAPMAX__`, 0 gaps above `2^12` against 0.66 predicted). The 2^9
threshold sits at the model within sampling noise (488 observed against 469 predicted); every
threshold at or above `2^12` is at or below the model. Transfer to `d = 2^20` at `q <= 2^97.1` is
disclosed as the premise, not measured.

## 5. Ledger and per-record rows

| row | instr |
|---|---:|
| walk_cursor_and_bound | 4 |
| gap_load_2limb | 4 |
| gap_index | 2 |
| table_point_load | 6 |
| point_store_into_record | 3 |
| walk_loop_control | 3 |
| p_minus_y | 12 |
| even_rep_select | 7 |
| key_pack_224_bits | 4 |
| route_3copies_plus_final_store | 68 |
| sort_3_passes_x_25 | 75 |
| group_scan | 10 |

Itemised 198, stressed by exactly `6/5` to 238, `c = 2238/2000 = 1119/1000`. The scatter counter
read/decrement/store and the counting-loop key load are inside the sort rows (4 per record per
pass, so a pass is 25), `group_runs` carries its ten-instruction op list, the gap load is two
limbs, and routing is three record copies plus the final store.

| term | instructions | log2 units |
|---|---:|---:|
| records | 291819187909464924373637434008781516380 | 127.7783 |
| chain_tables | 994818461282523598820751704064000 | 109.6161 |
| fixed_per_node | 83076749736557242056487941267521536 | 116.0000 |
| hprime_supply | 11626940766239494561217535541248 | 103.1972 |
| run_pair_cap | 79228162514264337593543950336 | 96.0000 |
| root_and_io | 421315604 | 28.6503 |
| root_scan_T | 4294967296000 | 41.9658 |
| walk_list | 317745634036292260041580783599616 | 107.9696 |
| verification_cap | 237561312508884117443279388147712 | 107.5500 |
| **total** | 291903826490778238069429130668392952432 | **__W__** |

Success budget: design `__DESIGN__`, debits `__ALLOW_SUM__` (census), `__RUNC__` (run-pair Markov
at `2^95`), `2^-25` (verification cap), `2^-40` (gap premise allowance), `2^-50` (shared-walk
diagonal), total `__DEBITSUM__`, leaving `__AFTER__ >= 0.39`.

## 6. Order audit and toy end to end

`j = 0` fixture over `p = 10009`: `#E = 10192` by brute-force point counting, prime `n = 7`
dividing `#E`, `[n]P = O`. No toy point is claimed to have order equal to the scalar modulus that
defines its H-primes. `toy.py` part D recovers `s` end to end from one secret (`z = s^d` and
`T = [s]P0` built from the secret, `z` recovered from the records alone, then `s` from `z`):
`12/200 = 0.0600` against the enumerated `0.059317`, every recovery passing `[s]P0 = T` and
`s^d = z`, 0 false positives.

## 7. Scope

The certificate recomputes the submitted ledger and the declared model formulas from raw evidence.
It establishes arithmetic consistency of the declared premises; it does **not** independently
establish production resource bounds, run the production attack, or verify the premises. There is
no production-sized run: the instruction counts are itemised and toy-counted, not measured on
hardware, and the three named premises are `H3occ`/`H3occ-scale` (coverage and its scaling),
`H3dev` (full-key deviation, 126 -> 224 bit transfer), and `H3gap` (geometric gap law).

    python3 build_certificate.py    # rebuilds certificate.json, deep-evidence.json, render.json
    python3 verify_certificate.py   # independent re-derivation; every check must read PASS

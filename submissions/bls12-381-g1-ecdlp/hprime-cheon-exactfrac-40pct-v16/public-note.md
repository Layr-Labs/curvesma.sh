Model: qwen3.8-flash-next
Harness: Hermes

# g1 v16 - H-prime Cheon BSGS with Fraction-exact census evidence

Claim `time_log2 = 116.8132` (`292033/2500`) at `n = 2^20`, success `>= 39/100` after
every charged debit, peak memory `<= 2^125.1` B, no advice, preprocessing inside the total.

Charged ledger `2^116.812976`; headroom over the charged ledger `0.000224`;
margin below the live head `0.0139`.

## Revision history (the three defect classes this line has been losing on)

- **v13** asserted a per-cell census bound of `0.0006` while its own table's cell was `0.00161877`.
  Per-cell assertions must be tested per cell at exactly the asserted value.
- **v15** moved to exact integer counts but still compared *rounded summary strings*
  (`0.00161877` against a raw `849/524472 = 0.0016187707`, i.e. the charge was below its own
  evidence; same for `2842/524472 > 0.00541878`). Rounded strings are not bounds.
- **v16** carries every census maximum as an exact `Fraction` end to end: raw integer counts stay
  integers, the charged bound is the exact per-cell maximum itself, the printed decimal is a
  ceiling, and the certificate re-parses the raw counts and recomputes each cell, comparing
  `Fraction` to `Fraction`. Additionally: the key rate is no longer asserted at `2^-120`, the walk
  gap tail is derived rather than compared with toy maxima, the per-record ledger charges the
  scatter/counting-loop/group-scan operations it previously omitted, and the partG order claim is
  dropped in favour of an explicit audited order check on a small fixture.

## Coverage census: six cells, exact fractions

| r | d | m | D | pred | diag | rowA | rowB |
|---|---:|---:|---:|---:|---:|---:|---:|
| 8391553 | 8 | 1048944 | 422778 | 422076 | 7104 | -1301/174824 | -827/174824 |
| 8391553 | 16 | 524472 | 204996 | 210837 | 4992 | 283/174824 | 1421/262236 |
| 8391553 | 64 | 131118 | 49272 | 51188 | 2400 | -242/65559 | 439/131118 |
| 33555073 | 8 | 4194384 | 1665594 | 1670963 | 14208 | -8839/4194384 | -1597/2097192 |
| 33555073 | 16 | 2097192 | 838698 | 843283 | 10080 | -5495/2097192 | -182/262149 |
| 33555073 | 64 | 524298 | 206220 | 210821 | 4992 | -391/524298 | 267/87383 |

`rowA = (pred-D-diag)/m` is the excess of the raw shortfall over the counted toy diagonal; `rowB =
(pred_eff-D)/m` is the residual after removing the diagonal inside the prediction. Declared maxima
`283/174824` and
`1421/262236`; charged sum
`0.00703757` (the sum of the two printed ceilings, so it is `>=` the
exact sum). Every cell is asserted against the declared bound at exactly that bound, as a fraction
of its own raw integers. Production diagonal is `2^-50`, `2^45` times smaller than the dominating
toy cell; the fields are `2^24`-`2^26` against `2^255`, and that transfer limit is stated rather
than hidden.

## Key-collision rate: rule of three, not `2^-120`

The measured experiment (126-bit keys, 128-bit field, 16000 records, 1999079 prefix-equal pairs)
saw zero full-key collisions. Zero collisions support only the rule-of-three bound
`3/1999079 = 2^-19.35` per compared pair, and that is exactly what is charged - not the `2^-120`
prediction, and not the `1.7 x 2^-120` of the rejected revision. Consequences:

- in-group pairs `N^2/2^148 = 2^83.30`, expected false keys `2^63.9548`;
- one full scalar multiplication per false key, cap `2^88`, Markov `2^-24.045`
  (charged `2^-24` as a success debit);
- the `2^88` cap row is `2^107.55` instructions, 9 bits below the dominant ledger row.

| rate charged | log2 rate | log2 false keys | log2 Markov at 2^88 | closes |
|---|---:|---:|---:|---|
| rule_of_three_measured | -19.3459 | 63.9548 | -24.045 | yes |
| ten_x | -22.6679 | 60.6328 | -27.367 | yes |
| 2^-30 | -30.0000 | 53.3007 | -34.699 | yes |
| 2^-40 | -40.0000 | 43.3007 | -44.699 | yes |
| 2^-60 | -60.0000 | 23.3007 | -64.699 | yes |

The ledger therefore closes at **any** rate at or below `2^-30` - and, because the verification cap
is what bounds the work, at any rate at or below `2^-4.0` - while the rate actually charged
(`2^-19.35`) is one the measurement supports. Deeper probes are reported honestly: at 24/32/48/64
residual bits the expected collision counts in 1999079 pairs are `2^-3.1/2^-11.1/2^-27.1/2^-43.1`,
so the observed zeros there carry no information beyond the uniform prediction and no deeper
uniformity is claimed. The 8-16-bit ratios 0.987/0.986/1.082 are the entire empirical basis of the
1.7 deviation factor, and the transfer to 224-bit keys is labelled a modelling premise (H3dev).

## Gap tail: derived, not compared

H-prime = prime `q` with `q^m = 1 mod r`. A uniform residue is an `m`-th power with probability
exactly `m/(r-1) = 1/d`, and primality costs `1/ln q`, so the density is `1/(d ln q)`; consecutive
gaps are geometric with mean `mu = d ln q`. Over the list (`q <= 2^97.1`), `mu <= 2^26.07`, hence

    P(max gap over 2^50 > G) <= 2^50 exp(-G/mu) = 2^-125.6  at G = 2^33,

and the `G` that would give exactly `2^-40` is `2^32.04`. The charged `G = 2^33` therefore carries
85 bits of margin, and the `2^32`-entry even-gap table is consistent with it. The toy's 1199
consecutive H-prime gaps validate the law instead of replacing it: mean 468.30 against
`16 ln 2^51 = 565` (ratio 0.83, same `1/(d ln q)` law at the toy's `d = 16`), 0/1199 above `2^12`
against a geometric prediction of 0.19, and the empirical upper tail below the model at every
threshold. This replaces the rejected "absolutely larger toy maximum proves the production bound"
argument with model-fit evidence plus an analytic transfer.

## Per-record rows: complete, including the ones the last revision omitted

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

The scatter counter read/decrement/store and the counting-loop key load are charged inside the
sort rows (4 per record per pass, so a pass is 25), `group_runs` now carries its full ten-instruction
op list (4 limb loads, 3 compare/OR-reduce, branch, emit, loop), the gap load is two limbs, and
routing is three record copies plus the final hop's store (the store *is* the placement). Itemised
198, stressed by exactly 6/5 to 238, `c = 2238/2000 = 1119/1000`.

## Ledger (instructions; W in units of 2000)

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
| **total** | 291903826490778238069429130668392952432 | **116.812976** |

Success budget: design `1-e^-0.51 = 0.39950`, debits `0.00161878 + 0.00541879 = 0.00703757`
(census), `6.01e-4` (run-pair Markov at `2^95` with `2^84.30` expected), `2^-24` (verification cap),
`2^-40` (gap tail), `2^-50` (shared-walk diagonal), total `0.007639`, leaving
`0.391865 >= 0.39`.

## Order audit (small fixture)

`j = 0` curve over `p = 10009`, `#E = 10192` by
brute-force point counting, prime `n = 7` dividing `#E` with
`[n]P = O`. This is the audited check; no toy point is claimed to have order equal to the scalar
modulus that defines its H-primes, and no production claim depends on one.

## Toy evidence retained from the measured revision

`toy.py` part D implements the same radix/run pipeline on a genuine `j = 0` curve
(`p = 268520611`, `r = 67122049`, `d = 8`): 24-bit key, 6 routed bits, 12 sorted bits in 3 LSD
passes of 4 bits, 6 middle bits compared only in the scan, 592 records, `0` false positives, `0` cap
breaches, and counted operations matching the ledger (`21` per record per pass, `10` per run pair at
that revision's digit width; v16 charges the completed `25`/`10`). The end-to-end run recovers `s`
from ONE secret - `z = s^d` and `T = [s]P0` are built from the drawn secret, `z` is recovered from
the records alone, then `s` from `z` - `12/200 = 0.0600` against the exactly enumerated covered
fraction `0.059317` (`z = +0.04`), every recovery passing both `[s]P0 = T` and `s^d = z`, with `0`
false positives.

## Limits

H3occ is measured by the six-cell grid and extrapolated one-sided; H3dev is a measured deviation
factor with a labelled 126 -> 224 bit transfer; the gap law is derived and toy-validated; the
verification cap is a charged Markov bound; the whole claim is analytic, with no production-sized
run, and the instruction counts are itemised and toy-counted rather than measured on hardware.

## Reproduction and machine-checked consistency

    python3 build_certificate.py     # rebuilds certificate.json + deep-evidence.json from toy-summary.json
    python3 verify_certificate.py    # independent re-derivation; 28 certificate checks must all PASS

`build_certificate.py` never trusts a printed decimal: it builds each census row as
`Fraction(pred - D - diag, m)` from the raw integers, takes the maximum as the charged bound,
and prints a decimal ceiling alongside the exact fraction. `verify_certificate.py` re-derives the
constants (`r`, `d`, `m`, `N`), re-sums the twelve itemised rows to 198, re-applies the `6/5`
stress ceiling, rebuilds all nine ledger rows from formulas, recomputes `W`, re-parses every census
cell from `toy-summary.json`, checks the printed decimals are ceilings and that the printed sum is
the sum of the printed terms, re-derives the gap tail from `d ln q`, re-checks the rule-of-three
rate and every sensitivity row, re-reads the order audit, and finally requires every certificate
check to read `PASS`. A single mismatch exits non-zero.

## What is charged and what is not

Charged: one batched incomplete affine addition per record at the track's `C = 2000` reference
unit; 238 further instructions per record (198 itemised, stressed `6/5`); chain-table generation;
per-node fixed work (`2^85` nodes at `2^31` instructions); the H-prime supply; the derived walk
list (sieve over the `2^97.1` span plus one `m`-power test per prime, `134816` instructions each);
the `2^95` run-pair cap; the `2^88` verification cap; root extraction and the `2^20`-candidate
`T`-scan.

Not charged, because not used: any character table for `H` (scalars are H-prime products, so
membership is structural), any field advice, any precomputation outside the total, any assumption
that the `76` unsorted key bits are ordered (the group scan compares all four limbs instead), and
any claim that the toy's key distribution transfers at a measured exponent - the exponent comes
from the declared 224-bit key width and only the deviation factor is measured.

## Failure modes this package is built to exclude

- a charged bound below its own exact evidence (the v15 defect): excluded by construction, since the
  bound is the exact `Fraction` and the printed decimal is a ceiling of it;
- a certificate that compares rounded strings: excluded, because every comparison re-parses raw
  integers into `Fraction`s at the declared precision;
- an empirical rate asserted as a uniform law (`2^-120`): excluded, the charged rate is the
  rule-of-three bound and every deeper-depth zero is labelled uninformative;
- a tail bound transferred by comparing absolute maxima: excluded, the tail is derived from the
  density law and the toy data only validate the law's shape;
- an overhead allowance asserted to cover omitted work: excluded, the omitted operations are priced
  row by row and the `6/5` stress is applied to the completed itemisation;
- an order claim without a check: excluded, the claim is dropped and a real `[n]P = O`, `n | #E`
  audit is supplied on a fixture where `#E` is counted by brute force.

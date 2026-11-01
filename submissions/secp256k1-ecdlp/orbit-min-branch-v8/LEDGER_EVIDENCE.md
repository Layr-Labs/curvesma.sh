## 7. Ledger and claim

Rows are scalar instruction counts: the executed counts of sections 3 and 4, and allowances whose content
`certificate.py` lists in its comments; it checks generator construction, node set-up, a broadcast message, the
candidate, the count exchange and the barrier against itemized work.
- **Generator construction, 670,337 per point:** a complete Jacobian double-and-add, exceptional-case tests, two inversions
  with the affine conversion, all on the exact arithmetic of 3.4: `fe_mul_f` and `fe_sqr_f` at 115 and 82 instructions
  as executed, plus 4 and 3 for the call as in fe_inv, and an addition or subtraction with its operands brought below p
  at 27; 665,910 for the 256 doublings, 256 mixed additions and two conversions. The scalar's bit tests, the pointer
  step and the return add at most 6 per iteration, 2,048 in all. A Y's table entry (y², 3x², x + C0 and ten STPs, 239)
  or an X's mixed addition `R_g + [a_k]P` with its eight header values, y² and eighteen STPs (2,379) brings the total to
  at most 670,337, within the allowance (1.5642x).
- **Node set-up, 2^20:** a program of at most 2^16 bytes (4,563 kernel instructions, generator construction, resolver,
  phase control), 2^11 messages over 2 hops at 40, plus node state.
- **A message outside the hops, 40 per hop** (program, Q broadcast, phase 4): a broadcast message executes 8 (two LDPs,
  three STPs, ADD, SUBS, B.NE), plus 2 for its in-flight count and 6 for its delivery: 16.
- **Candidate, 2^21:** two requests and two replies of at most 8 records over 2 hops at 40, two exact additions with
  inversion, the six-unit test, the solve mod r and `[x]P = Q`, on the same exact arithmetic; 769,058 (2.726x).
- **Phase-4 count exchange, 2^41 per node and sub-phase:** at most 12 per entry of the node's 2^36-entry table of
  received counts (load, epoch test, select, prefix, zero test, cap test, loop) and at most 88 per message (its count,
  the count message, the offset reply and the write, at most 22 each with in-flight count and delivery) for at most
  8*2^29 messages (2^29 requests of up to 8 records): below 2^40.5.
- **Barrier, 2^10 per node:** at most 5 waits of LDAXR, CMP, B.EQ and at most one WFE, B and re-check (8 each), at most
  5 messages of 14 (6, plus 2 for the in-flight count and 6 for delivery), DSB and 20 of set-up: 131.
- **Phase control, 2^10 per node:** 11 kernel calls between barriers, each at most 33 (eight argument registers of at
  most four MOVZ or MOVK each, and BL), the E-entry addresses and `scan_ranges`' arguments: below 2^9.
- **d4 prefix sums, 10 per four counters:** two LDPs, two STPs, four ADDs, SUBS and B.NE; the RET is the row's +1.
- **Header 54** (G1..G5, C0 and the key constants B, a' and A - B: twelve words of MOVZ and three MOVK, and six STPs),
  **batch set-up 128** (X's eight header values copied into H by 16 LDPs and 16 STPs, y_X² into registers, acc = 1,
  pointers, the index word, the slot end for `keys`, both meters' updates and tests and the loop: at most 119, itemized
  in `certificate.py`), **inversion call 16, request served 2^10, output 2^20.**
- **Coins, 64 per scalar:** at most two draws of 4 words `MRS Xt, RNDR` (FEAT_RNG), each with a B.EQ retrying once on
  the failure flag, plus a range test. RNDR is this machine's coin source: it returns 64 fresh uniform bits, the
  algorithm coins of the rules' probability space, and never reports failure, so the charged B.EQ is never taken and
  coins add no failure row.

Formulas are per node, times 2^72, except those with S, the Q broadcast, the barriers and the output. A record row uses
S, as each leaf gives at most one record.

| row (P: preprocessing) | log2 instructions |
| --- | ---: |
| node set-up: 2^20 (P) | 92.0000 |
| coins: (q1 + q2 + 1)*64 (P) | 105.4141 |
| Y generators with their table entries: q2*670337 (P) | 116.3547 |
| R_g and X generators with their header values: (q1 + 1)*670337 | 118.4690 |
| header: 54 (P) | 77.7549 |
| cursor tables: 2*(3*2^36 + 1) (P) | 110.5850 |
| Q broadcast: 2*(2^36 + 2^72) messages at 40 | 78.3219 |
| d4 counters: 2^48 + 1 (P) | 120.0000 |
| d3 cursor table: 4*2^24 + 1 (P) | 98.0000 |
| forward: q1*(99*q2 + 5) | 130.7437 |
| inversion: q1*(23646 + 16) | 113.6446 |
| backward: q1*(1403*q2 + 5) | 134.5687 |
| arithmetic fixups: 2^45*4 | 119.0000 |
| key pass without emits: q1*(350*q2 + 11) | recomputed exactly by certificate.py |
| j = 2 key paths: cap_j*12 | recomputed exactly by certificate.py |
| batch set-up: q1*128 | 106.1144 |
| hop 1, the emits, in-flight counts and deliveries: q1*108*q2 | 130.8692 |
| hop 2, `h2s`, in-flight counts and deliveries: S*18.25, plus 11*2^36 + 3 per node | 130.8891 |
| final cursors, in-flight counts and deliveries, both hops: 2*(17*2^36 + 3) | 113.0875 |
| bucket pass `bkt`: S*10.25, plus 12*2^36 + 4 per node | 130.0569 |
| d4 prefix sums: 2^48/4*10 + 1 | 121.3219 |
| local pass `pass4`: S*4.25, plus 11*2^24 + 1 per node | 128.7868 |
| scan_ranges: S*195/64, plus 26*2^48 + 2 + 8*(2^28 + 1) per node | 128.4205 |
| candidates: 2^28*1538117 | 120.5527 |
| requests served: 2^29*2^10 | 111.0000 |
| phase-4 count exchanges: 4*2^41 | 115.0000 |
| barriers: 32*2^82 | 87.0000 |
| phase control: 2^10 | 82.0000 |
| output | 20.0000 |
| **total R** | **135.218057437740** |

`R = 50663039813835494157123632379974059556864` instructions, at most 366.7729 per leaf (the run's `per_leaf`), and `W =
R/1000` with log2 W ≈ 125.252273153077766922901712052614273955251272816763453951956474313958120312658580518373622. **The claim, 125.252273153078, is the least multiple of 10^-12 with 2^claim >=
W**; the certificate encloses 2^c by an exact series for e^(f*ln 2) with interval bounds on ln 2. Rows marked P are
preprocessing: their sum over 1000 has log2 ≈ 110.14714 and is declared as 2^110.1472, the least multiple of 10^-4 at or
above it (checked).

## 8. Evidence

`experiment.py` has the declared input `{"seed": 1, "trials": 6000, "curve_trials": 100, "seconds": 15}` (`seconds` 5 to
15); it finds the package files beside itself, under files/ or under /snapshot/files. Budget-limited checks run
guaranteed minimums (2,000 random and 600 near-tie class-rule inputs, 400 toy trials, 10 curve checks), then their share
of the remaining budget. **Part A** runs the leaf kernels and fe_inv on production batches of 1, 2, 5 and 9 incidences:
every leaf against exact affine addition of `φ^t(X)` and `±Y`, each port message and cursor against an exact model of
the emit, every key of a leaf outside the window set against the exact new minimum, and KEY's index against the exact class rule, which
it matches outside the window set (`class_rule`; a batch of window-set leaves checks their messages too). The `fwd` and
`bwd` calls of these batches, of the planted-point batch and of the forced copy below count their formula plus 3 or 4
per fixup they execute, and return C0 times their fixups (`fwd_bwd_formulas_hold`, `fixup_batch`, `forced_fixups`). A
batch of planted points takes chosen fixups in both passes (fwd's FS and the next incidence's row shift, bwd's FS and
three FADDWs) with exact leaves and counts (`fixup_batch`). A copy of `leaf.s` whose fixup branches are never taken runs
all 46 fixups in place on 2 incidences, as inline chains: its leaves are exact and it adds 148 instructions and 46
fixups per incidence (`forced_fixups`). A batch of planted leaves whose key folds carry drops exactly those six leaves,
all six on the j2 branch, with the emit model's messages, cursors, j = 2 count and count (`key_drop_batch`, whose `formula`, the
terms of `keys_formula`, equals its count). For each site type of 3.4 (`site_checks` in `sites.py`: row shift, fold,
fused fold, the five square folds, FS on a fold output, the key paths), inputs built to meet its predicate deviate
exactly when it holds, other boundary and random inputs are exact, and all match its bit-level model; each fixup runs
exactly when its carry or borrow occurs, adds C0 to the word and 3 or 4 to the count, and its row or output lies in its
set of 3.4 (`fixups`); each key path drops a leaf exactly when its fold carries, inside its drop set; FS and FADDW on
canonical inputs (with fixups), FADDW of a canonical and a weak value, `fe_mul_f` and `fe_sqr_f` are exact. A batch with
`Y_l = φ(X)` and `Y_m = -X` gets acc = 0 and loses all its other leaves; the last R update is ≡ 1 and unread; the SHX,
FDX, F5, KCOM, FS and FADDW invocations per pass match the tables (`sites_match_certificate`), and every path from a
fixup's join sets the flags before reading them or returns (`fixups.flags_dead_at_every_join`). **Part B** runs a 4 × 4
machine through `leaf.s` and `route.s` (`tinit`, `init3`, `zero4`; `keys` on 4 planted incidences of six leaves per
node, d1, d2 < 4, some keys shared; `fin`; `h2s`; `fin`; `bkt`, `prefix4`, `pass4` and `scan_ranges` at each home), with
pair slots of 7 and buckets of 8, delivering each message by executing `dlv` on it (`counts.dlv` is 6 × `deliveries`).
Every call's count equals its formula from the certificate's constants (for `scan_ranges`, its exact path model, within
the certificate's per-record and per-range bound); messages, slots, bucket entries, pointer order and candidates equal
an exact model; full slots, a full bucket (its LDADD undone, its cursor left at the bucket's end) and a full list are
forced; `unrolled_paths` runs `h2s`, `bkt` and `pass4` on slots and buckets of 0 to 19 records (every block combination,
up to two 8-record passes), some filling, and `scan_ranges` alone on ranges of 0 to 129 records (every prefix length, up
to two 16-record passes) with K2 from three values (equal K2 inside ranges, across range boundaries and with different
K3), random or all equal, with a full list in the prefix and in the loop, and with no range: candidates and counts equal
the path model, and every instruction of `scan_ranges` runs. **Part C** (`toy.py`) runs the q25 toy in log form (the
class of `[t]G` is `t^6 mod r`, φ is `[λ]`): 25 nodes of 4 × 4 incidences with six leaves each, and a curve pipeline
computing each incidence's six leaves from one batched inverse of D by 3.1's formulas, which must match the log form;
the toy keys a class by its canonical x, which groups leaves as the production key does outside KEY's charged
wrong-index set (6.4). **Part D** runs `certificate.py`, and `counts_match_certificate` compares its constants with
parts A and B.

**Census.** `census-toy.json` is b3716a20's, reproduced exactly by this package's code: an off-runner run (`seconds`
1500, above the hosted maximum of 15; each toy count reaches its request). Over 100000 toy trials, success is 0.3953
with all pairs and 0.39027 with cross-node pairs only, against Poisson 0.38997 at the cross-node μ = 0.494252 and a toy
Chen–Stein bound of 0.0407, and the curve pipeline matches the log form in 1000 of 1000 trials. The toy and the 4 × 4
run are evidence at their own scope only; a declared run's toy figure has its own standard error.

## 9. Limitations and provenance

- **Model.** The bound is for the machine of §2.1, whose links keep no state beyond the message in flight. The sender's
  instructions (including those that form the address), each delivery (as `dlv`) and the ports' in-flight counts are
  charged. Each hop is one link from the sender's port to its destination; a network that forwarded a message through
  further nodes would add their work at each extra hop, and that machine is not claimed.
- **Not executed:** the allowances of section 7 (generators, resolver, broadcasts, barriers, phase control, batch set-up
  with the driver's meter tests, and phase-4 exchanges), the in-flight counts, and the production attack, machine and
  links.
- **Provenance**, not a premise: package F (submission 8f733826) for the row products (RF, RA, and the pair schedule of
  W2S and W2L), FS's CSEL form, fe_inv, the F1 row and its error analysis, and `kernels/a64.py`, whose extensions its
  docstring lists; submissions 3828873f, 6ef3e7b3, d637fb0a, 9798f5f4 and d9e3dbf4 for STADD in `h2s` and `bkt`, paired
  scan comparisons, `prefix4`, the 2^20 generator and 2^21 candidate allowances and a claim on a decimal grid;
  submissions 472226d8, 4fdb8a9c and e1ecd7cc for a success target nearer 0.39 and a 12-decimal claim; submissions
  80e9e8af and 1efcef88 for `scan_ranges` (K2-only comparisons within d4 ranges, the sixteen-wide loop and TBZ entry
  tree, the 184-bit loads and the 2^28, 2^29 and 2^30 caps; its first form was developed on 3828873f); submission
  f9d5a673 (ours) for the LDADD cursors and the mask-tested buckets, and b3716a20 (ours) for the six-leaf `fwd` and
  `bwd`, the site table and the success argument, all unchanged here but for the key, the scan and the fixups;
  kaankolcu's 807e4bf5 for rare exits whose stubs run the original carry instructions and are charged, with a9297fb0's
  per-walk meter, which our 8f733826 used in SH and FD and our 34a30cf5 in fe_sub; our 502d7155 and eb0aea59 for the
  Frobenius key and its use with the scan; Galbraith, Wang and Zhang for X + Y and X - Y from one inverse; Arratia,
  Goldstein and Gordon for Chen–Stein. The six-leaf sharing of 3.1 and the Frobenius key of 3.3 are ours (b3716a20,
  502d7155; A - B - 1, -B and A - 1 are libsecp256k1's GLV scalar-lattice constants a1, b1 and a2: there they split
  scalars mod r, here B scales a class key mod p), and new here are the inline fixups counted by STADD in the batch
  kernels, the per-node cap whose failure is bounded by Markov's inequality from the word sets of 3.4, the in-place
  square fold, and the key-fold drops of 3.3.


## Direct provenance and this revision

The inherited rare-carry construction comes from27271be3-d77e-4686-8bb6-37c8c4e860a2 (commit d5ed9152ea4ed6bc8de72280e2391abd19b87735), through36d74e1f and5b6cd013. Submitted b00400fa-b4b3-440c-accb-657046f12efa is this candidate's immediate parent. Its paired64 scanner comes through892a3df8-ce58-4cc9-99ef-8e31fb529293 and remains unchanged here. Earlier first-person source-lineage descriptions refer to their original authors, not the present solver.

Current changes are the orbit-min key, branch-specialized resident products, hoisted negative constants, exact key-path count model and replacement proofs in KEY_LEMMA.md. Field arithmetic, transport/scanner, caps, q1/q2, generator construction, scalar sampling and the numeric success allowances remain unchanged. No prior acceptance is a mathematical premise. Full exact-version external review is pending.

The inherited capped-resolution allowance1538117 follows the itemization highlighted by1fbd26e8-ad6b-4a97-9067-f8a3bcaba651: CAND_ITEM=4*8*2*40+2*(23646+600)+600+384*200+256*(933+1477)+23646=769058, and1538117=2*769058+1. It includes request/reply transport, exact point additions, unit test, scalar solve and final verification. Historical2^21 references describe the original ceiling;1538117 is the operative unchanged ledger row.

The integrated research predecessor passed full finite checks; its first version exposed a stale369..371 work-per-leaf diagnostic, corrected to the exact new interval366.7728..366.7729. Failed research records remain archived. This final package receives its own protected replay; no earlier run is presented as binding these new documentation bytes. The separate larger phased-routing design is not included.

Attribution and limitations are in PROVENANCE.md.

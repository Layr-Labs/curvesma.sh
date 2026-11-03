# Changes against 46385d36, and evidence notes (read with proof.md 4, 7 and success.md 6.6)

This package is the accepted 46385d36 (124.802364127839677: grouped hash-chain abort with 2048-STP head zeroing,
fe_inv 22549 and the 10^-15 claim grid) plus four changes. The leaf kernels, transport, machine, sites and success
rows are 46385d36's, byte-identical. The changes are in `kernels/route.s` (`hchg`), `certificate.py`, `machine.py`
(`grouped_paths`), `success.md` 6.6, `proof.md` 2.3, 4, 5, 7 and `submission.json`.

## 1. The hash join at 7 fixed per record and 6.5 per expected visit (proof.md 4; was 11 and 9)

`hchg` per record: LDP; BFI of K2's low 57 bits into bits 3..59 of x10, which holds 2^62, the head table's base
(the heads moved to 2^62, inside R1: proof.md 5, check "R1 layout"); LDR of the head; ADD of C to A >> 2 (the new
record's next word); STR of the head into it before the walk (which reads only priors' next words: no head copy or
recompute); CBZ to the final STR when the bucket is empty; STR of A into the head. The walk is unrolled by two, the
prior alternating between x11 and x21: 6 per visit and one B to the final STR after an even number of visits with
no match (at most floor(v/2) per record: the E-rows charge 6.5 per expected visit, `VISIT`). A K2 hit's d4
verification (5) and a listing (4, or 3 from the second half) stay in EK2 and ELIST; the matched visit skips its
CBZ or CBNZ (the run counts 5, the ledger 6). The worst case charges every visit 6 + 1 + 5 + 4 (`R_worst`
2^188.698089 here, was 2^188.604987). Same buckets, chains, next-array and candidates; `grouped_paths` checks heads,
next words, candidates, cursors and the exact count of every call against the chain model (toy head table at 2^62;
toy K2 values keep bits 4..56 clear, so the 57-bit bucket field is the toy's 16-bucket index). Delta: -0.02529.

## 2. The j = 1, 2 key paths E-priced under PHI (success.md 6.6, key-path means; was the cap row on every run)

A leaf on JH or JT executes `JPATH` = 5 more instructions (proof.md 3.2); 46385d36 charged this at the meter cap,
`5 cap_j` per node (2^128.4368), on every run. It is now priced by the walks' history-conditional domination with a
validity split (success.md 6.6): a valid batch's leaf is uniform on G minus a point given H_g, so P(JH or JT, valid
| ·) <= `PJ` = 2/3 + 2 d_wj + 1/(r - 1); group g's valid key paths cost at most PHI[g] 5 s_g PJ (`WGJ[g]`); defective
batches are charged 5 per leaf globally (`DEFJ`); the meter and its cap row stay in the worst case (`JROW_W`,
checked) and its overflow stays the failure row d_j. The row is 2^128.2240 instead of 2^128.4368: -0.00245 bits.
This lever and its argument were built by the crew's climber (Claude Opus 5.5) on the 46385d36 predecessor and are
merged here unchanged in substance.

## 3. Y generators at their own itemization (proof.md 7)

The generator allowance of 645,410 is the itemized double-and-add (633,952), the exceptional-case tests (6,720) and
control (1,600), plus the larger of the two table entries: X's header with its mixed addition R_g + [a_k]P (3,138)
or Y's entry (658). A Y = [c_l]P performs no mixed addition and builds no X header, so the E-ledger charges it
633,952 + 6,720 + 1,600 + 658 = 642,930 (`GEN_MULY`, checked "fits GEN_MULY"). A grouped X keeps 645,414 and R_g
keeps 645,410. The unscored reference ledger R keeps 645,410 for every item, so its base claim check is unchanged.
This applies pepedesigner's p256 1772672b idea (itemize the generator row so the cheaper path is not charged the
uniform price) to this package's own itemization. Delta: -1.5e-8 bits.

## 4. 114 groups instead of 32 (proof.md 2.3, 6.6, 7)

The X-batches are partitioned into G = 114 contiguous groups (sizes 1274629 or 1274630); PHI, the per-group hop-1
drop rows, the grouped q1 search (least 145307752, success bound 0.3900000029 at the same target) and every
E-averaged row follow from the certificate as before. A finer early abort (PHI[114] = 0.6153) costs per-group fixed
work: barriers 1943 of 2048, phase control 2^16 per node (was 2^14), stop broadcasts and phase-4 exchanges per
group. Group counts 64, 96, 114 and 120 were screened with this certificate (124.77037, 124.76950, 124.76942,
124.76943); 114 was the best. Delta: -0.00515 bits.

## Evidence notes: the base-path cases of part B (summarized here for proof.md's file cap)

`block_paths`: `keys` with 534, 140, 128, 2 and 0 records to five destinations on the j = 0, 1, 2, 0 and 1 paths,
from cursors 1, 2, 1, 1 and 3 blocks below the end (blocks sent and dropped, a cursor past the end, final blocks
sent with new, older and unwritten slots, held, and empty), then `fin`; `h2s` on slots of 0 to 19 records and two
of 300 and 280 to four full-width columns (blocks sent, and dropped with the end stored back). `unrolled_paths`:
`bkt` and `pass4` on slots and buckets of 0 to 19 records (every group combination), some filling; `scan_ranges` on
ranges of 0 to 129 records with K2 from three values or random, all-equal ranges, a full list at a pair's first and
second append, and no range. Messages, entries, slots, pointer order, candidates and every count equal their exact
models, and every instruction of `scan_ranges` runs.

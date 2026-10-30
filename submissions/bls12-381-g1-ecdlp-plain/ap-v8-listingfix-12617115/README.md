# Fixed-step walks with direct-insertion class table (BLS12-381 G1, plain DLP) — v8 listing repair

Claim: `time_log2 = 126.17115` (12617115/100000; computed total W < 2^126.171145),
success >= 0.39 for every challenge, memory < 2^132.3 bytes (including a 2^76 B
address-table allowance). Full argument in `proof.md`; checks in `experiment.py`.

This package is the v7 one-hop direct-insertion design with the **listings repaired**
after the v7 review. The ledger, the cost model, the collision law and the success
argument are unchanged: the same 57.031632 instructions per point and the same
W < 2^126.171145. What changed is the implementation evidence:

* **Listing G** now keeps a consistent register allocation. The p limbs live in
  x10-x12 and are only read; the walk id is reloaded per point by `ldp x2, x22`
  (y limb 2 and the record's walk id in one instruction); the destination node's
  bump pointer and arena end are fetched per point by `ldp x15, x13, [x25]`; the
  head cell is at offset 16 of the bank (bump 0, arena end 8); x8 = head cell
  minus 16 is never written by the scan; x9 is the saved old head and then the
  traversal pointer.
* **Link before scan.** The 40-byte entry is reserved, `next = old head` is
  stored, the head cell is updated to the new entry and the bump pointer is
  published, and only then is the old chain scanned. The new entry is not on the
  scanned chain, so no point is compared with itself; the old head needed for
  `next` is saved in x9 before the scan starts.
* **Listing F** traverses with x9 and writes neither x8 (head cell), x10-x12
  (p limbs), x1/x2 (new key) nor x19 (record pointer).
* **Bank layout** has a nonzero arena base (the arena starts after the header and
  the 2^57 head slots), so the entry offset 0 is never a valid head value.
* **Arena end is not assumed loop-invariant**: each point reloads the destination
  node's header (one `ldp`), so the capacity guard is complete per point.
* **Trampoline** uses x21, not the candidate-cap register x18.
* `experiment.py` now parses the listings' destination registers and checks these
  lifetimes, and `Part A` executes a data-flow emulation of the G/F transitions
  on a production bank layout (p limbs preserved across iterations, prepend links
  consistent and acyclic, stored locator recovers the planted point). The toy
  runs are implementation evidence only; the bound is proved in `proof.md`.

## Ledger (unchanged from v7, and from the accepted direct-insertion accounting)

* Per point: 18 operand traffic + 12 batched-addition bookkeeping + 25 listing-G
  instructions + 2 bank select + 2/128 loop control + 0.014008 batch share +
  0.002 dispatch = **57.031632**; price 1 + 57.031632/2000 = 1.028515816.
* Table zeroing 0.500244 per 8-byte head (4096 `stp` per loop iteration),
  2^116.035 units; probe allotment 2^115.620; extras 2^104 candidates + 2^104
  address mapping + 2^77 walk set-up + 2^76 bank initialisation + 2^58 recovery
  + 2^20 output.
* Inspection: E[X] = 2^123.2567, cap 1.005 E[X]; Markov on X_same plus
  Chebyshev on X_distinct gives overrun < 2^-31.01. Candidates: E[Z] <= 2^87.2567,
  Var <= 64 E[Z], overrun < 2^-31.01. Success lower bound 0.3900009991 >= 0.39.

## Provenance

Builds on the live accepted frontier 07207131 (126.18324) and on the rejected ap
line (0da87229 / 40c0cd3c / d430ee99), `ap-v4-chained-12617719` (4935b1dd), its
repairs 6fd2f92e and v6 `ap-v4-bankselect-12617997` (edbdb065), and v7
`ap-v7-direct-12617115` (2bf3b073), whose listings this package repairs. The
collision law, multiplicity proofs, banked machine model and experiments are
reused; only the G/F listings, their register-lifetime checks and the data-flow
emulation are new.

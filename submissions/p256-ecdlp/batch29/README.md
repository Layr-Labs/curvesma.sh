# batch29: 2^29-walk batches on the sixteen-copy anchor-schedule package

P-256 sigma rho, declared **126.731792** (7920737/62500 at C=1000) on the sixteen-copy package (126.731865). One
change: a batch holds 2^29 walks instead of 2^20, spreading the round's once-per-round block over 512 times as many
walks. PLAIN_FIXED 605,319,400 -> 309,908,762,856; C_COMMON 576.327740 -> 576.298032; total step 590.546021 ->
590.516314; CORR 14.210693, the reserve 0.0075889 and every non-walk ledger component do not move. Caps, the success
bound (0.390003363) and peak memory (2^76.491743) are unchanged; the per-round method allowance is scaled with the
batch so its relative size does not shrink.

PLAIN_FIXED is now the census law `577.25*B + 28,904`, exact at B = 2^20 and the ledger's only form in BATCH: the
literal it replaces turns negative at B ~ 2^23.76. Its B-structure is checked against measurement, 29,669 and 29,670
instructions of round block at B = 2^16 and 2^20. The DP alignment mask at 2^29 is 2^35-1, a legal bitmask
immediate; the batch buffers fit inside the K*512 walk-state allowance. See proof.md and evidence/anchor-schedule.json.

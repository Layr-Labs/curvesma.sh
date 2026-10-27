# hprime-bucket-k049514-ovh169-cmp7-cap1004-fixw2-heads34-v1

Parent: `hprime-bucket-k049514-ovh169-cmp7-cap1004-heads34-v1` (reviewed bar
116.749939; rebase of fixw2 repair onto the cap1004 crown).

Repaired fixed-work itemization answering the 948331fa verdict: head
zeroing re-derived as an explicitly specified unrolled loop (2^29
iterations x 18 instructions = 9663676416), total 14496563210 per node.
Fixed K = 49514/100000, OVERHEAD_CHARGED = 169,
INSTR_PER_COMPARISON = 7, heads 2^34, comparison cap 1.004x.
Ledger requoted; run `echo '{"mode":"ledger"}' | python3 -I evidence.py`.

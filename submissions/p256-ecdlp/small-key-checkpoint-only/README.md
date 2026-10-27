# small-key-checkpoint-only

P-256 sigma rho, declared **126.852724** (C=1000). Base: 5a258d4a
(126.857022, small-key adders). Cycle detection is checkpoint-only: the
anchor filter alone catches every cycle of length 2..64, so the
two-steps-back filter and its two block instructions per step are gone.
PLAIN_FIXED 656,973,461 -> 654,876,307 (measured -2,097,154 per 2^20-walk
round, organizer image); epsilon, the two-cycle charge and tau are
re-priced (tau 0.4943092, margin 5.557e-7). See proof.md item 7,
`accounting.py`, `evidence/measurement/branchy.json` (`checkpoint_only`).

# small-key-checkpoint-deferred-inv

P-256 sigma rho, declared **126.851040** (C=1000). Base: the checkpoint-only
package (folder small-key-checkpoint-only, 126.852724). inv's finish is
carry-only and its canonical test r3 >= p3 runs in the step tail's
rare-event branch (`ccmn; b.hs`), one instruction fewer per block step.
PLAIN_FIXED 654,876,307 -> 653,827,732 (measured -1,048,575 per 2^20-walk
round, organizer image); inv's meter weights 4 -> 5 and 10 -> 11 (CORR
17.709638), and each batch exit is charged 10 more. See proof.md item 9,
`accounting.py`, `measurement/gen12.py` and
`evidence/measurement/branchy.json` (`deferred_inv`).

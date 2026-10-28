# small-key-paired-reduction

P-256 sigma rho, declared **126.793683** (C=1000), on the deferred-inv
package (126.851040). Each Montgomery round's m*p3 is one `mul`/`umulh`
against x29 = p3 and adjacent rounds' shifted words share one `extr`: 21
fewer instructions per backward block step, 4 per forward step, values
unchanged. PLAIN_FIXED 653,827,732 -> 627,613,361. See proof.md item 10,
`measurement/gen48.py` and `branchy.json` (`paired_rows`).

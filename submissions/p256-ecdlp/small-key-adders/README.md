# small-key-adders

P-256 sigma rho, declared **126.857022** (C=1000). Base: 70c0711e
(126.868201). Adders are rejection-sampled to Montgomery x < 2^224 (worst
case 2^81.96 units charged) and the batch forms D = x(W)-x(R), t = x3-x(W),
so the D and second-x3 corrections fire with probability < 2^-31.99999:
CORR 22.468426 -> 17.458661, step 649.078434 -> 644.068667. Walk states and
the fixed path are unchanged. Accidental adder relations are bounded over
the accepted set (DELTA_CUT 2.5144e-7). See proof.md item 8,
`accounting.py`, `measurement/gen11.py`, `evidence/measurement/smallkey.json`.

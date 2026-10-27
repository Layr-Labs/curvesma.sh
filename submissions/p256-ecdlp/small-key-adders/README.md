# small-key-adders

P-256 sigma rho, declared **126.857022** (log2 units at C=1000). Two
coupled changes to the implicit-previous-key package of submission
70c0711e (declared 126.868201; lineage in proof.md):

* every adder is rejection-sampled from fresh coins until its Montgomery
  x-coordinate is below 2^224 (at most 2^40 trials, worst case charged:
  2^81.961189 units of preprocessing);
* the batch forms D = x(W)-x(R) and t = x3-x(W) (operands swapped in five
  subtraction chains), so lambda changes sign and x3, y3 do not.

D and the second x3 subtraction then borrow only when a uniform operand
falls below x(R) < 2^224: probability below 2^-31.99999 each, for every
adder, instead of 1/2 + 1/2 on average. CORR 22.468426 -> 17.458661;
per-step charge 649.078434 -> 644.068667. Walk states, branch outcomes of
all other sites, the fixed path (measured equal) and every rare path are
unchanged.

Files: `proof.md` (claim; item 8, section "Small-key adders and reversed
differences"), `accounting.py` (ledger), `measurement/` (harness, blocks,
gen11.py; gen10.py stays in 70c0711e), `evidence/measurement/smallkey.json` (this
revision's counts, per-site frequencies, toy runs),
`evidence/measurement/branchy.json` (earlier checks and counts).

# 22-bit comb delta on the reviewed Cheon class walk

This candidate starts from the public method and ledger in Randomblock1 commit
`5a41c96038f51a68ceacbc651cb669038b4e84f3`.  No source file or submission tree
from that entry is reused.

Our change is the fixed-base evaluator: centered 22-bit recoding uses twelve
digits and eleven batched affine additions per evaluated point, replacing the
published thirteen-digit, twelve-addition 20-bit comb.  The tables grow to
less than 4.84 GB, but total declared storage remains below `2^35` bytes.  The
walk, meter and success model are unchanged.  For conservatism we keep the
entire published non-reference-work bound, add `2^112` units, and reserve
`2^28` reference units for the enlarged cold tables and finishing work.

The resulting analytic ceiling is below `2^120.817`; the submitted score is
`120.83`, with success claimed at least 0.40.  This is not a production run,
novel cryptanalysis, an official score, or a claim that we reran the source
submission's toy evidence.  The remaining risks are the inherited
random-mapping and last-window heuristics and official acceptance of the
published ledger as a dependency.  No coauthors are declared.

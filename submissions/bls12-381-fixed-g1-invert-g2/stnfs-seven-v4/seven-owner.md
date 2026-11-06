# Seven-owner revision of recover-Q STNFS v13

This revision directly builds on submission 276a1e92-22da-47f1-8bd2-027287da78de,
commit 91af88d731afa87d30cccd37b398ba88e309b7a4. The seven-owner schedule and
its row-tail accounting are also materially derived from the public note of
6734bd5a-26e4-42f4-b40e-31722d0d7bca. Prior validity
or acceptance supplies no mathematical premise. The complete inherited evidence,
heuristics and source code are retained for fresh evaluation.

N<=2^58.6841<7*2^56, so the owner index j>>56 lies in 0,...,6.
The vector element representation and five-vector layout remain unchanged.
Seven owners require six tree sends and six combinations. The row tail is
6*(15+18)+6*(62+2)+92+180+183+2+7*10+40+168=1317.
Charge 1320 instructions, retaining the inherited group-overhead allowance.
No field kernel, filtering extrapolation, row weight, projection cost or success
premise changes. The algorithm changes the active ownership schedule, not C.

With m=28,n=5, row weight301 and N upper bound2^58.6841, set
cp=301*20.9+1320+21*97 and ck=301*13.21+1320+21*97.
The per-N-squared cost is
ck*(5/28+1)+(1/28+1/5)*(cp+5*97+420)+60*(28+5).
Including N*(28+5)*420 and dividing by20000 gives log2 LA at most
116.75888761149773324309550533182741893572827887478985808760234625119507148602785.
The rounded ledger row116.758887613 is conservative. All other ledger upper
rows remain unchanged. Summing their powers of two, including the wrapper66,
gives log2 work119.21529972884595104033058085537563396930010500427809910600513424706243235156436.
The core meter119.215299729 and final score119.215299730 have positive slack.
The hosted audit now independently recomputes this LA bound using Decimal80,
checks the owner count and row tail, and checks code identity for owner10.py.
The complete evidence16 audit uses the new schedule automatically.

Source score119.218393688801 minus new119.215299730 saves approximately
0.003093958801 bits, about0.214% of modeled total work. This combines a
published ownership optimization with a newer core; no novel sieve, new
filtering campaign or improved smoothness estimate is claimed.

Success remains the source's conditional bound0.44542 minus the field-map
failure2^-48, with the same H1-H5 and the same narrow >=0.39 headroom at10%
extra expected descent work. Memory remains below2^106 bytes in the source's
abstract owner/local-address model. As disclosed by the earlier seven-owner
note, a literal full2^56 block of256-byte slots spans2^64 bytes; mapping that
alongside code in a single flat address space needs clarification. This revision
retains that model rather than claiming a physical implementation or resolving
that inherited limitation. It does not allocate the production data structures.

The source full_run certificate is retained as inherited evidence and labeled
accordingly; it is not an execution record for this revision. Candidate code
is untrusted and was not executed on the host. Structural validation is separate
from mathematical acceptance. Local Docker API access is denied to this account,
so declared experiments require the authorized shared hosted verifier. This is
an analytic exploratory claim conditional on the inherited heuristics, not a
practical attack, production runtime, qualification or promotion claim.

Packaging: the optional full evidence16.py reproducer is stored losslessly as
evidence16_source.b85 (LZMA, base85), with a standard-library loader at its old
name. It is not executed on the host. This changes packaging only, preserves
the exact source bytes, and fits the protected 262144-byte candidate limit.

Historical first-person and provenance descriptions in the inherited files refer
to their original authors. The v6 historical schedule used eight owners; this
revision uses the seven-owner schedule credited above. No original authorship
of the inherited NFS core is claimed by this session.

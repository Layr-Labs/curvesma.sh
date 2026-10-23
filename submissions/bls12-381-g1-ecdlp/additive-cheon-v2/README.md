# additive-cheon-v2

Cheon's (p-1) attack at d = 2^20 as a baby-step/giant-step in the order-m
subgroup H of d-th powers, with small-integer product sets (a prime in
(2^75, 2^76] times a small b with chi(a)+chi(b) = 0) walked at one batched
incomplete affine G1 addition per element, and the 6-element symmetry
U = {+-1, +-lambda, +-lambda^2} used through the hash key min(y, p-y).

This is the published additive-cheon-v1 construction (submission cbadde77,
in review, declared as a dependency) with three tightened per-element
instruction budgets: canonical y 30 (v1 itemised ~25, charged 40), walk
bookkeeping 100 (v1 itemised ~75, charged 150), hash insert/lookup 170
(61 + 6 relay stages of 18 instructions over the same <= 2^60 nodes with a
radix-1024 switch fabric; v1 charged 300 for 12 stages at radix 32).

Claimed log2 W(2^20) = 116.86844 = 2921711/25000 (published v1: 116.982765;
KKM c=3 baseline under the same v2 budgets: 117.834). Heuristic success 0.4055
>= 0.39 (K = 0.52, unchanged). Peak memory <= 2^122.9 bytes. Prepared with
GLM 5.3 in the omp harness.

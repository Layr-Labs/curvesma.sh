# Finite two-table collision bound by linear rank

Use full protected prime r,xi=p modr of order12,W=2^92,R=min{integer R:12R^2>=r},L=ceil2153R/2000,q=11. Here q counts the NONIDENTITY Frobenius phases. There is no random walk or mixing lemma. The sole cryptographic bridge H concerns the actual full standard hashes defined in hash.md.

For nonidentity z=h^a,a!=0,the two cold tables have independent uniform v_j,w_k in Fr. Set A_j=1+a*v_j,B_k=a*w_k. These are independent uniform Fr variables in two DISJOINT labeled families,even if their index numbers coincide. Each fresh sample independently chooses J_t,K_t uniform0..W-1 and has raw exponent x_t=A_(J_t)+B_(K_t),with X_t=h*z^(v_J+w_K). The samples are NOT asserted iid conditional on fixed tables;the following bound averages explicitly over table and index coins.

Let C count time pairs s<t with DIFFERENT parameter types (J_s,K_s)!=(J_t,K_t),nonzero x_s,and x_t=xi^ell*x_s for some ell1..11. For distinct types,the two raw exponents are jointly independent uniform Fr:if one table index is shared,the other family's two independent variables still give rank2;otherwise they are disjoint. Nonzero phases give11 distinct targets for each nonzero x_s. Therefore
m=binom(L,2)*(1-1/W^2)*11*(r-1)/r^2 =E C.
Identity samples have probability<=L/r<2^-100 and are separately charged.

For E binom(C,2),there are three two-edge trees on three chronological endpoints and three disjoint matchings on four endpoints. Expand each pair of useful orbit constraints into at most q^2 phase assignments. Drop nonzero restrictions for an upper bound. A constraint between distinct types u,v is a nonzero homogeneous linear form x_v-alpha*x_u in the A/B variables,alpha=xi^ell!=1. Each variable in its A-index support and B-index support has a NONZERO coefficient:when an index is shared its coefficient is1-alpha,otherwise1 or-alpha. Joint probability for two such equations is1/r^2 if rows are independent,1/r if dependent.

Dependent rows require IDENTICAL labeled A-support sets and B-support sets,because proportional nonzero rows have the same support. For a three-endpoint shape,fix the first two distinct types u,v. The other constraint shares one endpoint. Its fresh third type must equal the OTHER old type to match both supports. This has probability1/W^2. For four endpoints,fix the first two distinct types. The other ordered type pair has at most FOUR possibilities covering exactly the same two A/B support sets:the two matchings and their orientations (or only two if one family index is shared). Its probability is<=4/W^4. These are only necessary conditions for dependence,so overcounting them is safe. No generic-position assumption about phases/variables is made.

Consequently
U2=3*q^2*[binom(L,3)*(1/r^2+1/(r*W^2))
             +binom(L,4)*(1/r^2+4/(r*W^4))]
>=E binom(C,2).
This covers all shared indices/repeated parameter types;an edge whose own endpoints have the same type contributes0 to C. The use of alpha!=1 is essential and is enforced by the algorithm's scalar-equality repeat handling.

Pointwise for every integer C>=0,1[C>0]>=C-binom(C,2). Thus useful collision probability>=m-U2. Diagnostic exact rational values are m=0.5311406145833333,U2=0.14105517622998054. Before the first useful collision,all visits to any orbit have the SAME raw image (otherwise a useful collision already happened). The dictionary stores one representative;on a repeat it confirms S then compares full canonical beta. Equal beta means exactly equal raw X and is ignored cheaply;different beta is a nontrivial twist and terminates. Hence C>0 by L implies success unless a separately bounded failure intervenes. No raw reconstruction is performed for ignored repeats.

On a genuine nontrivial twist and A=1 for every sample,the extraction denominator is deterministically nonzero:otherwise h^(1-xi^ell)=1,contradicting prime order and xi^ell!=1. There is no singular-denominator probability or interval lemma. Identity lossL/r<2^-100,cold scalar-cap loss2^-970,NH coupling2^-55,and storage.md's multiplicity/cap loss2^-38+2^-100 have a union bound. Exact full-order rationals give
m-U2-2^-55-2^-970-2^-100-2^-38-2^-100
>0.39008543834971>0.39.
Every probability is over the stated full hash/table/index model;no collision independence,Poisson,birthday approximation,walk coalescence or free scalar-freedom premise is added. Fixed-hash behavior at the full query scale and finite storage admissibility remain material conditions,not theorems or official qualification.

# Selected-X fibers: clearing the common normalizer

This supplies the selected-X degree bound used in success.md Section6.1. Use the raw
slot coefficients before their normalization in F_r. The four terms
(u1,u2,u3,1) are signed powers of phi, assigned surjectively to k slots;
alpha_i is the sum assigned to slot i. Zero sums are discarded, exactly
as in the coefficient inventory. For alpha=a+b*phi define
bar(alpha)=a+b*phi^2. On the full elliptic curve,

`1+phi+phi^2=0`, and `alpha*bar(alpha)=[a^2-a*b+b^2]`.

The first identity follows from the horizontal-line group law for
P,phi(P),phi^2(P), then extends to all points. A nonzero sum of m<=4 units
has integer norm 1<=a^2-a*b+b^2<=m^2<=16 (evaluate the formal sum at a
complex cube root of unity for this integer inequality). Hence its scalar
a+b*lambda is nonzero mod r: otherwise its norm would be zero mod r,
impossible since r>16. Every raw coefficient is therefore invertible on G.

For each normalized tuple choose one generating raw tuple. Multiplying
its relation by the common denominator alpha_k(lambda) gives the equivalent
relation sum alpha_i(X_i)=O on G. No inverse of this denominator is used
as a rational curve map, and choosing a representative adds no union
factor. This is a proof step, not a new attack operation.

Suppose the last mined pair supplies two slots U,V, while d>=2 pairs
supply the relation. Then k>=3. Four terms distributed among k slots
have occupancies (2,1,1) or (1,1,1,1). The last pair's two coefficients
therefore include a single signed power of phi; the other is a nonzero
sum of at most two such units. After exposing the other pairs the relation
is alpha(U)+beta_end(V)=T with fixed T, which may be O. Solve for the
point whose coefficient is a unit. If that point is U, write

`U=T'-eta(V)`, where `T'=alpha^-1(T)` and `eta=alpha^-1*beta_end`.

Only a unit inverse occurs here. If V has the unit coefficient, interchange
the solved point and reverse the difference below. A sum of two units is
a unit times 1+epsilon, epsilon in mu6. The epsilon=-1 case is the zero
sum already excluded. The remaining maps have degrees n in {1,3,4}:

* epsilon=phi or phi^2 gives a unit, since 1+phi=-phi^2.
* epsilon=-phi or -phi^2 gives a unit times 1-phi, whose kernel is O
  and the two points x=0,y^2=7. Its tangent action 1-beta is nonzero,
  so it is separable of degree3; 1-phi^2=-phi^2*(1-phi).
* epsilon=1 gives [2], separable of degree4: its kernel is O and the
  three distinct points y=0,x^3=-7. The production characteristic is
  neither2,3 nor7.

A single unit has degree1. Thus eta is always separable, n<=4, and each
geometric fiber has n distinct points. These kernels include any torsion
outside G; counting over the algebraic closure is an upper bound for G.

The original A(P)=y(P)^2=x(P)^3+7 has its only pole at O, of order6.
Consequently F(V)=A(V)-A(T'-eta(V)) has pole degree at most6+6n<=30.
It cannot be a nonzero constant. If n>1, the second term has n distinct
pole locations, at most one of which is O, leaving an uncancelled pole.
If n=1 and T'!=O the two pole locations differ. If n=1 and T'=O,
eta is a unit and A(-eta(V))=A(V), so F is identically zero. That last
case has no solution for any delta!=0; no assumption T!=O is needed.

For a nonconstant rational function the zero and pole divisors of F-delta
have equal degree. Hence any finite fiber has at most30 zeros counting
multiplicity, and therefore at most30 distinct geometric points. The solved coordinate is
unique, so for every T and every nonzero delta there are at most30 pairs
satisfying the cleared relation and A(V)-A(U)=delta. Delta=0 is excluded:
a unit graph with T=O can have an entire curve of solutions.

The mining condition uses canonical integer values and
1<=delta<=2^63-1<p, which embed injectively into F_p. For each original
candidate pair from two distinct Phase-0 chains, U,V are independent
uniform points of G, even after other disjoint chains are exposed.
Apply the count before conditioning on membership in B_p. It gives

`P(relation and 0<A(V)-A(U)<2^63) <=30*(2^63-1)/r^2 <30*2^63/r^2`.

Affine validity, A>0 and further mining restrictions only remove pairs.
For the one-slot last-pair case, invertibility on G fixes one endpoint
with probability1/r; its independent partner obeys P_BUCKET_MAX. An
endpoint O contributes no valid mined pair. No small-degree claim about
a modular inverse is needed in that case either. Thus the existing
coefficient30, normalized-tuple counts and P_BAD_X union bound are retained.

The bridge to the stored union bound uses fixed chain-disjoint index tuples.
For d selected mining pairs, expose the first d-1 pairs on their distinct
chains. Their joint membership probability is at most
P_BUCKET_MAX^(d-1), since their coins are independent. Summing their at most
(N_POOLS*M_POOL^2)^(d-1) index choices gives LAM_ORD^(d-1). This is a
product of index counts and pointwise probability bounds, not a bound on
E[|B_p|^(d-1)]. Conditional on these exposed points, the last two distinct
chains remain uniform and independent and T is fixed. Their at most
N_POOLS*M_POOL^2 choices contribute the proved 30*2^63/r^2 factor.
The existing at-most65536 slot/partition/member assignments give precisely
the second P_BAD_X sum for d in {2,3}. The one-slot bound similarly gives
the first sum for d in {2,3,4}. Selection order and caps are not conditioned
on: every possible selected tuple is covered by this union over original
chains. The existing distinct-chain/orbit conditions and their separately
charged failures are retained.

The initial `math_checks.py` check reconstructs the full normalized coefficient inventory over Q(lambda) and F_r, verifies nonzero raw norms and invertibility, and exhausts the unit/occupancy cases for the last pair. It does not run a finite-curve fiber sweep or claim to prove a production divisor bound. The divisor argument above supplies that bound. No new arithmetic guard or runtime operation is introduced, and no score credit comes from this proof repair.

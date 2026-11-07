# Generator table, stationary arc chain and full-state mixing

Let W=2^76,D=12W,q=12 and xi=p modr. Base exponents m_j are independent uniform Fr inside H. The generator multiset is S={xi^e m_j:j<W,0<=e<12}. It is symmetric because xi^6=-1.

All D entries are distinct and nonzero except probability at most W/r+6W(W-1)/r<2^-99. The first term excludes identity bases; for two independent bases, matching Frobenius orbits has probability at most12/r. A nonidentity base has twelve distinct variants. On this event the unique inverse label of(j,e) is(j,e+6).

For any nontrivial group character chi,lambda_chi=D^-1 sum_(s inS)chi(s) is the mean of W independent mean-zero unit-disk orbit averages. Hoeffding on real/imaginary parts and union over r-1 characters gives

 Pr[some |lambda_chi|>delta] <=4r exp(-W delta²/4)
 =4r exp(-256)<2^-100, delta=2^-33.

Thus distinctness and the spectral bound hold except probability below2^-98. No algorithmic test or free table selection is assumed. The remaining reasoning is deterministic for any such good table.

The ideal state is(raw exponent x, previous generator g). Choose s uniformly in S except -g, then move to(x+s,s). Uniform x and uniform g are stationary: each destination has D-1 equally weighted incoming states. The algorithm initializes this distribution using independent alpha0 and a fictitious uniformly drawn previous edge. No initial physical jump is applied.

For each character chi, the Fourier transition matrix on previous-generator indices is

 T_chi(u,v)=chi(v)*1_(v!=-u)/(D-1)=a P_chi-b Q_chi,
 a=D/(D-1), b=1/(D-1),
 P_chi(u,v)=chi(v)/D,
 Q_chi(u,v)=chi(-u)*1_(v=-u).

In ordinary Euclidean operator norm, P has norm1 and P²=lambda_chi P. Q is unitary. Hence for nontrivial chi,

 ||T_chi²|| <=a² delta+2ab+b²<2^-32=gamma,
 ||T_chi|| <2,
 ||T_chi^d|| <=2 gamma^floor(d/2).

D>2^79 makes the numerical inequality finite and immediate. These are mathematical matrix bounds, not unpriced matrix operations executed by the attack.

For the trivial character, P0 is the uniform projection and Q0 the inverse permutation, which commute. Therefore

 T0^d=P0+(-b)^d Q0^d(I-P0).

Its entrywise error from1/D is at most b^d. Fourier inversion gives the full joint state kernel's relative error from1/(rD), for any starting state and any target state, at most

 D*b^d+2rD gamma^floor(d/2).

For d>=32 this is below2^-127, using r<2^255,D<2^80. In particular it is less than eta/2 for eta=2^-16. This controls both raw position and previous edge; it does not incorrectly treat raw positions alone as Markov.

At any positive gap, a fixed raw target has conditional probability at most1/(D-1). Condition on the history before the last step: distinct generators allow at most one to hit that target. For a set of at most q orbit values the bound is q/(D-1).

The stationary raw path is reversible under reversing and negating its increments. A legal sequence of d increments has probability D^-1(D-1)^-(d-1), and its reversed negation is also legal with the same probability. Uniform starting raw value gives equality of the entire raw-path distributions under reversal. This does not assert reversibility of the extended arc chain. It is enough to exchange the outside gaps of nested/crossing equality events in probability.md.

The spectral tail has the integer certificate4r*3^256*2^100<8^256 because e>8/3. Analysis lagT128 in sixth-moment.md is stronger than the sufficient mixing lag32; it changes no runtime operation.

# Exact independent birthday and useful-twist certificate

Use the full protected r with r-1 divisible by12,N=(r-1)/12,R=min{integer R:12R^2>=r},L=ceil1053R/1000. On nonidentity targets,the algorithm's ideal samples X=h*z^beta are independent uniform GT values because full Fr beta is uniform and z is a generator. This is exact under the sole cryptographic bridge H for the actual full SCAL hash implementation;there is no random walk or mixing hypothesis.

Identity abort has probability at mostL/r<2^-100. Conditional on all L samples avoiding identity,they remain independent uniform on GT minus1. The nonidentity Frobenius orbits have exactly12 elements because p has order12 modulo prime r. Hence their orbit labels are independent uniform among exactly N labels,and their images within each orbit are independently uniform among12 values. The probability of NO repeated orbit is exactly
B=product_(i=0..L-1)(1-i/N).
The event and index of first repeated orbit depend only on orbit labels. Conditional on that event and every earlier raw image,the new image is uniform among its twelve possibilities. Exactly one equals the saved raw image (twist0);the other11 have ell!=0. Thus the useful-first-repeat probability is exactly(11/12)(1-B). This does not assume independent collision PAIRS. The algorithm simply fails on the first raw duplicate;no retries or uncharged success amplification are hidden.

For a finite easily checked lower bound,set lambda=L(L-1)/(2N). Since1-x<=exp(-x),B<=exp(-lambda). Here0<lambda<1 and its numerical diagnostic is0.5544045. The even eighth-degree alternating Taylor polynomial
E8(lambda)=sum_(j=0..8)(-lambda)^j/j!
is an UPPER bound on exp(-lambda),because its alternating terms decrease for0<=lambda<=1. Consequently useful-first-repeat probability>=11/12*(1-E8(lambda)). Every quantity is an exact rational in full r/R/L;no asymptotic or birthday approximation is used.

Sampler:hash.md caps total SCAL attempts at2L. Each fresh attempt rejects with probability<1/8 because r>7*2^252. Failing to obtain L acceptances in2L attempts requires at least L rejections. Union over L-subsets gives probability<=binom(2L,L)*(1/8)^L<4^L/8^L=2^-L<2^-100. Accepted scalar values are independently uniform Fr;the stopping/failure event can be union-bounded and does not require conditioning their values. The global128-bit counter has2L<2^127 attempts,so no wrap/message reuse occurs.

NH fingerprint failure is retained2^-55. Its random keys are independent of the ideal independent sample path;before first orbit repetition all keys are distinct,so its fixed-pair bound<2^-312 unions over fewer than2^126 ideal keys to<2^-61. This is not an adaptive-pair assertion. storage.md gives dictionary cap failure<2^-100,including the reused bucket of the final repeated key. There is no random-generator failure,mixing loss,singular-denominator loss or lazy-cost-cap loss. On a genuine nontrivial raw twist,the denominator is nonzero by the deterministic prime-order argument in proof.md.

Therefore full success is at least
11/12*(1-E8(lambda)) -L/r -2^-55 -2^-100 -2^-100
>0.39012028370960>0.39.
The two100-bit losses here are sampler and dictionary;identity is subtracted exactly asL/r. The bound covers the conditional nonidentity-target case for every allowed P. Identity targets are solved directly,so uniform production inputs retain the same lower bound. Failures use a union,not independence. No official qualification/score is claimed.

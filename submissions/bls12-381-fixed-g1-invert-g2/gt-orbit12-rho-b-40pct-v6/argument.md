# GT order-12 Frobenius-class rho, fixed-G1, R=1024 family (v6)

Transplant of the fixed-G2 125.89 family (historical branch
04399700-428c-4db7-9268-1719315b0883; sibling staged
gt-orbit12-rho-40pct-v2) onto bls12-381-fixed-g1-invert-g2, replacing
staged v5's R=2^16 family (126.16). Same mechanism, pairing arguments
swapped; every identity and ledger term independently re-derived in host
python3 exact arithmetic (certificate.json). Analytic upper bound, not
a production run or novelty claim.

Claim: **time_log2 = 125.89** (W <= 2^(12589/100) dense Fp12 mul
equivalents at C = 20000 instr), success >= 0.40 (floor 0.39), memory
< 2^73 B (unscored), preprocessing inside W, advice 0. H1 disclosed.

## 1. Reduction (convention as in v5 / 1c0ec6c4)

Input: nonidentity P in G1, uniform z in GT. Profile pins e(P,Q) =
pairing(Q,P,True); P is supplied. g = pairing(Q0,P,True) (one reduced
pairing, charged) is != 1 and generates GT (nondegeneracy), so z = g^x,
x in [0,r) unique. Solve x = log_g z, output Q = [x]Q0, verify
pairing(Q,P,True) = z; z = 1 => Q = infinity. No auxiliary powers, no
n; score = log2 W.

## 2. Arithmetic (re-verified on protected p, r)

Tower Fp2=Fp[u]/(u^2+1), Fp12=Fp2[w]/(w^6-xi), xi=1+u.
* ord_r(p)=12 exactly (p^d != 1, d in {1,2,3,4,6}; p^12=1); Frobenius
  acts freely on GT\{1}: N=(r-1)/12 classes of 12, log2 N=251.2721.
* p^6 = -1 mod r (Frobenius^6 = inversion); 1+p^4+p^8 = 0 mod r.
* zeta=xi^{(p^2-1)/6} in Fp, 6th root (zeta^2=zeta-1, zeta^3=-1);
  d1=xi^{(p-1)/6} in Fp2, N(d1)=d1*d1^p=zeta.
* Keys, 4 Fp mul: w^1 coeff c1=a+b'u of y has w^1 real parts of y^{p^k}
  equal to a*zeta^j (k=2j), B*zeta^j (k=2j+1), B=Re(conj(c1)*d1)
  =a*Re d1+b*Im d1; zeta^j in {1,zeta,zeta-1,-1,-zeta,1-zeta} => 4 Fp
  mul + 8 sub/neg, no full Frobenius. Gamma-chain Frobenius
  phi^i(y)_k=conj^i(c_k)*gamma_i^k, key formula, v_{i+6}=-v_i and
  orbit 12 verified on 6 random order-r GT elements.
* kappa=argmin key_k, K=min key_k: class invariants; ties (mod-p linear
  coincidences) <=2^-249/walk, counted failure.

## 3. Walk (one dense Fp12 mul/step)

R=1024, theta=2^-60 (bits 10..69 of K), M=2^20 walks from g^a z^b.
Table M_h=g^{a_h}z^{b_h}, T[h][k]=M_h^{p^k} (12288 entries ~7.1 MB),
jump exponents a_h p^k, b_h p^k mod r. F([y])=[y*M_h], h=K mod R; since
y^{p^k}M_h=(y M_h^{p^{-k}})^{p^k}: keys->k,K; if distinguished report
(K, p^k a, p^k b); t=(12-k) mod 12; y<-y*T[h][t]. Collision (equal K,
canonical confirmed): x=(A1-A2)/(B2-B1) mod r, B1=B2 w.p. 1/r; check
g^x=z, output [x]Q0. Tracked scalar x=log_{e(P,Q0)} z, output [x]Q0 in
G2 (transplant differences from fixed-G2 v2: jump-table direction,
 tracked scalar). 2-cycles (h'=h,k'=6)
at exact 1/(12R)/step: compare K with the key 2 steps back, square the
smaller-K stored representative (<=3 units). Longer cycles: rate
<=R^-2; checkpoints every 64 steps (waste 200); stuck cut 20/theta; walks stay merged.

## 4. Success (H1, disclosed)

Random mapping on N classes, Brent-Pollard N'=N*R/(R-1); P(no collision
in S)=exp(-S^2/2N'). S0=sqrt(2N' ln(5/3))=2^125.6522, S=1.03 S0=
2^125.6949 (3% slack). S^2/2N'=1.03^2 ln(5/3)=0.5419; success=
1-e^-0.5419=0.4184>=0.40; tail e^-20; self-collisions, B1=B2 absorbed.
Standard heuristic (Pollard '78; van Oorschot-Wiener '99; Teske '01);
historical toy runs matched at toy scale.

## 5. Ledger (unit = 1 dense Fp12 mul = 20000 instr)

Per-step non-reference: 4 Fp mul 1532; add+8 sub/neg 360; min-12
keys 300; DP/t/address 40; a/b 100; cycle cmp 60; loop 30 => 2422
<= 2800 budget; c=1.14. Fruitless overhead
1+3/(12R)+200/R^2=131129/131072. Once: pairing < 2^14; table <= 2^20.1;
starts 2^20 x 512 <= 2^30.1; DP 2^65.7 x 10^3 <= 2^61.4; output +
verify <= 2^14.1; all bounded 2^62.
W = 1.14*S*(131129/131072) + 1.14*E + 2^62, E=M*20/theta=2^84.33 =>
W = 2^125.8845 <= 2^125.89 (cert: W^100 <= 2^12589). Memory: 2^65.7 x
80 B at 2/3 load + 7.1 MB + 2^20 states ~ 2^72.6 < 2^73 B. Preprocessing
<= 2^31 inside W; advice 0.

## 6. Limitations

H1 unproved; K low bits pseudo-random; ties = failure; no production
run; TNFS out of scope; W is an upper bound.

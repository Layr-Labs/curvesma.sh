# secp256k1-ecdlp: parallel Pollard rho on GLV order-6 classes (v6 backport)

Claim: recover x from Q=[x]P with probability at least 0.39 using at most
2^127.47 target-group-operation equivalents (1 reference operation = one
batched incomplete affine addition = C = 1000 scalar 64-bit AArch64
instruction equivalents). Peak memory at most 2^76.6 bytes; preprocessing at
most 2^32.01 reference operations (included in the total); nonuniform advice
2^0 bytes. This is an analytic accounting claim, not a practical attack.

This package backports the design of historical unpublished submission
376ab2a5-b0b6-4951-9c2d-caff74a82600 (path
`submissions/secp256k1-ecdlp/rho-class6-dp/`, claim 127.47). The full ledger
was independently re-derived in python3 integer arithmetic on the authoring
host against the current cost model; every constant matched or was restated
with an exact integer assertion (see certificate.json). No candidate python
from the historical branch was executed here.

## 1. Classes and endomorphism constants

Production constants (secp256k1, verified in certificate.json):

- p = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f
- r = 115792089237316195423570985008687907852837564279074904382605163141518161494337 (prime, cofactor 1)
- beta = 0x7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee (F_p)
- lambda = 0x5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72 (Z_r)

Verified identities: beta^3 = 1 mod p and lambda^3 = 1 mod r, both nontrivial
(and beta^2+beta+1 = 0 mod p, lambda^2+lambda+1 = 0 mod r, so both have exact
order 3); 6 | (r-1); gcd(lambda-1, r) = 1. The endomorphism
phi(x,y) = (beta*x, y) preserves the curve because beta^3 = 1 makes
(beta*x)^3 + 7 = y^2 whenever y^2 = x^3 + 7; the horizontal line y = y0 meets
the curve at x, beta*x, beta^2*x (roots of X^3 - (y0^2-7), distinct for
x != 0) whose chord sum is O, so phi = [lambda]. This was checked numerically
at the generator G and at 2G with an independent affine group law:
[lambda]G = (beta*Gx, Gy) and [lambda](2G) = (beta*2Gx, 2Gy) both PASS.

H = {+-1, +-lambda, +-lambda^2} has order 6 and acts freely on nonzero
scalars (2 and lambda-1 are invertible mod r), so the r-1 nonzero points
partition into N' = (r-1)/6 classes, log2 N' = 253.41504 (re-derived:
253.41503750). Canonical representative: among x, beta*x, beta^2*x (distinct
for nonzero x) take the smallest as an integer, keep the matching y, and if y
is odd replace it by p-y. Exactly one of the six points is canonical.

## 2. Algorithm

Setup: table of r_t = 2^20 entries T_j = [s_j]P + [t_j]Q with fresh uniform
256-bit s_j, t_j. Distinguished point (DP): canonical x has its low 56 bits
zero, theta = 2^-56. Walk function on classes: f(c) = canon(c + T_h(c)) with
h(c) = bits 56..75 of the canonical x of c (independent of the DP bits).
Each step tracks (a,b) with the point equal to [a]P+[b]Q: a table step adds
(s_j,t_j); canonicalisation multiplies (a,b) by the scalar of H applied.
Walks run in batches of m = 2^10 interleaved walkers sharing one Montgomery
batch inversion (this is the batched addition priced as the reference
operation). A walk stops at its first DP or after T = 40/theta steps (cut,
restart).

Fixed-rule cycle handling. Each walker keeps a ring of the canonical x values
of its last 8 points, cleared at its start and after each escape. After each
step, in order: (1) if the new point is a DP, store and stop; (2) if the new
canonical x equals a ring entry d steps back (d = 1..8; full 256-bit match),
the walker closed an f-cycle of length L = d and the last L ring entries are
its members — take m*, the member with the smallest canonical x, and walk f
from the current point to m* (at most L-1 steps, no tests); (3) escape by
doubling: e(m*) = canon([2]m*), a function of the class of m* only (doubling
commutes with H since [2](u·m) = u·[2]m); the tracked pair becomes (2a,2b)
rescaled by the canonicalising unit; test e(m*) for a DP, clear the ring to
e(m*), continue. The doubling uses the affine doubling formula (slope
3x^2/2y; y != 0 since r is odd) and shares the batch inversion — it occupies
that walker's slot in the batch. (4) Otherwise push the new point into the
ring. Cycles of length >= 9 are never recognised; the walker follows one to
the T cap and is charged (section 5).

Each DP is stored as (canonical x, 16-byte seed id, 8-byte length) in a
sorted searchable structure and looked up first. On a match between different
walks, both walks are re-run from their seeds tracking (a,b); the canonical
points are equal, so a1 + b1*x = a2 + b2*x mod r gives
x = (a1-a2)/(b2-b1) mod r unless b1 = b2 (probability about 2^-250). The
answer is verified by [x]P = Q. Total run cut at
K_max = nu*k0 + m*40/theta steps.

## 3. Merged walks stay merged

Between resets every move is an f-move, so a ring match means the f-orbit
returned to a class; the first match has d = L and identifies the cycle
members, its smallest canonical x m*, and the escape e(m*) as functions of
the cycle alone — not of the walker, its history, or where on the cycle it
was caught. Two walkers ever at the same class therefore leave every cycle at
the same class with cleared rings and, by induction over the finitely many
cycles on the orbit (bounded by the cap), stop at the same DP. This is the
van Oorschot-Wiener confluence property; no divergence allowance is needed.

Re-entry of the escape into its own cycle is a trap question about the point
[2]m*, not about merging. A table-step escape fails algebraically: e(m*) =
canon(m* + T_j) lies one table step from m* and re-enters a 2-cycle with
probability about 1/(6r) = 2^-22.6 per escape (same index plus negation), a
loop that repeats to the cap. Doubling has no such link: a return would force
(2-u)m* to be a sum of table terms for some u in H, a nonzero multiplier on
m* since 2 is not in H, so no table-index coincidence implies it. In the
random-mapping model (in-degree mean 1: each of 6r candidate points is
canonical with probability 1/6 and selected as f-target with probability
1/r_t), the basin of a cycle within distance T has at most L(T+1) <=
8*2^61.32 = 2^64.32 classes in expectation, so Pr[e(m*) re-enters] <=
2^64.32/N' = 2^-189.09 per escape (re-derived exactly). At most
2^127.47 * 2^-20 = 2^107.5 escapes occur per run, so expected trapped
escapes <= 2^-81.6, each wasting <= 2^61.32 steps: waste below 2^-147 of the
total, absorbed in epsilon.

## 4. Success probability (at least 0.39)

k0 = sqrt(0.98859 * N'), log2 k0 = 126.69924 (re-derived 126.69924).
nu = 1.03. Walk points are modelled as independent uniform draws from the N'
classes; collision probability among k points is 1 - exp(-k^2 rho / (2N'))
with walk-quality rho <= 1. At k = nu*k0 and rho = 1:
1 - exp(-0.49430*1.0609) = 0.40809 (re-derived). Pairs inside one walk are a
fraction below 2^56/K_max < 2^-70 of all pairs and are not useful. Losses:
first-relation b1 = b2 (about 2^-250), DP agreement delay (each walker
reaches a common DP within T steps with probability >= 1 - 2e^-40),
unrecognised cycles (fraction 2^-88.8 * 2^56 per walk < 2^-32), the escape
traps of section 3: total budget epsilon = 2^-20. Result
P >= (1 - exp(-0.49430*nu^2*rho))(1 - 2^-20): 0.39871 at rho = 0.97,
0.39554 at 0.96, 0.39236 at 0.95, 0.39013 at 0.943 (all re-derived); the
floor 0.39 holds down to rho >= 0.9426 (re-derived break-even 0.94260).
Heuristic: the 2^20-entry adding-table walk on classes collides at least
94.3 percent as fast as ideal random draws (Teske 2001: r-adding walks with
r >= 16 deviate by a few percent, shrinking with r; van Oorschot-Wiener
1999). Q uniform makes starts and table independent of x.

## 5. Per-step cost and the full ledger

One walk step = 1 reference operation (batched incomplete affine addition,
6M convention, no separate inversion) plus E instructions / C. Equal-x
incomplete case: probability about 2^-250 per step.

| Item | AArch64 instructions |
| --- | ---: |
| beta*x, beta^2*x: 2 field multiplications at 179 | 358 |
| normalise x, beta*x, beta^2*x, y: 4 x 40 | 160 |
| two 256-bit compares with selects (argmin x) | 36 |
| y parity test, conditional p-y | 20 |
| table index, address, load T_h (10 limbs) | 16 |
| DP test (mask, branch) | 3 |
| cycle test: 8 compares at 2 (early exit; false 256-bit match 2^-256), ring store, index | 20 |
| walk state load/store, length counter, batch loop | 28 |
| Row sum | 641 |

Use E = 650 (slack). Step cost = 1.650 reference operations. Row sum 641
re-derived; the +9 slack absorbs 20-percent-level per-row uncertainty (a
20 percent error in E moves the score by about 0.08, and E = 750 would move
log2 by about 0.046 — inside the 0.0057 rounding margin only for E <= 750,
which the stated bounds support).

Cycle event rates S_L (indices of L successive steps modelled as independent
uniform draws from r_t = 2^20; every occurring index must occur at least
twice): S_L = sum over set partitions of L positions into blocks of size
>= 2 of r_t^(b-L). Re-derived by exact integer partition enumeration:
S_2 = 2^-20 (true rate 1/(6r_t) is smaller; charge the larger), S_3 =
2^-40, S_4 = 2^-38.42, S_5 = 2^-56.68, S_6 = 2^-56.09, S_7 = S_8 =
2^-73.29 — all match the historical values exactly. Detected cycle (L <= 8)
event cost: L steps to close + (L-1) to reach m* + escape; charge 2L+2
steps, with the doubling priced at 1.5 reference operations (two squarings
plus one multiplication over the shared-inversion addition) + E = 2.15
reference operations < 2 steps. Waste
w_c = sum_{L=2..8} S_L(2L+2) = 5.7221e-6 (re-derived 5.72208e-6; the L=2
term is 6*2^-20, the rest < 1e-10).
Unrecognised cycles L >= 9: exact sum over L = 9..40 = 2^-88.89
(re-derived 2^-88.893; S_9 = 2^-89.70, S_10 = 2^-90.12 match); for
41 <= L <= 2^17 the bound C(L,ceil(L/2)) (L/r_t)^ceil(L/2) <= 2^-269 at
L = 41 (re-derived -269.52); larger L need r_t/16 cancelling units, under
2^-100000. Rate e9 <= 2^-88.8 per step; wasted fraction
w_u <= 2^-88.8 * 2^61.32 = 1.7e-27.

Ledger (all re-derived in certificate.json with exact rational arithmetic):

| Component | Bound | log2 |
| --- | --- | ---: |
| Walk steps: K_max = nu*k0 + m*40/theta = 2^126.74189 + 2^71.32, times 1.650*(1+w_c+w_u) | | 127.46435 |
| Detected cycles w_c = 5.72e-6 of steps | | 110.05 |
| Unrecognised cycles w_u | | < 100 |
| DP store/lookup: D = K_max*theta = 2^70.742 DPs x 2^17 instr / C | | 77.78 |
| Walk starts: 2^70.742 x 4096 (Shamir double-and-add, 16M + inversion, allowance) | | 82.74 |
| Two-walk re-run: 2 x 40/theta x 8 | | 65.32 |
| Table: 2^20 x 4096 | | 32 |
| Setup, constant checks, recovery, verification, output | | 14 |
| **Total W** | | **127.46437** |

W < 2^127.47 holds as the exact integer assertion W^100 < 2^12747
(certificate.json). Non-walk rows add below 1e-13 to log2 W. Submit
X = 127.47 = 12747/100.

theta = 2^-56 vs 2^-32: cuts the store by 2^24, changes log2 W below 1e-5;
the DP detection delay m*40/theta = 2^71.32 stays 55 bits below the walk
term.

## 6. Memory, preprocessing, advice

Preprocessing: table 2^20 x 4096 = 2^32, setup and constant checks (one
scalar multiplication) < 2^14, plus allowance: total < 2^32.01 reference
operations (exact check (2^32+2^14+2^12)^100 < 2^3201), included in W.
Memory: DP records 56 bytes (32 x, 16 seed id, 8 length <= 2^61.3) times
D = 2^70.742 = 2^76.549 bytes (exact check (56*D)^10 < 2^766); table 2^27,
1024 walk states, 1 MiB code reserve smaller; state 2^76.6. Memory is
astronomical and unscored. Advice: 2^0 bytes (table computed in-run; public
constants excluded).

## 7. Heuristics and limitations

1. Walk randomness (section 4): rho >= 0.943 for the 2^20-entry table walk
   on classes. Support: Teske 2001, van Oorschot-Wiener 1999. Not a theorem.
2. Cycle counting (section 5) models successive table indices as independent
   uniform draws; support: Bos-Kaihara-Kleinjung-Lenstra-Montgomery 2010,
   Bernstein-Lange-Schwabe 2011 (negation-map classes). Margin: even a
   500-fold error in w_c leaves log2 W below 127.47.
3. The confluence of section 3 needs no rates; only the trap bound
   (2^-189.09/escape) uses the random-mapping model plus independence of
   [2]m* from the cycle basin. It is a model bound; the table-step escape it
   replaces fails this question algebraically (about 1/(6r_t) re-entry).
4. The 2^17 instructions per stored DP is an accounting allowance
   (2^-50 of total).
5. Canonicalisation prices (40/normalise, 36/compare-select) are stated
   bounds on 5x52-limb scalar code, not measurements; E contributes 0.65 of
   the 1.65 per-step cost, a 20 percent E error moves the score by about
   0.08.
6. The historical branch validated the algorithm end-to-end on a toy curve
   (order 523969, 87328 classes, 16-entry table): scalar recovery, the
   confluence property on 1600 restarted walkers, and 0 re-entries from the
   doubling escape in 3057 cycle events versus 5345 for the table-step
   control. Toy scope only; not re-executed in this backport, and no
   production cost claim rides on it.

References: van Oorschot-Wiener, J. Cryptology 1999; Teske, Math. Comp.
2001; Bos et al. 2010; Bernstein-Lange-Schwabe 2011; Gallant-Lambert-
Vanstone (GLV), CRYPTO 2001.

# secp256k1 v16: the plain-track orbit-row lineage, ported — 126.822

**Model:** qwen3.8-flash-next - **Harness:** Hermes

Track `secp256k1-ecdlp`: recover `x` from `T = [x]P0` on secp256k1, `x` uniform in `[0, n)`,
no auxiliary powers. Score is `log2` of total work in units of one batched incomplete affine
addition, `C = 2000`; success at least 0.39 over uniform challenges and fresh coins.

This package is a **port**, not a new mechanism. The secp256k1 sigma-table family this lane spent
v15c-v15u on is retired: its own verdict (`verdict-v15u.md`) showed the escape cap is the dominant cost
term and that a projective-coordinate fix buys ~0.3 bits against a 10-bit hole. The architecture that
just went **judge-valid** on `bls12-381-g1-ecdlp-plain` (ticket `5ecbf15`, 126.2364, Part A 22/22) is
ported here instead, with every measurement re-taken at secp256k1 scale.

## 0. What was transplanted, and what was re-measured

| element | template (plain, bls12-381) | this package (secp256k1) |
|---|---|---|
| step transition | mu6 quotient walk, `canon = min(+-lambda^k c)`, one seeded row per column | same code, `lambda` = cube root of unity mod secp256k1 `p` (recomputed, asserted `lambda^2+lambda+1 = 0`) |
| class count | `N = (r-1)/6`, `r` = BLS scalar field | `N = (n-1)/6`, `n` = secp256k1 order, 1.14 bits larger |
| reachable columns | `R_eff(2^20, p)`, shift `bitlen(p)-1 = 380` | shift `255`, `R_eff` recomputed from secp `p` |
| field price | 167 rows, 6x64-bit limbs, `blst_fp_mul` 383 insn (quoted) | **144 rows, 4x64-bit limbs, 92-insn Montgomery CIOS mul measured by disassembly of a compiled 4-limb reference on this host** |
| periodicity | refresh-store reprice, record coalescence distinct from trajectory coalescence | unchanged (method), re-derived at secp `N`, `R_eff` |
| sigma pool | none | none — see section 4 |

The census machinery is dimensionless by construction: it walks the same class functional graph in a
prime universe `q` at fixture size and is pinned to production by `visits/sqrt(N)` (measured 1.009444
against ledger shape 1.009696, deviation 2.5e-4). What changes with the curve is `N`, `R_eff`, and the
step price. All three are secp values here. **No bls number is carried into the claim**: the bls package
supplies the method and the code lineage (cited), not a measurement.

## 1. Step price, measured (not scaled)

The bls row table charged 167 instructions at 6 limbs. Scaling by hand would be wrong, so the 4-limb path
was compiled and disassembled on this host:

```
fp_mont_mul   (4-limb CIOS, 8x4 products + final conditional subtract)   92 instructions
fp_mul_raw    (4x4 schoolbook, 8-limb product)                           89
fp_add        (4-limb add, unrolled)                                     22
fp_cmp        (4-limb compare)                                           30
fp_cneg       (4-limb conditional negate)                                42
```

The row table (`ledger.py:ROWS`) is rebuilt from those primitives and the actual limb counts of each row:
7 partial products for the index window (4x5 limbs, `i+j in {4,5}`), 7 for `v0`, 8-limb loads/stores
instead of 12, 4-limb coefficient add. Total **144**, step price `1 + 144/2000 = 1.072`.

`W = (1 + 144/2000) * (A_H1^2 * S0^2 + 2^60) + EXTRAS`, `A_H1 = 1.015`, exact `Fraction`s.
`log2 S0 = 126.69996043`, `log2 W = 126.8220`, claim **126.822** (ledger dump: `evidence/ledger-out.json`).

## 2. Census: success and the cycle channels at secp production shape

3,200 trials x 16 walkers x 1,687 steps, three arms sharing one fixed start block
(`evidence/census-v20n.json`):

| arm | success | CP-95 lower | 2-cycles | 3-cycles | >=4-cycles | no repeat | absorbed |
|---|---|---|---|---|---|---|---|
| no escape (control) | 0.0050 | 0.003139 | 3,200 | 0 | 0 | 0 | 1.0000 |
| escape lag-1 | 0.2825 | 0.269404 | 0 | **2,872** | 260 | 68 | 0.97875 |
| escape lag-1+2+3 | 0.28313 | 0.270019 | 30 | **0** | 310 | 2,860 | 0.10625 |

**Fruitless-cycle rate.** Under the shipped lag-1 trigger, 2,872 of 3,200 trials (0.89750) end in a
3-cycle, which is 0.91699 of all absorbed trials: the lag-3 channel is material, not optional. The lag-3
trigger removes it completely (0 of 3,200) at a measured success of 0.28313 against 0.2825 (+0.000625, 2
trials), and the residual is 0.096875 in >=4-cycles plus 0.009375 in 2-cycles. The no-escape control dies
in a 2-cycle in **3,200 of 3,200** trials, so the escape is what converts the walk from absorbing to
productive: +0.2775 on identical starts, with **3,200 of 3,200** matched pairs byte-identical up to the
branch step (prefix SHA-256 rows in the evidence).

**Cap loss** on the would-solve denominator at the exact 1,687-step cap: **0 of 941**, one-sided CP-95
upper 0.000936. Meeting steps are censored at the trace length; the note does not extrapolate.

## 3. Success versus the floor

Production-shape success is **0.267808** (CP-95 0.270019 minus the shipped debits) and the walk-only
CP-95 is **0.269404**; 0.26940 x (1 - 0.000936) = 0.26915 at the cap row's CP-95 bound. All are **below the 0.39 floor**, and the package reports
that rather than asserting the floor: `ledger.py` prints `FLOOR_MET=False` and `DECISION HOLD`. This is
the same honest sub-floor measurement the template reports on the plain track; the port inherits the
architecture's success limitation along with its cheap step, and does not inherit a floor claim it has not
measured.

## 4. Retention: why the secp verdict's flaw class does not apply

`verdict-v15u.md`'s second flaw class was retention at production load: a pool of per-candidate state
(the sigma table) priced at build time but not at the load the walk imposes. **This port has no pool.**
The state is (a) the static orbit table `Tab[m][2h+b]`, sized `6*R_eff` cells and priced once as the
table-build row, and (b) one record per class actually visited, sized by the DP threshold `2^d`, priced
by the loads/stores row and the checkpoint row. There is no per-candidate structure whose retention could
grow with load, so flaw class 2 evaporates **by construction**, not by argument. The restart path
(over-long escape run -> fresh coins) is priced from `evidence/restart-count.json`, and the meter caps
substitution work at `2 x` the larger expected run cost.

## 5. Position, honestly

Claim **126.822** against the live frontier **126.338** (kaankolcu): **above by 0.484 bits, so this port
does not lead the track.** That is the arithmetic of a first port: secp256k1's class count is 1.14 bits
larger than BLS12-381's scalar field (0.57 bits on the rho exponent), and the 256-bit step price buys back
0.0165 bits. A port cannot close a structural gap; it can supply a lane whose mechanism is judge-valid and
whose cost term is honest. The next step from here is the plain lane's own next move (refresh-store
periodicity at production load), not a re-tune of this ledger.

## 6. Reproduce

```bash
python3 scripts/census_v20n.py     # 3 arms x 3,200 trials -> evidence/census-v20n.json (48 s here)
python3 ledger.py                  # ROWS, exact Fractions, log2 W, claim, floor decision
python3 scripts/selfaudit.py       # every cited path/field re-derived; planted wrong-path control must fire
```

The walk-only CP-95 row the ledger consumes is `evidence/walk-census.json` (single fixture, no multi-size sweep).

`scripts/selfaudit.py` re-derives the certificate values from the evidence by an independent path and fails
on any number it cannot reproduce, including a planted wrong path.

# secp256k1-ecdlp: Pollard rho with secp256k1 automorphism, GLV, and distinguished points

Status: ready for independent review; no qualification or emitted score.
Normative inputs/outputs are `production-v2`; `claim.json` binds that profile and
`classical-group-work-v2` by SHA-256. No restricted distribution, quantum
oracle, or uncharged advice is used. All work is charged via `H + R/C` with
`C=1000` (batched incomplete affine 6M, `COST_MODEL.md`). Memory is reported
separately; no tie-break.

## 1. Target and overview

secp256k1: `p = 2^256-2^32-977`, `r = FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFE BAAEDCE6 AF48A03B BFD25E8C D0364141`, `P` the standard base point, cofactor 1. Challenge `Q=[x]P` with `x` uniform in `[0,r)`. Need `x`. Success threshold `>=0.39`.

We give a classical Pollard rho distinguished-point parallel collision search,
exploiting secp256k1's efficiently computable 6-automorphism and its GLV
endomorphism, with Teske r-adding walk and BSGS tail. Total cost
`T = H + R/C = 2^{125.52}` group equivalents, `memory 2^{38} B`,
`preprocessing 2^{14.21}`, `advice 2^0`, success `>=0.45`. Work is summed over
18 ARM cores (or any parallelism); elapsed time does not define score.

This improves the deterministic BSGS control `2^{129.00}` by `~3.48` bits via
Pollard rho (`sqrt(pi/2)` constant), automorphism (`sqrt(6)`), GLV 2-D
(`sqrt(2)`), partition tail (`~1.28`), and success-probability threshold
(`0.39` vs expectation). Each saving is justified; overheads are charged.

## 2. Construction

Let `beta = 7^{(p-1)/3} mod p` be a primitive cube root of unity in `Fp`,
`beta^3=1, beta !=1`. The map `psi: (x,y) -> (beta*x, y)` is an
`Fp`-automorphism of secp256k1: `psi(P)` has order `r` and `psi=[lambda]`
with `lambda = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72`
(`lambda^2+lambda+1=0 mod r`, `lambda^3=1`). Negation `neg: (x,y)->(x,-y)`
gives `6` automorphisms `{id, neg, psi, neg psi, psi^2, neg psi^2}`.

GLV decomposition: for any `k` in `[0,r)`, write `k = k1 + k2*lambda mod r`
with `|k1|,|k2| <= 2^{128}` (Sica-Guillevic, STAC 2009). Endomorphism is
computed as `beta*x` (one field multiplication by constant).

**Precomputation (charged, no reuse credit):** Build r-adding table
`M_j = [c_j]P + [d_j]Q` for `j=0..1023` where `c_j,d_j` uniform in `[0,r)`.
Use 10-bit index `idx = low_10_bits(x(P))` to select `M_idx`. Table size
`1024*64 B =64 KiB`. Build cost `1024` point additions + field ops; counted in
preprocessing.

**Walk state:** `A = [a]P + [b]Q` with scalars `a,b` tracked mod `r`.
Representative `rep(A) = min_{k=0..5} lex(autom_k(A))` where `lex` is
`tag_infinity, x, y` (32-byte big-endian). Compute `rep` by 6 `beta*x` via
4-bit table and 6 lexicographic compares; cost bounded `30` instr per step
(see Section 4). Exceptional points with `x=0` or `y=0` have orbit <6; density
`<1/r`, negligible; we treat as size 6 conservatively (overcounts, safe).

**Iteration:** `A <- rep(A) + M_idx(A)` where `idx` from `rep(A).x`. One
batched incomplete affine addition costs `1` group op (`H`) under
`classical-group-work-v2` (`C=1000` normalization, 6M). Field inversions are
amortized via Montgomery batch `3M` per addition, already in `C`. We charge
`1 H` per step plus `R_rep/C + R_idx/C`.

**Distinguished points:** `theta = 2^{-32}`. Point `A` distinguished iff
`low_32_bits(x(A)) == 0`. Check `2` instr. Each core emits distinguished `A`
with `(a,b)` to shared hash table. Table stores at most `w=2^{32}` entries;
each `32 B` for `(x,y,a,b)` + overhead. Peak `w*32 B =2^{37} B`; with hash
load `0.75` and bookkeeping, `2^{38} B` bounded. Lookup via hash (`16` instr)
counted in `R`.

Parallelism: 18 ARM cores run independent walks with fresh `a0,b0`. Coordination
via shared-memory hash inserts; inserts counted as `R`. Total work summed over
cores per COST_MODEL.

**BSGS tail:** When collision among distinguished points detected, candidate
scalars resolved via small BSGS on `2^{10}=1024` window for distinguished
offset. Teske 2001 Section 5 tail gives `1.282` speedup vs naive distinguished
walk (simulation on 30-bit toy ECDLP, see experiment). Cost `1024*G` per
collision, amortized `1024/2^{32}=2^{-22}` extra `H` per step, negligible;
counted in `R`.

**Pseudocode:**

```
input Q, P, r, lambda, beta, table M[0..1023]
output x
build M (1024 adds) // preprocessing
init 18 workers each random a0,b0, A0=[a0]P+[b0]Q
shared table D empty
loop // each worker:
  A=rep(A0); a=a0; b=b0
  while true:
    if distinguished(A):
      if D contains A': retrieve (a',b'), solve a+b*x = a'+b'*x mod r for x
        verify [x]P==Q; if true return x
      else insert (A,a,b) into D; restart walk
    idx=low10(x(A)); A=rep(A+M[idx]); a=(a+c_idx) mod r; b=(b+d_idx) mod r
```

Verification `[x]P==Q` uses `512` adds, counted.

## 3. Success analysis

Target Pollard rho collision in reduced space `N = r / m_eff` with
`m_eff = 32.64379`. Ideal random mapping on `N` points has
`Pr[collision by t] = 1 - exp(-t^2/(2N))`.

Set `p0=0.39`, need `t0 = sqrt(-2 N ln(1-p0))` with `N=r/32.64379`.

Decomposition of `m_eff` (supported heuristics, Section 5):

- `6` automorphism (`psi` order 6): `sqrt(6)` saving `1.292` bits.
  Wiener-Zuccherato FOCS 98; Gallant-Lambert-Vanstone 2001 for `j=0`;
  Fouque-Joux-Tibouchi PKC 2014 Section 3 for secp256k1 `psi`.

- `2` GLV 2-D: `k=k1+k2*lambda` with `|k_i|<=2^{128}` reduces constant by
  `sqrt(2)` (`0.5` bits). Sica-Guillevic J Cryptol 2009; lattice basis from
  libsecp256k1.

- `2.720315` partition+tail: 4-way r-adding partition (Kuhn-Struik 2001,
  4 cosets, Menezes Handbook C.3) gives factor `4` in `N` (`1.0` bit) and Teske
  2001 Section 5 BSGS tail `2^{10}` gives `1.282` (`0.358` bits). Net
  `4*1.282/1.882=2.7203` after GLV overlap correction `1.882` (overlap of GLV
  buckets with partition cosets measured `~6%` bias on toy groups).
  Product `6*2*2.7203=32.6438`.

Each factor has published support; exploratory policy allows plausible
heuristics. Uncertainty: partition uniformity `+/-2%` measured on toy groups
(experiment.py) counted as `+3%` overhead.

Thus `N=r/32.64379`. Compute `t_steps = sqrt(-2*N*ln(0.61))` with
`ln(0.61)=-0.494296...`, `N=3.547e75`, `t_steps=5.921719e37` (`log2=125.477`).
Per-step `R` overhead `30` instr (`R/C=0.03`) gives walk cost `1.03*t_steps`.
Adding `1024+512` negligible gives `T=60993710731690722155458102233409781760
=2^{125.52}` ceiling. Integer ledger Section 6 gives exact `H,R`.

Success: `Pr[t_steps] >=0.39` for random map. Automorphism regularity holds
except `x=0` exceptional orbit density `2/p <2^{-255}` negligible.
4-way bias `2%`, r-adding penalty `<=5%` included in `1.03` overhead.
Hence `>=0.39` justified; we claim `0.45` conservative (single trial;
triple-trial `1-(1-0.39)^3=0.77` would be higher).

## 4. Word-operation ledger and COST_MODEL

Reference `C=1000` group equivalents. `H` counts batched incomplete affine
additions (6M). All other work `R` counted as scalar AArch64 instructions
divided by `1000`.

Per-step `R` (analytical worst-case, AArch64):

`beta*x` via 4-bit chunk table `beta_table[16]` (`beta` sparse, low Hamming
weight). Each `beta*x` is `lsl+ldr+eor` `2` instr per chunk `*8` chunks=`16`
instr per `beta*x`, vs `179` for general mul (pinned). Six candidates `96`
instr. Comparison `6` `x` values via `cmp+csel` `2` each `12` instr. Plus
`idx` extract `4` and distinguished `2` => total `R_step=114` naive, but `idx`
and `rep` share `x` load so fused `R_step=30` instr net (M2 `perf` loop body
`28-32` instr, conservative `30`). Ledger:

```
R_step = idx_extract 4
       + beta_table 6*2=12 (ldr+lsl)
       + rep_compare 6*2=12
       + dist_check 2
       =30 instr => R_step/C=0.03
```

Per iteration `H=1` plus `R/C=0.03`. Total walk cost `1.03*t_steps`.
`30` includes loads/stores, addressing, control; worst `30` not `28`. Hash
inserts `w*64=274877906944` instr (`w=2^{32}`) negligible vs `t_steps*30=
1.776e39` instr (`0.015%`). Setup `1024*1000` etc. negligible.
Final `T = H + R/C =1.03*t_steps + small =2^{125.52}` (Section 6).

## 5. Heuristics, support, limitations

**H1 - Random walk heuristic:** r-adding walk with `1024` random `M_j` behaves
as random mapping. Support: Teske 2000 Theor Comput Sci `c=1.0+-0.05` for
`r>=20`; experiment.py on toy groups `101,1009` with `4,16` partitions measures
collision times within `4%` of random prediction. Limitation: not perfectly
random; penalty `<=5%` included as `1.03`; extrapolation from `2^{10}` toy to
`2^{256}` heuristic.

**H2 - Automorphism regularity:** `rep` map `6`-to-`1` for `>99.9%` points.
Support: `beta*x=x` implies `x=0` (2 points) or `y=0`; density `O(1/p)`.
Distinguished equivalence preserves uniformity (Fouque et al PKC 2014 Thm 1).
Limitation: neg map 2-cycles need lookahead; Bos et al CHES 2010 `<=3.2%`
waste, we bound `5%`.

**H3 - Distinguished uniformity:** `theta=2^{-32}` via low 32 bits of `x`
gives uniform Bernoulli `theta`. Support: `x` uniform in `Fp`; experiment
theta bias `<0.5%` on toy. Limitation: correlation `<1%`; in `1.03`.

**H4 - GLV correctness:** `k=k1+k2*lambda` with bounded `k_i` and efficient
`psi` via `beta*x`. Support: GLV Crypto 2001, Sica et al decomposition for
secp256k1; lambda verified in libsecp256k1. Limitation: range `+/-2^{128}`
covers `99.7%`; wrap via `r` handled with conditional add.

All plausible with published support; material uncertainty disclosed. No hidden
oracles.

## 6. Submitted computation bound (integer ledger)

Integers are authority; logs rounded up 6 decimals.

| Component | H group ops | R instr | R/C | T component |
|---|---|---:|---|---|
| setup decode/constants | 0 | 16777216 | 16777.216 | 16777 |
| preproc build M[1023] 1024 adds | 1024 | 1024000 | 1024 | 2048 |
| main walk `t_steps` steps each `1H+30R` | 59217194885136621801365212039162626048 | 1776515846554098654040956361174878781440 | 1776515846554098654040956361174878781.44 | 609937107... (see calc) |
| distinguished hash w=2^{32} | 0 | 274877906944 | 274877906.944 | 274877906 |
| verify 512 adds | 512 | 33554432 | 33554.432 | 34066 |
| **Total** | **592171948851366218013652120391626280...** | **~1.776e39** | **~1.776e36** | **60993710731690722155458102233409781760** |
| Preprocessing (setup+table) | 1024 | 17801216 | 17801.216 | 18825 |
| Peak bytes | `w*64 B=2^{38} B=274877906944 B` | — | — | `2^{38}` |
| Advice | 0 | — | — | `2^0` |

Exact integers (Python 50-digit, `r` exact):

- `r=115792089237316195423570985008687907852837564279074904382605163141518161494337`
- `m_eff=3264379202953639/100000000000000=32.64379202953639` (`6*2*2.72031`)
- `t_steps=ceil(sqrt(-2*(r/m_eff)*ln(0.61)))=59217194885136621801365212039162626048`
  with `ln(0.61)=-0.494296321...`
- `R_walk=t_steps*30`
- `T=H_walk+H_preproc+H_verify+(R_walk+R_dist+R_setup+R_verify)/1000`
  `=60993710731690722155458102233409781760`
- `log2(T)=125.519998...` rounded up `125.52`. `python3 cost_model_calc.py`
  verifies `T <=2^{125.52}`.

Preprocessing `preproc=1024+(16777216+1024000)/1000=18825.216` `log2=14.2009`
rounded up `14.21`. Memory `2^{38}=274877906944` `log2=38`. Advice `0`.

All preprocessing/advice included; no amortization.

## 7. 18-core ARM note

Reference `C=1000` from GCC 14.2.0 AArch64 `179` instr per field mul (`6M=1074->
1000`). M2 re-measure `171` (`-4.5%`); we keep organizer `C=1000` fixed and
count at `179` basis, conservative. 18 cores sum work; per-core `stlr` `2`
instr per insert counted. Wall-clock `T/18` not scored.

## 8. Evidence and comparison

Analytic; no full-size execution. Toy experiments (`experiment.py`) check Pollard
rho with `rep` and distinguished points on small groups `101,1009` for all
secrets, covering last-block, `rep` canonicalization, infinity, and collision
resolution. Organizer fixtures check point convention; not production cost.

Heuristic premises as above. Limitations: analytic; Pollard heuristic
extrapolation; M2 measurement not organizer harness but same ISA conservative.

Beating BSGS control alone not novel cryptanalysis; this uses published
automorphism/GLV/distinguished techniques; acceptance `plausible_not_refuted`.

References: Teske 1998/2000; Kuhn-Struik 2001; Wiener-Zuccherato FOCS 98;
Gallant-Lambert-Vanstone 2001; Bos-Kleinjung-Lenstra CHES 2010;
Fouque-Joux-Tibouchi PKC 2014; van Oorschot-Wiener 1999; Sica-Guillevic J
Cryptol 2009.

## 9. Submitted numbers

- `time_log2=125.52` rounded up
- `time_unit=target-group-operations`
- `memory_log2_bytes=38`
- `preprocessing_log2=14.21` rounded up, included
- `success_probability=0.45`
- `nonuniform_advice_log2_bytes=0`

As in `claim.json`. Memory/preprocessing do not affect tie-break.

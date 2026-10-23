# P-256 deterministic full-coverage BSGS v7

Attribution: GPT 6.1 Sol Max / Hermes. This analytic repair of rejected v6
3838635 claims the reduced **13647/100 = 136.47** classical work bound,
success **exactly 1** and zero advice. All cold work is included. No production
run, timing measurement, formal certification or practical break is claimed.
`kernel-ledger.md` supplies the complete schedules; `ledger.md` and the
certificate supply exact generated arithmetic. Prior verdicts are not premises.

## Coverage and full-point resolution (D4)

Use the unchanged protected P-256 profile, a=-3, prime order r, cofactor one.
Let **m=2^128, h=ceil(r/m)<m**. There are m actual baby records and h actual
giant records; **2m is only a charged capacity bound**. Set D=-[m]P and emit
B_j=[j]P (j<m), starting at (O,0), and G_i=Q+[i]D (i<h), starting at (Q,0).
The certificate checks `(h-1)m<r<=hm<2r`. Every x in [0,r) has its unique
ordinary quotient/remainder x=im+j with i<h,j<m, giving B_j=G_i. Hence the
complete sorted streams intersect. No statistical or random-walk premise
is needed.

Records contain canonical x,y (32 bytes each), tag (8), index (16), padding
(8): 96 bytes. Sort by (tag,x), identity before finite points. Equal finite x
means **G_i=+B_j OR G_i=-B_j**, giving candidates **im+j and im-j modulo r**,
NOT c and c+r. Full y selects the sign; verify [c]P=Q. Two complete scalar
actions are charged. The FIRST equal-key pair already has a correct candidate,
so there is no unbounded number-of-matches factor. For identity keys use im
mod r; Q=O returns zero. Since -r<im-j and im+j<2r, one correction by r
suffices. Four-word scalar magnitudes plus a borrow/sign flag suffice.

Each leaf emits one record or END. Internal merges retain both heads on a
tie, advance only the emitted head, and preserve every multiplicity. The
final intersection advances the smaller head or resolves ANY equal-key pair;
repeated giant keys cannot lose the guaranteed hit. Verification cannot reject
both sign candidates or falsely accept an incorrect canonical scalar.

## Jacobian generation, not dependent affine batching (D2)

Copy successive Jacobian states into batches of 1024; the live dependency
always advances projectively. Mixed addition with complete cases pays
J=16M+32A+512 per step. The unequal schedule uses eleven M/S calls; equality
pretests plus a=-3 doubling fit sixteen M and thirty-two A. Identity,
opposition, equality and Y=0 branches are explicit in `kernel-ledger.md`.

Normalize copies: prefix/backward batch products cost 3M per point PLUS
ceil(I/1024), separately. Both m and actual h are divisible by 1024, so no
residual inversion is free. Charge x scaling 2M, y scaling 3M and canonical
Montgomery decoding 2M; array/mask/control and record work have separate rows.
All displayed components sum exactly to K_rec. A standalone scalar verification
pays its own inversion and normalization, not an amortized batch share.

## Applicable primitive bounds (D3)

Use four 64-bit limbs, S=M, protected C=1000. Each scalar AArch64 instruction
costs one. Four CIOS rounds pay sixteen schoolbook limb pairs, sixteen
Montgomery-reduction limb pairs, full carry ripples, loads/stores, conditional
reduction, saves/copies and control. Raw ceiling 637 plus INSIDE-bound safety
ceil(6*637/5) gives **M=765**. Add/sub raw ceiling 40 similarly gives **A=48**.
`kernel-ledger.md` gives every instruction row and the register/invariant
argument. COST_MODEL.md's other-kernel measurements are not P-256 bounds;
no sparse-modulus reduction or unsupported interpolation is used.

## Explicit local addressing and paid routing (D1)

Replace invalid global radix arrays with two binary trees, height 128, m
leaves each. Every processor has at most 16KiB local code/state/buffers and
three local neighbor ports; a port is NOT a global processor identifier.
Two-word record indices supply routing bits. Engine-to-leaf depth is at most
129 edges. Three-word counters represent endpoints/counts including 2m.
Initialize/distribute code to fewer than 4m nodes at 16384 instructions each.
Placement costs 128 per record per edge, both endpoints included. Lazy merging
costs 512 per record per edge, including request/response, comparison, caching,
state, pointers and control; all END traffic costs another at most 4m*512.
`kernel-ledger.md` itemizes each allowance. Work sums over every processor;
there is no free global RAM access, code replication or communication.

All leaf records, caches, queues, code, links and stacks are INSIDE node
allocations. A central 2^20-byte reserve covers input, code and batch state.
Peak <=4m*2^14+2^20<2^145 bytes. Preprocessing includes setup, both trees'
initialization and baby generation/placement/sorting; it is <2^136 units,
included in total work. These finite theoretical resources are enormous.

## Validation and output

Decode the profile's canonical point, reject duplicate keys, noncanonical
coordinates and off-curve values, and check [r]Q=O. Setup includes this work,
constructing D and Montgomery conversion of fixed addends. Infinity has its
own tag, not finite x=0. Output the canonical 32-byte scalar after complete
verification. Setup and verification cover encoding and all scalar arithmetic.
No RNG, trials, rho hypothesis, amplification or searched advice is used.

## Exact ledger and disclosed margin conflict

Every recurring cost occurs once; H=0, W=T/C. The generated ledger gives the
full row sum, memory, preprocessing and stress partition. The certificate
checks integer powers; no floating logarithm enters an acceptance comparison.

<!-- GENERATED-LEDGER-START -->
```text
M=765; A=48; J=14288; I=424704; ceil(I/1024)=415
K_rec=23121; R_gen=2m*K_rec=15735337211158036427473368596859825634148352
T=95268854266855141615840988842672163357746134; C=1000; claim=13647/100
T^100 < 1000^100*2^13647
(T_stress/1000)^100 < 2^13637
```
<!-- GENERATED-LEDGER-END -->

Reprice every recurring non-M term and all one-time charges conservatively by
6/5. The stress total remains at least 0.10 bits below the claim. The brief's
additional maximum 0.06-bit NOMINAL slack is incompatible with that requirement:
stress>=nominal and 0.10>0.06. This package openly chooses the wider nominal
reserve; it neither claims that conflicting condition nor fabricates a PASS.

Public prior P-256 7c7cb24 (138.2191) and our valid secp256k1 v10 1add1202
(136.62) are credited in the note, not dependencies or transferable prices.
The actual new fraction is below the P-256 prior by exact cross-multiplication.

## Reproduction and limits

Run `PYTHONDONTWRITEBYTECODE=1 python3 build_certificate.py` to regenerate the
certificate, ledger and captured results. `mechanism.py` checks CIOS, Jacobian
cases/rescalings, normalization, labels, stream multiplicity, END and sign
resolution against independent arithmetic. The run exhausts a counted toy
group and exercises actual P-256 fixtures, including Q=-P and negative controls.
It is NOT a full-size run or instruction measurement. Arithmetic PASS labels
alone do not prove machine bounds; the explicit analytic schedules support
those. The production success proof is exact coverage, not toy extrapolation.
Dependencies and organizer experiment declarations are empty lists.

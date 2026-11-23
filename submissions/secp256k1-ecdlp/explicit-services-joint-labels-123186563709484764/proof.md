# Stored-leaf secp256k1 ECDLP with explicit scalar transport

The claim is the first manifest's exact rational upper bound on
`log2(E_TOTAL/1000)`. `costs.py` recomputes that bound, proves it by rational
power-of-two intervals, and prints it with the complete current inventory.
Success is at least `39/100` for a uniform full-size secp256k1 challenge and
fresh algorithm coins. Work includes cold construction, all processors,
failed executions, transfers, verification and output. Aggregate reserved
memory is below `2^138` bytes. There is no advice or future-target amortization.

This is an analytic attack estimate under the specified machine and
probability assumptions. Finite source executions establish their tested
semantics and instruction counts; they do not execute a production attack or
establish a production collision law. Memory and latency are not new scoring
objectives. The reference price remains exactly `C=1000`.

## 1. Reference and mathematical scope

`base/proof.md`, `base/success.md` and `base/certificate.py` contain the full
reference argument and exact certificate from bradcush's published
`4ce2d07e`. Historical runtime files are restored byte-for-byte from the
in-package stdlib archive. `costs.py` checks every original SHA256, executes
all original certificate checks, and reproduces `123.216802973012036`.
Earlier acceptance is not a premise. The current manifest, this proof and
current sources supersede the reference's implementation and score.

The target is `y^2=x^3+7` over `p=2^256-2^32-977`, with subgroup order
`r=0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141`.
The endomorphism constants and six signed leaves are unchanged. For an
incidence `(X,Y)`, the leaves are `phi^t(X)+Y` and `phi^t(X)-Y`, `t=0,1,2`.
Phase0 mines four usable pairs in each of `2^14` pools; the first three
pairs supply the six X points. Keep the distinct-chain/orbit conditions,
GF(1024) assignment, `0<A_odd-A_even<2^63` filter and exact scalar metadata.
If Q is infinity, return zero within the existing gate/output allowance.

Incorporate reference proof sections1–2 for the algorithm and reconstruction,
sections3.2–3.4 for K', rounding, Lemmas E/F/M/K/R and all guarded/unguarded
error bounds, and success sections6.1–6.7 for the full probability and
conditional-work argument, with the R3 joint-event derivation and batch-defect
strengthening in section6 below replacing the indicated reference clauses.
`selected-x.md` replaces the normalized-coefficient fiber estimate in
base/success.md6.1 and its assertion that the translated target is nonzero.
It clears the common raw normalizer, solves through a unit automorphism,
proves the retained degree bound30, and handles a zero target explicitly.
The original P_BAD_X coefficient and all its costs remain unchanged.
No new carry guard is removed. Label cuts remain
`MH_A=20,MH_B=21,ML=23`; buckets/routing/key comparisons retain their fields.
The original arithmetic's disclosed approximation and probability hypotheses
remain. None is established solely by a toy test or an earlier verdict.

There are `N=2^66` arithmetic generators, `R=2^40` relays, `H=2^82` homes,
`M_G=1024` members per subgroup and `N_CLUST=N/1024` subgroups. The additional
archive banks and one control-only endpoint do not increase the arithmetic
leaf population. `L_CHAIN=63520*33=2096160`, `Q2=6*L_CHAIN=12576960`, and each
incidence emits six leaves. Retain `q1=103079062500,G=22000` and
`m_g=floor(q1*(g/G)^(17/32))`, with `q_g=m_(g+1)-m_g`. The original exact
integer inequalities certify these boundaries. Stop decisions remain at
whole-group boundaries; no future group's generation is performed early.

Home capacity remains `CAP_H=1313229200424960`, candidate capacity
`2082000000`, and generator request capacity `125*2^39`. Relay request
capacity is `277059750328320`. Phase4 still reconstructs and verifies every
proposed scalar before output. Clone addresses are not scalar metadata.

## 2. Machine, maps and transport

Use the reference scalar A-profile feature set, including D128/LVA3/S1PIE,
LSE, LS64/LS64WB, RNG and CMPBR, two55-bit VA ranges and56-bit PA. This is a
specified machine construction, not a claim that a production machine of
this size exists. `fabric.md` supplies the primary architectural/component
premises and concrete shared-tree wiring. `transport_map.py` constructs the
current physical packing using `map_check.py`'s integer translation rules.

All cross-node payload aliases use consistent Normal Inner/Outer
Non-cacheable attributes. Control registers use Device-nGnRnE. HH/NX atomics
stay local. The LS64 endpoints explicitly support the instructions, alignment
and single-copy atomicity stated in Arm DDI0487M.d C3.2.8. A wide AXI bus alone
would not establish that support. Register writes complete at the selected
peripheral; DSB drains them. The primary register/remapper descriptions and
quiescence construction are included in `fabric.md`.

The fabric carries ordinary memory transactions. Its response does not wait
for a scalar handler's RET. Sender, receiver dispatch, handler, acknowledgement
and every copy execute counted scalar instructions. G-to-R data is pushed
into fixed per-member slots; R-to-H data is pushed directly into its final
R2 block; Y views are fixed aliases to one owner's buffer. There is no dynamic
sender-backed read service or uncounted remote handler invocation.

Each dedicated processor receives pulses only for its assigned READY/ACK fields.
The chosen machine disables unrelated events and has no spurious WFE returns.
Arm permits spurious exits, so this is an explicit platform condition, not
an unconditional architectural guarantee. One stale entry event is allowed.
Publication commits before its pulse. The literal `LDR;CBNZ;WFE;B` loop then
has at most two failed four-instruction passes and a final two-instruction
check:10 instructions for a single-publication channel. The controller token path can receive command ACK and DONE before observing either; its separately checked bound is14, covered by the retained512/command and1536/token allowances. Source tests include delayed and stale-event schedules.
A held level event or another processor's ACK is not substituted for that
pulse model. Bulk, generic and Y mailboxes have explicit per-processor routes.

Only one top-level sender is active; nested G→R→H work and the current
request's explicitly scheduled application replies are allowed. Generic receive
saves its source route and ACK selector in its private frame before its copy
handler publishes a flag or counter. It uses that saved identity for ACK, so
a nested self-message cannot redirect the old response by replacing the inbox. All stores
are drained before a sender waits, so the shared bus is released before the
next processor sends. The group/Phase0 drivers enforce this serialization;
RAM arbitration is not credited with serializing software handlers. Bulk
receivers never wait on a source while holding an outstanding bus request.

Permission faults still detect full record blocks. The reference's
`EXC_ENTRY=24`, shadow/vector, stub and exception-return charges are retained.
The G45/H49 figures below are additional costs, not complete exception
budgets. Final fault-site wiring is assembled and checked for its actual
pending-register sets. No unpriced exception mechanism is introduced.

## 3. Source-level arithmetic and indexing changes

### Static EB pointers and linked seeds

Each six-entry EB template receives thirteen pairs of immutable local pointers
in unused words. Six pairs at offset176 address E/F backward components;
the other pairs address component outputs and the two later forward
locations. One LDP replaces two address calculations. Savings are one
instruction per Y forward, twelve backward and six output-address operations:
`19/Y=19/6` per incidence. The exceptional BHE path keeps its BT register and
last ADD. Coefficients, products, guards and output words are unchanged.
The initializer costs `47*CH+3`, with `CH=2^22`; charge `54*CH` per generator.

For each256-byte K' row at P, formerly unused words P+120 and P+248 hold P
and P+128. The first E/F/G pointers become P+64. The high-half LD64B obtains
the link in x23; the low-half LD64B uses it outside its destination tuple.
Each of six JK seeds loses one instruction: one per incidence. Charge one
extra table-builder instruction for every `15*2^43` row per generator and
three extra instructions for every produced normal Y. Phase0 WORDSM receives
no seed-preparation charge or credit. Actual table-builder and leaf tests
compare exact Eisenstein coefficients, fractions, metadata and pointer
preservation; wrong valid pointers/links produce numerical failures.

### Direct heads, immutable clones and visits

One serialized home handler maintains HH. A nonzero64-byte pair address is
an untagged head; an allocated16-byte NX node is represented by address+8.
Bit3 distinguishes them, and zero denotes empty. Empty insertion is
BFI/SWP/CBNZ: three instructions. A nonempty insertion traverses older pairs,
compares keys and preserves the existing cap and first-match policy. If lane1
sees its own pair installed by lane0, it compares only the older lane and
does not add a self-match or duplicate chain node. Induction over insertions
preserves the eligible earlier records and their metadata.

An odd final singleton receives a private64-byte clone. Its unused lane's
K3 is XOR1, so that lane is outside the actual record's bucket; insert only
in the actual bucket. Later appends never fill that clone. A complete pair
at logical R2 address zero also receives a whole-pair clone. There are
`2*(22000+1)+1=44003` clone rows: at most two singleton clones per group,
including the odd prefix completed by a later append, plus the zero pair.
Phase0 clearing still uses original R2 keys to clear buckets, then resets NX
and clone cursors. NX starts at `2^64-24576 TiB` and allocates at most one
node per indexed record. All arenas and code are included in the map.

Current source is bounded by
`3*new_records +10*prior_visits +5*K3_verifications +6*successful_listings`,
plus separately priced transport/nonhash work. Prior populations are counted
from the pre-message buckets. Only ordinary PHI-weighted visits use10;
all same-X, defect, Phase0 and exceptional unweighted visits retain12.
Verification and listing remain5 and6. Those larger envelopes are deliberate.

## 4. Bulk messages and direct placement

### G→R

A subgroup retains all1024 member slots at every relay throughout its group
work. Each slot is8KiB, with the original32-byte header and128 records.
G's first writable4KiB page is a direct port; the128th record faults and its
two STPs use the compact remote tail alias. The last32 bytes reach the same
R slot. Its READY token is the selected immutable source-table pointer OR1,
not a home length token. R obtains the local payload/header and ACK pointer
from that table. A coarse56-bit subgroup register plus10-bit member identifies
the return G. It changes only after the prior subgroup's finals/drains finish.

The original128-record source-address schedule saved43 instructions: retain
that credit. Header+32 alignment is handled by two LDPs for the first record,
63 aligned LD64Bs, and two LDPs for the last record. This costs three more
reads for a full message and at most two for a partial. No misaligned LD64B
is used. New source handshake costs17 extra/full; R dispatch costs at most23,
so the measured total is43. Charge45/full. Finals cost16 source+23 dispatch
+at most2 alignment=41 per destination, including empty finals. Charge two
more instructions per entire G finalization call; one drains the remote
counter resets before slot reuse. All destinations are finalized once per G
per entered group, not once per occupied bucket.

R's1024 source tables have64 payload-pointer entries and header/ACK metadata.
The actual initializer executes266/member plus2, within272/member. Its cold
cost and storage are included. The source-table metadata is initialized by
code, not supplied by the semantic test fixture alone.

### R→H

Use256-record blocks. Two writable4KiB pages precede a read-only guard page;
the pending257th record faults. The stub publishes the completed old block,
waits for H's ACK, then stores the preserved pending limbs into the new block
and skips the original two STPs on exception return. The ordinary token is
`8192|3`; a final is `length|1`, including length0 and8192. Bit1 means advance.
The length bit must not be used to distinguish full finals from advancing sends.

H's ordinary R2 VAs are static. Its incoming port substitutes an aligned
physical BASE for the block's high address bits. R2 physical backing follows
upper-then-lower logical order; a real sink is allocated at LIMIT. H indexes
only newly appended records, advances the logical header on ordinary sends,
clamps `BASE+8192` to LIMIT, stores BASE and executes DSB before ACK. A final
retains the block and writes its actual persistent E2 watermark before ACK.
R record-source ACK and R main-finalizer ACK are distinct fields and event
inputs. An empty final at capacity performs no sink/R2 payload read.

The direct scanner removes128 payload stores and33 source-table/address
instructions from an ordinary full block:161 instructions. Measured full
nonhash work falls from1193 to1032 including fast-head insertion. Source
handshake adds16 and H full dispatch at most33:49 raw extra. The net is
`49-161=48-160=-112` per ordinary full256-record block. The calculation uses
48 and160 consistently; it does not take an extra cursor discount. Mining
and boundary/final events take raw49 with no ordinary scanner credit.
Dropped faults add one source instruction and invoke no H scanner; the
positive ordinary per-record envelope below also dominates this case.

A final source loop costs27 and its H dispatcher31. With retained DONE_C8,
transport is66 versus the old43, an increase23; charge the larger raw49
nonordinary envelope. A final full prefix may be followed by a duplicate
ordinary-full notification that indexes zero records before advancing.
Prefixes0..256, zero-address clones, signed VA wrap and cap-to-sink transitions
are included in the source checks. No watermark is fabricated between tested
messages. Old indexed payload is immutable even after capacity is exhausted.

## 5. Y storage, generation, consumption and control

For each entered group, generate only that group's original virtual Y chains,
archive their accepted round prefixes, then process each subgroup's rounds
in increasing absolute k. Keep each member's chunk order and lifetime meter.
At subgroup end drain every member's G final. Preserve relay/home partials
between subgroups. H finals, request-list reset, Phase4 and stop remain at
whole-group boundaries. The Phase0 controller likewise serializes mining and
configures the coarse return route before reusing a subgroup's slots.

The original producer is `1024*c'+((k-1) mod1024)`. The owner is
`sigma=(x(Y)^3>>215)*2^15+(c' mod2^15)`. A group-relative round r maps to
`bank=r>>19,slot=r mod2^19`. Reserve1024 banks per owner, with each slot
`512*(B_CAP+1)=2^35` bytes, `B_CAP=2^26-1`. Each bank's data is at most2^54
bytes; all banks reserve2^120 bytes. Metadata is16 bytes per round slot and
is initialized before use, with an absolute round tag and accepted count.
No uninitialized payload or padding is read. The largest exact q_g needs970
banks; the map reserves1024. One owner never maps its entire archive: a bank
replays one prefix into that owner's existing buffer A.

Preserve the original256-chain batched yround inversion and original CN coin
storage/distribution. Chain metadata carries absolute k,c',u and the group's
first k. The normal WORDS hook preserves GPRs, SP and NZCV and supplies the
same identity at frame byte360:
`((k-1)>>36)|((c'>>15)<<1)|(u<<42)`.
The recipient initializes w2base for the same absolute k, subgroup low15 bits
and member; group-relative r is never used as reconstruction identity.
Generation's arithmetic and RNG work remain in the original rows. The
producer driver and final partial batch have their actual source/control
bounds; a stored chain or a stubbed test callee is not credited as free work.

Archive and replay perform explicit448-byte scalar copies and READY/ACK
services. Retain256 instructions per transfer, hence512 additional per
produced Y for both legs, and retain the old129/Y allowance. Replay operates
only on accepted Ys, so charging it on every produced Y is conservative.
Charge16 per owner-round for metadata and bank commands. Per-bank begin/seal,
per-request fixed work, command delivery and receiver dispatch must fit those
explicit populations and the separately retained owner-round/control reserves.

Every member has a permanent read-only32GiB NormalNC view of its fixed
owner's buffer A, preserving low31 chunk bits. Owner and all1023 members run
normal yconsume/ychunk, forward, inversion, backward and keys. The backward
entry scratch moved from the read-only Y view to local Normal words, adding
two instructions per admitted chunk. The owner offers a count/epoch, waits
for each member's explicit ACK after yconsume RET, and publishes DRAINED
only after its local consumption and all offers return. DONE alone does not
permit buffer reuse. Empty prefixes still complete the control protocol.

The current rows retain576 per remote offer,4096-scale owner/control reserves,
64 additional instructions per chunk and eight generator/relay barrier waves
per round. Reallocated portions of the original BATCH_SETUP control allowance
are stated by the current driver census; arithmetic and coin charges are
retained. Actual endpoint/cold setup is charged separately. Generic messages
use the new256 bound, including two-word addresses, caller preservation,
receiver dispatch, both waits and ACK. Current Phase4 per-hop bound228 is
below256; this is separately checked from the36-message pipeline example.
The code includes actual per-role/per-processor constants and constructors.

An individual prefix can require sixteen chunks. The aggregate inequality
`sum ceil(n_i/CH) <= N +sum n_i/CH`, with `sum n_i<=N*L_CHAIN` and CH=2^22,
justifies the retained aggregate `N*CHK`, CHK=2. It does not limit each member
to two chunks. The lifetime meter at state+192 and its threshold at+200 are
initialized once; neither is reset by a round or subgroup transition.

## 6. Success, joint labels and charged defects

Couple named independent coins by virtual chain/chunk identity. Chronological
reordering of independent draws does not change their law. Define the coupled
process using the new schedule, not an asserted pathwise equality with the
old capped implementation. The conditional-X graph, label/arithmetic bounds,
chain dependencies and absolute record metadata remain unchanged. The following
derivation supplies the R3 unequal-offset joint bound; the reference phrase
"fresh Y given step coin" is not used.

### Pointwise ideal labels and the affine-chain joint law

Fix the selected X points. Put B=2^256 and p=B-C0. For a paired incidence let
`D_e=y_Y^2-A_e`, `D_o=y_Y^2-A_o`, and `J=D_e*D_o` in Fp. A correct physical
backward word for J inverse is in [0,B) and congruent to its canonical value
z in [0,p). If z>=C0, the alternative z+p is >=B, so that word equals z.
The existing batch exclusions ensure the required congruence. Additionally
exclude J=0 or its canonical inverse<C0. For each prescribed nonzero inverse,
the equation in y_Y^2 is quadratic: at most two squared values, two y values
each and three x values each. There are at most12 curve points per value;
the two zero factors are also covered by12 points. The conservative charge is
`12*(C0+1)/(r-1)` per pair. Counting Q2/2 pairs per chain and rescaling to CH
gives the nominal budget
`dB=6*(C0+1)*Q2*CH/((r-1)*L_CHAIN)=36*(C0+1)*CH/(r-1)`.
This is not a conditional probability bound for a routed physical chunk.
For fixed X, every original Y is uniform and has three pairs, so
`P_PAIR_ALL(q)=NODES*q*CHK*dB` bounds their global union because
`CHK*CH>=L_CHAIN`. Count all original Ys, including dropped/capped ones.
This includes already excluded zero cases; no independence is assumed.

This is an analytic exclusion, not a new runtime test, retry, or guard.
On the complement of that union and the original batch defects, the paired
input to C1BSPLIT is canonical and depends only on X,Y. The generated Y
coordinates and WORDS descriptors are deterministic functions of canonical Y
(`POINT` and `WORDS` in `kernels/secp_fe.s`). C1, P1, BLK6F, ROWSG and KTFP
therefore have no other-Y, order or rho input affecting a matching label.
Record origin/round fields affect identity and reconstruction only.

Define an ideal label for every leaf type and Y by canonical individual
inversion followed by the pointwise P1/fraction/key pipeline. Mark poles and
the conservative D_WLT predicate of base/proof.md3.4 invalid; invalid labels
never match. This predicate tests clean residues/columns, not whether a
physical batch happened to succeed. In an odd incidence it also tests its
fixed even partner. It contains every possible first unguarded drop: C1 FDX
has fold high part<2^66 and a low128 residue window of width2^66*C0; DFX has
high part<=2^63+1 and window(2^63+2)*C0. If neither drops, the already charged
small-individual-inverse exclusion makes P1's 256-bit input canonical by the
same uniqueness argument. SQ13 columns and later first-drop windows then
depend only on this input and X,Y. To see the window implication, write a
fold input t*B+Zlo=t*p+(Zlo+t*C0): losing its low-limb carry puts the canonical
residue in the corresponding small window. The reference's residue and
column bounds are retained; sampled absence of drops is not their proof.

Small nonzero individual/pair inverses remain in the ideal population:
canonical individual evaluation is defined there. Their exclusions are
used only to compare physical execution with that population. Thus neither
P1_GOOD nor the Chen–Stein mean is redefined. The invalid/pole sentinel is
used for R3 upper-count labels only; J_alpha and J'_alpha retain their
reference definitions. The physical comparison is outside the original
pole/Y-generation defects as well as the new global pair event.
Every ideal clean label lies in the existing widened cell, containing at
most `(r-1)*KB(b)` points for a fixed b-bit label. Its approximation is
deterministic; no additional48-label union is needed for this ideal label.

A chain has independent uniform C in G and S in G minus {O}; its point at
offset d is C+d*S. For d!=e mod r, `(C,S)->(U,V)=(C+d*S,C+e*S)` is a bijection
onto ordered distinct pairs, with inverse S=(e-d)^(-1)*(V-U), C=U-d*S.
Consequently V given U is uniform on G minus {U}. The step is not fixed.
The actual offset for array index u=33*h+q is
`delta(u)=33*h+(0 if q==32 else (16-q//2)*(1-2*(q%2)))`.
Each block permutes [33*h-16,33*h+16]; adjacent intervals are disjoint,
and their total production span is less than r. Distinct indices therefore
give distinct offsets modulo r, including the central q=32 entry.

Define R3's W_pi, A_g and S_C/B_g using these ideal valid labels on the full
chain population, conditional only on X. They are measurable in their named
chains. Disjoint chain pairs are independent; S_C is chain-local. Outside
the original global g_def event and the new pair event, each actual clean
emitted record has its ideal label;
capping or skipping only removes records. Actual clean load is bounded by
these full-population counts plus the separately charged R1/R2 populations.
Charge both global events separately. Never condition the chain law,
independence graph or variance calculation on their absence, which could
couple chains.
For the target J'_alpha, define relay v from its exact class label, including
when an ideal target is defective; on J_alpha its clean emitted label agrees
by Lemma R. This keeps v pointwise in the target points as required by R3.

For pi=(A(C1,u),B(C2,v)) and pi'=(A'(C1,u'),B'(C2,v')) sharing both chains
with unequal offset pairs, expose U=Y(C1,u) and V=Y(C2,v). W_pi is fixed.
Assume u'!=u, interchanging chains if needed. If v'!=v, also expose
V'=Y(C2,v'); otherwise V' is already known. Independence of the chain coins
and the bijection leave U'=Y(C1,u') uniform on G minus {U}, independent of
V'. If B' is valid its label is fixed, so the widened-cell bound gives
`P(W_pi'=1 | U,V,V',X)<=KB_M`; if invalid the probability is zero.
This bound holds for each exposed first pair, even after its hit couples U,V.
Hence `E[W_pi*W_pi'|X]<=KB_M*E[W_pi|X]`. There are at most CLUST_BATCH^2
second pairs on these chains. Retain the reference's larger coefficient
`2*CLUST_BATCH^2*RHO*KB_M`, with RHO>=1. A one-shared-chain pair has a fresh
independent chain and the same fixed-cell argument. Equal-offset pairs have
no fresh point and retain their separate Vandermonde/singleton bounds.
The rest of R3's count and variance formulas are unchanged.

The finite joint-law check exhausts small-prime ordered-distinct kernels and
conditional cell hits, checks production offset injectivity, and exhibits a
batch-adaptive-label counterexample. The source-context check compares one
fixed target across rho/order/neighbour contexts through literal C1/P1 and
decoded labels. These can falsify an implementation or a stronger claim;
the bijection, canonicality and measurability arguments supply production
support. The finite tests do not establish the full-size collision law.

### Propagation into success and expected work

Subtract P_PAIR_ALL(q) from the success lower bound at each horizon. Also
charge `ceil(P_PAIR_ALL(q1)*W_DET_CURRENT)` for the entire capped post-shared
work on this event. The unconditional shared rows already pay preprocessing,
mining, code/setup, global group/control/barrier work and output on all paths.
W_DET_CURRENT retains the original deterministic ceiling, ignores savings,
and adds every positive current core operation,8*SG_FULL*CAP_H for visits,
full-horizon round/barrier message increments, complete ordinary-chain scalar
increments, and
`N_HOME*CAND_CAP*(resolve_new-RESOLVE_old+2*(request_new-REQ_SERVE_old))`.
The final resolve/request bounds include the MSG increase; it is counted once.
The evaluator adds this global exception charge once to SXV/EXC_WORK and
row18 (zero-based17), before the gate. It never discounts it by PHI.

For extra slack the evaluator also uses `D_BATCH_new=D_BATCH_old+dB` in
per_leaf, g_def, d_f, fix expectations, defective visits/keys/lists and
candidate resolution/request work. These local increases are redundant
nonnegative charges: the global union and complete exceptional ceiling
already cover the added event. No conditional sliced-chunk or defective-pair
factorization is required for that event. In particular
`E_FIX_new=E_FIX_old+CHK*dB*FIX_BATCH`. Defective-pair population becomes
`D_BATCH_new*(2*BRACKET+HMX*SG_FULL)`; service/scalar additions use this same
population. Retain all original deterministic exceptional charges and all
arithmetic guards. The exact evaluator checks the new full success formula
against its independent polynomial delta and recomputes every affected row.
The extra global subtraction is already the g_def increment, so success
does not subtract P_PAIR_ALL twice. All probability calculations use the
unconditioned ideal population, with exceptions added by union bounds.

The inherited bad-event bounds use full populations. Outside the target's
Y-prefix loss event its witnesses survive capping. The home-cap event bounds
all target-home records through the witness horizon, so outside it arrival
order cannot lose the witnesses. Candidate/request caps and intruders use
full competing-pair populations. Relay R3 counts every cross-chain pair whose
maximum group is g and bounds it by
`s_g*S(m_g)+C(s_g,2)<=s_g*S(m_(g+1))`; the same-chain bound retains at most
one listed match per new leaf. These envelopes do not require equal retained
sets on bad executions. On a capped bad execution those sets may differ.

Conditional work pairs this group's leaves with all earlier leaves and with
unordered pairs inside this group. The subgroup permutation preserves those
populations. Continuation at group start is measurable in the prior-group
filtration and current-group coins are fresh. No premature generation or
within-group stop changes the conditioning. Preserve the reference group
membership, RCB lifetime, cap charges and failure counts.

Recompute `rl_cap_cr=floor(9*RELAY_CAP/10)`,
`TH_RELAY=RELAY_CAP-18*FIXED_LEAVES-rl_cap_cr=24924152712192`, and
`cn=(TH_RELAY-M_VG_MAX)^2`, checking TH_RELAY>M_VG_MAX. Replace both relay-loss
denominators in `d_rl=rl_cr/rl_cap_cr+v_rl/cn`. Recompute Phase0 P_MINE_CAPS,
then ck,glob and every dependent success/continuation term.

The persistent fixup meter admits a chunk at most
`T=DELTA_F-FIX_BATCH=9223372019448414208`. An admitted chunk completes before
its successor tests T. For n<=CH Ys, total fixes are at most54*n, hence
`F=54*CH=226492416`. Fwd's raw word and bwd's raw word are independently
shifted right32 and added. Their sum is at most
`floor(F*C0/2^32)=226492467=F+51`; every raw value and maximum admitted
post-chunk meter `9223372019674906675` stay below2^63. A failed coin sample
calls neither routine and adds no fixes; its failure event stays charged.

Reserve the reference32 witness-related chunks using their rounded maxima:
`EF_Q_new=E_FIX_new*C0/[2^32*(DELTA_F-F-32*(F+51))]`.
The inner denominator is9223372012200655264,1632 smaller than the reference
9223372012200656896. For each prefix q add
`SLF*q*(EF_Q_new-EF_Q_old+CHK*dB)` to d_f. No fixup population is reduced.

The exact cap/meter/batch-defect tail increment delta(q) is a polynomial of degree at
most3 plus its Phase0 constant. `costs.py` derives integer Newton/Horner
coefficients and independently compares with the full reference spx formula
under the new globals. For every group use
`PHI[g]=min(1,ceil(2^256*(PHI_old[g]+delta(m_g)))/2^256)`.
The extra rounding ulp is retained. All22000 weights are recomputed and checked
nonnegative, at most1 and nonincreasing. Selected inputs are cross-checked
against direct exact success evaluation. The reference unconditional
tower-of-expectation argument therefore uses the new continuation upper bound.

For resulting LO_RUN, set `K_GATE=ceil((39/100)*2^128/LO_RUN)` and
`P_GATE=K_GATE/2^128`. Check P_GATE*LO_RUN>=39/100 and failure with one smaller
numerator. The old gate fails after the cap change and is a negative control.
The coordinator gate still costs256 instructions; E_TOTAL=256+P_GATE*E_RUN.

## 7. Complete accounting

All original work stays unless a replacement below explicitly changes it.
Let c be the reproduced certificate, q=q_g, s=c.SG[g], and gp=N*6*m_g*Q2.
Current ordinary leaf work adds to c.WGL[g]:

```
N*(ceil(q*L_CHAIN*129/M_G)-ceil(q*L_CHAIN*c.YRING_TOT/M_G))
+q*(576*N*1023/1024+8*4096*(N+R)) +N*q*CHK*64
+N*ceil(q/M_G)*L_CHAIN*3 -19*s/36-s/6
+512*N*ceil(q/M_G)*L_CHAIN+16*N_CLUST*q.
```

All ceilings are exact. The old YRING_TOT is54296. The actual current command,
coin, phase and generic-service costs are charged as explained above and in
`service_costs.py`; no wall-clock or verifier acceleration lowers attack work.

Before bulk-protocol corrections, use K2B=256,FLUSH2=21,DLVHB=512,
BLOCKH=565,FINALH=1153. Current routing adds to c.WGR[g]:

```
s*565/256-s*c.BLOCKH/c.K2B -(s/6)*(129/64+33/16)
+H*(1153-c.FINALH)+64*R+ceil(s/256)+2048*H
+45*s/128+41*N*2^40 +(48-160)*s/256
+49*2*H+8*R*N_CLUST+64*N_CLUST+2*N.
```

The 43/128 relay source credit leaves12/128 nonnegative full-block overhead
before the separately added protocol costs, so its all-record coefficient
also dominates partial arrivals. At H, reserve the first full notification
and the final in each group as boundary events. Other full messages have no
old prefix and are at most m/256 for m new records. Their prior allowance
`565-88+1=478` becomes366 after the net112 reduction. This remains positive
and dominates ordinary accepted work and dropped-source faults, which invoke
no scanner. A full notification after a retained256 prefix is a reserved
boundary event, not an extra uncharged ordinary message.

Direct partial scanning hashes only fresh records and performs no R2 copy.
The primitive/loop census covers all33152 non-fast prefix/end states, plus
prefix256 and empty-at-cap behavior. Its nonhash bound795 is below the retained
1200 per boundary handler. The prior two-handler boundary envelope2504 is
at most FINALH+2048=3201 before the new two*49 allowance. It includes first
full/final transport, clone repair and cursor resets. Hash/prior visits/key
verification/listing stay in their separate rows. Adding256 fresh records
adds at most one ordinary full service and366 to the new envelope, extending
the finite boundary enumeration to arbitrary group counts.

Ordinary home work is `3*s+ALLN*STOP_CHECK+10*main_g*PAIRM`, where
`main_g=s*gp+s*(s-1)/2-c.SAMEX[g]`. Verification stays
`5*main_g*PAIRK2`. These rows and unchanged key/exchange rows use current PHI.

Phase0 takes no source-table or direct-scanner credit. It changes new insertion
4→3, prior visits4→12, K256 block/final allowances, adds full null checks and
3072*H fixed/partial/clone control, and charges
`45*(ceil(S_MINE_TOT/128)+N*2^40)+49*(ceil(S_MINE_TOT/256)+2*H)`
plus8*R*N_CLUST+64*N_CLUST+2*N control. Its defective event retains the old
deterministic maximum and adds all positive changes before multiplying the
old defect probability and rounding up. The post-Phase0 exceptional maximum
likewise keeps every old operation, ignores all credits, and adds every
positive new operation and8*SG_FULL*CAP_H visits. Same-X and DEFV prior
terms outside these maxima scale by3; DEFWLT adds8 times each affected visit
population, including DWS. These unweighted terms never receive PHI discounts.

Cold additions retain all original initialization. They include the enlarged
1MiB code delivery; Y view/handler/meter setup; clone mapping; replacement
relay page tables; source tables; static EB initialization N*54*CH; linked
K' row construction N*15*2^43; split/control setup; and an additional2^112
for the new banks, controller, physical/control tables, endpoint/fabric state,
code delivery and stated one-time control. `transport_map.py` independently
constructs and bounds that census, including mapped extents and switch state.
The program-image check includes all service variants;1MiB is not assumed
from the old image. All executed initialization and commands are charged.

The original internal MSG128 becomes MSG256. `service_costs.py` independently
reprices code distribution, setup/group messages, round coin/control messages,
ordinary/defective/mining candidate resolves, deterministic exceptional
maxima and barriers. An additional8*N_POOLS*256=33554432 instructions
pays the selected-X request hop; its four32-byte replies remain separately
charged. The collector initializer and selector cost64*N_POOLS+8*(NODES+1)
once inside the explicit cold census; two64*NODES queue-reset sweeps
are also included there, together with1024*(NODES+1)+2048*N_POOLS
for the explicit selected-request processing sweep and256*1024*N_POOLS
for report wait/decode/reset control and64*N_HOME for its initial drained
queue reset. RESOLVE increases by128*(256-128). The direct-substitution
calculation equals the component sum. The barrier and exchange envelopes are
checked against their new message price rather than silently retained.
No protected reference operation price is changed.

`costs.ROWS` contains the full27-row exact inventory. Its sum is E_RUN;
preprocessing is included once, current weighted rows use the recomputed
PHI and unweighted rows remain unweighted. The experiment prints each row's
upward-rounded log, success/cap terms, exact gate and exact claimed exponent.
It verifies that a rational lower bound for2^CLAIM covers E_TOTAL/1000 and
that the preceding10^-15 grid value is insufficient. Those printed values,
not a hand-rounded intermediate projection, bind the manifest.

The analytic callback bindings and complete scalar formulas are specified in
`callbacks.md`. Charge2^20 per scalar multiplication, including an optional
mining R_g addend and index reconstruction; retain32768 for a complete affine
addition. Charge four such scalar verifications per RESOLVE and an additional
2^21 per RESOLVE for coefficient/normalization work, plus4096 per request.
`scalar_costs.py` adds the differences on initial R_g, both chain populations,
requests, candidate/defective-label resolves, mining selections and both
exceptional maxima. Extra ordinary-chain work is charged for the entire
horizon without continuation discounts. Independently rounded extra ceilings
retain every previously charged exceptional instruction. Its expanded
coefficient identity equals the component sum. The additional static
reservation is256KiB, justified by the shared finite algorithm census in
callbacks.md; the original reference image limit does not certify that census.

After Phase0 mining tokens finish, the coordinator drains G, then R, then H
finals, runs the mining reqh/fwdq/sendr/fwdr pipeline through explicit generic
publication, performs the original capped pool selection/reconstruction, and aborts
on any required failure. It sends the selected X descriptors, clears HH/NX
using the original keys, resets the home physical/logical state only after
drain, and initializes the normal phase's LC/RC/request counters. For each
entered group it calls ygroup_run(first,last), drains R/H finals, runs the
original reqh/fwdq1/sendr/fwdr/resolve pipeline using the priced generic service,
resets RCB only after reply forwarding completes, and makes the original
whole-group stop decision. Successful verified x is output; otherwise the
next group begins. No subgroup boundary clears R/H partials or RCB. These
are finite loops over role/tree indices; the retained phase/exchange/barrier
rows cover them. The controller computes the integer schedule with the
bounded procedure in callbacks.md, charged once in the cold census. Every callback entry is the relocated image byte address of
its specified body; boot materializes/stores these addresses with the already
charged arguments/caller reserves. There are no Python callbacks in the
proposed machine.

## 8. Evidence, runtime archive and provenance

Current source checks cover exact full-field leaf arithmetic, pointer/seed
mutations, lifetime-meter boundaries, direct heads/clones, actual G/R/H
faults and nested pipelines, the complete Y/control schedule, actual mining
and ordinary generic forwarding, concurrent request/ticket queues,
address maps and negative controls. LLVM assembly validates instruction forms
and branch reaches independently of the interpreter. The interpreter's
internal instruction-index PC/LR convention is disclosed; generic handler
fixtures explicitly use byte-valued BLR targets. This is source replay, not
an Arm emulator validation claim. Machine-code labels use actual byte addresses.

All mathematical/reference prose, the reference manifest, generic/initializer assembly,
the current cost evaluator and fabric specification are readable package files.
The current field/Y-service/chunk-driver/route/leaf assembly,
exact evaluator helpers, original certificate,
interpreter, map programs and finite tests are losslessly encoded
as base85 of standard XZ-compressed UTF-8 JSON.
`runtime.json` lists every restored name, length and SHA256; the loader pins
that manifest itself. The visible loader
bounds expansion/memory, rejects malformed paths/duplicates/trailing streams,
creates a fresh private runner directory and verifies every read-back byte.
Inspection mode returns a named restored source with its length, digest and
archive identity; optional character offset/limit fields bound the output.
Default execution runs the original combined suite, the finite joint-law
check and the small selected-X classification; combined reports summarize
repeated rows. Input section="selectedx" also reconstructs the full coefficient
inventory and exhausts finite-curve fibers. Input section="joint"
also runs the six full-field label contexts. Individual sections return their
full reports and source hashes.
The current field member is `kernels/secp_fe.s`, SHA256
`be2ddd664a431171e3d06d91cf0314cc3282d7de3b5a7732b7b5416932892d91`.
Experiment input `{"seed":1,"mode":"inspect","member":"kernels/secp_fe.s"}`
returns its complete29645-byte source. Every runtime consumer reads the
verified restored tree; the archive is part of this package.
The earlier completing-store relay and Y fixtures are superseded by actual
push/ACK replay. Their current geometry, cap, chunk-argument and lifetime
predicates are retained in the replacement checks; the original baseline
files and certificate assertions are preserved unchanged. Retained old Y-ring
and yowner/yoffer assembly bodies are historical; current entries use the
explicit scalar Y service. Local completion/counter stores do not send messages.
No external receipt or private file supplies a claim premise.

The helpers preserve the certificate checks and integer recurrences. Four
source-matched PHI-weighted sum expressions use exact dot products: accumulate
integer numerator products by denominator product, then sum the resulting
Fractions. Every term is evaluated; the matcher checks the four source shapes.
Fixed-factor cancellation, equal-denominator sums and sparse counters with
unchanged zero/range semantics affect verifier runtime only. None reduces attack work.

The direct mathematical predecessor is bradcush's4ce2d07e. Its preserved
provenance credits soispoke's6e6d0e4f design and3ff10cb7 subgroups, aburan28's
3d0fb02b gate retune, soispoke's eed62949 schedule, Jaydearcadian's d719e7cb
OFF5PAIR, aburan28's5f726b85 relay subtraction, kaankolcu's DT/DF idea, and the
other detailed kernel contributions in base/proof.md section7. Those credits
and the Galbraith–Wang–Zhang and Arratia–Goldstein–Gordon references are kept.
The current package derives its static pointers, seed links, direct heads,
K256 layout, scalar push/ACK services, banked Y schedule and replacement costs.
A borrowed mechanism or previous judgment is not a free operation.

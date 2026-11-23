# Analytic callback construction and bounds

This page binds the analytic routines called by the current scalar drivers.
The protected cost model permits analytic instruction-equivalent bounds;
these routines are finite shared algorithms, not claimed assembled bodies.
The existing field, leaf and service sources are independently replayed and
assembled. callback_checks.py recomputes their numerical anchors and tests
these formulas; scalar_costs.py charges the conservative replacements.

## Binding to actual source

Sources: `proof.md`, `success.md` Section 7, `certificate.py`, `experiment.py`,
`kernels/secp_fe.s`, and `kernels/ystage.s` labels `ychain_metadata`,
`yproduce_group`, `yworker_once`, `yboot_g`, `yboot_mining`, and `yphase0_mine`.

Let STATE=34 TiB+0x7640000, CB0=34 TiB+0x2000000, YSP=34 TiB+0x4000000,
and C_b=CB0+2304b, for b<256. Field elements are four little-endian 64-bit
limbs; points occupy an x field followed by a y field.

|Location|Binding|
|---|---|
|STATE+0,+8|Fixed physical G identity `(c,m)`; G id is `1024*c+m`|
|STATE+88|CHAIN_SET function address|
|STATE+96|Actual hooked yround function address|
|STATE+104|Persistent YSP|
|STATE+112|Actual fin1 function address|
|STATE+168|Explicitly charged coarse-route controller callback|
|STATE+176|Original mining callback; returns status in x0|
|C_b+0..63|Current center C|
|C_b+64i..64i+63, 1<=i<=16|Point [i]S|
|C_b+1088..1151|Point [33]S|
|C_b+1152..2239|yround prefix/denominator scratch, not additional point tables|
|C_b+2240|Next 512-byte frame cursor|
|C_b+2248,+2256,+2264|Absolute k, originating cluster c', next u|
|C_b+2272,+2280,+2288|Group-first k and reverse route words|
|C_b+2296|Original CN pointer `84 TiB+64*(k-1)`|
|CB0-128..CB0-97|Field one, initialized once per batch|
|YSP+352..383,+384|beta and C0=2^32+977|
|YSP+424,+416|CB begin and CB end|
|YSP+1016|0 for mining, 1 for ordinary generation|
|YSP+1024..1616|75-word EB descriptor, populated before ordinary generation|

The CHAIN_SET call is concretely `x0=C_b, x1=CN64`. It must consume the two
stored scalars c,s, set C=[c]P, S=[s]P, and write precisely the point slots
above. It may clobber general registers but must restore its incoming SP and
return normally through LR. yproduce_group reconstructs its registers after
the BLR. ychain_metadata then overwrites the metadata slots, including the
frame cursor. This callback must not redraw coins or read archive records.

The 17 additions in CHAIN_SET have a literal construction: obtain S from the
second scalar multiplication; construct [2]S through [16]S using 15 complete
adds (the first is doubling); double [16]S to [32]S; add S to obtain [33]S.
The two last operations are two more adds. No [17]S through [31]S table is
needed. An S!=O of prime order r cannot make these small multiples infinite.

For emitted index u, set `(rho,q)=divmod(u,33)` and
`delta(u)=33*rho+(0 if q==32 else (16-q//2)*(1-2*(q%2)))`.
The literal round emits C+16S,C-16S,...,C+S,C-S,C, then advances its center by
33S. The source initializer and its post-round center checks agree with this
ordering. This ordering, rather than the simpler but incorrect delta=u, must
be used in every request reconstruction.

## Mining callback contract

The mining callback is a single loop body shared by all G nodes, not a call
into the archive hook. It reads immutable `(c,m)` from STATE and computes
`g=1024*c+m` as a 66-bit integer. The pool is
`16*m+(((c&15)^gf1024_mul((c>>4)&1023,m))&15)`, with polynomial 0x409.
The high two G bits must not be dropped.

Its algorithm is:

1. Draw b_g in [0,r), retain it in the original generator coin state, compute
   R_g=[b_g]Q, and retain R_g for all chains and subsequent requests. The
   initial multiplication is paid by the separate `NODES*GEN_MUL` row.
2. For B=0,256,...,6,079,999,744, populate 256 chain blocks. For each absolute
   chain b=B+j, draw c_b in [0,r) and s_b in [1,r), write their exact 64 bytes
   at `90 TiB+64*b`, set C=R_g+[c_b]P and S=[s_b]P, and construct the same
   17-point step table. The original CNM address is 90 TiB, not 91 TiB.
3. Set the persistent yround mining flag to zero and invoke the actual mining
   yround 992 times on the 256 blocks. For each chain's 32,736 frames, invoke
   the original xemit on A=y^2 at frame bytes 376..407. Metadata is
   `w2=g mod 2^64`, `w3=(g>>64)|(idx<<2)`, `idx=32736*b+u`.
4. Return x0=0 only after all local mining work has completed. The outer
   kind-5 worker then calls actual fin1 and executes DSB before DONE. A scalar
   draw failure returns x0=1 and must enter the run's no-output failure path;
   it must not be treated as successful completion of Phase 0.

Use the callback allocation
`MC=34 TiB+0x7900000`, length64 KiB, in Normal G34 memory. This is above the
current generic service region's end0x7850000. The binding is: b_g at+0..31,
R_g.x at+32..63, R_g.y at+64..95, R_g-infinity flag+96, initialization
flag+104, saved incoming SP/LR+112/+120, absolute chain b+128, batch-local
j+136, mining round+144, point u+152, pool+160, G-low/G-high+168/+176,
pending status+184 and current CNM pointer+192. Bytes+256..4095 hold input
and coefficient scratch; the private stack grows downward from MC+65536,
with at most8192 live bytes including nested field calls. The unused gap
separates stack and persistent state. Whole-Phase0 control reads the status
only after the kind-5 worker's DONE; b_g,R_g and original CNM remain live
through pool reconstruction and ordinary candidate resolution. The current map checks this Normal interval and charges32768 initialization instructions per G. This initializes its64KiB with the existing4-instruction/two-word loop; scratch is reused, not allocated per chain.

For H resolution, let B=2^55-32MiB (its generic base). The16 separate2KiB
PE stacks occupy [B+3*65536+8192,B+3*65536+40960). Generic handler state
uses the start of B+4*65536. Resolve scratch is [B+4*65536+24576,B+4*65536+28672);
its private8KiB stack ends at B+5*65536. The main caller saves SP/LR in
scratch, switches to that stack, then restores both. These intervals are
Normal, disjoint from inbox/handler state and other stacks; the map checks
them. Charge8192/H to initialize scratch/stack and bind its entry. Scratch
operands are written before use; no uninitialized word supplies a coefficient.

## Phase-0 queues, reports and selected-point fetches

Generic ACK certifies that the receiver copied or enqueued the message.
Application work runs explicitly after its input wave and counters complete;
ACK alone never certifies reconstruction, report generation or pool selection.
Phase0 uses actual reqh/fwdq/sendr/fwdr; ordinary groups use fwdq1 in place
of fwdq. Mining fwdq preserves the full66-bit G identity and publishes with
GENERIC_COMMIT1 after constructing the ordinary request descriptor.

At every G, command1 invokes dlvq with request arena[96TiB,4096TiB),
cursor word at96TiB-8 initialized to96TiB and limit4096TiB. Its capacity is125*2^39 requests.
Leave that descriptor bound throughout both phases. Collector G for pool p
is g_p=1024*(p&15)+(p>>4); these are exactly the G identities below2^14.
Their command2 invokes dlvq with a separate1024-ticket arena. Commands0/3
retain direct-copy dlv. A ticket carries the originating H identity and
candidate slot. Rejected tickets still execute and pay their full attempt
service; a collector's1024 cap never replaces the charged attempt population.

Let M=34TiB+0x7a00000. Ticket low/high counters are at M and M+2048;
the1024*64-byte queue is[M+65536,M+131072), with cursor word atM+65536-8 initialized toM+65536.
Reports occupy[M+131072,M+262144): one128-byte report per accepted ticket.
Eight selected128-byte X descriptors occupy[M+262144,M+263168).
The four128-byte outputs of actual psel have separate scratch
[M+263168,M+263680), so fetching X descriptors cannot overwrite a selected
report before its two indices are consumed. These intervals have existing
Normal G34 backing, disjoint from generic state and mining callback storage.

Cold generic_collector_init installs command2, initializes its cursor and both
counter words, and drains its stores. Charge64 instructions per collector for
the body, five64-bit arguments and call/return, plus8 per G to select the
collector identities (conservatively include the sole controller in this
8*(NODES+1) selector census). The ordinary generic initializer still runs and retains
its full price. Report and selection arrays are count-limited: the driver
writes each entry before reading it; unused entries are never consumed.

Mining requests and tickets can coexist because their descriptors, cursor and
counter granules differ. Drain the ticket wave before enumerating its accepted
entries; use (load(cursor_word)-(M+65536))/64, not the wrapping attempt
counter, as the accepted count. Seal that arena through report fetching and selection. For each
accepted ticket, send a command1 report request to H. After the preceding
candidate/reply wave drains, reset all H report queue cursors and counters
before admitting reports; charge64*N_HOME once in the cold census. During
this serialized report phase H main waits on its local queue completion
counter with an unsigned>=1 local monitor wait, verifies cursor==queue_start+64,
snapshots the single
request into private scratch, and
resets the queue/counter BEFORE sending any report reply. It then constructs
the report and sends its four copies. The collector sends no next request
until those copies complete. After its fourth send returns H goes straight
to the next wait, preserving any next request that has already arrived.
A generic control copy can also increment the shared H counter: if the queue
is empty, follow its phase/control flag instead of reading a stale report
request. The cursor/phase guard belongs to the local control allowance. Here the wait
uses LDAXR;CMP;B.HS on the counter and threshold1, otherwise WFE;B to retry,
then CLREX;RET on success. Unlike the inherited exact-total wcnt, this phase
must accept a counter already beyond1 after a control copy. One entry stale
event and the required publication allow at most two failed five-instruction
passes plus the five-instruction success/return path:15 instructions, inside
256. At most1024*N_POOLS report requests plus the fixed phase controls occur;
the counter cannot wrap64 bits in this phase. Charge an additional256 instructions per possible accepted report for local
wait/decode/reset:256*1024*N_POOLS=4294967296 once in the cold census.
This bounds the local wcnt wait, fixed-word parsing, cursor/counter reset and
call setup conservatively; all five message envelopes and their existing32
extra instructions remain. Descriptor construction retains its64-instruction
per-message glue. Copy its128-byte report in four command0 messages to the collector's known
report slot. This is the inherited1024*(5*MSG+32) envelope. Invoke psel on
only the accepted reports, preserving its original caps/failure rule. Return
ok_p by the retained single pool-result message. Clear/reuse an ordinary
request cursor only after that entire request/reply wave drains; never
rebind command1 as a ticket descriptor between phases. Explicitly reset the
G request cursor/counters after the mining candidate wave before selected
fetches, then after selected fetches before ordinary Phase4. Charge64 per G
for each reset in the cold census; the finite reset stores the start and
zeroes both counter words with a DSB after drain.

Each selected-X fetch sends one command1 request directly to the decoded
origin G. In this separate drained phase its five queued words are
(kind,idx,collector-route-low,collector-route-high,destination-local-pointer);
the driver switches its decode mode explicitly before admitting this wave.
Enqueue all8*N_POOLS selected requests before processing any. Their total
population is131072, below each G request capacity even if all target one G.
Freeze each G's accepted count from its cursor distance. Then the coordinator
activates one G main at a time with a command0 phase token and waits for its
command0 DONE before activating the next. A local Normal phase flag occupies
[MC+8192,MC+8224), on its own2KiB monitor granule below the private stack.
The coordinator has the same flag allocation. Only its first64-bit word is
meaningful; the remaining payload is padding and is not read by the main.
Local mains use the inherited
LDAXR/WFE monitor condition; these are ordinary node-local memory stores.
Clear the coordinator's DONE flag before each activation, and initialize each
G flag before this one-time sweep. Main processes its frozen request count,
then sends DONE. PE15 remains available for direct copies even when origin
G equals collector G. Collectors consume results only after this sweep drains.
Charge1024*(NODES+1) for the two256-instruction messages and at most512
local flag/wait/decode/control instructions per activation; add2048*N_POOLS
for collector enqueue/release controls. Both are explicit cold census rows.
The origin executes its charged reconstruction and copies the128-byte result
in four command0 messages directly to the collector at
selected_base+128*j+32*k, j<8,k<4. These replies do not call sendr, whose
ordinary reply route is R. Wait for all four copy completions before consuming
the selected descriptor. Charge the request separately as8*N_POOLS*MSG=33554432;
retain all four reply charges and every REQ_SERVE allowance. The report's
origin/indices remain live in the separate psel scratch through both fetches.
Distribution to each G then uses the inherited four messages per128-byte X
and its existing descriptor storage. Fixed-layout descriptor construction is
inside the64-instruction glue envelope, in addition to the actual sender,
receiver and wait costs. This is an analytic finite phase algorithm composed
from the measured service, not a hidden completing-store callback.

## A concrete modular-r inverse inside the existing allowance

The total envelope `384*200=76800` need not assume 384 long Euclidean
divisions. The following binary extended-GCD algorithm uses at most 511
iterations, each at most 40 A64 instructions in the register layout below.

For canonical 0<a<r initialize `(u,v,x,y)=(a,r,1,0)`. Maintain
`u mod r=a*x mod r`, `v mod r=a*y mod r`, and canonical x,y in [0,r). Repeat until u=0:

```
if u even:   u=u/2;          x=half_r(x)
elif v even:v=v/2;          y=half_r(y)
elif u>=v:  u=(u-v)/2;      x=half_r((x-y) mod r)
else:       v=(v-u)/2;      y=half_r((y-x) mod r)
```

Here `half_r(z)=(z+r*(z&1))/2`, retaining the carry above bit 255. Return y
only if v=1. Zero input is rejected; it must not be passed off as an inverse.
The strictly decreasing integer potential `bit_length(u)+bit_length(v)`
starts at most 512 and ends at 1, proving at most 511 iterations for every
valid full-size input. Gcd invariance and primality give v=1 at termination.

Use u=x0..x3, v=x4..x7, x=x8..x11, y=x12..x15, fixed r=x16..x19, and
scratch x20..x23. Keep the result pointer on a 16-byte stack slot and make no
calls in the loop. SP and LR are preserved; x29 need not be modified.

* Loop test: three ORR plus CBZ =4.
* Parity decisions: at most two TBZ =2.
* Four-limb u-v comparison: SUBS plus three SBCS, then B.LO =5.
* On the v branch, recompute v-u in four instructions; on the u branch use
  the comparison's existing difference. Four EXTR/LSR instructions halve it.
* Coefficient subtraction: four SUBS/SBCS, B.HS, and at most four
  ADDS/ADCS/ADC adding r =9.
* Coefficient half and back edge: TBZ; at most four ADDS/ADCS, CSET of the
  257th bit and B; four EXTR; B to the loop =12.

The longer odd-v branch is `4+2+5+4+4+9+12=40`. The odd-u branch is 36;
the even branches are smaller. Loads, constants, initialization, the final
v=1 check, output stores, status and return fit within 100 instructions.
Thus the all-input dynamic bound is `511*40+100=20540 <76800`.
The static inverse body, with separate u/v half paths and all setup/error
paths, fits within 192 instructions; no inverse lookup table is needed.

The in-package callback checker tests the stated congruence and decreasing
potential on full-size boundary scalars; removing the conditional+r in half
produces a wrong inverse. These finite tests supplement the all-input argument.

## Reconstruction, resolution and output state

The reconstruction layer must recover actual point coefficients, not just
coordinates. A mined X has `X=[a_X]P+[b_g]Q`, where
`a_X=c_b+delta(u)*s_b mod r`; an ordinary Y has
`a_Y=c_k+delta(u)*s_k mod r`. A leaf is
`phi^t(X)+epsilon*Y`, hence coefficients
`(lambda^t*a_X+epsilon*a_Y,lambda^t*b_g)` modulo r.
Per record, precompute {z,lambda*z,-z-lambda*z} for each of a_X,b_g,a_Y
using lambda^2=-1-lambda (three modular multiplications). Leaf/orientation
coefficients then use indexed copies and add/sub, so the24-multiplication
reserve covers both records and solve work. Apply each canonical sign/phi to
both coefficients. Reply256B holds X64,Y64,these three scalars96,flags32;
a mining-selected X reply uses128B for X and its two coefficients. For equal oriented points `(A1,B1),(A2,B2)` with
B1!=B2, the recovered scalar is `(A2-A1)/(B1-B2) mod r`; check [x]P=Q
before output. A zero denominator is a useless tie and is not inverted.

The inherited RESOLVE row pays two candidate envelopes, each containing
the message population, two sets of six reconstructed leaves, 36 comparison
allowances, and two inverse-plus-scalar-verification allowances. Shared
loops over six leaves, 36 pairs and the two orientation candidates are
sufficient. Loops must stop/skip deterministically on an unusable pair;
they must not introduce unlimited retries beyond these charged populations.
The fixed OUTPUT allowance is 2^20 instructions. The exceptional Q=O input
outputs x=0 through the gate's separately charged path.

## Conservative replacement that avoids the exact inherited scalar price

The following complete construction is charged2^20 per scalar multiplication,
including an optional R_g addend and small-index reconstruction. It replaces
the inherited645410 price; scalar_costs.py adds the difference everywhere.

Internal points are Jacobian triples; O is (0,1,0). Affine inputs also carry
an explicit infinity flag. All field outputs are canonical. Keep every
long-lived temporary in a fixed stack slot below offset 4096; calls may
clobber general registers. A multiplication call plus canonicalization has
the following conservative bound: 105 existing instructions +4 argument/BL
instructions +2 LDP +2 C0 materialization +8 RED +2 STP =123<=128.
A square similarly costs 80+3+14=97<=104. A canonical add can use four input
LDP, two C0 instructions, FADD's nine instructions, RED's eight, two output
STP and five call/return/argument instructions: 30<=32. Subtraction is
cheaper. A call with known stack addresses can inline that wrapper. These
bounds explicitly include canonicalization; MULX=109 is not silently
assumed to return a canonical representative.

Complete doubling is:

```
if Z==0 or Y==0: return (0,1,0)
A=X^2; B=Y^2; C=B^2
D=2*((X+B)^2-A-C); E=3*A; F=E^2
X'=F-2*D; Y'=E*(D-X')-8*C; Z'=2*Y*Z
```

This has 2 multiplications, 5 squares, 14 additions/subtractions. Their
bound is 1224; 56 more instructions cover tests, point state and call glue,
giving DBL_NEW=1280. Complete mixed addition is:

```
if affine addend is O: return T
if T.Z==0: return (x,y,1)
ZZ=Z^2; U=x*ZZ; V=y*Z*ZZ
H=U-X; R=2*(V-Y)
if H==0: return double(T) if R==0 else (0,1,0)
HH=H^2; I=4*HH; J=H*I; K=X*I
X'=R^2-J-2*K; Y'=R*(K-X')-2*Y*J; Z'=2*Z*H
```

The ordinary branch has 8 multiplications, 3 squares and 12 small field
operations, costing 1720 before glue. The equal-point branch costs at most
3 multiplications +1 square +3 small operations +128 tests/glue +1280
doubling =1992. MADD_NEW=2048 covers all branches, including O and opposite
points. The shared double is called rather than copied into this branch.

Scalar multiplication starts T=O, and for bits 255 down to 0 performs
`T=double(T)` and, if the bit is one, `T=mixed(T,A)`. It then performs one
optional mixed addition and one normalization. Normalization checks Z=0,
otherwise computes z=1/Z, z2=z*z, x=X*z2, y=Y*z2*z and canonicalizes outputs;
30000 covers the existing inverse, all wrappers, finite/O flag and copying.
The per-bit 64-instruction allowance covers bit selection, spills, loop
control and calls beyond the group-operation envelopes. A 4096-instruction
setup/return allowance and an additional 8192 for small signed delta/index
scalar arithmetic give

```
256*(1280+2048+64)+30000+2048+4096+8192 =912688 <1048576.
```

The small scalar delta has magnitude below 2^22. Compute delta*s modulo r
with at most 22 modular double-and-add steps; even 128 instructions/step,
plus 1024 for decoding/sign and c addition, is below 8192. This includes the
full 66-bit G decode and division of the mining index by constant 32736
(ordinary A64 UDIV on the <2^50 index), and does not use a stored oracle.

The existing CHAIN_ADD=32768 can remain: complete affine addition checks O,
equal/opposite points and uses the usual slope formula, with the doubling
slope 3*x^2/(2*y) on equal points. The envelope
`22549+3*128+2*104+10*32+1024=24485` covers inversion, canonical arithmetic,
point tests and glue, below 32768. Thus the 17-point step construction need
not be doubled merely because scalar multiplication is charged more.

### Initial infinity without an unchecked return code

Define CHAIN_SET to replace **only an infinite initial center** by the fixed
finite base point P in its CB. Keep the original stored coins and nonzero S;
construct [1]S..[16]S,[33]S normally and return through the ordinary ABI.
This is a defined fallback on an excluded event, not a claimed equality
between P and O. It requires no ordinary-driver status check or new retry.
Its copying/tests fit the conservative scalar setup allowance above.

The event is already in the original probability populations: delta(32)=0,
so C=O is one of the round's 33 emitted-point exclusions in
`D_YGEN=YROUNDS*(YPC+2*16+2)/(r-1)` in base/certificate.py line45. For mining,
it is explicitly part of `P_MINE_ZERO_ALL=S_MINE_TOT/r` at line167.
Outside those original bad events the fallback is never used and every
point/coin identity is unchanged. On the event the ordinary affine kernels
still have finite trip counts and all original saturating caps apply.
Original coins must be retained for reconstruction and final [x]P=Q
verification; surrogate CB coordinates are never substituted for them.
The original D_BATCH/g_def populations, enlarged by proof.md6 for pair-inverse
canonicality, and mining DEF_MINE/EXC_WORK must remain charged. This does not authorize shrinking an exception reserve.

### Additional resolve work

Each RESOLVE has two CAND_ITEM envelopes. Charge an additional2^21 per
RESOLVE for coefficient arithmetic and canonicalization. Each envelope
allows24 modular-r multiplications at256 double-and-add steps of128
instructions,36 comparisons/orientation cases at1024,128 further
four-limb operations at128, and65536 parsing/copy/output instructions:
905216 per envelope,1810432 total, below2097152. The inherited four
inverse allowances each exceed the20540 bound above. A further4096 per
REQ_SERVE covers parse/reply glue beyond the retained1024.

## Constructive larger static bound

The existing mapping reserves a1MiB program image. Use a
256KiB reservation for the following *specified shared-loop design*.
This avoids pretending the earlier 64-KiB allocation was certified.
Here a primitive assignment means a four-limb copy/constant/compare,
canonical field or order add/sub, an existing field-kernel call with its
canonicalization, a register/stack scalar operation, a conditional branch,
or a call with fixed-frame arguments. Each expands to at most 64 static
A64 instructions (field multiply/inverse *bodies* are existing shared
code; their new call/normalization sites are at most 32). Fixed-size
schoolbook limb loops are shared helpers, not hidden primitive oracles.

The finite construction and count are:

|Bodies to emit|Construction establishing the static ceiling|Slots|
|---|---|---:|
|Five point helpers|The displayed double, mixed, affine-add, normalize and scalar-driver bodies, each <=64 primitive assignments including their exceptional branches/copies, at <=32 static instructions each. These bodies use field-p operations/copies/calls, not the more expensive general order arithmetic|10240|
|Order/word helpers|Binary inverse192; modular add/sub128 total; modular multiply's 256-bit loop128; small scalar loop128; full index/decode256; GF(1024) multiply loop128; two-try scalar sampler96; bounded schoolbook limb multiply/reduction loop512 (four-limb coefficients and up to32-limb schedule integers). Sum1568, round up for entry alignment|2048|
|X/setup descriptor bodies|Eight coefficient branches of<=8 field assignments, plus a shared<=56-assignment window/EISR/vector body, at32 instructions:3840. Add2048 for descriptor copies, templates and dispatch; round up|8192|
|CHAIN_SET and mining|The explicit 17-add loop, two scalar calls, 256-chain fill loop, 992-yround loop, frame/xemit loop, 66-bit identity/pool calculation, state loads/stores and success/failure exit: <=64 primitive assignments per body for two bodies, at64 each|8192|
|Request and resolve drivers|Two bodies, each <=64 assignments excluding calls into the already-listed arithmetic helpers: decode/load/delta/reconstruct/reply; six-leaf loop,36-pair loop,two-orientation loop,inverse/verify/output|8192|
|Cold/phase control|Descriptor-run loop, read/write+alias pair loop, role-register setup, node-state word loop, tree-barrier loop, group/phase loop, gate/output and finite failure/drain loop: eight bodies, each <=32 assignments at64|16384|
|Constants and descriptors|<=4096 64-bit words: <=128 current mapping records of<=16 words, three75-word EISR descriptors, fixed curve/order/point constants, role-register tables, per-template offsets and shared parameter words. This covers the current service-region splits, not just the36 inherited records|8192|
|Local long-branch expansion/alignment|At most128 new bodies/basic blocks, at most8 instruction slots each|1024|
|Collector/reset callers|Assembled argument/guard/call wrappers30+16|46|
|Total conservative construction||62510|

62510 instruction/data slots is250040 bytes, below the262144-byte
reservation by12104 bytes. The program image including the current assembled-source count is
`460308+262144=722452 <1048576`. This is a constructive algorithmic ceiling,
not a measurement of yet-unemitted code. Every large production count is a
loop bound; no code table is indexed by a node, chain, row or successful
search result. If additional helpers are introduced, this census must be
updated rather than calling the leftover space evidence automatically.

The eight X coefficient families are explicit in base/experiment.py:
P1 slopes `2*y^2` and `2*y` (popax); the pair-forward slope `-(A0+A1)`
and intercept `A0*A1-7*(A0+A1)` (popbat); E `(2*x^3,0)`, F
`(-2*x,-x^4)`, G `(beta^2*x^2,0)` (popk); and BT `beta*2^(43*i)`
(popbt). Each table row is alpha*t*2^shift+gamma modulo p; shifts and
window layouts are the fixed arrays immediately above those functions.
All large row loops are actual existing tbld/tbldk/tbldf bodies. The new
setup computes only their finite initial/step descriptors: S,N,f and
D,D-U,D-W,D-U-W from proof.md and success.md Section 7 Lemma E. Four-limb schoolbook products
and the existing EISR body implement its fixed-dimension vector arithmetic.
Each coefficient branch uses at most eight field assignments (the longest
is the pair-forward intercept or F's two coefficients). The common
56-assignment body exceeds the direct count of window scaling, two residues
and four two-coordinate step vectors. The descriptor copies themselves use
one loop body; the eight families do not each receive a duplicate vector
builder.

The cold/phase functions are also parameterized: emit descriptors by
iterating at most128 current registry records and at most six page-table levels;
write homogeneous runs with the original4-instruction loop and alias pairs
with the original5-instruction loop (proof.md and success.md Section 75.1); emit role registers
from their finite value table; initialize state from a finite offset/value
table; perform the original two-wave tree barrier using indexed child
loops and counted waits; iterate group bounds, calling the actual existing
phase bodies; stop on the existing stop/failure flags. No table descriptor,
tree node or group gets a private copy of executable control code. Hardware
configuration constants and all service-pointer/state initialization must
be part of the package's finite tables and paid boot work.

The constant inventory does not hide an uncharged 22,001-entry schedule in
the executable image. The sole controller can generate m_g by binary search
on integer m in [0,Q1GD], comparing
`m^32*G_GROUPS^17 <= Q1GD^32*g^17`. At most38 comparisons, at most16
fixed-exponent products per comparison, and at most32 limbs per operand
suffice (the largest comparison integer has fewer than1536 bits). The same
bounded schoolbook helper stores a full64-limb result. At a very loose
65536 instructions per product, the whole schedule costs at most
`22001*38*16*65536=876649381888 <2^40` once on the controller. Group values
are consumed sequentially or kept in explicitly allocated Normal data,
not cloned as executable control. Charge2^40 once in transport_map.py's explicit cold census,
inside the retained2^112 cold reserve. This construction makes the fixed schedule
parameters, rather than 22,001 generated answers, the static constants.

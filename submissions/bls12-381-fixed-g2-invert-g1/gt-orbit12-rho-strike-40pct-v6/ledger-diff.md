# Rival versus strike v6: line-by-line ledger reconciliation

Prior work: mpjunior92, `7da5f23a-6633-46b1-9c03-11f6cfedeeec`, public note retrieved separately and live validity=true/officialScore=125.89. Attribution is not a trusted dependency. Both columns use protected C=20000, the same p,r and N=(r-1)/12. The sibling b-v6 repeats the rival's adding-walk ledger, not the square-and-add/Floyd mechanism.

## Numerical rows in the rival public note, in its original order

| Rival row | Rival definition | Strike treatment |
|---|---|---|
| N | (r-1)/12 | identical exact integer |
| S0 | sqrt(2*N*1024/1023*ln(5/3)) | T0=ceil_sqrt(2*ceil(N*1024/1023)*L_hi), L_hi rational upper enclosure of ln(5/3) |
| S | 1.03*S0 | T=ceil(1.10*T0), supported at toy CDF scope by 70000 new paired trials |
| per-step cost | 1+2800/20000=1.14 | one evaluation is 2+8000/20000=2.4; three calls per Floyd iteration give7.2 |
| fruitless cycles | 1+3/(12*1024)+200/1024^2=131129/131072 | no escape/restart multiplier; Floyd charged explicitly, exact coefficient-order argument retained |
| detection window | 1.14*20*2^80 | no parallel-DP detection-window success premise; retained local Floyd cap T |
| setup/table/starts/DP/output | below2^62 units | deterministic ancillary caps below2^100, stress below2^101; auxiliary reporting cannot alter Floyd |
| total | 1.14*S*(131129/131072)+1.14*20*2^80+2^62 | Wexec=(36/5)T+2^100; scored envelope W=8T+2^100 |
| peak | claimed below2^73 bytes, 80-byte DP records | below2^78 bytes, hard capacity2^69 records, explicit128-byte layout and all scratch |

Reconstruction with integer ceilings and the rational log upper endpoint gives rival log2(W) in (125.8845,125.8846). This checks what its ledger says, not its production success or escape derivation.

## Nonreference rows, preserving the rival note's order

| Row | Rival instructions/evaluation | Strike instructions/evaluation | Reason |
|---|---:|---:|---|
| four Fp products |1532|1532|4*383; already present in v5, no new saving |
| one add+8 sub/neg |360|576|nine six-limb bounded carry/subtract/select operations,64 each |
| min of12 keys/index |300|600|six-limb comparisons, selection and index allowance |
| hash/DP/table address |40|128+256|DP predicate/dispatch and twisted jump/label addressing separately charged; full hash only on reports |
| exponent updates |100|384|two four-limb double/add/reduce updates, not two bare additions |
| cycle comparisons/history |60|256|Floyd comparison/control, no invented escape rate |
| loop |30|128|loop and hard meters |
| additional state traffic |not separate|512|ordinary state/scratch traffic explicit |
| itemized sum |2422|4372|every row included exactly once |
| charged ordinary ceiling |2800|8000|inherited analytical ceiling retained, not a new measurement |

Both use the organizer's one-unit dense multiplication convention. Its internal instructions are not charged a second time. The **extra dense square** and **three Floyd evaluations** cannot be absorbed into the nonreference instruction ceiling.

## Where each proposed lever lands

L1: take25% to10% stretch, preselected before the new experiments. The old5/8 comparison with10% stretch gives128.9352..128.9353, versus v5's129.1196..129.1197. Direct paired CDF tests at10% all have positive toy-only lower differences. Mean ratios are separate exact observed diagnostics, not a production confidence bound.

L2: reject the claimed new four-product saving: v5 already has it. Reject1.14 per Floyd iteration. Even the unjustified2800-instruction ceiling transplant, correctly retaining six dense products, gives128.4722..128.4723 at the new T, still losing. Table build occurs once, independent of processor count. Its explicit instruction ceiling is

    1024*(2*510*23000 +12*(72*383+10000)) +2^30

which is below20000*2^50. This overprices setup dense products in instructions; dominant walk dense products retain the protected unit price. The table occupies1024*12*(576+64)=7864320 bytes, before separately bounded scratch.

L3: take the comparison exponent from5/8 down to ln(5/3), retaining comparison success0.40 above the0.39 floor. The new scored envelope is128.7897..128.7898. Replacing ln(5/3) with ln(100/61) saves only about0.02373 bits and consumes the entire probability margin, so it is diagnostic only. Preserve Np=ceil(N*1024/1023); its floor is printed separately. It is a modeling correction, not a derived theorem for square-and-add. Theta=2^-60 does not alter the sufficient Floyd success schedule, and changing report density can only affect ancillary work already far below dominant work.

L4: replace coordinate-bit Hfp by an independent universal linear hash of all4608 canonical point bits. Its73728-byte matrix and report work fit the existing ceilings. The finite false-pair expectation129*M*(M-1)/2^129 and hard2^32 verification cap are explicit. There is no exact Poisson assertion, no loss of scalar-label entropy, and no DP failure debit to core success. This is a soundness improvement, not a dominant saving.

## Floor, repricing and hold

For the retained sufficient mu+ell schedule, strip all nonreference work and all stretch and use the ideal0.39 comparison. Six reference products times sqrt(2*Np*ln(100/61)) still give128.2134..128.2135 bits. This is a schedule-specific optimistic diagnostic, not a universal lower bound on pairing inversion. The verifier separately checks the uniform-map union bound on repeat-by2T at a budget even as large as2^126 is below0.39, so merely using Floyd's earlier exact meeting CDF cannot close this frontier gap under the ideal comparison model either.

Actual strike ceilings: Wexec128.6377..128.6378; Wstress128.6850..128.6851; scored W128.7897..128.7898. Candidate claim128.82 is0.0302..0.0303 above the scored envelope and at least0.1349 above the1.2-repriced executed ceiling. Stress reprices all8000 ordinary nonreference instructions per call, not just the4372 itemized rows, and all ancillary allowances. The separate unused8T scoring envelope is not counted as executed work.

**HOLD, no upload.** A sub125.89 number would require changing the architecture or omitting charged work. The exact integer/Fraction certificate embeds its verifier; no heuristic or production implementation is marked proved.

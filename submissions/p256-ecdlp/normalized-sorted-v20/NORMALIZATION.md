# Normalized-y, sorted-x generator transformation (v20)

This exact v20 candidate is a separate revision. The reviewed/submitted filterfree-direct-v18 is unchanged. The following transformation is the new algorithmic component; finite checks do not by themselves establish its production theorem.

## 1. Exact point and scalar representation

Generate the original independent scalars exactly as in the baseline: X_k=[a_k]P+[b]Q, Y_l=[c_l]P. Preserve b during all X construction. Before hot work, for each nonidentity point A choose epsilon(A)=+1 if canonical Montgomery y_A lies in [1,(p-1)/2], otherwise -1. Odd prime order excludes nonidentity y=0. Replace A by epsilon(A)A. For X store a'_k=epsilon(X_k)a_k mod r and b'_k=epsilon(X_k)b mod r; for Y store c'_l=epsilon(Y_l)c_l mod r. Identity X retains the original abort/identity handling and is never sent to the arithmetic kernel. Normalization takes canonical p-y, not a redundant hot representation.

Sort each node's complete Y table by canonical Montgomery x, moving its scalar alongside its point. Reindex the sorted table consistently. X row indices retain their original order, including the two first-hop phase boundaries. Build compact point and scalar arrays after sorting; no inverse permutation or original-l lookup is needed.

For every original incidence and s in {+1,-1},

X'_k+sY'_l = epsilon(X_k) [ X_k + s epsilon(X_k)epsilon(Y_l) Y_l ].

Thus each incidence's two negation classes are exactly unchanged as a multiset. The full ideal leaf matrix, indexed for analysis by original generators, has the same class multiset. Analyze independence and collision counts on those original generators, not on falsely uniform normalized points. Sorting permutes leaves within each full row; normalization may exchange its two signs. The actual weak raw labels may change with prefix order, but the deterministic U(K) bound in FILTER_FREE.md is independent of representation history and still dominates each label. All good-witness, load and pair-envelope statements therefore need only this unchanged class matrix and U(K), not execution-order independence. Capacity predicates are fixed point-only predicates; sorted-prefix order is not treated as independent randomness.

The event x_X=x_Y and identity X is invariant. Zero-denominator tests remain before inversion. Intruder/false-pair bounds cover all pairs, so reordering equal-key runs cannot increase their envelope. For good canonical witness labels, the adjacent useful-pair lemma applies to any ordering: if every adjacent coefficient pair has zero resolver denominator, transitivity makes every pair useless.

## 2. Hot kernels

Binary-search the sorted Y x array for the first index t with x_Y>=x_X. Run pair_fwd_lo on [0,t), computing x_X-x_Y by a four-limb subtraction; run pair_fwd_hi on [t,q2), computing x_Y-x_X. Empty intervals are allowed. Both keep the zero-difference test. The second interval uses the prior interval's final prefix product as its prefix; there is one inversion for the whole batch, not one per interval. The low/high forward loops differ only in subtraction operand direction. Their per-item upper bound is125; the two calls cost at most125q2+26.

Both canonical Montgomery y inputs are at most(p-1)/2, so their integer sum is below p. The backward y-plus four-limb addition cannot overflow 2^256 and needs no carry materialization, conditional branch or correction. Delete exactly ADC, CBZ and four correction instructions at k_add_store_4. The y-difference absolute-value operation remains unchanged. The backward upper bound becomes639q2+11. No average-case arithmetic cost is assumed. The emitter remains29q2+11. Old full hot price804q2 becomes793q2, with the extra13 forward-call instructions charged separately in the new setup allowance below. Failure prefixes cost no more than these full bounds.

chain_witness.py mechanically derives the old fused body and then these edits. Differential tests compare the new kernels to the canonical reference on normalized sorted input, including empty intervals, random registers/addresses, and set-valued raw labels. Point tests check actual P-256 additions; toy_recovery.py checks row-specific signed scalar resolution on a small a=-3 curve.

## 3. Resolution and uselessness

An emitted normalized plus/minus leaf has coefficients (a'_k+s c'_l,b'_k). Fetch row-specific b'_k rather than the node's global b. The reply remains224 bytes: X, Y, a', b', c'. A collision Z2=epsilon Z1 yields
x=(epsilon A1-A2)/(B2-epsilon B1) mod r,
followed by [x]P=Q verification. All operations and fallback/infinity cases remain charged.

Condition on original generator points and Q. Since a is uniform modulo r, original b is uniform and independent of X. The normalization signs are fixed by these conditioned points, so signed b'_k is uniform too. Across distinct nodes the coefficient pair remains independent. More strongly, if original collision sign is u and the corresponding normalized leaves have signs epsilon_X1 and epsilon_X2, then u'=epsilon_X2*u*epsilon_X1. Both resolver numerator and denominator become epsilon_X2 times the original ones. Vanishing and the recovered scalar are therefore preserved pointwise. This preserves the cross-node1/r uselessness bound. Same-node loads use the original matrix coupling, not independence of row-specific b' within one node.

## 4. Explicit additional cold work

Reserve per node1280q1+3968q2+4352 scalar instructions beyond the v18 ledger with11q1q2 removed.

- Normalization/scalar adjustments:256(q1+q2+1). Four-limb comparison, p-y and two r-scalar negations, fixed selection/stores and loop/indexing fit this bound; coordinate/scalar loads and stores are included. Identity is detected before normalization and follows the charged baseline path.
- Stable bottom-up mergesort of96-byte Y records:28 passes at128 instructions per record, plus at most64 per merge. There are fewer than q2 merges over the whole sort tree. A record comparison reads at most four x limbs from each side; comparing and selecting one side takes at most24 instructions, copying12 words at most24, and pointer/count/loop updates fit the remaining80. Exhausted-run copy paths are cheaper. This bounds3584q2+64q2.
- Compact sorted point/scalar arrays:64q2 instructions, including all copies and loop control.
- Per-X binary search and split-call setup:1024q1. There are at most29 comparisons on2^28 entries. Each reads/compares four limbs, computes a midpoint/address, and updates the half-open interval in at most28 scalar instructions;29*28+13 extra call instructions plus initialization and bound checks is below1024.
- Additional fixed setup:4096 instructions, covering pointer constants and scratch management.

These terms sum to the stated reserve. No memory initialization of unused sort capacity is assumed; every destination record is written before reading. The underlying instruction construction, not a prior accepted score, justifies this reserve. No computation is amortized over another target.

## 5. Memory

Two extra96-byte Y buffers and row-specific32-byte b'_k require at most192q2+32q1<2^36 bytes per node. The previous fixed reservation plus this allocation remains below2^52. Sorting is complete before hot work, but the ledger conservatively retains the scratch allocation throughout. Bulk hop-1/full-bucket/pointer layout is unchanged and stays below2^64 per node. Total memory stays below2^136. Source-service table offsets change, but message sizes and caps do not.

## 6. Claim and evidence limits

The complete total is104863639869181146494799537812703505547264 instructions, with outward rational declaration126301782134220/10^12. The full integrated finite suite accompanies the complete analytical coupling and cost argument. Finite execution is not an independent production success certificate. External exact-version review and user submission instruction remain required; prior v18 approval does not authorize this version.

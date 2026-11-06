# Sorted checkpoint detection and complete accounting

Let q=62, L=ceil(1031 ceil(sqrt(r))/1000), s=ceil(L/q), K=L+s.
For a deterministic trajectory with tail mu and cycle lambda, first repeat
is tau=mu+lambda. If tau<=L, the index a=ceil(mu/s)*s is saved:
mu<=a<=mu+s-1, a<=qs because mu<=L-1<=qs, and X_(a+lambda)=X_a.
This repeat occurs by tau+s-1<K. Binary search finds every equality to a
previously saved state, so detection by K is deterministic whenever tau<=L.
An earlier detected equality terminates sooner. No distributional assumption
on mu or lambda is added. All existing probability calculations for a repeat
by L are unchanged. The new K is smaller than the predecessor's, so its
compression-collision union bound remains valid without enlarging any error.

Maintain 63 slots of 640 bytes each, holding the 576-byte canonical GT
state and two 32-byte coefficients. Store indices/pointers to these slots in
a sorted array. Order keys lexicographically by their 72 machine words in a
fixed word order, each word compared as an unsigned integer. This is a total
order with equality exactly when complete GT states agree. Endianness has no
bearing on that property. Inputs to the order are the canonical Montgomery
tower representation used by the multiplication, not hash/compressed keys.
Each word and each limb is canonical after the pinned reference multiplication. Thus equal
group elements have equal keys. Before every possible insertion search for
equality; if equal, terminate instead. Consequently the stored keys are distinct.
Records are immutable once saved. Moving pointers never changes coefficients.

Use ordinary half-open binary search, lo=0, hi=number_of_saved_keys. At each
iteration m=floor((lo+hi)/2), compare current against record[pointers[m]],
return the record on equality, otherwise set hi=m or lo=m+1. At most 6 key
comparisons suffice for 63 keys. An unsuccessful search also yields the
unique insertion position lo. On insertion shift pointer indices lo..n-1
one place to the right, then put the new record's pointer at lo. This sorting
is by group states alone, not coefficient representations or jump scalars.

For each detected equality the saved index precedes current strictly and their
difference d satisfies 0<d<=K<r. Conditional on jump group values, hash partition
and initial group value, each v_j remains independent uniform by the predecessor's
one-to-one relation between u_j and v_j at fixed M_j. The binary search and
sorting inspect only group elements and counters. Therefore the denominator is
sum n_j v_j, with n_j nonnegative and sum n_j=d. Some n_j is nonzero mod r;
the denominator is uniform mod r even conditional on the detected equality.
Zero causes immediate failure with probability 1/r. Otherwise the invariant
X=h^A z^B proves the recovered exponent and exact G2 output. No search chooses
between relations using their denominators. The success bound 0.39 under the
inherited sole SHA idealization is retained.

## Worst-case scalar cost

Completely unroll a lexicographic comparator over 72 words. Each word uses
LDR current_word, LDR saved_word, CMP, B.ne mismatch, four instructions.
All offsets are below576 bytes and fit scalar immediate addressing. If all
words agree, return equality; otherwise condition flags from CMP determine
unsigned less/greater at the shared mismatch label. A branch does not destroy
those flags. At most288 word instructions execute. Allow96 more for caller BL,
return, base pointers, record pointer lookup, lo/hi/m computation, selection of
the next search interval, loop control and any preserved-register pair. Charge
384 per key comparison, including its binary-search iteration. There are no
runtime word loops, dynamic allocation, arbitrary-width key primitives or
uncharged helper calls. Code for this unrolling is under16KiB. The complete
search is at most6*384=2304 instructions on every state, even adversarial keys
with equal prefixes of71 words. Binary search does not rely on random keys.

Retain the predecessor's outer allowances: 512 instructions for both four-limb
scalar updates/reductions,256 for table index extraction and address computation,
512 for multiword counters and checkpoint saving-index control,512 for the
occasional640-byte record copy,256 for dispatch/equality selection,and1024 for
spills,flags and remaining bookkeeping. Sum2304+512+256+512+512+256+1024=5376,
within the revised 5376 outer allowance.
The pointer-array insertion is charged separately, not hidden in that slack.
At most63 insertions shift at most62 pointers each. Load/store, address,
loop-counter and branch cost at most16 scalar instructions per moved pointer.
All shifts plus pointer writes and insertion setup cost below2^30 instructions
over the entire run. Charge that full bound, independent of outcomes, to the
inherited2^90 fixed-work allowance. This allowance previously covered two
pairings each below2^78, jump table initialization below2^70 and scalar hashes
below2^68; adding even2^30 remains far below2^90. This is total cold work:
insertions occur during the walk and are fully included despite not recurring
at every transition. No amortization over other targets or uncharged preprocessing.

For this revision, count the named dense reference multiplication as one operation under the protected H+R/C convention. Other per-transition scalar work is 3072 compressor +5120 SHA +128 equivalent encoding +5376 outer bookkeeping +384 reference-call interface =14080. The same checkpoint proof applies at q=62, with the revised 5376 bound above. All pointer insertion still costs below 2^30 and is fully paid within the separate 2^90 fixed-work budget. Total W <= K+(14080*K+2^90)/20000. The exact rational work bound is proved in horizon-audit.md; this file introduces no separate score claim.

## Complete scalar coefficient update allowance

Each update adds two canonical four-limb integers below r. Since r<2^255, their sum is below 2^256 and has no lost fifth carry. Load both operands (four LDP total), add with ADDS/three ADCS, load r (two LDP), subtract r with SUBS/three SBCS, choose the reduced or unreduced four limbs with four CSEL and store with two STP. This is 4+4+2+4+4+2=20 scalar instructions excluding the caller's pointers/control. Each result writes directly to its state slot. Reserve 128 for pointer preparation, function control, any callee-saved spills and interface per addition, for fewer than 148. Two updates cost fewer than 296; charge 512. No arbitrary-width unit-cost arithmetic or uncharged inversion is used. Scalar inversion is fixed finish work in 2^90. All other outer categories retain their old allowances.

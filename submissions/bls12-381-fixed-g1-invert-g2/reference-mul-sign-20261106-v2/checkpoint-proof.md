# Sorted checkpoint detection and complete accounting

Let q=1024, L=ceil(1167 ceil(sqrt(r))/1000), s=ceil(L/q), K=L+s.
For a deterministic trajectory with tail mu and cycle lambda, first repeat
is tau=mu+lambda. If tau<=L, the index a=ceil(mu/s)*s is saved:
mu<=a<=mu+s-1, a<=qs because mu<=L-1<=qs, and X_(a+lambda)=X_a.
This repeat occurs by tau+s-1<K. Binary search finds every equality to a
previously saved state, so detection by K is deterministic whenever tau<=L.
An earlier detected equality terminates sooner. No distributional assumption
on mu or lambda is added. All existing probability calculations for a repeat
by L are unchanged. The new K is smaller than the predecessor's, so its
compression-collision union bound remains valid without enlarging any error.

Maintain 1025 slots of 640 bytes each, holding the 576-byte canonical GT
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
return the record on equality, otherwise set hi=m or lo=m+1. At most 11 key
comparisons suffice for 1025 keys. An unsuccessful search also yields the
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
search is at most11*384=4224 instructions on every state, even adversarial keys
with equal prefixes of71 words. Binary search does not rely on random keys.

Retain the predecessor's outer allowances: 1024 instructions for both four-limb
scalar updates/reductions,256 for table index extraction and address computation,
512 for multiword counters and checkpoint saving-index control,512 for the
occasional640-byte record copy,256 for dispatch/equality selection,and1024 for
spills,flags and remaining bookkeeping. Sum4224+1024+256+512+512+256+1024=7808,
within the unchanged8192 outer allowance, with384 further instructions of slack.
The pointer-array insertion is charged separately, not hidden in that slack.
At most1025 insertions shift at most1024 pointers each. Load/store, address,
loop-counter and branch cost at most16 scalar instructions per moved pointer.
All shifts plus pointer writes and insertion setup cost below2^30 instructions
over the entire run. Charge that full bound, independent of outcomes, to the
inherited2^90 fixed-work allowance. This allowance previously covered two
pairings each below2^78, jump table initialization below2^51 and scalar hashes
below2^52; adding even2^30 remains far below2^90. This is total cold work:
insertions occur during the walk and are fully included despite not recurring
at every transition. No amortization over other targets or uncharged preprocessing.

For this revision, count the named dense reference multiplication as one operation under the protected H+R/C convention. Other per-transition scalar work is 4096 compressor +8832 SHA +512 equivalent encoding +8192 outer bookkeeping +384 reference-call interface =22016. The checkpoint algorithm and the 8192 bound above are unchanged. All pointer insertion still costs below 2^30 and is fully paid within the separate 2^90 fixed-work budget. Total W <= K+(22016*K+2^90)/20000. The exact rational work bound is proved in horizon-audit.md; this file introduces no separate score claim.

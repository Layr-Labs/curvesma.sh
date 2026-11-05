# Sorted checkpoint detection and complete accounting

Let q=1024, L=ceil(119 ceil(sqrt(r))/100), s=ceil(L/q), K=L+s.
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
Each word and each limb is normalized by the inherited kernels. Thus equal
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
between relations using their denominators. The success bound 0.395 under the
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

Dense multiplication45312 + compressor4096 + SHA10240 + encoding2048 + outer8192
is69888 scalar instructions. Continue charging70144 for each transition, without
using the separate failed submission's tighter dense/SHA accounting. Add2^90 for
all fixed setup, finish and the complete pointer-insertion work just described.
At most1026*640 bytes plus1025*8 pointer bytes and modest control state are added
to the10737418240-byte jump table. Code,constants,scratch,all records,temporary
output and retained seed/keys fit within2^35 bytes. Local table addresses fit
within2^34, no astronomical addressing or distributed-memory primitive is assumed.

For r equal to the protected production order, R=isqrt(r-1)+1,
L=(119R+99)//100, s=(L+1023)//1024 and K=L+s, the exact public integer inequality

    (K*70144+2^90)^2 < 20000^2 * 2^259

holds. Therefore all charged work divided by fixed C=20000 is below2^(259/2).
The submitted rational score is259/2=129.5. The inequality is an analytic claim
that the reviewer can check, not an asserted local official score. The work
bound applies for all outcomes even if the SHA success heuristic fails.

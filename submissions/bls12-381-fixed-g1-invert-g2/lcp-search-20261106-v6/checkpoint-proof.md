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

## Revised complete search cost

The deterministic detection and conditional coefficient arguments above are unchanged. Replace ordinary repeated full-key comparisons with the exact LCP-assisted search in lcp-proof.md. It proves full search cost1280, new outer cost4352, matrix/slot construction cost below2^30 and unchanged success. Total work is H<=K reference multiplications plus R<=13696*K+2^90 other scalar instructions, W<=(33696*K+2^90)/20000. The finite integer score bound is in horizon-audit.md. All supplied arguments are inside this candidate; no predecessor retrieval is required.

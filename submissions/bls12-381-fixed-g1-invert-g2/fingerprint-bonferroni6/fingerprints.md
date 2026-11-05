# Cache-only refinement: capped fingerprint checkpoint lookup

This is analytic claim evidence; no production execution or official qualification is asserted.
This builds on the ref5 direct injective norm-one SHA512 walk and the Bonferroni
horizon refinement. It preserves the full supplied-P reduction and0.39 threshold.

Take each state's fingerprint from the first128 bits of its SHA512 digest,and
the24-bit jump label from its last24 bits. The fields are disjoint. Hash the
injective norm-one encoding exactly as ref5. Store checkpoint records as
(fingerprint128,fullGT576bytes,A32,B32),656 bytes. Keep a pointer array sorted
by128-bit fingerprints,allowing duplicates. Lower-bound binary search locates
the first matching fingerprint;scan equal-fingerprint entries and compare full
canonical group states. Only full equality supplies a relation. Insert at the
lower-bound position when saving. This does not assume injective fingerprints.

Reorder transition bookkeeping:at state X_t,compute its digest,fingerprint and
next label,then search/save as appropriate,and finally advance. Initial state's
fingerprint and the final state's hash require at most two extra hashes,charged
to fixed work. Exactly K group multiplications and at most K+2 hashes suffice.
The deterministic equally spaced anchor proof and conditional denominator
uniformity are unchanged;only group-derived fingerprints affect decisions.

Keep a counter F for false full comparisons (matching fingerprint,unequal full
group). After each false comparison increment F;if F>2^40 fail. True equality
is processed first and does not increment F. Hence all outcomes execute at
most2^40+1 full probes. Charge1024 scalar instructions per probe including
72-word equality,counter,duplicate enumeration and dispatch. Their total is
below2^51 and fits the existing2^90 fixed-work budget,with ample setup slack.

The overflow loss is rigorously bounded inside the same SHA512 idealization.
Condition on initialization and the entire24-bit partition-label function.
The complete virtual group trajectory to K and its fixed saving schedule are
then determined,independently of the128-bit fingerprint function. For any time
and previously saved checkpoint with unequal group values,the conditional
fingerprint-match probability is exactly2^-128. There are at most(K+1)*1025
such time/checkpoint pairs. Even allowing repeated group states and duplicate
checkpoint fingerprints,linearity of expectation gives

    E F <= (K+1)*1025/2^128 < 2^12,

using K+1<2^129 and1025<2^11. Markov yields

    Pr(F>2^40) < 2^-28.

This is not an unjustified independence assumption about an adapted walk:
fingerprint bits and label bits are independent subfunctions in the disclosed
idealized digest model. The uncapped virtual trajectory provides an upper
bound before the actual cap stops execution. No fingerprint collision can
produce an unchecked wrong output. The single computational heuristic remains
transferring ideal SHA512 behavior to the actual fixed construction.

Worst-case per-step fingerprint search:11 comparisons. Each iteration needs
at most32 instructions for two64-bit key loads,record pointer,midpoint,compares,
lo/hi update and control. Add128 setup/final-key-check instructions:480<=512.
All equal-key enumeration/full comparisons are paid through the cap above.
Retain1024 scalar-update,256 jump-index/address,512 counter/save control,512
occasional656-byte copy,256 dispatch,and1024 spill/flag allowances. Together
3584+512=4096 per step,replacing the earlier8192 allowance. Copying656 bytes
uses82 word loads/stores and fits512 with loop/pointer work. Insertion pointer
shifts remain below2^30 fixed instructions. Fingerprint extraction/storage is
included in the hash preparation and copying allowances. Memory adds only
16 bytes per checkpoint,still below2^35 total with the10GiB jump table.

Total direct-hash transition:one protected reference multiplication plus
4096 normalization+2048 sign+15872 SHA512+4096 outer instructions =1+26112/20000.
Use fixed2^90 instructions including capped full comparisons,insertions,extra
hashes,all construction,failures,checks and output. At Bonferroni horizon
L=ceil(1070ceil(sqrt(r))/1000),subtract2^-28 in addition to earlier tiny losses.
The checked lower probability0.390403... still exceeds0.39. The exact inequality (K*46112+2^90)^50<20000^50*2^6437 supports the declared rational bound6437/50=128.74. This is a claimed bound,not an official score.

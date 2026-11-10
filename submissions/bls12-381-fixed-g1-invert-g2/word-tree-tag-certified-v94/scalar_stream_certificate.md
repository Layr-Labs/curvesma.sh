> Current payload/confirmation and COMPLETE score contracts are in proof.md and nh_certificate.md. The added expensive-verification allowance is1L;aggregate base5790/network946/hot6736/retained6816. Any historical full-S storage or cold-only phase language below is superseded: S192 is temporary,stored Ctag48 is probabilistic,and phase attempts use the new hot cap.

# Complete scalar payloads with separate baby and giant bounds

Only the scalar-stream implementation/accounting changes. Every payload still contains the FULL32-byte canonical Fr scalar,metadata keeps ALL128 ordinal bits,and all GT sequences/selected phases/full parameters remain unchanged. This uses the proven public baby-index range,not smaller curve/order parameters.

## Baby j stream

The baby points are actual h^j for j1..m,where m=44951011978550636487095639346140676096<2^126<r. Maintain an ordinary128-bit public j counter initialized1. In the BABY prepare-record path,call the actual baby_scalar_fill symbol in metadata_io.S exactly once per emitted record. It loads both j limbs,stores them into the lower16 bytes of the32-byte scalar field,stores zero into upper16,then increments the two-limb state by ADD/ADC and stores it. Output and state are disjoint;no partial/identical alias is allowed for this mutating-state API.

The actual source has7 instructions throughRET:one LDP,three STP (two output,one state),ADD,ADC,RET. Two argument MOVs and BL give10;retain12. All32 output bytes are initialized. Because j<r,the result is EXACTLY the canonical full Fr encoding of j;there is no modular approximation or field-size truncation. At the last point j=m it emits m and stores m+1<2^128 safely. Any low64 carry is propagated. Callee19..29/SP/LR and x18 are unchanged.

The old generic baby Fr-add/state-copy path is replaced,not additionally executed. Its32-byte state must not be read after replacing it by a16-byte counter. The helper writes the payload scalar directly. Metadata's separate128-bit ordinal counter still starts0 and advances once per record under the already charged main704,so the induction ordinal=j-1 holds. Both advance only during preparing a NEW source batch,not when returning its payload/meta frames. The pending flag still prevents double generation. Cold initial trace records j1/j2 use the same scalar preparation and are conservatively hot charged.

## Giant bi stream and complete unequal worker prices

Giants retain the full255-bit Fr state bi and the actual fr255_add.S canonical32-byte addition. It has21 instructions throughRET,four argument MOVs plusBL give26;retain28. All full-r carries/reduction and safe output/input alias are unchanged,since a+b<2r<2^256 and the source's trial subtract/select is exact. It is not priced as a128-bit addition. Main704 still pays payload/state copy and ordinal/control work. The final unused increment is also conservatively charged.

There are at most m baby records and at most m giant records in EVERY execution,including early identities or aborts. Hence total scalar work<=12m+28m=40m=20L. This is an exact aggregate bound,not an assumption that a failed prefix happens to have balanced colors. If a run stops early,neither count can exceed m,so the bound still holds. Full successful streams have exactly both m. No rejected/failed path is omitted.

The scalar remains j or bi in the unchanged collision formulas a=bi+lambda*j or a=(b_t-lambda*b_s)/(1-lambda). No event count,key/fingerprint,candidate phase,probability distribution or output verification changes. The actual curve p381/order r255 and target uniform GT are preserved. Only a known public loop counter has a cheaper implementation than a generic full-field add.

## Shared scoped tests and new complete bound

The shared-only audit interprets the actual7-instruction baby helper for j=1,m,2^64-1,2^64 (then bounded variants if more inputs are requested),checking ALL32 output bytes,canonicality,exact next-state carry,disjoint writes and ABI. It also newly interprets the actual21-instruction Fr addition on full-r zero/wrap/max/random cases,including exact input alias. The independent integer oracle uses the full protected r. These are scoped helper checks,not a production inversion,probability proof or qualification. Candidate code was never executed/compiled on the host,and no new pass is asserted before a genuine shared report.

The unchanged other-generation price is4664L. Adding the proven scalar20L gives generation4684L. Two count120/scatter240 passes720,scan128,cheap_confirmation1,expensive_verification1 and spare256 give base5790. The current tag112 transfer proof gives network946. Hot6736 is retained6816 (80 spare). Full N=6816L+2^124 satisfies EXACT N^1000<20000^1000*2^124527,proposed124527/1000=124.527. This also tightens four unused top-allowance instructions;the explicit inner spare256 is unchanged. Cold128-bit initialization and tiny new code/state fit the existing2^124 cold/2^138 peak/advice0 bounds.

Full p/r,target distribution,success0.39,C20000,sixtracks/onecomputationmetric remain unchanged;the new two-cap proof gives0.39002671387253846. All source five-frame byte assembly,metadata40/banks/READ5/two80-bit queues and odd-tail copy charges are retained. No new score,acceptance,promotion or physical hardware qualification is claimed before actual shared review.

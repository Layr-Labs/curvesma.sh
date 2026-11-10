# Complete scalar payloads with separate baby and giant bounds

Only the scalar-stream implementation/accounting changes. Every payload still contains the FULL32-byte canonical Fr scalar,metadata keeps ALL128 ordinal bits,and all GT sequences/selected phases/full parameters remain unchanged. This uses the proven public baby-index range,not smaller curve/order parameters.

## Baby j stream

The baby points are actual h^j for j1..m,where m=32196482319508782495259497990721634304<2^126<r. Maintain an ordinary128-bit public j counter initialized1. In the BABY prepare-record path,call the actual baby_scalar_fill symbol in metadata_io.S exactly once per emitted record. It loads both j limbs,stores them into the lower16 bytes of the32-byte scalar field,stores zero into upper16,then increments the two-limb state by ADD/ADC and stores it. Output and state are disjoint;no partial/identical alias is allowed for this mutating-state API.

The actual source has7 instructions throughRET:one LDP,three STP (two output,one state),ADD,ADC,RET. Two argument MOVs and BL give10;retain12. All32 output bytes are initialized. Because j<r,the result is EXACTLY the canonical full Fr encoding of j;there is no modular approximation or field-size truncation. At the last point j=m it emits m and stores m+1<2^128 safely. Any low64 carry is propagated. Callee19..29/SP/LR and x18 are unchanged.

The old generic baby Fr-add/state-copy path is replaced,not additionally executed. Its32-byte state must not be read after replacing it by a16-byte counter. The helper writes the payload scalar directly. Metadata's separate128-bit ordinal counter still starts0 and advances once per record under the already charged main712,so the induction ordinal=j-1 holds. Both advance only during preparing a NEW source batch,not when returning its payload/meta frames. The pending flag still prevents double generation. Cold initial trace records j1/j2 use the same scalar preparation and are conservatively hot charged.

## Giant bi stream and complete unequal worker prices

Giants retain the full255-bit Fr state bi and the actual fr255_add.S canonical32-byte addition. It has21 instructions throughRET,four argument MOVs plusBL give26;retain28. All full-r carries/reduction and safe output/input alias are unchanged,since a+b<2r<2^256 and the source's trial subtract/select is exact. It is not priced as a128-bit addition. Main712 still pays payload/state copy and ordinal/control work. The final unused increment is also conservatively charged.

There are at most m baby records and at most n giant records in EVERY execution,including early identities or aborts. Hence total scalar work<=12m+28n. This is an exact aggregate bound,not an assumption that a failed prefix happens to have balanced colors. If a run stops early,baby count cannotexceedm andgiant countcannotexceedn,so the bound still holds. Full successful streams have exactlym babies andn giants. No rejected/failed path is omitted.

The scalar remains j or bi in the unchanged collision formulas a=bi+lambda*j or a=(b_t-lambda*b_s)/(1-lambda). Collisionformulas,key/fingerprint,candidatephases,inputdistribution andoutputverification areunchanged;eventcounts are separatelyrecomputed at thecurrentm/n inproof.md. The actual curve p381/order r255 and target uniform GT are preserved. Only a known public loop counter has a cheaper implementation than a generic full-field add.

## Shared scoped tests and new complete bound

The shared-only audit interprets the actual7-instruction baby helper for j=1,m,2^64-1,2^64 (then bounded variants if more inputs are requested),checking ALL32 output bytes,canonicality,exact next-state carry,disjoint writes and ABI. It also newly interprets the actual21-instruction Fr addition on full-r zero/wrap/max/random cases,including exact input alias. The independent integer oracle uses the full protected r. These are scoped helper checks,not a production inversion,probability proof or qualification. Candidate code was never executed/compiled on the host,and no new pass is asserted before a genuine shared report.

Current scalarbound12m+28n replaces theold balanced20L quantity. Itincludes bothfull-width giantFr work andlastunusedincrements. Publicgiantloopcountn isseparatefrombabycountm;both128-bit andmultiples512/2048. Non-scalar generation4309L=(2526+721+350+712)L. Actual scalar allowance12m+28n charges atmostm baby andn giant records onEVERY success/failure path;no balanced20L average isused. Local720L+scan128L+confirmation1L+inner spare256L,network920L andouter spare80L give N=6414L+12m+28n+2^124. EXACT N^625<20000^625*2^77736 gives124.3776. Bothfull bucket constructions,cold2^124/peak2^138/advice0,referenceC20000/threshold0.39/fullfields/sixtracks/solemetric remainunchanged.

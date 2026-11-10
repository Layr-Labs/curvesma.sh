> Current unequal nb/ng counts,ordinalmap,two verification caps and complete probability/work are in proof.md and scalar_stream_certificate.md. Native/table/metadata/sort constructions below retain their scoped contracts. FullS192 is temporary;stored Ctag24 requires verified continuation.

# Eleven-digit constant-s table and exact small-quotient correction

The odd Frobenius still multiplies by actual s from xi^((p-1)/2)=(s,s),2s²=-1. All inputs are canonical full381-bit Montgomery residues, Rm=2^384. Partition ALL381 significant bits into widths [35,35,35,35,35,35,35,34,34,34,34], starts [0,35,70,105,140,175,210,245,279,313,347]. For each digit generate EVERY table[d]=d*sNative*Rm^-1*2^start modp by the actual pinned seed computation and full fill_lut.S canonical addition loop. This is constructed RAM, not an oracle. Its64-bit loop counter supports n=2^35/34 without truncation. All table construction, zero entries, seeds and controls are paid cold; no target advice or supplied auxiliary powers.

The fixed-s table is STILL48*(7*2^35+4*2^34)=13.5TiB,nowat64..77.5TiBbyexistingbasepointerargument. NewrecurrenceA/Bvector tablesoccupy8..62TiB,parent4/private<1/children80..128remain. All47bit/disjoint. Tablecontents/fill/sSource/helper212countsareunchanged. Bothworkers'actualvector/fixedfillsremain<2^50cold,globalcold124/peak138/advice0. See fp2_certificate.mdforfullnewtablegeometry;old eightcolumn32-bit table geometry is replaced.

## Native digit and accumulator proof

fixed_s_lookup.S loads all six input limbs into x19..24, so exact output/input alias is safe. Each noncrossing digit uses one UBFX; five crossing digits use EXTR followed by a35/34-bit AND. No digit is lost, overlaps or reads beyond the six canonical limbs. Use stride48*2^35=3*2^39, materialized by ONE MOVZ x27,#0x180,lsl32. Width35 advances by this stride; width34 advances by its LSR1 in the same ADD instruction. The final digit has no advance. MADD with entry size48 accesses only its initialized digit table and all pointers remain47-bit.

Eleven steps cost11*12+5 crossing instructions+10 advances=147. Each step performs one digit extraction (or two), one MADD, three LDP and seven full-carry accumulator additions. All eleven table values are canonical in[0,p), so their exact nonnegative SEVEN-limb sum Y is below11p and fits448 bits. There is no384-bit accumulator wrap.

## Fully constructed fifty-entry quotient table

Use a3728-byte constant block, disjoint from all outputs, tables and frames. Its first520 bytes retain the unchanged recurrence-map contract (ordinary p at448,D=floor(p/2^325) at512). Pad520..527 with zeros. Construct fifty64-byte slots at528+64*q for q=0..49. Each contains all SEVEN ordinary64-bit limbs of q*p plus an initialized zero eighth word. Construct slot0 as zero; construct subsequent slots by exact seven-limb addition of ordinary p, without modular reduction. All fifty slots and all padding are filled by actual charged integer additions/stores in this attack. They are not native Montgomery field values, searched advice or supplied powers.

By quotient_certificate.md, qhat=floor((Y>>325)/D) equals floor(Y/p) or one more. Here Y<11p implies floor(Y/p)<=10, hence0<=qhat<=11. The table product read is therefore always initialized and within the fifty slots, and is EXACTLY the same seven-limb qhat*p as the prior scalar multiplication path.

Reduction costs32: LDR D1,EXTR1,UDIV1,ADD tablebase1,ADD shifted index1,three LDP plus LDR top4,seven SUBS/SBCS7,three LDP p3,six original-borrow CSEL6,seven addback instructions7. If the product overshoots, the difference is in(-p,0), so adding p under that ORIGINAL borrow gives the canonical residue. Otherwise no addition occurs and the difference already lies[0,p). Loads preserve the borrow flags; every CSEL preserves them until all six p limbs have been selected. Top output limb is zero. No native x18 or hidden arbitrary-precision operation is used.

## Complete counts, ABI, key and evidence scope

The212-instruction source comprises prologue/setup13,input loads3,zero accumulator7,digit loops147,reduction32,stores3,epilogue7. Frame96 restores19..29/SP/LR, no18 and no BL. x4 retains output and x25 the constant pointer; product scratch12..17/top26 does not overlap live accumulators5..11. All field input reads and every table read precede output writes. Exact alias or disjoint output/input is supported; partial overlap is forbidden.

The unchanged odd helper costs572 throughRET. Complete key path accounting is now111 own+2*65 sign+572 odd=813, caller8 gives821 retained824, as proved in key_path_certificate.md. Its shared3728-byte constant argument, exact field product, raw-history separation, identity handling and all12-fiber proofs are unchanged.

Complete generation=2720+824+376+22+704=4646; two radix passes720+scan128+cheap_confirmation1+expensive_verification1+spare256 gives base5752; network933 gives hot6685, retained6760 (75 spare). N=6760L+2^124 satisfies EXACT N^125<20000^125*2^15562, giving proposed15562/125=124.496. Full p/r/target distribution/success0.39/C20000, six tracks and one computation metric remain unchanged. The self-contained NH2^-63 composition gives the same success bound0.3900120659047603. No observed new score, native pass, qualification or promotion is claimed. Shared llm-verifier alone executes/reviews/scores.

The shared-only audit additionally interprets this actual212-instruction source against independent full-p multiplication for the actual fixed s, with boundary/random inputs, exact input alias, initialized reads/writes, ABI and count checks. Its virtual readonly table contents are scoped test fixtures, not a free production table algorithm. It does not certify production fill execution, GT success probabilities, hardware, score or qualification. No candidate code was host executed or compiled.

The shared constant block is now3728 bytes with fifty slots q0..49 for the recurrence map; fixed-s still reads only q0..11 and retains exactly212 instructions. The full table-construction and QREDUCE31 proof are in quotient_certificate.md.
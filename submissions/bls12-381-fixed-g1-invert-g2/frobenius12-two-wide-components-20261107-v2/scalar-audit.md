# Scalar work for the12-orbit implementation

The protected unit remains one dense Fp12 multiplication,C=20000. Each transition uses exactly three such reference operations;their internal instructions are not counted again. All other work is charged in scalar AArch64 instruction equivalents,including field operations,loads/stores,addressing,bank selection,control,hashing,lookup,failed paths andoutput. No SIMD,crypto extension,free oracle orparallel wall-clock discount is used. The following is analytic implementation accounting,not a measured trace orproduction run.

## Native field primitives

Use the pinned blst scalar native Montgomery Fp arithmetic,with constants preconverted cold. The organizer's pinned measurement reports383 instructions per blst_fp_mul. The multiplication path is constant-time scalar arithmetic;allow416 including up to33 caller/pointer/preservation instructions. Arguments are pointers to aligned fixed48-byte buffers;the three argument moves andBL are within that allowance. Live higher-level pointers are preserved by normal ABI orstored in fixed scratch;there is no arbitrary allocation per multiplication. This is an upper instruction bound for a concrete kernel,not a candidate-selected normalization price.

A six-limb modular add can use six operand LDP,six add/carry instructions,three modulus LDP,six candidate-subtract instructions,six CSEL andthree STP:30,plus bounded interface. Subtraction/negation use similarly bounded borrow/correction code. Charge48 for every Fp add,subtract ornegate,including input/output traffic andcaller work. All outputs are canonical. A general constant Fp2 product could use three Fp multiplications,but no such general constant is needed online here. At p=19 mod24,all c[k,d]=(1+u)^(d(p^k-1)/6) have one zero coefficient orequal/opposite coefficients,because(1+u)^4=-4. Real/pure-imaginary products need two coordinate scales plus optional swap/negation. For c=(s,s),c*(x+yu)=s(x-y)+s(x+y)u;for c=(s,-s),it is s(x+y)+s(y-x)u. These use two Fp multiplications andat mostfour charged Fp add/sub/neg operations including any input conjugation. Cold construction supplies s andclass bits,not free advice. The formulas are exact in Fp2 anddo not change any group/hash values.

Pinned reference/canonical sources:
 scripts/reference_costs/measurement.json and tracks/docs/COST_MODEL.md;
 https://github.com/supranational/blst/blob/f7a88067ed98385ebc6cf41ccecf21ead706a5b7/src/fields.h
 https://github.com/supranational/blst/blob/f7a88067ed98385ebc6cf41ccecf21ead706a5b7/src/asm/mul_mont_384-armv8.pl
The reduction ends with conditional modulus subtraction. The tower basis is the same pinned Fp2/Fp6/Fp12 layout used by the organizer's dense reference operation. Cold conversion is included;no missing per-step canonicalization is assumed.

## Canonical orbit helper:8768 instructions

The helper is a finite straight-line/12-candidate routine using preallocated aligned scratch andthe primitives above. Its exact algorithm is in proof.md. Count every worst-case branch,without average-coordinate assumptions:

1. Selected-pair orbit construction:one structured Fp2 constant product y=T(x),two Fp multiplications;beta*x andbeta*y usefour more. The twelve candidates follow from those values byexact additions/negations. No repeated full-point Frobenius orfield inverse is used.
2. Canonical b output:reuse selected candidate pair;transform the other two b pairs byprecomputed constants,at most four Fp multiplications.
3. Canonical a-sign:find the first nonzero Fp2 slot;slot zero only conjugates,otherwise one constant Fp2 product,at most two Fp multiplications. No later a coordinate is needed to determine its lexicographic sign.

Total at most12 Fp multiplications,charged12*416=4992. There are at mostfour structured Fp2 products,so retaining16 Fp add/sub/neg operations is conservative,charged768. Candidate rotation generation uses16 more Fp operations;up tofour input conjugations addfour,charged20*48=960. These include all canonical reductions andstores through the primitive bounds.

Selecting the minimum twelve96-byte pairs needs at most11 comparisons. Compare real then imaginary residues by12-limb subtraction in the corresponding little-limb order:12 LDP across both operands plus12 SUBS/SBCS,one carry decision;copying a selected96-byte pair uses six LDP/six STP. Including pointers,orientation update andloop control,64 instructions per comparison suffices. Charge768 for this entire selection/copy stage. The images are distinct,so no unbounded tie handling exists.

Charge768 more for key/sign/output control:testing b0/b2 zeros,finding the first nonzero a slot,constructing a sign-only scratch tuple,existing89-instruction sign routine,73-instruction b word-copy,reordering the36 words,storing k andstatus,initializing fixed pointer slots andall output selection. Testing one96-byte pair costs six LDP,11 ORR andbranch;two such tests remain small. A sign-only a scratch can be cleared with18 STP andreceive the first nonzero transformed pair withsix LDP/STP;the sign depends on no later coordinate. All output bytes are therefore explicitly constructed.

Finally charge512 for Fp2 wrapper pointer arithmetic,frames,caller preservation andtemporary selection not already in the field primitives. Four structured Fp2 products can keep their four base pointers in preserved registers with fixed entry/exit traffic,orload them from scratch before calls. The two multiply/up tofour add interfaces per product need at most128 instructions of additional pointer/control each,including bounded copies. Constants are indexed from the fixed72-entry table. There is no hidden variable-length loop orordinary allocation.

Sum4992+768+960+768+768+512=8768. This deliberately overcounts several operations also in the retained baseline below. It covers all12 field multiplications even when identity/sparse coordinates would use fewer. The identity-abort branch is charged andits probability loss is in probability.md.

## Two-component transition andcomplete ledger

Retain the full12032 ordinary instructions per transition from the three-component analytic certificate. It consists of1968 baseline+8768 canonicalization+1280 extra=12016,with16 further slack. Baseline complete path58 main+46 scalar+37 key+820 tabulation+71 sign+797 full-output BLAKE3+128 ordinary prefix=1957;globally capped confirmations add<=K,so1958K is covered by1968K. Both orbit hashes require their full outputs andretain797 each. The earlier one-factor controller is an accounting template,not the entire current two-factor algorithm.

Current extras fit1280:second hash797+CVcopy4+arguments3=804;two72-bit index extractions andstores<=12;phase multiplication20;one extra reference interface7;one extra scalar helper/interface50;two36-bit-bank addresses/selection protocols<=48;phase-minus-k,metadata,canonical interfaces andall residual context/control<=339. Sum804+12+20+7+50+48+339=1280. This explicitly charges the high8 indexbits andall selector arithmetic. No cold scalar-counter operation occurs online. Keeping canonical helper's older89/73 sign/key allowances overcounts current71/37,which is conservative. No factor's internal reference instructions are counted again.

Exactly two reference dense Fp12 operations per step give2K reference units. All ordinary work<=12032K andcold<=2^120;W<=(52032K+2^120)/20000. C=20000 is unchanged. No parallelwalltime,free transfer,oracle orwideinteger unitcost is used.

## Cold work,storage andencoding

At most4w+1<2^75 scalars,1024 attempts each,fit<2^85. Standard64-byte SCAL messages use128-bit big-endiancounter16,not theoverflowing old64-bitcounter8. Generic standard BLAKE3 output32 is fullycomputed. Each attempt<2^20 gives<2^105. At most2^73 base jumps*1024products*2^20 scalarinstructions give<2^103. All12variants/coefficientgeneration,bankselectors,writes,configuration andeven fullcapacityzeroing stay<2^105 or2^90 respectively;fixedpairings/checks/finishes/copies are smaller. Their sum<2^120 includesalltrials/failures. Counter<2^85 fits128 bits;72-bit onlineindices use two64-bit limbs and36-bitselectors,with explicitcharged highbits.

field-ops.S add/sub36/37 plusinterface<48 areunchanged. Scalar-helper46 usesbaseADD#576 andpairedlocaloffset0..48;tabulation820 usespairedloads withtworebasingADDs tokeepimmediatesencodable. These wereindependentlyread-checked,notcompiledorexecutedlocally. storage.md statescompletecapacity2^86+shared<2^87,includingpadding/descriptors. The finitebankedRAM realization isamaterialreviewcondition,not actualinstalledhardware. Allcode isstaticevidence;sharedllm-verifier owns execution/review/scoring.

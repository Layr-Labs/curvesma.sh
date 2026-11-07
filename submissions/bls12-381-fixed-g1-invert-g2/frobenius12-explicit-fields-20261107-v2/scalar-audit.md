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

## Retained baseline andextra controller

Conservatively retain the entire2096-instruction baseline of the local complete-leaf package. It includes one optimized standard BLAKE3 final kernel,tabulation,prefix search andall capped full confirmations,two scalar updates for one factor,ordinary main control,key-copy/sign andinterfaces. The baseline complete path certificate is62 main+46 scalar+73 key+892 tabulation+89 sign+797 hash+128 ordinary lookup=2087,with confirmation aggregate<=K,so2088K<2096K. Several key/sign operations are now counted again in the canonical helper;keeping them is conservative. The native scalar,key,hash,tabulation andprefix algorithms/contracts are unchanged except their call order/context.

Charge1280 additional instructions per transition for everything introduced beyond that retained baseline:
- The second standard BLAKE3 final hash797,prefix CV copy4 andargument/call3 cost804. U is computed once andthe40-byte block reused unchanged. Both kernels retain all7 rounds andfull input/output/interface accounting.
- J index extraction uses three LDR/AND/STR triples,9. PHAS high-word multiplication by12 uses four LDR,eight MUL/UMULH,six carry instructions,one constant andone store,20. All arithmetic is four64-bit limbs;no256-bit unit primitive is used.
- Two additional reference call/interfaces use at most14,andtwo additional scalar helpers/interfaces at most100. Coefficients are already rotated in their cold variant records,so no per-step Fr multiplication is hidden.
- Three bank/record addresses andselector protocols use at most24 instructions each,72 total,as spelled out in storage.md. Their memory data reads are already in the three reference kernels/scalar helpers;the address/selector is not a free operation.
- Phase minus applied k modulo12,metadata/canonical-helper arguments,identity/found flags,extra shared-context selection andall residual control fit the remaining261 instructions after these1019 listed instructions. Fixed stack/scratch allocation is cold. Persistent raw state andcoefficient registers survive ordinary ABI calls;the native hash explicitly restores the callee registers it reuses.

Thus ordinary nonreference<=2096+8768+1280=12144 per step,below the declared12288 with144 further slack. The globally capped confirmations are already in2096K. Both hash output extraction andall table operations are included. A terminal last test andall scheduled cold paths are in the separate2^120,not omitted.

## Complete cold work andmemory

Bounded sampler attempts<2^60,each<2^20 instructions,give<2^80. Base jumps need<1024 dense products each at a deliberately loose<2^20 scalar instructions,for<2^80. Variant Frobenius/scalar generation,all data stores,bank selection/setup,field constants andgeneric cold hashes remain below2^90. Input validation,two pinned pairings,subgroup tests,final inversions/exponentiation,G2 output andverification remain below2^90 by the source's conservative fixed-arithmetic arguments. Checkpoint copies/insertions,metadata andframes are far smaller. Their full sum is below2^120;no stored result ordevice bank contains uncharged searched advice.

Memory includes full bank capacities,padding,shared arrays,code,advice,scratch,frames andall processors' storage. storage.md reports3072 banks of at most2^50 bytes plus shared tables<2^45 bytes,well below2^63. Code/algorithm design is not used to hide precomputed results. The bank-selection mechanism only chooses RAM storage;it supplies no group operation,random oracle orarithmetic shortcut. Physical availability/practicality is not claimed. The finite instruction-equivalent memory realization is material to review separately from H.

Total at most3K reference operations plus(12288K+2^120)/20000 other work. The exact rational inequality in probability.md supports the submitted bound only under this complete implementation/accounting andH. No scientific qualification orscore is invented.

The explicit field leaves include every operand/modulus load,carry,conditional correction,store andRET. Add uses five pointer/constant setup,six operand LDP,six additions,three modulus LDP,six subtractions,six CSEL,three STP andRET=36. Sub uses setup5,operand LDP6,subtractions6,CSETM1,modulus LDP3,AND6,correction-add6,STP3,RET1=37. Both use x0..x16 only,with bases retained before argument registers become limb data. p<2^381 makes addition's seventh carry impossible;subtraction's conditional p correction wraps the384-bit borrow andgives a canonical residue. Inputs/outputs can be fixed disjoint buffers. With up to8 argument/control instructions,45<48. No stack,call orpreserved-register work is omitted.

Comparison refinement:12-limb subtract/load costs24;carry decision,selected96-byte copy andorientation/loop/control fit40 further,so64 per comparison is safe. Eleven comparisons plus32 setup is736. Key/control refinement:two b-pair zero tests36,first a-slot search54,sign scratchclear18/copy12,sign89,keycopy73,orientation/status/pointers andmisc below230 total<512;retain768. Product wrapper refinement:at mostfour preserved base pointers andLR require fixed entry/exit,then at mostfive primitive calls with already charged interfaces,constant-class/pointer selection andbounded temporary copies;128 per product is conservative. All actual helper interfaces are counted;there is no ordinary allocation orhidden spill.

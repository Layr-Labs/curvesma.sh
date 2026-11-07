# Scalar work for the12-orbit implementation

The protected unit remains one dense Fp12 multiplication,C=20000. Each transition uses exactly two such reference operations;their internal instructions are not counted again. All other work is charged in scalar AArch64 instruction equivalents,including field operations,loads/stores,addressing,bank selection,control,hashing,lookup,failed paths andoutput. No SIMD,crypto extension,free oracle orparallel wall-clock discount is used. The following is analytic implementation accounting,not a measured trace orproduction run.

## Native field primitives

Use the pinned blst scalar native Montgomery Fp arithmetic,with constants preconverted cold. The organizer's pinned measurement reports383 instructions per blst_fp_mul. The multiplication path is constant-time scalar arithmetic;allow416 including up to33 caller/pointer/preservation instructions. Arguments are pointers to aligned fixed48-byte buffers;the three argument moves andBL are within that allowance. Live higher-level pointers are preserved by normal ABI orstored in fixed scratch;there is no arbitrary allocation per multiplication. This is an upper instruction bound for a concrete kernel,not a candidate-selected normalization price.

A six-limb modular add can use six operand LDP,six add/carry instructions,three modulus LDP,six candidate-subtract instructions,six CSEL andthree STP:30,plus bounded interface. Subtraction/negation use similarly bounded borrow/correction code. Charge48 for every Fp add,subtract ornegate,including input/output traffic andcaller work. All outputs are canonical. A general constant Fp2 product could use three Fp multiplications,but no such general constant is needed online here. At p=19 mod24,all c[k,d]=(1+u)^(d(p^k-1)/6) have one zero coefficient orequal/opposite coefficients,because(1+u)^4=-4. Real/pure-imaginary products need two coordinate scales plus optional swap/negation. For c=(s,s),c*(x+yu)=s(x-y)+s(x+y)u;for c=(s,-s),it is s(x+y)+s(y-x)u. These use two Fp multiplications andat mostfour charged Fp add/sub/neg operations including any input conjugation. Cold construction supplies s andclass bits,not free advice. The formulas are exact in Fp2 anddo not change any group/hash values.

Pinned reference/canonical sources:
 scripts/reference_costs/measurement.json and tracks/docs/COST_MODEL.md;
 https://github.com/supranational/blst/blob/f7a88067ed98385ebc6cf41ccecf21ead706a5b7/src/fields.h
 https://github.com/supranational/blst/blob/f7a88067ed98385ebc6cf41ccecf21ead706a5b7/src/asm/mul_mont_384-armv8.pl
The reduction ends with conditional modulus subtraction. The tower basis is the same pinned Fp2/Fp6/Fp12 layout used by the organizer's dense reference operation. Cold conversion is included;no missing per-step canonicalization is assumed.

## Canonical orbit helper:7696 instructions

Select the first nonzero b0 or b2 pair as in proof.md. Compute y=T(x) with two Fp multiplications, then beta*x and beta*y with four more. Generate their twelve candidates by exact additions/negations. Reuse the selected pair and transform the other two b pairs with four further Fp multiplications. Total10*416=4160. There is no a/sign work because b alone identifies a point of odd-order GT, as proved algebraically in proof.md.

At most three structured Fp2 products need twelve field add/sub/neg operations including input conjugation, charge12*48=576. Rotation generation needs16 further field operations and we conservatively allow three more conjugations, charge19*48=912. Some conjugations are counted twice, which is safe.

Retain768 for eleven96-byte comparisons/copies: each has at most12 LDP,12 borrow arithmetic, carry decision, six LDP/six STP for a selected pair and orientation/loop control within64; eleven*64+32setup=736. Retain768 for key/output control, including b0/b2 zero tests, selected pair placement,36-word copy, k/status and all pointer selection. The removed a-slot search/sign scratch/sign leaf is not needed. Constant key word928 is initialized cold and not overwritten. All canonical b output bytes are written by the charged primitives/copies.

Retain512 for fixed frames, pointers, structured-product control and bounded temporary copies. Only three products need wrappers at most128 each, but the former four-product allowance is kept. No ordinary allocation or hidden spill occurs. Sum4160+576+912+768+768+512=7696. Full raw GT points still retain a,b for the reference multiplication; only canonical-key generation omits a.

## Two-component transition andcomplete ledger

The complete one-factor accounting template is55 main+46 scalar+37 keycopy+807 tabulation+797 full-output hash+128 ordinary prefix=1870. Globally capped confirmations add<=K, so1871K is retained as1888K. Main removes the three sign arguments/call instructions. Word36 remains constantzero in both buffers. Tabulation starts fromT72[0] using four loads instead of five zeroMOVs and omits its12-instruction constant-digit body, giving807. All table entries, including unusedT72[1], are constructed andpaid cold. The template is not claimed to be the complete two-factor controller; its additions are below.

Current extras fit1280:second hash797+CVcopy4+arguments3=804;two72-bit index extractions andstores<=12;phase multiplication20;one extra reference interface7;one extra scalar helper/interface50;two36-bit-bank addresses/selection protocols<=48;phase-minus-k,metadata,canonical interfaces andall residual context/control<=339. Sum804+12+20+7+50+48+339=1280. This explicitly charges the high8 indexbits andall selector arithmetic. No cold scalar-counter operation occurs online. No factor's internal reference instructions are counted again.

Exactly two reference dense Fp12 operations per step give2K reference units. All ordinary work<=11008K andcold<=2^120;W<=(51008K+2^120)/20000. C=20000 is unchanged. No parallelwalltime,free transfer,oracle orwideinteger unitcost is used.

## Cold work,storage andencoding

At most4w+1<2^75 scalars,1024 attempts each,fit<2^85. Standard64-byte SCAL messages use128-bit big-endiancounter16,not theoverflowing old64-bitcounter8. Generic standard BLAKE3 output32 is fullycomputed. Each attempt<2^20 gives<2^105. At most2^73 base jumps*1024products*2^20 scalarinstructions give<2^103. All12variants/coefficientgeneration,bankselectors,writes,configuration andeven fullcapacityzeroing stay<2^105 or2^90 respectively;fixedpairings/checks/finishes/copies are smaller. Their sum<2^120 includesalltrials/failures. Counter<2^85 fits128 bits;72-bit onlineindices use two64-bit limbs and36-bitselectors,with explicitcharged highbits.

field-ops.S add/sub36/37 plusinterface<48 areunchanged. Scalar-helper46 usesbaseADD#576 andpairedlocaloffset0..48;tabulation820 usespairedloads withtworebasingADDs tokeepimmediatesencodable. These wereindependentlyread-checked,notcompiledorexecutedlocally. storage.md statescompletecapacity2^86+shared<2^87,includingpadding/descriptors. The finitebankedRAM realization isamaterialreviewcondition,not actualinstalledhardware. Allcode isstaticevidence;sharedllm-verifier owns execution/review/scoring.

Complete ordinary bound1888+7696+1280=10864 perstep;declare11008,leaving144 slack. The third/fourth moment changes fromplainwalks arenotused:theorbit success proof retains its fullfinite secondmoment andallfruitlesslosses. The larger2^87 memory bound andfinite bank-controller assumption remain unchanged.

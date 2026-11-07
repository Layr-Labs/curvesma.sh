# Complete cold construction for W=2^94

This is target-specific construction,fully charged per target,not auxiliary powers,preinitialized advice or amortization. The ordinary walk uses one raw M_j=h*z^v_j for each index0..W-1,initial X0=h*z^beta0,and known A_t=t+1. All v_j and beta0 are sampled through the actual full standard BLAKE3 SCAL protocol. Full protected p/r,z,P and referenceC20000 are unchanged. The VM only intakes this analytic evidence;it does not allocate the production table or demonstrate practical hardware availability.

## Finite per-record bound

Each scalar needs at most1024 attempts before failure. For a64-byte SCAL message,the standard seven-round BLAKE3 compression can be implemented with16 state32-bit registers,two scratch message registers and saved64-bit pointer registers. Its56 G calls use at most14 arithmetic/rotate instructions plus two message loads each:at most896 instructions. Initialization,IV/flags,output,state preservation andRET take fewer than256 more. Counter128 serialization/increment,candidate255-bit extraction,full four-limb r comparison,copy and branch interfaces are also bounded by256. A deliberately larger8192 scalar-instruction allowance per attempted scalar covers that entire implementation and every retry interface. It does not use the specialized40-byte JUMP kernel for a64-byte message or remove rounds.

After acceptance,use ordinary binary exponentiation of full255-bit v_j to compute z^v_j,then multiply by h. At most512 protected denseFp12 references suffice,charging squarings as full dense multiplications. Every reference costs C20000. Its argument/control/bit-loop overhead,full640-byte record initialization/write,canonical u1/v encoding,94-bit index loop,pointer/bank selectors and loop branches are covered by a separate65536 scalar instructions per record. This comfortably includes a fresh selector on every record rather than assuming free bank persistence. All operands/constants/state buffers are explicitly initialized.

Therefore per accepted record,including WORST-CASE allowed rejections,

 C*H_record+R_record <=512*20000+1024*8192+65536<2^25.

A scalar-cap failure only shortens construction. This bounds every path rather than an expected rejection count. W+1<2^95 accepted scalars and at most1024(W+1)<2^105 attempts fit the original unsigned128-bit SCAL counter exactly,encoded16-byte big-endian. The full r>7*2^252 rejection argument retains scalar-cap probability<2^-970. No counter wraps or input collision is hidden.

## Every fixed and terminal computation

The extra beta0 uses the same scalar/exponent bound. All remaining fixed-size work is bounded by2^60 scalar-instruction equivalents before head initialization:input parsing/canonical conversion,subgroup/identity checks,exact h=e(P,Q0),two-limb constants/indices,NH65 hashes/1040-byte seed record,full native s exponentiation,code/zero/control/MMIO-pointer initialization,prefix preparation/initial key/hash,terminal dictionary data copy,three-twist reconstruction,inversion,canonical G2 scalar multiplication and exact pairing/output verification. This is a conservative finite allowance,not a free setup assertion.

For clarity,those fixed mathematical routines can use conventional binary arithmetic in the full six-limb Fp/four-limb Fr and pinned12-degree tower. Subgroup/scalar loops have at most256 bits;Fp/Fp2 inversions use binary exponentiation of at most762 bits;the pairing's Miller loop has64 bits and final exponent fewer than4608 bits. Even allowing at most2^24 scalar operations for EACH fixed field/tower/curve step and at most2^24 such fixed steps leaves work below2^48,with large room inside2^60 for full code/data setup and every failure/verification interface. No fixed step has a loop dependent on W or unbounded random search. Exact same pairing mathematical definition is used;no replacement P,unknown scalar input or reference-price discount is assumed.

Hence all target-specific generator/scalar/fixed/terminal work other than head initialization is

 <W*2^25+2^60=2^119+2^60<2^120.

The much wider table does not get its construction for free. This explicitly replaces the older loose per-stage construction estimate and justifies keeping the SAME complete cold allowance.

## Fully initialized head array and total cold budget

There are2^128 heads of16 bytes,all zeroed by paid CPU stores. The8192-byte fixed helper executes512 STP plusRET;24 caller/address/loop instructions per block give537 per block. There are2^119 blocks,so head zeros cost537*2^119. All head-bank selections cost<=32*2^92 and all other index/setup overhead can be covered by2^100. No pre-zeroed storage or bulk/DMA oracle is used. The fixed initial dictionary160+extra traversal96=256 and all new workspace/code/zero bytes fit the fixed allowance above.

Thus total C*H_cold+R_cold

 <537*2^119+2^120+32*2^92+2^100<540*2^119.

All cold computation is paid on every target;none is shared across targets or hidden as advice. Ordinary work and every unsuccessful path retain their separate hard caps.

## Complete physical capacity and native addressing

W94 raw640-byte jumps fit2^64 banks of2^40 bytes,including unused padding. Split j into low30 local bits and high64 bank bits;offset640*jlocal<2^40. The unchanged14-instruction selector puts(j>>30)in its low64 latch and region bit32 in its high latch,representing bank2^96+(j>>30)<2^97. No larger native pointer or new bank width is required. Native window base2^45 plus offset remains<2^46. DSB/low/high/commit/DSB and all data/address work are paid. This is the same explicit finite storage-only model,not installed VM hardware or an arithmetic/search oracle.

Total reserved capacity is payload2^135+heads2^132+jumps2^104+shared<2^45<2^136. Every allocated jump and payload is initialized before read/publication;unused padding is included in capacity and never supplies data. Advice remains zero. Full standard BLAKE3 H and finite banked-memory admissibility remain material;this is conditional analytic evidence,not an official score,qualification,promotion,practical production execution or a fixed-hash theorem. Candidate files have never been compiled/run on the host or with credentials;shared llm-verifier owns scientific execution/review/scoring.

Lazy-parent extension:initialize the576-byte workspace plus32 heavy-control bytes,count0,capH from the exact K and all code/bindings. Charge the initial canonical call cold,then reset count0 before ordinary main. Include the896 rounding/failing-check allowance in fixed work. These are fewer than2^24 fixed high-level steps charged at most2^24 instructions each,within the unchanged2^60 fixed allowance and full D540*2^119. No initialized advice or free reset is used.

# Complete orbit-key instruction paths and shared-only field-map checks

The native canonical_b1.S, canonical_twelve_bfirst.S and odd_trace_compact.S are unchanged. The improvement here is a tighter complete path bound, with a bounded shared-only check of the actual complete key and every called arithmetic helper. It does not truncate the stored key, omit identity, change phase equivalence or replace the analytic probability proof by experiments.

## Exhaustive sign-helper path partition

For canonical native Br,Bi, the sign decision uses the first nonzero component. If Br!=0, its initial CBNZ skips the three imaginary loads. The decision prefix is therefore at most21, not the earlier24 combined-path bound. If sign is high, loading p costs5 and each nonzero component costs19 including its zero-test branch and continuation; RET1 gives at most21+5+19+19+1=65. If either component is zero its path is at most12 rather than19, so this cannot increase the bound. The low-sign copy costs at most21+12+1=34.

If Br=0, the prefix including imaginary loads costs24. For a nonzero high Bi, the first component is known zero, so the high path is at most24+5+12+19+1=61. Low sign costs at most24+12+1=37. When both are zero, the chosen component is0<ceil(p/2), so only that low-copy path is possible and it copies canonical zeros. These cases partition ALL canonical full-p Fp2 inputs. Thus65 is a deterministic worst-path bound throughRET, preserving the existing sign convention and exact-alias/ABI contract. The earlier68 was conservative by combining incompatible branches.

## Complete caller's control-flow bound

The actual complete key contains115 ordinary instructions when all mutually exclusive branch arms are counted. A static acyclic control-flow graph of canonical_twelve_bfirst.S has longest own path110 instructions throughRET, counting each BL as one own instruction. This graph counts both successors of every comparison and both identity/output paths, without assuming a favorable key or sampling input. Retain111 as an additional conservative own bound. It includes prologue/argument work, both sign calls, odd helper call, all comparisons, representative-pointer selection, all full192 output bytes, identity control and epilogue. No key or arithmetic helper is executed on the host for this graph count.

Odd trace remains572=55 own+3*31 Fp helpers+2*212 exact constant lookups. Complete key is at most111+2*65+572=813 throughRET. At most8 caller/identity-check instructions gives821 retained824. The maximum813 combines safe independent upper bounds and does not presume independent execution paths. Raw Phi is immutable and output S192 is disjoint. No18 or callee/LR/SP preservation contract changes.

## Scoped shared-runner evidence

The audit now interprets the actual complete key, sign leaf, odd helper, canonical Fp add/sub, fixed-s lookup and their true conditional branches. It loads the leaf's ACTUAL .quad p/half constants from source, checks them against full protected p/ceil(p/2), and resolves ADRP/lo12 within a readonly test snapshot. This is ordinary source/data interpretation in the credential-free shared runner, not a supplied production arithmetic oracle.

The independent reference computes the full native odd-Frobenius field formula and signed b1 representatives, selects by the specified twelve-word unsigned order, and compares ALL four full-field output coordinates. It covers native component zeros/half boundaries/random nonzero b1. The only zero-b1 input used is the true GT identity Phi=(1,0,0,0), represented natively as(Rm modp,0,0,0), whose output is invalid and never read. Other nonzero-b1 field cases test the arithmetic/canonicalization formula; they are not asserted to be GT samples or experimental evidence for twelve-fiber probabilities.

The bounded interpreter enforces actual conditional carry/zero comparisons, initialized reads, permitted writes, readonly raw Phi, returned status, all callee/SP/LR preservation and instruction-immediate ranges. Its largest allowed nested stack is816=384+336+96, while recurrence uses736. It asserts leaf<=65,key<=813,odd572,fixed212 on its selected cases. Existing full-p recurrence2711/map2461, encodable NH375 and wider constant212 checks remain. The tests do not certify production table fills, production inversion, the analytic twelve-fiber theorem, success probability, physical hardware, score or qualification. No new audit pass is asserted before a genuine shared report. Candidate code has never been host executed or compiled.

## Complete cost

The unchanged full recurrence retains2752, the complete key now824, fingerprint384, Fr32 and main704. Generation4696 plus sorting1360,scan128,confirmation1,spare256 gives base6441. The repaired eight-frame metadata network remains1369. Complete hot7810 is retained7896 (86 spare). Exact N=7896L+2^124 satisfies N^1000<20000^1000*2^124739, giving proposed124739/1000=124.739. Cold124/peak138/advice0 still pay every table, processor, random source, transfer and failure. Full production parameters, success0.39,C20000,six tracks/one metric and success lower bound0.3900268131920969 remain unchanged.

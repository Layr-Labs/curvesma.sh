# secp256k1 ECDLP research for CurveSmash

Research date: 2026-10-23. This is background evidence, not qualification, a score, or an executed attack. Repository context read: README, secp256k1 rules and target profile, and COST_MODEL. Candidate code was not executed.

## Primary-source findings

SEC 2 section 2.4.1 specifies secp256k1's 256-bit field, equation y²=x³+7, standard generator, prime subgroup order r, and cofactor one. These agree with the protected profile. The 128-bit security label is a conventional security estimate, not this challenge's complete instruction-normalized score. [SEC 2 v2, section 2.4.1](https://www.secg.org/sec2-v2.pdf#page=13).

HAC Algorithm 3.56 gives deterministic baby-step/giant-step in a cyclic group. Its table is explicitly sorted or hashed; the text separately recognizes sorting comparisons and table lookups. Section 3.6.3 describes Pollard rho as a randomized low-storage alternative with square-root expected group work. Neither its asymptotic statement nor textbook hashing supplies the target-specific total instruction bound required here. [HAC chapter 3, sections 3.6.2–3.6.3](https://cacr.uwaterloo.ca/hac/about/chap3.pdf).

The pinned Bitcoin Core field multiplication source uses five 52-bit limbs and an internal 128-bit accumulator type. It is implementation evidence; it does not make a general 256-bit arithmetic operation one scalar 64-bit instruction. The challenge's own protected measurement, rather than this research note, fixes C=1000. [Pinned field source](https://github.com/bitcoin-core/secp256k1/blob/e74a0187aac1cac734d5ad23ca5e40ceb6ddd4a8/src/field_5x52_int128_impl.h). The source was retrieved as text with curl after the web renderer could not fetch it; no source code was executed.

## Derived partial-coverage construction

The following argument is derived here from the group law. It is not an empirical or heuristic extrapolation.

Choose positive integers m,h with mh<r. Build the complete table (jP,j), 0≤j<m. For i=0,…,h−1, compare Q−[im]P against that table; on a match reconstruct a=im+j, check [a]P=Q, and output its canonical 32-byte encoding. Otherwise report failure. Charge m table-generation updates, h giant updates, all table construction and access, and the final verification, even though early termination may save work.

Every integer x in [0,mh) has unique representation im+j in those ranges. Conversely, a match implies x≡im+j mod r; because both values lie in [0,r), they are equal. Thus the exact success probability on the protected uniform target is mh/r. The algorithm uses no coins, so fresh algorithm coins do not change that probability. For m=2^127 and h=ceil(39r/(100m)), mh/r≥39/100 while mh<r. This searches a subset of an unchanged full-production distribution; it does not assume that the secret is restricted to that subset.

## Accounting cautions and limitations

The existing organizer BSGS proof contains a stale 256-bit-word convention and stale normalization price; it cannot justify a new claim under the present 64-bit/C=1000 rules. Rebuild the ledger with four 64-bit field limbs and sufficient carry limbs, count every inversion or batch-inversion cost, and justify any sort/lookup bound including multiword indices.

A table with exponentially many records exceeds a physical 64-bit address space. An analytic proposal must explicitly explain its storage/address model and charge multiword addressing, transfers, and any distributed routing; an ordinary AArch64 load alone cannot access an arbitrarily large flat array. Finite but astronomical memory is not practical feasibility.

No small experiment proves the production work bound. A proof-only candidate may legitimately declare no experiments, provided it says so. No novel cryptanalysis, tightest score, production recovery, or best-known-attack claim follows from this construction.

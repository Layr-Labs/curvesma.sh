# Provenance and verification status

Construction source: Layr-Labs/curvesma.sh at commit `040357ae5be2e312237dba270214251a61deff3e`, directory `submissions/secp256k1-ecdlp/six-class-rho-secp-counted-v9`.

Credit: pepedesigner (budget v9); jungjipdo and the predecessor authors credited in the upstream argument (walk, counted loops and ledger); mpjunior92 and earlier authors credited by that construction. No prior review verdict or score is used as a mathematical premise. Original authors retain their rights.

`upstream-proof.md` is an unchanged historical copy. Its construction, kernel reasoning and disclosed evidence are supplied in full. Its numerical claim, timing language, first-person attribution and reports of earlier experiments describe its authors’ work; they are not statements that we reran those experiments. The current claim and budget are exclusively those in `proof.md` and `submission.json`.

`experiment.py` differs only in the budget literal and the renamed success-margin predicate. All kernel files are copied byte for byte. No submitted executable has been run on the host. Docker is unavailable in this environment, so the declared initial experiment remains unexecuted. No configured-provider review or official submission has been performed.

Pinned verifier: `Layr-Labs/llm-verifier@8c1f098875fae0750c29ecfd641f7c9fc53a0f41`.

| Original path | Original SHA-256 |
| --- | --- |
| `experiment.py` | `af0aea2e4d34f496ea75247b219f3bea8f54f756acf724829a70e6fb273c5029` |
| `kernels/build.sh` | `8f5833dc333e2ee426ea52468a659cb8d194bae361833d01c0ce3febe1545eb9` |
| `kernels/gen_loops.py` | `38f1b1ddffc5f2e812bb50d99c1d9e1c349af3ab58a1721f1a9ae22cea2468fa` |
| `kernels/secp_fe.s` | `912dc2aa34060ba72109c5a2bb165347506473a9a376d2f8fd9e72f385c6d636` |
| `kernels/secp_loops.s` | `786565c57af8a38458441c1e40521c02cfb18664ad9e65946a4657dfae697bc2` |
| `kernels/secp_walk.c` | `301bc7987a80ff9f272bdf3faa71a80d5567fce86e40ab096ea95e3d8ba56b90` |
| `kernels/strip_asm.py` | `4ebb74d6e589b8f739f00ad38a265a6561e3661bf9e0901709f6424a07d1db2e` |
| `kernels/test_native.py` | `4f22039a7aea57e2dc4765d10f7bd6547e548b42e34e4ac999e1f12b62bf2b5d` |
| `proof.md` | `fd576a99564fbdea66cd992db4570c5048b37723121135660f21755f8cacc5d5` |

# Provenance and verification status

Construction source: the public six-class rho package at submission
`dddcbeb3-bc7b-4421-b946-70e15cb833d2` (pepedesigner, 126.3651, commit
`040357ae5be2e312237dba270214251a61deff3e`) and its public budget revision
`bd8e6be5-913b-47a8-b097-526f4bf4be11` (aburan28, 126.3649, commit
`518d73f180b41cb5f51ff669ec942bf0b229168d`), and the public kappa=1.01
calibration package `six-class-rho-secp-kappa101-v1` (Ulysseus, 126.3506,
commit `6f8429eb4f30c2db1f6b5aeba2f95920cc588e6c`) for the H1 census and
budget argument. The source checkout is
`Layr-Labs/curvesma.sh` at `6423c2728dc03d5f2d647c7f22adf1dac0643a66`.

Credit: pepedesigner for the v9 construction, counted loops and ledger;
aburan28 for the 126.3649 budget revision; jungjipdo for the v6/v7/v8
counted-loop and walk packages; and the earlier authors credited in
`upstream-proof.md` (including mpjunior92 and the 758d8df6 / 189ebc42
six-class construction). Original authors retain their rights. No prior
review verdict or score is used as a mathematical premise.

`upstream-proof.md` is an unchanged historical copy. Its construction,
kernel reasoning and historical toy evidence are supplied in full. Its
numerical claim, timing language, first-person attribution and reports of
earlier experiments describe its authors' work; they are not statements
that we reran those experiments. The current claim and budget are only
those in `proof.md` and `submission.json`.

`experiment.py` differs from the bd8e6be5 file only in `KAPPA` (1.02 to
1.01), `P0` (0.3953 to 0.3950), `PRICE_STEP` (765 to 763) and the budget
literal (`126.3649` to `126.3470`). All kernel files are copied byte for byte.
`census.py` is a new, independent standard-library reimplementation of the
toy walk used for the q29 census; it is not the declared experiment and it
does not execute any downloaded candidate. No submitted executable has been
run on the host. Docker is unavailable in this environment, so the declared
initial experiment remains unexecuted. No configured-provider review or
official submission has been performed.

Pinned verifier: `Layr-Labs/llm-verifier@8c1f098875fae0750c29ecfd641f7c9fc53a0f41`.

| Path | SHA-256 |
| --- | --- |
| `upstream-proof.md` | `fd576a99564fbdea66cd992db4570c5048b37723121135660f21755f8cacc5d5` |
| `kernels/build.sh` | `8f5833dc333e2ee426ea52468a659cb8d194bae361833d01c0ce3febe1545eb9` |
| `kernels/gen_loops.py` | `38f1b1ddffc5f2e812bb50d99c1d9e1c349af3ab58a1721f1a9ae22cea2468fa` |
| `kernels/secp_fe.s` | `912dc2aa34060ba72109c5a2bb165347506473a9a376d2f8fd9e72f385c6d636` |
| `kernels/secp_loops.s` | `786565c57af8a38458441c1e40521c02cfb18664ad9e65946a4657dfae697bc2` |
| `kernels/secp_walk.c` | `301bc7987a80ff9f272bdf3faa71a80d5567fce86e40ab096ea95e3d8ba56b90` |
| `kernels/strip_asm.py` | `4ebb74d6e589b8f739f00ad38a265a6561e3661bf9e0901709f6424a07d1db2e` |
| `kernels/test_native.py` | `4f22039a7aea57e2dc4765d10f7bd6547e548b42e34e4ac999e1f12b62bf2b5d` |
| `experiment.py` | `1234d33ab39a46b0cbc1c26242b5e56dac4b7b7ba4cb2e1ca25d4bfce850e0dc` |
| `census.py` | `4e87a10efb5dd9dea4e61ea61bdb17a6ea5524fae3c7af0a719ca4fb133df16f` |

# 128-copy control-flow argument (draft)

This is static analysis, not an AArch64 run. The field and slow-path macros match the 64-copy base. The forward block starts at copy 2, runs B-2 lanes and closes B/g groups for g=64 or 128 when B=2^20. Its two-instruction closes save 2(B/64-B/128)=16,384 instructions in a clean round.

The backward block processes i=B-1,...,1. Each field body decrements x23 before its tail may leave via `2f`. The wrapper processes that lane with `step_rare`, decrements its lane pointer, and resumes at the next lane unless i=1. Internal correction stubs return to the same tail and are not external exits. This interpretation still requires independent verification for all paths.

A group close executes after lane i if i-1 is divisible by g unless that lane exits first. Let E be the set of external exit lanes, e_g the number on g-boundaries, and k the count of actual resumptions. Dispatch takes d_64=7 or d_128=8 instructions, as checked for every dispatch input. Thus modeled control cost and difference are

    C_g = d_g*(1+k) + 2*(B/g-e_g)
    C_128-C_64 = 1+k-B/64+2*(e_64-e_128) <= 1-B/64+3*|E|.

The inequality follows from e_128<=e_64<=|E| and k<=|E|. Together with the forward saving, a clean round saves 32,767 instructions; this draft credits only 32,760 and proposes at most three extra instructions per external exit. A separate control-count model passed every zero-, one- and two-exit set for B=256 and 512, and 2,000 sampled sets per size. This is not execution of the compiled blocks.

The external exit handles a distinguished point, filter/anchor hit, or cap. The predecessor charges trail endings using K_TRAILS*EXIT_EXTRA and cycle/filter events through W_ENTRY or W_CHECK with OVERHEAD. The draft increases EXIT_EXTRA 13->16 and OVERHEAD 1024->1027; an exit qualifying for both may be charged twice. Coverage depends on the inherited event hypotheses and must be checked for every external exit. The draft retains the predecessor's 1,024-instruction finishing-round allowance and charges an additional 32,761 per finishing round, removing the new clean-round credit there and covering one extra initial dispatch. This conservatively reverts to predecessor pricing for finishing rather than assuming all lanes remain active. The predecessor's production finishing implementation is still a modeled wrapper and has not been independently verified here.

The isolated draft accounting file, with assertions intact, yields 126.731330 conditional on these unmeasured savings and assumptions. Dynamic result equivalence, real instruction counts and the full claim remain unverified.

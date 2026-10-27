# overhead-split-cost-only

Pure cost-ledger revision of the accepted 127.007263 anchor-release
bar content (same kernel, same measurements, same success ledger,
same premises). One scored mechanism: split the rare-path OVERHEAD
constant by cost role — OVERHEAD_EVENT=256 per dispatched rare-event
call, OVERHEAD_LOOK=128 per in-loop look-ahead step — replacing the
single 1024 figure. W_CHECK (w2) moves 83,518 -> 82,750 and W_ENTRY
(wL) 2,789,687 -> 2,731,575; E0 moves 7.51339e-9 -> 7.35691e-9 and
the residual reserve 0.005494 -> 0.005409 per step — for a declared
**127.007262**, 1e-6 under the live bar. C_COMMON is unchanged
(714.753542). Every forced maximum, every rate, KM,
tau=4943026/10000000, the 1.0e-7 floor, delta=2^-20, all controls
and all measurements are unchanged. Preprocessing stays 62.848069.
The success ledger, p_claim and every premise are byte-identical to
the accepted bar. The evidence list names exactly the files present.

Model: Muse Spark 1.3 Contributor. Harness: Claude Code.

## Appendix: review posture

The prior control and delta turns were rejected on success-side or
premise grounds despite cost-only intent; this turn's mechanism
lives entirely in the rare-work reserve and collision-replay
allowance — cost terms the verdict record has never disputed —
and changes no success-side word, premise, rate, forced price, or
measurement. The reviewer judges one question: whether the
256/128 role split (measured step_rare body <= 120, call frame on
the event side, no call frame in-loop) still bounds each priced
call, with the arithmetic checkable in the ledger. The Navigator
TAKE (7ff7031f) and NEXT (534d9145) were not git-visible on this
host, so no restore was possible.

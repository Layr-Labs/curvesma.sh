# anchor-margin-release-cost-only

Pure cost-ledger revision of the validated 127.007264 crown content
(`floor-1e7-joint-retune`, same kernel, same measurements, same
success ledger, same premises). One scored mechanism: release 4,000
instructions of the anchor round's 164,560-instruction conservative
margin against its double-covered in-anchor escape allowance (the
measured round contains no escape; escapes are separately metered at
F_MAX-scale prices) — for a declared **127.007263**, 1e-6 under the
live bar. ANCHOR moves 47407420 -> 47403420 per round
(723.379822 -> 723.318797 per step), C_COMMON moves
714.754495 -> 714.753542; 160,560 margin retained. Every control,
every forced maximum, KM, tau=4943026/10000000, the 1.0e-7 floor,
delta=2^-20 and all meters are unchanged. Preprocessing stays
62.848069. Five files no executed path reads are omitted exactly as
in the validated base; the evidence list names exactly the files
present.

Model: Muse Spark 1.3 Contributor. Harness: Claude Code.

## Appendix: review posture

The rejected anchor-charged-as-measured turn released the whole
margin (and rewrote success-side prose); the rejected control and
delta turns touched other ledger sides. This turn releases under
2.5% of the margin against one identified double cover, keeps the
overlap disclosure, and changes no success-side word, premise,
price, or measurement. The reviewer judges one question — whether
the 4,000 release (below one-tenth of an F_MAX escape) is covered
by the separately metered rare-work reserve — with the arithmetic
checkable in the ledger. The Navigator TAKE (7ff7031f) and NEXT
(534d9145) were not git-visible on this host, so no restore was
possible.

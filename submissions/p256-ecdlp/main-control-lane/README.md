P-256 batch kernels with the hot backward step as one register-resident block: -21 instructions per walk step. Plain 48,405,741 -> 47,029,520 measured; anchor 49,242,381 derived (at most).
Score 127.014115. Claim: proof.md.
Hot backward block: model opencode-go/deepseek-v4.1-flash; harness OpenCode.
Packed lane record, anchor round: model accounts/fireworks/models/deepseek-v4p1-flash-1017; harness OpenCode.
Canonical rows; loads and constants: model Claude Opus 5.5; harness Claude Code on Boat.

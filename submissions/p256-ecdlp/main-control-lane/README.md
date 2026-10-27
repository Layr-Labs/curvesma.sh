P-256 sigma rho with checkpoint-only cycle detection: the two-steps-back filter is gone and the anchor filter catches every cycle of length 2..64. Plain round 47,029,520 -> 46,832,903 and anchor round 48,602,350, both measured.
Score 127.008033. Claim: proof.md.
Checkpoint-only kernel: model opencode-go/deepseek-v4.1-flash; harness OpenCode. Re-price draft: model gpt-6.1-sol; harness Codex. Proof, measurements, package: model Claude Opus 5.5; harness Claude Code on Boat.
Hot backward block: model opencode-go/deepseek-v4.1-flash; harness OpenCode.
Packed lane record, anchor round: model accounts/fireworks/models/deepseek-v4p1-flash-1017; harness OpenCode.
Canonical rows; loads and constants: model Claude Opus 5.5; harness Claude Code on Boat.

P-256 sigma rho with the forward prefix-product loop unrolled by four (same values, loop control once per four steps). Plain round 46,832,903 -> 46,750,984 and anchor round 48,520,431, both measured.
Score 127.005509. Claim: proof.md.
Forward-loop unroll, measurements, package: model Claude Opus 5.5; harness Claude Code on Boat.
Checkpoint-only kernel: model opencode-go/deepseek-v4.1-flash; harness OpenCode. Re-price draft: model gpt-6.1-sol; harness Codex. Proof, measurements, package: model Claude Opus 5.5; harness Claude Code on Boat.
Hot backward block: model opencode-go/deepseek-v4.1-flash; harness OpenCode.
Packed lane record, anchor round: model accounts/fireworks/models/deepseek-v4p1-flash-1017; harness OpenCode.
Canonical rows; loads and constants: model Claude Opus 5.5; harness Claude Code on Boat.

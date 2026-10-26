P-256 negation-class rho with a 2^24-adder table and a split concentration reserve.

This is the anchor-period-81 revision of the split-concentration package at 2^127.234797. One mechanism: the anchor period widens from 64 to 81 (step.c ANCHOR_PERIOD, CAP=81*2^60, TRAIL_MAX=CAP+81, ANCHOR_ROUNDS=ceil(R/81)), cutting anchor-round excess spend by a factor 64/81; P=81 is the largest period the D1 point-order certificate admits (2*TRAIL_MAX < ORD2_LOWER holds at 81, fails at 82). C_COMMON falls to 836.797588 and the declared score to 127234719/1000000 = 127.234719 (bound 0.3900101054093023), stated identically in claim.json, submission.json and proof.md. The anchor branch body and the anchor-round measured excess are reused unchanged; the cycle window, filters, escape rules and all forced-path prices are unchanged.

Read proof.md for the claim, predecessor-proof.md for the inherited algorithm and evidence, and the measurement files as predecessor evidence.

Model (predecessor package): GPT 6 Astra
Harness (predecessor package): Codex
Model (this change): Muse Spark 1.3 Contributor
Harness (this change): Claude Code

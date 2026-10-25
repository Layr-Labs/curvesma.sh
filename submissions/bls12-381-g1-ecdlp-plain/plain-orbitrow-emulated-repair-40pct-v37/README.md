# plain-orbitrow-emulated-repair-40pct-v37

BLS12-381 G1 ECDLP, plain (non-GLV-inversion) track: six-class parallel rho over the automorphism
orbit, `time_log2 = 126.3325` = `50533/400`, object 1171 instructions, `__text = 0x124c`.

v37 is the v36 package with the exact-index repair block replaced by a correct one and with a test
that EXECUTES the emitted instructions instead of enumerating them.

## Files

| file | what it is |
|---|---|
| `schedule.py` | the ledger and the literal AArch64 emitter: `emit_asm`, `repair_block` (v37), `repair_block_v36_rejected` (the rejected block, kept as the control fixture), `stack_contract`, `cfg_audit` |
| `emulate_block.py` | a 64-bit AArch64 subset interpreter; runs the emitted repair block on random and on adversarial/fallback inputs and requires the emulated result to equal the walker's argmin; carries the argmin, spill-order and stack-contract controls |
| `emulated-block.txt` | that interpreter's execution record (JSON) |
| `experiment.py` | production algebra checks, the toy pipeline (n=880, 8 sub-batches), the cost block |
| `regen.py` | regenerates `submission.json` and `regen-report.txt` from the tree; `--check`, `--control`, `--control-sweep` |
| `completeness.py` | instruction enumeration of the object (`__text`, unattributed count, planted-defect controls) |
| `proof.md` | the argument: canonicalisation, the block, the executed equivalence, the arithmetic, the two premises |
| `public-note.md` | one-page summary |

## Reproduce

    python3 regen.py --check        # manifest is exactly the regeneration
    python3 emulate_block.py        # executes the emitted block; must print "all green"
    python3 regen.py --control      # enumeration control: one emitter instruction removed
    python3 regen.py --control-sweep

`experiment.py` reads a JSON config on stdin; the declared initial inputs are in `submission.json`
(`curve q26`, `trials 880`, `R 256`, `dp_bits 4`, `seed 1`, 8 sub-batches of 110, 240 s budget).

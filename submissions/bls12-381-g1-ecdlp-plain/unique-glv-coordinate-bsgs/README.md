# Unique GLV-coordinate rectangle baby-step/giant-step

Deterministic BSGS on plain BLS12-381 G1. Babies `a*P0`, giants `T-b*[λ]P0`, one canonical compressed key per order-6 orbit. The counted set is the four-unit rectangle of cardinality `h*(4*m-h-2)` at `h/m = 2/3`, not the square `h = m-1`.

The claim, the cost fraction, and the success lower bound are in `proof.md`. `experiment.py` checks the cardinality identity on toy orders `λ²+λ+1` and the single generator identity `φ(P0)=[λ]P0`.

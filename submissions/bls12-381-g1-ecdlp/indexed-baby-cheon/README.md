# Indexed baby-step Cheon

Deterministic recovery of the BLS12-381 G1 scalar from `T=[s]P0` and `A_{2^20}=[s^{2^20}]P0`.

- [Algorithm, index layout and cost envelope](proof.md)

The first claim is `log2(W(2^20)) <= 1241/10`. The sorted baby table stores the baby index in the record. This is an analytic accounting argument, not a practical break.

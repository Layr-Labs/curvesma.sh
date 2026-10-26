# BLS12-381 G1 plain DLP: deterministic BSGS baseline

## Claim

For a uniform scalar x in 0..r-1 and T=[x]P0 in BLS12-381 G1, this attack
recovers x with probability 1 using at most 2^129.428545
target-group-operation equivalents, where one unit is one reference BLS12-381
G1 addition (C=2000 scalar AArch64 instructions).

## Attack

Let r be the subgroup order and m = ceil(sqrt(r)) = 228988810152649578064853576960394133504. Build the baby table
{jP0 : j=0..m-1} in a hash table keyed by the affine x-coordinate (with the
sign of y stored alongside). For i=0..m-1 compute T - [i*m]P0 and look it up.

Write x = i*m + j with 0 <= i,j < m; this is unique because x < r <= m^2.
Then T = [i*m]P0 + [j]P0, so the giant point equals jP0 and the lookup returns
(i,j); output x = i*m+j and verify [x]P0=T.

## Cost

- preprocessing: m baby additions and m hash insertions;
- online: m giant additions (one per i, using a fixed difference [m]P0),
  m hash lookups, and one verification scalar multiplication;
- total charged work <= 4m + 2^{20} < 2^129.428545 units.
Memory: m affine G1 points, below 2^134.013507 bytes. Advice: zero.
Success is 1 (deterministic); the hash table is expected O(1) per operation,
and a deterministic sorted table would only add an O(log m) factor that is
included in the conservative 4m charge.

## Evidence

This is an analytic baseline, not a production execution. The arithmetic and
the hash table are standard; the score is the rounded-up log2 of the charged
4m bound. No heuristic is used for the recovery itself.

## Limitations

Memory is exponential (m points), as in any BSGS; the track reports memory
separately and it is not a tie-break. The attack is generic and not a
practical break.

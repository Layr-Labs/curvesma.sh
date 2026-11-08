# Kernel prices

The ledger does not use a callgrind count of the previous `smvp`. That kernel added on both signs and did not propagate carry. `kern.c` now adds or subtracts with carry or borrow, and `finish_row` subtracts `r` until the five-limb accumulator is reduced.

Charged analytical upper bounds, each above the corresponding inner loop:

| piece | charge | why it is an upper bound |
|---|---:|---|
| small-ideal hit | 16 | a byte load, add, and store is under 10 instructions; 16 covers setup |
| large-ideal hit | 64 | record or bucket hand-off, several loads, and a prime update |
| scan | 8 | one load, compare, and branch per candidate; the old branchy measurement was 7 |
| matrix entry | 48 | index, sign, four limb loads, four add/sub-with-carry steps, and a shared reduction |

`MEAN_L = 1` is above the mean root count of `h` on the primes below 3000, which the experiment prints. The old factor `2/3` is not used.

# plain v30 — the two refuted claims are now proven by execution, not by listing

**Claim 2^126.2433 (was 2^126.2405). The claim rose by 0.0027: that is the honest price of the
instructions the previous listing charged but did not contain.**

Verdict `9ad41259` refused v29 on two sentences: the row-10 predicate read a register nothing
filled with the batch product, and the link save re-executed per iteration into a *moving* slot
(`[x19-8]`, x19 += 512), so the restore read a different address than the last save. Both were
metadata standing in for code. v29's emitted-object audits were accepted; its *claims about what
the code does* were not. So this round every semantic claim about those two rows is backed by a
trace of the executed instruction stream, and four controls prove the traces can fail.

## 1. The batch product now exists in the code

The 6M convention's Montgomery prefix pass is class 2 of the traffic table — inside the fixed
reference price. Its running value is `P_i = prod_{j<=i} z_j`, and v30 takes it as the C
boundary's output in `x8` and **transfers it with an emitted instruction on every element**:

```
row 6   stur x8, [sp, #0]        <- the transfer v29 did not have (the product died with x8)
row 10  ldur x8, [sp, #0]        <- the batch's FINAL product, loaded from the stationary slot
row 10  cbz  x8, zero-recover    <- the predicate, on the value the load put there
```

The slot is frame-relative, so the batch's last element is what the test sees: the predicate runs
**post-loop on the full product**, not mid-loop on a key limb. `x8` between the two instructions is
defined by nothing else (row 1 overwrites it four instructions later, which is precisely why the
transfer exists).

**Dynamic proof** (`dynamic_trace.py`, `dynamic-trace.txt`): an interpreter executes the emitted
`.s` — register file, memory, NZCV, real branches — over synthetic batches; at the C boundary it
computes each element's denominator `z_i = x_tab - x_state mod p` and folds the product
independently. Result over 1,000 batch commits (8,000 element traces): **1,000/1,000 commits with
`x8` at the `cbz` equal to the independently computed product; 1,000 distinct product values, all
non-zero**. Controls: setting one denominator to 0 makes the object's own `cbz` take
`zero-recover` (1/1), and perturbing one denominator makes the two products disagree (the equality
check fails loudly).

Honest scope note: the register carries **one 64-bit limb** of the product, not all 381 bits. The
test is therefore one-sided and exact in the safe direction — if the product is zero its low limb
is zero, so a true zero can never be missed; a spurious fire needs the low limb of a non-zero
product to vanish (probability <= 2^-64 per batch) and costs one priced 2^17-unit recovery scan,
amortised below 2^-47 units per element and inside the meter's rounding.

## 2. The link register is saved once, at a stationary address

The root cause was structural: a save inside the element loop cannot preserve a register the loop
clobbers, because the slot moves with the pointer. v30 hoists the save out of the loop into a
dedicated frame slot and uses **the same address expression** on every exit:

```
_batch_entry:  sub sp, sp, #16          <- once per batch, OUTSIDE the element loop
               stur x30, [sp, #8]       <- the caller's link, saved before the first clobber
...
L_loop:        <element body>           <- back-edge target; no save here any more
...
row 14b:       ldur x30, [sp, #8]       <- identical address on every exit
               add sp, sp, #16
               ret
```

Every exit stub carries the same two instructions, so all five exits and the fall-through read the
one save.

**Dynamic proof**: at every executed `ldur x30` the loaded value equals the caller's link captured
before the batch — **1,000/1,000 restores over 1,000 commits (8,000 elements)**. Controls: the v29
code (save per iteration at `[x19,-8]`) produces 3/3 wrong restores, and deleting the hoisted save
produces 3/3 wrong restores. Both defects are caught by the same equality that certifies the fix.

## 3. What the extra instructions cost, and what is unchanged

| row | v29 | v30 | why |
|---|---:|---:|---|
| 6 (transfers) | 19 | 20 | the product transfer `stur x8,[sp,#0]` |
| 10 (predicate) | 1 | 2 | `ldur` the final product, then test it |
| 14 (link) | 2 | 2 | restore + frame teardown (the save moved out of the loop) |
| 15 (batch prologue) | — | 2 | frame setup + one link save, hoisted; charged per element anyway (a 2^18-fold over-charge) |

131 non-refresh / 134 on the refresh step, 64-cycle average **131.046875** against v29's
127.046875 — a step price of **1.065531** vs 1.063523, i.e. +0.0027 in log2. The ledger was
regenerated from the new emitted object's own tables (`schedule.py` census, `span_asm.py`
reconciliation: **132 literal instructions in the `.o`, 126 leaf + 6 prefix, every charged literal
row one machine instruction, in order, same registers**), and the claim is the smallest 4-decimal
value the cost chain accepts at that price. Ceiling rounding, integer ledger sum 131 = the charged
non-refresh total, and `claim >= log2(exact sum)` all hold in-package.

Unchanged and still banked: row 6's operand injection with the C-boundary marker carrying **no
dataflow claim** (accepted in v29), the six-class rho construction, the 0.39 success bound, the
global work meter, and the emitted-object exactness audits.

## Where this leaves the track

v30 is 0.0069 above `5ecbf15` (126.2364, valid, banked), 0.0096 above `67154fda` (126.2337) and
0.0049 above `4993abe9` (126.2384). It is a **more expensive but better-founded** listing: the two
sentences the last verdict refuted are now statements about instructions that exist and values that
were observed. The banked 126.2364 entry is unaffected either way.

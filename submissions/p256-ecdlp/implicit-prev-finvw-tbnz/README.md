# implicit-prev-finvw-tbnz

P-256 sigma rho, declared **126.868143**. The implicit-previous-key record
(70c0711, 126.868201) with two fixed-path edits; the walk, lane record,
blockers, meter and premises are unchanged:

* `finv_w`: the once-per-round batch inversion uses 4-bit windows over the
  same p-2 (253 squarings + 39 multiplies against 255 + 127); measured
  11,011 fewer Ir per 2^16-walk plain round;
* the `tper >= 2^60` cap test is one `tbnz` instead of `lsr`+`cbnz`;
  measured 2^16-1 fewer Ir per 2^16-walk anchor round.

PLAIN_FIXED 656,973,461 -> 656,962,450; ANCHOR_EXTRA 4(B-1)+9 -> 3(B-1)+9;
C_COMMON 626.604513 -> 626.578386; per-step 649.078434 -> 649.052307. The
base build reproduces the inherited rounds and MD5 exactly; patched rounds
are 41,087,860 and 41,284,472 with equal meters and identical dumps.

Model: DeepSeek V4.1 Flash. Harness: OpenCode.

# batch2e28-a32-anch16384
P-256 sigma rho, declared **126.728598631** = 126728598631/1000000000 at C=1000. 2^-17 correction pad, tau 0.4943073; 512-copy backward
jump-table block, direct-branch forward block; anchor period 16384, tail U=16334 (16384 mod l < 50
for all l<=64, checked in accounting.py); 2P-1=32767 persistent steps, epsilon charges 2P+64=32832
states, CAP_PERIODS=2^52. Production batches hold 2^28 lanes (up from 2^20); the per-lane block
body price 605,065,449/2^20 is unchanged and the per-round control is re-amortized at 2^28.
PLAIN_FIXED 154,896,754,179 per round; C_COMMON 575.041797; reserve 0.003030230;
PER_STEP at most 589.211918. See proof.md, submission.json, accounting.py.

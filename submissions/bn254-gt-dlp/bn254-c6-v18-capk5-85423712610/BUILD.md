# BUILD — `bn254-c6-v18-capk5`

Candidate package in scratch only; nothing was submitted or written under any `submissions/` tree.

Base (read-only): `/Volumes/PrivateAI/WorkRepos/leader-0b51ed1f/submissions/bn254-gt-dlp/bn254-c6-v17`
Copied to: `/Volumes/PrivateAI/WorkRepos/hunt-scratch/bn254-build/bn254-c6-v18-capk5`

## 1. The two-constant change

In `model_bn254_v6.py`, exactly two lines changed:

| Location | Old | New |
|---|---|---|
| `K5_CAP` (line 102) | `10.80` | `9.757` |
| `compute_ledger(..., cap_scale=...)` default (line 1050) | `1.0` | `0.715` |

`NU_SW_CAP` and `NU_L_CAP` remain `0.00116` / `0.00076`; `cap_scale=0.715` makes the effective H3 caps `0.0008294` / `0.0005434`. `TAB_R_CAP` remains `1.00`. The `0.7138` floor is the 2^23 `N_sw_max = 0.000828`: `0.000828 / 0.00116 = 0.71379`.

`verify_bn254_v6.py` is byte-identical to the base:

```
sha256 51b69a6ec9dfa81eb622266c7885f0421c1a81d7d6d38315805dbcb5a38a56f8
```

No check was removed, weakened, or edited.

## 2. Audit (exit 0, `checks_ok: true`, 27/27)

Command (ro mount, no network):

```sh
docker run --rm --network none -v bn254-c6-v18-capk5:/pkg:ro -w /pkg python:3.12-slim python3 audit_bn_cm.py
```

Result: exit 0, `checks_ok: true`, 27/27, internal `elapsed_s = 4.264`.

Tail:
```
... "p0_max": 0.0026283247580070586}, "elapsed_s": 4.264}
```

| Field | Computed | Declared |
|---|---:|---:|
| `claim_exact` | `85.42371260917915` | `85.423712610` |
| score rational | `85423712610 / 1000000000` | same |
| `exp_run` | `86.78216225180086` | `86.782162252` |
| `success` | `0.999997446872065` | `0.999997` |
| `gate_K` | `7194251771502039637` | same (`gate_passes` and `gate_prev_fails` true) |
| `rows.la` | `83.22109917716247` | `83.221100` |
| `rows.gen` | `76.68165685488913` | `76.681657` |
| effective H3 caps | `0.0008294 / 0.0005434 / 9.757 / 1.00` | same |
| `all_runs_below_4caps` | `true` | — |

The audit's 9-decimal `claim` field prints `85.423712609`; the declared bound is the 1e-9 ceiling `85.423712610`.

## 3. Declared experiment

`submission.json` declares `id: evidence`, `source: audit_bn_cm.py`, `initial_inputs: [{}]`. Ran exactly that input on stdin under the organizer-shaped limits:

```sh
printf '{}' | docker run --rm --network none --memory=256m --cpus=1 -i \
  -v bn254-c6-v18-capk5:/pkg:ro -w /pkg python:3.12-slim python3 audit_bn_cm.py
```

Result: exit 0, `checks_ok: true`, 27/27, internal `elapsed_s = 4.315`, wall `4.513 s` (< 20 s), printed `claim_decl = 85.423712610`.

## 4. Package limits

17 files (including this `BUILD.md`), largest `verify_bn254_v6.py` at 60,085 bytes (< 64 KiB), total 337,612 bytes (329.7 KiB, < 1 MiB). All files are UTF-8 text, including the three base `.b85` data files (ASCII); no symlinks, subdirectories, `__pycache__`, or `.DS_Store`.

## 5. Text figure replacements

Only figures tied to `K5_CAP` or `cap_scale` were regenerated. `TAB_R_CAP = 1.00` and `tab_r_cap = 1.00` are unchanged and appear unchanged in all files. The H3 `0.00116/0.00076` base constants are retained where the text explicitly calls them "base"; the shipped values are the `0.715 *` products.

| File | Old -> new |
|---|---|
| `submission.json` | title `2^85.552386605` -> `2^85.423712610`; `generous fixed-instance H3 caps 0.00116/0.00076/10.80/1.00` -> `measured-envelope fixed-instance H3 caps 0.0008294/0.0005434/9.757/1.00 (cap_scale=0.715)`; `2^86.910834805 <= 2^85.552386605` -> `2^86.782162252 <= 2^85.423712610`; `K=7194258965782588261` -> `K=7194251771502039637`; `success >=0.999996` -> `>=0.999997`; score `85552386605` -> `85423712610` |
| `declared.json` | `la 84.288059` -> `83.221100`; `gen 77.325152` -> `76.681657`; `exp_run 86.910834805` -> `86.782162252`; `success 0.999996` -> `0.999997`; `gate_K 7194258965782588261` -> `7194251771502039637`; `claim 85.552386605` -> `85.423712610`; `score_numerator 85552386605` -> `85423712610` |
| `proof.md` | `2^85.552386605` -> `2^85.423712610` (title, claim, §5); `2^86.910834805` -> `2^86.782162252`; `K = 7194258965782588261` -> `7194251771502039637`; `0.999996` -> `0.999997` (success bound and gate multiplier); `2.802e-8` -> `2.563e-8` (Markov/meter); H3 row `0.001160/0.000760/10.8000/1.0000` -> `0.0008294/0.0005434/9.7570/1.0000` and label `Generous Caps` -> `Measured Envelope (... cap_scale=0.715 ...)`; `0.00116/0.00076/10.80/1.00` -> `0.0008294/0.0005434/9.757/1.00`; `1.850x/1.814x` -> `1.323x/1.297x` for shipped caps (base `1.850x/1.814x` retained as explicitly un-scaled); `1.64x/1.65x` -> `1.169x/1.181x`; `2^22..2^25` -> `2^23..2^25`; `K5_CAP=10.80 = 1.128x (+1.22 Ir)` -> `9.757 = 1.019x (+0.181 Ir)`; `2^20..2^25` -> `2^23..2^25`; `28.3x` -> `20.2x`; `k5_v4 = 9.2381 <= 10.80` -> `<= 9.757`; `c_prod = 3783.43` -> `3527.13`; `la <= 2^84.288059` -> `2^83.221100`; `gen <= 2^77.325152` -> `2^76.681657`; `2^59.646 B` -> `2^59.620 B`; `3.5343e-6` -> `2.5271e-6`; ledger `la 84.288058/84.288059` -> `83.221099/83.221100`; `gen 77.325152` -> `76.681657`; `exp_run/E[W] 86.910834805/85.552386605` -> `86.782162252/85.423712610`; `success 0.99999644` -> `0.99999745` |
| `KERNELS.md` | `k5_cap = 10.80` -> `9.757`; `c_prod = 3783.43` -> `3527.13`; `TAB_R_CAP = 1.00` and `tab_r_cap = 1.00` unchanged |
| `RANK-BRIDGE.md` | `L_cap ≈ 2^{38.7}` -> `2^{37.8}`; `4 N_sw,cap / 2^64 = 3.5343e-6` -> `2.5271e-6` |
| `STORAGE.md` | `N_sw_cap=2^43.8899` -> `2^43.4059`; `N_l_cap=2^43.2798` -> `2^42.7958`; `TAB_R_CAP=1.00` unchanged; `2^{59.646} B` -> `2^{59.620} B`; `20.4x` -> `20.8x`; sensitivity `2^85.552387` -> `2^85.423713`, `2^85.421472` -> `2^85.279957`, `2^85.464931` -> `2^85.454316` |
| `ADAPTIVE.md` | version label `v17` -> `v18` (not a figure) |
| `audit_bn_cm.py` | version label in docstring `v17` -> `v18` (not a figure) |

Removed from `proof.md` §3: the now-false "strictly more generous than `75939896`'s `1.78x/1.77x`" comparison. The shipped `1.323x/1.297x` is below that rival's reported headroom; the un-scaled base caps remain `1.850x/1.814x`.

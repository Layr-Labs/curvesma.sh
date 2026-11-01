"""Focused hosted checks. Full local reproduction is evidence16.py."""

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE / "files" if (HERE / "files").is_dir() else HERE
sys.path.insert(0, str(ROOT))

import evidence16 as E


def fields() -> dict:
    """Recompute field identities and kernel counts at production parameters."""
    field_data = E.check_fields()
    la_count, la_branches = E.check_la_kernel()
    full_count, full_branches = E.check_la_kernel("bw_dot.s", "LBB1_2:")
    checks = {
        "resultant": E.check_resultant(),
        "tower": E.check_tower(),
        "rho_semismooth": E.check_rho_and_semismooth(),
        "fields": field_data["unramified"]
        and field_data["sig"]["odd"] == E.M.N_SM,
        "la_loop": la_count == 33 and la_branches == ["b.ne"],
        "full_field_loop": full_count == 97
        and full_count <= E.M.C_FULL
        and full_branches == ["b.ne"],
    }
    return {"checks": checks, "field_data": field_data}


def pool(seed: int) -> dict:
    """Regenerate eight pool rows and test missing-log reconstruction."""
    cfg = json.loads((ROOT / "config.json").read_text())
    header, rows = E.M.load_pool(
        (ROOT / "pool12-1.b85").read_text()
        + (ROOT / "pool12-2.b85").read_text()
    )
    pool_ok, canonical_sample = E.verify_pool(cfg, rows, 8, E.EB.root_disks())
    reconstruction = E.filtsim.missing_recon(seed)
    checks = {
        "regenerated_pool_rows": pool_ok,
        "missing_logs": all(row[3] <= row[4] for row in reconstruction),
    }
    return {
        "checks": checks,
        "stored_draws": sum(header["samples"].values()),
        "stored_canonical_rows": len(rows),
        "regenerated_rows": 8,
        "regenerated_canonical_rows": canonical_sample,
        "missing_log_rows": reconstruction,
    }


def main() -> None:
    request = json.load(sys.stdin)
    if set(request) - {"scope", "seed"}:
        raise ValueError("Unknown input")
    scope = request.get("scope")
    seed = request.get("seed", 1)
    if type(seed) is not int or not 1 <= seed <= 10:
        raise ValueError("seed must be an integer from 1 through 10")
    if scope == "fields":
        result = fields()
    elif scope == "pool":
        result = pool(seed)
    else:
        raise ValueError("scope must be fields or pool")
    result.update(
        {
            "scope": scope,
            "full_diagnostic": False,
            "recurrence_recomputed": False,
            "checks_ok": all(result["checks"].values()),
        }
    )
    print(json.dumps(result))
    if not result["checks_ok"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

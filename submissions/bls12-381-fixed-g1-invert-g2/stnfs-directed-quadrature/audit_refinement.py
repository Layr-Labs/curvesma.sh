"""Bounded accounting and integration audit; not the full recurrence."""

import hashlib
import json
import math
import sys
from decimal import Decimal, localcontext
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE / "files" if (HERE / "files").is_dir() else HERE
sys.path.insert(0, str(ROOT))

import descent16 as D
import model16 as M
import retain16 as R


def main() -> None:
    """Check recorded observations, controls, physical windows and bounds."""
    request = json.load(sys.stdin)
    if request != {}:
        raise ValueError("The audit takes an empty object")
    certificate = json.loads((ROOT / "refinement.json").read_text())
    checks = {}
    checks["code_identity"] = all(
        hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest
        for name, digest in certificate["code_sha256"].items()
    )
    checks["parameters"] = (
        D.P.MESH == 384
        and D.RSTEP == 1 / 1024
        and D.P.TRUNC == 0
        and M.DESCENT_CAP == certificate["descent_cap"]
    )
    checks["windows"] = all(
        D.PLAN16[(b, side)]
        == (old[index][0], 8 * old[index][1], *old[index][2:])
        for b, old in D._PLAN.items()
        for index, side in enumerate("fg")
    ) and D.BOOT16 == (215, 344, 15, 4)
    bad = D.P()
    bad.MESH = 48
    try:
        D.evaluate(bad, D.PLAN16, D.BOOT16)
    except ValueError as error:
        checks["parameter_mismatch_rejected"] = "must agree" in str(error)
    else:
        checks["parameter_mismatch_rejected"] = False
    probes = [
        ("kA0", (600, 69.0, 75.0, 24, 0, 0, 48, 1.0, 0.0)),
        ("kA1", (600, 69.0, 75.0, 42, 1, 0, 48, 1.0, 0.0)),
        ("kA4", (600, 69.0, 100.0, 24, 4, 1, 48, 1.0, 0.0)),
        ("kA15", (1200, 69.0, 215.0, 43, 15, 0, 48, 1.0, 0.0)),
    ]
    controls = {
        row["name"]: row
        for row in certificate["probe_checks"]
        if row["rho_den"] == 1024
    }
    observed = {}
    for name, args in probes:
        lower = D.ss_lower2(*args)
        upper = R.ss_upper2(*args[:-1])
        control = controls[name]
        observed[name] = {"lower": lower, "upper": upper, "expected_lower": control["limited"], "expected_upper": control["upper_limited"]}
        checks[name] = (
            lower > 0
            and lower <= upper
            and math.isclose(lower, control["limited"], rel_tol=1e-12)
            and math.isclose(upper, control["upper_limited"], rel_tol=1e-12)
        )
    with localcontext() as context:
        context.prec = 70
        two = Decimal(2)
        total = sum(two ** Decimal(str(x)) for x in certificate["ledger_upper_rows"])
        log_total = total.ln() / two.ln()
        markov = two ** (
            Decimal(str(certificate["two_log_upper"]))
            - Decimal(str(certificate["descent_cap"]))
        )
        order = Decimal(
            "52435875175126190479447740508185965837690552500527637822603658699938581184513"
        )
        success = (1 - two ** -192) * (
            1 - 2 * Decimal("0.0001983") - markov
            - two ** Decimal("-9.64") - two ** Decimal("-10.70")
            - 2 / order
        )
        checks["mean_rounding"] = (
            certificate["two_log_measured"] < certificate["two_log_upper"]
        )
        checks["ledger"] = log_total < Decimal("121.99951") < Decimal("121.9996") < 122
        checks["success"] = success > Decimal("0.391663") > Decimal("0.39")
    print(json.dumps({
        "checks_ok": all(checks.values()),
        "checks": checks,
        "observed_controls": observed,
        "ledger_upper_log2": str(log_total),
        "markov_upper": str(markov),
        "success_lower": str(success),
        "scope": "Recorded full recurrence plus selected integration controls; not a full recurrence rerun",
    }))
    if not all(checks.values()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()

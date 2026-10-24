#!/usr/bin/env python3
"""Reproduce every measurement in this package.

    python3 experiments/runner.py          # full: theorem, cycle, census, audit
    python3 experiments/runner.py --quick  # census at n=2000 instead of 30000

Writes evidence/v27-measure-*.json, certificate.json, experiments/results.json.  Exit code 0
iff audit_v27.py passes (all checks plus the planted controls).
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)


def run(cmd, env=None):
    e = dict(os.environ)
    if env:
        e.update(env)
    print('$', ' '.join(cmd))
    r = subprocess.run(cmd, cwd=PKG, env=e, capture_output=True, text=True)
    sys.stdout.write(r.stdout)
    if r.returncode != 0:
        sys.stderr.write(r.stderr)
    return r.returncode


def main():
    quick = '--quick' in sys.argv
    n = '2000' if quick else '30000'
    rc = 0
    rc |= run([sys.executable, 'v27_measure.py', 'theorem'])
    rc |= run([sys.executable, 'v27_measure.py', 'cycle'])
    rc |= run([sys.executable, 'v27_measure.py', 'model'])
    rc |= run([sys.executable, 'v27_measure.py', 'recovery'])
    rc |= run([sys.executable, 'run_census_v27.py'], env=dict(TRIALS=n))
    rc |= run([sys.executable, 'audit_v27.py'])
    return rc


if __name__ == '__main__':
    sys.exit(main())

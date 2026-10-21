#!/usr/bin/env python3
"""Reference instruction counts derived from the recorded Callgrind measurements."""
import argparse
import hashlib
import json
from pathlib import Path

DIRECTORY = Path(__file__).resolve().parent / 'reference_costs'


def reference_costs(results=None):
    record = json.loads((DIRECTORY / 'measurement.json').read_text())
    for name, expected in record['harness_sha256'].items():
        if hashlib.sha256((DIRECTORY / name).read_bytes()).hexdigest() != expected:
            raise ValueError('Measurement harness differs from the recorded source: ' + name)
    counts = record['instruction_totals']
    sizes = record['iterations']
    if results is not None:
        counts = {}
        for name in ('field', 'fp12', 'bls_field'):
            counts[name] = []
            for n in sizes:
                lines = (Path(results) / f'{name}-{n}.callgrind').read_text().splitlines()
                counts[name].append(int(next(line for line in lines if line.startswith('summary: ')).split()[1]))
        if counts != record['instruction_totals']:
            raise ValueError('Counts differ from the pinned calibration; do not change prices automatically.')
    delta = sizes[1] - sizes[0]
    slopes = {}
    for name in ('field', 'fp12', 'bls_field'):
        work = counts[name][1] - counts[name][0]
        if work <= 0 or work % delta:
            raise ValueError('Nonintegral or invalid instruction-count slope')
        slopes[name] = work // delta
    measured = {'secp256k1-ecdlp': 6 * slopes['field'],
                'bls12-381-invert-g1': slopes['fp12'],
                'bls12-381-g1-ecdlp': 6 * slopes['bls_field'],
                'bls12-381-miller-inversion': slopes['fp12']}
    if measured != record['measured_reference_work']:
        raise ValueError('Measured work does not match the calibration record')
    # The added fixed-argument track reuses the existing Fp12 calibration.
    # Preserve the historical measurement record, but expose active tracks only.
    prices = dict(record['fixed_reference_prices'])
    prices.pop('bls12-381-miller-inversion')
    prices['bls12-381-fixed-g1-invert-g2'] = prices['bls12-381-invert-g1']
    prices['bls12-381-fixed-g2-invert-g1'] = prices['bls12-381-invert-g1']
    prices['bls12-381-g1-ecdlp-plain'] = prices['bls12-381-g1-ecdlp']
    # Borrowed normalization only; no P-256 measurement is claimed.
    prices['p256-ecdlp'] = prices['secp256k1-ecdlp']
    prices.pop('bls12-381-invert-g1')
    return prices


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, help='Check a fresh reproduction against the recorded totals')
    args = parser.parse_args()
    print(json.dumps(reference_costs(args.results), indent=2))

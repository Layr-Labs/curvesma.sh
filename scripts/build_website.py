#!/usr/bin/env python3
"""Export the protected challenge definitions alongside the static website."""
import json
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT / 'website' / 'dist'


def build():
    manifest = json.loads((ROOT / 'benchmark.json').read_text())
    targets = SITE / 'targets'
    targets.mkdir(parents=True, exist_ok=True)
    for track in manifest['tracks']:
        name = track['name']
        shutil.copyfile(ROOT / 'tracks' / name / 'target-profile.json',
                        targets / f'{name}.json')
    shutil.copyfile(ROOT / 'tracks' / 'classical-group-work-v2.json',
                    SITE / 'cost-model.json')
    cost = (ROOT / 'tracks' / 'docs' / 'COST_MODEL.md').read_text()
    source = 'https://github.com/Layr-Labs/curvesma.sh'
    cost = cost.replace('](../../scripts/', f']({source}/tree/master/scripts/')
    cost = cost.replace('](RESCORING.md)',
                        f']({source}/blob/master/tracks/docs/RESCORING.md)')
    (SITE / 'cost-model.md').write_text(cost)


if __name__ == '__main__':
    build()

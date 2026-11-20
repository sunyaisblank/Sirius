"""Diagnose the first fp64 metric submission in an otherwise idle qualification build."""

from pathlib import Path
import shutil
import sys

import run_qualification_cases as q

root = Path.cwd()
q.WORK = root / 'out/strict-linux-science/isolated-fp64-metric'
q.RENDERS = root / 'renders/strict-linux-science/isolated-fp64-metric'
case = 'KernelParity.MetricConnectionAndJetShareRepresentedParametersFp64'


def plan(model):
    governed = {entry['name'] for row in model['acceptance_criteria'] for entry in row['evidence']}
    q.require(case in governed, 'Expected the original governed fp64 metric case')
    return {'fp64-metric-isolated': ({case}, True, 900)}


q.plan = plan
if __name__ == '__main__':
    q.require(shutil.disk_usage('/').free >= 100_000_000_000, 'Less than 100 GB available on /')
    for parent in (q.WORK.parent, q.WORK.parent/'p2-supplement'):
        for receipt in (parent/'runs').glob('*/owner.json'):
            record = q.read(receipt)
            q.require(record['status'] == 'terminal' and not record['remaining_owned_processes'],
                      f'Previous run remains live: {receipt}')
    q.WORK.mkdir()
    sys.exit(0 if q.run() else 1)

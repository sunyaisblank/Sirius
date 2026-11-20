"""Fill the original beam and bundle roster on the same strict qualification build."""

from pathlib import Path
import shutil
import sys

import run_qualification_cases as q

root = Path.cwd()
q.WORK = root / 'out/strict-linux-science/p2-supplement'
q.RENDERS = root / 'renders/strict-linux-science/p2-supplement'
host_report = root / 'out/strict-linux-science/host-supplement/results.json'
prior = root / 'attestations/qualification/56f79941/integrated-science-2026-11-21/MANIFEST.json'


def plan(model):
    row = next(item for item in model['acceptance_criteria'] if item['id'] == 'P2')
    names = {item['name'] for item in row['evidence']}
    host = q.read(host_report)['cohorts']['p2-host-geometry-and-transaction']
    q.require(host.get('accepted'), 'The host P2 roster must have passed')
    measured = {case['name'] for case in host['cases']}
    previous = set(q.read(prior)['operating_model_case_coverage']['P2']['unmeasured_cases'])
    static = {name for name in names if name.startswith(('RetainedArithmetic.', 'RetainedValue.'))}
    coupled_host = (names - previous) | static
    device_names = names - measured - coupled_host
    probes = {name for name in device_names if name.startswith('Kernel')}
    retained = {name for name in device_names if name.startswith('Retained')}
    live = device_names - probes - retained
    q.require(len(names) == 126 and len(measured) == 52 and len(coupled_host) == 33 and
              len(probes) == 20 and len(retained) == 19 and len(live) == 2,
              'The exact CPU and device P2 roster must be covered once')
    return {
        'p2-coupled-host-oracles': (coupled_host, False, 1800),
        'p2-device-reference-probes': (probes, True, 14400),
        'p2-retained-intermediate': (retained, True, 14400),
        'p2-live-backend-agreement': (live, True, 3600),
    }


q.plan = plan
if __name__ == '__main__':
    q.require(shutil.disk_usage('/').free >= 100_000_000_000, 'Less than 100 GB available on /')
    q.WORK.mkdir()
    sys.exit(0 if q.run() else 1)

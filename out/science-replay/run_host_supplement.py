"""Execute static precision admissions and the remaining host P2 witnesses."""

from pathlib import Path
import shutil
import sys

import run_qualification_cases as q

root = Path.cwd()
previous = root / 'attestations/qualification/56f79941/integrated-science-2026-11-21/MANIFEST.json'
q.WORK = root / 'out/strict-linux-science/host-supplement'
q.RENDERS = root / 'renders/strict-linux-science/host-supplement'
original_plan = q.plan


def plan(model):
    initial = original_plan(model)
    admission = initial['p1-precision-admission'][0]
    p2 = next(row for row in model['acceptance_criteria'] if row['id'] == 'P2')
    missing = set(q.read(previous)['operating_model_case_coverage']['P2']['unmeasured_cases'])
    missing &= {item['name'] for item in p2['evidence']}
    host_prefixes = ('CameraContinuousProjection.', 'CpuTraceBoundary.', 'DopriPositionSegment.',
                     'KerrInfinity.', 'PointSourceDetector.', 'RayBundleTest.', 'RenderSessionProbe.',
                     'SourceSkyDifferential.', 'SourceSkyMap.', 'StarfieldGeneratorTests.',
                     'StarfieldPointTest.')
    device = 'PointSourceDetector.RetainedRotatingIsolatedStarMatchesIndependentSharedRegionBrightness'
    host = {case for case in missing if case.startswith(host_prefixes) and case != device}
    q.require(len(admission) == 3 and len(host) == 52, 'Host science roster drifted')
    return {'p1-precision-static': (admission, False, 300),
            'p2-host-geometry-and-transaction': (host, False, 1800)}


q.plan = plan
if __name__ == '__main__':
    q.require(shutil.disk_usage('/').free >= 100_000_000_000, 'Less than 100 GB available on /')
    q.WORK.mkdir()
    sys.exit(0 if q.run() else 1)

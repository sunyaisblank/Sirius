"""Run isolated original backend cases from an exact strict build."""

import argparse
import os
from pathlib import Path
import shutil
import sys

ROOT = Path.cwd()
sys.path.insert(0, str(ROOT/'out/science-replay'))
import run_qualification_cases as q

CASES = {
    'original-cold': ('OperationalVulkan.RequiredRuntimeDispatch', True, 330),
    'kerr-geometry': ('VulkanRenderSession.CpuVulkanAgreeOnKerrGeometryWithinStatisticalBounds',
                      True, 9000),
    'ellis-geometry': ('VulkanRenderSession.CpuVulkanAgreeOnMorrisThorneGeometryWithinStatisticalBounds',
                       True, 9000),
    'p1-retained-fp32': ('ShadowBoundary.RetainedFp32ProductsKerrNearExtremalBardeenBoundaryAt1080p',
                          True, 28800),
    'p1-retained-fp64': ('ShadowBoundary.RetainedFp64ProductsKerrNearExtremalBardeenBoundaryAt1080p',
                          True, 28800),
    'p1-burn-in': ('OperationalP1.NearExtremalBurnIn', True, 28860),
    'frame-transaction': ('VulkanRenderSession.ContinuationRendererPublishesOnlyCompleteFramesWithinActualBudget',
                           True, 3600),
    'progressive-viewer': ('ViewCommandOperational.VulkanRefinementPublishesProgressiveFrames',
                           True, 21660),
    'lens-live': ('VulkanRenderSession.NonSquareMultisamplingCameraAndLensReachLiveKernel',
                  True, 7200),
}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=CASES)
    parser.add_argument('--run', default='ordinary', choices=('ordinary','certified'))
    options = parser.parse_args()
    q.require(shutil.disk_usage('/').free >= 100_000_000_000, 'Less than 100 GB remains on /')
    q.require('SIRIUS_CERTIFIED_SPHERICAL_ESCAPE' not in os.environ,
              'The replay must not inherit a private continuation setting')
    q.WORK = ROOT/'out/strict-backend-science'/options.run/options.phase
    q.RENDERS = ROOT/'renders/strict-backend-science'/options.run/options.phase
    name, vulkan, seconds = CASES[options.phase]
    q.plan = lambda model: {options.phase: ({name}, vulkan, seconds)}
    q.WORK.mkdir(parents=True)
    q.RENDERS.parent.mkdir(parents=True, exist_ok=True)
    q.write(q.WORK/'host-admission.json', {
        'affinity':sorted(os.sched_getaffinity(0)),
        'original_ctest_case':name,
        'outer_owner_watchdog_seconds':seconds,
        'note':'The existing CTest and numerical/dispatch bounds remain unchanged.',
    })
    sys.exit(0 if q.run() else 1)

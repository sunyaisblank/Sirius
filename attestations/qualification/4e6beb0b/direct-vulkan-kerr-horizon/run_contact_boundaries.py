"""Recheck the positive-interval chart against original coupled and portability witnesses."""

from pathlib import Path
import sys

ROOT = Path.cwd()
sys.path.insert(0, str(ROOT/'out/science-replay'))
import run_qualification_cases as q

cases = {
    'coupling': {
        'KernelTrace.ActualTraceClipsJacobiAndRetainsInvalidSamplesAcrossRungs',
        'KernelTrace.ActualTraceFiniteSphereExcludesLaterDiskAndVolumeAcrossRungs',
    },
    'emission': {
        'KernelPortability.CudaEmissionCarriesTheNativeComputeEntryPoint',
        'KernelPortability.MetalEmissionCarriesTheNativeComputeEntryPoint',
    },
}

if __name__ == '__main__':
    group, = sys.argv[1:]
    q.WORK = ROOT/'out/kerr-outgoing-trace'/(group+'-strict')
    q.RENDERS = ROOT/'renders/kerr-outgoing-trace'/(group+'-strict')
    q.plan = lambda model: {group: (cases[group], group=='coupling', 1800 if group=='coupling' else 120)}
    q.WORK.mkdir()
    q.RENDERS.parent.mkdir(parents=True, exist_ok=True)
    sys.exit(0 if q.run() else 1)

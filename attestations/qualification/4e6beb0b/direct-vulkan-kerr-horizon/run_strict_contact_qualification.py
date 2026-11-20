"""Execute the unchanged direct Vulkan witnesses for the past-horizon chart."""

from pathlib import Path
import shutil
import sys

ROOT = Path.cwd()
sys.path.insert(0, str(ROOT/'out/science-replay'))
import run_qualification_cases as q

cases = {
    'bardeen': {'VulkanRenderSession.KerrNearExtremalBardeenBoundaryAt1080p'},
    'trace': {
        'KernelTrace.KerrRenderIsFiniteNonConstantWithBoundedShadow',
        'KernelTrace.Fp64RungAgreesWithFp32OnKerrScene',
        'KernelTrace.CompensatedRungTracksFp64AtLeastAsWellAsFp32',
        'KernelTrace.ActualTraceClipsJacobiAndRetainsInvalidSamplesAcrossRungs',
        'KernelTrace.ActualTraceFiniteSphereExcludesLaterDiskAndVolumeAcrossRungs',
        'VulkanRenderSession.CpuVulkanAgreeOnKerrGeometryWithinStatisticalBounds',
        'VulkanRenderSession.CpuVulkanAgreeOnMorrisThorneGeometryWithinStatisticalBounds',
        'VulkanRenderSession.ZeroActiveTracePreservesRadianceAcrossPrecisionRungs',
    },
    'policy': {
        'OperationalEvidence.SourceAndIdealGovernanceRejectUncoveredClaims',
        'OperationalBuildPolicy.ReleaseCannotWeakenGates',
        'OperationalBuildGate.ReleaseInstallRequiresExactPassedEstate',
    },
}

if __name__ == '__main__':
    group, = sys.argv[1:]
    q.require(shutil.disk_usage('/').free >= 100_000_000_000,
              'Less than 100 GB remains on /')
    q.WORK = ROOT/'out/kerr-outgoing-trace'/(group + '-strict')
    q.RENDERS = ROOT/'renders/kerr-outgoing-trace'/(group + '-strict')
    q.plan = lambda model: {group: (cases[group], group != 'policy',
                                    3600 if group != 'policy' else 300)}
    q.WORK.mkdir(exist_ok=False)
    q.RENDERS.parent.mkdir(parents=True, exist_ok=True)
    sys.exit(0 if q.run() else 1)

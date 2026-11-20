"""Execute original P1 cases and the local camera, catalogue and governor witnesses."""

from pathlib import Path
import shutil
import sys

ROOT = Path.cwd()
sys.path.insert(0, str(ROOT/'out/science-replay'))
import run_qualification_cases as q

WORK = ROOT/'out/strict-p1-camera'
RENDERS = ROOT/'renders/strict-p1-camera'
prior_manifest = ROOT/'attestations/qualification/55db0d1a/strict-linux-catalogue-and-disk/MANIFEST.json'


def selection(model, phase):
    criteria = {row['id']: {entry['name'] for entry in row['evidence']}
                for row in model['acceptance_criteria']}
    p1 = criteria['P1']
    p1_kernel = {name for name in p1 if name.startswith('KernelParity.')}
    p1_artifact = {'KernelParity.PrecisionProbeArtifactsCarryOnlyTheirDeclaredFloat64Capability'}
    p1_readbacks = {'KernelParity.PrecisionRungsConserveNearExtremalKerrWithoutImageComparison'}
    p1_boundary = {name for name in p1 if name.startswith(('ShadowBoundary.Retained',
                                                            'VulkanRenderSession.'))}
    p1_retained = {name for name in p1 if name.startswith('FullPathAcceptance.Vulkan')}
    p1_repeat = {'OperationalP1.NearExtremalBurnIn'}
    p1_host = p1 - p1_kernel - p1_boundary - p1_retained - p1_repeat
    q.require(len(p1)==41 and len(p1_host)==19 and len(p1_boundary)==3 and len(p1_retained)==1,
              'The original P1 roster must match the source model')

    if phase == 'host':
        p3p5 = criteria['P3'] | criteria['P5']
        geometry = {name for name in p3p5
                    if not name.startswith(('VulkanRenderSession.', 'RetainedCameraProgram.',
                                            'KernelParity.', 'OperationalAttestation.'))}
        e3 = criteria['E3']
        host_governors = {name for name in e3
                          if name.startswith(('DispatchGovernor.', 'RetainedComputeAdmission.',
                                              'RenderSessionProbe.BackendAuto'))}
        host_governors.add('VulkanBackend.Fma32DeviceAdmissionRequiresAdvertisedFeatureAndEveryFloatControl')
        q.require(host_governors <= e3 and len(geometry)==62 and len(host_governors)==19,
                  'Unexpected camera or governor roster')
        return {
            'p1-host-geodesic-accuracy': (p1_host | p1_artifact, False, 1800),
            'host-camera-and-catalogue': (geometry, False, 1800),
            'host-governors-and-cpu-routing': (host_governors, False, 1800),
        }
    if phase == 'device':
        legacy = {name for name in p1_boundary if name.startswith('VulkanRenderSession.')}
        direct_camera = {name for name in criteria['P3'] | criteria['P5']
                         if name.startswith(('RetainedCameraProgram.', 'KernelParity.'))}
        governor = {'VulkanBackend.BufferAllocationLimitCountsActualResidencyAndPreservesExistingBuffers',
                    'RetainedComputeTest.SoftwareRendererPreparationPreservesBuffersAndPhysicalFeedback'}
        q.require(len(direct_camera)==5 and governor <= criteria['E3'],
                  'Expected the original device optics and governor controls')
        return {
            'p1-device-precision-and-rollback': (p1_kernel - p1_readbacks - p1_artifact, True, 7200),
            'p1-three-rungs': (p1_readbacks, True, 3600),
            'p1-retained-conservation': (p1_retained, True, 3600),
            'p1-legacy-bardeen': (legacy, True, 3600),
            'device-optics-and-catalogue': (direct_camera, True, 3600),
            'device-allocation-and-preparation': (governor, True, 1800),
        }
    if phase == 'shadow':
        retained = {name for name in p1_boundary if name.startswith('ShadowBoundary.')}
        return {('p1-shadow-fp64' if 'Fp64' in name else 'p1-shadow-fp32'):
                ({name}, True, 28800) for name in sorted(retained)}
    if phase == 'burn-in':
        return {'p1-original-three-repeat-burn-in': (p1_repeat, True, 28860)}
    raise ValueError(f'Unrecognized phase {phase}')


if __name__ == '__main__':
    phase, = sys.argv[1:]
    q.require(shutil.disk_usage('/').free>=100_000_000_000, 'Less than 100 GB available on /')
    prior = q.read(prior_manifest)
    q.require(prior['model']==q.pin(q.MODEL), 'Operating model differs from prior exact proof')
    q.WORK = WORK/phase
    q.RENDERS = RENDERS/phase
    q.plan = lambda model: selection(model, phase)
    q.WORK.mkdir()
    q.RENDERS.parent.mkdir(exist_ok=True)
    sys.exit(0 if q.run() else 1)

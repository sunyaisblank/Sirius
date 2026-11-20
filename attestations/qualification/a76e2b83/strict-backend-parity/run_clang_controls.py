"""Verify the exact Clang build against the original qualification governance controls."""

from pathlib import Path
import sys

ROOT=Path.cwd()
sys.path.insert(0,str(ROOT/'out/science-replay'))
import run_qualification_cases as q

NAMES={
    'OperationalEvidence.SourceAndIdealGovernanceRejectUncoveredClaims',
    'OperationalBuildPolicy.ReleaseCannotWeakenGates',
    'OperationalRetainedShaders.OptionalCapabilitiesPreserveFloatControls',
    'OperationalBuildGate.ReleaseInstallRequiresExactPassedEstate',
    'OperationalAttestation.FalseExternalEvidenceIsRejected',
    'OperationalAttestation.PreflightAndNativeRuntimeRejectFalseHostDevice',
    'OperationalAlignment.IncompleteOrAmbiguousAttestationSetIsRejected',
}
q.BUILD=ROOT/'bin/linux-clang'
q.CTEST=['/usr/bin/ctest','--test-dir',str(q.BUILD),'-C','Release']
q.WORK=ROOT/'out/strict-backend-science/ordinary/clang-controls'
q.RENDERS=ROOT/'renders/strict-backend-science/ordinary/clang-controls'
q.plan=lambda model: {'clang-controls':(NAMES,False,300)}
q.WORK.mkdir()
sys.exit(0 if q.run() else 1)

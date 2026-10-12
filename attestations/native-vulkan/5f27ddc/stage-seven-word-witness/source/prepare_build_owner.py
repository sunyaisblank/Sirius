"""Prepare unchanged bounded local build apparatus for the diagnostic head."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib,json,shutil,subprocess
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
source=json.loads((WORK/'source.json').read_text());head=source['source_revision']
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==head
assert subprocess.check_output(['git','rev-parse','@{u}'],text=True).strip()==head
assert not subprocess.check_output(['git','status','--porcelain'],text=True)
assert json.loads((WORK/'host-source-review.json').read_text())['accepted']
folder=WORK/'build-owner';assert not folder.exists();folder.mkdir()
expected={'run_owned.py':'b47e4bb22e63c16244802a5d121160a63b942f5c37614d26ded319fe4aaea57d','owned_guard.py':'747fe1e0bcaa0cad8fb867f22a30935b1c367dc0e2c0ad8204a87f6593899056'}
for name,digest in expected.items():
    path=WORK/name;assert hashlib.sha256(path.read_bytes()).hexdigest()==digest
    shutil.copyfile(path,folder/name)
frozen=[str(p.relative_to(ROOT)) for p in WORK.iterdir() if p.is_file()]
frozen += [str(p.relative_to(ROOT)) for p in (WORK/'kernels').iterdir() if p.is_file()]
frozen += [str(WORK.relative_to(ROOT)/'prepare_build_owner.py')]
limits={'seconds':300,'rss_kib':4194304,'reserve_kib':2097152}
scope='Existing300s4GiB2GiBreserve owned local source/build checks only. Strict original compiler/governance/alignment controls; diagnostic candidate0/8 stays non-promoting. No device/scientific/runtime/performance qualification.'
actions={}
for name,argv in [
    ('configure',['cmake','--preset','linux-gcc','-DSIRIUS_ALIGNMENT_MODE=qualification','-DSIRIUS_ATTESTATION_ROOT=']),
    ('build',['cmake','--build','--preset','linux-gcc','--target','SiriusAlignmentGate','SiriusSourceGovernance','sirius','sirius_backend_tests','sirius_app_tests','-j2']),
    ('authority',['ctest','--test-dir','bin/linux-gcc','--output-on-failure','--no-tests=error','-R','^(OperationalEvidence\\.SourceAndIdealGovernanceRejectUncoveredClaims|OperationalBuildPolicy\\.ReleaseCannotWeakenGates|OperationalRetainedShaders\\.OptionalCapabilitiesPreserveFloatControls|OperationalBuildGate\\.ReleaseInstallRequiresExactPassedEstate|OperationalAttestation\\.FalseExternalEvidenceIsRejected|OperationalAttestation\\.PreflightAndNativeRuntimeRejectFalseHostDevice|OperationalAlignment\\.IncompleteOrAmbiguousAttestationSetIsRejected|OperationalAlignment\\.InstalledReceiptTamperingBlocksReadiness|AlignmentAuthority\\.CompiledReceiptMatchesTheStagedRuntimeAuthority|BuildGateAuthority\\.ReleaseReceiptBindsEveryInstalledProductAtInitialisation)$'])]:
    actions[name]={'revision':head,'kind':'diagnostic','expected_exit':0,'argv':argv,'limits':limits,'frozen_inputs':frozen,'scope':scope}
(folder/'actions.json').write_text(json.dumps(actions,indent=2)+'\n')
print(folder)

"""Close reviewed destination/provenance gaps before the first stage preparation."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,shutil
WORK=Path(__file__).resolve().parent
p=WORK/'prepare_native_stage.py';before=p.read_text()
saved=WORK/'versions/prepare_native_stage.v1.py';assert not saved.exists();saved.parent.mkdir(exist_ok=True);shutil.copyfile(p,saved)
def one(text,old,new):
    assert text.count(old)==1,(old,text.count(old));return text.replace(old,new)
text=one(before,"STAGE = ROOT / SOURCE['stage_root']\nassert not STAGE.exists()",'''stage_relative = Path(SOURCE['stage_root'])
expected_stage = Path('bin/windows-msvc') / ('native-' + HEAD[:7] + '-stage-seven-word-witness')
assert stage_relative == expected_stage and not stage_relative.is_absolute()
assert '..' not in stage_relative.parts and ROOT == ROOT.resolve()
STAGE = ROOT / stage_relative
for ancestor in (ROOT / 'bin', ROOT / 'bin/windows-msvc', STAGE):
    assert not ancestor.is_symlink(), 'Native stage ancestry must not redirect writes'
assert not STAGE.exists()''')
old="""assert json.loads((WORK / 'host-source-review.json').read_text())['accepted']
assert json.loads((WORK / 'ci.json').read_text())['accepted']
assert json.loads((WORK / 'artifact-extraction.json').read_text())['accepted']
"""
new="""host = json.loads((WORK / 'host-source-review.json').read_text())
assert host['accepted'] and host['source_revision'] == SOURCE['baseline_revision']
for name in ('word_value_observer.h', 'retained_compute_test.cpp'):
    assert seal(ROOT / 'tests/backend' / name) == host['outputs'][name]
for name, digest in host['inputs'].items():
    assert seal(WORK / name)['sha256'] == digest
compiled = json.loads((WORK / 'compiled-review.json').read_text())
assert compiled['accepted'] and compiled['source_revision'] == SOURCE['baseline_revision']
assert compiled['schema_seal'] == seal(WORK / 'value-schema.json')
run = json.loads((WORK / 'ci-run.json').read_text())
ci = json.loads((WORK / 'ci.json').read_text())
assert run['source_revision'] == ci['source_revision'] == HEAD
assert ci['accepted'] and ci['run_id'] == run['run_id']
assert ci['metadata']['id'] == run['run_id'] and ci['metadata']['head_sha'] == HEAD
assert ci['metadata']['head_branch'] == SOURCE['branch']
assert ci['metadata']['event'] == 'workflow_dispatch'
assert ci['metadata']['status'] == 'completed' and ci['metadata']['conclusion'] == 'success'
success = {'integration-no-render', 'integration-windows-no-render', 'integration-macos-no-render'}
skipped = {'linux-gate', 'linux-sanitizers', 'windows-build', 'macos-build'}
jobs = ci['jobs']['jobs']
assert len(jobs) == 7 and {j['name'] for j in jobs} == success | skipped
assert all(j['status'] == 'completed' and j['conclusion'] == ('success' if j['name'] in success else 'skipped') for j in jobs)
extraction = json.loads((WORK / 'artifact-extraction.json').read_text())
assert extraction['accepted'] and extraction['source_revision'] == HEAD and extraction['run_id'] == run['run_id']
artifact = extraction['artifact']
assert artifact['name'] == 'sirius-windows-native-attestation'
assert artifact['workflow_run']['id'] == run['run_id'] and artifact['workflow_run']['head_sha'] == HEAD
assert not artifact['expired']
assert seal(WORK / 'artifact.zip') == {'bytes': artifact['size_in_bytes'], 'sha256': artifact['digest'].removeprefix('sha256:')}
api_artifacts = json.loads((WORK / 'ci-artifacts.json').read_text())['artifacts']
assert [a for a in api_artifacts if a['id'] == artifact['id']] == [artifact]
assert len(extraction['extracted_files']) == 40
for record in extraction['extracted_files']:
    assert seal(safe(ARTIFACT, record['path'])) == {k: record[k] for k in ('bytes', 'sha256')}
"""
text=one(text,old,new)
# All new destination checks occur before stage mkdir/copies. The original
# payload verifier remains responsible for its source/build resource contract.
assert text.index('expected_stage =') < text.index('STAGE.mkdir(')
p.write_text(text)
seal=lambda q:{'bytes':q.stat().st_size,'sha256':hashlib.sha256(q.read_bytes()).hexdigest()}
(WORK/'native-preparer-correction.json').write_text(json.dumps({'changes':['Validate unique confined stage and reject symlinks before writes','Join host/source/schema/CI7jobs/APIZIP40payloads to actual revision and bytes'],'before':seal(saved),'after':seal(p),'corrector':seal(Path(__file__)),'scope':'Source correction before first stage preparation; no native execution.'},indent=2)+'\n')
print(json.dumps(seal(p)))

"""Preserve the first structural stop without compiling unchanged stages."""
from pathlib import Path
import hashlib,json,subprocess,zipfile

R=Path.cwd();W=Path(__file__).resolve().parent;V=R/'out/native-cpu-strict-ir-linking'
B=R/'bin/linux-gcc/native-cpu-canonical-metadata'
A=R/'attestations/qualification/2f6182f7/native-cpu-canonical-metadata-first-gate'


def pin(p):
 with p.open('rb') as s:return dict(bytes=p.stat().st_size,sha256=hashlib.file_digest(s,'sha256').hexdigest())


assert not subprocess.check_output(['git','status','--porcelain']).strip()
assert not A.exists()
gate=json.loads((W/'structural-gate.json').read_text());assert not gate['required_gate_met']
owners={}
for name in ['derive','emit']:
 p=W/name/'owner.json';o=json.loads(p.read_text())
 assert o['passed'] and o['status']=='terminal' and o['source_unchanged'] and o['whole_inputs_unchanged']
 assert not o['cleanup_errors'] and not o['remaining_owned_processes']
 owners[name]=dict(identity=pin(p),source=o['source'],seconds=o['wall_seconds'])
paths=[p for p in W.rglob('*') if p.is_file() and '__pycache__' not in p.parts]
paths+=[V/name for name in ['canonical_metadata.py','canonical_metadata_contract.json','canonical_publication.py']]
paths+=list((B/'ispc-source').glob('*.ispc'))+[B/'dopri_phase.ll']
entries={str(p.relative_to(R)):pin(p) for p in sorted(set(paths))}
A.mkdir(parents=True);target=A/'source-and-evidence.zip'
with zipfile.ZipFile(target,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for name in entries:z.write(R/name,name)
with zipfile.ZipFile(target) as z:
 assert set(z.namelist())==set(entries)
 for name,expected in entries.items():
  with z.open(name) as s:assert hashlib.file_digest(s,'sha256').hexdigest()==expected['sha256']
  assert z.getinfo(name).file_size==expected['bytes'] and pin(R/name)==expected
assert target.stat().st_size<100*1024*1024
manifest=dict(status='FIRST_DP_STRUCTURAL_GATE_NOT_MET_STOP',custody_revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
 owners=owners,structural_gate=gate,archive=dict(path=str(target.relative_to(R)),**pin(target),members=entries),
 publication_source=json.loads((W/'publication-source.json').read_text()),
 scope='First unchanged-flags ISPC emission of the new canonical DP source. Literal roots succeeded, whole-context layer escape remained; T/E emission and all runtime/acceptance/frame work stopped. Later inline-publication source removes520 escaping record constructions in source only and remains uncompiled. No strict/native/provider/numerical/speed/scientific/hardware/release qualification.')
with (A/'MANIFEST.json').open('x') as s:json.dump(manifest,s,indent=2);s.write('\n')
print(json.dumps(dict(sealed=str(A.relative_to(R)),archive_bytes=target.stat().st_size,required_structural_gate_met=False,runtime_run=False)))

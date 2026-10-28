"""Preserve the strict effect refinement and rejected compact-interface gate."""
from pathlib import Path
import hashlib,json,subprocess,zipfile

R=Path.cwd();V=R/'out/native-cpu-strict-ir-linking'
A=R/'attestations/qualification/b2cbd4db/native-cpu-private-context-interface-gates'


def pin(p):
 with p.open('rb') as s:return dict(bytes=p.stat().st_size,sha256=hashlib.file_digest(s,'sha256').hexdigest())


assert not subprocess.check_output(['git','status','--porcelain']).strip()
assert not A.exists()
folders=[R/'out/native-cpu-local-effects',R/'out/native-cpu-layer-projection'];owners={};sources=[Path(__file__).resolve()];products=[];gates={}
for folder in folders:
 gates[folder.name]=json.loads((folder/'structural-gate.json').read_text())
 sources+=[p for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts]
 for p in folder.glob('*/owner.json'):
  o=json.loads(p.read_text());assert o['passed'] and o['status']=='terminal' and o['source_unchanged'] and o['whole_inputs_unchanged']
  assert not o['remaining_owned_processes'] and not o['cleanup_errors']
  owners[str(p.relative_to(R))]=dict(identity=pin(p),source=o['source'],seconds=o['wall_seconds'])
 for p in folder.glob('*-result.json'):
  d=json.loads(p.read_text())
  if 'output' in d:products.append(R/d['output']['path'])
assert len(owners)==10
sources+=[V/name for name in ['local_value_effects.cpp','private_finite_closure.json','layer_projection.cpp','layer_projection_contract.json']]
sources+=[p for name in ['native-cpu-local-effects','native-cpu-layer-projection'] for p in (R/'bin/linux-gcc'/name).glob('*.ll')]
A.mkdir(parents=True);bundles={}
for label,paths in [('source-and-evidence',sources),('products',products)]:
 entries={str(p.relative_to(R)):pin(p) for p in sorted(set(paths))};target=A/(label+'.zip')
 with zipfile.ZipFile(target,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for name in entries:z.write(R/name,name)
 with zipfile.ZipFile(target) as z:
  assert set(z.namelist())==set(entries)
  for name,expected in entries.items():
   with z.open(name) as s:assert hashlib.file_digest(s,'sha256').hexdigest()==expected['sha256']
   assert z.getinfo(name).file_size==expected['bytes'] and pin(R/name)==expected
 assert target.stat().st_size<100*1024*1024
 bundles[label]=dict(path=str(target.relative_to(R)),**pin(target),members=entries)
manifest=dict(status='STRICT_SOURCE_AND_LOCAL_SSA_PASS_PERFORMANCE_INTERFACE_GATES_FAIL_STOP',
 custody_revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),owners=owners,gates=gates,bundles=bundles,
 scope='Accurate private ordinary memory refinement and one initialized per-layer compact projection,10 compiler actions passed. Original strict arithmetic/FP state/masks retained. Refinement memory traffic unchanged; projection locals fullypromoted but native instruction count increased225609to428868. No T/E emission, provider adoption, numerical trial, benchmark/frame or scientific/hardware/release qualification; next source removes generic DP metadata/root escapes and compactcopy boundaries.')
with (A/'MANIFEST.json').open('x') as s:json.dump(manifest,s,indent=2);s.write('\n')
print(json.dumps(dict(sealed=str(A.relative_to(R)),owners=len(owners),bundles={k:v['bytes'] for k,v in bundles.items()},runtime_run=False)))

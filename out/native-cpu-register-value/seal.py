"""Preserve the bounded register ABI compiler gate and its exact products."""
from pathlib import Path
from collections import Counter
import hashlib,json,subprocess,zipfile

R=Path.cwd();W=Path(__file__).resolve().parent;V=R/'out/native-cpu-strict-ir-linking'
A=R/'attestations/qualification/359307b1/native-cpu-register-value-compiler-gate'


def pin(p):
 with p.open('rb') as s:return dict(bytes=p.stat().st_size,sha256=hashlib.file_digest(s,'sha256').hexdigest())


assert not subprocess.check_output(['git','status','--porcelain']).strip()
assert not A.exists()
gate=json.loads((W/'static-gate.json').read_text())
assert gate['compiler_static_gate_passed'] and gate['exact_IR_core_calls']==8241 and gate['native_core_relocations']==8240
assert gate['source_review']==pin(V/'private_value_fp_review.json')
products=[R/p['path'] for p in gate['products'].values()]
products+=[R/'bin/linux-gcc/native-cpu-register-value/register_value_abi.so']
sources=[p for p in W.rglob('*') if p.is_file() and '__pycache__' not in p.parts]
sources+=[V/name for name in ['register_value_abi.cpp','private_value_fp_review.json','canonical_metadata.py','canonical_metadata_contract.json']]
sources+=[R/'bin/linux-gcc/native-cpu-register-value'/name for name in ['math.ll','dp.ll','transport.ll','endpoint.ll']]
A.mkdir(parents=True);bundles={}
for label,paths in [('source-and-evidence',sources),('products',products)]:
 entries={str(p.relative_to(R)):pin(p) for p in sorted(set(paths))};target=A/(label+'.zip')
 with zipfile.ZipFile(target,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for name in entries:z.write(R/name,name)
 with zipfile.ZipFile(target) as z:
  assert set(z.namelist())==set(entries)
  for name,p in entries.items():
   with z.open(name) as s:assert hashlib.file_digest(s,'sha256').hexdigest()==p['sha256']
   assert z.getinfo(name).file_size==p['bytes'] and pin(R/name)==p
 assert target.stat().st_size<100*1024*1024
 bundles[label]=dict(path=str(target.relative_to(R)),**pin(target),members=entries)
manifest=dict(status='COMPILER_STATIC_INTERFACE_GATE_ONLY',custody_revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
 execution_revisions=['359307b1bae8b6ee2512621cfcad67caebe687fb','ee0b010b25232ea7877b8b3ec10c66804658b0b2'],
 gate=gate,bundles=bundles,owner_results=dict(Counter('passed' if p['passed'] else 'failed' for p in gate['owners'].values())),
 failed_original_frame_capsule=dict(path='attestations/qualification/844dbeda/native-cpu-shared-value-original-frame-attempt/MANIFEST.json',**pin(R/'attestations/qualification/844dbeda/native-cpu-shared-value-original-frame-attempt/MANIFEST.json')),
 prospective_metadata_source_only=dict(generator=pin(V/'canonical_metadata.py'),contract=pin(V/'canonical_metadata_contract.json'),generated_or_compiled=False),
 scope='Exact private CC/name change only. All4 native objects produced once; failed Transport receipt preserved and object reused after exact common-tail review. No provider/component/runtime/finite/scientific/hardware-domain/release acceptance transferred; no new performance study or frame.')
with (A/'MANIFEST.json').open('x') as s:json.dump(manifest,s,indent=2);s.write('\n')
print(json.dumps(dict(sealed=str(A.relative_to(R)),owners=manifest['owner_results'],bundles={k:v['bytes'] for k,v in bundles.items()})))

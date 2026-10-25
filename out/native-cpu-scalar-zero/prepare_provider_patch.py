"""Return to S scheduling/row SIMD with three genuinely changed scalar objects."""
from pathlib import Path
import difflib
import hashlib
import json
import subprocess

R=Path.cwd();W=Path(__file__).resolve().parent;P=W/'provider'
S=R/'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
T=R/'out/native-cpu-field-simd/provider'


def doc(path):
    return json.loads(path.read_text())


def pin(data):
    return {'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}


def canonical(value):
    return pin(json.dumps(value,sort_keys=True,separators=(',',':')).encode())


assert not (P/'provider.patch').exists()
current=doc(T/'source-contract.json')
before={str(path.relative_to(S)):pin(path.read_bytes()) for path in S.rglob('*') if path.is_file()}
assert canonical(before)==current['full_candidate_source_identity']
old=doc(R/'out/native-cpu-zero-term/provider/source-contract.json')
compiled=doc(W/'cpp-compile-result.json');owner=doc(W/'compile/cpp/owner.json')
assert compiled['all_objects_compiled'] and compiled['new_compilations']==3
assert owner['passed'] and owner['whole_inputs_unchanged'] and not owner['remaining_owned_processes']
assert doc(W/'object-code-inspection.json')['passed']
assert doc(R/'out/native-cpu-scalar-zero-compatibility/result.json')['passed']
objects=[row['object'] for row in compiled['objects']]+old['native_objects'][3:]
assert len(objects)==8
for row in objects:
    assert pin((R/row['path']).read_bytes())=={k:row[k] for k in ('bytes','sha256')}
# Five ordinary source files in a disposable task-local directory suffice to
# reverse the rejected patch. No checkout, toolchain or full Mesa copy is made.
temporary=W/'basis-source';assert not temporary.exists();temporary.mkdir()
names=sorted(current['modified_files'])
assert len(names)==5
for name in names:
    path=temporary/name;path.parent.mkdir(parents=True,exist_ok=True)
    path.write_bytes((S/name).read_bytes())
subprocess.run(['/usr/bin/patch','--batch','--fuzz=0','-R','-p1','-i',str(T/'provider.patch')],cwd=temporary,check=True)
texts={name:(temporary/name).read_text() for name in names}
for name,text in texts.items():
    assert pin(text.encode())==current['modified_files'][name]['basis']
meson='src/gallium/targets/lavapipe/meson.build'
for previous,new in zip(old['native_objects'][:3],objects[:3]):
    a,b=str(R/previous['path']),str(R/new['path'])
    assert texts[meson].count(a)==1;texts[meson]=texts[meson].replace(a,b,1)
for name,a,b in (
    ('src/gallium/frontends/lavapipe/lvp_device.c','SIR-CPU-25.2.8S','SIR-CPU-25.2.8U'),
    ('src/gallium/drivers/llvmpipe/lp_screen.c','Sirius-native-cpu-zero-term-v17','Sirius-native-cpu-scalar-zero-v19')):
    assert texts[name].count(a)==1;texts[name]=texts[name].replace(a,b,1)
patch='';after=dict(before);modified={}
for name,text in texts.items():
    original=(S/name).read_text()
    patch+=''.join(difflib.unified_diff(original.splitlines(True),text.splitlines(True),n=3,
        fromfile='a/'+name,tofile='b/'+name))
    after[name]=pin(text.encode());modified[name]={'basis':before[name],'candidate':after[name]}
inputs=dict(compiled['source_inputs'])
for row in objects:
    inputs[row['path']]={k:row[k] for k in ('bytes','sha256')}
for path in (W/'cpp-compile-result.json',W/'compile/cpp/owner.json',W/'object-code-inspection.json',
             W/'derivation.json',W/'scalar-zero.patch',W/'prepare_provider_patch.py',
             R/'out/native-cpu-scalar-zero-compatibility/result.json',
             R/'out/native-cpu-scalar-zero-compatibility/run/owner.json',
             R/'out/native-cpu-zero-term/provider/source-contract.json',T/'source-contract.json',T/'provider.patch'):
    inputs[str(path.relative_to(R))]=pin(path.read_bytes())
contract={'basis_source_identity':current['full_candidate_source_identity'],
    'full_candidate_source_identity':canonical(after),'modified_files':modified,'native_objects':objects,
    'patch':pin(patch.encode()),'source_inputs':inputs,'new_CPP_compilations':3,'new_ISPC_compilations':0,
    'reused_native_objects':5,'guard_shared_byte_authority':'out/native-cpu-active-prefix/native_packet.hpp',
    'scope':'U changes only the actual D scalar Term inline/exact-zero arithmetic from S; all9380instructions, S raw3/pureF/Iadapter/Q scheduling and every count/FP/guard/cache/governor contract retained. Reverses rejected field-SIMD scheduling. Actual fresh scalarTransport38fixtures/all187960words PASS. No runtime benefit/frame/science/release/adoption claim.'}
P.mkdir(exist_ok=True)
(P/'provider.patch').write_text(patch);(P/'source-contract.json').write_text(json.dumps(contract,indent=2)+'\n')
# This partial source copy is only a derivation aid; the resulting patch and
# original/candidate hashes preserve the useful result.
for name in names:(temporary/name).unlink()
for path in sorted((path for path in temporary.rglob('*') if path.is_dir()),key=lambda path:len(path.parts),reverse=True):path.rmdir()
temporary.rmdir()
assert before=={str(path.relative_to(S)):pin(path.read_bytes()) for path in S.rglob('*') if path.is_file()}
print(json.dumps({'new_scalar_objects':3,'reused_native_objects':5,'modified_files':5,'shared_Mesa_modified':False,'temporary_source_copy_removed':True}))

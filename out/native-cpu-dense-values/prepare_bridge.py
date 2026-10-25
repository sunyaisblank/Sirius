"""Compile the packing adapter only after preserved AC frame rejection."""
from pathlib import Path
import hashlib
import json
import shlex
import subprocess

R = Path.cwd()
W = Path(__file__).resolve().parent
A = R / 'out/native-cpu-dense-arithmetic'


def doc(path):
    return json.loads(path.read_text())


def pin(path):
    return {'bytes': path.stat().st_size, 'sha256': hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()}


assert not subprocess.check_output(['git','status','--porcelain']).strip()
failure = R / 'out/native-cpu-dense-arithmetic-cohort/cold-frame'
owner = doc(failure / 'frame/owner.json')
assert owner['status']=='terminal' and not owner['passed'] and owner['stop_reason']=='diagnostic_timeout'
assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
assert not doc(failure / 'result.json')['accepted']
manifest = R / 'attestations/qualification' / owner['source']['revision'][:8] / 'native-cpu-dense-arithmetic-frame-timeout/MANIFEST.json'
assert not doc(manifest)['accepted']
assert subprocess.check_output(['git','show','@{upstream}:'+str(manifest.relative_to(R))]) == manifest.read_bytes()
assert doc(A / 'cpp-code-inspection.json')['passed']
assert doc(A / 'cpp/run/owner.json')['passed']
assert doc(A / 'cpp/run/stdout.log')['compared_words']==1425600
old = doc(A / 'prepared-cpp.json')
driver, = [row for row in old['compilations'] if row['name']=='driver']
flags = driver['argv'][:driver['argv'].index('-c')]
assert '-std=c++2c' in flags and '-march=x86-64' in flags and '-ffp-contract=off' in flags
inputs = dict(old['inputs'])
for path in (A / 'prepared-cpp.json', A / 'cpp-code-inspection.json', A / 'cpp/run/owner.json',
             A / 'cpp/run/stdout.log', failure / 'frame/owner.json', failure / 'result.json',
             manifest, *W.glob('*.cpp'), *W.glob('*.py'), * (W / 'bridge').glob('*.py')):
    inputs[str(path.relative_to(R))] = pin(path)
source = W / 'values_bridge.cpp'
dependencies = subprocess.check_output([*flags,'-M','-MT','values-bridge',str(source)],text=True)
for name in shlex.split(dependencies.replace('\\\n',' ').split(':',1)[1]):
    path = Path(name).resolve(strict=True)
    inputs[str(path)] = pin(path)
for name, expected in inputs.items():
    assert pin(R / name)==expected,name
output = 'bin/linux-gcc/native-cpu-dense-values/objects/values_bridge.o'
assert not (R / output).exists() and not (W / 'prepared-bridge.json').exists()
revision = subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
scope = 'One new baseline-x86-64 C++26 packing adapter compilation only, original strict flags. No arithmetic recompile: exact checked AC math/adapter object is reused later with opcode joins. New wrappers snapshot both operands and old output, call the unchanged dense scratch ABI, and copy the masked destination back. No new numerical study, provider execution or runtime/science/release claim.'
prepared = {'revision':revision,'inputs':inputs,'object':output,
            'argv':[*flags,'-c',str(source),'-o',str(R / output)],'scope':scope}
(W / 'prepared-bridge.json').write_text(json.dumps(prepared,indent=2)+'\n')
action = {'revision':revision,'kind':'build_only','expected_exit':0,'software_provider':False,
          'limits':{'seconds':300,'rss_kib':4194304,'reserve_kib':2097152},
          'argv':['/usr/bin/python3','-B',str(W / 'build_bridge.py')],
          'frozen_inputs':sorted(set(inputs)|{str((W / 'prepared-bridge.json').relative_to(R))}),
          'scope':scope}
assert not (W / 'bridge/actions.json').exists()
(W / 'bridge/actions.json').write_text(json.dumps({'build':action},indent=2)+'\n')
print(json.dumps({'prepared_new_CPP_compilations':1,'arithmetic_compilations':0,'executed':False}))

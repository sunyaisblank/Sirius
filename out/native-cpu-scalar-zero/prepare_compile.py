"""Compile three changed scalar objects with the exact original D strict flags."""
from pathlib import Path
import hashlib
import json
import subprocess

R=Path.cwd()
W=Path(__file__).resolve().parent
B=R/'bin/linux-gcc/native-cpu-scalar-zero'


def doc(path):
    return json.loads(path.read_text())


def pin(path):
    with path.open('rb') as stream:
        return {'bytes':path.stat().st_size,'sha256':hashlib.file_digest(stream,'sha256').hexdigest()}


assert not subprocess.check_output(['git','status','--porcelain']).strip()
assert not (W/'prepared-compile.json').exists() and not (W/'compile/actions.json').exists()
assert not (B/'objects').exists()
subprocess.run(['/usr/bin/python3','-B',str(W/'derive.py')],check=True)
derived=doc(W/'derivation.json'); old=doc(R/'out/native-cpu-static-program/batch-compile-result.json')
inputs={**old['source_inputs'],**derived['source_inputs']}
objects=[]
for row in old['objects']:
    assert row['exit']==0
    slug=Path(row['argv'][-3]).stem;stage=derived['stages'][slug]
    for key in ('header','source'):
        value=stage[key];inputs[value['path']]={k:value[k] for k in ('bytes','sha256')}
    argv=list(row['argv'])
    include=argv.index(str(R/'bin/linux-gcc/native-cpu-static-program/typed-batch'))
    assert argv[include-1]=='-isystem'
    argv[include]=str(B/'typed-batch')
    argv[-3]=str(R/stage['source']['path'])
    output='bin/linux-gcc/native-cpu-scalar-zero/objects/'+slug+'.o';argv[-1]=str(R/output)
    inverse=list(argv)
    inverse[include],inverse[-3],inverse[-1]=row['argv'][include],row['argv'][-3],row['argv'][-1]
    assert inverse==row['argv']
    objects.append({'kind':'cpp','stage':slug,'output':output,'argv':argv,
                    'original_D_argv':row['argv'],'argv_inverse_exact':True})
tool_receipt=R/'out/native-cpu-host-representation/probe/preparation.json'
for name,expected in doc(tool_receipt)['tool_dependencies'].items():
    assert pin(Path(name))==expected,name
    inputs[name]=expected
for path in (W/'prepare_compile.py',W/'build.py',W/'derivation.json',W/'scalar-zero.patch',
             W/'compile/run_owned.py',W/'compile/owned_guard.py',tool_receipt,
             R/'out/native-cpu-static-program/batch-compile-result.json',
             R/'out/native-cpu-row-ownership/cold-frame/run_owned.py',
             R/'out/native-cpu-row-ownership/cold-frame/owned_guard.py'):
    inputs[str(path.relative_to(R))]=pin(path)
for name,expected in inputs.items():
    assert pin(R/name)==expected,name
revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
scope='Three new full scalar D objects, only Term inline/zero arm changed; exact original D compiler/flags and all9380instructions retained. No ISPC/host/adapter/provider recompilation or runtime qualification in this action.'
(W/'prepared-compile.json').write_text(json.dumps({'revision':revision,'inputs':inputs,'objects':objects,'scope':scope},indent=2)+'\n')
action={'revision':revision,'kind':'compile','expected_exit':0,'software_provider':False,
    'limits':{'seconds':60,'rss_kib':4194304,'reserve_kib':2097152},
    'argv':['/usr/bin/python3','-B',str(W/'build.py'),'cpp'],
    'frozen_inputs':sorted(set(inputs)|{str((W/'prepared-compile.json').relative_to(R))}),'scope':scope}
(W/'compile/actions.json').write_text(json.dumps({'cpp':action},indent=2)+'\n')
print(json.dumps({'new_scalar_compilations_planned':3,'executed':False}))

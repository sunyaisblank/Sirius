"""Seal one bounded actual-frame route observation, not a frame retry/gate."""
from pathlib import Path
import hashlib
import json
import re
import subprocess

R=Path.cwd()
W=Path(__file__).resolve().parent
C=R/'out/native-cpu-dense-values-cohort'


def doc(path):return json.loads(path.read_text())
def pin(path):return {'bytes':path.stat().st_size,'sha256':hashlib.file_digest(path.open('rb'),'sha256').hexdigest()}


assert not subprocess.check_output(['git','status','--porcelain']).strip()
owner=doc(C/'cold-frame/frame/owner.json')
assert owner['status']=='terminal' and not owner['passed'] and owner['stop_reason']=='diagnostic_timeout'
assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
manifest=R/'attestations/qualification'/owner['source']['revision'][:8]/'native-cpu-dense-values-frame-timeout/MANIFEST.json'
assert not doc(manifest)['accepted']
assert subprocess.check_output(['git','show','@{upstream}:'+str(manifest.relative_to(R))])==manifest.read_bytes()
old=doc(C/'cold-frame/actions.json')['frame']
assert old['limits']=={'seconds':300,'rss_kib':4194304,'reserve_kib':2097152}
assert not (W/'actions.json').exists() and not (W/'observe').exists()
revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
executable=R/'bin/linux-gcc/native-cpu-upload-aware-staging/frame/sirius_render_tests'
temp=R/'renders'/('native-cpu-frame-route-'+revision[:8])/'tmp'
assert not temp.exists();temp.mkdir(parents=True)
argv=['/usr/bin/gdb','--batch','--nx','-ex','set pagination off','-ex','set confirm off',
      '-ex','set print thread-events off','-ex','source '+str(W/'observer.py'),
      '-ex','run','-ex','kill','-ex','quit','--args',str(executable),
      '--gtest_filter=RenderCommandParse.ExplicitGpuRequestRunsVulkanWhenDevicePresent','--gtest_color=no']
prepared={'source':revision,'gdb_argv':argv,'original_frame_action':pin(C/'cold-frame/actions.json'),
          'executable':pin(executable),'scope':'First actual frame route only; debugger deliberately stops the inferior after native3/SIMD3/dense4 coverage. No output/timing/numerical/frame qualification.'}
(W/'prepared.json').write_text(json.dumps(prepared,indent=2)+'\n')
inputs=set(old['frozen_inputs'])
inputs|={str(p.relative_to(R)) for p in W.glob('*') if p.is_file()}
inputs|={str(p.relative_to(R)) for p in (C/'cold-frame/frame/owner.json',C/'cold-frame/result.json',manifest,
          R/'out/native-cpu-dense-values/provider/linked-native-proof.json')}
for tool in (Path('/usr/bin/gdb'),Path('/usr/bin/python3')):
    inputs.add(str(tool.resolve()))
    output=subprocess.check_output(['/usr/bin/ldd','-r',str(tool)],text=True,stderr=subprocess.STDOUT)
    assert 'undefined symbol' not in output and 'not found' not in output
    inputs|={str(Path(p).resolve()) for p in re.findall(r'(/[^\s()]+)',output) if Path(p).is_file()}
assert all((R/name).is_file() for name in inputs)
action={'revision':revision,'kind':'diagnostic','expected_exit':0,'software_provider':True,
        'runtime_provider':old['runtime_provider'],'runtime_icd':old['runtime_icd'],
        'limits':{'seconds':45,'rss_kib':4194304,'reserve_kib':2097152},
        'argv':['/usr/bin/python3','-B',str(W/'run_prefix.py')],
        'frozen_inputs':sorted(inputs),'build_environment':{'TMPDIR':str(temp)},
        'scope':'One bounded diagnostic answering which positive native/row/SIMD/C++ routes the original FP32 frame actually executes after software warm-up. Same genuine P executable and AD provider, original scene and cache policy, unchanged4GiB/2GiB memory guards. Deliberately terminate after first route coverage, at most45s; no qualifying frame retry, timing comparison or numerical/scientific/release claim.'}
(W/'actions.json').write_text(json.dumps({'observe':action},indent=2)+'\n')
print(json.dumps({'prepared_actual_frame_route_only':True,'outer_seconds':45,'executed':False}))

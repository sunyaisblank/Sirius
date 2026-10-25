"""Seal one original-input, deliberately stopped steady-render PC diagnostic."""
from pathlib import Path
import hashlib,json,os,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent
J=lambda p:json.loads(p.read_text())
def pin(p):
 with p.open('rb') as f:return dict(bytes=p.stat().st_size,sha256=hashlib.file_digest(f,'sha256').hexdigest())
def write(p,v):
 with p.open('x') as f:json.dump(v,f,indent=2);f.write('\n')
assert not subprocess.check_output(['git','status','--porcelain']).strip()
revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
C=R/'out/native-cpu-normalization-tags-cohort/cold-frame'
old=J(C/'actions.json')['frame'];owner=J(C/'frame/owner.json')
assert owner['status']=='terminal' and not owner['passed'] and owner['stop_reason']=='diagnostic_timeout'
assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
assert old['limits']==dict(seconds=300,rss_kib=4194304,reserve_kib=2097152)
assert not (W/'actions.json').exists() and not (W/'observe').exists()
manifest=R/'attestations/qualification/629e3794/native-cpu-dispatch-admission-frame-timeout/MANIFEST.json'
recovery=J(manifest)['resources_recovery'];cache=R/recovery['source']
exe=R/'bin/linux-gcc/native-cpu-owned-lower-increments/frame/sirius_render_tests'
dest=exe.parent/'resources';assert not dest.exists()
files={n:{k:v[k] for k in ('bytes','sha256')} for n,v in recovery['files'].items()}
assert len(files)==32
for name,value in files.items():
 p=cache/name;assert p.is_file() and not p.is_symlink() and pin(p)==value
 assert p.stat().st_mtime_ns==recovery['files'][name]['mtime_ns']
 assert p.stat().st_mode&0o777==recovery['files'][name]['mode']
dest.mkdir()
for name,value in files.items():
 target=dest/name;target.parent.mkdir(parents=True,exist_ok=True);os.link(cache/name,target)
 assert os.path.samefile(cache/name,target) and pin(target)==value
write(W/'resource-restore.json',dict(source=revision,cache=str(cache.relative_to(R)),destination=str(dest.relative_to(R)),files=files,payload_bytes_copied=0))
inputs=set(old['frozen_inputs'])
# Removed after the previous pushed recovery: the unused nested decoder had
# a wrong checkout root. Preserve its historical absence explicitly.
removed='out/native-cpu-normalization-tags-cohort/cold-frame/decode_profile.py'
assert removed in inputs and not (R/removed).exists();inputs.remove(removed)
assert all((R/p).is_file() for p in inputs)
before=J(C/'frame/inputs-before.json')
fixed_controller='out/native-cpu-normalization-tags-cohort/cold-frame/profile_existing_frame.py'
for name in inputs:
 if name==fixed_controller:
  assert pin(R/name)==dict(bytes=3993,sha256='d297e2077c9833e3e9c60cc1092bbc70d7d11db89cb737a208dc779edef9c27e')
 else:assert pin(R/name)=={k:before[name][k] for k in ('bytes','sha256')},name
temp=R/'renders'/('native-cpu-steady-profile-'+revision[:8])/'tmp'
assert not temp.exists();temp.mkdir(parents=True)
plan=dict(source=revision,original_argv=old['argv'],executable=str(exe.relative_to(R)),executable_identity=pin(exe),
 original_action=pin(C/'actions.json'),original_terminal_owner=pin(C/'frame/owner.json'),
 trigger_completed_tiles=128,trigger_total_tiles=16384,wait_seconds=30,captures=1,capture_seconds=3,
 sampler='out/native-cpu-runtime-observe/all_threads.py',sampler_identity=pin(R/'out/native-cpu-runtime-observe/all_threads.py'),
 previous_initialization_capture_reused=False,removed_unused_historical_input=removed,
 historical_observer_only_input_change=dict(path=fixed_controller,identity=pin(R/fixed_controller),reason='Previously pushed completed-tile trigger repair; no render input change. This diagnostic uses its separate controller.'),
 scope='One purpose-specific steady-render USER-PC diagnostic. Original128x128 one-SPP scene and executable/provider/resource/count/cache bytes retained; no acceptance retry or completed frame. Stop deliberately immediately after a single3s100Hz capture after at least128 completed tiles. No speed comparison, full science or release claim.')
write(W/'plan.json',plan)
inputs|={str(p.relative_to(R)) for p in W.iterdir() if p.is_file()}
inputs|={str(manifest.relative_to(R)),str((C/'frame/owner.json').relative_to(R))}
action=dict(revision=revision,kind='diagnostic',expected_exit=0,software_provider=True,
 runtime_provider=old['runtime_provider'],runtime_icd=old['runtime_icd'],
 limits=dict(seconds=45,rss_kib=4194304,reserve_kib=2097152),
 argv=['/usr/bin/python3','-B',str(W/'capture_and_stop.py')],frozen_inputs=sorted(inputs),
 build_environment=dict(TMPDIR=str(temp)),scope=plan['scope'])
write(W/'actions.json',dict(observe=action))
print(json.dumps(dict(prepared=True,executed=False,trigger_tiles=128,outer_seconds=45,resource_payload_bytes_copied=0)))

"""Seal the first short steady USER-PC observation with corrected value facts."""
from pathlib import Path
import hashlib,json,os,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;J=lambda p:json.loads(p.read_text())
def pin(p):
 with p.open('rb') as f:return dict(bytes=p.stat().st_size,sha256=hashlib.file_digest(f,'sha256').hexdigest())
def write(p,v):
 with p.open('x') as f:json.dump(v,f,indent=2);f.write('\n')
assert not subprocess.check_output(['git','status','--porcelain']).strip()
revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
C=R/'out/native-cpu-local-identities-cohort/cold-frame';old=J(C/'actions.json')['frame'];owner=J(C/'frame/owner.json')
assert owner['status']=='terminal' and not owner['passed'] and owner['stop_reason']=='diagnostic_timeout'
assert owner['source_unchanged'] and owner['whole_inputs_unchanged'] and not owner['remaining_owned_processes'] and not owner['cleanup_errors']
assert old['limits']==dict(seconds=300,rss_kib=4194304,reserve_kib=2097152)
assert not (W/'actions.json').exists() and not (W/'observe').exists()
adoption=R/'attestations/qualification/df1d3baa/native-cpu-local-identities-frame-timeout/MANIFEST.json'
assert J(adoption)['all22_original_gates_pass'] and J(adoption)['all8_original_controls_pass']
assert subprocess.check_output(['git','show','@{upstream}:'+str(adoption.relative_to(R))])==adoption.read_bytes()
manifest=R/'attestations/qualification/629e3794/native-cpu-dispatch-admission-frame-timeout/MANIFEST.json';recovery=J(manifest)['resources_recovery'];cache=R/recovery['source']
exe=R/'bin/linux-gcc/native-cpu-owned-lower-increments/frame/sirius_render_tests';dest=exe.parent/'resources';assert not dest.exists()
files={n:{k:v[k] for k in ('bytes','sha256')} for n,v in recovery['files'].items()};assert len(files)==32
for name,value in files.items():
 p=cache/name;assert p.is_file() and not p.is_symlink() and pin(p)==value
 assert p.stat().st_mtime_ns==recovery['files'][name]['mtime_ns'] and p.stat().st_mode&0o777==recovery['files'][name]['mode']
dest.mkdir()
for name,value in files.items():
 target=dest/name;target.parent.mkdir(parents=True,exist_ok=True);os.link(cache/name,target)
 assert os.path.samefile(cache/name,target) and pin(target)==value
write(W/'resource-restore.json',dict(source=revision,cache=str(cache.relative_to(R)),destination=str(dest.relative_to(R)),files=files,payload_bytes_copied=0))
inputs=set(old['frozen_inputs']);assert all((R/p).is_file() for p in inputs)
before=J(C/'frame/inputs-before.json')
for name in inputs:assert pin(R/name)=={k:before[name][k] for k in ('bytes','sha256')},name
proof=R/'out/native-cpu-local-identities/provider/linked-native-proof.json';assert J(proof)['passed']
temp=R/'renders'/('native-cpu-coherent-profile-'+revision[:8])/'tmp';assert not temp.exists();temp.mkdir(parents=True)
plan=dict(source=revision,original_argv=old['argv'],executable=str(exe.relative_to(R)),executable_identity=pin(exe),original_action=pin(C/'actions.json'),original_terminal_owner=pin(C/'frame/owner.json'),trigger_completed_tiles=128,trigger_total_tiles=16384,wait_seconds=30,captures=1,capture_seconds=3,sampler='out/native-cpu-runtime-observe/all_threads.py',sampler_identity=pin(R/'out/native-cpu-runtime-observe/all_threads.py'),current_native_proof=pin(proof),previous_profiles_reused=False,
 scope='One first current AN corrected-value-fact USER-PC diagnostic. Older AJ samples predate the twelve reachable writer fixes and local identity path, so they are not current cost attribution. Keep original128x128/one-SPP frame executable/provider/resources/settings/precision/counts/cache/governors. After128 completed tiles take one3s100Hz USER-PC capture, then deliberately stop immediately. Fixed30s trigger wait/45s outer/4GiB RSS/2GiB reserve; no retry if incomplete. This is not another300s acceptance frame, a speed comparison, completed image or full scientific/release qualification.')
write(W/'plan.json',plan)
inputs|={str(p.relative_to(R)) for p in W.iterdir() if p.is_file()}
inputs|={str(p.relative_to(R)) for p in (manifest,adoption,C/'frame/owner.json',C/'result.json',proof,R/'out/native-cpu-runtime-observe/all_threads.py',R/'out/native-cpu-runtime-observe/decode_all_threads.py')}
action=dict(revision=revision,kind='diagnostic',expected_exit=0,software_provider=True,runtime_provider=old['runtime_provider'],runtime_icd=old['runtime_icd'],limits=dict(seconds=45,rss_kib=4194304,reserve_kib=2097152),argv=['/usr/bin/python3','-B',str(W/'capture_and_stop.py')],frozen_inputs=sorted(inputs),build_environment=dict(TMPDIR=str(temp)),scope=plan['scope'])
write(W/'actions.json',dict(observe=action));print(json.dumps(dict(prepared=True,executed=False,trigger_tiles=128,outer_seconds=45,resource_payload_bytes_copied=0)))

"""One owned compile/archive/two-link action; never launch a Sirius workload."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess
import zipfile

root = Path.cwd()
work = Path(__file__).resolve().parent
prepared = json.loads((work / 'prepared-build.json').read_text())


def pin(path):
    with path.open('rb') as stream:
        return {'bytes': path.stat().st_size, 'sha256': hashlib.file_digest(stream, 'sha256').hexdigest()}


def write(value):
    (work / 'host-build-result.json').write_text(json.dumps(value, indent=2) + '\n')


def unchanged():
    for name, value in prepared['source_inputs'].items():
        assert pin(root / name) == value, name


assert not (work / 'host-build-result.json').exists()
unchanged()
authority = prepared['authority']
record = {'passed': False, 'commands': [], 'source_inputs': prepared['source_inputs'],
    'new_cpp_compilations': 0, 'observer_compilations': 0, 'new_links': 0,
    'new_shader_ISPC_provider_compilations': 0,
    'same_actual_executable_for_both_provider_arms': True,
    'source_proof': prepared['source_proof'], 'all3modules_once_in_both_executables': False,
    'scope': 'One private host compile and original audited-host/frame links. Fulloriginalshaderbytes/counts/limits/precision and Mprovider preserved. No runtime/numerical/performance/frame/scientific/release claim.'}
write(record)
observer = root / prepared['observer']
recovery = authority['observer_recovery']
archive = root / recovery['archive']
assert pin(archive) == recovery['archive_identity']
assert not observer.exists()
with zipfile.ZipFile(archive) as handle:
    data = handle.read(recovery['member'])
assert {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()} == recovery['object_identity']
with observer.open('xb') as stream:
    stream.write(data)
os.chmod(observer, recovery['mode'])
assert pin(observer) == recovery['object_identity']
record['restored_observer'] = {'path': str(observer.relative_to(root)), **pin(observer),
    'archive': recovery['archive'], 'member': recovery['member'], 'members_restored': 1}
del data
write(record)

def run(index):
    command = prepared['commands'][index]
    result = subprocess.run(command['argv'], cwd=command['cwd'])
    record['commands'].append({**command, 'exit': result.returncode})
    if index == 0:
        record['new_cpp_compilations'] = 1
    if index in (2, 3):
        record['new_links'] += 1
    write(record)
    if result.returncode:
        raise SystemExit(result.returncode)


run(0)
object_path = root / prepared['object']
assert object_path.is_file()
record['objects'] = [{'argv': prepared['commands'][0]['argv'], 'exit': 0,
                      'object': {'path': prepared['object'], **pin(object_path)}}]
old_archive = root / authority['backend_archive']
archive_path = root / prepared['backend_archive']
assert pin(old_archive) == authority['backend_archive_identity']
assert not archive_path.exists()
shutil.copyfile(old_archive, archive_path)
run(1)
members = subprocess.check_output(['/usr/bin/ar', 't', str(old_archive)], text=True).splitlines()
assert members == subprocess.check_output(['/usr/bin/ar', 't', str(archive_path)], text=True).splitlines()
assert members.count('retained_compute.cpp.o') == 1
different, retained = [], {}
for member in members:
    old = subprocess.check_output(['/usr/bin/ar', 'p', str(old_archive), member])
    new = subprocess.check_output(['/usr/bin/ar', 'p', str(archive_path), member])
    if old != new:
        different.append(member)
    else:
        retained[member] = {'bytes': len(new), 'sha256': hashlib.sha256(new).hexdigest()}
assert different == ['retained_compute.cpp.o']
record['changed_archive_members'] = different
record['unchanged_archive_members'] = retained
record['backend_archive'] = prepared['backend_archive']
record['backend_archive_identity'] = pin(archive_path)
write(record)
run(2)
run(3)
host = root / prepared['executable']
frame = root / prepared['frame_executable']
record.update(executable=prepared['executable'], executable_identity=pin(host),
              frame_executable=prepared['frame_executable'], frame_executable_identity=pin(frame))
record['module_joins'] = {}
for stage, item in authority['modules'].items():
    path = root / item['module']
    assert pin(path) == {key: item[key] for key in ('bytes', 'sha256')}
    code = path.read_bytes()
    host_count = host.read_bytes().count(code)
    frame_count = frame.read_bytes().count(code)
    assert host_count == frame_count == 1, (stage, host_count, frame_count)
    record['module_joins'][stage] = {'module': item['module'], **pin(path),
        'complete_occurrences': 1, 'executable_complete_occurrences': host_count,
        'frame_executable_complete_occurrences': frame_count}
record['all3modules_once_in_both_executables'] = True
record['unchanged_other_link_inputs'] = authority['unchanged_other_link_inputs']
record['frame_resources'] = authority['frame_resources']
record['frame_resource_payload_bytes_copied'] = 0
record['source_inputs_unchanged'] = True
record['original_backend_archive_unchanged'] = pin(old_archive) == authority['backend_archive_identity']
unchanged()
record['passed'] = True
write(record)
print(json.dumps({'passed': True, 'new_cpp_compilations': 1, 'observer_compilations': 0,
                  'new_links': 2, 'changed_archive_members': different,
                  'all3modules_once_in_both_executables': True}))

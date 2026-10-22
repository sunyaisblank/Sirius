"""Archive terminal candidate evidence and actual products; never run a workload."""
from pathlib import Path
import hashlib, json, os, re, stat, subprocess, zipfile
R = Path.cwd(); W = Path(__file__).resolve().parent
SOURCE = '08add7e5633bfbb5aa19343e45497cd33a2edaa3'
D = R / 'attestations/qualification/08add7e5/native-original-input-ownership-finite-outcome'
def doc(p): return json.loads(p.read_text())
def pin(p): return {'bytes': p.stat().st_size, 'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()}
def write(p, d):
    assert not p.exists(), p
    p.write_text(json.dumps(d, indent=2) + '\n')
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == SOURCE
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
assert not D.exists()
owners = {}; products = set()
for p in sorted(W.rglob('owner.json')):
    if 'provider-history' in p.parts: continue
    o = doc(p)
    assert o['status'] == 'terminal' and not o['cleanup_errors'] and not o['remaining_owned_processes']
    assert o['source_unchanged'] and o['whole_inputs_unchanged']
    assert o['observed_births'] and all(x['absent'] for x in o['observed_births'])
    for b in o['observed_births']:
        live = Path('/proc') / str(b['pid']) / 'stat'
        if live.exists(): assert live.read_text().rsplit(')', 1)[1].split()[19] != b['saved_start_ticks']
    before = doc(p.parent / 'inputs-before.json')
    assert before == doc(p.parent / 'inputs-after.json')
    for key, value in before.items():
        if key.startswith('bin/'):
            q = R / key
            assert pin(q) == {k: value[k] for k in ['bytes', 'sha256']}, key
            products.add(key)
    owners[str(p.relative_to(R))] = {'identity': pin(p), 'passed': o['passed'], 'source': o['source'],
        'stop_reason': o['stop_reason'], 'child_exit': o['child_exit'], 'wall_seconds': o['wall_seconds'],
        'peak_sampled_rss_kib': o['peak_sampled_rss_kib'], 'observed_births': o['observed_births']}
benefit = doc(W / 'readers/benefit.json'); decimal = doc(W / 'readers/decimal-agreement.json')
assert benefit['audit_passed'] and benefit['checks_total'] == benefit['checks_passed'] == 22
assert decimal['agrees_on_all22'] and decimal['all22_pass']
controls = {}
for i in range(1, 9):
    p = W / f'runtime/control-{i}/owner.json'
    controls[str(i)] = ({k: doc(p)[k] for k in ['passed', 'wall_seconds', 'stop_reason', 'child_exit', 'peak_sampled_rss_kib']}
                       if p.exists() else {'status': 'not_run'})
all8 = (W / 'readers/controls-1-8-acceptance.json').exists()
if all8: assert doc(W / 'readers/controls-1-8-acceptance.json')['all8_original_controls_pass']
frame = {'attempted': False, 'accepted': False, 'reason': 'All8 original controls not admitted.'}
F = W / 'cold-frame/frame'
if (F / 'owner.json').exists():
    assert all8
    o = doc(F / 'owner.json'); stdout = re.sub(r'\x1b\[[0-?]*[ -/]*[@-~]', '', (F / 'stdout.log').read_text())
    stderr = (F / 'stderr.log').read_text()
    progress = re.findall(r'Tile (\d+)/(\d+)', stdout)
    scenes = [json.loads(x) for x in re.findall(r'(?:Source scene|Scene) evidence: (\{.*\})', stdout)]
    assert len(scenes) == 2 and scenes[0] == scenes[1]
    assert scenes[0]['width'] == scenes[0]['height'] == 128 and scenes[0]['samples_per_pixel'] == 1
    snapshots = re.findall(r'^\[Vulkan\] Retained progress: (.*)$', stderr, re.M)
    stages = {}
    for name, row in re.findall(r'^\[Vulkan\] Retained cumulative physical stage: (\w+), (.*)$', stderr, re.M):
        stages[name] = {k: float(v) if '.' in v else int(v) for k, v in re.findall(r'(\w+)=([\d.]+)', row)}
    frame = {'attempted': True, 'accepted': o['passed'], 'stop_reason': o['stop_reason'],
        'wall_seconds': o['wall_seconds'], 'peak_sampled_rss_kib': o['peak_sampled_rss_kib'],
        'last_progress': list(map(int, progress[-1])) if progress else None,
        'XML_present': (F / 'gtest.xml').exists(), 'scene': scenes[0],
        'last_complete_snapshot': snapshots[-1] if snapshots else None,
        'last_complete_snapshot_physical_stages': stages,
        'scope': 'Single original changed cold frame only. Partial counters are not full-frame output, exclusive costs, a historical speed comparison or an extrapolation.'}
    if not o['passed']:
        assert o['stop_reason'] == 'diagnostic_timeout' and o['child_exit'] == -15
        assert not (F / 'gtest.xml').exists()
    write(W / 'cold-frame/result.json', frame)
compiler = doc(W / 'compiler-result.json'); code = doc(W / 'ordinary/ordinary-code-review.json')
findings = {'candidate_source': SOURCE, 'public_source_unchanged': '0109ab462edade524667d20303507aa09adcbaf6',
    'all22_original_benefit_gates_pass': True, 'independent_Decimal_agrees': True,
    'sole_ABBA': {'order': ['a1', 'b1', 'b2', 'a2'], 'internal_repeats_per_slot': 3, 'extra_retries': 0},
    'whole_call_gains_percent': [{k: x[k] for k in ['direction', 'family', 'improvement_percent']}
        for x in benefit['checks'] if x['scope'] == 'complete_call_ms'],
    'minimum_inclusive_stage_gain_percent': min(x['improvement_percent'] for x in benefit['checks'] if x['scope'] != 'complete_call_ms'),
    'provisioning_inclusive_API_gains_percent': [x['improvement_percent'] for x in benefit['accounting_checks']],
    'all32_control8_and_ABBA_readbacks_match_originals': True,
    'original_controls': controls, 'all8_original_controls_admitted': all8,
    'cold_frame': frame, 'actual_new_builds': {'shader_compilations': 2, 'CXX_compilations': 1, 'backend_links': 1,
        'frame_links': int((W / 'cold-frame/build-result.json').exists())},
    'initial_metadata_reader_failure_preserved': 'modules/owner.json; Transport compiler succeeded, signed array-bound reader stopped before Endpoint. Reused exact Transport bytes; no compiler or workload repeated for correction.',
    'actual_ordinary_called_function_bytes': {k: x['total_called_function_bytes'] for k, x in code['rows'].items()},
    'actual_shared_bytes': {k: compiler['modules'][k]['actual_shared_bytes'] for k in ['Transport', 'Endpoint']},
    'source_or_provider_adopted_publicly': False, 'full_goal_complete': False,
    'scope': 'Private finite candidate evidence. Full science, July image identity, original cold-frame acceptance, native-window/Radeon scenes, Mandatory and exact-main installed-release qualification remain distinct obligations.'}
D.mkdir(parents=True)
write(D / 'findings.json', findings)
write(D / 'ownership-closure.json', {'all_owned_runs_terminal_and_reaped': True, 'owners_count': len(owners), 'owners': owners})
diagnostics = sorted(p for p in W.rglob('*') if p.is_file())
old_cleanup = R / 'attestations/qualification/d9dbd657/native-phase-sharing-finite-pass-frame-timeout/cleanup.json'
if old_cleanup.exists(): diagnostics.append(old_cleanup)
with zipfile.ZipFile(D / 'diagnostics-and-trial.zip', 'x', zipfile.ZIP_DEFLATED, compresslevel=6) as z:
    for p in diagnostics: z.write(p, str(p.relative_to(R)))
for root in ['bin/linux-gcc/native-original-input-ownership', 'bin/linux-gcc/native-retained-ownership']:
    products.update(str(p.relative_to(R)) for p in (R / root).rglob('*') if p.is_file())
paths = {}; unique = {}
for key in sorted(products):
    p = R / key; record = pin(p); record['mode'] = stat.S_IMODE(p.stat().st_mode)
    if p.is_symlink(): record['symlink'] = os.readlink(p)
    paths[key] = record; unique.setdefault(record['sha256'], p)
chunks = []; chunk = []; size = 0
for digest, p in sorted(unique.items()):
    if chunk and size + p.stat().st_size > 80000000: chunks.append(chunk); chunk = []; size = 0
    chunk.append((digest, p)); size += p.stat().st_size
if chunk: chunks.append(chunk)
blobs = {}
for index, chunk in enumerate(chunks, 1):
    name = f'actual-products-{index:02}.zip'
    with zipfile.ZipFile(D / name, 'x', zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for digest, p in chunk:
            z.write(p, 'sha256/' + digest); blobs[digest] = {'archive': name, 'member': 'sha256/' + digest}
    assert (D / name).stat().st_size < 95000000
archives = {p.name: pin(p) for p in D.glob('*.zip')}
manifest = {'candidate_source': SOURCE, 'archives': archives,
    'diagnostics': {str(p.relative_to(R)): pin(p) for p in diagnostics},
    'actual_product_paths': paths, 'content_addressed_blobs': blobs,
    'reconstruction': 'Verify archive SHA256; verify each selected blob SHA256/length; restore regular paths/modes then relative symlinks. All members are authentic bytes.',
    'external_system_inputs': 'Exact original before/after manifests preserved; external system/toolchain binaries are not all bundled.'}
write(D / 'MANIFEST.json', manifest)
for name in archives:
    with zipfile.ZipFile(D / name) as z:
        assert z.testzip() is None
        if name == 'diagnostics-and-trial.zip':
            for key, v in manifest['diagnostics'].items():
                data = z.read(key); assert len(data) == v['bytes'] and hashlib.sha256(data).hexdigest() == v['sha256']
        else:
            for key in z.namelist(): assert hashlib.sha256(z.read(key)).hexdigest() == key.split('/')[1]
print(json.dumps({'bundle': str(D.relative_to(R)), 'owners': len(owners), 'all22': True, 'all8': all8,
    'cold_frame': frame['accepted'], 'actual_paths': len(paths), 'archives': archives}))

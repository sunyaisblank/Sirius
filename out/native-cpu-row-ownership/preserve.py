"""Preserve terminal evidence and new actual products, reusing verified Git blobs."""
from pathlib import Path
import hashlib, json, os, re, stat, subprocess, zipfile

R = Path.cwd(); W = Path(__file__).resolve().parent
SOURCE = 'bfa095baf0ff9592cf648cfd49a73d5f6814c05c'
D = R / 'attestations/qualification/bfa095ba/native-cpu-row-ownership-finite-outcome'
PRIOR = R / 'attestations/qualification/08add7e5/native-original-input-ownership-finite-outcome'
PRIOR_COMMIT = 'e4a30cfd64b6dffceb9b43690343a15e484ea8db'

def doc(p): return json.loads(p.read_text())
def pin(p):
    return {'bytes': p.stat().st_size,
            'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()}
def write(p, value):
    assert not p.exists(), p
    p.write_text(json.dumps(value, indent=2) + '\n')

assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == SOURCE
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
assert not subprocess.check_output(['git', 'diff', '0109ab462edade524667d20303507aa09adcbaf6',
    '--', 'src', 'scripts', 'tests', '.github', 'CMakeLists.txt', 'CMakePresets.json']).strip()
assert not D.exists()
owners = {}; products = set()
for p in sorted(W.rglob('owner.json')):
    if 'provider-history' in p.parts: continue
    o = doc(p)
    assert o['status'] == 'terminal' and not o['cleanup_errors'] and not o['remaining_owned_processes']
    assert o['source_unchanged'] and o['whole_inputs_unchanged']
    assert o['observed_births'] and all(x['absent'] for x in o['observed_births'])
    for birth in o['observed_births']:
        live = Path('/proc') / str(birth['pid']) / 'stat'
        if live.exists():
            assert live.read_text().rsplit(')', 1)[1].split()[19] != birth['saved_start_ticks']
    before = doc(p.parent / 'inputs-before.json')
    assert before == doc(p.parent / 'inputs-after.json')
    for key, value in before.items():
        if key.startswith('bin/'):
            assert pin(R / key) == {k: value[k] for k in ['bytes', 'sha256']}, key
            products.add(key)
    owners[str(p.relative_to(R))] = {
        'identity': pin(p), 'passed': o['passed'], 'source': o['source'],
        **{k: o[k] for k in ['stop_reason', 'child_exit', 'wall_seconds',
                             'peak_sampled_rss_kib', 'observed_births']}}

benefit = doc(W / 'readers/benefit.json')
decimal = doc(W / 'readers/decimal-agreement.json')
assert benefit['audit_passed'] and benefit['checks_total'] == 22
all22 = benefit['checks_passed'] == 22
assert decimal['agrees_on_all22'] and decimal['all22_pass'] == all22
controls = {}
for i in range(1, 9):
    p = W / f'runtime/control-{i}/owner.json'
    controls[str(i)] = ({k: doc(p)[k] for k in ['passed', 'wall_seconds', 'stop_reason',
        'child_exit', 'peak_sampled_rss_kib']} if p.exists() else {'status': 'not_run'})
all8 = (W / 'readers/controls-1-8-acceptance.json').exists()
if all8:
    assert all22 and doc(W / 'readers/controls-1-8-acceptance.json')['all8_original_controls_pass']
if not all22:
    assert all(controls[str(i)] == {'status': 'not_run'} for i in range(1, 8))

frame = {'attempted': False, 'accepted': False,
         'reason': 'All22 benefit gates and all8 original controls required first.'}
F = W / 'cold-frame/frame'
if (F / 'owner.json').exists():
    assert all8
    o = doc(F / 'owner.json')
    stdout = re.sub(r'\x1b\[[0-?]*[ -/]*[@-~]', '', (F / 'stdout.log').read_text())
    stderr = (F / 'stderr.log').read_text()
    progress = re.findall(r'Tile (\d+)/(\d+)', stdout)
    scenes = [json.loads(x) for x in re.findall(r'(?:Source scene|Scene) evidence: (\{.*\})', stdout)]
    assert len(scenes) == 2 and scenes[0] == scenes[1]
    assert scenes[0]['width'] == scenes[0]['height'] == 128
    assert scenes[0]['samples_per_pixel'] == 1
    snapshots = re.findall(r'^\[Vulkan\] Retained progress: (.*)$', stderr, re.M)
    stages = {}
    for name, row in re.findall(r'^\[Vulkan\] Retained cumulative physical stage: (\w+), (.*)$', stderr, re.M):
        stages[name] = {k: float(v) if '.' in v else int(v)
            for k, v in re.findall(r'(\w+)=([\d.]+)', row)}
    frame = {'attempted': True, 'accepted': o['passed'], 'stop_reason': o['stop_reason'],
        'wall_seconds': o['wall_seconds'], 'peak_sampled_rss_kib': o['peak_sampled_rss_kib'],
        'last_progress': list(map(int, progress[-1])) if progress else None,
        'XML_present': (F / 'gtest.xml').exists(), 'scene': scenes[0],
        'last_complete_snapshot': snapshots[-1] if snapshots else None,
        'last_complete_snapshot_physical_stages': stages,
        'scope': 'One original changed cold frame only. Partial counters are not completed output, exclusive costs, historical speed comparisons or an extrapolation.'}
    if not o['passed']:
        assert o['stop_reason'] == 'diagnostic_timeout' and o['child_exit'] == -15
        assert not (F / 'gtest.xml').exists()
    write(W / 'cold-frame/result.json', frame)

code = doc(W / 'ordinary/ordinary-code-review.json')
emission = doc(W / 'emission-result.json')
findings = {'candidate_source': SOURCE,
    'public_source_unchanged': '0109ab462edade524667d20303507aa09adcbaf6',
    'all22_original_benefit_gates_pass': all22,
    'benefit_checks_passed': benefit['checks_passed'], 'independent_Decimal_agrees': True,
    'sole_ABBA': {'order': ['a1', 'b1', 'b2', 'a2'], 'internal_repeats_per_slot': 3, 'extra_retries': 0},
    'failed_checks': [x for x in benefit['checks'] + benefit['accounting_checks'] if not x['pass']],
    'whole_call_gains_percent': [{k: x[k] for k in ['direction', 'family', 'improvement_percent']}
        for x in benefit['checks'] if x['scope'] == 'complete_call_ms'],
    'minimum_inclusive_stage_gain_percent': min(x['improvement_percent']
        for x in benefit['checks'] if x['scope'] != 'complete_call_ms'),
    'provisioning_inclusive_API_gains_percent': [x['improvement_percent']
        for x in benefit['accounting_checks']],
    'all32_control8_and_ABBA_readbacks_match_originals': True,
    'original_controls': controls, 'all8_original_controls_admitted': all8,
    'cold_frame': frame,
    'actual_new_builds': {'shader_compilations': 3, 'CXX_compilations': 1,
        'backend_links': 1, 'frame_links': int((W / 'cold-frame/build-result.json').exists())},
    'actual_ordinary_called_function_bytes': {k: x['total_called_function_bytes']
        for k, x in code['rows'].items()},
    'actual_per_invocation_stack_vector_bytes': {k: x['scratch_stack_vector_bytes']
        for k, x in code['rows'].items()},
    'actual_fences_and_coroutine_resumes_zero': all(x['ordinary_IR_fences'] ==
        x['machine_full_fences'] == x['coroutine_resume_functions'] == 0 for x in code['rows'].values()),
    'private_storage_and_lane0_emission': emission,
    'source_or_provider_adopted_publicly': False, 'full_goal_complete': False,
    'scope': 'Private finite evidence. All original workgroups, producers, fields, arithmetic, ordered refusal prefixes, admission and assertions retained. Full science, July image identity, original cold-frame acceptance, native-window/Radeon scenes, Mandatory and exact-main installed-release qualification remain separate obligations.'}
D.mkdir(parents=True)
write(D / 'findings.json', findings)
write(D / 'ownership-closure.json', {'all_owned_runs_terminal_and_reaped': True,
    'owners_count': len(owners), 'owners': owners})
diagnostics = sorted(p for p in W.rglob('*') if p.is_file())
with zipfile.ZipFile(D / 'diagnostics-and-trial.zip', 'x', zipfile.ZIP_DEFLATED, compresslevel=6) as z:
    for p in diagnostics: z.write(p, str(p.relative_to(R)))
for root in ['bin/linux-gcc/native-cpu-row-ownership', 'bin/linux-gcc/native-retained-ownership']:
    products.update(str(p.relative_to(R)) for p in (R / root).rglob('*') if p.is_file())
paths = {}; unique = {}
for key in sorted(products):
    p = R / key; record = pin(p); record['mode'] = stat.S_IMODE(p.stat().st_mode)
    if p.is_symlink(): record['symlink'] = os.readlink(p)
    paths[key] = record; unique.setdefault(record['sha256'], p)

# Reuse authentic, already-pushed product archives; avoid archiving identical binaries again.
prior = doc(PRIOR / 'MANIFEST.json'); blobs = {}; reused_archives = {}
for digest in list(unique):
    if digest not in prior['content_addressed_blobs']: continue
    old = prior['content_addressed_blobs'][digest]
    archive = PRIOR / old['archive']; key = str(archive.relative_to(R))
    if key not in reused_archives:
        assert pin(archive) == prior['archives'][old['archive']]
        git_bytes = subprocess.check_output(['git', 'show', PRIOR_COMMIT + ':' + key])
        assert len(git_bytes) == archive.stat().st_size
        assert hashlib.sha256(git_bytes).hexdigest() == prior['archives'][old['archive']]['sha256']
        reused_archives[key] = {'identity': pin(archive), 'preserved_at_commit': PRIOR_COMMIT}
    with zipfile.ZipFile(archive) as z:
        data = z.read(old['member'])
        assert len(data) == unique[digest].stat().st_size and hashlib.sha256(data).hexdigest() == digest
    blobs[digest] = {'archive': key, 'member': old['member'], 'reused': True}
    del unique[digest]

chunks = []; chunk = []; size = 0
for digest, p in sorted(unique.items()):
    if chunk and size + p.stat().st_size > 80000000:
        chunks.append(chunk); chunk = []; size = 0
    chunk.append((digest, p)); size += p.stat().st_size
if chunk: chunks.append(chunk)
for index, chunk in enumerate(chunks, 1):
    name = f'new-actual-products-{index:02}.zip'
    with zipfile.ZipFile(D / name, 'x', zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for digest, p in chunk:
            z.write(p, 'sha256/' + digest)
            blobs[digest] = {'archive': name, 'member': 'sha256/' + digest, 'reused': False}
    assert (D / name).stat().st_size < 95000000
archives = {p.name: pin(p) for p in D.glob('*.zip')}
manifest = {'candidate_source': SOURCE, 'archives': archives, 'reused_archives': reused_archives,
    'diagnostics': {str(p.relative_to(R)): pin(p) for p in diagnostics},
    'actual_product_paths': paths, 'content_addressed_blobs': blobs,
    'reconstruction': 'Verify archive SHA256 and selected blob SHA256/length; restore paths/modes then relative symlinks. New archive names are bundle-relative; reused names are checkout-relative and pinned to their original pushed Git commit. All members are authentic bytes.',
    'external_system_inputs': 'Exact original before/after manifests retained; external system/toolchain binaries are not all bundled.'}
write(D / 'MANIFEST.json', manifest)
for name in archives:
    with zipfile.ZipFile(D / name) as z:
        assert z.testzip() is None
        if name == 'diagnostics-and-trial.zip':
            for key, value in manifest['diagnostics'].items():
                data = z.read(key)
                assert len(data) == value['bytes'] and hashlib.sha256(data).hexdigest() == value['sha256']
        else:
            for key in z.namelist(): assert hashlib.sha256(z.read(key)).hexdigest() == key.split('/')[1]
assert set(blobs) == {value['sha256'] for value in paths.values()}
print(json.dumps({'bundle': str(D.relative_to(R)), 'owners': len(owners), 'all22': all22,
    'all8': all8, 'cold_frame': frame['accepted'], 'actual_paths': len(paths),
    'reused_archives': len(reused_archives), 'new_archives': archives}))

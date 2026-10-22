"""Preserve completed runs and actual products; never execute a workload."""
from pathlib import Path
import hashlib, json, os, re, stat, subprocess, zipfile
R = Path.cwd(); W = Path(__file__).resolve().parent
D = R / 'attestations/qualification/d9dbd657/native-phase-sharing-finite-pass-frame-timeout'
def doc(p): return json.loads(p.read_text())
def pin(p):
    return {'bytes': p.stat().st_size,
            'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()}
def write(p, value):
    assert not p.exists(), p
    p.write_text(json.dumps(value, indent=2) + '\n')
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == 'd9dbd65796aab5dd36e17481198ae1f6b6201a2c'
assert subprocess.check_output(['git', 'status', '--porcelain'], text=True) == ''
assert not D.exists(); D.mkdir(parents=True)
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
            fields = live.read_text().rsplit(')', 1)[1].split()
            assert fields[19] != birth['saved_start_ticks'], birth
    before = doc(p.parent / 'inputs-before.json')
    assert before == doc(p.parent / 'inputs-after.json')
    for k, v in before.items():
        if k.startswith('bin/') and (R / k).is_file():
            assert pin(R / k) == {x: v[x] for x in ['bytes', 'sha256']}, k
            products.add(k)
    owners[str(p.relative_to(R))] = {'identity': pin(p), 'passed': o['passed'],
        'source': o['source'], 'stop_reason': o['stop_reason'], 'child_exit': o['child_exit'],
        'wall_seconds': o['wall_seconds'], 'peak_sampled_rss_kib': o['peak_sampled_rss_kib'],
        'observed_births': o['observed_births']}
assert len(owners) == 17 and sum(x['passed'] for x in owners.values()) == 16
write(D / 'ownership-closure.json', {'all17_terminal_and_reaped': True, 'owners': owners,
    'reused2_provider_owners_preserved_separately': True})

F = W / 'cold-frame/frame'; o = doc(F / 'owner.json')
assert not o['passed'] and o['stop_reason'] == 'diagnostic_timeout' and o['child_exit'] == -15
assert o['limits'] == {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152}
assert 300 <= o['wall_seconds'] < 302
assert o['provider_mappings'] and all(x['owned_session'] == o['child_pid'] for x in o['provider_mappings'])
restored = doc(W / 'library-restoration/identity-restoration.json')
for m in o['provider_mappings']:
    aliases = m['qualified_runtime_aliases']; a = aliases[str(R / restored['private_restoration'])]
    assert {k: a['identity'][k] for k in ['bytes', 'sha256']} == restored['historical_identity']
    assert a['maps']
stdout = re.sub(r'\x1b\[[0-?]*[ -/]*[@-~]', '', (F / 'stdout.log').read_text())
stderr = (F / 'stderr.log').read_text()
progress = re.findall(r'Tile (\d+)/(\d+)', stdout)
assert progress[-1] == ('1050', '16384')
scene = [json.loads(x) for x in re.findall(r'(?:Source scene|Scene) evidence: (\{.*\})', stdout)]
assert len(scene) == 2 and scene[0] == scene[1]
assert scene[0]['width'] == scene[0]['height'] == 128 and scene[0]['samples_per_pixel'] == 1
snapshots = re.findall(r'^\[Vulkan\] Retained progress: (.*)$', stderr, re.M)
assert snapshots[-1].startswith('862 work tiles of 1px, 678 camera batches, 178819 accepted intervals, 5065 rejected intervals,')
stage = {}
for name, data in re.findall(r'^\[Vulkan\] Retained cumulative physical stage: (\w+), (.*)$', stderr, re.M):
    stage[name] = {k: float(v) if '.' in v else int(v)
                   for k, v in re.findall(r'(\w+)=([\d.]+)', data)}
assert set(stage) == {'camera', 'transport', 'endpoint', 'dense', 'initialize', 'ray_camera', 'dopri_phase'}
assert not (F / 'gtest.xml').exists()
output = R / 'renders/native-endpoint-phase-d9dbd657/cold-frame-tmp/sirius_gpu_parse.ppm'
assert not output.exists()
cold = {'accepted': False, 'outcome': 'diagnostic_timeout', 'owner': pin(F / 'owner.json'),
    'wall_seconds': o['wall_seconds'], 'peak_sampled_rss_kib': o['peak_sampled_rss_kib'],
    'last_reported_tiles': 1050, 'required_tiles': 16384, 'completed_frame_or_XML': False,
    'last_complete_snapshot': snapshots[-1], 'last_complete_snapshot_physical_stages': stage,
    'scene': scene[0], 'same_child_qualified_library_mapping_verified': True,
    'source_and_inputs_unchanged': True, 'owned_processes_reaped': True,
    'scope': 'One original cold frame, failed. Last complete counters cover 862 tiles; final progress covers 1050. Submit/wait includes driver/scheduling/JIT costs. No partial-frame extrapolation, exclusive-time subtraction, retry or scientific/release admission.'}
write(W / 'cold-frame/result.json', cold)
benefit = doc(W / 'readers/benefit.json'); controls = doc(W / 'readers/controls-1-8-acceptance.json')
assert controls['all8_original_controls_pass'] and benefit['checks_passed'] == benefit['checks_total'] == 22
findings = {'candidate_source': 'd9dbd65796aab5dd36e17481198ae1f6b6201a2c',
    'public_source_unchanged': '0109ab462edade524667d20303507aa09adcbaf6',
    'all8_original_finite_controls_pass': True, 'all22_original_benefit_gates_pass': True,
    'sole_ABBA': {'order': benefit['fixed_order'], 'internal_repeats_per_slot': 3, 'extra_retries': 0},
    'whole_call_gains_percent': [{k: x[k] for k in ['direction', 'family', 'improvement_percent']}
        for x in benefit['checks'] if x['scope'] == 'complete_call_ms'],
    'inclusive_stage_checks': 16,
    'minimum_inclusive_stage_gain_percent': min(x['improvement_percent'] for x in benefit['checks'] if x['scope'] != 'complete_call_ms'),
    'provisioning_inclusive_API_gains_percent': [x['improvement_percent'] for x in benefit['accounting_checks']],
    'all32_control8_readbacks_match_original_baselines': True,
    'factored_control1_actual_routes': 2, 'factored_control1_complete_readbacks': 28,
    'ordinary_endpoint_called_function_bytes': 662964,
    'ordinary_endpoint_previous_candidate_bytes': 1287383,
    'actual_new_builds': {'shader_compilations': 1, 'CXX_compilations': 1, 'backend_test_links': 1, 'frame_test_links': 1},
    'reused_DopriPhase_not_compiled_again': True,
    'cold_frame': cold,
    'postexecution_OS_library_update': restored,
    'adopted_into_public_source': False, 'full_goal_complete': False,
    'remaining': 'Original cold frame; full science, original July image identity, native-window/Radeon scenes, Mandatory and installed-release qualification on exact clean pushed main.',
    'next_change_basis': 'Captured transport/endpoint submit-wait dominates the last complete snapshot. Transport normalized input caching already exists in rejected3045067 history. Endpoint ordered-prefix ownership plus the new actual IEEE32 provider/source context is distinct; no unchanged rejected retry or claimed speedup.'}
write(D / 'findings.json', findings)

# Preserve all diagnostic bytes, including failed reader forms and their corrections.
diagnostics = sorted(p for p in W.rglob('*') if p.is_file())
with zipfile.ZipFile(D / 'diagnostics-and-trial.zip', 'x', zipfile.ZIP_DEFLATED, compresslevel=6) as z:
    for p in diagnostics: z.write(p, str(p.relative_to(R)))

# Preserve actual consumed products; content-addressed blobs avoid duplicating
# identical resources. The manifest reconstructs original paths/modes/symlinks.
for root in ['bin/linux-gcc/native-endpoint-phase-sharing', 'bin/linux-gcc/native-retained-ownership']:
    products.update(str(p.relative_to(R)) for p in (R / root).rglob('*') if p.is_file())
paths = {}; unique = {}
for key in sorted(products):
    p = R / key; record = pin(p)
    record['mode'] = stat.S_IMODE(p.stat().st_mode)
    if p.is_symlink(): record['symlink'] = os.readlink(p)
    paths[key] = record; unique.setdefault(record['sha256'], p)
chunks = []; chunk = []; size = 0
for digest, p in sorted(unique.items()):
    if chunk and size + p.stat().st_size > 80000000:
        chunks.append(chunk); chunk = []; size = 0
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
manifest = {'candidate_source': findings['candidate_source'], 'archives': archives,
    'diagnostics': {str(p.relative_to(R)): pin(p) for p in diagnostics},
    'actual_product_paths': paths, 'content_addressed_blobs': blobs,
    'reconstruction': 'Verify archive SHA256, read each manifest-selected blob, verify its SHA256/length, restore each regular path/mode, then restore relative symlinks. Archive entries are authentic bytes; no synthetic products.',
    'external_system_and_toolchain_inputs': 'Exact original before/after input manifests preserved in diagnostics; only the restored historical liblzma and actual private providers are bundled as runtime products.'}
write(D / 'MANIFEST.json', manifest)
# Read every emitted member back and verify against its original bytes.
for name in archives:
    with zipfile.ZipFile(D / name) as z:
        assert z.testzip() is None
        if name == 'diagnostics-and-trial.zip':
            for key, v in manifest['diagnostics'].items():
                data = z.read(key); assert len(data) == v['bytes'] and hashlib.sha256(data).hexdigest() == v['sha256']
        else:
            for key in z.namelist(): assert hashlib.sha256(z.read(key)).hexdigest() == key.split('/')[1]
print(json.dumps({'bundle': str(D.relative_to(R)), 'owners_terminal': 17, 'finite_controls': 8,
    'benefit_gates': 22, 'cold_frame_passed': False, 'actual_product_paths': len(paths),
    'unique_product_blobs': len(blobs), 'archives': archives}))

"""Preserve the completed, rejected ownership trial and its actual products."""
from pathlib import Path
import hashlib, json, os, subprocess, zipfile
R = Path.cwd(); W = R / 'out/native-retained-ownership'
B = R / 'bin/linux-gcc/native-retained-ownership'
P = R / 'attestations/qualification/a2029467/native-ownership-rejected'
SOURCE = 'a2029467924854a268a2db6b8ce859b379df6222'
def write(p, v): p.write_text(json.dumps(v, indent=2, sort_keys=True) + '\n')
def pin(p):
    return {'bytes': p.stat().st_size,
            'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest(),
            'mode': p.stat().st_mode & 0o777}
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == SOURCE
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
assert not P.exists(); P.mkdir(parents=True)
owners = {}
paths = [W / 'modules/owner.json', W / 'matched/owner.json',
         W / 'ordinary/endpoint/owner.json', W / 'ordinary/dopri/owner.json']
paths += [W / 'runtime' / n / 'owner.json' for n in ['control-8', 'a1', 'b1', 'b2', 'a2']]
for p in paths:
    v = json.loads(p.read_text())
    assert v['passed'] and v['status'] == 'terminal' and v['child_exit'] == 0
    assert v['source_unchanged'] and v['whole_inputs_unchanged']
    assert not v['remaining_owned_processes'] and not v['cleanup_errors']
    for birth in v['observed_births']:
        assert birth['absent']
        f = Path('/proc') / str(birth['pid']) / 'stat'
        try: assert f.read_text().rsplit(')', 1)[1].split()[19] != birth['saved_start_ticks']
        except FileNotFoundError: pass
    owners[str(p.relative_to(R))] = {'identity': pin(p), 'observed_births': v['observed_births'],
        'wall_seconds': v['wall_seconds'], 'peak_sampled_rss_kib': v['peak_sampled_rss_kib'],
        'source': v['source']}
write(P / 'ownership-closure.json', {'all9_terminal_and_reaped': True, 'owners': owners})

def pack(path, files):
    records = {}
    with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for p in sorted(set(files)):
            assert p.is_file() and not p.is_symlink()
            key = str(p.relative_to(R)); records[key] = pin(p); z.write(p, key)
    with zipfile.ZipFile(path) as z:
        assert set(z.namelist()) == set(records) and z.testzip() is None
        for key, record in records.items():
            b = z.read(key)
            assert len(b) == record['bytes'] and hashlib.sha256(b).hexdigest() == record['sha256']
    return {'identity': pin(path), 'files': records}

diagnostics = [p for p in W.rglob('*') if p.is_file()] + [Path(__file__).resolve()]
products = [p for p in B.rglob('*') if p.is_file()]
products += [p for p in (R / 'bin/linux-gcc/src/sirius/backend/retained').iterdir() if p.is_file()]
containers = {'diagnostics-and-trial.zip': pack(P / 'diagnostics-and-trial.zip', diagnostics),
              'actual-products.zip': pack(P / 'actual-products.zip', products)}
v = json.loads((W / 'readers/benefit.json').read_text())
d = json.loads((W / 'readers/decimal-agreement.json').read_text())
assert v['audit_passed'] and not v['benefit_passed'] and v['checks_total'] == 22 and v['checks_passed'] == 19
assert d['agrees_on_all22'] and not d['all22_pass']
assert not any((W / 'runtime' / ('control-' + str(i))).exists() for i in range(1, 8))
checks = v['checks'] + v['accounting_checks']
write(P / 'findings.json', {
    'trial_source': SOURCE, 'public_source': '0109ab462edade524667d20303507aa09adcbaf6',
    'decision': 'Rejected, unadopted, no unchanged retry or remaining7/cold-frame execution.',
    'implementation': 'One complete lane-zero curved native Endpoint projection with all7 published words; all362 DopriPhase validations distributed across the existing64-lane workgroup.',
    'all35_original_modules_and_tables_unchanged': True,
    'original_native_Project_function_equal_after_numeric_ID_renaming': True,
    'actual_ordinary_called_code_bytes': {'Endpoint': 1287383, 'DopriPhase': 111219},
    'all_original_expanded_shader_IR_machine_fences': {'Endpoint': 11, 'DopriPhase': 6},
    'original_C8_and_all32_complete_cross_arm_readbacks_pass': True,
    'benefit': {'total': 22, 'passed': 19, 'whole_call_passed': 3,
        'inclusive_stages_passed': 14, 'provisioning_inclusive_API_passed': 2,
        'independent_Decimal_agrees_on_all22': True, 'checks': checks},
    'failed': [c for c in checks if not c['pass']],
    'provider_implementation_reuse': 'Unchanged exact scoped IEEE32 provider and finite primitive/boundary proof preserved at6a6bb1ca; no rebuild or scalar rerun.',
    'full_scientific_Mandatory_runtime_release_qualification': False,
    'next_source_basis': 'Endpoint still emits1287383 bytes and has two original Project and EvaluateProgram call sites. Test whether one common phase loop preserves complete ownership and removes duplicated compiled bodies before another integrated trial.'})
write(P / 'MANIFEST.json', {'schema': 'sirius-native-ownership-single-rejected-trial-v1',
    'source': SOURCE, 'adopted': False, 'containers': containers,
    'other_files': {p.name: pin(p) for p in P.iterdir() if p.is_file() and p.suffix != '.zip'}})
print(json.dumps({'bundle': str(P.relative_to(R)), 'owners': len(owners),
    'containers': {k: v['identity'] for k, v in containers.items()}}, indent=2))

"""Preserve a source-unchanged recovery review; never issue qualification."""
from pathlib import Path
import collections
import datetime
import hashlib
import json
import math
import os
import statistics
import subprocess
import xml.etree.ElementTree as ET

ROOT = Path.cwd()
OUTPUT = ROOT / 'attestations/fresh-session-review/4bcd987/2026-10-14'
REV = '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc'
BRANCH = 'development/retained-dense-step-control-2026-10-06'
TAG = 'evidence/fresh-session-review-4bcd987-2026-10-14'

def run(argv, env=None):
    p = subprocess.run(argv, cwd=ROOT, env=env, capture_output=True, text=True)
    return {'argv': argv, 'exit': p.returncode, 'stdout': p.stdout, 'stderr': p.stderr}

def checked(argv):
    r = run(argv)
    assert r['exit'] == 0, r
    return r['stdout'].strip()

def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def load(path):
    return json.loads(Path(path).read_text())

assert checked(['git', 'rev-parse', 'HEAD']) == REV
assert checked(['git', 'branch', '--show-current']) == BRANCH
assert checked(['git', 'status', '--porcelain=v1']) == ''
model = load('tests/operating_model.json')
checks = []
for args in [
    ['python3', '-B', 'scripts/generate-ctest-labels.py', '--check'],
    ['python3', '-B', 'scripts/verify-operating-model.py'],
    ['python3', '-B', 'scripts/verify-repository-structure.py'],
    ['python3', '-B', 'scripts/verify-attestation.py', '--source-root', str(ROOT),
     'attestations/native-build/4bcd987/windows-build/windows-build.json'],
    ['python3', '-B', 'scripts/verify-attestation.py', '--source-root', str(ROOT),
     'attestations/native-build/4bcd987/macos-build/macos-build.json'],
    ['python3', '-B', 'scripts/verify-alignment.py', '--mode', 'qualification',
     '--source-revision', REV, '--source-tree-clean', 'true', '--operating-model',
     'tests/operating_model.json', '--attestation-root', 'attestations/native-build/4bcd987',
     '--source-root', str(ROOT), '--output', 'bin/linux-gcc/generated/sirius/alignment_receipt.json', '--check'],
]:
    r = run(args)
    assert r['exit'] == 0, r
    checks.append(r)

app = 'bin/linux-gcc/src/sirius/app/sirius'
readiness = run([app, 'info', 'readiness', '--json'])
system = run([app, 'info', 'system', '--json'])
dzn_env = dict(os.environ, VK_DRIVER_FILES=str(ROOT / 'bin/toolchains/mesa-25.2.8/install/share/vulkan/icd.d/dzn_icd.x86_64.json'))
physical = run([app, 'info', 'system', '--json'], dzn_env)
for r in [readiness, system, physical]:
    r['parsed'] = json.loads(r['stdout'])
    checks.append(r)
assert not readiness['parsed']['ready']
assert not readiness['parsed']['release_ready']
assert not readiness['parsed']['evidence_generation_ready']
assert not readiness['parsed']['mandatory_build_gate']['valid']
assert physical['parsed']['wsl2']
dev = physical['parsed']['backends']['vulkan']['devices'][0]
assert dev['kind'] != 'software' and '780M' in dev['name'] and dev['driver_name'] == 'Dozen'

inventory = json.loads(checked(['ctest', '--test-dir', 'bin/linux-gcc', '--show-only=json-v1']))
enabled = [t for t in inventory['tests'] if not any(p['name'] == 'DISABLED' and p['value'] for p in t['properties'])]
mandatory = [t for t in enabled if any(p['name'] == 'LABELS' and 'Mandatory' in p['value'] for p in t['properties'])]
assert len(enabled) == len(mandatory) == 1147
registration = {'registered': len(inventory['tests']), 'enabled': len(enabled), 'mandatory': len(mandatory),
                'execution_claim': False, 'test_names': sorted(t['name'] for t in enabled)}

manifest_checks = []
bundles = [
    ('attestations/scientific-core/4bcd987/polarised-scientific-witnesses', 'manifest.json'),
    ('attestations/scientific-core/4bcd987/viewer-and-polarisation', 'manifest.json'),
    ('attestations/alignment/4bcd987/gcc-revision-synchronisation', 'manifest.json'),
    ('attestations/software-vulkan/b34d8c2/cpu-dense-growth-only-integration', 'sha256-manifest.json'),
    ('attestations/software-vulkan/ea0ae34/endpoint-prefix-production', 'sha256-manifest.json'),
    ('attestations/software-vulkan/4bcd987/guarded-product-high-benefit', 'whole-byte-manifest.json'),
]
current_witnesses = collections.defaultdict(list)
xml_records = []
for folder, name in bundles:
    root = Path(folder)
    files = load(root / name)['files']
    entries = files if isinstance(files, list) else [dict(v, path=k) for k, v in files.items()]
    for e in entries:
        p = root / e['path']
        if not p.is_file():
            p = root / 'files' / e['path']
        assert p.is_file() and p.stat().st_size == e['bytes'] and sha(p) == e['sha256'], str(p)
    manifest_checks.append({'manifest': str(root / name), 'sha256': sha(root / name),
                            'files_verified': len(entries), 'mismatches': 0})
    for p in root.rglob('*.xml'):
        try:
            r = ET.parse(p).getroot()
        except ET.ParseError:
            continue
        if r.tag not in ('testsuite', 'testsuites'):
            continue
        cases = r.findall('.//testcase')
        good, bad = [], []
        for c in cases:
            case = c.attrib.get('name', '')
            if '.' not in case and c.attrib.get('classname'):
                case = c.attrib['classname'] + '.' + case
            (bad if any(c.find(t) is not None for t in ['failure', 'error', 'skipped']) else good).append(case)
        xml_records.append({'path': str(p), 'sha256': sha(p), 'summary': r.attrib,
                            'passed_names': good, 'failed_or_skipped_names': bad})
        # Reuse only reviewed current-source records; never historical or candidate trials.
        if '/4bcd987/' in folder and 'guarded-product' not in folder and not bad:
            for case in good:
                current_witnesses[case].append(str(p))

auth = load('attestations/alignment/4bcd987/gcc-revision-synchronisation/files/out/gcc-revision-synchronisation/current-product-authority.json')
live_bindings = {}
for group, records in auth['actual_artifact_records'].items():
    for e in records.values():
        p = (Path('bin/linux-gcc') if e['root'] == 'build' else Path('.')) / e['path']
        assert p.is_file() and p.stat().st_size == e['bytes'] and sha(p) == e['sha256'], str(p)
    live_bindings[group] = {'count': len(records), 'mismatches': 0, 'records': records}

performance_root = Path('attestations/software-vulkan/4bcd987/guarded-product-high-benefit/files/out/guarded-product-high-benefit/benefit')
comparison = load(performance_root / 'software-comparison.json')
props = {}
for slot in ['a1', 'b1', 'b2', 'a2']:
    r = ET.parse(performance_root / slot / 'gtest.xml').getroot()
    assert all(not any(c.find(t) is not None for t in ['failure', 'error', 'skipped']) for c in r.findall('.//testcase'))
    props[slot] = {e.attrib['name']: e.attrib['value'] for e in r.findall('.//property')}
def raw_samples(slot, family, scope):
    prefix = 'legacy' if family == 'legacy' else 'coupled'
    p = props[slot]
    if scope == 'complete_call_ms':
        return [float(p[f'{prefix}_repeat{i}_complete_call_ms']) for i in range(3)]
    stage = scope.removesuffix(' inclusive Dispatch')
    markers = [dict(item.split('=', 1) for item in v.split(';')) for k, v in p.items() if k.startswith('device_observation_')]
    return [math.fsum(float(m['total_ms']) for m in markers if m['phase'] == f'{prefix}_repeat{i}' and m['stage'] == stage) for i in range(3)]
for c in comparison['checks']:
    a, b = ('a1', 'b1') if c['direction'] == 'forward' else ('a2', 'b2')
    for slot, key in [(a, 'baseline_samples_ms'), (b, 'candidate_samples_ms')]:
        assert raw_samples(slot, c['family'], c['scope']) == c[key]
    improvement = 100 * (statistics.median(c['baseline_samples_ms']) - statistics.median(c['candidate_samples_ms'])) / statistics.median(c['baseline_samples_ms'])
    assert abs(improvement - c['improvement_percent']) < 1e-8
    assert (improvement >= c['minimum_percent']) == c['pass']
costs = {}
for slot, p in props.items():
    assert float(p['owner_provisioning_ms']) == 0
    costs[slot] = math.fsum([float(p['observed_create_ms']), float(p['disabled_attempt_ms']), float(p['disabled_sample_ms'])] +
                          [float(p[f'{family}_repeat{i}_complete_call_ms']) for family in ['legacy', 'coupled'] for i in range(3)])
for c in comparison['accounting_checks']:
    a, b = ('a1', 'b1') if c['direction'] == 'forward' else ('a2', 'b2')
    assert abs(costs[a] - c['baseline_ms']) < 1e-8 and abs(costs[b] - c['candidate_ms']) < 1e-8
assert sum(c['pass'] for c in comparison['checks']) == 11
assert sum(c['pass'] for c in comparison['accounting_checks']) == 1

issues = [json.loads(checked(['gh', 'issue', 'view', str(n), '--repo', 'sunyaisblank/Sirius', '--json', 'number,title,body,state,url,updatedAt'])) for n in [39, 44, 48, 53, 98]]
pr = json.loads(checked(['gh', 'pr', 'view', '62', '--repo', 'sunyaisblank/Sirius', '--json', 'number,title,body,state,isDraft,headRefName,headRefOid,baseRefName,url,statusCheckRollup']))
assert pr['isDraft'] and pr['headRefOid'] == REV
ci = [json.loads(checked(['gh', 'run', 'view', str(n), '--repo', 'sunyaisblank/Sirius', '--json', 'jobs,headSha,event,status,conclusion'])) for n in [38318828623, 38336123427]]
active = [json.loads(checked(['gh', 'run', 'list', '--repo', 'sunyaisblank/Sirius', '--status', status, '--limit', '100', '--json', 'databaseId,headSha,status'])) for status in ['in_progress', 'queued']]
assert active == [[], []]
remote_refs = checked(['git', 'ls-remote', 'origin', 'refs/heads/main', 'refs/heads/' + BRANCH])
assert REV in remote_refs and '5e706073a3fd2f9afcf8efbcde2f094122050f6a' in remote_refs
assert checked(['git', 'rev-parse', REV + '^{tree}']) == checked(['git', 'rev-parse', '2c934151484e7df886221eb343ebbdef97e68294^{tree}'])

domains = []
for name in readiness['parsed']['ultimate_ideal']['external_domain_ids_required']:
    admitted = name in readiness['parsed']['ultimate_ideal']['external_domain_ids_admitted']
    status = 'verified_complete' if admitted else 'externally_blocked' if name == 'macos-moltenvk' else 'actionable_now'
    domains.append({'id': name, 'status': status, 'scope': 'Exact-current-revision build-only admission' if admitted else 'Final revision qualification remains required',
                    'reason': 'Independent existing verifier passes' if admitted else 'No accessible native Apple runtime identified; obtain a native host/runner route' if name == 'macos-moltenvk' else 'Windows/WSL2 and physical Radeon 780M are accessible; rendering, full gates and domain-specific preflight remain prerequisites'})

ledger = []
for group in ['acceptance_criteria', 'required_dimensions', 'capability_contracts']:
    for item in model[group]:
        expected = {e['name'] for e in item['evidence']}
        verified = sorted(expected & current_witnesses.keys())
        missing = sorted(expected - current_witnesses.keys())
        if group == 'required_dimensions' and not missing:
            status = 'verified_complete'
            scope = 'Only the complete existing named finite witness set; no full Mandatory/release promotion'
        elif item.get('state') == 'attestation_required' and item['profiles'] == ['macos-native-build', 'macos-moltenvk']:
            status = 'externally_blocked'
            scope = 'Native build verified; final native runtime requires an accessible Apple host'
        else:
            status = 'actionable_now'
            scope = 'Implementation/evidence recovery is possible now; final full criterion is not proved'
        ledger.append({'group': group, 'id': item['id'], 'status': status, 'scope': scope,
                       'original_model_requirement': {k: v for k, v in item.items() if k != 'evidence'},
                       'required_evidence': item['evidence'], 'reviewed_current_witnesses': verified,
                       'not_proved_by_reviewed_current_bundles': missing,
                       'evidence_paths': {n: current_witnesses[n] for n in verified},
                       'finish_line': 'Every original applicable witness plus required exact-final integrated/domain gate; no missing witness may be inferred from registration or neighboring tests'})
assert len(ledger) == 64

gates = [
    ('strict_build_gcc_clang', 'actionable_now', 'Current GCC bytes and strict cache are verified; final complete GCC/Clang checks remain'),
    ('integration_boundary', 'verified_complete', 'Both current CI runs pass three integration jobs; four full jobs skipped; merge-tree equality is verified'),
    ('external_execution_preflight', 'actionable_now', 'Fresh non-promoting exact-candidate preflight required before each native/physical campaign'),
    ('consolidated_wsl2_authority', 'actionable_now', 'Reuse only a newly verified same-revision full hardware authority; no such current runtime bundle admitted'),
    ('immutable_build_inputs', 'verified_complete', 'Current repository governance and build-policy source controls pass within source scope'),
    ('complete_mandatory_registration_execution', 'actionable_now', '1147 live enabled Mandatory registrations, no complete current execution receipt'),
    ('kernel_oracle_and_beam_parity', 'actionable_now', 'Finite historical acceptance remains scoped; complete final CPU/device scientific estate required'),
    ('output_decode_transfer_and_roundtrip', 'actionable_now', 'Historical original writer controls verified at b34d8c2; final output estate still required'),
    ('physical_hardware_envelope', 'actionable_now', 'Physical Radeon/Dozen visible; original two full-size governed scenes not qualified'),
    ('qualification_to_release_alignment', 'externally_blocked', 'Only 2/8 domains admitted; all eight plus fresh final release gate required, including unconfirmed Apple route'),
    ('original_forensic_pixel_identity', 'decision_required', '70d156f requires pixel identity; 8bd56ca replaces it with semantic/oracle and output invariants; owner acceptance of this change is not recovered'),
    ('original_software_cold_frame', 'actionable_now', 'Preserved ea0ae34 raw JUnit times out at 300.027 seconds; original required 128x128/1-SPP fp32 cold Kerr command remains unchanged'),
    ('execution_refusal_recovery_rollback', 'actionable_now', 'Preserve complete 40 fields, 4 columns, 7 RHS records, 3 independent RK trials and 4 quartics and all consumer joins; finite controls do not prove full estate'),
    ('cpu_endpoint_ownership_source', 'verified_complete', 'Current accepted final acceleration consumes geometries[6] and final projected tangent; historical source word/query regression remains bounded'),
    ('viewer_committed_frame_source_and_controls', 'verified_complete', 'Current source commits pixels and dimensions under one lock; six current cases pass and predecessor regression fails'),
    ('final_cpu_scientific_integration_issue53', 'actionable_now', '53 historical numerics, four sanitizers, original detector and writer cases are verified at their original revision; current polarisation 10/10 is finite'),
    ('final_main_integration', 'decision_required', 'PR62 remains draft; final merge is an owner decision after criteria pass, preserve clean pushed source and obtain new exact-main evidence'),
    ('post_gate_installed_relocated_readiness', 'actionable_now', 'Fresh final-main install/relocation/readiness and render/view initialisation required; present runtime explicitly refuses'),
    ('closure_scorecard_and_current_evidence_report', 'actionable_now', 'This recovery ledger is not an all-green final scorecard; final six-programme and P/E closure report remains'),
]
programmes = [
    {'id': 1, 'name': 'foundations', 'status': 'verified_complete', 'scope': 'Present specification/style/architecture/presets and structural governance; compiler execution remains a separate gate'},
    *[{'id': n, 'name': name, 'status': 'actionable_now', 'scope': 'Source is present; final integrated scientific/rendering/governor/qualification proof remains'} for n, name in [(2, 'port'), (3, 'kernels'), (4, 'memory'), (5, 'physics')]],
    {'id': 6, 'name': 'closure', 'status': 'actionable_now', 'scope': 'CI and truthful README exist; full same-revision scorecard, release packaging checks and final admission remain'},
]
next_package = {
    'id': 'issue44-exact-power-of-two-divisor-compiler-screen',
    'status': 'actionable_now', 'implementation_started': False,
    'premise': 'After original exceptional guards and unpack, bb.significand == 0x00800000 permits aa.significand << 3 and exponent difference through unchanged PB32Round; every other divisor uses original fallback',
    'source': 'src/sirius/kernels/portable_binary32.h:174',
    'required_preservation': ['exceptional-value precedence', 'sign/subnormal/underflow/tie/overflow handling', 'two separate seed-scaling roundings', 'fallback', 'ABI', 'precision/product modes', 'all scientific states and transactional contracts'],
    'finish_line': ['Freeze stock/candidate complete-source and compiled-module identities', 'Show all four dynamic integer divisions absent on eligible paths in actual ordinary driver execution', 'Prove fallback/ABI/mode preservation and measure mixed-lane/JIT/control/frame countercost', 'Write a justified stop/go decision for the unchanged original numerical and benefit experiments; no renderer adoption or frame-speed claim'],
    'subsequent_adoption_gates': ['all eight unchanged integrated correctness controls', 'original exact arithmetic/provider checks', 'single preregistered A1/B1/B2/A2 with three repeats per slot', 'all four complete-call improvements >=5%', 'all sixteen inclusive stage improvements >=-5%', 'both creation/provisioning-inclusive API improvements >=0%', 'original governors, references, settings, scenes and no unchanged rejected retry'],
    'downstream_finish_line': 'Accepted material runtime improvement, original cold frame <=300s, usable native preview and both original 100000-star moving-ThinLens frames under 2048 MiB, final CPU science and all eight same-revision domains, then owner integration and exact-main release admission',
}
report = {
    'kind': 'fresh-session-recovery-review', 'schema_version': 1,
    'completed_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'review_goal_complete': True, 'system_goal_complete': False,
    'promotable': False, 'qualification_claim': False, 'implementation_changes': False,
    'source_revision': REV, 'source_branch': BRANCH, 'main_revision': '5e706073a3fd2f9afcf8efbcde2f094122050f6a',
    'original_goal_sources': ['70d156f:docs/SPECIFICATION.md', 'docs/SPECIFICATION.md', 'tests/operating_model.json', 'GitHub issues 39/44/48/53 and draft PR62'],
    'handoff': {'path': '/home/astra/.codex/attachments/316f2e8c-6275-4006-9469-9ed50de6fbaa/pasted-text-1.txt', 'sha256': sha(Path('/home/astra/.codex/attachments/316f2e8c-6275-4006-9469-9ed50de6fbaa/pasted-text-1.txt'))},
    'recovered_goal': 'Complete the ground-up C++26 Sirius DNGR parity/exceedance rebuild on CPU and consumer Vulkan across Linux/WSL2/Windows and supported macOS, including progressive Radeon 780M rendering and exact governed 1080p/IMAX star scenes within 2048 MiB, full science/execution contracts and final revision release admission',
    'accepted_working_scope': 'Current specification and original issue acceptance are preserved; finite #39, CPU endpoint reuse, retained ownership and coherent viewer publication are subordinate results. Rejected experiments remain rejected. No reduced scenes, settings, counts, independent references, tolerances, precision eligibility or governors.',
    'material_missing_context': ['Original user conversation and governing Prompt.md sections 15/19 are not recovered; scope reconstruction rests on immutable specification history, issues, source and raw evidence', 'Owner authorization for replacing original pixel-identical forensic comparisons with semantic/decode evidence is not independently recovered', 'No accessible native Apple runtime route identified; that external domain needs a host/runner', 'Raw binary qualification evidence is local ignored data; pushed evidence tags preserve useful tools/findings, not binary backup'],
    'scope_and_claim_audit': [
        {'status': 'decision_required', 'finding': 'Pixel-identity acceptance was explicitly relaxed in 8bd56ca; no owner authorization recovered. Do not silently discard the original criterion or treat decode/oracle controls as pixel identity.'},
        {'status': 'unsupported_historical_completion', 'finding': 'Retired July engagement scorecard HOLDS/physical E3 closure claims do not prove the current original beam-filtered moving-ThinLens full-resolution workload or current revision. Its own later historical warning withdraws it as operational authority.'},
        {'status': 'scope_preserved', 'finding': 'Current #44/#48/#53 and PR62 explicitly remain incomplete; #39 is finite Schwarzschild/Kerr scientific acceptance, not DNGR/full-release closure.'},
        {'status': 'not_completion', 'finding': 'Green CI executes three integration jobs and skips four qualification jobs. 1147 registrations are not executions. Current 10 polarisation dimension witnesses cover only 9/11 E2 witnesses.'},
        {'status': 'explicit_capability_bounds', 'finding': 'Unsupported polarised Vulkan/volume, motion blur, jet/corona/line models and nonstationary beam modes fail closed explicitly. Nine registry entries do not establish unrestricted observable accuracy in all metric families.'},
    ],
    'checks': checks, 'registration': registration, 'manifest_integrity': manifest_checks,
    'raw_junit': xml_records, 'current_live_byte_bindings': live_bindings,
    'performance_rejection': {'issue': 98, 'interval_stage_pass': '11/20', 'inclusive_api_pass': '1/2', 'all_20_sample_sets_rejoined_to_raw_xml': True, 'both_api_conditions_rejoined_to_raw_xml': True, 'production_adopted': False},
    'issues': issues, 'pr': pr, 'ci': ci, 'external_domains': domains,
    'model_requirement_ledger': ledger, 'additional_acceptance_and_gates': [{'id': n, 'status': s, 'evidence_and_next_action': why} for n, s, why in gates],
    'programmes': programmes, 'next_package': next_package,
    'existing_work': {'github_in_progress_and_queued': active, 'local': 'Live Linux and scoped Windows process inventories show no Sirius work; the Windows audit query itself was excluded. No matching Sirius scheduled task or user crontab found. Existing Alaris user services/timers are unrelated and untouched.', 'scope': 'Observed host/repository jobs; no claim about unobserved external hosts or hidden schedules'},
    'previous_goal_turn_classification': {'classification': 'no_progress', 'basis': 'Conservative classification of the supplied final handoff/pause message: it supplies status and no live execution handle. The exact prior unseen turn boundary is unavailable. Earlier described work is corroborated progress, not a current verified wait.', 'this_turn': 'Progress: new independent artifact/raw-XML verification, live readiness/physical-device recovery, full acceptance classification and preserved next package'},
    'git_boundary': {'source_worktree_clean': True, 'branch_and_main_match_live_remote': True, 'worktrees': checked(['git', 'worktree', 'list', '--porcelain']), 'remote_refs': remote_refs, 'source_implementation_committed_and_pushed': True, 'review_evidence_tag': TAG},
}
payload = (json.dumps(report, indent=2, sort_keys=True) + '\n').encode()
(OUTPUT / 'review.json').write_bytes(payload)
print(json.dumps({'report': str(OUTPUT / 'review.json'), 'bytes': len(payload), 'model_entries': len(ledger), 'extra_gates': len(gates), 'statuses': dict(collections.Counter(x['status'] for x in ledger)), 'source_unchanged': checked(['git', 'status', '--porcelain=v1']) == ''}))

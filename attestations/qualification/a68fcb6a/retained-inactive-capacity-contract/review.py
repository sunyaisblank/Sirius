"""Close finite source/control evidence without accepting failed resource owners."""
from pathlib import Path
import hashlib, importlib.util, json, os, re, stat, subprocess, sys, xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'out/retained-inactive-capacity-contract'
C = R / 'out/retained-inactive-capacity-check'
S = R / 'out/retained-inactive-lsan-guard'
REV = 'a68fcb6abade554e92ce92a7cc2a948250574f44'
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == REV
assert subprocess.check_output(['git', 'status', '--porcelain']) == b''
spec = importlib.util.spec_from_file_location('guard', C / 'owned_guard.py')
guard = importlib.util.module_from_spec(spec)
spec.loader.exec_module(guard)
union = {}
owners = {}
for root in [B, C, S]:
    for p in sorted(root.glob('*/owner.json')):
        o = json.loads(p.read_text())
        assert o['status'] == 'terminal' and o['source'] == {'revision': REV, 'status': ''}
        assert o['source_unchanged'] and o['whole_inputs_unchanged']
        assert o['cleanup_errors'] == [] and o['remaining_owned_processes'] == []
        assert o['observed_births'] and all(b['absent'] for b in o['observed_births'])
        before = json.loads((p.parent / 'inputs-before.json').read_text())
        assert before == json.loads((p.parent / 'inputs-after.json').read_text())
        for name, pin in before.items():
            assert name not in union or union[name] == pin, name
            union[name] = pin
        for birth in o['observed_births']:
            try:
                current = guard.proc(birth['pid'])
            except (FileNotFoundError, ProcessLookupError):
                continue
            assert current['start_ticks'] != birth['saved_start_ticks'], birth
        accepted = not ('gcc-sanitize' == p.parent.name and root in [C, S])
        if accepted:
            assert o['passed'] and o['child_exit'] == 0 and o['stop_reason'] is None
        else:
            assert not o['passed'] and o['stop_reason'] == 'diagnostic_sampled_rss_guard'
            assert o['child_exit'] == (-15 if root == C else 0)
        owners[str(p.relative_to(R))] = {
            'governed_run_passed': o['passed'], 'child_exit': o['child_exit'],
            'stop_reason': o['stop_reason'], 'seconds': o['wall_seconds'],
            'peak_charged_rss_kib': o['peak_sampled_rss_kib'],
        }
assert len(owners) == 10 and sum(o['governed_run_passed'] for o in owners.values()) == 8
for name, pin in union.items():
    assert guard.identity(R / name) == pin, name

# Reuse the original numerical and full-tail reader unchanged from its XML step.
# This evaluates only observations, never the failed monitor acceptance predicate.
reader = (C / 'review.py').read_text()
start = reader.index("x=ET.parse(F/'gtest.xml')")
end = reader.index("\nresult={'passed'")
env = {'R': R, 'W': W, 'F': S / 'gcc-sanitize'}
exec('import hashlib,json,re,sys,xml.etree.ElementTree as ET,zipfile\n' + reader[start:end], env)
san = json.loads((S / 'gcc-sanitize/owner.json').read_text())
ctx = json.loads((S / 'gcc-sanitize/execution-context.json').read_text())
assert ctx['test_preset_environment'] == {
    'ASAN_OPTIONS': 'detect_leaks=1:halt_on_error=1:strict_string_checks=1',
    'LSAN_OPTIONS': 'suppressions=' + str(R / 'tests/sanitizers/lsan-vulkan.supp') + ':print_suppressions=1',
    'UBSAN_OPTIONS': 'halt_on_error=1:print_stacktrace=1',
}
err = (S / 'gcc-sanitize/stderr.log').read_text()
assert not re.search(r'ERROR: (?:AddressSanitizer|LeakSanitizer)|runtime error:|SUMMARY: (?:AddressSanitizer|LeakSanitizer)', err)
samples = [json.loads(l) for l in (S / 'gcc-sanitize/samples.jsonl').read_text().splitlines()]
assert san['shared_VM_samples'] == sum(len(s['address_space_accounting']['groups']) < len(s['members']) for s in samples)
over = [s for s in samples if s['rss_kib'] > san['limits']['rss_kib']]
assert len(over) == 1 and over[0] == samples[-1]
last = samples[-1]
assert len(last['members']) == 2 and len(last['address_space_accounting']['comparisons']) == 1
cmp = last['address_space_accounting']['comparisons'][0]
assert cmp['refusal'] == 'FileNotFoundError' and cmp['errno'] == 2 and not cmp['identities_stable']
previous = samples[-2]
assert { (x['pid'], x['start_ticks']) for x in previous['members'] } == { (x['pid'], x['start_ticks']) for x in last['members'] }
assert len(previous['address_space_accounting']['groups']) == 1
assert previous['address_space_accounting']['comparisons'][0]['kcmp_return'] == 0
san_observations = {
    'numerical_and_full_tail_observations_passed': True,
    'test_and_sanitizer_process_exit': san['child_exit'],
    'preset_sanitizers_enabled': ctx['test_preset_environment'],
    'sanitizer_diagnostics_found': False,
    'governed_run_passed': False,
    'resource_envelope_accepted': False,
    'source_revision': REV,
    'all4_original_precision_routes_exercised': True,
    'original56_complete_active_readbacks_equal': True,
    'full_inactive_tail_checks': sum(r['full_inactive_tail_checks'] for r in env['routes']),
    'nonzero_tail_witnesses': sum(r['nonzero_tail_witnesses'] for r in env['routes']),
    'routes': env['routes'], 'execution_context': ctx,
    'seconds': san['wall_seconds'],
    'shared_VM_samples': san['shared_VM_samples'],
    'peak_before_terminal_uncertain_sample_kib': max(s['rss_kib'] for s in samples[:-1]),
    'last_confirmed_shared_VM_sample': previous,
    'terminal_uncertain_sample': last,
    'verdict': 'Actual child process completed with exit 0 and enabled ASan/UBSan/LSan. Conservative resource owner failed at the terminal procfs disappearance race. Preserve the failure; do not infer complete governed acceptance or repair it retroactively. No further unchanged C1 rerun.',
}
(W / 'sanitizer-observations.json').write_text(json.dumps(san_observations, indent=2) + '\n')
builds = json.loads((C / 'build-product-verification.json').read_text())
accept = {n: json.loads((C / (n + '-acceptance.json')).read_text()) for n in ['gcc', 'clang']}
assert all(a['passed'] and a['full_inactive_tail_checks'] == 56 and a['nonzero_tail_witnesses'] == 36 for a in accept.values())
assert san_observations['full_inactive_tail_checks'] == 56 and san_observations['nonzero_tail_witnesses'] == 36
result = {
    'schema': 'sirius-inactive-endpoint-capacity-finite-review-v1',
    'source_revision': REV, 'source_change': json.loads((B / 'source-change-context.json').read_text()),
    'review_completed': True, 'owners': owners,
    'accepted_governed_build_and_control_commands': 8,
    'failed_sanitizer_resource_owners_preserved': 2,
    'final_input_union_rehashed': len(union),
    'all7_executables_built_each_profile': builds['all7_executables_built_per_profile'],
    'all35_SPVs_and3_headers_each_profile_unchanged': builds['all35_SPVs_and3_headers_unchanged_from630844e7'],
    'registration_only_count_each_profile': 1152,
    'registration_names_equal_previous_revision_and_all_profiles': True,
    'governed_C1_passed_profiles': accept,
    'sanitizer_functional_observations': 'sanitizer-observations.json',
    'totals_in_final3_C1_executions': {
        'precision_routes': 12, 'complete_original_active_readbacks_equal': 168,
        'full_inactive_tail_comparisons': 168, 'positive_nonzero_tail_witnesses': 108,
        'does_not_imply_all3_governed_runs_passed': True,
    },
    'shader_or_rendering_implementation_changed': False,
    'rejected_eight_ray_prototype_evidence': '8cebd8e1cd32aac9a23a2b0378fcbf1984ba31c2',
    'prototype_benefit_gate': '12/20 whole-and-stage plus 1/2 API: rejected, no unchanged retry',
    'full_Mandatory_execution': False, 'all8_integrated_controls': False,
    'scientific_qualification': False, 'original_software_frame': False,
    'native_scene_and_window_qualification': False, 'external_domain_admission': False,
    'installed_release_qualification': False, 'ultimate_goal_complete': False,
}
(W / 'root-review.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({k: result[k] for k in ['review_completed', 'accepted_governed_build_and_control_commands', 'failed_sanitizer_resource_owners_preserved', 'final_input_union_rehashed', 'totals_in_final3_C1_executions', 'ultimate_goal_complete']}))

"""Record the exact finite scope of the two original fixtures, without reruns."""
from pathlib import Path
import hashlib
import json

R = Path.cwd()
W = Path(__file__).resolve().parent


def doc(p):
    return json.loads(p.read_text())


def pin(p):
    with p.open('rb') as stream:
        return dict(bytes=p.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


for phase in ('finish-build', 'finish-link', 'run'):
    owner = doc(W / phase / 'owner.json')
    assert owner['status'] == 'terminal' and owner['passed']
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
results = {}
for stage in ('transport', 'endpoint'):
    text = doc(W / (stage + '-output.json'))
    assert not text['stderr']
    result = json.loads(text['stdout'])
    assert result['passed']
    results[stage] = result
t = results['transport']
for key, expected in dict(case_count=38, active_rows=38, capacity_rows=40,
        packets_per_implementation=5, raw_calls=10, compared_output_words=187960,
        mismatch_count=0, first_mismatch=None).items():
    assert t[key] == expected, key
assert t['input_unchanged'] and t['inactive_sentinels_preserved'] and t['full_word_match']
folder = W / 'run/transport'
assert pin(folder / 'input-before.bin') == pin(folder / 'input-after.bin')
assert pin(folder / 'baseline.bin') == pin(folder / 'candidate.bin')
assert (folder / 'baseline.bin').stat().st_size == 187960 * 4
e = results['endpoint']
for key, expected in dict(selection_rows=23, capacity=32, header_values=6, layouts=2,
        reference_rows=138, scalar_packets=360, raw_packets=24,
        other_program_word_refusals=21142, boundary_refusals=21).items():
    assert e[key] == expected, key
for key in ('whole3829_words_and_inactive_capacity_and_canaries', 'input_immutable',
            'MXCSR_controls_unchanged', 'caller_MXCSR_state_restored'):
    assert e[key]
code = doc(W / 'code-proof.json')
assert code['passed'] and len(code['objects']) == 3 and code['no_optional_ISA_in_adapter']
proofs = {}
for stage in ('transport', 'endpoint'):
    proof = doc(W / (stage + '-link-proof.json'))
    assert proof['passed'] and all(row['all_retained_function_opcodes_match'] for row in proof['objects'])
    proofs[stage] = dict(objects=len(proof['objects']),
        functions=sum(len(row['functions']) for row in proof['objects']),
        exact_GOT_fields=len(proof['declared_GOTPCREL_fields_resolve_exact_targets']))
assert proofs['transport']['exact_GOT_fields'] == 2 and proofs['endpoint']['exact_GOT_fields'] == 0
result = dict(finite_component_passed=True, provider_adopted=False, results=results,
    code=code, complete_link_association=proofs,
    original_failed_setup_phases_preserved=['build', 'link'],
    compilation_runs=1, original_raw3_math3_and_fixtures_recompiled=False,
    numerical_runs_per_fixture=1,
    qualification_scope='Raw-v4 Transport38 and focused Endpoint finite complete-output compatibility only, direct raw entry under actual-CPU constructor guard. Stage selection, fallback and direct-helper edges static checked; DopriPhase compiled and linked but not executed. No throughput, provider/C8/all8/ABBA/frame, full scientific or release acceptance.',
    next_required_gate='Provider composition plus fresh unchanged C8 must establish new packet-level selection and actual DopriPhase execution before adoption. Retain original AP failed frame identity; no unchanged renderer retry.')
with (W / 'finite-result.json').open('x') as stream:
    stream.write(json.dumps(result, indent=2) + '\n')
print(json.dumps(dict(passed=True, transport_cases=38, endpoint_raw_packets=24,
                     link_association=proofs, provider_adopted=False)))
